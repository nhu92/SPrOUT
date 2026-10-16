"""Merge frozen C, gene rescue/stability and ML evidence without truth access.

The auxiliary route must retain all three checks: existing C membership, stable
gene rescue, and at least one independently counted ML-supporting gene. Review
candidates never silently enter the final prediction set.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED


DEFAULT_POLICY = Path(__file__).with_name("fusion_policy.json")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    return rows


def keyed(rows, key):
    result = {r[key]: r for r in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate {key}")
    return result


def boolean(value):
    if value not in {"True", "False"}:
        raise ValueError(f"Expected serialized boolean, got {value!r}")
    return value == "True"


def decide(row, policy):
    votes = row["ml_gene_votes"]
    if type(votes) is not int or votes < 0:
        raise ValueError("ML votes must be a nonnegative gene count")
    core = row["c_selected"] and votes >= policy["core_ml_genes"]
    auxiliary = (row["c_selected"] and row["rescue_selected"]
                 and row["stability_lower"] > policy["auxiliary_stability_lower_exclusive"]
                 and votes >= policy["auxiliary_ml_genes"])
    merged = core or auxiliary
    review = (row["b_selected"] or row["c_selected"]
              or votes >= policy["standalone_ml_genes"]) and not merged
    route = "core_and_auxiliary" if core and auxiliary else "core" if core else "auxiliary" if auxiliary else "review" if review else "not_selected"
    return {"core": core, "auxiliary": auxiliary, "merged": merged, "review": review, "route": route}


def csv_bytes(rows):
    stream = io.StringIO(newline="")
    w = csv.DictWriter(stream, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
    return stream.getvalue().encode("utf-8-sig")


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def collect(transfer, legacy, panels, policy):
    """Validate source agreement and join evidence; no sample truth is read."""
    transfer, legacy, panels = map(Path, (transfer, legacy, panels))
    sources = {}

    def record(path):
        sources[str(path.resolve())] = digest(path)
        return path

    def names(path):
        values = record(path).read_text(encoding="utf-8").splitlines()
        if len(set(values)) != len(values) or any(not v.strip() for v in values):
            raise ValueError(f"Duplicate/blank candidate: {path}")
        return set(values)

    source_audit = json.loads(record(transfer / "AUDIT.json").read_text())
    if not all(source_audit[k] for k in ("all_100_complete", "all_reused_input_hashes_unchanged", "code_model_protocol_unchanged")):
        raise ValueError("Transfer input audit did not pass")
    record(transfer / "PREDICTIONS_COMPLETE")
    selection = json.loads(record(transfer / "selection.json").read_text())["selected"]
    old_selection = json.loads(record(legacy / "selection.json").read_text())["selected"]
    if len(selection) != 100 or {r["project"] for r in selection} != {r["project"] for r in old_selection}:
        raise ValueError("The source sample sets do not match")
    old_by_project = keyed(old_selection, "project")
    rows, calls = [], []
    for sample in selection:
        project = sample["project"]
        if sample["sample"] != old_by_project[project]["sample"]:
            raise ValueError("Source sample names do not match")
        ml_path = transfer / "samples" / project
        old_path = legacy / "samples" / project
        stat_path = panels / project
        record(ml_path / "SUCCESS")
        audit = json.loads(record(ml_path / "audit.json").read_text())
        if not math.isclose(audit["threshold"], policy["ml_match_cutoff"], rel_tol=0, abs_tol=1e-12) or audit["truth_used"]:
            raise ValueError("Unexpected ML model threshold/provenance")
        ml = keyed(read_csv(record(ml_path / "family_evidence.csv")), "family")
        old = keyed(read_csv(record(old_path / "family_evidence.csv")), "taxon")
        stat = keyed(read_csv(record(stat_path / "evidence.csv")), "family")
        c = names(old_path / "family_candidates.txt")
        if c != names(ml_path / "baseline_confirmed.txt") or c != names(stat_path / "candidates.txt"):
            raise ValueError(f"Frozen C lists disagree: {project}")
        stable = names(stat_path / "stable_majority_candidates.txt")
        method_sets = {m: set() for m in ("baseline_c", "stable_c", "core", "auxiliary", "merged", "review", "standalone_ml")}
        if set(ml) != set(old) or set(ml) != set(stat):
            raise ValueError(f"Family panels disagree: {project}")
        for family in sorted(ml):
            m, o, s = ml[family], old[family], stat[family]
            lower, upper = json.loads(s["selection_frequency_mc_interval"])
            if not 0 <= lower <= upper <= 1:
                raise ValueError("Invalid stability interval")
            in_c, rescue, b = boolean(o["selected"]), boolean(o["rescue_selected"]), boolean(o["hybrid_selected"])
            if (in_c != (family in c) or in_c != boolean(m["baseline_confirmed"])
                    or in_c != boolean(s["selected"]) or rescue != boolean(s["rescue_selected"])
                    or b != boolean(s["hybrid_selected"]) or (rescue and not in_c)):
                raise ValueError(f"Admission evidence disagrees: {project}/{family}")
            if (family in stable) != (in_c and lower > 0.5):
                raise ValueError("Stable-C list disagrees with evidence")
            if m["parent_order"] != o["parent_order"] or m["parent_order"] != s["parent_order"]:
                raise ValueError("Parent labels disagree")
            row = {"sample": sample["sample"], "project": project, "family": family,
                   "parent_order": m["parent_order"], "c_selected": in_c,
                   "b_selected": b, "rescue_selected": rescue,
                   "native_z_selected": boolean(o["native_z_selected"]),
                   "ml_gene_votes": int(m["profile_gene_votes"]),
                   "nearest_gene_votes": int(m["nearest_gene_votes"]),
                   "valid_genes": int(m["valid_genes"]), "max_ml_match_score": float(m["max_match_score"]),
                   "original_z": float(o["original_z"]), "discriminating_genes": int(o["discriminating_genes"]),
                   "comparable_genes": int(o["comparable_genes"]), "bootstrap_lower": float(o["bootstrap_lower"]),
                   "stability_frequency": float(s["selection_frequency"]), "stability_lower": lower,
                   "max_gene_share": float(s["max_gene_share"]), "stable_c": family in stable}
            row.update(decide(row, policy))
            if row["core"] != boolean(m["profile_confirmed_baseline"]):
                raise ValueError("Core does not reproduce previous ML confirmation")
            for method, chosen in (("baseline_c", in_c), ("stable_c", family in stable),
                                   ("standalone_ml", row["ml_gene_votes"] >= policy["standalone_ml_genes"]),
                                   *[(k, row[k]) for k in ("core", "auxiliary", "merged", "review")]):
                if chosen:
                    method_sets[method].add(family)
            rows.append(row)
        calls.append({"sample": sample["sample"], "project": project,
                      "historical30": sample["previously_seen"],
                      **{k: sorted(v) for k, v in method_sets.items()}})
    for path, expected in sources.items():
        if digest(path) != expected:
            raise ValueError(f"Source changed during merge: {path}")
    return rows, calls, sources


def build(transfer, legacy, panels, output, policy_path=DEFAULT_POLICY):
    output = Path(output)
    if output.exists():
        raise ValueError("Use a new output directory; previous predictions are immutable")
    policy_hash = digest(policy_path)
    policy = json.loads(Path(policy_path).read_text())
    rows, calls, sources = collect(transfer, legacy, panels, policy)
    if digest(policy_path) != policy_hash:
        raise ValueError("Policy changed during merge")
    output.mkdir(parents=True)
    write_json(output / "POLICY.json", policy)
    write_json(output / "predictions.json", calls)
    (output / "evidence.csv.gz").write_bytes(gzip.compress(csv_bytes(rows), mtime=0))
    active = [r for r in rows if r["merged"] or r["review"]]
    (output / "candidate_evidence.csv").write_bytes(csv_bytes(active))
    lists = [{"sample": c["sample"], "family": f, "route": next(r["route"] for r in rows if r["sample"] == c["sample"] and r["family"] == f)} for c in calls for f in c["merged"]]
    (output / "merged_predictions.csv").write_bytes(csv_bytes(lists))
    with ZipFile(output / "prediction_sets.zip", "w", compression=ZIP_DEFLATED) as archive:
        for sample in calls:
            for method in ("core", "auxiliary", "merged", "review"):
                archive.writestr(f'{sample["project"]}/{method}_families.txt', "".join(f + "\n" for f in sample[method]))
        archive.writestr("POLICY.json", (output / "POLICY.json").read_bytes())
    artifact_hashes = {p.name: digest(p) for p in output.iterdir() if p.is_file()}
    write_json(output / "MANIFEST.json", {"schema": policy["schema"], "samples": len(calls),
        "evidence_rows": len(rows), "source_hashes": sources, "code_sha256": digest(__file__),
        "policy_sha256": policy_hash, "artifact_sha256": artifact_hashes,
        "source_hashes_unchanged": True, "truth_read_by_predictor": False,
        "final_prediction_set": "merged", "review_is_part_of_final": False,
        "retrospective": True, "candidate_count_cap": None})
    (output / "PREDICTIONS_COMPLETE").write_text("100/100 complete; no truth input\n")
    print(json.dumps({"samples": len(calls), "merged_events": sum(len(c["merged"]) for c in calls),
                      "auxiliary_only_events": sum(len(set(c["merged"]) - set(c["core"])) for c in calls),
                      "review_events": sum(len(c["review"]) for c in calls)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transfer-batch", required=True)
    parser.add_argument("--legacy-batch", required=True)
    parser.add_argument("--statistical-panels", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--policy", default=str(DEFAULT_POLICY))
    args = parser.parse_args()
    build(args.transfer_batch, args.legacy_batch, args.statistical_panels, args.output, args.policy)
