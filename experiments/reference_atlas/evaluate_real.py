"""Score sealed real-mixture calls against a separately audited S7 mapping CSV.

This local evaluator is never installed on the prediction host. The mapping must
contain sample, source_cell, species, family and order, with one species per row.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path

from .real_truth import verify_complete_batch
from .real_evaluation import score, family_scores, audit_exported_counts, prediction
from .export_real_predictions import sha


def evaluate(batch, mapping, output, sealed_output=None):
    batch, mapping, output = Path(batch), Path(mapping), Path(output)
    calls = verify_complete_batch(batch)
    sealed_output = Path(sealed_output) if sealed_output else output
    seal = json.loads((sealed_output / "PREDICTION_EXPORT.json").read_text())
    for name, digest in seal["artifact_sha256"].items():
        if sha(sealed_output / name) != digest:
            raise ValueError("Sealed prediction artifact changed: " + name)
    if sha(batch / "predictions.json") != sha(sealed_output / "predictions.json"):
        raise ValueError("Batch differs from sealed predictions")
    with mapping.open(encoding="utf-8-sig", newline="") as f:
        mapped = list(csv.DictReader(f))
    truth = {}
    for row in mapped:
        if not all(row.get(k) for k in ("sample", "source_cell", "species", "family", "order")):
            raise ValueError("Incomplete truth mapping row")
        truth.setdefault(row["sample"], set()).add(row["family"])
    with (batch / "evidence.csv").open(encoding="utf-8-sig", newline="") as f:
        evidence = list(csv.DictReader(f))
    lookup = {(r["sample"], r["family"]): r for r in evidence}
    references, coverage, excluded = set(), [], []
    # Include every family actually represented in scoring trees/matrices, as
    # well as selected panel recoveries that later lacked usable query evidence.
    references.update(r["family"] for r in evidence)
    for sample in calls:
        folder = batch / "samples" / sample["project"]
        panel = json.loads((folder / "02_family/family_panel_report.json").read_text())
        families = {r["family"] for r in panel["selection"]}
        references.update(families)
        trees = json.loads((folder / "02_family/family_tree_validation.json").read_text())
        audit = json.loads((folder / "AUDIT.json").read_text())
        parents = set(panel["parent_taxa"])
        for family in sorted(truth[sample["sample"]]):
            expected_orders = {r["order"] for r in mapped if r["sample"] == sample["sample"] and r["family"] == family}
            row = lookup.get((sample["sample"], family), {})
            coverage.append({"sample": sample["sample"], "family": family,
                "truth_orders": ";".join(sorted(expected_orders)), "parent_admitted": bool(expected_orders & parents),
                "reference_present": family in families, "scoring_evidence_present": bool(row),
                "statistical_coverage": row.get("statistical_coverage", "False"),
                "ml_valid_genes": row.get("valid_genes", 0), "valid_genes_sample": audit["valid_genes"],
                "valid_trees_sample": audit["valid_trees"], "merged_selected": family in sample["merged"]})
        excluded.append({"sample": sample["sample"], "valid_genes": audit["valid_genes"],
            "valid_trees": audit["valid_trees"], "tree_validation": trees})
    by_sample, summary, universe = score(calls, truth, references)
    audited = audit_exported_counts(batch, truth, universe, by_sample, summary)
    changes, errors = [], []
    for sample in calls:
        sid, expected = sample["sample"], truth[sample["sample"]]
        if not set(sample["merged"]) <= set(sample["baseline_c"]):
            raise ValueError("Merged calls must be a subset of C")
        for comparison, source, dest in (("C_to_merged", "baseline_c", "merged"), ("core_to_merged", "core", "merged")):
            before, after = prediction(sample, source), prediction(sample, dest)
            for direction, families in (("removed", before - after), ("added", after - before)):
                for family in sorted(families):
                    changes.append({"sample": sid, "comparison": comparison, "direction": direction,
                        "family": family, "true_family": family in expected, **lookup.get((sid, family), {})})
        for method in ("baseline_c", "merged", "standalone_ml", "expanded_review_union"):
            chosen = prediction(sample, method)
            for error, families in (("FP", chosen - expected), ("FN", expected - chosen)):
                for family in sorted(families):
                    errors.append({**lookup.get((sid, family), {}), "sample": sid, "method": method,
                        "error": error, "family": family, "in_review": family in sample["review"]})
    output.mkdir(parents=True, exist_ok=True)
    for filename, rows in (("summary.csv", summary), ("by_sample.csv", by_sample),
                           ("truth_coverage.csv", coverage), ("method_changes.csv", changes), ("errors.csv", errors)):
        if rows:
            fields = list(dict.fromkeys(k for row in rows for k in row))
            with (output / filename).open("w", encoding="utf-8-sig", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
    with gzip.open(output / "by_family.csv.gz", "wt", encoding="utf-8", newline="") as f:
        rows = family_scores(calls, truth, universe)
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for name, value in (("truth.json", {k: sorted(v) for k, v in truth.items()}),
                        ("tree_coverage.json", excluded), ("EVALUATION_AUDIT.json", {
                            **audited, "sealed_predictions_unchanged": True,
                            "truth_mapping_sha256": sha(mapping), "truth_events": sum(map(len, truth.values())),
                            "truth_species_rows": len(mapped), "negative_universe_families": sorted(universe)})):
        (output / name).write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sealed-output", help="Original sealed predictions when writing a separate diagnostic view")
    args = parser.parse_args()
    evaluate(args.batch, args.mapping, args.output, args.sealed_output)
