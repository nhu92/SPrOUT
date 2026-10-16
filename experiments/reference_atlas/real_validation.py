"""Frozen family-fusion validation on UMD's eleven real mixtures.

Only prediction-side code belongs on the cluster. S7 and the evaluator remain
local. Query assemblies and complete per-gene order evidence are reused; order
admission is replayed without the legacy capacity limit, using distinct genes.
References and full family trees are rebuilt for all admitted parents.
"""
import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import traceback
from collections import defaultdict


MODEL_HASH = "ae8673bd3a42a4220abcd35fed00f60c5e93c7ad69c9ae77cb2ae021ee3247d8"
LEGACY_ROOT = Path("/scratch/zt1/project/jmeng-prj/user/nhu/mrknowitall_full_20261014/realmix_runs")
SAMPLE_ROOT = Path("/scratch/zt1/project/jmeng-prj/user/nhu/202610/realmixes")
BACKBONE = Path("/home/nhu/mrknowitall_test_20261013/backbone_corrected")
CACHE = Path("/scratch/zt1/project/jmeng-prj/user/nhu/mrknowitall_hierarchy_20261013/kew_cache")
PACKAGES = ("numpy", "scipy", "scikit-learn", "biopython", "joblib", "pandas")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def write_names(path, values):
    Path(path).write_text("".join(f + "\n" for f in sorted(values)))


def verify_files(hashes):
    for path, expected in hashes.items():
        if sha(path) != expected:
            raise RuntimeError("Frozen input changed: " + path)


def source_hashes(root):
    paths = [p for folder in (root / "code", root / "sources")
             for p in sorted(folder.rglob("*")) if p.is_file() and "__pycache__" not in p.parts]
    paths += sorted(root.glob("*.sbatch")) + [root / "REQUEST.json"]
    return {str(p.resolve()): sha(p) for p in paths}


def check_order_rule(evidence, candidates):
    """Legacy order admission must equal the existing uncapped 2-gene/Z/top rule."""
    admitted = set()
    for r in evidence:
        z = float(r["aggregate_z"]) if r.get("aggregate_z") else float("-inf")
        rank = int(r["aggregate_rank"]) if r.get("aggregate_rank") else None
        if int(r["locus_votes"]) >= 2 or z >= .9 or rank == 1:
            admitted.add(r["taxon"])
    if admitted != candidates:
        raise ValueError("Legacy order calls do not match the frozen uncapped rule")
    if any("capacity_excluded" in str(r) for r in evidence):
        raise ValueError("Legacy order evidence contains capacity exclusions")


def uncapped_orders(stage):
    """Recover distinct-gene votes from the complete, unfiltered PEA evidence.

    Legacy aggregate locus_votes may count exons and legacy selected_by fields
    may be capacity-truncated. Neither is used to decide the new parent set.
    """
    genes = [g["gene"] for g in stage["genes"]]
    if len(genes) != len(set(genes)) or len(genes) < 10:
        raise ValueError("Invalid full gene evidence")
    votes = defaultdict(set)
    for gene in stage["genes"]:
        for match in gene["matches"]:
            if int(match["query_votes"]) > 0:
                votes[match["taxon"]].add(gene["gene"])
    predicted = {r["taxon"]: r for r in stage["predictions"]}
    rows = []
    for taxon in sorted(set(votes) | set(predicted)):
        r = predicted.get(taxon, {})
        z, rank = r.get("z_score"), r.get("rank")
        selected = len(votes[taxon]) >= 2 or z is not None and z >= .9 or rank == 1
        rows.append({"taxon": taxon, "distinct_gene_votes": len(votes[taxon]),
            "aggregate_z": z, "aggregate_rank": rank, "selected": bool(selected)})
    return {r["taxon"] for r in rows if r["selected"]}, rows


def prepare(root):
    if (root / "PROTOCOL.json").exists():
        raise ValueError("Refusing to replace a frozen protocol")
    # Fail on missing runtime dependencies before reading thousands of inputs.
    runtime = {p: importlib.metadata.version(p) for p in PACKAGES}
    protocol = read_json(root / "REQUEST.json")
    if sha(root / "code/model.joblib") != MODEL_HASH:
        raise ValueError("The original frozen TTU model is required")
    selection_path = root / "selection.json"
    cached = read_json(selection_path)["selected"] if selection_path.exists() else None
    selection = cached if cached is not None else []
    expected_samples = [f"EMseq{i:03d}" for i in range(1, 12)]
    if sorted(p.name for p in SAMPLE_ROOT.iterdir() if p.is_dir()) != expected_samples:
        raise ValueError("Unexpected real-mixture sample directories")
    if cached is not None:
        if [s["sample"] for s in cached] != expected_samples or any(
            s["index"] != i or s["project"] != s["sample"].casefold()
            or s["run"] != str(root / "samples" / s["project"])
            or s["source_order"] != str(LEGACY_ROOT / s["sample"] / "order/results")
            for i, s in enumerate(cached)):
            raise ValueError("The partial selection cache belongs to a different run")
        # Every source hash is independently rechecked by execute() before its
        # sample is consumed, and again by collect(). This cache only avoids
        # repeating discovery after a preflight-only environment failure.
    for i, sample in enumerate(expected_samples):
        if cached is not None:
            break
        old = LEGACY_ROOT / sample
        order, family = old / "order/results", old / "family/results"
        assert (order / "SUCCESS").exists() and (family / "SUCCESS").exists()
        project = sample.casefold()
        exons = sorted((order / "exon").glob("*.fasta"))
        if len(exons) < 10:
            raise ValueError("Too few assembled exon files: " + sample)
        document = read_json(order / "stage.pea.json")
        order_stage = next(s for s in document["stages"] if s["level"] == "order")
        candidates, order_evidence = uncapped_orders(order_stage)
        old_candidates = set((order / "candidates.txt").read_text().splitlines())
        genes = set((order / "gene_list.txt").read_text().splitlines())
        available = genes & {p.stem for p in BACKBONE.glob("*.fasta")}
        if len(available) < 10:
            raise ValueError("Fewer than ten reference/query genes: " + sample)
        files = exons + [order / name for name in
            ("gene_list.txt", "candidates.txt", "candidate_evidence.csv", "predictions.csv", "stage.pea.json", "RUN_SUMMARY.txt")]
        reuse = {str(p.resolve()): sha(p) for p in files}
        # Store size/mtime only for large reads; no FASTQ is copied or altered.
        read_metadata = [{"path": str(p), "size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns}
                         for p in sorted((SAMPLE_ROOT / sample).iterdir()) if p.is_file()]
        selection.append({"index": i, "sample": sample, "project": project,
            "run": str(root / "samples" / project), "source_order": str(order), "source_family": str(family),
            "source_genes": sorted(available), "orders": sorted(candidates), "reuse_sha256": reuse,
            "order_evidence": order_evidence, "legacy_orders": sorted(old_candidates),
            "legacy_order_limit_recovered": sorted(candidates - old_candidates),
            "reads_metadata": read_metadata, "previously_seen": True})
    write_json(root / "selection.json", {"selected": selection, "truth_used": False})
    external = {str(p.resolve()): sha(p) for p in sorted(BACKBONE.glob("*.fasta"))}
    external.update({str((CACHE / name).resolve()): sha(CACHE / name) for name in
                     ("release.json", "by_recovery_index.json", "sequence_manifest.tsv", "treeoflife.current.tree")})
    protocol.update(selection_sha256=sha(root / "selection.json"), code_and_model_sha256=source_hashes(root),
        runtime=runtime, python=sys.version, partial_selection_cache_reused=cached is not None,
        reference_metadata_and_backbone_sha256=external,
        frozen_model_sha256=MODEL_HASH, truth_uploaded=False, truth_used=False)
    write_json(root / "PROTOCOL.json", protocol)
    print(json.dumps({"samples": len(selection), "reused_exon_and_reference_files": sum(len(s["reuse_sha256"]) for s in selection),
                      "runtime": protocol["runtime"]}))


def prefetch(root):
    """Metadata/download orchestration only; no alignment or scoring on login."""
    sys.path.insert(0, str(root / "sources/MrKnowItAll/src"))
    from sprout_refbuilder.clade import BalancedKewPanelBuilder
    from sprout_refbuilder.kew import KewClient
    protocol = read_json(root / "PROTOCOL.json")
    verify_files(protocol["code_and_model_sha256"])
    selected = read_json(root / "selection.json")["selected"]
    parents = sorted({p for s in selected for p in s["orders"]})
    with KewClient(CACHE) as client:
        report = BalancedKewPanelBuilder(BACKBONE, root / "prefetch_only", representatives_per_taxon=4,
            pool_multiplier=2, minimum_recovered_loci=50, download_workers=8).prefetch("family", parents, client)
    verify_files(protocol["reference_metadata_and_backbone_sha256"])
    write_json(root / "PREFETCH.json", report)
    (root / "PREFETCH_COMPLETE").write_text("download orchestration complete\n")
    print(json.dumps(report), flush=True)


def validate_trees(source, output, project):
    from Bio import Phylo
    output.mkdir(exist_ok=True)
    valid, excluded = [], []
    for path in sorted(source.glob("*.tre")):
        try:
            tree = Phylo.read(path, "newick")
            names = [t.name or "" for t in tree.get_terminals()]
            if len(set(names)) != len(names):
                raise ValueError("duplicate terminal names")
            if not any(project in n for n in names) or not any(project not in n and "NODE" not in n for n in names):
                raise ValueError("query/reference leaves missing")
            link = output / path.name
            if not link.exists():
                shutil.copy2(path, link)
            valid.append(path.name)
        except Exception as e:
            excluded.append({"tree": path.name, "error": str(e)})
    result = {"valid_trees": valid, "excluded_trees": excluded,
              "distinct_genes": len({name.split("_exon_")[0] for name in valid})}
    write_json(output.parent / "family_tree_validation.json", result)
    if result["distinct_genes"] < 10:
        raise ValueError("Fewer than ten valid genes")
    return result


def execute(root, index):
    protocol = read_json(root / "PROTOCOL.json")
    if sha(root / "selection.json") != protocol["selection_sha256"]:
        raise ValueError("Selection changed")
    verify_files(protocol["code_and_model_sha256"])
    verify_files(protocol["reference_metadata_and_backbone_sha256"])
    if not (root / "PREFETCH_COMPLETE").exists():
        raise ValueError("Complete metadata/download orchestration first")
    sample = read_json(root / "selection.json")["selected"][index]
    run = Path(sample["run"])
    run.mkdir(parents=True, exist_ok=True)
    if (run / "SUCCESS").exists():
        print("Already complete:", sample["sample"])
        return
    if (run / "FAILED").exists():
        raise ValueError("Preserve and resolve the failed attempt before a controlled retry")
    verify_files(sample["reuse_sha256"])
    try:
        execute_stages(root, run, sample)
        verify_files(sample["reuse_sha256"])
        verify_files(protocol["code_and_model_sha256"])
        verify_files(protocol["reference_metadata_and_backbone_sha256"])
        (run / "SUCCESS").write_text("complete\n")
        print("COMPLETE", sample["sample"], flush=True)
    except Exception:
        (run / "FAILED").write_text(traceback.format_exc())
        raise


def execute_stages(root, run, sample):
    sys.path[:0] = [str(root / "sources"), str(root / "code")]
    from sprout_evidence.features import extract, save_json
    from sprout_evidence.statistics import analyze, Settings
    from sprout_evidence.frozen.adaptive import extract_family_features, read_predictions
    from sprout_evidence.frozen.precision import confirmed_family_selection
    from transfer import predict_sample
    from fuse import decide
    project = sample["project"]
    stage = run / "02_family"
    stage.mkdir(exist_ok=True)
    order_stage = run / "01_order"
    order_stage.mkdir(exist_ok=True)
    write_names(order_stage / "order_candidates.txt", sample["orders"])
    write_csv(order_stage / "order_evidence.csv", sample["order_evidence"])

    def status(value):
        (run / "CURRENT_STAGE").write_text(value + "\n")
        print(sample["sample"], value, flush=True)

    def command(args):
        subprocess.run([sys.executable, *map(str, args)], cwd=stage, check=True)

    if not (stage / "REFERENCE_COMPLETE").exists():
        status("uncapped_family_reference")
        os.environ["PYTHONPATH"] = str(root / "sources/MrKnowItAll/src")
        os.environ["SPROUT_REF_OFFLINE"] = "1"
        command(["-m", "sprout_refbuilder.cli", "--cache", CACHE, "build-kew-panel", "--level", "family",
            "--parent-file", order_stage / "order_candidates.txt", "--backbone", BACKBONE,
            "--alignment-gene-list", Path(sample["source_order"]) / "gene_list.txt",
            "--output", stage / "family_reference", "--representatives", 4, "--pool-multiplier", 2,
            "--minimum-recovered-loci", 50, "--download-workers", 8, "--threads", 8])
        available = set((stage / "family_reference/gene.list.txt").read_text().splitlines()) & set(sample["source_genes"])
        if len(available) < 10:
            raise ValueError("Reference reconstruction yielded fewer than ten query genes")
        write_names(stage / "gene_list.txt", available)
        shutil.copy2(stage / "family_reference/panel_report.json", stage / "family_panel_report.json")
        (stage / "REFERENCE_COMPLETE").write_text("complete\n")
    reference_hashes = {str(p): sha(p) for p in (stage / "family_reference/ref").glob("*.fasta")}
    if not (stage / "TREES_COMPLETE").exists():
        status("family_trees")
        command([root / "sources/02_exon_trees.py", "-t", 4, "-p", project,
                 "-g", stage / "gene_list.txt", "-e", Path(sample["source_order"]) / "exon",
                 "-r", stage / "family_reference/ref", "--output_phylo", stage / "phylo"])
        validate_trees(stage / "phylo", stage / "valid_phylo", project)
        (stage / "TREES_COMPLETE").write_text("complete\n")
    tree_hashes = {str(p): sha(p) for p in (stage / "valid_phylo").glob("*.tre")}
    if not (stage / "DISTANCES_COMPLETE").exists():
        status("family_distances")
        command([root / "sources/03_distance_matrices.py", "-t", 4, "-p", project,
                 "-g", stage / "gene_list.txt", "--query-only", "--input_phylo", stage / "valid_phylo",
                 "--output_tree", stage / "family_matrix"])
        command([root / "sources/04_prediction.py", "-i", stage / f"{project}.cumulative_dist.csv",
                 "-o", stage / "native_predictions.csv", "-tl", "f", "-z", .9, "-to", stage / "native_z_candidates.txt"])
        matrices = sorted((stage / "family_matrix").glob("*.cleaned.csv"))
        if len(matrices) != len(tree_hashes):
            raise ValueError("Missing distance matrices for valid trees")
        (stage / "DISTANCES_COMPLETE").write_text("complete\n")
    status("frozen_C_and_stability")
    result = confirmed_family_selection(extract_family_features(stage / "family_matrix", project),
        read_predictions(stage / "native_predictions.csv"), z_threshold=.9)
    write_csv(stage / "family_evidence.csv", result["rows"])
    write_names(stage / "family_candidates.txt", result["selected"])
    write_names(stage / "prefilter_candidates.txt", result["prefilter_selected"])
    write_names(stage / "unconfirmed_candidates.txt", result["unconfirmed"])
    features = extract(stage / "family_matrix", stage / "native_predictions.csv", project)
    save_json(stage / "features.json.gz", features)
    statistics = analyze(features, Settings())
    save_json(stage / "statistics.json", statistics)
    raw_c = {r["taxon"]: r for r in result["rows"]}
    stat = {r["family"]: r for r in statistics["rows"]}
    if set(raw_c) != set(stat):
        raise ValueError("C/statistical family panels differ")
    for family in stat:
        if raw_c[family]["selected"] != stat[family]["selected"] or raw_c[family]["rescue_selected"] != stat[family]["rescue_selected"]:
            raise ValueError("Frozen C and statistical reconstruction differ: " + family)
    status("frozen_ML")
    predict_sample(root, root / "code", root / "ml", sample["index"])
    ml = {r["family"]: r for r in read_csv(root / "ml/samples" / project / "family_evidence.csv")}
    # Trees and raw matrices can have different reference coverage. Never fill
    # a missing statistical family with fabricated support or biological absence.
    if not set(raw_c) <= set(ml):
        raise ValueError("ML evidence missing a C reference family")
    policy = read_json(root / "code/fusion_policy.json")
    rows = []
    calls = {k: [] for k in ("original_z", "baseline_c", "stable_c", "core", "auxiliary", "merged", "review", "standalone_ml")}
    for family in sorted(set(ml) | set(raw_c)):
        m, c, s = ml.get(family, {}), raw_c.get(family, {}), stat.get(family, {})
        record = {"sample": sample["sample"], "project": project, "family": family,
            "parent_order": m.get("parent_order", c.get("parent_order", "")),
            "c_selected": bool(c.get("selected", False)), "b_selected": bool(c.get("hybrid_selected", False)),
            "rescue_selected": bool(c.get("rescue_selected", False)), "ml_gene_votes": int(m.get("profile_gene_votes", 0)),
            "nearest_gene_votes": int(m.get("nearest_gene_votes", 0)), "valid_genes": int(m.get("valid_genes", 0)),
            "original_z": c.get("original_z"), "discriminating_genes": c.get("discriminating_genes"),
            "comparable_genes": c.get("comparable_genes"), "bootstrap_lower": c.get("bootstrap_lower"),
            "stability_lower": s.get("selection_frequency_mc_interval", [0])[0],
            "stability_frequency": s.get("selection_frequency"), "max_gene_share": s.get("max_gene_share"),
            "statistical_coverage": family in stat}
        record.update(decide(record, policy))
        rows.append(record)
        flags = {k: record[k] for k in ("core", "auxiliary", "merged", "review")}
        flags.update(original_z=bool(c.get("native_z_selected", False)), baseline_c=record["c_selected"],
            stable_c=record["c_selected"] and record["stability_lower"] > .5,
            standalone_ml=record["ml_gene_votes"] >= 3)
        for name, selected in flags.items():
            if selected:
                calls[name].append(family)
    for name, values in calls.items():
        write_names(run / f"{name}_families.txt", values)
    write_csv(run / "merged_evidence.csv", rows)
    write_json(run / "predictions.json", {"sample": sample["sample"], "project": project, **calls})
    verify_files(tree_hashes)
    verify_files(reference_hashes)
    verify_files(features["input_sha256"])
    status("complete")
    write_json(run / "AUDIT.json", {"truth_used": False, "code_model_hashes_verified": True,
        "reused_assemblies_and_references_unchanged": True, "scoring_tree_hashes_unchanged": True,
        "tree_sha256": tree_hashes, "distance_sha256": features["input_sha256"],
        "reference_sha256": reference_hashes,
        "valid_genes": len(features["genes"]), "valid_trees": len(tree_hashes), "candidate_count_cap": None,
        "rebuilt_family_trees": True, "original_C_matches_shared_gene_reconstruction": True})


def collect(root):
    protocol = read_json(root / "PROTOCOL.json")
    verify_files(protocol["code_and_model_sha256"])
    verify_files(protocol["reference_metadata_and_backbone_sha256"])
    samples = read_json(root / "selection.json")["selected"]
    if not all((Path(s["run"]) / "SUCCESS").exists() and not (Path(s["run"]) / "FAILED").exists() for s in samples):
        raise ValueError("All eleven must succeed before collection or evaluation")
    calls, rows, audits = [], [], []
    for sample in samples:
        verify_files(sample["reuse_sha256"])
        run = Path(sample["run"])
        audit = read_json(run / "AUDIT.json")
        verify_files(audit["tree_sha256"])
        verify_files(audit["distance_sha256"])
        verify_files(audit["reference_sha256"])
        audits.append({"sample": sample["sample"], **audit})
        calls.append(read_json(run / "predictions.json"))
        rows.extend(read_csv(run / "merged_evidence.csv"))
    write_json(root / "predictions.json", calls)
    write_csv(root / "evidence.csv", rows)
    write_json(root / "AUDIT.json", {"all_11_complete": True, "truth_used": False,
        "frozen_code_model_policy_and_inputs_unchanged": True, "samples": audits})
    (root / "PREDICTIONS_COMPLETE").write_text("11/11 complete before S7 evaluation\n")
    with tarfile.open(root / "compact_results.tar.gz", "w:gz") as archive:
        for name in ("PROTOCOL.json", "REQUEST.json", "selection.json", "predictions.json", "evidence.csv", "AUDIT.json", "PREDICTIONS_COMPLETE"):
            archive.add(root / name, arcname=name)
        for sample in samples:
            run = Path(sample["run"])
            for p in sorted(run.iterdir()):
                if p.is_file():
                    archive.add(p, arcname=str(p.relative_to(root)))
            for name in ("family_evidence.csv", "family_tree_validation.json", "family_panel_report.json", "native_predictions.csv", "statistics.json"):
                archive.add(run / "02_family" / name, arcname=str((run / "02_family" / name).relative_to(root)))
            for name in ("family_evidence.csv", "exon_predictions.csv", "audit.json"):
                p = root / "ml/samples" / sample["project"] / name
                archive.add(p, arcname=str(p.relative_to(root)))
    print(json.dumps({"complete": 11, "artifact": str(root / "compact_results.tar.gz")}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "prefetch", "run", "collect"))
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.command == "prepare":
        prepare(root)
    elif args.command == "prefetch":
        prefetch(root)
    elif args.command == "collect":
        collect(root)
    else:
        if args.index is None:
            parser.error("run requires --index")
        execute(root, args.index)
