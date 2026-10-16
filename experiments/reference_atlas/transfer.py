"""Reference-trained full-profile family scoring on immutable exon trees.

The transfer model learns a shared candidate-family match score from known
reference pseudo-queries. It never consumes mixture identities or truth labels.
All references of a pseudo-query species are removed from its anchor library.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import joblib
import numpy as np
from Bio import Phylo
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

FEATURES = ["nearest", "mean_distance", "sibling_margin", "centered_profile_rmse",
            "profile_correlation", "profile_offset", "distance_rank"]
METHODS = ("profile_gene3", "nearest_gene3", "profile_confirmed_baseline")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stable(text):
    return int(hashlib.sha256(text.encode()).hexdigest()[:12], 16)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows, fields=None):
    with Path(path).open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def taxonomy(name):
    parts = name.strip("'").removeprefix("_R_").split("_")
    if len(parts) < 4 or not parts[3] or not parts[3][0].islower():
        return None
    return parts[0], parts[1], (parts[2] + "_" + parts[3]).casefold()


def tree_distances(path):
    """All terminal patristic distances in quadratic work, without changing tree."""
    tree = Phylo.read(path, "newick")
    tips = tree.get_terminals()
    names = [tip.name for tip in tips]
    if len(set(names)) != len(names):
        raise ValueError(f"Duplicate terminal names: {path}")
    index = {id(t): i for i, t in enumerate(tips)}
    depths = tree.depths()
    d = np.zeros((len(tips), len(tips)))
    leaf_depth = np.array([depths[t] for t in tips])
    descendants = {}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            descendants[id(node)] = np.array([index[id(node)]])
        else:
            children = [descendants[id(c)] for c in node.clades]
            for i, a in enumerate(children):
                for b in children[i + 1:]:
                    block = leaf_depth[a, None] + leaf_depth[None, b] - 2 * depths[node]
                    d[np.ix_(a, b)] = block
                    d[np.ix_(b, a)] = block.T
            descendants[id(node)] = np.concatenate(children)
    if (d < -1e-9).any() or not np.isfinite(d).all():
        raise ValueError(f"Invalid branch distances: {path}")
    return names, np.maximum(d, 0)


def reference_block(path, project):
    names, matrix = tree_distances(path)
    refs, taxa, queries = [], [], []
    for i, name in enumerate(names):
        if project in name:
            queries.append(i)
        else:
            tax = taxonomy(name)
            if tax:
                refs.append(i)
                taxa.append(tax)
    if len(refs) < 4:
        raise ValueError(f"Too few reference leaves: {path}")
    # One accession per species, selected deterministically, for equal representation.
    chosen = {}
    for i, tax in enumerate(taxa):
        chosen.setdefault(tax[2], i)
    unique = list(chosen.values())
    refs = np.array(refs)[unique]
    taxa = np.array(taxa)[unique]
    return (matrix[np.ix_(refs, refs)], taxa,
            matrix[np.ix_(queries, refs)], [names[i] for i in queries])


def profiles(q, reference, taxa, allowed_orders=None):
    """Full distance-shape features for every candidate, with no label-specific fit.

    Centering the entire query-to-anchor fingerprint removes an additive terminal
    branch offset. Family means (rather than sums) avoid reference-count inflation.
    Only reference distances set the common scale. No zero self-hit is permitted
    when this function is used to construct training pseudo-queries.
    """
    families = sorted(set(taxa[:, 1]))
    groups = {f: np.flatnonzero(taxa[:, 1] == f) for f in families}
    cross = taxa[:, 1, None] != taxa[None, :, 1]
    positive = reference[cross & (reference > 0)]
    scale = float(np.median(positive)) if len(positive) else 1.0
    scale = max(scale, 1e-8)
    near = {f: float(q[g].min()) for f, g in groups.items()}
    parents = {f: str(taxa[groups[f][0], 0]) for f in families}
    features, output = [], []
    for family in families:
        if allowed_orders is not None and parents[family] not in allowed_orders:
            continue
        ids = groups[family]
        others = np.flatnonzero(taxa[:, 1] != family)
        if not len(others):
            continue
        siblings = [near[f] for f in families if f != family and parents[f] == parents[family]]
        competitor = min(siblings) if siblings else min(near[f] for f in families if f != family)
        proto = reference[np.ix_(ids, others)].mean(axis=0)
        query = q[others]
        delta = query - proto
        offset = float(np.median(delta))
        rmse = float(np.sqrt(np.mean((delta - offset) ** 2))) / scale
        denominator = float(np.linalg.norm(query - query.mean()) * np.linalg.norm(proto - proto.mean()))
        correlation = float(np.dot(query - query.mean(), proto - proto.mean()) / denominator) if denominator > 1e-12 else 0
        features.append([near[family] / scale, float(q[ids].mean()) / scale,
                         (competitor - near[family]) / max(competitor + near[family], 1e-8),
                         rmse, correlation, offset / scale,
                         sum(v < near[family] for v in near.values()) / len(near)])
        output.append((family, parents[family]))
    return np.array(features).reshape(-1, len(FEATURES)), output


def choose_training_trees(source_root):
    selection = json.loads((Path(source_root) / "selection.json").read_text())["selected"]
    chosen = {}
    # One tree per gene, deterministic across all 100 contexts, independent of truth.
    for row in selection:
        stage = Path(row["run"]) / "02_family"
        valid = json.loads((stage / "family_tree_validation.json").read_text())["valid_trees"]
        for name in valid:
            gene = name.split("_")[0]
            key = stable(row["project"] + "/" + name)
            if gene not in chosen or key < chosen[gene][0]:
                chosen[gene] = (key, row["project"], stage / "valid_phylo" / name)
    return [(gene, *chosen[gene][1:]) for gene in sorted(chosen)]


def train(source_root, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    rows, audit = [], []
    training_trees = choose_training_trees(source_root)
    for gene, project, path in training_trees:
        before = digest(path)
        d, taxa, _, _ = reference_block(path, project)
        ids = sorted(range(len(taxa)), key=lambda i: stable(gene + taxa[i, 2]))[:64]
        for i in ids:
            keep = taxa[:, 2] != taxa[i, 2]
            if not (taxa[keep, 1] == taxa[i, 1]).any():
                continue
            x, labels = profiles(d[i, keep], d[np.ix_(keep, keep)], taxa[keep])
            order = sorted(range(len(labels)), key=lambda j: x[j, 0])
            selected = [j for j in order if labels[j][0] == taxa[i, 1]]
            selected += [j for j in order if labels[j][0] != taxa[i, 1]][:4]
            partition = "calibration" if stable(taxa[i, 2]) % 5 == 0 else "train"
            for j in selected:
                rows.append({"x": x[j].tolist(), "y": int(labels[j][0] == taxa[i, 1]),
                             "species": str(taxa[i, 2]), "family": str(taxa[i, 1]),
                             "partition": partition})
        assert digest(path) == before
        audit.append({"gene": gene, "project": project, "path": str(path), "sha256": before})
    fit = [r for r in rows if r["partition"] == "train"]
    calibration = [r for r in rows if r["partition"] == "calibration"]
    if len(fit) < 100 or len(calibration) < 30:
        raise ValueError("Insufficient independent reference pseudo-queries")
    # Equal total influence per known reference species, then class-balanced loss.
    counts = Counter(r["species"] for r in fit)
    weights = np.array([1 / counts[r["species"]] for r in fit])
    weights *= len(weights) / weights.sum()
    model = make_pipeline(StandardScaler(), LogisticRegression(C=0.1, class_weight="balanced", max_iter=2000))
    model.fit([r["x"] for r in fit], [r["y"] for r in fit], logisticregression__sample_weight=weights)
    # Unknown reference controls remove the entire pseudo-query family, not just
    # its best hit, before scoring the most convincing remaining family.
    null_scores, known_scores = [], []
    cal_species = set()
    for gene, project, path in training_trees:
        d, taxa, _, _ = reference_block(path, project)
        ids = sorted(range(len(taxa)), key=lambda i: stable(gene + taxa[i, 2]))[:64]
        for i in ids:
            if stable(taxa[i, 2]) % 5 != 0:
                continue
            keep = taxa[:, 2] != taxa[i, 2]
            if not (taxa[keep, 1] == taxa[i, 1]).any():
                continue
            x, labels = profiles(d[i, keep], d[np.ix_(keep, keep)], taxa[keep])
            p = model.predict_proba(x)[:, 1]
            # The binary model was trained on nearby hard negatives. Use it as
            # a match verifier, never extrapolate its score to rank far families.
            winner = int(x[:, 0].argmin())
            known_scores.append({"species": str(taxa[i, 2]), "score": float(p[winner]),
                                 "correct": labels[winner][0] == taxa[i, 1]})
            absent = taxa[:, 1] != taxa[i, 1]
            if absent.sum() < 4 or len(set(taxa[absent, 1])) < 2:
                continue
            nx, _ = profiles(d[i, absent], d[np.ix_(absent, absent)], taxa[absent])
            if len(nx):
                nearest = int(nx[:, 0].argmin())
                null_scores.append({"species": str(taxa[i, 2]), "score": float(model.predict_proba(nx)[nearest, 1])})
                cal_species.add(str(taxa[i, 2]))
    # Equal total influence per species. Select a threshold on held-out reference
    # queries only, balancing correct known acceptance against family-absent rejection.
    if len(cal_species) < 10:
        raise ValueError("Insufficient family-absent reference calibration species")
    def species_mean(rows, predicate):
        groups = defaultdict(list)
        for row in rows:
            groups[row["species"]].append(float(predicate(row)))
        return float(np.mean([np.mean(v) for v in groups.values()]))
    curve = []
    for cutoff in np.linspace(0, 1, 101):
        correct_accept = species_mean(known_scores, lambda r: r["correct"] and r["score"] >= cutoff)
        absent_reject = species_mean(null_scores, lambda r: r["score"] < cutoff)
        curve.append({"threshold": float(cutoff), "balanced_reference_utility": (correct_accept + absent_reject) / 2,
                      "known_correct_acceptance": correct_accept, "family_absent_rejection": absent_reject})
    chosen = max(curve, key=lambda r: (r["balanced_reference_utility"], r["threshold"]))
    threshold = chosen["threshold"]
    bundle = {"model": model, "threshold": threshold, "features": FEATURES, "schema": "full-profile-transfer-v1"}
    joblib.dump(bundle, output / "model.joblib", compress=3)
    write_json(output / "training.json", {"schema": bundle["schema"], "truth_used": False,
        "source_root": str(source_root), "training_trees": audit, "training_rows": len(fit),
        "training_species": sorted(set(r["species"] for r in fit)), "calibration_species": sorted(cal_species),
        "calibration_rows": len(calibration), "threshold": threshold,
        "null_calibration_rule": "species-balanced known correct acceptance / family-absent rejection utility, 101 fixed cutoffs; nearest-candidate verification only",
        "calibration_choice": chosen, "calibration_curve": curve,
        "known_top1_accuracy": float(np.mean([r["correct"] for r in known_scores])),
        "known_accepted_accuracy": float(np.mean([r["correct"] for r in known_scores if r["score"] >= threshold]))
            if any(r["score"] >= threshold for r in known_scores) else None,
        "known_acceptance": float(np.mean([r["score"] >= threshold for r in known_scores])),
        "features": FEATURES, "code_sha256": digest(__file__)})
    print(json.dumps({"threshold": threshold, "fit_rows": len(fit), "calibration_species": len(cal_species)}), flush=True)


def predict_sample(source_root, model_dir, output_root, index):
    source_root, model_dir, output_root = map(Path, (source_root, model_dir, output_root))
    row = json.loads((source_root / "selection.json").read_text())["selected"][index]
    source = Path(row["run"])
    output = output_root / "samples" / row["project"]
    output.mkdir(parents=True, exist_ok=True)
    if (output / "SUCCESS").exists():
        return
    bundle = joblib.load(model_dir / "model.joblib")
    model, threshold = bundle["model"], bundle["threshold"]
    order_file = source / "01_order/order_candidates.txt"
    allowed = set(order_file.read_text().splitlines())
    stage = source / "02_family"
    validation = stage / "family_tree_validation.json"
    tree_files = json.loads(validation.read_text())["valid_trees"]
    inputs = {str(p): digest(p) for p in [order_file, validation, stage / "family_panel_report.json",
                                         stage / "family_candidates.txt", stage / "native_z_candidates.txt"]}
    candidates = set((stage / "family_candidates.txt").read_text().splitlines())
    native = set((stage / "native_z_candidates.txt").read_text().splitlines())
    genes, full_genes, nearest_genes, scores, parents = set(), defaultdict(set), defaultdict(set), defaultdict(list), {}
    exon_rows = []
    feature_rows = []
    for name in sorted(tree_files):
        path = stage / "valid_phylo" / name
        inputs[str(path)] = digest(path)
        d, taxa, queries, query_names = reference_block(path, row["project"])
        gene = name.split("_")[0]
        genes.add(gene)
        for q, query_name in zip(queries, query_names):
            x, labels = profiles(q, d, taxa, allowed)
            if not len(x):
                continue
            p = model.predict_proba(x)[:, 1]
            nearest = int(x[:, 0].argmin())
            winner = nearest
            family, parent = labels[winner]
            parents.update(dict(labels))
            accepted = bool(p[winner] >= threshold)
            if accepted:
                full_genes[family].add(gene)
            nearest_genes[labels[nearest][0]].add(gene)
            scores[family].append(float(p[winner]))
            exon_rows.append({"gene": gene, "exon": name, "query": query_name,
                              "family": family, "parent_order": parent, "score": float(p[winner]),
                              "accepted": accepted, "nearest_family": labels[nearest][0]})
            for label, value, features in zip(labels, p, x):
                feature_rows.append({"gene": gene, "exon": name, "query": query_name,
                                     "family": label[0], "score": float(value),
                                     **dict(zip(FEATURES, features.tolist()))})
    if len(genes) < 10:
        raise ValueError(f"{row['project']}: fewer than ten valid genes")
    predicted = {"profile_gene3": {f for f, g in full_genes.items() if len(g) >= 3},
                 "nearest_gene3": {f for f, g in nearest_genes.items() if len(g) >= 3},
                 "profile_confirmed_baseline": {f for f in candidates if len(full_genes[f]) >= 2}}
    for method, selected in predicted.items():
        (output / (method + ".txt")).write_text("".join(f + "\n" for f in sorted(selected)))
    for name, selected in [("baseline_confirmed", candidates), ("baseline_z09", native)]:
        (output / (name + ".txt")).write_text("".join(f + "\n" for f in sorted(selected)))
    evidence = []
    for family in sorted(parents):
        evidence.append({"family": family, "parent_order": parents[family],
            "profile_gene_votes": len(full_genes[family]), "nearest_gene_votes": len(nearest_genes[family]),
            "valid_genes": len(genes), "max_match_score": max(scores[family], default=0),
            **{m: family in predicted[m] for m in METHODS}, "baseline_confirmed": family in candidates})
    write_csv(output / "family_evidence.csv", evidence)
    write_csv(output / "exon_predictions.csv", exon_rows,
              ["gene", "exon", "query", "family", "parent_order", "score", "accepted", "nearest_family"])
    write_csv(output / "profile_features.csv", feature_rows,
              ["gene", "exon", "query", "family", "score"] + FEATURES)
    for path, expected in inputs.items():
        if digest(path) != expected:
            raise RuntimeError(f"Reused input changed during prediction: {path}")
    write_json(output / "audit.json", {"sample": row["sample"], "project": row["project"],
        "previously_seen": row["previously_seen"], "input_sha256": inputs,
        "orders": sorted(allowed), "valid_trees": len(tree_files), "valid_genes": len(genes),
        "model_sha256": digest(model_dir / "model.joblib"), "code_sha256": digest(__file__),
        "threshold": threshold, "truth_used": False, "candidate_count_cap": None,
        "reuse": "assembly, order predictions, references, family topology and branch lengths; no recomputation"})
    (output / "SUCCESS").write_text("complete\n")
    print(row["project"], {m: len(v) for m, v in predicted.items()}, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    t = sub.add_parser("train")
    t.add_argument("--source-root", required=True)
    t.add_argument("--output", required=True)
    t = sub.add_parser("predict")
    t.add_argument("--source-root", required=True)
    t.add_argument("--model-dir", required=True)
    t.add_argument("--output-root", required=True)
    t.add_argument("--index", type=int, required=True)
    args = p.parse_args()
    with threadpool_limits(limits=2):
        if args.command == "train":
            train(args.source_root, args.output)
        else:
            predict_sample(args.source_root, args.model_dir, args.output_root, args.index)


if __name__ == "__main__":
    main()
