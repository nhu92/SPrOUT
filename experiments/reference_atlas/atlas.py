"""Fixed reference atlas with supervised panel selection and pattern baselines.

Run `python -m experiments.reference_atlas.atlas --help` from the repository root.
Distances are retained in full dimension for prediction. PCoA is fitted to anchors
only. Species used as anchors, training queries, validation queries and test queries
are disjoint, including all their loci and accessions.
"""
import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
from collections import Counter
from importlib.metadata import version
from pathlib import Path

import joblib
import numpy as np
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

RANKS = {"order": 0, "family": 1, "genus": 2, "species": 3}


def taxon_id(header):
    """SPrOUT's Order_Family_Genus_species convention, dropping source suffixes."""
    parts = header.split()[0].lstrip(">").split("_")
    if len(parts) < 4 or any(not p for p in parts[:4]) or not parts[3][0].islower():
        raise ValueError(f"Expected Order_Family_Genus_species: {header}")
    return "_".join(parts[:4])


def read_fasta(path, skip_invalid=False, diagnostics=None):
    records, name, seq = {}, None, []
    for line in Path(path).read_text(encoding="utf-8").splitlines() + [">"]:
        if line.startswith(">"):
            if name is not None:
                s = "".join(seq).upper()
                try:
                    key = taxon_id(name)
                except ValueError:
                    if not skip_invalid:
                        raise
                    if diagnostics is not None:
                        diagnostics.append({"file": Path(path).name, "header": name,
                                            "reason": "missing/invalid species identity"})
                    name, seq = line[1:], []
                    continue
                # One representative per species/locus; never split accessions.
                if key not in records or sum(c in "ACGT" for c in s) > sum(
                    c in "ACGT" for c in records[key]
                ):
                    records[key] = s
            name, seq = line[1:], []
        elif line.strip():
            if name is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            seq.append(line.strip())
    if not records or len({len(s) for s in records.values()}) != 1:
        raise ValueError(f"Empty or unaligned FASTA: {path}")
    return records


def p_distances(sequences, min_overlap=80):
    """Pairwise deletion p-distance. Gaps/ambiguity never count as matches."""
    if not sequences or len({len(s) for s in sequences}) != 1:
        raise ValueError("Sequences must have the same alignment length")
    x = np.array([list(s) for s in sequences])
    valid = np.isin(x, list("ACGT")).astype(np.float32)
    overlap = valid @ valid.T
    matches = np.zeros_like(overlap)
    for base in "ACGT":
        b = (x == base).astype(np.float32)
        matches += b @ b.T
    distance = np.divide(overlap - matches, overlap,
                         out=np.full_like(overlap, np.nan), where=overlap >= min_overlap)
    return distance, overlap / x.shape[1]


def prepare(alignments, output, min_overlap=80, max_loci=None, gene_list=None, provenance=None,
            skip_invalid=False):
    files = sorted(Path(alignments).glob("*.fasta"))
    if gene_list is not None:
        requested = {Path(name).stem for name in Path(gene_list).read_text(encoding="utf-8").split()}
        files = [p for p in files if p.stem in requested]
        missing = requested - {p.stem for p in files}
        if missing:
            raise ValueError(f"Requested loci missing from reference: {sorted(missing)}")
    if max_loci is not None:
        files = files[:max_loci]
    if not files:
        raise ValueError("No *.fasta alignments found")
    excluded = []
    records = [read_fasta(p, skip_invalid, excluded) for p in files]
    names = sorted(set().union(*(r.keys() for r in records)))
    index = {s: i for i, s in enumerate(names)}
    distances = np.full((len(files), len(names), len(names)), np.nan, dtype=np.float32)
    coverage = np.zeros_like(distances)
    source_files = []
    for k, (p, r) in enumerate(zip(files, records)):
        ids = [index[s] for s in r]
        d, c = p_distances(list(r.values()), min_overlap)
        distances[k][np.ix_(ids, ids)] = d
        coverage[k][np.ix_(ids, ids)] = c
        source_files.append({"file": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()})
        print(f"Prepared {p.stem}: {len(r)} species", flush=True)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, distances=distances, coverage=coverage,
                        names=np.array(names), loci=np.array([p.stem for p in files]),
                        metadata=json.dumps({"distance_kind": "pairwise-deletion-p",
                                             "min_overlap": min_overlap, "files": source_files,
                                             "reference_builder": provenance,
                                             "excluded_invalid_identities": excluded}))


def prepare_refbuilder(args):
    """Use MrKnowItAll's real installer/corrections, not copied taxonomy rules."""
    source = Path(args.refbuilder_source).resolve()
    if not (source / "src/sprout_refbuilder/backbone.py").is_file():
        raise ValueError("--refbuilder-source must be a MrKnowItAll checkout")
    destination = Path(args.backbone_dir)
    if destination.exists() and any(destination.iterdir()) and not args.reuse_installed:
        raise FileExistsError("Use a new empty --backbone-dir for a reproducible import")
    sys.path.insert(0, str(source / "src"))
    from sprout_refbuilder.backbone import install_backbone_archive
    if args.reuse_installed:
        install = json.loads((destination / "backbone_install.json").read_text(encoding="utf-8"))
        archive_digest = hashlib.sha256(Path(args.archive).read_bytes()).hexdigest()
        corrections_digest = hashlib.sha256((source / "src/sprout_refbuilder/resources/taxonomy_corrections.tsv").read_bytes()).hexdigest()
        if install["source_sha256"] != archive_digest or install["taxonomy_corrections_sha256"] != corrections_digest:
            raise ValueError("Existing install differs from requested archive/corrections")
    else:
        install = install_backbone_archive(Path(args.archive), destination)
    revision = subprocess.run(["git", "-C", str(source), "rev-parse", "HEAD"],
                              check=True, capture_output=True, text=True).stdout.strip()
    corrections = source / "src/sprout_refbuilder/resources/taxonomy_corrections.tsv"
    provenance = {"repository": "https://github.com/nhu92/MrKnowItAll", "revision": revision,
                  "installer_sha256": hashlib.sha256((source / "src/sprout_refbuilder/backbone.py").read_bytes()).hexdigest(),
                  "corrections_sha256": hashlib.sha256(corrections.read_bytes()).hexdigest(),
                  "installation": install}
    prepare(destination, args.output, args.min_overlap, None, args.gene_list, provenance,
            skip_invalid=True)


def read_square(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        names = next(reader)[1:]
        rows = list(reader)
    row_names = [r[0] for r in rows]
    if len(set(names)) != len(names) or set(row_names) != set(names) or len(row_names) != len(names):
        raise ValueError(f"Expected unique identical row/column labels in {path}")
    row_map = {r[0]: r[1:] for r in rows}
    matrix = np.array([row_map[n] for n in names], dtype=float)
    if matrix.shape != (len(names), len(names)) or not np.allclose(matrix, matrix.T, equal_nan=True):
        raise ValueError(f"Expected a symmetric square matrix in {path}")
    if np.isinf(matrix).any() or (matrix[np.isfinite(matrix)] < 0).any() or (matrix == 999).any():
        raise ValueError(f"Invalid/flagged distances in {path}")
    return names, matrix


def prepare_matrices(args):
    """Reference-only raw stage-03 matrices plus matched site-overlap fractions."""
    base = Path(args.manifest).parent
    with Path(args.manifest).open(newline="", encoding="utf-8") as f:
        entries = list(csv.DictReader(f))
    if not entries or len({e["locus"] for e in entries}) != len(entries):
        raise ValueError("Manifest must contain distinct nonempty loci")
    matrices, coverages, identifiers, sources = [], [], [], []
    for entry in entries:
        names, d = read_square(base / entry["path"])
        cn, c = read_square(base / entry["coverage_path"])
        if names != cn or (c[np.isfinite(c)] > 1).any() or not np.isfinite(c).all():
            raise ValueError("Coverage must match distance labels and contain fractions in 0..1")
        ids = [taxon_id(n) for n in names]
        if len(set(ids)) != len(ids):
            raise ValueError("Use one accession per species/locus when preparing reference matrices")
        if not np.allclose(np.diag(d)[np.isfinite(np.diag(d))], 0):
            raise ValueError("Distance matrix diagonal must be zero")
        matrices.append(d)
        coverages.append(c)
        identifiers.append(ids)
        sources.append({"locus": entry["locus"], "distance_sha256": hashlib.sha256(
            (base / entry["path"]).read_bytes()).hexdigest(), "coverage_sha256": hashlib.sha256(
            (base / entry["coverage_path"]).read_bytes()).hexdigest()})
    names = sorted(set().union(*map(set, identifiers)))
    lookup = {n: i for i, n in enumerate(names)}
    distances = np.full((len(entries), len(names), len(names)), np.nan, dtype=np.float32)
    coverage = np.zeros_like(distances)
    for k, ids in enumerate(identifiers):
        idx = [lookup[n] for n in ids]
        distances[k][np.ix_(idx, idx)] = matrices[k]
        coverage[k][np.ix_(idx, idx)] = coverages[k]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, distances=distances, coverage=coverage, names=np.array(names),
                        loci=np.array([e["locus"] for e in entries]),
                        metadata=json.dumps({"distance_kind": args.distance_kind, "files": sources}))


def split_species(names, rank, seed):
    """Four partitions within taxa with >=4 species; smaller taxa are unknown."""
    identities = ["_".join(n.split("_")[2:4]).casefold() for n in names]
    if len(set(identities)) != len(identities):
        raise ValueError("A binomial occurs under conflicting taxonomy; curate identities before splitting")
    labels = np.array([s.split("_")[RANKS[rank]] if rank != "species" else s for s in names])
    rng = np.random.default_rng(seed)
    split = {s: [] for s in ("anchor", "train", "validation", "test", "unknown")}
    for label in sorted(set(labels)):
        ids = rng.permutation(np.flatnonzero(labels == label))
        if len(ids) < 4:
            split["unknown"].extend(ids.tolist())
            continue
        n_anchor = max(1, len(ids) // 4)
        n_eval = max(1, len(ids) // 6)
        split["anchor"].extend(ids[:n_anchor].tolist())
        split["validation"].extend(ids[n_anchor:n_anchor + n_eval].tolist())
        split["test"].extend(ids[n_anchor + n_eval:n_anchor + 2 * n_eval].tolist())
        split["train"].extend(ids[n_anchor + 2 * n_eval:].tolist())
    if len(set(labels[split["train"]])) < 2:
        raise ValueError(f"Rank {rank} needs >=2 taxa with >=4 distinct species each. "
                         "This reference cannot validate that rank on unseen species.")
    return labels, {k: np.array(v, dtype=int) for k, v in split.items()}


def normalize_distances(distances, anchors):
    """Each locus uses a positive median computed exclusively between anchors."""
    scales = []
    for matrix in distances:
        aa = matrix[np.ix_(anchors, anchors)]
        positive = aa[np.isfinite(aa) & (aa > 0)]
        scales.append(float(np.median(positive)) if len(positive) else 1.0)
    return distances / np.array(scales)[:, None, None], np.array(scales)


def aggregate_profiles(distances):
    count = np.isfinite(distances).sum(axis=0)
    return np.divide(np.nansum(distances, axis=0), count,
                     out=np.full(count.shape, np.nan, dtype=float), where=count > 0)


def select_panel(distances, labels, split, budget, seed):
    """Supervised between/within class separation of distance-to-anchor features.

    Scores use training query species only. One anchor per known taxon is retained,
    then remaining slots follow the discriminative score. No test/validation labels
    or distances affect selection.
    """
    candidates, queries = split["anchor"], split["train"]
    x = aggregate_profiles(distances[:, queries][:, :, candidates])
    median = np.array([np.median(col[np.isfinite(col)]) if np.isfinite(col).any() else 1
                       for col in x.T])
    x = np.where(np.isfinite(x), x, median)
    y = labels[queries]
    means = np.array([x[y == c].mean(axis=0) for c in sorted(set(y))])
    within = np.mean([np.mean((x[y == c] - means[i]) ** 2, axis=0)
                      for i, c in enumerate(sorted(set(y)))], axis=0)
    scores = means.var(axis=0) / (within + 0.01)
    classes = sorted(set(labels[candidates]))
    if budget < len(classes):
        raise ValueError(f"Panel budget {budget} must cover {len(classes)} known taxa")
    budget = min(budget, len(candidates))
    chosen = [max(np.flatnonzero(labels[candidates] == c), key=lambda j: scores[j])
              for c in classes]
    chosen += [j for j in np.argsort(-scores) if j not in chosen][:budget - len(chosen)]
    rng = np.random.default_rng(seed)
    random = [int(rng.choice(np.flatnonzero(labels[candidates] == c))) for c in classes]
    random += rng.permutation([j for j in range(len(candidates)) if j not in random])[
        :budget - len(random)].tolist()
    return candidates[chosen], candidates[random], scores


def fit_pcoa(distance):
    """Classical PCoA. Report discarded negative eigenvalue mass explicitly."""
    d = np.asarray(distance, dtype=float)
    if d.ndim != 2 or d.shape[0] != d.shape[1] or not np.isfinite(d).all():
        raise ValueError("Atlas needs a finite square anchor distance matrix")
    if (d < 0).any() or not np.allclose(d, d.T, atol=1e-6):
        raise ValueError("Atlas distances must be nonnegative and symmetric")
    np.fill_diagonal(d, 0)
    d2 = d ** 2
    row, grand = d2.mean(axis=1), float(d2.mean())
    b = -0.5 * (d2 - row[:, None] - row[None, :] + grand)
    values, vectors = np.linalg.eigh(b)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    keep = values > max(1e-12, abs(values).max() * 1e-10)
    if keep.sum() < 2:
        raise ValueError("Anchor matrix has fewer than two positive PCoA dimensions")
    basis = vectors[:, keep][:, :2]
    eig = values[keep][:2]
    # Deterministic sign convention for portable image coordinates.
    for j in range(2):
        if basis[np.argmax(abs(basis[:, j])), j] < 0:
            basis[:, j] *= -1
    return {"xy": basis * np.sqrt(eig), "basis": basis, "eigenvalues": eig,
            "row_mean_d2": row, "grand_mean_d2": grand,
            "positive_mass_2d": float(eig.sum() / values[keep].sum()),
            "negative_mass_fraction": float(abs(values[values < 0]).sum() /
                                             max(abs(values).sum(), 1e-12))}


def complete_anchor_distances(distance, method="error"):
    """Optional layout-only completion from the observed anchor graph.

    Observed distances remain unchanged. Query features are never graph-imputed.
    Report the imputed pair count; disconnection is a hard failure.
    """
    d = np.array(distance, dtype=float, copy=True)
    np.fill_diagonal(d, 0)
    missing = ~np.isfinite(d)
    count = int(np.triu(missing, 1).sum())
    if count and method == "error":
        raise ValueError(f"Atlas lacks {count} anchor pairs. Use --atlas-missing shortest-path "
                         "for explicit layout-only graph completion, or improve reference coverage")
    if count:
        paths = np.where(missing, np.inf, d)
        for k in range(len(paths)):
            paths = np.minimum(paths, paths[:, k, None] + paths[None, k, :])
        if not np.isfinite(paths).all():
            raise ValueError("Anchor distance graph is disconnected; cannot build a shared atlas")
        d[missing] = paths[missing]
    return d, count


def project_pcoa(distance_to_anchor, atlas):
    d2 = distance_to_anchor ** 2
    b = -0.5 * (d2 - d2.mean(axis=1, keepdims=True) - atlas["row_mean_d2"] +
                atlas["grand_mean_d2"])
    return b @ atlas["basis"] / np.sqrt(atlas["eigenvalues"])


def raster_kernel(xy, size=12, bandwidth=0.12):
    center, span = xy.min(axis=0), np.ptp(xy, axis=0)
    unit = (xy - center) / np.maximum(span, 1e-12)
    grid = np.stack(np.meshgrid(np.linspace(0, 1, size), np.linspace(0, 1, size)), -1)
    kernel = np.exp(-((unit[:, None, :] - grid.reshape(-1, 2)) ** 2).sum(axis=2) /
                    (2 * bandwidth ** 2))
    return kernel


def feature_rows(distance, coverage):
    available = np.isfinite(distance)
    # Zero affinity is used only alongside explicit availability/overlap channels.
    affinity = np.where(available, np.exp(-np.nan_to_num(distance)), 0)
    return np.stack([affinity, np.where(available, coverage, 0), available.astype(float)], axis=1)


def rasterize(channels, kernel):
    # Average splatting keeps anchor density from dominating image intensities.
    return (channels @ kernel / np.maximum(kernel.sum(axis=0), 1e-12)).reshape(len(channels), -1)


def classifier():
    return make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=1500,
                                                              class_weight="balanced"))


def fuse(probabilities, ids, targets):
    """Equal locus mean per query species/contig; never a mixture-wide softmax."""
    result = []
    for target in targets:
        rows = probabilities[ids == target]
        result.append(rows.mean(axis=0) if len(rows) else np.full(probabilities.shape[1], np.nan))
    return np.array(result)


def choose_threshold(probabilities, truth, classes, target_precision=0.9):
    """Validation-only selective classification cutoff, without calibration claims."""
    confidence = probabilities.max(axis=1)
    correct = classes[probabilities.argmax(axis=1)] == truth
    candidates = []
    for cutoff in np.unique(confidence):
        selected = confidence >= cutoff
        # A single lucky validation query cannot establish a usable cutoff.
        if selected.sum() >= 5 and correct[selected].mean() >= target_precision:
            candidates.append(float(cutoff))
    return min(candidates) if candidates else 1.01


def evaluate(probabilities, truth, classes, threshold=None):
    finite = np.isfinite(probabilities).all(axis=1)
    p, y = probabilities[finite], truth[finite]
    if not len(y):
        return {"n": 0, "no_evidence": int((~finite).sum())}
    predicted = classes[p.argmax(axis=1)]
    result = {"n": len(y), "no_evidence": int((~finite).sum()),
              "accuracy": float(accuracy_score(y, predicted)),
              "macro_f1": float(f1_score(y, predicted, labels=classes,
                                         average="macro", zero_division=0))}
    if threshold is not None:
        selected = p.max(axis=1) >= threshold
        result.update({"threshold": threshold, "accepted": int(selected.sum()),
                       "coverage": float(selected.mean()),
                       "selective_accuracy": float((predicted[selected] == y[selected]).mean())
                       if selected.any() else None})
    return result


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                          encoding="utf-8")


def benchmark(args):
    with np.load(args.dataset, allow_pickle=False) as archive:
        data = {key: archive[key] for key in archive.files}
    if getattr(args, "within_order", None):
        selected = np.array([n.split("_")[0] == args.within_order for n in data["names"]])
        if not selected.any():
            raise ValueError(f"No reference species in order {args.within_order}")
        data["names"] = data["names"][selected]
        for key in ("distances", "coverage"):
            data[key] = data[key][:, selected][:, :, selected]
    names, loci = data["names"], data["loci"]
    labels, split = split_species(names, args.rank, args.seed)
    normalized, scales = normalize_distances(data["distances"], split["anchor"])
    panel, random_panel, scores = select_panel(normalized, labels, split, args.panel_size, args.seed)
    # All candidate anchors use the same training-only normalization for fair ablations.
    aa = aggregate_profiles(normalized[:, panel][:, :, panel])
    aa, imputed_pairs = complete_anchor_distances(aa, getattr(args, "atlas_missing", "error"))
    atlas = fit_pcoa(aa)
    kernel = raster_kernel(atlas["xy"], args.grid_size)
    all_ids = np.tile(np.arange(len(names)), len(loci))
    locus_ids = np.repeat(np.arange(len(loci)), len(names))
    panel_d = normalized[:, :, panel].reshape(-1, len(panel))
    panel_c = data["coverage"][:, :, panel].reshape(-1, len(panel))
    channels = feature_rows(panel_d, panel_c)
    usable = np.isfinite(panel_d).sum(axis=1) >= max(2, len(panel) // 4)
    train = np.isin(all_ids, split["train"]) & usable
    if set(labels[all_ids[train]]) != set(labels[split["train"]]):
        raise ValueError("At least one training class has no usable query loci")
    # Train-only fill for PCoA placement and its two-dimensional baselines.
    fill = np.array([np.median(c[np.isfinite(c)]) if np.isfinite(c).any() else 1.0
                     for c in panel_d[train].T])
    xy_queries = project_pcoa(np.where(np.isfinite(panel_d), panel_d, fill), atlas)
    random_d = normalized[:, :, random_panel].reshape(-1, len(random_panel))
    random_c = data["coverage"][:, :, random_panel].reshape(-1, len(random_panel))
    all_d = normalized[:, :, split["anchor"]].reshape(-1, len(split["anchor"]))
    all_c = data["coverage"][:, :, split["anchor"]].reshape(-1, len(split["anchor"]))
    features = {"distance_learned_panel": channels.reshape(len(channels), -1),
                "distance_random_panel": feature_rows(random_d, random_c).reshape(len(channels), -1),
                "distance_all_anchors": feature_rows(all_d, all_c).reshape(len(channels), -1),
                "raster_pattern": rasterize(channels, kernel), "pcoa2_logistic": xy_queries}
    classes = np.array(sorted(set(labels[split["train"]])))
    # Give species equal total influence even when their locus recovery differs.
    counts = Counter(all_ids[train])
    weights = np.array([1 / counts[i] for i in all_ids[train]])
    weights *= len(weights) / weights.sum()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    results, fitted, fused_predictions = {}, {}, {}
    for name, x in features.items():
        print(f"Fitting {name}: {train.sum()} locus profiles, {x.shape[1]} features", flush=True)
        model = classifier()
        model.fit(x[train], labels[all_ids[train]], logisticregression__sample_weight=weights)
        probabilities = model.predict_proba(x)
        probabilities[~usable] = np.nan
        valid_rows = usable
        val = fuse(probabilities[valid_rows], all_ids[valid_rows], split["validation"])
        if not np.isfinite(val).all():
            raise ValueError("Validation species without evidence; increase coverage or change inputs")
        threshold = choose_threshold(val, labels[split["validation"]], classes)
        test = fuse(probabilities[valid_rows], all_ids[valid_rows], split["test"])
        unknown = fuse(probabilities[valid_rows], all_ids[valid_rows], split["unknown"])
        result = {"validation": evaluate(val, labels[split["validation"]], classes),
                  "test": evaluate(test, labels[split["test"]], classes, threshold)}
        finite_unknown = unknown[np.isfinite(unknown).all(axis=1)]
        if len(finite_unknown):
            result["unknown"] = {"n": len(finite_unknown),
                "false_accept_rate": float((finite_unknown.max(axis=1) >= threshold).mean()),
                "known_vs_unknown_confidence_auroc": float(roc_auc_score(
                    [1] * len(test[np.isfinite(test).all(axis=1)]) + [0] * len(finite_unknown),
                    np.r_[test[np.isfinite(test).all(axis=1)].max(axis=1), finite_unknown.max(axis=1)]))}
        results[name], fitted[name] = result, model
        fused_predictions[name] = (test, threshold)
    # A fixed PCoA + K-means analogue of the abandoned approach. Cluster-to-label
    # mapping uses anchor labels only; repeated labels in separate clusters allowed.
    km = KMeans(n_clusters=len(classes), random_state=args.seed, n_init=10).fit(atlas["xy"])
    cluster_labels = {}
    for cluster in range(len(classes)):
        members = labels[panel][km.labels_ == cluster]
        cluster_labels[cluster] = Counter(members).most_common(1)[0][0] if len(members) else classes[0]
    guessed = np.array([cluster_labels[c] for c in km.predict(xy_queries)])
    kp = (guessed[:, None] == classes[None, :]).astype(float)
    kt = fuse(kp[usable], all_ids[usable], split["test"])
    results["pcoa2_kmeans"] = {"test": evaluate(kt, labels[split["test"]], classes)}
    # Nearest reference in the full normalized distance profile (same panel).
    nearest = np.argmin(np.where(np.isfinite(panel_d), panel_d, np.inf), axis=1)
    npred = (labels[panel][nearest, None] == classes[None, :]).astype(float)
    nt = fuse(npred[usable], all_ids[usable], split["test"])
    results["nearest_anchor"] = {"test": evaluate(nt, labels[split["test"]], classes)}
    report = {"schema": "sprout-reference-atlas-v1", "seed": args.seed, "rank": args.rank,
              "within_order": getattr(args, "within_order", None),
              "software": {"python": platform.python_version(), **{p: version(p) for p in
                           ("numpy", "scipy", "scikit-learn", "matplotlib", "joblib")}},
              "experiment_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "dataset_metadata": json.loads(str(data["metadata"])),
              "n_species": len(names), "n_loci": len(loci), "n_classes": len(classes),
              "panel_size": len(panel), "anchor_candidates": len(split["anchor"]),
              "split": {k: names[v].tolist() for k, v in split.items()},
              "atlas": {k: atlas[k] for k in ("positive_mass_2d", "negative_mass_fraction")},
              "atlas_completion": {"method": getattr(args, "atlas_missing", "error"),
                                   "imputed_pairs": imputed_pairs},
              "results": results,
              "limitations": ["Reference pseudo-queries, no raw-read or assembled-mixture validation",
                "Unknown taxa are too rare to train, not a random holdout of entire well-sampled clades",
                "Confidence is an uncalibrated model score; cutoff uses known validation queries only",
                "Loci are genes in supplied sample data, not independently identified exons",
                "No species-level generalization claim or DNA abundance estimate"]}
    write_json(out / "metrics.json", report)
    bundle = {"schema": report["schema"], "rank": args.rank, "panel": names[panel],
              "within_order": report["within_order"],
              "loci": loci, "scales": scales, "atlas": atlas, "kernel": kernel,
              "grid_size": args.grid_size, "classes": classes, "fill": fill,
              "models": {k: fitted[k] for k in ("distance_learned_panel", "raster_pattern")},
              "thresholds": {k: fused_predictions[k][1] for k in ("distance_learned_panel", "raster_pattern")},
              "distance_kind": report["dataset_metadata"]["distance_kind"]}
    joblib.dump(bundle, out / "model.joblib", compress=3)
    np.savez_compressed(out / "locus_patterns.npz", channels=rasterize(channels[usable], kernel).reshape(
        -1, 3, args.grid_size, args.grid_size), species=names[all_ids[usable]],
        split=np.array([next(k for k, v in split.items() if i in v) for i in all_ids[usable]]),
        loci=loci[locus_ids[usable]])
    with (out / "panel.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["species", "label", "x", "y", "discrimination_score"])
        candidate_lookup = {v: i for i, v in enumerate(split["anchor"])}
        for j, idx in enumerate(panel):
            writer.writerow([names[idx], labels[idx], *atlas["xy"][j], scores[candidate_lookup[idx]]])
    with (out / "predictions.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["species", "truth", "method", "prediction", "score", "accepted"])
        for name, (p, threshold) in fused_predictions.items():
            for idx, row in zip(split["test"], p):
                if np.isfinite(row).all():
                    writer.writerow([names[idx], labels[idx], name, classes[row.argmax()],
                                     row.max(), bool(row.max() >= threshold)])
    # Fixed, declared synthetic contig pool, never selected based on model success.
    mix_ids = np.array([next(i for i in split["test"] if labels[i] == c)
                        for c in classes[:3]])
    for name in ("distance_learned_panel", "raster_pattern"):
        probs = fitted[name].predict_proba(features[name])
        rows = usable & np.isin(all_ids, mix_ids)
        # One contribution per represented gene, avoiding raw contig-count inflation.
        evidence = np.array([probs[rows & (locus_ids == k)].mean(axis=0)
                             for k in range(len(loci)) if (rows & (locus_ids == k)).any()]).mean(axis=0)
        write_json(out / f"synthetic_mixture_{name}.json", {
            "kind": "in-silico pool of reference locus profiles, no assembly simulation",
            "species": names[mix_ids].tolist(), "truth": labels[mix_ids].tolist(),
            "class_evidence": dict(zip(classes.tolist(), evidence.tolist()))})
        save_surface(bundle, evidence, out / f"surface_{name}.png", name,
                     "Synthetic reference pool: " + ", ".join(labels[mix_ids]))
    save_comparison(results, out / "comparison.png", args.rank)
    save_patterns(channels, all_ids, locus_ids, split, names, loci, kernel, args.grid_size,
                  out / "patterns.png")
    print(json.dumps({k: v["test"] for k, v in results.items()}, indent=2))


def save_surface(bundle, evidence, path, method, subtitle):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xy = bundle["atlas"]["xy"]
    panel_labels = np.array([s.split("_")[RANKS[bundle["rank"]]] for s in bundle["panel"]])
    counts = Counter(panel_labels)
    # Share each class score among its anchors, preserving total evidence mass.
    z = np.array([evidence[np.flatnonzero(bundle["classes"] == c)[0]] / counts[c] for c in panel_labels])
    gx, gy = np.meshgrid(np.linspace(xy[:, 0].min(), xy[:, 0].max(), 65),
                         np.linspace(xy[:, 1].min(), xy[:, 1].max(), 65))
    grid = np.stack([gx, gy], axis=-1)
    bw = max(float(np.ptp(xy, axis=0).max()) * 0.10, 1e-8)
    k = np.exp(-((grid[:, :, None, :] - xy) ** 2).sum(axis=-1) / (2 * bw ** 2))
    zz = k @ z
    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(gx, gy, zz, cmap="viridis", alpha=0.9, linewidth=0)
    ax.scatter(xy[:, 0], xy[:, 1], z, c=z, cmap="viridis", edgecolor="black", s=18)
    for cls in bundle["classes"][np.argsort(-evidence)[:5]]:
        idx = np.flatnonzero(panel_labels == cls)[0]
        ax.text(xy[idx, 0], xy[idx, 1], z[idx] + 0.01, cls, fontsize=8)
    ax.set(xlabel="Reference PCoA 1", ylabel="Reference PCoA 2")
    fig.text(0.88, 0.52, "Relative evidence", rotation=90, va="center", fontsize=10)
    fig.suptitle("Reference atlas: learned evidence surface", fontsize=15)
    ax.set_title(method.replace("_", " "), fontsize=11)
    fig.text(0.5, 0.04, subtitle + "\nHeight is model support, not DNA abundance or a calibrated probability.",
             ha="center", fontsize=10)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def save_patterns(channels, species_ids, locus_ids, split, names, loci, kernel, size, path):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 3, figsize=(10, 9), layout="constrained")
    for row, species in enumerate(split["test"][:3]):
        candidates = np.flatnonzero((species_ids == species) & (channels[:, 2].sum(axis=1) > 0))
        if not len(candidates):
            continue
        idx = candidates[0]
        images = rasterize(channels[idx:idx + 1], kernel).reshape(3, size, size)
        for col, channel in enumerate(("Affinity", "Pairwise site overlap", "Availability mask")):
            axes[row, col].imshow(images[col], origin="lower", cmap="viridis", vmin=0, vmax=1)
            axes[row, col].set(xticks=[], yticks=[])
            if row == 0:
                axes[row, col].set_title(channel)
            if col == 0:
                parts = names[species].split("_")
                axes[row, col].set_ylabel(f"{parts[0]}\n{parts[2]} {parts[3]}\nLocus {loci[locus_ids[idx]]}")
    fig.suptitle("Actual pattern inputs on one frozen reference grid", fontsize=14)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def save_comparison(results, path, rank):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = list(results)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.barh(names, [results[k]["test"]["macro_f1"] for k in names], color="#217c85")
    ax.set(xlim=(0, 1.06), xlabel=f"Macro F1 on held-out species (known {rank} taxa)",
           title="Reference pseudo-query benchmark")
    for i, k in enumerate(names):
        ax.text(results[k]["test"]["macro_f1"] + 0.01, i,
                f'{results[k]["test"]["macro_f1"]:.3f}', va="center")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def read_query_matrix(path, panel, prefix):
    """Adapt raw square stage-03 CSVs, not filtered/flagged summary distances."""
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader)[1:]
        rows = list(reader)
    columns = {}
    panel_set = set(panel)
    for j, name in enumerate(header):
        if name.startswith(prefix):
            continue
        try:
            key = taxon_id(name)
        except ValueError:
            continue
        if key in panel_set:
            if key in columns:
                raise ValueError(f"Duplicate anchor {key} in {path}")
            columns[key] = j
    queries, values = [], []
    for row in rows:
        if row[0].startswith(prefix):
            queries.append(row[0])
            values.append([float(row[1 + columns[s]]) if s in columns else np.nan for s in panel])
    d = np.array(values, dtype=float).reshape(-1, len(panel))
    if np.isinf(d).any() or (d[np.isfinite(d)] < 0).any() or (d[np.isfinite(d)] == 999).any():
        raise ValueError(f"Invalid/flagged distances in {path}; use raw .matrix files")
    return queries, d


def predict_matrices(args):
    model = joblib.load(args.model)
    if args.distance_kind != model["distance_kind"]:
        raise ValueError("Prediction distance kind differs from training. Do not apply a p-distance "
                         "model to tree branch lengths; prepare a compatible training distance dataset.")
    locus_map = {s: i for i, s in enumerate(model["loci"])}
    probabilities, names, locus_names = [], [], []
    # Manifest makes gene/exon identity explicit; matching filenames heuristically
    # can silently join the wrong loci (e.g. 4691.1 versus 4691_exon_1).
    with Path(args.manifest).open(newline="", encoding="utf-8") as f:
        entries = list(csv.DictReader(f))
    seen = set()
    for entry in entries:
        locus = entry["locus"]
        if locus not in locus_map:
            raise ValueError(f"Untrained locus {locus}")
        path = Path(args.manifest).parent / entry["path"]
        query, d = read_query_matrix(path, model["panel"], args.query_prefix)
        if not len(query):
            continue
        d /= model["scales"][locus_map[locus]]
        usable = np.isfinite(d).sum(axis=1) >= max(2, len(model["panel"]) // 4)
        # Stage-03 matrices have no site-overlap metadata. Require the user to
        # supply compatible coverage, rather than treating unknown overlap as 100%.
        if "coverage_path" not in entry or not entry["coverage_path"]:
            raise ValueError("Manifest needs coverage_path: same query/anchor CSV layout, "
                             "with pairwise valid-site fractions (0..1)")
        cov_names, cov = read_query_matrix(Path(args.manifest).parent / entry["coverage_path"],
                                           model["panel"], args.query_prefix)
        if cov_names != query or (cov[np.isfinite(cov)] > 1).any():
            raise ValueError("Coverage rows must match distances, and fractions must be in 0..1")
        if not np.isfinite(cov[np.isfinite(d)]).all():
            raise ValueError("Every observed distance needs a coverage fraction")
        ch = feature_rows(d, np.nan_to_num(cov))
        x = ch.reshape(len(ch), -1) if args.method == "distance_learned_panel" else rasterize(ch, model["kernel"])
        p = model["models"][args.method].predict_proba(x)
        for q, row, ok in zip(query, p, usable):
            if (q, locus) in seen:
                raise ValueError(f"Duplicate query/locus {q}, {locus}")
            seen.add((q, locus))
            if ok:
                names.append(q)
                locus_names.append(locus)
                probabilities.append(row)
    if not names:
        raise ValueError("No usable query profiles found")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    threshold = model["thresholds"][args.method]
    with (out / "query_predictions.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["query", "locus", "rank", "prediction", "score", "accepted"])
        for q, locus, p in zip(names, locus_names, probabilities):
            w.writerow([q, locus, model["rank"], model["classes"][np.argmax(p)],
                        max(p), bool(max(p) >= threshold)])
    # Every locus contributes equally after averaging its contig profiles.
    evidence = np.mean([np.mean([p for p, loc in zip(probabilities, locus_names) if loc == locus], axis=0)
                        for locus in sorted(set(locus_names))], axis=0)
    write_json(out / "mixture_evidence.json", {"rank": model["rank"], "method": args.method,
        "height_meaning": "relative evidence, not abundance or calibrated presence probability",
        "class_evidence": dict(zip(model["classes"].tolist(), evidence.tolist()))})
    save_surface(model, evidence, out / "surface.png", args.method, "Query locus profiles")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare", help="Build distances from aligned reference FASTAs")
    p.add_argument("--alignments", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-overlap", type=int, default=80)
    p.add_argument("--max-loci", type=int)
    p.add_argument("--gene-list", help="Optional explicit set of reference loci")
    p = sub.add_parser("prepare-refbuilder", help="Install/correct a backbone using MrKnowItAll")
    p.add_argument("--refbuilder-source", required=True)
    p.add_argument("--archive", required=True)
    p.add_argument("--backbone-dir", required=True, help="New empty directory; no source checkout is changed")
    p.add_argument("--reuse-installed", action="store_true", help="Reuse an installation with matching source/correction hashes")
    p.add_argument("--gene-list", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-overlap", type=int, default=80)
    p = sub.add_parser("prepare-matrices", help="Build compatible reference dataset from raw CSV matrices")
    p.add_argument("--manifest", required=True, help="CSV: locus,path,coverage_path; reference-only squares")
    p.add_argument("--distance-kind", choices=["pairwise-deletion-p", "patristic"], required=True)
    p.add_argument("--output", required=True)
    p = sub.add_parser("benchmark", help="Four species partitions and paired method comparisons")
    p.add_argument("--dataset", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--rank", choices=RANKS, default="order")
    p.add_argument("--within-order", help="Optional parent-order scope for a local family/genus atlas")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--panel-size", type=int, default=48)
    p.add_argument("--grid-size", type=int, default=12)
    p.add_argument("--atlas-missing", choices=["error", "shortest-path"], default="error",
                   help="Explicit layout-only completion for missing anchor pairs")
    p = sub.add_parser("predict-matrices", help="Predict compatible query-to-anchor CSV profiles")
    p.add_argument("--model", required=True)
    p.add_argument("--manifest", required=True, help="CSV: locus,path,coverage_path")
    p.add_argument("--query-prefix", required=True)
    p.add_argument("--distance-kind", choices=["pairwise-deletion-p", "patristic"], required=True)
    p.add_argument("--method", choices=["distance_learned_panel", "raster_pattern"], default="distance_learned_panel")
    p.add_argument("--output", required=True)
    args = parser.parse_args()
    if getattr(args, "grid_size", 12) < 2 or getattr(args, "panel_size", 48) < 2:
        parser.error("grid-size and panel-size must be >=2")
    if getattr(args, "min_overlap", 80) < 1 or (getattr(args, "max_loci", None) is not None and args.max_loci < 1):
        parser.error("min-overlap and max-loci must be positive")
    with threadpool_limits(limits=2):
        if args.command == "prepare":
            prepare(args.alignments, args.output, args.min_overlap, args.max_loci, args.gene_list)
        elif args.command == "prepare-refbuilder":
            prepare_refbuilder(args)
        elif args.command == "prepare-matrices":
            prepare_matrices(args)
        elif args.command == "benchmark":
            benchmark(args)
        else:
            predict_matrices(args)


if __name__ == "__main__":
    main()
