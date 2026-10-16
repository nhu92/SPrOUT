"""Numerical, leakage and integration checks for the experimental atlas."""
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import joblib
import numpy as np
import pytest
from scipy.spatial.distance import cdist
from threadpoolctl import threadpool_limits

from experiments.reference_atlas.atlas import (
    aggregate_profiles, benchmark, choose_threshold, complete_anchor_distances, feature_rows, fit_pcoa,
    normalize_distances, p_distances, predict_matrices, prepare_matrices,
    project_pcoa, raster_kernel, rasterize, read_fasta, read_query_matrix,
    select_panel, split_species,
)


def test_pairwise_deletion_never_treats_gaps_or_ambiguity_as_matches():
    d, overlap = p_distances(["ACGT-N", "ATGTAA", "NNNNNN"], min_overlap=4)
    assert d[0, 1] == pytest.approx(0.25)
    assert overlap[0, 1] == pytest.approx(4 / 6)
    assert np.isnan(d[0, 2])
    assert np.isnan(d[2, 2])


def test_unaligned_fasta_fails_and_accessions_collapse(tmp_path):
    path = tmp_path / "gene.fasta"
    path.write_text(">O_F_G_s_sourceA\nA---\n>O_F_G_s_sourceB\nACGT\n")
    assert read_fasta(path) == {"O_F_G_s": "ACGT"}
    path.write_text(">O_F_G_s\nAAA\n>P_F_H_t\nAAAA\n")
    with pytest.raises(ValueError, match="unaligned"):
        read_fasta(path)


def test_missing_species_cannot_turn_gene_suffix_into_independent_species(tmp_path):
    path = tmp_path / "gene.fasta"
    path.write_text(">Asparagales_Asphodelaceae_Hemerocallis_1KP-4691\nACGT\n"
                    ">Asparagales_Asphodelaceae_Hemerocallis_sp_1KP-4691\nACGT\n")
    with pytest.raises(ValueError, match="species"):
        read_fasta(path)
    excluded = []
    records = read_fasta(path, skip_invalid=True, diagnostics=excluded)
    assert list(records) == ["Asparagales_Asphodelaceae_Hemerocallis_sp"]
    assert len(excluded) == 1


def synthetic_dataset():
    rng = np.random.default_rng(7)
    names, points = [], []
    for k, center in enumerate(([0, 0, 0], [3, 0, 0], [0, 3, 0])):
        for i in range(6):
            names.append(f"Order{k}_Family{k}_Genus{k}_species{i}")
            points.append(np.array(center) + rng.normal(0, 0.06, 3))
    names.append("Novel_Fam_Gen_species")
    points.append([5, 5, 5])
    d = cdist(points, points)
    distances = np.stack([d, d * 2, d * 0.5])
    coverage = np.ones_like(distances)
    return np.array(names), distances, coverage


def test_species_splits_are_disjoint_and_sparse_classes_are_unknown():
    names, _, _ = synthetic_dataset()
    labels, split = split_species(names, "order", 42)
    sets = [set(v) for v in split.values()]
    assert set.union(*sets) == set(range(len(names)))
    assert sum(map(len, sets)) == len(set.union(*sets))
    assert names[split["unknown"]].tolist() == ["Novel_Fam_Gen_species"]
    with pytest.raises(ValueError, match="cannot validate"):
        split_species(names, "species", 42)


def test_selection_scales_and_geometry_do_not_see_validation_or_test():
    names, d, _ = synthetic_dataset()
    labels, split = split_species(names, "order", 42)
    normalized, scales = normalize_distances(d, split["anchor"])
    panel, random, scores = select_panel(normalized, labels, split, 3, 42)
    changed = d.copy()
    for part in ("test", "validation", "unknown"):
        ids = split[part]
        changed[:, ids, :] = 1000
        changed[:, :, ids] = 1000
    n2, s2 = normalize_distances(changed, split["anchor"])
    p2, r2, sc2 = select_panel(n2, labels, split, 3, 42)
    np.testing.assert_array_equal(scales, s2)
    np.testing.assert_array_equal(panel, p2)
    np.testing.assert_array_equal(random, r2)
    np.testing.assert_allclose(scores, sc2)
    a1 = fit_pcoa(aggregate_profiles(normalized[:, panel][:, :, panel]))
    a2 = fit_pcoa(aggregate_profiles(n2[:, panel][:, :, panel]))
    np.testing.assert_allclose(a1["xy"], a2["xy"])
    assert not set(panel) & set(split["train"])


def test_pcoa_out_of_sample_reproduces_reference_and_euclidean_points():
    points = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    d = cdist(points, points)
    atlas = fit_pcoa(d)
    np.testing.assert_allclose(project_pcoa(d, atlas), atlas["xy"], atol=1e-10)
    q = np.array([[0.2, 0.6]])
    xy = project_pcoa(cdist(q, points), atlas)
    np.testing.assert_allclose(cdist(xy, atlas["xy"]), cdist(q, points), atol=1e-10)


def test_negative_eigenvalue_mass_is_reported():
    d = np.array([[0, 1, 1, 5], [1, 0, 1, 1], [1, 1, 0, 1], [5, 1, 1, 0]], dtype=float)
    atlas = fit_pcoa(d)
    assert atlas["negative_mass_fraction"] > 0


def test_layout_completion_is_explicit_preserves_observed_and_rejects_disconnection():
    d = np.array([[0, 1, np.nan], [1, 0, 2], [np.nan, 2, 0]])
    with pytest.raises(ValueError, match="lacks"):
        complete_anchor_distances(d)
    completed, count = complete_anchor_distances(d, "shortest-path")
    assert count == 1
    assert completed[0, 2] == 3
    np.testing.assert_array_equal(completed[np.isfinite(d)], d[np.isfinite(d)])
    d[0, 1] = d[1, 0] = np.nan
    with pytest.raises(ValueError, match="disconnected"):
        complete_anchor_distances(d, "shortest-path")


def test_raster_is_invariant_to_anchor_order_and_retains_missingness():
    d = np.array([[0.2, np.nan, 0.7]])
    c = np.array([[0.8, 0, 0.9]])
    xy = np.array([[0, 0], [0, 1], [1, 0]])
    ch = feature_rows(d, c)
    assert ch[0, 0, 1] == 0 and ch[0, 2, 1] == 0
    order = [2, 0, 1]
    kernel = raster_kernel(xy, 5)
    np.testing.assert_allclose(rasterize(ch, kernel), rasterize(ch[:, :, order], kernel[order]))
    assert rasterize(ch, kernel).shape == (1, 75)


def test_validation_threshold_rejects_when_evidence_insufficient():
    p = np.array([[0.9, 0.1], [0.8, 0.2], [0.2, 0.8]])
    assert choose_threshold(p, np.array(["A", "A", "B"]), np.array(["A", "B"])) > 1


def write_matrix(path, names, d):
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([""] + list(names))
        for name, row in zip(names, d):
            writer.writerow([name] + list(row))


def test_adapter_rejects_filtered_sentinel_and_duplicate_anchor(tmp_path):
    path = tmp_path / "q.matrix"
    anchor = "O_F_G_s"
    write_matrix(path, [anchor, "query1"], [[0, 999], [999, 0]])
    with pytest.raises(ValueError, match="flagged"):
        read_query_matrix(path, [anchor], "query")
    write_matrix(path, [anchor, anchor + "_accession", "query1"], np.zeros((3, 3)))
    with pytest.raises(ValueError, match="Duplicate anchor"):
        read_query_matrix(path, [anchor], "query")


def test_reference_matrix_preparation(tmp_path):
    names, d, c = synthetic_dataset()
    write_matrix(tmp_path / "d.csv", names, d[0])
    write_matrix(tmp_path / "c.csv", names, c[0])
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("locus,path,coverage_path\ngene1,d.csv,c.csv\n")
    output = tmp_path / "prepared.npz"
    prepare_matrices(SimpleNamespace(manifest=manifest, distance_kind="patristic", output=output))
    data = np.load(output)
    assert data["distances"].shape == (1, len(names), len(names))
    assert json.loads(str(data["metadata"]))["distance_kind"] == "patristic"


def test_benchmark_model_reload_and_query_inference(tmp_path):
    names, d, c = synthetic_dataset()
    dataset = tmp_path / "data.npz"
    np.savez(dataset, names=names, distances=d, coverage=c, loci=np.array(["a", "b", "c"]),
             metadata=json.dumps({"distance_kind": "patristic"}))
    output = tmp_path / "bench"
    with threadpool_limits(limits=2):
        benchmark(SimpleNamespace(dataset=dataset, output=output, rank="order", seed=42,
                                  panel_size=3, grid_size=4))
    metrics = json.loads((output / "metrics.json").read_text())
    assert metrics["results"]["distance_learned_panel"]["test"]["accuracy"] == 1
    bundle = joblib.load(output / "model.joblib")
    panel = list(bundle["panel"])
    indices = [list(names).index(n) for n in panel]
    query_species = metrics["split"]["test"][0]
    idx = list(names).index(query_species)
    square_ids = indices + [idx]
    write_matrix(tmp_path / "query.csv", panel + ["query1"], d[0][np.ix_(square_ids, square_ids)])
    write_matrix(tmp_path / "coverage.csv", panel + ["query1"], c[0][np.ix_(square_ids, square_ids)])
    manifest = tmp_path / "query_manifest.csv"
    manifest.write_text("locus,path,coverage_path\na,query.csv,coverage.csv\n")
    for method in ("distance_learned_panel", "raster_pattern"):
        args = SimpleNamespace(model=output / "model.joblib", manifest=manifest, query_prefix="query",
                               distance_kind="patristic", method=method, output=tmp_path / method)
        predict_matrices(args)
        rows = list(csv.DictReader((args.output / "query_predictions.csv").open()))
        assert rows[0]["prediction"] == query_species.split("_")[0]
        args.distance_kind = "pairwise-deletion-p"
        with pytest.raises(ValueError, match="distance kind"):
            predict_matrices(args)
