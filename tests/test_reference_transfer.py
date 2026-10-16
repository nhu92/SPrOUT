from io import StringIO
import json

import joblib
import numpy as np
from Bio import Phylo
from sklearn.dummy import DummyClassifier
from experiments.reference_atlas.transfer import profiles, tree_distances, reference_block, predict_sample
from experiments.reference_atlas.evaluate_transfer import measures


def test_tree_matrix_matches_biopython_and_does_not_edit_source(tmp_path):
    path = tmp_path / "tree.tre"
    content = "((O_F_G_a:0.1,O_F_H_b:0.2):0.3,(O_X_J_c:0.4,O_X_K_d:0.5):0.6);"
    path.write_text(content)
    names, distances = tree_distances(path)
    tree = Phylo.read(StringIO(content), "newick")
    expected = [[tree.distance(a, b) for b in names] for a in names]
    np.testing.assert_allclose(distances, expected)
    assert path.read_text() == content


def test_full_fingerprint_removes_additive_query_branch_offset():
    ref = np.array([[0, 1, 5, 6], [1, 0, 6, 5], [5, 6, 0, 1], [6, 5, 1, 0]], dtype=float)
    taxa = np.array([["O", "F", "a"], ["O", "F", "b"], ["O", "X", "c"], ["O", "X", "d"]])
    q = np.array([0.5, 0.7, 5.2, 5.5])
    x, labels = profiles(q, ref, taxa)
    shifted, _ = profiles(q + 2, ref, taxa)
    np.testing.assert_allclose(x[:, 3], shifted[:, 3])
    np.testing.assert_allclose(x[:, 4], shifted[:, 4])
    assert labels == [("F", "O"), ("X", "O")]
    assert x[0, 0] < x[1, 0]


def test_anchors_are_species_deduplicated_and_parent_scope_is_respected(tmp_path):
    path = tmp_path / "tree.tre"
    path.write_text("((O_F_G_a_source1:0.1,O_F_G_a_source2:0.2):0.1,"
                    "O_F_H_b:0.2,P_X_J_c:0.3,P_X_K_d:0.4,sample_NODE_1:0.1);")
    d, taxa, q, names = reference_block(path, "sample")
    assert len(taxa) == 4
    assert names == ["sample_NODE_1"]
    x, labels = profiles(q[0], d, taxa, {"O"})
    assert labels == [("F", "O")]
    assert x.shape == (1, 7)


def test_confusion_metrics_count_false_negatives_and_undefined_empty_precision():
    result = measures(8, 2, 4, 86)
    assert result["precision"] == 0.8
    assert result["recall"] == 8 / 12
    assert result["F1"] == 16 / 22
    assert result["NPV"] == 86 / 90
    assert measures(0, 0, 4, 96)["precision"] is None
    assert measures(0, 0, 4, 96)["recall"] == 0


def test_reuse_prediction_keeps_inputs_and_counts_genes_once(tmp_path):
    source = tmp_path / "source"
    sample = source / "samples/sample"
    stage = sample / "02_family"
    phylo = stage / "valid_phylo"
    phylo.mkdir(parents=True)
    (sample / "01_order").mkdir()
    (sample / "01_order/order_candidates.txt").write_text("O\n")
    trees = [f"g{i}_exon_1.tre" for i in range(10)] + ["g0_exon_2.tre"]
    content = "((O_F_G_a:0.1,O_F_H_b:0.2,sample_NODE_1:0.1):0.2,(O_X_J_c:0.2,O_X_K_d:0.1):0.5);"
    for name in trees:
        (phylo / name).write_text(content)
    (stage / "family_tree_validation.json").write_text(json.dumps({"valid_trees": trees}))
    (stage / "family_panel_report.json").write_text("{}")
    for name in ("family_candidates.txt", "native_z_candidates.txt"):
        (stage / name).write_text("F\n")
    (source / "selection.json").write_text(json.dumps({"selected": [
        {"run": str(sample), "sample": "test", "project": "sample", "previously_seen": False}]}))
    model_dir = tmp_path / "model"
    model_dir.mkdir()
    model = DummyClassifier(strategy="constant", constant=1).fit(np.zeros((2, 7)), [0, 1])
    joblib.dump({"model": model, "threshold": 0.5}, model_dir / "model.joblib")
    output = tmp_path / "output"
    predict_sample(source, model_dir, output, 0)
    out = output / "samples/sample"
    assert (out / "profile_gene3.txt").read_text() == "F\n"
    audit = json.loads((out / "audit.json").read_text())
    assert audit["valid_genes"] == 10 and audit["valid_trees"] == 11
    assert all((phylo / name).read_text() == content for name in trees)
    import csv
    with (out / "family_evidence.csv").open() as f:
        evidence = {r["family"]: r for r in csv.DictReader(f)}
    assert evidence["F"]["profile_gene_votes"] == "10"
