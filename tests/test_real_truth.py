import pytest

from experiments.reference_atlas.real_truth import resolve_taxonomy, taxonomy_index, verify_complete_batch
from experiments.reference_atlas.real_evaluation import METHODS, score


def test_taxonomy_resolves_only_without_cross_family_ambiguity(tmp_path):
    tree = tmp_path / "taxonomy.tree"
    tree.write_text("(O_F_G_a_1,O_F_G_b_2,P_H_J_c_3,Q_K_J_d_4);")
    species, genera = taxonomy_index(tree)
    assert resolve_taxonomy("G a", species, genera)["mapping_rule"].startswith("exact")
    assert resolve_taxonomy("G absent", species, genera)["family"] == "F"
    with pytest.raises(ValueError):
        resolve_taxonomy("J absent", species, genera)
    with pytest.raises(ValueError):
        resolve_taxonomy("Missing species", species, genera)


def test_truth_access_guard_rejects_incomplete_predictions(tmp_path):
    with pytest.raises(ValueError, match="before all predictions"):
        verify_complete_batch(tmp_path)


def test_family_truth_is_deduplicated_and_deferred_missing_family_stays_fn():
    sample = {method: ["A"] for method in METHODS if method != "expanded_review_union"}
    sample.update(sample="real1", review=["B"])
    rows, summary, universe = score([sample], {"real1": ["A", "A", "B", "MissingReference"]}, {"A", "B", "Negative"})
    merged = next(r for r in rows if r["method"] == "merged")
    assert (merged["TP"], merged["FP"], merged["FN"], merged["TN"]) == (1, 0, 2, 1)
    assert "MissingReference" in universe
    expanded = next(r for r in rows if r["method"] == "expanded_review_union")
    assert expanded["TP"] == 2 and expanded["FN"] == 1
    total = next(r for r in summary if r["method"] == "merged")
    assert total["macro_recall_present_families"] == 1 / 3
