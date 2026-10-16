import pytest

from experiments.reference_atlas.real_truth import resolve_taxonomy, taxonomy_index, verify_complete_batch
from experiments.reference_atlas.real_evaluation import METHODS, score
from experiments.reference_atlas.parse_real_truth import expected_rows


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


def test_duplicate_sample_rows_cannot_inflate_scores():
    sample = {method: ["A"] for method in METHODS if method != "expanded_review_union"}
    sample.update(sample="real1", review=[])
    with pytest.raises(ValueError, match="sample IDs"):
        score([sample, sample], {"real1": ["A"]}, {"A"})


def test_s7_uses_expected_columns_and_library_ids_not_legacy_predictions():
    cells = {"B3": "Library ID", "D3": "Species", "I3": "Expected Orders", "J3": "Expected Families"}
    for i in range(1, 12):
        row = i + 3
        cells.update({f"B{row}": f"EM-{i}", f"C{row}": 1, f"D{row}": "G species",
                      f"E{row}": .1, f"H{row}": 100, f"I{row}": "TrueOrder", f"J{row}": "TrueFamily",
                      f"K{row}": "WrongOrder", f"L{row}": "WrongFamily"})
    rows = expected_rows(cells)
    assert len(rows) == 11
    assert rows[-1]["sample"] == "EMseq011"
    assert {r["family"] for r in rows} == {"TrueFamily"}
    cells["B14"] = "EM-12"
    with pytest.raises(ValueError, match="Unexpected library"):
        expected_rows(cells)
