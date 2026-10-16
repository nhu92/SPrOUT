import pytest

from experiments.reference_atlas.real_validation import check_order_rule, uncapped_orders, validate_trees


def test_reused_order_selection_must_preserve_the_uncapped_union():
    rows = [{"taxon": "A", "locus_votes": "0", "aggregate_z": ".8", "aggregate_rank": "1"},
            {"taxon": "B", "locus_votes": "2", "aggregate_z": "-.4", "aggregate_rank": "2"},
            {"taxon": "C", "locus_votes": "1", "aggregate_z": ".9", "aggregate_rank": "3"},
            {"taxon": "D", "locus_votes": "0", "aggregate_z": "", "aggregate_rank": ""}]
    check_order_rule(rows, {"A", "B", "C"})
    with pytest.raises(ValueError):
        check_order_rule(rows, {"A", "C"})


def test_full_tree_validation_rejects_duplicate_labels_and_preserves_sources(tmp_path):
    source = tmp_path / "trees"
    source.mkdir()
    for i in range(10):
        (source / f"g{i}_exon_1.tre").write_text("(O_F_G_a:0.1,emseq001_NODE_1:0.2);")
    invalid = source / "bad_exon_1.tre"
    content = "(O_F_G_a:0.1,O_F_G_a:0.2,emseq001_NODE_1:0.2);"
    invalid.write_text(content)
    result = validate_trees(source, tmp_path / "valid", "emseq001")
    assert result["distinct_genes"] == 10
    assert result["excluded_trees"] == [{"tree": invalid.name, "error": "duplicate terminal names"}]
    assert invalid.read_text() == content


def test_legacy_cap_and_exon_vote_counts_do_not_change_distinct_gene_admission():
    stage = {"genes": [{"gene": str(i), "matches": []} for i in range(10)],
             "predictions": [{"taxon": "A", "z_score": 4, "rank": 1},
                             {"taxon": "B", "z_score": -.3, "rank": 2},
                             {"taxon": "C", "z_score": -.4, "rank": 3}]}
    stage["genes"][0]["matches"] = [{"taxon": "B", "query_votes": 200}, {"taxon": "C", "query_votes": 1}]
    stage["genes"][1]["matches"] = [{"taxon": "C", "query_votes": 1}]
    called, rows = uncapped_orders(stage)
    assert called == {"A", "C"}
    assert next(r for r in rows if r["taxon"] == "B")["distinct_gene_votes"] == 1
