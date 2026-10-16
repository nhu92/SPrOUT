import json

import pytest

from experiments.reference_atlas.fuse import DEFAULT_POLICY, boolean, decide, keyed


POLICY = json.loads(DEFAULT_POLICY.read_text())


def row(**changes):
    return {"c_selected": True, "b_selected": True, "rescue_selected": True,
            "stability_lower": 0.7, "ml_gene_votes": 2, **changes}


def test_stable_rescue_can_restore_one_ml_gene_but_requires_all_channels():
    assert decide(row(ml_gene_votes=1), POLICY)["route"] == "auxiliary"
    for changes in ({"c_selected": False}, {"rescue_selected": False},
                    {"stability_lower": 0.5}, {"ml_gene_votes": 0}):
        candidate = row(ml_gene_votes=1)
        candidate.update(changes)
        assert not decide(candidate, POLICY)["merged"]


def test_ml_only_candidates_remain_review_and_route_overlap_is_deduplicated():
    result = decide(row(c_selected=False, b_selected=False, ml_gene_votes=8), POLICY)
    assert result["review"] and not result["merged"]
    both = decide(row(), POLICY)
    assert both["route"] == "core_and_auxiliary" and both["merged"] and not both["review"]
    absent = decide(row(c_selected=False, b_selected=False, ml_gene_votes=0), POLICY)
    assert absent["route"] == "not_selected"


def test_low_stability_does_not_override_two_gene_core_confirmation():
    result = decide(row(stability_lower=0.1, rescue_selected=False), POLICY)
    assert result["core"] and result["merged"] and not result["auxiliary"]


def test_bad_evidence_cannot_silently_merge():
    with pytest.raises(ValueError):
        boolean("true")
    with pytest.raises(ValueError):
        keyed([{"family": "A"}, {"family": "A"}], "family")
    with pytest.raises(ValueError):
        decide(row(ml_gene_votes=-1), POLICY)
