import json
from pathlib import Path

import pytest

from experiments.hierarchy100.evaluation import complete,genus_loss
from experiments.hierarchy100.runner import validate_cohort


def cohort():
    source=Path(__file__).resolve().parents[1]/'experiments/reference_atlas/results/transfer100_20261016/selection.json'
    return [{**r,'genus_development_ten':i<10} for i,r in enumerate(json.loads(source.read_text())['selected'])]


def put(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value))


def test_incomplete_hundred_cannot_be_scored_even_if_top_level_audit_claims_success(tmp_path):
    rows=cohort();put(tmp_path/'selection.json',{'selected':rows})
    put(tmp_path/'AUDIT.json',{'all_100_complete':True,'source_models_inputs_products_verified':True})
    (tmp_path/'PREDICTIONS_COMPLETE').write_text('100/100')
    put(tmp_path/'predictions.json',[{'sample':r['sample']} for r in rows[:-1]])
    with pytest.raises(ValueError,match='Full cohort prediction mismatch'):complete(tmp_path)


def test_hundred_sample_directories_must_not_collide():
    rows=cohort();rows[1]['project']=rows[0]['project']
    with pytest.raises(ValueError,match='unique'):validate_cohort(rows)


def test_sealing_rechecks_genus_text_instead_of_trusting_combined_json(tmp_path):
    rows=cohort();put(tmp_path/'selection.json',{'selected':rows})
    put(tmp_path/'AUDIT.json',{'all_100_complete':True,'source_models_inputs_products_verified':True})
    (tmp_path/'PREDICTIONS_COMPLETE').write_text('100/100')
    calls=[{'sample':r['sample'],'orders':[],'family':{},'genus':{'auxiliary':['Citrus']}} for r in rows]
    put(tmp_path/'predictions.json',calls)
    first=tmp_path/'samples'/rows[0]['project'];put(first/'predictions.json',calls[0])
    (first/'01_order').mkdir();(first/'01_order/order_candidates.txt').write_text('')
    directory=first/'03_genus/result';directory.mkdir(parents=True)
    (directory/'auxiliary_genera.txt').write_text('Alepidea\n')
    with pytest.raises(ValueError,match='Genus raw set mismatch'):complete(tmp_path)


def test_unresolved_true_genus_is_a_miss_before_any_classifier_claim(tmp_path):
    (tmp_path/'initial_genera.txt').write_text('Citrus\n')
    (tmp_path/'unresolved_genera.txt').write_text('Citrus\n')
    assert genus_loss(tmp_path,'Citrus',set())=='unresolved_search_region'
