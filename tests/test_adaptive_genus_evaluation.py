"""Evaluation gates must not turn incomplete searches into successful calls."""
import json

import pytest

from experiments.adaptive_genus.evaluation import complete, counts, export, loss_stage, verify_export


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def test_incomplete_batch_is_rejected_before_truth_access(tmp_path):
    put(tmp_path/'AUDIT.json', {'all_10_complete':False})
    with pytest.raises(ValueError, match='All ten'):
        complete(tmp_path)


def test_absent_calls_still_count_as_false_negatives():
    assert counts(0, 0, 7)=={'TP':0, 'FP':0, 'FN':7, 'precision':None, 'recall':0, 'F1':0}


def test_first_loss_distinguishes_parent_routing_and_unresolved(tmp_path):
    row={'families':['F']}
    (tmp_path/'initial_genera.txt').write_text('A\nB\n')
    (tmp_path/'unresolved_genera.txt').write_text('B\n')
    round_path=tmp_path/'round_01';round_path.mkdir()
    (round_path/'retained_genera.txt').write_text('B\n')
    assert loss_stage('A','missing',row,tmp_path,set())=='absent_parent_family'
    assert loss_stage('C','F',row,tmp_path,set())=='absent_recoverable_initial_reference'
    assert loss_stage('A','F',row,tmp_path,set())=='round_01_routing_exclusion'
    assert loss_stage('B','F',row,tmp_path,set())=='unresolved_search_region'


def test_sealed_export_has_ten_samples_and_detects_tampering(tmp_path):
    from experiments.adaptive_genus.evaluation import METHODS
    batch=tmp_path/'batch';batch.mkdir()
    put(batch/'AUDIT.json',{'all_10_complete':True,'code_model_inputs_products_verified':True})
    (batch/'PREDICTIONS_COMPLETE').write_text('done')
    put(batch/'PROTOCOL.json',{'synthetic':True})
    put(batch/'SEARCHES.json',[])
    rows=[];calls=[]
    for i in range(10):
        row={'sample':str(i),'project':'mix_'+str(i)};rows.append(row)
        call={**row,**{m:[] for m in METHODS}};calls.append(call)
        path=batch/'samples'/row['project'];path.mkdir(parents=True)
        (path/'SUCCESS').write_text('done')
        put(path/'predictions.json',call)
        put(path/'SEARCH.json',{'search_status':'unresolved_partial'})
        for m in METHODS:
            (path/(m+'_genera.txt')).write_text('')
        for name in ('initial_genera','final_search_genera','unresolved_genera'):
            (path/(name+'.txt')).write_text('')
    put(batch/'selection.json',{'selected':rows})
    put(batch/'predictions.json',calls)
    output=tmp_path/'export'
    export(batch,output)
    verify_export(output)
    (output/'predictions.json').write_text('[]')
    with pytest.raises(ValueError,match='artifact changed'):
        verify_export(output)
