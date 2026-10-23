import json
from pathlib import Path
import subprocess
import sys

import pytest

from experiments.hierarchy100.evaluation import complete,genus_loss
from experiments.hierarchy100.runner import validate_cohort


def test_local_evaluator_import_does_not_require_inference_scientific_packages():
    repo=Path(__file__).resolve().parents[1]
    result=subprocess.run([sys.executable,'-S','-c','import experiments.hierarchy100.evaluation'],
                          cwd=repo,capture_output=True,text=True)
    assert result.returncode==0,result.stderr


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


def test_all_hundred_unresolved_samples_still_contribute_full_recall_denominator(tmp_path,monkeypatch):
    from experiments.hierarchy100 import evaluation as ev
    rows=cohort();batch=tmp_path/'batch';output=tmp_path/'export'
    put(batch/'selection.json',{'selected':rows})
    put(batch/'AUDIT.json',{'all_100_complete':True,'source_models_inputs_products_verified':True})
    (batch/'PREDICTIONS_COMPLETE').write_text('100/100')
    for name in ('PROTOCOL.json','PAYLOAD.json'):put(batch/name,{})
    put(batch/'STATUS.json',[])
    (batch/'PLAN.md').write_bytes(Path(ev.__file__).with_name('PLAN.md').read_bytes())
    sources={i:{'genus':'Genus'+str(i),'family':'Family','order':'Order'} for i in range(1,31)}
    calls=[];timings=[]
    for row in rows:
        directory=batch/'samples'/row['project'];directory.mkdir(parents=True)
        call={'sample':row['sample'],'orders':['Order'],'family':{'merged':[]},'genus':{'auxiliary':[]}}
        calls.append(call);put(directory/'predictions.json',call)
        timings.append({'sample':row['sample'],'inference_seconds':1,'phases':[{'phase':'synthetic_test','seconds':1}]})
        put(directory/'TIMING.json',{'phases':[]})
        (directory/'01_order').mkdir();(directory/'01_order/order_candidates.txt').write_text('Order\n')
        (directory/'merged_families.txt').write_text('')
        gd=directory/'03_genus/result';gd.mkdir(parents=True)
        for name in ('initial_genera.txt','unresolved_genera.txt'):
            (gd/name).write_text(''.join(r['genus']+'\n' for r in sources.values()))
        for name in ('final_search_genera.txt','auxiliary_genera.txt'):(gd/name).write_text('')
        put(gd/'SEARCH.json',{'rounds':[],'initial_genera':30,'remaining_search_genera':30,'unresolved_genera':30,'final':{},
             'search_status':'unresolved_partial','seconds':1,'fetch_seconds':0,'qc_seconds':0,'reference_alignment_seconds':0})
    put(batch/'TIMINGS.json',timings)
    put(batch/'predictions.json',calls)
    ev.seal(batch,output)
    monkeypatch.setattr(ev,'truth',lambda path:(sources,{}))
    ev.evaluate(batch,output,tmp_path/'never-opened.xlsx')
    summary=ev.read(output/'EVALUATION_AUDIT.json')
    assert summary['genus_first_losses']=={'unresolved_search_region':560}
    assert len(ev.read(output/'predictions.json'))==100
    from experiments.hierarchy100.report import report
    report(batch,output)
    assert (output/'validation.png').stat().st_size>1000
    assert (output/'REPORT.md').exists()
