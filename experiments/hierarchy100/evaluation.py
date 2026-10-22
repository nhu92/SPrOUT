"""Local-only sealing and literal Table S1 scoring of the hundred fresh runs.

Never upload the truth workbook or this evaluation's derived labels for inference.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import shutil

from experiments.adaptive_genus.evaluation import read,sha,dump,names,write_csv,counts,table_truth,TRUTH_SHA256
from experiments.hierarchy100.runner import validate_cohort


def complete(batch):
    batch=Path(batch);rows=read(batch/'selection.json')['selected'];validate_cohort(rows)
    audit=read(batch/'AUDIT.json')
    if not (batch/'PREDICTIONS_COMPLETE').exists() or not audit.get('all_100_complete') or not audit.get('source_models_inputs_products_verified'):
        raise ValueError('All 100 audited fresh pipelines must finish before scoring')
    calls=read(batch/'predictions.json');by_id={r['sample']:r for r in calls}
    if len(calls)!=100 or len(by_id)!=100 or set(by_id)!={r['sample'] for r in rows}:
        raise ValueError('Full cohort prediction mismatch')
    for row in rows:
        directory=batch/'samples'/row['project'];r=by_id[row['sample']]
        if read(directory/'predictions.json')!=r:
            raise ValueError('Per-sample and combined outputs differ')
        if set(r['orders'])!=names(directory/'01_order/order_candidates.txt'):
            raise ValueError('Order raw set mismatch')
        for method,values in r['family'].items():
            if method not in ('sample','project') and set(values)!=names(directory/(method+'_families.txt')):
                raise ValueError('Family raw set mismatch')
        for method,values in r['genus'].items():
            if method not in ('sample','project') and set(values)!=names(directory/'03_genus/result'/(method+'_genera.txt')):
                raise ValueError('Genus raw set mismatch')
    return rows,by_id


def seal(batch,output):
    batch,output=Path(batch),Path(output)
    rows,calls=complete(batch)
    if output.exists():raise FileExistsError('Use a fresh result export')
    output.mkdir(parents=True)
    (output/'.gitattributes').write_text('* -text whitespace=cr-at-eol\n')
    for name in ('PROTOCOL.json','selection.json','predictions.json','TIMINGS.json','STATUS.json','AUDIT.json','PREDICTIONS_COMPLETE','PLAN.md','PAYLOAD.json'):
        shutil.copy2(batch/name,output/name)
    predictions=[]
    for row in rows:
        r=calls[row['sample']]
        for rank,methods in [('order',{'admission':r['orders']}),('family',r['family']),('genus',r['genus'])]:
            for method,taxa in methods.items():
                if method in ('sample','project'):continue
                predictions.extend({'sample':row['sample'],'rank':rank,'method':method,'taxon':taxon,
                     'genus_development_ten':row['genus_development_ten']} for taxon in sorted(taxa))
    write_csv(output/'all_method_predictions.csv',predictions)
    primary={('order','admission'),('family','merged'),('genus','auxiliary')}
    write_csv(output/'primary_predictions.csv',[r for r in predictions if (r['rank'],r['method']) in primary])
    dump(output/'PREDICTION_SEAL.json',{'samples':100,'sealed_utc':datetime.now(timezone.utc).isoformat(),
         'truth_read_by_exporter':False,'retrospective_development':True,
         'plan_sha256':sha(Path(__file__).with_name('PLAN.md')),
         'source_batch_sha256':{p.relative_to(batch).as_posix():sha(p) for p in batch.rglob('*') if p.is_file()},
         'artifact_sha256':{p.name:sha(p) for p in output.iterdir() if p.is_file()}})


def verify_seal(output,batch=None):
    record=read(output/'PREDICTION_SEAL.json')
    if record['plan_sha256']!=sha(Path(__file__).with_name('PLAN.md')):raise ValueError('Changed evaluation plan')
    for name,digest in record['artifact_sha256'].items():
        if sha(output/name)!=digest:raise ValueError('Changed sealed artifact: '+name)
    if batch is not None:
        for name,digest in record['source_batch_sha256'].items():
            if sha(batch/name)!=digest:raise ValueError('Changed raw batch artifact: '+name)


def truth(workbook):
    from openpyxl import load_workbook
    from experiments.reference_atlas.real_truth import xml_sheet
    sources,cells=table_truth(workbook,'B',True)
    book=load_workbook(workbook,read_only=True,data_only=True);sheet=book['Table S1'];xml=xml_sheet(workbook,'Table S1')
    for identity,row in sources.items():
        address='G'+row['cells'][0][1:];value=str(sheet[address].value).strip()
        if value!=str(xml[address]).strip():raise ValueError('Order truth XML disagreement: '+address)
        row['order']=value;cells[address]=sheet[address].value
    book.close()
    if sha(workbook)!=TRUTH_SHA256:raise ValueError('Truth workbook changed')
    return sources,cells


def genus_loss(directory,genus,called):
    if genus not in names(directory/'initial_genera.txt'):return 'absent_initial_reference_or_parent'
    for path in sorted(directory.glob('round_*')):
        if genus not in names(path/'retained_genera.txt'):return path.name+'_routing_exclusion'
    if genus in names(directory/'unresolved_genera.txt'):return 'unresolved_search_region'
    if genus not in names(directory/'final_search_genera.txt'):return 'outside_final_joint_search'
    path=directory/'final_reference/panel.json';panel=read(path) if path.exists() else {'selection':[]}
    if genus not in {r['genus'] for r in panel['selection']}:return 'final_reference_QC_unavailable'
    if genus not in called:return 'final_classifier_rejection'
    return 'recovered'


def score(expected,called):
    return counts(len(expected&called),len(called-expected),len(expected-called))


def evaluate(batch,output,workbook):
    batch,output,workbook=map(Path,(batch,output,workbook));rows,calls=complete(batch);verify_seal(output,batch)
    for name in ('PROTOCOL.json','selection.json','predictions.json'):
        if sha(batch/name)!=sha(output/name):raise ValueError('Batch differs from sealed export')
    sources,cells=truth(workbook)
    dump(output/'truth.json',{'sheet':'Table S1','view':'literal','workbook_sha256':TRUTH_SHA256,
         'sources':list(sources.values()),'independently_verified_cells':cells})
    metrics=[];losses=[];timings=[];rounds=[]
    for row in rows:
        sample=row['sample'];directory=batch/'samples'/row['project'];call=calls[sample]
        source=[sources[int(i)] for i in sample.split('x')]
        for rank,methods in [('order',{'admission':call['orders']}),('family',call['family']),('genus',call['genus'])]:
            expected={r[rank] for r in source}
            for method,values in methods.items():
                if method in ('sample','project'):continue
                called=set(values)
                metrics.append({'sample':sample,'rank':rank,'method':method,'components':len(source),
                     'genus_development_ten':row['genus_development_ten'],**score(expected,called),
                     'expected':';'.join(sorted(expected)),'predicted':';'.join(sorted(called)),
                     'false_positives':';'.join(sorted(called-expected)),'false_negatives':';'.join(sorted(expected-called))})
        gd=directory/'03_genus/result';primary=set(call['genus']['auxiliary'])
        for genus in sorted({r['genus'] for r in source}):
            losses.append({'sample':sample,'genus':genus,'first_loss':genus_loss(gd,genus,primary),
                           'genus_development_ten':row['genus_development_ten']})
        for phase in read(directory/'TIMING.json')['phases']:
            timings.append({'sample':sample,**phase})
        search=read(gd/'SEARCH.json')
        for r in search['rounds']:rounds.append({'sample':sample,**r})
    summaries=[]
    groups=[('all100',lambda r:True),('genus_development10',lambda r:r['genus_development_ten']),
            ('other90',lambda r:not r['genus_development_ten'])]
    groups += [('components'+str(n),lambda r,n=n:r['components']==n) for n in (3,6,10)]
    for group,predicate in groups:
        for rank,method in dict.fromkeys((r['rank'],r['method']) for r in metrics):
            sub=[r for r in metrics if r['rank']==rank and r['method']==method and predicate(r)]
            summaries.append({'subset':group,'rank':rank,'method':method,'samples':len(sub),
                 **counts(*(sum(r[k] for r in sub) for k in ('TP','FP','FN'))),
                 'exact_sets':sum(not r['FP'] and not r['FN'] for r in sub)})
    for name,data in [('summary.csv',summaries),('by_sample.csv',metrics),('genus_first_losses.csv',losses),
                       ('timings.csv',timings),('genus_rounds.csv',rounds)]:write_csv(output/name,data)
    primary=next(r for r in summaries if r['subset']=='all100' and r['rank']=='genus' and r['method']=='auxiliary')
    if sum(r['first_loss']=='recovered' for r in losses)!=primary['TP']:raise ValueError('Loss audit disagrees with primary recall')
    dump(output/'EVALUATION_AUDIT.json',{'all_100_sealed_before_this_evaluation':True,'literal_truth_openpyxl_XML_agree':True,
         'raw_per_sample_sets_match_combined_predictions':True,'workbook_sha256':TRUTH_SHA256,
         'genus_first_losses':dict(Counter(r['first_loss'] for r in losses)),
         'previous_genus_development_ten':sum(r['genus_development_ten'] for r in rows),
         'independent_biological_validation':False})
    print([r for r in summaries if r['subset']=='all100' and (r['rank'],r['method']) in [('order','admission'),('family','merged'),('genus','auxiliary')]])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('seal','evaluate'))
    parser.add_argument('--batch',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workbook',type=Path);args=parser.parse_args()
    if args.command=='seal':seal(args.batch,args.output)
    else:evaluate(args.batch,args.output,args.workbook)
