"""Rerun order, family fusion and adaptive genus searches without old calls."""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import time
import traceback

from experiments.adaptive_genus.collect_parallel import sha_stream, verify_parallel
from experiments.adaptive_genus.evaluation import dump


def read(path):
    return json.loads(Path(path).read_text())


def names(path, values):
    Path(path).write_text(''.join(v+'\n' for v in sorted(set(values))))


def hashes(paths):
    from concurrent.futures import ThreadPoolExecutor
    paths=list(paths)
    with ThreadPoolExecutor(max_workers=8) as pool:
        return {str(p.resolve()):v for p,v in zip(paths,pool.map(sha_stream,paths))}


def verify(mapping):
    verify_parallel(mapping,workers=8)


def family_parents(calls):
    return sorted(set(calls['original_z'])|set(calls['merged'])|set(calls['review']))


def validate_cohort(rows):
    if len(rows)!=100 or len({r['sample'] for r in rows})!=100 or len({r['project'] for r in rows})!=100:
        raise ValueError('Expected exactly 100 unique mixture IDs and project outputs')
    for r in rows:
        if r['project']!='mix_'+r['sample'] or any(not p.isdigit() or not 1<=int(p)<=30 for p in r['sample'].split('x')):
            raise ValueError('Unexpected sample identity')


def prepare(root):
    if (root/'PROTOCOL.json').exists():
        raise FileExistsError('Use a fresh frozen run root')
    request=read(root/'REQUEST.json');cohort=read(root/'cohort.json')
    validate_cohort(cohort['selected'])
    rows=[]
    for i,r in enumerate(cohort['selected']):
        query=Path(request['assembly_source'])/'samples'/r['project']/'00_query'
        genes=sorted((query/'gene_list.txt').read_text().split())
        if len(genes)!=50 or len(set(genes))!=50:
            raise ValueError('Expected the frozen fifty target genes: '+r['sample'])
        files=sorted((query/'exon').glob('*.fasta'))+[query/'gene_list.txt']
        if len(files)<11:
            raise ValueError('Missing query assembly: '+r['sample'])
        rows.append({'index':i,'sample':r['sample'],'project':r['project'],'run':str(root/'samples'/r['project']),
                     'source_exon':str(query/'exon'),'genes':genes,'input_sha256':hashes(files),
                     'previously_seen':True,'genus_development_ten':r['genus_development_ten']})
    dump(root/'selection.json',{'selected':rows,'previous_predictions_imported':False})
    external=[Path(request['cache'])/n for n in ('release.json','sequence_manifest.tsv','by_recovery_index.json','treeoflife.current.tree')]
    external += [Path(request[k])/(g+'.fasta') for k in ('backbone','order_reference') for g in rows[0]['genes']]
    code=[p for folder in ('code','sources','models') for p in (root/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    code+=list(root.glob('*.sbatch'))+[root/n for n in ('PLAN.md','REQUEST.json','cohort.json','PAYLOAD.json')]
    for rank,expected in request['model_sha256'].items():
        if sha_stream(root/'models'/rank/'model.joblib')!=expected:
            raise ValueError('Changed '+rank+' checkpoint')
    dump(root/'PROTOCOL.json',{**request,'selection_sha256':sha_stream(root/'selection.json'),
         'code_sha256':hashes(code),'reference_metadata_and_backbone_sha256':hashes(external),
         'runtime':{'python':sys.version,**{p:importlib.metadata.version(p) for p in ('numpy','scipy','scikit-learn','biopython','joblib','pandas','httpx')}},
         'frozen_utc':datetime.now(timezone.utc).isoformat(),'truth_uploaded':False,'truth_used_by_pipeline':False})
    print('Frozen 100 order-to-genus tasks; only query assemblies reused',flush=True)


def verify_root(root):
    protocol=read(root/'PROTOCOL.json')
    if sha_stream(root/'selection.json')!=protocol['selection_sha256']:
        raise ValueError('Changed sample selection')
    verify(protocol['code_sha256']);verify(protocol['reference_metadata_and_backbone_sha256'])
    return protocol


class Timer:
    def __init__(self,output):
        self.output=output;self.started=time.monotonic();self.current=None;self.tick=self.started;self.rows=[]

    def stage(self,name):
        now=time.monotonic()
        if self.current is not None:
            self.rows.append({'phase':self.current,'seconds':now-self.tick})
        self.current=name;self.tick=now
        dump(self.output/'TIMING.json',{'completed_phases':self.rows,'current':name,'elapsed_seconds':now-self.started})
        (self.output/'CURRENT_STAGE').write_text(name+'\n')
        print(name,flush=True)

    def finish(self):
        self.stage('complete')
        dump(self.output/'TIMING.json',{'phases':self.rows,'inference_seconds':sum(r['seconds'] for r in self.rows)})


def command(args,cwd,log):
    with log.open('a') as handle:
        subprocess.run([sys.executable,*map(str,args)],cwd=cwd,stdout=handle,stderr=subprocess.STDOUT,check=True)


def order_stage(root,sample,protocol,timer):
    from experiments.reference_atlas.real_validation import validate_trees,uncapped_orders,write_csv
    from export_pea import gene_evidence,read_predictions
    run=Path(sample['run']);stage=run/'01_order';stage.mkdir(exist_ok=True)
    (stage/'exon').symlink_to(sample['source_exon'],target_is_directory=True)
    names(stage/'gene_list.txt',sample['genes'])
    timer.stage('order_align_tree')
    command([root/'sources/02_exon_trees.py','-t',8,'-p',sample['project'],'-g',stage/'gene_list.txt',
             '-e',sample['source_exon'],'-r',protocol['order_reference'],'--output_phylo',stage/'phylo'],stage,stage/'COMMANDS.log')
    validate_trees(stage/'phylo',stage/'valid_phylo',sample['project'])
    timer.stage('order_distances')
    command([root/'sources/03_distance_matrices.py','-t',8,'-p',sample['project'],'-g',stage/'gene_list.txt',
             '--query-only','--input_phylo',stage/'valid_phylo','--output_tree',stage/'order_matrix'],stage,stage/'COMMANDS.log')
    timer.stage('order_prediction')
    command([root/'sources/04_prediction.py','-i',stage/(sample['project']+'.cumulative_dist.csv'),
             '-o',stage/'native_predictions.csv','-tl','o','-z',.9,'-to',stage/'native_z_candidates.txt'],stage,stage/'COMMANDS.log')
    genes,coverage=gene_evidence(stage/'order_matrix',sample['project'],'order')
    document={'genes':genes,'predictions':read_predictions(stage/'native_predictions.csv',None,None)}
    orders,evidence=uncapped_orders(document)
    if not orders:
        raise ValueError('No order parent admitted')
    names(stage/'order_candidates.txt',orders);write_csv(stage/'order_evidence.csv',evidence)
    dump(stage/'ORDER.json',{'orders':sorted(orders),'coverage':coverage,'evidence':document})
    return {**sample,'orders':sorted(orders),'order_evidence':evidence,'source_order':str(stage),'source_genes':sample['genes']}


def genus_stage(root,sample,protocol,families,timer):
    from experiments.adaptive_genus import runner as genus
    base=Path(sample['run'])/'03_genus';base.mkdir()
    (base/'model').symlink_to(root/'models/genus',target_is_directory=True)
    row={**sample,'run':str(base/'result'),'families':family_parents(families),
         'family_search_routes':{k:families[k] for k in ('original_z','merged','review')}}
    dump(base/'selection.json',{'selected':[row],'truth_used':False})
    policy=read(root/'code/experiments/adaptive_genus/policy.json')
    dump(base/'PROTOCOL.json',{**protocol,'policy':policy,'selection_sha256':sha_stream(base/'selection.json'),
         'model_sha256':sha_stream(root/'models/genus/model.joblib'),'training_sha256':sha_stream(root/'models/genus/TRAINING.json')})
    genus.verify=verify
    timer.stage('genus_adaptive_search_and_final')
    genus.execute(base,0)
    return read(base/'result/predictions.json')


def run(root,index):
    protocol=verify_root(root);sample=read(root/'selection.json')['selected'][index]
    output=Path(sample['run']);output.mkdir(parents=True,exist_ok=True)
    if (output/'SUCCESS').exists():
        return
    if (output/'ATTEMPT_STARTED').exists():
        raise FileExistsError('Preserve incomplete attempt before a controlled retry: '+str(output))
    verify(sample['input_sha256']);(output/'ATTEMPT_STARTED').write_text(datetime.now(timezone.utc).isoformat())
    timer=Timer(output)
    try:
        runtime=order_stage(root,sample,protocol,timer)
        from sprout_refbuilder.clade import BalancedKewPanelBuilder
        from sprout_refbuilder.kew import KewClient
        timer.stage('family_prefetch')
        with KewClient(Path(protocol['cache'])) as client:
            prefetch=BalancedKewPanelBuilder(Path(protocol['backbone']),output/'prefetch_only',
                representatives_per_taxon=4,pool_multiplier=2,minimum_recovered_loci=50,
                download_workers=4).prefetch('family',runtime['orders'],client)
        dump(output/'FAMILY_PREFETCH.json',prefetch)
        from experiments.reference_atlas import real_validation as family
        family.CACHE=Path(protocol['cache']);family.BACKBONE=Path(protocol['backbone'])
        family.execute_stages(root,output,runtime,status_callback=lambda n:timer.stage('family_'+n),tree_workers=8)
        families=read(output/'predictions.json')
        dump(output/'FAMILY_AUDIT.json',read(output/'AUDIT.json'))
        genera=genus_stage(root,sample,protocol,families,timer)
        dump(output/'predictions.json',{'sample':sample['sample'],'project':sample['project'],
             'orders':runtime['orders'],'family':families,'genus':genera,
             'primary_family':'merged','primary_genus':'auxiliary'})
        timer.finish()
        products=[p for d in ('01_order','02_family') for p in (output/d).rglob('*')
                  if p.is_file() and p.suffix in ('.tre','.fasta','.csv')]
        dump(output/'ORDER_FAMILY_PRODUCTS.json',hashes(products))
        verify(sample['input_sha256']);verify_root(root)
        dump(output/'AUDIT.json',{'fresh_order_trees_and_predictions':True,'fresh_family_references_trees_predictions':True,
             'fresh_adaptive_genus_search':True,'previous_taxonomic_predictions_imported':False,
             'truth_used_by_pipeline':False,'model_sources_and_input_hashes_verified':True,
             'order_family_products_manifest_sha256':sha_stream(output/'ORDER_FAMILY_PRODUCTS.json'),
             'family_audit_sha256':sha_stream(output/'FAMILY_AUDIT.json'),
             'genus_audit_sha256':sha_stream(output/'03_genus/result/AUDIT.json')})
        (output/'SUCCESS').write_text('Fresh order/family/genus inference complete; inspect unresolved statuses\n')
    except Exception:
        (output/'FAILED').write_text(traceback.format_exc());raise


def collect(root):
    protocol=verify_root(root);rows=read(root/'selection.json')['selected'];validate_cohort(rows)
    status=[{'sample':r['sample'],'success':(Path(r['run'])/'SUCCESS').exists(),
             'failed':(Path(r['run'])/'FAILED').exists(),
             'stage':(Path(r['run'])/'CURRENT_STAGE').read_text().strip() if (Path(r['run'])/'CURRENT_STAGE').exists() else 'not_started'} for r in rows]
    dump(root/'STATUS.json',status)
    if not all(r['success'] and not r['failed'] for r in status):
        raise ValueError('All 100 must finish before sealing/evaluation; see STATUS.json')
    predictions=[];timings=[]
    for row in rows:
        output=Path(row['run']);verify(row['input_sha256']);verify(read(output/'ORDER_FAMILY_PRODUCTS.json'))
        ga=read(output/'03_genus/result/AUDIT.json');verify(ga['products_sha256']);verify(ga['recovery_sha256'])
        predictions.append(read(output/'predictions.json'));timings.append({'sample':row['sample'],**read(output/'TIMING.json')})
    dump(root/'predictions.json',predictions);dump(root/'TIMINGS.json',timings)
    dump(root/'AUDIT.json',{'all_100_complete':True,'source_models_inputs_products_verified':True,'truth_used_by_pipeline':False})
    (root/'PREDICTIONS_COMPLETE').write_text('100/100 fresh order-to-genus predictions audited\n')
    with tarfile.open(root/'compact_results.tar.gz','w:gz') as archive:
        for name in ('PROTOCOL.json','selection.json','predictions.json','TIMINGS.json','STATUS.json','AUDIT.json','PREDICTIONS_COMPLETE','PLAN.md','PAYLOAD.json'):
            archive.add(root/name,arcname=name)
        for row in rows:
            output=Path(row['run'])
            for p in output.rglob('*'):
                if not p.is_file() or p.is_symlink() or 'exon' in p.parts or 'family_reference' in p.parts or 'reference_bank' in p.parts:
                    continue
                if p.suffix in ('.json','.txt','.gz') or p.name in ('genus_evidence.csv','exon_predictions.csv','merged_evidence.csv','order_evidence.csv','family_evidence.csv','routing_evidence.csv'):
                    archive.add(p,arcname=str(p.relative_to(root)))
    print('100/100 results collected and audited',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','run','collect'))
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--index',type=int)
    args=parser.parse_args();root=args.root.resolve()
    if args.command=='run':run(root,args.index)
    elif args.command=='prepare':prepare(root)
    else:collect(root)
