"""Complete a frozen strict search without dropping supported unresolved genera.

This post-evaluation development run consumes frozen candidate pools, not truth.
The genus model, reference QC, tree inference and final scoring are unchanged.
"""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
from pathlib import Path
import shutil
import time
import traceback

from . import runner
from .collect_parallel import sha_stream, verify_parallel
from .references import ReferenceBank, dump
from .phylogeny import TreeIndex
from .scoring import final_scores
from .trees import infer_exons

runner.verify = verify_parallel
METHODS = ('balanced_z','c_style','core','auxiliary','merged','standalone_ml','review')


def names(path):
    return set(Path(path).read_text().splitlines())


def hashes(paths):
    paths=list(paths)
    with ThreadPoolExecutor(max_workers=16) as pool:
        return {str(p):value for p,value in zip(paths,pool.map(sha_stream,paths))}


def final_pool(remaining,budget=512):
    remaining=set(remaining)
    return (remaining,set()) if len(remaining)<=budget else (set(),remaining)


def prepare(root):
    root=Path(root).resolve()
    if (root/'PROTOCOL.json').exists():
        raise FileExistsError('Use a fresh completion root')
    reuse=runner.read(root/'coarse_reuse_manifest.json')
    verify_parallel(reuse['global_sha256'])
    source=Path(reuse['source_root'])
    parent=runner.read(source/'PROTOCOL.json')
    verify_parallel(parent['code_sha256'])
    verify_parallel(parent['reference_metadata_and_backbone_sha256'])
    if sha_stream(root/'model/model.joblib')!=parent['model_sha256'] or sha_stream(root/'model/TRAINING.json')!=parent['training_sha256']:
        raise ValueError('Completion must use the unchanged genus model')
    rows=[]
    for row in runner.read(source/'selection.json')['selected']:
        metadata=reuse['samples'][row['sample']]
        verify_parallel(metadata['metadata_sha256'])
        rows.append({**row,'run':str(root/'samples'/row['project']),
                     'source_search_run':row['run'],**metadata})
    dump(root/'selection.json',{'selected':rows,'truth_used_by_pipeline':False,
                               'human_development_used_previous_evaluation':True})
    code=[p for folder in (root/'code',root/'sources') for p in folder.rglob('*')
          if p.is_file() and '__pycache__' not in p.parts]
    code+=list(root.glob('*.sbatch'))+[root/'coarse_reuse_manifest.json',root/'POSTHOC_COMPLETION_PLAN.md']
    protocol={**parent,'schema':'adaptive-genus-complete-joint-v3',
              'policy':runner.read(root/'code/experiments/adaptive_genus/policy.json'),
              'code_sha256':hashes(code),'selection_sha256':sha_stream(root/'selection.json'),
              'source_strict_root':str(source),'source_protocol_sha256':sha_stream(source/'PROTOCOL.json'),
              'truth_used_by_pipeline':False,'human_development_used_previous_evaluation':True,
              'scope':'Post-evaluation development correction on the same ten mixtures; unchanged strict search and genus model',
              'completion_rule':'Jointly verify every retained genus at <=512; no partial small-cell substitution',
              'cost_scope':'SEARCH.seconds includes the reused source pipeline plus incremental completion, including the original partial final run'}
    dump(root/'PROTOCOL.json',protocol)
    print('Frozen complete-pool correction for all ten mixtures',flush=True)


def reuse_reference_cache(bank,source):
    """Use fixed-coordinate positive records from audited source panels."""
    from sprout_refbuilder.fasta import read_fasta
    reports=sorted(source.glob('round_*/reference/panel.json'))
    if (source/'final_reference/panel.json').exists():
        reports.append(source/'final_reference/panel.json')
    positive={}
    for path in reports:
        report=runner.read(path)
        for row in report['selection']:
            key=row['recovery']
            if key not in bank.entries or bank.header(key)!=row['header']:
                raise ValueError('Changed recovery identity in reference cache')
            if row['accepted_loci']<bank.minimum:
                raise ValueError('Source reference failed the unchanged target-locus QC')
            positive[key]=row
        for gene in report['genes']:
            for sequence in read_fasta(path.parent/'ref'/(gene+'.fasta')):
                previous=bank.aligned[gene].get(sequence.id)
                if previous is not None and previous.sequence!=sequence.sequence:
                    raise ValueError('Conflicting fixed coordinates in reused reference cache')
                bank.aligned[gene][sequence.id]=sequence
    for key,row in positive.items():
        if any(row['header'] in bank.aligned[g] for g in bank.genes):
            bank.quality[key]={'accepted_loci':row['accepted_loci'],'minimum_loci':bank.minimum,
                               'passed':True,'accepted':{}}
    for gene,rows in bank.aligned.items():
        if len({len(row.sequence) for row in rows.values()})>1:
            raise ValueError('Reused projections do not share fixed coordinates: '+gene)
    positive_count=len(bank.quality)
    qc_path=source/'reference_bank/qc.jsonl'
    if qc_path.exists():
        observed=set();accepted=defaultdict(set)
        with qc_path.open() as handle:
            for line in handle:
                row=json.loads(line);identity=row['sequence_id'];observed.add(identity)
                if row['accepted']:
                    accepted[identity].add(row['gene_id'])
        by_id=defaultdict(list)
        for key,entry in bank.entries.items():
            by_id[entry['record'].sequence_id].append(key)
        for identity in observed:
            if len(by_id[identity])==1 and len(accepted[identity])<bank.minimum:
                key=by_id[identity][0]
                if key not in bank.quality:
                    bank.quality[key]={'accepted_loci':len(accepted[identity]),'minimum_loci':bank.minimum,
                                       'passed':False,'accepted':{}}
    return positive_count


def link(source,destination):
    destination.symlink_to(source,target_is_directory=True)


def run(root,number):
    root=Path(root).resolve();protocol=runner.verify_root(root)
    sample=runner.read(root/'selection.json')['selected'][number]
    verify_parallel(sample['input_sha256']);verify_parallel(sample['metadata_sha256'])
    verify_parallel(sample['reference_sha256'])
    source=Path(sample['source_search_run']);output=Path(sample['run'])
    output.mkdir(parents=True,exist_ok=True)
    if (output/'SUCCESS').exists():
        return
    started=time.monotonic()
    def status(value):
        (output/'CURRENT_STAGE').write_text(value+'\n')
        print(sample['sample'],value,flush=True)
    try:
        previous=runner.read(source/'SEARCH.json');audit=runner.read(source/'AUDIT.json')
        remaining=names(source/('round_%02d/retained_genera.txt'%previous['rounds'][-1]['round'])) if previous['rounds'] else names(source/'initial_genera.txt')
        if len(remaining)!=previous['remaining_search_genera']:
            raise ValueError('Frozen retained pool count mismatch')
        for directory in sorted(source.glob('round_*')):
            link(directory,output/directory.name)
        shutil.copy2(source/'initial_genera.txt',output/'initial_genera.txt')
        products=dict(audit['products_sha256']);recoveries=dict(audit['recovery_sha256'])
        if remaining==names(source/'final_search_genera.txt') and not names(source/'unresolved_genera.txt'):
            status('reuse_existing_complete_pool_final')
            for path in source.iterdir():
                if path.is_file() and path.name not in ('AUDIT.json','SEARCH.json','SUCCESS','CURRENT_STAGE'):
                    shutil.copy2(path,output/path.name)
            for name in ('final_reference','final_trees'):
                link(source/name,output/name)
            record={**previous,'completion_final_reused':True,'incremental_completion_seconds':time.monotonic()-started,
                    'source_pipeline_seconds':previous['seconds'],'development_after_truth_evaluation':True}
            record['seconds']=previous['seconds']+record['incremental_completion_seconds']
        else:
            selected,unresolved=final_pool(remaining,protocol['policy']['completion_joint_genus_budget'])
            runner.names(output/'final_search_genera.txt',selected)
            runner.names(output/'unresolved_genera.txt',unresolved)
            calls={method:[] for method in METHODS};final={};omitted=sorted(selected)
            cache_count=cache_total=0;fetch=qc=align=downloaded=0;checked=0
            if selected:
                status('complete_joint_panel_%d_genera'%len(selected))
                from sprout_refbuilder.kew import KewClient
                index=TreeIndex(Path(protocol['cache'])/'treeoflife.current.tree')
                with KewClient(Path(protocol['cache'])) as client:
                    bank=ReferenceBank(index,client,protocol['backbone'],sample['genes'],sample['families'],
                                       output/'reference_bank',workers=8,
                                       minimum_fraction=protocol['policy']['minimum_accepted_target_locus_fraction'])
                    cache_count=reuse_reference_cache(bank,source)
                    cache_total=len(bank.quality)
                    status('cached_%d_positive_recoveries'%cache_count)
                    keys,omitted,counts=bank.genus_panel(selected,protocol['policy']['final_species_per_genus'])
                    if len(keys)>=4:
                        reference,panel=bank.panel(keys,output/'final_reference')
                        status('infer_joint_trees_%d_references'%len(keys))
                        trees,info=infer_exons(reference,Path(sample['source_exon']),panel['genes'],output/'final_trees',
                                               workers=8,minimum_bases=protocol['policy']['minimum_query_bases'])
                        status('unchanged_genus_C_ML_fusion')
                        calls,evidence,exons,statistics,features=final_scores(trees,root/'model/model.joblib',sample['project'])
                        runner.csv_file(output/'genus_evidence.csv',evidence)
                        runner.csv_file(output/'exon_predictions.csv',exons)
                        dump(output/'STATISTICS.json',statistics)
                        with gzip.open(output/'GENE_FEATURES.json.gz','wt') as handle:
                            json.dump(features,handle)
                        final={'reference_genera':sum(v>0 for v in counts.values()),'reference_species':len(keys),
                               'omitted_qc_genera':omitted,'valid_genes':info['valid_genes'],
                               'valid_trees':len(trees),'tree_seconds':info['seconds']}
                        products.update(hashes([*trees,*reference.glob('*.fasta')]))
                    for path,digest in bank.recovery_hashes.items():
                        if path in recoveries and recoveries[path]!=digest:
                            raise ValueError('Previously frozen recovery changed: '+path)
                        recoveries[path]=digest
                    fetch,qc,align=bank.fetch_seconds,bank.qc_seconds,bank.align_seconds
                    downloaded,checked=bank.downloaded_bytes,len(bank.quality)
            for method,values in calls.items():
                if not set(values)<=selected:
                    raise ValueError('Final prediction outside the retained search pool')
                runner.names(output/(method+'_genera.txt'),values)
            dump(output/'predictions.json',{'sample':sample['sample'],'project':sample['project'],**calls})
            elapsed=time.monotonic()-started
            record={**previous,'final_search_genera':len(selected),'unresolved_genera':len(unresolved),
                    'search_status':'unresolved_partial' if unresolved or omitted or not final else 'refined',
                    'low_qc_final_genera':omitted,'final':final,
                    'seconds':previous['seconds']+elapsed,'source_pipeline_seconds':previous['seconds'],
                    'incremental_completion_seconds':elapsed,'completion_final_reused':False,
                    'source_partial_final_tree_seconds':previous['final'].get('tree_seconds',0),
                    'cached_positive_recoveries':cache_count,'new_download_bytes':previous['new_download_bytes']+downloaded,
                    'cached_negative_recoveries':cache_total-cache_count,
                    'fetch_seconds':previous['fetch_seconds']+fetch,'qc_seconds':previous['qc_seconds']+qc,
                    'reference_alignment_seconds':previous['reference_alignment_seconds']+align,
                    'checked_recoveries':previous['checked_recoveries']+checked-cache_total,
                    'development_after_truth_evaluation':True}
        dump(output/'SEARCH.json',record)
        verify_parallel(products);verify_parallel(recoveries)
        verify_parallel(sample['metadata_sha256']);verify_parallel(sample['input_sha256'])
        runner.verify_root(root)
        dump(output/'AUDIT.json',{'truth_used_by_pipeline':False,'human_development_used_previous_evaluation':True,
                                 'input_model_code_unchanged':True,'products_sha256':products,'recovery_sha256':recoveries,
                                 'source_search':str(source),'source_metadata_sha256':sample['metadata_sha256'],
                                 'cached_reference_projections_verified':True})
        status('complete');(output/'SUCCESS').write_text('Complete retained-pool run; see SEARCH.json for QC gaps\n')
    except Exception:
        (output/'FAILED').write_text(traceback.format_exc())
        raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','run'))
    parser.add_argument('--root',required=True,type=Path)
    parser.add_argument('--index',type=int)
    args=parser.parse_args()
    if args.command=='prepare':
        prepare(args.root)
    else:
        run(args.root,args.index)
