"""Freeze and run the ten-sample adaptive genus experiment on TTU."""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import tarfile
import time
import traceback
from importlib.metadata import version

from .phylogeny import Cell,TreeIndex
from .references import ReferenceBank,dump
from .scoring import route,final_scores
from .trees import infer_exons


BASE = Path('/lustre/scratch/nhu/202609')
CACHE = BASE/'ttu_nocap100_z09_20261015/reference_cache'
BACKBONE = BASE/'871ref/backbone_corrected'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def csv_file(path,rows):
    if not rows:
        Path(path).write_text('')
        return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def names(path,values):
    Path(path).write_text(''.join(v+'\n' for v in sorted(set(values))))


def verify(hashes):
    for path,digest in hashes.items():
        if sha(path)!=digest:
            raise ValueError('Frozen file changed: '+path)


def prepare(root,backbone=None):
    root=Path(root).resolve()
    if (root/'PROTOCOL.json').exists():
        raise FileExistsError('This root is already frozen')
    backbone=Path(backbone) if backbone else BACKBONE
    if not backbone.is_dir():
        raise FileNotFoundError(backbone)
    policy=read(root/'code/experiments/adaptive_genus/policy.json')
    cohort=read(root/'code/experiments/adaptive_genus/cohort.json')
    if not (root/'model/SUCCESS').exists():
        raise ValueError('Genus-specific reference model must finish before freezing samples')
    rows=[]
    for row in cohort['selected']:
        source=Path(row['source_run'])/'00_query'
        genes=sorted(g for g in (source/'gene_list.txt').read_text().split() if (backbone/(g+'.fasta')).exists())
        files=sorted((source/'exon').glob('*.fasta'))+[source/'gene_list.txt']
        if len(genes)<10 or len(files)<10:
            raise ValueError('Insufficient reusable input: '+row['sample'])
        rows.append({**row,'source_exon':str(source/'exon'),'genes':genes,
                     'run':str(root/'samples'/row['project']),
                     'input_sha256':{str(p):sha(p) for p in files}})
    dump(root/'selection.json',{'selected':rows,'truth_used':False})
    code_files=[p for directory in (root/'code',root/'sources') for p in directory.rglob('*')
                if p.is_file() and '__pycache__' not in p.parts]
    code_files+=sorted(root.glob('*.sbatch'))
    verify(read(root/'model/TRAINING.json')['source_alignment_sha256'])
    metadata=[CACHE/n for n in ('release.json','sequence_manifest.tsv','by_recovery_index.json','treeoflife.current.tree')]
    all_genes=sorted({g for r in rows for g in r['genes']})
    refs=[backbone/(g+'.fasta') for g in all_genes]
    protocol={'schema':'adaptive-phylogenetic-genus-v1','rank':'genus','samples':10,
              'parent_source':'Union of existing original Z, frozen family merged and review; never expected taxa',
              'truth_used':False,'policy':policy,'backbone':str(backbone),'cache':str(CACHE),
              'code_sha256':{str(p):sha(p) for p in code_files},
              'model_sha256':sha(root/'model/model.joblib'),'training_sha256':sha(root/'model/TRAINING.json'),
              'selection_sha256':sha(root/'selection.json'),
              'reference_metadata_and_backbone_sha256':{str(p):sha(p) for p in metadata+refs},
              'software':{'python':platform.python_version(),**{p:version(p) for p in ('numpy','scipy','biopython','scikit-learn','joblib','httpx')}},
              'scope':'New genus experiment on previously analyzed mixtures; reference-only calibration, not a new independent biological cohort',
              'score_changes':'Genus-balanced contrast Z; same C confirmation/rescue and core/auxiliary forms; genus-specific model. Not the old frozen family rule.',
              'alignment_changes':'Cached Kew projection; reference-only trimming; reduced reference MAFFT addfragments/keeplength; fresh full exon trees each round'}
    dump(root/'PROTOCOL.json',protocol)
    print(json.dumps({'frozen_samples':10,'model_sha256':protocol['model_sha256'],'backbone':str(backbone)}),flush=True)


def verify_root(root):
    protocol=read(root/'PROTOCOL.json')
    verify(protocol['code_sha256'])
    verify(protocol['reference_metadata_and_backbone_sha256'])
    if (sha(root/'selection.json')!=protocol['selection_sha256']
        or sha(root/'model/model.joblib')!=protocol['model_sha256']
        or sha(root/'model/TRAINING.json')!=protocol['training_sha256']):
        raise ValueError('Selection/model changed')
    return protocol


def regroup(cells,index):
    groups=defaultdict(set)
    for cell in cells:
        groups[cell.family].update(cell.tips)
    return [Cell(f,tuple(sorted(tips))) for f,tips in sorted(groups.items())]


def execute(root,number):
    root=Path(root).resolve()
    protocol=verify_root(root);policy=protocol['policy']
    sample=read(root/'selection.json')['selected'][number]
    verify(sample['input_sha256'])
    run=Path(sample['run']);run.mkdir(parents=True,exist_ok=True)
    if (run/'SUCCESS').exists():
        return
    started=time.monotonic()
    def status(message):
        (run/'CURRENT_STAGE').write_text(message+'\n')
        print(sample['sample'],message,flush=True)
    try:
        from sprout_refbuilder.kew import KewClient
        status('load_tree_and_index')
        index=TreeIndex(Path(protocol['cache'])/'treeoflife.current.tree')
        with KewClient(Path(protocol['cache'])) as client:
            bank=ReferenceBank(index,client,protocol['backbone'],sample['genes'],sample['families'],
                               run/'reference_bank',workers=8,
                               minimum_fraction=policy['minimum_accepted_target_locus_fraction'])
            initial=[Cell(f,tuple(sorted(set(index.by_family[f])&set(bank.by_tip)))) for f in sample['families']]
            empty_families=[c.family for c in initial if not c.tips]
            active=[c for c in initial if c.tips]
            original={g for c in active for g in index.genera(c)}
            names(run/'initial_genera.txt',original)
            history=[];frozen_products={};stalls=0
            cells=active
            for round_number in range(1,policy['maximum_rounds']+1):
                genera={g for c in active for g in index.genera(c)}
                if len(genera)<=policy['final_joint_genus_budget']:
                    break
                directory=run/f'round_{round_number:02d}';directory.mkdir(exist_ok=True)
                status(f'round_{round_number:02d}_sample_{len(genera)}_genera')
                phase_start=time.monotonic()
                budget=policy['maximum_cells'] if stalls else policy['initial_cells']
                cells=index.partition(regroup(active,index),budget=max(budget,len(regroup(active,index))))
                selected_keys,anchor_cells,details=bank.anchors(cells,policy['anchor_qc_attempts'])
                if not selected_keys:
                    break
                names(directory/'candidate_genera.txt',genera)
                dump(directory/'plan.json',{'distance_kind':index.distance_kind,
                     'cells':[{'id':c.key,'family':c.family,'tips':list(c.tips),'genera':sorted(index.genera(c))} for c in cells],
                     'representatives':details,'anchor_cells':anchor_cells})
                ref,panel=bank.panel(selected_keys,directory/'reference')
                status(f'round_{round_number:02d}_trees_{len(selected_keys)}_references')
                trees,tree_info=infer_exons(ref,Path(sample['source_exon']),panel['genes'],directory/'trees',workers=8,
                                            minimum_bases=policy['minimum_query_bases'])
                chosen,evidence,anchor_votes=route(trees,anchor_cells,cells,index,
                    min_genes=policy['route_minimum_genes'],ambiguous_genes=policy['ambiguous_second_minimum_genes'],
                    gap=policy['ambiguous_distance_gap'])
                supported=[r['tip'] for r in details if not r['quality_unresolved']
                           and len(anchor_votes.get(bank.header(r['recovery']),set()))>=policy['route_minimum_genes']]
                kept=index.neighborhood_guard(cells,chosen,supported,policy['nearby_guard_genera'])
                remaining={g for c in kept for g in index.genera(c)}
                if not remaining<=genera:
                    raise ValueError('Search reintroduced pruned genera')
                csv_file(directory/'routing_evidence.csv',evidence)
                names(directory/'retained_genera.txt',remaining)
                record={'round':round_number,'initial_genera':len(genera),'retained_genera':len(remaining),
                        'reference_species':len(selected_keys),'sampled_genera':len({bank.entries[k]['taxonomy'].genus for k in selected_keys}),
                        'cells':len(cells),'selected_cells':len(chosen),'valid_genes':tree_info['valid_genes'],
                        'trees':len(trees),'quality_unresolved_cells':sum(r['quality_unresolved'] for r in evidence),
                        'seconds':time.monotonic()-phase_start,'tree_seconds':tree_info['seconds'],
                        'new_download_bytes_cumulative':bank.downloaded_bytes,
                        'qc_seconds_cumulative':bank.qc_seconds,'reference_alignment_seconds_cumulative':bank.align_seconds,
                        'reduction_fraction':1-len(remaining)/max(1,len(genera))}
                dump(directory/'ROUND.json',record);history.append(record)
                frozen_products.update({str(p):sha(p) for p in [*trees,*ref.glob('*.fasta')]})
                active=regroup(kept,index)
                if len(remaining)/max(1,len(genera))>=policy['stalled_retained_fraction']:
                    stalls+=1
                else:
                    stalls=0
                print(json.dumps(record),flush=True)
                if stalls>=policy['stalled_rounds']:
                    break
            remaining={g for c in active for g in index.genera(c)}
            unresolved=set()
            if len(remaining)>policy['final_joint_genus_budget']:
                # Only fully local neighborhoods are eligible for a partial
                # joint run. Coarse anchors are never promoted into genus calls.
                final_set={g for c in cells if len(index.genera(c))<=2
                           for g in index.genera(c)}&remaining
                unresolved=remaining-final_set
            else:
                final_set=remaining
            names(run/'unresolved_genera.txt',unresolved)
            names(run/'final_search_genera.txt',final_set)
            status(f'final_joint_panel_{len(final_set)}_genera')
            calls={name:[] for name in ('balanced_z','c_style','core','auxiliary','merged','standalone_ml','review')}
            final_record={};omitted=sorted(final_set)
            if final_set:
                keys,omitted,counts=bank.genus_panel(final_set,policy['final_species_per_genus'])
                if len(keys)>=4:
                    ref,panel=bank.panel(keys,run/'final_reference')
                    trees,tree_info=infer_exons(ref,Path(sample['source_exon']),panel['genes'],run/'final_trees',workers=8,
                                                minimum_bases=policy['minimum_query_bases'])
                    status('final_genus_C_ML_fusion')
                    calls,evidence,exons,statistics,features=final_scores(trees,root/'model/model.joblib',sample['project'])
                    csv_file(run/'genus_evidence.csv',evidence)
                    csv_file(run/'exon_predictions.csv',exons)
                    dump(run/'STATISTICS.json',statistics)
                    with gzip.open(run/'GENE_FEATURES.json.gz','wt') as f:
                        json.dump(features,f)
                    final_record={'reference_genera':sum(v>0 for v in counts.values()),'reference_species':len(keys),'omitted_qc_genera':omitted,
                                  'valid_genes':tree_info['valid_genes'],'valid_trees':len(trees),'tree_seconds':tree_info['seconds']}
                    frozen_products.update({str(p):sha(p) for p in [*trees,*ref.glob('*.fasta')]})
            for method,values in calls.items():
                if not set(values)<=final_set:
                    raise ValueError('Final call outside tested genus pool')
                names(run/(method+'_genera.txt'),values)
            dump(run/'predictions.json',{'sample':sample['sample'],'project':sample['project'],**calls})
            record={'sample':sample['sample'],'initial_genera':len(original),'remaining_search_genera':len(remaining),
                    'final_search_genera':len(final_set),'unresolved_genera':len(unresolved),
                    'search_status':'unresolved_partial' if unresolved or omitted or not final_record else 'refined',
                    'empty_parent_families':empty_families,'low_qc_final_genera':omitted,
                    'rounds':history,'final':final_record,'seconds':time.monotonic()-started,
                    'new_download_bytes':bank.downloaded_bytes,'fetch_seconds':bank.fetch_seconds,
                    'qc_seconds':bank.qc_seconds,'reference_alignment_seconds':bank.align_seconds,
                    'checked_recoveries':len(bank.quality),'model_threshold':read(root/'model/TRAINING.json')['threshold'],
                    'truth_used':False}
            dump(run/'SEARCH.json',record)
            verify(frozen_products);verify(bank.recovery_hashes)
            verify(sample['input_sha256']);verify_root(root)
            dump(run/'AUDIT.json',{'truth_used':False,'input_model_code_unchanged':True,
                                  'products_sha256':frozen_products,'recovery_sha256':bank.recovery_hashes})
            status('complete');(run/'SUCCESS').write_text('pipeline complete; see SEARCH.json for unresolved regions\n')
    except Exception:
        (run/'FAILED').write_text(traceback.format_exc())
        raise


def collect(root):
    root=Path(root).resolve();protocol=verify_root(root)
    verify(read(root/'model/TRAINING.json')['source_alignment_sha256'])
    rows=read(root/'selection.json')['selected']
    predictions,searches=[],[]
    for row in rows:
        run=Path(row['run'])
        if not (run/'SUCCESS').exists() or (run/'FAILED').exists():
            raise ValueError('All ten pipelines must finish, including unresolved cases')
        verify(row['input_sha256'])
        audit=read(run/'AUDIT.json');verify(audit['products_sha256']);verify(audit['recovery_sha256'])
        predictions.append(read(run/'predictions.json'));searches.append(read(run/'SEARCH.json'))
    dump(root/'predictions.json',predictions);dump(root/'SEARCHES.json',searches)
    dump(root/'AUDIT.json',{'all_10_complete':True,'code_model_inputs_products_verified':True,'truth_used':False,
                          'unresolved_samples':[r['sample'] for r in searches if r['search_status']!='refined']})
    (root/'PREDICTIONS_COMPLETE').write_text('10/10 pipelines finished before truth evaluation\n')
    with tarfile.open(root/'compact_results.tar.gz','w:gz') as archive:
        for name in ('PROTOCOL.json','selection.json','predictions.json','SEARCHES.json','AUDIT.json','PREDICTIONS_COMPLETE','model/TRAINING.json'):
            archive.add(root/name,arcname=name)
        for row in rows:
            run=Path(row['run'])
            for p in run.iterdir():
                if p.is_file():
                    archive.add(p,arcname=str(p.relative_to(root)))
            for directory in sorted(run.glob('round_*')):
                for name in ('candidate_genera.txt','retained_genera.txt','plan.json','routing_evidence.csv','ROUND.json','reference/panel.json','trees/trees.json'):
                    p=directory/name
                    if p.exists():
                        archive.add(p,arcname=str(p.relative_to(root)))
            for name in ('final_reference/panel.json','final_trees/trees.json'):
                p=run/name
                if p.exists():
                    archive.add(p,arcname=str(p.relative_to(root)))
    print('All 10 completed and audited',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','run','collect'))
    parser.add_argument('--root',required=True,type=Path)
    parser.add_argument('--backbone',type=Path)
    parser.add_argument('--index',type=int)
    args=parser.parse_args()
    if args.command=='prepare':
        prepare(args.root,args.backbone)
    elif args.command=='collect':
        collect(args.root)
    else:
        execute(args.root,args.index)
