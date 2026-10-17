"""Matched no-verifier diagnostic using the same recorded nearest candidates.

The quorum is the unchanged three-gene ML quorum. A second diagnostic adds the
existing C-style requirement for three uniquely discriminating genes (margin .1).
Plan/code were added during
post-evaluation development, before opening complete-pool scores. No truth is
read when constructing these predictions. Replaying accepted fragments must
exactly reproduce the sealed standalone-ML set.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .evaluation import complete, counts, dump, read, read_csv, sha, write_csv


def export(batch,output):
    batch,output=Path(batch),Path(output)
    selection,calls=complete(batch)
    rows=[];predictions=[];sources={}
    for sample in selection:
        path=batch/'samples'/sample['project']
        raw,accepted=defaultdict(set),defaultdict(set)
        source=path/'exon_predictions.csv'
        if source.exists():
            sources[str(source)]=sha(source)
            for r in read_csv(source):
                raw[r['genus']].add(r['gene'])
                if r['accepted']=='True':
                    accepted[r['genus']].add(r['gene'])
        ml={g for g,genes in accepted.items() if len(genes)>=3}
        if ml!=set(calls[sample['sample']]['standalone_ml']):
            raise ValueError('Accepted-fragment replay differs from frozen ML calls')
        selected=sorted(g for g,genes in raw.items() if len(genes)>=3)
        evidence=path/'genus_evidence.csv'
        info={r['genus']:r for r in read_csv(evidence)} if evidence.exists() else {}
        if evidence.exists():
            sources[str(evidence)]=sha(evidence)
        diagnostic=sorted(g for g in ml if int(info[g]['discriminating_genes'])>=3)
        predictions.append({'sample':sample['sample'],'nearest_gene3':selected,'ml_diagnostic_gene3':diagnostic})
        rows.extend({'sample':sample['sample'],'genus':g,'nearest_winner_gene_votes':len(raw[g])} for g in selected)
    output.mkdir(parents=True,exist_ok=False)
    write_csv(output/'predictions.csv',rows,['sample','genus','nearest_winner_gene_votes'])
    dump(output/'predictions.json',predictions)
    dump(output/'SEAL.json',{'sealed_utc':datetime.now(timezone.utc).isoformat(),
         'constructor_code_sha256':sha(__file__),
         'method':'nearest_gene3','same_pre_verifier_nearest_winners_and_gene_quorum':True,
         'additional_fixed_rule':'ml_diagnostic_gene3: ML >=3 genes AND >=3 distinct genes with unique nearest genus and sibling-relative margin >=.1',
         'accepted_fragment_replay_matches_standalone_ML':True,'truth_read_by_constructor':False,
         'scope':'Post-evaluation diagnostic; whole-method recommendation remains in-sample development',
         'source_exon_prediction_sha256':sources,
         'artifact_sha256':{p.name:sha(p) for p in output.iterdir() if p.is_file()}})


def evaluate(output,truth_dir):
    output,truth_dir=Path(output),Path(truth_dir)
    for name,digest in read(output/'SEAL.json')['artifact_sha256'].items():
        if sha(output/name)!=digest:
            raise ValueError('Nearest control predictions changed after sealing')
    expected={r['sample']:set(r['expected'].split(';')) for r in read_csv(truth_dir/'by_sample.csv') if r['method']=='merged'}
    rows=[]
    for r in read(output/'predictions.json'):
        for method in ('nearest_gene3','ml_diagnostic_gene3'):
            truth=expected[r['sample']];called=set(r[method])
            rows.append({'sample':r['sample'],'method':method,'components':len(r['sample'].split('x')),
                         **counts(len(truth&called),len(called-truth),len(truth-called)),
                         'expected':';'.join(sorted(truth)),'predicted':';'.join(sorted(called)),
                         'false_positives':';'.join(sorted(called-truth)),'false_negatives':';'.join(sorted(truth-called))})
    totals=[]
    for method in ('nearest_gene3','ml_diagnostic_gene3'):
        sub=[r for r in rows if r['method']==method]
        totals.append({'variant':'complete_joint','view':'table_literal','method':method,'samples':len(sub),
                       **counts(*(sum(r[k] for r in sub) for k in ('TP','FP','FN'))),
                       'exact_sets':sum(r['FP']==r['FN']==0 for r in sub)})
    write_csv(output/'by_sample.csv',rows);write_csv(output/'summary.csv',totals)
    return totals


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('export','evaluate'))
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--batch',type=Path)
    parser.add_argument('--truth-dir',type=Path)
    args=parser.parse_args()
    if args.command=='export':
        export(args.batch,args.output)
    else:
        print(evaluate(args.output,args.truth_dir))
