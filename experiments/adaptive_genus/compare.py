"""Compare the two predeclared, fully sealed routing policies on identical data."""
import argparse
from pathlib import Path
from statistics import median

from .evaluation import dump, read, read_csv, verify_export, write_csv


def normalized_code(protocol):
    return {key.split('/code/',1)[-1] if '/code/' in key else key.split('/sources/',1)[-1]:value
            for key,value in protocol['code_sha256'].items()
            if not key.endswith('.sbatch') and not key.endswith('/policy.json')}


def compare(root):
    root=Path(root);paths={'v1':root/'v1','strict_route':root/'strict_route'}
    for p in paths.values():
        verify_export(p)
    a,b=(read(paths[k]/'PROTOCOL.json') for k in paths)
    if normalized_code(a)!=normalized_code(b):
        raise ValueError('Inference code differs beyond policy and submission wrappers')
    for key in ('model_sha256','training_sha256','reference_metadata_and_backbone_sha256'):
        if a[key]!=b[key]:
            raise ValueError('Frozen input/model mismatch: '+key)
    rows={key:read(p/'selection.json')['selected'] for key,p in paths.items()}
    same_keys=('sample','project','families','family_search_routes','source_exon','genes','input_sha256')
    if [{k:r[k] for k in same_keys} for r in rows['v1']] != [{k:r[k] for k in same_keys} for r in rows['strict_route']]:
        raise ValueError('Cohort, parent evidence or query inputs differ')
    differences={key:{'v1':a['policy'][key],'strict_route':b['policy'][key]} for key in a['policy'] if a['policy'][key]!=b['policy'][key]}
    if set(differences)!={'schema','route_minimum_genes','ambiguous_second_minimum_genes','nearby_guard_genera','evaluation'}:
        raise ValueError('Unexpected policy change')
    if read(paths['v1']/'truth.json') != read(paths['strict_route']/'truth.json'):
        raise ValueError('Truth views differ')
    results=[];summary=[]
    for variant,p in paths.items():
        timing={r['sample']:r for r in read_csv(p/'timings.csv')}
        for result in read_csv(p/'by_sample.csv'):
            if result['method']=='merged':
                t=timing[result['sample']]
                results.append({'variant':variant,**result,**{k:t[k] for k in ('initial_genera','remaining_search_genera','final_search_genera',
                                'final_reference_genera','unresolved_genera','search_status','rounds','seconds','new_download_bytes')}})
        for result in read_csv(p/'summary.csv'):
            summary.append({'variant':variant,**result})
    write_csv(root/'comparison_by_sample.csv',results)
    write_csv(root/'comparison_summary.csv',summary)
    dump(root/'COMPARISON_AUDIT.json',{'same_ten_query_inputs_and_parent_evidence':True,'same_model_and_training':True,
         'same_backbone_and_reference_metadata':True,'inference_code_identical_except_policy_and_job_wrappers':True,
         'same_literal_truth_and_used_cells':True,'policy_differences':differences,
         'v2_chosen_before_truth_from_v1_search_stagnation':True,'cache_order':'v1 launched first; both use shared raw Kew cache',
         'scope':'Retrospective algorithm development; no independent validation or exhaustive speed baseline'})
    pct=lambda value:'NA' if value=='' else f'{float(value):.1%}'
    lines=['# Adaptive genus search: ten-mixture TTU experiment','',
           'Completed both predeclared routing policies on the same ten in-silico mixtures: **20/20 pipelines audited**. '
           'Both batches were sealed before the original Table S1 genus truth was opened. '
           'The stricter policy was specified from truth-blinded search stagnation, with the same model and query data.', '',
           '| Routing policy; merged classifier | TP | FP | FN | Precision | Recall | F1 | Exact sets | Median retained genera | Median search minutes | ≤64 genera / full QC |',
           '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    stats={}
    for variant,p in paths.items():
        r=next(s for s in summary if s['variant']==variant and s['method']=='merged')
        t=read_csv(p/'timings.csv')
        stats[variant]={'accuracy':r,'median_pool':median(int(x['remaining_search_genera']) for x in t),
                        'median_minutes':median(float(x['seconds'])/60 for x in t),
                        'median_reduction':median(1-int(x['remaining_search_genera'])/int(x['initial_genera']) for x in t),
                        'bounded':sum(int(x['remaining_search_genera'])<=64 for x in t),
                        'refined':sum(x['search_status']=='refined' for x in t)}
        s=stats[variant]
        lines.append(f'| {variant} | {r["TP"]} | {r["FP"]} | {r["FN"]} | {pct(r["precision"])} | {pct(r["recall"])} | '
                     f'{float(r["F1"]):.3f} | {r["exact_sets"]}/10 | {s["median_pool"]:g} | {s["median_minutes"]:.1f} | {s["bounded"]}/10 / {s["refined"]}/10 |')
    dump(root/'COMPARISON_METRICS.json',stats)
    lines+=['','v1 requires two routing genes (three for ambiguous second routes) and four nearby guard genera per '
            'supported anchor. The strict policy requires five genes for both routes and one guard genus. '
            'Both use 32 initial tree cells, 64 after a stalled round, at most eight rounds and a 64-genus '
            'final joint-panel budget. These constrain reference work, not the number of positive predictions.', '',
            'Each round chooses central, target-locus-QC-passing Kew representatives of tree neighborhoods, '
            'aligns the same query exons against the reduced references, infers new exon trees and retains '
            'regions with repeated independent-gene support. Surviving regions are sampled again. '
            'The Kew planning tree supplies topology distances, not dated evolutionary branch lengths. '
            'Final calls combine the genus C-style rule and a reference-trained genus ML verifier. '
            'Coarse anchor genera are never promoted directly into final calls.', '',
            '## Results and reusable outputs','',
            '- [Strict routing: detailed results, trajectories and errors](strict_route/REPORT.md)',
            '- [Original routing: detailed results, trajectories and errors](v1/REPORT.md)',
            '- [Strict routing merged genus predictions](strict_route/merged_predictions.csv)',
            '- [Original routing merged genus predictions](v1/merged_predictions.csv)',
            '- [Paired sample comparison](comparison_by_sample.csv)',
            '- [Input/model/policy comparison audit](COMPARISON_AUDIT.json)','',
            'Review sets, unresolved genera, alternative final rules and all per-round candidates are retained '
            'in each result directory; none are silently added to the merged predictions. '
            'Each detailed report identifies the first stage losing a true genus, including upstream parent '
            'misses and reference QC failures.', '',
            '## Interpretation limits','',
            'This is retrospective development on previously analyzed mixtures drawn from the same 30 source species. '
            'The ten selected mixtures contain four 3-component, three 6-component and three 10-component samples. '
            'A preference based on these ten results needs validation on a separate cohort. '
            'The old progressive ten-mixture cohort is different, so it is not a matched accuracy/runtime control.', '',
            'Reference-pool reduction is measured; speedup versus exhaustive genus inference is not. '
            'The second policy can benefit from earlier downloads into the shared cache. '
            'Concurrency was a throughput choice rather than a controlled timing variable: '
            'v1 allowed four tasks, while strict routing increased from two to four to six. '
            'Per-sample search minutes include reference sampling, fetching/QC, reference projection, query alignment, '
            'trees and scoring; they exclude prior assembly, initial/final audits, Slurm setup/cleanup and queueing. '
            'Allocated job durations and CPU accounting are archived separately. '
            'Shared model training costs 5.9 minutes and is not counted again for the second policy.', '',
            'The original Table S1 is the literal scoring source. Prior family-level source discrepancies for '
            'IDs 1 and 25 are disclosed; no genus corrections were inferred from the new predictions. '
            'Species presence, abundance, negative-only mixtures and independent ML generalization are not evaluated.', '',
            'Branch: `experiment/adaptive-genus-search-20261017`. All implementation and evaluation code is separate '
            'from previous production defaults. The exact frozen inference payloads, shared trained model, '
            'source hashes, Slurm records, raw result archives and independent audits are included. '
            'See [design](../../DESIGN.md) and [evaluation plan](../../EVALUATION_PLAN.md).','']
    (root/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True,type=Path)
    compare(parser.parse_args().root)
