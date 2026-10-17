"""Compare the sealed baselines with the disclosed complete-pool correction."""
import argparse
from pathlib import Path
import shutil
from statistics import median
from zipfile import ZipFile, ZIP_DEFLATED

from .compare import compare, normalized_code
from .evaluation import dump, read, read_csv, sha, verify_export, write_csv
from .nearest_control import evaluate as evaluate_nearest


def region_sets(path):
    result={}
    for row in read_csv(path/'search_regions.csv'):
        result.setdefault((row['sample'],row['stage']),set()).add(row['genus'])
    return result


def plot(root,summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    variants=('v1','strict_route','complete_joint')
    labels=('Partial v1','Partial strict','Complete pool\n(post-evaluation)')
    fig,axes=plt.subplots(1,2,figsize=(11,4.7),constrained_layout=True)
    x=np.arange(3);width=.24
    for offset,field,color,label in [(-1,'precision','#4477aa','Precision'),(0,'recall','#228833','Recall'),(1,'F1','#aa3377','F1')]:
        values=[float(next(r for r in summary if r['variant']==v and r['method']=='merged')[field]) for v in variants]
        axes[0].bar(x+offset*width,values,width,color=color,label=label)
    axes[0].set_xticks(x,labels);axes[0].set(ylim=(0,1.13),ylabel='Micro metric',title='A  Same merged classifier')
    axes[0].legend(fontsize=8,ncol=3)
    stages=('initial','after_routing','final_search','final_reference','merged')
    all_values={}
    for variant in variants:
        rows=read_csv(root/variant/'truth_survival.csv');searches=read(root/variant/'SEARCHES.json')
        lookup={(r['sample'],r['stage']):int(r['true_retained']) for r in rows}
        values=[0]*5
        for search in searches:
            sample=search['sample'];last=f'round_{search["rounds"][-1]["round"]:02d}' if search['rounds'] else 'initial'
            for i,stage in enumerate(('initial',last,'final_search','final_reference','merged')):
                values[i]+=lookup[sample,stage]
        all_values[variant]=values
    x=np.arange(5)
    for offset,variant,color,label in [(-1,'v1','#999999','Partial v1'),(0,'strict_route','#4477aa','Partial strict'),(1,'complete_joint','#228833','Complete pool')]:
        axes[1].bar(x+offset*width,all_values[variant],width,color=color,label=label)
    axes[1].set_xticks(x,['Initial','Routed','Joint pool','QC refs','Merged'])
    axes[1].set(ylim=(0,72),ylabel='True sample–genus events',title='B  Where the original stop rule lost recall')
    axes[1].legend(fontsize=8,ncol=3)
    for ax in axes:
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
        ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    for suffix in ('png','svg'):
        fig.savefig(root/('comparison.'+suffix),dpi=220)
    plt.close(fig)
    path=root/'comparison.svg'
    path.write_text('\n'.join(line.rstrip() for line in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    return all_values


def final_report(root):
    root=Path(root)
    compare(root)
    shutil.copy2(root/'REPORT.md',root/'BASELINE_COMPARISON.md')
    paths={v:root/v for v in ('v1','strict_route','complete_joint')}
    for path in paths.values():
        verify_export(path)
    strict=read(paths['strict_route']/'PROTOCOL.json');completed=read(paths['complete_joint']/'PROTOCOL.json')
    if not completed.get('human_development_used_previous_evaluation') or completed['source_protocol_sha256']!=sha(paths['strict_route']/'PROTOCOL.json'):
        raise ValueError('Completion must identify the unchanged strict source and disclose previous evaluation')
    for key in ('model_sha256','training_sha256','reference_metadata_and_backbone_sha256'):
        if strict[key]!=completed[key]:
            raise ValueError('Final model or input changed: '+key)
    old_code,new_code=normalized_code(strict),normalized_code(completed)
    for key,value in old_code.items():
        if key!='experiments/adaptive_genus/runner.py' and new_code.get(key)!=value:
            raise ValueError('Changed QC, alignment, classifier or dependency: '+key)
    same_keys=('sample','project','families','family_search_routes','source_exon','genes','input_sha256')
    rows={v:read(path/'selection.json')['selected'] for v,path in paths.items()}
    if [{k:r[k] for k in same_keys} for r in rows['strict_route']] != [{k:r[k] for k in same_keys} for r in rows['complete_joint']]:
        raise ValueError('Completion changed the mixture cohort or parent evidence')
    old_regions,new_regions=region_sets(paths['strict_route']),region_sets(paths['complete_joint'])
    for key,values in old_regions.items():
        if key[1]=='initial' or key[1].startswith('round_'):
            if new_regions.get(key,set())!=values:
                raise ValueError('Completion changed a frozen routing pool')
    for search in read(paths['strict_route']/'SEARCHES.json'):
        sample=search['sample'];last=f'round_{search["rounds"][-1]["round"]:02d}_after' if search['rounds'] else 'initial'
        if new_regions.get((sample,'final_search'),set())!=old_regions.get((sample,last),set()):
            raise ValueError('Completion did not test the entire frozen retained pool')
    if read(paths['strict_route']/'truth.json')!=read(paths['complete_joint']/'truth.json'):
        raise ValueError('Literal truth changed')
    all_summary=[];by_sample=[];stats={}
    for variant,path in paths.items():
        all_summary.extend({'variant':variant,**r} for r in read_csv(path/'summary.csv'))
        by_sample.extend({'variant':variant,**r} for r in read_csv(path/'by_sample.csv'))
        timing=read_csv(path/'timings.csv')
        stats[variant]={'median_retained_genera':median(int(r['remaining_search_genera']) for r in timing),
                        'median_reduction':median(1-int(r['remaining_search_genera'])/int(r['initial_genera']) for r in timing),
                        'median_combined_search_minutes':median(float(r['seconds'])/60 for r in timing),
                        'target64_reached':sum(int(r['remaining_search_genera'])<=64 for r in timing),
                        'no_unresolved_search_regions':sum(int(r['unresolved_genera'])==0 for r in timing)}
    control=root/'nearest_control'
    if (control/'SEAL.json').exists():
        all_summary.extend(evaluate_nearest(control,paths['complete_joint']))
        by_sample.extend({'variant':'complete_joint',**r} for r in read_csv(control/'by_sample.csv'))
    write_csv(root/'all_variant_summary.csv',all_summary)
    write_csv(root/'all_variant_by_sample.csv',by_sample)
    options=[r for r in all_summary if r['variant']=='complete_joint' and r['method']!='review']
    best=max(options,key=lambda r:(float(r['F1']),float(r['precision'] or 0),float(r['recall']),r['method']=='merged'))
    predictions=[r for r in read_csv(paths['complete_joint']/'all_method_predictions.csv') if r['method']==best['method']]
    if best['method'] in ('nearest_gene3','ml_diagnostic_gene3'):
        evidence={(r['sample'],r['genus']):r for r in read_csv(paths['complete_joint']/'candidate_evidence.csv')}
        status={r['sample']:r['search_status'] for r in read_csv(paths['complete_joint']/'timings.csv')}
        predictions=[]
        raw={(r['sample'],r['genus']):r['nearest_winner_gene_votes'] for r in read_csv(control/'predictions.csv')}
        for record in read(control/'predictions.json'):
            for genus in record[best['method']]:
                e=evidence[record['sample'],genus]
                predictions.append({'sample':record['sample'],'genus':genus,'family':e['family'],'method':best['method'],
                                    'ml_gene_votes':e['ml_gene_votes'],'nearest_gene_votes':e['nearest_gene_votes'],
                                    'nearest_winner_gene_votes':raw.get((record['sample'],genus),''),'search_status':status[record['sample']]})
    if len(predictions)!=int(best['TP'])+int(best['FP']):
        raise ValueError('Recommended whole-rule predictions disagree with evaluation counts')
    write_csv(root/'recommended_predictions.csv',predictions,
              ['sample','genus','family','method','ml_gene_votes','nearest_gene_votes','nearest_winner_gene_votes','search_status'])
    selected_rows=[r for r in by_sample if r['variant']=='complete_joint' and r['method']==best['method']]
    if len(selected_rows)!=10:
        raise ValueError('Recommendation must account for every selected mixture')
    write_csv(root/'recommended_by_sample.csv',selected_rows)
    with ZipFile(root/'recommended_sets.zip','w',compression=ZIP_DEFLATED) as archive:
        for row in rows['complete_joint']:
            genera=sorted(r['genus'] for r in predictions if r['sample']==row['sample'])
            archive.writestr(row['sample']+'/genera.txt',''.join(g+'\n' for g in genera))
    recommendation={'variant':'complete_joint','method':best['method'],'metrics':best,
                    'selection':'Highest observed micro F1 among all fixed positive-call paths and matched diagnostics; review excluded; ties favor precision, recall, then merged',
                    'scope':'Post-evaluation whole-method development selection on these same ten mixtures; not independent validation',
                    'manual_genus_additions_or_removals':False,'sample_genus_events':len(predictions)}
    dump(root/'RECOMMENDATION.json',recommendation)
    retention=plot(root,all_summary)
    dump(root/'COMPLETION_COMPARISON_AUDIT.json',{'same_ten_inputs_and_parents':True,'same_model_and_training':True,
         'same_QC_alignment_scoring_and_dependencies':True,'same_initial_and_every_routing_pool':True,
         'every_retained_genus_in_complete_joint_search':True,'same_literal_truth':True,
         'human_development_used_previous_evaluation':True,'true_genus_retention_by_stage':retention,
         'stage_order':['initial','after_routing','joint_pool','QC_references','merged'],
         'search_work':stats,'recommendation_is_posthoc':True})
    pct=lambda x:'NA' if x in ('',None) else f'{float(x):.1%}'
    distinct_genera={g for r in by_sample if r['variant']=='complete_joint' and r['method']=='merged'
                    for g in r['expected'].split(';') if g}
    primary=next(r for r in all_summary if r['variant']=='complete_joint' and r['method']=='merged')
    lines=['# Hierarchical genus search: ten in-silico mixtures on TTU','',
           f'Best observed fixed final rule after the complete-pool correction: **{best["method"]}**, '
           f'precision **{pct(best["precision"])}**, recall **{pct(best["recall"])}**, F1 **{float(best["F1"]):.3f}** '
           f'({best["TP"]} TP, {best["FP"]} FP, {best["FN"]} FN; {best["exact_sets"]}/10 exact genus sets). '
           '[Download the recommended genus set](recommended_predictions.csv) or [one text list per mixture](recommended_sets.zip).', '',
           '**Development scope:** two ten-mixture baselines were frozen, predicted, audited and sealed before genus truth was opened. '
           'Their evaluation exposed a stopping-rule defect. The complete-pool correction and the recommendation above use the same '
           'ten mixtures after that evaluation; they are in-sample development results, not a new blind or independent validation. '
           'The genus model, its 0.83 threshold, QC, alignment and final-rule thresholds were unchanged.', '',
           '## What worked and what failed','',
           'Both searches retained **all 60 true sample–genus events through every routing round**. '
           f'The stricter policy removed a median {stats["strict_route"]["median_reduction"]:.1%} of the initial genus pool. '
           'Central Kew representatives and repeated independent-gene support provided useful narrowing.', '',
           'The initial final stage tested only locally resolved small cells when more than 64 genera remained. '
           'That excluded 54 true events in v1 and 40 in strict routing; query fragments from those untested genera '
           'could then match related substitutes. The low end-to-end recall was caused by this final-panel fallback, '
           'rather than loss of the true genera during routing.', '',
           'The correction keeps the exact frozen strict candidate lists and jointly verifies every retained genus '
           'when the pool is at most 512. The largest pool here is 447 genera. Five existing complete-pool final '
           'tests are reused; five partial tests are completed. No expected genera are added, no individual '
           'predictions are edited using truth, and mixture component counts never cap the answer.', '',
           '| Same merged final classifier | TP | FP | FN | Precision | Recall | F1 | Median retained genera |',
           '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for variant in paths:
        r=next(r for r in all_summary if r['variant']==variant and r['method']=='merged')
        lines.append(f'| {variant} | {r["TP"]} | {r["FP"]} | {r["FN"]} | {pct(r["precision"])} | {pct(r["recall"])} | '
                     f'{float(r["F1"]):.3f} | {stats[variant]["median_retained_genera"]:g} |')
    lines+=['','![Genus accuracy and retention through the pipeline](comparison.png)','',
            '## Completed-pool method comparison','',
            '| Fixed rule | TP | FP | FN | Precision | Recall | F1 | Exact sets |',
            '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in all_summary:
        if r['variant']=='complete_joint':
            lines.append(f'| {r["method"]} | {r["TP"]} | {r["FP"]} | {r["FN"]} | {pct(r["precision"])} | {pct(r["recall"])} | '
                         f'{float(r["F1"]):.3f} | {r["exact_sets"]}/10 |')
    lines+=['','The merged rule remains the original primary comparison. Review is a separate follow-up set. '
            'The recommendation compares every fixed positive-call path, including core and auxiliary component paths, '
            'and the matched diagnostics. Including auxiliary in this choice was decided after inspecting these results. '
            'The complete output of one rule is used for every mixture; no genus-specific edits or truth-based '
            'top-N truncation are applied. This observed-F1 choice adds development selection bias.', '',
            'The auxiliary path requires a C-style selected call, gene rescue, a Monte Carlo stability lower bound '
            '>0.5, and at least one ML-supported gene. It is narrower than the original core-or-auxiliary union. '
            'Its stability score is not a calibrated biological presence probability.', '',
            'The nearest_gene3 control uses exactly the same pre-verifier nearest candidates and three-gene quorum '
            'as standalone ML. Replaying its accepted fragments must reproduce the sealed ML set. This diagnostic '
            'was added after the baseline evaluations, before inspecting complete-pool metrics; its prediction '
            'constructor does not read truth. It separates the verifier effect from sampling and reference coverage.', '',
            'The ml_diagnostic_gene3 diagnostic additionally requires three genes with a unique nearest genus '
            'and sibling-relative margin ≥0.1, using the existing C-style diagnostic-gene evidence. Its quorum '
            'and margin were fixed before complete-pool scores were inspected, after the baseline evaluations. '
            'It addresses ambiguous conserved fragments and is eligible only as a disclosed development choice.', '',
            'At the same 54 recovered events and three-gene quorum, the matched nearest-distance control produces '
            '226 false positives, versus 20 with the learned verifier. Adding the three-diagnostic-gene condition '
            'does not change the standalone ML calls in this cohort. This supports a verifier benefit on these '
            'frozen panels; it does not isolate the sampler effect from all other pipeline changes.', '',
            '## Remaining errors','',
            'All six merged false negatives first fail final-reference coverage, rather than routing or the final '
            'classifier: Gynostemma occurs in three missed events with a best checked reference of 11 accepted '
            'target loci, and Agropyropsis in three with 19. Both fall below the fixed 30-locus reference gate. '
            'The [coverage audit](REFERENCE_COVERAGE_DIAGNOSIS.json) records the source QC hashes and accepted '
            'locus counts. No coverage threshold was lowered using these known misses.', '',
            'The recommended set still has 16 false positives and only one of ten mixtures has an exact genus set. '
            'Treat it as a candidate-screening result. A future generic low-coverage reference tier and any '
            'stronger false-positive filter should be fixed before testing a fresh held-out mixture cohort.', '',
            '## Data, time and limitations','',
            'The cohort contains four 3-component, three 6-component and three 10-component mixtures, selected '
            f'by sample-ID hashes from the completed 100-mixture dataset. They represent {len(distinct_genera)} distinct source genera '
            'and 60 sample–genus presence events. Reused assemblies cover 50 target genes. '
            'Search parents come from existing raw Z, merged and review family evidence. Target-locus QC requires '
            'at least 30 accepted genes; up to two species represent each final genus. Unavailable genera stay '
            'explicit and count as misses.', '',
            'Kew planning distances count topology edges because the supplied tree has no branch lengths. '
            'Each coarse round builds fresh exon trees against reduced references. The stricter policy uses '
            'five routing genes and one nearby-genus guard, with 32 initial cells and 64 after stagnation. '
            'Final genus features use full patristic distance profiles and a separately reference-trained verifier.', '',
            'Reference projection is a substantial runtime cost. Completion reuses only hashed, QC-passing '
            'fixed-coordinate references and caches negative QC checks. New references receive the same QC. '
            'Completion timings include both the reused strict pipeline and the incremental joint work, including '
            'the already-spent partial final run. Incremental times and reused-final flags are reported separately. '
            'These warm-cache observations do not establish speedup against a matched exhaustive genus run. '
            'Queueing, audit and Slurm overhead are separate from inference phases.', '',
            'The final panels can contain the source species; this is a closed-reference experiment. '
            'Mixture genera were not required to be absent from model fitting. The model was fitted without mixture labels, '
            'but the cohort and reference library have been used '
            'in previous work. Results do not establish accuracy on new species, new mixtures, absent-reference '
            'genera, negative-only samples, species-level presence or abundance.', '',
            'Truth is the first word of Table S1 scientific names, using source IDs in the fixed mixture IDs. '
            'All used A/B/H cells agree between openpyxl and independent ZIP/XML reads; raw method text lists '
            'independently reproduce TP/FP/FN. Source family annotations for IDs 1 and 25 have known discrepancies, '
            'but no genus substitutions are made. All 60 expected genera were already in the initial pools.', '',
            '## Files','',
            '- [Complete-pool correction: per-sample results, timing and QC](complete_joint/REPORT.md)',
            '- [Original strict policy](strict_route/REPORT.md) and [v1](v1/REPORT.md)',
            '- [Original two-policy comparison](BASELINE_COMPARISON.md)',
            '- [All variant and method scores](all_variant_summary.csv)',
            '- [Recommended whole-rule selection record](RECOMMENDATION.json)',
            '- [Recommended per-sample accuracy and error lists](recommended_by_sample.csv)',
            '- [Source/model/routing reuse audit](COMPLETION_COMPARISON_AUDIT.json)',
            '- [Shared trained genus model](shared_model/model.joblib)',
            '- [Post-evaluation correction plan](../../POSTHOC_COMPLETION_PLAN.md)','',
            'All source payloads, job records, hash audits and raw compact results are in artifacts/. '
            'The original workbook and DNA read data are not redistributed. '
            'Branch: `experiment/adaptive-genus-search-20261017`; production defaults and previous family results are unchanged.','']
    (root/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
    print(recommendation)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True,type=Path)
    final_report(parser.parse_args().root)
