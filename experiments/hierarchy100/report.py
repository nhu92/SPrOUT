"""Report the fixed primary order/family/genus rules without reselection."""
import argparse
from collections import defaultdict
from pathlib import Path
from statistics import median

from experiments.adaptive_genus.evaluation import read,read_csv,write_csv
from experiments.hierarchy100.evaluation import verify_seal


PRIMARY={'order':'admission','family':'merged','genus':'auxiliary'}


def pct(value):
    return 'NA' if value in ('',None) else f'{float(value):.1%}'


def report(batch,output):
    verify_seal(output,batch)
    summary=read_csv(output/'summary.csv');by_sample=read_csv(output/'by_sample.csv')
    rows=read(output/'selection.json')['selected'];timing=read(output/'TIMINGS.json')
    phases=defaultdict(list);sample_seconds={}
    for row in timing:
        sample_seconds[row['sample']]=row['inference_seconds']
        for r in row['phases']:phases[r['phase']].append(float(r['seconds']))
    searches={r['sample']:read(batch/'samples'/r['project']/'03_genus/result/SEARCH.json') for r in rows}
    search_rows=[{'sample':s,'initial_genera':r['initial_genera'],'remaining_genera':r['remaining_search_genera'],
         'final_reference_genera':r['final'].get('reference_genera',0),'status':r['search_status'],
         'rounds':len(r['rounds']),'genus_seconds':r['seconds'],
         'fetch_seconds':r['fetch_seconds'],'QC_seconds':r['qc_seconds'],
         'reference_projection_seconds':r['reference_alignment_seconds']} for s,r in searches.items()]
    write_csv(output/'genus_search_summary.csv',search_rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    x=np.arange(3);width=.24;ranks=list(PRIMARY)
    for delta,key,color in [(-1,'precision','#4477aa'),(0,'recall','#228833'),(1,'F1','#aa3377')]:
        values=[next(r for r in summary if r['subset']=='all100' and r['rank']==rank and r['method']==PRIMARY[rank])[key] for rank in ranks]
        axes[0].bar(x+delta*width,[float(v) if v!='' else 0 for v in values],width,label=key,color=color)
    axes[0].set_xticks(x,['Order admission','Family merged','Genus auxiliary']);axes[0].set(ylim=(0,1.17),ylabel='Micro metric',title='A  Fixed primary rules, all 100')
    axes[0].legend(ncol=3,fontsize=8)
    for row in rows:
        r=searches[row['sample']];values=[r['initial_genera']]+[p['retained_genera'] for p in r['rounds']]
        axes[1].plot(range(len(values)),[v/max(1,values[0]) for v in values],color='#aa3377' if row['genus_development_ten'] else '#4477aa',alpha=.35,lw=.9)
    axes[1].set(xlabel='Completed routing rounds',ylabel='Fraction of initial genera',ylim=(-.02,1.05),title='B  Genus candidate pool trajectories')
    for ax in axes:
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    for ext in ('png','svg'):fig.savefig(output/('validation.'+ext),dpi=220)
    plt.close(fig)
    svg=output/'validation.svg';svg.write_text('\n'.join(s.rstrip() for s in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    lines=['# 100-mixture fresh order-to-genus validation on TTU','',
         'All 100 samples completed fresh order/family/genus inference and source/model/input/product checks. '
         'Only query assemblies and reference resources were reused. Fixed primary rules are order admission, '
         'family merged and genus auxiliary; the best-performing rule is not reselected using this run.', '',
         '| Subset | Rank / fixed rule | Samples | TP | FP | FN | Precision | Recall | F1 | Exact sets |',
         '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in summary:
        if r['subset'] in ('all100','genus_development10','other90') and r['method']==PRIMARY.get(r['rank']):
            lines.append(f'| {r["subset"]} | {r["rank"]} / {r["method"]} | {r["samples"]} | {r["TP"]} | {r["FP"]} | {r["FN"]} | {pct(r["precision"])} | {pct(r["recall"])} | {float(r["F1"]):.3f} | {r["exact_sets"]}/{r["samples"]} |')
    lines+=['','![Primary accuracy and genus pool trajectories](validation.png)','',
         'Pink trajectories are the previous ten genus-development samples; blue trajectories are the other ninety. '
         'The cohort contains 40 three-component, 40 six-component and 20 ten-component mixtures. '
         'All samples and their source species were analyzed previously. The other ninety are additional genus '
         'evaluation mixtures, not independent biological samples or unseen-source-species validation.', '',
         '## Search and errors','']
    reduction=[1-r['remaining_search_genera']/r['initial_genera'] for r in searches.values() if r['initial_genera']]
    if reduction:
        lines.append(f'Median within-genus-search candidate reduction: {median(reduction):.1%}. '
                     f'{sum(r["search_status"]=="refined" for r in searches.values())}/100 final panels are fully QC represented. '
                     'Missing parents, missing references and unresolved searches remain in recall denominators. '
                     'The normal routing target is 64 genera and the final complete-pool limit is 512; larger pools abstain.')
    lines+=['','| First genus loss / recovery | Occurrences |','|---|---:|']
    for stage,count in read(output/'EVALUATION_AUDIT.json')['genus_first_losses'].items():lines.append(f'| {stage} | {count} |')
    unresolved=[r for r in rows if searches[r['sample']]['unresolved_genera']>0]
    lines+=['',f'{len(unresolved)}/100 searches finish with an unresolved final pool. '
         'The fully QC-represented count above requires a reference for every candidate genus; '
         'a missing reference for an irrelevant candidate also marks a panel as partial. '
         'This status is separate from successful computational completion.', '',
         '| Genus auxiliary by mixture size | Samples | TP | FP | FN | Precision | Recall | F1 |',
         '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in summary:
        if r['rank']=='genus' and r['method']=='auxiliary' and r['subset'].startswith('components'):
            lines.append(f'| {r["subset"]} | {r["samples"]} | {r["TP"]} | {r["FP"]} | {r["FN"]} | {pct(r["precision"])} | {pct(r["recall"])} | {float(r["F1"]):.3f} |')
    lines+=['','Per-sample truth and all error lists are in [by_sample.csv](by_sample.csv). '
         'Every fixed alternate method and mixture-size stratum is in [summary.csv](summary.csv). '
         'Primary outputs are in [primary_predictions.csv](primary_predictions.csv). '
         'Order calls and the broad family Z/merged/review union are search admissions; a genus can be '
         'supported through a review family rather than the final family merged list. '
         'No ancestral closure or truth-based manual replacement is applied.', '',
         '## Time','',f'Fresh inference per sample: {min(sample_seconds.values())/60:.2f}–{max(sample_seconds.values())/60:.2f} minutes; '
         f'median {median(sample_seconds.values())/60:.2f} minutes. This includes per-genus in-task audits. '
         'Final parent-product auditing, preparation, final collection, queueing and prior assembly are separate.', '',
         '| Phase | Median minutes | Sum of sample-phase minutes |','|---|---:|---:|']
    for phase,values in phases.items():lines.append(f'| {phase} | {median(values)/60:.2f} | {sum(values)/60:.2f} |')
    lines+=['','Tasks execute concurrently according to the scheduler; phase sums are not cohort wall time or CPU time. '
         'No workflow array concurrency throttle is applied. Training is not repeated. '
         'Cached Kew DNA and pretrained checkpoints make this a warm-cache benchmark. '
         'Fine genus rounds are in [genus_rounds.csv](genus_rounds.csv); QC/projection/final coverage summaries '
         'are in [genus_search_summary.csv](genus_search_summary.csv).', '',
         '## Interpretation and provenance','',
         'The plan, policies and models were fixed before this inference batch. The genus auxiliary choice '
         'and complete-pool stopping correction were selected in the previous ten-mixture development. '
         'The family and genus verifier thresholds remain 0.94 and 0.83. Reference coverage remains 30 of 50 '
         'target genes for genus representatives. Source species are allowed in reference panels. '
         'Scores and selection stability are not calibrated biological presence probabilities.', '',
         'Truth uses literal Table S1 source IDs, G order, H family and the first word of B scientific names. '
         'All used cells are checked with openpyxl and independent ZIP/XML. Known family-name discrepancies '
         'for source IDs 1 and 25 are retained in this literal view; genus names are not substituted. '
         'The workbook stays local and unchanged. Outputs are sealed and raw text lists checked before scoring.', '',
         'Branch: `experiment/order-to-genus100-20261022`. See [PLAN.md](PLAN.md), '
         '[PROTOCOL.json](PROTOCOL.json), [PREDICTION_SEAL.json](PREDICTION_SEAL.json) and '
         '[EVALUATION_AUDIT.json](EVALUATION_AUDIT.json).','']
    (output/'REPORT.md').write_text('\n'.join(s.rstrip() for s in lines).rstrip()+'\n',encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--batch',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path);args=parser.parse_args();report(args.batch,args.output)
