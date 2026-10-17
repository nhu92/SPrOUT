"""Render the sealed ten-mixture evaluation into a report and scientific figures."""
import argparse
from collections import Counter
from pathlib import Path
from statistics import median

from .evaluation import read, read_csv, verify_export


def pct(value):
    return 'NA' if value is None or value=='' else f'{float(value):.1%}'


def figures(output, selection, timings, by_sample, searches):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','savefig.dpi':220})
    colors=plt.cm.tab10(np.arange(10))
    fig, axes=plt.subplots(1,2,figsize=(11,4.8),constrained_layout=True)
    for i,row in enumerate(selection):
        search=next(s for s in searches if s['sample']==row['sample'])
        pool=[search['initial_genera']]+[r['retained_genera'] for r in search['rounds']]
        axes[0].plot(range(len(pool)),pool,'o-',color=colors[i],label=f'S{i+1}')
        axes[1].plot(range(len(pool)),[n/pool[0] if pool[0] else 0 for n in pool], 'o-',color=colors[i])
    axes[0].set(xlabel='Completed routing rounds',ylabel='Remaining search genera',yscale='log',title='A  Search pool size')
    axes[0].axhline(64,color='0.4',ls='--',lw=1,label='Joint-panel budget')
    axes[1].set(xlabel='Completed routing rounds',ylabel='Fraction of initial genera retained',ylim=(-.02,1.03),title='B  Relative reduction')
    axes[1].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    for ax in axes:
        ax.grid(alpha=.2)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    axes[0].legend(ncol=4,fontsize=8,loc='upper right')
    for suffix in ('png','svg'):
        fig.savefig(output/('search_reduction.'+suffix))
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),constrained_layout=True)
    primary=[next(r for r in by_sample if r['sample']==s['sample'] and r['method']=='merged') for s in selection]
    x=np.arange(10);width=.25
    for offset,field,color,label in [(-1,'TP','#228833','Recovered'),(0,'FP','#cc6677','False positive'),(1,'FN','#ddaa33','Missed')]:
        axes[0].bar(x+offset*width,[int(r[field]) for r in primary],width,color=color,label=label)
    axes[0].set(ylabel='Sample–genus events',title='A  Frozen merged classifier')
    axes[0].legend(fontsize=8)
    axes[0].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    axes[1].bar(x,[float(r['seconds'])/60 for r in timings],color='#4477aa')
    axes[1].set(ylabel='Elapsed inference minutes',title='B  Per-mixture elapsed time')
    for ax in axes:
        ax.set_xticks(x,[f'S{i+1}' for i in range(10)])
        ax.set_xlabel('Mixture (mapping in report)')
        ax.grid(axis='y',alpha=.2)
        ax.set_axisbelow(True)
    for suffix in ('png','svg'):
        fig.savefig(output/('accuracy_runtime.'+suffix))
    plt.close(fig)
    # Matplotlib emits trailing spaces in SVG style blocks; preserve clean diffs.
    for path in output.glob('*.svg'):
        path.write_text('\n'.join(line.rstrip() for line in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')


def report(output):
    output=Path(output);verify_export(output)
    summary=read_csv(output/'summary.csv');by_sample=read_csv(output/'by_sample.csv')
    timings=read_csv(output/'timings.csv');rounds=read_csv(output/'rounds.csv')
    selection=read(output/'selection.json')['selected'];searches=read(output/'SEARCHES.json')
    training=read(output/'model/TRAINING.json');audit=read(output/'EVALUATION_AUDIT.json')
    strict = 'strict-routing' in read(output/'PROTOCOL.json')['policy']['schema']
    primary=next(r for r in summary if r['method']=='merged')
    reductions=[1-float(t['remaining_search_genera'])/float(t['initial_genera']) for t in timings if float(t['initial_genera'])]
    seconds=[float(r['seconds']) for r in timings]
    refined=sum(r['search_status']=='refined' for r in timings)
    bounded=sum(int(r['remaining_search_genera'])<=64 for r in timings)
    figures(output,selection,timings,by_sample,searches)
    lines=['# Ten in-silico mixtures: adaptive phylogenetic genus search','',
           f'All ten TTU pipelines completed and passed frozen input/model/code/product checks. '
           f'The predeclared merged genus method recovered {primary["TP"]} true sample–genus events, '
           f'with {primary["FP"]} false positives and {primary["FN"]} missed events: '
           f'precision **{pct(primary["precision"])}**, recall **{pct(primary["recall"])}**, F1 **{float(primary["F1"]):.3f}**.', '',
           f'The search removed a median {median(reductions):.1%} of the initial recoverable genus pool. '
           f'{bounded}/10 searches reached the ≤64-genus search target. '
           f'{refined}/10 searches reached a fully QC-represented final joint panel; '
           f'{10-refined}/10 retained explicit unresolved regions or final reference omissions. '
           f'Per-mixture elapsed inference was {min(seconds)/60:.1f}–{max(seconds)/60:.1f} minutes '
           f'(median {median(seconds)/60:.1f}); shared reference-only model training took {training["seconds"]/60:.1f} minutes.', '',
           'This is a retrospective method experiment on ten mixtures from the previously analyzed 100. '
           'The 4/3/3 selection across 3-, 6- and 10-component strata was fixed by sample-ID hashes, '
           'without genus truth or prior success scores. It is not an independent biological validation. '
           'Mixture component counts never determine the number of predicted genera.', '',
           '## Search behavior','',
           'The existing family search parents are the union of raw Z, frozen merged and review families. '
           'Inside those parents, partition the recoverable Kew tree into 32 neighborhoods, subdivide '
           'each into representative regions, and choose central genera using genus-balanced minimax/mean '
           'distance to a farthest-point landmark cover. The supplied Kew tree has no branch lengths; '
           'these planning distances count topology edges, not evolutionary time. All recoverable placements '
           'of polyphyletic genera remain eligible.', '',
           'Representative sequences must pass MrKnowItAll QC on at least 30 of the actual 50 target genes. '
           'Download only proposed recoveries, replace failed representatives locally, and preserve uncovered '
           'regions as unresolved. Cache QC and projected reference coordinates per sample. '
           'At each round, align the existing query exons to the reduced panel and build fresh full exon trees. '
           'Two independent genes retain a region; a close second route needs three genes. '
           'Nearby guards add four genera per supported anchor, within the previous active pool. '
           'The anchor genus labels its region until the final joint test.', '',
           'Repartition and repeat. A stalled round permits 64 cells; two consecutive ≥95% retention rounds '
           'or eight total rounds stop the search as unresolved. At ≤64 surviving genera, jointly test up to '
           'two QC-passing species per genus. Larger unresolved searches can test fully local neighborhoods '
           'while retaining all other genera explicitly as unresolved. No biological top-N prediction cap is used.', '',
           '![Search pool trajectories](search_reduction.png)','',
           '| ID | Mixture | Search genera: initial → retained | Final tested genera / species | Rounds | Minutes | Status | TP / FP / FN |',
           '| --- | --- | ---: | ---: | ---: | ---: | --- | ---: |']
    for i,row in enumerate(selection):
        t=next(r for r in timings if r['sample']==row['sample'])
        p=next(r for r in by_sample if r['sample']==row['sample'] and r['method']=='merged')
        lines.append(f'| S{i+1} | {row["sample"]} | {t["initial_genera"]} → {t["remaining_search_genera"]} | '
                     f'{t["final_reference_genera"]} / {t["final_reference_species"]} | {t["rounds"]} | '
                     f'{float(t["seconds"])/60:.1f} | {t["search_status"]} | {p["TP"]} / {p["FP"]} / {p["FN"]} |')
    lines+=['',f'Across {len(rounds)} routing rounds, median per-round reduction was '
            f'{median(float(r["reduction_fraction"]) for r in rounds):.1%}. '
            'Full per-round panels, counts, evidence, QC and timing are in the sealed batch archive and rounds.csv. '
            'The pilot and training costs are separate from the ten sample times. '
            f'Slurm used eight CPUs per mixture with up to {6 if strict else 4} mixtures of this policy concurrently; '
            'the strict-policy throttle increased from two to four to six as earlier work finished. '
            'Reported inference times exclude queueing, prior assembly, initial/final audits and Slurm setup/cleanup. '
            'Allocated job elapsed times are recorded separately in Slurm accounting. '
            'Fetch/QC/alignment durations overlap worker tasks and are not CPU-time totals.', '',
            'Reference-count reduction measures avoided reference work. There is no exhaustive genus run '
            'on these same ten samples, so these results do not establish a measured speedup over the old '
            'pipeline. The previous progressive ten-mixture cohort differs from this cohort.', '',
            '## Genus predictions','',
            'Primary output: [merged_predictions.csv](merged_predictions.csv). '
            'The [complete method sets](prediction_sets.zip), [candidate evidence](candidate_evidence.csv) '
            'and [unresolved search regions](search_regions.csv) are separate exports.', '',
            '| Method | TP | FP | FN | Precision | Recall | F1 | Exact sets |',
            '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in summary:
        lines.append(f'| {r["method"]} | {r["TP"]} | {r["FP"]} | {r["FN"]} | {pct(r["precision"])} | '
                     f'{pct(r["recall"])} | {float(r["F1"]):.3f} | {r["exact_sets"]}/10 |')
    lines+=['','Review is an unmerged follow-up set, not another primary call set. '
            'C-style uses a new genus-balanced distance contrast and the previous gene-bootstrap confirmation/rescue formulas. '
            'Merged keeps C-style calls with at least two ML-supported genes, or the predeclared stable auxiliary rescue '
            'with one ML-supported gene. Balanced Z is a within-panel standardized contrast; it is not the old raw family Z.', '',
            '![Merged genus accuracy and elapsed time](accuracy_runtime.png)', '',
            '## Where true genera were lost','',
            '| First loss / recovery stage | Sample–genus events |','| --- | ---: |']
    for stage,n in sorted(audit['first_losses'].items()):
        lines.append(f'| {stage} | {n} |')
    lines+=['','[truth_survival.csv](truth_survival.csv) records every true genus surviving each round; '
            '[first_losses.csv](first_losses.csv) identifies the first exclusion. All ten mixtures contribute to recall, '
            'including absent parents, unavailable references, unresolved searches and final rejections. '
            'The original Table S1 is the primary literal truth. Older family reports documented source-identity '
            'discrepancies for IDs 1 and 25; no genus substitution was inferred from these predictions.', '',
            '## Model and reproducibility','',
            f'The new genus verifier used {training["reference_species"]} Kew species across '
            f'{training["reference_genera"]} genera and {training["reference_only_trees"]} reference-only gene trees. '
            'Entire genera are partitioned for fitting, calibration and testing. Known-genus and genus-absent '
            f'pseudoqueries fixed the threshold at {training["threshold"]:.2f} before mixture prediction. '
            f'Reference-held-out correct known-genus acceptance was {training["test_known_correct_acceptance"]:.1%}; '
            f'genus-absent rejection was {training["test_genus_absent_rejection"]:.1%}. '
            'These are reference fragment checks, not ten-mixture accuracy or biological presence probabilities. '
            'The reference library came from previously constructed panels and is not an independent library.', '',
            'The final query workflow also changes reference-column trimming and alignment projection. '
            'Accuracy/runtime differences cannot be attributed to the sampler alone. Sparse loci and cross-genus '
            'similarity can defeat routing; unresolved regions are visible instead of forced into representative labels. '
            'Species-level presence and abundance are not evaluated.', '',
            'Inference implementation: `1fbf38d`; evaluation plan/export implementation: `502e1e1`. '
            'Branch: `experiment/adaptive-genus-search-20261017`. '
            'The previous pipeline defaults and family results are unchanged. '
            'No mixture-truth-based threshold tuning was performed in this frozen experiment.', '',
            'The [protocol](PROTOCOL.json), [prediction seal](PREDICTION_EXPORT.json), '
            '[independent evaluation audit](EVALUATION_AUDIT.json), [training record](model/TRAINING.json), '
            'runtime source snapshot, trained model, job records and batch archive preserve the run. '
            'The model training used its separately archived source snapshot; technical pilot failure records are retained. '
            'The original workbook and raw reads are not redistributed.', '',
            f'Truth source: original 30-source Table S1, SHA-256 `{audit["workbook_sha256"]}`. '
            'Every used source-ID, genus/name and family cell was checked with both openpyxl and ZIP/XML. '
            'Per-method plain text prediction files independently reproduce TP/FP/FN. '
            'The truth workbook was opened only after all ten outputs were complete, downloaded and sealed locally.', '']
    (output/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    report(parser.parse_args().output)
