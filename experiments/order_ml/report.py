"""Human-readable comparison of all prespecified order operating points."""
import argparse
import csv
from pathlib import Path
import statistics

from .common import read


def fmt(row, key, spec):
    value = row[key]
    return format(float(value), spec) if value not in ('', None) else '—'


def report(output):
    output = Path(output)
    with (output/'summary.csv').open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    training = read(output/'TRAINING.json'); timings = read(output/'TIMINGS.json')
    all100 = [r for r in rows if r['subset'] == 'all100']
    z = next(r for r in all100 if r['method'] == 'z_only')
    primary = next(r for r in all100 if r['method'] == 'z_plus_ml_specificity')
    metric = lambda row, key: float(row[key])
    if metric(primary, 'F1') > metric(z, 'F1') and metric(primary, 'precision') >= metric(z, 'precision'):
        verdict = 'The prespecified specificity union improved F1 without reducing observed precision.'
    elif metric(primary, 'F1') > metric(z, 'F1'):
        verdict = 'The prespecified specificity union improved F1, with an observed precision tradeoff.'
    elif primary['TP'] == z['TP'] and primary['FP'] == z['FP']:
        verdict = 'The prespecified specificity union did not change the baseline calls.'
    else:
        verdict = 'The prespecified specificity union did not improve F1; retain the original Z output for formal order prediction.'
    text = f'''# Reference-trained order ML on the frozen 100 mixtures

{verdict}

This is the order-only experiment on `experiment/order-ml100-20261023`.
Primary new method: `z_plus_ml_specificity`, fixed before mixture inference.
Existing family/genus outputs and their broad search admissions were not rerun or changed.

| Method | TP | FP | FN | Precision | Recall | F1 | Exact sets /100 | Rescued Z FN | New FP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
'''
    for r in all100:
        text += f"| {r['method']} | {r['TP']} | {r['FP']} | {r['FN']} | {fmt(r, 'precision', '.3%')} | {fmt(r, 'recall', '.3%')} | {fmt(r, 'F1', '.4f')} | {r['exact_sets']} | {r['rescued_z_false_negatives']} | {r['new_false_positives']} |\n"
    text += '\nThe union cannot lose existing Z true positives. ML-only rows can remove old calls. '
    text += 'A rescue is a sample–order presence event; repeated species from the same order count once. '
    text += 'Precision measures correctness of positive calls, whereas recall measures recovery of true orders.\n\n'
    text += '| Subset | Method | Samples | Precision | Recall | F1 | Exact sets |\n|---|---|---:|---:|---:|---:|---:|\n'
    for r in rows:
        if r['subset'] != 'all100' and r['method'] in ('z_only', 'z_plus_ml_specificity'):
            text += f"| {r['subset']} | {r['method']} | {r['samples']} | {float(r['precision']):.3%} | {float(r['recall']):.3%} | {float(r['F1']):.4f} | {r['exact_sets']} |\n"
    text += '\n## Reference fitting and held-out checks\n\n'
    split = training['split']
    text += f"Fresh reference-only trees: {len(training['trees'])}. Whole-order split: {len(split['fit'])} fitting, {len(split['calibration'])} calibration, {len(split['test'])} test, {len(split['anchor_only'])} anchor-only orders. "
    text += f"Fitting: {training['fit_species']} species, {training['fit_rows']} rows including {training['fit_order_absent_negatives']} nearest-candidate negatives after removing the entire true order. "
    text += f"Held-out calibration/test species: {training['calibration_species']}/{training['test_species']}.\n\n"
    text += '| Fixed point | Threshold | Calibration known recall | Calibration absent rejection | Test known recall | Test absent rejection |\n|---|---:|---:|---:|---:|---:|\n'
    for key in ('balanced', 'specificity'):
        cal = training['thresholds'][key]; test = training['reference_test'][key]
        text += f"| {key} | {cal['threshold']:.2f} | {cal['known_correct_acceptance']:.1%} | {cal['absent_species_rejection']:.1%} | {test['known_correct_acceptance']:.1%} | {test['absent_species_rejection']:.1%} |\n"
    text += '\nThresholds use only the reference calibration orders. Both reference and mixture acceptance require at least three distinct genes for the same order. '
    text += 'Reference absent rejection is the fraction of species for which no wrong order reaches that quorum. '
    text += 'It is not a per-fragment average and is not a guarantee of mixture precision. '
    text += 'Reference profiles use seven distance features and a species-weighted standardized logistic verifier; nearest-order selection precedes verification. '
    text += 'The model score is a discriminative verifier score, not a calibrated probability of sample presence.\n\n'
    seconds = [r['seconds'] for r in timings]
    text += '## Computation and provenance\n\n'
    text += f"Reference tree building: {training['tree_seconds']/60:.2f} minutes; entire fitting/calibration/test stage: {training['total_seconds']/60:.2f} minutes. "
    text += f"Per-mixture inference/audit wall time: median {statistics.median(seconds):.2f}s, range {min(seconds):.2f}–{max(seconds):.2f}s; summed sample wall times {sum(seconds)/3600:.3f} hours. "
    text += 'Sample wall times overlap; they must not be added to claim elapsed cohort time. Slurm accounting records queue/start/end times.\n\n'
    text += f"Checkpoint SHA256: `{training['model_sha256']}`. "
    text += 'Input tree hashes were checked against the original sealed run and rechecked before/after inference; all raw exon scores reproduce every output set. '
    text += 'All 100 outputs were audited and sealed before this evaluation. Literal order labels are Table S1 column G, independently checked with openpyxl and workbook ZIP/XML locally; the workbook was not uploaded. '
    text += 'See `PLAN.md`, `PROTOCOL.json`, `MODEL_FROZEN.json`, `PREDICTION_SEAL.json`, `EVALUATION_AUDIT.json`, and `artifacts/` for source, checkpoint, raw results and provenance.\n\n'
    text += '## Interpretation limits\n\n'
    text += 'This is a retrospective experiment on a known cohort whose earlier truth had already been inspected. '
    text += 'The prior genus-development ten and other ninety are reported separately; neither is a new independent biological cohort. '
    text += 'Whole-order holdout tests verifier transfer but shares a reference-only phylogenetic scaffold; source species may be in the reference panel. '
    text += 'Filtering unnamed species and collapsing accessions changes the ML anchor library relative to Z; therefore this comparison tests the entire fixed ML path, not an isolated causal effect of logistic regression. '
    text += 'A single known/absent reference pseudoquery has fewer opportunities for false quorums than a complex mixed sample, so high reference specificity can fail to transfer. '
    text += 'The two thresholds were fixed before mixture scoring; no best-performing mixture threshold is selected here. '
    text += 'These results do not establish family/genus accuracy after narrowing their candidate pools.\n'
    (output/'REPORT.md').write_text(text, encoding='utf-8')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    report(p.parse_args().output)
