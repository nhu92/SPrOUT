"""Local-only audit/seal/evaluation. Importable without NumPy, Bio or sklearn."""
import argparse
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import shutil

from .common import METHODS, read, dump, sha, combine
from experiments.adaptive_genus.evaluation import counts, names, write_csv


def verify_local_mapping(batch, mapping):
    remote = read(batch/'REQUEST.json')['root']
    for name, digest in mapping.items():
        path = Path(name)
        relative = path.relative_to(remote)
        if sha(batch/relative) != digest:
            raise ValueError('Changed downloaded frozen file '+str(relative))


def complete(batch):
    batch = Path(batch); audit = read(batch/'AUDIT.json')
    if not (batch/'PREDICTIONS_COMPLETE').exists() or not audit.get('all_100_complete') or not audit.get('code_model_inputs_outputs_verified'):
        raise ValueError('All hundred samples and immutable-input audits must finish')
    selected = read(batch/'selection.json')['selected']; calls = read(batch/'predictions.json')
    by_id = {r['sample']: r for r in calls}
    if len(selected) != 100 or len(calls) != 100 or len(by_id) != 100 or set(by_id) != {r['sample'] for r in selected}:
        raise ValueError('Full cohort prediction mismatch')
    if len({r['project'] for r in selected}) != 100:
        raise ValueError('Duplicate project directory')
    protocol = read(batch/'PROTOCOL.json'); frozen = read(batch/'MODEL_FROZEN.json')
    if sha(batch/'PLAN.md') != protocol['plan_sha256'] or sha(batch/'PLAN.md') != sha(Path(__file__).with_name('PLAN.md')):
        raise ValueError('Frozen/local plan mismatch')
    if sha(batch/'selection.json') != protocol['selection_sha256']:
        raise ValueError('Cohort manifest changed')
    if sha(batch/'PROTOCOL.json') != frozen['protocol_sha256']:
        raise ValueError('Frozen protocol changed')
    verify_local_mapping(batch, frozen['sha256'])
    if sha(batch/'model/model.joblib') != frozen['model_sha256']:
        raise ValueError('Checkpoint mismatch')
    thresholds = read(batch/'model/THRESHOLDS.json')
    for row in selected:
        directory = batch/'samples'/row['project']
        if (directory/'FAILED').exists() or not (directory/'SUCCESS').exists():
            raise ValueError('Incomplete sample '+row['sample'])
        success = read(directory/'SUCCESS'); verify_local_mapping(batch, success['sha256'])
        sample = by_id[row['sample']]
        if read(directory/'predictions.json') != sample:
            raise ValueError('Combined/per-sample mismatch')
        recomputed = combine(row['z_only'], read(directory/'HITS.json'), thresholds)
        if sample['methods'] != recomputed:
            raise ValueError('Raw scores do not reproduce frozen prediction sets')
        for method in METHODS:
            values = sample['methods'][method]
            if len(values) != len(set(values)) or set(values) != names(directory/(method+'.txt')):
                raise ValueError('Raw text prediction mismatch')
        sample_audit = read(directory/'AUDIT.json')
        if sample_audit['model_sha256'] != frozen['model_sha256'] or sample_audit['model_freeze_sha256'] != sha(batch/'MODEL_FROZEN.json'):
            raise ValueError('Per-sample model freeze mismatch')
    return selected, by_id


def verify_prior_baseline(batch, original):
    """Link each reused tree/baseline to the already sealed original run."""
    old_seal = read(original/'PREDICTION_SEAL.json')['source_batch_sha256']
    selected = read(batch/'selection.json')['selected']
    for row in selected:
        prior_relative = 'samples/'+row['project']+'/ORDER_FAMILY_PRODUCTS.json'
        if old_seal[prior_relative] != row['prior_products_sha256']:
            raise ValueError('Reused-tree manifest differs from original seal')
        for file, key in (('native_z_candidates.txt', 'baseline_sha256'), ('ORDER.json', 'order_json_sha256')):
            relative = 'samples/'+row['project']+'/01_order/'+file
            if old_seal[relative] != row[key]:
                raise ValueError('Baseline differs from original seal: '+relative)
    return True


def seal(batch, output, original):
    batch, output, original = map(Path, (batch, output, original))
    rows, calls = complete(batch); verify_prior_baseline(batch, original)
    # The exact source payload is retained beside the downloaded batch. Its
    # bytes must match every inference module frozen on TTU before fitting.
    import hashlib
    import tarfile
    protocol = read(batch/'PROTOCOL.json')
    with tarfile.open(batch.parent/'payload.tar.gz') as payload:
        for absolute, digest in protocol['code_sha256'].items():
            relative = Path(absolute).relative_to(protocol['root']).as_posix()
            content = payload.extractfile(relative).read()
            if hashlib.sha256(content).hexdigest() != digest:
                raise ValueError('Frozen code differs from retained payload: '+relative)
    if output.exists():
        raise FileExistsError('Preserve earlier export; choose a new directory')
    output.mkdir(parents=True)
    (output/'.gitattributes').write_text('* -text whitespace=cr-at-eol\n')
    files = ('PLAN.md', 'REQUEST.json', 'PROTOCOL.json', 'selection.json', 'MODEL_FROZEN.json',
             'predictions.json', 'TIMINGS.json', 'AUDIT.json', 'PREDICTIONS_COMPLETE')
    for name in files:
        shutil.copy2(batch/name, output/name)
    for name in ('TRAINING.json', 'THRESHOLDS.json', 'CATALOG.json'):
        shutil.copy2(batch/'model'/name, output/name)
    export = [{'sample': r['sample'], 'method': method, 'order': order}
              for r in rows for method in METHODS for order in calls[r['sample']]['methods'][method]]
    write_csv(output/'all_method_predictions.csv', export)
    dump(output/'PREDICTION_SEAL.json', {'sealed_utc': datetime.now(timezone.utc).isoformat(),
         'all100_complete_before_truth_evaluation': True, 'truth_read_by_exporter': False,
         'retrospective_known_cohort': True, 'baseline_matches_original_seal': True,
         'source_sha256': {p.relative_to(batch).as_posix(): sha(p) for p in batch.rglob('*') if p.is_file()},
         'artifacts_sha256': {p.name: sha(p) for p in output.iterdir() if p.is_file()},
         'evaluation_source_sha256': sha(Path(__file__)),
         'local_plan_sha256': sha(Path(__file__).with_name('PLAN.md'))})


def verify_seal(batch, output):
    record = read(output/'PREDICTION_SEAL.json')
    for name, digest in record['source_sha256'].items():
        if sha(batch/name) != digest:
            raise ValueError('Changed sealed source '+name)
    for name, digest in record['artifacts_sha256'].items():
        if sha(output/name) != digest:
            raise ValueError('Changed sealed export '+name)
    if sha(Path(__file__)) != record['evaluation_source_sha256']:
        raise ValueError('Changed evaluator after seal')
    if sha(Path(__file__).with_name('PLAN.md')) != record['local_plan_sha256']:
        raise ValueError('Changed local evaluation plan')


def evaluate(batch, output, workbook):
    from experiments.hierarchy100.evaluation import truth
    batch, output = Path(batch), Path(output)
    selected, calls = complete(batch); verify_seal(batch, output)
    sources, cells = truth(Path(workbook))
    dump(output/'truth.json', {'sources': list(sources.values()), 'independently_verified_cells': cells,
                              'workbook_sha256': sha(workbook), 'sheet': 'Table S1', 'order_column': 'G'})
    metrics, changes = [], []
    for row in selected:
        sample = row['sample']; expected = {sources[int(s)]['order'] for s in sample.split('x')}
        baseline = set(calls[sample]['methods']['z_only'])
        for method in METHODS:
            called = set(calls[sample]['methods'][method])
            metrics.append({'sample': sample, 'method': method, 'components': len(sample.split('x')),
                 'genus_development_ten': row['genus_development_ten'],
                 **counts(len(expected & called), len(called-expected), len(expected-called)),
                 'expected': ';'.join(sorted(expected)), 'predicted': ';'.join(sorted(called)),
                 'false_positives': ';'.join(sorted(called-expected)), 'false_negatives': ';'.join(sorted(expected-called))})
            changes.append({'sample': sample, 'method': method,
                 'rescued_z_false_negatives': len((called-baseline) & expected),
                 'new_false_positives': len(called-baseline-expected),
                 'lost_z_true_positives': len((baseline-called) & expected),
                 'removed_z_false_positives': len(baseline-called-expected),
                 'rescued_orders': ';'.join(sorted((called-baseline) & expected)),
                 'new_wrong_orders': ';'.join(sorted(called-baseline-expected))})
    summaries = []
    subsets = [('all100', lambda r: True), ('genus_development10', lambda r: r['genus_development_ten']),
               ('other90', lambda r: not r['genus_development_ten'])]
    subsets += [('components'+str(n), lambda r, n=n: r['components'] == n) for n in (3, 6, 10)]
    changes_by_key = {(r['sample'], r['method']): r for r in changes}
    for subset, predicate in subsets:
        for method in METHODS:
            sub = [r for r in metrics if r['method'] == method and predicate(r)]
            delta = [changes_by_key[(r['sample'], method)] for r in sub]
            summaries.append({'subset': subset, 'method': method, 'samples': len(sub),
                 **counts(*(sum(r[k] for r in sub) for k in ('TP', 'FP', 'FN'))),
                 'exact_sets': sum(not r['FP'] and not r['FN'] for r in sub),
                 **{k: sum(r[k] for r in delta) for k in ('rescued_z_false_negatives', 'new_false_positives',
                      'lost_z_true_positives', 'removed_z_false_positives')}})
    write_csv(output/'summary.csv', summaries); write_csv(output/'by_sample.csv', metrics)
    write_csv(output/'changes_from_z.csv', changes)
    dump(output/'EVALUATION_AUDIT.json', {'all100_sealed_before_evaluation': True,
         'literal_order_truth_openpyxl_XML_agree': True, 'all_raw_scores_reproduce_predictions': True,
         'baseline_matches_prior_seal': True, 'truth_uploaded': False,
         'independent_biological_validation': False, 'mixture_accuracy_used_to_select_threshold': False})
    print([r for r in summaries if r['subset'] == 'all100'])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('command', choices=('seal', 'evaluate'))
    p.add_argument('--batch', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--original', type=Path); p.add_argument('--workbook', type=Path)
    a = p.parse_args()
    if a.command == 'seal': seal(a.batch, a.output, a.original)
    else: evaluate(a.batch, a.output, a.workbook)
