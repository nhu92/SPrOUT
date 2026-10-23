"""Freeze, batch-fit, score and collect the fixed order ML experiment on TTU."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import importlib.metadata
from pathlib import Path
import shutil
import sys
import tarfile
import time
import traceback

from .common import read, dump, sha, combine


def hashes(paths):
    paths = list(paths)
    with ThreadPoolExecutor(max_workers=8) as pool:
        return dict(zip(map(str, paths), pool.map(sha, paths)))


def verify(mapping):
    actual = hashes(map(Path, mapping))
    bad = [p for p in mapping if actual[p] != mapping[p]]
    if bad:
        raise ValueError('Immutable input changed: '+str(bad[:5]))


def utc():
    return datetime.now(timezone.utc).isoformat()


def freeze(root):
    request = read(root/'REQUEST.json'); source = Path(request['source'])
    old = read(source/'selection.json')['selected']
    if len(old) != 100 or len({r['sample'] for r in old}) != 100:
        raise ValueError('Need original hundred unique samples')
    if not read(source/'AUDIT.json')['all_100_complete']:
        raise ValueError('Original source audit incomplete')
    selected = []
    for row in old:
        stage = Path(row['run'])/'01_order'
        valid = read(stage/'family_tree_validation.json')['valid_trees']
        if not valid or len(valid) != len(set(valid)):
            raise ValueError('Missing/duplicated valid order trees')
        paths = [stage/'valid_phylo'/name for name in valid]
        tree_hash = hashes(paths)
        prior_path = Path(row['run'])/'ORDER_FAMILY_PRODUCTS.json'
        prior = read(prior_path)
        if any(prior.get(p) != digest for p, digest in tree_hash.items()):
            raise ValueError('Order tree differs from original frozen run')
        zpath = stage/'native_z_candidates.txt'; orderpath = stage/'ORDER.json'
        z = set(zpath.read_text().splitlines())
        expected = {r['taxon'] for r in read(orderpath)['evidence']['predictions'] if float(r['z_score']) > .9}
        if z != expected:
            raise ValueError('Strict Z baseline mismatch')
        selected.append({k: row[k] for k in ('index', 'sample', 'project', 'genus_development_ten')})
        selected[-1].update(trees=tree_hash, baseline_path=str(zpath), baseline_sha256=sha(zpath),
                            order_json_path=str(orderpath), order_json_sha256=sha(orderpath), z_only=sorted(z),
                            prior_products_path=str(prior_path), prior_products_sha256=sha(prior_path))
    dump(root/'selection.json', {'selected': selected})
    code = list((root/'code').rglob('*.py'))+[root/'PLAN.md']
    alignments = list(Path(request['reference']).glob('*.fasta'))
    if len(alignments) != 50:
        raise ValueError('Expected fifty reference alignments')
    dump(root/'PROTOCOL.json', {**request, 'schema': 'order-ml-100-v1', 'frozen_utc': utc(),
         'source_selection_sha256': sha(source/'selection.json'),
         'source_audit_sha256': sha(source/'AUDIT.json'), 'selection_sha256': sha(root/'selection.json'),
         'code_sha256': hashes(code), 'reference_sha256': hashes(alignments),
         'runtime': {'python': sys.version, **{p: importlib.metadata.version(p)
                      for p in ('numpy', 'scipy', 'scikit-learn', 'biopython', 'joblib')}},
         'truth_uploaded': False, 'mixture_truth_consumed_for_training': False})


def verify_root(root, reference=False):
    protocol = read(root/'PROTOCOL.json')
    if sha(root/'selection.json') != protocol['selection_sha256']:
        raise ValueError('Changed cohort/tree/baseline manifest')
    verify(protocol['code_sha256'])
    if reference:
        verify(protocol['reference_sha256'])
    return protocol


def prepare(root):
    from .model import train
    freeze(root)
    protocol = verify_root(root, True)
    result = train(protocol['reference'], root/'model')
    verify_root(root, True)
    files = [p for p in (root/'model').rglob('*') if p.is_file()]
    dump(root/'MODEL_FROZEN.json', {'frozen_utc': utc(), 'sha256': hashes(files),
         'protocol_sha256': sha(root/'PROTOCOL.json'), 'model_sha256': result['model_sha256'],
         'thresholds': result['thresholds'], 'all_fitting_complete_before_any_mixture_inference': True})
    print('Frozen order-specific checkpoint and reference-only thresholds', flush=True)


def run(root, index):
    import joblib
    from .model import block, profiles
    start = time.monotonic()
    verify_root(root)
    frozen = read(root/'MODEL_FROZEN.json')
    verify(frozen['sha256'])
    if sha(root/'PROTOCOL.json') != frozen['protocol_sha256']:
        raise ValueError('Changed protocol')
    row = read(root/'selection.json')['selected'][index]
    if row['index'] != index:
        raise ValueError('Array index mismatch')
    directory = root/'samples'/row['project']; directory.mkdir(parents=True, exist_ok=True)
    if (directory/'SUCCESS').exists() or (directory/'FAILED').exists():
        raise FileExistsError('Preserve previous result; do not overwrite')
    try:
        verify(row['trees'])
        verify({row['baseline_path']: row['baseline_sha256'], row['order_json_path']: row['order_json_sha256']})
        bundle = joblib.load(root/'model/model.joblib')
        def one(tree):
            path = Path(tree); gene = path.name.split('_')[0]
            d, taxa, queries, names = block(path, bundle['catalog']['genes'][gene], row['project'])
            hits = []
            for q, name in zip(queries, names):
                x, labels = profiles(q, d, taxa)
                winner = int(x[:, 0].argmin())
                score = float(bundle['model'].predict_proba(x[winner:winner+1])[0, 1])
                hits.append({'gene': gene, 'tree': path.name, 'query': name,
                             'order': labels[winner][0], 'score': score, 'features': x[winner].tolist()})
            return hits
        hits = []
        with ThreadPoolExecutor(max_workers=8) as pool:
            for rows in pool.map(one, sorted(row['trees'])):
                hits.extend(rows)
        calls = {'sample': row['sample'], 'project': row['project'],
                 'methods': combine(row['z_only'], hits, bundle['thresholds'])}
        dump(directory/'HITS.json', hits); dump(directory/'predictions.json', calls)
        for method, values in calls['methods'].items():
            (directory/(method+'.txt')).write_text(''.join(v+'\n' for v in values))
        verify(row['trees'])
        verify({row['baseline_path']: row['baseline_sha256'], row['order_json_path']: row['order_json_sha256']})
        verify_root(root)
        dump(directory/'AUDIT.json', {'sample': row['sample'], 'index': index, 'trees': len(row['trees']),
             'genes': len({Path(p).name.split('_')[0] for p in row['trees']}), 'query_hits': len(hits),
             'model_sha256': frozen['model_sha256'], 'model_freeze_sha256': sha(root/'MODEL_FROZEN.json'),
             'input_trees_verified_before_and_after': True, 'seconds': time.monotonic()-start})
        outputs = [p for p in directory.iterdir() if p.is_file()]
        dump(directory/'SUCCESS', {'sha256': hashes(outputs), 'finished_utc': utc()})
    except Exception:
        (directory/'FAILED').write_text(traceback.format_exc())
        raise


def collect(root):
    verify_root(root, True); frozen = read(root/'MODEL_FROZEN.json'); verify(frozen['sha256'])
    selected = read(root/'selection.json')['selected']; calls = []; timings = []
    all_inputs = {}; output_hashes = {}
    for row in selected:
        directory = root/'samples'/row['project']
        if (directory/'FAILED').exists() or not (directory/'SUCCESS').exists():
            raise ValueError('Incomplete sample '+row['sample'])
        record = read(directory/'SUCCESS'); verify(record['sha256']); output_hashes.update(record['sha256'])
        all_inputs.update(row['trees']); all_inputs[row['baseline_path']] = row['baseline_sha256']
        all_inputs[row['order_json_path']] = row['order_json_sha256']
        all_inputs[row['prior_products_path']] = row['prior_products_sha256']
        calls.append(read(directory/'predictions.json')); timings.append(read(directory/'AUDIT.json'))
    verify(all_inputs)
    dump(root/'predictions.json', calls); dump(root/'TIMINGS.json', timings)
    dump(root/'AUDIT.json', {'all_100_complete': len(calls) == 100, 'code_model_inputs_outputs_verified': True,
         'samples': len(calls), 'reused_order_trees': sum(len(r['trees']) for r in selected),
         'output_sha256': output_hashes, 'collected_utc': utc()})
    (root/'PREDICTIONS_COMPLETE').write_text(utc()+'\n')
    # Include checkpoint/reference-only trees and raw query scores, but no truth workbook.
    with tarfile.open(root/'results.tar.gz', 'w:gz') as archive:
        for name in ('PLAN.md', 'REQUEST.json', 'PROTOCOL.json', 'selection.json', 'MODEL_FROZEN.json',
                     'predictions.json', 'TIMINGS.json', 'AUDIT.json', 'PREDICTIONS_COMPLETE', 'model', 'samples', 'logs'):
            archive.add(root/name, arcname=name)
    print('All 100 sealed; archive complete', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=('prepare', 'run', 'collect'))
    p.add_argument('--root', type=Path, required=True); p.add_argument('--index', type=int)
    args = p.parse_args()
    if args.command == 'run': run(args.root, args.index)
    else: globals()[args.command](args.root)
