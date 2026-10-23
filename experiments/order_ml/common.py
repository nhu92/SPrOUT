"""Order verifier protocol primitives; safe to import without scientific packages."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path

METHODS = ('z_only', 'ml_balanced', 'z_plus_ml_balanced',
           'ml_specificity', 'z_plus_ml_specificity')


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def dump(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def stable(value):
    return hashlib.sha256(value.encode()).hexdigest()


def quorum(rows, threshold, minimum_genes=3):
    support = defaultdict(set)
    for row in rows:
        if row['score'] >= threshold:
            support[row['order']].add(row['gene'])
    return {order for order, genes in support.items() if len(genes) >= minimum_genes}


def reference_metrics(rows, species_orders, threshold):
    groups = defaultdict(list)
    for row in rows:
        groups[(row['species'], row['context'])].append(row)
    n = len(species_orders)
    if not n:
        raise ValueError('No held-out reference species')
    known = sum(order in quorum(groups[(sp, 'known')], threshold)
                for sp, order in species_orders.items()) / n
    rejected = sum(not quorum(groups[(sp, 'absent')], threshold)
                   for sp in species_orders) / n
    # Also disclose incorrect orders alongside correct calls in known context.
    known_fp = sum(bool(quorum(groups[(sp, 'known')], threshold)-{order})
                   for sp, order in species_orders.items()) / n
    return {'threshold': threshold, 'species': n, 'known_correct_acceptance': known,
            'absent_species_rejection': rejected, 'known_species_any_wrong_order': known_fp,
            'balanced_reference_utility': (known+rejected)/2}


def choose_thresholds(rows, species_orders):
    curve = [reference_metrics(rows, species_orders, i/100) for i in range(101)]
    balanced = max(curve, key=lambda r: (r['balanced_reference_utility'], r['threshold']))
    candidates = [r for r in curve if r['absent_species_rejection'] >= .99]
    specificity = max(candidates, key=lambda r: (r['known_correct_acceptance'],
                       r['absent_species_rejection'], r['threshold'])) if candidates else curve[-1]
    return {'balanced': balanced, 'specificity': specificity,
            'specificity_target_met': bool(candidates), 'curve': curve}


def combine(z, hits, thresholds):
    result = {'z_only': sorted(z)}
    for label in ('balanced', 'specificity'):
        failed_specificity = label == 'specificity' and not thresholds.get('specificity_target_met', True)
        calls = set() if failed_specificity else quorum(hits, thresholds[label]['threshold'])
        result['ml_'+label] = sorted(calls)
        result['z_plus_ml_'+label] = sorted(set(z) | calls)
    return {key: result[key] for key in METHODS}
