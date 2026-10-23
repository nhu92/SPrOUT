"""Reference-only, whole-order-held-out training; no mixture truth inputs."""
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil
import subprocess
import time

from Bio import SeqIO
import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from experiments.reference_atlas.transfer import taxonomy, tree_distances, profiles, FEATURES
from .common import dump, read, sha, stable, choose_thresholds, reference_metrics


def order_taxa(taxa):
    """The generic candidate axis is column 1; all orders share one root."""
    return np.array([('Angiosperms', tax[0], tax[2]) for tax in taxa])


def eligible_catalog(reference):
    genes = sorted(Path(reference).glob('*.fasta'))
    if len(genes) != 50:
        raise ValueError('Expected fifty original order alignments')
    catalog, coverage, order_by_species = {}, Counter(), {}
    for path in genes:
        best = {}
        for rec in SeqIO.parse(path, 'fasta'):
            tax = taxonomy(rec.id)
            if not tax or tax[2].split('_', 1)[1] in ('sp', 'spp', 'cf', 'aff', 'indet'):
                continue
            bases = sum(c in 'ACGT' for c in str(rec.seq).upper())
            if bases < 80:
                continue
            species = tax[2]
            if species in order_by_species and order_by_species[species] != tax[0]:
                raise ValueError('Species assigned to multiple orders: '+species)
            order_by_species[species] = tax[0]
            key = (-bases, rec.id)
            if species not in best or key < best[species]['choice']:
                best[species] = {'choice': key, 'name': rec.id, 'order': tax[0], 'bases': bases}
        catalog[path.stem] = best
        coverage.update(best.keys())
    eligible = {s for s, n in coverage.items() if n >= 30}
    by_order = defaultdict(list)
    for s in sorted(eligible):
        by_order[order_by_species[s]].append(s)
    multi = sorted((o for o, sp in by_order.items() if len(sp) >= 2),
                   key=lambda o: stable('order-ml-split/'+o))
    if len(multi) < 30:
        raise ValueError(f'Need at least 30 eligible multi-species orders; found {len(multi)}')
    split = {'calibration': multi[:10], 'test': multi[10:20], 'fit': multi[20:],
             'anchor_only': sorted(o for o, sp in by_order.items() if len(sp) < 2)}
    clean = {g: {s: {k: v for k, v in r.items() if k != 'choice'}
                 for s, r in rows.items() if s in eligible} for g, rows in catalog.items()}
    return {'genes': clean, 'coverage': dict(coverage), 'eligible_species': sorted(eligible),
            'order_species': dict(by_order), 'split': split}


def build_trees(reference, output, workers=8):
    output = Path(output); output.mkdir()
    executable = shutil.which('FastTree') or shutil.which('fasttree')
    if not executable:
        raise FileNotFoundError('FastTree')
    def one(path):
        trimmed, tree = output/(path.stem+'.fasta'), output/(path.stem+'.tre')
        subprocess.run(['trimal', '-in', str(path), '-out', str(trimmed), '-gt', '.5'],
                       check=True, capture_output=True, timeout=1800)
        with tree.open('w') as handle:
            r = subprocess.run([executable, '-gtr', '-gamma', '-nt', str(trimmed)],
                               stdout=handle, stderr=subprocess.PIPE, text=True, timeout=1800)
        (output/(path.stem+'.log')).write_text(r.stderr)
        if r.returncode:
            raise RuntimeError(r.stderr[-2000:])
        tree_distances(tree)
        return {'gene': path.stem, 'tree_sha256': sha(tree), 'alignment_sha256': sha(trimmed)}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(one, sorted(Path(reference).glob('*.fasta'))))


def block(tree, catalog, project=None):
    names, distances = tree_distances(tree)
    normalized = {name.strip("'").removeprefix('_R_'): i for i, name in enumerate(names)}
    if len(normalized) != len(names):
        raise ValueError('Duplicate normalized leaf IDs')
    refs, taxa = [], []
    for species, rec in sorted(catalog.items()):
        if rec['name'] not in normalized:
            raise ValueError('Frozen eligible reference missing: '+rec['name'])
        refs.append(normalized[rec['name']])
        taxa.append((rec['order'], '', species))
    refs = np.array(refs, dtype=int)
    queries = [i for i, n in enumerate(names) if project and project in n]
    if project and not queries:
        raise ValueError('No mixture query leaves: '+str(tree))
    return (distances[np.ix_(refs, refs)], order_taxa(taxa),
            distances[np.ix_(queries, refs)], [names[i] for i in queries])


def pseudo_features(d, taxa, i, allowed, absent=False):
    keep = np.array([o in allowed for o in taxa[:, 1]])
    keep &= (taxa[:, 1] != taxa[i, 1]) if absent else (taxa[:, 2] != taxa[i, 2])
    if not absent and not (taxa[keep, 1] == taxa[i, 1]).any():
        return None
    if len(set(taxa[keep, 1])) < 2:
        return None
    return profiles(d[i, keep], d[np.ix_(keep, keep)], taxa[keep])


def heldout(blocks, catalog, partition, model):
    split = catalog['split']
    allowed = set(split['fit']+split['anchor_only']+split[partition])
    eligible = {s: order for order in split[partition] for s in catalog['order_species'][order]}
    rows = []
    represented = defaultdict(set)
    for gene, (d, taxa) in blocks.items():
        for i, tax in enumerate(taxa):
            species = str(tax[2])
            if species not in eligible:
                continue
            for absent in (False, True):
                value = pseudo_features(d, taxa, i, allowed, absent)
                if value is None:
                    continue
                x, labels = value
                winner = int(x[:, 0].argmin())
                probability = float(model.predict_proba(x[winner:winner+1])[0, 1])
                rows.append({'gene': gene, 'species': species, 'truth_order': str(tax[1]),
                             'context': 'absent' if absent else 'known',
                             'order': labels[winner][0], 'score': probability})
                if not absent:
                    represented[species].add(gene)
    # Require every held-out species to have a measurable multi-gene known control.
    insufficient = {sp: len(represented[sp]) for sp in eligible if len(represented[sp]) < 3}
    if insufficient:
        raise ValueError('Insufficient held-out controls: '+str(insufficient))
    return rows, eligible


def train(reference, output):
    output = Path(output); output.mkdir()
    start = time.monotonic()
    catalog = eligible_catalog(reference)
    dump(output/'CATALOG.json', catalog)
    trees = build_trees(reference, output/'trees')
    tree_end = time.monotonic()
    blocks = {gene: block(output/'trees'/(gene+'.tre'), refs)[:2]
              for gene, refs in catalog['genes'].items()}
    split = catalog['split']; allowed = set(split['fit']+split['anchor_only'])
    rows = []
    for gene, (d, taxa) in blocks.items():
        ids = sorted((i for i, t in enumerate(taxa) if t[1] in split['fit']),
                     key=lambda i: stable(gene+'/'+taxa[i, 2]))[:64]
        for i in ids:
            species, true_order = str(taxa[i, 2]), str(taxa[i, 1])
            value = pseudo_features(d, taxa, i, allowed)
            if value is None:
                continue
            x, labels = value
            ranked = sorted(range(len(labels)), key=lambda j: x[j, 0])
            selected = [j for j in ranked if labels[j][0] == true_order]
            selected += [j for j in ranked if labels[j][0] != true_order][:3]
            for j in selected:
                rows.append({'x': x[j].tolist(), 'y': int(labels[j][0] == true_order),
                             'species': species, 'gene': gene, 'absent': False})
            value = pseudo_features(d, taxa, i, allowed, True)
            if value:
                nx, _ = value; j = int(nx[:, 0].argmin())
                rows.append({'x': nx[j].tolist(), 'y': 0, 'species': species,
                             'gene': gene, 'absent': True})
    if len(rows) < 100:
        raise ValueError('Insufficient fitting rows')
    count = Counter(r['species'] for r in rows)
    weights = np.array([1/count[r['species']] for r in rows]); weights *= len(weights)/weights.sum()
    model = make_pipeline(StandardScaler(), LogisticRegression(C=.1, class_weight='balanced', max_iter=2000))
    with threadpool_limits(limits=1):
        model.fit([r['x'] for r in rows], [r['y'] for r in rows],
                  logisticregression__sample_weight=weights)
        cal_rows, cal_species = heldout(blocks, catalog, 'calibration', model)
        thresholds = choose_thresholds(cal_rows, cal_species)
        test_rows, test_species = heldout(blocks, catalog, 'test', model)
    if int(model[-1].n_iter_.max()) >= 2000:
        raise ValueError('Logistic regression failed to converge')
    bundle = {'schema': 'order-full-profile-v1', 'model': model, 'features': FEATURES,
              'thresholds': thresholds, 'catalog': catalog, 'minimum_genes': 3}
    joblib.dump(bundle, output/'model.joblib')
    dump(output/'CALIBRATION_ROWS.json', cal_rows)
    dump(output/'REFERENCE_TEST_ROWS.json', test_rows)
    dump(output/'THRESHOLDS.json', thresholds)
    audit = {'schema': bundle['schema'], 'model_sha256': sha(output/'model.joblib'),
             'features': FEATURES, 'split': split, 'fit_rows': len(rows), 'fit_species': len(count),
             'fit_order_absent_negatives': sum(r['absent'] for r in rows),
             'fit_labels': dict(Counter(r['y'] for r in rows)),
             'calibration_species': len(cal_species), 'test_species': len(test_species),
             'thresholds': {k: thresholds[k] for k in ('balanced', 'specificity', 'specificity_target_met')},
             'reference_test': {k: reference_metrics(test_rows, test_species, thresholds[k]['threshold'])
                                for k in ('balanced', 'specificity')},
             'trees': trees, 'tree_seconds': tree_end-start, 'total_seconds': time.monotonic()-start,
             'no_mixture_queries_or_truth_used_for_fitting_or_thresholds': True,
             'shared_reference_only_phylogenetic_scaffold': True}
    dump(output/'TRAINING.json', audit)
    return audit
