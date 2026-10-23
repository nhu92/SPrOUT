import subprocess
import sys

import numpy as np

from experiments.order_ml.common import quorum, combine, choose_thresholds, reference_metrics
from experiments.order_ml.model import order_taxa, pseudo_features


def test_candidate_axis_is_order_even_when_families_differ():
    taxa = order_taxa([('A', 'fam1', 'sp1'), ('A', 'fam2', 'sp2'), ('B', 'fam3', 'sp3')])
    assert list(taxa[:, 1]) == ['A', 'A', 'B']
    assert set(taxa[:, 0]) == {'Angiosperms'}


def test_species_leaveout_and_entire_order_absence(monkeypatch):
    import experiments.order_ml.model as module
    taxa = np.array([['root', 'A', 'sp1'], ['root', 'A', 'sp1'], ['root', 'A', 'sp2'],
                     ['root', 'B', 'sp3'], ['root', 'C', 'sp4'], ['root', 'D', 'sp5']])
    d = np.ones((6, 6)); np.fill_diagonal(d, 0)
    captured = []
    def probe(q, reference, anchors):
        captured.append(anchors.tolist()); return [], []
    monkeypatch.setattr(module, 'profiles', probe)
    pseudo_features(d, taxa, 0, {'A', 'B', 'C'})
    assert [t[2] for t in captured[-1]] == ['sp2', 'sp3', 'sp4']
    pseudo_features(d, taxa, 0, {'A', 'B', 'C'}, True)
    assert [t[1] for t in captured[-1]] == ['B', 'C']


def test_exons_from_one_gene_cannot_create_a_quorum():
    hits = [{'gene': 'g1', 'order': 'A', 'score': .99} for _ in range(99)]
    assert not quorum(hits, .9)
    hits += [{'gene': g, 'order': 'A', 'score': .95} for g in ('g2', 'g3')]
    assert quorum(hits, .9) == {'A'}


def test_absent_species_false_quorum_is_not_diluted_by_other_loci():
    hits = [{'species': 's', 'context': 'absent', 'gene': str(g), 'order': 'WRONG', 'score': .95}
            for g in range(3)]
    hits += [{'species': 's', 'context': 'absent', 'gene': str(g), 'order': 'OTHER', 'score': .1}
             for g in range(3, 50)]
    result = reference_metrics(hits, {'s': 'A'}, .9)
    assert result['absent_species_rejection'] == 0


def test_two_reference_only_operating_points_are_fixed_by_quorum():
    hits = []
    for g in range(3):
        hits += [{'species': 's', 'context': 'known', 'gene': str(g), 'order': 'A', 'score': .9},
                 {'species': 's', 'context': 'absent', 'gene': str(g), 'order': 'B', 'score': .8}]
    thresholds = choose_thresholds(hits, {'s': 'A'})
    assert thresholds['balanced']['threshold'] == .9
    assert thresholds['specificity']['threshold'] == .9
    assert thresholds['specificity_target_met']


def test_union_preserves_baseline_even_when_ml_abstains():
    thresholds = {k: {'threshold': .99} for k in ('balanced', 'specificity')}
    calls = combine({'A', 'B'}, [], thresholds)
    assert calls['z_plus_ml_specificity'] == ['A', 'B']
    assert calls['z_plus_ml_balanced'] == ['A', 'B']
    assert calls['ml_specificity'] == []


def test_failure_to_calibrate_specificity_forces_abstention_even_at_score_one():
    thresholds = {k: {'threshold': 1.} for k in ('balanced', 'specificity')}
    thresholds['specificity_target_met'] = False
    hits = [{'gene': str(g), 'order': 'WRONG', 'score': 1.} for g in range(3)]
    assert combine({'A'}, hits, thresholds)['z_plus_ml_specificity'] == ['A']


def test_standard_library_import_for_local_postprocessing():
    r = subprocess.run([sys.executable, '-S', '-c', 'import experiments.order_ml.evaluation'],
                       text=True, capture_output=True)
    assert r.returncode == 0, r.stderr


def test_catalog_coverage_duplicate_accessions_and_whole_order_split(tmp_path):
    from experiments.order_ml.model import eligible_catalog
    for gene in range(50):
        records = [(f'O{i:02d}_F_G{i:02d}{j}_alpha_acc', 'A'*80)
                   for i in range(32) for j in range(2)]
        records.append(('O00_F_G000_alpha_better', 'A'*90))
        if gene < 29:
            records.append(('O00_F_Rare_alpha_acc', 'A'*80))
        (tmp_path/(str(gene)+'.fasta')).write_text(''.join('>'+name+'\n'+seq+'\n' for name, seq in records))
    catalog = eligible_catalog(tmp_path)
    assert catalog['coverage']['g000_alpha'] == 50
    assert 'rare_alpha' not in catalog['eligible_species']
    assert catalog['genes']['0']['g000_alpha']['name'].endswith('_better')
    split = catalog['split']
    assert len(split['fit']) == 12
    assert not set(split['fit']) & set(split['calibration']+split['test'])
    assert not set(split['calibration']) & set(split['test'])
