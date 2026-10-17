from collections import Counter
import json

import joblib
import numpy as np
import pytest
from dataclasses import dataclass,asdict

from experiments.adaptive_genus.phylogeny import TreeIndex,Cell
from experiments.adaptive_genus.scoring import blocks,route,final_scores
from experiments.adaptive_genus.select_samples import select
from experiments.adaptive_genus.references import ReferenceBank


def tree_index(tmp_path):
    path=tmp_path/'species.tree'
    path.write_text('((O_F_A_a_r1:1,O_F_C_c_r2:1):1,(O_F_B_b_r3:1,O_F_D_d_r4:1):1);')
    return TreeIndex(path)


def test_partition_preserves_every_placement_and_monotonic_guard(tmp_path):
    index=tree_index(tmp_path)
    initial=Cell('F',tuple(index.taxonomy))
    cells=index.partition([initial],budget=2,target_genera=1)
    assert len(cells)==2
    assert set().union(*(set(c.tips) for c in cells))==set(initial.tips)
    singleton=Cell('F',(cells[0].tips[0],))
    guarded=index.neighborhood_guard(cells,[singleton],singleton.tips,neighbors=1)
    assert len({g for c in guarded for g in index.genera(c)})==2
    assert set().union(*(set(c.tips) for c in guarded))<=set(initial.tips)
    assert index.neighborhood_guard(cells,[singleton],singleton.tips,neighbors=0)==[singleton]


def test_polyphyletic_genus_survives_any_remaining_placement(tmp_path):
    path=tmp_path/'poly.tree'
    path.write_text('((O_F_G_a_1,O_F_H_a_2),(O_F_G_b_3,O_F_J_a_4));')
    index=TreeIndex(path)
    cells=index.partition([Cell('F',tuple(index.taxonomy))],2,1)
    assert all('G' in index.genera(c) for c in cells)
    assert len(index.by_genus['G'])==2


def test_quality_center_is_interior_of_a_long_branch(tmp_path):
    path=tmp_path/'medoid.tree'
    path.write_text('((O_F_A_a_1:1,O_F_B_b_2:1):1,(O_F_C_c_3:1,O_F_Z_z_4:15):1);')
    index=TreeIndex(path);cell=Cell('F',tuple(index.taxonomy))
    ranked=index.medoid_order(cell)
    assert index.taxonomy[ranked[0]][2]!='Z'
    assert set(ranked)==set(cell.tips)


def gene_trees(tmp_path,count=10,same_gene=False):
    paths=[]
    for i in range(count):
        gene='100' if same_gene else str(100+i)
        path=tmp_path/f'{gene}_exon_{i}.tre'
        path.write_text('(((O_F_A_a_r1:.05,O_F_A_b_r2:.05):.01,QUERY_q:.01):.1,'
                        '(O_F_B_c_r3:.05,O_F_B_d_r4:.05):.2);')
        paths.append(path)
    return paths


def test_route_counts_genes_not_exons_and_does_not_call_anchors(tmp_path):
    index=tree_index(tmp_path)
    cells=index.partition([Cell('F',tuple(index.taxonomy))],2,1)
    mapping={'O_F_A_a_r1':[cells[0].key],'O_F_A_b_r2':[cells[0].key],
             'O_F_B_c_r3':[cells[1].key],'O_F_B_d_r4':[cells[1].key]}
    trees=gene_trees(tmp_path)
    chosen,rows,_=route(trees,mapping,cells,index)
    assert len(chosen)==1
    assert rows[0]['nearest_genes'] in (0,10)
    assert all('genus_prediction' not in r for r in rows)
    with pytest.raises(ValueError,match='ten different genes'):
        route(gene_trees(tmp_path,10,same_gene=True),mapping,cells,index)


def test_reference_failure_retains_unresolved_region(tmp_path):
    index=tree_index(tmp_path)
    cells=index.partition([Cell('F',tuple(index.taxonomy))],2,1)
    mapping={'O_F_A_a_r1':[cells[0].key],'O_F_A_b_r2':[cells[0].key]}
    selected,rows,_=route(gene_trees(tmp_path),mapping,cells,index)
    assert cells[1] in selected
    assert next(r for r in rows if r['cell']==cells[1].key)['quality_unresolved']


class AcceptVerifier:
    def predict_proba(self,x):
        return np.tile([.01,.99],(len(x),1))


def test_genus_rank_adapter_preserves_gene_budget_and_primary_set(tmp_path):
    model=tmp_path/'model.joblib'
    joblib.dump({'model':AcceptVerifier(),'threshold':.9},model)
    calls,evidence,_,statistics,features=final_scores(gene_trees(tmp_path),model,'synthetic')
    assert calls['merged']==['A']
    assert calls['merged']==calls['core']
    assert set(calls['merged'])<=set(calls['c_style'])
    assert {r['genus'] for r in statistics['rows']}=={'A','B'}
    assert np.max(features['acs'])<=1
    assert next(r for r in evidence if r['genus']=='A')['ml_gene_votes']==10


def test_cohort_is_fixed_and_stratified_without_taxon_labels():
    rows=[]
    for n in (3,6,10):
        for i in range(12):
            rows.append({'sample':'x'.join(str(i*20+j) for j in range(n))})
    chosen=select(rows)
    assert Counter(len(r['sample'].split('x')) for r in chosen)=={3:4,6:3,10:3}
    assert select(list(reversed(rows)))==chosen


def test_panel_metadata_keeps_full_species_name_and_epithet_separate():
    @dataclass
    class Taxon:
        order:str='O'
        family:str='F'
        genus:str='G'
        species:str='epithet'
        sprout_prefix:str='O_F_G_epithet'
    @dataclass
    class Recovery:
        species:str='G epithet'
        sequence_id:str='X'
        def to_dict(self):
            return asdict(self)
    bank=ReferenceBank.__new__(ReferenceBank)
    bank.entries={'record':{'record':Recovery(),'taxonomy':Taxon()}}
    bank.quality={'record':{'accepted_loci':40}}
    row=bank.metadata('record')
    assert row['species']=='G epithet' and row['species_epithet']=='epithet'
