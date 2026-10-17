"""Full patristic profiles, gene-balanced routing and genus C-style fusion.

The existing statistics engine's internal `families` axis is adapted to genera
and its parent axis to families. Tree labels and distance values are not edited.
The bounded genus contrast replaces reference-count-dependent ACS in this new
experiment; these are C-style genus rules, not the unchanged family predictor.
"""
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np

from experiments.reference_atlas.transfer import tree_distances, profiles


def stable(text):
    return int(hashlib.sha256(text.encode()).hexdigest()[:12], 16)


def blocks(path):
    names, distances = tree_distances(path)
    refs, queries, taxa, reference_names = [], [], [], []
    seen = set()
    for i, name in enumerate(names):
        if name.startswith('QUERY_'):
            queries.append(i)
            continue
        fields = name.removeprefix('_R_').split('_', 4)
        if len(fields) != 5:
            raise ValueError('Unrecognized tree leaf: ' + name)
        identity = (fields[2]+'_'+fields[3]).casefold()
        if identity in seen:
            continue
        seen.add(identity)
        refs.append(i)
        taxa.append((fields[1], fields[2], identity))
        reference_names.append(name.removeprefix('_R_'))
    return (distances[np.ix_(refs, refs)], np.array(taxa),
            distances[np.ix_(queries, refs)], [names[i] for i in queries], reference_names)


def route(trees, anchor_cells, cells, index, min_genes=2, ambiguous_genes=3, gap=.10):
    primary, secondary, anchor_votes = defaultdict(set), defaultdict(set), defaultdict(set)
    valid = set()
    cell_keys = {c.key for c in cells}
    for path in trees:
        gene = path.name.split('_')[0]
        _, _, queries, _, names = blocks(path)
        if not len(queries):
            continue
        valid.add(gene)
        groups = {key: [i for i, name in enumerate(names) if key in anchor_cells.get(name, [])] for key in cell_keys}
        groups = {key: ids for key, ids in groups.items() if ids}
        for query in queries:
            nearest = {key:float(query[ids].min()) for key, ids in groups.items()}
            if not nearest:
                continue
            ranked = sorted(nearest, key=lambda key:(nearest[key],key))
            best = nearest[ranked[0]]
            winners = [key for key in ranked if np.isclose(nearest[key],best,rtol=1e-9,atol=1e-12)]
            for key in winners:
                primary[key].add(gene)
            for i, name in enumerate(names):
                if np.isclose(query[i], best, rtol=1e-9, atol=1e-12):
                    anchor_votes[name].add(gene)
            if len(ranked)>1 and (nearest[ranked[1]]-best)/max(nearest[ranked[1]],1e-12)<=gap:
                secondary[ranked[1]].add(gene)
    if len(valid)<10:
        raise ValueError('Routing requires at least ten different genes')
    represented = set().union(*anchor_cells.values()) if anchor_cells else set()
    selected = [c for c in cells if len(primary[c.key])>=min_genes
                or len(secondary[c.key])>=ambiguous_genes or c.key not in represented]
    rows = [{'cell':c.key,'family':c.family,'genera':len(index.genera(c)),
             'nearest_genes':len(primary[c.key]),'ambiguous_second_genes':len(secondary[c.key]),
             'quality_unresolved':c.key not in represented,'selected':c in selected} for c in cells]
    return selected, rows, anchor_votes


def final_scores(trees, model_path, project):
    from .gene_statistics import analyze, Settings
    bundle = joblib.load(model_path)
    model, threshold = bundle['model'], bundle['threshold']
    per_exon = []
    parent = {}
    genes = set()
    ml_votes, nearest_votes = defaultdict(set), defaultdict(set)
    ml_rows = []
    for path in trees:
        gene = path.name.split('_')[0]
        reference, taxa, queries, query_names, _ = blocks(path)
        if len(taxa)<4 or not len(queries):
            continue
        genes.add(gene)
        genera = sorted(set(taxa[:,1]))
        groups = {genus:np.flatnonzero(taxa[:,1]==genus) for genus in genera}
        parent.update({genus:str(taxa[ids[0],0]) for genus,ids in groups.items()})
        contrast, margins = defaultdict(float), defaultdict(float)
        present = defaultdict(set)
        for genus in genera:
            present[parent[genus]].add(genus)
        comparable = {genus for genus in genera if len(present[parent[genus]])>=2
                      or len([f for f in genera if parent[f]!=parent[genus]])>0}
        winners_seen = set()
        for query, query_name in zip(queries,query_names):
            near = {genus:float(query[ids].min()) for genus,ids in groups.items()}
            means = {genus:float(query[ids].mean()) for genus,ids in groups.items()}
            background = float(np.median(list(means.values())))
            for genus in genera:
                contrast[genus] = max(contrast[genus],max(0.,min(1.,(background-means[genus])/background)) if background>0 else 0.)
            best = min(near.values())
            winners = [genus for genus in genera if np.isclose(near[genus],best,rtol=1e-9,atol=1e-12)]
            for genus in winners:
                nearest_votes[genus].add(gene)
                winners_seen.add(genus)
            if len(winners)==1:
                genus = winners[0]
                peers = [near[g] for g in genera if g!=genus and parent[g]==parent[genus]]
                if not peers:
                    peers = [near[g] for g in genera if g!=genus]
                if peers and min(peers)>0:
                    margins[genus] = max(margins[genus],max(0.,min(1.,(min(peers)-best)/min(peers))))
            x, labels = profiles(query,reference,taxa)
            if len(x):
                winner = int(x[:,0].argmin())
                genus = labels[winner][0]
                probability = float(model.predict_proba(x[winner:winner+1])[0,1])
                accepted = probability>=threshold
                if accepted:
                    ml_votes[genus].add(gene)
                ml_rows.append({'gene':gene,'exon':path.name,'query':query_name,'genus':genus,
                                'score':probability,'accepted':accepted})
        per_exon.append((gene,genera,dict(contrast),dict(margins),comparable))
    if len(genes)<10:
        raise ValueError('Final scoring requires at least ten different genes')
    genes, genera = sorted(genes), sorted(parent)
    shape = (len(genes),len(genera))
    a,w,e,k,v = [np.zeros(shape) for _ in range(5)]
    for i,gene in enumerate(genes):
        exons = [row for row in per_exon if row[0]==gene]
        for j,genus in enumerate(genera):
            available = [row for row in exons if genus in row[1]]
            comparable = [row for row in available if genus in row[4]]
            a[i,j] = sum(row[2].get(genus,0.) for row in available)/len(exons)
            e[i,j] = bool(comparable)
            w[i,j] = sum(row[3].get(genus,0.) if row[3].get(genus,0.)>=.1 else 0. for row in comparable)/len(comparable) if comparable else 0.
            k[i,j] = any(row[3].get(genus,0.)>=.1 for row in comparable)
            v[i,j] = gene in nearest_votes[genus]
    data = {'project':project,'genes':genes,'families':genera,'family_parent':parent,
            'acs':a.tolist(),'weights':w.tolist(),'eligible':e.tolist(),'informative':k.tolist(),'votes':v.tolist(),
            'reference_species':{},'missing_sibling_genes':{},'input_sha256':{}}
    statistics = analyze(data,Settings())
    calls = {name:[] for name in ('balanced_z','c_style','core','auxiliary','merged','standalone_ml','review')}
    evidence = []
    for row in statistics['rows']:
        genus = row['family']
        votes = len(ml_votes[genus])
        c = row['selected']
        core = c and votes>=2
        auxiliary = c and row['rescue_selected'] and row['selection_frequency_mc_interval'][0]>.5 and votes>=1
        merged = core or auxiliary
        review = not merged and (row['hybrid_selected'] or votes>=3)
        flags = {'balanced_z':row['native_selected'],'c_style':c,'core':core,'auxiliary':auxiliary,
                 'merged':merged,'standalone_ml':votes>=3,'review':review}
        for method,selected in flags.items():
            if selected:
                calls[method].append(genus)
        evidence.append({'genus':genus,'family':parent[genus],'ml_gene_votes':votes,
                         'nearest_gene_votes':len(nearest_votes[genus]),'z':row['z'],
                         'discriminating_genes':row['discriminating_genes'],'comparable_genes':row['comparable_genes'],
                         'stability_lower':row['selection_frequency_mc_interval'][0],**flags})
    statistics['rank']='genus';statistics['parent_rank']='family'
    statistics['z_origin']='Genus-balanced contrast standardized in the final panel; Z is not a calibrated genus presence probability'
    statistics['rows']=[{**{k:v for k,v in r.items() if k not in ('family','parent_order')},
                         'genus':r['family'],'family':r['parent_order']} for r in statistics['rows']]
    return {name:sorted(values) for name,values in calls.items()}, evidence, ml_rows, statistics, data
