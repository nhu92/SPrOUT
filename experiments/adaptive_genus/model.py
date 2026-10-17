"""Genus-held-out reference-only full-profile match verifier."""
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil
import time

import joblib
import numpy as np
from Bio import SeqIO
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from experiments.reference_atlas.transfer import FEATURES, profiles
from .scoring import blocks, stable
from .trees import command


def train(source_panels, genes, output, workers=8):
    output = Path(output);output.mkdir(parents=True,exist_ok=True)
    started = time.monotonic()
    species_genes = defaultdict(set)
    per_gene = {}
    source_hashes = {}
    import hashlib
    for gene in genes:
        records = {}
        widths = set()
        for panel in source_panels:
            path = Path(panel)/'ref'/(gene+'.fasta')
            if not path.exists():
                continue
            source_hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
            for record in SeqIO.parse(path,'fasta'):
                fields = record.id.split('_',4)
                if len(fields)!=5 or fields[3].casefold().strip('.') in {'sp','spp','x','cf','aff','indet'}:
                    continue
                key = (fields[2]+'_'+fields[3]).casefold()
                width = len(record.seq);widths.add(width)
                bases = sum(c.upper() in 'ACGT' for c in str(record.seq))
                if bases<80:
                    continue
                if key not in records or bases>records[key][0]:
                    records[key]=(bases,record)
        if len(widths)>1:
            raise ValueError('Training alignments have incompatible fixed coordinates: '+gene)
        per_gene[gene] = records
        for key in records:
            species_genes[key].add(gene)
    eligible = {s for s,gs in species_genes.items() if len(gs)>=max(10,int(.6*len(genes)))}
    genera = defaultdict(list)
    for species in eligible:
        genera[species.split('_')[0]].append(species)
    chosen = {s for g,ss in genera.items() for s in sorted(ss,key=lambda x:stable('reference-library/'+x))[:2]}
    split = {g:('calibration' if stable('genus-split/'+g)%5==0 else 'test' if stable('genus-split/'+g)%5==1 else 'fit') for g in genera}
    multiplets = {g for g,ss in genera.items() if len(set(ss)&chosen)>=2}
    counts = Counter(split[g] for g in multiplets)
    if counts['fit']<10 or counts['calibration']<10 or counts['test']<5:
        raise ValueError('Insufficient independent multi-species genera for reference training: '+str(counts))
    reference_dir = output/'reference_only';reference_dir.mkdir(exist_ok=True)
    fasttree = shutil.which('fasttree') or shutil.which('FastTree')
    def build(gene):
        raw = reference_dir/(gene+'.fasta');trim = reference_dir/(gene+'.trim.fasta');tree=reference_dir/(gene+'.tre')
        rows = [v[1] for s,v in per_gene[gene].items() if s in chosen]
        if len(rows)<20:
            return None
        SeqIO.write(rows,raw,'fasta')
        command(['trimal','-in',raw,'-out',trim,'-gt','.5'])
        command([fasttree,'-gtr','-gamma','-nt',trim],tree)
        return tree
    with ThreadPoolExecutor(max_workers=workers) as pool:
        trees = [p for p in pool.map(build,genes) if p]
    if len(trees)<10:
        raise ValueError('Too few reference-only training genes')
    fit,calibration,test = [],[],[]
    cache = []
    for path in trees:
        d,taxa,_,_,_ = blocks(path)
        species = taxa[:,2]
        # Split entire genera, including their use as anchors in fitting.
        kind = np.array([split[s.split('_')[0]] for s in species])
        cache.append((path.stem,d,taxa,kind))
        candidates = [i for i,s in enumerate(species) if kind[i]=='fit' and s.split('_')[0] in multiplets]
        candidates = sorted(candidates,key=lambda i:stable(path.stem+'/'+species[i]))[:64]
        for i in candidates:
            keep = (kind=='fit') & (species!=species[i])
            x,labels = profiles(d[i,keep],d[np.ix_(keep,keep)],taxa[keep])
            if not len(x):
                continue
            correct = [j for j,(g,_) in enumerate(labels) if g==taxa[i,1]]
            if not correct:
                continue
            negatives = [j for j in np.argsort(x[:,0]) if j not in correct][:3]
            for j in correct+negatives:
                fit.append({'species':str(species[i]),'features':x[j],'label':int(j in correct)})
    species_counts = Counter(r['species'] for r in fit)
    model = make_pipeline(StandardScaler(),LogisticRegression(C=.3,class_weight='balanced',max_iter=1000))
    model.fit(np.array([r['features'] for r in fit]),np.array([r['label'] for r in fit]),
              logisticregression__sample_weight=np.array([1/species_counts[r['species']] for r in fit]))
    for gene,d,taxa,kind in cache:
        for partition,store in (('calibration',calibration),('test',test)):
            candidates = [i for i,s in enumerate(taxa[:,2]) if kind[i]==partition and s.split('_')[0] in multiplets]
            candidates = sorted(candidates,key=lambda i:stable(gene+'/'+taxa[i,2]))[:48]
            for i in candidates:
                allowed = (kind=='fit')|(kind==partition)
                keep = allowed & (taxa[:,2]!=taxa[i,2])
                x,labels = profiles(d[i,keep],d[np.ix_(keep,keep)],taxa[keep])
                if not len(x):
                    continue
                winner = int(x[:,0].argmin())
                known = float(model.predict_proba(x[winner:winner+1])[0,1])
                absent = allowed & (taxa[:,1]!=taxa[i,1])
                nx,_ = profiles(d[i,absent],d[np.ix_(absent,absent)],taxa[absent])
                nw = int(nx[:,0].argmin())
                null = float(model.predict_proba(nx[nw:nw+1])[0,1])
                store.append({'species':str(taxa[i,2]),'genus':str(taxa[i,1]),'known_score':known,
                              'known_correct':labels[winner][0]==taxa[i,1],'absent_score':null})
    def mean(rows,fn):
        groups = defaultdict(list)
        for row in rows:
            groups[row['species']].append(float(fn(row)))
        return float(np.mean([np.mean(v) for v in groups.values()]))
    curve = []
    for threshold in np.linspace(0,1,101):
        accept = mean(calibration,lambda r:r['known_correct'] and r['known_score']>=threshold)
        reject = mean(calibration,lambda r:r['absent_score']<threshold)
        curve.append({'threshold':float(threshold),'known_correct_acceptance':accept,'genus_absent_rejection':reject,'utility':(accept+reject)/2})
    selected = max(curve,key=lambda r:(r['utility'],r['threshold']))
    threshold = selected['threshold']
    bundle = {'schema':'reference-genus-full-profile-v1','model':model,'threshold':threshold,'features':FEATURES}
    joblib.dump(bundle,output/'model.joblib',compress=3)
    report = {'truth_used':False,'reference_only_trees':len(trees),'source_alignment_sha256':source_hashes,
              'reference_species':len(chosen),'reference_genera':len(genera),'multi_species_genera_by_split':dict(counts),
              'whole_genus_partitions':split,'fit_species':sorted(species_counts),'fit_rows':len(fit),
              'calibration_species':sorted({r['species'] for r in calibration}),
              'test_species':sorted({r['species'] for r in test}),
              'threshold':threshold,'calibration_choice':selected,'calibration_curve':curve,
              'test_known_correct_acceptance':mean(test,lambda r:r['known_correct'] and r['known_score']>=threshold),
              'test_genus_absent_rejection':mean(test,lambda r:r['absent_score']<threshold),
              'seconds':time.monotonic()-started,
              'scope':'Reference genus-held-out pseudoqueries; previously constructed Kew panels condition the library. Mixture truth is never read.'}
    (output/'TRAINING.json').write_text(json.dumps(report,indent=2)+'\n')
    (output/'SUCCESS').write_text('reference-only genus model complete\n')
    return report
