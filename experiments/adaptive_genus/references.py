"""Lazy Kew downloads, target-locus QC, and reusable fixed reference coordinates."""
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

from .phylogeny import Cell


def dump(path, value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


class ReferenceBank:
    def __init__(self, index, client, backbone, genes, families, output, workers=8, minimum_fraction=.6):
        from sprout_refbuilder.models import Taxonomy
        self.index,self.client,self.backbone,self.genes,self.families,self.output = index,client,Path(backbone),sorted(genes),sorted(families),Path(output)
        self.workers = workers
        self.orders = {tax[0] for tax in index.taxonomy.values() if tax[1] in families}
        self.minimum = max(10,int(len(genes)*minimum_fraction))
        self.output.mkdir(parents=True,exist_ok=True)
        self.entries,self.by_tip,self.by_genus = {},defaultdict(list),defaultdict(list)
        by_species = defaultdict(list)
        for tip,tax in index.taxonomy.items():
            if tax[1] in families:
                by_species[(tax[2]+' '+tax[3]).casefold()].append(tip)
        for record in client.records():
            if not record.recovery_file:
                continue
            tips = by_species.get(' '.join(record.species.replace('_',' ').split()[:2]).casefold(),[])
            if not tips:
                continue
            exact = [t for t in tips if index.taxonomy[t][4].casefold()==record.sequence_id.casefold()]
            points = exact or tips
            for tip in points:
                tax = Taxonomy(*index.taxonomy[tip][:4])
                key = record.recovery_file
                if key not in self.entries:
                    self.entries[key] = {'record':record,'taxonomy':tax}
                    self.by_genus[tax.genus].append(key)
                self.by_tip[tip].append(key)
        for mapping in (self.by_tip,self.by_genus):
            for name,keys in mapping.items():
                mapping[name] = sorted(set(keys),key=self.priority)
        self.quality = {}
        self.aligned = {g:{} for g in self.genes}
        self.fetch_seconds,self.qc_seconds,self.align_seconds = 0.,0.,0.
        self.downloaded_bytes = 0
        self.recovery_hashes = {}

    def priority(self,key):
        row = self.entries[key]['record']
        kind = row.sequence_type.casefold()
        return (0 if kind=='annotated_genome' else 1 if 'genome' in kind else 2,
                0 if 'paftol' in row.project.casefold() else 1,row.species.casefold(),row.sequence_id)

    def load(self, keys):
        from sprout_refbuilder.fasta import read_fasta,parse_kew_gene_id
        from sprout_refbuilder.clade import CladeRecovery,_evaluate_locus
        keys = sorted(set(keys)-set(self.quality))
        if not keys:
            return
        start = time.monotonic()
        def fetch(key):
            entry = self.entries[key]
            expected = self.client.cache_dir/'recoveries'/key
            existed = expected.exists()
            path = self.client.fetch_recovery(entry['record'])
            groups = defaultdict(list)
            for seq in read_fasta(path):
                gene = parse_kew_gene_id(seq)
                if gene in self.genes:
                    groups[gene].append(seq)
            return key,path,dict(groups),existed
        with ThreadPoolExecutor(max_workers=4) as pool:
            fetched = list(pool.map(fetch,keys))
        self.fetch_seconds += time.monotonic()-start
        records = []
        for key,path,groups,existed in fetched:
            entry = self.entries[key]
            self.recovery_hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
            self.downloaded_bytes += 0 if existed else path.stat().st_size
            recovery = CladeRecovery(entry['record'],entry['taxonomy'],path,len(groups))
            records.append((key,recovery,groups))
        # Only the loci used by this predictor are QC'd; every locus still uses
        # its complete corrected backbone for taxonomic/sequence checks.
        jobs = [(g,self.backbone/(g+'.fasta'),'genus',.5,
                 [(r,groups.get(g,[])) for _,r,groups in records]) for g in self.genes]
        start = time.monotonic()
        with ProcessPoolExecutor(max_workers=self.workers) as pool:
            results = list(pool.map(_evaluate_locus,jobs))
        self.qc_seconds += time.monotonic()-start
        accepted,all_qc = defaultdict(dict),[]
        ids = {r.record.sequence_id:key for key,r,_ in records}
        if len(ids)!=len(records):
            raise ValueError('Duplicate recovery sequence IDs need explicit disambiguation')
        for gene,rows,qc in results:
            all_qc.extend(qc)
            for seq,identifier in rows:
                accepted[ids[identifier]][gene] = seq
        for key,_,_ in records:
            self.quality[key] = {'accepted_loci':len(accepted[key]),'minimum_loci':self.minimum,
                                 'passed':len(accepted[key])>=self.minimum,
                                 'accepted':accepted[key]}
        with (self.output/'qc.jsonl').open('a') as f:
            for row in all_qc:
                f.write(json.dumps(row)+'\n')

    def usable(self,key):
        return self.quality.get(key,{}).get('passed',False)

    def anchors(self,cells,max_attempts=6):
        requests = []
        for cell in cells:
            for region in self.index.representative_regions(cell):
                proposals,seen = [],set()
                ranked = self.index.medoid_order(region,self.by_tip)
                seen_genera,first,later = set(),[],[]
                for tip in ranked:
                    genus = self.index.taxonomy[tip][2]
                    if genus not in seen_genera:
                        first.append(tip);seen_genera.add(genus)
                    else:
                        later.append(tip)
                ordered = first+later
                for backup in range(2):
                    for tip in ordered:
                        keys = self.by_tip[tip]
                        if backup>=len(keys):
                            continue
                        key = keys[backup]
                        if key not in seen:
                            proposals.append((tip,key));seen.add(key)
                requests.append({'cell':cell,'proposals':proposals[:max_attempts],'chosen':None})
        for attempt in range(max_attempts):
            missing = [r for r in requests if r['chosen'] is None and len(r['proposals'])>attempt]
            self.load([r['proposals'][attempt][1] for r in missing])
            for request in missing:
                tip,key = request['proposals'][attempt]
                if self.usable(key):
                    request['chosen'] = (tip,key)
        mapping,chosen,details = defaultdict(set),set(),[]
        for request in requests:
            value = request['chosen']
            detail = {'cell':request['cell'].key,'proposed_records':len(request['proposals']),
                      'quality_unresolved':value is None}
            if value:
                tip,key = value
                chosen.add(key)
                header = self.header(key)
                mapping[header].add(request['cell'].key)
                detail.update(tip=tip,recovery=key,genus=self.entries[key]['taxonomy'].genus,
                              species=self.entries[key]['record'].species,
                              accepted_loci=self.quality[key]['accepted_loci'])
            details.append(detail)
        return sorted(chosen),{h:sorted(v) for h,v in mapping.items()},details

    def genus_panel(self,genera,representatives=2,max_attempts=8):
        chosen,omitted = set(),[]
        requests = {g:list(self.by_genus.get(g,[])) for g in sorted(genera)}
        selected,species = defaultdict(list),defaultdict(set)
        for attempt in range(max_attempts):
            trials = {}
            for genus,keys in requests.items():
                if len(selected[genus])>=representatives:
                    continue
                options = [k for k in keys if k not in selected[genus]
                           and self.entries[k]['record'].species not in species[genus]
                           and not (k in self.quality and not self.usable(k))]
                if options:
                    known = sorted((k for k in options if self.usable(k)),
                                   key=lambda k:(-self.quality[k]['accepted_loci'],self.priority(k)))
                    trials[genus] = (known or options)[0]
            self.load(list(trials.values()))
            for genus,key in trials.items():
                if self.usable(key):
                    selected[genus].append(key)
                    species[genus].add(self.entries[key]['record'].species)
                    chosen.add(key)
        for genus in genera:
            if not selected[genus]:
                omitted.append(genus)
        return sorted(chosen),sorted(omitted),{g:len(v) for g,v in selected.items()}

    def header(self,key):
        row = self.entries[key]
        return row['taxonomy'].sprout_prefix+'_KEW_'+row['record'].sequence_id

    def metadata(self,key):
        entry=self.entries[key]
        tax=asdict(entry['taxonomy'])
        epithet=tax.pop('species')
        return {'recovery':key,'header':self.header(key),
                'accepted_loci':self.quality[key]['accepted_loci'],
                **entry['record'].to_dict(),**tax,'species_epithet':epithet}

    def panel(self,keys,output):
        from sprout_refbuilder.clade import _align_locus
        from sprout_refbuilder.fasta import write_fasta
        keys = sorted(set(keys))
        unknown = [key for key in keys if not self.usable(key)]
        if unknown:
            raise ValueError('Panel contains failed/unchecked references')
        output = Path(output);output.mkdir(parents=True,exist_ok=True)
        todo = [key for key in keys if not any(self.header(key) in self.aligned[g] for g in self.genes)]
        cache = output/'new_alignment';cache.mkdir(exist_ok=True)
        jobs = []
        for gene in self.genes:
            seqs = [self.quality[k]['accepted'][gene] for k in todo if gene in self.quality[k]['accepted']]
            if seqs:
                # Sample parent context stays fixed across rounds, coordinates
                # stay fixed, and each downloaded reference is projected once.
                jobs.append((gene,self.backbone/(gene+'.fasta'),seqs,True,cache/(gene+'.fasta'),
                             'genus',{order.casefold() for order in self.orders},
                             {f.casefold() for f in self.families},
                             {genus.casefold() for genus in self.by_genus},False))
        start = time.monotonic()
        with ProcessPoolExecutor(max_workers=self.workers) as pool:
            list(pool.map(_align_locus,jobs))
        self.align_seconds += time.monotonic()-start
        from sprout_refbuilder.fasta import read_fasta
        for gene,*_ in jobs:
            self.aligned[gene].update({s.id:s for s in read_fasta(cache/(gene+'.fasta'))})
        allowed = {self.header(k) for k in keys}
        reference_dir = output/'ref';reference_dir.mkdir(exist_ok=True)
        valid = []
        for gene in self.genes:
            rows = [s for h,s in self.aligned[gene].items() if h in allowed]
            if len(rows)>=4:
                write_fasta(reference_dir/(gene+'.fasta'),rows)
                valid.append(gene)
        report = {'genes':valid,'references':len(keys),
                  'selection':[self.metadata(k) for k in keys],
                  'new_references_aligned':len(todo),'fetch_seconds_cumulative':self.fetch_seconds,
                  'qc_seconds_cumulative':self.qc_seconds,'reference_alignment_seconds_cumulative':self.align_seconds,
                  'new_download_bytes_cumulative':self.downloaded_bytes,
                  'quality_records_checked':len(self.quality),'minimum_accepted_target_genes':self.minimum}
        dump(output/'panel.json',report)
        return reference_dir,report
