"""Reduced-panel alignments and full exon trees, with explicit exclusions."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord


def command(args, output=None):
    if output:
        with Path(output).open('w') as handle:
            result = subprocess.run(list(map(str,args)),stdout=handle,stderr=subprocess.PIPE,text=True,timeout=1800)
    else:
        result = subprocess.run(list(map(str,args)),stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=1800)
    if result.returncode:
        raise RuntimeError(' '.join(map(str,args[:3]))+': '+result.stderr[-1500:])
    return result


def infer_exons(reference_dir, exon_dir, genes, output, workers=8, minimum_bases=80):
    reference_dir,exon_dir,output = map(Path,(reference_dir,exon_dir,output))
    output.mkdir(parents=True,exist_ok=True)
    paths,excluded,query_map = [],[],[]
    fasttree = shutil.which('fasttree') or shutil.which('FastTree')
    if not fasttree:
        raise FileNotFoundError('FastTree is not on PATH')
    start = time.monotonic()
    files_by_gene = {g:sorted(p for p in exon_dir.glob('*.fasta')
                             if re.search(r'(?:^|_)'+re.escape(g)+r'(?:_|$)',p.stem)) for g in genes}
    def process(gene):
        good,bad,mapping = [],[],[]
        ref = reference_dir/(gene+'.fasta')
        if not ref.exists():
            return good,bad,mapping
        trimmed = output/(gene+'_reference_trimmed.fasta')
        command(['trimal','-in',ref,'-out',trimmed,'-gt','.5'])
        for number,source in enumerate(files_by_gene[gene],1):
            name = gene+'_exon_'+str(number)
            query = []
            for i,record in enumerate(SeqIO.parse(source,'fasta')):
                if sum(c.upper() in 'ACGT' for c in str(record.seq))<minimum_bases:
                    continue
                identity = 'QUERY_'+gene+'_'+str(number)+'_'+str(i)
                query.append(SeqRecord(record.seq,id=identity,description=''))
                mapping.append({'query':identity,'source':str(source),'original_id':record.id})
            if not query:
                continue
            qpath,alignment,tree = (output/(name+suffix) for suffix in ('_query.fasta','_aligned.fasta','.tre'))
            SeqIO.write(query,qpath,'fasta')
            try:
                command(['mafft','--quiet','--thread','1','--preservecase','--adjustdirection','--addfragments',qpath,'--keeplength',trimmed],alignment)
                rows = list(SeqIO.parse(alignment,'fasta'))
                for row in rows:
                    row.id = row.id.removeprefix('_R_');row.description = ''
                projected = [row for row in rows if not row.id.startswith('QUERY_')
                             or sum(c.upper() in 'ACGT' for c in str(row.seq))>=minimum_bases]
                if not any(row.id.startswith('QUERY_') for row in projected):
                    raise ValueError('No query retains minimum aligned ACGT overlap')
                SeqIO.write(projected,alignment,'fasta')
                command([fasttree,'-gtr','-gamma','-nt',alignment],tree)
                parsed = Phylo.read(tree,'newick')
                labels = [n.name for n in parsed.get_terminals()]
                if len(labels)!=len(set(labels)) or not any(n.startswith('QUERY_') for n in labels):
                    raise ValueError('Missing query or duplicate names')
                good.append(tree)
            except (RuntimeError,ValueError) as error:
                bad.append({'tree':str(tree),'gene':gene,'error':str(error)})
        return good,bad,mapping
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for good,bad,mapping in pool.map(process,genes):
            paths.extend(good);excluded.extend(bad);query_map.extend(mapping)
    record = {'trees':[p.name for p in paths],'valid_genes':len({p.name.split('_')[0] for p in paths}),
              'excluded_trees':excluded,'seconds':time.monotonic()-start,
              'query_contigs':len(query_map),'reference_trim_rule':'trimAl -gt .5 on references before adding queries',
              'query_alignment':'mafft --addfragments --keeplength, one thread per independent gene worker'}
    (output/'trees.json').write_text(json.dumps(record,indent=2)+'\n')
    (output/'query_ids.json').write_text(json.dumps(query_map,indent=2)+'\n')
    return sorted(paths),record
