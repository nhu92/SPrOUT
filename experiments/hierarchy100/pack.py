"""Build a self-contained, auditable TTU batch payload from existing sources."""
from pathlib import Path
import hashlib
import io
import json
import tarfile


def pack():
    repo=Path(__file__).resolve().parents[2];base=repo/'experiments/hierarchy100'
    run=base/'runs/ttu_20261022';run.mkdir(parents=True,exist_ok=True)
    sources=repo/'experiments/reference_atlas/runs/umd_real11_20261016/payload_v2/sources'
    if not (sources/'sprout_evidence').is_dir():
        raise FileNotFoundError('Frozen shared family/QC implementation is required')
    files={}
    for package in ('adaptive_genus','reference_atlas','hierarchy100'):
        for p in (repo/'experiments'/package).glob('*.py'):
            files['code/experiments/'+package+'/'+p.name]=p.read_bytes()
    files['code/experiments/__init__.py']=b''
    for p in (repo/'experiments/adaptive_genus').glob('policy*.json'):
        files['code/experiments/adaptive_genus/'+p.name]=p.read_bytes()
    for name in ('transfer.py','fuse.py','fusion_policy.json'):
        files['code/'+name]=(repo/'experiments/reference_atlas'/name).read_bytes()
    for p in sources.rglob('*'):
        if p.is_file() and '__pycache__' not in p.parts:
            files['sources/'+p.relative_to(sources).as_posix()]=p.read_bytes()
    modifications=[]
    path='sources/02_exon_trees.py';old=files[path]
    pattern=b'gene, threads, input_exon_dir, ref_dir, output_phylo, log_file, min_size,'
    if old.count(pattern)!=1:
        raise ValueError('Unexpected tree worker implementation')
    files[path]=old.replace(pattern,b'gene, 1, input_exon_dir, ref_dir, output_phylo, log_file, min_size,')
    modifications.append({'path':path,'original_sha256':hashlib.sha256(old).hexdigest(),
                          'change':'One MAFFT thread per independent gene worker; no nested CPU oversubscription'})
    path='sources/MrKnowItAll/src/sprout_refbuilder/kew.py';old=files[path]
    files[path]+=b'''\n# Operational guard: only duplicate writers of the same recovery wait.\n_original_fetch_recovery = KewClient.fetch_recovery\ndef _locked_fetch_recovery(self, record):\n    import hashlib, fcntl\n    token = hashlib.sha256(str(record.recovery_file).encode()).hexdigest()\n    locks = self.cache_dir / '.recovery_locks'\n    locks.mkdir(exist_ok=True)\n    with (locks / token).open('a+b') as handle:\n        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)\n        try:\n            return _original_fetch_recovery(self, record)\n        finally:\n            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)\nKewClient.fetch_recovery = _locked_fetch_recovery\n'''
    modifications.append({'path':path,'original_sha256':hashlib.sha256(old).hexdigest(),
                          'change':'Advisory lock per recovery file; independent records remain concurrent'})
    for rank,prefix in [('family',repo/'experiments/reference_atlas/models/full_profile_transfer_v1'),
                        ('genus',repo/'experiments/adaptive_genus/results/ttu10_20261017/shared_model')]:
        files['models/'+rank+'/model.joblib']=(prefix/'model.joblib').read_bytes()
        if (prefix/'TRAINING.json').exists():files['models/'+rank+'/TRAINING.json']=(prefix/'TRAINING.json').read_bytes()
    files['code/model.joblib']=files['models/family/model.joblib']
    for name in ('PLAN.md','REQUEST.json','cohort.json','prepare.sbatch','array.sbatch','collect.sbatch'):
        files[name]=(base/name if name=='PLAN.md' else run/name).read_bytes()
        if name.endswith('.sbatch'):
            files[name]=files[name].replace(b'\r\n',b'\n')
    manifest={'operational_source_changes':modifications,
              'artifact_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in sorted(files.items())}}
    files['PAYLOAD.json']=(json.dumps(manifest,indent=2)+'\n').encode()
    with tarfile.open(run/'payload.tar.gz','w:gz') as archive:
        for name,data in sorted(files.items()):
            info=tarfile.TarInfo(name);info.size=len(data);info.mode=0o644
            archive.addfile(info,io.BytesIO(data))
    print('Packed',len(files),'files;',run/'payload.tar.gz')


if __name__=='__main__':pack()
