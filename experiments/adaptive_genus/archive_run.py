"""Package completed HPC metadata without changing frozen inference files.

Run after the collector. Large sequence libraries, raw reads and inferred trees
stay on the host; product/recovery hashes and exact software sources are kept.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile


def package(root, job_ids):
    root=Path(root).resolve();base=root.parent
    if not (root/'PREDICTIONS_COMPLETE').exists():
        raise ValueError('Collector has not finished')
    audit=json.loads((root/'AUDIT.json').read_text())
    if not audit['all_10_complete'] or not audit['code_model_inputs_products_verified']:
        raise ValueError('Incomplete audited batch')
    output=base/('reproducibility_'+root.name+'.tar.gz')
    if output.exists():
        raise FileExistsError(output)
    files={}
    for name in ('train.py','train_code.tar.gz','TRAIN_CODE_SHA.json','RUN_STATE.json','PROBE.json','TOOLCHAIN.json','toolchain.py','EVALUATION_PLAN.md','collect_parallel.py','collect_parallel_v1.sbatch','collect_parallel_strict_route.sbatch'):
        p=base/name
        if p.exists():
            files[name]=p
    for name in ('sources.tar.gz','model/model.joblib','model/TRAINING.json','prepare.sbatch','array.sbatch','collect.sbatch'):
        files['inference/'+name]=root/name
    if (root/'COLLECTION_IMPLEMENTATION.json').exists():
        files['inference/COLLECTION_IMPLEMENTATION.json']=root/'COLLECTION_IMPLEMENTATION.json'
    for name in ('coarse_reuse_manifest.json','POSTHOC_COMPLETION_PLAN.md'):
        if (root/name).exists():
            files['inference/'+name]=root/name
    for name in ('PILOT.json','pilot.py'):
        p=base/'pilot_02'/name
        if p.exists():
            files['pilot_02/'+name]=p
    for p in (base/'logs').glob('*'):
        if p.is_file():
            files['earlier_logs/'+p.name]=p
    for p in (root/'logs').glob('*'):
        if p.is_file():
            files['inference/logs/'+p.name]=p
    for directory in sorted((root/'samples').glob('mix_*')):
        q=directory/'reference_bank/qc.jsonl'
        if q.exists():
            files['inference/'+str(q.relative_to(root))]=q
        for q in directory.glob('*/trees/query_ids.json'):
            files['inference/'+str(q.relative_to(root))]=q
        q=directory/'final_trees/query_ids.json'
        if q.exists():
            files['inference/'+str(q.relative_to(root))]=q
    # Preserve Slurm's actual elapsed/resource status, including failed pilots.
    result=subprocess.run(['sacct','-P','-j',','.join(map(str,job_ids)),
                           '--format=JobID,JobName,State,Submit,Start,End,Elapsed,TotalCPU,AllocCPUS,MaxRSS,ExitCode'],
                          check=True,text=True,capture_output=True)
    accounting=base/('SLURM_ACCOUNTING_'+root.name+'.psv');accounting.write_text(result.stdout)
    files[accounting.name]=accounting
    manifest={'packaged_utc':datetime.now(timezone.utc).isoformat(),'inference_root':str(root),
              'job_ids':job_ids,'truth_included':False,
              'artifact_sha256':{name:hashlib.sha256(p.read_bytes()).hexdigest() for name,p in sorted(files.items())}}
    record=base/('REPRODUCIBILITY_'+root.name+'.json');record.write_text(json.dumps(manifest,indent=2)+'\n')
    files[record.name]=record
    with tarfile.open(output,'w:gz') as archive:
        for name,p in sorted(files.items()):
            archive.add(p,arcname=name)
    for p in (root/'compact_results.tar.gz',output):
        print(hashlib.sha256(p.read_bytes()).hexdigest(),p)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True,type=Path)
    parser.add_argument('--jobs',required=True,nargs='+',type=int)
    args=parser.parse_args()
    package(args.root,args.jobs)
