"""One-shot continuation of the authorized HPC workflow; no server truth upload.

Wait for this batch's collector, verify/download its archive, seal, evaluate and
report locally. Never repair scientific failures or change frozen parameters.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import time


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def dump(path,value):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
    temporary.replace(path)


def execute(args,timeout=120,check=True):
    r=subprocess.run(list(map(str,args)),text=True,capture_output=True,timeout=timeout)
    if check and r.returncode:raise RuntimeError(r.stderr[-2500:] or r.stdout[-2500:])
    return r


def ssh(command):
    return execute(['ssh','-q','-o','BatchMode=yes','-o','ConnectTimeout=15',
                    '-o','ServerAliveInterval=10','-o','ServerAliveCountMax=2','ttu',command]).stdout.strip()


def finish(run):
    repo=Path(__file__).resolve().parents[2];run=run.resolve();os.chdir(repo)
    os.environ['GIT_TERMINAL_PROMPT']='0';os.environ['GCM_INTERACTIVE']='Never'
    state=json.loads((run/'RUN_STATE.json').read_text());root=state['root'];collector=state['jobs']['collector']
    branch=state['branch'];started=time.monotonic()
    modules=[Path(__file__),Path(__file__).with_name('evaluation.py'),Path(__file__).with_name('PLAN.md'),
             Path(__file__).with_name('report.py'),Path(__file__).with_name('runner.py')]
    snapshot={str(p):sha(p) for p in modules}
    status={'state':'waiting_for_collector','pid':os.getpid(),'jobs':state['jobs'],'root':root,
            'local_source_sha256':snapshot,'started_utc':datetime.now(timezone.utc).isoformat(),
            'truth_workbook_uploaded':False}
    status_path=run/'CONTINUATION.json';dump(status_path,status)
    while True:
        if time.monotonic()-started>48*3600:raise TimeoutError('48-hour local continuation limit; HPC jobs are not cancelled')
        try:
            monitor=json.loads(ssh(f'/usr/bin/python3 -S {root}/monitor.py {root}'))
            job_state=ssh(f'sacct -X -n -P -j {collector} --format=State').splitlines()
            job_state=job_state[0].split('|')[0].strip() if job_state else 'UNKNOWN'
            status.update(state='waiting_for_collector',collector_state=job_state,remote=monitor,
                          last_check_utc=datetime.now(timezone.utc).isoformat())
            dump(status_path,status)
            if monitor['failed']:raise ValueError('Inference failed; preserve outputs and inspect sample failures')
            if job_state in ('FAILED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL') or job_state.startswith('CANCELLED'):
                raise ValueError('Collector did not complete: '+job_state)
            if job_state=='COMPLETED' and monitor['collected'] and monitor['complete']==100:break
        except ValueError:
            raise
        except Exception as exc:
            status.update(state='connection_retry',last_connection_error=str(exc)[-1200:]);dump(status_path,status)
        time.sleep(60)
    for path,digest in snapshot.items():
        if sha(path)!=digest:raise ValueError('Local evaluation code changed during inference: '+path)
    status['state']='download_verify_seal_evaluate';dump(status_path,status)
    remote_archive=root+'/compact_results.tar.gz'
    expected=ssh('sha256sum '+remote_archive).split()[0]
    archive=run/'compact_results.tar.gz';temporary=run/'compact_results.download'
    execute(['scp','-q','-B','-o','ConnectTimeout=15','-o','ServerAliveInterval=10',
             '-o','ServerAliveCountMax=2','ttu:'+remote_archive,temporary],timeout=3600)
    if sha(temporary)!=expected:raise ValueError('TTU/local archive digest mismatch')
    temporary.replace(archive)
    batch=run/'received_batch'
    if batch.exists():raise FileExistsError('Preserve existing received batch before restarting')
    batch.mkdir()
    with tarfile.open(archive) as tar:tar.extractall(batch,filter='data')
    from experiments.hierarchy100.evaluation import seal,evaluate,verify_seal
    output=repo/'experiments/hierarchy100/results/ttu100_20261022'
    seal(batch,output)
    workbook=Path('C:/Users/kudos/OneDrive - Texas Tech University/2024Fall/CodFather(TM)/SubmissionAPPS/STable.xlsx')
    evaluate(batch,output,workbook)
    # XLSX truth extraction uses this bundled runtime; scientific plotting uses
    # the existing project environment, which has Matplotlib installed.
    plot_python=repo.parent/'.venv312/Scripts/python.exe'
    execute([plot_python,'-m','experiments.hierarchy100.report','--batch',batch,'--output',output],timeout=1800)
    verify_seal(output,batch)
    artifacts=output/'artifacts';artifacts.mkdir()
    for name in ('payload.tar.gz','tests_prelaunch.log','JOBS.json'):
        shutil.copy2(run/name,artifacts/name)
    if archive.stat().st_size<=40*1024*1024:
        shutil.copy2(archive,artifacts/'compact_results.tar.gz')
    else:
        parts=[]
        with archive.open('rb') as source:
            while True:
                block=source.read(20*1024*1024)
                if not block:break
                p=artifacts/f'compact_results.tar.gz.part{len(parts)+1:03d}';p.write_bytes(block)
                parts.append({'file':p.name,'sha256':sha(p)})
        dump(artifacts/'ARCHIVE_PARTS.json',{'concatenate_in_order':parts,'complete_archive_sha256':expected})
    accounting=ssh('sacct -P -j '+','.join(str(v) for v in state['jobs'].values())+' --format=JobID,JobName,State,Submit,Start,End,Elapsed,TotalCPU,AllocCPUS,MaxRSS,ExitCode')
    (artifacts/'SLURM_ACCOUNTING.psv').write_text(accounting+'\n',encoding='utf-8')
    dump(output/'DOWNLOAD_AUDIT.json',{'remote_archive_sha256':expected,'local_archive_sha256':sha(archive),
         'safe_extraction_filter':'data','all_100_raw_sets_verified':True,'local_sources_unchanged':True,
         'truth_workbook_uploaded':False,'visual_inspection_by_agent_pending':True})
    files=[p for p in output.rglob('*') if p.is_file()]
    dump(output/'ARTIFACTS_SHA256.json',{'files':{p.relative_to(output).as_posix():sha(p) for p in sorted(files)}})
    current=execute(['git','branch','--show-current']).stdout.strip()
    if current!=branch:raise ValueError('Results produced; branch changed, so automatic Git publication is skipped')
    # Stage only this run's result package; never include unrelated working files.
    execute(['git','add',output.relative_to(repo)])
    execute(['git','diff','--cached','--check'])
    execute(['git','commit','-q','--only','-m','Record all 100 fresh TTU order-to-genus predictions and fixed-rule validation',output.relative_to(repo)])
    execute(['git','push','origin','HEAD:refs/heads/'+branch],timeout=1800)
    head=execute(['git','rev-parse','HEAD']).stdout.strip()
    remote_head=execute(['git','ls-remote','--heads','origin','refs/heads/'+branch]).stdout.split()[0]
    if head!=remote_head:raise ValueError('Remote result commit not confirmed')
    status.update(state='complete',report=str(output/'REPORT.md'),commit=head,
                  finished_utc=datetime.now(timezone.utc).isoformat(),visual_inspection_by_agent_pending=True)
    dump(status_path,status)
    state.update(state='all 100 evaluated and published',remaining=[],final_commit=head)
    dump(run/'RUN_STATE.json',state)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);args=p.parse_args()
    try:finish(args.run)
    except Exception as exc:
        path=args.run.resolve()/'CONTINUATION.json'
        value=json.loads(path.read_text()) if path.exists() else {}
        value.update(state='needs_attention',error=repr(exc));dump(path,value)
        raise
