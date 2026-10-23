"""One-shot local continuation of this authorized batch, with sealed evaluation."""
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import time

from .common import read, dump, sha
from .launch import execute, ssh


def finish(run):
    repo = Path(__file__).resolve().parents[2]; os.chdir(repo); run = run.resolve()
    state = read(run/'RUN_STATE.json'); root = state['root']; jobs = state['jobs']
    branch = state['branch']; started = time.monotonic()
    files = list(Path(__file__).parent.glob('*.py'))+[Path(__file__).with_name('PLAN.md')]
    snapshot = {str(p): sha(p) for p in files}
    status_path = run/'CONTINUATION.json'
    status = {'state': 'waiting_for_slurm', 'pid': os.getpid(), 'root': root, 'jobs': jobs,
              'started_utc': datetime.now(timezone.utc).isoformat(), 'truth_uploaded': False,
              'local_source_sha256': snapshot}
    dump(status_path, status)
    while True:
        if time.monotonic()-started > 48*3600:
            raise TimeoutError('48-hour continuation limit; cluster jobs are preserved')
        try:
            accounting = ssh('sacct -X -n -P -j '+','.join(jobs.values())+' --format=JobIDRaw,State,ExitCode')
            states = {r.split('|')[0]: r.split('|')[1] for r in accounting.splitlines() if '|' in r}
            status.update(job_states=states, last_check_utc=datetime.now(timezone.utc).isoformat())
            dump(status_path, status)
            failed = {job: value for job, value in states.items()
                      if value in ('FAILED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL', 'BOOT_FAIL') or value.startswith('CANCELLED')}
            if failed:
                raise ValueError('Batch failed; preserve outputs and inspect: '+str(failed))
            if states.get(jobs['collect']) == 'COMPLETED':
                break
        except ValueError:
            raise
        except Exception as exc:
            status.update(state='connection_retry', last_connection_error=str(exc)[-1500:]); dump(status_path, status)
        time.sleep(60)
    for path, digest in snapshot.items():
        if sha(path) != digest:
            raise ValueError('Local evaluation source changed while waiting: '+path)
    status['state'] = 'download_verify_seal_evaluate'; dump(status_path, status)
    remote_archive = root+'/results.tar.gz'
    expected = ssh('sha256sum '+remote_archive).split()[0]
    archive = run/'results.tar.gz'; temporary = run/'results.download'
    result = subprocess.run(['scp', '-q', '-B', '-o', 'ConnectTimeout=15', '-o', 'ServerAliveInterval=10',
                             '-o', 'ServerAliveCountMax=2', 'ttu:'+remote_archive, str(temporary)],
                            capture_output=True, text=True, timeout=3600)
    if result.returncode:
        raise RuntimeError(result.stderr[-1500:])
    if sha(temporary) != expected:
        raise ValueError('Remote/local archive digest mismatch')
    temporary.replace(archive)
    batch = run/'received_batch'
    if batch.exists():
        raise FileExistsError('Preserve previous download before restarting')
    batch.mkdir()
    with tarfile.open(archive) as tar:
        tar.extractall(batch, filter='data')
    from .evaluation import seal, evaluate, verify_seal
    from .report import report
    output = repo/'experiments/order_ml/results/ttu100_20261023'
    original = repo/'experiments/hierarchy100/results/ttu100_20261022'
    seal(batch, output, original)
    workbook = Path('C:/Users/kudos/OneDrive - Texas Tech University/2024Fall/CodFather(TM)/SubmissionAPPS/STable.xlsx')
    evaluate(batch, output, workbook); report(output); verify_seal(batch, output)
    artifacts = output/'artifacts'; artifacts.mkdir()
    for name in ('payload.tar.gz', 'PAYLOAD.json', 'JOBS.json', 'tests_prelaunch.log'):
        shutil.copy2(run/name, artifacts/name)
    if archive.stat().st_size <= 40*1024*1024:
        shutil.copy2(archive, artifacts/'results.tar.gz')
    else:
        parts = []
        with archive.open('rb') as f:
            while True:
                content = f.read(20*1024*1024)
                if not content:
                    break
                path = artifacts/f'results.tar.gz.part{len(parts)+1:03d}'; path.write_bytes(content)
                parts.append({'file': path.name, 'sha256': sha(path)})
        dump(artifacts/'ARCHIVE_PARTS.json', {'concatenate_in_order': parts, 'complete_archive_sha256': expected})
    accounting = ssh('sacct -P -j '+','.join(jobs.values())+' --format=JobID,JobName,State,Submit,Start,End,Elapsed,TotalCPU,AllocCPUS,MaxRSS,ExitCode')
    (artifacts/'SLURM_ACCOUNTING.psv').write_text(accounting+'\n', encoding='utf-8')
    dump(output/'DOWNLOAD_AUDIT.json', {'remote_archive_sha256': expected, 'local_archive_sha256': sha(archive),
         'safe_extraction_filter': 'data', 'all100_raw_scores_reproduced': True, 'truth_uploaded': False,
         'source_sha256_unchanged': True})
    dump(output/'ARTIFACTS_SHA256.json', {'files': {p.relative_to(output).as_posix(): sha(p)
         for p in sorted(output.rglob('*')) if p.is_file()}})
    status.update(state='evaluated', report=str(output/'REPORT.md')); dump(status_path, status)
    if execute(['git', 'branch', '--show-current']) != branch:
        raise ValueError('Evaluation complete; branch changed, automatic publication skipped')
    os.environ['GIT_TERMINAL_PROMPT'] = '0'; os.environ['GCM_INTERACTIVE'] = 'Never'
    execute(['git', 'add', output.relative_to(repo)])
    execute(['git', 'diff', '--cached', '--check'])
    execute(['git', 'commit', '-q', '--only', '-m', 'Record fixed order ML and Z rescue results on 100 mixtures', output.relative_to(repo)])
    execute(['git', 'push', 'origin', 'HEAD:refs/heads/'+branch])
    head = execute(['git', 'rev-parse', 'HEAD'])
    remote_head = execute(['git', 'ls-remote', '--heads', 'origin', 'refs/heads/'+branch]).split()[0]
    if head != remote_head:
        raise ValueError('Remote result commit not confirmed')
    status.update(state='complete', commit=head, finished_utc=datetime.now(timezone.utc).isoformat())
    dump(status_path, status)
    state.update(state='all100 evaluated and published', final_commit=head, remaining=[]); dump(run/'RUN_STATE.json', state)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--run', type=Path, required=True)
    a = p.parse_args()
    try:
        finish(a.run)
    except Exception as exc:
        path = a.run.resolve()/'CONTINUATION.json'; record = read(path) if path.exists() else {}
        record.update(state='needs_attention', error=repr(exc), failed_utc=datetime.now(timezone.utc).isoformat())
        dump(path, record)
        raise
