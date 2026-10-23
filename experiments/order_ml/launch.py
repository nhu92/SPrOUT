"""Local packaging and one-time Slurm submission; never send truth to TTU."""
import argparse
import io
from pathlib import Path
import subprocess
import tarfile

from .common import dump, sha

ROOT = '/lustre/scratch/nhu/202609/ttu_order_ml100_20261023'
SOURCE = '/lustre/scratch/nhu/202609/ttu_hierarchy100_20261022'
REFERENCE = '/lustre/scratch/nhu/202609/ttu_100mixes10_20261013/order_reference'


def execute(args):
    r = subprocess.run(list(map(str, args)), text=True, capture_output=True, timeout=300)
    if r.returncode:
        raise RuntimeError(r.stderr or r.stdout)
    return r.stdout.strip()


def ssh(command):
    return execute(['ssh', '-q', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15', 'ttu', command])


def launch(directory):
    repo = Path(__file__).resolve().parents[2]; directory = directory.resolve(); directory.mkdir(parents=True, exist_ok=True)
    if (directory/'JOBS.json').exists():
        raise FileExistsError('Jobs already submitted; inspect rather than duplicate')
    files = [Path('experiments/order_ml')/name for name in ('__init__.py', 'common.py', 'model.py', 'runner.py')]
    files += [Path('experiments/reference_atlas/transfer.py')]
    request = {'source': SOURCE, 'reference': REFERENCE, 'root': ROOT, 'primary': 'z_plus_ml_specificity',
               'plan_sha256': sha(repo/'experiments/order_ml/PLAN.md'),
               'branch': execute(['git', '-C', repo, 'branch', '--show-current'])}
    dump(directory/'REQUEST.json', request)
    scripts = {}
    for name, command, minutes, memory, array in (
        ('prepare', 'prepare', 120, '16G', ''),
        ('array', 'run --index "$SLURM_ARRAY_TASK_ID"', 60, '16G', '#SBATCH --array=0-99\n'),
        ('collect', 'collect', 60, '16G', '')):
        scripts[name] = f'''#!/bin/bash
#SBATCH --job-name=oml_{name}
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem={memory}
#SBATCH --time={minutes}
{array}#SBATCH --output={ROOT}/logs/{name}_%A_%a.out
#SBATCH --error={ROOT}/logs/{name}_%A_%a.err
set -euo pipefail
export PATH=/home/nhu/miniconda3/envs/sprout/bin:$PATH
export PYTHONPATH={ROOT}/code
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python -m experiments.order_ml.runner {command} --root {ROOT}
'''
        (directory/(name+'.sbatch')).write_text(scripts[name], encoding='utf-8', newline='\n')
    payload = directory/'payload.tar.gz'
    with tarfile.open(payload, 'w:gz') as tar:
        for path in files:
            tar.add(repo/path, arcname='code/'+path.as_posix())
        tar.add(repo/'experiments/order_ml/PLAN.md', arcname='PLAN.md')
        tar.add(directory/'REQUEST.json', arcname='REQUEST.json')
        for name, script in scripts.items():
            content = script.encode(); info = tarfile.TarInfo(name+'.sbatch'); info.size = len(content)
            tar.addfile(info, io.BytesIO(content))
    dump(directory/'PAYLOAD.json', {'sha256': sha(payload), 'files': {p.as_posix(): sha(repo/p) for p in files},
                                  'truth_uploaded': False})
    ssh(f'test ! -e {ROOT} && mkdir -p {ROOT}/logs')
    execute(['scp', '-q', '-B', '-o', 'ConnectTimeout=15', payload, 'ttu:'+ROOT+'/payload.tar.gz'])
    remote = ssh(f'sha256sum {ROOT}/payload.tar.gz').split()[0]
    if remote != sha(payload):
        raise ValueError('Payload transfer mismatch')
    ssh(f'tar -xzf {ROOT}/payload.tar.gz -C {ROOT}')
    jobs = {}
    def submit(name, dependency=''):
        answer = ssh(f'sbatch --parsable {dependency} {ROOT}/{name}.sbatch')
        job = answer.split(';')[0]
        if not job.isdigit():
            raise ValueError('Unexpected job ID '+answer)
        jobs[name] = job
        dump(directory/'JOBS.json', {'root': ROOT, 'jobs': jobs})
        return job
    prepare = submit('prepare')
    array = submit('array', '--dependency=afterok:'+prepare)
    submit('collect', '--dependency=afterany:'+array)
    dump(directory/'RUN_STATE.json', {'root': ROOT, 'jobs': jobs, 'branch': request['branch'],
         'state': 'submitted', 'remaining': ['reference fitting', '100 predictions', 'sealed local evaluation']})
    print(jobs)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--run', type=Path, required=True)
    launch(p.parse_args().run)
