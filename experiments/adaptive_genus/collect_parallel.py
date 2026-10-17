"""Operational collector: parallelize hash reads of completed frozen runs.

Inference sources and predictions are unchanged. The original frozen collector
still assembles the results; only its read-only verifier is replaced with an
equivalent streaming SHA-256 verifier. Archive this wrapper with its audit.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys


def sha_stream(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def verify_parallel(hashes,workers=16):
    def check(item):
        path,expected=item
        if sha_stream(path)!=expected:
            raise ValueError('Frozen file changed: '+str(path))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for _ in pool.map(check,hashes.items()):
            pass


def collect(root,workers):
    root=Path(root).resolve()
    sys.path[:0]=[str(root/'code'),str(root/'sources/MrKnowItAll/src'),str(root/'sources')]
    from experiments.adaptive_genus import runner
    if Path(runner.__file__).resolve()!=root/'code/experiments/adaptive_genus/runner.py':
        raise ValueError('Must import the frozen inference collector')
    record={'started_utc':datetime.now(timezone.utc).isoformat(),'workers':workers,
            'wrapper_sha256':sha_stream(__file__),'original_collector':str(Path(runner.__file__).resolve()),
            'change':'Streaming SHA-256 checks run concurrently; original frozen collector and all inference files remain intact',
            'inference_rerun':False,'truth_read':False}
    (root/'COLLECTION_IMPLEMENTATION.json').write_text(json.dumps(record,indent=2)+'\n')
    runner.verify=lambda hashes:verify_parallel(hashes,workers)
    runner.collect(root)
    print(json.dumps(record),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True,type=Path)
    parser.add_argument('--workers',type=int,default=16)
    args=parser.parse_args()
    collect(args.root,args.workers)
