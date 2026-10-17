"""Verify downloaded archives and their runtime payload before local unblinding."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile

from .evaluation import complete, dump, read, sha


def receive(batch_archive, reproducibility_archive, batch_sha256, reproducibility_sha256, output):
    batch_archive,reproducibility_archive,output=map(Path,(batch_archive,reproducibility_archive,output))
    if sha(batch_archive)!=batch_sha256 or sha(reproducibility_archive)!=reproducibility_sha256:
        raise ValueError('Remote/local archive SHA-256 mismatch')
    if output.exists():
        raise FileExistsError('Use a new verified extraction directory')
    output.mkdir(parents=True)
    batch,repro=output/'batch',output/'reproducibility'
    for source,destination in ((batch_archive,batch),(reproducibility_archive,repro)):
        destination.mkdir()
        with tarfile.open(source) as archive:
            archive.extractall(destination,filter='data')
    selection,calls=complete(batch)
    protocol=read(batch/'PROTOCOL.json')
    manifests=list(repro.glob('REPRODUCIBILITY_*.json'))
    if len(manifests)!=1:
        raise ValueError('Ambiguous archive manifest')
    manifest=read(manifests[0])
    for name,digest in manifest['artifact_sha256'].items():
        if sha(repro/name)!=digest:
            raise ValueError('Reproducibility member changed: '+name)
    inference=repro/'inference'
    if sha(inference/'model/model.joblib')!=protocol['model_sha256'] or sha(inference/'model/TRAINING.json')!=protocol['training_sha256']:
        raise ValueError('Checkpoint/training differs from frozen protocol')
    if (inference/'model/TRAINING.json').read_bytes()!=(batch/'model/TRAINING.json').read_bytes():
        raise ValueError('Training copies disagree')
    remote_root=manifest['inference_root'].rstrip('/')+'/'
    with tarfile.open(inference/'sources.tar.gz') as source:
        members={m.name:hashlib.sha256(source.extractfile(m).read()).hexdigest() for m in source if m.isfile()}
    for remote,digest in protocol['code_sha256'].items():
        if not remote.startswith(remote_root):
            raise ValueError('Code path outside archived inference root')
        relative=remote[len(remote_root):]
        actual=sha(inference/relative) if relative.endswith('.sbatch') else members.get(relative)
        if actual!=digest:
            raise ValueError('Frozen runtime payload member mismatch: '+relative)
    training_hashes=read(repro/'TRAIN_CODE_SHA.json')
    with tarfile.open(repro/'train_code.tar.gz') as archive:
        for name,digest in training_hashes.items():
            if hashlib.sha256(archive.extractfile(name).read()).hexdigest()!=digest:
                raise ValueError('Original model-training source mismatch: '+name)
    record={'remote_archive_sha256':{'batch':batch_sha256,'reproducibility':reproducibility_sha256},
            'samples':len(selection),'plain_prediction_sets_match_combined_JSON':True,
            'frozen_code_sources_and_submission_scripts_verified':True,'frozen_model_and_training_verified':True,
            'original_model_training_code_verified':True,'reproducibility_member_hashes_verified':True,
            'truth_read':False,'inference_root':manifest['inference_root']}
    dump(output/'DOWNLOAD_AUDIT.json',record)
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch-archive',required=True,type=Path)
    parser.add_argument('--reproducibility-archive',required=True,type=Path)
    parser.add_argument('--batch-sha256',required=True)
    parser.add_argument('--reproducibility-sha256',required=True)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    receive(args.batch_archive,args.reproducibility_archive,args.batch_sha256,args.reproducibility_sha256,args.output)
