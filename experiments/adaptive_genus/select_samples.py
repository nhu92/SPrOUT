"""A reproducible mixture-size sample, independent of taxon truth or accuracy."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path


def select(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[len(row['sample'].split('x'))].append(row)
    chosen,compositions = [],set()
    for size,quota in ((3,4),(6,3),(10,3)):
        ordered = sorted(grouped[size],key=lambda r:hashlib.sha256(('adaptive-genus-20261017/'+r['sample']).encode()).hexdigest())
        count = 0
        for row in ordered:
            identity = tuple(sorted(row['sample'].split('x')))
            if identity in compositions:
                continue
            chosen.append(row);compositions.add(identity);count+=1
            if count==quota:
                break
        if count!=quota:
            raise ValueError('Not enough distinct mixtures in stratum '+str(size))
    return chosen


def main():
    root = Path(__file__).resolve().parents[2]
    source = root/'experiments/reference_atlas/results/transfer100_20261016/selection.json'
    fusion = root/'experiments/reference_atlas/results/fusion100_20261016/predictions.json'
    native = root/'experiments/reference_atlas/runs/transfer100_20261016/batch/samples'
    rows = select(json.loads(source.read_text())['selected'])
    merged = {r['sample']:r for r in json.loads(fusion.read_text())}
    result,hashes = [],{str(source.relative_to(root)):hashlib.sha256(source.read_bytes()).hexdigest(),
                       str(fusion.relative_to(root)):hashlib.sha256(fusion.read_bytes()).hexdigest()}
    for i,row in enumerate(rows):
        zfile = native/row['project']/'baseline_z09.txt'
        z = set(zfile.read_text().split())
        hashes[str(zfile.relative_to(root))] = hashlib.sha256(zfile.read_bytes()).hexdigest()
        prior = merged[row['sample']]
        families = sorted(z|set(prior['merged'])|set(prior['review']))
        result.append({'index':i,'sample':row['sample'],'project':row['project'],
                       'source_run':row['run'],'families':families,
                       'family_search_routes':{'original_z':sorted(z),'merged':prior['merged'],'review':prior['review']}})
    output = Path(__file__).with_name('cohort.json')
    output.write_text(json.dumps({'selection_rule':'SHA256 adaptive-genus-20261017/<sample>; 4/3/3 at 3/6/10 components; unique compositions',
                                 'truth_used':False,'source_sha256':hashes,'selected':result},indent=2)+'\n')
    print(json.dumps([{'sample':r['sample'],'family_search_pool':len(r['families'])} for r in result],indent=2))


if __name__=='__main__':
    main()
