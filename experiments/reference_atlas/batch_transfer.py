"""Freeze, run and collect a separate 100-sample transfer experiment on TTU."""
import argparse
import hashlib
import json
import platform
import subprocess
import sys
import tarfile
from importlib.metadata import version
from pathlib import Path

from transfer import digest, predict_sample, write_json

PYTHON = "/home/nhu/miniconda3/envs/sprout/bin/python"


def freeze(root, source):
    root, source = Path(root).resolve(), Path(source).resolve()
    if (root / "PROTOCOL.json").exists():
        raise FileExistsError("Protocol already frozen")
    original = json.loads((source / "selection.json").read_text())
    rows = original["selected"]
    assert len(rows) == 100 and (source / "PREDICTIONS_COMPLETE").exists()
    assert all((Path(s["run"]) / "SUCCESS").exists() for s in rows)
    selection = [{k: s[k] for k in ("index", "sample", "project", "run", "previously_seen")} for s in rows]
    write_json(root / "selection.json", {"selected": selection, "unique_compositions": original["unique_compositions"]})
    code = {p.name: digest(p) for p in (root / "transfer.py", root / "batch_transfer.py")}
    write_json(root / "PROTOCOL.json", {"schema": 1, "source_root": str(source),
        "software": {"python": platform.python_version(), **{p: version(p) for p in
                     ("numpy", "scikit-learn", "biopython", "joblib")}},
        "source_selection_sha256": digest(source / "selection.json"),
        "selection_sha256": digest(root / "selection.json"), "code_sha256": code,
        "model_sha256": digest(root / "model/model.joblib"),
        "training_sha256": digest(root / "model/training.json"), "truth_used_for_training": False,
        "parameters": {"primary": "profile_gene3", "independent_gene_votes": 3,
                       "min_valid_genes": 10, "baseline_confirmation_gene_votes": 2,
                       "candidate_count_cap": None, "reference_training_trees_per_gene": 1,
                       "reference_training_queries_per_tree": 64,
                       "calibration": "reference-only species-balanced known acceptance / absent-family rejection utility; 101 fixed cutoffs",
                       "prediction": "nearest family by full distances, verified by learned full-profile features"},
        "reuse": "identical order predictions, assembly, family reference labels, topology and branch lengths",
        "evaluation_status": "retrospective reanalysis of previously evaluated samples; not a new independent test"})
    (root / "logs").mkdir(exist_ok=True)


def verify(root):
    root = Path(root)
    protocol = json.loads((root / "PROTOCOL.json").read_text())
    assert all(digest(root / name) == sha for name, sha in protocol["code_sha256"].items())
    assert digest(root / "model/model.joblib") == protocol["model_sha256"]
    assert digest(root / "model/training.json") == protocol["training_sha256"]
    assert digest(root / "selection.json") == protocol["selection_sha256"]
    assert digest(Path(protocol["source_root"]) / "selection.json") == protocol["source_selection_sha256"]
    return protocol


def submit(root):
    root = Path(root)
    verify(root)
    if (root / "JOBS.json").exists():
        raise FileExistsError("Jobs already submitted")
    wrap = f'{PYTHON} {root}/batch_transfer.py predict --root {root} --index "$SLURM_ARRAY_TASK_ID"'
    job = subprocess.check_output(["sbatch", "--parsable", "--partition=nocona", "--nodes=1",
        "--ntasks=1", "--cpus-per-task=2", "--mem=8G", "--time=04:00:00", "--array=0-99%8",
        "--job-name=atlas_family100", "--output=" + str(root / "logs/%A_%a.out"),
        "--error=" + str(root / "logs/%A_%a.err"), "--wrap=" + wrap], text=True).strip()
    write_json(root / "JOBS.json", {"array": job, "collector": None})
    collect = subprocess.check_output(["sbatch", "--parsable", "--partition=nocona", "--nodes=1", "--ntasks=1", "--cpus-per-task=1",
        "--mem=4G", "--time=00:30:00", "--job-name=atlas_collect", "--dependency=afterany:" + job,
        "--output=" + str(root / "logs/collect_%j.out"), "--error=" + str(root / "logs/collect_%j.err"),
        "--wrap=" + f"{PYTHON} {root}/batch_transfer.py collect --root {root}"], text=True).strip()
    write_json(root / "JOBS.json", {"array": job, "collector": collect})
    print(json.dumps({"array": job, "collector": collect}), flush=True)


def predict(root, index):
    root = Path(root)
    protocol = verify(root)
    row = json.loads((root / "selection.json").read_text())["selected"][index]
    out = root / "samples" / row["project"]
    out.mkdir(parents=True, exist_ok=True)
    try:
        predict_sample(protocol["source_root"], root / "model", root, index)
    except Exception as exc:
        (out / "FAILED").write_text(repr(exc))
        raise


def collect(root):
    root = Path(root)
    protocol = verify(root)
    rows = json.loads((root / "selection.json").read_text())["selected"]
    status = []
    for row in rows:
        p = root / "samples" / row["project"]
        status.append({"sample": row["sample"], "success": (p / "SUCCESS").exists(),
                       "failed": (p / "FAILED").exists()})
    write_json(root / "STATUS.json", status)
    assert all(s["success"] and not s["failed"] for s in status), "Not all 100 samples completed"
    for row in rows:
        audit = json.loads((root / "samples" / row["project"] / "audit.json").read_text())
        assert audit["model_sha256"] == protocol["model_sha256"]
        assert audit["code_sha256"] == protocol["code_sha256"]["transfer.py"]
        assert all(digest(path) == sha for path, sha in audit["input_sha256"].items())
    write_json(root / "AUDIT.json", {"all_100_complete": True, "all_reused_input_hashes_unchanged": True,
                                   "code_model_protocol_unchanged": True, "truth_used": False})
    (root / "PREDICTIONS_COMPLETE").write_text("100/100\n")
    with tarfile.open(root / "compact_results.tar.gz", "w:gz") as tar:
        for name in ("selection.json", "PROTOCOL.json", "STATUS.json", "AUDIT.json", "PREDICTIONS_COMPLETE", "JOBS.json", "model/training.json"):
            tar.add(root / name, arcname=name)
        for row in rows:
            p = root / "samples" / row["project"]
            for file in p.iterdir():
                if file.is_file() and file.name != "profile_features.csv":
                    tar.add(file, arcname=str(file.relative_to(root)))
    print("100/100 complete; immutable-input audit passed", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["freeze", "submit", "predict", "collect", "status"])
    parser.add_argument("--root", required=True)
    parser.add_argument("--source")
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    if args.command == "freeze":
        freeze(args.root, args.source)
    elif args.command == "submit":
        submit(args.root)
    elif args.command == "predict":
        predict(args.root, args.index)
    elif args.command == "collect":
        collect(args.root)
    else:
        root = Path(args.root)
        rows = json.loads((root / "selection.json").read_text())["selected"]
        complete = sum((root / "samples" / r["project"] / "SUCCESS").exists() for r in rows)
        failed = [r["sample"] for r in rows if (root / "samples" / r["project"] / "FAILED").exists()]
        print(json.dumps({"complete": complete, "failed": failed, "collected": (root / "PREDICTIONS_COMPLETE").exists()}))


if __name__ == "__main__":
    main()
