"""Export completed, audited real-mixture predictions without S7 access."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
from zipfile import ZipFile, ZIP_DEFLATED

from .real_truth import verify_complete_batch


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export(batch, output):
    batch, output = Path(batch), Path(output)
    calls = verify_complete_batch(batch)
    if output.exists():
        raise ValueError("Use a new export directory")
    protocol = json.loads((batch / "PROTOCOL.json").read_text())
    policy = Path(__file__).with_name("fusion_policy.json")
    expected = next(v for k, v in protocol["code_and_model_sha256"].items() if k.endswith("/code/fusion_policy.json"))
    if sha(policy) != expected:
        raise ValueError("Local policy differs from the frozen prediction policy")
    with (batch / "evidence.csv").open(encoding="utf-8-sig", newline="") as f:
        evidence = list(csv.DictReader(f))
    lookup = {(r["sample"], r["family"]): r for r in evidence}
    if len(lookup) != len(evidence):
        raise ValueError("Duplicate sample/family evidence")
    flat = []
    for sample in calls:
        for family in sample["merged"]:
            r = lookup[sample["sample"], family]
            if r["merged"] != "True":
                raise ValueError("Evidence and final candidates disagree")
            flat.append({"sample": sample["sample"], "family": family, "parent_order": r["parent_order"],
                "route": r["route"], "ml_gene_votes": r["ml_gene_votes"],
                "stability_lower": r["stability_lower"], "discriminating_genes": r["discriminating_genes"],
                "comparable_genes": r["comparable_genes"]})
    if sum(r["merged"] == "True" for r in evidence) != len(flat):
        raise ValueError("Final candidate/evidence cardinalities differ")
    output.mkdir(parents=True)
    with (output / "merged_predictions.csv").open("w", encoding="utf-8-sig", newline="") as f:
        fields = ["sample", "family", "parent_order", "route", "ml_gene_votes", "stability_lower", "discriminating_genes", "comparable_genes"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flat)
    for source, name in ((batch / "predictions.json", "predictions.json"),
                         (batch / "evidence.csv", "candidate_evidence.csv"), (policy, "POLICY.json")):
        shutil.copy2(source, output / name)
    with ZipFile(output / "prediction_sets.zip", "w", compression=ZIP_DEFLATED) as archive:
        archive.write(policy, "POLICY.json")
        for sample in calls:
            for method in ("merged", "core", "auxiliary", "review"):
                path = batch / "samples" / sample["project"] / f"{method}_families.txt"
                if path.read_text().splitlines() != sample[method]:
                    raise ValueError("Exported source list mismatch")
                archive.write(path, f'{sample["sample"]}/{method}_families.txt')
    record = {"samples": 11, "rank": "family", "primary": "merged", "truth_read_by_exporter": False,
        "sealed_utc": datetime.now(timezone.utc).isoformat(),
        "prediction_audit_sha256": sha(batch / "AUDIT.json"), "frozen_protocol_sha256": sha(batch / "PROTOCOL.json"),
        "artifact_sha256": {p.name: sha(p) for p in output.iterdir() if p.is_file()}}
    (output / "PREDICTION_EXPORT.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"samples": 11, "merged_sample_family_events": len(flat), "truth_read": False}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    export(args.batch, args.output)
