#!/usr/bin/env python3
"""Export small, readable SPrOUT reports for SPrOUT-PEA.

Only rankings and the nearest reference per query exon are exported. FASTQ,
alignments, complete trees, and distance matrices never enter the web report.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import median

SCHEMA = "sprout-pea/v1"
LEVEL_INDEX = {"order": 0, "family": 1, "genus": 2, "species": 3}
STAGE_ORDER = {"order": 0, "family": 1, "branch": 2, "genus": 3, "species": 4}
MAX_REPORT_BYTES = 20 * 1024 * 1024


def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def taxon_from_leaf(leaf: str, level: str) -> str | None:
    parts = leaf.strip().strip("'\"").removeprefix("_R_").split("_")
    if len(parts) < 4 or parts[0] == "NODE":
        return None
    if level == "species":
        return " ".join(parts[2:4])
    return parts[LEVEL_INDEX[level]] or None


def read_names(path: Path | None) -> set[str]:
    if path is None:
        return set()
    return {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


def read_predictions(path: Path, evidence: Path | None, candidates: Path | None) -> list[dict]:
    extra = {}
    if evidence is not None:
        with evidence.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                extra[row["taxon"]] = row
    chosen = read_names(candidates)
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for rank, row in enumerate(csv.DictReader(handle), 1):
            name = (row.get("row_name") or "").strip()
            if not name:
                continue
            data = extra.get(name, {})
            rows.append({
                "taxon": name,
                "rank": rank,
                "score": number(row.get("sum_of_total_value")),
                "z_score": number(row.get("z_score")),
                "pipeline_selected": name in chosen,
                "selected_by": [x for x in data.get("selected_by", "").split(";") if x],
                "locus_votes": int(data.get("locus_votes") or 0),
                "query_votes": int(data.get("query_votes") or 0),
                "median_winning_distance": number(data.get("median_winning_distance")),
            })
    if not rows:
        raise ValueError(f"No predictions in {path}")
    return rows


def gene_evidence(matrix_dir: Path | None, project: str, level: str) -> tuple[list[dict], dict]:
    if matrix_dir is None:
        return [], {"matrices": 0, "query_exons": 0, "reference_taxa": 0}
    files = sorted(matrix_dir.glob("*.cleaned.csv"))
    if not files:
        raise ValueError(f"No *.cleaned.csv matrices in {matrix_dir}")
    per_gene = defaultdict(lambda: {"query_exons": 0, "trees": 0, "taxa": defaultdict(list)})
    reference_taxa = set()
    total_queries = 0
    for file in files:
        gene = file.name.split(".", 1)[0]
        with file.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            if len(fields) < 2:
                continue
            query_columns = [name for name in fields[1:] if project in name]
            records = []
            for row in reader:
                leaf = (row.get(fields[0]) or "").strip()
                if not leaf or project in leaf or "NODE" in leaf:
                    continue
                taxon = taxon_from_leaf(leaf, level)
                if taxon:
                    reference_taxa.add(taxon)
                    records.append((leaf, taxon, row))
        bucket = per_gene[gene]
        bucket["trees"] += 1
        for query in query_columns:
            nearest = []
            for leaf, taxon, row in records:
                distance = number(row.get(query))
                if distance is not None and distance >= 0:
                    nearest.append((distance, taxon, leaf))
            if not nearest:
                continue
            best = min(item[0] for item in nearest)
            winners = [item for item in nearest if math.isclose(item[0], best, rel_tol=1e-9, abs_tol=1e-12)]
            bucket["query_exons"] += 1
            total_queries += 1
            # One vote per taxon and query exon, even when several references tie.
            for taxon in sorted({item[1] for item in winners}):
                leaf = next(item[2] for item in winners if item[1] == taxon)
                bucket["taxa"][taxon].append((best, leaf))
    genes = []
    for gene, data in sorted(per_gene.items()):
        matches = []
        for taxon, hits in data["taxa"].items():
            matches.append({
                "taxon": taxon,
                "query_votes": len(hits),
                "median_distance": median(distance for distance, _ in hits),
                "example_reference": hits[0][1],
            })
        matches.sort(key=lambda item: (-item["query_votes"], item["median_distance"], item["taxon"]))
        genes.append({"gene": gene, "trees": data["trees"], "query_exons": data["query_exons"], "matches": matches})
    return genes, {"matrices": len(files), "query_exons": total_queries, "reference_taxa": len(reference_taxa)}


def build_stage(args) -> dict:
    genes, counts = gene_evidence(args.matrix_dir, args.project, args.level)
    gene_list = read_names(args.gene_list)
    return {
        "id": args.stage_id or args.level,
        "level": args.level,
        "default_z": args.z_threshold,
        "predictions": read_predictions(args.predictions, args.evidence, args.candidates),
        "genes": genes,
        "counts": {**counts, "requested_genes": len(gene_list)},
        "evidence_available": bool(genes),
        "note": (
            "Per-gene nearest-reference evidence is available. Distances are tree distances, not read proportions."
            if genes else
            "Per-gene matrices were not supplied; only aggregate rankings are available."
        ),
    }


def validate_report(report: dict) -> None:
    if report.get("schema_version") != SCHEMA:
        raise ValueError("Unsupported SPrOUT-PEA schema")
    if not isinstance(report.get("sample_id"), str) or not report["sample_id"]:
        raise ValueError("Missing sample_id")
    if not isinstance(report.get("stages"), list) or not report["stages"]:
        raise ValueError("No stages in report")
    ids = [stage["id"] for stage in report["stages"]]
    if len(set(ids)) != len(ids):
        raise ValueError(f"Duplicate stage IDs: {ids}")


def markdown(report: dict) -> str:
    lines = [f"# SPrOUT results: {report['sample_id']}", "", "This is a research result, not a confirmed species identification.", ""]
    for stage in report["stages"]:
        lines.extend([
            f"## {stage['id'].title()} ({stage['level']})", "",
            f"Default Z threshold: {stage['default_z']:g}; "
            f"genes with evidence: {len(stage['genes'])}; "
            f"query exons: {stage['counts']['query_exons']}.", "",
            "| Rank | Taxon | Score | Z | Locus votes | Pipeline selected |",
            "| ---: | --- | ---: | ---: | ---: | :---: |",
        ])
        for item in stage["predictions"][:30]:
            score = "—" if item["score"] is None else f"{item['score']:.3f}"
            z = "—" if item["z_score"] is None else f"{item['z_score']:.3f}"
            lines.append(f"| {item['rank']} | {item['taxon']} | {score} | {z} | {item['locus_votes']} | {'yes' if item['pipeline_selected'] else 'no'} |")
        lines.extend(["", stage["note"], ""])
        if stage["genes"]:
            lines.extend([
                "### Gene-by-gene nearest references", "",
                "Each vote is a query exon whose closest reference is in that taxon; ties may vote for multiple taxa.", "",
                "| Gene | Trees | Query exons | Leading taxon | Votes | Median tree distance |",
                "| --- | ---: | ---: | --- | ---: | ---: |",
            ])
            for gene in stage["genes"]:
                leader = gene["matches"][0] if gene["matches"] else None
                taxon = leader["taxon"] if leader else "—"
                votes = leader["query_votes"] if leader else "—"
                distance = f"{leader['median_distance']:.4f}" if leader else "—"
                lines.append(
                    f"| {gene['gene']} | {gene['trees']} | {gene['query_exons']} | "
                    f"{taxon} | {votes} | {distance} |"
                )
            lines.append("")
    lines.extend([
        "## Interpretation", "",
        "Changing a Z threshold in SPrOUT-PEA only explores the existing ranking. "
        "It does not add references or recompute the trees. "
        "Low-abundance taxa may be missed; review gene evidence and reference coverage before rerunning.", "",
    ])
    return "\n".join(lines)


def write_report(report: dict, output: Path, markdown_path: Path | None) -> None:
    validate_report(report)
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(report, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")
    if len(payload) > MAX_REPORT_BYTES:
        raise ValueError(f"Report is {len(payload):,} bytes; limit is {MAX_REPORT_BYTES:,}. Reduce per-gene evidence.")
    if output.name.endswith(".gz"):
        with gzip.open(output, "wb", compresslevel=9) as handle:
            handle.write(payload)
    else:
        output.write_bytes(payload)
    if markdown_path:
        markdown_path.parent.mkdir(parents=True, exist_ok=True)
        markdown_path.write_text(markdown(report), encoding="utf-8")
    print(f"SPrOUT-PEA report: {output} ({output.stat().st_size:,} bytes)")
    if markdown_path:
        print(f"Readable summary: {markdown_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    stage = sub.add_parser("stage", help="Export one SPrOUT order/family/genus stage")
    stage.add_argument("--project", required=True)
    stage.add_argument("--level", required=True, choices=LEVEL_INDEX)
    stage.add_argument("--stage-id")
    stage.add_argument("--predictions", required=True, type=Path)
    stage.add_argument("--matrix-dir", type=Path)
    stage.add_argument("--evidence", type=Path, help="MrKnowItAll candidate_evidence.csv")
    stage.add_argument("--candidates", type=Path)
    stage.add_argument("--gene-list", type=Path)
    stage.add_argument("--z-threshold", type=float, default=0.0)
    stage.add_argument("--output", required=True, type=Path)
    stage.add_argument("--markdown", type=Path)
    merge = sub.add_parser("merge", help="Combine stage reports into one upload file")
    merge.add_argument("--project", required=True)
    merge.add_argument("--input", required=True, type=Path, nargs="+")
    merge.add_argument("--output", required=True, type=Path)
    merge.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.project):
        parser.error("Project must contain only letters, numbers, dots, underscores, or hyphens")
    if args.command == "stage":
        stages = [build_stage(args)]
    else:
        stages = []
        for path in args.input:
            raw = gzip.open(path, "rt", encoding="utf-8") if path.name.endswith(".gz") else path.open(encoding="utf-8")
            with raw as handle:
                item = json.load(handle)
            validate_report(item)
            if item["sample_id"] != args.project:
                parser.error(f"Sample mismatch in {path}: {item['sample_id']}")
            stages.extend(item["stages"])
        stages.sort(key=lambda item: (STAGE_ORDER.get(item["id"], 99), item["id"]))
    report = {
        "schema_version": SCHEMA,
        "sample_id": args.project,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "producer": "SPrOUT / MrKnowItAll",
        "stages": stages,
    }
    write_report(report, args.output, args.markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
