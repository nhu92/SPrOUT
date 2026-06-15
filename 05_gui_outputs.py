#!/usr/bin/env python3
"""
05_gui_outputs.py – Build GUI-ready SPrOUT result bundle.

Collects outputs from prior SPrOUT steps into a compact, documented directory that
can be consumed by an external GUI (for example an R Shiny app). It creates:
- run_metadata.json: basic run and file-count information
- exon_metrics.csv: per-exon sequence length and mapping-coverage summaries
- all_exon_trees.nwk: all Newick exon trees concatenated with tree labels
- tree_inventory.csv: one row per exon tree for GUI tree selection
- tree_contributions.csv: per-tree ACS contribution to final prediction taxa
- result_manifest.json: machine-readable list of generated GUI files
- <project>.sprout_results.zip: portable archive of the GUI-ready directory
"""
import argparse
import json
import re
import shutil
import statistics
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from pipeline_utils import is_valid_project_name, load_config


def _safe_read_csv(path):
    return pd.read_csv(path) if path and Path(path).exists() else None


def parse_tree_label(tree_path):
    stem = Path(tree_path).stem
    match = re.match(r"(?P<gene>.+)_exon_(?P<exon>\d+)$", stem)
    return {
        "tree_id": stem,
        "gene": match.group("gene") if match else stem,
        "exon_index": int(match.group("exon")) if match else None,
        "tree_file": str(tree_path),
    }


def build_tree_inventory(input_phylo, output_file, merged_tree_file):
    rows = []
    tree_paths = sorted(Path(input_phylo).glob("*.tre")) if input_phylo else []
    with open(merged_tree_file, "w", encoding="utf-8") as merged:
        for tree_path in tree_paths:
            row = parse_tree_label(tree_path)
            newick = tree_path.read_text(encoding="utf-8").strip()
            row["newick_length"] = len(newick)
            rows.append(row)
            merged.write(f"[{row['tree_id']}] {newick}\n")
    pd.DataFrame(rows).to_csv(output_file, index=False)
    return rows


def exon_metrics(input_exon, exon_split_dir, output_file):
    rows = []
    exon_files = sorted(Path(input_exon).glob("*.fasta")) if input_exon else []
    for fasta in exon_files:
        records = list(SeqIO.parse(str(fasta), "fasta"))
        lengths = [len(record.seq) for record in records]
        parts = fasta.stem.split("_exon_")
        rows.append({
            "exon_id": fasta.stem,
            "gene": parts[0] if parts else fasta.stem,
            "exon_index": int(parts[1]) if len(parts) == 2 and parts[1].isdigit() else None,
            "sequence_count": len(records),
            "min_length": min(lengths) if lengths else 0,
            "mean_length": round(statistics.mean(lengths), 2) if lengths else 0,
            "max_length": max(lengths) if lengths else 0,
            "total_bases": sum(lengths),
            "fasta_file": str(fasta),
        })

    split_rows = []
    if exon_split_dir:
        for tsv in sorted(Path(exon_split_dir).glob("*_exon_split.tsv")):
            try:
                df = pd.read_csv(tsv, sep="\t")
            except Exception:
                continue
            if "exon_names" in df.columns:
                split_rows.append({"exon_split_file": str(tsv), "alignment_hit_rows": len(df)})

    metrics = pd.DataFrame(rows)
    if not metrics.empty:
        metrics.to_csv(output_file, index=False)
    else:
        pd.DataFrame(columns=[
            "exon_id", "gene", "exon_index", "sequence_count", "min_length",
            "mean_length", "max_length", "total_bases", "fasta_file"
        ]).to_csv(output_file, index=False)
    return rows, split_rows


def infer_tree_id_from_matrix(filename):
    stem = Path(filename).name.replace(".cleaned.csv", "")
    if "." in stem:
        gene, idx = stem.rsplit(".", 1)
        return gene, int(idx) if idx.isdigit() else None, stem
    return stem, None, stem


def build_tree_contributions(matrix_dir, project, final_summary_file, output_file):
    final = _safe_read_csv(final_summary_file)
    selected_taxa = set(final["row_name"].astype(str)) if final is not None and "row_name" in final else set()
    rows = []
    for matrix_path in sorted(Path(matrix_dir).glob("*.cleaned.csv")) if matrix_dir else []:
        df = pd.read_csv(matrix_path)
        if df.empty:
            continue
        label_col = df.columns[0]
        project_cols = [col for col in df.columns[1:] if project in col]
        if not project_cols:
            project_cols = list(df.columns[1:])
        numeric = df[project_cols].apply(pd.to_numeric, errors="coerce")
        similarity = 1 / (1 + numeric)
        taxa = df[label_col].astype(str).str.replace(r"\d+", "", regex=True).str.rstrip("_")
        gene, exon_index, matrix_id = infer_tree_id_from_matrix(matrix_path)
        for taxon in sorted(set(taxa)):
            mask = taxa == taxon
            acs = float(similarity.loc[mask].sum().sum())
            rows.append({
                "tree_id": matrix_id,
                "gene": gene,
                "exon_index": exon_index,
                "taxon": taxon,
                "acs_contribution": round(acs, 6),
                "included_in_final_threshold_result": taxon in selected_taxa,
                "matrix_file": str(matrix_path),
            })
    pd.DataFrame(rows).to_csv(output_file, index=False)
    return rows


def write_metadata(args, output_dir, tree_rows, exon_rows, split_rows, contribution_rows):
    metadata = {
        "project": args.proj_name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": vars(args),
        "counts": {
            "exon_fastas": len(exon_rows),
            "exon_trees": len(tree_rows),
            "exon_split_tables": len(split_rows),
            "tree_contribution_rows": len(contribution_rows),
        },
    }
    path = Path(output_dir) / "run_metadata.json"
    path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def zip_output(output_dir, project):
    archive_base = Path(output_dir).with_name(f"{project}.sprout_results")
    archive = shutil.make_archive(str(archive_base), "zip", root_dir=output_dir)
    return archive


def main():
    parser = argparse.ArgumentParser(description="Create GUI-ready SPrOUT result files and ZIP archive.")
    parser.add_argument("-c", "--config", help="Path to config file (YAML/JSON/TOML)")
    parser.add_argument("-p", "--proj_name", help="Project name identifier")
    parser.add_argument("--input_exon", default=None, help="Directory with extracted exon FASTA files")
    parser.add_argument("--input_phylo", default=None, help="Directory with exon .tre files")
    parser.add_argument("--matrix_dir", default=None, help="Directory with per-tree cleaned matrices")
    parser.add_argument("--prediction_summary", default=None, help="Final prediction summary CSV from 04_prediction.py")
    parser.add_argument("--exon_split_dir", default=None, help="Directory containing *_exon_split.tsv files")
    parser.add_argument("-o", "--output_dir", default=None, help="Directory for GUI-ready outputs")
    args = parser.parse_args()

    config = load_config(args.config) if args.config else {}
    args.proj_name = args.proj_name or config.get("proj_name")
    args.input_exon = args.input_exon or config.get("output_exon", "02_exon_extracted")
    args.input_phylo = args.input_phylo or config.get("output_phylo", "03_phylo_results")
    args.matrix_dir = args.matrix_dir or config.get("output_tree", "04_all_trees")
    args.prediction_summary = args.prediction_summary or config.get("output_file")
    args.exon_split_dir = args.exon_split_dir or args.input_exon
    args.output_dir = args.output_dir or config.get("gui_output_dir", "05_gui_results")

    if not args.proj_name:
        parser.error("proj_name is required.")
    if not is_valid_project_name(args.proj_name):
        parser.error(f"Project name '{args.proj_name}' contains invalid characters.")

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    exon_rows, split_rows = exon_metrics(args.input_exon, args.exon_split_dir, out / "exon_metrics.csv")
    tree_rows = build_tree_inventory(args.input_phylo, out / "tree_inventory.csv", out / "all_exon_trees.nwk")
    contribution_rows = build_tree_contributions(args.matrix_dir, args.proj_name, args.prediction_summary, out / "tree_contributions.csv")
    metadata = write_metadata(args, out, tree_rows, exon_rows, split_rows, contribution_rows)
    manifest = {
        "metadata": "run_metadata.json",
        "exon_metrics": "exon_metrics.csv",
        "tree_inventory": "tree_inventory.csv",
        "merged_exon_trees": "all_exon_trees.nwk",
        "tree_contributions": "tree_contributions.csv",
        "notes": "Use tree_id/gene/exon_index as stable keys for GUI tree views and slider-driven prediction summaries.",
    }
    (out / "result_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    archive = zip_output(out, args.proj_name)
    print(json.dumps({"output_dir": str(out), "archive": archive, "counts": metadata["counts"]}, indent=2))


if __name__ == "__main__":
    main()
