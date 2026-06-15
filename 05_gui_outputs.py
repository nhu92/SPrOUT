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
<<<<<<< ours
=======
import ast
import importlib.util
>>>>>>> theirs
import json
import re
import shutil
import statistics
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from pipeline_utils import is_valid_project_name, load_config

<<<<<<< ours

def _safe_read_csv(path):
    return pd.read_csv(path) if path and Path(path).exists() else None


def parse_tree_label(tree_path):
    stem = Path(tree_path).stem
    match = re.match(r"(?P<gene>.+)_exon_(?P<exon>\d+)$", stem)
=======
EMPTY_EXON_COLUMNS = [
    "exon_id", "gene", "exon_index", "sequence_count", "min_length",
    "mean_length", "max_length", "total_bases", "mapped_bases",
    "target_bases", "mapping_coverage", "alignment_hit_rows", "fasta_file",
]
EMPTY_TREE_COLUMNS = ["tree_id", "gene", "exon_index", "tree_file", "newick_length"]
EMPTY_CONTRIBUTION_COLUMNS = [
    "tree_id", "gene", "exon_index", "taxon", "acs_contribution",
    "fraction_of_tree_acs", "included_in_final_threshold_result", "matrix_file",
]


def load_distance_helpers():
    """Load helpers from 03_distance_matrices.py, whose filename is not importable normally."""
    module_path = Path(__file__).with_name("03_distance_matrices.py")
    spec = importlib.util.spec_from_file_location("sprout_distance_matrices", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.clean_up_matrix, module.distance_to_similarity


CLEAN_UP_MATRIX, DISTANCE_TO_SIMILARITY = load_distance_helpers()


def safe_read_csv(path):
    return pd.read_csv(path) if path and Path(path).exists() else None


def safe_literal_list(value):
    if pd.isna(value):
        return []
    try:
        parsed = ast.literal_eval(str(value))
    except (ValueError, SyntaxError):
        return []
    return parsed if isinstance(parsed, list) else []


def interval_total(intervals):
    total = 0
    for item in intervals:
        if isinstance(item, (list, tuple)) and len(item) >= 2:
            try:
                start, end = int(item[0]), int(item[1])
            except (TypeError, ValueError):
                continue
            total += max(0, end - start)
    return total


def parse_tree_label(tree_path):
    stem = Path(tree_path).stem
    underscore_match = re.match(r"(?P<gene>.+)_exon_(?P<exon>\d+)$", stem)
    dotted_match = re.match(r"(?P<gene>.+)\.(?P<exon>\d+)$", stem)
    match = underscore_match or dotted_match
>>>>>>> theirs
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
<<<<<<< ours
    pd.DataFrame(rows).to_csv(output_file, index=False)
    return rows


def exon_metrics(input_exon, exon_split_dir, output_file):
    rows = []
=======
    pd.DataFrame(rows, columns=EMPTY_TREE_COLUMNS).to_csv(output_file, index=False)
    return rows


def summarize_exon_split_tables(exon_split_dir):
    summaries = {}
    split_rows = []
    if not exon_split_dir:
        return summaries, split_rows

    for tsv in sorted(Path(exon_split_dir).glob("*_exon_split.tsv")):
        try:
            df = pd.read_csv(tsv, sep="\t")
        except Exception:
            continue
        split_rows.append({"exon_split_file": str(tsv), "alignment_hit_rows": len(df)})
        if "exon_names" not in df.columns:
            continue
        for _, row in df.iterrows():
            exon_names = safe_literal_list(row.get("exon_names"))
            mapped_ranges = safe_literal_list(row.iloc[13]) if len(row) > 13 else []
            target_ranges = safe_literal_list(row.iloc[6]) if len(row) > 6 else []
            for idx, exon_id in enumerate(exon_names):
                mapped_bases = interval_total([mapped_ranges[idx]]) if idx < len(mapped_ranges) else 0
                target_bases = interval_total([target_ranges[idx]]) if idx < len(target_ranges) else mapped_bases
                current = summaries.setdefault(exon_id, {
                    "mapped_bases": 0,
                    "target_bases": 0,
                    "alignment_hit_rows": 0,
                })
                current["mapped_bases"] += mapped_bases
                current["target_bases"] += target_bases
                current["alignment_hit_rows"] += 1
    return summaries, split_rows


def exon_metrics(input_exon, exon_split_dir, output_file):
    rows = []
    split_summaries, split_rows = summarize_exon_split_tables(exon_split_dir)
>>>>>>> theirs
    exon_files = sorted(Path(input_exon).glob("*.fasta")) if input_exon else []
    for fasta in exon_files:
        records = list(SeqIO.parse(str(fasta), "fasta"))
        lengths = [len(record.seq) for record in records]
        parts = fasta.stem.split("_exon_")
<<<<<<< ours
=======
        split_info = split_summaries.get(fasta.stem, {})
        mapped_bases = split_info.get("mapped_bases", sum(lengths))
        target_bases = split_info.get("target_bases", sum(lengths))
        coverage = round(mapped_bases / target_bases, 6) if target_bases else None
>>>>>>> theirs
        rows.append({
            "exon_id": fasta.stem,
            "gene": parts[0] if parts else fasta.stem,
            "exon_index": int(parts[1]) if len(parts) == 2 and parts[1].isdigit() else None,
            "sequence_count": len(records),
            "min_length": min(lengths) if lengths else 0,
            "mean_length": round(statistics.mean(lengths), 2) if lengths else 0,
            "max_length": max(lengths) if lengths else 0,
            "total_bases": sum(lengths),
<<<<<<< ours
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
=======
            "mapped_bases": mapped_bases,
            "target_bases": target_bases,
            "mapping_coverage": coverage,
            "alignment_hit_rows": split_info.get("alignment_hit_rows", len(records)),
            "fasta_file": str(fasta),
        })

    pd.DataFrame(rows, columns=EMPTY_EXON_COLUMNS).to_csv(output_file, index=False)
>>>>>>> theirs
    return rows, split_rows


def infer_tree_id_from_matrix(filename):
    stem = Path(filename).name.replace(".cleaned.csv", "")
    if "." in stem:
        gene, idx = stem.rsplit(".", 1)
        return gene, int(idx) if idx.isdigit() else None, stem
    return stem, None, stem


<<<<<<< ours
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
=======
def load_selected_taxa(taxonomy_output_file, prediction_summary_file, zscore_threshold):
    if taxonomy_output_file and Path(taxonomy_output_file).exists():
        with open(taxonomy_output_file, encoding="utf-8") as handle:
            return {line.strip() for line in handle if line.strip()}
    final = safe_read_csv(prediction_summary_file)
    if final is None or "row_name" not in final:
        return set()
    if zscore_threshold is not None and "z_score" in final:
        return set(final.loc[final["z_score"] > zscore_threshold, "row_name"].astype(str))
    return set(final["row_name"].astype(str))


def build_tree_contributions(
    matrix_dir, project, prediction_summary_file, taxonomy_output_file,
    zscore_threshold, threshold, use_flag, use_threshold, output_file,
):
    selected_taxa = load_selected_taxa(taxonomy_output_file, prediction_summary_file, zscore_threshold)
    rows = []
    for matrix_path in sorted(Path(matrix_dir).glob("*.cleaned.csv")) if matrix_dir else []:
        raw_df = pd.read_csv(matrix_path)
        if raw_df.empty:
            continue
        prefix = matrix_path.name.split("cleaned.csv")[0]
        taxa_file = matrix_path.with_name(f"{prefix}list.txt")
        cleaned = CLEAN_UP_MATRIX(
            raw_df, project, threshold,
            str(taxa_file) if taxa_file.exists() else None,
            use_flag, use_threshold,
        )
        if cleaned.empty:
            continue
        similarity = DISTANCE_TO_SIMILARITY(cleaned)
        label_col = similarity.columns[0]
        value_cols = [col for col in similarity.columns[1:] if project in col]
        if not value_cols:
            value_cols = list(similarity.columns[1:])
        gene, exon_index, matrix_id = infer_tree_id_from_matrix(matrix_path)
        tree_total = float(similarity[value_cols].sum().sum()) if value_cols else 0.0
        for taxon in sorted(set(similarity[label_col].astype(str))):
            mask = similarity[label_col].astype(str) == taxon
            acs = float(similarity.loc[mask, value_cols].sum().sum()) if value_cols else 0.0
>>>>>>> theirs
            rows.append({
                "tree_id": matrix_id,
                "gene": gene,
                "exon_index": exon_index,
                "taxon": taxon,
                "acs_contribution": round(acs, 6),
<<<<<<< ours
                "included_in_final_threshold_result": taxon in selected_taxa,
                "matrix_file": str(matrix_path),
            })
    pd.DataFrame(rows).to_csv(output_file, index=False)
=======
                "fraction_of_tree_acs": round(acs / tree_total, 6) if tree_total else 0.0,
                "included_in_final_threshold_result": taxon in selected_taxa,
                "matrix_file": str(matrix_path),
            })
    pd.DataFrame(rows, columns=EMPTY_CONTRIBUTION_COLUMNS).to_csv(output_file, index=False)
>>>>>>> theirs
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
<<<<<<< ours
=======
    parser.add_argument("--taxonomy_output_file", default=None, help="Selected taxa text file from 04_prediction.py")
    parser.add_argument("--zscore_threshold", type=float, default=None, help="Threshold used by 04_prediction.py")
    parser.add_argument("--threshold", type=float, default=None, help="Distance filtering threshold used by 03_distance_matrices.py")
    parser.add_argument("--use_flag", action="store_true", help="Use flag method when recomputing per-tree ACS contributions")
    parser.add_argument("--use_threshold", action="store_true", help="Use threshold filtering when recomputing per-tree ACS contributions")
>>>>>>> theirs
    parser.add_argument("--exon_split_dir", default=None, help="Directory containing *_exon_split.tsv files")
    parser.add_argument("-o", "--output_dir", default=None, help="Directory for GUI-ready outputs")
    args = parser.parse_args()

    config = load_config(args.config) if args.config else {}
    args.proj_name = args.proj_name or config.get("proj_name")
    args.input_exon = args.input_exon or config.get("output_exon", "02_exon_extracted")
    args.input_phylo = args.input_phylo or config.get("output_phylo", "03_phylo_results")
    args.matrix_dir = args.matrix_dir or config.get("output_tree", "04_all_trees")
    args.prediction_summary = args.prediction_summary or config.get("output_file")
<<<<<<< ours
=======
    args.taxonomy_output_file = args.taxonomy_output_file or config.get("taxonomy_output_file")
    args.zscore_threshold = args.zscore_threshold if args.zscore_threshold is not None else config.get("zscore_threshold")
    args.threshold = args.threshold if args.threshold is not None else config.get("threshold", 1.96)
    args.use_flag = args.use_flag or bool(config.get("use_flag", False))
    args.use_threshold = args.use_threshold or bool(config.get("use_threshold", False))
>>>>>>> theirs
    args.exon_split_dir = args.exon_split_dir or args.input_exon
    args.output_dir = args.output_dir or config.get("gui_output_dir", "05_gui_results")

    if not args.proj_name:
        parser.error("proj_name is required.")
    if not is_valid_project_name(args.proj_name):
        parser.error(f"Project name '{args.proj_name}' contains invalid characters.")
<<<<<<< ours
=======
    if args.use_flag and args.use_threshold:
        parser.error("You cannot enable both --use_flag and --use_threshold at the same time.")

    args.threshold = float(args.threshold)
    args.zscore_threshold = float(args.zscore_threshold) if args.zscore_threshold is not None else None
>>>>>>> theirs

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    exon_rows, split_rows = exon_metrics(args.input_exon, args.exon_split_dir, out / "exon_metrics.csv")
    tree_rows = build_tree_inventory(args.input_phylo, out / "tree_inventory.csv", out / "all_exon_trees.nwk")
<<<<<<< ours
    contribution_rows = build_tree_contributions(args.matrix_dir, args.proj_name, args.prediction_summary, out / "tree_contributions.csv")
=======
    contribution_rows = build_tree_contributions(
        args.matrix_dir, args.proj_name, args.prediction_summary, args.taxonomy_output_file,
        args.zscore_threshold, args.threshold, args.use_flag, args.use_threshold,
        out / "tree_contributions.csv",
    )
>>>>>>> theirs
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
