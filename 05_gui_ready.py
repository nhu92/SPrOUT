#!/usr/bin/env python3
"""
05_gui_ready.py – Post-processing utility to package SPrOUT outputs for GUI consumption.

This script does not build a GUI. Instead, it organizes run artifacts into a
standard bundle, emits summary JSON/CSV files for quick loading, and produces a
single ZIP archive that can be handed to a downstream interface.

Key features:
- Merge all exon trees into one multi-tree Newick file.
- Generate a tree index table with tip counts and branch-length statistics.
- Consolidate exon length/coverage metrics generated during assembly.
- Save per-tree contribution tables from distance calculations (if present).
- Build a run summary JSON and manifest describing included files.
- Package everything into a zip (sprout_<project>_results.zip).
"""
import argparse
import glob
import json
import os
import shutil
import zipfile
from datetime import datetime
from statistics import median

import pandas as pd
from Bio import Phylo

from pipeline_utils import is_valid_project_name, load_config


def merge_tree_files(tree_dir, merged_path):
    """Concatenate all .tre files into a single multi-tree Newick file."""
    tree_files = sorted(glob.glob(os.path.join(tree_dir, "*.tre")))
    if not tree_files:
        return 0
    with open(merged_path, "w") as outfile:
        for path in tree_files:
            with open(path) as infile:
                contents = infile.read().strip()
                if contents:
                    outfile.write(contents.rstrip(";"))
                    outfile.write(";\n")
    return len(tree_files)


def build_tree_index(tree_dir, output_path):
    """Create a CSV summarizing each tree for GUI browsing."""
    rows = []
    for path in sorted(glob.glob(os.path.join(tree_dir, "*.tre"))):
        try:
            tree = Phylo.read(path, "newick")
        except Exception:
            continue
        tips = tree.get_terminals()
        branches = [clade.branch_length for clade in tree.get_nonterminals() if clade.branch_length]
        rows.append(
            {
                "tree_file": os.path.basename(path),
                "gene": os.path.basename(path).split("_exon_")[0],
                "tip_count": len(tips),
                "branch_min": min(branches) if branches else None,
                "branch_median": median(branches) if branches else None,
                "branch_max": max(branches) if branches else None,
            }
        )
    if rows:
        pd.DataFrame(rows).to_csv(output_path, index=False)
    return len(rows)


def collect_exon_metrics(exon_dir, project, output_path):
    metrics_files = sorted(glob.glob(os.path.join(exon_dir, f"{project}_*_exon_metrics.tsv")))
    if not metrics_files:
        return 0
    frames = [pd.read_csv(path, sep="\t") for path in metrics_files]
    pd.concat(frames, ignore_index=True).to_csv(output_path, index=False)
    return sum(len(frame) for frame in frames)


def copy_if_exists(paths, destination_dir):
    os.makedirs(destination_dir, exist_ok=True)
    for path in paths:
        for match in glob.glob(path):
            if os.path.isfile(match):
                shutil.copy2(match, destination_dir)


def summarize_run(project, gene_list, tree_count, matrix_dir, exon_metric_rows, threshold, use_flag, use_threshold, output_path):
    summary = {
        "project": project,
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "gene_count": len(gene_list),
        "tree_count": tree_count,
        "matrix_files": len(glob.glob(os.path.join(matrix_dir, "*cleaned.csv"))),
        "exon_metric_rows": exon_metric_rows,
        "threshold": threshold,
        "use_flag": use_flag,
        "use_threshold": use_threshold,
    }
    with open(output_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


def write_manifest(staging_dir, manifest_path):
    entries = []
    for root, _, files in os.walk(staging_dir):
        for name in files:
            full_path = os.path.join(root, name)
            rel_path = os.path.relpath(full_path, staging_dir)
            entries.append(rel_path)
    entries.sort()
    with open(manifest_path, "w") as fh:
        for rel in entries:
            fh.write(rel + "\n")
    return entries


def zip_staging(staging_dir, zip_path):
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for root, _, files in os.walk(staging_dir):
            for file in files:
                full_path = os.path.join(root, file)
                arcname = os.path.relpath(full_path, staging_dir)
                zf.write(full_path, arcname)


def load_gene_names(gene_list_path):
    if not gene_list_path or not os.path.exists(gene_list_path):
        return []
    with open(gene_list_path) as handle:
        return [line.strip().replace('.fasta', '') for line in handle if line.strip()]


def main():
    parser = argparse.ArgumentParser(description="Package SPrOUT outputs for GUI usage (no GUI built).")
    parser.add_argument("-c", "--config", help="Path to config file (YAML/JSON/TOML)")
    parser.add_argument("-p", "--proj_name", help="Project name identifier")
    parser.add_argument("-g", "--gene_list", help="Path to gene list file", default="gene_list.txt")
    parser.add_argument("-e", "--exon_dir", help="Directory with exon FASTAs", default="02_exon_extracted")
    parser.add_argument("--tree_dir", help="Directory with exon trees", default="03_phylo_results")
    parser.add_argument("--matrix_dir", help="Directory with distance matrices", default="04_all_trees")
    parser.add_argument("--threshold", type=float, help="Threshold used for distance filtering")
    parser.add_argument("--use_flag", action="store_true", help="Flag method used in distance filtering")
    parser.add_argument("--use_threshold", action="store_true", help="Threshold filtering enabled")
    parser.add_argument("--package_dir", help="Directory for GUI-ready bundle", default="gui_ready")
    parser.add_argument("--zip_name", help="Custom zip file name")
    args = parser.parse_args()

    config = {}
    if args.config:
        config = load_config(args.config)
    proj_name = args.proj_name or config.get('proj_name')
    gene_list_path = args.gene_list if args.gene_list != parser.get_default('gene_list') else config.get('gene_list', "gene_list.txt")
    exon_dir = args.exon_dir if args.exon_dir != parser.get_default('exon_dir') else config.get('exon_dir', "02_exon_extracted")
    tree_dir = args.tree_dir if args.tree_dir != parser.get_default('tree_dir') else config.get('tree_dir', "03_phylo_results")
    matrix_dir = args.matrix_dir if args.matrix_dir != parser.get_default('matrix_dir') else config.get('matrix_dir', "04_all_trees")
    package_dir = args.package_dir if args.package_dir != parser.get_default('package_dir') else config.get('package_dir', "gui_ready")
    zip_name = args.zip_name or config.get('zip_name')
    threshold = args.threshold if args.threshold is not None else config.get('threshold')
    use_flag = args.use_flag or bool(config.get('use_flag', False))
    use_threshold = args.use_threshold or bool(config.get('use_threshold', False))

    gene_names = load_gene_names(gene_list_path)
    if not proj_name:
        parser.error("Project name is required for packaging.")
    if not is_valid_project_name(proj_name):
        parser.error(f"Project name '{proj_name}' contains invalid characters.")
    staging_dir = os.path.join(package_dir, f"{proj_name}_gui_bundle")
    os.makedirs(staging_dir, exist_ok=True)

    # Merge trees and build index
    trees_dir = os.path.join(staging_dir, "trees")
    os.makedirs(trees_dir, exist_ok=True)
    merged_tree_path = os.path.join(trees_dir, "merged_exon_trees.nwk")
    tree_count = merge_tree_files(tree_dir, merged_tree_path)
    tree_index_path = os.path.join(trees_dir, "tree_index.csv")
    build_tree_index(tree_dir, tree_index_path)

    # Consolidate exon metrics
    metrics_dir = os.path.join(staging_dir, "exon_metrics")
    os.makedirs(metrics_dir, exist_ok=True)
    metrics_path = os.path.join(metrics_dir, "exon_metrics.csv")
    exon_metric_rows = collect_exon_metrics(exon_dir, proj_name, metrics_path)

    # Copy predictions and matrices
    matrices_out = os.path.join(staging_dir, "matrices")
    copy_if_exists(
        [
            os.path.join(matrix_dir, "*summary_dist.csv"),
            os.path.join(matrix_dir, "*cumulative_dist.csv"),
            os.path.join(matrix_dir, "*_tree_contributions.csv"),
            os.path.join(matrix_dir, "*similarity.csv"),
        ],
        matrices_out,
    )

    # Copy logs and predictions if present
    logs_out = os.path.join(staging_dir, "logs")
    copy_if_exists([f"{proj_name}_*.log", f"{proj_name}_*.out"], logs_out)

    # Save run summary and manifest
    summaries_dir = os.path.join(staging_dir, "summaries")
    os.makedirs(summaries_dir, exist_ok=True)
    summary_path = os.path.join(summaries_dir, "run_summary.json")
    summarize_run(proj_name, gene_names, tree_count, matrix_dir, exon_metric_rows, threshold, use_flag, use_threshold, summary_path)
    manifest_path = os.path.join(staging_dir, "manifest.txt")
    write_manifest(staging_dir, manifest_path)

    # Create ZIP archive
    os.makedirs(package_dir, exist_ok=True)
    archive_name = zip_name or f"sprout_{proj_name}_results.zip"
    zip_path = os.path.join(package_dir, archive_name)
    zip_staging(staging_dir, zip_path)
    print(f"GUI-ready bundle created at {zip_path}")


if __name__ == "__main__":
    main()
