#!/usr/bin/env python3
"""
03_distance_matrices.py – Compute genetic distance matrices from gene trees and aggregate results.
This script processes phylogenetic trees for each gene, calculating pairwise genetic distances
and generating a summary distance matrix for all taxa.

It can be run with command-line arguments or a configuration file.

Arguments:
- `-c`, `--config`: Path to configuration file (YAML/JSON/TOML).
- `-t`, `--threads`: Number of threads for parallel processing. Defaults to 1 if not specified.
- `-p`, `--proj_name`: Project name identifier for output files.
- `-g`, `--gene_list`: Path to file containing list of gene names.  Defaults to "gene_list.txt".
- `--threshold`: Threshold for distance filtering (default: 1.96). 
- `--use_flag`: Use flag method for filtering (min=0, others=999).
- `--use_threshold`: Enable threshold-based filtering (default: off).
- `--input_phylo`: Directory containing input .tre files (default: "03_phylo_results").
- `--output_dir`: Directory for output matrices (default: "04_all_trees").

Usage:
python 03_distance_matrices.py -c config.yaml -t 4 -p my_project -g gene_list.txt --threshold 1.96 --use_flag --input_phylo 03_phylo_results --output_dir 04_all_trees
or
python 03_distance_matrices.py --config config.yaml --threads 4 --proj_name my_project --gene_list gene_list.txt --threshold 1.96 --use_flag --input_phylo 03_phylo_results --output_dir 04_all_trees
"""
import os
import glob
import argparse
from concurrent.futures import ThreadPoolExecutor
from Bio import Phylo
import pandas as pd
import numpy as np
from pipeline_utils import log_status, load_config, is_valid_project_name

def find_clade_and_move(tree, node_name):
    """
    Find sister taxa for a collapsed node (e.g., "NODE_x"): return names of all non-NODE tips 
    in the sister clade of the smallest clade containing node_name, if support > 0.7.
    """
    target = None
    for leaf in tree.get_terminals():
        if leaf.name == node_name:
            target = leaf
            break
    related_taxa = []
    if target:
        # Traverse from target leaf up toward root
        path = tree.get_path(target)
        for clade in reversed(path):
            if clade.clades:  # if not a terminal
                # Check each sibling clade for real (non-"NODE") taxa
                for sister in clade.clades:
                    if sister is not clade and any("NODE" not in leaf.name for leaf in sister.get_terminals()):
                        related_taxa.extend([leaf.name for leaf in sister.get_terminals() if "NODE" not in leaf.name])
                        break
                # If this clade had any sister taxa, and clade support is >0.7, stop climbing up
                if related_taxa and (clade.confidence is None or clade.confidence > 0.7):
                    break
    return related_taxa

def calculate_genetic_distance(tree):
    """
    Calculate pairwise distances between all leaves in the tree.
    Returns a tuple of (list_of_taxa, distance_matrix_numpy).
    """
    taxa = [leaf.name for leaf in tree.get_terminals()]
    n = len(taxa)
    distances = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = tree.distance(taxa[i], taxa[j])
            distances[i, j] = distances[j, i] = d
    return taxa, distances

def distance_to_similarity(dist_df):
    """Transform a distance DataFrame into similarity values using 1/(1+d)."""
    sim_df = dist_df.copy()
    numeric_cols = sim_df.select_dtypes(include=['float', 'int']).columns
    sim_df[numeric_cols] = 1 / (1 + sim_df[numeric_cols])
    return sim_df

def clean_up_matrix(df, project, threshold, taxa_file=None, use_flag=False, use_threshold=True):
    """
    Clean a distance matrix DataFrame by filtering out irrelevant entries:
    - Remove rows where the taxon name contains the project name (i.e., self-hits).
    - Keep only columns (species) that belong to the project.
    - Apply either standard deviation threshold filtering or "flag" method to mark outliers as 999.
    - If a taxa mapping file is provided, mask distances for species that do not co-occur with expected taxa.
    """
    # Exclude any rows that correspond to the project’s own sequences
    df = df[~df.iloc[:, 0].str.contains(project)]
    # Keep only columns where header contains the project name (plus the first column for row labels)
    cols_to_keep = [df.columns[0]] + [col for col in df.columns[1:] if project in col]
    df = df[cols_to_keep]
    # Apply filtering to each project column
    for col in df.columns[1:]:
        if use_flag:
            # Flag method: set the minimum value in this column to 0 (best match) and all others to 999
            min_val = df[col].min()
            df[col] = df[col].apply(lambda x: 0 if x == min_val else 999)
        elif use_threshold:
            # Standard deviation method: mark as 999 any value higher than (mean - threshold*std) for the column
            mean = df[col].mean()
            sd = df[col].std()
            df[col] = df[col].apply(lambda x: 999 if x > (mean - threshold * sd) else x)
        # If neither flag nor threshold, do not filter
    # If a taxa mapping file is available, use it to further clean the matrix
    if taxa_file and os.path.exists(taxa_file):
        species_to_taxa = {}
        with open(taxa_file, 'r') as tf:
            for line in tf:
                parts = line.strip().split(':')
                if len(parts) == 2:
                    species, taxa_list = parts
                    species_to_taxa[species.strip()] = [tax.strip() for tax in taxa_list.split(';')]
        for header in df.columns[1:]:
            for species, taxa_list in species_to_taxa.items():
                if header == species:
                    # For each row, if its taxon (row_name) is not in the expected taxa list for this species, mark distance as 999
                    for idx, row in df.iterrows():
                        if row[df.columns[0]] not in taxa_list:
                            df.at[idx, header] = 999
    # Simplify row names by removing any trailing numbers/underscores (from original sample IDs)
    df.iloc[:, 0] = df.iloc[:, 0].str.replace(r'\d+', '', regex=True).str.rstrip('_')
    return df

def process_matrices(matrix_dir, project, threshold, use_flag, use_threshold):
    """
    Combine all per-gene distance matrices in `matrix_dir` into one summary DataFrame.
    Converts distances to similarities and computes a total sum per taxon.
    """
    all_dfs = []
    for filename in os.listdir(matrix_dir):
        if filename.endswith('cleaned.csv'):
            file_path = os.path.join(matrix_dir, filename)
            df = pd.read_csv(file_path)
            # Identify corresponding taxa list file (if exists) for further filtering
            prefix = filename.split('cleaned.csv')[0]
            taxa_file = os.path.join(matrix_dir, f"{prefix}list.txt")
            df = clean_up_matrix(df, project, threshold, taxa_file if os.path.exists(taxa_file) else None, use_flag, use_threshold)
            df = distance_to_similarity(df)
            all_dfs.append(df)
    # Concatenate all gene similarity data
    total_df = pd.concat(all_dfs, ignore_index=True)
    total_df.fillna(0, inplace=True)
    # Sum similarity scores across all genes for each taxon
    value_cols = [col for col in total_df.columns if col != 'Unnamed: 0']
    total_df['total_value'] = total_df[value_cols].sum(axis=1)
    # Return a DataFrame with taxon (row_name) and its aggregated total value
    result = total_df[['Unnamed: 0', 'total_value']].rename(columns={'Unnamed: 0': 'row_name'})
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute distance matrices from exon trees and aggregate them.")
    parser.add_argument("-c", "--config", help="Path to config file (YAML/JSON/TOML)")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads for parallel processing")
    parser.add_argument("-p", "--proj_name", help="Project name identifier")
    parser.add_argument("-g", "--gene_list", help="Path to gene list file", default="gene_list.txt")
    parser.add_argument("--threshold", type=float, help="Threshold for distance filtering", default=1.96)
    parser.add_argument("--use_flag", action="store_true", help="Use flag method (min=0, others=999) for filtering")
    parser.add_argument("--use_threshold", action="store_true", help="Enable threshold-based filtering (default: off)")
    parser.add_argument("--input_phylo", help="Directory with input .tre files", default="03_phylo_results")
    parser.add_argument("--output_dir", help="Directory for output matrices", default="04_all_trees")
    args = parser.parse_args()

    # Load config if provided
    config = {}
    if args.config:
        config = load_config(args.config)
    threads = args.threads if args.threads is not None else config.get('threads')
    proj_name = args.proj_name or config.get('proj_name')
    gene_list_path = args.gene_list if args.gene_list != parser.get_default('gene_list') else config.get('gene_list', "gene_list.txt")
    threshold = args.threshold if args.threshold is not None else config.get('threshold', 1.96)
    use_flag = args.use_flag or bool(config.get('use_flag', False))
    use_threshold = args.use_threshold or bool(config.get('use_threshold', False))
    input_phylo = args.input_phylo if args.input_phylo != parser.get_default('input_phylo') else config.get('input_phylo', "03_phylo_results")
    output_dir = args.output_dir if args.output_dir != parser.get_default('output_dir') else config.get('output_dir', "04_all_trees")

    # Conflict check: use_flag and use_threshold cannot both be True
    if use_flag and use_threshold:
        parser.error("You cannot enable both --use_flag and --use_threshold at the same time.")

    if threads is None or not proj_name:
        parser.error("Required parameters missing: threads and proj_name must be specified.")
    if not is_valid_project_name(proj_name):
        parser.error(f"Project name '{proj_name}' contains invalid characters.")
    threads = int(threads)
    threshold = float(threshold)
    # Initialize log
    log_file = f"{proj_name}_03_distance_calc.log"
    if os.path.exists(log_file):
        os.remove(log_file)
    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {threads}")
    log_status(log_file, f"  Project Name: {proj_name}")
    log_status(log_file, f"  Gene List: {gene_list_path}")
    log_status(log_file, f"  Threshold: {threshold} (enabled: {use_threshold})")
    log_status(log_file, f"  Use Flag: {use_flag}")
    log_status(log_file, f"  Input Directory: {input_phylo}")
    log_status(log_file, f"  Output Directory: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    log_status(log_file, f"Created directory {output_dir}")
    # Load gene names and process each gene's trees in parallel
    with open(gene_list_path, 'r') as f:
        gene_names = [line.strip() for line in f if line.strip()]
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_gene, gene, input_phylo, output_dir, log_file) for gene in gene_names]
        for future in futures:
            future.result()  # raise exceptions if any
    # Combine all matrices and output summary
    summary_df = process_matrices(output_dir, proj_name, threshold, use_flag, use_threshold)
    summary_csv = os.path.join(output_dir, f"{proj_name}.summary_dist.csv")
    summary_df.to_csv(summary_csv, index=False)
    log_status(log_file, f"Processed matrices saved to {summary_csv}")
    # Generate cumulative distance summary across all genes
    cumulative_csv = f"{proj_name}.cumulative_dist.csv"
    group_and_sum(summary_csv, cumulative_csv)
    log_status(log_file, f"Generated cumulative distance file: {cumulative_csv}")
    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for details.")
