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
- `--output_tree`: Directory for output matrices (default: "04_all_trees").

Usage:
python 03_distance_matrices.py -c config.yaml -t 4 -p my_project -g gene_list.txt --threshold 1.96 --use_flag --input_phylo 03_phylo_results --output_tree 04_all_trees
or
python 03_distance_matrices.py --config config.yaml --threads 4 --proj_name my_project --gene_list gene_list.txt --threshold 1.96 --use_flag --input_phylo 03_phylo_results --output_tree 04_all_trees
"""
import os
import glob
import json
import argparse
import shutil
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
            # Flag method: keep only the maximum value (best match) as 0, set others to 999
            max_val = df[col].max()
            df[col] = df[col].apply(lambda x: 0 if x == max_val else 999)
        elif use_threshold:
            # Threshold method: mark as 999 any value higher than (mean - threshold*std)
            mean = df[col].mean()
            sd = df[col].std()
            if sd > 0:
                # For similarities, keep values above (mean + threshold*std) to identify good matches
                threshold_val = mean + threshold * sd
                df[col] = df[col].apply(lambda x: x if x > threshold_val else 999)
            else:
                # If std is 0, keep values above mean
                df[col] = df[col].apply(lambda x: x if x >= mean else 999)
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
        
        row_name_col = df.columns[0]
        for header in df.columns[1:]:
            if header in species_to_taxa:
                taxa_list = species_to_taxa[header]
                # Mark as 999 any row whose taxon is not in the expected taxa list for this species
                mask = ~df[row_name_col].astype(str).isin(taxa_list)
                df.loc[mask, header] = 999
    # Simplify row names by removing any trailing numbers/underscores (from original sample IDs)
    df.iloc[:, 0] = df.iloc[:, 0].str.replace(r'\d+', '', regex=True).str.rstrip('_')
    return df

def _sum_reindexed_matrices(matrices):
    """Return the element-wise sum of similarity matrices, aligning indexes/columns."""
    if not matrices:
        return pd.DataFrame()

    all_index = sorted({idx for df in matrices for idx in df.index})
    all_columns = sorted({col for df in matrices for col in df.columns})
    combined = pd.DataFrame(0.0, index=all_index, columns=all_columns)
    for df in matrices:
        reindexed = df.reindex(index=all_index, columns=all_columns, fill_value=0.0)
        combined += reindexed
    return combined


def _aggregate_manifest_entries(entries, aggregate_by, matrix_dir):
    """Aggregate manifest entries according to the requested level."""
    if aggregate_by == "tree":
        return entries, []

    aggregated = []
    aggregated_paths = []
    key_field = "gene" if aggregate_by == "gene" else "exon"
    for key, group in sorted(((k, v) for k, v in _group_entries(entries, key_field).items()), key=lambda item: item[0]):
        matrices = [entry.pop("_dataframe") for entry in group]
        combined_df = _sum_reindexed_matrices(matrices)
        if combined_df.empty:
            total = 0.0
        else:
            total = float(combined_df.to_numpy().sum())

        if aggregate_by == "gene":
            gene = key
            exon = "all"
            identifier = f"{gene}.all"
            filename = f"{gene}.aggregated.similarity.json"
        else:
            gene = "all"
            exon = key
            identifier = f"exon_{key}"
            filename = f"{identifier}.aggregated.similarity.json"

        aggregate_path = os.path.join(matrix_dir, filename)
        combined_df.to_json(aggregate_path, orient='split')
        aggregated_paths.append(os.path.abspath(aggregate_path))
        aggregated.append({
            "id": identifier,
            "gene": gene,
            "exon": exon,
            "matrix": filename,
            "total": total,
            "member_count": len(group),
            "_dataframe": combined_df,
        })
    return aggregated, aggregated_paths


def _group_entries(entries, field):
    groups = {}
    for entry in entries:
        key = entry[field]
        groups.setdefault(key, []).append(entry)
    return groups


def process_matrices(matrix_dir, project, threshold, use_flag, use_threshold, aggregate_by="tree"):
    """
    Combine all per-gene distance matrices in `matrix_dir` into one summary DataFrame.
    Converts distances to similarities, persists per-tree similarity matrices, and
    generates a manifest with normalized weights for downstream visualization.

    Returns
    -------
    tuple
        (summary_df, manifest_path, similarity_jsons, aggregated_similarity_jsons)
    """
    all_dfs = []
    manifest_entries = []
    similarity_paths = []
    aggregated_paths = []
    for filename in os.listdir(matrix_dir):
        if filename.endswith('cleaned.csv'):
            file_path = os.path.join(matrix_dir, filename)
            df = pd.read_csv(file_path)
            # Identify corresponding taxa list file (if exists) for further filtering
            prefix = filename.split('cleaned.csv')[0]
            taxa_file = os.path.join(matrix_dir, f"{prefix}list.txt")
            df = clean_up_matrix(df, project, threshold, taxa_file if os.path.exists(taxa_file) else None, use_flag, use_threshold)
            df = distance_to_similarity(df)
            row_name_col = df.columns[0]
            df = df.rename(columns={row_name_col: 'row_name'})
            similarity_df = df.set_index('row_name')

            base = prefix.rstrip('.')
            if '.' in base:
                gene_name, exon_id = base.split('.', 1)
            else:
                gene_name, exon_id = base, '1'
            similarity_filename = f"{base}.similarity.json"
            similarity_path = os.path.join(matrix_dir, similarity_filename)
            similarity_df.to_json(similarity_path, orient='split')
            similarity_paths.append(os.path.abspath(similarity_path))

            tree_total = float(similarity_df.to_numpy().sum())
            manifest_entries.append({
                "id": f"{gene_name}.{exon_id}",
                "gene": gene_name,
                "exon": exon_id,
                "matrix": similarity_filename,
                "total": tree_total,
                "member_count": 1,
                "_dataframe": similarity_df,
            })
            all_dfs.append(df)

    if manifest_entries:
        aggregated_entries, aggregated_paths = _aggregate_manifest_entries(manifest_entries, aggregate_by, matrix_dir)
        if aggregate_by == "tree":
            aggregated_entries = manifest_entries
            aggregated_paths = []
        grand_total = sum(entry["total"] for entry in aggregated_entries)
        for entry in aggregated_entries:
            entry["weight"] = entry["total"] / grand_total if grand_total else 0.0
            entry.pop("_dataframe", None)
        manifest_entries_sorted = sorted(aggregated_entries, key=lambda item: item.get("id", ""))
        manifest_data = {
            "project": project,
            "aggregate_by": aggregate_by,
            "entries": [{
                **{k: v for k, v in entry.items() if k != "_dataframe"},
                "total": float(entry["total"]),
                "weight": float(entry["weight"]),
            } for entry in manifest_entries_sorted],
        }
    else:
        manifest_data = {
            "project": project,
            "aggregate_by": aggregate_by,
            "entries": [],
        }

    manifest_path = os.path.join(matrix_dir, f"{project}.manifest.json")
    with open(manifest_path, 'w') as manifest_file:
        json.dump(manifest_data, manifest_file, indent=2)

    manifest_path_abs = os.path.abspath(manifest_path)
    if not all_dfs:
        return pd.DataFrame(columns=['row_name', 'total_value']), manifest_path_abs, sorted(similarity_paths), sorted(aggregated_paths)

    # Concatenate all gene similarity data
    total_df = pd.concat(all_dfs, ignore_index=True)
    total_df.fillna(0, inplace=True)
    # Sum similarity scores across all genes for each taxon
    value_cols = [col for col in total_df.columns if col != 'row_name']
    total_df['total_value'] = total_df[value_cols].sum(axis=1)
    # Return a DataFrame with taxon (row_name) and its aggregated total value
    result = total_df[['row_name', 'total_value']]
    return result, manifest_path_abs, sorted(similarity_paths), sorted(aggregated_paths)

def genetic_distance_matrix(tree_file, node_output_file, output_file):
    """
    Given a Newick tree file, compute the genetic distance matrix and save to output_file.
    Also, record sister taxa for any collapsed nodes in node_output_file.
    Reroot the tree using known outgroups if present; otherwise, use midpoint rooting.

    Parameters:
    - tree_file: Path to the input Newick tree file.
    - node_output_file: Path to output file for recording sister taxa of collapsed nodes.
    - output_file: Path to output CSV file for the genetic distance matrix.
    Returns:
    None
    """
    tree = Phylo.read(tree_file, 'newick')
    reroot_taxa = ["Amborella", "Nymphaea", "Austrobaileya"]
    
    # Try to root with known outgroups
    rooted = False
    for taxa in reroot_taxa:
        for clade in tree.find_clades():
            if clade.name and taxa in clade.name:
                try:
                    tree.root_with_outgroup(clade)
                    rooted = True
                    break
                except Exception as e:
                    print(f"Warning: Could not root with {taxa}: {e}")
        if rooted:
            break
    
    # If no outgroup found, use midpoint rooting
    if not rooted:
        try:
            tree.root_at_midpoint()
        except Exception as e:
            print(f"Warning: Midpoint rooting failed: {e}")
    
    node_records = []
    for tip in tree.get_terminals():
        if tip.name and "NODE" in tip.name:
            recorded_taxa = find_clade_and_move(tree, tip.name)
            if recorded_taxa:
                node_records.append(f'{tip.name}: {"; ".join(recorded_taxa)}')
    with open(node_output_file, 'w') as file:
        for record in node_records:
            file.write(f'{record}\n')
    clades, distance_matrix = calculate_genetic_distance(tree)
    df = pd.DataFrame(distance_matrix)
    df.columns = clades
    df.index = clades
    df.to_csv(output_file)

def group_and_sum(input_file, output_file):
    """
    Group the input CSV by 'row_name', summing 'total_value' for each unique taxon.
    Save the aggregated results to output_file.
    Parameters:
    - input_file: Path to the input CSV file with 'row_name' and 'total_value' columns.
    - output_file: Path to the output CSV file for aggregated results.
    Returns:
    None
    """
    data = pd.read_csv(input_file)
    filtered_data = data[~data['row_name'].str.contains("NODE")]
    grouped_sum = filtered_data.groupby('row_name')['total_value'].sum().reset_index()
    grouped_sum.to_csv(output_file, index=False)

def process_gene(gene_name_shorter, input_dir, output_dir, log_file):
    """Generate per-tree matrices for one gene and return created files."""
    artifacts = {
        "gene": gene_name_shorter,
        "tree_files": [],
        "node_lists": [],
        "raw_matrices": [],
        "cleaned_matrices": [],
        "errors": None,
    }
    try:
        # list all .tre files for this gene inside input_dir
        pattern = os.path.join(input_dir, f"{gene_name_shorter}*tre")
        tree_files = sorted(glob.glob(pattern))
        if not tree_files:
            log_status(log_file, f"[WARN] No trees found for {gene_name_shorter} in {input_dir}")
            return artifacts

        for i, filename in enumerate(tree_files, start=1):
            base = f"{gene_name_shorter}.{i}"
            node_output_file = os.path.join(output_dir, f"{base}.list.txt")
            matrix_file      = os.path.join(output_dir, f"{base}.matrix")
            cleaned_csv      = os.path.join(output_dir, f"{base}.cleaned.csv")

            genetic_distance_matrix(filename, node_output_file, matrix_file)
            log_status(log_file, f"Generated matrix for {gene_name_shorter} tree {i}")
            artifacts["tree_files"].append(os.path.abspath(filename))
            artifacts["node_lists"].append(os.path.abspath(node_output_file))
            artifacts["raw_matrices"].append(os.path.abspath(matrix_file))

            # The genetic_distance_matrix writes a square CSV; copy/rename to *.cleaned.csv
            # (keep a separate 'cleaned' name because downstream code searches for it)
            shutil.copy2(matrix_file, cleaned_csv)
            log_status(log_file, f"Copied matrix to cleaned CSV for {gene_name_shorter} tree {i}")
            artifacts["cleaned_matrices"].append(os.path.abspath(cleaned_csv))

        log_status(log_file, f"Finished {gene_name_shorter}")
    except Exception as e:
        log_status(log_file, f"Failed processing {gene_name_shorter}: {e}")
        print(f"Failed processing {gene_name_shorter}: {e}")
        artifacts["errors"] = str(e)
    return artifacts

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
    parser.add_argument("--output_tree", help="Directory for output matrices", default="04_all_trees")
    parser.add_argument(
        "--aggregate_by",
        choices=["tree", "gene", "exon"],
        default="tree",
        help="Aggregation level for GUI data manifest (default: tree)",
    )
    args = parser.parse_args()

    # Load config if provided
    config = {}
    if args.config:
        config = load_config(args.config)
    threads = args.threads if args.threads is not None else config.get('threads')
    proj_name = args.proj_name or config.get('proj_name')
    gene_list_path = args.gene_list or config.get('gene_list', "gene_list.txt")
    threshold = args.threshold if args.threshold is not None else config.get('threshold', 1.96)
    use_flag = args.use_flag or bool(config.get('use_flag', False))
    use_threshold = args.use_threshold or bool(config.get('use_threshold', False))
    input_phylo = args.input_phylo or config.get('input_phylo', "03_phylo_results")
    output_tree = args.output_tree or config.get('output_tree', "04_all_trees")
    aggregate_by = args.aggregate_by or config.get('aggregate_by', "tree")

    if aggregate_by not in {"tree", "gene", "exon"}:
        parser.error("--aggregate_by must be one of: tree, gene, exon")

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
    log_status(log_file, f"  Output Directory: {output_tree}")
    log_status(log_file, f"  Aggregate Level: {aggregate_by}")
    os.makedirs(output_tree, exist_ok=True)
    log_status(log_file, f"Created directory {output_tree}")
    # Load gene names and process each gene's trees in parallel
    with open(gene_list_path, 'r') as f:
        gene_names = [line.strip() for line in f if line.strip()]
    gene_results = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_gene, gene, input_phylo, output_tree, log_file) for gene in gene_names]
        for future in futures:
            result = future.result()
            if result:
                gene_results.append(result)
    # Combine all matrices and output summary
    summary_df, similarity_manifest, similarity_files, aggregated_similarity = process_matrices(
        output_tree,
        proj_name,
        threshold,
        use_flag,
        use_threshold,
        aggregate_by,
    )
    summary_csv = os.path.join(output_tree, f"{proj_name}.summary_dist.csv")
    summary_df.to_csv(summary_csv, index=False)
    log_status(log_file, f"Processed matrices saved to {summary_csv}")
    # Generate cumulative distance summary across all genes
    cumulative_csv = f"{proj_name}.cumulative_dist.csv"
    group_and_sum(summary_csv, cumulative_csv)
    log_status(log_file, f"Generated cumulative distance file: {cumulative_csv}")
    manifest_artifacts = {
        "tree_files": sorted({path for res in gene_results for path in res.get("tree_files", [])}),
        "node_lists": sorted({path for res in gene_results for path in res.get("node_lists", [])}),
        "raw_matrices": sorted({path for res in gene_results for path in res.get("raw_matrices", [])}),
        "cleaned_matrices": sorted({path for res in gene_results for path in res.get("cleaned_matrices", [])}),
        "summary_matrix": os.path.abspath(summary_csv),
        "cumulative_matrix": os.path.abspath(cumulative_csv),
        "similarity_manifest": similarity_manifest,
        "similarity_matrices": similarity_files,
        "aggregated_similarity_matrices": aggregated_similarity,
    }
    stage_manifest = {
        "stage": 3,
        "project": proj_name,
        "log_file": os.path.abspath(log_file),
        "parameters": {
            "threads": threads,
            "gene_list": os.path.abspath(gene_list_path),
            "threshold": threshold,
            "use_flag": use_flag,
            "use_threshold": use_threshold,
            "input_phylo_dir": os.path.abspath(input_phylo),
            "output_tree_dir": os.path.abspath(output_tree),
        },
        "artifacts": manifest_artifacts,
        "downstream_inputs": {
            "stage4": {
                "cumulative_matrix": manifest_artifacts["cumulative_matrix"],
                "summary_matrix": manifest_artifacts["summary_matrix"],
            }
        },
    }
    manifest_path = f"{proj_name}_stage3_manifest.json"
    with open(manifest_path, "w") as manifest_file:
        json.dump(stage_manifest, manifest_file, indent=2)
    log_status(log_file, f"Stage 3 manifest saved to {manifest_path}")
    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for details.")
