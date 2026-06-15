#!/usr/bin/env python3
"""
04_prediction.py – Summarize total scores by taxonomic level and filter by significance.
This script processes cumulative scores from a CSV file, summarizes them by taxonomic level,
and filters taxa based on a z-score threshold.
It supports command-line arguments or a configuration file for flexibility.

Arguments:
    -c, --config: Path to a configuration file (YAML/JSON/TOML).
    -i, --input_file: Input CSV file containing cumulative scores.
    -o, --output_file: Output CSV file for summarized scores by taxonomy.
    -tl, --taxonomic_level: Taxonomic level to summarize (o = Order, f = Family, g = Genus, s = Species).
    -z, --zscore_threshold: Z-score threshold for filtering significant taxa.
    -to, --taxonomy_output_file: Output file for selected taxonomy names based on z-score.

Usage:
    python 04_prediction.py -i input_scores.csv -o summary_scores.csv -tl g -z 2.0 -to selected_taxa.txt
    or
    python 04_prediction.py --config config.yaml
"""
import argparse
import pandas as pd
from scipy.stats import zscore
from pathlib import Path
from pipeline_utils import load_config, package_outputs
import os
import json

def process_column(column, level):
    """
    Extract the specified taxonomic level from each string in the column.
    Levels: 'o' = Order, 'f' = Family, 'g' = Genus, 's' = Species (Genus + species).
    """
    level_map = {'o': 0, 'f': 1, 'g': 2, 's': [2, 3]}
    processed = []
    for item in column:
        parts = item.split('_')
        if level == 's':
            # Species: join genus and species epithet
            processed.append('_'.join(parts[2:4]) if len(parts) >= 4 else item)
        else:
            idx = level_map[level]
            if isinstance(idx, list):
                # Not used in current levels (only 's' uses list)
                name = '_'.join(parts[idx[0]:idx[-1]+1])
            else:
                name = parts[idx] if len(parts) > idx else item
            processed.append(name)
    return processed

def select_taxonomy_by_zscore(summary_df, z_threshold):
    """Return list of taxon names from summary_df where z_score > z_threshold."""
    if 'z_score' not in summary_df.columns:
        raise ValueError("summary_df must contain 'z_score' column")
    if 'row_name' not in summary_df.columns:
        raise ValueError("summary_df must contain 'row_name' column")
    
    selected = summary_df.loc[summary_df['z_score'] > z_threshold, 'row_name'].tolist()
    return selected

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize and filter taxa by total scores and z-score.')
    parser.add_argument('-c', '--config', help='Path to config file (YAML/JSON/TOML)')
    parser.add_argument('-i', '--input_file', help='Input CSV file (cumulative scores)')
    parser.add_argument('-o', '--output_file', help='Output CSV for summarized scores by taxonomy')
    parser.add_argument('-tl', '--taxonomic_level', choices=['o', 'f', 'g', 's'], help='Taxonomic level (o, f, g, s)')
    parser.add_argument('-z', '--zscore_threshold', type=float, help='Z-score threshold for significance')
    parser.add_argument('-to', '--taxonomy_output_file', help='Output file for selected taxonomy names')
    parser.add_argument('-p', '--project_name', help='Project identifier used for packaging metadata.')
    parser.add_argument('--bundle-output', help='Destination directory or archive for the results bundle.')
    parser.add_argument('--bundle-format', choices=['directory', 'zip'], help='Bundle format: directory structure or zip archive.')
    parser.add_argument('--skip-bundle', action='store_true', help='Skip packaging outputs into a bundle.')
    parser.add_argument('--bundle-overwrite', action='store_true', help='Overwrite the bundle destination if it already exists.')
    args = parser.parse_args()

    # Load config if given
    config = {}
    if args.config:
        config = load_config(args.config)
    # Gather parameters (CLI or config)
    input_file = args.input_file or config.get('input_file')
    output_file = args.output_file or config.get('output_file')
    taxonomic_level = args.taxonomic_level or config.get('taxonomic_level')
    z_threshold = args.zscore_threshold if args.zscore_threshold is not None else config.get('zscore_threshold')
    taxonomy_output = args.taxonomy_output_file or config.get('taxonomy_output_file')
    project_name = args.project_name or config.get('project_name')
    bundle_format = args.bundle_format or config.get('bundle_format', 'directory')
    bundle_output = args.bundle_output or config.get('bundle_output')
    skip_bundle = args.skip_bundle or config.get('skip_bundle', False)
    bundle_overwrite = args.bundle_overwrite or config.get('bundle_overwrite', False)

    if not input_file or not output_file or not taxonomic_level or z_threshold is None or not taxonomy_output:
        parser.error("Parameters missing: input_file, output_file, taxonomic_level, zscore_threshold, taxonomy_output_file are required.")
    z_threshold = float(z_threshold)

    # Read input data and compute summary by taxonomic level
    df = pd.read_csv(input_file)
    
    # Ensure total_value column is numeric
    if 'total_value' in df.columns:
        df['total_value'] = pd.to_numeric(df['total_value'], errors='coerce')
    else:
        # If total_value doesn't exist, use the second column
        df.rename(columns={df.columns[1]: 'total_value'}, inplace=True)
        df['total_value'] = pd.to_numeric(df['total_value'], errors='coerce')
    
    # Remove rows with NaN total_value (from failed conversions)
    if df['total_value'].isnull().any():
        print(f"Warning: {df['total_value'].isnull().sum()} rows had non-numeric total_value and were excluded")
        df = df.dropna(subset=['total_value'])
    
    df['taxon_level'] = process_column(df.iloc[:, 0], taxonomic_level)
    summary = df.groupby('taxon_level')['total_value'].sum().reset_index()
<<<<<<< ours
<<<<<<< ours
<<<<<<< ours
<<<<<<< ours
<<<<<<< ours
    
    # Ensure summary['total_value'] is numeric before zscore
    summary['total_value'] = pd.to_numeric(summary['total_value'], errors='coerce')
    
    # Only calculate zscore if we have numeric data
    if summary['total_value'].dtype in [float, int, 'float64', 'int64']:
        summary['z_score'] = zscore(summary['total_value'])
    else:
        print(f"Error: total_value column has dtype {summary['total_value'].dtype}, expected numeric")
        raise ValueError(f"Cannot calculate z-score on non-numeric data: {summary['total_value'].dtype}")
    
=======
=======
>>>>>>> theirs
=======
>>>>>>> theirs
=======
>>>>>>> theirs
=======
>>>>>>> theirs
    if len(summary) > 1 and summary['total_value'].std(ddof=0) != 0:
        summary['z_score'] = zscore(summary['total_value'])
    else:
        # scipy.stats.zscore returns NaN for a single value or zero variance,
        # which breaks threshold-driven downstream reporting. Treat all tied
        # summaries as neutral rather than significant.
        summary['z_score'] = 0.0
<<<<<<< ours
<<<<<<< ours
<<<<<<< ours
<<<<<<< ours
>>>>>>> theirs
=======
>>>>>>> theirs
=======
>>>>>>> theirs
=======
>>>>>>> theirs
=======
>>>>>>> theirs
    summary = summary.rename(columns={'taxon_level': 'row_name', 'total_value': 'sum_of_total_value'})
    summary = summary.sort_values(by='sum_of_total_value', ascending=False)
    summary.to_csv(output_file, index=False)
    print(f"Summary has been written to {output_file}")

    # Filter by z-score threshold and save selected taxa
    significant_taxa = select_taxonomy_by_zscore(summary, z_threshold)
    
    if not significant_taxa:
        print(f"Warning: No taxa found with z_score > {z_threshold}")
    
    with open(taxonomy_output, 'w') as fout:
        for name in significant_taxa:
            fout.write(str(name) + "\n")
    print(f"Selected {len(significant_taxa)} taxa (z_score > {z_threshold}) written to {taxonomy_output}")
    inferred_project = config.get('proj_name') if isinstance(config, dict) else None
    if not inferred_project:
        base_name = os.path.basename(input_file)
        inferred_project = base_name.split('.')[0] if base_name else "project"
    stage_manifest = {
        "stage": 4,
        "project": inferred_project,
        "parameters": {
            "input_file": os.path.abspath(input_file),
            "output_file": os.path.abspath(output_file),
            "taxonomy_output_file": os.path.abspath(taxonomy_output),
            "taxonomic_level": taxonomic_level,
            "zscore_threshold": z_threshold,
        },
        "artifacts": {
            "taxonomy_summary": os.path.abspath(output_file),
            "selected_taxa": os.path.abspath(taxonomy_output),
        },
        "source_inputs": {
            "cumulative_matrix": os.path.abspath(input_file),
        },
        "downstream_inputs": {},
    }
    manifest_path = f"{inferred_project}_stage4_manifest.json"
    with open(manifest_path, "w") as manifest_file:
        json.dump(stage_manifest, manifest_file, indent=2)
    print(f"Stage 4 manifest saved to {manifest_path}")
