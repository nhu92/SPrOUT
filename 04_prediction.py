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
from pipeline_utils import load_config

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
    return summary_df.loc[summary_df['z_score'] > z_threshold, 'row_name'].tolist()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize and filter taxa by total scores and z-score.')
    parser.add_argument('-c', '--config', help='Path to config file (YAML/JSON/TOML)')
    parser.add_argument('-i', '--input_file', help='Input CSV file (cumulative scores)')
    parser.add_argument('-o', '--output_file', help='Output CSV for summarized scores by taxonomy')
    parser.add_argument('-tl', '--taxonomic_level', choices=['o', 'f', 'g', 's'], help='Taxonomic level (o, f, g, s)')
    parser.add_argument('-z', '--zscore_threshold', type=float, help='Z-score threshold for significance')
    parser.add_argument('-to', '--taxonomy_output_file', help='Output file for selected taxonomy names')
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
    if not input_file or not output_file or not taxonomic_level or z_threshold is None or not taxonomy_output:
        parser.error("Parameters missing: input_file, output_file, taxonomic_level, zscore_threshold, taxonomy_output_file are required.")
    z_threshold = float(z_threshold)

    # Read input data and compute summary by taxonomic level
    df = pd.read_csv(input_file)
    df['taxon_level'] = process_column(df.iloc[:, 0], taxonomic_level)
    summary = df.groupby('taxon_level')['total_value'].sum().reset_index()
    summary['z_score'] = zscore(summary['total_value'])
    summary = summary.rename(columns={'taxon_level': 'row_name', 'total_value': 'sum_of_total_value'})
    summary = summary.sort_values(by='sum_of_total_value', ascending=False)
    summary.to_csv(output_file, index=False)
    print(f"Summary has been written to {output_file}")

    # Filter by z-score threshold and save selected taxa
    significant_taxa = select_taxonomy_by_zscore(summary, z_threshold)
    with open(taxonomy_output, 'w') as fout:
        for name in significant_taxa:
            fout.write(name + "\n")
    print(f"Selected taxonomy names (z_score > {z_threshold}) have been written to {taxonomy_output}")
