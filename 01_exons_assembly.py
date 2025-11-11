#!/usr/bin/env python3
"""
01_exons_assembly.py – Assemble raw reads and extract exon sequences.
This script performs the following steps:
1. Quality trimming of raw reads using fastp.
2. Assembly of target sequences using HybPiper.
3. For each gene in the provided gene list:
    - Processes exon data from exonerate results.
    - Extracts exon sequences from assembled contigs.
    - Writes exon assignments to a TSV file and creates FASTA files for each exon.

It requires a config file for parameters, or command-line arguments can be used to override defaults.  
It is expected to be run in an environment with the necessary dependencies installed, including Biopython, pandas, and HybPiper.
It also assumes the presence of a shared utilities module (`pipeline_utils`) for logging and command execution
Arguments:
- `-c`, `--config`: Path to a configuration file (YAML/JSON/TOML).
- `-t`, `--threads`: Number of threads to use for assembly.
- `-r1`, `--read1`: Path to the first reads file (FASTQ).
- `-r2`, `--read2`: Path to the second reads file (FASTQ).
- `-m`, `--mega353`: Path to the target FASTA file (default is "angiosperms353_v2_interim_targetfile.fasta").
- `-p`, `--proj_name`: Project name identifier.
- `-g`, `--gene_list`: Path to the gene list file.
- `-ov`, `--overlap`: Overlap ratio to consider the same exon (default is 0.8).
- `--output_hyb`: Output folder for HybPiper results (default is "01_hyb_output").
- `--output_exon`: Output folder for exon FASTAs (default is "02_exon_extracted").

Example usage:
python 01_exons_assembly.py -c config.yaml -t 8 -r1 reads_1.fastq -r2 reads_2.fastq -p my_project -g gene_list.txt
or
python 01_exons_assembly.py --threads 8 --read1 reads_1.fastq --read2 reads_2.fastq --proj_name my_project --gene_list gene_list.txt
"""
import os
import sys
import argparse
import ast
import json
import pandas as pd
from Bio import SeqIO
# Import shared utilities
from pipeline_utils import log_status, run_command, is_valid_project_name, load_config


def parse_ranges_field(value):
    """Safely parse exon coordinate ranges stored as strings."""
    if isinstance(value, list):
        ranges = value
    elif isinstance(value, str):
        value = value.strip()
        if not value:
            return []
        try:
            ranges = ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return []
    else:
        return []
    parsed = []
    for item in ranges:
        if isinstance(item, (list, tuple)) and len(item) >= 2:
            try:
                start = int(float(item[0]))
                end = int(float(item[1]))
                parsed.append((start, end))
            except (TypeError, ValueError):
                continue
    return parsed


def coerce_numeric(value):
    """Convert a value or collection of values to a representative float if possible."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        try:
            return float(stripped)
        except ValueError:
            try:
                literal = ast.literal_eval(stripped)
            except (ValueError, SyntaxError):
                return None
            return coerce_numeric(literal)
    if isinstance(value, (list, tuple)):
        numeric_values = [coerce_numeric(v) for v in value]
        numeric_values = [v for v in numeric_values if v is not None]
        if numeric_values:
            return sum(numeric_values) / len(numeric_values)
        return None
    return None

def extract_contigs(row, fasta_sequences, output_dir):
    """Append contig segments for each exon in ``row`` to the appropriate FASTA file."""
    ranges = row.get('parsed_ranges', [])
    exon_names = row.get('exon_names', [])
    sequence_id = row.get('sequence_id')
    if sequence_id is None and len(row) > 3:
        sequence_id = row.iloc[3]
    # Find the full sequence record corresponding to this contig ID
    sequence = next((seq for seq in fasta_sequences if seq.id == sequence_id), None)
    if sequence is None:
        print(f"Warning: Sequence {sequence_id} not found in FASTA.")
        return
    for i, (start, end) in enumerate(ranges):
        if i >= len(exon_names):
            continue
        exon_name = exon_names[i]
        contig = sequence[start:end]
        contig.id = f"{exon_name}_{sequence_id}_{i+1}"
        contig.description = ""
        exon_file = os.path.join(output_dir, f"{exon_name}.fasta")
        with open(exon_file, "a") as fh:
            SeqIO.write(contig, fh, "fasta")


def clean_fasta(row, fasta_sequences, output_dir):
    """Remove stale exon FASTA files prior to writing updated contig sequences."""
    ranges = row.get('parsed_ranges', [])
    exon_names = row.get('exon_names', [])
    for i in range(min(len(ranges), len(exon_names))):
        exon_name = exon_names[i]
        exon_file = os.path.join(output_dir, f"{exon_name}.fasta")
        if os.path.exists(exon_file):
            os.remove(exon_file)

def check_overlap(exon_ranges, start, end, overlap_threshold):
    """
    Check if the interval (start, end) overlaps significantly (>= overlap_threshold) with any existing exon interval.
    Returns the name of an overlapping exon if found, otherwise None.
    """
    for (existing_start, existing_end), exon_name in exon_ranges:
        # Compute overlap length between intervals
        overlap_len = min(end, existing_end) - max(start, existing_start)
        # Check if overlap is at least the given fraction of the combined interval length
        if overlap_len > 0 and overlap_len / (max(end, existing_end) - min(start, existing_start)) >= overlap_threshold:
            return exon_name
    return None

def process_exon_data(input_dir, gene_name, output_dir, overlap_threshold):
    """Process HybPiper exon statistics for ``gene_name`` and extract contig FASTAs."""
    data_label = os.path.basename(input_dir)
    stats_file = os.path.join(input_dir, gene_name, data_label, 'exonerate_stats.tsv')
    if not os.path.isfile(stats_file):
        raise FileNotFoundError(f"Expected stats file not found: {stats_file}")

    df = pd.read_csv(stats_file, sep='\t')
    sentinel_idx = df[df.iloc[:, 0] == 'Hits filtered to remove hits with frameshifts'].index
    if len(sentinel_idx) > 0:
        df = df.loc[:sentinel_idx[0] - 1]
    df = df.dropna(how='all')

    output_tsv = os.path.join(output_dir, f"{data_label}_{gene_name}_exon_split.tsv")

    if df.empty:
        df.to_csv(output_tsv, sep='\t', index=False)
        return {
            "gene": gene_name,
            "assignment_table": output_tsv,
            "exon_fastas": [],
            "metrics": [],
        }

    range_col = None
    for col in df.columns:
        non_null = df[col].dropna().head(5)
        for value in non_null:
            if isinstance(value, (list, tuple)):
                range_col = col
                break
            if isinstance(value, str) and value.strip().startswith('['):
                range_col = col
                break
        if range_col:
            break
    if range_col is None:
        raise ValueError(f"Could not locate exon coordinate column in {stats_file}")

    df['parsed_ranges'] = df[range_col].apply(parse_ranges_field)

    candidate_sequence_cols = [
        col for col in df.columns
        if any(token in str(col).lower() for token in ('contig', 'sequence', 'seq', 'name', 'id'))
    ]
    candidate_sequence_cols += [col for col in df.columns if col not in candidate_sequence_cols]

    sequence_col = None
    for col in candidate_sequence_cols:
        if col in df.columns and df[col].notna().any():
            sequence_col = col
            break
    if sequence_col is None:
        sequence_col = df.columns[0]
    df['sequence_id'] = df[sequence_col].astype(str)

    depth_columns = [col for col in df.columns if 'depth' in str(col).lower()]
    score_columns = [col for col in df.columns if 'score' in str(col).lower()]

    exon_ranges = []
    exon_names_per_row = []
    metrics_lookup = {}
    exon_counter = 0

    for _, row in df.iterrows():
        ranges = row['parsed_ranges'] or []
        depth_value = None
        for col in depth_columns:
            depth_value = coerce_numeric(row.get(col))
            if depth_value is not None:
                break
        score_value = None
        for col in score_columns:
            score_value = coerce_numeric(row.get(col))
            if score_value is not None:
                break

        row_exon_names = []
        for start, end in ranges:
            overlap_exon = check_overlap(exon_ranges, start, end, overlap_threshold)
            if overlap_exon is None:
                exon_counter += 1
                exon_name = f"{data_label}_{gene_name}_exon_{exon_counter}"
                exon_ranges.append(((start, end), exon_name))
            else:
                exon_name = overlap_exon

            entry = metrics_lookup.setdefault(
                exon_name,
                {
                    'project': data_label,
                    'gene': gene_name,
                    'exon_name': exon_name,
                    'lengths': [],
                    'depth_values': [],
                    'score_values': [],
                },
            )
            length_bp = max(0, int(end) - int(start))
            if length_bp:
                entry['lengths'].append(length_bp)
            if depth_value is not None:
                entry['depth_values'].append(depth_value)
            if score_value is not None:
                entry['score_values'].append(score_value)
            row_exon_names.append(exon_name)

        exon_names_per_row.append(row_exon_names)
        print(f"{gene_name}: identified {len(row_exon_names)} exons in one alignment hit.")

    df['exon_names'] = exon_names_per_row
    df.to_csv(output_tsv, sep='\t', index=False)

    contigs_fasta = os.path.join(input_dir, gene_name, f"{gene_name}_contigs.fasta")
    if not os.path.isfile(contigs_fasta):
        raise FileNotFoundError(f"Expected FASTA file not found: {contigs_fasta}")
    fasta_sequences = list(SeqIO.parse(contigs_fasta, "fasta"))

    unique_exons = sorted({name for names in exon_names_per_row for name in names})
    for exon_name in unique_exons:
        exon_file = os.path.join(output_dir, f"{exon_name}.fasta")
        if os.path.exists(exon_file):
            os.remove(exon_file)

    for _, row in df.iterrows():
        if row['parsed_ranges']:
            extract_contigs(row, fasta_sequences, output_dir)

    exon_files = [
        os.path.join(output_dir, f"{exon}.fasta")
        for exon in unique_exons
        if os.path.exists(os.path.join(output_dir, f"{exon}.fasta"))
    ]

    metrics_records = []
    for entry in metrics_lookup.values():
        lengths = entry['lengths']
        depth_values = entry['depth_values']
        score_values = entry['score_values']
        metrics_records.append(
            {
                'project': entry['project'],
                'gene': entry['gene'],
                'exon_name': entry['exon_name'],
                'length_bp': max(lengths) if lengths else None,
                'mean_depth': (sum(depth_values) / len(depth_values)) if depth_values else None,
                'alignment_score': (sum(score_values) / len(score_values)) if score_values else None,
            }
        )

    return {
        "gene": gene_name,
        "assignment_table": output_tsv,
        "exon_fastas": exon_files,
        "metrics": metrics_records,
    }

def sequence_assembly(num_threads, read1, read2, target_fasta, project, log_file, output_hyb_dir):
    """
    Step 1: Run quality trimming and assembly:
    - Uses fastp for read trimming.
    - Runs HybPiper to assemble target sequences from trimmed reads.
    """
    # 1A. Read trimming with fastp
    fastp_cmd = (
        f"fastp -i {read1} -I {read2} "
        f"-o {read1}.trimmed.fastq.gz -O {read2}.trimmed.fastq.gz "
        f"-j fastp.json -h fastp.html"
    )
    run_command(fastp_cmd, "Sequence Trimming (fastp)", log_file, critical=True)
    # 1B. Run HybPiper assembly
    os.makedirs(output_hyb_dir, exist_ok=True)
    log_status(log_file, f"Create Output Directory ({output_hyb_dir}): SUCCESS")
    hybpiper_cmd = (
        f"hybpiper assemble -t_dna {target_fasta} "
        f"-r {read1}.trimmed.fastq.gz {read2}.trimmed.fastq.gz "
        f"--prefix {project} --bwa --cpu {num_threads} -o {output_hyb_dir}"
    )
    run_command(hybpiper_cmd, "Sequence Assembly (HybPiper)", log_file, critical=True)
    trimmed_read1 = os.path.abspath(f"{read1}.trimmed.fastq.gz")
    trimmed_read2 = os.path.abspath(f"{read2}.trimmed.fastq.gz")
    return {
        "trimmed_reads": [trimmed_read1, trimmed_read2],
        "fastp_reports": [os.path.abspath("fastp.json"), os.path.abspath("fastp.html")],
        "hybpiper_output": os.path.abspath(output_hyb_dir),
    }

def exon_extraction(gene_list_path, overlap_threshold, project, log_file, input_hyb_dir, output_exon_dir):
    """
    Step 2: For each gene in the gene list, process exons and extract contigs.
    """
    # Read gene names (strip any “.fasta” extension in the list if present)
    with open(gene_list_path, 'r') as f:
        gene_names = [line.strip().replace('.fasta', '') for line in f if line.strip()]
    # Save a cleaned gene list for compatibility with downstream steps
    with open('gene_list.txt', 'w') as f_out:
        for gene in gene_names:
            f_out.write(gene + '\n')
    log_status(log_file, "Modify Gene List: SUCCESS")
    os.makedirs(output_exon_dir, exist_ok=True)
    log_status(log_file, f"Create Output Directory ({output_exon_dir}): SUCCESS")
    input_project_dir = os.path.join(input_hyb_dir, project)
    assignment_tables = []
    exon_fastas = []
    all_metrics = []
    metrics_csv = None
    metrics_json = None
    metrics_manifest = None
    for gene in gene_names:
        try:
            gene_outputs = process_exon_data(input_project_dir, gene, output_exon_dir, overlap_threshold)
            if gene_outputs:
                assignment_table = gene_outputs.get("assignment_table")
                if assignment_table and os.path.exists(assignment_table):
                    assignment_tables.append(os.path.abspath(assignment_table))
                exon_fastas.extend(
                    os.path.abspath(path)
                    for path in gene_outputs.get("exon_fastas", [])
                    if os.path.exists(path)
                )
                all_metrics.extend(gene_outputs.get("metrics", []))
            log_status(log_file, f"Processed Exons for Gene {gene}: SUCCESS")
        except Exception as e:
            log_status(log_file, f"Failed to Process Exons for Gene {gene}: {e}: FAILURE")
            print(f"Error processing exons for gene {gene}: {e}")
    if all_metrics:
        metrics_df = pd.DataFrame(all_metrics)
        if not metrics_df.empty:
            metrics_df.sort_values(['gene', 'exon_name'], inplace=True)
        metrics_csv = os.path.join(output_exon_dir, f"{project}_exon_metrics.csv")
        metrics_df.to_csv(metrics_csv, index=False)
        records = []
        for record in metrics_df.to_dict(orient='records'):
            cleaned_record = {key: (None if pd.isna(value) else value) for key, value in record.items()}
            records.append(cleaned_record)
        metrics_json = os.path.join(output_exon_dir, f"{project}_exon_metrics.json")
        with open(metrics_json, 'w') as fh:
            json.dump(records, fh, indent=2)
        metrics_manifest = os.path.join(output_exon_dir, f"{project}_exon_metrics_manifest.json")
        manifest = {
            'project': project,
            'metrics_csv': metrics_csv,
            'metrics_json': metrics_json,
            'gene_count': int(metrics_df['gene'].nunique()),
            'exon_count': int(len(metrics_df)),
        }
        with open(metrics_manifest, 'w') as fh:
            json.dump(manifest, fh, indent=2)
        log_status(log_file, f"Exon metrics written to {metrics_csv}")
        log_status(log_file, f"Exon metrics JSON written to {metrics_json}")
        log_status(log_file, f"Exon metrics manifest written to {metrics_manifest}")
    else:
        log_status(log_file, "No exon metrics were generated.")
    return {
        "assignment_tables": sorted(set(assignment_tables)),
        "exon_fastas": sorted(set(exon_fastas)),
        "gene_list": os.path.abspath('gene_list.txt'),
        "metrics_csv": os.path.abspath(metrics_csv) if metrics_csv else None,
        "metrics_json": os.path.abspath(metrics_json) if metrics_json else None,
        "metrics_manifest": os.path.abspath(metrics_manifest) if metrics_manifest else None,
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assemble reads and extract exons from a mixed sample.")
    parser.add_argument("-c", "--config", help="Path to config file (YAML/JSON/TOML)")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to use")
    parser.add_argument("-r1", "--read1", help="Path to first reads file (FASTQ)")
    parser.add_argument("-r2", "--read2", help="Path to second reads file (FASTQ)")
    parser.add_argument("-m", "--mega353", help="Path to target FASTA (default Angiosperms353)", 
                        default="angiosperms353_v2_interim_targetfile.fasta")
    parser.add_argument("-p", "--proj_name", help="Project name identifier")
    parser.add_argument("-g", "--gene_list", help="Path to gene list file")
    parser.add_argument("-ov", "--overlap", type=float, help="Overlap ratio to consider same exon (0-1)",
                        default=0.8)
    parser.add_argument("--output_hyb", help="Output folder for HybPiper results", default="01_hyb_output")
    parser.add_argument("--output_exon", help="Output folder for exon FASTAs", default="02_exon_extracted")
    args = parser.parse_args()

    # Load config file if provided
    config = {}
    if args.config:
        config = load_config(args.config)
    # Merge CLI args and config, giving precedence to CLI
    threads = args.threads if args.threads is not None else config.get('threads')
    read1 = args.read1 or config.get('read1')
    read2 = args.read2 or config.get('read2')
    mega353 = args.mega353 or config.get('mega353', "angiosperms353_v2_interim_targetfile.fasta")
    proj_name = args.proj_name or config.get('proj_name')
    gene_list = args.gene_list or config.get('gene_list')
    overlap = args.overlap if args.overlap is not None else config.get('overlap', 0.8)
    output_hyb = args.output_hyb or config.get('output_hyb', "01_hyb_output")
    output_exon = args.output_exon or config.get('output_exon', "02_exon_extracted")
    # Validate required params
    if threads is None or not read1 or not read2 or not proj_name or not gene_list:
        parser.error("Missing required parameters (threads, read1, read2, proj_name, gene_list).")
    if not is_valid_project_name(proj_name):
        sys.exit(f"Error: Project name '{proj_name}' is invalid (contains disallowed characters).")
    threads = int(threads)
    overlap = float(overlap)
    # Initialize log file
    log_file = f"{proj_name}_01_exon_assembly.out"
    if os.path.exists(log_file):
        os.remove(log_file)
    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {threads}")
    log_status(log_file, f"  Read1: {read1}")
    log_status(log_file, f"  Read2: {read2}")
    log_status(log_file, f"  Mega353: {mega353}")
    log_status(log_file, f"  Project Name: {proj_name}")
    log_status(log_file, f"  Gene List: {gene_list}")
    log_status(log_file, f"  Overlap Threshold: {overlap}")
    log_status(log_file, f"  Output Hyb: {output_hyb}")
    log_status(log_file, f"  Output Exon: {output_exon}")
    # Run steps 1 and 2
    assembly_outputs = sequence_assembly(threads, read1, read2, mega353, proj_name, log_file, output_hyb)
    exon_outputs = exon_extraction(gene_list, overlap, proj_name, log_file, output_hyb, output_exon)
    stage_manifest = {
        "stage": 1,
        "project": proj_name,
        "log_file": os.path.abspath(log_file),
        "parameters": {
            "threads": threads,
            "read1": os.path.abspath(read1),
            "read2": os.path.abspath(read2),
            "target_fasta": os.path.abspath(mega353),
            "gene_list_source": os.path.abspath(gene_list),
            "overlap_threshold": overlap,
            "output_hyb": os.path.abspath(output_hyb),
            "output_exon": os.path.abspath(output_exon),
        },
        "artifacts": {
            "trimmed_reads": assembly_outputs.get("trimmed_reads", []),
            "fastp_reports": assembly_outputs.get("fastp_reports", []),
            "hybpiper_output": assembly_outputs.get("hybpiper_output"),
            "exon_assignment_tables": exon_outputs.get("assignment_tables", []),
            "exon_fastas": exon_outputs.get("exon_fastas", []),
            "clean_gene_list": exon_outputs.get("gene_list"),
            "exon_metrics_csv": exon_outputs.get("metrics_csv"),
            "exon_metrics_json": exon_outputs.get("metrics_json"),
            "exon_metrics_manifest": exon_outputs.get("metrics_manifest"),
        },
        "downstream_inputs": {
            "stage2": {
                "exon_fasta_directory": os.path.abspath(output_exon),
                "gene_list": exon_outputs.get("gene_list"),
            }
        },
    }
    manifest_path = f"{proj_name}_stage1_manifest.json"
    with open(manifest_path, "w") as manifest_file:
        json.dump(stage_manifest, manifest_file, indent=2)
    log_status(log_file, f"Stage 1 manifest saved to {manifest_path}")
    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for details.")
