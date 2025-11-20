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
from Bio.SeqRecord import SeqRecord
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
        try:
            if isinstance(item, (list, tuple)) and len(item) == 2:
                start = int(item[0])
                end = int(item[1])
                if start < end:
                    normalized.append((start, end))
        except (ValueError, TypeError):
            continue
    return normalized


def check_overlap(existing_exons, start, end, threshold):
    """
    Determine if exon (start, end) overlaps with any existing exon above a given threshold.
    `existing_exons` is a list of tuples: [((start, end), exon_name), ...].
    Returns the exon_name if overlapping exon is found, otherwise None.
    """
    for (exon_start, exon_end), exon_name in existing_exons:
        overlap_start = max(start, exon_start)
        overlap_end = min(end, exon_end)
        overlap_len = max(0, overlap_end - overlap_start)
        exon_len = exon_end - exon_start
        if exon_len > 0 and (overlap_len / exon_len) >= threshold:
            return exon_name
    return None


def coerce_numeric(value):
    """
    Try to coerce a value (possibly a string or list) to a float.
    Returns None if this isn't possible.
    """
    if value is None:
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
    # Handle both dict and Series objects from pandas
    def get_value(row, key, default=None):
        try:
            if hasattr(row, 'get'):
                return row.get(key, default)
            else:
                return row[key] if key in row else default
        except (KeyError, TypeError):
            return default
    
    ranges = get_value(row, 'parsed_ranges', [])
    exon_names = get_value(row, 'exon_names', [])
    sequence_id = get_value(row, 'sequence_id')
    
    # Debug output
    if not ranges or not exon_names:
        print(f"Debug: ranges={ranges}, exon_names={exon_names}, sequence_id={sequence_id}")
        return
    
    # Find the sequence with the matching ID
    sequence = next((seq for seq in fasta_sequences if seq.id == str(sequence_id)), None)
    if sequence is None:
        print(f"Warning: Sequence {sequence_id} not found in FASTA. Available sequences: {[seq.id for seq in fasta_sequences[:5]]}...")
        return

    # Write each exon segment to its respective FASTA file
    for i, (start, end) in enumerate(ranges):
        try:
            exon_name = exon_names[i]
        except IndexError:
            print(f"Warning: exon_names index mismatch for sequence_id={sequence_id}, ranges={ranges}, exon_names={exon_names}")
            continue
        contig = sequence[start:end]
        contig.id = f"{exon_name}_{sequence_id}_{i+1}"
        contig.description = ""
        exon_file = os.path.join(output_dir, f"{exon_name}.fasta")
        with open(exon_file, "a") as fh:
            SeqIO.write(contig, fh, "fasta")


def clean_fasta(row, fasta_sequences, output_dir):
    """
    Remove existing exon FASTA files for the exons present in the given row.
    This prevents old data from previous genes from accumulating in the files.
    """
    exon_names = row.get('exon_names', []) if hasattr(row, 'get') else row['exon_names']
    if not exon_names:
        return
    unique_exons = set(exon_names)
    for exon in unique_exons:
        exon_file = os.path.join(output_dir, f"{exon}.fasta")
        if os.path.exists(exon_file):
            os.remove(exon_file)


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

    # --------- FIXED COLUMN SELECTION (matches the original working script) ----------
    # Use fixed column indices to match the original behaviour.
    # Column index 6 (7th column) holds the exon coordinate ranges, typically
    # something like "query_HSPFragment_ranges", containing a Python-like list
    # of (start, end) tuples.
    df['parsed_ranges'] = df.iloc[:, 6].apply(parse_ranges_field)

    # Column index 3 (4th column) holds the contig / hit identifier that
    # matches the record IDs in <gene>_contigs.fasta (usually the 'hit_id').
    # We cast to string to be safe and to ensure matching against FASTA ids.
    df['sequence_id'] = df.iloc[:, 3].astype(str)
    # -------------------------------------------------------------------------------

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

            row_exon_names.append(exon_name)

            if exon_name not in metrics_lookup:
                metrics_lookup[exon_name] = {
                    'project': data_label,
                    'gene': gene_name,
                    'exon_name': exon_name,
                    'lengths': [],
                    'depth_values': [],
                    'score_values': [],
                }
            metrics_lookup[exon_name]['lengths'].append(end - start)
            if depth_value is not None:
                metrics_lookup[exon_name]['depth_values'].append(depth_value)
            if score_value is not None:
                metrics_lookup[exon_name]['score_values'].append(score_value)

        exon_names_per_row.append(row_exon_names)
        print(f"{gene_name}: identified {len(row_exon_names)} exons in one alignment hit.")

    df['exon_names'] = exon_names_per_row
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"{gene_name}: Wrote assignment table with {len(df)} rows to {output_tsv}")
    print(f"{gene_name}: Columns in dataframe: {list(df.columns)}")
    print(f"{gene_name}: Sample exon_names: {exon_names_per_row[:2] if exon_names_per_row else 'EMPTY'}")

    contigs_fasta = os.path.join(input_dir, gene_name, f"{gene_name}_contigs.fasta")
    if not os.path.isfile(contigs_fasta):
        raise FileNotFoundError(f"Expected FASTA file not found: {contigs_fasta}")
    fasta_sequences = list(SeqIO.parse(contigs_fasta, "fasta"))
    print(f"{gene_name}: Loaded {len(fasta_sequences)} contigs from {contigs_fasta}")

    for _, row in df.iterrows():
        clean_fasta(row, fasta_sequences, output_dir)
    for _, row in df.iterrows():
        extract_contigs(row, fasta_sequences, output_dir)

    unique_exons = sorted({exon for exon_list in exon_names_per_row for exon in exon_list})
    exon_files = [
        os.path.join(output_dir, f"{exon}.fasta")
        for exon in unique_exons
        if os.path.exists(os.path.join(output_dir, f"{exon}.fasta"))
    ]
    
    print(f"Gene {gene_name}: Created {len(exon_files)} exon FASTA files out of {len(unique_exons)} expected exons")
    if len(exon_files) < len(unique_exons):
        missing_exons = [e for e in unique_exons if not os.path.exists(os.path.join(output_dir, f"{e}.fasta"))]
        print(f"Warning: Missing FASTA files for exons: {missing_exons}")
    else:
        print(f"Success: All {len(exon_files)} exon FASTA files created successfully!")

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
    Step 1: Run fastp for quality trimming, then HybPiper for assembly.
    """
    log_status(log_file, "Starting sequence assembly")
    # 1A. Quality trimming with fastp
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


def exon_extraction(gene_list_path, overlap_threshold, project, log_file, input_hyb_dir, output_exon_dir):
    """
    Step 2: For each gene in the gene list, process exons and extract contigs.
    """
    # Read gene names
    with open(gene_list_path, "r") as f:
        gene_names = [line.strip() for line in f if line.strip()]
    # Save gene list in output folder
    gene_list_output = os.path.join(output_exon_dir, "gene_list.txt")
    with open(gene_list_output, "w") as f_out:
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
        metrics_csv = os.path.join(output_exon_dir, f"{project}_exon_metrics.csv")
        metrics_df.to_csv(metrics_csv, index=False)
        metrics_json = os.path.join(output_exon_dir, f"{project}_exon_metrics.json")
        metrics_df.to_json(metrics_json, orient="records", indent=2)
        metrics_manifest = os.path.join(output_exon_dir, f"{project}_exon_metrics_manifest.json")
        manifest_data = {
            "project": project,
            "metric_table": os.path.abspath(metrics_csv),
            "metric_json": os.path.abspath(metrics_json),
        }
        with open(metrics_manifest, "w") as mf:
            json.dump(manifest_data, mf, indent=2)
        log_status(log_file, f"Saved exon metrics to {metrics_csv} and {metrics_json}")
    else:
        log_status(log_file, "No exon metrics to save.")
    return {
        "gene_list": gene_list_output,
        "assignment_tables": assignment_tables,
        "exon_fastas": exon_fastas,
        "metrics_csv": metrics_csv,
        "metrics_json": metrics_json,
        "metrics_manifest": metrics_manifest,
    }


def main():
    parser = argparse.ArgumentParser(description="01_exons_assembly – Assemble reads and extract exons.")
    parser.add_argument(
        "-c", "--config", type=str, required=False,
        help="Path to configuration file (YAML/JSON/TOML)."
    )
    parser.add_argument("--threads", "-t", type=int, default=None, help="Number of CPU threads.")
    parser.add_argument("--read1", "-r1", type=str, default=None, help="Path to read1 FASTQ.")
    parser.add_argument("--read2", "-r2", type=str, default=None, help="Path to read2 FASTQ.")
    parser.add_argument("--target_fasta", "-tf", type=str, default=None, help="Target FASTA for HybPiper.")
    parser.add_argument("--proj_name", "-p", type=str, default=None, help="Project name (sample identifier).")
    parser.add_argument("--gene_list", "-g", type=str, default=None, help="Gene list file.")
    parser.add_argument("--overlap_threshold", "-o", type=float, default=None,
                        help="Overlap threshold for exons (0–1).")
    parser.add_argument("--output_hyb_dir", type=str, default=None, help="Output directory for HybPiper.")
    parser.add_argument("--output_exon_dir", type=str, default=None, help="Output directory for exon FASTAs.")
    parser.add_argument("--log_file", type=str, default="01_exons_assembly.log", help="Log file.")

    args = parser.parse_args()

    config = load_config(args.config) if args.config else {}
    num_threads = args.threads or config.get("threads", 4)
    read1 = args.read1 or config.get("read1")
    read2 = args.read2 or config.get("read2")
    target_fasta = args.target_fasta or config.get("target_fasta")
    project = args.proj_name or config.get("proj_name")
    gene_list = args.gene_list or config.get("gene_list")
    overlap_threshold = args.overlap_threshold or config.get("overlap_threshold", 0.8)
    output_hyb = args.output_hyb_dir or config.get("output_hyb_dir", "01_HybPiper_output")
    output_exon = args.output_exon_dir or config.get("output_exon_dir", "02_Exons_output")
    log_file = args.log_file or config.get("log_file", "01_exons_assembly.log")

    if not is_valid_project_name(project):
        print(f"Error: Invalid project name '{project}'.")
        sys.exit(1)
    if not read1 or not read2 or not target_fasta or not gene_list:
        print("Error: Missing required inputs (read1, read2, target_fasta, gene_list).")
        sys.exit(1)

    log_status(log_file, "Starting 01_exons_assembly stage.")
    os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)

    sequence_assembly(num_threads, read1, read2, target_fasta, project, log_file, output_hyb)

    exon_outputs = exon_extraction(
        gene_list_path=gene_list,
        overlap_threshold=overlap_threshold,
        project=project,
        log_file=log_file,
        input_hyb_dir=output_hyb,
        output_exon_dir=output_exon,
    )

    stage_manifest = {
        "stage": "01_exons_assembly",
        "project": project,
        "log_file": os.path.abspath(log_file),
        "outputs": {
            "hybpiper_dir": os.path.abspath(output_hyb),
            "exon_dir": os.path.abspath(output_exon),
            "gene_list": exon_outputs.get("gene_list"),
            "assignment_tables": exon_outputs.get("assignment_tables", []),
            "exon_fastas": exon_outputs.get("exon_fastas", []),
            "metrics_csv": exon_outputs.get("metrics_csv"),
            "metrics_json": exon_outputs.get("metrics_json"),
            "metrics_manifest": exon_outputs.get("metrics_manifest"),
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


if __name__ == "__main__":
    main()
