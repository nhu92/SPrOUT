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
import pandas as pd
from Bio import SeqIO
# Import shared utilities
from pipeline_utils import log_status, run_command, is_valid_project_name, load_config

def extract_contigs(row, fasta_sequences, output_dir):
    """
    Extract contig sequences for each exon in the DataFrame row and append to exon-specific FASTA files.
    """
    ranges = eval(row.iloc[13])  # exon coordinate ranges
    sequence_id = row.iloc[3]
    # Find the full sequence record corresponding to this contig ID
    sequence = next((seq for seq in fasta_sequences if seq.id == sequence_id), None)
    if sequence is None:
        print(f"Warning: Sequence {sequence_id} not found in FASTA.")
        return
    # Write each exon segment to its respective FASTA file
    for i, (start, end) in enumerate(ranges):
        exon_name = row['exon_names'][i]
        contig = sequence[start:end]              # extract subsequence for the exon
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
    ranges = eval(row.iloc[13])
    sequence_id = row.iloc[3]
    sequence = next((seq for seq in fasta_sequences if seq.id == sequence_id), None)
    if sequence is None:
        return
    for i, _ in enumerate(ranges):
        exon_name = row['exon_names'][i]
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
    """
    Process exonerate results for a single gene to define unique exons and extract corresponding contigs.
    - Reads the exonerate_stats.tsv for the gene.
    - Determines exon names (ensuring overlapping hits get the same name).
    - Writes a TSV of exon assignments and creates FASTA files for each exon.
    """
    data_label = os.path.basename(input_dir)  # project or sample identifier (folder name)
    stats_file = os.path.join(input_dir, gene_name, data_label, 'exonerate_stats.tsv')
    if not os.path.isfile(stats_file):
        raise FileNotFoundError(f"Expected stats file not found: {stats_file}")
    df = pd.read_csv(stats_file, sep='\t')
    # Truncate DataFrame at the line indicating frameshift filtering (if present)
    end_idx = df[df.iloc[:, 0] == 'Hits filtered to remove hits with frameshifts'].index
    if len(end_idx) > 0:
        df = df.loc[:end_idx[0]-1]
    exon_ranges = []  # list of ((start, end), exon_name) for discovered exons
    exon_names_per_row = []  # exon name list for each alignment hit (row)
    for _, row in df.iterrows():
        ranges = eval(row.iloc[6])  # parse stringified list of exon coordinates
        row_exon_names = []
        for i, (start, end) in enumerate(ranges):
            overlap_exon = check_overlap(exon_ranges, start, end, overlap_threshold)
            if overlap_exon is None:
                # Define a new exon
                exon_name = f"{data_label}_{gene_name}_exon_{i+1}"
                exon_ranges.append(((start, end), exon_name))
            else:
                # Use the existing exon name for overlapping region
                exon_name = overlap_exon
            row_exon_names.append(exon_name)
        exon_names_per_row.append(row_exon_names)
        print(f"{gene_name}: identified {len(row_exon_names)} exons in one alignment hit.")
    df['exon_names'] = exon_names_per_row
    # Save exon assignment table for reference
    output_tsv = os.path.join(output_dir, f"{data_label}_{gene_name}_exon_split.tsv")
    df.to_csv(output_tsv, sep='\t', index=False)
    # Load assembled contigs FASTA for this gene (from HybPiper output)
    contigs_fasta = os.path.join(input_dir, gene_name, f"{gene_name}_contigs.fasta")
    fasta_sequences = list(SeqIO.parse(contigs_fasta, "fasta"))
    # Remove any old entries in exon FASTA files for this gene, then extract new sequences
    for _, row in df.iterrows():
        clean_fasta(row, fasta_sequences, output_dir)
    for _, row in df.iterrows():
        extract_contigs(row, fasta_sequences, output_dir)

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
    for gene in gene_names:
        try:
            process_exon_data(input_project_dir, gene, output_exon_dir, overlap_threshold)
            log_status(log_file, f"Processed Exons for Gene {gene}: SUCCESS")
        except Exception as e:
            log_status(log_file, f"Failed to Process Exons for Gene {gene}: {e}: FAILURE")
            print(f"Error processing exons for gene {gene}: {e}")

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
    sequence_assembly(threads, read1, read2, mega353, proj_name, log_file, output_hyb)
    exon_extraction(gene_list, overlap, proj_name, log_file, output_hyb, output_exon)
    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for details.")
