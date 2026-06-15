<<<<<<< ours
<<<<<<< ours
<<<<<<< ours
#!/usr/bin/env python3
"""
02_exon_trees.py – Align exon sequences and build gene trees.
This script processes exon sequences for each gene, aligning them to a reference
alignment, trimming the alignment, and constructing a phylogenetic tree.

It supports parallel processing of multiple genes and logs the status of each step.
It requires the following tools:
- MAFFT for sequence alignment
- trimAl for trimming alignments
- FastTree for phylogenetic tree construction
It can be run with command-line arguments or a configuration file.

Arguments:
- -c, --config: Path to configuration file (YAML/JSON/TOML)
- -t, --threads: Number of CPU threads to use for parallel processing
- -e, --input_exon: Directory of extracted exon FASTA files
- -r, --ref_alignment: Directory of reference alignments
- -g, --gene_list: Path to gene list file
- -p, --proj_name: Project name identifier
- -o, --output_phylo: Directory for output trees
- -m, --min_exon_size: Minimum exon length to include in analysis
- --tree_method: Phylogenetic tree construction method (fasttree or iqtree)
- --iqtree_mode: IQ-TREE mode (fixed, fixed+gamma, or mfp)

Usage:
python 02_exon_trees.py -c config.yaml -t 8 -e 02_exon_extracted -r ref -g gene_list.txt -p my_project -o 03_phylo_results -m 80 --tree_method fasttree --iqtree_mode fixed
or
python 02_exon_trees.py --threads 8 --input_exon 02_exon_extracted --ref_alignment ref --gene_list gene_list.txt --proj_name my_project --output_phylo 03_phylo_results --min_exon_size 80 --tree_method fasttree --iqtree_mode fixed
"""
import os
import glob
import json
import argparse
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from pipeline_utils import log_status, run_command, load_config, is_valid_project_name

def all_sequences_meet_minimum_length(fasta_path, min_length=80):
    """Check if all sequences in the FASTA file are at least `min_length` bases long."""
    for record in SeqIO.parse(fasta_path, "fasta"):
        if len(record.seq) < min_length:
            return False
    return True

def process_gene_exon_alignment(
    gene_name, threads, input_dir, ref_dir, output_phylo, log_file, min_size,
    tree_method="fasttree", iqtree_mode="fixed"
):
    """
    Align and build tree for all exon contigs of a given gene.
    Supports FastTree and IQ-TREE (fixed or MFP mode).
    Args:
        gene_name: Name of the gene to process.
        threads: Number of threads to use for parallel processing.
        input_dir: Directory containing exon FASTA files.
        ref_dir: Directory containing reference alignments.
        output_phylo: Directory to save output trees.
        log_file: Log file to record status messages.
        min_size: Minimum length of exon sequences to include.
        tree_method: Phylogenetic tree construction method ('fasttree' or 'iqtree').
        iqtree_mode: IQ-TREE mode ('fixed', 'fixed+gamma', or 'mfp').
 
        """
    # Find all exon FASTA files for this gene
    gene_outputs = {
        "gene": gene_name,
        "aligned_fastas": [],
        "trimmed_fastas": [],
        "tree_files": [],
        "skipped_exons": [],
        "errors": None,
    }
    try:
        pattern = os.path.join(input_dir, f"*{gene_name}*.fasta")
        exon_files = glob.glob(pattern)
    except Exception as e:
        log_status(log_file, f"List exons for {gene_name}: FAILURE")
        print(f"Error listing exons for {gene_name}: {e}")
        gene_outputs["errors"] = str(e)
        return gene_outputs
    log_status(log_file, f"List exons for {gene_name}: SUCCESS")
    if not exon_files:
        log_status(log_file, f"No exon files found for {gene_name}, skipping.")
        return gene_outputs
    exon_files.sort()
    for i, exon_path in enumerate(exon_files, start=1):
        # Enforce minimum exon length
        if not all_sequences_meet_minimum_length(exon_path, min_size):
            log_status(log_file, f"Skipping {os.path.basename(exon_path)} (sequences < {min_size} bp)")
            gene_outputs["skipped_exons"].append(os.path.abspath(exon_path))
            continue
        # Alignment with MAFFT
        ref_alignment = os.path.join(ref_dir, f"{gene_name}.fasta")
        if not os.path.isfile(ref_alignment):
            log_status(log_file, f"Reference alignment not found: {ref_alignment}. Skipping {gene_name} exon {i}")
            gene_outputs["skipped_exons"].append(os.path.abspath(exon_path))
            continue
        aligned_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}_aligned.fasta")
        mafft_cmd = (
            f"mafft --preservecase --maxiterate 1000 --localpair --adjustdirection "
            f"--thread {threads} --addfragments {exon_path} {ref_alignment} > {aligned_out}"
        )
        run_command(mafft_cmd, f"MAFFT alignment for {gene_name} exon {i}", log_file)
        if not os.path.isfile(aligned_out):
            log_status(log_file, f"MAFFT failed to produce output for {gene_name} exon {i}")
            continue
        gene_outputs["aligned_fastas"].append(os.path.abspath(aligned_out))
        # Trim alignment with trimAl
        trimmed_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}_trimmed.fasta")
        trimal_cmd = f"trimal -in {aligned_out} -out {trimmed_out} -gt 0.5"
        run_command(trimal_cmd, f"Trim alignment for {gene_name} exon {i}", log_file)
        if not os.path.isfile(trimmed_out):
            log_status(log_file, f"trimAl failed to produce output for {gene_name} exon {i}")
            continue
        gene_outputs["trimmed_fastas"].append(os.path.abspath(trimmed_out))
        # Build tree with selected method
        tree_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}.tre")
        if tree_method == "fasttree":
            fasttree_cmd = f"fasttree -gtr -gamma -nt {trimmed_out} > {tree_out}"
            run_command(fasttree_cmd, f"Tree construction for {gene_name} exon {i} (FastTree)", log_file)
            if os.path.isfile(tree_out):
                gene_outputs["tree_files"].append(os.path.abspath(tree_out))
            else:
                log_status(log_file, f"FastTree failed to produce tree for {gene_name} exon {i}")
        elif tree_method == "iqtree":
            if iqtree_mode == "mfp":
                iqtree_cmd = (
                    f"iqtree2 -s {trimmed_out} -nt {threads} -m MFP -pre {tree_out.replace('.tre','')} "
                    f"--quiet"
                )
            elif iqtree_mode == "fixed+gamma":
                # GTR with gamma distribution
                iqtree_cmd = (
                    f"iqtree2 -s {trimmed_out} -nt {threads} -m GTR+G "
                    f"-pre {tree_out.replace('.tre','')} --quiet"
                )
            else:  # fixed model (GTR without gamma)
                iqtree_cmd = (
                    f"iqtree2 -s {trimmed_out} -nt {threads} -m GTR "
                    f"-pre {tree_out.replace('.tre','')} --quiet"
                )
            run_command(iqtree_cmd, f"Tree construction for {gene_name} exon {i} (IQ-TREE)", log_file)
            # IQ-TREE outputs .treefile, so rename/move to .tre for consistency
            iqtree_treefile = tree_out.replace('.tre', '.treefile')
            if os.path.exists(iqtree_treefile):
                os.replace(iqtree_treefile, tree_out)
                gene_outputs["tree_files"].append(os.path.abspath(tree_out))
            else:
                log_status(log_file, f"IQ-TREE failed to produce treefile for {gene_name} exon {i}")
        else:
            log_status(log_file, f"Unknown tree method: {tree_method}")
    return gene_outputs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align exons and build exon trees for each gene.")
    parser.add_argument("-c", "--config", help="Path to config file (YAML/JSON/TOML)")
    parser.add_argument("-t", "--threads", type=int, help="Number of CPU threads for parallel processing")
    parser.add_argument("-e", "--input_exon", help="Directory of extracted exon FASTA files", default="02_exon_extracted")
    parser.add_argument("-r", "--ref_alignment", help="Directory of reference alignments", default="ref")
    parser.add_argument("-g", "--gene_list", help="Path to gene list file", default="gene_list.txt")
    parser.add_argument("-p", "--proj_name", help="Project name identifier")
    parser.add_argument("-o", "--output_phylo", help="Directory for output trees", default="03_phylo_results")
    parser.add_argument("-m", "--min_exon_size", type=int, help="Minimum exon length to include", default=80)
    parser.add_argument("--tree_method", choices=["fasttree", "iqtree"], default="fasttree",
                        help="Phylogeny method: fasttree or iqtree (default: fasttree)")
    parser.add_argument("--iqtree_mode", choices=["fixed", "fixed+gamma", "mfp"], default="fixed",
                        help="IQ-TREE mode: fixed (GTR), fixed+gamma (GTR+G), or mfp (ModelFinder Plus)")
    args = parser.parse_args()

    # Load config if provided
    config = {}
    if args.config:
        config = load_config(args.config)
    # Determine parameters (CLI overrides config)
    threads = args.threads if args.threads is not None else config.get('threads')
    proj_name = args.proj_name or config.get('proj_name')
    input_exon_dir = args.input_exon or config.get('input_exon', "02_exon_extracted")
    ref_dir = args.ref_alignment or config.get('ref_alignment', "ref")
    gene_list_path = args.gene_list or config.get('gene_list', "gene_list.txt")
    output_phylo = args.output_phylo or config.get('output_phylo', "03_phylo_results")
    min_size = args.min_exon_size if args.min_exon_size is not None else config.get('min_exon_size', 80)
    tree_method = args.tree_method or config.get('tree_method', 'fasttree')
    iqtree_mode = args.iqtree_mode or config.get('iqtree_mode', 'fixed')
    # Disable iqtree_mode if tree_method is fasttree
    if tree_method == "fasttree":
        iqtree_mode = None
    if threads is None or not proj_name:
        parser.error("Required parameters missing: threads and proj_name must be specified.")
    if not is_valid_project_name(proj_name):
        parser.error(f"Project name '{proj_name}' contains invalid characters.")
    threads = int(threads)
    min_size = int(min_size)
    # Initialize log file
    log_file = f"{proj_name}_02_exons_phylo.log"
    if os.path.exists(log_file):
        os.remove(log_file)
    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {threads}")
    log_status(log_file, f"  Input Exon Directory: {input_exon_dir}")
    log_status(log_file, f"  Reference Alignment Directory: {ref_dir}")
    log_status(log_file, f"  Gene List: {gene_list_path}")
    log_status(log_file, f"  Project Name: {proj_name}")
    log_status(log_file, f"  Output Directory: {output_phylo}")
    log_status(log_file, f"  Minimum Exon Size: {min_size}")
    log_status(log_file, f"  Tree Method: {tree_method}")
    if tree_method == "iqtree":
        log_status(log_file, f"  IQ-TREE Mode: {iqtree_mode}")
    os.makedirs(output_phylo, exist_ok=True)
    log_status(log_file, f"Created directory {output_phylo}")
    # Read gene list and execute alignments/trees in parallel
    with open(gene_list_path, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    gene_results = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                process_gene_exon_alignment,
                gene, threads, input_exon_dir, ref_dir, output_phylo, log_file, min_size,
                tree_method, iqtree_mode
            )
            for gene in genes
        ]
        for future in futures:
            result = future.result()
            if result:
                gene_results.append(result)
    manifest_artifacts = {
        "aligned_fastas": sorted({path for res in gene_results for path in res.get("aligned_fastas", [])}),
        "trimmed_fastas": sorted({path for res in gene_results for path in res.get("trimmed_fastas", [])}),
        "tree_files": sorted({path for res in gene_results for path in res.get("tree_files", [])}),
        "skipped_exons": sorted({path for res in gene_results for path in res.get("skipped_exons", [])}),
    }
    stage_manifest = {
        "stage": 2,
        "project": proj_name,
        "log_file": os.path.abspath(log_file),
        "parameters": {
            "threads": threads,
            "input_exon_dir": os.path.abspath(input_exon_dir),
            "reference_alignment_dir": os.path.abspath(ref_dir),
            "gene_list": os.path.abspath(gene_list_path),
            "output_phylo_dir": os.path.abspath(output_phylo),
            "min_exon_size": min_size,
            "tree_method": tree_method,
            "iqtree_mode": iqtree_mode,
        },
        "artifacts": manifest_artifacts,
        "downstream_inputs": {
            "stage3": {
                "tree_directory": os.path.abspath(output_phylo),
                "tree_files": manifest_artifacts["tree_files"],
            }
        },
    }
    manifest_path = f"{proj_name}_stage2_manifest.json"
    with open(manifest_path, "w") as manifest_file:
        json.dump(stage_manifest, manifest_file, indent=2)
    log_status(log_file, f"Stage 2 manifest saved to {manifest_path}")
    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for details.")
=======
#!/usr/bin/env python3
"""
02_exon_trees.py – Align exon sequences and build gene trees.
This script processes exon sequences for each gene, aligning them to a reference
alignment, trimming the alignment, and constructing a phylogenetic tree.

It supports parallel processing of multiple genes and logs the status of each step.
It requires the following tools:
- MAFFT for sequence alignment
- trimAl for trimming alignments
- FastTree for phylogenetic tree construction
It can be run with command-line arguments or a configuration file.

Arguments:
- -c, --config: Path to configuration file (YAML/JSON/TOML)
- -t, --threads: Number of CPU threads to use for parallel processing
- -e, --input_exon: Directory of extracted exon FASTA files
- -r, --ref_alignment: Directory of reference alignments
- -g, --gene_list: Path to gene list file
- -p, --proj_name: Project name identifier
- -o, --output_phylo: Directory for output trees
- -m, --min_exon_size: Minimum exon length to include in analysis
- --tree_method: Phylogenetic tree construction method (fasttree or iqtree)
- --iqtree_mode: IQ-TREE mode (fixed, fixed+gamma, or mfp)

Usage:
python 02_exon_trees.py -c config.yaml -t 8 -e 02_exon_extracted -r ref -g gene_list.txt -p my_project -o 03_phylo_results -m 80 --tree_method fasttree --iqtree_mode fixed
or
python 02_exon_trees.py --threads 8 --input_exon 02_exon_extracted --ref_alignment ref --gene_list gene_list.txt --proj_name my_project --output_phylo 03_phylo_results --min_exon_size 80 --tree_method fasttree --iqtree_mode fixed
"""
import os
import glob
import argparse
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from pipeline_utils import log_status, run_command, load_config, is_valid_project_name

def all_sequences_meet_minimum_length(fasta_path, min_length=80):
    """Check if all sequences in the FASTA file are at least `min_length` bases long."""
    for record in SeqIO.parse(fasta_path, "fasta"):
        if len(record.seq) < min_length:
            return False
    return True

def process_gene_exon_alignment(
    gene_name, threads, input_dir, ref_dir, output_phylo, log_file, min_size,
    tree_method="fasttree", iqtree_mode="fixed"
):
    """
    Align and build tree for all exon contigs of a given gene.
    Supports FastTree and IQ-TREE (fixed or MFP mode).
    Args:
        gene_name: Name of the gene to process.
        threads: Number of threads to use for parallel processing.
        input_dir: Directory containing exon FASTA files.
        ref_dir: Directory containing reference alignments.
        output_phylo: Directory to save output trees.
        log_file: Log file to record status messages.
        min_size: Minimum length of exon sequences to include.
        tree_method: Phylogenetic tree construction method ('fasttree' or 'iqtree').
        iqtree_mode: IQ-TREE mode ('fixed', 'fixed+gamma', or 'mfp').
 
        """
    # Find all exon FASTA files for this gene
    try:
        pattern = os.path.join(input_dir, f"*{gene_name}*.fasta")
        exon_files = glob.glob(pattern)
    except Exception as e:
        log_status(log_file, f"List exons for {gene_name}: FAILURE")
        print(f"Error listing exons for {gene_name}: {e}")
        return
    log_status(log_file, f"List exons for {gene_name}: SUCCESS")
    if not exon_files:
        log_status(log_file, f"No exon files found for {gene_name}, skipping.")
        return
    exon_files.sort()
    for i, exon_path in enumerate(exon_files, start=1):
        # Enforce minimum exon length
        if not all_sequences_meet_minimum_length(exon_path, min_size):
            log_status(log_file, f"Skipping {os.path.basename(exon_path)} (sequences < {min_size} bp)")
            continue
        # Alignment with MAFFT
        ref_alignment = os.path.join(ref_dir, f"{gene_name}.fasta")
        aligned_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}_aligned.fasta")
        mafft_cmd = (
            f"mafft --preservecase --maxiterate 1000 --localpair --adjustdirection "
            f"--thread {threads} --addfragments {exon_path} {ref_alignment} > {aligned_out}"
        )
        run_command(mafft_cmd, f"MAFFT alignment for {gene_name} exon {i}", log_file)
        # Trim alignment with trimAl
        trimmed_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}_trimmed.fasta")
        trimal_cmd = f"trimal -in {aligned_out} -out {trimmed_out} -gt 0.5"
        run_command(trimal_cmd, f"Trim alignment for {gene_name} exon {i}", log_file)
        # Build tree with selected method
        tree_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}.tre")
        if tree_method == "fasttree":
            fasttree_cmd = f"fasttree -gtr -gamma -nt {trimmed_out} > {tree_out}"
            run_command(fasttree_cmd, f"Tree construction for {gene_name} exon {i} (FastTree)", log_file)
        elif tree_method == "iqtree":
            if iqtree_mode == "mfp":
                iqtree_cmd = (
                    f"iqtree2 -s {trimmed_out} -nt {threads} -m MFP -pre {tree_out.replace('.tre','')} "
                    f"--quiet"
                )
=======
#!/usr/bin/env python3
"""
02_exon_trees.py – Align exon sequences and build gene trees.
This script processes exon sequences for each gene, aligning them to a reference
alignment, trimming the alignment, and constructing a phylogenetic tree.

It supports parallel processing of multiple genes and logs the status of each step.
It requires the following tools:
- MAFFT for sequence alignment
- trimAl for trimming alignments
- FastTree for phylogenetic tree construction
It can be run with command-line arguments or a configuration file.

Arguments:
- -c, --config: Path to configuration file (YAML/JSON/TOML)
- -t, --threads: Number of CPU threads to use for parallel processing
- -e, --input_exon: Directory of extracted exon FASTA files
- -r, --ref_alignment: Directory of reference alignments
- -g, --gene_list: Path to gene list file
- -p, --proj_name: Project name identifier
- -o, --output_phylo: Directory for output trees
- -m, --min_exon_size: Minimum exon length to include in analysis
- --tree_method: Phylogenetic tree construction method (fasttree or iqtree)
- --iqtree_mode: IQ-TREE mode (fixed, fixed+gamma, or mfp)

Usage:
python 02_exon_trees.py -c config.yaml -t 8 -e 02_exon_extracted -r ref -g gene_list.txt -p my_project -o 03_phylo_results -m 80 --tree_method fasttree --iqtree_mode fixed
or
python 02_exon_trees.py --threads 8 --input_exon 02_exon_extracted --ref_alignment ref --gene_list gene_list.txt --proj_name my_project --output_phylo 03_phylo_results --min_exon_size 80 --tree_method fasttree --iqtree_mode fixed
"""
import os
import glob
import argparse
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from pipeline_utils import log_status, run_command, load_config, is_valid_project_name

def all_sequences_meet_minimum_length(fasta_path, min_length=80):
    """Check if all sequences in the FASTA file are at least `min_length` bases long."""
    for record in SeqIO.parse(fasta_path, "fasta"):
        if len(record.seq) < min_length:
            return False
    return True

def process_gene_exon_alignment(
    gene_name, threads, input_dir, ref_dir, output_phylo, log_file, min_size,
    tree_method="fasttree", iqtree_mode="fixed"
):
    """
    Align and build tree for all exon contigs of a given gene.
    Supports FastTree and IQ-TREE (fixed or MFP mode).
    Args:
        gene_name: Name of the gene to process.
        threads: Number of threads to use for parallel processing.
        input_dir: Directory containing exon FASTA files.
        ref_dir: Directory containing reference alignments.
        output_phylo: Directory to save output trees.
        log_file: Log file to record status messages.
        min_size: Minimum length of exon sequences to include.
        tree_method: Phylogenetic tree construction method ('fasttree' or 'iqtree').
        iqtree_mode: IQ-TREE mode ('fixed', 'fixed+gamma', or 'mfp').
 
        """
    # Find all exon FASTA files for this gene
    try:
        pattern = os.path.join(input_dir, f"*{gene_name}*.fasta")
        exon_files = glob.glob(pattern)
    except Exception as e:
        log_status(log_file, f"List exons for {gene_name}: FAILURE")
        print(f"Error listing exons for {gene_name}: {e}")
        return
    log_status(log_file, f"List exons for {gene_name}: SUCCESS")
    if not exon_files:
        log_status(log_file, f"No exon files found for {gene_name}, skipping.")
        return
    exon_files.sort()
    for i, exon_path in enumerate(exon_files, start=1):
        # Enforce minimum exon length
        if not all_sequences_meet_minimum_length(exon_path, min_size):
            log_status(log_file, f"Skipping {os.path.basename(exon_path)} (sequences < {min_size} bp)")
            continue
        # Alignment with MAFFT
        ref_alignment = os.path.join(ref_dir, f"{gene_name}.fasta")
        aligned_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}_aligned.fasta")
        mafft_cmd = (
            f"mafft --preservecase --maxiterate 1000 --localpair --adjustdirection "
            f"--thread {threads} --addfragments {exon_path} {ref_alignment} > {aligned_out}"
        )
        run_command(mafft_cmd, f"MAFFT alignment for {gene_name} exon {i}", log_file)
        # Trim alignment with trimAl
        trimmed_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}_trimmed.fasta")
        trimal_cmd = f"trimal -in {aligned_out} -out {trimmed_out} -gt 0.5"
        run_command(trimal_cmd, f"Trim alignment for {gene_name} exon {i}", log_file)
        # Build tree with selected method
        tree_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}.tre")
        if tree_method == "fasttree":
            fasttree_cmd = f"fasttree -gtr -gamma -nt {trimmed_out} > {tree_out}"
            run_command(fasttree_cmd, f"Tree construction for {gene_name} exon {i} (FastTree)", log_file)
        elif tree_method == "iqtree":
            if iqtree_mode == "mfp":
                iqtree_cmd = (
                    f"iqtree2 -s {trimmed_out} -nt {threads} -m MFP -pre {tree_out.replace('.tre','')} "
                    f"--quiet"
                )
>>>>>>> theirs
=======
#!/usr/bin/env python3
"""
02_exon_trees.py – Align exon sequences and build gene trees.
This script processes exon sequences for each gene, aligning them to a reference
alignment, trimming the alignment, and constructing a phylogenetic tree.

It supports parallel processing of multiple genes and logs the status of each step.
It requires the following tools:
- MAFFT for sequence alignment
- trimAl for trimming alignments
- FastTree for phylogenetic tree construction
It can be run with command-line arguments or a configuration file.

Arguments:
- -c, --config: Path to configuration file (YAML/JSON/TOML)
- -t, --threads: Number of CPU threads to use for parallel processing
- -e, --input_exon: Directory of extracted exon FASTA files
- -r, --ref_alignment: Directory of reference alignments
- -g, --gene_list: Path to gene list file
- -p, --proj_name: Project name identifier
- -o, --output_phylo: Directory for output trees
- -m, --min_exon_size: Minimum exon length to include in analysis
- --tree_method: Phylogenetic tree construction method (fasttree or iqtree)
- --iqtree_mode: IQ-TREE mode (fixed, fixed+gamma, or mfp)

Usage:
python 02_exon_trees.py -c config.yaml -t 8 -e 02_exon_extracted -r ref -g gene_list.txt -p my_project -o 03_phylo_results -m 80 --tree_method fasttree --iqtree_mode fixed
or
python 02_exon_trees.py --threads 8 --input_exon 02_exon_extracted --ref_alignment ref --gene_list gene_list.txt --proj_name my_project --output_phylo 03_phylo_results --min_exon_size 80 --tree_method fasttree --iqtree_mode fixed
"""
import os
import glob
import argparse
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from pipeline_utils import log_status, run_command, load_config, is_valid_project_name

def all_sequences_meet_minimum_length(fasta_path, min_length=80):
    """Check if all sequences in the FASTA file are at least `min_length` bases long."""
    for record in SeqIO.parse(fasta_path, "fasta"):
        if len(record.seq) < min_length:
            return False
    return True

def process_gene_exon_alignment(
    gene_name, threads, input_dir, ref_dir, output_phylo, log_file, min_size,
    tree_method="fasttree", iqtree_mode="fixed"
):
    """
    Align and build tree for all exon contigs of a given gene.
    Supports FastTree and IQ-TREE (fixed or MFP mode).
    Args:
        gene_name: Name of the gene to process.
        threads: Number of threads to use for parallel processing.
        input_dir: Directory containing exon FASTA files.
        ref_dir: Directory containing reference alignments.
        output_phylo: Directory to save output trees.
        log_file: Log file to record status messages.
        min_size: Minimum length of exon sequences to include.
        tree_method: Phylogenetic tree construction method ('fasttree' or 'iqtree').
        iqtree_mode: IQ-TREE mode ('fixed', 'fixed+gamma', or 'mfp').
 
        """
    # Find all exon FASTA files for this gene
    try:
        pattern = os.path.join(input_dir, f"*{gene_name}*.fasta")
        exon_files = glob.glob(pattern)
    except Exception as e:
        log_status(log_file, f"List exons for {gene_name}: FAILURE")
        print(f"Error listing exons for {gene_name}: {e}")
        return
    log_status(log_file, f"List exons for {gene_name}: SUCCESS")
    if not exon_files:
        log_status(log_file, f"No exon files found for {gene_name}, skipping.")
        return
    exon_files.sort()
    for i, exon_path in enumerate(exon_files, start=1):
        # Enforce minimum exon length
        if not all_sequences_meet_minimum_length(exon_path, min_size):
            log_status(log_file, f"Skipping {os.path.basename(exon_path)} (sequences < {min_size} bp)")
            continue
        # Alignment with MAFFT
        ref_alignment = os.path.join(ref_dir, f"{gene_name}.fasta")
        aligned_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}_aligned.fasta")
        mafft_cmd = (
            f"mafft --preservecase --maxiterate 1000 --localpair --adjustdirection "
            f"--thread {threads} --addfragments {exon_path} {ref_alignment} > {aligned_out}"
        )
        run_command(mafft_cmd, f"MAFFT alignment for {gene_name} exon {i}", log_file)
        # Trim alignment with trimAl
        trimmed_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}_trimmed.fasta")
        trimal_cmd = f"trimal -in {aligned_out} -out {trimmed_out} -gt 0.5"
        run_command(trimal_cmd, f"Trim alignment for {gene_name} exon {i}", log_file)
        # Build tree with selected method
        tree_out = os.path.join(output_phylo, f"{gene_name}_exon_{i}.tre")
        if tree_method == "fasttree":
            fasttree_cmd = f"fasttree -gtr -gamma -nt {trimmed_out} > {tree_out}"
            run_command(fasttree_cmd, f"Tree construction for {gene_name} exon {i} (FastTree)", log_file)
        elif tree_method == "iqtree":
            if iqtree_mode == "mfp":
                iqtree_cmd = (
                    f"iqtree2 -s {trimmed_out} -nt {threads} -m MFP -pre {tree_out.replace('.tre','')} "
                    f"--quiet"
                )
>>>>>>> theirs
            else:  # fixed model (GTR+G or GTR)
                model = "GTR+G" if iqtree_mode == "fixed+gamma" else "GTR"
                iqtree_cmd = (
                    f"iqtree2 -s {trimmed_out} -nt {threads} -m {model} "
                    f"-pre {tree_out.replace('.tre','')} --quiet"
                )
<<<<<<< ours
<<<<<<< ours
            run_command(iqtree_cmd, f"Tree construction for {gene_name} exon {i} (IQ-TREE)", log_file)
            # IQ-TREE outputs .treefile, so rename/move to .tre for consistency
            iqtree_treefile = tree_out.replace('.tre', '.treefile')
            if os.path.exists(iqtree_treefile):
                os.replace(iqtree_treefile, tree_out)
        else:
            log_status(log_file, f"Unknown tree method: {tree_method}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align exons and build exon trees for each gene.")
    parser.add_argument("-c", "--config", help="Path to config file (YAML/JSON/TOML)")
    parser.add_argument("-t", "--threads", type=int, help="Number of CPU threads for parallel processing")
    parser.add_argument("-e", "--input_exon", help="Directory of extracted exon FASTA files", default="02_exon_extracted")
    parser.add_argument("-r", "--ref_alignment", help="Directory of reference alignments", default="ref")
    parser.add_argument("-g", "--gene_list", help="Path to gene list file", default="gene_list.txt")
    parser.add_argument("-p", "--proj_name", help="Project name identifier")
    parser.add_argument("-o", "--output_phylo", help="Directory for output trees", default="03_phylo_results")
    parser.add_argument("-m", "--min_exon_size", type=int, help="Minimum exon length to include", default=80)
=======
            run_command(iqtree_cmd, f"Tree construction for {gene_name} exon {i} (IQ-TREE)", log_file)
            # IQ-TREE outputs .treefile, so rename/move to .tre for consistency
            iqtree_treefile = tree_out.replace('.tre', '.treefile')
            if os.path.exists(iqtree_treefile):
                os.replace(iqtree_treefile, tree_out)
        else:
            log_status(log_file, f"Unknown tree method: {tree_method}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align exons and build exon trees for each gene.")
    parser.add_argument("-c", "--config", help="Path to config file (YAML/JSON/TOML)")
    parser.add_argument("-t", "--threads", type=int, help="Number of CPU threads for parallel processing")
    parser.add_argument("-e", "--input_exon", help="Directory of extracted exon FASTA files", default="02_exon_extracted")
    parser.add_argument("-r", "--ref_alignment", help="Directory of reference alignments", default="ref")
    parser.add_argument("-g", "--gene_list", help="Path to gene list file", default="gene_list.txt")
    parser.add_argument("-p", "--proj_name", help="Project name identifier")
    parser.add_argument("-o", "--output_phylo", help="Directory for output trees", default="03_phylo_results")
    parser.add_argument("-m", "--min_exon_size", type=int, help="Minimum exon length to include", default=80)
>>>>>>> theirs
=======
            run_command(iqtree_cmd, f"Tree construction for {gene_name} exon {i} (IQ-TREE)", log_file)
            # IQ-TREE outputs .treefile, so rename/move to .tre for consistency
            iqtree_treefile = tree_out.replace('.tre', '.treefile')
            if os.path.exists(iqtree_treefile):
                os.replace(iqtree_treefile, tree_out)
        else:
            log_status(log_file, f"Unknown tree method: {tree_method}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align exons and build exon trees for each gene.")
    parser.add_argument("-c", "--config", help="Path to config file (YAML/JSON/TOML)")
    parser.add_argument("-t", "--threads", type=int, help="Number of CPU threads for parallel processing")
    parser.add_argument("-e", "--input_exon", help="Directory of extracted exon FASTA files", default="02_exon_extracted")
    parser.add_argument("-r", "--ref_alignment", help="Directory of reference alignments", default="ref")
    parser.add_argument("-g", "--gene_list", help="Path to gene list file", default="gene_list.txt")
    parser.add_argument("-p", "--proj_name", help="Project name identifier")
    parser.add_argument("-o", "--output_phylo", help="Directory for output trees", default="03_phylo_results")
    parser.add_argument("-m", "--min_exon_size", type=int, help="Minimum exon length to include", default=80)
>>>>>>> theirs
    parser.add_argument("--tree_method", choices=["fasttree", "iqtree"], default=None,
                        help="Phylogeny method: fasttree or iqtree (default: fasttree)")
    parser.add_argument("--iqtree_mode", choices=["fixed", "fixed+gamma", "mfp"], default=None,
                        help="IQ-TREE mode: fixed (GTR), fixed+gamma (GTR+G), or mfp (ModelFinder Plus)")
<<<<<<< ours
<<<<<<< ours
    args = parser.parse_args()

    # Load config if provided
    config = {}
    if args.config:
        config = load_config(args.config)
    # Determine parameters (CLI overrides config)
    threads = args.threads if args.threads is not None else config.get('threads')
    proj_name = args.proj_name or config.get('proj_name')
    input_exon_dir = args.input_exon if args.input_exon != parser.get_default('input_exon') else config.get('input_exon', "02_exon_extracted")
    ref_dir = args.ref_alignment if args.ref_alignment != parser.get_default('ref_alignment') else config.get('ref_alignment', "ref")
    gene_list_path = args.gene_list if args.gene_list != parser.get_default('gene_list') else config.get('gene_list', "gene_list.txt")
    output_phylo = args.output_phylo if args.output_phylo != parser.get_default('output_phylo') else config.get('output_phylo', "03_phylo_results")
    min_size = args.min_exon_size if args.min_exon_size != parser.get_default('min_exon_size') else config.get('min_exon_size', 80)
    tree_method = args.tree_method or config.get('tree_method', 'fasttree')
    iqtree_mode = args.iqtree_mode or config.get('iqtree_mode', 'fixed')
    # Disable iqtree_mode if tree_method is fasttree
    if tree_method == "fasttree":
        iqtree_mode = None
    if threads is None or not proj_name:
        parser.error("Required parameters missing: threads and proj_name must be specified.")
    if not is_valid_project_name(proj_name):
        parser.error(f"Project name '{proj_name}' contains invalid characters.")
    threads = int(threads)
    min_size = int(min_size)
    # Initialize log file
    log_file = f"{proj_name}_02_exons_phylo.log"
    if os.path.exists(log_file):
        os.remove(log_file)
    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {threads}")
    log_status(log_file, f"  Input Exon Directory: {input_exon_dir}")
    log_status(log_file, f"  Reference Alignment Directory: {ref_dir}")
    log_status(log_file, f"  Gene List: {gene_list_path}")
    log_status(log_file, f"  Project Name: {proj_name}")
    log_status(log_file, f"  Output Directory: {output_phylo}")
    log_status(log_file, f"  Minimum Exon Size: {min_size}")
    log_status(log_file, f"  Tree Method: {tree_method}")
    if tree_method == "iqtree":
        log_status(log_file, f"  IQ-TREE Mode: {iqtree_mode}")
    os.makedirs(output_phylo, exist_ok=True)
    log_status(log_file, f"Created directory {output_phylo}")
    # Read gene list and execute alignments/trees in parallel
    with open(gene_list_path, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                process_gene_exon_alignment,
                gene, threads, input_exon_dir, ref_dir, output_phylo, log_file, min_size,
                tree_method, iqtree_mode
            )
            for gene in genes
        ]
        for future in futures:
            future.result()
    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for details.")
>>>>>>> theirs
=======
    args = parser.parse_args()

    # Load config if provided
    config = {}
    if args.config:
        config = load_config(args.config)
    # Determine parameters (CLI overrides config)
    threads = args.threads if args.threads is not None else config.get('threads')
    proj_name = args.proj_name or config.get('proj_name')
    input_exon_dir = args.input_exon if args.input_exon != parser.get_default('input_exon') else config.get('input_exon', "02_exon_extracted")
    ref_dir = args.ref_alignment if args.ref_alignment != parser.get_default('ref_alignment') else config.get('ref_alignment', "ref")
    gene_list_path = args.gene_list if args.gene_list != parser.get_default('gene_list') else config.get('gene_list', "gene_list.txt")
    output_phylo = args.output_phylo if args.output_phylo != parser.get_default('output_phylo') else config.get('output_phylo', "03_phylo_results")
    min_size = args.min_exon_size if args.min_exon_size != parser.get_default('min_exon_size') else config.get('min_exon_size', 80)
    tree_method = args.tree_method or config.get('tree_method', 'fasttree')
    iqtree_mode = args.iqtree_mode or config.get('iqtree_mode', 'fixed')
    # Disable iqtree_mode if tree_method is fasttree
    if tree_method == "fasttree":
        iqtree_mode = None
    if threads is None or not proj_name:
        parser.error("Required parameters missing: threads and proj_name must be specified.")
    if not is_valid_project_name(proj_name):
        parser.error(f"Project name '{proj_name}' contains invalid characters.")
    threads = int(threads)
    min_size = int(min_size)
    # Initialize log file
    log_file = f"{proj_name}_02_exons_phylo.log"
    if os.path.exists(log_file):
        os.remove(log_file)
    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {threads}")
    log_status(log_file, f"  Input Exon Directory: {input_exon_dir}")
    log_status(log_file, f"  Reference Alignment Directory: {ref_dir}")
    log_status(log_file, f"  Gene List: {gene_list_path}")
    log_status(log_file, f"  Project Name: {proj_name}")
    log_status(log_file, f"  Output Directory: {output_phylo}")
    log_status(log_file, f"  Minimum Exon Size: {min_size}")
    log_status(log_file, f"  Tree Method: {tree_method}")
    if tree_method == "iqtree":
        log_status(log_file, f"  IQ-TREE Mode: {iqtree_mode}")
    os.makedirs(output_phylo, exist_ok=True)
    log_status(log_file, f"Created directory {output_phylo}")
    # Read gene list and execute alignments/trees in parallel
    with open(gene_list_path, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                process_gene_exon_alignment,
                gene, threads, input_exon_dir, ref_dir, output_phylo, log_file, min_size,
                tree_method, iqtree_mode
            )
            for gene in genes
        ]
        for future in futures:
            future.result()
    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for details.")
>>>>>>> theirs
=======
    args = parser.parse_args()

    # Load config if provided
    config = {}
    if args.config:
        config = load_config(args.config)
    # Determine parameters (CLI overrides config)
    threads = args.threads if args.threads is not None else config.get('threads')
    proj_name = args.proj_name or config.get('proj_name')
    input_exon_dir = args.input_exon if args.input_exon != parser.get_default('input_exon') else config.get('input_exon', "02_exon_extracted")
    ref_dir = args.ref_alignment if args.ref_alignment != parser.get_default('ref_alignment') else config.get('ref_alignment', "ref")
    gene_list_path = args.gene_list if args.gene_list != parser.get_default('gene_list') else config.get('gene_list', "gene_list.txt")
    output_phylo = args.output_phylo if args.output_phylo != parser.get_default('output_phylo') else config.get('output_phylo', "03_phylo_results")
    min_size = args.min_exon_size if args.min_exon_size != parser.get_default('min_exon_size') else config.get('min_exon_size', 80)
    tree_method = args.tree_method or config.get('tree_method', 'fasttree')
    iqtree_mode = args.iqtree_mode or config.get('iqtree_mode', 'fixed')
    # Disable iqtree_mode if tree_method is fasttree
    if tree_method == "fasttree":
        iqtree_mode = None
    if threads is None or not proj_name:
        parser.error("Required parameters missing: threads and proj_name must be specified.")
    if not is_valid_project_name(proj_name):
        parser.error(f"Project name '{proj_name}' contains invalid characters.")
    threads = int(threads)
    min_size = int(min_size)
    # Initialize log file
    log_file = f"{proj_name}_02_exons_phylo.log"
    if os.path.exists(log_file):
        os.remove(log_file)
    log_status(log_file, "Pipeline started with the following parameters:")
    log_status(log_file, f"  Threads: {threads}")
    log_status(log_file, f"  Input Exon Directory: {input_exon_dir}")
    log_status(log_file, f"  Reference Alignment Directory: {ref_dir}")
    log_status(log_file, f"  Gene List: {gene_list_path}")
    log_status(log_file, f"  Project Name: {proj_name}")
    log_status(log_file, f"  Output Directory: {output_phylo}")
    log_status(log_file, f"  Minimum Exon Size: {min_size}")
    log_status(log_file, f"  Tree Method: {tree_method}")
    if tree_method == "iqtree":
        log_status(log_file, f"  IQ-TREE Mode: {iqtree_mode}")
    os.makedirs(output_phylo, exist_ok=True)
    log_status(log_file, f"Created directory {output_phylo}")
    # Read gene list and execute alignments/trees in parallel
    with open(gene_list_path, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                process_gene_exon_alignment,
                gene, threads, input_exon_dir, ref_dir, output_phylo, log_file, min_size,
                tree_method, iqtree_mode
            )
            for gene in genes
        ]
        for future in futures:
            future.result()
    log_status(log_file, "Pipeline completed successfully.")
    print(f"Pipeline completed. Check {log_file} for details.")
>>>>>>> theirs
