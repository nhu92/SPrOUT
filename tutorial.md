# Pipeline Tutorial for Plant Species Identification from Mixed Samples

Nan Hu, 2024 August

---

This tutorial provides step-by-step instructions to run a pipeline designed for identifying plant species from mixed DNA samples using the Angiosperms353 target sequencing kit and the HybPiper workflow. The pipeline is efficient and cost-effective, making it a valuable tool in various scientific and practical applications.

## Prerequisites

Before running the pipeline, ensure that you have the following:

- **Software and Tools**:
  - Python 3.9+
  - `HybPiper 2.2.0`, `fastp 0.23.4`, `mafft 7.526`, `fasttree 2.1.11`, `seqkit 2.8.2`, `trimal 1.5.0`. Suggest using Conda to install
  - Required Python libraries: `pandas 2.2.2`, `argparse 1.4.0`, `scipy 1.14.0`, `scikit-learn 1.5.1`, `numpy 2.0.1`, `biopython 1.84`

```bash
  # Create environment for HybPiper
  conda config --add channels defaults
  conda config --add channels bioconda
  conda config --add channels conda-forge
  conda create -n pacificod hybpiper
  conda activate pacificod
  # Install dependencies
  conda install seqkit fasttree fastp trimal
  pip install numpy pandas scipy scikit-learn biopython
```

- **Data**:
  - Paired-end reads from mixed plant DNA samples.
  - Reference database of Angiosperms353 sequences.
  > Reference sequences name should follow Order_Family_Genus_Species format for prediction. For example: >Rosales_Rosaceae_Rose_rosa

### Quick Example Run

Run this [Preparing commands](https://github.com/nhu92/Mixed-DNA-Project/blob/main/release/v1.1_underconstruction/test_run.sh) to clone the package and download the sample input files.

To run the entire pipeline, execute the script [here](https://github.com/nhu92/Mixed-DNA-Project/blob/main/release/v1.1_underconstruction/sample_command.sh) (Using SLURM job submission system as an example). These steps require running on a node/job submission systems, preferring high performance computer clusters. You need to modify the job names and replace the placeholder "RENAME" with your sample names. Then, submit the script with:
```bash
sbatch sample_command.sh
```

## Background

Plant species identification from mixed samples is critical in fields such as ecology, conservation, and agriculture. Traditional methods face challenges due to high costs, computational inefficiencies, and incomplete databases. This pipeline addresses these issues by leveraging a novel approach that combines target sequencing and phylogenetic inference.

## 1. Data Processing

- **Input**: The process starts with raw mixed plant samples.
- **Steps**:
  - **DNA Extraction**: DNA is extracted from the mixed plant sample following established protocols (e.g., Hale et al., 2019).
  - **Target Sequencing Library Preparation**: The extracted DNA undergoes library preparation for target sequencing (as Slimp et al., 2021).
  - **Illumina Sequencing**: The prepared libraries are sequenced using Illumina technology to produce short reads, which serve as the input for the next stage.

## 2. Target Assembly

- **Input**: Short reads generated from the Illumina sequencing step.
- **Steps**:
  - **FastP Trimming & QC**: The short reads undergo quality control and trimming using FastP to produce clean reads.
  - **HybPiper**: The clean reads are assembled using HybPiper, generating BAM files, contigs in FASTA format, and exons in FASTA format. These outputs are critical for the subsequent phylogenetic inference.

## 3. Phylogenetic Inference

- **Input**: Contigs and exon FASTA files from the Target Assembly stage.
- **Steps**:
  - **MAFFT Alignment**: The contigs and exons are aligned using MAFFT against a reference panel.
  - **Trimal**: The aligned exons are trimmed to ensure consistency and remove low-quality alignments.
  - **FastTree/IQ-TREE**: The trimmed alignments are used to infer phylogenetic trees, which illustrate the evolutionary relationships among the sequences.

## 4. Prediction

- **Input**: Exon trees generated during Phylogenetic Inference.
- **Steps**:
  - **Genetic Distance Calculation**: Distance matrices are calculated from the phylogenetic trees to quantify the genetic differences among the taxa.
  - **ACS Evaluation**: An ACS (Adaptive Cluster Sampling) array is created based on the distance matrices to evaluate the species' relationships.
  - **Species Prediction**: The final species prediction is made using the ACS array, identifying the plant species present in the mixed sample.

## Pipeline Overview

The pipeline consists of four main steps:

1. **Sequence Assembly**
2. **Exon Tree Creation**
3. **Distance Matrix Calculation**
4. **Prediction and Identification**

### Step 1: Sequence Assembly (`01_exons_assembly.py`)

This step involves trimming and assembling short paired-end reads using `fastp` and `HybPiper`.

#### Usage

```bash
python 01_exons_assembly.py -t <threads> -r1 <read1.fastq> -r2 <read2.fastq> -m <mega353.fasta> -p <project_name> -g <gene_list>
```

**Arguments**:

- `-t` or `--threads`: Number of CPU threads to use.
- `-r1` or `--read1`: Path to the first read file (FASTQ format).
- `-r2` or `--read2`: Path to the second read file (FASTQ format).
- `-m` or `--mega353`: Path to the MEGA353 reference file (FASTA format). Default is `angiosperms353_v2_interim_targetfile.fasta`.
- `-p` or `--project_name`: Project name for output files.
- `-g` or `--gene_list`: Path to the selected gene list.
- `-ov` or `--overlapping_rate`: Overlapping to consider the same exon (0-1, default 0.8).
- `--output_hyb`: Output directory for HybPiper results (default is '01_hyb_output').
- `--output_exon`: Output directory for extracted exons (default is '02_exon_extracted').

### Step 2: Exon Tree Creation (`02_exon_trees.py`)

In this step, the pipeline extracts exons from the assembled sequences and builds exon trees for phylogenetic inference.

#### Usage

```bash
python 02_exon_trees.py -t <threads> -e <exon_dir> -r <ref_dir> -p <project_name> -g <gene_list>
```

**Arguments**:

- `-t` or `--threads`: Number of CPU threads to use.
- `-e` or `--input_exon`: Directory of extracted exon sequences. Default is `02_exon_extracted`.
- `-r` or `--ref_alignment`: Directory of reference alignments. Default is `ref`.
- `-p` or `--project_name`: Project name for output files.
- `-g` or `--gene_list`: Path to the list of gene names. Default is generated from Step 1 `gene_list.txt`.
- `--exon_min_size`: Minimum exon size to include in analysis (default is 80bp)
- `--output_dir`: Output directory for phylogenetic trees (default is 03_phylo_results)

### Step 3: Distance Matrix Calculation (`03_distance_matrices.py`)

This step calculates distance matrices based on the exon trees created in the previous step. These matrices are used in phylogenetic analysis to infer relationships among the taxa in the mixed sample.

#### Usage

```bash
python 03_distance_matrices.py -t <threads> -p <project_name> -g <gene_list> --threshold <float_num> [--use_flag] 
```

**Arguments**:

- `-t` or `--threads`: Number of CPU threads to use.
- `-p` or `--project_name`: Project name for output files.
- `-g` or `--gene_list`: Path to the list of gene names. Default is generated from Step 1 `gene_list.txt`.
- `--threshold`: Parameter to control specificity of candidate branches on phylogeny. Default is `1.96`, which means only genetic similarity is over `1.96*SD+Mean` will be counted.
- `--use_flag`: A boolean parameter to mark if all the non-outliers' similarity set to 0.
- `--input_dir`: Input directory containing the tree files (default is 03_phylo_results)
- `--output_dir`: Output directory for storing results (default is 04_all_trees)

### Step 4: Prediction and Identification (`04_prediction.py`)

The final step involves predicting and identifying the plant species present in the mixed sample by applying phylogenetic inference on the computed distance matrices.

#### Usage

```bash
python 04_prediction.py -i <input_ACS_file> -o <output_file> -tl <taxonomic_level>
```

**Arguments**:

- `-i` or `--input_distance_matrix`: Path to the input cumulative distances file (*cumulative_dist.csv).
- `-o` or `--output_file`: Path to save the prediction results.
- `-tl` or `--taxonomic_level`: Taxonomic level to process (choices: `o` for Order, `f` for Family, `g` for Genus, `s` for Species).
- `-z` or `--zscore_threshold`: Z-score threshold to select taxonomy names.
- `-to` or `--taxonomy_output_file`: Path to the output file for selected taxonomy names.
