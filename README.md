# SPrOUT

Nan Hu, 2025 November

---
For citation, please refer to this:

Nan Hu, Madison Bullock, Chris Jackson, Courtney Miller, Elizabeth Sage Hunter, Charles Huff, Yanni Chen, Sara M. Handy, Matthew G. Johnson. SPrOUT: A computational and targeted sequencing approach for mixed plant DNA identification with Angiosperms353. *In Review*.

---
SPrOUT is a computational pipeline designed for predicting species taxonomic information from single and mixed plant samples, using target sequencing arrays. We provides step-by-step instructions to run a pipeline designed for identifying plant species from mixed DNA samples using the Angiosperms353 target sequencing kit and the HybPiper workflow. SPrOUT is designed to be computational efficient and lab cost-effective for complicated mixed species identification with nuclear gene involved.

## Prerequisites

Before running the pipeline, ensure that you have the following:

- **Software and Tools**:
  - Python 3.11+
  - `HybPiper 2.2.0`, `fastp 0.23.4`, `mafft 7.526`, `fasttree 2.1.11`, `seqkit 2.8.2`, `trimal 1.5.0`. Suggest using Conda to install
  - Required Python libraries: `pandas 2.2.2`, `argparse 1.4.0`, `scipy 1.14.0`, `scikit-learn 1.5.1`, `numpy 2.0.1`, `biopython 1.84`

```bash
  # Create environment for HybPiper
  conda config --add channels defaults
  conda config --add channels bioconda
  conda config --add channels conda-forge
  conda create -n sprout hybpiper -y
  conda activate sprout
  # Install dependencies
  conda install -y seqkit fasttree fastp trimal
  pip install numpy pandas scipy scikit-learn biopython
```

- **Data**:
  - Paired-end reads from mixed plant DNA samples.
  - Reference database of Angiosperms353 sequences. (Actually, these can be any target sequencing data but requires similar  well curated gene alignments as references)
  > Reference sequences name should follow Order_Family_Genus_Species format for prediction. For example: >Rosales_Rosaceae_Rose_rosa

## Installation

```bash
  mkdir SPrOUT
  git clone https://github.com/nhu92/SPrOUT.git
  cd SPrOUT
```

### Quick Example Run

Run this [Preparing commands](https://github.com/nhu92/SPrOUT/blob/main/test_run.sh) to clone the package and download the sample input files.

To run the entire pipeline, execute the script [here](https://github.com/nhu92/SPrOUT/blob/main/sample_command.sh) (Using SLURM job submission system as an example). These steps require running on a node/job submission systems, preferring high performance computer clusters. You need to modify the input arguments if you aim to test your own sample. After then, submit the script with:
```bash
sbatch sample_command.sh
```

Detailed instructions please refer to [Wiki page](https://github.com/nhu92/SPrOUT/wiki)

<<<<<<< ours
<<<<<<< ours
<<<<<<< ours
<<<<<<< ours
<<<<<<< ours
<<<<<<< ours
## Stage 4 Output Bundle Contract

Starting with this release, Stage 4 automatically packages the prediction artefacts into a portable bundle that downstream tooling (including the upcoming GUI) can consume without re-discovering file locations.

### CLI options

`04_prediction.py` accepts additional arguments to control the bundle:

- `-p` / `--project_name` – Optional project identifier recorded in the manifest.
- `--bundle-format` – Either `directory` (default) or `zip`.
- `--bundle-output` – Custom destination for the generated bundle.
- `--skip-bundle` – Disable bundling if you only need the raw CSV/TXT outputs.
- `--bundle-overwrite` – Replace an existing bundle at the target path.

### Bundle layout

By default, a run that processes project `01x02x03` produces `01x02x03_bundle/` with the following layout (the same hierarchy is zipped when `--bundle-format zip` is chosen):

```
01x02x03_bundle/
├── manifest.json
├── inputs/
│   └── 01x02x03.cumulative_dist.csv
└── reports/
    ├── 01x02x03.predictions.csv
    └── 01x02x03.order_candidates.txt
```

- `manifest.json` contains metadata (project name, taxonomic level, z-score threshold, script version) and SHA-256 checksums for every payload file.
- `inputs/` stores the cumulative distance matrix consumed by Stage 4.
- `reports/` stores the summarized predictions and the taxonomy candidate list.

This contract is shared with the Wiki so automated services can rely on the structure when ingesting completed Stage 4 runs.
=======
=======
>>>>>>> theirs
=======
>>>>>>> theirs
=======
>>>>>>> theirs
=======
>>>>>>> theirs
=======
>>>>>>> theirs
## Repository overview

SPrOUT is organized as a four-step command-line pipeline plus a GUI-export helper:

1. `01_exons_assembly.py` trims reads, runs HybPiper, parses exonerate output, and writes per-exon FASTA files plus exon split tables.
2. `02_exon_trees.py` aligns exon FASTAs to reference alignments, trims alignments, and builds one Newick tree per exon with FastTree or IQ-TREE.
3. `03_distance_matrices.py` converts exon trees to per-tree distance matrices, applies optional filtering, transforms distances to ACS-style similarity scores, and writes cumulative taxon scores.
4. `04_prediction.py` summarizes cumulative scores at order, family, genus, or species level and applies a z-score threshold to produce final predicted taxa.
5. `05_gui_outputs.py` packages the outputs above into files intended for a separate visualization repository.

Shared helpers live in `pipeline_utils.py`; example commands and small test data live in `sample_command.sh`, `test_run.sh`, and `sample_data/`.

## GUI-ready outputs for SPrOUT 1.2

To support a future R Shiny app or equivalent GUI, run the GUI export step after generating distance matrices and predictions:

```bash
python 05_gui_outputs.py \
  --proj_name my_project \
  --input_exon 02_exon_extracted \
  --input_phylo 03_phylo_results \
  --matrix_dir 04_all_trees \
  --prediction_summary my_project.summary_scores.csv \
<<<<<<< ours
<<<<<<< ours
=======
  --taxonomy_output_file selected_taxa.txt \
  --zscore_threshold 1.96 \
  --threshold 1.96 \
>>>>>>> theirs
=======
  --taxonomy_output_file selected_taxa.txt \
  --zscore_threshold 1.96 \
  --threshold 1.96 \
>>>>>>> theirs
  --output_dir 05_gui_results
```

The exporter creates:

- `run_metadata.json`: project name, input locations, creation time, and counts for exon FASTAs, trees, split tables, and contribution rows.
<<<<<<< ours
<<<<<<< ours
- `exon_metrics.csv`: per-exon sequence count, minimum/mean/maximum exon length, total bases, and FASTA path. These fields are intended for basic run-quality cards in a GUI.
- `tree_inventory.csv`: one row per exon tree, with stable `tree_id`, `gene`, `exon_index`, source path, and Newick length for a tree-selector panel.
- `all_exon_trees.nwk`: all exon trees concatenated into one labeled Newick file for review in external tree viewers.
- `tree_contributions.csv`: per-tree, per-taxon ACS contribution values and whether the taxon appears in the thresholded final result. This table is intended to drive contribution bar plots and per-tree drill-downs.
=======
=======
>>>>>>> theirs
- `exon_metrics.csv`: per-exon sequence count, minimum/mean/maximum exon length, total bases, mapped bases, target bases, mapping coverage, alignment-hit count, and FASTA path. These fields are intended for basic run-quality cards in a GUI.
- `tree_inventory.csv`: one row per exon tree, with stable `tree_id`, `gene`, `exon_index`, source path, and Newick length for a tree-selector panel.
- `all_exon_trees.nwk`: all exon trees concatenated into one labeled Newick file for review in external tree viewers.
- `tree_contributions.csv`: per-tree, per-taxon ACS contribution values, fraction of each tree ACS, and whether the taxon appears in the thresholded final result. This table is intended to drive contribution bar plots and per-tree drill-downs.
<<<<<<< ours
>>>>>>> theirs
=======
>>>>>>> theirs
- `result_manifest.json`: a small machine-readable index of the GUI bundle.
- `<project>.sprout_results.zip`: a compressed archive for sharing or upload into a separate GUI repository.

Suggested GUI views for the separate application:

- **Run overview**: read `run_metadata.json` and show total exon FASTAs, exon trees, and contribution rows.
- **Exon QC table**: read `exon_metrics.csv` to display exon length distributions and sequence counts; highlight unusually short or sparse exons.
- **Exon tree browser**: read `tree_inventory.csv` and `all_exon_trees.nwk`; let users select a gene/exon and render the corresponding Newick tree.
- **Prediction threshold explorer**: read the final prediction CSV from `04_prediction.py`; use a slider for z-score threshold and update the selected taxa table interactively.
- **ACS contribution view**: read `tree_contributions.csv`; aggregate `acs_contribution` by taxon, gene, or exon tree to explain how individual trees support final predictions.
- **Download panel**: provide the ZIP bundle and selected CSV/Newick files as direct downloads.
<<<<<<< ours
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
=======
>>>>>>> theirs
