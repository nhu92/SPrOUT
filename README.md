# SPrOUT

Nan Hu, 2025 August

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

## Preparing a GUI-ready results bundle

After completing the pipeline, you can package outputs for a downstream GUI (e.g., R Shiny or Dash) without building the GUI itself. Run:

```bash
python 05_gui_ready.py -p <project_name> \
  --exon_dir 02_exon_extracted --tree_dir 03_phylo_results --matrix_dir 04_all_trees
```

This creates `gui_ready/<project>_gui_bundle/` containing:
- `trees/merged_exon_trees.nwk` and `trees/tree_index.csv` for quick browsing.
- `exon_metrics/exon_metrics.csv` with per-exon length and coverage summaries.
- `matrices/` with similarity tables, cumulative/summary distances, and tree contribution metrics.
- `summaries/run_summary.json` plus a manifest listing all packaged files.
- A compressed archive `sprout_<project>_results.zip` for download or transfer to the GUI repository.
