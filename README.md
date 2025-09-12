# SPrOUT

Nan Hu, 2025 August

---

This tutorial provides step-by-step instructions to run a pipeline designed for identifying plant species from mixed DNA samples using the Angiosperms353 target sequencing kit and the HybPiper workflow. The pipeline is efficient and cost-effective, making it a valuable tool in various scientific and practical applications.


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
  conda create -n sprout hybpiper
  conda activate sprout
  # Install dependencies
  conda install seqkit fasttree fastp trimal
  pip install numpy pandas scipy scikit-learn biopython
```

## Installation

```bash
  mkdir SPrOut_test
  git clone https://github.com/nhu92/SPrOut.git
  cd SPrOut_test
```

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
  conda create -n sprout hybpiper
  conda activate sprout
  # Install dependencies
  conda install seqkit fasttree fastp trimal
  pip install numpy pandas scipy scikit-learn biopython
```

- **Data**:
  - Paired-end reads from mixed plant DNA samples.
  - Reference database of Angiosperms353 sequences.
  > Reference sequences name should follow Order_Family_Genus_Species format for prediction. For example: >Rosales_Rosaceae_Rose_rosa

### Quick Example Run

Run this [Preparing commands](https://github.com/nhu92/SPrOUT/blob/main/test_run.sh) to clone the package and download the sample input files.

To run the entire pipeline, execute the script [here](https://github.com/nhu92/SPrOUT/blob/main/sample_command.sh) (Using SLURM job submission system as an example). These steps require running on a node/job submission systems, preferring high performance computer clusters. You need to modify the job names and replace the placeholder "RENAME" with your sample names. Then, submit the script with:
```bash
sbatch sample_command.sh
```

Detailed instructions please refer to [Wiki page](https://github.com/nhu92/SPrOUT/wiki)
