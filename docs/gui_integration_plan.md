# SPrOUT Overview and GUI Integration Plan

## Repository Overview
- **Purpose**: SPrOUT is a computational pipeline for identifying plant species from mixed DNA using target sequencing (e.g., Angiosperms353). It automates read processing, exon extraction, phylogenetic tree construction, distance summarization, and prediction scoring.
- **Core pipeline stages**:
  1. **Exon assembly and extraction** (`01_exons_assembly.py`): trims reads, runs HybPiper assembly, and extracts exon contigs per gene with overlap-aware naming and per-gene logs.
  2. **Exon alignment and tree building** (`02_exon_trees.py`): aligns exon contigs to reference alignments, trims them, and builds exon-level gene trees via FastTree or IQ-TREE with parallel execution and minimum-length filters.
  3. **Distance matrices** (`03_distance_matrices.py`): roots trees, computes pairwise genetic distances, filters matrices, converts to similarities, and aggregates similarity totals per taxon for downstream scoring; records sister taxa for collapsed nodes.
  4. **Prediction aggregation** (`04_prediction.py`): combines similarity totals with thresholds/flags to report candidate species predictions and summary tables.
- **Shared utilities**: `pipeline_utils.py` provides logging, shell command helpers, config loading, and project-name validation.
- **Configuration and examples**: `config.yaml`, `sample_command.sh`, `test_run.sh`, and `sample_data/` illustrate end-to-end execution.

## Goals for a GUI-Ready v1.2
Create a companion GUI (e.g., R Shiny, Dash) to visualize SPrOUT outputs and streamline downloads while keeping the CLI pipeline intact.

### Packaging and Download Experience
- Produce a **single compressed archive** per run (e.g., `sprout_<project>_results.zip`) bundling logs, exon FASTAs, trees, distance matrices, and predictions, with a manifest and README for the GUI.
- Include a compact **run summary JSON** (inputs, versions, gene counts, timing) for quick GUI ingestion.

### Visualization-Ready Exon Trees
- **Merge exon trees into one file** (multi-tree Newick/NeXML) so users can load all exon trees into external viewers; keep per-exon identifiers for filtering.
- Add a **tree index table** (tree file, gene, exon label, sequence count, min/median/ max branch length) to drive GUI navigation.

### Feature Metrics for Users
- Expose **per-exon statistics**: assembled length distribution, alignment trimming retention, mapping coverage, and overlap notes from `exonerate_stats` parsing.
- Provide **run-level QC**: read trimming stats (fastp), total contigs per gene, failed/filtered exons, and reference alignment coverage.

### Prediction Transparency
- Emit a **contribution table** mapping each exon tree to its weight toward the final similarity or ACS-style score (e.g., percentage of total similarity, count of votes after thresholding), and surface the threshold used.
- Allow a **GUI slider to reapply thresholds** by saving intermediate similarity matrices and a small recomputation script/notebook the GUI can call.

### Data Model for the GUI Repository
- Input: zipped run artifact + summary JSON.
- Backend helpers: lightweight APIs to unpack artifacts, load tree bundle, compute thresholded predictions, and summarize contributions.
- Frontend views: run overview, exon-tree gallery with search/filter, interactive threshold slider with live prediction table, and QC/download section.

### Additional Output Improvements
- Normalize file naming and directory structure (exons → alignments → trees → matrices → predictions) and capture it in the manifest.
- Bundle **co-occurrence and sister-taxa notes** (from collapsed nodes) alongside distance matrices to aid interpretation.
- Provide **export shortcuts**: combined TSV/CSV tables for exon stats, tree index, contribution metrics, and final predictions for downstream sharing.

### Implementation Roadmap (CLI side)
1. Add a post-processing step to build the tree bundle, contribution table, and summary JSON.
2. Enhance distance/prediction steps to emit per-exon metrics and intermediate similarity tables for threshold replays.
3. Add a packaging script to assemble the standardized archive and manifest for GUI consumption.
