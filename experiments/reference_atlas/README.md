# Learned reference atlas experiment

This branch revisits the reference-panel / PCA / PCoA idea in
[Mixed-DNA-Project](https://github.com/nhu92/Mixed-DNA-Project), with a runnable,
species-disjoint benchmark and actual image-pattern inputs. It is separate from
SPrOUT steps 01–04. See [the research design in Chinese](DESIGN.zh-CN.md) and
[the recorded results](results/REPORT.zh-CN.md).

The experiment asks whether full query-to-reference distance profiles retain
predictive information that a two-dimensional plot loses, and whether a spatial
encoding can preserve that information. A regularized linear classifier on image
pixels provides a first pattern-learning baseline. No CNN or pretrained AI model
is implemented in this version.

## Run from the repository root

Use Python 3.11+ in the SPrOUT environment. Dependencies are isolated to this
experiment. `requirements.txt` gives supported ranges; `requirements-tested.txt`
records the direct versions used for the checked-in results.

```bash
python -m pip install -r experiments/reference_atlas/requirements-tested.txt

python -m experiments.reference_atlas.atlas prepare \
  --alignments sample_data/family298_refs_50genes \
  --output experiments/reference_atlas/runs/family298_distances.npz

python -m experiments.reference_atlas.atlas benchmark \
  --dataset experiments/reference_atlas/runs/family298_distances.npz \
  --rank order --seed 42 --panel-size 48 --grid-size 12 \
  --output experiments/reference_atlas/runs/seed42

# Repeat with --seed 43 and --seed 44, changing the output directory.
python -m pytest tests/test_reference_atlas.py -q
```

The commands use shell line continuations for Bash. In PowerShell, put each
command on one line or replace the trailing backslash with a backtick.

## Use MrKnowItAll for a larger family training library

The user-supplied Angiosperms353 interim ZIP and MrKnowItAll's versioned installer
provide a corrected, denser reference library. This adapter calls the builder's
actual installer and correction code without changing its checkout. The recorded
run used MrKnowItAll commit `c3d7df29aed81071e18b32401fd21b1d0bd541be`.

```bash
git clone https://github.com/nhu92/MrKnowItAll.git ../MrKnowItAll-reference
git -C ../MrKnowItAll-reference checkout c3d7df29aed81071e18b32401fd21b1d0bd541be

python -m experiments.reference_atlas.atlas prepare-refbuilder \
  --refbuilder-source ../MrKnowItAll-reference \
  --archive /path/to/angiosperms_353_v2_interim_targetfile_gene_alignments.zip \
  --backbone-dir experiments/reference_atlas/runs/refbuilder_backbone \
  --gene-list sample_data/gene.list.txt \
  --output experiments/reference_atlas/runs/refbuilder_distances.npz

python -m experiments.reference_atlas.atlas benchmark \
  --dataset experiments/reference_atlas/runs/refbuilder_distances.npz \
  --rank family --panel-size 64 --atlas-missing shortest-path --seed 42 \
  --output experiments/reference_atlas/runs/refbuilder_family_seed42

python -m experiments.reference_atlas.atlas benchmark \
  --dataset experiments/reference_atlas/runs/refbuilder_distances.npz \
  --rank family --within-order Lamiales --panel-size 8 \
  --atlas-missing shortest-path --seed 42 \
  --output experiments/reference_atlas/runs/refbuilder_lamiales_seed42
```

Repeat both benchmarks for seeds 43 and 44. `--reuse-installed` can reuse a prior
installation after matching archive/correction hashes. Import records all source
checksums and installer metadata. The corrected resource has 797 named identity
groups, 59 orders and 291 families; 47 families have sufficient species for the
four-partition benchmark. Headers without a species identity are excluded with
an explicit audit list. In particular, a genus-only `Hemerocallis_1KP-4691` header
must not become a different pseudo-species for each gene suffix. Species placeholders
such as `Genus_sp` are conservatively grouped together; scientific species-resolution
claims require further specimen/name curation.

For an existing MrKnowItAll Kew panel, use its aligned `ref/` directory with
`prepare --alignments ... --gene-list ...`. Build a sufficiently dense **training
library** before selecting a small deployment panel. Two references per family
cannot support four independent species partitions.

`prepare` reads aligned nucleotide `*.fasta` files. It calculates pairwise-deletion
p-distances on A/C/G/T sites, requiring 80 shared valid sites by default, and
records pairwise site overlap. Gaps and ambiguous bases are missing observations.
Headers follow SPrOUT's `Order_Family_Genus_species_source` convention. Different
accessions collapse to one species per locus, choosing the sequence with most
valid bases. Header taxonomy requires curation before scientific use; embedded
underscores within a taxonomic name are unsupported.

The input FASTAs in the repository are **gene alignments**. This first experiment
does not infer exon boundaries. `--max-loci 5` is available for a smoke run, but
those results should not be compared to the recorded 50-gene results.

## What is fitted

For each taxon with at least four species, all accessions/loci of a species are
assigned together to exactly one partition:

1. Anchor candidates, roughly one quarter of species.
2. Training query species, all remaining species after validation/test allocation.
3. Validation query species, at least one per class.
4. Test query species, at least one per class.

Taxa with fewer than four species form an additional unknown-class test set. Its
composition is determined by sparse reference sampling, so it is not a random
holdout of well-sampled clades. Unknowns never tune the model or its threshold.

Locus normalization uses only positive anchor-to-anchor distances. Supervised
panel selection ranks candidate anchors by between-class / within-class separation
of training query distance profiles, retaining at least one anchor per known
class. The default selects 48 of 54 candidates on the included dataset. This is
a modest reduction, not a demonstrated optimal panel. No fitted component uses
test species. Each training species receives equal total weight across its
available loci.

PCoA fits the mean normalized anchor distances, reports negative eigenvalue mass
and positive eigenvalue mass retained in two dimensions, and freezes coordinates.
Out-of-sample PCoA projection never recomputes the reference map around a query.
Missing projection distances use training-query column medians. The two-dimensional
K-means baseline uses anchor labels to map clusters to taxa. It is an analogue of
the abandoned strategy, not an exact reproduction of the 2023 implementation.

Every query locus has three channels at every anchor: `exp(-normalized_distance)`,
valid-site overlap fraction, and availability. Missing similarity is zero with
explicit coverage/mask channels. Gaussian splatting produces a fixed 12×12×3
tensor. The raster classifier and full-profile classifiers use the same regularized
logistic regression (C=0.1). The seven paired comparisons are learned panel,
stratified random panel, all candidate anchors, raster pixels, projected 2D
logistic regression, projected 2D K-means, and nearest anchor.

Missing anchor pairs fail by default. `--atlas-missing shortest-path` explicitly
fills only missing layout distances through the observed anchor graph, preserving
all observed distances and reporting the imputed pair count. A disconnected graph
still fails. Query profiles retain their original missingness masks. In the
MrKnowItAll global-family runs, seeds 42/43/44 impute 0/7/0 pairs, respectively.
This approximate geometry is a limitation of the image representation, not an
imputation of query genetic distances.

Query-species predictions average available locus scores equally. A validation-only
confidence threshold targets at least 90% empirical precision over at least five
validation species. If no cutoff qualifies, all queries abstain. This threshold
does **not** guarantee test precision or detect unknowns reliably. Model scores
are uncalibrated. The report records test accuracy, macro F1, abstention coverage,
selective accuracy, unknown false acceptance and confidence AUROC.

## Outputs

Each benchmark run writes:

- `metrics.json`: all species partitions, results, software versions, source-file
  checksums and experiment-code checksum.
- `panel.csv`, `predictions.csv`: selected anchors, fixed coordinates, held-out
  species predictions for the five supervised comparisons.
- `model.joblib`: fitted learned-panel and raster classifiers, frozen atlas,
  normalization, class labels and validation cutoffs.
- `locus_patterns.npz`: image tensors for **all usable partitions**, with species,
  partition and locus labels. Split on the provided partition field before any
  further ML training. This file is not an independent test set.
- `patterns.png`, `comparison.png`, `surface_*.png`: actual input patterns,
  comparisons and evidence surfaces.
- `synthetic_mixture_*.json`: a deterministic pool of three held-out reference
  species from distinct classes. This is not a FASTQ/assembly mixture experiment.

For the surface, each class score is divided equally among its anchors and
displayed with a Gaussian kernel. Heights show relative model evidence. Nearby
classes may merge visually. No surface peak count is interpreted as species count,
and no height is interpreted as DNA abundance. The renderer smooths the display;
classification does not read the smoothed evidence surface.

Large generated data/models live under ignored `runs/`. Small metrics, predictions,
panels and representative figures for seeds 42/43/44 are checked into `results/`.

## Train on exon-tree distances and predict new profiles

To test SPrOUT patristic distances, first build a **compatible training dataset**
from reference-only raw square matrices. Supply a CSV manifest:

```csv
locus,path,coverage_path
4691_exon_1,reference/4691_exon_1.matrix,reference/4691_exon_1.coverage.csv
4691_exon_2,reference/4691_exon_2.matrix,reference/4691_exon_2.coverage.csv
```

Each distance CSV has leaf names as the first column and header row, exactly as
stage 03 raw `.matrix` files. Coverage CSVs have identical labels/layout, containing
shared A/C/G/T sites divided by homologous alignment length, in [0,1]. Calculate
these from the corresponding alignments. Reference preparation requires a single
accession per species/locus. Query leaves belong only in prediction inputs.

```bash
python -m experiments.reference_atlas.atlas prepare-matrices \
  --manifest reference_manifest.csv --distance-kind patristic \
  --output experiments/reference_atlas/runs/exon_tree_distances.npz

python -m experiments.reference_atlas.atlas benchmark \
  --dataset experiments/reference_atlas/runs/exon_tree_distances.npz \
  --output experiments/reference_atlas/runs/exon_tree_model

python -m experiments.reference_atlas.atlas predict-matrices \
  --model experiments/reference_atlas/runs/exon_tree_model/model.joblib \
  --manifest query_manifest.csv --query-prefix my_sample \
  --distance-kind patristic --method distance_learned_panel \
  --output experiments/reference_atlas/runs/query
```

The query manifest uses the same columns with query+anchor raw matrices and
matching coverage CSVs. Paths are relative to the manifest. Training and query
locus IDs must match exactly and represent the same homologous sites. This avoids
silently interpreting stage-03 enumeration (`4691.1`) as exon identity. Every raw
query leaf must start with `--query-prefix`; unknown names are never parsed as
reference labels. Missing anchors are masked, duplicate anchors/query-locus pairs
and flagged distances (`999`) fail explicitly. At least one quarter of panel
distances, and at least two anchors, must be observed for a usable query locus.
Unusable query loci are omitted; no usable query profiles causes an error.

Prediction writes per-query/locus scores and a mixture evidence surface. Different
exons of one contig need externally consistent IDs for subsequent contig-level
fusion. Mixture evidence averages contigs within locus, then averages loci. It
still needs calibrated presence thresholds and correlated-evidence handling.

The sample p-distance model cannot consume patristic distances, and the CLI rejects
that mismatch. Distance-kind agreement alone does not ensure compatibility: branch
length estimation, locus identities, alignment coordinates and reference processing
must also agree. A frozen map does not freeze tree branch lengths when new queries
are added; this potential domain shift needs an independent placement benchmark.

The generic `--rank` switch accepts order/family/genus/species. Ranks with too few
independent species fail explicitly. Under this unseen-species benchmark, exact
species is intentionally unidentifiable from its own training class: it needs
multiple specimens/accessions and a different, specimen-disjoint evaluation.
Automatic hierarchical routing and open-set calibration remain research work,
specified in the design document.

## Iteration on the 100 existing TTU mixtures

[Completed 100-sample results and diagnostics](results/transfer100_20261016/REPORT.zh-CN.md)
include both the literal supplementary-table evaluation and the previously
established two-taxon harmonization. All 100 predictions succeeded; the input
hash audit and independent ZIP/XML truth/count audit passed.

`transfer.py` fits a shared match verifier from reference-only leave-species-out
pseudo-queries on existing family exon trees. This adapts the full-profile idea
to patristic distances without applying the previous p-distance model to a new
distance domain. The closest family by full distances proposes a fragment match;
the learned verifier reads the entire distance-shape residual, correlation,
terminal offset and competing-family evidence. It is not permitted to rank distant
families by an extrapolated binary score.

`batch_transfer.py` freezes the 100-sample manifest, model and code, then scores
immutable existing family trees under identical reused order predictions. It
rechecks input hashes, produces per-fragment and independent-gene evidence, and
collects results only after all 100 succeed. No assembly, alignment, reference
selection or tree inference is repeated. The primary rule requires three distinct
genes supporting a verified family; there is no sample-wide family count cap.
Nearest-reference gene voting and verification of the old final candidates are
predeclared comparisons. Thresholds use species-balanced known-reference versus
whole-family-absent reference controls, never mixture truth.

`evaluate_transfer.py` runs locally after completion. It needs `openpyxl` to read
the original STable and outputs literal-table and previously established two-taxon
harmonization views, along with TP/FP/FN/TN, precision, recall, F1, NPV and MCC.
The truth workbook and evaluator are not uploaded to the prediction host. This
is a retrospective iteration on previously evaluated mixtures, not an independent
new blind validation. The HPC runtime versions are recorded in its protocol.

The recorded HPC model used Python 3.13.7, NumPy 2.3.3, SciPy 1.16.2,
scikit-learn 1.7.2, Biopython 1.85 and joblib 1.5.2. To repeat the run, copy
`transfer.py` and `batch_transfer.py` together into a new empty HPC directory,
then invoke the scripts directly:

```bash
python transfer.py train --source-root /path/to/completed100 --output /path/to/new-run/model
python batch_transfer.py freeze --root /path/to/new-run --source /path/to/completed100
python batch_transfer.py submit --root /path/to/new-run
```

The batch helper targets TTU's `nocona` partition and the interpreter path stated
in the script. It reuses the completed batch's `selection.json`, `01_order` calls,
validated `02_family` trees and existing final/native candidate lists. Use a new
run directory to preserve previous results. `AUDIT.json` and `compact_results.tar.gz`
are produced after all 100 succeed. In the recorded run, resource-only scheduling
adjustments and completion of file verification/packaging are documented separately
in `RESOURCE_UPDATE.json`; no prediction code or scientific parameter changed.

After downloading and safely extracting the compact archive, run the evaluator
and the independent count auditor locally from the SPrOUT checkout:

```bash
python -m experiments.reference_atlas.evaluate_transfer \
  --batch /path/to/extracted-batch --truth /path/to/STable.xlsx --output /path/to/evaluation
python -m experiments.reference_atlas.audit_transfer \
  --batch /path/to/extracted-batch --truth /path/to/STable.xlsx --evaluation /path/to/evaluation
```

The auditor reads workbook ZIP/XML separately from openpyxl and recomputes calls
directly from prediction files. The checked-in result directory compresses the
large per-family table as `by_family.csv.gz`, and includes all five candidate sets
per sample in `predictions.json`. Full per-fragment features remain in the HPC run.

## Combined C, auxiliary gene rescue and ML prediction sets

[Completed fusion results and exported 100-sample lists](results/fusion100_20261016/REPORT.zh-CN.md)
retain the two-gene ML-confirmed C core, with an auxiliary route for C gene-rescue
candidates having a stability Monte Carlo lower bound >0.5 and at least one
ML-verified gene. `merged` is the deduplicated union. B/C/three-gene ML candidates
outside that union remain in `review`, which is not final presence. Channels
share trees and do not yield independent biological presence probabilities.

The fusion rule was developed after seeing earlier errors on these mixtures;
it is a retrospective improvement, not independent validation. The literal-table
and established two-taxon-harmonization views both remain in the report, including
the independent ML comparator that has better literal-table F1. No taxon-specific
or sample-specific exception is used.

```bash
python -m experiments.reference_atlas.fuse \
  --transfer-batch /path/to/transfer100/batch \
  --legacy-batch /path/to/nocap100/batch \
  --statistical-panels /path/to/statistical_evidence/panels \
  --output /path/to/new-fusion
python -m experiments.reference_atlas.evaluate_fusion \
  --predictions /path/to/new-fusion --truth /path/to/STable.xlsx
```

`fuse` refuses an existing output directory and verifies source C memberships,
family panels, parent labels, stability membership and ML threshold before joining.
Its inputs are the completed original C batch, audited ML batch and existing
statistical panels, not sample truth. It writes merged/core/auxiliary/review files
for all 100 samples inside `prediction_sets.zip`, a flat `merged_predictions.csv`,
the full gzipped evidence snapshot, the frozen policy and a hash manifest.
`evaluate_fusion` subsequently reads truth, computes both views and independently
recounts ZIP/XML truth and exported lists. It also reports the cost of adding all
review candidates, so deferral does not hide false negatives.
