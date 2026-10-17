# Adaptive phylogenetic genus search

Experimental branch: `experiment/adaptive-genus-search-20261017`. The cohort,
search policy and scientific choices are recorded in [DESIGN](DESIGN.md),
[policy.json](policy.json) and [cohort.json](cohort.json).

This replaces alphabetical branch anchors and whole-sister-branch expansion with
genus-balanced central representatives, target-gene sequence QC, repeated
independent-gene routing and small nearby-genus guards. A representative genus
stands for a search region until a final joint genus run. Polyphyletic genera
retain multiple recoverable placements. Unresolved regions stay explicit.

Each coarse panel has up to two representatives per tree cell. The initial
budget is 32 cells; one stalled round permits 64 cells. Surviving regions are
subdivided again. The joint stage starts at at most 64 candidate genera, using
up to two different species per genus. These are reference-work budgets, not a
cap on the number of predicted genera or an expected mixture size.

Kew sequences are retrieved lazily, checked against the complete corrected
backbone for each of the actual target genes, and required to pass at least 60%
of those genes (minimum ten). The same reference's QC and projected coordinates
are cached within a sample. Every round still aligns queries to the reduced
panel and infers new full exon trees. References alone determine trimAl columns,
so additional query contigs cannot change the reference trimming threshold.

The final classifier uses the seven full-distance profile features from the
family experiment, with a **new genus-specific reference-trained verifier**.
Entire genera are kept out of fitting for calibration and reference testing.
Known-reference acceptance and genus-absent rejection choose the threshold;
mixture truth does not. The C-style comparison uses a bounded genus-balanced
contrast and the existing gene-bootstrap/confirmation formulas, followed by
the same core/auxiliary merge structure. It is a new genus adaptation, not the
unchanged family checkpoint or unchanged raw ACS statistic.

The shared gene-statistics implementation is copied unchanged from the archived
SPrOUT family experiment (`gene_statistics.py`, original SHA-256
`e46c4d014b1621f2fd972b902fa951903f407eedfca930ca9861f1afd04b036e`).
The adapter changes the input axis to genera with family parents and labels its
public outputs accordingly. Stability is not biological presence probability.

## Run

Python scientific dependencies match the existing SPrOUT environment. Kew/QC
integration additionally imports `sprout_refbuilder` from the recorded
MrKnowItAll source snapshot. MAFFT, trimAl and FastTree must be on PATH. HPC
payloads include the exact source snapshot and its license. No global Conda
environment is modified.

```bash
python -m pytest tests/test_adaptive_genus.py -q
python -m experiments.adaptive_genus.select_samples
# Prepare a new HPC root with source, policy/cohort and completed reference model.
python -m experiments.adaptive_genus.runner prepare --root /path/to/new-run \
  --backbone /path/to/corrected-backbone
python -m experiments.adaptive_genus.runner run --root /path/to/new-run --index 0
# Indices 0–9 may run independently. Collect only after all ten complete.
python -m experiments.adaptive_genus.runner collect --root /path/to/new-run
```

The TTU experiment reuses the completed 50-target-locus exon assemblies and
family search evidence, not FASTQ assembly. Outputs retain every round's
candidate genera, sampled species, QC, routing evidence and timings. A completed
pipeline may still be taxonomically unresolved: inspect `SEARCH.json` rather
than treating `SUCCESS` as a biological result. The original 100-mixture data
have been analyzed before; this is retrospective algorithm development, not a
new independent cohort. No species-level presence or abundance claim is made.
