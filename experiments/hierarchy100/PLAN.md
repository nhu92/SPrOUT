# Order → family → genus: 100-mixture TTU rerun, 2026-10-22

Run all 100 existing paired-read mixtures. Reuse only completed query exon
assemblies and frozen reference resources. Rebuild order and family trees,
distance matrices, predictions and genus searches in fresh output directories.
Never import previous order/family/genus prediction lists. This starts at order
inference, rather than repeating FASTQ assembly. Cached public Kew DNA is reusable.

Use the existing 50 target genes, order reference panel and corrected backbone.
Order admission: at least two **distinct genes**, aggregate Z >= 0.9 or aggregate
rank 1, without a taxon count cap. Family reference construction and C, stability,
ML fusion use the existing real-validation code and unchanged family verifier
(threshold 0.94). Genus parents are fresh family original-Z ∪ merged ∪ review.
Genus uses the complete-pool policy from the ten-mixture experiment, unchanged
reference coverage gate (30/50), eight rounds, target 64, completion limit 512,
and unchanged verifier threshold 0.83. Larger pools abstain explicitly.
Primary final genus output is the existing stable auxiliary path, selected
during the previous ten-mixture development. Preserve all other output paths.

Freeze source, models, input hashes, policies, sample list and this plan before
inference. No truth workbook is uploaded. Seal all hundred sample results before
local scoring. Count unresolved, coverage-excluded and absent-parent taxa as
misses, rather than removing them from denominators. Incomplete jobs block the
all-100 evaluation and remain visible in STATUS.json.

Score order admission, family merged and genus auxiliary against literal Table S1
G/H/first-word B. Verify every used cell independently through ZIP/XML. Preserve
known taxonomy annotation discrepancies; do not derive name replacements from
new predictions. Report all-100, previous genus-development ten, other ninety,
mixture-size strata, per-sample errors, first-loss stage and phase timings.
The cohort and thirty source species have been analyzed before. The other ninety
were not used for the previous genus development, but are not an independent
biological cohort or a previously unseen species validation.

Submit substantial work as TTU Slurm batch jobs. One array element processes its
own order→family→genus dependencies; all 100 elements are submitted without a
workflow concurrency throttle. Scheduler/account limits determine concurrency.
Each task uses eight allocated CPUs; the legacy tree function is invoked with
one MAFFT thread per independent gene worker to avoid nested oversubscription.
Outputs and projected reference caches are sample-local. Shared raw Kew cache
downloads use per-recovery advisory locks to avoid concurrent writes to the same
file; unrelated references and samples remain independent. Immutable cache
metadata and backbone hashes are checked before and after each task.

Training is not repeated. Report fresh inference times separately from shared
preflight/audit jobs, queueing and historical assembly. Cached reference DNA and
models make this a warm-cache inference benchmark, not a cold-install timing.
