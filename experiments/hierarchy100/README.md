# Fresh order-to-genus inference on 100 mixtures

Branch: `experiment/order-to-genus100-20261022`.
The [frozen plan](PLAN.md) defines the 100-mixture retrospective validation.
All order/family/genus inference is rerun. Only query exon assemblies, public
reference resources and the two existing verifier checkpoints are reused.
The previous ten genus-development mixtures and other ninety are reported
separately. Neither subset is a new biological cohort.

All 100 tasks and the collector completed successfully. The [final report](results/ttu100_20261022/REPORT.md)
records primary family precision/recall/F1 of 87.5%/88.5%/0.880 and genus
74.2%/72.0%/0.731. Ten searches remain above the 512-genus work limit, accounting
for 100 missed genus occurrences; reference coverage excludes another 50, and
routing excludes seven. Those misses stay in all-100 denominators. The original
ten development genus call sets reproduce exactly; other-ninety genus F1 is 0.717.
Primary outputs are in [primary_predictions.csv](results/ttu100_20261022/primary_predictions.csv).
All 54 local tests pass. A local completion-only import dependency failure was
repaired and preserved in the artifact audit; remote prediction code, models,
policy and outputs were not changed.

TTU root: `/lustre/scratch/nhu/202609/ttu_hierarchy100_20261022`.
Preparation job `26230805` completed successfully. Array `26230808` contains
all indices 0–99, without a workflow concurrency cap. The array and collector
`26230810` both completed. The initial submission records and final result
records are kept separately.

Every sample has independent output directories and executes order → family →
genus in sequence. Each array task requests eight CPUs, 32 GB and twelve hours.
The scheduler controls simultaneous allocations. Per-reference file locks
prevent duplicate concurrent Kew writes, without locking unrelated references.
The archived legacy tree wrapper uses one MAFFT thread per gene worker to stay
within allocated CPUs. These operational source adaptations are recorded in
the payload manifest.

Order admission uses two distinct genes, Z ≥ 0.9 or rank 1, without a candidate
count cap. Family merged calls use the frozen C/stability/ML fusion (ML threshold
0.94). Fresh original-Z, merged and review families define the genus search pool.
The genus policy is unchanged from the complete-pool ten-mixture development:
30-of-50-locus reference coverage, threshold 0.83, up to eight rounds, 64-genus
routing target and 512-genus final work limit. Larger unresolved pools abstain.
Primary genus calls use the already selected auxiliary output path; every other
fixed method is preserved. No parameters change in response to new truth scores.

## Reproduce the submitted payload

The checked-in launch files specify this run's exact paths and cohort. Use a
new output root when adapting them; do not overwrite a frozen run.

```powershell
New-Item -ItemType Directory -Force experiments/hierarchy100/runs/ttu_20261022
Copy-Item experiments/hierarchy100/launch/* experiments/hierarchy100/runs/ttu_20261022/
python -m experiments.hierarchy100.pack
```

Packing requires the recorded family/QC source snapshot and model checkpoints;
it fails if they are absent. Runtime sources, both models, configuration and
batch scripts are hash-frozen by the preparation batch job. The exact uploaded
payload is retained separately from later local evaluation/report helpers.

After the uncapped inference array and collector finish, verify the remote
archive digest, download and extract the compact result archive safely. The
truth workbook stays local. Seal first, then evaluate:

```bash
python -m experiments.hierarchy100.evaluation seal --batch /path/to/received-batch --output /path/to/new-export
python -m experiments.hierarchy100.evaluation evaluate --batch /path/to/received-batch --output /path/to/new-export --workbook /path/to/original-STable.xlsx
```

The gate requires all 100 audited predictions and checks per-sample text lists
against combined JSON, including empty calls. Literal Table S1 source IDs,
scientific names, orders and families are checked with independent spreadsheet
and ZIP/XML readers. Parent/reference omissions and unresolved genera remain in
the recall denominator. Missing jobs are reported separately and block scoring.
The supplied workbook is never modified or uploaded.

Per-phase timing separates order alignment/trees, distances and admission;
family reference construction, trees, distances, C/stability and ML; and the full
genus search. Genus SEARCH/ROUND records retain its finer phase/round breakdown.
Shared preparation/audit, queueing and historical FASTQ assembly are outside
fresh inference timings. This is a warm-cache benchmark; it does not establish
a cold-install or matched-exhaustive speedup.

## One-shot local completion

`finish_local.py` can continue this authorized batch after the chat turn ends.
It waits for the collector's Slurm COMPLETED state, verifies and safely downloads
the archive, seals all 100 outputs, scores the local truth workbook, renders the
report and publishes only this run's result package on the experimental branch.
It does not submit recurring jobs, upload truth, edit thresholds or skip failed
samples. It stops for attention on scientific failures or changed local evaluator
sources. SSH interruptions are retried for up to 48 hours. The TTU VPN must remain
available on the host for transfer; disconnecting the VPN does not stop submitted
Slurm computation. State is recorded in `runs/ttu_20261022/CONTINUATION.json`.
Plots are generated automatically; human/agent visual inspection is recorded as
pending until it is actually performed.

Use the bundled Python runtime for the local XLSX truth reader. Scientific plots
use the existing `.venv312` environment, which supplies Matplotlib. No runtime
installation or shared environment modification is performed. A Windows
background launch uses a hidden window and separate stdout/stderr logs.
