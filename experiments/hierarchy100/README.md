# Fresh order-to-genus inference on 100 mixtures

Branch: `experiment/order-to-genus100-20261022`.
The [frozen plan](PLAN.md) defines the 100-mixture retrospective validation.
All order/family/genus inference is rerun. Only query exon assemblies, public
reference resources and the two existing verifier checkpoints are reused.
The previous ten genus-development mixtures and other ninety are reported
separately. Neither subset is a new biological cohort.

TTU root: `/lustre/scratch/nhu/202609/ttu_hierarchy100_20261022`.
Preparation job `26230805` completed successfully. Array `26230808` contains
all indices 0–99, without a workflow concurrency cap. Collector `26230810`
depends on array completion and blocks the all-100 result seal if a task fails.
Submission records indicate launched work, not completed validation results.

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
