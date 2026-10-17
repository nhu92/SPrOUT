# Adaptive genus search: ten-mixture TTU experiment

Completed both predeclared routing policies on the same ten in-silico mixtures: **20/20 pipelines audited**. Both batches were sealed before the original Table S1 genus truth was opened. The stricter policy was specified from truth-blinded search stagnation, with the same model and query data.

| Routing policy; merged classifier | TP | FP | FN | Precision | Recall | F1 | Exact sets | Median retained genera | Median search minutes | ≤64 genera / full QC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| v1 | 6 | 41 | 54 | 12.8% | 10.0% | 0.112 | 0/10 | 381 | 14.8 | 0/10 / 0/10 |
| strict_route | 20 | 41 | 40 | 32.8% | 33.3% | 0.331 | 0/10 | 74.5 | 13.7 | 5/10 / 0/10 |

| Policy | Components | Samples | TP | FP | FN | Precision | Recall | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| v1 | 3 | 4 | 4 | 10 | 8 | 28.6% | 33.3% | 0.308 |
| v1 | 6 | 3 | 1 | 11 | 17 | 8.3% | 5.6% | 0.067 |
| v1 | 10 | 3 | 1 | 20 | 29 | 4.8% | 3.3% | 0.039 |
| strict_route | 3 | 4 | 12 | 8 | 0 | 60.0% | 100.0% | 0.750 |
| strict_route | 6 | 3 | 6 | 13 | 12 | 31.6% | 33.3% | 0.324 |
| strict_route | 10 | 3 | 2 | 20 | 28 | 9.1% | 6.7% | 0.077 |

v1 requires two routing genes (three for ambiguous second routes) and four nearby guard genera per supported anchor. The strict policy requires five genes for both routes and one guard genus. Both use 32 initial tree cells, 64 after a stalled round, at most eight rounds and a 64-genus final joint-panel budget. These constrain reference work, not the number of positive predictions.

Each round chooses central, target-locus-QC-passing Kew representatives of tree neighborhoods, aligns the same query exons against the reduced references, infers new exon trees and retains regions with repeated independent-gene support. Surviving regions are sampled again. The Kew planning tree supplies topology distances, not dated evolutionary branch lengths. Final calls combine the genus C-style rule and a reference-trained genus ML verifier. Coarse anchor genera are never promoted directly into final calls.

## Results and reusable outputs

- [Strict routing: detailed results, trajectories and errors](strict_route/REPORT.md)
- [Original routing: detailed results, trajectories and errors](v1/REPORT.md)
- [Strict routing merged genus predictions](strict_route/merged_predictions.csv)
- [Original routing merged genus predictions](v1/merged_predictions.csv)
- [Paired sample comparison](comparison_by_sample.csv)
- [Input/model/policy comparison audit](COMPARISON_AUDIT.json)

Review sets, unresolved genera, alternative final rules and all per-round candidates are retained in each result directory; none are silently added to the merged predictions. Each detailed report identifies the first stage losing a true genus, including upstream parent misses and reference QC failures.

## Interpretation limits

This is retrospective development on previously analyzed mixtures drawn from the same 30 source species. The ten selected mixtures contain four 3-component, three 6-component and three 10-component samples. A preference based on these ten results needs validation on a separate cohort. The old progressive ten-mixture cohort is different, so it is not a matched accuracy/runtime control.

Reference-pool reduction is measured; speedup versus exhaustive genus inference is not. The second policy can benefit from earlier downloads into the shared cache. Concurrency was a throughput choice rather than a controlled timing variable: v1 allowed four tasks, while strict routing increased from two to four to six. Per-sample search minutes include reference sampling, fetching/QC, reference projection, query alignment, trees and scoring; they exclude prior assembly, initial/final audits, Slurm setup/cleanup and queueing. Allocated job durations and CPU accounting are archived separately. Shared model training costs 5.9 minutes and is not counted again for the second policy.

The original Table S1 is the literal scoring source. Prior family-level source discrepancies for IDs 1 and 25 are disclosed; no genus corrections were inferred from the new predictions. Species presence, abundance, negative-only mixtures and independent ML generalization are not evaluated.

Branch: `experiment/adaptive-genus-search-20261017`. All implementation and evaluation code is separate from previous production defaults. The exact frozen inference payloads, shared trained model, source hashes, Slurm records, raw result archives and independent audits are included. See [design](../../DESIGN.md) and [evaluation plan](../../EVALUATION_PLAN.md).
