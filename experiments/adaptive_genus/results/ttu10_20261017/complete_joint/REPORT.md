# Ten in-silico mixtures: adaptive phylogenetic genus search

All ten TTU pipelines completed and passed frozen input/model/code/product checks. The predeclared merged genus method recovered 54 true sample–genus events, with 30 false positives and 6 missed events: precision **64.3%**, recall **90.0%**, F1 **0.750**.

The search removed a median 91.6% of the initial recoverable genus pool. 5/10 searches reached the ≤64-genus search target. 0/10 searches reached a fully QC-represented final joint panel; 10/10 retained explicit unresolved regions or final reference omissions. Combined source-plus-completion inference was 5.8–47.8 minutes (median 14.8); shared reference-only model training took 5.9 minutes.

This is a retrospective method experiment on ten mixtures from the previously analyzed 100. The 4/3/3 selection across 3-, 6- and 10-component strata was fixed by sample-ID hashes, without genus truth or prior success scores. It is not an independent biological validation. Mixture component counts never determine the number of predicted genera.

This is the post-evaluation complete-pool correction. The original two batches were sealed and evaluated first; their failure diagnosis motivated this stopping-rule change. The same ten mixtures and unchanged final model are reused. It is an in-sample development result, not a new blind validation.

## Search behavior

The existing family search parents are the union of raw Z, frozen merged and review families. Inside those parents, partition the recoverable Kew tree into 32 neighborhoods, subdivide each into representative regions, and choose central genera using genus-balanced minimax/mean distance to a farthest-point landmark cover. The supplied Kew tree has no branch lengths; these planning distances count topology edges, not evolutionary time. All recoverable placements of polyphyletic genera remain eligible.

Representative sequences must pass MrKnowItAll QC on at least 30 of the actual 50 target genes. Download only proposed recoveries, replace failed representatives locally, and preserve uncovered regions as unresolved. Cache QC and projected reference coordinates per sample. At each round, align the existing query exons to the reduced panel and build fresh full exon trees. 5 independent genes retain a region; a close second route needs 5 genes. Nearby guards add 1 genera per supported anchor, within the previous active pool. The anchor genus labels its region until the final joint test.

Repartition and repeat. A stalled round permits 64 cells; two consecutive ≥95% retention rounds or eight total rounds stop routing. After the normal ≤64 target or eight-round limit, jointly verify the entire retained pool when it contains ≤512 genera. Larger pools remain wholly unresolved; no small-cell subset is substituted. Five original full-pool final tests are reused; the five formerly partial tests are completed with unchanged cached projections and newly QC-screened references. No biological top-N prediction cap is used.

![Search pool trajectories](search_reduction.png)

| ID | Mixture | Search genera: initial → retained | Final tested genera / species | Rounds | Minutes | Status | TP / FP / FN |
| --- | --- | ---: | ---: | ---: | ---: | --- | ---: |
| S1 | 07x26x30 | 1234 → 37 | 33 / 41 | 3 | 6.9 | unresolved_partial | 3 / 2 / 0 |
| S2 | 18x25x12 | 316 → 56 | 35 / 53 | 2 | 5.8 | unresolved_partial | 3 / 1 / 0 |
| S3 | 30x04x03 | 956 → 19 | 17 / 22 | 3 | 7.2 | unresolved_partial | 3 / 3 / 0 |
| S4 | 26x01x24 | 432 → 39 | 36 / 47 | 2 | 5.8 | unresolved_partial | 3 / 2 / 0 |
| S5 | 30x04x19x07x03x16 | 2224 → 187 | 154 / 196 | 8 | 23.2 | unresolved_partial | 5 / 0 / 1 |
| S6 | 30x04x13x26x19x28 | 1098 → 93 | 58 / 79 | 8 | 16.3 | unresolved_partial | 4 / 6 / 2 |
| S7 | 04x07x08x15x16x01 | 1361 → 56 | 44 / 55 | 8 | 13.3 | unresolved_partial | 6 / 3 / 0 |
| S8 | 04x30x25x16x10x03x28x01x13x08 | 3719 → 289 | 260 / 341 | 8 | 47.8 | unresolved_partial | 9 / 2 / 1 |
| S9 | 25x28x30x07x22x24x12x04x03x18 | 2143 → 355 | 288 / 373 | 8 | 41.4 | unresolved_partial | 9 / 4 / 1 |
| S10 | 07x22x30x08x18x13x03x19x04x02 | 1861 → 447 | 391 / 503 | 8 | 43.8 | unresolved_partial | 9 / 7 / 1 |

Across 58 routing rounds, median per-round reduction was 27.0%. Full per-round panels, counts, evidence, QC and timing are in the sealed batch archive and rounds.csv. The pilot and training costs are separate from the ten sample times. Slurm used eight CPUs per mixture with up to 3 mixtures of this phase concurrently; the strict-policy throttle increased from two to four to six as earlier work finished. Reported inference times exclude queueing, prior assembly, initial/final audits and Slurm setup/cleanup. Allocated job elapsed times are recorded separately in Slurm accounting. Fetch/QC/alignment durations overlap worker tasks and are not CPU-time totals.

Reference-count reduction measures avoided reference work. There is no exhaustive genus run on these same ten samples, so these results do not establish a measured speedup over the old pipeline. The previous progressive ten-mixture cohort differs from this cohort.

| Runtime phase | Median seconds per mixture | Sum of sample-phase minutes |
| --- | ---: | ---: |
| fetch_parse | 8.9 | 1.6 |
| reference_QC | 68.3 | 11.7 |
| reference_projection | 571.9 | 149.8 |
| query_align_tree | 172.4 | 32.2 |
| planning_routing_scoring_other | 66.1 | 16.2 |

Phase sums combine overlapping Slurm tasks and are not cohort wall time or CPU time. The residual phase includes tree planning, index loading, routing and final scoring. For this completion run, totals include the reused source pipeline and its original partial final test, plus new completion work. Incremental completion seconds and reused-final flags are in timings.csv; these are warm-cache observations.

## Genus predictions

Primary output: [merged_predictions.csv](merged_predictions.csv). The [complete method sets](prediction_sets.zip), [candidate evidence](candidate_evidence.csv) and [unresolved search regions](search_regions.csv) are separate exports.

| Method | TP | FP | FN | Precision | Recall | F1 | Exact sets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| merged | 54 | 30 | 6 | 64.3% | 90.0% | 0.750 | 0/10 |
| c_style | 54 | 141 | 6 | 27.7% | 90.0% | 0.424 | 0/10 |
| balanced_z | 31 | 239 | 29 | 11.5% | 51.7% | 0.188 | 0/10 |
| standalone_ml | 54 | 20 | 6 | 73.0% | 90.0% | 0.806 | 0/10 |
| core | 54 | 28 | 6 | 65.9% | 90.0% | 0.761 | 0/10 |
| auxiliary | 54 | 16 | 6 | 77.1% | 90.0% | 0.831 | 1/10 |
| review | 0 | 221 | 60 | 0.0% | 0.0% | 0.000 | 0/10 |

| Merged, by mixture complexity | Samples | TP | FP | FN | Precision | Recall | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 components | 4 | 12 | 8 | 0 | 60.0% | 100.0% | 0.750 |
| 6 components | 3 | 15 | 9 | 3 | 62.5% | 83.3% | 0.714 |
| 10 components | 3 | 27 | 13 | 3 | 67.5% | 90.0% | 0.771 |

Review is an unmerged follow-up set, not another primary call set. C-style uses a new genus-balanced distance contrast and the previous gene-bootstrap confirmation/rescue formulas. Merged keeps C-style calls with at least two ML-supported genes, or the predeclared stable auxiliary rescue with one ML-supported gene. Balanced Z is a within-panel standardized contrast; it is not the old raw family Z.

![Merged genus accuracy and elapsed time](accuracy_runtime.png)

## Where true genera were lost

| First loss / recovery stage | Sample–genus events |
| --- | ---: |
| final_reference_QC_unavailable | 6 |
| recovered | 54 |

![Retention of true genera through the pipeline](truth_survival.png)

[truth_survival.csv](truth_survival.csv) records every true genus surviving each round; [first_losses.csv](first_losses.csv) identifies the first exclusion. All ten mixtures contribute to recall, including absent parents, unavailable references, unresolved searches and final rejections. The original Table S1 is the primary literal truth. Older family reports documented source-identity discrepancies for IDs 1 and 25; no genus substitution was inferred from these predictions.

## Model and reproducibility

The new genus verifier used 397 Kew species across 297 genera and 50 reference-only gene trees. Entire genera are partitioned for fitting, calibration and testing. Known-genus and genus-absent pseudoqueries fixed the threshold at 0.83 before mixture prediction. Reference-held-out correct known-genus acceptance was 55.6%; genus-absent rejection was 96.5%. These are reference fragment checks, not ten-mixture accuracy or biological presence probabilities. The reference library came from previously constructed panels and is not an independent library. Source species were not explicitly excluded from the final Kew panels; this is a closed-reference experiment, not a test of genus recognition for wholly unrepresented source species.

The final query workflow also changes reference-column trimming and alignment projection. Accuracy/runtime differences cannot be attributed to the sampler alone. Sparse loci and cross-genus similarity can defeat routing; unresolved regions are visible instead of forced into representative labels. Species-level presence and abundance are not evaluated.

Initial inference implementation: `1fbf38d`; truth-blinded strict policy: `3f302fb`. Later commits add local evaluation gates, report figures and artifact audits. Branch: `experiment/adaptive-genus-search-20261017`. The previous pipeline defaults and family results are unchanged. Final ML/C-style/vote thresholds are unchanged; the stopping-rule correction was chosen after the baseline evaluation.

The [protocol](PROTOCOL.json), [prediction seal](PREDICTION_EXPORT.json), [independent evaluation audit](EVALUATION_AUDIT.json), [training record](model/TRAINING.json), runtime source snapshot, trained model, job records and batch archive preserve the run. The model training used its separately archived source snapshot; technical pilot failure records are retained. The original workbook and raw reads are not redistributed.

Truth source: original 30-source Table S1, SHA-256 `421cca29377ab01e57cb6aeded7bf446ba4faf06282120a4ed602e8c67bbb01f`. Every used source-ID, genus/name and family cell was checked with both openpyxl and ZIP/XML. Per-method plain text prediction files independently reproduce TP/FP/FN. The baseline truth had already been opened when this development correction was designed; new outputs were sealed before their rescoring.
