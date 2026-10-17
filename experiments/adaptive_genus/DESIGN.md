# Adaptive genus search experiment (2026-10-17)

This branch tests repeated phylogenetic sampling after the existing order/family
stages. It does not change the previous family experiment or production defaults.

* Choose 10 of the completed 100 in-silico mixtures deterministically: four
  3-component, three 6-component and three 10-component mixtures, ranked within
  each stratum by SHA-256 of `adaptive-genus-20261017/<sample>`. Do not inspect
  their genus truth to select samples, representatives or stopping thresholds.
* Reuse the existing assemblies and search families from original Z, frozen
  merged family calls and review combined. These are search parents, not new
  final family-presence calls. Missing true parents remain end-to-end misses.
* Initially partition the recoverable Kew tree inside these families into at
  most 32 neighborhoods; allow 64 after a stalled round. Split the largest
  active regions along actual tree structure.
  Polytomies remain forests. Keep all recoverable placements of a genus.
* Download and QC central representatives of each neighborhood. Centers minimize
  maximum/mean distance to a genus-balanced landmark cover; they are not the
  first alphabetical genus. Prefer two different representative genera per
  region and enough accepted query-target loci. Failed representatives are
  replaced locally; an uncovered region is unresolved, never negative evidence.
* Each round builds a reduced Kew panel, aligns the same query exons, infers
  FastTree GTR+Gamma trees and aggregates support by distinct gene. A sampled
  anchor genus represents its neighborhood, not an asserted true genus.
* Retain regions with repeated-gene evidence and evidence-based ambiguity paths.
  Add only a small set of nearby genera around supported anchors, rather than
  whole sister branches. Repartition surviving tree placements and repeat.
  Record original pool, selected references, retained genera, quality failures,
  per-gene evidence, discarded regions, elapsed time and work at every round.
* Stop for a joint genus panel when the active pool is at most 64 genera. Stop
  as unresolved if repeated routing does not shrink it or the round budget is
  exhausted. Do not silently force a final top-N taxonomic answer.
* Reuse the seven full-distance profile features, but fit a genus-rank verifier
  from reference-only trees and calibrate it with held-out reference genera and
  genus-absent controls. The family checkpoint is not a calibrated genus model.
  Final comparisons include balanced genus Z, C-style gene rescue/confirmation,
  standalone ML, the analogous core/auxiliary merge and review. Their genus
  application is a new experiment, not an unchanged family validation.
* Build reference-only training from already-QC'd Kew genus alignments where
  available. Cache recovery QC and aligned reference coordinates; record this
  one-off cost separately from sample inference, along with download and queue
  costs. Do not claim measured end-to-end acceleration from reference-count
  reductions alone.
* Freeze the test cohort, policy, model and code before the full ten-sample
  run. Mixture truth is evaluated locally after prediction outputs are sealed.
  Report per-round true-genus survival after unblinding, full-set precision,
  recall/F1, parent/reference misses, final candidate counts and actual timing.
  Development/reference-only checks are not biological validation.
