# Complete the retained genus pool

This is a development correction after inspecting the two sealed ten-mixture
evaluations, not a new blind or independent validation.

The strict search retained all 60 expected sample–genus events through every
routing round. Its partial final-panel fallback excluded 40 of those events as
unresolved; merged recall was 20/60, with 41 false positives. v1 excluded 54/60.
Scoring all query fragments against the small locally resolved subset also
offered close relatives as substitutes for untested true genera.

Keep both baseline batches unchanged. Reuse the exact strict search, assemblies,
family parents, candidate lists, QC rules, fixed reference coordinates, genus
model and final classifiers. Test the entire retained pool jointly when it is
at most 512 genera, after the existing ≤64 target or eight-round routing limit.
The five existing full-pool final runs can be reused; complete the five partial
runs with retained counts 93, 187, 289, 355 and 447. These are the whole frozen
candidate pools, not truth-selected additions. Genus count is never capped by
the expected number of components. If a future retained pool exceeds 512,
report the whole region unresolved rather than classify its fragments against
an arbitrary small subset.

Reuse only QC-passing projected reference sequences whose original files and
metadata match the sealed strict-run hashes. New recoveries undergo the same
target-gene QC and alignment procedure. Prefer up to two species per genus as
before; unavailable genera remain explicit and count as misses. Do not lower
QC, add named expected genera or change ML/C-style thresholds using truth.

Report this correction alongside both baselines. Separate reused source cost
from incremental completion cost and disclose warm caches. Reusing these same
ten mixtures can establish that this implementation defect was corrected, but
cannot establish independent accuracy on new mixtures.
