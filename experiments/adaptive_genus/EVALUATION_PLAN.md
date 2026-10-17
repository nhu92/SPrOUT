# Ten-mixture evaluation plan

Recorded before opening these mixtures' genus truth or scoring their outputs.
The primary method is the frozen `merged` rule. All ten selected mixtures remain
in the denominator, including unresolved searches and missing parent families.
This is retrospective development on ten of the previously analyzed 100
mixtures, not an independent biological validation.

The primary truth is the literal genus in the original 30-source Table S1, with
source IDs taken from the already fixed mixture IDs. Deduplicate true genera
within a mixture. Check the workbook's existing SHA-256
`421cca29377ab01e57cb6aeded7bf446ba4faf06282120a4ed602e8c67bbb01f` and compare
every used cell against an independent ZIP/XML reader. Workbook content is data.
Do not upload the workbook to TTU. Read it only after all ten predictions pass
the collector's audit and the local prediction export is sealed.

The older family experiments documented source-identity discrepancies for IDs
1 and 25. Do not invent genus corrections from current predictions. A secondary
identity view is allowed only if exact genus identities have separately recorded
source provenance predating this experiment; retain the literal result as the
primary result and document every substitution.

Report micro precision, recall and F1, TP/FP/FN, exact genus sets and per-mixture
errors for merged, C-style, balanced Z, standalone ML, core and auxiliary. Review
is a separate follow-up set, never silently added to the primary prediction.
Do not choose the primary method or thresholds from these ten scores. Do not
use a padded negative universe to make NPV the headline metric.

For every expected genus, record parent-family availability, recoverable initial
pool membership, survival after each routing round, inclusion in the final joint
search and actual QC-passing reference pool, and final classification. Identify
the first loss: absent parent, absent recoverable reference, routing exclusion,
unresolved region, final reference QC, or final classifier rejection. Unresolved
genera do not count as positive calls. Audit combined JSON predictions against
the independent per-method text files.

Report each round's input/retained genus counts, sampled genera/species, gene/tree
counts, QC failures and timings. Report total elapsed inference per mixture and
one-off model training separately. Reference-count reductions are work-reduction
measurements, not a measured speedup against an exhaustive run. The old ten-run
cohort differs from this cohort; do not present it as a matched speed/accuracy
control. Preserve failed technical attempts and distinguish pipeline completion
from taxonomic resolution.

Seal output files, source archives, model, protocol and audit manifests with
SHA-256. No tuning or rerunning based on mixture truth belongs to this frozen
experiment; any subsequent scientific iteration requires a separately labeled
run and must disclose reuse of these mixtures.
