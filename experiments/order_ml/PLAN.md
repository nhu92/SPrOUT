# Order ML, 2026-10-23: fixed reference-only training and hundred-mixture scoring

Use the same hundred mixtures and frozen order trees from hierarchy100_20261022.
Do not rerun assembly, order alignment, order query trees or downstream family /
genus. Do not alter the old candidate lists, models, reference cache or outputs.
All substantial training and inference use TTU Slurm batch jobs. Submit every
independent sample array element without a workflow concurrency throttle;
training/freezing must precede prediction because predictions require its model.

Build fifty new reference-only gene trees from the original order reference
alignments, trimAl -gt .5 and FastTree GTR+Gamma, eight independent one-thread
workers. These training trees contain no mixture queries. Keep named species
with >=80 ACGT bases per used locus and >=30 of fifty loci. Collapse accessions
by species, choosing most valid bases and then lexicographic header. An order
with fewer than two eligible species can remain an anchor but cannot supply a
known-order held-out pseudoquery.

Split entire eligible multi-species orders by SHA256(order-ml-split/<order>):
first ten calibration orders, next ten test orders, all others fitting; require
at least ten fitting orders. Calibration/test orders do not appear in fitting
anchors. Reference phylogenetic scaffold geometry is shared and uses reference
DNA only; this is not an independent tree-scaffold test. Fitting pseudoqueries
remove every accession of their own species. Known rows include the correct
order and three nearest wrong orders. Additional negative rows remove the whole
pseudoquery order and use its nearest remaining candidate. Limit fitting
pseudoqueries per locus to 64 by SHA256(gene/species). All species have equal
total fitting weight; logistic regression uses standardized existing seven
full-profile features, C=.1 and class_weight=balanced. Candidate taxa are orders;
all share the Angiosperms parent. No order ID or mixture label is a feature.

For calibration, use fitting+calibration orders; remove own species for known
pseudoqueries and the entire true order for order-absent controls. Keep fitting
and test orders for the independent reference-test pass, excluding calibration.
For each threshold in the fixed grid 0,.01,...,1 aggregate nearest-verified hits
by species and distinct genes. Known acceptance requires >=3 genes supporting
the correct order. Absent rejection requires no wrong order reaching >=3 genes.
Choose balanced threshold by maximal mean known acceptance and absent rejection,
ties favor higher thresholds. Choose specificity threshold among those with
>=99% empirical absent-species rejection by greatest known acceptance, then
greater rejection and higher threshold. If no threshold qualifies, abstain at 1
and explicitly report failure to reach that target. Small reference calibration
counts cannot guarantee mixture precision or 99% rejection on future samples.

Freeze code, this plan, fifty reference alignments and all reused order tree /
baseline hashes before training. After reference-only threshold selection freeze
the checkpoint/thresholds before releasing all hundred prediction tasks. Score
each query exon from the existing full tree: closest order by normalized distance
is verified by the learned score. Count each supporting gene once per order.
Report original Z>0.9; balanced ML>=3 genes; balanced Z-union-ML; specificity ML;
specificity Z-union-ML. Primary new choice is specificity Z-union-ML, compared
with original Z. Preserve all outputs, including abstentions and no-improvement
results. The specificity point and quorum are fixed before mixture scoring, not
selected from mixture accuracy. No borrowed family/genus threshold is used.

Keep the truth workbook local. After all hundred predictions finish, audit and
seal before evaluation against literal Table S1 column G, with independent
openpyxl and ZIP/XML checks. Report all100, prior genus-development10, other90,
and component-size strata. Count recovered Z false negatives, newly introduced
false positives and full per-sample sets. Do not silently replace the prior
primary hierarchy100 report or change downstream search admissions. These are
retrospective results on known mixtures whose earlier truth has been inspected;
reference species can also occur in mixture sources. This is no new independent
biological/species validation. A new order output does not establish accuracy
after pruning and rerunning the downstream family/genus workflow.
