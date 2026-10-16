"""Local family-level scoring for completed real-mixture predictions.

The S7 parser is kept separate from prediction. No simulated source-ID corrections
are inherited. A deferred true family is still a false negative for final calls.
"""
import csv
import json
from pathlib import Path

from .evaluate_transfer import measures, write_csv


METHODS = ("original_z", "baseline_c", "stable_c", "core", "merged", "standalone_ml", "expanded_review_union")


def prediction(sample, method):
    if method == "expanded_review_union":
        return set(sample["merged"]) | set(sample["review"])
    return set(sample[method])


def score(calls, truth, reference_families):
    by_sample = []
    if len(calls) != len(truth) or {r["sample"] for r in calls} != set(truth):
        raise ValueError("Truth/prediction sample IDs must match exactly")
    universe = set(reference_families) | {f for families in truth.values() for f in families}
    for sample in calls:
        expected = set(truth[sample["sample"]])
        if not expected:
            raise ValueError("An empty truth set requires explicit negative-control handling")
        for method in METHODS:
            chosen = prediction(sample, method)
            if not chosen <= universe:
                raise ValueError("Predictions outside the reference/truth universe")
            counts = map(len, (expected & chosen, chosen - expected, expected - chosen, universe - expected - chosen))
            by_sample.append({"sample": sample["sample"], "method": method, **measures(*counts),
                "expected": ";".join(sorted(expected)), "predicted": ";".join(sorted(chosen)),
                "false_positives": ";".join(sorted(chosen - expected)), "false_negatives": ";".join(sorted(expected - chosen))})
    summary = []
    family_rows = family_scores(calls, truth, universe)
    for method in METHODS:
        rows = [r for r in by_sample if r["method"] == method]
        present = [r for r in family_rows if r["method"] == method and r["true_sample_count"] > 0]
        summary.append({"method": method, "samples": len(rows), "universe_families": len(universe),
            **measures(*[sum(r[k] for r in rows) for k in ("TP", "FP", "FN", "TN")]),
            "exact_sets": sum(r["FP"] == r["FN"] == 0 for r in rows),
            "macro_recall_present_families": sum(r["recall"] for r in present) / len(present),
            "macro_F1_present_families": sum(r["F1"] for r in present) / len(present)})
    return by_sample, summary, universe


def family_scores(calls, truth, universe):
    rows = []
    for method in METHODS:
        for family in sorted(universe):
            counts = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
            for sample in calls:
                expected = family in truth[sample["sample"]]
                selected = family in prediction(sample, method)
                counts["TP" if expected and selected else "FN" if expected else "FP" if selected else "TN"] += 1
            rows.append({"method": method, "family": family, "true_sample_count": counts["TP"] + counts["FN"],
                         **measures(*(counts[k] for k in ("TP", "FP", "FN", "TN")))})
    return rows


def audit_exported_counts(batch, truth, universe, by_sample, summary):
    """Read plain exported lists independently from combined predictions.json."""
    batch = Path(batch)
    samples = {r["sample"]: r for r in json.loads((batch / "predictions.json").read_text())}
    for row in by_sample:
        sample = samples[row["sample"]]
        folder = batch / "samples" / sample["project"]
        if row["method"] == "expanded_review_union":
            names = set((folder / "merged_families.txt").read_text().splitlines())
            names |= set((folder / "review_families.txt").read_text().splitlines())
        else:
            names = set((folder / (row["method"] + "_families.txt")).read_text().splitlines())
        expected = set(truth[row["sample"]])
        counts = (len(expected & names), len(names - expected), len(expected - names), len(universe - expected - names))
        if counts != tuple(row[k] for k in ("TP", "FP", "FN", "TN")):
            raise ValueError("Candidate file/summary mismatch")
    for total in summary:
        rows = [r for r in by_sample if r["method"] == total["method"]]
        if any(total[k] != sum(r[k] for r in rows) for k in ("TP", "FP", "FN", "TN")):
            raise ValueError("Aggregate count mismatch")
    return {"sample_method_counts_checked": len(by_sample), "summary_rows_checked": len(summary),
            "all_exported_list_counts_match": True}
