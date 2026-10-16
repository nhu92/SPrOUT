"""Render family-presence errors and fixed-method metrics from real11 exports."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np


def plot(folder):
    folder = Path(folder)
    calls = json.loads((folder / "predictions.json").read_text())
    truth = json.loads((folder / "truth.json").read_text())
    families = sorted({f for row in calls for f in row["merged"]} | {f for values in truth.values() for f in values})
    values = np.zeros((len(families), len(calls)), dtype=int)
    for j, sample in enumerate(calls):
        for i, family in enumerate(families):
            present = family in truth[sample["sample"]]
            chosen = family in sample["merged"]
            values[i, j] = 1 if present and chosen else 2 if present else 3 if chosen else 0
    with (folder / "summary.csv").open(encoding="utf-8-sig", newline="") as f:
        summary = list(csv.DictReader(f))
    names = {"original_z": "Original Z", "baseline_c": "Complete C", "stable_c": "Stable C",
             "core": "C + ML core", "merged": "Merged (primary)", "standalone_ml": "ML alone",
             "expanded_review_union": "Merged + all review"}
    fig, axes = plt.subplots(1, 2, figsize=(13.5, max(6.5, len(families) * .26)),
                             gridspec_kw={"width_ratios": [1.1, 1]}, layout="constrained")
    colors = ["#f0f1f3", "#26846b", "#dd8a2b", "#b44766"]
    axes[0].imshow(values, aspect="auto", cmap=ListedColormap(colors), norm=BoundaryNorm(np.arange(-.5, 4), 4))
    axes[0].set_xticks(range(len(calls)), [s["sample"].replace("EMseq", "") for s in calls], rotation=0)
    axes[0].set_yticks(range(len(families)), families, fontsize=8)
    axes[0].set_xlabel("Real mixture (EMseq)")
    axes[0].set_title("Frozen merged calls versus S7 families", loc="left", fontsize=11)
    from matplotlib.patches import Patch
    axes[0].legend([Patch(facecolor=c) for c in colors[1:]], ["TP", "FN", "FP"],
                   loc="upper center", bbox_to_anchor=(.5, -.08), ncol=3, frameon=False)
    y = np.arange(len(summary))
    for offset, metric, color in ((-.23, "precision", "#286b9b"), (0, "recall", "#dc9b32"), (.23, "F1", "#26846b")):
        axes[1].barh(y + offset, [float(r[metric]) for r in summary], height=.21, label=metric, color=color)
    axes[1].set_yticks(y, [names[r["method"]] for r in summary], fontsize=9)
    axes[1].invert_yaxis()
    axes[1].set_xlim(0, 1.04)
    axes[1].set_xlabel("Micro score across all 11 mixtures")
    axes[1].set_title("Predeclared methods, identical samples", loc="left", fontsize=11)
    axes[1].legend(loc="upper center", bbox_to_anchor=(.5, -.08), ncol=3, frameon=False)
    axes[1].grid(axis="x", alpha=.15)
    axes[1].set_axisbelow(True)
    fig.savefig(folder / "real11_validation.png", dpi=180)
    svg = folder / "real11_validation.svg"
    fig.savefig(svg)
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines()) + "\n", encoding="utf-8")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True)
    plot(parser.parse_args().results)
