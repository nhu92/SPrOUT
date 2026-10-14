"""Compact, explicitly pruned phylogenetic context for portable PEA reports.

No reads or sequences are written to the report. The scored evidence is derived
from SPrOUT's cleaned query-to-reference distance matrices and trimmed exon
alignments; the topology is copied from the corresponding FastTree/IQ-TREE file.
"""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import median

from Bio import Phylo, SeqIO


TREE_NAME = re.compile(r"^(.+)_exon_(\d+)\.tre$")


def _usable_bp(sequence: str) -> int:
    return sum(base in "ACGT" for base in sequence.upper())


def _alignment_lengths(path: Path, project: str) -> tuple[dict[str, int], float | None, int | None]:
    if not path.is_file():
        return {}, None, None
    with path.open(encoding="utf-8") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    lengths = {record.id: _usable_bp(str(record.seq)) for record in records}
    references = [length for name, length in lengths.items() if project not in name and length]
    return lengths, median(references) if references else None, len(records[0].seq) if records else None


def _matrix_support(
    matrix: Path, project: str, level: str, lengths: dict[str, int],
    reference_bp: float | None, taxon_from_leaf,
) -> tuple[list[dict], int]:
    if not matrix.is_file():
        return [], 0
    scores: dict[str, float] = defaultdict(float)
    bp_total = 0
    with matrix.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(fields) < 2:
            return [], 0
        rows = [row for row in reader if project not in (row.get(fields[0]) or "")]
        for query in (field for field in fields[1:] if project in field):
            query_bp = lengths.get(query, lengths.get(query.removeprefix("_R_"), 0))
            bp_total += query_bp
            if not query_bp or not reference_bp:
                continue
            distances = []
            for row in rows:
                taxon = taxon_from_leaf(row.get(fields[0]) or "", level)
                try:
                    distance = float(row.get(query) or "nan")
                except ValueError:
                    continue
                if taxon and math.isfinite(distance) and distance >= 0:
                    distances.append((distance, taxon))
            if not distances:
                continue
            best = min(distance for distance, _ in distances)
            winners = sorted({taxon for distance, taxon in distances if math.isclose(
                distance, best, rel_tol=1e-9, abs_tol=1e-12,
            )})
            # Coverage is bounded at one: extra-long contigs cannot dominate.
            # Similarity is the same monotone 1/(1+d) transform used by SPrOUT.
            score = min(1.0, query_bp / reference_bp) / (1.0 + best) / len(winners)
            for taxon in winners:
                scores[taxon] += score
    return [
        {"taxon": taxon, "score": round(score, 6)}
        for taxon, score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    ], bp_total


def _parents(root) -> dict:
    parents = {}
    stack = [root]
    while stack:
        node = stack.pop()
        for child in node.clades:
            parents[child] = node
            stack.append(child)
    return parents


def _nearest_sibling(tip, parents: dict, keep: set[str], classify) -> str | None:
    original = classify(tip.name or "")
    current = tip
    while current in parents:
        parent = parents[current]
        for sibling in parent.clades:
            if sibling is current:
                continue
            for leaf in sibling.get_terminals():
                if (leaf.name and leaf.name not in keep and classify(leaf.name)
                        and classify(leaf.name) != original):
                    return leaf.name
        current = parent
    return None


def _serialize_pruned(node, keep: set[str], classify, query_name: str | None = None):
    length = round(max(0.0, node.branch_length or 0.0), 8)
    if node.is_terminal():
        if node.name not in keep:
            return None
        return {
            "l": length, "n": node.name or "", "t": classify(node.name or ""),
            "q": bool(query_name and query_name in (node.name or "")),
        }
    children = [item for child in node.clades if (item := _serialize_pruned(
        child, keep, classify, query_name,
    )) is not None]
    if not children:
        return None
    if len(children) == 1:
        children[0]["l"] = round(children[0]["l"] + length, 8)
        return children[0]
    return {"l": length, "c": children}


def _exon_context(tree, project: str, level: str, candidates: list[str],
                  taxon_from_leaf, max_tips: int) -> dict:
    tips = tree.get_terminals()
    query = [tip for tip in tips if tip.name and project in tip.name]
    keep = {tip.name for tip in query[:max_tips]}
    classify = lambda name: taxon_from_leaf(name, level)
    by_taxon = defaultdict(list)
    for tip in tips:
        if tip.name and project not in tip.name:
            taxon = classify(tip.name)
            if taxon:
                by_taxon[taxon].append(tip)
    anchors = []
    for taxon in candidates:
        choices = by_taxon.get(taxon, [])
        if choices and len(keep) < max_tips:
            tip = choices[0]
            keep.add(tip.name)
            anchors.append(tip)
    parents = _parents(tree.root)
    for tip in query[:12] + anchors[:30]:
        if len(keep) >= max_tips:
            break
        neighbor = _nearest_sibling(tip, parents, keep, classify)
        if neighbor:
            keep.add(neighbor)
    if len(keep) < 2:
        for tip in tips:
            if tip.name and classify(tip.name):
                keep.add(tip.name)
            if len(keep) >= min(2, max_tips):
                break
    return {
        "tree": _serialize_pruned(tree.root, keep, classify, project),
        "source_tips": len(tips), "displayed_tips": len(keep),
        "omitted_queries": max(0, len(query) - max_tips),
    }


def export_exon_trees(phylo_dir: Path, matrix_dir: Path | None, project: str,
                      level: str, predictions: list[dict], taxon_from_leaf,
                      max_tips: int = 80) -> dict[str, list[dict]]:
    if max_tips < 2:
        raise ValueError("max_tips must be at least 2")
    candidates = list(dict.fromkeys(
        [row["taxon"] for row in predictions if row.get("pipeline_selected")]
        + [row["taxon"] for row in predictions[:40]]
    ))
    by_gene = defaultdict(list)
    trees = defaultdict(list)
    for path in sorted(phylo_dir.glob("*.tre")):
        match = TREE_NAME.fullmatch(path.name)
        if match:
            trees[match.group(1)].append(path)
    for gene, files in sorted(trees.items()):
        for matrix_index, path in enumerate(files, 1):
            try:
                tree = Phylo.read(path, "newick")
            except ValueError as exc:
                if "There are no trees in this file" in str(exc):
                    continue
                raise
            alignment = phylo_dir / path.name.replace(".tre", "_trimmed.fasta")
            lengths, reference_bp, alignment_bp = _alignment_lengths(alignment, project)
            matrix = matrix_dir / f"{gene}.{matrix_index}.cleaned.csv" if matrix_dir else Path()
            support, query_bp = _matrix_support(
                matrix, project, level, lengths, reference_bp, taxon_from_leaf,
            ) if matrix_dir else ([], 0)
            context = _exon_context(
                tree, project, level, candidates, taxon_from_leaf, max_tips,
            )
            by_gene[gene].append({
                "id": path.stem, "alignment_bp": alignment_bp,
                "reference_median_bp": round(reference_bp, 1) if reference_bp else None,
                "query_usable_bp": query_bp, "support": support,
                **context,
            })
    return dict(by_gene)


def export_summary_tree(context_path: Path, level: str, predictions: list[dict],
                        taxon_from_leaf, max_tips: int = 80) -> dict:
    """A Kew reference topology, using one representative species per rank.

    This does not infer a mixture/species tree. Missing or non-monophyletic
    groups are deliberately not claimed to have a single clade boundary.
    """
    tree = Phylo.read(context_path, "newick")
    classify = lambda name: taxon_from_leaf(name, level)
    by_taxon = {}
    for tip in tree.get_terminals():
        if tip.name:
            taxon = classify(tip.name)
            if taxon:
                by_taxon.setdefault(taxon, tip)
    requested = list(dict.fromkeys(
        [row["taxon"] for row in predictions if row.get("pipeline_selected")]
        + [row["taxon"] for row in predictions[:40]]
    ))
    keep = set()
    anchors = []
    for taxon in requested:
        tip = by_taxon.get(taxon)
        if tip and len(keep) < max_tips - 1:
            keep.add(tip.name)
            anchors.append(tip)
    parents = _parents(tree.root)
    for tip in anchors[:30]:
        if len(keep) >= max_tips:
            break
        neighbor = _nearest_sibling(tip, parents, keep, classify)
        if neighbor:
            keep.add(neighbor)
    represented = {classify(name) for name in keep}
    return {
        "source": "Kew Tree of Life species-tree topology; one representative specimen per ranked taxon",
        "tree": _serialize_pruned(tree.root, keep, classify),
        "displayed_tips": len(keep), "requested_taxa": len(requested),
        "represented_taxa": sorted(represented & set(requested)),
        "missing_taxa": [taxon for taxon in requested if taxon not in represented],
    }
