"""Tree neighborhoods are search regions, never labels for unseen genera.

All recoverable placements of a polyphyletic genus are retained. Removing one
cell does not remove that genus when another placement remains active.
"""
from dataclasses import dataclass
import hashlib
from collections import defaultdict

import numpy as np
from Bio import Phylo


@dataclass(frozen=True)
class Cell:
    family: str
    tips: tuple[int, ...]

    @property
    def key(self):
        return hashlib.sha256((self.family + ':' + ','.join(map(str, self.tips))).encode()).hexdigest()[:16]


class TreeIndex:
    def __init__(self, path):
        tree = Phylo.read(path, 'newick')
        self.nodes = list(tree.find_clades(order='preorder'))
        index = {id(node): i for i, node in enumerate(self.nodes)}
        self.parent = [-1] * len(self.nodes)
        self.children = [[] for _ in self.nodes]
        self.edges = [[] for _ in self.nodes]
        self.depth = [0] * len(self.nodes)
        self.taxonomy = {}
        self.by_genus = defaultdict(list)
        self.by_family = defaultdict(list)
        self.by_name = {}
        supplied_lengths = 0
        for i, node in enumerate(self.nodes):
            for child in node.clades:
                j = index[id(child)]
                length = 1.0 if child.branch_length is None else float(child.branch_length)
                if not np.isfinite(length) or length < 0:
                    raise ValueError('Invalid Kew tree edge')
                supplied_lengths += child.branch_length is not None
                self.parent[j] = i
                self.depth[j] = self.depth[i] + 1
                self.children[i].append(j)
                self.edges[i].append((j, length))
                self.edges[j].append((i, length))
            if node.is_terminal():
                fields = str(node.name or '').split('_', 4)
                if len(fields) == 5:
                    self.taxonomy[i] = tuple(fields)
                    self.by_genus[fields[2]].append(i)
                    self.by_family[fields[1]].append(i)
                    self.by_name[node.name] = i
        self.distance_kind = 'topological_edges' if not supplied_lengths else 'patristic_with_unit_missing_edges'
        self._distance_cache = {}

    def distances(self, source):
        if source not in self._distance_cache:
            if len(self._distance_cache) >= 256:
                self._distance_cache.clear()
            result = np.zeros(len(self.nodes))
            stack = [(source, -1)]
            while stack:
                node, previous = stack.pop()
                for child, length in self.edges[node]:
                    if child != previous:
                        result[child] = result[node] + length
                        stack.append((child, node))
            self._distance_cache[source] = result
        return self._distance_cache[source]

    def lca(self, tips):
        left = tips[0]
        for right in tips[1:]:
            a, b = left, right
            while a != b:
                if self.depth[a] >= self.depth[b]:
                    a = self.parent[a]
                else:
                    b = self.parent[b]
            left = a
        return left

    def genera(self, cell):
        return {self.taxonomy[t][2] for t in cell.tips}

    def split(self, cell):
        if len(cell.tips) < 2 or len(self.genera(cell)) < 2:
            return None
        root = self.lca(cell.tips)
        groups = defaultdict(list)
        for tip in cell.tips:
            node = tip
            while self.parent[node] != root:
                node = self.parent[node]
            groups[node].append(tip)
        if len(groups) < 2:
            return None
        # A polytomy is represented as a forest; never invent a resolved edge.
        bins, weights = [[], []], [0, 0]
        for group in sorted(groups.values(), key=lambda g: (-len({self.taxonomy[t][2] for t in g}), min(g))):
            side = 0 if weights[0] <= weights[1] else 1
            bins[side].extend(group)
            weights[side] += len({self.taxonomy[t][2] for t in group})
        return tuple(Cell(cell.family, tuple(sorted(group))) for group in bins)

    def partition(self, cells, budget=32, target_genera=2):
        cells = list(cells)
        if len(cells) > budget:
            raise ValueError('Budget cannot represent all active families/cells')
        while len(cells) < budget:
            options = [(len(self.genera(c)), len(c.tips), c.key, i) for i, c in enumerate(cells)
                       if len(self.genera(c)) > target_genera and self.split(c) is not None]
            if not options:
                break
            _, _, _, i = max(options)
            cells[i:i+1] = self.split(cells[i])
        return sorted(cells, key=lambda c: (c.family, c.key))

    def medoid_order(self, cell, eligible=None):
        """Genus-balanced minimax centers, then mean distance; not alphabetic ends.

        A farthest-point cover of at most 32 landmarks bounds the planning cost.
        The full cell, not the downloaded anchors, remains the search hypothesis.
        """
        candidates = sorted(set(cell.tips) & set(eligible if eligible is not None else cell.tips))
        if not candidates:
            return []
        tips = np.array(cell.tips)
        genus_ids = {}
        for t in tips:
            genus_ids.setdefault(self.taxonomy[int(t)][2], []).append(int(t))
        # One deterministic landmark per genus avoids weighting species-rich genera.
        pool = [min(v) for _, v in sorted(genus_ids.items())]
        landmarks = [pool[0]]
        nearest = self.distances(pool[0])[pool].copy()
        for _ in range(min(31, len(pool)-1)):
            next_tip = pool[int(np.argmax(nearest))]
            if next_tip in landmarks:
                break
            landmarks.append(next_tip)
            nearest = np.minimum(nearest, self.distances(next_tip)[pool])
        distances = np.array([self.distances(t)[candidates] for t in landmarks])
        ranked = sorted(range(len(candidates)), key=lambda i: (float(distances[:, i].max()),
                         float(distances[:, i].mean()), self.nodes[candidates[i]].name))
        return [candidates[i] for i in ranked]

    def neighborhood_guard(self, previous, chosen, anchor_tips, neighbors=4):
        """Add only nearby genera, rather than every genus in a sister branch."""
        if neighbors<=0:
            return list(chosen)
        universe = {t for cell in previous for t in cell.tips}
        kept = {t for cell in chosen for t in cell.tips}
        additions = set()
        for anchor in anchor_tips:
            family = self.taxonomy[anchor][1]
            present = {self.taxonomy[t][2] for t in kept}
            ranked = sorted((t for t in universe if self.taxonomy[t][1] == family
                             and self.taxonomy[t][2] not in present),
                            key=lambda t: (self.distances(anchor)[t], self.nodes[t].name))
            genera = []
            for tip in ranked:
                genus = self.taxonomy[tip][2]
                if genus not in genera:
                    genera.append(genus)
                if len(genera) >= neighbors:
                    break
            additions.update(t for t in universe if self.taxonomy[t][1] == family
                             and self.taxonomy[t][2] in genera)
        by_family = defaultdict(list)
        for tip in additions - kept:
            by_family[self.taxonomy[tip][1]].append(tip)
        return [*chosen, *(Cell(family, tuple(sorted(tips))) for family, tips in sorted(by_family.items()))]

    def representative_regions(self, cell):
        if len(self.genera(cell))<2:
            return [cell]
        first = self.medoid_order(cell)[0]
        different = [t for t in cell.tips if self.taxonomy[t][2]!=self.taxonomy[first][2]]
        second = max(different,key=lambda t:(self.distances(first)[t],self.nodes[t].name))
        halves = [[],[]]
        for tip in cell.tips:
            side = 0 if self.distances(first)[tip]<=self.distances(second)[tip] else 1
            halves[side].append(tip)
        return [Cell(cell.family,tuple(sorted(tips))) for tips in halves if tips]
