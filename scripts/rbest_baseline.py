#!/usr/bin/env python3
"""Reconstruct orthogroups from reciprocal-best chain edges alone.

This is a reference implementation, not a pipeline step. It exists to answer one
question: how good an ortholog table can you get from the pairwise alignment
evidence by itself, with no graph and no projections?

The answer matters because the shipped WF-E table disagrees with it wildly. On
the 2026-06-12 run (invocation cc39af39a106fd9e) the consensus tool labelled 21
of 5,817 orthogroups CORE-1:1 -- 0.4% -- for eight conspecific P. vivax strains,
where the same rbest edges alone give 3,979 of 5,804, or 68.6%. The difference is
that the consensus tool admits a gene's native id and its anchor-derived
projected id as two separate nodes, so nearly every group carries two "copies"
per strain and the labels, which key on max_copies, are mostly artifacts.

So this script is the acceptance test for fixing phase_e_consensus: once the
aliasing is right, the consensus table should land near this baseline, and any
remaining difference should be attributable to evidence the baseline ignores.

The rbest edges are 1:1 by construction -- on the run above, 100.0% of
(gene, target strain) lookups returned exactly one partner (24 exceptions in
154,463) -- which is what makes connected components a defensible grouping AND
what makes the clique test below meaningful.

A caveat the numbers above inherit: Sal-I contributes no rbest edges at all, so
the baseline covers 7 of the 8 strains and "all strains" in the CORE-1:1 test
means 7. The cause is a chromosome-naming mismatch, not the alignment -- Sal-I's
assembly uses GenBank accessions (CM000442.1) while its gene BED uses PlasmoDB
internal names (PVAD80_MIT), so phase_e_rbest_overlap can never intersect the two.
Its annotation was never chrom-reconciled to its assembly. Fixing that is a
data-prep job, and until it is done Sal-I's column in any ortholog table is
populated only by projections.

Usage:

    python scripts/rbest_baseline.py --edges rbest_edges.tsv \\
        [--table ortholog_table.tsv] [--out baseline.tsv]

`--table` is optional: give it the WF-E ortholog_table.tsv and the script prints
a side-by-side label comparison.
"""
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


def load_edges(path: Path):
    """Return undirected gene-gene edges as {frozenset({node_a, node_b})}.

    rbest ships both directions of most pairs. Deduplicating here is what makes
    the clique ratio a true 0..1 fraction -- counting directed edges against an
    undirected expectation of k*(k-1)/2 inflates it past 1.0.
    """
    edges = set()
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            a = f'{row["strain_a"]}#{row["gene_a"]}'
            b = f'{row["strain_b"]}#{row["gene_b"]}'
            if a != b:
                edges.add(frozenset((a, b)))
    return edges


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:      # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def label_of(n_strains, max_copies, n_all):
    """The thresholds phase_e_consensus.py uses, applied unchanged."""
    if n_strains == n_all and max_copies == 1:
        return "CORE-1:1"
    if n_strains == n_all and max_copies >= 2:
        return "CORE-VAR"
    if max_copies >= 3:
        return "FAMILY"
    if n_strains <= 2:
        return "LINEAGE-SPECIFIC"
    return "PARTIAL"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # ⚠ type=Path, and it is load-bearing: the bodies below call `.open()` on all three. A bare
    # string here is an AttributeError on every invocation, not a style preference.
    ap.add_argument("--edges", required=True, type=Path, help="rbest_edges.tsv from WF-E")
    ap.add_argument("--table", type=Path, help="WF-E ortholog_table.tsv, to compare against")
    ap.add_argument("--out", type=Path, help="write the reconstructed table here")
    ap.add_argument("--min-clique", type=float, default=0.9,
                    help="report groups below this clique completeness (default 0.9)")
    a = ap.parse_args()

    edges = load_edges(a.edges)
    uf = UnionFind()
    for e in edges:
        x, y = tuple(e)
        uf.union(x, y)

    comps = defaultdict(set)
    for node in uf.parent:
        comps[uf.find(node)].add(node)

    edge_count = Counter()
    for e in edges:
        edge_count[uf.find(next(iter(e)))] += 1

    strains = sorted({n.split("#", 1)[0] for n in uf.parent})
    n_all = len(strains)
    print(f"edges: {len(edges):,} undirected   strains: {n_all} ({', '.join(strains)})")
    print(f"orthogroups: {len(comps):,}\n")

    rows, labels, cliques, ragged = [], Counter(), [], []
    for i, (root, nodes) in enumerate(sorted(comps.items(), key=lambda kv: -len(kv[1])), 1):
        per = defaultdict(list)
        for n in nodes:
            s, g = n.split("#", 1)
            per[s].append(g)
        k, mx = len(per), max(len(v) for v in per.values())
        lab = label_of(k, mx, n_all)
        labels[lab] += 1
        expected = k * (k - 1) // 2
        clique = edge_count[root] / expected if expected else 1.0
        cliques.append(clique)
        if clique < a.min_clique:
            ragged.append((f"OG{i:06d}", k, mx, len(nodes), round(clique, 3)))
        rows.append({"orthogroup_id": f"OG{i:06d}", "label": lab, "n_strains": k,
                     "max_copies": mx, "clique": round(clique, 3),
                     **{s: ",".join(sorted(per.get(s, []))) or "-" for s in strains}})

    print("labels:")
    for lab, n in labels.most_common():
        print(f"   {lab:18} {n:>6,}  {100 * n / len(comps):5.1f}%")

    cliques.sort()
    below = sum(1 for c in cliques if c < a.min_clique)
    print("\nclique completeness (undirected edges / k*(k-1)/2):")
    print(f"   median {cliques[len(cliques) // 2]:.3f}   "
          f"10th pct {cliques[len(cliques) // 10]:.3f}")
    print(f"   below {a.min_clique}: {below:,} groups ({100 * below / len(comps):.1f}%) "
          f"-- these are chained, and a correct consensus should split them")
    for og, k, mx, n, c in sorted(ragged, key=lambda r: r[4])[:10]:
        print(f"      {og}  {k} strains, {mx} max copies, {n} genes, clique {c}")

    if a.table:
        # SIM115/PTH123: the reader was built over a bare open() whose handle was never closed.
        with a.table.open(newline="") as fh:
            shipped = Counter(r["label"] for r in csv.DictReader(fh, delimiter="\t"))
        total = sum(shipped.values())
        print("\nshipped WF-E table vs this baseline:")
        print(f"   {'label':18} {'shipped':>16} {'baseline':>16}")
        for lab in sorted(set(shipped) | set(labels)):
            print(f"   {lab:18} {shipped[lab]:>7,} {100*shipped[lab]/total:5.1f}% "
                  f"{labels[lab]:>8,} {100*labels[lab]/len(comps):5.1f}%")

    if a.out:
        with a.out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {a.out} ({len(rows):,} orthogroups)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
