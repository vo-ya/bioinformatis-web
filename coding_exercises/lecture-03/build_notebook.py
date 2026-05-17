"""Build exercise.ipynb for L03 — DNA Sequence Assembly: From Reads to Genomes.

Run from this directory:
    python3 build_notebook.py

Emits exercise.ipynb. Re-running overwrites the file.
"""

from __future__ import annotations

import os
# Make the shared Colab-form helper importable from the parent dir.
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from apply_colab_form import apply_colab_form  # noqa: E402

import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def md(text: str):
    """Markdown cell."""
    return new_markdown_cell(text)


def code(source: str):
    """Plain (visible) code cell."""
    return new_code_cell(source)


def hidden(source: str):
    """Code cell whose source is collapsed by default in Colab / Jupyter."""
    cell = new_code_cell(source)
    cell.metadata = {
        "jupyter": {"source_hidden": True},
        "cellView": "form",
    }
    return cell


# ---------------------------------------------------------------------------
# Cell sources
# ---------------------------------------------------------------------------


TITLE_MD = """# L03 — DNA Sequence Assembly: From Reads to Genomes

In this exercise you build a tiny de Bruijn assembler from scratch. You
simulate a 10 kb random genome, draw 500 short reads with sequencing errors,
count k-mers, build the directed graph whose nodes are (k-1)-mers and edges
are k-mers, walk it to extract linear contigs, and finally measure assembly
quality (N50, coverage) and how that quality depends on the choice of k.
"""


AHA_MD = """> **Aha.** Assembly is the overlap problem **recast as graph traversal**.
> k-mer counting is a hash-table sketch of the read set; de Bruijn edges
> stitch overlapping k-mers into a graph whose **non-branching paths are
> contigs**. Low k gives long contigs but the graph collapses at any
> repeat; high k gives a clean graph but introduces coverage gaps. The
> sweet spot is whatever maximises **N50** for your read length and depth.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 matplotlib==3.8.4 -q
"""


IMPORTS = """import math
import random
from collections import Counter, defaultdict

import numpy as np
import matplotlib.pyplot as plt

# Deterministic for the whole notebook.
SEED = 42
np.random.seed(SEED)
random.seed(SEED)

ALPHABET = "ACGT"


def random_genome(length: int, seed: int = SEED) -> str:
    \"\"\"A reproducible random ACGT string.\"\"\"
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list(ALPHABET), size=length))


def simulate_reads(genome: str, n_reads: int, read_len: int,
                   error_rate: float = 0.01, seed: int = SEED) -> list:
    \"\"\"Sample n_reads of length read_len uniformly from genome with i.i.d. base errors.\"\"\"
    rng = np.random.default_rng(seed)
    reads = []
    L = len(genome)
    starts = rng.integers(0, L - read_len + 1, size=n_reads)
    err_mask = rng.random((n_reads, read_len)) < error_rate
    alt = rng.integers(1, 4, size=(n_reads, read_len))  # 1..3 offset -> guaranteed substitution
    base_to_int = {b: i for i, b in enumerate(ALPHABET)}
    int_to_base = ALPHABET
    for r in range(n_reads):
        s = starts[r]
        seq_ints = np.array([base_to_int[c] for c in genome[s:s + read_len]])
        seq_ints = np.where(err_mask[r], (seq_ints + alt[r]) % 4, seq_ints)
        reads.append("".join(int_to_base[i] for i in seq_ints))
    return reads


print("Imports OK; numpy:", np.__version__)
"""


STEP1_MD = """## Step 1 (8 min) — Simulate the genome, draw reads, count k-mers

We need a deterministic toy problem. Generate a 10 000 bp random ACGT
genome, sample 500 reads of length 50 with a 1 % per-base substitution
error rate (so ~100× nominal coverage), then count every length-k substring
across the read set. Singleton k-mers are almost always sequencing errors —
drop them. The surviving k-mer multiset is the **sketch** of the read pile.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Simulate the genome and reads, count k-mers, drop singletons.
# ----------------------------------------------------------------------

GENOME_LEN  = 10_000
N_READS     = 500
READ_LEN    = 50
ERROR_RATE  = 0.01
K_DEFAULT   = 25


def count_kmers(reads, k):
    """Return a Counter mapping kmer -> occurrence count over all reads."""
    # TODO: slide a length-k window over every read, count occurrences.
    raise NotImplementedError


def drop_singletons(kmer_counts):
    """Return a dict of k-mers with count >= 2."""
    # TODO
    raise NotImplementedError


# Build the data and run k-mer counting.
# genome = random_genome(GENOME_LEN)
# reads  = simulate_reads(genome, N_READS, READ_LEN, ERROR_RATE)
# kmers_raw = count_kmers(reads, K_DEFAULT)
# kmers     = drop_singletons(kmers_raw)
# print(f"k={K_DEFAULT}: {len(kmers_raw):,} distinct k-mers, {len(kmers):,} survive singleton drop")
'''


CLICK_TO_EXPAND = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

GENOME_LEN  = 10_000
N_READS     = 500
READ_LEN    = 50
ERROR_RATE  = 0.01
K_DEFAULT   = 25


def count_kmers(reads, k):
    c = Counter()
    for read in reads:
        for i in range(len(read) - k + 1):
            c[read[i:i + k]] += 1
    return c


def drop_singletons(kmer_counts):
    return {km: n for km, n in kmer_counts.items() if n >= 2}


genome = random_genome(GENOME_LEN)
reads  = simulate_reads(genome, N_READS, READ_LEN, ERROR_RATE)

kmers_raw = count_kmers(reads, K_DEFAULT)
kmers     = drop_singletons(kmers_raw)

print(f"Genome length:        {len(genome):,} bp")
print(f"Reads:                {len(reads):,} x {READ_LEN} bp  "
      f"(nominal coverage = {len(reads) * READ_LEN / len(genome):.1f}x)")
print(f"k = {K_DEFAULT}")
print(f"  distinct k-mers in reads:     {len(kmers_raw):,}")
print(f"  surviving singleton drop:     {len(kmers):,}")
genome_kmers = {genome[i:i + K_DEFAULT] for i in range(len(genome) - K_DEFAULT + 1)}
print(f"  distinct k-mers in genome:    {len(genome_kmers):,}")

# Quick coverage histogram — error k-mers cluster at count=1, true k-mers higher.
counts = np.array(list(kmers_raw.values()))
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.hist(counts, bins=np.arange(0, counts.max() + 2) - 0.5, edgecolor="black")
ax.set_xlabel("k-mer occurrence count")
ax.set_ylabel("number of distinct k-mers")
ax.set_title(f"k-mer coverage spectrum (k = {K_DEFAULT})")
ax.axvline(1.5, color="red", linestyle="--", label="singleton threshold")
ax.set_yscale("log")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP2_MD = """## Step 2 (10 min) — Build the de Bruijn graph

For each surviving k-mer `w[0..k-1]`, draw a directed edge from `w[0..k-2]`
(its (k-1)-prefix) to `w[1..k-1]` (its (k-1)-suffix). The resulting graph
has at most `4^(k-1)` nodes; in practice only the (k-1)-mers that actually
appear get added. We will represent it as `dict[node, list[node]]` for
outgoing edges plus an indegree counter — that is everything Step 3 needs.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Build the de Bruijn graph from the surviving k-mers.
# ----------------------------------------------------------------------


def build_debruijn(kmers, k):
    """Return (out_adj, in_deg) where:
        out_adj[node] = list of successor (k-1)-mer nodes
        in_deg[node]  = number of incoming edges
    """
    # TODO: for each k-mer, split into prefix + suffix and add an edge.
    raise NotImplementedError


# out_adj, in_deg = build_debruijn(kmers, K_DEFAULT)
# print(f"de Bruijn graph: {len(out_adj):,} nodes, "
#       f"{sum(len(v) for v in out_adj.values()):,} edges")
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.


def build_debruijn(kmers, k):
    out_adj = defaultdict(list)
    in_deg  = defaultdict(int)
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        out_adj[prefix].append(suffix)
        in_deg[suffix] += 1
        # Make sure the suffix node exists in out_adj (even if it has no successors yet).
        out_adj.setdefault(suffix, [])
        in_deg.setdefault(prefix, 0)
    return dict(out_adj), dict(in_deg)


out_adj, in_deg = build_debruijn(kmers, K_DEFAULT)

n_nodes = len(out_adj)
n_edges = sum(len(v) for v in out_adj.values())
n_branching = sum(1 for v in out_adj.values() if len(v) > 1)

print(f"de Bruijn graph (k = {K_DEFAULT}):")
print(f"  nodes (distinct (k-1)-mers): {n_nodes:,}")
print(f"  edges (distinct k-mers):     {n_edges:,}")
print(f"  branching nodes (out > 1):   {n_branching:,}")

# Degree spectrum.
out_deg_counts = Counter(len(v) for v in out_adj.values())
in_deg_counts  = Counter(in_deg.values())
print("  out-degree distribution:", dict(sorted(out_deg_counts.items())))
print("  in-degree  distribution:", dict(sorted(in_deg_counts.items())))
'''


STEP3_MD = """## Step 3 (12 min) — Greedy walk: extract linear contigs

A **non-branching path** is a maximal walk through nodes whose in- and
out-degree are both exactly 1, anchored at a node that is either a source
(in-degree ≠ 1 or out-degree > 1) or has no incoming edges. Each such path
spells out one contig: start with the first node's (k-1)-mer and append the
last base of every successor edge. Halt at branch points (out-degree ≠ 1)
or at dead ends.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Walk every non-branching path; emit one contig per walk.
# ----------------------------------------------------------------------


def is_branching(node, out_adj, in_deg):
    """True if node has in-deg != 1 or out-deg != 1 (a contig start/end)."""
    # TODO
    raise NotImplementedError


def extract_contigs(out_adj, in_deg, k):
    """Return a list of contig strings.

    Walk every edge once. Start a new contig at any branching node that has
    outgoing edges; for cycles of unbranching nodes, walk once around.
    """
    # TODO
    raise NotImplementedError


# contigs = extract_contigs(out_adj, in_deg, K_DEFAULT)
# contigs.sort(key=len, reverse=True)
# print(f"emitted {len(contigs)} contigs;  longest = {len(contigs[0])} bp;  "
#       f"shortest = {len(contigs[-1])} bp")
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.


def is_branching(node, out_adj, in_deg):
    return in_deg.get(node, 0) != 1 or len(out_adj.get(node, [])) != 1


def extract_contigs(out_adj, in_deg, k):
    contigs = []
    # Use multi-edge counts to consume each edge exactly once.
    edge_iter = {n: list(succs) for n, succs in out_adj.items()}

    # Pass 1: start at every branching node that has outgoing edges.
    for node in list(out_adj):
        if is_branching(node, out_adj, in_deg) and edge_iter.get(node):
            for _ in range(len(edge_iter[node])):
                # Pop one outgoing edge and walk forward until we hit a branch / dead end.
                nxt = edge_iter[node].pop()
                path_chars = [node, nxt[-1]]
                cur = nxt
                while not is_branching(cur, out_adj, in_deg) and edge_iter.get(cur):
                    nxt = edge_iter[cur].pop()
                    path_chars.append(nxt[-1])
                    cur = nxt
                contigs.append(path_chars[0] + "".join(path_chars[1:]))

    # Pass 2: anything left over lives on isolated cycles. Walk them once around.
    for node in list(edge_iter):
        while edge_iter.get(node):
            nxt = edge_iter[node].pop()
            path_chars = [node, nxt[-1]]
            cur = nxt
            steps = 0
            while edge_iter.get(cur) and cur != node and steps < len(out_adj) + 5:
                nxt = edge_iter[cur].pop()
                path_chars.append(nxt[-1])
                cur = nxt
                steps += 1
            contigs.append(path_chars[0] + "".join(path_chars[1:]))

    return contigs


contigs = extract_contigs(out_adj, in_deg, K_DEFAULT)
contigs.sort(key=len, reverse=True)
print(f"Extracted {len(contigs)} contigs")
print(f"  longest:  {len(contigs[0]):,} bp")
print(f"  median:   {int(np.median([len(c) for c in contigs])):,} bp")
print(f"  shortest: {len(contigs[-1]):,} bp")
print(f"  total bp: {sum(len(c) for c in contigs):,}")

# Show the top-5 contig lengths.
print("\\nTop-5 contig lengths:",
      [len(c) for c in contigs[:5]])
'''


STEP4_MD = """## Step 4 (15 min) — Map contigs back to the genome; compute N50 + coverage

To grade the assembly we need ground truth — and we have it, because we
simulated the genome. For each contig, locate its position on the genome
by a **naive substring search** (it is only a 10 kb haystack, so this is
fine). Then compute the assembly's headline metrics:

- **N50**: sort contigs by length descending; the N50 is the length of the
  contig at which the **cumulative** length first reaches half of the total
  assembled bp.
- **Coverage fraction**: the fraction of genome positions covered by at
  least one contig hit.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Align contigs to the genome; compute N50 + coverage.
# ----------------------------------------------------------------------


def find_first_occurrence(haystack, needle):
    """Return the 0-based start position of needle in haystack, or -1."""
    # TODO: builtin .find(...) is fine here.
    raise NotImplementedError


def assembly_n50(lengths):
    """Standard N50: smallest L such that contigs >= L cover >= 50 % of total bp."""
    # TODO
    raise NotImplementedError


def coverage_fraction(contigs, genome):
    """Fraction of genome positions covered by at least one mapped contig."""
    # TODO: walk each contig, find its hit, mark a boolean mask.
    raise NotImplementedError


# n50 = assembly_n50([len(c) for c in contigs])
# cov = coverage_fraction(contigs, genome)
# print(f"N50 = {n50:,} bp;  coverage = {cov:.1%}")
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.


def find_first_occurrence(haystack, needle):
    return haystack.find(needle)


def assembly_n50(lengths):
    if not lengths:
        return 0
    lens = sorted(lengths, reverse=True)
    total = sum(lens)
    cum = 0
    for L in lens:
        cum += L
        if cum >= total / 2:
            return L
    return lens[-1]


def coverage_fraction(contigs, genome):
    mask = np.zeros(len(genome), dtype=bool)
    for c in contigs:
        if len(c) < 10:  # ignore microscopic contigs
            continue
        pos = find_first_occurrence(genome, c)
        if pos >= 0:
            mask[pos:pos + len(c)] = True
    return mask.mean()


lengths = [len(c) for c in contigs]
n50 = assembly_n50(lengths)
cov = coverage_fraction(contigs, genome)
total_bp = sum(lengths)

print(f"Assembly metrics for k = {K_DEFAULT}:")
print(f"  contigs:            {len(contigs):,}")
print(f"  total assembled bp: {total_bp:,}")
print(f"  longest contig:     {max(lengths):,} bp")
print(f"  N50:                {n50:,} bp")
print(f"  coverage of truth:  {cov:.1%}")

# How many contigs map at all?
mapped = sum(1 for c in contigs if len(c) >= 10 and find_first_occurrence(genome, c) >= 0)
print(f"  mapped contigs:     {mapped} / {len(contigs)}")

# Sanity: cumulative-length curve.
sorted_lens = sorted(lengths, reverse=True)
cumlen = np.cumsum(sorted_lens)
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(np.arange(1, len(sorted_lens) + 1), cumlen / cumlen[-1], marker="o", ms=3)
ax.axhline(0.5, color="red", linestyle="--", label=f"N50 = {n50:,} bp")
ax.set_xlabel("rank (longest first)")
ax.set_ylabel("cumulative fraction of assembled bp")
ax.set_title(f"N50 curve (k = {K_DEFAULT})")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP5_MD = """## Step 5 (15 min) — Sweep k ∈ {15, 25, 35}; plot N50 vs k

The choice of k trades off two failure modes:

- **k too small:** every short genomic repeat collapses to the same node;
  the graph branches; contigs shatter at every repeat boundary. Short
  contigs, low N50.
- **k too large:** any k-mer that contains an error is unique and gets
  dropped; the graph develops gaps; contigs end at every error. Many
  short contigs again.

Run the full pipeline at k ∈ {15, 25, 35} and plot N50 + coverage to see
the sweet spot.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Sweep k; plot N50 + coverage vs k.
# ----------------------------------------------------------------------


def assemble(reads, genome, k):
    """Glue together Steps 1-4. Return dict(k, n_contigs, longest, n50, coverage)."""
    # TODO
    raise NotImplementedError


# K_SWEEP = [15, 25, 35]
# results = [assemble(reads, genome, k) for k in K_SWEEP]
# for r in results:
#     print(r)
'''

STEP5_SOLUTION = '''# Reference solution — Step 5.


def assemble(reads, genome, k):
    raw = count_kmers(reads, k)
    surv = drop_singletons(raw)
    if not surv:
        return {"k": k, "n_contigs": 0, "longest": 0, "n50": 0, "coverage": 0.0}
    oa, ind = build_debruijn(surv, k)
    cs = extract_contigs(oa, ind, k)
    lens = [len(c) for c in cs]
    return {
        "k":         k,
        "n_contigs": len(cs),
        "longest":   max(lens) if lens else 0,
        "n50":       assembly_n50(lens),
        "coverage":  coverage_fraction(cs, genome),
    }


K_SWEEP = [15, 25, 35]
results = [assemble(reads, genome, k) for k in K_SWEEP]

print(f"{'k':>4s}  {'n_contigs':>10s}  {'longest':>9s}  {'N50':>8s}  {'coverage':>9s}")
for r in results:
    print(f"{r['k']:>4d}  {r['n_contigs']:>10d}  {r['longest']:>9d}  "
          f"{r['n50']:>8d}  {r['coverage']:>9.1%}")

ks = [r["k"] for r in results]
n50s = [r["n50"] for r in results]
covs = [r["coverage"] for r in results]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
ax1.plot(ks, n50s, marker="o", ms=8, lw=2, color="steelblue")
ax1.set_xlabel("k")
ax1.set_ylabel("N50 (bp)")
ax1.set_title("N50 vs k")
ax1.grid(alpha=0.3)

ax2.plot(ks, covs, marker="o", ms=8, lw=2, color="darkorange")
ax2.set_xlabel("k")
ax2.set_ylabel("coverage fraction")
ax2.set_title("Coverage vs k")
ax2.set_ylim(0, 1.05)
ax2.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# Find the best k.
best = max(results, key=lambda r: r["n50"])
print(f"\\nBest k by N50: k = {best['k']}, N50 = {best['n50']:,} bp, coverage = {best['coverage']:.1%}")
'''


SELFCHECK_MD = """## Self-check

These asserts pin down the load-bearing numerical pieces of the pipeline.
If you ran the reference solutions above they should pass; if you wrote
your own and an assert fires, revisit the corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. count_kmers + drop_singletons sanity on a tiny worked example.
toy_reads = ["AAACCCG", "AACCCGT", "ACCCGTA"]
toy_counts = count_kmers(toy_reads, k=3)
# k-mer "CCC" appears in all three reads, so its count must be 3.
assert toy_counts["CCC"] == 3, f"CCC count expected 3, got {toy_counts['CCC']}"
# k-mer "AAA" appears once (only in the first read).
assert toy_counts["AAA"] == 1, f"AAA count expected 1, got {toy_counts['AAA']}"
toy_surv = drop_singletons(toy_counts)
assert "AAA" not in toy_surv, "singleton AAA should have been dropped"

# 2. build_debruijn on a deterministic toy: 4-mers from "ACGTAC" with k=4 -> 3 edges.
toy_kmers = {"ACGT": 1, "CGTA": 1, "GTAC": 1}
oa, ind = build_debruijn(toy_kmers, k=4)
# The three edges are ACG -> CGT, CGT -> GTA, GTA -> TAC.
assert "ACG" in oa and "CGT" in oa["ACG"], "edge ACG -> CGT missing"
assert "CGT" in oa and "GTA" in oa["CGT"], "edge CGT -> GTA missing"
assert "GTA" in oa and "TAC" in oa["GTA"], "edge GTA -> TAC missing"

# 3. extract_contigs on the same toy reconstructs the full 6-mer.
toy_contigs = extract_contigs(oa, ind, k=4)
assert "ACGTAC" in toy_contigs, f"expected ACGTAC, got {toy_contigs}"

# 4. assembly_n50 sanity: lengths [10, 5, 3, 2, 1] (sum=21, half=10.5)
#    cumulative [10, 15, ...] -> N50 is the contig at which we cross 10.5, length 5.
assert assembly_n50([10, 5, 3, 2, 1]) == 5

# 5. End-to-end pipeline: at k=25 with 50 bp reads and 1% error, the maximum
#    possible contig length is bounded by read length and how often error-free
#    overlaps chain. We should still recover at least 30% of the genome and
#    produce contigs longer than the read length somewhere.
final_n50 = assembly_n50([len(c) for c in contigs])
final_cov = coverage_fraction(contigs, genome)
longest = max(len(c) for c in contigs) if contigs else 0
print(f"Final assembly:  N50 = {final_n50:,} bp,  longest = {longest:,} bp,  coverage = {final_cov:.1%}")
assert final_n50 >= 30, f"N50 = {final_n50} suspiciously low for k = {K_DEFAULT}"
assert longest >= READ_LEN, f"longest contig = {longest} is below read length"
assert final_cov >= 0.30, f"coverage = {final_cov:.1%} suspiciously low"

print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — hash sketch, graph traversal, trellis decoding

The pipeline you just built is three classic EE objects in costume:

1. **k-mer counting = hash-table sketch.** A read pile of 500 × 50 bp is
   25 000 sliding windows; you compress it to a count vector indexed by
   k-mers. Singleton-dropping is **denoising** — a per-bin amplitude
   threshold on the sketch.
2. **Contig extraction = path-finding on a DAG.** Non-branching paths are
   maximum-likelihood walks through a trellis whose nodes are (k-1)-mer
   states and whose edges are observed k-mer transitions. The hard step
   in real assemblers is **disambiguating repeats** — same problem as
   Viterbi decoding through a trellis when two paths share a stretch of
   identical state.
3. **k = the receptive-field knob.** Small k = small receptive field, the
   detector confuses repeats (model bias, low resolution). Large k =
   large receptive field, the detector is sensitive to **single-base
   errors** (model variance, no robustness). N50 vs k is your standard
   bias-variance curve, and you found its peak empirically.

This is exactly why long reads (L11) change the game: a 10 kb HiFi read
**is** a contig, no graph traversal needed for any unique region. The
graph-vs-overlap choice in real assemblers is essentially a sample-rate
decision: short reads sample the genome densely (de Bruijn wins), long
reads sample it sparsely with high SNR (overlap-layout-consensus wins).
"""


# ---------------------------------------------------------------------------
# Assemble + write
# ---------------------------------------------------------------------------


def build():
    nb = new_notebook()
    nb.cells = [
        md(TITLE_MD),
        md(AHA_MD),
        code(PREAMBLE),
        code(IMPORTS),

        md(STEP1_MD),
        code(STEP1_TODO),
        md(CLICK_TO_EXPAND),
        hidden(STEP1_SOLUTION),

        md(STEP2_MD),
        code(STEP2_TODO),
        md(CLICK_TO_EXPAND),
        hidden(STEP2_SOLUTION),

        md(STEP3_MD),
        code(STEP3_TODO),
        md(CLICK_TO_EXPAND),
        hidden(STEP3_SOLUTION),

        md(STEP4_MD),
        code(STEP4_TODO),
        md(CLICK_TO_EXPAND),
        hidden(STEP4_SOLUTION),

        md(STEP5_MD),
        code(STEP5_TODO),
        md(CLICK_TO_EXPAND),
        hidden(STEP5_SOLUTION),

        md(SELFCHECK_MD),
        code(SELFCHECK),

        md(EE_MD),
    ]
    nb.metadata = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.11"},
        "colab": {"provenance": [], "toc_visible": True},
    }
    return nb


if __name__ == "__main__":
    nb = build()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "exercise.ipynb")
    with open(out_path, "w", encoding="utf-8") as f:
        apply_colab_form(nb)
        nbformat.write(nb, f)
    nb2 = nbformat.read(out_path, as_version=4)
    n_md = sum(1 for c in nb2.cells if c.cell_type == "markdown")
    n_code = sum(1 for c in nb2.cells if c.cell_type == "code")
    n_hidden = sum(
        1 for c in nb2.cells
        if c.cell_type == "code"
        and c.metadata.get("jupyter", {}).get("source_hidden")
    )
    print(f"Wrote {out_path}")
    print(f"  cells: {len(nb2.cells)} total  ({n_md} md, {n_code} code, {n_hidden} hidden-solution)")
