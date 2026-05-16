"""Build the L11 exercise notebook — Long Reads and the Pangenome.

Run:  python3 build_notebook.py
Emits: exercise.ipynb in the same directory.

The notebook is structured as:
  1. Title markdown
  2. Aha callout markdown
  3. Preamble (!pip install)
  4. Imports + seeds
  5..9. Step N markdown + TODO cell + hidden solution cell
  10. Final self-check assert cell
  11. EE framing markdown

Solution cells use metadata.jupyter.source_hidden = True so they collapse
in Colab. Each solution cell is preceded by a short prompt line that tells
the student how to expand it.

All data is generated in-notebook from deterministic NumPy seeds — no
external fetches, per the L11 spec.
"""

from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell


HERE = Path(__file__).parent
OUT = HERE / "exercise.ipynb"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def md(text: str):
    return new_markdown_cell(text)


def code(text: str):
    return new_code_cell(text)


def hidden_solution(text: str):
    """Code cell flagged source_hidden in Jupyter metadata.

    Colab respects `metadata.jupyter.source_hidden`; classic Jupyter does too.
    The cell still executes; only the source is collapsed.
    """
    cell = new_code_cell(text)
    cell.metadata = {
        "jupyter": {"source_hidden": True},
        "cellView": "form",
    }
    return cell


# ---------------------------------------------------------------------------
# Cell bodies
# ---------------------------------------------------------------------------

TITLE_MD = """# L11 — Long Reads and the Pangenome

Long reads (PacBio HiFi, ONT) span structural variants outright and unlock two
techniques that short reads cannot match: (1) **direct SV detection from
coverage and split alignments**, and (2) **alignment to a pangenome graph**,
where the reference is no longer a single linear string but a DAG over
alternative haplotypes.

In this exercise you will:
1. Simulate a 100 kb reference, a homozygous 5 kb deletion, and 50 HiFi-style
   reads that span the breakpoints.
2. Recover the deletion from coverage and split-read evidence.
3. Build a 5-node pangenome GFA and visualise the DAG.
4. Generalise the linear-chain Viterbi from L02/L21 to a **Viterbi on a DAG**
   over a minimizer-seeded alignment problem.
5. Phase the reads into two haplotype sets by allele support across the
   deletion locus."""


AHA_MD = """> **Aha.** A **pangenome is a DAG** over the linear reference: variant branches
> attach to backbone nodes, and an alignment is now a *path through the graph*
> rather than a position on a string. Long reads can phase haplotypes
> **directly**, because a single read spans both the SV and the heterozygous
> SNVs nearby — the source-separation problem becomes trivial once your
> measurements are long enough.
>
> **Algorithmically:** Viterbi on a DAG is the *same* dynamic program as the
> linear-chain Viterbi from L02/L21 — we just process nodes in topological
> order and take the max over **all predecessors** instead of one. Same
> log-space arithmetic, same backtrace pointers."""


PREAMBLE_CODE = """# Colab preamble — pinned installs. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 matplotlib==3.8.4 networkx==3.2.1 -q"""


IMPORTS_CODE = """import math
import time
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx

# Deterministic for the whole notebook.
np.random.seed(42)
RNG = np.random.default_rng(42)

# Canonical DNA alphabet.
NT = "ACGT"
NT_SET = set(NT)
print("Numpy:", np.__version__, " | NetworkX:", nx.__version__)"""


# ---------------------------------------------------------------------------
# Step 1 — Generate ref + reads + spot SV in coverage
# ---------------------------------------------------------------------------

STEP1_MD = """## Step 1 (8 min) — Simulate a 100 kb reference with a 5 kb deletion; map reads; spot the SV in coverage

We build a synthetic 100 kb reference `REF` and a "sample" sequence
`SAMPLE` that has a clean 5 kb deletion at positions `[DEL_START, DEL_END)`
(default: 40 000 - 45 000). Then we draw 50 HiFi-style long reads of
~12 kb each from `SAMPLE` with a low (~1%) substitution error rate.

A naive aligner that maps each read to `REF` and records the matched-region
will leave a **coverage trough** at the deletion: reads that span the
breakpoints split into two segments, and the deleted interval simply has
nothing mapped to it.

Your job: implement the simulator, run a toy aligner, and plot the
coverage track. The deletion should jump out as a clean zero-coverage gap."""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Reference + sample with a 5 kb deletion; reads; coverage track.
# ----------------------------------------------------------------------

REF_LEN     = 100_000
DEL_START   = 40_000
DEL_END     = 45_000      # half-open; deletion size = DEL_END - DEL_START
N_READS     = 50
READ_LEN    = 12_000
ERR_RATE    = 0.01        # HiFi-like substitution rate


def simulate_reference(L: int, seed: int = 0) -> str:
    """Return a length-L DNA string drawn uniformly from ACGT."""
    # TODO: use np.random.default_rng(seed).choice(...)
    raise NotImplementedError


def apply_deletion(ref: str, start: int, end: int) -> str:
    """Return ref with positions [start, end) excised — the 'sample' haplotype."""
    # TODO
    raise NotImplementedError


def simulate_reads(sample: str, n: int, read_len: int, err_rate: float,
                   seed: int = 1) -> list[tuple[int, str]]:
    """Draw n reads from sample. Return list of (start_in_sample, read_seq).

    Each read should:
      - Have a uniformly random start position in [0, len(sample) - read_len].
      - Carry per-base substitution errors at rate err_rate.
    """
    # TODO
    raise NotImplementedError


def naive_align_to_ref(read_seq: str, read_start_in_sample: int,
                       del_start: int, del_end: int) -> list[tuple[int, int]]:
    """Map a read drawn from `sample` back to `ref` coordinates.

    Returns a list of (ref_start, ref_end) match intervals — one if the read
    sits entirely outside the deletion (i.e. fully on one side), two if it
    spans the breakpoint (a split alignment).

    The trick: positions in `sample` to the right of `del_start` correspond
    to positions in `ref` shifted by + (del_end - del_start).
    """
    # TODO
    raise NotImplementedError


# Build and visualise.
ref = simulate_reference(REF_LEN, seed=0)
sample = apply_deletion(ref, DEL_START, DEL_END)
reads = simulate_reads(sample, N_READS, READ_LEN, ERR_RATE, seed=1)

# Coverage track over REF coordinates.
coverage = np.zeros(REF_LEN, dtype=np.int32)
for start_in_sample, rseq in reads:
    for r_start, r_end in naive_align_to_ref(rseq, start_in_sample, DEL_START, DEL_END):
        coverage[r_start:r_end] += 1
'''


STEP1_SOL_MD = """*Click ▶ to expand the reference solution.*"""


STEP1_SOL = '''# Reference solution — Step 1.

REF_LEN     = 100_000
DEL_START   = 40_000
DEL_END     = 45_000
N_READS     = 50
READ_LEN    = 12_000
ERR_RATE    = 0.01


def simulate_reference(L: int, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    arr = rng.choice(list(NT), size=L)
    return "".join(arr)


def apply_deletion(ref: str, start: int, end: int) -> str:
    return ref[:start] + ref[end:]


def simulate_reads(sample: str, n: int, read_len: int, err_rate: float,
                   seed: int = 1):
    rng = np.random.default_rng(seed)
    max_start = len(sample) - read_len
    reads = []
    for _ in range(n):
        s = int(rng.integers(0, max_start + 1))
        rseq = list(sample[s:s + read_len])
        # Position-independent substitution errors.
        n_err = rng.binomial(read_len, err_rate)
        if n_err > 0:
            positions = rng.choice(read_len, size=n_err, replace=False)
            for p in positions:
                cur = rseq[p]
                alts = [c for c in NT if c != cur]
                rseq[p] = alts[int(rng.integers(0, 3))]
        reads.append((s, "".join(rseq)))
    return reads


def naive_align_to_ref(read_seq: str, read_start_in_sample: int,
                       del_start: int, del_end: int):
    """Map a read from sample back to ref. Returns 1 or 2 (ref_start, ref_end) intervals."""
    L = len(read_seq)
    s_start = read_start_in_sample
    s_end = s_start + L
    del_len = del_end - del_start

    if s_end <= del_start:
        # Entirely left of the breakpoint — direct mapping.
        return [(s_start, s_end)]
    if s_start >= del_start:
        # Entirely right of the breakpoint in `sample` — shift by deletion size.
        return [(s_start + del_len, s_end + del_len)]
    # Spans the breakpoint: split alignment.
    left_len = del_start - s_start
    left  = (s_start, s_start + left_len)
    right = (del_end, del_end + (L - left_len))
    return [left, right]


# Build and visualise.
ref = simulate_reference(REF_LEN, seed=0)
sample = apply_deletion(ref, DEL_START, DEL_END)
reads = simulate_reads(sample, N_READS, READ_LEN, ERR_RATE, seed=1)

coverage = np.zeros(REF_LEN, dtype=np.int32)
split_breakpoints = 0
for start_in_sample, rseq in reads:
    intervals = naive_align_to_ref(rseq, start_in_sample, DEL_START, DEL_END)
    if len(intervals) == 2:
        split_breakpoints += 1
    for r_start, r_end in intervals:
        coverage[r_start:r_end] += 1

print(f"Reference: {len(ref):,} bp")
print(f"Sample:    {len(sample):,} bp (deletion of {DEL_END-DEL_START} bp)")
print(f"Reads:     {len(reads)} x {READ_LEN} bp (error rate {ERR_RATE})")
print(f"Split (breakpoint-spanning) alignments: {split_breakpoints}")

fig, ax = plt.subplots(figsize=(11, 3.5))
ax.fill_between(np.arange(REF_LEN), coverage, color="#3aa07b", alpha=0.7)
ax.axvspan(DEL_START, DEL_END, color="#d36b3b", alpha=0.25,
           label=f"deletion [{DEL_START:,}, {DEL_END:,})")
ax.set_xlabel("Reference position (bp)")
ax.set_ylabel("Coverage")
ax.set_title("Coverage track — the deletion is a clean zero-coverage gap")
ax.legend(loc="upper right")
plt.tight_layout()
plt.show()
'''


# ---------------------------------------------------------------------------
# Step 2 — Pangenome GFA
# ---------------------------------------------------------------------------

STEP2_MD = """## Step 2 (10 min) — Build a pangenome GFA and visualise the graph

A pangenome encodes alternative haplotypes as **branches off a linear
backbone**. The simplest non-trivial case for our 100 kb region: model the
deletion locus as two paths through a 4-node bubble (`N1 → {N2, N3} → N4`).

To make the bubble discriminative for graph alignment, each branch needs to
carry **real, allele-specific sequence**. We give the bubble a `BRIDGE`-bp
shoulder on each side so junction-spanning reads have something to anchor on:

| Node | Span on REF                                                | What it represents |
|------|------------------------------------------------------------|--------------------|
| N1   | `ref[0 : del_start - BRIDGE]`                              | left flank backbone |
| N2   | `ref[del_start - BRIDGE : del_end + BRIDGE]`               | **reference allele in context** (5 kb deletion *present* in REF + BRIDGE on each side) |
| N3   | `ref[del_start - BRIDGE : del_start] + ref[del_end : del_end + BRIDGE]` | **deletion-junction allele** — the 2·BRIDGE-bp bridge across the excised interval |
| N4   | `ref[del_end + BRIDGE :]`                                  | right flank backbone |

The trick is N3: it is **not the empty allele**. It is the
breakpoint-junction sequence — the only place in the graph where the last
BRIDGE bp before `del_start` are immediately followed by the first BRIDGE bp
after `del_end`. A read drawn from the deletion haplotype that spans the
breakpoint contains this exact sequence; its minimizers therefore land
**only** on N3.

In **GFA-1 syntax**:
- `S <id> <sequence>` declares a segment.
- `L <from> <from_orient> <to> <to_orient> <overlap>` declares a link.
- The reference path is `N1 -> N2 -> N4`; the deletion path is `N1 -> N3 -> N4`.

Build the GFA, parse it into a directed graph with NetworkX, and visualise."""


STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Build a pangenome GFA; parse to networkx DiGraph.
# ----------------------------------------------------------------------

BRIDGE = 200  # bp of shoulder context flanking the bubble on each side

# Segments (see table in the markdown above):
#   N1 = ref[:del_start - BRIDGE]
#   N2 = ref[del_start - BRIDGE : del_end + BRIDGE]             (REF allele)
#   N3 = ref[del_start - BRIDGE : del_start]
#      + ref[del_end : del_end + BRIDGE]                        (DEL junction)
#   N4 = ref[del_end + BRIDGE :]


def build_pangenome_gfa(ref_seq: str,
                        del_start: int, del_end: int,
                        bridge: int = BRIDGE) -> str:
    """Return a GFA-1 string with 4 segments and 4 links.

    Layout:
      N1 = ref_seq[:del_start - bridge]
      N2 = ref_seq[del_start - bridge : del_end + bridge]           (REF allele)
      N3 = ref_seq[del_start - bridge : del_start]
           + ref_seq[del_end : del_end + bridge]                    (DEL junction)
      N4 = ref_seq[del_end + bridge :]

    Links: N1->N2, N1->N3, N2->N4, N3->N4.
    Paths: REF = N1,N2,N4; DEL = N1,N3,N4.
    """
    # TODO: assemble the GFA string with S/L/P lines and return it.
    raise NotImplementedError


def parse_gfa(gfa_text: str) -> tuple[dict, "nx.DiGraph", dict]:
    """Parse a GFA-1 string. Return (nodes, graph, paths) where:
        nodes  = {segment_id: sequence_str}
        graph  = nx.DiGraph with directed edges from L lines
        paths  = {path_name: [segment_id, ...]}
    """
    # TODO: split on '\\n', dispatch on the first token (S/L/P), build structures.
    raise NotImplementedError


def draw_pangenome(graph: "nx.DiGraph", nodes: dict, paths: dict):
    """Render the DAG with one row per path; colour ref vs del paths."""
    # TODO: networkx draw_networkx with a positions dict; legend by path.
    raise NotImplementedError


# gfa_text = build_pangenome_gfa(ref, DEL_START, DEL_END)
# nodes, graph, paths = parse_gfa(gfa_text)
# draw_pangenome(graph, nodes, paths)
'''


STEP2_SOL = '''# Reference solution — Step 2.

BRIDGE = 200


def build_pangenome_gfa(ref_seq: str, del_start: int, del_end: int,
                        bridge: int = BRIDGE) -> str:
    n1 = ref_seq[:del_start - bridge]
    n2 = ref_seq[del_start - bridge : del_end + bridge]
    n3 = (ref_seq[del_start - bridge : del_start]
          + ref_seq[del_end : del_end + bridge])
    n4 = ref_seq[del_end + bridge :]

    lines = [
        "H\\tVN:Z:1.0",
        f"S\\tN1\\t{n1}",
        f"S\\tN2\\t{n2}",
        f"S\\tN3\\t{n3}",
        f"S\\tN4\\t{n4}",
        "L\\tN1\\t+\\tN2\\t+\\t0M",
        "L\\tN1\\t+\\tN3\\t+\\t0M",
        "L\\tN2\\t+\\tN4\\t+\\t0M",
        "L\\tN3\\t+\\tN4\\t+\\t0M",
        "P\\tREF\\tN1+,N2+,N4+\\t*",
        "P\\tDEL\\tN1+,N3+,N4+\\t*",
    ]
    return "\\n".join(lines) + "\\n"


def parse_gfa(gfa_text: str):
    nodes = {}
    graph = nx.DiGraph()
    paths = {}
    for raw in gfa_text.splitlines():
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split("\\t")
        tag = parts[0]
        if tag == "S":
            seg_id, seq = parts[1], parts[2]
            nodes[seg_id] = seq
            graph.add_node(seg_id, length=len(seq))
        elif tag == "L":
            frm, _frm_orient, to, _to_orient, _overlap = parts[1:6]
            graph.add_edge(frm, to)
        elif tag == "P":
            name, segs, _ovl = parts[1], parts[2], parts[3]
            path = [s.rstrip("+-") for s in segs.split(",")]
            paths[name] = path
    return nodes, graph, paths


def draw_pangenome(graph: nx.DiGraph, nodes: dict, paths: dict):
    # Lay out nodes left-to-right by topological order, two rows for the bubble.
    pos = {
        "N1": (0, 0),
        "N2": (1, 0.5),   # reference allele — upper branch
        "N3": (1, -0.5),  # deletion allele — lower branch
        "N4": (2, 0),
    }
    fig, ax = plt.subplots(figsize=(9, 3.5))

    # Edges first, coloured by which path uses them.
    ref_edges = list(zip(paths["REF"][:-1], paths["REF"][1:]))
    del_edges = list(zip(paths["DEL"][:-1], paths["DEL"][1:]))
    nx.draw_networkx_edges(graph, pos, edgelist=ref_edges, ax=ax,
                           edge_color="#3aa07b", width=3.0, arrows=True,
                           arrowsize=20, connectionstyle="arc3,rad=0.0")
    nx.draw_networkx_edges(graph, pos, edgelist=del_edges, ax=ax,
                           edge_color="#d36b3b", width=3.0, arrows=True,
                           arrowsize=20, connectionstyle="arc3,rad=0.0")

    # Nodes.
    nx.draw_networkx_nodes(graph, pos, node_size=2200,
                           node_color="#e8eaf0", edgecolors="#333",
                           linewidths=1.5, ax=ax)
    labels = {nid: f"{nid}\\n{len(nodes[nid])} bp" for nid in graph.nodes}
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=10, ax=ax)

    ax.set_title("Pangenome DAG — reference allele (green) vs deletion allele (orange)")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(-1.2, 1.2)

    legend_handles = [
        mpatches.Patch(color="#3aa07b", label="REF path: N1 → N2 → N4"),
        mpatches.Patch(color="#d36b3b", label="DEL path: N1 → N3 → N4"),
    ]
    ax.legend(handles=legend_handles, loc="lower center",
              bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False)
    plt.tight_layout()
    plt.show()


gfa_text = build_pangenome_gfa(ref, DEL_START, DEL_END)
nodes, graph, paths = parse_gfa(gfa_text)
print(f"Parsed GFA: {len(nodes)} segments, {graph.number_of_edges()} links, {len(paths)} paths")
for name, p in paths.items():
    total_len = sum(len(nodes[n]) for n in p)
    print(f"  Path {name}: {' -> '.join(p)}  (total {total_len:,} bp)")

draw_pangenome(graph, nodes, paths)
'''


# ---------------------------------------------------------------------------
# Step 3 — Minimizer seeding across nodes
# ---------------------------------------------------------------------------

STEP3_MD = """## Step 3 (12 min) — Minimizer seeding across the graph

A **minimizer** of a window is the lexicographically smallest k-mer in that
window. Minimizer sketching reduces the seed count to roughly `1/(w-k+1)`
of all k-mers, while *guaranteeing* that overlapping windows share seeds.
This is the same locality-sensitive hashing trick that powers
minimap2 / Winnowmap.

For our graph aligner we collect minimizers **per node** of the GFA, then
combine them into a single `{minimizer: [(node, offset), ...]}` index. A
read's minimizers then resolve to candidate (node, offset) hits — possibly
spanning multiple nodes for breakpoint-crossing reads.

Implement minimizer extraction and build the per-graph index."""


STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Minimizer seeding across graph nodes.
# ----------------------------------------------------------------------

K_MIN = 15
W_MIN = 10  # window size; minimizer per (w-k+1)-mer sliding window


def minimizers(seq: str, k: int = K_MIN, w: int = W_MIN) -> list[tuple[str, int]]:
    """Return all (kmer, position) minimizers from `seq`.

    For each window of length (w + k - 1) starting at i in [0, len-w-k+2):
      - enumerate the w k-mers seq[i:i+k], seq[i+1:i+k+1], ..., seq[i+w-1:i+w+k-1]
      - pick the lexicographically smallest; record (kmer, absolute_position)
    Deduplicate consecutive duplicate hits (same position).
    """
    # TODO: sliding window + min selection.
    raise NotImplementedError


def build_graph_minimizer_index(nodes: dict, k: int = K_MIN, w: int = W_MIN) -> dict:
    """Return {kmer: [(node_id, offset), ...]} across all node sequences.

    Important: when two branch nodes (e.g. the REF allele N2 and the DEL
    junction N3) share flanking shoulder sequence with the backbone,
    minimizers in those shoulders appear in BOTH branches and carry no
    discriminative signal. **Drop minimizers that appear in more than one
    branch node** (N2, N3) so each branch's index contains only its
    branch-unique k-mers. Backbone nodes (N1, N4) keep all their seeds.
    """
    # TODO: minimizers() per node; then post-process so shared-branch kmers
    # are removed from the branch-node lists.
    raise NotImplementedError


def find_read_seeds(read_seq: str, graph_index: dict,
                    k: int = K_MIN, w: int = W_MIN) -> list[tuple[int, str, int]]:
    """For a read, compute its minimizers and look them up in graph_index.

    Return list of (read_pos, node_id, node_offset) hits.
    """
    # TODO
    raise NotImplementedError


# graph_index = build_graph_minimizer_index(nodes)
# seeds_for_read0 = find_read_seeds(reads[0][1], graph_index)
'''


STEP3_SOL = '''# Reference solution — Step 3.

K_MIN = 15
W_MIN = 10


def minimizers(seq: str, k: int = K_MIN, w: int = W_MIN):
    """Sliding-window minimizers."""
    L = len(seq)
    out = []
    win_kmers = w  # number of k-mers per window
    n_windows = L - k - w + 2
    if n_windows <= 0:
        return out
    last_pos = -1
    for i in range(n_windows):
        best_kmer = None
        best_pos = -1
        for j in range(win_kmers):
            kmer = seq[i + j : i + j + k]
            # Skip k-mers containing non-canonical bases.
            if not all(c in NT_SET for c in kmer):
                continue
            if best_kmer is None or kmer < best_kmer:
                best_kmer = kmer
                best_pos = i + j
        if best_kmer is not None and best_pos != last_pos:
            out.append((best_kmer, best_pos))
            last_pos = best_pos
    return out


BRANCH_NODES = {"N2", "N3"}


def build_graph_minimizer_index(nodes: dict, k: int = K_MIN, w: int = W_MIN) -> dict:
    """Per-node minimizer index, with branch-shared k-mers removed from N2/N3.

    Step 1 builds a raw `{kmer: [(node_id, pos), ...]}` map from per-node
    minimizers. Step 2 walks the map: for each k-mer, look at how many
    *branch* nodes (N2 / N3) it appears in. If more than one, that k-mer
    is non-discriminative for the bubble — drop its branch-node entries
    (but keep entries on backbone nodes like N1 / N4 if any).
    """
    raw = defaultdict(list)
    for node_id, seq in nodes.items():
        for kmer, pos in minimizers(seq, k, w):
            raw[kmer].append((node_id, pos))

    idx = {}
    for kmer, hits in raw.items():
        branch_hits = [h for h in hits if h[0] in BRANCH_NODES]
        if len(branch_hits) > 1:
            # Non-discriminative for the bubble — keep only non-branch hits.
            kept = [h for h in hits if h[0] not in BRANCH_NODES]
            if kept:
                idx[kmer] = kept
        else:
            idx[kmer] = hits
    return idx


def find_read_seeds(read_seq: str, graph_index: dict,
                    k: int = K_MIN, w: int = W_MIN):
    hits = []
    for kmer, read_pos in minimizers(read_seq, k, w):
        for node_id, node_off in graph_index.get(kmer, []):
            hits.append((read_pos, node_id, node_off))
    return hits


t0 = time.time()
graph_index = build_graph_minimizer_index(nodes)
print(f"Graph index: {len(graph_index):,} unique minimizers "
      f"({sum(len(v) for v in graph_index.values()):,} total positions) "
      f"in {time.time()-t0:.2f}s")

# Diagnostic: where do read 0's seeds land across the graph?
seeds_for_read0 = find_read_seeds(reads[0][1], graph_index)
node_hits = defaultdict(int)
for _, node_id, _ in seeds_for_read0:
    node_hits[node_id] += 1
print(f"\\nRead 0 (drawn from sample pos {reads[0][0]:,}): "
      f"{len(seeds_for_read0)} total seeds")
for node_id in sorted(node_hits):
    print(f"  on node {node_id}: {node_hits[node_id]} seeds")

# Aggregate across all 50 reads.
all_seeds = [find_read_seeds(rseq, graph_index) for _, rseq in reads]
print(f"\\nAcross 50 reads: median {int(np.median([len(s) for s in all_seeds]))} "
      f"seeds per read")
'''


# ---------------------------------------------------------------------------
# Step 4 — Viterbi on a DAG
# ---------------------------------------------------------------------------

STEP4_MD = """## Step 4 (15 min) — Viterbi on a DAG: pick the best path per read

This is the heart of the exercise. Given a read and its seeds across the
graph, we want the **single best path** the read takes through the DAG.
We score each candidate (node, alignment) by the number of supporting
minimizer seeds on that node, then propagate scores along graph edges.

**The recurrence** (in log-space, $-\\infty$ for "no support"):

$$
\\delta[v] = \\log s(v) \\;+\\; \\max_{u \\in \\text{pred}(v)} \\,\\delta[u]
$$

where $s(v)$ is the seed support for node $v$ from this read, and the max
is taken over all predecessors $u$ of $v$ in the DAG. This is *exactly*
the linear-chain Viterbi from L02/L21 — the only change is that nodes are
processed in **topological order** and the max ranges over multiple
predecessors instead of one. Backtrace pointers store the argmax
predecessor.

Decoding a read into REF vs DEL is then trivial: which path through the
bubble (N2 or N3) does the best alignment go through?"""


STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Viterbi on a DAG over graph nodes.
# ----------------------------------------------------------------------

LOG_ZERO = -1e18


def node_support(seeds: list[tuple[int, str, int]]) -> dict:
    """Return {node_id: count} — number of seed hits to each node from this read."""
    # TODO
    raise NotImplementedError


def viterbi_on_dag(graph: nx.DiGraph, node_support_counts: dict,
                   topo_order: list[str]) -> tuple[list[str], float]:
    """Find the maximum-score path through the DAG.

    Score of a node v = log(1 + support_count[v])  (additive, monotone in support)
    Path score        = sum over nodes on path of node score
    Recurrence:
        delta[v] = score(v) + max over u in pred(v) of delta[u]    (or 0 if v is a source)
        ptr[v]   = argmax u   (or None if v is a source)

    Return (best_path_node_ids, best_path_score).
    """
    # TODO
    raise NotImplementedError


def classify_read(path_nodes: list[str], support: dict) -> str:
    """Return 'REF' if path traverses N2, 'DEL' if it traverses N3, else 'AMBIG'.

    Important: if BOTH N2 and N3 have zero seed support, the read carries
    no information about the bubble (it sits entirely in flanking N1 / N4)
    and the call is AMBIG. Don't be fooled by an arbitrary Viterbi
    tie-break in that case.
    """
    # TODO
    raise NotImplementedError


# topo = list(nx.topological_sort(graph))
# decoded = []
# for (_, rseq), s in zip(reads, all_seeds):
#     support = node_support(s)
#     path, score = viterbi_on_dag(graph, support, topo)
#     decoded.append((classify_read(path, support), path, score))
'''


STEP4_SOL = '''# Reference solution — Step 4.

LOG_ZERO = -1e18


def node_support(seeds):
    counts = defaultdict(int)
    for _read_pos, node_id, _node_off in seeds:
        counts[node_id] += 1
    return dict(counts)


def viterbi_on_dag(graph: nx.DiGraph, node_support_counts: dict,
                   topo_order):
    """Max-score path through the DAG by Viterbi-on-DAG.

    Node score: log(1 + support_count[v]).  Sources start at their own score.
    """
    delta = {}
    ptr = {}
    for v in topo_order:
        score_v = math.log1p(node_support_counts.get(v, 0))
        preds = list(graph.predecessors(v))
        if not preds:
            delta[v] = score_v
            ptr[v] = None
            continue
        best_u = None
        best_pred_score = LOG_ZERO
        for u in preds:
            if delta[u] > best_pred_score:
                best_pred_score = delta[u]
                best_u = u
        delta[v] = score_v + best_pred_score
        ptr[v] = best_u

    # Best terminal node = argmax over all sinks (nodes with no successors)
    # AND, for safety, over all nodes (the optimal path always ends at a sink
    # for a connected DAG with non-negative scores).
    sinks = [v for v in topo_order if graph.out_degree(v) == 0]
    if not sinks:
        sinks = list(topo_order)
    end = max(sinks, key=lambda v: delta[v])
    score = float(delta[end])

    # Backtrace.
    path = [end]
    while ptr[path[-1]] is not None:
        path.append(ptr[path[-1]])
    path.reverse()
    return path, score


def classify_read(path_nodes, support):
    # If neither branch has any seed support, the read is not informative
    # for the SV — it lies entirely in the flanking backbone.
    if support.get("N2", 0) == 0 and support.get("N3", 0) == 0:
        return "AMBIG"
    if "N2" in path_nodes and "N3" not in path_nodes:
        return "REF"
    if "N3" in path_nodes and "N2" not in path_nodes:
        return "DEL"
    return "AMBIG"


topo = list(nx.topological_sort(graph))
print(f"Topological order: {topo}")

decoded = []
for (read_start, rseq), s in zip(reads, all_seeds):
    support = node_support(s)
    path, score = viterbi_on_dag(graph, support, topo)
    call = classify_read(path, support)
    decoded.append({"start": read_start, "seq": rseq, "path": path,
                    "score": score, "call": call, "support": support})

# Tally.
call_counts = defaultdict(int)
for d in decoded:
    call_counts[d["call"]] += 1
print(f"\\nCalls across 50 reads: REF={call_counts['REF']}  "
      f"DEL={call_counts['DEL']}  AMBIG={call_counts['AMBIG']}")
print(f"  (AMBIG = read covers neither N2 nor N3 — sits entirely in flanking backbone.)")

# Show 5 example reads.
print("\\nExample reads (first 5):")
for i, d in enumerate(decoded[:5]):
    p = " -> ".join(d["path"])
    print(f"  read {i:2d}  start_in_sample={d['start']:>6,}  "
          f"call={d['call']:<5s}  path={p}  score={d['score']:.2f}")

# Aggregate SV genotype call: majority over informative (non-AMBIG) reads.
informative = [d for d in decoded if d["call"] != "AMBIG"]
inf_ref = sum(1 for d in informative if d["call"] == "REF")
inf_del = sum(1 for d in informative if d["call"] == "DEL")
print(f"\\nInformative reads (those crossing the bubble): {len(informative)}/50")
print(f"  REF support: {inf_ref}  |  DEL support: {inf_del}")
sample_genotype = "DEL" if inf_del > inf_ref else ("REF" if inf_ref > inf_del else "AMBIG")
print(f"  Aggregated SV call: {sample_genotype}  (truth = DEL)")
'''


# ---------------------------------------------------------------------------
# Step 5 — Phasing
# ---------------------------------------------------------------------------

STEP5_MD = """## Step 5 (15 min) — Phase the reads by allele support

In a real sample with a heterozygous deletion, half the reads carry the
reference allele and half carry the deletion. Phasing then becomes
**source separation**: partition the reads into two haplotype sets based
on which allele they support.

Two important caveats up front:

- **Only reads that span the bubble are phasable.** A read entirely inside
  N1 or N4 carries no allele information — we call those AMBIG. With 50
  random 12 kb reads on a 95-100 kb sample, only ~10-15 reads will cross
  the breakpoint. That handful is enough.
- **Our Step 1 sample was homozygous DEL**, so direct phasing on it gives
  one set; the source-separation aha is more interesting on a
  **heterozygous mock** where half the reads come from REF, half from DEL.

Implement the het-mixed simulator, re-run the graph aligner on it, and
verify that informative reads partition cleanly by their Viterbi path."""


STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Phase the reads by allele support (heterozygous mock).
# ----------------------------------------------------------------------

def simulate_het_reads(ref_seq: str, del_start: int, del_end: int,
                       n_reads: int, read_len: int, err_rate: float,
                       seed: int = 7):
    """Half the reads from REF haplotype (ref_seq), half from DEL haplotype.

    Returns list of dicts with keys:
        truth   : 'REF' or 'DEL' (ground-truth source)
        seq     : the read sequence
    """
    # TODO: simulate n_reads/2 reads from ref_seq, n_reads/2 from the deletion
    # haplotype = apply_deletion(ref_seq, del_start, del_end). Add per-base
    # errors.
    raise NotImplementedError


def phase_reads(het_reads, graph, graph_index, topo):
    """For each het read: run Viterbi on the DAG; classify as REF/DEL/AMBIG.

    Return a list of dicts adding 'call' = REF/DEL/AMBIG.
    """
    # TODO
    raise NotImplementedError


def phasing_accuracy(phased_reads):
    """Fraction of reads whose call matches truth (AMBIG counts as wrong)."""
    # TODO
    raise NotImplementedError


# het_reads = simulate_het_reads(ref, DEL_START, DEL_END, n_reads=50,
#                                 read_len=READ_LEN, err_rate=ERR_RATE)
# phased = phase_reads(het_reads, graph, graph_index, topo)
# acc = phasing_accuracy(phased)
'''


STEP5_SOL = '''# Reference solution — Step 5.

def simulate_het_reads(ref_seq, del_start, del_end, n_reads,
                       read_len, err_rate, seed=7):
    rng = np.random.default_rng(seed)
    n_ref = n_reads // 2
    n_del = n_reads - n_ref
    out = []

    # REF haplotype = ref_seq.
    max_start_ref = len(ref_seq) - read_len
    for _ in range(n_ref):
        s = int(rng.integers(0, max_start_ref + 1))
        rseq = list(ref_seq[s:s + read_len])
        n_err = rng.binomial(read_len, err_rate)
        if n_err > 0:
            positions = rng.choice(read_len, size=n_err, replace=False)
            for p in positions:
                cur = rseq[p]
                alts = [c for c in NT if c != cur]
                rseq[p] = alts[int(rng.integers(0, 3))]
        out.append({"truth": "REF", "seq": "".join(rseq)})

    # DEL haplotype.
    del_hap = apply_deletion(ref_seq, del_start, del_end)
    max_start_del = len(del_hap) - read_len
    for _ in range(n_del):
        s = int(rng.integers(0, max_start_del + 1))
        rseq = list(del_hap[s:s + read_len])
        n_err = rng.binomial(read_len, err_rate)
        if n_err > 0:
            positions = rng.choice(read_len, size=n_err, replace=False)
            for p in positions:
                cur = rseq[p]
                alts = [c for c in NT if c != cur]
                rseq[p] = alts[int(rng.integers(0, 3))]
        out.append({"truth": "DEL", "seq": "".join(rseq)})

    rng.shuffle(out)
    return out


def phase_reads(het_reads, graph, graph_index, topo):
    out = []
    for r in het_reads:
        seeds = find_read_seeds(r["seq"], graph_index)
        support = node_support(seeds)
        path, score = viterbi_on_dag(graph, support, topo)
        call = classify_read(path, support)
        out.append({**r, "call": call, "path": path, "score": score,
                    "support": support})
    return out


def phasing_accuracy(phased, informative_only: bool = True):
    """Phasing accuracy.

    informative_only=True restricts to reads whose call is REF or DEL
    (i.e., reads that actually crossed the bubble — non-AMBIG). Reads
    that sit entirely in flanking backbone cannot be phased; their
    "AMBIG" call is a truthful refusal, not a misclassification.
    """
    pool = [r for r in phased if r["call"] != "AMBIG"] if informative_only else phased
    if not pool:
        return 0.0
    correct = sum(1 for r in pool if r["call"] == r["truth"])
    return correct / len(pool)


het_reads = simulate_het_reads(ref, DEL_START, DEL_END, n_reads=50,
                                read_len=READ_LEN, err_rate=ERR_RATE,
                                seed=7)
phased = phase_reads(het_reads, graph, graph_index, topo)

# Confusion matrix.
truth_calls = defaultdict(lambda: defaultdict(int))
for r in phased:
    truth_calls[r["truth"]][r["call"]] += 1

acc_inf = phasing_accuracy(phased, informative_only=True)
n_inf = sum(1 for r in phased if r["call"] != "AMBIG")
n_correct_inf = sum(1 for r in phased if r["call"] != "AMBIG" and r["call"] == r["truth"])
print(f"Phasing accuracy on informative (bubble-spanning) reads: "
      f"{acc_inf:.3f} ({n_correct_inf}/{n_inf})")
print(f"  ({len(phased) - n_inf}/{len(phased)} reads are AMBIG — entirely in flank, "
      "no allele information.)")
print()
print("Confusion matrix (rows = truth, cols = call):")
print(f"{'':>8s} {'REF':>6s} {'DEL':>6s} {'AMBIG':>6s}")
for truth in ("REF", "DEL"):
    row = truth_calls[truth]
    print(f"{truth:>8s} {row['REF']:>6d} {row['DEL']:>6d} {row['AMBIG']:>6d}")

# Source-separated read sets.
ref_set = [r for r in phased if r["call"] == "REF"]
del_set = [r for r in phased if r["call"] == "DEL"]
print(f"\\nSeparated into {len(ref_set)} REF reads and {len(del_set)} DEL reads "
      f"(ambig: {sum(1 for r in phased if r['call'] == 'AMBIG')}).")

# Visualisation: stack of reads coloured by truth, marked by call.
fig, ax = plt.subplots(figsize=(10, 4))
for i, r in enumerate(phased):
    color = "#3aa07b" if r["truth"] == "REF" else "#d36b3b"
    marker = "o" if r["call"] == r["truth"] else "x"
    ax.scatter([i], [r["score"]], c=color, marker=marker, s=50,
               edgecolors="black", linewidths=0.5)
ax.set_xlabel("Read index (shuffled)")
ax.set_ylabel("Viterbi log-score")
ax.set_title("Per-read Viterbi scores — colour = truth haplotype, x = misclassified")

legend_handles = [
    mpatches.Patch(color="#3aa07b", label="truth REF"),
    mpatches.Patch(color="#d36b3b", label="truth DEL"),
]
ax.legend(handles=legend_handles, loc="upper right")
plt.tight_layout()
plt.show()
'''


# ---------------------------------------------------------------------------
# Self-check
# ---------------------------------------------------------------------------

SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing pieces of the pipeline. If you ran
the reference solutions above they will pass; if you wrote your own and an
assert fails, revisit the corresponding step."""


SELFCHECK_CODE = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1) Reference / sample dimensions.
assert len(ref) == REF_LEN, f"Reference length wrong: {len(ref)}"
assert len(sample) == REF_LEN - (DEL_END - DEL_START), (
    f"Sample length wrong: {len(sample)} vs expected {REF_LEN - (DEL_END - DEL_START)}"
)
print(f"[1] Reference: {len(ref):,} bp; sample: {len(sample):,} bp  OK")

# 2) Coverage track has a zero-coverage gap exactly at the deletion.
gap = coverage[DEL_START:DEL_END]
flank_left = coverage[max(0, DEL_START - 5000):DEL_START]
flank_right = coverage[DEL_END:DEL_END + 5000]
assert gap.max() == 0, f"Coverage inside deletion should be 0, got max={gap.max()}"
assert flank_left.mean() > 1.0, (
    f"Left flank mean coverage too low: {flank_left.mean():.2f}"
)
assert flank_right.mean() > 1.0, (
    f"Right flank mean coverage too low: {flank_right.mean():.2f}"
)
print(f"[2] Coverage gap clean: 0 inside [{DEL_START:,}, {DEL_END:,}); "
      f"flank mean ~ {0.5 * (flank_left.mean() + flank_right.mean()):.1f}  OK")

# 3) Pangenome graph topology: 4 nodes, 4 directed edges, single source N1,
#    single sink N4.
assert set(graph.nodes) == {"N1", "N2", "N3", "N4"}, f"Bad node set: {set(graph.nodes)}"
assert graph.number_of_edges() == 4, f"Bad edge count: {graph.number_of_edges()}"
assert nx.is_directed_acyclic_graph(graph), "Graph must be a DAG"
sources = [v for v in graph.nodes if graph.in_degree(v) == 0]
sinks   = [v for v in graph.nodes if graph.out_degree(v) == 0]
assert sources == ["N1"], f"Bad sources: {sources}"
assert sinks == ["N4"], f"Bad sinks: {sinks}"
print(f"[3] DAG: 4 nodes / 4 edges / DAG / single source N1 / single sink N4  OK")

# 4) Viterbi on a tiny hand-constructed DAG: 3-node line A -> B -> C with
#    weights {A:1, B:1, C:1}: best path = A,B,C with score = 3 * log1p(1).
test_graph = nx.DiGraph()
test_graph.add_edges_from([("A", "B"), ("B", "C")])
test_support = {"A": 1, "B": 1, "C": 1}
test_topo = ["A", "B", "C"]
test_path, test_score = viterbi_on_dag(test_graph, test_support, test_topo)
expected = 3 * math.log1p(1)
assert test_path == ["A", "B", "C"], f"Hand-test path wrong: {test_path}"
assert abs(test_score - expected) < 1e-9, (
    f"Hand-test score wrong: {test_score} vs {expected}"
)
print(f"[4] Viterbi-on-DAG hand test: path A->B->C, score {test_score:.4f}  OK")

# 5) Viterbi on a 3-node BUBBLE A -> {B, C} -> D with support {A:0, B:5, C:0, D:0}
#    should pick the A,B,D branch.
bubble = nx.DiGraph()
bubble.add_edges_from([("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")])
bubble_support = {"A": 0, "B": 5, "C": 0, "D": 0}
bubble_topo = ["A", "B", "C", "D"]
bp, _bs = viterbi_on_dag(bubble, bubble_support, bubble_topo)
assert bp == ["A", "B", "D"], f"Bubble Viterbi went wrong: {bp}"
print(f"[5] Bubble Viterbi prefers high-support branch: {bp}  OK")

# 6) On the homozygous-deletion sample (Step 4 decoded all 50 reads):
#    among reads that actually cross the bubble (call != AMBIG), virtually
#    all should be called DEL. Reads outside the bubble can't be phased.
informative_calls = [d["call"] for d in decoded if d["call"] != "AMBIG"]
assert len(informative_calls) >= 5, (
    f"Only {len(informative_calls)} informative reads — sample too sparse"
)
del_frac = informative_calls.count("DEL") / len(informative_calls)
assert del_frac >= 0.85, (
    f"Of {len(informative_calls)} informative reads, only {del_frac:.0%} "
    "called DEL on homozygous-deletion sample"
)
print(f"[6] Homozygous-deletion: {informative_calls.count('DEL')}/{len(informative_calls)} "
      f"informative reads called DEL  OK")

# 7) Heterozygous mock phasing: accuracy on informative reads >= 0.85.
het_acc = phasing_accuracy(phased, informative_only=True)
n_inf_het = sum(1 for r in phased if r["call"] != "AMBIG")
assert n_inf_het >= 5, f"Only {n_inf_het} informative het reads — too few"
assert het_acc >= 0.85, (
    f"Het phasing accuracy {het_acc:.3f} < 0.85 on {n_inf_het} informative reads"
)
print(f"[7] Heterozygous phasing accuracy: {het_acc:.3f} on {n_inf_het} informative reads  OK")

print()
print("✅ Self-check passed.")
'''


# ---------------------------------------------------------------------------
# EE framing closing
# ---------------------------------------------------------------------------

EE_MD = """## EE framing — DAG Viterbi, LSH seeding, source separation

You implemented the three load-bearing ideas in graph genomics:

- **Pangenome = DAG over the linear reference.** Variant alleles are
  branches off a backbone. An alignment is now a *path* through the graph;
  positions are `(node_id, offset)` pairs, not a single integer. The
  algorithmic upgrade from L02 (FM-index on a string) is conceptual, not
  enormous.

- **Minimizer seeding = locality-sensitive hashing / sketching.** Picking
  the smallest k-mer in every window guarantees that overlapping windows
  share at least one minimizer (the sketch is *locality-preserving*),
  while reducing seed counts by ~ `1/(w-k+1)`. This is the same idea
  behind MinHash and the count-sketch sublinear estimators.

- **Viterbi on a DAG is the same DP as Viterbi on a chain.** Process
  nodes in topological order; max over **all predecessors** instead of
  one; everything stays in log-space; backtrace pointers store the
  argmax predecessor. The L02/L21 max-product belief-propagation
  intuition transfers verbatim. The only place the chain assumption was
  load-bearing was the *predecessor count*; the rest is identical.

- **Phasing = source separation, but only on informative measurements.**
  Half the reads come from one haplotype, half from the other; the
  graph alignment naturally partitions reads that **span the bubble**
  by which allele node their best path traversed. Reads sitting entirely
  in the flanking backbone carry no allele information and must be
  declared AMBIG — refusing to call is the right move when the
  measurement has zero discriminating signal. Long reads make this
  trivial: a single 12 kb read spans the SV and reveals its haplotype.
  With short reads you would need per-SNV statistics and an iterative
  haplotype-assembly step (HapCUT, WhatsHap)."""


# ---------------------------------------------------------------------------
# Assemble notebook
# ---------------------------------------------------------------------------

def build():
    nb = new_notebook()
    cells = []

    cells.append(md(TITLE_MD))
    cells.append(md(AHA_MD))
    cells.append(code(PREAMBLE_CODE))
    cells.append(code(IMPORTS_CODE))

    # Step 1
    cells.append(md(STEP1_MD))
    cells.append(code(STEP1_TODO))
    cells.append(md(STEP1_SOL_MD))
    cells.append(hidden_solution(STEP1_SOL))

    # Step 2
    cells.append(md(STEP2_MD))
    cells.append(code(STEP2_TODO))
    cells.append(md("*Click ▶ to expand the reference solution.*"))
    cells.append(hidden_solution(STEP2_SOL))

    # Step 3
    cells.append(md(STEP3_MD))
    cells.append(code(STEP3_TODO))
    cells.append(md("*Click ▶ to expand the reference solution.*"))
    cells.append(hidden_solution(STEP3_SOL))

    # Step 4
    cells.append(md(STEP4_MD))
    cells.append(code(STEP4_TODO))
    cells.append(md("*Click ▶ to expand the reference solution.*"))
    cells.append(hidden_solution(STEP4_SOL))

    # Step 5
    cells.append(md(STEP5_MD))
    cells.append(code(STEP5_TODO))
    cells.append(md("*Click ▶ to expand the reference solution.*"))
    cells.append(hidden_solution(STEP5_SOL))

    # Self-check + closing
    cells.append(md(SELFCHECK_MD))
    cells.append(code(SELFCHECK_CODE))
    cells.append(md(EE_MD))

    nb.cells = cells
    nb.metadata = {
        "kernelspec": {
            "display_name": "Python 3",
            "language":     "python",
            "name":         "python3",
        },
        "language_info": {"name": "python", "version": "3.11"},
        "colab": {"provenance": [], "toc_visible": True},
    }

    nbf.write(nb, OUT)
    return nb


if __name__ == "__main__":
    nb = build()
    nb2 = nbf.read(OUT, as_version=4)
    n_md = sum(1 for c in nb2.cells if c.cell_type == "markdown")
    n_code = sum(1 for c in nb2.cells if c.cell_type == "code")
    n_hidden = sum(
        1 for c in nb2.cells
        if c.cell_type == "code"
        and c.metadata.get("jupyter", {}).get("source_hidden")
    )
    print(f"Wrote {OUT}")
    print(f"  cells: {len(nb2.cells)} total  ({n_md} md, {n_code} code, {n_hidden} hidden-solution)")
