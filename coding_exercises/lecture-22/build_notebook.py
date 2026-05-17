"""Build exercise.ipynb for L22 — Network Biology and Pathway Analysis.

Run from this directory:
    python3 build_notebook.py

Emits exercise.ipynb. Re-running overwrites the file.
"""

from __future__ import annotations

import os
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


TITLE_MD = """# L22 — Network Biology and Pathway Analysis

In this exercise you treat a protein-protein interaction (PPI) network as a
**graph signal-processing** problem. You will:

1. Build a 500-node synthetic PPI with planted communities.
2. Run **random walk with restart (RWR)** from a 20-gene seed set and sweep
   the restart probability.
3. Find communities with the **Louvain** algorithm and compare modularity
   against the ground-truth partition.
4. Run **hypergeometric (ORA)** pathway enrichment on RWR's top-30 list
   against a tiny toy gene-set database.
5. Re-interpret RWR as a **low-pass filter on the graph Laplacian**: the
   restart probability is the cutoff frequency.

Everything is synthetic and deterministic; no network fetches.
"""


AHA_MD = """> **Aha.** Network propagation is a **low-pass filter** in the graph Fourier
> basis: the graph Laplacian's small-eigenvalue eigenvectors are slowly
> varying across edges, and RWR damps the high-frequency components by
> $\\lambda_i / (\\lambda_i + \\alpha)$ where $\\alpha = r / (1-r)$. Louvain
> partitions the same Laplacian's low-frequency subspace — it is essentially
> k-means in eigenspace. Hypergeometric enrichment is a two-sample binomial
> on a finite urn. Three apparently-distinct ideas, one spectral story.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 networkx==3.3 matplotlib==3.8.4 -q
"""


IMPORTS = """import math
import random
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import networkx as nx
import matplotlib.pyplot as plt

# Deterministic for the whole notebook.
SEED = 42
np.random.seed(SEED)
random.seed(SEED)

print(f"networkx {nx.__version__}, scipy {sp.__name__.split('.')[0]} loaded")
"""


STEP1_MD = """## Step 1 (8 min) — Build a 500-node synthetic PPI with planted communities

A real PPI from STRING has thousands of proteins and a thick tail of
high-degree hubs. We replicate the structure — but not the size — with a
**stochastic block model**: 500 nodes split into 5 communities of 100 nodes
each, with a high within-block edge probability and a low between-block
probability. The result is a sparse graph with a clear (but noisy) modular
structure, which is exactly what we want for the rest of the exercise.

We label communities `C0`–`C4` and treat them as **ground-truth modules**
for the modularity comparison in Step 3. Nodes are named `g000`–`g499`;
genes `g000`–`g019` (i.e. the first 20 nodes, all in community C0) form
the **seed set** for the RWR.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Build the synthetic PPI graph.
# ----------------------------------------------------------------------

N_NODES        = 500
N_COMMUNITIES  = 5
NODES_PER_COMM = N_NODES // N_COMMUNITIES   # 100
P_WITHIN       = 0.08    # edge probability within a community
P_BETWEEN      = 0.002   # edge probability between communities
SEED_SIZE      = 20      # first 20 nodes form the RWR seed (all in C0)


def build_ppi(seed: int = SEED):
    """Return (G, ground_truth) where G is a 500-node nx.Graph and
    ground_truth is a dict {node_name: community_label}."""
    # TODO: use nx.stochastic_block_model with sizes=[100]*5 and a 5x5
    # probability matrix with P_WITHIN on the diagonal and P_BETWEEN
    # off-diagonal. Rename nodes to "g000"..."g499". Drop any isolates so
    # the graph is one big connected component (use the largest component).
    raise NotImplementedError


# Build the graph and inspect it.
# G, ground_truth = build_ppi()
# print(f"nodes={G.number_of_nodes()}, edges={G.number_of_edges()}")
'''


STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

N_NODES        = 500
N_COMMUNITIES  = 5
NODES_PER_COMM = N_NODES // N_COMMUNITIES
P_WITHIN       = 0.08
P_BETWEEN      = 0.002
SEED_SIZE      = 20


def build_ppi(seed: int = SEED):
    sizes = [NODES_PER_COMM] * N_COMMUNITIES
    P = np.full((N_COMMUNITIES, N_COMMUNITIES), P_BETWEEN)
    np.fill_diagonal(P, P_WITHIN)

    G_raw = nx.stochastic_block_model(sizes, P.tolist(), seed=seed)

    # Rename 0..499 -> g000..g499 and record community labels.
    mapping = {i: f"g{i:03d}" for i in G_raw.nodes()}
    G = nx.relabel_nodes(G_raw, mapping)
    ground_truth = {}
    for i in G_raw.nodes():
        # nx assigns the block label as a node attribute "block".
        ground_truth[mapping[i]] = f"C{G_raw.nodes[i]['block']}"

    # Keep the largest connected component so RWR has no surprises with isolates.
    largest_cc = max(nx.connected_components(G), key=len)
    G = G.subgraph(largest_cc).copy()
    ground_truth = {n: ground_truth[n] for n in G.nodes()}
    return G, ground_truth


G, ground_truth = build_ppi()
print(f"nodes      = {G.number_of_nodes()}")
print(f"edges      = {G.number_of_edges()}")
print(f"density    = {nx.density(G):.4f}")
print(f"avg degree = {2 * G.number_of_edges() / G.number_of_nodes():.2f}")

# Validate: degree distribution + connectivity.
degrees = [d for _, d in G.degree()]
print(f"degree min/median/max = {min(degrees)} / {int(np.median(degrees))} / {max(degrees)}")
print(f"connected? {nx.is_connected(G)}")

# Communities are ~equal-sized.
print("community sizes:", Counter(ground_truth.values()))

# Plot the degree distribution.
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.hist(degrees, bins=30, color="steelblue", alpha=0.8)
ax.set_xlabel("degree")
ax.set_ylabel("count")
ax.set_title("Degree distribution of the synthetic PPI")
plt.tight_layout()
plt.show()

# The first 20 nodes (g000..g019) all live in community C0 — they are our seed.
SEED_NODES = [f"g{i:03d}" for i in range(SEED_SIZE)]
assert all(ground_truth[n] == "C0" for n in SEED_NODES), \\
    "Seed nodes should all be in community C0"
print(f"\\nSeed set: {len(SEED_NODES)} nodes, all in C0.")
'''


STEP2_MD = """## Step 2 (12 min) — RWR from a 20-gene seed, sweep restart probability

**Random walk with restart** simulates a random walker that, at every step,
either follows a random edge (with probability $1 - r$) or jumps back to one
of the seed nodes (with probability $r$). The stationary distribution
$\\mathbf{p}$ ranks every node by its "topological closeness" to the seeds.

The closed-form linear-system view is:
$$\\mathbf{p} = r \\, (I - (1-r) \\, W)^{-1} \\, \\mathbf{p}_0$$
where $W$ is the column-stochastic transition matrix (each column sums to 1)
and $\\mathbf{p}_0$ is the uniform seed distribution. We solve this directly
with a sparse linear solve.

You will sweep $r \\in \\{0.3, 0.5, 0.7\\}$ and observe how the **diffusion
radius** shrinks as $r$ grows — small $r$ means the walker wanders far
before restarting, so the propagation spreads across more of the network.
"""


STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Random walk with restart.
# ----------------------------------------------------------------------

RESTART_PROBS = [0.3, 0.5, 0.7]


def column_stochastic_W(G):
    """Return (W, nodes) where W is the column-stochastic transition matrix
    of G as a scipy.sparse matrix and nodes is the ordered node list."""
    # TODO: A = nx.adjacency_matrix; D_inv = 1/degree; W = A @ D_inv
    # (each column of W sums to 1).
    raise NotImplementedError


def rwr(G, seeds, r: float):
    """Return a pandas.Series {node: RWR score} solving p = r*(I - (1-r)*W)^{-1} * p0."""
    # TODO:
    # 1. Build W and the seed distribution p0 (uniform over `seeds`, zero elsewhere).
    # 2. Solve (I - (1-r) * W) x = r * p0 with scipy.sparse.linalg.spsolve.
    # 3. Return the result as a Series indexed by node names.
    raise NotImplementedError


# Sweep RESTART_PROBS and tabulate the top-30 nodes.
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.
from scipy.sparse.linalg import spsolve

RESTART_PROBS = [0.3, 0.5, 0.7]


def column_stochastic_W(G):
    nodes = list(G.nodes())
    A = nx.to_scipy_sparse_array(G, nodelist=nodes, format="csr", dtype=float)
    deg = np.asarray(A.sum(axis=0)).flatten()
    deg_inv = np.where(deg > 0, 1.0 / deg, 0.0)
    D_inv = sp.diags(deg_inv)
    W = A @ D_inv   # column-stochastic
    return W.tocsr(), nodes


def rwr(G, seeds, r: float):
    W, nodes = column_stochastic_W(G)
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)

    p0 = np.zeros(n)
    for s in seeds:
        p0[idx[s]] = 1.0
    p0 /= p0.sum()

    I = sp.eye(n, format="csr")
    M = (I - (1.0 - r) * W).tocsc()
    p = spsolve(M, r * p0)
    return pd.Series(p, index=nodes).sort_values(ascending=False)


# Run for each restart probability and collect the top-30 hits.
rwr_results = {}
for r in RESTART_PROBS:
    scores = rwr(G, SEED_NODES, r)
    rwr_results[r] = scores
    top30 = scores.head(30).index.tolist()
    c0_in_top30 = sum(1 for n in top30 if ground_truth[n] == "C0")
    print(f"r={r}: top-30 contains {c0_in_top30}/30 community-C0 nodes; "
          f"score range [{scores.min():.4f}, {scores.max():.4f}]")

# Visual: distribution of RWR score per community for r=0.5 (the canonical value).
r_pick = 0.5
scores = rwr_results[r_pick]
df = pd.DataFrame({
    "score": scores.values,
    "community": [ground_truth[n] for n in scores.index],
})

fig, ax = plt.subplots(figsize=(7, 3.8))
communities = sorted(df["community"].unique())
for i, c in enumerate(communities):
    vals = df.loc[df["community"] == c, "score"].values
    ax.scatter(np.full_like(vals, i) + np.random.uniform(-0.15, 0.15, size=len(vals)),
               vals, alpha=0.6, s=14, label=c)
ax.set_yscale("log")
ax.set_xticks(range(len(communities)))
ax.set_xticklabels(communities)
ax.set_xlabel("ground-truth community")
ax.set_ylabel("RWR score (log)")
ax.set_title(f"RWR scores by community (r = {r_pick}; seeds in C0)")
ax.legend(loc="upper right", fontsize=8)
plt.tight_layout()
plt.show()

# Sanity check: the median RWR score in C0 should be much higher than other communities.
median_by_comm = df.groupby("community")["score"].median()
print("\\nMedian RWR score by community:")
print(median_by_comm.round(5))
'''


STEP3_MD = """## Step 3 (15 min) — Louvain community detection vs ground truth

Louvain greedily moves nodes between communities to maximise the
**modularity**
$$Q = \\frac{1}{2m} \\sum_{ij} \\left[ A_{ij} - \\frac{k_i k_j}{2m} \\right] \\delta(c_i, c_j),$$
where $A$ is the adjacency, $k_i$ is the degree of node $i$, $m$ is the
total number of edges, and $\\delta(c_i, c_j) = 1$ when nodes $i$ and $j$
are in the same community. The maximum $Q$ in a planted-partition graph
identifies the planted communities (provided $P_{within} \\gg P_{between}$).

You will run `nx.community.louvain_communities`, compute $Q$ for the
returned partition and for the ground-truth partition, and compare via the
**adjusted Rand index** — a clustering-comparison metric that is 1 when
two partitions agree perfectly and ≈ 0 for random partitions.
"""


STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Louvain community detection.
# ----------------------------------------------------------------------


def run_louvain(G, seed: int = SEED) -> list[set]:
    """Run Louvain on G; return a list of sets, one per detected community."""
    # TODO: use nx.community.louvain_communities(G, seed=seed).
    raise NotImplementedError


def modularity(G, partition: list[set]) -> float:
    """Compute modularity Q for a given partition (list of node sets)."""
    # TODO: use nx.community.modularity.
    raise NotImplementedError


def adjusted_rand_index(labels_a: dict, labels_b: dict) -> float:
    """ARI between two node->label dicts (same key set)."""
    # TODO: build paired label arrays in a fixed node order; use
    # sklearn.metrics.adjusted_rand_score, OR implement the closed-form
    # ARI from the contingency table (more illustrative).
    raise NotImplementedError


# Run Louvain; compute Q on both partitions; print ARI.
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.
from itertools import combinations


def run_louvain(G, seed: int = SEED):
    return nx.community.louvain_communities(G, seed=seed)


def modularity(G, partition):
    return nx.community.modularity(G, partition)


def _labels_to_partition(labels: dict, all_nodes):
    """Turn {node: label} into [set(...), set(...), ...]."""
    groups = defaultdict(set)
    for n in all_nodes:
        groups[labels[n]].add(n)
    return list(groups.values())


def adjusted_rand_index(labels_a: dict, labels_b: dict) -> float:
    """Closed-form ARI from a contingency table.

    n_ij = |a_i intersect b_j|; sums over (n_ij choose 2) give the index.
    """
    nodes = list(labels_a.keys())
    a = [labels_a[n] for n in nodes]
    b = [labels_b[n] for n in nodes]

    # Contingency table.
    a_keys = sorted(set(a))
    b_keys = sorted(set(b))
    ai = {k: i for i, k in enumerate(a_keys)}
    bi = {k: i for i, k in enumerate(b_keys)}
    table = np.zeros((len(a_keys), len(b_keys)), dtype=int)
    for x, y in zip(a, b):
        table[ai[x], bi[y]] += 1

    def comb2(x):
        return x * (x - 1) // 2

    n = len(nodes)
    sum_comb_table = sum(comb2(int(v)) for v in table.flatten())
    sum_comb_a = sum(comb2(int(v)) for v in table.sum(axis=1))
    sum_comb_b = sum(comb2(int(v)) for v in table.sum(axis=0))
    expected_index = sum_comb_a * sum_comb_b / comb2(n)
    max_index = 0.5 * (sum_comb_a + sum_comb_b)
    if max_index == expected_index:
        return 1.0
    return (sum_comb_table - expected_index) / (max_index - expected_index)


# Run Louvain.
detected = run_louvain(G)
print(f"Louvain detected {len(detected)} communities; sizes = "
      f"{sorted((len(c) for c in detected), reverse=True)}")

# Ground-truth partition.
truth_partition = _labels_to_partition(ground_truth, G.nodes())

Q_detected = modularity(G, detected)
Q_truth    = modularity(G, truth_partition)

print(f"Modularity Q (Louvain)     = {Q_detected:.4f}")
print(f"Modularity Q (ground truth)= {Q_truth:.4f}")

# Convert detected partition to a {node: label} dict and compute ARI.
detected_labels = {}
for i, comm in enumerate(detected):
    for n in comm:
        detected_labels[n] = f"L{i}"

ari = adjusted_rand_index(ground_truth, detected_labels)
print(f"Adjusted Rand index (Louvain vs ground truth) = {ari:.4f}")

# Visualise the agreement as a confusion matrix.
truth_labels_list = sorted(set(ground_truth.values()))
detected_labels_list = sorted(set(detected_labels.values()))
M = np.zeros((len(truth_labels_list), len(detected_labels_list)), dtype=int)
ti = {k: i for i, k in enumerate(truth_labels_list)}
di = {k: i for i, k in enumerate(detected_labels_list)}
for n in G.nodes():
    M[ti[ground_truth[n]], di[detected_labels[n]]] += 1

fig, ax = plt.subplots(figsize=(6, 4))
im = ax.imshow(M, cmap="Blues", aspect="auto")
ax.set_xticks(range(len(detected_labels_list)))
ax.set_xticklabels(detected_labels_list, rotation=45, ha="right")
ax.set_yticks(range(len(truth_labels_list)))
ax.set_yticklabels(truth_labels_list)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if M[i, j] > 0:
            ax.text(j, i, str(M[i, j]),
                    ha="center", va="center",
                    color="white" if M[i, j] > M.max() / 2 else "black",
                    fontsize=9)
ax.set_xlabel("Louvain community")
ax.set_ylabel("Ground-truth community")
ax.set_title(f"Confusion matrix  (Q_louvain={Q_detected:.3f}, ARI={ari:.3f})")
plt.colorbar(im, ax=ax, fraction=0.04)
plt.tight_layout()
plt.show()
'''


STEP4_MD = """## Step 4 (15 min) — Hypergeometric pathway enrichment on RWR top-30

A **gene set** is a curated collection of genes that share a biological
property (a pathway, a complex, a Gene Ontology term). Given (a) the
universe $U$ of all genes (in our case, all 500 nodes), (b) the
candidate list $L$ (RWR's top 30 with $r = 0.5$), and (c) one gene set
$G$, the question is: *is the overlap $|L \\cap G|$ bigger than chance
would predict?*

Under the null that $L$ is a uniform random subset of $U$ of size
$|L|$, the count $|L \\cap G|$ follows a **hypergeometric distribution**
with parameters
- $N = |U|$ (universe size),
- $K = |G|$ (gene-set size),
- $n = |L|$ (list size).

The one-sided p-value is $P[\\text{Hyp}(N, K, n) \\geq |L \\cap G|]$. You
will test the RWR list against a small toy database of 4 gene sets — one
of which is engineered to overlap heavily with community C0 (where the
seeds live) and should pop out as the top enrichment.
"""


STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Hypergeometric enrichment on RWR top-30.
# ----------------------------------------------------------------------

# Toy gene-set database. Each pathway is a list of node names.
# - "C0_module"     : 30 nodes from community C0 (should be HIGHLY enriched).
# - "C2_module"     : 30 nodes from community C2 (should NOT be enriched).
# - "ribosome_like" : 25 random nodes drawn from the full graph (null-ish).
# - "wide_pathway"  : 60 nodes mixed across communities (mild signal at best).

def build_gene_sets(ground_truth, rng):
    """Return {pathway_name: list[node_name]}."""
    # TODO: build the 4 pathway lists described above using the rng for
    # reproducibility. C0 / C2 pathways should be drawn from their
    # respective communities; the others should be sampled from all nodes.
    raise NotImplementedError


def hypergeom_enrichment(L: set, gene_set: set, universe: set) -> dict:
    """Return {'k': k, 'K': K, 'n': n, 'N': N, 'p_value': p}.

    k = |L intersect G|, K = |G|, n = |L|, N = |U|; p = P[X >= k] under
    Hypergeometric(N, K, n).
    """
    # TODO: use scipy.stats.hypergeom.sf(k - 1, N, K, n).
    raise NotImplementedError


# Run enrichment for every pathway in the toy DB; rank by p-value.
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.


def build_gene_sets(ground_truth, rng):
    by_comm = defaultdict(list)
    for n, c in ground_truth.items():
        by_comm[c].append(n)
    for c in by_comm:
        by_comm[c].sort()

    all_nodes = sorted(ground_truth.keys())

    gene_sets = {
        # 30 nodes drawn from community C0 (where the seeds live)
        "C0_module":     list(rng.choice(by_comm["C0"], size=30, replace=False)),
        # 30 nodes from a different community — should NOT light up
        "C2_module":     list(rng.choice(by_comm["C2"], size=30, replace=False)),
        # 25 random nodes — null-ish
        "ribosome_like": list(rng.choice(all_nodes, size=25, replace=False)),
        # 60-node "wide" pathway sampling across communities — weak at best
        "wide_pathway":  list(rng.choice(all_nodes, size=60, replace=False)),
    }
    return gene_sets


def hypergeom_enrichment(L: set, gene_set: set, universe: set) -> dict:
    N = len(universe)
    K = len(gene_set & universe)
    n = len(L & universe)
    k = len(L & gene_set & universe)
    # P[X >= k] = 1 - CDF(k-1) = SF(k-1)
    p = float(stats.hypergeom.sf(k - 1, N, K, n))
    return {"k": k, "K": K, "n": n, "N": N, "p_value": p}


rng = np.random.default_rng(SEED)
gene_sets = build_gene_sets(ground_truth, rng)

# Take the RWR top-30 with r=0.5.
L = set(rwr_results[0.5].head(30).index)
universe = set(G.nodes())

rows = []
for name, gs in gene_sets.items():
    gs_set = set(gs)
    res = hypergeom_enrichment(L, gs_set, universe)
    res["pathway"] = name
    res["fold_enrichment"] = (
        (res["k"] / res["n"]) / (res["K"] / res["N"]) if res["K"] > 0 else float("nan")
    )
    rows.append(res)

enrich = pd.DataFrame(rows)[
    ["pathway", "k", "K", "n", "N", "fold_enrichment", "p_value"]
].sort_values("p_value")

# BH FDR (small list — Benjamini-Hochberg by rank).
m = len(enrich)
enrich = enrich.reset_index(drop=True)
enrich["rank"] = enrich.index + 1
enrich["p_adj"] = (enrich["p_value"] * m / enrich["rank"]).clip(upper=1.0)

print("ORA hypergeometric enrichment on RWR top-30:")
print(enrich.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

# Bar plot of -log10(p) per pathway.
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.barh(
    enrich["pathway"],
    -np.log10(np.maximum(enrich["p_value"], 1e-300)),
    color=["crimson" if p < 0.01 else "steelblue" for p in enrich["p_value"]],
)
ax.axvline(-np.log10(0.05), color="grey", linestyle="--", label="p = 0.05")
ax.axvline(-np.log10(0.01), color="black", linestyle="--", label="p = 0.01")
ax.set_xlabel("-log10(p-value)")
ax.set_title("Hypergeometric enrichment on RWR top-30 (r=0.5)")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP5_MD = """## Step 5 (10 min) — RWR is a low-pass filter on the graph Laplacian

The graph Laplacian $L = D - A$ has eigenvalues $0 = \\lambda_0 \\leq
\\lambda_1 \\leq \\dots \\leq \\lambda_{N-1}$ and eigenvectors $u_i$.
Small-$\\lambda$ eigenvectors are **smooth** across edges (slowly varying);
large-$\\lambda$ ones **oscillate** rapidly between neighbours. The graph
Fourier transform of a signal $f$ is $\\hat f_i = u_i^T f$.

For the (symmetrically normalised) Laplacian $\\mathcal{L} = I - W_{sym}$,
the closed-form RWR solution
$$\\mathbf{p} = r \\, (I - (1-r) W_{sym})^{-1} \\, \\mathbf{p}_0$$
acts in the Fourier basis as the diagonal multiplier
$$\\hat{p}_i = \\frac{r}{r + (1 - r)\\,\\lambda_i} \\, \\hat{p}_{0,i}.$$
This is a **low-pass filter** with cutoff controlled by $r$:
- $r \\to 1$: $\\hat p_i \\to \\hat p_{0,i}$ — no smoothing.
- $r \\to 0$: only $\\lambda_0 = 0$ survives — total smoothing to the
  stationary distribution.

You will compute the actual transfer function $r / (r + (1-r) \\lambda_i)$
on this graph's Laplacian spectrum and overlay it for $r \\in \\{0.3, 0.5,
0.7\\}$ — concretely showing that **restart probability = cutoff frequency**.
"""


STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — RWR as a graph low-pass filter.
# ----------------------------------------------------------------------

# TODO:
# 1. Build the symmetrically normalised Laplacian
#    L_sym = I - D^{-1/2} A D^{-1/2} as a dense numpy array (500x500 is fine).
# 2. Eigendecompose with np.linalg.eigh; sort eigenvalues ascending.
# 3. For r in {0.3, 0.5, 0.7}, plot the transfer function
#    H(lambda) = r / (r + (1-r) * lambda) vs lambda.
# 4. Mark the half-power point H(lambda) = 1/2: lambda_half = r / (1 - r).
'''


STEP5_SOLUTION = '''# Reference solution — Step 5.

nodes = list(G.nodes())
A = nx.to_numpy_array(G, nodelist=nodes)
deg = A.sum(axis=1)
d_inv_sqrt = np.where(deg > 0, 1.0 / np.sqrt(deg), 0.0)
L_sym = np.eye(len(nodes)) - (d_inv_sqrt[:, None] * A * d_inv_sqrt[None, :])

# Eigendecomposition of the symmetric Laplacian.
eigvals, _ = np.linalg.eigh(L_sym)
eigvals = np.sort(eigvals)
print(f"L_sym eigenvalues: min={eigvals[0]:.4f}, max={eigvals[-1]:.4f} "
      f"(theory: 0 <= lambda <= 2)")

# Transfer function for each restart probability.
lambdas = np.linspace(0, 2, 400)

fig, ax = plt.subplots(figsize=(7.5, 4))
for r in RESTART_PROBS:
    H = r / (r + (1.0 - r) * lambdas)
    ax.plot(lambdas, H, lw=2, label=f"r = {r}")
    lam_half = r / (1.0 - r)
    if lam_half <= 2:
        ax.axvline(lam_half, color="grey", linestyle=":", alpha=0.5)
        ax.text(lam_half, 0.55, rf"$\\lambda_{{1/2}}={lam_half:.2f}$",
                rotation=90, va="bottom", ha="right", fontsize=8)

# Overlay the actual eigenvalue spectrum as a rug.
ax.scatter(eigvals, np.full_like(eigvals, -0.04), s=3, color="black", alpha=0.4,
           clip_on=False, label="Laplacian eigenvalues")

ax.axhline(0.5, color="black", linestyle="--", alpha=0.4)
ax.set_xlim(0, 2)
ax.set_ylim(-0.05, 1.05)
ax.set_xlabel(r"Laplacian eigenvalue $\\lambda$  (graph frequency)")
ax.set_ylabel(r"$H(\\lambda) = r / (r + (1-r)\\lambda)$")
ax.set_title("RWR as a low-pass filter on the graph Laplacian")
ax.legend(loc="upper right")
plt.tight_layout()
plt.show()

# Print the half-power cutoff for each restart probability.
print("\\nHalf-power cutoff (lambda where H = 1/2):")
for r in RESTART_PROBS:
    print(f"  r = {r:.1f}  ->  lambda_1/2 = r / (1 - r) = {r / (1 - r):.3f}")

print(
    "\\nSmall r  -> low cutoff -> smoothing aggressive -> RWR spreads far from seeds.\\n"
    "Large r  -> high cutoff -> little smoothing      -> RWR sticks to seeds.\\n"
    "Restart probability r is literally the cutoff frequency of the graph filter."
)
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing pieces of the exercise. If you ran
the reference solutions above they should all pass; if you wrote your own
and an assert fails, revisit the corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. The graph has 500 nodes (or very close — we drop tiny components, never C0).
assert 480 <= G.number_of_nodes() <= 500, \\
    f"unexpected node count {G.number_of_nodes()}"
assert nx.is_connected(G), "graph should be connected after taking the largest component"

# 2. RWR returns a proper probability vector for r=0.5.
scores = rwr_results[0.5]
assert abs(scores.sum() - 1.0) < 1e-6, f"RWR scores sum to {scores.sum():.6f}, not 1"
assert (scores >= -1e-12).all(), "RWR scores should be non-negative"

# 3. Top-30 RWR list is enriched for community C0.
top30 = scores.head(30).index.tolist()
c0_in_top30 = sum(1 for n in top30 if ground_truth[n] == "C0")
assert c0_in_top30 >= 20, f"top-30 has only {c0_in_top30} C0 nodes; RWR is broken"

# 4. Louvain modularity is high (Q > 0.5 is "strong community structure").
assert Q_detected > 0.5, f"Louvain modularity {Q_detected:.3f} too low for SBM"
assert ari > 0.7, f"Louvain vs ground truth ARI {ari:.3f} too low — Louvain failed"

# 5. The C0_module pathway is the top enrichment hit by p-value.
top_pathway = enrich.iloc[0]["pathway"]
assert top_pathway == "C0_module", \\
    f"expected C0_module to top the enrichment, got {top_pathway}"
assert enrich.iloc[0]["p_value"] < 1e-3, \\
    f"C0_module p-value {enrich.iloc[0]['p_value']:.2e} larger than expected"

# 6. C2_module — same size as C0_module but distant from the seed — should NOT
#    be significant after Bonferroni or even at p < 0.01.
c2_row = enrich.loc[enrich["pathway"] == "C2_module"].iloc[0]
assert c2_row["p_value"] > 0.01, \\
    f"C2_module unexpectedly significant (p={c2_row['p_value']:.3g})"

# 7. RWR's low-pass interpretation: at r=0.5 the half-power cutoff is lambda=1.
#    Compute the closed-form transfer function and check the math.
H_at_lambda_1 = 0.5 / (0.5 + 0.5 * 1.0)
assert abs(H_at_lambda_1 - 0.5) < 1e-9, "low-pass transfer function math is off"

print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — graph signal processing in three movements

You just implemented three apparently-distinct ideas; they are all the same
story told in the graph Laplacian's eigenbasis.

1. **RWR = low-pass filter.** The closed-form solution is a diagonal
   matrix in the Fourier basis with entries $H(\\lambda_i) = r / (r + (1-r)
   \\lambda_i)$. That is a first-order low-pass filter with cutoff
   $\\lambda_{1/2} = r / (1-r)$ — exactly the role $\\omega_c$ plays in
   an RC filter. Network propagation is **heat diffusion** on a resistor
   network; the restart is a grounding resistor that bleeds the signal back
   to the seed at rate $r$.

2. **Louvain = k-means in eigenspace.** Modularity maximisation is
   equivalent (in the limit of small graphs) to projecting onto the smallest
   non-zero eigenvectors of the Laplacian and clustering there — exactly
   the spectral-clustering recipe. Both algorithms find the slow-varying
   eigenvectors and cut where they change sign.

3. **Hypergeometric enrichment = two-sample binomial.** ORA tests the null
   that the RWR top-$n$ is a uniform random subset of the universe; the
   answer is the same as the urn problem that produces every "draws
   without replacement" question in introductory probability. GSEA (which
   we did not implement) replaces "is the count high?" with "is the rank
   distribution shifted?" — a two-sample Mann-Whitney or KS test.

The unifying observation: a biological network is just a finite-dimensional
signal-processing domain whose Fourier basis is the Laplacian's
eigenvectors. Every standard EE filter (low-pass, high-pass, matched filter)
has a graph analogue; every clustering method you know already has a
graph-spectrum implementation. Bioinformatics calls this "network biology";
the EE department calls it **graph signal processing**.
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
        md(STEP1_SOLUTION_HEADER),
        hidden(STEP1_SOLUTION),

        md(STEP2_MD),
        code(STEP2_TODO),
        md(STEP1_SOLUTION_HEADER),
        hidden(STEP2_SOLUTION),

        md(STEP3_MD),
        code(STEP3_TODO),
        md(STEP1_SOLUTION_HEADER),
        hidden(STEP3_SOLUTION),

        md(STEP4_MD),
        code(STEP4_TODO),
        md(STEP1_SOLUTION_HEADER),
        hidden(STEP4_SOLUTION),

        md(STEP5_MD),
        code(STEP5_TODO),
        md(STEP1_SOLUTION_HEADER),
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
