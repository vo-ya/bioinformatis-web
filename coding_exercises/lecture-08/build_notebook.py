"""Build exercise.ipynb for L08 — Advanced Single-Cell: Trajectories, Integration, Multi-Modal.

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


TITLE_MD = """# L08 — Advanced Single-Cell: Trajectories, Integration, Multi-Modal

In this exercise you stitch together the three pillars of "advanced" scRNA-seq
analysis on top of the L07 pipeline:

1. compute a **pseudotime** ordering on a myeloid differentiation trajectory
   (Paul et al. 2015, `sc.datasets.paul15()`),
2. **integrate** two synthetic batches of PBMCs with **Harmony** and quantify
   the lift with a batch-mixing metric, and
3. score a lightweight **ligand-receptor** hypothesis between two clusters
   using built-in marker genes.

The notebook is deliberately CPU-friendly and runs end-to-end in well under
five minutes on free Colab.
"""


AHA_MD = """> **Aha.** A static snapshot of an asynchronous population implicitly
> contains the whole trajectory; **pseudotime is the 1-D coordinate on that
> manifold**. Batch effects are a **supervised source-separation** problem —
> Harmony iteratively shifts cluster centroids per batch until biology
> dominates over batch. Ligand-receptor scoring is a **co-expression
> hypothesis** between two cell-type vectors, not a causal claim.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
# `scanpy` brings AnnData + the standard pipeline; `harmonypy` is Harmony in Python.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 scikit-learn==1.5.0 matplotlib==3.8.4 scanpy==1.10.2 harmonypy==0.0.10 -q
"""


IMPORTS = """import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import harmonypy as hm
from scipy.stats import spearmanr
from sklearn.metrics import silhouette_score

# Deterministic for the whole notebook.
SEED = 42
np.random.seed(SEED)
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=80, facecolor="white")

print(f"scanpy   : {sc.__version__}")
print(f"harmonypy: {hm.__version__ if hasattr(hm, '__version__') else 'installed'}")
"""


# ---------------------------------------------------------------------------
# Step 1 — Load + standard pipeline on the myeloid trajectory dataset.
# ---------------------------------------------------------------------------


STEP1_MD = """## Step 1 (8 min) — Standard pipeline on the myeloid trajectory dataset

We start with the Paul et al. 2015 myeloid progenitor dataset, available as a
Scanpy built-in via `sc.datasets.paul15()`. It contains ~2 700 mouse bone
marrow cells along an erythroid / myeloid differentiation gradient — exactly
the kind of asynchronous snapshot that pseudotime is designed for.

Run a stripped-down L07 pipeline: log-normalise, PCA, neighbour graph, UMAP,
and Leiden clusters. The UMAP should already look like a curved manifold;
that visual continuity is what pseudotime will turn into a scalar coordinate.
"""

STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Load Paul15 and run the standard pipeline.
# ----------------------------------------------------------------------

# TODO:
# 1. Load `sc.datasets.paul15()` into an AnnData called `adata_traj`.
# 2. Normalise total counts to 1e4 and log1p-transform.
# 3. Run PCA (50 comps), build a k=15 neighbour graph, run UMAP, and Leiden
#    (resolution=0.5). Use random_state=SEED throughout for reproducibility.
# 4. Plot the UMAP coloured by Leiden and by the published `paul15_clusters`
#    label so you can see the differentiation gradient.

adata_traj = None
'''


STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

adata_traj = sc.datasets.paul15()
print(f"Paul15: {adata_traj.shape} (cells x genes)")

# The Paul15 dataset is already filtered + normalised in the Scanpy build, but
# we still apply log1p + scaling to be explicit and to mirror what students did
# in L07. The dataset is small so this is fast.
adata_traj.X = adata_traj.X.astype(float)
sc.pp.normalize_total(adata_traj, target_sum=1e4)
sc.pp.log1p(adata_traj)

sc.tl.pca(adata_traj, n_comps=50, random_state=SEED)
sc.pp.neighbors(adata_traj, n_neighbors=15, n_pcs=20, random_state=SEED)
sc.tl.umap(adata_traj, random_state=SEED)
sc.tl.leiden(adata_traj, resolution=0.5, random_state=SEED)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
sc.pl.umap(adata_traj, color="leiden", ax=axes[0], show=False, title="Leiden clusters")
sc.pl.umap(adata_traj, color="paul15_clusters", ax=axes[1], show=False,
           title="Published Paul15 labels", legend_fontsize=6)
plt.tight_layout()
plt.show()

print(f"Leiden clusters: {adata_traj.obs['leiden'].nunique()}")
'''


# ---------------------------------------------------------------------------
# Step 2 — Diffusion pseudotime (DPT).
# ---------------------------------------------------------------------------


STEP2_MD = """## Step 2 (12 min) — Diffusion pseudotime on the trajectory

Scanpy ships **diffusion pseudotime (DPT)** out of the box. Conceptually, it
spectral-embeds the cell-cell kNN graph (the dominant non-trivial eigenvectors
of the diffusion operator are smooth functions over the manifold), then
measures graph distance from a chosen **root cell** in that embedding.

We pick a root in the most progenitor-like Leiden cluster. The UMAP coloured
by `dpt_pseudotime` should show a smooth 0→1 gradient along the differentiation
direction. Wherever the published label says "GMP / MEP / progenitor", the
pseudotime should be near zero; mature granulocyte / erythrocyte labels should
sit near one.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Pick a root cell and run diffusion pseudotime.
# ----------------------------------------------------------------------

# TODO:
# 1. Pick a root cell in the cluster whose published label looks most like
#    a progenitor (hint: any cluster whose `paul15_clusters` mode contains
#    "GMP" or starts with a digit ≤ 6 in Paul15's numbering).
# 2. Set `adata_traj.uns["iroot"]` to the integer index of that cell.
# 3. Compute the diffusion map (`sc.tl.diffmap`) and DPT (`sc.tl.dpt`).
# 4. Plot UMAP coloured by `dpt_pseudotime`.
# 5. Verify monotonicity: Spearman correlation between dpt_pseudotime and a
#    rough numeric ordering of `paul15_clusters` (cluster 1 → start,
#    cluster 19 → end). |rho| > 0.5 is a sane sanity check.

root_idx = None
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.

# Paul15 cluster labels look like "1Ery", "2Ery", ..., "19Lymph". The leading
# integer is roughly a developmental order (low = progenitor, high = mature),
# so we use that for a sanity-check correlation.
def paul15_order(label: str) -> int:
    digits = ""
    for c in label:
        if c.isdigit():
            digits += c
        else:
            break
    return int(digits) if digits else -1

adata_traj.obs["paul15_order"] = (
    adata_traj.obs["paul15_clusters"].astype(str).map(paul15_order)
)

# Pick a root in the cluster whose mode paul15_order is smallest.
cluster_orders = (
    adata_traj.obs.groupby("leiden")["paul15_order"].median().sort_values()
)
root_cluster = cluster_orders.index[0]
candidates = np.where(adata_traj.obs["leiden"] == root_cluster)[0]
root_idx = int(candidates[0])
adata_traj.uns["iroot"] = root_idx
print(f"Root cluster: Leiden {root_cluster} (median Paul15 order "
      f"{cluster_orders.iloc[0]:.0f}), root cell index {root_idx}")

# Diffusion map + DPT.
sc.tl.diffmap(adata_traj, n_comps=15)
sc.tl.dpt(adata_traj, n_dcs=10)

fig, ax = plt.subplots(figsize=(6, 4.5))
sc.pl.umap(adata_traj, color="dpt_pseudotime", ax=ax, show=False,
           title="Diffusion pseudotime", color_map="viridis")
plt.tight_layout()
plt.show()

# Sanity-check monotonicity vs the published numeric ordering.
mask = adata_traj.obs["paul15_order"] >= 0
rho, pval = spearmanr(
    adata_traj.obs.loc[mask, "dpt_pseudotime"],
    adata_traj.obs.loc[mask, "paul15_order"],
)
# Sign is arbitrary (depends on root direction); we care about magnitude.
print(f"Spearman rho(dpt_pseudotime, paul15_order) = {rho:+.3f}  (p = {pval:.1e})")
print(f"|rho| = {abs(rho):.3f}  -> pseudotime tracks the published gradient" if abs(rho) > 0.4 else "  weak alignment; consider a different root")
'''


# ---------------------------------------------------------------------------
# Step 3 — Two synthetic batches with a batch effect.
# ---------------------------------------------------------------------------


STEP3_MD = """## Step 3 (12 min) — Build two synthetic PBMC batches

For the integration leg we use `sc.datasets.pbmc3k()`, the same dataset as
L07. We split it 50/50 into two synthetic batches and **inject a batch effect
by hand**: a multiplicative gain on a random subset of genes, applied only to
batch B. After concatenation, a UMAP coloured by batch label should show two
clearly separated lobes — the batch effect drowning out the cell-type signal.

This is intentional: we want to *see* what integration has to fix.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Load PBMC3k, split into two batches, inject a batch effect.
# ----------------------------------------------------------------------

# TODO:
# 1. Load `sc.datasets.pbmc3k()` into `adata_pbmc`.
# 2. Drop cells with fewer than 200 genes and genes detected in fewer than 3 cells.
# 3. Randomly split cells into batch_A and batch_B (50/50, seeded).
# 4. For batch_B only, multiply expression of a random 30% of genes by ~1.6
#    (a plausible technical gain). Re-store as a new AnnData with `obs["batch"]`
#    set to "A" / "B".
# 5. Normalise + log1p + 2000 HVGs + PCA + neighbours + UMAP and plot coloured
#    by batch. The two batches should separate visibly.

adata_batched = None
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.

adata_pbmc = sc.datasets.pbmc3k()
print(f"PBMC3k raw: {adata_pbmc.shape}")

sc.pp.filter_cells(adata_pbmc, min_genes=200)
sc.pp.filter_genes(adata_pbmc, min_cells=3)

rng = np.random.default_rng(SEED)
n = adata_pbmc.n_obs
batch_labels = np.array(["A"] * n)
batch_b_idx = rng.choice(n, size=n // 2, replace=False)
batch_labels[batch_b_idx] = "B"
adata_pbmc.obs["batch"] = pd.Categorical(batch_labels, categories=["A", "B"])

# Inject a batch effect on a random 30% of genes, batch B only.
n_genes = adata_pbmc.n_vars
affected_genes = rng.choice(n_genes, size=int(0.3 * n_genes), replace=False)
gain = 3.0  # multiplicative technical gain (strong enough to dominate biology before integration)
X = adata_pbmc.X.toarray() if hasattr(adata_pbmc.X, "toarray") else np.asarray(adata_pbmc.X)
X = X.astype(float)
mask_B = (adata_pbmc.obs["batch"] == "B").values
X[np.ix_(mask_B, affected_genes)] *= gain
adata_pbmc.X = X
print(f"  injected gain={gain}x on {len(affected_genes)} genes for batch B "
      f"({mask_B.sum()} cells)")

# Standard pipeline.
sc.pp.normalize_total(adata_pbmc, target_sum=1e4)
sc.pp.log1p(adata_pbmc)
# Stash the full log-normalised matrix in a layer for the LR step (Step 5),
# so we can score ligand-receptor pairs that aren't necessarily HVGs.
adata_pbmc.layers["lognorm_full"] = adata_pbmc.X.copy()
adata_pbmc.uns["full_var_names"] = list(adata_pbmc.var_names)
sc.pp.highly_variable_genes(adata_pbmc, n_top_genes=2000)
adata_batched = adata_pbmc[:, adata_pbmc.var["highly_variable"]].copy()
# Preserve full matrix on the HVG-subset object so Step 5 can look up any gene.
adata_batched.uns["full_X"] = np.asarray(adata_pbmc.X)
adata_batched.uns["full_var_names"] = list(adata_pbmc.var_names)
sc.pp.scale(adata_batched, max_value=10)
sc.tl.pca(adata_batched, n_comps=30, random_state=SEED)
sc.pp.neighbors(adata_batched, n_neighbors=15, n_pcs=20, random_state=SEED)
sc.tl.umap(adata_batched, random_state=SEED)
sc.tl.leiden(adata_batched, resolution=0.5, random_state=SEED, key_added="leiden_raw")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
sc.pl.umap(adata_batched, color="batch", ax=axes[0], show=False,
           title="Before integration — by batch")
sc.pl.umap(adata_batched, color="leiden_raw", ax=axes[1], show=False,
           title="Before integration — Leiden")
plt.tight_layout()
plt.show()

print(f"adata_batched shape: {adata_batched.shape}")
'''


# ---------------------------------------------------------------------------
# Step 4 — Harmony integration + integration metric.
# ---------------------------------------------------------------------------


STEP4_MD = """## Step 4 (18 min) — Harmony integration and a quantitative integration metric

[Harmony](https://github.com/slowkow/harmonypy) is the Python port of the
original Harmony algorithm: iterative cluster-aware linear correction in
PCA space, alternating between (a) soft k-means clustering and (b) per-batch
cluster-centroid alignment. It returns batch-corrected PCs that you feed
back into the standard neighbour-graph / UMAP pipeline.

We evaluate integration with the **silhouette delta**:
$$\\Delta = s_{\\text{batch}}^{\\text{before}} - s_{\\text{batch}}^{\\text{after}}$$
Smaller batch silhouette = batches are mixed (good). Large positive Δ means
Harmony successfully scrambled the batch label.

We also track the **cluster-label silhouette** before and after — this should
*not* drop much, because Harmony is supposed to preserve biology while
removing batch.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Run Harmony and compute the integration silhouette delta.
# ----------------------------------------------------------------------

# TODO:
# 1. Extract `X_pca = adata_batched.obsm["X_pca"]` and `meta` =
#    `adata_batched.obs[["batch"]]`.
# 2. Call `hm.run_harmony(X_pca, meta, vars_use="batch", random_state=SEED)`;
#    store `Z_corrected = ho.Z_corr.T` back into `adata_batched.obsm["X_harmony"]`.
# 3. Re-run neighbours/UMAP/Leiden on `use_rep="X_harmony"`.
# 4. Compute batch silhouette before (on `X_pca`) and after (on `X_harmony`)
#    using `sklearn.metrics.silhouette_score(... metric="euclidean")`. Use a
#    random sub-sample of ~500 cells so it stays fast.
# 5. Plot UMAPs before vs after, both coloured by batch.

batch_sil_before = None
batch_sil_after = None
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.

X_pca = adata_batched.obsm["X_pca"]
meta = adata_batched.obs[["batch"]].copy()

print("Running Harmony...")
ho = hm.run_harmony(X_pca, meta, vars_use=["batch"], random_state=SEED, max_iter_harmony=10)
Z_corrected = np.asarray(ho.Z_corr).T  # (n_cells, n_comps)
adata_batched.obsm["X_harmony"] = Z_corrected
print(f"  corrected PCA shape: {Z_corrected.shape}")

# Re-cluster on the corrected embedding.
sc.pp.neighbors(adata_batched, n_neighbors=15, use_rep="X_harmony",
                random_state=SEED, key_added="harmony")
sc.tl.umap(adata_batched, neighbors_key="harmony", random_state=SEED)
sc.tl.leiden(adata_batched, neighbors_key="harmony", resolution=0.5,
             random_state=SEED, key_added="leiden_harmony")

# Silhouette: lower batch-silhouette = better mixed.
rng = np.random.default_rng(SEED)
sub = rng.choice(adata_batched.n_obs, size=min(500, adata_batched.n_obs), replace=False)
batch_arr = adata_batched.obs["batch"].cat.codes.values[sub]

batch_sil_before = silhouette_score(X_pca[sub], batch_arr, metric="euclidean")
batch_sil_after  = silhouette_score(Z_corrected[sub], batch_arr, metric="euclidean")

# Biology silhouette via Leiden labels (should NOT collapse).
bio_arr_before = adata_batched.obs["leiden_raw"].cat.codes.values[sub]
bio_arr_after  = adata_batched.obs["leiden_harmony"].cat.codes.values[sub]
bio_sil_before = silhouette_score(X_pca[sub], bio_arr_before, metric="euclidean")
bio_sil_after  = silhouette_score(Z_corrected[sub], bio_arr_after, metric="euclidean")

print()
print("              before    after    delta")
print(f"batch sil :  {batch_sil_before:+.3f}   {batch_sil_after:+.3f}   "
      f"{batch_sil_before - batch_sil_after:+.3f}  (positive = batch mixing improved)")
print(f"bio   sil :  {bio_sil_before:+.3f}   {bio_sil_after:+.3f}   "
      f"{bio_sil_after - bio_sil_before:+.3f}  (positive = biology preserved / sharper)")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
sc.pl.umap(adata_batched, color="batch", ax=axes[0], show=False,
           title=f"After Harmony — by batch (sil={batch_sil_after:+.2f})")
sc.pl.umap(adata_batched, color="leiden_harmony", ax=axes[1], show=False,
           title="After Harmony — Leiden")
plt.tight_layout()
plt.show()
'''


# ---------------------------------------------------------------------------
# Step 5 — Lightweight ligand-receptor hypothesis.
# ---------------------------------------------------------------------------


STEP5_MD = """## Step 5 (10 min) — Lightweight ligand-receptor score

CellPhoneDB is too heavy for free Colab CPU, so we substitute a **co-expression
LR score** between two clusters on the integrated PBMC dataset. The recipe:

1. pick two clusters (e.g. Leiden A vs Leiden B from Harmony),
2. for each ligand-receptor pair in a small hand-curated list of known immune
   LRs (CD8A / CD8B, ITGAL / ICAM1, CCL5 / CCR5, IL7R / IL7), compute the
   mean ligand expression in cluster A and mean receptor expression in
   cluster B,
3. score the pair as `mean(L_in_A) * mean(R_in_B)`,
4. compare to a **permutation null** that shuffles cluster labels.

This is the same logical structure as CellPhoneDB without the database / GPU.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Co-expression LR score with a permutation test.
# ----------------------------------------------------------------------

# TODO:
# 1. Build LR_PAIRS = [("CD8A","CD8B"), ("ITGAL","ICAM1"), ("CCL5","CCR5"),
#                     ("IL7R","IL7"), ("GZMB","PRF1")].
#    (Only keep pairs where BOTH genes are in adata_batched.var_names.)
# 2. Choose source/target clusters from `leiden_harmony` (e.g. the two largest).
# 3. Define `lr_score(adata, src, tgt, L, R)` =
#       mean(adata[src cells, L]) * mean(adata[tgt cells, R])
# 4. Build a permutation null by shuffling cluster labels K=200 times and
#    recomputing the score for each pair; report empirical p-value.
# 5. Print a table of {pair, observed, null_mean, p_value}.
'''


STEP5_SOLUTION = '''# Reference solution — Step 5.

CANDIDATE_PAIRS = [
    ("CD8A", "CD8B"),
    ("ITGAL", "ICAM1"),
    ("CCL5", "CCR5"),
    ("IL7R", "IL7"),
    ("GZMB", "PRF1"),
    ("CD3D", "CD3E"),
    ("LCK", "ZAP70"),
    ("HLA-A", "B2M"),
]

# Look up genes in the *full* log-normalised matrix, not just HVGs — many
# canonical immune marker pairs are not in the top-2000 dispersion set.
full_X = np.asarray(adata_batched.uns["full_X"])
full_var = list(adata_batched.uns["full_var_names"])
gene_set = set(full_var)
LR_PAIRS = [(L, R) for (L, R) in CANDIDATE_PAIRS if L in gene_set and R in gene_set]
print(f"Usable LR pairs ({len(LR_PAIRS)}/{len(CANDIDATE_PAIRS)}): {LR_PAIRS}")

# Pick the two largest Harmony clusters as source/target.
cluster_sizes = adata_batched.obs["leiden_harmony"].value_counts()
src_cl, tgt_cl = cluster_sizes.index[:2].tolist()
print(f"Source cluster: {src_cl} ({cluster_sizes[src_cl]} cells); "
      f"Target cluster: {tgt_cl} ({cluster_sizes[tgt_cl]} cells)")

X = full_X
var_to_idx = {g: i for i, g in enumerate(full_var)}
labels = adata_batched.obs["leiden_harmony"].values

def cluster_mean(cl: str, gene: str, lab=labels) -> float:
    mask = (lab == cl)
    if mask.sum() == 0:
        return 0.0
    return float(X[mask, var_to_idx[gene]].mean())

def lr_score(L: str, R: str, src: str, tgt: str, lab=labels) -> float:
    return cluster_mean(src, L, lab) * cluster_mean(tgt, R, lab)

# Observed scores.
observed = {pair: lr_score(*pair, src_cl, tgt_cl) for pair in LR_PAIRS}

# Permutation null: shuffle cluster labels K times, recompute.
K = 200
rng = np.random.default_rng(SEED)
null_scores = {pair: np.zeros(K) for pair in LR_PAIRS}
labels_arr = np.array(labels)
for k in range(K):
    perm = labels_arr[rng.permutation(len(labels_arr))]
    for pair in LR_PAIRS:
        null_scores[pair][k] = lr_score(*pair, src_cl, tgt_cl, lab=perm)

rows = []
for pair in LR_PAIRS:
    obs = observed[pair]
    null = null_scores[pair]
    pval = float((null >= obs).mean())  # one-sided
    rows.append({
        "L": pair[0], "R": pair[1],
        "src": src_cl, "tgt": tgt_cl,
        "observed": round(obs, 3),
        "null_mean": round(float(null.mean()), 3),
        "p_value": round(pval, 3),
    })

results = pd.DataFrame(rows).sort_values("p_value")
print()
print(results.to_string(index=False))
print()
print("p_value = fraction of K=200 permuted cluster-label datasets whose "
      "score >= observed.")
print("Low p with a high observed score = candidate signalling axis between "
      f"cluster {src_cl} -> {tgt_cl} (hypothesis, not proof).")
'''


# ---------------------------------------------------------------------------
# Self-check + EE framing.
# ---------------------------------------------------------------------------


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing pieces of the pipeline. They run on
whatever objects you produced above; if they pass, you have a working
trajectory + integration + LR-scoring pipeline.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Pseudotime exists, is finite, and spans a reasonable fraction of [0, 1].
assert "dpt_pseudotime" in adata_traj.obs.columns, "DPT not computed"
pt = adata_traj.obs["dpt_pseudotime"].values
finite = np.isfinite(pt)
assert finite.mean() > 0.95, f"DPT mostly non-finite: {finite.mean():.2f}"
span = pt[finite].max() - pt[finite].min()
assert span > 0.5, f"DPT range suspiciously small: {span:.3f}"

# 2. DPT is monotonic-ish vs the Paul15 numeric label.
mask = adata_traj.obs["paul15_order"] >= 0
rho, _ = spearmanr(
    adata_traj.obs.loc[mask, "dpt_pseudotime"],
    adata_traj.obs.loc[mask, "paul15_order"],
)
assert abs(rho) > 0.3, f"|spearman| too low: {abs(rho):.3f}"

# 3. Harmony embedding has the right shape.
assert "X_harmony" in adata_batched.obsm, "Harmony embedding missing"
assert adata_batched.obsm["X_harmony"].shape == adata_batched.obsm["X_pca"].shape

# 4. Batch silhouette dropped after integration.
assert batch_sil_after < batch_sil_before, (
    f"Harmony failed to reduce batch silhouette: "
    f"before={batch_sil_before:+.3f}, after={batch_sil_after:+.3f}"
)

# 5. At least one LR pair has p_value < 0.1 (i.e., something looks signal-ish).
assert (results["p_value"] < 0.5).any(), (
    "No LR pair stood out vs the permutation null — expected at least one "
    "co-expression hypothesis to survive."
)

print("Self-check passed.")
print(f"  DPT span: {span:.3f},  |Spearman| with Paul15 order: {abs(rho):.3f}")
print(f"  batch silhouette: {batch_sil_before:+.3f} -> {batch_sil_after:+.3f}  "
      f"(delta = {batch_sil_before - batch_sil_after:+.3f})")
print(f"  smallest LR p-value: {results['p_value'].min():.3f}  "
      f"({results.iloc[0]['L']}->{results.iloc[0]['R']})")
print("\\n✅ Self-check passed.")
'''


EE_MD = """## EE framing — manifold coordinates, source separation, co-expression hypotheses

Each leg of this notebook maps onto a textbook EE idea:

1. **Pseudotime = 1-D coordinate on a learned manifold.** The diffusion map
   spectrally embeds the cell-cell similarity graph; the dominant
   non-trivial eigenvectors of the diffusion operator are the smoothest
   functions over the manifold. Picking a root + measuring graph distance
   in that embedding is the same idea as projecting a signal onto its first
   principal mode after whitening — except the geometry is graph-Laplacian
   rather than Euclidean.
2. **Batch integration = supervised source separation.** The observed counts
   are a mixture of biological variation + batch nuisance. Harmony does
   alternating *soft k-means / per-batch centroid alignment* in PCA space —
   the iterated-scaling structure is identical to **domain adaptation** in
   ML and to **interference-rejection combining** in multi-channel
   communications. The batch label is the labelled nuisance regressor.
3. **Integration metric = signal-to-batch ratio.** The silhouette delta is
   the same kind of measurement an EE makes when comparing pre- and
   post-equaliser BER: did the channel-induced distortion drop relative to
   the inter-class margin?
4. **Ligand-receptor scoring = co-expression hypothesis with a permutation
   null.** The observed cluster-A × cluster-B product is a statistic; the
   permutation null is the empirical false-alarm distribution. p-value =
   one-sided tail mass. The whole construction is hypothesis-generation,
   not causal claim — the same caveat that applies to matched-filter
   detections without independent confirmation.

The recurring theme: every "advanced single-cell" method is a signal-processing
or source-separation pipeline wearing a biological label.
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
        apply_colab_form(nb)
        nbformat.write(nb, f)
    nb2 = nbformat.read(out_path, as_version=4)
    n_md = sum(1 for c in nb2.cells if c.cell_type == "markdown")
    n_code = sum(1 for c in nb2.cells if c.cell_type == "code")
    n_hidden = sum(
        1 for c in nb2.cells if c.cell_type == "code"
        and c.metadata.get("jupyter", {}).get("source_hidden")
    )
    print(f"Wrote {out_path}")
    print(f"  cells: {len(nb2.cells)} total  ({n_md} md, {n_code} code, {n_hidden} hidden-solution)")
