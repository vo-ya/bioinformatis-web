"""Build exercise.ipynb for L07 — Single-Cell RNA-seq Fundamentals.

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


TITLE_MD = """# L07 — Single-Cell RNA-seq Fundamentals

In this exercise you run the canonical scRNA-seq pipeline end-to-end on the
3k-cell PBMC dataset: QC filtering, log-normalisation, highly variable gene
(HVG) selection, PCA, UMAP, Leiden community detection, and Wilcoxon marker
gene discovery. The same five steps are what every published scRNA-seq paper
runs as its first figure — and they are nothing more than a linear-algebra
stack glued together with a graph clusterer.
"""


AHA_MD = """> **Aha.** A scRNA-seq cell is a **sparse vector in 20 000-dim gene
> space**. Most of the variance lives on a low-dimensional manifold; PCA
> gets you the best linear approximation, UMAP wrings out the residual
> nonlinear structure, and Leiden carves the resulting kNN graph into
> communities. Discrete cell types fall out as clusters because biology
> built them that way — distinct transcriptional programmes encode
> distinct cell states.
"""


PREAMBLE = """# Install pinned scientific stack on first run. scanpy pulls anndata + leiden + umap.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 scikit-learn==1.5.0 matplotlib==3.8.4 scanpy==1.10.2 leidenalg==0.10.2 igraph==0.11.6 -q
"""


IMPORTS = """import io
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp

import scanpy as sc

# Deterministic for the whole notebook.
SEED = 42
np.random.seed(SEED)
sc.settings.verbosity = 1            # 0=silent, 1=brief, 3=very chatty
sc.settings.set_figure_params(dpi=80, frameon=False)

# Quiet down a few cosmetic warnings from scanpy / anndata / numba.
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

print(f"scanpy {sc.__version__}")
"""


STEP1_MD = """## Step 1 (8 min) — Load PBMC3k and inspect QC metrics

Scanpy ships the canonical 3k-cell PBMC dataset (10x Genomics, healthy
donor) as `sc.datasets.pbmc3k()`. It is a sparse count matrix of shape
`(n_cells, n_genes)` — about 2 700 cells × 32 738 genes. We compute three
per-cell QC metrics:

- `n_counts` (or `total_counts`) — total UMI count per cell. Empty
  droplets / dying cells have low counts; doublets have high counts.
- `n_genes_by_counts` — number of genes with at least one UMI. Tightly
  correlated with total counts; the slope deviates for ambient-RNA cells.
- `pct_counts_mt` — fraction of UMIs mapping to mitochondrial genes
  (gene names start with `MT-`). High mt% = stressed or dying cell.

If `sc.datasets.pbmc3k()` cannot reach its CDN, we fall back to a tiny
deterministic synthetic count matrix that exercises the exact same
pipeline (just with smaller cluster separation).
"""

STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Load PBMC3k and compute QC metrics.
# ----------------------------------------------------------------------

def load_pbmc3k():
    """Return an AnnData object. Falls back to a synthetic matrix on failure."""
    # TODO:
    #  - Try sc.datasets.pbmc3k(); return the AnnData on success.
    #  - On any exception, build a tiny synthetic AnnData with two cell types
    #    so the rest of the pipeline still runs.
    raise NotImplementedError


def annotate_qc(adata):
    """Add var['mt'] flag and per-cell QC metrics in place."""
    # TODO:
    #  - Flag mitochondrial genes: adata.var["mt"] = var_names.str.startswith("MT-").
    #  - Run sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True).
    raise NotImplementedError


# Plot the three QC distributions as side-by-side violins.
def plot_qc(adata):
    """Three QC violin plots: total_counts, n_genes_by_counts, pct_counts_mt."""
    # TODO: use plt.subplots(1, 3) and ax.violinplot(...) for each metric.
    raise NotImplementedError


# adata = load_pbmc3k()
# annotate_qc(adata)
# plot_qc(adata)
'''


STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.
import anndata as ad


def _synthetic_pbmc(n_cells: int = 600, n_genes: int = 1500, seed: int = 0):
    """Tiny deterministic two-cluster count matrix. Same schema as 10x output."""
    rng = np.random.default_rng(seed)
    # Two cell types of equal size.
    half = n_cells // 2
    labels = np.array(["A"] * half + ["B"] * (n_cells - half))

    # Background expression: NB-ish counts with mean ~ 0.3 per gene.
    bg = rng.poisson(0.3, size=(n_cells, n_genes)).astype(np.float32)

    # 30 cluster-A marker genes, 30 cluster-B marker genes.
    a_markers = rng.choice(n_genes, size=30, replace=False)
    b_markers = rng.choice([g for g in range(n_genes) if g not in a_markers],
                           size=30, replace=False)
    for g in a_markers:
        bg[labels == "A", g] += rng.poisson(8, size=half)
    for g in b_markers:
        bg[labels == "B", g] += rng.poisson(8, size=n_cells - half)

    # Add 10 MT- genes with a few percent of counts each.
    mt_idx = np.arange(n_genes - 10, n_genes)
    var_names = [f"GENE{i:05d}" for i in range(n_genes - 10)] + \\
                [f"MT-{i:02d}" for i in range(10)]
    for g in mt_idx:
        bg[:, g] += rng.poisson(0.5, size=n_cells)

    obs = pd.DataFrame({"true_label": labels},
                       index=[f"cell{i:04d}" for i in range(n_cells)])
    var = pd.DataFrame(index=var_names)
    adata = ad.AnnData(X=sp.csr_matrix(bg), obs=obs, var=var)
    return adata


def load_pbmc3k():
    try:
        adata = sc.datasets.pbmc3k()
        print(f"  pbmc3k from scanpy CDN: {adata.shape}")
        return adata
    except Exception as exc:
        print(f"  pbmc3k fetch failed ({exc!r}); using synthetic fallback.")
        return _synthetic_pbmc()


def annotate_qc(adata):
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True,
    )


def plot_qc(adata):
    metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
    for ax, m in zip(axes, metrics):
        vals = np.asarray(adata.obs[m].values, dtype=float)
        parts = ax.violinplot(vals, showmedians=True, widths=0.7)
        for pc in parts["bodies"]:
            pc.set_facecolor("#4c72b0")
            pc.set_alpha(0.7)
        ax.set_title(m)
        ax.set_xticks([])
        ax.set_ylabel(m)
    fig.suptitle("Per-cell QC distributions (pre-filter)")
    plt.tight_layout()
    plt.show()


adata = load_pbmc3k()
annotate_qc(adata)
print(f"adata: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"  median total_counts        = {np.median(adata.obs['total_counts']):.0f}")
print(f"  median n_genes_by_counts   = {np.median(adata.obs['n_genes_by_counts']):.0f}")
print(f"  median pct_counts_mt       = {np.median(adata.obs['pct_counts_mt']):.2f}%")
plot_qc(adata)
'''


STEP2_MD = """## Step 2 (10 min) — Drop low-quality cells

We trim cells with extreme QC values. The cutoffs in the spec
(`pct_mt > 15`, `total_counts < 500` or `> 6000`) are deliberate:

- The lower count floor (500) removes empty / dying droplets.
- The upper count ceiling (6000) trims likely doublets.
- The mt% ceiling (15 %) flags stressed cells whose cytoplasmic mRNA has
  leaked out, leaving mitochondrial mRNA over-represented.

Lab-grade pipelines often tune these per-experiment from the QC violins
above; for PBMC3k these defaults work fine.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Filter low-quality cells (and ultra-rare genes).
# ----------------------------------------------------------------------

MT_MAX_PCT     = 15.0   # cells with > 15% mt counts are dying — drop.
COUNT_MIN      = 500    # < 500 UMIs = ambient / empty droplet.
COUNT_MAX      = 6000   # > 6000 UMIs = likely doublet.
MIN_CELLS_GENE = 3      # drop genes seen in < 3 cells.


def filter_cells(adata):
    """Apply the four cutoffs above. Returns a *new* AnnData (subset)."""
    # TODO:
    #  - mask = (pct_counts_mt < MT_MAX_PCT) & (total_counts > COUNT_MIN) & (...)
    #  - subset adata in place or return adata[mask].copy()
    #  - also drop genes in < MIN_CELLS_GENE cells via sc.pp.filter_genes.
    raise NotImplementedError


# Print before/after sizes so you can see how many cells were dropped.
'''

STEP2_SOLUTION = '''# Reference solution — Step 2.

MT_MAX_PCT     = 15.0
COUNT_MIN      = 500
COUNT_MAX      = 6000
MIN_CELLS_GENE = 3


def filter_cells(adata):
    n0 = adata.shape[0]
    mask = (
        (adata.obs["pct_counts_mt"] < MT_MAX_PCT)
        & (adata.obs["total_counts"] > COUNT_MIN)
        & (adata.obs["total_counts"] < COUNT_MAX)
    )
    adata = adata[mask].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_GENE)
    print(f"  cells: {n0} -> {adata.shape[0]} (kept {adata.shape[0]/n0:.1%})")
    print(f"  genes: kept {adata.shape[1]} with >= {MIN_CELLS_GENE} expressing cells")
    return adata


adata = filter_cells(adata)
'''


STEP3_MD = """## Step 3 (15 min) — Log-normalise and pick highly variable genes

Raw UMI counts are sequencing-depth-dependent. The standard transform is:

1. **CP10k**: divide each cell by its total count, multiply by 10 000 — so
   every cell sums to 10 000 transcripts. Equivalent to TPM at single-cell
   scale.
2. **log1p**: take `log(x + 1)` to stabilise variance and squash the long
   right tail of highly expressed genes.

After normalisation we pick the top **2 000 highly variable genes** by
normalised dispersion (variance / mean after binning by expression
strength). HVGs capture the biological signal; the other ~28 000 genes are
mostly housekeeping or noise — keeping them in PCA would dilute the signal
with thousands of constant directions.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Log-normalise and select 2000 HVGs.
# ----------------------------------------------------------------------

N_HVG = 2000


def normalise_and_log(adata):
    """In-place CP10k normalisation followed by log1p."""
    # TODO:
    #  - sc.pp.normalize_total(adata, target_sum=1e4)
    #  - sc.pp.log1p(adata)
    raise NotImplementedError


def select_hvgs(adata, n_top: int = N_HVG):
    """Annotate adata.var['highly_variable']; subset to HVGs."""
    # TODO:
    #  - sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat")
    #  - keep a *copy* of the un-subset adata for marker tests later
    #  - return adata[:, adata.var['highly_variable']].copy()
    raise NotImplementedError


# adata_full = adata.copy()
# normalise_and_log(adata)
# adata_hvg  = select_hvgs(adata)
'''

STEP3_SOLUTION = '''# Reference solution — Step 3.

N_HVG = 2000


def normalise_and_log(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)


def select_hvgs(adata, n_top: int = N_HVG):
    # 'seurat' flavor works on log-normalised data (which we just produced).
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top, flavor="seurat",
    )
    sub = adata[:, adata.var["highly_variable"]].copy()
    return sub


# Keep the full log-normalised matrix around so we can run marker tests on
# all genes (Wilcoxon at Step 5 wants more than just HVGs).
normalise_and_log(adata)
adata_full = adata.copy()
adata_hvg  = select_hvgs(adata)

n_hvg = int(adata.var["highly_variable"].sum())
print(f"  HVGs flagged: {n_hvg} / {adata.shape[1]} genes")
print(f"  HVG matrix shape: {adata_hvg.shape}")

# Visualise the dispersion vs mean plot — HVGs sit above the trend line.
means = np.asarray(adata.var["means"])
disps = np.asarray(adata.var["dispersions_norm"])
hvg   = np.asarray(adata.var["highly_variable"])
fig, ax = plt.subplots(figsize=(6, 4))
ax.scatter(means[~hvg], disps[~hvg], s=4, alpha=0.4, c="grey", label="other")
ax.scatter(means[hvg],  disps[hvg],  s=4, alpha=0.8, c="C3",   label="HVG")
ax.set_xscale("log")
ax.set_xlabel("mean expression (CP10k log1p)")
ax.set_ylabel("normalised dispersion")
ax.set_title(f"Highly variable gene selection (top {N_HVG})")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP4_MD = """## Step 4 (18 min) — PCA → UMAP → Leiden

The PCA on HVG-restricted cells gives us the best linear rank-30
approximation of the data. Those 30 PCs are what UMAP and the kNN graph
both consume:

- **kNN graph** on PC-space (k=15 neighbours) — the local structure of
  the manifold.
- **UMAP** lays the manifold out in 2-D for visualisation.
- **Leiden** finds graph communities that maximise modularity; resolution
  0.5 yields a reasonable cluster count for PBMC3k (~7-10 clusters
  corresponding to T / B / NK / monocyte / DC subsets).
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — PCA, kNN graph, UMAP, Leiden clustering.
# ----------------------------------------------------------------------

N_PCS   = 30
N_NEIGH = 15
LEIDEN_RES = 0.5


def reduce_and_cluster(adata):
    """Run PCA, build kNN graph, embed with UMAP, run Leiden."""
    # TODO (in order):
    #  - sc.pp.scale(adata, max_value=10)
    #  - sc.tl.pca(adata, n_comps=N_PCS)
    #  - sc.pp.neighbors(adata, n_neighbors=N_NEIGH, n_pcs=N_PCS)
    #  - sc.tl.umap(adata)
    #  - sc.tl.leiden(adata, resolution=LEIDEN_RES, random_state=SEED)
    raise NotImplementedError


def plot_umap(adata):
    """2-D UMAP coloured by Leiden cluster."""
    # TODO: scatter adata.obsm['X_umap'] coloured by adata.obs['leiden'].
    raise NotImplementedError


# reduce_and_cluster(adata_hvg)
# plot_umap(adata_hvg)
'''

STEP4_SOLUTION = '''# Reference solution — Step 4.

N_PCS      = 30
N_NEIGH    = 15
LEIDEN_RES = 0.5


def reduce_and_cluster(adata):
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=N_PCS, random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGH, n_pcs=N_PCS, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, random_state=SEED)


def plot_umap(adata):
    coords = adata.obsm["X_umap"]
    labels = adata.obs["leiden"].astype(str).values
    cats = sorted(set(labels), key=lambda s: int(s) if s.isdigit() else s)
    cmap = plt.cm.tab20(np.linspace(0, 1, max(len(cats), 1)))

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for i, c in enumerate(cats):
        mask = labels == c
        ax.scatter(coords[mask, 0], coords[mask, 1], s=4, alpha=0.8,
                   color=cmap[i % len(cmap)], label=c)
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.set_title(f"UMAP coloured by Leiden cluster (res={LEIDEN_RES})")
    ax.legend(title="cluster", bbox_to_anchor=(1.02, 1), loc="upper left",
              fontsize=8, markerscale=2)
    plt.tight_layout()
    plt.show()


reduce_and_cluster(adata_hvg)
plot_umap(adata_hvg)
n_clusters = adata_hvg.obs["leiden"].nunique()
print(f"  Leiden produced {n_clusters} clusters at res={LEIDEN_RES}")

# Variance explained by the first 30 PCs — a quick sanity-check on the
# dimensionality-reduction quality.
var_ratio = adata_hvg.uns["pca"]["variance_ratio"]
print(f"  total variance explained by {N_PCS} PCs: {var_ratio.sum():.1%}")
'''


STEP5_MD = """## Step 5 (9 min) — Wilcoxon marker test + heatmap

For one cluster we run a Wilcoxon rank-sum test on each gene (cluster vs
all-other-cells) and keep the top 5 markers by score. We then visualise
their expression across **all** clusters as a heatmap (mean log-CP10k per
cluster, z-scored across rows). Real markers light up *only* in the home
cluster — that is the visual signature of a clean cluster.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Wilcoxon marker test for one cluster, top-5 heatmap.
# ----------------------------------------------------------------------

TARGET_CLUSTER = "0"   # we look at the largest cluster by default
N_TOP_MARKERS  = 5


def find_markers(adata_hvg, adata_full, target: str, n_top: int):
    """Return (marker_genes, df) where df is a tidy long-form table."""
    # TODO:
    #  - Carry over Leiden labels from adata_hvg to adata_full.
    #  - sc.tl.rank_genes_groups(adata_full, "leiden", method="wilcoxon")
    #  - Pull the top-n gene names for the target group.
    raise NotImplementedError


def plot_marker_heatmap(adata_full, markers):
    """Heatmap of marker mean log-expression per Leiden cluster (z-scored)."""
    # TODO:
    #  - Loop over clusters; compute mean log-CP10k for each marker.
    #  - Z-score rows; imshow with seaborn-style colormap.
    raise NotImplementedError


# markers, df = find_markers(adata_hvg, adata_full, TARGET_CLUSTER, N_TOP_MARKERS)
# plot_marker_heatmap(adata_full, markers)
'''

STEP5_SOLUTION = '''# Reference solution — Step 5.

TARGET_CLUSTER = sorted(adata_hvg.obs["leiden"].astype(str).unique(),
                         key=lambda s: int(s) if s.isdigit() else s)[0]
N_TOP_MARKERS  = 5


def find_markers(adata_hvg, adata_full, target: str, n_top: int):
    # Propagate Leiden labels onto the full (all-gene) matrix.
    adata_full.obs["leiden"] = adata_hvg.obs["leiden"].astype("category").values

    sc.tl.rank_genes_groups(
        adata_full, groupby="leiden", method="wilcoxon", n_genes=50,
    )
    names = pd.DataFrame(adata_full.uns["rank_genes_groups"]["names"])
    scores = pd.DataFrame(adata_full.uns["rank_genes_groups"]["scores"])
    top = names[target].iloc[:n_top].tolist()
    return top, names


def plot_marker_heatmap(adata_full, markers):
    clusters = sorted(adata_full.obs["leiden"].astype(str).unique(),
                       key=lambda s: int(s) if s.isdigit() else s)
    mat = np.zeros((len(markers), len(clusters)))
    X = adata_full.X
    is_sparse = sp.issparse(X)

    gene_to_col = {g: i for i, g in enumerate(adata_full.var_names)}

    for j, c in enumerate(clusters):
        mask = (adata_full.obs["leiden"].astype(str) == c).values
        for i, gene in enumerate(markers):
            col_idx = gene_to_col.get(gene)
            if col_idx is None:
                continue
            col = X[:, col_idx]
            if is_sparse:
                col = np.asarray(col.todense()).ravel()
            else:
                col = np.asarray(col).ravel()
            mat[i, j] = col[mask].mean()

    # Z-score per row.
    row_mu  = mat.mean(axis=1, keepdims=True)
    row_sd  = mat.std(axis=1, keepdims=True) + 1e-9
    z = (mat - row_mu) / row_sd

    fig, ax = plt.subplots(figsize=(max(5, len(clusters) * 0.8),
                                     max(2.5, len(markers) * 0.5)))
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels(clusters)
    ax.set_yticks(range(len(markers)))
    ax.set_yticklabels(markers)
    ax.set_xlabel("Leiden cluster")
    ax.set_ylabel("marker gene")
    ax.set_title(f"Top-{len(markers)} markers for cluster {TARGET_CLUSTER} "
                  f"(row-zscored mean log-CP10k)")
    fig.colorbar(im, ax=ax, shrink=0.7, label="z-score")
    plt.tight_layout()
    plt.show()


markers, all_names = find_markers(adata_hvg, adata_full, TARGET_CLUSTER, N_TOP_MARKERS)
print(f"Top-{N_TOP_MARKERS} markers for cluster {TARGET_CLUSTER}: {markers}")
plot_marker_heatmap(adata_full, markers)
'''


SELFCHECK_MD = """## Self-check

These asserts validate that each pipeline stage produced sensible numerical
output. If you ran the reference solutions above they should all pass; if
you wrote your own and an assert fails, revisit the corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Cell count is in the expected band: pbmc3k drops from ~2700 to ~2500
#    after our QC; synthetic fallback starts at 600 and ends close to it.
n_cells = adata_hvg.shape[0]
assert 100 < n_cells < 4000, f"unexpected cell count after QC: {n_cells}"

# 2. We asked for 2000 HVGs.
n_hvg = adata_hvg.shape[1]
assert 1000 <= n_hvg <= 2500, f"HVG count out of band: {n_hvg}"

# 3. PCA produced exactly N_PCS components and a reasonable variance share.
assert "X_pca" in adata_hvg.obsm, "PCA output missing"
assert adata_hvg.obsm["X_pca"].shape[1] == N_PCS
total_var = adata_hvg.uns["pca"]["variance_ratio"].sum()
# PBMC3k legitimately reaches ~10% in 30 PCs (highly sparse, ~2k HVGs);
# synthetic fallback concentrates more variance, so the band is wide.
assert total_var > 0.05, f"PC variance ratio suspiciously low: {total_var:.3f}"

# 4. UMAP coordinates exist and have finite range.
assert "X_umap" in adata_hvg.obsm
umap_range = adata_hvg.obsm["X_umap"].ptp(axis=0)
assert np.all(umap_range > 1), f"UMAP collapsed: range={umap_range}"

# 5. Leiden produced more than one cluster.
n_clusters = adata_hvg.obs["leiden"].nunique()
assert n_clusters >= 2, f"Leiden produced only {n_clusters} cluster(s)"

# 6. Markers came back as strings present in the var_names index.
assert len(markers) == N_TOP_MARKERS, f"expected {N_TOP_MARKERS} markers, got {len(markers)}"
for g in markers:
    assert g in adata_full.var_names, f"marker {g!r} not in var_names"

print(f"Pipeline summary:")
print(f"  cells (post-QC):   {n_cells}")
print(f"  HVGs:              {n_hvg}")
print(f"  Leiden clusters:   {n_clusters}")
print(f"  top-cluster markers: {markers}")
print()
print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — manifold learning + graph clustering on a sparse signal

The five steps you just ran are a linear-algebra stack with one nonlinear
embedding and one graph-clustering pass on top:

1. **Normalisation = sensor calibration.** CP10k removes the per-cell
   sequencing-depth gain; log1p reshapes the heavy-tailed multiplicative
   noise into roughly additive noise. The same idea as taking the dB of
   a power measurement before any further analysis.
2. **HVG selection = feature gating.** Genes whose dispersion sits on
   the trend line (housekeeping, noise) carry no discriminative
   information — keeping them adds rank but no signal-to-noise. Picking
   the top 2 000 by normalised dispersion is just a SNR-ordered feature
   selection step.
3. **PCA = optimal rank-30 linear approximation.** Under a Gaussian noise
   model the top-K PCs are the **minimum-variance unbiased estimator**
   of the K-dimensional signal subspace. The PC variance ratio you
   printed is the explained-variance fraction — the same metric a comms
   engineer uses to set a Karhunen-Loeve transform.
4. **kNN + UMAP = nonlinear manifold embedding.** UMAP optimises a
   cross-entropy between local kNN distances and an embedded 2-D layout.
   It is a nonlinear cousin of multidimensional scaling, and a special
   case of stochastic-neighbour embedding (t-SNE is its sibling).
5. **Leiden = modularity-maximising graph partition.** Modularity is the
   sum over edges of `A_ij − k_i k_j / 2m`, the deviation of observed
   from expected (random-graph) edge density. Maximising it is **spectral
   clustering in disguise**: the top eigenvectors of the modularity
   matrix span the same low-frequency subspace as the bottom eigenvectors
   of the graph Laplacian.

The whole pipeline is a textbook signal-processing chain — sensor
calibration → feature selection → dimensionality reduction → manifold
embedding → clustering — with one biological twist: the manifold has
**discrete clumps** because cell types are quasi-discrete attractors of
the gene-regulatory dynamical system. That clumpiness is what makes the
pipeline work.
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
        1 for c in nb2.cells
        if c.cell_type == "code"
        and c.metadata.get("jupyter", {}).get("source_hidden")
    )
    print(f"Wrote {out_path}")
    print(f"  cells: {len(nb2.cells)} total  ({n_md} md, {n_code} code, {n_hidden} hidden-solution)")
