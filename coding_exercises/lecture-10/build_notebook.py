"""Build exercise.ipynb for L10 — Methylation, Hi-C, and 3D Genome Organisation.

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


TITLE_MD = """# L10 — Methylation, Hi-C, and 3D Genome Organisation

In this exercise you treat two epigenomic data modalities as **linear-algebra
problems**:

1. **Methylation** as a *proportion* with quantified uncertainty — you'll fit
   the conjugate **Beta posterior** to bisulfite count data and watch how
   sampling depth tightens the credible interval.
2. **Hi-C** as a *covariance-like contact matrix* — you'll **ICE-normalise**
   away coverage bias, score **TAD boundaries** as change-points on a
   1-D insulation signal, and recover **A/B compartments** as the sign of
   the first eigenvector after observed-over-expected correction.

Everything is generated synthetically inside the notebook from a fixed seed,
so it runs in under 5 minutes on free Colab CPU with no external fetches.
"""


AHA_MD = """> **Aha.** Methylation proportions are **Beta-distributed**: a conjugate
> prior + binomial likelihood gives you a closed-form posterior, which is the
> right way to handle low-coverage CpGs. Hi-C contact matrices are a
> **covariance** of physical contact frequencies; **TADs** are blocks on the
> diagonal that you find with an **edge detector** (insulation score = 1-D
> change-point detection), and **A/B compartments** are the **sign of the
> first principal component** of the contact correlation matrix.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 matplotlib==3.8.4 -q
"""


IMPORTS = """import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats
from scipy.signal import find_peaks

# Deterministic for the whole notebook.
RNG = np.random.default_rng(42)
np.random.seed(42)

print(f"numpy {np.__version__}, pandas {pd.__version__}")
"""


STEP1_MD = """## Step 1 (8 min) — Generate the methylation landscape and spot a CpG island

We simulate a 100-CpG × 5-sample bisulfite count table. The first ~80 CpGs sit
in **open sea** (global background) and are **mostly methylated** (β ≈ 0.85),
while CpGs 30–50 sit in a **CpG island** that is **mostly unmethylated**
(β ≈ 0.05). Each CpG-sample pair stores `(M, U)` = (methylated reads,
unmethylated reads); the observed proportion is `β̂ = M / (M + U)`.

A CpG **island** shows up as a dark vertical stripe across all samples in the
β̂ heatmap — that's the regulatory signature of a typically-unmethylated
promoter region.
"""

STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Simulate methylation counts and plot the landscape.
# ----------------------------------------------------------------------

N_CPG     = 100   # CpGs along the locus
N_SAMPLE  = 5     # biological replicates
ISLAND    = (30, 50)  # half-open: CpGs 30..49 inclusive form the island
MEAN_COV  = 20    # mean read depth per (cpg, sample)


def simulate_methylation(n_cpg=N_CPG, n_sample=N_SAMPLE, island=ISLAND,
                         mean_cov=MEAN_COV, seed=0):
    """Return three (n_cpg, n_sample) arrays: M, U, beta_hat = M/(M+U)."""
    # TODO:
    # 1. Pick a per-CpG true methylation proportion (0.85 outside the island,
    #    0.05 inside the island). Add a small per-CpG noise (~Beta(80,20) etc.)
    #    so each CpG isn't identical.
    # 2. Sample a Poisson depth per (cpg, sample) with mean = mean_cov.
    # 3. M ~ Binomial(depth, true_beta); U = depth - M.
    # 4. beta_hat = M / max(1, M+U).
    raise NotImplementedError


# Run it; produce M, U, beta_hat with shape (N_CPG, N_SAMPLE).
# Visualise beta_hat as a heatmap. The island should be a dark vertical band.
'''


STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""


STEP1_SOLUTION = '''# Reference solution — Step 1.

N_CPG     = 100
N_SAMPLE  = 5
ISLAND    = (30, 50)
MEAN_COV  = 20


def simulate_methylation(n_cpg=N_CPG, n_sample=N_SAMPLE, island=ISLAND,
                         mean_cov=MEAN_COV, seed=0):
    rng = np.random.default_rng(seed)
    # Per-CpG true beta: 0.85 background, 0.05 inside island.
    true_beta = np.full(n_cpg, 0.85)
    true_beta[island[0]:island[1]] = 0.05
    # A little per-CpG dispersion around its centre so the heatmap is not flat.
    true_beta = np.clip(true_beta + rng.normal(0, 0.05, size=n_cpg), 0.01, 0.99)

    depth = rng.poisson(lam=mean_cov, size=(n_cpg, n_sample))
    depth = np.maximum(depth, 1)  # avoid divide-by-zero
    M = rng.binomial(n=depth, p=true_beta[:, None])
    U = depth - M
    beta_hat = M / (M + U)
    return M, U, beta_hat, true_beta


M, U, beta_hat, true_beta = simulate_methylation(seed=42)
print(f"M shape: {M.shape}, mean depth: {(M+U).mean():.1f}x")
print(f"Background mean beta_hat (outside island): {beta_hat[np.r_[:ISLAND[0], ISLAND[1]:N_CPG]].mean():.3f}")
print(f"Island     mean beta_hat (CpGs 30..49):    {beta_hat[ISLAND[0]:ISLAND[1]].mean():.3f}")

fig, ax = plt.subplots(figsize=(9, 3))
im = ax.imshow(beta_hat.T, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
ax.set_xlabel("CpG index")
ax.set_ylabel("sample")
ax.set_title(r"observed methylation $\\hat\\beta = M/(M+U)$  —  dark band = CpG island")
ax.axvspan(ISLAND[0] - 0.5, ISLAND[1] - 0.5, fill=False, edgecolor="black", lw=1.5, linestyle="--")
plt.colorbar(im, ax=ax, label=r"$\\hat\\beta$")
plt.tight_layout()
plt.show()
'''


STEP2_MD = """## Step 2 (10 min) — Beta posterior for one sample

A methylation proportion estimated from `(M, U)` reads is a **Binomial
proportion** estimate. The conjugate prior is `Beta(α, β)`. With a flat
prior `Beta(1, 1)` (Laplace's rule of succession), the posterior is

$$\\text{Beta}(1 + M,\\; 1 + U)$$

— a closed form, no MCMC needed. The 95 % credible interval is the
`[0.025, 0.975]` quantile range of that Beta. Crucially, the **width**
shrinks as `1/√(M+U)`: a CpG covered by 5 reads has a much wider CI than
one covered by 50.

You will plot the **per-CpG 95 % CI** for one sample and verify the
shrinkage-with-depth scaling.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Beta posterior + 95% credible intervals.
# ----------------------------------------------------------------------

def beta_posterior(m_arr, u_arr, alpha_prior=1.0, beta_prior=1.0):
    """Return posterior alpha and beta arrays (shape = m_arr.shape)."""
    # TODO: alpha_post = alpha_prior + m_arr ; beta_post = beta_prior + u_arr
    raise NotImplementedError


def credible_interval(alpha_post, beta_post, level=0.95):
    """Return (lo, hi) arrays using the Beta distribution quantile function."""
    # TODO: use scipy.stats.beta.ppf at (1-level)/2 and 1-(1-level)/2.
    raise NotImplementedError


# Pick sample 0; compute alpha_post, beta_post, lo, hi (each length N_CPG).
# Plot mean (alpha/(alpha+beta)) with error bars = [mean-lo, hi-mean].
# Then verify CI-width ~ 1/sqrt(depth) by scatter.
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.

def beta_posterior(m_arr, u_arr, alpha_prior=1.0, beta_prior=1.0):
    return alpha_prior + m_arr, beta_prior + u_arr


def credible_interval(alpha_post, beta_post, level=0.95):
    lo_q = (1.0 - level) / 2.0
    hi_q = 1.0 - lo_q
    return stats.beta.ppf(lo_q, alpha_post, beta_post), stats.beta.ppf(hi_q, alpha_post, beta_post)


sample = 0
alpha_post, beta_post = beta_posterior(M[:, sample], U[:, sample])
mean_post = alpha_post / (alpha_post + beta_post)
lo, hi = credible_interval(alpha_post, beta_post, level=0.95)

fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))

axes[0].errorbar(np.arange(N_CPG), mean_post,
                 yerr=[mean_post - lo, hi - mean_post],
                 fmt="o", ms=3, lw=0.8, capsize=0, alpha=0.7)
axes[0].axvspan(ISLAND[0] - 0.5, ISLAND[1] - 0.5, color="orange", alpha=0.15,
                label="CpG island")
axes[0].set_xlabel("CpG index")
axes[0].set_ylabel(r"posterior mean $\\beta$ (95% CI)")
axes[0].set_title(f"Beta posterior — sample {sample}")
axes[0].set_ylim(-0.02, 1.02)
axes[0].legend(loc="upper right")

# CI width vs depth: should scale as ~ 1/sqrt(depth).
depth = M[:, sample] + U[:, sample]
width = hi - lo
order = np.argsort(depth)
axes[1].scatter(depth, width, s=12, alpha=0.6)
xs = np.linspace(depth.min(), depth.max(), 100)
# Heuristic 1/sqrt scaling: width ~ 2 * z * sqrt(p(1-p)/n) ~ c / sqrt(n).
c = np.median(width * np.sqrt(depth))
axes[1].plot(xs, c / np.sqrt(xs), "r--", lw=2, label=r"$c/\\sqrt{n}$ guide")
axes[1].set_xlabel("read depth (M+U)")
axes[1].set_ylabel("95% CI width")
axes[1].set_title(r"Posterior CI shrinks as $1/\\sqrt{n}$")
axes[1].legend()

plt.tight_layout()
plt.show()

print(f"Median CI width: {np.median(width):.3f}")
print(f"Shrinkage constant c (width * sqrt(depth)): {c:.2f}")
'''


STEP3_MD = """## Step 3 (12 min) — Simulate a Hi-C matrix and ICE-normalise it

A Hi-C contact matrix `C[i, j]` is the count of paired-end reads where one
mate falls in bin `i` and the other in bin `j`. Two structural features
dominate:

1. **Diagonal decay**: bins close in genomic distance contact each other more
   often. `C[i,j] ∝ |i-j|⁻¹` is a common power-law approximation.
2. **TAD blocks**: bins inside the same Topologically Associating Domain
   contact each other ~3-5× more often than expected from distance alone.

On top of that, real data has **per-bin coverage bias** (mappability,
restriction-site density, GC content). **ICE** (Iterative Correction and
Eigenvector decomposition) removes this bias by repeatedly rescaling rows and
columns to a common marginal sum. After ICE, every row/column sums to (close
to) the same value.

We simulate a 50×50 matrix for a 2 Mb region binned at 40 kb, inject a known
per-bin bias vector, and watch ICE recover the unbiased matrix.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Simulate Hi-C with TADs + bias, then ICE-normalise.
# ----------------------------------------------------------------------

N_BIN = 50
TAD_BOUNDARIES = [0, 12, 25, 38, 50]  # 4 TADs: [0,12), [12,25), [25,38), [38,50)


def simulate_hic(n_bin=N_BIN, tad_boundaries=TAD_BOUNDARIES, seed=0):
    """Return (C_raw, bias_true) where C_raw has per-bin bias applied."""
    # TODO:
    # 1. Build the *unbiased* expected matrix:
    #       E[i,j] = base_rate / (1 + |i-j|)   (diagonal decay)
    #       + tad_boost if i and j sit in the same TAD
    # 2. Sample C_clean[i,j] ~ Poisson(E[i,j]). Symmetrise.
    # 3. Draw bias_true ~ LogNormal(0, 0.4); inject as C_raw[i,j] = b[i]*b[j]*C_clean[i,j].
    raise NotImplementedError


def ice_normalize(C, n_iter=30, tol=1e-5):
    """Iterative Correction — return C_norm and the per-bin bias estimate.

    Each iteration computes row sums s_i, divides row/column i by sqrt(s_i),
    then re-symmetrises. Convergence when max(|s - mean(s)|) < tol.
    """
    # TODO: standard ICE loop. Mask zero-coverage rows so you do not divide by 0.
    raise NotImplementedError


# C_raw, bias_true = simulate_hic(seed=7)
# C_ice, bias_est = ice_normalize(C_raw)
# Plot raw vs ICE side-by-side (log scale) and check that ICE row sums are flat.
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.

N_BIN = 50
TAD_BOUNDARIES = [0, 12, 25, 38, 50]


def _tad_of(i, bounds):
    for t, (a, b) in enumerate(zip(bounds[:-1], bounds[1:])):
        if a <= i < b:
            return t
    return -1


def simulate_hic(n_bin=N_BIN, tad_boundaries=TAD_BOUNDARIES, seed=0):
    rng = np.random.default_rng(seed)
    tad_id = np.array([_tad_of(i, tad_boundaries) for i in range(n_bin)])

    # Unbiased expectation: diagonal power-law decay + TAD-block boost.
    E = np.zeros((n_bin, n_bin))
    base = 50.0
    for i in range(n_bin):
        for j in range(n_bin):
            d = abs(i - j)
            decay = base / (1.0 + d)               # power-law-ish
            boost = 3.0 if tad_id[i] == tad_id[j] else 1.0
            E[i, j] = decay * boost

    # Sample symmetric counts.
    C_clean = rng.poisson(E)
    C_clean = (C_clean + C_clean.T) // 2  # symmetrise

    # Per-bin bias.
    bias_true = rng.lognormal(mean=0.0, sigma=0.4, size=n_bin)
    bias_true = bias_true / np.exp(np.log(bias_true).mean())  # geom-mean 1
    C_raw = bias_true[:, None] * bias_true[None, :] * C_clean
    return C_raw.astype(float), bias_true, C_clean.astype(float)


def ice_normalize(C, n_iter=30, tol=1e-5):
    C = C.astype(float).copy()
    n = C.shape[0]
    bias_est = np.ones(n)

    # Mask zero-coverage rows.
    mask = C.sum(axis=1) > 0
    for it in range(n_iter):
        s = C.sum(axis=1)
        # Target: every nonzero row sums to the mean of nonzero row sums.
        target = s[mask].mean()
        scale = np.where(mask & (s > 0), s / target, 1.0)
        sqrt_scale = np.sqrt(scale)
        # Avoid division by zero
        sqrt_scale[sqrt_scale == 0] = 1.0
        C = C / (sqrt_scale[:, None] * sqrt_scale[None, :])
        bias_est *= sqrt_scale
        # Check convergence.
        s_new = C.sum(axis=1)
        if mask.any() and np.max(np.abs(s_new[mask] - s_new[mask].mean())) < tol * target:
            break
    return C, bias_est


C_raw, bias_true, C_clean = simulate_hic(seed=7)
C_ice, bias_est = ice_normalize(C_raw)

fig, axes = plt.subplots(1, 3, figsize=(13, 4))
for ax, mat, title in zip(
    axes,
    [C_clean, C_raw, C_ice],
    ["unbiased truth", "raw (biased)", "after ICE"],
):
    im = ax.imshow(np.log1p(mat), cmap="hot", aspect="equal")
    ax.set_title(title)
    ax.set_xlabel("bin")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
axes[0].set_ylabel("bin")
plt.tight_layout()
plt.show()

# Row sums should be approximately flat after ICE.
print(f"Row-sum SD before ICE: {C_raw.sum(axis=1).std():.2f}")
print(f"Row-sum SD after  ICE: {C_ice.sum(axis=1).std():.2f}")

# Correlate recovered bias to truth.
r_bias = np.corrcoef(np.log(bias_true), np.log(bias_est))[0, 1]
print(f"Pearson(log bias_true, log bias_est) = {r_bias:.3f}")
'''


STEP4_MD = """## Step 4 (18 min) — Insulation score and TAD-boundary detection

Once the matrix is ICE-normalised, a **TAD boundary** is a position where
upstream-vs-downstream contacts drop sharply. The standard summary statistic
is the **insulation score** at bin `i`:

$$I_i = \\text{mean}\\bigl( C[i-w:i, \\; i:i+w] \\bigr)$$

— the average contact frequency in the off-diagonal `w × w` square centred
on `i`. Low `I_i` ⇒ few contacts crossing position `i` ⇒ boundary. Local
**minima** of the insulation score mark TAD boundaries; this is equivalent
to a **1-D edge detector** convolved with the diagonal of the matrix.

You'll compute `I` with a window of 5 bins (200 kb at 40 kb resolution),
plot it alongside the matrix, and recover the four planted boundary
positions with `scipy.signal.find_peaks` applied to `-I`.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Insulation score → TAD boundary calls.
# ----------------------------------------------------------------------

WINDOW = 5  # bins; 5 * 40 kb = 200 kb


def insulation_score(C, window=WINDOW):
    """Return I, an array of length n. At each i, average C[i-w:i, i:i+w]."""
    # TODO: handle edges by setting I[i] = NaN where the square would run off
    # the matrix. Otherwise just mean over the off-diagonal block.
    raise NotImplementedError


def call_boundaries(I, min_dist=4):
    """Return integer indices of local minima of I (= TAD boundary calls)."""
    # TODO: scipy.signal.find_peaks on -I with distance=min_dist.
    raise NotImplementedError


# Compute I, plot below the C_ice heatmap, mark called boundaries vs truth.
# Compare set(called) to set(TAD_BOUNDARIES[1:-1]).
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.

WINDOW = 5  # bins


def insulation_score(C, window=WINDOW):
    n = C.shape[0]
    I = np.full(n, np.nan)
    for i in range(window, n - window):
        block = C[i - window:i, i:i + window]
        I[i] = block.mean()
    return I


def call_boundaries(I, min_dist=4):
    # Mask NaNs by replacing with the max so they aren't picked as minima.
    I_clean = np.where(np.isnan(I), np.nanmax(I) + 1, I)
    peaks, _ = find_peaks(-I_clean, distance=min_dist)
    return peaks


I = insulation_score(C_ice, window=WINDOW)
called = call_boundaries(I)
true_inner = TAD_BOUNDARIES[1:-1]  # interior boundaries only

fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                          gridspec_kw={"height_ratios": [3, 1]})

im = axes[0].imshow(np.log1p(C_ice), cmap="hot", aspect="auto",
                    extent=[0, N_BIN, N_BIN, 0])
for b in true_inner:
    axes[0].axvline(b, color="cyan", lw=1.2, alpha=0.8, linestyle=":")
    axes[0].axhline(b, color="cyan", lw=1.2, alpha=0.8, linestyle=":")
axes[0].set_ylabel("bin")
axes[0].set_title("ICE-normalised Hi-C  +  true TAD boundaries (cyan)")
plt.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04)

axes[1].plot(np.arange(N_BIN), I, lw=1.5, color="navy")
for b in true_inner:
    axes[1].axvline(b, color="cyan", lw=1.2, alpha=0.7, linestyle=":", label="truth" if b == true_inner[0] else None)
for c in called:
    axes[1].axvline(c, color="red", lw=1.2, alpha=0.7, linestyle="--", label="called" if c == called[0] else None)
axes[1].set_xlabel("bin")
axes[1].set_ylabel("insulation score")
axes[1].legend(loc="upper right", fontsize=9)

plt.tight_layout()
plt.show()

print(f"True boundaries:   {sorted(true_inner)}")
print(f"Called boundaries: {sorted(called.tolist())}")
matched = sum(1 for b in true_inner if any(abs(b - c) <= 2 for c in called))
print(f"Recovered {matched}/{len(true_inner)} TAD boundaries within ±2 bins.")
'''


STEP5_MD = """## Step 5 (12 min) — A/B compartments via the first eigenvector

A/B compartments are a **megabase-scale** chromatin organisation: the genome
splits into two sets of bins, `A` (active, open) and `B` (inactive, closed),
where bins in the same compartment contact each other more often *even at
long range*. Empirically, you can recover them like this:

1. Compute the **observed/expected (O/E)** matrix: divide each anti-diagonal
   by its mean so distance-decay is normalised away.
2. Compute the Pearson correlation of the O/E matrix (`corr(C_oe)`).
3. Take the **first eigenvector**. **Sign of PC1 ⇒ compartment label**: bins
   with positive PC1 are in one compartment, negative in the other. The
   eigenvalue magnitude reflects how strong the compartmentalisation is.

This is the same operation as **PCA on a covariance matrix** — exactly what
you would do in EE to find the principal direction of a multi-channel
signal. Here the "signal" is the bin-by-bin contact pattern.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Observed/expected, correlation, PC1 = A/B compartments.
# ----------------------------------------------------------------------

# We rebuild a *larger* simulation that has a clear A/B split: every other TAD
# is in compartment A vs B (long-range A-A and B-B contacts boosted).

def simulate_hic_with_compartments(n_bin=N_BIN, seed=0):
    """Return C with both diagonal decay and a long-range A/B boost.

    Compartment label for bin i: A if i // 12 is even, else B. A-A and B-B
    distal contacts get a 1.5x boost; A-B contacts are baseline.
    """
    # TODO: same as simulate_hic but add a compartment_boost when same compartment
    # and |i-j| > some long-range threshold.
    raise NotImplementedError


def observed_over_expected(C):
    """Divide each anti-diagonal of C by its mean."""
    # TODO
    raise NotImplementedError


def compartments_pc1(C):
    """Return PC1 of corrcoef(observed/expected). Sign convention: pin sign so
    PC1[0] >= 0 for reproducibility."""
    # TODO
    raise NotImplementedError


# Compare PC1 sign to the true A/B label; report accuracy.
'''


STEP5_SOLUTION = '''# Reference solution — Step 5.

def simulate_hic_with_compartments(n_bin=N_BIN, seed=0):
    rng = np.random.default_rng(seed)
    # Compartment label: alternate every 12 bins (matches TAD_BOUNDARIES roughly).
    comp = np.array([0 if (i // 12) % 2 == 0 else 1 for i in range(n_bin)])

    base = 50.0
    long_range = 8  # bins
    E = np.zeros((n_bin, n_bin))
    for i in range(n_bin):
        for j in range(n_bin):
            d = abs(i - j)
            decay = base / (1.0 + d)
            boost = 1.0
            # TAD boost (short range, same TAD block).
            if (i // 12) == (j // 12):
                boost *= 2.5
            # A/B compartment boost (long range, same compartment).
            if d >= long_range and comp[i] == comp[j]:
                boost *= 1.6
            E[i, j] = decay * boost
    C_clean = rng.poisson(E)
    C_clean = (C_clean + C_clean.T) // 2
    return C_clean.astype(float), comp


def observed_over_expected(C):
    n = C.shape[0]
    oe = np.zeros_like(C, dtype=float)
    for d in range(n):
        diag_vals = np.array([C[i, i + d] for i in range(n - d)])
        mean_d = diag_vals.mean() if len(diag_vals) > 0 else 1.0
        if mean_d <= 0:
            mean_d = 1.0
        for i in range(n - d):
            oe[i, i + d] = C[i, i + d] / mean_d
            oe[i + d, i] = oe[i, i + d]
    return oe


def compartments_pc1(C):
    oe = observed_over_expected(C)
    # Correlation matrix of the O/E pattern.
    corr = np.corrcoef(oe)
    # Replace NaNs (from zero-variance rows) with 0.
    corr = np.nan_to_num(corr, nan=0.0)
    vals, vecs = np.linalg.eigh(corr)
    # eigh returns ascending; take the largest eigenvalue's eigenvector.
    pc1 = vecs[:, -1]
    # Sign convention.
    if pc1[0] < 0:
        pc1 = -pc1
    return pc1, vals[-1]


C_comp, comp_true = simulate_hic_with_compartments(seed=11)
C_comp_ice, _ = ice_normalize(C_comp)
pc1, lam = compartments_pc1(C_comp_ice)

# Pick the sign that maximises agreement with the truth (compartment label is
# unsigned; the eigenvector has an arbitrary global sign).
pred_a = (pc1 > 0).astype(int)
acc_pos = (pred_a == comp_true).mean()
acc_neg = (1 - pred_a == comp_true).mean()
acc = max(acc_pos, acc_neg)

fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                         gridspec_kw={"height_ratios": [3, 1]})
oe = observed_over_expected(C_comp_ice)
corr = np.corrcoef(oe)
corr = np.nan_to_num(corr, nan=0.0)
im = axes[0].imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
axes[0].set_title("corr(observed/expected) — plaid pattern = A/B compartments")
axes[0].set_ylabel("bin")
plt.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, label="Pearson r")

bar_pc1 = pc1 if acc_pos >= acc_neg else -pc1
colors = ["tab:blue" if v >= 0 else "tab:red" for v in bar_pc1]
axes[1].bar(np.arange(N_BIN), bar_pc1, color=colors, width=1.0, edgecolor="none")
axes[1].axhline(0, color="black", lw=0.5)
axes[1].set_xlabel("bin")
axes[1].set_ylabel("PC1")
axes[1].set_title(f"PC1 sign = compartment label   (accuracy vs truth: {acc:.0%}, λ₁={lam:.2f})")

plt.tight_layout()
plt.show()

print(f"PC1-based A/B compartment recovery accuracy: {acc:.0%}")
print(f"Leading eigenvalue λ₁ = {lam:.2f}  (larger ⇒ stronger compartmentalisation)")
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing pieces of the pipeline. If you ran
the reference solutions above they should all pass; if you wrote your own
and an assert fails, revisit the corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Beta-posterior CI shrinks with depth.
a_lo, b_lo = beta_posterior(np.array([2]), np.array([3]))    # depth 5
a_hi, b_hi = beta_posterior(np.array([20]), np.array([30]))  # depth 50
lo1, hi1 = credible_interval(a_lo, b_lo)
lo2, hi2 = credible_interval(a_hi, b_hi)
w_low_depth = float(np.atleast_1d(hi1 - lo1)[0])
w_high_depth = float(np.atleast_1d(hi2 - lo2)[0])
assert w_high_depth < w_low_depth, (
    f"Higher depth should give a narrower CI; got {w_high_depth:.3f} vs {w_low_depth:.3f}"
)

# 2. CpG-island band: mean beta_hat inside ISLAND is much lower than outside.
inside = beta_hat[ISLAND[0]:ISLAND[1]].mean()
outside = beta_hat[np.r_[:ISLAND[0], ISLAND[1]:N_CPG]].mean()
assert outside - inside > 0.5, (
    f"CpG-island contrast too small: outside={outside:.3f}, inside={inside:.3f}"
)

# 3. ICE flattens row sums.
row_sd_raw = float(C_raw.sum(axis=1).std())
row_sd_ice = float(C_ice.sum(axis=1).std())
assert row_sd_ice < 0.25 * row_sd_raw, (
    f"ICE did not flatten row sums enough: {row_sd_ice:.2f} vs raw {row_sd_raw:.2f}"
)

# 4. TAD-boundary recovery: ≥2 of 3 inner boundaries within ±2 bins.
true_inner = TAD_BOUNDARIES[1:-1]
recovered = sum(1 for b in true_inner if any(abs(b - c) <= 2 for c in called))
assert recovered >= 2, f"Recovered only {recovered}/{len(true_inner)} TAD boundaries"

# 5. A/B compartment recovery accuracy ≥ 0.80.
acc_check = max((pred_a == comp_true).mean(), (1 - pred_a == comp_true).mean())
assert acc_check >= 0.80, f"Compartment recovery too low: {acc_check:.0%}"

print(f"Beta CI width  depth-5  vs depth-50:  {w_low_depth:.3f} -> {w_high_depth:.3f}")
print(f"CpG island contrast (outside - inside): {outside - inside:.3f}")
print(f"ICE row-sum SD reduction:  {row_sd_raw:.2f} -> {row_sd_ice:.2f}")
print(f"TAD boundaries recovered:  {recovered}/{len(true_inner)}")
print(f"A/B compartment accuracy:  {acc_check:.0%}")
print()
print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — Bayesian inference + covariance / PCA / change-point

You used four canonical EE / DSP tools in this exercise:

1. **Beta posterior = conjugate Bayesian inference on a Bernoulli channel.**
   Each bisulfite read is a noisy 1-bit measurement of the CpG's methylation
   state. The Beta posterior is the closed-form Kalman update for a Bernoulli
   observation model — `Beta(α + M, β + U)` is *exactly* the posterior over the
   channel parameter. CI width ∝ 1/√n is the same √n scaling that governs
   Cramér-Rao bounds on any estimator of a binomial proportion.

2. **ICE = iterative bias scaling / proportional fitting.** The same algorithm
   used to balance contingency tables in classical statistics
   (Sinkhorn-Knopp / Iterative Proportional Fitting). It estimates a
   diagonal-rank-1 multiplicative bias `b[i] b[j]` and divides it out. After
   convergence every row marginal equals the geometric mean — a
   *whitening* operation in the row-marginal sense.

3. **Insulation score = 1-D edge detector.** The off-diagonal block mean is
   exactly the cross-correlation of upstream-vs-downstream contact vectors at
   lag 0. A local minimum is the discrete analogue of a derivative
   zero-crossing in a Laplacian-of-Gaussian edge filter. TAD calling = peak
   detection on the negated insulation track.

4. **A/B compartments = PC1 of contact correlation.** The O/E matrix removes
   the deterministic distance trend; the correlation matrix is then the
   sample covariance of the residual contact pattern across bins; PC1 picks
   out the dominant mode of co-variation. This is identical to the spatial
   PCA you would run on an array of correlated sensors to find the
   strongest common-mode signal.

Hi-C teaches us that a contact matrix isn't just an image — it is a
**covariance matrix in disguise**, and every standard covariance-matrix tool
(whitening, eigen-decomposition, change-point detection) maps onto a
biological question (bias correction, compartments, TADs).
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
