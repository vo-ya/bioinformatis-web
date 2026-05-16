"""Build exercise.ipynb for L09 — ChIP-seq, ATAC-seq, and Peak Calling.

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


TITLE_MD = """# L09 — ChIP-seq, ATAC-seq, and Peak Calling

In this exercise you implement the analysis stack for a chromatin-accessibility
experiment: peak calling with a local-Poisson test (MACS2-style CFAR), replicate
concordance, a negative-binomial GLM for differential accessibility (callback
to L06), and a PWM motif scan against the canonical JASPAR MA0139 CTCF matrix.
All data is generated in-notebook from a deterministic seed — no external
fetches.
"""


AHA_MD = """> **Aha.** Peak calling is **constant-false-alarm-rate (CFAR) detection**:
> the threshold is set adaptively from the *local* background so that wide
> regions of elevated noise do not blow up the false-positive rate. Differential
> accessibility reuses the **NB GLM** from L06 — chromatin is just another count
> matrix. A PWM scan is a **matched filter**: the inner product of one-hot DNA
> against a log-odds template peaks at the motif's preferred sites.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 statsmodels==0.14.2 matplotlib==3.8.4 -q
"""


IMPORTS = """import math
import time
from dataclasses import dataclass

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm

# Deterministic for the whole notebook.
RNG = np.random.default_rng(seed=42)
np.random.seed(42)

# 1 Mb synthetic chromosome at 1 bp resolution would be 1M values; we work at
# 10 bp resolution to keep arrays small and computations under 1 s per step.
BIN_SIZE       = 10           # bp per coverage bin
CHROM_LEN_BP   = 1_000_000    # 1 Mb synthetic chromosome
N_BINS         = CHROM_LEN_BP // BIN_SIZE   # 100,000 bins
ALPHABET       = "ACGT"

print(f"Chromosome: {CHROM_LEN_BP/1e6:.1f} Mb at {BIN_SIZE} bp resolution -> {N_BINS:,} bins")
"""


STEP1_MD = """## Step 1 (6 min) — Generate a synthetic ATAC bedgraph and inspect signal vs background

Real ATAC coverage is a low-rate Poisson background with a small fraction of
genomic positions showing sharp pile-ups at open-chromatin sites. We simulate
this directly:

- **Background**: Poisson with mean `lambda_bg = 0.6` reads per 10 bp bin (so
  ~6 reads / 100 bp on average — realistic for 30M-read ATAC).
- **Peaks**: 200 true binding sites at random positions, each modelled as a
  ~250 bp Gaussian pile-up with a per-peak amplitude drawn from a heavy-tail
  log-normal so weak and strong sites coexist.
- **Treated replicate**: the same chromosome with a fraction of peaks
  up-regulated (gain of accessibility) and a different fraction down-regulated.

The aha here is visual: a coverage histogram is a **mixture** of a tight
Poisson lump (background) and a long tail (peaks). The peak-calling task is
to learn the threshold that separates them — adaptively, because the
background rate varies across the genome (mappability, sonication bias,
local depth).
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Generate the synthetic ATAC coverage tracks.
# ----------------------------------------------------------------------

N_PEAKS         = 200
PEAK_WIDTH_BP   = 250         # full-width-ish at half-max
LAMBDA_BG       = 0.6         # background reads per 10 bp bin

# Treatment design: 30% of true peaks get 1.8x stronger, 30% get 0.55x weaker,
# 40% stay the same. This is what differential-binding tests must recover.
FRAC_UP         = 0.30
FRAC_DOWN       = 0.30
UP_FOLD         = 1.80
DOWN_FOLD       = 0.55


def make_peaks(n: int, n_bins: int, rng: np.random.Generator):
    """Return arrays (centres_bin, sigma_bin, amplitudes) for n peaks."""
    # TODO:
    # - sample n peak centres uniformly along the chromosome (in bin units)
    # - convert PEAK_WIDTH_BP to sigma in bin units (sigma_bp ~ PEAK_WIDTH_BP/2.355)
    # - sample amplitudes from a log-normal so a few peaks are very tall
    raise NotImplementedError


def synth_coverage(centres, sigmas, amps, n_bins: int, lam_bg: float,
                   rng: np.random.Generator) -> np.ndarray:
    """Background Poisson + Gaussian peak envelopes, sampled as integer counts."""
    # TODO:
    # - start with the Poisson background array (size n_bins)
    # - for each peak, add a Gaussian-shaped expected pile-up to the rate, then
    #   resample those bins from a Poisson with the local rate
    # Return an int array of read counts per bin.
    raise NotImplementedError


# Generate control + treated tracks.
# centres, sigmas, amps = make_peaks(N_PEAKS, N_BINS, RNG)
# cov_ctrl    = synth_coverage(centres, sigmas, amps, N_BINS, LAMBDA_BG, RNG)
# cov_treated = synth_coverage(centres, sigmas, amps_treated, N_BINS, LAMBDA_BG, RNG)
'''


STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""


STEP1_SOLUTION = '''# Reference solution — Step 1.

N_PEAKS         = 200
PEAK_WIDTH_BP   = 250
LAMBDA_BG       = 0.6
FRAC_UP         = 0.30
FRAC_DOWN       = 0.30
UP_FOLD         = 1.80
DOWN_FOLD       = 0.55


def make_peaks(n: int, n_bins: int, rng: np.random.Generator):
    # Centres in bin units. Avoid the very edges so the Gaussian fits.
    margin = 50
    centres = rng.integers(margin, n_bins - margin, size=n)
    sigma_bp = PEAK_WIDTH_BP / 2.355  # FWHM -> sigma
    sigmas = np.full(n, sigma_bp / BIN_SIZE)
    # log-normal amplitudes: mean ~5 expected pile-up reads at the peak centre,
    # heavy tail up to ~25.
    amps = rng.lognormal(mean=np.log(5.0), sigma=0.6, size=n)
    return centres, sigmas, amps


def synth_coverage(centres, sigmas, amps, n_bins: int, lam_bg: float,
                   rng: np.random.Generator) -> np.ndarray:
    # Local Poisson rate: background + sum of Gaussian peak envelopes.
    rate = np.full(n_bins, lam_bg, dtype=float)
    x = np.arange(n_bins)
    for c, sig, A in zip(centres, sigmas, amps):
        lo = max(0, int(c - 4 * sig))
        hi = min(n_bins, int(c + 4 * sig) + 1)
        rate[lo:hi] += A * np.exp(-0.5 * ((x[lo:hi] - c) / sig) ** 2)
    return rng.poisson(rate).astype(np.int32)


centres, sigmas, amps = make_peaks(N_PEAKS, N_BINS, RNG)

# Assign treatment effect per peak.
peak_class = np.zeros(N_PEAKS, dtype=int)  # 0=null, 1=up, -1=down
idx = RNG.permutation(N_PEAKS)
n_up = int(round(FRAC_UP * N_PEAKS))
n_down = int(round(FRAC_DOWN * N_PEAKS))
peak_class[idx[:n_up]] = 1
peak_class[idx[n_up:n_up + n_down]] = -1

amps_treated = amps.copy()
amps_treated[peak_class == 1]  *= UP_FOLD
amps_treated[peak_class == -1] *= DOWN_FOLD

# Two control replicates + two treated replicates (for the NB GLM in Step 4).
cov_ctrl_a    = synth_coverage(centres, sigmas, amps,         N_BINS, LAMBDA_BG, RNG)
cov_ctrl_b    = synth_coverage(centres, sigmas, amps,         N_BINS, LAMBDA_BG, RNG)
cov_treated_a = synth_coverage(centres, sigmas, amps_treated, N_BINS, LAMBDA_BG, RNG)
cov_treated_b = synth_coverage(centres, sigmas, amps_treated, N_BINS, LAMBDA_BG, RNG)

# Inspect: histogram of bin counts.
fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
axes[0].plot(np.arange(N_BINS) * BIN_SIZE / 1e3, cov_ctrl_a, lw=0.4)
axes[0].set_xlabel("position (kb)")
axes[0].set_ylabel("reads / 10 bp bin")
axes[0].set_title("Control replicate A — coverage track")
axes[0].set_xlim(0, 200)  # zoom into first 200 kb so peaks are visible

bins = np.arange(0, 30)
axes[1].hist(cov_ctrl_a, bins=bins, alpha=0.7, label="observed", color="#4C72B0")
# Overlay pure background Poisson.
bg_pmf = stats.poisson.pmf(bins, mu=LAMBDA_BG) * N_BINS
axes[1].plot(bins + 0.5, bg_pmf, "r-", lw=2, label=f"Poisson({LAMBDA_BG}) [background only]")
axes[1].set_yscale("log")
axes[1].set_xlabel("reads / bin")
axes[1].set_ylabel("count (log)")
axes[1].set_title("Bin-count histogram = background + peak tail")
axes[1].legend()
plt.tight_layout()
plt.show()

print(f"Total reads (control A):  {cov_ctrl_a.sum():,}")
print(f"Total reads (control B):  {cov_ctrl_b.sum():,}")
print(f"Total reads (treated A):  {cov_treated_a.sum():,}")
print(f"Total reads (treated B):  {cov_treated_b.sum():,}")
print(f"True peak classes -> up: {n_up}, down: {n_down}, null: {N_PEAKS - n_up - n_down}")
'''


STEP2_MD = """## Step 2 (12 min) — Local-Poisson peak calling at q ∈ {0.01, 0.05, 0.1}

MACS2's core idea: at every candidate position, compare the observed pile-up to
a **local background** estimate. The local lambda is the maximum of several
window sizes (e.g. 1 kb, 5 kb, 10 kb plus a chromosome-wide floor) — taking the
max means that if any nearby region has elevated coverage (mappability, copy
number, sonication bias) the threshold rises with it. This is the same logic
as a radar's CFAR detector: never let local noise raise the false-alarm rate.

You will:

1. Slide a **window** along the bin track to get smoothed bin counts.
2. At each candidate bin, compute `lambda_local = max(lambda_chrom, lambda_5kb, lambda_10kb)`.
3. p-value = `1 - Poisson(lambda_local).cdf(obs - 1)`.
4. Benjamini-Hochberg FDR; threshold at q ∈ {0.01, 0.05, 0.1}.
5. Merge adjacent significant bins into peaks.
"""


STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Local-Poisson peak calling (CFAR-style).
# ----------------------------------------------------------------------

WIN_BINS         = 25          # ~250 bp window for candidate pile-up estimate
LOCAL_WIN_5KB    = 500         # 5 kb in 10 bp bins
LOCAL_WIN_10KB   = 1000        # 10 kb in 10 bp bins
Q_LEVELS         = [0.01, 0.05, 0.10]


def sliding_sum(x: np.ndarray, w: int) -> np.ndarray:
    """Sum of every contiguous w-length window, output length len(x)-w+1."""
    # TODO: use np.cumsum or np.convolve with np.ones(w).
    raise NotImplementedError


def local_lambda(cov: np.ndarray, window_bins: int) -> float:
    """Chromosome-wide rate scaled to a window of `window_bins` bins."""
    # TODO: total reads / n_bins * window_bins  (a single scalar).
    raise NotImplementedError


def call_peaks_poisson(cov: np.ndarray, q: float) -> pd.DataFrame:
    """Return a DataFrame of called peaks with columns ['start_bin','end_bin','score'].

    score = -log10(q-value at the peak's most significant bin).
    """
    # TODO:
    # 1. Compute the WIN_BINS-wide pile-up at every position.
    # 2. For each position, compute lambda_local = max(chrom, 5kb-window, 10kb-window),
    #    each scaled so they refer to the same WIN_BINS-wide window.
    # 3. p-value at each position = 1 - poisson.cdf(obs - 1, mu=lambda_local).
    # 4. BH-FDR -> q-values; positions with q <= threshold are significant.
    # 5. Merge adjacent significant positions into peaks; gap <= 1 window allowed.
    # Return one row per merged peak.
    raise NotImplementedError


# peaks_q01 = call_peaks_poisson(cov_ctrl_a, 0.01)
# peaks_q05 = call_peaks_poisson(cov_ctrl_a, 0.05)
# peaks_q10 = call_peaks_poisson(cov_ctrl_a, 0.10)
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.

WIN_BINS         = 25          # 25 bins x 10 bp = 250 bp candidate pile-up window
LOCAL_WIN_5KB    = 500         # 5 kb
LOCAL_WIN_10KB   = 1000        # 10 kb
Q_LEVELS         = [0.01, 0.05, 0.10]


def sliding_sum(x: np.ndarray, w: int) -> np.ndarray:
    cs = np.concatenate([[0], np.cumsum(x, dtype=np.int64)])
    return (cs[w:] - cs[:-w]).astype(np.int64)


def local_lambda_scalar(cov: np.ndarray) -> float:
    return cov.mean()


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values for an array of p-values."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    qvals = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p down
    qvals = np.minimum.accumulate(qvals[::-1])[::-1]
    qvals = np.clip(qvals, 0, 1)
    out = np.empty_like(qvals)
    out[order] = qvals
    return out


def call_peaks_poisson(cov: np.ndarray, q: float) -> pd.DataFrame:
    n_bins = len(cov)

    # 1. Pile-up across WIN_BINS-wide windows. obs[i] sums bins [i, i+WIN_BINS).
    obs = sliding_sum(cov, WIN_BINS)

    # 2. Local lambdas — each scaled to a WIN_BINS-wide window so they are comparable.
    #    a) chromosome-wide
    lam_chrom = cov.mean() * WIN_BINS
    #    b) 5 kb window centred on each candidate (clip the sliding sum at the edges)
    sum5 = sliding_sum(cov, LOCAL_WIN_5KB)
    lam5 = sum5 / LOCAL_WIN_5KB * WIN_BINS  # length n_bins - LOCAL_WIN_5KB + 1
    #    c) 10 kb window
    sum10 = sliding_sum(cov, LOCAL_WIN_10KB)
    lam10 = sum10 / LOCAL_WIN_10KB * WIN_BINS

    # Align all three to the same indexing (the WIN_BINS pile-up coordinate).
    n_pos = len(obs)  # n_bins - WIN_BINS + 1
    lam_local = np.full(n_pos, lam_chrom)
    # For the 5kb window: position i refers to bins [i, i+WIN_BINS). Centre the
    # 5kb window on i+WIN_BINS/2. Index into lam5 at i+WIN_BINS//2 - LOCAL_WIN_5KB//2.
    for arr, win in [(lam5, LOCAL_WIN_5KB), (lam10, LOCAL_WIN_10KB)]:
        offset = WIN_BINS // 2 - win // 2
        idx = np.clip(np.arange(n_pos) + offset, 0, len(arr) - 1)
        lam_local = np.maximum(lam_local, arr[idx])

    # 3. Poisson p-values: P(X >= obs | lambda_local) = 1 - cdf(obs-1).
    pvals = stats.poisson.sf(obs - 1, mu=lam_local)
    # 4. BH-FDR q-values.
    qvals = bh_fdr(pvals)

    sig = qvals <= q

    # 5. Merge adjacent significant positions into peaks. Vectorised run-finder:
    #    transitions in `sig` mark peak start/end indices.
    sig_int = sig.astype(np.int8)
    if not sig_int.any():
        return pd.DataFrame(columns=["start_bin", "end_bin", "score"])
    diff = np.diff(np.concatenate([[0], sig_int, [0]]))
    starts = np.where(diff == 1)[0]
    ends   = np.where(diff == -1)[0]  # exclusive

    # -log10 q per significant position; score = max within the run.
    mlq = -np.log10(np.clip(qvals, 1e-300, 1.0))
    scores = np.array([mlq[s:e].max() for s, e in zip(starts, ends)])
    end_bins = ends - 1 + WIN_BINS  # last sig position's window goes up to e-1+WIN_BINS

    return pd.DataFrame({"start_bin": starts, "end_bin": end_bins, "score": scores})


def precompute_qvals(cov: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (qvals, obs) so multiple q-thresholds can be applied without refit."""
    obs = sliding_sum(cov, WIN_BINS)
    lam_chrom = cov.mean() * WIN_BINS
    sum5 = sliding_sum(cov, LOCAL_WIN_5KB)
    lam5 = sum5 / LOCAL_WIN_5KB * WIN_BINS
    sum10 = sliding_sum(cov, LOCAL_WIN_10KB)
    lam10 = sum10 / LOCAL_WIN_10KB * WIN_BINS
    n_pos = len(obs)
    lam_local = np.full(n_pos, lam_chrom)
    for arr, win in [(lam5, LOCAL_WIN_5KB), (lam10, LOCAL_WIN_10KB)]:
        offset = WIN_BINS // 2 - win // 2
        idx = np.clip(np.arange(n_pos) + offset, 0, len(arr) - 1)
        lam_local = np.maximum(lam_local, arr[idx])
    pvals = stats.poisson.sf(obs - 1, mu=lam_local)
    qvals = bh_fdr(pvals)
    return qvals, obs


def merge_significant(qvals: np.ndarray, q: float) -> pd.DataFrame:
    sig = qvals <= q
    sig_int = sig.astype(np.int8)
    if not sig_int.any():
        return pd.DataFrame(columns=["start_bin", "end_bin", "score"])
    diff = np.diff(np.concatenate([[0], sig_int, [0]]))
    starts = np.where(diff == 1)[0]
    ends   = np.where(diff == -1)[0]
    mlq = -np.log10(np.clip(qvals, 1e-300, 1.0))
    scores = np.array([mlq[s:e].max() for s, e in zip(starts, ends)])
    end_bins = ends - 1 + WIN_BINS
    return pd.DataFrame({"start_bin": starts, "end_bin": end_bins, "score": scores})


# Precompute once; threshold at each q level.
t0 = time.time()
qvals_ctrl_a, _ = precompute_qvals(cov_ctrl_a)
print(f"q-value computation (control A): {time.time() - t0:.1f}s")

peaks_by_q = {}
for q in Q_LEVELS:
    peaks_by_q[q] = merge_significant(qvals_ctrl_a, q)
    print(f"q <= {q:.2f}: {len(peaks_by_q[q])} peaks called")

# Sanity vs ground truth: how many called peaks overlap a true peak centre?
def count_recovered(called: pd.DataFrame, true_centres: np.ndarray,
                    slack_bins: int = 30) -> int:
    n_hit = 0
    for c in true_centres:
        mask = (called["start_bin"] - slack_bins <= c) & (c <= called["end_bin"] + slack_bins)
        if mask.any():
            n_hit += 1
    return n_hit

n_true = len(centres)
for q, df in peaks_by_q.items():
    rec = count_recovered(df, centres)
    print(f"  q={q:.2f}: recovered {rec}/{n_true} true peaks (recall={rec/n_true:.0%}); called {len(df)} total")
'''


STEP3_MD = """## Step 3 (12 min) — Replicate concordance

A peak called in one replicate but absent in the other is suspicious. Real
pipelines (IDR, ENCODE) require replicate overlap. We compute it the simple
way: a control-A peak counts as **reproduced** if any control-B peak overlaps
it within a 500 bp window (50 bins).

Build a contingency table and an overlap-recall curve as a function of the
slack window.
"""


STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Replicate concordance between control A and control B.
# ----------------------------------------------------------------------

OVERLAP_SLACK_BP = 500
SLACK_BINS       = OVERLAP_SLACK_BP // BIN_SIZE


def overlap_count(peaks_a: pd.DataFrame, peaks_b: pd.DataFrame,
                  slack_bins: int) -> int:
    """How many rows of peaks_a have at least one peaks_b row within slack_bins?"""
    # TODO: a simple O(N*M) loop is fine at our scale; a sort-and-sweep is nicer.
    raise NotImplementedError


# Re-call peaks on replicate B at q=0.05; report:
#   - peaks called in each replicate
#   - reproduced peaks (slack 500 bp)
#   - sweep the slack from 0 -> 2000 bp and plot recall
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.

OVERLAP_SLACK_BP = 500
SLACK_BINS       = OVERLAP_SLACK_BP // BIN_SIZE


def overlap_count(peaks_a: pd.DataFrame, peaks_b: pd.DataFrame,
                  slack_bins: int) -> int:
    if len(peaks_b) == 0 or len(peaks_a) == 0:
        return 0
    # Sort B by start; for each A peak, binary-search a candidate window.
    b_starts = peaks_b["start_bin"].to_numpy()
    b_ends   = peaks_b["end_bin"].to_numpy()
    order = np.argsort(b_starts)
    b_starts = b_starts[order]
    b_ends   = b_ends[order]

    n_hit = 0
    for s, e in zip(peaks_a["start_bin"], peaks_a["end_bin"]):
        # any b interval where b_start <= e + slack and b_end >= s - slack
        # use searchsorted to skip far-away rows
        i_hi = np.searchsorted(b_starts, e + slack_bins + 1, side="right")
        # look back from i_hi for rows whose end >= s - slack
        if i_hi == 0:
            continue
        sub_ends = b_ends[:i_hi]
        if (sub_ends >= s - slack_bins).any():
            n_hit += 1
    return n_hit


qvals_ctrl_b, _ = precompute_qvals(cov_ctrl_b)
peaks_b_q05 = merge_significant(qvals_ctrl_b, q=0.05)
peaks_a_q05 = peaks_by_q[0.05]

n_a = len(peaks_a_q05)
n_b = len(peaks_b_q05)
n_rep = overlap_count(peaks_a_q05, peaks_b_q05, SLACK_BINS)

print(f"Replicate A peaks (q<=0.05): {n_a}")
print(f"Replicate B peaks (q<=0.05): {n_b}")
print(f"Reproduced (A peak within {OVERLAP_SLACK_BP} bp of any B peak): {n_rep}")
print(f"Reproducibility rate A->B: {n_rep / max(1, n_a):.0%}")

# Sweep slack to show the trade-off.
slacks_bp = np.array([0, 100, 200, 500, 1000, 2000])
recall = []
for s_bp in slacks_bp:
    recall.append(overlap_count(peaks_a_q05, peaks_b_q05, s_bp // BIN_SIZE) / max(1, n_a))

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(slacks_bp, recall, "o-", color="#55A868")
ax.set_xlabel("overlap slack (bp)")
ax.set_ylabel("fraction of A peaks reproduced in B")
ax.set_ylim(0, 1.05)
ax.axhline(0.8, ls="--", color="grey", lw=1, label="ENCODE-style 80% target")
ax.set_title("Replicate-overlap recall vs slack")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP4_MD = """## Step 4 (18 min) — Differential accessibility via an NB GLM

The peak set is fixed (use the replicate-A peaks at q=0.05). For each peak,
count reads in the four samples — two controls + two treated — and fit a
per-peak **negative-binomial** GLM with a single condition coefficient. This is
the same NB GLM as L06 RNA-seq, just on chromatin counts.

We use a fixed dispersion `alpha = 0.02`. (For real ATAC a higher value is
typical, but our synthetic counts have only Poisson-around-mean noise, so a
small dispersion is appropriate. A production pipeline estimates dispersion
per peak with empirical-Bayes shrinkage exactly as in L06.) The log-fold-change
Wald statistic + BH FDR gives a volcano plot.
"""


STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — NB GLM on peak counts, control vs treated.
# ----------------------------------------------------------------------

NB_ALPHA = 0.02  # Fixed dispersion; small because our synthetic data has only Poisson noise around a fixed mean.


def peak_counts(cov: np.ndarray, peaks: pd.DataFrame) -> np.ndarray:
    """Sum cov over each peak's [start_bin, end_bin) range."""
    # TODO
    raise NotImplementedError


def diff_test_nb(counts_matrix: np.ndarray, condition: np.ndarray,
                 alpha: float = NB_ALPHA) -> pd.DataFrame:
    """One NB GLM per peak with a binary condition predictor.

    counts_matrix is (n_peaks, n_samples); condition is (n_samples,) binary.
    Returns a DataFrame with columns: log2FC, lfcSE, wald_z, pvalue, padj.
    Use statsmodels.GLM with NegativeBinomial family and the supplied alpha.
    """
    # TODO
    # For each peak:
    #   y      = counts_matrix[i]                            (n_samples,)
    #   X      = [[1, condition_0], [1, condition_1], ...]   (intercept + condition)
    #   offset = log(library_size) - mean(log(library_size)) (library-size norm)
    # Fit GLM; extract beta[1] (condition), se, z, two-sided p; BH-FDR.
    raise NotImplementedError


# Use the q=0.05 peakset from replicate A as the universe.
# counts_ca = peak_counts(cov_ctrl_a,    peaks_a_q05)
# counts_cb = peak_counts(cov_ctrl_b,    peaks_a_q05)
# counts_ta = peak_counts(cov_treated_a, peaks_a_q05)
# counts_tb = peak_counts(cov_treated_b, peaks_a_q05)
# counts_mat = np.stack([counts_ca, counts_cb, counts_ta, counts_tb], axis=1)
# condition  = np.array([0, 0, 1, 1])
# diff = diff_test_nb(counts_mat, condition)
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.

NB_ALPHA = 0.02


def peak_counts(cov: np.ndarray, peaks: pd.DataFrame) -> np.ndarray:
    out = np.empty(len(peaks), dtype=np.int64)
    for i, (s, e) in enumerate(zip(peaks["start_bin"], peaks["end_bin"])):
        s_clip = max(0, int(s))
        e_clip = min(len(cov), int(e))
        out[i] = int(cov[s_clip:e_clip].sum())
    return out


def diff_test_nb(counts_matrix: np.ndarray, condition: np.ndarray,
                 alpha: float = NB_ALPHA) -> pd.DataFrame:
    """counts_matrix is (n_peaks, n_samples); condition is (n_samples,) binary."""
    n_peaks, n_samples = counts_matrix.shape
    # Library sizes for offset.
    lib = counts_matrix.sum(axis=0).astype(float)
    log_lib = np.log(lib)
    log_lib -= log_lib.mean()

    X = np.column_stack([np.ones(n_samples), condition.astype(float)])

    results = []
    for i in range(n_peaks):
        y = counts_matrix[i].astype(float)
        if y.sum() < n_samples:  # too few reads to fit
            results.append((np.nan, np.nan, np.nan, 1.0))
            continue
        try:
            fam = sm.families.NegativeBinomial(alpha=alpha)
            mod = sm.GLM(y, X, family=fam, offset=log_lib)
            fit = mod.fit(method="bfgs", disp=False)
            beta = fit.params[1]
            se   = fit.bse[1]
            z    = beta / max(se, 1e-12)
            p    = 2.0 * stats.norm.sf(abs(z))
            log2fc = beta / math.log(2.0)
            results.append((log2fc, se / math.log(2.0), z, p))
        except Exception:
            results.append((np.nan, np.nan, np.nan, 1.0))

    df = pd.DataFrame(results, columns=["log2FC", "lfcSE", "wald_z", "pvalue"])
    df["padj"] = bh_fdr(df["pvalue"].fillna(1.0).to_numpy())
    return df


counts_ca = peak_counts(cov_ctrl_a,    peaks_a_q05)
counts_cb = peak_counts(cov_ctrl_b,    peaks_a_q05)
counts_ta = peak_counts(cov_treated_a, peaks_a_q05)
counts_tb = peak_counts(cov_treated_b, peaks_a_q05)
counts_mat = np.stack([counts_ca, counts_cb, counts_ta, counts_tb], axis=1)  # n_peaks x 4
condition  = np.array([0, 0, 1, 1])

diff = diff_test_nb(counts_mat, condition)
diff = pd.concat([peaks_a_q05.reset_index(drop=True), diff], axis=1)

# Volcano plot.
fig, ax = plt.subplots(figsize=(7, 4.5))
sig = (diff["padj"] < 0.05) & (np.abs(diff["log2FC"]) > 0.5)
ax.scatter(diff.loc[~sig, "log2FC"], -np.log10(diff.loc[~sig, "padj"].clip(lower=1e-300)),
           s=10, alpha=0.5, color="grey", label="ns")
ax.scatter(diff.loc[sig & (diff["log2FC"] > 0), "log2FC"],
           -np.log10(diff.loc[sig & (diff["log2FC"] > 0), "padj"].clip(lower=1e-300)),
           s=14, color="#C44E52", label="up in treated")
ax.scatter(diff.loc[sig & (diff["log2FC"] < 0), "log2FC"],
           -np.log10(diff.loc[sig & (diff["log2FC"] < 0), "padj"].clip(lower=1e-300)),
           s=14, color="#4C72B0", label="down in treated")
ax.axhline(-np.log10(0.05), ls="--", color="grey", lw=1)
ax.axvline( 0.5, ls="--", color="grey", lw=1)
ax.axvline(-0.5, ls="--", color="grey", lw=1)
ax.set_xlabel("log2 fold change (treated / control)")
ax.set_ylabel("-log10 adjusted p-value")
ax.set_title("Differential accessibility — NB GLM")
ax.legend()
plt.tight_layout()
plt.show()

n_up_called   = int(((diff["padj"] < 0.05) & (diff["log2FC"] >  0.5)).sum())
n_down_called = int(((diff["padj"] < 0.05) & (diff["log2FC"] < -0.5)).sum())
print(f"Up-regulated peaks called   (padj<0.05, |log2FC|>0.5): {n_up_called}")
print(f"Down-regulated peaks called (padj<0.05, |log2FC|>0.5): {n_down_called}")
print(f"Ground truth: {n_up} true-up + {n_down} true-down peaks "
      f"(plus {N_PEAKS - n_up - n_down} null peaks; you only see the subset that got called).")
'''


STEP5_MD = """## Step 5 (12 min) — PWM matched-filter scan with the canonical CTCF motif

We have peak locations but not their underlying sequence — so we generate a
synthetic chromosome with **planted CTCF motifs** at the true peak centres
(canonical consensus `CCCTCNNNNNGGTGG`) on a 25% uniform-AT/GC background.
Then we scan the differentially-accessible peaks with the JASPAR MA0139.1
CTCF position-weight matrix (counts shown below) and compare hit frequency
inside peak windows vs random background windows.

The PWM scan is a **matched filter**: at every position you compute the inner
product of a one-hot DNA window with the log-odds template. Sites where the
window matches the template's preferred letters score high; the score's
empirical distribution under random DNA gives a calibrated threshold.
"""


STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — PWM scan: canonical CTCF MA0139.1, matched-filter style.
# ----------------------------------------------------------------------

# JASPAR MA0139.1 CTCF position-frequency matrix (rows A, C, G, T; 19 columns).
# Published in JASPAR 2024 from CTCF ChIP-seq peaks, public domain. The consensus
# letter is the per-column argmax and reads as the canonical CTCF binding pattern.
CTCF_COUNTS = np.array([
    # A,  C,  G,  T  ordered as four rows below
    [ 87, 167, 281,  56,   8, 744,  40, 107, 851,   5, 333,  54,  12,  56, 104, 372,  82, 117, 402],  # A
    [291, 145,  49, 800, 903,  13, 528, 433,  11,   0,   3,  12,   0,   8, 733,  13, 482, 322, 181],  # C
    [ 76, 414, 449,  21,   0,  65, 334,  48,  32, 903, 566, 504, 890, 775,   5, 507, 307,  73, 266],  # G
    [459, 187, 134,  36,   2,  91,  11, 324,  18,   3,  17, 343,  11,  74,  72,  21,  42, 401,  64],  # T
])  # shape (4, 19)

PWM_PSEUDO   = 0.5
MOTIF_LEN    = CTCF_COUNTS.shape[1]
# Consensus letters per column (per-column argmax). We plant this exact 19-mer
# at each true peak so the matched-filter scan finds it. This consensus is a
# faithful rendering of MA0139.1's most-frequent base at each position.
ACGT_ORDER   = "ACGT"


def counts_to_log_odds(counts: np.ndarray, pseudo: float = PWM_PSEUDO,
                       bg: np.ndarray | None = None) -> np.ndarray:
    """Convert a 4xL count matrix to a 4xL log2-odds PWM."""
    # TODO:
    # - add pseudo to counts
    # - convert each column to a probability
    # - log2(prob / bg) where bg defaults to uniform 0.25
    raise NotImplementedError


def synth_chromosome(centres_bp: np.ndarray, motif: str,
                     length_bp: int = CHROM_LEN_BP, seed: int = 11) -> str:
    """Random ACGT background with `motif` planted at each centre (truth set)."""
    # TODO:
    # - sample length_bp ACGT bases uniformly
    # - for each centre, overwrite the window [centre - m/2, centre + m/2) with motif
    # - return the joined string
    raise NotImplementedError


def pwm_scan(seq: str, pwm: np.ndarray) -> np.ndarray:
    """Best log-odds score over both strands at every start position."""
    # TODO
    # - one-hot encode the sequence as a 4xN matrix
    # - for each start s in [0, len(seq) - L + 1):
    #     score_fwd = sum over j of pwm[base, j] where base = seq[s+j]
    #     score_rc  = same with reverse-complement of the window
    # - return the per-position max(fwd, rc).
    raise NotImplementedError


# pwm = counts_to_log_odds(CTCF_COUNTS)
# ctcf_consensus = "".join(ACGT_ORDER[i] for i in CTCF_COUNTS.argmax(axis=0))
# chrom_seq = synth_chromosome(centres * BIN_SIZE, motif=ctcf_consensus)
# pwm_scores = pwm_scan(chrom_seq, pwm)
'''


STEP5_SOLUTION = '''# Reference solution — Step 5.

CTCF_COUNTS = np.array([
    [ 87, 167, 281,  56,   8, 744,  40, 107, 851,   5, 333,  54,  12,  56, 104, 372,  82, 117, 402],
    [291, 145,  49, 800, 903,  13, 528, 433,  11,   0,   3,  12,   0,   8, 733,  13, 482, 322, 181],
    [ 76, 414, 449,  21,   0,  65, 334,  48,  32, 903, 566, 504, 890, 775,   5, 507, 307,  73, 266],
    [459, 187, 134,  36,   2,  91,  11, 324,  18,   3,  17, 343,  11,  74,  72,  21,  42, 401,  64],
])

PWM_PSEUDO   = 0.5
MOTIF_LEN    = CTCF_COUNTS.shape[1]
COMPLEMENT   = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
ACGT_ORDER   = "ACGT"


def counts_to_log_odds(counts: np.ndarray, pseudo: float = PWM_PSEUDO,
                       bg: np.ndarray | None = None) -> np.ndarray:
    if bg is None:
        bg = np.full(4, 0.25)
    counts = counts + pseudo
    probs = counts / counts.sum(axis=0, keepdims=True)
    return np.log2(probs / bg[:, None])


def synth_chromosome(centres_bp: np.ndarray, motif: str,
                     length_bp: int = CHROM_LEN_BP, seed: int = 11) -> str:
    """Random ACGT background with `motif` planted at each centre (truth set)."""
    rng = np.random.default_rng(seed)
    base_arr = rng.choice(list(ALPHABET), size=length_bp)
    m = len(motif)
    for c in centres_bp:
        s = int(c) - m // 2
        if 0 <= s and s + m <= length_bp:
            base_arr[s:s + m] = list(motif)
    return "".join(base_arr.tolist())


def _seq_to_index(seq: str) -> np.ndarray:
    idx = np.full(len(seq), -1, dtype=np.int8)
    for k, b in enumerate("ACGT"):
        idx[np.frombuffer(seq.encode("ascii"), dtype=np.uint8) == ord(b)] = k
    return idx


def pwm_scan(seq: str, pwm: np.ndarray) -> np.ndarray:
    L = pwm.shape[1]
    n = len(seq) - L + 1
    idx = _seq_to_index(seq)
    # Forward strand: scores[s] = sum_j pwm[idx[s+j], j]
    # Build a length-n score by accumulating per-column contributions vectorised.
    scores_fwd = np.zeros(n, dtype=np.float32)
    for j in range(L):
        col = pwm[:, j]                       # length 4
        contrib = np.where(idx >= 0, col[idx.clip(min=0)], -100.0)
        scores_fwd += contrib[j:j + n]

    # Reverse-complement strand: PWM reversed left-right, base index 0<->3, 1<->2.
    rc_swap = np.array([3, 2, 1, 0], dtype=np.int8)
    pwm_rc = pwm[rc_swap][:, ::-1]
    scores_rc = np.zeros(n, dtype=np.float32)
    for j in range(L):
        col = pwm_rc[:, j]
        contrib = np.where(idx >= 0, col[idx.clip(min=0)], -100.0)
        scores_rc += contrib[j:j + n]

    return np.maximum(scores_fwd, scores_rc)


pwm = counts_to_log_odds(CTCF_COUNTS)
# Consensus = per-column argmax of the count matrix.
ctcf_consensus = "".join(ACGT_ORDER[i] for i in CTCF_COUNTS.argmax(axis=0))
print(f"PWM shape: {pwm.shape}")
print(f"CTCF MA0139.1 consensus (per-column argmax): {ctcf_consensus}")

# Plant the actual PWM consensus at each peak centre so the scan can find it.
chrom_seq = synth_chromosome(centres * BIN_SIZE, motif=ctcf_consensus)
pwm_scores = pwm_scan(chrom_seq, pwm)

# Compare PWM scores inside diff-accessible peaks vs random background windows.
diff_up_peaks = diff[(diff["padj"] < 0.05) & (diff["log2FC"] > 0.5)]
print(f"Scanning {len(diff_up_peaks)} up-regulated peaks for CTCF...")

def best_pwm_in_peak(peak_row, scores):
    s_bp = int(peak_row["start_bin"]) * BIN_SIZE
    e_bp = int(peak_row["end_bin"])   * BIN_SIZE
    s_idx = max(0, s_bp)
    e_idx = min(len(scores), max(s_idx + 1, e_bp - MOTIF_LEN + 1))
    if e_idx <= s_idx:
        return np.nan
    return float(scores[s_idx:e_idx].max())

peak_best = diff_up_peaks.apply(best_pwm_in_peak, axis=1, scores=pwm_scores).dropna().to_numpy()

# Background: 500 random windows the same width as the median peak.
median_w = int((diff_up_peaks["end_bin"] - diff_up_peaks["start_bin"]).median() * BIN_SIZE)
bg_rng = np.random.default_rng(7)
n_bg = 500
bg_best = []
for _ in range(n_bg):
    s = bg_rng.integers(0, len(pwm_scores) - median_w + 1)
    bg_best.append(float(pwm_scores[s:s + median_w].max()))
bg_best = np.array(bg_best)

threshold = np.quantile(bg_best, 0.95)  # 5% false-alarm rate on background windows
hit_rate_peaks = float((peak_best >= threshold).mean()) if len(peak_best) else 0.0
hit_rate_bg    = float((bg_best   >= threshold).mean())
enrichment     = hit_rate_peaks / max(hit_rate_bg, 1e-6)

fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(bg_best,  bins=30, alpha=0.6, color="grey",   label=f"random windows (n={n_bg})", density=True)
ax.hist(peak_best, bins=30, alpha=0.6, color="#C44E52",
        label=f"up-regulated peaks (n={len(peak_best)})", density=True)
ax.axvline(threshold, ls="--", color="black", lw=1, label=f"5% FA threshold = {threshold:.1f}")
ax.set_xlabel("best CTCF PWM score in window (log2-odds)")
ax.set_ylabel("density")
ax.set_title("CTCF MA0139.1 matched-filter scan")
ax.legend()
plt.tight_layout()
plt.show()

print(f"Hit rate in up-regulated peaks: {hit_rate_peaks:.0%}")
print(f"Hit rate in random background : {hit_rate_bg:.0%}")
print(f"Enrichment factor             : {enrichment:.1f}x")
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces of the pipeline. If
you ran the reference solutions above they should all pass; if you wrote your
own and an assert fails, revisit the corresponding step.
"""


SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Synthetic coverage actually has a heavy upper tail (peaks).
assert cov_ctrl_a.max() > 5 * cov_ctrl_a.mean(), "no peak tail in control coverage"

# 2. Peak calling at q=0.05 recovers a strong majority of true peaks.
peaks_q05_calls = peaks_by_q[0.05]
recovered = count_recovered(peaks_q05_calls, centres)
assert recovered >= 0.7 * N_PEAKS, f"recall too low: {recovered}/{N_PEAKS}"

# 3. q=0.01 should be at least as stringent as q=0.10.
assert len(peaks_by_q[0.01]) <= len(peaks_by_q[0.10]), "FDR levels inverted"

# 4. Replicate concordance is reasonable.
rep_rate = overlap_count(peaks_a_q05, peaks_b_q05, SLACK_BINS) / max(1, len(peaks_a_q05))
assert rep_rate >= 0.6, f"replicate concordance too low: {rep_rate:.0%}"

# 5. NB GLM recovers a majority of the planted up/down peaks at padj<0.05.
n_diff_called = int(((diff["padj"] < 0.05) & (np.abs(diff["log2FC"]) > 0.5)).sum())
true_diff_expected = n_up + n_down  # at most this many of the called peakset
# We only fit on peaks that got called in rep A, so accept >= 30% of expected.
assert n_diff_called >= 0.3 * true_diff_expected, (
    f"diff accessibility recovery suspiciously low: {n_diff_called}/{true_diff_expected}"
)

# 6. PWM hit-rate is enriched at least 2x in peaks vs background.
assert enrichment >= 2.0, f"CTCF PWM enrichment too low: {enrichment:.1f}x"

print("All self-checks passed.")
print(f"  recall @ q=0.05         : {recovered}/{N_PEAKS} = {recovered/N_PEAKS:.0%}")
print(f"  replicate concordance   : {rep_rate:.0%}")
print(f"  diff-acc peaks called   : {n_diff_called} (vs {true_diff_expected} planted)")
print(f"  CTCF enrichment factor  : {enrichment:.1f}x")
'''


EE_MD = """## EE framing — CFAR detection, matched filter, NB GLM as count-domain regression

You implemented the three load-bearing detection problems in chromatin genomics:

1. **Peak calling = CFAR detection.** MACS2's `lambda_local = max(chrom, 5kb, 10kb)`
   is a textbook **constant-false-alarm-rate** detector. The threshold tracks the
   local noise floor; a regional spike in noise pulls the threshold up with it,
   preserving the per-bin false-alarm probability. The same algorithm pattern
   shows up in radar (range cells around the cell under test set the threshold),
   in cosmic-ray spectroscopy, and in network-anomaly detection.

2. **PWM scan = matched filter.** Convolving a log-odds template with the one-hot
   sequence is **literally the matched filter** for a known waveform in additive
   noise. The PWM column entries are log-odds — i.e. **log-likelihood
   contributions** — so the per-position score is a log-likelihood ratio, and
   the empirical null distribution (from random windows) calibrates the false-alarm
   threshold. The MAP detector says: declare a hit at any position where the LLR
   exceeds the threshold.

3. **NB GLM = L06 callback in the count domain.** The same family of models that
   handled RNA-seq differential expression handles differential binding. Counts
   are over-dispersed Poisson → negative binomial; the condition coefficient is
   a log-fold-change estimate; the Wald-test p-value plus BH-FDR is the standard
   multi-test machinery. The point is that **chromatin is just another high-dim
   count matrix** — peaks are the rows, samples are the columns, and the GLM
   stack from L06 transfers verbatim.

Footprinting (the inverse of peak-calling, recovering TF binding sites *inside*
open regions from coverage troughs) is a √N signal-averaging problem: stack
the coverage around N motif instances and the per-bp noise falls as 1/√N, so
the footprint trough emerges from the background. We left that for the
artefact viewers — what you have here is the four-step detection pipeline
that every chromatin paper runs.
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
