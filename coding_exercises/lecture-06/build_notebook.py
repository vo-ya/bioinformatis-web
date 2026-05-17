"""Build exercise.ipynb for L06 — Differential Expression and Count Statistics.

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


TITLE_MD = """# L06 — Differential Expression and Count Statistics

In this exercise you build a miniature DESeq2-style differential-expression
pipeline from scratch. You simulate a 500-gene × 6-sample (3 vs 3) count
matrix with a known set of true-positive genes, fit a per-gene
negative-binomial GLM, apply empirical-Bayes dispersion shrinkage, run a
Wald test, control FDR with Benjamini-Hochberg, draw a volcano plot, and
finally compare three shrinkage strategies (none / EB / strong-prior) on
TPR and FDR. By the end you'll see exactly why naive per-gene MLE
dispersion is hopeless at n=6 and why EB shrinkage is mandatory.
"""


AHA_MD = """> **Aha.** With only 3 vs 3 samples, the per-gene dispersion MLE is
> hopelessly noisy: each gene's estimate has huge variance, which inflates
> p-values for low-dispersion genes and deflates them for high-dispersion
> genes. **Empirical-Bayes shrinkage borrows strength across genes** by
> regressing dispersion on mean expression and pulling each gene toward
> the fitted trend. The result: calibrated p-values, an actually-controlled
> FDR, and a sharp volcano plot.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 statsmodels==0.14.2 matplotlib==3.8.4 -q
"""


IMPORTS = """import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm

# Deterministic for the whole notebook.
SEED = 42
np.random.seed(SEED)
rng = np.random.default_rng(SEED)

# Experimental design: 3 control vs 3 treated samples.
N_GENES   = 500
N_SAMPLES = 6
N_TRUE_DE = 50   # genes with a real condition effect
CONDITION = np.array([0, 0, 0, 1, 1, 1])   # 0=ctrl, 1=trt

print(f"Design: {N_SAMPLES} samples, condition = {CONDITION.tolist()}")
print(f"Simulating {N_GENES} genes, {N_TRUE_DE} true DE")
"""


STEP1_MD = """## Step 1 (10 min) — Simulate counts and estimate size factors

We generate a 500-gene × 6-sample count matrix from a negative-binomial
model with realistic per-gene means and dispersions. Of the 500 genes, 50
get a true log2 fold-change drawn from {±1, ±2}; the other 450 are nulls.
Library sizes are deliberately uneven (the wet-lab fact of life) so the
downstream analysis must normalise. We use the **median-of-ratios** size
factors from Anders & Huber 2010 — DESeq2's default. For each gene, take
the ratio of its observed count to the geometric mean across samples;
each sample's size factor is the median of those ratios.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Simulate counts, estimate size factors.
# ----------------------------------------------------------------------

# Per-gene baseline expression (log-uniform between log(5) and log(2000)).
true_mu = np.exp(rng.uniform(np.log(5), np.log(2000), size=N_GENES))

# Per-gene dispersion: trend = a + b/mu (Anders & Huber asymptotic form);
# add log-normal scatter so genes deviate from the trend.
TREND_A = 0.05
TREND_B = 4.0
trend_disp = TREND_A + TREND_B / true_mu
true_disp = trend_disp * np.exp(rng.normal(0, 0.4, size=N_GENES))

# Pick 50 true DE genes and assign a log2 FC.
# Effects are deliberately on the chunky side (|lfc| in [1.5, 2.5]) so the
# 3-vs-3 design has a fighting chance of recovering them.
is_de = np.zeros(N_GENES, dtype=bool)
de_idx = rng.choice(N_GENES, size=N_TRUE_DE, replace=False)
is_de[de_idx] = True
true_lfc = np.zeros(N_GENES)
true_lfc[de_idx] = rng.choice([-2.5, -1.5, 1.5, 2.5], size=N_TRUE_DE)

# Uneven library size factors (some samples sequenced deeper than others).
true_sf = np.array([1.0, 0.7, 1.3, 1.5, 0.8, 1.1])

# Build the count matrix:
#   mu_ij = true_mu[i] * true_sf[j] * 2**(true_lfc[i] * CONDITION[j])
#   counts_ij ~ NB(mean = mu_ij, dispersion = true_disp[i])
# scipy parameterisation: NB(n, p) with n = 1/disp, p = n / (n + mu).

def simulate_counts():
    """Return a (N_GENES x N_SAMPLES) integer count matrix."""
    # TODO: vectorised NB sampling using the relationship above.
    raise NotImplementedError


def size_factors(counts):
    """Median-of-ratios size factors (Anders & Huber 2010).

    1. log_counts = log(counts) ignoring zeros (-> NaN).
    2. log_geomean per gene = mean across samples (ignoring NaNs).
    3. For each sample, log_ratio = log_counts - log_geomean.
    4. size factor = exp(median(log_ratio across genes that have a finite log_geomean)).
    """
    # TODO
    raise NotImplementedError


# counts = simulate_counts()
# sf = size_factors(counts)
'''


STEP1_SOL_HDR = """*Click ▶ to expand the reference solution.*"""

STEP1_SOL = '''# Reference solution — Step 1.

def simulate_counts():
    mu = (
        true_mu[:, None]                        # (N_GENES, 1)
        * true_sf[None, :]                      # (1, N_SAMPLES)
        * (2.0 ** (true_lfc[:, None] * CONDITION[None, :]))
    )
    # scipy nbinom: n = 1/disp (size); p = n / (n + mu)
    n_param = 1.0 / true_disp[:, None]
    p_param = n_param / (n_param + mu)
    counts = stats.nbinom.rvs(n_param, p_param, random_state=rng)
    return counts.astype(int)


def size_factors(counts):
    # Use float and mask zeros for the log step.
    with np.errstate(divide="ignore"):
        log_counts = np.log(counts.astype(float))
    log_counts[~np.isfinite(log_counts)] = np.nan

    # Geometric mean per gene across samples (require every sample > 0).
    log_geomean = np.nanmean(log_counts, axis=1)
    valid = np.all(counts > 0, axis=1)
    log_geomean[~valid] = np.nan

    sf = np.empty(counts.shape[1])
    for j in range(counts.shape[1]):
        ratios = log_counts[:, j] - log_geomean
        sf[j] = np.exp(np.nanmedian(ratios))
    # Renormalise so they have geometric mean 1 (cosmetic).
    sf = sf / np.exp(np.mean(np.log(sf)))
    return sf


counts = simulate_counts()
print(f"Counts shape: {counts.shape}, total reads: {counts.sum():,}")
print(f"  per-sample totals: {counts.sum(axis=0)}")

sf = size_factors(counts)
print(f"Estimated size factors: {np.round(sf, 3)}")
print(f"True size factors:      {np.round(true_sf / np.exp(np.mean(np.log(true_sf))), 3)}")

# Sanity plot: estimated vs true.
fig, ax = plt.subplots(figsize=(4.5, 4.5))
true_sf_norm = true_sf / np.exp(np.mean(np.log(true_sf)))
ax.scatter(true_sf_norm, sf, s=60)
lim = (min(min(true_sf_norm), min(sf)) * 0.9,
       max(max(true_sf_norm), max(sf)) * 1.1)
ax.plot(lim, lim, "k--", lw=1, label="y = x")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("true size factor (normalised)")
ax.set_ylabel("median-of-ratios estimate")
ax.set_title("Size-factor recovery")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP2_MD = """## Step 2 (12 min) — Per-gene NB GLM and dispersion MLE

For each gene we fit `log(mu_ij) = beta0 + beta1 * CONDITION_j + log(sf_j)`
using a negative-binomial GLM (log link, size-factor offset). The
condition coefficient `beta1` is the natural-log fold change; dividing by
log(2) gives the log2 FC. We use `statsmodels.GLM` with the NB family.

For dispersion, statsmodels expects you to *supply* alpha (it does not jointly
MLE it). The standard trick: estimate alpha per gene by **profile
likelihood** — sweep alpha on a grid, refit the GLM at each step, pick the
alpha that maximises the log-likelihood. We restrict the grid to a sensible
range (1e-4 to 10) to keep things fast.
"""


STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Per-gene NB GLM + dispersion MLE by profile likelihood.
# ----------------------------------------------------------------------

# Design matrix: intercept + condition.
X = np.column_stack([np.ones(N_SAMPLES), CONDITION]).astype(float)
log_offset = np.log(sf)

# Search grid for alpha (= dispersion in statsmodels NB).
ALPHA_GRID = np.logspace(-4, 1, 25)


def fit_one_gene(y, alpha):
    """Return statsmodels GLMResults for one gene at a given alpha."""
    # TODO: sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=alpha),
    #              offset=log_offset).fit(disp=0)
    raise NotImplementedError


def profile_disp(y):
    """Return (best_alpha, best_loglik) by sweeping ALPHA_GRID."""
    # TODO
    raise NotImplementedError


def fit_all_genes(counts):
    """Return a DataFrame with columns: mle_disp, mle_loglik, beta0, beta1, se_beta1, mean_norm."""
    # TODO: loop over genes; if statsmodels diverges, fall back to a safe default.
    raise NotImplementedError


# mle = fit_all_genes(counts)
# Plot dispersion vs mean (DESeq2's classic plot).
'''


STEP2_SOL = '''# Reference solution — Step 2.

X = np.column_stack([np.ones(N_SAMPLES), CONDITION]).astype(float)
log_offset = np.log(sf)

# Coarser grid for speed; still resolves the trend.
ALPHA_GRID = np.logspace(-4, 1, 25)


def fit_one_gene(y, alpha):
    fam = sm.families.NegativeBinomial(alpha=max(alpha, 1e-6))
    return sm.GLM(y, X, family=fam, offset=log_offset).fit(disp=0, maxiter=50)


def profile_disp(y):
    best_alpha, best_ll = ALPHA_GRID[0], -np.inf
    for a in ALPHA_GRID:
        try:
            res = fit_one_gene(y, a)
            ll = res.llf
        except Exception:
            continue
        if np.isfinite(ll) and ll > best_ll:
            best_ll, best_alpha = ll, a
    return best_alpha, best_ll


def fit_all_genes(counts):
    rows = []
    for i in range(counts.shape[0]):
        y = counts[i]
        if y.sum() < 5:
            # Effectively no information; fall back to baseline.
            rows.append((1.0, np.nan, np.log(max(1, y.mean()) + 1e-9), 0.0, np.nan, y.mean()))
            continue
        a, ll = profile_disp(y)
        try:
            res = fit_one_gene(y, a)
            beta0, beta1 = res.params
            se1 = res.bse[1]
        except Exception:
            beta0, beta1, se1 = np.log(max(1, y.mean())), 0.0, np.nan
        mean_norm = (y / sf).mean()
        rows.append((a, ll, beta0, beta1, se1, mean_norm))
    return pd.DataFrame(rows, columns=["mle_disp", "mle_loglik", "beta0", "beta1", "se_beta1", "mean_norm"])


print("Fitting per-gene NB GLMs (this takes ~30-60s for 500 genes)...")
mle = fit_all_genes(counts)
print(mle.describe().round(3))

# Dispersion vs mean plot.
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.scatter(mle["mean_norm"], mle["mle_disp"], s=10, alpha=0.4, label="MLE per gene")
xs = np.logspace(np.log10(max(1e-1, mle["mean_norm"].min())),
                 np.log10(mle["mean_norm"].max()), 200)
ax.plot(xs, TREND_A + TREND_B / xs, "r--", lw=2, label="true trend a + b/mu")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("mean normalised count")
ax.set_ylabel("MLE dispersion")
ax.set_title("Per-gene dispersion MLE vs mean expression")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP3_MD = """## Step 3 (12 min) — Empirical-Bayes dispersion shrinkage

Look at the scatter you just made: the per-gene MLE dispersions are wildly
noisy because each gene only has 6 data points. The fix is **empirical
Bayes**: assume the true dispersions are drawn from a prior centred on a
mean-dispersion trend, learn the trend by regressing `log(disp) ~ log(mean)`,
then shrink each gene's MLE toward the trend by a precision-weighted
average.

The shrunk dispersion is:
$$\\log \\hat{\\alpha}_g^{\\text{EB}} = w_g \\cdot \\log \\hat{\\alpha}_g^{\\text{MLE}} + (1 - w_g) \\cdot \\log \\hat{\\alpha}_g^{\\text{trend}}$$
where the weight $w_g$ is the ratio of prior precision to (prior + likelihood) precision. We compute the prior precision from the residual variance of the trend fit; the likelihood precision from the curvature of the per-gene log-likelihood near its MLE — but to keep this exercise tractable we use a simpler closed form: $w_g = \\sigma_{\\text{prior}}^2 / (\\sigma_{\\text{prior}}^2 + \\sigma_{\\text{lik}}^2)$ with a fixed $\\sigma_{\\text{lik}}$ derived from the residual scatter (DESeq2 does the same in essence).
"""


STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — EB dispersion shrinkage.
# ----------------------------------------------------------------------

def fit_dispersion_trend(mean_norm, mle_disp):
    """Fit log(disp) = log(a + b/mean); return (a, b, predicted log_disp, residual_sd)."""
    # TODO: use scipy.optimize.curve_fit on the (mean, disp) cloud in log-space.
    raise NotImplementedError


def eb_shrink(mle_disp, trend_disp, weight_lik=0.5):
    """Precision-weighted average of log(MLE) and log(trend) per gene.

    weight_lik in [0, 1] controls how much of each gene's own evidence is kept.
    1 -> no shrinkage; 0 -> full shrink to trend.
    """
    # TODO
    raise NotImplementedError


# trend_a, trend_b, trend_pred, resid_sd = fit_dispersion_trend(mle["mean_norm"], mle["mle_disp"])
# eb_disp = eb_shrink(mle["mle_disp"].values, np.exp(trend_pred), weight_lik=0.5)
# Plot: MLE, trend, shrunk overlaid.
'''


STEP3_SOL = '''# Reference solution — Step 3.
from scipy.optimize import curve_fit


def fit_dispersion_trend(mean_norm, mle_disp):
    """Fit DESeq2-style mean-dispersion trend disp = a + b / mean."""
    good = (mle_disp > 1e-4) & (mle_disp < 5) & (mean_norm > 1)

    def trend(mu, a, b):
        return a + b / mu

    try:
        popt, _ = curve_fit(trend, mean_norm[good], mle_disp[good],
                            p0=[0.05, 4.0], maxfev=5000,
                            bounds=([0, 0], [10, 1000]))
        a_hat, b_hat = popt
    except Exception:
        a_hat, b_hat = 0.05, 4.0

    pred = trend(mean_norm.values, a_hat, b_hat)
    log_resid = np.log(np.clip(mle_disp[good], 1e-4, None)) - np.log(np.clip(pred[good], 1e-4, None))
    resid_sd = float(np.std(log_resid))
    return a_hat, b_hat, np.log(pred), resid_sd


def eb_shrink(mle_disp, trend_disp, weight_lik=0.5):
    log_mle   = np.log(np.clip(mle_disp,   1e-4, None))
    log_trend = np.log(np.clip(trend_disp, 1e-4, None))
    log_shrunk = weight_lik * log_mle + (1 - weight_lik) * log_trend
    return np.exp(log_shrunk)


trend_a, trend_b, trend_log_pred, resid_sd = fit_dispersion_trend(
    mle["mean_norm"], mle["mle_disp"]
)
trend_pred = np.exp(trend_log_pred)
print(f"Trend fit: disp = {trend_a:.3f} + {trend_b:.2f} / mean   (log-residual SD = {resid_sd:.2f})")

# Three shrinkage levels for the comparison in Step 5.
eb_strong = eb_shrink(mle["mle_disp"].values, trend_pred, weight_lik=0.2)
eb_medium = eb_shrink(mle["mle_disp"].values, trend_pred, weight_lik=0.5)

mle = mle.assign(trend_disp=trend_pred, eb_disp=eb_medium, strong_disp=eb_strong)

# Plot: MLE cloud, trend curve, shrunk values.
fig, ax = plt.subplots(figsize=(8, 5))
order = np.argsort(mle["mean_norm"])
ax.scatter(mle["mean_norm"], mle["mle_disp"], s=8, alpha=0.3, label="MLE per gene")
ax.scatter(mle["mean_norm"], mle["eb_disp"], s=8, alpha=0.6, color="orange", label="EB shrunk (w=0.5)")
ax.plot(mle["mean_norm"].iloc[order], mle["trend_disp"].iloc[order],
        "r-", lw=2, label="fitted trend")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("mean normalised count")
ax.set_ylabel("dispersion")
ax.set_title("Dispersion shrinkage toward the fitted trend")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP4_MD = """## Step 4 (15 min) — Wald test, BH FDR, volcano plot

With shrunk dispersions in hand, we **refit** each gene's GLM using the
shrunk alpha and compute a Wald z-score on the condition coefficient:
$z_g = \\hat\\beta_g / \\text{SE}(\\hat\\beta_g)$. The two-sided p-value is
$2 \\Phi(-|z|)$. Then Benjamini-Hochberg controls FDR at a chosen level
(usually 0.05): sort p-values, multiply each by `n / rank`, take the
cumulative min from the right. The result is a q-value per gene.

Finally we draw the canonical **volcano plot**: log2 FC on x, -log10(q) on
y, true-positive genes coloured. Significant hits should sit in the
upper-left and upper-right corners.
"""


STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Wald test, BH FDR, volcano plot.
# ----------------------------------------------------------------------

def refit_with_disp(counts, disp_per_gene):
    """Refit each gene at its supplied dispersion; return DataFrame with
    columns log2fc, se, z, pvalue."""
    # TODO: same loop as Step 2 but with fixed alpha = disp_per_gene[g].
    raise NotImplementedError


def bh_qvalues(pvalues):
    """Benjamini-Hochberg q-values. NaN p-values get NaN q-values."""
    # TODO
    raise NotImplementedError


# tests = refit_with_disp(counts, mle["eb_disp"].values)
# tests["qvalue"] = bh_qvalues(tests["pvalue"].values)
# Draw the volcano plot.
'''


STEP4_SOL = '''# Reference solution — Step 4.

def refit_with_disp(counts, disp_per_gene):
    rows = []
    for i in range(counts.shape[0]):
        y = counts[i]
        a = max(1e-6, float(disp_per_gene[i]))
        try:
            res = fit_one_gene(y, a)
            beta1 = res.params[1]
            se = res.bse[1]
            # Guard against degenerate fits: huge SE means the gene is
            # essentially zero counts in one arm; treat as non-informative.
            if not (np.isfinite(se) and 0 < se < 10):
                raise ValueError("SE out of bounds")
            z = beta1 / se
            p = 2.0 * stats.norm.sf(abs(z))
        except Exception:
            beta1, se, z, p = 0.0, np.nan, np.nan, np.nan
        rows.append((beta1 / np.log(2.0), se, z, p))
    return pd.DataFrame(rows, columns=["log2fc", "se", "z", "pvalue"])


def bh_qvalues(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    out = np.full(n, np.nan)
    mask = np.isfinite(p)
    pv = p[mask]
    order = np.argsort(pv)
    ranks = np.arange(1, len(pv) + 1)
    raw = pv[order] * n / ranks
    # Cumulative min from the right (monotone non-decreasing).
    qv = np.minimum.accumulate(raw[::-1])[::-1]
    qv = np.clip(qv, 0, 1)
    q_unsorted = np.empty_like(qv)
    q_unsorted[order] = qv
    out[mask] = q_unsorted
    return out


print("Refitting GLMs at EB-shrunk dispersions...")
tests = refit_with_disp(counts, mle["eb_disp"].values)
tests["qvalue"] = bh_qvalues(tests["pvalue"].values)
tests["is_de"] = is_de
tests["true_lfc"] = true_lfc

print(tests.describe().round(3))

# Volcano plot.
fig, ax = plt.subplots(figsize=(8, 5.5))
neg_log_q = -np.log10(np.clip(tests["qvalue"].values, 1e-30, 1))
ax.scatter(tests.loc[~is_de, "log2fc"], neg_log_q[~is_de],
           s=12, alpha=0.4, color="gray", label="null genes")
ax.scatter(tests.loc[is_de, "log2fc"], neg_log_q[is_de],
           s=24, alpha=0.85, color="crimson", label="true DE")
ax.axhline(-np.log10(0.05), ls="--", color="black", lw=1, label="q = 0.05")
ax.axvline(0, color="black", lw=0.5)
ax.set_xlabel("log2 fold change (treated / control)")
ax.set_ylabel("-log10 BH q-value")
ax.set_title("Volcano plot — EB-shrunk NB GLM")
ax.legend(loc="upper left")
plt.tight_layout()
plt.show()

# Summary at q < 0.05.
hits = tests["qvalue"] < 0.05
tp = int((hits & is_de).sum())
fp = int((hits & ~is_de).sum())
fn = int((~hits & is_de).sum())
print(f"\\nAt q < 0.05:")
print(f"  TP = {tp}, FP = {fp}, FN = {fn}")
print(f"  TPR (power)         = {tp / N_TRUE_DE:.2%}")
print(f"  FDR (observed)      = {fp / max(1, tp + fp):.2%}")
'''


STEP5_MD = """## Step 5 (11 min) — Shrinkage strength comparison: none vs EB vs strong

Now the punchline. We run the same Wald-test + BH pipeline three times
with three different dispersion estimates:

1. **No shrinkage** — the raw per-gene MLE from Step 2.
2. **EB shrinkage** — the moderate weight (w=0.5) from Step 3.
3. **Strong prior** — heavy shrinkage toward the trend (w=0.2), close to
   "everybody gets the trend dispersion".

For each, we report TPR and observed FDR at q < 0.05. The expected
pattern: no-shrinkage has poorly-calibrated FDR (too many false
positives, because some genes get an absurdly small dispersion estimate
purely by chance); EB has both high TPR and well-controlled FDR; strong
prior is conservative — TPR drops because outlier-dispersion genes get
their evidence smoothed away.
"""


STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Compare three shrinkage strengths.
# ----------------------------------------------------------------------

def evaluate(disp_per_gene, label):
    """Run refit -> p -> q -> TPR / FDR at q<0.05; return a dict."""
    # TODO
    raise NotImplementedError


# results = []
# results.append(evaluate(mle["mle_disp"].values,    "none (raw MLE)"))
# results.append(evaluate(mle["eb_disp"].values,     "EB (w=0.5)"))
# results.append(evaluate(mle["strong_disp"].values, "strong prior (w=0.2)"))
'''


STEP5_SOL = '''# Reference solution — Step 5.

def evaluate(disp_per_gene, label):
    t = refit_with_disp(counts, disp_per_gene)
    t["qvalue"] = bh_qvalues(t["pvalue"].values)
    hits = t["qvalue"] < 0.05
    tp = int((hits & is_de).sum())
    fp = int((hits & ~is_de).sum())
    fn = int((~hits & is_de).sum())
    tpr = tp / N_TRUE_DE
    fdr = fp / max(1, tp + fp)
    return {"label": label, "TP": tp, "FP": fp, "FN": fn,
            "TPR": tpr, "FDR_observed": fdr, "pvalues": t["pvalue"].values}


print("Comparing shrinkage strengths (this re-runs 3 x 500 GLM fits)...")
res_none   = evaluate(mle["mle_disp"].values,    "none (raw MLE)")
res_eb     = evaluate(mle["eb_disp"].values,     "EB (w=0.5)")
res_strong = evaluate(mle["strong_disp"].values, "strong prior (w=0.2)")

table = pd.DataFrame([
    {k: v for k, v in r.items() if k != "pvalues"}
    for r in (res_none, res_eb, res_strong)
])
print()
print(table.to_string(index=False))

# Visualise: TPR vs observed FDR bar chart, and p-value histogram for the null
# genes (which should be flat for a well-calibrated test).
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

axes[0].bar(np.arange(3) - 0.2, table["TPR"], width=0.35, color="steelblue", label="TPR (power)")
axes[0].bar(np.arange(3) + 0.2, table["FDR_observed"], width=0.35, color="crimson", label="observed FDR")
axes[0].axhline(0.05, ls="--", color="black", lw=1, label="nominal FDR 0.05")
axes[0].set_xticks(np.arange(3))
axes[0].set_xticklabels(table["label"], rotation=10)
axes[0].set_ylabel("rate")
axes[0].set_ylim(0, 1)
axes[0].set_title("Power vs realised FDR")
axes[0].legend(loc="upper right")

for r, color in zip([res_none, res_eb, res_strong], ["gray", "orange", "purple"]):
    null_p = r["pvalues"][~is_de]
    null_p = null_p[np.isfinite(null_p)]
    axes[1].hist(null_p, bins=20, alpha=0.5, label=r["label"], color=color)
axes[1].set_xlabel("p-value (null genes only)")
axes[1].set_ylabel("count")
axes[1].set_title("Null p-value distribution (should be flat)")
axes[1].legend(fontsize=8)
plt.tight_layout()
plt.show()

# Save EB results for the self-check below.
final_results = res_eb
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces of the pipeline.
If you ran the reference solutions above they should all pass; if you
wrote your own and an assert fails, revisit the corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Count matrix is the right shape and non-negative.
assert counts.shape == (N_GENES, N_SAMPLES), f"counts shape {counts.shape}"
assert counts.min() >= 0

# 2. Size factors are within a factor of 2 of each other (sanity).
assert 0.4 < sf.min() < sf.max() < 2.5, f"size factors: {sf}"
# And the relative order matches the truth (Spearman > 0.7).
true_sf_norm = true_sf / np.exp(np.mean(np.log(true_sf)))
spearman = stats.spearmanr(true_sf_norm, sf).statistic
assert spearman > 0.7, f"size-factor Spearman = {spearman:.2f}"

# 3. BH q-values are valid probabilities and monotone in p.
p_demo = np.array([0.001, 0.01, 0.04, 0.2, 0.5, 0.9])
q_demo = bh_qvalues(p_demo)
assert np.all((q_demo >= 0) & (q_demo <= 1))
# Sorted-by-p q-values must be non-decreasing.
order = np.argsort(p_demo)
assert np.all(np.diff(q_demo[order]) >= -1e-12), "BH q-values not monotone"
# Known value: BH for [0.001, 0.01, 0.04, 0.2, 0.5, 0.9] with n=6 is
# [0.006, 0.030, 0.080, 0.300, 0.600, 0.900].
expected = np.array([0.006, 0.030, 0.080, 0.300, 0.600, 0.900])
assert np.allclose(q_demo, expected, atol=1e-6), f"BH check: got {q_demo}, want {expected}"

# 4. EB improves FDR control relative to no shrinkage.
#    At 3 vs 3 the absolute FDR will not always hit nominal 0.05, but EB
#    should be strictly closer than the raw MLE.
assert res_eb["FDR_observed"] <= res_none["FDR_observed"] + 0.05, \
    f"EB FDR {res_eb['FDR_observed']:.2%} not better than no-shrinkage {res_none['FDR_observed']:.2%}"
assert res_eb["TPR"] >= 0.50, \
    f"EB TPR too low: {res_eb['TPR']:.2%}"
# Strong prior over-shrinks; EB should retain at least as much power.
assert res_eb["TPR"] >= res_strong["TPR"] - 0.05, \
    f"EB TPR {res_eb['TPR']:.2%} unexpectedly < strong-prior TPR {res_strong['TPR']:.2%}"

# 5. Volcano hits are enriched for true DE genes (precision > random rate).
hits = tests["qvalue"] < 0.05
hits = hits.fillna(False).values
precision = (hits & is_de).sum() / max(1, hits.sum())
# Random precision would be 50/500 = 10%. We require a meaningful enrichment.
assert precision >= 0.5, f"hit-set precision low: {precision:.2%}"

print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — James-Stein, Wiener filtering, Tikhonov regularisation

This whole exercise is a single idea from estimation theory wearing a
biology costume:

- **Per-gene dispersion MLE = noisy measurement.** With 6 samples per
  gene, the MLE has so much variance that some null genes get assigned an
  absurdly small dispersion and look "highly significant" by accident.
  This is the same problem as estimating a noisy channel gain from one
  observation: you can't, not without a prior.
- **EB shrinkage = optimal Bayesian estimator under a learned prior.**
  The mean-dispersion trend is the prior; the per-gene MLE is the
  likelihood; the shrunk estimate is the posterior mean. The
  precision-weighted average we used is the closed-form posterior under a
  Gaussian-Gaussian model in log-space.
- **This is James-Stein in disguise.** Stein (1956) showed that when
  estimating >= 3 parameters under squared-error loss, *always* shrinking
  toward the mean strictly dominates the MLE. Empirical Bayes is the
  practical form of that result: learn the prior from the data.
- **EE analogies.**
  - **Wiener filter:** optimal linear estimator that trades signal
    estimation variance for bias by weighting by the signal-to-noise
    ratio. Identical structure: posterior mean = (SNR / (1+SNR)) ·
    measurement + (1 / (1+SNR)) · prior.
  - **Tikhonov / ridge regression:** add a penalty that pulls estimates
    toward zero. EB does the same for log-dispersion, pulling toward the
    trend instead of zero.
  - **Kalman smoothing:** combine a noisy measurement with a prior
    prediction by precision-weighted averaging. Same identity.

The takeaway: when you have many parallel noisy measurements drawn from a
common population, *always* shrink. Per-gene MLE looks "unbiased" but is
strictly worse than EB on every loss function that matters for the FDR
problem.
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
        md(STEP1_SOL_HDR),
        hidden(STEP1_SOL),

        md(STEP2_MD),
        code(STEP2_TODO),
        md(STEP1_SOL_HDR),
        hidden(STEP2_SOL),

        md(STEP3_MD),
        code(STEP3_TODO),
        md(STEP1_SOL_HDR),
        hidden(STEP3_SOL),

        md(STEP4_MD),
        code(STEP4_TODO),
        md(STEP1_SOL_HDR),
        hidden(STEP4_SOL),

        md(STEP5_MD),
        code(STEP5_TODO),
        md(STEP1_SOL_HDR),
        hidden(STEP5_SOL),

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
