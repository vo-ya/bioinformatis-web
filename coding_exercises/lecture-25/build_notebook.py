"""Build exercise.ipynb for L25 — Causal Inference and Mendelian Randomisation.

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


TITLE_MD = """# L25 — Causal Inference and Mendelian Randomisation

In this exercise you implement two-sample Mendelian randomisation from
scratch on 100 simulated SNP instruments. You will fit the three
workhorse estimators — **IVW**, **MR-Egger**, **weighted median** —
identify outliers via Cook's distance, run a Steiger directionality
check, and read off how triangulation flags hidden horizontal
pleiotropy.
"""


AHA_MD = """> **Aha.** Mendelian randomisation is an **errors-in-variables**
> regression where the instrument-exposure (β_ZX) and instrument-outcome
> (β_ZY) coefficients live on the axes. Under the three IV assumptions
> the regression slope **is** the causal effect — no covariates needed.
> IVW forces the intercept to 0; MR-Egger lets it float and absorbs
> directional pleiotropy at the cost of power; the weighted median
> stays consistent so long as ≥ 50% of instruments are valid. When all
> three agree you triangulate; when they diverge, something is leaking
> through a back door.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 matplotlib==3.8.4 -q
"""


IMPORTS = """import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Deterministic for the whole notebook.
RNG_SEED = 25
np.random.seed(RNG_SEED)

# True causal effect of exposure X on outcome Y -- the number we want to recover.
TRUE_BETA = 0.5

# Instrument counts.
N_VALID       = 80   # SNPs satisfying all three IV assumptions
N_PLEIOTROPIC = 20   # SNPs with horizontal pleiotropy (violate exclusion)
N_SNPS        = N_VALID + N_PLEIOTROPIC

# Sample sizes for the two GWAS (two-sample MR).
N_X = 20_000        # exposure GWAS
N_Y = 25_000        # outcome GWAS

print(f"Target: recover TRUE_BETA = {TRUE_BETA} from {N_SNPS} instruments "
      f"({N_VALID} valid + {N_PLEIOTROPIC} pleiotropic).")
"""


STEP1_MD = """## Step 1 (10 min) — Simulate two-sample summary statistics

Real MR papers receive **summary statistics** from two independent
GWAS: per-SNP marginal regression coefficients β̂ and standard errors
SE for the exposure (X) and the outcome (Y). We will fabricate them
under a known causal structure, then forget the truth and recover it.

For each SNP $j$ with minor-allele frequency $p_j$:

- True effect on exposure: $\\gamma_j$ drawn from $\\mathcal{N}(0, 0.05^2)$.
- True direct effect on outcome via pleiotropy: $\\alpha_j = 0$ for valid SNPs;
  drawn from $\\mathcal{N}(0.1, 0.04^2)$ for the 20 pleiotropic SNPs
  (directional pleiotropy — the MR-Egger intercept should pick it up).
- True total effect on outcome: $\\Gamma_j = \\beta \\gamma_j + \\alpha_j$.

The reported $\\hat{\\beta}_{ZX}$ and $\\hat{\\beta}_{ZY}$ are these true
quantities plus Gaussian noise with SE $\\propto 1 / \\sqrt{N \\cdot 2 p (1-p)}$.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Simulate (beta_ZX, SE_ZX) and (beta_ZY, SE_ZY) per SNP.
# ----------------------------------------------------------------------

def simulate_mr_summary(
    n_valid: int = N_VALID,
    n_pleio: int = N_PLEIOTROPIC,
    beta: float = TRUE_BETA,
    n_x: int = N_X,
    n_y: int = N_Y,
    seed: int = RNG_SEED,
) -> pd.DataFrame:
    """Return a DataFrame with columns:
        snp, maf, gamma, alpha, valid,
        beta_zx, se_zx, beta_zy, se_zy
    where gamma/alpha are the TRUE per-SNP exposure and pleiotropy effects,
    and the beta_z*/se_z* columns are the observed summary stats.

    Valid SNPs (alpha == 0) must occupy the first n_valid rows; pleiotropic SNPs
    (alpha drawn from N(0.1, 0.04^2)) the remaining rows.
    """
    # TODO:
    # 1. Draw MAFs uniformly in [0.05, 0.5].
    # 2. Draw gamma_j ~ N(0, 0.05^2).
    # 3. Build alpha vector (zeros for valid, N(0.1, 0.04^2) for pleiotropic).
    # 4. Gamma_j = beta * gamma_j + alpha_j.
    # 5. SE for a marginal allele-dosage regression: 1 / sqrt(N * 2 * p * (1 - p)).
    # 6. Observed beta_z* = true + N(0, SE^2).
    raise NotImplementedError


df = None  # df = simulate_mr_summary()
'''


STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution -- Step 1.

def simulate_mr_summary(
    n_valid: int = N_VALID,
    n_pleio: int = N_PLEIOTROPIC,
    beta: float = TRUE_BETA,
    n_x: int = N_X,
    n_y: int = N_Y,
    seed: int = RNG_SEED,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n_snps = n_valid + n_pleio

    maf = rng.uniform(0.05, 0.5, size=n_snps)
    gamma = rng.normal(0.0, 0.05, size=n_snps)

    alpha = np.zeros(n_snps)
    alpha[n_valid:] = rng.normal(0.10, 0.04, size=n_pleio)  # directional pleiotropy
    valid = np.zeros(n_snps, dtype=bool)
    valid[:n_valid] = True

    Gamma = beta * gamma + alpha  # true total effect on outcome

    # Per-allele-dosage regression: Var(b) = sigma^2 / (N * Var(G)); Var(G) = 2 p (1 - p).
    # Take sigma = 1 (standardised traits) so SE = 1 / sqrt(N * 2 p (1 - p)).
    se_zx = 1.0 / np.sqrt(n_x * 2.0 * maf * (1.0 - maf))
    se_zy = 1.0 / np.sqrt(n_y * 2.0 * maf * (1.0 - maf))

    beta_zx = gamma + rng.normal(0.0, se_zx)
    beta_zy = Gamma + rng.normal(0.0, se_zy)

    df = pd.DataFrame({
        "snp":    [f"rs{1_000_000 + i}" for i in range(n_snps)],
        "maf":    maf,
        "gamma":  gamma,
        "alpha":  alpha,
        "valid":  valid,
        "beta_zx": beta_zx,
        "se_zx":   se_zx,
        "beta_zy": beta_zy,
        "se_zy":   se_zy,
    })
    return df


df = simulate_mr_summary()

# Sanity print: show a few rows and the F-statistic (instrument strength).
print(df.head().to_string(index=False))
F = (df["beta_zx"] / df["se_zx"]) ** 2
print(f"\\nMean F-statistic across instruments: {F.mean():6.1f} "
      f"(>>10 = strong instruments)")
print(f"  n valid:        {df['valid'].sum()}")
print(f"  n pleiotropic:  {(~df['valid']).sum()}")

# Wald scatter: beta_zy vs beta_zx. Slope ~ TRUE_BETA on valid SNPs.
fig, ax = plt.subplots(figsize=(6.5, 5))
ax.errorbar(df.loc[df["valid"], "beta_zx"], df.loc[df["valid"], "beta_zy"],
            xerr=df.loc[df["valid"], "se_zx"], yerr=df.loc[df["valid"], "se_zy"],
            fmt="o", ms=4, alpha=0.55, color="#3367d6", label="valid (n=80)")
ax.errorbar(df.loc[~df["valid"], "beta_zx"], df.loc[~df["valid"], "beta_zy"],
            xerr=df.loc[~df["valid"], "se_zx"], yerr=df.loc[~df["valid"], "se_zy"],
            fmt="s", ms=5, alpha=0.85, color="#d62728", label="pleiotropic (n=20)")
xs = np.linspace(df["beta_zx"].min(), df["beta_zx"].max(), 50)
ax.plot(xs, TRUE_BETA * xs, "k--", lw=1.5, label=f"truth (slope = {TRUE_BETA})")
ax.axhline(0, color="grey", lw=0.5)
ax.axvline(0, color="grey", lw=0.5)
ax.set_xlabel(r"$\\hat\\beta_{ZX}$ (instrument -> exposure)")
ax.set_ylabel(r"$\\hat\\beta_{ZY}$ (instrument -> outcome)")
ax.set_title("Wald scatter -- two-sample MR summary stats")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP2_MD = """## Step 2 (12 min) — Inverse-variance-weighted estimator

The per-SNP **Wald ratio** $\\hat\\beta_j = \\hat\\beta_{ZY,j} / \\hat\\beta_{ZX,j}$
is a noisy estimator of the causal effect. Pool them with weights
$w_j = \\hat\\beta_{ZX,j}^2 / \\mathrm{SE}(\\hat\\beta_{ZY,j})^2$ (delta-method
approximation):

$$\\hat\\beta_{\\text{IVW}} = \\frac{\\sum_j w_j \\hat\\beta_j}{\\sum_j w_j}, \\qquad
\\mathrm{SE}(\\hat\\beta_{\\text{IVW}}) = \\sqrt{\\frac{1}{\\sum_j w_j}}.$$

Equivalently, IVW is the **weighted least squares** slope of
$\\hat\\beta_{ZY}$ on $\\hat\\beta_{ZX}$ **through the origin**, with weights
$1 / \\mathrm{SE}(\\hat\\beta_{ZY})^2$.

When *all* instruments are valid the IVW estimate is consistent. When 20%
of SNPs leak pleiotropic effects (as in our cohort), IVW becomes biased
because it has no parameter to absorb the bias.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 -- IVW estimator.
# ----------------------------------------------------------------------

def ivw_estimate(df: pd.DataFrame) -> dict:
    """Return {'beta': beta_IVW, 'se': SE, 'ci_low', 'ci_high', 'z', 'p'}.

    Weights are 1 / SE(beta_ZY)^2 in the WLS-through-origin formulation,
    which is algebraically the same as the delta-method per-SNP weighting
    when the delta-method neglects uncertainty in beta_ZX (the NO Measurement
    Error assumption -- "NOME").
    """
    # TODO: implement IVW via WLS through the origin.
    # beta_IVW = sum(w_j * beta_zx_j * beta_zy_j) / sum(w_j * beta_zx_j^2)
    # SE_IVW   = sqrt(1 / sum(w_j * beta_zx_j^2))
    raise NotImplementedError


ivw = None  # ivw = ivw_estimate(df)
'''


STEP2_SOLUTION = '''# Reference solution -- Step 2.

def ivw_estimate(df: pd.DataFrame) -> dict:
    bx = df["beta_zx"].to_numpy()
    by = df["beta_zy"].to_numpy()
    sy = df["se_zy"].to_numpy()
    w  = 1.0 / sy ** 2

    num = np.sum(w * bx * by)
    den = np.sum(w * bx * bx)
    beta = num / den
    se   = np.sqrt(1.0 / den)
    z    = beta / se
    p    = 2.0 * (1.0 - stats.norm.cdf(abs(z)))
    return {
        "method":  "IVW",
        "beta":    beta,
        "se":      se,
        "ci_low":  beta - 1.96 * se,
        "ci_high": beta + 1.96 * se,
        "z":       z,
        "p":       p,
    }


ivw = ivw_estimate(df)
print(f"IVW    beta = {ivw['beta']:.3f}  "
      f"SE = {ivw['se']:.3f}  "
      f"95% CI = [{ivw['ci_low']:.3f}, {ivw['ci_high']:.3f}]  "
      f"p = {ivw['p']:.2e}")
print(f"  bias vs truth: {ivw['beta'] - TRUE_BETA:+.3f}  "
      f"(directional pleiotropy pulls IVW off {TRUE_BETA} -- sign depends on the draw)")

# IVW restricted to the *valid* subset should be unbiased -- a useful diagnostic.
ivw_clean = ivw_estimate(df[df["valid"]].reset_index(drop=True))
print(f"IVW (valid-only oracle) beta = {ivw_clean['beta']:.3f}  "
      f"-- close to {TRUE_BETA}, confirming the simulator is well-specified.")
'''


STEP3_MD = """## Step 3 (13 min) — MR-Egger with free intercept + outlier detection

MR-Egger relaxes the **intercept-through-origin** constraint:

$$\\hat\\beta_{ZY,j} = \\alpha_0 + \\beta_{\\text{Egger}}\\,\\hat\\beta_{ZX,j} + \\epsilon_j,
\\quad \\epsilon_j \\sim \\mathcal{N}(0, \\mathrm{SE}(\\hat\\beta_{ZY,j})^2).$$

The intercept $\\alpha_0$ absorbs **average directional pleiotropy**: under
the InSIDE assumption (Instrument Strength Independent of Direct Effect)
the slope $\\beta_{\\text{Egger}}$ stays consistent for the causal effect
even with pleiotropic SNPs. A non-zero $\\alpha_0$ is itself a diagnostic
("Egger intercept test").

You will then flag outliers by **Cook's distance** on the Egger WLS fit
— SNPs with $D > 4/N$ are conventional candidates for removal. Compare
to the ground-truth `valid` flag (you won't recover all 20 pleiotropic
SNPs, but you should catch the worst offenders).
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 -- MR-Egger + Cook's distance outlier detection.
# ----------------------------------------------------------------------

def mr_egger(df: pd.DataFrame) -> dict:
    """Return {'beta', 'se', 'intercept', 'intercept_se', 'intercept_p', 'cook'}.

    Fit a weighted least squares regression
        beta_zy ~ intercept + slope * beta_zx,   weights w_j = 1 / SE(beta_zy)^2.
    Return per-SNP Cook's distances and the standard Egger intercept test.
    """
    # TODO: implement WLS by hand (or via np.linalg.lstsq on whitened design).
    # Cook's distance: D_j = (r_j^2 / (p * MSE)) * h_jj / (1 - h_jj)^2,
    # where r_j is the standardized residual, h_jj is the hat-matrix diagonal,
    # p = 2 parameters (intercept + slope).
    raise NotImplementedError


egger = None
'''


STEP3_SOLUTION = '''# Reference solution -- Step 3.

def mr_egger(df: pd.DataFrame) -> dict:
    bx = df["beta_zx"].to_numpy()
    by = df["beta_zy"].to_numpy()
    sy = df["se_zy"].to_numpy()
    n  = len(bx)

    # Weighted least squares: whiten by 1 / sy, then ordinary lstsq.
    w_sqrt = 1.0 / sy
    X = np.column_stack([np.ones(n), bx])                # (n, 2)
    Xw = X * w_sqrt[:, None]                              # whitened design
    yw = by * w_sqrt                                      # whitened response

    # Coefficients via normal equations.
    XtX = Xw.T @ Xw
    XtX_inv = np.linalg.inv(XtX)
    coef = XtX_inv @ (Xw.T @ yw)
    intercept, beta = coef

    # Residuals + MSE.
    resid_w = yw - Xw @ coef
    p_params = 2
    dof = n - p_params
    mse = (resid_w @ resid_w) / dof

    # Coefficient SEs.
    se_coef = np.sqrt(np.diag(XtX_inv) * mse)
    intercept_se, beta_se = se_coef
    intercept_z = intercept / intercept_se
    intercept_p = 2.0 * (1.0 - stats.norm.cdf(abs(intercept_z)))

    # Cook's distance (per-observation influence).
    # Hat matrix H = Xw (Xw'Xw)^-1 Xw'; we only need its diagonal.
    H_diag = np.einsum("ij,jk,ik->i", Xw, XtX_inv, Xw)
    # Studentized residuals.
    r_std = resid_w / np.sqrt(mse * (1.0 - H_diag))
    cook = (r_std ** 2) / p_params * (H_diag / (1.0 - H_diag))

    return {
        "method":       "MR-Egger",
        "beta":         beta,
        "se":           beta_se,
        "ci_low":       beta - 1.96 * beta_se,
        "ci_high":      beta + 1.96 * beta_se,
        "intercept":    intercept,
        "intercept_se": intercept_se,
        "intercept_p":  intercept_p,
        "cook":         cook,
    }


egger = mr_egger(df)
print(f"MR-Egger  beta      = {egger['beta']:.3f}  "
      f"SE = {egger['se']:.3f}  "
      f"95% CI = [{egger['ci_low']:.3f}, {egger['ci_high']:.3f}]")
print(f"          intercept = {egger['intercept']:+.4f}  "
      f"SE = {egger['intercept_se']:.4f}  "
      f"p = {egger['intercept_p']:.2e}")
print(f"  Non-zero intercept = signal of directional pleiotropy "
      f"(true mean alpha = {df.loc[~df['valid'], 'alpha'].mean():+.3f}).")

# Cook's distance: flag SNPs with D > 4/N.
threshold = 4.0 / len(df)
flagged = np.where(egger["cook"] > threshold)[0]
print(f"\\nCook's distance threshold (4/N): {threshold:.4f}")
print(f"Flagged {len(flagged)} outlier SNPs.")

# How many of the flagged SNPs are *actually* pleiotropic?
pleio_idx = np.where(~df["valid"].to_numpy())[0]
true_positives = np.intersect1d(flagged, pleio_idx)
false_positives = np.setdiff1d(flagged, pleio_idx)
print(f"  True positives (truly pleiotropic):  {len(true_positives):2d} / {len(pleio_idx)} pleiotropic SNPs")
print(f"  False positives (valid but flagged): {len(false_positives):2d} / {df['valid'].sum()} valid SNPs")

fig, ax = plt.subplots(figsize=(8, 4))
colors = np.where(df["valid"], "#3367d6", "#d62728")
ax.stem(np.arange(len(df)), egger["cook"], linefmt="grey", markerfmt=" ", basefmt=" ")
ax.scatter(np.arange(len(df)), egger["cook"], c=colors, s=18, zorder=3)
ax.axhline(threshold, color="black", lw=1.2, ls="--", label=f"4/N = {threshold:.4f}")
ax.set_xlabel("SNP index")
ax.set_ylabel("Cook's distance")
ax.set_title("Egger residual influence -- pleiotropic SNPs in red")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP4_MD = """## Step 4 (12 min) — Weighted median + estimator comparison

The weighted median is the **robust** alternative. If
$\\hat\\beta_1 \\le \\hat\\beta_2 \\le \\dots \\le \\hat\\beta_N$ are the
per-SNP Wald ratios sorted ascending with normalised weights
$\\tilde w_j$ summing to 1, the **weighted median** is the value at the
point where the cumulative weight crosses 0.5 — i.e., a linear
interpolation between the two SNPs whose cumulative weights bracket the
median.

Bowden et al. showed the weighted median is consistent as long as **at
least 50 % of the weight comes from valid instruments** — strictly
weaker than IVW (which needs 100 %).

Standard errors come from a parametric bootstrap: resample
$\\hat\\beta_{Z*} \\sim \\mathcal{N}(\\hat\\beta_{Z*}, \\mathrm{SE}_*^2)$,
recompute the weighted median on each draw, take the empirical SD.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 -- weighted median + bootstrap SE.
# ----------------------------------------------------------------------

def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    """Linearly interpolated weighted median (Bowden 2016 definition)."""
    # TODO: sort by value, build cumulative weight, find the index where cumweight
    # crosses 0.5 (after normalisation), and linearly interpolate.
    raise NotImplementedError


def weighted_median_mr(df: pd.DataFrame, n_boot: int = 1000, seed: int = RNG_SEED) -> dict:
    """Return {'beta', 'se', 'ci_low', 'ci_high'} for the weighted-median MR estimator."""
    # TODO:
    # 1. Compute per-SNP Wald ratios beta_j = beta_zy_j / beta_zx_j.
    # 2. Weights w_j = beta_zx_j^2 / SE(beta_zy_j)^2  (same as IVW weights).
    # 3. Point estimate = weighted_median(beta_j, w_j).
    # 4. Bootstrap SE: for b in range(n_boot), resample beta_z* ~ N(beta_z*, SE_z*^2)
    #    and recompute the weighted median.
    raise NotImplementedError


wmed = None
'''


STEP4_SOLUTION = '''# Reference solution -- Step 4.

def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    v_sorted = values[order]
    w_sorted = weights[order]
    cw = np.cumsum(w_sorted) / w_sorted.sum()
    # Linear interpolation between the two SNPs straddling cumulative weight = 0.5.
    idx_hi = np.searchsorted(cw, 0.5)
    if idx_hi == 0:
        return float(v_sorted[0])
    if idx_hi >= len(v_sorted):
        return float(v_sorted[-1])
    cw_lo, cw_hi = cw[idx_hi - 1], cw[idx_hi]
    v_lo, v_hi   = v_sorted[idx_hi - 1], v_sorted[idx_hi]
    if cw_hi == cw_lo:
        return float(v_lo)
    # Interpolate so that the function crosses 0.5 exactly.
    frac = (0.5 - cw_lo) / (cw_hi - cw_lo)
    return float(v_lo + frac * (v_hi - v_lo))


def weighted_median_mr(df: pd.DataFrame, n_boot: int = 1000, seed: int = RNG_SEED) -> dict:
    rng = np.random.default_rng(seed)
    bx = df["beta_zx"].to_numpy()
    by = df["beta_zy"].to_numpy()
    sx = df["se_zx"].to_numpy()
    sy = df["se_zy"].to_numpy()

    wald = by / bx
    weights = bx ** 2 / sy ** 2  # delta-method approx (NOME)

    point = weighted_median(wald, weights)

    # Bootstrap: redraw the summary stats and recompute on each iteration.
    boot = np.empty(n_boot)
    for b in range(n_boot):
        bx_b = bx + rng.normal(0.0, sx)
        by_b = by + rng.normal(0.0, sy)
        wald_b = by_b / bx_b
        w_b = bx_b ** 2 / sy ** 2
        boot[b] = weighted_median(wald_b, w_b)

    se = float(boot.std(ddof=1))
    return {
        "method":  "Weighted median",
        "beta":    point,
        "se":      se,
        "ci_low":  point - 1.96 * se,
        "ci_high": point + 1.96 * se,
    }


wmed = weighted_median_mr(df, n_boot=1000)
print(f"Weighted median  beta = {wmed['beta']:.3f}  "
      f"SE = {wmed['se']:.3f}  "
      f"95% CI = [{wmed['ci_low']:.3f}, {wmed['ci_high']:.3f}]")
print(f"  Robust to up to 50% invalid instruments (we have 20% -- comfortably under).")

# Assemble the three-estimator comparison table.
results = pd.DataFrame([
    {"method": "IVW",             "beta": ivw["beta"],   "se": ivw["se"],   "ci_low": ivw["ci_low"],   "ci_high": ivw["ci_high"]},
    {"method": "MR-Egger",        "beta": egger["beta"], "se": egger["se"], "ci_low": egger["ci_low"], "ci_high": egger["ci_high"]},
    {"method": "Weighted median", "beta": wmed["beta"],  "se": wmed["se"],  "ci_low": wmed["ci_low"],  "ci_high": wmed["ci_high"]},
])
print("\\nEstimator triangulation:")
print(results.to_string(index=False, float_format=lambda x: f"{x:6.3f}"))

# Forest plot.
fig, ax = plt.subplots(figsize=(7, 3.2))
y_pos = np.arange(len(results))[::-1]
ax.errorbar(results["beta"], y_pos,
            xerr=[results["beta"] - results["ci_low"], results["ci_high"] - results["beta"]],
            fmt="o", ms=8, color="#3367d6", capsize=4)
ax.axvline(TRUE_BETA, color="black", lw=1.2, ls="--", label=f"truth = {TRUE_BETA}")
ax.axvline(0, color="grey", lw=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels(results["method"])
ax.set_xlabel(r"causal effect estimate $\\hat\\beta$")
ax.set_title("MR triangulation -- 3 estimators, 95% CIs")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP5_MD = """## Step 5 (13 min) — Steiger directionality + sensitivity triangulation

The third IV assumption — **exclusion**: $Z$ acts on $Y$ only through
$X$ — is not testable from data. But we *can* test **directionality**:
which trait is the SNP a stronger instrument for? The **Steiger test**
compares $R^2_{ZX}$ vs $R^2_{ZY}$ per SNP; if $R^2_{ZY} > R^2_{ZX}$,
the SNP is plausibly acting on $Y$ first (reverse causation), and the
SNP should be excluded.

We approximate $R^2 \\approx \\hat\\beta^2 \\cdot 2p(1-p)$ (the variance
explained by an additive allele-dosage regressor on a standardised
trait). Then we:

1. Drop any SNP that fails the Steiger test (forward direction).
2. Drop any SNP flagged by Cook's distance from the Egger fit.
3. Re-run all three estimators on the cleaned cohort.
4. Read off the triangulation: do IVW, Egger, weighted median now agree?
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 -- Steiger directionality filter + post-filter triangulation.
# ----------------------------------------------------------------------

def steiger_filter(df: pd.DataFrame) -> np.ndarray:
    """Return a boolean mask: True for SNPs passing the Steiger forward test
    (R^2_ZX > R^2_ZY, where R^2 ~= 2 p (1-p) beta^2)."""
    # TODO
    raise NotImplementedError


# After applying both filters (Steiger + Cook's > 4/N from `egger` above),
# recompute IVW, MR-Egger, weighted median. Display the comparison.
'''


STEP5_SOLUTION = '''# Reference solution -- Step 5.

def steiger_filter(df: pd.DataFrame) -> np.ndarray:
    maf = df["maf"].to_numpy()
    bx  = df["beta_zx"].to_numpy()
    by  = df["beta_zy"].to_numpy()
    r2_zx = 2.0 * maf * (1.0 - maf) * bx ** 2
    r2_zy = 2.0 * maf * (1.0 - maf) * by ** 2
    return r2_zx > r2_zy


steiger_ok = steiger_filter(df)
cook_ok    = egger["cook"] <= 4.0 / len(df)
keep       = steiger_ok & cook_ok
print(f"Steiger pass:       {steiger_ok.sum():3d} / {len(df)}")
print(f"Cook's-distance ok: {cook_ok.sum():3d} / {len(df)}")
print(f"Both filters:       {keep.sum():3d} / {len(df)}  (dropped {(~keep).sum()})")

df_clean = df[keep].reset_index(drop=True)

ivw_c   = ivw_estimate(df_clean)
egger_c = mr_egger(df_clean)
wmed_c  = weighted_median_mr(df_clean, n_boot=1000)

results_clean = pd.DataFrame([
    {"method": "IVW",             "beta": ivw_c["beta"],   "se": ivw_c["se"],   "ci_low": ivw_c["ci_low"],   "ci_high": ivw_c["ci_high"]},
    {"method": "MR-Egger",        "beta": egger_c["beta"], "se": egger_c["se"], "ci_low": egger_c["ci_low"], "ci_high": egger_c["ci_high"]},
    {"method": "Weighted median", "beta": wmed_c["beta"],  "se": wmed_c["se"],  "ci_low": wmed_c["ci_low"],  "ci_high": wmed_c["ci_high"]},
])
print("\\nPost-filter triangulation:")
print(results_clean.to_string(index=False, float_format=lambda x: f"{x:6.3f}"))

# Before / after forest plot.
fig, ax = plt.subplots(figsize=(8, 4.2))
methods = results["method"].tolist()
y_before = np.arange(len(methods)) + 0.18
y_after  = np.arange(len(methods)) - 0.18
ax.errorbar(results["beta"], y_before,
            xerr=[results["beta"] - results["ci_low"], results["ci_high"] - results["beta"]],
            fmt="o", ms=7, color="#d62728", capsize=4, label="raw (100 SNPs)")
ax.errorbar(results_clean["beta"], y_after,
            xerr=[results_clean["beta"] - results_clean["ci_low"], results_clean["ci_high"] - results_clean["beta"]],
            fmt="s", ms=7, color="#3367d6", capsize=4,
            label=f"post-filter ({keep.sum()} SNPs)")
ax.axvline(TRUE_BETA, color="black", lw=1.2, ls="--", label=f"truth = {TRUE_BETA}")
ax.axvline(0, color="grey", lw=0.5)
ax.set_yticks(np.arange(len(methods)))
ax.set_yticklabels(methods)
ax.set_xlabel(r"causal effect estimate $\\hat\\beta$")
ax.set_title("Sensitivity triangulation: raw vs Steiger+Cook-filtered")
ax.legend()
plt.tight_layout()
plt.show()

# Spread of the three estimators -- agreement is what we want.
spread_before = results["beta"].max() - results["beta"].min()
spread_after  = results_clean["beta"].max() - results_clean["beta"].min()
print(f"\\nEstimator spread (max - min):  before = {spread_before:.3f}  after = {spread_after:.3f}")
print("Tighter spread + all three CIs covering TRUE_BETA = triangulated causal claim.")
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing pieces of the pipeline against
the known ground truth. If a check fails, walk back through the
corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check -- ties every estimator back to TRUE_BETA = 0.5.
# ----------------------------------------------------------------------

# 1. Simulator returns the right shape and the F-statistic is comfortably > 10.
assert df is not None and len(df) == N_SNPS, "df missing or wrong size"
F_mean = ((df["beta_zx"] / df["se_zx"]) ** 2).mean()
assert F_mean > 10.0, f"Weak instruments: mean F = {F_mean:.1f}"

# 2. IVW is biased upwards by directional pleiotropy but the CI is at least
#    in the right neighbourhood.
assert 0.30 <= ivw["beta"] <= 0.80, f"IVW beta out of plausible band: {ivw['beta']:.3f}"

# 3. MR-Egger intercept is significantly positive (we injected mean alpha = 0.1).
assert egger["intercept"] > 0.0, f"Egger intercept not positive: {egger['intercept']:.4f}"
assert egger["intercept_p"] < 0.05, (
    f"Egger intercept test p={egger['intercept_p']:.2e} -- expected significant "
    "given simulated directional pleiotropy"
)

# 4. Weighted median is consistent: should recover TRUE_BETA within tolerance.
assert abs(wmed["beta"] - TRUE_BETA) < 0.15, (
    f"Weighted median far from truth: {wmed['beta']:.3f} vs {TRUE_BETA}"
)

# 5. Cook's distance flags at least *some* of the truly pleiotropic SNPs.
threshold = 4.0 / len(df)
flagged = np.where(egger["cook"] > threshold)[0]
pleio_idx = np.where(~df["valid"].to_numpy())[0]
hits = np.intersect1d(flagged, pleio_idx)
assert len(hits) >= 3, (
    f"Cook's distance recovered only {len(hits)} / {len(pleio_idx)} pleiotropic SNPs"
)

# 6. After Steiger + Cook filtering, *all three* estimators recover TRUE_BETA
#    within a generous tolerance, AND the IVW estimator improves.
TOL = 0.20
for r in results_clean.itertuples():
    assert abs(r.beta - TRUE_BETA) < TOL, (
        f"{r.method} post-filter beta = {r.beta:.3f} too far from {TRUE_BETA}"
    )
assert abs(ivw_c["beta"] - TRUE_BETA) < abs(ivw["beta"] - TRUE_BETA), (
    "Filtering did not reduce IVW bias"
)

print()
print(f"Recovered TRUE_BETA = {TRUE_BETA} via triangulation:")
print(f"  IVW             post-filter: {ivw_c['beta']:.3f}  (raw: {ivw['beta']:.3f})")
print(f"  MR-Egger        post-filter: {egger_c['beta']:.3f}  (raw: {egger['beta']:.3f})")
print(f"  Weighted median post-filter: {wmed_c['beta']:.3f}  (raw: {wmed['beta']:.3f})")
print()
print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — instrumental variables = orthogonal regression

You implemented the canonical two-sample MR pipeline:

1. **Two-stage least squares as orthogonal regression.** The Wald scatter
   places $\\hat\\beta_{ZX}$ and $\\hat\\beta_{ZY}$ on the two axes. Under
   the three IV assumptions the slope **through the origin** is the
   causal effect — exactly the **errors-in-variables** regression an EE
   would fit when both axes are noisy measurements of underlying truth.
   IVW is just weighted least squares with the intercept clamped to 0;
   MR-Egger lifts that constraint to absorb a bias term.
2. **L₁ vs L₂ robustness.** IVW minimises a sum of squared residuals →
   one bad SNP can pull the line. The weighted median minimises a sum
   of weighted absolute deviations → up to 50 % of weight can be
   corrupted before the estimator breaks down. Same story as median vs
   mean, or L₁ vs L₂ filter design: trading some efficiency under the
   null for robustness to outliers.
3. **Cook's distance = leverage × residual.** A single observation can
   move a least-squares fit if both (a) its residual is large and
   (b) it sits at high leverage (far from the centroid of $X$). Cook's
   distance multiplies the two. Same idea as the hat matrix in adaptive
   beamforming — knowing which sensor has the loudest say in the
   output.
4. **Triangulation = sensor fusion.** Each MR estimator makes a
   *different* assumption about the noise. When their estimates agree
   you have a consistent reading across several independent assumption
   sets — the same logic as combining radar + lidar + camera estimates
   of position. When they disagree, somebody is wrong, and the
   disagreement tells you in which direction the bias points.
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
