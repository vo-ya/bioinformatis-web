"""Build exercise.ipynb for L13 — GWAS and Statistical Genetics.

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


TITLE_MD = """# L13 — GWAS and Statistical Genetics

In this exercise you simulate a Genome-Wide Association Study from scratch:
a 10 000-individual × 50 000-SNP cohort with 100 causal SNPs **and** a hidden
ancestry axis that is mildly correlated (r = 0.15) with the phenotype. You
will run 50 000 parallel per-SNP regressions, build Manhattan and QQ plots,
compute the genomic inflation factor `lambda_GC`, then re-run the analysis
after regressing out the ancestry axis and watch the QQ tail collapse.

The whole pipeline runs in well under 5 minutes on free Colab CPU.
"""


AHA_MD = """> **Aha.** Population **stratification** lifts the *whole* QQ tail —
> every test statistic gets inflated because every SNP's allele frequency
> co-varies with ancestry, which co-varies with phenotype. **Polygenicity**
> lifts only the *upper* QQ tail — most SNPs are still null and follow the
> diagonal; only the causal ~100 escape. The *shape* of the QQ plot tells
> the two apart; `lambda_GC` alone cannot.
"""


PREAMBLE = """# Install the pinned scientific stack. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 matplotlib==3.8.4 -q
"""


IMPORTS = """import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Deterministic for the whole notebook.
RNG_SEED = 13
rng = np.random.default_rng(RNG_SEED)

# Cohort and genotype-matrix dimensions.
N_INDIV   = 10_000   # individuals
N_SNPS    = 50_000   # variants tested
N_CAUSAL  = 100      # SNPs with true non-zero effect
ANC_RHO   = 0.15     # correlation between latent ancestry and phenotype

# Genome-wide significance threshold (Bonferroni-style; Manhattan reference line).
GW_SIG    = 5e-8

print(f"Simulation knobs:")
print(f"  N_INDIV  = {N_INDIV:,}")
print(f"  N_SNPS   = {N_SNPS:,}")
print(f"  N_CAUSAL = {N_CAUSAL}")
print(f"  ANC_RHO  = {ANC_RHO}")
print(f"  GW_SIG   = {GW_SIG:.0e}  (Manhattan reference line)")
"""


STEP1_MD = """## Step 1 (10 min) — Simulate the cohort

We need three ingredients:

1. **A latent ancestry axis** `a ~ N(0, 1)` per individual. Think of it as the
   leading PC of the genotype matrix in a real cohort: a continuous coordinate
   along which allele frequencies drift.
2. **A genotype matrix** `G ∈ {0,1,2}^{N×M}`. Most SNPs are *neutral* — their
   minor-allele frequency does not depend on ancestry. A subset of **stratified**
   SNPs has an ancestry-dependent MAF; this is the structure that creates
   spurious associations.
3. **A binary phenotype** `y` whose log-odds depend on (a) the 100 causal SNP
   genotypes with small per-SNP effects, plus (b) a contribution from the
   ancestry axis itself — a non-genetic confounder. We tune the ancestry weight
   so that `corr(a, y) ≈ 0.15`.

Memory note: 10 000 × 50 000 in `int8` is ~500 MB. We use `uint8` and avoid
copies; the whole matrix fits comfortably in Colab RAM.
"""

STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Simulate ancestry, genotypes, causal effects, phenotype.
# ----------------------------------------------------------------------

def simulate_cohort(
    n_indiv: int = N_INDIV,
    n_snps: int = N_SNPS,
    n_causal: int = N_CAUSAL,
    target_anc_rho: float = ANC_RHO,
    seed: int = RNG_SEED,
):
    """Return a dict with keys:
        ancestry        (n_indiv,) float
        genotypes       (n_indiv, n_snps) uint8 in {0,1,2}
        causal_idx      (n_causal,) int      — indices of causal SNPs
        causal_beta     (n_causal,) float    — true effects (log-odds)
        phenotype       (n_indiv,) int       — 0/1 binary trait
        is_stratified   (n_snps,) bool       — SNPs with ancestry-dependent MAF
    """
    # TODO:
    #   1. Draw ancestry ~ N(0,1).
    #   2. For each SNP, draw a base MAF ~ Uniform(0.05, 0.5).
    #      Mark roughly half of all SNPs as "stratified": their per-individual
    #      MAF shifts linearly with ancestry (clamped to [0.02, 0.98]).
    #   3. Sample genotypes Binomial(2, p_individual) per (SNP, individual).
    #   4. Pick `n_causal` SNP indices; draw small effects (e.g. N(0, 0.1)).
    #      Compute the genetic log-odds contribution.
    #   5. Add an ancestry contribution to the log-odds; tune its scale so the
    #      observed corr(ancestry, phenotype) is close to target_anc_rho.
    #   6. Pass log-odds through sigmoid -> Bernoulli sample for the phenotype.
    raise NotImplementedError


cohort = None  # cohort = simulate_cohort()
'''

STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def simulate_cohort(
    n_indiv: int = N_INDIV,
    n_snps: int = N_SNPS,
    n_causal: int = N_CAUSAL,
    target_anc_rho: float = ANC_RHO,
    seed: int = RNG_SEED,
):
    rng = np.random.default_rng(seed)

    # 1. Latent ancestry axis, standardised.
    ancestry = rng.standard_normal(n_indiv).astype(np.float32)

    # 2. Per-SNP MAFs. ~50% of SNPs are stratified (ancestry-dependent freq).
    base_maf = rng.uniform(0.05, 0.5, size=n_snps).astype(np.float32)
    is_stratified = rng.random(n_snps) < 0.5
    # Stratification slope: each stratified SNP has its own random allele-frequency
    # gradient w.r.t. ancestry. The 0.10 scale is enough to cause inflation but
    # not so much that frequencies routinely hit the bounds.
    strat_slope = np.where(
        is_stratified,
        rng.normal(0.0, 0.10, size=n_snps).astype(np.float32),
        0.0,
    )

    # 3. Sample genotypes one SNP at a time so we never materialise a full
    #    (N x M) float matrix of probabilities.
    genotypes = np.empty((n_indiv, n_snps), dtype=np.uint8)
    for j in range(n_snps):
        p = base_maf[j] + strat_slope[j] * ancestry
        np.clip(p, 0.02, 0.98, out=p)
        genotypes[:, j] = rng.binomial(2, p).astype(np.uint8)

    # 4. Causal SNPs: small additive log-odds effects.
    causal_idx = rng.choice(n_snps, size=n_causal, replace=False)
    causal_idx.sort()
    causal_beta = rng.normal(0.0, 0.10, size=n_causal).astype(np.float32)
    # Standardise causal genotypes column-wise before applying betas so the
    # phenotypic variance contributed by causal SNPs is in a sensible range.
    G_caus = genotypes[:, causal_idx].astype(np.float32)
    G_caus -= G_caus.mean(axis=0, keepdims=True)
    sd = G_caus.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    G_caus /= sd
    genetic_logodds = G_caus @ causal_beta

    # 5. Ancestry contribution to the log-odds; pick its weight so that
    #    corr(ancestry, phenotype) ~= target_anc_rho.
    # The target correlation between a continuous predictor and a Bernoulli
    # outcome whose logit is alpha + w * ancestry is roughly w / sqrt(w^2 + c),
    # for some constant c that depends on the link. We solve for w empirically
    # at a coarse grid and pick the closest match.
    best_w, best_diff = 0.0, np.inf
    for w in np.linspace(0.05, 1.5, 30):
        logits = -0.2 + genetic_logodds + w * ancestry
        p_y = _sigmoid(logits)
        y_try = (rng.random(n_indiv) < p_y).astype(np.int8)
        r = float(np.corrcoef(ancestry, y_try)[0, 1])
        diff = abs(r - target_anc_rho)
        if diff < best_diff:
            best_diff, best_w = diff, w
    anc_weight = float(best_w)

    # 6. Generate the phenotype with the chosen ancestry weight.
    logits = -0.2 + genetic_logodds + anc_weight * ancestry
    phenotype = (rng.random(n_indiv) < _sigmoid(logits)).astype(np.int8)

    obs_rho = float(np.corrcoef(ancestry, phenotype)[0, 1])
    print(f"Ancestry weight (chosen) = {anc_weight:.3f}")
    print(f"Observed corr(ancestry, phenotype) = {obs_rho:+.3f}  (target {target_anc_rho:+.2f})")
    print(f"Case rate = {phenotype.mean():.3f}")

    return {
        "ancestry":      ancestry,
        "genotypes":     genotypes,
        "causal_idx":    causal_idx,
        "causal_beta":   causal_beta,
        "phenotype":     phenotype,
        "is_stratified": is_stratified,
        "anc_weight":    anc_weight,
    }


t0 = time.time()
cohort = simulate_cohort()
print(f"Simulated cohort in {time.time() - t0:.1f}s")
print(f"  genotypes:  shape={cohort['genotypes'].shape}, dtype={cohort['genotypes'].dtype}, "
      f"size={cohort['genotypes'].nbytes / 1e6:.1f} MB")
print(f"  causal SNPs: {len(cohort['causal_idx'])}  "
      f"(stratified among them: {int(cohort['is_stratified'][cohort['causal_idx']].sum())})")
'''


STEP2_MD = """## Step 2 (12 min) — Per-SNP association test (vectorised)

The textbook test for a binary phenotype is a per-SNP logistic regression. For
a 50 000-SNP scan we approximate it with a much faster **correlation z-test**
that gives essentially the same z-scores when effects are small (which they
always are in real GWAS):

$$r_j = \\frac{\\mathrm{cov}(G_j, y)}{\\sqrt{\\mathrm{var}(G_j)\\cdot\\mathrm{var}(y)}}\\qquad z_j = r_j\\sqrt{N-2}$$

This is one matrix-vector multiply against the column-centred genotype matrix.
The two-sided p-value is `2 * Φ(-|z|)`. We will compute z-scores and p-values
for every SNP **with stratification present** — the noisy baseline. Step 4
will repeat the scan after whitening out the ancestry axis.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Vectorised per-SNP association scan.
# ----------------------------------------------------------------------

def run_scan(
    genotypes: np.ndarray,
    y: np.ndarray,
    covariates: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (z_scores, p_values), one per SNP.

    Parameters
    ----------
    genotypes : (N, M) uint8 in {0,1,2}
    y         : (N,)  int   in {0,1}
    covariates: (N, K) float, optional. If given, regress each genotype column
        AND y on the covariates and run the scan on the residuals (Step 4).
    """
    # TODO:
    #   1. Cast genotypes to float32 (batch by columns if needed to fit memory).
    #   2. Residualise on covariates if provided (project them out of G and y).
    #   3. Compute centred y and centred G per column.
    #   4. cov_j = (G_j_centered . y_centered) / (N-1)
    #      var_j = var(G_j); var_y = var(y)
    #   5. r_j = cov_j / sqrt(var_j * var_y);  z_j = r_j * sqrt(N - 2)
    #   6. p_j = 2 * survival_function(|z_j|) for the standard normal.
    raise NotImplementedError


z_unc, p_unc = (np.zeros(N_SNPS), np.zeros(N_SNPS))   # uncorrected scan
'''

STEP2_SOLUTION = '''# Reference solution — Step 2.

def _residualise(X: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Regress columns of X on C (with intercept absorbed into C if you want).
    Returns the residual; works on float32 to keep memory down.
    """
    # Add an intercept column so we also remove the mean.
    C1 = np.column_stack([np.ones(C.shape[0], dtype=np.float32), C.astype(np.float32)])
    # Beta = (C1^T C1)^-1 C1^T X
    G_inv = np.linalg.pinv(C1.T @ C1)
    beta = G_inv @ (C1.T @ X)
    return X - C1 @ beta


def run_scan(
    genotypes: np.ndarray,
    y: np.ndarray,
    covariates: np.ndarray | None = None,
    batch: int = 5000,
) -> tuple[np.ndarray, np.ndarray]:
    n = genotypes.shape[0]
    m = genotypes.shape[1]
    y = y.astype(np.float32)

    if covariates is not None:
        # Residualise y once.
        y_res = _residualise(y.reshape(-1, 1), covariates).ravel()
    else:
        y_res = y - y.mean()

    var_y = float((y_res * y_res).sum() / (n - 1))

    z_out = np.empty(m, dtype=np.float32)

    # Process columns in batches so we never need an (N, M) float32 matrix in RAM.
    for s in range(0, m, batch):
        e = min(s + batch, m)
        G = genotypes[:, s:e].astype(np.float32)
        if covariates is not None:
            G = _residualise(G, covariates)
        else:
            G = G - G.mean(axis=0, keepdims=True)

        # Per-column variance and covariance with residualised y.
        var_G = (G * G).sum(axis=0) / (n - 1)
        cov_Gy = (G * y_res[:, None]).sum(axis=0) / (n - 1)
        # Pearson correlation per SNP, then z = r * sqrt(N-2) (large-N asymptote of the t-test).
        denom = np.sqrt(var_G * var_y)
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(denom > 0, cov_Gy / np.where(denom > 0, denom, 1.0), 0.0)
        z = r * np.sqrt(n - 2)
        z_out[s:e] = z.astype(np.float32)

    p_out = 2.0 * stats.norm.sf(np.abs(z_out))
    # Clip away exact 0 / 1 so log transforms remain finite.
    p_out = np.clip(p_out, 1e-300, 1.0)
    return z_out, p_out


t0 = time.time()
z_unc, p_unc = run_scan(cohort["genotypes"], cohort["phenotype"])
print(f"Uncorrected scan: {len(z_unc):,} SNPs in {time.time() - t0:.1f}s")
print(f"  median |z| = {np.median(np.abs(z_unc)):.3f}")
print(f"  #SNPs with p < 5e-8: {int((p_unc < GW_SIG).sum())}")
print(f"  #true-causal with p < 5e-8: "
      f"{int((p_unc[cohort['causal_idx']] < GW_SIG).sum())} / {N_CAUSAL}")
'''


STEP3_MD = """## Step 3 (15 min) — Manhattan plot and the 5e-8 threshold

The canonical GWAS visualisation: one dot per SNP, x-axis = SNP index (proxy
for genomic position), y-axis = `-log10(p)`. The horizontal line at
`-log10(5e-8) ≈ 7.30` is the genome-wide significance threshold (Bonferroni at
α = 0.05 over ~10⁶ effectively independent tests).

Highlight the true-causal SNPs in a different colour so we can see how many
escape the threshold versus how many spurious hits the stratification creates.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Manhattan plot.
# ----------------------------------------------------------------------

def manhattan(
    p_values: np.ndarray,
    causal_idx: np.ndarray,
    title: str = "GWAS — Manhattan plot",
    threshold: float = GW_SIG,
):
    """Render a Manhattan plot with the genome-wide significance line and the
    true-causal SNPs highlighted.
    """
    # TODO: scatter all SNPs faintly; overlay causal SNPs in a bold colour;
    # add a horizontal line at -log10(threshold).
    raise NotImplementedError


# manhattan(p_unc, cohort["causal_idx"], title="Uncorrected GWAS — stratification present")
'''

STEP3_SOLUTION = '''# Reference solution — Step 3.

def manhattan(
    p_values: np.ndarray,
    causal_idx: np.ndarray,
    title: str = "GWAS — Manhattan plot",
    threshold: float = GW_SIG,
):
    neglog = -np.log10(p_values)
    x = np.arange(len(p_values))

    fig, ax = plt.subplots(figsize=(10, 3.5))
    # Faint cloud for all SNPs.
    ax.scatter(x, neglog, s=2, c="#888", alpha=0.4, rasterized=True)
    # Bold dots for true-causal SNPs.
    ax.scatter(causal_idx, neglog[causal_idx], s=14, c="#d62728",
               label="true-causal SNPs", zorder=3)
    # Threshold.
    ax.axhline(-np.log10(threshold), color="#1f77b4", linestyle="--",
               label=f"genome-wide sig. ({threshold:.0e})")
    ax.set_xlabel("SNP index (proxy for genomic position)")
    ax.set_ylabel(r"$-\\log_{10}(p)$")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.show()


manhattan(p_unc, cohort["causal_idx"],
          title="Uncorrected GWAS — stratification present")

n_sig_unc       = int((p_unc < GW_SIG).sum())
n_sig_causal    = int((p_unc[cohort["causal_idx"]] < GW_SIG).sum())
n_sig_spurious  = n_sig_unc - n_sig_causal
print(f"Genome-wide significant SNPs (uncorrected): {n_sig_unc:,}")
print(f"  of which true-causal: {n_sig_causal}  ({100*n_sig_causal/max(1,n_sig_unc):.1f}%)")
print(f"  of which spurious   : {n_sig_spurious:,}")
'''


STEP4_MD = """## Step 4 (15 min) — QQ plot, `lambda_GC`, and the whitened re-scan

Two diagnostics together tell you whether your GWAS is calibrated:

1. **QQ plot.** Sort the observed `-log10(p)` ascending; plot against the
   theoretical quantiles of a uniform-on-(0,1) p-value distribution. Under a
   clean null the points hug the diagonal until the upper tail, where genuine
   hits depart. **Stratification** lifts the *whole* tail above the diagonal;
   **polygenicity** lifts only the *upper* tail.
2. **`lambda_GC`** — the genomic inflation factor. Take the median χ² statistic
   across SNPs and divide by the median of a χ²₁ (= 0.456). `lambda_GC = 1.0`
   means calibrated; `> 1.05` is suspicious; `> 1.2` is a red flag.

You will compute both for the uncorrected scan, then re-run the scan with the
ancestry axis as a covariate (this is the EE story: **PCA correction = noise
whitening**), and overlay the two QQs. The whitened QQ should drop onto the
diagonal until the very top, where the 100 causal SNPs poke out.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — QQ plot, lambda_GC, whitened re-scan.
# ----------------------------------------------------------------------

def lambda_gc(p_values: np.ndarray) -> float:
    """Genomic inflation factor = median(chi^2) / 0.4549.
    chi^2 = norm.ppf(1 - p/2)**2 for a 1-df test.
    """
    # TODO
    raise NotImplementedError


def qq_plot(
    p_values_by_label: dict,
    title: str = "QQ plot",
):
    """For each {label: p_array} entry, plot expected vs observed -log10(p)."""
    # TODO: sort each p_array, compute expected quantiles
    #   exp = -log10((rank + 0.5) / N), plot exp vs obs; overlay y=x reference.
    raise NotImplementedError


# 1. Compute lambda_GC for the uncorrected scan.
# 2. Re-run the scan including cohort["ancestry"] as a single-column covariate.
# 3. Compute lambda_GC for the corrected scan.
# 4. Overlay the two QQs and the Manhattan after correction.
z_corr, p_corr = (np.zeros(N_SNPS), np.zeros(N_SNPS))
'''

STEP4_SOLUTION = '''# Reference solution — Step 4.

CHI2_1_MEDIAN = stats.chi2.ppf(0.5, df=1)  # ~0.4549


def lambda_gc(p_values: np.ndarray) -> float:
    # Convert p to chi^2 (1 df) via the inverse normal CDF.
    z = stats.norm.isf(p_values / 2.0)
    chi2 = z * z
    return float(np.median(chi2) / CHI2_1_MEDIAN)


def qq_plot(p_values_by_label: dict, title: str = "QQ plot"):
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    colours = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd"]
    max_y = 0.0
    for (label, pvals), col in zip(p_values_by_label.items(), colours):
        n = len(pvals)
        obs = -np.log10(np.sort(pvals))
        exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
        # Downsample for plotting speed; keep the upper tail intact.
        if n > 20000:
            tail = 5000
            stride = max(1, (n - tail) // 15000)
            idx_body = np.arange(0, n - tail, stride)
            idx_tail = np.arange(n - tail, n)
            idx = np.concatenate([idx_body, idx_tail])
        else:
            idx = np.arange(n)
        ax.scatter(exp[idx], obs[idx], s=4, alpha=0.6, c=col,
                   label=f"{label}  (λ_GC = {lambda_gc(pvals):.3f})")
        max_y = max(max_y, float(obs[idx].max()))
    lim = max(max_y, 8.0)
    ax.plot([0, lim], [0, lim], "k-", lw=1, label="y = x (null)")
    ax.set_xlabel(r"expected $-\\log_{10}(p)$")
    ax.set_ylabel(r"observed $-\\log_{10}(p)$")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.show()


# 1. lambda_GC for the uncorrected scan.
lam_unc = lambda_gc(p_unc)

# 2. Re-run the scan with ancestry as a covariate (PCA correction = noise whitening).
covariates = cohort["ancestry"].reshape(-1, 1)
t0 = time.time()
z_corr, p_corr = run_scan(cohort["genotypes"], cohort["phenotype"],
                          covariates=covariates)
print(f"Corrected scan (ancestry covariate): {time.time() - t0:.1f}s")

lam_corr = lambda_gc(p_corr)

print(f"\\nlambda_GC")
print(f"  uncorrected (stratification + polygenicity): {lam_unc:.3f}")
print(f"  corrected   (ancestry projected out)        : {lam_corr:.3f}")

# 3. Overlay both QQs.
qq_plot(
    {
        "uncorrected": p_unc,
        "corrected (ancestry covariate)": p_corr,
    },
    title="QQ plot — stratification lifts the whole tail",
)

# 4. Manhattan after correction.
manhattan(p_corr, cohort["causal_idx"],
          title="Corrected GWAS — ancestry regressed out")

n_sig_corr = int((p_corr < GW_SIG).sum())
n_sig_causal_corr = int((p_corr[cohort["causal_idx"]] < GW_SIG).sum())
print(f"\\nGenome-wide significant SNPs (corrected): {n_sig_corr:,}")
print(f"  of which true-causal: {n_sig_causal_corr}  "
      f"({100*n_sig_causal_corr/max(1,n_sig_corr):.1f}%)")
'''


STEP5_MD = """## Step 5 (8 min) — Top-10 SNP table and the final λ summary

The deliverable from a GWAS scan is a sorted summary-statistics table. Build a
DataFrame with the top-10 SNPs by p-value from each scan (uncorrected and
corrected), flagging whether each was actually causal in the simulation. Print
the two λ_GC values side-by-side as the final EE-style calibration report.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Top-10 SNP tables and side-by-side lambda_GC.
# ----------------------------------------------------------------------

def top_n_table(z: np.ndarray, p: np.ndarray, causal_idx: np.ndarray,
                is_strat: np.ndarray, n: int = 10) -> pd.DataFrame:
    """Top-n SNPs by p-value with causal / stratified flags."""
    # TODO: order = np.argsort(p)[:n]; build a DataFrame.
    raise NotImplementedError


# top_unc  = top_n_table(z_unc,  p_unc,  cohort["causal_idx"], cohort["is_stratified"])
# top_corr = top_n_table(z_corr, p_corr, cohort["causal_idx"], cohort["is_stratified"])
'''

STEP5_SOLUTION = '''# Reference solution — Step 5.

def top_n_table(z: np.ndarray, p: np.ndarray, causal_idx: np.ndarray,
                is_strat: np.ndarray, n: int = 10) -> pd.DataFrame:
    order = np.argsort(p)[:n]
    causal_set = set(int(i) for i in causal_idx)
    rows = []
    for snp in order:
        rows.append({
            "snp_idx":      int(snp),
            "z":            float(z[snp]),
            "p":            float(p[snp]),
            "neglog10_p":   float(-np.log10(max(p[snp], 1e-300))),
            "causal":       int(snp) in causal_set,
            "stratified":   bool(is_strat[snp]),
        })
    return pd.DataFrame(rows)


top_unc  = top_n_table(z_unc,  p_unc,  cohort["causal_idx"], cohort["is_stratified"])
top_corr = top_n_table(z_corr, p_corr, cohort["causal_idx"], cohort["is_stratified"])

print("Top-10 SNPs — uncorrected scan (stratification present):")
print(top_unc.to_string(index=False))
print()
print("Top-10 SNPs — corrected scan (ancestry covariate):")
print(top_corr.to_string(index=False))
print()
print("lambda_GC summary")
print(f"  uncorrected: {lambda_gc(p_unc):.3f}")
print(f"  corrected  : {lambda_gc(p_corr):.3f}")
print()
print("Interpretation:")
print("  Uncorrected: the whole QQ tail is lifted; many top hits are spurious")
print("    (stratified=True but causal=False). lambda_GC well above 1.0.")
print("  Corrected:  ancestry projected out (= noise whitening). QQ drops to the")
print("    diagonal except in the upper tail where the real ~100 causal SNPs poke")
print("    through. lambda_GC near 1 with a (slight) upward bump from polygenicity.")
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces of the pipeline. If
you ran the reference solutions above they should all pass; if you wrote your
own and an assert fails, revisit the corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Cohort shape is exactly what we asked for.
G = cohort["genotypes"]
assert G.shape == (N_INDIV, N_SNPS), f"Bad genotype shape: {G.shape}"
assert G.dtype == np.uint8, f"Genotypes should be uint8, got {G.dtype}"
assert set(np.unique(G[:100, :100]).tolist()).issubset({0, 1, 2}), \\
    "Genotypes must lie in {0, 1, 2}"

# 2. Ancestry-phenotype correlation is in the intended ballpark (target = 0.15).
obs_rho = float(np.corrcoef(cohort["ancestry"], cohort["phenotype"])[0, 1])
assert 0.05 < obs_rho < 0.25, \\
    f"corr(ancestry, phenotype) = {obs_rho:.3f}, expected ~0.15"

# 3. Scan returned a z and p per SNP and the p-values are valid.
assert z_unc.shape == (N_SNPS,) and p_unc.shape == (N_SNPS,)
assert np.all((p_unc > 0) & (p_unc <= 1.0)), "p-values out of (0, 1]"

# 4. lambda_GC moves the right way: uncorrected is meaningfully > 1, corrected drops.
lam_u = lambda_gc(p_unc)
lam_c = lambda_gc(p_corr)
assert lam_u > 1.05, \\
    f"Uncorrected lambda_GC = {lam_u:.3f} — stratification should lift it above 1.05"
assert lam_c < lam_u, \\
    f"Correction did not reduce lambda_GC: {lam_u:.3f} -> {lam_c:.3f}"
assert lam_c < 1.20, \\
    f"Corrected lambda_GC = {lam_c:.3f} — should drop below 1.20 after whitening"

# 5. Causal SNPs are over-represented in the corrected top-N — the real signal
#    survives whitening.
n_top = 200
top_corr_idx = np.argsort(p_corr)[:n_top]
causal_set = set(int(i) for i in cohort["causal_idx"])
hits_top = sum(int(i) in causal_set for i in top_corr_idx)
# Random expectation = n_top * N_CAUSAL / N_SNPS = 200 * 100/50000 = 0.4
assert hits_top >= 5, \\
    f"Only {hits_top} of top-{n_top} corrected SNPs are causal — signal is missing"

print(f"lambda_GC uncorrected = {lam_u:.3f}")
print(f"lambda_GC corrected   = {lam_c:.3f}")
print(f"causal SNPs in corrected top-{n_top}: {hits_top}  (random ≈ 0.4)")
print()
print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — multi-channel detection and noise whitening

You ran **50 000 parallel single-channel hypothesis tests**. Each SNP is one
detector with its own z-score; the genome is the channel axis. Two clean
analogies fall out:

1. **Genome-wide significance = controlled false-alarm rate at scale.**
   At M = 5 × 10⁴ tests we used the Bonferroni-style threshold 5 × 10⁻⁸ that
   modern GWAS uses across ~10⁶ effectively independent SNPs. This is the
   same multi-hypothesis problem as multi-band radar (one detector per
   frequency cell), multi-target tracking (one detector per range gate), or
   per-pixel anomaly scoring in imaging.

2. **PCA correction = noise whitening.** The ancestry axis induces a
   *correlated* component in every SNP's null distribution — its test
   statistic depends on the same hidden variable as the phenotype. Projecting
   that direction out of the genotype matrix is exactly what an adaptive
   array does to suppress jamming or what a Wiener filter does to convert
   coloured noise into white noise. Once the noise is white, per-SNP z
   really is N(0, 1) under the null, and the QQ plot snaps back to the
   diagonal.

The QQ plot itself is the **inverse-CDF probability transform** — it asks
"do my observed p-values match a Uniform(0, 1) reference?" — exactly the
same diagnostic an EE student would run on a noise model.

If you want to keep going, the natural next step is to compute the
**leading principal components of the genotype matrix** directly (instead
of using the known ancestry axis) and use those as covariates. That is the
PCA-correction recipe that became the field standard in Price et al. 2006.
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
