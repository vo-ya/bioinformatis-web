"""Build exercise.ipynb for L24 — CRISPR Functional Screens and DepMap.

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


TITLE_MD = """# L24 — CRISPR Functional Screens and DepMap

In this exercise you implement the core of **MAGeCK**: a per-sgRNA
negative-binomial test followed by a rank-based gene-level aggregation
(α-RRA). You will analyse a simulated pooled CRISPR knockout screen —
10 000 sgRNAs × 6 samples (3 control + 3 treatment) — recover the 50
planted essential genes and 20 resistance genes from 430 distractors,
and inspect the false-discovery behaviour.
"""


AHA_MD = """> **Aha.** A pooled CRISPR screen is **compressed sensing**: a sparse
> signal (≲ 70 hit genes out of 500) is recovered from 10 000 noisy
> count measurements. The negative-binomial test is a **matched filter**
> in the count domain (overdispersion-aware). Robust rank aggregation
> turns 20 noisy per-sgRNA p-values per gene into a single gene-level
> call — built-in replication. Multiple sgRNAs per gene are the
> repetition code that makes the inverse problem well-posed.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 matplotlib==3.8.4 -q
"""


IMPORTS = """import math
import time
from dataclasses import dataclass

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Deterministic for the whole notebook.
RNG = np.random.default_rng(42)
np.random.seed(42)

# Library and experiment shape.
N_SGRNAS_PER_GENE = 20      # CRISPR screens use 4-6; we use 20 to make rank-aggregation behaviour vivid
N_GENES           = 500     # 50 essential + 20 resistant + 430 nulls
N_ESSENTIAL       = 50
N_RESISTANT       = 20
N_NULL            = N_GENES - N_ESSENTIAL - N_RESISTANT
N_SGRNAS          = N_GENES * N_SGRNAS_PER_GENE   # 10 000
N_CONTROL         = 3
N_TREATMENT       = 3
N_SAMPLES         = N_CONTROL + N_TREATMENT

# Truth: which genes have which planted effect (in log2-fold-change units).
LFC_ESSENTIAL = -2.0    # depletion under treatment (~4x dropout)
LFC_RESISTANT = +1.0    # enrichment under treatment (~2x gain)

print(f"Screen shape: {N_SGRNAS} sgRNAs x {N_SAMPLES} samples")
print(f"  {N_ESSENTIAL} essential genes (planted LFC = {LFC_ESSENTIAL})")
print(f"  {N_RESISTANT} resistant genes (planted LFC = {LFC_RESISTANT})")
print(f"  {N_NULL} null genes (planted LFC = 0)")
print(f"  {N_SGRNAS_PER_GENE} sgRNAs per gene")
"""


STEP1_MD = """## Step 1 (8 min) — Generate the synthetic count matrix

Real screens count sgRNA amplicons via NGS: counts are integer, **overdispersed**
(Var > Mean), and have **size-factor** differences between samples (sequencing
depth varies). We simulate a negative-binomial count for each sgRNA × sample
cell, with a planted log-fold-change per gene that applies only to the
treatment samples.

The NB(μ, α) parametrisation we use:
$$
\\mathrm{Var}(X) = \\mu + \\alpha \\mu^2
$$
α is the dispersion. SciPy's `nbinom(n, p)` parameters relate via
`n = 1/α`, `p = n / (n + μ)`.

You will (a) draw baseline mean counts from a realistic log-normal sgRNA
abundance distribution, (b) apply per-sample size factors (≈ 0.8–1.25), (c)
multiply treatment means by `2**LFC` for the planted genes, (d) draw NB
counts. You will also write a **median-of-ratios** size-factor estimator
(DESeq2-style) and verify it recovers the size factors you put in.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Simulate counts and write a size-factor normaliser.
# ----------------------------------------------------------------------

DISPERSION = 0.05        # alpha; biological replicate variability in real screens
MEAN_BASELINE = 500      # log-normal median sgRNA count
LIBRARY_SIZE_FACTORS_TRUE = np.array([0.85, 1.00, 1.15, 0.90, 1.10, 1.20])


def make_gene_truth(rng: np.random.Generator) -> pd.DataFrame:
    """Return a DataFrame with columns gene, category, true_lfc."""
    # TODO: build gene names like "GENE_0001"; assign first N_ESSENTIAL as
    # "essential" with true_lfc = LFC_ESSENTIAL, next N_RESISTANT as "resistant"
    # with true_lfc = LFC_RESISTANT, remaining N_NULL as "null" with true_lfc = 0.
    raise NotImplementedError


def make_sgrna_table(genes: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Return a DataFrame with columns sgrna, gene, base_mean, true_lfc."""
    # TODO: for each gene, generate N_SGRNAS_PER_GENE sgRNAs. base_mean is
    # log-normal: exp(N(log(MEAN_BASELINE), 0.6)).
    raise NotImplementedError


def simulate_counts(sgrnas: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    """Return an (N_SGRNAS, N_SAMPLES) integer count matrix.

    Control samples (cols 0..N_CONTROL-1) use base_mean; treatment samples
    (cols N_CONTROL..) use base_mean * 2**true_lfc. Per-sample multiply by
    LIBRARY_SIZE_FACTORS_TRUE. Draw NB with dispersion DISPERSION.
    """
    # TODO: vectorise with scipy.stats.nbinom or numpy.random.negative_binomial.
    raise NotImplementedError


def estimate_size_factors(counts: np.ndarray) -> np.ndarray:
    """DESeq2 median-of-ratios estimator.

    1. Compute the geometric mean of each row (sgRNA) across samples.
    2. Divide each column by that pseudo-reference; take the median (over
       sgRNAs with positive ratio) — that is the size factor for that sample.
    """
    # TODO
    raise NotImplementedError


def normalise(counts: np.ndarray, sf: np.ndarray) -> np.ndarray:
    """Return counts / size_factors (broadcast along columns)."""
    # TODO
    raise NotImplementedError


# Populate `genes`, `sgrnas`, `counts`, `sf`, `norm`.
genes = None
sgrnas = None
counts = None
sf = None
norm = None
'''


STEP_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

DISPERSION = 0.05
MEAN_BASELINE = 500
LIBRARY_SIZE_FACTORS_TRUE = np.array([0.85, 1.00, 1.15, 0.90, 1.10, 1.20])


def make_gene_truth(rng: np.random.Generator) -> pd.DataFrame:
    categories = (
        ["essential"] * N_ESSENTIAL
        + ["resistant"] * N_RESISTANT
        + ["null"] * N_NULL
    )
    true_lfc = (
        [LFC_ESSENTIAL] * N_ESSENTIAL
        + [LFC_RESISTANT] * N_RESISTANT
        + [0.0] * N_NULL
    )
    return pd.DataFrame({
        "gene":     [f"GENE_{i:04d}" for i in range(N_GENES)],
        "category": categories,
        "true_lfc": true_lfc,
    })


def make_sgrna_table(genes: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for gi, row in genes.iterrows():
        base_means = np.exp(rng.normal(np.log(MEAN_BASELINE), 0.6, size=N_SGRNAS_PER_GENE))
        for j in range(N_SGRNAS_PER_GENE):
            rows.append({
                "sgrna":     f"{row['gene']}_sg{j+1:02d}",
                "gene":      row["gene"],
                "base_mean": float(base_means[j]),
                "true_lfc":  row["true_lfc"],
            })
    return pd.DataFrame(rows)


def simulate_counts(sgrnas: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    base = sgrnas["base_mean"].to_numpy()                # (N_SGRNAS,)
    lfc  = sgrnas["true_lfc"].to_numpy()                 # (N_SGRNAS,)
    fold = np.power(2.0, lfc)                            # (N_SGRNAS,)

    # Per-sample mean matrix.
    mu = np.empty((N_SGRNAS, N_SAMPLES), dtype=float)
    for s in range(N_SAMPLES):
        is_treatment = s >= N_CONTROL
        per_sgrna = base * (fold if is_treatment else 1.0)
        mu[:, s] = per_sgrna * LIBRARY_SIZE_FACTORS_TRUE[s]

    # NB sampling: Var = mu + alpha * mu^2. scipy nbinom(n, p) with n = 1/alpha.
    alpha = DISPERSION
    n_par = 1.0 / alpha
    p_par = n_par / (n_par + mu)                         # (N_SGRNAS, N_SAMPLES)
    # Sample via numpy.
    counts = rng.negative_binomial(n_par, p_par)
    return counts.astype(np.int64)


def estimate_size_factors(counts: np.ndarray) -> np.ndarray:
    # Geometric mean across samples per sgRNA. Use log domain to be safe.
    with np.errstate(divide="ignore"):
        log_counts = np.log(counts.astype(float))
    # Only rows with no zero count contribute to the geometric mean.
    finite_rows = np.all(np.isfinite(log_counts), axis=1)
    log_geo = log_counts[finite_rows].mean(axis=1)       # (n_finite,)
    log_ratios = log_counts[finite_rows] - log_geo[:, None]
    return np.exp(np.median(log_ratios, axis=0))


def normalise(counts: np.ndarray, sf: np.ndarray) -> np.ndarray:
    return counts.astype(float) / sf[None, :]


# Build everything.
genes  = make_gene_truth(RNG)
sgrnas = make_sgrna_table(genes, RNG)
counts = simulate_counts(sgrnas, RNG)
sf     = estimate_size_factors(counts)
norm   = normalise(counts, sf)

print(f"counts matrix: shape={counts.shape}, dtype={counts.dtype}")
print(f"  median per column:    {np.median(counts, axis=0).round(1)}")
print(f"  estimated size factors: {sf.round(3)}")
print(f"  true       size factors: {LIBRARY_SIZE_FACTORS_TRUE}")
print(f"  recovery error (relative): "
      f"{np.abs(sf / sf.mean() - LIBRARY_SIZE_FACTORS_TRUE / LIBRARY_SIZE_FACTORS_TRUE.mean()).max():.3f}")

# Quick visualisation: count distributions before / after normalisation.
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
for s in range(N_SAMPLES):
    label = f"ctrl_{s+1}" if s < N_CONTROL else f"trt_{s-N_CONTROL+1}"
    axes[0].hist(np.log10(counts[:, s] + 1), bins=50, histtype="step", label=label)
    axes[1].hist(np.log10(norm[:, s] + 1), bins=50, histtype="step", label=label)
axes[0].set_title("raw counts (log10)")
axes[1].set_title("size-factor normalised (log10)")
for ax in axes:
    ax.set_xlabel("log10 (count + 1)")
    ax.set_ylabel("sgRNAs")
axes[0].legend(fontsize=8, ncol=2)
plt.tight_layout()
plt.show()
'''


STEP2_MD = """## Step 2 (12 min) — Per-sgRNA negative-binomial test

For each sgRNA we test the null that mean(control) == mean(treatment) under
a NB count model with a shared dispersion. We use the **likelihood-ratio
test** between (a) one common mean across all 6 samples and (b) two
separate means per condition. The test statistic
$$
\\Lambda = 2 (\\ell_{\\text{full}} - \\ell_{\\text{null}})
$$
is approximately χ²₁ under the null (one extra parameter). The signed
log₂-fold-change uses the size-factor-normalised counts.

Then apply **Benjamini-Hochberg FDR** to the 10 000 p-values.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Per-sgRNA NB likelihood-ratio test + BH FDR.
# ----------------------------------------------------------------------

ALPHA_NB = DISPERSION    # we know the dispersion (use a moment-of-method estimate in real life)


def nb_log_likelihood(counts_row: np.ndarray, mu_row: np.ndarray, alpha: float) -> float:
    """Sum of log NB PMFs across samples.

    Var = mu + alpha * mu**2, so n = 1/alpha, p = n / (n + mu).
    """
    # TODO: use scipy.stats.nbinom.logpmf with the (n, p) parametrisation.
    raise NotImplementedError


def sgrna_test(counts_row: np.ndarray, sf: np.ndarray, alpha: float = ALPHA_NB):
    """Return (lfc, pvalue) for one sgRNA's row of raw counts.

    Use size-factor-normalised counts to estimate the per-condition mean;
    feed (mean_per_condition * sf) into nb_log_likelihood as mu_row.
    """
    # TODO:
    # 1. Compute normalised = counts_row / sf.
    # 2. Estimate mu_ctrl, mu_trt as means of normalised over each condition.
    # 3. Build mu under full model: [mu_ctrl]*3 + [mu_trt]*3, then multiply by sf.
    # 4. Build mu under null: pooled mean across all 6, multiply by sf.
    # 5. lfc = log2((mu_trt + 0.5) / (mu_ctrl + 0.5)).
    # 6. LLR -> chi^2_1 -> p-value.
    raise NotImplementedError


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (i.e. q-values)."""
    # TODO
    raise NotImplementedError


# Loop over rows; collect lfc + pvalue.
sgrna_lfc = None
sgrna_p   = None
sgrna_q   = None
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.

ALPHA_NB = DISPERSION


def nb_log_likelihood(counts_row: np.ndarray, mu_row: np.ndarray, alpha: float) -> float:
    mu = np.clip(mu_row, 1e-6, None)
    n_par = 1.0 / alpha
    p_par = n_par / (n_par + mu)
    return float(stats.nbinom.logpmf(counts_row, n_par, p_par).sum())


def sgrna_test(counts_row: np.ndarray, sf: np.ndarray, alpha: float = ALPHA_NB):
    norm_row = counts_row / sf
    mu_ctrl = norm_row[:N_CONTROL].mean()
    mu_trt  = norm_row[N_CONTROL:].mean()
    mu_pool = norm_row.mean()

    # Pseudocount in LFC to handle zeros.
    lfc = math.log2((mu_trt + 0.5) / (mu_ctrl + 0.5))

    # Build mu_row under each model (in raw-count space, scaled by size factors).
    mu_full = np.empty(N_SAMPLES)
    mu_full[:N_CONTROL] = mu_ctrl
    mu_full[N_CONTROL:] = mu_trt
    mu_full = mu_full * sf

    mu_null = np.full(N_SAMPLES, mu_pool) * sf

    ll_full = nb_log_likelihood(counts_row, mu_full, alpha)
    ll_null = nb_log_likelihood(counts_row, mu_null, alpha)
    llr = 2.0 * (ll_full - ll_null)
    # Negative LLR can occur from numerical noise on near-null sgRNAs; clip to 0.
    llr = max(llr, 0.0)
    pval = float(stats.chi2.sf(llr, df=1))
    return lfc, pval


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1.0)
    # Enforce monotonicity from the largest down.
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


t0 = time.time()
lfc_arr = np.empty(N_SGRNAS)
p_arr   = np.empty(N_SGRNAS)
for i in range(N_SGRNAS):
    lfc_arr[i], p_arr[i] = sgrna_test(counts[i], sf)
sgrna_lfc = lfc_arr
sgrna_p   = p_arr
sgrna_q   = bh_fdr(sgrna_p)
print(f"sgRNA NB tests: {N_SGRNAS} sgRNAs in {time.time()-t0:.1f}s")

# Attach to sgrnas frame.
sgrnas = sgrnas.assign(lfc=sgrna_lfc, pval=sgrna_p, qval=sgrna_q)

# Diagnostics: median LFC per category should reflect what we planted.
for cat in ["essential", "resistant", "null"]:
    sel = sgrnas["gene"].isin(genes.loc[genes.category == cat, "gene"])
    print(f"  {cat:>9s}  n={sel.sum():5d}  median LFC = {sgrnas.loc[sel, 'lfc'].median():+.2f}  "
          f"frac q<0.05 = {(sgrnas.loc[sel, 'qval'] < 0.05).mean():.2%}")

# p-value histogram — null should be uniform, hits give an enriched left bin.
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.hist(sgrnas.loc[sgrnas["gene"].isin(genes.loc[genes.category=="null", "gene"]), "pval"],
        bins=20, alpha=0.6, label="null sgRNAs", density=True)
ax.hist(sgrnas.loc[~sgrnas["gene"].isin(genes.loc[genes.category=="null", "gene"]), "pval"],
        bins=20, alpha=0.6, label="hit sgRNAs", density=True)
ax.set_xlabel("p-value")
ax.set_ylabel("density")
ax.set_title("Per-sgRNA p-value distribution")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP3_MD = """## Step 3 (14 min) — Gene-level α-RRA aggregation and volcano plot

Each gene has 20 sgRNAs, hence 20 p-values. Robust Rank Aggregation (RRA) —
the heart of MAGeCK — works as follows:

1. **Rank** all sgRNAs globally by p-value (sorted in the direction of
   interest: smallest p first for depletion / enrichment as appropriate).
2. For each gene, convert its sgRNAs' ranks to **normalised ranks**
   `r_i = rank_i / N`, all uniform on (0, 1] under the null.
3. **α-cut:** keep only sgRNAs with `r_i ≤ α` (typically 0.25). This is the
   robustness knob — a couple of bad sgRNAs in a gene's set should not
   dominate.
4. For each kept rank, the order statistic
   `β_{i,n} = Beta(i, n-i+1).cdf(r_(i))` is the p-value for "this is the
   i-th smallest of n uniforms". The gene's RRA score is `min_i β_{i,n}`.
5. Permute (or here: use the analytic Beta-min approximation) to get a
   p-value for the gene, then BH-FDR.

You will run this twice: once with sgRNAs sorted ascending by signed LFC
(depletion / essentiality), once descending (enrichment / resistance), and
combine via min FDR with a sign tag. Then plot a volcano (gene LFC vs
−log₁₀ FDR) and colour by ground truth.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — alpha-RRA gene-level aggregation and volcano plot.
# ----------------------------------------------------------------------

ALPHA_RRA = 0.25
N_PERM = 2000   # permutation count for the null distribution


def alpha_rra_score(norm_ranks: np.ndarray, alpha: float = ALPHA_RRA) -> float:
    """Compute rho = min_i Beta(i, n-i+1).cdf(r_(i)) over kept ranks.

    norm_ranks is the array of normalised ranks (in (0, 1]) for one gene's
    sgRNAs. Sort ascending, keep those <= alpha, compute the order-statistic
    p-values, return the minimum (or 1.0 if nothing survives the alpha cut).
    """
    # TODO
    raise NotImplementedError


def rra_null(n_per_gene: int, alpha: float, n_perm: int, rng: np.random.Generator) -> np.ndarray:
    """Permutation null: sample n_perm gene-sized batches of uniform ranks,
    compute the RRA score on each, return the array of null scores."""
    # TODO
    raise NotImplementedError


def gene_rra(sgrnas: pd.DataFrame,
             genes: pd.DataFrame,
             rank_by: str,
             ascending: bool,
             alpha: float = ALPHA_RRA,
             n_perm: int = N_PERM,
             rng: np.random.Generator = RNG) -> pd.DataFrame:
    """Return a DataFrame with one row per gene: rra_score, rra_p, rra_q, mean_lfc."""
    # TODO:
    # 1. Globally rank sgRNAs by the chosen column (e.g. "lfc" ascending for depletion).
    # 2. Convert to normalised ranks.
    # 3. Per gene, compute alpha-RRA score.
    # 4. Build a permutation null at n_per_gene = sgRNAs per gene (assume uniform across genes).
    # 5. Empirical p-value for each gene = (1 + #(null <= score)) / (n_perm + 1).
    # 6. BH-FDR adjust; also return mean LFC per gene.
    raise NotImplementedError


# Compute gene-level results for depletion (ascending LFC) and enrichment (descending LFC).
gene_depletion = None
gene_enrichment = None
gene_results = None
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.

ALPHA_RRA = 0.25
N_PERM = 2000


def alpha_rra_score(norm_ranks: np.ndarray, alpha: float = ALPHA_RRA) -> float:
    r_sorted = np.sort(norm_ranks)
    kept = r_sorted[r_sorted <= alpha]
    if kept.size == 0:
        return 1.0
    n = len(r_sorted)
    # Order-statistic p-values: P[U_(i) <= r_(i)] = Beta(i, n - i + 1).cdf(r_(i)).
    i = np.arange(1, kept.size + 1)
    beta_p = stats.beta.cdf(kept, i, n - i + 1)
    return float(beta_p.min())


def rra_null(n_per_gene: int, alpha: float, n_perm: int,
             rng: np.random.Generator) -> np.ndarray:
    # Vectorised: draw (n_perm, n_per_gene) uniforms, sort each row, evaluate
    # Beta order-statistic CDFs in one pass, take the row-wise min over the
    # kept (r <= alpha) entries.
    u = rng.uniform(size=(n_perm, n_per_gene))
    u_sorted = np.sort(u, axis=1)
    i_idx = np.arange(1, n_per_gene + 1)
    # Beta(i, n - i + 1).cdf evaluated per element.
    beta_cdf = stats.beta.cdf(u_sorted, i_idx, n_per_gene - i_idx + 1)
    # Mask out r > alpha (do not count those order stats).
    beta_cdf = np.where(u_sorted <= alpha, beta_cdf, 1.0)
    return beta_cdf.min(axis=1)


def gene_rra(sgrnas: pd.DataFrame,
             genes: pd.DataFrame,
             rank_by: str,
             ascending: bool,
             alpha: float = ALPHA_RRA,
             n_perm: int = N_PERM,
             rng: np.random.Generator = RNG) -> pd.DataFrame:
    n_total = len(sgrnas)
    order = sgrnas[rank_by].rank(method="average", ascending=ascending).to_numpy()
    norm_rank_all = order / n_total
    df = sgrnas.assign(_nrank=norm_rank_all)

    scores = {}
    mean_lfc = {}
    for gene_name, sub in df.groupby("gene", sort=False):
        scores[gene_name] = alpha_rra_score(sub["_nrank"].to_numpy(), alpha)
        mean_lfc[gene_name] = sub["lfc"].mean()

    score_arr = np.array([scores[g] for g in genes["gene"]])
    null = rra_null(N_SGRNAS_PER_GENE, alpha, n_perm, rng)

    # Empirical p-value: fraction of null scores at least as small as the gene's score.
    null_sorted = np.sort(null)
    pvals = (np.searchsorted(null_sorted, score_arr, side="right") + 1) / (n_perm + 1)
    pvals = np.clip(pvals, 1.0 / (n_perm + 1), 1.0)

    qvals = bh_fdr(pvals)
    return pd.DataFrame({
        "gene":      genes["gene"].to_numpy(),
        "category":  genes["category"].to_numpy(),
        "true_lfc":  genes["true_lfc"].to_numpy(),
        "rra_score": score_arr,
        "rra_p":     pvals,
        "rra_q":     qvals,
        "mean_lfc":  np.array([mean_lfc[g] for g in genes["gene"]]),
    })


t0 = time.time()
gene_depletion  = gene_rra(sgrnas, genes, rank_by="lfc", ascending=True)
gene_enrichment = gene_rra(sgrnas, genes, rank_by="lfc", ascending=False)
print(f"alpha-RRA gene aggregation: {2*N_GENES} gene tests in {time.time()-t0:.1f}s")

# Combine: a gene is called depleted or enriched depending on which side wins.
# Take min q from the two tails, with a "direction" tag.
combined = gene_depletion.merge(
    gene_enrichment[["gene", "rra_p", "rra_q", "rra_score"]],
    on="gene", suffixes=("_dep", "_enr"),
)
combined["direction"] = np.where(
    combined["rra_q_dep"] <= combined["rra_q_enr"], "depleted", "enriched"
)
combined["best_q"] = np.minimum(combined["rra_q_dep"], combined["rra_q_enr"])
combined["best_p"] = np.where(
    combined["direction"] == "depleted", combined["rra_p_dep"], combined["rra_p_enr"]
)
gene_results = combined

# Per-category recovery.
print()
print("Per-category recovery (FDR < 0.10, correct sign):")
for cat, expected_dir in [("essential", "depleted"), ("resistant", "enriched"), ("null", None)]:
    sub = gene_results[gene_results["category"] == cat]
    if expected_dir is None:
        called = sub["best_q"] < 0.10
    else:
        called = (sub["best_q"] < 0.10) & (sub["direction"] == expected_dir)
    print(f"  {cat:>9s}  {called.sum():4d} / {len(sub):4d}  ({called.mean():.0%})")

# Volcano plot.
fig, ax = plt.subplots(figsize=(8, 5))
neg_log_q = -np.log10(np.clip(gene_results["best_q"], 1e-6, 1.0))
colours = {"essential": "#d62728", "resistant": "#2ca02c", "null": "#777777"}
for cat in ["null", "essential", "resistant"]:
    sel = gene_results["category"] == cat
    ax.scatter(gene_results.loc[sel, "mean_lfc"], neg_log_q[sel],
               s=18, alpha=0.7, color=colours[cat], label=cat,
               edgecolors="none")
ax.axhline(-np.log10(0.10), color="k", linestyle="--", lw=0.8, label="FDR = 0.10")
ax.axvline(0, color="k", lw=0.5)
ax.set_xlabel("gene mean log2 fold-change (treatment / control)")
ax.set_ylabel("-log10 (RRA FDR)")
ax.set_title("Volcano: alpha-RRA gene calls vs ground truth")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP4_MD = """## Step 4 (14 min) — sgRNA concordance and observed-vs-expected FDR

A hit gene should have **concordant** sgRNAs: most of its 20 sgRNAs point
in the same direction. A null gene should have ~10 / 10. The concordance
fraction (max of "n sgRNAs LFC < 0" / 20 and "n sgRNAs LFC > 0" / 20) is a
classic MAGeCK sanity diagnostic.

You will also compare the **claimed FDR** (BH-adjusted q-value) to the
**observed FDR** in your called set. They should match in expectation; if
they diverge it usually means the NB dispersion was misestimated (we
used the true value here, so they should agree).
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Concordance and observed-vs-expected FDR.
# ----------------------------------------------------------------------


def gene_concordance(sgrnas: pd.DataFrame) -> pd.Series:
    """Return Series indexed by gene with concordance fraction in [0.5, 1.0]."""
    # TODO: per gene compute max(frac sgRNAs LFC<0, frac sgRNAs LFC>0).
    raise NotImplementedError


def fdr_curve(gene_results: pd.DataFrame) -> pd.DataFrame:
    """Sort genes by best_q; for each threshold compute claimed vs observed FDR.

    A "true positive" = (category != 'null') AND direction matches truth.
    A "false positive" = category == 'null'.
    """
    # TODO
    raise NotImplementedError


concordance = None
fdr_df = None
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.


def gene_concordance(sgrnas: pd.DataFrame) -> pd.Series:
    grouped = sgrnas.groupby("gene")["lfc"]
    frac_neg = grouped.apply(lambda x: (x < 0).mean())
    frac_pos = grouped.apply(lambda x: (x > 0).mean())
    return pd.concat([frac_neg, frac_pos], axis=1).max(axis=1)


concordance = gene_concordance(sgrnas)
gene_results = gene_results.merge(
    concordance.rename("concordance").reset_index(), on="gene"
)

# Diagnostic: median concordance by category.
print("Median sgRNA concordance per category:")
for cat in ["essential", "resistant", "null"]:
    sel = gene_results["category"] == cat
    print(f"  {cat:>9s}  median = {gene_results.loc[sel, 'concordance'].median():.2f}  "
          f"(min {gene_results.loc[sel, 'concordance'].min():.2f}, "
          f"max {gene_results.loc[sel, 'concordance'].max():.2f})")

# Observed vs claimed FDR.
def fdr_curve(gene_results: pd.DataFrame) -> pd.DataFrame:
    g = gene_results.sort_values("best_q").reset_index(drop=True)
    # "True hit" definition: non-null and direction matches.
    correct_direction = (
        ((g["category"] == "essential") & (g["direction"] == "depleted"))
        | ((g["category"] == "resistant") & (g["direction"] == "enriched"))
    )
    is_null_or_wrongdir = ~correct_direction
    n_called = np.arange(1, len(g) + 1)
    observed_fp = is_null_or_wrongdir.cumsum().to_numpy()
    return pd.DataFrame({
        "rank":          n_called,
        "claimed_q":     g["best_q"].to_numpy(),
        "observed_fdr":  observed_fp / n_called,
        "gene":          g["gene"].to_numpy(),
        "category":      g["category"].to_numpy(),
    })


fdr_df = fdr_curve(gene_results)
print()
print("First 5 / last few rows of FDR curve:")
print(fdr_df.head().to_string(index=False))

# Plot.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))

# Concordance distributions.
for cat, colour in [("null", "#777777"), ("essential", "#d62728"), ("resistant", "#2ca02c")]:
    sel = gene_results["category"] == cat
    axes[0].hist(gene_results.loc[sel, "concordance"], bins=20,
                 alpha=0.6, color=colour, label=cat, density=True)
axes[0].axvline(0.5, color="k", lw=0.5)
axes[0].set_xlabel("max(frac neg, frac pos) sgRNAs per gene")
axes[0].set_ylabel("density")
axes[0].set_title("sgRNA concordance")
axes[0].legend()

# Observed vs claimed FDR.
axes[1].plot(fdr_df["claimed_q"], fdr_df["observed_fdr"], lw=1.5, label="observed")
axes[1].plot([0, 1], [0, 1], "k--", lw=0.8, label="y = x (ideal)")
axes[1].set_xlim(0, 0.5)
axes[1].set_ylim(0, 0.5)
axes[1].set_xlabel("claimed FDR (BH q-value)")
axes[1].set_ylabel("observed FDR (in called set)")
axes[1].set_title("Calibration of BH FDR")
axes[1].legend()

plt.tight_layout()
plt.show()

# Summary at FDR = 0.10.
mask = fdr_df["claimed_q"] < 0.10
if mask.any():
    last = fdr_df[mask].iloc[-1]
    print(f"\\nAt FDR < 0.10: called {int(last['rank'])} genes, "
          f"observed FDR = {last['observed_fdr']:.1%}")
'''


STEP5_MD = """## Step 5 (12 min) — Screen as compressed sensing

We end with the EE punchline. The whole screen can be written as
$$
\\mathbf{y} = \\Phi \\mathbf{x} + \\boldsymbol\\varepsilon
$$
where

- `x ∈ ℝ^G` is the **gene fitness vector** — sparse (only ~70 of 500 genes
  are non-zero in this screen).
- `Φ ∈ {0,1}^{S×G}` is the **sgRNA × gene perturbation matrix** — exactly
  20 ones per row, all in the same column (each sgRNA targets one gene).
- `y ∈ ℝ^S` is the **measured per-sgRNA log-fold-change** — noisy.

That makes recovery a sparse-inverse-sensing problem. Multiple sgRNAs per
gene give the **redundancy** (repetition coding) that lets us beat the
per-measurement noise.

You will:

1. Build Φ, set `x` from the ground truth, simulate `y = Φx + ε`, and
   solve via least squares (closed-form because Φ has orthogonal columns
   here).
2. Compare to the **`mean LFC per gene`** estimator we used in Step 3 —
   they should be identical.
3. Sweep the **number of sgRNAs per gene** (k ∈ {2, 4, 10, 20}) by
   subsampling and report the recovery RMSE — it falls as 1/√k, the
   matched-filter / repetition-code rate.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Screen as a sparse-recovery problem; sgRNAs as a repetition code.
# ----------------------------------------------------------------------


def build_perturbation_matrix(sgrnas: pd.DataFrame, genes: pd.DataFrame) -> np.ndarray:
    """Return Phi of shape (N_SGRNAS, N_GENES)."""
    # TODO: one-hot encode the sgRNA -> gene mapping.
    raise NotImplementedError


def least_squares_solve(Phi: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Closed-form OLS: x_hat = (Phi^T Phi)^-1 Phi^T y."""
    # TODO: use np.linalg.lstsq.
    raise NotImplementedError


def sweep_redundancy(sgrnas: pd.DataFrame, genes: pd.DataFrame,
                     ks: list[int], rng: np.random.Generator) -> pd.DataFrame:
    """For each k in ks, sub-sample to k sgRNAs per gene, refit, return RMSE."""
    # TODO
    raise NotImplementedError


Phi = None
x_hat = None
rmse_curve = None
'''


STEP5_SOLUTION = '''# Reference solution — Step 5.


def build_perturbation_matrix(sgrnas: pd.DataFrame, genes: pd.DataFrame) -> np.ndarray:
    gene_index = {g: i for i, g in enumerate(genes["gene"])}
    Phi = np.zeros((len(sgrnas), len(genes)), dtype=np.float32)
    for s_idx, row in enumerate(sgrnas.itertuples(index=False)):
        Phi[s_idx, gene_index[row.gene]] = 1.0
    return Phi


def least_squares_solve(Phi: np.ndarray, y: np.ndarray) -> np.ndarray:
    # With orthogonal columns this reduces to per-gene mean.
    x_hat, *_ = np.linalg.lstsq(Phi, y, rcond=None)
    return x_hat


Phi = build_perturbation_matrix(sgrnas, genes)
y = sgrnas["lfc"].to_numpy()
x_hat = least_squares_solve(Phi, y)

# x_hat should match the gene mean LFC we computed earlier.
mean_lfc_by_gene = gene_results.set_index("gene")["mean_lfc"].reindex(genes["gene"]).to_numpy()
assert np.allclose(x_hat, mean_lfc_by_gene, atol=1e-5), "LS disagrees with per-gene mean"
print(f"OLS recovery == mean-per-gene (max abs diff = "
      f"{np.abs(x_hat - mean_lfc_by_gene).max():.2e})")
print(f"  true sparsity: {(genes['true_lfc'] != 0).sum()} / {N_GENES} non-zero genes")


def sweep_redundancy(sgrnas: pd.DataFrame, genes: pd.DataFrame,
                     ks: list[int], rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    truth = genes["true_lfc"].to_numpy()
    for k in ks:
        # Subsample k sgRNAs per gene.
        keep = sgrnas.groupby("gene", group_keys=False).apply(
            lambda g: g.sample(n=k, random_state=int(rng.integers(0, 1_000_000)))
        )
        Phi_k = build_perturbation_matrix(keep, genes)
        y_k = keep["lfc"].to_numpy()
        x_hat_k = least_squares_solve(Phi_k, y_k)
        rmse = float(np.sqrt(np.mean((x_hat_k - truth) ** 2)))
        rows.append({"k_sgRNAs_per_gene": k, "rmse": rmse})
    return pd.DataFrame(rows)


rmse_curve = sweep_redundancy(sgrnas, genes, ks=[2, 4, 10, 20], rng=RNG)
print()
print(rmse_curve.to_string(index=False))

# Theoretical: at k sgRNAs per gene, RMSE per gene-mean drops like sigma / sqrt(k).
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(rmse_curve["k_sgRNAs_per_gene"], rmse_curve["rmse"], "o-",
        label="empirical RMSE", lw=2)
sigma = rmse_curve["rmse"].iloc[0] * np.sqrt(rmse_curve["k_sgRNAs_per_gene"].iloc[0])
ax.plot(rmse_curve["k_sgRNAs_per_gene"],
        sigma / np.sqrt(rmse_curve["k_sgRNAs_per_gene"]),
        "k--", lw=1, label=r"$\\sigma / \\sqrt{k}$")
ax.set_xlabel("sgRNAs per gene (k)")
ax.set_ylabel("RMSE of recovered gene LFC")
ax.set_title("Repetition coding: noise floor falls as 1 / sqrt(k)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()
plt.tight_layout()
plt.show()

print()
print("Reads: more sgRNAs per gene = lower noise floor, exactly as a repetition")
print("code stacks N independent observations to gain sqrt(N) SNR.")
'''


SELFCHECK_MD = """## Self-check

If the reference solutions ran, these asserts validate the load-bearing
numerical pieces of the analysis. If you wrote your own code and one
fails, revisit the relevant step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Count matrix has the right shape and is non-negative integer.
assert counts.shape == (N_SGRNAS, N_SAMPLES), f"bad counts shape: {counts.shape}"
assert counts.dtype.kind in "iu", "counts should be integer"
assert (counts >= 0).all(), "negative counts?"

# 2. Size factors recover (up to common rescale) the true values to within 10%.
sf_norm = sf / sf.mean()
true_norm = LIBRARY_SIZE_FACTORS_TRUE / LIBRARY_SIZE_FACTORS_TRUE.mean()
rel_err = np.abs(sf_norm - true_norm).max()
assert rel_err < 0.10, f"size factors off by {rel_err:.2%}"

# 3. The median LFC for essential / resistant / null sgRNAs is roughly right.
for cat, expected_sign, target in [
    ("essential", "neg", LFC_ESSENTIAL),
    ("resistant", "pos", LFC_RESISTANT),
    ("null",      "zero", 0.0),
]:
    sel = sgrnas["gene"].isin(genes.loc[genes.category == cat, "gene"])
    med = float(sgrnas.loc[sel, "lfc"].median())
    if expected_sign == "neg":
        assert med < -0.5, f"essential median LFC = {med:.2f}, expected << 0"
        assert abs(med - target) < 1.0, f"essential median LFC {med:.2f} far from {target}"
    elif expected_sign == "pos":
        assert med > 0.3, f"resistant median LFC = {med:.2f}, expected > 0"
        assert abs(med - target) < 0.7, f"resistant median LFC {med:.2f} far from {target}"
    else:
        assert abs(med) < 0.15, f"null median LFC = {med:.2f}, expected ~ 0"

# 4. The RRA pipeline recovers >= 80% of essentials and >= 70% of resistants at FDR < 0.10
#    with correct sign, and labels <= 5% of nulls.
def hit_rate(cat: str, direction: str) -> float:
    sub = gene_results[gene_results["category"] == cat]
    called = (sub["best_q"] < 0.10) & (sub["direction"] == direction)
    return called.mean()

ess_rate = hit_rate("essential", "depleted")
res_rate = hit_rate("resistant", "enriched")
null_called = ((gene_results["category"] == "null") & (gene_results["best_q"] < 0.10)).mean()
print(f"Essential recovered: {ess_rate:.0%}")
print(f"Resistant recovered: {res_rate:.0%}")
print(f"Null mis-called:     {null_called:.1%}")
assert ess_rate >= 0.80, f"essential recovery only {ess_rate:.0%}"
assert res_rate >= 0.70, f"resistant recovery only {res_rate:.0%}"
assert null_called <= 0.05, f"null mis-call rate {null_called:.1%} too high"

# 5. Observed FDR in called set should be close to claimed at q=0.10 (within 5pp).
mask = fdr_df["claimed_q"] < 0.10
if mask.any():
    obs = fdr_df.loc[mask, "observed_fdr"].iloc[-1]
    assert obs <= 0.15, f"observed FDR {obs:.1%} much higher than claimed 10%"

# 6. Repetition-code scaling: RMSE at k=20 << RMSE at k=2.
r2 = rmse_curve.loc[rmse_curve["k_sgRNAs_per_gene"] == 2, "rmse"].iloc[0]
r20 = rmse_curve.loc[rmse_curve["k_sgRNAs_per_gene"] == 20, "rmse"].iloc[0]
ratio = r2 / r20
# Theory says sqrt(20/2) ~ 3.16; allow a wide band for the finite-sample empirical.
assert 2.0 <= ratio <= 5.0, f"RMSE ratio k=2 / k=20 = {ratio:.2f}, expected ~3"

print("\\n✅ Self-check passed.")
'''


EE_MD = """## EE framing — repetition coding, matched filter, sparse recovery

You have just built three EE ideas back-to-back:

1. **Matched filter in the count domain.** The NB likelihood-ratio test
   integrates the per-sample counts against the alternative-hypothesis
   mean (treatment-vs-control), weighted by the variance model. Replace
   "samples" with "time bins" and "NB" with "Gaussian" and you have the
   discrete matched filter from your DSP class. The size-factor
   normalisation is the gain-control step that makes the filter
   shift-invariant across sample depths.

2. **Repetition code → matched-filter SNR gain.** A gene's 20 sgRNAs are
   20 independent noisy estimates of the same scalar (the gene's
   fitness contribution). Averaging them drops noise by √20. The
   final RMSE-vs-k plot in Step 5 reproduces the classic
   `1 / √N` repetition-code curve.

3. **Compressed sensing on a sparse perturbation matrix.** The
   sgRNA-to-gene matrix Φ is binary, super-sparse, and has orthogonal
   columns by construction (each sgRNA targets one gene). The fitness
   vector x is itself sparse (~14% non-zero here). Recovery is OLS in
   the easy case; in real screens, where sgRNAs have off-target effects,
   Φ has small but non-zero off-diagonal entries and Lasso / RRA-style
   robust aggregation handles the dependent measurements.

That last bit is why MAGeCK, BAGEL2, DrugZ all exist as separate tools:
they make different assumptions about Φ's off-target structure and the
noise distribution, but the **inverse-problem framing is shared**.
DepMap is then just this analysis run on ~1 000 cell lines, with the
**recovered x-vectors stacked** to give a (gene × cell line) dependency
matrix — the cancer-vulnerability atlas.
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
        md(STEP_SOLUTION_HEADER),
        hidden(STEP1_SOLUTION),

        md(STEP2_MD),
        code(STEP2_TODO),
        md(STEP_SOLUTION_HEADER),
        hidden(STEP2_SOLUTION),

        md(STEP3_MD),
        code(STEP3_TODO),
        md(STEP_SOLUTION_HEADER),
        hidden(STEP3_SOLUTION),

        md(STEP4_MD),
        code(STEP4_TODO),
        md(STEP_SOLUTION_HEADER),
        hidden(STEP4_SOLUTION),

        md(STEP5_MD),
        code(STEP5_TODO),
        md(STEP_SOLUTION_HEADER),
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
