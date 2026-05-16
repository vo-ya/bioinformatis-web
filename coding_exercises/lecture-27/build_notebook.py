"""Build exercise.ipynb for L27 — Mass-Spectrometry Proteomics + Metabolomics.

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


SOL_HEADER = "*Click ▶ to expand the reference solution.*"


# ---------------------------------------------------------------------------
# Cell sources
# ---------------------------------------------------------------------------


TITLE_MD = """# L27 — Mass-Spectrometry Proteomics + Metabolomics

In this exercise you implement the **target-decoy FDR** machinery that sits at
the heart of every shotgun-proteomics pipeline, then stack two more layers on
top: **parsimonious protein inference** from shared peptides, and **label-free
quantification** with a volcano plot. The whole pipeline runs on 20 000
synthetic PSMs and 100 proteins × 6 samples — all generated in-notebook from a
deterministic seed.
"""


AHA_MD = """> **Aha.** A reverse / shuffled **decoy** database gives you an
> *empirical null* for your PSM scoring function. The fraction of decoys above
> any threshold is a direct, unbiased estimate of the false-discovery rate —
> no parametric distribution required. Everything downstream (protein
> inference, LFQ, volcano) is multiple-testing built on top of that one
> empirical null.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 matplotlib==3.8.4 -q
"""


IMPORTS = """import math
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Deterministic for the whole notebook.
SEED = 42
rng = np.random.default_rng(SEED)
np.random.seed(SEED)

print("Imports OK.  numpy", np.__version__, " pandas", pd.__version__)
"""


# ---------------------------------------------------------------------------
# Step 1 — simulate PSMs
# ---------------------------------------------------------------------------


STEP1_MD = """## Step 1 (8 min) — Simulate the PSM list

In a real LC-MS/MS pipeline, each MS2 spectrum is searched against (a) a
**target** protein database and (b) a **decoy** database built by reversing
or shuffling the targets. Each search yields a **PSM** (peptide-spectrum
match) with a search-engine score — here we use a Comet-style **XCorr**.
Targets are a mix of true positives (drawn from a high-scoring distribution)
and false positives (which behave statistically like decoys). Decoys are
all-false by construction.

The model used in this notebook:

- **20 000 PSMs total**: 10 000 target + 10 000 decoy.
- **Targets**: 60% true matches drawn from `N(35, 5)`, 40% false matches drawn from `N(15, 5)`.
- **Decoys**: all from `N(15, 5)` — same false-match distribution as the bad targets.

Each PSM also carries a **peptide sequence** and the **protein** it maps to;
we generate 100 target proteins and randomly assign peptides, sometimes
sharing a peptide between proteins so that protein inference has something
to think about.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Simulate 20 000 PSMs with XCorr scores + peptide/protein metadata.
# ----------------------------------------------------------------------

N_TARGETS  = 10_000
N_DECOYS   = 10_000
N_PROTEINS = 100
TRUE_FRAC  = 0.60   # fraction of target PSMs that are real matches


def simulate_psms(rng) -> pd.DataFrame:
    """Return a DataFrame with columns:
       psm_id, label ('target'/'decoy'), xcorr, peptide, protein
    """
    # TODO:
    #   1. Generate target XCorr scores: a TRUE_FRAC mixture from N(35,5) + N(15,5).
    #   2. Generate decoy XCorr scores: all from N(15, 5).
    #   3. Make peptide identifiers (e.g. "PEP00001") and assign each peptide
    #      to one of N_PROTEINS proteins. Inject ~5% shared peptides that map
    #      to two proteins (record the *first* of the two as primary; we'll
    #      handle the second via a side-table in Step 3).
    #   4. Concatenate target + decoy rows; shuffle; reset index.
    raise NotImplementedError


# psms = simulate_psms(rng)
# print(psms.head())
# print(psms["label"].value_counts())
'''


STEP1_SOLUTION = '''# Reference solution — Step 1.

N_TARGETS  = 10_000
N_DECOYS   = 10_000
N_PROTEINS = 100
TRUE_FRAC  = 0.60


def simulate_psms(rng) -> pd.DataFrame:
    # --- Target scores: TRUE_FRAC from the "good" distribution, rest from null.
    n_true = int(round(N_TARGETS * TRUE_FRAC))
    n_bad  = N_TARGETS - n_true
    target_scores = np.concatenate([
        rng.normal(loc=35.0, scale=5.0, size=n_true),
        rng.normal(loc=15.0, scale=5.0, size=n_bad),
    ])
    rng.shuffle(target_scores)

    # --- Decoy scores: all from the empirical null.
    decoy_scores = rng.normal(loc=15.0, scale=5.0, size=N_DECOYS)

    # --- Peptide / protein metadata. ---------------------------------------
    # 8000 distinct target peptides assigned uniformly across 100 proteins.
    # Some peptides will be hit by multiple PSMs (the same peptide can match
    # several MS2 spectra) — that is how protein-level evidence accrues.
    n_unique_peps = 8000
    pep_ids       = np.array([f"PEP{i:05d}" for i in range(n_unique_peps)])
    pep_protein   = rng.integers(0, N_PROTEINS, size=n_unique_peps)

    # ~5% of peptides are "shared" — they map to a second protein as well.
    n_shared = int(0.05 * n_unique_peps)
    shared_mask = np.zeros(n_unique_peps, dtype=bool)
    shared_mask[:n_shared] = True
    rng.shuffle(shared_mask)
    pep_protein_b = np.where(
        shared_mask,
        rng.integers(0, N_PROTEINS, size=n_unique_peps),
        -1,
    )

    target_pep_idx = rng.integers(0, n_unique_peps, size=N_TARGETS)
    target_peps    = pep_ids[target_pep_idx]
    target_prots   = np.array([f"PROT{p:03d}" for p in pep_protein[target_pep_idx]])
    target_prots_b = np.array([
        f"PROT{p:03d}" if p >= 0 else ""
        for p in pep_protein_b[target_pep_idx]
    ])

    # Decoys map to "decoy-only" proteins so they cannot be confused with targets.
    decoy_pep_ids  = np.array([f"DECPEP{i:05d}" for i in range(N_DECOYS)])
    decoy_prots    = np.array([f"DECPROT{rng.integers(0, N_PROTEINS):03d}" for _ in range(N_DECOYS)])

    target_df = pd.DataFrame({
        "psm_id":    [f"PSM{i:06d}" for i in range(N_TARGETS)],
        "label":     "target",
        "xcorr":     target_scores,
        "peptide":   target_peps,
        "protein":   target_prots,
        "protein_b": target_prots_b,
    })
    decoy_df = pd.DataFrame({
        "psm_id":    [f"PSM{i:06d}" for i in range(N_TARGETS, N_TARGETS + N_DECOYS)],
        "label":     "decoy",
        "xcorr":     decoy_scores,
        "peptide":   decoy_pep_ids,
        "protein":   decoy_prots,
        "protein_b": "",
    })

    psms = pd.concat([target_df, decoy_df], ignore_index=True)
    psms = psms.sample(frac=1.0, random_state=SEED).reset_index(drop=True)
    return psms


psms = simulate_psms(rng)
print(f"Simulated {len(psms):,} PSMs  "
      f"({(psms['label'] == 'target').sum():,} target, "
      f"{(psms['label'] == 'decoy').sum():,} decoy)")
print()
print(psms.head())
print()
print("XCorr summary by label:")
print(psms.groupby("label")["xcorr"].agg(["mean", "std", "min", "max"]).round(2))

# Quick visual: the two distributions.
fig, ax = plt.subplots(figsize=(7, 3.5))
bins = np.linspace(0, 60, 60)
ax.hist(psms.loc[psms["label"] == "target", "xcorr"], bins=bins,
        alpha=0.6, label="target", color="C0")
ax.hist(psms.loc[psms["label"] == "decoy",  "xcorr"], bins=bins,
        alpha=0.6, label="decoy",  color="C3")
ax.set_xlabel("XCorr score")
ax.set_ylabel("PSMs")
ax.set_title("Target vs decoy score distributions")
ax.legend()
plt.tight_layout()
plt.show()
'''


# ---------------------------------------------------------------------------
# Step 2 — FDR sweep and ROC
# ---------------------------------------------------------------------------


STEP2_MD = """## Step 2 (10 min) — Sweep score thresholds → FDR curve

For any threshold $t$, define

$$\\widehat{\\text{FDR}}(t) = \\frac{\\#\\{\\text{decoys with score} \\ge t\\}}{\\#\\{\\text{targets with score} \\ge t\\}}.$$

Sweeping $t$ from high to low and plotting FDR vs target count traces out the
empirical operating curve. The classical cutoff is **1% FDR**. We will also
**monotonise** the raw curve so FDR never decreases as we add hits (a
standard "step-down" / q-value trick), and compute an **ROC** by treating
decoys as negatives.
"""


STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Threshold sweep, FDR curve, q-values, ROC.
# ----------------------------------------------------------------------


def fdr_sweep(psms: pd.DataFrame) -> pd.DataFrame:
    """For every distinct threshold, return DataFrame with:
       threshold, n_target, n_decoy, fdr, q_value
    Sort descending in `threshold` and accumulate counts.
    """
    # TODO:
    #   1. Sort PSMs by xcorr descending.
    #   2. cumulative count of targets and decoys.
    #   3. fdr = n_decoy / n_target  (clip n_target to >= 1).
    #   4. q_value = monotone-from-bottom minimum of fdr (i.e. q[i] = min(fdr[i:])).
    raise NotImplementedError


def threshold_for_fdr(sweep: pd.DataFrame, target_fdr: float = 0.01) -> float:
    """Return the *lowest* threshold whose q-value is still <= target_fdr."""
    # TODO
    raise NotImplementedError


# sweep = fdr_sweep(psms)
# t1pct = threshold_for_fdr(sweep, 0.01)
# Then plot fdr vs n_target and a target-vs-decoy ROC.
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.


def fdr_sweep(psms: pd.DataFrame) -> pd.DataFrame:
    df = psms.sort_values("xcorr", ascending=False).reset_index(drop=True)
    is_target = (df["label"] == "target").to_numpy()
    is_decoy  = (df["label"] == "decoy").to_numpy()

    n_target = np.cumsum(is_target)
    n_decoy  = np.cumsum(is_decoy)

    fdr = n_decoy / np.maximum(n_target, 1)
    # Monotone q-values: q[i] = min(fdr[i:])  (running minimum from the right).
    q_value = np.minimum.accumulate(fdr[::-1])[::-1]

    return pd.DataFrame({
        "threshold": df["xcorr"].to_numpy(),
        "n_target":  n_target,
        "n_decoy":   n_decoy,
        "fdr":       fdr,
        "q_value":   q_value,
    })


def threshold_for_fdr(sweep: pd.DataFrame, target_fdr: float = 0.01) -> float:
    ok = sweep[sweep["q_value"] <= target_fdr]
    if len(ok) == 0:
        # Nothing passes the cutoff — return +inf (no PSMs survive).
        return float("inf")
    # Largest n_target row with q <= target_fdr — that's the deepest cutoff we
    # can pick while still satisfying FDR.
    return float(ok["threshold"].iloc[-1])


sweep = fdr_sweep(psms)

t1pct = threshold_for_fdr(sweep, 0.01)
t5pct = threshold_for_fdr(sweep, 0.05)

n_at_1pct = int(sweep[sweep["q_value"] <= 0.01]["n_target"].max())
n_at_5pct = int(sweep[sweep["q_value"] <= 0.05]["n_target"].max())

print(f"XCorr threshold at 1% FDR: {t1pct:.2f}  ->  {n_at_1pct:,} target PSMs survive")
print(f"XCorr threshold at 5% FDR: {t5pct:.2f}  ->  {n_at_5pct:,} target PSMs survive")

# --- Plot 1: FDR curve (n_target on x, FDR on y) ----------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4))

ax = axes[0]
ax.plot(sweep["n_target"], sweep["fdr"],    color="C7", alpha=0.5, label="raw FDR")
ax.plot(sweep["n_target"], sweep["q_value"], color="C0",            label="q-value (monotone)")
ax.axhline(0.01, color="red", lw=1, ls="--", label="1% cutoff")
ax.set_xlabel("# target PSMs accepted (rank)")
ax.set_ylabel("FDR  =  #decoy / #target")
ax.set_title("Target-decoy FDR curve")
ax.set_ylim(0, 0.5)
ax.legend(loc="upper left")

# --- Plot 2: ROC (decoys as negatives) -------------------------------------
ax = axes[1]
ax.plot(sweep["n_decoy"], sweep["n_target"], color="C0")
ax.plot([0, sweep["n_decoy"].max()],
        [0, sweep["n_decoy"].max()], "k--", alpha=0.5, label="random")
ax.set_xlabel("# decoy PSMs accepted")
ax.set_ylabel("# target PSMs accepted")
ax.set_title("Empirical ROC")
ax.legend(loc="lower right")

plt.tight_layout()
plt.show()
'''


# ---------------------------------------------------------------------------
# Step 3 — Protein inference (parsimony + >=2 unique peptides)
# ---------------------------------------------------------------------------


STEP3_MD = """## Step 3 (12 min) — Protein inference: parsimony + ≥ 2 unique peptides

We have a list of high-confidence target PSMs. To answer the biologists'
actual question — "which **proteins** are in the sample?" — we collapse
peptides up to proteins with two rules:

1. **Parsimony.** When a peptide is shared between several proteins, assign
   it to the protein with the most *unique* peptide evidence. This is
   Occam's razor: don't postulate two proteins when one already covers the
   peptide.
2. **≥ 2 unique peptides per protein.** Any inferred protein must have at
   least two peptides that map *only* to it; one-hit-wonders are dropped.
   This is the standard MaxQuant / FragPipe filter and roughly enforces
   protein-level FDR control on top of peptide-level FDR.

Operate only on PSMs that survive the 1% FDR cutoff from Step 2.
"""


STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Parsimonious protein inference.
# ----------------------------------------------------------------------


def infer_proteins(passing: pd.DataFrame, min_unique: int = 2) -> pd.DataFrame:
    """Given the surviving target PSMs (with `peptide`, `protein`, `protein_b`),
    return a DataFrame of inferred proteins with columns:
       protein, n_unique_peptides, n_total_psms
    Apply parsimony for shared peptides, then filter by `min_unique`.
    """
    # TODO:
    #   1. Build {peptide: set(proteins_it_maps_to)} from the passing PSMs.
    #      (Use both protein and protein_b — protein_b is empty for non-shared.)
    #   2. For each protein, count its *unique* peptides (those that map only to it).
    #   3. For each shared peptide, assign it to the protein with the most unique
    #      peptides (parsimony). Break ties by lexicographic protein id.
    #   4. Aggregate to per-protein counts of unique peptides + total PSMs.
    #   5. Drop any protein with < min_unique unique peptides.
    raise NotImplementedError


# passing = psms[(psms['label'] == 'target') & (psms['xcorr'] >= t1pct)].copy()
# proteins = infer_proteins(passing, min_unique=2)
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.


def infer_proteins(passing: pd.DataFrame, min_unique: int = 2) -> pd.DataFrame:
    # ---- 1. Build peptide -> set(proteins) ---------------------------------
    pep_to_prots: dict[str, set[str]] = defaultdict(set)
    for pep, a, b in zip(passing["peptide"], passing["protein"], passing["protein_b"]):
        pep_to_prots[pep].add(a)
        if isinstance(b, str) and b:
            pep_to_prots[pep].add(b)

    # ---- 2. Count uniques (peptides mapping to exactly one protein) -------
    unique_count = defaultdict(int)
    for pep, prots in pep_to_prots.items():
        if len(prots) == 1:
            unique_count[next(iter(prots))] += 1

    # ---- 3. Parsimony for shared peptides ---------------------------------
    #     Assign each shared peptide to its "best" protein.
    pep_assigned: dict[str, str] = {}
    for pep, prots in pep_to_prots.items():
        if len(prots) == 1:
            pep_assigned[pep] = next(iter(prots))
        else:
            # Pick protein with most unique support; tiebreak alphabetical.
            pep_assigned[pep] = max(prots, key=lambda p: (unique_count[p], -ord(p[-1])))

    # ---- 4. Aggregate to per-protein counts -------------------------------
    psms_per_protein = defaultdict(int)
    peps_per_protein = defaultdict(set)
    for pep, a in zip(passing["peptide"], passing["protein"]):
        # PSMs are counted once per row.
        prot = pep_assigned.get(pep, a)
        psms_per_protein[prot] += 1
        peps_per_protein[prot].add(pep)

    rows = []
    for prot, peps in peps_per_protein.items():
        # Re-count uniques *under the assignment* (a shared peptide assigned to
        # this protein still counts as a peptide of evidence — but for the
        # "unique" filter we keep the original strict-uniqueness rule).
        strict_uniques = sum(1 for p in peps if len(pep_to_prots[p]) == 1)
        rows.append({
            "protein":            prot,
            "n_unique_peptides":  strict_uniques,
            "n_peptides":         len(peps),
            "n_psms":             psms_per_protein[prot],
        })

    out = pd.DataFrame(rows).sort_values("n_psms", ascending=False).reset_index(drop=True)
    # ---- 5. >= min_unique filter ------------------------------------------
    return out[out["n_unique_peptides"] >= min_unique].reset_index(drop=True)


passing = psms[(psms["label"] == "target") & (psms["xcorr"] >= t1pct)].copy()
print(f"PSMs surviving 1% FDR: {len(passing):,}")

proteins = infer_proteins(passing, min_unique=2)
print(f"Inferred proteins (>= 2 unique peptides): {len(proteins)}")
print()
print(proteins.head(10))
'''


# ---------------------------------------------------------------------------
# Step 4 — LFQ + volcano
# ---------------------------------------------------------------------------


STEP4_MD = """## Step 4 (15 min) — Label-free quantification + volcano plot

For each inferred protein, simulate a **6-sample** experiment (3 control,
3 treated). Quantification follows the standard LFQ recipe:

1. **TopN MS1 intensity per protein** — sum the three most intense peptide
   features. Real tools (MaxLFQ, Top3, iBAQ) all share this top-of-stack
   flavour.
2. **Median-of-ratios normalisation** — for each sample, divide every
   protein's intensity by the geometric mean across samples, take the
   median of those ratios, and rescale the column. This is DESeq2's
   normalisation on a different substrate.
3. **Log₂ fold change** between treated and control means.
4. **Welch's t-test** per protein; cap p-values at the resolution limit.
5. **Volcano plot**: log₂ FC on the x-axis, $-\\log_{10}p$ on the y-axis.
   Twenty pre-spiked-in upregulated proteins (4× effect) should pop out.

Note: we are *simulating MS1 intensity directly* rather than going through
peptide-level integration, because (a) Step 2's PSMs are scores not
intensities, and (b) the statistical structure of the volcano is what the
exercise is teaching.
"""


STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Simulate intensity, normalise, t-test, volcano.
# ----------------------------------------------------------------------

N_UPREGULATED = 20


def simulate_intensity_table(proteins: pd.DataFrame, rng) -> pd.DataFrame:
    """Return wide DataFrame indexed by protein, with 6 columns:
       C1, C2, C3, T1, T2, T3
    Each entry is a positive intensity.  Pick 20 proteins to be upregulated
    in the treated samples by ~2-3x.
    """
    # TODO:
    #   1. Base intensity per protein: lognormal(mu=12, sigma=1.0).
    #   2. Sample-to-sample noise per (protein, sample): lognormal(0, 0.15).
    #   3. Pick `N_UPREGULATED` proteins at random; multiply their T1..T3 columns by 4.0.
    #   4. Stash the spiked-in protein set on `df.attrs['upregulated']` for the
    #      self-check at the end. Return the wide DataFrame.
    raise NotImplementedError


def median_of_ratios(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise each column so the median (column / geometric-mean reference) = 1."""
    # TODO:
    #   1. Compute the geometric mean of each protein row (over all samples).
    #   2. For each column, divide by that geometric mean; take the median ratio.
    #   3. Divide each column by its median ratio.
    raise NotImplementedError


def per_protein_ttest(df: pd.DataFrame, controls: list, treated: list) -> pd.DataFrame:
    """For each row, Welch t-test on log2(intensity); return DataFrame with:
       protein, log2fc, p_value, neglog10_p
    """
    # TODO: use scipy.stats.ttest_ind(equal_var=False).
    raise NotImplementedError


# intensity = simulate_intensity_table(proteins, rng)
# intensity_norm = median_of_ratios(intensity)
# stats_df = per_protein_ttest(intensity_norm,
#                              controls=['C1','C2','C3'], treated=['T1','T2','T3'])
# Then make a volcano plot.
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.

N_UPREGULATED = 20
SAMPLES_CTRL  = ["C1", "C2", "C3"]
SAMPLES_TRT   = ["T1", "T2", "T3"]
ALL_SAMPLES   = SAMPLES_CTRL + SAMPLES_TRT


def simulate_intensity_table(proteins: pd.DataFrame, rng) -> pd.DataFrame:
    p = proteins["protein"].to_numpy()
    n = len(p)

    base = rng.lognormal(mean=12.0, sigma=1.0, size=n)  # per-protein abundance
    noise = rng.lognormal(mean=0.0, sigma=0.15, size=(n, 6))

    intensity = (base[:, None] * noise)

    # Spike in upregulation in treated samples.
    if n > 0:
        up_idx = rng.choice(n, size=min(N_UPREGULATED, n), replace=False)
        intensity[up_idx, 3:] *= 4.0
    else:
        up_idx = np.array([], dtype=int)

    df = pd.DataFrame(intensity, columns=ALL_SAMPLES)
    df.insert(0, "protein", p)
    df.attrs["upregulated"] = set(p[up_idx])
    return df


def median_of_ratios(df: pd.DataFrame) -> pd.DataFrame:
    samples = [c for c in df.columns if c != "protein"]
    X = df[samples].to_numpy().astype(float)

    # Geometric mean per protein row (over samples). Use log-mean for stability.
    log_X = np.log(np.maximum(X, 1e-12))
    geo_mean = np.exp(log_X.mean(axis=1, keepdims=True))   # shape (n_protein, 1)

    ratios = X / geo_mean                                  # (n_protein, n_sample)
    size_factor = np.median(ratios, axis=0)                # per-sample

    norm = X / size_factor
    out = df[["protein"]].copy()
    out[samples] = norm
    out.attrs["upregulated"] = df.attrs.get("upregulated", set())
    out.attrs["size_factors"] = dict(zip(samples, size_factor))
    return out


def per_protein_ttest(df: pd.DataFrame, controls: list, treated: list) -> pd.DataFrame:
    log_X = np.log2(np.maximum(df[controls + treated].to_numpy(), 1e-12))
    c_idx = list(range(len(controls)))
    t_idx = list(range(len(controls), len(controls) + len(treated)))

    log2fc = log_X[:, t_idx].mean(axis=1) - log_X[:, c_idx].mean(axis=1)

    tres = stats.ttest_ind(
        log_X[:, t_idx], log_X[:, c_idx], axis=1, equal_var=False,
    )
    p = np.where(np.isfinite(tres.pvalue), tres.pvalue, 1.0)
    p = np.clip(p, 1e-12, 1.0)  # cap resolution

    out = pd.DataFrame({
        "protein":    df["protein"].to_numpy(),
        "log2fc":     log2fc,
        "p_value":    p,
        "neglog10_p": -np.log10(p),
    })
    return out


intensity = simulate_intensity_table(proteins, rng)
intensity_norm = median_of_ratios(intensity)
stats_df = per_protein_ttest(intensity_norm, controls=SAMPLES_CTRL, treated=SAMPLES_TRT)

print("Size factors after median-of-ratios normalisation:")
for s, sf in intensity_norm.attrs["size_factors"].items():
    print(f"  {s}: {sf:.3f}")

# Top up-/down-regulated by |log2fc| with p < 0.05.
sig = stats_df[(stats_df["p_value"] < 0.05) & (stats_df["log2fc"].abs() > 1.0)]
print(f"\\nProteins with p<0.05 and |log2FC|>1.0: {len(sig)}  "
      f"(spiked-in upregulated: {len(intensity.attrs['upregulated'])})")

# ----- Volcano plot --------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
up_set = intensity.attrs["upregulated"]
mask_up = stats_df["protein"].isin(up_set)

ax.scatter(stats_df.loc[~mask_up, "log2fc"],
           stats_df.loc[~mask_up, "neglog10_p"],
           s=14, alpha=0.5, color="C7", label="other")
ax.scatter(stats_df.loc[mask_up, "log2fc"],
           stats_df.loc[mask_up, "neglog10_p"],
           s=30, alpha=0.9, color="C3", edgecolors="black", linewidths=0.5,
           label="spiked-in upregulated")
ax.axhline(-math.log10(0.05), color="grey", ls="--", lw=1)
ax.axvline( 1.0, color="grey", ls="--", lw=1)
ax.axvline(-1.0, color="grey", ls="--", lw=1)
ax.set_xlabel("log2 fold change  (treated / control)")
ax.set_ylabel(r"$-\\log_{10}$ p-value")
ax.set_title("Volcano plot — LFQ proteomics")
ax.legend(loc="upper left")
plt.tight_layout()
plt.show()
'''


# ---------------------------------------------------------------------------
# Step 5 — Missing values + isotopic labels discussion
# ---------------------------------------------------------------------------


STEP5_MD = """## Step 5 (15 min) — Missing values + SILAC / TMT alternatives

Two practical issues that hit every real-world LFQ pipeline:

- **Missing values.** Up to 30–50% of peptide-feature cells are blank in a
  real LFQ matrix — sometimes "missing at random" (MAR, instrument
  dropouts) and sometimes "missing not at random" (MNAR, true low
  abundance below detection limit). We simulate 30% MNAR missingness and
  compare three imputation strategies:
  - **Drop rows with any missing value** (baseline; you lose statistical
    power proportionally to dropout).
  - **Half-minimum imputation** — replace missing with `0.5 ×
    column-minimum`. Standard MNAR proxy.
  - **K-nearest-neighbours** on the protein direction. Standard MAR proxy.

- **Isotopic labelling.** SILAC encodes treated vs control in heavy/light
  amino acids; TMT in chemical tags on the N-terminus. Both reduce CV
  dramatically vs LFQ because every condition is measured in *the same
  MS1 scan*. We simulate the per-protein **coefficient of variation** for
  LFQ vs SILAC vs TMT to show the cost/precision tradeoff.

This step is short on TODOs and long on the comparison; tweak the
parameters to feel the trade-offs.
"""


STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Missing values + isotopic-label CV comparison.
# ----------------------------------------------------------------------


def inject_mnar(intensity: pd.DataFrame, frac: float, rng) -> pd.DataFrame:
    """Knock out `frac` of the lowest-intensity cells (MNAR pattern)."""
    # TODO:
    #   - Pick the `frac` smallest entries across the whole matrix (excluding
    #     the `protein` column) and set them to NaN.
    raise NotImplementedError


def impute_half_min(df: pd.DataFrame) -> pd.DataFrame:
    """Per-column half-minimum imputation."""
    # TODO
    raise NotImplementedError


def impute_knn(df: pd.DataFrame, k: int = 5) -> pd.DataFrame:
    """Naive KNN: for each missing cell, average the k nearest rows
    (by Euclidean distance over the non-missing columns) at that column."""
    # TODO
    raise NotImplementedError


def cv_per_protein(df: pd.DataFrame, samples: list) -> np.ndarray:
    """Coefficient of variation across `samples` for each row."""
    # TODO: std / mean over the selected columns.
    raise NotImplementedError


# Use these to:
#   1. Compute how many "spiked-in upregulated" proteins are recovered (q < 0.05)
#      under each imputation strategy.
#   2. Simulate per-protein CV under LFQ (sigma=0.30), SILAC (sigma=0.10),
#      TMT (sigma=0.05) and plot the distributions.
'''


STEP5_SOLUTION = '''# Reference solution — Step 5.


def inject_mnar(intensity: pd.DataFrame, frac: float, rng) -> pd.DataFrame:
    samples = [c for c in intensity.columns if c != "protein"]
    X = intensity[samples].to_numpy().astype(float).copy()

    flat = X.flatten()
    n_drop = int(frac * flat.size)
    threshold = np.partition(flat, n_drop)[n_drop - 1] if n_drop > 0 else -np.inf
    X[X <= threshold] = np.nan

    out = intensity[["protein"]].copy()
    out[samples] = X
    out.attrs.update(intensity.attrs)
    return out


def impute_half_min(df: pd.DataFrame) -> pd.DataFrame:
    samples = [c for c in df.columns if c != "protein"]
    X = df[samples].to_numpy().astype(float).copy()
    for j in range(X.shape[1]):
        col = X[:, j]
        col_min = np.nanmin(col)
        col[np.isnan(col)] = 0.5 * col_min
        X[:, j] = col
    out = df[["protein"]].copy()
    out[samples] = X
    out.attrs.update(df.attrs)
    return out


def impute_knn(df: pd.DataFrame, k: int = 5) -> pd.DataFrame:
    samples = [c for c in df.columns if c != "protein"]
    X = df[samples].to_numpy().astype(float).copy()
    n_row = X.shape[0]

    for i in range(n_row):
        if not np.any(np.isnan(X[i])):
            continue
        missing_cols = np.where(np.isnan(X[i]))[0]
        # Find rows with no NaN in the *observed* columns of row i.
        obs_cols = np.where(~np.isnan(X[i]))[0]
        if len(obs_cols) == 0:
            # Pathological — impute with global column means.
            for j in missing_cols:
                X[i, j] = np.nanmean(X[:, j])
            continue
        # Candidate rows: those fully observed at obs_cols *and* missing_cols.
        candidates = np.where(
            ~np.isnan(X[:, missing_cols]).any(axis=1)
            & ~np.isnan(X[:, obs_cols]).any(axis=1)
        )[0]
        candidates = candidates[candidates != i]
        if len(candidates) == 0:
            for j in missing_cols:
                X[i, j] = np.nanmean(X[:, j])
            continue
        # Euclidean distance over obs_cols (log-space for stability).
        dists = np.sqrt(np.sum(
            (np.log(np.maximum(X[candidates][:, obs_cols], 1e-12))
             - np.log(np.maximum(X[i, obs_cols], 1e-12))) ** 2,
            axis=1,
        ))
        nearest = candidates[np.argsort(dists)[:k]]
        for j in missing_cols:
            X[i, j] = np.mean(X[nearest, j])
    out = df[["protein"]].copy()
    out[samples] = X
    out.attrs.update(df.attrs)
    return out


def cv_per_protein(df: pd.DataFrame, samples: list) -> np.ndarray:
    X = df[samples].to_numpy().astype(float)
    mu = np.nanmean(X, axis=1)
    sd = np.nanstd(X, axis=1)
    return sd / np.maximum(mu, 1e-12)


# ---- Compare imputation strategies on recovery of upregulated proteins ----

def recovery_after_imputation(intensity_missing, strategy_name, imputed) -> int:
    imputed_norm = median_of_ratios(imputed)
    sdf = per_protein_ttest(imputed_norm, controls=SAMPLES_CTRL, treated=SAMPLES_TRT)
    up_set = intensity.attrs["upregulated"]
    sig_up = sdf[
        (sdf["protein"].isin(up_set))
        & (sdf["p_value"] < 0.05)
        & (sdf["log2fc"] > 1.0)
    ]
    print(f"  {strategy_name:<20s}: recovered {len(sig_up):2d}/{len(up_set)} upregulated")
    return len(sig_up)


print("Injecting 30% MNAR missingness...")
intensity_missing = inject_mnar(intensity, frac=0.30, rng=np.random.default_rng(SEED + 1))
n_missing = intensity_missing[ALL_SAMPLES].isna().sum().sum()
print(f"  Total missing cells: {n_missing} / {intensity_missing[ALL_SAMPLES].size}\\n")

# Strategy 1: drop incomplete rows.
intensity_drop = intensity_missing.dropna(subset=ALL_SAMPLES).copy()
intensity_drop.attrs.update(intensity_missing.attrs)
print(f"  drop-incomplete-rows : {len(intensity_drop):3d} proteins remain")
recovery_after_imputation(intensity_missing, "drop-incomplete", intensity_drop)

# Strategy 2: half-minimum imputation.
recovery_after_imputation(
    intensity_missing,
    "half-min impute",
    impute_half_min(intensity_missing),
)

# Strategy 3: KNN imputation.
recovery_after_imputation(
    intensity_missing,
    "KNN impute (k=5)",
    impute_knn(intensity_missing, k=5),
)


# ---- CV comparison: LFQ vs SILAC vs TMT -----------------------------------
print("\\nSimulating per-protein CV under three labelling strategies...")

def sim_cv(sigma_noise, n=len(proteins), seed=0):
    r = np.random.default_rng(seed)
    base = r.lognormal(mean=12.0, sigma=1.0, size=n)
    noise = r.lognormal(mean=0.0, sigma=sigma_noise, size=(n, 6))
    X = base[:, None] * noise
    return np.std(X, axis=1) / np.mean(X, axis=1)

cv_lfq   = sim_cv(0.30, seed=SEED + 10)
cv_silac = sim_cv(0.10, seed=SEED + 11)
cv_tmt   = sim_cv(0.05, seed=SEED + 12)

print(f"  median CV  LFQ   = {np.median(cv_lfq)*100:5.1f}%")
print(f"  median CV  SILAC = {np.median(cv_silac)*100:5.1f}%")
print(f"  median CV  TMT   = {np.median(cv_tmt)*100:5.1f}%")

fig, ax = plt.subplots(figsize=(8, 4))
bins = np.linspace(0, 0.6, 40)
ax.hist(cv_lfq,   bins=bins, alpha=0.6, label="LFQ   (σ_noise=0.30)")
ax.hist(cv_silac, bins=bins, alpha=0.6, label="SILAC (σ_noise=0.10)")
ax.hist(cv_tmt,   bins=bins, alpha=0.6, label="TMT   (σ_noise=0.05)")
ax.set_xlabel("coefficient of variation across 6 samples")
ax.set_ylabel("# proteins")
ax.set_title("Quantification precision by labelling strategy")
ax.legend()
plt.tight_layout()
plt.show()
'''


# ---------------------------------------------------------------------------
# Self-check
# ---------------------------------------------------------------------------


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces of the pipeline.
"""


SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Total PSM count and target/decoy balance.
assert len(psms) == 20_000, f"Expected 20000 PSMs, got {len(psms)}"
assert (psms["label"] == "target").sum() == 10_000, "Target count wrong"
assert (psms["label"] == "decoy").sum()  == 10_000, "Decoy count wrong"

# 2. Mean XCorr score: targets clearly higher than decoys (the whole point).
mean_t = psms.loc[psms["label"] == "target", "xcorr"].mean()
mean_d = psms.loc[psms["label"] == "decoy",  "xcorr"].mean()
assert mean_t > mean_d + 5, (
    f"Target mean XCorr ({mean_t:.2f}) should clearly exceed decoy mean ({mean_d:.2f})"
)

# 3. 1% FDR cutoff: 4000-6500 target PSMs should survive given the 60/40 mixture.
n_pass = int(((psms["label"] == "target") & (psms["xcorr"] >= t1pct)).sum())
print(f"PSMs surviving 1% FDR: {n_pass}")
assert 4_000 <= n_pass <= 6_500, (
    f"PSMs at 1% FDR is {n_pass}, expected ~4000-6500 (the TRUE_FRAC=0.6 should pass cleanly)"
)

# 4. q-value is monotone non-decreasing as you walk down the ranked PSM list
#    (because at each step we add either a target — which can lower fdr — or a
#    decoy — which raises it; the running-min-from-the-right enforces this).
diffs = np.diff(sweep["q_value"].to_numpy())
assert np.all(diffs >= -1e-9), (
    f"q-value is not monotone non-decreasing as expected; min diff = {diffs.min()}"
)

# 5. Inferred proteins: between 50 and 100 with >= 2 unique peptides.
assert 50 <= len(proteins) <= 100, (
    f"Inferred {len(proteins)} proteins; expected 50-100 from 100 simulated"
)
assert (proteins["n_unique_peptides"] >= 2).all(), "Min-2-unique filter violated"

# 6. Volcano: at least 15 of the 20 spiked-in proteins should land in the
#    "significant + up" quadrant (p < 0.05 and log2fc > 1).
up_set = intensity.attrs["upregulated"]
recovered = stats_df[
    (stats_df["protein"].isin(up_set))
    & (stats_df["p_value"] < 0.05)
    & (stats_df["log2fc"] > 1.0)
]
print(f"Recovered spiked-in upregulated: {len(recovered)}/{len(up_set)}")
assert len(recovered) >= 15, (
    f"Volcano recovered only {len(recovered)}/{len(up_set)} upregulated proteins; "
    "check normalisation + t-test."
)

# 7. CV ordering: TMT < SILAC < LFQ. If you skipped Step 5's CV comparison,
#    this assert will be skipped.
try:
    assert np.median(cv_tmt) < np.median(cv_silac) < np.median(cv_lfq), (
        "Expected median CV ordering TMT < SILAC < LFQ"
    )
except NameError:
    print("(skipping CV ordering check — Step 5 not run)")

print("✅ Self-check passed.")
'''


# ---------------------------------------------------------------------------
# Closing EE framing
# ---------------------------------------------------------------------------


EE_MD = """## EE framing — empirical-null detection at scale

You implemented the canonical proteomics statistics pipeline. Three EE
ideas thread the whole thing:

1. **Matched-filter scoring in $m/z$ space.** XCorr (and SEQUEST-style
   scores in general) is a normalised dot product between an observed
   MS2 spectrum and the theoretical b/y ions of a candidate peptide.
   That is a one-dimensional matched filter, and the search engine just
   runs the filter bank across every tryptic peptide in the database.

2. **Decoys as an empirical null distribution.** Instead of assuming
   the score-under-null follows some parametric form (Gumbel, Weibull,
   χ²-ish, etc.), you generate decoys that match the target null *by
   construction*. Decoy-count / target-count is then an unbiased
   estimator of the local FDR. This is the proteomics version of
   constant-false-alarm-rate (CFAR) detection — calibrate the threshold
   off measured noise, not theory.

3. **Multiple-testing the right way.** With $10^5$ parallel hypothesis
   tests per LC-MS run, controlling FWER (Bonferroni) is suicidal. FDR
   / q-values trade a quantifiable false-positive rate for *vastly*
   more power; the volcano plot lives in this regime. Same logic in
   genomics (Benjamini-Hochberg on RNA-seq tests) and in any other
   "tens-of-thousands-of-tests" pipeline.

4. **Labelling vs LFQ is signal-conditioning vs raw acquisition.**
   SILAC and TMT modulate every condition into a single MS1 scan
   (think DSB → coherent demodulation), so the per-condition noise is
   common-mode and cancels in ratios. LFQ measures each condition
   independently (think envelope detection on separate channels) and
   suffers run-to-run gain drift — higher CV, more missing values, but
   no labelling cost. Pick your noise model deliberately.
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
        md(SOL_HEADER),
        hidden(STEP1_SOLUTION),

        md(STEP2_MD),
        code(STEP2_TODO),
        md(SOL_HEADER),
        hidden(STEP2_SOLUTION),

        md(STEP3_MD),
        code(STEP3_TODO),
        md(SOL_HEADER),
        hidden(STEP3_SOLUTION),

        md(STEP4_MD),
        code(STEP4_TODO),
        md(SOL_HEADER),
        hidden(STEP4_SOLUTION),

        md(STEP5_MD),
        code(STEP5_TODO),
        md(SOL_HEADER),
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
