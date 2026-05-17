"""Build exercise.ipynb for L18 — Cancer Genomics: Integrated Capstone.

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


TITLE_MD = """# L18 — Cancer Genomics: Integrated Capstone

In this exercise you decompose a synthetic 50-tumour cohort of 96-bin
mutational spectra into its underlying mutagenic processes. You implement
the Lee-Seung multiplicative-update **non-negative matrix factorisation
(NMF)** algorithm, recover SBS1 (clock), SBS4 (tobacco), SBS7a (UV), and
SBS3 (HRD) up to permutation, and sweep the model order `K` to see
under-fit vs over-fit behaviour. The aha is that mutational signatures
are exactly **non-negative blind source separation** — the same maths
used for hyperspectral unmixing and topic modelling.
"""


AHA_MD = """> **Aha.** A tumour's 96-bin mutation spectrum is a non-negative
> linear combination of a handful of mutagenic-process signatures.
> Recovering the signatures **and** the per-tumour exposures from the
> raw spectra alone is the textbook **NMF** problem: positivity-constrained
> matrix factorisation `V approximately = W H` solvable by
> multiplicative-update projected gradient descent. Pick `K` too small
> and recovered components blur multiple signatures together; pick `K`
> too large and the extra components fit noise.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scikit-learn==1.5.0 matplotlib==3.8.4 -q
"""


IMPORTS = """import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import permutations

# Deterministic for the whole notebook.
np.random.seed(42)

# Canonical 96 trinucleotide-substitution context labels.
# Order: 6 substitution types x 16 (5'-base, 3'-base) contexts.
SUB_TYPES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
BASES = ["A", "C", "G", "T"]


def context_labels():
    \"\"\"Return the 96 standard trinucleotide-context labels (e.g. 'A[C>A]A').\"\"\"
    labels = []
    for sub in SUB_TYPES:
        ref = sub[0]
        for left in BASES:
            for right in BASES:
                labels.append(f"{left}[{sub}]{right}")
    return labels


CONTEXTS = context_labels()
assert len(CONTEXTS) == 96
"""


STEP1_MD = """## Step 1 (8 min) — Build the 50-sample cohort and barplot one spectrum

The four canonical signatures we will recover are simplified, normalised
approximations of the published COSMIC v3.4 profiles for SBS1, SBS4,
SBS7a, and SBS3. Each is a length-96 probability vector over the standard
trinucleotide-substitution alphabet. We mix them with random non-negative
exposures plus Poisson sampling noise to synthesise 50 tumour spectra,
each summing to ~1500 mutations. The mixing matrix `H_true` (4 x 50)
and the basis `W_true` (96 x 4) are the **ground truth** we will try
to recover blind.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Build COSMIC-like signatures, simulate a 50-sample cohort.
# ----------------------------------------------------------------------

# Simplified COSMIC v3.4 SBS1/4/7a/3 profiles (96-bin, normalised).
# Order: 6 substitution types (C>A, C>G, C>T, T>A, T>C, T>G) x 16 contexts.
SIG_NAMES = ["SBS1", "SBS4", "SBS7a", "SBS3"]

def build_signatures() -> np.ndarray:
    """Return a (96, 4) matrix whose columns are SBS1, SBS4, SBS7a, SBS3."""
    # TODO: hardcode the four 96-bin patterns (see hidden solution for the
    # canonical-shape numbers); normalise each column to sum to 1.
    raise NotImplementedError


def simulate_cohort(W_true: np.ndarray, n_samples: int = 50,
                    mean_count: int = 1500, seed: int = 42) -> tuple[np.ndarray, np.ndarray]:
    """Return (V, H_true): V is (96, n_samples) Poisson-noised count matrix,
    H_true is (4, n_samples) ground-truth exposures (columns sum to 1)."""
    # TODO:
    # 1. Draw H_true ~ Dirichlet over 4 signatures (so each sample is a
    #    non-negative mixture summing to 1).
    # 2. Total mutation count per sample: Poisson(mean_count).
    # 3. V_noisefree = W_true @ (H_true * total_counts).
    # 4. V = Poisson(V_noisefree).
    raise NotImplementedError


# W_true = build_signatures()
# V, H_true = simulate_cohort(W_true)
# print("W_true shape:", W_true.shape, "  V shape:", V.shape)

# Barplot one sample's spectrum, coloured by substitution type.
'''


STEP1_SOLUTION_HEADER = """*Click to expand the reference solution.*"""


STEP1_SOLUTION = '''# Reference solution -- Step 1.
# Simplified COSMIC v3.4 SBS1/4/7a/3 profiles; shapes follow Alexandrov et al.
# and the COSMIC public release. Each pattern is a length-16 vector over the
# (5\'-base, 3\'-base) context grid for one substitution type.

SIG_NAMES = ["SBS1", "SBS4", "SBS7a", "SBS3"]


def _make_sig(per_sub: list[list[float]]) -> np.ndarray:
    """Stack 6 length-16 patterns into a 96-vector and normalise to sum=1."""
    sig = np.array(per_sub, dtype=float).ravel()
    return sig / sig.sum()


def build_signatures() -> np.ndarray:
    # SBS1 -- clock-like; spontaneous deamination of 5mC -> C>T at NCG.
    # The NCG trinucleotide contexts sit at positions [2, 6, 10, 14] within
    # each substitution\'s 16-context block (5\'-base x 3\'-base lexical order).
    low = [0.3 + 0.1 * (j % 3) for j in range(16)]
    sbs1_ct = [0.5, 0.6, 8.0, 0.6,  0.7, 0.6, 6.5, 0.7,
               0.5, 0.7, 5.0, 0.7,  0.5, 0.7, 7.0, 0.7]
    SBS1 = _make_sig([low, low, sbs1_ct, low, low, low])

    # SBS4 -- tobacco smoking; broad C>A peaks.
    sbs4_ca = [3.0, 4.0, 1.5, 4.5,  5.5, 6.0, 2.5, 6.5,
               4.5, 5.0, 2.0, 5.5,  3.5, 4.5, 2.0, 4.5]
    flat = lambda h: [h + 0.4 * (j % 4) for j in range(16)]
    SBS4 = _make_sig([sbs4_ca, flat(0.6), flat(0.8), flat(0.6), flat(0.9), flat(0.5)])

    # SBS7a -- UV; high C>T at CC and TC dimers (positions 4-7 = CCN; 12-15 = TCN).
    sbs7_ct = [0.6, 0.5, 0.7, 0.6,  5.0, 6.5, 5.2, 4.5,
               0.6, 0.7, 0.8, 0.6,  4.8, 5.8, 5.0, 4.0]
    SBS7a = _make_sig([low, low, sbs7_ct, low, low, low])

    # SBS3 -- HRD; flat / diffuse across all contexts.
    flat3 = lambda h: [h + 0.2 * np.sin(j * 0.7) for j in range(16)]
    SBS3 = _make_sig([flat3(1.0), flat3(1.0), flat3(1.4),
                      flat3(0.9), flat3(1.2), flat3(1.0)])

    return np.stack([SBS1, SBS4, SBS7a, SBS3], axis=1)  # (96, 4)


def simulate_cohort(W_true: np.ndarray, n_samples: int = 50,
                    mean_count: int = 1500, seed: int = 42):
    rng = np.random.default_rng(seed)
    K = W_true.shape[1]
    # Dirichlet exposures -- each column is a non-negative mixture summing to 1.
    # Use uneven concentrations so some signatures dominate some samples.
    alpha = np.array([0.7, 0.7, 0.7, 0.7])
    H_true = rng.dirichlet(alpha, size=n_samples).T  # (K, n_samples)
    totals = rng.poisson(mean_count, size=n_samples).astype(float)
    V_noisefree = W_true @ (H_true * totals[None, :])
    V = rng.poisson(V_noisefree).astype(float)
    return V, H_true


W_true = build_signatures()
V, H_true = simulate_cohort(W_true, n_samples=50, mean_count=1500)
print(f"W_true shape: {W_true.shape}   V shape: {V.shape}")
print(f"per-sample totals: min={V.sum(0).min():.0f}  median={np.median(V.sum(0)):.0f}  max={V.sum(0).max():.0f}")
print(f"H_true column sums (should all be 1): min={H_true.sum(0).min():.3f} max={H_true.sum(0).max():.3f}")

# Barplot the first sample\'s spectrum, coloured by substitution type.
SUB_COLORS = ["#1e40af", "#b45309", "#c4342c", "#525252", "#0d7377", "#7c3aed"]
fig, ax = plt.subplots(figsize=(11, 3.5))
bins = np.arange(96)
for s, color in enumerate(SUB_COLORS):
    sl = slice(s * 16, (s + 1) * 16)
    ax.bar(bins[sl], V[sl, 0], color=color, label=SUB_TYPES[s], width=0.9)
# Light vertical separators between substitution-type blocks.
for s in range(1, 6):
    ax.axvline(s * 16 - 0.5, color="black", alpha=0.1, lw=0.5)
ax.set_xlim(-0.5, 95.5)
ax.set_xlabel("trinucleotide context (96 bins)")
ax.set_ylabel("mutation count")
ax.set_title(f"Sample 0 spectrum -- dominant exposures: " +
             ", ".join(f"{n}={H_true[i, 0]:.2f}" for i, n in enumerate(SIG_NAMES)))
ax.legend(ncol=6, loc="upper right", fontsize=8)
plt.tight_layout()
plt.show()

# Label the top trinucleotide contexts in sample 0.
top_idx = np.argsort(V[:, 0])[::-1][:5]
print("\\nTop-5 trinucleotide bins in sample 0:")
for i in top_idx:
    print(f"  {CONTEXTS[i]:>10s}  count={V[i, 0]:6.0f}")
'''


STEP2_MD = """## Step 2 (12 min) — Implement Lee-Seung multiplicative-update NMF

Lee & Seung (1999, 2001) showed that **multiplicative updates** preserve
non-negativity for free and monotonically decrease the Frobenius error
`||V - W H||_F^2`. The update rules are:

$$ H \\leftarrow H \\odot \\frac{W^\\top V}{W^\\top W H} \\qquad
   W \\leftarrow W \\odot \\frac{V H^\\top}{W H H^\\top} $$

where `odot` is element-wise multiplication and division is element-wise.
A small epsilon in the denominators prevents division-by-zero. This is
**projected gradient descent** with a clever step size that keeps `W` and
`H` non-negative without an explicit projection step.
"""


STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 -- multiplicative-update NMF.
# ----------------------------------------------------------------------

EPS = 1e-10


def nmf_lee_seung(V: np.ndarray, K: int, n_iter: int = 200,
                  seed: int = 0) -> tuple[np.ndarray, np.ndarray, list[float]]:
    """Lee-Seung multiplicative-update NMF.

    Returns (W, H, errors) where W is (n_features, K), H is (K, n_samples),
    and errors is the Frobenius reconstruction error per iteration.
    """
    # TODO:
    # 1. Random non-negative init of W (shape (96, K)) and H (shape (K, n_samples)).
    # 2. For n_iter iterations:
    #     H *= (W.T @ V) / (W.T @ W @ H + EPS)
    #     W *= (V @ H.T) / (W @ H @ H.T + EPS)
    #     record ||V - W H||_F
    # 3. Normalise columns of W to sum to 1 (and absorb the scale into H) at the end.
    raise NotImplementedError


# W_hat, H_hat, errs = nmf_lee_seung(V, K=4, n_iter=200, seed=0)
'''


STEP2_SOLUTION = '''# Reference solution -- Step 2.

EPS = 1e-10


def nmf_lee_seung(V: np.ndarray, K: int, n_iter: int = 200,
                  seed: int = 0) -> tuple[np.ndarray, np.ndarray, list[float]]:
    rng = np.random.default_rng(seed)
    n_features, n_samples = V.shape
    # Init scaled so that W @ H has the same order of magnitude as V.
    scale = np.sqrt(V.mean() / K)
    W = rng.uniform(EPS, 1.0, size=(n_features, K)) * scale
    H = rng.uniform(EPS, 1.0, size=(K, n_samples)) * scale

    errors = []
    for it in range(n_iter):
        # Update H first, then W (any order is fine; this matches sklearn).
        WtV  = W.T @ V
        WtWH = W.T @ W @ H + EPS
        H *= WtV / WtWH

        VHt  = V @ H.T
        WHHt = W @ H @ H.T + EPS
        W *= VHt / WHHt

        err = np.linalg.norm(V - W @ H, "fro")
        errors.append(err)

    # Normalise columns of W to sum to 1 (signatures are probability vectors);
    # absorb the scale into H so V \\approx W @ H still holds.
    col_sums = W.sum(axis=0, keepdims=True) + EPS
    W = W / col_sums
    H = H * col_sums.T

    return W, H, errors


W_hat, H_hat, errs = nmf_lee_seung(V, K=4, n_iter=200, seed=0)
print(f"Final Frobenius error: {errs[-1]:.2f}   (initial: {errs[0]:.2f})")

# Plot the error curve -- should drop fast then plateau.
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(errs, lw=1.5)
ax.set_xlabel("iteration")
ax.set_ylabel("||V - W H||_F")
ax.set_title("Lee-Seung NMF reconstruction error vs iteration")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()
'''


STEP3_MD = """## Step 3 (15 min) — Cosine-similarity match to ground truth

NMF only recovers `W` and `H` up to a permutation (and a non-negative
rescaling — we already fixed the scale by normalising `W`'s columns to
sum to 1). To grade the fit we compute the **cosine similarity** between
each recovered column and each true signature, then find the assignment
that maximises the mean similarity. With four signatures there are
only `4! = 24` permutations, so brute-force is fine.
"""


STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 -- cosine similarity + optimal permutation match.
# ----------------------------------------------------------------------


def cosine_sim(a: np.ndarray, b: np.ndarray) -> float:
    """Cosine similarity between two non-negative vectors."""
    # TODO
    raise NotImplementedError


def match_signatures(W_hat: np.ndarray, W_true: np.ndarray) -> tuple[tuple[int, ...], np.ndarray]:
    """Return (best_perm, per_signature_cosines).

    best_perm[i] is the index in W_hat that best matches W_true column i.
    """
    # TODO: enumerate all permutations of range(K) and pick the one that
    # maximises the average cosine sim of (W_hat[:, perm[i]], W_true[:, i]).
    raise NotImplementedError


# perm, sims = match_signatures(W_hat, W_true)
# print("Per-signature cosines:", dict(zip(SIG_NAMES, sims)))
'''


STEP3_SOLUTION = '''# Reference solution -- Step 3.


def cosine_sim(a: np.ndarray, b: np.ndarray) -> float:
    na = np.linalg.norm(a) + EPS
    nb = np.linalg.norm(b) + EPS
    return float(a @ b / (na * nb))


def match_signatures(W_hat: np.ndarray, W_true: np.ndarray):
    K = W_true.shape[1]
    best_perm, best_mean = None, -1.0
    for perm in permutations(range(K)):
        sims = [cosine_sim(W_hat[:, perm[i]], W_true[:, i]) for i in range(K)]
        m = float(np.mean(sims))
        if m > best_mean:
            best_mean, best_perm = m, perm
    sims = np.array([cosine_sim(W_hat[:, best_perm[i]], W_true[:, i]) for i in range(K)])
    return best_perm, sims


perm, sims = match_signatures(W_hat, W_true)
print("Best permutation (truth -> recovered column index):", perm)
print("Per-signature cosine similarities:")
for name, s in zip(SIG_NAMES, sims):
    print(f"  {name:>6s}: {s:.3f}")
print(f"Mean cosine similarity: {sims.mean():.3f}")

# Side-by-side barplots: true (left) vs recovered (right) for each signature.
fig, axes = plt.subplots(4, 2, figsize=(13, 9), sharex=True)
for i, name in enumerate(SIG_NAMES):
    j = perm[i]
    for s, color in enumerate(SUB_COLORS):
        sl = slice(s * 16, (s + 1) * 16)
        axes[i, 0].bar(np.arange(16) + s * 16, W_true[sl, i], color=color, width=0.9)
        axes[i, 1].bar(np.arange(16) + s * 16, W_hat[sl, j],  color=color, width=0.9)
    axes[i, 0].set_ylabel(name)
    axes[i, 0].set_title(f"truth -- {name}" if i == 0 else "")
    axes[i, 1].set_title(f"recovered (cos={sims[i]:.3f})" if i == 0 else f"cos={sims[i]:.3f}")
for ax in axes[-1, :]:
    ax.set_xlabel("trinucleotide context (96 bins)")
plt.tight_layout()
plt.show()
'''


STEP4_MD = """## Step 4 (15 min) — Model-order selection: sweep K in {2, 3, 4, 5}

How would you choose `K = 4` if you didn't know the truth? The classic
NMF / PCA / filter-order trick: plot the **reconstruction error** vs `K`
and look for the elbow. Under-fitting (`K = 2` or `3`) leaves a lot
of variance unexplained; over-fitting (`K = 5`) only marginally
improves the error because the fifth component soaks up Poisson noise.

We track three quantities per `K`:

1. **Frobenius reconstruction error** -- should drop sharply up to the
   true `K`, then plateau.
2. **Mean cosine similarity to the true signatures** -- only defined
   when at least 4 components are available; we report the best-match
   subset.
3. **Component stability** -- check whether the recovered fifth column
   at `K = 5` is *interpretable* or just noise.
"""


STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 -- model-order sweep.
# ----------------------------------------------------------------------

K_values = [2, 3, 4, 5]

# TODO:
# 1. For each K, run nmf_lee_seung(V, K, n_iter=200, seed=0).
# 2. Record final Frobenius error.
# 3. For K >= 4, compute mean cosine sim to W_true (only the top-4 matching
#    columns count -- enumerate combinations(range(K), 4) and pick the best).
# 4. Plot error vs K; annotate the elbow at K=4.
'''


STEP4_SOLUTION = '''# Reference solution -- Step 4.
from itertools import combinations

K_values = [2, 3, 4, 5]

errs_by_K = {}
sims_by_K = {}
W_by_K = {}
for K in K_values:
    W_k, H_k, errs_k = nmf_lee_seung(V, K=K, n_iter=200, seed=0)
    errs_by_K[K] = errs_k[-1]
    W_by_K[K] = W_k

    if K >= 4:
        # Brute force: pick which 4 of the K recovered columns best explain the truth.
        best_mean = -1.0
        best_subset_sims = None
        for subset in combinations(range(K), 4):
            sub_W = W_k[:, list(subset)]
            _, s = match_signatures(sub_W, W_true)
            if s.mean() > best_mean:
                best_mean = s.mean()
                best_subset_sims = s
        sims_by_K[K] = best_subset_sims
    else:
        sims_by_K[K] = None

# Summary table.
print(f"{'K':>3s}  {'Frob. error':>14s}  {'mean cos sim (top-4)':>22s}")
for K in K_values:
    s = sims_by_K[K]
    s_str = f"{s.mean():.3f}" if s is not None else "n/a (K < 4)"
    print(f"{K:>3d}  {errs_by_K[K]:>14.2f}  {s_str:>22s}")

# Plot reconstruction error vs K with the elbow at K=4.
fig, axes = plt.subplots(1, 2, figsize=(13, 4))

ax = axes[0]
ax.plot(K_values, [errs_by_K[K] for K in K_values], "o-", lw=2, ms=8)
ax.axvline(4, color="green", linestyle="--", alpha=0.5, label="truth K=4")
ax.set_xlabel("model order K")
ax.set_ylabel("final ||V - W H||_F")
ax.set_title("Reconstruction error vs K (elbow at K=4)")
ax.set_xticks(K_values)
ax.legend()
ax.grid(alpha=0.3)

# Plot the spurious 5th signature at K=5 -- expect noise / split-signature.
ax = axes[1]
W5 = W_by_K[5]
# Identify the 5th column (the one not picked by the best-matching subset).
best_subset = None
best_mean = -1.0
for subset in combinations(range(5), 4):
    sub_W = W5[:, list(subset)]
    _, s = match_signatures(sub_W, W_true)
    if s.mean() > best_mean:
        best_mean = s.mean()
        best_subset = subset
spurious = [c for c in range(5) if c not in best_subset][0]
for s, color in enumerate(SUB_COLORS):
    sl = slice(s * 16, (s + 1) * 16)
    ax.bar(np.arange(16) + s * 16, W5[sl, spurious], color=color, width=0.9)
ax.set_xlabel("trinucleotide context (96 bins)")
ax.set_ylabel("normalised intensity")
ax.set_title(f"K=5 spurious 5th component (col {spurious}) -- soaks up Poisson noise")
plt.tight_layout()
plt.show()
'''


STEP5_MD = """## Step 5 (10 min) — Therapeutic interpretation per sample

Now that you have recovered exposures `H_hat`, you can stratify patients
by their dominant mutagenic process and translate that into a candidate
treatment hypothesis. From the lecture:

- **SBS3 (HRD)** dominant -> PARP-inhibitor candidate (olaparib, niraparib).
- **SBS7a (UV)** dominant -> consistent with melanoma; mostly diagnostic.
- **SBS4 (tobacco)** dominant -> smoking-related; consistent with NSCLC.
- **SBS1 (clock)** dominant -> age-related background; no specific therapy.

You will tag every one of the 50 simulated tumours with its top
signature (using the recovered, label-permuted `H_hat`) and identify
the PARP-inhibitor cohort.
"""


STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 -- interpret recovered exposures.
# ----------------------------------------------------------------------

# TODO:
# 1. Permute the rows of H_hat (from the K=4 run) to match the truth order
#    using `perm` from Step 3.
# 2. Renormalise each column of H_hat to sum to 1 -- so each sample\'s
#    exposures are a probability distribution over signatures.
# 3. For each sample, take argmax over signatures -> dominant_signature[i].
# 4. Tabulate sample, dominant signature, candidate therapy.
# 5. Compare the recovered dominant-signature labels to the truth from H_true
#    (also argmax-reduced): report the agreement rate.
'''


STEP5_SOLUTION = '''# Reference solution -- Step 5.

# Permute recovered H rows to match SIG_NAMES order.
H_hat_aligned = H_hat[list(perm), :]
# Renormalise each column so it sums to 1 -- exposures as a probability.
H_hat_norm = H_hat_aligned / (H_hat_aligned.sum(axis=0, keepdims=True) + EPS)

dominant_recovered = np.argmax(H_hat_norm, axis=0)
dominant_truth     = np.argmax(H_true, axis=0)
agreement = (dominant_recovered == dominant_truth).mean()
print(f"Dominant-signature agreement (recovered vs truth): {agreement * 100:.0f}% ({(dominant_recovered == dominant_truth).sum()}/50)")

THERAPY = {
    "SBS1":  "no specific therapy (age-related background)",
    "SBS4":  "smoking-related; standard NSCLC workup",
    "SBS7a": "UV; melanoma workup (no signature-specific therapy)",
    "SBS3":  "PARP inhibitor candidate (olaparib / niraparib)",
}

rows = []
for i in range(V.shape[1]):
    rows.append({
        "sample":      i,
        "dominant":    SIG_NAMES[dominant_recovered[i]],
        "exposure":    round(float(H_hat_norm[dominant_recovered[i], i]), 2),
        "truth":       SIG_NAMES[dominant_truth[i]],
        "candidate":   THERAPY[SIG_NAMES[dominant_recovered[i]]],
    })
cohort = pd.DataFrame(rows)
print()
print("First 10 samples:")
print(cohort.head(10).to_string(index=False))

# How many PARP candidates did we identify?
parp = cohort[cohort["dominant"] == "SBS3"]
print(f"\\nPARP-inhibitor candidates (SBS3 dominant): {len(parp)} / 50")
print(parp[["sample", "exposure", "truth"]].to_string(index=False))

# Stacked bar of recovered exposures across the cohort.
fig, ax = plt.subplots(figsize=(11, 3.8))
order = np.argsort(dominant_recovered)  # sort samples by dominant signature
bottom = np.zeros(V.shape[1])
for i, name in enumerate(SIG_NAMES):
    color = ["#c4342c", "#1e40af", "#b45309", "#0d7377"][i]
    ax.bar(np.arange(V.shape[1]), H_hat_norm[i, order], bottom=bottom, color=color, label=name)
    bottom += H_hat_norm[i, order]
ax.set_xlabel("sample (sorted by dominant signature)")
ax.set_ylabel("recovered exposure")
ax.set_title("Recovered per-sample signature exposures (cohort, n=50)")
ax.legend(ncol=4)
plt.tight_layout()
plt.show()
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces of the
decomposition. If you ran the reference solutions above they should all
pass; if you wrote your own and an assert trips, revisit the
corresponding step.
"""


SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check -- runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. The four signatures are 96-long probability vectors.
assert W_true.shape == (96, 4)
for k in range(4):
    assert abs(W_true[:, k].sum() - 1.0) < 1e-9, f"signature {k} not normalised"

# 2. The cohort is non-negative, 96 x 50, and column sums are positive.
assert V.shape == (96, 50)
assert (V >= 0).all()
assert (V.sum(axis=0) > 0).all()

# 3. Reconstruction error decreases monotonically (multiplicative-update guarantee).
assert errs[0] > errs[-1], "NMF error did not decrease"
# Strict monotonicity isn\'t guaranteed every iteration due to floating-point;
# allow tiny upticks but require an overall ~10x drop.
assert errs[-1] < 0.5 * errs[0], "NMF error barely decreased"

# 4. All four signatures are recovered with cosine similarity >= 0.85.
for name, s in zip(SIG_NAMES, sims):
    assert s >= 0.85, f"recovered {name} cosine = {s:.3f} < 0.85 -- decomposition failed"

# 5. Mean cosine similarity rises sharply from K=3 to K=4 (the true model order).
# At K=3 we can\'t cover four signatures, so we expect a much lower mean
# (computed by enumerating combinations of 3 truth columns vs 3 recovered).
def _mean_sim_at_K(K):
    s = sims_by_K[K]
    return s.mean() if s is not None else None

assert _mean_sim_at_K(4) >= 0.90, f"K=4 mean cosine = {_mean_sim_at_K(4):.3f}; expected >= 0.90"

# 6. Dominant-signature agreement is at least 80%.
assert agreement >= 0.80, f"dominant-signature agreement only {agreement * 100:.0f}%"

print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — non-negative blind source separation

You implemented the canonical **non-negative blind source separation**
pipeline. The mathematics is identical to:

- **Hyperspectral unmixing.** Pixels are mixtures of pure-material
  spectra; spectra are non-negative; abundances are non-negative.
- **Speech separation under non-negative constraints.** Magnitude
  spectrograms of co-occurring speakers add up rather than cancel.
- **Topic modelling.** A document is a non-negative mixture of topics;
  each topic is a non-negative distribution over words.

The Lee-Seung multiplicative-update rule is a **projected gradient
descent** with a step size chosen so non-negativity is preserved without
an explicit projection. The convergence guarantee uses an auxiliary
function (a tight upper bound on the loss); the proof is in the original
1999 NIPS paper.

**Model-order selection** is the same problem you meet in PCA (`k`),
linear-prediction filter order (`p`), Gaussian mixture components (`K`),
or radar target enumeration: pick the smallest model that explains the
data well, and watch the residual elbow. In cancer genomics the
production tool (SigProfilerExtractor) wraps this with a stability
analysis -- run NMF many times with different inits and keep only `K`
values whose components recur across runs.

The signal-processing intuition is exact, not metaphorical: a tumour's
mutation spectrum is the sum of mutagenic-process "sources" weighted by
their activity in that tumour, and recovering both the sources and the
weights from many noisy spectra is exactly the problem NMF was designed
to solve.
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
