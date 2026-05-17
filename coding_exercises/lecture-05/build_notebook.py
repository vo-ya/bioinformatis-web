"""Build exercise.ipynb for L05 — Bulk RNA-seq: From Reads to Transcript Abundances.

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


TITLE_MD = """# L05 — Bulk RNA-seq: From Reads to Transcript Abundances

In this exercise you build the **compatibility-class + EM** core of a
modern RNA-seq quantifier (Kallisto / Salmon / RSEM lineage). Three
overlapping isoforms share exons, so most reads are **ambiguous** — they
are compatible with more than one transcript. You will index isoform
k-mers, classify each read into its compatibility class, and run
Expectation-Maximization to recover the ground-truth abundances.
"""


AHA_MD = """> **Aha.** Ambiguous reads are not garbage. **Soft EM assignments**
> use the current abundance estimate as a prior, and within a handful of
> iterations the fractional read weights converge on the true mixture.
> Hard-assigning each ambiguous read to its first match would throw away
> the very signal that distinguishes overlapping isoforms.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 matplotlib==3.8.4 -q
"""


IMPORTS = '''import math
import random
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Deterministic for the whole notebook.
SEED = 42
np.random.seed(SEED)
random.seed(SEED)

# DNA alphabet for the simulated isoforms.
NT = "ACGT"

print("Imports loaded; numpy", np.__version__)
'''


STEP1_MD = """## Step 1 (8 min) — Simulate three overlapping isoforms, then index their k-mers

We build a tiny "gene" with four exons (A, B, D, C). Three isoforms each
skip one of A / B / D and all share the last exon C:

- **t1 = A-B-C**   (skips D)
- **t2 = A-D-C**   (skips B)
- **t3 = B-D-C**   (skips A)

Reads from the shared exon C are compatible with all three isoforms.
Exon A distinguishes {t1, t2} from {t3}; exon B distinguishes {t1, t3}
from {t2}; exon D distinguishes {t2, t3} from {t1}. All isoforms have
the same length, which keeps the EM model simple — we will come back to
length normalisation in the EE framing at the end. We then build a
`kmer -> {isoform set}` index: this is the compatibility-class table
that the Kallisto / Salmon quantifiers use in place of position-level
alignment.
"""

STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Simulate isoforms and build the k-mer compatibility index.
# ----------------------------------------------------------------------

K = 25                  # k-mer length used by the index
EXON_LEN = 600          # bp per exon
ISOFORM_NAMES = ["t1", "t2", "t3"]
TRUE_ABUNDANCE = np.array([0.40, 0.35, 0.25])   # ground truth mixture


def random_exon(length: int, rng: np.random.Generator) -> str:
    """Return a random ACGT string of the given length."""
    # TODO: sample `length` independent nucleotides from NT with a fixed rng.
    raise NotImplementedError


def build_isoforms(seed: int = 0) -> dict[str, str]:
    """Return {isoform_name: sequence} for t1=ABC, t2=ADC, t3=BDC.

    All three isoforms share exon C. Exon A distinguishes {t1,t2} from t3;
    exon B distinguishes {t1,t3} from t2; exon D distinguishes {t2,t3} from t1.
    """
    # TODO:
    # 1. Use np.random.default_rng(seed) to draw four 600-bp exons A, B, C, D.
    # 2. Assemble t1 = A+B+C, t2 = A+D+C, t3 = B+D+C  (all 1800 bp).
    # 3. Return a dict keyed by isoform name.
    raise NotImplementedError


def build_kmer_index(isoforms: dict[str, str], k: int = K) -> dict[str, frozenset[str]]:
    """For every k-mer that appears in any isoform, list the set of isoforms
    it is found in. The set IS the compatibility class.
    """
    # TODO: slide a length-k window over each isoform; accumulate {isoform names}.
    # Return {kmer: frozenset(isoform_names)}.
    raise NotImplementedError


# isoforms = build_isoforms(seed=SEED)
# kmer_index = build_kmer_index(isoforms, k=K)
isoforms = {}
kmer_index = {}
'''


STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

K = 25
EXON_LEN = 600
ISOFORM_NAMES = ["t1", "t2", "t3"]
TRUE_ABUNDANCE = np.array([0.40, 0.35, 0.25])


def random_exon(length: int, rng: np.random.Generator) -> str:
    return "".join(rng.choice(list(NT), size=length))


def build_isoforms(seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    A = random_exon(EXON_LEN, rng)
    B = random_exon(EXON_LEN, rng)
    C = random_exon(EXON_LEN, rng)
    D = random_exon(EXON_LEN, rng)
    return {
        "t1": A + B + C,   # skip D
        "t2": A + D + C,   # skip B
        "t3": B + D + C,   # skip A
    }


def build_kmer_index(isoforms: dict, k: int = K) -> dict:
    idx: dict[str, set] = defaultdict(set)
    for name, seq in isoforms.items():
        for i in range(len(seq) - k + 1):
            idx[seq[i:i + k]].add(name)
    return {km: frozenset(s) for km, s in idx.items()}


isoforms = build_isoforms(seed=SEED)
kmer_index = build_kmer_index(isoforms, k=K)

# Report sizes and compatibility-class composition.
print("Isoform lengths:")
for n, s in isoforms.items():
    print(f"  {n}: {len(s)} bp")

class_counts = Counter(kmer_index.values())
print(f"\\nTotal distinct k-mers: {len(kmer_index):,}")
print("k-mers per compatibility class:")
for cls, n in sorted(class_counts.items(), key=lambda x: -x[1]):
    print(f"  {sorted(cls)}: {n:,}")
'''


STEP2_MD = """## Step 2 (10 min) — Sample reads and assign each to its compatibility class

Sample 2 000 reads of length 50 bp from the three isoforms, weighted by the
ground-truth abundance vector (40 / 35 / 25 %). For each read we look up
**all** of its k-mers in the index; the intersection of the per-k-mer
isoform sets is the read's compatibility class. With clean reads and a
long enough k, every read maps to exactly the isoform(s) that contain
the exon it came from.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Simulate reads and classify each into a compatibility class.
# ----------------------------------------------------------------------

N_READS = 2000
READ_LEN = 50


def sample_reads(isoforms: dict[str, str],
                 abundances: np.ndarray,
                 n_reads: int = N_READS,
                 read_len: int = READ_LEN,
                 seed: int = 1) -> list[tuple[str, str]]:
    """Return [(source_isoform, read_seq), ...] with source isoform sampled
    in proportion to `abundances` and start position uniform within the
    chosen isoform.
    """
    # TODO:
    # 1. rng = np.random.default_rng(seed)
    # 2. For each of n_reads:
    #    a. Pick an isoform name with probabilities = abundances.
    #    b. Pick a random start in [0, len(seq) - read_len].
    #    c. Append (isoform_name, seq[start:start+read_len]).
    raise NotImplementedError


def classify_read(read: str, kmer_index: dict, k: int = K) -> frozenset[str]:
    """Intersect the compatibility classes of every k-mer in `read`.

    A k-mer that is absent from the index has class = empty; the read's
    overall class is the intersection over all its k-mers.
    """
    # TODO: slide a length-k window over the read; intersect per-window classes.
    raise NotImplementedError


# reads = sample_reads(isoforms, TRUE_ABUNDANCE)
# read_classes = [classify_read(r, kmer_index) for _, r in reads]
reads = []
read_classes = []
'''

STEP2_SOLUTION = '''# Reference solution — Step 2.


def sample_reads(isoforms: dict, abundances: np.ndarray,
                 n_reads: int = N_READS, read_len: int = READ_LEN,
                 seed: int = 1) -> list:
    rng = np.random.default_rng(seed)
    names = list(isoforms.keys())
    probs = np.asarray(abundances, dtype=float)
    probs = probs / probs.sum()
    out = []
    for _ in range(n_reads):
        name = rng.choice(names, p=probs)
        seq = isoforms[name]
        start = int(rng.integers(0, len(seq) - read_len + 1))
        out.append((name, seq[start:start + read_len]))
    return out


def classify_read(read: str, kmer_index: dict, k: int = K) -> frozenset:
    if len(read) < k:
        return frozenset()
    cls = None
    for i in range(len(read) - k + 1):
        kmer = read[i:i + k]
        c = kmer_index.get(kmer, frozenset())
        if cls is None:
            cls = c
        else:
            cls = cls & c
        if not cls:
            return frozenset()
    return cls if cls is not None else frozenset()


reads = sample_reads(isoforms, TRUE_ABUNDANCE, seed=1)
read_classes = [classify_read(r, kmer_index) for _, r in reads]

# Tally compatibility-class membership.
class_to_count = Counter(read_classes)
print(f"Sampled {len(reads)} reads ({READ_LEN} bp). Compatibility-class tally:")
for cls, n in sorted(class_to_count.items(), key=lambda x: -x[1]):
    if cls:
        print(f"  {sorted(cls)}: {n:3d}")
    else:
        print(f"  (unmapped / empty class): {n:3d}")

# Sanity: fraction of unambiguous reads (class size == 1).
unique = sum(n for cls, n in class_to_count.items() if len(cls) == 1)
ambig  = sum(n for cls, n in class_to_count.items() if len(cls) > 1)
print(f"\\nUnambiguous reads: {unique}/{len(reads)}")
print(f"Ambiguous reads:   {ambig}/{len(reads)}")
'''


STEP3_MD = """## Step 3 (15 min) — E-step: fractional read assignment under current abundances

The EM algorithm treats each ambiguous read as a soft mixture over its
compatibility class. We parameterise the EM in terms of *expected reads
per isoform* η_t (because all three isoforms have the same length here,
η_t is also proportional to the molecular abundance — Kallisto / Salmon
divide by transcript length to convert; we revisit that in the EE
framing). The **posterior probability** that read `r` came from isoform
`t` is:

$$P(t \\mid r,\\eta) = \\frac{\\eta_t \\cdot \\mathbb{1}[t \\in C_r]}{\\sum_{t' \\in C_r} \\eta_{t'}}$$

— renormalised expected-read-count within the read's compatibility class.
Reads outside the class contribute zero. A uniquely-assigned read gets
all its weight on the one isoform it is compatible with.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — E-step: compute fractional responsibilities for every read.
# ----------------------------------------------------------------------

# We work with read *classes* (frozensets) rather than the reads themselves —
# all reads in the same compatibility class get identical responsibilities,
# so we can collapse to per-class counts and multiply through.

NAMES = ["t1", "t2", "t3"]
NAME_TO_IDX = {n: i for i, n in enumerate(NAMES)}


def e_step(class_counts: dict[frozenset, int], eta: np.ndarray) -> np.ndarray:
    """Return a (3,) vector of expected reads-from-each-isoform under eta.

    `eta` is the current expected-reads-per-isoform vector (the natural
    EM parameter in read units; length normalisation comes later).

    For each compatibility class C with n_C reads:
      responsibility_t = n_C * eta_t / sum_{t' in C} eta_{t'}    for t in C
      responsibility_t = 0                                       otherwise
    Sum responsibilities across classes to get expected counts per isoform.
    """
    # TODO: iterate classes, distribute n_C across the isoforms in C in proportion to eta.
    raise NotImplementedError


# class_counts = Counter(read_classes); see solution for usage
'''

STEP3_SOLUTION = '''# Reference solution — Step 3.

NAMES = ["t1", "t2", "t3"]
NAME_TO_IDX = {n: i for i, n in enumerate(NAMES)}


def e_step(class_counts: dict, eta: np.ndarray) -> np.ndarray:
    """Expected reads per isoform under current eta (read-unit parameter)."""
    exp = np.zeros(len(NAMES))
    for cls, n in class_counts.items():
        if not cls:
            continue  # drop unmapped reads
        idx = [NAME_TO_IDX[t] for t in cls if t in NAME_TO_IDX]
        if not idx:
            continue
        weights = eta[idx]
        s = weights.sum()
        if s <= 0:
            # Degenerate eta — distribute uniformly within the class.
            weights = np.ones_like(weights)
            s = weights.sum()
        share = n * weights / s
        for j, i in enumerate(idx):
            exp[i] += share[j]
    return exp


# Show one E-step from a uniform read-share prior to make the soft-assignment idea concrete.
class_counts = Counter(read_classes)
eta0 = np.full(3, N_READS / 3.0)  # naive: split reads evenly across isoforms
expected_counts = e_step(class_counts, eta0)
print("After ONE E-step from uniform read prior eta = (N/3, N/3, N/3):")
for n, e in zip(NAMES, expected_counts):
    print(f"  expected reads from {n}: {e:6.2f}")
print(f"  total reads mapped to a class: {expected_counts.sum():.0f}")
'''


STEP4_MD = """## Step 4 (12 min) — M-step + iterate to convergence

The **M-step** is trivial: the expected read counts *are* the next
estimate of η. No renormalisation needed in read units:

$$\\eta_t^{(i+1)} = E[n_t \\mid \\eta^{(i)}]$$

At the very end, convert from read counts to **molecular abundance**
(transcripts per cell) by dividing by transcript length and renormalising:

$$\\theta_t = \\frac{\\eta_t / L_t}{\\sum_{t'} \\eta_{t'} / L_{t'}}$$

Iterate E and M for ~20 rounds; track the L1 distance between consecutive
η vectors as the convergence diagnostic. EM is guaranteed to
**monotonically increase the log-likelihood**, so the iterate sequence
is stable — no oscillation, no overshoot.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Alternate E and M; track convergence; convert to molecular abundance.
# ----------------------------------------------------------------------

N_ITERS = 20


def run_em(class_counts: dict[frozenset, int], n_iters: int = N_ITERS,
           eta_init: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Run EM in read units. Return (eta_final, history) where history is
    (n_iters+1, 3) holding eta after each M-step (history[0] = eta_init).
    """
    # TODO:
    # 1. Start from eta_init or split total reads evenly across isoforms.
    # 2. For i in range(n_iters):
    #       eta = e_step(class_counts, eta)
    #       (M-step is identity in read units — expected counts ARE the new estimate)
    #       record eta into history.
    # 3. Return (eta, history).
    raise NotImplementedError


def reads_to_abundance(eta: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Convert read counts to molecular abundance (transcripts per cell).

    Divide by transcript length, then renormalise to sum to 1.
    """
    # TODO
    raise NotImplementedError


# eta_hat, history = run_em(Counter(read_classes), n_iters=N_ITERS)
# lengths = np.array([len(isoforms[n]) for n in NAMES])
# theta_hat = reads_to_abundance(eta_hat, lengths)
'''

STEP4_SOLUTION = '''# Reference solution — Step 4.

N_ITERS = 20


def run_em(class_counts: dict, n_iters: int = N_ITERS,
           eta_init: np.ndarray | None = None) -> tuple:
    total_mapped = sum(n for cls, n in class_counts.items() if cls)
    if eta_init is None:
        eta = np.full(len(NAMES), total_mapped / len(NAMES))
    else:
        eta = eta_init.copy()
    history = [eta.copy()]
    for _ in range(n_iters):
        eta = e_step(class_counts, eta)
        history.append(eta.copy())
    return eta, np.asarray(history)


def reads_to_abundance(eta: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    rate = eta / lengths
    total = rate.sum()
    if total <= 0:
        return np.ones_like(rate) / len(rate)
    return rate / total


eta_hat, history = run_em(Counter(read_classes), n_iters=N_ITERS)
lengths = np.array([len(isoforms[n]) for n in NAMES])
theta_hat = reads_to_abundance(eta_hat, lengths)

# Convert ground truth read-shares to molecular abundances for fair comparison.
# Ground-truth abundance vector means: probability a read came from this isoform.
# So expected reads per isoform is N_READS * TRUE_ABUNDANCE; convert to molecular
# abundance via the same length divide.
true_molecular = reads_to_abundance(N_READS * TRUE_ABUNDANCE, lengths)

print("EM trajectory (eta = expected reads per isoform, every other iter):")
print(pd.DataFrame(history, columns=NAMES).round(2).iloc[::2])

print("\\nRead-share recovery (compare to TRUE_ABUNDANCE):")
print(f"  Truth read shares : {np.round(TRUE_ABUNDANCE, 4)}")
print(f"  EM read shares    : {np.round(eta_hat / eta_hat.sum(), 4)}")

print("\\nMolecular abundance recovery (length-normalised, sums to 1):")
print(f"  Truth molecular   : {np.round(true_molecular, 4)}")
print(f"  EM molecular      : {np.round(theta_hat, 4)}")
print(f"  Abs error         : {np.round(np.abs(theta_hat - true_molecular), 4)}")
'''


STEP5_MD = """## Step 5 (15 min) — Plot convergence and stress-test on a different mixture

Two plots:

1. **Convergence trace** — *read share* of each isoform (η_t / Σ η) vs
   EM iteration, with the ground-truth dashed line. The trace should
   asymptote within ~5 iterations.
2. **Perturbed ground truth** — rerun the whole pipeline (sample reads,
   classify, EM) with a different ground-truth mixture (e.g. 0.7 / 0.2 / 0.1)
   to confirm the algorithm tracks the underlying mixture rather than
   memorising the first answer.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Plot convergence; rerun with a different ground-truth mixture.
# ----------------------------------------------------------------------

# TODO:
# 1. Plot history[:, 0..2] vs iteration. Overlay TRUE_ABUNDANCE as dashed lines.
# 2. Rerun the simulation with TRUE_ABUNDANCE_2 = [0.70, 0.20, 0.10]. Re-sample
#    reads, classify, run EM, and report the estimate alongside the new truth.
# 3. Print the absolute error in both cases.

TRUE_ABUNDANCE_2 = np.array([0.70, 0.20, 0.10])
'''

STEP5_SOLUTION = '''# Reference solution — Step 5.

# 1. Convergence trace — read share per iteration.
fig, ax = plt.subplots(figsize=(8, 4.5))
iters = np.arange(history.shape[0])
row_sums = history.sum(axis=1, keepdims=True)
row_sums[row_sums == 0] = 1.0
read_share = history / row_sums
colors = ["tab:blue", "tab:orange", "tab:green"]
for j, name in enumerate(NAMES):
    ax.plot(iters, read_share[:, j], color=colors[j], lw=2, label=f"{name} (EM)")
    ax.axhline(TRUE_ABUNDANCE[j], color=colors[j], linestyle="--", alpha=0.6,
               label=f"{name} truth = {TRUE_ABUNDANCE[j]:.2f}")
ax.set_xlabel("EM iteration")
ax.set_ylabel("read share estimate (eta_t / sum eta)")
ax.set_title("EM convergence on three overlapping isoforms")
ax.legend(loc="center right", ncol=2, fontsize=9)
ax.set_ylim(0, 1)
plt.tight_layout()
plt.show()

# 2. Stress-test on a perturbed mixture.
TRUE_ABUNDANCE_2 = np.array([0.70, 0.20, 0.10])

reads_2 = sample_reads(isoforms, TRUE_ABUNDANCE_2, seed=7)
classes_2 = [classify_read(r, kmer_index) for _, r in reads_2]
eta_hat_2, history_2 = run_em(Counter(classes_2), n_iters=N_ITERS)
read_share_2 = eta_hat_2 / eta_hat_2.sum()
theta_hat_2 = reads_to_abundance(eta_hat_2, lengths)
true_molecular_2 = reads_to_abundance(N_READS * TRUE_ABUNDANCE_2, lengths)

print("Mixture 1 (truth 0.40 / 0.35 / 0.25):")
print(f"  EM read share: {np.round(eta_hat / eta_hat.sum(), 3)}   abs err: "
      f"{np.round(np.abs(eta_hat / eta_hat.sum() - TRUE_ABUNDANCE), 3)}")
print(f"  EM molecular : {np.round(theta_hat, 3)}   truth molecular: {np.round(true_molecular, 3)}")
print()
print("Mixture 2 (truth 0.70 / 0.20 / 0.10):")
print(f"  EM read share: {np.round(read_share_2, 3)}   abs err: "
      f"{np.round(np.abs(read_share_2 - TRUE_ABUNDANCE_2), 3)}")
print(f"  EM molecular : {np.round(theta_hat_2, 3)}   truth molecular: {np.round(true_molecular_2, 3)}")

# Bar plot comparing both mixtures (read-share view).
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
x = np.arange(3)
width = 0.4
for ax, truth, est, title in [
    (axes[0], TRUE_ABUNDANCE,   eta_hat / eta_hat.sum(), "Mixture 1: 40/35/25"),
    (axes[1], TRUE_ABUNDANCE_2, read_share_2,            "Mixture 2: 70/20/10"),
]:
    ax.bar(x - width / 2, truth, width=width, label="truth", color="#888")
    ax.bar(x + width / 2, est,   width=width, label="EM",    color="tab:blue")
    ax.set_xticks(x)
    ax.set_xticklabels(NAMES)
    ax.set_title(title)
    ax.set_ylim(0, 1)
    ax.legend()
axes[0].set_ylabel("read share")
plt.tight_layout()
plt.show()
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces of the pipeline.
If you ran the reference solutions above they should all pass; if you
wrote your own and an assert fails, revisit the corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Isoforms have the expected lengths (all three are 3-exon, each skips a different exon).
assert set(isoforms.keys()) == {"t1", "t2", "t3"}, "isoform names wrong"
for name in ("t1", "t2", "t3"):
    assert len(isoforms[name]) == 3 * EXON_LEN, f"{name} should be 3 exons long"

# 2. The k-mer index returns reasonable compatibility classes.
#    The shared exon C makes many k-mers map to all three isoforms.
all_three = sum(1 for cls in kmer_index.values() if cls == frozenset({"t1", "t2", "t3"}))
assert all_three > 100, f"expected many k-mers in shared exon C; got {all_three}"

# 3. The E-step preserves the total reads-mapped-to-a-class count.
cc = Counter(read_classes)
mapped = sum(n for cls, n in cc.items() if cls)
exp = e_step(cc, np.full(3, mapped / 3.0))
assert math.isclose(exp.sum(), mapped, rel_tol=1e-9), \
    f"E-step total {exp.sum()} != mapped reads {mapped}"

# 4. reads_to_abundance output is a valid probability simplex.
ra = reads_to_abundance(np.array([700.0, 300.0, 500.0]), np.array([1800, 1800, 1800]))
assert math.isclose(ra.sum(), 1.0, abs_tol=1e-9)
assert np.all(ra >= 0)

# 5. EM recovers the true READ SHARE to within 0.05 absolute error per isoform.
em_read_share = eta_hat / eta_hat.sum()
abs_err = np.abs(em_read_share - TRUE_ABUNDANCE)
print(f"EM read share : {np.round(em_read_share, 3)}")
print(f"Ground truth  : {np.round(TRUE_ABUNDANCE, 3)}")
print(f"Abs error/iso : {np.round(abs_err, 3)}")
assert abs_err.max() < 0.05, f"EM did not recover read shares: max abs err {abs_err.max():.3f}"

# 6. Length-normalised molecular abundance also matches truth within 0.05.
abs_err_mol = np.abs(theta_hat - true_molecular)
print(f"EM molecular  : {np.round(theta_hat, 3)}")
print(f"Truth molec.  : {np.round(true_molecular, 3)}")
assert abs_err_mol.max() < 0.05, f"molecular abundance off: max abs err {abs_err_mol.max():.3f}"

# 7. Convergence: by iteration 15, the iterate has stabilised
#    (L1 change in read-share < 0.005).
ten_share = history[14] / max(1e-9, history[14].sum())
eleven_share = history[15] / max(1e-9, history[15].sum())
late_diff = np.abs(eleven_share - ten_share).sum()
assert late_diff < 0.005, f"EM still moving at iter 15: L1 step = {late_diff:.4f}"

print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — EM as iterative soft-decision decoding

You implemented the canonical **iterative soft-decoder**:

1. **Latent variable**: which isoform did each read come from? Unknown.
2. **Soft assignment (E-step)**: given the current abundance estimate
   `θ`, compute the **posterior probability** that read *r* came from
   each isoform in its compatibility class. This is the same arithmetic
   as the *extrinsic information* exchanged between component decoders in
   a turbo code, or the *messages* passed in belief propagation on a
   factor graph.
3. **Parameter update (M-step)**: re-estimate `θ` as the renormalised
   sum of soft assignments. Equivalent to maximum-likelihood given the
   currently-believed responsibilities.
4. **Iterate**: EM monotonically increases the log-likelihood, so each
   pass is a guaranteed-non-regressive step toward a local optimum.
   Within 5–10 iterations the iterate has converged.

The same recipe appears across signal processing:

- **Turbo codes / LDPC decoders** — soft information about each bit
  exchanged between component decoders until they agree.
- **Kalman smoothers** — forward / backward soft estimates of hidden
  state combined into a posterior.
- **GMM clustering** — soft-assign points to Gaussians (E), re-estimate
  means and covariances (M), iterate.

So when Kallisto reports a TPM for an isoform with 80 % of its sequence
shared with a sibling isoform, it is **not pretending the reads were
unambiguous**: it is using exactly the kind of iterative soft-decoder
that lets a deep-space probe extract a 10-bit/s signal from -30 dB SNR.
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
