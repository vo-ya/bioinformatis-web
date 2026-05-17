"""Build exercise.ipynb for L04 — Variant Calling: From Aligned Reads to Called Differences.

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


TITLE_MD = """# L04 — Variant Calling: From Aligned Reads to Called Differences

In this exercise you simulate a small 2 kb diploid region with a single
heterozygous SNV, generate noisy reads with Phred-coded errors, build a
pileup at the candidate site, compute genotype log-likelihoods, normalise
with a uniform prior to get a posterior, call the MAP genotype, and then
sweep the per-base error rate to see how robust the caller is.
"""


AHA_MD = """> **Aha.** Variant calling is not thresholding ("≥ 20% reads carry the
> alt allele, call het"). It is **Bayesian MAP inference**: the pileup
> column is a sufficient statistic, each read contributes a log-likelihood
> term computed from its Phred quality, and the call is the
> argmax-of-posterior over `{AA, AG, GG}`. Once you have the math, you
> get a calibrated genotype quality for free.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 matplotlib==3.8.4 -q
"""


IMPORTS = """import math
from dataclasses import dataclass
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Deterministic for the whole notebook.
RNG_SEED = 42
np.random.seed(RNG_SEED)

# Canonical DNA alphabet.
BASES = "ACGT"
BASES_SET = set(BASES)


def phred_to_perr(q: int) -> float:
    \"\"\"Phred quality -> per-base error probability.\"\"\"
    return 10.0 ** (-q / 10.0)


def perr_to_phred(p: float) -> int:
    \"\"\"Per-base error probability -> Phred quality (rounded).\"\"\"
    return int(round(-10.0 * math.log10(max(p, 1e-12))))


print(f"Q20 -> p_err = {phred_to_perr(20):.4f}")
print(f"Q30 -> p_err = {phred_to_perr(30):.4f}")
"""


STEP1_MD = """## Step 1 (8 min) — Simulate a diploid region and reads

We build the toy genome ourselves: a 2 kb random reference plus a diploid
"true" genotype that is identical to the reference everywhere except at
position 1000, where one chromosome carries `A` (reference) and the other
carries `G` (alternate). That makes the ground-truth genotype at the
variant site **0/1 (het A/G)**.

Reads are 100 bp, drawn uniformly from random positions, randomly assigned
to one of the two chromosomes, and corrupted by Phred-coded per-base
errors. With ~50× depth and 500 reads we get plenty of coverage at the
variant site.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Simulate the diploid genome and the reads.
# ----------------------------------------------------------------------

REF_LEN     = 2000          # 2 kb reference
VARIANT_POS = 1000          # 0-based variant position
REF_ALLELE  = "A"
ALT_ALLELE  = "G"
N_READS     = 500
READ_LEN    = 100
MEAN_Q      = 30            # per-base Phred quality; p_err ~= 10^-3


def simulate_reference(length: int, seed: int = 0) -> str:
    """Random ACGT reference of given length."""
    # TODO: draw `length` letters from BASES with numpy and return as a string.
    raise NotImplementedError


def build_diploid(ref: str, var_pos: int, ref_allele: str, alt_allele: str) -> tuple[str, str]:
    """Return (chrom_A, chrom_B): one is reference, one carries the ALT at var_pos."""
    # TODO: chrom_A = ref (carries ref allele at var_pos)
    #       chrom_B = ref with position var_pos overwritten by alt_allele
    raise NotImplementedError


def simulate_reads(chroms: tuple[str, str], n_reads: int, read_len: int,
                   mean_q: int = MEAN_Q, seed: int = 0):
    """Return list of (read_seq, quals, ref_start) tuples.

    Each read is drawn from a random chromosome and starting position;
    every base has its own quality drawn near mean_q; the observed base is
    the chromosome's base, flipped to a different base with probability
    10^(-q/10).
    """
    # TODO: implement read simulation as described.
    raise NotImplementedError


# ref:       str of length 2000
# chroms:    (chrom_A, chrom_B)
# reads:     list[(seq, quals, ref_start)]
ref = None
chroms = None
reads = None
'''


STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

REF_LEN     = 2000
VARIANT_POS = 1000
REF_ALLELE  = "A"
ALT_ALLELE  = "G"
N_READS     = 500
READ_LEN    = 100
MEAN_Q      = 30


def simulate_reference(length: int, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    arr = rng.choice(list(BASES), size=length)
    return "".join(arr)


def build_diploid(ref: str, var_pos: int, ref_allele: str, alt_allele: str):
    assert ref[var_pos] == ref_allele, (
        f"reference base at {var_pos} is {ref[var_pos]!r}, not {ref_allele!r}"
    )
    chrom_A = ref
    chrom_B = ref[:var_pos] + alt_allele + ref[var_pos + 1:]
    return chrom_A, chrom_B


def _flip_base(true_base: str, rng) -> str:
    """Return a random base != true_base."""
    others = [b for b in BASES if b != true_base]
    return rng.choice(others)


def simulate_reads(chroms, n_reads: int, read_len: int,
                   mean_q: int = MEAN_Q, seed: int = 0,
                   error_scale: float = 1.0):
    """error_scale multiplies every per-base error probability (>=1 worsens, <=1 improves)."""
    rng = np.random.default_rng(seed)
    L = len(chroms[0])
    reads = []
    for _ in range(n_reads):
        chrom = chroms[rng.integers(0, 2)]
        start = int(rng.integers(0, L - read_len + 1))
        true_window = chrom[start:start + read_len]
        # Quality scores: small jitter around mean_q so the likelihood
        # downstream actually depends on the value.
        quals = np.clip(rng.normal(mean_q, 3, size=read_len), 5, 41).astype(int)
        seq_chars = []
        for true_base, q in zip(true_window, quals):
            p_err = phred_to_perr(int(q)) * error_scale
            p_err = min(p_err, 0.75)
            if rng.random() < p_err:
                seq_chars.append(_flip_base(true_base, rng))
            else:
                seq_chars.append(true_base)
        reads.append(("".join(seq_chars), quals.tolist(), start))
    return reads


ref = simulate_reference(REF_LEN, seed=RNG_SEED)
# Force the reference base at VARIANT_POS to REF_ALLELE so the diploid is well-defined.
ref = ref[:VARIANT_POS] + REF_ALLELE + ref[VARIANT_POS + 1:]
chroms = build_diploid(ref, VARIANT_POS, REF_ALLELE, ALT_ALLELE)
reads = simulate_reads(chroms, N_READS, READ_LEN, MEAN_Q, seed=RNG_SEED)

print(f"Reference: {REF_LEN} bp, base at variant pos {VARIANT_POS} = {ref[VARIANT_POS]!r}")
print(f"Chromosome A at {VARIANT_POS}: {chroms[0][VARIANT_POS]!r}")
print(f"Chromosome B at {VARIANT_POS}: {chroms[1][VARIANT_POS]!r}")
print(f"Simulated {len(reads)} reads of {READ_LEN} bp each (~{N_READS * READ_LEN / REF_LEN:.0f}x depth)")
'''


STEP2_MD = """## Step 2 (10 min) — Build the pileup at the variant position

The pileup is the central data structure for per-position variant calling:
for one reference column, list every read that overlaps that column, the
base each read carries there, and that base's quality. It is a sufficient
statistic — every per-position caller, from `bcftools` to GATK, eats a
pileup as input.

Build the pileup column at `VARIANT_POS`, count base composition, and
estimate the variant allele fraction (VAF). The true het VAF is 0.5; with
500 reads and ~25 covering this site you should land near that, plus or
minus sampling noise.
"""


STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Pileup at the variant position.
# ----------------------------------------------------------------------


@dataclass
class PileupBase:
    base: str   # observed base in the read
    qual: int   # Phred quality of that base


def build_pileup(reads, position: int, read_len: int = READ_LEN) -> list[PileupBase]:
    """Return all (base, qual) observations covering `position`."""
    # TODO: for each (seq, quals, start) in reads, if start <= position < start + read_len,
    # append PileupBase(seq[offset], quals[offset]) with offset = position - start.
    raise NotImplementedError


def pileup_summary(pileup):
    """Return dict {base: count} and the VAF for the ALT allele."""
    # TODO: count letters; VAF = count[ALT_ALLELE] / total.
    raise NotImplementedError


# pileup_at_variant = build_pileup(reads, VARIANT_POS)
# counts, vaf = pileup_summary(pileup_at_variant)
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.


@dataclass
class PileupBase:
    base: str
    qual: int


def build_pileup(reads, position: int, read_len: int = READ_LEN):
    out = []
    for seq, quals, start in reads:
        if start <= position < start + read_len:
            offset = position - start
            out.append(PileupBase(seq[offset], int(quals[offset])))
    return out


def pileup_summary(pileup):
    counts = Counter(p.base for p in pileup)
    total = sum(counts.values())
    vaf = counts.get(ALT_ALLELE, 0) / total if total else 0.0
    return dict(counts), vaf


pileup_at_variant = build_pileup(reads, VARIANT_POS)
counts, vaf = pileup_summary(pileup_at_variant)

print(f"Pileup at position {VARIANT_POS}: depth = {len(pileup_at_variant)}")
print(f"  base composition: {counts}")
print(f"  VAF({ALT_ALLELE}) = {vaf:.3f}   (truth for a het: 0.500)")
print()
print("First 12 (base, qual) observations:")
for pb in pileup_at_variant[:12]:
    print(f"  {pb.base}  Q{pb.qual}  (p_err={phred_to_perr(pb.qual):.4f})")
'''


STEP3_MD = """## Step 3 (12 min) — Genotype log-likelihoods

For a diploid site with two candidate alleles `A` and `G`, the candidate
genotypes are `AA`, `AG`, `GG`. The per-read likelihood under genotype
`G = (a1, a2)` with per-base error probability `e` is:

- if the read's observed base equals `a1` or `a2`: `P(b | G) = 0.5·(1 − e) + 0.5·(e/3)`
- if it equals neither: `P(b | G) = e/3`

(The 0.5 weights each of the two parental alleles equally because we have
no idea which chromosome the read came from.) For a **homozygous** site
the two parental alleles are the same letter so the per-read likelihood
collapses to `(1 − e)` for a match and `e/3` for a mismatch.

We work in **log-space** because likelihoods underflow once you multiply
~25 reads together:

`log L(G) = Σ_reads log P(b_r | G, e_r)`

Implement this for `AA, AG, GG` and return the three log-likelihoods.
"""


STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Genotype log-likelihoods.
# ----------------------------------------------------------------------

GENOTYPES = ["AA", "AG", "GG"]   # ordered: hom-ref, het, hom-alt


def read_log_likelihood(observed_base: str, qual: int, genotype: str) -> float:
    """log P(observed_base | genotype, per-base-error from qual).

    genotype is a 2-character string like "AG"; alleles are unordered.
    """
    # TODO:
    # e = phred_to_perr(qual)
    # if observed_base equals either allele in genotype:
    #     allele match probability is 0.5*(1-e) + 0.5*(e/3) when alleles differ,
    #     or simply (1-e) when both alleles are the same letter.
    # else: probability is e/3.
    # Return math.log(probability).
    raise NotImplementedError


def genotype_log_likelihoods(pileup, genotypes=GENOTYPES) -> dict[str, float]:
    """Return {genotype: log L(genotype | pileup)}."""
    # TODO: sum read_log_likelihood over every (base, qual) in the pileup.
    raise NotImplementedError


# logL = genotype_log_likelihoods(pileup_at_variant)
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.

GENOTYPES = ["AA", "AG", "GG"]


def read_log_likelihood(observed_base: str, qual: int, genotype: str) -> float:
    e = phred_to_perr(qual)
    a1, a2 = genotype[0], genotype[1]
    if a1 == a2:
        # Homozygous: match -> 1-e; mismatch -> e/3.
        p = (1.0 - e) if observed_base == a1 else (e / 3.0)
    else:
        # Heterozygous: each parental allele weighted 0.5.
        if observed_base == a1 or observed_base == a2:
            p = 0.5 * (1.0 - e) + 0.5 * (e / 3.0)
        else:
            p = e / 3.0
    return math.log(max(p, 1e-300))


def genotype_log_likelihoods(pileup, genotypes=GENOTYPES):
    return {
        g: sum(read_log_likelihood(pb.base, pb.qual, g) for pb in pileup)
        for g in genotypes
    }


logL = genotype_log_likelihoods(pileup_at_variant)
print("Genotype log-likelihoods at the variant position:")
for g, ll in logL.items():
    print(f"  {g}: log L = {ll:10.3f}")

# Sanity-plot: bar chart of log-likelihoods.
fig, ax = plt.subplots(figsize=(5, 3))
ax.bar(list(logL.keys()), list(logL.values()), color=["#6c8ebf", "#82b366", "#b85450"])
ax.set_ylabel("log L(G | data)")
ax.set_title(f"Genotype log-likelihoods at site {VARIANT_POS}")
ax.axhline(0, color="black", lw=0.5)
for i, (g, ll) in enumerate(logL.items()):
    ax.text(i, ll, f"{ll:.1f}", ha="center", va="bottom" if ll < 0 else "top")
plt.tight_layout()
plt.show()
'''


STEP4_MD = """## Step 4 (15 min) — Bayes-rule posterior, MAP call, genotype quality

With a uniform prior `P(G) = 1/3` for each of `AA, AG, GG`, the posterior
is simply the normalised likelihood. To avoid numerical underflow we work
in log-space and **log-sum-exp** to normalise:

`log P(G | D) = log L(G) − logsumexp_{G'} log L(G')`

The **MAP genotype** is `argmax_G P(G | D)`. The **genotype quality** is
the Phred-scaled probability of error of the MAP call:

`GQ = −10 log₁₀ (1 − P(MAP | D))`

A GQ of 30 means the caller is 99.9% sure of the MAP call. Implement the
posterior, identify the MAP genotype, and compute its GQ.
"""


STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Posterior, MAP genotype, genotype quality.
# ----------------------------------------------------------------------


def logsumexp(values):
    """Numerically stable log(sum(exp(v) for v in values))."""
    # TODO: subtract max first, then exp+sum+log+add back.
    raise NotImplementedError


def posterior(log_likelihoods: dict[str, float],
              prior: dict[str, float] | None = None) -> dict[str, float]:
    """Return {genotype: posterior probability}. prior defaults to uniform."""
    # TODO: log_post = log_L + log_prior - logsumexp(log_L + log_prior).
    # Return exp(log_post).
    raise NotImplementedError


def map_genotype(post: dict[str, float]) -> tuple[str, float]:
    """Return (best_genotype, posterior probability)."""
    # TODO: argmax over post.
    raise NotImplementedError


def genotype_quality(post: dict[str, float]) -> float:
    """Phred-scaled probability that the MAP genotype is wrong."""
    # TODO: 1 - P(MAP) is the error probability; convert to Phred.
    raise NotImplementedError


# post = posterior(logL)
# gt, p_gt = map_genotype(post)
# gq = genotype_quality(post)
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.


def logsumexp(values):
    vals = list(values)
    m = max(vals)
    return m + math.log(sum(math.exp(v - m) for v in vals))


def posterior(log_likelihoods, prior=None):
    genotypes = list(log_likelihoods.keys())
    if prior is None:
        prior = {g: 1.0 / len(genotypes) for g in genotypes}
    log_unnorm = {g: log_likelihoods[g] + math.log(prior[g]) for g in genotypes}
    Z = logsumexp(log_unnorm.values())
    return {g: math.exp(log_unnorm[g] - Z) for g in genotypes}


def map_genotype(post):
    g_best = max(post, key=post.get)
    return g_best, post[g_best]


def genotype_quality(post):
    g_best, p_best = map_genotype(post)
    p_err = max(1.0 - p_best, 1e-12)
    return -10.0 * math.log10(p_err)


post = posterior(logL)
gt, p_gt = map_genotype(post)
gq = genotype_quality(post)

print("Posterior P(G | data) with uniform prior:")
for g, p in post.items():
    print(f"  {g}: P = {p:.6f}")
print()
print(f"MAP genotype: {gt}   posterior = {p_gt:.6f}")
print(f"Genotype quality GQ = {gq:.1f}  (Phred-scaled)")
print()
print(f"Truth: het {REF_ALLELE}/{ALT_ALLELE}   Called: {gt}   "
      f"{'CORRECT' if gt == REF_ALLELE + ALT_ALLELE else 'WRONG'}")

# Posterior bar chart.
fig, ax = plt.subplots(figsize=(5, 3))
ax.bar(list(post.keys()), list(post.values()), color=["#6c8ebf", "#82b366", "#b85450"])
ax.set_ylabel("P(G | data)")
ax.set_title(f"Posterior over genotypes at site {VARIANT_POS} (uniform prior)")
ax.set_ylim(0, 1.05)
for i, (g, p) in enumerate(post.items()):
    ax.text(i, p, f"{p:.3f}", ha="center", va="bottom")
plt.tight_layout()
plt.show()
'''


STEP5_MD = """## Step 5 (15 min) — Sweep the per-base error rate

The whole point of the Bayesian framing is that it degrades gracefully:
as base-quality drops, the likelihood functions get flatter, the
posterior gets less peaky, and eventually the caller starts to make
mistakes. Quantify that — and notice that the failure mode at a true
hom-ref site is different from at a true het site:

- **Het site (truth = AG).** Both alleles are abundant; even at high error
  the het signal is hard to wash out (errors mostly push toward `C`/`T`,
  not toward the other het allele).
- **Hom-ref site (truth = AA).** Only one allele is present. Every error
  is a fake "ALT" observation. As the error rate climbs, the caller
  starts to spuriously prefer `AG` — a **false-positive variant call**.

Sweep four error scales, repeat each setting with many seeds, and report
the genotype accuracy at both a het site (pos 1000) and a hom-ref site
(pos 500, where we force the reference to be `A`).
"""


STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Sweep error rate; measure caller accuracy at het and hom-ref.
# ----------------------------------------------------------------------

# Multipliers on the Q30 base-error rate (Q30 -> 0.001 per-base error).
# At 100x we reach ~10% effective per-base error, which breaks hom-ref calls.
ERROR_SCALES = [0.0, 1.0, 10.0, 100.0]
N_REPLICATES = 30                # independent seeds per setting

# Pick a hom-ref site: any position not at VARIANT_POS where the reference is 'A'.
HOMREF_POS = 500   # we will force ref[HOMREF_POS] = 'A' below
TRUTH_HET    = REF_ALLELE + ALT_ALLELE   # "AG" at VARIANT_POS
TRUTH_HOMREF = REF_ALLELE + REF_ALLELE   # "AA" at HOMREF_POS


def call_site(reads, position: int) -> str:
    """Given a reads list and a position, run the full caller and return MAP genotype."""
    # TODO:
    # 1. pu = build_pileup(reads, position)
    # 2. compute log-likelihoods over GENOTYPES on pu
    # 3. compute posterior
    # 4. return map_genotype(post)[0]   (handle empty pileup -> "??")
    raise NotImplementedError


# Loop ERROR_SCALES x N_REPLICATES; collect het_accuracy and homref_accuracy.
'''


STEP5_SOLUTION = '''# Reference solution — Step 5.

ERROR_SCALES = [0.0, 1.0, 10.0, 100.0]
N_REPLICATES = 30
HOMREF_POS = 500
TRUTH_HET    = REF_ALLELE + ALT_ALLELE   # "AG"
TRUTH_HOMREF = REF_ALLELE + REF_ALLELE   # "AA"

# Force ref[HOMREF_POS] = 'A' so the hom-ref truth is well-defined; rebuild chroms.
ref_local = ref[:HOMREF_POS] + REF_ALLELE + ref[HOMREF_POS + 1:]
chroms_local = build_diploid(ref_local, VARIANT_POS, REF_ALLELE, ALT_ALLELE)


def call_site(reads_list, position: int) -> str:
    pu = build_pileup(reads_list, position)
    if not pu:
        return "??"
    ll = genotype_log_likelihoods(pu)
    post_local = posterior(ll)
    gt_local, _ = map_genotype(post_local)
    return gt_local


het_acc = {}
homref_acc = {}
for es in ERROR_SCALES:
    het_calls    = []
    homref_calls = []
    for i in range(N_REPLICATES):
        sim = simulate_reads(chroms_local, N_READS, READ_LEN, MEAN_Q,
                             seed=1000 + i, error_scale=es)
        het_calls.append(call_site(sim, VARIANT_POS))
        homref_calls.append(call_site(sim, HOMREF_POS))
    het_acc[es]    = sum(1 for c in het_calls    if c == TRUTH_HET)    / N_REPLICATES
    homref_acc[es] = sum(1 for c in homref_calls if c == TRUTH_HOMREF) / N_REPLICATES
    eff = phred_to_perr(MEAN_Q) * es if es > 0 else 0.0
    print(f"  error_scale={es:>5.1f}  (eff p_err~={eff:.4f})  "
          f"het_acc={het_acc[es]:.2f}  homref_acc={homref_acc[es]:.2f}  "
          f"hom-ref calls: {Counter(homref_calls)}")

accuracy = het_acc          # kept as an alias for self-check below

# Plot accuracy vs effective error rate.
xs = [max(phred_to_perr(MEAN_Q) * es, 1e-5) for es in ERROR_SCALES]
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(xs, [het_acc[es]    for es in ERROR_SCALES], "o-",
        color="#82b366", lw=2, markersize=10, label="het site (truth AG)")
ax.plot(xs, [homref_acc[es] for es in ERROR_SCALES], "s-",
        color="#b85450", lw=2, markersize=10, label="hom-ref site (truth AA)")
for x, es in zip(xs, ERROR_SCALES):
    ax.annotate(f"x{es}", (x, homref_acc[es]),
                textcoords="offset points", xytext=(8, -16), fontsize=9)
ax.set_xlabel("effective per-base error probability")
ax.set_ylabel(f"genotype accuracy ({N_REPLICATES} replicates)")
ax.set_ylim(-0.05, 1.05)
ax.set_xscale("log")
ax.set_title("Caller robustness vs sequencing error rate")
ax.legend(loc="lower left")
ax.grid(True, alpha=0.3, which="both")
plt.tight_layout()
plt.show()
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces. If you ran the
reference solutions above they should all pass; if any fails, revisit the
corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Phred conversion matches the canonical Q20 -> 1% error.
assert abs(phred_to_perr(20) - 0.01) < 1e-9, "Q20 should be 1% error"

# 2. Pileup at variant site has reasonable depth and a non-trivial alt allele count.
assert len(pileup_at_variant) >= 15, (
    f"variant-site depth suspiciously low: {len(pileup_at_variant)}"
)
counts_check, vaf_check = pileup_summary(pileup_at_variant)
assert counts_check.get(REF_ALLELE, 0) >= 3, "ref allele count too low"
assert counts_check.get(ALT_ALLELE, 0) >= 3, "alt allele count too low"

# 3. VAF for a true het should sit near 0.5 (allow wide band for sampling noise).
assert 0.25 <= vaf_check <= 0.75, (
    f"het VAF should be near 0.5; got {vaf_check:.3f}"
)

# 4. MAP genotype at the variant site is the true het AG.
assert gt == "AG", f"MAP genotype expected 'AG', got {gt!r}"

# 5. Heterozygous log-likelihood beats both homozygous alternatives.
assert logL["AG"] > logL["AA"], "het should beat hom-ref"
assert logL["AG"] > logL["GG"], "het should beat hom-alt"

# 6. Posterior probabilities sum to 1.
assert abs(sum(post.values()) - 1.0) < 1e-9, "posterior must normalise"

# 7. Genotype quality is high at default Q30 (caller should be very confident).
assert gq >= 20, f"GQ at default Q30 expected >=20; got {gq:.1f}"

# 8. Het accuracy at zero error is 1.0 (clean signal, easy call).
assert het_acc[0.0] == 1.0, f"het accuracy at error_scale=0 must be 1; got {het_acc[0.0]:.2f}"

# 9. Hom-ref accuracy degrades monotonically (with some noise) as error rises.
#    At 100x the caller should make false-positive variant calls at the hom-ref site.
assert homref_acc[0.0] == 1.0, (
    f"hom-ref accuracy at error_scale=0 must be 1; got {homref_acc[0.0]:.2f}"
)
assert homref_acc[100.0] <= homref_acc[1.0], (
    "hom-ref accuracy at 100x error should not exceed accuracy at 1x error"
)
assert homref_acc[100.0] < 1.0, (
    f"hom-ref accuracy at 100x error should drop below 1 (caller emits false positives); "
    f"got {homref_acc[100.0]:.2f}"
)

# 10. logsumexp on a singleton equals the value itself.
assert abs(logsumexp([0.7]) - 0.7) < 1e-12, "logsumexp([x]) must equal x"

print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — signal detection with unequal priors

What you implemented is a **maximum a posteriori (MAP) detector** with a
finite hypothesis set `{AA, AG, GG}`:

1. **Sufficient statistic.** The pileup column is a data-compression
   step. Whatever else the read carries about other positions, only the
   `(base, qual)` pair at this column matters for this column's call.
   That is the bioinformatics analogue of a matched-filter projection
   onto the signal subspace before the decision stage.
2. **Log-likelihood = inner product.** `log L(G) = Σ_reads log P(b_r | G, e_r)`
   is a sum over independent observations — i.e. an inner product
   between the data and the log-emission template of each genotype. The
   decision rule is "pick the genotype whose template best aligns with
   the data."
3. **Unequal priors.** Real callers don't use a uniform prior. Population
   allele frequencies (gnomAD, 1000 Genomes) give `P(AG)` ≈ 2pq and
   `P(GG)` ≈ q². That is just the **a priori** information a Bayes
   detector folds in. With a strong prior toward `AA` the same data needs
   higher evidence to call a het; that is exactly what a population
   prior buys you for rare variants.
4. **Genotype quality = Phred(P_err).** The genotype quality is Phred-scaled
   `1 − P(MAP | D)`, the false-alarm probability of the MAP decision.
   Same math as a CFAR detector reporting `P_FA` in dB.
5. **Graceful degradation.** As `e` rises, every per-read likelihood gets
   pushed toward `1/4` for all genotypes. The likelihood ratio collapses,
   the posterior flattens, and the GQ drops. Calling does not break
   abruptly — it loses confidence the way any well-posed detector should.

Thresholding by VAF ("≥ 20% reads carry the alt allele, call het") has
none of these properties: it has no notion of per-base confidence, no
prior, no calibrated quality. That is why every modern caller is
Bayesian.
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
