"""Build the L01 Colab notebook: exercise.ipynb.

Re-run with `python3 build_notebook.py` after editing this script.
The output is a single notebook with hidden solution cells (Colab collapses
them via `metadata.jupyter.source_hidden = true`).
"""

from __future__ import annotations

import os

# Make the shared Colab-form helper importable from the parent dir.
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from apply_colab_form import apply_colab_form  # noqa: E402

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook


def md(source: str) -> nbformat.NotebookNode:
    return new_markdown_cell(source)


def code(source: str, hidden: bool = False) -> nbformat.NotebookNode:
    cell = new_code_cell(source)
    if hidden:
        cell.metadata["jupyter"] = {"source_hidden": True}
        # Colab also respects `cellView: form` for collapse-style hides,
        # but `source_hidden` is the portable Jupyter standard.
    return cell


def step_block(
    n: int,
    minutes: int,
    title: str,
    why: str,
    todo_code: str,
    solution_code: str,
) -> list[nbformat.NotebookNode]:
    """Emit the 4-cell block for one of the 5 steps."""
    cells: list[nbformat.NotebookNode] = []
    cells.append(md(f"## Step {n} ({minutes} min) — {title}\n\n{why}"))
    cells.append(code(todo_code))
    cells.append(md("*Click ▶ to expand the reference solution.*"))
    cells.append(code(solution_code, hidden=True))
    return cells


# ---------------------------------------------------------------------------
# Cell content
# ---------------------------------------------------------------------------

TITLE_MD = """# L01 — Foundations: From Cells to Sequences to FASTQ

In this exercise you will generate a synthetic Illumina-style FASTQ from
scratch, parse it, and fingerprint the sequencer's error profile. The aim
is to make the link between *symbol on disk* (Phred character) and
*physics of measurement* (per-base error probability) concrete in code."""

AHA_MD = """> **Aha.** A FASTQ file is a lossy compression of raw trace
> intensity into two parallel strings: a base call and a per-base
> log-confidence (`Q = −10 log₁₀ P_err`). The position-wise quality
> decay you'll plot in Step 4 is a *fingerprint of the instrument* —
> the same fingerprint a QC tool like FastQC looks for."""

PREAMBLE_CODE = """# Pinned install — L01 only needs numpy + matplotlib (Biopython kept
# light because the FASTQ work is pure-Python). Quiet flag keeps Colab tidy.
!pip install numpy==1.26.4 matplotlib==3.8.4 -q"""

IMPORTS_CODE = """import io
import math
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt

# Reproducibility — fix every random draw downstream.
SEED = 42
rng = np.random.default_rng(SEED)
np.random.seed(SEED)

# Constants used across steps.
N_READS = 10_000
READ_LEN = 150
PHRED_OFFSET = 33  # Sanger / Illumina 1.8+ encoding."""

# ----- Step 1 ---------------------------------------------------------------

STEP1_WHY = (
    "Real Illumina reads degrade as the cycle count grows — dephasing, "
    "fluorophore wear, and phasing/prephasing errors cause Phred to drift "
    "down from ~Q38 at cycle 1 to ~Q20 at cycle 150. We *simulate* that "
    "physics in three lines of NumPy so the rest of the notebook has a "
    "deterministic FASTQ to chew on. Goal: write `fastq_text` (a single "
    "string) containing 10 000 four-line records."
)

STEP1_TODO = '''# TODO: build a synthetic FASTQ as one big string `fastq_text`.
#
# Specs:
#   - N_READS reads of length READ_LEN.
#   - Per-cycle quality decays roughly linearly from Q≈38 at cycle 0
#     to Q≈20 at cycle READ_LEN-1, with Gaussian jitter of ~2 Phred units.
#   - Convert Q -> error probability P = 10**(-Q/10).
#   - At each base, with probability P emit a *wrong* nucleotide
#     (uniform over the other three); otherwise emit a uniform random
#     correct base. We have no "true" reference here — the point is the
#     error model, not the genome.
#   - Encode each Q as chr(Q + PHRED_OFFSET).
#   - Record format (4 lines):
#         @read_{i}
#         <seq>
#         +
#         <qualstring>
#
# Hint: use the module-level `rng` (numpy default_rng) for all draws.

def build_synthetic_fastq(n_reads: int = N_READS,
                          read_len: int = READ_LEN) -> str:
    # TODO replace this stub
    return ""

fastq_text = build_synthetic_fastq()
print(f"fastq_text length: {len(fastq_text):,} chars")
print("First record preview:\\n" + "\\n".join(fastq_text.splitlines()[:4]))
'''

STEP1_SOL = '''def build_synthetic_fastq(n_reads: int = N_READS,
                          read_len: int = READ_LEN) -> str:
    bases = np.array(list("ACGT"))
    # Quality curve: linear decay 38 -> 20 with 2 Phred jitter.
    base_curve = np.linspace(38.0, 20.0, read_len)
    # Shape (n_reads, read_len): one jittered Phred per base.
    q = base_curve + rng.normal(0.0, 2.0, size=(n_reads, read_len))
    q = np.clip(q, 2, 41).round().astype(np.int8)
    # Per-base error probability under the Phred definition.
    p_err = 10.0 ** (-q / 10.0)
    # Draw 'wrong-base?' Bernoulli per base.
    is_err = rng.random(q.shape) < p_err
    # Choose a base uniformly; if it's an error, pick a different one.
    base_idx = rng.integers(0, 4, size=q.shape)
    wrong_offset = rng.integers(1, 4, size=q.shape)
    base_idx = np.where(is_err, (base_idx + wrong_offset) % 4, base_idx)
    seq_arr = bases[base_idx]

    # Vectorised string build: join each row, format the 4-line record.
    lines: list[str] = []
    for i in range(n_reads):
        seq = "".join(seq_arr[i].tolist())
        qual = "".join(chr(int(v) + PHRED_OFFSET) for v in q[i])
        lines.append(f"@read_{i}")
        lines.append(seq)
        lines.append("+")
        lines.append(qual)
    return "\\n".join(lines) + "\\n"

fastq_text = build_synthetic_fastq()
print(f"fastq_text length: {len(fastq_text):,} chars")
print("First record preview:\\n" + "\\n".join(fastq_text.splitlines()[:4]))
'''

# ----- Step 2 ---------------------------------------------------------------

STEP2_WHY = (
    "Now parse it back. A FASTQ record is exactly four lines: header, "
    "sequence, `+`, quality. The quality string is ASCII-encoded Phred "
    "(`Q = ord(c) − 33`). Convert each character to its error "
    "probability `P_err = 10^(−Q/10)`. You should end up with a "
    "(N_READS, READ_LEN) int array of Phred scores and a matching "
    "float array of error probabilities."
)

STEP2_TODO = '''# TODO: parse fastq_text into arrays.
#
# Produce:
#   - seqs:   list[str] of length N_READS
#   - phred:  np.ndarray, shape (N_READS, READ_LEN), dtype int8
#   - p_err:  np.ndarray, shape (N_READS, READ_LEN), dtype float64
#
# Then verify by printing the mean Phred across all bases.

def parse_fastq(text: str):
    seqs: list[str] = []
    phred_rows: list[np.ndarray] = []
    # TODO: split text into 4-line chunks; pull out the seq + qual lines.
    return seqs, np.zeros((0, 0), dtype=np.int8)

seqs, phred = parse_fastq(fastq_text)
p_err = 10.0 ** (-phred.astype(float) / 10.0) if phred.size else np.zeros(0)
print(f"Parsed {len(seqs):,} reads")
print(f"Mean Phred: {phred.mean() if phred.size else float('nan'):.2f}")
print(f"Mean P_err: {p_err.mean() if p_err.size else float('nan'):.4f}")
'''

STEP2_SOL = '''def parse_fastq(text: str):
    lines = text.rstrip("\\n").split("\\n")
    assert len(lines) % 4 == 0, "FASTQ must have a multiple of 4 lines"
    n = len(lines) // 4
    seqs: list[str] = []
    phred_rows = np.empty((n, READ_LEN), dtype=np.int8)
    for i in range(n):
        header = lines[4 * i]
        seq = lines[4 * i + 1]
        plus = lines[4 * i + 2]
        qual = lines[4 * i + 3]
        assert header.startswith("@") and plus.startswith("+"), \\
            f"Malformed record at read {i}"
        assert len(seq) == len(qual) == READ_LEN, \\
            f"Length mismatch at read {i}"
        seqs.append(seq)
        phred_rows[i] = np.frombuffer(qual.encode("ascii"),
                                      dtype=np.uint8) - PHRED_OFFSET
    return seqs, phred_rows

seqs, phred = parse_fastq(fastq_text)
p_err = 10.0 ** (-phred.astype(float) / 10.0)
print(f"Parsed {len(seqs):,} reads")
print(f"Mean Phred: {phred.mean():.2f}")
print(f"Mean P_err: {p_err.mean():.4f}")
'''

# ----- Step 3 ---------------------------------------------------------------

STEP3_WHY = (
    "Per-read summary: the *mean Phred* of a read is a one-number "
    "QC score. Plot its histogram across the 10 000 reads. The expected "
    "number of errors per read is `λ = Σ P_err ≈ READ_LEN · ⟨P_err⟩`, "
    "which under independence is approximately Poisson(λ). Overlay the "
    "Poisson PMF on the observed *expected-error-count* histogram so "
    "you can see how the symbol-level model maps to the discrete-event "
    "model. (Heads-up: this instrument is good — λ ends up < 1, so the "
    "Poisson is dominated by the zero-error mode. That *is* the "
    "signature; it's why most reads are clean.)"
)

STEP3_TODO = '''# TODO: build two arrays and two plots.
#
#  (a) mean_phred_per_read: shape (N_READS,) — mean Phred per read.
#      Plot histogram.
#  (b) errors_per_read: shape (N_READS,) — count of *actually wrong*
#      bases per read is unknown (we don't carry a reference), so use
#      the *expected* count from the error model: round(p_err.sum(axis=1)).
#      Plot histogram and overlay the Poisson PMF with rate
#      λ = expected_errors_per_read.mean().

mean_phred_per_read = np.zeros(0)       # TODO
expected_errors_per_read = np.zeros(0)  # TODO
lam = 0.0                               # TODO

# Plot 1: mean Phred distribution.
plt.figure()
plt.title("Mean Phred per read (TODO)")
plt.xlabel("Mean Phred")
plt.ylabel("Count")
plt.show()

# Plot 2: expected errors per read with Poisson overlay.
plt.figure()
plt.title("Expected errors per read vs Poisson(λ) (TODO)")
plt.xlabel("Errors per read")
plt.ylabel("Density")
plt.show()
'''

STEP3_SOL = '''mean_phred_per_read = phred.mean(axis=1)
expected_errors_per_read = p_err.sum(axis=1)
lam = float(expected_errors_per_read.mean())

# Plot 1: mean Phred distribution.
plt.figure(figsize=(6, 4))
plt.hist(mean_phred_per_read, bins=40, color="#3b82f6",
         edgecolor="white")
plt.title("Mean Phred per read")
plt.xlabel("Mean Phred")
plt.ylabel("Count")
plt.tight_layout()
plt.show()

# Plot 2: expected errors per read with Poisson overlay.
obs = expected_errors_per_read.round().astype(int)
k_max = int(obs.max()) + 1
ks = np.arange(0, k_max + 1)
# Poisson PMF, computed in log space for stability.
log_pmf = (ks * math.log(lam) - lam
           - np.array([math.lgamma(k + 1) for k in ks]))
poisson_pmf = np.exp(log_pmf)

plt.figure(figsize=(6, 4))
plt.hist(obs, bins=np.arange(-0.5, k_max + 1.5, 1.0),
         density=True, color="#10b981", edgecolor="white",
         label="Observed (expected-error count)")
plt.plot(ks, poisson_pmf, "o-", color="#ef4444",
         label=f"Poisson(λ={lam:.2f})")
plt.title("Errors per read vs Poisson")
plt.xlabel("Errors per read")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.show()

print(f"λ (mean expected errors per read) = {lam:.3f}")
'''

# ----- Step 4 ---------------------------------------------------------------

STEP4_WHY = (
    "This is the money plot. For each cycle 0..149, average the Phred "
    "score across all reads and plot it as a curve. You should see the "
    "characteristic monotone decay from ~Q38 to ~Q20 — the instrument "
    "fingerprint. Overlay a horizontal line at Q=20 (the conventional "
    "trimming threshold, `P_err = 0.01`) to make the cutoff visible."
)

STEP4_TODO = '''# TODO: compute and plot the position-wise quality curve.
#
#   mean_phred_per_cycle: shape (READ_LEN,)
#   plot it vs cycle index 0..READ_LEN-1
#   draw a horizontal reference at Q = 20.

mean_phred_per_cycle = np.zeros(READ_LEN)  # TODO

plt.figure()
plt.title("Position-wise quality (TODO)")
plt.xlabel("Cycle")
plt.ylabel("Mean Phred")
plt.show()
'''

STEP4_SOL = '''mean_phred_per_cycle = phred.mean(axis=0)

plt.figure(figsize=(7, 4))
plt.plot(np.arange(READ_LEN), mean_phred_per_cycle,
         color="#6366f1", lw=2.0)
plt.axhline(20, color="#ef4444", ls="--", lw=1.0,
            label="Q=20 (P_err=0.01)")
plt.title("Position-wise mean Phred across 10,000 reads")
plt.xlabel("Cycle (position in read)")
plt.ylabel("Mean Phred")
plt.ylim(15, 42)
plt.legend()
plt.tight_layout()
plt.show()

print(f"Cycle 0  mean Phred: {mean_phred_per_cycle[0]:.2f}")
print(f"Cycle {READ_LEN - 1} mean Phred: {mean_phred_per_cycle[-1]:.2f}")
'''

# ----- Step 5 ---------------------------------------------------------------

STEP5_WHY = (
    "Quality trimming is a sliding-window threshold operation. Apply "
    "the classic *trim-from-3'-end at first base below Q10* policy "
    "(Q10 ⇒ P_err = 0.1, a generous threshold) per read. Drop reads "
    "shorter than 30 bp after trimming. Report:\n\n"
    "  - how many reads survive,\n"
    "  - the new median length,\n"
    "  - the new mean Phred across surviving bases."
)

STEP5_TODO = '''# TODO: implement Q10 3'-trim, drop reads < 30 bp post-trim.
#
# For each read: find the first cycle (scanning from the 3' end) where
# Phred drops below 10. Trim everything from that cycle to the end.
# Equivalent: trim at the *last* index from the left where Phred >= 10
# (i.e., new length = position of last good base + 1).
#
# Outputs:
#   trimmed_lengths: np.ndarray of length N_READS
#   survivors_mask: np.ndarray[bool], True where trimmed length >= 30
#   surviving_reads: int
#   median_len_after_trim: int
#   mean_phred_surviving: float

Q_TRIM = 10
MIN_LEN = 30

trimmed_lengths = np.zeros(N_READS, dtype=int)  # TODO
survivors_mask = np.ones(N_READS, dtype=bool)   # TODO
surviving_reads = 0                              # TODO
median_len_after_trim = 0                        # TODO
mean_phred_surviving = 0.0                       # TODO

print(f"Surviving reads:        {surviving_reads:,} / {N_READS:,}")
print(f"Median length post-trim: {median_len_after_trim}")
print(f"Mean Phred (surviving):  {mean_phred_surviving:.2f}")
'''

STEP5_SOL = '''Q_TRIM = 10
MIN_LEN = 30

# `good[i, j]` = True if base j of read i is at or above Q_TRIM.
good = phred >= Q_TRIM
# For each read, the trimmed length is (index of the last good base) + 1.
# If no base is good, the read collapses to length 0.
# Vectorised: reverse along axis=1, argmax finds the first True from the right.
has_good = good.any(axis=1)
reversed_good = good[:, ::-1]
last_good_from_end = reversed_good.argmax(axis=1)
trimmed_lengths = np.where(has_good, READ_LEN - last_good_from_end, 0)

survivors_mask = trimmed_lengths >= MIN_LEN
surviving_reads = int(survivors_mask.sum())
median_len_after_trim = int(np.median(trimmed_lengths[survivors_mask]))

# Mean Phred over surviving *bases* only.
total_q = 0.0
total_n = 0
for i in np.where(survivors_mask)[0]:
    L = trimmed_lengths[i]
    total_q += float(phred[i, :L].sum())
    total_n += int(L)
mean_phred_surviving = total_q / total_n

print(f"Surviving reads:        {surviving_reads:,} / {N_READS:,}")
print(f"Median length post-trim: {median_len_after_trim}")
print(f"Mean Phred (surviving):  {mean_phred_surviving:.2f}")

# Bonus: histogram of post-trim lengths.
plt.figure(figsize=(6, 4))
plt.hist(trimmed_lengths, bins=50, color="#f59e0b", edgecolor="white")
plt.axvline(MIN_LEN, color="#ef4444", ls="--", label=f"MIN_LEN={MIN_LEN}")
plt.title("Read length after Q10 3' trim")
plt.xlabel("Length")
plt.ylabel("Count")
plt.legend()
plt.tight_layout()
plt.show()
'''

# ----- Self-check -----------------------------------------------------------

SELFCHECK_CODE = '''# Self-check: ties together all five steps.
#
# Expected (deterministic given SEED=42, linear Q 38→20, σ=2 jitter):
#   - Step 1: fastq_text parses to 10,000 reads of length 150.
#   - Step 2: ⟨Phred⟩ ≈ 29 (midpoint of 38 and 20), so
#             ⟨P_err⟩ ≈ 10^(-29/10) ≈ 0.0013, but Jensen's inequality
#             pulls the *average* P_err up — empirically ~0.0027.
#   - Step 3: λ (expected errors per read) = READ_LEN · ⟨P_err⟩ ≈ 0.4.
#             Most reads have zero expected errors; that *is* the
#             Poisson signature.
#   - Step 4: cycle 0 mean Phred ≈ 38; cycle 149 mean Phred ≈ 20.
#   - Step 5: every read should survive — the curve bottoms out at
#             Q≈20, well above the Q10 cutoff. Median length = 150.

assert len(seqs) == N_READS, f"Expected {N_READS} reads, got {len(seqs)}"
assert phred.shape == (N_READS, READ_LEN), f"phred shape {phred.shape}"

mean_p_err = float(p_err.mean())
assert 0.0015 < mean_p_err < 0.005, (
    f"Mean P_err {mean_p_err:.4f} outside [0.0015, 0.005] — "
    "did the Q-curve change?"
)

q_cycle0 = float(mean_phred_per_cycle[0])
q_cycle_last = float(mean_phred_per_cycle[-1])
assert 36.0 < q_cycle0 < 40.0, f"Cycle 0 Phred {q_cycle0:.2f} not near 38"
assert 18.0 < q_cycle_last < 22.0, (
    f"Cycle {READ_LEN - 1} Phred {q_cycle_last:.2f} not near 20"
)

assert 0.15 <= lam <= 0.80, f"λ={lam:.2f} outside expected [0.15, 0.80]"

assert surviving_reads >= int(0.99 * N_READS), (
    f"Only {surviving_reads}/{N_READS} survived Q10 trim — too aggressive"
)
assert median_len_after_trim >= 140, (
    f"Median post-trim length {median_len_after_trim} < 140 — "
    "instrument curve degraded too far below Q10"
)

print("✅ Self-check passed.")
print(f"   ⟨P_err⟩ = {mean_p_err:.4f}")
print(f"   Cycle 0 / cycle {READ_LEN - 1} Phred = {q_cycle0:.2f} / {q_cycle_last:.2f}")
print(f"   λ (mean expected errors/read) = {lam:.2f}")
print(f"   {surviving_reads:,} / {N_READS:,} reads survived Q10 trim "
      f"(median len {median_len_after_trim})")
'''

# ----- EE framing -----------------------------------------------------------

EE_FRAMING_MD = """## EE framing

| FASTQ concept | EE analogue |
|---|---|
| Phred score `Q = −10 log₁₀ P_err` | **SNR in decibels.** Same definition, different domain. |
| Position-wise Q-curve | **Channel-response measurement** — flat at first, droops with cycle count due to dephasing (think filter rolloff). |
| Quality trimming at Q10 | **Pre-processing a noisy signal** before downstream correlation (alignment). |
| The Phred byte itself | **Lossy compression**: the full intensity trace → 1 base call + 1 confidence symbol. |

**Bottom line.** A sequencer is a noisy channel. FASTQ stores the
demodulated symbol *and* a per-symbol log-likelihood. Every downstream
step (alignment, variant calling, expression estimation) is a decoder
that consumes that confidence — which is why throwing Phred away is
throwing away a free Bayesian prior."""


# ---------------------------------------------------------------------------
# Assemble notebook
# ---------------------------------------------------------------------------

def build() -> nbformat.NotebookNode:
    nb = new_notebook()
    cells: list[nbformat.NotebookNode] = []

    cells.append(md(TITLE_MD))
    cells.append(md(AHA_MD))
    cells.append(code(PREAMBLE_CODE))
    cells.append(code(IMPORTS_CODE))

    cells.extend(step_block(1, 8, "Generate the synthetic FASTQ",
                            STEP1_WHY, STEP1_TODO, STEP1_SOL))
    cells.extend(step_block(2, 10, "Parse + decode Phred",
                            STEP2_WHY, STEP2_TODO, STEP2_SOL))
    cells.extend(step_block(3, 12, "Per-read quality histogram + Poisson",
                            STEP3_WHY, STEP3_TODO, STEP3_SOL))
    cells.extend(step_block(4, 15, "Position-wise quality curve",
                            STEP4_WHY, STEP4_TODO, STEP4_SOL))
    cells.extend(step_block(5, 15, "Quality-trim at Q10",
                            STEP5_WHY, STEP5_TODO, STEP5_SOL))

    cells.append(code(SELFCHECK_CODE))
    cells.append(md(EE_FRAMING_MD))

    nb["cells"] = cells

    # Notebook-level metadata: kernel + language info so Colab opens cleanly.
    nb["metadata"] = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python"},
        "colab": {"provenance": []},
    }
    return nb


def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(here, "exercise.ipynb")
    nb = build()
    with open(out_path, "w", encoding="utf-8") as fh:
        apply_colab_form(nb)
        nbformat.write(nb, fh)
    print(f"Wrote {out_path}")

    # Re-read and report a tiny structural summary.
    nb2 = nbformat.read(out_path, as_version=4)
    n_md = sum(1 for c in nb2.cells if c.cell_type == "markdown")
    n_code = sum(1 for c in nb2.cells if c.cell_type == "code")
    hidden = sum(
        1 for c in nb2.cells
        if c.cell_type == "code"
        and c.metadata.get("jupyter", {}).get("source_hidden") is True
    )
    print(f"Total cells: {len(nb2.cells)} ({n_md} markdown, {n_code} code, "
          f"{hidden} hidden code)")


if __name__ == "__main__":
    main()
