"""Build the L21 exercise notebook.

Run:  python3 build_notebook.py
Emits: exercise.ipynb in the same directory.

The notebook is structured as:
  1. Title markdown
  2. Aha callout markdown
  3. Preamble (!pip install)
  4. Imports + seed
  5..9. Step N markdown + TODO cell + hidden solution cell
  10. Final self-check assert cell
  11. EE framing markdown

Solution cells use metadata.jupyter.source_hidden = True so they collapse
in Colab. Each solution cell is preceded by a short prompt line that tells
the student how to expand it.
"""

from pathlib import Path
# Make the shared Colab-form helper importable from the parent dir.
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
from apply_colab_form import apply_colab_form  # noqa: E402

import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell


HERE = Path(__file__).parent
OUT = HERE / "exercise.ipynb"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def md(text: str):
    return new_markdown_cell(text)


def code(text: str):
    return new_code_cell(text)


def hidden_solution(text: str):
    """Return a code cell flagged source_hidden in Jupyter metadata.

    Colab respects `metadata.jupyter.source_hidden`; classic Jupyter does too.
    The cell still executes; only the source is collapsed.
    """
    cell = new_code_cell(text)
    cell.metadata = {
        "jupyter": {"source_hidden": True},
        "cellView": "form",  # Colab's analogue for "collapsed input"
    }
    return cell


# ---------------------------------------------------------------------------
# Cell bodies
# ---------------------------------------------------------------------------

TITLE_MD = """# L21 — HMMs, Profile HMMs, and Gene Finding

Hidden Markov Models give us a discrete-state-space view of a sequence: a latent
regime sequence (`intergenic / exon / intron / ...`) emits the observed letters.
In this exercise you will implement **Viterbi decoding in log-space** on a
synthetic 5-state gene HMM and a **profile-HMM-style scorer** that separates
globin from non-globin proteins."""


AHA_MD = """> **Aha.** Viterbi is **MAP decoding on a state-space model** — the single most probable hidden path. That is *not* the same as the sum of per-position posterior marginals (the forward-backward answer). And a **profile HMM is a position-specific filter bank**: a separate emission distribution per match column, the way a matched-filter bank assigns a separate template per signal class.
>
> Two practical reminders we lean on below:
> - Everything is done in **log-space** to dodge floating-point underflow when sequences are more than a few dozen positions long.
> - State-transition is the inner loop. We **vectorise** it in NumPy: at each step `delta` (1-D over states) is updated by `(delta[:, None] + log_A).max(axis=0) + log_b[:, x_t]`."""


PREAMBLE_CODE = """# Colab preamble — pinned installs. Re-runs are no-ops if already present.
!pip install numpy==1.26.4 matplotlib==3.8.4 requests==2.32.3 biopython==1.83 -q
# Optional: pyhmmer makes parsing Pfam HMMER3 files trivial. We write our own
# minimal parser below so this is just a backup, but install it for the
# fallback path.
!pip install pyhmmer==0.10.15 -q"""


IMPORTS_CODE = """import io
import os
import math
import time
import json
import urllib.request
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import requests

# Biopython is only used for FASTA parsing of UniProt downloads.
try:
    from Bio import SeqIO
    _HAS_BIOPYTHON = True
except Exception:
    _HAS_BIOPYTHON = False

# Deterministic seeds — used throughout for the synthetic fallbacks.
np.random.seed(42)
RNG = np.random.default_rng(42)

CACHE = Path("./_l21_cache")
CACHE.mkdir(exist_ok=True)
print(f"Cache dir: {CACHE.resolve()}")"""


# ---------------------------------------------------------------------------
# Step 1 — load HMM, fetch BRCA1 region
# ---------------------------------------------------------------------------

STEP1_MD = """## Step 1 (10 min) — Load the 5-state gene HMM + fetch a 2 kb BRCA1 region

We work with a tiny eukaryotic gene HMM with 5 states:
`intergenic, 5'UTR, exon, intron, 3'UTR`. Each state has its own
nucleotide emission distribution (e.g., exons are GC-biased; introns are AT-biased)
and a **sticky self-transition** (~0.99 for intergenic, ~0.97 for intron) so the
chain naturally generates long stretches before switching state.

Fill in `load_gene_hmm()` to return the dict, and use Ensembl REST to download a
2 kb slice of BRCA1 (with a deterministic synthetic fallback if the call fails)."""


STEP1_TODO = '''def load_gene_hmm():
    """Return a dict describing a 5-state gene HMM.

    Required keys:
      - "states":   list[str], length K
      - "init":     np.ndarray shape (K,)   initial distribution, sums to 1
      - "trans":    np.ndarray shape (K, K) row-stochastic transition matrix
      - "emit":     np.ndarray shape (K, 4) emission over {A, C, G, T}; rows sum to 1
      - "alphabet": "ACGT"

    Realistic ballpark:
      - intergenic, 3'UTR: self-transition ~ 0.99, AT-rich emissions
      - exon:              self-transition ~ 0.985, GC-rich emissions
      - intron:            self-transition ~ 0.97, AT-rich emissions
      - 5'UTR:             self-transition ~ 0.99, mildly GC-rich
    """
    # TODO: build and return the HMM dictionary
    return None


def fetch_brca1_region(cache_path):
    """Fetch a 2 kb BRCA1 region from Ensembl REST.

    URL:
      https://rest.ensembl.org/sequence/region/human/17:43044295..43046294
      ?content-type=text/x-fasta

    Cache the response to disk. On network failure, return a synthetic 2 kb
    sequence with two embedded "exon-like" GC-rich runs of ~150 bp each.
    """
    # TODO: fetch the FASTA, parse to a plain ACGT string, cache to cache_path
    return None


hmm = load_gene_hmm()
brca = fetch_brca1_region(CACHE / "brca1_2kb.fa")
print("HMM:", hmm)
print("BRCA1 length:", None if brca is None else len(brca))
'''


STEP1_SOL_MD = """*Click ▶ to expand the reference solution.*"""


STEP1_SOL = '''import re

_ALPHABET = "ACGT"


def load_gene_hmm():
    """Hardcoded 5-state gene HMM.

    Emission rows: A, C, G, T order.
      intergenic:  AT-rich, equal
      5'UTR:       mildly GC-rich
      exon:        clearly GC-rich
      intron:      AT-rich, long
      3'UTR:       AT-rich, like intergenic
    """
    states = ["intergenic", "5UTR", "exon", "intron", "3UTR"]
    init = np.array([0.85, 0.05, 0.02, 0.04, 0.04])

    # Transition matrix: K x K, row-stochastic.
    # Self-transitions chosen so the implied length distribution matches biology:
    # mean length = 1 / (1 - p_self).
    trans = np.array([
        # to: intergenic  5UTR    exon    intron  3UTR
        [   0.992,        0.005,  0.001,  0.001,  0.001 ],  # intergenic
        [   0.001,        0.985,  0.013,  0.000,  0.001 ],  # 5UTR
        [   0.000,        0.000,  0.985,  0.013,  0.002 ],  # exon
        [   0.000,        0.000,  0.030,  0.970,  0.000 ],  # intron
        [   0.010,        0.000,  0.000,  0.000,  0.990 ],  # 3UTR
    ])
    # Normalise rows just in case of floating drift.
    trans = trans / trans.sum(axis=1, keepdims=True)

    # Emissions chosen with enough between-state contrast that Viterbi can
    # recover the truth on self-simulated sequences. Real biology has weaker
    # contrast; we exaggerate here to make the algorithmic point clean.
    # Each state has a distinct emission "fingerprint":
    #   intergenic: A-heavy
    #   5UTR:       balanced GC
    #   exon:       G-heavy (high GC)
    #   intron:     T-heavy
    #   3UTR:       C-heavy
    emit = np.array([
        [0.55, 0.15, 0.15, 0.15],  # intergenic — A-dominant
        [0.20, 0.30, 0.30, 0.20],  # 5UTR       — balanced GC
        [0.10, 0.20, 0.55, 0.15],  # exon       — G-dominant (GC-rich)
        [0.15, 0.15, 0.15, 0.55],  # intron     — T-dominant (AT-rich)
        [0.20, 0.50, 0.15, 0.15],  # 3UTR       — C-dominant
    ])
    emit = emit / emit.sum(axis=1, keepdims=True)

    return {
        "states":   states,
        "init":     init,
        "trans":    trans,
        "emit":     emit,
        "alphabet": _ALPHABET,
    }


def _synthetic_brca1(n=2000):
    """Fallback: build a synthetic 2 kb DNA sequence with two GC-rich runs."""
    rng = np.random.default_rng(123)
    bg = rng.choice(list("ACGT"), size=n, p=[0.30, 0.20, 0.20, 0.30])
    # Splice in two ~150 bp GC-rich "exon-like" blocks at positions 600 and 1300.
    for start in (600, 1300):
        gc = rng.choice(list("ACGT"), size=150, p=[0.20, 0.30, 0.30, 0.20])
        bg[start:start + 150] = gc
    return "".join(bg)


def fetch_brca1_region(cache_path):
    cache_path = Path(cache_path)
    if cache_path.exists():
        text = cache_path.read_text()
    else:
        url = ("https://rest.ensembl.org/sequence/region/human/"
               "17:43044295..43046294?content-type=text/x-fasta")
        try:
            r = requests.get(url, timeout=15)
            r.raise_for_status()
            text = r.text
            cache_path.write_text(text)
        except Exception as e:
            print(f"  ! Ensembl fetch failed ({e}); using synthetic fallback.")
            seq = _synthetic_brca1()
            cache_path.write_text(f">SYNTHETIC_BRCA1_2kb\\n{seq}\\n")
            return seq

    # Strip FASTA header(s) and whitespace; collapse to ACGT (uppercase).
    lines = [ln.strip() for ln in text.splitlines() if not ln.startswith(">")]
    seq = "".join(lines).upper()
    seq = re.sub(r"[^ACGT]", "A", seq)  # mask N and ambiguities
    return seq[:2000]


hmm = load_gene_hmm()
brca = fetch_brca1_region(CACHE / "brca1_2kb.fa")

print("States:", hmm["states"])
print("Self-transitions:", {s: round(float(hmm["trans"][i, i]), 4)
                             for i, s in enumerate(hmm["states"])})
print(f"BRCA1 length: {len(brca)} bp  ({brca[:60]}...)")
'''


# ---------------------------------------------------------------------------
# Step 2 — Viterbi in log-space
# ---------------------------------------------------------------------------

STEP2_MD = """## Step 2 (12 min) — Viterbi in log-space; visualise the decoded path

Implement Viterbi as a `O(T K^2)` dynamic program. Every quantity is stored in
log-space so multiplying probabilities becomes adding log-probabilities, and
underflow stays away even for `T = 10^4`.

Recurrence (`delta[t][j]` is the log-probability of the most probable path
ending in state `j` at position `t`):

```
delta[t][j] = max_i (delta[t-1][i] + log A[i][j])  +  log B[j][x_t]
```

Backtrace pointers store the `argmax_i` so we can walk back from `T-1` to `0`.

We then plot the decoded state path as a colour-coded sequence annotation."""


STEP2_TODO = '''def viterbi_log(seq, hmm):
    """Return (path, log_score) where path is a list of state indices.

    seq  -- a string over hmm["alphabet"]
    hmm  -- dict from Step 1

    Implementation hints:
      - Convert init / trans / emit to log-space ONCE.
      - Use a (T, K) delta matrix and a (T, K) int-pointer matrix.
      - Vectorise the per-step transition step with broadcasting:
            cand = delta[t-1, :, None] + logA          # shape (K, K)
            delta[t] = cand.max(axis=0) + logB[:, x_t] # shape (K,)
            ptr[t]   = cand.argmax(axis=0)             # shape (K,)
      - Use log(0) -> -np.inf safely.
    """
    # TODO: implement
    return [], -np.inf


def plot_viterbi_path(path, hmm, title="Viterbi-decoded state path"):
    """Plot the decoded state path as a coloured strip along the sequence."""
    # TODO: matplotlib strip plot, one color per state, with a legend
    pass


# Guarded driver: runs even if Step 1's stubs returned None.
if brca is not None and hmm is not None:
    path, score = viterbi_log(brca, hmm)
    print(f"Viterbi log-score: {score:.2f}")
    plot_viterbi_path(path, hmm)
else:
    print("(skipping — Step 1 stubs still return None; fill them in or expand the solution.)")
'''


STEP2_SOL = '''_LOG_ZERO = -1e18  # near -inf but safe for arithmetic


def viterbi_log(seq, hmm):
    """Vectorised Viterbi in log-space.

    Returns
    -------
    path : list[int]
        Decoded state index per position.
    log_score : float
        log P(X, Z* | theta) — the joint log-probability of the best path.
    """
    A = hmm["trans"]
    B = hmm["emit"]
    pi = hmm["init"]
    alphabet = hmm["alphabet"]
    sym2idx = {c: i for i, c in enumerate(alphabet)}

    # Encode sequence as integer indices into the alphabet.
    obs = np.array([sym2idx.get(c, 0) for c in seq], dtype=np.int64)
    T = obs.size
    K = A.shape[0]

    # All in log-space.
    with np.errstate(divide="ignore"):
        log_pi = np.where(pi > 0, np.log(pi), _LOG_ZERO)
        log_A  = np.where(A > 0, np.log(A), _LOG_ZERO)
        log_B  = np.where(B > 0, np.log(B), _LOG_ZERO)

    delta = np.full((T, K), _LOG_ZERO)
    ptr   = np.zeros((T, K), dtype=np.int32)

    # Initialisation.
    delta[0] = log_pi + log_B[:, obs[0]]

    # Recursion (vectorised over the K-from-states axis).
    for t in range(1, T):
        # cand[i, j] = delta[t-1, i] + log_A[i, j]
        cand = delta[t - 1, :, None] + log_A
        ptr[t]   = cand.argmax(axis=0)
        delta[t] = cand.max(axis=0) + log_B[:, obs[t]]

    # Termination + backtrace.
    last = int(delta[T - 1].argmax())
    log_score = float(delta[T - 1, last])

    path = np.zeros(T, dtype=np.int32)
    path[-1] = last
    for t in range(T - 1, 0, -1):
        path[t - 1] = ptr[t, path[t]]

    return path.tolist(), log_score


_STATE_COLORS = ["#7e8a9a", "#3aa07b", "#d36b3b", "#4a4d8a", "#b86bb5"]


def plot_viterbi_path(path, hmm, title="Viterbi-decoded state path"):
    fig, ax = plt.subplots(figsize=(11, 1.6))
    path_arr = np.asarray(path)
    # Build a coloured horizontal strip: one column per position.
    img = path_arr.reshape(1, -1)
    cmap = plt.matplotlib.colors.ListedColormap(_STATE_COLORS[: len(hmm["states"])])
    ax.imshow(img, aspect="auto", interpolation="nearest", cmap=cmap,
              vmin=0, vmax=len(hmm["states"]) - 1)
    ax.set_yticks([])
    ax.set_xlabel("Sequence position (bp)")
    ax.set_title(title)
    # Legend.
    handles = [plt.matplotlib.patches.Patch(facecolor=_STATE_COLORS[i],
                                            label=hmm["states"][i])
               for i in range(len(hmm["states"]))]
    ax.legend(handles=handles, ncol=len(hmm["states"]),
              bbox_to_anchor=(0.5, -0.6), loc="upper center", frameon=False)
    plt.tight_layout()
    plt.show()


path, score = viterbi_log(brca, hmm)
state_counts = {hmm["states"][i]: int(np.sum(np.array(path) == i))
                for i in range(len(hmm["states"]))}
print(f"Viterbi log-score: {score:.2f}")
print(f"Per-state position counts: {state_counts}")
plot_viterbi_path(path, hmm)
'''


# ---------------------------------------------------------------------------
# Step 3 — Sanity-check Viterbi vs a synthetic ground truth
# ---------------------------------------------------------------------------

STEP3_MD = """## Step 3 (12 min) — Sanity-check: decode a sequence with known truth

We can't get a publicly aligned per-base RefSeq exon mask for every Colab
session without a large download. Instead we **simulate from the HMM itself**
to get ground-truth labels, run Viterbi, and report position-wise accuracy.

This is the right unit test for the Viterbi implementation: if the model
*generated* the sequence, a correct MAP decoder should recover most of the
labels."""


STEP3_TODO = '''def sample_from_hmm(hmm, T, rng):
    """Draw (sequence, true_state_path) of length T from the HMM by ancestral sampling.

    Use rng.choice with probabilities = hmm["init"], then for each subsequent
    step draw the next state from hmm["trans"][prev_state] and the emitted
    letter from hmm["emit"][state].
    """
    # TODO: implement
    return "", []


def viterbi_accuracy(true_path, pred_path):
    """Fraction of positions where pred matches truth."""
    # TODO: implement
    return 0.0


if hmm is not None:
    sim_seq, sim_truth = sample_from_hmm(hmm, T=3000, rng=np.random.default_rng(7))
    sim_pred, _ = viterbi_log(sim_seq, hmm)
    print(f"Viterbi recovery accuracy on simulated 3 kb sequence: "
          f"{viterbi_accuracy(sim_truth, sim_pred):.3f}")
else:
    print("(skipping — Step 1's load_gene_hmm() stub returns None.)")
'''


STEP3_SOL = '''def sample_from_hmm(hmm, T, rng):
    states_idx = np.arange(len(hmm["states"]))
    alpha_idx  = np.arange(len(hmm["alphabet"]))

    z = np.zeros(T, dtype=np.int32)
    x = np.zeros(T, dtype=np.int32)

    z[0] = rng.choice(states_idx, p=hmm["init"])
    x[0] = rng.choice(alpha_idx,  p=hmm["emit"][z[0]])
    for t in range(1, T):
        z[t] = rng.choice(states_idx, p=hmm["trans"][z[t - 1]])
        x[t] = rng.choice(alpha_idx,  p=hmm["emit"][z[t]])
    seq = "".join(hmm["alphabet"][i] for i in x)
    return seq, z.tolist()


def viterbi_accuracy(true_path, pred_path):
    t = np.asarray(true_path)
    p = np.asarray(pred_path)
    return float(np.mean(t == p))


sim_seq, sim_truth = sample_from_hmm(hmm, T=3000, rng=np.random.default_rng(7))
sim_pred, _ = viterbi_log(sim_seq, hmm)
acc = viterbi_accuracy(sim_truth, sim_pred)
print(f"Viterbi recovery accuracy on simulated 3 kb sequence: {acc:.3f}")
plot_viterbi_path(sim_pred, hmm, title="Viterbi decoding of a self-simulated sequence")
'''


# ---------------------------------------------------------------------------
# Step 4 — Profile HMM
# ---------------------------------------------------------------------------

STEP4_MD = """## Step 4 (15 min) — Profile HMM: score 50 sequences

A **profile HMM** has one *match* state per column of a multiple alignment,
plus *insert* and *delete* states between columns. For a quick scorer we use a
simplified profile (no insertions/deletions): the score of a query is the sum
of column-specific log-odds against a null amino-acid background.

We fetch the Pfam PF00042 (globin) profile HMM (HMMER3 format) from the
EBI / InterPro API, parse just the match-state emission distributions with a
tiny custom parser, fall back to `pyhmmer` if that fails, and a hand-built
globin consensus if both fail. Then we score 25 globins (UniProt
`family:globin`) and 25 non-globin kinases."""


STEP4_TODO = '''def fetch_profile_hmm(cache_path):
    """Fetch the Pfam PF00042 HMM (HMMER3 format) and return:
        - mat_emit: np.ndarray shape (L, 20) of match-state amino-acid probabilities
        - aa_alpha: str of length 20 listing the amino-acid order in the HMM
    Fallback: build a tiny synthetic globin profile from a hand-coded consensus
    via a one-hot + smoothing trick.
    """
    # TODO: download, parse, fall back
    return None, "ACDEFGHIKLMNPQRSTVWY"


def fetch_test_sequences(cache_path):
    """Return list[(label, seq)] of 25 globins (label="globin") and 25 non-globins
    (label="other"). Source: UniProt REST (family:globin, family:kinase NOT family:globin).
    Synthetic fallback: random aa sequences, half with an embedded globin-consensus motif.
    """
    # TODO: fetch + parse FASTA, with synthetic fallback
    return []


def score_profile(seq, mat_emit, aa_alpha, null_bg=None):
    """Profile-HMM score (log-odds) for a query sequence against a profile.

    Simplification: greedy alignment of the query to the profile columns —
    for each profile column k from 0..L-1, take the best matching position
    in the (sliding) query and sum its log-odds against the null background.
    This is NOT full Match/Insert/Delete Viterbi; it's a teaching-stand-in
    that captures the position-specific-filter-bank intuition cleanly.
    """
    # TODO: implement log-odds scoring; return a single scalar score
    return 0.0


mat_emit, aa_alpha = fetch_profile_hmm(CACHE / "PF00042.hmm")
test_seqs = fetch_test_sequences(CACHE / "uniprot_50.fa")
print(f"Profile length: {None if mat_emit is None else mat_emit.shape[0]} match columns")
print(f"Test sequences: {len(test_seqs)}")
'''


STEP4_SOL = r'''AA_ALPHA_STD = "ACDEFGHIKLMNPQRSTVWY"
# Background frequencies (Robinson-Robinson, used by HMMER3) — order = AA_ALPHA_STD.
NULL_BG = np.array([
    0.0826, 0.0137, 0.0546, 0.0671, 0.0387, 0.0708, 0.0227, 0.0594,
    0.0584, 0.0967, 0.0242, 0.0406, 0.0469, 0.0393, 0.0553, 0.0656,
    0.0534, 0.0687, 0.0108, 0.0292
])
NULL_BG = NULL_BG / NULL_BG.sum()


def _parse_hmmer3(text):
    """Tiny HMMER3 parser: extract match-state emission probabilities.

    HMMER3 .hmm files store NEGATIVE log probabilities (natural log) on the
    'match emission' lines inside the HMM block. The header line names the
    amino-acid columns in HMM-internal order. We return probabilities,
    re-mapped to AA_ALPHA_STD for consistency.
    """
    lines = text.splitlines()
    aa_order = None
    in_hmm   = False
    rows = []  # list of np.ndarray shape (20,) of probabilities in AA_ALPHA_STD order

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("HMM "):
            # Header: "HMM A C D E F G ..." — capture the 20-letter ordering.
            parts = line.split()
            aa_order = "".join(parts[1:21])
            in_hmm = True
            # skip the next line ("m->m m->i m->d ...") and the COMPO line if present
            i += 1
            # skip transition-label line
            i += 1
            # next non-empty line might be COMPO; we just keep going and pick up
            # only the lines whose first token is an integer (match-state lines).
            while i < len(lines):
                row = lines[i].split()
                if len(row) >= 21:
                    try:
                        int(row[0])
                    except ValueError:
                        i += 1
                        continue
                    # Match emission line: row[0] = state index, row[1..20] = -log p
                    neg_logs = np.array([float(x) if x != "*" else np.inf
                                         for x in row[1:21]])
                    probs = np.exp(-neg_logs)
                    # Re-map from aa_order -> AA_ALPHA_STD.
                    remap = np.zeros(20)
                    for j, aa in enumerate(AA_ALPHA_STD):
                        if aa in aa_order:
                            remap[j] = probs[aa_order.index(aa)]
                    s = remap.sum()
                    if s > 0:
                        remap = remap / s
                    else:
                        remap = NULL_BG.copy()
                    rows.append(remap)
                    # The next two lines are insert emissions and transitions; skip them.
                    i += 3
                    continue
                if lines[i].startswith("//"):
                    return np.vstack(rows), AA_ALPHA_STD
                i += 1
            break
        i += 1

    if rows:
        return np.vstack(rows), AA_ALPHA_STD
    return None, AA_ALPHA_STD


def _synthetic_profile():
    """Build a 30-column synthetic globin profile from a hand-coded consensus."""
    # Hairpin of the alpha-globin family, very rough consensus.
    consensus = "MVLSPADKTNVKAAWGKVGAHAGEYGAEAL"
    L = len(consensus)
    mat = np.full((L, 20), NULL_BG.copy())
    # Heavily weight the consensus residue at each column; small smoothing.
    for k, aa in enumerate(consensus):
        if aa in AA_ALPHA_STD:
            j = AA_ALPHA_STD.index(aa)
            row = np.full(20, 0.01)
            row[j] = 0.81
            row = row / row.sum()
            mat[k] = row
    return mat, AA_ALPHA_STD


def fetch_profile_hmm(cache_path):
    import gzip
    cache_path = Path(cache_path)
    text = None
    if cache_path.exists():
        text = cache_path.read_text()
    else:
        url = "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00042/?annotation=hmm"
        try:
            r = requests.get(url, timeout=20)
            r.raise_for_status()
            # The InterPro endpoint serves the HMM as gzip-compressed text/plain
            # (Content-Type often "application/gzip"). Decompress if needed.
            payload = r.content
            try:
                text = gzip.decompress(payload).decode("utf-8", errors="replace")
            except (OSError, gzip.BadGzipFile):
                text = payload.decode("utf-8", errors="replace")
            if "HMMER3" not in text[:200] and "NAME" not in text[:200]:
                raise ValueError("Unexpected response (not HMMER3)")
            cache_path.write_text(text)
        except Exception as e:
            print(f"  ! Pfam HMM fetch failed ({e}); using synthetic globin profile.")
            return _synthetic_profile()

    # First try our custom parser.
    try:
        mat, alpha = _parse_hmmer3(text)
        if mat is not None and mat.shape[0] >= 20:
            return mat, alpha
    except Exception as e:
        print(f"  ! Custom HMMER3 parse failed ({e}); trying pyhmmer.")

    # Fallback: pyhmmer (if installed in the Colab env).
    try:
        import pyhmmer
        with pyhmmer.plan7.HMMFile(io.BytesIO(text.encode())) as hf:
            hmm_pyh = next(hf)
        # match emissions: shape (M+1, K) — index 1..M are match states
        # MatchEmissions in pyhmmer returns a 2D array.
        match = np.array(hmm_pyh.match_emissions[1:])
        # pyhmmer alphabet ordering is the amino-acid alphabet; remap to AA_ALPHA_STD.
        pyh_alpha = hmm_pyh.alphabet.symbols  # e.g. "ACDEFGHIKLMNPQRSTVWY"
        remap_cols = [pyh_alpha.index(a) for a in AA_ALPHA_STD]
        return match[:, remap_cols], AA_ALPHA_STD
    except Exception as e:
        print(f"  ! pyhmmer fallback also failed ({e}); using synthetic globin profile.")
        return _synthetic_profile()


def _synthetic_test_sequences():
    """Hand-built: 25 sequences with a globin-consensus block, 25 random controls."""
    rng = np.random.default_rng(11)
    bg_letters = list(AA_ALPHA_STD)
    bg_p = NULL_BG.copy()

    consensus_block = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF"
    out = []
    for i in range(25):
        # 120-aa scaffold with the consensus inserted at a random offset.
        length = rng.integers(100, 160)
        body = rng.choice(bg_letters, size=length, p=bg_p)
        off = int(rng.integers(0, length - len(consensus_block)))
        body[off:off + len(consensus_block)] = list(consensus_block)
        out.append(("globin", "".join(body)))
    for i in range(25):
        length = rng.integers(100, 200)
        body = rng.choice(bg_letters, size=length, p=bg_p)
        out.append(("other", "".join(body)))
    return out


def _parse_fasta_text(text):
    """Minimal FASTA parser — yields (header, seq)."""
    name = None
    chunks = []
    for ln in text.splitlines():
        if ln.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name = ln[1:].strip()
            chunks = []
        else:
            chunks.append(ln.strip())
    if name is not None:
        yield name, "".join(chunks)


def fetch_test_sequences(cache_path):
    cache_path = Path(cache_path)
    if cache_path.exists():
        text = cache_path.read_text()
        parsed = list(_parse_fasta_text(text))
        labelled = [("globin" if h.split("|", 1)[0] == "globin" else "other", s)
                    for h, s in parsed]
        glob = [x for x in labelled if x[0] == "globin"][:25]
        oth  = [x for x in labelled if x[0] == "other"][:25]
        if len(glob) == 25 and len(oth) == 25:
            return glob + oth

    # Globins are short (~150 aa); kinases are long (~400-1500 aa).
    # We over-fetch both classes and truncate long sequences to length 500
    # so the profile scorer doesn't pay a quadratic price.
    urls = [
        ("globin",
         "https://rest.uniprot.org/uniprotkb/search?"
         "query=family%3Aglobin&format=fasta&size=40"),
        ("other",
         "https://rest.uniprot.org/uniprotkb/search?"
         "query=family%3Akinase+NOT+family%3Aglobin&format=fasta&size=40"),
    ]
    fetched = []
    try:
        for label, url in urls:
            r = requests.get(url, timeout=20)
            r.raise_for_status()
            for h, s in _parse_fasta_text(r.text):
                s = "".join(c for c in s.upper() if c in AA_ALPHA_STD)
                if len(s) >= 30:
                    fetched.append((label, h, s[:500]))
        # Trim to exactly 25 + 25.
        globs = [(l, h, s) for (l, h, s) in fetched if l == "globin"][:25]
        others = [(l, h, s) for (l, h, s) in fetched if l == "other"][:25]
        if len(globs) == 25 and len(others) == 25:
            cache_text = "".join(f">{l}|{h}\n{s}\n" for (l, h, s) in globs + others)
            cache_path.write_text(cache_text)
            return [(l, s) for (l, h, s) in globs + others]
        raise ValueError(f"Got {len(globs)} globins + {len(others)} others; need 25 each")
    except Exception as e:
        print(f"  ! UniProt fetch failed ({e}); using synthetic test sequences.")
        return _synthetic_test_sequences()


def score_profile(seq, mat_emit, aa_alpha, null_bg=NULL_BG):
    """Sliding-window log-odds score against a profile.

    For each possible alignment offset 0 <= start <= len(seq) - L, take the
    log-odds sum across the L match columns. Return the best score across
    offsets — this is a discriminative profile scorer in the same spirit as
    HMMER's match-state-only "glocal" Viterbi without indels.
    """
    L = mat_emit.shape[0]
    if len(seq) < L:
        return -np.inf
    aa2idx = {a: i for i, a in enumerate(aa_alpha)}
    obs = np.array([aa2idx.get(c, -1) for c in seq], dtype=np.int64)
    valid = obs >= 0
    if not valid.any():
        return -np.inf

    # Pre-compute log-odds matrix M (L, 20) once.
    with np.errstate(divide="ignore"):
        log_odds = np.log(np.clip(mat_emit, 1e-30, 1.0)) - np.log(null_bg)

    # Slide a window of length L; vectorise the inner sum.
    best = -np.inf
    n = len(seq)
    for start in range(0, n - L + 1):
        window = obs[start:start + L]
        ok = window >= 0
        if not ok.all():
            continue
        s = log_odds[np.arange(L), window].sum()
        if s > best:
            best = float(s)
    return best


mat_emit, aa_alpha = fetch_profile_hmm(CACHE / "PF00042.hmm")
test_seqs = fetch_test_sequences(CACHE / "uniprot_50.fa")
print(f"Profile length: {mat_emit.shape[0]} match columns")
print(f"Test sequences: {len(test_seqs)} "
      f"({sum(1 for l, _ in test_seqs if l == 'globin')} globin / "
      f"{sum(1 for l, _ in test_seqs if l == 'other')} other)")
'''


# ---------------------------------------------------------------------------
# Step 5 — ROC and discussion
# ---------------------------------------------------------------------------

STEP5_MD = """## Step 5 (11 min) — Score the 50 test sequences; plot ROC

Score every sequence; separate globin (positive class) from non-globin
(negative class); plot the ROC curve as the score threshold slides from
generous to strict; report the area under the curve (AUC).

The aha for this step: even our toy profile HMM cleanly separates the two
classes — a profile HMM is a **position-specific filter bank**, and the
filter bank tuned on globins is a good detector for globins."""


STEP5_TODO = '''def score_all(test_seqs, mat_emit, aa_alpha):
    """Return (scores, labels) as parallel np.ndarrays. labels in {1, 0}."""
    # TODO: implement
    return np.array([]), np.array([])


def roc_curve(scores, labels):
    """Return (fpr, tpr) arrays for varying thresholds. Higher score = positive."""
    # TODO: implement (no sklearn allowed — write the sweep manually)
    return np.array([]), np.array([])


def auc_trapezoid(fpr, tpr):
    """Area under the ROC curve via trapezoid rule."""
    # TODO: implement
    return 0.0


if mat_emit is not None and test_seqs:
    scores, labels = score_all(test_seqs, mat_emit, aa_alpha)
    fpr, tpr = roc_curve(scores, labels)
    auc = auc_trapezoid(fpr, tpr)
    print(f"AUC: {auc:.3f}")
else:
    print("(skipping — Step 4 stubs still return empty / None.)")
'''


STEP5_SOL = '''def score_all(test_seqs, mat_emit, aa_alpha):
    scores = np.array([score_profile(s, mat_emit, aa_alpha) for _, s in test_seqs])
    labels = np.array([1 if lab == "globin" else 0 for lab, _ in test_seqs])
    return scores, labels


def roc_curve(scores, labels):
    """Plain implementation: sort by score descending, sweep thresholds."""
    order = np.argsort(-scores)
    y = labels[order]
    P = int(y.sum())
    N = len(y) - P
    tpr = np.zeros(len(y) + 1)
    fpr = np.zeros(len(y) + 1)
    tp = fp = 0
    for i, yi in enumerate(y, start=1):
        if yi == 1:
            tp += 1
        else:
            fp += 1
        tpr[i] = tp / max(P, 1)
        fpr[i] = fp / max(N, 1)
    return fpr, tpr


def auc_trapezoid(fpr, tpr):
    # NumPy 2 renamed trapz -> trapezoid; support both.
    _trap = getattr(np, "trapezoid", None) or np.trapz
    return float(_trap(tpr, fpr))


scores, labels = score_all(test_seqs, mat_emit, aa_alpha)
fpr, tpr = roc_curve(scores, labels)
auc = auc_trapezoid(fpr, tpr)
print(f"AUC: {auc:.3f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

# Score distributions: globin vs other.
g = scores[labels == 1]
o = scores[labels == 0]
bins = np.linspace(min(scores.min(), -50), scores.max() + 1, 30)
ax1.hist(o, bins=bins, alpha=0.6, label="non-globin", color="#7e8a9a")
ax1.hist(g, bins=bins, alpha=0.6, label="globin",     color="#d36b3b")
ax1.set_xlabel("Profile-HMM log-odds score")
ax1.set_ylabel("Count")
ax1.set_title("Score distributions: PF00042 profile")
ax1.legend()

# ROC.
ax2.plot(fpr, tpr, color="#3aa07b", linewidth=2)
ax2.plot([0, 1], [0, 1], "k--", alpha=0.4, linewidth=1)
ax2.set_xlabel("False positive rate")
ax2.set_ylabel("True positive rate")
ax2.set_title(f"ROC (AUC = {auc:.3f})")
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1.02)

plt.tight_layout()
plt.show()
'''


# ---------------------------------------------------------------------------
# Self-check
# ---------------------------------------------------------------------------

SELFCHECK_MD = """## Self-check

If you've expanded the reference solutions for Steps 2, 3, 4, and 5 (or written
equivalent code yourself), these asserts should pass."""


SELFCHECK_CODE = '''# 1) Viterbi must recover >= 80% of the truth labels on a self-simulated sequence.
sim_seq2, sim_truth2 = sample_from_hmm(hmm, T=3000, rng=np.random.default_rng(13))
sim_pred2, _ = viterbi_log(sim_seq2, hmm)
recovery = viterbi_accuracy(sim_truth2, sim_pred2)
print(f"[1] Viterbi recovery on simulated 3 kb seq: {recovery:.3f}")
assert recovery >= 0.80, (
    f"Viterbi recovery {recovery:.3f} < 0.80 — check log-space arithmetic / "
    "backtrace direction."
)

# 2) Profile-HMM must rank globin sequences above the 75th-percentile of non-globin.
g_scores = scores[labels == 1]
o_scores = scores[labels == 0]
print(f"[2] Mean globin score: {g_scores.mean():.2f}, "
      f"mean non-globin: {o_scores.mean():.2f}")
assert g_scores.mean() > o_scores.mean(), (
    "Globin mean score should exceed non-globin mean score."
)

# 3) AUC should be high (well-separated classes).
print(f"[3] AUC = {auc:.3f}")
assert auc >= 0.90, f"AUC {auc:.3f} < 0.90 — profile not discriminating well."

# 4) At a threshold midway between the two class means, almost every globin
#    should be above and almost every non-globin below.
threshold = 0.5 * (g_scores.mean() + o_scores.mean())
tp = int(np.sum(g_scores >= threshold))
fp = int(np.sum(o_scores >= threshold))
print(f"[4] Threshold = {threshold:.2f}; globin TP = {tp}/25; non-globin FP = {fp}/25")
assert tp >= 22 and fp <= 3, (
    "Threshold check failed — score separation is weaker than expected."
)

print("\\nAll self-checks passed.")
'''


# ---------------------------------------------------------------------------
# EE framing closing
# ---------------------------------------------------------------------------

EE_MD = """## EE framing

- **HMM = discrete-state-space dynamical system.** Latent state evolves by a
  Markov chain; observations are emitted via a state-conditional distribution.
  The continuous-state analogue you already know is the Kalman filter.
- **Viterbi = MAP decoding via DP on the trellis.** Identical machinery to
  ML-decoding a convolutional code: max-product belief propagation along a
  chain. Implement in log-space because the path probability collapses
  exponentially with sequence length.
- **Forward = sum-product BP.** Replace `max` with `logsumexp` and you get
  the sequence likelihood; combined with backward you get per-position
  posterior marginals — and **the marginal-max sequence is not the
  Viterbi path** (this is the same gap as MAP-vs-MMSE in EE).
- **Profile HMM = position-specific matched-filter bank.** One emission
  distribution per match column makes the model a *sequence-shaped* filter
  bank: the score is exactly the sum of column-by-column log-odds — a
  log-likelihood-ratio detector — against a null background.
- **Why HMMs still win on sparse-data segmentation:** strong inductive bias
  (Markov state + position-specific emissions), interpretable parameters,
  zero-training-data flavour (specify by hand), and `O(T K^2)` exact
  inference that you can audit. Where deep models win: dense data and the
  budget to fine-tune."""


# ---------------------------------------------------------------------------
# Assemble notebook
# ---------------------------------------------------------------------------

def build():
    nb = new_notebook()
    cells = []

    cells.append(md(TITLE_MD))
    cells.append(md(AHA_MD))
    cells.append(code(PREAMBLE_CODE))
    cells.append(code(IMPORTS_CODE))

    # Step 1
    cells.append(md(STEP1_MD))
    cells.append(code(STEP1_TODO))
    cells.append(md(STEP1_SOL_MD))
    cells.append(hidden_solution(STEP1_SOL))

    # Step 2
    cells.append(md(STEP2_MD))
    cells.append(code(STEP2_TODO))
    cells.append(md("*Click ▶ to expand the reference solution.*"))
    cells.append(hidden_solution(STEP2_SOL))

    # Step 3
    cells.append(md(STEP3_MD))
    cells.append(code(STEP3_TODO))
    cells.append(md("*Click ▶ to expand the reference solution.*"))
    cells.append(hidden_solution(STEP3_SOL))

    # Step 4
    cells.append(md(STEP4_MD))
    cells.append(code(STEP4_TODO))
    cells.append(md("*Click ▶ to expand the reference solution.*"))
    cells.append(hidden_solution(STEP4_SOL))

    # Step 5
    cells.append(md(STEP5_MD))
    cells.append(code(STEP5_TODO))
    cells.append(md("*Click ▶ to expand the reference solution.*"))
    cells.append(hidden_solution(STEP5_SOL))

    # Self-check + closing
    cells.append(md(SELFCHECK_MD))
    cells.append(code(SELFCHECK_CODE))
    cells.append(md(EE_MD))

    nb.cells = cells
    nb.metadata = {
        "kernelspec": {
            "display_name": "Python 3",
            "language":     "python",
            "name":         "python3",
        },
        "language_info": {"name": "python"},
        "colab": {"provenance": [], "toc_visible": True},
    }

    apply_colab_form(nb)
    nbf.write(nb, OUT)
    return nb


if __name__ == "__main__":
    nb = build()
    n_md = sum(1 for c in nb.cells if c.cell_type == "markdown")
    n_code = sum(1 for c in nb.cells if c.cell_type == "code")
    print(f"Wrote {OUT}  ({len(nb.cells)} cells: {n_md} markdown, {n_code} code)")
