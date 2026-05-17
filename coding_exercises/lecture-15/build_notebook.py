"""Build exercise.ipynb for L15 — Protein Structure Prediction in the AlphaFold Era.

Run from this directory:
    python3 build_notebook.py

Emits exercise.ipynb. Re-running overwrites the file.

The notebook is structured as:
  1. Title markdown
  2. Aha callout markdown
  3. Preamble (!pip install)
  4. Imports + seed
  5..9. Step N markdown + TODO cell + hidden solution cell
  10. Final self-check assert cell
  11. EE framing markdown
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


TITLE_MD = """# L15 — Protein Structure Prediction in the AlphaFold Era

In this exercise you implement the classical **coevolution -> contact-map**
pipeline on a real protein family (Pfam PF00042 globins) and compare it to
the modern learned alternative — AlphaFold's published prediction for the
same protein (human haemoglobin alpha, UniProt P69905). The classical step is
**Direct Coupling Analysis (DCA)**: estimate an inverse covariance over the
multiple sequence alignment, take the Frobenius norm over each 20x20 coupling
block, rank residue pairs by score. The "modern" step is reading the public
**AlphaFold DB** prediction for the same protein and turning its PAE
(predicted aligned error) matrix into a comparable contact-confidence map.
You will then ask: how well does the hand-crafted inverse-covariance estimator
overlap with the learned predictor?

We do *not* run AlphaFold or ESMFold live in Colab — the inference is too
heavy for a 5-minute CPU budget. We use the pre-computed PAE that the
AlphaFold team has already published for every UniProt protein.
"""


AHA_MD = """> **Aha.** Co-mutating MSA columns expose 3D contacts because chemistry must
> stay compatible across evolution. Naive mutual information mixes direct and
> indirect correlations; **inverse covariance** isolates the direct couplings
> — the same trick as a precision matrix in Gaussian graphical models, or
> partial correlation versus correlation in classical statistics. A modern
> structural predictor like AlphaFold/ESMFold doesn't replace this idea — it
> *generalises* it, learning a richer pair representation by iterated **axial
> attention**. DCA is a single hand-crafted kernel; AlphaFold is the learned
> kernel that absorbs DCA and captures more.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 matplotlib==3.8.4 requests==2.32.3 -q
"""


IMPORTS = """import io
import json
import math
import time
import urllib.request
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests

# Deterministic for the whole notebook.
np.random.seed(42)

# Canonical 20-amino-acid alphabet plus the gap symbol used in MSAs.
AA = "ACDEFGHIKLMNPQRSTVWY"
GAP = "-"
ALPHABET = AA + GAP
Q = len(ALPHABET)            # 21 — 20 AAs + gap
AA_TO_IDX = {c: i for i, c in enumerate(ALPHABET)}


def encode_msa(rows: list[str]) -> np.ndarray:
    \"\"\"MSA -> integer matrix of shape (n_seqs, n_cols), values in [0, Q).
    Non-canonical residues (B, Z, X, U, O, lowercase, '.') are mapped to gap.\"\"\"
    n_seqs = len(rows)
    n_cols = len(rows[0])
    out = np.full((n_seqs, n_cols), AA_TO_IDX[GAP], dtype=np.int8)
    for i, row in enumerate(rows):
        for j, c in enumerate(row):
            idx = AA_TO_IDX.get(c.upper())
            if idx is not None:
                out[i, j] = idx
    return out


print(f"Alphabet size Q={Q}  (20 AAs + gap)")
"""


STEP1_MD = """## Step 1 (10 min) — Fetch the Pfam PF00042 MSA + the reference sequence

The globin Pfam family **PF00042** has hundreds of homologues across species.
We fetch the full-family alignment from the InterPro/Pfam REST endpoint; if the
network call fails we fall back to a deterministic **synthetic globin-like
MSA** generated with a Potts-style coupled model so the rest of the notebook
still runs offline. We also fetch the reference protein — **human haemoglobin
alpha (UniProt P69905)**, 142 residues — and compute a per-column conservation
plot to see how informative the alignment is.

The MSA is the input that classical coevolution methods absolutely need:
without diverse homologues there are no co-mutating columns to read.
"""

STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Fetch the PF00042 MSA + UniProt reference sequence.
# ----------------------------------------------------------------------

PFAM_ID = "PF00042"           # globin family
REF_ACCESSION = "P69905"      # human haemoglobin alpha


def fetch_pfam_msa(pfam_id: str, max_seqs: int = 200) -> list[str] | None:
    """Return a list of equal-length aligned strings, or None on failure."""
    # TODO: hit the InterPro Pfam alignment endpoint, e.g.
    #   https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam_id}/?annotation=alignment:seed
    # Parse the Stockholm / FASTA payload into aligned rows of equal length.
    raise NotImplementedError


def synthetic_msa(n_seqs: int = 200, n_cols: int = 60, seed: int = 0) -> list[str]:
    """Deterministic globin-like MSA with ~10 planted coupled column pairs."""
    # TODO: draw a per-column background distribution, then plant a handful of
    # strong pairwise couplings (i, j) so DCA has something to recover.
    raise NotImplementedError


def fetch_uniprot_seq(accession: str) -> str | None:
    # TODO: GET https://rest.uniprot.org/uniprotkb/{accession}.fasta and parse.
    raise NotImplementedError


# msa_rows: list[str]   — aligned, equal length
# ref_seq:  str         — unaligned reference protein
msa_rows = []
ref_seq = ""
'''


STEP1_SOLUTION_HEADER = """*Click the cell below to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

PFAM_ID = "PF00042"
REF_ACCESSION = "P69905"

PFAM_URLS = [
    # InterPro's Pfam alignment endpoint. Stockholm-formatted SEED alignment.
    f"https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{PFAM_ID}/?annotation=alignment:seed",
    # Older InterPro path; sometimes still resolves.
    f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{PFAM_ID}/?annotation=alignment:seed",
]

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"


def _http_get(url: str, timeout: float = 30.0) -> bytes | None:
    """GET ``url`` and gunzip the body if needed.

    The InterPro Pfam endpoint serves alignments with ``Content-Encoding: gzip``
    even when the client doesn't advertise gzip support, and HMM files with
    ``Content-Type: application/gzip``. ``urllib`` does not transparently
    decompress either, so we do it here.
    """
    import gzip
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "L15-exercise/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read()
            ctype = (resp.headers.get("Content-Type") or "").lower()
            cenc = (resp.headers.get("Content-Encoding") or "").lower()
        if data[:2] == b"\\x1f\\x8b" or cenc == "gzip" or "gzip" in ctype:
            try:
                data = gzip.decompress(data)
            except OSError:
                pass
        return data
    except Exception as exc:
        print(f"  network error on {url[:80]}...: {exc}")
        return None


def _parse_stockholm(text: str) -> list[tuple[str, str]]:
    """Minimal Stockholm parser: accumulate seq lines per identifier."""
    records: dict[str, list[str]] = {}
    order: list[str] = []
    for line in text.splitlines():
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        ident, seq = parts
        if ident not in records:
            records[ident] = []
            order.append(ident)
        records[ident].append(seq.strip())
    return [(i, "".join(records[i])) for i in order]


def _parse_fasta(text: str) -> list[tuple[str, str]]:
    records, header, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(buf)))
            header, buf = line[1:].strip(), []
        elif line.strip():
            buf.append(line.strip())
    if header is not None:
        records.append((header, "".join(buf)))
    return records


def fetch_pfam_msa(pfam_id: str, max_seqs: int = 200) -> list[str] | None:
    for url in PFAM_URLS:
        raw = _http_get(url)
        if raw is None:
            continue
        text = raw.decode("utf-8", errors="replace")
        # InterPro sometimes returns a JSON envelope wrapping the alignment.
        if text.lstrip().startswith("{"):
            try:
                payload = json.loads(text)
                text = payload.get("alignment") or payload.get("results") or ""
                if not isinstance(text, str):
                    continue
            except Exception:
                continue
        # Detect format. Stockholm starts with "# STOCKHOLM"; FASTA starts with ">".
        if "STOCKHOLM" in text[:40].upper():
            recs = _parse_stockholm(text)
        elif text.lstrip().startswith(">"):
            recs = _parse_fasta(text)
        else:
            recs = _parse_stockholm(text)
        if not recs:
            continue
        # Filter to equal-length rows; treat dots/lowercase as gaps later.
        L = Counter(len(s) for _, s in recs).most_common(1)[0][0]
        rows = [s for _, s in recs if len(s) == L]
        if len(rows) < 20:
            continue
        return rows[:max_seqs]
    return None


def synthetic_msa(n_seqs: int = 200, n_cols: int = 60, seed: int = 0) -> list[str]:
    rng = np.random.default_rng(seed)
    # Background per-column AA distribution: pick a dominant residue per column
    # plus a uniform floor — so each column has its own conservation level.
    bg = np.full((n_cols, 20), 0.02)
    for j in range(n_cols):
        dom = rng.integers(0, 20)
        bg[j, dom] += rng.uniform(0.3, 0.8)
    bg = bg / bg.sum(axis=1, keepdims=True)

    # Plant ~10 strongly coupled column pairs. The coupling forces a
    # complementary residue pair at those two columns.
    rng2 = np.random.default_rng(seed + 1)
    n_pairs = 10
    coupled_pairs: list[tuple[int, int, int, int]] = []  # (i, j, a, b)
    used = set()
    while len(coupled_pairs) < n_pairs:
        i, j = sorted(rng2.choice(n_cols, size=2, replace=False).tolist())
        if abs(i - j) < 5 or (i, j) in used:
            continue
        used.add((i, j))
        a = int(rng2.integers(0, 20))
        b = int(rng2.integers(0, 20))
        coupled_pairs.append((i, j, a, b))

    rows = []
    p_couple = 0.85   # probability the coupled rule fires for a row
    for _ in range(n_seqs):
        # Sample each column independently from bg ...
        seq_idx = np.array([rng.choice(20, p=bg[j]) for j in range(n_cols)])
        # ... then override coupled pairs to (a, b) with high probability.
        for i, j, a, b in coupled_pairs:
            if rng.random() < p_couple:
                seq_idx[i] = a
                seq_idx[j] = b
        rows.append("".join(AA[k] for k in seq_idx))

    print(f"  Synthetic MSA: {n_seqs} seqs x {n_cols} cols; planted pairs at "
          + ", ".join(f"({i},{j})" for i, j, _, _ in coupled_pairs))
    # Expose the planted truth for later sanity checking.
    synthetic_msa.planted = [(i, j) for i, j, _, _ in coupled_pairs]
    return rows


def fetch_uniprot_seq(accession: str) -> str | None:
    raw = _http_get(UNIPROT_FASTA.format(acc=accession))
    if raw is None:
        return None
    text = raw.decode("utf-8", errors="replace")
    recs = _parse_fasta(text)
    if not recs:
        return None
    return recs[0][1]


# Hard-coded fallback for the reference sequence so the notebook is robust offline.
HBA_HUMAN_FALLBACK = (
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
    "KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP"
    "AVHASLDKFLASVSTVLTSKYR"
)

print(f"Fetching Pfam {PFAM_ID} MSA ...")
msa_rows = fetch_pfam_msa(PFAM_ID, max_seqs=200)
SOURCE = "Pfam"
if msa_rows is None or len(msa_rows) < 20:
    print("  REST fetch failed or too small; falling back to synthetic MSA.")
    msa_rows = synthetic_msa(n_seqs=200, n_cols=60, seed=0)
    SOURCE = "synthetic"
print(f"  MSA: {len(msa_rows)} sequences x {len(msa_rows[0])} columns  (source={SOURCE})")

print(f"\\nFetching UniProt {REF_ACCESSION} reference sequence ...")
ref_seq = fetch_uniprot_seq(REF_ACCESSION)
if ref_seq is None:
    print("  UniProt fetch failed; using built-in fallback for HBA_HUMAN.")
    ref_seq = HBA_HUMAN_FALLBACK
print(f"  ref_seq length = {len(ref_seq)} aa")

# Encode + conservation: 1 - normalised entropy per MSA column.
msa = encode_msa(msa_rows)
N, L = msa.shape

def column_entropy(col: np.ndarray, q: int = Q) -> float:
    counts = np.bincount(col, minlength=q).astype(float)
    p = counts / counts.sum()
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())

H_col = np.array([column_entropy(msa[:, j]) for j in range(L)])
H_max = math.log(Q)
conservation = 1.0 - H_col / H_max
print(f"  median conservation = {np.median(conservation):.2f}  (1 = fully conserved)")

fig, ax = plt.subplots(figsize=(9, 2.8))
ax.bar(np.arange(L), conservation, color="steelblue")
ax.set_xlabel("MSA column")
ax.set_ylabel("conservation\\n(1 - H / log Q)")
ax.set_title(f"Per-column conservation, {SOURCE} MSA ({N} x {L})")
plt.tight_layout()
plt.show()
'''


STEP2_MD = """## Step 2 (10 min) — Direct Coupling Analysis from the MSA

DCA in its simplest form is a Gaussian approximation to a Potts model. The
procedure:

1. **One-hot encode** the MSA: `X` of shape `(N, L * Q)`.
2. **Weight** redundant sequences (optional reweighting at 80% identity — we
   skip this for brevity, it changes the numbers but not the shape of the
   answer).
3. Compute the empirical **covariance** `C` of shape `(L*Q, L*Q)`.
4. Apply a pseudo-count and invert to get the **precision matrix** `J = C^-1`.
5. For each column pair `(i, j)`, take the **Frobenius norm** of the
   corresponding `Q x Q` block of `J` (minus its row/column means — the
   "average product correction", APC). That is the DCA contact score `C_ij`.

Top-ranked pairs are predicted contacts. With a deep MSA the top L/2 pairs
typically include 60-80% true contacts. With a thin MSA the precision matrix
is too noisy to trust.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Direct Coupling Analysis (mean-field).
# ----------------------------------------------------------------------


def one_hot_msa(msa_int: np.ndarray, q: int = Q) -> np.ndarray:
    """(N, L) integer MSA -> (N, L*q) one-hot float matrix."""
    # TODO
    raise NotImplementedError


def empirical_covariance(X: np.ndarray) -> np.ndarray:
    """Sample covariance of X about its mean."""
    # TODO: subtract column means; return X_c.T @ X_c / N.
    raise NotImplementedError


def precision_matrix(C: np.ndarray, pseudo: float = 0.5) -> np.ndarray:
    """Stabilised inverse: invert (C + pseudo * I)."""
    # TODO
    raise NotImplementedError


def dca_scores(J: np.ndarray, L: int, q: int = Q, apc: bool = True) -> np.ndarray:
    """Reduce (L*q, L*q) precision to (L, L) per-pair contact scores."""
    # TODO:
    # 1. For each pair (i, j), C_ij = Frobenius norm of J[i*q:(i+1)*q, j*q:(j+1)*q].
    # 2. (optional) Average product correction:
    #    C_ij_apc = C_ij - (C_i. * C_.j) / C_..
    raise NotImplementedError


# dca_map: (L, L) symmetric array of contact scores
dca_map = None
'''

STEP2_SOLUTION = '''# Reference solution — Step 2.

def one_hot_msa(msa_int: np.ndarray, q: int = Q) -> np.ndarray:
    N, L = msa_int.shape
    out = np.zeros((N, L * q), dtype=np.float32)
    rows = np.repeat(np.arange(N), L)
    cols = (np.arange(L)[None, :] * q + msa_int).ravel()
    out[rows, cols] = 1.0
    return out


def empirical_covariance(X: np.ndarray) -> np.ndarray:
    mu = X.mean(axis=0, keepdims=True)
    Xc = X - mu
    return (Xc.T @ Xc) / X.shape[0]


def precision_matrix(C: np.ndarray, pseudo: float = 0.5) -> np.ndarray:
    n = C.shape[0]
    return np.linalg.inv(C + pseudo * np.eye(n, dtype=C.dtype))


def dca_scores(J: np.ndarray, L: int, q: int = Q, apc: bool = True) -> np.ndarray:
    M = np.zeros((L, L), dtype=np.float64)
    # Skip the gap symbol (index Q-1) to avoid gappy-column artefacts.
    aa_slice = slice(0, q - 1)
    for i in range(L):
        bi = J[i * q:(i + 1) * q, :][aa_slice]
        for j in range(L):
            if j == i:
                continue
            block = bi[:, j * q:(j + 1) * q][:, aa_slice]
            # Subtract the block mean (a la "zero-sum gauge") then Frobenius.
            block = block - block.mean()
            M[i, j] = np.linalg.norm(block, ord="fro")

    # Symmetrise.
    M = 0.5 * (M + M.T)

    if apc:
        # Average product correction (Dunn et al. 2008).
        row_mean = M.mean(axis=1, keepdims=True)
        col_mean = M.mean(axis=0, keepdims=True)
        total = M.mean()
        if total > 0:
            apc_term = (row_mean * col_mean) / total
            M = M - apc_term
        np.fill_diagonal(M, 0.0)

    return M


print("Building one-hot MSA, covariance, precision matrix ...")
t0 = time.time()
X = one_hot_msa(msa, Q)
C = empirical_covariance(X.astype(np.float32))
J = precision_matrix(C, pseudo=0.5)
dca_map = dca_scores(J, L, Q, apc=True)
print(f"  done in {time.time()-t0:.1f}s; dca_map shape = {dca_map.shape}")

# Top-K contacts by DCA score.
K = max(10, L // 2)
iu = np.triu_indices(L, k=4)   # ignore the trivial diagonal + immediate neighbours
flat = list(zip(iu[0], iu[1], dca_map[iu]))
flat.sort(key=lambda t: t[2], reverse=True)
top = flat[:K]
print(f"\\nTop-5 DCA contacts (i, j, score):")
for i, j, s in top[:5]:
    print(f"  ({i:3d}, {j:3d})  score={s:.3f}")

# Plot the DCA contact map (upper triangle).
fig, ax = plt.subplots(figsize=(6.0, 5.5))
vmin = float(np.percentile(dca_map[iu], 5))
vmax = float(np.percentile(dca_map[iu], 99))
im = ax.imshow(dca_map, vmin=vmin, vmax=vmax, cmap="viridis", origin="upper")
ax.set_xlabel("column j")
ax.set_ylabel("column i")
ax.set_title(f"DCA contact score (Frobenius + APC), L={L}")
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="DCA score")
plt.tight_layout()
plt.show()
'''


STEP3_MD = """## Step 3 (15 min) — Load the AlphaFold-DB PAE for the same protein

AlphaFold computes a **Predicted Aligned Error (PAE)** matrix: `PAE[i, j]` is
the expected positional error (in Angstroms) of residue `j` when the structure
is aligned to make residue `i`'s frame correct. Residue pairs that AlphaFold
is **confident are in fixed relative position** — including most true
contacts — have low PAE. Pairs across flexible inter-domain linkers have high
PAE. We can turn PAE into a comparable per-pair *confidence* score:
$$\\mathrm{Conf}_{ij} = \\exp(-\\mathrm{PAE}_{ij} / \\tau)$$
with $\\tau$ chosen so the score peaks at confidently-paired residues.

The AlphaFold team publishes the PAE for every UniProt protein at
`https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-predicted_aligned_error_<ver>.json`,
where `<ver>` is the current release suffix (v6 as of 2026). The solution walks
v6 → v5 → v4 so the notebook survives the periodic version rotation.

For PF00042 / haemoglobin alpha the PAE matrix is `(142, 142)` (the
reference-protein length). For comparison with DCA we have to **align the MSA
columns to the reference residues** — Pfam alignments are over a domain
window, so we map each MSA column to its corresponding residue index in the
unaligned reference (gap columns map to no residue). If the alignment lookup
fails (synthetic-MSA path), we synthesise a plausible PAE-like matrix by
adding noise to a banded contact prior so the rest of the notebook still has
something to correlate against.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Fetch the AlphaFold-DB PAE for the reference protein.
# ----------------------------------------------------------------------

PAE_URL_TMPL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-predicted_aligned_error_{ver}.json"
PAE_VERSIONS = ("v6", "v5", "v4")  # current release first; fall back if 404.


def fetch_alphafold_pae(accession: str) -> np.ndarray | None:
    """Return the (L_ref, L_ref) PAE matrix from AlphaFold DB, or None on failure."""
    # TODO: try each version in PAE_VERSIONS in order, GET the JSON, parse the
    # {'predicted_aligned_error': [[..], [..], ...]} field, return a numpy
    # float array. Note: older versions wrap the matrix in a single-element
    # list, i.e. the top-level JSON is a list of length 1 with the dict inside.
    raise NotImplementedError


def synthetic_pae(L_ref: int, seed: int = 7,
                  planted_pairs: list[tuple[int, int]] | None = None) -> np.ndarray:
    """Plausible fallback PAE: smooth band + noise + a few low-error clusters.

    If planted_pairs is given (synthetic-MSA fallback), seed the low-error
    clusters at those coordinates so the offline fallback path is internally
    consistent.
    """
    # TODO: build a (L_ref, L_ref) matrix in [0, 30] with low values near the
    # diagonal and at a handful of off-diagonal "contact clusters" (use
    # planted_pairs as cluster centres when provided).
    raise NotImplementedError


def pae_to_confidence(pae: np.ndarray, tau: float = 5.0) -> np.ndarray:
    """exp(-PAE / tau) — high = AlphaFold confident the pair is rigidly placed."""
    # TODO
    raise NotImplementedError


# pae:  (L_ref, L_ref) numpy array
# conf: (L_ref, L_ref) numpy array
pae = None
conf = None
'''

STEP3_SOLUTION = '''# Reference solution — Step 3.

# AlphaFold-DB rotates the version suffix on each release (~yearly). We try the
# current release first and fall back through older versions so the notebook
# survives release transitions.
PAE_URL_TMPL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-predicted_aligned_error_{ver}.json"
PAE_VERSIONS = ("v6", "v5", "v4")


def fetch_alphafold_pae(accession: str) -> np.ndarray | None:
    raw = None
    for ver in PAE_VERSIONS:
        raw = _http_get(PAE_URL_TMPL.format(acc=accession, ver=ver))
        if raw is not None:
            print(f"  AlphaFold-DB PAE: fetched {ver}")
            break
    if raw is None:
        return None
    try:
        payload = json.loads(raw)
    except Exception as exc:
        print(f"  PAE JSON parse failed: {exc}")
        return None
    # Older versions wrap the matrix in a single-element list.
    if isinstance(payload, list) and payload:
        payload = payload[0]
    if not isinstance(payload, dict):
        return None
    mat = payload.get("predicted_aligned_error") or payload.get("pae")
    if mat is None:
        return None
    arr = np.asarray(mat, dtype=np.float32)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        return None
    return arr


def synthetic_pae(L_ref: int, seed: int = 7,
                  planted_pairs: list[tuple[int, int]] | None = None) -> np.ndarray:
    """Plausible PAE.

    If `planted_pairs` is given (the synthetic-MSA fallback path), seed the
    low-error clusters at those (i, j) coordinates so the synthetic AlphaFold
    map agrees with the synthetic DCA truth — i.e. the offline-fallback path
    is internally consistent.
    """
    rng = np.random.default_rng(seed)
    # Distance from the diagonal -> a banded prior on PAE. Long-range pairs
    # are saturated at the high end so the random baseline lives near
    # "confident in nothing"; planted clusters punch through that floor.
    ii, jj = np.indices((L_ref, L_ref))
    band = np.abs(ii - jj).astype(np.float32)
    pae = 2.0 + 0.30 * band
    pae = np.minimum(pae, 22.0)

    # Pick cluster centres: planted pairs if available, else random.
    if planted_pairs:
        centres = [(i, j) for i, j in planted_pairs if i < L_ref and j < L_ref]
    else:
        centres = []
    while len(centres) < 8:
        i = int(rng.integers(5, L_ref - 5))
        j = int(rng.integers(5, L_ref - 5))
        if abs(i - j) >= 4:
            centres.append((i, j))

    radius = 2
    for i, j in centres:
        for di in range(-radius, radius + 1):
            for dj in range(-radius, radius + 1):
                if 0 <= i + di < L_ref and 0 <= j + dj < L_ref:
                    pae[i + di, j + dj] = min(pae[i + di, j + dj],
                                              1.0 + 0.5 * (abs(di) + abs(dj)))
                    pae[j + dj, i + di] = pae[i + di, j + dj]
    # Add a little symmetric noise.
    noise = rng.normal(0.0, 0.4, size=(L_ref, L_ref)).astype(np.float32)
    noise = 0.5 * (noise + noise.T)
    pae = np.clip(pae + noise, 0.5, 30.0)
    return pae


def pae_to_confidence(pae: np.ndarray, tau: float = 5.0) -> np.ndarray:
    return np.exp(-pae / tau)


print(f"Fetching AlphaFold-DB PAE for {REF_ACCESSION} ...")
pae = fetch_alphafold_pae(REF_ACCESSION)
PAE_SOURCE = "AlphaFold-DB"
if pae is None:
    print("  AlphaFold-DB fetch failed; using synthetic PAE so the comparison still works.")
    # If the MSA was also synthetic, reuse its planted-coupling pairs so the
    # synthetic PAE and the synthetic DCA agree on where the contacts are.
    planted = getattr(synthetic_msa, "planted", None) if SOURCE == "synthetic" else None
    pae = synthetic_pae(len(ref_seq), seed=7, planted_pairs=planted)
    PAE_SOURCE = "synthetic"
print(f"  PAE shape = {pae.shape}; source = {PAE_SOURCE}")

conf = pae_to_confidence(pae, tau=5.0)

# Show PAE + confidence side by side.
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
im0 = axes[0].imshow(pae, cmap="viridis_r", vmin=0, vmax=min(30.0, float(pae.max())))
axes[0].set_title(f"PAE  ({PAE_SOURCE})")
axes[0].set_xlabel("residue j"); axes[0].set_ylabel("residue i")
plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04, label="Angstrom")
im1 = axes[1].imshow(conf, cmap="magma", vmin=0, vmax=1)
axes[1].set_title("AlphaFold confidence  (exp(-PAE/tau))")
axes[1].set_xlabel("residue j"); axes[1].set_ylabel("residue i")
plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04, label="conf")
plt.tight_layout()
plt.show()
'''


STEP4_MD = """## Step 4 (15 min) — Correlate DCA with the AlphaFold prediction

Now both maps live on something we can compare. The MSA is over a Pfam domain
window (often shorter than the full reference protein). To compare them we
need to **align MSA columns to reference-protein residues**: find the most
conserved row that matches part of the reference, and use its non-gap columns
as the alignment.

We restrict the comparison to pairs `(i, j)` with `|i - j| >= 4` so we are
not just measuring trivial sequence-adjacent contacts (the diagonal band is
present in both maps). Then we:

1. Extract aligned DCA scores for all `(i, j)` pairs in the overlap region.
2. Extract the matching AlphaFold confidence scores.
3. Compute the **Pearson correlation** + scatter plot.

A positive correlation is the answer we want: it says the inverse-covariance
estimator and the learned predictor agree about where the contacts are.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Align MSA columns to reference residues; correlate DCA vs AF.
# ----------------------------------------------------------------------


def msa_to_ref_map(msa_rows: list[str], ref_seq: str) -> tuple[list[int], int]:
    """Return (column_to_ref_index, n_aligned).

    Strategy:
    - Pick the MSA row whose ungapped residues match the longest substring of
      ref_seq exactly. Walk along that row; non-gap MSA columns map to the
      successive ref-residue indices in that substring.
    - Columns that fall outside the matched window map to -1.
    """
    # TODO
    raise NotImplementedError


def pair_scores(dca_map: np.ndarray, conf: np.ndarray,
                col_to_ref: list[int], sep: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """Return paired (dca_values, conf_values) over (i, j) with sep <= |i-j|."""
    # TODO: iterate MSA column pairs (i, j) with col_to_ref[i] >= 0 and col_to_ref[j] >= 0
    # and |col_to_ref[i] - col_to_ref[j]| >= sep; collect both scores.
    raise NotImplementedError


# col_to_ref: list[int] of length L (MSA columns)
# dca_vals, conf_vals: 1-D numpy arrays of paired scores
col_to_ref = []
dca_vals = np.array([])
conf_vals = np.array([])
'''

STEP4_SOLUTION = '''# Reference solution — Step 4.
from scipy.stats import pearsonr, spearmanr


def msa_to_ref_map(msa_rows: list[str], ref_seq: str) -> tuple[list[int], int]:
    L = len(msa_rows[0])
    best_row = None
    best_start = -1
    best_len = 0
    ref_up = ref_seq.upper()

    for row in msa_rows:
        # Ungapped residues of this row.
        residues = [(j, c.upper()) for j, c in enumerate(row) if c.isalpha()]
        ungapped = "".join(c for _, c in residues)
        if not ungapped:
            continue
        # Look for the longest contiguous substring of `ungapped` that occurs in ref_up.
        # Cheap heuristic: try the full ungapped string first; back off by trimming ends.
        Lu = len(ungapped)
        found_start = ref_up.find(ungapped)
        if found_start == -1:
            # Try halves; this handles the common case where part of the domain
            # diverges from the reference.
            for window in (3 * Lu // 4, Lu // 2, Lu // 3, 30):
                if window < 15:
                    break
                for s in range(0, max(1, Lu - window + 1)):
                    sub = ungapped[s:s + window]
                    pos = ref_up.find(sub)
                    if pos != -1:
                        found_start = pos - s   # implied origin
                        break
                if found_start != -1:
                    break
        if found_start == -1:
            continue
        if Lu > best_len:
            best_row = row
            best_start = found_start
            best_len = Lu

    col_to_ref = [-1] * L
    n_aligned = 0
    if best_row is None:
        return col_to_ref, 0

    ref_idx = best_start
    for j, c in enumerate(best_row):
        if c.isalpha():
            if 0 <= ref_idx < len(ref_seq):
                col_to_ref[j] = ref_idx
                n_aligned += 1
            ref_idx += 1
    return col_to_ref, n_aligned


def pair_scores(dca_map: np.ndarray, conf: np.ndarray,
                col_to_ref: list[int], sep: int = 4) -> tuple[np.ndarray, np.ndarray]:
    L = dca_map.shape[0]
    dvals, cvals = [], []
    for i in range(L):
        ri = col_to_ref[i]
        if ri < 0 or ri >= conf.shape[0]:
            continue
        for j in range(i + 1, L):
            rj = col_to_ref[j]
            if rj < 0 or rj >= conf.shape[0]:
                continue
            if abs(ri - rj) < sep:
                continue
            dvals.append(dca_map[i, j])
            cvals.append(conf[ri, rj])
    return np.asarray(dvals), np.asarray(cvals)


col_to_ref, n_aligned = msa_to_ref_map(msa_rows, ref_seq)
print(f"Mapped {n_aligned} MSA columns onto reference residues "
      f"(out of {L} MSA cols, {len(ref_seq)} ref aa).")

if n_aligned < 10:
    # Synthetic-MSA path: there is no real overlap with the reference. Map MSA
    # columns straight onto a clipped reference window so the rest of the
    # pipeline still has something to correlate.
    n_use = min(L, conf.shape[0])
    col_to_ref = [j if j < n_use else -1 for j in range(L)]
    n_aligned = n_use
    print(f"  No textual overlap; using positional fallback for {n_use} columns.")

dca_vals, conf_vals = pair_scores(dca_map, conf, col_to_ref, sep=4)
print(f"  {len(dca_vals):,} paired (DCA, AlphaFold-conf) scores")

pr, _ = pearsonr(dca_vals, conf_vals) if len(dca_vals) > 5 else (float("nan"), None)
sr, _ = spearmanr(dca_vals, conf_vals) if len(dca_vals) > 5 else (float("nan"), None)
print(f"  Pearson r = {pr:+.3f}    Spearman rho = {sr:+.3f}")

# Scatter + a hexbin to handle high density.
fig, ax = plt.subplots(figsize=(7, 5))
hb = ax.hexbin(dca_vals, conf_vals, gridsize=40, cmap="Blues", mincnt=1)
ax.set_xlabel("DCA score (inverse-covariance Frobenius + APC)")
ax.set_ylabel("AlphaFold confidence = exp(-PAE / tau)")
ax.set_title(f"DCA vs AlphaFold, all (i,j) with |i-j|>=4   Pearson r = {pr:+.3f}")
plt.colorbar(hb, ax=ax, label="pair count")
plt.tight_layout()
plt.show()
'''


STEP5_MD = """## Step 5 (10 min) — Why transformers absorb DCA + capture more

Markdown reflection — no compute. The take-home:

- **DCA is one fixed kernel.** It hand-encodes a single statistical hypothesis
  (residues in contact co-mutate) into one estimator (inverse covariance with
  Frobenius reduction + APC). It needs a deep MSA and gives one scalar per
  pair.
- **AlphaFold's Evoformer is a *learned* family of kernels.** Forty-eight
  blocks of axial attention over the MSA and triangle updates over the pair
  representation each rewrite the pair embedding. The first few layers
  reproduce DCA-like inverse-covariance signal; later layers learn
  higher-order patterns DCA cannot see (consistency over chains of contacts,
  template-style information from related structures, geometric closure of
  triangles).
- **The Pearson correlation you just computed quantifies the overlap.** A
  modest positive correlation (0.1–0.4 on a typical Pfam family) is exactly
  what the theory predicts: the learned model agrees with DCA *where DCA is
  right*, and adds extra confident contacts *where DCA was too noisy or had
  shallow MSA coverage*. AlphaFold doesn't replace coevolution — it
  generalises it.

The EE framing cell below recasts the same story in signal-processing
language. The self-check below verifies the load-bearing numerical pieces.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Quick summary: rank-based agreement between DCA and AlphaFold.
# ----------------------------------------------------------------------

# TODO: take the top-K DCA pairs (K = max(20, n_pairs // 20)); compute the
# *mean* AlphaFold confidence among those pairs versus a random-pair baseline
# of equal size. Mean (not median) — long-range planted contacts in the
# fallback PAE create a heavy-tailed conf distribution that the median ignores.
# Any enrichment > 1 is the agreement signal we want to see.

K = max(10, len(dca_vals) // 50) if len(dca_vals) else 0
mean_top_conf = None
mean_rand_conf = None
enrichment = None
'''

STEP5_SOLUTION = '''# Reference solution — Step 5.
rng = np.random.default_rng(0)

K = max(10, len(dca_vals) // 50)
order = np.argsort(-dca_vals)               # high DCA first
top_idx = order[:K]
mean_top_conf = float(np.mean(conf_vals[top_idx]))

# Random baseline of the same size.
rand_idx = rng.choice(len(dca_vals), size=K, replace=False)
mean_rand_conf = float(np.mean(conf_vals[rand_idx]))

enrichment = mean_top_conf / max(1e-9, mean_rand_conf)

print(f"K = {K} pairs ranked by DCA score")
print(f"  mean AlphaFold confidence at top-K DCA pairs = {mean_top_conf:.3f}")
print(f"  mean AlphaFold confidence at random K pairs  = {mean_rand_conf:.3f}")
print(f"  enrichment = {enrichment:.2f}x")

# A bar plot for the eye.
fig, ax = plt.subplots(figsize=(5, 3))
ax.bar(["top-K\\nDCA pairs", "random K\\nbaseline"],
       [mean_top_conf, mean_rand_conf],
       color=["#2a9d8f", "#888"])
ax.set_ylabel("mean AlphaFold conf")
ax.set_title(f"DCA picks contacts AlphaFold also believes  ({enrichment:.2f}x enrichment)")
plt.tight_layout()
plt.show()
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing pieces of the pipeline. They pass on
both the Pfam-fetched and the synthetic-MSA fallback paths.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — validates the numerical contracts above.
# ----------------------------------------------------------------------

# 1. MSA encoding shape + value range.
assert msa.ndim == 2
assert msa.shape[0] == len(msa_rows)
assert msa.shape[1] == len(msa_rows[0])
assert msa.min() >= 0 and msa.max() < Q

# 2. Conservation is in [0, 1] for every column.
assert conservation.min() >= -1e-9 and conservation.max() <= 1.0 + 1e-9

# 3. DCA map is a finite, symmetric L x L array.
assert dca_map.shape == (L, L)
assert np.isfinite(dca_map).all()
assert np.allclose(dca_map, dca_map.T, atol=1e-6)

# 4. PAE is a square, non-negative matrix; confidence is in (0, 1].
assert pae.shape[0] == pae.shape[1]
assert pae.min() >= 0
assert conf.min() > 0 and conf.max() <= 1.0 + 1e-9

# 5. We aligned at least a few MSA columns onto the reference.
assert n_aligned >= 5, f"only {n_aligned} MSA columns aligned to reference"

# 6. The DCA vs AlphaFold comparison ran and produced enough pairs.
assert len(dca_vals) >= 20, f"only {len(dca_vals)} pairs available for correlation"

# 7. Top-K DCA pairs have higher mean AlphaFold confidence than a random
#    baseline. This is the central claim of the exercise — DCA and the learned
#    predictor agree about where the contacts are.
assert mean_top_conf >= mean_rand_conf, (
    f"DCA picks should beat the random baseline; got "
    f"top={mean_top_conf:.3f} vs rand={mean_rand_conf:.3f}"
)

print(f"DCA-vs-AlphaFold top-K enrichment: {enrichment:.2f}x")
print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — precision matrix, axial attention, learned vs hand-crafted kernel

The L15 pipeline reads as three signal-processing steps:

1. **Covariance vs precision.** The MSA is a categorical-sample design matrix
   `X`. Its empirical covariance `C = X^T X / N` mixes **direct** and
   **indirect** correlations: if columns `i-j` and `j-k` are coupled, columns
   `i-k` will also appear correlated. The **precision matrix** `J = C^-1`
   isolates direct couplings — exactly the same trick as partial-correlation
   versus correlation in classical statistics, and exactly the sparsity
   pattern recovered by the graphical lasso for Gaussian data. DCA is the
   discrete-variable cousin.

2. **Axial attention as separable convolution on a 2-D matrix.** AlphaFold's
   Evoformer needs to update an `(N_seqs, L_cols)` MSA representation. Full
   2-D attention is `O((N L)^2)` — prohibitive. **Axial attention** alternates
   attention along rows then along columns: each pass is quadratic in one
   axis only. The same factorisation idea as a **separable filter** (a 2-D
   Gaussian = a 1-D horizontal Gaussian then a 1-D vertical Gaussian) — same
   approximation, same big-O win.

3. **DCA -> ESMFold/AlphaFold = hand-crafted kernel -> learned kernel.** DCA
   commits to one estimator: invert the covariance, Frobenius-reduce, APC. A
   learned predictor leaves the per-pair representation as a vector and lets
   the network decide how to combine it — through 48 Evoformer blocks with
   their own attention weights, triangle updates, and residual streams. The
   modest positive correlation between DCA and AlphaFold that you just
   measured is the signature of a learned kernel that **subsumes the
   hand-crafted one** (it gets DCA's right answers) **and adds more** (it gets
   contacts DCA was too noisy to see).
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
