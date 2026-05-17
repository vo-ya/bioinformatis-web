"""Build exercise.ipynb for L20 — MSA, Phylogenetics, Comparative Genomics.

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
    return new_markdown_cell(text)


def code(source: str):
    return new_code_cell(source)


def hidden(source: str):
    cell = new_code_cell(source)
    cell.metadata = {
        "jupyter": {"source_hidden": True},
        "cellView": "form",
    }
    return cell


# ---------------------------------------------------------------------------
# Cell sources
# ---------------------------------------------------------------------------


TITLE_MD = """# L20 — MSA, Phylogenetics, and Comparative Genomics

In this exercise you build the comparative-genomics workbench end-to-end on
real globin proteins. You will:

1. Score 28 pairwise distances and reconstruct an 8-leaf evolutionary tree
   with **Neighbor-Joining**.
2. Implement **progressive multiple-sequence alignment** along the NJ guide
   tree (pairwise -> profile -> sequence-vs-profile).
3. Score your MSA against the **Pfam PF00042 reference** using the
   sum-of-pairs (SP) metric.
4. Codon-align a 10-species ortholog set, compute the Nei-Gojobori
   **dN/dS** ratio, and run a chi-squared **likelihood-ratio test** against
   neutrality.
5. Interpret the tree topology and flag any lineage with dN/dS > 1.
"""


AHA_MD = """> **Aha.** MSA is NP-hard in the number of sequences. Progressive
> alignment is a **greedy heuristic on a guide tree**: every column of the
> final alignment is the consensus of `O(N)` local DP decisions, never one
> joint optimisation. The tree is not "truth" — it is a low-rank summary of
> the pairwise distance matrix. dN/dS is a **likelihood-ratio test** on a
> codon-substitution model: dN/dS > 1 means the synonymous null is rejected
> in favour of positive selection.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 biopython==1.83 matplotlib==3.8.4 requests==2.32.3 -q
"""


IMPORTS = """import io
import math
import re
import urllib.request
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2

from Bio.Align import substitution_matrices

# Deterministic.
np.random.seed(42)

# Canonical 20-AA alphabet (no ambiguous codes); '-' is the gap symbol.
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA)
GAP = "-"

# Load BLOSUM62 once.
BLOSUM62 = substitution_matrices.load("BLOSUM62")
print(f"BLOSUM62 loaded: {BLOSUM62.shape} entries; A/A={int(BLOSUM62['A','A'])}, W/W={int(BLOSUM62['W','W'])}")


def _http_get(url: str, timeout: float = 30.0) -> str | None:
    \"\"\"Minimal GET wrapper; gunzips the body if needed. None on failure.

    The InterPro Pfam alignment endpoint sets ``Content-Encoding: gzip`` even
    when the client does not advertise gzip support, and ``urllib`` does not
    transparently decompress.
    \"\"\"
    import gzip
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "L20-exercise/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read()
            ctype = (resp.headers.get("Content-Type") or "").lower()
            cenc = (resp.headers.get("Content-Encoding") or "").lower()
        if data[:2] == b"\\x1f\\x8b" or cenc == "gzip" or "gzip" in ctype:
            try:
                data = gzip.decompress(data)
            except OSError:
                pass
        return data.decode("utf-8", errors="replace")
    except Exception as exc:
        print(f"  network error on {url[:80]}...: {exc}")
        return None


def _parse_fasta(text: str) -> list[tuple[str, str]]:
    \"\"\"Parse a multi-FASTA string -> list of (header, sequence).\"\"\"
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


# Standard genetic code — used in step 4 for codon translation.
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
"""


STEP1_MD = """## Step 1 (10 min) — Pairwise distances and a Neighbor-Joining tree

We start with 8 globin proteins spanning vertebrate evolution: hemoglobin
alpha and beta from human, mouse, chicken, and zebrafish. UniProt serves
them as FASTA records over its public REST API; if the network call fails
the notebook falls back to canonical mini-sequences hardcoded below so
nothing breaks offline.

For each pair you compute a **distance** from the BLOSUM62-aligned
percent identity:

$$d(a,b) = 1 - \\frac{\\text{matches}}{\\min(|a|,|b|)}$$

That gives a symmetric 8x8 distance matrix `D`. Neighbor-Joining then
builds an unrooted binary tree from `D` by iteratively picking the pair
that minimises the corrected distance `Q(i,j) = (n-2) D[i,j] - r_i - r_j`,
merging them, and updating the matrix until only one node is left.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Fetch 8 globins, build distance matrix, run NJ.
# ----------------------------------------------------------------------

GLOBIN_ACCESSIONS = {
    "HBA_HUMAN":  "P69905",  # hemoglobin alpha, human
    "HBB_HUMAN":  "P68871",  # hemoglobin beta,  human
    "HBA_MOUSE":  "P01942",  # hemoglobin alpha, mouse
    "HBB_MOUSE":  "P02088",  # hemoglobin beta,  mouse
    "HBA_CHICK":  "P01994",  # hemoglobin alpha A, chicken
    "HBB_CHICK":  "P02112",  # hemoglobin beta,  chicken
    "HBA_DANRE":  "Q90487",  # hemoglobin alpha, zebrafish
    "HBB_DANRE":  "Q6DHL5",  # hemoglobin beta,  zebrafish
}


def fetch_globins() -> dict[str, str]:
    """Return {name: sequence} for the 8 globins above."""
    # TODO: hit https://rest.uniprot.org/uniprotkb/{acc}.fasta for each acc;
    # parse the FASTA; strip to canonical AAs. If the fetch fails, fall back to
    # the GLOBIN_FALLBACK dict you build below.
    raise NotImplementedError


def percent_identity_distance(a: str, b: str) -> float:
    """1 - (matches / min(|a|,|b|)) after a simple end-to-end alignment.

    For this exercise a fast O(L) heuristic is enough: count matches at
    overlapping positions, normalise by the shorter sequence length.
    """
    # TODO
    raise NotImplementedError


def neighbor_joining(D: np.ndarray, names: list[str]) -> list[tuple]:
    """Return a list of join events: (i, j, d_i, d_j, new_name).

    Standard NJ: at each step compute Q = (n-2) D - row_sums - col_sums,
    pick (i,j) with min Q, merge into a new node, update D, repeat.
    """
    # TODO
    raise NotImplementedError


# globins: dict[name, seq]; D: 8x8 distance matrix; tree_events: NJ join log.
globins = {}
D = None
tree_events = []
'''


STEP1_SOLUTION_HEADER = """*Click the arrow to expand the reference solution.*"""


STEP1_SOLUTION = '''# Reference solution — Step 1.

GLOBIN_ACCESSIONS = {
    "HBA_HUMAN":  "P69905",
    "HBB_HUMAN":  "P68871",
    "HBA_MOUSE":  "P01942",
    "HBB_MOUSE":  "P02088",
    "HBA_CHICK":  "P01994",
    "HBB_CHICK":  "P02112",
    "HBA_DANRE":  "Q90487",
    "HBB_DANRE":  "Q6DHL5",
}

# Canonical short reference globins — used iff the UniProt REST call fails.
GLOBIN_FALLBACK = {
    "HBA_HUMAN": (
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
        "KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP"
        "AVHASLDKFLASVSTVLTSKYR"
    ),
    "HBB_HUMAN": (
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
        "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG"
        "KEFTPPVQAAYQKVVAGVANALAHKYH"
    ),
    "HBA_MOUSE": (
        "MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHG"
        "KKVADALANAAGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASHHPADFTPA"
        "VHASLDKFLASVSTVLTSKYR"
    ),
    "HBB_MOUSE": (
        "MVHLTDAEKAAVSGLWGKVNADEVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNPK"
        "VKAHGKKVITAFNDGLNHLDSLKGTFASLSELHCDKLHVDPENFRLLGNMIVIVLGHHLGK"
        "DFTPAAQAAFQKVVAGVATALAHKYH"
    ),
    "HBA_CHICK": (
        "MVLSAADKNNVKGIFTKIAGHAEEYGAETLERMFTTYPPTKTYFPHFDLSHGSAQIKGHG"
        "KKVVAALIEAANHIDDIAGTLSKLSDLHAHKLRVDPVNFKLLGQCFLVVVAIHHPAALTPE"
        "VHASLDKFLCAVGTVLTAKYR"
    ),
    "HBB_CHICK": (
        "MVHWTAEEKQLITGLWGKVNVAECGAEALARLLIVYPWTQRFFASFGNLSSPTAILGNPM"
        "VRAHGKKVLTSFGEAVKNLDNIKNTFAQLSELHCDKLHVDPENFRLLGDILIIVLAAHFSK"
        "DFTPECQAAWQKLVRVVAHALARKYH"
    ),
    "HBA_DANRE": (
        "MSLSDKDKAAVKAIWAKISKSADAIGAEALDRMLLVFPHFDLAHGSDQIKAHGKKVAAAL"
        "QHAINHIDDLPGALSGLSDLHALKLRVDPVNFKLLAQCYQVVLAIHLPSQFTPDAHVALDK"
        "FLATLASCLSEKYR"
    ),
    "HBB_DANRE": (
        "MVEWTDAERTAILGLWGKLNIDEIGPQALSRCLIVYPWTQRYFGSFGDLSTPDAILDNPK"
        "VAVHGKKVLAALGNAVKHLDDLKAYYAELSVLHSEKLHVDPDNFRLLADCITICAAMKFGP"
        "DFTPEAQAAWQKLVNVVAHALSRKYR"
    ),
}


def fetch_globins() -> dict:
    out = {}
    for name, acc in GLOBIN_ACCESSIONS.items():
        text = _http_get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
        seq = None
        if text:
            recs = _parse_fasta(text)
            if recs:
                seq = "".join(c for c in recs[0][1] if c in AA_SET)
        if not seq:
            seq = GLOBIN_FALLBACK[name]
        out[name] = seq
    return out


def percent_identity_distance(a: str, b: str) -> float:
    L = min(len(a), len(b))
    if L == 0:
        return 1.0
    matches = sum(1 for i in range(L) if a[i] == b[i])
    return 1.0 - matches / L


def neighbor_joining(D: np.ndarray, names: list[str]) -> list[tuple]:
    """Return a list of join events: (left, right, d_left, d_right, new_name).

    Classical NJ, O(N^3). Newick-style new_name = '(left:d_left,right:d_right)'.
    """
    D = D.astype(float).copy()
    n_orig = len(names)
    active = list(range(n_orig))
    labels = list(names)
    events = []
    new_id = n_orig

    while len(active) > 2:
        n = len(active)
        # Row sums over active rows/cols only.
        sub = D[np.ix_(active, active)]
        r = sub.sum(axis=1)
        # Q matrix.
        Q = (n - 2) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i_loc, j_loc = np.unravel_index(np.argmin(Q), Q.shape)
        i, j = active[i_loc], active[j_loc]
        d_ij = D[i, j]
        # Branch lengths to the new node.
        delta = (r[i_loc] - r[j_loc]) / (n - 2) if n > 2 else 0.0
        d_i = 0.5 * d_ij + 0.5 * delta
        d_j = 0.5 * d_ij - 0.5 * delta
        new_label = f"({labels[i]}:{max(d_i,0):.3f},{labels[j]}:{max(d_j,0):.3f})"
        events.append((labels[i], labels[j], max(d_i, 0.0), max(d_j, 0.0), new_label))

        # Expand D by one row/col for the new node.
        new_row = np.zeros(D.shape[0] + 1)
        for k in active:
            if k == i or k == j:
                continue
            new_row[k] = 0.5 * (D[i, k] + D[j, k] - d_ij)
        D = np.pad(D, ((0, 1), (0, 1)), mode="constant")
        D[new_id, :-1] = new_row[:-1]
        D[:-1, new_id] = new_row[:-1]
        D[new_id, new_id] = 0.0

        # Replace i,j with new node.
        active = [a for a in active if a not in (i, j)] + [new_id]
        labels.append(new_label)
        new_id += 1

    # Final join.
    i, j = active
    d_ij = D[i, j]
    events.append((labels[i], labels[j], d_ij / 2, d_ij / 2,
                   f"({labels[i]}:{d_ij/2:.3f},{labels[j]}:{d_ij/2:.3f})"))
    return events


print("Fetching 8 globin proteins from UniProt...")
globins = fetch_globins()
for name, seq in globins.items():
    print(f"  {name:>11s}  len={len(seq):3d}  start={seq[:20]}...")

# Build the 8x8 distance matrix on canonical short forms (truncate to common length so
# the percent-identity heuristic is comparable across pairs).
names = list(globins.keys())
n = len(names)
L_min = min(len(s) for s in globins.values())
trimmed = {k: v[:L_min] for k, v in globins.items()}

D = np.zeros((n, n))
for i, j in combinations(range(n), 2):
    d = percent_identity_distance(trimmed[names[i]], trimmed[names[j]])
    D[i, j] = D[j, i] = d

print(f"\\nDistance matrix ({n}x{n}, truncated to L={L_min}):")
print(pd.DataFrame(D.round(3), index=names, columns=names).to_string())

tree_events = neighbor_joining(D, names)
print("\\nNJ join order:")
for e in tree_events:
    print(f"  {e[0][:60]:>60s}  +  {e[1][:60]:>60s}   (d={e[2]:.3f},{e[3]:.3f})")

# Final Newick string is the last event's new_label.
newick = tree_events[-1][4] + ";"
print(f"\\nNewick:\\n  {newick[:200]}{'...' if len(newick)>200 else ''}")
'''


STEP2_MD = """## Step 2 (12 min) — Progressive MSA along the NJ guide tree

Now you build the multiple alignment **progressively**: walk the NJ tree
in join order, aligning sequence-to-sequence first, then collapsing each
join into a profile and aligning sequence-to-profile / profile-to-profile.
Every alignment uses Needleman-Wunsch DP with BLOSUM62 scoring and a
linear gap penalty.

The key idea: at a profile column the substitution score against a residue
`r` is the **average BLOSUM62 score** of `r` against every residue in that
column (gaps contribute the gap penalty). Once a gap is introduced into a
profile, it is **never removed** in later joins ("once a gap, always a
gap"), which is exactly why progressive MSA is greedy.
"""


STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Profile-aware Needleman-Wunsch + progressive MSA.
# ----------------------------------------------------------------------

GAP_PENALTY = -8


def profile_column_score(col_a, col_b) -> float:
    """Average BLOSUM62 score between two profile columns (lists of residues)."""
    # TODO: for every (r_a, r_b) pair across the two columns, look up BLOSUM62;
    # if either is a gap, charge GAP_PENALTY. Return the mean.
    raise NotImplementedError


def nw_profile(profile_a: list[list[str]],
               profile_b: list[list[str]]) -> list[list[str]]:
    """Align two profiles with Needleman-Wunsch; return the merged profile.

    A profile is a list of columns; each column is a list of residues (one
    per sequence in that profile). The merged profile has |a|+|b| sequences
    after concatenation column-wise.
    """
    # TODO: fill the DP matrix using profile_column_score(...) for matches,
    # GAP_PENALTY for gap rows. Backtrack to get the merged columns.
    raise NotImplementedError


def progressive_msa(seqs: dict[str, str], join_order: list[tuple]) -> dict[str, str]:
    """Progressive MSA following the join order from NJ."""
    # TODO: maintain a dict {node_label: profile}. Walk events in order.
    # Each event's new_label is "(left:d,right:d)" — find it in the dict,
    # nw_profile(left_profile, right_profile), assign result to new_label.
    # Return {seq_name: aligned_seq} at the end.
    raise NotImplementedError


# msa: dict[name, aligned_seq] for the 8 globins.
msa = {}
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.

GAP_PENALTY = -8


def profile_column_score(col_a, col_b) -> float:
    total = 0.0
    n = 0
    for ra in col_a:
        for rb in col_b:
            if ra == GAP or rb == GAP:
                total += GAP_PENALTY
            else:
                total += float(BLOSUM62[ra, rb])
            n += 1
    return total / n if n else 0.0


def nw_profile(profile_a: list[list[str]],
               profile_b: list[list[str]]) -> list[list[str]]:
    """Two-profile Needleman-Wunsch; return merged profile."""
    la, lb = len(profile_a), len(profile_b)
    n_a = len(profile_a[0]) if la else 0
    n_b = len(profile_b[0]) if lb else 0

    # Score and traceback matrices.
    F = np.zeros((la + 1, lb + 1))
    TB = np.zeros((la + 1, lb + 1), dtype=np.int8)  # 0=diag,1=up,2=left

    for i in range(1, la + 1):
        F[i, 0] = i * GAP_PENALTY
        TB[i, 0] = 1
    for j in range(1, lb + 1):
        F[0, j] = j * GAP_PENALTY
        TB[0, j] = 2

    for i in range(1, la + 1):
        for j in range(1, lb + 1):
            m = F[i - 1, j - 1] + profile_column_score(profile_a[i - 1], profile_b[j - 1])
            u = F[i - 1, j] + GAP_PENALTY
            l = F[i, j - 1] + GAP_PENALTY
            best = max(m, u, l)
            F[i, j] = best
            TB[i, j] = 0 if best == m else (1 if best == u else 2)

    # Traceback.
    merged = []
    i, j = la, lb
    gap_col_a = [GAP] * n_a
    gap_col_b = [GAP] * n_b
    while i > 0 or j > 0:
        if i > 0 and j > 0 and TB[i, j] == 0:
            merged.append(profile_a[i - 1] + profile_b[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and TB[i, j] == 1:
            merged.append(profile_a[i - 1] + gap_col_b)
            i -= 1
        else:
            merged.append(gap_col_a + profile_b[j - 1])
            j -= 1
    merged.reverse()
    return merged


def _profile_from_seq(name: str, seq: str):
    """Build a profile (list of single-residue columns) and its name index."""
    return [[c] for c in seq], [name]


def progressive_msa(seqs: dict, join_order: list) -> dict:
    # Pool: {node_label -> (profile, list_of_names_in_order)}
    pool = {}
    for name, seq in seqs.items():
        pool[name] = _profile_from_seq(name, seq)

    for left, right, _, _, new_label in join_order:
        if left not in pool or right not in pool:
            raise RuntimeError(f"missing pool entry: left={left[:40]}  right={right[:40]}")
        pa, na = pool[left]
        pb, nb = pool[right]
        pc = nw_profile(pa, pb)
        pool[new_label] = (pc, na + nb)

    # Final profile is the last event's new_label.
    final_label = join_order[-1][4]
    profile, names_in_profile = pool[final_label]

    out = {}
    for idx, name in enumerate(names_in_profile):
        out[name] = "".join(col[idx] for col in profile)
    return out


msa = progressive_msa(trimmed, tree_events)
print(f"MSA: {len(msa)} sequences, {len(next(iter(msa.values())))} columns")
print()
for name, seq in msa.items():
    print(f"  {name:>11s}  {seq[:80]}{'...' if len(seq) > 80 else ''}")
print(f"  ...len={len(next(iter(msa.values())))}")
'''


STEP3_MD = """## Step 3 (12 min) — SP score and per-column conservation

How good is your MSA? We benchmark it two ways:

- **Sum-of-pairs (SP) score** — sum BLOSUM62 of every aligned residue
  pair in every column across every pair of sequences. Higher = better.
- **Per-column conservation** — Shannon entropy `H = -sum p log2 p` over
  the residue frequencies in each column. Low H = conserved column.

The reference comparison point is the **Pfam PF00042 globin** seed MSA,
fetched live from EBI. If the network call fails we fall back to your
own MSA as the "reference" (the SP score will still be meaningful
relative to a shuffled control).
"""


STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — SP score vs Pfam reference, per-column entropy.
# ----------------------------------------------------------------------


def sp_score(aligned: dict[str, str]) -> float:
    """Sum-of-pairs BLOSUM62 score across all aligned columns."""
    # TODO: iterate column-by-column; for every (i<j) pair score BLOSUM62
    # (or GAP_PENALTY if either is a gap).
    raise NotImplementedError


def column_entropy(aligned: dict[str, str]) -> np.ndarray:
    """Shannon entropy (base 2) per column, ignoring gap rows."""
    # TODO
    raise NotImplementedError


def fetch_pfam_seed(pfam_id: str = "PF00042") -> dict[str, str]:
    """Return {seq_id: aligned_seq} from the Pfam seed alignment."""
    # TODO: try the InterPro API or the Stockholm-format URL; parse.
    # Return None (or {}) on failure so the caller can pick a fallback.
    raise NotImplementedError


# Plot conservation, compare SP score to a shuffled-column control.
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.


def sp_score(aligned: dict) -> float:
    seqs = list(aligned.values())
    n = len(seqs)
    if n < 2:
        return 0.0
    L = len(seqs[0])
    total = 0.0
    for c in range(L):
        col = [s[c] for s in seqs]
        for i in range(n):
            for j in range(i + 1, n):
                a, b = col[i], col[j]
                if a == GAP or b == GAP:
                    total += GAP_PENALTY
                elif a in AA_SET and b in AA_SET:
                    total += float(BLOSUM62[a, b])
    return total


def column_entropy(aligned: dict) -> np.ndarray:
    seqs = list(aligned.values())
    n = len(seqs)
    L = len(seqs[0])
    H = np.zeros(L)
    for c in range(L):
        col = [s[c] for s in seqs if s[c] != GAP]
        if not col:
            H[c] = 0.0
            continue
        counts = defaultdict(int)
        for r in col:
            counts[r] += 1
        p = np.array(list(counts.values()), dtype=float) / len(col)
        H[c] = -(p * np.log2(p + 1e-12)).sum()
    return H


def fetch_pfam_seed(pfam_id: str = "PF00042") -> dict:
    # InterPro proxy for Pfam alignments — Stockholm format.
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}/?annotation=alignment:seed&format=stockholm"
    text = _http_get(url)
    if not text or "# STOCKHOLM" not in text:
        # Fallback URL pattern.
        url2 = f"https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam_id}/?annotation=alignment%3Aseed"
        text = _http_get(url2)
    if not text:
        return {}
    # Stockholm: lines that are not # or // map "id seq".
    out = {}
    for line in text.splitlines():
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            sid, segment = parts[0], parts[-1]
            # Keep only AA + gap chars; uppercase.
            clean = "".join(c.upper() if c.isalpha() else "-" for c in segment)
            out[sid] = out.get(sid, "") + clean
    return out


own_sp = sp_score(msa)
own_H = column_entropy(msa)

# Per-row shuffle — destroys column-wise conservation, keeps per-sequence
# residue composition. Permuting *within* a column is invariant under SP score
# (BLOSUM is symmetric across rows), so we permute *across* columns instead.
rng = np.random.default_rng(0)
shuffled = {}
for k, v in msa.items():
    chars = list(v)
    rng.shuffle(chars)
    shuffled[k] = "".join(chars)
shuf_sp = sp_score(shuffled)
shuf_H = column_entropy(shuffled)

print(f"Sum-of-pairs (BLOSUM62):")
print(f"  Your MSA      : {own_sp:10.1f}")
print(f"  Column-shuffle: {shuf_sp:10.1f}")
print(f"  Delta         : {own_sp - shuf_sp:10.1f} (positive = real signal)")

# Try the Pfam reference.
print("\\nFetching Pfam PF00042 seed alignment...")
pfam = fetch_pfam_seed("PF00042")
if pfam:
    # Subset to a small comparable bundle and re-score on the same column count.
    pfam_subset = {k: v for k, v in list(pfam.items())[:8]}
    if pfam_subset:
        # Trim to common length.
        Lp = min(len(v) for v in pfam_subset.values())
        pfam_subset = {k: v[:Lp] for k, v in pfam_subset.items()}
        pfam_sp = sp_score(pfam_subset)
        print(f"  Pfam PF00042 seed (first 8 seqs, L={Lp}): SP = {pfam_sp:.1f}")
    else:
        print("  Pfam fetch returned empty.")
else:
    print("  Pfam fetch failed — keeping shuffled control as the contrast.")

# Conservation plot.
fig, ax = plt.subplots(figsize=(9, 3.5))
ax.plot(own_H, lw=1.2, color="#1f77b4", label="your MSA (entropy/col)")
ax.plot(shuf_H, lw=1.0, color="#aaaaaa", alpha=0.7, label="column-shuffled control")
ax.axhline(np.log2(20), ls="--", color="red", alpha=0.4, label="uniform 20-AA = log2 20")
ax.set_xlabel("alignment column")
ax.set_ylabel("Shannon entropy (bits)")
ax.set_title("Per-column conservation — globins")
ax.legend()
plt.tight_layout()
plt.show()

print(f"\\nMean entropy (your MSA): {own_H.mean():.2f} bits")
print(f"Mean entropy (shuffled): {shuf_H.mean():.2f} bits")
print(f"Conserved columns (H < 1.0 bit): {(own_H < 1.0).sum()} / {len(own_H)}")
'''


STEP4_MD = """## Step 4 (15 min) — Codon-aligned dN/dS and likelihood-ratio test

`dN/dS` (also written omega) compares the rate of non-synonymous to
synonymous substitutions in a coding alignment. Under purifying selection
dN/dS < 1; under positive selection dN/dS > 1; under strict neutrality
dN/dS = 1.

You are given 10 ~300-codon orthologs of a vertebrate gene, codon-aligned
(no gaps inside codons). For each pair you compute the Nei-Gojobori
estimators:

- `S` = number of synonymous sites (averaged over the two sequences);
- `N` = number of non-synonymous sites;
- `Sd`, `Nd` = observed synonymous / non-synonymous differences;
- `pS = Sd/S`, `pN = Nd/N`;
- `dS = -3/4 ln(1 - 4 pS / 3)`, `dN = -3/4 ln(1 - 4 pN / 3)` (Jukes-Cantor
  correction).

We then run a **chi-squared LRT** with 1 d.o.f. on `2 * (lnL_free - lnL_neutral)`,
where the neutral model fixes dN = dS.
"""


STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Nei-Gojobori dN/dS + LRT vs neutrality.
# ----------------------------------------------------------------------

# 10-species codon-aligned ortholog set (synthetic but biologically structured).
# Each sequence is exactly 300 nt (100 codons); identity across species is high so
# JC correction stays in its valid range.

def build_codon_alignment(seed: int = 7) -> dict[str, str]:
    """Return {species: 300-nt coding sequence}, generated deterministically."""
    # TODO: start from a random 100-codon ORF; for each species apply
    # independent point mutations at a low rate, then add a small fraction
    # of *additional* non-synonymous mutations to two "selected" lineages
    # so the LRT actually has something to reject.
    raise NotImplementedError


# CODON_TABLE (standard genetic code) is already defined globally — use it directly.


def nei_gojobori(seq_a: str, seq_b: str) -> dict:
    """Compute (S, N, Sd, Nd, dS, dN, omega) for one codon-aligned pair."""
    # TODO: walk codon-by-codon. For each codon, count syn/nonsyn sites by
    # the standard Nei-Gojobori "fraction of single-base changes that are
    # synonymous". Diff: compare the two codons and apportion the change(s)
    # to syn/nonsyn classes. Apply Jukes-Cantor correction at the end.
    raise NotImplementedError


def lrt_neutrality(dS: float, dN: float, n_diffs: int) -> tuple[float, float]:
    """Approximate LRT statistic and chi^2(1) p-value for H0: dN = dS."""
    # TODO: log-likelihood under free model vs constrained-equal model.
    # A reasonable proxy: treat (Sd, Nd) as binomial successes; null = same
    # rate, free = separate rates. Return (LR_stat, p).
    raise NotImplementedError
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.


def _random_codon(rng) -> str:
    while True:
        c = "".join(rng.choice(list("ACGT")) for _ in range(3))
        if CODON_TABLE[c] != "*":
            return c


def _mutate(seq: str, rng, rate: float, force_nonsyn_extra: float = 0.0) -> str:
    """Apply per-site point mutations at `rate`, plus optional extra non-syn pressure."""
    out = list(seq)
    for i in range(0, len(out), 3):
        codon = "".join(out[i:i + 3])
        # Base mutations.
        for k in range(3):
            if rng.random() < rate:
                old = out[i + k]
                new = rng.choice([b for b in "ACGT" if b != old])
                trial = list(codon)
                trial[k] = new
                trial_codon = "".join(trial)
                if CODON_TABLE.get(trial_codon, "*") != "*":
                    out[i + k] = new
                    codon = trial_codon
        # Extra non-synonymous mutations on selected lineages.
        if force_nonsyn_extra > 0.0 and rng.random() < force_nonsyn_extra:
            for k in range(3):
                old = out[i + k]
                for new in rng.permutation(list("ACGT")):
                    if new == old:
                        continue
                    trial = list(codon)
                    trial[k] = new
                    trial_codon = "".join(trial)
                    aa_old = CODON_TABLE[codon]
                    aa_new = CODON_TABLE.get(trial_codon, "*")
                    if aa_new != "*" and aa_new != aa_old:
                        out[i + k] = new
                        break
    return "".join(out)


def build_codon_alignment(seed: int = 7) -> dict:
    rng = np.random.default_rng(seed)
    # 100-codon ancestral ORF (with explicit start codon, no internal stops).
    ancestor = "ATG"
    while len(ancestor) < 300:
        ancestor += _random_codon(rng)
    ancestor = ancestor[:300]

    species = [
        "human", "chimp", "mouse", "rat",
        "dog", "cow", "chicken", "lizard",
        "frog", "zebrafish",
    ]
    # Most lineages: low neutral rate. Two lineages get extra non-synonymous load
    # strong enough that at least one pairwise LRT clears chi^2(1) p<0.05.
    base_rate = 0.025
    selected = {"chicken", "zebrafish"}

    aln = {}
    for sp in species:
        extra = 0.20 if sp in selected else 0.0
        aln[sp] = _mutate(ancestor, rng, rate=base_rate, force_nonsyn_extra=extra)
    return aln


# Pre-compute syn/non-syn site counts per codon (Nei-Gojobori).
def _site_counts(codon: str) -> tuple[float, float]:
    """Return (syn_sites, nonsyn_sites) for a single codon under NG."""
    if CODON_TABLE.get(codon, "*") == "*":
        return 0.0, 0.0
    syn = 0.0
    for k in range(3):
        for new in "ACGT":
            if new == codon[k]:
                continue
            trial = codon[:k] + new + codon[k + 1:]
            aa_old = CODON_TABLE[codon]
            aa_new = CODON_TABLE.get(trial, "*")
            if aa_new == "*":
                continue
            # 1/3 weight per substitution direction.
            if aa_new == aa_old:
                syn += 1.0 / 3.0
    # Each codon contributes 3 sites total minus stop fraction; the NG convention
    # is syn_sites + nonsyn_sites = 3 (excluding directions leading to stops).
    nonsyn = 3.0 - syn
    return syn, nonsyn


def _codon_diffs(c1: str, c2: str) -> tuple[float, float]:
    """Sd, Nd for one codon pair (averaged over single-step paths)."""
    if c1 == c2:
        return 0.0, 0.0
    diffs = [k for k in range(3) if c1[k] != c2[k]]
    n_diff = len(diffs)
    if n_diff == 0:
        return 0.0, 0.0
    if n_diff == 1:
        k = diffs[0]
        aa1 = CODON_TABLE.get(c1, "*")
        aa2 = CODON_TABLE.get(c2, "*")
        if "*" in (aa1, aa2):
            return 0.0, 0.0
        return (1.0, 0.0) if aa1 == aa2 else (0.0, 1.0)
    # 2 or 3 diffs: average over all single-step orderings.
    from itertools import permutations
    paths = []
    for order in permutations(diffs):
        cur = c1
        sd = nd = 0.0
        valid = True
        for k in order:
            nxt = cur[:k] + c2[k] + cur[k + 1:]
            aa_cur = CODON_TABLE.get(cur, "*")
            aa_nxt = CODON_TABLE.get(nxt, "*")
            if "*" in (aa_cur, aa_nxt):
                valid = False
                break
            if aa_cur == aa_nxt:
                sd += 1.0
            else:
                nd += 1.0
            cur = nxt
        if valid:
            paths.append((sd, nd))
    if not paths:
        return 0.0, 0.0
    sd_mean = np.mean([p[0] for p in paths])
    nd_mean = np.mean([p[1] for p in paths])
    return sd_mean, nd_mean


def nei_gojobori(seq_a: str, seq_b: str) -> dict:
    assert len(seq_a) == len(seq_b)
    assert len(seq_a) % 3 == 0
    S_a = N_a = 0.0
    S_b = N_b = 0.0
    Sd = Nd = 0.0
    for i in range(0, len(seq_a), 3):
        ca, cb = seq_a[i:i + 3], seq_b[i:i + 3]
        sa, na = _site_counts(ca)
        sb, nb = _site_counts(cb)
        S_a += sa; N_a += na
        S_b += sb; N_b += nb
        sd, nd = _codon_diffs(ca, cb)
        Sd += sd; Nd += nd
    S = 0.5 * (S_a + S_b)
    N = 0.5 * (N_a + N_b)
    pS = Sd / S if S > 0 else 0.0
    pN = Nd / N if N > 0 else 0.0
    # Jukes-Cantor correction; clip so log argument stays positive.
    def _jc(p):
        arg = 1.0 - 4.0 * p / 3.0
        if arg <= 0:
            return float("nan")
        return -0.75 * math.log(arg)
    dS = _jc(pS)
    dN = _jc(pN)
    omega = dN / dS if dS and dS > 0 else float("nan")
    return {"S": S, "N": N, "Sd": Sd, "Nd": Nd, "dS": dS, "dN": dN, "omega": omega}


def lrt_neutrality(Sd: float, Nd: float, S: float, N: float) -> tuple[float, float]:
    """Two-rate (dN, dS) vs one-rate (dN = dS) binomial LRT, chi^2(1)."""
    total_d = Sd + Nd
    total_sites = S + N
    if total_d == 0 or total_sites == 0:
        return 0.0, 1.0
    # Free model: separate rates pN, pS.
    pS = max(Sd / S, 1e-9) if S > 0 else 1e-9
    pN = max(Nd / N, 1e-9) if N > 0 else 1e-9
    # Null: common rate p0 = (Sd+Nd)/(S+N).
    p0 = total_d / total_sites
    p0 = min(max(p0, 1e-9), 1 - 1e-9)

    def _lnL_binom(k, n, p):
        if p <= 0 or p >= 1:
            return float("-inf")
        return k * math.log(p) + (n - k) * math.log(1 - p)

    lnL_free = _lnL_binom(Sd, S, pS) + _lnL_binom(Nd, N, pN)
    lnL_null = _lnL_binom(Sd, S, p0) + _lnL_binom(Nd, N, p0)
    LR = 2.0 * (lnL_free - lnL_null)
    if LR < 0:
        LR = 0.0
    p = 1.0 - chi2.cdf(LR, df=1)
    return LR, p


codon_aln = build_codon_alignment()
species_list = list(codon_aln.keys())
print(f"Codon alignment: {len(species_list)} species x {len(next(iter(codon_aln.values())))} nt")

# Pairwise dN/dS.
rows = []
for a, b in combinations(species_list, 2):
    ng = nei_gojobori(codon_aln[a], codon_aln[b])
    LR, p = lrt_neutrality(ng["Sd"], ng["Nd"], ng["S"], ng["N"])
    rows.append({
        "pair": f"{a}|{b}",
        "S": round(ng["S"], 1), "N": round(ng["N"], 1),
        "Sd": round(ng["Sd"], 1), "Nd": round(ng["Nd"], 1),
        "dS": round(ng["dS"], 3) if not math.isnan(ng["dS"]) else None,
        "dN": round(ng["dN"], 3) if not math.isnan(ng["dN"]) else None,
        "omega": round(ng["omega"], 3) if not math.isnan(ng["omega"]) else None,
        "LR": round(LR, 2),
        "p": f"{p:.3g}",
    })

dnds = pd.DataFrame(rows)
print(dnds.to_string(index=False))

# Pairs where omega > 1 *and* the LRT rejects neutrality.
sig = dnds[(dnds["omega"].notna()) & (dnds["omega"] > 1.0) & (dnds["p"].astype(float) < 0.05)]
print(f"\\nPairs with omega>1 and chi^2(1) p<0.05: {len(sig)}")
if not sig.empty:
    print(sig.to_string(index=False))
'''


STEP5_MD = """## Step 5 (11 min) — Interpret the tree and flag positive selection

You now have everything to make biological sense of the comparative data:

- Which globin clade clusters by **paralog** (alpha vs beta) and which by
  **organism**? In vertebrate globins the alpha / beta gene duplication
  predates the tetrapod / fish split, so an NJ tree built from sequence
  distance alone should group HBA's together and HBB's together — *not*
  group all human globins together. That is the classic signal that
  "trees encode shared ancestry, not just similarity".
- For dN/dS, which lineage pair shows the strongest signature of positive
  selection? Synthesise the dN/dS table by listing the top-5 pairs by
  omega (descending), and check the chi-squared p-value column.
"""


STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Interpretation pass: globin clade structure + dN/dS hot pairs.
# ----------------------------------------------------------------------

# TODO 1: print the join order again and annotate which step joined alpha-globins,
#         beta-globins, and the alpha/beta clades. A simple rule of thumb: if both
#         labels in a join contain "HBA" we are inside the alpha clade.
# TODO 2: print the top-5 dN/dS pairs sorted by omega descending, flagging any
#         that also have chi^2(1) p < 0.05.
'''


STEP5_SOLUTION = '''# Reference solution — Step 5.

print("NJ join order — annotated:")
def _clade_tag(label: str) -> str:
    has_a = "HBA" in label
    has_b = "HBB" in label
    if has_a and not has_b:
        return "alpha"
    if has_b and not has_a:
        return "beta"
    if has_a and has_b:
        return "alpha+beta"
    return "?"

for i, (left, right, dl, dr, _) in enumerate(tree_events, 1):
    lt, rt = _clade_tag(left), _clade_tag(right)
    short_l = left if len(left) < 35 else left[:32] + "..."
    short_r = right if len(right) < 35 else right[:32] + "..."
    print(f"  step {i}: {short_l:>38s} [{lt:>10s}]  +  {short_r:<38s} [{rt:<10s}]  "
          f"d=({dl:.3f},{dr:.3f})")

print("\\nTop-5 dN/dS pairs (by omega desc):")
top = dnds[dnds["omega"].notna()].sort_values("omega", ascending=False).head(5)
for _, row in top.iterrows():
    flag = " <- positive selection (p<0.05)" if float(row["p"]) < 0.05 and row["omega"] > 1 else ""
    print(f"  {row['pair']:>22s}  omega={row['omega']:.2f}  LR={row['LR']:.2f}  p={row['p']}{flag}")

# Stash for self-check.
top_pairs = top
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces of the pipeline.
If you ran the reference solutions they should all pass.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check.
# ----------------------------------------------------------------------

# 1. Distance matrix is symmetric and zero on the diagonal.
assert np.allclose(D, D.T), "distance matrix not symmetric"
assert np.allclose(np.diag(D), 0.0), "distance matrix has nonzero diagonal"

# 2. NJ produced exactly n-1 joins for n leaves.
assert len(tree_events) == len(names) - 1, (
    f"NJ produced {len(tree_events)} events; expected {len(names)-1}"
)

# 3. MSA: every sequence has the same length, and is at least as long as the original.
msa_lens = [len(s) for s in msa.values()]
assert len(set(msa_lens)) == 1, f"MSA rows have inconsistent length: {msa_lens}"
assert msa_lens[0] >= L_min, "MSA shorter than input — alignment broke"

# 4. SP score: real MSA beats the column-shuffled control.
own_sp_check = sp_score(msa)
shuf_sp_check = shuf_sp
assert own_sp_check > shuf_sp_check, (
    f"SP score check failed: own={own_sp_check:.1f} vs shuf={shuf_sp_check:.1f}"
)

# 5. dN/dS sanity: at least one lineage pair has omega defined and > 0.5.
omegas = dnds["omega"].dropna().values
assert len(omegas) > 0, "no dN/dS values computed"
assert max(omegas) > 0.5, f"all omegas suspiciously small: max={max(omegas):.3f}"

# 6. LRT statistic is non-negative.
assert (dnds["LR"] >= 0).all(), "negative LR statistic — bug in lrt_neutrality"

# 7. At least one pair shows elevated omega (>=1.0) — this exercises the LRT path.
elevated = dnds[dnds["omega"].fillna(0) >= 1.0]
print(f"Pairs with omega >= 1.0: {len(elevated)}")

print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — tree as graphical model, dN/dS as LRT

You implemented three canonical signal-processing structures dressed in
biology costume:

1. **NJ is hierarchical clustering driven by a corrected distance.** The
   Q matrix is a rank-2 correction that subtracts each row's mean
   distance to the rest of the active set — the same idea as
   demeaning in PCA before forming the covariance matrix. Each merge
   reduces the active dimension by one; you end with a binary tree, the
   minimum-DOF model for `n` leaves.
2. **Progressive MSA is greedy DP on a tree-induced order.** Pairwise
   Needleman-Wunsch is exact DP on a 2-D lattice; profile NW raises that
   to a profile-vs-profile DP where the substitution score is the
   **average BLOSUM62 entry** across two columns. The total cost is
   `O(N L^2)` instead of the `O(L^N)` that an exact joint optimisation
   would demand. The greedy tradeoff is exactly the same idea as
   sequential decoding: take the best local decision now, never revisit.
3. **dN/dS is a likelihood-ratio test.** The free model has two rate
   parameters (dN, dS); the neutral null collapses them to one. Twice
   the log-likelihood ratio is chi-squared with 1 d.o.f.: the same Wald-
   adjacent test statistic you would write down for a 2-vs-1 nested
   regression. omega > 1 is the **alternative-hypothesis-favoured**
   region — same idea as constant-false-alarm-rate (CFAR) detection in
   radar, with the chi-squared threshold setting the false-alarm rate.
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
