"""Build exercise.ipynb for L19 — BLAST and Sequence Search Statistics.

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


TITLE_MD = """# L19 — BLAST and Sequence Search Statistics

In this exercise you implement the core of protein BLAST: a seed-and-extend
cascade that scores candidate alignments with BLOSUM62, then converts raw
scores into bit scores and Karlin-Altschul E-values. You will run your
pipeline on five real query proteins against a ~20k-sequence UniProt subset
and compare hit lists to what NCBI BLAST returns.
"""


AHA_MD = """> **Aha.** BLAST is a **matched-filter cascade**: a cheap k-mer prefilter
> hands a small set of candidate diagonals to an expensive BLOSUM-scored
> ungapped extension. The bit score is a **log-likelihood ratio** test
> statistic comparing the target model to the random background; the
> E-value is its **false-alarm rate** scaled by the search space `m·n`.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 biopython==1.83 matplotlib==3.8.4 requests==2.32.3 -q
"""


IMPORTS = """import io
import math
import time
import random
import urllib.request
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests

from Bio.Align import substitution_matrices

# Deterministic for the whole notebook.
np.random.seed(42)
random.seed(42)

# Canonical 20-amino-acid alphabet (no ambiguous codes).
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA)

# Load BLOSUM62 once; the matrix is a Bio.Align.substitution_matrices.Array.
BLOSUM62 = substitution_matrices.load("BLOSUM62")
print(f"BLOSUM62 loaded: {BLOSUM62.shape} entries; A/A={int(BLOSUM62['A','A'])}, W/W={int(BLOSUM62['W','W'])}")


def score_pair(a: str, b: str) -> int:
    \"\"\"BLOSUM62 score for two equal-length residue strings (X / unknown -> 0).\"\"\"
    s = 0
    for x, y in zip(a, b):
        if x in AA_SET and y in AA_SET:
            s += int(BLOSUM62[x, y])
    return s
"""


STEP1_MD = """## Step 1 (10 min) — Build the database and inspect BLOSUM62

We need (a) a small protein search space and (b) the substitution matrix that
turns residue pairs into log-odds scores. We fetch ~5 000 reviewed human
proteins from UniProt; if the network call fails we fall back to a deterministic
synthetic database built from amino-acid background frequencies so the rest of
the notebook still works offline. You will also hand-score one alignment with
BLOSUM62 to make sure the matrix lookup matches a hand calculation.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Load / build the protein database and the query set.
# ----------------------------------------------------------------------

QUERY_ACCESSIONS = {
    "INS_HUMAN":   "P01308",  # human insulin
    "HBA_HUMAN":   "P69905",  # human hemoglobin alpha
    "LYZ_CHICK":   "P00698",  # hen egg-white lysozyme
    "BRCA1_HUMAN": "P38398",  # human BRCA1
    "MYG_HUMAN":   "P02144",  # human myoglobin
}

def fetch_uniprot_fasta(accession: str) -> tuple[str, str]:
    """Return (header, sequence) for one UniProt accession."""
    # TODO: GET https://rest.uniprot.org/uniprotkb/{accession}.fasta
    # Parse out the header line and concatenate the remaining lines.
    raise NotImplementedError


def fetch_uniprot_db(n: int = 5000) -> dict[str, str]:
    """Return {accession: sequence} for ~n reviewed human proteins."""
    # TODO: hit https://rest.uniprot.org/uniprotkb/search with
    # query=reviewed:true AND organism_id:9606, format=fasta, size=n
    # Parse the multi-FASTA stream; return a dict.
    raise NotImplementedError


def synthetic_db(n: int = 5000, mean_len: int = 350, seed: int = 0) -> dict[str, str]:
    """Deterministic fallback: n random proteins drawn from AA background freqs."""
    # TODO: sample lengths around mean_len; sample residues from a fixed frequency vector.
    raise NotImplementedError


# Hand-score one alignment to sanity-check BLOSUM62.
def hand_score_example() -> int:
    """Score the alignment 'WKAA' vs 'WRAA' under BLOSUM62."""
    # TODO: use score_pair(...) and assert the answer matches your by-hand sum.
    raise NotImplementedError


# queries: dict[name, seq]; db: dict[acc, seq]
queries = {}
db = {}
'''


STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.
import urllib.request

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
UNIPROT_SEARCH = (
    "https://rest.uniprot.org/uniprotkb/search?"
    "query=reviewed:true+AND+organism_id:9606&format=fasta&size={n}"
)

QUERY_ACCESSIONS = {
    "INS_HUMAN":   "P01308",
    "HBA_HUMAN":   "P69905",
    "LYZ_CHICK":   "P00698",
    "BRCA1_HUMAN": "P38398",
    "MYG_HUMAN":   "P02144",
}


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


def _http_get(url: str, timeout: float = 30.0) -> str | None:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "L19-exercise/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except Exception as exc:
        print(f"  network error on {url[:60]}...: {exc}")
        return None


def fetch_uniprot_fasta(accession: str):
    text = _http_get(UNIPROT_FASTA.format(acc=accession))
    if text is None:
        return None
    recs = _parse_fasta(text)
    if not recs:
        return None
    return recs[0]


def fetch_uniprot_db(n: int = 5000):
    text = _http_get(UNIPROT_SEARCH.format(n=n))
    if text is None:
        return None
    out = {}
    for header, seq in _parse_fasta(text):
        # header looks like "sp|P12345|GENE_HUMAN Description ..."
        acc = header.split("|")[1] if "|" in header else header.split()[0]
        # Keep only canonical AAs; drop very short fragments.
        clean = "".join(c for c in seq if c in AA_SET)
        if len(clean) >= 30:
            out[acc] = clean
    return out


# Synthetic fallback — deterministic, no network.
# Background amino-acid frequencies (Robinson & Robinson 1991 / UniProt averages).
_AA_FREQ = np.array([
    0.074, 0.025, 0.054, 0.054, 0.047, 0.074, 0.026, 0.068,
    0.058, 0.099, 0.025, 0.045, 0.039, 0.034, 0.052, 0.057,
    0.051, 0.073, 0.013, 0.032,
])
_AA_FREQ = _AA_FREQ / _AA_FREQ.sum()


def synthetic_db(n: int = 5000, mean_len: int = 350, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    db = {}
    for i in range(n):
        L = max(50, int(rng.normal(mean_len, 80)))
        seq = "".join(rng.choice(list(AA), size=L, p=_AA_FREQ))
        db[f"SYN{i:05d}"] = seq
    return db


def hand_score_example() -> int:
    # WKAA vs WRAA -> W/W + K/R + A/A + A/A
    s = score_pair("WKAA", "WRAA")
    # By-hand: W/W=11, K/R=2, A/A=4, A/A=4 -> 21
    assert s == 21, f"Hand-scored alignment expected 21, got {s}"
    return s


print("Fetching 5 query proteins from UniProt...")
queries = {}
for name, acc in QUERY_ACCESSIONS.items():
    rec = fetch_uniprot_fasta(acc)
    if rec is None:
        # Tiny built-in fallback for the queries — minimal canonical sequences.
        FALLBACK = {
            "P01308": (
                "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
                "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
            ),
            "P69905": (
                "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
                "KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPA"
                "VHASLDKFLASVSTVLTSKYR"
            ),
            "P00698": (
                "MKALIVLGLVLLSVTVQGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQA"
                "TNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDG"
                "NGMNAWVAWRNRCKGTDVQAWIRGCRL"
            ),
            "P38398": (
                "MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQ"
                "CPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKEN"
            ),
            "P02144": (
                "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASE"
                "DLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHP"
                "GDFGADAQGAMNKALELFRKDMASNYKELGFQG"
            ),
        }
        rec = (f"sp|{acc}|{name}", FALLBACK.get(acc, ""))
    header, seq = rec
    queries[name] = "".join(c for c in seq if c in AA_SET)
    print(f"  {name:>12s}  {acc}  len={len(queries[name])}")

print("Fetching ~5000 reviewed human proteins...")
db = fetch_uniprot_db(5000)
if db is None or len(db) < 500:
    print("  REST fetch failed or too small; using synthetic fallback.")
    db = synthetic_db(5000)
print(f"Database: {len(db)} sequences, {sum(len(s) for s in db.values()):,} residues total")

print("Hand-scored alignment WKAA vs WRAA:", hand_score_example())
'''


STEP2_MD = """## Step 2 (12 min) — Query 3-mer neighbourhoods

BLAST does not scan the database with the query's own 3-mers; it builds a
**neighbourhood**: for each query word, every 3-mer in the AA³ space whose
BLOSUM62 score against the query word is at least the seeding threshold `T`.
A higher `T` shrinks the neighbourhood (fewer seeds, faster, less sensitive);
a lower `T` grows it. You will enumerate neighbourhoods at `T = 11` (BLAST
default) and `T = 5` (very permissive) and plot the size distribution.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — k-mer neighbourhoods for a query at two thresholds.
# ----------------------------------------------------------------------

K = 3


def neighbourhood(word: str, T: int) -> list[tuple[str, int]]:
    """All (k-mer, score) pairs with BLOSUM62 score >= T against `word`."""
    # TODO: iterate all 20**K possible k-mers and keep those that beat T.
    # For K=3 there are 8000 — fine.
    raise NotImplementedError


def query_neighbourhoods(seq: str, T: int) -> dict[int, list[tuple[str, int]]]:
    """For each window-start position in `seq`, the neighbourhood list."""
    # TODO: slide a length-K window over seq; call neighbourhood(...) at each pos.
    raise NotImplementedError


# Compare neighbourhood sizes at T=11 vs T=5 for one query (e.g. INS_HUMAN).
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.
from itertools import product

K = 3


def _all_kmers(k: int = K):
    # Cached on first call.
    if not hasattr(_all_kmers, "_cache"):
        _all_kmers._cache = ["".join(p) for p in product(AA, repeat=k)]
    return _all_kmers._cache


def neighbourhood(word: str, T: int) -> list[tuple[str, int]]:
    out = []
    for kmer in _all_kmers(len(word)):
        s = score_pair(word, kmer)
        if s >= T:
            out.append((kmer, s))
    return out


def query_neighbourhoods(seq: str, T: int) -> dict:
    return {
        i: neighbourhood(seq[i:i + K], T)
        for i in range(len(seq) - K + 1)
    }


# Compare T=11 vs T=5 on insulin.
ins = queries["INS_HUMAN"]
nbhd_high = query_neighbourhoods(ins, T=11)
nbhd_low  = query_neighbourhoods(ins, T=5)

sizes_high = [len(v) for v in nbhd_high.values()]
sizes_low  = [len(v) for v in nbhd_low.values()]

print(f"INS_HUMAN ({len(ins)} aa)")
print(f"  T=11: mean nbhd size = {np.mean(sizes_high):5.1f}, total seeds = {sum(sizes_high):,}")
print(f"  T= 5: mean nbhd size = {np.mean(sizes_low):5.1f}, total seeds = {sum(sizes_low):,}")

fig, ax = plt.subplots(figsize=(7, 3.5))
ax.hist(sizes_high, bins=30, alpha=0.7, label="T = 11 (default)")
ax.hist(sizes_low,  bins=30, alpha=0.5, label="T =  5 (permissive)")
ax.set_xlabel("neighbourhood size per query 3-mer")
ax.set_ylabel("count of query positions")
ax.set_title("Seed neighbourhood size — insulin")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP3_MD = """## Step 3 (12 min) — Scan database and X-drop ungapped extension

The database is indexed by its 3-mers; each query-side neighbourhood word
hits a list of (db-seq, db-position) pairs in O(1). At each hit we anchor an
**ungapped extension**: walk outward in both directions, adding BLOSUM62
column scores, and stop when the running score has dropped more than `X`
below the best score seen so far. The maximum score over a seed's
extension is one **HSP** (high-scoring segment pair). Report the top-10
HSPs per query.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — DB k-mer index, seeding scan, X-drop extension.
# ----------------------------------------------------------------------

T_SEED  = 11   # neighbourhood threshold
X_DROP  = 7    # bits-ish ungapped extension drop
TOP_N   = 10   # HSPs to keep per query


def build_db_index(db: dict[str, str], k: int = K) -> dict[str, list[tuple[str, int]]]:
    """{kmer: [(db_acc, pos), ...]}"""
    # TODO: for every db sequence, slide a length-k window and append to the index.
    raise NotImplementedError


def xdrop_extend(query: str, target: str, q_pos: int, t_pos: int, X: int = X_DROP) -> tuple[int, int, int]:
    """Return (best_score, q_start, t_start) for an ungapped extension.

    Both endpoints inclusive of the seed; extend left, then right; stop when
    running score has fallen more than X below the best running score.
    """
    # TODO: implement two-sided X-drop walk.
    raise NotImplementedError


def blast_one_query(query: str, db: dict, db_index: dict,
                    T: int = T_SEED, X: int = X_DROP, top_n: int = TOP_N) -> list[dict]:
    """Return a list of HSP dicts sorted by score desc."""
    # TODO:
    # 1. Build query neighbourhoods at threshold T.
    # 2. For each (q_pos, neighbour_word) lookup db_index[neighbour_word].
    # 3. For each hit, X-drop extend; remember best HSP per (db_acc, diagonal).
    # 4. Sort and return top_n.
    raise NotImplementedError


# db_index = build_db_index(db)
# hsps_per_query = {name: blast_one_query(seq, db, db_index) for name, seq in queries.items()}
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.

T_SEED  = 11
X_DROP  = 7
TOP_N   = 10


def build_db_index(db: dict, k: int = K) -> dict:
    idx = defaultdict(list)
    for acc, seq in db.items():
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            # Skip windows containing non-canonical residues.
            if all(c in AA_SET for c in kmer):
                idx[kmer].append((acc, i))
    return dict(idx)


def xdrop_extend(query: str, target: str, q_pos: int, t_pos: int,
                 X: int = X_DROP) -> tuple[int, int, int, int]:
    """Two-sided X-drop ungapped extension.

    Returns (best_score, q_start, t_start, length).
    """
    # Seed score: 3-mer match — start with the K-residue seed.
    seed_score = score_pair(query[q_pos:q_pos + K], target[t_pos:t_pos + K])
    best = seed_score
    q_start, t_start, length = q_pos, t_pos, K

    # Extend left.
    score = best
    qi, ti = q_pos - 1, t_pos - 1
    cur = best
    cur_q_start, cur_t_start, cur_len = q_pos, t_pos, K
    while qi >= 0 and ti >= 0:
        a, b = query[qi], target[ti]
        if a not in AA_SET or b not in AA_SET:
            break
        cur += int(BLOSUM62[a, b])
        cur_len += 1
        cur_q_start, cur_t_start = qi, ti
        if cur > best:
            best = cur
            q_start, t_start, length = cur_q_start, cur_t_start, cur_len
        if best - cur > X:
            break
        qi -= 1
        ti -= 1

    # Extend right.
    cur = best
    qi = q_start + length
    ti = t_start + length
    while qi < len(query) and ti < len(target):
        a, b = query[qi], target[ti]
        if a not in AA_SET or b not in AA_SET:
            break
        cur += int(BLOSUM62[a, b])
        if cur > best:
            best = cur
            length = qi - q_start + 1
        if best - cur > X:
            break
        qi += 1
        ti += 1

    return best, q_start, t_start, length


def blast_one_query(query: str, db: dict, db_index: dict,
                    T: int = T_SEED, X: int = X_DROP,
                    top_n: int = TOP_N) -> list[dict]:
    # 1. Query neighbourhoods.
    nbhd = query_neighbourhoods(query, T)

    # 2. Track best HSP per (db_acc, diagonal) to deduplicate seeds on same diagonal.
    best_by_diag: dict[tuple[str, int], dict] = {}

    for q_pos, words in nbhd.items():
        for word, _ in words:
            hits = db_index.get(word)
            if not hits:
                continue
            for acc, t_pos in hits:
                target = db[acc]
                score, qs, ts, length = xdrop_extend(query, target, q_pos, t_pos, X)
                diag = ts - qs
                key = (acc, diag)
                prev = best_by_diag.get(key)
                if prev is None or score > prev["score"]:
                    best_by_diag[key] = {
                        "acc": acc, "score": score,
                        "q_start": qs, "t_start": ts, "length": length,
                        "diag": diag,
                    }

    hsps = sorted(best_by_diag.values(), key=lambda d: d["score"], reverse=True)
    return hsps[:top_n]


t0 = time.time()
db_index = build_db_index(db)
print(f"DB k-mer index: {len(db_index):,} unique 3-mers, built in {time.time()-t0:.1f}s")

hsps_per_query = {}
for name, seq in queries.items():
    t1 = time.time()
    hsps_per_query[name] = blast_one_query(seq, db, db_index)
    print(f"  {name}: top score = {hsps_per_query[name][0]['score']:4d}  "
          f"({time.time()-t1:.1f}s, {len(hsps_per_query[name])} HSPs)")

# Show the top-5 HSPs for one query.
print("\\nTop-5 HSPs for INS_HUMAN:")
for h in hsps_per_query["INS_HUMAN"][:5]:
    print(f"  acc={h['acc']}  score={h['score']:4d}  "
          f"q[{h['q_start']}:{h['q_start']+h['length']}]  "
          f"t[{h['t_start']}:{h['t_start']+h['length']}]")
'''


STEP4_MD = """## Step 4 (15 min) — Bit scores and Karlin-Altschul E-values

The raw alignment score `S` depends on the chosen substitution matrix. The
**bit score** normalises it so that one bit means a 2× likelihood ratio:
$$S' = \\frac{\\lambda S - \\ln K}{\\ln 2}$$
For ungapped BLOSUM62 the canonical parameters are λ ≈ 0.318 and K ≈ 0.13.
The **E-value** is the expected number of HSPs with bit score ≥ S' under the
random-string null:
$$E = m \\cdot n \\cdot 2^{-S'}$$
with `m` = query length and `n` = total DB length. Compute both for every
top-10 HSP and assemble the results into a pandas DataFrame.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Karlin-Altschul bit scores and E-values.
# ----------------------------------------------------------------------

# Canonical ungapped BLOSUM62 parameters.
LAMBDA = 0.318
K_KA   = 0.13


def bit_score(raw_score: int, lam: float = LAMBDA, K: float = K_KA) -> float:
    """Karlin-Altschul bit-score conversion."""
    # TODO
    raise NotImplementedError


def e_value(bits: float, m: int, n: int) -> float:
    """Expected number of chance HSPs >= these bits in an m vs n search."""
    # TODO
    raise NotImplementedError


# Assemble a DataFrame with columns:
#   query, acc, raw_score, bits, E, q_start, t_start, length
# rows = union of top-10 HSPs across the 5 queries.
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.

LAMBDA = 0.318
K_KA   = 0.13


def bit_score(raw_score: int, lam: float = LAMBDA, K: float = K_KA) -> float:
    return (lam * raw_score - math.log(K)) / math.log(2.0)


def e_value(bits: float, m: int, n: int) -> float:
    return m * n * (2.0 ** -bits)


n_db_residues = sum(len(s) for s in db.values())

rows = []
for qname, hsps in hsps_per_query.items():
    m = len(queries[qname])
    for h in hsps:
        b = bit_score(h["score"])
        E = e_value(b, m, n_db_residues)
        rows.append({
            "query":     qname,
            "acc":       h["acc"],
            "raw_score": h["score"],
            "bits":      round(b, 1),
            "E":         f"{E:.2e}",
            "q_start":   h["q_start"],
            "t_start":   h["t_start"],
            "length":    h["length"],
        })

results = pd.DataFrame(rows)
print(f"Search space: m varies per query, n = {n_db_residues:,} residues")
print()
print(results.to_string(index=False))

# Sanity check 1: top hit for each query should be the query itself (when the
# query is in the human DB) or at least have a tiny E-value.
print("\\nBest E-value per query:")
for q, grp in results.groupby("query"):
    best = grp.iloc[grp["bits"].astype(float).idxmax() - grp.index[0]]
    print(f"  {q:>12s}  bits={best['bits']:7.1f}  E={best['E']}")
'''


STEP5_MD = """## Step 5 (11 min) — Compare to NCBI BLAST and visualise the cascade

NCBI BLAST is a black-box version of what you just built (with gapped
extension and composition correction layered on top). Open
https://blast.ncbi.nlm.nih.gov in another tab, paste **one** of your query
sequences into BLASTP against `nr`, and compare your top-5 hits against
theirs.

In the notebook we plot the **score histogram with the Karlin-Altschul
expected tail overlaid**: the chance HSP score distribution is
approximately Gumbel, and on a log-y axis its tail is a straight line.
Your real hits should sit far to the right of that tail.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Visualise the matched-filter cascade.
# ----------------------------------------------------------------------

# TODO:
# 1. Run a quick reseed against a synthetic ("random") database the same size
#    as `db` (use synthetic_db(n=len(db))). These are nulls — their best HSPs
#    define the empirical false-alarm distribution.
# 2. Plot a histogram of best HSP bit scores from the synthetic search and
#    overlay the Karlin-Altschul expected-tail curve y = m*n*2^(-x).
# 3. Mark the real-query top hits as vertical lines and verify they are far
#    to the right of the null tail.
'''


STEP5_SOLUTION = '''# Reference solution — Step 5.

# Build a *small* random database of the same total size for a quick null distribution.
print("Building null DB (synthetic) for empirical false-alarm comparison...")
null_db = synthetic_db(n=min(1000, len(db)), mean_len=350, seed=99)
null_index = build_db_index(null_db)

null_top_bits = []
# Use the shortest real query (insulin) — fastest, and still gives us a tail.
qname = "INS_HUMAN"
qseq = queries[qname]
m = len(qseq)
n_null = sum(len(s) for s in null_db.values())

null_hsps = blast_one_query(qseq, null_db, null_index, top_n=50)
for h in null_hsps:
    null_top_bits.append(bit_score(h["score"]))

real_top_bits = [bit_score(h["score"]) for h in hsps_per_query[qname][:10]]

# Plot.
fig, ax = plt.subplots(figsize=(8, 4.5))
bins = np.linspace(0, max(max(real_top_bits), max(null_top_bits)) + 5, 40)
ax.hist(null_top_bits, bins=bins, alpha=0.6, label=f"null DB ({len(null_db)} random seqs)")

# Karlin-Altschul expected-tail curve, plotted as a count using bin width.
bin_width = bins[1] - bins[0]
xs = np.linspace(bins[0], bins[-1], 200)
expected_count = m * n_null * (2.0 ** -xs) * bin_width / max(1, len(null_top_bits))
# Scale to the histogram by total null count.
ax.plot(xs, expected_count * len(null_top_bits) / max(1e-9, expected_count.max()),
        "r--", lw=2, label="Karlin-Altschul tail (scaled)")

for b in real_top_bits:
    ax.axvline(b, color="green", alpha=0.5, lw=1)
ax.axvline(real_top_bits[0], color="green", lw=2, label=f"real {qname} top hits")

ax.set_yscale("log")
ax.set_xlabel("bit score")
ax.set_ylabel("count (log)")
ax.set_title(f"Score distribution: random null vs real hits ({qname})")
ax.legend()
plt.tight_layout()
plt.show()

print()
print(f"Real {qname} top-10 bit scores: " + ", ".join(f"{b:.1f}" for b in real_top_bits))
print(f"Null top bit scores (min/median/max): "
      f"{min(null_top_bits):.1f} / {np.median(null_top_bits):.1f} / {max(null_top_bits):.1f}")
print()
print("Open https://blast.ncbi.nlm.nih.gov, paste INS_HUMAN into BLASTP vs nr,")
print("and compare your top-5 hits qualitatively. Bit-score and E-value ordering")
print("should agree, even though NCBI also does gapped extension + composition")
print("correction (we did not).")
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces of the pipeline. If
you ran the reference solutions above they should all pass; if you wrote
your own and an assert fails, revisit the corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. BLOSUM62 lookup matches the hand calculation.
assert score_pair("WKAA", "WRAA") == 21, "BLOSUM62 hand-score failed"

# 2. The bit-score formula matches the canonical value at a known raw score.
# With lambda=0.318, K=0.13: (0.318*50 - ln 0.13) / ln 2 ~= 25.88
b50 = bit_score(50)
assert 24.5 < b50 < 27.5, f"bit_score(50) out of expected band: {b50:.2f}"

# 3. E-value monotonically decreases with bit score.
assert e_value(50.0, 1000, 1_000_000) > e_value(70.0, 1000, 1_000_000)

# 4. Self-vs-self alignment for INS_HUMAN reaches a very high score.
# Build a tiny one-entry "database" containing just the query and re-search.
self_db = {"SELF": queries["INS_HUMAN"]}
self_index = build_db_index(self_db)
self_hsps = blast_one_query(queries["INS_HUMAN"], self_db, self_index, top_n=1)
assert self_hsps, "self-search returned no HSPs"
self_bits = bit_score(self_hsps[0]["score"])
assert self_bits > 80, f"self-vs-self bit score suspiciously low: {self_bits:.1f}"

# 5. Best E-value for INS_HUMAN against the full DB is tiny if the DB really is human
#    (else, with the synthetic fallback, the best HSP still beats E=1).
ins_top = hsps_per_query["INS_HUMAN"][0]
ins_bits = bit_score(ins_top["score"])
ins_E = e_value(ins_bits, len(queries["INS_HUMAN"]), sum(len(s) for s in db.values()))
print(f"INS_HUMAN top hit  acc={ins_top['acc']}  bits={ins_bits:.1f}  E={ins_E:.2e}")
assert ins_E < 1.0, f"INS_HUMAN top-hit E={ins_E:.2e} is uncomfortably large"

print("All self-checks passed.")
'''


EE_MD = """## EE framing — cascade detection, LLR test statistic, false-alarm rate

You implemented the canonical **matched-filter cascade**:

1. **Stage 1 — cheap prefilter.** The query 3-mer neighbourhood is a hash
   lookup. It throws away almost everything but is sensitive enough that the
   true matches keep at least one seed. Threshold `T` is the sensitivity /
   speed knob — the same knob you would tune in a CFAR radar detector.
2. **Stage 2 — expensive scoring.** X-drop ungapped extension is a 1-D
   matched-filter walk along the diagonal: accumulate BLOSUM62 column scores
   until the running sum drops `X` below its peak. BLOSUM62 entries are
   **log-odds**, so the running sum is a **log-likelihood ratio** comparing
   the target (matched) model to a random-protein null.
3. **Stage 3 — calibrate the statistic.** The Karlin-Altschul theory says
   that the maximum HSP score in a random search is approximately
   Gumbel-distributed; the bit score is the LLR rescaled into a unit where
   chance doubles per bit. The E-value is then the **false-alarm rate**
   scaled by the size of the search space `m·n`.

So when a BLAST result reports `E < 1e-50` it is saying: under the random
background, you would expect 10⁻⁵⁰ chance HSPs this good in a database this
big. That is the same false-alarm-rate language a radar engineer uses for a
P_FA on a constant-false-alarm-rate detector — the analogy isn't metaphor,
it is the math.
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
