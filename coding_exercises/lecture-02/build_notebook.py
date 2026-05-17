"""Build exercise.ipynb for L02 — Read Alignment: From Brute Force to FM Index and Back.

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


TITLE_MD = """# L02 — Read Alignment: From Brute Force to FM Index and Back

In this exercise you build the data structures behind every modern short-read
aligner: a **suffix array** over the *E. coli* K-12 reference, the
**Burrows-Wheeler transform** derived from it, and the **FM-index backward
search** that turns exact pattern matching into binary search. Then you
benchmark exact-match queries against the brute-force scan and watch the
scaling collapse from O(N) per query to O(L log N).
"""


AHA_MD = """> **Aha.** Preprocess once, query forever. Building the suffix array costs
> O(N log N) up front; afterwards every exact-match query is a sequence of
> rank lookups whose cost depends only on the pattern length `L`, **not** on
> the reference size `N`. The BWT is a data-dependent transform that
> concentrates like symbols into runs, making those rank lookups cheap.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 matplotlib==3.8.4 -q
"""


IMPORTS = """import random
import time
import urllib.request

import numpy as np
import matplotlib.pyplot as plt

# Deterministic for the whole notebook.
np.random.seed(42)
random.seed(42)

# DNA alphabet plus the BWT sentinel '$' (lexicographically smaller than ACGT).
ALPHABET = "$ACGT"
ALPHABET_IDX = {c: i for i, c in enumerate(ALPHABET)}
"""


STEP1_MD = """## Step 1 (8 min) — Load the reference and build a suffix array

We need a single ~4.6 Mb DNA string. The *E. coli* K-12 MG1655 reference
(GenBank `U00096.3`) is the standard teaching genome. We fetch it from the
NCBI EFetch REST endpoint; if the network call fails the notebook falls back
to a deterministic synthetic 200 kb sequence so the rest of the pipeline still
runs. We then build a **suffix array** — the permutation of indices `0..N-1`
that sorts the suffixes of `T + '$'` lexicographically.

The textbook construction is to sort `N` suffix pointers by their slice; that
is `O(N log² N)` and good enough for a teaching exercise on ~200 kb of text.
For the full 4.6 Mb *E. coli* genome we use NumPy's `argsort` on the integer
recoding of the text, which is fast enough on free Colab.
"""

STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Load the reference and build the suffix array.
# ----------------------------------------------------------------------

ECOLI_FASTA_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id=U00096.3&rettype=fasta&retmode=text"
)

def fetch_ecoli_reference(url: str = ECOLI_FASTA_URL) -> str | None:
    """Return the E. coli K-12 MG1655 reference as an uppercase ACGT string,
    or None on network failure."""
    # TODO: urlopen the URL, drop the FASTA header (lines starting with '>'),
    # concatenate the rest, uppercase, and keep only ACGT characters.
    raise NotImplementedError


def synthetic_reference(length: int = 200_000, seed: int = 0) -> str:
    """Deterministic random ACGT string for offline use."""
    # TODO: rng = np.random.default_rng(seed); sample ACGT with uniform freq.
    raise NotImplementedError


def build_suffix_array(text: str) -> np.ndarray:
    """Return the suffix array of text + '$' as an int32 numpy array of length len(text)+1.

    Hint: a clean approach is to sort range(N+1) by the suffix slice text[i:] + '$'.
    For longer text the slice-based key blows memory; recoding the text as a
    short-integer array and using np.argsort on increasing window lengths
    (prefix-doubling, optional) is faster.
    """
    # TODO: build and return the suffix-array permutation.
    raise NotImplementedError


# Placeholder values so downstream cells can run even with stubs.
T = ""          # the reference + sentinel will be assigned in the solution
SA = np.array([], dtype=np.int32)
'''


SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

ECOLI_FASTA_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id=U00096.3&rettype=fasta&retmode=text"
)


def fetch_ecoli_reference(url: str = ECOLI_FASTA_URL, timeout: float = 30.0):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "L02-exercise/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read().decode("utf-8", errors="replace")
    except Exception as exc:
        print(f"  network error: {exc}")
        return None
    # Drop the FASTA header line, keep only ACGT.
    seq_lines = [ln.strip() for ln in raw.splitlines() if ln and not ln.startswith(">")]
    seq = "".join(seq_lines).upper()
    seq = "".join(c for c in seq if c in "ACGT")
    return seq if len(seq) > 1_000_000 else None


def synthetic_reference(length: int = 200_000, seed: int = 0) -> str:
    rng = np.random.default_rng(seed)
    arr = rng.integers(0, 4, size=length, dtype=np.int8)
    return "".join("ACGT"[i] for i in arr)


def build_suffix_array(text: str) -> np.ndarray:
    """Suffix array of text + '$' using prefix-doubling on integer recoding.

    Runs in O(N log N) time using numpy sorts. Handles N up to a few million
    in a few seconds on free Colab.
    """
    s = text + "$"
    n = len(s)
    # Map each character to its rank (use ALPHABET ordering: $ < A < C < G < T).
    rank = np.array([ALPHABET_IDX.get(c, 0) for c in s], dtype=np.int64)
    sa = np.argsort(rank, kind="stable").astype(np.int64)
    tmp = np.empty(n, dtype=np.int64)

    k = 1
    while k < n:
        # Composite key: (rank[i], rank[i+k] if in-range else -1).
        right = np.full(n, -1, dtype=np.int64)
        idx = np.arange(n) + k
        valid = idx < n
        right[valid] = rank[idx[valid]]
        # Sort sa by (rank[sa], right[sa]).
        keys = rank * (n + 2) + (right + 1)
        sa = sa[np.argsort(keys[sa], kind="stable")]
        # Rebuild rank.
        tmp[sa[0]] = 0
        prev = keys[sa[0]]
        r = 0
        for i in range(1, n):
            cur = keys[sa[i]]
            if cur != prev:
                r += 1
                prev = cur
            tmp[sa[i]] = r
        rank = tmp.copy()
        if rank[sa[-1]] == n - 1:
            break
        k *= 2
    return sa.astype(np.int32)


print("Fetching E. coli K-12 MG1655 reference from NCBI EFetch...")
genome = fetch_ecoli_reference()
if genome is None:
    print("  EFetch unavailable; falling back to a 200 kb synthetic reference.")
    genome = synthetic_reference(length=200_000, seed=0)
else:
    # Keep the build under 5 min on free Colab — sub-sample the first ~500 kb
    # for the SA build. Backward search and brute-force scan see this same T.
    genome = genome[:500_000]

T = genome
N = len(T)
print(f"Reference length N = {N:,} bp")

t0 = time.time()
SA = build_suffix_array(T)
print(f"Suffix array built in {time.time() - t0:.2f}s, |SA| = {len(SA):,}")

# Sanity-check: the first suffix in sorted order ends with '$' alone.
text_plus = T + "$"
assert SA[0] == N, "first sorted suffix should be the lone sentinel"
# Pairwise check: adjacent suffixes are lexicographically increasing.
for i in range(1, 5):
    a, b = SA[i - 1], SA[i]
    assert text_plus[a:a + 20] <= text_plus[b:b + 20], f"SA not sorted at i={i}"
print("First 5 suffix-array entries:", SA[:5].tolist())
print("Their leading 12 chars:")
for i in SA[:5]:
    print(f"  SA[{int(i):>6d}] -> {text_plus[i:i+12]!r}")
'''


STEP2_MD = """## Step 2 (8 min) — Burrows-Wheeler transform and LF mapping

The **BWT** of `T$` is the last column `L` of the sorted-rotations matrix.
With the suffix array in hand it is trivial: `L[i] = T$[(SA[i] - 1) mod (N+1)]`.

For backward search we also need two precomputed arrays:

- `C[c]` = number of characters in `T$` that are **strictly smaller** than
  symbol `c` in alphabet order. `C[A]` is the count of `$`, `C[C]` is the
  count of `$` + `A`, and so on.
- `Occ[c, i]` = number of occurrences of symbol `c` in `L[:i]`. With these
  two arrays the **LF mapping** that walks the BWT backwards is
  `LF(i) = C[L[i]] + Occ[L[i], i]`.

You will also verify the inverse-BWT on a hand example (`"banana$"`) to make
sure your LF mapping is correct.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — BWT, C array, Occ table, LF mapping.
# ----------------------------------------------------------------------

def build_bwt(text: str, sa: np.ndarray) -> str:
    """Return BWT(T + '$') as a string of length N+1.

    L[i] = text_plus[(SA[i] - 1) mod (N+1)]   with text_plus = text + '$'.
    """
    # TODO
    raise NotImplementedError


def build_c_array(bwt: str) -> dict:
    """Return {c: count of symbols in bwt strictly less than c} for c in ALPHABET."""
    # TODO
    raise NotImplementedError


def build_occ(bwt: str) -> dict:
    """Return {c: int32 array of length N+2} so Occ[c, i] = count of c in bwt[:i]."""
    # TODO: a column-wise cumulative count works (one int32 array per alphabet symbol).
    raise NotImplementedError


def lf(c_arr: dict, occ: dict, bwt: str, i: int) -> int:
    """LF mapping: position i in BWT -> position of the same row in F column."""
    # TODO: return c_arr[bwt[i]] + occ[bwt[i]][i]
    raise NotImplementedError


def inverse_bwt(bwt: str) -> str:
    """Reconstruct text + '$' from its BWT using repeated LF walks from row 0."""
    # TODO
    raise NotImplementedError


# Verify on the classical 'banana' example.
# bwt('banana$') should be 'annb$aa'  (sorted rotations).
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.

def build_bwt(text: str, sa: np.ndarray) -> str:
    s = text + "$"
    n = len(s)
    out = []
    for i in sa:
        j = int(i) - 1
        if j < 0:
            j = n - 1
        out.append(s[j])
    return "".join(out)


def build_c_array(bwt: str) -> dict:
    # Use whichever symbols actually occur; '$' should always sort first.
    alphabet = sorted(set(bwt))
    counts = {c: 0 for c in alphabet}
    for ch in bwt:
        counts[ch] += 1
    # Cumulative: C[c] = sum of counts[c'] for c' < c.
    C = {}
    running = 0
    for c in alphabet:
        C[c] = running
        running += counts[c]
    return C


def build_occ(bwt: str) -> dict:
    n = len(bwt)
    alphabet = sorted(set(bwt))
    occ = {c: np.zeros(n + 1, dtype=np.int32) for c in alphabet}
    # Precompute per-position cumulative counts.
    arr = np.frombuffer(bwt.encode("ascii"), dtype=np.uint8)
    for c in alphabet:
        mask = (arr == ord(c)).astype(np.int32)
        occ[c][1:] = np.cumsum(mask)
    return occ


def lf(c_arr: dict, occ: dict, bwt: str, i: int) -> int:
    ch = bwt[i]
    return c_arr[ch] + int(occ[ch][i])


def inverse_bwt(bwt: str) -> str:
    c_arr = build_c_array(bwt)
    occ = build_occ(bwt)
    n = len(bwt)
    # Standard inversion: start at the row whose L column is '$' (so that
    # walking LF reads T+'$' character by character, ending with the first
    # char of T). Reverse at the end to get the natural order.
    i = bwt.index("$")
    out = []
    for _ in range(n):
        out.append(bwt[i])
        i = lf(c_arr, occ, bwt, i)
    return "".join(reversed(out))


# ----- Hand example: banana -----
def _ban_sa(text: str) -> list:
    s = text + "$"
    return sorted(range(len(s)), key=lambda k: s[k:])


banana = "banana"
ban_sa = np.array(_ban_sa(banana), dtype=np.int32)
ban_bwt = build_bwt(banana, ban_sa)
print(f"BWT('banana$') = {ban_bwt!r}  (expected 'annb$aa')")
assert ban_bwt == "annb$aa", f"banana BWT mismatch: {ban_bwt!r}"
recovered = inverse_bwt(ban_bwt)
print(f"inverse_bwt -> {recovered!r}")
assert recovered == "banana$", f"inverse BWT failed: {recovered!r}"

# ----- Build on the real reference -----
t0 = time.time()
BWT = build_bwt(T, SA)
C = build_c_array(BWT)
OCC = build_occ(BWT)
print(f"BWT + C + Occ built in {time.time() - t0:.2f}s")
print(f"C array: {{ {', '.join(f'{c!r}:{v:,}' for c, v in C.items())} }}")
print(f"BWT first 60 chars: {BWT[:60]!r}")
'''


STEP3_MD = """## Step 3 (12 min) — FM-index backward search

Backward search reads the pattern `P` **right to left**. Maintain an interval
`[lo, hi)` in the F column that contains exactly the rows whose suffix begins
with the pattern processed so far. For each preceding character `c`:

```
lo' = C[c] + Occ[c, lo]
hi' = C[c] + Occ[c, hi]
```

If `lo' >= hi'` the pattern has zero matches. Otherwise the surviving SA
entries `SA[lo':hi']` are the genomic positions of all exact matches. That is
`O(L)` rank lookups per query, with `L = len(P)`, independent of the
reference size.

You will also write a brute-force scanner for comparison — the textbook O(N·L)
naive matcher.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Backward search (FM-index) and brute-force scan.
# ----------------------------------------------------------------------

def backward_search(pattern: str, bwt: str, c_arr: dict, occ: dict) -> tuple[int, int]:
    """Return (lo, hi) — half-open interval in SA of rows matching `pattern`.

    Returns (0, 0) if the pattern does not occur.
    """
    # TODO: walk pattern right-to-left, updating lo and hi by the LF rule.
    raise NotImplementedError


def fm_occurrences(pattern: str, bwt: str, sa: np.ndarray,
                   c_arr: dict, occ: dict) -> np.ndarray:
    """Return all genomic positions where pattern occurs exactly."""
    # TODO: call backward_search; slice SA[lo:hi]; return as a sorted numpy array.
    raise NotImplementedError


def brute_force_occurrences(pattern: str, text: str) -> np.ndarray:
    """Naive O(N*L) scan over text for exact matches of pattern."""
    # TODO: slide a window across text; record start positions where text[i:i+L] == pattern.
    raise NotImplementedError


# Smoke test: pick a pattern that definitely occurs (a slice of T itself).
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.

def backward_search(pattern: str, bwt: str, c_arr: dict, occ: dict) -> tuple:
    if not pattern:
        return 0, len(bwt)
    lo, hi = 0, len(bwt)
    for ch in reversed(pattern):
        if ch not in c_arr:
            return 0, 0
        lo = c_arr[ch] + int(occ[ch][lo])
        hi = c_arr[ch] + int(occ[ch][hi])
        if lo >= hi:
            return 0, 0
    return lo, hi


def fm_occurrences(pattern: str, bwt: str, sa: np.ndarray,
                   c_arr: dict, occ: dict) -> np.ndarray:
    lo, hi = backward_search(pattern, bwt, c_arr, occ)
    if lo >= hi:
        return np.array([], dtype=np.int32)
    return np.sort(sa[lo:hi].astype(np.int32))


def brute_force_occurrences(pattern: str, text: str) -> np.ndarray:
    L = len(pattern)
    n = len(text)
    out = []
    for i in range(n - L + 1):
        if text[i:i + L] == pattern:
            out.append(i)
    return np.array(out, dtype=np.int32)


# Smoke test: a 50-mer pulled directly from the reference must be found.
test_pos = 10_000
test_pat = T[test_pos:test_pos + 50]
fm_hits = fm_occurrences(test_pat, BWT, SA, C, OCC)
bf_hits = brute_force_occurrences(test_pat, T)
print(f"Smoke test pattern at position {test_pos}, length {len(test_pat)}")
print(f"  FM hits:          {fm_hits.tolist()}")
print(f"  Brute-force hits: {bf_hits.tolist()}")
assert test_pos in fm_hits and test_pos in bf_hits
assert np.array_equal(fm_hits, bf_hits), "FM and brute-force disagree"
print("FM and brute-force agree.")

# Negative test: a pattern that is unlikely to occur.
ghost = "ZZZZZ"  # contains a non-ACGT character
ghost_hits = fm_occurrences(ghost, BWT, SA, C, OCC)
assert len(ghost_hits) == 0, "expected no hits for non-ACGT pattern"
print(f"Negative test ({ghost!r}): {len(ghost_hits)} hits — OK.")
'''


STEP4_MD = """## Step 4 (12 min) — 100 simulated reads: FM vs brute force

Now we measure the speedup. Sample 100 read-positions uniformly from the
reference; extract 50-bp slices as exact-match queries. Time both the FM
backward search and the brute-force scan over all 100 queries; report the
mean time per query and the speedup ratio. Verify that the two methods
return the same set of hits for every query.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Time 100 queries; FM vs brute force.
# ----------------------------------------------------------------------

NUM_READS = 100
READ_LEN  = 50

def sample_reads(text: str, n: int, L: int, seed: int = 7) -> list[tuple[int, str]]:
    """Return n (pos, slice) pairs sampled uniformly from text."""
    # TODO
    raise NotImplementedError


def time_queries(patterns, query_fn) -> tuple[float, list]:
    """Call query_fn(pattern) for each pattern; return (total_seconds, list_of_results)."""
    # TODO
    raise NotImplementedError


# Sample reads, time both methods, assert hit-set agreement, print summary.
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.

NUM_READS = 100
READ_LEN  = 50


def sample_reads(text: str, n: int, L: int, seed: int = 7) -> list:
    rng = np.random.default_rng(seed)
    starts = rng.integers(0, len(text) - L, size=n)
    return [(int(s), text[s:s + L]) for s in starts]


def time_queries(patterns, query_fn):
    t0 = time.time()
    results = [query_fn(p) for p in patterns]
    return time.time() - t0, results


reads = sample_reads(T, NUM_READS, READ_LEN, seed=7)
patterns = [p for _, p in reads]
true_positions = [pos for pos, _ in reads]

t_fm, fm_results = time_queries(patterns,
                                lambda p: fm_occurrences(p, BWT, SA, C, OCC))
t_bf, bf_results = time_queries(patterns,
                                lambda p: brute_force_occurrences(p, T))

# Hit-set agreement (sorted arrays already).
for i, (fm, bf) in enumerate(zip(fm_results, bf_results)):
    assert np.array_equal(fm, bf), f"hit-set disagreement on read {i}"
    assert true_positions[i] in fm, f"true position missing from FM hits at read {i}"

print(f"{NUM_READS} queries of length {READ_LEN}:")
print(f"  FM backward search:  total {t_fm*1000:8.1f} ms"
      f"   mean {t_fm/NUM_READS*1000:7.3f} ms/query")
print(f"  Brute-force scan:    total {t_bf*1000:8.1f} ms"
      f"   mean {t_bf/NUM_READS*1000:7.3f} ms/query")
print(f"  Speedup: {t_bf/t_fm:.0f}x")

# Histogram: number of occurrences per read (uniformly-sampled 50-mers in a
# 500 kb reference are mostly unique).
n_hits = np.array([len(r) for r in fm_results])
print(f"\\nHit-count summary across {NUM_READS} reads:")
print(f"  min/median/max occurrences: "
      f"{n_hits.min()} / {int(np.median(n_hits))} / {n_hits.max()}")
'''


STEP5_MD = """## Step 5 (20 min) — Scaling: query time vs pattern length

Vary the pattern length `L ∈ {20, 30, 50, 100, 200}` and time both methods.
The brute-force scan is roughly `O(N · L)` per query (every position requires
up to `L` comparisons); the FM backward search is `O(L)` rank lookups per
query, independent of `N`. We plot the measured mean-time-per-query against
`L` for both methods, with `L log N / log N₀` as a sanity reference for the
FM curve. The horizontal-ish FM line next to the linearly-growing brute-force
line is the punchline.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Scaling with pattern length.
# ----------------------------------------------------------------------

LENGTHS = [20, 30, 50, 100, 200]
N_PER_LEN = 30  # queries per length for a stable mean

def scan_lengths(text: str, lengths, n_per_len) -> dict:
    """Return {length: (t_fm_mean_s, t_bf_mean_s)} aggregated over n_per_len queries."""
    # TODO: for each L, sample n_per_len reads, time both methods, store means.
    raise NotImplementedError


# Build a log-log plot of mean time vs L, one curve per method.
'''


STEP5_SOLUTION = '''# Reference solution — Step 5.

LENGTHS = [20, 30, 50, 100, 200]
N_PER_LEN = 30


def scan_lengths(text: str, lengths, n_per_len) -> dict:
    results = {}
    for L in lengths:
        reads = sample_reads(text, n_per_len, L, seed=100 + L)
        pats = [p for _, p in reads]

        t_fm, fm_res = time_queries(pats,
                                    lambda p: fm_occurrences(p, BWT, SA, C, OCC))
        t_bf, bf_res = time_queries(pats,
                                    lambda p: brute_force_occurrences(p, T))
        # Sanity: methods agree.
        for fm, bf in zip(fm_res, bf_res):
            assert np.array_equal(fm, bf), f"disagreement at L={L}"
        results[L] = (t_fm / n_per_len, t_bf / n_per_len)
    return results


timings = scan_lengths(T, LENGTHS, N_PER_LEN)

# Tabulate.
print(f"{'L':>5s}  {'FM ms/query':>14s}  {'BF ms/query':>14s}  {'speedup':>8s}")
for L in LENGTHS:
    tf, tb = timings[L]
    print(f"{L:5d}  {tf*1000:14.4f}  {tb*1000:14.4f}  {tb/max(tf, 1e-9):8.1f}x")

# Plot.
Ls = np.array(LENGTHS, dtype=float)
fm_ms = np.array([timings[L][0] * 1000 for L in LENGTHS])
bf_ms = np.array([timings[L][1] * 1000 for L in LENGTHS])

# Theoretical reference: FM scales like L; BF scales like L (the cost is
# O((N - L + 1)) but the inner loop early-exits, so it is roughly constant
# in L for small L over a fixed N — we plot a flat reference at the mean BF).
fig, ax = plt.subplots(figsize=(7.5, 4.5))
ax.loglog(Ls, fm_ms, "o-", label="FM backward search")
ax.loglog(Ls, bf_ms, "s-", label="Brute-force scan")
# Theoretical O(L) reference for FM, anchored at L=50.
ref_anchor = fm_ms[LENGTHS.index(50)]
ax.loglog(Ls, ref_anchor * Ls / 50.0, "k--", alpha=0.4, label="O(L) reference")
ax.set_xlabel("pattern length L (bp)")
ax.set_ylabel("mean time per query (ms)")
ax.set_title(f"Exact-match query scaling — reference length N = {len(T):,} bp")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.show()

print()
print("Punchline:")
print("  Brute force scans N positions per query — runtime tracks N, not L.")
print("  FM backward search does L rank lookups per query — runtime tracks L, not N.")
print(f"  Net speedup at L=50: {bf_ms[LENGTHS.index(50)] / fm_ms[LENGTHS.index(50)]:.0f}x")
'''


SELFCHECK_MD = """## Self-check

These asserts pin down the load-bearing pieces of the pipeline. If your own
implementations replaced the reference solutions, an assert failure here
points at the step you need to revisit.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Suffix array is a permutation of [0, N+1).
assert sorted(SA.tolist()) == list(range(len(T) + 1)), "SA is not a permutation"

# 2. Adjacent suffixes are lexicographically increasing on a small window.
text_plus = T + "$"
for i in range(1, 20):
    a, b = int(SA[i - 1]), int(SA[i])
    assert text_plus[a:a + 30] <= text_plus[b:b + 30], f"SA sort broken at i={i}"

# 3. BWT has the right length and only canonical symbols.
assert len(BWT) == len(T) + 1
assert set(BWT) <= set(ALPHABET), f"BWT has unexpected symbols: {set(BWT) - set(ALPHABET)}"

# 4. Inverse BWT round-trips on the classical banana example.
assert inverse_bwt("annb$aa") == "banana$"

# 5. Hand example: backward search for 'ana' in 'banana' returns rows whose SA
#    entries are exactly {1, 3} (the two 'ana' starts).
_ban_t = "banana"
_ban_sa = np.array(sorted(range(len(_ban_t) + 1),
                          key=lambda k: (_ban_t + "$")[k:]), dtype=np.int32)
_ban_bwt = build_bwt(_ban_t, _ban_sa)
_ban_c = build_c_array(_ban_bwt)
_ban_occ = build_occ(_ban_bwt)
ana_hits = fm_occurrences("ana", _ban_bwt, _ban_sa, _ban_c, _ban_occ)
assert sorted(ana_hits.tolist()) == [1, 3], f"banana 'ana' hits wrong: {ana_hits}"

# 6. FM and brute-force agree on a fresh random read sampled from T.
rng = np.random.default_rng(2024)
spot = int(rng.integers(0, len(T) - 80))
probe = T[spot:spot + 80]
fm_h = fm_occurrences(probe, BWT, SA, C, OCC)
bf_h = brute_force_occurrences(probe, T)
assert np.array_equal(fm_h, bf_h), "FM vs brute-force disagree on random probe"
assert spot in fm_h, "true position missing from FM hits"

# 7. Scaling sanity: FM mean time at L=50 is strictly less than brute force.
tf_50, tb_50 = timings[50]
assert tf_50 < tb_50, f"FM ({tf_50*1000:.3f} ms) not faster than BF ({tb_50*1000:.3f} ms)"

print("All self-checks passed.")
'''


EE_MD = """## EE framing — successive approximation on a sorted list

The suffix array reduces exact-pattern matching to **binary search** in a
sorted list of suffixes; the BWT then collapses that binary search into a
sequence of constant-time rank lookups via the LF mapping. The whole thing is
a textbook example of trading **storage for query latency**:

- The brute-force scan is an O(N·L) **direct comparison** at every offset —
  the genomic analogue of correlating a template against every shift of a
  long signal.
- The suffix array is a **successive-approximation lookup** — like a SAR
  ADC, each comparison halves the candidate set; `log N` comparisons pin the
  pattern in the sorted suffix list.
- The BWT is a **data-dependent transform** that concentrates equal symbols
  into runs, dropping the entropy of the last column. Low-entropy strings
  compress, and they support fast cumulative-count queries — which is
  exactly what `Occ[c, i]` needs.
- Backward search is the dual of that compression: read the pattern right to
  left and watch the `[lo, hi)` interval **whiten** down to the matching
  rows. Each step is one cheap lookup, independent of N.

So when bwa or bowtie2 indexes a 3 Gb human reference once and serves
billion-read alignments forever, it is reaping the same engineering trade
you just measured: heavy preprocessing buys cheap, query-length-only matching
forever after.
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
        md(SOLUTION_HEADER),
        hidden(STEP1_SOLUTION),

        md(STEP2_MD),
        code(STEP2_TODO),
        md(SOLUTION_HEADER),
        hidden(STEP2_SOLUTION),

        md(STEP3_MD),
        code(STEP3_TODO),
        md(SOLUTION_HEADER),
        hidden(STEP3_SOLUTION),

        md(STEP4_MD),
        code(STEP4_TODO),
        md(SOLUTION_HEADER),
        hidden(STEP4_SOLUTION),

        md(STEP5_MD),
        code(STEP5_TODO),
        md(SOLUTION_HEADER),
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
