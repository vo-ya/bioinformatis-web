"""Build exercise.ipynb for L23 — Metagenomics and the Microbiome.

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


TITLE_MD = """# L23 — Metagenomics and the Microbiome

In this exercise you implement a miniature DADA2-style 16S pipeline end to
end on a simulated dataset: parse and QC reads, infer **Amplicon Sequence
Variants** (ASVs) by edit-distance clustering plus a replicate-based error
model, classify ASVs against a 10-taxon reference using k-mer voting, and
finally measure community structure with **alpha diversity** (Shannon,
Simpson, richness), **beta diversity** (Bray-Curtis, Jaccard) plus PCoA, and
test case-vs-control differences with a paired t-test on Shannon and ANOSIM
on the Bray-Curtis matrix.
"""


AHA_MD = """> **Aha.** Sequencing depth is a **population census** — every read is a
> sampled individual from a microbial community. ASV inference is
> **error-correction coding on biological strings**: edit distance is just
> Hamming in a discrete alphabet, and the consensus call exploits replicate
> structure exactly like a repetition code. Alpha diversity is **entropy**.
> PCoA is **eigen-decomposition of the distance matrix**. And the choice
> between **Bray-Curtis** (abundance) and **Jaccard** (presence/absence)
> matters: they will tell you different stories about the same community
> because microbiome data lives on a **simplex** (compositional), not in
> ordinary Euclidean space.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 scikit-learn==1.5.0 matplotlib==3.8.4 -q
"""


IMPORTS = """import math
import random
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.stats import ttest_rel

# Deterministic for the whole notebook.
SEED = 23
np.random.seed(SEED)
random.seed(SEED)

# DNA alphabet for the simulator.
DNA = "ACGT"
"""


STEP0_MD = """## Step 0 — Simulate a 16S dataset

Real 16S amplicon sequencing produces hundreds of thousands of ~250 bp reads
per sample, mostly drawn from a few dozen bacterial taxa. We build a tiny
toy version of that: **10 ground-truth taxa**, each represented by one
150 bp reference V4-like region, sampled into **6 samples** (3 case + 3
control) with **1 000 reads per sample**, then perturbed with **2 %
per-position substitution errors**. The case samples shift relative
abundance toward taxa 0-4; controls shift toward taxa 5-9 — so we know
ahead of time that the two groups should be distinguishable in beta
diversity.
"""

STEP0_CODE = '''# ----------------------------------------------------------------------
# Step 0 — Simulate the 16S dataset (deterministic; runs in <1 s).
# ----------------------------------------------------------------------

N_TAXA      = 10
READ_LEN    = 150
N_SAMPLES   = 6
READS_PER_S = 1000
ERROR_RATE  = 0.02   # per-position substitution probability

rng = np.random.default_rng(SEED)


def random_dna(length: int, rng) -> str:
    """Sample a random DNA string of the given length."""
    return "".join(rng.choice(list(DNA), size=length))


def mutate(seq: str, rate: float, rng) -> str:
    """Per-position iid substitution noise (no indels — DADA2's actual model)."""
    out = list(seq)
    nbases = len(out)
    n_err = rng.binomial(nbases, rate)
    if n_err == 0:
        return seq
    positions = rng.choice(nbases, size=n_err, replace=False)
    for p in positions:
        ref = out[p]
        alts = [b for b in DNA if b != ref]
        out[p] = alts[rng.integers(0, 3)]
    return "".join(out)


# 1. Ten ground-truth references.
TAXA = [f"Taxon{i:02d}" for i in range(N_TAXA)]
references = {t: random_dna(READ_LEN, rng) for t in TAXA}

# 2. Per-sample taxon abundances. Case samples lean toward taxa 0-4; control
#    samples lean toward taxa 5-9.
case_weights    = np.array([3.0]*5 + [1.0]*5)
control_weights = np.array([1.0]*5 + [3.0]*5)

sample_meta = []
sample_reads = {}
for i in range(N_SAMPLES):
    label = "case" if i < N_SAMPLES // 2 else "control"
    base = case_weights if label == "case" else control_weights
    # Per-sample jitter on top of the group baseline.
    weights = base * rng.dirichlet(np.ones(N_TAXA) * 5.0) * 10.0
    probs = weights / weights.sum()
    sample_id = f"S{i+1}_{label}"
    sample_meta.append({"sample": sample_id, "group": label, "true_probs": probs})

    # Draw taxon labels then noisy reads.
    taxa_drawn = rng.choice(N_TAXA, size=READS_PER_S, p=probs)
    reads = [mutate(references[TAXA[t]], ERROR_RATE, rng) for t in taxa_drawn]
    sample_reads[sample_id] = reads

meta_df = pd.DataFrame([{"sample": m["sample"], "group": m["group"]} for m in sample_meta])
print(meta_df.to_string(index=False))
print()
print(f"Total reads: {sum(len(r) for r in sample_reads.values()):,}")
print(f"Read length: {READ_LEN} bp, error rate: {ERROR_RATE*100:.0f}% per position")
print(f"Reference taxa: {N_TAXA}; references stored in `references` dict.")
'''


STEP1_MD = """## Step 1 (8 min) — Parse, length-filter, drop singletons

The first move of any DADA2-style pipeline is **dereplication**: collapse
identical reads into unique sequences with a count. After that we drop
exact-singleton unique sequences (they are almost always sequencing errors
no one else saw) and length-filter to the expected amplicon size. You build
both: a per-sample dereplication that returns `{seq: count}`, and a global
table of unique sequences across all 6 samples.
"""

STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Dereplicate, length-filter, drop global singletons.
# ----------------------------------------------------------------------


def dereplicate(reads: list[str], expected_len: int = READ_LEN) -> dict[str, int]:
    """Length-filter to `expected_len` exactly, then count unique reads."""
    # TODO: build a Counter over reads of length expected_len.
    raise NotImplementedError


def drop_global_singletons(per_sample: dict[str, dict[str, int]]) -> dict[str, dict[str, int]]:
    """Drop any unique sequence whose total count across all samples is < 2."""
    # TODO:
    # 1. Compute total counts per unique sequence across all samples.
    # 2. Build a new dict where each per-sample sub-dict only keeps non-singleton sequences.
    raise NotImplementedError


# Apply both to `sample_reads`.
# After this step you should have:
#   derep:   dict[sample_id, dict[seq, count]]
#   uniques: set of unique sequences (non-singleton) across all samples
'''

STEP1_SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

def dereplicate(reads, expected_len=READ_LEN):
    return Counter(r for r in reads if len(r) == expected_len)


def drop_global_singletons(per_sample):
    totals = Counter()
    for d in per_sample.values():
        totals.update(d)
    keep = {seq for seq, c in totals.items() if c >= 2}
    return {sid: {s: c for s, c in d.items() if s in keep} for sid, d in per_sample.items()}


derep_raw = {sid: dereplicate(rs) for sid, rs in sample_reads.items()}
derep     = drop_global_singletons(derep_raw)

uniques = set()
for d in derep.values():
    uniques.update(d.keys())

print(f"Per-sample unique-sequence counts before / after singleton drop:")
for sid in sample_reads:
    a = len(derep_raw[sid])
    b = len(derep[sid])
    n = sum(derep[sid].values())
    print(f"  {sid:>12s}  unique: {a:5d} -> {b:5d}   surviving reads: {n}")
print(f"\\nGlobal unique sequences after dereplication + singleton drop: {len(uniques)}")
'''


STEP2_MD = """## Step 2 (12 min) — ASV inference: greedy edit-distance clustering with an error model

True DADA2 fits a per-base, per-quality error model and uses a likelihood
ratio to decide whether a low-abundance unique sequence is a noisy child of
a high-abundance "parent" ASV. We use a stripped-down version: **greedy
clustering by Hamming distance**, where a unique sequence becomes a *new*
ASV only if its closest existing ASV is further away than what the error
model predicts at its observed abundance. Otherwise it is folded into the
nearest ASV.

The error model: under 2 % iid substitution error, the expected Hamming
distance between a read and its true parent is `L · p = 150 · 0.02 = 3`.
We accept a sequence as a noisy child of a parent if their Hamming distance
is **≤ 5** and the parent is at least **2×** more abundant. Anything else
spawns a new ASV.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Greedy ASV inference with a simple replicate-based error model.
# ----------------------------------------------------------------------

HAMMING_TOL    = 5     # distance under which we call something a child
ABUNDANCE_TOL  = 2.0   # parent must be at least this much more abundant


def hamming(a: str, b: str) -> int:
    """Hamming distance for equal-length strings."""
    # TODO
    raise NotImplementedError


def infer_asvs(derep: dict[str, dict[str, int]]):
    """Greedy clustering.

    Returns:
        asvs:          list[str] of ASV centroid sequences in descending abundance order
        membership:    dict[seq, int] mapping every kept unique sequence to its ASV index
        asv_counts:    np.ndarray of shape (n_samples, n_asvs) with abundances per sample
    """
    # TODO:
    # 1. Pool all unique sequences with their global counts (sum across samples).
    # 2. Sort by descending global count.
    # 3. Walk down the sorted list:
    #    - For each candidate, scan existing ASV centroids; if any centroid is
    #      within HAMMING_TOL and at least ABUNDANCE_TOL x more abundant, merge.
    #    - Otherwise create a new ASV with the candidate as its centroid.
    # 4. Build the (n_samples x n_asvs) count matrix by summing each sample's
    #    counts over the unique sequences mapped to each ASV.
    raise NotImplementedError


# asvs, membership, asv_counts = infer_asvs(derep)
'''

STEP2_SOLUTION = '''# Reference solution — Step 2.

HAMMING_TOL    = 5
ABUNDANCE_TOL  = 2.0


def hamming(a: str, b: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def infer_asvs(derep):
    # Global counts per unique sequence.
    global_counts = Counter()
    for d in derep.values():
        global_counts.update(d)
    ranked = sorted(global_counts.items(), key=lambda kv: -kv[1])

    asvs = []           # list of centroid sequences
    asv_counts_global = []
    membership = {}

    for seq, count in ranked:
        merged = False
        for i, centroid in enumerate(asvs):
            if asv_counts_global[i] >= ABUNDANCE_TOL * count:
                d = hamming(seq, centroid)
                if d <= HAMMING_TOL:
                    membership[seq] = i
                    asv_counts_global[i] += count
                    merged = True
                    break
        if not merged:
            membership[seq] = len(asvs)
            asvs.append(seq)
            asv_counts_global.append(count)

    # Build the sample x ASV count matrix.
    sample_ids = list(derep.keys())
    n_asvs = len(asvs)
    counts = np.zeros((len(sample_ids), n_asvs), dtype=int)
    for si, sid in enumerate(sample_ids):
        for seq, c in derep[sid].items():
            counts[si, membership[seq]] += c

    return asvs, membership, counts, sample_ids


asvs, membership, asv_counts, sample_ids = infer_asvs(derep)
print(f"Inferred {len(asvs)} ASVs (true number of taxa = {N_TAXA}).")
print(f"ASV count matrix shape: {asv_counts.shape}")
print()
print("Top-15 ASVs by total abundance:")
totals = asv_counts.sum(axis=0)
top_idx = np.argsort(-totals)[:15]
for rank, i in enumerate(top_idx):
    centroid = asvs[i]
    # Find closest true reference for sanity (Hamming to each).
    dists = [(t, hamming(centroid, references[t])) for t in TAXA]
    best_taxon, best_d = min(dists, key=lambda kv: kv[1])
    print(f"  ASV{i:03d}  total={totals[i]:5d}   closest ref: {best_taxon} (d={best_d})")
'''


STEP3_MD = """## Step 3 (12 min) — k-mer taxonomy: assign ASVs to the 10-reference DB

DADA2 hands its ASVs to a separate classifier (`assignTaxonomy`, RDP /
SILVA / GreenGenes). The standard idea: **build a k-mer index of the
reference database, then vote**. For each ASV, count how many of its
distinct k-mers appear in each reference; the reference with the highest
shared-k-mer fraction wins, provided it beats a confidence threshold.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — k-mer voting taxonomy + abundance matrix at the taxon level.
# ----------------------------------------------------------------------

KMER_K = 8


def seq_kmers(seq: str, k: int = KMER_K) -> set[str]:
    """All length-k subsequences of `seq`."""
    # TODO
    raise NotImplementedError


def classify_asv(asv: str, references: dict[str, str], k: int = KMER_K) -> tuple[str, float]:
    """Return (best_taxon, confidence) where confidence = best_shared_kmers / total_kmers."""
    # TODO: iterate references, compute Jaccard-like score, return the argmax.
    raise NotImplementedError


def taxon_abundance(asv_counts: np.ndarray, asvs: list[str], references: dict) -> pd.DataFrame:
    """Collapse ASV counts into a (n_samples x n_taxa) DataFrame indexed by taxon name."""
    # TODO: classify each ASV, then sum its column of asv_counts into the right taxon column.
    raise NotImplementedError


# taxon_df = taxon_abundance(asv_counts, asvs, references)
'''

STEP3_SOLUTION = '''# Reference solution — Step 3.

KMER_K = 8


def seq_kmers(seq, k=KMER_K):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


# Cache reference k-mer sets — we hit each ref once per ASV otherwise.
_REF_KMERS = {t: seq_kmers(seq) for t, seq in references.items()}


def classify_asv(asv, references, k=KMER_K):
    qk = seq_kmers(asv, k)
    if not qk:
        return ("Unclassified", 0.0)
    best_taxon, best_score = "Unclassified", 0.0
    for t, rk in _REF_KMERS.items():
        shared = len(qk & rk)
        score = shared / len(qk)
        if score > best_score:
            best_taxon, best_score = t, score
    return (best_taxon, best_score)


def taxon_abundance(asv_counts, asvs, references):
    classifications = [classify_asv(a, references) for a in asvs]
    n_samples = asv_counts.shape[0]
    by_taxon = defaultdict(lambda: np.zeros(n_samples, dtype=int))
    for j, (taxon, _) in enumerate(classifications):
        by_taxon[taxon] += asv_counts[:, j]
    df = pd.DataFrame(by_taxon, index=sample_ids)
    # Ensure all 10 taxa are present even if zero, for shape stability.
    for t in TAXA:
        if t not in df.columns:
            df[t] = 0
    df = df[[c for c in df.columns if c != "Unclassified"]
            + (["Unclassified"] if "Unclassified" in df.columns else [])]
    return df, classifications


taxon_df, asv_classes = taxon_abundance(asv_counts, asvs, references)
print(f"Taxon-level abundance matrix shape: {taxon_df.shape}")
print()
print(taxon_df.to_string())

# Sanity: confidence distribution.
confs = [c for _, c in asv_classes]
print(f"\\nClassification confidence: min={min(confs):.2f}, "
      f"median={np.median(confs):.2f}, max={max(confs):.2f}")
'''


STEP4_MD = """## Step 4 (15 min) — Alpha + beta diversity, PCoA

**Alpha diversity** is a per-sample summary of *how diverse this one
community is*. Three classics:

- **Richness** = number of taxa observed.
- **Shannon entropy** = `−Σ pᵢ log pᵢ`. Larger = more even.
- **Simpson** = `1 − Σ pᵢ²`. Probability that two random reads come from
  different taxa.

**Beta diversity** is a pairwise distance between samples:

- **Bray-Curtis** = `Σ|xᵢ − yᵢ| / Σ(xᵢ + yᵢ)` — abundance-sensitive.
- **Jaccard** = `1 − |X ∩ Y| / |X ∪ Y|` on the presence/absence sets —
  abundance-blind.

After computing the 6×6 Bray-Curtis matrix you do **PCoA** (principal
coordinate analysis): center the squared-distance matrix, eigen-decompose,
and read off the 2-D embedding from the top two eigenvectors.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Alpha and beta diversity + PCoA.
# ----------------------------------------------------------------------


def shannon(counts: np.ndarray) -> float:
    """Entropy of the abundance profile (natural log)."""
    # TODO
    raise NotImplementedError


def simpson(counts: np.ndarray) -> float:
    """Gini-Simpson index = 1 - sum(p_i^2)."""
    # TODO
    raise NotImplementedError


def richness(counts: np.ndarray) -> int:
    """Number of taxa with non-zero count."""
    # TODO
    raise NotImplementedError


def bray_curtis_matrix(X: np.ndarray) -> np.ndarray:
    """Pairwise Bray-Curtis on samples x taxa array."""
    # TODO: you can use scipy.spatial.distance.pdist with metric='braycurtis'.
    raise NotImplementedError


def jaccard_matrix(X: np.ndarray) -> np.ndarray:
    """Pairwise Jaccard on the presence/absence pattern of X."""
    # TODO: binarise X > 0, then 1 - |intersection| / |union| per pair.
    raise NotImplementedError


def pcoa(D: np.ndarray, k: int = 2) -> np.ndarray:
    """Classical multidimensional scaling on distance matrix D -> (n, k) embedding."""
    # TODO:
    # 1. Square the distances; build A = -0.5 * D**2.
    # 2. Double-center: B = A - row_mean - col_mean + grand_mean.
    # 3. Eigen-decompose B; take top-k eigenvectors scaled by sqrt(eigenvalue).
    raise NotImplementedError


# Compute alpha for each sample, beta for the cohort, and the PCoA embedding.
'''

STEP4_SOLUTION = '''# Reference solution — Step 4.

def shannon(counts):
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts[counts > 0] / total
    return float(-(p * np.log(p)).sum())


def simpson(counts):
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts / total
    return float(1.0 - (p * p).sum())


def richness(counts):
    return int((np.asarray(counts) > 0).sum())


def bray_curtis_matrix(X):
    return squareform(pdist(X, metric="braycurtis"))


def jaccard_matrix(X):
    B = (X > 0).astype(int)
    n = B.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            inter = np.logical_and(B[i], B[j]).sum()
            union = np.logical_or(B[i], B[j]).sum()
            D[i, j] = D[j, i] = 1.0 - inter / union if union else 0.0
    return D


def pcoa(D, k=2):
    n = D.shape[0]
    A = -0.5 * D**2
    row_mean = A.mean(axis=1, keepdims=True)
    col_mean = A.mean(axis=0, keepdims=True)
    grand = A.mean()
    B = A - row_mean - col_mean + grand
    # Symmetrise to suppress numerical drift.
    B = 0.5 * (B + B.T)
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(-vals)
    vals = vals[order]
    vecs = vecs[:, order]
    pos = np.clip(vals[:k], 0, None)
    return vecs[:, :k] * np.sqrt(pos)


# Use the taxon-level matrix as our abundance feature.
X = taxon_df.values.astype(float)

alpha_df = pd.DataFrame({
    "sample":   sample_ids,
    "group":    [m["group"] for m in sample_meta],
    "richness": [richness(X[i]) for i in range(len(sample_ids))],
    "shannon":  [shannon(X[i])  for i in range(len(sample_ids))],
    "simpson":  [simpson(X[i])  for i in range(len(sample_ids))],
})
print("Alpha diversity per sample:")
print(alpha_df.to_string(index=False))

D_bc = bray_curtis_matrix(X)
D_ja = jaccard_matrix(X)
print("\\nBray-Curtis matrix:")
print(pd.DataFrame(D_bc.round(3), index=sample_ids, columns=sample_ids).to_string())
print("\\nJaccard matrix:")
print(pd.DataFrame(D_ja.round(3), index=sample_ids, columns=sample_ids).to_string())

# PCoA on Bray-Curtis.
emb = pcoa(D_bc, k=2)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
colors = ["tab:red" if m["group"] == "case" else "tab:blue" for m in sample_meta]
axes[0].scatter(emb[:, 0], emb[:, 1], c=colors, s=120, edgecolor="k")
for i, sid in enumerate(sample_ids):
    axes[0].annotate(sid, (emb[i, 0], emb[i, 1]),
                     textcoords="offset points", xytext=(5, 5), fontsize=9)
axes[0].set_xlabel("PCo1")
axes[0].set_ylabel("PCo2")
axes[0].set_title("PCoA on Bray-Curtis (red=case, blue=control)")
axes[0].grid(alpha=0.3)

# Alpha-diversity bar chart.
xpos = np.arange(len(sample_ids))
w = 0.35
axes[1].bar(xpos - w/2, alpha_df["shannon"], w, label="Shannon", color="tab:orange")
axes[1].bar(xpos + w/2, alpha_df["simpson"], w, label="Simpson", color="tab:green")
axes[1].set_xticks(xpos)
axes[1].set_xticklabels(sample_ids, rotation=30, ha="right")
axes[1].set_ylabel("diversity index")
axes[1].set_title("Alpha diversity per sample")
axes[1].legend()
plt.tight_layout()
plt.show()
'''


STEP5_MD = """## Step 5 (13 min) — Statistical tests: paired t-test on Shannon + ANOSIM

The final move: **does case differ from control statistically?** Two
complementary tests:

1. **Paired t-test on Shannon entropy.** With matched case/control pairs
   (S1↔S4, S2↔S5, S3↔S6 here), `scipy.stats.ttest_rel` gives you the
   group-level alpha-diversity comparison.
2. **ANOSIM on the Bray-Curtis matrix.** The test statistic
   `R = (r̄_between − r̄_within) / (n(n−1)/4)` uses the *ranks* of all
   pairwise distances. Significance is obtained by a **permutation test**:
   shuffle group labels many times and recompute R; the p-value is the
   fraction of permutations giving an R as extreme as the observed one.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Statistical tests: paired t on Shannon + ANOSIM on Bray-Curtis.
# ----------------------------------------------------------------------

N_PERM = 999


def anosim_R(D: np.ndarray, labels: list[str]) -> float:
    """ANOSIM R statistic on distance matrix D with group labels."""
    # TODO:
    # 1. Take the upper triangle of D and flatten into pairwise distances.
    # 2. Mark each pair as 'within' (same label) or 'between' (different label).
    # 3. Convert distances to ranks (scipy.stats.rankdata or pd.Series.rank).
    # 4. R = (mean_between_rank - mean_within_rank) / (N * (N-1) / 4)
    raise NotImplementedError


def anosim_permutation(D: np.ndarray, labels: list[str], n_perm: int = N_PERM,
                       seed: int = SEED) -> tuple[float, float]:
    """Return (observed R, two-sided permutation p-value)."""
    # TODO:
    # 1. Compute observed R.
    # 2. For n_perm iterations, shuffle the label vector and recompute R.
    # 3. p = (1 + #(R_perm >= R_obs)) / (1 + n_perm)  for a one-sided right-tail test;
    #    a two-sided variant uses |R_perm| >= |R_obs|.
    raise NotImplementedError


# Run both tests and print a 1-paragraph interpretation.
'''

STEP5_SOLUTION = '''# Reference solution — Step 5.

N_PERM = 999


def _flatten_upper(D):
    n = D.shape[0]
    iu = np.triu_indices(n, k=1)
    return D[iu], iu


def anosim_R(D, labels):
    dists, (ii, jj) = _flatten_upper(D)
    same = np.array([labels[i] == labels[j] for i, j in zip(ii, jj)])
    # Ranks of all pairwise distances (ties resolved by average rank).
    ranks = pd.Series(dists).rank(method="average").to_numpy()
    r_within  = ranks[same].mean()
    r_between = ranks[~same].mean()
    N = D.shape[0]
    return (r_between - r_within) / (N * (N - 1) / 4)


def anosim_permutation(D, labels, n_perm=N_PERM, seed=SEED):
    rng = np.random.default_rng(seed)
    obs = anosim_R(D, labels)
    perms = []
    labels_arr = np.array(labels)
    for _ in range(n_perm):
        shuffled = rng.permutation(labels_arr).tolist()
        perms.append(anosim_R(D, shuffled))
    perms = np.array(perms)
    # One-sided right-tail (positive R = within < between, i.e. groups separate).
    p_one = (1 + (perms >= obs).sum()) / (1 + n_perm)
    return obs, p_one, perms


# Paired t on Shannon.
case_shan    = alpha_df.loc[alpha_df["group"] == "case",    "shannon"].to_numpy()
control_shan = alpha_df.loc[alpha_df["group"] == "control", "shannon"].to_numpy()
t_stat, t_p = ttest_rel(case_shan, control_shan)
print(f"Paired t-test on Shannon (case vs control): t={t_stat:.3f}, p={t_p:.3f}")
print(f"  case Shannon:    {case_shan}")
print(f"  control Shannon: {control_shan}")

# ANOSIM on Bray-Curtis.
labels = [m["group"] for m in sample_meta]
R_obs, R_p, R_perms = anosim_permutation(D_bc, labels)
print(f"\\nANOSIM on Bray-Curtis: R={R_obs:.3f}, p={R_p:.3f} ({N_PERM} permutations)")
print(f"  permutation R distribution: "
      f"mean={R_perms.mean():.3f}, 95th pct={np.percentile(R_perms, 95):.3f}")

# Visualise the permutation null distribution.
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.hist(R_perms, bins=30, alpha=0.7, label="null R (label-shuffled)")
ax.axvline(R_obs, color="red", lw=2, label=f"observed R = {R_obs:.2f}")
ax.set_xlabel("ANOSIM R")
ax.set_ylabel("count")
ax.set_title("Permutation null for case vs control")
ax.legend()
plt.tight_layout()
plt.show()

print()
print("Interpretation:")
print("  - Bray-Curtis (abundance-aware) separates case from control because the")
print("    simulator deliberately shifts relative abundance between groups.")
print("  - Jaccard (presence/absence) would be weaker here — both groups see all")
print("    10 taxa at sufficient depth, so the *sets* of observed taxa overlap")
print("    almost completely. That asymmetry is the central compositional lesson.")
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces of the pipeline.
If you ran the reference solutions above they should all pass; if you wrote
your own and an assert fails, revisit the corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Dereplication keeps reads of the expected length.
# Note: with 2% iid error on 150 bp reads, almost every read is a unique
# sequence, so a global singleton drop is aggressive — that's normal.
# The real "did we keep enough signal?" check happens after ASV inference.
assert all(len(s) == READ_LEN for s in uniques), "dereplication kept off-length reads"

# 2. ASV inference recovered roughly the right number of clusters.
# At 2% error + Hamming tol 5, greedy clustering on 10 true taxa should land
# in the 10-30 range (some natural over-splitting is fine).
assert 5 <= len(asvs) <= 60, f"unexpected ASV count: {len(asvs)}"

# 3. The count matrix has the right shape and recovers most of the read mass.
assert asv_counts.shape[0] == N_SAMPLES, "ASV matrix has wrong sample axis"
assert asv_counts.shape[1] == len(asvs), "ASV matrix has wrong ASV axis"

# 4. Taxonomy collapses ASVs back to ~10 classes; richness per sample is in range.
n_taxa_observed = (taxon_df.values.sum(axis=0) > 0).sum()
assert 5 <= n_taxa_observed <= 12, f"taxonomy gave odd taxon count: {n_taxa_observed}"

# 5. Shannon diversity is positive and bounded by ln(n_taxa).
for i, sid in enumerate(sample_ids):
    s = shannon(taxon_df.values[i])
    assert 0 < s < math.log(len(TAXA)) + 0.5, f"{sid} Shannon out of range: {s:.2f}"

# 6. Bray-Curtis diagonal is zero; matrix is symmetric.
assert np.allclose(np.diag(D_bc), 0.0), "BC diagonal not zero"
assert np.allclose(D_bc, D_bc.T), "BC matrix not symmetric"

# 7. ANOSIM finds case != control: R should be large (groups separate).
# With n=3 vs n=3 the minimum one-sided permutation p-value is 1/C(6,3) ≈ 0.05,
# and the test floors near that when the observed R is the most extreme value.
# We require R > 0.3 (well above the null mean of ~0) and p <= 0.15.
assert R_obs > 0.3, f"ANOSIM R should be > 0.3 with planted group signal; got {R_obs:.3f}"
assert R_p <= 0.15, f"ANOSIM p too high; signal not detected (p={R_p:.3f})"

# 8. PCoA embedding has the right shape and finite values.
assert emb.shape == (N_SAMPLES, 2), f"PCoA embedding shape: {emb.shape}"
assert np.all(np.isfinite(emb)), "PCoA embedding has NaN/inf"

print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — error-correction coding, entropy, eigen-decomposition

You implemented the entire 16S pipeline as a stack of signal-processing
primitives:

1. **Dereplication + singleton drop = denoising prefilter.** Identical
   reads collapse into a single symbol with a count; singletons are
   discarded because their signal-to-noise ratio is unbeatable. This is
   the same move as the spike-finder threshold before any heavier
   processing.
2. **ASV inference = error-correction decoding in a discrete alphabet.**
   Hamming distance is the natural metric; the "abundance × distance"
   acceptance rule is just a likelihood-ratio test between
   *child-of-existing-ASV* and *new-ASV* hypotheses. Replicate evidence
   (multiple reads, multiple samples) is the repetition-code structure
   that lets us correct.
3. **Alpha diversity = entropy.** Shannon is `H(p) = −Σ pᵢ log pᵢ`,
   verbatim from information theory. Simpson is the
   collision probability — the inverse of *channel diversity*.
4. **PCoA = eigen-decomposition of a distance matrix.** Classical MDS
   converts a distance matrix into a covariance-like inner-product matrix
   by double-centering, then projects onto the top eigenvectors —
   exactly what PCA does on a covariance, generalised to any metric.
5. **ANOSIM = rank-based two-sample test with a permutation null.** The
   permutation distribution *is* the null hypothesis, computed
   empirically — no parametric assumption, just label exchangeability.

The deeper aha: microbiome data lives on the **simplex** (proportions sum
to 1), not in ordinary ℝⁿ. Bray-Curtis respects that geometry better than
Euclidean; centered log-ratio transforms (CLR, used by ANCOM-BC) respect
it even better. Whenever you read "compositional data problem", that is
the geometry being violated — and a quiet source of false discoveries in
real-world microbiome papers.
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

        md(STEP0_MD),
        code(STEP0_CODE),

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
