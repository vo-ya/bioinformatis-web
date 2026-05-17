"""Build exercise.ipynb for L26 — Drug Discovery and Chemoinformatics.

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


TITLE_MD = """# L26 — Drug Discovery and Chemoinformatics

In this exercise you implement the chemoinformatics core that every
modern lead-discovery pipeline rests on: parse SMILES with RDKit,
compute Morgan (ECFP) fingerprints, build a Tanimoto similarity matrix
on 100 drug-like molecules, cluster them hierarchically, and apply
Lipinski + QED drug-likeness filters. The dataset comes from the
public ChEMBL REST API; if the network call fails the notebook drops
back to a curated 20-drug anchor set so the self-check still passes
offline.
"""


AHA_MD = """> **Aha.** A Morgan fingerprint is a **local feature extractor with a
> bounded receptive field** — each on-bit records the presence of a
> circular substructure of radius `r` around some atom. Tanimoto
> similarity is the **Jaccard index on this binary feature vector**,
> the same metric used by min-hash / locality-sensitive hashing. Once
> molecules live in this bit-vector space, hierarchical clustering is
> ordinary vector quantisation, and Lipinski's Rule of 5 is just a
> cheap rectangular filter on (MW, logP, HBD, HBA) — useful, but no
> substitute for actual assay data.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 matplotlib==3.8.4 requests==2.32.3 rdkit-pypi==2022.9.5 -q
"""


IMPORTS = """import io
import json
import math
import random
import urllib.request
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski
from rdkit.Chem.QED import qed

# Silence RDKit's chatty C++ warnings — we'll handle parse errors in Python.
RDLogger.DisableLog("rdApp.*")

# Deterministic for the whole notebook.
np.random.seed(42)
random.seed(42)

print("RDKit:", Chem.__name__, "imported successfully")
print("Test parse — aspirin:", Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O") is not None)
"""


STEP1_MD = """## Step 1 (8 min) — Fetch SMILES and compute Lipinski descriptors

We need 100 drug-like molecules. The cleanest source is the **ChEMBL
REST API**, filtered to approved drugs with MW between 150 and 500 Da.
That's the same molecular-weight window the Lipinski Rule of 5 cares
about, and it screens out salts / counter-ions that would crash the
RDKit parser. If the network call fails (sandboxed Colab, ChEMBL
maintenance, ...) we fall back to a hard-coded list of 20 anchor drugs
that everyone recognises — aspirin, ibuprofen, caffeine, metformin,
atorvastatin, dabigatran, ezetimibe, sertraline, paracetamol, and
friends. The pipeline then runs end-to-end on whichever set we get.

For each molecule you compute the four **Lipinski descriptors**:
molecular weight, octanol-water partition coefficient (logP),
hydrogen-bond donor count (HBD), and hydrogen-bond acceptor count
(HBA). These are the inputs to the Rule of 5 in step 4.
"""


STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Fetch SMILES (ChEMBL or fallback) and compute Lipinski descriptors.
# ----------------------------------------------------------------------

CHEMBL_URL = (
    "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
    "?molecule_properties__mw_freebase__range=150,500"
    "&max_phase=4&limit=100"
)

# Curated anchor-drug fallback (>= 20 SMILES of well-known approved drugs).
# Used if the ChEMBL REST call fails — the self-check must still pass.
ANCHOR_DRUGS = {
    "aspirin":         "CC(=O)Oc1ccccc1C(=O)O",
    "ibuprofen":       "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "caffeine":        "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "paracetamol":     "CC(=O)Nc1ccc(O)cc1",
    "metformin":       "CN(C)C(=N)NC(=N)N",
    "sertraline":      "CNC1CCc2ccc(Cl)c(Cl)c2C1c1ccccc1",
    "atorvastatin":    "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
    "dabigatran":      "CCCCCCOC(=O)N=C(N)c1ccc(NCc2nc3cc(C(=O)N(CCC(=O)OCC)c4ccccn4)ccc3n2C)cc1",
    "ezetimibe":       "OC(CCC1C(=O)N(c2ccc(F)cc2)C1c1ccc(O)cc1)c1ccc(F)cc1",
    "warfarin":        "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
    "diazepam":        "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "fluoxetine":      "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",
    "omeprazole":      "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",
    "simvastatin":     "CCC(C)(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12",
    "amoxicillin":     "CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O",
    "morphine":        "CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5",
    "codeine":         "COc1ccc2CC3C(C=CC(O)C34CCN(C)C24)O1",
    "nicotine":        "CN1CCCC1c1cccnc1",
    "acetaminophen":   "CC(=O)Nc1ccc(O)cc1",   # duplicate of paracetamol — both names common
    "ranitidine":      "CNC(=Cn1ccnc1C)NCCSCc1ccc(CN(C)C)o1",  # canonical-ish
    "loratadine":      "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1",
    "rosuvastatin":    "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1C=CC(O)CC(O)CC(=O)O",
    "ciprofloxacin":   "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
    "metoprolol":      "COCCc1ccc(OCC(O)CNC(C)C)cc1",
    "lisinopril":      "NCCCCC(NC(CCc1ccccc1)C(=O)N1CCCC1C(=O)O)C(=O)O",
    "tamoxifen":       "CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",
}


def fetch_chembl_smiles(n: int = 100) -> list[tuple[str, str]] | None:
    """Return [(chembl_id, smiles)] for ~n approved drug-like molecules."""
    # TODO: GET CHEMBL_URL, parse JSON, extract canonical_smiles per molecule,
    # filter out None / salts. Return None if the call fails.
    raise NotImplementedError


def parse_molecules(records: list[tuple[str, str]]) -> pd.DataFrame:
    """Parse each SMILES into an RDKit Mol; compute Lipinski descriptors."""
    # TODO: For each (id, smiles):
    #   - Chem.MolFromSmiles(smiles)
    #   - skip on None
    #   - compute MW (Descriptors.MolWt), logP (Crippen.MolLogP),
    #     HBD (Lipinski.NumHDonors), HBA (Lipinski.NumHAcceptors)
    # Return a DataFrame with columns: id, smiles, mol, MW, logP, HBD, HBA.
    raise NotImplementedError


# Set up the working dataframe. Try ChEMBL first; fall back to anchor drugs.
smiles_records = None  # list[(id, smiles)]
df = None              # DataFrame with descriptors
'''


SOLUTION_HEADER = """*Click ▶ to expand the reference solution.*"""


STEP1_SOLUTION = '''# Reference solution — Step 1.

def _http_get_json(url: str, timeout: float = 30.0):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "L26-exercise/1.0",
                                                   "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8", errors="replace"))
    except Exception as exc:
        print(f"  network error on {url[:60]}...: {exc}")
        return None


def fetch_chembl_smiles(n: int = 100):
    payload = _http_get_json(CHEMBL_URL)
    if payload is None:
        return None
    out = []
    for m in payload.get("molecules", []):
        chembl_id = m.get("molecule_chembl_id")
        structures = m.get("molecule_structures") or {}
        smi = structures.get("canonical_smiles")
        if chembl_id and smi:
            # Drop multi-fragment SMILES (salts) — keep largest fragment heuristic.
            frags = smi.split(".")
            smi = max(frags, key=len)
            out.append((chembl_id, smi))
        if len(out) >= n:
            break
    return out if out else None


def parse_molecules(records):
    rows = []
    for ident, smi in records:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        rows.append({
            "id":     ident,
            "smiles": Chem.MolToSmiles(mol),
            "mol":    mol,
            "MW":     Descriptors.MolWt(mol),
            "logP":   Crippen.MolLogP(mol),
            "HBD":    Lipinski.NumHDonors(mol),
            "HBA":    Lipinski.NumHAcceptors(mol),
        })
    return pd.DataFrame(rows)


print("Trying ChEMBL REST API...")
smiles_records = fetch_chembl_smiles(100)
if smiles_records is None or len(smiles_records) < 20:
    print("  ChEMBL fetch failed or too small; using anchor-drug fallback.")
    smiles_records = list(ANCHOR_DRUGS.items())

print(f"Got {len(smiles_records)} SMILES records; parsing with RDKit...")
df = parse_molecules(smiles_records)
print(f"Parsed {len(df)} valid molecules.")
print()
print(df.drop(columns=["mol"]).head(10).to_string(index=False))
print()
print("Descriptor summary:")
print(df[["MW", "logP", "HBD", "HBA"]].describe().round(2).to_string())
'''


STEP2_MD = """## Step 2 (12 min) — Morgan fingerprints and Tanimoto similarity matrix

A **Morgan fingerprint** (also called ECFP) assigns each atom an
initial integer hash from its element / charge / hybridisation, then
iteratively folds in the hashes of its neighbours out to radius `r`.
The set of hashes is then projected into a fixed-size bit vector
(2048 bits is the field default) by `hash % 2048`. Two molecules with
many shared substructures will share many on-bits.

**Tanimoto similarity** of two bit vectors `a` and `b` is the Jaccard
index:
$$T(a, b) = \\frac{|a \\cap b|}{|a \\cup b|} = \\frac{\\text{popcount}(a \\,\\&\\, b)}{\\text{popcount}(a \\,|\\, b)}$$
Range [0, 1]; 1.0 means identical fingerprints (not necessarily
identical molecules — collisions happen at 2048 bits).
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Morgan fingerprints + Tanimoto similarity matrix.
# ----------------------------------------------------------------------

FP_RADIUS = 2
FP_BITS   = 2048


def morgan_fp(mol, radius: int = FP_RADIUS, nbits: int = FP_BITS) -> np.ndarray:
    """Return a uint8 NumPy array of length nbits (0/1)."""
    # TODO: AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    # then convert to np.array of uint8 via DataStructs / ToBitString.
    raise NotImplementedError


def tanimoto_matrix(fps: np.ndarray) -> np.ndarray:
    """Pairwise Tanimoto similarity for an (N, B) uint8 matrix."""
    # TODO: vectorised intersection / union over bits.
    # |a & b| / |a | b| with no division-by-zero (default to 1 on the diagonal).
    raise NotImplementedError


# Compute fingerprints for the full set, assemble matrix.
'''


STEP2_SOLUTION = '''# Reference solution — Step 2.
from rdkit import DataStructs


def morgan_fp(mol, radius: int = FP_RADIUS, nbits: int = FP_BITS) -> np.ndarray:
    bv = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    arr = np.zeros((nbits,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(bv, arr)
    return arr


def tanimoto_matrix(fps: np.ndarray) -> np.ndarray:
    # fps: (N, B) uint8 / 0-1.
    fps = fps.astype(np.uint8)
    # popcount per row
    popcount = fps.sum(axis=1)
    # intersection via dot product (since values are 0/1)
    inter = fps @ fps.T
    # |a ∪ b| = |a| + |b| - |a ∩ b|
    union = popcount[:, None] + popcount[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        sim = np.where(union > 0, inter / np.maximum(union, 1), 1.0)
    # Diagonal should be exactly 1.0 by definition.
    np.fill_diagonal(sim, 1.0)
    return sim


fps = np.stack([morgan_fp(m) for m in df["mol"].tolist()])
print(f"Fingerprint matrix shape: {fps.shape}  (mean on-bit count = {fps.sum(axis=1).mean():.1f})")

sim = tanimoto_matrix(fps)
print(f"Similarity matrix shape: {sim.shape}")
print(f"Off-diagonal mean Tanimoto: {sim[np.triu_indices_from(sim, k=1)].mean():.3f}")
print(f"Off-diagonal max  Tanimoto: {sim[np.triu_indices_from(sim, k=1)].max():.3f}")

# Heatmap of similarity.
fig, ax = plt.subplots(figsize=(6, 5.2))
im = ax.imshow(sim, cmap="viridis", vmin=0.0, vmax=1.0, aspect="auto")
ax.set_title("Tanimoto similarity matrix (Morgan r=2, 2048 bits)")
ax.set_xlabel("molecule index")
ax.set_ylabel("molecule index")
fig.colorbar(im, ax=ax, label="Tanimoto")
plt.tight_layout()
plt.show()
'''


STEP3_MD = """## Step 3 (12 min) — Hierarchical clustering on 1 - Tanimoto

Tanimoto similarity has a natural distance form `d = 1 - T`. We feed
this distance matrix to `scipy.cluster.hierarchy.linkage` with the
**average-linkage** rule (every chemoinformatics paper's default) and
plot a dendrogram. The leaf order rearranges molecules so structural
neighbours sit next to each other — chemists call this a
**scaffold-clustered view**. Cutting the dendrogram at a fixed height
(say `d = 0.7`) yields discrete clusters; we'll colour the most
populous one and inspect a couple of its members.
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Hierarchical clustering on 1 - Tanimoto.
# ----------------------------------------------------------------------

LINKAGE_METHOD = "average"
CLUSTER_CUTOFF = 0.7  # in distance units = 1 - Tanimoto


def cluster_molecules(sim: np.ndarray,
                      method: str = LINKAGE_METHOD,
                      cutoff: float = CLUSTER_CUTOFF):
    """Return (Z, cluster_labels) — Z is the linkage matrix, labels are int IDs."""
    # TODO:
    # 1. d = 1 - sim, zero on the diagonal, then squareform(d) -> condensed.
    # 2. Z = linkage(condensed, method=method).
    # 3. labels = fcluster(Z, t=cutoff, criterion="distance").
    raise NotImplementedError


# Plot a dendrogram and report cluster sizes.
'''


STEP3_SOLUTION = '''# Reference solution — Step 3.

def cluster_molecules(sim: np.ndarray,
                      method: str = LINKAGE_METHOD,
                      cutoff: float = CLUSTER_CUTOFF):
    d = 1.0 - sim
    np.fill_diagonal(d, 0.0)
    # Numerical noise can push entries just below 0 — clip.
    d = np.clip(d, 0.0, 1.0)
    condensed = squareform(d, checks=False)
    Z = linkage(condensed, method=method)
    labels = fcluster(Z, t=cutoff, criterion="distance")
    return Z, labels


Z, labels = cluster_molecules(sim)
df["cluster"] = labels

cluster_sizes = Counter(labels)
top_clusters = cluster_sizes.most_common(5)
print(f"Found {len(cluster_sizes)} clusters at cutoff d = {CLUSTER_CUTOFF}")
print("Top 5 cluster sizes:", top_clusters)

# Dendrogram.
n = len(df)
fig, ax = plt.subplots(figsize=(10, 4))
ddata = dendrogram(
    Z,
    labels=[str(i) for i in range(n)],
    leaf_font_size=7,
    color_threshold=CLUSTER_CUTOFF,
    ax=ax,
)
ax.axhline(CLUSTER_CUTOFF, color="gray", lw=0.8, ls="--",
           label=f"cut at d = {CLUSTER_CUTOFF}")
ax.set_ylabel("distance (1 - Tanimoto)")
ax.set_title(f"Hierarchical clustering of {n} molecules ({LINKAGE_METHOD} linkage)")
ax.legend()
plt.tight_layout()
plt.show()

# Inspect the largest cluster.
if top_clusters:
    biggest_id, biggest_size = top_clusters[0]
    members = df[df["cluster"] == biggest_id].head(8)
    print(f"\\nLargest cluster (id={biggest_id}, size={biggest_size}) — first members:")
    print(members[["id", "smiles", "MW", "logP"]].to_string(index=False))
'''


STEP4_MD = """## Step 4 (15 min) — Lipinski Rule of 5 and QED

Christopher Lipinski's 1997 rule of thumb says oral drugs typically
satisfy:

- **MW** ≤ 500 Da
- **logP** ≤ 5
- **HBD** ≤ 5
- **HBA** ≤ 10

Violating any one of those gives one "Lipinski strike". A molecule
with ≥ 2 strikes is *probably* not a great oral drug, although several
approved drugs (atorvastatin, dabigatran, ezetimibe, ...) cheerfully
break the rule and work anyway.

A more graded score is **QED — Quantitative Estimate of Drug-likeness**
(Bickerton et al. 2012). QED combines 8 descriptors (MW, logP, HBD,
HBA, PSA, rotatable bonds, aromatic rings, structural alerts) into a
single number in [0, 1] using empirical desirability functions; higher
is more drug-like. RDKit ships an implementation.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Lipinski strikes + QED, scatter on (logP, MW).
# ----------------------------------------------------------------------

def lipinski_strikes(row) -> int:
    """Number of Rule-of-5 violations (0 = perfect, 4 = all four broken)."""
    # TODO: count violations of MW>500, logP>5, HBD>5, HBA>10.
    raise NotImplementedError


def add_qed(df: pd.DataFrame) -> pd.DataFrame:
    """Add a 'QED' column using rdkit.Chem.QED.qed."""
    # TODO
    raise NotImplementedError


# Scatter on (logP, MW); colour by Lipinski strikes; size by QED.
'''


STEP4_SOLUTION = '''# Reference solution — Step 4.

def lipinski_strikes(row) -> int:
    s = 0
    s += int(row["MW"]   > 500)
    s += int(row["logP"] > 5)
    s += int(row["HBD"]  > 5)
    s += int(row["HBA"]  > 10)
    return s


def add_qed(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["QED"] = [qed(m) for m in df["mol"].tolist()]
    return df


df["strikes"] = df.apply(lipinski_strikes, axis=1)
df = add_qed(df)

n_pass = (df["strikes"] == 0).sum()
n_fail2 = (df["strikes"] >= 2).sum()
print(f"Lipinski Rule of 5: {n_pass}/{len(df)} pass with 0 strikes; "
      f"{n_fail2} have >= 2 strikes")
print(f"QED — mean {df['QED'].mean():.3f}, median {df['QED'].median():.3f}, "
      f"min {df['QED'].min():.3f}, max {df['QED'].max():.3f}")

# Show worst Lipinski offenders, if any.
worst = df.sort_values("strikes", ascending=False).head(5)
print("\\nTop Lipinski violators:")
print(worst[["id", "smiles", "MW", "logP", "HBD", "HBA", "strikes", "QED"]].to_string(index=False))

# Scatter (logP, MW) coloured by strikes, sized by QED.
fig, ax = plt.subplots(figsize=(7.5, 5.5))
sc = ax.scatter(
    df["logP"], df["MW"],
    c=df["strikes"], cmap="RdYlGn_r", vmin=0, vmax=3,
    s=30 + 250 * df["QED"],
    alpha=0.75, edgecolors="k", linewidths=0.4,
)
ax.axhline(500, color="gray", lw=0.7, ls="--", label="MW = 500")
ax.axvline(5,   color="gray", lw=0.7, ls=":",  label="logP = 5")
ax.set_xlabel("logP (octanol-water)")
ax.set_ylabel("Molecular weight (Da)")
ax.set_title("Lipinski-strikes colour, QED size — Rule-of-5 boundaries dashed")
ax.legend(loc="upper left")
cbar = fig.colorbar(sc, ax=ax, label="Lipinski strikes")
cbar.set_ticks([0, 1, 2, 3])
plt.tight_layout()
plt.show()
'''


STEP5_MD = """## Step 5 (13 min) — Interpretation: where Lipinski fails

Some approved drugs **knowingly** break Lipinski's rule because their
target biology demanded it:

- **Dabigatran** (oral thrombin inhibitor): high HBD + HBA count and
  high MW because the molecule must mimic peptide hydrogen-bonding to
  fit thrombin's active site. The prodrug strategy (dabigatran
  etexilate) was the engineering response to the resulting poor oral
  absorption.
- **Ezetimibe** (cholesterol-absorption inhibitor): high logP / MW
  combination, again because it has to embed in the lipid layer of
  the intestinal brush border.
- **Atorvastatin** (statin): high MW + many H-bond acceptors; the
  formulation team had to do the heavy lifting that the chemistry
  could not.

Below we identify any anchor drugs in our dataset that fail Lipinski,
and quantify the **structure-similarity ↔ property-similarity**
relationship: do nearby molecules in fingerprint space have nearby
logP / MW / QED? This is the **molecular-similarity principle** —
the empirical observation that justifies fingerprint-based virtual
screening in the first place.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Activity-cliff sanity check + similarity-vs-property scatter.
# ----------------------------------------------------------------------

# TODO:
# 1. Print rows where strikes >= 2; comment on why these can still be drugs.
# 2. For every off-diagonal pair (i, j), compute pair Tanimoto sim
#    and pair |QED_i - QED_j|.  Scatter sim vs |dQED|; overlay a moving median
#    binned by sim. The trend should be: higher sim -> lower |dQED|.
# 3. Flag "activity-cliff candidates": pairs with sim > 0.5 but |dQED| > 0.3.
#    Print the count.
'''


STEP5_SOLUTION = '''# Reference solution — Step 5.

# 1. Lipinski violators in the dataset.
violators = df[df["strikes"] >= 2]
if len(violators):
    print(f"{len(violators)} molecules break >= 2 Lipinski rules:")
    print(violators[["id", "smiles", "MW", "logP", "HBD", "HBA", "strikes"]].to_string(index=False))
else:
    print("No molecules break >= 2 Lipinski rules in this dataset.")

# 2. Pairwise (similarity, |dQED|) scatter.
n = len(df)
iu = np.triu_indices(n, k=1)
pair_sim  = sim[iu]
pair_dqed = np.abs(df["QED"].values[iu[0]] - df["QED"].values[iu[1]])

# Moving-median trend by similarity bin.
bins = np.linspace(0, 1, 11)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
med = []
for lo, hi in zip(bins[:-1], bins[1:]):
    mask = (pair_sim >= lo) & (pair_sim < hi)
    med.append(np.median(pair_dqed[mask]) if mask.any() else np.nan)

fig, ax = plt.subplots(figsize=(7.5, 5))
ax.scatter(pair_sim, pair_dqed, s=6, alpha=0.25, label="pair")
ax.plot(bin_centers, med, "ro-", lw=2, label="median per sim bin")
ax.set_xlabel("Pairwise Tanimoto similarity")
ax.set_ylabel("|dQED|")
ax.set_title("Molecular-similarity principle: similar fingerprint -> similar QED")
ax.set_xlim(0, 1)
ax.set_ylim(0, max(0.55, pair_dqed.max() * 1.05))
ax.legend()
plt.tight_layout()
plt.show()

# Sanity check: lowest-similarity bins should have larger |dQED| than highest-similarity bins.
low_med  = np.nanmedian(pair_dqed[pair_sim < 0.2])
high_med = np.nanmedian(pair_dqed[pair_sim > 0.5]) if (pair_sim > 0.5).any() else low_med
print(f"\\nMedian |dQED| for low-similarity pairs  (sim < 0.2): {low_med:.3f}")
print(f"Median |dQED| for high-similarity pairs (sim > 0.5): {high_med:.3f}")

# 3. Activity-cliff candidates.
cliff_mask = (pair_sim > 0.5) & (pair_dqed > 0.3)
print(f"\\nActivity-cliff candidates (sim > 0.5 AND |dQED| > 0.3): {int(cliff_mask.sum())}")
if cliff_mask.any():
    i_arr, j_arr = iu[0][cliff_mask], iu[1][cliff_mask]
    print("Examples:")
    for i, j in zip(i_arr[:5], j_arr[:5]):
        print(f"  sim={sim[i, j]:.2f}  "
              f"{df.iloc[i]['id']} (QED={df.iloc[i]['QED']:.2f}) "
              f"vs {df.iloc[j]['id']} (QED={df.iloc[j]['QED']:.2f})")
'''


SELFCHECK_MD = """## Self-check

These asserts validate the load-bearing numerical pieces of the
pipeline. If you ran the reference solutions above they should all
pass; if you wrote your own and an assert fails, revisit the
corresponding step.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. The dataset has at least 20 valid molecules.
assert len(df) >= 20, f"only {len(df)} molecules after parsing — need >= 20"

# 2. Lipinski descriptors are in plausible ranges.
assert (df["MW"] > 0).all()
assert df["MW"].max() < 2000, "molecular weight implausibly large"
assert df["HBD"].min() >= 0 and df["HBA"].min() >= 0

# 3. Fingerprint matrix has the right shape and uses both 0 and 1.
assert fps.shape == (len(df), FP_BITS), f"fingerprint shape mismatch: {fps.shape}"
assert fps.sum() > 0 and fps.sum() < fps.size, "fingerprint matrix is constant"

# 4. Tanimoto matrix is symmetric, in [0, 1], with unit diagonal.
assert sim.shape == (len(df), len(df))
assert np.allclose(np.diag(sim), 1.0), "Tanimoto diagonal is not 1.0"
assert np.allclose(sim, sim.T, atol=1e-9), "Tanimoto matrix is not symmetric"
assert sim.min() >= -1e-9 and sim.max() <= 1.0 + 1e-9, "Tanimoto out of [0,1]"

# 5. A molecule must be perfectly similar to itself by Tanimoto.
i = 0
assert sim[i, i] == 1.0, "self-similarity != 1"

# 6. Clustering produced at least one cluster and labels for every row.
assert "cluster" in df.columns
assert df["cluster"].nunique() >= 1
assert len(df["cluster"]) == len(df)

# 7. QED is in [0, 1] for every molecule.
assert (df["QED"] >= 0).all() and (df["QED"] <= 1).all(), "QED out of [0, 1]"

# 8. The similarity principle holds: median |dQED| is lower in the high-similarity
#    bin than the low-similarity bin (when there's enough data).
iu = np.triu_indices(len(df), k=1)
pair_sim_chk  = sim[iu]
pair_dqed_chk = np.abs(df["QED"].values[iu[0]] - df["QED"].values[iu[1]])
low  = np.median(pair_dqed_chk[pair_sim_chk < 0.2]) if (pair_sim_chk < 0.2).any() else 1.0
high = np.median(pair_dqed_chk[pair_sim_chk > 0.5]) if (pair_sim_chk > 0.5).any() else 0.0
assert high <= low + 0.05, (
    f"similarity principle violated: high-sim median |dQED|={high:.3f} "
    f"vs low-sim median |dQED|={low:.3f}"
)

# 9. A hand-checked similarity sanity test: aspirin vs salicylic acid (parent
#    metabolite) should be highly similar. We compute this independently of
#    the dataset, so this assertion does not depend on which records ChEMBL
#    returned.
m_aspirin = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
m_salic   = Chem.MolFromSmiles("Oc1ccccc1C(=O)O")
fp_a = morgan_fp(m_aspirin)
fp_s = morgan_fp(m_salic)
t_as = (fp_a & fp_s).sum() / max(1, (fp_a | fp_s).sum())
assert t_as > 0.30, f"aspirin/salicylic-acid Tanimoto unexpectedly low: {t_as:.2f}"

print(f"\\nDataset: {len(df)} molecules; {df['cluster'].nunique()} clusters; "
      f"mean QED = {df['QED'].mean():.3f}")
print(f"Aspirin / salicylic acid Tanimoto = {t_as:.3f}")
print()
print("\\u2705 Self-check passed.")
'''


EE_MD = """## EE framing — local feature extraction, LSH, vector quantisation

You implemented a **bit-vector pipeline** with three classic
signal-processing analogues:

1. **Morgan fingerprint = local feature extractor.** Each on-bit
   records the presence of a circular substructure of radius `r`
   around some atom — a finite receptive field, just like a
   convolutional filter on an image. Folding 32-bit hashes into 2048
   bits is **feature hashing** (the same trick scikit-learn's
   `HashingVectorizer` uses for text).
2. **Tanimoto = Jaccard on binary feature vectors.** This is the
   distance metric that **MinHash / locality-sensitive hashing** is
   designed around — Tanimoto >= t can be approximated in
   sub-linear time on millions of fingerprints, which is how every
   real-world virtual screen at billion-molecule scale (Enamine REAL,
   ZINC22) actually works.
3. **Hierarchical clustering = vector quantisation.** Average-linkage
   on `1 - Tanimoto` is a standard agglomerative VQ scheme, and the
   dendrogram cut threshold sets the codebook size. **Lipinski + QED
   are the dimensionality-reduction step** from molecular graph down
   to a 4-D / 1-D summary that's amenable to threshold-based filtering.

So the chemoinformatics stack is signal processing in disguise:
extract local features, hash them into a bit-vector representation,
compare with Jaccard, cluster with VQ, then filter on cheap summary
descriptors. The same playbook works for documents, images, and
chemical structures.
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
