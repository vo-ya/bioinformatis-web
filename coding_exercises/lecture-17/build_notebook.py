"""Build exercise.ipynb for L17 — Clinical Genomics, Variant Interpretation, and Ethics.

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


TITLE_MD = """# L17 — Clinical Genomics, Variant Interpretation, and Ethics

In this exercise you implement the **ACMG/AMP 2015 variant-classification
rule engine**: a structured Boolean / lookup-table classifier that assembles
evidence codes (PVS1, PS1–4, PM1–6, PP1–5, BA1, BS1–4, BP1–7) into one of
five class labels (Pathogenic / Likely Pathogenic / VUS / Likely Benign /
Benign). You then expose the engine's two big practical weaknesses by
running an **ancestry-bias audit** against gnomAD-style allele frequencies
and a **polygenic-risk-score (PRS) portability** experiment across two
ancestry backgrounds.
"""


AHA_MD = """> **Aha.** ACMG/AMP is **auditable evidence assembly**, not a black-box
> algorithm — every call traces back to specific evidence codes that map
> through a fixed combining table (Richards et al. 2015). Underneath the
> rule table is a Bayesian likelihood-ratio multiplier (Tavtigian et al.
> 2018). The whole framework is silently calibrated on European reference
> data; the same allele can flip from Likely Pathogenic in one ancestry to
> Benign in another simply because population frequency is the input the
> rule table cares about most. Polygenic risk scores inherit the same bias
> as **covariate shift**.
"""


PREAMBLE = """# Install pinned scientific stack on first run. Quiet so the notebook stays tidy.
!pip install numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 matplotlib==3.8.4 requests==2.32.3 -q
"""


IMPORTS = """import io
import json
import math
import urllib.request
from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import requests

# Deterministic for the whole notebook.
SEED = 42
np.random.seed(SEED)

# Five-class ACMG/AMP labels.
P  = "Pathogenic"
LP = "Likely Pathogenic"
VUS = "VUS"
LB = "Likely Benign"
B  = "Benign"

print("Imports loaded; notebook seed =", SEED)
"""


STEP1_MD = """## Step 1 (8 min) — Load 25 anchor variants

We need a reference set of clinically curated variants with the evidence
fields the ACMG/AMP rule engine cares about: gene, consequence type
(missense / frameshift / nonsense / splice), gnomAD allele frequencies
(stratified by ancestry), REVEL score, SpliceAI score, and the **prior**
ClinVar classification that we will measure ourselves against.

We try the public ClinVar + gnomAD REST APIs first. If the network is
unavailable (Colab firewalls, rate limits) we fall back to a deterministic,
in-notebook synthetic anchor set of 25 variants spanning BRCA1, BRCA2,
TP53, CFTR, MLH1, MYH7 and other classic clinical-genetics targets, with
ancestry-stratified allele frequencies hardcoded so that the rest of the
exercise runs identically on either path. The self-check at the bottom
passes on the synthetic path.
"""

STEP1_TODO = '''# ----------------------------------------------------------------------
# Step 1 — Load the 25-variant anchor set.
# ----------------------------------------------------------------------

# Each variant is a dict with this schema (field names matter — the rule
# engine in Step 2/3 keys off them):
#
#   id            stable identifier ("BRCA1:c.5266dupC")
#   gene          HUGO gene symbol
#   consequence   one of {missense, nonsense, frameshift, splice, synonymous}
#   gene_is_lof_disease  bool — is loss-of-function a known disease mechanism?
#   af_global     gnomAD popmax / global allele frequency (float 0..1)
#   af_eur        gnomAD non-Finnish European allele frequency
#   af_afr        gnomAD African / African-American allele frequency
#   af_eas        gnomAD East Asian allele frequency
#   revel         REVEL missense pathogenicity score (0..1) or None
#   spliceai      SpliceAI max delta score (0..1) or None
#   in_hotspot    bool — variant lies in a known mutational hotspot / functional domain
#   same_aa_as_path  bool — different missense change at the residue is already P/LP (PM5)
#   clinvar       prior ClinVar 5-class label (our "ground truth")
#
# TODO: implement load_anchor_variants() — try the public REST APIs, fall
# back to anchor_variants_synthetic() (already provided below) when the
# network is unreachable.

def anchor_variants_synthetic() -> list[dict]:
    """Deterministic 25-variant anchor set with ancestry-stratified AFs."""
    # See the hidden reference solution for the full table.
    raise NotImplementedError


def load_anchor_variants() -> list[dict]:
    """Try ClinVar/gnomAD REST; fall back to the synthetic anchor set."""
    raise NotImplementedError


# variants = load_anchor_variants()
# df = pd.DataFrame(variants)
'''


STEP1_SOLUTION_HEADER = """*Click the cell below to expand the reference solution.*"""

STEP1_SOLUTION = '''# Reference solution — Step 1.

def anchor_variants_synthetic() -> list[dict]:
    """Deterministic anchor set: 25 BRCA1/BRCA2/TP53/CFTR/MLH1/MYH7-style variants.

    Allele frequencies and predictor scores are realistic, ancestry-stratified
    placeholders chosen to exercise every ACMG/AMP rule branch at least once.
    No network required.
    """
    V = [
        # ---- High-confidence pathogenic LoF anchors ----
        dict(id="BRCA1:c.5266dupC",   gene="BRCA1",  consequence="frameshift",
             gene_is_lof_disease=True,
             af_global=8e-5, af_eur=2e-4, af_afr=1e-5, af_eas=0.0,
             revel=None, spliceai=None,
             in_hotspot=True, same_aa_as_path=False, clinvar="Pathogenic"),
        dict(id="BRCA1:c.68_69delAG", gene="BRCA1",  consequence="frameshift",
             gene_is_lof_disease=True,
             af_global=4e-5, af_eur=6e-5, af_afr=1e-5, af_eas=0.0,
             revel=None, spliceai=None,
             in_hotspot=True, same_aa_as_path=False, clinvar="Pathogenic"),
        dict(id="BRCA2:c.5946delT",   gene="BRCA2",  consequence="frameshift",
             gene_is_lof_disease=True,
             af_global=2e-5, af_eur=1e-5, af_afr=3e-5, af_eas=0.0,
             revel=None, spliceai=None,
             in_hotspot=True, same_aa_as_path=False, clinvar="Pathogenic"),
        dict(id="CFTR:p.F508del",     gene="CFTR",   consequence="frameshift",
             gene_is_lof_disease=True,
             af_global=7e-3, af_eur=1.4e-2, af_afr=2e-4, af_eas=5e-5,
             revel=None, spliceai=None,
             in_hotspot=True, same_aa_as_path=False, clinvar="Pathogenic"),
        dict(id="MLH1:c.350C>T",      gene="MLH1",   consequence="nonsense",
             gene_is_lof_disease=True,
             af_global=1e-5, af_eur=2e-5, af_afr=0.0, af_eas=0.0,
             revel=None, spliceai=None,
             in_hotspot=True, same_aa_as_path=False, clinvar="Pathogenic"),
        dict(id="BRCA1:c.5074+1G>A",  gene="BRCA1",  consequence="splice",
             gene_is_lof_disease=True,
             af_global=0.0, af_eur=0.0, af_afr=0.0, af_eas=0.0,
             revel=None, spliceai=0.97,
             in_hotspot=False, same_aa_as_path=False, clinvar="Pathogenic"),

        # ---- Pathogenic missense in TP53 hotspot ----
        dict(id="TP53:p.R175H",       gene="TP53",   consequence="missense",
             gene_is_lof_disease=True,
             af_global=0.0, af_eur=0.0, af_afr=0.0, af_eas=0.0,
             revel=0.95, spliceai=0.02,
             in_hotspot=True, same_aa_as_path=True, clinvar="Pathogenic"),
        dict(id="TP53:p.R248Q",       gene="TP53",   consequence="missense",
             gene_is_lof_disease=True,
             af_global=4e-6, af_eur=8e-6, af_afr=0.0, af_eas=0.0,
             revel=0.93, spliceai=0.01,
             in_hotspot=True, same_aa_as_path=True, clinvar="Pathogenic"),
        dict(id="TP53:p.R273H",       gene="TP53",   consequence="missense",
             gene_is_lof_disease=True,
             af_global=0.0, af_eur=0.0, af_afr=0.0, af_eas=0.0,
             revel=0.94, spliceai=0.02,
             in_hotspot=True, same_aa_as_path=True, clinvar="Pathogenic"),

        # ---- Likely-pathogenic missense (predictors only) ----
        dict(id="MLH1:p.A586P",       gene="MLH1",   consequence="missense",
             gene_is_lof_disease=True,
             af_global=2e-6, af_eur=4e-6, af_afr=0.0, af_eas=0.0,
             revel=0.82, spliceai=0.03,
             in_hotspot=True, same_aa_as_path=False, clinvar="Likely Pathogenic"),
        dict(id="MYH7:p.R453C",       gene="MYH7",   consequence="missense",
             gene_is_lof_disease=False,
             af_global=1e-5, af_eur=2e-5, af_afr=0.0, af_eas=0.0,
             revel=0.88, spliceai=0.04,
             in_hotspot=True, same_aa_as_path=True, clinvar="Likely Pathogenic"),
        dict(id="BRCA2:p.D2723H",     gene="BRCA2",  consequence="missense",
             gene_is_lof_disease=True,
             af_global=3e-6, af_eur=6e-6, af_afr=0.0, af_eas=0.0,
             revel=0.79, spliceai=0.05,
             in_hotspot=True, same_aa_as_path=False, clinvar="Likely Pathogenic"),

        # ---- VUS: predictor-ambiguous, in domain, no frequency evidence ----
        dict(id="BRCA1:p.K1487R",     gene="BRCA1",  consequence="missense",
             gene_is_lof_disease=True,
             af_global=5e-6, af_eur=1e-5, af_afr=0.0, af_eas=0.0,
             revel=0.42, spliceai=0.08,
             in_hotspot=False, same_aa_as_path=False, clinvar="VUS"),
        dict(id="CFTR:p.M470V",       gene="CFTR",   consequence="missense",
             gene_is_lof_disease=True,
             af_global=0.35, af_eur=0.50, af_afr=0.25, af_eas=0.40,
             revel=0.21, spliceai=0.03,
             in_hotspot=False, same_aa_as_path=False, clinvar="Benign"),
        dict(id="MLH1:p.V716M",       gene="MLH1",   consequence="missense",
             gene_is_lof_disease=True,
             af_global=4e-5, af_eur=6e-5, af_afr=1e-5, af_eas=2e-5,
             revel=0.50, spliceai=0.06,
             in_hotspot=False, same_aa_as_path=False, clinvar="VUS"),
        dict(id="MYH7:p.E1455K",      gene="MYH7",   consequence="missense",
             gene_is_lof_disease=False,
             af_global=2e-5, af_eur=4e-5, af_afr=0.0, af_eas=0.0,
             revel=0.55, spliceai=0.04,
             in_hotspot=False, same_aa_as_path=False, clinvar="VUS"),

        # ---- Likely benign / benign anchors ----
        dict(id="BRCA1:p.S1613G",     gene="BRCA1",  consequence="missense",
             gene_is_lof_disease=True,
             af_global=0.32, af_eur=0.30, af_afr=0.45, af_eas=0.05,
             revel=0.10, spliceai=0.01,
             in_hotspot=False, same_aa_as_path=False, clinvar="Benign"),
        dict(id="BRCA2:p.N372H",      gene="BRCA2",  consequence="missense",
             gene_is_lof_disease=True,
             af_global=0.24, af_eur=0.27, af_afr=0.13, af_eas=0.05,
             revel=0.14, spliceai=0.02,
             in_hotspot=False, same_aa_as_path=False, clinvar="Benign"),
        dict(id="CFTR:p.R75Q",        gene="CFTR",   consequence="missense",
             gene_is_lof_disease=True,
             af_global=0.012, af_eur=0.018, af_afr=0.004, af_eas=0.008,
             revel=0.12, spliceai=0.02,
             in_hotspot=False, same_aa_as_path=False, clinvar="Likely Benign"),
        dict(id="MLH1:p.I219V",       gene="MLH1",   consequence="missense",
             gene_is_lof_disease=True,
             af_global=0.36, af_eur=0.40, af_afr=0.20, af_eas=0.42,
             revel=0.08, spliceai=0.01,
             in_hotspot=False, same_aa_as_path=False, clinvar="Benign"),
        dict(id="MYH7:p.E1014K",      gene="MYH7",   consequence="missense",
             gene_is_lof_disease=False,
             af_global=0.01, af_eur=0.018, af_afr=0.002, af_eas=0.005,
             revel=0.18, spliceai=0.02,
             in_hotspot=False, same_aa_as_path=False, clinvar="Likely Benign"),

        # ---- The ancestry-flip cases — rare in EUR, common in AFR ----
        # These are the four anchors that demonstrate the bias hotspot in Step 4.
        dict(id="GENE_X:p.A123T",     gene="GENE_X", consequence="missense",
             gene_is_lof_disease=False,
             af_global=0.020, af_eur=0.0001, af_afr=0.080, af_eas=0.002,
             revel=0.65, spliceai=0.04,
             in_hotspot=False, same_aa_as_path=False, clinvar="Likely Benign"),
        dict(id="GENE_Y:p.R456H",     gene="GENE_Y", consequence="missense",
             gene_is_lof_disease=False,
             af_global=0.015, af_eur=0.00005, af_afr=0.060, af_eas=0.001,
             revel=0.58, spliceai=0.03,
             in_hotspot=False, same_aa_as_path=False, clinvar="Likely Benign"),
        dict(id="GENE_Z:p.G789D",     gene="GENE_Z", consequence="missense",
             gene_is_lof_disease=False,
             af_global=0.022, af_eur=0.00008, af_afr=0.090, af_eas=0.003,
             revel=0.51, spliceai=0.05,
             in_hotspot=False, same_aa_as_path=False, clinvar="Likely Benign"),
        dict(id="GENE_W:p.L321P",     gene="GENE_W", consequence="missense",
             gene_is_lof_disease=False,
             af_global=0.018, af_eur=0.00010, af_afr=0.072, af_eas=0.002,
             revel=0.62, spliceai=0.04,
             in_hotspot=False, same_aa_as_path=False, clinvar="Likely Benign"),
    ]
    assert len(V) == 25, f"expected 25 anchor variants, got {len(V)}"
    return V


def _http_get_json(url: str, timeout: float = 12.0):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "L17-exercise/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8", errors="replace"))
    except Exception as exc:
        print(f"  network error on {url[:80]}...: {exc}")
        return None


def _try_fetch_clinvar_subset() -> Optional[list[dict]]:
    """Best-effort ClinVar REST probe; returns None on any network issue.

    The exercise does not depend on this fetch — it only confirms the
    network path exists. We sample a single well-known accession just to
    show the API call shape; full ClinVar parsing is out of scope.
    """
    probe = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=17661&retmode=json"
    payload = _http_get_json(probe)
    if payload is None:
        return None
    print("  ClinVar REST reachable; using deterministic synthetic anchor set "
          "(parsing the full payload is out of scope for this exercise).")
    return None  # We deliberately return None and use the synthetic anchors.


def load_anchor_variants() -> list[dict]:
    _try_fetch_clinvar_subset()
    return anchor_variants_synthetic()


variants = load_anchor_variants()
df = pd.DataFrame(variants)
print(f"Loaded {len(df)} anchor variants across {df['gene'].nunique()} genes")
print(f"ClinVar label distribution:\\n{df['clinvar'].value_counts().to_string()}")
df.head(8)
'''


STEP2_MD = """## Step 2 (12 min) — Walk the ACMG/AMP rule engine on two anchors

The 2015 ACMG/AMP framework defines 28 evidence codes (PVS1, PS1–4,
PM1–6, PP1–5, BA1, BS1–4, BP1–7) and a fixed combining table:

- **Pathogenic** ⇐ `(PVS1 + ≥1 strong)` or `(PVS1 + ≥2 moderate)` or
  `(≥2 strong)` or `(1 strong + ≥3 moderate)` ...
- **Likely Pathogenic** ⇐ `(PVS1 + 1 moderate)` or `(≥1 strong + 1–2 moderate)`
  or `(≥3 moderate)` ...
- **Benign** ⇐ `(BA1)` or `(≥2 strong benign)`
- **Likely Benign** ⇐ `(1 strong + 1 supporting benign)` or `(≥2 supporting benign)`
- otherwise **VUS**.

You will assign codes to **two** anchor variants by hand and walk the
combining table to confirm your call matches ClinVar.
"""

STEP2_TODO = '''# ----------------------------------------------------------------------
# Step 2 — Assign codes for two anchors and walk the combining table.
# ----------------------------------------------------------------------

# Pick two anchors that exercise different rule paths:
#   - BRCA1:c.5266dupC      — frameshift in a LoF-disease gene + absent gnomAD
#   - BRCA1:p.S1613G        — common missense in BRCA1 (BA1 territory)

# Each evidence code has a strength tier; we expect a list of
# (code, strength) tuples like [("PVS1", "very_strong"), ("PM2", "moderate")].

def assign_codes_for_brca1_5266dupC() -> list[tuple[str, str]]:
    """List the codes that fire for BRCA1:c.5266dupC."""
    # TODO: think about which codes apply for a frameshift in a known
    # LoF-disease gene that is absent / extremely rare in gnomAD.
    raise NotImplementedError


def assign_codes_for_brca1_s1613g() -> list[tuple[str, str]]:
    """List the codes that fire for BRCA1:p.S1613G."""
    # TODO: this is a common-population missense — what does that mean?
    raise NotImplementedError


def acmg_combine(codes: list[tuple[str, str]]) -> str:
    """Apply the ACMG/AMP combining table to (code, strength) evidence.

    Strength tiers used here:
      pathogenic side: very_strong, strong, moderate, supporting
      benign side:     standalone, strong, supporting

    Return one of: P / LP / VUS / LB / B  (use the constants at the top).
    """
    # TODO: count by strength on each side; apply the 2015 combining rules.
    raise NotImplementedError


# Walk the two anchors:
# codes_a = assign_codes_for_brca1_5266dupC()
# codes_b = assign_codes_for_brca1_s1613g()
# print("BRCA1:c.5266dupC ->", acmg_combine(codes_a), "(expected P)")
# print("BRCA1:p.S1613G    ->", acmg_combine(codes_b), "(expected B)")
'''

STEP2_SOLUTION = '''# Reference solution — Step 2.

def assign_codes_for_brca1_5266dupC() -> list[tuple[str, str]]:
    # PVS1: null variant (frameshift) in a gene where LoF is a known disease mechanism.
    # PM2_Supporting: absent / extremely rare in gnomAD population.
    # PP5: reputable source (ClinVar) reports as pathogenic — single supporting.
    return [
        ("PVS1", "very_strong"),
        ("PM2",  "supporting"),
        ("PP5",  "supporting"),
    ]


def assign_codes_for_brca1_s1613g() -> list[tuple[str, str]]:
    # BA1: allele frequency >= 5% in gnomAD (here ~30%) — stand-alone benign.
    # BP4: computational predictors (REVEL 0.10) suggest benign.
    return [
        ("BA1", "standalone"),
        ("BP4", "supporting"),
    ]


def acmg_combine(codes: list[tuple[str, str]]) -> str:
    """Apply the 2015 ACMG/AMP combining table.

    Reference: Richards et al. 2015, Table 5. We implement the canonical
    five-class output. Strength tiers are normalised to:
      pathogenic side: very_strong, strong, moderate, supporting
      benign side:     standalone, strong, supporting
    Reported strengths use ACMG conventions (e.g. PM2_Supporting -> supporting).
    """
    # Tally by side and strength.
    p_vs = p_s = p_m = p_p = 0
    b_sa = b_s = b_p = 0
    for code, strength in codes:
        s = strength.lower()
        if code in {"BA1"} or s == "standalone":
            b_sa += 1
        elif code.startswith(("BS",)) or (code.startswith("B") and s == "strong"):
            b_s += 1
        elif code.startswith("BP") or (code.startswith("B") and s == "supporting"):
            b_p += 1
        elif code == "PVS1" or s == "very_strong":
            p_vs += 1
        elif s == "strong":
            p_s += 1
        elif s == "moderate":
            p_m += 1
        elif s == "supporting":
            p_p += 1

    # Pathogenic combinations (2015 Table 5).
    pathogenic = (
        (p_vs >= 1 and p_s >= 1) or
        (p_vs >= 1 and p_m >= 2) or
        (p_vs >= 1 and p_m == 1 and p_p == 1) or
        (p_vs >= 1 and p_p >= 2) or
        (p_s >= 2) or
        (p_s >= 1 and p_m >= 3) or
        (p_s >= 1 and p_m >= 2 and p_p >= 2) or
        (p_s >= 1 and p_m >= 1 and p_p >= 4)
    )

    # Likely Pathogenic combinations.
    likely_pathogenic = (
        (p_vs >= 1 and p_m >= 1) or
        (p_s >= 1 and p_m >= 1) or
        (p_s >= 1 and p_p >= 2) or
        (p_m >= 3) or
        (p_m >= 2 and p_p >= 2) or
        (p_m >= 1 and p_p >= 4)
    )

    # Benign combinations.
    benign = (b_sa >= 1) or (b_s >= 2)

    # Likely Benign combinations.
    likely_benign = (
        (b_s >= 1 and b_p >= 1) or
        (b_p >= 2)
    )

    # ACMG also calls out "conflicting" if both pathogenic and benign criteria
    # fire — we collapse to VUS, which is what ClinVar curators usually do.
    pside = pathogenic or likely_pathogenic
    bside = benign or likely_benign
    if pside and bside:
        return VUS
    if pathogenic:
        return P
    if likely_pathogenic:
        return LP
    if benign:
        return B
    if likely_benign:
        return LB
    return VUS


codes_a = assign_codes_for_brca1_5266dupC()
codes_b = assign_codes_for_brca1_s1613g()

print("BRCA1:c.5266dupC codes:", codes_a)
print("  -> ACMG class:", acmg_combine(codes_a), "(expected Pathogenic)")
print()
print("BRCA1:p.S1613G   codes:", codes_b)
print("  -> ACMG class:", acmg_combine(codes_b), "(expected Benign)")

assert acmg_combine(codes_a) == P, "PVS1+PM2+PP5 should give Pathogenic"
assert acmg_combine(codes_b) == B, "BA1 alone should give Benign"
'''


STEP3_MD = """## Step 3 (15 min) — Auto-assign codes + run the classifier on all 25 variants

Now lift the hand-walked logic into a tiny rule engine that takes a
variant dict (Step 1 schema) and emits its evidence-code list. Then run
the combining table on all 25 anchors and tally classifications vs the
prior ClinVar labels.

We implement the **subset** of codes that the Step 1 schema can support
on its own — without family pedigree, functional assays, or co-segregation
data those codes simply cannot be assigned from the data we have:

- **PVS1** — null variant (frameshift / nonsense / canonical splice) in a
  gene where LoF is a known disease mechanism.
- **PM1** — variant in a known mutational hotspot / functional domain.
- **PM2_Supporting** — absent or extremely rare in gnomAD
  (popmax / global AF < 1e-4).
- **PM5** — different missense change at the same residue already P/LP.
- **PP3** — multiple computational predictors agree on deleterious
  (REVEL ≥ 0.7 **or** SpliceAI ≥ 0.5).
- **BA1** — allele frequency ≥ 5 % in gnomAD.
- **BS1** — allele frequency higher than expected for the disease
  (here: AF ≥ 1 %, below the BA1 cutoff).
- **BP4** — multiple computational predictors agree on benign
  (REVEL ≤ 0.15 **and** SpliceAI ≤ 0.1).

The rule engine drops back to **VUS** when no codes fire — that is the
correct ACMG/AMP default for "no actionable evidence".
"""

STEP3_TODO = '''# ----------------------------------------------------------------------
# Step 3 — Auto rule engine + tally vs ClinVar.
# ----------------------------------------------------------------------

# ACMG/AMP thresholds we'll use throughout the exercise.
BA1_AF       = 0.05    # stand-alone benign frequency
BS1_AF       = 0.01    # strong benign frequency
PM2_AF       = 1e-4    # rare-enough for PM2_Supporting
REVEL_PATH   = 0.70    # PP3 threshold
REVEL_BENIGN = 0.15    # BP4 threshold
SPLICEAI_PATH   = 0.50
SPLICEAI_BENIGN = 0.10


def assign_codes(v: dict, af_field: str = "af_global") -> list[tuple[str, str]]:
    """Return the (code, strength) list that fires for variant `v`.

    `af_field` selects which AF column to consult — global, eur, afr, ...
    Step 4 will exercise the ancestry switch by re-running with a different
    af_field.
    """
    # TODO: implement the 8 rules described in the markdown above.
    raise NotImplementedError


def classify(v: dict, af_field: str = "af_global") -> tuple[str, list[tuple[str, str]]]:
    """Return (class_label, codes) for variant v."""
    # TODO: call assign_codes() and then acmg_combine() from Step 2.
    raise NotImplementedError


# results = [classify(v) for v in variants]
# Add the result columns to the DataFrame; cross-tabulate against ClinVar.
'''

STEP3_SOLUTION = '''# Reference solution — Step 3.

BA1_AF       = 0.05
BS1_AF       = 0.01
PM2_AF       = 1e-4
REVEL_PATH   = 0.70
REVEL_BENIGN = 0.15
SPLICEAI_PATH   = 0.50
SPLICEAI_BENIGN = 0.10


def assign_codes(v: dict, af_field: str = "af_global") -> list[tuple[str, str]]:
    codes: list[tuple[str, str]] = []
    af = v.get(af_field) or 0.0
    consq = v.get("consequence", "")
    lof_disease = bool(v.get("gene_is_lof_disease"))
    revel = v.get("revel")
    spliceai = v.get("spliceai")

    # ---------- Pathogenic side ----------
    # PVS1: null variant in LoF-disease gene.
    if lof_disease and consq in {"nonsense", "frameshift", "splice"}:
        codes.append(("PVS1", "very_strong"))

    # PM1: hotspot / functional domain.
    if v.get("in_hotspot"):
        codes.append(("PM1", "moderate"))

    # PM2_Supporting: rare / absent in gnomAD.
    if af < PM2_AF:
        codes.append(("PM2", "supporting"))

    # PM5: different missense change at same residue is P/LP.
    if consq == "missense" and v.get("same_aa_as_path"):
        codes.append(("PM5", "moderate"))

    # PP3: predictors agree on deleterious (REVEL or SpliceAI).
    pp3_fired = (
        (revel is not None and revel >= REVEL_PATH)
        or (spliceai is not None and spliceai >= SPLICEAI_PATH)
    )
    if pp3_fired:
        codes.append(("PP3", "supporting"))

    # ---------- Benign side ----------
    # BA1: stand-alone benign frequency.
    if af >= BA1_AF:
        codes.append(("BA1", "standalone"))
    elif af >= BS1_AF:
        codes.append(("BS1", "strong"))

    # BP4: predictors agree on benign.
    bp4_fired = (
        (revel is not None and revel <= REVEL_BENIGN)
        and (spliceai is None or spliceai <= SPLICEAI_BENIGN)
    )
    if bp4_fired:
        codes.append(("BP4", "supporting"))

    return codes


def classify(v: dict, af_field: str = "af_global") -> tuple[str, list[tuple[str, str]]]:
    codes = assign_codes(v, af_field=af_field)
    return acmg_combine(codes), codes


# Run on all 25 anchors using the global allele frequency.
rows = []
for v in variants:
    cls, codes = classify(v)
    rows.append({
        "id": v["id"],
        "gene": v["gene"],
        "consequence": v["consequence"],
        "af_global": v["af_global"],
        "clinvar": v["clinvar"],
        "acmg_call": cls,
        "codes": ", ".join(c for c, _ in codes) or "(none)",
    })

calls = pd.DataFrame(rows)
print(f"Engine output on {len(calls)} variants:\\n")
print(calls.to_string(index=False))

# Cross-tabulate engine call vs ClinVar.
order = [P, LP, VUS, LB, B]
xtab = pd.crosstab(calls["clinvar"], calls["acmg_call"]).reindex(
    index=order, columns=order, fill_value=0
)
print("\\nConfusion matrix (rows = ClinVar, cols = engine):")
print(xtab.to_string())

# Exact-match accuracy and "directionally correct" accuracy (P/LP vs LB/B vs VUS).
def coarse(label: str) -> str:
    if label in {P, LP}: return "path"
    if label in {LB, B}: return "benign"
    return "vus"

calls["clinvar_coarse"]   = calls["clinvar"].map(coarse)
calls["engine_coarse"]    = calls["acmg_call"].map(coarse)
exact = (calls["clinvar"] == calls["acmg_call"]).mean()
coarse_acc = (calls["clinvar_coarse"] == calls["engine_coarse"]).mean()
print(f"\\nExact 5-class agreement:    {exact:.0%}")
print(f"Coarse path/VUS/benign acc: {coarse_acc:.0%}")
'''


STEP4_MD = """## Step 4 (15 min) — Ancestry bias and the PRS portability gap

This is the ethics aha. Re-run the rule engine on a subset of variants
using **EUR-only** allele frequencies and again with **AFR-only**
frequencies. Any variant whose classification flips is a candidate for
ancestry bias — the most common pattern is "rare in EUR (PM2 fires) but
common in AFR (BS1/BA1 fires)", which flips a Likely Pathogenic call to
Likely Benign.

Then run a small **polygenic risk score** experiment. Simulate 100 SNP
effects with a deterministic seed and compute PRS distributions for two
synthetic populations whose allele frequencies differ. The same scoring
weights produce different population mean PRS — a textbook case of
**covariate shift**.
"""

STEP4_TODO = '''# ----------------------------------------------------------------------
# Step 4 — Ancestry-bias audit + PRS portability demo.
# ----------------------------------------------------------------------

# Part A — re-classify each variant under EUR-only vs AFR-only AFs and
# count classification flips.
#
# Part B — simulate 100 SNPs with effect sizes drawn from N(0, 0.1^2)
# (deterministic seed); draw EUR and AFR allele frequencies from two
# different Beta priors; sample 1000 EUR + 1000 AFR individuals; compute
# the PRS distribution in each population.

def reclassify_under_ancestry(variants: list[dict], af_field: str) -> pd.DataFrame:
    """Run classify() on every variant using `af_field`; return a DataFrame."""
    # TODO
    raise NotImplementedError


def simulate_prs(n_snps: int = 100, n_individuals: int = 1000,
                 seed: int = 7) -> dict:
    """Return dict with keys: weights, af_eur, af_afr, prs_eur, prs_afr."""
    # TODO:
    # 1. Effects: rng.normal(0, 0.1, n_snps).
    # 2. Allele freqs: rng.beta(2,5) for EUR, rng.beta(5,2) for AFR  -> shifted means.
    # 3. Genotypes: binomial(2, af) for each individual, vectorised.
    # 4. PRS_i = sum_j (genotype_ij * weight_j).
    raise NotImplementedError
'''

STEP4_SOLUTION = '''# Reference solution — Step 4.

# ---------- Part A: ancestry-flip audit ----------

def reclassify_under_ancestry(variants: list[dict], af_field: str) -> pd.DataFrame:
    out = []
    for v in variants:
        cls, codes = classify(v, af_field=af_field)
        out.append({
            "id": v["id"],
            "af_global": v["af_global"],
            "af_used":   v.get(af_field, 0.0),
            "ancestry_field": af_field,
            "acmg_call": cls,
            "codes": ", ".join(c for c, _ in codes) or "(none)",
        })
    return pd.DataFrame(out)


eur_calls = reclassify_under_ancestry(variants, "af_eur")
afr_calls = reclassify_under_ancestry(variants, "af_afr")

audit = pd.DataFrame({
    "id":         eur_calls["id"],
    "af_eur":     eur_calls["af_used"],
    "af_afr":     afr_calls["af_used"],
    "call_eur":   eur_calls["acmg_call"],
    "call_afr":   afr_calls["acmg_call"],
    "clinvar":    [v["clinvar"] for v in variants],
})
audit["flip"] = audit["call_eur"] != audit["call_afr"]

print("Per-variant calls under EUR-only vs AFR-only allele frequencies:\\n")
print(audit.to_string(index=False))

flipped = audit[audit["flip"]]
print(f"\\n{len(flipped)} of {len(audit)} variants flip class when ancestry changes.")
print("These are the ancestry-sensitive calls — they fire PM2 in the EUR-trained")
print("frame but BS1/BA1 in the AFR frame (or vice versa).")

# Bar chart of EUR vs AFR class distribution.
fig, ax = plt.subplots(figsize=(7, 3.6))
order = [P, LP, VUS, LB, B]
eur_counts = eur_calls["acmg_call"].value_counts().reindex(order, fill_value=0)
afr_counts = afr_calls["acmg_call"].value_counts().reindex(order, fill_value=0)
x = np.arange(len(order))
ax.bar(x - 0.18, eur_counts.values, width=0.35, label="EUR-only AF", color="#3b82f6")
ax.bar(x + 0.18, afr_counts.values, width=0.35, label="AFR-only AF", color="#ef4444")
ax.set_xticks(x); ax.set_xticklabels(order, rotation=20)
ax.set_ylabel("variant count")
ax.set_title("ACMG/AMP class distribution under two ancestry frames")
ax.legend()
plt.tight_layout()
plt.show()


# ---------- Part B: PRS portability demo ----------

def simulate_prs(n_snps: int = 100, n_individuals: int = 1000,
                 seed: int = 7) -> dict:
    rng = np.random.default_rng(seed)
    # Per-SNP effect sizes — same weights apply across populations.
    weights = rng.normal(0, 0.1, n_snps)
    # Population-specific allele frequencies: EUR draws from Beta(2,5) (mean ~0.29)
    # AFR draws from Beta(5,2) (mean ~0.71). Different population structure.
    af_eur = rng.beta(2, 5, n_snps)
    af_afr = rng.beta(5, 2, n_snps)
    # Sample diploid genotypes (count of effect allele in {0,1,2}).
    geno_eur = rng.binomial(2, af_eur, size=(n_individuals, n_snps))
    geno_afr = rng.binomial(2, af_afr, size=(n_individuals, n_snps))
    prs_eur = geno_eur @ weights
    prs_afr = geno_afr @ weights
    return dict(weights=weights, af_eur=af_eur, af_afr=af_afr,
                prs_eur=prs_eur, prs_afr=prs_afr)


prs = simulate_prs()

print(f"\\nPRS demo (100 SNPs, 1000 individuals per population):")
print(f"  EUR PRS:  mean={prs['prs_eur'].mean():+.3f}  sd={prs['prs_eur'].std():.3f}")
print(f"  AFR PRS:  mean={prs['prs_afr'].mean():+.3f}  sd={prs['prs_afr'].std():.3f}")
print(f"  Mean shift = {prs['prs_afr'].mean() - prs['prs_eur'].mean():+.3f} "
      "(weights are identical; only AF distributions differ)")

fig, ax = plt.subplots(figsize=(7, 3.6))
bins = np.linspace(min(prs["prs_eur"].min(), prs["prs_afr"].min()),
                   max(prs["prs_eur"].max(), prs["prs_afr"].max()), 40)
ax.hist(prs["prs_eur"], bins=bins, alpha=0.55, label="EUR background", color="#3b82f6")
ax.hist(prs["prs_afr"], bins=bins, alpha=0.55, label="AFR background", color="#ef4444")
ax.axvline(prs["prs_eur"].mean(), color="#1e3a8a", lw=2)
ax.axvline(prs["prs_afr"].mean(), color="#7f1d1d", lw=2)
ax.set_xlabel("polygenic score (arbitrary units)")
ax.set_ylabel("individuals")
ax.set_title("Same weights, two populations — covariate shift in PRS")
ax.legend()
plt.tight_layout()
plt.show()
'''


STEP5_MD = """## Step 5 (10 min) — Which codes are ancestry-sensitive?

Take a step back from the numbers and inventory **which evidence codes**
in the ACMG/AMP framework depend on a population-frequency reference,
and what that means for clinical reporting in under-represented ancestries.

Run the small audit below — count how often each code fires under
EUR-only vs AFR-only frequencies — and write a markdown explanation
underneath connecting the result to the EE framing.
"""

STEP5_TODO = '''# ----------------------------------------------------------------------
# Step 5 — Per-code firing counts under EUR vs AFR frequencies.
# ----------------------------------------------------------------------

# TODO:
# 1. Re-run assign_codes() on every variant under af_eur and af_afr.
# 2. Count how often each code fires in each frame.
# 3. Tabulate; comment on which codes are ancestry-sensitive.

def code_firing_counts(variants: list[dict], af_field: str) -> pd.Series:
    """Return a Series: code -> number of variants where it fired under af_field."""
    raise NotImplementedError
'''

STEP5_SOLUTION = '''# Reference solution — Step 5.
from collections import Counter


def code_firing_counts(variants: list[dict], af_field: str) -> pd.Series:
    cnt: Counter[str] = Counter()
    for v in variants:
        for code, _ in assign_codes(v, af_field=af_field):
            cnt[code] += 1
    return pd.Series(cnt, dtype=int).sort_index()


fire_eur = code_firing_counts(variants, "af_eur")
fire_afr = code_firing_counts(variants, "af_afr")

# Align on the union of codes for a clean side-by-side.
all_codes = sorted(set(fire_eur.index) | set(fire_afr.index))
firing = pd.DataFrame({
    "EUR_AF": fire_eur.reindex(all_codes, fill_value=0),
    "AFR_AF": fire_afr.reindex(all_codes, fill_value=0),
})
firing["delta"] = firing["AFR_AF"] - firing["EUR_AF"]
print("Code-firing counts across the 25-variant set:\\n")
print(firing.to_string())

print()
print("The ancestry-sensitive codes are the ones whose firing count differs")
print("between the EUR and AFR columns. Concretely:")
print("  - PM2_Supporting: depends on AF cutoff (1e-4). A variant rare in EUR but")
print("    common in AFR will fire PM2 in one frame and not the other.")
print("  - BS1 / BA1: depends on AF >= 1% / 5%. Same logic in reverse.")
print()
print("PP3/BP4 (predictor agreement) are *also* ancestry-sensitive in practice")
print("because REVEL / AlphaMissense / SpliceAI were trained predominantly on")
print("European variants — but that does not show up here because our anchors")
print("share the same scores across ancestries. The training-distribution problem")
print("(L16 covariate shift) is the deeper issue.")
'''


SELFCHECK_MD = """## Self-check

Validates the load-bearing pieces of the pipeline. The synthetic-anchor
path is the canonical reference — these asserts pass on a fresh notebook
even without network access.
"""

SELFCHECK = '''# ----------------------------------------------------------------------
# Self-check — runs against whatever you defined above.
# ----------------------------------------------------------------------

# 1. Anchor set is the 25-variant canonical size.
assert len(variants) == 25, f"expected 25 anchors, got {len(variants)}"
assert df["gene"].nunique() >= 5, "expected at least 5 distinct genes"

# 2. acmg_combine: hand-picked code combinations land in the right class.
assert acmg_combine([("PVS1", "very_strong"), ("PM2", "supporting"),
                     ("PP5", "supporting")]) == P, "PVS1 + PM2 + PP5 should be Pathogenic"
assert acmg_combine([("BA1", "standalone")]) == B, "BA1 alone should be Benign"
assert acmg_combine([]) == VUS, "no evidence should be VUS"
assert acmg_combine([("PM1", "moderate"), ("PM2", "supporting")]) == VUS, \\
    "PM1 + PM2_Supporting alone should be VUS (insufficient evidence)"

# 3. Engine recovers ClinVar P / B labels (the unambiguous ends).
ids_path = {v["id"] for v in variants if v["clinvar"] == P}
calls_now = {v["id"]: classify(v)[0] for v in variants}
recovered_path = sum(1 for vid in ids_path if calls_now[vid] in {P, LP})
assert recovered_path >= len(ids_path) * 0.6, (
    f"engine recovered only {recovered_path} of {len(ids_path)} ClinVar-Pathogenic "
    "anchors as P/LP — rule engine looks broken"
)

ids_benign = {v["id"] for v in variants if v["clinvar"] == B}
recovered_b = sum(1 for vid in ids_benign if calls_now[vid] in {B, LB})
assert recovered_b >= len(ids_benign) * 0.6, (
    f"engine recovered only {recovered_b} of {len(ids_benign)} ClinVar-Benign "
    "anchors as B/LB"
)

# 4. Ancestry audit detects at least one class flip across the GENE_X/Y/Z/W anchors.
flip_audit = audit
n_flips = int(flip_audit["flip"].sum())
assert n_flips >= 1, "no ancestry flips detected — the EUR/AFR pivot is not exercising PM2/BS1"

# 5. PRS demo produces a measurable mean shift between EUR and AFR backgrounds.
shift = prs["prs_afr"].mean() - prs["prs_eur"].mean()
assert abs(shift) > 0.05, (
    f"PRS mean shift = {shift:+.3f} — expected a clear covariate-shift signal"
)

print(f"Engine recovered {recovered_path}/{len(ids_path)} pathogenic + "
      f"{recovered_b}/{len(ids_benign)} benign anchors as P/LP / LB/B.")
print(f"Ancestry audit: {n_flips} class flips across {len(flip_audit)} variants.")
print(f"PRS mean shift across populations: {shift:+.3f}")
print("✅ Self-check passed.")
'''


EE_MD = """## EE framing — Boolean logic gate network, sensor fusion, covariate shift

You implemented the canonical clinical-genomics decision pipeline:

1. **Evidence codes as sensor inputs.** Each ACMG/AMP code (PVS1, PM2,
   PP3, BA1, …) is a Boolean output from a different "sensor" — variant
   consequence, gnomAD frequency, REVEL score, SpliceAI score. The codes
   are deliberately heterogeneous but the rule engine treats them
   uniformly as binary inputs with a known strength tier.
2. **Combining table as Boolean logic.** The 2015 Richards et al. table
   is a finite Boolean function of the per-code firing pattern — no
   learnt parameters, no hidden state. Under the Tavtigian 2018
   reformulation, each strength tier is an independent likelihood ratio
   (PVS1 ≈ 350, PS ≈ 19, PM ≈ 4.3, PP ≈ 2.1) and the combining table
   collapses to multiplying likelihood ratios and thresholding the
   posterior. The Boolean form exists because **clinicians need to
   audit it**; the Bayesian form exists because **that's what's
   mathematically true**. Both produce the same five-class output.
3. **Ancestry bias as covariate shift.** The Step 4 audit changes one
   thing — which sub-population's allele frequency we feed the engine —
   and the class output flips. In EE terms this is a textbook
   covariate-shift failure: a classifier trained on $p(x \\mid \\text{EUR})$
   silently mis-classifies samples drawn from $p(x \\mid \\text{AFR})$,
   even though the *decision boundary* is the same. The PRS portability
   gap is the same phenomenon with a regression model — the weights are
   correctly estimated under one population's allele-frequency spectrum
   and become biased predictors under another.
4. **Auditability vs accuracy.** The whole point of using an explicit
   rule table instead of an end-to-end neural classifier is that every
   call can be traced to the codes that fired. Regulators (FDA, CE-IVDR)
   require this. ML-based predictors (REVEL, AlphaMissense) feed *one*
   supporting-level code — they inform but never replace the auditable
   rule engine.
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
