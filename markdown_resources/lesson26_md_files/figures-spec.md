# Lecture 26 (proposed L24) — Figures Specification

> **Scope**: Static diagrams for Lecture 26 (Drug Discovery and Chemoinformatics).
> **Companion files**: `lecture-26.md`, `artifacts-spec.md`.

## Conventions

- Filenames `NN-name-kebab.svg` zero-padded.
- Each figure legible at 720 px; scales to 1200 px.
- Atom colours (CPK convention): C grey, N cobalt, O red, S amber, H white, halogens green.
- Drug-likeness colours: drug-like cobalt, lead amber, hit grey.
- Pipeline-stage colours: discovery cobalt, lead-opt amber, pre-clinical green, clinical violet, approved red.
- Typography: Inter for chrome; JetBrains Mono for SMILES, IDs, IC50 values.

## Figure budget — 12 figures

| # | Title | Part |
|---|---|---|
| 1 | Drug discovery pipeline | Part 1 |
| 2 | SMILES → graph → fingerprint | Part 2 |
| 3 | Morgan fingerprint algorithm | Part 3 |
| 4 | Vina docking pose generation | Part 4 |
| 5 | ADMET radar chart | Part 5 |
| 6 | Drug-discovery DL workflow | Part 6 |
| 7 | Chemical space exploration | Part 1 |
| 8 | Tanimoto similarity matrix | Part 3 |
| 9 | Lipinski Rule of Five compliance | Part 5 |
| 10 | GNN message passing on a molecule | Part 6 |
| 11 | Hit-to-lead optimisation | Part 7 |
| 12 | AlphaFold-driven design pipeline | Part 6 |

---

## Figure 1 — Drug discovery pipeline

**File**: `diagrams/lecture-26/01-pipeline.svg`
**ViewBox**: `0 0 1200 600`

Horizontal flow chart with stages and approximate timelines:

- Target ID (1y): Genetics, GWAS, CRISPR.
- Hit discovery (2y): HTS / virtual screening.
- Hit-to-lead (1-2y): SAR + ADMET.
- Lead optimisation (2-3y): med chem.
- Pre-clinical (1-2y): PK / tox / efficacy.
- Phase 1-3 (5-7y).
- FDA review (1y).

Each stage with attrition rate (% advancing). Total timeline: 10-15 years; cost: $1-3B per approved drug.

Annotation: "computational methods accelerate every stage; pre-clinical and clinical remain the rate-limiting bottlenecks."

---

## Figure 2 — SMILES → graph → fingerprint

**File**: `diagrams/lecture-26/02-representations.svg`
**ViewBox**: `0 0 1200 540`

Three side-by-side panels showing the same molecule (aspirin) in three forms:

- Panel 1: SMILES string `CC(=O)Oc1ccccc1C(=O)O` with annotation of atom/bond grammar.
- Panel 2: 2D molecular graph drawn in CPK colours; atoms labelled.
- Panel 3: 1024-bit Morgan fingerprint visualised as a binary heat map.

Side annotation: "SMILES = serialisation; graph = working representation; fingerprint = ML feature vector."

---

## Figure 3 — Morgan fingerprint algorithm

**File**: `diagrams/lecture-26/03-morgan.svg`
**ViewBox**: `0 0 1200 600`

A small molecule (~6 atoms). Iteration steps:

- Iteration 0: each atom has its initial identifier (element + valence).
- Iteration 1: each atom incorporates immediate neighbours.
- Iteration 2: each atom incorporates 2-bond-distant neighbours.
- All unique atom-environment IDs hashed into a 1024-bit Morgan fingerprint.

Annotation: "ECFP4 = Morgan radius 2; ECFP6 = radius 3."

---

## Figure 4 — Vina docking pose generation

**File**: `diagrams/lecture-26/04-vina-docking.svg`
**ViewBox**: `0 0 1200 720`

Top: receptor pocket with surface coloured by hydrophobicity / charge.

Below: 5 candidate ligand poses drawn in different orientations within the pocket, each with its computed Vina score (kcal/mol).

Best pose highlighted; H-bond contacts annotated with dashed lines.

Side panel: schematic of Vina scoring function: VdW + electrostatic + H-bond + hydrophobic.

---

## Figure 5 — ADMET radar chart

**File**: `diagrams/lecture-26/05-admet-radar.svg`
**ViewBox**: `0 0 1200 720`

Spider/radar chart with 8 ADMET axes:

- Oral bioavailability.
- Plasma protein binding.
- BBB permeability.
- CYP450 inhibition.
- Hepatotoxicity.
- hERG blockade.
- AMES mutagenicity.
- Half-life.

Two example molecules overlaid:

- Aspirin (mostly green / acceptable across all axes).
- A toxicology-flagged compound (red on hERG and CYP3A4).

Annotation: "drug optimisation = pushing the radar outward across all axes simultaneously."

---

## Figure 6 — Drug-discovery DL workflow

**File**: `diagrams/lecture-26/06-dl-workflow.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom flowchart:

1. Target structure (AlphaFold or experimental).
2. Pocket detection.
3. Generative model (REINVENT / MolDiff / Pocket2Mol) → ~10⁵ candidate molecules.
4. ADMET filter (Chemprop ML predictor) → ~10⁴.
5. Docking + rescoring → ~10³.
6. Manual triage / synthetic accessibility → ~10² to synthesise.
7. Experimental testing → ~5-15% hit rate.

Annotation: "DL accelerates the upstream funnel; experimental validation remains the bottleneck."

---

## Figure 7 — Chemical space exploration

**File**: `diagrams/lecture-26/07-chemical-space.svg`
**ViewBox**: `0 0 1200 720`

A 2D t-SNE / UMAP embedding of ~10,000 ZINC molecules:

- Coloured by drug-likeness (QED).
- Annotated regions: drug-like central cluster, oral bioavailability tail, NCE cluster, fragments.

Several known drugs marked as labelled points (aspirin, imatinib, sotorasib, etc.).

Annotation: "drug-like chemical space is bounded; most synthesisable molecules don't satisfy ADMET; design = guided exploration."

---

## Figure 8 — Tanimoto similarity matrix

**File**: `diagrams/lecture-26/08-tanimoto.svg`
**ViewBox**: `0 0 1200 720`

Heat map of pairwise Tanimoto similarities for 30 molecules from a chemical series:

- Diagonal: 1.0 (self).
- Most off-diagonal: 0.3-0.6 (random pairs).
- Visible block of high-similarity (~0.85+) corresponding to a single SAR series.

Side panel: dendrogram from hierarchical clustering by Tanimoto distance.

---

## Figure 9 — Lipinski Rule of Five compliance

**File**: `diagrams/lecture-26/09-lipinski.svg`
**ViewBox**: `0 0 1200 540`

Bar chart of % of FDA-approved drugs vs Lipinski criteria:

- MW ≤ 500: ~85%.
- logP ≤ 5: ~88%.
- HBD ≤ 5: ~95%.
- HBA ≤ 10: ~93%.

Annotation: "Lipinski is empirical; ~10% of approved drugs violate ≥ 1 rule. Useful for triage, not strict filter."

---

## Figure 10 — GNN message passing on a molecule

**File**: `diagrams/lecture-26/10-gnn-mp.svg`
**ViewBox**: `0 0 1200 600`

A small molecule's graph + 3 iterations of message passing:

- Iteration 0: per-atom feature vectors (element, charge, hybridisation).
- Iteration 1: each atom updates from immediate neighbours.
- Iteration 2: each atom updates from 2-bond neighbourhood.
- Iteration 3: graph-level pooling produces a single molecule-level embedding.

Side panel: how the embedding is fed to a property-prediction head (regression for IC50 or binary for ADMET).

---

## Figure 11 — Hit-to-lead optimisation

**File**: `diagrams/lecture-26/11-hit-to-lead.svg`
**ViewBox**: `0 0 1200 720`

Iterative SAR optimisation:

- Hit (μM IC50, 1 modification away from a known active).
- Round 1 substituent variations: 6 analogues, IC50s shown.
- Round 2: best from round 1, further variations.
- Round 3: lead candidate (sub-μM IC50, improved ADMET).

Visual: tree of related molecules with potency annotations.

Annotation: "medicinal chemists run ~5-10 SAR rounds per program; modern DL accelerates analogue selection."

---

## Figure 12 — AlphaFold-driven design pipeline

**File**: `diagrams/lecture-26/12-alphafold-design.svg`
**ViewBox**: `0 0 1200 720`

Top-to-bottom integration:

1. Sequence of target.
2. AlphaFold 3 structure prediction.
3. Pocket detection (cryptic + canonical).
4. Pocket-conditioned generative model → candidate molecules.
5. Vina docking → predicted poses.
6. MM-GBSA rescoring → ranked hits.
7. Synthesise + experimentally validate.

Side panel: example success — Insilico's INS018_055 (target ID via AlphaFold; molecule via REINVENT; IND-stage).
