# Lecture 26 (proposed L24) — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 26 (Drug Discovery and Chemoinformatics).
> **Companion files**: `lecture-26.md`, `figures-spec.md`.

## Conventions (lecture-wide)

- Each artifact is a single self-contained HTML file in `artifacts/lecture-26/NN-name.html`.
- Vanilla HTML / CSS / JavaScript; no build step.
- Tokens via `../_shared/artifact-theme.css`.
- **`<script src="../_shared/resize.js" defer></script>` exactly once near `</body>`.**
- Atom colours (CPK): C grey, N cobalt, O red, S amber, halogens green.
- Drug-likeness colours: drug-like cobalt, lead amber, hit grey.
- Pipeline-stage colours: discovery cobalt, lead amber, pre-clinical green, clinical violet, approved red.
- Typography: Inter for chrome; JetBrains Mono for SMILES, IDs, IC50 values.
- Default state instructive; outcome banner; "Educational tool" disclaimer.

## Artifact budget — 7 interactive tools

| # | Title | Anchor |
|---|---|---|
| 1 | SMILES → fingerprint converter | §2 |
| 2 | Tanimoto similarity calculator | §3.4 |
| 3 | Lipinski / QED filter | §5.2 |
| 4 | Vina docking visualiser | §4 |
| 5 | ADMET predictor | §5 |
| 6 | GNN molecular property predictor | §6.1 |
| 7 | Drug discovery pipeline simulator | §7 |

---

## Artifact #1 — SMILES → Fingerprint Converter

**File**: `artifacts/lecture-26/01-smiles-fingerprint.html`
**Anchor**: §2

### Teaching purpose

Enter a SMILES; artifact parses to graph, computes Morgan fingerprint, displays as bit pattern.

### UI layout

- SMILES text input (preset to aspirin).
- Canonical SMILES output.
- 2D molecular graph rendered.
- 1024-bit Morgan fingerprint displayed as a binary heat map.
- Property summary: MW, logP, TPSA, HBD/HBA.
- Outcome banner: "aspirin: MW = 180.16, logP = 1.19, ECFP4 has 47 set bits."

### Target aha

Different representations encode different info; SMILES → graph → fingerprint is a deterministic chain.

### Acceptance criteria

- 5 preset molecules; manual SMILES input also works.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #2 — Tanimoto Similarity Calculator

**File**: `artifacts/lecture-26/02-tanimoto.html`
**Anchor**: §3.4

### Teaching purpose

Compute pairwise Tanimoto similarities for a set of molecules; cluster by similarity.

### UI layout

- Pre-loaded 20 molecules (mix of similar series + diverse).
- Heat map of pairwise similarities.
- Hierarchical clustering dendrogram on the side.
- Slider: similarity threshold for clustering (0.5 to 1.0).
- Output: number of clusters at each threshold.
- Outcome banner: "at threshold 0.6: 4 chemical series. At 0.85: 12 series (more granular)."

### Target aha

Tanimoto thresholds define chemical series; choice of threshold reflects how much variation you'll accept.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #3 — Lipinski / QED Filter

**File**: `artifacts/lecture-26/03-lipinski-filter.html`
**Anchor**: §5.2

### Teaching purpose

Apply Lipinski's Rule of 5 + QED score to a virtual library; report pass / fail.

### UI layout

- Pre-loaded 100-molecule library.
- Per-molecule property table.
- Visual radar chart per molecule showing 5 Lipinski axes.
- Color-coded pass / fail.
- Histogram of QED scores.
- Outcome banner: "from 100 molecules: 80 pass full Lipinski; 60 with QED ≥ 0.6 (drug-like)."

### Target aha

Lipinski filters out implausible drugs early; QED gives a single drug-likeness score.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #4 — Vina Docking Visualiser

**File**: `artifacts/lecture-26/04-vina-docking.html`
**Anchor**: §4

### Teaching purpose

Pre-computed Vina docking results for a target. Browse poses; rotate the receptor; identify key interactions.

### UI layout

- Receptor: a kinase (e.g., ABL1 1IEP).
- Library: 5 ligands docked, each with 5-10 poses.
- 3D viewer (using a simple SVG-based 2D approximation, or a small WebGL).
- Per-pose Vina score displayed.
- Highlight residues within 4 Å of ligand.
- Outcome banner: "best pose for compound A: ΔG = -10.5 kcal/mol; H-bonds with Glu286 and Met318."

### Target aha

Docking generates poses with predicted binding scores; visual inspection separates plausible from artefactual.

### Acceptance criteria

- 5 ligands × 5 poses each.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #5 — ADMET Predictor

**File**: `artifacts/lecture-26/05-admet.html`
**Anchor**: §5

### Teaching purpose

For a query molecule, run a panel of ADMET ML predictions; output radar chart.

### UI layout

- SMILES input or preset.
- 8-property radar chart: oral bioavailability, plasma protein binding, BBB, CYP3A4 inhib, hepatotoxicity, hERG, AMES, half-life.
- Per-property: predicted value + confidence.
- Pass / fail per property; overall druglikeness verdict.
- Outcome banner: "this compound: PASS oral bioavailability, FAIL hERG (predicted high blockade) — flag for cardiac risk."

### Target aha

ADMET is multi-dimensional; one fail can sink an otherwise-potent compound.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #6 — GNN Molecular Property Predictor

**File**: `artifacts/lecture-26/06-gnn-predictor.html`
**Anchor**: §6.1

### Teaching purpose

Visualise GNN message passing on a molecule; show how atom-level info aggregates to molecule-level prediction.

### UI layout

- Pre-loaded molecule.
- Message-passing iteration slider (0-3).
- Per-atom embedding visualised as a small histogram.
- Final molecule embedding → property prediction (e.g., logP).
- Side panel: prediction vs true value comparison for several test molecules.
- Outcome banner: "GNN predicts logP for aspirin = 1.20 (true 1.19); for ibuprofen = 3.97 (true 3.97)."

### Target aha

GNNs learn task-specific molecular embeddings; message passing builds neighbourhood-aware features.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #7 — Drug Discovery Pipeline Simulator

**File**: `artifacts/lecture-26/07-pipeline-sim.html`
**Anchor**: §7

### Teaching purpose

Run a simulated end-to-end drug discovery pipeline. Show attrition at each stage.

### UI layout

- Stage controls: virtual screen library size, hit-to-lead success rate, ADMET pass rate, etc.
- Output: number of compounds at each stage; total cost; total time.
- Bar chart of attrition.
- "Optimise pipeline parameters" button: run 100 simulations with varied parameters; visualise the optimum.
- Outcome banner: "from 10⁶ virtual screen → ~5 IND candidates → ~1 approved. Cost: $2B; time: 11 years."

### Target aha

Drug discovery is a funnel; each stage cuts ~10x; total attrition is multiplicative.

### Acceptance criteria

- Realistic attrition rates per stage.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Cross-artifact consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- CPK atom colour scheme consistent across #1, #4 and lecture figures.
- Pipeline-stage colours consistent across #7 and figure 1.

## Testing checklist (per artifact)

Standard checklist (renders standalone; controls function; acceptance criteria pass; legible 720px → 1200px; resize.js × 1; outcome banner; disclaimer).
