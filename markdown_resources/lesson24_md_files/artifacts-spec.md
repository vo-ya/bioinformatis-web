# Lecture 24 (proposed L21) — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 24 (CRISPR Screens and DepMap).
> **Companion files**: `lecture-24.md`, `figures-spec.md`.

## Conventions (lecture-wide)

- Each artifact is a single self-contained HTML file in `artifacts/lecture-24/NN-name.html`.
- Vanilla HTML / CSS / JavaScript; no build step.
- Tokens via `../_shared/artifact-theme.css`.
- **`<script src="../_shared/resize.js" defer></script>` exactly once near `</body>`.**
- Screen-direction colours: positive selection red, dropout cobalt.
- Cell-line lineages: hematopoietic red, breast cobalt, lung amber, colon green.
- Score-distribution colours: essential dropout cobalt, non-essential null grey, resistance positive red.
- Typography: Inter for chrome; JetBrains Mono for sgRNA, gene, score readouts.
- Default state instructive; outcome banner; "Educational tool" disclaimer.

## Artifact budget — 7 interactive tools

| # | Title | Anchor |
|---|---|---|
| 1 | sgRNA library design + scoring | §2 |
| 2 | Pooled-screen simulator | §1 |
| 3 | MAGeCK gene-level analysis | §3 |
| 4 | DepMap dependency browser | §4 |
| 5 | MAVE → variant interpretation | §5 |
| 6 | CRISPRko vs CRISPRi comparator | §1.3 |
| 7 | From screen to drug-target picker | §7 |

---

## Artifact #1 — sgRNA Library Design + Scoring

**File**: `artifacts/lecture-24/01-sgrna-design.html`
**Anchor**: §2

### Teaching purpose

Pick a target gene; artifact returns 4 best sgRNAs by Doench Rule Set 2 + CFD off-target score.

### UI layout

- Gene dropdown (50 examples).
- Computed sgRNA candidates with on-target / off-target scores.
- Top 4 selected by combined ranking.
- Visualisation of sgRNA position along gene.
- Outcome banner: "for gene TP53: 4 sgRNAs selected, mean on-target score 0.72, max off-target predicted matches 3."

### Target aha

sgRNA quality matters for screen performance; scoring algorithms quantify trade-offs.

### Acceptance criteria

- 50 genes with ~20 candidate sgRNAs each.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #2 — Pooled-Screen Simulator

**File**: `artifacts/lecture-24/02-screen-simulator.html`
**Anchor**: §1

### Teaching purpose

Simulate a pooled screen on a small library; user controls MOI, depth, replicates. Outputs realistic count tables.

### UI layout

- Slider: MOI (0.1–1.0).
- Slider: cell coverage per sgRNA (50–1000).
- Slider: replicate count (1–4).
- Pre-set library: 200 sgRNAs targeting 50 genes (mix of essential, non-essential, resistance).
- Output: count tables for T0 and T-final.
- QC dashboard: sgRNA recovery, replicate concordance, CEG dropout.
- Outcome banner: "MOI 0.3, cov 500, n=3 replicates: 96% sgRNA recovery, 0.94 concordance, 87% CEG dropout. Pass."

### Target aha

Screen design parameters (MOI, coverage, replicates) directly determine result quality.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #3 — MAGeCK Gene-Level Analysis

**File**: `artifacts/lecture-24/03-mageck.html`
**Anchor**: §3

### Teaching purpose

Run a simplified MAGeCK on the simulated screen output. Show per-sgRNA p-values, RRA aggregation, gene-level FDR.

### UI layout

- Pre-loaded count table (or use output from artifact #2).
- Step-by-step display: per-sgRNA LFC + p-value table; RRA computation; final gene-level rankings.
- Volcano plot of all genes.
- Outcome banner: "top 5 essential genes detected at FDR < 0.05: A, B, C, D, E. Top 3 resistance genes: X, Y, Z."

### Target aha

MAGeCK produces interpretable gene rankings via principled rank aggregation.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #4 — DepMap Dependency Browser

**File**: `artifacts/lecture-24/04-depmap-browser.html`
**Anchor**: §4

### Teaching purpose

Browse DepMap dependencies for a chosen cancer type; identify lineage-selective targets.

### UI layout

- Cancer type dropdown (~10: AML, ALL, breast, lung, colon, glioma, etc.).
- Heat map of essential genes × cell lines (subset of 100 genes most relevant).
- Lineage-selective filter slider: minimum effect-size difference.
- Output: top 10 lineage-selective genes per cancer.
- Drug-target overlay: which targets have approved drugs / clinical trials.
- Outcome banner: "AML: top selective dependency = MCL1; clinical-stage drug venetoclax targets pathway."

### Target aha

DepMap reveals lineage-specific drug targets; selectivity is the key.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #5 — MAVE → Variant Interpretation

**File**: `artifacts/lecture-24/05-mave-interpretation.html`
**Anchor**: §5

### Teaching purpose

Show how saturation mutagenesis MAVE data informs variant classification.

### UI layout

- Pre-loaded BRCA1 MAVE map (Findlay 2018 simulated).
- Heat map: AA position × possible AA change.
- Click on a variant: shows fitness score, ClinVar status, ACMG/AMP code triggered.
- Side panel: how MAVE data supports PS3 (functional study damaging) for clinical interpretation.
- Outcome banner: "BRCA1 R1443*: MAVE damaging score, ClinVar pathogenic — PS3 + PVS1 + PM2 → PATHOGENIC."

### Target aha

MAVEs provide the functional evidence for clinical variant interpretation; the bridge to L17.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #6 — CRISPRko vs CRISPRi Comparator

**File**: `artifacts/lecture-24/06-modality-comparator.html`
**Anchor**: §1.3

### Teaching purpose

Run the same screen with CRISPRko vs CRISPRi simulated; compare results on essential vs non-essential genes.

### UI layout

- Pre-loaded library hits (50 genes mix).
- Toggle: knockout vs interference modality.
- Side-by-side comparison of LFC distributions.
- Output: agreement table (essential genes detected in both, in only one, etc.).
- Outcome banner: "essential genes: CRISPRko has higher LFC magnitude but is fatal — can't capture some critical genes. CRISPRi reveals all without killing the population."

### Target aha

Modality choice depends on whether you can tolerate cell death; CRISPRi for essentials.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #7 — From Screen to Drug-Target Picker

**File**: `artifacts/lecture-24/07-screen-to-target.html`
**Anchor**: §7

### Teaching purpose

Triage CRISPR screen hits for drug-target potential. Apply druggability filters.

### UI layout

- Pre-loaded list of 50 screen hits.
- Filters: druggability score (small-molecule pocket), surface-accessibility (for biologics), normal-tissue dependency (cytotoxicity risk), patent landscape, current clinical-stage drugs.
- Output: ranked list with annotations.
- Outcome banner: "from 50 screen hits: 15 druggable, 5 surface-accessible, 3 with clean cytotoxicity profile, 1 already in clinical trials. Top novel target: gene X."

### Target aha

Screen → drug requires more than essentiality; druggability + safety filter dominate.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Cross-artifact consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Score-distribution colours consistent across #2, #3, #4 and figures.
- Cancer-lineage colours consistent across #4 and figure 4.

## Testing checklist (per artifact)

Standard checklist (renders standalone; controls function; acceptance criteria pass; legible 720px → 1200px; resize.js × 1; outcome banner; disclaimer).
