# Lecture 27 (proposed L11) — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 27 (Mass-Spectrometry Proteomics Primer).
> **Companion files**: `lecture-27.md`, `figures-spec.md`.

## Conventions (lecture-wide)

- Each artifact is a single self-contained HTML file in `artifacts/lecture-27/NN-name.html`.
- Vanilla HTML / CSS / JavaScript; no build step.
- Tokens via `../_shared/artifact-theme.css`.
- **`<script src="../_shared/resize.js" defer></script>` exactly once near `</body>`.**
- Spectrum colours: MS1 cobalt, MS2 amber, b-ions cobalt, y-ions red.
- Acquisition-mode colours: DDA amber, DIA cobalt, PRM green.
- Quantification colours: LFQ grey, SILAC paired (cobalt+amber), TMT multi-colour.
- Typography: Inter for chrome; JetBrains Mono for peptides, m/z, IDs.
- Default state instructive; outcome banner; "Educational tool" disclaimer.

## Artifact budget — 7 interactive tools

| # | Title | Anchor |
|---|---|---|
| 1 | Peptide fragment ion calculator | §1 |
| 2 | DDA vs DIA acquisition simulator | §2 |
| 3 | Target-decoy FDR explorer | §3.1 |
| 4 | Quantification method comparator | §4 |
| 5 | Proteomics differential abundance | §5 |
| 6 | PTM mass-shift detector | §3.3 |
| 7 | Plasma proteome biomarker explorer | §6 |

---

## Artifact #1 — Peptide Fragment Ion Calculator

**File**: `artifacts/lecture-27/01-fragment-calculator.html`
**Anchor**: §1

### Teaching purpose

Enter a peptide sequence; artifact computes b-ion and y-ion masses; shows predicted MS2 spectrum.

### UI layout

- Peptide input box (preset to LSDPYHRGSP; user can edit).
- Fragment table: b1, b2, ..., y1, y2, ... with masses.
- Predicted MS2 spectrum visualisation.
- Toggle: charge state (1+, 2+, 3+).
- Outcome banner: "LSDPYHRGSP: 9 b-ions + 9 y-ions; full sequence determinable from spectrum."

### Target aha

Peptide identification = matching observed peaks to predicted b/y masses; full series → unique sequence.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #2 — DDA vs DIA Acquisition Simulator

**File**: `artifacts/lecture-27/02-dda-vs-dia.html`
**Anchor**: §2

### Teaching purpose

Simulate a DDA and DIA acquisition on the same peptide library; compare coverage and reproducibility.

### UI layout

- Synthetic library: 1000 peptides with known intensities.
- Toggle: DDA top-N (slider 5-20) or DIA window count (slider 20-80).
- Peptide identification overlap across 3 simulated replicates.
- Output: identified peptide count + replicate concordance.
- Outcome banner: "DDA top-15 = 600 peptides identified, replicate concordance 70%. DIA 40-windows = 850 peptides, concordance 88%."

### Target aha

DIA's deterministic acquisition gives better reproducibility than DDA's stochastic top-N selection.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #3 — Target-Decoy FDR Explorer

**File**: `artifacts/lecture-27/03-target-decoy.html`
**Anchor**: §3.1

### Teaching purpose

Apply target-decoy FDR estimation to a synthetic identification set; explore score thresholds.

### UI layout

- Pre-loaded score distribution: 10,000 target hits + 10,000 decoy hits.
- Score-threshold slider.
- Live update: target hits accepted, decoy hits accepted, FDR.
- 1% / 5% FDR threshold markers.
- Outcome banner: "at score threshold 30: 7,200 target hits, 80 decoy hits → FDR 1.1%. Acceptable."

### Target aha

FDR = decoy/target at any threshold; modern proteomics tools select threshold giving exactly 1% FDR.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #4 — Quantification Method Comparator

**File**: `artifacts/lecture-27/04-quant-comparator.html`
**Anchor**: §4

### Teaching purpose

Run simulated experiment under each of LFQ, SILAC, TMT; compare CV, missing values, statistical power.

### UI layout

- Experimental design: 30 proteins, 6 samples (case n=3, control n=3), 30% true-DE.
- Tabs: LFQ | SILAC | TMT-6plex.
- Per-method: simulated abundances, missing-value pattern, computed CV, p-values for DE proteins.
- Power curves: true positive rate vs effect size for each method.
- Outcome banner: "LFQ: 60% TPR. SILAC: 85% TPR. TMT-6plex: 75% TPR (with 5% ratio compression bias)."

### Target aha

Quantification method affects power; choose by question + budget.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #5 — Proteomics Differential Abundance

**File**: `artifacts/lecture-27/05-diff-abundance.html`
**Anchor**: §5

### Teaching purpose

Run DEqMS-style differential expression analysis on a synthetic protein-abundance matrix with realistic missing values.

### UI layout

- 100 proteins × 20 samples synthetic dataset.
- Toggle: imputation strategy (none, half-min, KNN).
- t-test / limma / DEqMS comparison.
- Volcano plot per method.
- Outcome banner: "with KNN imputation: limma identifies 42 DE proteins at FDR < 0.05. With half-min: 35 (more conservative)."

### Target aha

Imputation affects DE power; missing-value handling is a real choice.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #6 — PTM Mass-Shift Detector

**File**: `artifacts/lecture-27/06-ptm-detector.html`
**Anchor**: §3.3

### Teaching purpose

Search a synthetic spectrum for PTM-modified peptides; identify the modification type from the mass shift.

### UI layout

- 5 query peptides; one is phosphorylated.
- Database with variable PTM search (phospho, acetyl, methyl).
- Identification table.
- Mass-shift histogram showing the +79.97 phospho-shift signal.
- Outcome banner: "1 of 5 peptides modified — phosphorylation at S/T/Y position. Mass shift +79.97 Da matches phospho variant."

### Target aha

PTMs are detected by characteristic mass shifts in the search expansion.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #7 — Plasma Proteome Biomarker Explorer

**File**: `artifacts/lecture-27/07-plasma-biomarkers.html`
**Anchor**: §6

### Teaching purpose

Browse the plasma proteome; pick a target abundance level; see which assays can detect it.

### UI layout

- Slider: target protein abundance (10⁻¹³ to 10⁻³ M).
- Dropdown: example proteins at that abundance.
- Output: which detection technologies (LC-MS/MS, Olink, SomaScan, Simoa) can quantify at this abundance.
- Cost / time / sensitivity comparison.
- Outcome banner: "10⁻¹⁰ M = troponin range. LC-MS/MS targeted: yes. Simoa: yes. Olink: yes. Routine LC-MS/MS untargeted: no — too low abundance."

### Target aha

Plasma proteome spans 10 orders of magnitude; assay selection depends on target abundance.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Cross-artifact consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Spectrum colour scheme consistent across #1, #2, #6 and lecture figures.
- Quantification method colours consistent across #4 and figure 4.

## Testing checklist (per artifact)

Standard checklist (renders standalone; controls function; acceptance criteria pass; legible 720px → 1200px; resize.js × 1; outcome banner; disclaimer).
