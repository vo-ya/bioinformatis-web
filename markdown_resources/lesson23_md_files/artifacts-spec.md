# Lecture 23 (proposed L18) — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 23 (Metagenomics and the Microbiome).
> **Companion files**: `lecture-23.md`, `figures-spec.md`.

## Conventions (lecture-wide)

- Each artifact is a single self-contained HTML file in `artifacts/lecture-23/NN-name.html`.
- Vanilla HTML / CSS / JavaScript; no build step.
- Tokens via `../_shared/artifact-theme.css`.
- **`<script src="../_shared/resize.js" defer></script>` exactly once near `</body>`.**
- Phylum colours: Firmicutes amber, Bacteroidetes cobalt, Proteobacteria red, Actinobacteria green, Verrucomicrobia teal.
- Diversity-metric colours: alpha violet, beta cobalt.
- Typography: Inter for chrome; JetBrains Mono for taxonomy, sequences, abundances.
- Default state instructive; outcome banner; "Educational tool — not a clinical microbiome diagnostic" disclaimer.

## Artifact budget — 7 interactive tools

| # | Title | Anchor |
|---|---|---|
| 1 | 16S vs shotgun design chooser | §1 |
| 2 | DADA2 ASV denoising stepper | §2.2 |
| 3 | Kraken2 k-mer classifier demo | §3.2 |
| 4 | Alpha + beta diversity calculator | §4 |
| 5 | Microbiome differential abundance explorer | §5 |
| 6 | UniFrac vs Bray-Curtis comparator | §4.2 |
| 7 | MAG quality assessment | §3.5 |

---

## Artifact #1 — 16S vs Shotgun Design Chooser

**File**: `artifacts/lecture-23/01-design-chooser.html`
**Anchor**: §1

### Teaching purpose

Pick study constraints (cost budget, biomass level, taxonomic resolution, functional content); artifact recommends 16S or shotgun.

### UI layout

- 4 toggle/slider inputs: budget, biomass, resolution, functional-needs.
- Output: recommendation card + reasoning.
- Cost/sample readout, expected resolution, expected analysis complexity.
- Outcome banner: "for low-biomass skin samples on a tight budget, 16S V4 is the right choice — strain-level resolution sacrificed for affordability."

### Target aha

Design choice is driven by question + constraints; no universally-best assay.

### Acceptance criteria

- 5 distinct presets demonstrate correct routing.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #2 — DADA2 ASV Denoising Stepper

**File**: `artifacts/lecture-23/02-dada2-stepper.html`
**Anchor**: §2.2

### Teaching purpose

Walk through DADA2's denoising on a small toy read set with intentional sequencing errors. Show how single-base differences are resolved as ASVs vs noise.

### UI layout

- Input: 100 simulated reads from 3 ground-truth sequences with 1% per-base error.
- Step 1: error-rate learning visualisation.
- Step 2: ASV inference with pairwise distance threshold.
- Output: 3 recovered ASVs vs ground truth.
- Toggle: 97% OTU clustering — shows fewer, less precise units.
- Outcome banner: "DADA2 recovers all 3 ASVs at single-base resolution. 97% OTU clustering merges them into 2 OTUs."

### Target aha

ASVs are exact sequences, not 97% clusters; DADA2's error model distinguishes biology from sequencer noise.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #3 — Kraken2 k-mer Classifier Demo

**File**: `artifacts/lecture-23/03-kraken-demo.html`
**Anchor**: §3.2

### Teaching purpose

Classify simulated reads against a tiny reference database via k-mer matching + LCA aggregation.

### UI layout

- Reference DB: 5 toy genomes (different species).
- Read input: 50 simulated reads from various sources + chimeric reads.
- Slider: k-mer size (15-35).
- Per-read assignment with confidence visualisation.
- Output: read counts per species; misclassification rate.
- Outcome banner: "k=31: 92% reads correctly classified. k=15: 65% (too many random matches). k=51: 75% (too few k-mers per read)."

### Target aha

k controls precision-sensitivity trade-off; LCA gracefully handles ambiguous reads.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #4 — Alpha + Beta Diversity Calculator

**File**: `artifacts/lecture-23/04-diversity-calc.html`
**Anchor**: §4

### Teaching purpose

Compute diversity metrics on user-defined community compositions; compare across samples.

### UI layout

- Input: 5 sample-composition sliders (5 species × 5 samples).
- Output: per-sample alpha (Shannon, Simpson, observed richness, Faith's PD).
- Beta-diversity heat map (Bray-Curtis dissimilarity matrix).
- PCoA scatter plot.
- Outcome banner: "Sample 1 is most diverse (Shannon 2.10); samples 2 and 3 are most similar (BC 0.15)."

### Target aha

Diversity metrics capture different aspects: richness, evenness, phylogenetic depth.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #5 — Microbiome Differential Abundance Explorer

**File**: `artifacts/lecture-23/05-diff-abundance.html`
**Anchor**: §5

### Teaching purpose

Run differential-abundance test on case vs control with compositional data correction.

### UI layout

- Pre-loaded 30-sample dataset (case n=15, control n=15) with 50 taxa.
- Tabs: t-test (naive) | DESeq2 (RNA-seq adapted) | ANCOM-BC (compositional).
- Volcano plot per method.
- FDR-significant taxa highlighted.
- Outcome banner: "naive t-test: 18 significant. ANCOM-BC: 5 significant — most naive hits are compositional artefacts."

### Target aha

Compositional bias inflates naive test significance; ANCOM-BC gives the right answer.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #6 — UniFrac vs Bray-Curtis Comparator

**File**: `artifacts/lecture-23/06-unifrac-vs-bc.html`
**Anchor**: §4.2

### Teaching purpose

Show how phylogenetically-aware UniFrac differs from compositional Bray-Curtis on cases where shared taxa are phylogenetically close vs far.

### UI layout

- Two example sample-pair comparisons:
  - Pair A: shared taxa cluster phylogenetically — same Bray-Curtis = 0.5, but UniFrac = 0.3 (close phylogenies).
  - Pair B: shared taxa span the tree — same Bray-Curtis = 0.5, but UniFrac = 0.7 (distant phylogenies).
- Phylogenetic-tree visualisation per pair.
- Outcome banner: "compositional dissimilarity (BC) is the same; phylogenetic dissimilarity (UniFrac) differs by 2x. Pick UniFrac when phylogeny matters."

### Target aha

Choice of metric encodes biology assumption.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #7 — MAG Quality Assessment

**File**: `artifacts/lecture-23/07-mag-quality.html`
**Anchor**: §3.5

### Teaching purpose

Score MAGs by completeness + contamination using a simulated single-copy-marker-gene assessment.

### UI layout

- Pre-loaded 5 MAGs from a simulated metagenome.
- Per-MAG: marker-gene presence table (~120 markers).
- Computed completeness + contamination.
- Quality category (high / medium / low).
- Outcome banner: "MAG 1: 95% / 1% — high quality. MAG 5: 45% / 12% — too fragmented + contaminated for downstream analysis."

### Target aha

MAG quality gates downstream usefulness; CheckM-style assessment is the standard.

### Acceptance criteria

- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Cross-artifact consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Phylum colours consistent across #4 and lecture figures.
- ASV / OTU representation consistent across #2 and figure 2.

## Testing checklist (per artifact)

Standard checklist (renders standalone; controls function; acceptance criteria pass; legible 720px → 1200px; resize.js × 1; outcome banner; disclaimer).
