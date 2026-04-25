# Lecture 20 (proposed L4) — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 20 (MSA, Phylogenetics, Comparative Genomics).
> **Companion files**: `lecture-20.md`, `figures-spec.md`.

## Conventions (lecture-wide)

- Each artifact is a single self-contained HTML file in `artifacts/lecture-20/NN-name.html`.
- Vanilla HTML / CSS / JavaScript; no build step.
- Tokens via `../_shared/artifact-theme.css`.
- **`<script src="../_shared/resize.js" defer></script>` exactly once near `</body>`.**
- Sequence colours: hydrophobic amber, polar teal, charged red.
- Selection-regime colours: positive-selection red, purifying cobalt, neutral grey.
- Synteny colours: conserved-block green, rearranged amber, inverted red.
- Typography: Inter for chrome; JetBrains Mono for sequences, codons, scores.
- Default state instructive; outcome banner on every artifact; "Educational tool" disclaimer.

## Artifact budget

Seven interactive tools.

| # | Title | Anchor |
|---|---|---|
| 1 | Progressive MSA stepper | §2.1 |
| 2 | NJ vs ML tree comparator | §3 |
| 3 | Molecular clock simulator | §4 |
| 4 | dN/dS calculator | §5 |
| 5 | MSA quality scorer | §2 |
| 6 | Synteny browser | §6.2 |
| 7 | Long-branch attraction demo | §3.5 |

---

## Artifact #1 — Progressive MSA Stepper

**File**: `artifacts/lecture-20/01-progressive-msa.html`
**Anchor**: §2.1

### Teaching purpose

Walk through progressive alignment on a small toy dataset. Each step adds one sequence to a profile; user controls speed.

### UI layout

- 5 short input sequences (preset; user-editable).
- Step-forward / step-back buttons.
- Step indicator (1 / 4).
- Current alignment display in monospace.
- Guide-tree mini-diagram on the side.
- Outcome banner: "step N: aligning profile of {A, B} with sequence C — 'once a gap, always a gap'."

### Target aha

Progressive alignment is heuristic; the early-step gaps freeze. Iterative refinement (next artifact, in spirit) gives them another chance.

### Acceptance criteria

- 4 steps to build the 5-sequence MSA.
- Visualise all-pair scores at each step.
- Highlight gap positions added at each step.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #2 — NJ vs ML Tree Comparator

**File**: `artifacts/lecture-20/02-nj-vs-ml.html`
**Anchor**: §3

### Teaching purpose

Build NJ and (simplified) ML trees on the same toy data; compare topologies; show bootstrap support.

### UI layout

- Input: a small distance matrix (preset; editable by user).
- Tabs: "NJ tree" | "ML tree" | "Bootstrap analysis".
- Each tab renders the tree as an SVG with branch lengths and bootstrap values.
- Slider: substitution-model strictness (Jukes-Cantor → K2P → GTR).
- Outcome banner: "topologies agree in X / N branches; ML and NJ disagree on the placement of {taxon Y}".

### Target aha

NJ is fast and usually correct; ML is more rigorous; both agree on well-resolved branches but can differ at hard branches.

### Acceptance criteria

- Default 6-taxon dataset.
- 100-iteration bootstrap simulator.
- Visual highlight of disagreement branches.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #3 — Molecular Clock Simulator

**File**: `artifacts/lecture-20/03-clock-simulator.html`
**Anchor**: §4

### Teaching purpose

Simulate substitution accumulation as a Poisson process. Show how branch length corresponds to time when clock holds; how relaxed clock generalises.

### UI layout

- Slider: per-site substitution rate $\mu$.
- Slider: clock relaxation level (0 = strict; 1 = strongly relaxed).
- Toggle: lineage-rate variation on/off.
- Visualisation: phylogenetic tree growing over time; substitutions appearing as ticks on each branch.
- Output: estimated divergence times with credible intervals.
- Outcome banner: "strict clock: time = branches × constant. Relaxed clock: time-credible intervals widen by ~40% under typical settings."

### Target aha

The molecular clock is a Poisson process; relaxation handles real-world rate variation; calibrated against fossils gives absolute dates.

### Acceptance criteria

- Animated tree growth.
- Toggle works; calibrated dates update.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #4 — dN/dS Calculator

**File**: `artifacts/lecture-20/04-dnds.html`
**Anchor**: §5

### Teaching purpose

Compute dN/dS for a simulated codon alignment. Detect positive selection.

### UI layout

- Pre-loaded codon alignment of 5 sequences (~100 codons).
- Sliders: rate of synonymous mutations, rate of non-synonymous mutations.
- Output: per-site $\omega$, gene-level $\omega$, likelihood-ratio test p-value.
- Visualisation: $\omega$ along the gene with positive-selection sites highlighted.
- Outcome banner: "ω = 1.4 with p < 0.001 → positive selection at sites 23, 67, 89. Real candidates for adaptive evolution."

### Target aha

dN/dS quantifies selection; LRT distinguishes neutral (ω = 1) from selected (ω ≠ 1).

### Acceptance criteria

- Slider-driven simulation; user can construct purifying / neutral / positive scenarios.
- LRT p-value computed by chi-squared on simulated likelihood differences.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #5 — MSA Quality Scorer

**File**: `artifacts/lecture-20/05-msa-quality.html`
**Anchor**: §2

### Teaching purpose

Score an MSA's quality (sum-of-pairs, column-wise conservation). Compare strict vs lenient gap costs.

### UI layout

- Pre-loaded 6-row MSA with intentional bad columns.
- Toggle: strict (BLOSUM62) vs lenient (PAM250) substitution.
- Visualisation: per-column conservation score; bad columns highlighted in red.
- Outcome banner: "MSA SP-score = X. Bad columns at positions [...] could be trimmed (or fixed by re-alignment).

### Target aha

MSA quality is column-by-column; bad columns inflate downstream tree errors. Trim before tree-building.

### Acceptance criteria

- 6-row 30-column alignment with 3 deliberately bad columns.
- Score updates with toggle.
- "Trim" button removes bad columns; new score displays.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #6 — Synteny Browser

**File**: `artifacts/lecture-20/06-synteny-browser.html`
**Anchor**: §6.2

### Teaching purpose

Visualise synteny between two genomes; identify rearrangements vs conserved blocks.

### UI layout

- Two-genome selection: human chr17 vs (mouse chr11 / rat chr10 / dog chr9).
- Dot-plot rendering with conserved-block green, rearranged amber, inverted red.
- Click on a feature: side panel showing the gene names and inferred rearrangement type.
- Outcome banner: "12 conserved synteny blocks; 3 inversions; 1 translocation. Most striking: HOXB cluster preserved in all 3."

### Target aha

Diagonal stripes = collinear blocks; off-diagonals = rearrangements; conservation reveals selection on gene order.

### Acceptance criteria

- 3 species pairs work.
- Synteny blocks identifiable.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Artifact #7 — Long-Branch Attraction Demo

**File**: `artifacts/lecture-20/07-lba-demo.html`
**Anchor**: §3.5

### Teaching purpose

Show how parsimony fails on long-branch-attraction scenarios; ML succeeds.

### UI layout

- 4-taxon configuration: A, B with short branches; C, D with long branches.
- Slider: substitution rate on long branches.
- Two trees: parsimony-inferred (often A-C close, B-D close — wrong) and ML-inferred (correct).
- Statistics: parsimony score, log-likelihood for each candidate topology.
- Outcome banner: "at high long-branch rate: parsimony places A and C as sister (incorrect); ML correctly groups A-B and C-D."

### Target aha

Parsimony has known failure modes; ML is more robust. Use ML.

### Acceptance criteria

- Default rate produces clear LBA failure.
- Slider drives rate; failure threshold crossed at ω substitution rate.
- HTML parses; JS passes `node --check`; resize.js × 1.

---

## Cross-artifact consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Selection-regime colour conventions consistent across #4 and lecture figures.
- Synteny colours consistent across #6 and lecture figures.
- All artifacts emit outcome banners and educational disclaimers.

## Testing checklist (per artifact)

Standard checklist (renders standalone; controls function; acceptance criteria pass; legible 720px → 1200px; resize.js × 1; outcome banner; disclaimer).
