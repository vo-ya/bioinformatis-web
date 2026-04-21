# Lecture 7 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 7 (Single-Cell RNA-seq Fundamentals).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-07.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

These conventions apply to every artifact in this lecture. Per-artifact sections below override them only when they need to.

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-07/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript. External libraries only if justified per-artifact.
- Must render standalone when opened directly in a browser.
- Embedded in the lecture page via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` near the end of `<body>`.** C6 smoke gate: exactly one occurrence.

### 1.2 Visual design

- Use tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- Cell-type and cluster colours: use the same palette across artifacts for the same concepts. Distinct accent-family hues for up to ~10 clusters.
- Typography: **Inter** for UI chrome; **JetBrains Mono** for barcodes, UMIs, gene symbols, count values, numerics.
- Default state is instructive: the artifact opens showing a meaningful example, no user interaction required.
- Controls grouped in a panel above or to the left of the visualisation.
- No animations longer than ~400 ms.

### 1.3 Interaction model

- **Sliders / dropdowns / inputs** — editable parameters validated against sensible ranges.
- **Step / Play / Pause / Reset** — for iterative or sequential processes (graph construction, clustering).
- **Re-simulate** — where stochastic data is involved.
- Illegal input shows a quiet inline message (`--fg-muted`), not a modal.

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question at the end:

- Deduplication → shows **how many reads collapsed to how many molecules**.
- Knee-plot threshold → shows **how many barcodes pass / fail**.
- Normalisation → shows **how the same gene looks across three representations**.
- Clustering → shows **the cluster count and dominant markers per cluster**.
- Annotation → shows **the cell-type label and the evidence (marker genes or reference similarity)**.

### 1.5 Feasibility gate on user input (required where input is free-form)

Artifacts accepting user input validate before running: check row formats, matching dimensions, positive-value constraints; report rejections inline with line numbers.

### 1.6 Pedagogical constraint

Every artifact produces a **specific realization** — the target aha moment named in its section. If the student plays with the artifact and doesn't land on that realization, the artifact has failed.

### 1.7 Out of scope

- No logins, accounts, or persistence between sessions.
- No telemetry.
- No external data files larger than ~100 KB (some scRNA-seq examples need enough data to show 3+ cell types with visible clusters).

---

## 2. Artifact #1 — UMI Deduplication Demo

**File**: `artifacts/lecture-07/01-umi-dedup.html`
**Lecture anchor**: §1.4 Cell barcodes and UMIs
**EE framing reinforced**: UMIs as error-correcting tags surviving the amplification channel.

### Teaching purpose

Visualise how reads with the same cell barcode + UMI collapse to a single molecule count, making PCR amplification bias irrelevant to the final count.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Configuration:                                              │
│   Cells: [──●── 3]                                          │
│   mRNAs per cell: [──●── 4]                                 │
│   PCR cycles: [──●── 5 ] (2^5 = 32× amplification)          │
│   Read sampling rate: [──●── 30% ]                          │
│ [Run pipeline]                                              │
├─────────────────────────────────────────────────────────────┤
│ Step 1 — Original mRNAs (12 total):                          │
│   Cell A: mRNA₁ mRNA₂ mRNA₃ mRNA₄                           │
│   Cell B: mRNA₅ mRNA₆ mRNA₇ mRNA₈                           │
│   Cell C: mRNA₉ mRNA₁₀ mRNA₁₁ mRNA₁₂                        │
│                                                             │
│ Step 2 — Tagged with UMI at reverse transcription:          │
│   (CB_A, UMI_1) · mRNA₁                                     │
│   (CB_A, UMI_2) · mRNA₂                                     │
│   ... 12 tagged molecules total                             │
│                                                             │
│ Step 3 — PCR amplification (×32):                            │
│   12 molecules → 384 reads, each carrying (CB, UMI) tag      │
│                                                             │
│ Step 4 — Sequencing (30% sample):                            │
│   ~115 reads sequenced                                      │
│                                                             │
│ Step 5 — Deduplication by (CB, UMI):                         │
│   115 reads → 12 distinct (CB, UMI) pairs                    │
│                                                             │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ ✓ Counting recovered 12 true molecules                │   │
│ │   From 115 reads across PCR amplification (32× factor)│   │
│ │   Counting reads would report ~38 per cell (wrong);   │   │
│ │   counting distinct UMIs reports 4 per cell (correct).│   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Cells slider: 1–8 (default 3).
- mRNAs per cell slider: 2–10 (default 4).
- PCR cycles slider: 1–12 (default 5 → 32× amplification).
- Read sampling rate slider: 5–100% (default 30%).
- Run pipeline button.

### What they see

- Each step rendered as a distinct panel.
- A colour-coded visualisation of how molecules multiply and then collapse back.
- Final outcome banner contrasting the "read count" (wrong) vs "UMI count" (correct).

### Target aha moment

Crank PCR cycles to 10. 12 original molecules become 12,288 reads. Sample 30% → 3,686 reads. Dedup on (CB, UMI) → back to 12. Student sees: the PCR factor is completely absorbed by the UMI-based dedup; the final count is invariant to amplification.

### Technical notes

- Pure JS. All simulation deterministic or seeded-stochastic for reproducibility.
- Reads displayed as small coloured boxes grouped by molecule.
- No feasibility gate needed beyond slider ranges.

### Acceptance criteria

- [ ] Default config produces correct 12-molecule recovery.
- [ ] PCR cycles slider changes the read count but not the final molecule count.
- [ ] Visual grouping clearly shows which reads collapse to which molecule.
- [ ] Opens standalone with default config pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — Knee Plot Explorer

**File**: `artifacts/lecture-07/02-knee-plot.html`
**Lecture anchor**: §2.3 Empty droplets and the knee plot
**EE framing reinforced**: bimodal distribution separation as detection problem.

### Teaching purpose

Interactive knee plot with adjustable threshold; student sees how cell count and ambient fraction change the plot's shape.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Parameters:                                                 │
│   True cell count: [──●── 10000 ]                           │
│   Ambient barcode count: [──●── 500000 ]                    │
│   Cell UMI range: [──●── 1000-30000 ]                       │
│   Ambient UMI range: [──●── 10-500 ]                        │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Knee plot (SVG):                                            │
│   [ log-log scatter of barcodes sorted by UMI count ]        │
│   dashed threshold line at current slider value              │
│   real-cells colour; ambient-barcodes colour                │
├─────────────────────────────────────────────────────────────┤
│ Threshold position: [──●── auto (EmptyDrops) / manual ]     │
│ Manual threshold: [──●── 100 UMIs ]                         │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 9,872 cells kept · 500,128 discarded                  │   │
│ │ True cells recovered: 98.7%                           │   │
│ │ Ambient included as "cells": 0.03%                    │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- True cell count, ambient barcode count, cell UMI range, ambient UMI range sliders.
- Auto-threshold (EmptyDrops simulation) vs manual threshold toggle.
- Manual threshold slider.
- Re-simulate button.

### What they see

- Log-log knee plot with clear separation between "cell" region and "ambient" region.
- Live update on threshold change: counts of kept vs discarded, true-positive rate, false-positive rate.

### Target aha moment

Start with default (10k cells, 500k ambient). Note the clean knee. Now raise ambient barcode count to 5M — the tail gets much longer but the knee stays at the same UMI count. The student sees that the threshold is about UMI-count separability, not about absolute barcode numbers.

### Technical notes

- Pure JS. Simulate cell UMIs from log-normal in the cell-UMI range; ambient from power-law or log-normal in ambient range.
- Knee-point auto-detection: simple algorithm (locate the point of maximum curvature in log-log space).
- Feasibility gate: cell count ≤ ambient count (warn if violated).

### Acceptance criteria

- [ ] Default simulation produces ~10k cells clearly separated.
- [ ] Threshold slider updates both the line position and the count tallies live.
- [ ] Auto-threshold places the line at/near the knee.
- [ ] Re-simulate produces visibly different but consistent distributions.
- [ ] Opens with default pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — QC Filter Playground

**File**: `artifacts/lecture-07/03-qc-filter.html`
**Lecture anchor**: §2.5 Mitochondrial fraction and ambient RNA
**EE framing reinforced**: multi-criterion outlier filtering.

### Teaching purpose

Multiple QC sliders (total UMIs, genes-per-cell, MT%) on a simulated dataset of ~2000 cells. Student watches the population change as thresholds move.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Dataset: 2000 simulated cells, 3 cell types + debris        │
│ [Re-simulate dataset]                                       │
├─────────────────────────────────────────────────────────────┤
│ Filter thresholds (all applied jointly):                    │
│   Min total UMIs per cell:   [──●── 500 ]                    │
│   Max total UMIs per cell:   [──●── 30000 ]                  │
│   Min genes per cell:        [──●── 200 ]                    │
│   Max MT fraction (%):       [──●── 15 ]                     │
├─────────────────────────────────────────────────────────────┤
│ QC scatter / violins (SVG):                                  │
│   Left: genes vs UMIs scatter, kept cells filled + discarded│
│     ones greyed                                             │
│   Right: MT% violin per cell type, threshold line marked   │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Kept: 1,612 / 2,000 cells · Removed: 388             │   │
│ │  Removed by: Min UMI 48 · Min genes 142 · MT% 198    │   │
│ │ Cell types in kept: 3 of 3 preserved                  │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- 4 threshold sliders (min/max UMI, min genes, max MT%).
- Re-simulate dataset button.

### What they see

- Scatter of genes-per-cell vs total-UMIs, with cells colour-coded by kept/removed status.
- Violin of MT% per cell type, with threshold line.
- Live count of kept vs removed, broken down by which filter removed each.

### Target aha moment

Start with permissive thresholds — 1900 of 2000 cells pass. Tighten MT% to 5% — 40% of cells get dropped, including healthy high-mito cells like neurons (the simulated dataset includes a "neuron-like" cell type with naturally elevated MT%). Student sees that aggressive filtering removes real biology, not just noise.

### Technical notes

- Pure JS. Simulate 3 cell types + ambient debris with realistic per-cell-type UMI/gene/MT distributions.
- Feasibility gate: min < max on each slider pair.

### Acceptance criteria

- [ ] Default thresholds retain ~80% of simulated cells.
- [ ] Each slider visibly changes the kept/removed split.
- [ ] Removed-by-filter breakdown sums to total removed.
- [ ] Re-simulate produces consistent but distinct distributions.
- [ ] Opens with default thresholds pre-applied.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — Normalisation Comparison

**File**: `artifacts/lecture-07/04-normalisation.html`
**Lecture anchor**: §3.2 Log-normalisation and SCTransform
**EE framing reinforced**: variance stabilisation on overdispersed count data.

### Teaching purpose

Given a simulated 3-cell-type scRNA-seq dataset at different total UMI levels, compare raw, log-normalised, and SCTransform-style (Pearson residual) representations for the same gene.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Dataset: 3 cell types, 500 cells, varying depth             │
│ Gene to display: [ MS4A1 ▾ / CD3D / ACTB / GAPDH ]          │
│ [Re-simulate dataset]                                       │
├─────────────────────────────────────────────────────────────┤
│ Three parallel panels:                                      │
│   (a) Raw counts — bar chart per cell                        │
│   (b) Log-normalised — scaled to 10^4 + log                  │
│   (c) Pearson residuals (SCTransform-style)                  │
├─────────────────────────────────────────────────────────────┤
│ Per-cell-type violins for each representation:              │
│   Raw:     [violins showing mean-variance scaling]           │
│   Log:     [violins showing more stable variance]            │
│   Pearson: [violins showing flat variance, clear type       │
│             separation]                                     │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Gene: MS4A1 (B-cell marker)                          │   │
│ │  Raw: high in B cells, low in T cells but noisy       │   │
│ │  Log-norm: clear B vs T separation, some variance     │   │
│ │  Pearson: stable variance, high B-cell SNR            │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Gene selector: 4 presets (marker gene, housekeeping, high-variance, low-count).
- Re-simulate dataset.

### What they see

- Three side-by-side representations of the same cells.
- Violins per cell-type per representation — allowing direct visual comparison.

### Target aha moment

Select a housekeeping gene (GAPDH). Raw counts look bimodal because of depth variation. Log-norm partially stabilises it. Pearson residuals show it as nearly uniform across cell types (which is correct — housekeeping *is* uniform). Student sees that normalisation *removes* the depth-confound and reveals the biology.

### Technical notes

- Pure JS. Simulate count matrix with known per-type expression + per-cell depth noise.
- Log-normalisation: standard size-factor + log(count/sf*1e4 + 1).
- Pearson residual: use a simplified NB model fit per gene, compute residuals.
- Feasibility gate: none (preset-driven).

### Acceptance criteria

- [ ] Default shows 3 distinct cell types across three representations.
- [ ] Switching genes updates all three panels.
- [ ] Pearson residual panel shows stable variance (the teaching point).
- [ ] Opens with marker gene pre-selected.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — PCA + UMAP Embedder

**File**: `artifacts/lecture-07/05-pca-umap-embedder.html`
**Lecture anchor**: §4.3 UMAP and t-SNE for visualisation
**EE framing reinforced**: linear vs nonlinear manifold learning.

### Teaching purpose

Given a simulated 3-cell-type dataset in "gene space" (100 genes), perform PCA, show the scree, then compute a UMAP-style embedding. Compare PCA-2D to UMAP-2D.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Dataset: 3 cell types, 500 cells × 100 genes                 │
│ UMAP neighbours (k): [──●── 15 ]                            │
│ UMAP min-dist: [──●── 0.1 ]                                  │
│ [Re-simulate] [Recompute UMAP]                              │
├─────────────────────────────────────────────────────────────┤
│ Panel 1 — PCA scree (top 20 PCs):                            │
│   [ bar chart showing variance per PC ]                      │
│                                                             │
│ Panel 2 — PCA 2D (PC1 × PC2):                                │
│   [ scatter, cells coloured by true type ]                   │
│                                                             │
│ Panel 3 — UMAP 2D:                                           │
│   [ scatter, cells coloured by true type ]                   │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ PCA: 3 types visible in PC1/PC2                       │   │
│ │ UMAP: 3 tight clusters, clear boundaries              │   │
│ │ PCA captures 58% variance in first 10 PCs              │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- UMAP k-neighbours slider (5–50).
- UMAP min-dist slider (0.01–1).
- Re-simulate, Recompute UMAP buttons.

### What they see

- Scree plot showing variance per PC.
- PCA 2D projection (PC1 × PC2).
- UMAP 2D embedding computed from the top 10 PCs.
- Both scatters colour-coded by true cell type.

### Target aha moment

Compare PCA-2D and UMAP-2D on the same data. PCA shows linear separation with cells bleeding into each other. UMAP shows tight, well-separated clusters. Student sees the difference between linear variance-maximisation (PCA) and topology-preserving embedding (UMAP).

### Technical notes

- Pure JS. PCA via power iteration on the covariance matrix (or SVD via Golub-Kahan bidiagonalisation).
- UMAP via simplified algorithm (kNN graph in PCA space, force-directed layout with attraction + repulsion).
- Feasibility gate: k ≤ cell count.

### Acceptance criteria

- [ ] Default shows 3 types in both PCA and UMAP.
- [ ] Scree plot decreases monotonically.
- [ ] UMAP is cleaner than PCA (visibly tighter clusters).
- [ ] Slider changes update UMAP live.
- [ ] Opens with default pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — Leiden Clustering Visualizer

**File**: `artifacts/lecture-07/06-leiden-clustering.html`
**Lecture anchor**: §4.4 Graph-based clustering
**EE framing reinforced**: graph community detection via modularity maximisation.

### Teaching purpose

Build a kNN graph on a simulated dataset, visualise it on UMAP coordinates, run Leiden at adjustable resolution. Watch clusters form, split, and merge as resolution changes.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Dataset: Simulated 500 cells in 2D (UMAP-like layout)        │
│ kNN neighbours (k): [──●── 15 ]                             │
│ Leiden resolution: [──●── 1.0 ]                             │
│ [Rebuild graph] [Recluster] [Reset]                         │
├─────────────────────────────────────────────────────────────┤
│ Main view (SVG):                                            │
│   Cells as coloured dots (colour = cluster assignment)       │
│   kNN edges drawn as faint lines                            │
│   Cluster boundaries (convex hulls) in cluster colour       │
├─────────────────────────────────────────────────────────────┤
│ Resolution slider effects:                                  │
│   Low (0.1): 2 clusters                                     │
│   Medium (1.0): 5 clusters                                  │
│   High (2.0): 10 clusters                                   │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Resolution: 1.0 · Clusters: 5                         │   │
│ │ Modularity score: 0.62                                │   │
│ │ Median cluster size: 98 cells                         │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- kNN k slider (5–30).
- Leiden resolution slider (0.1–3.0).
- Rebuild graph button (regenerates kNN from scratch).
- Recluster button (rerun Leiden with new resolution).

### What they see

- 500 cells as a 2D scatter.
- kNN edges drawn as thin grey lines (toggleable).
- Cluster assignments as colour-coded cell fills.
- Optional convex-hull outlines around clusters.

### Target aha moment

Resolution 0.5 → 3 clusters. Resolution 1.5 → 7 clusters (same dataset). Raise to 2.5 → fragmentation; 12 small clusters, several of size < 10. Student sees: resolution is a free parameter; more clusters ≠ better clusters; over-splitting produces clusters with no biological interpretation.

### Technical notes

- Pure JS. kNN computed via brute-force distance.
- Leiden simplified: Louvain-style modularity optimisation with iterative refinement until stable.
- Feasibility gate: k ≤ cell count; resolution > 0.

### Acceptance criteria

- [ ] Resolution 1.0 on default data yields ~5 clusters.
- [ ] Raising resolution increases cluster count.
- [ ] Lowering resolution decreases cluster count.
- [ ] Rebuild graph produces different (but consistent) result.
- [ ] Opens with default resolution 1.0 pre-clustered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Cross-Artifact Consistency

- All six artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Cell-type colours are consistent across artifacts 3, 4, 5, 6 (e.g., "B cell" is always the same hue).
- UMAP-style 2D scatters in artifacts 5 and 6 use the same axis conventions.
- Every artifact emits an **outcome banner** per convention §1.4.

## 9. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls.
- [ ] `<script src="../_shared/resize.js" defer></script>` embedded near `</body>`.
- [ ] Outcome banner or equivalent verdict line visible at the end of any user interaction.
- [ ] User-input artifacts pre-flight inputs with explicit pass/fail messaging.
