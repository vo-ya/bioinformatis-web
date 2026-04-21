# Lecture 8 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 8 (Advanced Single-Cell: Trajectories, Integration, Multi-Modal).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-08.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

These conventions apply to every artifact in this lecture. Per-artifact sections below override them only when they need to.

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-08/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript. External libraries only if justified per-artifact.
- Must render standalone when opened directly in a browser.
- Embedded in the lecture page via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` near the end of `<body>`.** C6 smoke gate: exactly one occurrence.

### 1.2 Visual design

- Use tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- Cell-type colours: keep consistent with Lecture 7 palette (T cell cobalt `#1e3a8a`, B cell red `#c4342c`, monocyte amber `#b45309`, NK pink `#ec4899`, DC green `#15803d`, neuron violet `#7c3aed`).
- Typography: **Inter** for UI chrome; **JetBrains Mono** for gene symbols, rate constants, numeric values.
- Default state is instructive: the artifact opens showing a meaningful example, no user interaction required.
- Controls grouped in a panel above or to the left of the visualisation.
- No animations longer than ~400 ms.

### 1.3 Interaction model

- **Sliders / dropdowns / inputs** — editable parameters validated against sensible ranges.
- **Play / Pause / Step / Reset** — for iterative processes (integration rounds, velocity iteration).
- **Re-simulate** — where stochastic data is involved.
- Illegal input shows a quiet inline message (`--fg-muted`), not a modal.

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question at the end:

- Pseudotime → Kendall τ between inferred and true ordering.
- RNA velocity → correlation between inferred and true velocity directions.
- Batch integration → LISI (batch mixing) + cell-type separation.
- CITE-seq joint → recovered discrepant cell population in the joint embedding.
- Spot deconvolution → recovered vs true cell-type proportions per spot (mean absolute error).
- Ligand-receptor → filtered network edge count + top pair by score.

### 1.5 Feasibility gate on user input (required where input is free-form)

Artifacts accepting user input validate before running: check row formats, matching dimensions, positive-value constraints; report rejections inline with line numbers.

### 1.6 Pedagogical constraint

Every artifact produces a **specific realization** — the target aha moment named in its section. If the student plays with the artifact and doesn't land on that realization, the artifact has failed.

### 1.7 Out of scope

- No logins, accounts, or persistence between sessions.
- No telemetry.
- No external data files larger than ~100 KB.

---

## 2. Artifact #1 — Pseudotime Ordering Demo

**File**: `artifacts/lecture-08/01-pseudotime-ordering.html`
**Lecture anchor**: §1.2 Monocle, Slingshot, and PAGA
**EE framing reinforced**: principal-curve fitting as smoothing in low-dimensional space.

### Teaching purpose

Given a 2D synthetic dataset where cells have a known true developmental time (from simulation), fit a principal curve through them, compute pseudotime as arc-length along the curve, and report the Kendall τ between inferred pseudotime and the hidden ground-truth time. Student sees how close the snapshot ordering tracks the true ordering.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Dataset configuration:                                      │
│   Number of cells: [──●── 300]                              │
│   Noise level (σ):  [──●── 0.30]                            │
│   Trajectory shape: [ arc ▾ / S-curve / branch ]            │
│   [Re-simulate]                                             │
├─────────────────────────────────────────────────────────────┤
│ Pseudotime curve parameters:                                │
│   Smoothing λ:      [──●── 1.0]                             │
│   Iterations:       [──●── 20 ]                             │
│   [Fit principal curve]                                     │
├─────────────────────────────────────────────────────────────┤
│ 2D embedding (SVG):                                         │
│   Left half  — cells coloured by TRUE time                  │
│   Right half — cells coloured by INFERRED pseudotime        │
│   Principal curve drawn on both panels                      │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Kendall τ = 0.89 (perfect would be 1.00)              │   │
│ │ Curve fits cleanly along the true trajectory          │   │
│ │ Noise beyond σ > 0.6 degrades τ below 0.5             │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Number of cells: 50–1000 (default 300).
- Noise level σ: 0.0–1.0 (default 0.3).
- Trajectory shape: arc / S-curve / branch (default arc).
- Smoothing λ: 0.1–5.0 (default 1.0).
- Principal-curve iterations: 5–100 (default 20).
- Re-simulate + Fit buttons.

### What they see

- Two side-by-side scatters of the same cells. Left panel: cells coloured by true latent time. Right panel: cells coloured by inferred pseudotime. The principal curve overlaid on both in `--fg` 2 px.
- Kendall τ computed on the two orderings and displayed prominently.

### Target aha moment

Start with 300 cells, noise 0.3, arc trajectory, λ=1.0 → τ ≈ 0.9. Crank noise to 0.8 → τ drops below 0.4 because the true ordering gets washed out by within-cell-type scatter. Lower smoothing λ=0.1 with high noise → curve overfits to the noise, chasing individual points; τ stays bad. The student sees that pseudotime quality depends on (a) population density along the trajectory, and (b) smoothing hyperparameter choice.

### Technical notes

- Pure JS. Seeded PRNG (mulberry32) for reproducibility.
- Trajectory simulation: sample t uniform on [0, 1]; map (x, y) = shape_fn(t) + Gaussian noise.
- Principal-curve fit: start with PCA first component as initialisation. Iterate: (1) project each point to nearest point on current curve; (2) fit a smoothing spline or local-linear regression through the projected points ordered by arc length; (3) update curve. Use λ as smoothing kernel bandwidth.
- Kendall τ: standard rank-correlation formula over all cell pairs; report ± sign.
- Feasibility gate: warn if iterations × cells > 30,000 (slow).

### Acceptance criteria

- [ ] Default config (300 cells, σ=0.3, arc) yields τ ≥ 0.85.
- [ ] Noise slider monotonically degrades τ (on average).
- [ ] S-curve shape produces τ lower than arc due to 2D ambiguity.
- [ ] Branching shape (if implemented) shows τ degrading where the branch splits.
- [ ] Opens with default pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — RNA Velocity Vector Field

**File**: `artifacts/lecture-08/02-rna-velocity.html`
**Lecture anchor**: §2.2 The dynamical model and state-space estimator
**EE framing reinforced**: state-space estimator, linear-ODE fit on the u–s plane.

### Teaching purpose

Given a synthetic dataset of cells with simulated (spliced, unspliced) counts per gene, fit the steady-state γ per gene, compute per-cell velocity vectors, and project onto a 2D embedding. Student sees velocity arrows forming a flow field pointing from "immature" cells to "mature" cells.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Dataset:                                                    │
│   Cell count:    [──●── 300]                                │
│   Gene count:    [──●── 20 ]                                │
│   Noise level:   [──●── 0.2]                                │
│   Trajectory:    [ linear ▾ / branch ]                      │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Velocity model:                                             │
│   [ steady-state ▾ / stochastic ]                           │
│   Include cell-cycle genes: [ off ▾ / on ]                  │
├─────────────────────────────────────────────────────────────┤
│ Main panel:                                                 │
│   Left:  u-s phase diagram (one gene selected)              │
│          (cells as dots with γ steady-state line)            │
│   Right: 2D embedding with velocity arrows overlaid         │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Fitted γ (mean over genes): 0.92 (true 1.0)           │   │
│ │ Arrow-direction accuracy: 84% of cells point correctly│   │
│ │ Cell-cycle contamination increases wrong arrows to 26%│   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Cell count: 50–1000 (default 300).
- Gene count: 5–50 (default 20).
- Noise: 0.0–0.5 (default 0.2).
- Trajectory: linear / branch (default linear).
- Velocity model: steady-state / stochastic (default steady-state).
- Include cell-cycle genes: toggle to introduce oscillating genes that degrade the fit.
- Re-simulate button.

### What they see

- Left panel: u–s phase scatter for one selected gene. The steady-state line y = γ·x drawn. Cells on / above / below the line colour-coded.
- Right panel: 2D embedding (PCA-2 of spliced counts). Each cell gets a velocity arrow projected from high-dim to 2D via cosine similarity with neighbour displacements. Arrows `--accent`, 5–8 px long.
- Outcome banner: mean γ over all genes (inferred vs true = 1.0), arrow-direction accuracy (fraction of cells whose inferred velocity has positive cosine similarity with the true trajectory tangent), degradation with cell-cycle genes included.

### Target aha moment

With cell-cycle genes off, arrows form a clean flow field along the simulated trajectory; arrow-direction accuracy ≥ 80%. Turn cell-cycle genes on: arrows become noisy, some point backwards, accuracy drops to ~60%. Student sees the documented failure mode — cell-cycle genes contaminate velocity and need to be regressed out.

### Technical notes

- Pure JS. Seeded PRNG.
- Simulation: sample true time t per cell; for each gene, pick rate parameters (α_j, β_j=1, γ_j); simulate (u, s) under the ODE at time t with noise.
- Steady-state γ fit: linear regression of u on s through origin using the top quantile (say, top 95% of s values) — standard velocyto heuristic.
- Velocity: v_s = β·u − γ·s per gene per cell.
- Projection: for each cell c, find its 20 nearest neighbours in gene space; project velocity as weighted sum of neighbour-displacements cosine-aligned with v(c).
- Feasibility gate: warn if gene count < 5 (too few to fit a meaningful trajectory).

### Acceptance criteria

- [ ] Default linear trajectory with cell-cycle off produces ≥80% arrow-direction accuracy.
- [ ] Turning cell-cycle on visibly reduces accuracy (>10% drop).
- [ ] u-s phase diagram shows the steady-state line and off-diagonal cells.
- [ ] Gene selector switches the u-s panel to a different gene.
- [ ] Opens with default pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — Harmony vs scVI Batch Integration

**File**: `artifacts/lecture-08/03-batch-integration.html`
**Lecture anchor**: §3.2–§3.3 Harmony and scVI
**EE framing reinforced**: source separation — same biological source, two batch conditions.

### Teaching purpose

Given two simulated batches of the same 3 cell types with batch-specific mean shifts, compare three embeddings: (a) naive concatenation, (b) Harmony-style linear correction in PCA space, (c) a simplified scVI-style encoder that conditions on batch. Student sees why naive fails and how each method recovers cell-type structure.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Dataset:                                                    │
│   Cells per batch:    [──●── 200 ]                          │
│   Batch effect size:  [──●── 2.0 ]                          │
│   Cell-type separation: [──●── 1.5 ]                        │
│   Number of batches:   [──●── 2 ]                           │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Three parallel UMAPs:                                       │
│   (a) naive concat — batches separate                       │
│   (b) Harmony      — iterative cluster-centroid correction  │
│   (c) scVI-style   — VAE with batch conditioning            │
├─────────────────────────────────────────────────────────────┤
│ Metrics:                                                    │
│   LISI (batch mixing, higher is better)                     │
│   ASW (cell-type separation, higher is better)              │
│   Table of all three across methods                         │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ naive:   LISI 1.08  |  ASW 0.34 (batch dominates)     │   │
│ │ Harmony: LISI 1.83  |  ASW 0.56 (well mixed)          │   │
│ │ scVI:    LISI 1.91  |  ASW 0.62 (best here)           │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Cells per batch: 50–500 (default 200).
- Batch effect size: 0.0–5.0 in σ units (default 2.0).
- Cell-type separation: 0.5–5.0 (default 1.5).
- Number of batches: 2–4 (default 2).
- Re-simulate button.

### What they see

- Three UMAPs side-by-side. Cells coloured by cell type in all three (batches distinguishable by marker shape — circles vs squares vs triangles).
- A metrics table below: LISI(batch), LISI(cell-type), ASW, per method.
- Outcome banner contrasting the three.

### Target aha moment

At default batch effect = 2.0, naive concatenation shows 2 clear batch lobes, LISI ~1.0 (one batch per neighbourhood, the failure mode). Harmony brings LISI to ~1.8, batches mixed; but if you crank batch effect to 5.0, Harmony's linear correction starts over-smoothing cell-type boundaries (ASW drops). scVI-style continues to handle it because the decoder can fit nonlinear batch effects through the conditional. Student sees the limits of linear correction.

### Technical notes

- Pure JS. Seeded PRNG.
- Data simulation: 3 cell types, each with a mean gene-expression profile (Gaussian in 50D); per-batch add a shared additive shift in a random direction; per-cell small noise.
- PCA + UMAP: standard implementations from Lecture 7's artifact (power iteration + force-directed).
- Harmony simulation: iterative correction. Per iteration: soft-cluster cells (k-means with k=3–5), compute per-cluster per-batch centroid, shift each cell toward the cross-batch centroid weighted by cluster membership. ~20 iterations.
- scVI simulation (simplified): train a small 2-layer encoder-decoder with batch one-hot concatenated to input/output. Keep architecture tiny (latent dim 10, hidden 64) for runtime. Optimise NB-like reconstruction loss for a few hundred gradient steps; extract latent.
- LISI: for each cell, compute effective number of batches in k-NN (k=30); harmonic mean.
- ASW: standard silhouette width using cell-type labels.

### Acceptance criteria

- [ ] Default config shows naive<Harmony<scVI on LISI.
- [ ] ASW ordering reflects biological separation preservation.
- [ ] Cranking batch effect to extremes breaks Harmony before scVI.
- [ ] Adjusting cell-type separation down makes all methods harder.
- [ ] Opens with all three pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — CITE-seq Joint Embedding

**File**: `artifacts/lecture-08/04-cite-seq-joint.html`
**Lecture anchor**: §4.3 Joint embeddings
**EE framing reinforced**: sensor fusion — two modalities with different noise characteristics.

### Teaching purpose

Given simulated CITE-seq data where RNA and protein *mostly* agree but one cell population shows modality-specific structure (post-transcriptional regulation), compute RNA-only, protein-only, and WNN-style joint embeddings. Student sees which modality catches the discrepant population.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Dataset:                                                    │
│   Cells total:      [──●── 400 ]                            │
│   RNA noise:        [──●── 0.3 ]                            │
│   Protein noise:    [──●── 0.5 ]                            │
│   Discrepancy cell population size: [──●── 80 ]             │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Three UMAPs side-by-side:                                   │
│   (a) RNA-only      — 5 clusters visible                    │
│   (b) Protein-only  — 5 clusters, but discrepant population │
│                        lands differently                    │
│   (c) WNN joint     — discrepant cluster shows as distinct  │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Discrepant cells are lumped with B cells in RNA-only  │   │
│ │ WNN embedding separates them as a new cluster         │   │
│ │ ARI(true vs RNA): 0.82 | WNN: 0.94 | Protein: 0.76   │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Cells total: 100–1000 (default 400).
- RNA noise, protein noise: 0.0–1.0 (defaults 0.3, 0.5).
- Discrepancy population size: 0–150 (default 80) — cells where RNA looks like cell type A but protein looks like cell type B.
- Re-simulate.

### What they see

- Three small UMAPs coloured by *true* cell-type label.
- ARI (adjusted Rand index) computed between Leiden-style clustering of each UMAP and the true labels. Higher ARI = better recovery.
- Outcome banner states which modality revealed the discrepant population.

### Target aha moment

Raise discrepancy population size to 100. RNA-only UMAP lumps them with their "RNA-type" cluster; protein-only lumps them with their "protein-type" cluster; WNN puts them in a separate cluster that respects both modalities. ARI for WNN highest. Student sees that joint embeddings find structure neither single modality could.

### Technical notes

- Pure JS. Seeded PRNG.
- Simulation: 5 cell types in both modalities (RNA 50D, protein 20D). Discrepant cells have RNA profile = type A but protein profile = type B.
- WNN simulation: compute kNN in each modality separately, compute per-cell weight = local signal / (local signal + noise) per modality, combine into weighted graph.
- PCA + UMAP: from Lecture 7's toolkit.
- Clustering: simple k-means or Leiden on each embedding.
- ARI: standard adjusted Rand index against the true labels.

### Acceptance criteria

- [ ] Default config produces WNN ARI > RNA-only ARI.
- [ ] Raising discrepancy size widens the WNN advantage.
- [ ] Zero discrepancy makes all three UMAPs equivalent.
- [ ] Opens with default pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — Visium Spot Deconvolution

**File**: `artifacts/lecture-08/05-spot-deconvolution.html`
**Lecture anchor**: §5.2 The mixed-pixel problem and deconvolution
**EE framing reinforced**: constrained linear unmixing — same math as multispectral remote-sensing endmember unmixing.

### Teaching purpose

Given a simulated Visium dataset where each spot contains a known mixture of 4 cell types, fit per-spot proportions via constrained non-negative least squares against reference scRNA-seq profiles. Student sees inferred proportions vs ground truth and builds intuition for when deconvolution succeeds or fails.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Dataset:                                                    │
│   Number of spots:     [──●── 50 ]                          │
│   Reference signature noise: [──●── 0.2 ]                   │
│   Spot noise:           [──●── 0.3 ]                        │
│   Cells per spot:       [──●── 10 ]                         │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Left: Spatial grid of spots (coloured by dominant type)     │
│ Right: Selected spot — proportion bar chart                 │
│         (inferred + ground-truth ticks overlaid)            │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Mean abs error (proportion): 0.08                     │   │
│ │ Dominant-type accuracy: 94% (47 of 50 spots)          │   │
│ │ Two spots mis-assigned where min cell count is 2      │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Number of spots: 10–100 (default 50).
- Reference signature noise: 0.0–0.5 (default 0.2) — how noisy is the scRNA-seq reference you're using.
- Spot noise: 0.0–0.5 (default 0.3) — measurement noise in Visium counts.
- Cells per spot: 2–20 (default 10).
- Re-simulate.

### What they see

- Left panel: a grid of 50 circular spots, coloured by inferred dominant cell type. Hover/click a spot to select it.
- Right panel: bar chart for the selected spot — 4 bars (T cell, B cell, monocyte, fibroblast), inferred proportions filled, ground-truth proportions as tick marks on top of each bar.
- Outcome banner: mean absolute error of proportion estimates across all spots; dominant-type accuracy; number of mis-classified spots.

### Target aha moment

Set reference noise to 0.0 → near-perfect deconvolution, MAE < 0.05. Crank reference noise to 0.5 → deconvolution degrades fast; MAE > 0.2. Student sees that deconvolution quality is bounded by reference quality — garbage in, garbage out. Similarly, cells-per-spot = 2 produces high variance; cells-per-spot = 20 averages noise out.

### Technical notes

- Pure JS.
- Simulation: 4 cell-type reference profiles (50-gene vectors). Per-spot: randomly pick 2–20 cells by type according to a ground-truth proportion vector. Sum their profiles. Add Poisson/Gaussian noise.
- Deconvolution: constrained NNLS. For each spot x_i, minimise ||x_i − M w||² s.t. w ≥ 0, Σw = 1. Use standard Lawson-Hanson active-set NNLS (or simpler projected gradient for tutorial simplicity); reproject to simplex each step.
- Reference matrix M (genes × cell-types) from the reference profiles, optionally with added noise.
- MAE: mean of |w_inferred − w_true| across all spots and cell types.
- Feasibility gate: warn if cells per spot < 2.

### Acceptance criteria

- [ ] Default config produces MAE < 0.15.
- [ ] Reference-noise slider monotonically degrades MAE.
- [ ] Cells-per-spot slider visibly changes variance of inferred proportions.
- [ ] Clicking a spot updates the right-panel bar chart.
- [ ] Opens with default pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — Ligand-Receptor Network Browser

**File**: `artifacts/lecture-08/06-ligand-receptor.html`
**Lecture anchor**: §6.1 Ligand-receptor inference
**EE framing reinforced**: graph-based network inference from co-expression.

### Teaching purpose

Given simulated expression data for 5 cell types and a small curated ligand-receptor database, compute communication scores for every cell-type pair × L-R pair. Student explores the resulting network, filters by ligand family, thresholds edges, and sees how the network changes.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Dataset:                                                    │
│   Cells per type: [──●── 200 ]                              │
│   L-R database:  [ toy (10 pairs) ▾ / extended (30 pairs) ] │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Filters:                                                    │
│   Ligand family:       [ all ▾ | cytokine | chemokine ]     │
│   Score threshold:     [──●── 0.5 ]                         │
│   Minimum fraction expressing: [──●── 0.2 ]                 │
├─────────────────────────────────────────────────────────────┤
│ Main panel:                                                 │
│   Directed-graph SVG:                                       │
│     nodes = cell types                                      │
│     edges = ligand→receptor axes                            │
│     thickness = score; colour = family                      │
│   Side panel:                                               │
│     top 10 L-R pairs ranked by score                        │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Network shown: 7 edges after threshold 0.5            │   │
│ │ Top pair: IFNG (T cell) → IFNGR1 (macrophage)         │   │
│ │ Score 0.83; background null 0.12                      │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Cells per type: 50–500 (default 200).
- L-R database: toy (10 curated pairs) / extended (30 pairs).
- Ligand family filter dropdown.
- Score threshold slider: 0.0–1.0 (default 0.5).
- Minimum fraction expressing: 0.0–1.0 (default 0.2) — only consider pairs where ≥fraction of cells in a type express the ligand/receptor.
- Re-simulate.

### What they see

- Main SVG: 5 cell-type nodes arranged in a circle. Directed arrows between them representing active L-R axes. Arrow thickness = communication score; colour = ligand family.
- Side panel: ranked table of the top 10 L-R pairs with columns (ligand, receptor, sender cell type, receiver cell type, score, permutation p-value).

### Target aha moment

Start with score threshold 0.1 → dense network, ~20 edges, hard to read. Raise to 0.5 → 7 clear edges showing the core signalling. Apply ligand-family filter "cytokine" only → 3 edges remain. Now raise minimum-fraction-expressing to 0.5 → some edges disappear (the signalling was based on a rare subpopulation). Student sees that communication networks are thresholding-sensitive: the same data produces wildly different "mechanistic" stories at different filter settings.

### Technical notes

- Pure JS.
- Simulation: 5 cell types. Each has a characteristic ligand / receptor expression profile. Add noise.
- Toy database: 10 ligand-receptor pairs spanning 3 families (cytokine / chemokine / growth factor). Extended: 30 pairs.
- Score: geometric mean of (ligand mean in sender type) × (receptor mean in receiver type). Normalise so max = 1.
- Permutation p-value: shuffle cell-type labels; recompute; report fraction of shuffles ≥ observed. 100 permutations — precompute at dataset load.
- Network layout: fixed circular positions for cell types.
- Feasibility gate: warn if fraction-expressing threshold > all cells (empty network).

### Acceptance criteria

- [ ] Default config shows a non-empty network (≥5 edges).
- [ ] Threshold slider visibly prunes edges.
- [ ] Family filter restricts colours shown.
- [ ] Top 10 table updates when filters change.
- [ ] Opens with default pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Cross-Artifact Consistency

- All six artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Cell-type colours are consistent with Lecture 7 palette across artifacts that show clustered cells (T cell cobalt, B cell red, monocyte amber, NK pink, DC green, neuron violet).
- Any UMAP-style 2D scatter uses the same axis conventions (minimal axes, coordinates unlabelled).
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
