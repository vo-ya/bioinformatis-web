# Lecture 10 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 10 (Methylation, Hi-C, and 3D Genome Organisation).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-10.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-10/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript.
- Must render standalone.
- Embedded in the lecture via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` exactly once near the end of `<body>`.** C6 smoke gate.

### 1.2 Visual design

- Tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- Contact matrices: heatmap using `--accent` gradient; log-scaled intensity.
- DNA letters: base-palette colours (A `#c4342c`, C `#1e3a8a`, G `#b45309`, T `#2d7a3e`).
- Methylation state markers: filled circle = methylated; open circle = unmethylated.
- Typography: Inter for UI chrome; JetBrains Mono for DNA sequences, matrix coordinates, numerical values.
- Default state is instructive: opens with a meaningful example pre-rendered.
- Controls grouped in a panel above or to the left of the visualisation.
- Animations ≤ 400 ms.

### 1.3 Interaction model

- **Sliders / toggles / dropdowns** — editable parameters with sensible ranges.
- **Step / Run / Reset** — for iterative algorithms (ICE, eigenvector power iteration, insulation-score sweep).
- **Re-simulate** — for stochastic simulations.
- Illegal input → quiet inline message (`--fg-muted`).

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question:

- Bisulfite conversion → number of CpGs correctly recoverable in the output.
- Three-letter alignment → best mapping position, methylation calls.
- DMR fit → posterior probability of group A > group B, credible intervals.
- Contact matrix viewer → identified features (TADs, compartment stripes) under current view mode.
- ICE normalisation → drop in coefficient of variation across rows between raw and normalised matrix.
- TAD insulation → precision/recall of called boundaries vs ground truth.
- A/B compartment PCA → agreement of inferred sign with ground truth (accuracy).

### 1.5 Feasibility gate on user input (required where input is free-form)

- Three-letter alignment: pasted DNA must be ACGTN only; length ≥ template length; inline error with line number if invalid.
- Other artifacts use preset + slider inputs only; feasibility gates on slider ranges where relevant.

### 1.6 Pedagogical constraint

Every artifact produces its named aha moment. If the student plays with sliders and doesn't land on it, the artifact has failed.

### 1.7 Out of scope

- No accounts, no telemetry, no network calls beyond declared CDN libraries (none here).
- No external data files > 100 KB.

---

## 2. Artifact #1 — Bisulfite Conversion Simulator

**File**: `artifacts/lecture-10/01-bisulfite-conversion.html`
**Lecture anchor**: §2.1 Bisulfite conversion chemistry
**EE framing reinforced**: known, input-dependent symbol substitution on a DNA channel.

### Teaching purpose

Let the student build intuition for the bisulfite transformation by toggling per-CpG methylation state and watching each unmethylated C flip to T, while methylated Cs remain as C. Output is the methylation call per CpG.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Reference DNA (editable or preset):                         │
│   preset: [ exon + CpG island ▾ / repeat region / random ]  │
│   [AATCGTCGAAGCGTCGAAGCG ... ]  (read-only or paste)        │
├─────────────────────────────────────────────────────────────┤
│ Methylation state per CpG (click to toggle):                │
│   [●○●●○●●●○●]  10 CpGs marked above reference              │
├─────────────────────────────────────────────────────────────┤
│ Apply chemistry:                                            │
│   [Show: before → bisulfite → PCR ]                         │
│   Stages rendered top-to-bottom:                            │
│     (1) reference + methylation state                       │
│     (2) after bisulfite (unmethylated C → U, highlighted)   │
│     (3) after PCR (U → T)                                   │
├─────────────────────────────────────────────────────────────┤
│ Conversion efficiency: [──●── 99.0% ]                       │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Recovered: 10/10 CpGs correctly called                │   │
│ │ Conversion efficiency < 98% → 0.3% false methylation  │   │
│ │ meth calls match input toggle state: ✓                │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Sequence preset dropdown or paste textarea (ACGTN only, length ≥ 40 bp).
- Per-CpG methylation-state toggles (click each circle to flip).
- Conversion efficiency slider: 90–100% (default 99%).
- Also show a "simulate 3 reads" button that draws 3 independent conversion outcomes so student can see stochastic variation from incomplete conversion.

### What they see

- Three stacked DNA strand renderings showing each stage.
- Methylated/unmethylated annotation above each CpG.
- Highlighting of base changes at each stage (warning colour for converted bases).
- Outcome banner with the recovery rate and any spurious methylation from incomplete conversion.

### Target aha

Set conversion efficiency to 95% → ~5% of unmethylated Cs retain as C spuriously, looking like methylation. Student sees how non-conversion inflates the apparent methylation rate. Also: toggling methylation state produces the expected downstream pattern in the PCR product.

### Technical notes

- Pure JS, seeded mulberry32 for the stochastic conversion.
- Preset sequences: hand-crafted 60–80 bp strings with 8–12 CpGs each.
- At each unmethylated CpG, flip C → T with probability = conversion_efficiency.
- Methylated CpGs always retain as C.
- Output methylation call per CpG = "1 if C in output, 0 if T".
- Feasibility gate on pasted sequence: ACGTN only, length ≥ 40, at least 3 CpGs.

### Acceptance criteria

- [ ] Default config produces 10/10 correct calls.
- [ ] Lowering conversion efficiency introduces false-methylation calls.
- [ ] Toggling methylation state updates the output bases.
- [ ] All three stages render cleanly.
- [ ] Pasted invalid sequence rejected with inline message.
- [ ] Opens with preset pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — Three-Letter Alignment Demo

**File**: `artifacts/lecture-10/02-three-letter-alignment.html`
**Lecture anchor**: §2.3 Three-letter alignment
**EE framing reinforced**: joint decoding under a known symbol-flip channel.

### Teaching purpose

Show why naive alignment fails and how three-letter conversion fixes it. Student pastes a bisulfite read, attempts naive alignment (high-mismatch result), toggles three-letter mode (clean alignment), then sees methylation called in original coordinates.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Reference (preset or paste): [ promoter CpG island ▾ ]       │
│ Bisulfite read: [ pre-built read from promoter ▾ / custom ] │
├─────────────────────────────────────────────────────────────┤
│ Algorithm mode:                                             │
│   [ naive 4-letter alignment ▾ / three-letter (C→T) ]       │
├─────────────────────────────────────────────────────────────┤
│ Alignment output (SVG):                                     │
│   Reference row + read row, with mismatches highlighted      │
│   Best match position, # mismatches, # aligned bases         │
├─────────────────────────────────────────────────────────────┤
│ Methylation calls (if alignment succeeded):                  │
│   CpG 1: methylated   (read C at ref C)                      │
│   CpG 2: unmethylated (read T at ref C)                      │
│   ... N CpG calls                                           │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Naive mode: best map 8 mismatches, alignment fails.   │   │
│ │ Three-letter mode: 0 mismatches, aligned successfully │   │
│ │ Decoded: 4 methylated / 3 unmethylated CpGs           │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Reference preset dropdown.
- Read preset dropdown + custom paste textarea.
- Algorithm mode toggle (naive / three-letter).
- Mismatch-tolerance slider (0–10, default 3): filters how many mismatches an aligner accepts before calling "no alignment".

### What they see

- Two-row SVG alignment visualisation (reference + read), base colours per letter, mismatches highlighted in warning.
- A table of per-CpG methylation calls below.
- Outcome banner comparing naive vs three-letter regimes.

### Target aha

Load a read that has 6 unmethylated CpGs. Naive alignment fails (too many mismatches). Switch to three-letter mode → alignment succeeds with 0 mismatches. Student sees the "decode in reduced space" pattern: hard in 4-letter DNA, trivial in 3-letter.

### Technical notes

- Pure JS.
- Naive alignment: O(n·m) Smith-Waterman-like search with scoring matrix penalising any mismatch. Find best position; report # mismatches there.
- Three-letter conversion: build a C→T version of both reference and read; re-run the alignment; same logic.
- Methylation decoding: at each original reference-C position, check the read base in original (pre-conversion) coordinates: C = methylated, T = unmethylated.
- Feasibility: pasted read ACGTN only, length ≥ 10 bp and ≤ reference length.

### Acceptance criteria

- [ ] Default config: naive mode fails, three-letter mode succeeds.
- [ ] Switching modes updates alignment visualisation and mismatch count.
- [ ] Methylation-call table updates when three-letter alignment succeeds.
- [ ] Invalid custom read rejected with inline message.
- [ ] Opens pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — DMR Beta Fit Explorer

**File**: `artifacts/lecture-10/03-dmr-beta-fit.html`
**Lecture anchor**: §3.3 methylKit, BSmooth, and DSS
**EE framing reinforced**: Beta-Binomial conjugate inference for proportions.

### Teaching purpose

Fit a Beta posterior per group at a single CpG, plot both posteriors, and compute the posterior probability that group A has higher methylation than group B. Student sees how read depth tightens confidence intervals and how to express differential significance without a classical p-value.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulated CpG counts (2 groups × 3 replicates):             │
│   Group A:  [──●── 70% ] true methylation                   │
│   Group B:  [──●── 30% ] true methylation                   │
│   Reads per sample: [──●── 30 ]                             │
│   Prior: Beta(α, β): α [──●── 1] β [──●── 1]                │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Posterior plot (SVG):                                       │
│   Beta posterior curve for group A (pooled across reps)     │
│   Beta posterior curve for group B                          │
│   Observed group means as vertical ticks                    │
│   95% credible interval shading for each                    │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ P(μ_A > μ_B | data) = 0.983                           │   │
│ │ Group A posterior mean: 0.68  (95% CI: 0.57 – 0.78)   │   │
│ │ Group B posterior mean: 0.31  (95% CI: 0.22 – 0.41)   │   │
│ │ Raising reads/sample tightens both intervals          │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- True group-A methylation: 0.0–1.0 (default 0.7).
- True group-B methylation: 0.0–1.0 (default 0.3).
- Reads per sample (all replicates): 5–200 (default 30).
- Prior α, β sliders (each 0.5–20, default 1 = uniform).
- Re-simulate button.

### What they see

- Two Beta posterior density curves overlaid on [0, 1].
- Observed sample means marked as small ticks.
- 95% credible intervals shown as shaded regions under each curve.
- Outcome banner with `P(μ_A > μ_B | data)` and both posterior means + CIs.

### Target aha

Default setup: posteriors separate cleanly, P(A>B) > 0.98. Drop reads/sample to 5 → posteriors broaden; P(A>B) drops to ~0.6 even though true means still differ. Student sees read depth matters — and posterior probabilities capture that explicitly (the two curves overlap, reflecting uncertainty).

### Technical notes

- Pure JS, seeded mulberry32.
- Per-sample per-group: sample `m_ij ~ Binomial(reads, true_μ_group)`. Sum methylated and total across 3 replicates per group.
- Posterior for each group: `Beta(α + total_m, β + total_n - total_m)`.
- `P(μ_A > μ_B)` computed via Monte-Carlo: sample 20,000 pairs from the two posteriors, count fraction where A > B.
- Credible interval: quantile-based from the Beta density (use an inverse-regularised-incomplete-beta approximation — bisection on CDF is fine here).

### Acceptance criteria

- [ ] Default config yields P(A>B) ≥ 0.95.
- [ ] Lowering reads/sample → P(A>B) drops and CIs widen.
- [ ] Setting A = B → P(A>B) ≈ 0.5.
- [ ] Curves + CIs render correctly.
- [ ] Opens pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — Hi-C Contact Matrix Viewer

**File**: `artifacts/lecture-10/04-hic-contact-matrix.html`
**Lecture anchor**: §4.4 The contact matrix as raw output
**EE framing reinforced**: reading a 2D matrix as a distance/similarity structure.

### Teaching purpose

Render a simulated contact matrix and let the student toggle between raw / log-scaled / distance-detrended views. Goal: learn what the three canonical structures look like (diagonal, distance decay, TADs, compartments) by eye.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulated genome region:                                    │
│   Bins: [──●── 80 ]                                         │
│   TADs planted: [──●── 5 ]                                  │
│   Compartment checkerboard strength: [──●── 0.6 ]           │
│   Loop (focal) count: [──●── 3 ]                            │
│   Sequencing depth (read count): [──●── 500,000 ]            │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ View mode: [ raw ▾ / log-scaled / O/E detrended ]           │
├─────────────────────────────────────────────────────────────┤
│ Heatmap (SVG, 80×80 cells):                                 │
│   colour = contact intensity                                │
│   diagonal bright                                           │
│   TADs as triangular blocks                                 │
│   compartments as checkerboard                              │
│   loops as isolated bright dots                             │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Raw mode: diagonal dominates; structure hard to see   │   │
│ │ Log mode: off-diagonal detail emerges                 │   │
│ │ O/E mode: TADs + compartments + loops all visible     │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- # bins: 40–120 (default 80).
- # TADs: 0–10 (default 5).
- Compartment strength: 0.0–1.0 (default 0.6).
- # loops: 0–10 (default 3).
- Sequencing depth: 100k – 5M (default 500k).
- View mode: raw / log / O/E.
- Re-simulate.

### What they see

- SVG heatmap (one rectangle per cell).
- Cells coloured by current view mode.
- Annotations on the right showing "diagonal = adjacent contacts", "off-diagonal checkerboard = compartments", "triangular blocks = TADs", "focal spots = loops".
- Outcome banner comparing the three views.

### Target aha

Load defaults in raw mode → mostly see diagonal + distance decay; TADs barely visible. Switch to log mode → features emerge. Switch to O/E mode → clear TAD triangles and checkerboard. Student sees that the right transformation makes the structure readable.

### Technical notes

- Pure JS, seeded mulberry32.
- Simulate bin-bin contact counts using a power-law distance decay: `C[i,j] ∝ |i-j|^(-0.9)`. Add TAD block bonuses for same-TAD pairs. Add compartment checkerboard by A/B assignment. Add loop spots at specific (i, j) pairs.
- Sample total_reads reads; distribute according to probabilities proportional to the expected C matrix; Poisson-sample each cell.
- View modes:
  - raw: `C[i,j]` linearly coloured.
  - log: `log(C[i,j] + 1)` coloured.
  - O/E: divide by expected-by-distance profile, log-scale the ratio centred at 1.

### Acceptance criteria

- [ ] Default config produces clearly visible diagonal + TAD structure in O/E mode.
- [ ] Raw mode shows diagonal dominance; log mode shows intermediate; O/E shows cleanest structure.
- [ ] Re-simulate produces visually similar but different matrices.
- [ ] Slider changes update the view live.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — ICE Normalisation Demo

**File**: `artifacts/lecture-10/05-ice-normalisation.html`
**Lecture anchor**: §5.1 Normalisation — ICE and KR
**EE framing reinforced**: iterative matrix balancing to remove multiplicative biases.

### Teaching purpose

Run iterative correction on a contact matrix with planted bin-bias artefacts. Show the banding disappear across iterations; report the row-sum uniformity metric converging.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulated matrix:                                           │
│   Bins: [──●── 60 ]                                         │
│   Bin-bias strength: [──●── 2.5× ]                          │
│   # biased bins: [──●── 8 ]                                  │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Iteration: [──●── 0 ]  (0 = raw, 20 = converged)             │
│ [Step] [Run to convergence] [Reset]                         │
├─────────────────────────────────────────────────────────────┤
│ Two heatmaps side-by-side: current matrix + target          │
│   (raw → after iteration K)                                 │
├─────────────────────────────────────────────────────────────┤
│ Convergence metric:                                         │
│   Row-sum CV across iterations (line chart)                 │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Iter 0: row-sum CV = 0.42 (strong banding)            │   │
│ │ Iter 20: row-sum CV = 0.02 (converged)                │   │
│ │ Effective removal of bin-specific visibility biases    │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- # bins: 40–80 (default 60).
- Bin-bias strength: 1.0–5.0× (default 2.5×).
- # biased bins: 0–20 (default 8).
- Iteration slider: 0–30.
- Step / Run / Reset buttons.
- Re-simulate.

### What they see

- Left heatmap: the current (iteration-K) matrix.
- Right heatmap: the converged (iteration-∞) matrix as a target reference.
- Line chart below: row-sum coefficient of variation (CV) as a function of iteration number, showing monotone decrease.
- Outcome banner summarising CV drop and visible effect.

### Target aha

Default: raw matrix shows clear banding at biased bins → after 20 iterations, banding is gone. Student sees iterative balancing converge. Raise bin-bias strength to 5× → takes more iterations to converge; cv drops more slowly. Student sees the algorithm is robust but scales with severity of bias.

### Technical notes

- Pure JS, seeded mulberry32.
- Simulate a base contact matrix with smooth distance decay and TAD structure.
- Multiply rows/columns corresponding to biased bins by the bias factor (symmetrically, so matrix stays symmetric).
- ICE iterative correction: at each iteration, compute row sums, divide each row (and column, since symmetric) by its sum / target_mean. Repeat.
- Track row-sum CV across iterations.

### Acceptance criteria

- [ ] Default: raw matrix shows banding; after 20 iterations it's gone.
- [ ] CV metric decreases monotonically across iterations.
- [ ] Step / Run / Reset buttons work correctly.
- [ ] Re-simulate produces a new biased matrix.
- [ ] Opens pre-computed with default (iter=0).
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — TAD Insulation Score Calculator

**File**: `artifacts/lecture-10/06-tad-insulation.html`
**Lecture anchor**: §5.2 TADs
**EE framing reinforced**: 2D → 1D reduction, local-minima detection as change-point analysis.

### Teaching purpose

Compute the insulation score from a simulated contact matrix with planted TAD boundaries. Detect local minima as predicted boundaries; compare to ground truth.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulated matrix:                                           │
│   Bins: [──●── 80 ]                                         │
│   # planted TADs: [──●── 6 ]                                │
│   TAD strength: [──●── 2.0× baseline ]                      │
│   Noise level: [──●── 0.2 ]                                 │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Insulation-score parameters:                                │
│   Window size (bins): [──●── 10 ]                           │
│   Minimum-call threshold: [──●── -1.0 ] (log2 ratio)         │
├─────────────────────────────────────────────────────────────┤
│ Two-panel visualisation:                                    │
│   Top: contact matrix with planted boundaries marked         │
│         (green ticks) and called boundaries (amber triangles)│
│   Bottom: insulation-score track (line plot), threshold      │
│           line, local-minima-below-threshold highlighted     │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Planted: 6 · Called: 7 · Precision 0.86 · Recall 1.00 │   │
│ │ Reducing window size → more false positives           │   │
│ │ Raising threshold → miss real boundaries              │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- # bins: 50–120 (default 80).
- # planted TADs: 3–10 (default 6).
- TAD strength: 1.2–4.0× (default 2.0).
- Noise level: 0.0–0.5 (default 0.2).
- Window size: 3–25 bins (default 10).
- Call threshold: -3.0 to 0.0 (default -1.0, log2-ratio scale).
- Re-simulate + Run buttons.

### What they see

- SVG heatmap (top) with planted-boundary green ticks below the diagonal and called-boundary amber triangles above.
- Line plot (bottom) showing insulation score per bin; threshold as a horizontal dashed line; called minima as highlighted dots.
- Outcome banner with precision / recall / sensitivity hints.

### Target aha

Default config → precision ≥ 0.8, recall ≥ 0.9. Shrink window size to 3 → many spurious minima, precision drops. Raise threshold toward 0 → miss real boundaries. Student sees the two knobs — window size and threshold — trade off noise against sensitivity.

### Technical notes

- Pure JS, seeded mulberry32.
- Simulate a contact matrix with a base distance decay and K planted TADs of given strength (elevate in-TAD contacts proportionally). Add Poisson noise scaled by noise_level.
- Insulation score at bin `i`: `log2( mean_contacts_crossing_bin_i / mean_contacts_in_full_window_at_i )`, over a window of `W` bins.
- Detect local minima where IS < threshold AND IS < 0.9 × min(neighbours).
- Compare called boundaries to planted (±1 bin tolerance).

### Acceptance criteria

- [ ] Default config → precision ≥ 0.8, recall ≥ 0.9.
- [ ] Shrinking window increases false positives.
- [ ] Raising threshold decreases recall.
- [ ] Planted and called boundaries rendered on shared x-axis.
- [ ] Opens pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Artifact #7 — A/B Compartment Eigenvector

**File**: `artifacts/lecture-10/07-ab-compartment-pca.html`
**Lecture anchor**: §6.1 Hi-C contact matrix as covariance; compartments as PCA
**EE framing reinforced**: PCA / power iteration on a covariance matrix.

### Teaching purpose

Take a simulated O/E-normalised contact matrix, compute the correlation matrix, extract the first eigenvector by power iteration, and show how the eigenvector's sign recovers A/B compartment identity. This is the cleanest "genomics-is-linear-algebra" artifact in the course.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulated chromosome:                                       │
│   Bins: [──●── 60 ]                                         │
│   # A/B switches: [──●── 8 ]                                 │
│   Compartment-contact strength: [──●── 0.8 ]                 │
│   Noise: [──●── 0.3 ]                                       │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Power iteration:                                            │
│   Iteration: [──●── 0 ]  (0 = init, 20 = converged)          │
│   [Step] [Run to convergence]                               │
├─────────────────────────────────────────────────────────────┤
│ Three panels:                                               │
│   (a) O/E contact matrix                                    │
│   (b) Correlation matrix (computed from a)                   │
│   (c) Eigenvector 1D track (current estimate)                │
│        positive / negative regions shaded                    │
│   Bottom: binarised A/B call vs ground truth                │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Eigenvalue λ₁ = 11.3                                  │   │
│ │ A/B sign agreement with truth: 58/60 bins (97%)        │   │
│ │ Raising noise → agreement drops below 80%              │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- # bins: 40–80 (default 60).
- # A/B switches along the chromosome: 4–20 (default 8).
- Compartment-contact strength: 0.0–1.5 (default 0.8).
- Noise level: 0.0–1.0 (default 0.3).
- Iteration slider: 0–30.
- Step / Run / Reset.
- Re-simulate.

### What they see

- Panel (a): contact matrix as a heatmap.
- Panel (b): correlation matrix (of columns of a) — sharper checkerboard.
- Panel (c): 1D eigenvector track, sign-shaded.
- Bottom strip: binarised A/B call vs ground-truth track (two rows, with agreement/disagreement annotation).
- Outcome banner: eigenvalue, sign-agreement percentage, noise-sensitivity hint.

### Target aha

Default noise: eigenvector sign perfectly recovers A/B assignment. Raise noise to 0.8 → agreement drops; power iteration still converges but to a less informative eigenvector. Student sees PCA on a covariance matrix is the right framing — when the signal is strong, the first PC is the compartmentalisation direction; when noise dominates, the first PC is contaminated.

### Technical notes

- Pure JS, seeded mulberry32.
- Simulate A/B labels along `N` bins with `k` switches. Build O/E-like matrix where same-label bins have elevated expected contacts by `compartment_strength`, different-label bins have depressed contacts.
- Add Gaussian noise scaled by noise_level.
- Compute correlation matrix: for each pair (i,j), correlate row i with row j across all other bins.
- Power iteration: start with a random vector; iteratively `v ← C v`, normalise; 20 iterations suffices.
- Track convergence: angle between current and final eigenvectors.
- Binarise: sign of eigenvector at each bin → A (+) or B (−). Fix sign so A is the majority-positive label.

### Acceptance criteria

- [ ] Default config: agreement ≥ 95%.
- [ ] Raising noise degrades agreement.
- [ ] Power iteration converges within 20 iterations.
- [ ] Sign-agreement display correct.
- [ ] Panels update live on slider change.
- [ ] Opens pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 9. Artifact #8 — Methylation Landscape Viewer

**File**: `artifacts/lecture-10/08-methylation-landscape.html`
**Lecture anchor**: §1.2 CpG islands and the methylation landscape
**EE framing reinforced**: tissue-specific regulation localised to a small fraction of genomic loci (the CpG-island promoters).

### Teaching purpose

Show the hierarchical methylation pattern — dense unmethylated CpGs at CpG islands, sparse methylated CpGs everywhere else — and make it tissue-specific. Student flips between tissue presets and watches the promoter-methylation switch change while gene-body methylation stays nearly constant. A predicted expression state is derived from promoter methylation to close the loop.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Region: 50 kb window around GENE_X                          │
│ Tissue preset: [ active (liver) ▾ / silenced (lung) /        │
│                  partial (kidney) ]                          │
│ CpG-island parameters:                                      │
│   min length (bp): [──●── 200 ]                             │
│   min observed/expected CpG: [──●── 0.6 ]                   │
│   min GC content: [──●── 0.50 ]                             │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ SVG panel:                                                  │
│   (a) gene model (exons, introns, TSS arrow)                 │
│   (b) CpG density track (per-100 bp CpG count)              │
│   (c) CpG-island highlight (region passing thresholds)       │
│   (d) per-CpG methylation circles (filled = methylated)      │
├─────────────────────────────────────────────────────────────┤
│ Side readout:                                               │
│   Promoter methylation: 8% (island unmethylated)             │
│   Predicted state: EXPRESSED                                │
│   Switch to silenced tissue → promoter methyl 92% →          │
│   predicted state: SILENCED                                 │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Identified: 1 CpG island at promoter (0.8 kb)         │   │
│ │ Active tissue: island unmethylated → expressed         │   │
│ │ Silenced tissue: island methylated → silenced          │   │
│ │ Gene-body methylation differs &lt; 10% across tissues     │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Tissue preset dropdown: "active (liver)" / "silenced (lung)" / "partial (kidney)" — distinct per-tissue methylation patterns.
- CpG-island detector thresholds: min length (100–1000 bp, default 200), min observed/expected CpG ratio (0.4–0.8, default 0.6), min GC content (0.40–0.60, default 0.50).
- Re-simulate button (regenerates slight stochastic variation within the current tissue).

### What they see

- SVG with four horizontally-aligned tracks:
  - Gene model (exons, introns, TSS arrow).
  - CpG density curve (per-100-bp CpG count).
  - CpG-island highlight (regions passing all three threshold criteria, shaded).
  - Per-CpG methylation status (filled circle = methylated; open = unmethylated), positioned at the genomic coordinate.
- Right-side readout showing the promoter methylation fraction and a derived "expression state" based on a simple rule (promoter CpG-island methylation > 50% → SILENCED; else EXPRESSED).
- Outcome banner summarising per-tissue comparison and the finding that variation concentrates at the island.

### Target aha

Switch between active / silenced tissues. The gene-body methylation barely moves; the promoter CpG island flips from mostly-unmethylated to mostly-methylated, and the expression-state readout flips with it. Student sees that tissue-specific regulation localises to the CpG island while most of the landscape is structurally stable.

### Technical notes

- Pure JS, seeded mulberry32.
- Simulate a 50 kb region with a planted gene (TSS at ~15 kb, 3 exons, ending at ~25 kb).
- CpG positions: Poisson process with rate ~0.015 per bp outside the island, boosted to 0.08 per bp inside the island (~500 bp region overlapping the promoter).
- Methylation state per CpG per tissue:
  - Active tissue: island CpGs ~5% methylated, shore ~30%, gene body / intergenic ~80%.
  - Silenced tissue: island CpGs ~90% methylated, shore ~60%, gene body / intergenic ~82%.
  - Partial tissue: island CpGs ~50% methylated, shore ~45%, gene body / intergenic ~80%.
- Draw each per-CpG state from Bernoulli with the tissue-specific probability.
- CpG-island detection: slide a 500-bp window; compute length, obs/exp CpG ratio, GC% per window; call runs where all three thresholds are met. Highlight the detected island.
- Expression readout rule: mean methylation over CpGs inside the detected island → if > 50%, "SILENCED"; else "EXPRESSED".

### Acceptance criteria

- [ ] Active tissue preset → island unmethylated → state EXPRESSED.
- [ ] Silenced tissue preset → island methylated → state SILENCED.
- [ ] Partial tissue → intermediate island methylation → state shown with uncertainty.
- [ ] Switching tissues preserves gene-body methylation pattern approximately.
- [ ] CpG-island detector correctly identifies the planted island.
- [ ] Adjusting thresholds may cause the island to fall out of detection.
- [ ] Opens with active-tissue preset pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 10. Cross-Artifact Consistency

- All eight artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Contact-matrix visualisations use the same `--accent` gradient colour scheme.
- Methylation-state markers (filled/open circles) consistent across Artifacts #1, #2, and #8.
- Eigenvector sign shading: `--accent` for +, `--warning` for − (Artifact #7 and Figure 11).
- Every artifact emits an **outcome banner** per convention §1.4.

## 11. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls.
- [ ] `<script src="../_shared/resize.js" defer></script>` embedded near `</body>`.
- [ ] Outcome banner or equivalent verdict line visible at end of any user interaction.
- [ ] User-input artifacts pre-flight inputs with explicit pass/fail messaging.
