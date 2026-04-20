# Lecture 6 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 6 (Differential Expression and Count Statistics).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-06.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

These conventions apply to every artifact in this lecture. Per-artifact sections below override them only when they need to.

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-06/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript. External libraries only if justified in the per-artifact section.
- The file must render standalone when opened directly in a browser.
- Artifact is embedded in the lecture page via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` near the end of `<body>`.** Verified by the C6 smoke check; missing this silently breaks iframe auto-resize.

### 1.2 Visual design

- Use the design tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- No nucleotide colours in this lecture — DE content is statistical, not sequence-level.
- Typography: **Inter** for UI chrome; **JetBrains Mono** for counts, p-values, coefficients, formulas.
- Default state is instructive: the artifact opens showing a meaningful example, no user input required.
- Controls grouped in a panel above or to the left of the visualization.
- No animations longer than ~400 ms. Motion only when it carries information.

### 1.3 Interaction model

- **Input** — count matrices, gene lists, thresholds, dispersion values — editable text fields or sliders, validated.
- **Step / Play / Pause / Reset** — where the artifact shows an iterative or sequential process (BH procedure, GSEA running sum).
- **Speed** — optional slider, 0.25×–4×. Default 1×.
- Illegal input shows a quiet inline message (`--fg-muted`), not a modal.

### 1.4 Explicit outcome reporting (required)

Every artifact in this lecture answers its own question at the end. Concretely:

- If the artifact fits a model, it shows the **fitted parameters** with standard errors and a verdict ("✓ coefficient β̂₁ = 1.15 · z = 3.2 · p = 0.001").
- If the artifact computes a test statistic, it shows the **value and the reference distribution** with the tail.
- If the artifact runs a multiple-testing adjustment, it shows **before-and-after**: raw rejections at α, adjusted rejections at the same α.
- If the artifact computes an enrichment score, it shows **the score and its permutation-null p-value**.

### 1.5 Feasibility gate on user input (required where input is free-form)

Artifacts accepting user input (count matrices, p-value lists, gene rankings) must pre-flight the input and report *why* it is or isn't runnable before the computation runs. Rejected inputs get an inline explanation.

### 1.6 Pedagogical constraint

Every artifact must produce a **specific realization** — the "target aha moment" named in its section. If the student plays with the artifact and doesn't land on that realization, the artifact has failed and should be revised.

### 1.7 Out of scope

- No logins, accounts, or persistence between sessions.
- No telemetry or analytics.
- No external data files larger than ~50 KB.

---

## 2. Artifact #1 — Volcano Plot Explorer

**File**: `artifacts/lecture-06/01-volcano-explorer.html`
**Lecture anchor**: §1.1 What we're asking
**EE framing reinforced**: detection thresholds on a 2D (effect × confidence) map.

### Teaching purpose

Turn the abstract "DE gene list" into an interactive volcano. Student drags thresholds, watches the hit count update, and gets a visual feel for the trade between effect size and significance cuts.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Dataset: [Preset ▾] [Paste your own DE results]             │
│ Preset: [ airway smooth-muscle + dex · HCC vs normal liver  │
│          · simulated null · simulated heavy-signal ▾ ]      │
├─────────────────────────────────────────────────────────────┤
│ Thresholds:                                                 │
│   padj  ≤ [──●── 0.05]        (drag to 0.01, 0.1, ...)       │
│   |log₂FC| ≥ [──●── 1.0]      (drag to 0.5, 2.0, ...)        │
│   Highlight only PASS:  [✓]                                 │
├─────────────────────────────────────────────────────────────┤
│ Volcano plot (SVG):                                         │
│   • horizontal dashed line at y = -log10(padj-threshold)    │
│   • vertical dashed lines at ±log2FC threshold              │
│   • points coloured by significance quadrant                │
│   • highlighted genes labeled with gene symbols             │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 487 genes meet both thresholds · 312 up, 175 down     │   │
│ │ Most extreme up: IL6 (log₂FC 4.2, padj 1e-12)         │   │
│ │ Most extreme down: FOXO3 (log₂FC -3.1, padj 5e-9)     │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Dataset preset**: 4 scenarios.
- **padj threshold slider**: 0.001 to 0.5 (log scale).
- **log₂FC threshold slider**: 0 to 4.
- **Highlight-only-PASS toggle**: greys out non-passing points to emphasise the hit set.
- **Paste your own**: textarea accepting a DE table (CSV: gene, log₂FC, padj).

### What they see

- SVG scatter of ~2000 gene points.
- Dashed threshold lines that move as the sliders change.
- Points coloured by quadrant (up-significant = `--accent`, down-significant = `--error`, not-significant = `--fg-muted`).
- Up to 20 top-significant genes labeled with gene symbols in JetBrains Mono.
- Hit count + most-extreme genes in the outcome banner.

### Target aha moment

Start with the default (airway dex). Tighten padj to 0.001 — the hit count drops from ~500 to ~100. Now raise log₂FC to 2 — down to ~30. Student sees that each threshold independently reshapes the gene list, and "significant" is an interval of choice.

### Technical notes

- Pure JS. Presets hardcode ~2000 gene rows.
- Points drawn as SVG circles; thresholds as dashed lines; re-render on slider change (debounced via requestAnimationFrame).
- Feasibility gate: pasted data parsed; malformed rows reported with line numbers.

### Acceptance criteria

- [ ] Default preset shows ~500 hits at the default thresholds.
- [ ] Moving padj and log₂FC sliders updates the scatter live.
- [ ] At least 4 preset scenarios, each with distinct signal distributions.
- [ ] Custom input with malformed rows shows inline errors.
- [ ] Opens standalone at default preset with volcano already drawn.
- [ ] HTML parses; JS passes `node --check`; contains exactly 1 `_shared/resize.js`.

---

## 3. Artifact #2 — Dispersion Shrinkage Visualizer

**File**: `artifacts/lecture-06/02-dispersion-shrinkage.html`
**Lecture anchor**: §3.3 Empirical-Bayes shrinkage toward the fitted trend
**EE framing reinforced**: regularisation toward a learned prior.

### Teaching purpose

Let the student watch empirical-Bayes shrinkage operate on a simulated dispersion cloud. Change the shrinkage strength; watch gene-wise estimates move toward the fitted trend. See why shrinkage calibrates p-values.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Scenario: [Preset ▾ · realistic RNA-seq]                    │
│ Number of genes: [──●── 500 ] (more = slower)               │
│ Shrinkage strength: [──●── 0.6 ] (0 = no shrinkage)          │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Log-log scatter: mean (x) vs dispersion (y)                 │
│   • raw MLE circles (unfilled, grey)                        │
│   • shrunk estimates (filled, accent)                       │
│   • arrows from raw to shrunk                               │
│   • fitted trend curve overlaid                             │
│                                                             │
│ Click a gene to see its raw, shrunk, and trend values       │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 500 genes · shrinkage strength 0.6                    │   │
│ │ Average shrinkage: raw SD 0.48 → shrunk SD 0.15       │   │
│ │ Genes >3σ from trend: 8 (resist shrinkage most)       │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Scenario preset**: "realistic RNA-seq", "high-noise (small n)", "sparse signal".
- **Number of genes slider**: 100 to 2000.
- **Shrinkage strength slider**: 0 (no shrinkage, raw MLEs) to 1 (full shrinkage to trend).
- **Re-simulate** button: generates a new sample.

### What they see

- Log-log scatter with raw-MLE circles, shrunk dots, connecting arrows, and the fitted trend.
- Live update on slider changes (some debouncing for re-simulate).
- Statistics: how much each gene moved, the average shrinkage magnitude, and how many "outliers" resist.

### Target aha moment

Slide shrinkage strength from 0 to 1. Genes cluster onto the trend; the scatter tightens dramatically. But a handful of "outlier" genes (far from trend) resist — they don't snap to the line. Student sees that shrinkage is informed compromise, not blind pull-to-mean: outliers stay outliers.

### Technical notes

- Pure JS. Simulate dispersions from a log-normal distribution centred on a parametric trend function. Add per-gene noise.
- Shrinkage formula: `shrunk = raw * (1 - s) + trend * s` where s is the strength slider.
- SVG scatter with live redraw on slider change.

### Acceptance criteria

- [ ] Default scenario shows visible shrinkage trajectories.
- [ ] Strength slider 0 → unchanged raw; 1 → all points on trend.
- [ ] Outlier genes resist shrinkage (stay far from trend even at high strength).
- [ ] Re-simulate produces visibly different but consistent distributions.
- [ ] Opens standalone with default scenario pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — NB GLM Fitter

**File**: `artifacts/lecture-06/03-nb-glm-fitter.html`
**Lecture anchor**: §3.4 Hypothesis tests — Wald and LRT
**EE framing reinforced**: maximum-likelihood inference with a log-link.

### Teaching purpose

Fit the NB GLM to a toy 6-sample, 2-condition dataset. Show β̂₀, β̂₁, their standard errors. Compute Wald and LRT p-values. Change counts and watch the fit update.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Toy dataset (editable):                                     │
│   Sample  Condition  Size factor  Count                     │
│   S1      control    1.00         12                        │
│   S2      control    1.05         10                        │
│   S3      control    0.95         14                        │
│   S4      treated    1.02         42                        │
│   S5      treated    0.98         38                        │
│   S6      treated    1.01         45                        │
│ Dispersion α: [──●── 0.05 ]                                  │
│ [Refit GLM]                                                 │
├─────────────────────────────────────────────────────────────┤
│ Fitted model:                                               │
│   log μ = log s + β₀ + β₁ · condition                       │
│   β̂₀ = 2.53  SE 0.16                                         │
│   β̂₁ = 1.26  SE 0.21    (log2 FC = 1.82)                    │
│                                                             │
│ Tests:                                                      │
│   Wald:  z = 1.26 / 0.21 = 6.0 → p = 2.0e-9                 │
│   LRT:   λ = 2·(ℓ_full − ℓ_red) = 35.2 → p = 3.0e-9         │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ ✓ β̂₁ = 1.26 (log2 FC = 1.82, ~3.5× up)                │   │
│ │ Wald p = 2.0e-9 · LRT p = 3.0e-9                     │   │
│ │ Strong evidence of upregulation.                      │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Editable count matrix**: six cells, integer inputs.
- **Size factor slider per sample** (optional): show offset effect.
- **Dispersion α slider**: 0.001 to 1.0.
- **Refit GLM** button.

### What they see

- Fitted coefficients with standard errors.
- Both Wald and LRT p-values side by side, usually agreeing.
- Log2 fold change and a one-line verdict.

### Target aha moment

Start with the default — clear up-regulation. Reduce the treated counts so they're close to control counts. The β̂₁ shrinks, Wald and LRT p-values rise above 0.05, and the verdict flips from "✓ upregulation" to "✗ no significant difference." Now crank the dispersion α up to 0.5 — even with clear mean difference, the p-value becomes non-significant because the noise model now allows big mean differences by chance.

### Technical notes

- Pure JS. Implement NB-GLM fitting via iteratively reweighted least squares (IRLS) with the specified dispersion held fixed. Full fit + reduced fit for LRT.
- Log-likelihood computation for NB: use `lgamma` approximation (Stirling).
- Feasibility gate: non-integer / negative counts rejected with inline error.

### Acceptance criteria

- [ ] Default preset produces clear upregulation with Wald p &lt; 0.001.
- [ ] Editing counts updates fit live (or on Refit click).
- [ ] Wald and LRT p-values agree within reasonable tolerance (~10% for this sample size).
- [ ] Dispersion slider visibly changes p-value without changing the mean estimate much.
- [ ] Opens with default preset already fit.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — FDR Simulator

**File**: `artifacts/lecture-06/04-fdr-simulator.html`
**Lecture anchor**: §4.3 FDR interpretation and q-values
**EE framing reinforced**: detection with controlled false-alarm rate.

### Teaching purpose

Simulate 20,000 gene tests with a controlled fraction of true positives. Run Bonferroni and Benjamini–Hochberg at multiple thresholds. Report empirical FDR and power. Show how the two procedures trade off differently.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulation parameters:                                      │
│   Total genes: [──●── 20000 ]                               │
│   Fraction true positives: [──●── 5% ] = 1000 true DE        │
│   Effect size (log₂FC under H₁): [──●── 1.0 ]                │
│   Sample size per condition: [──●── 3 ]                      │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ p-value histogram (simulated):                              │
│   [ flat null distribution + spike near 0 from TPs ]         │
├─────────────────────────────────────────────────────────────┤
│ Correction: [ Bonferroni │ Benjamini–Hochberg (BH) ]         │
│ Threshold: [──●── 0.05 ]                                     │
│                                                             │
│ Results:                                                    │
│   Rejections: 487                                           │
│   True positives among rejections: 441 (90.6% precision)    │
│   False positives: 46 (9.4% of rejections)                  │
│   Power: 441 / 1000 = 44.1%                                 │
│   Empirical FDR: 9.4% (target ≤ 5.0%)                       │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ BH: 487 rejections · 9.4% false positive rate         │   │
│ │ Would-have been Bonferroni: 128 rejections · 0.2% FPR │   │
│ │ BH gives 3.8× more hits at the cost of 47× FPR        │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Total genes slider**: 1000 to 30000.
- **Fraction true positives slider**: 0 to 50%.
- **Effect size slider**: 0 to 3 log₂FC under H₁.
- **Sample size slider**: 2 to 20 per condition.
- **Correction choice**: Bonferroni or BH.
- **Threshold slider**: 0.001 to 0.2.
- **Re-simulate** button.

### What they see

- p-value histogram live — gets sharper spike with more signal, flatter with less.
- Rejection count, true-positive count, power, empirical FDR — all live.
- Side-by-side comparison of Bonferroni and BH at the same threshold.

### Target aha moment

Start with the default (5% true positives, effect 1, n = 3). Toggle Bonferroni vs BH at threshold = 0.05. Bonferroni finds ~100 genes; BH finds ~500. Power goes from 10% to 50%. Empirical FDR goes from 0.1% to 5%. Student sees BH's trade: 5× more power for a specific, controlled increase in false-positive rate.

### Technical notes

- Pure JS. Simulate p-values: null genes uniform, TP genes from beta(a, 1) distribution biased toward 0 based on effect size and sample size.
- Bonferroni: reject if p &lt; α / m.
- BH: standard Benjamini-Hochberg, implemented in ~10 lines.
- Histogram: 50-bin fixed.
- Feasibility gate: none needed (sliders constrained).

### Acceptance criteria

- [ ] Default shows clear difference between Bonferroni (~100) and BH (~500) rejections.
- [ ] Empirical FDR matches the BH threshold (within ~20%).
- [ ] Zero-signal scenario (0% TPs) produces ~0 hits from both corrections.
- [ ] Re-simulate produces visibly different but consistent results.
- [ ] Opens with default scenario rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — ORA Contingency Explorer

**File**: `artifacts/lecture-06/05-ora-contingency.html`
**Lecture anchor**: §5.1 Over-representation analysis (ORA)
**EE framing reinforced**: the hypergeometric test as detection of population enrichment.

### Teaching purpose

Compute Fisher's exact test for gene-set enrichment, showing the 2×2 contingency table and the resulting hypergeometric p-value.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Parameters:                                                 │
│   DE gene list size:   [──●── 500 ]                          │
│   Gene-set size:       [──●── 127 ] (G2/M cell cycle)        │
│   Overlap (DE ∩ set):  [──●── 42 ]                           │
│   Background (total):  [──●── 20000 ]                        │
│ [Compute]                                                   │
├─────────────────────────────────────────────────────────────┤
│ Contingency table:                                          │
│             in set G  not in G    total                     │
│   DE           42       458        500                      │
│   not DE       85      19415     19500                      │
│   total       127      19873     20000                      │
│                                                             │
│ Observed vs expected:                                       │
│   Observed DE in G: 42                                      │
│   Expected by chance: 500 × 127/20000 = 3.18                │
│   Ratio: 13.2× enrichment                                   │
├─────────────────────────────────────────────────────────────┤
│ Test results:                                               │
│   Fisher's exact p (one-sided, over): 3.2e-8                │
│   Chi-square p:                       2.0e-8                │
│   Odds ratio:                         13.4                  │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ ✓ Significant enrichment                              │   │
│ │ G2/M genes are 13× over-represented in DE list        │   │
│ │ Fisher's exact p = 3.2e-8                             │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **DE gene list size slider**: 10 to 2000.
- **Gene-set size slider**: 10 to 1000.
- **Overlap slider**: 0 to min(DE, set).
- **Background slider**: 10000 to 50000.
- **Compute** button (or live update).

### What they see

- Live-updated 2×2 table with all margins.
- Observed / expected ratio.
- Fisher's exact p-value, chi-square p-value (for comparison), odds ratio.
- Side-by-side: a hypergeometric PMF plot with the observed value marked.

### Target aha moment

Start with the default (42 observed, 3.18 expected, 13× enrichment, p ≈ 1e-8). Drop the overlap to 5 — now expected is still ~3, observed is 5, barely enriched, p rises to 0.2. Student sees that "enrichment" is relative to what's expected by chance, and small absolute numbers with small expected values are surprisingly insignificant.

### Technical notes

- Pure JS. Fisher's exact: sum over hypergeometric PMF from observed to max. Use `lgamma` for numerical stability.
- Chi-square: `(O - E)² / E` summed over the 4 cells, compared to χ²(1).
- Feasibility gate: overlap must be ≤ min(DE, set); inline error otherwise.

### Acceptance criteria

- [ ] Default produces p ≈ 1e-8 with 13× enrichment.
- [ ] Equal-to-expected overlap produces p ≈ 1.
- [ ] Small overlap (2-3 genes) with small expected produces intermediate p.
- [ ] Sliders update the table and p-values live.
- [ ] Opens with default parameters computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — GSEA Running-Sum Simulator

**File**: `artifacts/lecture-06/06-gsea-running-sum.html`
**Lecture anchor**: §5.2 Gene-set enrichment analysis (GSEA)
**EE framing reinforced**: CUSUM-style change-point detection on a ranked list.

### Teaching purpose

Let the student watch the GSEA running-sum statistic walk down a ranked gene list, see the enrichment score emerge at the peak, and compare to a permutation null.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Scenario: [Preset ▾ · concentrated gene set at top]          │
│ Rank list length: [──●── 2000 ]                              │
│ Gene set size: [──●── 40 ]                                   │
│ Enrichment at top: [──●── 70% ] (of set genes in top 20%)    │
│ [Re-shuffle] [Step along] [Run full walk]                    │
│ Step: 0 / 2000                                              │
├─────────────────────────────────────────────────────────────┤
│ Running sum plot:                                           │
│   [ curve walking +/- along rank axis ]                     │
│   current position marked                                   │
│ Ranked gene list (bottom strip): set-genes as ticks          │
├─────────────────────────────────────────────────────────────┤
│ Enrichment score: ES = 0.48 (peak at rank 283)              │
│ Permutation null (N=1000 shuffles): observed > 99.2%        │
│ p-value: 0.008                                              │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ ✓ Significant enrichment                              │   │
│ │ Set of 40 genes enriched at top of ranking            │   │
│ │ ES = 0.48 · permutation p = 0.008                    │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Scenario preset**: "concentrated top", "spread", "anti-correlated (concentrated at bottom)", "random".
- **Rank list length slider**: 100 to 5000.
- **Gene set size slider**: 5 to 200.
- **Enrichment-at-top slider**: 0% (random placement) to 100% (all set genes at top).
- **Re-shuffle** button.
- **Step** button: advance one gene in the walk.
- **Run full walk** button: animate to completion.

### What they see

- A ranked gene list shown as a horizontal strip with set-genes marked as ticks.
- A running-sum curve above, animating upward when a set-gene is encountered, downward otherwise.
- The enrichment score (peak height) labeled.
- A permutation null histogram showing where the observed ES falls.

### Target aha moment

Start with the default "concentrated top" preset. Watch the running sum climb steeply in the first 200 genes (dense with set-genes), then slowly decline. Peak at ES = 0.48. Switch to "spread" preset — running sum barely deviates from zero, ES ≈ 0.1, p ≈ 0.3. Switch to "anti-correlated" — running sum goes negative (enrichment at bottom), ES of opposite sign. Student sees GSEA isn't just "count genes in set" — it's "how concentrated are they at the top or bottom?"

### Technical notes

- Pure JS. Ranked list simulated based on scenario. Running sum updated step-by-step with weighted increment (per Subramanian 2005 formula using signed scores).
- Permutation null: shuffle set membership labels, recompute ES, repeat 1000 times.
- Feasibility gate: gene set size must be positive; list length must be positive.

### Acceptance criteria

- [ ] "Concentrated top" produces ES ≈ 0.4-0.6, p &lt; 0.01.
- [ ] "Random" produces ES ≈ 0.1, p &gt; 0.1.
- [ ] "Anti-correlated" produces negative ES.
- [ ] Step button advances the walk one gene at a time; run-full animates to completion.
- [ ] Permutation null histogram updates on re-simulate.
- [ ] Opens with default preset showing initial state.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Cross-Artifact Consistency

- All six artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- p-values and z-scores are rendered in JetBrains Mono consistently; gene symbols also in JetBrains Mono.
- The Volcano Explorer (#1) and FDR Simulator (#4) should use the same threshold-line conventions.
- The Dispersion Shrinkage Visualizer (#2) and GSEA Running-Sum Simulator (#6) share the scatter-plus-curve visual vocabulary.
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
