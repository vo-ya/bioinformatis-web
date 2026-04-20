# Lecture 6 — Figures Specification

> **Scope**: Static diagrams for Lecture 6 (Differential Expression and Count Statistics).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow the parent `diagram-style-guide.md` for all visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-06.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG. DE content is statistical and algorithmic; no photographs.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure must be legible at 720 px and scale cleanly up to 1200 px.
- Line-art first; fill only where the fill earns the ink.
- Base colours (`--base-a/t/g/c`) are unused this lecture — no nucleotide characters are rendered.
- Variable names, coefficients, and numeric values use **JetBrains Mono**; all other labels use **Inter**.
- Arrows follow the shared `<marker id="arrow-accent">` pattern.
- Escape `&`, `<`, `>` as XML entities in text content.

## Figure Budget

Eleven figures for a ~3h 30min lecture. Placement by part:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | Anatomy of a volcano plot | Part 1 | Custom SVG |
| 2 | Why t-tests fail on counts | Part 1 | Custom SVG |
| 3 | Mean–dispersion relationship | Part 2 | Custom SVG |
| 4 | NB GLM — design matrix and link | Part 2 | Custom SVG |
| 5 | Per-gene dispersion MLE on few samples | Part 3 | Custom SVG |
| 6 | Empirical-Bayes dispersion shrinkage | Part 3 | Custom SVG |
| 7 | Wald vs LRT — when each is right | Part 3 | Custom SVG |
| 8 | The p-value histogram diagnostic | Part 4 | Custom SVG |
| 9 | The Benjamini–Hochberg procedure | Part 4 | Custom SVG |
| 10 | ORA — Fisher's exact | Part 5 | Custom SVG |
| 11 | GSEA — the running-sum statistic | Part 5 | Custom SVG |

---

## Figure 1 — Anatomy of a volcano plot

**File**: `diagrams/lecture-06/01-volcano-anatomy.svg`
**Lecture anchor**: §1.1 What we're asking
**ViewBox**: `0 0 900 460`

### Purpose

Set the lecture's visual language in one image. Every DE result is summarised as a volcano; readers should read one at a glance by the end of the lecture.

### Content

- Standard volcano plot: x = log₂ fold change (range roughly −5 to 5), y = −log₁₀(padj) (range 0 to ~15).
- Scatter of ~2,000 gene points in `--fg-muted` at low-significance / low-effect, darkening to `--fg` and `--accent` at significant positions.
- Vertical dashed lines at log₂FC = ±1; horizontal dashed line at padj = 0.05 (y ≈ 1.3). All three in `--accent` 1px dashed.
- Label the four quadrants: "significant up-regulated", "not significant up", "not significant down", "significant down-regulated".
- Highlight 5 named example genes with labels and callouts (e.g., "IL6 · log₂FC 3.2, padj 10⁻⁸", "MKI67 · log₂FC 1.1, padj 0.003"). Accent-bg halo around each.
- Small panel at top right: "significance threshold: padj &lt; 0.05 · effect threshold: |log₂FC| &gt; 1".

### Style notes

- Axes: `--fg-muted` 1px stroke.
- Points: 2px radius circles with graded opacity by significance.
- Highlighted genes: `--accent` 4px radius circles, labeled in JetBrains Mono 10.
- Threshold lines: `--accent` 1px dashed with small labels "padj = 0.05", "log₂FC = ±1".

---

## Figure 2 — Why t-tests fail on counts

**File**: `diagrams/lecture-06/02-t-test-failure.svg`
**Lecture anchor**: §1.2 Why t-tests fail on counts
**ViewBox**: `0 0 900 400`

### Purpose

Show in one picture both why count data violates t-test assumptions and what the downstream p-value cost is.

### Content

**Left panel — count histogram vs normal.** A histogram of a low-count gene (mean ~5, heavy right skew, spike at 0) in solid `--accent` fill. Overlaid with a normal density of the same mean and variance in `--fg-muted` dashed line. Visible mismatch: the normal goes below zero (impossible for counts), the data has a zero-spike (impossible for normal). Annotations point out both problems.

**Right panel — p-value distribution under the null.** The empirical p-value histogram from running a t-test per gene on a dataset where every gene is truly null. Should be uniform. It isn't — there's a hump near 0 (false positives) and a deficit near 1. Annotation: "under the null, p-values should be uniform — this t-test produces p &lt; 0.05 ~12% of the time, not 5%."

### Style notes

- Histograms in `--accent-bg` fill with `--accent` stroke.
- Normal overlay: `--fg-muted` dashed 1.5px.
- Annotations in Inter 10 with `--warning` arrows.
- Footer caption: "The t-test assumes normal data with constant variance. Counts are neither."

---

## Figure 3 — Mean–dispersion relationship

**File**: `diagrams/lecture-06/03-mean-dispersion.svg`
**Lecture anchor**: §2.1 The NB model in more depth
**ViewBox**: `0 0 900 420`

### Purpose

Show the empirical trend that makes empirical-Bayes shrinkage (§3.3) possible. Points at low mean have high dispersion; points at high mean have low dispersion; the relationship is noisy but structured.

### Content

- Log-log scatter: x = mean count per gene (10⁰ to 10⁴), y = per-gene dispersion α (10⁻³ to 10¹).
- ~2,000 gene points in `--fg` at 1px radius.
- Overlaid fitted trend line — smooth decreasing curve from top-left to bottom-right — in `--accent` 2px.
- Annotations at three example points: "low-count gene, high dispersion", "medium-expression gene, moderate dispersion", "highly-expressed gene, tight dispersion".
- Marginal histograms along each axis (thin) showing the distribution of mean and dispersion values.

### Style notes

- Axes: `--fg-muted` with `--fg-subtle` grid lines.
- Points: 1px `--fg` circles.
- Fitted trend: `--accent` 2px solid.
- Annotations in Inter 10.

---

## Figure 4 — NB GLM, design matrix and link

**File**: `diagrams/lecture-06/04-nb-glm-design.svg`
**Lecture anchor**: §2.2 The NB generalised linear model
**ViewBox**: `0 0 980 440`

### Purpose

Show the GLM's moving pieces concretely: a design matrix, a gene's counts, the log-link equation, and the fitted coefficients. Readers should be able to set up the model for a new experimental design after seeing this.

### Content

**Left — sample table.** 6 samples with columns: Sample ID (S1–S6), Condition (control/treat), Batch (A/B), Size factor (0.9–1.1). The first 3 are control (Batch A+B), the last 3 are treated (Batch A+B).

**Middle — design matrix.** A matrix with rows = samples, columns = (Intercept, Condition, Batch). Values 0/1 per cell in JetBrains Mono, with column headers in `--accent`.

**Right — one gene's fit.** Gene counts (one column, 6 rows). Arrow down to: the GLM equation: `log μ = log s + β₀ + β₁·Condition + β₂·Batch` (in JetBrains Mono). Below the equation: the fitted coefficients β̂ = (4.2, 0.8, -0.1) with standard errors. At the bottom: log₂FC = β̂₁ / log(2) = 1.15.

### Style notes

- Design matrix: `--bg-muted` fill, `--border-strong` stroke, values in JetBrains Mono 11.
- GLM equation: emphasized box with `--accent-bg` fill.
- Arrows between panels: `--accent` 1.5px with arrowheads.

---

## Figure 5 — Per-gene dispersion MLE on few samples

**File**: `diagrams/lecture-06/05-dispersion-mle.svg`
**Lecture anchor**: §3.2 Per-gene dispersion estimation
**ViewBox**: `0 0 900 440`

### Purpose

Make the "3-replicate MLE is terrible" point visceral. Side-by-side: the log-likelihood surface for a single gene with 3 vs 20 replicates.

### Content

**Left — 3 replicates per condition.** Six data points shown as a small strip plot (3 control counts, 3 treated). The NB log-likelihood as a function of α, plotted from α = 10⁻³ to 10¹. The curve is very broad, almost flat — essentially no information about the true α. The MLE (maximum of the curve) is marked, with a 95% CI bracket spanning more than two orders of magnitude.

**Right — 20 replicates per condition.** Forty data points. The same log-likelihood plot. Now the curve has a clearly visible maximum with a CI spanning ~30% of the MLE value. Much more information.

**Bottom caption**: "Per-gene MLE dispersion from 3 replicates is statistically unreliable. The next step — empirical-Bayes shrinkage — fixes this by sharing information across the other 19,999 genes in the dataset."

### Style notes

- Strip plots: small circles per sample, coloured by condition.
- Log-likelihood curves: `--accent` 2px solid.
- MLE marker: `--accent` 4px circle.
- CI bracket: horizontal line at MLE with whiskers to CI endpoints.

---

## Figure 6 — Empirical-Bayes dispersion shrinkage

**File**: `diagrams/lecture-06/06-eb-shrinkage.svg`
**Lecture anchor**: §3.3 Empirical-Bayes shrinkage toward the fitted trend
**ViewBox**: `0 0 900 440`

### Purpose

Show the central algorithmic idea of DESeq2/edgeR: per-gene dispersion is shrunk toward the fitted mean-dispersion trend. Genes near the trend get strongly shrunk; genes far from the trend resist.

### Content

- The same mean-dispersion scatter as Figure 3, but with two layers:
  - **Raw MLEs**: unfilled circles in `--border-strong`, ~50 selected example genes.
  - **Shrunk estimates**: filled `--accent` circles, same 50 genes.
  - Each raw MLE is connected to its shrunk counterpart by a thin `--fg-muted` arrow — showing the shrinkage trajectory.
- The fitted trend passes through the cloud in `--accent` 2px.
- Highlight three genes by name:
  - Gene A (near the trend, short arrow to trend) — "strongly shrunk, trend-consistent"
  - Gene B (far above trend, medium arrow) — "moderately shrunk, partially resisted"
  - Gene C (very far above trend, short arrow) — "almost no shrinkage, likely real outlier"

### Style notes

- Raw MLE circles: unfilled with `--border-strong` 1.5px stroke.
- Shrunk circles: filled `--accent`.
- Shrinkage arrows: `--fg-muted` 1px.
- Named genes: labeled in JetBrains Mono 10, positioned near their shrunk circles.

---

## Figure 7 — Wald vs LRT: when each is right

**File**: `diagrams/lecture-06/07-wald-vs-lrt.svg`
**Lecture anchor**: §3.4 Hypothesis tests — Wald and LRT
**ViewBox**: `0 0 900 440`

### Purpose

Side-by-side geometric intuition for the two tests, with a small decision table for which to pick.

### Content

**Top half — Wald test.** A normal density with x-axis labeled "β̂ / SE(β̂) — z-score". Marked point at z = 2.4 with shading under the tails beyond ±2.4; p-value = 0.016 annotated. Caption: "asks: is β̂ far from 0 relative to its uncertainty?"

**Bottom half — LRT.** Two overlapping log-likelihood surfaces (full model with β₁, reduced without) with maxima labeled. The vertical difference λ = 2·(ℓ_full − ℓ_reduced) is shown as a double-arrow. The chi-square reference distribution on the right, with the observed λ value marked. p-value = 0.014 annotated. Caption: "asks: does removing β₁ materially hurt the fit?"

**Decision table (right side, small)**:
| Hypothesis | Preferred test |
|---|---|
| 2-condition comparison | Wald |
| 3+ level factor (ANOVA-like) | LRT |
| Time-course interaction | LRT |
| Small n, nonlinear GLM | LRT |

### Style notes

- Normal curve: `--fg` stroke.
- Chi-square curve: `--fg` stroke.
- Shaded tails: `--accent-bg` fill.
- Decision table: `--bg-muted` with `--border-strong` grid.

---

## Figure 8 — The p-value histogram diagnostic

**File**: `diagrams/lecture-06/08-pvalue-histogram.svg`
**Lecture anchor**: §4.1 The problem at 20,000 tests
**ViewBox**: `0 0 900 440`

### Purpose

Three-panel diagnostic reference: what p-value histograms look like for (a) healthy DE, (b) no signal, (c) miscalibrated model.

### Content

Three side-by-side histograms, each labeled:

**(a) Healthy DE result.** Mostly flat bulk (null genes) with a sharp spike near 0 (true-positive genes). ~6% of genes in the leftmost bin. Caption: "what a real analysis should look like".

**(b) No signal.** Entirely flat, uniform distribution. Caption: "every gene is null — expected under H₀".

**(c) Miscalibrated model.** U-shape with a spike near 0 AND a hump near 1. Caption: "something is wrong — model misspecification, outliers, batch effects". Red `--warning` tint on the bars to emphasise "bad".

### Style notes

- Bars: `--accent-bg` fill, `--accent` stroke (for a, b); `--warning-bg` fill, `--warning` stroke (for c).
- Each panel with axis labels "p-value" (x, 0 to 1) and "count" (y).
- Callout: "plot this first, before reporting any results".

---

## Figure 9 — The Benjamini–Hochberg procedure

**File**: `diagrams/lecture-06/09-bh-procedure.svg`
**Lecture anchor**: §4.2 Bonferroni vs Benjamini–Hochberg
**ViewBox**: `0 0 900 440`

### Purpose

Visualise the BH procedure as a line-crossing operation on sorted p-values, so readers can run it by hand.

### Content

Two panels:

**Left — high-signal scenario.** x = rank k (1 to 500, log scale), y = p-value. Sorted p-values plotted as dots. The BH line `p = (k/m) · α` overlaid for α = 0.05 — a straight line from (1, 5·10⁻⁶) to (m=20000, 0.05). The largest k where p₍ₖ₎ is below the line is marked with a vertical drop line; k* ≈ 450 here. Rejections = 450, expected false positives = 450·0.05 = ~22.

**Right — low-signal scenario.** Same axes. The p-values now sit near or above the BH line for all ranks. Only k* ≈ 12 below the line. Rejections = 12. A lot less power.

Side annotation: "BH interpretation: if I reject all genes with padj ≤ 0.05, ~5% of my rejections are false positives in expectation."

### Style notes

- Sorted p-values: `--fg` 2px circles.
- BH line: `--accent` 2px solid.
- Crossover point: `--accent` vertical drop line, labeled "k*".
- Log scales on both axes.

---

## Figure 10 — ORA: Fisher's exact

**File**: `diagrams/lecture-06/10-ora-fisher.svg`
**Lecture anchor**: §5.1 Over-representation analysis
**ViewBox**: `0 0 900 400`

### Purpose

Walk through the ORA logic with a concrete example so readers can compute enrichment by hand.

### Content

**Left — 2×2 contingency table.** Labeled rows: "DE gene" / "Not DE". Labeled columns: "in gene set G2/M" / "not in G2/M". Numbers filled in: (42, 458, 500) / (85, 19415, 19500). Margins (column sums) = (127, 19873). Total = 20000.

**Middle — the hypergeometric intuition.** "If you draw 500 genes at random from a bag of 20,000 (127 white, 19,873 black), what's the probability of getting ≥42 white?" A small hypergeometric PMF plot with the tail beyond 42 shaded.

**Right — the Fisher's exact p-value.** Numeric result: p = 3.2 × 10⁻⁸. Small interpretation: "DE genes are ~6.6× enriched for G2/M. Fisher says this would happen &lt; 1 in 30 million times by chance."

### Style notes

- Table: `--bg-muted` fill, `--border-strong` grid.
- Hypergeometric PMF: `--accent-bg` bars with tail shaded `--warning-bg`.
- Numbers in JetBrains Mono 11.
- p-value in large JetBrains Mono 14.

---

## Figure 11 — GSEA: the running-sum statistic

**File**: `diagrams/lecture-06/11-gsea-running-sum.svg`
**Lecture anchor**: §5.2 Gene-set enrichment analysis
**ViewBox**: `0 0 900 460`

### Purpose

Show GSEA's core mechanism: the running sum walking along a ranked gene list, peaking at the point of maximum enrichment.

### Content

**Bottom band — ranked gene list.** A long horizontal strip showing 20,000 genes ranked left-to-right by some signed DE statistic. Genes in the target gene set G are marked as short vertical ticks along the strip; other genes are white. Visibly more ticks near the left (top of the ranking).

**Middle band — running sum.** A curve that rises when a G-gene is encountered and falls when a non-G gene is. The curve peaks somewhere in the top third (where G-gene density is highest), then declines. The peak height = ES (enrichment score). A horizontal line at zero for reference.

**Top band — enrichment score annotation.** Label the peak "ES = 0.48 · p = 0.003 (permutation-based)".

**Right inset — permutation null.** A small histogram of ES values from 1000 label-shufflings; the observed value marked as an arrow above the right tail.

### Style notes

- Ticks for G-genes: `--accent` 1px vertical lines.
- Running sum curve: `--fg` 2px solid; area above zero shaded `--accent-bg`.
- Peak marker: `--accent` filled circle with label.
- Permutation histogram: `--fg-muted` bars with the observed value arrow in `--accent`.

---

## Cross-Figure Consistency Notes

- **Dispersion scatters** in Figures 3 and 6 use the same log-log axes and the same cloud of ~2,000 points.
- **p-value displays** in Figures 8 and 9 share axis conventions — 0 to 1 for raw p-values, log scale for sorted ranks.
- **Gene-set operations** in Figures 10 and 11 should use the same gene-set colour (`--accent`) and non-set colour (`--fg-muted`).
- **GLM equations** in Figures 4, 5, 7 should use identical notation (μ, s_i, β₀, β₁, α) so that readers can trace variables across figures.

## Pre-Submission Checklist (Lecture-Wide)

- [ ] All eleven figures render standalone in the browser with no external dependencies.
- [ ] No figure uses a gradient, drop shadow, glow, or 3D effect.
- [ ] All numeric values, variables, and coefficients are in JetBrains Mono; all other labels in Inter.
- [ ] Base colours (`--base-*`) are NOT used in this lecture — no nucleotides rendered.
- [ ] Every figure has `role="img"`, `<title>`, and `<desc>`.
- [ ] Every figure is legible at 720 px.
- [ ] Filenames follow `NN-name-kebab.svg` with zero-padded numbering.
- [ ] All `&`, `<`, `>` in text content are XML-escaped.
