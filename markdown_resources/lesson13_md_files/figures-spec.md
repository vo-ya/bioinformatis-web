# Lecture 13 — Figures Specification

> **Scope**: Static diagrams for Lecture 13 (GWAS and statistical genetics).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow `diagram-style-guide.md` for visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-13.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG; content is plot-heavy (Manhattan, QQ, PCA, locus zoom, heatmap).
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure legible at 720 px; scales to 1200 px.
- Plots use standard x/y-axis conventions; use Source Serif italic for math ($p$, $\beta$, $\lambda$, $r^2$, $h^2$).
- Monospace (JetBrains Mono) for: SNP labels (rs…), −log10(p) values, λ_GC, and numerical annotations.
- Inter for UI labels and annotations.
- Ancestry palette consistent with Lecture 12: African cobalt, European amber, East Asian red, South Asian green, Admixed violet.
- Test-statistic palette: null-grey `#525252`; significant-tail `--accent` cobalt; genome-wide threshold amber `#b45309`.
- Escape `&`, `<`, `>` as XML entities (`&amp;`, `&lt;`, `&gt;`).

## Figure Budget

Twelve figures for a ~3h 30min lecture:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | GWAS cohort design | Part 1 | Custom SVG |
| 2 | Single-SNP regression model | Part 1 | Custom SVG |
| 3 | Manhattan plot of a well-powered GWAS | Part 2 | Custom SVG |
| 4 | QQ plot with inflation diagnosis | Part 2 | Custom SVG |
| 5 | Effective number of independent tests | Part 3 | Custom SVG |
| 6 | Population structure in PCA space | Part 4 | Custom SVG |
| 7 | Kinship matrix heatmap | Part 4 | Custom SVG |
| 8 | Mixed-model vs OLS inflation | Part 4 | Custom SVG |
| 9 | PRS portability across ancestries | Part 5 | Custom SVG |
| 10 | Fine-mapping credible set | Part 6 | Custom SVG |
| 11 | Colocalization of GWAS and eQTL | Part 6 | Custom SVG |
| 12 | Burden test vs SKAT at a gene | Part 7 | Custom SVG |

---

## Figure 1 — GWAS cohort design

**File**: `diagrams/lecture-13/01-cohort-design.svg`
**Lecture anchor**: §1.1 Why association testing
**ViewBox**: `0 0 1080 440`

### Purpose

Introduce the two standard GWAS designs (case-control and quantitative-trait) and show how each individual provides one row in both a phenotype vector and a genotype matrix.

### Content

**Left panel — case-control design.** Two stacked cohorts: "Cases" (10 individuals with a red affected-marker icon) and "Controls" (10 individuals, no marker). Each individual has a row of 6 SNP cells coloured by genotype (0/1/2 count). To the right of the two cohorts, a comparison arrow labelled "per-SNP test: Δ allele frequency, χ²".

**Right panel — quantitative-trait design.** A single cohort of 20 individuals, each with a continuous phenotype bar (e.g. height in cm) and the same genotype row. Comparison arrow: "per-SNP test: linear regression β".

**Shared caption bar across bottom.** "Both designs: one regression per SNP; M ≈ 10⁶ tests; genome-wide significance at 5×10⁻⁸."

### Style notes

- Genotype cells: 0 copies white, 1 copy `--accent` cobalt (light), 2 copies cobalt (dark).
- Affected marker: red circle.
- Phenotype bars: horizontal, `--accent-soft` fill.
- Two panels separated by a thin vertical divider.

---

## Figure 2 — Single-SNP regression model

**File**: `diagrams/lecture-13/02-single-snp-regression.svg`
**Lecture anchor**: §1.2 Per-SNP regression
**ViewBox**: `0 0 1080 420`

### Purpose

Show the additive genetic regression model: genotype (0, 1, 2) on the x-axis, phenotype on the y-axis, a straight line through the three mean values, with effect-size $\beta$ annotated as the slope.

### Content

**Main panel (left ~65%).** Scatter plot. X-axis: genotype $G \in \{0, 1, 2\}$. Y-axis: phenotype (e.g. LDL cholesterol, mg/dL). Three clusters of ~30 points each centred near genotype values; vertical jitter within each cluster represents inter-individual noise. Mean horizontal line per cluster. Regression line: $Y = \beta_0 + \beta G$ drawn through the three means; $\beta$ labelled as the slope; $\beta_0$ labelled as the intercept.

**Right inset — hypothesis panel.** Two overlapping Gaussian curves for the sampling distribution of $\hat{\beta}$ under $H_0: \beta = 0$ (grey) and $H_1: \beta > 0$ (cobalt, shifted right). Observed $\hat{\beta}$ shown as a vertical tick. Shaded tail region: "p-value". $z = \hat{\beta} / \text{SE}(\hat{\beta})$ labelled.

### Style notes

- Scatter points: small dots at `--accent-soft`; means as slightly-larger dots in `--fg`.
- Regression line: bold `--accent` cobalt.
- Hypothesis-panel Gaussians: grey for null, cobalt for alternative.
- Tail shading: amber.

---

## Figure 3 — Manhattan plot of a well-powered GWAS

**File**: `diagrams/lecture-13/03-manhattan-plot.svg`
**Lecture anchor**: §2.1 Manhattan plots
**ViewBox**: `0 0 1200 440`

### Purpose

The canonical GWAS display. Reader leaves able to read peak positions, peak heights, and interpret the genome-wide significance threshold.

### Content

**Plot area.** X-axis spans 22 autosomes concatenated left-to-right, labelled "chr1 … chr22". Y-axis: $-\log_{10}(p)$ from 0 to 30. ~1000 dots scattered at varying heights; most between 0 and 5. Alternating chromosome tints (alternating panels of `--bg-subtle` vs `--bg-paper`).

**Three or four tall peaks** clearly visible at chr 2, chr 10, chr 16, chr 18. Each peak has neighbouring dots also elevated (LD propagation — a peak is a cluster, not a spike). Annotations above the peaks with nearby gene names: *FTO*, *TCF7L2*, *LCT*, *HLA*.

**Horizontal dashed line at $y = 7.3$** labelled "Genome-wide significance (5×10⁻⁸)". Minor dotted line at $y = 5$ labelled "Suggestive threshold".

**Axis label at top-left**: "N = 300,000 participants; M = 10⁶ SNPs tested".

### Style notes

- Dots: 1.5 px circles; alternating chromosome tints for readability.
- Peaks: dots colourised darker cobalt when above significance.
- Annotation text in italic.
- Threshold lines: dashed amber.

---

## Figure 4 — QQ plot with inflation diagnosis

**File**: `diagrams/lecture-13/04-qq-plot.svg`
**Lecture anchor**: §2.2 QQ plots and inflation
**ViewBox**: `0 0 1080 440`

### Purpose

Teach the reader how to diagnose inflation vs genuine polygenic signal from a QQ plot.

### Content

**Two side-by-side panels**, each a QQ plot.

**Left — well-controlled.** X-axis: expected $-\log_{10}(p)$ from 0 to 8. Y-axis: observed $-\log_{10}(p)$ from 0 to 30. Diagonal 45° reference line. Dot cloud follows the diagonal closely for most of the range, then departs up-and-right near the top — a few genome-wide hits. Annotation: "λ_GC = 1.02". Title: "Well-controlled GWAS: late departure = real signal".

**Right — inflated.** Same axes. Dot cloud lifts above the diagonal starting near the origin and stays above it throughout. Annotation: "λ_GC = 1.25". Title: "Inflated GWAS: early departure = confounding".

**Bottom caption bar.** "λ_GC > 1.05 warrants investigation. Distinguish polygenic signal (LDSC intercept near 1) from true confounding (LDSC intercept > 1)."

### Style notes

- Reference line: dashed grey.
- Dot cloud: small cobalt dots.
- Lift-off region (top-right of left panel): highlighted with a small curved arrow.
- Early lift (right panel): shaded region beneath the cloud.

---

## Figure 5 — Effective number of independent tests

**File**: `diagrams/lecture-13/05-effective-tests.svg`
**Lecture anchor**: §3.1 Why 5×10⁻⁸
**ViewBox**: `0 0 1080 440`

### Purpose

Explain how the Bonferroni threshold of $5 \times 10^{-8}$ arises from the effective (LD-independent) count of tests, not the raw SNP count.

### Content

**Top row — raw SNPs.** A horizontal bar along the page (representing 1 Mb of chromosome) with ~300 short tick marks crowded together. Label: "~10 M SNPs genome-wide (imputed)". Small count badge: "M_raw ≈ 10⁷".

**Middle row — LD blocks.** Same bar, now tick marks grouped into ~30 wider blocks (representing LD blocks). Each block drawn as a light-grey filled rectangle covering its SNPs. Label: "~1 M LD-independent blocks (Pe'er et al. 2008)". Count badge: "M_eff ≈ 10⁶".

**Bottom — Bonferroni math.** Large equation: "α_adjusted = 0.05 / M_eff = 0.05 / 10⁶ = 5 × 10⁻⁸". Note: "⟵ the genome-wide significance threshold".

### Style notes

- Tick marks: thin vertical lines.
- LD blocks: filled rectangles in `--bg-subtle`.
- Equation: large JetBrains Mono text, centred.
- Arrow from "M_eff" badge to the denominator in the equation.

---

## Figure 6 — Population structure in PCA space

**File**: `diagrams/lecture-13/06-pca-structure.svg`
**Lecture anchor**: §4.1 Population stratification
**ViewBox**: `0 0 1080 440`

### Purpose

Show how the top two PCs of a genotype matrix reveal ancestry groups and how phenotype can track ancestry, producing spurious associations if uncorrected.

### Content

**Two side-by-side scatter panels**, both PC1 vs PC2.

**Left — coloured by ancestry.** Five visible clusters: cobalt (European, large), red (East Asian), amber (African), green (South Asian), violet (admixed American; stretched between European and African clusters). Axis labels "PC1 (explains 4.2%)", "PC2 (explains 1.8%)".

**Right — coloured by phenotype.** Same point positions. Colour now encodes phenotype value (e.g. fasting glucose) on a grey-to-amber continuous scale. A gradient is visible: lower values on the European side, higher values on the African side. Annotation: "phenotype tracks ancestry → uncorrected GWAS will call ancestry-differentiated SNPs significant".

### Style notes

- Scatter dots: 3 px filled circles.
- Ancestry palette as in Lecture 12 (consistency).
- Right-panel colour scale: sequential amber gradient.
- Small legend strip at the bottom of each panel.

---

## Figure 7 — Kinship matrix heatmap

**File**: `diagrams/lecture-13/07-kinship-matrix.svg`
**Lecture anchor**: §4.3 Cryptic relatedness and kinship
**ViewBox**: `0 0 1000 520`

### Purpose

Visualise the pairwise-relatedness structure of a cohort of ~100 individuals.

### Content

**Main heatmap.** A square 100 × 100 matrix of coloured cells. Colour scale: white (kinship ≈ 0, unrelated) to deep cobalt (kinship = 0.5, first-degree). Diagonal: all cells dark cobalt (self-kinship).

**Two small blocks off the diagonal** representing related pairs: a 3×3 block near (20, 20) (a sibling trio) and a 2×2 block near (70, 70) (a parent–child pair).

**Colour-scale legend on the right**: 0.0 white, 0.05 pale cobalt, 0.25 cobalt, 0.5 deep cobalt. Labels: "unrelated", "2nd cousin", "sib/parent-child", "MZ twin/self".

**Text annotation below the plot**: "BOLT-LMM / SAIGE / REGENIE all consume this matrix as a random-effect covariance component."

### Style notes

- Axes: individual index 1..100 along both x and y.
- Colour scale: sequential cobalt.
- Annotated block outlines: thin amber rectangles.

---

## Figure 8 — Mixed-model vs OLS inflation

**File**: `diagrams/lecture-13/08-lmm-vs-ols.svg`
**Lecture anchor**: §4.4 Mixed-model regression
**ViewBox**: `0 0 1080 440`

### Purpose

Show, on the **same** simulated cohort, that OLS regression inflates and mixed-model regression (LMM) restores calibration.

### Content

**Two side-by-side QQ panels** (same axes as Figure 4).

**Left — OLS.** Dot cloud lifts above the diagonal starting at low-$-\log_{10}(p)$. Annotation: "OLS: λ_GC = 1.15. No kinship adjustment → inflated.".

**Right — LMM (BOLT-LMM / REGENIE).** Same dots, now close to the diagonal. Only the top-right tail lifts off. Annotation: "LMM: λ_GC = 1.02. Kinship-structured covariance absorbs relatedness and residual structure."

**Shared caption at the bottom**: "Same cohort, two methods. OLS: GLS with V=I. LMM: GLS with V = σ_g² K + σ_e² I."

### Style notes

- Matched palette to Figure 4.
- Small arrow between panels labelled "noise whitening via V⁻¹".

---

## Figure 9 — PRS portability across ancestries

**File**: `diagrams/lecture-13/09-prs-portability.svg`
**Lecture anchor**: §5.4 Polygenic Risk Scores
**ViewBox**: `0 0 1080 440`

### Purpose

Show the portability problem: a PRS trained on European GWAS data loses 3–5× predictive power when applied to non-European ancestries.

### Content

**Left — violin plots of PRS distribution per ancestry.** Three stacked violin plots for EUR / EAS / AFR ancestries. All centered near 0 with similar spread (the PRS is zero-meaned per-ancestry by convention). Phenotype correlation annotated above each violin.

**Right — bar chart of $R^2$ (PRS vs phenotype) per ancestry.** Three bars: EUR (training cohort) $R^2 = 0.12$ (tallest, cobalt), EAS $R^2 = 0.05$ (amber), AFR $R^2 = 0.03$ (red). Baseline dashed line at $R^2 = 0$.

**Caption bar at the bottom.** "PRS trained on UK-Biobank-style European GWAS. Portability drops with genetic distance from training population. Remedies: diverse GWAS + ancestry-aware PRS methods."

### Style notes

- Violins: filled with each ancestry's colour at low opacity.
- Bars: solid fills.
- Annotation callouts use italic.

---

## Figure 10 — Fine-mapping credible set

**File**: `diagrams/lecture-13/10-fine-mapping.svg`
**Lecture anchor**: §6.2 SuSiE, FINEMAP, CAVIAR
**ViewBox**: `0 0 1200 520`

### Purpose

Demonstrate how fine-mapping narrows a GWAS peak from hundreds of correlated SNPs down to a small credible set.

### Content

**Three stacked tracks** sharing a genomic coordinate axis (a single 300 kb locus).

**Top track — GWAS $-\log_{10}(p)$.** ~200 SNP dots; a clear peak in the middle rising to $-\log_{10}(p) = 25$. Smooth envelope over the peak region.

**Middle track — LD (r²) to lead SNP.** Same x-axis. Each SNP coloured by its $r^2$ to the lead SNP: deep red for $r^2 > 0.8$, orange for $0.5$–$0.8$, yellow for $0.2$–$0.5$, grey otherwise. A classic LD-block appearance.

**Bottom track — SuSiE Posterior Inclusion Probability (PIP).** Same x-axis. PIP bars, most near zero, five SNPs with PIP between 0.1 and 0.5 highlighted as the credible set. Credible-set annotation: "95% credible set: 5 SNPs, total PIP = 0.97". One SNP in the credible set has PIP 0.48 — the leading candidate.

**Shared caption.** "Fine-mapping: SuSiE narrows 200 peak SNPs to a 5-SNP credible set containing the causal variant with 95% posterior probability."

### Style notes

- Three tracks aligned vertically, sharing an x-axis label ("Genomic position (kb)").
- LD track: standard locus-zoom colour scheme.
- PIP bars: vertical bars in cobalt; credible-set members in amber.

---

## Figure 11 — Colocalization of GWAS and eQTL

**File**: `diagrams/lecture-13/11-colocalization.svg`
**Lecture anchor**: §6.3 Colocalization
**ViewBox**: `0 0 1200 440`

### Purpose

Show the co-occurring peaks of a GWAS trait and a tissue-specific eQTL at the same locus — the basis of coloc's $P(H_4)$ test.

### Content

**Two stacked locus-zoom panels** sharing a genomic coordinate axis (~200 kb).

**Top — GWAS LDL cholesterol signal.** ~200 SNP dots; clear peak at position ~100 kb rising to $-\log_{10}(p) = 18$. Lead SNP annotated.

**Bottom — eQTL for *SORT1* in liver tissue.** Same x-axis. ~200 SNP dots; peak at the **same position** rising to $-\log_{10}(p) = 15$. Lead SNP annotated (same or very close to the GWAS lead).

**Annotation box on the right.** "coloc output: P(H_0) = 0.00, P(H_1) = 0.02, P(H_2) = 0.00, P(H_3) = 0.05, P(H_4) = 0.93. → Colocalizes. SORT1 is a candidate causal gene for LDL regulation."

### Style notes

- Two locus-zoom plots sharing the x-axis.
- GWAS panel dots cobalt; eQTL panel dots amber.
- Vertical dashed line through both panels at the shared peak position.
- Annotation box: `--bg-subtle` fill, `--accent` border.

---

## Figure 12 — Burden test vs SKAT at a gene

**File**: `diagrams/lecture-13/12-burden-skat.svg`
**Lecture anchor**: §7.3 SKAT and SKAT-O
**ViewBox**: `0 0 1200 480`

### Purpose

Illustrate when burden tests vs SKAT are each optimal, using three cartoon genes with different rare-variant directional patterns.

### Content

**Three side-by-side gene cartoons**, each a long horizontal rectangle representing a gene with ~20 rare-variant tick marks along its body.

**Left gene — "all deleterious".** All 20 tick marks red (deleterious direction). Effect directions point the same way.

**Middle gene — "mixed directions".** 8 red, 7 green (protective), 5 grey (neutral). Effects point in different directions.

**Right gene — "few large effects".** 18 grey, 2 very red (large deleterious effect).

**Below each gene**, two small test-statistic bars:

- Burden statistic (height).
- SKAT statistic (height).

**Summary table at the bottom.** Columns: Gene / Directional pattern / Burden test significant / SKAT significant / SKAT-O significant. Entries illustrate:

- Left ("all deleterious"): Burden **✓**, SKAT **✓**, SKAT-O **✓**.
- Middle ("mixed directions"): Burden **✗**, SKAT **✓**, SKAT-O **✓**.
- Right ("few large effects"): Burden weak, SKAT **✓**, SKAT-O **✓**.

### Style notes

- Gene body: light grey rectangle.
- Tick marks: thin vertical lines coloured by variant direction.
- Bars: small, labelled.
- Summary table: minimal Inter-rendered with cells separated by thin lines.
