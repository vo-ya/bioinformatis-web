# Lecture 6 — Differential Expression and Count Statistics

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-06.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Explain why a t-test on RNA-seq counts produces unreliable p-values, and state the two properties of the count distribution that break it.
2. Write down the negative-binomial generalised linear model (GLM) that DESeq2 and edgeR fit, and identify the role of the size-factor offset.
3. Describe the dispersion-estimation pipeline: gene-wise maximum likelihood, fitted mean–dispersion trend, empirical-Bayes shrinkage.
4. Pick the right test — Wald vs likelihood-ratio test (LRT) — for a given DE question, and explain the difference.
5. Explain why the Benjamini–Hochberg procedure controls the false-discovery rate at 20k-gene scale, and contrast it with Bonferroni.
6. Interpret a p-value histogram: what a healthy RNA-seq DE analysis should look like, and what the failure modes are.
7. Pick between over-representation analysis (ORA) and gene-set enrichment analysis (GSEA) given a DE result, and describe what each actually tests.
8. Sketch a DE workflow end-to-end, from the count matrix produced by Lecture 5 to an annotated gene list suitable for biological interpretation.

---

## Part 1 — The DE Problem (≈25 min)

### 1.1 What we're asking (≈8 min)

"Differential expression" is the question a biologist usually wants to answer: given RNA-seq counts from samples under two or more conditions, which genes are expressed at measurably different levels between the conditions? The output is a list of genes, each with a log-fold-change estimate and a confidence statement.

Mechanically, differential expression is a gene-by-gene hypothesis test with thousands of tests run in parallel. For each gene:

- H₀: expression is the same in both conditions.
- H₁: expression differs by some non-zero amount.

The output you give a biologist is a **volcano plot**: log₂ fold change on the x-axis, −log₁₀ adjusted p-value on the y-axis. Genes far from the origin are "differentially expressed" (large effect, high confidence). Genes near the origin are not.

**FIGURE — Figure #1: Anatomy of a volcano plot** → `diagrams/lecture-06/01-volcano-anatomy.svg`
*A volcano plot with ~2,000 gene points, the significance / fold-change quadrants labeled, and the canonical thresholds (padj &lt; 0.05, |log₂FC| &gt; 1) drawn as dashed lines.*

**EMBED — Artifact #1: Volcano Plot Explorer** → `artifacts/lecture-06/01-volcano-explorer.html`
*Load a preset DE result (or paste your own); drag the significance / fold-change thresholds and watch the number of "hits" update live.*

> **Intuition box**: A volcano plot answers two questions for every gene at once: how big is the effect (horizontal position) and how certain are we (vertical position)? The whole of this lecture is about how the vertical axis — the adjusted p-value — gets computed honestly. Getting the horizontal axis wrong is a small problem; getting the vertical axis wrong turns every DE analysis into a random-number generator.

### 1.2 Why t-tests fail on counts (≈10 min)

A reasonable first instinct: "I have counts in two conditions; run a t-test per gene." This fails for three reasons rooted in the data.

**Reason 1: counts are integers, not normals.** A Student's t-test assumes data are drawn from a normal distribution with a common variance parameter. Counts violate both — they're non-negative integers, they have a hard floor at zero, and their variance is a *function* of their mean (bigger counts mean bigger variance, §1.3). Low-count genes are the ones most affected, and low-count genes are the ones most often relevant to biology.

**Reason 2: too few replicates.** A typical bulk RNA-seq study has 3–6 replicates per condition. That's not enough to estimate the per-gene variance reliably with a t-test — the sample standard deviation from 3 numbers has a roughly 40% relative error. With 20,000 genes tested, thousands of them will have a *spuriously small* sample variance and will look "highly significant" purely by chance. A test that is correct only on average becomes spectacularly wrong when you run it 20,000 times and pick the "winners."

**Reason 3: the mean–variance relationship.** In real RNA-seq data, the variance of a gene's counts grows roughly quadratically with its mean (the NB model from Lecture 5 §4.4). A t-test assumes homoscedastic variance. On RNA-seq data, the t-test will underestimate variance for highly expressed genes and overestimate it for lowly expressed genes — producing p-values that are systematically wrong in opposite directions depending on gene expression level.

**FIGURE — Figure #2: Why t-tests fail on counts** → `diagrams/lecture-06/02-t-test-failure.svg`
*Left: a typical low-count gene histogram vs a normal of the same mean — visible skew and floor at zero. Right: the resulting p-value distribution under the null for a t-test on count data — non-uniform, biased toward false positives.*

> **EE framing**: RNA-seq counts are photon-counting-style observations with additional multiplicative noise. The right detection statistic is not the same as for Gaussian-noise channels — you need a likelihood-ratio test built against the correct noise distribution, not a fixed-variance t-test. Using a t-test on counts is like using a matched-filter detector assuming white noise when the actual noise is 1/f shot noise: the detector runs, it produces numbers, and the numbers are wrong in systematic ways.

### 1.3 What the noise model must capture (≈7 min)

The DE noise model has to represent three empirical facts:

1. **Discrete integer counts with a zero floor.** A gene with true mean 0.3 reads per sample will have many samples with exactly 0 observed reads. A Gaussian model can't represent that.
2. **Variance grows with the mean.** Low-count genes have variance near their mean (Poisson-like at very low counts); high-count genes have variance much larger than their mean (overdispersion).
3. **The overdispersion is structured.** Genes of similar expression level tend to have similar overdispersion. This is the hook that lets empirical Bayes work (§3.3).

The **negative binomial** distribution captures all three at the cost of one extra parameter beyond the Poisson: a per-gene **dispersion** α that sets the quadratic-variance term. That's the model DESeq2 and edgeR fit. The rest of the lecture is how to fit it well with only a few replicates.

> **Warning box**: "Just use more replicates" is sometimes the right answer, but it rarely rescues a study whose replicate count was set by budget. At n = 2–3 per condition you have almost no power to estimate dispersion per gene; you're going to have to share information across genes (§3.3). Any method that doesn't share information will be over-confident on some genes and under-confident on others, with no principled way to predict which.

---

## Part 2 — Negative Binomial and the GLM (≈45 min)

### 2.1 The NB model in more depth (≈12 min)

Lecture 5 §4.4 introduced the NB as an overdispersed Poisson — a Poisson whose rate parameter is itself drawn from a gamma distribution. Writing that out:

$$ Y_{ig} \sim \text{NB}(\mu_{ig}, \alpha_g) $$

with mean E[Y] = μ and variance Var[Y] = μ + α·μ². The dispersion α is per-gene (not per-sample, not per-condition) because it's a property of how noisily that particular gene is expressed, and it's the thing we'll spend most of §3 estimating.

Three important special cases:

- **α = 0** gives the Poisson. Observed at very low counts where there isn't much to be overdispersed about.
- **Small α** (~0.01) gives near-Poisson behaviour. Typical for highly-expressed housekeeping genes with tight biological control.
- **Large α** (~0.3) gives heavy overdispersion. Typical for genes with bursty expression or strong biological variability.

A typical RNA-seq dataset has dispersions spanning two orders of magnitude across the 20,000 genes. The mean-dispersion relationship follows a predictable pattern: low-mean genes have high dispersion (the α·μ² term is small but Poisson-like floor dominates the noise), medium-mean genes have moderate dispersion, high-mean genes have low dispersion (tight biological control on the highly-expressed genes).

**FIGURE — Figure #3: Mean–dispersion relationship in RNA-seq** → `diagrams/lecture-06/03-mean-dispersion.svg`
*Log-log scatter of per-gene mean (x) vs per-gene dispersion (y) from a real RNA-seq dataset (~2,000 points); the characteristic downward-and-flat trend visible; a fitted trend line overlaid.*

> **EE framing**: The dispersion parameter is the noise-floor characterisation of each gene's expression channel. High-dispersion genes are noisy signal sources where even a big effect is hard to distinguish from fluctuation; low-dispersion genes are clean channels where small effects are visible. The DE problem is detection on a channel whose noise level you have to estimate from the same few samples you're trying to detect in — with 20,000 channels to estimate jointly.

### 2.2 The NB generalised linear model (≈15 min)

DESeq2 and edgeR fit a **generalised linear model** (GLM) per gene. The canonical form, for gene g and sample i:

$$ \log \mu_{ig} = \log s_i + \beta_{0,g} + \beta_{1,g} \cdot x_i $$

where:

- μ_{ig} is the expected count for gene g in sample i.
- s_i is the size factor for sample i (from DESeq2's median-of-ratios, Lecture 5 §4.3).
- x_i is the condition indicator (0 for control, 1 for treatment).
- β₀,g is the gene's baseline log-expression.
- β₁,g is the **log-fold-change** for gene g between conditions — the thing we want to test against zero.

More generally, x_i can be a full design matrix with multiple coefficients (batch, sex, condition, treatment × time interactions, etc.). The GLM machinery is exactly the same as in linear regression, but with the log-link function and the NB likelihood.

**FIGURE — Figure #4: NB GLM — design matrix and link function** → `diagrams/lecture-06/04-nb-glm-design.svg`
*A 6-sample toy dataset: design matrix (intercept + condition + batch), gene counts, log-link equation, fitted coefficients. Arrows connect design-matrix columns to their role in the prediction.*

> **EE framing**: The NB GLM is maximum-likelihood inference with a log-link — structurally identical to logistic regression (Bernoulli likelihood + logit link) or Poisson regression (Poisson likelihood + log link). The dispersion α replaces the scale parameter you'd have in OLS. Every DE tool that ends in "-seq2" or "-R" ultimately fits this model; they differ in how they estimate α and which likelihood-ratio or Wald form they use for testing.

The practical consequence: once the GLM is fit, every DE question becomes "is some linear combination of β coefficients significantly different from zero?" — which is just a standard hypothesis test.

### 2.3 Size factors and log-fold change (≈10 min)

The size factor s_i in the GLM is an **offset**, not a coefficient. We treat it as known (it was estimated from the full count matrix in Lecture 5 §4.3). In log space:

$$ \log \mu_{ig} = \log s_i + \beta_{0,g} + \beta_{1,g} \cdot x_i $$

the log s_i term shifts the prediction up for samples with deeper libraries and down for samples with shallower ones, leaving the β coefficients to carry the *composition-normalised* expression change. Without the offset, β₁ would conflate "more reads in this sample" with "more expression of this gene" — exactly the composition confounding that §4.3 of Lecture 5 was designed to remove.

The coefficient β₁,g, after the offset has absorbed depth, is the **log-fold change** between conditions:

$$ \beta_{1,g} = \log\left(\frac{\mu_{g,\text{treated}}}{\mu_{g,\text{control}}}\right) $$

converting to log₂ base (the biologist's convention) by dividing by log(2): `log2FC = β₁ / log(2)`.

> **Intuition box**: A size factor is not a per-gene multiplier. It's a per-sample multiplier applied to *all* genes in that sample. Think of it as each sample's "sequencing depth in composition-adjusted units." A sample with size factor 1.2 sequenced 20% more reads than the population average; every gene's prediction in that sample is scaled up by 1.2× before comparing to observations. Size factors are the bridge between a count and a comparable-across-samples abundance.

### 2.4 Biological versus technical variance (≈8 min)

Overdispersion in RNA-seq is not a single thing. It has two sources, and distinguishing them matters for experimental design.

**Technical variance** comes from the library-prep and sequencing pipeline: pipetting error, polyadenylation bias, fragmentation noise, PCR duplicates that weren't marked, sequencer calibration drift. This is the component that shrinks with more reads per sample (higher depth).

**Biological variance** comes from real expression differences between replicate samples — in animal experiments, slight differences in genetic background, circadian phase, cage environment. Even in "identical" cell-line samples, there's cell-cycle phase heterogeneity and clone-to-clone drift. This is the component that shrinks only with more biological replicates, never with more reads.

In practice, biological variance dominates for medium-to-highly-expressed genes, and technical variance dominates for low-count genes. The consequence for experiment design: if your effect is large and your genes are highly expressed, 3 replicates suffice; if your effect is small or you care about low-count genes (transcription factors, regulators), you need 6+ replicates. More reads per sample doesn't fix the problem.

> **Warning box**: The replicate question is unforgiving. Published DE analyses with 2 replicates per condition are common and mostly unreliable. DESeq2 and edgeR still run on n=2, but the dispersion estimates are noise, the empirical-Bayes shrinkage does its best, and the p-values come out with wide and often miscalibrated confidence intervals. If you see "n=2 per group" in a methods section, treat the gene list as exploratory and prioritise follow-up validation.

---

## Part 3 — DESeq2 and edgeR Mechanics (≈60 min)

### 3.1 The DE pipeline, end to end (≈5 min)

Both DESeq2 and edgeR follow the same conceptual flow, differing mostly in the details of dispersion estimation:

```
count matrix (from L5)
    ↓
estimate size factors                 (callback: L5 §4.3, median of ratios)
    ↓
estimate gene-wise dispersion         (§3.2)
    ↓
fit a mean–dispersion trend           (§3.3)
    ↓
shrink gene-wise dispersions toward the trend (empirical Bayes)
    ↓
fit NB GLM per gene using shrunk dispersion   (§2.2)
    ↓
test β coefficients (Wald or LRT)     (§3.4)
    ↓
raw p-values (one per gene)
    ↓
multiple-testing adjustment           (BH, §4)
    ↓
filtered / annotated gene list        (§5)
```

Every substantive disagreement between DESeq2 and edgeR is about how steps 3–5 are done. The rest is bookkeeping.

### 3.2 Per-gene dispersion estimation (≈15 min)

Given the counts for gene g across n samples and the fitted intercept β₀,g, the per-gene maximum-likelihood estimate of α_g is the value that makes the observed counts most probable under the NB model. Computed by a couple of iterations of Newton-Raphson on the NB log-likelihood; in code, DESeq2 uses `nbinomDispEstimate`, edgeR uses `estimateGLMCommonDisp`.

At n = 3–6 samples per condition, this estimator is essentially garbage. A single-gene dispersion MLE with 3–6 observations has a confidence interval that spans more than two orders of magnitude. That's why the next step is critical.

**FIGURE — Figure #5: Per-gene dispersion — MLE on few samples** → `diagrams/lecture-06/05-dispersion-mle.svg`
*Left: three replicates per condition, gene counts plotted. Center: the NB log-likelihood surface for α — broad, flat, highly uncertain. Right: the MLE from 3 replicates (one point with a massive CI bar) vs the MLE from 20 replicates (narrow CI) for comparison.*

Two details that matter:

- **Pooling across conditions matters**: for dispersion, we pool residuals from *all* samples of a gene (not just one condition), because we assume the dispersion is the same in both conditions. Only the mean differs.
- **Outliers bite hard**: a single aberrant count in a gene with 3 replicates per condition will inflate the gene-wise MLE by an order of magnitude. DESeq2 replaces such outliers with an imputed count using Cook's distance.

> **Warning box**: The gene-wise dispersion estimate from 3 replicates is not a reliable number on its own. If you see a per-gene dispersion value reported in a DE pipeline output, it's probably the *shrunk* estimate (post §3.3), not the raw MLE. If it's the raw MLE, treat it with deep suspicion.

### 3.3 Empirical-Bayes shrinkage toward the fitted trend (≈15 min)

Here is the central statistical idea of DESeq2 and edgeR, and the reason they work on 3-replicate studies.

**The empirical observation**: if you plot per-gene mean against per-gene (raw, MLE) dispersion for all 20,000 genes in a dataset, the points form a clear downward-sloping cloud. Low-expression genes have high dispersion on average; high-expression genes have low dispersion on average. This trend is real and it's consistent across datasets.

**The algorithmic move**: fit a parametric trend through the cloud (DESeq2 uses a regression of the form `dispersion = α₀ + α₁/mean`; edgeR uses a cubic-smoothing-spline variant). Then, for each gene, "shrink" its raw MLE toward the trend line. Genes that disagree strongly with the trend (real outliers) are shrunk less; genes close to the trend are shrunk more. The result is a per-gene dispersion estimate that is informed by both its own data and the global trend — an **empirical-Bayes** estimator.

**FIGURE — Figure #6: Empirical-Bayes dispersion shrinkage** → `diagrams/lecture-06/06-eb-shrinkage.svg`
*The same mean-vs-dispersion scatter as Figure 3, but with raw MLEs as outlined circles and shrunk estimates as filled dots — each outlined circle connected to its shrunk counterpart by a grey arrow. The fitted trend line passes through the cloud. Genes far from the trend get only a little shrinkage; close-to-trend genes move to the trend.*

**EMBED — Artifact #2: Dispersion Shrinkage Visualizer** → `artifacts/lecture-06/02-dispersion-shrinkage.html`
*A mean-vs-dispersion scatter with draggable genes. Toggle shrinkage on/off; see the raw MLE vs the fitted-and-shrunk estimate; change the strength parameter and watch how much each gene moves.*

The statistical payoff: the shrunk estimator has much smaller variance than the raw MLE for genes near the trend, at the cost of a small bias toward the trend for genes that are genuinely atypical. In large-scale testing that tradeoff is strongly favourable — the variance reduction buys calibrated p-values, the small bias doesn't qualitatively change which genes are significant.

> **EE framing**: Empirical-Bayes shrinkage is regularisation toward a learned prior. The "prior" is the fitted mean–dispersion trend, learned from the data. Shrinking each gene's MLE toward the prior is the frequentist analog of ridge regression: trade a tiny amount of bias for a big reduction in variance. The same structural move shows up in James-Stein estimation (shrink k≥3 means toward their grand mean), BLUP estimation in genetics, and every hierarchical Bayesian model ever written. What makes the RNA-seq application distinctive is that the prior is itself empirical — fitted from the data being tested, not specified up front.

> **Historical pointer**: The DESeq2 paper (Love, Huber, Anders 2014) packaged this approach into a single R function call. It wasn't the first EB-dispersion method — edgeR had a version in 2010, DSS in 2013 — but it was the first that became the community default. Its p-value calibration on realistic data was the thing that won out. As of 2024, the vast majority of published bulk RNA-seq DE analyses are DESeq2.

### 3.4 Hypothesis tests — Wald and LRT (≈15 min)

With the GLM fit (β̂₁,g) and its standard error (SE(β̂₁,g)) in hand, the DE test for gene g against the null H₀: β₁,g = 0 has two canonical forms.

**Wald test.** Compute the z-statistic and its p-value:

$$ z_g = \frac{\hat\beta_{1,g}}{\text{SE}(\hat\beta_{1,g})}, \qquad p_g = 2 \cdot \Phi(-|z_g|) $$

Simple, fast, asymptotically valid. Works well for single-coefficient tests (two-condition comparisons). The default in DESeq2.

**Likelihood-ratio test (LRT).** Fit the full model (with β₁,g) and the reduced model (without β₁,g). Compute:

$$ \lambda_g = 2 \cdot (\ell_g^{\text{full}} - \ell_g^{\text{reduced}}) \sim \chi^2_1 $$

where ℓ is the log-likelihood. The LRT is more robust for **complex** hypotheses — testing multiple coefficients jointly (a condition × time interaction, a factor with three levels), and for small-sample cases where the Wald approximation can be off. Used in DESeq2 with `test="LRT"` and in edgeR's `glmQLFTest`.

**FIGURE — Figure #7: Wald vs LRT — when each is right** → `diagrams/lecture-06/07-wald-vs-lrt.svg`
*Top: Wald test geometry — β̂ / SE as a z-score on a normal density. Bottom: LRT geometry — full vs reduced model's log-likelihoods, their difference as a chi-square statistic. A side-by-side comparison table: simple 2-condition / Wald, multi-level factor / LRT, time-course interaction / LRT.*

**EMBED — Artifact #3: NB GLM Fitter** → `artifacts/lecture-06/03-nb-glm-fitter.html`
*A toy 6-sample 2-condition dataset with one gene. Fit the NB GLM; see β₀, β₁, their standard errors; run both Wald and LRT and compare the resulting p-values. Drag counts to change the fit.*

> **Intuition box**: The Wald test asks "is this coefficient far from zero relative to its own uncertainty?" The LRT asks "does removing this coefficient materially hurt the fit?" They agree for simple cases and big datasets. They diverge when the model is curved near the maximum (Wald's linear approximation breaks), or when multiple coefficients should be tested jointly (Wald would need a quadratic form; LRT handles it directly).

> **Discussion prompt**: You are analysing a time-course RNA-seq experiment with 4 time points and want to test "any difference across time" per gene. Would you use Wald or LRT, and why? (LRT, applied to the time-factor as a whole. Wald would test each pairwise comparison separately and you'd have to correct for those; LRT tests the joint null "all β_t = 0" in one shot with one p-value per gene.)

### 3.5 limma-voom: the variance-stabilising alternative (≈10 min)

**limma-voom** (Law, Chen, Shi, Smyth 2014) is a different approach that ends at the same place: a per-gene p-value and log-fold change.

The core idea: transform the count data so it *looks* Gaussian-enough for standard limma linear models to work. The transformation has two steps:

1. **Log-CPM** the counts (Lecture 5 §4.2). Now the values are on an approximately additive scale.
2. Voom weighting — compute a precision weight per (sample, gene) cell, estimated from the empirical mean–variance relationship. Small counts get small weights; high counts get large weights. Feed weighted log-CPMs into limma's standard linear-model + empirical-Bayes testing framework.

The output: linear-model β coefficients, moderated t-statistics, and p-values — all the machinery of classical microarray analysis, applied to transformed RNA-seq counts. Benchmarks consistently show limma-voom is competitive with DESeq2 and edgeR on most datasets; it's especially fast on large cohorts (thousands of samples) where the GLM-based tools slow down.

> **EE framing**: Voom's variance stabilisation is a *preprocessing* trick rather than a new likelihood. It transforms the data onto a scale where a constant-variance Gaussian model is approximately right, then applies constant-variance Gaussian inference. This is the same philosophy as the Anscombe transform for Poisson data or the Box-Cox family of power transforms — when the noise has a known variance function, transform the data to stabilise variance, then use off-the-shelf tools. The tradeoff: you lose some efficiency (the transform isn't perfect; the subsequent test is using an approximate distribution), but you gain huge computational and methodological simplicity.

Pragmatic rule for 2024: for small experiments (n ≤ 20 per condition), DESeq2 is the default. For large cohorts (n ≥ 50 per condition) or genome-scale studies (GWAS-like designs on expression), limma-voom scales better and gives equivalent answers.

---

## Part 4 — Multiple Testing at Gene Scale (≈40 min)

### 4.1 The problem at 20,000 tests (≈15 min)

A standard DESeq2 run tests around 20,000 genes in parallel. Even if every gene were truly null (no differential expression at all), an α = 0.05 per-test threshold would produce 20,000 × 0.05 = 1,000 "significant" genes purely by chance. No biologist will believe that any single raw p-value &lt; 0.05 result means anything in this regime — and rightly so.

A useful way to visualise what's going on: **the p-value histogram.** Under the null hypothesis, p-values are uniformly distributed on [0, 1]. In a real DE analysis most genes are null, but some fraction has real effects. The histogram should look like: a flat uniform region (the nulls), with a spike near zero (the real effects). If the histogram has a spike near 1 (U-shape) or a systematic slope, the p-values are miscalibrated — a sign of model misspecification.

**FIGURE — Figure #8: The p-value histogram diagnostic** → `diagrams/lecture-06/08-pvalue-histogram.svg`
*Three panels: (a) healthy — flat bulk with spike near 0, ~5% of genes with effect. (b) no signal — entirely flat. (c) miscalibrated — U-shape or slope, indicating model problems.*

> **Intuition box**: The p-value histogram tells you everything about a DE analysis in one glance. Healthy histograms mean the testing framework is well-calibrated and there's real signal. U-shapes or slopes mean the p-values themselves are wrong — and no amount of multiple-testing correction will save you. Always plot the histogram before reporting results.

### 4.2 Bonferroni vs Benjamini–Hochberg (≈10 min)

Two canonical ways to control false positives when running many tests.

**Bonferroni** controls the **family-wise error rate (FWER)**: the probability of making *any* false discovery at all. With m tests at per-test level α, Bonferroni rejects only those with `p ≤ α/m`. For 20,000 genes at α = 0.05 Bonferroni, the per-gene threshold is 2.5 × 10⁻⁶. Highly conservative — you'll get a handful of very-strong hits and miss most real effects.

**Benjamini–Hochberg (BH)** controls the **false-discovery rate (FDR)**: the expected *fraction* of false discoveries among those rejected. Sort p-values ascending, p₍₁₎ ≤ p₍₂₎ ≤ … ≤ p₍ₘ₎. Find the largest k such that `p₍ₖ₎ ≤ (k/m)·α`. Reject the k hypotheses corresponding to the smallest k p-values.

BH's interpretation: "of the ~500 genes on my gene list, I'm willing to accept that up to 5% are false positives in exchange for not missing the 95%." That's an acceptable trade for biological follow-up where wet-lab validation will filter false hits.

**FIGURE — Figure #9: The Benjamini–Hochberg procedure** → `diagrams/lecture-06/09-bh-procedure.svg`
*Sorted p-values plotted against rank k; the BH line p = (k/m)·α overlaid; the largest k where p₍ₖ₎ is below the line marked. Left panel: many null genes. Right panel: many discoveries.*

**EMBED — Artifact #4: FDR Simulator** → `artifacts/lecture-06/04-fdr-simulator.html`
*Simulate 20,000 genes with a user-controlled fraction of true positives and effect size. Run both Bonferroni and BH at several FDR thresholds. Watch empirical FDR and power update live as you slide the parameters.*

> **EE framing**: FDR control is detection theory with a controlled false-alarm *rate* rather than a controlled false-alarm *probability*. Bonferroni is Neyman–Pearson with a single-test-level threshold adapted for m independent tests — guarantees that the probability of any false positive across the whole experiment is ≤ α. BH's FDR is closer to the radar-engineering concept of "keep the fraction of false alarms in the detection list below a controlled rate" — more useful in regimes where you're going to follow up on the detections anyway and a few false positives are tolerable. Both have their applications; DE analysis almost always uses FDR.

> **Historical pointer**: The Benjamini–Hochberg procedure was published in 1995 as a theoretical proposal for controlling FDR. It wasn't used much in biology until high-throughput data (microarrays, genotyping arrays) made m-scale testing common around 2000. By 2010 BH was ubiquitous in genomics; by 2024 it's the default in every DE tool without exception.

### 4.3 FDR interpretation and q-values (≈10 min)

A single number from a DE tool's output matters: **padj** (DESeq2) or **FDR** (edgeR) — the BH-adjusted p-value. It has a specific meaning:

> *If you reject all genes with padj ≤ q, then (in expectation) q · (number of rejections) of them are false positives.*

So at padj ≤ 0.05, if the gene list has 500 entries, ~25 are expected false positives.

The **q-value** (Storey 2003) is a refinement: it estimates the proportion of null hypotheses π₀ directly from the p-value distribution (the uniform-like "flat" region), which lets it give a slightly less conservative estimate than BH. For most RNA-seq DE work the BH/padj number and the Storey q-value agree to within a few percent; the distinction matters more in dense-signal regimes than typical bulk RNA-seq.

**Practical rule**: report padj (BH) everywhere. Report q-value in addition if you're doing methodological benchmarking. Don't report raw p-values without an FDR-adjusted companion — raw p-values at large m are meaningless without the adjustment.

### 4.4 Independent filtering and weighted FDR (≈5 min)

One last refinement that DESeq2 applies by default: **independent filtering.**

The idea: genes with extremely low total counts have *zero* chance of being significant even if their effect is large — the variance is so high that no test can discriminate the effect from noise. Including them in the multiple-testing adjustment wastes power. So DESeq2 filters them out *before* running BH, using a data-dependent threshold that maximises the number of rejections.

The filter criterion must be independent of the test statistic under the null, to avoid biasing p-values. Mean count (or more robust variants) satisfies this.

The result: the effective m for BH is smaller — sometimes by 30–50%. The BH threshold is less conservative. More true positives survive.

> **Warning box**: If you're reading DESeq2 output and see genes with NA in the padj column — they were filtered out by independent filtering. This is usually fine; they were almost certainly going to be non-significant anyway. But check: the "filtered" genes column should be mostly low-count, and if you see high-count genes with NA padj, something has gone wrong.

More sophisticated approaches (**Independent Hypothesis Weighting**, IHW; Ignatiadis et al. 2016) generalise this: weight each p-value by an independently-estimated signal covariate (such as gene length or mean expression). Used sparingly in RNA-seq; common in large-cohort studies where the extra power matters.

---

## Part 5 — From Gene List to Biology (≈30 min)

### 5.1 Over-representation analysis (ORA) (≈8 min)

A DE analysis produces a list of "significant" genes. That list is rarely meaningful on its own — a biologist wants to know which *pathways* or *processes* those genes point to. The simplest approach: **over-representation analysis (ORA)**.

Given a curated gene set (say "G2/M cell-cycle genes", 127 genes annotated from MSigDB), and your list of DE genes (say 500 genes), ask: "Is the overlap between these two sets larger than expected by chance?"

Mechanically, this is a 2×2 contingency test:

|  | in gene set | not in gene set | total |
|---|---|---|---|
| DE gene | 42 | 458 | 500 |
| not DE | 85 | 19,415 | 19,500 |
| total | 127 | 19,873 | 20,000 |

Fisher's exact test (or the asymptotic chi-square equivalent) gives a p-value for over-representation. Do this for thousands of gene sets (all of MSigDB, all GO terms), apply BH correction to the resulting p-values, and report the ones with significant enrichment.

**FIGURE — Figure #10: Over-representation analysis — Fisher's exact** → `diagrams/lecture-06/10-ora-fisher.svg`
*A 2×2 contingency table with numbers; the hypergeometric distribution the test is drawn from; the resulting p-value for "G2/M cell-cycle" gene set enrichment in a DE result.*

**EMBED — Artifact #5: ORA Contingency Explorer** → `artifacts/lecture-06/05-ora-contingency.html`
*Adjust DE list size, gene-set size, observed overlap, background gene count; see Fisher's exact p-value and the expected vs observed numbers update live.*

ORA is cheap and intuitive. It has two big drawbacks: (1) it requires a hard threshold for "DE" (loses information from sub-threshold genes), (2) it treats all genes above threshold as equivalent (a gene with log₂FC = 6 counts the same as one with log₂FC = 0.6).

### 5.2 Gene-set enrichment analysis (GSEA) (≈12 min)

**GSEA** (Subramanian et al. 2005) solves both ORA problems by using the full ranked gene list, not a thresholded subset.

Algorithm sketch:
1. Rank all 20,000 genes by some continuous DE statistic (signed log₂FC / -log₁₀p, or a signed Wald z-score).
2. For each gene set S of interest, walk down the ranked list computing a running-sum statistic: +1 weight when you hit a gene in S, −1 weight when you don't (normalised by set size). The running sum goes up when you're enriched for set genes, down when you're not.
3. The **enrichment score (ES)** is the maximum deviation from zero along the walk.
4. Estimate statistical significance by permutation — shuffle gene labels, recompute ES many times, build the null distribution.

**FIGURE — Figure #11: GSEA — the running-sum statistic** → `diagrams/lecture-06/11-gsea-running-sum.svg`
*A ranked gene list at the bottom; above it, the running-sum curve. Peak at the point of maximum enrichment; the ES is the peak height. Inset: permutation null distribution of ES with the observed value marked.*

**EMBED — Artifact #6: GSEA Running-Sum Simulator** → `artifacts/lecture-06/06-gsea-running-sum.html`
*Load a ranked gene list, pick a gene set; step through the walk and watch the running sum build. Run permutations to get an empirical p-value.*

> **EE framing**: GSEA's running-sum statistic is a CUSUM (cumulative sum) detector — the same test used for change-point detection in time-series monitoring. The "change" it detects is a sub-population of labeled elements clustered somewhere in the ranked list. The permutation test makes it nonparametric: no assumption on the null distribution of the ES, just empirical calibration via shuffling. The same algorithmic idea shows up in signal-processing literature as the "rank-sum runs test."

### 5.3 Pathway databases and interpretation caveats (≈10 min)

The quality of any enrichment analysis is bounded by the quality of the gene sets it runs against. Major sources:

- **MSigDB** (Molecular Signatures Database, Broad Institute). The single most widely used collection. Hallmark gene sets (H), curated gene sets (C1–C7), motif-based (C3), cancer-focused (C6), immune-focused (C7). ~25,000 gene sets total.
- **GO (Gene Ontology).** A controlled vocabulary of biological processes, molecular functions, and cellular components, each with its associated gene set. Universal; hierarchical (e.g. "immune response" has many child categories).
- **KEGG** (Kyoto Encyclopedia of Genes and Genomes). Manually curated metabolic and signaling pathways, each with a gene list. Smaller but higher-quality than GO.
- **Reactome.** Human-curated pathways, more detailed biochemistry than KEGG. Popular in European labs.
- **WikiPathways.** Community-curated pathways, edited like Wikipedia.

Pragmatic rules:

- Start with MSigDB Hallmarks (50 high-level gene sets) for a quick overview.
- Use the larger C2 / C7 collections for specific biological questions.
- Always report which database version you used — they update annually, and results shift.
- Be wary of hits with fewer than 5 genes overlapping — even if significant, the biology is thin.

> **Discussion prompt**: A DESeq2 analysis returns ~800 DE genes (FDR &lt; 0.05). GSEA reports strong enrichment (FDR &lt; 0.01) for both "inflammatory response" and "interferon-γ response" — two related but distinct immune gene sets. What's the right way to report this? (Report both. Don't collapse them — they describe different molecular programs, and their distinction may be biologically meaningful. Also look at whether the same genes are driving both hits, or different ones. If a single core set of cytokine genes is responsible for both enrichments, the "story" is about those genes, not about the pathways.)

---

## Wrap-up (≈10 min)

### What you should take away

- **DE is gene-by-gene hypothesis testing at 20,000-test scale with few replicates.** Every step of the pipeline exists because naive approaches (t-test per gene, Bonferroni correction, no dispersion sharing) fail at this scale.
- **The negative-binomial GLM is the backbone.** DESeq2 and edgeR differ in implementation details but fit the same model: log μ = offset + β·x, with NB likelihood and gene-specific dispersion α.
- **Empirical-Bayes shrinkage of dispersion is why DE works on 3 replicates.** Sharing information across genes via the mean-dispersion trend turns a noisy per-gene estimate into a calibrated one. Without this shrinkage, the analysis is uncalibrated.
- **FDR via Benjamini–Hochberg is the right multiple-testing framework.** Bonferroni is too conservative; FDR gives you a list with a controlled fraction of false positives, which is what biological follow-up needs.
- **The p-value histogram is the universal diagnostic.** Flat bulk with near-zero spike = healthy. Any other shape = something is wrong with your model or design.
- **Interpretation requires gene sets, not individual genes.** A list of 500 DE genes tells you almost nothing; an ORA or GSEA run against MSigDB tells you which molecular programs are shifting. Use the enrichment layer as the reporting unit, not the individual-gene layer.

### Next lecture

Single-cell RNA-seq fundamentals: droplet protocols, UMIs, the transition from bulk to per-cell resolution, dimensionality reduction with PCA and UMAP, clustering and cell-type annotation.

### Homework

1. Download a published bulk RNA-seq dataset with at least 3 replicates per condition (GEO accession GSE52778 works — airway smooth muscle treated with dexamethasone, 4 vs 4). Run DESeq2 end-to-end. Report: number of DE genes at padj &lt; 0.05, distribution of log₂FC values, and a screenshot of the p-value histogram.
2. Take the DE result from problem 1. Plot the volcano. What fraction of DE genes has |log₂FC| &gt; 1? Do those genes have smaller or larger raw p-values on average than the rest of the DE set?
3. Implement the BH procedure by hand on a sorted p-value list from your DESeq2 output. Verify your hand-computed padj matches DESeq2's output. At what rank does the BH line touch the sorted p-values for FDR = 0.05?
4. Run GSEA against MSigDB Hallmarks on the DE result. Report the top 5 hallmark gene sets by enrichment score. Compare ORA (Fisher's exact on DE-genes-only) against GSEA (running-sum on the full ranked list). Which method identifies more hits? Are they the same hits?
5. Analyse a simulated dataset where you know the ground truth. Use the FDR Simulator artifact with 20,000 genes, 500 true positives, modest effect size, and 3 replicates. At FDR = 0.05, what fraction of rejected genes are true positives (precision)? What fraction of true positives are rejected (power)? How does this compare to Bonferroni at α = 0.05?

### Recommended reading

- Love, M. I., Huber, W., &amp; Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15, 550. (The DESeq2 paper.)
- Robinson, M. D., McCarthy, D. J., &amp; Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* 26, 139–140.
- Law, C. W., Chen, Y., Shi, W., &amp; Smyth, G. K. (2014). voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology* 15, R29.
- Benjamini, Y., &amp; Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B* 57, 289–300.
- Subramanian, A., Tamayo, P., Mootha, V. K., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102, 15545–15550. (The GSEA paper.)
- Liberzon, A., Birger, C., Thorvaldsdóttir, H., et al. (2015). The Molecular Signatures Database (MSigDB) Hallmark gene set collection. *Cell Systems* 1, 417–425.
- DESeq2 vignette: <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — The DE Problem                            | 25&nbsp;min | 0:25 |
| Part 2 — Negative Binomial and the GLM              | 45&nbsp;min | 1:10 |
| Part 3 — DESeq2 and edgeR Mechanics                 | 60&nbsp;min | 2:10 |
| Part 4 — Multiple Testing at Gene Scale             | 40&nbsp;min | 2:50 |
| Part 5 — From Gene List to Biology                  | 30&nbsp;min | 3:20 |
| Wrap-up                                             | 10&nbsp;min | 3:30 |

**Total:** ~3h 30min of content.
