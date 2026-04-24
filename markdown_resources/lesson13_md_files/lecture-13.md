# Lecture 13 — GWAS and Statistical Genetics

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-13.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Describe the GWAS experimental design (cases vs controls, quantitative trait cohorts) and write down the per-SNP regression model used to test association.
2. Read a Manhattan plot and a QQ plot, interpret genome-wide significance (5×10⁻⁸), and diagnose inflation and LD-driven peak structure.
3. Explain why the GWAS significance threshold is 5×10⁻⁸ in terms of the number of effectively independent tests across the genome.
4. Diagnose population-stratification confounding and correct for it using principal-component covariates and mixed-model regression (BOLT-LMM, SAIGE, REGENIE).
5. Compute a basic polygenic risk score (PRS) from GWAS summary statistics and explain the portability problem across ancestries.
6. Describe the fine-mapping problem as a sparse inverse problem, and apply SuSiE-style credible-set reasoning to a locus.
7. Distinguish single-variant and rare-variant association tests (burden and SKAT), and explain when each is appropriate.

---

## Part 1 — The GWAS Framework (≈25 min)

### 1.1 Why association testing (≈6 min)

By Lecture 12 you know how allele frequencies are distributed in populations and how they evolve. A fundamentally different question motivates this lecture: **which genetic variants affect a phenotype?** Is there a SNP whose allele count correlates with disease risk, height, LDL cholesterol, or schizophrenia susceptibility?

The simple answer (ignoring confounders for one more subsection): for each SNP in the genome, test whether its allele count is associated with the phenotype. That's a regression problem — one per SNP.

Genome-Wide Association Studies (GWAS) scale this to **millions of SNPs across tens of thousands to millions of individuals**. First landmark: the Wellcome Trust Case Control Consortium (WTCCC) 2007 paper testing ~500,000 SNPs in 17,000 individuals for seven common diseases. Most modern studies are 10–1000× larger: UK Biobank (~500k participants), FinnGen (~500k), All of Us (~1M), Million Veteran Program (~1M).

GWAS is now the default first-pass analysis for any genetic contribution to any trait — disease, anthropometric, behavioural, molecular. The mechanics are the same; the scale is enormous.

**FIGURE — Figure #1: GWAS cohort design** → `diagrams/lecture-13/01-cohort-design.svg`
*Two side-by-side panels. Left: case-control design — cases (disease) and controls (no disease) each with genotype arrays (SNP1..SNPk) and the comparison statistic. Right: quantitative-trait design — a cohort of individuals each with a continuous phenotype (e.g. height in cm) and genotypes. Each SNP is tested independently in either design.*

### 1.2 Per-SNP regression (≈10 min)

The model at a single SNP is a simple regression. For a quantitative trait $Y$ (e.g. height), individual $i$, genotype $G_i \in \{0, 1, 2\}$ (count of the "effect" allele under the additive model):

$$Y_i = \beta_0 + \beta \cdot G_i + \boldsymbol{\gamma}^T \mathbf{X}_i + \varepsilon_i$$

where $\mathbf{X}_i$ is a vector of covariates (age, sex, PCs — more in §4), $\boldsymbol{\gamma}$ is the covariate effect vector, and $\varepsilon_i$ is noise. The test is: $H_0 : \beta = 0$ vs $H_1 : \beta \ne 0$.

For a binary trait (case / control), swap linear for logistic:

$$\log \frac{P(Y_i = 1)}{P(Y_i = 0)} = \beta_0 + \beta \cdot G_i + \boldsymbol{\gamma}^T \mathbf{X}_i$$

Three genetic models are standard:

- **Additive** ($G \in \{0, 1, 2\}$ as a numeric predictor). Default in most GWAS. Assumes each extra allele copy has the same effect.
- **Dominant** ($G \in \{0, 1\}$: has at least one effect allele vs none). For traits where a single copy is sufficient.
- **Recessive** ($G \in \{0, 1\}$: has two copies vs not). For traits where both copies matter.

The additive model is the workhorse. Most published GWAS assume it by default.

**FIGURE — Figure #2: Single-SNP regression model** → `diagrams/lecture-13/02-single-snp-regression.svg`
*Scatter plot: genotype (0, 1, 2 along x-axis) vs phenotype (y-axis, e.g. LDL). Three clusters of points at x=0, x=1, x=2 with regression line through the means. Side panel: null and alternative hypothesis distributions with the test statistic marked; effect size β annotated as the slope.*

**EMBED — Artifact #1: GWAS Single-SNP Association Explorer** → `artifacts/lecture-13/01-single-snp-test.html`
*Slide allele frequency, effect size β, sample size N, noise level. See genotype scatter, null and alternative sampling distributions for β̂, p-value, and power. Target aha: GWAS effect sizes are tiny (β typically 0.01–0.05 SD for common variants) — you need huge N to resolve them.*

### 1.3 Scale: millions of SNPs, millions of people (≈5 min)

GWAS data objects:

- **Genotype matrix**: $N$ individuals × $M$ SNPs. Modern arrays genotype ~1M SNPs directly; imputation (using a reference panel like TOPMed, HRC, or the HPRC pangenome) expands this to ~10M SNPs with imputed dosages.
- **Phenotype vector**: $N$-long vector of trait values.
- **Covariate matrix**: age, sex, principal components, cohort indicators — typically 10–30 columns.

The per-SNP regression runs $M$ times, producing $M$ $\beta$-hats, $M$ p-values, $M$ effect-size standard errors. The output is a **GWAS summary statistic table** — arguably the dominant data product of the whole field, since many downstream analyses (PRS, fine-mapping, LDSC) operate on summary stats rather than raw genotypes.

> **EE framing — GWAS as massively-multichannel detection**: A GWAS is a bank of $M \approx 10^6$ parallel single-channel detectors, one per SNP, running a classical likelihood-ratio test (t-statistic / Wald / logistic z-score) and reporting the p-value at each channel. The genome is the frequency axis of a coarse filter bank; each SNP is one filter; the phenotype is the signal being correlated against each channel. You then have to control the genome-wide false-alarm rate across all $M$ channels — a multi-hypothesis detection problem, identical in structure to multi-band radar, multi-target tracking, or multi-pixel anomaly detection in imaging. The two defining problems — how to control genome-wide FAR (Part 3), and how to handle correlated channels via LD (Parts 2, 6) — are familiar problems in disguise.

### 1.4 What a successful GWAS looks like (≈4 min)

When everything works:

- Hundreds to thousands of SNPs exceed the genome-wide significance threshold.
- Peaks cluster into dozens to hundreds of independent loci (because nearby SNPs are correlated by LD).
- Effect sizes are uniformly small: $\beta$ in units of the trait's SD is typically 0.01–0.05 for common variants; odds ratios (OR) for binary traits are usually 1.02–1.2. The days of Mendelian-sized effects are gone — modern GWAS finds **many small contributors** to polygenic traits.
- Heritability explained (Part 5) is often modest — only 10–30% of the genetic variance attributable to the detected SNPs, the rest being "missing heritability" (Part 5).

Key cohorts to know:

- **UK Biobank** (~500k participants; deep phenotyping; permissive access).
- **FinnGen** (~500k Finns; founder-population-specific variants).
- **Biobank Japan, China Kadoorie**, **All of Us**, **Million Veteran Program**, **H3Africa** — filling out global representation.

> **Historical pointer**: The first widely-cited GWAS was Klein et al. 2005, testing 116,204 SNPs in 96 cases (age-related macular degeneration) + 50 controls, finding a hit in the complement factor H (CFH) gene. The WTCCC 2007 paper is the field's "big bang" moment — 500k SNPs × 17,000 individuals × 7 diseases, explicitly demonstrating that common diseases have many common-variant contributions and that genome-wide scans find them. Nearly every methodological standard used today — genome-wide significance at 5×10⁻⁸, Manhattan plot presentation, QQ plot inflation diagnostics, PCA-based stratification correction — was established or consolidated in the 2005–2010 window.

---

## Part 2 — The Manhattan Plot and What It Shows (≈25 min)

### 2.1 Manhattan plots (≈8 min)

A **Manhattan plot** is the canonical GWAS visualisation: one dot per SNP, x-axis is genomic position across all 22 autosomes plus X, y-axis is $-\log_{10}(p)$. Alternating chromosome colours make the boundaries visible. Peaks rising far above the horizontal significance line are the "hits."

Why it's called "Manhattan": the skyscrapers.

What to read off it:

- **Number of peaks above $5 \times 10^{-8}$**: count of independent hits.
- **Peak heights**: strength of the association. Very tall peaks mean either very large effects or very large sample sizes.
- **Peak shapes**: a tall column with neighbouring SNPs also elevated signals LD-driven structure — the peak is driven by a few causal variants propagating their signal to correlated neighbours (Part 6, fine-mapping).
- **Chromosome patterns**: clustering of hits on specific chromosomes for polygenic traits often reflects chromosome-length and gene-density effects.

**FIGURE — Figure #3: Manhattan plot of a well-powered GWAS** → `diagrams/lecture-13/03-manhattan-plot.svg`
*Manhattan-style plot. X-axis: 22 autosomes alternating colour. Y-axis: −log10(p) from 0 to 30. Horizontal dashed line at 7.3 (5×10⁻⁸). Three or four clear peaks above the line annotated with nearby gene names (e.g. FTO, TCF7L2). Several near-significant bumps. Axis label suggests ~1 million SNPs tested, ~300k participants.*

### 2.2 QQ plots and inflation (≈8 min)

A **QQ (quantile-quantile) plot** compares observed vs expected $-\log_{10}(p)$ values. Under the null (no real signals, no confounding), the observed quantiles should match the expected ones — a straight 45° line.

Two kinds of deviation:

- **Late departure (upper right)**: the top-right of the plot lifts above the diagonal. This is what you want to see — the strongest hits are real signals, exceeding what you'd see under pure null.
- **Early departure (origin outwards)**: even small-p SNPs are above the diagonal. This is **inflation**: the whole distribution of test statistics is shifted up. Causes: population stratification (Part 4), cryptic relatedness, genotyping-batch effects, incorrect covariance modelling.

The **genomic control inflation factor** $\lambda_{GC}$ quantifies this:

$$\lambda_{GC} = \frac{\text{median}(\chi^2_{\text{observed}})}{\text{median}(\chi^2_{\text{null}})} = \frac{\text{median}(\chi^2_{\text{observed}})}{0.456}$$

where the denominator is the median of a $\chi^2_1$ distribution. $\lambda_{GC} = 1$ means no inflation. $\lambda_{GC} > 1.05$ suggests confounding. $\lambda_{GC} > 1.2$ is a red flag. But note: for polygenic traits with real widespread signal, $\lambda_{GC}$ also goes up above 1 because of genuine genetic signal, not just confounding. Modern practice: use **LD Score regression** (Part 5) to partition $\lambda_{GC}$ into "real polygenic signal" and "true confounding" components.

**FIGURE — Figure #4: QQ plot with inflation** → `diagrams/lecture-13/04-qq-plot.svg`
*Two stacked QQ panels. Top: observed vs expected −log10(p); straight 45° line plus observed dots. Left panel shows a well-controlled GWAS — dots follow the line until upper right, then lift. Right panel shows an inflated GWAS — dots lift above the line starting near the origin. λ_GC annotated: 1.02 (good) vs 1.25 (inflated). Caption: "inflation diagnoses confounding".*

### 2.3 Interpreting peak structure via LD (≈6 min)

A single SNP with a true effect isn't alone on a Manhattan plot. It sits in an LD block (Lecture 12) with many correlated SNPs; each of those neighbours also exceeds significance because their genotypes are correlated with the causal variant's. The "peak" around a true hit typically spans 10 kb to 1 Mb, consisting of tens to hundreds of sub-significant SNPs with only one or two causal.

This is both a blessing and a curse:

- **Blessing**: robust hits stand out — the peak shape (cluster of hits) is harder to fake than a single isolated tall SNP.
- **Curse**: you can't tell from the Manhattan plot which SNP is causal. The "lead SNP" (highest −log10(p)) is often not the causal variant; it's the best tag in your genotyping / imputation panel. Fine-mapping (Part 6) tries to narrow down.

Practical consequence: the first-round output of a GWAS isn't usually "this SNP causes this disease." It's "this locus (a ~100kb window) harbours a variant affecting this disease." Going from locus to causal variant is a separate analysis.

### 2.4 Manhattan plots as spectra (≈3 min)

> **EE framing — Manhattan plot as genome-wide periodogram**: Lay the genome along an axis; at each position compute a test statistic and its p-value. The Manhattan plot is a 1D "spectrum" of signal strength along the genome coordinate. Peaks above threshold are detections; the clutter below is the null. LD-driven peak width is analogous to **filter bandwidth** — each causal variant's signal leaks into correlated neighbours, broadening its footprint on the spectrum. Wider LD blocks = broader peaks = fewer effective independent tests (Part 3). Just as a coarse periodogram trades resolution for noise floor, a GWAS trades per-SNP precision for genome-wide coverage.

> **Intuition box**: The Manhattan plot is to GWAS what a spectrogram is to an audio signal. You're looking for peaks above the noise floor. A well-powered study with genuine hits shows a handful of clear peaks over a flat background. A study that's too small or too confounded shows either nothing clear (too small) or a noisy haze with ambiguous bumps (confounded). Learning to read Manhattans is most of the skill needed to assess a GWAS paper.

---

## Part 3 — Multiple Testing at SNP Scale (≈25 min)

### 3.1 Why $5 \times 10^{-8}$ (≈8 min)

If you run $M$ independent tests at per-test significance level $\alpha$, the expected number of false positives is $M \alpha$. Under **Bonferroni correction**, you demand $\alpha_{\text{adjusted}} = \alpha / M$ to keep the family-wise false-positive rate below $\alpha$.

For a 0.05 FWER threshold and $M = 1{,}000{,}000$ independent tests:

$$\alpha_{\text{adjusted}} = 0.05 / 10^6 = 5 \times 10^{-8}$$

That's the origin of the $5 \times 10^{-8}$ genome-wide threshold. The key number is $10^6$ — the **effective number of independent tests** across the genome. Empirically (Risch & Merikangas 1996; Pe'er et al. 2008), the number of LD-independent SNPs in European populations is about $10^6$, despite there being ~10M SNPs after imputation. Similar exercises in African populations give ~1.5×10⁶ because of shorter LD blocks; LD-based test-count estimation is population-specific.

Why not use all 10M imputed SNPs as the test count? Because they're **not independent** — adjacent SNPs are in LD, so testing each of them gives correlated p-values. Using the physical SNP count in Bonferroni would over-correct dramatically; using 1 million (the LD-independent count) is right.

Modern practice occasionally adjusts the threshold for specific populations or whole-genome sequencing studies (where rare variants introduce more effectively independent tests), but $5 \times 10^{-8}$ remains the de facto standard.

**FIGURE — Figure #5: Multiple-testing and effective independent tests** → `diagrams/lecture-13/05-effective-tests.svg`
*Three stacked panels. Top: a cartoon genome coordinate with SNP tick marks very densely packed (~10 million). Middle: LD-merged blocks (~1 million) — same SNPs grouped by LD. Bottom: Bonferroni threshold line for α=0.05 against M=10⁶ → p=5×10⁻⁸. Annotation: "effective count is LD-independent blocks, not raw SNP count".*

### 3.2 FWER vs FDR (≈7 min)

Two different error-control philosophies:

- **Family-Wise Error Rate (FWER)**: probability of any false positive across all tests. Controlled by Bonferroni or its variants (Holm, Hochberg). Strict and conservative.
- **False Discovery Rate (FDR)**: expected fraction of discoveries that are false. Controlled by the Benjamini-Hochberg procedure. Lets more real hits through at the cost of a known contamination rate.

GWAS traditionally uses FWER (Bonferroni at $5 \times 10^{-8}$). Reason: the field wants near-certainty at each reported hit, because each hit becomes a biological hypothesis that downstream work will spend years investigating. A 5% FDR would mean 1 in 20 hits is bogus — too noisy for the downstream fine-mapping, functional validation, and clinical interpretation to absorb.

Some sub-fields (e.g. transcriptome-wide association studies of gene expression) use FDR instead. The choice depends on what downstream tolerance for false positives looks like.

### 3.3 Power calculations (≈6 min)

For any GWAS, **power** depends on four quantities:

1. Effect size $\beta$ (smaller = harder to detect).
2. Allele frequency $p$ (lower = harder because fewer informative individuals).
3. Sample size $N$.
4. Significance threshold $\alpha$.

The test statistic for an additive model (in the simple case) is roughly:

$$\chi^2 \approx N \cdot 2p(1-p) \cdot \beta^2 / \sigma^2$$

(derivation: $\beta^2$ / variance of estimated $\beta$; variance decays as $1/(N \cdot 2p(1-p))$.)

To detect $\beta = 0.02$ SD at $p = 0.2$ (common variant) with $\alpha = 5 \times 10^{-8}$ and power 0.8, you need:

$$N \approx \frac{\chi^2_{\alpha, \text{power}}}{2p(1-p) \beta^2 / \sigma^2} \approx \frac{29}{2 \cdot 0.2 \cdot 0.8 \cdot 0.02^2} \approx 226{,}000$$

Hundreds of thousands of samples for a single $\beta = 0.02$ variant. This is why GWAS cohort sizes have been racing upward and why the field has converged on biobank-scale consortia.

For rare variants ($p = 0.001$, $\beta = 0.1$), power at $N = 10^5$ is ~3%. Single-variant testing on rare variants is underpowered; Part 7 covers the alternative.

### 3.4 Callback thread (≈4 min)

> **Intuition box**: This is the fourth time in the course you've seen controlled false-alarm rate at scale. Lecture 4 (variant calling): FDR on per-site variant calls. Lecture 6 (expression): FDR on differential-expression genes. Lecture 9 (peak calling): FDR on ChIP/ATAC peaks. Lecture 12 (selection scans): FDR on neutrality test statistics. Now: GWAS at 10⁶ tests. Each domain has the same mathematical structure (null distribution + test statistic + multi-hypothesis correction) but different null distributions, different notions of "independent tests," and different tolerances. Pay attention to *what defines the effective number of independent tests* in each — LD blocks (GWAS), gene co-expression (differential expression), peak spacing (ChIP), coalescent null (selection). The statistical framework is universal; the per-domain null characterisation is the expertise.

> **Warning box**: Genome-wide significant does **not** mean "true." It means "unlikely under the pure-null hypothesis of zero effect and no confounding." In practice, confounding from population stratification, genotyping batch effects, phenotype measurement error, and survivorship bias all produce genome-wide significant hits that aren't causal. Replication in an independent cohort, PheWAS checks, and mechanistic follow-up are required to call a hit real. $5 \times 10^{-8}$ is a necessary condition, not a sufficient one.

---

## Part 4 — Confounders and Corrections (≈40 min)

### 4.1 Population stratification (≈10 min)

The biggest confounder in GWAS is **population stratification**: systematic allele-frequency differences between groups that also differ in phenotype for non-genetic reasons.

Classic toy example: if your cohort is 60% North European and 40% South European, and the trait is something cultural (say, a dietary preference), then SNPs whose allele frequencies differ between the two groups will "predict" the trait — through ancestry, not genetics. The association is real (allele counts correlate with phenotype) but causally spurious (the actual driver is ancestry).

Stratification is the default state of any non-randomised cohort. It can never be eliminated; it can only be corrected. Without correction, a GWAS inflates across the board (every SNP with inter-group allele-frequency difference gets a bump in the test statistic) — visible as QQ-plot early departure and $\lambda_{GC} > 1$.

**FIGURE — Figure #6: Population structure in PCA space** → `diagrams/lecture-13/06-pca-structure.svg`
*Two-panel scatter plot of PC1 vs PC2 of a GWAS cohort's genotype matrix. Left: ancestry-coloured points forming visible population clusters (European, East Asian, African, admixed American). Right: the same PCs with phenotype overlay — showing how phenotype tracks ancestry. Annotation: "PC1 and PC2 capture most population structure; including them as covariates corrects for ancestry confounding".*

### 4.2 PCA correction (≈10 min)

The standard first-line correction: compute principal components of the genotype matrix (SNP × individual) and include the top $k$ PCs as covariates in every per-SNP regression.

Algorithm:

1. Centre and scale the genotype matrix. Optionally prune it to common, unlinked SNPs.
2. Compute the top $k$ eigenvectors of the genotype covariance matrix (typically $k = 5$–$20$).
3. For each individual, extract their loadings on these PCs.
4. Include them as covariates in the per-SNP regression.

The intuition: the top PCs capture population structure (inter-individual variation in allele frequencies that has a low-rank structure). Conditioning on them removes the between-group signal and keeps the within-group genetic variation.

How many PCs? Typically 10 for homogeneous cohorts, 20 for heterogeneous ones (UK Biobank uses 40 for some analyses). Too few leaves residual stratification; too many eats statistical power (each PC costs one degree of freedom).

Modern tools compute PCs efficiently from millions of SNPs × hundreds of thousands of individuals (PLINK2, FastPCA, TeraPCA). For biobank scales, specialised randomised-SVD approaches are standard.

> **EE framing — PCA correction as noise whitening**: A GWAS with population structure is a detection problem where the "noise" (per-SNP test statistic under the null) has a non-white covariance — test statistics at SNPs with similar allele-frequency differentiation patterns are correlated. Top PCA components of the genotype matrix are the principal eigen-directions of this noise covariance. Including them as covariates is equivalent to projecting out those directions — exactly **noise whitening** from adaptive-array signal processing or matched-filter detection against non-white noise. The residual test statistic is then closer to white noise, restoring valid per-SNP inference. Without whitening, the detector (Wald test) is miscalibrated; with whitening, it matches the assumed null distribution.

### 4.3 Cryptic relatedness and kinship (≈7 min)

A second confounder: **cryptic relatedness**. Even in nominally-unrelated cohorts, by chance some pairs of participants are cousins or more distant relatives. Relatives share alleles beyond what's expected by chance; failure to account for this inflates test statistics at correlated SNPs.

Solution: compute a **kinship matrix** $\mathbf{K}$ (or Genetic Relationship Matrix, GRM) capturing pairwise genetic similarity across all participants. Entry $K_{ij}$ is roughly the fraction of alleles shared identical-by-descent between individuals $i$ and $j$. For unrelated pairs, $K_{ij} \approx 0$. For parent-child, $K_{ij} \approx 0.5$. For second cousins, $K_{ij} \approx 0.03$.

**FIGURE — Figure #7: Kinship matrix heatmap** → `diagrams/lecture-13/07-kinship-matrix.svg`
*A kinship-matrix heatmap (N individuals × N individuals) with diagonal ~0.5 (self-kinship), small off-diagonal blocks showing clusters of related pairs (family groups). Colour scale from white (unrelated) to deep cobalt (first-degree relative). Annotation: "first-degree relatives at 0.5; unrelateds at ~0".*

### 4.4 Mixed-model regression (≈10 min)

The modern workhorse for GWAS at biobank scale is **mixed-model regression**: fit each SNP with both fixed (per-SNP effect, covariates) and random (individual-level residuals with kinship-structured covariance) effects.

Model for a quantitative trait:

$$\mathbf{Y} = \mathbf{X} \boldsymbol{\gamma} + G \beta + \mathbf{u} + \boldsymbol{\varepsilon}$$

where $\mathbf{u} \sim \mathcal{N}(0, \sigma_g^2 \mathbf{K})$ is a random effect with **kinship-structured covariance** (captures relatedness and ancestry together), and $\boldsymbol{\varepsilon} \sim \mathcal{N}(0, \sigma_e^2 \mathbf{I})$ is residual noise.

This is a **generalised least squares** problem. The per-SNP estimate of $\beta$ is:

$$\hat{\beta} = (G^T \mathbf{V}^{-1} G)^{-1} G^T \mathbf{V}^{-1} \mathbf{Y}$$

where $\mathbf{V} = \sigma_g^2 \mathbf{K} + \sigma_e^2 \mathbf{I}$ is the total-covariance matrix. Inverting $\mathbf{V}$ for $N = 500{,}000$ is a $500{,}000^2$-entry matrix inversion — impractical directly. Modern tools use algorithmic shortcuts:

- **BOLT-LMM** (Loh et al. 2015): iterative methods + spectral tricks; scales to N > 500k. Default for quantitative traits.
- **SAIGE** (Zhou et al. 2018): saddle-point approximation for case-control imbalance. Default for imbalanced binary traits (rare diseases).
- **REGENIE** (Mbatchou et al. 2021): two-stage ridge regression + single-SNP test; scales to UK Biobank's 500k × 10M SNPs. Emerging default.

Mixed models are often shown to reduce inflation ($\lambda_{GC} \to 1$) while preserving real signal, achieve higher power than PCA-only correction (by capturing more structure), and handle imbalanced case-control cohorts correctly.

> **EE framing — mixed models as GLS with structured covariance**: A mixed-model GWAS is generalised least squares with a specific covariance structure: total-residual covariance = genetic-kinship component + environmental-noise component. The kinship matrix $\mathbf{K}$ encodes who is related to whom (and how much). The variance components $\sigma_g^2, \sigma_e^2$ weight the two components. Standard LS assumes $\mathbf{I}$ covariance — white noise; the test statistic is then wrong whenever residuals are correlated across individuals. GLS with $\mathbf{V}^{-1}$ as the inverse-covariance weight is the optimal linear unbiased estimator in this setting (Gauss-Markov). BOLT-LMM, SAIGE, REGENIE are efficient approximations to this GLS for biobank-scale data.

**FIGURE — Figure #8: Mixed-model vs OLS inflation comparison** → `diagrams/lecture-13/08-lmm-vs-ols.svg`
*Two QQ panels side-by-side. Left: OLS regression, λ_GC = 1.15, clear early inflation. Right: mixed-model (LMM) regression on the same data, λ_GC = 1.02, inflation resolved. Annotation highlights the same GWAS run both ways: LMM restores proper calibration.*

### 4.5 Batch effects and technical confounders (≈3 min)

Beyond genetics, technical confounders are widespread:

- **Genotyping batch**: different plates, machines, or time periods produce systematically different genotype calls. Include batch indicator as a covariate.
- **DNA source** (blood vs saliva): occasionally produces SNP-specific biases.
- **Imputation quality**: imputed SNPs with low confidence scores inflate the test.
- **Phenotype measurement heterogeneity** (multi-centre studies): include centre indicators.

Every modern GWAS includes several technical covariates beyond the genetic ones. Documentation is ritualistic — published GWAS protocols list 10–40 covariates.

> **Warning box**: PCA correction is not always sufficient. For recently-admixed populations (Hispanic/Latino, African American cohorts), a handful of PCs miss fine-scale structure. Mixed models with the full kinship matrix usually correct better. For populations with strong founder effects (Finns, Ashkenazi), a smaller number of PCs plus mixed models is standard. Always inspect the QQ plot and λ_GC **after** correction; if still inflated, you're probably missing structure.

---

## Part 5 — Beyond Single-SNP: Heritability and Polygenic Risk (≈30 min)

### 5.1 Heritability (≈7 min)

**Heritability** $h^2$ is the fraction of phenotypic variance explained by genetic variation:

$$h^2 = \sigma_g^2 / \sigma_p^2$$

where $\sigma_g^2$ is genetic variance, $\sigma_p^2$ is total phenotypic variance. Two flavours:

- **Broad-sense heritability** $H^2$: all genetic effects (additive + dominance + epistasis).
- **Narrow-sense heritability** $h^2$: additive genetic effects only.

GWAS typically targets $h^2$ (additive). Classical twin studies estimate $h^2$ from the comparison of MZ and DZ twin correlations. Modern SNP-heritability estimates use GWAS data directly:

- **GCTA** (Yang et al. 2011): computes $h^2$ from the kinship matrix + phenotype via REML.
- **LDSC** (Bulik-Sullivan et al. 2015): computes $h^2$ from summary statistics alone, exploiting the relationship between a SNP's LD score and its expected χ² under polygenic signal.

Both are widely used; LDSC is the favourite when raw genotypes aren't available.

### 5.2 LD Score regression (≈8 min)

**LDSC** is the most influential summary-statistic method of the 2010s. The idea: under polygenic signal, a SNP's expected $\chi^2$ statistic should scale with how many other SNPs it's in LD with (its "LD score"). Specifically:

$$\mathbb{E}[\chi^2_j] = 1 + N h^2 \ell_j / M + N a$$

where $\ell_j$ is the LD score of SNP $j$ (sum of $r^2$ with all neighbouring SNPs), $N$ is sample size, $M$ is total SNP count, and $a$ captures confounding.

Regressing $\chi^2$ on LD score across all SNPs gives:

- **Slope** $\propto h^2$ (the polygenic signal).
- **Intercept** $\approx 1 + a$ (the confounding offset).

This separates real polygenic signal from confounding inflation — the key diagnostic $\lambda_{GC}$ alone can't do. A GWAS with high $\lambda_{GC}$ and high LDSC slope but intercept near 1 is polygenic, not confounded. A GWAS with high $\lambda_{GC}$ and intercept well above 1 is confounded.

LDSC also estimates **genetic correlation** $r_g$ between pairs of traits from cross-product summary statistics, and has been extended to partition heritability across functional categories (stratified LDSC, baseline-LD v2 model).

> **EE framing — LDSC as spectral decomposition**: LDSC is a linear regression of per-SNP test statistic on LD score. The LD score $\ell_j$ is the $j$th diagonal entry of $\mathbf{R} \mathbf{R}^T$ where $\mathbf{R}$ is the SNP-SNP LD matrix. The regression slope is a scalar summary of how polygenic signal spreads across the LD-eigenvalue spectrum — loosely, a spectral decomposition of the kinship-driven covariance. In the EE vernacular: "spectral density of the phenotypic signal projected on the genotype's covariance eigenbasis." The intercept isolates out confounding (a uniform offset across the spectrum) from the slope (the spectral footprint of true polygenic signal).

### 5.3 Missing heritability (≈5 min)

Heritability estimates for polygenic traits (height, IQ, schizophrenia) are typically $h^2 \approx 0.5$–$0.8$ from twin studies. But summing up the effects of GWAS-significant hits explains only a fraction — perhaps 10–40% — of the variance. The rest is "missing."

Where is it?

- **Many small effects below genome-wide significance**. Individual SNPs don't cross $5 \times 10^{-8}$ but contribute collectively. PRS (§5.4) captures these implicitly; SNP-heritability estimates (GCTA, LDSC) show that the "chip heritability" (total additive variance captured by common SNPs) is ~50–80% of twin-study heritability.
- **Rare variants**. Common-SNP arrays miss low-frequency variation; sequencing + burden tests (Part 7) recover some.
- **Non-additive effects** (dominance, epistasis). Small in practice; don't explain much.
- **Structural variants**. Under-captured by SNP arrays; better with long-read sequencing (Lecture 11) and pangenome references.
- **Gene × environment interactions**. Hard to estimate; some missing heritability hides here.

Current consensus: for most common diseases, genuine common-variant polygenic contribution is large (70–80% of heritability is capturable in principle) but individual variants are so weak that power-limited GWAS hasn't found them all yet. Larger cohorts progressively recover more.

### 5.4 Polygenic Risk Scores (≈10 min)

A **polygenic risk score (PRS)** aggregates signal across many SNPs into a single per-individual score. For individual $i$:

$$\text{PRS}_i = \sum_j \hat{\beta}_j G_{ij}$$

where $\hat{\beta}_j$ is the GWAS-estimated effect of SNP $j$ and $G_{ij}$ is the individual's allele count. Variants include "clumping + thresholding" (C+T; only include SNPs at $p < \alpha$, pruned for LD), Bayesian shrinkage methods (LDpred, PRS-CS, SBayesR), and deep-learning extensions.

PRS is becoming clinically relevant:

- Coronary artery disease: top 1% of PRS has ~3× the population risk — comparable to rare Mendelian variants for this trait.
- Breast cancer: PRS identifies a ~20% of women in whom lifetime risk approaches BRCA1-carrier risk.
- Psychiatric traits: PRS-schizophrenia accounts for a meaningful fraction of within-family variance.

**The portability problem**: A PRS trained in European individuals performs ~3–5× worse in non-European ancestry. Reasons:

1. **Allele frequencies differ** — SNPs common in Europeans are rare or absent in Africans; those SNPs contribute nothing in the wrong population.
2. **LD patterns differ** — the "tag SNPs" in the training cohort aren't tagging the same causal variants in a different ancestry.
3. **Effect-size heterogeneity** — real causal variants have different effects across ancestries (gene-environment interactions, epistasis with background).

The practical and ethical consequences are severe. PRS-based clinical tools (e.g. breast-cancer screening risk calculators) trained exclusively on European data underperform in minority populations, potentially worsening health disparities. Resolving this requires **equally-sized GWAS in each ancestry**, which doesn't yet exist for most traits. The H3Africa consortium, Biobank Japan, PAGE (US), and TOPMed are building this inventory.

**FIGURE — Figure #9: PRS portability across ancestries** → `diagrams/lecture-13/09-prs-portability.svg`
*Three stacked violin plots of PRS distribution separated by population (EUR / AFR / EAS) for a PRS trained in EUR. X-axis: PRS value. Annotated: R² of PRS vs trait in each population. EUR: R² = 0.12 (the training case). AFR: R² = 0.03. EAS: R² = 0.05. Annotation: "PRS loses ~3–5× predictive power in non-training ancestry".*

**EMBED — Artifact #2: Manhattan + QQ Plot Inspector** → `artifacts/lecture-13/02-manhattan-qq.html`
*Simulate a GWAS with tunable polygenicity, stratification strength, and N. Render Manhattan + QQ; compute λ_GC; diagnose inflation vs polygenic signal. Target aha: the QQ plot distinguishes inflation from real polygenic signal — you need it alongside the Manhattan.*

**EMBED — Artifact #3: PRS Portability Simulator** → `artifacts/lecture-13/03-prs-portability.html`
*Train a PRS on one simulated ancestry; test on another with different LD and allele-frequency structure. See R² drop. Target aha: portability loss is dominated by differing LD patterns, not by allele-frequency absence alone.*

---

## Part 6 — Fine-Mapping and Colocalization (≈30 min)

### 6.1 From peak to causal variant (≈5 min)

A GWAS peak is not a variant; it's a locus containing tens to hundreds of correlated SNPs, one or a few of which is causal. **Fine-mapping** is the inverse problem of identifying the causal variant(s) from the peak.

This is hard because:

- SNPs in LD have highly correlated test statistics. The most-significant SNP often isn't the causal one; it's the closest tag in the array / imputation panel.
- Multiple causal variants can coexist in a locus — the sum of their LD-propagated signals is the visible peak.
- Rare causal variants in LD with more common tagging SNPs mis-attribute signal to the common tag.

Classical fine-mapping used conditional analysis: regress on the top SNP, re-test all others, identify whether any remaining peak survives. Modern approaches are statistical (Bayesian) and return a **credible set** — a set of SNPs guaranteed to contain the causal variant with specified posterior probability.

### 6.2 SuSiE, FINEMAP, CAVIAR (≈10 min)

**SuSiE** (Sum of Single Effects; Wang et al. 2020) reformulates fine-mapping as sparse Bayesian regression:

$$\mathbf{Y} = \sum_{l=1}^{L} \sum_{j=1}^{M} \gamma_{lj} X_j \beta_{lj} + \boldsymbol{\varepsilon}$$

where $L$ is the assumed number of causal signals per locus, and $\gamma_{lj}$ is a one-hot indicator of which SNP drives signal $l$. Fit via variational EM. Output: a **posterior inclusion probability (PIP)** for each SNP, and one **credible set per signal** (a set of SNPs whose total PIP > 0.95 for that signal).

**FINEMAP** (Benner et al. 2016): shotgun stochastic search over configurations of causal SNPs (0, 1, 2, … simultaneously causal). Posterior over configurations; credible sets derived.

**CAVIAR** (Hormozdiari et al. 2014): early, elegant Bayesian formulation. Assumes a maximum number of causal variants and enumerates configurations.

All three require the LD matrix $\mathbf{R}$ as input — typically computed from the GWAS cohort itself or a reference panel (1000 Genomes, UK Biobank). They need summary statistics $(\beta_j, \sigma_j^2)$ and LD; they don't need individual-level genotype data.

Typical output for a well-fine-mapped locus:

- 1–5 independent signals.
- Each signal's credible set contains 1–20 SNPs.
- The lead SNP (highest PIP) is often the causal variant, especially at the 1-SNP credible sets.
- When the credible set is large (say 100 SNPs), the locus is under-fine-mapped — more data is needed, or the causal variant isn't in the imputation panel.

**FIGURE — Figure #10: Fine-mapping credible set** → `diagrams/lecture-13/10-fine-mapping.svg`
*Locus zoom showing ~200 SNPs at a GWAS hit. Top track: −log10(p) per SNP with the peak visible. Middle track: LD (r²) to the lead SNP encoded as colour gradient. Bottom track: SuSiE posterior inclusion probability per SNP; credible set (5 SNPs with total PIP > 0.95) highlighted. Annotation: "credible set narrows 200 peak SNPs down to 5 candidates".*

> **EE framing — fine-mapping as sparse inverse problem**: At a GWAS locus, the observed test-statistic vector is the LD matrix times the sparse causal-effect vector, plus noise: $\mathbf{z} = \mathbf{R} \boldsymbol{\beta}_c / \sqrt{N} + \boldsymbol{\varepsilon}$ where $\boldsymbol{\beta}_c$ has a handful of non-zero entries (the causal variants). Recovering $\boldsymbol{\beta}_c$ from $\mathbf{z}$ given $\mathbf{R}$ is **sparse linear inverse problem** — exactly compressed sensing / basis pursuit / LASSO / L0 regularised reconstruction. SuSiE's "sum of single effects" is an explicit sparse decomposition; FINEMAP's model-averaging is Bayesian sparse recovery. LD-block shape (autocorrelation of $\mathbf{R}$) determines whether the inverse is well-posed. The genomics community discovered this framework independently of the signal-processing community; the underlying math is identical.

### 6.3 Colocalization (≈8 min)

A GWAS signal says *a* variant in the locus affects the trait. It doesn't say *which* gene or *how*. Colocalization adds a mechanistic layer: **does the GWAS signal match the signal from an expression QTL (eQTL)?**

If yes, the candidate gene is likely the one whose expression mediates the GWAS signal.

**coloc** (Giambartolomei et al. 2014) tests five hypotheses at each locus:

- $H_0$: no association with either trait.
- $H_1$: association with trait 1 only.
- $H_2$: association with trait 2 only.
- $H_3$: both associated, but different causal variants.
- $H_4$: both associated, same causal variant.

Posterior probabilities of each are computed from summary statistics + LD. $P(H_4) > 0.8$ is the standard threshold for claiming colocalization.

Applications:

- GWAS + eQTL → gene mediating disease.
- GWAS + pQTL (protein QTL) → protein intermediate.
- GWAS + sQTL (splicing QTL) → splicing isoform causal.
- GWAS + metabolite QTL → metabolic mediator.

Large eQTL datasets (GTEx, eQTLGen, CommonMind) make colocalization routine. Typical success rate: ~20–40% of GWAS loci colocalise with at least one eQTL tissue.

**FIGURE — Figure #11: Colocalization of GWAS and eQTL** → `diagrams/lecture-13/11-colocalization.svg`
*Two locus-zoom panels stacked at the same genomic coordinates. Top: GWAS −log10(p) for a trait (e.g. LDL). Bottom: eQTL −log10(p) for a gene (e.g. SORT1) in liver tissue. Both panels show peaks at the same position. Annotation: "P(H4) = 0.93 — GWAS signal colocalises with SORT1 eQTL → SORT1 is a strong candidate causal gene".*

**EMBED — Artifact #4: Fine-Mapping Credible Set Viewer** → `artifacts/lecture-13/04-fine-mapping.html`
*Simulate a locus with tunable LD structure and 1–3 planted causal variants. Run SuSiE-like credible-set inference. Show how LD block structure and sample size determine credible-set size. Target aha: strong LD merges signals into wide credible sets; breaking LD with more recombination or larger N shrinks them.*

### 6.4 Fine-mapping as spectral inversion (≈3 min)

> **Intuition box**: Colocalization asks a visually simple question — do the GWAS trace and the eQTL trace for the candidate gene peak at the same place along the genome? If yes, a single variant plausibly drives both; the gene whose eQTL co-peaks is the causal candidate. The coloc machinery formalises this with Bayesian posteriors because two traces peaking near each other can either share a causal variant ($H_4$) or have nearby-but-distinct causal variants ($H_3$), and LD makes these visually indistinguishable. The mental picture is: two 1D signals lined up on the genome axis — you want peak coincidence, not just peak proximity.

### 6.5 When fine-mapping fails (≈4 min)

Fine-mapping fails when:

- **Credible sets remain huge** (>50 SNPs): locus has many equally-good candidates, often because the causal variant isn't in the panel or because LD is very tight. Improving the panel (WGS imputation, denser reference) shrinks the set.
- **Multiple signals present but poorly separated**: SuSiE with $L = 5$ finds 5 signals, but credible sets overlap heavily. Usually needs more data.
- **Causal variant is structural**: fine-mapping assumes SNP-based signal. If the causal is a CNV or SV not captured by the array, fine-mapping converges to a tagging SNP that's nearby the SV but isn't causal.

Realistic success rate: for well-powered loci (lead $p < 10^{-20}$ in large European cohorts), single-variant resolution is achieved at maybe 30% of loci. The rest remain candidate-set-level.

> **Warning box**: A 95% credible set contains the causal variant with 95% posterior probability only if the prior assumptions (sparsity $L$, LD matrix accuracy, causal-effect prior) are correct. In practice, the prior on $L$ is a hyperparameter; the LD matrix may be mis-estimated from a reference panel that doesn't match the GWAS cohort; and effect priors assume a functional form that may not hold. Credible sets should be read as "strong candidates," not "the causal variant is certainly in here."

---

## Part 7 — Rare-Variant Association (≈25 min)

### 7.1 Why single-SNP fails on rare variants (≈5 min)

For very rare variants ($\text{MAF} < 0.01$), single-variant GWAS is underpowered:

- Variance of $\hat{\beta}$ scales as $1 / [N \cdot 2p(1-p)]$.
- At $p = 0.001$, $2p(1-p) \approx 0.002$ — 250× smaller than at $p = 0.5$.
- Detecting $\beta = 0.1$ at $\alpha = 5 \times 10^{-8}$ with 80% power needs $N \approx 5 \times 10^{5}$ for common variants, but $N \approx 10^8$ for $p = 0.001$.

No cohort is that large. Single-variant testing of rare variants is futile.

**The workaround**: aggregate rare variants within some grouping (usually a gene) and test the group collectively for association. If a gene harbours many rare loss-of-function variants in cases but few in controls, there's signal even though no individual variant is individually significant.

### 7.2 Burden tests (≈8 min)

**Burden tests** collapse all rare variants in a gene into a single summary statistic per individual, then regress phenotype on that summary:

$$\text{BurdenScore}_i = \sum_{j \in \text{gene}} w_j G_{ij}$$

where $G_{ij}$ is genotype at rare variant $j$ in the gene and $w_j$ is a weight (often 1, or inverse-allele-frequency weighted, or functional-impact weighted).

The regression is then:

$$Y_i = \beta_0 + \beta_{\text{burden}} \cdot \text{BurdenScore}_i + \boldsymbol{\gamma}^T \mathbf{X}_i + \varepsilon_i$$

One test per gene → ~20,000 tests, not millions. Bonferroni threshold is $0.05 / 20{,}000 = 2.5 \times 10^{-6}$ — much less stringent than single-SNP.

Works great when: most rare variants in the gene have effects in the same direction (e.g. all loss-of-function). Fails when variants have mixed directions (some gain of function, some loss) — they cancel in the sum.

### 7.3 SKAT and SKAT-O (≈7 min)

**SKAT** (Sequence Kernel Association Test; Wu et al. 2011) fixes the "mixed directions" problem:

$$Q = (\mathbf{Y} - \boldsymbol{\mu})^T \mathbf{G} \mathbf{W} \mathbf{G}^T (\mathbf{Y} - \boldsymbol{\mu})$$

where $\mathbf{G}$ is the gene's genotype matrix (rare variants only) and $\mathbf{W}$ is a diagonal weight matrix. Under the null, $Q$ follows a mixture of chi-squared distributions (Davies approximation).

Intuition: SKAT is a **kernel test** — it compares the phenotype variance explained by the gene's genotype kernel against a null. Unlike the burden test, SKAT doesn't assume directional concordance; it just asks "does this gene's genotype collectively explain phenotype variance?"

**SKAT-O** (Lee et al. 2012) combines burden and SKAT adaptively: an optimally-weighted hybrid that gets the best of both depending on the underlying directional pattern of variants. Default choice for modern rare-variant GWAS.

**FIGURE — Figure #12: Burden test vs SKAT at a gene** → `diagrams/lecture-13/12-burden-skat.svg`
*A cartoon gene with ~20 rare variants, each coloured by direction (red = deleterious, green = protective, grey = neutral). Left case: "all same direction" → burden test fires. Middle case: "mixed directions" → burden cancels, SKAT detects. Right case: "mostly neutral, few large effects" → SKAT detects through kernel weighting. Table below summarising which test wins in each regime.*

### 7.4 Data sources: gnomAD and friends (≈3 min)

Rare-variant testing needs sequencing, not arrays. Exome sequencing captures ~1–2% of the genome (all protein-coding regions); whole-genome sequencing captures everything.

Key resources:

- **gnomAD** (Karczewski et al. 2020): >800k exomes + >150k genomes from healthy populations. Standard source for rare-variant allele frequencies. Every clinical genetics pipeline uses gnomAD as the "is this variant rare?" reference.
- **UK Biobank exome / WGS releases**: exomes for all 500k; WGS for ~500k (released 2022–23). Powers rare-variant burden tests across hundreds of traits.
- **TOPMed** (~150k whole genomes, US multi-ethnic): rare variants in diverse populations.

Current scale of rare-variant GWAS: UK Biobank exome BMI burden tests ran 18,000 genes × 80 phenotypes on 200k exomes; found dozens of gene-trait associations, most novel.

### 7.5 Rare-variant caveats (≈2 min)

- **Annotation-dependent**: burden tests require a definition of "rare functional variants in this gene." Definitions vary (LoF only vs LoF + missense, allele-frequency cutoffs at 0.01 vs 0.001, etc.). Sensitivity analyses across definitions are standard.
- **Genotyping / calling errors are amplified**: false-positive rare variant calls directly inflate the burden. QC is more demanding than for common variants.
- **Case-control imbalance is severe**: many cohorts have few cases; single-variant cases have few affected individuals. SAIGE-GENE / REGENIE-GENE handle this.

**EMBED — Artifact #5: Population Stratification + PCA Correction** → `artifacts/lecture-13/05-pca-correction.html`
*Simulate a two-subpopulation cohort with an ancestry-correlated phenotype and no real SNP effects. Run GWAS without correction → massive inflation. Include top PCs → λ_GC returns to 1. Target aha: PCA correction is noise whitening; without it, every SNP with inter-group allele-frequency difference fires.*

**EMBED — Artifact #6: Burden and SKAT Simulator** → `artifacts/lecture-13/06-burden-skat.html`
*Simulate a gene with user-controlled mix of LoF / missense / neutral variants at user-controlled frequencies + effect sizes. Compute burden statistic and SKAT; show which wins in which regime. Target aha: burden dominates when directions agree; SKAT dominates when mixed; SKAT-O approximates whichever is better per locus.*

**EMBED — Artifact #7: Kinship-Corrected GWAS Mini-Simulator** → `artifacts/lecture-13/07-mixed-model.html`
*Simulate a small cohort with relatedness (siblings, cousins). Run OLS GWAS and a mini mixed-model GWAS; compare QQ plots and effect estimates. Target aha: with cryptic relatedness, OLS inflates and the fix is to model the kinship covariance explicitly.*

> **EE framing — burden testing as matched filter on rare-variant count**: A burden test aggregates rare variants in a gene into one linear summary, then regresses phenotype on that summary. Equivalent to a **matched filter** where the template is "sum of rare damaging variants in gene $G$" and the signal is phenotype. SKAT is the kernel version, using a quadratic form $\mathbf{Y}^T \mathbf{K} \mathbf{Y}$ with $\mathbf{K} = \mathbf{G} \mathbf{W} \mathbf{G}^T$ — variance-based detection that doesn't require directional alignment. SKAT-O is adaptive combination of the two detectors, analogous to a **composite-hypothesis matched filter** that doesn't assume the signal template's exact form.

> **Discussion prompt**: You run a burden test at $N = 10{,}000$ exomes with 80% power to detect OR = 3 for a gene with MAF-sum 1%. No genes are significant at $\alpha = 2.5 \times 10^{-6}$. Three explanations: (a) no rare-variant signal exists for this trait, (b) the trait has rare-variant signal but your cohort is too small, (c) rare variants exist but with mixed directional effects that burden cancelled. How would you distinguish these? What additional analysis would you run, and how does scaling to $N = 100{,}000$ exomes change the picture?

---

## Wrap-up (≈10 min)

### What you should take away

- **GWAS is a bank of per-SNP regressions** with the phenotype as outcome, genotype as predictor, covariates including PCs, age, sex, and batch. At biobank scale this is a massively-multichannel detection problem.
- **Genome-wide significance is $5 \times 10^{-8}$** — Bonferroni correction for ~10⁶ LD-independent tests. The Manhattan plot displays the signal; the QQ plot checks calibration.
- **Population stratification inflates every SNP uniformly.** PCA correction (top 10–40 PCs as covariates) is the first line of defence; mixed-model regression with kinship (BOLT-LMM, SAIGE, REGENIE) is the modern workhorse.
- **LD Score regression separates polygenic signal from confounding** using summary statistics alone. Slope ∝ heritability; intercept ∝ confounding.
- **Polygenic Risk Scores aggregate signal across the genome into a per-individual score.** PRS has real clinical utility for coronary disease, breast cancer, schizophrenia — but portability across ancestries is ~3–5× worse when training data is European-dominated.
- **Fine-mapping is a sparse inverse problem.** SuSiE, FINEMAP, CAVIAR return credible sets — small groups of SNPs that collectively explain the locus signal with specified posterior probability.
- **Colocalization links GWAS loci to genes** via matched eQTL signal. coloc's $P(H_4) > 0.8$ is the standard threshold.
- **Rare-variant GWAS uses burden tests and SKAT.** Single-variant testing fails below MAF ~0.01; aggregation across a gene recovers power. gnomAD is the standard rare-variant reference.

### Next lecture

Data engineering, file formats, and reproducibility. FASTQ → BAM/CRAM; VCF/BCF; HDF5/Zarr for single-cell; Parquet for tabular omics. SRA/ENA/GEO/dbGaP data repositories. Containerization (Docker, Singularity) and the dependency-hell problem. Workflow languages (Nextflow, Snakemake, WDL). GIAB truth sets and benchmarking. The culture of bioinformatics: preprints, nf-core, open-source norms.

### Homework

1. Download the 1000 Genomes chr22 genotype data + a simulated phenotype. Run a per-SNP linear regression. Plot the Manhattan and QQ plots. Report $\lambda_{GC}$.
2. Re-run the same GWAS including 10 principal components as covariates. Report the new QQ plot and $\lambda_{GC}$. Did inflation go down?
3. Compute a PRS from a published GWAS summary statistic (e.g. for height from the Yengo et al. meta-analysis). Apply it to a held-out cohort. Report $R^2$ of PRS vs phenotype.
4. At a single well-known GWAS locus (e.g. *FTO* for BMI), run SuSiE on the summary statistics + LD. Report the credible set(s). What's the smallest credible set you can achieve?
5. Use gnomAD to check the allele frequency of three variants: rs334 (sickle cell), rs80357906 (BRCA1), and rs1815739 (ACTN3). Comment on how MAF and population structure affect each.

### Recommended reading

- Klein, R. J. et al. (2005). Complement factor H polymorphism in age-related macular degeneration. *Science* 308, 385–389. (First widely-cited GWAS.)
- Wellcome Trust Case Control Consortium (2007). Genome-wide association study of 14,000 cases of seven common diseases and 3,000 shared controls. *Nature* 447, 661–678.
- Visscher, P. M., Wray, N. R., Zhang, Q., et al. (2017). 10 years of GWAS discovery: biology, function, and translation. *American Journal of Human Genetics* 101, 5–22.
- Price, A. L., Patterson, N. J., Plenge, R. M., et al. (2006). Principal components analysis corrects for stratification in genome-wide association studies. *Nature Genetics* 38, 904–909.
- Loh, P. R., Tucker, G., Bulik-Sullivan, B. K., et al. (2015). Efficient Bayesian mixed-model analysis increases association power in large cohorts. *Nature Genetics* 47, 284–290. (BOLT-LMM.)
- Zhou, W., Nielsen, J. B., Fritsche, L. G., et al. (2018). Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. *Nature Genetics* 50, 1335–1341. (SAIGE.)
- Mbatchou, J., Barnard, L., Backman, J., et al. (2021). Computationally efficient whole-genome regression for quantitative and binary traits. *Nature Genetics* 53, 1097–1103. (REGENIE.)
- Bulik-Sullivan, B. K., Loh, P. R., Finucane, H. K., et al. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nature Genetics* 47, 291–295. (LDSC.)
- Wang, G., Sarkar, A., Carbonetto, P., &amp; Stephens, M. (2020). A simple new approach to variable selection in regression, with application to genetic fine-mapping. *Journal of the Royal Statistical Society, Series B* 82, 1273–1300. (SuSiE.)
- Giambartolomei, C., Vukcevic, D., Schadt, E. E., et al. (2014). Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genetics* 10, e1004383. (coloc.)
- Wu, M. C., Lee, S., Cai, T., et al. (2011). Rare-variant association testing for sequencing data with the sequence kernel association test. *American Journal of Human Genetics* 89, 82–93. (SKAT.)
- Lee, S., Emond, M. J., Bamshad, M. J., et al. (2012). Optimal unified approach for rare-variant association testing with application to small-sample case-control whole-exome sequencing studies. *American Journal of Human Genetics* 91, 224–237. (SKAT-O.)
- Karczewski, K. J., Francioli, L. C., Tiao, G., et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature* 581, 434–443. (gnomAD.)
- Martin, A. R., Kanai, M., Kamatani, Y., et al. (2019). Clinical use of current polygenic risk scores may exacerbate health disparities. *Nature Genetics* 51, 584–591. (PRS portability.)
- UK Biobank: <https://www.ukbiobank.ac.uk/>
- GTEx portal: <https://gtexportal.org/>
- gnomAD browser: <https://gnomad.broadinstitute.org/>
- SuSiE documentation: <https://stephenslab.github.io/susieR/>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — The GWAS Framework                                   | 25&nbsp;min | 0:25 |
| Part 2 — The Manhattan Plot and What It Shows                 | 25&nbsp;min | 0:50 |
| Part 3 — Multiple Testing at SNP Scale                         | 25&nbsp;min | 1:15 |
| Part 4 — Confounders and Corrections                            | 40&nbsp;min | 1:55 |
| Part 5 — Beyond Single-SNP: Heritability and PRS                | 30&nbsp;min | 2:25 |
| Part 6 — Fine-Mapping and Colocalization                        | 30&nbsp;min | 2:55 |
| Part 7 — Rare-Variant Association                                | 25&nbsp;min | 3:20 |
| Wrap-up                                                           | 10&nbsp;min | 3:30 |

**Total:** ~3h 30min of content.
