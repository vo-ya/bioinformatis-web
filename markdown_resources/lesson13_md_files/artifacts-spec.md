# Lecture 13 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 13 (GWAS and statistical genetics).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-13.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-13/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript.
- Must render standalone.
- Embedded in the lecture via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` exactly once near the end of `<body>`.** C6 smoke gate.

### 1.2 Visual design

- Tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- Manhattan / QQ dots: `--accent` cobalt for significant, `--fg-muted` grey for null.
- Genome-wide significance line: amber dashed at $y = 7.3$.
- Ancestry palette (consistent with Lecture 12): EUR cobalt, AFR amber, EAS red, SAS green.
- LD heatmap: sequential red from white to deep red.
- Kinship heatmap: sequential cobalt.
- Typography: Inter for UI chrome; JetBrains Mono for p-values, λ_GC values, test statistics, rs-IDs, allele counts.
- Default state is instructive: opens with a pre-computed example.
- Controls grouped in a panel above or to the left of the visualisation.
- Animations ≤ 400 ms.

### 1.3 Interaction model

- **Sliders / toggles / dropdowns** — editable parameters with sensible ranges.
- **Run / Re-simulate / Reset** — for stochastic simulations.
- **Threshold slider** — where applicable, user-adjustable significance threshold.
- Illegal input → quiet inline message (`--fg-muted`); never an `alert()`.

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question:

- Single-SNP Association Explorer → observed $\hat{\beta}$, z-statistic, p-value, power at the current $(N, p, \beta)$.
- Manhattan + QQ Inspector → computed $\lambda_{GC}$ and an interpretation label (calibrated / inflated / polygenic-with-inflation).
- PRS Portability Simulator → $R^2$ of PRS vs phenotype in training vs test ancestry.
- Fine-Mapping Credible Set Viewer → credible-set membership list + total PIP + lead-SNP PIP; comparison to ground-truth causal variant when known.
- PCA-Correction Simulator → Manhattan + $\lambda_{GC}$ before/after including top-$k$ PCs as covariates.
- Burden / SKAT Simulator → burden and SKAT p-values, directional-agreement diagnostic, winner.
- Kinship-Corrected GWAS Mini-Simulator → OLS vs mixed-model QQ plot comparison + $\lambda_{GC}$ values.

### 1.5 Feasibility gate on user input (required where input is free-form)

- All seven artifacts use slider / dropdown / numeric-input controls only (no free-form text or uploads).
- Numeric inputs clamp to their valid range and display the clamped value.

### 1.6 Pedagogical constraint

Every artifact produces its named aha moment. If the student plays with the controls and doesn't land on it, the artifact has failed.

### 1.7 Out of scope

- No accounts, no telemetry, no network calls beyond declared CDN libraries (KaTeX permitted for in-page math; otherwise none).
- No external data files > 100 KB (per-artifact hardcoded simulated cohorts).

---

## 2. Artifact #1 — GWAS Single-SNP Association Explorer

**File**: `artifacts/lecture-13/01-single-snp-test.html`
**Lecture anchor**: §1.2 Per-SNP regression
**EE framing reinforced**: per-SNP Wald test as matched filter; detectability governed by $N \cdot 2p(1-p) \cdot \beta^2$.

### Teaching purpose

Slide allele frequency $p$, effect size $\beta$, sample size $N$, and noise $\sigma$. Watch the simulated scatter of $(G, Y)$ and the null / alternative sampling distribution of $\hat{\beta}$. Compute the observed $z$-statistic, p-value, and power at $5 \times 10^{-8}$. Show why GWAS effect sizes require biobank-scale $N$.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Allele frequency p:    [──●── 0.30 ]                        │
│ Effect size β (SD):    [──●── 0.05 ]                        │
│ Sample size N:         [──●── 100000 (log) ]                │
│ Residual σ:            [──●── 1.0 ]                         │
│ Significance α:        [──●── 5×10⁻⁸ ]                      │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Scatter panel:                                              │
│   G ∈ {0,1,2} vs Y (jittered)                              │
│   regression line                                           │
├─────────────────────────────────────────────────────────────┤
│ Hypothesis panel:                                           │
│   null & alternative Gaussians for β̂                        │
│   observed β̂ marker; tail region shaded                     │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ β̂ = 0.048  SE = 0.008  z = 6.0  p = 1.9×10⁻⁹         │   │
│ │ Power at α=5×10⁻⁸: 89% (analytical)                    │   │
│ │ Try N=20,000: power drops to ~25%                       │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Allele frequency $p$: 0.01 – 0.50, default 0.30.
- Effect size $\beta$: 0 – 0.2 (in SD units), default 0.05.
- Sample size $N$: 1000 – 1,000,000 (log slider), default 100,000.
- Residual $\sigma$: 0.5 – 2.0, default 1.0.
- Significance threshold: 5×10⁻⁸, 5×10⁻⁵, or 0.05 (dropdown).
- Re-simulate.

### What they see

- Scatter of genotype vs phenotype.
- Overlay of regression line and its 95% CI.
- Null / alternative sampling-distribution panel for $\hat{\beta}$.
- Outcome banner with statistics and analytical power.

### Target aha

Start at $p = 0.30, \beta = 0.05, N = 100{,}000$ → z ≈ 6, genome-wide significant. Drop $N$ to 20,000 → z ≈ 2.7, p ≈ 0.007, **not** significant at genome-wide level despite the same effect. Drop $p$ to 0.02 at large $N$ → effect gets harder to detect. Student sees: you need **both** large $N$ and adequately common alleles, with effect-size that overcomes the $1/\sqrt{N \cdot 2p(1-p)}$ standard error.

### Technical notes

- Pure JS, seeded PRNG.
- Simulate genotypes from Binomial(2, p); phenotype $Y = \beta G + \varepsilon$, $\varepsilon \sim \mathcal{N}(0, \sigma^2)$.
- OLS fit: $\hat{\beta} = \text{cov}(G, Y) / \text{var}(G)$; $\text{SE} = \sigma / \sqrt{N \cdot \text{var}(G)}$; $z = \hat{\beta}/\text{SE}$.
- Analytical power: normal-approximation formula for two-sided test.
- KaTeX defer-loaded for math.

### Acceptance criteria

- [ ] Default (p=0.3, β=0.05, N=10⁵) → z ≈ 6, p ≈ 10⁻⁹, power > 0.8.
- [ ] Setting β=0 → z hovers near 0; p uniformly distributed on re-simulation.
- [ ] Halving N → SE increases by √2.
- [ ] Power calculation matches simulated frequency of genome-wide significance over 100 resamples.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — Manhattan + QQ Plot Inspector

**File**: `artifacts/lecture-13/02-manhattan-qq.html`
**Lecture anchor**: §2.1 Manhattan plots; §2.2 QQ plots and inflation
**EE framing reinforced**: QQ plot as a null-vs-observed calibration check; inflation as uniform channel-gain error.

### Teaching purpose

Simulate a GWAS with user-controlled (a) polygenicity (number of truly causal SNPs), (b) stratification strength (a uniform inflation added to $\chi^2$), and (c) sample size $N$. Render the Manhattan plot and the QQ plot. Compute $\lambda_{GC}$. Show that the QQ plot distinguishes real polygenic signal (late departure) from confounding (early departure).

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ # causal SNPs:          [──●── 50 ]  (polygenicity)         │
│ Average effect:         [──●── 0.04 SD ]                    │
│ Stratification σ_strat: [──●── 0.0 ]  (0 = clean)            │
│ Sample size N:          [──●── 200,000 (log) ]              │
│ SNP count M:            [──●── 500,000 ]                    │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Manhattan panel (top):                                      │
│   M dots across 22 chrom bands with threshold line          │
├─────────────────────────────────────────────────────────────┤
│ QQ panel (bottom-left):                                     │
│   observed vs expected −log10(p); 45° reference              │
├─────────────────────────────────────────────────────────────┤
│ Stats panel (bottom-right):                                 │
│   λ_GC, # genome-wide hits, # causal recovered               │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ λ_GC = 1.08.  QQ lift-off near the top (real signal).  │   │
│ │ 28 / 50 causal SNPs genome-wide significant             │   │
│ │ Add stratification σ=0.3 → λ_GC = 1.35, early inflation │   │
│ │ flag; λ_GC alone can't distinguish from polygenicity    │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Number of causal SNPs: 0 – 500, default 50.
- Average effect size: 0.01 – 0.1 SD, default 0.04.
- Stratification strength $\sigma_{\text{strat}}$: 0 – 1.0, default 0. Adds a per-SNP inflation: $z_j \to z_j + \sigma_{\text{strat}} \cdot \delta_j$ where $\delta_j$ is a per-SNP ancestry-structure term.
- Sample size $N$: 1,000 – 1,000,000 (log).
- SNP count $M$: 100,000 – 1,000,000.
- Re-simulate.

### What they see

- Manhattan plot (top).
- QQ plot (bottom-left).
- Statistics panel (bottom-right): $\lambda_{GC}$, number of hits, number of causal-SNPs recovered.
- Outcome banner with interpretation.

### Target aha

With stratification = 0 and polygenicity = 50 → QQ lifts only at the top-right; $\lambda_{GC}$ slightly above 1. Add stratification = 0.3 → QQ lifts from the origin; $\lambda_{GC}$ climbs to 1.3+. The Manhattan plot looks similar in both cases. Student sees: **the QQ plot**, not $\lambda_{GC}$ alone, tells you whether inflation is real signal or confounding.

### Technical notes

- Pure JS, seeded PRNG.
- Assign SNPs uniformly to 22 virtual chromosomes.
- Causal SNPs: pick $k$ at random; assign effects $\beta_j$; compute $z_j$ from an approximate non-central $\chi^2$ model with per-SNP SE.
- Stratification: add $\delta_j \sim \mathcal{N}(0, \sigma_{\text{strat}}^2)$ correlated across SNPs in chunks of 1000 (simulates LD-correlated stratification effect).
- $\lambda_{GC} = \text{median}(\chi^2) / 0.456$.
- Manhattan: use SVG or canvas; dots coloured by alternating chromosome.

### Acceptance criteria

- [ ] Default (no stratification, polygenic) → λ_GC near 1; QQ late lift-off.
- [ ] Setting polygenicity=0, stratification=0.3 → λ_GC > 1.2; QQ early lift-off.
- [ ] Manhattan threshold line at 7.3 visible.
- [ ] Significance hits count updates on re-simulate.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — PRS Portability Simulator

**File**: `artifacts/lecture-13/03-prs-portability.html`
**Lecture anchor**: §5.4 Polygenic Risk Scores
**EE framing reinforced**: PRS as linear combiner trained in one channel and deployed in another; portability loss as training-vs-test distribution mismatch.

### Teaching purpose

Simulate a GWAS in one ancestry (training population) and derive a PRS. Apply it to a test cohort of a second ancestry. Show the $R^2$ drop due to LD and allele-frequency divergence between populations. Student can toggle different match/mismatch scenarios.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Training population:     [ EUR ▾ ]                          │
│ Test population:         [ AFR ▾ / EAS / SAS / EUR ]        │
│ Ancestry-divergence factor: [──●── 1.0 ]                    │
│   (1 = full divergence; 0 = identical populations)         │
│ Training N:              [──●── 200,000 ]                   │
│ Test N:                  [──●── 10,000 ]                    │
│ PRS method:              [ C+T ▾ / PRS-CS / SBayesR ]       │
│ [Re-train and test]                                          │
├─────────────────────────────────────────────────────────────┤
│ Allele-frequency-difference panel:                          │
│   Per-SNP bar: |p_train - p_test|                           │
├─────────────────────────────────────────────────────────────┤
│ PRS-distribution panel:                                     │
│   Violin of PRS values per population                       │
├─────────────────────────────────────────────────────────────┤
│ Performance table:                                          │
│   R² in training pop (EUR): 0.12                           │
│   R² in test pop (AFR):     0.03                            │
│   Portability ratio: 0.25                                   │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ PRS trained EUR → tested AFR: R² drops 4×              │   │
│ │ Dominant cause: different LD patterns (85% of loss)     │   │
│ │ Allele-frequency absence contributes 15%                │   │
│ │ Try divergence=0.3 → R² matches; 1.0 → maximum drop    │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Training population (EUR / AFR / EAS / SAS), default EUR.
- Test population (defaults to different ancestry from training).
- Ancestry-divergence factor: 0 – 1.0, default 1.0 (1.0 = full divergence, 0 = identical populations).
- Training $N$: 10,000 – 500,000 (log).
- Test $N$: 1,000 – 50,000 (log).
- PRS method dropdown: C+T / PRS-CS / SBayesR (affects shrinkage but not fundamental portability pattern).
- Re-train and test.

### What they see

- Allele-frequency-difference panel per SNP.
- PRS distribution per population (violin plot).
- Performance table: $R^2$ in training vs test populations.
- Outcome banner with decomposition of portability loss.

### Target aha

Start with full EUR → AFR portability test → $R^2$ drops ~4×. Drop divergence factor to 0 (identical populations) → $R^2$ matches. Compare the two decompositions: when LD and MAF are matched, portability is preserved; when they diverge, PRS fails. Student sees that **fair deployment requires training data from the deployment population.**

### Technical notes

- Pure JS, seeded PRNG.
- Simulate training cohort with population-specific allele frequencies (draw from Dirichlet around a shared prior); LD structure generated from a simple block model with ancestry-specific block structure.
- Test cohort: same causal SNPs but with the test-population's allele frequencies and LD.
- PRS: C+T (clump top SNPs at threshold); PRS-CS (Bayesian shrinkage, simplified).
- $R^2$ computed as Pearson correlation-squared between PRS and simulated phenotype.

### Acceptance criteria

- [ ] Default EUR → AFR with full divergence → R² drops 3–5×.
- [ ] Divergence=0 → R² preserved.
- [ ] Training-cohort N scaling up improves training $R^2$ but doesn't close the portability gap.
- [ ] Decomposition (LD vs MAF) varies with divergence setting.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — Fine-Mapping Credible Set Viewer

**File**: `artifacts/lecture-13/04-fine-mapping.html`
**Lecture anchor**: §6.2 SuSiE, FINEMAP, CAVIAR
**EE framing reinforced**: fine-mapping as sparse inverse problem; credible-set size determined by LD-matrix conditioning.

### Teaching purpose

Simulate a GWAS locus with 1–3 **planted** causal variants and user-controlled LD structure. Run a simplified SuSiE-style credible-set inference. Show how LD block tightness and sample size determine credible-set size. Student toggles between strong-LD (wide credible set) and weak-LD (narrow credible set) scenarios.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ # causal variants L_true:  [──●── 2 ]                       │
│ LD block tightness:         [──●── 0.8 ]  (0=breaks,1=solid)│
│ Sample size N:              [──●── 200,000 ]                │
│ Effect size β:              [──●── 0.05 ]                   │
│ Noise σ:                    [──●── 1.0 ]                    │
│ # SNPs in locus:            [──●── 100 ]                    │
│ Fitted L_susie:             [──●── 3 ]                      │
│ [Re-simulate]  [Run SuSiE]                                  │
├─────────────────────────────────────────────────────────────┤
│ −log10(p) track (top):                                      │
│   per-SNP p-value; causal SNPs marked                        │
├─────────────────────────────────────────────────────────────┤
│ LD matrix (middle):                                          │
│   100×100 r² heatmap; block structure visible                │
├─────────────────────────────────────────────────────────────┤
│ PIP track (bottom):                                         │
│   per-SNP PIP from SuSiE; credible-set members highlighted   │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 2 credible sets found.                                  │   │
│ │ Set 1: 4 SNPs, total PIP 0.96, contains true causal 1   │   │
│ │ Set 2: 9 SNPs, total PIP 0.93, contains true causal 2   │   │
│ │ Drop block-tightness to 0.3 → each credible set shrinks │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- $L_{\text{true}}$: 1 – 3, default 2.
- LD block tightness: 0 – 1, default 0.8. Controls pairwise $r^2$ within blocks.
- Sample size $N$: 10,000 – 1,000,000 (log).
- Effect size $\beta$: 0.02 – 0.15.
- Noise $\sigma$: 0.5 – 2.0.
- SNPs in locus: 50 – 300.
- Fitted $L_{\text{susie}}$: 1 – 5, default 3.
- Re-simulate / Run.

### What they see

- $-\log_{10}(p)$ track per SNP with causals marked.
- LD heatmap.
- PIP track per SNP with credible sets highlighted.
- Outcome banner listing credible sets and their content.

### Target aha

With tight LD (0.9) → credible sets are wide (10–20 SNPs). With loose LD (0.2) → credible sets shrink to 1–3 SNPs per causal. Increase $N$ → credible sets shrink but only modestly vs reducing LD. Student sees: LD tightness fundamentally limits resolution; more data helps less than breaking LD.

### Technical notes

- Pure JS, seeded PRNG.
- Build LD matrix: diagonal block structure with within-block $r^2$ = tightness; between-block $r^2$ ≈ 0.
- Simulate summary stats: $z_j = \sum_c R_{jc} \beta_c \sqrt{N} / \sigma + \text{noise}$, where $R$ is LD and $c$ indexes causals.
- Simplified SuSiE: coordinate-ascent variational EM over $L_{\text{susie}}$ one-hot causal distributions.
- Credible set per signal: smallest set of SNPs whose cumulative PIP for that signal ≥ 0.95.

### Acceptance criteria

- [ ] Default (L_true=2, tightness=0.8) → SuSiE finds 2 credible sets; each contains one of the true causals.
- [ ] Tightness=0.2 → credible sets are narrower (typically 1–2 SNPs).
- [ ] Tightness=0.95 → credible sets merge (SuSiE can't separate causals).
- [ ] Fitted L_susie < L_true → misses one causal; L_susie > L_true → extra spurious set.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — Population Stratification + PCA Correction

**File**: `artifacts/lecture-13/05-pca-correction.html`
**Lecture anchor**: §4.1 Population stratification; §4.2 PCA correction
**EE framing reinforced**: PCA correction as noise whitening; top PCs capture the low-rank ancestry covariance.

### Teaching purpose

Simulate a two-subpopulation cohort with an ancestry-correlated phenotype and zero true SNP effects. Run GWAS without correction → massive genome-wide inflation. Include top $k$ PCs as covariates → $\lambda_{GC}$ returns to 1. Student controls the admixture fraction and phenotype-ancestry correlation.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Subpop A fraction:        [──●── 0.60 ]                     │
│ Subpop B fraction:        [──●── 0.40 ]                     │
│ F_ST (between-pop diff):  [──●── 0.10 ]                     │
│ Phenotype-ancestry corr: [──●── 0.30 ]                     │
│ N individuals:            [──●── 10,000 ]                   │
│ M SNPs:                   [──●── 20,000 ]                   │
│ # PCs to include:         [──●── 0 → 20 ]                   │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ PCA scatter (top):                                          │
│   PC1 vs PC2 of the genotype matrix                         │
│   coloured by subpopulation                                 │
├─────────────────────────────────────────────────────────────┤
│ Manhattan + QQ (middle, side-by-side):                     │
│   without PC correction                                     │
│   with top k PCs as covariates                              │
├─────────────────────────────────────────────────────────────┤
│ λ_GC readout:                                              │
│   no correction: 1.42                                       │
│   with PCs:      1.03                                       │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Zero causal SNPs planted.  Without correction, λ_GC    │   │
│ │  = 1.42 and several SNPs exceed genome-wide threshold. │   │
│ │ Including 2 PCs → λ_GC = 1.03; spurious hits vanish.   │   │
│ │ Turn correlation to 0 → no inflation even without PCs. │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Subpop A fraction: 0.1 – 0.9, default 0.6.
- $F_{ST}$ between subpops: 0 – 0.3, default 0.1.
- Phenotype-ancestry correlation: 0 – 0.6, default 0.3.
- $N$: 500 – 50,000 (log).
- $M$: 1,000 – 50,000 (log).
- PCs to include: 0 – 20, default 0 (no correction).
- Re-simulate.

### What they see

- PCA scatter with subpopulation colouring.
- Side-by-side Manhattan + QQ with vs without PC correction.
- $\lambda_{GC}$ readout.
- Outcome banner.

### Target aha

Start with default → $\lambda_{GC}$ large, spurious hits everywhere. Slide PCs to 2 → $\lambda_{GC} \to 1$, hits disappear. Turn phenotype-ancestry correlation to 0 → no inflation regardless of correction. Student sees: **the confounding was real, and PCA correction resolved it by projecting out the ancestry direction.**

### Technical notes

- Pure JS, seeded PRNG.
- Generate per-SNP allele frequencies per subpopulation (Balding-Nichols model with $F_{ST}$).
- Generate genotypes.
- Phenotype: $Y = \rho \cdot \text{ancestry-indicator} + \varepsilon$ (no SNP effects; all signals are ancestry artefacts).
- Compute top PCs via simplified randomised SVD.
- Run GWAS with $k$ PCs as additional covariates in each regression.
- $\lambda_{GC}$ computed for each version.

### Acceptance criteria

- [ ] Default → no-correction $\lambda_{GC}$ ≥ 1.3; with 2 PCs $\lambda_{GC}$ ≤ 1.1.
- [ ] Phenotype-correlation = 0 → $\lambda_{GC}$ ≤ 1.05 even without PCs.
- [ ] Higher $F_{ST}$ → larger inflation without correction.
- [ ] PC scatter correctly shows subpopulation structure.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — Burden and SKAT Simulator

**File**: `artifacts/lecture-13/06-burden-skat.html`
**Lecture anchor**: §7.3 SKAT and SKAT-O
**EE framing reinforced**: burden = sum statistic (matched filter on count); SKAT = kernel statistic (quadratic-form variance detector); SKAT-O = adaptive combiner.

### Teaching purpose

Simulate a gene containing a user-specified mix of rare functional variants. Student sets the directional agreement pattern. Compute burden, SKAT, and SKAT-O statistics and their p-values. Demonstrate when each test wins.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Rare variants in gene:    [──●── 20 ]                       │
│ Fraction deleterious:     [──●── 1.0 ]                      │
│ Fraction protective:       [──●── 0.0 ]                      │
│ Fraction neutral:          [──●── 0.0 ]                      │
│ Average effect size:       [──●── 0.5 SD ]                  │
│ Mean MAF:                  [──●── 0.003 ]                   │
│ Sample size N:             [──●── 20,000 ]                  │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Gene cartoon:                                               │
│   horizontal bar with variant tick marks coloured by effect │
├─────────────────────────────────────────────────────────────┤
│ Test statistics panel:                                      │
│   Burden:   T_burden = 4.2   p = 2e-5                       │
│   SKAT:     Q_SKAT   = 30.5  p = 1e-7                       │
│   SKAT-O:                      p = 3e-7                     │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Current pattern: all deleterious. All three tests fire.│   │
│ │ Flip half to protective → burden p ≈ 0.7 (cancels);    │   │
│ │   SKAT still p ≈ 1e-5 (detects variance contribution). │   │
│ │ SKAT-O ≈ min(burden, SKAT) — adaptive winner.          │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Rare variants in gene: 5 – 50, default 20.
- Fractions (deleterious / protective / neutral) summing to 1, with sliders enforcing the constraint.
- Average effect size: 0.1 – 1.0 SD, default 0.5.
- Mean MAF: 0.0001 – 0.01, default 0.003.
- $N$: 1,000 – 200,000 (log).
- Re-simulate.

### What they see

- Gene cartoon with coloured variant markers.
- Test statistics panel: burden, SKAT, SKAT-O with p-values.
- Outcome banner comparing patterns.

### Target aha

All deleterious → burden wins (smallest p). Half-half directional → burden cancels (p ≈ 1), SKAT still fires (small p). SKAT-O adaptively picks the winner. Student sees: the **right test depends on the underlying directional pattern**, which is unknown a priori — SKAT-O hedges.

### Technical notes

- Pure JS, seeded PRNG.
- Simulate genotypes: Binomial(2, MAF_j) per individual per variant.
- Phenotype: $Y_i = \sum_j \beta_j G_{ij} + \varepsilon_i$ where $\beta_j$ has the direction chosen by the mixture.
- Burden statistic: regress $Y$ on $\sum_j G_{ij}$; Wald test.
- SKAT statistic: $Q = (Y - \hat{\mu})^T G W G^T (Y - \hat{\mu})$ with $W$ MAF-inverse-weighted.
- SKAT p-value: Davies approximation (scaled chi-squared mixture).
- SKAT-O: compute for a grid of mixing parameters and take the minimum p-value after perturbation-based correction.

### Acceptance criteria

- [ ] Default all-deleterious → burden, SKAT, and SKAT-O all fire.
- [ ] Mixed directions → burden fails (p > 0.05); SKAT still fires.
- [ ] Neutral-only gene → all three fail.
- [ ] SKAT-O p is approximately the minimum of burden and SKAT p-values.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Artifact #7 — Kinship-Corrected GWAS Mini-Simulator

**File**: `artifacts/lecture-13/07-mixed-model.html`
**Lecture anchor**: §4.3 Cryptic relatedness and kinship; §4.4 Mixed-model regression
**EE framing reinforced**: mixed models as GLS with structured covariance; kinship as the variance-component weight matrix.

### Teaching purpose

Simulate a small cohort with user-controlled relatedness structure (randomly placed sibling groups, cousin pairs). Run both OLS GWAS and a mini-mixed-model GWAS. Compare QQ plots and $\lambda_{GC}$. Demonstrate that OLS inflates under cryptic relatedness and the mixed model corrects it.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ N individuals:          [──●── 1000 ]                       │
│ # sibling groups:       [──●── 10 ]  (each 3–5 sibs)        │
│ # cousin pairs:         [──●── 20 ]                         │
│ # unrelated singletons: [──●── remainder ]                  │
│ M SNPs (null):          [──●── 5,000 ]                      │
│ σ_g² (genetic variance):[──●── 0.3 ]                        │
│ σ_e² (residual):         [──●── 0.7 ]                        │
│ [Re-simulate]  [Run OLS + LMM]                              │
├─────────────────────────────────────────────────────────────┤
│ Kinship matrix heatmap:                                     │
│   N × N matrix; blocks visible for sibling groups           │
├─────────────────────────────────────────────────────────────┤
│ Side-by-side QQ plots:                                      │
│   OLS (left):  λ_GC shown                                   │
│   LMM (right): λ_GC shown                                   │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ OLS: λ_GC = 1.18 (inflated — cryptic relatedness)      │   │
│ │ LMM: λ_GC = 1.03 (corrected by kinship covariance)     │   │
│ │ Halve sibling groups → OLS λ_GC drops to 1.08           │   │
│ │ Raise σ_g² → more polygenic signal; both still match   │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- $N$: 200 – 5,000, default 1000.
- Number of sibling groups: 0 – 50, default 10.
- Number of cousin pairs: 0 – 100, default 20.
- $M$: 500 – 20,000, default 5000 (null SNPs — no true effects in this simulator).
- $\sigma_g^2$: 0 – 1, default 0.3.
- $\sigma_e^2$: 0 – 1, default 0.7.
- Re-simulate / Run.

### What they see

- Kinship matrix heatmap with sibling blocks visible.
- Side-by-side QQ plots (OLS vs LMM) with $\lambda_{GC}$ annotated.
- Outcome banner.

### Target aha

More relatives → larger OLS inflation. Apply the kinship-aware LMM → $\lambda_{GC}$ returns to ~1. Student sees: the **structured noise covariance** (relatedness) has to go into the model, not be ignored.

### Technical notes

- Pure JS, seeded PRNG.
- Build kinship matrix: unrelated pairs $K_{ij} = 0$; full siblings $K_{ij} = 0.5$; cousins $K_{ij} = 0.125$; self $K_{ii} = 0.5$ (1.0 = diploid self-kinship in strict convention).
- Simulate phenotype: $Y = \mathbf{u} + \boldsymbol{\varepsilon}$ with $\mathbf{u} \sim \mathcal{N}(0, \sigma_g^2 \mathbf{K})$, $\boldsymbol{\varepsilon} \sim \mathcal{N}(0, \sigma_e^2 \mathbf{I})$. No SNP effects.
- Simulate genotypes independently of phenotype.
- OLS: per-SNP linear regression.
- Mini-LMM: precompute $\mathbf{V}^{-1/2}$ via eigendecomposition of $\mathbf{K}$; transform $Y$ and $G_j$ into "whitened" space; OLS there.
- $\lambda_{GC}$ per approach.

### Acceptance criteria

- [ ] Default → OLS $\lambda_{GC}$ > 1.15; LMM $\lambda_{GC}$ < 1.08.
- [ ] Zero siblings / cousins → $\lambda_{GC}$ matches 1 in both.
- [ ] Increasing $\sigma_g^2$ → OLS inflation grows; LMM still matches 1.
- [ ] Kinship heatmap shows sibling blocks in the correct locations.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 9. Cross-Artifact Consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Manhattan / QQ / ancestry palettes identical across artifacts (§1.2 above).
- Log-scale sliders (for $N$, $M$, MAF) behave consistently.
- Every artifact emits an **outcome banner** per convention §1.4.
- KaTeX loaded (when used) from the same CDN script tag across artifacts.

## 10. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls beyond KaTeX (where used).
- [ ] `<script src="../_shared/resize.js" defer></script>` embedded near `</body>`.
- [ ] Outcome banner or equivalent verdict line visible at end of any user interaction.
