#import "../theme/book-theme.typ": *

= #idx("GWAS")GWAS and Statistical Genetics: A Million Parallel Detectors <ch:gwas>

#matters[
  By 2023, the published GWAS catalogue contained more than half a
  million reported variant-trait associations across roughly five
  thousand human traits. Every one of those associations is, mechanically,
  a regression coefficient at a single #idx("SNP")SNP in some cohort, surviving a
  p-value threshold of $5 times 10^(-8)$. The shape of that pipeline
  is so uniform that you can read most GWAS papers in fifteen minutes
  once you know what to look at: how many cases, how many controls,
  how many SNPs imputed, the #idx("Manhattan plot")Manhattan plot, the #idx("QQ plot")QQ plot, the inflation
  factor, the lead-SNP table. What the field has actually been doing
  for twenty years is running a very large multi-hypothesis detection
  problem against a noisy, structured background. Almost every
  controversy in human genetics over that period — missing heritability,
  the polygenic score portability crisis, the long argument about
  whether GWAS findings are mostly tagging or mostly causal — is a
  consequence of mis-specifying one of the boxes in that pipeline.
  This chapter unpacks each box.
]

A genome-wide association study answers one question: which positions
in the genome carry an allele whose count correlates with a phenotype
across a population? The question sounds modest, and the model used
to answer it at each position is one a first-year statistics student
already knows — a linear or logistic regression. What makes GWAS
different is not the per-position model. It is the fact that the
model runs in parallel at $10^6$ effectively-independent positions on
the same cohort, that the regressors are correlated by linkage
disequilibrium, that the residuals are correlated by relatedness and
ancestry, and that the effects being detected are small enough that
the whole exercise lives or dies on the calibration of its tail
probabilities.

This chapter walks through the pipeline in seven sections.
@sec:framework lays out the per-SNP regression and the experimental
designs (case-control and quantitative trait) that feed it.
@sec:manhattan covers the two plots every GWAS paper publishes — the
Manhattan and the QQ — and what they tell a trained reader.
@sec:gwas-multiple-testing derives the famous $5 times 10^(-8)$ threshold
and explains why it is the right number despite the obvious
arithmetic mismatch between $5 times 10^(-8)$ and the ten million
SNPs in a modern imputation panel. @sec:confounders is the longest
section, because confounding is what kills naive GWAS: population
stratification, cryptic relatedness, mixed-model regression, batch
effects. @sec:heritability turns from per-SNP detection to
genome-wide summaries — heritability estimation, LD-score
regression, polygenic risk scores, and the portability problem that
has become the defining ethical crisis of the field.
@sec:finemapping treats #idx("fine-mapping")fine-mapping as the sparse inverse problem
it is, and shows how SuSiE-style credible sets fall out.
@sec:rare-variants closes with the rare-variant regime, where
single-SNP testing collapses and burden and #idx("SKAT")SKAT tests take over.

The chapter assumes you have read Chapter 12 on population genetics
— in particular the parts on #idx("allele frequency")allele frequency, #idx("linkage disequilibrium")linkage disequilibrium,
and Hardy–Weinberg equilibrium, all of which become load-bearing
inside the next ten pages.


== The GWAS Framework <sec:framework>

The motivating question is direct: pick a trait — height, LDL
cholesterol, schizophrenia susceptibility, age at menarche, COVID-19
hospitalization — and find the genetic variants whose allele counts
correlate with it across a cohort.

Two experimental designs dominate. *Case-control* studies sample
$N_1$ individuals affected by a disease and $N_0$ unaffected controls,
genotype both groups at $M$ SNPs, and test for allele-frequency
differences at each SNP. *Quantitative-trait* studies sample $N$
individuals from a population, measure a continuous phenotype on each,
and test whether allele count at each SNP predicts the measurement.
The two designs run different per-SNP statistical tests but produce
the same output object: a table with one row per SNP and at least three
columns — an effect-size estimate, a standard error, and a p-value.

#figure(
  image("../../diagrams/lecture-13/01-cohort-design.svg", width: 95%),
  caption: [
    The two canonical GWAS designs. Case-control compares allele
    frequencies between affected and unaffected groups; quantitative-trait
    regresses a continuous measurement on genotype. The output of either
    is one effect estimate, one standard error, and one p-value per SNP.
  ],
) <fig:cohort>

=== Per-SNP Regression

At a single SNP, the additive model for a quantitative trait $Y$ on
individual $i$ is

$ Y_i = beta_0 + beta dot G_i + bold(gamma)^T bold(X)_i + epsilon_i $

where $G_i in {0, 1, 2}$ counts copies of the "effect" allele, $bold(X)_i$
is a vector of covariates (age, sex, principal components — see
@sec:confounders), and $epsilon_i$ is residual noise assumed normal
with variance $sigma^2$. The null hypothesis at this SNP is
$H_0 : beta = 0$; the alternative is $H_1 : beta != 0$. The Wald
statistic is $z = hat(beta) / "SE"(hat(beta))$, and its two-sided p-value
is computed against a standard normal.

For a binary trait, swap linear regression for logistic:

$ log P(Y_i = 1) / P(Y_i = 0) = beta_0 + beta dot G_i + bold(gamma)^T bold(X)_i $

The likelihood-ratio or score test against $H_0 : beta = 0$ gives a
p-value the same way.

Three genetic models exist. Under the *additive* model, $G$ is treated
as a numeric predictor and each additional copy of the effect allele
contributes another $beta$ to the phenotype. Under the *dominant*
model, $G$ is collapsed to a binary indicator of carrying at least one
effect allele, which is right when one copy is sufficient to drive the
trait. Under the *recessive* model, $G$ is a binary indicator of
carrying two effect alleles, right when both copies matter. The
additive model is the workhorse — virtually all published GWAS use it
as the default and only run dominant or recessive as sensitivity
analyses on specific variants.

#figure(
  image("../../diagrams/lecture-13/02-single-snp-regression.svg", width: 95%),
  caption: [
    The additive regression model at a single SNP. Each additional copy
    of the effect allele shifts the phenotype mean by $beta$ — the slope
    of the regression line through the three genotype means. The inset
    shows the sampling distribution of $hat(beta)$ under the null and
    alternative.
  ],
) <fig:single-snp>

=== The Scale

GWAS data come in three pieces: a genotype matrix of $N$ individuals
by $M$ SNPs (entries in ${0, 1, 2}$ for direct genotypes, or
real-valued dosages in $[0, 2]$ for imputed SNPs), a phenotype vector
of length $N$, and a covariate matrix that is usually 10–30 columns
wide. Modern genotyping arrays type 500 thousand to a few million SNPs
directly; imputation against a reference panel — TOPMed, the HRC, or
the more recent #idx("pangenome")pangenome — expands this to roughly 10 million SNPs
in the working table, with the caveat that imputation quality (the
INFO score) varies and low-INFO sites get dropped before testing.

The regression then runs $M$ times — one fit per SNP, with the per-SNP
$G$ vector and the shared covariate matrix as inputs. The output is
the *GWAS summary statistic table*: $M$ rows, each carrying at least
a SNP identifier, #idx("chromosome")chromosome, position, reference and effect alleles,
allele frequency, effect-size estimate $hat(beta)$, standard error
$"SE"(hat(beta))$, and p-value. This summary table — not the underlying
genotype matrix — is what gets shared between research groups and is
the input to nearly every downstream analysis (#idx("PRS")PRS, #idx("LDSC")LDSC, fine-mapping,
#idx("colocalization")colocalization). Treating the summary statistics as the data product
of GWAS, rather than the raw genotypes, is most of what makes the
modern statistical-genetics ecosystem possible.

#note[
  A GWAS is a bank of $M approx 10^6$ parallel single-channel detectors,
  one per SNP, each running a classical likelihood-ratio or #idx("Wald test")Wald test
  and reporting a p-value. The genome is the frequency axis of a
  coarse filter bank; each SNP is one filter; the phenotype is the
  signal being correlated against each channel. You then have to
  control the false-alarm rate across all $M$ channels — a
  multi-hypothesis detection problem with the same structure as
  multi-band radar, multi-target tracking, or multi-pixel anomaly
  detection in imaging. The two defining engineering problems —
  controlling genome-wide false-alarm rate (@sec:gwas-multiple-testing),
  and handling correlated channels via LD (@sec:manhattan,
  @sec:finemapping) — are familiar problems in disguise.
]

=== A Short Historical Arc

The first GWAS to genuinely catch the field's #idx("attention")attention was Klein et al.
(2005), testing 116,204 SNPs in 96 cases of age-related macular
degeneration and 50 controls — laughably small by today's standards.
A common variant in the complement-factor-H gene (`CFH`) came out at
$p approx 10^(-7)$ and replicated. The result was strong enough to
make the case for hypothesis-free, genome-wide testing of common
variants in common disease.

The big bang of the field, two years later, was the Wellcome Trust
Case Control Consortium's 2007 paper — 500,000 SNPs typed across
14,000 cases of seven common diseases (bipolar, coronary artery
disease, Crohn's, hypertension, rheumatoid arthritis, type 1 and
type 2 diabetes) plus 3,000 shared controls. The WTCCC paper
established almost every methodological convention used today: the
$5 times 10^(-8)$ threshold, the Manhattan plot, the QQ plot
inflation diagnostic, principal-component covariates for stratification
control, and the practice of designating SNPs as "lead" within an
LD-defined locus. Two of the seven diseases — Crohn's and type 1
diabetes — yielded multiple genome-wide significant hits;
the others required larger meta-analyses to catch up.

The cohort sizes have been racing upward ever since. The big modern
biobanks — UK Biobank ($approx 500,000$ participants with array
genotypes, exomes, and now whole genomes), FinnGen ($approx 500,000$
Finns with founder-population enrichment), Biobank Japan, China
Kadoorie, All of Us, the US Million Veteran Program — provide the
sample sizes needed to detect the small effects that dominate common
disease. By 2020, a height GWAS at $N approx 5,000,000$ recovered
over twelve thousand independent associated SNPs; coronary disease at
$N approx 1,300,000$ recovered over two hundred. The hits at the top
of the Manhattan plot are by now an established biological catalogue;
the disagreements have moved downstream, to fine-mapping and to
mechanistic interpretation.


== The Manhattan Plot and the QQ Plot <sec:manhattan>

A GWAS paper has two diagnostic plots, and a trained reader looks at
them in order. The *Manhattan plot* shows the signal. The
*QQ plot* shows the calibration. A paper that publishes only the
Manhattan is hiding something; a paper that publishes only the QQ
plot has no findings to report.

=== Manhattan Plots

A Manhattan plot puts one dot per SNP, with genomic position on the
x-axis (concatenating chromosomes 1 through 22, X, and sometimes Y
left to right) and $-log_(10)(p)$ on the y-axis. Alternating
chromosome tints make the boundaries visible. A horizontal dashed
line at $y = 7.3$ marks the genome-wide significance threshold of
$5 times 10^(-8)$. Peaks rising above the line are the "hits."

#figure(
  image("../../diagrams/lecture-13/03-manhattan-plot.svg", width: 95%),
  caption: [
    A Manhattan plot from a well-powered GWAS. Each dot is one SNP;
    peaks above the dashed threshold at $-log_(10)(p) = 7.3$ are
    genome-wide significant. The cluster of elevated SNPs around each
    lead is LD-driven peak structure — neighbours of a causal variant
    inherit some of its signal through correlated allele counts.
  ],
) <fig:manhattan>

Four things a reader extracts from the plot:

- The number of peaks above $5 times 10^(-8)$ — the count of
  independent hits, after collapsing within-locus correlation.
- The peak heights — strength of association. Very tall peaks
  ($-log_(10)(p) > 30$) mean either very large effects or very large
  $N$; small-effect common-disease GWAS rarely produces them outside
  of mixture loci like `HLA` or `FTO`.
- The peak shapes — a tall column with neighbouring SNPs also
  elevated signals LD-driven structure around a causal variant. An
  isolated single tall SNP with no elevated neighbours is much
  more suspicious; it suggests a genotyping artifact or a private
  variant rather than a clean signal.
- The chromosome distribution — clustering of hits on specific
  chromosomes for polygenic traits is expected and reflects
  chromosome length and gene density, not biology per se.

=== QQ Plots and the Inflation Factor

The QQ plot compares the observed quantiles of $-log_(10)(p)$ across
all tested SNPs against what would be expected under the global null
(no true associations, no confounding). The expected distribution is
uniform on $[0,1]$ for the p-values, so $-log_(10)(p)$ follows an
exponential distribution with a known set of quantiles. Plotting
observed against expected produces a straight diagonal under the
null. Two kinds of deviation matter.

*Late departure* — the upper-right tail of the plot lifts above the
diagonal while the lower-left stays on it. This is what a successful
GWAS looks like. The strongest hits are real signals exceeding what
you would see under pure null; everything else is consistent with
calibrated noise.

*Early departure* — even moderately-small p-values are above the
diagonal, with the lift starting near the origin. This is *inflation*:
the whole distribution of test statistics is shifted up. The classic
cause is population stratification, but cryptic relatedness,
genotyping batch effects, miscoded covariates, and unmodelled
phenotype heterogeneity all produce the same shape.

The standard one-number summary is the *genomic control inflation factor*

$ lambda_("GC") = "median"(chi^2_"observed") / "median"(chi^2_"null") = "median"(chi^2_"observed") / 0.456 $

where $0.456$ is the median of a $chi^2_1$ distribution. $lambda_("GC") = 1$
means no inflation; values up to about $1.05$ are routinely tolerated;
$lambda_("GC") > 1.1$ raises eyebrows; $lambda_("GC") > 1.2$ is a red
flag. The complication: for highly polygenic traits, $lambda_("GC")$
also goes up because of genuine widespread signal, not just confounding.
A height GWAS at $N approx 700,000$ has $lambda_("GC") approx 2$ from
real polygenic signal alone. The modern fix — LD-score regression
(@sec:heritability) — partitions the inflation into a true
polygenic-signal component and a confounding component.

#figure(
  image("../../diagrams/lecture-13/04-qq-plot.svg", width: 95%),
  caption: [
    The QQ plot is the calibration check for a GWAS. A well-controlled
    study (left) follows the diagonal and lifts only at the top-right
    where real hits live. An inflated study (right) lifts early,
    indicating widespread test-statistic elevation from confounders or
    polygenic signal.
  ],
) <fig:qq>

=== Peak Structure and LD

A causal SNP with a true effect rarely shows up alone on the Manhattan
plot. It sits in a *linkage disequilibrium* block — a region of the
genome where genotype correlations between nearby SNPs are high
because too few historical #idx("recombination")recombination events have separated them.
Every SNP in LD with the causal one carries a fraction of the signal:
its allele frequency in cases vs controls tracks the causal SNP's
allele frequency in cases vs controls, scaled by $r^2$ between the two.

A typical European-ancestry LD block spans 10 to 100 kilobases. Inside
the block, a causal variant with a clean $-log_(10)(p) = 25$ produces
neighbours with $-log_(10)(p)$ scaled roughly by $r^2$ — the SNP at
$r^2 = 0.8$ shows $-log_(10)(p) approx 20$, the one at $r^2 = 0.4$
shows $-log_(10)(p) approx 10$, and so on, fading into the background
at the edges of the block. The visible "peak" is the LD shadow of one
or a few causal variants.

This is simultaneously a blessing and a curse. The blessing: a peak
shape — a tall lead with elevated neighbours falling off smoothly on
both sides — is harder to fake than an isolated single SNP. A
genotyping artifact tends to be a one-SNP spike; a real signal is a
hill. The curse: you cannot tell from the Manhattan plot which SNP in
the peak is causal. The "lead SNP" (highest $-log_(10)(p)$) is often
not the causal variant — it is the best-imputed tag for the causal
variant in your reference panel. Going from peak to causal variant
is a separate analysis (@sec:finemapping).

#note[
  The Manhattan plot is to GWAS what a spectrogram is to an audio
  signal. Lay the genome along an axis, compute a test statistic at
  each position, plot the strength of the response. Peaks above the
  noise floor are detections; the clutter below is the null. LD-driven
  peak width is filter bandwidth — each causal variant's signal leaks
  into correlated neighbours, broadening its #idx("footprint")footprint on the
  spectrum. Wider LD blocks mean broader peaks and fewer effective
  independent tests, which feeds directly into the multiple-testing
  arithmetic of the next section.
]


== #idx("multiple testing")Multiple Testing at SNP Scale <sec:gwas-multiple-testing>

If you run $M$ independent tests at per-test significance level
$alpha$, the expected number of false positives is $M alpha$. To
control the family-wise error rate at $alpha = 0.05$ across all $M$
tests, #idx("Bonferroni")Bonferroni demands the adjusted per-test level
$alpha_"adj" = alpha / M$. For $M = 10^6$ this gives
$alpha_"adj" = 5 times 10^(-8)$ — the genome-wide significance
threshold that has defined GWAS for almost twenty years.

The arithmetic is simple. The question that took years to settle is
where the $10^6$ comes from when modern imputed datasets contain
roughly $10^7$ SNPs.

=== Where $10^6$ Comes From

The answer is that the $10^7$ SNPs in an imputed dataset are not
independent. Two SNPs in tight LD ($r^2 > 0.8$) give correlated test
statistics — if one is significant, the other is very likely
significant too, simply because their genotypes are nearly the same.
Bonferroni-correcting at $M = 10^7$ would massively over-correct,
because most of those tests are not genuine extra opportunities to
make a false discovery.

The right count is the *effective number of independent tests* —
the number of pairwise-near-independent loci across the genome.
Several papers estimated it empirically in the late 2000s using
LD-pruning of HapMap and 1000 Genomes data; the consensus number for
European-ancestry cohorts came out to roughly $10^6$. African-ancestry
cohorts have shorter LD blocks and therefore more effectively
independent tests — about $1.5 times 10^6$ — and a stricter threshold
should in principle be used, although in practice most published GWAS
of African-ancestry cohorts still use $5 times 10^(-8)$.

#figure(
  image("../../diagrams/lecture-13/05-effective-tests.svg", width: 95%),
  caption: [
    The Bonferroni threshold of $5 times 10^(-8)$ comes from the
    effective LD-independent SNP count, not the raw imputed SNP
    count. The roughly $10^7$ SNPs in a modern panel collapse to
    roughly $10^6$ LD-independent tests; $0.05 / 10^6 = 5 times 10^(-8)$.
  ],
) <fig:effective-tests>

#figure(
  image("../figures/ch13/f1-bonferroni-arithmetic.svg", width: 95%),
  caption: [
    The standard $5 times 10^(-8)$ threshold survives roughly across
    array densities because the LD-independent test count is bounded
    above by the genome's natural LD-block structure, not by how many
    SNPs the array happens to type. Denser arrays add correlated tests
    inside the same blocks, which add power per locus but not
    additional independent comparisons.
  ],
) <fig:bonferroni>

=== FWER versus #idx("FDR")FDR

Two different error-control philosophies live in statistical genetics
and adjacent fields.

*Family-wise error rate (FWER).* The probability of any false
positive across all tests. Controlled by Bonferroni or its refinements
(Holm, Hochberg). Strict and conservative. The GWAS standard.

*#idx("false discovery rate")False discovery rate (FDR).* The expected fraction of declared
discoveries that are false. Controlled by the Benjamini–Hochberg
procedure. Lets more real hits through at the cost of a known
contamination rate.

GWAS uses FWER. The reason is practical: each genome-wide significant
hit becomes a biological hypothesis that some downstream team will
spend years working on — fine-mapping, eQTL colocalization, mouse
knockouts, drug-target validation. A 5 % FDR across two hundred hits
would mean ten of them are bogus, distributed among the rest at
random. The cost of chasing the false ones is too high for the field
to tolerate. eQTL and other molecular-trait studies, where downstream
follow-up is cheaper and the per-hit cost of being wrong is small,
sometimes use FDR instead.

=== Power and Why Cohorts Keep Growing

The power of a single-SNP test depends on four numbers: the effect
size $beta$, the allele frequency $p$, the sample size $N$, and the
significance threshold $alpha$. The non-centrality parameter of the
chi-squared statistic under the additive model is approximately

$ chi^2 approx N dot 2 p (1 - p) dot beta^2 / sigma^2 $

To detect $beta = 0.02$ standard deviations at allele frequency
$p = 0.2$ with $alpha = 5 times 10^(-8)$ and 80 % power, this comes
out to $N approx 230,000$. That is the order of magnitude that makes
biobanks necessary. For smaller effects — say $beta = 0.01$, common
for individual loci in a polygenic trait — the required $N$ scales as
$1 / beta^2$, so quadruples to roughly $900,000$.

The arithmetic explains the field's two-decade #idx("trajectory")trajectory. The 2007
WTCCC paper at $N = 17,000$ could detect $beta approx 0.1$ at common
frequency and 80 % power; today's $N = 500,000$–$5,000,000$ studies
detect down to $beta approx 0.01$ and lower. Most of the gain has gone
into picking up smaller and smaller effect sizes for the same
biologically interesting traits; the per-locus effects have not
gotten bigger, the detection floor has dropped.

#figure(
  image("../figures/ch13/f2-power-curves.svg", width: 95%),
  caption: [
    Sample-size requirements at $alpha = 5 times 10^(-8)$ and 80 %
    power. Each curve fixes a target effect size $beta$ in trait
    standard deviations; the y-axis reports the cohort size needed to
    detect it at the given allele frequency. The biobank scale of
    $10^5$–$10^6$ samples covers most of the common-variant common-disease
    landscape; rare-variant detection at $p < 0.001$ remains out of
    reach for single-SNP tests at any feasible cohort size.
  ],
) <fig:power>

#warn[
  Genome-wide significant does not mean true. It means the
  observed test statistic is unlikely under the pure-null hypothesis
  of zero effect and no confounding. In practice, confounding from
  population stratification, batch effects, phenotype measurement
  error, and survivorship bias all produce genome-wide significant
  hits that are not causal. Replication in an independent cohort,
  PheWAS scanning for off-trait associations, and direct mechanistic
  follow-up are required before any hit gets called real.
  $5 times 10^(-8)$ is a necessary condition, not a sufficient one.
]

This is the fourth time in this book that multi-hypothesis correction
at scale has appeared. Chapter 4 used FDR on per-site variant calls,
chapter 6 on differential-expression genes, chapter 9 on #idx("ChIP-seq")ChIP-seq
peaks, and now $10^6$ tests across the genome. The mathematical
structure — null distribution, test statistic, correction for the
effective number of tests — recurs in every domain. What changes is
the definition of "effective number of independent tests": LD blocks
for GWAS, gene co-expression structure for #idx("differential expression")differential expression,
peak spacing for ChIP-seq, the #idx("coalescent")coalescent null for #idx("selection")selection scans.
Knowing where the effective count comes from in each setting is most
of what separates a clean p-value from an inflated one.


== Confounders and Their Corrections <sec:confounders>

Population structure is the most reliable way to wreck a GWAS. It
manifests as systematic allele-frequency differences between groups
that also differ in phenotype for reasons that have nothing to do
with genetic causation. Any SNP whose allele frequency tracks the
group structure will appear associated with the trait, even though
the causal driver is whatever made the groups behave differently.
Twenty years of methods development on GWAS is largely the story of
controlling for this.

=== Population Stratification

A canonical toy example: suppose a cohort is 60 % from Northern
Europe and 40 % from Southern Europe, and the trait of interest is
something cultural — a dietary preference, a behavioural outcome,
something with environmental rather than genetic drivers. SNPs
whose allele frequencies differ between the two source populations
(and many SNPs do, even within Europe) will correlate with the trait
through ancestry, not biology. The association is real in the
arithmetic sense — allele counts do correlate with phenotype — but
causally spurious. The actual driver is geography, not genetics.

Stratification is the *default state* of any non-randomized cohort.
The only way to eliminate it would be to draw participants
independently from a panmictic population, which has not existed in
human history. The practical task is therefore not to remove
stratification but to model it well enough that the per-SNP test
statistics are calibrated. Without correction, a stratified GWAS
inflates across the board — visible as QQ-plot early departure and
$lambda_("GC") > 1$ — and produces false positives concentrated at
ancestry-differentiated SNPs.

#figure(
  image("../../diagrams/lecture-13/06-pca-structure.svg", width: 95%),
  caption: [
    Population structure visible in #idx("PCA")PCA space. The top two principal
    components of the genotype matrix separate ancestry groups cleanly
    (left). When phenotype also tracks ancestry (right), any SNP
    with allele-frequency differences between groups will spuriously
    associate with the trait unless the structure is conditioned on.
  ],
) <fig:pca-structure>

=== PCA Correction

The standard first-line correction, introduced by Price et al. (2006),
is to compute principal components of the genotype matrix and include
the top $k$ PCs as covariates in every per-SNP regression. The
intuition is that the top PCs of a population's genotype matrix
capture exactly the kind of large-scale, between-group allele-frequency
variation that defines population structure. Conditioning the per-SNP
regression on those PCs removes the between-group signal while leaving
within-group genetic variation in place.

The recipe:

1. Center and scale the genotype matrix. Optionally prune to a set of
   common, approximately unlinked SNPs (one per LD block).
2. Compute the top $k$ eigenvectors of the centered genotype matrix's
   covariance, typically with a randomized SVD on the order of a few
   minutes for biobank scale.
3. For each individual, extract the loadings on these $k$ PCs.
4. Include them as columns of the covariate matrix $bold(X)$ in
   every per-SNP regression.

How many PCs? Practice converges on 10 for homogeneous cohorts, 20 for
heterogeneous ones, and 40 for the most diverse — UK Biobank uses 40
for analyses spanning the full participant pool. Too few PCs leave
residual stratification visible in the QQ plot; too many eat
statistical power without commensurate benefit (each PC costs one
degree of freedom).

#figure(
  image("../figures/ch13/f3-pca-correction.svg", width: 95%),
  caption: [
    Population stratification and its correction. Without conditioning
    on the top genotype PCs, ancestry-differentiated SNPs across the
    genome appear associated with an ancestry-correlated phenotype,
    inflating the QQ plot uniformly. Adding the top ten PCs as
    covariates absorbs the ancestry component into the model's noise
    structure, restoring a calibrated QQ plot with $lambda_("GC") approx 1$.
  ],
) <fig:pca-correction>

#note[
  PCA correction is noise whitening. A GWAS with stratification is a
  detection problem where the per-SNP noise has a non-white covariance
  structure: test statistics at SNPs with similar allele-frequency
  differentiation patterns are correlated through ancestry. The top
  principal components of the genotype matrix are the principal
  eigen-directions of that noise covariance. Including them as
  covariates is mathematically equivalent to projecting them out of
  the residuals — the exact operation noise whitening performs in
  adaptive-array signal processing or matched-filter detection
  against coloured noise. The residual test statistic is then closer
  to white, restoring valid per-SNP inference.
]

=== Cryptic Relatedness and #idx("kinship")Kinship

A second confounder hides inside the cohort even after stratification
has been handled. In any nominally-unrelated sample of $N$ individuals
drawn from a real population, some pairs will turn out to be
biological relatives — cousins, half-siblings, even occasional
parent-child or sibling pairs missed at recruitment. Relatives share
alleles beyond random expectation, which means their phenotypes are
not statistically independent draws. Treating them as independent
inflates the test statistic at any SNP correlated with the kinship
structure.

The fix is to compute a *kinship matrix* $bold(K)$ (or genetic
relationship matrix, GRM) that summarizes pairwise genetic similarity
across all participants. The entry $K_(i j)$ is roughly twice the
probability that a randomly-chosen allele at a randomly-chosen SNP is
identical-by-descent between individuals $i$ and $j$. For unrelated
pairs $K_(i j) approx 0$; for second cousins, $approx 0.03$; for
half-siblings, $approx 0.25$; for parent-child, $approx 0.5$;
for monozygotic twins or self-pairs, $approx 1.0$.

#figure(
  image("../../diagrams/lecture-13/07-kinship-matrix.svg", width: 90%),
  caption: [
    The kinship matrix (genetic relationship matrix) for a cohort of
    100 individuals. The diagonal is uniformly high ($approx 1$ for
    self-pairs); family clusters appear as small off-diagonal blocks.
    Mixed-model tools incorporate this matrix as a structured
    covariance term to prevent inflation from cryptic relatedness.
  ],
) <fig:kinship>

=== Mixed-Model Regression

The modern workhorse for biobank-scale GWAS is *mixed-model regression*,
which folds the kinship matrix into the residual covariance and runs
the per-SNP regression as a generalized least-squares problem rather
than ordinary least squares. The full model for a quantitative trait is

$ bold(Y) = bold(X) bold(gamma) + bold(G) beta + bold(u) + bold(epsilon) $

where $bold(u) tilde cal(N)(0, sigma_g^2 bold(K))$ is a random effect
with *kinship-structured covariance* and $bold(epsilon) tilde
cal(N)(0, sigma_e^2 bold(I))$ is residual environmental noise.
Equivalently, the joint residual covariance is

$ bold(V) = sigma_g^2 bold(K) + sigma_e^2 bold(I) $

and the per-SNP regression is a GLS fit weighted by $bold(V)^(-1)$:

$ hat(beta) = (bold(G)^T bold(V)^(-1) bold(G))^(-1) bold(G)^T bold(V)^(-1) bold(Y) $

Solving this naively requires inverting an $N times N$ matrix at every
SNP, which is computationally hopeless at biobank scale.
Algorithmic shortcuts make it tractable:

- *#idx("BOLT-LMM")BOLT-LMM* (Loh et al., 2015). Iterative conjugate-gradient methods
  with spectral approximations to the kinship matrix. Scales to
  $N > 500,000$; the default choice for quantitative traits in UK
  Biobank.
- *#idx("SAIGE")SAIGE* (Zhou et al., 2018). Saddle-point approximation for binary
  traits with severe case-control imbalance — a regime where ordinary
  logistic regression's test statistic is poorly calibrated even
  before kinship enters the picture.
- *REGENIE* (Mbatchou et al., 2021). A two-stage approach: a first
  stage fits ridge regression in blocks to compute a per-individual
  polygenic prediction, then a second stage runs the per-SNP test
  conditional on that prediction. Scales to UK Biobank's full
  $500,000 times 10,000,000$ table on commodity hardware.

#figure(
  image("../../diagrams/lecture-13/08-lmm-vs-ols.svg", width: 95%),
  caption: [
    Mixed-model regression versus OLS on the same cohort with cryptic
    relatedness. OLS treats residuals as independent and inflates
    test statistics. The mixed model absorbs kinship into the
    covariance matrix $bold(V)$, restoring calibration — the
    noise-whitening interpretation.
  ],
) <fig:lmm>

#note[
  A mixed-model GWAS is generalized least squares with a specific
  structured covariance: the total residual covariance is a
  genetic-kinship component plus an environmental-noise component.
  The kinship matrix encodes who is related to whom. Ordinary least
  squares assumes the residual covariance is the identity, which is
  white noise; the test statistic is then wrong whenever residuals are
  correlated across individuals. GLS with the inverse-covariance
  weight is the optimal linear unbiased estimator under correlated
  residuals (Gauss–Markov). BOLT-LMM, SAIGE, and REGENIE are
  efficient approximations to this GLS for biobank-scale data.
]

=== Batch Effects and Technical Confounders

Beyond population structure and kinship, a long tail of technical
confounders waits for the inattentive analyst.

*Genotyping batch.* Different plates, machines, sample-prep
operators, or time periods produce systematically different genotype
calls. The standard fix is to include batch indicators as covariates,
and to drop SNPs whose call rate or #idx("Hardy-Weinberg")Hardy-Weinberg p-value varies
across batches.

*#idx("DNA")DNA source.* Saliva and blood occasionally produce different
genotype distributions at specific SNPs because of bacterial DNA
contamination or differential cell-type sampling. Worth checking;
rarely a dominant effect.

*Imputation quality.* Imputed SNPs have an INFO score (the
imputation-server's estimate of how well the SNP is tagged by
neighbouring directly-genotyped SNPs). Low-INFO SNPs ($"INFO" < 0.8$)
have noisy dosages and inflate the test if included naively.
Filtering on INFO before testing is standard.

*Phenotype heterogeneity in multi-centre studies.* Different
recruitment sites have different ascertainment criteria, measurement
protocols, and demographic compositions. Include centre as a
covariate.

#warn[
  PCA correction is not always sufficient. For recently-admixed
  populations — Hispanic / Latino, African American, South Asian
  diaspora cohorts — a handful of top PCs misses fine-scale local
  ancestry structure that varies across chromosomal segments within
  individuals. Mixed models with the full kinship matrix do better.
  For founder populations like Finns or Ashkenazi Jews, a small
  number of PCs plus mixed models is standard but the residual
  inflation can still be non-trivial because LD blocks are wider and
  variance is concentrated in fewer haplotypes. Always re-inspect the
  QQ plot and $lambda_("GC")$ after correction; if still inflated,
  the model is still mis-specified somewhere.
]


== Heritability and Polygenic Risk Scores <sec:heritability>

The per-SNP picture is local. Genome-wide summaries — what fraction
of the trait's variance is genetic at all, what fraction of that
genetic variance is being captured by the GWAS, and how to aggregate
genome-wide signal into a per-individual prediction — require
different tools.

=== Heritability

*Heritability* is the fraction of phenotypic variance explained by
genetic variance:

$ h^2 = sigma_g^2 / sigma_p^2 $

with $sigma_g^2$ the genetic variance and $sigma_p^2$ the total
phenotypic variance. Two flavours show up. *Broad-sense*
heritability $H^2$ includes additive, dominance, and epistatic genetic
effects. *Narrow-sense* heritability $h^2$ includes only additive
effects, which is what GWAS measures under the additive model.
Twin and family studies have estimated $h^2$ for hundreds of traits;
representative values: height $approx 0.8$, BMI $approx 0.4$–$0.6$,
schizophrenia $approx 0.7$, type 2 diabetes $approx 0.3$.

Two modern methods estimate heritability directly from GWAS data.
*#idx("GCTA")GCTA* (Yang et al., 2011) estimates $h^2$ from the kinship matrix
plus phenotype via restricted #idx("maximum likelihood")maximum likelihood (REML) on the same
mixed-model framework as the per-SNP tests. *LD-score regression*
(LDSC, Bulik-Sullivan et al., 2015) estimates $h^2$ from summary
statistics alone, without ever touching the raw genotype matrix.

=== LD-Score Regression

LDSC is the most influential summary-statistic method of the 2010s.
It exploits a single, elegant identity: under a polygenic signal model,
a SNP's expected $chi^2$ test statistic scales linearly with how many
other SNPs it is in LD with — its *LD score*, $ell_j = sum_(k) r^2_(j k)$,
summed over neighbouring SNPs $k$.

The model is

$ EE[chi^2_j] = 1 + N h^2 ell_j / M + N a $

with $N$ the GWAS sample size, $M$ the SNP count, and $a$ a
confounding offset. Regressing the per-SNP $chi^2$ statistic on
LD score across all SNPs gives a slope proportional to $h^2$ and an
intercept of $1 + N a$ — which separates polygenic signal from
confounding.

The diagnostic application is direct. A GWAS with $lambda_("GC") = 1.4$
could be inflated from polygenic signal or from confounding;
$lambda_("GC")$ alone cannot distinguish them. LDSC can: a slope much
larger than the intercept-minus-one says polygenic, a flat slope with
large intercept says confounding, and the typical real GWAS shows some
of each.

#note[
  LDSC is a linear regression of per-SNP test statistic on LD score.
  The LD score is the row sum of the squared LD matrix. The regression
  slope is a scalar summary of how polygenic signal spreads across
  the LD-eigenvalue spectrum of the genotype covariance — a spectral
  decomposition of the kinship-driven correlation structure. The
  intercept isolates a uniform offset across the spectrum (confounding)
  from the slope (the spectral footprint of true polygenic signal).
  Same machinery as separating coherent signal from white-noise floor
  in a spectrum analyzer.
]

=== Missing Heritability

Twin studies estimate $h^2 approx 0.8$ for height. The sum of effects
of GWAS-significant SNPs from a 2010 height GWAS explained about
$0.05$ — five percent of the trait's variance, against an expected
ceiling of eighty. The gap was the original *missing heritability*
problem, and it dominated methodological discussion in the field for
a decade.

By 2023 the gap had narrowed considerably. SNP-heritability from LDSC
on the largest height GWAS recovers $h^2_"SNP" approx 0.5$, against
twin $h^2 approx 0.8$ — most of the supposedly missing heritability
turns out to live in many SNPs with effects too small to clear
$5 times 10^(-8)$ individually. Polygenic risk scores (below) capture
this implicitly. The residual gap to twin estimates can be parcelled
into:

- *Effects too small to be individually significant*, but captured
  by a PRS or by LDSC. Most of the recovery.
- *Rare variants* missed by common-SNP arrays. Sequencing studies
  recover some additional fraction.
- *Structural variants* under-represented in SNP-array data and
  better captured by long reads (Chapter 11) and the pangenome.
- *Non-additive effects* (dominance, epistasis) — small but non-zero.
- *Gene–environment interactions* — hard to estimate, residually
  important for some traits.

The lesson: "missing" heritability was mostly an artifact of stopping
at the $5 times 10^(-8)$ threshold. The signal was there; the
threshold hid it.

=== Polygenic Risk Scores

A *#idx("polygenic risk score")polygenic risk score* (PRS) aggregates signal across many SNPs into
a single per-individual score, by taking a weighted sum of allele
counts with weights given by GWAS effect estimates:

$ "PRS"_i = sum_(j in S) hat(beta)_j G_(i j) $

where the sum is over some subset $S$ of SNPs. The choice of subset
distinguishes the PRS methods:

- *Clumping + thresholding* (C+T, PRSice). Include SNPs at
  $p < alpha$, after LD-pruning to keep only one SNP per locus. The
  $alpha$ threshold is a hyperparameter; the optimal value depends on
  the trait and is usually found by cross-validation against an
  independent target cohort.
- *LDpred / LDpred2*. Bayesian #idx("shrinkage")shrinkage with a spike-and-slab prior
  on the per-SNP effects, applied genome-wide. Outperforms C+T because
  it borrows strength across all SNPs rather than just the significant
  ones.
- *PRS-CS*. Continuous shrinkage prior (a horseshoe-like distribution)
  fit to summary statistics with reference LD; typically the
  best-performing summary-statistic method as of 2023.
- *SBayesR / SBayesS*. Bayesian regression with biological priors on
  effect distributions across allele-frequency strata. Comparable
  performance to PRS-CS on most traits.

PRS has become genuinely clinically useful for a handful of traits.
For coronary artery disease, individuals in the top one percent of
the PRS distribution have about three times the population risk —
comparable to inheriting a rare Mendelian variant in `LDLR` or `PCSK9`.
For breast cancer, a well-calibrated PRS identifies a top-tier
roughly 20 % of women whose lifetime risk approaches that of
`BRCA1`/`BRCA2` carriers without rare-variant disease in either gene.
For schizophrenia and bipolar disorder, PRS captures a meaningful
fraction of within-family variance and is being explored as a
risk-stratification tool in clinical settings.

=== The Portability Problem

A PRS trained on European-ancestry GWAS performs roughly three to
five times worse — measured as $R^2$ of PRS against phenotype — when
applied to non-European individuals. The drop is steepest for
African-ancestry cohorts and intermediate for East Asian and South
Asian cohorts. Three factors drive it.

*Allele-frequency differences.* SNPs common in Europeans can be rare
or absent in Africans (and vice versa). A SNP that contributes
$hat(beta) dot G$ to the PRS in Europeans where $G$ varies across
$[0, 2]$ contributes $hat(beta) dot 0 = 0$ everywhere in a population
where the allele has been lost — its slot in the PRS is wasted.

*LD-pattern differences.* The "tag SNPs" that carry the GWAS signal
in the training cohort are not always tagging the same causal
variant in a different ancestry. African-ancestry LD blocks are
shorter, so a European-trained tag SNP at $r^2 = 0.9$ with the
causal variant may sit at $r^2 = 0.3$ with the same causal variant
in an African cohort, slashing its contribution to the per-individual
score.

*Effect-size heterogeneity.* Some genuine causal variants have
different effect sizes across ancestries because of differing
gene-environment interactions or differing genetic backgrounds
(epistasis). True biological heterogeneity is hard to disentangle
from the LD and frequency effects above, but it accounts for some
non-zero fraction of the portability gap.

#figure(
  image("../../diagrams/lecture-13/09-prs-portability.svg", width: 95%),
  caption: [
    Polygenic-score portability across ancestries. A PRS trained on a
    European-ancestry GWAS loses three to five times its predictive
    accuracy when applied to non-European cohorts, primarily because
    of differing LD patterns and allele frequencies between training
    and test populations.
  ],
) <fig:portability>

The practical and ethical consequences are severe. A PRS-based
clinical tool trained exclusively on European data underperforms in
minority patient populations, which means using it in the clinic
risks worsening existing health disparities rather than narrowing
them. The technical fixes are known — *equally-sized GWAS in each
ancestry*, ancestry-aware PRS methods (PRS-CSx, PolyPred), and
diverse reference panels for imputation — but the sample-size
disparity is real and growing slowly. As of 2023, fewer than ten
percent of GWAS participants worldwide were of non-European ancestry,
and most of that fraction lives in East Asian biobanks. African,
Hispanic, and South Asian representation remains a small fraction
of total.

#warn[
  Do not deploy a European-trained PRS to a clinical cohort of
  different ancestry without recalibration and validation against a
  matching-ancestry test set. The clinical literature contains
  documented cases where the same PRS that is well-calibrated in
  Europeans systematically over- or under-predicts risk in
  African-ancestry patients, exactly because of the LD and
  allele-frequency mismatches described above. This is one of the
  small number of places in this book where a methodological mistake
  has unambiguous clinical-harm consequences.
]


== Fine-Mapping and Colocalization <sec:finemapping>

A GWAS peak is not a variant. It is a locus — a window of ten to a
few hundred kilobases — containing tens to hundreds of SNPs whose
allele frequencies in cases and controls are correlated through LD
with one or a few causal variants. *Fine-mapping* is the inverse
problem of recovering the causal variants from the peak.

=== From Peak to Causal Variant

The reason fine-mapping is hard is the same reason peaks have shape
in the first place: SNPs in LD have highly correlated test statistics.
The most-significant SNP in the peak — the "lead SNP" — is often not
the causal one. It is whichever SNP happens to be best tagged by the
genotyping array or imputation panel for the actual causal variant.
The causal SNP might not even be in the panel.

Classical fine-mapping used *stepwise conditional analysis*: take
the lead SNP, include it as a covariate, re-run the per-SNP regression
across the locus, and see whether any remaining SNP still exceeds
significance. Repeat until no SNP survives. This identifies the
number of approximately independent signals but does not localize any
of them sharply. Modern approaches reformulate fine-mapping as
Bayesian sparse regression and return a *credible set* per signal —
a set of SNPs guaranteed (under the model's assumptions) to contain
the causal variant with specified posterior probability.

=== #idx("SuSiE")SuSiE, #idx("FINEMAP")FINEMAP, and #idx("CAVIAR")CAVIAR

*SuSiE* (Sum of Single Effects; Wang et al., 2020) is the dominant
modern method. It models the phenotype at the locus as a sum of
$L$ single-effect components, each of which picks exactly one SNP:

$ bold(Y) = sum_(l=1)^L sum_(j=1)^M gamma_(l j) bold(X)_j beta_(l j) + bold(epsilon) $

The $gamma_(l j)$ is a one-hot indicator vector — within signal $l$,
exactly one SNP $j$ contributes. The posterior over which SNP that is
gives a *posterior inclusion probability* (PIP) per SNP per signal,
and the *credible set* for signal $l$ is the smallest set of SNPs
whose summed PIP exceeds a chosen threshold (typically $0.95$). The
optimization is a coordinate-ascent variational scheme that converges
in seconds per locus.

*FINEMAP* (Benner et al., 2016) takes a different route — a
stochastic-search Bayesian model selection over all possible
configurations of up to $K$ causal SNPs (typically $K = 5$). It
samples configurations weighted by their posterior probability and
reports per-SNP PIPs as marginals. Slower than SuSiE but more
flexible about the prior on the number of causal signals.

*CAVIAR* (Hormozdiari et al., 2014) was the first widely-used
Bayesian fine-mapper, with an exact enumeration over configurations
up to a small $K$. Computationally limited to small loci but
historically important.

All three methods take *summary statistics plus an LD matrix* as
input. They do not need individual-level genotypes. The LD matrix is
typically computed from the GWAS cohort itself if available, or from
a matched reference panel (1000 Genomes, UK Biobank) when only
summary statistics are public.

#figure(
  image("../../diagrams/lecture-13/10-fine-mapping.svg", width: 95%),
  caption: [
    Fine-mapping with SuSiE at a GWAS locus. The top track shows the
    $-log_(10)(p)$ signal; the middle track shows LD ($r^2$) to the
    lead SNP; the bottom track shows posterior inclusion probabilities.
    The 95 % credible set narrows the locus from about 200 peak SNPs
    down to a handful of candidates that collectively cover the
    causal variant with $0.95$ posterior probability.
  ],
) <fig:finemap>

#figure(
  image("../figures/ch13/f4-sparse-inverse.svg", width: 95%),
  caption: [
    Fine-mapping as a sparse linear inverse problem. The observed
    z-statistic vector $bold(z)$ at a locus equals the LD matrix
    $bold(R)$ times a sparse causal-effect vector $bold(beta)_c$,
    plus noise. Recovery of $bold(beta)_c$ from $bold(z)$ given
    $bold(R)$ is exactly compressed sensing — basis pursuit, LASSO,
    or $L_0$-regularized reconstruction depending on the prior.
    Wider LD blocks make the inverse problem more ill-conditioned.
  ],
) <fig:sparse-inverse>

#note[
  Fine-mapping is a sparse linear inverse problem. The observed
  per-SNP z-statistic vector at a locus is

  $ bold(z) = bold(R) bold(beta)_c sqrt(N) + bold(eta) $

  where $bold(R)$ is the LD matrix at the locus, $bold(beta)_c$ is
  the (mostly-zero) vector of true causal effects, and $bold(eta)$ is
  approximately white noise. Recovering $bold(beta)_c$ from $bold(z)$
  given $bold(R)$ is compressed sensing — basis pursuit, LASSO, or
  $L_0$-regularized reconstruction, depending on the prior. SuSiE's
  decomposition into a sum of single effects is an explicit sparse
  recovery; FINEMAP's posterior over configurations is Bayesian sparse
  recovery. The eigenstructure of $bold(R)$ — the LD block shape —
  determines whether the inverse is well-conditioned. The
  statistical-genetics community discovered this framework
  independently of the signal-processing community; the math is
  identical.
]

=== Colocalization

A GWAS signal at a locus tells you that *some* variant in the locus
affects the trait. It does not say which gene the variant acts
through, or how. *Colocalization* answers the second question: does
the GWAS signal share a causal variant with the signal from an
expression QTL (eQTL) for a nearby gene?

If yes, the gene whose eQTL co-peaks with the GWAS is a strong
candidate for being the mediating gene — the variant changes the
gene's expression, the changed expression affects the phenotype, and
the GWAS picks up the downstream phenotypic effect. If no, the
mediator is something else: a different gene, a protein-coding
change, a regulatory effect on a tissue not sampled by the eQTL
dataset.

The standard tool is *coloc* (Giambartolomei et al., 2014). At each
locus, it tests five mutually-exclusive hypotheses:

- $H_0$: no association with either trait.
- $H_1$: association with trait 1 only (GWAS).
- $H_2$: association with trait 2 only (eQTL).
- $H_3$: both associated, but with different causal variants.
- $H_4$: both associated, with the same causal variant.

Posterior probabilities of each hypothesis are computed from the two
sets of summary statistics under a sparsity prior on causal SNPs.
$P(H_4) > 0.8$ is the conventional threshold for calling
colocalization. Large eQTL resources — GTEx (Genotype-Tissue
Expression, ~50 tissues), eQTLGen (whole blood), CommonMind (brain) —
make colocalization a routine first step after a GWAS hit.
Empirically about 20–40 % of GWAS loci colocalize with at least one
tissue's eQTL.

#figure(
  image("../../diagrams/lecture-13/11-colocalization.svg", width: 95%),
  caption: [
    Colocalization of a GWAS signal with an eQTL. The GWAS peak for
    LDL cholesterol and the eQTL peak for `SORT1` expression in liver
    tissue coincide at the same locus. The coloc posterior
    $P(H_4) = 0.93$ indicates the two signals share a causal variant,
    implicating `SORT1` as the mediating gene.
  ],
) <fig:coloc>

=== When Fine-Mapping Fails

Three failure modes recur.

*Huge credible sets.* When the LD block is tight and the genuine
causal variant is one of many indistinguishable tags, the 95 %
credible set can contain tens to hundreds of SNPs. There is no
statistical resolution to be had at the current sample size and
imputation panel. Improving the panel — moving from array to whole-genome
sequencing, or using a denser reference like the pangenome —
sometimes shrinks the set; sometimes does not.

*Overlapping signals.* When the locus contains multiple causal
variants in moderate LD, SuSiE with $L = 5$ may report five signals
with credible sets that overlap heavily, and disentangling them
requires either more data or a stronger prior on the number of
signals. This is the regime where FINEMAP's stochastic search can
out-perform SuSiE's variational scheme.

*Structural causal.* Fine-mapping assumes a SNP-based causal model.
If the causal variant is actually a #idx("structural variant")structural variant — a #idx("CNV")CNV, a
short tandem repeat expansion, a transposable-element insertion —
that is not represented in the SNP panel, fine-mapping converges to
the SNP tagging the #idx("SV")SV and reports a credible set that does not
contain the actual cause. Empirically maybe 10–20 % of GWAS loci
have a plausible structural cause; long-read resequencing of GWAS
loci is the active area here.

#warn[
  A 95 % credible set contains the causal variant with 95 % posterior
  probability only if the prior assumptions hold — the assumed number
  of signals $L$, the LD matrix's accuracy, and the per-SNP effect
  prior. In practice $L$ is a hyperparameter, the LD matrix is
  estimated from a reference panel that may not match the GWAS cohort,
  and the effect prior is parametric. Credible sets are best read as
  "strong candidates" rather than as guaranteed #idx("coverage")coverage. Always
  inspect the locus plot — peak shape, LD pattern, and the relative
  PIPs — alongside the credible-set summary.
]


== Rare-Variant Association <sec:rare-variants>

Single-SNP GWAS works for common variants — minor allele frequencies
above roughly 1 %. Below that, the math breaks. The variance of
$hat(beta)$ at a single SNP scales as $1 / (N dot 2 p (1 - p))$, so a
variant at minor allele frequency $0.001$ delivers roughly $500$
times less information per individual than one at $0.5$. To detect a
$beta = 0.1$ effect at $p = 0.001$ with $alpha = 5 times 10^(-8)$ and
80 % power, the cohort needs to be of order $10^8$ individuals —
larger than any biobank that exists or is likely to exist.

The workaround is to give up on per-variant testing and aggregate
rare variants within some grouping — usually a gene — and test the
group collectively. The two dominant approaches are *burden tests*
and the *sequence kernel association test* (SKAT), with an adaptive
hybrid (*SKAT-O*) as the modern default.

=== Burden Tests

A #idx("burden test")burden test collapses all rare variants in a gene into a single
summary statistic per individual, then regresses the phenotype on
that summary. The simplest form is

$ "Burden"_i = sum_(j in "gene") w_j G_(i j) $

where $G_(i j)$ is the genotype of individual $i$ at rare variant $j$
inside the gene, and $w_j$ is a weight — often 1 for the simple
allele-count burden, or proportional to $1 / sqrt(p_j (1 - p_j))$ to
upweight rare variants, or chosen by functional-impact annotation
(only loss-of-function variants, only high-CADD-score missense, etc.).
The per-individual burden then enters a standard regression as a
single covariate:

$ Y_i = beta_0 + beta_"burden" dot "Burden"_i + bold(gamma)^T bold(X)_i + epsilon_i $

The hypothesis test is $H_0 : beta_"burden" = 0$, run once per gene.
For $approx 20,000$ protein-coding genes, the Bonferroni-adjusted
threshold is $0.05 / 20,000 = 2.5 times 10^(-6)$ — much less
stringent than the single-SNP genome-wide $5 times 10^(-8)$, because
there are far fewer hypothesis tests.

Burden tests work when most rare variants in the gene have effects
in the same direction — for example, when the gene is a
loss-of-function-intolerant tumour suppressor and the rare variants
are all damaging. They fail when effects are bidirectional: a gene
with both loss-of-function and gain-of-function rare variants will
see them cancel inside the sum, leaving no net signal.

=== SKAT and SKAT-O

*SKAT* (Wu et al., 2011) handles the bidirectional regime. Instead of
summing signed allele counts, SKAT computes a quadratic form in
genotypes:

$ Q = (bold(Y) - bold(mu))^T bold(G) bold(W) bold(G)^T (bold(Y) - bold(mu)) $

where $bold(G)$ is the $N times M$ rare-variant genotype matrix for
the gene and $bold(W)$ is a diagonal weight matrix. Under the null,
$Q$ follows a mixture of chi-squared distributions (the Davies
approximation gives the tail probability). The intuition: SKAT asks
whether the gene's genotype kernel explains more phenotype variance
than expected under the null, without requiring the variants to point
in the same direction. A gene with five protective and five damaging
rare variants would still register if both contributed real variance.

*SKAT-O* (Lee et al., 2012) is an adaptive combination of burden and
SKAT, parameterized by a mixing weight $rho in [0, 1]$ that
interpolates between the burden test ($rho = 1$) and the SKAT test
($rho = 0$). The optimal $rho$ is chosen by minimizing the resulting
p-value over a grid, with the multiple-testing penalty folded in. The
result is a test that approximates whichever of burden or SKAT is
better suited to the locus, without requiring the analyst to commit
in advance. SKAT-O is the default in most modern rare-variant
pipelines.

#figure(
  image("../../diagrams/lecture-13/12-burden-skat.svg", width: 95%),
  caption: [
    When burden, SKAT, and SKAT-O win. Burden tests aggregate
    directional signals and win when rare-variant effects are
    concordant (left). SKAT detects variance explained even when
    effects are bidirectional (middle). SKAT-O adapts to whichever
    regime dominates at each locus.
  ],
) <fig:burden-skat>

#note[
  A burden test is a #idx("matched filter")matched filter on rare-variant counts. The
  template is "sum of rare variants in gene $G$"; the signal is the
  phenotype; the test statistic is the regression coefficient on the
  templated burden. SKAT is the kernel version, using a quadratic
  form with kernel matrix $bold(K) = bold(G) bold(W) bold(G)^T$ —
  variance-based detection that does not require the signal to
  align with a specific template direction. SKAT-O is composite-
  hypothesis detection that does not commit in advance to whether
  the signal is matched-filter-shaped or kernel-shaped. Each of these
  is a standard pattern in classical detection theory; the
  bioinformatics labels are new, the underlying decision rules are
  old.
]

=== Data Sources

Rare-variant testing requires sequencing rather than genotyping
arrays, because arrays capture only the common variants their content
panel was designed to type. Three resources dominate.

*#idx("gnomAD")gnomAD* (Karczewski et al., 2020) aggregates exome and whole-genome
sequencing from public studies — over 800 thousand exomes and 150
thousand genomes by version 4. It is the reference for rare-variant
allele frequencies across continental ancestries; every clinical
genetics pipeline checks gnomAD for "how rare is this variant?"
before annotating it as a candidate.

*UK Biobank exome and WGS releases.* Exome sequences for all 500,000
UK Biobank participants were released in 2021; whole genomes finished
in 2023. The combination of phenotypic depth and genome-wide
sequencing makes UKB the largest source of rare-variant association
results across hundreds of traits.

*TOPMed.* The Trans-Omics for Precision Medicine programme
contributes about 150,000 whole genomes from US multi-ethnic cohorts.
Particularly useful for rare variants in non-European populations
that gnomAD and UKB undersample.

=== Caveats

A short list of failure modes specific to rare-variant tests:

- *Annotation dependence.* Burden tests require a definition of
  "qualifying" rare variants in each gene. Definitions vary in two
  axes: allele-frequency cutoff (0.01? 0.001? 0.0001?) and functional
  filter (loss-of-function only? Plus damaging missense? CADD threshold?).
  Sensitivity analyses across definitions are standard practice.
- *Calling error amplification.* False-positive rare-variant calls
  feed directly into the burden. A 0.1 % per-variant false-positive
  rate from sequencing can produce a noticeable inflation in the
  per-gene burden if not filtered out. Quality control on rare-variant
  calls is more demanding than on common variants.
- *Case-control imbalance.* For binary traits with severe imbalance
  (rare disease versus general-population controls), ordinary
  regression's saddle-point approximation breaks. SAIGE-GENE and
  REGENIE's gene-level extension fix this with the same
  saddle-point machinery the single-SNP versions use.
- *Cohort heterogeneity.* Aggregating rare variants across genuinely
  different ancestral populations dilutes signals; ancestry-stratified
  analysis plus meta-analysis is usually safer than pooled analysis.


== Summary <sec:summary>

- GWAS is a bank of per-SNP regressions running in parallel against a
  shared cohort. Linear regression for quantitative traits, logistic
  regression for case-control. At biobank scale this is a
  massively-multichannel detection problem on the order of $10^6$
  effectively independent tests.
- The genome-wide significance threshold of $5 times 10^(-8)$ is
  Bonferroni correction for the LD-independent count of tests, not
  the raw SNP count. The Manhattan plot displays the signal; the QQ
  plot checks calibration; the genomic-control inflation factor
  $lambda_("GC")$ is a one-number summary of the QQ plot.
- Population stratification inflates every SNP that differs in
  frequency between groups. PCA correction (top 10–40 PCs as
  covariates) is the first line of defence; mixed-model regression
  with kinship (BOLT-LMM, SAIGE, REGENIE) is the modern workhorse and
  the equivalent of GLS noise whitening.
- LD-score regression separates polygenic signal from confounding
  using summary statistics alone. Slope $prop$ heritability;
  intercept $prop$ confounding.
- Polygenic risk scores aggregate genome-wide signal into a
  per-individual prediction. They have clinical utility for several
  common diseases. Their portability across ancestries is poor when
  training data is European-dominated — a problem that is both
  technical and ethical.
- Fine-mapping is a sparse linear inverse problem. SuSiE, FINEMAP,
  and CAVIAR return credible sets of SNPs that collectively cover
  the causal variant with specified posterior probability.
- Colocalization links a GWAS locus to a candidate gene by matching
  the GWAS signal against an eQTL signal. $P(H_4) > 0.8$ in coloc is
  the standard threshold.
- Rare-variant association uses aggregation. Burden tests sum
  variants within a gene; SKAT runs a kernel test that handles
  bidirectional effects; SKAT-O adaptively combines them. gnomAD is
  the standard reference for rare-variant allele frequencies.


== Exercises <sec:exercises>

#strong[1.] #emph[Single-SNP power.]
Using the chi-squared non-centrality approximation,
$lambda_"ncp" = N dot 2 p (1 - p) dot beta^2$, compute the cohort
size $N$ required to detect $beta = 0.03$ (in trait-SD units) at
allele frequency $p = 0.2$ with 80 % power at the genome-wide
threshold $alpha = 5 times 10^(-8)$. Then repeat for $p = 0.05$,
$p = 0.01$, and $p = 0.001$. Plot $N$ against $log_(10)(p)$ on
log–log axes. At what allele frequency does $N$ exceed $10^7$, and
what does that say about single-SNP testing of rare variants?

#strong[2.] #emph[QQ-plot interpretation.]
You are given a QQ plot showing $lambda_("GC") = 1.18$ on a polygenic
trait at $N = 200,000$. LD-score regression reports slope
$0.6 / "Mb"$ and intercept $1.02$. Is the inflation primarily
polygenic signal or primarily confounding? Justify in two sentences.
What additional diagnostic would you run to decide whether the
study's results are trustworthy?

#strong[3.] #emph[Bonferroni arithmetic.]
A whole-exome rare-variant burden test scans 19,800 protein-coding
genes. Compute the gene-level Bonferroni threshold for FWER $alpha =
0.05$. The same study also runs single-SNP tests on the same exome,
which retain 600,000 variants after frequency filtering. Compute the
SNP-level threshold. Why might the gene-level threshold be the easier
one to clear despite testing on the same data?

#strong[4.] #emph[PCA covariate dimensionality.]
For a cohort of 50,000 individuals, you run a per-SNP regression
including 20 principal components as covariates. The effective sample
size for the per-SNP test drops by how much, in percentage terms?
At what number of PCs does the percentage drop exceed 5 %? Why is
this argument less important at $N = 500,000$?

#strong[5.] #emph[Reading a PRS curve.]
@fig:portability shows a PRS-vs-phenotype $R^2$ of 0.12 in the
European training cohort, 0.05 in East Asian, and 0.03 in African.
For a binary disease with population prevalence 5 %, translate the
$R^2$ values into approximate odds ratios for individuals in the top
10 % of PRS versus the bottom 10 %, separately by ancestry. (Use the
liability-threshold approximation; document any assumptions.)

#strong[6.] #emph[Credible set arithmetic.]
A SuSiE credible set has 5 SNPs with posterior inclusion probabilities
$(0.42, 0.31, 0.14, 0.06, 0.04)$. What is the credible-set coverage
probability? Is this a 95 % credible set? If you wanted to shrink the
set to 3 SNPs while maintaining $>= 0.85$ coverage, is it possible
without additional data? Justify.

#strong[7.] #emph[Burden versus SKAT.]
A gene contains ten rare variants in cases versus the same ten
variants in controls. In setup A, all ten variants have effects of
$+0.5$ in cases. In setup B, five variants have effect $+0.5$ and
five have effect $-0.5$ in cases. For each setup, sketch the expected
behaviour of a burden test and a SKAT test, and explain which test
wins where.

#strong[8.] #emph[(Open-ended.)]
Pick one of the modern mixed-model GWAS tools — BOLT-LMM, SAIGE, or
REGENIE — and read its primary publication. In one paragraph, describe
the single most consequential algorithmic shortcut the authors used to
make biobank-scale GWAS tractable, and identify one regime in which
the shortcut might fail.


== Further Reading <sec:further-reading>

- #strong[Wellcome Trust Case Control Consortium.] (2007). "Genome-wide
  Association Study of 14,000 Cases of Seven Common Diseases and
  3,000 Shared Controls." _Nature_ 447: 661–678. The big-bang paper.
- #strong[Price, A. L., Patterson, N. J., Plenge, R. M., et al.]
  (2006). "Principal Components Analysis Corrects for Stratification
  in Genome-Wide Association Studies." _Nature Genetics_ 38: 904–909.
  The PCA-correction paper.
- #strong[Yang, J., Lee, S. H., Goddard, M. E., and Visscher, P. M.]
  (2011). "GCTA: A Tool for Genome-Wide Complex Trait Analysis."
  _American Journal of Human Genetics_ 88: 76–82. SNP-heritability
  estimation from kinship.
- #strong[Loh, P. R., Tucker, G., Bulik-Sullivan, B. K., et al.]
  (2015). "Efficient Bayesian Mixed-Model Analysis Increases
  Association Power in Large Cohorts." _Nature Genetics_ 47: 284–290.
  BOLT-LMM.
- #strong[Zhou, W., Nielsen, J. B., Fritsche, L. G., et al.] (2018).
  "Efficiently Controlling for Case-Control Imbalance and Sample
  Relatedness in Large-Scale Genetic Association Studies."
  _Nature Genetics_ 50: 1335–1341. SAIGE.
- #strong[Mbatchou, J., Barnard, L., Backman, J., et al.] (2021).
  "Computationally Efficient Whole-Genome Regression for Quantitative
  and Binary Traits." _Nature Genetics_ 53: 1097–1103. REGENIE.
- #strong[Bulik-Sullivan, B. K., Loh, P. R., Finucane, H. K., et al.]
  (2015). "#idx("LD score regression")LD Score Regression Distinguishes Confounding from
  Polygenicity in Genome-Wide Association Studies." _Nature Genetics_
  47: 291–295. LDSC.
- #strong[Wang, G., Sarkar, A., Carbonetto, P., and Stephens, M.]
  (2020). "A Simple New Approach to Variable Selection in Regression,
  with Application to Genetic Fine-Mapping." _Journal of the Royal
  Statistical Society, Series B_ 82: 1273–1300. SuSiE.
- #strong[Giambartolomei, C., Vukcevic, D., Schadt, E. E., et al.]
  (2014). "Bayesian Test for Colocalisation Between Pairs of Genetic
  Association Studies Using Summary Statistics." _PLoS Genetics_ 10:
  e1004383. coloc.
- #strong[Wu, M. C., Lee, S., Cai, T., et al.] (2011). "Rare-Variant
  Association Testing for Sequencing Data with the Sequence Kernel
  Association Test." _American Journal of Human Genetics_ 89:
  82–93. SKAT.
- #strong[Karczewski, K. J., Francioli, L. C., Tiao, G., et al.]
  (2020). "The Mutational Constraint Spectrum Quantified from
  Variation in 141,456 Humans." _Nature_ 581: 434–443. gnomAD.
- #strong[Martin, A. R., Kanai, M., Kamatani, Y., et al.] (2019).
  "Clinical Use of Current Polygenic Risk Scores May Exacerbate
  Health Disparities." _Nature Genetics_ 51: 584–591. The PRS
  portability paper.
- #strong[Visscher, P. M., Wray, N. R., Zhang, Q., et al.] (2017).
  "10 Years of GWAS Discovery: Biology, Function, and #idx("translation")Translation."
  _American Journal of Human Genetics_ 101: 5–22. A decade in review.
