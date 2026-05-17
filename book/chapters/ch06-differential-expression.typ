#import "../theme/book-theme.typ": *

= #idx("differential expression")Differential Expression and Count Statistics <ch:differential-expression>

#matters[
  After Chapter 5 you have a count matrix: a rectangle of integers,
  twenty thousand rows by a handful of columns, with one row per gene
  and one column per sample. The matrix sits on disk. It is not yet
  biology. To turn it into biology you have to answer the question that
  every wet-lab collaborator asks first: _which genes are different
  between these conditions?_ The answer is a list of gene names with
  log-fold changes and adjusted p-values attached. Producing that list
  honestly — so the p-values mean what they claim to mean and the
  fold changes are not artefacts of library depth — is what this chapter
  is about. Get any part of it wrong and the list becomes a
  random-number generator dressed up in statistical clothing. Get it
  right and you have the substrate on which an entire pharmaceutical
  programme can be defended.
]

A modern differential-expression analysis runs twenty thousand
hypothesis tests in parallel with three samples per condition. Stated
that way, it should not work. The classical statistician's intuition
says you need many replicates to estimate variance, you need normally
distributed data to use a t-test, and you need a #idx("Bonferroni")Bonferroni penalty
that makes any single result implausible. None of those constraints
hold in modern #idx("RNA-seq")RNA-seq, and yet the field produces calibrated p-values
and reproducible gene lists routinely. The reason is a particular
combination of moves — a negative-binomial generalised linear model,
empirical-Bayes #idx("shrinkage")shrinkage of #idx("dispersion")dispersion across genes, and
Benjamini–Hochberg false-discovery-rate control — that together turn an
under-determined statistical problem into a tractable one.

This chapter walks each move. Section 6.1 sets up the question and
shows why the most natural first instinct — a per-gene t-test — fails
on count data. Section 6.2 introduces the negative-binomial
distribution and the generalised linear model that #idx("DESeq2")DESeq2 and #idx("edgeR")edgeR
fit. Section 6.3 covers the DE pipeline in detail: per-gene dispersion
estimation, the mean–dispersion trend, empirical-Bayes shrinkage, and
the choice between Wald and likelihood-ratio tests. Section 6.4 turns
to multiple-testing correction and the p-value-histogram diagnostic
every analyst should plot before reporting anything. Section 6.5 takes
a gene list and turns it into pathway-level claims via
over-representation analysis and the Subramanian gene-set-enrichment
running-sum statistic.

The chapter assumes you have a count matrix in hand and have computed
size factors — the upstream work covered in Chapter 5. Everything
downstream of the gene list — single-cell analyses, integrative
multi-omics, clinical interpretation — assumes you understand what is
and is not implied by an adjusted p-value of $5 times 10^(-3)$ in a
DESeq2 results table.


== The Differential-Expression Problem <sec:de-problem>

A biologist's question is straightforward. Two groups of samples,
typically three to six replicates each. One group received a treatment,
the other did not. Which genes responded? The answer is the list of
genes whose expression differs between the groups by an amount that
exceeds what you would expect from biological and technical noise
alone.

Stated as inference, the question becomes a per-gene hypothesis test.
For each gene $g$ in the genome:

- $H_0$: the gene's expression is the same in both conditions.
- $H_1$: the gene's expression differs by some non-zero amount.

The output is, conceptually, twenty thousand pairs of numbers: an
estimate of how different the expression is (the _log-fold change_),
and a confidence statement about that estimate (the _p-value_,
eventually adjusted for the fact that you ran the test twenty thousand
times). The canonical visual for this output is the *#idx("volcano plot")volcano plot* —
log-fold change on the horizontal axis, $-log_(10)$ of the adjusted
p-value on the vertical — which carves the gene cloud into a familiar
butterfly shape with the "hits" lifted far from the origin in both
directions.

#figure(
  image("../../diagrams/lecture-06/01-volcano-anatomy.svg", width: 95%),
  caption: [
    Anatomy of a volcano plot. Every DE result is a scatter in
    (effect, confidence) space; thresholds on both axes carve the
    plane into "hit" and "non-hit" regions. Most genes sit near the
    origin; the biologically interesting ones are pushed into the
    upper corners.
  ],
) <fig:volcano>

#note[
  A volcano plot answers two questions at once: how big is the effect
  (horizontal position) and how certain are we (vertical position).
  This chapter is mostly about how the vertical axis gets computed
  honestly. Getting the horizontal axis off by a factor of two is
  embarrassing; getting the vertical axis wrong turns the analysis into
  a random number generator.
]

=== Why the T-Test Fails on Counts <sec:why-t-fails>

A reasonable first instinct is to run a Student's t-test per gene. You
already trust t-tests in other domains; the test is one line in any
scripting language; the output is a p-value. Done. The instinct is
wrong, and it is worth understanding exactly why.

The first reason is _distributional_. A t-test assumes the data are
drawn from a normal distribution with some variance $sigma^2$. Counts
are not normal. They are non-negative integers, they have a hard floor
at zero, and their distribution at low expression levels is decidedly
asymmetric — a gene with true expression of half a read per sample will
have several samples with exactly zero observed reads and a long
right tail of occasional higher counts. A normal approximation centred
on that mean assigns positive density to negative values, which is
nonsense, and fails entirely to capture the zero spike.

The second reason is _sample size_. A typical bulk RNA-seq study has
three to six replicates per condition. The sample standard deviation
computed from three numbers has a relative error of roughly forty
percent — and that is best-case, assuming the data really are normal.
With twenty thousand genes tested, several thousand of them will have
spuriously small sample variances purely by chance, and a fixed
$p < 0.05$ threshold will flag them as "highly significant." A test
that is correct only on average becomes spectacularly wrong when you
run it twenty thousand times and pick the apparent winners.

The third reason is the _mean–variance relationship_. In real RNA-seq
data, the variance of a gene's counts grows roughly as
$sigma^2 approx mu + alpha mu^2$ with the gene's mean expression. A
t-test assumes homoscedastic variance — the same $sigma^2$ at every
mean. On RNA-seq data, the t-test will systematically underestimate
variance for highly expressed genes and overestimate it for lowly
expressed ones. The resulting p-values are wrong in opposite
directions at the two ends of the expression range.

#figure(
  image("../../diagrams/lecture-06/02-t-test-failure.svg", width: 95%),
  caption: [
    Why a t-test is the wrong detector for RNA-seq counts. The
    distribution is not normal, and the resulting p-values are
    systematically wrong under the null. Low-count genes — where
    the failure is worst — are the ones most often relevant to
    biology.
  ],
) <fig:t-test-failure>

#note[
  RNA-seq counts behave like photon-counting observations with extra
  multiplicative noise. The right detection statistic is not the same
  as for Gaussian-noise channels — you need a likelihood-ratio test
  built against the correct noise distribution, not a fixed-variance
  t-test. Using a t-test on counts is like using a #idx("matched filter")matched filter that
  assumes white noise when the actual noise is shot noise with a
  $1/f$ component: the detector runs, it produces numbers, and the
  numbers are wrong in systematic ways.
]

=== What the Noise Model Has to Capture <sec:noise-model-requirements>

Any honest DE model has to represent three empirical facts about
RNA-seq counts.

First, the data are _discrete and have a zero floor_. Genes with low
true expression produce many zero observations and a thin right tail.
A Gaussian model cannot represent this. The minimum required
distribution family is one parametrised on non-negative integers — the
Poisson, the #idx("negative binomial")negative binomial, the zero-inflated Poisson, and a few
exotic alternatives are the candidate set.

Second, _variance grows with the mean_. A gene at mean count one
behaves close to Poisson — variance equal to mean. A gene at mean
count one thousand has variance much larger than its mean; the excess
variance is what we call _overdispersion_. The variance-versus-mean
plot for any real RNA-seq dataset shows a clear quadratic component
at high expression levels.

Third, the overdispersion is _structured_. Genes at similar expression
levels have similar overdispersion. There is a smooth mean–dispersion
trend across the transcriptome — high-mean genes have low dispersion,
low-mean genes have high dispersion. This structure is the hook that
empirical-Bayes shrinkage hangs on in @sec:eb-shrinkage.

The negative binomial captures all three facts at the cost of one
extra parameter beyond the Poisson: a per-gene dispersion $alpha_g$
that sets the quadratic-variance term. That is the model DESeq2 and
edgeR fit. The rest of this chapter is how to fit it well when you
have only a few replicates per condition.


== Negative Binomial and the GLM <sec:nb-glm>

The *negative binomial* distribution emerges from a simple two-stage
generative story. Imagine that each gene has a true mean expression
that is itself a random draw from a gamma distribution; conditional on
that mean, the observed counts are Poisson. Marginalising over the
gamma gives a negative-binomial distribution for the counts. The
parameters of the result are the mean $mu$ and a dispersion $alpha$
that controls how spread the underlying gamma is. The variance comes
out as

$ "Var"[Y] = mu + alpha mu^2. $

The first term is the Poisson contribution — pure counting noise. The
second is the overdispersion — extra variability from real biological
or technical fluctuation in the rate. When $alpha = 0$ the distribution
collapses to a Poisson. When $alpha$ is large the variance dominates
and the distribution becomes heavy-tailed.

Three regimes are worth naming explicitly. Very-low-count genes have
mean below one, and Poisson noise dominates the variance — the
$alpha mu^2$ term is small no matter what $alpha$ is. Highly expressed
housekeeping genes have tight biological control and small $alpha$,
typically around $0.01$. Bursty regulatory genes, immune genes, or any
gene under stochastic control have $alpha$ closer to $0.3$ — heavy
overdispersion, low signal-to-noise ratio per replicate. A typical
RNA-seq dataset has dispersions spanning two orders of magnitude
across the twenty thousand genes.

#figure(
  image("../figures/ch06/f1-nb-vs-poisson.svg", width: 95%),
  caption: [
    Negative binomial and Poisson PMFs at three mean counts. At
    $mu = 0.5$ both distributions look almost identical; the zero
    spike dominates. At $mu = 5$ a slight excess in the right tail
    distinguishes the NB. At $mu = 50$ the NB has a visibly broader
    spread, with mass extending past 80 counts that the Poisson
    forbids.
  ],
) <fig:nb-vs-poisson>

#figure(
  image("../../diagrams/lecture-06/03-mean-dispersion.svg", width: 95%),
  caption: [
    The mean–dispersion relationship in real RNA-seq data. High-mean
    genes have low dispersion; low-mean genes have high dispersion.
    This structure is what makes empirical-Bayes shrinkage possible
    (@sec:eb-shrinkage).
  ],
) <fig:mean-dispersion>

#note[
  The dispersion parameter is the noise-floor characterisation of each
  gene's expression channel. High-dispersion genes are noisy signal
  sources where even a big effect is hard to distinguish from
  fluctuation; low-dispersion genes are clean channels where small
  effects are visible. The DE problem is detection on a channel whose
  noise level you have to estimate from the same few samples you are
  trying to detect in — with twenty thousand channels to estimate
  jointly.
]

=== The NB Generalised Linear Model <sec:nb-glm-model>

DESeq2 and edgeR fit a *generalised linear model* (GLM) per gene. For
gene $g$ and sample $i$, the canonical form is

$ log mu_(i g) = log s_i + beta_(0,g) + beta_(1,g) dot x_i $

where $mu_(i g)$ is the expected count for gene $g$ in sample $i$;
$s_i$ is the *#idx("size factor")size factor* for sample $i$ (estimated by the
median-of-ratios procedure from Chapter 5); $x_i$ is a condition
indicator (zero for control, one for treatment); $beta_(0,g)$ is the
gene's baseline log-expression; and $beta_(1,g)$ is the
*log-fold change* between conditions — the quantity we ultimately
test against zero.

Real designs go far beyond a single binary indicator. The right-hand
side of the GLM is a sum of any number of $beta_(k,g) dot x_(i,k)$
terms encoding batch, sex, age, treatment, time, and their
interactions. The machinery is unchanged. Once the GLM is fit, every
DE question reduces to "is some linear combination of $beta$
coefficients significantly different from zero?" — a standard
hypothesis test about an estimated parameter.

#figure(
  image("../../diagrams/lecture-06/04-nb-glm-design.svg", width: 95%),
  caption: [
    The NB GLM in practice. A design matrix encoding condition and
    batch, plus one gene's counts, plus the log-link equation,
    produces fitted coefficients whose condition column is the
    log-fold change. The same fit happens twenty thousand times,
    once per gene.
  ],
) <fig:nb-glm-design>

#note[
  The NB GLM is maximum-likelihood inference with a log link —
  structurally identical to logistic regression (Bernoulli likelihood,
  logit link) or Poisson regression (Poisson likelihood, log link).
  The dispersion $alpha$ replaces the scale parameter you would have
  in ordinary least squares. Every DE tool whose name ends in `-seq2`
  or `-R` ultimately fits this model; the tools differ in how they
  estimate $alpha$ and which form of test they use for $beta$.
]

=== Size Factors and Log-Fold Change <sec:size-factors>

The size factor $s_i$ enters the GLM as an *offset*, not a coefficient.
We treat it as known — it was estimated from the full count matrix in
Chapter 5. In log space, the $log s_i$ term shifts the prediction up
for samples with deeper libraries and down for samples with shallower
ones, leaving the $beta$ coefficients to carry the
composition-normalised expression change. Without the offset, the
condition coefficient would conflate "more reads in this sample" with
"more expression of this gene" — exactly the composition confounding
the median-of-ratios procedure was designed to remove.

After the offset has absorbed depth, the coefficient $beta_(1,g)$ is
the log-fold change between conditions:

$ beta_(1,g) = log(mu_(g, "treated") / mu_(g, "control")). $

Converting to the biologist's convention of $log_2$ #idx("fold change")fold change is one
division by $log(2)$. A $beta_(1,g)$ value of $0.693$ corresponds to a
two-fold increase; a value of $-1.386$ corresponds to a four-fold
decrease.

#tip[
  A size factor is not a per-gene multiplier. It is a per-sample
  multiplier applied to _all_ genes in that sample. Think of it as
  the sample's "sequencing depth in composition-adjusted units." A
  sample with size factor $1.2$ contributed twenty percent more
  reads than the population average; every gene's prediction in
  that sample is scaled up by $1.2 times$ before being compared with
  the observed count. Size factors are the bridge between a raw
  count and a comparable-across-samples abundance.
]

=== Biological versus Technical Variance <sec:bio-vs-tech>

The overdispersion captured by $alpha$ is not a single physical
quantity. It has two sources, and distinguishing them matters for
experimental design.

*Technical variance* comes from the library-prep and sequencing
pipeline. Pipetting error, polyadenylation bias, fragmentation noise,
unmarked PCR duplicates, sequencer calibration drift — all contribute.
This component shrinks with more reads per sample. If your problem is
technical, you can fix it by sequencing deeper.

*Biological variance* comes from genuine differences between replicate
samples that you were treating as nominally identical. In animal
experiments, slight differences in genetic background, circadian
phase, or cage environment. In nominally identical cell-line samples,
cell-cycle phase heterogeneity and clone-to-clone drift. This
component shrinks only with more biological replicates. More reads per
sample does not fix it.

In practice, biological variance dominates for medium-to-highly
expressed genes; technical variance dominates only at the bottom of
the dynamic range, for genes near the detection limit. The
consequence for design is unforgiving: if your effect is small or you
care about lowly expressed regulators (#idx("transcription")transcription factors, signalling
molecules, splicing factors), you need six or more biological
replicates. More sequencing depth on three samples will not save you.

#warn[
  Published DE analyses with two replicates per condition are common
  and mostly unreliable. DESeq2 and edgeR still _run_ on $n = 2$, but
  the dispersion estimates are noise, the empirical-Bayes shrinkage
  does its best, and the p-values come out with wide and often
  miscalibrated confidence intervals. If you see "n = 2 per group"
  in a methods section, treat the gene list as exploratory and
  prioritise wet-lab validation.
]


== Inside DESeq2 and edgeR <sec:deseq-edger>

Both DESeq2 and edgeR follow the same conceptual pipeline, differing
mostly in how they execute steps 3 through 5.

#figure(
  image("../figures/ch06/f2-de-pipeline.svg", width: 95%),
  caption: [
    The DE pipeline end to end, from the count matrix produced in
    Chapter 5 to the annotated gene list at the end of Section 6.5.
    Every substantive disagreement between DESeq2 and edgeR is about
    dispersion estimation and test choice; the rest is bookkeeping.
  ],
) <fig:de-pipeline>

The pipeline runs: estimate per-sample size factors (median of ratios,
Chapter 5); estimate per-gene dispersion as a maximum-likelihood
estimate from the residuals; fit a smooth mean–dispersion trend across
all genes; shrink the per-gene dispersions toward the trend via
#idx("empirical Bayes")empirical Bayes; refit the GLM per gene with the shrunk dispersion;
test each $beta$ coefficient against zero using a Wald or
likelihood-ratio test; apply Benjamini–Hochberg correction across the
genes; filter and annotate the survivors.

Two of those steps — dispersion shrinkage and #idx("FDR")FDR control — are where
the real action is. Everything else is well-understood textbook
machinery.

=== Per-Gene Dispersion #idx("maximum likelihood")Maximum Likelihood <sec:disp-mle>

Given the counts for gene $g$ across $n$ samples and the fitted
intercept $beta_(0,g)$, the per-gene maximum-likelihood estimate of
$alpha_g$ is the value that makes the observed counts most probable
under the NB model. Numerically this is found by a few iterations of
Newton–Raphson on the NB log-likelihood. DESeq2 implements this in
its `nbinomDispEstimate` routine; edgeR uses `estimateGLMTagwiseDisp`.

The trouble is that at three to six samples per condition the
likelihood surface for $alpha$ is broad and flat. The single-gene MLE
has a confidence interval that spans more than two orders of
magnitude. With twenty samples per condition the same estimator
converges to a tight peak. With three samples, you are essentially
guessing.

#figure(
  image("../../diagrams/lecture-06/05-dispersion-mle.svg", width: 95%),
  caption: [
    Per-gene dispersion MLE from few samples. Three replicates produce
    an essentially flat log-likelihood in $alpha$; twenty replicates
    produce a tight peak. Statistical sharing is what rescues small-$n$
    designs.
  ],
) <fig:dispersion-mle>

Two practical details matter even before shrinkage. First, dispersion
is estimated by _pooling residuals across both conditions_ for the
same gene, not separately within each condition. The assumption is
that the dispersion is the same in both conditions — only the mean
differs. This is the same intuition that pools the variance in a
two-sample t-test. Second, a single aberrant count in a gene with
three replicates per condition will inflate the gene-wise MLE by an
order of magnitude. DESeq2 detects such outliers using Cook's
distance — the same leverage-and-residual statistic from ordinary
least squares — and either replaces them with imputed counts or
excludes the affected gene from testing.

#warn[
  The gene-wise dispersion estimate from three replicates is not a
  reliable number on its own. If you see a per-gene dispersion value
  reported in a DE pipeline output, it is almost certainly the
  _shrunk_ estimate after the next subsection's procedure. If it is
  the raw MLE, treat it with deep suspicion.
]

=== Empirical-Bayes Shrinkage <sec:eb-shrinkage>

Here is the central statistical idea of DESeq2 and edgeR, and the
reason they work on three-replicate studies.

The empirical observation, repeated in essentially every RNA-seq
dataset ever published: if you plot per-gene mean against per-gene
MLE dispersion for all twenty thousand genes, the points form a
clear downward-sloping cloud. Low-expression genes have high
dispersion on average; high-expression genes have low dispersion on
average. The trend is real, smooth, and consistent across datasets.

The algorithmic move is to fit a parametric trend through the cloud
and then to pull each gene's raw MLE toward the trend line. DESeq2
fits a parametric form

$ alpha("trend")(mu) = alpha_0 + alpha_1 / mu $

with $alpha_0$ and $alpha_1$ estimated from the data; edgeR fits a
cubic-smoothing-spline variant. Once the trend is in hand, each gene's
dispersion estimate is a weighted average of its own raw MLE and the
trend value at its mean. The weight depends on how informative the
gene-wise estimate is — that is, how concentrated the per-gene
likelihood is around its maximum — and on how far the gene's MLE
sits from the trend. Genes whose MLE is close to the trend get pulled
all the way in; genes whose MLE is genuinely far from the trend resist
the pull.

#figure(
  image("../../diagrams/lecture-06/06-eb-shrinkage.svg", width: 95%),
  caption: [
    Empirical-Bayes dispersion shrinkage. Raw per-gene MLEs are pulled
    toward the fitted mean–dispersion trend to produce the shrunk
    estimates. Trend-consistent genes shrink a lot; outliers resist.
    The shrunk estimate has much lower variance than the raw MLE, at
    the cost of a small bias toward the trend.
  ],
) <fig:eb-shrinkage>

The statistical payoff is large. The shrunk estimator has much smaller
variance than the raw MLE for genes near the trend, at the cost of a
small bias toward the trend for genes that are genuinely atypical. In
large-scale testing that tradeoff is strongly favourable: the variance
reduction buys calibrated p-values; the small bias does not
qualitatively change which genes are significant. This is the move
that makes three-replicate studies tractable.

#note[
  Empirical-Bayes shrinkage is regularisation toward a learned prior.
  The "prior" is the fitted mean–dispersion trend, learned from the
  data being tested. Shrinking each gene's MLE toward the prior is the
  frequentist analog of ridge regression: trade a tiny amount of bias
  for a big reduction in variance. The same structural move shows up
  in James–Stein estimation (shrink $k >= 3$ means toward their grand
  mean), in best linear unbiased prediction in genetics, and in every
  hierarchical Bayesian model ever written. What makes the RNA-seq
  application distinctive is that the prior is itself empirical —
  fitted from the data being tested, not specified in advance.
]

The history is worth one paragraph. The 2010 edgeR paper (Robinson,
McCarthy, and Smyth) introduced shared-dispersion estimation for
RNA-seq, building on Smyth's earlier work on moderated t-statistics
for microarrays. The 2010 DESeq paper (Anders and Huber, then at EMBL)
proposed a slightly different shrinkage scheme and introduced the
median-of-ratios normalisation that the field still uses. The 2014
DESeq2 paper (Love, Huber, and Anders) unified the two approaches,
added empirical-Bayes shrinkage of fold-change estimates as well as
dispersion, and produced the implementation that, a decade later, is
still the community default. Its p-value calibration on realistic
data — visible in well-shaped histograms — is what won the field over.
As of the mid-2020s, the vast majority of published bulk RNA-seq DE
analyses use DESeq2.

=== Wald versus Likelihood-Ratio Tests <sec:wald-vs-lrt>

With the GLM fit and the shrunk dispersion in hand, the test against
$H_0 : beta_(1,g) = 0$ has two canonical forms.

The *#idx("Wald test")Wald test* computes a z-statistic from the fitted coefficient and
its standard error:

$ z_g = hat(beta)_(1,g) / "SE"(hat(beta)_(1,g)), quad p_g = 2 dot Phi(-|z_g|). $

It is simple, fast, and asymptotically valid. It works well for
single-coefficient tests — the two-condition comparison that
constitutes most DE work — and it is the default in DESeq2.

The *likelihood-ratio test* compares the full model fit (with
$beta_(1,g)$) to a reduced fit (without $beta_(1,g)$) by computing the
log-likelihood difference:

$ lambda_g = 2 dot (cal(l)_g^"full" - cal(l)_g^"reduced") tilde chi^2_1. $

The LRT is more robust for _complex_ hypotheses — testing several
coefficients jointly, as you would for a factor with three or more
levels or for a condition-by-time interaction — and for small-sample
cases where the Wald's linear approximation breaks down. DESeq2
exposes the LRT via the `test = "LRT"` argument; edgeR's quasi-likelihood
F-test (`glmQLFTest`) is a closely related construction.

#figure(
  image("../../diagrams/lecture-06/07-wald-vs-lrt.svg", width: 95%),
  caption: [
    Wald versus likelihood-ratio. The Wald test asks "is $hat(beta)$
    far from zero relative to its standard error?". The LRT asks
    "does removing $beta$ from the model materially hurt the fit?".
    The two agree for simple cases and large samples; the LRT wins for
    multi-coefficient hypotheses and curved likelihoods.
  ],
) <fig:wald-vs-lrt>

#note[
  The Wald and LRT tests are the same idea expressed at different
  points in the inference geometry. Wald looks at the curvature of the
  log-likelihood at the maximum and assumes the surface is well
  approximated by its quadratic Taylor expansion. The LRT looks at the
  actual height of the surface at two relevant points and asks how
  much falls off when you slide from one to the other. When the
  likelihood is well-approximated by a quadratic, the two agree. When
  the likelihood is curved or asymmetric — small samples, dispersions
  near the boundary, edge effects — only the LRT keeps its calibration.
]

=== Limma-Voom: The Variance-Stabilising Alternative <sec:limma-voom>

A third family of methods, headed by *limma-voom* (Law, Chen, Shi, and
Smyth, 2014), reaches the same destination by a different route. The
strategy is to transform the count data so it looks Gaussian-enough
for standard linear-model machinery to apply, rather than to fit an
explicit NB model.

The transformation has two pieces. First, log-counts-per-million
(log-CPM) the data — put it on an additive, approximately
homoscedastic scale. Second, apply a per-observation precision weight
estimated from the empirical mean–variance relationship: small counts
get small weights, large counts get large weights. Feed the weighted
log-CPMs into #idx("limma")limma's existing empirical-Bayes-moderated t-test, a
piece of machinery developed for microarrays in the early 2000s and
already well-tested.

The output is conceptually identical to DESeq2 and edgeR: per-gene
log-fold changes, moderated t-statistics, and p-values. Benchmarks
consistently show limma-voom is competitive with the GLM-based tools
on most datasets, and substantially faster on cohorts of thousands of
samples where the per-gene Newton–Raphson optimisation in DESeq2
begins to dominate runtime. The pragmatic rule for the mid-2020s: for
small experiments (twenty or fewer samples per condition), DESeq2 is
the default. For large cohorts (fifty or more per condition) or
genome-scale studies, limma-voom scales better and gives equivalent
answers.

#note[
  #idx("voom")Voom's variance stabilisation is a _preprocessing_ trick rather than
  a new likelihood. It transforms the data onto a scale where a
  constant-variance Gaussian model is approximately right, then
  applies constant-variance Gaussian inference. This is the same
  philosophy as the Anscombe transform for Poisson data or the
  Box–Cox family of power transforms. When the noise has a known
  variance function, transform the data to stabilise variance, then
  use off-the-shelf tools. The tradeoff is a small loss of efficiency
  in exchange for huge computational and methodological simplicity.
]


== #idx("multiple testing")Multiple Testing at Genome Scale <sec:multiple-testing>

A standard DE run tests around twenty thousand genes in parallel. Even
if every gene were truly null, the conventional $p < 0.05$ threshold
applied per gene would produce $20000 times 0.05 = 1000$ "significant"
genes by chance alone. No biologist will believe that any single raw
$p < 0.05$ means anything in this regime, and rightly so. The job of
multiple-testing correction is to control the rate of false positives
at the level of the gene _list_, not the individual gene.

Before getting into the procedures, there is one diagnostic plot every
analyst should produce: the *p-value histogram*. Under the null
hypothesis, p-values are distributed uniformly on $[0, 1]$. In a real
DE analysis most genes are null but some have real effects. The
histogram should show a flat bulk across most of the range (the
nulls), with a spike near zero (the real effects). If you see a
U-shape, a systematic slope, or a bump in the middle of the range,
the p-values themselves are miscalibrated — a sign of model
misspecification — and no amount of multiple-testing correction will
rescue them.

#figure(
  image("../../diagrams/lecture-06/08-pvalue-histogram.svg", width: 95%),
  caption: [
    The p-value histogram diagnostic. Healthy DE has a flat bulk with
    a spike near zero. Any other shape signals trouble in the model,
    not in the biology. Always plot this before reporting results.
  ],
) <fig:pvalue-histogram>

#tip[
  Plot the p-value histogram before you plot anything else. A flat
  null bulk with a near-zero spike means your testing framework is
  well-calibrated and there is real signal to report. A U-shape, a
  bump in the middle, or a systematic slope means the p-values are
  wrong — the model is misspecified, you have an unmodeled batch
  effect, or your size factors are off. Always plot the histogram
  before the volcano.
]

=== Bonferroni and Benjamini–Hochberg <sec:bonferroni-bh>

The classical correction is *Bonferroni's*. It controls the
*family-wise error rate* (FWER) — the probability of making any false
discovery at all across the whole experiment. With $m$ tests at
overall level $alpha$, Bonferroni rejects only those with
$p <= alpha / m$. For twenty thousand genes at $alpha = 0.05$, the
per-gene threshold becomes $2.5 times 10^(-6)$. The procedure is
provably valid, but it is _conservative to the point of uselessness_
in the regime where most tests are not null: you get a handful of
very-strong hits and miss most real effects.

*Benjamini and Hochberg's* 1995 procedure controls a different
quantity: the *false-discovery rate* (FDR), the expected fraction of
false discoveries among the rejected hypotheses. The procedure is
disarmingly simple. Sort the p-values in ascending order, $p_((1)) <=
p_((2)) <= dots <= p_((m))$. Find the largest rank $k$ such that

$ p_((k)) <= (k / m) dot alpha. $

Reject the $k$ hypotheses with the smallest $k$ p-values. The
adjusted p-value attached to gene at rank $k$ is the smallest level
$alpha$ at which that gene would have been rejected: $"padj"_((k)) =
min(p_((j)) dot m / j)$ over $j >= k$.

The interpretation is operational: a rejected list at $"padj" <= 0.05$
contains, in expectation, at most five percent false positives. If
your gene list has 500 entries at $"padj" <= 0.05$, then approximately
25 of them are expected false positives — an acceptable trade for
biological follow-up where wet-lab validation will filter false hits
later.

#figure(
  image("../../diagrams/lecture-06/09-bh-procedure.svg", width: 95%),
  caption: [
    The Benjamini–Hochberg procedure. Sorted p-values are plotted
    against rank; the BH line $p = (k/m) dot alpha$ is overlaid; the
    largest rank where the sorted p-value lies below the line gives
    the rejection cutoff. The same procedure produces wildly different
    rejection counts depending on signal strength in the data.
  ],
) <fig:bh-procedure>

#figure(
  image("../figures/ch06/f3-fdr-vs-fwer.svg", width: 95%),
  caption: [
    Bonferroni and Benjamini–Hochberg on the same simulated DE
    experiment. Both procedures control the rate of errors, but they
    control different rates and produce wildly different gene-list
    sizes. Bonferroni rejects only the bluest of the blue points; BH
    rejects everything below the line and produces a useful list
    whose false-positive fraction is the controlled quantity.
  ],
) <fig:fdr-vs-fwer>

#note[
  FDR control is detection theory with a controlled false-alarm _rate_
  rather than a controlled false-alarm _probability_. Bonferroni is
  the Neyman–Pearson worst-case bound — guaranteed that no false
  positive contaminates the whole experiment. BH's FDR is closer to
  the radar-engineering convention of "keep the fraction of false
  alarms in the detection list below a controlled rate" — more useful
  in regimes where you are going to follow up on the detections anyway
  and a small fraction of false positives is tolerable. Both have
  legitimate uses; DE analysis almost always wants FDR.
]

The 1995 Benjamini–Hochberg paper is one of the most-cited statistics
papers of the late twentieth century, but it sat largely unused in
biology for a decade. High-throughput microarrays in the early 2000s
made $m$-scale multiple testing routine, and BH became the default
correction by 2010. By the mid-2020s, every DE tool ships BH as the
default and there is no debate.

=== FDR Interpretation and q-Values <sec:fdr-qvalues>

The single number to read from a DE tool's output is the *adjusted
p-value* — `padj` in DESeq2, `FDR` in edgeR. Its operational meaning
is precise: if you reject all genes with $"padj" <= q$, then in
expectation $q$ times the number of rejections of them are false
positives.

A small refinement, the *q-value* (Storey 2003), estimates the
proportion of null hypotheses $pi_0$ directly from the p-value
distribution — specifically, from the height of the flat
"null bulk" of the histogram. The q-value is then $"padj" / pi_0$,
slightly less conservative than the raw BH adjustment because BH
implicitly assumes $pi_0 = 1$. For most bulk RNA-seq the two agree to
within a few percent. The distinction matters more in dense-signal
regimes — large-cohort phenotype-association studies — than in typical
DE work.

The practical rule is to report `padj` (BH) everywhere. Report a
q-value in addition if you are doing methodological benchmarking. Do
not report raw p-values without an FDR-adjusted companion — at large
$m$, raw p-values without correction are meaningless and routinely
misinterpreted by readers.

=== Independent Filtering and Weighted FDR <sec:indep-filtering>

DESeq2 applies one further refinement by default: *independent
filtering*. The observation is that genes with very low total counts
have essentially zero chance of producing a significant p-value even
if their true log-fold change is large — the variance is so high that
no test can discriminate the effect from noise. Including those genes
in the BH adjustment wastes power: they contribute to the denominator
$m$ but never to the numerator of rejections.

DESeq2 filters out low-count genes before applying BH, using a
data-dependent threshold that maximises the number of rejections. The
filter criterion (mean normalised count, or a more robust variant)
must be _independent of the test statistic under the null_ — otherwise
the procedure biases p-values. Mean count satisfies this requirement
because under $H_0$ the mean and the test statistic are independent.
The effective $m$ for BH shrinks by thirty to fifty percent in
typical datasets, the threshold becomes less conservative, and more
true positives survive.

#warn[
  If you see genes with `NA` in the `padj` column of a DESeq2 results
  table, they were filtered out by the independent-filtering step.
  This is usually fine: they were almost certainly going to be
  non-significant anyway. But check — the filtered set should be
  dominated by low-count genes. If high-count genes are getting
  `NA` padj values, something has gone wrong with the filter and you
  should investigate.
]

The more sophisticated *Independent Hypothesis Weighting* (IHW;
Ignatiadis et al. 2016) generalises the idea: weight each p-value by
an independently-estimated covariate (gene length, mean expression,
#idx("chromosome")chromosome) that predicts power without biasing the null. Used sparingly
in routine RNA-seq; common in large-cohort studies where the extra
power matters.


== From Gene List to Biology <sec:enrichment>

A DE result is a list of 500 to 2000 genes, each with a log-fold
change and an adjusted p-value. That list is rarely meaningful on
its own. The biological question is not "did `IFIT1` change?" — it is
"did the interferon-response programme turn on?". The transition from
individual genes to coordinated programmes is *enrichment analysis*.

Two families of methods dominate. *Over-representation analysis*
asks whether a curated gene set is over-represented in the
thresholded DE list. *Gene-set enrichment analysis* asks whether the
members of a gene set cluster toward one end of the full ranked
gene list, without thresholding. Both have uses; the right choice
depends on the data and the question.

=== Over-Representation Analysis <sec:ora>

Given a curated gene set — say MSigDB's "G2/M cell-cycle" collection,
127 genes annotated from the literature — and the DE-gene list from
your analysis — say 500 genes at $"padj" < 0.05$ — *over-representation
analysis* asks: is the overlap between these two sets larger than
chance would predict?

The contingency formulation is a 2-by-2 table. Each gene falls in one
of four cells defined by two indicators: is it in the curated gene
set or not, and is it in the DE list or not. With twenty thousand
total tested genes, 127 in the set, 500 in the DE list, and 42 in
both, the table looks like

#align(center)[
#table(
  columns: (auto, auto, auto, auto),
  align: (left, center, center, center),
  stroke: 0.5pt,
  table.header[][*in gene set*][*not in set*][*total*],
  [*DE gene*], [42], [458], [500],
  [*not DE*], [85], [19,415], [19,500],
  [*total*], [127], [19,873], [20,000],
)
]

Under the null hypothesis that DE genes are a random sample of the
twenty thousand tested, the number of DE genes falling in the gene
set follows a #idx("hypergeometric")hypergeometric distribution. #idx("Fisher's exact test")Fisher's exact test
computes the tail probability of seeing 42 or more set members in the
500-gene draw — the over-representation p-value. Equivalent
chi-square approximations exist; the exact form is preferred when any
table cell is small. Run the same test against thousands of gene sets
in MSigDB or the #idx("Gene Ontology")Gene Ontology, apply BH correction to the resulting
p-values, and report the enrichments that survive.

#figure(
  image("../../diagrams/lecture-06/10-ora-fisher.svg", width: 95%),
  caption: [
    Over-representation analysis via Fisher's exact test. The 2-by-2
    contingency of (DE, in-set) counts is compared to the
    hypergeometric null; the tail beyond the observed overlap gives
    the enrichment p-value.
  ],
) <fig:ora-fisher>

#idx("ORA")ORA is cheap, intuitive, and easy to explain to a wet-lab
collaborator. Its drawbacks are also straightforward. It requires a
hard threshold to define "DE" — sub-threshold genes contribute nothing
even when their fold changes are large and consistent. It treats every
gene above threshold as equivalent — a gene with $log_2$ fold change
of six counts the same as one with fold change of $0.6$. And it
discards the directional information: an enrichment driven by ten
upregulated genes and one downregulated gene looks the same as one
driven by five and six. Each of these problems has a fix; the
ranking-based #idx("GSEA")GSEA family addresses all three.

=== Gene-Set Enrichment Analysis <sec:gsea>

*GSEA* (Subramanian et al. 2005) abandons the threshold and uses the
full ranked gene list. The algorithm in four steps.

First, rank all twenty thousand genes by a continuous DE statistic —
the signed Wald z-score is a good default; signed
$-log_(10)(p) dot "sign"(log "FC")$ also works. Most-upregulated at
one end, most-downregulated at the other.

Second, for each gene set $S$ of interest, walk down the ranked list
and compute a running sum. When you encounter a gene in $S$, add a
positive weight proportional to the gene's ranking statistic. When you
encounter a gene not in $S$, subtract a small constant designed so
that the walk would return to zero if $S$ were a random subset.

Third, the *enrichment score* (ES) is the maximum absolute deviation of
the running sum from zero across the walk. A positive ES means the
set's members cluster at the top of the list (over-expressed in the
condition); a negative ES means they cluster at the bottom.

Fourth, the statistical significance of the ES is computed by
*permutation*. Shuffle the gene labels (in the gene-permutation
variant) or the sample labels (in the phenotype-permutation variant)
many times, recompute the ES for each shuffled dataset, and build the
null distribution. The empirical $p$ is the fraction of permutations
producing an ES as extreme as the observed one. With twenty thousand
genes and standard gene-set sizes, a thousand permutations is
typically sufficient.

#figure(
  image("../../diagrams/lecture-06/11-gsea-running-sum.svg", width: 95%),
  caption: [
    The GSEA running-sum statistic. Walk down the ranked list adding
    weight for set members and subtracting for non-members. The peak
    deviation from zero is the enrichment score; permutation gives
    the null distribution.
  ],
) <fig:gsea-running-sum>

#figure(
  image("../figures/ch06/f4-ora-vs-gsea-decision.svg", width: 95%),
  caption: [
    ORA and GSEA on the same DE result, side by side. ORA tests
    only thresholded genes against a contingency table; GSEA walks
    the full ranked list. The two methods often agree on strong
    enrichments and diverge on diffuse, low-amplitude programmes —
    where GSEA's pooled-evidence design wins.
  ],
) <fig:ch06-ora-vs-gsea>

#note[
  GSEA's running-sum statistic is a CUSUM (cumulative-sum) detector —
  the same test used for change-point detection in time-series
  monitoring. The "change" it detects is a sub-population of labelled
  elements clustered somewhere in the ranked list. Permutation makes
  the test nonparametric: no assumption on the null distribution of
  the enrichment score, just empirical calibration via shuffling. The
  algorithm appears in classical signal-processing literature as the
  rank-sum-runs test.
]

A related method, *CAMERA* (Wu and Smyth, 2012), addresses one
practical issue with naive permutation-based GSEA: genes in a real
biological pathway are not independent — they tend to be co-regulated
— and the permutation null treats them as if they were. Co-regulation
inflates Type I error in tests that ignore it. CAMERA estimates an
inter-gene correlation from the data and adjusts the variance of the
enrichment statistic accordingly. The resulting p-values are
better-calibrated than naive GSEA's at the cost of a more
elaborate computation. For routine analysis the difference is small;
for studies with strong co-regulation (developmental programmes,
synchronised cell cycle) CAMERA is worth the extra steps.

=== Pathway Databases <sec:pathway-databases>

The quality of any enrichment analysis is bounded by the quality of
the gene sets it runs against. The major sources are worth knowing.

*MSigDB* — the Molecular Signatures Database from the Broad Institute
— is the single most widely used collection. Its Hallmark gene sets
(`H`) are a hand-curated set of fifty high-level pathway signatures
distilled from the larger collections; the curated `C2` collection
draws gene sets from canonical pathway databases and published
expression signatures; `C7` is immunology-focused. Roughly twenty-five
thousand gene sets in total.

*Gene Ontology* (GO) is a controlled vocabulary of biological
processes, molecular functions, and cellular components, each with an
associated gene set. Universal across species and hierarchical — an
"immune response" annotation entails membership in many child
categories.

*#idx("KEGG")KEGG* — the Kyoto Encyclopedia of Genes and Genomes — provides
manually curated metabolic and signalling pathways, each with a
canonical gene list and a wiring diagram. Smaller and higher-quality
than GO; especially useful for metabolic-focused analyses.

*#idx("Reactome")Reactome* offers more detailed biochemistry than KEGG and is popular
in European labs and in cardiovascular and immunological research.
*WikiPathways* is community-curated and edited like Wikipedia.

The pragmatic rules: start with the Hallmark sets for a quick
overview; use the larger collections for specific biological
questions; always report the database version (these update annually
and results shift); be wary of any enrichment driven by fewer than
five overlapping genes — even when statistically significant, the
biology is thin.


== Summary <sec:summary>

- Differential expression is gene-by-gene hypothesis testing at
  twenty-thousand-test scale with three to six replicates per
  condition. Every step of the pipeline exists because naive
  approaches — t-tests per gene, Bonferroni correction across genes,
  no dispersion sharing — fail at this scale.
- The negative-binomial generalised linear model is the backbone of
  modern DE. DESeq2 and edgeR fit the same model, $log mu_(i g) =
  log s_i + sum_k beta_(k,g) x_(i,k)$, with NB likelihood and a
  per-gene dispersion $alpha_g$.
- Empirical-Bayes shrinkage of dispersion toward the
  mean–dispersion trend is why the model works on three replicates.
  Sharing information across genes turns a noisy per-gene estimate
  into a calibrated one; without shrinkage the analysis is
  uncalibrated and the p-values are unreliable.
- Benjamini–Hochberg controls the false-discovery rate rather than
  the family-wise error rate. Reporting `padj` from the BH procedure
  gives you a gene list with a controlled fraction of false
  positives — exactly what biological follow-up needs. Bonferroni is
  too conservative for DE work.
- The p-value histogram is the universal diagnostic. A flat null
  bulk with a near-zero spike means the test is well-calibrated and
  there is real signal. Any other shape signals model misspecification.
- Interpretation requires gene sets, not individual genes. A list of
  500 DE genes tells you little; an ORA or GSEA run against MSigDB
  tells you which molecular programmes are shifting. Use the
  enrichment layer as the reporting unit.


== Exercises <sec:exercises>

#strong[1.] #emph[Hand-fit a Poisson and an NB.]
Three samples have observed counts $5, 7, 11$ for the same gene.
Compute the maximum-likelihood Poisson rate $hat(lambda)$ and the
implied variance under the Poisson model. Compute the empirical
sample variance. By inspection, do these data look Poisson or
overdispersed? Give a back-of-envelope estimate of the NB dispersion
$alpha$ from $hat(sigma)^2 = mu + alpha mu^2$.

#strong[2.] #emph[Wald-test arithmetic by hand.]
A DE analysis reports for one gene $hat(beta)_1 = 0.62$ (log-fold
change in natural-log units) and $"SE"(hat(beta)_1) = 0.18$.
Compute the Wald $z$ statistic, convert to a two-sided p-value, and
convert $hat(beta)_1$ to a $log_2$ fold change. Would you call this
gene "significant" at a raw threshold of $p < 0.05$? Why is the raw
threshold the wrong question to ask in a twenty-thousand-gene
analysis?

#strong[3.] #emph[Implement Benjamini–Hochberg.]
Given the sorted p-values $0.001, 0.008, 0.04, 0.05, 0.20, 0.50$
from six tests, compute the BH-adjusted p-values at $alpha = 0.05$.
How many hypotheses are rejected? Show the BH line $(k/m) dot alpha$
at each rank and identify which is the largest $k$ at which the
sorted p-value lies below.

#strong[4.] #emph[Read a p-value histogram.]
Sketch the p-value histograms you would expect under (a) a healthy
DE analysis with strong signal, (b) a healthy analysis with no
signal, (c) an analysis whose model misspecification (an unmodeled
#idx("batch effect")batch effect) generates excess small p-values, and (d) an analysis
whose model is over-confident and generates inflated p-values
clustered near 1. For each, describe in one sentence what an analyst
should do next.

#strong[5.] #emph[ORA contingency by hand.]
A DE analysis identifies 800 genes from 20,000 tested. A curated
gene set has 200 members, of which 35 appear in the DE list. Build
the 2-by-2 contingency table and compute the expected number of
set-members in the DE list under the null. Is 35 over-represented?
Estimate the Fisher's-exact p-value to one decimal order of
magnitude (you may use the normal approximation to the
hypergeometric distribution).

#strong[6.] #emph[Size-factor offset.]
Sample A has size factor $1.5$; sample B has size factor $0.6$. For
a gene whose true expression is identical in both samples, what is
the ratio of expected counts $E[Y_A] / E[Y_B]$ predicted by the NB
GLM? Why is this ratio absorbed by the offset rather than appearing
in $hat(beta)_1$? What would happen to $hat(beta)_1$ if you omitted
the offset entirely?

#strong[7.] #emph[Wald vs LRT for a multi-level factor.]
You are analysing a time-course experiment with four time points,
coded as three indicator variables relative to baseline. You want a
single p-value per gene for "any difference across time." Which
test do you use and why? How many p-values would the Wald approach
produce per gene before correction, and what additional step would
you need?

#strong[8.] #emph[(Open-ended.)]
Pick one DE tool whose paper you have read — DESeq2, edgeR, limma-voom,
sleuth, or another — and describe in one paragraph the single most
distinctive design choice the authors made. What problem does the
choice solve, and where in the pipeline does it sit?


== Further Reading <sec:further-reading>

- *Anders, S., and Huber, W.* (2010). "Differential Expression Analysis
  for Sequence Count Data." _Genome Biology_ 11: R106. The original
  DESeq paper, introducing the median-of-ratios normalisation and
  shared-dispersion estimation.
- *Robinson, M. D., McCarthy, D. J., and Smyth, G. K.* (2010).
  "edgeR: A Bioconductor Package for Differential Expression Analysis
  of Digital Gene Expression Data." _Bioinformatics_ 26: 139–140. The
  edgeR paper.
- *Love, M. I., Huber, W., and Anders, S.* (2014). "Moderated
  Estimation of Fold Change and Dispersion for RNA-seq Data with
  DESeq2." _Genome Biology_ 15: 550. The DESeq2 paper. Probably the
  most-cited bulk RNA-seq paper of the decade.
- *Law, C. W., Chen, Y., Shi, W., and Smyth, G. K.* (2014). "Voom:
  Precision Weights Unlock Linear Model Analysis Tools for RNA-seq
  Read Counts." _Genome Biology_ 15: R29. The limma-voom paper.
- *Benjamini, Y., and Hochberg, Y.* (1995). "Controlling the False
  Discovery Rate: A Practical and Powerful Approach to Multiple
  Testing." _Journal of the Royal Statistical Society Series B_ 57:
  289–300. The foundational FDR paper.
- *Subramanian, A., Tamayo, P., Mootha, V. K., et al.* (2005). "Gene
  Set Enrichment Analysis: A Knowledge-Based Approach for Interpreting
  Genome-Wide Expression Profiles." _PNAS_ 102: 15545–15550. The GSEA
  paper.
- *Wu, D., and Smyth, G. K.* (2012). "CAMERA: A Competitive Gene Set
  Test Accounting for Inter-Gene Correlation." _Nucleic Acids
  Research_ 40: e133. The fix for correlated genes in
  permutation-based enrichment.
- *Liberzon, A., Birger, C., Thorvaldsdóttir, H., et al.* (2015). "The
  Molecular Signatures Database (MSigDB) Hallmark Gene Set
  Collection." _Cell Systems_ 1: 417–425. The Hallmarks paper.
- *DESeq2 vignette.* `bioconductor.org/packages/release/bioc/vignettes/DESeq2`
  Worked end-to-end examples on real datasets, kept in sync with each
  release.
