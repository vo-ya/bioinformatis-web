#import "../theme/book-theme.typ": *

= Causal Inference and Mendelian Randomisation <ch:causal-inference>

#matters[
  Every GWAS hit in this book — every locus you have learned to call, fine-map,
  meta-analyse, and partition — establishes one thing: a statistical
  association between a stretch of DNA and a phenotype. None of it, by itself,
  tells you what happens if you intervene. The difference between $P(Y mid X)$
  and $P(Y mid op("do")(X))$ is the difference between an observation and a
  recommendation. It is the difference between "people with low LDL have less
  heart disease" and "lower this patient's LDL and their heart-disease risk
  falls." Drug-development programmes worth billions of dollars rise and fall
  on which of those two claims the data actually supports. Mendelian
  randomisation is the move that lets us bridge the gap using nothing more
  than the GWAS summary statistics you have already been working with — and
  the laws Gregor Mendel published in 1866.
]

The previous twenty-four chapters have given you a long string of tools that
turn raw sequence into statistical statements: read alignment, variant calling,
expression quantification, GWAS, fine-mapping, polygenic scoring. By Chapter
13 you knew how to run a GWAS and read its Manhattan plot. By Chapter 14 you
had a pipeline that processed a biobank end-to-end. By this point in the
course a single SNP that crosses a $5 times 10^(-8)$ threshold is a routine
observation. The question this chapter takes seriously is what such an
observation _means_, and how — and when — you can turn it into a causal claim
without running a randomised trial.

The answer involves three pieces. The first is a vocabulary borrowed from
econometrics and statistics: instrumental variables, two-stage least squares,
and Judea Pearl's do-calculus. The second is a specific trick — using
germline SNPs as instrumental variables for exposures of interest — that
George Davey Smith and Shah Ebrahim christened _Mendelian randomisation_ in
2003. The third is a disciplined sensitivity workflow — MR-Egger, MR-PRESSO,
weighted median, the Steiger directionality test — that gives you a fighting
chance of detecting it when the trick goes wrong.

The chapter walks the three in that order. We start by sharpening the
distinction between association and causation, then introduce the
instrumental-variable estimator and derive the Wald ratio from first
principles. The middle of the chapter is the MR sensitivity quartet and the
practical workflow that ties it together. The last two sections cover the
mediation problem and the application most likely to put MR in your CV
within a pharmaceutical R&D group: drug-target validation by cis-MR. The
chapter closes with a brief tour of the wider causal-inference landscape —
difference-in-differences, regression discontinuity, doubly robust
estimators, causal forests — so you know what else lives in the toolbox when
MR is the wrong tool for the job.

This is the chapter that turns a Chapter 13 association into a Chapter 26
intervention. Read it before you give a drug-development team a target
recommendation, or before you accept one.


== Association, Confounding, and the Do-Operator <sec:ch25-do>

Take a SNP $Z$ associated with a phenotype $Y$ at $p = 4 times 10^(-12)$ in a
GWAS of 500,000 Europeans. What does the association tell you? Strictly:
$P(Y = 1 mid Z = 1) != P(Y = 1 mid Z = 0)$. Knowing $Z$ shifts the
probability distribution of $Y$. It is a _conditional_ statement about the
joint distribution of $(Y, Z)$ as observed in the wild. It is not a
statement about what happens when somebody changes $Z$ deliberately. The
gap between those two statements is the entire problem this chapter exists
to address.

There are exactly three causal stories consistent with a positive observed
correlation between two variables $X$ and $Y$. The first is _direct
causation_: $X$ acts on $Y$, so intervening on $X$ moves $Y$. The second is
_reverse causation_: $Y$ acts on $X$. The two variables move together
because $X$ responds to $Y$, not the other way around. The third is
_confounding_: a third variable $C$ acts on both $X$ and $Y$, and the
correlation is an indirect consequence of their shared parent — the
back-door path $X arrow.l C arrow Y$. The data alone cannot distinguish
the three. A regression of $Y$ on $X$ produces the same slope under any of
them; the same R-squared, the same residuals, the same standard error.

#figure(
  image("../../diagrams/lecture-25/01-three-explanations.svg", width: 92%),
  caption: [
    Three causal-graph topologies consistent with the same observed
    correlation between $X$ and $Y$: direct causation, reverse causation,
    confounding. Observational data cannot distinguish them; an intervention
    or an instrument can.
  ],
) <fig:three>

For germline GWAS hits, reverse causation is easy to rule out by
time-ordering — your genotype was fixed at conception, well before any
adult disease was diagnosed. (Somatic mutations break this argument and we
will return to them in @sec:ch25-pitfalls.) But confounding remains. The
classic genetic confounder is _population stratification_: the SNP and the
phenotype may both vary across ancestry strata for entirely independent
reasons. Adjusting GWAS for the top ten principal components — the move
you met in Chapter 13 — addresses the bulk of this, but not all of it, and
not at the long tail of fine-scale stratification that biobank-era cohorts
are starting to reveal.

The cleanest formalism for the gap is _Judea Pearl's do-calculus_, developed
over the late 1980s and 1990s and synthesised in the book _Causality_
(2000, second edition 2009). Pearl distinguishes two operations on a
probability distribution. The first is _conditioning_, written $P(Y mid X
= x)$, which selects the sub-population of individuals whose $X$ happens
to equal $x$. The second is _intervention_, written $P(Y mid op("do")(X =
x))$, which describes the population you _would_ observe if you reached
into the system and set $X$ to $x$ for everybody, regardless of their
upstream causes. The two distributions agree if and only if the graph
contains no _open back-door path_ from $X$ to $Y$ — that is, no path that
starts with an arrow into $X$ and ends in $Y$. Confounding is exactly the
existence of such a path.

#figure(
  image("../figures/ch25/f1-do-vs-see.svg", width: 95%),
  caption: [
    Seeing versus doing. On the left, conditioning on $X$ leaves the
    back-door path through the confounder $C$ open and the observed slope
    mixes the causal effect with the confounded path. On the right, the
    intervention $op("do")(X)$ severs the arrow into $X$, the back-door
    closes, and the slope is the causal effect $beta$ alone.
  ],
) <fig:do>

#note[
  The "surgery" metaphor matters. $op("do")$ is not a conditioning
  operation — it is a graph operation. To compute $P(Y mid op("do")(X =
  x))$ you literally delete every arrow into $X$ in the causal graph,
  then read off the distribution of $Y$ under the modified graph. The
  resulting object is a different probability distribution, defined on a
  hypothetical world. RCTs implement the surgery in the wet lab by
  randomly assigning $X$; MR implements it by Mendel's randomisation
  upstream of every confounder.
]

For an EE-trained reader the cleanest way to think about the difference is
that the causal effect is a _partial derivative under intervention_. If
$Y = beta X + gamma C + epsilon$ and $X$ also depends on $C$, then the
observational slope $partial E[Y mid X] / partial X$ picks up both $beta$
and a nuisance term proportional to $gamma$ times the dependence of $C$
on $X$. The causal slope $partial E[Y mid op("do")(X)] / partial X$ is
just $beta$. They
agree exactly when the back-door is closed and disagree when it is open.
This is the same thing as the bias term in an omitted-variable regression,
written in the language of graphs rather than the language of matrices.

The gold standard for closing the back door is the _randomised controlled
trial_. By tossing a coin to assign $X$, the experimenter guarantees that
$X$ is independent of every pre-treatment variable $C$, regardless of
whether $C$ was measured. The back-door arrow $C arrow X$ is gone by
construction. RCTs are expensive (typical phase-III trials cost a hundred
million dollars or more), slow (years), often unethical (you cannot
randomise people to smoking), and sometimes infeasible (you cannot
randomise lifetime LDL exposure). The trick of Mendelian randomisation
is to find an intervention nature has already run for you.


== Mendel's Randomisation as Nature's RCT <sec:ch25-mr-idea>

The premise is one sentence. At meiosis, each parent's pair of alleles at a
locus is randomly assigned to each gamete with probability one-half. This is
Mendel's second law — independent assortment, published in _Versuche über
Pflanzenhybriden_ in 1866 and rediscovered in 1900. From the child's point
of view, which of the parent's two alleles arrived is essentially a coin
toss. Across a large population, two randomly drawn individuals who carry
different alleles at a SNP differ at that SNP for reasons that have nothing
to do with their adult lifestyles, their childhood environments, or their
later disease risks. Their genotype was assigned upstream of all of those
confounders, in utero, before any choice they would ever make.

That coin toss is the randomisation an RCT achieves with an actual coin.
The argument runs: if individuals with different genotypes at a SNP have
different rates of some adult outcome $Y$, and that SNP is known (or known
to be highly likely) to influence $Y$ only through a specific intermediate
$X$ — say, the protein product of the gene the SNP sits in — then the SNP
"randomises" individuals to different levels of $X$. We can use it as a
substitute for a manipulator who randomly assigns $X$. The slope of $Y$
against the SNP-induced variation in $X$ is the causal effect of $X$ on
$Y$.

The idea did not arrive in 2003. It is older than Davey Smith and Ebrahim's
landmark _International Journal of Epidemiology_ paper, which named the
method and made it canonical. Martijn Katan published a remarkably modern
version in _The Lancet_ in 1986, asking whether low serum cholesterol
_caused_ cancer or merely correlated with it; he reasoned about the
apolipoprotein E genotype as a randomly assigned natural experiment that
should reveal the difference. Katan's argument was correct and went largely
ignored for fifteen years. The 2003 paper imported the econometric
machinery of instrumental variables, gave the method a name, and tied it
back to the GWAS era that was about to break. By 2010 the method had
software (`MendelianRandomization` for R, then `TwoSampleMR`); by 2015 it
had its own database (MR-Base); by 2020 it was a mandatory checkpoint in
many pharmaceutical target-validation pipelines.

#note[
  The instrumental-variable estimator is older than Mendelian randomisation
  by a lifetime. Sewall Wright introduced it in 1928 in a paper on the
  price of corn, with his father Philip Wright as co-author; the
  econometrics canon traces the modern form to Arthur Goldberger and the
  Cowles Commission in the 1950s and 1960s. Mendelian randomisation is
  the genomics-flavoured application of a tool that had been sitting in
  the econometrics toolbox for seventy years.
]

What MR buys you is a way to interrogate causal effects without running a
trial. What it does not buy you is freedom from assumptions. Every IV
estimator rests on assumptions that are themselves uncheckable in the same
data that produced the estimate. The art of doing MR honestly is the art of
defending those assumptions transparently, and the next two sections walk
the assumptions and the diagnostics one at a time.


== Instrumental Variables and the Wald Ratio <sec:ch25-iv>

Formally, an _instrumental variable_ for the causal effect of $X$ on $Y$ is a
variable $Z$ satisfying three conditions:

1. *Relevance.* $Z$ predicts $X$ — there is a real first-stage effect
   $pi = partial X slash partial Z$ that is not zero.
2. *Independence (exchangeability).* $Z$ is statistically independent of every
   confounder $C$ between $X$ and $Y$.
3. *Exclusion restriction.* $Z$ affects $Y$ only through $X$ — there is no
   direct path $Z arrow Y$ that does not pass through $X$, and no indirect
   path that bypasses $X$.

If all three hold, then a clever piece of arithmetic recovers the causal
effect $beta$ of $X$ on $Y$ from two regressions you can run on observable
data. Suppose the structural model is

$ X = pi Z + e, quad Y = beta X + epsilon $

with $E[e epsilon] != 0$ (so $X$ and $Y$ share a confounder buried in their
error terms — that is exactly what makes the naive OLS regression of $Y$
on $X$ biased). Substitute the first equation into the second:

$ Y = beta (pi Z + e) + epsilon = (beta pi) Z + (beta e + epsilon) $

Call the coefficient on $Z$ in this _reduced-form_ regression $Gamma$. We
just showed $Gamma = beta pi$. Hence

$ beta = Gamma slash pi = (op("Cov")(Z, Y)) / (op("Cov")(Z, X)) $

which is the _Wald ratio_ estimator. Two OLS slopes, one division. The
denominator is the first-stage slope $hat(pi)$ (the effect of $Z$ on $X$);
the numerator is the reduced-form slope $hat(Gamma)$ (the effect of $Z$ on
$Y$). Both are estimable from data, and the IV assumptions ensure that the
unobserved confounders, which contaminate the OLS of $Y$ on $X$, do
_not_ contaminate either of these two regressions: by independence and
exclusion, $Z$ is uncorrelated with the error in the structural equation
for $Y$.

#figure(
  image("../figures/ch25/f2-wald-derivation.svg", width: 95%),
  caption: [
    The Wald ratio in three regressions. The first stage gives $hat(pi)$ —
    the slope of $X$ on $Z$. The reduced form gives $hat(Gamma)$ — the
    slope of $Y$ on $Z$. The causal slope is their ratio. The substitution
    box on the right shows why $Gamma = beta pi$ under the structural
    model.
  ],
) <fig:wald-derivation>

#figure(
  image("../../diagrams/lecture-25/07-wald-ratio.svg", width: 90%),
  caption: [
    Per-SNP Wald ratios. Each SNP plotted at its $hat(beta)_("Z" arrow "X")$
    versus $hat(beta)_("Z" arrow "Y")$; the slope from the origin to each
    point is that SNP's per-instrument estimate of the causal effect.
    Consistency of those slopes across SNPs is the visual signature of
    valid IV behaviour.
  ],
) <fig:wald>

When you have many independent IVs, you can generalise the Wald ratio in
several ways. The simplest is the _inverse-variance-weighted_ (IVW) MR
estimator, which is exactly the weighted mean of the per-SNP Wald ratios:

$ hat(beta)_("IVW") = (sum_k w_k hat(beta)_("Z"_k arrow "Y") slash hat(beta)_("Z"_k arrow "X")) / (sum_k w_k), quad w_k = 1 / op("SE")^2 (hat(beta)_("Z"_k arrow "Y")) $

If every SNP is a valid instrument and the per-SNP errors are independent,
IVW is the minimum-variance unbiased estimator. Under the hood, IVW is
weighted least-squares regression of $hat(beta)_("Z"_k arrow "Y")$ on
$hat(beta)_("Z"_k arrow "X")$ through the origin, with weights equal to the
inverse variance of the outcome estimates.

The other classical generalisation is _two-stage least squares_ (2SLS), the
estimator econometricians have used for half a century. 2SLS replaces the
two univariate regressions with two multivariate ones. Stage 1 regresses
$X$ on the full matrix of instruments $bold(Z)$ (plus any control
covariates $bold(W)$ — ancestry PCs, sex, age, batch indicators) and
forms the fitted value $hat(X)$:

$ X = bold(Z) bold(pi) + bold(W) bold(mu) + e $

Stage 2 then regresses $Y$ on $hat(X)$ (still with the same controls):

$ Y = beta hat(X) + bold(W) bold(gamma) + u $

The coefficient $hat(beta)$ from stage 2 is the IV estimate. The intuition
is unchanged: stage 1 projects $X$ onto the column space spanned by $bold(Z)$
and $bold(W)$, isolating the part of $X$ that is "instrument-driven" and,
by the IV assumptions, uncorrelated with the unobserved confounders. Stage
2 regresses $Y$ on this purified component of $X$ and recovers an unbiased
estimate of the causal effect.

#figure(
  image("../../diagrams/lecture-25/02-2sls.svg", width: 92%),
  caption: [
    Two-stage least squares as projection followed by OLS. Stage 1 projects
    $X$ onto the column space of $bold(Z)$ to produce $hat(X)$. Stage 2
    regresses $Y$ on $hat(X)$. The projection is the debiasing step.
  ],
) <fig:tsls>

#note[
  Stage 1 of 2SLS is exactly a linear projection in the signal-processing
  sense. The fitted $hat(X)$ is the orthogonal projection of $X$ onto the
  subspace spanned by the instrument matrix; the residual $X - hat(X)$ is
  the part of $X$ orthogonal to that subspace and, under the IV
  assumptions, the part contaminated by confounders. Stage 2 keeps the
  projection and discards the residual. If you ever forget the algebra
  of 2SLS, the diagram is "instrument-driven component first, then OLS
  on the cleaned signal."
]

A useful variant is _two-sample MR_. Because the Wald ratio is a ratio of
two slopes, you do not need both slopes from the same individuals. You can
take $hat(beta)_("Z" arrow "X")$ from a published GWAS of the exposure (the
GIANT consortium for BMI, the GLGC for lipids, the IBSEN consortium for
inflammatory markers) and $hat(beta)_("Z" arrow "Y")$ from a separate
published GWAS of the outcome (CARDIoGRAMplusC4D for coronary disease, the
PGC for psychiatric phenotypes, FinnGen for hundreds of others). As long as
the two GWAS were performed on non-overlapping samples drawn from the same
underlying population, the ratio of slopes is consistent for the causal
effect. This recipe drove the explosion of MR in the late 2010s. You do not
need IRB access to individual-level data; you need a few thousand lines of R
and an internet connection.


== The Three IV Assumptions, One at a Time <sec:ch25-assumptions>

Stating "IV assumptions are satisfied" is the easy part. Defending each
assumption in a specific application takes most of the actual work. We walk
the three.

*Relevance* is the easiest to test. The diagnostic is the first-stage
$F$-statistic, the standard regression $F$ for the joint significance of the
instruments in the stage-1 regression. With a single instrument the
$F$-statistic equals $(hat(pi) slash op("SE")(hat(pi)))^2$, the squared
$t$-statistic. The rule of thumb that Staiger and Stock (1997) popularised is
$F > 10$; below that, the IV estimator is biased — in the one-sample case,
toward the OLS estimate that we were trying to avoid; in the two-sample
case, toward the null. Modern best practice tightens the threshold to $F >
20$ or even $F > 50$ for two-sample MR with many weak instruments. When
your candidate SNPs reach genome-wide significance ($p < 5 times 10^(-8)$)
in a well-powered GWAS, $F$ is usually comfortably above 30; weak-instrument
bias is rarely the binding constraint.

*Independence* asks that $Z$ shares no common cause with the confounders that
link $X$ and $Y$. For germline SNPs in a single, well-defined ancestry,
Mendel's randomisation makes this assumption very nearly automatic — the
genotype was assigned at conception, before lifestyle, before disease,
before measurement. The classical violation is population stratification:
allele frequencies vary across ancestries for evolutionary reasons that have
nothing to do with the phenotype, and those same ancestries differ in
disease prevalence for environmental and cultural reasons. Adjusting the
underlying GWAS for the top ten principal components closes most of this,
and within-ancestry MR is the default. Two more subtle violations are
_assortative mating_ — partners non-randomly choose each other on the
phenotype, inducing genotype-phenotype correlation at the population level
beyond what Mendel predicts — and _cryptic relatedness_, which biobank-scale
cohorts have shown to be more pervasive than expected (the UK Biobank
contains tens of thousands of close-cousin pairs). Both are real, both have
been quantified in recent literature, and both shift estimates by a
percentage point or two in expectations that depend on the trait.

*Exclusion restriction* is the assumption that does the work and the
assumption that fails. It says that $Z$ acts on $Y$ only through $X$ — no
direct path, no parallel path. For an MR study using a SNP $Z$ as an
instrument for a protein $X$, exclusion fails any time the SNP influences
$Y$ via a route that does not pass through $X$. The biological name for
that route is _horizontal pleiotropy_. There are two kinds. _Vertical
pleiotropy_ — the SNP affects multiple downstream traits, but only after the
exposure of interest sits on the path — is fine: it does not bias MR. SNP
$arrow X arrow Y_1$ and SNP $arrow X arrow Y_2$ is vertical and harmless.
_Horizontal pleiotropy_ — the SNP affects $Y$ through a pathway parallel to
$X$ — is exactly the violation. SNP $arrow X arrow Y$ together with SNP
$arrow T arrow Y$, where $T$ is some other trait, biases the Wald ratio.

#figure(
  image("../figures/ch25/f4-iv-assumptions-checklist.svg", width: 95%),
  caption: [
    The three IV assumptions, what each says, the diagnostic that
    interrogates it, and the typical Mendelian-randomisation failure mode
    for each. Exclusion is the only one without a clean within-data test —
    sensitivity analyses do most of the work.
  ],
) <fig:iv-check>

#figure(
  image("../../diagrams/lecture-25/11-pleiotropy.svg", width: 88%),
  caption: [
    Vertical (allowed) versus horizontal (forbidden) pleiotropy. Only
    horizontal pleiotropy invalidates an instrumental variable: the SNP
    must not have a path to $Y$ that bypasses $X$.
  ],
) <fig:pleio>

There is no within-data test that proves the exclusion restriction. There
are only sensitivity analyses that interrogate it from different angles. The
next section is about those analyses.

#warn[
  A common error in early MR papers was to declare exclusion "plausible"
  on biological grounds and stop there. Plausibility is necessary but not
  sufficient. Modern best practice is to report the full sensitivity
  quartet — IVW, MR-Egger, weighted median, MR-PRESSO — plus the Steiger
  directionality test, even when the biology is clean. Reviewers expect it;
  reproducibility benefits from it; the cost is one extra figure.
]


== The Sensitivity Quartet <sec:ch25-sensitivity>

If the IVW estimator were the only tool, MR would have a thin
methodological skeleton: one estimator, three assumptions, plenty of room
for invalid instruments to hide. The thickening of the skeleton over the
2015–2018 window was the introduction of three estimators that survive
specific assumption violations, plus a directionality test. Together they
form what working MR analysts call the _sensitivity quartet_. The point of
the quartet is not that any single estimator is universally best — none
of them is. The point is that the four make _different_ assumptions, and
convergence across the four is much stronger evidence than agreement with
any single one of them.

*MR-Egger* (Bowden, Davey Smith and Burgess, _IJE_ 2015) is the IVW
regression with an intercept term left free:

$ hat(beta)_("Z"_k arrow "Y") = alpha + beta dot hat(beta)_("Z"_k arrow "X") + epsilon_k $

The slope $hat(beta)$ remains an estimate of the causal effect; the
intercept $hat(alpha)$ estimates the average horizontal pleiotropy
across SNPs. If all instruments have a small, _systematic_ direct effect on
$Y$ — the so-called "InSIDE" violation, _Instrument Strength Independent of
Direct Effect_ — IVW absorbs that systematic bias into its slope. MR-Egger
absorbs it into the intercept, leaving the slope clean. A statistically
significant non-zero $hat(alpha)$ is a signal that systematic pleiotropy
is present; in that case, trust the Egger slope over the IVW slope. A
non-significant $hat(alpha)$ does not _prove_ exclusion holds, but it does
make systematic pleiotropy unlikely.

The price MR-Egger pays for the extra parameter is variance. Egger
estimates are typically two to three times wider than IVW estimates for the
same data, and the method has its own assumption (the InSIDE assumption)
that the SNP-exposure effects are independent of the per-SNP horizontal
pleiotropies. When InSIDE fails — for instance, when a single biological
pathway feeds both $X$ and $Y$ in correlated proportions across many SNPs —
Egger is biased too, just biased differently from IVW.

#figure(
  image("../../diagrams/lecture-25/08-mr-egger.svg", width: 92%),
  caption: [
    MR-Egger as IVW with a free intercept. A significantly non-zero
    intercept is the test for systematic horizontal pleiotropy. When the
    intercept is significant, the Egger slope is the corrected causal
    estimate.
  ],
) <fig:egger>

*MR-PRESSO* (Verbanck, Chen, Neale, and Do, _Nature Genetics_ 2018) takes a
different cut at the same problem. Instead of accommodating systematic
pleiotropy via an intercept, MR-PRESSO assumes the bulk of the SNPs are
valid and looks for individual outliers that cannot be reconciled with the
consensus. The procedure is iterative: fit IVW, compute a per-SNP residual
sum of squares, test each SNP's residual for significance under a global
null, drop any SNPs that fail at a chosen threshold, and re-fit IVW on the
survivors. Reported alongside the corrected estimate are a "global test" for
the presence of any outliers and a "distortion test" for whether outlier
removal changed the estimate significantly.

MR-PRESSO is sensitive to a handful of bad SNPs that bias IVW heavily but
do not show up as a systematic intercept in MR-Egger (because they are not
average — they are extreme). It is also model-free in the sense that it does
not require the InSIDE assumption.

#figure(
  image("../../diagrams/lecture-25/09-mr-presso.svg", width: 92%),
  caption: [
    MR-PRESSO outlier detection. Per-SNP Wald estimates are tested against
    the consensus; SNPs whose residuals are statistically extreme are
    typically pleiotropic and are removed before re-estimation.
  ],
) <fig:presso>

*Weighted median* (Bowden, Davey Smith, Haycock, and Burgess, _Genetic
Epidemiology_ 2016) sidesteps the assumption that _any specific subset_ of
instruments is valid. Instead, it requires only that at least half of the
total instrument weight come from valid SNPs. Under that "majority valid"
condition, the weighted median of the per-SNP Wald ratios is consistent for
the causal effect. The estimator is robust because the median is robust:
arbitrarily many outliers can be present, in either direction, as long as
they do not collectively command more than half the weight.

The two cousins are the _simple median_ (equally weighted; needs $> 50%$
valid SNPs by count) and the _modal estimator_ (consistent if a plurality
of SNPs share a single underlying causal effect — useful when valid
instruments are an _aligned_ minority).

The quartet's fourth member is the _Steiger directionality test_ (Hemani
et al., _PLOS Genetics_ 2017). Causation has a direction; in genomics, the
SNP must work via $X$ on its way to $Y$, not via $Y$ on its way to $X$.
Reverse causation is unusual for germline SNPs — the genotype is fixed at
birth — but the equivalent failure mode does occur when the candidate
exposure $X$ is itself a downstream consequence of the outcome $Y$. The
Steiger test compares, per SNP, how much variance the SNP explains in $X$
versus how much it explains in $Y$:

$ R^2_(Z, X) "vs" R^2_(Z, Y) $

If the SNP explains substantially more variance in $X$ than in $Y$, it
behaves like an instrument for the forward direction $X arrow Y$. If the
reverse, the directionality is suspect. The test is asymmetric and
disagreeing instruments are typical for traits that are biologically
intertwined (the canonical example is BMI and depression — does carrying
extra weight cause low mood, or does low mood cause weight gain?). Steiger
flags the controversy without resolving it; the resolution requires a
mechanistic argument outside the data.

#figure(
  image("../../diagrams/lecture-25/10-steiger.svg", width: 90%),
  caption: [
    Steiger directionality plot. Each SNP plotted at $R^2_(Z, X)$ versus
    $R^2_(Z, Y)$; SNPs above the diagonal support $X arrow Y$, SNPs below
    support $Y arrow X$. The fraction above the diagonal sets the
    directional confidence in the inferred causal direction.
  ],
) <fig:steiger>

A full sensitivity report on a modern MR analysis runs the four estimators
side by side and asks the reader to inspect the agreement. The book example
in @fig:quartet has been built so that three of the four estimators
converge on a true causal effect of $beta = 0.5$, while IVW is dragged
upward to 0.62 by a handful of pleiotropic outliers. Egger absorbs the
systematic bias into its intercept, weighted median rides above the outliers,
MR-PRESSO drops them and re-fits. The disagreement between IVW and the
other three is a signal that something is off; the agreement between the
other three is a signal that the corrected estimate is trustworthy.

#figure(
  image("../figures/ch25/f3-sensitivity-quartet.svg", width: 95%),
  caption: [
    Four estimators on a single instrument set. IVW is biased upward by
    three pleiotropic outliers; MR-Egger absorbs the systematic component
    into its intercept; weighted median is robust to the outliers by
    construction; MR-PRESSO flags and drops them. Convergence across the
    three corrected estimators is the convincing evidence.
  ],
) <fig:quartet>

#tip[
  The honest reporting convention, adopted by most working MR groups, is
  to publish IVW as the primary estimate, the other three estimators as
  sensitivity analyses, the MR-Egger intercept ($p$-value and magnitude),
  the MR-PRESSO global and distortion tests, the Steiger directionality
  fraction, and a heterogeneity test (Cochran's $Q$). If the primary
  estimate and the sensitivity estimates agree to within their confidence
  intervals, report the IVW. If they disagree, report the picture and
  explain.
]


== A Two-Sample MR in Practice <sec:ch25-twosample>

The full workflow for a two-sample MR study is short enough to write down
in one go and tedious enough that running it for the first time is a useful
exercise. The canonical example is _LDL cholesterol on coronary artery
disease_, the study that turned MR from a niche tool into a target-validation
standard.

#figure(
  image("../../diagrams/lecture-25/03-mr-flow.svg", width: 95%),
  caption: [
    The two-sample Mendelian-randomisation pipeline. Each stage has a
    canonical tool (TwoSampleMR, MR-Base, PLINK clumping) and a reportable
    diagnostic.
  ],
) <fig:flow>

The recipe has eight steps.

1. *Pick an exposure–outcome pair.* For our running example, exposure is
   LDL cholesterol; outcome is coronary artery disease (CAD).

2. *Pull GWAS summary statistics for the exposure.* For LDL, the Global
   Lipids Genetics Consortium (GLGC; Willer et al., _Nature Genetics_ 2013;
   updated in Klarin et al., _Nature_ 2018; Graham et al., _Nature_ 2021).
   Each row gives a SNP, an effect allele, a non-effect allele, an estimated
   effect size $hat(beta)_("Z" arrow "X")$, and a standard error.

3. *Select instruments.* Take the SNPs that reach genome-wide significance
   ($p < 5 times 10^(-8)$) for the exposure. Clump them so that no two
   selected SNPs are in strong linkage disequilibrium ($r^2 < 0.001$ within
   10 Mb is the standard). For LDL, this typically yields 80 to 200
   independent instruments depending on the GWAS vintage.

4. *Pull GWAS summary statistics for the outcome.* For CAD, the
   CARDIoGRAMplusC4D consortium (Nikpay et al., _Nature Genetics_ 2015;
   updated van der Harst and Verweij, _Circulation Research_ 2018).

5. *Harmonise.* Make sure the effect alleles match across the two studies.
   This is the step that produces silent bugs: a SNP coded with respect to
   the major allele in the exposure GWAS and the minor allele in the
   outcome GWAS has effects of opposite sign for spurious reasons. The
   `TwoSampleMR` package's `harmonise_data` function handles this for you
   if you let it; if you write your own pipeline, double-check it.

6. *Compute per-SNP Wald ratios* and run the IVW estimator. Compute MR-Egger,
   weighted median, and MR-PRESSO in parallel.

7. *Run Steiger.* Confirm that, for the bulk of instruments, the SNPs
   explain more variance in LDL than in CAD — i.e., the SNPs behave as
   LDL-instruments, not as CAD-instruments masquerading as LDL.

8. *Report.* For LDL → CAD, the modern result is roughly $hat(beta)_("IVW") =
  0.55$ in standardised units (95% CI 0.45–0.65), with an Egger intercept
   indistinguishable from zero, weighted median and MR-PRESSO agreeing
   within a few percent. Steiger directionality runs cleanly in the LDL →
   CAD direction. The evidence is about as clean as MR ever gets.

This is the study that backed the development of PCSK9 inhibitors and the
post-hoc validation of statins — and that we will return to in
@sec:ch25-cis-mr from the cis-MR angle.

#note[
  In EE language, the IVW estimator is weighted least-squares regression of
  $hat(beta)_("Z" arrow "Y")$ on $hat(beta)_("Z" arrow "X")$ through the
  origin. MR-Egger is the same regression with a free intercept. Weighted
  median is a robust regression (LAD-style) with the same weights.
  MR-PRESSO is robust regression with explicit outlier removal. All four
  are textbook statistical estimators applied to a particular structural
  setting — the value-add is the structural setting, not the statistics.
]


== Mediation Analysis <sec:ch25-mediation>

Establishing that $X$ causes $Y$ leaves an open question: through what
pathway? For drug development, the answer matters. If LDL causes CAD entirely
through atherosclerotic plaque buildup, then any LDL-lowering drug should
help, regardless of mechanism. If most of the effect runs through a
specific inflammatory pathway, then the right target may be downstream of
LDL, not LDL itself. _Mediation analysis_ decomposes the total causal effect
of $X$ on $Y$ into a _direct_ piece that bypasses any specified mediator and
an _indirect_ piece that flows through the mediator.

The classical formulation is Baron and Kenny's 1986 paper in the _Journal of
Personality and Social Psychology_, which proposed a three-regression test:
verify that $X$ predicts a candidate mediator $M$; verify that $X$ predicts
$Y$; then check that $X$'s coefficient attenuates when $M$ is added to the
regression of $Y$. If all three hold, $M$ "mediates" $X$ on $Y$. The Baron-
Kenny approach is widely taught in psychology and epidemiology and widely
known to be fragile. It assumes no unmeasured confounders of either the
$X$–$M$ relationship or the $M$–$Y$ relationship; under realistic conditions
it produces biased estimates of both direct and indirect effects, and the
"attenuation" criterion has no formal statistical meaning.

The modern alternative is the _counterfactual mediation framework_ of Robins,
Greenland, Pearl, and VanderWeele, codified in VanderWeele's _Explanation in
Causal Inference_ (Oxford, 2015). It decomposes the total effect into:

- *Natural Direct Effect* (NDE): the effect of $X$ on $Y$ that would remain
  if $M$ were held at its baseline value while $X$ moves.
- *Natural Indirect Effect* (NIE): the effect of $X$ on $Y$ that would flow
  if $X$ were held at its exposed value while $M$ moves from its baseline
  to its exposed value.

The two add to the total effect: $"total" = "NDE" + "NIE"$. The framework
makes the identifying assumptions explicit — no unmeasured $X$–$Y$
confounders, no unmeasured $M$–$Y$ confounders, no unmeasured $X$–$M$
confounders, and a fourth "no exposure-induced mediator-outcome confounder"
condition — and gives identifiability formulas that produce unbiased
estimates whenever those conditions hold.

#figure(
  image("../../diagrams/lecture-25/04-mediation.svg", width: 95%),
  caption: [
    Mediation decomposition. Total effect of $X$ on $Y$ splits into a
    natural direct effect (NDE) that bypasses the mediator $M$ and a
    natural indirect effect (NIE) that flows through it. Two-step MR
    estimates each leg with separate instruments.
  ],
) <fig:mediation>

What MR adds to mediation analysis is a way of estimating the two legs of
the decomposition without relying on the no-unmeasured-confounder
assumptions that Baron-Kenny needs. _Two-step MR_ (Burgess, Daniel,
Butterworth, and Thompson, _IJE_ 2015) does exactly that:

1. Estimate the causal effect of $X$ on $M$ via MR, using SNPs for $X$ as
   instruments. Call the estimate $hat(beta)_(X arrow M)$.
2. Estimate the causal effect of $M$ on $Y$ via MR, using SNPs for $M$ as
   instruments. Call it $hat(beta)_(M arrow Y)$.
3. The estimated indirect effect is the product:
   $"NIE" approx hat(beta)_(X arrow M) dot hat(beta)_(M arrow Y)$.
4. The estimated direct effect is the difference between the total effect
   (estimated separately by MR of $X$ on $Y$) and the indirect effect.

The product-of-coefficients formula is exact under linearity and no
exposure-induced confounding; for binary outcomes or strong non-linearities,
the better formula involves natural effects on the additive scale.

The worked BMI-on-CAD example from the lecture: MR for BMI → LDL gives a
positive but modest effect; MR for LDL → CAD gives the strong effect we
just established; the product is roughly $0.3$ times the total MR estimate
of BMI → CAD. The conclusion: about 30% of BMI's CAD risk runs through LDL;
about 70% runs through other pathways (probably inflammation, insulin
resistance, blood pressure, body fat distribution). Targeting LDL alone will
not neutralise BMI's effect on CAD; a broader intervention is required.

#tip[
  Mediation analysis is the right tool when you have a specific mediator
  in mind and a clean MR for both the $X$–$M$ and $M$–$Y$ legs. It is the
  wrong tool when you have a vague suspicion that "something downstream"
  matters — without naming the mediator and instrumenting it, two-step MR
  has no traction.
]


== Pitfalls and Limits of MR <sec:ch25-pitfalls>

MR is not a substitute for an RCT and the differences matter. We walk the
main ones.

*Weak instruments.* The first-stage $F$-statistic captures whether the
instruments actually predict the exposure. In the LDL → CAD example, $F$ is
typically in the thousands; weak-instrument bias is negligible. For
exposures where the per-SNP heritability is small and few SNPs reach
genome-wide significance — many quantitative traits, most behavioural
phenotypes, many psychiatric outcomes — $F$ may sit close to 10. The bias
is then non-trivial. The standard mitigation is to enlarge the GWAS until
genome-wide significant SNPs proliferate; failing that, to use a
polygenic-score-based MR with explicit weak-instrument-robust standard
errors (Burgess and Thompson, _Genetic Epidemiology_ 2017).

*Pleiotropy.* The biggest threat to MR validity, and the one the sensitivity
quartet was designed to interrogate. The HLA region is notoriously
pleiotropic — almost any SNP in it affects many immune and inflammatory
traits. Excluding the HLA region from instrument sets is standard. Several
other regions deserve similar caution: the $italic("FTO")$ locus
(adiposity, but also dozens of other traits via $italic("IRX3")$ and
$italic("IRX5")$); the $italic("APOE")$ locus (LDL, Alzheimer's disease,
longevity); the $italic("ABO")$ locus (blood type, but also CVD,
inflammation, infection risk). Modern MR pipelines either exclude these
loci preemptively or test results with and without them.

*Population stratification.* MR within a single ancestry, with the
underlying GWAS PCA-adjusted, is the default. Cross-ancestry MR is hard
because allele frequencies differ across ancestries for evolutionary
reasons unrelated to the exposure — running an MR with European
SNP-exposure betas and East-Asian SNP-outcome betas can produce serious
bias. The cleanest cross-ancestry analyses run within-ancestry MR
separately and meta-analyse only the final effect estimates.

*Selection bias and survival.* GWAS cohorts are alive at the time of
recruitment. SNPs that lower survival enough to deplete carriers from the
recruited cohort produce a distorted picture: people who survived to age
60 with the harmful allele are systematically different from those who
did not. For severely disease-associated SNPs, the recruited sample is
non-representative and MR estimates are biased — the formal name is
_collider bias_, and it is a real concern for any MR analysis of a
condition with strong selective mortality. The mitigation is to avoid MR
for traits where genotype-conditional survival is heavily skewed, or to
correct via inverse-probability-of-survival weights when those can be
estimated.

*Lifetime versus point intervention.* This is the limit MR readers most
often miss. An MR estimate is the average effect of a _lifetime_
difference in exposure between genetic groups. PCSK9 carriers spent 50
years with slightly different LDL than non-carriers. An MR slope estimates
$partial Y slash partial X$ averaged across that 50-year window. A drug that
lowers LDL for 2 years may produce a much smaller effect; conversely, a
drug that produces an effect 50% of MR's lifetime estimate may be the
right answer for clinical practice. MR validates the _direction_ and
_rough magnitude_ of the causal effect at lifetime exposure. It does not
substitute for an RCT's dose-response and temporal-window data.

*Reverse causation and dynamic feedback.* For germline SNPs against adult
disease, reverse causation is ruled out by time-ordering. But many MR
applications involve _intermediate_ exposures — molecular phenotypes that
both cause and are caused by downstream physiology. Inflammation causes
heart disease; heart disease causes inflammation. The Steiger test
catches the simplest cases but cannot resolve the full dynamic.

*Canalisation and developmental compensation.* A subtle objection: lifetime
genetic exposure may be partially compensated for by developmental
adaptation in ways that pharmacological intervention cannot mimic. The
classical example is a knockout that causes a lethal phenotype during
embryonic development but whose biochemical role can be safely targeted in
adults. MR for the gene's adult effects would systematically underestimate
the drug's potential because the cohort it studies has already developed
around the deficit. The argument has been made primarily by Davey Smith
himself and is the strongest case for treating MR as supportive rather than
definitive evidence.


== Drug-Target Validation by cis-MR <sec:ch25-cis-mr>

The most consequential application of MR, by financial impact, is _drug-
target validation by cis-MR_. The idea is narrow and clever. Suppose you are
considering a drug that inhibits the protein product of gene $G$ to lower
exposure $X$ and reduce risk of disease $Y$. Before committing to the
hundreds of millions of dollars an RCT will cost, you would like to know
whether reducing $G$'s product would lower disease risk in real humans.
General MR for $X arrow Y$ — using whatever genome-wide SNPs predict $X$ —
answers a question that is too broad: it tells you whether _any_ pathway
lowering $X$ would lower $Y$, including pathways that no drug can hit.

Cis-MR sharpens the question. Instead of using all genome-wide significant
SNPs for $X$ as instruments, use only the SNPs in or near $G$ itself —
specifically, the cis-eQTLs or pQTLs within roughly $plus.minus 100$ kb
of the transcription start site of $G$. These SNPs alter the expression or
activity of $G$'s product and are biologically local: by the very locality
of their genomic position, they are far less likely to be horizontally
pleiotropic for other pathways. A cis-MR estimate using these SNPs is a
"natural drug experiment" for the gene of interest — what happens to
disease risk when the protein's expression or activity is altered by
nature's genetic variation.

#figure(
  image("../../diagrams/lecture-25/05-drug-target-mr.svg", width: 95%),
  caption: [
    The drug-target MR pipeline. Cis-eQTLs or cis-pQTLs near the candidate
    gene serve as natural drug experiments; a PheWAS scan with the same
    instruments flags potential off-target effects before any RCT runs.
  ],
) <fig:drugmr>

The canonical example is _PCSK9_. In 2003 Marianne Abifadel and colleagues
(Boileau lab, Paris) identified PCSK9 gain-of-function mutations in a
French family with familial hypercholesterolaemia. Three years later, Helen
Hobbs's group at UT Southwestern found loss-of-function PCSK9 mutations in
African Americans associated with dramatically lower LDL and dramatically
lower coronary disease risk (Cohen et al., _NEJM_ 2006). Both findings
constituted, in retrospect, a small MR study: people with naturally lower
PCSK9 activity had less coronary disease. The pharmaceutical industry took
the cue. Amgen's evolocumab and Sanofi/Regeneron's alirocumab — both
monoclonal antibodies that bind PCSK9 — moved through clinical trials in
the early 2010s. The FOURIER trial (Sabatine et al., _NEJM_ 2017) and the
ODYSSEY trial (Schwartz et al., _NEJM_ 2018) confirmed that pharmacological
PCSK9 inhibition reduced LDL by roughly 60% and CAD events by roughly 15%
to 20% over 2.2-year median follow-up. Both drugs were approved in 2015 at
list prices that have since fallen substantially; combined annual revenue
exceeds 1 billion dollars. The MR prediction preceded the trials by a
decade.

#figure(
  image("../figures/ch25/f5-pcsk9-cis-mr.svg", width: 95%),
  caption: [
    PCSK9 cis-MR walk-through. The locus zoom selects cis-eQTL SNPs within
    $plus.minus 100$ kb of the PCSK9 gene as instruments. The MR scatter
    shows clean LDL → CAD signal with negligible Egger intercept. The
    phenome-wide scan finds no off-target adverse effects across sixteen
    diverse phenotypes.
  ],
) <fig:pcsk9>

#figure(
  image("../../diagrams/lecture-25/12-cis-mr.svg", width: 92%),
  caption: [
    Cis-MR locus view for PCSK9. The eQTL signal in the gene's cis-window
    is the basis for the instrument set; the per-SNP Wald ratios cluster
    cleanly; the PheWAS panel scans for off-target effects.
  ],
) <fig:cismr>

The pipeline runs in five steps.

1. *Pick the target gene.* PCSK9 in our running example; IL-6R, HMGCR,
   CETP, IL-1$beta$, TYK2, GLP-1R in others.
2. *Find cis-instruments.* Take SNPs within $plus.minus 100$ kb of the
   gene's transcription start site (sometimes $plus.minus 1$ Mb for genes
   in gene-poor regions) that are significant eQTLs in a relevant tissue
   (GTEx for general tissue, eQTLGen for blood, the GTEx single-cell atlases
   for cell-type-specific signal). Clump to LD-independence ($r^2 < 0.1$
   is common for cis-instruments because the cis-region is small and
   strict clumping leaves too few SNPs).
3. *Run MR for the on-target effect.* IVW + Egger + weighted median +
   MR-PRESSO using the cis-instruments and the disease outcome GWAS. For
   PCSK9, the on-target effect is clean: LDL lowering, CAD reduction.
4. *Run a phenome-wide scan.* Apply the same cis-instruments to GWAS for
   every available phenotype — typically the GWAS catalogue, the UK
   Biobank PheWAS panel, FinnGen, BBJ. The output is a forest plot of
   estimated effects of the target on every measurable trait. Off-target
   effects show up as significant hits in phenotypes the drug is _not_
   meant to treat — diabetes, cancer, autoimmune disease, kidney function,
   psychiatric outcomes. For PCSK9, the PheWAS is clean; for the IL-6R
   tocilizumab work, the PheWAS predicted both the cardiovascular benefit
   and a modest infection-risk signal that subsequent trials confirmed.
5. *Decide.* If the on-target effect is large enough and the off-target
   signals are tolerable, the drug-development decision is reinforced. If
   the on-target effect is null, the programme should pause. If the
   off-target signals are severe, the programme should pause regardless of
   the on-target effect.

The economics of this matter. A pre-trial cis-MR study costs roughly
hundred thousand dollars in compute and analyst time. A failed phase-III
trial costs hundreds of millions and depletes a company's pipeline of
opportunity cost on the order of a billion. Pharmaceutical companies that
have institutionalised cis-MR as a target-validation gate — Genentech,
GSK, Regeneron, AstraZeneca's Open Targets initiative — report that
roughly one in three target proposals that look promising on cell-line
data fail the cis-MR check and are deprioritised before any trial. That
hit rate is industry-relevant.

#tip[
  When evaluating a published cis-MR study or running one yourself, three
  questions are worth asking. First, are the cis-instruments truly cis —
  within the $plus.minus 100$ kb window — or has the analyst slipped in
  trans-acting variants masquerading as cis-instruments? Second, has the
  PheWAS been run at sufficient phenotype breadth? A clean on-target
  signal with a narrow PheWAS is not a clean signal at all. Third, do the
  cis-instruments work in the right tissue? An eQTL in liver tells you
  about the protein in liver; a drug that fails to reach liver may not
  reproduce the effect.
]

The cis-MR move generalises beyond classical drug development. Combined
with CRISPR functional screens (Chapter 26), single-cell perturbation
atlases (Chapter 27), and downstream pathway modelling, cis-MR forms the
backbone of the Open Targets platform — a public consortium that scores
every gene–disease pair by every available causal-evidence stream and
publishes a target-prioritisation score. By 2024, Open Targets covered
more than 60,000 gene–disease pairs and had become a routine starting
point for pharmaceutical target proposals across the industry.


== Beyond Mendelian Randomisation <sec:ch25-beyond-mr>

MR is one corner of the causal-inference landscape. Several other tools
deserve recognition because they fill niches MR cannot, and because the
working bioinformatician will eventually need at least one of them.

#figure(
  image("../../diagrams/lecture-25/06-toolkit.svg", width: 95%),
  caption: [
    The causal-inference toolkit at a glance. Each technique trades a
    different set of assumptions against the data structure it requires;
    the appropriate choice depends on whether you have an instrument, a
    sharp threshold, a policy change, or only observational covariates.
  ],
) <fig:toolkit>

*Difference-in-differences* compares the change in outcome between an
exposed group and an unexposed group across a policy shift or natural
experiment. Useful for epidemiological questions where a policy creates a
clean before-versus-after contrast — the 2010 New York sodium reduction
order, the 2017 Australian sugar-tax pilot. Limited use in genomics but
genuinely useful where it applies; under parallel-trends assumptions, it
removes time-invariant confounders by differencing them out.

*Regression discontinuity* exploits sharp assignment thresholds. If a
patient is treated when a test score crosses a cutoff and untreated when
it does not, comparing the two sides of the threshold within a narrow
band gives a quasi-randomised estimate of the treatment effect. The
genomics analogue is polygenic-score-based screening — if patients above
the 95th-percentile PRS receive an intervention and those below do not,
the just-above versus just-below comparison estimates the intervention's
local average causal effect.

*Pearl's full do-calculus* and the related _back-door / front-door
criteria_ provide a formal framework for deciding, given an assumed
causal graph, which conditional independencies must be invoked to
identify a particular causal effect. The toolkit is mostly used to
diagnose whether a proposed analysis can recover the effect of interest
under the assumed graph. It is rarely used as the estimator itself;
rather, it is the bookkeeping that tells you which estimator is valid.

*Doubly robust estimators* — augmented inverse-probability weighting,
TMLE (targeted maximum likelihood estimation) — combine an outcome model
with an exposure model in a way that is consistent if _either_ model is
correctly specified. They have become standard in observational
pharmacoepidemiology and are increasingly seen in biobank-scale
genetic-epidemiology pipelines.

*Causal forests* (Athey and Wager, _JASA_ 2019) extend the random-forest
framework to estimate _heterogeneous_ treatment effects — how the causal
effect varies across individuals as a function of covariates. The natural
genomics application is precision-medicine analytics: does a drug work
better in patients with one polygenic profile than another? The literature
is still young but growing.

*MR-CAUSE* (Morrison, Knoblauch, Marcus, Mukherjee, _Nature Genetics_
2020) is a Bayesian MR estimator that explicitly models correlated
horizontal pleiotropy — a violation that the InSIDE assumption of MR-Egger
forbids and that MR-PRESSO does not capture. For complex traits with
intertwined biology, MR-CAUSE and its successors are increasingly the
honest choice when standard MR cannot defend its assumptions.

The frontier moves quickly. _Multi-trait MR_ estimates the joint effect of
several correlated exposures on a single outcome, addressing the
network-pleiotropy problem head-on. _Non-linear MR_ relaxes the linearity
assumption that the Wald ratio implicitly imposes. _Bidirectional MR_
formalises the Steiger directionality test into an estimator that returns
both directions and a posterior over which dominates. _Within-family MR_
uses siblings as natural controls for population structure and
assortative mating, at the cost of much smaller effective sample size.
Each is a useful refinement; none is a wholesale replacement.

#warn[
  The genuine frontier issue, as of 2024, is the gap between MR and the
  rapidly growing single-cell atlases. Most public eQTLs are bulk-tissue;
  most disease mechanisms are cell-type-specific. Cell-type-specific
  cis-MR is methodologically straightforward but data-limited: the
  cell-type-specific GWAS sample sizes are too small for genome-wide
  significance at the per-cell-type level for most genes. Single-cell
  consortia like the Human Cell Atlas, the eQTL Catalogue, and recent
  preprints on single-cell-resolved cis-eQTLs are filling the gap, but
  the methodological work and the data still have several years of run
  ahead of them.
]


== Summary <sec:ch25-summary>

- Association is not causation. Three causal stories — direct, reverse,
  confounded — produce the same correlation in observational data. The
  distinction matters because translation into therapy depends on which
  one is true.
- Pearl's do-calculus separates the observational distribution
  $P(Y mid X)$ from the interventional distribution $P(Y mid op("do")(X))$.
  RCTs implement the intervention by random assignment. Mendelian
  randomisation borrows the intervention from meiosis.
- An instrumental variable for $X arrow Y$ is any variable $Z$ that (i)
  predicts $X$, (ii) is independent of confounders, and (iii) affects $Y$
  only through $X$. The Wald ratio $hat(beta)_("Z" arrow "Y") slash
  hat(beta)_("Z" arrow "X")$ recovers the causal effect under the three
  assumptions.
- The inverse-variance-weighted estimator is the multi-instrument
  generalisation of the Wald ratio. Two-stage least squares is the
  same idea in matrix form, with stage 1 a linear projection onto the
  instrument subspace and stage 2 an OLS on the projected exposure.
- The exclusion restriction is the assumption that fails. The sensitivity
  quartet — MR-Egger for systematic pleiotropy, MR-PRESSO for outlier
  pleiotropy, weighted median for majority-validity, Steiger for
  directionality — interrogates it from four complementary angles. Honest
  MR reports all four.
- Two-sample MR combines GWAS summary statistics from non-overlapping
  cohorts via the same Wald ratio. It is enormously efficient: hundreds
  of MR analyses can be run with no individual-level data access.
- Mediation analysis decomposes a total effect into direct and indirect
  components. Two-step MR estimates each leg with separate instruments
  and combines them as a product, robust to confounders that would
  invalidate Baron-Kenny mediation.
- Cis-MR uses cis-eQTLs near a candidate drug-target gene as instruments
  for the target's expression or activity. The PCSK9 cis-MR study
  preceded evolocumab and alirocumab by a decade and is the canonical
  template for pharmaceutical target validation.
- MR is not an RCT. It estimates the lifetime average effect of natural
  exposure, not the effect of a 2-year pharmacological intervention.
  Selection bias, weak instruments, population stratification, and
  developmental canalisation are the standard limits.
- Other causal-inference tools — difference-in-differences, regression
  discontinuity, doubly robust estimators, causal forests, MR-CAUSE —
  fill the niches MR cannot. The right tool depends on the data
  structure and which assumptions you can defend.


== Exercises <sec:ch25-exercises>

#strong[1.] #emph[ASCII to causal graph.] For each of the following
biological scenarios, draw the causal graph (boxes for variables, arrows
for directed influences). Identify any open back-door paths from $X$ to
$Y$. Indicate whether the naive OLS regression of $Y$ on $X$ would
overestimate, underestimate, or correctly estimate the causal effect, and
why. (a) $X$ = serum LDL; $Y$ = coronary artery disease; unmeasured
confounder = adult diet. (b) $X$ = vitamin-D level; $Y$ = mortality;
unmeasured confounder = sunlight exposure (which also relates to physical
activity). (c) $X$ = body-mass index; $Y$ = depression; bidirectional
causation suspected.

#strong[2.] #emph[Wald ratio by hand.] You have three SNPs as instruments
for a single exposure on a single outcome. The per-SNP coefficients
$(hat(beta)_("Z" arrow "X"), hat(beta)_("Z" arrow "Y"))$ with standard
errors are: SNP1: $(0.20 plus.minus 0.02, 0.10 plus.minus 0.03)$; SNP2:
$(0.30 plus.minus 0.03, 0.16 plus.minus 0.04)$; SNP3: $(0.15 plus.minus
0.02, 0.07 plus.minus 0.02)$. Compute the three per-SNP Wald ratios.
Compute the IVW estimate using inverse-variance-of-outcome weights.
Compare to a simple unweighted mean of the three ratios.

#strong[3.] #emph[Implement two-sample MR.] In Python, load the public
TwoSampleMR LDL → CAD example dataset (`TwoSampleMR::TwoSampleMR_data` in
R, or the matched files at MR-Base). Implement IVW, MR-Egger, and
weighted median MR estimators from scratch. Compare your estimates with
those produced by the `TwoSampleMR` R package on the same data. Report
agreement to three significant figures and explain any discrepancy.

#strong[4.] #emph[Weak-instrument bias.] Simulate two-sample MR with one
instrument. The true causal effect is $beta = 0.5$. The first-stage
$F$-statistic is varied across simulations from 5 to 100. Plot the bias
and standard error of the IVW estimate as functions of $F$. Confirm the
rule of thumb $F > 10$ is too lenient for two-sample settings when the
instrument is weak. Briefly explain why two-sample MR bias goes toward
the null while one-sample bias goes toward OLS.

#strong[5.] #emph[Egger intercept interpretation.] You run MR on 50
instruments for an exposure $X$ and outcome $Y$. IVW gives $hat(beta) = 0.40$
(SE 0.05); MR-Egger gives $hat(beta) = 0.18$ (SE 0.12) with intercept
$hat(alpha) = 0.02$ (SE 0.005, $p = 6 times 10^(-5)$). Weighted median
gives $hat(beta) = 0.20$ (SE 0.06). What do you conclude about the
exclusion restriction? What is the most defensible causal effect estimate
to report? What additional analyses would you run before publishing?

#strong[6.] #emph[Mediation product.] Two-step MR for BMI → LDL gives
$hat(beta)_(X arrow M) = 0.18$. Two-step MR for LDL → CAD gives
$hat(beta)_(M arrow Y) = 0.55$. Direct MR for BMI → CAD gives
$hat(beta)_(X arrow Y) = 0.30$. What fraction of BMI's effect on CAD runs
through LDL? What does the residual direct effect (the part not running
through LDL) plausibly correspond to biologically?

#strong[7.] #emph[Cis-MR design.] You are a target-validation analyst at
a pharmaceutical company. Marketing has proposed inhibiting interleukin-6
receptor (IL-6R) for atherosclerosis. Design the cis-MR study. Specify:
(a) the cis-window for instrument selection; (b) the tissues whose eQTLs
you would use; (c) the primary outcome GWAS; (d) the panel of PheWAS
phenotypes for safety; (e) what evidence would advance the programme to a
phase-I trial, and what would kill it. Cite at least one published study
on IL-6R MR.

#strong[8.] #emph[(Open-ended.)] Pick one published MR study from
2018–2024 that reports a controversial or contested causal effect (for
example: vitamin D on COVID-19 mortality; alcohol on dementia; coffee on
cardiovascular disease). Audit the study against the sensitivity-quartet
discipline of @sec:ch25-sensitivity. Identify the weakest assumption. Cite
one follow-up paper that supports the original conclusion and one that
challenges it. Write a one-paragraph assessment of how confident a
clinician should be in the headline number.


== Further Reading <sec:ch25-further-reading>

- *Davey Smith, G., and Ebrahim, S.* (2003). "'Mendelian randomization':
  can genetic epidemiology contribute to understanding environmental
  determinants of disease?" #emph[International Journal of Epidemiology]
  32: 1–22. The paper that named the field and is still the cleanest
  introduction.
- *Bowden, J., Davey Smith, G., and Burgess, S.* (2015). "Mendelian
  randomization with invalid instruments: effect estimation and bias
  detection through Egger regression." #emph[International Journal of
  Epidemiology] 44: 512–525. The MR-Egger paper; pair with Bowden,
  Davey Smith, Haycock and Burgess (2016) for the weighted median.
- *Verbanck, M., Chen, C. Y., Neale, B., and Do, R.* (2018). "Detection
  of widespread horizontal pleiotropy in causal relationships inferred
  from Mendelian randomization between complex traits and diseases."
  #emph[Nature Genetics] 50: 693–698. MR-PRESSO. Their supplementary
  walks through the simulation code.
- *Hemani, G., et al.* (2018). "The MR-Base platform supports systematic
  causal inference across the human phenome." #emph[eLife] 7: e34408.
  MR-Base and the `TwoSampleMR` package — start here for any new
  analysis.
- *Pearl, J.* (2009). #emph[Causality: Models, Reasoning, and Inference],
  2nd edition. Cambridge University Press. The book that built the
  do-calculus. Read chapters 3 and 11 first.
- *VanderWeele, T. J.* (2015). #emph[Explanation in Causal Inference:
  Methods for Mediation and Interaction]. Oxford University Press. The
  modern reference for mediation analysis in the counterfactual
  framework.
- *Holmes, M. V., Ala-Korpela, M., and Davey Smith, G.* (2017).
  "Mendelian randomization in cardiometabolic disease: challenges in
  evaluating causality." #emph[Nature Reviews Cardiology] 14: 577–590.
  A clinically grounded survey of where MR has and has not delivered.
- *Sanderson, E., Davey Smith, G., Windmeijer, F., and Bowden, J.* (2019).
  "An examination of multivariable Mendelian randomization in the
  single-sample and two-sample summary data settings."
  #emph[International Journal of Epidemiology] 48: 713–727. Multivariable
  MR — the natural next step once a single exposure is no longer enough.
- *Open Targets Platform.* https://platform.opentargets.org. The public
  consortium consolidating MR and other causal-evidence streams across
  every gene–disease pair.
