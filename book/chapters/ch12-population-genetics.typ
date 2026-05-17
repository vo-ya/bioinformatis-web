#import "../theme/book-theme.typ": *

= Population Genetics: Allele Frequencies and Their Fates <ch:population-genetics>

#matters[
  Every analysis in this book up to now has treated a single genome at a
  time — one FASTQ, one BAM, one VCF. Population genetics asks a
  different question. Given a collection of genomes from a population
  — a thousand individuals, a hundred thousand individuals, all of
  Iceland — what can you say about the _distribution_ of variation? The
  answer matters for GWAS power, for clinical-allele frequency cutoffs,
  for evolutionary inference, and for the parts of human history that
  the historical record never touched. The mathematical objects in this
  chapter — Hardy-Weinberg ratios, the Wright-Fisher Markov chain, the
  Kingman coalescent — are also the cleanest stochastic-process
  isomorphisms in the entire course. Read them as both biology and
  signal processing.
]

The shift from one-genome to many-genome analysis is the shift from
descriptive bioinformatics to inferential bioinformatics. A single
variant call in a VCF tells you what the patient has. A collection of
variant calls across ten thousand patients tells you what the species
has — and, with enough cleverness, what it had a hundred thousand years
ago. Both are useful. The first lets you diagnose; the second lets you
predict, and to do model-based reasoning under uncertainty about the
generative process that produced the data.

This chapter walks the path from raw allele counts to that
generative-process reasoning. Section 12.1 lays out the basic
vocabulary — allele frequencies, genotype frequencies, Hardy-Weinberg
equilibrium — and explains why the centrepiece of the field is a model
that nobody believes is literally true. Section 12.2 introduces the
Wright-Fisher Markov chain, the simplest forward-in-time model of
allele-frequency dynamics, and the four forces that perturb it:
drift, selection, mutation, migration. Section 12.3 picks up linkage
disequilibrium — the genome-scale autocorrelation structure that GWAS
exploits and that selection scans use to date sweeps. Section 12.4
flips time around and develops the coalescent: the same Wright-Fisher
process run backward from a present-day sample, which is dramatically
cheaper to simulate. Section 12.5 turns those tools onto real data
with PSMC, SMC++, and ADMIXTURE — methods that infer effective
population size and ancestry composition from one or many sequenced
genomes. Section 12.6 covers selection scans: detecting the genomic
footprints of recent adaptation with controlled false-alarm rates.

The chapter assumes you have a multi-sample VCF or its equivalent in
hand. Everything downstream of this — GWAS in Chapter 13, polygenic
prediction, evolutionary genomics, conservation genetics — uses the
machinery this chapter develops.


== Allele Frequencies and Hardy-Weinberg <sec:hwe>

A *biallelic SNP* in a population has two alleles, say $A$ and $a$.
Pool the genomes of $N$ diploid individuals and you have $2N$
chromosomes total. Count the copies of each allele; call the counts
$k_A$ and $k_a$. The *allele frequency* of $A$ is
$ p = k_A / (2N), quad q = 1 - p. $
This single number $p$ is the basic currency of population genetics.
Every more sophisticated quantity — genotype frequency, heterozygosity,
effective population size, $F$-statistics, LD measures — is built on
allele-frequency tallies.

Three operational conventions are worth pinning down before going
further. Minor allele frequency (MAF) is $min(p, q)$ and is the
quantity most filtering pipelines use; the convention "MAF $> 5%$ = common variant,
$1text("–")5%$ = low-frequency, $< 1%$ = rare" recurs everywhere and
splits real datasets into regimes that different tools handle
differently. Genotype frequencies — the fraction of individuals with
each diploid genotype $A A$, $A a$, $a a$ — are not the same object as
allele frequencies, and the link between the two is the topic of the
next subsection. And ploidy matters: the X and Y chromosomes and
the mitochondrion break the diploid assumption and have their own
conventions; pipelines treat them as separate categories.

=== The Hardy-Weinberg Principle

In 1908, two people independently published the same result. *G. H.
Hardy*, a mathematician at Cambridge, wrote a half-page letter to
_Science_; *Wilhelm Weinberg*, a physician in Stuttgart, presented a
slightly longer paper in German. Both were responding to a practical
confusion in the early-Mendelian literature: an objection that
dominant alleles should become more common over generations simply
because they are dominant. Hardy and Weinberg showed, in a couple of
algebraic lines, that this is wrong. Under a set of idealised
conditions, allele frequencies do not change across generations, and
genotype frequencies are a deterministic function of a single allele
frequency:
$ P(A A) = p^2, quad P(A a) = 2 p q, quad P(a a) = q^2. $

The result is the oldest equation in population genetics and the most
heavily used. It is the *null model* — the distribution expected when
none of the evolutionary forces operate. Five idealising conditions
have to hold for HWE to be exact:

1. *Random mating* (no preference for particular genotypes).
2. *No selection* (every genotype has equal fitness).
3. *No mutation* (alleles don't switch identity).
4. *No migration* (no gene flow from other populations).
5. *Infinite population size* (no drift).

None of these conditions hold in any real population. So why is HWE
the centrepiece of the field?

Because the interesting question is never "is this locus in HWE?" The
question is whether the locus deviates from HWE _more than expected by
sampling noise_, and if so, in what direction. Excess homozygotes
point at population structure or inbreeding. Excess heterozygotes
point at overdominance. A specific direction of deviation, repeated
across loci, points at admixture. HWE is the noise floor against
which biology is detected.

#figure(
  image("../../diagrams/lecture-12/02-hardy-weinberg.svg", width: 85%),
  caption: [
    The Hardy-Weinberg parabolas. Genotype frequencies are a
    deterministic function of a single allele frequency — the
    heterozygote class peaks at $p = 0.5$ and vanishes at fixation
    ($p = 0$ or $p = 1$).
  ],
) <fig:hwe-parabolas>

#note[
  HWE is the stationary distribution of a Markov chain on diploid
  genotype states under idealised assumptions. In signal-processing
  terms, it is what a memoryless channel outputs under symmetric
  input: given allele frequencies $(p, q)$, the genotype distribution
  in the next generation is determined and does not depend on past
  history. Deviations from HWE are memory effects — population
  structure (memory in space), inbreeding (memory in mating
  structure), recent selection (memory in fitness history). Testing
  for HWE is testing whether the allele frequency is a sufficient
  statistic.
]

=== Testing the Null

The canonical HWE test is a chi-squared goodness-of-fit comparison.
Given observed genotype counts $n_(A A)$, $n_(A a)$, $n_(a a)$, the
expected counts under HWE are $N p^2$, $2 N p q$, $N q^2$, where $p$
is estimated from the sample as
$p = (2 n_(A A) + n_(A a)) / (2 N)$. The test statistic
$ chi^2 = sum (O_i - E_i)^2 / E_i $
follows a $chi^2_1$ distribution under the null (the one degree of
freedom reflects estimating $p$ from the data). For large genome-wide
panels, exact tests (Wigginton et al. 2005) are preferred at low
counts; tools like `vcftools --hardy` and `plink --hwe` compute these
per SNP across millions of variants in one pass.

#figure(
  image("../figures/ch12/f1-hwe-test-cascade.svg", width: 95%),
  caption: [
    A worked Hardy-Weinberg goodness-of-fit test on a 200-individual
    cohort. The observed het count is below the expected; the
    chi-squared statistic of 8.69 rejects HWE at $p approx 0.003$.
    What the rejection _means_ depends on the direction and the
    cohort-wide pattern.
  ],
) <fig:hwe-test>

The interpretive panel of @fig:hwe-test is the part that matters in
practice. A single-locus rejection can mean genotyping error
(by far the most common cause in real data), local population
structure, real selection at the site, or a rare variant whose
sampling distribution is poorly approximated by the asymptotic
chi-squared. A pattern of rejection _in the same direction across
many loci_ is much more informative — that is what tells you about
admixture or cohort-wide structure.

#warn[
  None of the four HWE-violating forces ever truly switch off. In
  every real population, mutation, drift, selection, and migration
  are all operating at some level. The question is therefore not "is
  this population in equilibrium?" — answer: no — but "is the
  deviation at this locus or in this cohort larger than the sampling
  noise allows?" HWE testing is a screen, not a verdict. Be especially
  cautious before invoking biology from a single rejected locus —
  genotyping artefacts dominate at the strict-significance tail.
]

=== Diagnostic Deviations

Five canonical deviation patterns each map to a different biological
process. *Excess homozygosity* — too many $A A$ and $a a$, not enough
$A a$ — typically reflects either the *Wahlund effect* (two
subpopulations with different allele frequencies pooled as if they
were one) or *inbreeding* (relatives mating, offspring inheriting the
same allele from both parents). *Excess heterozygosity* points at
*overdominance* — the canonical example is the $H b S$ (sickle-cell)
allele, where heterozygotes resist falciparum malaria while $H b S\/H b S$
homozygotes suffer sickle-cell disease — or at negative assortative
mating, common in plant breeding systems with self-incompatibility.
*Genotype patterns correlated with environment* indicate population
stratification, which becomes a confounder in GWAS (Chapter 13).
*Allele-frequency shifts between generations* are the live signature
of selection, drift, or migration; long-term sampling (museum
specimens, ancient DNA) makes these directly observable. *Single-site
extreme rejections* are usually genotyping artefacts: a clustered
miscall pattern, a copy-number variant masquerading as a SNP, or a
mismapped read pile.

#tip[
  When a single SNP rejects HWE at $p < 10^(-6)$ in a large cohort,
  the priors say it is a genotyping artefact roughly nine times out of
  ten. Tools like Hail's `hwe_normalize` and the GATK Best-Practices
  filters use HWE rejection precisely _as_ a quality filter, not as a
  biology detector. The biology shows up when many loci near each
  other reject in the same direction.
]


== The Wright-Fisher Model and the Forces of Change <sec:wright-fisher>

If Hardy-Weinberg is the null, what is the simplest non-trivial
process that produces deviations from it? The answer, formalised
independently by R. A. Fisher (1930) and Sewall Wright (1931), is the
*Wright-Fisher model*: a discrete-generation, finite-population
Markov chain on allele counts. Every more elaborate generative model
in population genetics — selection, mutation, structured populations,
the coalescent — is a generalisation of this kernel.

The setup is minimal. A population has $N$ diploid individuals, so
$2 N$ chromosomes. Generations are discrete and non-overlapping. To
produce generation $t + 1$, sample $2 N$ chromosomes _with
replacement_ from generation $t$. There is no selection, no mutation,
no migration — just sampling. The allele count in the next generation
is then
$ X_(t + 1) | X_t tilde "Binomial"(2 N, p_t), $
where $p_t = X_t / (2 N)$. This defines a Markov chain on the state
space ${0, 1, 2, ..., 2 N}$. The states $X = 0$ and $X = 2 N$ are
absorbing: once an allele is lost or fixed, it stays that way. Given
enough generations, every finite Wright-Fisher population eventually
fixes one allele or the other.

#figure(
  image("../../diagrams/lecture-12/03-wright-fisher.svg", width: 92%),
  caption: [
    The Wright-Fisher generative scheme. Each chromosome in the next
    generation is an independent random draw with replacement from
    the current generation — allele frequency drifts with variance
    $p(1 - p) / (2 N)$ per step.
  ],
) <fig:wf-scheme>

#note[
  The Wright-Fisher model is a discrete-time, discrete-state Markov
  chain whose transition kernel is binomial. In stochastic-process
  terms it is a birth-death chain on ${0, 1, ..., 2 N}$ with drift
  (under selection) and a diffusion approximation (the Kimura 1955
  forward equation). All of stochastic-process theory applies:
  stationary distributions, mixing times, expected hitting times to
  absorbing states. The whole quantitative output of classical
  population genetics — expected heterozygosity decay, fixation
  probability and time, the allele-frequency spectrum — follows from
  this Markov chain plus its diffusion limit.
]

=== Genetic Drift

*Genetic drift* is the per-generation change in allele frequency
arising purely from finite-population sampling. The conditional mean
and variance are straightforward consequences of the binomial:
$ E[p_(t + 1) | p_t] = p_t, quad
  "Var"(p_(t + 1) | p_t) = p_t(1 - p_t) / (2 N). $
Drift is unbiased in expectation but has variance inversely
proportional to population size. Small populations drift fast and
lose variation quickly; large populations drift slowly.

Three quantitative consequences fall out of this kernel and are worth
internalising:

- *Expected heterozygosity decays* at rate $1 \/ (2 N_e)$ per
  generation. Starting from $H_0 = 2 p q$, we have
  $H_t = H_0 (1 - 1 \/ (2 N_e))^t approx H_0 exp(- t \/ (2 N_e))$.
  The half-life is roughly $1.39 N_e$ generations.
- *Probability of fixation* of a neutral allele equals its current
  frequency: $P("fix") = p$. A brand-new mutation (frequency
  $1 \/ (2 N_e)$) has that same probability of eventually fixing.
- *Expected fixation time* for an allele that does fix is approximately
  $4 N_e$ generations (starting from $p = 0.5$; longer starts from
  the boundary).

@fig:drift-math gathers the numerics into a single panel for three
representative population sizes. The lesson is uniform: every
population-genetics quantity that depends on drift scales with $N_e$.
This is why $N_e$ is, in practice, the single most consequential
parameter in the field.

#figure(
  image("../figures/ch12/f2-drift-variance-math.svg", width: 95%),
  caption: [
    Drift as binomial sampling noise at three effective population
    sizes. Per-generation step size, fixation probability, and
    heterozygosity half-life all scale with $N_e$. The bottom panel
    sketches trajectory shapes; the right panel collects the
    long-run consequences.
  ],
) <fig:drift-math>

#figure(
  image("../../diagrams/lecture-12/04-drift-trajectories.svg", width: 95%),
  caption: [
    Drift trajectories at two population sizes from the same starting
    frequency. The small-$N$ panel shows wild excursions and frequent
    fixation; the large-$N$ panel stays tightly clustered near $0.5$
    — variance scales as $1 \/ (2 N)$.
  ],
) <fig:drift-trajectories>

=== Selection

*Selection* introduces a deterministic bias on top of drift. Genotypes
differ in fitness — survival times reproductive success — and
higher-fitness genotypes contribute disproportionately to the next
generation. Parameterise the three diploid genotypes by relative
fitness $w_(A A) = 1$, $w_(A a) = 1 - h s$, $w_(a a) = 1 - s$, where
$s$ is the *selection coefficient* against $a$ and $h$ is the
*dominance coefficient* ($h = 0$ is recessive, $h = 0.5$ is additive,
$h = 1$ is dominant). The deterministic allele-frequency update is
$ p_(t + 1) = p_t (p_t w_(A A) + q_t w_(A a)) / overline(w), $
where $overline(w)$ is mean population fitness. Stack the WF binomial
sampling on top of this deterministic update and you have a Markov
chain that combines drift with directional bias.

The most consequential observation in selection-vs-drift theory is
that the answer is governed not by $s$ alone, nor by $N$ alone, but by
their product $2 N_e s$. Three regimes:

- *Strong selection* ($|2 N_e s| >> 1$): the deterministic term
  dominates. Allele frequencies move predictably toward whichever
  side fitness prefers.
- *Weak selection* ($|2 N_e s| << 1$): drift dominates. The allele
  behaves nearly as if neutral.
- *Moderate selection* ($|2 N_e s| approx 1$): both forces matter;
  stochastic dynamics with a bias.

A selection coefficient of $s = 10^(-3)$ is strongly effective in
$N_e = 10^6$ ($2 N_e s = 2000$) but essentially neutral in
$N_e = 100$ ($2 N_e s = 0.2$). This is why microbial populations
(huge $N_e$) feel even tiny advantages, while vertebrate populations
(modest $N_e$) tolerate substantial mildly deleterious variation. The
scaling explains why human genetic load — the burden of slightly
deleterious variants drifting around at moderate frequency — looks
nothing like _E.~coli_'s.

#figure(
  image("../../diagrams/lecture-12/05-selection.svg", width: 95%),
  caption: [
    Selection vs drift. The scaled quantity $2 N_e s$ — not $s$ alone
    — determines whether an allele behaves deterministically or
    near-neutrally. The same selection coefficient gives qualitatively
    different fates in different effective sizes.
  ],
) <fig:selection>

=== Mutation, Migration, and Effective Size

The two remaining forces in the WF kernel are easier to summarise.
*Mutation* introduces new alleles at per-base rate $mu approx 10^(-8)$
per generation in humans — about 70 new mutations per child genome.
Most are lost within a few generations to drift; the surviving rare
fraction is what builds the standing pool of variation.
*Migration* moves alleles between populations at per-generation rate
$m$. Even modest migration ($m approx 0.01$) homogenises allele
frequencies between populations on evolutionary timescales — a fact
that has shaped human population structure over the past tens of
thousands of years.

The deeper concept is *effective population size*, $N_e$. The census
size of a real population — the number of individuals you would count
walking through a forest — is almost never the number that controls
genetic-drift behaviour. The relationship $N_e <= N$ usually holds
strictly, and often $N_e << N$. Several mechanisms shrink it:

- *Unequal reproductive success.* If a few individuals produce most
  of the offspring, drift is dominated by that subset.
- *Skewed sex ratio.* $N_e = 4 N_f N_m / (N_f + N_m)$ for an
  organism with $N_f$ breeding females and $N_m$ breeding males.
- *Fluctuating population size.* Bottlenecks drop $N_e$ far more than
  a simple time-average would suggest — the harmonic mean dominates.
- *Overlapping generations* and *age structure* further reduce $N_e$
  compared to the simple Wright-Fisher idealisation.

For humans, $N approx 8 times 10^9$ today but $N_e approx 10^4$ from
genetic evidence — a four-orders-of-magnitude gap. The 10 000 reflects
the historical narrowness of our lineage: bottlenecks during the
out-of-Africa expansion, founder events in non-African populations,
small hunter-gatherer band sizes over most of our evolutionary past.
$N_e$ is the relevant parameter for drift, selection thresholds, and
essentially every population-genetic prediction — it is _the_ number
to reach for whenever the formula calls for "$N$."

#warn[
  Effective size $N_e$ and census size $N$ are different quantities,
  and the difference can be more than a thousand-fold. When a paper
  or a software default says "population size 10 000," check whether
  it means census, current effective, or historical effective — these
  are three different numbers and substituting one for the other
  produces qualitatively wrong predictions.
]


== Linkage Disequilibrium: Genome-Scale Autocorrelation <sec:ld>

The Wright-Fisher model in §12.2 tracks a single locus. Real genomes
have millions of loci on a few dozen chromosomes, and alleles at
nearby loci are *not independent*: they share ancestry, they
co-segregate during meiosis, and they reflect the recombinational
history of the population. *Linkage disequilibrium* (LD) is the name
for that non-independence, and it is the centre of gravity of modern
population genetics. GWAS works because LD lets a small set of "tag"
SNPs stand in for millions of nearby ones. Selection scans work
because the LD around a swept allele is the durable fingerprint of
the sweep.

=== Where LD Comes From and How It Goes Away

If alleles at two loci were independent, the haplotype frequency
$P(A B)$ would equal the product $p_A p_B$ of the marginal allele
frequencies. The deviation
$ D = P(A B) - p_A p_B $
is the *raw LD coefficient*. Three generative processes produce
positive $D$:

1. *Historical founding.* When a new mutation arises, it appears on
   one specific chromosome. Initially it is in perfect LD ($D$ at
   its maximum) with every flanking variant on that chromosome.
2. *Admixture.* When two populations with different allele
   frequencies mix, chromosomes carry population-specific haplotypes
   and the admixed population has LD between any pair of
   differing-frequency loci.
3. *Selection.* A positively selected allele drags its flanking
   variants up in frequency with it, generating a footprint of
   elevated LD around the selected site — the "selective sweep"
   pattern picked up by iHS and XP-EHH (§12.6).

Three destructive processes erode LD: *recombination* breaks
haplotypes apart at meiosis; *drift* perturbs the deterministic
co-occurrences stochastically; and *time* simply lets recombination
compound. The single most consequential prediction in LD theory is the
exponential decay of LD with recombination distance and time, derived
in @fig:ld-decay.

=== Quantifying LD: $r^2$ and $D'$

The raw $D$ depends on the marginal allele frequencies and so cannot
be compared across loci. Two normalisations dominate practice.

The *squared correlation* $r^2$ is
$ r^2 = D^2 / (p_A p_a p_B p_b). $
This is the squared Pearson correlation between the binary allele
indicators at the two loci. $r^2 in [0, 1]$. It is the statistic that
GWAS power depends on: the variance explained by a tag SNP for an
unobserved causal SNP is proportional to $r^2$ between them.

The *normalised D-prime* is
$ D' = D / D_max, $
where $D_max$ is the maximum value $|D|$ could take given the marginal
allele frequencies. $D' in [-1, +1]$. $D'$ captures structural LD —
which haplotypes exist — better than $r^2$, which also requires the
allele frequencies to be similar to reach high values. The two
statistics together carry complementary information.

#tip[
  In a typical genome-wide LD analysis, $D' = 1$ means "no
  recombination has separated these two loci since they last shared an
  ancestor," while $r^2 = 1$ is a much stronger claim: it also
  requires the two alleles to have the same frequency. You can have
  $D' = 1$ and $r^2 approx 0$ when one locus is rare and the other
  common. Real LD-block boundaries are most cleanly seen in $D'$
  matrices; tag-SNP design is driven by $r^2$.
]

=== LD Decay with Distance

In a randomly mating population without selection, LD decays
geometrically with recombination distance and time. Per generation,
a fraction $c$ of pairs gets broken apart by recombination (where $c$
is the recombination probability between the two loci, with
$c approx$ genetic distance in Morgans). The recursion is
$D_(t + 1) = (1 - c) D_t$, so
$ D_t = D_0 (1 - c)^t approx D_0 exp(-c t), $
and consequently
$ r^2(t) approx r^2(0) exp(-2 c t). $
The factor of two in the exponent appears because $r^2$ is the square
of $D / sqrt("variance")$. Under drift-recombination equilibrium, the
expected value is
$ E[r^2] approx 1 / (1 + 4 N_e c), $
which makes the LD-decay rate scale with the product $N_e c$ — bigger
effective size means slower decay for any given physical distance.

#figure(
  image("../figures/ch12/f3-ld-decay-derivation.svg", width: 95%),
  caption: [
    The LD-decay derivation, with worked numerics and three example
    curves. Younger populations carry strong LD over long distances;
    in modern human cohorts $r^2$ halves at roughly 30 kb.
  ],
) <fig:ld-decay>

@fig:ld-matrix shows the genome-wide consequence: an $r^2$ heatmap of
nearby SNPs reveals a *block-diagonal* structure. Stretches of high
LD ("LD blocks") are punctuated by narrow regions where LD drops
sharply ("recombination hotspots"). The 2002 International HapMap
Project, followed by the 1000 Genomes Project (2008–2015) and the
Human Genome Diversity Panel, mapped these blocks empirically; their
catalogue underpins every commercial GWAS array.

#figure(
  image("../../diagrams/lecture-12/06-ld-matrix.svg", width: 95%),
  caption: [
    A genomic LD matrix. The block-diagonal structure reveals LD
    blocks punctuated by recombination hotspots — the empirical basis
    for GWAS tag-SNP design.
  ],
) <fig:ld-matrix>

#figure(
  image("../../diagrams/lecture-12/07-ld-decay.svg", width: 95%),
  caption: [
    Average $r^2$ vs physical distance in a human cohort. The
    half-decay distance lies around 30 kb — the operational scale of
    tag-SNP coverage.
  ],
) <fig:ld-decay-empirical>

#note[
  Linkage disequilibrium is the autocorrelation structure of a 1D
  stochastic process along the chromosome. Each SNP is a binary
  random variable at a fixed coordinate; the genome is the sample
  path. LD between two SNPs is their covariance; $r^2$ is their
  squared Pearson correlation. Genomic LD is the autocorrelation of
  a _non-stationary_ process — recombination rate varies along the
  genome, so the autocorrelation length depends on position. Every
  autocorrelation-based signal-processing tool generalises:
  windowed estimators, spectral methods, matched filtering for
  sweep signatures. Tag-SNP GWAS is a sampling-rate argument — if you
  sample at the autocorrelation scale, the missing SNPs are
  approximately recoverable.
]


== The Coalescent: Wright-Fisher Run Backward <sec:coalescent>

The Wright-Fisher model in §12.2 runs forward in time. You start at
some historical generation, sample binomially, and watch allele
frequencies evolve. This is conceptually clean but operationally
wasteful: most of the simulated lineages leave no present-day
descendants, and the relevant signal — what an actual sample of
contemporary chromosomes looks like — is buried inside a much larger
simulation. *J.~F.~C. Kingman*, a Cambridge probabilist with no
biological agenda, published a four-page paper in 1982 in _Stochastic
Processes and their Applications_ that inverted the entire framework.

His insight was this. Take a sample of $n$ chromosomes from the
present. Look at the pedigree backward in time. At each step
backward, the lineages that contributed to the sample either stay
distinct or two of them happen to share a parent — a *coalescence*
event. The waiting time to the next coalescence depends only on the
number of remaining lineages, not on the rest of the population. The
"coalescent process" thereby simulates only the small set of
lineages that actually matter, ignoring the $2 N - n$ irrelevant ones.
For samples of $n approx 10^2$ from populations of
$N approx 10^4$, the speed-up is a hundredfold or more.

=== Waiting Times and TMRCA

The quantitative core of the Kingman coalescent is one equation. In a
Wright-Fisher population of size $2 N$, two randomly chosen lineages
have a per-generation probability of $1 / (2 N)$ of sharing a parent
in the previous generation. With $k$ lineages, the number of
distinct pairs is $binom(k, 2)$, and the probability that _any_ pair
coalesces in the previous generation is approximately $binom(k, 2) /
(2 N)$. The waiting time $T_k$ until the next coalescence (going
backward) is therefore approximately exponential:
$ T_k tilde "Exp"(binom(k, 2) / (2 N)), $
with expected value $E[T_k] = 2 N / binom(k, 2)$ generations. The
expected time to the most recent common ancestor (TMRCA) of the
entire sample is the sum:
$ E[T_("MRCA")] = sum_(k = 2)^n 2 N / binom(k, 2)
                = 4 N sum_(k = 2)^n 1 / (k(k - 1))
                = 4 N (1 - 1 / n). $
The expression has a striking property: as $n -> infinity$,
$E[T_("MRCA")] -> 4 N$. Going from 10 samples to 1000 samples barely
moves the TMRCA. The last coalescence — from two lineages down to one
— takes as long as the first $n - 2$ combined, because the rate
$binom(2, 2) / (2 N) = 1 / (2 N)$ is the smallest rate the process
ever encounters.

#figure(
  image("../figures/ch12/f4-coalescent-waiting-times.svg", width: 95%),
  caption: [
    Per-interval waiting times for a sample of ten lineages. The
    final interval (2 lineages collapsing to 1) is the longest by a
    wide margin and dominates the total tree height. The bottom
    panel shows TMRCA's plateau at $4 N_e$ as sample size grows.
  ],
) <fig:coal-waiting>

#figure(
  image("../../diagrams/lecture-12/08-coalescent-tree.svg", width: 90%),
  caption: [
    A coalescent tree for $n = 8$. Branch lengths shrink as more
    lineages remain; the last pairwise coalescence dominates total
    depth, which is why TMRCA plateaus near $4 N_e$.
  ],
) <fig:coal-tree>

=== From Trees to Summary Statistics

Given a coalescent tree, every observable summary statistic of
genetic variation in the sample has a tractable expectation. With
mutations dropped onto branches at rate $mu$ per site per generation
(a Poisson process), the *number of segregating sites* in a sample
of $n$ has expectation
$ E[S] = theta sum_(k = 1)^(n - 1) 1 / k, $
where $theta = 4 N_e mu$ is the population-scaled mutation rate. The
expected number of *pairwise differences* between two randomly
sampled chromosomes is simply $E[pi] = theta$.

The two estimators give the same population parameter $theta$ under
neutrality but weight branches differently — $pi$ weights internal
branches (where most heterozygosity sits), and $S / sum 1/k$ weights
total tree length. The difference $pi - S / sum 1/k$ is the
numerator of *Tajima's D* (Tajima 1989), the classical test for
departures from neutrality (§12.6).

#figure(
  image("../../diagrams/lecture-12/09-coalescent-wf-duality.svg", width: 95%),
  caption: [
    The same realised population, run forward as Wright-Fisher and
    backward as the coalescent. Both views describe identical
    statistics on the sample; the coalescent simulates only the
    sub-pedigree that actually reaches the present.
  ],
) <fig:wf-coal-duality>

#note[
  The coalescent is a continuous-time Markov chain on the number of
  active lineages, with exponential waiting times between transitions
  $k -> k - 1$ at rate $binom(k, 2) / (2 N_e)$. In signal-processing
  terms it is an inverse-time process — the analogue of running a
  linear time-invariant system backward to recover driving inputs
  from outputs, or Kalman-smoothing a state trajectory from the
  present into the past. The Wright-Fisher-to-coalescent duality is a
  specific instance of a deep principle: given a Markov forward
  process, the conditional distribution of the past given the
  present is another Markov process. Population geneticists
  discovered the duality independently and gave it their own name.
]

=== Recombination and the Ancestral Recombination Graph

The standard coalescent assumes a single genealogy for all sites in
the sample — the entire chromosome shares one tree. Recombination
violates this. A crossover in the history of a lineage means that
different segments of the chromosome have different ancestral paths.
The generalisation is the *ancestral recombination graph (ARG)*:
along the chromosome, the tree changes at recombination breakpoints,
with most adjacent trees sharing most of their structure but differing
in one or a few subtree topologies. Modern coalescent simulators
(`msprime`, `tskit`; Kelleher et al. 2016, 2022) represent ARGs as
*tree sequences* — a compact data structure where the differences
between adjacent trees are stored as a sequence of edits, allowing
genome-scale simulation of millions of samples in seconds.

ARG-aware inference — reconstructing the actual ancestral history from
present-day data — is a young and increasingly important field. Tools
like `Relate` (Speidel et al. 2019) and `tsinfer` plus `tsdate`
(Kelleher et al. 2019; Wohns et al. 2022) infer approximate ARGs from
modern variation panels, opening direct access to allele ages and
genome-wide genealogical structure. The ARG is, in a sense, the
ultimate summary statistic — every other statistic ($pi$, $S$,
Tajima's D, $F_("ST")$, LD) is a deterministic function of the ARG.


== Inferring Population History <sec:demography>

The machinery of §12.2–12.4 gives us a generative model: a
demographic history (sequence of population sizes, split times,
migration rates) plus a mutation rate produces an expected pattern of
variation in any sample. Inverting that map — reading the demographic
history off the observed data — is *demographic inference*. The methods
below are state-of-the-art in 2024 and likely to remain so for the
near future.

=== PSMC: Demographic History from One Genome

The most surprising result in modern demographic inference is that
you can recover an effective-size trajectory $N_e(t)$ over a million
years from a _single diploid genome_. Heng Li and Richard Durbin
published *PSMC* (Pairwise Sequentially Markovian Coalescent) in
_Nature_ in 2011 demonstrating this. The intuition is that a diploid
genome contains two chromosome copies — maternal and paternal —
which have their own coalescent time at every locus. Locally low
heterozygosity means a recent MRCA at that locus (so small $N_e$
around that time); locally high heterozygosity means an ancient
MRCA (so large $N_e$ then).

The PSMC algorithm models the genome as a hidden Markov model. The
hidden state is the local coalescent time between the two diploid
copies, discretised into time bins. The emission is the observed
heterozygosity in a window. The transition probabilities come from
the *sequentially Markovian coalescent* (SMC) approximation: as you
walk along the chromosome, the local TMRCA changes only at
recombination events, and the new TMRCA is correlated with the old.
Fitting the HMM by Baum-Welch (an EM variant) gives a step-function
$N_e(t)$ that explains the observed heterozygosity profile.

#figure(
  image("../../diagrams/lecture-12/10-psmc-ne.svg", width: 95%),
  caption: [
    PSMC $N_e(t)$ trajectories for three continental populations. A
    single diploid genome reconstructs roughly a million years of
    demographic history; the out-of-Africa bottleneck appears as a
    distinct dip in non-African curves around 50–100 kya.
  ],
) <fig:psmc>

PSMC has two well-known blind spots. It cannot resolve the very
recent past — the last few thousand years — because the two diploid
alleles in any one genome simply have not had time to differ much in
that window. And it cannot resolve the very ancient past beyond about
$10 N_e$ generations, because by then the two alleles have always
coalesced and the data carries no signal. *SMC++* (Terhorst, Kamm,
and Song 2017) extends PSMC by combining the within-individual
diploid coalescent with the site-frequency spectrum across many
individuals — restoring resolution in the recent past at the cost of
needing dozens or hundreds of samples. *Relate*'s ARG-based approach
(Speidel et al. 2019) extracts demographic history directly from
inferred genealogies, giving the richest output of all but at higher
computational cost.

=== Admixture and Population Structure

Human populations are not isolated lineages. Modern Africans, Europeans,
East Asians, and the others differentiate, exchange migrants, mix
historically, and admix recently. *Admixture inference* asks: given a
panel of genomes, what is the best decomposition into a small number
of ancestral source populations, and what fraction of each genome
comes from each source?

The dominant formal model is unsupervised mixture-model decomposition.
Assume each genome $i$ is generated as a mixture of $K$ ancestral
populations with proportions $Q_i = (Q_(i, 1), ..., Q_(i, K))$ summing
to 1. Each ancestral population $k$ has its own allele-frequency
vector $F_k = (F_(k, 1), ..., F_(k, M))$ over $M$ SNPs. The observed
allele dosage at SNP $j$ in individual $i$ is binomial with success
probability $sum_k Q_(i, k) F_(k, j)$. Inference recovers both $Q$
and $F$ simultaneously from the genotype matrix.

Three implementations cover the workflow space. *STRUCTURE*
(Pritchard, Stephens, and Donnelly 2000) was first — a Bayesian MCMC
sampler, slow but principled. *ADMIXTURE* (Alexander, Novembre, and
Lange 2009) replaced the MCMC with a quasi-Newton optimisation of the
likelihood, achieving 10–100$times$ speedup with comparable accuracy;
it became the standard tool of the GWAS era. *fastSTRUCTURE* (Raj
et al. 2014) uses variational inference for similar speed with a
slightly different bias profile. All three require the user to
specify $K$, the number of ancestral populations.

#figure(
  image("../figures/ch12/f5-admixture-mixture-model.svg", width: 95%),
  caption: [
    ADMIXTURE as EM on a mixture model. Each genome is a convex
    combination of $K$ source allele-frequency vectors; EM
    alternates posterior-responsibility computation with
    parameter updates. The $K$-selection curve shows why
    cross-validation, not training likelihood, is the right model-
    order criterion.
  ],
) <fig:admixture-model>

#figure(
  image("../../diagrams/lecture-12/11-admixture.svg", width: 95%),
  caption: [
    Stacked-bar admixture output for a multi-continental panel.
    Each column is one individual; the bar segments are the inferred
    ancestry proportions — a direct visualisation of the $Q$ matrix.
  ],
) <fig:admixture-bars>

#note[
  ADMIXTURE is unsupervised mixture-model decomposition — the
  genomics instance of Gaussian-mixture clustering, spectral mixture
  analysis in remote sensing, or independent-component analysis
  variants in signal processing. Each observation (genome) is a
  mixture of $K$ latent sources (ancestral populations); the task is
  to recover both the source "spectra" (allele frequencies per
  population) and the per-sample mixture coefficients (ancestry
  proportions). Standard machinery applies: EM for point estimates,
  variational inference for posteriors, BIC or cross-validation for
  model-order selection. Admixture in population genetics has the
  same mathematical structure as unmixing a multispectral image.
]

The *$K$-selection* problem is the standard model-order selection
problem dressed up in genomic vocabulary. Cross-validation error,
not training likelihood, is the right criterion — training likelihood
increases monotonically with $K$ and gives no useful answer.
Real cross-validation curves are often flat across several values of
$K$, especially in admixed cohorts where the "true number of source
populations" is not well-defined. The honest interpretation is that
$K$ is a hyperparameter chosen for descriptive convenience, not a
biological truth claim. Reporting results at several values of $K$ and
showing how the picture changes is best practice.

#tip[
  The first sanity check on an admixture run is reproducibility across
  random seeds. ADMIXTURE's quasi-Newton optimiser has multiple local
  optima for $K >= 3$; running with three different seeds and
  checking that the $Q$ matrices agree on the dominant components is
  the cheapest way to detect that something is wrong with the fit.
]

=== Joint Inference and the Demographic Model Zoo

Real demographic inference combines multiple data signals — pairwise
heterozygosity (PSMC), joint site-frequency spectrum (`δaδi` from
Gutenkunst et al. 2009; `momi2` from Kamm et al. 2020), inferred
ARGs (`Relate` plus `tsdate`), and explicit pairwise statistics like
$F_("ST")$ and the $f$-statistics ($f_2$, $f_3$, $f_4$ from Patterson
et al. 2012). Modern human demographic models have ten to twenty
parameters: out-of-Africa split time, ancestral African $N_e$,
European-Asian split, sequence of $N_e$ changes per branch, Neanderthal
introgression timing, and so on. Fitting them well is an active
research area; competing reconstructions disagree at the few-thousand-
year level but agree on the overall shape.

#warn[
  Demographic models are degenerate in ways that can fool naïve
  fitting. An ancient bottleneck and a recent expansion produce
  similar site-frequency spectra; high $N_e$ in the ancestral
  population and low $N_e$ in one branch can mimic many ancestral
  splits. Cross-validating against an independent statistic — usually
  $F_("ST")$ or pairwise heterozygosity in a region not used for
  fitting — is essential before publishing a demographic model.
]


== Selection Scans <sec:selection-scans>

Population-genetic theory tells you what genetic variation should look
like under neutrality. Departures from that expectation, distributed
non-randomly along the genome, are the signature of natural
selection. *Selection scans* are systematic searches across the
genome for these departures. They are detection problems in the
classical signal-processing sense: a known noise distribution (the
neutral coalescent under a fitted demography), a test statistic at
each location, and a multiple-testing correction over the millions of
locations tested.

=== Tajima's D and the Site-Frequency Spectrum

The oldest neutrality test in popular use is *Tajima's D* (Tajima
1989). It compares two estimators of $theta = 4 N_e mu$ that weight
branches of the coalescent tree differently. Watterson's estimator
$theta_W = S / a_n$, where $a_n = sum_(k = 1)^(n - 1) 1 / k$, weights
total tree length. Tajima's $pi$ — average pairwise differences —
weights internal branches. Under strict neutrality both estimate
$theta$ and their difference is approximately zero. The standardised
test statistic
$ D = (pi - theta_W) / sqrt("Var"(pi - theta_W)) $
follows an approximately standard-normal null under neutrality.

Departures from $D approx 0$ have specific interpretations:

- *D \> 0*: excess of intermediate-frequency variants. Suggests
  *balancing selection* (heterozygote advantage, frequency-dependent
  selection) or recent population contraction.
- *D \< 0*: excess of rare variants. Suggests a *recent positive
  selection sweep* (which wipes out variation around the swept site,
  then new mutations accumulate as rare singletons) or *population
  expansion* (recent growth introduces many low-frequency variants).

Tajima's D is computed in sliding windows along the genome and
ranked. Tools like `vcftools --TajimaD` produce window-level
statistics; the top and bottom percentiles are candidate regions for
balancing or directional selection respectively.

=== iHS, XP-EHH, and the Haplotype-Length Tests

Tajima's D is sensitive to changes in the allele-frequency spectrum
but not to *haplotype structure*. A second class of tests looks at the
length distribution of haplotypes carrying each allele. *Extended
haplotype homozygosity* (EHH, Sabeti et al. 2002) at a focal SNP
measures the probability that two chromosomes carrying the same
allele at the SNP remain identical out to a given distance. Under
neutrality, EHH decays with distance at a rate set by local
recombination. A recently swept allele still carries an unusually
long haplotype, because there has not been time for recombination to
break it up.

The *integrated haplotype score* (iHS, Voight, Kudaravalli, Wen,
and Pritchard 2006) integrates EHH over distance for the two alleles
at a focal SNP and compares them. Extreme iHS values flag SNPs whose
derived allele carries an unusually long haplotype relative to the
ancestral — the signature of an incomplete sweep, where the favoured
allele is at intermediate frequency but still on a long shared
chromosomal background. *XP-EHH* (Sabeti et al. 2007) compares EHH at
the same SNP across two populations; extreme values indicate
population-specific sweeps. iHS detects sweeps within a population;
XP-EHH detects sweeps differentiated between populations.

#figure(
  image("../../diagrams/lecture-12/12-selection-scan.svg", width: 95%),
  caption: [
    A genome-wide selection-scan Manhattan plot. Annotated peaks
    correspond to known sweep targets: $L C T$ (lactase persistence
    in Europeans), $S L C 24 A 5$ (skin pigmentation), $E D A R$
    (East-Asian hair thickness), $A B C C 11$ (earwax type). Even
    canonical hits appear as modest peaks; robust detection requires
    a properly calibrated null.
  ],
) <fig:selection-manhattan>

=== Interpretation and Failure Modes

Selection scans produce long lists of candidate loci. Most candidates
turn out to be wrong. Four failure modes dominate:

1. *Population-structure confounding.* If subpopulations have
   different allele frequencies, pooled scans flag loci that are
   structure rather than selection. The fix is per-subpopulation
   scans plus replication across populations.
2. *Demographic events mimic selection.* A bottleneck reshapes the
   site-frequency spectrum and haplotype structure in ways that look
   like selection. The fix is a neutral-demography null calibrated
   by coalescent simulation under a fitted demographic model.
3. *Background selection.* Removal of linked deleterious variants
   reduces diversity genome-wide in gene-rich regions. Without
   correction, this looks like positive-selection signal in
   precisely the regions where you most want to detect it.
4. *Multiple testing.* Across millions of SNPs, the tail of any null
   distribution contains many false positives at any fixed nominal
   significance.

Good practice combines several mitigations: use a simulated null
under the best-fit demography (`msprime` or `SLiM` produce these),
combine multiple scan statistics (Tajima's D, iHS, XP-EHH) and require
agreement, and validate top hits with orthogonal evidence —
GWAS associations, functional studies, ancient-DNA dating.

#note[
  A selection scan is detection across millions of genomic windows.
  The null is neutrality under the inferred demography; the test
  statistic is Tajima's D, iHS, or XP-EHH; the multiple-testing
  problem is the standard one. Specific wrinkles: the null
  distribution has no closed form, so coalescent simulations stand
  in for it; demographic mis-specification translates directly into
  inflated false-alarm rates; the test statistic correlates across
  nearby windows due to LD. This is the same noise-floor-
  characterisation problem as radar — the detector works as well as
  the estimate of ambient clutter.
]

#warn[
  Selection scans are the genomics method most prone to
  over-interpretation. A locus flagged at the 99.9th percentile of
  iHS is a _candidate_, not a finding. Several published "selection
  sweeps" have been retracted when more-careful demographic
  modelling explained the signal without invoking selection. Always
  triangulate: multiple statistics, replication across populations,
  functional or GWAS evidence, a coherent evolutionary story.
]


== Summary <sec:summary>

- Allele frequencies are the basic currency of population genetics.
  Hardy-Weinberg gives the null distribution of genotype frequencies
  from a single allele frequency; deviations are the biological
  signal. HWE is never exactly true and is best read as a screen.
- The Wright-Fisher model is a discrete-time, discrete-state Markov
  chain on allele counts with a binomial kernel. Two absorbing
  states (fixation of either allele) make every finite WF population
  eventually monomorphic. Drift, selection, mutation, and migration
  are the four forces that perturb the kernel.
- Effective population size $N_e$, not census $N$, controls drift
  and selection thresholds. The product $2 N_e s$ determines whether
  an allele behaves deterministically or near-neutrally. Humans have
  $N approx 8 times 10^9$ but $N_e approx 10^4$.
- Linkage disequilibrium is genome-scale autocorrelation. $r^2$
  decays exponentially with $2 c t$; $E[r^2] approx 1 / (1 + 4 N_e c)$
  at equilibrium. LD blocks are autocorrelation plateaus separated
  by recombination hotspots — the empirical substrate of GWAS.
- The coalescent runs Wright-Fisher backward in time. Waiting times
  to the next coalescence among $k$ lineages are exponential with
  rate $binom(k, 2) / (2 N_e)$. TMRCA plateaus at $4 N_e$ as sample
  size grows because the final 2-to-1 merger dominates.
- PSMC infers $N_e(t)$ over a million years from a single diploid
  genome by modelling local coalescent time as an HMM. SMC++,
  `Relate`, and `tsinfer` extend the approach to multi-genome and
  ARG-based inference.
- ADMIXTURE is unsupervised mixture-model decomposition over $K$
  ancestral populations. $K$ is a hyperparameter, not a biological
  truth — choose it by cross-validation and report sensitivity to
  the choice.
- Selection scans (Tajima's D, iHS, XP-EHH) are detection problems
  with controlled false-alarm rates. The null must be a simulated
  neutral coalescent under the fitted demography, not an analytical
  approximation.


== Exercises <sec:exercises>

#strong[1.] #emph[HWE chi-squared by hand.] A cohort of 500 individuals
at a biallelic SNP shows genotype counts $A A = 220$, $A a = 200$,
$a a = 80$. Compute the sample allele frequency $p$, the expected
counts under HWE, the chi-squared statistic, and the corresponding
$p$-value. In which direction does the deviation point, and what
biological hypotheses are consistent with that direction?

#strong[2.] #emph[Drift variance and fixation probability.] A new
neutral mutation appears in a diploid population of $N_e = 5000$.
(a) What is the probability that the mutation eventually fixes? (b)
What is the expected fixation time conditional on fixation, in
generations and in years (assume a 25-year generation)? (c) The same
mutation in an $N_e = 5 times 10^5$ population: what are the answers?

#strong[3.] #emph[$2 N_e s$ critical scaling.] A mildly deleterious
variant has selection coefficient $s = -10^(-4)$. In which of the
following populations does drift dominate selection: human
($N_e approx 10^4$), domestic dog ($N_e approx 10^3$), _E.~coli_
($N_e approx 10^9$)? Compute $2 N_e |s|$ in each case and
characterise the regime.

#strong[4.] #emph[LD decay arithmetic.] Two SNPs are 50 kb apart in a
region with the human-average recombination rate. (a) Compute the
per-generation recombination probability $c$. (b) Starting from
$r^2(0) = 1$ at the moment of mutation, how many generations are
required for $r^2$ to fall below 0.1? (c) What does that translate
to in years, assuming 25 yr per generation? Comment on whether the
LD pattern around a 1000-year-old founder mutation should be visible
genome-wide.

#strong[5.] #emph[Coalescent waiting times.] Compute $E[T_k]$ for
$k = 2, 3, ..., 10$ in units of $2 N_e$ generations. Sum them to get
$E[T_("MRCA")]$ for a sample of 10. What fraction of the total time
is accounted for by the final 2-to-1 coalescence? Sketch how the
fraction would change if $n$ were 100 instead of 10.

#strong[6.] #emph[Tajima's D direction.] For each of the following
scenarios, predict the sign of Tajima's D: (a) a population that
underwent a strong recent bottleneck; (b) a population that has been
expanding for the last 50 000 years; (c) a region near a recently
swept advantageous allele; (d) a region under long-term balancing
selection. Justify each in terms of the site-frequency spectrum.

#strong[7.] #emph[ADMIXTURE K-selection.] You run ADMIXTURE at
$K = 2, 3, ..., 8$ on a worldwide reference panel. Cross-validation
error decreases sharply from $K = 2$ to $K = 4$, is approximately
flat from $K = 4$ to $K = 6$, and rises slightly thereafter. Which
$K$ do you report, and what caveats should accompany the result?
What additional information would let you settle the question?

#strong[8.] #emph[(Open-ended.)] Pick one of the tools mentioned in
the chapter — `PSMC`, `SMC++`, `Relate`, `ADMIXTURE`, `msprime`, or
another — and read its primary publication. In one paragraph,
describe the single design choice the authors made that is most
specific to population-genetic data (as opposed to a generic ML or
statistical method), and explain why it works on the empirical data
they show.


== Further Reading <sec:further-reading>

- #strong[Hardy, G. H.] (1908). "Mendelian Proportions in a Mixed
  Population." #emph[Science] 28: 49–50. Half a page; still the
  cleanest derivation of HWE.
- #strong[Fisher, R. A.] (1930). #emph[The Genetical Theory of
  Natural Selection.] Oxford University Press. The foundational
  text of mathematical population genetics; Chapters 1–4 cover the
  WF model and selection.
- #strong[Wright, S.] (1931). "Evolution in Mendelian Populations."
  #emph[Genetics] 16: 97–159. Wright's parallel construction of the
  same framework with a much heavier emphasis on drift.
- #strong[Kingman, J. F. C.] (1982). "The Coalescent."
  #emph[Stochastic Processes and their Applications] 13: 235–248.
  Four pages of pure probability; the paper that founded coalescent
  theory.
- #strong[Tajima, F.] (1989). "Statistical Method for Testing the
  Neutral Mutation Hypothesis by DNA Polymorphism." #emph[Genetics]
  123: 585–595. The original Tajima's D paper.
- #strong[Li, H., and Durbin, R.] (2011). "Inference of Human
  Population History from Individual Whole-Genome Sequences."
  #emph[Nature] 475: 493–496. The PSMC paper — one diploid, a
  million years of demography.
- #strong[Alexander, D. H., Novembre, J., and Lange, K.] (2009).
  "Fast Model-Based Estimation of Ancestry in Unrelated
  Individuals." #emph[Genome Research] 19: 1655–1664. ADMIXTURE.
- #strong[Voight, B. F., Kudaravalli, S., Wen, X., and Pritchard,
  J. K.] (2006). "A Map of Recent Positive Selection in the Human
  Genome." #emph[PLoS Biology] 4: e72. iHS, introduced with a
  worldwide application.
- #strong[Kelleher, J., Etheridge, A. M., and McVean, G.] (2016).
  "Efficient Coalescent Simulation and Genealogical Analysis for
  Large Sample Sizes." #emph[PLoS Computational Biology] 12:
  e1004842. `msprime` and tree sequences.
- #strong[Slatkin, M.] (2008). "Linkage Disequilibrium —
  Understanding the Evolutionary Past and Mapping the Medical
  Future." #emph[Nature Reviews Genetics] 9: 477–485. A compact
  review of LD as both evolutionary signal and GWAS substrate.
- #strong[1000 Genomes Project Consortium] (2015). "A Global
  Reference for Human Genetic Variation." #emph[Nature] 526:
  68–74. The empirical resource on which most modern human
  population-genetic analysis stands.
