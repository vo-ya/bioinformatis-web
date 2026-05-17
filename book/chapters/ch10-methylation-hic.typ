#import "../theme/book-theme.typ": *

= Methylation, Hi-C, and 3D Genome Organisation <ch:methylation-hic>

#matters[
  The genome you have spent nine chapters reading is a one-dimensional
  string of letters. The genome inside a cell is something else: a
  two-metre polymer folded into a five-micron nucleus, decorated with a
  second alphabet of chemical marks that does not appear in the FASTQ,
  organised into chromatin neighbourhoods and territories that decide
  which gene gets to talk to which enhancer. Two assays — bisulfite
  sequencing for the chemical marks, Hi-C for the folding — read both
  of these missing dimensions back into text files an aligner can chew.
  Each comes with a statistical problem dressed up in unfamiliar
  vocabulary, and each of those problems is, when you look at it
  squarely, a piece of textbook engineering math: a channel-decoding
  problem, a Bayesian proportion estimator, a covariance matrix's
  principal component, a change-point detector. This chapter does the
  translation.
]

DNA methylation and three-dimensional chromosome conformation sit in
the same chapter because both are layers of biology that the linear
sequence cannot tell you about — and because both, despite being
biologically and chemically unrelated, end up posing remarkably similar
inference problems. Methylation is a chemical mark on individual
cytosines, read out by sequencing after a channel-induced symbol flip.
Hi-C is a population-averaged measurement of which loci touch which
other loci in space, read out as an $N times N$ contact matrix. The
methylation pipeline boils down to decoding a channel and estimating
proportions; the Hi-C pipeline boils down to balancing a matrix and
running an eigendecomposition. By the end of the chapter the same
piece of linear algebra you would apply to a covariance matrix in any
EE 205 lab will partition a human chromosome into A and B compartments.

The chapter has two halves. The first (sections 10.1 through 10.4)
covers methylation: what the mark is, how bisulfite chemistry reads
it, how the three-letter alignment trick rescues otherwise-broken
aligners, and how Beta–Binomial models turn per-CpG counts into
differential-methylation calls. The second (sections 10.5 through 10.9)
covers Hi-C: what the protocol measures, what the contact matrix looks
like after iterative balancing, how TADs fall out of an insulation
score, and how A/B compartments fall out of the first eigenvector. The
last section pulls together the three linear-algebra framings — Hi-C
as covariance, compartments as PCA, 3D reconstruction as MDS — that
make the whole assay legible to anyone with a signals-and-systems
background.


== What DNA Methylation Is <sec:methylation-biology>

DNA methylation is a small chemical modification: a methyl group ($"CH"_3$)
covalently attached to the fifth carbon of cytosine, producing
*5-methylcytosine* (5mC). In every mammalian genome the modification
occurs almost exclusively at cytosines that sit immediately upstream
of a guanine on the same strand — the *CpG dinucleotide*, where the
"p" is a literal phosphodiester bond, not a base. Plants methylate
more permissively; mammals are CpG-centric.

The DNA alphabet at the level of individual bases is still
$\{A, C, G, T\}$. Layer the methyl mark on top and each C position
acquires a hidden bit: methylated or not. An effectively five-symbol
alphabet, if you want to track epigenetic state alongside sequence.
Every assay you read about in Chapters 1 through 9 throws the hidden
bit away — base calling ignores it, alignment ignores it, variant
calling ignores it. To see methylation you need an assay that converts
the hidden bit into a sequencing-visible signal, which is the topic of
the next section.

Methylation matters because it is stable, heritable, and regulatory.
Once a CpG is methylated, it usually stays methylated through many
cell divisions; the maintenance methyltransferase DNMT1 copies the
mark from parent strand to daughter strand after every replication
fork passes. The mark is therefore inherited the way DNA itself is
inherited, but the rules are softer — methylation states can change,
new methylation can appear, and the TET family of enzymes can
actively erase the mark by oxidising 5mC to 5-hydroxymethylcytosine
(5hmC) and then to other intermediates.

#figure(
  image("../../diagrams/lecture-10/01-methylation-chemistry.svg", width: 95%),
  caption: [
    5-methylcytosine. The same cytosine base, plus a methyl group at
    the 5-carbon. DNMTs add the mark; TET enzymes erase it via 5hmC
    intermediates. The same sequence can sit in two epigenetic states.
  ],
) <fig:methylation-chemistry>

The regulatory consequence depends on where in the genome the
methylation lands. Methylation at a gene's promoter typically silences
the gene: methylated CpGs recruit methyl-CpG-binding proteins (the
MBD family), which in turn recruit chromatin-compacting machinery,
which in turn keeps RNA polymerase from getting a foothold. Demethylation
at the same promoter releases the suppression and allows transcription.
The same sequence, read at two different methylation states, expresses
different genes. This is why cells with identical genomes — every
nucleated cell in your body — can become very different cell types:
the wiring diagram is shared, but the switches laid over it are not.

#note[
  Methylation is a bit strapped to each CpG. The DNA sequence alone
  does not tell you whether a gene is on or off, only what protein it
  could in principle produce. The methylation pattern is a sparse
  overlay that decides whether that protein actually gets produced in
  this cell at this time. Treat the sequence as the schematic and the
  methylation pattern as the breaker panel: the schematic does not
  change between rooms in the house, but the breakers do.
]

=== CpG Islands and the Methylation Landscape

Most CpGs in the mammalian genome are methylated — somewhere between
70 % and 85 % of all CpGs in a typical somatic cell carry the mark.
This seems paradoxical at first, because methylation makes the affected
cytosine far more prone to deamination: 5mC → T happens at perhaps ten
times the rate that unmodified C → U deamination happens, and the T is
not detected as damage in the way that U is. Over evolutionary time,
CpGs in heavily-methylated regions are slowly grazed away. CpGs are
roughly $5 times$ depleted relative to the equiprobable-base expectation
across mammalian genomes — a direct fossil of historical methylation.

A minority of CpGs survive this deamination grinder, concentrated in
short (200 bp to 3 kb) regions with locally high GC content (above
50 %) and high observed-to-expected CpG ratio (above 0.6). These are
*CpG islands*. They survive because they are usually unmethylated;
no mark, no accelerated deamination, no loss. About 60 to 70 percent
of human gene promoters overlap a CpG island, and in normal cells
those CpG islands are kept unmethylated regardless of whether the
linked gene is actively transcribed.

The global picture has three tiers, ordered from rarest and most
sensitive to most common and least dynamic:

- *CpG islands* (\~1 % of all CpGs): usually unmethylated; methylation
  gain silences the linked gene. The dynamic, regulatory tier.
- *CpG island shores* (the \~2 kb flanking each island): intermediate
  methylation; this is where most tissue-specific differences live.
- *Gene bodies and intergenic regions* (the remaining bulk):
  uniformly methylated at \~80 % or more; changes here are small and
  mostly correlative rather than causal.

#figure(
  image("../../diagrams/lecture-10/02-cpg-islands.svg", width: 95%),
  caption: [
    The three-tier methylation landscape. Most CpGs are methylated;
    promoter-overlapping CpG islands are unmethylated; the shores
    between them carry the tissue-specific signal.
  ],
) <fig:cpg-islands>

For a few hundred *imprinted* genes the rules change. The allele
inherited from the father is methylated on its promoter while the
maternal allele is unmethylated, or vice versa; only one of the two
parental alleles ever speaks. The asymmetry is set during germ-cell
development by sex-specific methylation at *differentially methylated
regions* (DMRs) and faithfully maintained through embryogenesis and
the rest of the organism's life. Loss of imprinting underlies several
clinical syndromes — Prader–Willi and Angelman both involve a
particular cluster on chromosome 15q11 — and shows up routinely in
cancer biology, where imprinted regions lose their parent-of-origin
methylation patterns as part of the broader epigenomic disorganisation
of tumour cells.

#note[
  Methylation as a gene-regulatory mechanism was proposed in the
  mid-1970s — Riggs in 1975, Holliday and Pugh independently the same
  year — well before sequencing was cheap enough to measure it
  genome-wide. The field spent twenty-five years studying individual
  loci one at a time, with locus-specific PCR and Sanger sequencing
  after bisulfite treatment. The 2008–2012 explosion of
  whole-genome bisulfite sequencing turned methylation from a
  per-locus into a per-genome assay — the same kind of transition that
  ChIP-seq caused for transcription-factor binding, in roughly the
  same five-year span. The Roadmap Epigenomics project assembled the
  first pan-tissue reference methylomes by 2015; a modern WGBS
  experiment covers all \~28 million human CpGs in one run.
]


== Bisulfite Sequencing and the Three-Letter Alphabet <sec:bisulfite>

The problem is straightforward to state and harder than it looks. A
sequencer reports base identities, not chemical modifications. To
read methylation, the chemistry has to translate the hidden bit
("methylated or not") into something the base caller can see —
ideally, into a different base. *Bisulfite conversion* (Frommer et al.,
1992) does exactly that.

The trick rests on a difference in chemical reactivity between
modified and unmodified cytosine. Treat denatured single-stranded DNA
with sodium bisulfite ($"NaHSO"_3$) under mildly acidic conditions and
elevated temperature, and two things happen:

1. Unmethylated cytosines react readily. The bisulfite ion adds across
   the 5,6 double bond, then the resulting intermediate deaminates,
   then the bisulfite leaves. The net transformation is C → U.
2. 5-methylcytosine does not react. The methyl group at C5 sterically
   blocks the bisulfite addition. 5mC survives the treatment unchanged
   and is read by the sequencer as C.

Run a PCR after bisulfite treatment and uracil is templated as
thymine, because U pairs with A and that A on the new strand pairs
with T on the strand after that. By the time the library reaches the
sequencer, the chemistry has carved a clean per-position record into
the base alphabet itself: every position that was originally
methylated reads as C, and every position that was originally
unmethylated reads as T. A surviving C in a read is a methylation
call; a C → T mismatch is an unmethylation call. The methylation bit
has become a base identity.

#figure(
  image("../../diagrams/lecture-10/03-bisulfite-chemistry.svg", width: 95%),
  caption: [
    Bisulfite conversion. Unmethylated C converts to U and then,
    through PCR, to T; methylated C survives unchanged. Every C in
    the post-PCR read is a methylated C in the original molecule.
  ],
) <fig:bisulfite>

#note[
  Bisulfite conversion is a *known input-dependent symbol substitution*.
  At each C position in the input sequence, the substitution depends
  on a hidden bit. Recasting it in communication-theory terms: the
  DNA sequence is the transmitted signal, the methylation state is a
  side-channel modulation, and bisulfite conversion is a channel that
  deterministically maps (C, unmethylated) → T and (C, methylated) → C
  while leaving A, G, and T untouched. The receiver — the sequencer —
  observes the channel output. Your job downstream is to decode two
  things at once: the transmitted base sequence (easy if you have a
  reference genome) and the modulation per C position (the methylation
  call). The rest of the methylation pipeline is a long exercise in
  decoding this channel cleanly.
]

=== Why Standard Aligners Break

A vanilla short-read aligner — BWA, Bowtie2 — was designed for
Chapter 2's problem: place a read against a reference under a small
random-error model. Most positions match; mismatches are rare,
quasi-Gaussian, and concentrate at high-error positions (read ends,
problematic motifs). A handful of mismatches per 150 bp read is fine.
A read with thirty mismatches is rejected outright.

A bisulfite read violates the model violently. Take a 100 bp reference
region containing twenty Cs. After conversion, somewhere between zero
and twenty of those Cs have flipped to T. A typical read shows
fifteen to twenty C → T positions. BWA with default parameters either
fails to align the read at all (too many mismatches in its
ungapped-extension scoring) or aligns it to whichever copy of the
reference happens to have a similar C → T pattern by chance, which is
not the right location more often than it is.

The mismatches are not random and they are not symmetric. They appear
at a predictable subset of positions (the Cs in the reference) and
always in one direction (C to T). An aligner that models errors as a
sequencing-noise distribution gives every C → T mismatch the same low
probability as every other one, missing the structural regularity
that all of those mismatches share an underlying explanation.

=== Three-Letter Alignment

*Bismark* (Krueger & Andrews, 2011) introduced the trick that has
dominated the field ever since: project both the reference and the
read into a reduced three-letter alphabet where the bisulfite channel
acts as the identity map. Align in that reduced space. Then, once
the read is placed, go back to the original four-letter sequences and
read methylation off as base agreement at C positions.

In one cycle the algorithm:

1. Builds two converted versions of the reference: a *C-to-T converted*
   index in which every C is replaced by T, and a *G-to-A converted*
   index for reads originating from the reverse-bisulfite strand.
2. Builds two converted versions of each read: a C-to-T read for the
   forward case and a G-to-A read for the reverse case.
3. Aligns each converted read against each converted reference using
   an off-the-shelf short-read aligner (Bismark wraps Bowtie2). The
   read and the reference both live in $\{A, G, T\}$ or $\{A, C, T\}$,
   so any C → T mismatch caused by bisulfite is hidden by construction.
4. Picks the best of the four resulting placements as the call.
5. Re-reads the *original*, unconverted read against the original
   reference at the chosen location. Every C in the read that
   coincides with a C in the reference is a methylated CpG (or
   non-CpG C, if the genome studied has those). Every T in the read
   at a reference-C position is an unmethylated C.

The clever piece is the projection: the alignment difficulty disappears
in the reduced alphabet, and the modulation is decoded *conditional on*
a successful alignment. Splitting the joint estimation problem in two
sequential steps recovers tractability without giving up much accuracy.

#figure(
  image("../../diagrams/lecture-10/04-bismark-alignment.svg", width: 95%),
  caption: [
    Bismark's three-letter strategy. Project both the reference and
    the read into the reduced alphabet where the bisulfite channel
    is identity, align there, then decode methylation by comparing
    the original sequences at C positions.
  ],
) <fig:bismark>

Other major short-read methylation aligners — BSMAP (Xi & Li, 2009),
methylpy (Schultz et al., 2015), Walt (Chen et al., 2016) — differ in
how they implement the reduction (hashing with degenerate bases,
allele-specific calling, GPU acceleration) but share the same
projection principle. Bismark remains the most widely used in 2024;
methylpy is the default for the large-scale Roadmap Epigenomics
re-analyses.

#note[
  Three-letter alignment is *joint channel decoding by alphabet
  reduction*. The channel has two possible transition kernels at each
  C position — identity, or the deterministic flip — and you do not
  know which applies until you have decoded the message. The fix is
  to project both reference and observation into a sub-alphabet where
  the two kernels agree, decode the message in the sub-alphabet, and
  then revisit the C positions to recover the modulation bit. The
  same pattern shows up in soft-decision decoding (decode in a
  quantised constellation, refine soft estimates afterward) and in
  error-correcting codes whose minimum distance is dominated by one
  axis of the symbol space.
]

#warn[
  Bisulfite conversion is not 100 % efficient in the wet lab.
  Production protocols achieve 99.0 to 99.8 % conversion of
  unmethylated C, leaving a small population of unconverted Cs that
  the caller mistakes for methylated sites. The standard control is
  a fully-unmethylated spike-in — typically lambda phage DNA at
  0.5 % of input mass — whose residual C rate after the pipeline runs
  bounds the false-positive methylation rate. A 1 % failure to convert
  on a truly-unmethylated locus produces a spurious "1 %
  methylation" call that an inexperienced eye reads as a faint
  biological signal. Always check conversion before believing
  small effects.
]

=== Long-Read Direct Calling

A different approach emerged in 2019–2022 from the long-read
platforms: call methylation directly from the raw sequencing signal,
without a bisulfite step at all. *Oxford Nanopore* sequences DNA by
threading it through a protein pore and recording an electrical
current; the methyl group at 5mC alters the current signature subtly
but consistently enough that a methylation-aware base caller (Guppy's
`hac_modbases` model, then Dorado in later releases) emits both a base
identity and a methylation probability per position. *PacBio HiFi*
measures polymerase kinetics during base incorporation; 5mC produces a
small but reproducible interpulse-duration shift, called by a
specialised HiFi 5mC model.

Per-CpG calling accuracy on long-read direct methylation is now in
the 95 % range, against bisulfite's \~99 % at typical depths. The
gap is closing as base-calling models improve. The trade-off shifts
the comparison in long-read's favour for cases where bisulfite is
otherwise painful — cancer whole-genome sequencing where the same
library has to deliver SNVs, SVs, and methylation in one pass;
phased methylation across haplotypes; calls in long repeats that
short-read aligners cannot place; and other modifications (6mA, 4mC,
5hmC) that bisulfite cannot distinguish from each other.

#figure(
  image("../../diagrams/lecture-10/05-longread-methylation.svg", width: 95%),
  caption: [
    Long-read direct methylation calling. Both platforms read the
    mark off their native signal modality — current deflection for
    Nanopore, polymerase slowdown for PacBio — through
    methylation-aware base callers.
  ],
) <fig:longread-methylation>

The 2024 picture is a transition. Bisulfite remains the cheaper
per-base option and is still the default for dedicated methylation
panels. Long-read direct calling is rapidly displacing bisulfite in
labs that already have a long-read instrument running for other
reasons. Five years from now, the bisulfite chapter of a bioinformatics
textbook may well be historical.


== Differential Methylation: Beta-Binomial Posteriors <sec:differential-methylation>

A single methylation experiment produces, per sample, a count pair
$(m_i, t_i)$ at every CpG: $m_i$ reads in which the cytosine was
methylated, $t_i$ total reads covering the position. The naive
methylation estimate is $hat(p)_i = m_i / t_i$ — the fraction of
methylated reads — and at infinite depth this would be the answer.
At finite depth it is a noisy estimate of a bounded quantity, and
the standard frequentist tools that work on count data (Poisson,
negative binomial, log-normal) do not respect the $[0, 1]$ bound the
true methylation proportion lives in.

The first move in any sane methylation pipeline is to put a *Beta*
prior on the per-CpG proportion. The Beta distribution has density

$ f(p; alpha, beta) = (p^(alpha-1) (1-p)^(beta-1)) / "B"(alpha, beta), quad p in [0, 1] $

with $alpha > 0$ and $beta > 0$ as shape parameters and
$"B"(alpha, beta) = Gamma(alpha) Gamma(beta) / Gamma(alpha + beta)$ as the
normaliser. The Beta is the conjugate prior to the Binomial: if
$p tilde "Beta"(alpha, beta)$ a priori and you observe $m$ methylated
reads out of $t$ total, the posterior is

$ p | (m, t) tilde "Beta"(alpha + m, beta + t - m). $

The posterior mean is $(alpha + m) / (alpha + beta + t)$, a smoothed
version of the naive proportion estimate that shrinks toward the prior
mean when the count is small and toward the data when the count is
large. Posterior variance falls as $1/(alpha + beta + t + 1)$ and the
credible interval respects the boundary at every depth, which is
exactly the property we need.

A weak default prior is $"Beta"(1, 1)$ — the uniform distribution on
$[0, 1]$, encoding "I have no preference, all proportions equally
likely a priori." Stronger priors centred near 0 or 1 encode the
biological reality that most CpGs in mammalian genomes are nearly all
or nothing.

#figure(
  image("../../diagrams/lecture-10/06-beta-distribution.svg", width: 95%),
  caption: [
    Beta densities for four parameter regimes. The shape parameters
    encode prior beliefs about the proportion; the conjugacy property
    means the posterior is again Beta, with parameters updated by the
    observed counts.
  ],
) <fig:beta>

#note[
  The Beta-Binomial is the canonical Bayesian model for "what
  fraction of trials succeeded?" The Beta prior encodes belief about
  the proportion before observing data; the Binomial likelihood
  updates that belief by the observed counts; the posterior is again
  Beta by conjugacy. The same structure appears in any engineering
  problem that estimates a probability bounded in $[0, 1]$:
  bit-error rate on a noisy link, detection probability under
  additive noise, false-alarm rate of a classifier. Methylation
  inference is the genomic instance, with the per-CpG methylation
  fraction as the parameter and binomial reads as the data. The Beta
  framing keeps the estimator calibrated near the boundary, where a
  Gaussian approximation would assign positive probability to
  proportions outside $[0, 1]$.
]

=== Two Groups, One Region

The actual question in a methylation experiment is rarely "what is
the methylation level at this single CpG?" It is "where does the
methylation pattern differ between two groups of samples?" — tumour
versus matched normal, drug-treated versus vehicle, embryonic versus
adult. *Differentially methylated regions* (DMRs) are the unit of
output: contiguous stretches of CpGs at which the group-level
methylation differs systematically.

The framework follows naturally from the per-CpG Beta-Binomial. At
each CpG and in each group, fit a Beta posterior. Compare the
posteriors — if the credible intervals do not overlap, the groups
differ at that site; equivalently, perform a likelihood-ratio test on
the two-group Beta-Binomial model. Adjust for multiple testing across
the genome (Benjamini–Hochberg, the standard default). Cluster
adjacent significant CpGs into DMR calls, with a minimum-size
threshold to suppress single-site flukes.

Three R packages dominate practice:

- *methylKit* (Akalin et al., 2012). Per-site Fisher exact test or
  logistic regression with no information sharing between sites. The
  fastest of the three, the smallest in default sensitivity. Good
  for screening.
- *BSmooth* (Hansen et al., 2012). Smooths per-CpG proportions across
  a sliding genomic window before testing, on the rationale that
  methylation varies on scales of hundreds of base pairs rather than
  per-site. Uses a local-likelihood smoother and emits contiguous
  DMRs as the output.
- *DSS* (Feng et al., 2014). Fits a Beta-Binomial mixed model with
  per-site dispersion estimated by empirical-Bayes shrinkage —
  borrowing strength across sites the way DESeq2 borrows it across
  genes in the differential-expression chapter. Best power in
  typical WGBS comparisons; the modern default.

A worked example is easier than the prose. Suppose at one CpG you
observe $(m_A, t_A) = (8, 10)$ in group $A$ and $(m_B, t_B) = (2, 10)$
in group $B$, with a weak $"Beta"(1, 1)$ prior. The posterior in
group $A$ is $"Beta"(9, 3)$ with mean $0.75$ and 95 % credible
interval roughly $[0.43, 0.93]$; the posterior in group $B$ is
$"Beta"(3, 9)$ with mean $0.25$ and 95 % credible interval roughly
$[0.07, 0.57]$. The intervals overlap modestly, the posterior
difference $hat(p)_A - hat(p)_B$ is 0.50 with substantial uncertainty,
and the call is marginally significant at this single CpG. Five
adjacent CpGs all showing the same pattern would compound to a
confident DMR; a single isolated CpG with this much depth, taken
alone, is suggestive but not strong.

#figure(
  image("../figures/ch10/f1-beta-posterior-update.svg", width: 95%),
  caption: [
    Beta posterior update for two groups at one CpG. The shapes
    sharpen and shift as data accumulates; the overlap between the
    two posteriors quantifies whether the groups differ. Cumulative
    evidence across adjacent CpGs is what consolidates a single
    suggestive site into a confident DMR call.
  ],
) <fig:beta-posterior>

#warn[
  Bulk methylation is averaged over a mixed tissue. A sample with two
  cell types — one 80 % methylated, one 0 % methylated, in a 50/50
  mix — produces a clean 40 % bulk reading that corresponds to no
  real cell state. The same 40 % could arise from every cell being
  truly intermediate, which is biologically very different. Always
  ask whether your "intermediate methylation" sites are mixtures
  rather than genuine partial methylation, especially in tumour
  samples where infiltrating immune cells routinely sit alongside
  tumour cells at variable proportions. Single-cell methylation
  assays (scBS-seq, scNMT-seq) resolve the ambiguity but are
  expensive and sparse.
]


== The Hi-C Protocol and the Contact Matrix <sec:hic-protocol>

The first nine chapters of this book treated the genome as a 1D
string. The physical reality is 3D: two metres of DNA per diploid
human cell, folded into a nucleus five microns across. Two loci that
sit a megabase apart on the linear sequence may end up half a micron
apart in three dimensions, or fifty nanometres apart, or in direct
contact. Gene regulation routinely depends on the difference. An
enhancer that lives a megabase upstream of its target promoter can
only activate the gene if the two stretches of DNA are physically
close — and whether they are close is decided by chromatin folding,
not by sequence.

*Hi-C* (Lieberman-Aiden et al., 2009) measures, for every pair of
genomic loci, how often they are in physical contact across a
population of cells. The output is a *contact matrix* $C in
"NN"^(N times N)$: for a discretised genome with $N$ bins of fixed
size (10 kb is the modern standard, with deep experiments going as
fine as 1 kb), the entry $C[i, j]$ counts paired-end sequencing reads
that map with one mate inside bin $i$ and the other inside bin $j$.
Large $C[i, j]$ means the two bins are frequently in close
three-dimensional contact; small $C[i, j]$ means they are rarely
together.

Hi-C is intrinsically a population assay. Any single cell at any
single moment carries a sparse, stochastic contact map. The matrix
you analyse is the sum across the millions of cells in your input
sample — an ensemble average of the population's instantaneous
configurations. Per-cell variability is real (single-cell Hi-C
exists, at a much higher per-cell cost) and shows up in the data as
diffuse noise around the population structure, but the chapter
treats Hi-C as a population measurement throughout.

=== The 3C Family

Hi-C is the all-vs-all member of a family of progressively-more-parallel
proximity assays, all built on the same chemistry trick — formaldehyde
crosslink, restriction-enzyme digest, proximity ligation, sequence
the resulting chimeric fragments — but differing in how many pairs
they can measure per experiment.

- *3C* (Dekker et al., 2002). Chromosome Conformation Capture. PCR a
  single predefined locus pair from the ligation library. One pair
  per experiment. The foundational paper.
- *4C* (Simonis et al., 2006). Circular 3C. Self-circularise the
  library, inverse-PCR from a single "viewpoint" locus, sequence its
  partners. One viewpoint versus everything else.
- *5C* (Dostie et al., 2006). Ligation-mediated amplification from
  many primers, up to a few hundred thousand pairs. Many-vs-many but
  over a pre-selected set of loci.
- *Hi-C* (Lieberman-Aiden et al., 2009). Biotin-tagged ligation
  junctions plus streptavidin pulldown plus Illumina sequencing. The
  first all-vs-all genome-wide contact matrix. Changed everything.
- *Micro-C* (Hsieh et al., 2015). Hi-C with MNase digestion to
  nucleosomal resolution instead of a restriction enzyme. Resolves
  features down to \~200 bp at the cost of more sequencing per cell.

#figure(
  image("../../diagrams/lecture-10/07-3c-family-tree.svg", width: 95%),
  caption: [
    The 3C-family genealogy. Each generation adds a dimension of
    parallelism. Hi-C's jump from 5C was the first time the full
    all-versus-all contact matrix became measurable in one experiment.
  ],
) <fig:3c-family>

=== The Chemistry, Step by Step

The Hi-C wet-lab protocol has six steps. Each one removes a different
piece of the linear-sequence prior and substitutes a 3D-proximity
signal in its place.

1. *Crosslink.* Treat intact cells with formaldehyde, which forms
   covalent bonds between proteins and DNA, and between proteins and
   proteins, that are close in three-dimensional space at the moment
   of fixation. The crosslinked chromatin "freezes" the 3D contact
   pattern.
2. *Digest.* Cut the crosslinked chromatin with a restriction enzyme
   — MboI or DpnII (4-bp cutters, denser cuts) or HindIII (a 6-bp
   cutter used in the original 2009 paper). The cuts leave sticky
   ends that stick up from the cross-linked nuclear matrix without
   diffusing away.
3. *Fill-in and biotin-tag.* Klenow polymerase fills the sticky ends,
   incorporating biotinylated dCTP in the process. Every cut site is
   now flagged with a biotin handle.
4. *Proximity-ligate.* T4 DNA ligase rejoins free ends. Because the
   ends are still tethered to the cross-linked nuclear matrix, ligation
   happens preferentially between ends that are physically close in
   3D — neighbours on the chromatin polymer, not neighbours on the
   linear sequence. The spatial information enters the library here.
5. *Shear, pull down, sequence.* Reverse the crosslinks, shear DNA to
   \~300 bp fragments, and use streptavidin beads to enrich for
   fragments carrying a biotinylated ligation junction. The pulldown
   discards non-Hi-C ligation byproducts and concentrates the true
   proximity-ligation events. Sequence paired-end on Illumina.
6. *Align.* Map each read of a pair independently against the
   reference. The two mates land at two distinct genomic positions —
   the two loci that were close in 3D space at the moment of fixation
   — and the read pair adds one count to $C[i, j]$ where $i$ and $j$
   are the bins of the two mates.

#figure(
  image("../../diagrams/lecture-10/08-hic-protocol.svg", width: 95%),
  caption: [
    The Hi-C protocol. Crosslinking freezes 3D contacts; restriction
    plus biotinylation plus proximity ligation captures the contact
    information into chimeric fragments; streptavidin pulldown plus
    sequencing reads out which pairs of loci were close.
  ],
) <fig:hic-protocol>

A modern Hi-C experiment runs at about 500 million paired-end reads
per sample, yielding a 10 kb-resolution contact matrix at acceptable
signal. The deepest published Hi-C and Micro-C experiments push to
several billion reads and 1 kb resolution.

=== The Contact Matrix Up Close

For a human genome binned at 10 kb the matrix is roughly
$290{,}000 times 290{,}000$. Stored densely the matrix has on the order
of $10^11$ cells. In practice the matrix is stored sparse: HDF5 with
the `cooler` schema (Abdennur & Mirny, 2020) or the `.hic` binary
format (Durand et al., 2016). Both formats index by bin coordinate so
that arbitrary genomic regions can be loaded as small dense submatrices.

The raw matrix is dominated by three features, in roughly decreasing
intensity:

- *A strong diagonal.* Bins adjacent on the linear sequence are
  almost always physically adjacent in space too. The diagonal is the
  noisiest line in any Hi-C visualisation.
- *Distance decay off-diagonal.* Contact frequency falls with
  genomic distance as a power law, $C[i, j] tilde lr(|i - j|)^(-gamma)$
  with $gamma$ in the neighbourhood of 1, because random thermal
  fluctuations of the chromatin fiber couple nearby loci more easily
  than distant ones. The decay is monotonic and smooth.
- *Structured deviations from the decay.* Away from the diagonal and
  beyond what the smooth decay predicts, the matrix carries the
  interesting biology — TADs, loops, A/B compartments. These features
  are what the rest of the analysis pipeline tries to extract.

#figure(
  image("../../diagrams/lecture-10/09-contact-matrix.svg", width: 95%),
  caption: [
    Raw and ICE-normalised contact matrices for the same region. The
    banding in the raw matrix is bin-specific bias; iterative
    normalisation removes it and the underlying TAD blocks and
    compartment checkerboard become visible.
  ],
) <fig:contact-matrix>

A raw map also carries two kinds of artefact that must be cleaned up
before structural interpretation:

- *Bin-specific coverage bias.* Some bins generate systematically
  more or fewer mapped read pairs because of mappability (a bin
  containing repeats accepts fewer unique mappings), GC content
  (PCR-amplification efficiency is GC-dependent), and
  restriction-site density (a bin with very few cut sites cannot
  contribute many ligation events anywhere). These biases appear as
  bright or dark bands running across whole rows and columns of the
  raw matrix.
- *Distance-decay background.* The power-law decay along the diagonal
  is real but uninteresting for most analyses, where the question is
  "is this contact stronger or weaker than expected for two loci
  separated by this distance?" Removing the smooth background is a
  prerequisite for compartment and loop calling.

The next section handles both.


== ICE Balancing and Distance Detrending <sec:ice>

The dominant normalisation for Hi-C in 2024 is *Iterative Correction
and Eigenvector decomposition* (ICE, Imakaev et al., 2012). ICE makes
a single assumption about the biases: each bin has a "visibility"
factor $b_i$, and the observed contact is the product of the two bins'
visibilities times the true underlying contact:

$ C[i, j] = b_i b_j T[i, j]. $

If you know the $b_i$ you can invert to recover $T = C / (b b^top)$.
You do not know them, but you can fit them by demanding that the
balanced matrix have constant row and column sums — the
*equal-visibility* assumption. The iterative algorithm:

1. Set $w_i = 1$ for all bins.
2. Compute row sums $s_i = sum_j C[i, j]$ of the current matrix.
3. Set $w_i arrow.l w_i / s_i$ for every bin.
4. Multiply every entry: $C[i, j] arrow.l C[i, j] dot w_i w_j$.
5. Repeat until row sums are within tolerance of each other.

In practice the algorithm converges in twenty to thirty iterations.
The final $w_i$ are the inverse visibility factors; the balanced
matrix has uniform marginals, and the bin-specific banding that
dominated the raw map disappears.

*Knight-Ruiz balancing* (KR, Knight & Ruiz, 2013) is a faster
matrix-balancing algorithm from numerical linear algebra that solves
exactly the same problem — find a diagonal matrix $D$ such that
$D C D$ has constant row and column sums — by a Newton-type method.
KR converges in a handful of iterations rather than twenty and is the
default in the Juicer toolchain (Rao et al., 2014). The output of
ICE and KR agree up to numerical tolerance; the difference is
runtime, not output.

A separate normalisation step removes the distance-decay background.
Compute the empirical expected contact at each genomic separation
$E(|i - j|)$ — typically the mean of all on-diagonal entries at that
separation across a chromosome — and form the *observed-versus-expected*
matrix

$ "OE"[i, j] = C[i, j] / E(lr(|i - j|)). $

Entries above 1 indicate contacts above what distance alone would
predict; entries below 1 indicate suppression. The OE matrix flattens
the diagonal and brings out compartment- and loop-scale structure that
the raw decay obscured. A common visualisation choice is $log_2("OE")$,
which puts no-effect at zero and is symmetric around it on a colour
scale.

#figure(
  image("../figures/ch10/f2-ice-balancing.svg", width: 95%),
  caption: [
    Iterative correction of a Hi-C contact matrix. Each iteration
    rescales bin visibilities to drive row and column sums toward a
    uniform target; the visible banding in the raw matrix dissolves
    within twenty iterations into the underlying TAD-and-compartment
    structure.
  ],
) <fig:ice-balancing>

#note[
  Matrix balancing is a classical problem in numerical linear algebra
  with applications well outside genomics — input-output economic
  models, doubly-stochastic Markov chains, transportation matrices.
  Sinkhorn's theorem (1964) guarantees that any non-negative matrix
  with full support is balanceable by alternating row and column
  rescaling, and that the result is unique up to a positive scalar.
  Imakaev's ICE is Sinkhorn balancing applied to Hi-C; the eigenvalue
  decomposition in the name refers to a separate cleanup step the
  paper bundled in. Knowing the matrix-balancing literature, you can
  drop in any of half a dozen alternative balancing algorithms with
  no biology-specific reasoning required.
]


== TADs and the Insulation Score <sec:tads>

Zoom into a balanced contact matrix at the megabase scale and a new
pattern appears: triangular blocks of high contact density along the
diagonal, separated by sharper drops at their borders. These are
*topologically associating domains* (TADs), introduced as a concept
in two simultaneous papers (Dixon et al., 2012; Nora et al., 2012).
A TAD is a contiguous genomic region — typically 200 kb to 2 Mb —
inside which contacts are frequent and across whose boundaries
contacts are relatively rare. The boundary is biologically real:
genes inside the same TAD tend to be co-regulated, enhancers tend to
activate only promoters within the same TAD, and TAD boundaries are
enriched for CTCF binding sites that anchor the boundary on the
chromatin polymer (Chapter 9 covered CTCF in the context of ChIP-seq).

The structural question — given a balanced contact matrix, where are
the TAD boundaries? — is a 1D boundary-detection problem dressed up
in 2D clothing. Several algorithms exist; all of them reduce the
matrix to a 1D signal and then run change-point detection on the
result.

*Insulation score* (Crane et al., 2015) is the simplest of the
reductions and the most widely used in 2024. For every bin $i$, define
an insulation window of half-width $w$ bins centred at $i$. Sum the
contacts in the window that cross bin $i$ — that is, all
$(j, k)$ pairs with $j <= i < k$ inside the window — and divide by
the local mean to normalise. The score is low when bin $i$ sits at a
strong boundary, because few read pairs span across it; the score is
high when bin $i$ sits inside a TAD interior, because most pairs in
its neighbourhood span across it freely. Local minima of the
insulation score are TAD-boundary candidates; minima that fall below
a noise threshold are emitted as boundaries.

#figure(
  image("../../diagrams/lecture-10/12-tad-insulation.svg", width: 95%),
  caption: [
    TAD calling reduces a 2D contact matrix to a 1D insulation-score
    signal, then runs standard change-point detection on it. Local
    minima are boundary candidates; the threshold separates boundary
    minima from intra-TAD fluctuations.
  ],
) <fig:insulation>

*Directionality index* (Dixon et al., 2012) was the original
proposal: for each bin, compare upstream-direction contacts to
downstream-direction contacts and look for sign flips. *Arrowhead*
(Rao et al., 2014) transforms the matrix into a scale-independent
TAD-strength map by a sliding-window operator and reports peaks of
the resulting score. All three methods agree on the strong boundaries
in any given dataset; they disagree at the weak ones, where the
choice of algorithm and threshold matters more than the data does.

#note[
  TAD calling is *change-point detection in 1D after 2D-to-1D
  reduction*. The insulation score is the reduction: it collapses
  each row-and-column neighbourhood of the contact matrix into a
  single scalar per bin. Change-point detection on the resulting 1D
  signal is then off-the-shelf — PELT, binary segmentation, Bayesian
  online change-point detection all apply. The same problem appears
  in any signal where you suspect piecewise-constant means with
  abrupt transitions: GPS trajectory segmentation, speech-pause
  detection, control-loop fault detection. The genomic specificity
  is entirely in the reduction; the change-point machinery is
  identical.
]

#warn[
  TAD boundaries are *resolution-dependent*. The same dataset
  analysed at 10 kb bins shows a different set of boundaries from
  the same dataset analysed at 40 kb. A "TAD" at 40 kb often
  resolves into three to five sub-TADs at 10 kb; some weak
  boundaries at 10 kb vanish at 40 kb. Always report the bin size
  and the insulation-window width together with the call. "Conserved
  TAD boundaries between species" claims that compare different
  resolutions across the comparison are systematically misleading.
]

=== Loops and Focal Enrichment

At a finer scale than TADs sit *chromatin loops*: a specific
bin-pair $(i, j)$ contacts more frequently than its immediate
neighbourhood does, producing a bright focal spot in the contact
matrix away from the diagonal. Most loops span tens to hundreds of
kilobases and connect either two CTCF-bound anchor sites or an
enhancer to its target promoter. The loop catalogue from Rao et al.
(2014) — about 10,000 loops in human GM12878 cells — is the
field-standard reference.

The standard loop caller is *HICCUPS* (Rao et al., 2014). For each
candidate matrix entry $(i, j)$ at sufficient distance from the
diagonal, compare its count to four background windows in the immediate
neighbourhood (one in each diagonal direction). Test whether the focal
count is significantly brighter than its neighbours under a local
Poisson null. The statistical framing is — as elsewhere in this course
— a local-Poisson detection problem, structurally identical to MACS2's
peak detection in Chapter 9. Different domain, same algorithm.


== A/B Compartments as the First Eigenvector <sec:compartments>

Zoom out one further scale, to the whole-chromosome view, and the
balanced contact matrix reveals a checkerboard pattern: long stretches
of the chromosome contact each other preferentially across the whole
chromosome arm, partitioning into two interleaved sets. The two sets
correspond to the *A compartment* (euchromatin, gene-dense,
transcriptionally active, generally accessible chromatin) and the *B
compartment* (heterochromatin, gene-poor, late-replicating, generally
inaccessible). The bipartition is biological and large — A and B
compartments collectively organise the megabase-scale architecture of
the genome.

The algorithmic recipe to find them (Lieberman-Aiden et al., 2009) is
worth working through in full because every step is something you
have seen in linear algebra and the genomics-specific naming is a
thin veneer over standard math.

1. Start with the OE matrix $"OE"[i, j] = C[i, j] / E(|i - j|)$
   computed in section 10.5. Distance decay is gone; only structured
   deviations remain.
2. Compute the *correlation matrix* of the OE. For each pair of bins
   $(i, j)$, treat row $i$ and row $j$ as vectors in $RR^N$ and
   compute their Pearson correlation across all other bins. Call the
   resulting matrix $R in RR^(N times N)$: $R[i, j]$ is the
   correlation between the contact profiles of bins $i$ and $j$.
3. Compute the *first eigenvector* of $R$ — the unit vector $v_1$
   such that $R v_1 = lambda_1 v_1$ with $lambda_1$ the largest
   eigenvalue. Each bin gets one coordinate, $v_1[i]$.
4. The sign of $v_1[i]$ partitions bins into two sets. By convention
   the positive set is the A compartment and the negative set is the
   B compartment; fixing the sign requires an external reference
   (gene density per bin works well — A compartments are gene-dense),
   because the eigenvector is defined up to a sign flip.

That is the entire algorithm. The output is one scalar per bin,
positive or negative, with magnitude proportional to how strongly
that bin participates in its compartment. Plot $v_1[i]$ against
genomic coordinate and you see long alternating stretches of positive
and negative values — the chromosome's compartment structure.

#figure(
  image("../../diagrams/lecture-10/10-tads-compartments.svg", width: 95%),
  caption: [
    Three feature scales on one balanced Hi-C matrix. The
    Mb-scale plaid is the A/B compartment checkerboard;
    triangular blocks along the diagonal are TADs at the
    100 kb–Mb scale; focal bright spots near TAD corners are
    individual loops.
  ],
) <fig:tads-compartments>

#figure(
  image("../../diagrams/lecture-10/11-compartments-eigenvector.svg", width: 95%),
  caption: [
    Compartments as the first eigenvector. The observed-over-expected
    matrix produces a correlation matrix; the first eigenvector of
    that correlation matrix partitions the chromosome into A and B
    compartments by sign.
  ],
) <fig:compartments-eigenvector>

#note[
  The Hi-C balanced contact matrix is, structurally, a covariance
  matrix: symmetric, positive-semidefinite, with entries proportional
  to how correlated the contact-behaviour of two bins is across the
  ensemble of cells. The first eigenvector of its correlation matrix
  is its first principal component — the single direction in
  bin-space along which contact patterns vary most. The bipartite
  sign structure arises because compartment identity is the largest
  source of variance in chromatin contact patterns. This is PCA on a
  covariance matrix. The same math defines eigenfaces, the first mode
  of a distributed sensor array, and the principal axes of an inertia
  tensor. Compartmentalisation is not an exotic genomics technique; it
  is the most direct translation in this whole course between an
  engineering construct and a piece of biology.
]


== Hi-C as Covariance, 3D as MDS <sec:linalg-framings>

The Hi-C analysis pipeline, taken as a whole, is a sequence of
linear-algebra operations dressed up in epigenomics vocabulary. Three
framings make the translation explicit; each one is a piece of
textbook math.

=== Hi-C as a Covariance Matrix; Compartments as PCA

The balanced, distance-detrended contact matrix is a covariance
matrix of bin activities, with the ensemble of cells as the sampling
distribution. Its first principal component is the A/B
compartmentalisation. The previous section walked the algebra in
detail; the structural point is that nothing in the calculation is
specific to genomics. Power iteration, Lanczos tridiagonalisation, or
any other standard eigenvector method that handles large sparse
symmetric matrices works directly on the cooler-format input.

=== TAD Calling as Block-Diagonal Detection

The same matrix, read at finer scale, has a block-diagonal structure
along its principal diagonal. TAD calling is the detection of those
blocks. Three equivalent framings, depending on which kind of
engineering problem you have seen most recently:

- *Block-diagonal structure detection.* A matrix with strong
  on-diagonal blocks and weak off-diagonal entries is the canonical
  form for community detection on a graph (the adjacency matrix is
  block-diagonal if communities are well-separated), Gaussian-mixture
  covariance estimation, and image-segmentation problems on
  pixel-similarity graphs.
- *Change-point detection.* A TAD boundary is a change-point in the
  per-bin contact distribution. PELT, binary segmentation, and
  Bayesian online change-point detection apply directly to the
  insulation score signal, with the algorithmic choice mattering less
  than the choice of insulation-window half-width.
- *2D image segmentation.* The contact matrix read as an image has
  segments — contiguous pixel regions of similar intensity. Watershed
  and graph-cut segmentation algorithms have been applied to Hi-C
  matrices directly, and they recover much the same boundaries that
  insulation-score methods do.

All three framings are equivalent views of the same object. An
engineer who has handled any one of them already has the right
intuition for TAD calling; the genomic specifics are entirely in the
2D-to-1D reduction.

#figure(
  image("../figures/ch10/f3-hic-as-covariance.svg", width: 95%),
  caption: [
    The full algebraic chain from the balanced Hi-C contact matrix to
    the A/B compartment vector. Distance detrending produces the OE
    matrix; row-wise correlation produces the correlation matrix;
    the first eigenvector's sign labels each bin as A or B.
  ],
) <fig:hic-as-covariance>

=== 3D Reconstruction as MDS

If Hi-C gives you a matrix of pairwise contact frequencies, the next
natural question is whether you can invert it to recover the actual
3D coordinates of each bin in the nucleus. Treat high-contact pairs
as spatially close and low-contact pairs as spatially distant, and
the contact matrix transforms into a *distance matrix*. Given a
distance matrix, recovering Euclidean coordinates is the classical
*multidimensional scaling* (MDS) problem, with a closed-form solution
in terms of eigenvectors of the doubly-centred squared-distance
matrix.

The recipe is direct:

1. Convert contact frequencies to distances. The common choice is
   $d_(i j) = C[i, j]^(-alpha)$ for some $alpha in [0.5, 1]$ that
   the user calibrates against expected fiber dimensions. Other
   monotone-decreasing transforms work too.
2. Form the squared-distance matrix $D in RR^(N times N)$ with
   $D_(i j) = d_(i j)^2$.
3. Double-centre: $B = -1/2 H D H$ where $H = I - 1/N bold(1) bold(1)^top$
   is the centring matrix.
4. Compute the top $k$ eigenvectors and eigenvalues of $B$. The 3D
   coordinates are $X = U_k Lambda_k^(1/2)$ where $U_k$ stacks the
   top-three eigenvectors and $Lambda_k$ is the diagonal matrix of
   their eigenvalues.

The output is one 3D point per bin: a polymer trajectory through
nuclear space. *3DMax* (Oluwadare et al., 2018), *PASTIS*
(Varoquaux et al., 2014), and *GEM* (Zhu et al., 2018) all solve
variants of this MDS problem, differing in how they handle the
non-Euclidean constraints (the polymer is a chain, with physical
limits on how sharply it can bend) and the ensemble-averaging caveat.

#figure(
  image("../figures/ch10/f4-mds-3d-reconstruction.svg", width: 95%),
  caption: [
    3D reconstruction from a Hi-C contact matrix as classical MDS.
    Contact frequencies transform into distances; double-centring
    produces a Gram matrix; the top three eigenvectors are the 3D
    coordinates of each bin in nuclear space.
  ],
) <fig:mds-3d>

#note[
  The ensemble-averaging caveat matters. The "3D structure" recovered
  by MDS is a single best-fit conformation to a population-averaged
  contact matrix — the configuration that, if every cell adopted it
  simultaneously, would produce the observed bulk Hi-C. Real
  chromatin is heterogeneous: at any moment different cells in the
  population occupy different configurations, and the population
  average is not necessarily a state that any single cell ever
  visits. Single-cell Hi-C (Nagano et al., 2013) recovers per-cell
  conformations at much lower read depth per cell. For bulk Hi-C,
  treat the MDS reconstruction as a useful caricature rather than as
  the physical structure.
]


== Summary <sec:summary>

- Methylation adds a second alphabet on top of DNA sequence. Each
  cytosine carries a one-bit modulation — methylated or not — that is
  invisible to standard sequencing and visible only with a
  specialised assay.
- Bisulfite conversion is a known input-dependent symbol substitution.
  Unmethylated C becomes T; methylated C survives as C. Decoding it
  recovers the methylation pattern at single-base resolution.
- Three-letter alignment rescues otherwise-broken aligners by
  projecting both the reference and the read into a reduced alphabet
  where the bisulfite channel acts as the identity. The full
  four-letter sequences are revisited after alignment to decode the
  modulation.
- Differential methylation is Beta–Binomial. The per-CpG proportion
  is Beta-distributed; the per-read count is Binomial; conjugacy
  gives a closed-form posterior update that respects the $[0, 1]$
  boundary. methylKit, BSmooth, and DSS implement increasingly
  sophisticated information-sharing strategies on top of the same
  underlying model.
- Hi-C measures 3D proximity by proximity-ligation-plus-sequencing.
  The output is an $N times N$ contact matrix whose entries count
  read pairs spanning each bin pair.
- ICE and KR balancing remove bin-specific bias; distance detrending
  removes the power-law background. The result is a contact matrix
  ready for structural interpretation.
- TADs are block-diagonal structure detected by reducing the 2D
  matrix to a 1D insulation score and running change-point detection.
  Resolution-dependent.
- A/B compartments are the sign of the first eigenvector of the
  contact correlation matrix — PCA on a covariance matrix, with no
  genomics-specific math required.
- 3D reconstruction from Hi-C is classical multidimensional scaling:
  eigendecomposition of a centred squared-distance matrix, top three
  eigenvectors.


== Exercises <sec:exercises>

#strong[1.] #emph[Bisulfite decoding by hand.]
A reference region reads `AATCGTCGAAGCG`. Three sequencing reads from
the same region come back as `AATCGTCGAAGCG`, `AATTGTTGAAGTG`, and
`AATCGTTGAAGCG`. (a) For each read, list the methylation state at
every CpG in the reference (filled = methylated, open = unmethylated).
(b) What is the per-CpG methylation rate across the three reads?
(c) Suppose conversion efficiency is 99 %. What is the most likely
explanation for any C that "looks methylated" with only one
supporting read?

#strong[2.] #emph[Beta posterior arithmetic.]
A weak prior $"Beta"(1, 1)$ is updated by 12 methylated reads out of
20 total at one CpG. (a) Write the posterior parameters. (b) Compute
the posterior mean and 95 % credible interval (use the inverse Beta
CDF or, for an approximate answer, a normal approximation). (c) The
same CpG is sequenced to 200 reads in a separate sample, with 120
methylated. Without re-deriving from scratch, write the posterior for
the second sample under the same prior, and comment on how much the
credible interval narrows.

#strong[3.] #emph[Three-letter alignment by hand.]
Take the reference `AATCGTCGAAGCG`. Produce its C-to-T converted
version. Take the bisulfite read `AATTGTTGAAGTG` and produce its
C-to-T converted version. Verify that the two converted strings
match. Now use the original sequences to decode methylation at each
of the three CpG positions (counted from the left).

#strong[4.] #emph[DMR significance.]
At one CpG you observe $(m_A, t_A) = (40, 50)$ in group A and
$(m_B, t_B) = (10, 50)$ in group B. (a) Fit Beta posteriors under a
$"Beta"(1, 1)$ prior. (b) Compute an approximate 95 % credible
interval for $p_A - p_B$ by sampling. (c) The same effect size is
observed at 20 adjacent CpGs in a 10 kb window. Without running the
full DSS pipeline, sketch why the joint evidence is much stronger
than any single CpG.

#strong[5.] #emph[ICE balancing.]
A toy contact matrix is
$C = mat(10, 2, 1; 2, 30, 6; 1, 6, 20)$.
Run two iterations of ICE by hand: compute row sums, divide each
column $j$ by the row sum's geometric mean
$sqrt(s_i s_j)$, repeat. Compare the marginals of the result to the
marginals of the raw matrix. Do they become more equal?

#strong[6.] #emph[Insulation score by hand.]
A small $7 times 7$ contact matrix is given with TAD blocks at bins
1–3 and 4–7. Compute the insulation score at bin 4 using a window
half-width of 2 (so the window is bins 2 through 6). Compare the
score at bin 4 to the score at bin 5. Which is the boundary
candidate?

#strong[7.] #emph[Compartment eigenvector.]
Given a $4 times 4$ toy OE matrix (chosen so that bins 1, 3 form one
group and bins 2, 4 form another), compute the row-wise Pearson
correlation matrix by hand or in Python. Show that the first
eigenvector has alternating signs that recover the bipartition. Why
does the sign convention require an external reference (e.g. gene
density) to choose A versus B?

#strong[8.] #emph[(Open-ended.)]
Pick one Hi-C analysis package (`cooler`, `cooltools`, `Juicer`,
`HiCExplorer`, or another) and read its user guide for compartment
calling. In one paragraph, describe how it handles two practical
issues that the textbook recipe glosses over: missing bins (low
mappability, no data) and chromosome-arm sign conventions across
multiple chromosomes.


== Further Reading <sec:further-reading>

- *Frommer, M., McDonald, L. E., Millar, D. S., et al.* (1992).
  "A Genomic Sequencing Protocol That Yields a Positive Display of
  5-Methylcytosine Residues in Individual DNA Strands." _PNAS_ 89:
  1827–1831. The original bisulfite-sequencing paper.
- *Krueger, F., & Andrews, S. R.* (2011). "Bismark: A Flexible
  Aligner and Methylation Caller for Bisulfite-Seq Applications."
  _Bioinformatics_ 27: 1571–1572. The three-letter alignment paper
  every methylation pipeline cites.
- *Feng, H., Conneely, K. N., & Wu, H.* (2014). "A Bayesian
  Hierarchical Model to Detect Differentially Methylated Loci from
  Single Nucleotide Resolution Sequencing Data." _Nucleic Acids
  Research_ 42: e69. The DSS paper.
- *Lieberman-Aiden, E., van Berkum, N. L., Williams, L., et al.*
  (2009). "Comprehensive Mapping of Long-Range Interactions Reveals
  Folding Principles of the Human Genome." _Science_ 326: 289–293.
  The Hi-C paper that introduced the contact matrix, ICE, and the
  first-eigenvector compartment definition.
- *Imakaev, M., Fudenberg, G., McCord, R. P., et al.* (2012).
  "Iterative Correction of Hi-C Data Reveals Hallmarks of Chromosome
  Organization." _Nature Methods_ 9: 999–1003. The ICE paper.
- *Dixon, J. R., Selvaraj, S., Yue, F., et al.* (2012). "Topological
  Domains in Mammalian Genomes Identified by Analysis of Chromatin
  Interactions." _Nature_ 485: 376–380. TADs.
- *Rao, S. S. P., Huntley, M. H., Durand, N. C., et al.* (2014).
  "A 3D Map of the Human Genome at Kilobase Resolution Reveals
  Principles of Chromatin Looping." _Cell_ 159: 1665–1680. The
  high-resolution Hi-C paper plus HICCUPS and Arrowhead.
- *Crane, E., Bian, Q., McCord, R. P., et al.* (2015).
  "Condensin-Driven Remodelling of X Chromosome Topology During
  Dosage Compensation." _Nature_ 523: 240–244. The insulation-score
  paper.
- *Abdennur, N., & Mirny, L. A.* (2020). "Cooler: Scalable Storage
  for Hi-C Data and Other Genomically Labeled Arrays."
  _Bioinformatics_ 36: 311–316. The modern Hi-C storage format.
- *Cooler documentation.* `cooler.readthedocs.io`. The reference
  for working with modern Hi-C data structures.
- *Bismark user guide.* `felixkrueger.github.io/Bismark`. The
  reference for end-to-end bisulfite alignment.
