#import "../theme/book-theme.typ": *

= CRISPR Functional Screens and DepMap <ch:crispr-screens>

#matters[
  By the end of Chapter 4 you could read variants out of a genome. By
  the end of Chapter 13 you could associate them with phenotypes at
  biobank scale. Neither chapter told you what those variants *do* in
  a cell, and neither could tell you which of the genes they touch is
  worth the next twenty years of medicinal chemistry. CRISPR functional
  screens are how the field answers both questions experimentally and
  at genome scale. A single Cas9 ribonucleoprotein and a library of
  eighty thousand short RNAs are enough to switch every protein-coding
  gene off in a population of cells and read out which losses matter
  for fitness, drug resistance, or any other selectable phenotype. The
  DepMap project industrialised that recipe across roughly a thousand
  cancer cell lines, and most cancer drugs entering clinical trials
  after 2020 trace back to a dependency it surfaced. Multiplexed assays
  of variant effects do the same job one residue at a time, and they
  feed straight into the clinical variant-classification machinery you
  will meet in Chapter 25. This chapter walks the engineering — guide
  design, infection, sequencing, robust rank aggregation, copy-number
  correction — and pulls out the compressed-sensing structure that
  makes pooled screens work in the first place.
]

A pooled CRISPR screen is, mechanically, five operations. Build a
library of about eighty thousand single-guide RNAs covering every
protein-coding gene at roughly four guides apiece. Package the library
into lentivirus and infect a population of cells at a multiplicity of
infection low enough that almost every infected cell carries exactly
one guide. Split the culture into a sample taken at the start of the
selection and one or more samples taken after fourteen to twenty-one
days of growth, drug pressure, viral challenge, or whatever phenotype
the screen is built around. Amplify the integrated guide cassette by
PCR from genomic DNA and sequence it to a depth of about five hundred
reads per guide. Aggregate guide-level counts into gene-level scores
with a statistical solver. The output is a ranked list of genes whose
perturbation moved the cells through the selection.

Read at the right level of abstraction, those five operations are an
inverse problem. The vector of guide-level log fold changes you measure
is a noisy linear function of an underlying vector of per-gene effects.
The design matrix that connects them — guides to genes — is sparse,
with four non-zero entries per column, and that sparsity is exactly
what makes the system identifiable from far fewer measurements than
the naive count of unknowns would suggest. Compressed sensing tells
you so. The rest of the chapter is the working-out of that observation:
how the guides are picked, how the matrix is inverted, what goes wrong
when the sparsity assumption breaks, and what you get on the other end
once a few cell lines and a few drug arms have been thrown into the
mix.

The remaining chapters of the book use those outputs without
re-deriving them. Chapter 25 picks up MAVE-derived variant scores as
PS3 evidence for ACMG/AMP classification. Chapter 26 reaches DepMap
for drug-target prioritisation when discussing translational pipelines.
The compressed-sensing framing recurs in Chapter 16's pattern of
matched-prior architectures and in Chapter 21's discussion of pooled
metagenomics. This chapter is the spine.


== Cas9, Guide RNAs, and the Loss-of-Function Move <sec:mechanism>

The natural history of CRISPR is now well known, but the parts that
matter for screen engineering are narrow. CRISPR-Cas9 systems are an
adaptive immune mechanism in prokaryotes — bacteria and archaea cache
fragments of past phage genomes in their own DNA and use them as a
guide library to recognise and cut re-infecting viruses. *Jennifer
Doudna*, *Emmanuelle Charpentier*, and their groups demonstrated in
August 2012 (Jinek et al., _Science_ 337) that the
_Streptococcus pyogenes_ Cas9 protein could be reprogrammed in vitro
with a synthetic guide RNA. *Feng Zhang*'s group at the Broad
Institute (Cong et al., _Science_ 2013) and *George Church*'s group at
Harvard (Mali et al., _Science_ 2013) showed within six months that
the same system worked in human cells. The reagent and the recipe were
both there by mid-2013. The 2014 _Science_ papers from the *Shalem*
and *Wang* labs took the next step and showed that a *library* of
guides delivered by lentivirus could screen tens of thousands of genes
at once. The 2020 Nobel Prize cited Doudna and Charpentier for the
biochemistry; the screening literature traces from Shalem and Wang.

For screens, all you need to know about the molecule is what determines
whether a given guide cuts where you point it.

A *single-guide RNA*, abbreviated *sgRNA*, is about a hundred
nucleotides long. Twenty of those nucleotides — the *spacer* — are
the programmable region that base-pairs with the genomic target. The
remaining sequence is structural scaffold that holds the guide in the
Cas9 active site. The protein itself is a roughly 1,400-amino-acid
nuclease with two cleavage domains, *HNH* and *RuvC*, that cut the
two strands of the target DNA. Cleavage requires an additional
recognition motif on the DNA immediately downstream of the
spacer-matching sequence — the *protospacer-adjacent motif* or *PAM*.
For _S. pyogenes_ Cas9, the PAM is `NGG`. The cut is blunt and lands
exactly three base pairs upstream of the PAM. @fig:cas9 lays out the
geometry.

#figure(
  image("../figures/ch24/f1-cas9-sgrna-mechanism.svg", width: 95%),
  caption: [
    Cas9 with a single-guide RNA at a target site. Twenty nucleotides
    of spacer base-pair with the genomic target; an `NGG` PAM licenses
    cleavage three bases upstream. The double-strand break is repaired
    by non-homologous end joining, which is short and error-prone, so
    most edited alleles end up with frameshifting indels and a
    truncated, non-functional protein.
  ],
) <fig:cas9>

The reason all of this matters for a knockout screen is that mammalian
cells repair double-strand breaks with a strong bias toward
*non-homologous end joining*, NHEJ. NHEJ is fast and inaccurate. It
glues the cut ends back together with frequent small insertions and
deletions of one to ten bases. Roughly eighty-five to ninety-five
percent of edited alleles end up with an indel that is not a multiple
of three. The mRNA produced from that allele has a shifted reading
frame, hits a premature stop codon within a hundred bases, and gets
degraded by nonsense-mediated decay or translated into a truncated
non-functional protein. A cell whose two chromosomal copies of a gene
have both received frameshifting indels is, for practical purposes, a
knockout. Bulk-population screens average over this — most cells in
the population are knockouts at the gene the guide targets, and the
phenotype is the average behaviour of the knockout subpopulation.

#note[
  The "knockout" produced by CRISPR is statistical, not surgical. In
  any given cell, both alleles need to receive a frame-shifting edit
  before the gene is genuinely lost. At guide activities you would
  call competent — Doench Rule Set 2 score above 0.6, say — the
  fraction of cells with bi-allelic knockout in a diploid human line is
  typically 60 — 90 percent after a fortnight in culture. Pooled screens
  read out the population-mean phenotype across that mixture; the
  signal is therefore attenuated by the fraction of un-edited cells.
  Modern libraries pick four guides per gene partly to insure against
  the unlucky-low-editing-efficiency outcome at any one guide.
]

Three engineering variants of the same backbone are now standard.
*CRISPR knockout* (CRISPRko) is what we have just described — wild-type
Cas9 cuts and NHEJ scrambles the locus. *CRISPR interference* (CRISPRi),
introduced by *Luke Gilbert* and *Stanley Qi* working with *Jonathan
Weissman* in 2013, replaces the nuclease activity with a catalytically
dead Cas9 (dCas9) fused to a KRAB transcriptional repressor. The
fusion sits on the promoter region of the targeted gene and silences
transcription without making a permanent edit. Knockdown rather than
knockout. The signal is partial — typically 70 — 95 percent reduction
in mRNA — and reversible, which is the point: essential genes survive
the screen long enough to be seen as dropping out, because partial
knockdown does not kill cells instantly the way a clean knockout does.
*CRISPR activation* (CRISPRa), introduced by *Silvana Konermann* and
*Feng Zhang* in 2015, takes the opposite tack: dCas9 fused to a VP64
activator domain (or the more potent SunTag and SAM systems) drives
expression of normally silent genes. Gain-of-function screens read out
the opposite direction — which forced overexpressions confer
resistance to a drug, which drive a developmental transition. The
three modalities are complementary measurements of the same gene from
different sides: knockout removes it, interference reduces it,
activation amplifies it. Modern projects often run all three on the
same gene set and triangulate.

Two further classes — *base editors* (Komor et al., 2016; Gaudelli et
al., 2017) and *prime editors* (Anzalone et al., 2019) — turn the
machinery into precision instruments rather than knockout factories.
A base editor fuses a deaminase to a catalytically attenuated Cas9
and converts a single nucleotide (C → T for cytidine deaminases, A → G
for adenine deaminases) inside a narrow editing window without making
a double-strand break. A prime editor uses a reverse-transcriptase
fusion and an extended guide that templates an arbitrary edit. Both
matter for the MAVE half of this chapter because they let a screen
introduce all possible single-base variants in a gene rather than
knock the gene out entirely.

#tip[
  When you are picking a modality, the question is not which is
  fanciest. It is what kind of perturbation the biology you want to
  read out is sensitive to. For an essentiality dropout screen in
  cancer cell lines, CRISPRko is the default because clean loss is
  what kills the cells. For surveying which genes you can drug to
  treat an essential-function disease, CRISPRi gives you the partial
  knockdown that approximates a small-molecule inhibitor. For target
  identification in gain-of-function biology — drug resistance,
  developmental driver genes, oncogenic transformations — CRISPRa is
  the only modality that asks the right question. The three are not
  interchangeable.
]


== From a Library to a Read Count <sec:library-to-counts>

A screen is engineered from the library outward. A *library* is a pool
of plasmids, each carrying one sgRNA expression cassette and a
selectable marker. A typical human CRISPRko library has roughly four
guides per gene across the 19,000-odd protein-coding genes for a total
near 76,000 — *Brunello* (Doench et al., 2016) is the canonical
example, *TKOv3* (Mair et al., 2019) is the matched essentiality-screen
companion, and *GeCKO v2* (Sanjana et al., 2014) is the older legacy
library with six guides per gene at the cost of more low-activity ones.
For CRISPRi the corresponding choices are *Dolcetto* and *Horlbeck*; for
CRISPRa, *Calabrese* and *Caprano*.

The libraries are not interchangeable. Two metrics distinguish them.
*On-target activity* is the average per-guide cleavage efficiency at
the intended site, and modern libraries select for high-activity
guides using the *Doench Rule Set 2* score — a linear regression on
sequence features (position-specific nucleotides, GC content, PAM
context) trained on a set of guides whose cutting efficiency was
measured directly. Azimuth 2.0 (Doench, 2016) and the deep-learning
extensions that followed (DeepCRISPR, CRISPR-Net) report higher
predicted activity at the high-score end of the distribution; in
practice they buy you a few percentage points of recall on dropout
screens, not the order-of-magnitude differences seen in other ML
domains. *Off-target burden* is the second axis: how often does a
guide cut somewhere it should not? The *CFD score* (Cutting Frequency
Determination, Doench 2016) and the older *MIT specificity score*
estimate per-guide off-target activity by summing penalised matches
across the genome. A high-quality library combines high on-target
activity with low off-target burden and is uniformly tiled across the
gene's coding sequence so that the four guides are not all in one
exon — pick a guide too close to the C-terminus and a knockout-grade
indel can produce a near-full-length protein with intact function.

The choice of library is the screen's first inductive bias. Pick a
library with low-activity guides and you will see weak signal because
many of your knockouts never actually happen. Pick a library with high
off-target burden and you will see strong signal where you do not
want it — on the bystander gene with the homologous protospacer, not
the gene you targeted. Brunello and TKOv3 both fix this, in different
ways; legacy libraries do not.

#warn[
  Guide-design tools are not interchangeable across organisms or
  Cas9 orthologues. Doench Rule Set 2 was trained on _S. pyogenes_
  Cas9 (`NGG` PAM) data from human and mouse cells; its predictions
  degrade outside that domain. If you are running an _S. aureus_ Cas9
  screen with an `NNGRRT` PAM, switching to its purpose-trained
  model (Najm 2018, _Nature Biotechnology_) or to a more recent
  pan-PAM predictor like DeepCRISPR is the correct move. The same
  caveat applies to base-editor libraries — they need their own
  scoring system because the active window is two-to-four bases
  rather than the full twenty.
]

With a library picked, the experimental side is a sequence of yields
and bottlenecks.

*Lentiviral packaging.* The plasmid pool is co-transfected with viral
packaging genes into HEK293T cells, which produce virions carrying the
guide cassettes. The titre — virions per millilitre — sets how many
cells can be infected per batch.

*Transduction at low MOI.* Cells are exposed to the lentivirus at a
*multiplicity of infection* (MOI) of about 0.3 — three virions per ten
cells, on average. At MOI = 0.3 the Poisson distribution
$P(k) = e^(-0.3) 0.3^k / k!$ tells you that about 26 percent of cells
get exactly one virion and only 3 percent get two or more. A higher
MOI gives more infected cells but starts mixing perturbations — a cell
with two guides has a confounded phenotype. A lower MOI keeps the
mixing rate down but throws away cells. The 0.3 sweet spot is so
standard that the field treats it as a constant.

*Cell coverage.* For each guide to have enough independent cellular
founders to average over, you need about five hundred cells per guide
at the moment of infection. With seventy-five thousand guides and a
30 percent infection rate, that means a starting population of about
$75{,}000 times 500 / 0.3 approx 1.25 times 10^8$ cells per replicate.
Two or three replicates per condition.

*Selection.* Apply the phenotype that will separate the guides — a
drug, a virus, gravity (literally, in some screens), continued
proliferation under no selection. For dropout screens the selection
is just time; let the cells grow for fourteen to twenty-one days, and
guides targeting essential genes lose representation as their cellular
hosts fail to divide.

*Genomic-DNA extraction and PCR amplification.* The integrated guide
cassettes are amplified out of the cellular genome with primers
flanking the spacer region. The amplicon is short — about 200 — 300 bp
— so single-end Illumina reads at 75 — 150 bp are enough.

*Sequencing.* Read each amplicon to a depth of about five hundred
reads per guide. For a 75,000-guide library that is $3.75 times 10^7$
reads per sample, or about $2 times 10^8$ reads for a six-sample screen.
A single NovaSeq lane easily covers this.

*Counting.* Map each read's spacer region back to the library by
exact match (or one-mismatch). The output is a *count table*: rows
are guides, columns are samples, cells are read counts. For a typical
screen the table is 75,000 by six-to-eight and weighs a few megabytes.

The count table is everything that downstream analysis has to work
with. Everything past this point is statistics applied to a sparse
non-negative integer matrix.

#figure(
  image("../../diagrams/lecture-24/01-pooled-screen.svg", width: 92%),
  caption: [
    The pooled-screen workflow. A lentiviral library is delivered to a
    cell population at MOI ≈ 0.3; cells are split into T0 and T-final
    samples; integrated guide cassettes are amplified by PCR from
    genomic DNA and sequenced; counts are normalised and aggregated to
    gene-level fitness scores. A 75 k-guide screen in two replicates
    of two conditions runs at roughly USD 20,000 of consumables in 2024.
  ],
) <fig:pooled>


== The Screen as a Sparse Linear Inverse Problem <sec:linear-system>

Stop and re-read the count table the way an EE would read it. A guide
has a count $c_(g, "pre")$ in the pre-selection sample and a count
$c_(g, "post")$ in the post-selection sample. The log fold-change is

$ y_g = log_2 ((c_(g, "post") + 1) / (c_(g, "pre") + 1)). $

Stack the $y_g$ for all guides into a vector $bold(y) in RR^M$ with
$M approx 80{,}000$. Each guide targets one gene out of $N approx
20{,}000$. Define the design matrix $A in {0, 1}^(M times N)$ with
$A_(g, j) = 1$ if guide $g$ targets gene $j$ and $0$ otherwise. The
matrix is sparse — exactly one non-zero per row, four non-zero
entries per column for a Brunello-style library. The underlying
quantity you would like to recover is $bold(x) in RR^N$, the
log-scale effect of each gene's knockout on fitness.

If the only noise were a per-guide measurement error $bold(epsilon)$,
the system reads

$ bold(y) = A bold(x) + bold(epsilon). $

@fig:linear shows the geometry. The system has more equations than
unknowns by a factor of four, so least squares would solve it; the
real difficulty is that $bold(epsilon)$ is heavy-tailed (PCR jackpotting,
guide-activity variation, copy-number-driven outliers) and that not
every entry of $A$ that is set to one represents the same effective
perturbation strength.

#figure(
  image("../figures/ch24/f2-screen-as-linear-system.svg", width: 96%),
  caption: [
    A pooled screen as the linear system $bold(y) = A bold(x) +
    bold(epsilon)$. The design matrix $A$ has roughly eighty thousand
    rows (guides), twenty thousand columns (genes), and four non-zero
    entries per column. The sparsity and the redundancy are what make
    the system identifiable; the solver of choice depends on what
    distribution you put on the noise.
  ],
) <fig:linear>

This is the *compressed sensing* analogy. Donoho (2006) and
Candès-Romberg-Tao (2006) formalised the conditions under which a
sparse vector can be recovered from a small number of incoherent
linear measurements. The standard pooled-screen design is not quite
that case — $bold(x)$ itself is not sparse, since most genes have
some non-zero fitness effect — but the structural insight is the
same. Four redundant measurements per gene plus a fast statistical
test recovers gene effects far more reliably than a naive
per-measurement test would. The four-per-gene redundancy is the
matched-filter design choice that, combined with a robust aggregation
rule, makes the screen work.

#note[
  The compressed-sensing reading is more than ornamental. Pooled
  combinatorial screens — where each cell gets two or more guides at
  once — are precisely a compressed-sensing setup where $bold(y)$ is
  now under-determined (each cell senses a linear combination of
  perturbations) and $bold(x)$ is genuinely sparse (most gene
  combinations have no effect). The 2019 _Compressed Perturb-seq_
  paper from Adamson, Norman, and Weissman makes this explicit, with
  hundreds of perturbations recovered from a number of cells far
  below the naive one-cell-per-condition count, by exploiting
  $ell_1$-regularised sparse recovery in the same family as LASSO.
]

The textbook MAGeCK / DrugZ / BAGEL2 pipelines do not invoke the
compressed-sensing vocabulary, but every one of them is a particular
choice of noise model and prior on the recovery problem. MAGeCK uses
a rank-based, non-parametric test that is robust to heavy tails.
DrugZ Z-scores guides against a non-targeting control distribution
and is well-tuned for drug-vs-vehicle resistance screens. BAGEL2 puts
a Bayesian prior — informative core-essential and non-essential gene
sets — on top of the likelihood and computes per-gene log-Bayes-factors.
CERES and Chronos use a least-squares model with explicit corrections
for copy-number bias and per-guide activity, and they handle the
DepMap-scale problem of pooling across many cell lines simultaneously.
Pick the prior that matches the noise you actually have.


== MAGeCK and Robust Rank Aggregation <sec:mageck>

*MAGeCK*, the Model-based Analysis of Genome-wide CRISPR-Cas9
Knockout, was published by *Wei Li* and *Han Xu* in *Xiaole Shirley
Liu*'s lab in 2014 (_Genome Biology_ 15: 554), and it is the
analytical default for pooled screens. The recipe runs in three
stages.

*Stage one — guide-level testing.* For each guide $g$ and each
condition, MAGeCK fits a *negative binomial* model to the count
distribution. The negative-binomial is the same workhorse that
Chapter 6 used for bulk RNA-seq differential expression: it handles
the over-dispersion that count data carries when biological
variability adds variance on top of pure Poisson sampling. The
mean-variance relationship is

$ "Var"(c_g) = mu_g + alpha_g mu_g^2, $

where $alpha_g$ is the per-guide dispersion. MAGeCK estimates
$alpha_g$ by pooling across guides with similar mean counts (the same
shrinkage idea Chapter 6 traced through DESeq2). Per-guide one-sided
$p$-values are computed by tail probability under the negative
binomial fit.

*Stage two — normalisation.* Per-sample size factors correct for
library-prep depth differences. MAGeCK's default is *median
normalisation* — divide each sample's counts by the sample's median,
re-scaled so the geometric mean across guides is preserved. For
high-quality libraries the result is statistically indistinguishable
from DESeq's geometric-mean normalisation. Some screens benefit from
*median normalisation over non-targeting controls only* — the field
has converged on this for drug screens because non-targeting guides
should have, by construction, zero true fold-change.

*Stage three — gene-level aggregation by modified RRA.* The interesting
piece. For each gene with $k$ guides, MAGeCK computes the $k$
guide-level $p$-values, ranks them across the genome to produce $k$
rank-percentiles $r_1, dots, r_k in [0, 1]$, and reports the gene's
score as

$ rho_(j) = min_(i = 1, dots, k) "Beta"(r_((i)); a = i, b = k - i + 1), $

where $r_((1)) <= r_((2)) <= dots <= r_((k))$ are the order statistics
and the Beta CDF gives the tail probability of seeing the $i$-th
smallest of $k$ uniform draws below $r_((i))$ under the null. The
minimum over $i$ — i.e. the most-extreme order statistic relative to
its own null — is the gene's RRA statistic. Significance is computed
either analytically (Beta tail) or by permutation.

In one line: a gene with four guides whose $p$-values are uniformly
small is unlikely under the null; a gene with one strong guide and
three weak ones is also unlikely if the one strong guide is strong
*enough*; a gene with four random guides has a uniform $r_((i))$
distribution and a Beta CDF spread evenly over $[0, 1]$. The
*modified* in modified-RRA is that MAGeCK uses the median rank rather
than the strict minimum to absorb one outlier guide — a single
non-functional guide does not kill the call. @fig:rra walks the
shape.

#figure(
  image("../figures/ch24/f3-rra-beta-null.svg", width: 95%),
  caption: [
    Robust rank aggregation in MAGeCK. The left panel shows three
    genes whose $k = 4$ guide-level $p$-value ranks differ
    qualitatively — concordant in the low tail (pan-essential),
    three-of-four with one outlier (context-essential), uniform
    (null). The right panel is the Beta(1, $k$) density that
    MAGeCK's analytic significance test uses; a gene's $p$-value is
    the tail probability of seeing the observed minimum rank under
    that null.
  ],
) <fig:rra>

#note[
  The non-parametric move at the heart of RRA is the same one that
  drives *gene set enrichment analysis*. Both algorithms ask whether a
  set of measurements occupies the extreme tail of a ranked
  distribution more than chance would predict. The order-statistic
  theory is identical; the difference is bookkeeping. Wherever you
  have a fixed-size set of related measurements and a ranking, this
  family of tests is a uniformly minimum-variance unbiased way to
  pool them.
]

MAGeCK's output is a per-gene log fold-change estimate, a $p$-value,
and a Benjamini-Hochberg-corrected FDR. For a typical dropout screen
on a cancer cell line you would expect roughly 1,500 — 2,500 genes at
FDR $< 0.05$ — about half of these are the *core essentials* that
drop out of every cancer cell line, and the rest are cell-line- or
context-specific. The interesting biology is in the second category.

*DrugZ* (Wang et al., 2017) is the standard tool for drug screens —
where you compare guide abundance in a drug-treated arm against a
vehicle-treated control. It computes guide-level Z-scores using the
control as the null distribution and aggregates them across guides by
summing into a gene-level Z. Its strength is that it does not assume
the bulk of guides have zero effect; for a drug screen where many
gene knockouts genuinely shift the phenotype, that matters.

*BAGEL2* (Kim et al., 2019) is the Bayesian classifier in the family.
It takes two reference gene sets — the *core essential* (CEG)
collection (about 600 genes that drop out of nearly every screen,
curated by *Traver Hart*) and a non-essential set of similar size —
and computes per-gene log-Bayes factors for the hypothesis that the
gene behaves like a CEG. Output is a single Bayes-factor score per
gene; high scores are essential, low scores are not. BAGEL2 is the
right tool when the goal is essentiality calling — yes-or-no on the
core-essential category — rather than fold-change ranking.

The 2020s consensus practice is to run all three on the same data and
take the union of significant hits as a top-line list and the
intersection as a high-confidence list. Agreement is generally good
across the strongly-essential and strongly-resistant tails; the
middle of the distribution is where the tools diverge and where
careful interpretation is required.

#figure(
  image("../../diagrams/lecture-24/03-mageck.svg", width: 92%),
  caption: [
    MAGeCK's pipeline. Count tables → median normalisation →
    per-guide negative-binomial test → per-gene modified RRA →
    Benjamini-Hochberg FDR → ranked hit list. The same backbone
    drives DrugZ and BAGEL2 with different choices of guide-level
    test statistic and aggregation rule.
  ],
) <fig:mageck>

Quality control before interpretation is non-negotiable. Four
metrics decide whether a screen is interpretable at all.

*Guide recovery* — the fraction of library guides detected with five
or more reads at T0 — should exceed 90 percent. A drop below that
implies a bottleneck somewhere in lentiviral packaging or PCR
amplification, and the resulting analysis will under-represent the
true diversity of the population.

*Replicate concordance* — Pearson correlation of guide-level log
counts between biological replicates of the same condition — should
exceed 0.9 in the post-selection samples. Lower correlations suggest
batch effects or insufficient cell coverage.

*Core-essential dropout* — the fraction of Hart's CEG list showing
significant dropout at T-final — should exceed 80 percent for a
healthy dropout screen. Lower numbers indicate a screen that is
under-powered or that has not been carried long enough.

*Non-essential null distribution* — the log-fold-change distribution
of *non-targeting controls* (guides with no genomic target) should be
centred on zero with reasonable spread. Skew or shifted centre flags
normalisation problems.

If any of these fails, re-run the screen. The temptation to bandage
QC failures with statistical post-hoc corrections is strong and
almost always wrong; a screen with 70 percent CEG recovery is not
fixable downstream.

#warn[
  Copy-number bias is the screen-killer most-likely to bite you in
  cancer cell lines. Cas9 cuts every copy of its target; in a
  high-copy-number locus, that means many simultaneous double-strand
  breaks, which produces a cell-fitness penalty unrelated to the
  function of the targeted gene. *CRISPRcleanR* (Iorio et al., 2018)
  and the CERES correction (Meyers et al., 2017) handle this by
  modelling guide log fold-change as a function of regional
  copy-number and subtracting the predicted bias. Skip the correction
  and a HeLa-scale screen will return amplification regions as the
  most essential parts of the genome — a confidently false answer.
]


== A Worked Example <sec:worked-example>

The reader-facing way to make all of this concrete is to walk one
specific screen. Take *TKOv3* (Mair et al., 2019), a 71,000-guide
library with four guides per gene optimised for dropout screens, and
run it in *A375*, a melanoma cell line, against a vehicle control,
for 21 days, with two biological replicates. Sequence to about 400
reads per guide per sample. Run MAGeCK.

The QC pass: guide recovery 96 percent, replicate Pearson 0.93,
non-targeting control log-fold-change distribution centred on −0.02
with standard deviation 0.41, core-essential dropout 88 percent. Pass.

The gene-level output ranks about 1,800 genes at FDR < 0.05 with
negative log fold-change. The top of the list is dominated by
ribosomal subunits (`RPS6`, `RPL5`, `RPL11`), proteasome components
(`PSMA3`, `PSMB5`), translation initiation factors (`EIF3A`), and
mitochondrial complex-V subunits — the universal essentials. None of
this is interesting on its own; every cancer cell line returns the
same core list, which is a sanity check that the screen worked.

The interesting biology starts where A375 differs from the average
cell line. `BRAF` and `MAPK1` show strong dropout because A375 is
`BRAF`-V600E mutant and depends on the MAPK pathway. `MITF` shows
dropout because A375 is melanocyte-lineage and depends on the
melanocyte transcription factor. `MDM2` shows dropout because A375
has wild-type `TP53` and is therefore exquisitely dependent on the
`MDM2` ubiquitin ligase to keep p53 levels in check. None of these
are surprises in 2024; they are the validation that the screen sees
real lineage biology. The genuinely-unknown signal is in the next two
hundred genes, where context-specific essentialities lurk that no one
has characterised yet.

The positive-selection arm — same screen, with a BRAF inhibitor at
sub-lethal dose for the same 21 days — produces a complementary
output. Guides targeting `NF1`, `PTEN`, and a handful of MAPK-pathway
phosphatases come up *enriched* (positive log fold-change), meaning
cells whose knockout removes a negative regulator of MAPK signalling
acquire a survival advantage under BRAF inhibition. This is the
2014 — 2016 mechanism-of-resistance story for melanoma, replayed in a
single screen.

#tip[
  A useful sanity check whenever you run a new screen is to scan the
  top thirty hits for known biology. If a `BRAF`-mutant melanoma
  screen does not surface `BRAF` and `MITF`, something is wrong with
  the screen, the cell line, or the analysis. The universally
  essential genes plus a small number of cell-line-specific drivers
  are the screen's known-answer key.
]


== DepMap and Cancer Dependencies <sec:depmap>

The single screen above is a snapshot of one cell line. The DepMap
project is the same screen, run in roughly a thousand cancer cell
lines covering thirty cancer types, over the course of seven years.

*DepMap* is a Broad Institute and Sanger Institute collaboration that
emerged from the 2017 paper "Defining a Cancer Dependency Map" by
*Aviad Tsherniak* and colleagues (_Cell_ 170: 564). The premise is
simple. The Cancer Cell Line Encyclopedia (CCLE) and the Sanger
COSMIC Cell Lines Project together hold characterised, omics-profiled
cancer cell lines numbered in the low thousands. Run a genome-scale
CRISPR knockout screen on every one of them, cross-reference against
the matched RNA-seq, methylation, mutation, and copy-number data, and
you have a four-dimensional map: for every gene, for every cell line,
how essential is its knockout, and what about the cell line's
molecular profile predicts that essentiality? Released publicly
through *depmap.org* at quarterly intervals, the project has by 2024
sequenced over 1,100 cell lines and runs at the centre of pharma
target-discovery pipelines.

The scoring algorithm has gone through two generations. *CERES* (Meyers
et al., 2017) was the original — least-squares with corrections for
copy-number bias and per-guide activity. *Chronos* (Dempster et al.,
2021) is the current version. Chronos fits a time-resolved
exponential-growth model to the count data, treating the screen as a
dynamical system rather than a two-time-point regression, and outputs
a per-cell-line, per-gene *gene effect* score. By convention a
gene-effect score of $-1$ corresponds to the median of the
core-essential set, and $0$ to the median of the non-essential set;
the scale is calibrated to make per-cell-line comparisons direct.

@fig:depmap shows what the resulting matrix looks like. Three classes
of column emerge.

#figure(
  image("../../diagrams/lecture-24/04-depmap.svg", width: 92%),
  caption: [
    DepMap dependency landscape. Pan-essentials (solid blue columns —
    essential everywhere) form a "no-go" set for drug targeting,
    pan-non-essentials (white columns) are dispensable everywhere, and
    *lineage-selective* dependencies (banded blue) are the drug-target
    gold.
  ],
) <fig:depmap>

*Pan-essentials* are genes whose loss kills nearly every cell line —
the ribosomal subunits, the proteasome, mitochondrial complexes. These
are the bottom of the universal-essentials list. From a drug-discovery
standpoint they are useless. Killing the ribosome would kill the
patient as fast as it killed the tumour.

*Pan-non-essentials* are genes that nothing depends on under any
condition. Most of the genome, by gene count. Also uninteresting for
drug targeting.

*Lineage-selective* and *context-selective* dependencies are the
interesting third class — genes whose loss kills one cancer type but
not others, or one mutational background but not others. These are
what the project was built to find. By 2024 about fifty such
dependencies have entered active clinical development; about ten have
reached late-stage trials and three or four (depending on how you
count partial successes) have produced FDA-approved drugs.

The canonical examples are worth naming because they reappear
throughout the chapter. *MCL1* dependency in haematological cancers
versus solid tumours drives the MCL1-inhibitor programmes at AstraZeneca
and AbbVie. *EZH2* dependency in lymphomas carrying gain-of-function
EZH2 mutations led to Tazemetostat (Tazverik, FDA-approved 2020). *WRN
helicase* dependency in microsatellite-instability-high (MSI-H)
cancers — a 2019 _Nature_ finding from three groups simultaneously
(Chan et al., Kategaya et al., Behan et al.) — is now in Phase II with
multiple compounds in development.

The most-cited DepMap-driven discovery is the *MTAP / PRMT5* synthetic
lethality. The *MTAP* gene sits adjacent to the tumour-suppressor
*CDKN2A* on chromosome 9p21 — a region deleted in roughly 15 percent
of all cancers as collateral to the CDKN2A loss. MTAP encodes
methylthioadenosine phosphorylase, an enzyme in the methionine
salvage pathway. MTAP loss causes accumulation of the substrate MTA,
which partially inhibits the arginine methyltransferase PRMT5. The
DepMap signal — strong PRMT5 essentiality specifically in MTAP-deleted
cell lines — was first reported by Marjon, Mavrakis, and colleagues in
2016 and confirmed across the expanding DepMap corpus. MTA-cooperative
PRMT5 inhibitors (MRTX1719, AMG 193, AG-270) are now in Phase I/II in
patients whose tumours carry MTAP deletion. @fig:synth-lethal walks
the logic.

#figure(
  image("../figures/ch24/f4-synthetic-lethality-logic.svg", width: 96%),
  caption: [
    Synthetic lethality logic with two clinical case studies. The 2 × 2
    truth table on the left is the abstract definition: both single
    knockouts are viable, the double knockout is lethal. The BRCA /
    PARP case study established the framework in the clinic (olaparib,
    FDA 2014); the MTAP / PRMT5 case study, anchored in DepMap, brings
    the discovery loop full circle.
  ],
) <fig:synth-lethal>

#note[
  The conceptual ancestor of all of this is *Leland Hartwell*'s 1997
  paper "Integrating genetic approaches into the discovery of
  anticancer drugs" (_Science_ 278: 1064) and *William Kaelin*'s 2005
  formalisation of the synthetic-lethality concept for cancer therapy
  (_Nature Reviews Cancer_ 5: 689). Both were written before CRISPR
  existed; both predicted exactly the kind of systematic discovery
  pipeline that DepMap eventually became. The first clinical
  realisation, the *PARP* inhibitor *olaparib* approved by the FDA in
  2014 for BRCA1/2-deficient ovarian cancer, came out of the Bryant
  (2005) and Farmer (2005) yeast-screen tradition — pre-CRISPR but
  philosophically continuous with what DepMap does at scale.
]

DepMap also incorporates a second layer of data: *drug screens*.
*PRISM* (Profiling Relative Inhibition Simultaneously in Mixtures,
Yu et al., 2016) is a barcode-multiplexed small-molecule screen that
runs about 6,000 compounds across 800 cell lines simultaneously, with
each cell line carrying a unique barcode that allows pooled
phenotyping. *GDSC* (Genomics of Drug Sensitivity in Cancer) and
*CTRPv2* (Cancer Therapeutics Response Portal) provide complementary
drug-response data at narrower depth. Cross-referencing CRISPR
dependencies with drug-response patterns is a hypothesis generator
for mechanism of action: a drug that kills the same cell lines whose
gene-$X$ knockout kills is likely targeting gene $X$ or its pathway.
For drugs with no known target — a not-uncommon situation for
phenotype-derived compounds — this is the cleanest available
mechanism-discovery handle.


== MAVEs: Variant-Effect Maps at Saturation <sec:mave>

Single-gene functional measurement is a separate use of the same
technology and deserves its own treatment. A *multiplexed assay of
variant effects* — MAVE — measures the functional effect of every
possible variant in a defined region. The locus might be a protein
domain, a regulatory element, or a whole gene; the variants might be
single-amino-acid changes, single-base substitutions, or short indels.
The unifying recipe is the same as a screen: build a library of all
the variants, introduce them into a population of cells under a
phenotypic selection that depends on the locus, and read out how each
variant's frequency changes over the selection.

The methodology was named by *Doug Fowler* and *Stan Fields* in their
2014 _Nature Methods_ paper, but the experimental tradition runs back
to the 1980s biochemistry of cassette mutagenesis. What changed in the
2010s was the scale: synthesised libraries of tens of thousands of
variants, delivered to cells by CRISPR-mediated knock-in or by base
editing, and read out by deep amplicon sequencing of the edited locus.
Two flavours dominate.

*Deep mutational scanning (DMS)* is the protein-level variant. Take a
gene's coding sequence, synthesise a library in which every codon is
replaced by every other codon at one position, repeated across all
positions in the gene. Express the library in cells and apply a
selection that depends on the protein's function. Sequence the
surviving variants. Each amino-acid substitution gets a *function
score* — log of post-selection frequency over pre-selection frequency
— and the result is a map of size $L times 20$ for an $L$-residue
protein.

*Saturation genome editing* is the in-locus variant. Rather than
expressing variants from a plasmid, you edit them directly into the
chromosomal copy of the gene, in cells under selection for the gene's
endogenous function. *Greg Findlay* and colleagues' 2018 _Nature_
paper on *BRCA1* (Findlay et al., 562: 217) is the canonical example:
roughly 4,000 single-nucleotide variants across the BRCA1 RING and
BRCT functional domains, edited into HAP1 cells (a haploid human
cell line where a single-copy knockout suffices), grown for 11 days
under selection for BRCA1-dependent homologous recombination, and
read out by amplicon sequencing. The output is a per-variant function
score.

@fig:mave-clinical shows what comes next. The variant-level scores
are binned — non-functional, intermediate, functional — and the bins
feed directly into the *ACMG/AMP* clinical-variant-classification
framework that Chapter 25 will spend most of its time on. A variant
scored non-functional in a well-calibrated MAVE earns *PS3* evidence
("functional studies show damaging effect"); a variant
indistinguishable from wild-type earns *BS3* ("functional studies
show no damaging effect"). Findlay's BRCA1 paper reclassified roughly
96 percent of prior BRCA1 variants of uncertain significance with
high confidence — a number that should be seen as a transformative
addition to clinical molecular diagnostics rather than as an
incremental improvement.

#figure(
  image("../figures/ch24/f5-mave-clinical-bridge.svg", width: 96%),
  caption: [
    The MAVE-to-clinic bridge. A saturation genome editing experiment
    produces per-variant function scores; binned scores map onto
    ACMG/AMP PS3 / BS3 evidence codes; the resulting reclassification
    moves variants from the "uncertain significance" bucket into
    definite Pathogenic / Likely Pathogenic or Likely Benign / Benign
    calls. The illustrative numbers reflect the magnitude of the shift
    Findlay (2018) measured for BRCA1.
  ],
) <fig:mave-clinical>

The same approach has been applied to a growing catalogue of clinically
important genes. *MSH2*, *MLH1*, *MSH6*, *PMS2* — the Lynch-syndrome
mismatch-repair genes — have published saturation editing studies.
*TP53* has multiple DMS datasets covering both DNA-binding and
oligomerisation domains. *PTEN* has DMS for phosphatase activity and
for cellular stability separately. *MAVEdb* (mavedb.org) is the
community repository and aggregates published MAVE datasets across
roughly fifty genes as of mid-2024.

#figure(
  image("../../diagrams/lecture-24/05-mave-map.svg", width: 92%),
  caption: [
    A MAVE heat map. Each cell of the matrix shows the experimental
    function score of one single-amino-acid variant at one position.
    Functional sites — catalytic residues, binding interfaces — appear
    as vertical bands of damaging variants. Wild-type residues show as
    neutral cells along the diagonal of the map.
  ],
) <fig:mave-map>

The honest assessment is that MAVE coverage is still uneven. Roughly
fifty published genes is not the 19,000-gene scale that variant
interpretation would need to be fully covered. For unmeasured genes,
the computational complement — *AlphaMissense* (Cheng et al., 2023,
trained partly on MAVE data), *EVE* (Frazer et al., 2021), and the
older PolyPhen-2 and CADD — supplies *PP3* / *BP4* evidence at lower
strength. The 2024 frontier is integration: train a deep model on
existing MAVEs across a panel of genes, and use its cross-gene
generalisation to predict missing entries, then experimentally verify
the high-impact predictions. This is a working strategy, not a
production one, and the field is still calibrating what *PS3* evidence
strength should be assigned to model-predicted scores versus
experimentally-measured ones.

#warn[
  A MAVE score is only as good as the assay that produced it. The
  BRCA1 study reads out homologous-recombination function specifically;
  a variant scored "functional" in that assay might still disrupt a
  different BRCA1 activity (E3 ligase, transcriptional regulation)
  that the assay cannot see. The PS3 evidence weight should reflect
  the assay's coverage of the clinically-relevant function, not the
  assay's per-variant precision. *Calibration to clinical truth* —
  validating MAVE scores against known-pathogenic and known-benign
  control variants — is the standard practice in the field, and the
  community has converged on it for the existing MAVEdb panel.
]


== Pitfalls and the Reproducibility Story <sec:pitfalls>

The chapter has flagged most of the pitfalls in passing; this section
collects them in one place and adds two more.

The *library pitfalls* are the easiest to avoid: use a modern,
well-validated library (Brunello, TKOv3, Dolcetto, Calabrese), insist
on at least four guides per gene, and do not extend a legacy GeCKO v2
analysis past its expiry date.

The *experimental pitfalls* are MOI control, cell coverage, and
replication. MOI above 0.3 produces guide co-occurrence in single
cells and confounded phenotypes; correct it by titrating the virus on
a fresh batch of cells and counting puromycin-resistant colonies
before the screen. Cell coverage below 500 cells per guide adds
sampling noise on top of biology; the fix is more cells, not more
statistics. Replicate concordance below 0.9 typically signals batch
effects; investigate before interpretation, do not pool the offending
arms.

The *bioinformatic pitfalls* are copy-number bias, library plasmid
contamination, and the always-essential dominance.
*Copy-number-variation bias* in cancer cell lines is the most common
silent killer; CRISPRcleanR or CERES handles it and should be
default-on for cancer-screen analysis. *Library plasmid contamination*
arises when residual plasmid DNA from the lentiviral packaging step
gets co-amplified with genomic DNA; the symptom is a count
distribution at T0 that mirrors the plasmid distribution too closely
and a non-essential-gene null that is too narrow. Hash-check the count
distribution against the input library and re-run the protocol with
DpnI digestion of the gDNA preparation if the symptom appears. The
*always-essential dominance* problem is mathematical: a few thousand
core-essential genes drop out hard in every screen and dominate the
empirical $p$-value distribution if you do not standardise the
analysis to per-gene effects rather than raw fold-changes. MAGeCK's
rank-based test partially handles this; BAGEL2's CEG-prior handles it
explicitly.

The *statistical pitfalls* are multiple-testing, power, and
imbalanced designs. Twenty thousand genes tested simultaneously means
that Benjamini-Hochberg FDR control is essential; a Bonferroni
correction is too conservative and most screens use BH. Power
calculations should target an effect size based on the screen's
purpose — dropout screens want enough power for log fold-changes of
roughly $-1$ to be detectable, drug-resistance screens want enough
power for $+0.5$ to be detectable, and the read-coverage and replicate
choices should follow. Two replicates is the absolute minimum for any
quantitative claim; three is standard practice in 2024; four is
worthwhile for weak-phenotype screens.

The *reproducibility story* is, on net, encouraging. *Erik Hanson*'s
2022 cross-study reanalysis showed that CRISPR-screen replication
rates between independent groups working on the same cell line and
condition are roughly 80 percent — much higher than the 30 — 50 percent
seen for microbiome studies or for early GWAS hits. Three factors
drive this. *Standardised libraries* (everyone uses Brunello or TKOv3)
remove a major source of cross-study variance. *Cell-line uniformity*
— a cancer cell line is the same genetic background across labs, in
contrast with patient cohorts — removes biological-heterogeneity
variance. And *fewer biological covariates* across the lab-to-lab
transfer reduces the dimensions on which a result can vary.

But the same uniformity is itself a limitation. Lab-adapted cell lines
do not perfectly reflect primary patient tissue, and CRISPR-screen
dependencies in cell lines occasionally fail to translate. The
*CCLE-vs-PDX* discrepancies are well-documented for several
oncogenes (notably KRAS, where some KRAS-mutant cell lines show much
weaker KRAS dependency than the patient-derived xenograft analogues
do). The current best-practice mitigation is to validate dependency
hits in patient-derived organoids and PDX models before committing to
a drug-discovery programme — a substantial additional step that the
top of the pipeline rarely advertises.

#figure(
  image("../../diagrams/lecture-24/06-workflow.svg", width: 92%),
  caption: [
    A working CRISPR-screen analysis pipeline. Each band has its own
    canonical tool and its own reproducibility hazard; QC checks
    before analysis are non-negotiable, and tool diversity at the
    aggregation step (MAGeCK + DrugZ + BAGEL2) is the standard cross-check.
  ],
) <fig:ch24-workflow>


== From a Hit to a Drug <sec:translation>

A screen hit is not a drug. The translation pipeline from the
ranked-genes output of MAGeCK to an FDA-approved compound runs through
seven gates, and most candidates fail at one of them. Walking the
gates briefly is useful both for calibrating expectations and for
seeing where computational biology stops and pharmacology takes over.

*Gate one — hit replication.* The first thing anyone does with a
hit list is replicate it in an independent screen, ideally with a
non-overlapping library and a different cell line. Roughly half the
gene-level hits in a single screen fail to replicate; the half that
do are the working set.

*Gate two — orthogonal validation.* Re-introduce the gene knockout
in an arrayed format — one well, one guide, one phenotype — and
confirm that the loss-of-function phenotype matches the pooled-screen
signal. The same step typically includes *isogenic-pair experiments*:
two cell lines that differ only at the gene of interest, derived from
the same parental line by CRISPR knock-in or knockout. If the
phenotype tracks the isogenic genotype, the dependency is real.

*Gate three — druggability assessment.* Not every gene's protein
product can be drugged with a small molecule. The classical
druggable-genome estimate (Hopkins and Groom, 2002) is about 10 — 15
percent of protein-coding genes — those with concave binding pockets
that small molecules can occupy with high affinity. Transcription
factors, scaffolding proteins, and disordered regions are
traditionally hard targets, although molecular glues and PROTACs are
beginning to crack some of them. Surface antigens are addressable by
antibody-drug conjugates or by CAR-T cells. The hit triage at this
gate determines what therapeutic modality the eventual drug will use.

*Gate four — hit-to-lead chemistry.* Once a target is committed, the
medicinal-chemistry programme begins. Virtual screens, fragment-based
screens, HTS campaigns generate a starting compound; iterative
optimisation improves potency, selectivity, and pharmacokinetics. The
typical hit-to-lead timeline is twelve to thirty-six months.

*Gate five — pre-clinical pharmacology.* In-vivo efficacy in mouse
models (xenografts of the relevant cancer type), safety in two animal
species, ADME profiling. The fraction of candidates that fail at this
stage for toxicity, off-target activity, or insufficient exposure is
high — roughly half.

*Gates six and seven — clinical trials.* Phase I (safety, ~30
patients, one to two years), Phase II (efficacy, ~100 — 300 patients,
two to four years), Phase III (large randomised trials, hundreds to
thousands of patients, three to five years), then FDA review. The
median oncology-drug development timeline from target nomination to
FDA approval is eight to twelve years. The attrition through
clinical trials is severe — roughly one in ten Phase-I oncology
candidates reaches approval (Wong, Siah, Lo, 2019, _Biostatistics_).

#figure(
  image("../../diagrams/lecture-24/12-screen-to-drug.svg", width: 92%),
  caption: [
    From CRISPR-screen hit to FDA approval. Seven gates, eight to
    twelve years, roughly one-in-ten survival rate from Phase I.
    PARP inhibitors and the MTAP / PRMT5 programme are the textbook
    waypoints; the 2020s pipeline is dense enough that two or three
    DepMap-anchored drugs reach approval per year.
  ],
) <fig:screen-to-drug>

The 2020s pharma landscape has CRISPR screens embedded across the
discovery side. Roche, Genentech, Pfizer, and Novartis have internal
DepMap-style operations. AI-first companies (Recursion, Tempus,
Nimbus, In8) have screen-derived target portfolios at the centre of
their pipelines. The variant-interpretation half of the story has
reached the clinic too: Vertex Pharmaceuticals' *Casgevy*, the
2023-approved gene-therapy product for sickle-cell disease and
beta-thalassemia, descends directly from the CRISPR-Cas9 mechanism
walked at the start of this chapter (although it is in-vivo editing of
the BCL11A enhancer, not a screen).

The applications beyond oncology have lagged but are arriving.
*Anti-viral host-factor screens* — the SARS-CoV-2 entry-factor work
from Daniloski, Wei, and others in 2021 surfaced ACE2, furin, and
cathepsin L as required host factors and shaped early therapeutic
hypotheses. *Inflammation and autoimmunity* screens in primary immune
cells (Schmidt et al., 2022, in T-cells; Freimer et al., 2022, in
macrophages) have opened new target classes. *Neurodegeneration*
screens in iPSC-derived neurons are now standard for early-stage
target identification in Alzheimer's and Parkinson's, although the
in-vivo translation is harder than in oncology.

#tip[
  When you read a 2024-vintage pharma-pipeline announcement and it
  describes a target as "synthetic-lethal" or "lineage-selective", the
  evidence almost always traces back to a DepMap entry. The DepMap
  portal (depmap.org) lets you check the evidence directly — pick the
  gene, look at the gene-effect score across the cell-line panel, ask
  whether the claimed selectivity holds in independent data. Five
  minutes of due diligence on the portal is the difference between
  reading a press release and reading the data.
]


== Summary <sec:ch24-summary>

- A pooled CRISPR screen is an inverse problem on a sparse linear
  system. Roughly 80,000 guide-level log fold-changes feed a design
  matrix with four non-zero entries per column; standard solvers
  (MAGeCK, DrugZ, BAGEL2) invert it to recover roughly 20,000
  per-gene effects with appropriate noise models.

- CRISPR knockout, interference, and activation are three modalities
  of the same backbone. CRISPRko removes a gene by NHEJ-induced
  frameshift; CRISPRi knocks it down with a dCas9-KRAB repressor;
  CRISPRa drives it up with a dCas9-VP64 activator. The choice
  depends on whether you can tolerate cell death and which side of
  the gene's dosage response you want to interrogate.

- Library quality is the screen's first inductive bias. Brunello and
  TKOv3 are the modern defaults; Doench Rule Set 2 is the
  standard-of-care for guide picking; CFD score handles off-target
  burden. Legacy libraries have many low-activity guides and should
  be retired for new screens.

- MAGeCK's robust rank aggregation is a Beta(1, $k$) tail test on the
  best-of-$k$ guide-level rank. It is the same family of order-statistic
  tests that drives GSEA; the modification (median rather than
  minimum) absorbs one outlier guide gracefully.

- Quality control before analysis is non-negotiable. Four checks —
  guide recovery $> 90$ percent, replicate concordance $> 0.9$,
  CEG dropout $> 80$ percent, non-targeting null centred at zero —
  decide whether a screen is interpretable.

- DepMap industrialises the dropout screen across roughly 1,100
  cancer cell lines and produces, per gene, the cross-lineage
  essentiality pattern. Pan-essentials are toxic targets;
  pan-non-essentials are uninteresting; lineage-selective and
  context-selective dependencies (MCL1, EZH2, WRN, PRMT5 in
  MTAP-deleted cancers) are the drug-target gold.

- Synthetic lethality is the screening-era reformulation of the
  Hartwell-Kaelin idea. Olaparib in BRCA1/2-deficient cancers is the
  pre-CRISPR clinical archetype; MTAP / PRMT5 is the DepMap-anchored
  modern example.

- MAVEs measure every possible variant in a region at once. Saturation
  genome editing of BRCA1 (Findlay 2018) reclassified ~96 percent of
  prior variants of uncertain significance; the resulting per-variant
  scores feed ACMG/AMP PS3 / BS3 evidence and bridge directly into
  the clinical-variant-interpretation pipeline of Chapter 25.

- Pitfalls cluster in five families: library (use modern, validated
  ones), experimental (control MOI and cell coverage), bioinformatic
  (correct for copy-number bias, watch for plasmid contamination),
  statistical (use BH FDR and adequate replication), and biological
  (cell lines are not patients — validate in organoids and PDX before
  committing to a programme).

- The screen-to-drug pipeline runs through seven gates over 8 — 12
  years. Roughly one in ten Phase-I candidates reaches FDA approval;
  most attrition is at pre-clinical pharmacology. The 2020s pharma
  industry has CRISPR screens embedded across discovery, and two to
  three DepMap-anchored drugs are reaching approval per year.


== Exercises <sec:ch24-exercises>

#strong[1.] #emph[End-to-end MAGeCK.] Download a public CRISPR-screen
count table from DepMap (the public 23Q4 release contains the raw
counts behind every Chronos run). Pick one cell line. Run MAGeCK
with default parameters and produce a ranked list of essential genes.
Report the top ten essentials and identify which are core-essentials
versus cell-line-specific. Show the volcano plot of log fold change
versus $-log_(10)$ FDR.

#strong[2.] #emph[Tool comparison.] Run the same count table through
MAGeCK, DrugZ, and BAGEL2. Produce a 3-way Venn of significant hits
at FDR $< 0.05$. Where do the tools agree, and where do they
disagree? Inspect ten genes called by only one tool and try to
explain what about their guide-level pattern produces the discrepancy.

#strong[3.] #emph[Beta(1, k) by simulation.] Generate ten thousand
sets of four uniform-on-$[0, 1]$ samples. For each set, record the
minimum. Plot the empirical CDF of the minima and overlay the
analytic Beta(1, 4) CDF, $F(r) = 1 - (1 - r)^4$. Verify that they
agree. Repeat with $k = 6$ and confirm the analytic CDF is now
$1 - (1 - r)^6$.

#strong[4.] #emph[DepMap browsing.] On depmap.org, pick a cancer
lineage you find biologically interesting. Identify the top five
lineage-selective dependencies (highest gene-effect contrast between
that lineage and all others, at moderate-to-strong essentiality in
the chosen lineage). Cross-reference each with the PRISM drug-screen
data: is there a compound that mirrors the gene-knockout pattern? If
so, name it and check whether it has a clinical-development programme
in the relevant indication.

#strong[5.] #emph[MAVE alignment with ClinVar.] Download the BRCA1
saturation genome-editing data from MAVEdb (Findlay 2018). Cross-reference
the top 100 ClinVar Pathogenic BRCA1 variants with the MAVE function
scores. Compute the per-variant concordance (Pathogenic should be
non-functional in MAVE). Inspect the disagreements and propose
hypotheses for each: assay limitation? Domain not covered by the
study? Mis-annotated ClinVar entry?

#strong[6.] #emph[Library design.] Pick a 50-gene custom gene set
relevant to a phenotype of your choice. Use CRISPick or its API to
generate four optimal sgRNAs per gene under Doench Rule Set 2, with
off-target burden under CFD. Report the mean on-target activity and
mean predicted off-target score for your library. Compare against
the equivalent statistics for the matching subset of Brunello.

#strong[7.] #emph[Copy-number-bias diagnostic.] Generate a synthetic
count table with two cell lines: one with normal ploidy and one with
a 10 Mb high-copy-number amplification on chromosome 7. Run MAGeCK
without and with CRISPRcleanR copy-number correction. Show the
difference in the top-50 essential gene list. Compute the false
positive rate for the uncorrected run.

#strong[8.] #emph[(Open-ended.)] Pick a 2023 — 2024 oncology drug that
entered Phase I or Phase II clinical trial with a CRISPR-screen-anchored
rationale. Trace the discovery story: which screen first surfaced
the dependency, which DepMap entry confirms it, which independent
replication exists in the literature, and which patient stratification
the trial design uses. Write a one-page case-study summary linking
the upstream evidence to the downstream trial.


== Further Reading <sec:ch24-further-reading>

- *Jinek, M., Chylinski, K., Fonfara, I., Hauer, M., Doudna, J. A., &
  Charpentier, E.* (2012). "A programmable dual-RNA-guided DNA
  endonuclease in adaptive bacterial immunity." _Science_ 337: 816 — 821.
  The biochemistry paper that started the engineering revolution.

- *Shalem, O., Sanjana, N. E., Hartenian, E., et al.* (2014).
  "Genome-scale CRISPR-Cas9 knockout screening in human cells."
  _Science_ 343: 84 — 87. The first pooled-screen demonstration at
  genome scale.

- *Li, W., Xu, H., Xiao, T., et al.* (2014). "MAGeCK enables robust
  identification of essential genes from genome-scale CRISPR/Cas9
  knockout screens." _Genome Biology_ 15: 554. The canonical screen
  analyser.

- *Doench, J. G., Fusi, N., Sullender, M., et al.* (2016). "Optimized
  sgRNA design to maximize activity and minimize off-target effects of
  CRISPR-Cas9." _Nature Biotechnology_ 34: 184 — 191. The Rule Set 2
  paper and the basis for modern library design.

- *Tsherniak, A., Vazquez, F., Montgomery, P. G., et al.* (2017).
  "Defining a Cancer Dependency Map." _Cell_ 170: 564 — 576. The DepMap
  founding paper.

- *Dempster, J. M., Boehm, J. S., McFarland, J. M., et al.* (2021).
  "Chronos: a CRISPR cell fitness time-series model that improves
  identification of cancer dependencies." _Genome Biology_ 22: 343.
  The current DepMap scoring algorithm.

- *Findlay, G. M., Daza, R. M., Martin, B., et al.* (2018). "Accurate
  classification of BRCA1 variants with saturation genome editing."
  _Nature_ 562: 217 — 222. The textbook clinical-grade MAVE.

- *Iorio, F., Behan, F. M., Gonçalves, E., et al.* (2018). "Unsupervised
  correction of gene-independent cell responses to CRISPR-Cas9
  targeting." _BMC Genomics_ 19: 604. The CRISPRcleanR copy-number
  correction.

- *DepMap portal*. `depmap.org`. The data and the entry point for
  every dependency question raised in this chapter.

- *MAVEdb*. `mavedb.org`. The community repository of multiplexed
  assays of variant effects.
