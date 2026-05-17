#import "../theme/book-theme.typ": *

= ML in Genomics: Architectures, Pitfalls, Frontiers <ch:ml-genomics>

#matters[
  Most of the methods in the preceding fifteen chapters arrived in the
  same way. Somebody noticed that a genomic measurement had a particular
  geometric or statistical shape — a 2-D pileup, a 1-D regulatory
  sequence, an over-dispersed count vector, a molecular graph, a 3-D
  protein — and somebody else discovered that a particular neural
  network architecture had an inductive bias that matched that shape
  almost exactly. The match made the architecture much more
  sample-efficient than its general-purpose alternatives. Get the match
  right and a tractable training set is enough; get it wrong and no
  amount of compute will save you. This chapter pulls the pattern out of
  the individual case studies. The five canonical pairings — pileups to
  CNNs, long-range regulation to dilated convolutions plus transformers,
  count matrices to negative-binomial VAEs, molecules to graph networks,
  3-D structure to equivariant attention — share a single design move,
  and three pathologies — label scarcity, data leakage, and oversold
  foundation models — share three responses. The chapter is the
  blueprint, and it is also the warning label.
]

A genome is a long string drawn from a four-letter alphabet, but every
piece of biology you have learned to do with it has involved a different
geometric arrangement of bytes. Variant calling presented you with a 2-D
array of stacked reads. Differential expression handed you a vector of
counts per gene. Single-cell analysis stacked those vectors into a
matrix that was mostly zeros. Long-range regulation gave you a
hundred-kilobase stretch of one-hot-encoded DNA and a thousand
quantitative readouts. Protein structure presented an L-by-L pair
representation that wanted to be a 3-D arrangement. Each of those
shapes admitted a different machine-learning architecture, and the
match between shape and architecture is the most useful thing you can
learn about genomics ML.

The chapter is organised around that observation. The first three
sections walk the canonical architecture–data pairings and extract
their common move. The fourth deals with the label-scarcity problem
that every genomics task eventually hits, and with the
self-supervised pretraining that the field has converged on as the
answer. The fifth section is the one most working bioinformaticians
discover the hard way: genomic data is correlated in ways that break
naive cross-validation, and the resulting *data leakage* is the most
common reason a published number does not replicate. The remaining
sections survey the current state of DNA language models and cell
foundation models honestly — what they buy, what they do not buy,
and how to read their evaluation tables without being fooled — and
close on the multimodal frontier.

This is the synthesis chapter. The technique that ran through
Chapters 4, 8, 9, and 15 — match the architectural prior to the
structural property of the data — is here unfolded as an explicit
design discipline. The remaining chapters of the book assume you can
do it.


== Five Pairings, One Design Move <sec:pairings>

The five problems below recur across genomics. Each presents data in a
distinct geometric arrangement, and each is solved by a network whose
inductive bias matches that arrangement. The vocabulary differs from
chapter to chapter, but the move is the same.

The first pairing is *pileups and convolutional networks*. The variant
caller you met in Chapter 4 receives, at every candidate site, a small
2-D array — rows are reads, columns are positions around the
candidate, channels are base identity plus base quality plus strand
plus mapping quality plus a handful of auxiliary signals. A typical
window is on the order of a hundred rows by two hundred columns by six
channels. The classification target is one of four genotypes. The
right architecture is a 2-D convolutional neural network — exactly
the architecture Yann LeCun's group introduced for handwritten-digit
recognition (LeNet, 1998) and that AlexNet (Krizhevsky, Sutskever, and
Hinton, 2012) scaled to natural images. DeepVariant (Poplin et al.,
*Nature Biotechnology* 2018) was the moment this transfer became
explicit: a plain ImageNet-style CNN fed pileup images outperformed a
decade of hand-engineered statistical machinery in
GATK's HaplotypeCaller. The CNN's two priors — locality (a small
filter sees a small neighbourhood) and translation equivariance
(the same filter applies at every column) — line up with what variant
evidence looks like. Mismatches and deletions are local features.
Whether they sit at column eighty or column one-twenty is irrelevant.
Pool the local evidence across reads and you get a variant call.

#figure(
  image("../../diagrams/lecture-16/01-deepvariant-cnn.svg", width: 92%),
  caption: [
    DeepVariant's pileup-as-image trick. A 2-D CNN's translation
    equivariance and locality bias match the structural property of
    the pileup exactly — the column of the candidate does not matter,
    and the predictive signal is a local pattern of mismatches across
    a handful of reads.
  ],
) <fig:deepvariant>

DeepVariant was not invented in a vacuum. Earlier work — DeepBind
(Alipanahi et al., 2015) for transcription-factor binding from DNA
sequence, Basset (Kelley, Snoek, and Rinn, 2016) for chromatin
accessibility — had already shown that a 1-D CNN over one-hot-encoded
DNA could rediscover known TF motifs as filters in its first
convolutional layer and predict binding intensities better than
position-weight-matrix methods. The 2018 PrecisionFDA truth
challenges, where DeepVariant won the SNP and indel categories on
both Illumina and PacBio data, cemented the pattern: an ImageNet
backbone, retrained on a genomics task with a few thousand
high-confidence labels, beats hand-tuned heuristics. The remainder of
the genomics-CNN literature is variations on the same theme.

The second pairing is *long-range regulatory sequence and dilated
convolutions plus self-attention*. The Chapter 9 problem was to
predict a vector of regulatory readouts — chromatin accessibility in
several hundred cell types, CAGE expression, ChIP-seq tracks — from
a hundred kilobases of surrounding DNA. The input is a 100,000-by-4
one-hot tensor. Pure self-attention is quadratic in length and would
require $10^{10}$ pairwise computations at this scale, which is
prohibitive. Pure convolutions cannot reach a hundred-kilobase
receptive field without an absurd number of layers. The Enformer
architecture (Avsec et al., *Nature Methods* 2021), which descends
from Basenji (Kelley, 2018) and Basenji2 (Kelley, 2020), resolves the
tension in two stages. A stack of dilated convolutions with dilations
$1, 2, 4, dots, 512$ down-samples the sequence and grows the
receptive field exponentially to about two kilobases — enough to
capture motifs and small clusters of them. The down-sampled
representation then feeds a transformer whose self-attention routes
information across the remaining hundred-kilobase span. Multi-scale
hierarchy meets sparse long-range routing, and the receptive field
ends up at about 196 kb without paying full quadratic cost.

#figure(
  image("../../diagrams/lecture-16/02-enformer.svg", width: 95%),
  caption: [
    Enformer. Dilated convolutions handle local-to-mid scale; the
    transformer routes the down-sampled features across the long
    range. The two-stage split is what makes hundred-kilobase
    receptive fields tractable.
  ],
) <fig:enformer>

The third pairing is *count matrices and variational autoencoders with
the right likelihood*. Single-cell RNA-seq (Chapter 7) gives you
roughly twenty thousand gene counts per cell, with mean expression
well below variance and a heavy zero-fraction from technical dropout.
A vanilla autoencoder minimises mean-squared error, which is the
log-likelihood of a Gaussian noise model. Counts are not Gaussian.
They are negative binomial, often zero-inflated. *scVI* (Lopez,
Regier, Cole, Jordan, and Yosef, *Nature Methods* 2018) replaces the
MSE reconstruction with a negative-binomial log-likelihood and
conditions the decoder on a batch covariate. The encoder maps a
cell's count vector to a ten-to-thirty-dimensional latent; the
decoder reconstructs the counts under an NB likelihood with
batch-aware dispersion. The architecture is straightforward; the
likelihood is the entire game. *totalVI* (Gayoso et al., 2021) adds
surface-protein measurements; *MultiVI* (Ashuach et al., 2023) adds
ATAC. None of them are deep — they are shallow MLPs by current
standards. The depth was never the point.

The fourth pairing is *molecular graphs and graph neural networks*.
Cheminformatics tasks — predicting solubility, binding affinity,
toxicity, ADMET properties from molecular structure — present a
chemical compound as a graph. Nodes are atoms with features (element,
hybridisation, formal charge); edges are bonds with features (single
or double, aromatic, in-ring). The right network is a message-passing
graph neural network: at each layer every node aggregates a learned
function of its neighbours' features, and after several layers each
node's representation encodes a neighbourhood of bounded radius.
Variants include the graph convolutional network of Kipf and Welling
(2017), the message-passing neural network of Gilmer et al. (2017),
and the graph isomorphism network of Xu et al. (2019). Two priors
match the chemistry. Permutation equivariance — the network's output
is unchanged by renumbering the atoms — is what allows the same
network to handle a benzene ring whose atoms are labelled in any
order. Locality — chemical properties are dominated by functional
groups and their local environment — is what the message-passing
range encodes. Beyond cheminformatics, graph networks appear inside
AlphaFold-3's diffusion module for ligand atoms and in many models
that predict protein–ligand interaction.

The fifth pairing is *protein structure and equivariant attention*.
Chapter 15 walked the Evoformer and the Invariant Point Attention
module of AlphaFold-2 (Jumper et al., *Nature* 2021). The
key priors there are axial attention on the L-by-L pair
representation (which makes full 2-D attention's $O(L^4)$ cost
collapse to $O(L^3)$), triangle multiplicative updates and triangle
attention (which enforce a geometric consistency: if residues $i$ and
$j$ are close in 3-D and $j$ and $k$ are close, then $i$ and $k$
inherit a triangle-inequality constraint), and SE(3)-equivariance in
the structure module (so that rotating the input rotates the output by
the same group element, with no data augmentation needed). AlphaFold
is the canonical demonstration that the right inductive bias plus a
lot of data produces a step change in accuracy. The CASP14 jump from
roughly 40 to roughly 90 GDT_TS is one of the most consequential
benchmark moves in modern computational science.

#figure(
  image("../../diagrams/lecture-16/03-architecture-map.svg", width: 95%),
  caption: [
    Five problems, five architectures, one rule. The architecture's
    inductive bias is chosen to match the data's structural property
    — translation equivariance for pileups, multi-scale plus sparse
    long-range for regulation, NB-likelihood for counts, permutation
    equivariance for graphs, SE(3)-equivariance for 3-D structure.
  ],
) <fig:archmap>

Across the five pairings the same move recurs. Each architecture
commits to a structural prior — a symmetry, an invariance, a
likelihood, a routing assumption — and the choice of prior is
calibrated to the structural property of the data. The stronger the
match, the more sample-efficient the network. The weaker the match,
the more data is needed before the network can learn the prior itself.

#note[
  "Commit to a prior" is the engineering reading of *inductive bias*.
  Stronger priors require less data on matching problems and fail
  badly on mismatched ones. Weaker priors are universal but
  data-hungry. The interesting design question is rarely "which
  architecture is best in the abstract" — it is "which prior fits
  the geometric and statistical shape of my measurements." The five
  pairings above are the canonical answers; the rest of the chapter
  is about how to evaluate whether you got the answer right.
]


== Inductive Bias, in Detail <sec:inductive-bias>

The casual definition of an inductive bias is "any assumption an
algorithm makes that the training data does not force on it." The
useful definition is sharper. Every supervised learner is a map from
training data to a hypothesis in some hypothesis class; the class
itself is the bias. Linear regression's hypothesis class is the set of
affine functions of the features. A random forest's class is the set
of axis-aligned piecewise-constant functions. A CNN's class is the
set of functions that are translation-equivariant and locally
parameterised. An equivariant network's class is the set of functions
that respect a specified group of symmetries.

The no-free-lunch theorem (Wolpert, 1996) says that, averaged across
all possible learning problems, no learner beats any other. The
theorem is vacuous as a statement about realistic problems and
sharply useful as a reminder: a learner that performs well on the
problems you actually face does so because its hypothesis class
matches those problems. Genomics is one of the most extreme examples
of this principle in modern ML, because the relevant priors —
translation equivariance for sliding-window features, NB statistics
for counts, SE(3) equivariance for 3-D geometry — are so specific to
the data that picking a network with the wrong prior produces an
order-of-magnitude penalty.

Two comparisons make the trade-off concrete. The first is on pileup-style
data. Train a fully-connected MLP, a CNN, and a transformer with
absolute positional embeddings on the DeepVariant task at matched
budget. The CNN reaches roughly 99.5 % F1; the transformer with
sinusoidal positions reaches about 98 % at the same budget and
roughly matches the CNN only when given ten times more training data;
the MLP underperforms both by several points and never closes the gap.
Add a 2-D positional bias plus local-attention windowing to the
transformer — that is, hand-engineer the CNN's priors back in — and
the transformer reaches CNN accuracy with comparable data. The
point is not that CNNs are better than transformers; it is that the
prior matters more than the family, and a "more general" architecture
needs to learn the prior from data when it could have been written down
for free.

#figure(
  image("../../diagrams/lecture-16/10-architecture-comparison.svg", width: 92%),
  caption: [
    Architecture comparison on a synthetic pileup task. The CNN's
    explicit translation-equivariance bias gives it a roughly 10 — 100 ×
    data-efficiency advantage on small data; the gap closes at very
    large data, where the transformer learns the equivalent prior
    from examples. The MLP, which has no useful prior, never closes
    the gap.
  ],
) <fig:archcompare>

The second comparison is the one Lopez and colleagues drew when
introducing scVI. Train a deep VAE on a sc-RNA-seq dataset under
three reconstruction likelihoods — Gaussian (MSE), Poisson, and
negative binomial. The Gaussian VAE predicts negative counts, drags
its estimates toward a mean that cannot accommodate the empirical
variance, and produces nonsense under any log-transform. The Poisson
VAE handles the integer constraint but assumes mean equals variance
— and single-cell data has variance much greater than mean, so the
fit under-estimates uncertainty everywhere. The NB VAE has a free
dispersion parameter that captures the over-dispersion, and the
ZINB extension handles 10× Chromium-era dropout zeros. Three
architectures, identical depth and width; the right likelihood buys
orders of magnitude in statistical efficiency. The architecture
depth was scaffolding. The likelihood was the algorithm.

#tip[
  When a deep-learning paper in genomics reports a large improvement
  over a previous method, the first question to ask is whether the
  improvement comes from the architecture or from a changed
  likelihood, a changed positional encoding, or a changed
  pretraining corpus. The famous improvements — DeepVariant's CNN,
  scVI's NB, AlphaFold-2's triangle updates — all came from a
  specific structural change with a clear interpretation. Generic
  "we used a transformer" is rarely the answer.
]

The third move in the inductive-bias family is *equivariance*. If
your problem has a known symmetry — rotation in 3-D, permutation of
items in a set, time-shift in a signal — you can write the symmetry
into the network and it will hold by construction. Translation
equivariance is what makes a CNN a CNN; group convolutions (Cohen and
Welling, 2016) generalise the idea to 2-D rotations; SE(3)-equivariant
networks (NequIP, e3nn, SE(3)-Transformer) generalise it to 3-D
rotations and translations together; permutation-equivariant networks
generalise it to unordered sets and graphs. Each equivariance
constrains the function class to those that respect the group action,
and the constraint typically buys an order of magnitude or two in
sample efficiency on problems where the symmetry genuinely holds.
AlphaFold-2's Invariant Point Attention is the canonical example: it
treats coordinates not as raw $(x, y, z)$ but as frames attached to
each residue, and every attention update is computed in a way that is
exactly invariant to rotations and translations of the global frame.
No data augmentation, no rotation jitter — the symmetry is built in.

#figure(
  image("../../diagrams/lecture-16/04-inductive-biases.svg", width: 95%),
  caption: [
    Architecture / inductive-bias / data-regime cheat sheet. Pick by
    the structural property of the data, not by familiarity with the
    architecture.
  ],
) <fig:inductive>

A useful framing for an EE-trained reader is the matched-filter
analogy. A matched filter in signal processing is optimal when both
the signal template and the noise covariance are known exactly. Get
either wrong — wrong template, wrong noise model — and performance
collapses. An ML architecture's inductive bias is the generalised
version: the symmetries it can efficiently represent play the role of
the template; the reconstruction or output likelihood plays the role
of the noise model. The matched filter for a rectangular pulse in
white Gaussian noise is the CNN with MSE on a denoising task; the
matched filter for over-dispersed counts is the NB-likelihood VAE;
the matched filter for rotation-invariant 3-D structure is an
SE(3)-equivariant network. Architecture design is matched-filter
design at scale, and the same intuition — get the template and the
noise model right — carries over directly.

#note[
  When you do have enough data, the inductive-bias advantage
  collapses. Dosovitskiy et al.'s Vision Transformer (2020) showed
  that on a hundred-million-image pretraining corpus, a plain
  transformer learns the equivalent of local translation
  equivariance from data and beats the strongest CNNs. AlphaFold-2's
  Evoformer learns enough about MSA statistics from
  150 million-sequence pretraining that it does not need the explicit
  contact priors that EVfold relied on. The general rule is: the
  more data, the weaker the prior you can afford to commit to. Most
  genomics tasks are still in the data-limited regime where the
  prior matters; the exceptions — protein structure, language
  modelling on bulk sequence — are the ones where transformers
  eat their CNN-flavoured ancestors.
]


== Information Theory in Genomics <sec:info-theory>

The inductive-bias question reduces, in the end, to one Shannon would
recognise: what is the minimum description length of the signal you
are trying to recover? The right architecture is the one whose
hypothesis class is small enough to generalise from limited data, but
rich enough to express the true signal. Several classical
information-theoretic ideas reappear in this chapter as the
quantitative language for that trade-off.

The *Shannon entropy* of a random variable $X$ over an alphabet
$Sigma$ is $H(X) = -sum_(x in Sigma) p(x) log_2 p(x)$. For a
uniform-random DNA sequence, $H = 2$ bits per base. Real genomes sit
below that ceiling because the empirical base frequencies are
non-uniform, regional GC content varies (gene-rich isochores are
typically 50 — 60 %, gene-poor regions 35 — 45 %), and longer-range
correlations — codon usage, CpG islands, repeat-element families —
shave off further bits. Generic compression of a human chromosome
with `gzip` reaches roughly 1.4 bits per base; reference-aware
genomic compressors (Genozip, Spring) reach about 0.4 bits per base
by exploiting redundancy with the reference and across reads. The
gap between 2 bits and 0.4 bits is where biology lives.

Per-position information content drives the position-specific
scoring matrices that have been the backbone of motif discovery for
thirty years and the sequence logos that visualise transcription-factor
preferences. For column $i$ of an aligned set of binding sites, the
information content is

$ "IC"_i = log_2 |Sigma| - H_i = 2 - (- sum_b p_(i, b) log_2 p_(i, b)) $

A fully conserved column has $"IC" = 2$ bits; a uniform column has
$"IC" = 0$. The total information of a motif, $sum_i "IC"_i$, is a
direct read on how easily it can be detected against a random-DNA
background — it is the matched-filter signal-to-noise ratio in
information units.

*Mutual information* between two positions, $I(X; Y) = sum p(x, y)
log frac(p(x, y), p(x) p(y))$, measures statistical dependence
without assuming a functional form. Three genomics applications use it
directly. Coevolution detection (Chapter 15) uses pairwise MI across
MSA columns as a proxy for spatial contact in 3-D. Regulatory-network
inference (ARACNe, CLR, GENIE3) uses MI between TF expression and
target expression to score putative regulatory edges. Feature
selection in GWAS uses MI between a SNP and a phenotype as a
model-free alternative to the linear-regression $p$-value, robust to
non-linear association forms. The DCA story in Chapter 15
is the most refined version of the MI move — adding an
inverse-covariance step to separate direct from indirect coupling —
but the underlying intuition is unchanged.

The genetic code itself is, on closer inspection, an
*error-correcting code*. The 64-codon to 20-amino-acid mapping is a
$(64, 20)$ code with structured redundancy. Mutations at codon
position 3 are usually synonymous because of the wobble degeneracy.
Mutations at position 1 typically change one amino acid for a
chemically similar one (hydrophobic for hydrophobic, charged for
charged). Mutations at position 2 are the most disruptive and tend to
change the chemical class outright. Freeland and Hurst (*Journal of
Molecular Evolution*, 1998) showed that the standard genetic code sits
in the top 0.02 % of randomly permuted code tables for minimising the
expected damage from random mutations. The code is error-resistant
*by selection*, without explicit parity bits. Few engineering systems
have it so good.

#figure(
  image("../../diagrams/lecture-16/13-info-theory.svg", width: 95%),
  caption: [
    Information theory in genomics. Sequence logos use per-column
    information content for motif detection; pairwise mutual
    information identifies coevolving sites; the genetic code's
    third-position wobble is a built-in error-correcting code that
    Freeland and Hurst showed is near-optimal among random
    permutations.
  ],
) <fig:info>

Several genomics measurements are *compressed-sensing* problems in
disguise. Single-cell RNA-seq observes a tiny fraction of a cell's
mRNA molecules; the recovery of a usable expression matrix from those
sparse observations is a low-rank matrix completion problem. DIA
proteomics co-fragments overlapping peptides and recovers individual
identities by sparse decomposition. Pooled CRISPR screens read a
population-level phenotype that is a linear measurement of the
underlying perturbation matrix. Hi-C contact recovery downsamples a
three-dimensional contact tensor and reconstructs it under sparsity
constraints. The pattern is the same in every case: an
underdetermined linear system $bold(y) = A bold(x)$ where $bold(x)$
is sparse, and recovery via $ell_1$-regularised convex optimisation
or low-rank factorisation. Donoho and Candès–Romberg–Tao formalised
this in 2006; modern genomics has rediscovered it independently in
each domain.

#note[
  The recurring pattern across these examples is the same one that
  motivates the inductive-bias choice in @sec:inductive-bias. In
  channel coding, the trade-off is between code rate and error
  correction — more parity bits give better recovery but lower
  information density. In inductive bias, the trade-off is between
  hypothesis-class size and generalisation — a stronger prior
  generalises better from less data but cannot represent signals
  outside its class. The matched filter, the constant-time decoder,
  the LDPC sparse mixing matrix — each is a point in a rate–distortion
  plane for recovering genomic signals. The right architecture is
  determined by the signal's source distribution, just as the right
  code is determined by the channel's noise distribution.
]


== Training Data Is the Bottleneck <sec:data>

Architecture choices are not free; they are constrained by how much
labelled data you have. In genomics the labelled-data budget is
almost always much smaller than it is in comparable ML domains, and
the shape of the data pyramid is inverted.

#figure(
  image("../../diagrams/lecture-16/05-data-pyramid.svg", width: 95%),
  caption: [
    The genomics labelled-data pyramid is inverted. Unlabelled
    sequence is abundant; gold-standard functional labels — DMS,
    MPRA, ClinVar pathogenicity, curated drug responses — are
    thousands of orders smaller. Foundation-model strategies exist
    largely to bridge that gap.
  ],
) <fig:pyramid>

Ranked by abundance, the genomics data hierarchy looks roughly like
this. At the top sit *unlabelled sequences*: tens of terabytes of DNA
and protein in UniProt, NCBI, the SRA, and Ensembl, with no labels
attached but enormous statistical regularity to learn from. Just
below are *bulk RNA-seq* samples — millions in GEO and the SRA, each
labelled imperfectly with a tissue or condition string and with
substantial cross-study batch effects. *GWAS associations* contribute
hundreds of thousands of SNP–trait pairs at biobank scale, noisy but
extensive. *Single-cell expression* now covers tens to a few hundred
million cells across the cellxgene and Human Cell Atlas archives,
with cluster-level labels that change every time the data is
reclustered. *Protein structures* number about 180,000 experimentally
determined entries in the PDB and about 214 million predicted entries
in the AlphaFold DB — and the predicted ones are usable as labels if
filtered by pLDDT. *Gene annotations* (RefSeq, Ensembl, GENCODE) are
in the tens of thousands per genome, mostly curated. *Variant
pathogenicity* in ClinVar contributes a few hundred thousand
expert-curated calls, growing slowly. *Functional experimental data*
from deep mutational scanning and massively parallel reporter assays
contributes hundreds of thousands of variant effects across about a
thousand genes — gold-standard signal, but small. *Clinical
outcomes* are mostly access-controlled (dbGaP, EGA), comparably tiny,
and rarely cross-study compatible.

The shape matters because supervised deep learning wants millions of
labelled examples and most genomics tasks meet the quality bar with
hundreds to thousands. The consequences are structural rather than
tactical. *Foundation-model strategies dominate*: pretrain on the
abundant unlabelled bases, fine-tune on the scarce labelled outcomes.
Self-supervised objectives — masked language modelling on DNA or
protein sequences, contrastive learning on paired biological views —
are the workhorses. Evaluation is perilous because the test sets are
small and noisy, and differences within a few percent are rarely
significant. Transfer between tasks is the only way most tasks get
solved at all.

*Self-supervised pretraining* is the move that made label-scarce
genomics tractable. You take unlabelled data, invent a label by
hiding part of the input, and train the network to predict the hidden
part. After enough such pretraining iterations, the network has
internalised the statistical structure of the data; you then fine-tune
on a small labelled task and the network does not start from scratch.
Two pretraining objectives dominate. *Masked language modelling*
hides about 15 % of tokens and asks the network to fill them in
from context — the BERT recipe (Devlin et al., 2018) for natural
language, adapted to DNA k-mers in DNABERT (Ji et al., 2021) and to
amino acids in ESM-1 (Rives et al., 2021) and its successors.
*Causal next-token prediction* asks the network to predict the next
token given the prefix — the GPT recipe, adapted to DNA in HyenaDNA
(Nguyen et al., 2023) and Evo (Nguyen et al., 2024), and to proteins
in ProGen (Madani et al., 2023). A few other objectives appear at
the edges: contrastive learning (SimCLR, MoCo) on paired biological
views, BYOL-style self-distillation, and protein-specific tricks like
inverse-folding pretraining (ESM-IF) where the input is a backbone
and the target is the sequence.

#figure(
  image("../../diagrams/lecture-16/06-ssl-flow.svg", width: 92%),
  caption: [
    Self-supervised pretraining. The ratio is what matters:
    hundreds of millions of unlabelled examples in pretraining, tens
    of thousands of labelled examples in fine-tuning. Pretraining
    bridges the gap that pure supervision cannot cross.
  ],
) <fig:ssl>

#figure(
  image("../figures/ch16/f3-ssl-data-efficiency.svg", width: 95%),
  caption: [
    Data efficiency of self-supervised pretraining on a downstream
    genomics task. Pretrained-and-fine-tuned models reach the same
    validation accuracy as scratch-trained models with roughly 10 — 100 ×
    fewer labels. The biggest gains arrive in the low-label regime
    that most genomics applications occupy.
  ],
) <fig:ssl-eff>

The data-efficiency gain is large enough to change the kinds of
projects that are feasible. A task with a thousand labelled examples
— a deep-mutational-scan-style variant-effect dataset for a single
gene, a small panel of regulatory annotations — was hopeless as a
supervised training set in 2017. After pretrained models became
available, the same thousand labels reach the accuracy that ten
thousand or a hundred thousand would have required from scratch.
@fig:ssl-eff sketches the typical shape.

Weak supervision matters at the margin. Distant supervision —
inferring labels from auxiliary data, as when GO terms are propagated
to uncharacterised proteins by sequence similarity — provides large
but noisy labelled sets that work when the noise is uncorrelated with
the signal of interest. Heuristic rule-based labelling, multi-task
learning across related tasks with more labels, and active learning
that prioritises expert labels for the samples the model is most
uncertain about are all standard ML moves that earn their keep in
genomics because the label budget is so tight.

#tip[
  In a fine-tuning setting, the first knob to turn is rarely the
  architecture. It is the pretraining corpus and the pretraining
  objective. A protein-LM pretrained on UniRef90 with masked-amino-acid
  prediction performs differently from one pretrained on the same
  corpus with inverse-folding, even at identical downstream
  fine-tuning. Audit the pretraining recipe before you blame the
  fine-tuning architecture.
]


== The Data Leakage Problem <sec:leakage>

Genomics data is correlated in ways that break naive cross-validation,
and the resulting *data leakage* — where train and test sets share
information they should not — is the most under-reported cause of
inflated benchmarks in the field. Shuffle a dataset and split it
80/20 in vanilla ML and you usually get a valid estimate of
generalisation. Do the same on genomic data and you measure
memorisation of correlated structure, not generalisation. The gap
between as-reported accuracy and after-leakage-correction accuracy is
routinely twenty to thirty percentage points. It is the single most
common reason a published genomics ML result fails to replicate.

The sources of leakage are structural to the biology. *Sequence
homology* puts protein sequences from related species at 60 — 95 %
identity; a random split puts close homologues on both sides of the
train/test boundary, and the test answers become trivially predictable
by nearest-neighbour lookup. *Genomic locality* puts nearby SNPs in
linkage disequilibrium; per-SNP random splits put correlated SNPs on
both sides. *Family relatedness* in biobank cohorts is enough that a
500 k-individual cohort like UK Biobank contains thousands of
unflagged cousin pairs; their shared haplotype segments allow a
model trained on one to score the other without learning the biology.
*Cross-study batch effects* in single-cell datasets — different
chemistry versions, different operator practices, different
ambient-RNA distributions — produce systematic technical signatures
that a classifier can latch onto in place of biological signal.

#figure(
  image("../../diagrams/lecture-16/07-leakage-splits.svg", width: 95%),
  caption: [
    Random versus chromosome-aware split. A random hold-out at the
    SNP level evaluates memorisation of linkage disequilibrium; a
    chromosome-block hold-out evaluates generalisation to unseen
    regions. The two protocols can disagree by twenty percentage
    points on the same model.
  ],
) <fig:splits>

#warn[
  Every train/test evaluation assumes the samples are independent
  and identically distributed and that the train and test partitions
  are independent draws. In genomics, neither assumption holds.
  Adjacent SNPs are autocorrelated along a chromosome; homologous
  sequences across species; samples from the same lab via batch
  effects. The right fix is *blocked cross-validation* — hold out
  whole regions of the correlation structure (chromosomes, families,
  studies, time windows) — so the test points are statistically
  independent of the training points. Get this wrong and your
  reported accuracy measures memorisation, not generalisation.
]

The remedies are domain-specific. For *protein-function or
protein-structure tasks*, the community standard is sequence-identity
clustering with CD-HIT or MMseqs2 at a 30 — 50 % identity threshold,
splitting at the cluster level so that no test cluster shares more
than the threshold identity with any training cluster. The CAFA
challenge (*Critical Assessment of Functional Annotation*, Radivojac
et al., 2013, running since 2010) takes a different cut at the same
problem with strict time-based splits: train on sequences whose
function was known by some date $T$; test on annotations added
between $T$ and $T + 2$ years. Both methods avoid the homology
leak; both occasionally still bleed via family-level correlations
that survive the clustering.

For *sequence-regulation tasks at hundred-kilobase scale*, the
standard is the *chromosome split*. Enformer (Avsec, 2021) holds out
chromosomes 8 and 9; Basenji and Basenji2 use similar partitions. The
Enformer authors did this carefully, and the model still leaked: a
2023 reanalysis showed that within-chromosome regulatory
auto-correlation meant the chromosome split was insufficient on its
tightest tasks. A stricter alternative is the *whole-species split*
— train on human, test on mouse — which is rarely done because most
labelled datasets are human-only but is the gold standard when
matched-species data exists.

For *biobank GWAS-style tasks*, the standard is *kinship-aware
splitting*. Compute the genome-wide kinship matrix between every pair
of individuals; reject any train/test split that contains pairs with
kinship above about 0.05 (third-cousin level or closer). GWAS
pipelines have implemented this for two decades. ML pipelines that
ingest biobank phenotypes often have not, with predictable results.

For *cross-study single-cell tasks*, the standard is
*leave-one-study-out* evaluation. Each held-out study evaluates
cross-study generalisation directly; aggregating across leave-one-out
runs gives a realistic estimate of how the model behaves on the next
unseen cohort. The Luecken et al. (2022) integration benchmark
enforces this protocol for atlas-level evaluations and produces very
different rankings from the within-study evaluations that preceded it.

#figure(
  image("../../diagrams/lecture-16/11-leakage-cases.svg", width: 95%),
  caption: [
    Three case studies of genomics ML benchmarks reported under naive
    splits versus tightened splits. Protein-function predictors from
    the mid-2010s, cancer-driver predictors from 2014 — 2018, and
    Enformer's 2021 chromosome split all lose substantial accuracy
    when the leakage is plugged.
  ],
) <fig:cases>

Three case studies put the magnitude on the table. Pre-2016
protein-function-prediction benchmarks routinely reported above 95 %
accuracy on naive splits and dropped to 40 — 60 % under CAFA's
time-based protocol. Cancer-driver-gene predictors from 2014 to 2018,
several of which trained and tested on overlapping known-driver
databases, lost twenty to thirty percentage points on independent
benchmarks (Tokheim et al., 2016). Even Enformer's initial-release
test passed its chromosome-split sanity check but lost meaningful
accuracy when a 2023 reanalysis tightened the protocol further.

#note[
  The right mental image for genomic data leakage is "adjacent
  pixels in an image." Hand a model a photograph and train it on a
  random 80 % of pixels, evaluating on the held-out 20 %, and the
  model barely needs to look at anything — it averages the
  neighbours. The exam was already seen. Genomic data is structurally
  identical: SNPs along a chromosome are the adjacent pixels;
  related individuals are adjacent patients; samples from one lab run
  are adjacent batches. Random hold-outs on any of these produce an
  exam you have already taken. Leakage-aware splits hold out whole
  regions of the correlation structure.
]

When a paper reports a striking accuracy, the first question to ask
is not whether the architecture is novel. It is how the train and
test sets were split. If the answer is "random 80/20", discount the
number heavily. If the answer is "stratified by gene", probe for
homology-cluster splits. If the answer is "chromosome split", ask
whether the model's effective receptive field bleeds across the split.
Reproducing a paper's result under a tighter split is a standard
sanity check; it frequently changes the headline number.


== DNA Language Models <sec:dna-lms>

A *DNA language model* is a large transformer or transformer variant
pretrained on genomic sequence with a self-supervised objective. The
recipe is direct: chop the genome into windows of 500 bp to 1 Mb,
tokenise (per-nucleotide, $k$-mer, or BPE-style), mask some fraction
of the tokens and predict them from context (or use causal next-token
prediction), and ship the resulting model for fine-tuning or
feature-extraction on downstream tasks. The hope, sometimes stated
explicitly, is that DNA pretraining will do for genomics what BERT and
GPT did for natural language.

The lineage runs from 2021 forward. *DNABERT* (Ji et al.,
*Bioinformatics* 2021) was the first scaled DNA-LM: BERT-style, $k$-mer
tokenisation, a 512-token context, and roughly a hundred million
parameters. *Nucleotide Transformer* (Dalla-Torre et al., 2024, from
InstaDeep and Cambridge) scaled to between 500 million and
2.5 billion parameters, used multi-species pretraining across 850
genomes, and extended the context to about twelve kilobases.
*HyenaDNA* (Nguyen et al., NeurIPS 2023) replaced quadratic attention
with the Hyena operator (Poli et al., 2023) — implicit long
convolutions parameterised by neural networks — and pushed the
context to one megabase at sub-quadratic cost. *Evo* (Nguyen et al.,
*Science* 2024, from Arc Institute) trained a Hyena-style model on
millions of microbial genomes and demonstrated convincing generative
quality at the megabase scale, including the design of plausible CRISPR
systems. *Caduceus* (Schiff et al., 2024) added reverse-complement
equivariance directly into the architecture, so the model treats a
strand and its complement as equivalent inputs by construction.

#figure(
  image("../../diagrams/lecture-16/08-dna-lm-timeline.svg", width: 95%),
  caption: [
    The DNA language-model landscape, 2021 — 2024. Context length is
    the dominant axis of competition because it determines what scale
    of regulatory interaction is learnable; parameter count is the
    secondary axis.
  ],
) <fig:dna-lm>

The tokenisation choice is itself an inductive-bias decision.
Per-nucleotide tokenisation gives the finest resolution but the
shortest effective context (a 512-token window covers 512 bases).
$k$-mer tokenisation packs more sequence into fewer tokens but loses
nucleotide-level positional precision and produces a vocabulary that
grows as $4^k$. Byte-pair encoding learns a vocabulary from
frequency statistics and lands somewhere in between. Reverse
complementarity adds a wrinkle that text does not have: the same
biological signal can sit on either strand, so the tokeniser and the
network together need to handle the symmetry, either by explicit
data augmentation or by built-in equivariance (Caduceus).

#figure(
  image("../figures/ch16/f1-dna-tokenisation.svg", width: 95%),
  caption: [
    Three DNA tokenisations of the same 24-base sequence —
    per-nucleotide, 6-mer (overlapping), and BPE-style — with
    resulting token counts and effective receptive field at a fixed
    512-token window. The tokeniser is a design choice; it sets the
    scale at which the model can attend.
  ],
) <fig:dna-tok>

An honest account of what DNA pretraining buys is shorter than the
marketing suggests. DNA LMs *do* recover motif-scale signal at no
extra training cost — TF binding motifs can be probed out of attention
heads or sequence-attribution gradients, in the spirit of DeepBind
(2015). They *do* deliver moderate, real, but not dramatic
improvements on downstream regulatory-prediction tasks — typically
one to five percentage points absolute over specialised baselines.
They *do* generalise across species better than single-species
specialised models; cross-species transfer is the clearest place where
the pretrained representation pays for itself. And they enable
*in-silico* mutagenesis at scale: scoring the difference between
reference and mutant context is essentially free once the model is
trained, and faster than SHAP-style attribution on specialised
networks.

What DNA pretraining does *not* yet deliver is a step change. Enformer
still beats DNA LMs on Enformer's own tasks. Long-range dependence is
nominal — HyenaDNA and Evo have megabase contexts in the name, but
distal-element interactions are still weakly learned by current
training recipes. Generative quality looks plausible at the
nucleotide level but is hard to validate experimentally — designed
sequences mostly have not been wet-lab tested at the depth that
would convince a critical reader. Species bias is real (most corpora
are dominated by human and a handful of model organisms). And
evaluation-set contamination is a recurring problem because
pretraining corpora often overlap published benchmarks.

#warn[
  DNA is not text in a four-letter alphabet, despite the convenient
  packaging of the analogy. The signal-to-noise per base is much
  lower than per English word; functional elements are sparse on
  the kilobase scale and absent on the megabase scale; tokenisation
  choices materially affect what scale is learnable; reverse-complement
  symmetry is real and has no analogue in text. The DNA-LM literature
  has occasionally overpromised on the strength of the analogy. The
  models are useful, but they are not the BERT moment for genomics
  — at least not yet.
]


== Foundation Models for Cells <sec:cell-fms>

A second wave of foundation models targets *cells* rather than
sequence. *Geneformer* (Theodoris et al., *Nature* 2023) was the
opening move: a BERT-style transformer pretrained on roughly thirty
million single-cell transcriptomes from the Human Cell Atlas, with a
rank-value tokenisation that orders genes by per-cell expression and
uses the rank as the token. *scGPT* (Cui et al., *Nature Methods*
2024) followed with binned expression values plus gene-identity tokens
as parallel channels, pretrained on tens of millions of cells.
*scFoundation* (Hao et al., *Nature Methods* 2024) used per-gene
real-valued expression via a custom embedding. *UCE* (Rosen et al.,
*bioRxiv* 2024) introduced a cross-species shared embedding so that
mouse and human cells live in the same latent space. The pitch is
uniform: a single model whose pretrained representations transfer to
any downstream single-cell task — cell-type classification,
perturbation-response prediction, gene-regulatory-network inference,
cross-species translation.

#figure(
  image("../../diagrams/lecture-16/09-cell-foundation-model.svg", width: 92%),
  caption: [
    A cell foundation model in the abstract: a cell's gene-expression
    vector is tokenised (rank-value, binned, or real-valued); a stack
    of transformer layers produces a cell-level embedding; downstream
    tasks consume the embedding. The architecture is conventional;
    the open question is whether tokenisation choices and pretraining
    corpora make the representation broadly useful.
  ],
) <fig:cell-fm>

The honest evaluation is mixed. *Cell-type classification under
leave-one-cell-out* is broadly easy — most cell types in the test set
are present in the training set under different labels, and the
foundation models reach 70 — 85 % zero-shot accuracy. *Cell-type
classification under leave-one-study-out* is the honest version of the
same protocol, and the foundation-model numbers drop substantially.
The gap is a frequent source of misleading marketing in the
literature. Perturbation-response prediction shows modest but real
success on scGPT. Cross-species cell-type alignment is the clearest
win for UCE, where the shared embedding genuinely improves
mouse-to-human transfer. Gene-gene interaction discovery via attention
is suggestive but hard to validate against ground truth.

#figure(
  image("../figures/ch16/f2-eval-protocol-gap.svg", width: 95%),
  caption: [
    The same model under two evaluation protocols. Leave-one-cell-out
    (LOCO) trains on a random fraction of cells and tests on the
    rest, including cells from the same studies. Leave-one-study-out
    (LOSO) holds out entire studies, removing the cross-study batch
    correlation. The accuracy difference between the two is what
    "zero-shot" reporting often obscures.
  ],
) <fig:eval-gap>

What is oversold. Many published numbers report zero-shot accuracy on
cell types that are present in pretraining under different ontology
labels — re-cluster the data, rename the clusters, and the
"zero-shot" task becomes nearest-neighbour lookup in the pretrained
embedding. Specialised methods still win on specialised tasks: scVI
for integration, Harmony for batch correction, CellTypist for
cell-type classification under standard ontologies. The pretraining
cost is high — thousands of GPU-hours per model — and the downstream
utility has not yet justified it across most tasks. The biological
interpretability of the learned representations is limited; they do
not straightforwardly recover gene regulatory networks or known
pathways. The GPT analogy is loose. Text tokens carry semantic
meaning individually; gene-expression bins do not. A cell is not a
sentence.

#warn[
  When you read a cell-foundation-model paper, the first three
  questions to ask are: which studies were in the pretraining corpus,
  what evaluation protocol was used (LOCO or LOSO), and what
  specialised baselines were compared against. If the answers are
  "many", "LOCO", and "a trivial baseline", discount the headline
  numbers. The strongest evaluations of cell foundation models so
  far have shown them as competitive with specialised tools on some
  tasks, not as replacements.
]

Where the models might still matter is in places where the
specialised tools also struggle. Cross-species alignment at scale is a
genuine UCE strength because no specialised method exists at the
same coverage. Patient-level predictions from single-cell data —
combining multi-patient data into a shared embedding for disease-state
classification — could parallel what radiomics has done for imaging.
*In-silico* perturbation screens at virtual scale, where a pretrained
representation plus a modest amount of experimental perturbation data
predicts CRISPR-screen outcomes across cell types, is a plausible
near-term application. Multi-modal integration — foundation models
that jointly cover RNA, ATAC, surface protein, and ideally lineage —
is the likely future. Current single-modality cell foundation models
are a stepping stone.


== Pitfalls, in a List <sec:pitfalls>

A working genomics ML project carries a recurring set of design
decisions, and the same set of mistakes recurs in published work. The
list below is not exhaustive, but it captures the failure modes that
the previous five sections of this chapter are organised to prevent.

#figure(
  image("../figures/ch16/f4-pitfall-checklist.svg", width: 95%),
  caption: [
    A design-review checklist for a genomics ML project. The
    five questions on the left map to the failure modes on the right;
    most published controversies in the field reduce to one of them.
  ],
) <fig:checklist>

*Architecture choice.* The first decision is whether the architecture's
inductive bias matches the structural property of the data. The five
canonical pairings in @sec:pairings are the starting point. If your
input is 2-D-image-shaped, the default is a CNN; if it is a long
1-D signal with sparse long-range dependence, the default is dilated
convolution plus attention; if it is a count vector, the default is a
VAE with the right discrete likelihood; if it is a graph, the default
is a message-passing GNN; if it is a 3-D arrangement, the default is
an equivariant transformer. Deviating from the default is fine, but
the deviation should be justified by a property of the data, not by
the convenience of a familiar library.

*Likelihood choice.* The reconstruction or output likelihood is, in
many genomics tasks, more important than the architecture depth. MSE
is the right loss for continuous, roughly Gaussian targets; cross-entropy
for categorical; negative binomial for over-dispersed counts;
zero-inflated negative binomial for over-dispersed counts with excess
zeros; Tweedie for compound Poisson-Gamma; ordinal regression for
ranked targets. Wrong likelihood is a frequent silent killer.

*Split protocol.* The most-discussed pitfall in this chapter. Random
splits leak via homology, locality, relatedness, or batch effects in
almost every genomics setting. Use blocked cross-validation matched
to the correlation structure: cluster-based splits for protein tasks,
chromosome splits for sequence-regulation tasks, kinship-aware
splits for biobank tasks, leave-one-study-out for cross-study
benchmarks.

*Baselines.* A strong specialised baseline is often missing from
foundation-model papers. The fair comparison for a cell foundation
model is not against an MLP from scratch but against scVI plus
CellTypist plus Harmony on the same data, with the same split. If the
specialised stack is comparable or stronger, the foundation model has
not earned the headline.

*Calibration and error analysis.* Reporting a single accuracy number
is not enough. A well-calibrated classifier should produce predicted
probabilities that match empirical frequencies. A model that produces
70 % accuracy on a task where the best baseline produces 65 % may
still be useless if its errors cluster in a clinically important
sub-population. Stratified evaluation across sub-groups, reliability
diagrams, expected calibration error — all standard ML hygiene and
all rarely seen in genomics papers.

*Test-set reuse.* Iterating on a fixed test set inflates the reported
accuracy through implicit hyper-parameter overfitting. The classic
remedy — keep a held-out set that is touched once — is hard to
enforce in a community benchmark, but the CAFA and CASP cultures, with
strict time-locked test sets, are the standard worth aspiring to.

#warn[
  Most published genomics ML controversies trace back to one of the
  five questions on the checklist. "Why doesn't this paper reproduce?"
  has a small number of common answers: leaky split, wrong baseline,
  inflated zero-shot setup, mis-specified likelihood, or
  test-set contamination. Run the checklist before the experiment, and
  rerun it before reading any paper's headline number.
]


== The Multimodal Frontier <sec:multimodal>

The chapter to this point has treated each modality in isolation. The
likely near future is *multimodal foundation models* that operate in
a shared latent space across more than one biological data type. The
engineering recipe is well understood — separate encoders per
modality, a shared transformer backbone, modality-specific heads — and
the architectural pieces are mature. The bottleneck is the data:
matched multi-modal measurements are far rarer than any individual
modality.

#figure(
  image("../../diagrams/lecture-16/12-multimodal.svg", width: 95%),
  caption: [
    The multimodal frontier. Existing models concentrate either on
    single modalities at high parameter counts (ESM-2, Evo) or
    moderate parameter counts at two-to-three modalities (AlphaFold-3,
    cell-state models). The unfilled region — many modalities at
    high parameter counts — is data-limited, not compute-limited.
  ],
) <fig:multi>

Three near-term directions look promising. *Protein-DNA co-design*
combines a protein-LM with a DNA-LM to design both a protein and its
target nucleic acid jointly — a transcription factor with its binding
site, an RNA-binding protein with its preferred RNA motif, a CRISPR
variant with its guide-RNA preference. RFdiffusion-All-Atom (Krishna
et al., 2024) and AlphaFold-3 (Abramson et al., *Nature* 2024)
already do partial versions of this. *Generative biology* moves past
prediction into design — designed proteins (RFdiffusion + ProteinMPNN
+ AlphaFold as critic, the pipeline Chapter 15 walked), designed DNA
and RNA (Evo and its successors), and, more speculatively, designed
cells. The dual-use questions raised in Chapter 15 become
sharper with each generation; the field has not yet settled on
clear norms. *Convergence with mainstream AI* is the longer-term
trajectory: genomics ML used to lag the mainstream by two or three
years, and that gap has closed. The same scaling laws, the same
optimisation tricks, the same architectural progress now arrives in
genomics within weeks of being published elsewhere.

#note[
  DeepVariant (2018) is a useful waypoint for the moment genomics ML
  became indistinguishable from mainstream ML. Before 2018, most
  bioinformatics tools were hand-engineered statistical methods —
  GATK's HaplotypeCaller is essentially an HMM with hand-crafted
  features. DeepVariant showed that an ImageNet-style CNN trained on
  pileup images could match a decade of careful engineering and
  retrain for a new sequencing technology in a few GPU-days. The
  2018 — 2024 window saw the pattern repeat across the field: scVI
  for single-cell, Basset and DeepSEA and then Enformer for
  regulation, AlphaFold for structure, ESM for protein language,
  scGPT and Geneformer for cells. In each case the trigger was the
  same — enough labelled or pretraining data had accumulated that
  the right deep-learning architecture beat the hand-engineered
  baseline.
]


== Summary <sec:ch16-summary>

- Five canonical architecture–data pairings recur across genomics:
  pileups to 2-D CNNs (DeepVariant), long-range regulation to dilated
  convolutions plus self-attention (Enformer), count matrices to
  negative-binomial VAEs (scVI), molecular graphs to message-passing
  GNNs, and protein 3-D structure to SE(3)-equivariant attention
  (AlphaFold-2). Each pairing is a match between the architecture's
  inductive bias and the data's structural property.
- The right likelihood usually matters more than the architecture
  depth. scVI's negative-binomial reconstruction is the move that
  made single-cell modelling work; the depth of the VAE is
  secondary. Wrong likelihoods are a frequent silent failure mode.
- Inductive bias is the matched-filter design problem at scale.
  Symmetries play the role of the signal template; the output
  likelihood plays the role of the noise model. Strong priors are
  sample-efficient on matching problems and fail badly elsewhere.
- The genomics data pyramid is inverted relative to mainstream ML.
  Unlabelled bases are abundant; gold-standard labels are scarce.
  Self-supervised pretraining followed by fine-tuning is the
  dominant strategy for traversing the gap, and the data-efficiency
  gain is one-to-two orders of magnitude.
- Data leakage is catastrophic and widely under-acknowledged.
  Random splits leak via homology, linkage disequilibrium, kinship,
  and cross-study batch effects. Use blocked cross-validation matched
  to the correlation structure: cluster splits for proteins,
  chromosome splits for regulation, kinship-aware splits for
  biobanks, leave-one-study-out for single-cell.
- DNA language models (DNABERT, Nucleotide Transformer, HyenaDNA,
  Evo, Caduceus) are useful but not transformative. They deliver
  modest improvements on downstream regulatory tasks and clear gains
  on cross-species transfer; they have not yet displaced
  specialised architectures on most genomics benchmarks.
- Cell foundation models (Geneformer, scGPT, scFoundation, UCE) are
  in their infancy. Honest evaluation — leave-one-study-out,
  comparison to specialised baselines — places them as competitive
  on some tasks rather than as replacements. The clearest win is
  cross-species cell-type alignment.
- The five recurring pitfalls — mismatched architecture, wrong
  likelihood, leaky split, weak baseline, naive zero-shot — account
  for most genomics ML controversies. Run the design-review checklist
  before the experiment and re-run it before trusting a paper's
  headline number.
- The multimodal frontier (protein-DNA co-design, generative biology,
  cross-modal foundation models) is data-limited rather than
  compute-limited. The architectural pieces are mature; the matched
  multi-modal corpora are not.


== Exercises <sec:ch16-exercises>

#strong[1.] #emph[Architecture matching.] For each of the following
genomics problems, identify the architecture with the closest match
between inductive bias and data: (a) predict the effect of a missense
variant on protein stability; (b) call structural variants from
long-read alignments; (c) impute missing genotypes from a sparsely
genotyped chip; (d) cluster cells in a CITE-seq dataset that has
both RNA counts and surface protein counts; (e) predict the binding
affinity of a small molecule to a protein, given the molecule's SMILES
and the protein's sequence. For each, write one sentence on why the
chosen architecture's prior fits the data.

#strong[2.] #emph[Likelihood swap.] Implement a small VAE (50 lines
of PyTorch) and train it on a subset of the PBMC 3 k single-cell
dataset under three likelihoods: Gaussian (MSE on log-transformed
counts), Poisson, and negative binomial with a per-gene dispersion
parameter. Compare reconstruction error on held-out cells under each
likelihood; produce UMAPs of the latent space and qualitatively
inspect whether known cell types remain separable. Write one
paragraph on why the NB likelihood produces the best results despite
identical architecture.

#strong[3.] #emph[Leakage diagnostic.] Take a publicly available
protein-function-prediction dataset (e.g., Pfam classification on a
held-out fraction of UniProt). Run two splits: a random 80/20 split,
and a CD-HIT 50 % identity cluster split. Train the same shallow
classifier on both. Report the accuracy gap. Predict, before running,
how large the gap will be; reflect on the difference between your
prediction and the measurement.

#strong[4.] #emph[Information content.] Compute the per-column
information content of the TBP-binding site (TATA box) from a curated
set of binding sites (JASPAR provides matrices and FASTA files of
aligned sites for each TF). Compute the total motif information
$sum_i "IC"_i$. Estimate the expected number of false positives if
this motif is scanned against a 3-billion-base human genome at the
threshold $sum_i log_2 (p_(i, b_i) / 0.25) > T$ for $T$ equal to
half the total information. Compare your estimate to a real
genome-wide TBP-motif scan if you can locate one in the literature.

#strong[5.] #emph[DNA-LM fine-tuning.] Take a small pretrained DNA-LM
(DNABERT or Nucleotide Transformer via HuggingFace) and fine-tune it
on the DeepProm promoter-classification dataset. Compare against a
CNN of comparable parameter count trained from scratch on the same
data. Report (a) validation accuracy under the standard random
split, (b) validation accuracy under a chromosome split, (c)
training time and GPU memory for each approach. Comment on whether
the DNA-LM's pretrained representation justifies its cost on this
specific task.

#strong[6.] #emph[Foundation-model audit.] Read the scGPT paper
(Cui et al., 2024, *Nature Methods*). Audit the evaluation: which
benchmarks use leave-one-study-out, which use leave-one-cell-out,
which compare against specialised baselines such as scVI or
CellTypist on the same split. Write a one-page assessment of which
reported numbers are likely to replicate and which are inflated by
the evaluation protocol. Cite specific figures, tables, or
supplementary sections.

#strong[7.] #emph[Architecture for Hi-C.] Sketch an architecture for
predicting Hi-C contact maps at 10 kb resolution from 1 Mb of
surrounding DNA sequence. Specify the inductive biases you would
build in, the output likelihood, the split protocol you would use
for evaluation, and the baselines you would compare against. Justify
each choice from first principles in two or three sentences each.

#strong[8.] #emph[(Open-ended.)] Pick a 2023 — 2025 genomics ML paper
that you find interesting and that reports a striking headline
number. Run the five-question design-review checklist
(@fig:checklist) against the paper. Identify the weakest of the five
links. Propose one experiment that would change your confidence in
the headline result, and predict whether it would raise or lower the
number. Cite the paper.


== Further Reading <sec:ch16-further-reading>

- *Poplin, R., Chang, P.-C., Alexander, D., et al.* (2018). "A universal
  SNP and small-indel variant caller using deep neural networks."
  _Nature Biotechnology_ 36: 983 — 987. The DeepVariant paper. The
  ImageNet-to-genomics transfer in one place.
- *Avsec, Ž., Agarwal, V., Visentin, D., et al.* (2021). "Effective
  gene expression prediction from sequence by integrating long-range
  interactions." _Nature Methods_ 18: 1196 — 1203. The Enformer paper.
  Pair with Kelley (2018, 2020) on Basenji to read the dilated-conv
  lineage.
- *Lopez, R., Regier, J., Cole, M. B., Jordan, M. I., & Yosef, N.*
  (2018). "Deep generative modeling for single-cell transcriptomics."
  _Nature Methods_ 15: 1053 — 1058. scVI; the right-likelihood
  argument made explicit.
- *Jumper, J., Evans, R., Pritzel, A., et al.* (2021). "Highly
  accurate protein structure prediction with AlphaFold." _Nature_
  596: 583 — 589. Read with Chapter 15 of this book for the
  inductive-bias view of the Evoformer and IPA.
- *Theodoris, C. V., Xiao, L., Chopra, A., et al.* (2023). "Transfer
  learning enables predictions in network biology." _Nature_ 618:
  616 — 624. The Geneformer paper. Pair with Cui et al. (2024) for
  scGPT and with Rosen et al. (2024) for UCE.
- *Nguyen, E., Poli, M., Faizi, M., et al.* (2024). "Sequence
  modeling and design from molecular to genome scale with Evo."
  _Science_ 386: eado9336. The current frontier of DNA language
  modelling at scale; useful as a stress-test for the
  honest-evaluation discipline of @sec:dna-lms.
- *Luecken, M. D., Büttner, M., Chaichoompu, K., et al.* (2022).
  "Benchmarking atlas-level data integration in single-cell genomics."
  _Nature Methods_ 19: 41 — 50. The community standard for honest
  cross-study evaluation in single-cell genomics. Read for the split
  protocols even if you skip the methods.
