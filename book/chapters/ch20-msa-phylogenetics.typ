#import "../theme/book-theme.typ": *

= #idx("multiple sequence alignment")Multiple Sequence Alignment, Phylogenetics, and Comparative Genomics <ch:msa-phylogenetics>

#matters[
  Until this chapter the unit of analysis has been a single genome, a
  single transcriptome, a single protein. The next move in
  bioinformatics is to compare. Once you can compare a gene across
  fifty species, you can ask which residues are under #idx("selection")selection, when
  two lineages diverged, whether a non-coding region is a regulatory
  element, and whether a gene family expanded through duplication or
  through speciation. The whole sub-field of evolutionary genomics
  rests on three computational pillars: a *multiple sequence
  alignment* that places homologous residues in the same column, a
  *phylogenetic tree* that organises the sequences by descent, and a
  *substitution model* that scores how plausible the alignment is
  given the tree. Each pillar comes with an honest tractability story
  — exact #idx("MSA")MSA is NP-hard, tree inference is a search over a
  super-exponential space, codon-level selection tests are
  likelihood-ratio statistics with subtle null distributions — and
  modern tooling is the set of heuristics that make the whole stack
  work in minutes rather than centuries.
]

A genome looks different the second time you see it. The first pass
through the assembly, the variant calls, the expression matrix, the
#idx("chromatin")chromatin tracks — the work of Chapters 1 through 19 — treats the
species in isolation. Real biological insight rarely lives inside a
single species. It lives in the differences. A glance at the
haemoglobin alpha chain in human, mouse, chicken, and lamprey tells
you which residues touch the haem iron, which residues touch the
neighbouring beta chain, and which residues are free to wander. A
glance at the same gene in five hundred mammals tells you which
positions accelerated in bats and which positions froze in cetaceans.
A glance across species in a syntenic interval tells you whether the
region you are staring at is a conserved #idx("regulatory element")regulatory element or a
neutral intergenic stretch. The unit of biological insight is not a
genome. It is a *comparison*.

This chapter is the toolbox for that comparison. It opens with
multiple sequence alignment — the unglamorous data structure on which
everything else hangs — and walks through why exact alignment is
intractable, why the field settled on progressive heuristics, and what
the modern tools (#idx("Clustal")Clustal Omega, #idx("MUSCLE")MUSCLE, #idx("MAFFT")MAFFT, #idx("T-Coffee")T-Coffee, ProbCons)
actually do under the hood. It then moves to phylogenetic inference
proper. Distance methods, parsimony, and #idx("maximum likelihood")maximum likelihood are three
distinct philosophies for turning an alignment into a tree, and a
good practitioner reaches for whichever one matches the problem.
#idx("Felsenstein")Felsenstein's pruning algorithm — belief propagation on a tree, seven
years before Pearl named graphical models — is the algorithmic heart
of the field. The chapter then turns calendrical: the #idx("molecular clock")molecular clock,
relaxed clocks, and divergence-time estimation. From there to selection:
the codon-substitution machinery of #idx("dN/dS")dN/dS, the likelihood-ratio
framework for detecting positive selection, and the worked examples
that made the technique famous (HA, MHC, sperm-egg recognition).
The closing sections leave the gene level entirely. #idx("synteny")Synteny aligns
whole genomes, conservation tracks score per-base evolutionary
constraint, and Infernal plus Rfam apply the same machinery to
non-coding RNAs whose function depends on structure rather than on
sequence.

A reader who has done Chapters 2 (#idx("read alignment")read alignment), 4 (variant
calling), and 9 (#idx("peak calling")peak calling) will recognise a recurring move in this
chapter. Each previous problem had an obvious objective function and a
tractable algorithm. This chapter's problems have obvious objective
functions and *intractable* algorithms; every working method is a
disciplined compromise. Phylogenetics is, in that sense, the maturity
test of bioinformatics. The compromises are the methodology.


== From Pairwise to Multiple <sec:pairwise-to-multiple>

The pairwise-alignment algorithms of Chapter 2 — Needleman–Wunsch for
global, Smith–Waterman for local — solved a clean two-string dynamic
programming problem in $O(L^2)$ time. The biological objective they
served was rarely about two strings alone. A biologist who has a
candidate gene almost always wants to align it against every known
homologue: all twenty haemoglobins, all fifty cytochrome $c$ orthologs,
all of #idx("Pfam")Pfam's globin domain. The right data structure for that comparison
is a *multiple sequence alignment* (MSA): a matrix of $N$ rows and
$M$ columns where every column claims a single evolutionary identity.
The residues in a column are either descended from a common ancestral
residue, or they are gaps marking an insertion or deletion in a
specific lineage.

A good MSA is the substrate for almost every comparative move in this
chapter. Phylogenetic trees are inferred from MSAs. Profile HMMs (the
Pfam domain models you met in Chapter 15) are estimated from MSAs.
#idx("AlphaFold-2")AlphaFold-2's MSA features power its co-evolutionary contact
predictions. Conserved-block detection, #idx("motif")motif discovery, codon-substitution
analysis, and ancestral-state reconstruction all start from a column
matrix where homologous residues line up. A bad MSA poisons all of
that downstream work in ways that rarely announce themselves loudly.

#figure(
  image("../../diagrams/lecture-20/01-pairwise-to-msa.svg", width: 95%),
  caption: [
    From pairwise to multiple. The exact-MSA tensor scales as $L^N$
    and explodes past four or five sequences. #idx("progressive alignment")Progressive alignment
    decomposes the tensor into a sequence of pairwise steps along a
    guide tree, paying $O(N^2 L^2)$ instead of $O(L^N)$.
  ],
) <fig:msa-tensor>

=== Why Exact MSA Is Hopeless

The natural generalisation of pairwise #idx("dynamic programming")dynamic programming is the
$N$-dimensional Carrillo–Lipman tensor. Each cell stores the best
alignment of a prefix combination drawn from all $N$ sequences;
each transition has up to $2^N - 1$ neighbours, one per non-empty
subset of sequences that advance. Time and space both scale as $O(L^N)$
in the alphabet's prefix space, which is unhelpfully written but
roughly correct. For five sequences of length three hundred the
state count is around $10^{12}$. For fifty sequences — a small family
— it is something like $10^{125}$. Wang and Jiang (1994) proved that
multiple alignment under the sum-of-pairs score is *NP-hard* in
general; no clever indexing trick rescues the exact problem.

The Carrillo–Lipman bound (Lipman, Altschul, and Kececioglu, 1989) does
permit exact alignment for very small $N$ and short $L$ by pruning
branches of the DP tensor that cannot beat the current best.
The technique still ships inside MSA, the tool, but rarely scales beyond five
sequences. The bound's main use is as a sanity check on heuristic
output for tiny problems, not as a production algorithm.

The standard scoring objective for the exact problem is *sum-of-pairs*:

$ "SP"(A) = sum_(i < j) "score"(A_i, A_j) $

— the sum, over all pairs of rows in the alignment, of the standard
pairwise alignment score with a #idx("substitution matrix")substitution matrix (#idx("BLOSUM")BLOSUM, #idx("PAM")PAM) and
gap penalties. SP is convenient because it factorises across pairs,
but it carries a known bias: identical columns count quadratically in
$N$, which over-weights conservation in dense families. T-Coffee
(Notredame, Higgins, and Heringa, 2000) attacks the same objective
with a *consistency-based* twist: collect pairwise local and global
alignments, build a library of weighted residue pairs, and let the
progressive aligner make its column choices consistent with as much
of the library as possible. Probalign (Roshan and Livesay, 2006) and
ProbCons (Do et al., 2005) take the consistency idea further by
working in posterior-probability space rather than with hard
alignments.

#note[
  Sum-of-pairs is one of several scoring choices the MSA literature
  argues over. *Maximum-weight trace*, *minimum entropy*, and
  *consensus-likelihood* objectives all have advocates. In practice
  the choice matters less than the heuristic search procedure; the
  same tool with the same SP-style score produces dramatically
  different alignments depending on guide-tree construction and
  iterative-refinement budget. Pick a tool, audit its alignment
  visually for the regions you care about, and accept that an MSA is
  always an approximate object.
]


== Progressive Alignment <sec:progressive>

Feng and Doolittle (1987) cut the Gordian knot with one move:
*align in tree order*. Build a rough guide tree from quick pairwise
distances, then align the two closest sequences first; at every
subsequent internal node align the existing profile against the next
sequence (or against another profile). The work shrinks from $O(L^N)$
to $O(N^2 L^2)$ — roughly the cost of $N$ pairwise alignments plus
the guide-tree construction. The tradeoff is that any gap placed
early in the process is locked into the final alignment. Feng and
Doolittle stated the principle as "*once a gap, always a gap*", and
it has shadowed the literature for forty years.

The CLUSTAL family (Higgins and Sharp, 1988; CLUSTAL W in 1994;
CLUSTAL Omega in 2011) was the first widely used implementation.
The recipe is direct. Compute all pairwise distances by $k$-tuple
counting — fast, but inaccurate at distant homology because shared
$k$-tuples become rare. Build a guide tree by UPGMA or neighbour-joining.
Walk the tree from leaves to root; at each node align the two child
profiles using a profile-profile dynamic programme with #idx("affine gap")affine gap
penalties. CLUSTAL Omega uses HHalign's profile-profile algorithm
underneath, which is part of why it scales to tens of thousands of
sequences. The result is a single alignment in column form, ready
for tree inference, #idx("profile HMM")profile #idx("HMM")HMM construction, or visualisation.

#figure(
  image("../../diagrams/lecture-20/02-progressive-msa.svg", width: 95%),
  caption: [
    Five globin sequences aligned progressively along a guide tree.
    Each internal node merges a profile with a sequence or with
    another profile. Conserved positions — including the
    heme-binding histidine — surface as fully populated columns by
    the final alignment.
  ],
) <fig:progressive>

=== Iterative Refinement and the MUSCLE Twist

The "once a gap, always a gap" lock-in led to a generation of tools
that revisited the alignment after the first pass. MUSCLE (Edgar,
*Nucleic Acids Research* 2004) added two innovations that have
become standard. First, MUSCLE *iterates the guide tree*: align,
derive a better distance matrix from the alignment itself, rebuild
the tree, and re-align. The improved tree usually fixes ordering
errors that bias the original alignment. Second, MUSCLE adds
*iterative refinement*: at each step pick a random internal edge of
the tree, bipartition the leaves into two profiles by cutting that
edge, re-align the two profiles, and accept the new alignment only if
the sum-of-pairs score improves. Roughly twenty rounds of this
escapes most of the local minima the first pass falls into.

#tip[
  Iterative refinement is *coordinate-descent on the alignment
  landscape*. Each bipartition picks a coordinate; re-alignment moves
  along it; the SP score is the objective. The landscape is heavily
  non-convex, gap-extension makes it noisy, and the best stopping
  rule is "until the score plateaus for ten consecutive iterations."
  This is the same algorithmic shape as Gibbs sampling, k-means
  re-initialisation, and other coordinate-wise heuristics across ML.
]

=== MAFFT and Fast Fourier Transform Distances

MAFFT (Katoh, Misawa, Kuma, and Miyata, *Nucleic Acids Research*
2002; substantially updated in 2013) is the de facto default for
large MSAs. Its central trick is to compute pairwise distances by
*fast Fourier transform* on physicochemical-property profiles of
each sequence rather than by full pairwise dynamic programming. Each
residue is mapped to a small vector of properties (volume, polarity,
charge); the cross-correlation of two such property sequences is a
proxy for alignment-shift detection, and a length-$L$ FFT computes it
in $O(L log L)$ rather than $O(L^2)$. The full pairwise DP is still
run, but only inside a narrow band around the best-shift peak.

MAFFT ships with several presets that trade speed for accuracy.
*FFT-NS-1* and *FFT-NS-2* are the fast modes appropriate for tens of
thousands of sequences. *L-INS-i* uses local pairwise alignments as
anchors and consistently scores highest on benchmark sets at the cost
of being slower. *G-INS-i* is the global-alignment counterpart for
sequences expected to be alignable end-to-end. The standard advice
in 2025 is: L-INS-i when accuracy matters; FFT-NS-2 when scale matters;
Clustal Omega when you need to align a hundred thousand sequences at
once.

#figure(
  image("../../diagrams/lecture-20/08-msa-benchmark.svg", width: 90%),
  caption: [
    Accuracy of common MSA tools on BAliBASE. Modern tools all reach
    sum-of-pairs scores above $0.85$ on the easier reference sets;
    the gap between fastest and most accurate widens on distant-homology
    benchmarks (BAliBASE references 4 and 5). MAFFT L-INS-i and
    ProbCons sit at the top; FFT-NS-2 and Clustal Omega are the
    practical compromises.
  ],
) <fig:msa-benchmark>

#warn[
  An MSA viewed as a single static object is misleading. Most working
  practitioners audit alignments visually with AliView, JalView, or
  MEGA before doing anything downstream. The first columns and last
  columns are usually noise — sequence-end artefacts that should be
  trimmed. Columns with more than thirty per cent gaps are usually
  trimmed too, with trimAl or Gblocks. The published convention is
  to report which trimming tool was used and at what aggressiveness;
  silent trimming makes reproducibility a nightmare.
]


== The Complexity Wall and Its Workarounds <sec:complexity>

The complexity story behind MSA is worth dwelling on because it
recurs throughout this chapter. Every interesting evolutionary
inference problem — alignment, tree search, divergence-time
estimation, selection detection — has an exact formulation whose
asymptotic cost is exponential or super-exponential in the size of
the input. The progress of the field has been the progress of
disciplined heuristics: which approximations are safe, which
approximations bias the answer, and which approximations are merely
slow.

#figure(
  image("../figures/ch20/f1-msa-complexity.svg", width: 95%),
  caption: [
    Three regimes of MSA complexity. Exact dynamic programming on the
    Carrillo--Lipman tensor scales as $L^N$ and is hopeless past five
    sequences. Progressive alignment is $O(N^2 L^2)$ and dominates
    practice. Profile-profile alignment with a fixed guide tree is
    $O(N L^2)$ amortised and underpins ultrafast methods such as
    Clustal Omega and FAMSA at the cost of accuracy on distant
    homology.
  ],
) <fig:complexity>

The same engineering pattern shows up in tree inference. Searching
all rooted bifurcating trees on $N$ taxa is $(2N - 3)!! / 2^(N - 1)$
topologies — about $7.9 times 10^4$ at $N = 8$, about $2.2 times 10^(20)$
at $N = 20$, astronomically large at $N = 50$. No exact tree search
is feasible past about ten taxa. Modern likelihood tree software
(#idx("RAxML")RAxML, #idx("IQ-TREE")IQ-TREE, FastTree) uses smart local-rearrangement
heuristics that find the optimum reliably in practice while making
no guarantees in theory. The honest framing of every modern
phylogenetics paper is *we found the best tree we could using a
heuristic search*; whether the tree we found is genuinely optimal
is rarely provable.


== Phylogenetic Trees: What They Are and How to Build Them <sec:trees>

A *phylogenetic tree* is a graph whose leaves are observed taxa
and whose internal nodes represent inferred common ancestors. The
edges carry weights — *branch lengths* — that measure
evolutionary distance, typically in expected substitutions per site.
A tree is *rooted* if it specifies which internal node is the
ancestor of everything else; *unrooted* otherwise. A tree of $N$
leaves has $2N - 2$ edges if rooted, $2N - 3$ if unrooted; the
combinatorics of tree space scale as the double factorial above.

Three philosophies dominate the inference of trees from sequence data.
*Distance methods* — UPGMA, neighbour-joining — collapse the MSA
into a pairwise distance matrix and reconstruct a tree from
distances alone. They are fast and frequently good enough.
*Character methods* — parsimony — keep the alignment columns
explicit and pick the tree that requires the fewest substitutions.
*Likelihood methods* — maximum likelihood and Bayesian inference —
assume an explicit substitution model and pick the tree under which
the data are most probable. Modern phylogenetics is dominated by
likelihood; distance methods survive as fast first-pass scaffolds
and parsimony survives as a sanity check.

#figure(
  image("../figures/ch20/f3-tree-methods.svg", width: 95%),
  caption: [
    Four tree-inference paradigms compared on the same axes. UPGMA
    assumes the molecular clock and is fast. Neighbour-joining drops
    the clock assumption and runs in $O(N^3)$. Maximum parsimony
    minimises total substitutions but is statistically inconsistent
    under #idx("long-branch attraction")long-branch attraction. Maximum likelihood and Bayesian
    inference are statistically consistent under their model
    assumptions and have become the modern defaults.
  ],
) <fig:tree-methods>

=== UPGMA: Hierarchical Clustering with the Clock

UPGMA — Unweighted Pair Group Method with Arithmetic Mean — is
agglomerative hierarchical clustering applied to a pairwise distance
matrix. Find the closest pair, merge them into an internal node with
branch length equal to half their distance, recompute distances from
the merged node to all remaining nodes by averaging, repeat until
one node remains. The algorithm runs in $O(N^2)$ time once the
distance matrix is in hand.

UPGMA's hidden assumption is the *molecular clock*: all lineages
accumulate substitutions at the same rate, so distance is
proportional to time. When the clock holds, UPGMA recovers the
correct rooted tree. When the clock fails — and it usually does —
UPGMA places long branches incorrectly because it cannot tell the
difference between an early divergence at the clock rate and a
recent divergence at an accelerated rate. The tree is biased toward
ultrametric topologies even when the underlying evolution is not.
Use UPGMA for fast scaffolds and for visualising distance matrices.
Do not use it as a publishable phylogeny.

=== Neighbour-Joining: Distance Without the Clock

Neighbour-joining (Saitou and Nei, 1987) drops the molecular-clock
assumption while keeping the distance-matrix simplicity. For each
pair $(i, j)$ it computes a corrected criterion

$ Q(i, j) = (N - 2) d(i, j) - sum_(k) d(i, k) - sum_(k) d(j, k) $

— the raw pairwise distance penalised by each taxon's average
distance to everything else. Pairs joined early are the ones whose
distance is small *relative to their distance to the rest*, which
removes the long-branch bias UPGMA suffers from. NJ runs in $O(N^3)$
time and is provably consistent when the input distances are
additive (i.e., when the distances genuinely came from some tree).
Real distances are rarely exactly additive — they are estimated from
finite-length alignments — but NJ is a robust heuristic that often
produces topologies indistinguishable from maximum-likelihood
results.

#figure(
  image("../../diagrams/lecture-20/03-nj-tree.svg", width: 90%),
  caption: [
    Neighbour-joining on a six-taxon distance matrix. The Q-matrix
    correction penalises each pair by the taxa's average distance to
    everything else, neutralising the long-branch bias that biases
    UPGMA. The resulting tree has lineage-specific branch lengths
    and matches the true topology when the inputs are additive.
  ],
) <fig:nj>

In practice neighbour-joining is the workhorse for "I just want a
tree". #idx("BLAST")BLAST followed by ClustalW followed by NJ is still the most
common workflow for the back-of-the-envelope phylogeny of a small
family.

=== Substitution Models: Jukes--Cantor to GTR

Before maximum likelihood we need a model of how substitutions
accumulate along a branch. The simplest model — Jukes--Cantor
(1969) — assumes all four bases are equiprobable at equilibrium and
all twelve directed substitution rates are equal. Under JC the
probability of a substitution after time $t$ at rate $mu$ is
$3/4 (1 - exp(-4 mu t / 3))$, and the corrected pairwise distance is

$ d_("JC") = -3/4 ln(1 - 4/3 p) $

where $p$ is the observed fraction of mismatched sites. The JC
correction inflates raw distances to account for multiple
substitutions at the same position, which become common past about
twenty per cent divergence.

Kimura's two-parameter model (K80, 1980) allows transitions
(purine-to-purine, pyrimidine-to-pyrimidine) and transversions
(purine-to-pyrimidine and back) to occur at different rates,
capturing the empirical observation that transitions outnumber
transversions roughly two-to-one. The Felsenstein 81 (F81), HKY85
(Hasegawa, Kishino, and Yano), and GTR (general time-reversible)
models add parameters for unequal base frequencies and unequal
substitution rates between specific base pairs. GTR is the most
permissive *reversible* model — twelve rate parameters constrained
by detailed balance leave six free — and is the default for
#idx("nucleotide")nucleotide phylogenetics. Protein analogs (JTT, WAG, LG) replace the
$4 times 4$ rate matrix with a $20 times 20$ matrix estimated from
empirical protein databases.

#figure(
  image("../figures/ch20/f2-substitution-models.svg", width: 95%),
  caption: [
    The nested ladder of substitution models. Jukes--Cantor sits at
    the bottom with one rate; Kimura's K80 splits transitions from
    transversions; F81 and HKY85 add unequal base frequencies; GTR
    adds six free rates under detailed balance. Each step adds
    parameters that better fit empirical data at the cost of
    estimation variance.
  ],
) <fig:sub-models>

#note[
  Site-rate heterogeneity matters as much as the choice of base
  model. Real positions evolve at very different rates: third #idx("codon")codon
  positions and rRNA loops fast, second codon positions and rRNA
  stems slow. The standard fix is to add a Gamma-distributed
  site-rate ($+G$) model, usually with four discrete rate categories,
  and a fraction of *invariant* sites ($+I$) for positions that are
  effectively frozen. Almost every modern phylogenetic analysis
  reports the model as "GTR+G+I" or similar; the model-selection
  step in IQ-TREE or jModelTest auto-picks the best fit by AIC or
  BIC.
]

=== Maximum Parsimony and the Long-Branch Trap

Parsimony picks the tree that requires the fewest substitutions to
explain the data. For a fixed tree, Fitch's algorithm (1971)
computes the minimum substitution count column by column in $O(N L)$
time; Sankoff's generalisation (Sankoff and Cedergren, 1983) handles
arbitrary substitution costs. Tree search is the hard part —
$(2N - 3)!! / 2^(N - 1)$ topologies is too many to enumerate past
ten taxa — so parsimony software uses local rearrangements
(nearest-neighbour interchange, subtree pruning and regrafting) to
hill-climb on parsimony score.

Parsimony was the dominant method through the 1980s and 1990s. It
fell out of favour after Felsenstein (1978) proved it is
*statistically inconsistent* under what came to be called *long-branch
attraction*. When two long branches have evolved fast enough to
accumulate convergent substitutions at random, parsimony sees the
shared substitutions as evidence of common ancestry and groups the
two long branches together — even when the true tree separates them.
The error gets *worse* with more data, because more sites accumulate
more convergent substitutions. Maximum-likelihood and Bayesian
methods, with proper substitution models, are statistically
consistent under the same conditions.

#figure(
  image("../../diagrams/lecture-20/11-lba.svg", width: 90%),
  caption: [
    Long-branch attraction. The true tree pairs $A$ with $B$ and $C$
    with $D$, but the long branches leading to $A$ and to $C$
    accumulate convergent substitutions by chance. Parsimony groups
    $A$ with $C$ — wrong — while likelihood methods using a proper
    substitution model recover the true topology.
  ],
) <fig:lba>

#warn[
  Long-branch attraction is the cautionary tale of phylogenetics.
  *Adding more sequence does not save you* — under LBA, more data
  drives the wrong answer with higher confidence. The fix is either
  to use likelihood with a model that captures the rate variation,
  or to break the long branches by sampling more intermediate taxa.
  When a textbook tree looks surprising — flatworms grouped with
  amoebozoans, microsporidia grouped with fungi or with amoebae —
  the first hypothesis to test is LBA.
]

=== Maximum Likelihood: Felsenstein's Pruning Algorithm

Maximum likelihood fits a tree by maximising the probability of the
observed alignment given the tree topology, branch lengths, and
substitution model. For a tree with $N$ leaves and $L$ alignment
columns, the data likelihood factorises across columns:

$ L("tree") = product_(c = 1)^L P("column"_c | "tree", "model") $

Each column's contribution is computed by Felsenstein's pruning
algorithm (Felsenstein, *Journal of Molecular Evolution* 1981). At
every leaf, set the per-state likelihood vector to a one-hot indicator
of the observed residue. At every internal node $v$ with children
$u$ and $w$ and edge lengths $t_u$ and $t_w$, compute

$ L_v(x) = (sum_y P(y | x, t_u) L_u(y)) (sum_z P(z | x, t_w) L_w(z)) $

— a sum over child states marginalised against the transition
probabilities given the edge lengths and the substitution model.
The recursion runs from leaves up to the root in $O(N |Sigma|^2)$
time per column, where $|Sigma|$ is the alphabet size. The root
likelihood is then averaged against the equilibrium distribution and
log-summed across columns.

#figure(
  image("../../diagrams/lecture-20/09-felsenstein.svg", width: 90%),
  caption: [
    Felsenstein's pruning algorithm on a five-leaf tree. Each
    internal node holds a vector of state-likelihoods marginalised
    over the subtree below it. The bottom-up sweep is exactly
    forward-pass belief propagation on a tree, seven years before
    Pearl named graphical models.
  ],
) <fig:felsenstein>

The algorithmic shape is familiar from electrical engineering. A
phylogenetic tree with substitution model is a *probabilistic
graphical model*; internal nodes are hidden, leaves are observed,
edges carry conditional transition probabilities, and Felsenstein's
pruning is *belief propagation* — the bottom-up forward pass that
computes the marginal likelihood at the root. The Bayesian network
literature describes the same algorithm as the *sum-product*
algorithm on a tree-structured graph. Felsenstein invented it for
the phylogenetic special case in 1981; Pearl gave it a general
treatment in *Probabilistic Reasoning* in 1988. The seven-year gap is
one of the more remarkable independent rediscoveries in modern
computational science.

Felsenstein's *pulley principle*, also from the 1981 paper, observes
that for any time-reversible substitution model the choice of root
on an unrooted tree does not affect the likelihood — sliding the
root along any branch shifts time forward in one subtree and back in
the other in a way that exactly cancels in the transition
probabilities. The principle means likelihood searches operate on
unrooted topologies; rooting is decided separately, usually by an
outgroup or by a relaxed-clock analysis.

Modern likelihood software — RAxML (Stamatakis, 2014), IQ-TREE
(Nguyen, Schmidt, von Haeseler, and Minh, 2015), FastTree (Price,
Dehal, and Arkin, 2010) — finds the maximum-likelihood tree by
combining Felsenstein's pruning for likelihood evaluation with smart
local-rearrangement heuristics for tree search. IQ-TREE is the
modern default for medium-sized problems (a few thousand taxa, a few
megabases of alignment) and includes a model-selection step
(ModelFinder) that auto-picks the best substitution model by BIC.

Bayesian phylogenetics — #idx("MrBayes")MrBayes (Huelsenbeck and Ronquist, 2001),
BEAST (Drummond and Rambaut, 2007) — replaces the maximum-likelihood
point estimate with a posterior distribution over trees. The MCMC
sampler explores tree space, branch-length space, and model-parameter
space simultaneously; the posterior on each clade is the fraction of
sampled trees that contain it. Bayesian posteriors give honest
uncertainty quantification on every internal branch, at the cost of
substantial extra compute. BEAST is the tool of choice when
divergence-time estimation is the goal; the relaxed-clock machinery
of @sec:clock lives there natively.

=== Bootstrapping and Branch Support

Felsenstein (*Evolution*, 1985) introduced *non-parametric
bootstrapping* as the standard way to express uncertainty about an
ML tree's branches. Resample the columns of the alignment with
replacement to form a pseudo-replicate of the same length; rebuild
the tree on the pseudo-replicate; repeat one hundred to one thousand
times; for each branch in the original tree, record the fraction of
bootstrap trees that contain the same clade. Branches with bootstrap
support of seventy or higher are usually considered reliable;
branches below fifty are usually collapsed into polytomies. RAxML
implements *rapid bootstrap* and IQ-TREE implements *ultrafast
bootstrap*, both of which exploit alignment-likelihood approximations
to give comparable results in a fraction of the time.

#tip[
  Bootstrap proportions are *not* posterior probabilities. The
  bootstrap measures resampling consistency, which is related to but
  not identical to clade confidence. A heavily supported clade
  under bootstrap may still be wrong if the model is misspecified.
  The Bayesian alternative — posterior clade probabilities from
  MrBayes or BEAST — provides a different (and arguably more
  interpretable) confidence measure, with its own quirks around
  prior choice and MCMC convergence.
]


== The Molecular Clock <sec:clock>

The branch lengths of a phylogenetic tree are measured in expected
substitutions per site. They are not, by default, calibrated in
years. The *molecular clock hypothesis* (Zuckerkandl and Pauling,
1965) says they ought to be: substitutions accumulate at a roughly
constant rate per lineage per unit time, so a tree with branch
lengths in substitutions can be re-scaled to a tree with branch
lengths in millions of years if we know the rate.

The hypothesis is correct in spirit and wrong in detail. Rates vary
across lineages (rodents evolve faster than primates), across
genes (mitochondrial #idx("DNA")DNA accumulates substitutions about ten times
faster than nuclear), across regions of a genome (synonymous codon
positions faster than non-synonymous), and across time within a
single lineage (rate slowdowns are common after large radiations).
A *strict* clock that assumes a single rate is rejected by every
real dataset on inspection.

=== The Poisson Process View

The cleanest mathematical statement of the clock concept is that
substitutions accumulate as a *Poisson process* with rate $mu$ per
site per unit time. Over time $t$, the number of substitutions at a
site is $X(t) tilde "Pois"(mu t)$, with mean and variance both equal
to $mu t$. The variance-equals-mean property is testable on real
data: under a strict clock, the variance of substitution counts
across sites in a single lineage should equal the mean. *Overdispersion*
— variance much larger than mean — is the empirical signature that
sites are evolving at different rates, which is the motivation for
the gamma-distributed site-rate corrections of @sec:trees.

#figure(
  image("../../diagrams/lecture-20/04-clock-dating.svg", width: 90%),
  caption: [
    Molecular clock dating. The chimpanzee--human split (six million
    years ago, calibrated from fossils) anchors the rate; branch
    lengths in substitutions translate to dates with credible
    intervals under a #idx("relaxed clock")relaxed clock. The mouse--rat split lands at
    about thirteen million years; the rodent--primate split at about
    eighty million years.
  ],
) <fig:clock>

=== Relaxed Clocks and Calibration

Because the strict clock fails, modern divergence-time estimation uses
*relaxed clock* models that allow rates to vary along the tree.
*Uncorrelated relaxed clocks* draw each branch's rate from a
lognormal or exponential distribution; *autocorrelated relaxed
clocks* let parent and child rates be correlated according to a
random-walk model. BEAST implements both flavours; the choice
between them is usually decided by Bayes-factor comparison or by
posterior-predictive checks.

Converting branch lengths from substitutions to years requires
*calibration points*: known divergence dates from the fossil record
or from independent biogeographic evidence. The chimpanzee–human
split at roughly six to seven million years ago is the canonical
mammalian calibration; the actinopterygian–sarcopterygian split at
around four hundred million years anchors deeper vertebrate dates.
Calibrations are usually specified as prior distributions on internal
node ages (lognormal, exponential with offset, uniform); BEAST
combines them with the relaxed clock and the substitution model into
a posterior over the entire dated tree.

A modern application of the dated phylogeny machinery — and the one
most readers will recognise from recent news — is the
*molecular-clock dating of pathogen outbreaks*. The dating of
SARS-CoV-2's most-recent common ancestor to roughly October to
December 2019 came from a relaxed-clock analysis of early genome
sequences combined with sampling-time calibrations. The same
machinery is the standard for influenza A pandemic dating, HIV
phylogeography, and the rate-calibrated outbreak reconstruction that
Nextstrain ships in production.

#note[
  The clock-dating literature has a recurring tension between
  *node-age priors* (specified by the analyst from fossil evidence)
  and *fossilised birth--death models* (which infer node ages from
  the diversification process itself). The latter are more honest
  about uncertainty but require more compute and more careful prior
  specification. Stadler's family of FBD models is the current
  state of the art for combining sparse fossil data with extant
  molecular sequences.
]


== dN/dS: Selection in the Codon Substitution Matrix <sec:dnds>

The chapter so far has treated sequences as strings of nucleotides
or amino acids. For protein-coding regions a finer-grained model
pays off: substitutions in coding DNA come in two flavours.
*Synonymous substitutions* preserve the encoded amino acid by
changing a codon to a synonym (most often at the third codon
position); they are nearly neutral with respect to protein function.
*Non-synonymous substitutions* change the encoded amino acid and
are subject to selection — purifying selection removes the
deleterious ones, positive selection preserves the beneficial ones.
The ratio of the two rates is a quantitative measure of selection.

=== Definitions and the $omega$ Ratio

Define $d_S$ as the rate of synonymous substitutions per synonymous
site, and $d_N$ as the rate of non-synonymous substitutions per
non-synonymous site. The two normalisations matter: a randomly
chosen codon position has roughly one synonymous site and roughly two
non-synonymous sites under the standard code, and the per-site rates
correct for that. The dN/dS ratio is

$ omega = d_N / d_S $

with three regimes:

- $omega = 1$: neutral evolution. Synonymous and non-synonymous
  substitutions accumulate at the same per-site rate, which is what
  you would see under no selection at all.
- $omega < 1$: *purifying selection*. Non-synonymous changes are
  removed faster than synonymous ones, which is the case for the
  overwhelming majority of genes most of the time. Typical
  whole-gene $omega$ in mammals is $0.1$ to $0.3$; conserved genes
  like histones reach $omega$ near $0$.
- $omega > 1$: *positive selection*. Non-synonymous changes are
  retained faster than synonymous ones, which requires an active
  selective advantage. Whole-gene $omega > 1$ is rare; site-specific
  $omega > 1$ at a handful of residues is common in fast-evolving
  immune and reproductive genes.

#figure(
  image("../../diagrams/lecture-20/05-dnds.svg", width: 95%),
  caption: [
    Distributions of $omega$ across genes under three regimes.
    Purifying selection (histones, ribosomal proteins) concentrates
    mass near zero. Neutral evolution (most pseudogenes) concentrates
    mass around one. Positive selection (influenza HA, MHC,
    sperm-egg-recognition proteins) puts a fraction of sites above
    one while the bulk of the gene remains under purifying selection.
  ],
) <fig:dnds>

=== Computing dN and dS

The Nei--Gojobori (1986) method counts synonymous and non-synonymous
sites and substitutions directly from a pairwise codon alignment
under the assumption that all substitution types are equally likely.
The method is fast and reasonable as a first pass but biased in
realistic data where transitions outnumber transversions. The
Yang--Nielsen (2000) family of codon-substitution models — implemented
in #idx("PAML")PAML's `codeml` — replace the counting with a full
maximum-likelihood model of codon substitution. The model is
$61 times 61$ (excluding stop codons) and parameterises codon
substitutions by base substitution rates, codon usage frequencies,
transition–transversion ratio, and the $omega$ ratio itself.
Likelihoods are computed by Felsenstein's pruning algorithm on the
codon-state space.

=== Likelihood-Ratio Tests for Selection

The standard way to detect positive selection is by *likelihood-ratio
test* (LRT). Fit two nested models to the same data:

- $M_0$: a single $omega$ for the whole gene across all sites.
- $M_1$: two site classes, one with $omega < 1$ (constrained) and
  one with $omega = 1$ (neutral).
- $M_2$: as $M_1$ but with a third class allowed $omega > 1$
  (positive selection).
- $M_(7), M_(8)$: a beta-distributed $omega$ over $(0, 1)$, with $M_8$
  adding a discrete class with $omega > 1$.

The LRT statistic is

$ Lambda = 2 (ln L_("alt") - ln L_("null")) tilde chi^2_(d f) $

with degrees of freedom equal to the difference in free parameters
between the two models. A significant LRT in $M_1$ vs $M_2$ (or
$M_7$ vs $M_8$) is the field's standard evidence for
positive selection. PAML's `codeml` ships with all of these as
preset model labels.

#note[
  The Bayes-Empirical-Bayes step in PAML reports posterior
  probabilities that specific sites belong to the positively-selected
  class. These per-site labels are the headline output of a
  selection scan; conservation track followups (mapping the
  positively selected residues onto a known structure) are the
  standard validation. The 2005 Yang--Wong--Nielsen paper is the
  reference for the BEB procedure.
]

=== The Canonical Examples

A short list of genes for which positive selection has been
repeatedly documented:

- *Influenza HA*. The haemagglutinin glycoprotein on the influenza
  virus surface has a small handful of antigenic-region residues
  with $omega gt.tilde 5$; these are the positions that change
  most rapidly between seasons and that drive the annual vaccine
  redesign. Bush et al. (*Science* 1999) is the textbook reference.
- *MHC class I and II*. The peptide-binding groove residues are
  highly polymorphic within populations and rapidly diverging
  between species. Hughes and Nei (1988) was the first paper to
  apply codon-level selection tests to MHC and is part of why dN/dS
  has its modern profile.
- *Sperm-egg-recognition proteins*. Bindin (sea urchin), ZP3
  (mammalian zona pellucida), and the lysin-VERL pair (abalone) all
  show elevated $omega$ on the gamete-interaction surfaces.
  Speciation-related selection on reproductive isolation is the
  standard interpretation.
- *Antiviral restriction factors*. TRIM5alpha, APOBEC3G, and other
  primate antiviral genes show $omega > 1$ at interface residues
  with viral capsid or genome proteins. Sawyer et al. (2005) is the
  reference for TRIM5alpha.
- *Histones, ribosomal proteins, actin*. The other end of the
  distribution: $omega$ near $0$ over hundreds of millions of years
  of evolution. These are the proteins where every amino acid
  matters.

#warn[
  dN/dS is not a magic positive-selection detector. The LRT can fire
  on technical artefacts: bad alignment around indels inflates
  $omega$ on a handful of sites; #idx("recombination")recombination breaks the
  single-tree assumption; saturation of $d_S$ at deep divergence
  makes the ratio unstable. The 2005--2015 wave of "every other
  gene is under positive selection" papers has been substantially
  walked back; current practice requires high-quality codon
  alignments, multiple substitution models, and conservative LRT
  thresholds.
]


== Comparative Genomics: Synteny, Orthology, Conservation <sec:comparative>

Phylogenetics has been gene-level so far. Comparative genomics asks
the same questions at *whole-genome* scale: how is gene order
preserved across species, which genomic regions are under purifying
selection, which species share a chromosomal rearrangement that
suggests a more recent common ancestor.

=== Synteny: Conserved Gene Order

*Synteny* is the preservation of gene order along a #idx("chromosome")chromosome
across species. Conserved synteny blocks are stretches of the
ancestral chromosome that have not been broken by inversion,
translocation, or fission in the time since the last common
ancestor. Human chromosome 17 and mouse chromosome 11 share a large
syntenic block — the HOXB cluster, p53 region, #idx("BRCA1")BRCA1 neighbourhood —
that has stayed intact across roughly ninety million years of
mammalian evolution. The block is interrupted by a handful of small
inversions and by lineage-specific gene gains and losses.

Synteny is detected by aligning two genomes (or many genomes
pairwise) at the level of homologous gene clusters. Mauve (Darling
et al., 2010) was the classic tool for bacterial-genome synteny via
*locally collinear blocks* (LCBs); MUMmer/nucmer (Kurtz et al., 2004)
handles arbitrary genome pairs with suffix-tree-based seed and
extend; modern tools (LAST, #idx("minimap2")minimap2, SyMAP) scale to gigabase
mammalian genomes. The output is a list of synteny blocks with
start, end, and orientation in each genome.

#figure(
  image("../../diagrams/lecture-20/06-synteny.svg", width: 95%),
  caption: [
    Synteny dot plot for human chromosome 17 against mouse chromosome
    11. The diagonal stripes are conserved synteny blocks — long
    runs of homologous genes that have stayed in the same order
    since the human--mouse common ancestor. The off-diagonal
    fragments mark inversions and small translocations.
  ],
) <fig:synteny>

A synteny dot plot is the standard visualisation. Place genome 1
along the $x$-axis, genome 2 along the $y$-axis, and put a dot at
every pair of homologous gene positions. Diagonal stripes are
collinear blocks. Off-diagonal stripes (inverted slope) are
inversions. Scattered dots are individual gene moves. Whole-genome
synteny browsers (Genomicus, Synima, the UCSC synteny tracks)
present the same information as ribbons rather than as dots.

=== Orthologs and Paralogs

The vocabulary of homology, made precise by Fitch (*Systematic
Zoology*, 1970), distinguishes two cases. *Orthologs* are
homologous genes in different species whose most recent common
ancestor was a *speciation* event — "the same gene in different
species". *Paralogs* are homologous genes whose most recent common
ancestor was a *duplication* event — "different genes in the same
family". The distinction matters for function: orthologs typically
retain the ancestral function across species, while paralogs often
sub-functionalise or neo-functionalise after duplication.

#figure(
  image("../../diagrams/lecture-20/10-ortho-paralog.svg", width: 95%),
  caption: [
    Orthologs versus paralogs in a small gene-family tree. A
    duplication node produced the $alpha$ and $beta$ sub-families.
    Subsequent speciation produced the human, mouse, and dog
    orthologs within each sub-family. Human-$alpha$ and mouse-$alpha$
    are orthologs; human-$alpha$ and human-$beta$ are paralogs;
    human-$alpha$ and mouse-$beta$ are sometimes called *out-paralogs*
    or, when the duplication preceded a whole-genome duplication
    event, *ohnologs*.
  ],
) <fig:ortho-paralog>

Practical #idx("ortholog")ortholog inference comes in three flavours. *Reciprocal
best hits* (RBH) is the simplest: BLAST each gene against the other
species; if gene $A$ in species 1 is the best hit of gene $B$ in
species 2 and vice versa, call them orthologs. RBH is fast and
reaches roughly eighty per cent accuracy on benchmarks; its main
failure modes are gene loss (an ortholog has been lost in one
species, so RBH picks the closest #idx("paralog")paralog) and recent duplication
(the closest paralog and the true ortholog have similar BLAST
scores). *Tree-based methods* — #idx("OrthoFinder")OrthoFinder (Emms and Kelly, 2015,
2019) being the canonical example — build a gene tree per family,
identify speciation and duplication nodes by reconciliation against
the species tree, and read orthologs and paralogs from the
reconciliation. OrthoFinder scales to hundreds of species and is the
modern default for *de novo* ortholog inference. *Synteny-based
methods* (OMA, MetaPhOrs) add gene-neighbourhood evidence and tend
to outperform sequence-only methods on recently duplicated families.

#tip[
  When in doubt, check OrthoDB or Ensembl Compara before doing your
  own ortholog inference. These pre-computed databases use the most
  expensive tree-based methods and curators have manually validated
  the boundaries of the major gene families. For a typical
  mammalian gene the right answer is already there.
]

=== Conservation Tracks

The other product of genome-scale alignment is a per-base
conservation score. *PhyloP* (Pollard et al., 2010) and
*PhastCons* (Siepel et al., 2005) are the two standard UCSC tracks.
PhyloP scores each base by the log of the likelihood ratio between
the neutral substitution rate and the observed rate across a
multi-species alignment; bases with high positive PhyloP have changed
slower than expected (purifying selection) and bases with high
negative PhyloP have changed faster (lineage-specific acceleration).
PhastCons identifies conserved *segments* under a phylo-HMM that
flips between fast and slow evolutionary regimes along the genome.

The tracks are workhorses for variant prioritisation in clinical
genomics — a missense variant at a PhyloP $>= 5$ position is
substantially more likely to be functional than a variant at a
PhyloP $approx 0$ position. They are also the foundation of
*phylogenetic footprinting*: a non-coding region that is conserved
across mammals is a strong candidate regulatory element even without
#idx("ChIP-seq")ChIP-seq data, and intersecting conservation tracks with ChIP-seq
peaks and motif scans gives the highest-confidence
transcription-factor binding-site inventory.

#figure(
  image("../../diagrams/lecture-20/12-tree-of-life.svg", width: 95%),
  caption: [
    A schematic excerpt of the tree of life. Modern tree-of-life
    work (Hug et al., *Nature Microbiology* 2016) uses concatenated
    marker genes across thousands of taxa and a relaxed clock for
    dating. The most consequential modern re-ordering is the
    placement of archaea as sister to eukaryotes rather than as a
    sister clade to bacteria.
  ],
) <fig:tol>


== Non-Coding #idx("RNA")RNA: Structure-Aware Alignment <sec:rna>

Most of this chapter has treated sequence as the unit of homology.
For non-coding RNAs that assumption breaks. The function of tRNA,
rRNA, riboswitches, microRNA precursors, and many long non-coding
RNAs depends on a *#idx("secondary structure")secondary structure* of base-paired stems and
loops that is conserved across enormous evolutionary distances even
as the underlying primary sequence diverges. Two tRNAs with twenty
per cent sequence identity may have identical cloverleaf structure.
Profile HMMs and standard MSA tools, which work on per-column
sequence statistics, cannot represent the constraints introduced by
non-local base-pairing. The right machinery is *stochastic
context-free grammars*.

=== RNA Secondary Structure Basics

RNA secondary structure is the set of Watson--Crick (A-U, G-C) and
#idx("wobble")wobble (G-U) base pairs formed by a single-stranded RNA molecule
folding back on itself. The most common structural motifs are
*hairpin loops* (a stem followed by a loop), *internal loops* (two
stems separated by an unpaired bulge), *multi-branch loops* (three
or more stems meeting at a common loop), and *pseudoknots* (two
stems whose base-pairs interleave). Pseudoknots are biologically
important but algorithmically painful: classical RNA-folding
algorithms exclude them because the minimum-free-energy problem with
pseudoknots is NP-hard (Lyngsø and Pedersen, 2000).

The standard objective for RNA folding is *minimum free energy*
under a thermodynamic parameter set (Turner energies). Zuker's
algorithm (Zuker and Stiegler, 1981) computes the MFE structure in
$O(L^3)$ time and $O(L^2)$ space by dynamic programming on substring
intervals. McCaskill (1990) extended the same machinery to compute
the *partition function* and per-base-pair probabilities, which gives
a probabilistic measure of structural reliability. ViennaRNA (Lorenz
et al., 2011) is the standard implementation; `RNAfold` produces
both the MFE structure (in dot-bracket notation) and the partition
function.

#figure(
  image("../../diagrams/lecture-20/13-rna-structure.svg", width: 95%),
  caption: [
    tRNA secondary structure as a cloverleaf, with the acceptor stem,
    D-arm, anticodon arm, and T$Psi$C arm marked. The same structure
    expressed in dot-bracket notation aligns each character of the
    sequence with a structural role. The bottom panel shows the
    profile covariance model — paired columns become paired states,
    loop columns become singlet states, and the SCFG carries the
    base-pair constraints into the alignment.
  ],
) <fig:rna>

=== Stochastic Context-Free Grammars and Infernal

A *stochastic context-free grammar* (SCFG) over RNA assigns
probabilities to derivations rather than just to strings. The grammar
contains rules of the form $W -> a W' b$ (a paired state emitting
base $a$ on the left, base $b$ on the right, with the rest of the
sub-derivation in between) and $W -> a W'$ (a singlet state emitting
a single base). SCFG parsing — the Cocke--Younger--Kasami algorithm
— computes the optimal derivation in $O(L^3)$ time, and the
inside--outside algorithm computes posterior pair probabilities in
the same complexity. Eddy and Durbin's 1994 paper formalised the
*covariance model* (CM) as the RNA analog of the profile HMM: each
column of a structure-annotated MSA becomes either a paired state
or a singlet state, with emission distributions estimated from the
training alignment.

*Infernal* (Eddy, *Bioinformatics* 2013) is the SCFG counterpart of
#idx("HMMER")HMMER. `cmbuild` constructs a CM from a structure-annotated seed
alignment; `cmsearch` scans a sequence database for hits to the
model with statistical significance reported as an #idx("E-value")E-value.
Crucially, Infernal scores both sequence similarity *and*
base-pairing covariation — pairs of positions in the model that
co-vary in a way consistent with maintaining base-pairing are
informative even when the individual positions are not conserved.
This is why CMs detect tRNA, rRNA, and snoRNA orthologs at sequence
identities (sometimes thirty per cent) where profile HMMs fail.

*Rfam* (Kalvari et al., *Nucleic Acids Research* 2021) is the Pfam
analog for RNA: roughly four thousand curated families covering
tRNAs, rRNAs, snoRNAs, miRNA precursors, long non-coding RNAs, and
bacterial riboswitches, each with a seed alignment, consensus
secondary structure, and an Infernal CM. Running `cmsearch` against
the Rfam library is the standard first pass for finding non-coding
RNAs in a new genome, and it is one of the few cases where the
right answer is genuinely a thirty-year-old algorithm: covariance
models from 1994 still outperform deep-learning RNA classifiers on
most tasks.

#tip[
  The CYK / inside--outside parsing of SCFGs is to context-free
  grammars what forward--backward is to HMMs. The complexity
  hierarchy — regular for HMMs, context-free for SCFGs,
  context-sensitive for pseudoknots — parallels the Chomsky hierarchy
  exactly. Pseudoknots require context-sensitive parsing, and the
  best polynomial-time pseudoknot algorithms (Akutsu, 2000; Rivas
  and Eddy, 1999) accept only restricted classes of pseudoknots in
  exchange for tractability.
]

=== #idx("tertiary structure")Tertiary Structure and the AlphaFold-3 Coda

Beyond secondary structure, predicting an RNA molecule's full
three-dimensional coordinates is the structural-RNA analog of the
protein-structure problem of Chapter 15. Pre-deep-learning tools
(RNAComposer, 3dRNA) assembled tertiary structure from
fragment libraries indexed by secondary-structure motifs and reached
reasonable accuracy on small structured RNAs. The deep-learning wave
arrived with RoseTTAFold-NA (Baek et al., 2024) and the RNA
extensions of #idx("AlphaFold")AlphaFold; AlphaFold-3 (Abramson et al., 2024) is the
current frontier and predicts protein, DNA, RNA, and small-molecule
ligand coordinates jointly. Accuracy on small structured RNAs (tRNAs,
riboswitches) is in the $2$–$4$ Å backbone-RMSD range; large
flexible lncRNAs remain hard.


== A Standard Workflow and Its Pitfalls <sec:workflow>

A working phylogenetic study in 2025 follows a standard pipeline.
The names of the tools change every five years; the sequence of
operations does not.

#figure(
  image("../../diagrams/lecture-20/07-workflow.svg", width: 95%),
  caption: [
    The standard phylogenetics workflow. Collect homologs by BLAST or
    OrthoFinder; align with MAFFT L-INS-i; trim with trimAl or
    Gblocks; run IQ-TREE with ModelFinder for the maximum-likelihood
    tree plus ultrafast bootstrap; visualise with iTOL or FigTree.
    For selection scans, branch off into PAML's `codeml`; for
    divergence dating, into BEAST with a relaxed clock.
  ],
) <fig:ch20-workflow>

The list of *pitfalls* is what separates a publishable analysis from
a misleading one. The recurring patterns:

- *Misaligned columns.* A handful of badly aligned positions can
  rotate a deep branch. Trim aggressively with trimAl in automated
  mode or Gblocks at a conservative threshold before tree-building.
  Visualise the alignment in AliView and sanity-check the regions
  flagged by the trimmer.
- *Long-branch attraction.* @sec:trees covered the algorithmic side.
  The empirical fix is to break long branches by adding intermediate
  taxa — if you have a tree where two distant taxa are grouping,
  find species that sit on the branches between them and add them.
- *Saturation.* At very high divergence, multiple substitutions
  per site cause the JC correction to blow up and the observed
  distance to plateau. The model-selection step (ModelFinder) will
  often refuse to fit GTR+G+I to saturated data; the practical
  remedy is to drop the saturated positions (typically third codon
  positions) or to switch to amino-acid-level analysis.
- *Missing data and gaps.* Indels are handled by Gappy columns being
  trimmed or by gap-aware likelihood. Reporting gap-handling
  explicitly is part of reproducible phylogenetics.
- *Wrong outgroup.* Rooting an unrooted ML tree requires an
  outgroup, and picking an outgroup that is itself misplaced inverts
  ancestor relationships. A common safeguard is to root by a
  relaxed-clock analysis instead, where the clock determines the
  root position.

#figure(
  image("../figures/ch20/f4-pitfalls-checklist.svg", width: 95%),
  caption: [
    Five-question design-review checklist for a phylogenetics
    project. Each question targets a documented failure mode; running
    them once before the experiment and once before trusting another
    paper's tree avoids most of the field's recurring mistakes.
  ],
) <fig:pitfalls>

The five questions cover most of the failure landscape: is the
alignment trimmed sensibly, is the substitution model chosen by a
formal criterion, has long-branch attraction been ruled out, are the
branch-support values bootstrap or posterior, and is the root
defensible. A published tree that survives all five questions is
likely to replicate; a published tree that survives none is not.


== Summary <sec:ch20-summary>

- Exact multiple sequence alignment is NP-hard; every working tool
  is a heuristic. Progressive alignment along a guide tree
  (CLUSTAL Omega, MUSCLE, MAFFT) is the dominant paradigm.
  Iterative refinement and consistency-based scoring (T-Coffee,
  ProbCons) trim the "once a gap, always a gap" lock-in. MAFFT
  L-INS-i is the modern accuracy default; FFT-NS-2 is the speed
  default; Clustal Omega scales to hundreds of thousands of
  sequences.
- Phylogenetic trees come from three philosophies. Distance methods
  (UPGMA, neighbour-joining) collapse the alignment into a distance
  matrix and build a tree from it. Parsimony minimises substitutions
  on the tree but is statistically inconsistent under long-branch
  attraction. Maximum likelihood and Bayesian inference are the
  modern defaults; RAxML, IQ-TREE, and MrBayes are the standard
  tools.
- Felsenstein's pruning algorithm is belief propagation on a
  tree-structured graphical model, computing the marginal likelihood
  at the root by a single bottom-up sweep. Felsenstein invented it
  for phylogenetics in 1981, seven years before Pearl gave the
  general graphical-model treatment.
- The molecular clock is a Poisson process. Strict clocks fail in
  practice; relaxed-clock models in BEAST allow rates to vary along
  the tree and combine with fossil calibration points to give dated
  trees with credible intervals. Pathogen-outbreak dating
  (SARS-CoV-2, influenza, HIV) is the modern showpiece of the
  technique.
- dN/dS quantifies codon-level selection. $omega < 1$ is purifying
  selection (most genes most of the time), $omega = 1$ is neutral,
  $omega > 1$ is positive selection. PAML's `codeml` implements the
  Yang--Nielsen codon model; the likelihood-ratio tests M1 vs M2 and
  M7 vs M8 are the standard for detecting site-specific positive
  selection. Influenza HA, MHC, and sperm-egg recognition proteins
  are the canonical examples.
- Comparative genomics works at whole-genome scale. Synteny tracks
  preserved gene order across species; OrthoFinder reconciles gene
  trees with species trees to call orthologs and paralogs; PhyloP
  and PhastCons score per-base evolutionary conservation across
  vertebrates. Conservation tracks are the foundation of
  phylogenetic footprinting and of clinical variant prioritisation.
- Non-coding RNA homology needs structure-aware alignment.
  Stochastic context-free grammars (SCFGs) generalise profile HMMs
  to the nested base-pairing structure that defines RNA function.
  Infernal and Rfam apply SCFGs at scale and remain the standard
  for ncRNA family detection; they routinely outperform sequence-only
  methods at distant homology because base-pair covariation carries
  more information than per-column conservation.
- The EE-framing distillation. MSA is constrained dynamic
  programming on a $D$-dimensional tensor. Phylogenetic inference is
  graphical-model inference (Felsenstein 1981 = Pearl 1988 for
  trees). The molecular clock is a Poisson process. dN/dS is a
  likelihood-ratio test. RNA folding is constrained context-free
  parsing. Most of evolutionary genomics is electrical-engineering
  formalism applied to biological substrates.


== Exercises <sec:ch20-exercises>

#strong[1.] #emph[Progressive alignment.] Pick five protein sequences
from the same family — five vertebrate myoglobins or five bacterial
50S ribosomal L1 proteins. Run MAFFT L-INS-i, Clustal Omega, and
MUSCLE on the same input. Quantify the differences with the column
score (number of columns identical between two alignments divided by
total columns). Identify one functionally important residue
(e.g., the proximal histidine of myoglobin) and verify all three
aligners place it in the same column. Where they disagree,
investigate the local context.

#strong[2.] #emph[Distance correction.] Compute the JC-corrected
pairwise distance between two sequences with $30%$ raw divergence
and the K2P-corrected distance assuming a transition/transversion
ratio of $2$. Compute the difference between the two as a function of
the underlying time. Where does the correction matter most? Where
does it matter least?

#strong[3.] #emph[Felsenstein's pruning by hand.] Take a four-leaf
tree with branch lengths $0.1, 0.2, 0.3, 0.4$ in substitutions per
site, with all four leaves observed as base $A$. Under Jukes--Cantor
($mu = 1$), compute the likelihood at the root by hand. Repeat for
the same tree with leaves $A, A, G, G$. Use the closed-form
$P(x | y, t) = 1/4 + 3/4 exp(-4 mu t / 3)$ for the diagonal and
$P(x | y, t) = 1/4 - 1/4 exp(-4 mu t / 3)$ for the off-diagonal.

#strong[4.] #emph[Long-branch attraction.] Simulate a four-taxon
dataset under the true topology $((A, B), (C, D))$ with branch
lengths $A: 1.5$, $B: 0.05$, $C: 1.5$, $D: 0.05$ — two long branches
($A$, $C$) and two short branches ($B$, $D$). Generate $1000$
alignment columns under Jukes--Cantor. Build a parsimony tree and a
maximum-likelihood tree under JC. Which method recovers the true
topology? At what point — as you shorten the long branches — does
parsimony stop failing?

#strong[5.] #emph[dN/dS from scratch.] Take a pair of codon-aligned
orthologous coding sequences (the chimpanzee and human FOXP2
sequences are a common teaching example). Compute the Nei--Gojobori
$d_N$ and $d_S$ directly from the alignment by counting synonymous
and non-synonymous sites and substitutions. Compare your $omega$
estimate with the value from PAML's `codeml`. Where do the two
methods disagree, and why?

#strong[6.] #emph[Synteny inspection.] Visit the UCSC Genome Browser
and open the synteny track for human chromosome 17 against mouse.
Identify one clear inversion and one clear translocation. Estimate
the size of each rearrangement in base pairs. Are there genes
spanning the breakpoints? What is the relevance for human disease if
any of those genes is broken by the inversion in mouse?

#strong[7.] #emph[Infernal in anger.] Take a $10$-kb intergenic
region from any sequenced bacterial genome and run `cmsearch` against
Rfam with an E-value threshold of $10^(-3)$. Tabulate the hits.
For one tRNA hit, extract the aligned sequence and the
secondary-structure annotation; verify the predicted cloverleaf with
`RNAfold`. Do the two predictions agree?

#strong[8.] #emph[(Open-ended.)] Pick a recently published
phylogenomic paper (any genome-scale tree paper from 2023--2025) and
run the five-question design-review checklist (@fig:pitfalls)
against it. Identify the weakest of the five links. Propose one
analysis that would change your confidence in the headline tree, and
predict whether the tree would survive. Cite the paper.


== Further Reading <sec:ch20-further-reading>

- *Felsenstein, J.* (1981). "Evolutionary trees from DNA sequences:
  a maximum likelihood approach." *Journal of Molecular Evolution*
  17: 368--376. The paper that founded modern phylogenetics. Read
  it for the pruning algorithm and the pulley principle.
- *Saitou, N., and Nei, M.* (1987). "The neighbor-joining method: a
  new method for reconstructing phylogenetic trees." *Molecular
  Biology and Evolution* 4: 406--425. The distance-method workhorse
  that still ships in every phylogenetics toolkit.
- *Edgar, R. C.* (2004). "MUSCLE: multiple sequence alignment with
  high accuracy and high throughput." *Nucleic Acids Research* 32:
  1792--1797. Iterative refinement formalised; pair with the MAFFT
  references for the modern aligner landscape.
- *Katoh, K., and Standley, D. M.* (2013). "MAFFT multiple sequence
  alignment software version 7: improvements in performance and
  usability." *Molecular Biology and Evolution* 30: 772--780.
  The current MAFFT reference; the FFT-based distance idea is in
  the 2002 paper, but this one covers the production tool.
- *Yang, Z., and Nielsen, R.* (2000). "Estimating synonymous and
  nonsynonymous substitution rates under realistic evolutionary
  models." *Molecular Biology and Evolution* 17: 32--43. The codon
  substitution model behind PAML's `codeml`.
- *Nguyen, L.-T., Schmidt, H. A., von Haeseler, A., and Minh, B. Q.*
  (2015). "IQ-TREE: a fast and effective stochastic algorithm for
  estimating maximum-likelihood phylogenies." *Molecular Biology
  and Evolution* 32: 268--274. The modern default for medium-scale
  ML tree inference.
- *Hug, L. A., et al.* (2016). "A new view of the tree of life."
  *Nature Microbiology* 1: 16048. The comparative-genomics paper
  that reorganised the deep tree, placing archaea as sister to
  eukaryotes.
- *Eddy, S. R.* (2013). "Computational analysis of conserved RNA
  secondary structure in transcriptomes and genomes." *Annual
  Review of Biophysics* 43: 433--456. The Infernal author's own
  review of the SCFG/CM machinery. Pair with the Kalvari et al.
  Rfam reference for the database side.
- *Lorenz, R., et al.* (2011). "ViennaRNA Package 2.0." *Algorithms
  for Molecular Biology* 6: 26. The reference for `RNAfold` and the
  rest of the ViennaRNA suite.
