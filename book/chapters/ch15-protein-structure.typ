#import "../theme/book-theme.typ": *

= #idx("protein structure")Protein Structure Prediction in the #idx("AlphaFold")AlphaFold Era <ch:protein-structure>

#matters[
  For fifty years, predicting a protein's three-dimensional structure
  from its amino-acid sequence was one of the canonical unsolved
  problems in computational biology. The fact that it is now largely
  solved — for single-chain folds, on the timescale of minutes per
  protein, by a freely available neural network — is the most concrete
  triumph #idx("machine learning")machine learning has delivered to the life sciences. Almost
  every downstream biology problem that involves a protein, which is
  to say most of them, has been quietly re-grounded on top of it.
  Drug discovery, enzyme engineering, functional annotation, evolutionary
  comparison, the entire field of protein design — none of these look
  the same after 2021 as they did before. This chapter is a tour of
  what changed, how it was made to change, and what is still hard.
]

A protein is not, in any useful sense, a sequence. It is a chain of
amino acids that has folded into a compact three-dimensional shape,
and the shape is what determines the function. The chain is a
combinatorial object; the fold is a geometric one. The chapter you
have just left behind in the table of contents — population genetics
and genome-wide association — treated proteins implicitly, as
side-effects of variant calls in coding regions. That works only as
long as the protein structure is someone else's problem. Open a Cancer
Genome Atlas record, scroll through twenty missense variants, and the
clinical question is no longer "what base changed" but "did the fold
survive." Until 2020 the second question had no general answer for
most proteins. It now does.

This chapter walks the arc from Christian Anfinsen's refolding
experiment in 1961, which established that the answer exists in
principle, to AlphaFold2 in 2020, which produced the answer in
practice. Along the way it covers the half-century of failed attempts,
the statistical-physics insight that opened the door in 2011, the
transformer-based architecture that walked through it in 2020, and
the inverse problem — protein design — that became tractable
immediately afterwards. The mathematics is dense in places, and the
chapter occasionally lingers on a derivation longer than the lecture
would. That is deliberate. The #idx("Evoformer")Evoformer and the #idx("structure module")structure module
are the most carefully reasoned pieces of architectural design in
applied machine learning; reading them quickly leaves their
elegance invisible.

If you want a single take-home for the chapter: the breakthrough was
not a new energy function or a faster sampler. It was a way of reading
information that had been hiding in plain sight for thirty years —
the pattern of co-mutation across homologous sequences — and folding
that information through a deep neural network that knew, by
construction, how rigid bodies move in three-dimensional space.


== Why Structure, and What It Is <sec:why-structure>

A typical protein is a polymer of fifty to two thousand amino acids
drawn from an alphabet of twenty standard residues. The chain is
synthesised by a #idx("ribosome")ribosome reading messenger #idx("RNA")RNA in groups of three
nucleotides (Chapter 1); the chemistry of the side chains decorates
the backbone with the full diversity of organic functional groups —
carboxylates, amines, aromatics, sulphurs, hydroxyls. Almost the
moment the chain leaves the ribosome it begins to fold, driven by the
hydrophobic effect to bury its greasy residues in a compact interior
and expose its polar residues to water. Within milliseconds to
seconds, the chain has settled into a near-unique three-dimensional
arrangement of atoms. That arrangement is the structure.

#figure(
  image("../../diagrams/lecture-15/01-structure-hierarchy.svg", width: 90%),
  caption: [
    The four structural levels. Primary is the linear amino-acid
    sequence; secondary is local hydrogen-bonded patterns (helices
    and sheets); tertiary is the global fold of a single chain;
    quaternary is multi-chain assembly. When people say "the
    structure" they almost always mean the tertiary fold.
  ],
) <fig:hierarchy>

The four hierarchical levels in @fig:hierarchy are conventional. The
#idx("primary structure")primary structure is the sequence itself, a #idx("STRING")string in a 20-letter
alphabet. The #idx("secondary structure")secondary structure is local backbone conformation —
the $alpha$-helix (a right-handed spiral with 3.6 residues per turn
and intra-helix hydrogen bonds) and the $beta$-strand (an extended
conformation that pairs with other strands to make a sheet). A
typical globular protein is thirty to fifty percent helix and
strand; the rest is loops and turns. The #idx("tertiary structure")tertiary structure is the
full three-dimensional arrangement of the chain — every backbone
torsion ($phi.alt$, $psi$, $omega$) and every side-chain torsion
($chi_1$, $chi_2$, $...$). Given those torsions and the standard
bond geometry, every atom's $(x, y, z)$ coordinate is determined.
The #idx("quaternary structure")quaternary structure is the assembly of multiple chains —
haemoglobin is four chains, the ribosome is eighty chains plus RNA.

In bioinformatics practice, "solving a structure" almost always means
solving the tertiary fold. Secondary structure is a side product —
classical tools like PSIPRED predicted it at ~85 % accuracy in the
1990s, and modern transformers do it essentially perfectly. Quaternary
structure is harder, and it is what AlphaFold-Multimer and AlphaFold-3
were specifically built for; we return to it in @sec:ch15-af3.

Three representations of a tertiary structure recur throughout this
chapter. *Cartesian coordinates*: an $(x, y, z)$ for every atom,
typically eight atoms per residue, written in PDB or mmCIF format.
*Internal coordinates*: backbone torsions plus side-chain torsions
plus standard bond geometry; more compact and the representation
AlphaFold actually outputs before converting back to Cartesian. And
the *contact map*: an $N times N$ binary or continuous matrix where
entry $(i, j)$ encodes the distance between residues $i$ and $j$.
Contact maps throw away the rotation and #idx("translation")translation of the whole
protein and keep only the relative geometry; this is what
coevolution-based methods predict directly, and what AlphaFold's
internal pair representation eventually becomes.

=== The Anfinsen Experiment

Christian Anfinsen, working at the National Institutes of Health in
the late 1950s, asked a deceptively simple question. If you take a
folded protein, unfold it completely by chemical denaturation, and
then remove the denaturant under conditions where the protein could
in principle refold — does it? And if so, does it refold to the same
structure it started with?

The model system was ribonuclease A, a small (124-residue) digestive
enzyme. Anfinsen denatured it with 8 M urea, reduced its four
disulphide bonds with $beta$-mercaptoethanol, then dialysed both
chemicals out under air. The enzyme refolded spontaneously, all four
disulphide bonds re-formed in their correct pairings, and the
catalytic activity returned to ninety-some percent of native. There
was no chaperone in the tube. There was no template. There was nothing
but the sequence, the solvent, and time.

The implication, set out fully in Anfinsen's 1973 Nobel lecture, is
what the field now calls the *thermodynamic hypothesis*: the native
fold of a protein is the global minimum of a free-energy landscape
over conformational space, and that minimum is uniquely determined by
the sequence. Computationally, this is the most important sentence
in the chapter. It says that the structure-prediction problem is
well-defined — given a sequence, there is in principle a correct
answer — and it transforms the problem from "what does this molecule
do" into "find this energy minimum."

#note[
  Anfinsen's conclusion has been refined over six decades but never
  overturned. Some proteins fold with help from molecular chaperones,
  some pass through transient kinetic traps, and a small minority
  ("intrinsically disordered proteins") never fold at all in
  isolation. But for the canonical case of a globular soluble
  protein, sequence determines structure. The fifty years from
  Anfinsen's paper to AlphaFold2's CASP14 result are essentially the
  history of trying to compute that determination.
]

The trouble, of course, is that the energy landscape is staggering.
*#idx("Levinthal")Levinthal's paradox* (Cyrus Levinthal, 1969) makes the staggering
concrete: a chain of 100 residues, each with three accessible
backbone conformations, has $3^100 approx 10^48$ conformations.
Even if the chain could sample one conformation per nanosecond, an
exhaustive search would take longer than the age of the universe.
Real proteins fold in milliseconds. The folding landscape is not
sampled randomly — it is *funnelled* toward the native state by the
hydrophobic effect and by a hierarchy of local stabilising
interactions that get committed early. Finding the algorithm that
captures this funnel, without exhaustively sampling, is the problem
the next fifty years tried to solve.

=== Domains and the Fold-Space Universe

Most natural proteins are mosaic objects. They are built from
*domains* — semi-independent structural units, fifty to two hundred
residues long, that fold autonomously and recur across many proteins
as mix-and-match modules. A kinase has a catalytic domain plus a
regulatory domain plus a targeting domain, and each is found, in
slightly different form, in dozens of other proteins.

The number of distinct *folds* — recurring three-dimensional
architectures of secondary-structure elements — is finite. SCOP and
CATH, the two classical fold classifications, list roughly one to
two thousand folds in nature. New folds were rare and exciting events
in the structural-biology literature before 2020. After AlphaFold
was applied to metagenomic sequences, a modest number of genuinely
new folds were discovered, but the headline result is that the fold
universe really is small — backbone physics and residue chemistry
constrain what can stably fold.

This matters for prediction because it sets the upper bound on what
the model needs to learn. If there were $10^9$ distinct folds in
nature, a model trained on the 200,000 structures in the PDB would
have no chance. With $approx 10^3$ folds and millions of protein
sequences in homologous families, the empirical statistics are
plentiful. AlphaFold's training set is finite; the universe it is
predicting over is also finite.


== The Classical Era and Why It Plateaued <sec:classical>

Before 2018, protein structure prediction had three main strategies,
each with a sharply different failure mode.

*Homology modelling* works when the query protein has a sequence-similar
homologue with a known experimental structure. The recipe is mechanical:
#idx("BLAST")BLAST the query against PDB sequences, find a "template" with high
identity, copy the template's backbone coordinates for aligned
residues, place the query's side chains using a rotamer library
(SCWRL, Rosetta), then energy-minimise. Tools: *Modeller* (Sali, 1993),
*SWISS-MODEL* (Schwede et al., 2003). Accuracy degrades with
sequence-identity-to-template: above 50 % identity the result is
typically usable (backbone RMSD $approx 1$ Å); between 30 and 50 % is
the "twilight zone" where the global fold is right but local details
are wrong; below 30 % the model is essentially unusable. For proteins
with no close PDB homologue — orphan proteins, or anything from
under-sequenced clades — homology modelling produces nothing.

*Threading* (Bowie, Lüthy, and Eisenberg, 1991; *I-TASSER* under Yang
Zhang from 2008) tries to break the homology requirement. For each
known fold in the PDB, score how well the query's sequence "fits"
that fold using sequence–structure compatibility metrics. Pick the
best fit, build the model. Threading recovers the case where the
query shares a fold with a known protein but has diverged at the
sequence level — it found cases homology modelling missed, but the
fits were always borderline.

*Ab initio* methods predict structure without any template. Two
flavours dominated. *Fragment assembly*, embodied by David Baker's
Rosetta (1997+), represented the protein chain as a sequence of short
(9-residue) structural fragments drawn from the PDB and assembled
them by Monte Carlo sampling under a learned energy function. *Molecular
dynamics* simulated the folding #idx("trajectory")trajectory atom by atom with a
classical force field; in 1998, Peter Kollman's group folded the
villin headpiece (35 residues) in months of supercomputer time.
Folding\@home and Anton scaled MD to milliseconds of simulated
trajectory by 2010, but anything beyond very small proteins remained
out of reach.

#figure(
  image("../../diagrams/lecture-15/02-casp-progression.svg", width: 90%),
  caption: [
    #idx("CASP")CASP free-modelling accuracy (median best GDT_TS on the hardest
    target category) from CASP1 in 1994 through CASP16 in 2024.
    Decades of incremental progress, then two step changes —
    CASP13 (AlphaFold-1) and CASP14 (#idx("AlphaFold-2")AlphaFold-2) — and a high
    plateau thereafter.
  ],
) <fig:casp>

The community benchmark for structure prediction is *CASP* — the
Critical Assessment of Structure Prediction, run biennially since
1994 by John Moult and colleagues. The format is straightforward:
participants are given the sequences of soon-to-be-released
experimental structures; they submit predictions; after the
experimental structures are released, predictions are scored against
them. CASP tracks several categories — template-based modelling,
free modelling, multimers, ligand #idx("docking")docking, contact prediction — and
its primary metric is *GDT_TS* (Global Distance Test, Total Score),
a percentage from zero to one hundred where ~50 means "recognisable
fold," ~70 means "good prediction," and ~90 corresponds roughly to
experimental quality.

@fig:casp tells the story compactly. From CASP1 (1994) through CASP12
(2016), progress on free-modelling targets was steady but slow.
Rosetta and I-TASSER traded the lead. Median best scores climbed from
the twenties to the forties over twenty years. At CASP13 (2018),
DeepMind's AlphaFold-1 entered, used a deep neural network to predict
inter-residue distance distributions, fed those into a gradient
descent folder, and jumped the free-modelling score by roughly twenty
GDT_TS points in one edition — an unprecedented leap. At CASP14
(2020), AlphaFold-2 reached a median best of ~90 GDT_TS on free
modelling: indistinguishable, on most targets, from a real
experimental structure. The room at the CASP14 conference, by all
accounts, fell quiet.

#note[
  The plateau in @fig:casp before 2018 is not a story about
  insufficient compute. Rosetta-era methods ran on supercomputers for
  weeks per protein and still topped out at ~40 GDT_TS on hard
  targets. The plateau is about *information* — the methods were not
  exploiting the rich constraint signal that lives in homologous
  sequence families. Two changes had to happen before that signal
  could be read: a statistical-physics method (#idx("DCA")DCA, 2011) that knew
  how to extract direct from indirect coupling, and a neural
  architecture (the Evoformer, 2020) that could fold the extracted
  signal into atomic-resolution geometry. The decade from 2011 to
  2020 was the story of those two pieces converging.
]

Two structural reasons explain the plateau. First, template
dependence: in the absence of a PDB homologue, classical methods had
no anchor and produced essentially random structures. Second,
imperfect energy functions: Rosetta's energy is an approximation to
the real physics, and its minimum does not always correspond to the
native state. Improving the energy function by hand-tuning terms gave
diminishing returns through CASP6–CASP12. The break came not from
improving the energy function but from changing the input. The
information was always in the sequences; the trick was to read many
sequences at once.


== #idx("coevolution")Coevolution: Reading Contact From the #idx("MSA")MSA <sec:coevolution>

Given a protein family — say, all known homologues of a kinase
catalytic domain — align the sequences into a *multiple sequence
alignment* (MSA). Each column is one residue position; each row is
one species' version of the protein. A globin family alignment from a
typical sequence database is several thousand rows by ~150 columns.

The empirical observation that drives the whole field is this. Some
pairs of columns *vary together*. If column $i$ mutates from leucine
to valine in a particular species, column $j$ — half the protein away
in the linear sequence — frequently also mutates, in the same row.
Across thousands of homologous sequences from thousands of species,
this co-mutation pattern is robust and reproducible.

The physical explanation is the punchline. Residues $i$ and $j$ are
in three-dimensional *contact* in the folded structure. A mutation at
$i$ that changes its size or polarity destabilises the contact; for
the fold to survive, $j$ must co-mutate to a compatible partner.
Evolution propagates many such double mutations across the family.
The MSA records the cumulative imprint. Reading which column pairs
co-mutate recovers the contact map.

#figure(
  image("../../diagrams/lecture-15/03-coevolution.svg", width: 90%),
  caption: [
    The mechanism. Pairs of residues in 3D contact must co-mutate to
    maintain stability; the MSA accumulates that constraint across
    thousands of homologues. Reading column-pair co-mutation
    recovers contacts.
  ],
) <fig:coevolution>

=== Mutual Information and Its Failure

The naïve estimator is mutual information. For each column pair
$(i, j)$, treat each column as a categorical distribution over the 20
amino acids and compute

$ "MI"(i, j) = sum_(a, b) P(a, b) log (P(a, b)) / (P(a) P(b)) $

where $P(a)$ is the marginal probability of amino acid $a$ at position
$i$, and $P(a, b)$ is the joint probability over columns $i$ and $j$.
High MI means the two columns carry information about each other.

The problem is *indirect correlations*. If columns $i$–$j$ and $j$–$k$
are both in contact, then columns $i$–$k$ will also show elevated MI
even if $i$ and $k$ are not in direct contact — the chain of
correlations propagates. Worse, if the family is shallow, sampling
noise alone produces enough joint-distribution variation to throw up
spurious MI signals. Through the 1990s and 2000s, MI-based contact
predictors were known to find roughly forty percent of true contacts
in the top-$L$ predictions — better than chance, but nowhere near the
accuracy needed to seed a folding calculation.

=== #idx("direct coupling analysis")Direct Coupling Analysis

The insight that unlocked the field came from statistical physics. If
you model the joint distribution over sequences as a Potts model —

$ P(sigma) prop exp ( sum_i h_i (sigma_i) + sum_(i < j) J_(i j) (sigma_i, sigma_j) ) $

where $sigma = (sigma_1, dots, sigma_L)$ is a sequence, the $h_i$ are
per-position biases capturing conservation, and the $J_(i j)$ are
pairwise couplings — then the maximum-entropy distribution
consistent with the observed pairwise marginals has *direct*
couplings $J_(i j)$ whose magnitude reveals which pairs are
genuinely in contact. Indirect correlations are "explained away" by
the model's chain of intermediate couplings.

Two algorithms compete on the inverse problem of estimating $J_(i j)$
from an MSA. *Mean-field DCA* (mfDCA; Faruck Morcos and collaborators
in Martin Weigt's group, 2011) gives a closed-form approximation that
runs in seconds per family. *Pseudo-likelihood DCA* (plmDCA; Magnus
Ekeberg and colleagues, 2013) is more accurate, runs in minutes, and
was the standard until #idx("deep learning")deep learning displaced it. Both produce, for
each column pair $(i, j)$, a 20 × 20 coupling matrix; the *contact
score* is the Frobenius norm

$ C_(i j) = || J_(i j) ||_F = sqrt(sum_(a, b) J_(i j) (a, b)^2). $

Large $C_(i j)$ means predicted contact. On proteins with deep MSAs,
the top-$L$ contacts from plmDCA are roughly 80 % correct — enough,
combined with fragment assembly in Rosetta, to fold the protein
approximately. *EVfold* (Debora Marks and colleagues, 2011) was the
first widely used DCA-based contact predictor; it was rapidly joined
by *PSICOV*, *GREMLIN*, and *CCMpred*. By CASP12 (2016) DCA + Rosetta
was the state of the art on hard targets, beating template-free
predictors that didn't use coevolution by a comfortable margin.

#figure(
  image("../../diagrams/lecture-15/04-dca-contacts.svg", width: 90%),
  caption: [
    DCA recovers sparse true contacts by inverting covariance.
    Mutual information tangles direct and indirect correlations;
    the inverse-covariance estimator separates them. Left: raw MI
    matrix, noisy. Middle: DCA inverse-covariance score, sharp.
    Right: native contact map from the experimental structure.
  ],
) <fig:dca>

#note[
  Reading co-mutations from an MSA is precisely the *graphical-model
  inverse-covariance-estimation* problem. The MSA is a matrix of
  categorical samples from a joint distribution; the true contact
  graph is the inverse-covariance sparsity pattern. Naïve mutual
  information is the analogue of raw covariance — it sums direct
  and indirect dependencies indiscriminately. DCA is the analogue of
  the precision matrix — it isolates direct dependencies by
  conditioning on every other variable. The Gaussian version of the
  same problem is Friedman, Hastie, and Tibshirani's graphical lasso
  (2008); plmDCA is the discrete-variable, Potts-model version. The
  DCA papers cite statistical-mechanics literature rather than
  machine-learning statistics, but the mathematical content is
  identical.
]

The catch is that DCA's accuracy depends critically on MSA depth.
@fig:msa-depth makes the dependence explicit. Below fifty effective
sequences, the estimator is too noisy to be useful — sampling
variance dominates the signal. Above a thousand effective sequences,
top-$L$ contact precision climbs to ~80 % and saturates around
50,000 sequences at ~85 %. Orphan proteins — those with no
detectable homologues — are invisible to DCA, and roughly 20 % of
human proteins fall in or near that category. The other curve in
@fig:msa-depth shows what changed with deep learning: AlphaFold-2's
contact precision on the same benchmark *starts* at ~0.45 for shallow
MSAs and saturates near 0.95 by a thousand sequences. The improvement
on deep MSAs is real but modest; the improvement on shallow ones is
the qualitative breakthrough.

#figure(
  image("../figures/ch15/f1-msa-depth-vs-accuracy.svg", width: 95%),
  caption: [
    Top-$L$ contact precision as a function of effective MSA depth,
    for classical DCA versus AlphaFold-2 on a CAMEO-style benchmark.
    DCA needs ~1,000 sequences for reasonable accuracy; AlphaFold-2
    extracts useful signal from shallow alignments because it has
    learned strong structural priors at training time.
  ],
) <fig:msa-depth>

=== From Contacts to Coordinates

A contact map by itself is not a structure. To go from predicted
contacts to a 3D fold, the pre-AlphaFold pipeline treated each
predicted contact $(i, j)$ as a soft distance constraint $d_(i j) < 8$
Å and ran fragment assembly (Rosetta) or distance-geometry embedding
(CNS) subject to those constraints. Pick the lowest-energy resulting
model. This was *Marks–Onuchic–Hopf 2011* on the toy case, *EVfold*
generalised, and Yang Zhang's I-TASSER folded into a workflow that
dominated CASP10–CASP12.

AlphaFold-1 (CASP13, 2018) replaced the Rosetta step with a deep
convolutional network that predicted distance *distributions*
($p(d | i, j)$) rather than binary contacts, and a gradient-descent
folder that minimised an aggregate distance loss. AlphaFold-2 (CASP14,
2020) abandoned the contact-then-fold two-step pipeline entirely.
It is end-to-end: sequence and MSA in, atomic coordinates out, with
the contact map appearing only as an emergent internal representation.
That is the architecture we walk through next.


== AlphaFold-2: Architecture <sec:af2>

AlphaFold-2 was published by John Jumper and 18 colleagues in
_Nature_ in July 2021, after the CASP14 result had already become the
biggest story in the field. The paper describes a network with about
93 million trainable parameters that takes a protein sequence and its
MSA as input and produces atomic coordinates plus per-residue
confidence as output. The architecture has three major components —
the Evoformer, the structure module, and a recycling loop — plus
several confidence heads and a template-embedding side channel.
@fig:af2 lays out the overall flow.

#figure(
  image("../../diagrams/lecture-15/05-af2-overview.svg", width: 95%),
  caption: [
    AlphaFold-2 end to end. Sequence and MSA in; coordinates and
    confidence out. The Evoformer refines two representations
    (MSA and pair) jointly; the structure module decodes them into
    3D coordinates; recycling feeds the structure back as input
    for another pass.
  ],
) <fig:af2>

=== The Two Representations

Two data objects are maintained throughout the network. The *MSA
representation* $bold(M) in bb(R)^(s times r times c)$ is a tensor of
shape (sequences $times$ residues $times$ channels), initially $c = 256$.
Each row is one MSA sequence, each column is one residue position, and
each cell holds a learned $c$-dimensional embedding. The *pair
representation* $bold(Z) in bb(R)^(r times r times c)$ is a residue
$times$ residue $times$ channels tensor — a learned, richly-featured
$N times N$ matrix where entry $(i, j)$ encodes the relationship
between residues $i$ and $j$. Both representations are updated
iteratively by the Evoformer. The structure module reads from $bold(Z)$
(and a "single" representation distilled from $bold(M)$) to produce
coordinates.

You can think of $bold(Z)$ as a generalised contact map: not a binary
matrix saying "in contact / not in contact," but a $c$-channel
embedding that encodes "how far apart, in what orientation, with what
backbone conformation." The whole point of the Evoformer is to refine
this embedding until the structure module can decode it into
coordinates.

=== The Evoformer Block

The Evoformer is 48 repeated blocks. Each block updates both $bold(M)$
and $bold(Z)$ through six sub-operations, run in a specific order:

1. *Row-wise MSA #idx("attention")attention with pair bias.* Within each MSA row (one
   sequence), every residue position attends to every other position.
   The attention map is biased by features pulled from the pair
   representation $bold(Z)$, so contact information already in
   $bold(Z)$ shapes which residues talk to which.
2. *Column-wise MSA attention.* Within each column (one residue
   position across all species), every sequence attends to every
   other. This is how each residue learns about its homologous
   versions.
3. *MSA transition* (a position-wise feedforward layer).
4. *Outer-product mean.* The MSA is summarised into an outer-product-
   mean statistic computed across species and fed into the pair
   representation. This is the channel through which MSA-derived
   information enters $bold(Z)$.
5. *Triangle multiplicative update.* The pair representation updates
   itself using triangle geometry: if $(i, k)$ and $(j, k)$ are both
   "contact-like" in $bold(Z)$, the update pushes $(i, j)$ toward
   contact-like as well. This imposes the triangle-inequality prior
   that any consistent set of pairwise distances must satisfy.
6. *Triangle attention.* An attention-based analogue of the
   multiplicative update: each pair $(i, j)$ attends over all pairs
   sharing a third residue $k$.

@fig:evoformer shows one block as a directed graph of these six
operations; the block stacks 48 deep, with separate weights at each
layer.

#figure(
  image("../../diagrams/lecture-15/06-evoformer-block.svg", width: 95%),
  caption: [
    One Evoformer block as a DAG. Six operations refine MSA and pair
    jointly; the block stacks 48 deep with separate weights per
    layer. The dotted "outer-product mean" arrow is the channel
    through which MSA-derived evolutionary signal enters the pair
    representation.
  ],
) <fig:evoformer>

The most striking feature of the design is *axial attention* on the
MSA. Full 2-D attention over an MSA of $s$ sequences and $r$ residues
would cost $O(s^2 r^2)$ per layer, which at $s = r = 512$ comes to
$approx 7 times 10^10$ operations per layer per attention head — at
48 blocks, the total is intractable. Row-wise attention costs
$O(s r^2)$; column-wise attention costs $O(r s^2)$. The combination
recovers most of the expressive power of full attention at a fraction
of the cost (@fig:axial-attention). Axial attention was popularised
by Jonathan Ho and colleagues in image transformers in 2019 and
adapted by Jumper's group for the MSA case.

#figure(
  image("../figures/ch15/f2-axial-attention.svg", width: 95%),
  caption: [
    Full 2-D attention on the MSA versus axial attention. The two
    alternating sweeps — row attention across residues, column
    attention across species — produce roughly the same effective
    receptive field as the full operation at $approx 250 times$
    lower cost for the canonical $s = r = 512$ case.
  ],
) <fig:axial-attention>

#note[
  The Evoformer is a two-axis #idx("transformer")transformer. One axis is residues
  (tokens along the chain); the other is species (tokens down the
  MSA). The triangle updates in the pair representation are
  effectively a graph-neural-network layer where edges update
  themselves based on two-step neighbours, with the triangle
  inequality as a hard-wired geometric prior. Reading the paper, the
  most surprising design move is just how *much* of the standard
  transformer toolkit Jumper's group threw at the problem — six
  distinct sub-operations per block, 48 blocks, three flavours of
  attention — and how each piece earns its keep in the ablations.
  Removing the triangle updates costs about ten GDT_TS points.
  Reducing depth to 24 blocks costs five.
]

=== Why It Works

An intuitive reading of the Evoformer: its job is to construct a
high-quality pair representation $bold(Z)$ whose internal embedding
encodes not just contact-or-not but distance and relative
orientation — enough geometric information for a downstream 3D
embedder to decode. The MSA is the input signal, a noisy compilation
of evolutionary co-mutation patterns. The pair representation is the
output signal, a clean per-residue-pair geometric description. The
48 blocks iteratively distil the gossip in the MSA into the structured
answer in $bold(Z)$.

Why 48 blocks? Empirical. Fewer is worse on benchmark accuracy; more
gives diminishing returns at higher memory cost. The same applies to
many of the architectural choices in AlphaFold — they survived because
they worked, not because anyone has a clean theoretical argument for
the exact depth.


== AlphaFold-2: The Structure Module <sec:structure-module>

After 48 blocks of Evoformer, the network has a pair representation
$bold(Z) in bb(R)^(r times r times c)$ and a "single" representation
$bold(s) in bb(R)^(r times c)$ distilled from the MSA. The *structure
module* takes these and produces 3D coordinates. It is 8 iterations
of the same block, applied with shared weights.

At each iteration, every residue maintains a *backbone frame*
$bold(T)_i in "SE"(3)$ — a 3D position plus a 3D orientation, six
degrees of freedom per residue. The block has four parts:

1. *#idx("invariant point attention")Invariant point attention (IPA).* The novel attention mechanism
   at the heart of the module, described below.
2. *Frame update.* Each residue gets translated and rotated based on
   the IPA output.
3. *Compute per-residue points.* The $C_alpha$ position and the
   side-chain torsion angles are read off from the updated frame.
4. *Predict per-atom coordinates.* The standard amino-acid geometry
   library (idealised bond lengths and angles) is composed with the
   torsions to recover all-atom coordinates.

The structure is iteratively refined: at each iteration every residue
knows its current 3D position and the positions of all others (via
IPA); it then updates its own frame based on local geometric signals
and the pair representation.

=== Invariant Point Attention

The crucial property of IPA is that it is *equivariant* under global
rotations and translations of the whole protein. If you take the
input — sequence, MSA, pair representation, all the same — and apply
some rotation $R in "SO"(3)$ to the structure, IPA's output is the
input rotated by the same $R$. This respects the physics: a protein's
energy, contacts, and stability are invariant to its orientation in
the lab; the predictor should be too.

IPA achieves equivariance by construction. Each residue $i$ projects
its query to a set of *points in 3D space*, expressed in $i$'s own
local frame $bold(T)_i$. Each residue $j$ projects its key and value
to points in $j$'s local frame $bold(T)_j$. To compute the attention
weight between $i$ and $j$, the network re-expresses both sets of
points in a common reference (either's frame works, by equivariance)
and computes the squared distance between $i$'s query points and
$j$'s key points. The attention weight is a softmax over the negative
distance. Because the distance is computed between points expressed
in their respective local frames, and because the local frames
rotate with the protein, applying a global rotation $R$ rotates
every per-residue frame by $R$ — leaving all pairwise distances
unchanged. The attention map is invariant; the output, which is a
geometric quantity, is equivariant.

#figure(
  image("../../diagrams/lecture-15/07-structure-module.svg", width: 95%),
  caption: [
    The structure module folds 8 IPA iterations deep. Each step
    every residue moves based on $"SE"(3)$-equivariant attention to
    all others; rotation of the input rotates the output the same
    way. The bottom inset shows the key geometric quantity — points
    expressed in each residue's local frame, mapped to a common
    reference for the distance computation.
  ],
) <fig:structure-module>

#note[
  IPA is one of a small family of *$"SE"(3)$-equivariant networks*
  developed in parallel by Mario Geiger and Tess Smidt (e3nn, NequIP),
  Fabian Fuchs and Daniel Worrall (SE(3)-Transformers), and others.
  The shared principle: treat per-residue features as tensors with
  well-defined rotation behaviour (scalars, vectors, higher-order
  tensors), and construct operations that preserve the rotation
  class. Without built-in equivariance, the network would have to
  *learn* rotational invariance from data, which is wasteful and
  imperfect. The same equivariance principle shows up in particle
  physics (gauge invariance), in image processing (rotation-equivariant
  CNNs), and in molecular property prediction (DimeNet, GemNet).
]

=== The FAPE Loss

The structure module's training signal is the *Frame Aligned Point
Error* (FAPE), an unusual but effective loss function. For every pair
$(i, j)$ of residues, transform the predicted $C_alpha$ position of
residue $j$ into residue $i$'s predicted local frame, and the true
$C_alpha$ of $j$ into $i$'s true local frame, then take the Euclidean
distance between the two transformed points. Average over all pairs.

The benefit over a naïve coordinate-RMSD loss is that FAPE is
*frame-aware*: it does not penalise a globally-rotated prediction;
it penalises only *relative* geometric error between residues. It is
also bounded above (clamped at 10 Å in the AlphaFold paper), which
prevents a single very-wrong residue from dominating the gradient.
FAPE gives the structure module a smooth, geometrically meaningful
training signal.

=== Recycling and Confidence

After the structure module produces coordinates, AlphaFold-2 does
something a little surprising: it feeds the final pair representation
*and* the predicted structure back into the Evoformer as input, then
reruns the whole stack. Typically three to four such recycling
iterations are used in inference. Why? Because the Evoformer can be
too aggressive in early passes and reach a weird intermediate state.
Recycling lets the network correct itself, with each cycle benefiting
from the structure information produced previously. The recycling
loop is what lets a single network with $approx 93$M parameters
deliver atomic-accuracy structures without an explicit ensembling
step.

The other thing AlphaFold-2 produces is *confidence*. Three signals
matter:

- *#idx("pLDDT")pLDDT* (predicted local distance difference test). Per-residue,
  scaled 0–100. Operationally: pLDDT $> 90$ means "very confident
  this residue is positioned within $approx 1$ Å of correct"; $70$–$90$
  means "confident, backbone usable"; $50$–$70$ means "low confidence,
  fold-class only"; $< 50$ means "probably disordered or wrong."
- *#idx("PAE")PAE* (#idx("predicted aligned error")Predicted Aligned Error). Per-pair, in Ångstroms.
  $"PAE"_(i j)$ is the expected error in residue $j$'s position when
  the structure is aligned to residue $i$'s frame. Low PAE means
  "$j$'s position is well-defined relative to $i$"; high PAE means
  "I don't know where $j$ is relative to $i$."
- *pTM* (predicted TM-score). A single scalar, 0–1, that summarises
  overall predicted correctness.

#figure(
  image("../../diagrams/lecture-15/08-plddt-pae.svg", width: 95%),
  caption: [
    pLDDT colours the structure (per-residue local confidence); PAE
    shows the matrix (per-pair relative confidence). A two-domain
    protein typically has two dark blocks on the PAE diagonal and
    bright off-diagonal regions, indicating confident individual
    domains and uncertain inter-domain orientation.
  ],
) <fig:plddt-pae>

The combination of pLDDT and PAE is critical for interpreting
multi-domain proteins. Residues within a single domain typically have
low intra-domain PAE — the domain is confidently folded. Residues
across domains often have high inter-domain PAE — each domain is
fine in isolation, but their relative orientation is uncertain
(@fig:plddt-pae). A reader who looks only at pLDDT and sees uniformly
high values can be misled into trusting the full multi-domain
geometry when the PAE matrix is screaming that the inter-domain
hinge is floppy.

#figure(
  image("../figures/ch15/f3-confidence-cheatsheet.svg", width: 95%),
  caption: [
    A reading guide for AlphaFold confidence metrics. The pLDDT
    bands tell you which residues to trust at what level of
    detail; canonical PAE patterns distinguish single domains
    from multi-domain proteins, flexible linkers, and intrinsic
    disorder; the decision tree at the bottom maps confidence
    signatures to recommended downstream actions.
  ],
) <fig:confidence>

#warn[
  pLDDT is calibrated in the sense that "pLDDT = 80" averaged over
  many residues really does correspond to "$approx 1$ Å backbone
  error." But it is a *local* metric. A protein where every residue
  has pLDDT 90 can still have its two halves rotated relative to each
  other if the inter-domain PAE is high. Reading the confidence
  metrics together is non-negotiable. @fig:confidence summarises the
  decision logic the rest of this chapter assumes.
]

=== Training Regime

AlphaFold-2 was trained on:

- The PDB ($approx 170,000$ experimental structures at training time).
- UniRef90 plus BFD (hundreds of millions of sequences, used for MSA
  generation by HHblits and JackHMMER).
- *Self-distillation*: predict structures for UniRef sequences without
  experimental data, take high-confidence predictions, add them to
  the training set. This roughly doubled the effective training
  corpus and was crucial for shallow-MSA generalisation.

Training ran for about 11 days on 128 TPUv3 cores. Inference on a
single protein of typical length takes one to five minutes on a
single A100 GPU. The compute asymmetry between training and inference
is favourable: train once, predict billions.


== AlphaFold-3 and Its Successors <sec:ch15-af3>

AlphaFold-2's domain was single-chain protein structures. The
chapter's next three years filled in the remaining biology.

*AlphaFold-Multimer* (Richard Evans and colleagues, 2022) extended
AlphaFold-2 to multi-chain assemblies by simply concatenating the
input sequences with a separator token and retraining on PDB
biological assemblies. It works reasonably for small complexes (2–5
chains) but struggles with large assemblies and antibody–antigen
interactions.

*AlphaFold-3* (Joshua Abramson and 26 co-authors, _Nature_ 2024)
rewrote the architecture more radically. Three changes mattered.
First, the structure module is replaced by a *diffusion model*:
instead of iterative IPA frame updates, AF3 generates coordinates by
denoising from Gaussian noise, conditioned on a refined pair
representation (@fig:af3). This makes the architecture more
flexible — it handles variable-size outputs (with or without a
ligand), multimodal distributions (multiple conformations), and
heterogeneous input types in a single pass.

#figure(
  image("../../diagrams/lecture-15/12-af3-diffusion.svg", width: 95%),
  caption: [
    AlphaFold-3's diffusion structure module. Forward process: add
    Gaussian noise to true coordinates over $T$ timesteps until
    coordinates are pure noise. Learned reverse process: denoise
    iteratively, conditioned on the pair representation. Multi-
    component inputs (protein + #idx("DNA")DNA + ligand) denoise together in
    a single pass.
  ],
) <fig:af3>

Second, AF3 supports *multi-modal inputs*: a single inference call
can contain protein chains, nucleic-acid chains (DNA or RNA), small-
molecule ligands, and ions. The output places all components in their
predicted bound state — opening up protein–drug docking, RNA-binding
proteins with their bound RNA, and antibody–antigen complexes as
routine inference tasks. Third, the Evoformer is replaced by a
"Pairformer" that drops the explicit MSA axis once early
representations are formed, making inference much faster.

#note[
  Diffusion in AF3 is a textbook *denoising diffusion probabilistic
  model* (DDPM) applied to 3D coordinates. Forward process: add
  Gaussian noise to true coordinates over $T$ timesteps until they
  are pure noise. Learned reverse process: a neural network predicts
  the noise (or the clean signal) at each timestep, trained across
  many proteins and noise levels. At inference, start from Gaussian
  noise, denoise iteratively conditioned on the pair representation,
  end with clean coordinates. The same mathematical framework
  underlies image diffusion (Ho et al., 2020), score-matching (Song
  et al., 2021), and small-molecule conformation generation
  (#idx("DiffDock")DiffDock, 2022). The advantage over IPA is generative flexibility;
  the cost is sharper sensitivity to the conditioning signal.
]

AlphaFold-3 reports accuracy at roughly 60–80 % of experimental
quality on multi-chain and nucleic-acid tasks, depending on the
sub-task. The protein-only accuracy is comparable to AlphaFold-2's
ceiling. The drug-docking accuracy is the eye-catching result — for
proteins with sufficient experimental data on similar systems,
AF3 reproduces ligand poses at near-crystallographic resolution.

One controversy attached to the AF3 release: the model weights were
not publicly released at launch. DeepMind made the system available
through a web interface (the AlphaFold Server at EMBL-EBI) but kept
the weights proprietary. The community pushback was loud — losing
reproducibility was a real cost — and partial relaxations followed
through 2024 and 2025. The episode is a foretaste of how the
open-science-versus-dual-use tension is likely to play out for
frontier biology models.

=== RoseTTAFold and the Open-Source Ecosystem

While DeepMind shipped AlphaFold-2 weights with a non-commercial
licence and AlphaFold-3 weights with restrictions, David Baker's group
at the University of Washington built a parallel architecture and
released it openly. *RoseTTAFold* (Minkyung Baek and colleagues,
2021) is a three-track architecture — sequence, pair, and 3D tracks
updating each other — that lags AF2 by a few GDT_TS points but is
fully available with weights. RoseTTAFold2 (2023) added multi-chain
support; RoseTTAFold-All-Atom (2024) added ligands and nucleic acids.
For labs that could not run AlphaFold-2 commercially, RoseTTAFold was
the workhorse.

*OpenFold* (Gustaf Ahdritz and colleagues, 2022) is an open-source
PyTorch reimplementation of AlphaFold-2 with retrained weights. The
architecture is identical; the training corpus and procedure are
documented; the released code supports #idx("fine-tuning")fine-tuning and feature
probing. OpenFold has become the canonical platform for
academic-research extensions of AlphaFold-2 — anything that requires
gradient access to the model or modifications to the training loop.

*ColabFold* (Milot Mirdita and colleagues, 2022) is a usability
wrapper. It replaces AlphaFold-2's slow MSA generation (HHblits and
JackHMMER against very large databases) with MMseqs2-based search
against a pre-computed cluster, cutting inference time from hours to
minutes. For most laboratories, "running AlphaFold" actually means
running ColabFold in a Google Colab notebook. The democratisation
ColabFold provides cannot be overstated — it was the difference
between AlphaFold being a flagship-lab tool and AlphaFold being a
graduate-student tool.

=== #idx("ESMFold")ESMFold and Single-Sequence Prediction

*ESMFold* (Zeming Lin and colleagues, _Science_ 2023) takes a
different bet. Instead of using an MSA as input, it passes the single
query sequence through *ESM-2*, a large protein-language transformer
trained by masked-residue prediction on 65 million UniRef sequences.
The resulting per-residue embeddings carry most of the evolutionary
information an MSA would supply, learned implicitly during language-
model #idx("pretraining")pretraining. ESMFold projects the embeddings into a pair
representation, then runs an AlphaFold-style structure module on top.

The crucial property is *no MSA search at inference time*. The
expensive HHblits / JackHMMER step disappears. ESMFold runs about
60× faster than AlphaFold-2 on the same hardware. Accuracy is lower
than AF2 on proteins with deep MSAs, comparable or better on
orphan and metagenomic proteins where MSAs are shallow or absent,
and competitive on most designed proteins. The Meta-AI team used
ESMFold to predict structures for $approx 600$ million metagenomic
proteins, released as the *ESM Metagenomic Atlas* — the first
structural view of a large fraction of the unexplored microbial
proteome.

The successor, *ESM-3* (2024), is a multi-modal protein foundation
model that handles sequence, structure, and function tokens jointly
in a single transformer, with conditioning on partial inputs (e.g.,
"design a sequence with this active site"). ESM-3 is the most
direct attempt yet to make protein modelling resemble large language
modelling: pre-train one big model, prompt it for specific tasks at
inference time.

The practical post-2024 workflow has settled into a few patterns.
ColabFold (AlphaFold-2) is the default for most protein-only
prediction; ESMFold is the choice when speed matters or the MSA is
shallow; AlphaFold-3 (via the EMBL-EBI server) or RoseTTAFold-All-Atom
handles ligands and nucleic acids; OpenFold is the option when the
research question requires gradient access. Five years after CASP14
the field is broadly saturated near $approx 90$ GDT_TS on single-
domain folds. The remaining frontiers — multi-state proteins,
conformational ensembles, intrinsic disorder, and very large
assemblies — are where the next decade of work will land.


== Inverse Folding and Protein Design <sec:design>

Structure prediction runs the map *sequence $arrow.r$ structure*.
Inverse folding runs *structure $arrow.r$ sequence*. Given a desired
three-dimensional fold, what amino-acid sequence would produce it?

The problem is ill-posed in a specific way: many sequences can fold
to the same structure (the forward map is many-to-one), so the
inverse has infinite solutions. The interesting direction is to find
sequences that fold *stably and accurately* to the target. Three
canonical applications motivate the field. *Stabilisation* — given a
natural protein, find a sequence with the same fold but a higher
melting temperature, lower aggregation propensity, or longer shelf
life. *Functional engineering* — given a de-novo designed structure
with a specific binding pocket, find a sequence that folds to it.
*Repair* — given a protein with a deleterious mutation, find a
near-native sequence that fixes the defect.

=== #idx("ProteinMPNN")ProteinMPNN

*ProteinMPNN* (Justas Dauparas and colleagues, _Science_ 2022) is the
Baker-lab inverse-folding network. The input is a *backbone-only*
3D structure — the positions of N, $C_alpha$, $C$, and $O$ atoms for
every residue, with side chains and sequence identity stripped away.
The output is a probability distribution over the 20 amino acids at
every position; sequences are drawn by autoregressive sampling.

The architecture is a *message-passing graph neural network*. Each
residue is a node, with its backbone frame as the node feature. Edges
connect each residue to its $k$ nearest neighbours in 3D space.
Message-passing iterations exchange geometric information between
neighbours; a final per-residue softmax gives the amino-acid
distribution. The autoregressive decoder fills in residues one by
one, conditioning each prediction on the residues already decoded.

#figure(
  image("../../diagrams/lecture-15/09-proteinmpnn.svg", width: 95%),
  caption: [
    ProteinMPNN as a backbone-conditioned amino-acid GNN. Backbone
    geometry in, sequence out. Sampled sequences fold back to the
    target ~50 % of the time in experimental tests — well above
    Rosetta's $approx 10$ % baseline.
  ],
) <fig:proteinmpnn>

The empirical result that drove adoption: in head-to-head experimental
tests in the Baker lab, $approx 50$ % of ProteinMPNN-designed
sequences expressed solubly and folded correctly, against
$approx 10$ % for Rosetta-designed sequences targeting the same
backbones. The five-fold uplift in wet-lab success made protein
design a routine enterprise for the first time.

#note[
  Inverse folding is a classical ill-posed inverse problem — the
  forward map (sequence to structure) is many-to-one, so the
  inverse has infinite solutions. ProteinMPNN regularises by
  incorporating strong *structural priors* learned from millions of
  known backbone–sequence pairs. The autoregressive decoder adds
  *compositional regularisation*: residue identities are not sampled
  independently, each conditions on its already-decoded neighbours.
  The structure parallels image super-resolution (ill-posed; prior =
  natural-image statistics) and compressed sensing (ill-posed; prior =
  sparsity). The novelty is that the prior is learned end-to-end on
  the structural data, not hand-crafted.
]

=== #idx("RFdiffusion")RFdiffusion and the Design Loop

Pure inverse folding still requires a *target backbone* — you give
ProteinMPNN a shape and ask for a sequence. *De novo* protein design
is more ambitious: invent a shape that nature has never seen, then
find a sequence for it. The Rosetta era achieved this for restricted
topologies (the Koga et al. 2012 papers built de-novo four-helix
bundles and $beta$-barrels from idealised geometry), but expanding
the design space required a generative model over backbones.

*RFdiffusion* (Joseph Watson and colleagues, _Nature_ 2023) is that
generative model. It fine-tunes the RoseTTAFold network as a
diffusion model trained to denoise corrupted PDB structures. At
inference, you specify constraints — "a four-helix bundle binding
this ligand," "an active site that holds these three catalytic
residues in this geometry," "a binder for this antibody Fc region" —
and RFdiffusion generates novel backbone structures satisfying the
constraints. The design pipeline composes three stages
(@fig:rfdiffusion):

1. *RFdiffusion*: given constraints, sample a backbone.
2. *ProteinMPNN*: given the backbone, design a sequence.
3. *AlphaFold-2 / ESMFold validation*: fold the designed sequence
   and check it converges to the same backbone (self-consistency).

#figure(
  image("../../diagrams/lecture-15/10-rfdiffusion.svg", width: 95%),
  caption: [
    The design pipeline as sculptor (RFdiffusion) $arrow.r$
    scriptwriter (ProteinMPNN) $arrow.r$ critic (AlphaFold)
    $arrow.r$ wet lab. The pipeline produces de-novo proteins at
    10–40 % experimental-success rates against $< 1$ % with the
    pre-AlphaFold Rosetta-only flow.
  ],
) <fig:rfdiffusion>

Only candidates that survive self-consistency are synthesised. In the
2023 Watson paper, $approx 10$–$40$ % of synthesised candidates
expressed, folded, and exhibited the desired function. Pre-AlphaFold,
equivalent success rates were below 1 %. The decisive change is the
fast in-silico critic: before AlphaFold, there was no quick way to
check whether a designed sequence would fold as intended, so designs
mostly failed at the wet-lab step. Adding a reliable folding critic
inside the design loop is what flipped the success rate by a factor
of fifty.

#tip[
  The design loop is a useful structural pattern beyond proteins. It
  is a generator-critic relay: a flexible generative model proposes
  candidates, a separately-trained validator checks them, only the
  validated candidates are spent. The same shape appears in
  generative chemistry (a graph-VAE proposes molecules, DFT predicts
  binding energy), in code generation (an LLM proposes patches, a
  test suite validates), and in adversarial training (a generator
  proposes images, a discriminator scores). The lesson the protein-
  design field learned ahead of others is that the validator's
  accuracy is the bottleneck. Weak critic = expensive wet lab. Strong
  critic = fast iteration.
]

=== Dual Use

Protein design is a dual-use technology. Beneficial applications are
already well-established: therapeutic binders (Baker-lab anti-SARS-CoV-2
miniproteins, 2020; anti-influenza designs, 2024), novel enzymes for
green chemistry (PETase variants for plastic degradation), de-novo
vaccine antigens. Concerning applications run in parallel: a
sufficiently powerful generative model could in principle design
novel toxins or pathogen proteins.

The field has begun to grapple with the implications. The Baker lab
adopted a screening policy that excludes designs with predicted
toxic activity. Several major model releases (Meta's ESM-3, DeepMind's
AlphaFold-3) underwent some form of biosecurity review before public
release. The Anthropic ASL-3 framework, drafted in 2024, explicitly
cites bioweapon uplift as a frontier-AI threshold. The US Executive
Order on AI (October 2023) requires major model developers to notify
the government of frontier biological capabilities. Whether weights
are publicly released is now a deliberate policy decision rather than
a default; the decision is unlikely to converge in a single direction
across labs and countries.

If you find yourself building one of these models — and many readers
of this book will — the obligation is not to wait for the policy
framework to mature before thinking about it. Know which classes of
designs your model might enable that didn't exist before. Know who
benefits from open release and who is harmed. Decide deliberately.

#note[
  Concretely: before releasing weights for a new design model, ask
  what *known* dangerous proteins your model can improve upon
  relative to the baseline; whether your training data was screened
  for explicitly dual-use sequences; whether you can detect misuse
  after release; and whether the marginal scientific benefit of open
  weights outweighs the marginal uplift to bad actors. There is no
  formula. There is the obligation to ask.
]


== The AlphaFold Database and What Changed <sec:afdb>

By the time AlphaFold-2's weights and code were released in July
2021, the field had been preparing for a year. Within twelve months
DeepMind and EMBL-EBI had jointly released the *AlphaFold Protein
Structure Database*: predicted structures for the entire human
proteome, then the proteomes of twenty major model organisms, then
all of UniRef50, then nearly all of UniProt. By 2024 the database
held about 214 million predicted structures (@fig:afdb-scale). For
comparison, the Protein Data Bank — the cumulative output of fifty
years of crystallography and cryo-EM — held about 180,000
experimentally determined structures. The ratio flipped from less
than 0.1 % of proteins having any structural information to more
than 99 %.

#figure(
  image("../figures/ch15/f4-afdb-scale.svg", width: 95%),
  caption: [
    The AlphaFold Database in scale and confidence. About 1,200
    times more predicted structures than experimentally determined
    ones; roughly two-thirds at confident-or-better pLDDT; four
    routine workflows transformed from "rarely possible" to
    "routine" between 2020 and 2024.
  ],
) <fig:afdb-scale>

The database is searchable by UniProt accession, by pLDDT
distribution, by length, and (with companion tools like *Foldseek*,
Martin Steinegger's structural-search system) by 3D similarity to a
query structure. Each entry comes with the structure coloured by
pLDDT and the PAE matrix.

#figure(
  image("../../diagrams/lecture-15/11-af-database.svg", width: 85%),
  caption: [
    PDB versus AlphaFold DB at scale. The bar chart is log-scaled
    because the linear comparison fails to fit on a page — the AFDB
    is roughly three orders of magnitude larger.
  ],
) <fig:afdb>

Four downstream workflows were transformed within twenty-four months
of the release.

*Drug discovery.* Pre-AlphaFold, structure-based drug design required
an experimentally determined target structure — X-ray crystallography
or cryo-EM, years of work per target. Most drug targets had no
structure. Post-AlphaFold, every human protein has a predicted
structure ready for inspection. Druggability prediction (does this
protein have a viable binding pocket?), virtual screening (which of
ten million library compounds binds this pocket?), and structure-based
lead optimisation now begin from an AlphaFold model as a routine
input step. The pharmaceutical industry absorbed this rapidly; every
major pharma now has an AlphaFold-anchored workflow for target
validation.

*Functional annotation of dark proteomes.* Roughly a third of human
proteins remained functionally uncharacterised in 2020 — the
"dark proteome." A predicted structure plus structural homology
search (Foldseek, Dali) allows function transfer from
structurally-similar characterised proteins even when sequence
similarity is too low for traditional homology methods. The fraction
of human proteins with at least a putative function climbed
noticeably in the 2022–2024 window because of this.

*Cross-species structural comparison*. Asking how an enzyme has
evolved across a clade used to require either a curated set of
crystal structures (rare) or homology models (unreliable on diverged
sequences). Now every species' version of the enzyme has a
high-confidence predicted structure, and comparing them is routine.
Structural conservation invisible at the sequence level becomes
visible.

*De-novo design has a critic*. As laid out in @sec:design,
RFdiffusion plus ProteinMPNN plus AlphaFold-as-validator is the
modern design pipeline. Without a fast accurate folding critic the
loop falls apart; AlphaFold is that critic.

=== What's Still Hard

Three classes of biology still resist the AlphaFold workflow.

*Conformational dynamics.* A protein is not a single static structure;
it is a thermodynamic ensemble fluctuating around one or more energy
minima. AlphaFold predicts a dominant state but not the distribution.
Some drugs (kinase inhibitors, allosteric modulators) bind transient
conformations the predictor doesn't surface. Recent work — AlphaFold-2
with MSA subsampling (Wayment-Steele et al., 2024), *AlphaFlow* (Jing
et al., 2024), Distributional Graphormer — addresses this partially
by sampling diverse predictions or by training generative models over
conformational distributions, but no method is yet a full
substitute for experimental dynamics data.

*Intrinsically disordered regions.* Roughly 30 % of eukaryotic
proteins contain significant disordered content — stretches that
never fold into a stable structure in isolation. AlphaFold correctly
identifies these as low-pLDDT regions, but the prediction is the
*absence* of structure rather than a distribution over the states
the disordered region actually visits. Many disordered regions are
biologically functional (transcription-factor activation domains,
signalling motifs). Predicting their behaviour requires methods that
treat the ensemble as the object of prediction, not the static state.

*Multi-state functional proteins.* Kinases, transporters, and most
membrane proteins cycle through distinct functional states — active
and inactive conformations, open and closed channel states.
AlphaFold predicts one state, usually the one most represented in
its training data, and is silent on the others. For drug-discovery
work where the relevant target is an inactive kinase that the
network has predicted in its active form, the prediction can be
worse than useless.

#warn[
  AlphaFold predictions are *not experimental structures*. For
  publication, regulatory, or clinical work that requires real
  structural data — atomic drug-binding modes, mutagenesis
  experiments targeting specific atoms, structure–function causation
  claims — the prediction is a starting hypothesis, not a definitive
  answer. Always cross-check with experiment when the stakes are
  real, and read the pLDDT and PAE before accepting a prediction at
  face value.
]


== Summary <sec:ch15-summary>

- Protein function follows from three-dimensional structure, and
  Anfinsen's 1961 refolding experiment established that the structure
  is determined by the sequence. The fifty years from Anfinsen to
  AlphaFold-2 were spent searching for the algorithm.
- The classical era (homology modelling, threading, fragment assembly)
  plateaued at ~40 GDT_TS on hard CASP targets because it could not
  exploit the coevolutionary signal in homologous sequence families.
- Direct Coupling Analysis (Morcos 2011, Marks 2011) read contacts
  from MSAs by inverse-covariance estimation, the discrete-variable
  cousin of the graphical lasso. DCA + Rosetta dominated CASP12
  before deep learning displaced it.
- AlphaFold-2 (Jumper 2021) replaced the contact-then-fold two-step
  pipeline with an end-to-end network. The Evoformer applies 48 blocks
  of axial attention and triangle updates to refine joint (MSA, pair)
  representations; the structure module decodes coordinates with
  $"SE"(3)$-equivariant Invariant Point Attention; recycling lets
  the network correct itself across multiple forward passes.
- Per-residue pLDDT and per-pair PAE are calibrated confidence metrics.
  Reading them together is non-negotiable — a high-pLDDT protein can
  still have its two halves mis-oriented if PAE is high across the
  inter-domain hinge.
- AlphaFold-3 (Abramson 2024) replaces IPA with a diffusion-based
  structure module and adds nucleic-acid and ligand inputs.
  ESMFold (Lin 2023) replaces the MSA with a protein-language-model
  embedding, trading some accuracy on deep families for a 60×
  speedup and competitive accuracy on orphan proteins.
- ProteinMPNN (Dauparas 2022) solves inverse folding — backbone in,
  sequence out — by a graph neural network with autoregressive
  decoding. RFdiffusion (Watson 2023) generates novel backbones from
  constraint specifications by denoising diffusion. Together with
  AlphaFold as a validation critic, they pushed de-novo design
  experimental-success rates from $< 1$ % to 10–40 %.
- The AlphaFold Protein Structure Database holds about 214 million
  predicted structures — roughly 1,200 times more than the cumulative
  experimental PDB. Drug discovery, functional annotation,
  cross-species comparison, and de-novo design all assume an
  available predicted structure as input.
- Open frontiers: conformational dynamics, intrinsic disorder,
  multi-state functional proteins, large multi-chain assemblies, and
  the open-science-versus-dual-use tension in releasing model
  weights.


== Exercises <sec:ch15-exercises>

#strong[1.] #emph[Reading pLDDT and PAE.] Download an #idx("AlphaFold-DB")AlphaFold-DB
prediction for a multi-domain human protein (kinase, ABC transporter,
or your favourite). Open the structure in a viewer that colours by
pLDDT (PyMOL, ChimeraX). Examine the PAE matrix. Identify (a) the
domain boundaries from the PAE block structure, (b) the highest-confidence
residue stretch, (c) the lowest-confidence stretch and whether it
corresponds to a likely disordered region (cross-check DisProt or
IUPred). Write one paragraph summarising which parts of the prediction
you would and would not trust for a downstream task.

#strong[2.] #emph[Counting Levinthal.] A 150-residue protein has, for
each backbone, three accessible $phi.alt$ and three accessible $psi$
combinations. (a) How many backbone conformations are possible if
torsions are sampled independently? (b) At one conformation per
nanosecond, how long would an exhaustive search take? (c) Real proteins
of this size fold in 10 ms to 1 s. Estimate the maximum effective
branching factor per residue that is consistent with that folding time,
assuming a sequential search. Compare to your answer in (a).

#strong[3.] #emph[Mutual information versus DCA.] Generate a synthetic
MSA from a Potts model with three known contacts: $(5, 10)$, $(10, 15)$,
$(20, 30)$. Use 200 sequences with mutation rates that mix the contacts
about evenly. Compute the mutual-information matrix and an estimate of
the DCA contact scores (the `pydca` package implements plmDCA; or
write a mean-field DCA in 50 lines of numpy). Plot both matrices. Which
contacts does MI find? Which extra (indirect) contacts does it report?
Which contacts does DCA find?

#strong[4.] #emph[Axial-attention cost.] An MSA of 4096 sequences and
2048 residues passes through one Evoformer block per inference. (a)
Compute the per-block FLOP cost of full 2-D attention. (b) Compute the
cost of axial attention (row plus column). (c) The Evoformer is 48
blocks; what is the inference-time savings, in TFLOPs, of axial over
full attention?

#strong[5.] #emph[FAPE versus RMSD.] Two predicted structures of a
40-residue $alpha$-helical peptide both have backbone RMSD = 2.0 Å
versus the experimental structure when superimposed globally. Structure
A has uniform 2 Å error spread along the helix; structure B has zero
error in the first 30 residues and 8 Å error concentrated in the last
10. Compute, by hand or with a short Python script, an approximate FAPE
for each (use $epsilon = 0$ in the clamping). Which structure is the
better prediction by FAPE? Which by RMSD? Discuss when each is the
more useful objective.

#strong[6.] #emph[Inverse-folding by template-matching.] Pick a small
PDB structure (50–100 residues) and strip the sequence (keep only
backbone atoms). Use ProteinMPNN (via the ColabDesign notebook) to
sample 20 sequences for the stripped backbone. (a) Compute the
sequence-recovery rate (fraction of positions matching the native).
(b) Fold each designed sequence with AlphaFold-2 (or ESMFold) and
compute self-consistency: the backbone RMSD of the prediction against
the original. (c) Plot the relationship between sequence-recovery
and self-consistency RMSD. Are higher-recovery sequences more
self-consistent?

#strong[7.] #emph[Reading the AlphaFold-2 paper.] Pick one Evoformer
sub-operation from @fig:evoformer — row attention with pair bias,
column attention, outer-product mean, triangle multiplicative update,
or triangle attention. From the supplementary methods of Jumper et al.
2021, write a one-page explanation in your own words, including the
tensor shapes at input and output and a complexity estimate per
forward pass.

#strong[8.] #emph[(Open-ended.)] Pick one of the open frontiers from
@sec:afdb (conformational dynamics, intrinsic disorder, multi-state
function, large assemblies, biosecurity). Find one paper from 2023–2025
that addresses it, read the abstract and figure 1, and write a single
paragraph describing the central trick the paper introduces and what
empirical evidence supports it. Cite the paper.


== Further Reading <sec:ch15-further-reading>

- *Anfinsen, C. B.* (1973). "Principles that govern the folding of
  protein chains." _Science_ 181: 223–230. The original statement
  of the thermodynamic hypothesis. Still readable in an afternoon.
- *Marks, D. S., Colwell, L. J., Sheridan, R., et al.* (2011).
  "Protein 3D structure computed from evolutionary sequence
  variation." _PLoS ONE_ 6: e28766. The EVfold paper — coevolution
  as a predictor of 3D contact, before deep learning.
- *Morcos, F., Pagnani, A., Lunt, B., et al.* (2011). "Direct-coupling
  analysis of residue coevolution captures native contacts across
  many protein families." _PNAS_ 108: E1293–E1301. The companion DCA
  paper; the inverse-covariance formulation laid out cleanly.
- *Jumper, J., Evans, R., Pritzel, A., et al.* (2021). "Highly
  accurate protein structure prediction with AlphaFold." _Nature_
  596: 583–589. The AlphaFold-2 paper. The main text is approachable;
  the supplementary methods reward careful reading.
- *Abramson, J., Adler, J., Dunger, J., et al.* (2024). "Accurate
  structure prediction of biomolecular interactions with AlphaFold-3."
  _Nature_ 630: 493–500. The AlphaFold-3 paper, including the
  diffusion structure module and multi-modal input handling.
- *Dauparas, J., Anishchenko, I., Bennett, N., et al.* (2022).
  "Robust deep-learning-based protein sequence design using
  ProteinMPNN." _Science_ 378: 49–56. The inverse-folding paper that
  reset the design-success baseline.
- *Watson, J. L., Juergens, D., Bennett, N. R., et al.* (2023). "De
  novo design of protein structure and function with RFdiffusion."
  _Nature_ 620: 1089–1100. The generative-design paper; pair with
  the Dauparas paper to read the design pipeline end to end.
- *Lin, Z., Akin, H., Rao, R., et al.* (2023). "Evolutionary-scale
  prediction of atomic-level protein structure." _Science_ 379:
  1123–1130. ESMFold and the Metagenomic Atlas — single-sequence
  prediction with a protein-language-model front end.
- *Varadi, M., Anyango, S., Deshpande, M., et al.* (2022). "AlphaFold
  Protein Structure Database: massively expanding the structural
  #idx("coverage")coverage of protein-sequence space with high-accuracy models."
  _Nucleic Acids Research_ 50: D439–D444. The database paper.
- *Mirdita, M., Schütze, K., Moriwaki, Y., et al.* (2022). "ColabFold:
  making protein folding accessible to all." _Nature Methods_ 19:
  679–682. The practical workhorse most labs actually use.
