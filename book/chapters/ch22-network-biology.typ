#import "../theme/book-theme.typ": *

= Network Biology and Pathway Analysis <ch:network-biology>

#matters[
  By the time you reach this chapter, the previous twenty-one have
  taught you to produce gene lists. #idx("differential expression")Differential expression hands you
  a few hundred up- and down-regulated transcripts. #idx("ChIP-seq")ChIP-seq returns
  TF-bound regions; #idx("GWAS")GWAS returns a set of significant SNPs and the
  genes nearest to them; somatic-variant calling on tumour cohorts
  returns a tail of mutated drivers. In every case the gene list is the
  output of the assay, not the answer to the biological question. The
  question is what the list _means_, and the framework that answers it
  is built on graphs. Pathway databases attach biological vocabulary to
  the list. Enrichment tests ask whether one of those vocabularies is
  over-represented. Network propagation spreads the signal across known
  interactions so the analysis sees a connected mechanism rather than
  isolated genes. Module detection partitions the resulting subgraph
  into biological units. All of it — every algorithm in this chapter —
  is graph signal processing applied to biology, and an EE student
  already knows the math.
]

You arrive at network biology with a list of genes and a question that
the list alone cannot answer. The genes might be the DEGs from a
two-condition #idx("RNA-seq")RNA-seq comparison, the marker set for a single-cell
cluster, the candidate region from a clinical exome, or the
significant SNPs from a biobank-scale GWAS. They are unsorted by
biological context, and the biological-context question — _what
pathway is broken? which complex is perturbed? which mechanism is
implicated?_ — is the one your collaborator will ask you over coffee
on the day you deliver the table.

This chapter is the framework that answers that question. The first
move is to recognise that the gene is not a self-contained unit;
proteins interact, #idx("transcription")transcription factors regulate, metabolites
participate in cascades. The natural representation is a graph. The
second move is to attach a controlled vocabulary to the graph: a
_pathway_ in #idx("KEGG")KEGG, a _term_ in #idx("Gene Ontology")Gene Ontology, a _gene set_ in MSigDB.
The third move is to ask which terms are statistically enriched in
your list — by counting (over-representation analysis) or by
following the full ranked list (gene-set enrichment analysis). The
fourth move is to propagate signal across the network so a gene that
is connected to many DEGs scores high even if its own expression did
not change. The fifth and final move is to cut the resulting
subgraph into modules — communities — that correspond to functional
units.

Across all five moves the same EE-flavoured mathematics keeps
recurring. The graph #idx("Laplacian")Laplacian is a spectral operator with a Fourier
basis. Network propagation is a low-pass filter in that basis.
#idx("spectral clustering")Spectral clustering operates in the same basis on a different
slice. #idx("GSEA")GSEA is a Kolmogorov–Smirnov statistic on a ranked list.
#idx("PageRank")PageRank — the founding algorithm of modern web search — is exactly
the random-walk-with-restart that drug-target prediction has been
running since the 2000s. The bioinformatics literature occasionally
re-invents these tools and gives them new names; the EE-trained
reader can usually recognise the underlying machinery on first sight.


== From a Gene List to a Network <sec:from-list-to-graph>

The gene list arrives in different shapes from different upstream
analyses, but the downstream question is uniform. Five places it
shows up in this book are worth listing because they look
superficially different and converge on the same network step.

- Differential expression on bulk RNA-seq (Chapter 6) returns a set
  of genes whose #idx("transcript abundance")transcript abundance differs significantly between
  conditions, typically a few hundred entries at #idx("FDR")FDR &lt; 0.05.
- ChIP-seq #idx("peak calling")peak calling (Chapter 9) returns transcription-factor
  binding sites which, after annotation, become a set of target
  genes.
- scRNA-seq clustering (Chapter 7) labels each cluster with a set of
  marker genes that distinguish it from its neighbours.
- GWAS (Chapter 13) returns a set of SNPs above genome-wide
  significance, and a standard interpretation pass maps each #idx("SNP")SNP to
  the nearest protein-coding gene.
- Tumour driver discovery on cancer cohorts returns a set of
  significantly mutated genes — the canonical drivers and the long
  tail of rare ones.

In each case the gene names are biologically meaningful but the list
itself is structurally flat. The genes are not annotated by which
pathway they belong to, which complex they participate in, which
disease they are associated with. The first move is therefore to add
that structure, and the natural object that carries it is a graph.

There are five distinct kinds of biological network you will
encounter, each with its own data sources and uses, and the
difference between them is more than vocabulary.

*Protein–protein interaction (PPI) networks* are the workhorse. Nodes
are proteins; edges are physical interactions detected by
affinity-purification #idx("mass spectrometry")mass spectrometry (AP-MS), yeast two-hybrid
(Y2H), proximity labelling (BioID, APEX), or structural prediction
(AlphaFold-Multimer). The major databases — #idx("STRING")STRING (von Mering et
al., 2003), #idx("BioGRID")BioGRID, IntAct, the Human Reference Interactome (HuRI) —
each curate the same underlying experiments with different policies
on what to include. STRING is broadest because it includes predicted
interactions weighted by confidence; BioGRID is narrower because it
restricts to curated experimental evidence. PPI networks are used to
infer protein-complex membership and to apply _guilt by association_
— a protein of unknown function whose neighbours are all DNA-repair
factors is probably itself a DNA-repair factor.

*Gene regulatory networks (GRNs)* are directed. Nodes are
transcription factors and their target genes; edges represent the
regulation of one by the other. The evidence comes from ChIP-seq #idx("motif")motif
calls, perturbation experiments — #idx("CRISPRi")CRISPRi screens that knock down
each TF in turn and measure the downstream cascade — and from
#idx("ATAC-seq")ATAC-seq combined with chromatin-state annotation. GRNs are
inherently noisier than PPI networks because regulation is contextual:
a TF binds #idx("chromatin")chromatin in some cell types and not others, and the same
binding event activates one target while repressing another.

*Metabolic networks* are bipartite. Nodes are either metabolites or
enzymes; edges represent participation of one in a reaction catalysed
by the other. KEGG, MetaCyc, and Recon3D are the standard databases.
The metabolic network is what makes flux-balance analysis possible:
treat the network as a linear system, write the steady-state
constraint that mass is conserved at every metabolite, and solve a
linear program to predict reaction fluxes from growth conditions.

*Signalling networks* are also directed and overlap with GRNs at the
TF-binding end. Nodes are signalling proteins — kinases, receptors,
adapters; edges are #idx("phosphorylation")phosphorylation events, allosteric activations,
or stable complex formation. PhosphoSitePlus, OmniPath, and #idx("Reactome")Reactome's
signalling-pathway portion are the main sources. Signalling networks
are where drug-mechanism inference and pathway-perturbation modelling
live.

*Disease networks* are the most abstract. Nodes are diseases; edges
encode shared genes, comorbidity in clinical records, or drug-target
overlap. DisGeNET, OpenTargets, and OMIM are the standard catalogues.
These networks are how _drug repositioning_ work begins — find two
diseases connected in the network, ask whether a drug effective on
one might be relevant to the other.

#figure(
  image("../../diagrams/lecture-22/01-network-types.svg", width: 95%),
  caption: [
    Biological network taxonomy. Five major network types — PPI, GRN,
    metabolic, signalling, disease — each with their own data
    sources, scale, and applications. The same mathematical
    framework, the graph Laplacian and its spectrum, applies to all
    of them.
  ],
) <fig:net-types>

#note[
  The five network types are not exclusive. Most working analyses
  build a _heterogeneous_ network that fuses two or more — drugs and
  their targets glued onto a PPI, GWAS hits projected onto a regulatory
  network, ligand–receptor pairs joining single cells into a
  communication graph. The graph algorithms below all extend
  naturally to heterogeneous edges; the only book-keeping change is
  that edge weights carry a type.
]

=== Topology Shared with the Web and the Citation Graph

Real biological networks share three topological features that are
not unique to biology and that the network-science literature has
characterised carefully. The first is *scale-free degree
distribution*: a handful of hubs have very high degree and most
nodes have low degree, with the distribution following a power law
$P(k) tilde k^(-gamma)$ with $gamma$ typically in the range two to
three. Barabási and Albert (1999) introduced the preferential-attachment
model that generates such distributions — _the rich get richer_,
because a new node is more likely to attach to an already
well-connected one. Human PPI networks have $gamma$ near $2.3$;
hubs include TP53, EGFR, and ACTB, each with several hundred
interaction partners.

The second is *small-world structure*: the average shortest-path
length grows logarithmically with the network size. The Watts–Strogatz
model (1998) showed that a small fraction of long-range "shortcut"
edges in an otherwise locally clustered graph is enough to produce
this regime. Human PPI networks are small-world with an average
path length around four to six — any two proteins are typically
connected by a chain of four to six interactions. The implication
for analysis is that propagation spreads quickly. A signal seeded at
ten DEGs reaches most of the proteome within a handful of hops; the
analytical challenge is _which_ neighbours matter, not how to reach
them.

The third is *#idx("modularity")modularity*: dense functional clusters connected to
each other by sparse inter-cluster edges. The cluster structure is
itself nested — protein complexes inside pathways inside larger
functional systems. Module detection in Section 22.5 is the
machinery for recovering these clusters from the graph alone.

None of these properties are unique to biology. Social networks,
the Web, citation graphs, the protein interactome — all show
scale-free degree distributions, small-world path lengths, and
modular structure. The toolkit of network science transfers
directly, which is why so many of the algorithms in this chapter
have their origin outside biology and arrived via #idx("translation")translation.

=== Networks as Graph Signal Processing

The mathematics that organises this chapter is the spectral theory of
the graph Laplacian. A network $G = (V, E)$ with $|V| = n$ nodes has
an adjacency matrix $A in {0, 1}^(n times n)$ where $A_(i j) = 1$
iff $(i, j) in E$, and a degree matrix $D = "diag"(d_1, dots, d_n)$
with $d_i$ the number of neighbours of node $i$. The *graph
Laplacian* is

$ L = D - A. $

The normalised variants $L_("sym") = I - D^(-1/2) A D^(-1/2)$ and
$L_("rw") = I - D^(-1) A$ are used in practice because they handle
degree variation better, but $L$ itself is the cleanest object to
think with.

$L$ is symmetric and positive semi-definite. Its eigenvalues are
real and non-negative, $0 = lambda_1 lt.eq lambda_2 lt.eq dots lt.eq
lambda_n$, with $lambda_1 = 0$ corresponding to the constant
eigenvector. The remaining eigenvectors $bold(v)_2, dots, bold(v)_n$
are orthogonal and form a basis for $RR^n$. A signal on the graph
is a function $f: V -> RR$ — for example, expression level at each
node — and any such signal has a unique decomposition

$ f = sum_(i = 1)^n hat(f)_i bold(v)_i, quad hat(f)_i = angle.l f,
bold(v)_i angle.r. $

This is the *graph Fourier transform*. Low-eigenvalue eigenvectors
are smooth across edges — neighbours have similar values; the
extreme case is the constant eigenvector at $lambda = 0$.
High-eigenvalue eigenvectors are rough — values change sign between
adjacent nodes; the extreme case at $lambda_n$ alternates as much
as the graph topology allows.

#figure(
  image("../figures/ch22/f1-graph-laplacian.svg", width: 95%),
  caption: [
    The graph Laplacian as a spectral operator. A small network with
    two communities (left), its adjacency and Laplacian matrices
    (middle), and three of its eigenvectors plotted as signals on
    the node index (right). The smallest non-zero eigenvector $bold(v)_2$
    splits the graph along the bridge edge — the canonical "Fiedler
    vector" cut. High-eigenvalue eigenvectors oscillate rapidly between
    neighbours and have no community-level interpretation.
  ],
) <fig:laplacian>

This decomposition is the engineering #idx("scaffold")scaffold for the rest of the
chapter. Network propagation in Section 22.4 is a low-pass filter:
multiply each Fourier coefficient $hat(f)_i$ by a decay factor that
shrinks the high-eigenvalue components and keeps the low-eigenvalue
ones. Spectral clustering in Section 22.5 is the dual move: take
the lowest-eigenvalue eigenvectors and use them as coordinates for
k-means. Both algorithms are operations on the graph Fourier
spectrum, dressed differently.

#note[
  The same mathematics is the foundation of the recently mature
  field of _graph signal processing_ (Shuman et al., 2013). The
  bioinformatics community arrived at the moves independently in the
  2000s, motivated by network propagation rather than by an explicit
  Fourier analogy. Reading the original Vanunu et al. (2010) PRINCE
  paper alongside Shuman's GSP survey makes the two-cultures
  convergence vivid: same equations, different vocabulary,
  developed in parallel.
]


== Pathway Databases <sec:pathway-databases>

Before any algorithm can attach biological meaning to a gene list, it
needs a vocabulary of biological _categories_ — pathways, complexes,
processes, locations, phenotypes — and a mapping from each gene to
the categories it belongs to. Four databases have become the
practical standard, each curated under different assumptions about
what a "pathway" is.

*KEGG*, the Kyoto Encyclopedia of Genes and Genomes (Kanehisa,
*Nucleic Acids Research*, 2000), is the oldest. Established in 1995
at the Kyoto University Bioinformatics Center, it provides roughly
550 manually-curated pathway maps — the iconic ball-and-stick
diagrams that decorate biology textbooks. KEGG pathway IDs (e.g.,
`hsa04110` for "cell cycle, human") are the de-facto standard for
pathway-level annotation. KEGG is strongest in metabolism, where the
ball-and-stick maps directly encode the chemistry; it is weaker in
signalling, where the level of abstraction does not always match
modern molecular detail. Licensing has become restrictive in the
last decade — academic use is free, commercial use is not — and
this has driven some downstream pipelines toward open alternatives.

#figure(
  image("../../diagrams/lecture-22/07-kegg-glycolysis.svg", width: 90%),
  caption: [
    KEGG glycolysis pathway. The ball-and-stick representation is the
    standard visual language for metabolic pathways: nodes are
    metabolites, edges are enzymes labelled by EC number, and
    branches show inputs to the pentose-phosphate pathway and
    outputs to lactate or the #idx("mitochondrion")mitochondrion.
  ],
) <fig:kegg>

*Reactome* (Joshi-Tope et al., *Nucleic Acids Research*, 2005) is the
younger and more granular competitor. Established in 2003 as an
open-source successor to the older Genome-Net pathway maps, it
provides roughly 2,500 human pathways, hierarchically organised so
that broad pathways like "#idx("DNA")DNA repair" nest into specific
sub-pathways like "translesion synthesis" and "mismatch repair".
Each pathway is decomposed into reactions; each reaction lists
inputs, outputs, catalysts, and the regulators that promote or
inhibit it. Reactome's hierarchy means an enrichment analysis can
run at the granularity that fits the question — broad terms for
exploratory work, specific terms when the biology is well understood.

*Gene Ontology* (Ashburner et al., *Nature Genetics*, 2000) takes a
different approach entirely. Rather than enumerate pathways, GO
defines a controlled vocabulary of biological concepts organised into
three directed acyclic graphs: _biological process_ (e.g., "cell
cycle progression"), _molecular function_ (e.g., "kinase activity"),
and _cellular component_ (e.g., "nuclear envelope"). The current
release contains about 50,000 terms, and each gene is annotated to
multiple terms — often dozens — by curators reading the primary
literature. The DAG structure encodes _is-a_ and _part-of_
relationships, so "cell-cycle phase transition" is a child of "cell
cycle" is a child of "cellular process" is a child of "biological
process". The _true-path rule_ says that a gene annotated to a deep
term is implicitly annotated to all of its ancestors, which has
consequences both for enrichment statistics and for redundancy
in the output.

#figure(
  image("../../diagrams/lecture-22/08-go-dag.svg", width: 90%),
  caption: [
    A subset of the Gene Ontology biological-process DAG. Each term
    has parent–child relationships; under the true-path rule, a gene
    annotated to a deep term is implicitly annotated to all of its
    ancestors. The DAG structure makes parent and child enrichment
    p-values correlated, which Section 22.3 addresses.
  ],
) <fig:go-dag>

*MSigDB*, the Molecular Signatures Database from the Broad Institute,
is an aggregator rather than a primary source. It collects roughly
30,000 gene sets from KEGG, Reactome, BioCarta, curated cancer
signatures, and the curated _Hallmark_ collection of 50 well-defined
functional sets. For a typical RNA-seq differential-expression
follow-up, MSigDB Hallmark plus KEGG plus Reactome is the standard
first pass.

#figure(
  image("../../diagrams/lecture-22/02-pathway-databases.svg", width: 95%),
  caption: [
    The pathway-database landscape. Each database has its own
    granularity and curation philosophy; the same gene appears in
    different "pathways" across them. Modern analyses run several
    in parallel and report consensus rather than a single source.
  ],
) <fig:databases>

#tip[
  When a database disagreement matters for interpretation — say,
  KEGG implicates one pathway and Reactome implicates two
  different ones — the right move is to look at the gene-level
  overlap, not at the pathway names. The pathways labelled "DNA
  repair" in different databases share most of their genes; the
  apparent disagreement is usually a vocabulary mismatch over a
  shared mechanism.
]

=== The Gene-Set Universe Has Internal Structure

A subtle point that catches new analysts is that the gene-set
universe has _internal structure_ which biases naive enrichment
statistics. Two issues recur. The first is overlap: pathways share
genes, and GO terms inherit annotations from their descendants. A
significant enrichment of "cell-cycle phase transition" plus a
significant enrichment of "cell cycle" is one biological finding
counted twice. The second is curation bias: better-studied genes
are annotated to more pathways than poorly-studied ones, simply
because more papers have been written about them. TP53 sits in 14
KEGG pathways, 87 Reactome pathways, 198 GO biological-process
terms, and 6 Hallmark sets. A novel uncharacterised protein may
appear in none. Enrichment tests that ignore this asymmetry
systematically over-call well-studied genes.

Modern analyses confront both problems explicitly. Redundancy
reduction tools — REVIGO (Supek et al., 2011), simplifyEnrichment
(Gu & Hübschmann, 2023) — cluster significant terms by semantic
similarity and report cluster representatives rather than raw
lists. Hub correction at the network step (Section 22.7) compensates
for the annotation asymmetry by computing a null distribution that
preserves degree and annotation density. Neither remedy is perfect;
both make the difference between a publishable interpretation and
a noise-driven artefact.


== Pathway Enrichment <sec:enrichment>

Pathway enrichment is the simplest network-aware analysis, and the
one you will run first on every gene list. Two flavours dominate:
*over-representation analysis* (#idx("ORA")ORA), which discretises the gene
list at a significance threshold and counts overlap, and *gene-set
enrichment analysis* (GSEA), which uses the full ranked list of all
tested genes.

=== Over-Representation Analysis

Set up the problem with three sets. $U$ is the *universe* — the set
of genes that were tested by your experiment; for an RNA-seq run,
that is the set of genes that passed your expression filter and
entered the differential test. $L$ is your candidate list — the
DEGs at FDR &lt; 0.05, the GWAS-significant SNPs' nearest genes,
the cluster's marker genes. $G$ is a gene set you want to test —
a KEGG pathway, a GO term, a Hallmark set.

The null hypothesis is that $L$ is a uniform random subset of $U$
of size $|L|$. Under that null, the overlap count $|L sect G|$ is
*#idx("hypergeometric")hypergeometric* — the same distribution as drawing $|L|$ balls
without replacement from an urn of $|U|$ balls of which $|G|$ are
marked. The one-sided p-value is the probability of seeing at least
the observed overlap by chance,

$ p = sum_(k = |L sect G|)^(min(|L|, |G|)) frac(binom(|G|, k)
binom(|U| - |G|, |L| - k), binom(|U|, |L|)). $

Equivalently, this is #idx("Fisher's exact test")Fisher's exact test applied to the $2 times
2$ contingency table of "in / not in $L$" against "in / not in
$G$", and the equivalence is exact rather than approximate.

#figure(
  image("../figures/ch22/f2-ora-worked.svg", width: 95%),
  caption: [
    A worked ORA example. With 200 DEGs from a universe of 20,000
    tested genes and a 500-gene pathway, the expected overlap under
    the null is 5; the observed overlap is 20, four times the
    expected, and the hypergeometric p-value is approximately
    $4 times 10^(-10)$. The footer lists the four ways the same
    calculation is most often misused.
  ],
) <fig:ora>

#figure(
  image("../../diagrams/lecture-22/03-ora-vs-gsea.svg", width: 95%),
  caption: [
    ORA versus GSEA. ORA discretises the gene list at the
    significance threshold and counts the overlap with the pathway;
    GSEA uses the full ranked list and asks whether pathway members
    cluster toward the extremes of the ranking. GSEA's strength is
    catching coherent pathway-wide shifts that ORA misses when most
    pathway genes are just below the cutoff.
  ],
) <fig:ch22-ora-vs-gsea>

ORA is fast, simple to implement, and easy to misuse. The most
common failure is the *wrong-universe* error: setting $|U|$ to "all
human genes" when only a fraction were actually tested in your
experiment. With $|U|$ inflated, every gene set's expected overlap
shrinks toward zero, every observed overlap looks improbable, and
every pathway you test comes back significant. The fix is to set
$|U|$ to the set of genes that survived the expression filter and
entered the test. A second failure is *unused background*: if your
candidate list comes from a targeted panel (a #idx("CRISPR")CRISPR screen against
500 cancer genes, say), the universe must be restricted to those
500 genes, not extended to the whole transcriptome.

The third failure is *small expected counts*. When the expected
overlap drops below about three, the asymptotic chi-squared
approximation to Fisher's test breaks down and the asymptotic
p-value can be wrong by an order of magnitude. The exact
hypergeometric is the safe default; modern enrichment packages
(`clusterProfiler`, `g:Profiler`, `enrichR`) use it by default and
the chi-squared shortcut survives only in older code.

==== #idx("multiple testing")Multiple Testing

Testing 1,000 pathways against the same list produces ~50 nominally
"significant" results at $p < 0.05$ purely by chance. The standard
remedy is Benjamini–Hochberg false-discovery-rate control on the
per-pathway p-values, with a cutoff of FDR &lt; 0.05 reported in
practice. The BH procedure is valid under positive dependence,
which holds for the gene-set universe because parent-child and
overlap-induced correlations between pathway p-values are positive.

#warn[
  The conservative-but-valid status of BH under positive dependence
  is not a license to ignore the dependence structure. Two
  significantly-enriched GO terms that share most of their genes are
  one finding, not two. Redundancy-reduction tools (REVIGO,
  simplifyEnrichment) cluster significant terms by semantic
  similarity; modern enrichment workflows report cluster
  representatives, not raw tables of correlated hits.
]

=== GSEA: Using the Full Ranked List

ORA's biggest weakness is that it throws away the rank information.
A pathway whose every gene is just below the FDR cutoff produces no
ORA enrichment because the overlap count is zero. *Gene-Set
Enrichment Analysis* (Subramanian et al., *PNAS*, 2005) was
designed to fix that. Rank every tested gene by a statistic — most
commonly the signed log fold-change or the signed t-statistic from
the differential test — and ask whether the genes of a pathway
cluster toward one end of the ranked list.

The algorithm walks down the ranked list, incrementing a running
*enrichment score* (ES) each time it encounters a gene in the
target set and decrementing it each time it encounters one outside
the set. The increments and decrements are weighted so the random
walk has zero expected drift under the null. The reported ES is the
maximum-magnitude excursion of the running statistic. The
statistical significance is computed by permutation: shuffle the
gene labels and recompute the ES; the empirical p-value is the
fraction of permutations whose ES exceeds the observed one.

#figure(
  image("../../diagrams/lecture-22/10-gsea-running.svg", width: 95%),
  caption: [
    The GSEA running enrichment statistic. The curve walks from left
    to right along the ranked list; it climbs at gene-set hits and
    falls at non-hits. The maximum-magnitude excursion is the
    enrichment score, and the genes contributing to it are the
    leading-edge subset — the candidates for downstream
    interpretation.
  ],
) <fig:gsea-running>

GSEA's contribution was not the underlying statistic — the
maximum-deviation running-sum trick is the *Kolmogorov–Smirnov
two-sample statistic*, a textbook tool in classical statistics. The
contribution was the visualisation (the running-score plot is
instantly readable), the permutation-based significance test that
is robust to gene–gene correlation, and the recognition that the
right input is the ranked statistic on every gene, not the
thresholded list.

#note[
  GSEA is mathematically equivalent to a signed Mann–Whitney test
  on the rank of pathway genes against the rank of non-pathway
  genes, with the rank-walk being a convenient bookkeeping device.
  An EE student already knows the underlying tool from non-parametric
  statistics courses; the bioinformatics packaging adds the
  visualisation and the gene-set machinery on top.
]

The signed-versus-unsigned distinction matters for interpretation. A
*signed* ranking statistic (log fold-change) detects _directional_
enrichment — pathway up or pathway down. An *unsigned* statistic
(absolute log fold-change or a $-log_10 p$) detects _deregulation_
in either direction. The choice depends on the question: "is
glycolysis up in tumour versus normal" wants signed; "is glycolysis
perturbed by drug X" wants unsigned.

Modern GSEA variants address speed and inter-gene correlation.
*fgsea* (Sergushichev, 2016) provides a fast preranked
implementation that is the workhorse for routine RNA-seq follow-up.
*camera* (Wu & Smyth, 2012) applies an explicit variance-inflation
correction for inter-gene correlation within pathways, which makes
the calibration more conservative on real data. *roast* and
*mroast* (Wu et al., 2010) use rotation-based testing to avoid
permutation entirely. For typical workflows the recipe is: fgsea
plus MSigDB Hallmark plus KEGG plus Reactome, FDR &lt; 0.05.


== Network Propagation <sec:propagation>

Pathway enrichment treats the gene list as a set; it sees no edges.
The next move is to embed the list in a network and let the
information about which genes interact with which inform the
analysis. The technique is *network propagation*: start with a
signal on a small set of seed nodes (the DEGs, the GWAS hits, the
known drug targets) and let it spread across the network. After
spreading, every node has a score that reflects its proximity to
the seeds, and the score ranks novel candidates by their network
distance to known biology.

=== #idx("random walk with restart")Random Walk With Restart

The most common propagation algorithm is *random walk with restart*
(#idx("RWR")RWR). Define the column-normalised adjacency $W = A D^(-1)$, so
$W_(i j)$ is the probability of stepping from $j$ to $i$. Pick a
restart probability $r in (0, 1)$, typically around $0.5$. Initialise
$bold(p)^((0))$ to a seed vector — $1$ at known-disease genes or
known drug targets, $0$ elsewhere, normalised so the entries sum to
one. The propagation iteration is

$ bold(p)^((t + 1)) = (1 - r) W bold(p)^((t)) + r bold(p)^((0)). $

At each step, with probability $1 - r$ the walker moves to a random
neighbour weighted by the column-normalised adjacency; with
probability $r$ it teleports back to the seed distribution. Iterate
for about fifty steps and the distribution converges. The
stationary distribution $bold(p)^(infinity)$ satisfies the linear
system $bold(p)^(infinity) = (1 - r) W bold(p)^(infinity) + r
bold(p)^((0))$, which has the closed-form solution

$ bold(p)^(infinity) = r (I - (1 - r) W)^(-1) bold(p)^((0)). $

#figure(
  image("../figures/ch22/f3-rwr-convergence.svg", width: 95%),
  caption: [
    RWR convergence and the restart-probability knob. Iterations 0,
    5, 20, and the stationary distribution (left) show mass
    redistributing from a single seed across the graph. Sweeping
    $r$ from 0.1 to 0.9 (right) controls how concentrated the
    stationary distribution stays around the seed: high $r$ keeps
    mass near the seed, low $r$ spreads it broadly. The literature
    default is $r = 0.5$ because it balances the two regimes; the
    answer is rarely sensitive to the exact choice within $[0.3,
    0.7]$.
  ],
) <fig:rwr-conv>

#figure(
  image("../../diagrams/lecture-22/04-network-propagation.svg", width: 95%),
  caption: [
    Network propagation as low-pass filtering. The signal starts as
    impulses at the seed nodes and smooths across the graph; the
    steady-state distribution decays with graph distance from the
    seeds. The right panel makes the analogy to Gaussian smoothing
    on a regular grid explicit — same low-pass logic, different
    underlying geometry.
  ],
) <fig:propagation>

#note[
  RWR with a uniform restart vector is *PageRank*: the algorithm
  Brin and Page (1998) introduced at Stanford to rank web pages by
  topological importance and used to launch Google. Personalised
  PageRank puts non-zero restart probability only on a small set of
  user-of-interest pages; biological RWR puts it only on seed
  genes. Same algebra, different seed. The bioinformatics literature
  rediscovered the algorithm independently around 2008 (PRINCE,
  Vanunu et al., 2010) and the convergence was made explicit in
  Cowen et al. (2017).
]

=== Heat Diffusion and the Spectral View

A continuous-time analogue replaces the discrete iteration with the
heat-diffusion equation on the graph,

$ frac(d bold(p), d t) = -L bold(p), $

with $L$ the (normalised) Laplacian. The solution is $bold(p)(t) =
e^(-L t) bold(p)(0)$, where $e^(-L t)$ is the *heat kernel*. In the
Laplacian eigenbasis,

$ bold(p)(t) = sum_i e^(-lambda_i t) angle.l bold(p)(0), bold(v)_i
angle.r bold(v)_i. $

The high-eigenvalue components decay quickly; only the low-eigenvalue
components persist. This is the *low-pass filter* statement of
network propagation: the high frequencies of the seed signal are
attenuated, and the smoothed output is dominated by the
slowly-varying components — exactly the smooth-across-edges
patterns that biological pathways correspond to.

The closed-form RWR solution and the heat kernel are two
parametrisations of the same low-pass filter on the graph spectrum.
RWR uses the spectral filter $r / (1 - (1 - r) mu_i)$ where
$mu_i$ are the eigenvalues of $W$; heat diffusion uses
$exp(-lambda_i t)$. The qualitative behaviour is identical —
neighbours of seeds accumulate score, distant nodes do not — and
the practical recipe (run for 50 iterations or invert the matrix)
is the same.

The HotNet2 method (Leiserson et al., *Nature Genetics*, 2015) uses
the heat-diffusion variant to discover "hot subnetworks" of cancer
mutations: seed the diffusion at mutated genes, run to a chosen
time scale, threshold the diffused score, and extract the connected
components of the high-score subgraph. The resulting modules
recover both well-known pathways (PI3K/AKT in glioblastoma,
RTK/RAS/MAPK in melanoma) and rare composite-driver modules where
no single gene is significantly mutated but the module as a whole
is.

PRINCE (Vanunu et al., *PLoS Computational Biology*, 2010) is the
canonical RWR application to disease-gene prioritisation. Seed the
walk at known disease genes; rank candidate genes in a chromosomal
region by their RWR scores; the top-ranked candidate is the most
likely culprit. The same recipe powers NetWAS for GWAS
interpretation: seed at GWAS-hit nearest genes; spread; the
high-score neighbours are the implicated mechanism even when they
themselves are not significant.

#figure(
  image("../../diagrams/lecture-22/09-drug-target-rwr.svg", width: 95%),
  caption: [
    RWR in a heterogeneous drug–target–disease network. Heat
    propagates from a drug seed through known targets to PPI
    neighbours and to drug-associated diseases; novel candidates
    rank high in the steady state. The same algorithm with
    different seeds answers different questions.
  ],
) <fig:drug-rwr>

=== What Propagation Is Doing Formally

The smoothing interpretation is the right one to internalise. In
the Fourier basis of $L$, propagation multiplies the seed signal's
coefficients by a monotonically decreasing function of the
eigenvalue. High-frequency components — values that change sign
between neighbours — are damped; low-frequency components — values
that vary smoothly across edges — are preserved. The result is the
graph analogue of Gaussian smoothing on an image grid, with the
graph Laplacian playing the role of the (negative) image Laplacian
in the heat equation.

The biological prior the smoothing encodes is _guilt by
association_: neighbours in the PPI network are more likely to share
function than random pairs of proteins. Network propagation
operationalises that prior. A gene whose own measurement did not
reach significance but whose neighbours did acquires a non-trivial
score; a gene with no relevant neighbours stays at zero. The signal
is amplified along the smooth low-frequency directions of the graph
spectrum and damped along the rough high-frequency directions.


== Module Detection <sec:modules>

The third major class of network operations partitions the graph
into communities — densely-connected subgraphs sparsely connected
to each other. In biology, communities correspond to functional
units: a protein complex, a signalling pathway, a single-cell
cluster. Module detection turns a structurally uniform graph into a
labelled one.

=== Modularity and Its Spectral Twin

The standard quality function for community partitions is *Newman
modularity* (Newman, 2004; Newman, *PNAS*, 2006). Given a partition
$\{C_1, dots, C_K\}$ of the nodes, the modularity is

$ Q = frac(1, 2 m) sum_(i, j) [A_(i j) - frac(k_i k_j, 2 m)]
delta(c_i, c_j), $

where $m = |E|$, $k_i$ is the degree of node $i$, and $delta(c_i,
c_j) = 1$ if nodes $i$ and $j$ are in the same community and 0
otherwise. The interpretation is that $Q$ compares the within-community
edge fraction to what you would expect under a
*degree-preserving random graph* — the configuration model.
High $Q$ means the partition has more within-community edges than
expected by chance given the node degrees.

Maximising $Q$ is NP-hard. Two practical algorithms have become the
standard.

*#idx("Louvain")Louvain* (Blondel et al., *Journal of Statistical Mechanics*, 2008)
is a greedy local-move heuristic. Initialise each node as its own
community. For each node in turn, move it to whichever neighbouring
community produces the largest increase in $Q$, breaking ties
arbitrarily. Iterate until no move increases $Q$. Then aggregate:
each community becomes a super-node, the inter-community edges
become weighted super-edges, and the algorithm restarts on the
aggregated graph. Repeat until the aggregation step finds no further
improvement. The whole procedure runs in $O(n log n)$ time in
practice and scales to millions of nodes.

*#idx("Leiden")Leiden* (Traag, Waltman, and van Eck, *Scientific Reports*, 2019)
improves Louvain by fixing a documented pathology: Louvain
occasionally produces disconnected communities, where a community's
member nodes do not form a connected subgraph. Leiden inserts a
refinement step between the local-move and aggregation phases that
checks community connectivity and guarantees, by construction, that
every community in the output is internally connected. The
single-cell community has converged on Leiden as the default
clustering algorithm; the same Leiden that runs inside #idx("Scanpy")scanpy and
#idx("Seurat")Seurat for cell-type discovery is the algorithm here.

=== Spectral Clustering

The alternative is *spectral clustering* (Shi and Malik, 1997; Ng,
Jordan, and Weiss, 2001), and it is mathematically related to
modularity maximisation through the graph Laplacian.

The recipe is short. Compute the Laplacian $L = D - A$. Find the
$k$ smallest non-zero eigenvalues with their eigenvectors $bold(v)_2,
dots, bold(v)_(k + 1)$. Stack them into an $n times k$ matrix $V$.
Treat each row of $V$ as a $k$-dimensional embedding of the
corresponding node, and cluster the rows with k-means.

What the algorithm is doing is using the low-frequency directions of
the graph spectrum as coordinates. The Fiedler vector $bold(v)_2$ —
the smallest non-zero eigenvector — already encodes the best
two-way cut of the graph in the sense of minimising the
normalised-cut criterion. Adding $bold(v)_3, dots$ extends the
embedding to more communities. K-means in this low-dimensional
spectral space recovers community labels that closely match what
modularity maximisation finds, because both algorithms are
operating on the same low-eigenvalue subspace.

#figure(
  image("../../diagrams/lecture-22/05-spectral-clustering.svg", width: 95%),
  caption: [
    Spectral clustering. Project the network onto the smallest
    non-zero eigenvectors of the Laplacian; communities separate
    cleanly in this low-frequency embedding. K-means in the
    eigenspace recovers the community assignment.
  ],
) <fig:spectral>

#note[
  The connection between modularity and spectral clustering was
  made rigorous in the early 2010s. Both algorithms produce similar
  partitions on biological networks in practice; the difference is
  speed (Louvain and Leiden scale better to millions of nodes) and
  theoretical guarantees (spectral clustering has a tighter
  cut-quality bound). For a working bioinformatics pipeline, Leiden
  is the right default; spectral clustering is the right baseline
  for theoretical comparison.
]

The resolution-limit question is worth flagging. Modularity has a
known _resolution limit_ (Fortunato and Barthélemy, 2007): below a
characteristic size set by the total number of edges, communities
become invisible to $Q$-maximisation no matter how dense they are.
The fix is the _resolution parameter_ in modern Leiden
implementations: it scales the expected-edge baseline and lets the
analyst dial the granularity of the partition between many small
modules and few large ones. Single-cell scRNA-seq analyses often
sweep the resolution parameter and pick the value at which the
recovered clusters match the expected biology.


== Applications: Drugs, Disease, Cells <sec:applications>

The four operations of the previous sections — pathway enrichment,
propagation, module detection, and the underlying graph algebra —
combine into a small number of recurring analyses that practising
bioinformaticians run repeatedly. The shape of each is the same:
seeds in, network operations applied, ranked candidates out.

=== Drug-Target Prediction

Build a heterogeneous network of drugs, protein targets, diseases,
and side effects. Edges encode drug-target binding (DrugBank,
#idx("ChEMBL")ChEMBL), protein–protein interactions (STRING), drug-disease
indications (Open Targets), and protein-disease associations
(DisGeNET). Seed the propagation at a known drug; let the heat
spread through known targets to their PPI neighbours and onward to
drug-associated diseases; rank candidate new targets by their
stationary score.

The pattern is *guilt by association* operationalised at scale.
The clinical applications fall into two categories. *Drug
repositioning* finds new uses for existing drugs: a drug approved
for one indication that propagates strongly to genes implicated in
another disease is a candidate for the second. The most-celebrated
example is sildenafil, originally developed for angina and
repositioned for erectile dysfunction; modern network-medicine
methods aim to make such discoveries systematic rather than serendipitous.
*Off-target prediction* anticipates safety issues before clinical
trials: a drug whose RWR distribution accumulates on proteins
unrelated to its intended target may have side effects the trial
will eventually surface, and the network analysis can flag them
months earlier than a phenotypic screen.

Methods in this space include *NeoDTI* (a graph-neural-network
extension of guilt-by-association), *DrugBank-based RWR*, and the
network-proximity framework of Cheng et al. (*Nature
Communications*, 2018), which formalises the proximity between
drug-target sets and disease-gene sets in the human interactome and
shows that small proximity predicts drug efficacy across hundreds
of drug-disease pairs.

=== Disease-Gene Prioritisation

A clinical-genetics scenario. A child presents with an undiagnosed
syndrome. Trio whole-exome sequencing identifies a chromosomal
region in which several candidate genes carry rare variants. Which
is the disease-causing gene?

The network-based approach takes a curated set of seed disease
genes for similar syndromes — from OMIM, DisGeNET, or the literature
— and propagates the signal across the STRING (or BioGRID) PPI
network. Among the candidates in the patient's region, the gene
with the highest propagation score is the most likely culprit, and
becomes the priority for variant interpretation and clinical
validation.

PRINCE (Vanunu et al., 2010) is the canonical implementation;
DOMINO (Levi, 2021) refines the prioritisation by enforcing module
structure on the result; GeneMANIA (Warde-Farley, 2010) provides an
interactive web front-end. Routine in modern clinical genomics for
novel disease-gene discovery, with growing integration into
hospital diagnostic pipelines.

#figure(
  image("../../diagrams/lecture-22/11-disease-module.svg", width: 95%),
  caption: [
    A disease-gene module on the PPI. Known disease genes form a
    connected module via high-confidence PPI edges; novel
    candidates that rank high by network proximity can be
    prioritised for clinical validation. The pathway-enrichment
    annotation of the module (right panel) supplies the mechanistic
    interpretation.
  ],
) <fig:disease-module>

The *disease-module hypothesis* (Barabási, Menche, Goh, and others,
2011 onward) is the conceptual underpinning. Empirical work on the
human interactome shows that genes associated with the same disease
tend to cluster in a localised neighbourhood of the network — the
#idx("disease module")disease module — rather than scattering uniformly. Diseases whose
modules overlap in the network tend to share clinical features and
respond to overlapping drug interventions. The module hypothesis
turns clinical comorbidity and drug repositioning into a network
question, and the answers it gives have been progressively
validated against clinical-outcome data over the last decade.

=== Pathway-Context Interpretation of DEGs

When the gene list is the output of a differential-expression
analysis, the standard follow-up is pathway enrichment, but
network methods take the analysis further. Project the DEGs onto
the PPI network; run propagation to find connected modules among
DEGs; identify "core" hubs that connect multiple DEGs and that may
explain the cascade; hypothesise a causal mechanism in which the
hub's perturbation propagates to its neighbours.

This is the workflow behind OmicsNet, NetworkAnalyst, and the
STRING-driven web interface that most experimentalists use without
realising it is running RWR under the hood. The output is a
publication-ready figure showing the DEG-induced subgraph annotated
by pathway membership, plus a short list of candidate driver hubs.
The interpretation is, as ever, an exercise in caution: hubs are
hubs in every analysis, and the next section returns to how to
correct for that.

=== Cell-Cell Communication

The most recent application is *cell–cell communication inference*
from #idx("single-cell RNA-seq")single-cell RNA-seq data. Single-cell experiments measure
expression at the level of individual cells; the cells partition
into types; and the question is which cell types are signalling to
which others via known #idx("ligand-receptor")ligand-receptor pairs.

#idx("CellChat")CellChat (Jin et al., *Nature Communications*, 2021) and
CellPhoneDB (Efremova et al., *Nature Protocols*, 2020) are the
standard tools. The recipe builds a directed network whose nodes
are cell-type pairs and whose edges are weighted by the joint
expression of curated ligand-receptor pairs from databases like
OmniPath. Integrative analysis identifies the dominant communication
pathways in a tissue, and the results overlay onto histology to
locate signalling neighbourhoods. Single-cell perturbation
experiments — CRISPR-Perturb-seq, drug screens — then validate the
inferred edges by knocking out the ligand and measuring the
downstream cascade in the receptor-expressing cells.

=== The Graph-Neural-Network Frontier

The 2024 frontier is *graph neural networks* (GNNs) replacing
classical RWR for some prediction tasks. Where RWR uses a fixed
propagation operator — the column-normalised adjacency or the heat
kernel — GNNs learn a task-specific node embedding via end-to-end
training on labelled examples. Variants include the graph
convolutional network of Kipf and Welling (2017), the
message-passing neural network of Gilmer et al. (2017), and the
graph #idx("attention")attention network of Veličković et al. (2018). Multi-modal
network medicine — PPI plus GRN plus drug-target plus GWAS in one
heterogeneous graph, with a single GNN learning unified embeddings —
is the active research area.

The honest reading is that GNNs do not yet uniformly replace RWR
for the canonical tasks. The Cowen et al. (2017) review concluded
that network propagation is _the_ universal amplifier of genetic
associations; the more recent benchmarks (Huang et al., 2022;
Mostafavi et al., 2023) show that GNNs match classical propagation
on most disease-gene-prioritisation tasks and beat it modestly on
heterogeneous-network tasks with large labelled training sets. The
inductive-bias discussion of Chapter 16 applies directly:
classical propagation embeds the smoothness prior explicitly;
GNNs have to learn it from data. When the training set is small —
the typical genomics setting — the explicit prior wins.


== Pitfalls and Practice <sec:pitfalls>

A working network-biology project carries a small set of recurring
design decisions, and the same handful of mistakes recurs in
published analyses. The list below is not exhaustive, but it
captures the failure modes the previous five sections of this
chapter are organised to prevent.

#figure(
  image("../figures/ch22/f4-pitfall-checklist.svg", width: 95%),
  caption: [
    A design-review checklist for a network-biology project. Five
    failure modes on the left, five diagnostic fixes on the right.
    Most published controversies in network medicine reduce to one
    of them.
  ],
) <fig:ch22-checklist>

*Hub bias.* The single most common failure mode. Hub proteins —
TP53, ACTB, EGFR — appear in every analysis as "significant"
because their high degree, high pathway-membership count, and
heavy literature curation conspire to make them score high under
any reasonable scheme. Reporting that TP53 is enriched in your
cancer DEGs is not biology; it is a statement about how much TP53
has been studied. The fix is a *degree-preserving null*. Compute
the score against permuted graphs that preserve every node's
degree, and report the z-score relative to the null rather than
the raw rank. The configuration-model permutation is computationally
cheap and changes top-ten lists dramatically.

*False edges.* High-throughput interaction screens have
false-positive rates in the 30 to 50 percent range, even at the
manufacturer's recommended cutoff. STRING's "predicted-only"
edges can be entirely spurious in poorly-studied gene families
where the prediction is driven by genomic-context heuristics rather
than physical evidence. The fix is to filter to high-confidence
experimental edges before running any analysis: STRING combined
score above 0.7, or restrict to BioGRID curated. The right
sensitivity analysis is to recompute the analysis under both a
permissive and a strict threshold; if the top-ten list changes
qualitatively, the noise is doing the work.

*Tissue mismatch.* PPI databases are species-aggregate, pooled
over every cell type and condition. Neural-specific signalling
interactions show up as edges in the #idx("MUSCLE")muscle PPI, where they may
have no biological relevance to a cardiomyopathy study. The fix
is to filter edges by the tissue expression of both endpoints:
use GTEx or the Human Protein Atlas to retrieve the
tissue-specific expression vector and keep only edges where both
proteins are expressed above a threshold. Tissue-specific PPI
networks (HumanNet-XC, GIANT) integrate this filtering at the
database level.

*Wrong gene-universe.* The ORA universe $|U|$ must be the set of
genes that were _tested_, not the set of all human genes. Inflating
$U$ shrinks the expected overlap and makes every gene set look
enriched. For RNA-seq, drop genes below the expression filter from
$U$ before running the test. For targeted panels — CRISPR screens,
exome panels — restrict $U$ to the panel. Almost every published
enrichment table this section is teaching you to read sceptically
contains this error somewhere.

*Overlapping pathways.* GO parent-and-child terms are heavily
correlated; KEGG super-pathways share genes; MSigDB collections
overlap. A significant parent plus a significant child is one
biological finding double-counted. Redundancy reduction — REVIGO
for GO, simplifyEnrichment for general pathway sets — clusters
significant terms by semantic similarity and reports cluster
representatives. Modern enrichment workflows use redundancy
reduction by default; check whether your pipeline does.

#figure(
  image("../../diagrams/lecture-22/12-hub-bias.svg", width: 95%),
  caption: [
    Hub-centrality bias. Hubs like TP53 (degree ~250) appear
    significant in every analysis by virtue of their high degree
    alone. The naive p-value (left) and the degree-corrected
    p-value (right) routinely disagree by orders of magnitude.
    Correct before interpreting.
  ],
) <fig:hub-bias>

=== The Standard Workflow

A canonical end-to-end network-analysis workflow has the following
shape, and most published analyses can be read as instances of it
with minor variation.

1. *Define seeds.* DEGs at FDR &lt; 0.05, GWAS-significant SNPs'
   nearest genes, tumour drivers from a cancer cohort, or known
   disease genes from a curated database.
2. *Pick a network.* STRING (broad), BioGRID (curated experimental),
   Reactome FI (pathway-based), or a tissue-specific subset.
3. *Filter the network.* High-confidence edges (STRING &gt; 0.7),
   tissue-expression filter, removal of multi-degree promiscuous
   hubs if appropriate.
4. *Run propagation or community detection.* RWR for seed-driven
   ranking; Leiden for unsupervised module discovery; HotNet2 for
   hot-subnetwork extraction from cancer mutation data.
5. *Annotate.* Pathway enrichment on each module via fgsea or
   clusterProfiler; map back to biological interpretation.
6. *Visualise.* Cytoscape for figures; networkx or igraph for
   programmatic analysis.

#figure(
  image("../../diagrams/lecture-22/06-workflow.svg", width: 95%),
  caption: [
    Standard network-analysis workflow. Modular: seeds in, biology
    out; each step has its own tool and runtime budget. The
    GWAS-driven variant in the side panel substitutes NetWAS or
    DOMINO for the propagation step.
  ],
) <fig:ch22-workflow>

The tools cluster around two ecosystems. *Cytoscape* is the
desktop standard for visualisation, with hundreds of plugins for
layout, statistics, and database integration. *Gephi* is its main
desktop alternative. For programmatic analysis, *networkx* (Python)
is the de-facto standard, with *igraph* (R / Python) the
performance-oriented alternative. The R *Bioconductor* stack —
clusterProfiler, fgsea, ReactomePA, STRINGdb — provides
domain-specific wrappers that the canonical RNA-seq follow-up runs
on. A typical published RNA-seq paper uses igraph or networkx for
analysis and Cytoscape for figure generation.

#warn[
  When you read a network-medicine paper that reports a striking
  novel drug target or disease-gene prioritisation, the questions
  to ask in order are: (1) was the universe correct for enrichment?
  (2) was hub bias corrected? (3) was the network tissue-filtered?
  (4) was the result robust to the edge-confidence threshold? Most
  unreplicated network-medicine findings trace back to one of
  these four. Reproducing the analysis under a tightened protocol
  is a standard sanity check; it frequently changes the headline
  result.
]


== Summary <sec:ch22-summary>

- Five biological network types — PPI, GRN, metabolic, signalling,
  and disease — share scale-free degree distributions, small-world
  path lengths, and modular structure with social and citation
  networks. The toolkit of network science transfers directly.
- The graph Laplacian $L = D - A$ is the spectral operator that
  organises this chapter. Its eigenvectors form the graph Fourier
  basis; low-eigenvalue eigenvectors are smooth across edges,
  high-eigenvalue eigenvectors are rough. Network propagation
  is a low-pass filter in this basis; spectral clustering uses
  the same basis as coordinates.
- Pathway enrichment comes in two flavours. ORA (hypergeometric or
  Fisher's exact) discretises the gene list at a significance
  threshold and counts overlap with a gene set; GSEA
  (Kolmogorov–Smirnov on a signed-rank statistic) uses the full
  ranked list and catches coherent shifts ORA misses.
- Network propagation — random walk with restart, heat diffusion
  — spreads signal from seed nodes across the network. The
  closed-form solution $bold(p)^infinity = r (I - (1 - r) W)^(-1)
  bold(p)^((0))$ is the same algebra as personalised PageRank.
  Used for disease-gene prioritisation (PRINCE), cancer-driver
  module discovery (HotNet2), GWAS interpretation (NetWAS), and
  drug repositioning.
- Module detection partitions the graph into communities.
  Modularity maximisation (Louvain, Leiden) and spectral
  clustering operate on the same low-eigenvalue subspace and
  produce similar partitions on biological networks. Leiden's
  connectivity guarantee makes it the default in single-cell
  analysis.
- The five recurring pitfalls — hub bias, false edges, tissue
  mismatch, wrong gene-universe, overlapping pathways — account
  for most network-medicine reproducibility failures. Run the
  design-review checklist before the analysis and re-run it
  before trusting another paper's headline result.
- The graph-neural-network frontier is active but does not yet
  uniformly replace classical RWR. The classical methods embed
  the smoothness prior explicitly; GNNs learn it from data. In
  the label-limited regime where most genomics applications
  live, the explicit prior usually wins.


== Exercises <sec:ch22-exercises>

#strong[1.] #emph[Hypergeometric ORA.] Given a universe of $|U| =
20{,}000$ tested genes, a candidate list of $|L| = 200$ DEGs, and
a target gene set of $|G| = 500$, compute the expected overlap
under the null and the hypergeometric p-value for an observed
overlap of $a = 20$. Verify your computation against
`scipy.stats.hypergeom.sf` and Fisher's exact test
(`scipy.stats.fisher_exact`). Then re-run with $|U| = 60{,}000$
(all human genes including non-coding) and explain in two
sentences why the p-value moves and which value is correct.

#strong[2.] #emph[GSEA implementation.] Implement the running
enrichment statistic of GSEA in fewer than 30 lines of Python.
Test it on a synthetic ranked list of 5,000 genes with a 100-gene
target set, where 50 of the target genes are placed in the top
500 ranks and 50 are scattered uniformly. Plot the running ES
and compute the maximum-magnitude excursion. Compare against the
fgsea package output on the same data; the two ES values should
agree to four decimal places.

#strong[3.] #emph[RWR on a small PPI.] Download a 100-node
subgraph of the human PPI from STRING (any seed gene of your
choice plus its two-hop neighbourhood). Implement RWR in
networkx; iterate from a single seed gene with $r = 0.5$ until
$norm(bold(p)^((t + 1)) - bold(p)^((t)))_1 < 10^(-6)$. Report
the top-10 ranked nodes. Compare against the closed-form
solution $bold(p)^infinity = r (I - (1 - r) W)^(-1)
bold(p)^((0))$ — the two should agree to numerical precision.

#strong[4.] #emph[Spectral clustering by hand.] Construct a
40-node synthetic graph with three planted communities of 12, 14,
and 14 nodes respectively, with intra-community edge probability
0.4 and inter-community edge probability 0.05. Compute the
Laplacian, extract the second and third smallest non-zero
eigenvectors, plot the 2-D embedding, and cluster with k-means
($k = 3$). Compute the normalised mutual information between
your inferred clusters and the planted ground truth. Then run
Leiden on the same graph (via `python-igraph`) and compare.

#strong[5.] #emph[Hub-bias correction.] Take any published GO
biological-process enrichment table from an RNA-seq paper of
your choice. Identify the three most-enriched terms and the
three highest-degree genes contributing to each. Compute the
degree-preserving null p-value for one of those terms: generate
1,000 random gene lists of the same size by sampling without
replacement from genes binned by degree, re-run the enrichment
test against each, and report the empirical p-value. How does it
compare with the paper's reported p-value?

#strong[6.] #emph[Universe matters.] Take a public RNA-seq
differential expression result (any GEO study). Run ORA against
KEGG with two universes: (a) all genes with an Ensembl ID, and
(b) the subset of genes that passed the expression filter and
entered the differential test. Compare the top-20 enriched
pathway lists. How many pathways disagree between the two? Which
of the two universes is correct?

#strong[7.] #emph[Disease-module hypothesis.] Pick a disease with
at least 20 known disease genes in DisGeNET (or pick a different
curated source). Extract the induced subgraph on the STRING PPI
and compute the disease module's _shortest-path distance_ to
1,000 random gene sets of the same size. Plot the distribution
of distances and locate where the disease module sits. The
disease-module hypothesis predicts the disease genes will sit at
a percentile far below the random distribution; verify whether
your disease confirms or violates it.

#strong[8.] #emph[(Open-ended.)] Find a published network-medicine
paper from 2022 — 2025 that reports a novel drug-repositioning
hypothesis. Run the five-question pitfall checklist
(@fig:ch22-checklist) against the paper. Identify the weakest of the
five links. Propose one robustness experiment that would change
your confidence in the headline hypothesis, and predict whether
it would strengthen or weaken the claim. Cite the paper.


== Further Reading <sec:ch22-further-reading>

- *Subramanian, A., Tamayo, P., Mootha, V. K., et al.* (2005).
  "Gene set enrichment analysis: a knowledge-based approach for
  interpreting genome-wide expression profiles." _PNAS_ 102:
  15545–15550. The GSEA paper. One of the most-cited
  bioinformatics papers ever.
- *Vanunu, O., Magger, O., Ruppin, E., et al.* (2010).
  "Associating genes and protein complexes with disease via
  network propagation." _PLoS Computational Biology_ 6: e1000641.
  The PRINCE paper — RWR for disease-gene prioritisation,
  rediscovering personalised PageRank in a biological context.
- *Cowen, L., Ideker, T., Raphael, B. J., and Sharan, R.* (2017).
  "Network propagation: a universal amplifier of genetic
  associations." _Nature Reviews Genetics_ 18: 551–562. The
  authoritative review of network propagation and its applications
  across the field.
- *Newman, M. E.* (2006). "Modularity and community structure in
  networks." _PNAS_ 103: 8577–8582. The modularity-spectral
  connection and the original modularity-maximisation framework.
- *Traag, V. A., Waltman, L., and van Eck, N. J.* (2019). "From
  Louvain to Leiden: guaranteeing well-connected communities."
  _Scientific Reports_ 9: 5233. The Leiden paper; read for the
  connectivity guarantee that Louvain lacks.
- *Barabási, A. L., Gulbahce, N., and Loscalzo, J.* (2011).
  "Network medicine: a network-based approach to human disease."
  _Nature Reviews Genetics_ 12: 56–68. The disease-module
  hypothesis and the foundational essay of network medicine.
- *Cheng, F., Desai, R. J., Handy, D. E., et al.* (2018). "Network-based
  approach to prediction and population-based validation of in
  silico drug repurposing." _Nature Communications_ 9: 2691. The
  network-proximity framework for drug repositioning, with
  population-scale clinical validation.
- *Shuman, D. I., Narang, S. K., Frossard, P., Ortega, A., and
  Vandergheynst, P.* (2013). "The emerging field of signal
  processing on graphs." _IEEE Signal Processing Magazine_ 30:
  83–98. The graph-signal-processing manifesto from the EE side
  of the same mathematics. Pair with Cowen et al. to see the
  two-cultures convergence.
