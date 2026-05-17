#import "../theme/book-theme.typ": *

= #idx("DNA")DNA Sequence Assembly: From Reads to Genomes <ch:assembly>

#matters[
  An assembled genome is the foundation on which every other analysis in
  this book quietly rests. Alignment in Chapter 2 only works against a
  reference somebody has already assembled. #idx("variant calling")Variant calling in Chapter 4
  is the difference between your reads and an assembled coordinate
  system. Expression analysis, regulatory annotation, the #idx("pangenome")pangenome —
  every one of them presumes a sequence to map onto. For most of the
  genomes you will encounter in practice — a soil microbe nobody has
  sequenced before, a tumour rearranged past recognition, a viral
  isolate from a clinical swab — that sequence does not yet exist. The
  only way to get one is to assemble it from the reads themselves.
]

The story of #idx("de novo assembly")de novo assembly is the story of one engineering decision
made well in 2001 and a long sequence of consequences. The decision was
to formulate the assembly problem on a different graph from the
intuitive one: not "a node per read" — that turns assembly into the
NP-hard Hamiltonian-path problem — but "a node per fixed-length
substring," which turns assembly into the linear-time Eulerian-path
problem. The same input data, the same desired output, a different
choice of state variable. Almost every short-read assembler since
#idx("Velvet")Velvet (Zerbino and Birney, 2008) is a refinement of that one
substitution.

This chapter walks you through the substitution and its consequences.
We start with the assembly pipeline as an end-to-end object: six stages
between raw reads and a polished sequence, each of which can fail
independently. We then climb into the two stages that carry most of the
algorithmic weight — #idx("error correction")error correction in the k-mer-frequency domain,
and graph traversal — before stepping back out to look at what really
limits assembly in practice (repeats) and how the field has learned to
measure assemblies honestly (#idx("N50")N50 and its sharper-edged cousins). Long
reads, mate-pair libraries, and #idx("Hi-C")Hi-C show up wherever they earn a
mention; they are not the protagonist of this chapter but they are the
reason any modern mammalian assembly exists at all.

The chapter is long because the topic earns it. Assembly is the place
where graph theory, signal processing, and biology all collide at once,
and the engineer who understands the collision is the engineer who can
spot a bad assembly without running QUAST.


== The Assembly Problem <sec:problem>

Chapter 2 closed with reads aligned against a reference genome. That
chapter took the reference for granted. For a few hundred well-studied
organisms — humans, mice, _E. coli_, _Drosophila_, _Arabidopsis_, and
maybe a hundred more — the assumption is valid: a high-quality reference
exists, and the right question is "where does this new read map?"
For everything else, the question is "what _is_ the sequence to begin
with?" There is no reference. There is no template. There is only a
tube of dissolved DNA, an instrument, and the #idx("FASTQ")FASTQ file the instrument
produces.

That reconstruction problem is _de novo_ assembly. Given millions of
short reads drawn — with errors — from an unknown genome, the assembler
returns its best guess at the original sequence: a set of long strings
called *contigs* that together cover as much of the genome as the data
allows, with explicit gaps wherever the data does not.

#note[
  The intuition is a jigsaw puzzle with no picture on the box. You have
  millions of pieces, each about 150 pixels wide, drawn from a picture
  three billion pixels long. The pieces overlap; each position in the
  picture appears in roughly thirty pieces. About one stroke in a
  hundred on each piece is wrong. Reconstruct the picture. That is
  what the assembler is doing — and the question of whether the
  resulting picture is correct, complete, and contiguous is what the
  rest of this chapter is about.
]

=== De novo versus resequencing

Two problems sound similar enough to be confused. *Resequencing*, the
subject of Chapter 2, takes a known reference and a sample expected to
be nearly identical to it; alignment is the bottleneck, the biology is
assumed, and the output is a list of differences (variants). This is
how 23andMe genotyping works, how clinical diagnostic panels work, and
how most of human population genetics works.

*De novo assembly* takes the same reads but no reference. Either none
exists or the sample diverges from the available reference far enough
that alignment would discard everything interesting — a cancer with
megabase-scale rearrangements, a new bacterial species, a plant that
has duplicated its entire genome twice in the last million years.
Reconstruction itself is the bottleneck. There is nothing to match
against.

The two problems converge in cost: resequencing is cheaper because
alignment is cheaper than reconstruction. They diverge in honesty:
resequencing inherits whatever biases were baked into the reference,
while a de novo assembly's biases are at least its own.

#figure(
  image("../../diagrams/lecture-03/01-denovo-vs-resequencing.svg", width: 95%),
  caption: [
    The same reads, two problems. Resequencing aligns to a known
    template; de novo assembly reconstructs the template from the reads.
    The pile of reads on the right has no horizontal axis to land on.
  ],
) <fig:denovo-vs-resequencing>

=== The pipeline

An assembly run is never a single program. It is a pipeline of six
stages, and the production assemblers you have heard of — #idx("SPAdes")SPAdes, Flye,
Canu, hifiasm — are orchestrators that run roughly the same sequence in
the same order:

1. *Input quality control.* Adapter trimming, length filtering,
   per-read quality trimming. Same tools as Chapter 1 — `fastp`,
   `trim_galore`. A pre-flight check, not the main event.
2. *Error correction.* Clean systematic read errors before they poison
   the graph. Tools: `BFC`, `Lighter`, SPAdes's built-in BayesHammer.
   Covered in @sec:errors.
3. *Graph construction.* Build a #idx("de Bruijn graph")de Bruijn graph (short reads) or an
   #idx("overlap graph")overlap graph (long reads). Covered in @sec:debruijn and @sec:overlap.
4. *Graph cleanup.* Remove tips, collapse bubbles, simplify or report
   tangles. Covered in @sec:topology.
5. *#idx("contig")Contig construction.* Walk the cleaned graph; emit one contig per
   non-branching path. Covered in @sec:topology.
6. *Scaffolding.* Use #idx("paired-end")paired-end, mate-pair, Hi-C, or optical-map
   information to link contigs across the gaps the reads could not
   bridge. Covered in @sec:scaffolding.

A seventh stage, *polishing*, is sometimes folded in: align the reads
back to the assembly and correct any residual single-base errors with
`Pilon` (short reads) or `Medaka`/`Racon` (long reads). It is optional
but standard for publication-quality assemblies.

#figure(
  image("../../diagrams/lecture-03/02-assembly-pipeline.svg", width: 100%),
  caption: [
    The assembly pipeline. Each stage can fail independently. Most
    post-mortems on a bad assembly are about figuring out which stage
    let a problem through.
  ],
) <fig:ch03-pipeline>

The flow is linear and the stages compose, but the same physical
assembler may iterate across stages — SPAdes, for example, builds graphs
at several values of _k_, merges them, and re-runs cleanup after the
merge. Treat the pipeline as a logical ordering, not a literal one.

=== Reads and read types

Not all reads are equally useful for assembly. Three axes matter:
length, accuracy, and whether the reads come in pairs whose insert size
is known.

- *Short accurate reads* (#idx("Illumina")Illumina, 150–300 bp, ≤0.5% per-base error).
  Cheap, deep, low error. Useless on repeats longer than the read.
  Dominant for bacterial and small-eukaryote assembly, and for the
  short-read half of hybrid assemblies.
- *Long noisy reads* (#idx("PacBio")PacBio #idx("CLR")CLR and #idx("Oxford Nanopore")Oxford Nanopore raw, 10–100 kb,
  5–15% raw error). Span most repeats. The noise corrects out under
  multi-fold consensus. Essential for mammalian-scale assembly.
- *Long accurate reads* (PacBio #idx("HiFi")HiFi, 10–25 kb, &lt;0.1% error). The
  best of both worlds, at a price premium. The current gold standard
  for de novo mammalian assembly; the #idx("T2T")T2T human genome was a HiFi-and-
  #idx("ONT")ONT project.
- *Paired-end and mate-pair reads.* Two reads from opposite ends of a
  single fragment of known approximate length. Paired-end inserts are
  300–600 bp; mate-pair (older terminology: "jumping libraries")
  inserts are 2–40 kb. The paired information bridges short repeats;
  the mate-pair information scaffolds.

#figure(
  image("../../diagrams/lecture-03/03-read-types.svg", width: 100%),
  caption: [
    Read types used for assembly, drawn on a common scale. The dashed
    outline on Oxford Nanopore indicates higher raw error; HiFi's solid
    outline indicates consensus-grade accuracy at long length.
  ],
) <fig:read-types>

The cheapest modern assemblies mix read types. Short reads bring per-
base accuracy; long reads bring span; mate-pairs and Hi-C bring long-
range order. The assembler that uses one type only is leaving
information on the table.


== #idx("coverage")Coverage: A Poisson Process with Uneven Bias <sec:coverage>

Before any assembly algorithm runs, the raw reads are characterised by
one summary number more than any other: *coverage*. Coverage at a
position is the number of reads that span that position. Averaged
across the genome it is the *average coverage* or *depth*:

$ "average coverage" = (N_"reads" dot |r|) / L_"genome" $

For a 5-megabase bacterial genome sequenced with 100,000 Illumina
150 bp reads, average coverage is $100{,}000 times 150 / 5{,}000{,}000 = 3 times$.
For a typical whole-genome shotgun run it is 30× to 60×.

The number is easy to compute and easy to misunderstand. Coverage is a
distribution, not a scalar.

#note[
  A bioinformatician used to thinking of "30× coverage" as a property
  of the assembly is closer in spirit to a power-supply engineer
  thinking of "12 V" as a property of a rail. Both abstractions are
  fine until something fails, at which point the engineer who has the
  full V-I curve at every node — or the full coverage histogram across
  every position — finds the failure and the other one does not.
]

=== The Poisson model

If reads were drawn uniformly at random across an _L_-base genome,
the depth at any single position would follow a Poisson distribution
with mean $lambda = N dot |r| / L$. The Poisson model is the same
distribution that photon-counting detectors, spike-sorting
neuroscientists, and twentieth-century telephone-switching engineers
have lived with for a century, and it predicts the coverage statistics
remarkably well over most of a real genome.

Two consequences of the Poisson model are worth carrying around as
back-of-envelope numbers.

First, at $lambda = 30$, the probability that a particular position
sees zero reads is $e^(-30) approx 9 times 10^(-14)$. For a 3-gigabase
human genome, the _expected_ number of fully uncovered positions under
the Poisson model is therefore $3 times 10^9 times 9 times 10^(-14)
approx 3 times 10^(-4)$ — effectively none. At $lambda = 5$, the same
expected count is $3 times 10^9 times e^(-5) approx 2 times 10^7$,
or 20 million uncovered bases. The Poisson model is steep, and the
practical implication is that doubling coverage from 5× to 30× moves
you out of "everything is missing" and into "everything is present."

Second, the standard deviation of the Poisson distribution is
$sqrt(lambda)$. At $lambda = 30$ the standard deviation is roughly 5.5,
so the per-position coverage in a Poisson world fluctuates within
something like ±2σ between 19 and 41. The graph cleanup heuristics in
@sec:topology lean on this fact: at 30× coverage you expect _true_
edges to carry weights in the 20–40 range, and a single-supporting-
read edge is overwhelmingly likely to be wrong.

#figure(
  image("../../diagrams/lecture-03/04-coverage-distribution.svg", width: 100%),
  caption: [
    Coverage as spatial pile-up (top) and statistical distribution
    (bottom). A Poisson curve fits the bulk of the histogram; the
    small spike at zero is the part the Poisson model misses — real
    dropouts.
  ],
) <fig:ch03-coverage>

=== Where the Poisson model breaks

The Poisson model assumes reads are drawn uniformly at random. They
are not. Three sources of non-uniformity matter in practice:

- *GC bias.* PCR amplification during library preparation is more
  efficient on moderate-GC fragments than on very high- or very low-GC
  regions. Extreme-GC regions end up under-represented by factors of
  two to ten. This is the dominant cause of short-read assembly gaps
  in real genomes.
- *Repeat collapse.* Reads from repetitive elements map equivocally to
  multiple positions. Some pipelines weight them fractionally, some
  drop them, some assign them randomly. The resulting coverage at
  repeat positions can look artificially low or artificially flat.
- *Strand bias and read-end effects.* The first and last few bases of
  a read have elevated error rates and are usually trimmed aggressively.
  Coverage piles up in the middle of reads and thins at their edges.

The cleanup for GC bias is a second PCR-free library preparation, or a
PCR-free chemistry such as #idx("10x Genomics")10X Genomics's linked reads or PCR-free
Illumina kits. The cleanup for repeats is long reads. There is no
cleanup for strand bias — you tolerate it and you report it.

#tip[
  Before you spend a day debugging an assembly, plot the coverage
  histogram. A clean Poisson-shaped peak around the expected depth
  means the data are usable. A bimodal histogram with a spike at zero
  means you have systematic dropouts and the assembly will be
  fragmented. A long right tail means you have a collapsed repeat
  inflating one or two contigs. Three commands — `samtools depth`,
  `awk` for the histogram, `gnuplot` or matplotlib for the picture —
  and you know more about your data than any single QC report will
  tell you.
]


== Error Correction in the k-mer-Frequency Domain <sec:errors>

Raw short reads have per-base error rates of roughly 0.1% to 1%. A
150 bp read therefore carries, on average, 0.15 to 1.5 errors. Long
raw reads from ONT have historically had 5% to 15% error rates. Feed
those reads directly into a de Bruijn graph and the graph becomes a
forest of spurious edges — each error creates a wrong _k_-mer that
branches off the true path and either terminates (if no other read
extends it) or rejoins the path further down (forming a bubble). Either
way, the graph balloons and the downstream traversal slows to a crawl.

Every production assembler corrects the reads before graph
construction. The dominant technique is *#idx("k-mer spectrum")#idx("k-mer")k-mer spectrum analysis*:
treat each read as a window over the genome, count every _k_-mer in
every read, and use the resulting frequency distribution to separate
true _k_-mers from errors.

=== The k-mer spectrum

Pick a value of _k_ — typically 21 or 31 for Illumina data. For each
read, slide a window of length _k_ and count how many times each
distinct _k_-mer appears across all reads.

The histogram of those counts has two recognisable populations.

First, a *true-k-mer peak*. A _k_-mer that occurs once in the genome
appears in about $lambda$ reads on average (where $lambda$ is the
coverage), with Poisson-shaped spread around that mean. At 30×
coverage this peak sits near 30. A _k_-mer that occurs twice in the
genome appears around 60. A _k_-mer that occurs $n$ times occurs at
$n lambda$. The peaks at multiples of $lambda$ correspond to repeat
copy numbers — a useful side-benefit of the spectrum.

Second, an *error-k-mer spike at frequency 1*. Each sequencing error
in a read produces (up to) _k_ distinct _k_-mers, each of which is
almost certainly unique (the probability that the _same_ random error
produces the _same_ erroneous _k_-mer in two different reads is
vanishingly small). So errors pile up at frequency 1, with a long thin
tail at frequencies 2 and 3 that comes from very rare repeated errors
or from heterozygous variants at low coverage.

#figure(
  image("../../diagrams/lecture-03/05-kmer-spectrum.svg", width: 100%),
  caption: [
    The k-mer spectrum at 30× coverage. The error peak at frequency 1
    is tall and narrow; the true-k-mer peak sits around the coverage
    depth with Poisson-shaped width. The threshold line separates them.
  ],
) <fig:kmer-spectrum>

The threshold is the heart of the trick. _k_-mers below it (say, count
&lt; 3 for moderate coverage) are presumed errors; _k_-mers above it
are presumed true. Setting the threshold is the same problem as setting
a detection threshold on a #idx("matched filter")matched filter — you pick the false-alarm
versus miss-detection tradeoff that matches your downstream noise
tolerance. Modern correctors do this Bayesianly, fitting both
populations and choosing the maximum-likelihood threshold.

=== The correction step

Once the threshold is in place, correction is a per-read walk. For each
position in the read, look at the _k_-mer ending at that position. If
the _k_-mer has high count, leave the base alone. If it has low count
— and especially if it sits inside a window of consecutive low-count
_k_-mers, which is the signature of a single base error — try the
three alternative bases at the suspect position and ask which (if any)
restores all the affected _k_-mers to the high-count peak. If exactly
one substitution does, accept it. If two do, the position is
ambiguous; leave it for the graph stage.

The geometry is worth pausing on. A single base error in the middle of
a read corrupts exactly _k_ consecutive _k_-mers (the ones whose window
spans the erroneous position). The error correction algorithm exploits
that signature: a long run of low-count _k_-mers is much more likely
to come from a real sequencing error than from a stretch of genuinely
rare _k_-mers, because rare-but-real _k_-mers do not show up in
contiguous runs of length _k_.

#figure(
  image("../figures/ch03/f1-kmer-correction-walkthrough.svg", width: 100%),
  caption: [
    A worked example of single-base correction at k = 9. The error
    base corrupts three consecutive 9-mers (the position-8 window).
    Of the four candidate substitutions, only one (C→G) restores all
    three to the high-count peak; the corrector accepts it and moves
    on.
  ],
) <fig:kmer-correction>

#note[
  Error correction is signal-noise separation in the k-mer-frequency
  domain. True _k_-mers sit in a Poisson peak around the coverage
  depth; erroneous _k_-mers sit in a long tail near frequency 1. The
  two populations are well separated at 30× coverage and overlap
  catastrophically at 5×, which is exactly why low-coverage assembly
  is hard — the corrector cannot tell signal from noise at the
  histogram level.
]

Modern assemblers (SPAdes, Canu, Flye) run a Bayesian corrector like
BayesHammer automatically. You do not usually run it as a separate
step; you benefit from it implicitly. The choice of $k$ is the only
parameter you might want to tune: smaller $k$ correct more aggressively
but lose discrimination at repeats; larger $k$ are conservative but
require deeper coverage to populate the true-k-mer peak.


== A Detour through Graph Theory <sec:graph-theory>

Before de Bruijn graphs, a minimum viable vocabulary.

A *graph* is a set of *nodes* connected by *edges*. Edges can be
*directed* ($A arrow.r B$) or *undirected* ($A - B$). Graphs can have
*cycles* (paths that return to their starting node) or be *acyclic*.
Nodes can have *labels*; edges can have *weights*.

The *degree* of a node is the number of edges touching it. For directed
graphs you have *in-degree* and *out-degree* separately. Most graph
algorithms look at degree as the first signal — a node with in-degree
one and out-degree one is on a chain; a node with multiple in-edges or
out-edges is a branch point.

Two kinds of path show up constantly. An *#idx("Eulerian path")Eulerian path* visits every
edge exactly once. A *Hamiltonian path* visits every node exactly
once. The two sound similar and the difference is fatal.

Euler proved in 1735, studying whether one could walk all seven bridges
of Königsberg without crossing any of them twice, that an Eulerian
path exists in a connected graph if and only if every node has even
degree, or exactly two nodes have odd degree (the start and end). For
directed graphs the condition becomes: in-degree equals out-degree
everywhere, or exactly two nodes are imbalanced by one (one source,
one sink). When the path exists, *#idx("Hierholzer")Hierholzer's algorithm* constructs
it in $O(|E|)$ time: walk until you get stuck, then for any remaining
unvisited edges, find an intermediate node on your current path that
still has out-degree, and splice a new cycle in there.

Hamiltonian-path detection, by contrast, is NP-complete. There is no
known polynomial-time algorithm and there is good reason to believe
none exists. The gap between the two problems is the gap between
"runs in seconds on a bacterial genome" and "infeasible for any genome
of practical size."

#note[
  Graphs are bookkeeping for "A connects to B." In a circuit, $A$ and
  $B$ are nodes and the edge is a wire with some impedance. In a
  finite-state machine the nodes are states and the edges are
  transitions. In a #idx("Viterbi")Viterbi decoder the trellis is a graph and the
  best decoded sequence is a path through it. In short-read assembly
  the nodes are short substrings and the edges are the single-letter
  shifts that take one substring to the next. The algorithms transfer.
]


== De Bruijn Graphs <sec:debruijn>

Here is the substitution that made modern short-read assembly
feasible.

Given a set of reads and a parameter $k$, the *de Bruijn graph* has:

- one *node* per distinct $(k - 1)$-mer that appears as a prefix or
  suffix of some _k_-mer in the reads;
- one *directed edge* per distinct $k$-mer in the reads, connecting
  the node labelled with its first $(k-1)$ characters to the node
  labelled with its last $(k-1)$ characters.

That is the whole definition. Two lines. Implementations vary in how
they store it — succinct hashing for memory efficiency, Bloom filters
for approximate construction, FM-indices for compressed lookup — but
the concept is this simple.

#figure(
  image("../../diagrams/lecture-03/06-debruijn-graph.svg", width: 100%),
  caption: [
    Constructing a de Bruijn graph from three short reads at k = 5.
    One node per distinct (k−1)-mer, one edge per distinct k-mer.
    Reading the genome is walking the graph.
  ],
) <fig:debruijn>

The remarkable property is what happens next. To recover the original
genome, walk an Eulerian path through the graph — every edge exactly
once. Write the label of the first node, then append the last
character of each subsequent node. The reconstructed #idx("STRING")string has length
$"(number of edges)" + (k - 1)$ and, in the absence of errors and
repeats, is the original genome.

#note[
  A de Bruijn graph is the state-transition diagram of a shift register
  whose state is the current $(k-1)$-mer. The edges are exactly the
  shift transitions the register performs as it processes a string,
  which makes assembly the problem of recovering a string from
  observations of its state transitions — the same problem a Viterbi
  decoder solves against a convolutional-code trellis. Same algebraic
  structure; the cost function and the target representation are the
  only things that differ.
]

=== A brief history of the idea

The graph is named after *Nicolaas de Bruijn*, the Dutch mathematician
who in 1946 studied universal cyclic sequences with special
combinatorial properties — sequences of length $4^n$ over a four-letter
alphabet in which every length-$n$ window appears exactly once. The
"de Bruijn graph" was the combinatorial object he used to enumerate
them. The graph predates DNA sequencing by three decades; it predates
the Watson-Crick model by seven years.

The application to assembly came in two waves. *Idury and Waterman*
applied the graph to fragment assembly in 1995, in a paper that was
ahead of the hardware: the genomes were small, the sequencing data
were small, and the practical advantage over overlap-layout-consensus
was not yet decisive. *Pavel Pevzner* and colleagues at UCSD in 2001
reformulated the assembly problem on the de Bruijn graph in a way that
explicitly traded the Hamiltonian formulation for the Eulerian one,
and proved that the resulting algorithm scaled to genomes that the
overlap approach could not handle. The Eulerian reformulation was the
intellectual centerpiece — once you had it, the rest was implementation.

Implementation arrived in 2008 with *Velvet*, by Daniel Zerbino and
Ewan Birney at the European Bioinformatics Institute. Velvet was the
first de Bruijn-graph assembler that ran on a desktop computer
overnight on a bacterial genome. It set the template for everything
since: build the graph, walk it, clean it, walk it again, emit contigs.
Within four years SPAdes (Bankevich, Nurk, Antipov, Pevzner and
others, 2012) had refined the recipe with iterative multi-_k_ graphs,
Bayesian error correction, and small-genome polishing, and by 2015 the
de Bruijn graph was the universal short-read substrate. Velvet, ABySS,
SOAPdenovo, Minia, MEGAHIT, and SPAdes are all variants of the same
basic object.

=== Eulerian traversal in practice

Hierholzer's algorithm — the linear-time Eulerian-path construction —
is conceptually three lines. Start at the source. Walk edges. When you
get stuck (in a connected graph that has an Eulerian path, this can
only happen at the sink), find an intermediate node on your current
path that still has unvisited out-edges, and splice a sub-tour in
there.

#figure(
  image("../figures/ch03/f2-eulerian-traversal.svg", width: 100%),
  caption: [
    Hierholzer's algorithm on a small de Bruijn graph. The first walk
    runs from source to sink and stalls with edges still unvisited.
    Splicing in a missed cycle at one of the intermediate nodes
    consumes the rest of the graph. The bottom row reconstructs the
    underlying sequence from the path.
  ],
) <fig:eulerian>

For a bacterial genome at 30× coverage with $k = 25$, the graph has
roughly $1.5 times 10^8$ edges and Hierholzer's algorithm finishes in
a few seconds on a modern laptop. For a mammalian genome at the same
coverage you are looking at billions of edges, but the algorithm is
still linear and still runs in minutes once the graph is built. The
expensive step in modern assemblers is graph _construction_, not
graph _traversal_.

#note[
  The 1735 result is the reason assembly is polynomial. If we still
  formulated assembly as a Hamiltonian-path problem over reads-as-nodes
  — which is what overlap-layout-consensus does — we would be paying
  NP-hard costs to assemble. Pre-2001 assemblers did pay those costs,
  and the resulting runtimes were a major obstacle to the first
  human-genome efforts. The de Bruijn shift is, in retrospect, the
  single most important algorithmic move in the history of the field.
]


== Overlap Graphs and Long Reads <sec:overlap>

The older, more intuitive way to formalise assembly is the *overlap
graph*: one node per read, one edge from read $a$ to read $b$ if the
suffix of $a$ matches the prefix of $b$ for at least some threshold
length $w$. To assemble, find a Hamiltonian path through the graph,
concatenate the reads along the path with overlaps merged, and you
have the genome.

This is the *Overlap-Layout-Consensus (OLC)* approach. *Celera
Genomics* used it to assemble the human genome in 2000-2001 (the
landmark paper is Myers, Sutton, Delcher, and others, 2000), and it is
still the substrate of every modern long-read assembler — *Canu*,
*Flye*, *miniasm*, *hifiasm*, *HiCanu*. The reason long-read assemblers
have not switched to de Bruijn graphs is that, with reads ten or fifty
or a hundred kilobases long, there is no redundancy for the de Bruijn
machinery to exploit. Picking $k = 1000$ on long reads produces a
graph where almost every $k$-mer is unique; there is nothing to gain
by paying the substitution cost.

#figure(
  image("../../diagrams/lecture-03/07-overlap-graph.svg", width: 100%),
  caption: [
    An overlap graph. Each read is a node; each suffix-prefix overlap
    above the threshold is a directed edge. Assembly reduces to finding
    a Hamiltonian path — NP-hard in general but tractable when the
    graph is nearly linear, as it is for clean long-read data.
  ],
) <fig:overlap>

Hamiltonian-path search is NP-hard in the worst case, but the overlap
graphs that arise from real long-read data are almost-linear: each
read has at most one or two high-confidence overlap neighbours on each
side, and a greedy traversal does the right thing. When the graph
diverges — at repeats, at structural variants, at homozygous-to-
heterozygous transitions in a diploid — the assembler reports the
ambiguity rather than guessing. The honesty is built in.

=== Choosing the graph

The choice is almost never free. Read length forces it.

#figure(
  align(center)[
    #table(
      columns: (auto, auto, auto),
      stroke: 0.5pt + rgb("#c4c1b6"),
      align: (left, left, left),
      inset: 7pt,
      table.header(
        [*Attribute*], [*De Bruijn*], [*Overlap*],
      ),
      [Node represents], [$(k-1)$-mer], [Read],
      [Edge represents], [Single $k$-mer], [Read–read suffix-prefix overlap],
      [Construction cost], [$O$(total read length)], [$O$(reads²) naive],
      [Search], [Eulerian (polynomial)], [Hamiltonian (NP-hard)],
      [Memory scales with], [Graph size], [Read count],
      [Best for], [Short accurate reads], [Long reads],
    )
  ],
  caption: [
    De Bruijn versus overlap graphs at a glance. The choice tracks
    read length: short reads are unworkable as overlap-graph nodes;
    long reads are unworkable as de Bruijn-graph edges.
  ],
) <fig:graph-choice>

For reads under 500 bp you effectively have to use de Bruijn — an
overlap graph would be quadratic in the number of reads, and at
$10^9$ reads that is $10^(18)$ comparisons. For reads over 5 kb you
effectively have to use overlap — de Bruijn with $k approx 1000$ has
no redundancy to leverage.

#warn[
  Hybrid assemblers — SPAdes-hybrid, Unicycler, MaSuRCA — use both.
  De Bruijn builds a clean short-read backbone; long-read overlaps
  bridge the repeats and gaps the short reads cannot resolve. If you
  have both data types, use a hybrid assembler. Running short-read
  and long-read assemblies separately and "merging" them afterwards
  almost always gives worse results than letting the assembler use
  the two signals jointly during graph construction.
]


== Contigs and Topology <sec:topology>

The graph is not the assembly. The assembly is a set of strings — the
*contigs* — each representing a region of the genome the assembler is
confident in.

A contig is produced by walking a *non-branching path*: a maximal
path where every internal node has in-degree one and out-degree one.
The walk terminates at a branch point (a node with multiple in-edges
or multiple out-edges) or at a dead end. Each non-branching path
becomes one contig.

Why break at forks? Because forks are exactly the places where more
than one reconstruction is consistent with the reads. The assembler
could guess, but guessing hides uncertainty. Emitting two separate
contigs at a fork is the honest output.

For a clean bacterial genome at 30× coverage with no repeats longer
than the read length, the de Bruijn graph is close to a single linear
path and the assembly is one contig covering the whole #idx("chromosome")chromosome.
Real genomes are not like that. Real genomes have repeats, errors, and
structural variation, and the graph has forks.

=== Three signatures, three responses <sec:topology-signatures>

Three patterns of "not-a-linear-path" show up in every assembly graph.
They are the assembler's basic vocabulary.

- *Tips* are short dead-end branches a few edges long. They come from
  single sequencing errors near the end of a read — the error creates
  a unique $k$-mer that branches off the true path and terminates
  because no other read extends it. The cleanup is *tip clipping*:
  remove any tip shorter than a threshold, typically about $2k$.
- *Bubbles* are parallel paths between the same pair of branch nodes.
  They come from SNPs (one read carries the major allele, the other
  the minor), from small indels, or from sequencing errors that happen
  to rejoin the true path further down. The cleanup is *bubble
  popping*: identify the two paths, compare their coverage, keep the
  higher-coverage path, drop the other.
- *Tangles* are densely interconnected subgraphs with many in- and
  out-edges. They come from repeats: every copy of a repeat collapses
  onto the same path in the graph, and reads coming into the repeat
  can leave through any of the copies. There is no clean cleanup.
  Tangles get reported as unresolved regions; long reads or paired-end
  information are needed to walk them.

#figure(
  image("../../diagrams/lecture-03/08-topology-signatures.svg", width: 100%),
  caption: [
    The three canonical topology signatures. Tips come from single-read
    errors; bubbles from small variants; tangles from repeats.
    Recognising the pattern is half of cleanup.
  ],
) <fig:topology-sig>

#figure(
  image("../figures/ch03/f3-topology-cleanup.svg", width: 100%),
  caption: [
    Before-and-after the two main cleanup operations. Tip clipping
    removes the dead-end branch off the main path; bubble popping
    collapses the two parallel paths into the higher-coverage one.
    Coverage is the discriminator in both cases.
  ],
) <fig:topology-cleanup>

#note[
  The topology of the graph is the signal pattern of different error
  types. A sequencing error shows up as a tip one or two edges long.
  A biological variant shows up as a bubble exactly two paths wide.
  A repeat of length $R$ and copy number $N$ shows up as a tangle
  with $N$ entries and $N$ exits. Each has a distinct signature, and
  the assembler's cleanup stage is essentially pattern-matched
  remediation — matched filter, then take the action the match
  implies.
]

=== Why bubble popping is harder than it sounds

Tip clipping is unambiguous. A two-edge dead-end branch with coverage
one in a graph where everything else has coverage twenty-five is, with
overwhelming probability, an error; the cleanup is to delete it.

Bubble popping has more subtlety. If the two parallel paths have
coverages 32 and 2, the high-coverage path is real and the low-coverage
path is the error; pop the low-coverage one. If the two paths have
coverages 16 and 14, you are almost certainly looking at a heterozygous
variant — the genome carries _both_ paths, one on each homologous
chromosome, and collapsing them throws away real biology. A
diploid-aware assembler keeps the bubble and emits two contigs, one
per #idx("haplotype")haplotype. A haploid-only assembler collapses it and notes the
ambiguity. The choice depends on what you intend to do with the
assembly downstream.

Modern long-read assemblers (hifiasm in particular) make the diploid
distinction first-class: they emit a *primary* assembly plus an
*alternate* haplotype assembly, with bubbles separated rather than
collapsed. For a heterozygous human sample at 30× HiFi, that is the
right answer; for a haploid bacterial sample it is a waste of memory.


== Repeats: The Fundamental Obstacle <sec:repeats>

Repeats are the limit. Everything else in assembly is engineering;
repeats are physics.

A *repeat* is a sequence that appears at more than one position in the
genome. The human genome is roughly half repeat by content. The
catalogue is worth knowing in numbers:

- *Alu elements*, primate-specific SINEs, about 300 bp long, about
  1.1 million copies, occupying about 10% of the genome;
- *LINE-1 elements*, about 6 kb long, about 500,000 copies, about 17%
  of the genome;
- *Centromeric and pericentromeric satellite arrays*, megabase-scale
  arrays of short motifs (171 bp for the human alpha-satellite),
  resolved only with the most recent ultra-long-read data;
- *Segmental duplications*, large blocks of 10 kb or more with
  &gt; 90% identity to another copy elsewhere, totalling about 5% of
  the genome;
- *Ribosomal #idx("RNA")RNA arrays*, hundreds of tandem copies of the rDNA
  cluster on five different chromosomes.

Here is why they break assembly. Take two identical copies of a 500 bp
repeat at positions $X$ and $Y$. For any $k <= 500$, every $k$-mer in
the repeat appears in both copies of the repeat. In the de Bruijn
graph the two copies *collapse into a single path* — the graph has
one path for the repeat, with the two upstream regions converging into
it and the two downstream regions diverging out of it. The assembler
sees the repeat. It cannot tell which upstream flank belongs to which
downstream flank.

#figure(
  image("../../diagrams/lecture-03/09-repeat-collapse.svg", width: 100%),
  caption: [
    Two identical repeat copies at positions X and Y collapse into one
    graph node, with the two upstream flanks converging in and the two
    downstream flanks diverging out. The assembler cannot pair upstreams
    with downstreams from the reads alone.
  ],
) <fig:repeat>

#warn[
  Any repeat longer than the read length cannot be resolved from the
  reads alone — the reads do not span it, so there is no information
  to distinguish the copies. Long reads help (if the read spans the
  repeat, the repeat is resolved). Mate-pair libraries help (if the
  insert spans the repeat, scaffolding can bridge it). But no
  algorithm resolves _all_ repeats. Published genomes with "complete"
  in the title usually have gaps inside the centromeric or telomeric
  repeat arrays — the T2T human assembly in 2022 was the first
  publication to genuinely close them, twenty years after the
  "complete" draft.
]

Assemblers handle repeats in three ways, in increasing order of
honesty:

1. *Collapse and report.* The repeat appears once in the assembly with
   elevated coverage. This is what short-read-only assemblers do for
   long repeats. Copy number is lost; the assembler hopes someone
   downstream notices the inflated coverage.
2. *Fork and report both.* The assembler emits two contigs at the
   repeat boundary and marks them as alternate paths. This requires
   reasoning about coverage and works best when the copies have
   accumulated enough divergent mutations to be slightly different.
3. *Span with long reads.* A read that is longer than the repeat
   unambiguously resolves which upstream flank belongs to which
   downstream flank. This is why long-read assembly is necessary for
   any mammalian genome — short reads simply cannot bridge the longest
   classes of repeats.

The third approach is the only one that gives you a complete sequence
through a repeat. Everything else is bookkeeping over the ambiguity.


== Scaffolding: Linking Contigs Across the Gap <sec:scaffolding>

After contig construction and graph cleanup, the assembler has a set
of contigs. Their internal order — which contig is to the left of which
— is unknown from the contigs alone. The reads inside each contig pin
down its sequence, but the reads inside contig $C_1$ have no information
about contig $C_2$ unless a read or read-pair touches both.

*Scaffolding* is the post-step of ordering and orienting contigs
across the unresolved regions, using information that was set aside
during graph construction.

The two workhorses:

- *Paired-end reads.* A pair of reads from the same approximately
  500 bp fragment. If one mate maps to contig $A$ and the other to
  contig $B$, then $A$ and $B$ are within roughly 500 bp of each other
  on the chromosome, and you know their relative orientation from the
  pair geometry.
- *Mate-pair libraries.* The same idea with larger inserts — 2 kb to
  40 kb. The long inserts jump over short repeats that paired-end
  cannot bridge. They are more expensive to prepare (the protocol
  circularises the fragment before sequencing, which adds molecular-
  biology steps and reduces complexity) but they are the difference
  between an assembly that is a hundred fragments and one that is a
  dozen.

#figure(
  image("../../diagrams/lecture-03/10-scaffolding-mate-pairs.svg", width: 100%),
  caption: [
    Scaffolding with mate pairs. Long-insert read pairs link contigs
    across unresolved regions; the resulting #idx("scaffold")scaffold preserves order
    and orientation with explicit N-padded gaps.
  ],
) <fig:scaffolding>

Modern equivalents include *Hi-C* (chromosome-conformation capture,
which links contigs that sit near each other in three-dimensional
nuclear space even when they are megabases apart in the linear
sequence) and *optical maps* (Bionano Genomics), which provide long-
range order independent of the reads.

The output of scaffolding is a *scaffold*: a sequence consisting of
contigs linked by gaps. Gaps are filled with `N` characters whose
length is the estimated gap size from the linking evidence. A scaffold
of length 2 Mb might be 1.8 Mb of real sequence and 0.2 Mb of `N`s
distributed across a dozen gaps. The `N`s carry no signal — they are
a placeholder that downstream tools can recognise and skip over. A
guess at what is in the gap would silently corrupt variant calling,
expression quantification, and any other analysis that maps onto the
assembly. Gaps are honest about what the data does not resolve.

#note[
  Mate-pair ("jumping") libraries were introduced by the *Celera
  Genomics* human-genome assembly (Myers, Sutton, and others, 2000)
  precisely to solve the repeat problem. The whole-genome shotgun
  approach without mate-pairs would have produced a million unordered
  fragments. With 2 kb, 10 kb, and 50 kb insert libraries layered
  together, Celera delivered an assembly competitive with the public
  HGP's clone-by-clone approach at roughly one-fifth the cost. The
  technique is dated by Hi-C and ultra-long reads today, but the
  insight — use insert-size diversity to bridge what any single insert
  cannot — is still the bedrock of every modern assembly.
]


== Measuring an Assembly <sec:metrics>

How do you tell a good assembly from a bad one?

There is no single number, but if there were it would be *N50*.

=== N50 and its variants

Sort the contigs by length, largest first. Walk down the sorted list
adding up lengths. *N50 is the length of the contig at which the
running total first exceeds half of the total assembly length.*
Equivalently: N50 is the largest length $L$ such that contigs of
length $>= L$ together cover at least 50% of the total assembly.

#note[
  N50 is the median of the contig-length distribution _weighted by
  length_. The usual median is the middle contig by count; N50 is
  the middle contig by base. That weighting is what makes the metric
  biologically meaningful — you care more about where half the genome
  sits than about where the middle contig out of ten thousand sits.
]

A worked example helps the definition land. Suppose your assembly
contains ten contigs with lengths (in kilobases) 800, 600, 500, 400,
300, 250, 200, 150, 100, 50. The total assembly is 3,350 kb. The
running cumulative sum is 800, 1400, 1900, 2300, 2600, 2850, 3050,
3200, 3300, 3350. Half of 3,350 is 1,675. The cumulative sum first
crosses 1,675 at the third contig (cumulative 1,900). Therefore
*N50 = 500*. The companion metric *L50 = 3* tells you how many contigs
were needed to cover half the assembly.

#figure(
  image("../figures/ch03/f4-n50-walkthrough.svg", width: 100%),
  caption: [
    N50 and #idx("NG50")NG50 computed on a ten-contig assembly. The cumulative-length
    curve crosses the 50%-of-total reference inside contig 3
    (N50 = 500 kb) and crosses the 50%-of-expected-genome reference
    inside contig 4 (NG50 = 400 kb). The summary table at lower left
    is the five-number report that always accompanies these metrics.
  ],
) <fig:n50>

=== NG50 and the contig curve

*NG50* is the same idea but using *expected genome size* as the
denominator rather than total assembly length. If your 3,350 kb
assembly is really half of a 6,700 kb genome (i.e. 50% incomplete),
N50 will still report a contig from inside the assembled half, but
NG50 will report a much smaller contig — possibly zero, if the
assembly does not even cover half of the genome. NG50 is the metric
to use when comparing assemblers on the same genome, because
incomplete-but-fragmented assemblies will not look misleadingly good.

The full curve of $N_x$ as $x$ varies from 0 to 100 — the
*assembly curve* — is more informative than any single point on it.
$N_{10}$ tells you about your longest few contigs; $N_{50}$ tells
you about the middle; $N_{90}$ tells you how long the small-contig
tail is. The shape of the curve diagnoses common failure modes:
a steep $N_{10}$ followed by a sharp drop into a long flat tail of
small contigs is the classic "one or two big collapsed-repeat contigs
plus a long tail of unresolved fragments."

#warn[
  N50 is easy to game. Aggressive joining through uncertain regions
  raises N50 while lowering accuracy. Always report N50 alongside
  *total length*, *number of contigs*, *largest contig*, and
  *genome-fraction covered*. A five-number summary tells you more than
  any single metric and is much harder to cheat — you cannot inflate
  one without deflating another.
]

=== Beyond length: consistency, completeness, and misassemblies

Length is not correctness. Three further metrics are routinely reported.

- *Assembly consistency.* The fraction of the assembly that aligns
  correctly to a trusted reference, when one exists. Requires a
  reference; useful for benchmarking on known organisms.
- *Genome completeness.* The fraction of the expected genome that the
  assembly covers. Tools like *BUSCO* estimate completeness without
  needing the full reference, by checking how many of a curated set of
  single-copy orthologs from related organisms appear once in the
  assembly. A BUSCO completeness score of 98% with 1% duplication
  rate is the kind of summary that ships in modern publications.
- *Misassembly count.* The number of places the assembly disagrees
  with a reference in ways not explainable by real variation —
  inverted segments, translocated segments, spurious junctions.
  *QUAST* computes them by aligning the assembly back to a reference
  and looking for break points.

The honest report on an assembly looks something like this: total
length 4.6 Mb across 87 contigs, largest contig 530 kb, N50 145 kb,
NG50 138 kb (expected genome 4.8 Mb), BUSCO completeness 97%, 0
misassemblies against the type strain. That is more useful than any
of its constituents in isolation, and it is a five-line summary that
anyone can produce from a single `quast` invocation.


== A Whole-Genome Shotgun, End to End <sec:wgs>

Putting it all together: how a bacterial genome gets assembled from a
tube of DNA to a polished FASTA.

1. Extract DNA. Fragment to 300–500 bp. Build a sequencing library;
   sequence on Illumina to get roughly 30× of 150 bp reads. Total cost
   for a small bacterial genome in 2024: about \$100, dominated by
   reagent and instrument time.
2. Run `fastp` on the reads for adapter trimming and quality
   filtering. Drop reads shorter than 50 bp. This usually removes a
   percent or two of the raw data.
3. Run *SPAdes* in `--isolate` (or `--careful`) mode. SPAdes
   internally runs BayesHammer error correction, builds de Bruijn
   graphs at several values of $k$, merges them, cleans tips and
   bubbles, and emits contigs.
4. (Optional) Sequence a second library: mate-pairs, or a PacBio HiFi
   run, for scaffolding.
5. Run a scaffolder (`SSPACE`, `BESST`, or SPAdes's built-in) to order
   and orient contigs across the unresolved regions.
6. Run `Pilon` (or `Racon` for long-read polishing) to correct any
   residual single-base errors.
7. Run `QUAST` to compute the metrics from @sec:metrics.

Total compute for a 5 Mb bacterial genome is roughly 30 minutes on a
modern laptop. The same pipeline on a 3 Gb mammalian genome requires
roughly a hundred times the compute, multi-day runtime, and long reads
throughout. Modern mammalian-genome assemblies are PacBio HiFi with
*hifiasm* or *HiCanu*; they reach contig N50 in the tens of megabases
and, with the addition of ultra-long ONT reads and Hi-C scaffolding,
reach chromosome-arm N50 values for the first time in 2022.

#note[
  The shotgun-versus-clone-based debate that drove the 1998-2001
  human-genome race is settled: shotgun won. The hierarchical-clone
  approach used by the public Human Genome Project produced a slightly
  higher-quality reference but cost roughly ten times more. Celera's
  shotgun approach — the technique described in this section — became
  the universal standard for every genome sequenced since, from
  microbes to mammals. The unresolved repeats that broke the original
  shotgun assemblies are now resolved by long reads, not by going back
  to clones.
]

=== Working with the BAM after assembly

The first thing you do with a finished assembly is align the reads
back to it. The motivation is mostly self-check: estimate per-contig
coverage, detect residual errors, look for contamination, verify the
assembler did not miss anything obvious. Those alignments land in a
*BAM* file — the binary indexed form of SAM, whose field layout was
introduced in Chapter 1. The standard tools all work:

- `samtools sort` / `index` — prepare the BAM for fast per-region
  access.
- `samtools depth` — per-position coverage; the input to every
  assembly-QC plot.
- `samtools flagstat` — sanity-check mapping rates and pair concordance.
- `bcftools mpileup` or `DeepVariant` — call variants from the pileup
  if you intend to compare the same reads against the assembled
  sequence (a useful internal-consistency check).

#warn[
  An assembly BAM is not interchangeable with a reference BAM. The
  reads are aligned to the assembly itself — the contigs — not to a
  reference genome. The reference column carries contig names like
  `NODE_1_length_12345`, not chromosome names like `chr3`. Tools that
  assume reference-aligned BAMs will silently produce nonsense on an
  assembly BAM. Always confirm what was used as the reference before
  feeding a BAM downstream.
]


== Summary <sec:summary>

- *De novo assembly is template-free reconstruction.* When no reference
  exists, assembly is the only way to get a sequence. The problem is
  fundamentally harder than resequencing because there is no signal
  to match against.
- *The de Bruijn graph is the central algorithmic idea of short-read
  assembly.* It turns reconstruction into an Eulerian-path problem,
  which is polynomial-time and scales to mammalian genomes. Every
  serious short-read assembler since 2008 is a de Bruijn variant.
- *Overlap graphs remain the right tool for long reads.* No redundancy
  to exploit at $k = 1000$; the assembler instead works at the read
  level and tolerates the Hamiltonian-path search because the graph is
  almost-linear in practice.
- *Coverage is a Poisson process with systematic non-uniformity.*
  Mean coverage is a weak summary; the distribution — especially its
  left tail at low depth — tells you where the assembly is going to
  fail.
- *Repeats are the limit.* No algorithm fully resolves repeats longer
  than the read. Long reads, mate-pairs, and Hi-C all help; none of
  them solve the problem. Assembly-with-gaps is the honest output.
- *N50 is the median contig length weighted by length.* Useful but
  gameable; always pair it with total length, contig count,
  largest contig, and genome-fraction-covered.


== Exercises <sec:exercises>

#strong[1.] _Computational._ Download a small bacterial genome —
_E. coli_ K-12 MG1655 is the canonical choice at roughly 4.6 Mb —
from NCBI. Simulate 30× Illumina-style reads at 150 bp with
`wgsim` using its default error model. Assemble the simulated reads
with SPAdes in `--isolate` mode. Report number of contigs, largest
contig, N50, NG50 (against the known true genome length), and total
assembly length. Compare these to the true genome and account for
the difference.

#strong[2.] _Pencil-and-paper._ By hand, construct the de Bruijn graph
with $k = 4$ from these three reads: `ACAGACGT`, `CAGACGTA`,
`AGACGTAC`. Draw the graph. Find the Eulerian path. Read off the
contig. Show your work for each step.

#strong[3.] _Conceptual._ Explain in one paragraph why a perfect
tandem repeat of length exactly $k - 1$ base pairs will _not_
collapse in a de Bruijn graph with parameter $k$. (Hint: think about
what happens to the $k$-mer that straddles the boundary between the
two repeat copies.)

#strong[4.] _Computational._ Given the Phred quality string and a
coverage histogram from a real sequencing run (any FASTQ from SRA),
write a Python program that fits a two-component model (true peak as
Poisson around mean $lambda$, error peak as exponential decay from 1)
to the $k$-mer histogram, and reports the maximum-likelihood threshold.
Compare your threshold to the one BayesHammer selects on the same data.

#strong[5.] _Quantitative._ A 3 Gb diploid genome is sequenced to 30×
coverage on an instrument with 1% per-base error and 150 bp reads.
At $k = 31$, estimate (a) the expected total number of distinct
$k$-mers in the read set, (b) the expected fraction of those that
contain at least one sequencing error, and (c) the expected count of
the error-k-mer peak at frequency 1. State your assumptions.

#strong[6.] _Conceptual._ Bubble popping discards the lower-coverage
path. In a diploid genome at 30× coverage, a heterozygous variant
produces a bubble with two paths at roughly 15× each. What happens to
the variant if you run a haploid assembler? What does a diploid-aware
assembler do instead, and how does it tell heterozygous bubbles from
sequencing-error bubbles?

#strong[7.] _Open-ended._ Pick one published assembler — SPAdes,
Flye, hifiasm, Canu, MEGAHIT, or another of your choice — and read
its main paper. In one page, identify (a) the specific design
decisions it makes that differ from the textbook de Bruijn / overlap
recipes, (b) the failure modes those decisions trade off, and (c) the
class of input data for which it is currently the best available
tool. Cite the paper.


== Further Reading <sec:further-reading>

- *de Bruijn, N. G.* (1946). "A Combinatorial Problem." _Proceedings
  of the Koninklijke Nederlandse Akademie van Wetenschappen_ 49:
  758–764. The original paper, six pages, all combinatorics. Worth
  reading once to see how far the application has travelled.
- *Idury, R. M., and Waterman, M. S.* (1995). "A New Algorithm for
  DNA Sequence Assembly." _Journal of Computational Biology_ 2:
  291–306. The first application of de Bruijn graphs to genome
  assembly.
- *Pevzner, P. A., Tang, H., and Waterman, M. S.* (2001). "An
  Eulerian Path Approach to DNA Fragment Assembly." _PNAS_ 98:
  9748–9753. The reformulation that made assembly polynomial.
- *Myers, E. W., Sutton, G. G., Delcher, A. L., et al.* (2000). "A
  Whole-Genome Assembly of _Drosophila_." _Science_ 287: 2196–2204.
  The Celera paper. The first OLC assembly of a multicellular genome
  and the methodological foundation of every modern long-read
  assembler.
- *Zerbino, D. R., and Birney, E.* (2008). "Velvet: Algorithms for
  De Novo Short Read Assembly Using De Bruijn Graphs." _Genome
  Research_ 18: 821–829. The first practical de Bruijn assembler.
- *Bankevich, A., Nurk, S., Antipov, D., et al.* (2012). "SPAdes: A
  New Genome Assembly Algorithm and Its Applications to Single-Cell
  Sequencing." _Journal of Computational Biology_ 19: 455–477. The
  modern short-read workhorse.
- *Koren, S., Walenz, B. P., Berlin, K., et al.* (2017). "Canu:
  Scalable and Accurate Long-Read Assembly via Adaptive k-mer
  Weighting and Repeat Separation." _Genome Research_ 27: 722–736.
  The long-read counterpart.
- *Compeau, P. E. C., Pevzner, P. A., and Tesler, G.* (2011). "How
  to Apply de Bruijn Graphs to Genome Assembly." _Nature
  Biotechnology_ 29: 987–991. A two-page primer; the cleanest
  short summary in print.
