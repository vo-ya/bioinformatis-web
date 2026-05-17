#import "../theme/book-theme.typ": *

= Bulk RNA-seq: From Reads to Transcript Abundances <ch:bulk-rna-seq>

#matters[
  A genome lists the parts a cell could in principle build. A
  transcriptome reports the parts the cell is actually building right
  now. The same DNA sits in every nucleus of your body; a liver cell and
  a neuron differ because they transcribe different subsets of it, at
  different rates, into different splice variants. RNA-seq is the
  measurement that turns that activity into numbers — and almost every
  paper on disease mechanism, drug response, developmental biology, or
  cancer prognosis that you will read in the next ten years opens with
  some flavour of it. The numerics carry consequences: a botched
  normalisation can flip the sign of an effect; a misconfigured strand
  setting can silently halve every count in the matrix. This chapter
  builds the pipeline from FASTQ to a clean count matrix and explains
  what each stage is buying you.
]

The output of a sequencing instrument is the same whether the input was
DNA or RNA — a FASTQ file of short reads with PHRED quality scores. The
divergence starts the moment you try to align those reads back to a
genome. A read from a piece of genomic DNA maps to a contiguous stretch
of reference. A read from a piece of mature mRNA does not: the cell has
already spliced its introns out, so the read corresponds to a
discontinuous slice of genome with multi-kilobase gaps inside it. The
rest of the chapter follows from that one fact.

Four moves carry the rest. Section 5.1 lays out the biology — splicing,
isoforms, and why naive read-per-gene counting is wrong on two
independent axes. Section 5.2 explains splice-aware alignment as a
hidden-Markov-model extension of the seed-and-extend framework from
Chapter 3, comparing the two dominant short-read aligners STAR
(Dobin 2013) and HISAT2 (Kim 2015, 2019). Section 5.3 takes the radical
alternative — pseudoalignment, introduced by Kallisto (Bray 2016) and
refined by Salmon (Patro 2017) — which throws away the base-level
alignment that quantification didn't need and runs 10–100× faster as a
consequence. Section 5.4 covers the expectation-maximisation algorithm
that resolves the ambiguity in compatibility classes, from RSEM (Li &
Dewey 2011) through to its modern descendants. Section 5.5 turns to the
count matrix itself: CPM, FPKM, TPM, and the median-of-ratios size
factor of DESeq2. Section 5.6 explains why count data is overdispersed
Poisson and what the negative-binomial fix buys you. Section 5.7 surveys
the broader assay family — ribo-seq, CLIP-seq, m^6A-seq, CAGE — which
all reuse this pipeline downstream of FASTQ.

The chapter assumes a working knowledge of short-read alignment from
Chapter 3 (seed-and-extend, BWT/FM-index) and the Phred-quality formalism
from Chapter 1. Differential expression — the statistical machinery that
takes the count matrix out the other side of this chapter and decides
which genes change between conditions — is the subject of Chapter 6.


== Why RNA Is Not DNA <sec:rna-vs-dna>

A short Illumina read coming off the instrument is 150 bp of A/C/G/T
with a PHRED string. There is no flag on the FASTQ record that says
"this came from genomic DNA" versus "this came from polyA-selected RNA."
The biology is hiding entirely in what the reference looks like, which is
why the same alignment tool you used in Chapter 3 will silently produce
nonsense if you point it at an RNA-seq library and ask for a contiguous
genomic alignment.

Three biological facts drive everything in this chapter.

*Transcription is selective.* Only a fraction of the genome is
transcribed in any given cell at any given time. Genes have promoters
that recruit RNA polymerase II; transcription factors bind regulatory
regions and turn genes on or off; the abundance of a transcript in a
cell is approximately proportional to its production rate times its
half-life. A liver hepatocyte transcribes the ~12,000 genes most
relevant to liver function at non-trivial levels; the other ~8,000
protein-coding genes in the genome are mostly silent.

*Eukaryotic genes are split.* Human genes are encoded as a series of
*exons* — the parts that end up in the mature mRNA — separated by
*introns* — the parts that get cut out. The spliceosome, a large
RNA-protein complex, recognises canonical splice-donor (GT) and
splice-acceptor (AG) dinucleotides at exon boundaries, excises the
intron, and joins the exons. A gene spanning 40 kb of genomic locus
typically produces a 2 kb mature mRNA: most of the locus is intron.
The largest human gene, *DMD*, spans 2.4 Mb on the X chromosome and
encodes a 14 kb mRNA. About 99.4 % of *DMD*'s genomic length is intron.

*One gene, many transcripts.* The same gene can be spliced in multiple
ways — different combinations of exons assembled into different final
mRNAs. These are called *isoforms* or splice variants. Human genes
average roughly four isoforms each; a small fraction have dozens.
Isoforms can differ in their protein product (an exon-skipped form may
lose a regulatory domain), their stability (alternative 3′ UTRs change
miRNA susceptibility), or just their localisation. From the alignment
algorithm's perspective, the same gene's two isoforms are two distinct
reference sequences that share most of their bases but disagree at the
junctions.

#figure(
  image("../../diagrams/lecture-05/01-splicing-isoforms.svg", width: 92%),
  caption: [
    Transcription and splicing. The same genomic locus produces
    different mature mRNAs depending on which exons survive splicing.
    This is the move that breaks the DNA-style alignment from Chapter 3.
  ],
) <fig:splicing>

These three facts make RNA-seq a fundamentally different alignment
problem from DNA-seq. The sample is not a copy of the reference; it is
a _transformed_ version of it. A short read drawn from a mature mRNA
may straddle an exon-exon junction; in the genome those exons are
separated by a gap of hundreds to hundreds of thousands of bases. A
contiguous aligner cannot place that read. Worse, the same read could
be _equally consistent_ with three different isoforms of the same gene,
because all three share the exonic content the read covers. The basics
of counting break down at every step.

=== Why Naive Counting Fails

Given the biology in @fig:splicing, the naive pipeline — "align reads
with BWA from Chapter 3 and count how many overlap each gene" — breaks
in three independent ways.

First, *reads span splice junctions*. A 150 bp read can sit across the
boundary between two exons. In the genome those exons are separated by
an intron — possibly 50 kb long. A DNA-style aligner that tries to map
the read contiguously to the reference either marks it unmapped or
splits it into two low-quality alignments. STAR and HISAT2 exist
specifically to handle these junction-spanning reads as a first-class
alignment case.

Second, *reads from one gene can't always be assigned to one isoform*.
A read that lands in an exon shared by three isoforms of the same gene
is, at the read level, ambiguous: we know the gene but not the
transcript. The naive solution — assign to a primary transcript and
move on — is wrong in a principled way, because the same gene's
isoforms can have radically different functional consequences in
disease. Salmon and Kallisto resolve the ambiguity with
expectation-maximisation, treating each read's contribution to each
compatible transcript as a fractional vote.

Third, *raw read counts conflate abundance with transcript length*. A
10 kb transcript drawn at the same molar abundance as a 1 kb transcript
collects ten times more reads, because each read is sampled from a
random position along the transcript. Raw read counts are a biased
estimator of molar abundance — you must normalise by length to compare
genes within a sample.

#warn[
  "Reads per gene" and "transcript abundance" are different quantities,
  in different units, with different rankings. Plotting one and
  labelling it as the other is the single most common bug in published
  RNA-seq analyses. A gene with high read count but very long length is
  not necessarily more expressed than a short gene with fewer reads. The
  units matter and this chapter is partly about getting them right.
]

=== The Pipeline at a Glance

Every bulk RNA-seq analysis follows one of two canonical paths. They
disagree on the middle stage but agree on the endpoints.

#figure(
  image("../../diagrams/lecture-05/02-rnaseq-pipeline.svg", width: 95%),
  caption: [
    The bulk RNA-seq pipeline. Two paths differ on the middle stage —
    align-then-count via STAR/HISAT2 plus `featureCounts`, or
    pseudoalign via Kallisto/Salmon — but both converge at a gene-level
    count matrix that Chapter 6 will consume.
  ],
) <fig:pipeline>

The *align-then-count* path produces a sorted BAM and then counts reads
overlapping each feature (gene, exon, or transcript) using
`featureCounts` or `htseq-count`. The *pseudoalign* path produces
transcript-level abundance estimates directly and never writes a BAM at
all. Both end at the same destination: a gene-level count matrix that
differential expression tools consume.

The pragmatic 2024 default is Salmon followed by `tximport` into
DESeq2; the legacy default — still right in many cases, especially when
you need the BAM for other purposes such as variant calling from RNA —
is STAR followed by `featureCounts`. The two paths land within a few
percent of each other on standard benchmarks.

#note[
  The historical timeline of bulk RNA-seq quantification is short and
  worth knowing. Tophat (Trapnell 2009) was the first widely-used
  splice-aware aligner. RSEM (Li & Dewey 2011) introduced
  expectation-maximisation on top of Bowtie alignments — the first
  principled treatment of multi-mapping reads. STAR (Dobin 2013) and
  HISAT (Kim 2015) made splice-aware alignment fast. Kallisto (Bray
  2016) replaced alignment with pseudoalignment and brought
  quantification down to a few minutes per sample. Salmon (Patro 2017)
  added bias-aware quantification on top of pseudoalignment. By 2018,
  the workflow most labs run today had stabilised.
]


== Splice-Aware Alignment <sec:splice-aware>

The splice-aware aligner is the workhorse that bridges short reads and a
genomic reference when the input is mRNA. Two tools dominate: STAR and
HISAT2. They solve the same problem with different data structures.

=== The Exon-Intron Problem

Concretely, imagine an exon of 100 bp followed by an intron of 5 kb
followed by the next exon. A read that straddles the exon-exon junction
will have its first 80 bp aligning to the end of exon 1 and its last 70
bp aligning to the start of exon 2. In between, on the reference, sits
5 kb of intron that the read does not contain. The aligner has to
recognise that 5 kb gap as a *splice* rather than as a deletion. A
deletion that large would be either a structural variant (Chapter 4) or
a complete misalignment; a splice is a routine biological event.

#figure(
  image("../../diagrams/lecture-05/03-junction-reads.svg", width: 90%),
  caption: [
    Three categories of reads on a gene with three exons. Reads that
    fall inside one exon are easy; junction reads need a single
    long-skip move; multi-junction reads need two.
  ],
) <fig:junction-reads>

The naive solution — let the aligner open arbitrarily large gaps — fails
because a 150 bp read against a 3 Gb genome will spuriously match many
random 150-bp slices with a few edits. Allowing 5 kb gaps would produce
millions of false-positive alignments per sample. What rescues the
problem is that *legitimate splice junctions are constrained*: they
occur almost exclusively at canonical GT-AG dinucleotides (with the
remaining few percent split between the GC-AG and AT-AC alternatives),
and the vast majority are pre-annotated in reference transcriptomes
like GENCODE or Ensembl. Splice-aware aligners exploit both: they allow
intron-sized gaps in the alignment, but only at sites that look like
canonical splice motifs and ideally already appear in the annotation.

#note[
  A junction-spanning read is, mechanically, an insertion-deletion of
  length 1–500 kb whose position is anchored to known splice-site
  motifs. The DNA-style aligners of Chapter 3 allow indels up to ~50 bp
  and assume they can appear anywhere. RNA-seq aligners allow indels up
  to half a megabase but only at biologically sanctioned positions.
  Same algorithmic shape, very different constraint geometry.
]

=== STAR

STAR — *S*pliced *T*ranscripts *A*lignment to a *R*eference — is the
dominant short-read RNA aligner. It extends the seed-and-extend
framework from Chapter 3 with three modifications, tuned for long
intron-sized jumps.

The algorithm has three phases per read.

+ *Seed search.* Instead of looking for one long exact match, STAR
  searches for the *maximal mappable prefix* (MMP) starting from the
  first base of the read. The MMP is the longest prefix of the read
  that occurs exactly somewhere in the genome. This is done with a
  suffix-array index over the reference. When the prefix grows beyond
  what matches uniquely, STAR records the MMP, jumps forward in the
  read by the matched length, and starts a new MMP from the next
  unaligned base. A 150 bp read might produce two or three MMPs.

+ *Seed clustering.* Multiple MMPs from one read are clustered by
  genomic location. Two MMPs within ~50 kb of each other on the same
  chromosome probably belong to the same alignment, separated by an
  intron. MMPs that land further apart, or on different chromosomes,
  are either chimeric (a structural-variant-like RNA fusion) or
  multi-mapping, and are demoted to non-primary alignments.

+ *Extend with splice awareness.* For each cluster, STAR extends the
  MMPs using a modified Smith-Waterman alignment that allows an unusual
  move: at any canonical GT…AG motif, the alignment can "jump" forward
  by the intron length at *zero* alignment cost. The splice move is a
  different operation from an ordinary insertion — insertions have a
  per-base gap-open and gap-extend penalty; splices are free if the
  boundary lands at a canonical splice motif.

#figure(
  image("../../diagrams/lecture-05/04-star-mechanics.svg", width: 92%),
  caption: [
    STAR's three phases: three MMPs are recovered by suffix-array
    search; two cluster within fifty kilobases on the same chromosome;
    the third is discarded; the surviving cluster is extended across
    the intron at a canonical GT-AG motif.
  ],
) <fig:star-mechanics>

@fig:star-trellis spells the same algorithm out as an explicit trellis
walk for a single 150-bp junction-spanning read. The take-home is that
the splice edge is one move in a state machine that otherwise looks
identical to the DNA aligner of Chapter 3 — the algorithmic machinery
is unchanged, the state space is enriched with one biology-aware
transition.

#figure(
  image("../figures/ch05/f1-star-trellis-walk.svg", width: 95%),
  caption: [
    The Viterbi-best alignment of a 150-bp junction-spanning read
    through STAR's seed cluster and extension trellis. The splice edge
    costs zero because the boundary lands on a canonical GT-AG motif;
    the alternative gap-extension path costs five thousand
    gap-extensions and loses.
  ],
) <fig:star-trellis>

#tip[
  Splice-aware alignment is a hidden Markov model whose state-transition
  diagram includes long-skip transitions at canonical splice sites. A
  junction alignment is the Viterbi path through a trellis where one
  edge can traverse 100 kb of reference in a single time step — but only
  if that edge is labelled "canonical splice." Drop the long-skip edges
  and you recover the DNA aligner. Add them and you recover STAR. The
  algorithmic vocabulary an EE student already has — Viterbi decoding,
  trellis search, constrained state transitions — carries over almost
  unchanged.
]

STAR's practical envelope: a human-genome index is about 30 GB of RAM
(the suffix array over the whole genome plus a junction-database
overlay). The aligner processes ~100 million reads per hour on 16
cores. The output is a coordinate-sorted BAM with `XS` tags indicating
splice strand and, optionally, a tab-separated table of detected splice
junctions per sample.

A second-pass mode is worth knowing. STAR's *two-pass mapping* first
runs the aligner with only the GENCODE-annotated junctions; the
detected novel junctions from pass one are added to the junction
database; the aligner runs again. Two-pass mapping recovers an
additional ~5–10 % of junction reads at the cost of roughly twice the
runtime. For RNA variant calling and for fusion detection, two-pass is
the right default.

=== HISAT2

HISAT2 (Kim, Paggi, Park, Bennett & Salzberg 2019) solves the same
problem as STAR with a different data structure. Instead of a
suffix array, HISAT2 uses a *hierarchical graph FM-index* (HGFM), which
extends the FM-index from Chapter 3 to a reference *graph* whose nodes
are genomic positions and whose edges include known splice junctions
and common variants from dbSNP.

The construction is layered.

- *Global index.* An FM-index over the whole reference genome, exactly
  as in Chapter 3. Handles the majority of within-exon alignment.
- *Local indexes.* Every ~57 kb of the genome has its own small
  FM-index, optimised for short extensions. A single read's alignment
  can hand off from the global index to a local index for the tricky
  last few bases without paying the cost of a full-genome backward
  search.
- *Graph edges.* Known splice junctions and common SNPs are encoded as
  alternative edges in the graph. The FM-index is extended to traverse
  these edges during search, so a read that crosses an annotated
  junction sees a contiguous query path even though the underlying
  reference is discontinuous.

#figure(
  image("../../diagrams/lecture-05/05-hisat2-hgfm.svg", width: 92%),
  caption: [
    HISAT2's hierarchical graph FM-index. A global FM-index covers the
    whole genome; hundreds of local indexes give cheap short
    extensions; annotated splice junctions appear as graph edges that
    the FM-index traverses during backward search.
  ],
) <fig:hisat2>

The payoff is memory: HISAT2's human index is about 8 GB, roughly four
times smaller than STAR's. The tradeoff is somewhat less sensitivity to
*novel* splice junctions not in the annotation, although in practice
the gap is small for libraries from well-annotated species.

#note[
  HISAT2 is FM-index alignment on a directed acyclic graph rather than
  on a linear string. The same Last-First mapping and backward-search
  procedure as Chapter 3 works — the FM-index alphabet is extended so
  that transitions across graph edges appear as valid moves during the
  search. This is the same generalisation that returns in the pangenome
  alignment work of recent years (vg, minigraph), where the reference
  is a population-scale variation graph rather than a single haplotype.
]

=== From BAM to Gene Counts

Once you have a sorted, indexed BAM, counting is mechanical. Two tools
dominate. *featureCounts* (part of the Subread suite) is a fast C
implementation that takes a BAM plus a GTF annotation file and emits a
gene-by-sample count matrix. *htseq-count* is the older Python reference
implementation — slower, but its `intersection-strict` mode is the
semantic ground truth everyone eventually calibrates against.

Both tools must make a decision for each read about which feature it
counts toward. The non-obvious cases are the ones that bite people.

- A read overlapping two genes (common in compact regions or near gene
  boundaries) — counted toward both, toward neither, or toward an
  "ambiguous" bucket?
- A read partly in an exon and partly in an intron — retained intron
  (biologically meaningful) or genomic-DNA contamination?
- A multi-mapping read (one primary plus several supplementary
  alignments) — split the count fractionally, assign to the primary
  alignment only, or drop the read entirely?

Every tool has slightly different defaults. The discipline: document
which tool and which counting mode produced the matrix, because
differential-expression results are not invariant to that choice.

#warn[
  RNA-seq libraries are usually *strand-specific* — the protocol
  preserves the strand of the original mRNA, and you can tell which
  strand each read came from. Strandedness affects which reads count
  toward which gene: a read on the forward strand that overlaps a gene
  encoded on the reverse strand should *not* count toward that gene.
  Running counting tools with `--stranded=no` on a stranded library
  silently double-counts antisense reads; running with
  `--stranded=reverse` on an unstranded library silently drops half the
  signal. Always verify strandedness with `infer_experiment.py` from
  RSeQC before counting. Strandedness misconfiguration is the most
  common bug in RNA-seq pipelines, and the only reliable symptom is a
  count matrix that doesn't match a public-data reference.
]


== Pseudoalignment <sec:pseudo>

Here is the observation that reorganised RNA-seq quantification in 2016
and ate the previous decade's tooling.

If your goal is transcript *abundance* — not base-level alignment
coordinates — you don't need to know *where* in each transcript a read
came from. You only need to know *which transcripts the read is
compatible with.* A read's "compatibility class" is the set of
transcripts whose sequence is consistent with the read.

The construction is straightforward. Pick a k-mer length (k = 31 is
typical). Extract every k-mer from every transcript in your reference
transcriptome and build a hash table that maps each k-mer to the set of
transcripts containing it. Some k-mers will be unique to one transcript;
some will appear in two or three closely related isoforms; some will
appear in dozens (think rRNA, or a paralog family). When a new read
arrives, chop it into k-mers, look each one up, intersect the resulting
transcript sets, and the intersection is the read's compatibility class.

#figure(
  image("../../diagrams/lecture-05/06-compatibility-classes.svg", width: 92%),
  caption: [
    Compatibility classes. Each k-mer carries a tag listing the set of
    transcripts that contain it; a read's class is the intersection of
    tags across all its k-mers. No Smith-Waterman, no seed extension,
    no base-level alignment.
  ],
) <fig:compat>

That is it. No Smith-Waterman, no seed extension, no base-level
alignment. The read's contribution to quantification is entirely
captured by the set of transcripts in its compatibility class. The
ambiguity in that class — a read that lands in a shared exon will have
a class of size three or four — is what the EM step of §5.4 resolves.

#tip[
  Pseudoalignment is sketch-based set membership. The transcriptome
  index is a hash-addressable summary of which transcripts contain
  which k-mers; a read is a query against that summary; the answer is a
  subset of the transcriptome universe. The data-flow is the same
  pattern as Bloom filters, MinHash, and HyperLogLog — query-time set
  operations against precomputed sketches. The difference from
  Chapter 3 indexing is the *task*: Chapter 3's indexes had to support
  coordinate retrieval, which is expensive. Pseudoalignment indexes
  support set membership only, which is cheap. Dropping the coordinate
  requirement is exactly what makes pseudoalignment 10–100× faster.
]

What you give up is exact base-level alignment coordinates, which means
you cannot directly genotype from RNA-seq reads with a pseudoaligner.
What you keep is everything you need for quantification: which
transcripts each read is compatible with, and how those compatibilities
distribute across the transcriptome.

=== Kallisto and the Target de Bruijn Graph

Kallisto (Bray, Pimentel, Melsted & Pachter 2016) was the first
pseudoalignment tool to ship. Its data structure is a *Target de Bruijn
Graph* (T-DBG): a de Bruijn graph (Chapter 7 builds these up
properly) constructed from the transcriptome, with each k-mer node
tagged by the set of transcripts that contain it.

For each incoming read, Kallisto walks its k-mers along the T-DBG and
intersects the transcript-set tags along the walk. The intersection
gives the compatibility class. A read whose walk produces an empty
class — no transcripts compatible with every observed k-mer — is
discarded as contamination or sequencing artifact.

#figure(
  image("../../diagrams/lecture-05/07-tdbg-kallisto.svg", width: 92%),
  caption: [
    The Target de Bruijn Graph. A read's k-mers trace a path through
    the graph; the intersection of transcript-set tags along that path
    is the compatibility class. Walk-and-intersect replaces alignment.
  ],
) <fig:tdbg>

The construction is linear in the transcriptome size (a few minutes for
a human transcriptome at k = 31). The query is amortised $O("read length")$.
A Kallisto quantification of 30 million RNA-seq reads against the human
transcriptome runs in under ten minutes on a laptop — an order of
magnitude faster than any aligner-based tool.

#note[
  Kallisto's release in May 2016 was a culture shock. Before Kallisto,
  RNA-seq quantification took CPU-hours per sample: STAR plus RSEM
  routinely ran for half a day on a 30 million-read library. Kallisto
  did the same job in minutes, with comparable accuracy on standard
  benchmarks. The paper argued — with hard benchmark numbers — that
  the explicit base-level alignment RSEM was spending all its time on
  was essentially wasted computation for the quantification task. The
  community took about two years to digest the argument. By 2018,
  pseudoalignment was the default.
]

=== Salmon and Selective Alignment

Salmon (Patro, Duggal, Love, Irizarry & Kingsford 2017) is
pseudoalignment with two refinements. The first, called *quasi-mapping*
in the original paper, replaced Kallisto's T-DBG with an enhanced suffix
array that explicitly tracked the approximate position of each k-mer
match within the transcript. Knowing position lets Salmon model
positional biases — reads do not come from a uniform distribution along
a transcript, because of 3′ bias from polyA selection, 5′ bias from
degradation, and position-dependent GC biases. Salmon corrects all
three explicitly in its likelihood model.

The second, introduced in Salmon 1.0 (2019), is *selective alignment*: a
hybrid mode where fast k-mer matching identifies candidate transcripts
and a fast, ungapped, banded Smith-Waterman alignment is then run only
on those candidates to produce a proper score. The result is the speed
of pseudoalignment with the accuracy of true alignment on the cases
where the difference matters — primarily, distinguishing reads from
highly similar paralogs.

The net effect: Salmon 1.0 on its default settings is *nearly* as fast
as Kallisto, *nearly* as accurate as STAR + RSEM, and handles
positional, GC, and sequence-context biases explicitly. For bulk
RNA-seq quantification in 2024 it is the pragmatic default. Kallisto
remains the right tool when speed is paramount (single-sample
turnaround on a laptop, very large cohort runs) and the bias model
matters less.

=== Pseudoalignment vs Alignment: What Each Gives Up

Pseudoalignment is not strictly better than alignment; it gives up some
information and keeps the rest. The tradeoffs:

- *Speed.* 10–100× faster than alignment-based quantification. A 30M
  read sample quantifies in 5–10 minutes on a laptop versus 4–8 CPU
  hours for STAR + RSEM.
- *Memory.* About 3–4 GB resident for a human transcriptome index, vs
  30 GB for STAR's genome index.
- *No BAM.* You cannot run downstream tools that need read coordinates
  — variant calling, splice-junction discovery, coverage tracks for
  visualisation. If you need a BAM for any reason, run the aligner.
- *Comparable accuracy on standard benchmarks.* Within 1–2 % of
  alignment-based tools on gene-level counts; 2–5 % on transcript-level
  counts; the gap is mostly in highly paralogous regions and in
  reads that span multiple junctions.

The decision tree most labs follow: if the only downstream goal is
differential expression, use Salmon. If you also need RNA variant
calling, splice-junction discovery, fusion detection, or coverage
plots, use STAR (often two-pass) plus `featureCounts`, and let Salmon
or RSEM consume the STAR alignments only if you have a specific reason.


== Expectation-Maximisation for Read Assignment <sec:em>

Compatibility classes are inherently ambiguous. A read that lands in
an exon shared by three transcripts has compatibility class
${T_1, T_2, T_3}$ — we know the read came from one of those three but
not which. The EM algorithm resolves this probabilistically.

EM was first formalised by Dempster, Laird, and Rubin in 1977 as a
general technique for maximum-likelihood estimation when some of the
data is missing. In RNA-seq the "missing data" is the true
read-to-transcript assignment: we observe only the compatibility class,
not the source transcript. The treatment that applies EM specifically
to transcript quantification is in Li and Dewey's RSEM paper (2011),
which Salmon and Kallisto's quantification step inherits almost
unchanged.

The model. Assume each read was generated in two steps: first pick a
transcript with probability proportional to its abundance times its
length (because longer transcripts produce more reads at the same molar
abundance), then pick a position within the transcript. The likelihood
of observing read $r$ with compatibility class $C(r)$, given transcript
abundances $theta = (theta_1, dots, theta_n)$, is

$ P(r | theta) = sum_(t in C(r)) (theta_t / ell_t) / (sum_(t') theta_(t') / ell_(t')) $

where $ell_t$ is the length of transcript $t$. Summing $log P(r | theta)$
across all reads gives the full data log-likelihood $L(theta)$.

We want to maximise $L(theta)$ over the simplex $theta_t >= 0$,
$sum_t theta_t = 1$. The log-sum structure couples the transcript
abundances, so no closed-form solution exists. EM handles the coupling
exactly by alternating two steps.

*E-step.* Compute the posterior probability that each read came from
each compatible transcript, given the current estimate of abundances:

$ z_(r,t)^((k+1)) = (theta_t^((k)) / ell_t) / (sum_(t' in C(r)) theta_(t')^((k)) / ell_(t')) $

The $z_(r,t)$ values are *fractional votes* — real numbers between 0 and
1, summing to 1 over the read's compatibility class. A read with a
unique compatible transcript votes 1.0 for that transcript and 0 for
everything else.

*M-step.* Update the abundance estimate by aggregating the fractional
votes and renormalising:

$ theta_t^((k+1)) = (sum_r z_(r,t)^((k+1))) / N $

where $N$ is the total number of reads. Iterate. Each iteration is
linear in the number of reads times the average compatibility class
size — typically a millisecond per iteration on a modern laptop.
Convergence to within $10^(-3)$ relative tolerance typically takes
100–1000 iterations, which is why a Salmon or Kallisto run can finish
in minutes.

#figure(
  image("../../diagrams/lecture-05/08-em-iterations.svg", width: 92%),
  caption: [
    Five EM iterations on a small example: the fractional-vote matrix
    redistributes each E-step, the abundance bars update each M-step,
    and the log-likelihood climbs monotonically to convergence.
  ],
) <fig:em-iterations>

To make the iteration concrete, @fig:em-worked carries the full
numerical evolution for a five-read, three-transcript toy problem.
Transcript T1 is much shorter than the others (500 vs 1500 vs 2000 nt),
but its high abundance is anchored by two reads that uniquely match it
(r5) or are dominated by it (r3, after a few iterations). Watch how
read r3, which starts the iteration with its vote split three ways,
ends up assigning 82 % of its weight to T1 by iteration ten — the
collective evidence from the unique-mapping reads has tilted the
posterior toward T1, and the ambiguous read follows.

#figure(
  image("../figures/ch05/f2-em-worked-iterations.svg", width: 96%),
  caption: [
    A full numerical EM walk for five reads, three transcripts of
    lengths 500, 1500, and 2000 nt, starting from a uniform abundance
    prior. By iteration ten the relative changes in $theta$ are below
    $10^(-3)$ and the algorithm has converged. Note that the
    most-ambiguous read (r3) eventually puts more than 80 % of its
    vote on T1 — its compatibility class is informed by the
    unambiguous reads.
  ],
) <fig:em-worked>

#tip[
  EM for transcript quantification is iterative soft-decision decoding.
  Each read is a noisy observation of a codeword (the transcript it
  came from); the compatibility class is the set of codewords
  consistent with the observation; the fractional votes are soft bits
  — probabilistic guesses at which codeword was sent. Iterating those
  soft bits against a model of codeword probabilities is exactly what
  a turbo decoder does on a convolutional code. The same algorithmic
  pattern returns in Chapter 8 for the variational EM used by scVI on
  single-cell counts, and it sits at the heart of every
  detection-theory curriculum.
]

A few practical points about the EM in production. Salmon and Kallisto
run a few thousand iterations of *bootstrapped* EM, resampling the
reads with replacement to produce a posterior distribution over
abundances per transcript rather than a single point estimate. The
bootstrap variance is what downstream tools like `sleuth` and
`tximport`'s offset machinery use as a per-gene uncertainty estimate
when feeding the negative-binomial GLM of Chapter 6. Don't throw it
away.

=== RSEM and the Alignment-Based Alternative

RSEM (Li & Dewey 2011) is the alignment-based predecessor to Salmon and
Kallisto. It runs Bowtie to produce base-level alignments against the
transcriptome and then applies the same EM as Salmon does to those
alignments. The conceptual logic is identical; the difference is where
compatibility classes come from — from actual Smith-Waterman alignments
rather than from k-mer intersections.

RSEM is slower (CPU-hours per sample) but produces base-level
information you can audit. For 2024+ work, Salmon or Kallisto handle 99
percent of bulk-RNA-seq use cases. RSEM remains in use when someone
already has STAR alignments and wants quantification without re-running
the alignment step, or when paranoid scrutiny of the alignment quality
matters more than runtime.


== Counts, Normalisation, and What the Units Mean <sec:norm>

After alignment-plus-counting or pseudoalignment, the output is a
*count matrix*: rows are genes (or transcripts), columns are samples,
and each cell is the number of reads assigned to that gene in that
sample. The count is not expression. It is expression × library depth
× transcript length × technical biases. To compare anything to anything
— gene to gene within a sample, the same gene across samples, or your
samples to a published reference — you need to strip out the systematic
factors and keep the biological signal.

Three categories of bias matter.

- *Library depth.* A sample sequenced to 50 M reads produces counts
  about twice as large as the same sample sequenced to 25 M. Any
  cross-sample comparison must correct for depth.
- *Transcript length.* At equal molar abundance, a 10 kb transcript
  collects ten times more reads than a 1 kb transcript. Any cross-gene
  comparison within or across samples must correct for length.
- *Library composition.* If one highly-expressed gene soaks up 50 % of
  reads in sample A but only 10 % in sample B, every *other* gene in
  sample A appears deflated relative to sample B even when its absolute
  expression is identical. This is the composition bias the
  median-of-ratios size factor of DESeq2 corrects for, and it is *not*
  the same as depth.

=== CPM, FPKM, TPM

Three units dominate practice. They look superficially interchangeable
and are not.

*CPM (counts per million)* normalises for library depth only:

$ "CPM"_(s, g) = c_(s,g) / (D_s / 10^6) $

where $c_(s,g)$ is the raw count of gene $g$ in sample $s$ and $D_s$ is
the total library depth. CPM makes samples of different depth
comparable gene by gene, but it does *not* correct for transcript
length — a long gene still has a larger CPM than a short gene at the
same molar abundance.

*FPKM (fragments per kilobase per million)* normalises for both depth
and gene length, in the order depth then length:

$ "FPKM"_(s,g) = c_(s,g) / ((D_s / 10^6) (ell_g / 1000)) $

FPKM is CPM divided by gene length in kilobases. RPKM is the same
quantity for single-end reads; FPKM is the paired-end version, where
the fragment (not the read) is the unit of accounting.

*TPM (transcripts per million)* normalises for length *first*, then for
depth:

$ "TPM"_(s,g) = (c_(s,g) / ell_g) / (sum_(g') c_(s,g') / ell_(g')) dot 10^6 $

The key difference: TPM divides each gene by its length before
normalising the sum across genes; FPKM normalises by total depth
before dividing by length. The order of operations matters. TPM values
sum to exactly $10^6$ across genes in every sample — they are a
proportion. FPKM values sum to some sample-dependent number that
depends on the gene-length composition of the library — they are
*not* a proportion and cannot be directly averaged across samples.

#figure(
  image("../../diagrams/lecture-05/09-normalization-units.svg", width: 92%),
  caption: [
    CPM, FPKM, and TPM on a small worked matrix. The three units rank
    genes differently depending on transcript length and sample
    composition.
  ],
) <fig:norm-units>

The arithmetic that makes the difference concrete is in
@fig:tpm-fpkm-arith. Four genes of different lengths and counts produce
four different rankings under CPM, FPKM, and TPM — and the rank flip
between gene B and gene D is large enough to mislead a casual reader of
the table.

#figure(
  image("../figures/ch05/f3-tpm-fpkm-cpm-arithmetic.svg", width: 96%),
  caption: [
    Four genes, the same input counts, three different normalisation
    units. CPM ranks Gene C highest because it has the most reads;
    TPM ranks Gene A highest once length is corrected. The unit you
    report is part of the result.
  ],
) <fig:tpm-fpkm-arith>

#warn[
  TPM is the right unit for *reporting* expression because it is a
  proportion and behaves cleanly across samples. CPM is the right unit
  to feed into differential expression tools like edgeR and limma-voom,
  because those tools model length as a regression offset internally
  and want depth-corrected counts as input. FPKM is legacy. Never
  compare a gene's FPKM between two samples without converting to TPM
  first — the result can be off by 50 % from composition effects alone.
]

=== Size Factors and Composition Bias

Depth normalisation by total library size assumes two samples with
different total counts differ only in depth. That assumption is usually
false in practice. Suppose sample A is a blood sample with hemoglobin
transcripts soaking up 40 % of reads, and sample B is a similar blood
sample but with globin depletion performed (a common protocol). Both
are sequenced to the same nominal depth. In sample A the globins eat 40
% of reads; every other gene's CPM is deflated. In sample B the
globins eat 2 %; the other genes' CPMs look 40 % higher. Just comparing
CPMs makes every non-globin gene appear up-regulated in B when in truth
the underlying expression is identical.

The DESeq2 fix, due to Anders and Huber (2010), is the *median-of-ratios*
size-factor estimator. The procedure is four steps.

+ For each gene $g$, compute the geometric mean of its counts across
  all samples in the experiment.
+ For each sample $s$, compute the ratio $c_(s,g) /
  "geom_mean"_g$ for every gene with non-zero count.
+ The sample's size factor $s_s$ is the *median* of those ratios
  across genes.
+ Divide every count in that sample by $s_s$ to get
  "size-factor-adjusted counts."

#figure(
  image("../../diagrams/lecture-05/10-size-factors.svg", width: 92%),
  caption: [
    DESeq2's median-of-ratios size factor in five panels — raw counts,
    geometric mean per gene, per-sample ratios, the median that
    becomes the size factor, and the normalised counts.
  ],
) <fig:size-factors>

The reasoning: if a sample is at population-average depth, the
distribution of its per-gene ratios is centred on 1. If its depth is
2× average, the ratios centre on 2. The *median* of the ratios is
robust to a handful of very highly (or very lowly) expressed genes
hijacking the estimator — exactly the failure mode that breaks naive
depth scaling on a sample with one runaway transcript. The estimator
behaves cleanly even when 30 % of reads come from a small number of
genes, which is the operating regime of every blood sample, every
muscle sample, and every cancer line.

#note[
  Median-of-ratios is the same robust-estimator move that returns
  throughout statistics — trimmed means in robust regression, the
  median absolute deviation as a scale estimator, the Hodges-Lehmann
  median of pairwise differences in nonparametric tests. The story is
  always the same: a small number of outliers shifts the mean
  arbitrarily but the median budges only when half the data moves.
  RNA-seq libraries always have outliers (the top 100 genes routinely
  account for 30–50 % of reads), so the right estimator is robust by
  construction.
]

edgeR uses a closely related method called *TMM* (trimmed mean of
M-values, Robinson & Oshlack 2010) that additionally trims out genes
with extreme log-fold-changes before averaging. For bulk RNA-seq with
standard biological replicates, DESeq2's median-of-ratios and edgeR's
TMM give nearly identical size factors; the differences appear only at
extreme library composition.

A reader who only takes one corollary from this section should take this:
*size factors correct for composition; they do not correct for gene
length*. You need both corrections, just for different questions. TPM
handles length-and-depth for intra-sample comparison; size factors
handle composition-and-depth for inter-sample comparison; the
differential-expression model in Chapter 6 handles all three together as
regression offsets and intercepts.


== The Negative Binomial <sec:nb>

Chapter 1 introduced the Poisson distribution as the natural model for
read counts: if reads are sampled uniformly at random from a long
transcriptome, the count of reads on any given feature is
Poisson-distributed with mean proportional to the feature's molar
abundance. For technical replicates — same RNA, re-sequenced — the
Poisson fits beautifully. For biological replicates — different mice,
different patients, different culture days — it doesn't.

The reason is variance. A Poisson distribution has the property that
$"var"[X] = mu$ when $X tilde "Poisson"(mu)$. Empirically, RNA-seq
biological-replicate data have variance much *larger* than the mean,
especially at high counts. The extra variance comes from biological
sources — cell-to-cell expression heterogeneity, subtle protocol
variation between library preps, slight chemistry drift between
sequencing runs. None of it is sampling noise; all of it is real
biology that the Poisson model has no way to express.

The fix is the *negative binomial* (NB) distribution. The NB has two
parameters, a mean $mu$ and a dispersion $alpha$, with

$ "var"[X] = mu + alpha mu^2 $

The $alpha mu^2$ term is the overdispersion: at $alpha = 0$ the NB
collapses to the Poisson; at $alpha > 0$ it has heavier tails. The NB
can be derived two equivalent ways: as a Poisson whose rate parameter
is itself drawn from a gamma distribution, or as a count of independent
Bernoulli failures before the $r$-th success. Both derivations land at
the same density; differential-expression tools mostly use the
gamma-Poisson form because the gamma prior is a natural way to model
biological variation in rate.

For real RNA-seq data, $alpha$ depends on the gene's mean expression.
Highly expressed genes have lower dispersion (the law of large numbers
helps); lowly expressed genes have higher dispersion (small counts are
noisier). Typical values are $alpha approx 0.01$ for technical
replicates, $alpha approx 0.1$ for clean biological replicates of cell
lines, and $alpha approx 0.3$ or above for clinical samples with
realistic between-patient variation.

#figure(
  image("../../diagrams/lecture-05/11-poisson-vs-nb.svg", width: 92%),
  caption: [
    Poisson and Negative Binomial at the same mean. Top: the two
    distributions side by side; the NB has a much longer right tail.
    Bottom: empirical mean-variance scatter from real biological
    replicates, sitting systematically above the Poisson line.
  ],
) <fig:nb-poisson>

@fig:nb-meanvar puts the empirical evidence on log-log axes. Each
point is one gene, plotted as its sample-wise mean count against its
sample-wise variance. The Poisson identity line is the diagonal.
Almost every gene with a non-trivial mean sits *above* the line; the
deviation grows with the mean exactly as $alpha mu^2$ predicts. The
three overlaid NB curves bracket the bulk of the data.

#figure(
  image("../figures/ch05/f4-meanvar-poisson-nb.svg", width: 96%),
  caption: [
    Empirical mean-variance scatter on log-log axes. The Poisson
    identity line lies systematically below the data; three negative
    binomial fits at dispersions 0.01, 0.1, and 0.3 fan out from the
    line as the mean grows. At $mu = 100$ a Poisson predicts standard
    deviation 10; real samples typically show 30 to 60.
  ],
) <fig:nb-meanvar>

#tip[
  The negative binomial is photon counting with multiplicative
  intensity noise. A pure Poisson process counts independent arrivals
  at a known rate; an NB process counts arrivals at a rate that itself
  fluctuates. The same shape appears in photon-limited imaging where
  the illumination flickers, in queuing theory with bursty arrivals,
  and in over-dispersed packet-loss models in networking. The
  dispersion parameter $alpha$ is exactly the relative variance of the
  intensity. Treating $alpha$ as zero — that is, using Poisson — gives
  spurious significance everywhere the data is heavier than the Poisson
  predicts.
]

Practically, the negative binomial is the noise model that DESeq2,
edgeR, and every modern differential-expression tool assume. Chapter 6
spends a section on how to estimate $alpha$ from a handful of biological
replicates (it shrinks toward a global mean using empirical-Bayes
inference, which trades a little bias for a lot of variance reduction).
The thing to take away from this chapter is just the shape: variance
grows faster than the mean, the NB captures that, and pretending data
is Poisson when it isn't inflates type-I error by a factor of three to
ten.


== The RNA-seq Variant Zoo <sec:variants>

Bulk RNA-seq measures one quantity: transcript abundance, averaged over
the cells in the tube. But cells regulate gene expression through
multiple layers — transcription rate, splicing choice, RNA-protein
binding, RNA modification, translation rate, RNA degradation. Each
layer has its own dedicated sequencing assay. They all reuse the
alignment, counting, normalisation, and NB-modelling machinery you have
just learned. What changes is the library-prep step upstream of the
sequencer and the biological interpretation of the resulting counts.

#figure(
  image("../../diagrams/lecture-05/12-ribo-seq.svg", width: 95%),
  caption: [
    The RNA-seq variant family as multi-channel signal acquisition.
    Bulk RNA-seq covers the whole transcript; ribo-seq footprints
    concentrate on protein-coding regions; CLIP-seq peaks mark
    RNA-binding-protein sites. Each channel uses the same downstream
    pipeline; joint analysis (e.g., translational efficiency =
    ribo / RNA) extracts derived quantities no single channel sees.
  ],
) <fig:ribo>

*Ribo-seq* (ribosome profiling, Ingolia, Ghaemmaghami, Newman & Weissman
2009) measures *translation*, not transcription. The trick: RNase-digest
a cell extract; everything is degraded except the ~28-nucleotide
footprints of mRNA that translating ribosomes were protecting at the
moment of lysis; sequence those footprints. The result is per-codon
ribosome occupancy across the transcriptome. The killer derived
quantity is *translational efficiency*:

$ "TE"_g = ("ribo-seq counts on gene " g) / ("RNA-seq counts on gene " g) $

A gene with high mRNA abundance but low ribosome occupancy is
transcribed but not translated — invisible to RNA-seq alone. TE is the
single most informative diagnostic for translational control. Stress
responses, viral infection, and the hallmark cancer translation programs
all reshape TE without corresponding mRNA changes. Tools: Plastid,
Ribotaper, RiboCode for ORF discovery; xtail and ribodiff for
differential TE testing.

*CLIP-seq* (cross-linking immunoprecipitation, Ule et al. 2003; iCLIP
and eCLIP are the modern variants) measures *RNA-protein binding*.
UV-crosslink an RNA-binding protein (RBP) to whatever RNAs it is
gripping; immunoprecipitate the RBP; sequence the protected RNA
fragments. The result is a transcriptome-wide footprint of one
specific RBP's binding sites. The ENCODE eCLIP atlas contains tracks
for about 150 human RBPs in two cell lines, including all the major
splicing factors (SF1, U2AF1/2, SRSF family), 3′ UTR regulators (HuR,
AUF1), and m^6^A readers/writers (METTL3, YTHDF1-3). Tools: PureCLIP,
omniCLIP, and CLIPper for peak calling.

Other variants in the family worth recognising:

- *m^6^A-seq / MeRIP-seq.* Maps N6-methyladenosine — the most
  abundant internal mRNA modification — at peak resolution. The
  derived "m^6^A epitranscriptome" is now a routine layer in stress and
  cancer biology.
- *CAGE* and *5′-RACE.* Map transcription start sites at single-base
  resolution; the FANTOM5 atlas of human CAGE data is the canonical
  reference for promoter usage.
- *3′-end seq* (PAS-Seq, QuantSeq). Map polyadenylation sites; cheaper
  than full RNA-seq for differential-expression studies where you don't
  need transcript-body coverage.
- *NET-seq / GRO-seq / PRO-seq.* Map active RNA Polymerase II
  positions; measure nascent transcription rather than steady-state
  mRNA.
- *SLAM-seq / TimeLapse-seq.* Metabolic labelling with a thiol-tagged
  uridine, then chemical conversion that turns labelled positions into
  apparent T-to-C mismatches. The mismatch rate per gene is the
  proportion of recently-synthesised RNA. Distinguishes synthesis from
  degradation rates without time-course sampling.

#tip[
  Treat the assay family as multi-channel signal acquisition. Each
  protocol gives you a different channel on the same underlying RNA
  biology; the downstream stack is identical (align, count, normalise,
  NB-test); the biology lives in what the upstream protocol selected
  for. Joint analyses — TE = ribo / RNA, m^6^A peaks vs CLIP peaks vs
  RNA-seq differential expression — extract quantities that no single
  channel can produce. The same multi-channel logic generalises to
  ATAC + RNA + ChIP joint inference in Chapter 9 and to multimodal
  single-cell in Chapter 8.
]


== Which Tool, When <sec:which-tool>

The RNA-seq tool landscape is rich, but the choice narrows once the
scientific goal is on the table.

For *standard bulk DE on gene-level counts*: Salmon → tximport → DESeq2
is the right default in 2024. Fast, accurate, handles biases, integrates
cleanly with the downstream stack. If you don't know what to use, use
that.

For *the same plus a BAM*: STAR (two-pass) → `featureCounts` → DESeq2.
Slower, more memory, but produces alignments you can use for variant
calling, fusion detection, splice-junction discovery, and coverage
plots.

For *transcript-level (isoform) quantification*: Salmon with
`--validateMappings`. The transcript-level output is what `tximport`
collapses to gene-level counts for DE; if isoform-level DE is the
question, hand the transcript counts directly to `sleuth` or to DRIMSeq.

For *long-read RNA-seq* (Iso-Seq, ONT direct RNA): Minimap2 in spliced
mode plus IsoQuant or Bambu. The long-read tools handle full-length
transcript identification natively; pseudoalignment doesn't apply.

For *fusion transcript detection in cancer*: STAR-Fusion or Arriba.
Both need STAR alignments as input; both look for chimeric reads
spanning gene-gene junctions.

For *RNA variant calling*: STAR two-pass → GATK HaplotypeCaller in RNA
mode. Calls germline variants from RNA-seq, with the caveats that you
only see variants in expressed transcripts and that allele-specific
expression biases the apparent VAF.

For *quick-look quantification on a laptop with no compute budget*:
Kallisto. Still hard to beat for raw throughput per CPU-second.


== Summary <sec:summary>

- RNA differs from DNA at the alignment step: splicing means reads
  don't map contiguously, and isoforms mean reads from one gene don't
  always pin to one transcript. The rest of the chapter follows from
  that one fact.
- Splice-aware aligners (STAR, HISAT2) extend the seed-and-extend
  framework of Chapter 3 with long-skip transitions at canonical
  GT-AG splice motifs. Both are hidden-Markov-model-flavoured; the
  state space is the only thing that changes.
- Pseudoalignment dispenses with base-level alignment for the
  quantification task. A read's compatibility class — the set of
  transcripts whose k-mers are consistent with the read — is sufficient
  for EM-based quantification, and the speedup is 10–100×.
- The EM algorithm resolves the ambiguity in compatibility classes by
  iterating fractional read-to-transcript assignments to convergence.
  The shape is iterative soft-decision decoding; the math is identical
  to RSEM's 2011 formulation.
- Normalisation has three layers. CPM corrects for depth; TPM corrects
  for length and depth; the DESeq2 size factor (median-of-ratios)
  corrects for composition. Different questions need different layers.
- Counts are overdispersed Poisson. The negative-binomial dispersion
  parameter $alpha$ captures excess variance from biological sources.
  Modelling counts as Poisson — pretending $alpha = 0$ — inflates
  apparent significance by an order of magnitude on real data.
- The same downstream pipeline serves ribo-seq, CLIP-seq, m^6^A-seq,
  CAGE, and the rest of the RNA-seq variant zoo. The differences live
  in what the library prep selected for upstream of the sequencer.


== Exercises <sec:exercises>

#strong[1.] #emph[Splice motif recognition.]
A 150-bp paired-end read on chromosome 7 has its first 80 bp matching
positions 117,479,920–117,479,999 in the reference and its last 70 bp
matching positions 117,484,580–117,484,649. The 4,580-bp gap on the
reference between these two regions is bounded by the dinucleotides
`GT…AG`. Should a splice-aware aligner report this as one alignment
with a splice, or as two alignments? Justify in one sentence with
reference to the canonical-motif rule. If the gap dinucleotides had
been `CA…TT` instead, what would change?

#strong[2.] #emph[Compatibility class arithmetic.]
You have three transcripts $T_1$, $T_2$, $T_3$ of equal length 1 kb,
sharing exons as follows: $T_1$ and $T_2$ share their first 500 bp;
$T_2$ and $T_3$ share their last 500 bp; $T_1$ and $T_3$ have no
overlap. Five reads arrive, each 100 bp long, sampled at positions
(a) 250 of $T_1$; (b) 750 of $T_1$; (c) 250 of $T_2$; (d) 750 of $T_2$;
(e) 250 of $T_3$. For each read, write down its compatibility class.
Which read is unambiguous?

#strong[3.] #emph[Hand EM.]
Implement EM for the toy problem in @fig:em-worked: three transcripts
of length $(500, 1500, 2000)$ nt and five reads with compatibility
classes $({T_1, T_2}, {T_2}, {T_1, T_2, T_3}, {T_3}, {T_1})$.
Initialise $theta = (1"/"3, 1"/"3, 1"/"3)$, run ten iterations, and
report the converged $theta$. Verify against the figure. What single
read carries the most information about the relative abundances of
$T_1$ and $T_2$, and why?

#strong[4.] #emph[CPM/TPM rank flip.]
Build a small four-gene count matrix that produces a different rank
ordering under CPM than under TPM, using only gene lengths in the
range 0.5–10 kb and counts in the range 1,000–5,000. Verify the
arithmetic by hand. Which gene's apparent rank changes the most, and
which normalisation does the lay reader most often default to?

#strong[5.] #emph[Median-of-ratios size factor.]
A four-sample, six-gene count matrix has the following row
geometric means: 10, 100, 1000, 10000, 50, 500. Sample 1's raw counts
for those genes are (8, 90, 1100, 9500, 60, 540). Compute the
per-gene ratios, then the median, then the size factor. Now suppose
gene 4's count in sample 1 was actually 95,000 instead of 9,500 (a
clerical mis-recording inflating the highest-expressed gene tenfold).
Recompute the size factor. By how much does the inflated outlier shift
the size factor? Compare to what would happen if you used the *mean*
ratio instead of the median.

#strong[6.] #emph[Mean-variance check.]
Pretend you have a count matrix from six biological replicates of a
control condition, with 18,000 genes. For each gene compute the
sample mean and sample variance across replicates. Plot the
mean-variance scatter on log-log axes; overlay the Poisson line
(slope 1, intercept 0). For genes with mean 100, what fraction of
points lies above the Poisson line in a typical bulk RNA-seq dataset?
What does that tell you about using a Poisson GLM directly for
differential expression?

#strong[7.] #emph[The strandedness check.]
You have a BAM from a bulk RNA-seq library and you don't know whether
the library was stranded. Describe in one paragraph the procedure
`infer_experiment.py` follows: which features of the read alignments
tell you the strand protocol, what the three possible answers are
(`unstranded`, `forward`, `reverse`), and what the consequence is of
running `featureCounts` with the wrong setting on a gene encoded on
the negative strand.

#strong[8.] #emph[(Open-ended.)]
Pick one tool from the chapter — Salmon, Kallisto, STAR, HISAT2, RSEM,
or DESeq2's size-factor estimation — and read its primary publication.
In one paragraph, describe the single most surprising design decision
the authors made and why it works on the empirical data they show. Pay
particular attention to anything in the paper that an EE student
recognises as a well-known technique from a different domain.


== Further Reading <sec:further-reading>

- *Dobin, A., Davis, C. A., Schlesinger, F., et al.* (2013). "STAR:
  Ultrafast Universal RNA-seq Aligner." _Bioinformatics_ 29: 15–21.
  The canonical STAR paper.

- *Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L.*
  (2019). "Graph-based Genome Alignment and Genotyping with HISAT2
  and HISAT-genotype." _Nature Biotechnology_ 37: 907–915. HISAT2
  with the graph FM-index.

- *Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L.* (2016).
  "Near-optimal Probabilistic RNA-seq Quantification." _Nature
  Biotechnology_ 34: 525–527. The Kallisto paper that introduced
  pseudoalignment.

- *Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford,
  C.* (2017). "Salmon Provides Fast and Bias-aware Quantification of
  Transcript Expression." _Nature Methods_ 14: 417–419. Salmon and
  the bias-aware quantification model.

- *Li, B., & Dewey, C. N.* (2011). "RSEM: Accurate Transcript
  Quantification from RNA-Seq Data with or without a Reference Genome."
  _BMC Bioinformatics_ 12: 323. The EM-for-RNA-seq paper.

- *Dempster, A. P., Laird, N. M., & Rubin, D. B.* (1977). "Maximum
  Likelihood from Incomplete Data via the EM Algorithm." _Journal of
  the Royal Statistical Society B_ 39: 1–38. The original EM paper.

- *Anders, S., & Huber, W.* (2010). "Differential Expression Analysis
  for Sequence Count Data." _Genome Biology_ 11: R106. The
  median-of-ratios size factor and the negative-binomial GLM.

- *Robinson, M. D., & Oshlack, A.* (2010). "A Scaling Normalization
  Method for Differential Expression Analysis of RNA-seq Data."
  _Genome Biology_ 11: R25. The TMM normalisation used in edgeR.

- *Ingolia, N. T., Ghaemmaghami, S., Newman, J. R. S., & Weissman, J.
  S.* (2009). "Genome-Wide Analysis in Vivo of Translation with
  Nucleotide Resolution Using Ribosome Profiling." _Science_ 324:
  218–223. The ribo-seq paper.

- *Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L.*
  (2016). "Transcript-level Expression Analysis of RNA-seq
  Experiments with HISAT, StringTie and Ballgown." _Nature Protocols_
  11: 1650–1667. A worked-end-to-end-pipeline reference.
