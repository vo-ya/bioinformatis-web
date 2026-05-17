#import "../theme/book-theme.typ": *

= Long Reads and the #idx("pangenome")Pangenome: From One #idx("STRING")String to a Graph <ch:long-reads>

#matters[
  Two pieces of bioinformatics infrastructure quietly changed between
  2019 and 2024, and almost every downstream pipeline in this book now
  inherits both. The first change is that long-read sequencing crossed
  the accuracy threshold where it is no longer just "long but noisy" — a
  #idx("PacBio")PacBio #idx("HiFi")HiFi read is now as accurate as a Sanger read, sixty times its
  length, at one-thousandth the cost; #idx("Oxford Nanopore")Oxford Nanopore has tracked
  alongside and is now within reach of the same quadrant. The second
  change is that the reference genome stopped being a string. The
  Human Pangenome Reference Consortium's first release, in May 2023,
  replaced #idx("GRCh38")GRCh38 with a graph of forty-seven haplotypes that encodes
  variation as branching structure rather than as a column in a #idx("VCF")VCF.
  Either change on its own would be infrastructural; together they
  retire most of the assumptions Chapters 2 through 4 were forced to
  make. This chapter is the formal reset.
]

For twenty years a genome assembly was a string, and a sequencing
project was a #idx("translation")translation problem: every read had to be placed on that
string, every difference reported as a delta. The string was good
enough when reads were short and variation was small, which is to say,
for almost everything covered in Chapters 2 through 4. The two
discontinuities described above broke that frame from both sides — reads
long enough to span whole structural events, and a reference graph
expressive enough to represent variation as first-class topology rather
than as second-class annotation. The result is a different stack: long
reads at the input, a directed acyclic graph as the reference, graph
alignment as the new map operation, and Viterbi-on-a-DAG as the math
that makes all of it tractable.

This chapter walks the stack from input to output. Section 11.1 takes
stock of where the long-read platforms sit in 2024 and explains why
crossing the Q20 line was the change worth paying #idx("attention")attention to. Section
11.2 catalogues what long reads make easy — structural variants, tandem
repeats, #idx("haplotype")haplotype #idx("phasing")#idx("phasing")phasing, #idx("telomere-to-telomere")telomere-to-telomere assembly — each of
which was genuinely hard or actually impossible with short reads alone.
Sections 11.3 and 11.4 introduce the pangenome graph and the file
format that made it portable. Sections 11.5 and 11.6 work through the
algorithms: #idx("seed-and-extend")seed-and-extend generalised to a DAG, #idx("Viterbi")Viterbi dynamic
programming on a DAG, and the complexity argument for when the graph
formulation pays for itself. Section 11.7 covers #idx("T2T-CHM13")T2T-CHM13 and #idx("HPRC")HPRC v1
as the two concrete artifacts of this transition, and Section 11.8
treats phasing — already mentioned in passing — as a stand-alone
source-separation problem.

The chapter assumes you understand sequencing platforms at the level of
Chapter 1, short-read alignment at the level of Chapter 2, and variant
calling at the level of Chapter 4. The two #idx("structural variant")structural variant sections
of Chapter 4 (§4.7 in particular) are the most-load-bearing reference
points; the algorithmic machinery from Chapter 2 (#idx("k-mer")k-mer indexing, seed
chaining, banded #idx("dynamic programming")dynamic programming) is the spine of Section 11.5.


== Long-Read Tech in 2024 <sec:long-reads-tech>

Chapter 1 sketched the long-read platforms in their pre-2019 form, as a
counterweight to #idx("Illumina")Illumina's reign. PacBio's single-molecule real-time
(#idx("SMRT")SMRT) instrument watched a #idx("polymerase")polymerase work in a zero-mode waveguide;
Oxford Nanopore threaded #idx("DNA")DNA through a protein pore and watched the
ionic current modulate. Read lengths reached ten to a hundred
kilobases; per-base error rates, in the early platforms, sat between
five and fifteen percent. The mental model the field operated under was
simple and durable: long reads are long but noisy; short reads are
short but accurate; pick one. That mental model is now obsolete.

=== HiFi and the Crossing of Q30

PacBio's *HiFi* (High-Fidelity) chemistry, released in 2019 alongside
the work of Wenger and colleagues, made one decisive change. The
template molecule was circularised, the polymerase ran around it eight
to fifteen times in a single SMRT well, and the per-pass reads were
combined into a circular-consensus read. Each pass had its own
~10–15 % error rate; the consensus across eight passes brought the
expected per-base error to roughly $10^(-3)$ — Q30. The 2024 generation
of HiFi, running on the Revio instrument, sits at Q33 with read
lengths around 20 kb. Two thousand gigabases per day per instrument;
ten to twenty US dollars per gigabase.

Q30 at 20 kb is a regime no sequencing platform had previously
occupied. #idx("Sanger sequencing")Sanger sequencing was Q30 at 800 bp. Illumina is Q30 at 150
bp. HiFi is Q30 at twenty thousand bases — Sanger-grade accuracy at one
hundred and twenty-fifth the cost per base, with reads that span things
no Sanger read or Illumina read could ever span.

=== Oxford Nanopore and the Q20+ Era

Oxford Nanopore took longer. The fundamental physics is different —
there is no consensus pass; each base is called once, from the
modulation it imparts on the pore current as it translocates through.
Accuracy depends entirely on the #idx("basecaller")basecaller's ability to decode that
modulation, and the basecaller is a neural network that has been
re-trained on every chemistry update.

Three chemistry generations carry the story. *R9.4* (deployed widely
2018–2021) sat at roughly Q12 raw — about 6 % per-base error — and
required heavy alignment-based consensus to be useful for anything
beyond structural variants and ultra-long #idx("contig")contig scaffolding. *R10.4*
(2022) moved the pore design to a dual-reader chemistry that sampled
the current at two points along the molecule, jumping raw accuracy to
Q18–Q20. *R10.4.1 + Q20+* (2023, paired with the modern Dorado
basecaller) crossed the Q20 line for raw reads, with Q25+ after
polishing. Read lengths range from 10 kb to over a megabase; the
"ultra-long" tail of the distribution is what makes #idx("ONT")ONT uniquely
capable for centromeric assembly.

ONT is not identical to HiFi. The error profile is structured —
specific k-mer contexts basecall worse than others, homopolymer runs
are still where most errors concentrate, and the per-base confidence
score from the basecaller is less calibrated than HiFi's. But the gap
that mattered for downstream tools has closed.

#figure(
  image("../figures/ch11/f1-quality-length-quadrants.svg", width: 95%),
  caption: [
    The read-length × per-base-accuracy plane, with platform positions
    plotted for 2018 and 2024. The "reference-grade" quadrant — length
    at least 10 kb, Q at least 20 — was empty until 2019. HiFi and
    modern ONT are now both in it.
  ],
) <fig:quadrants>

#note[
  The Q20 line is not arbitrary. Below it, every read is a noisy
  vote and #idx("variant calling")variant calling requires deep alignment-based consensus —
  many reads voting on each base. Above it, each read is a
  mostly-correct witness; small-variant calling becomes feasible from
  a single read, methylation-modification calling becomes feasible at
  per-read resolution, and phasing becomes a problem of grouping
  internally consistent reads rather than of inferring linkage from
  population statistics. Most of the analytical techniques in this
  chapter assume Q20+ data.
]

=== EE Framing: Read Length as Observation-Window Size

Worth one paragraph as an aside. In signal processing, the
detectability of a transient event scales with the observation window
relative to the event's duration. Events shorter than the window are
detectable by matched-filter correlation; events longer than the
window are partially observable at best, and recovering them requires
inference over multiple weak partial views. The same trade applies in
genomics: a structural variant shorter than the read length is
resolved directly by a single alignment, a structural variant longer
than the read length is visible only through indirect inference on
fragment geometry — split reads, discordant pair orientations, #idx("coverage")coverage
anomalies. The 50 bp–10 kb regime, essentially invisible at 150 bp
read length, becomes trivial at 20 kb read length. This is the
biophysical analogue of pulse-compression in radar: enlarge the
observation window and previously sub-threshold events come into the
detection regime.

#warn[
  Long reads still have error modes. HiFi's residual errors
  concentrate in homopolymer runs and at sites where the polymerase
  stutters; ONT's residual errors are bias-structured, with certain
  k-mers reproducibly basecalled worse than others. The "reads are
  perfect now" overclaim breaks at predictable loci. Use a modern
  long-read-aware variant caller — #idx("DeepVariant")DeepVariant in long-read mode,
  Clair3, or PEPPER-Margin-DeepVariant — rather than a short-read
  caller pointed at a #idx("BAM")BAM of long reads.
]


== What Long Reads Unlock <sec:unlock>

The story of the long-read transition is best told through the four
problem classes that flipped from hard or impossible under short reads
to routine under long. Each class had its own short-read workaround —
indirect #idx("SV")SV inference, PCR fragment sizing, statistical phasing, gap
masking — and each workaround is now obsolete.

=== Structural Variants, Resolved Directly

Chapter 4 covered structural variant calling under short reads
(@sec:sv-recap below). The recap: a 150 bp short read cannot span
anything longer than itself, so SVs are visible only through the three
indirect channels Chapter 4 catalogued — #idx("discordant read")discordant read pairs straddling
the breakpoint, split reads cut across it, and depth anomalies inside
deleted or duplicated regions. Sensitivity below 60 % for mid-size SVs
was typical even on high-coverage short-read whole-genome data; precision
was a constant negotiation between false positives in repeat-rich regions
and false negatives near segmental duplications.

Long reads change the inference problem to an alignment problem. A
single 20 kb HiFi read that spans an entire 5 kb tandem duplication
shows the duplication directly inside its own alignment — the gap, the
repeat, the flanks. No inference is needed; one or two reads at the
locus are enough to call the event.

#figure(
  image("../../diagrams/lecture-11/02-sv-long-vs-short.svg", width: 95%),
  caption: [
    Short reads can only infer a structural variant from indirect
    evidence — discordant pairs, split reads, coverage shifts. A single
    long read spans the entire event.
  ],
) <fig:sv-vs>

The dominant 2024 long-read SV callers are *Sniffles2*, *CuteSV*, and
*Severus*. All three operate on the same principle: scan the alignments
for read-internal evidence of insertions, deletions, inversions, and
translocations; cluster supporting reads; emit a VCF with single-base
breakpoint coordinates and inserted sequence where it can be assembled.
On the #idx("GIAB")GIAB HG002 benchmark for SVs between 50 bp and 10 kb, modern
long-read pipelines report both precision and recall above 95 %.
Short-read pipelines on the same sample sit closer to 70–80 %.

=== Spanning Repeats

A subset of tandem-repeat loci are clinically important because the
number of repeat units determines disease. *Huntington's disease* —
CAG repeats in `HTT`, #idx("pathogenic")pathogenic at ≥ 36 units, severity scaling with
expansion size. *Fragile X syndrome* — CGG repeats in `FMR1`,
pathogenic at ≥ 200 units, with mosaic premutation alleles in the
55–200 range. *Myotonic dystrophy* — CTG repeats in `DMPK`, pathogenic
above ~50 units. Two dozen more disorders follow the same pattern: a
short #idx("motif")motif, repeated many times, expanding across generations, with
disease determined by repeat count.

Short reads cannot count repeat units when the expansion exceeds the
read length. Clinical labs sized these alleles with PCR-based
fragment-length assays for thirty years because sequencing could not
do it. Long reads sidestep the problem entirely: a single 3 kb HiFi
read can span 100 CAG repeats plus a hundred bases of flanking
sequence on each side, and the units are counted directly. Tools —
*TRGT* (PacBio), *Straglr*, *LongTR*, *NanoRepeat* — wrap this into
clinical-grade genotype calls.

#figure(
  image("../../diagrams/lecture-11/03-tandem-repeats.svg", width: 95%),
  caption: [
    Long reads count tandem-repeat units directly by spanning the
    expansion. The HTT CAG locus, pathogenic above 36 units, is the
    canonical clinical example.
  ],
) <fig:repeats>

The same principle resolves segmental duplications (multi-kilobase
blocks duplicated at greater than 90 % identity across the genome) and
the long alpha-satellite arrays at centromeres. These were the
regions short-read pipelines masked out by default. Long reads
de-mask them.

=== Phasing from a Single Molecule

A diploid human carries two near-identical copies of each autosome,
one from each parent. At every heterozygous site, the sample has
both alleles; at every homozygous site, only one. *Phasing* is the
problem of deciding which heterozygous alleles share a parental
#idx("chromosome")chromosome — that is, recovering the two parental haplotypes from the
mixture of reads.

Short-read phasing is statistical. A 150 bp read overlaps at most one
or two heterozygous sites, and the combinations across many such reads
do not uniquely determine phase. The fix is to bring in a population
panel — 1000 Genomes, HGDP, #idx("gnomAD")gnomAD — and use *#idx("linkage disequilibrium")linkage disequilibrium*
(LD), the population-level correlation between nearby variants, as the
prior. The result is probabilistic, accurate for common variants in
well-sampled populations, and progressively worse for rare variants
and for populations under-represented in the reference panel.

Long reads break the ambiguity directly. A 20 kb HiFi read covers
five to ten heterozygous sites on average; an ultra-long ONT read
covers dozens. The phase among those sites is read off the single
molecule — no population, no statistics, no prior. The block length
of long-read phasing is set by the read length convolved with
coverage; for 30× HiFi, chromosome-scale phase blocks are routine.

#figure(
  image("../../diagrams/lecture-11/04-phasing-short-vs-long.svg", width: 95%),
  caption: [
    Short reads cover one or two heterozygous sites at a time; phase
    requires inference against a population LD panel. A single long
    read carries phase across all the het sites it spans.
  ],
) <fig:phasing>

We will treat phasing in more depth in @sec:phasing, including the
ML-natural framing of phasing as supervised source separation.

=== Telomere-to-Telomere Assemblies

The 2004 human-genome-project assembly was widely reported as
complete; it had known gaps. GRCh38, released in 2013, still
contained about 200 Mb of gaps and about 100 Mb of unplaced contigs.
Both gap sets concentrated in the same regions: centromeric
alpha-satellite arrays, ribosomal-DNA clusters, acrocentric short
arms, and segmental duplications. Short reads could not assemble these
regions because the repeats were longer than the reads.

The first actually complete human-genome assembly arrived in April
2022, when the *#idx("T2T")T2T Consortium* published the *T2T-CHM13* assembly
(Nurk et al.). It used the hydatidiform-mole cell line CHM13 — a cell
line that carries two copies of a single paternal genome and is
therefore effectively haploid, side-stepping the phasing problem — and
combined ultra-long ONT reads (200× coverage, with a sub-population
above 1 Mb), HiFi data at 40× for accurate consensus, optical mapping
for long-range scaffolding, and extensive manual curation. The
assemblers *Verkko* and *hifiasm* did most of the algorithmic work.

#figure(
  image("../../diagrams/lecture-11/05-t2t-before-after.svg", width: 95%),
  caption: [
    A chromosome before and after T2T. GRCh38 left the #idx("centromere")centromere and
    the acrocentric arm as a masked gap; T2T-CHM13 fills both in with
    resolved alpha-satellite higher-order repeats.
  ],
) <fig:t2t>

T2T-CHM13 added roughly 200 Mb of sequence to GRCh38, concentrated in
~70 Mb of centromeric alpha-satellite, ~40 Mb of acrocentric short
arms (with their rDNA gene clusters), and ~50 Mb of segmental
duplications. Concrete biology fell out immediately: the first
per-chromosome catalogue of centromere dynamics, complete rDNA arrays
that finally enabled rRNA-biology research that had been limited to
gene-copy counting, and resolved sequence for medically important
duplicated genes like `SMN1`/`SMN2` (spinal muscular atrophy) and the
`NBPF` family. The *T2T-Y* assembly followed in 2023, using a separate
donor, and the field now treats *T2T-CHM13 v2.0* (incorporating Y) as a
plausible default reference for production pipelines.


== The Linear-Reference Problem <sec:reference-bias>

T2T-CHM13 closed the gaps in the reference, but it did not solve the
deeper problem with using a single reference at all. GRCh38 — and
T2T-CHM13 v2.0, in the same way — is one specific combination of
alleles drawn from a small number of donors. Humans differ from each
other by about one base in a thousand at the single-nucleotide level
and by half a percent to one percent at the structural-variant level.
Every variant call against a single reference reports a delta against
*one combination*, and the choice of which combination matters.

=== #idx("reference bias")Reference Bias

The pathology is called *reference bias*. Variants common in the
population but absent from the reference are systematically
under-called. A read carrying a non-reference allele has a higher
probability of misalignment (it has fewer matching bases against the
reference around the variant site); it lands in low-mapping-quality
territory; the variant caller weights its evidence lower or drops it
entirely. The bias compounds at structural variants — a 5 kb insertion
common in a population but absent in the reference fundamentally cannot
be aligned, only inferred.

The 2019 work of Sherman and colleagues quantified this on African
sequencing data. Aligning African short-read whole-genome data to a
GRCh38 augmented with 296 Mb of African-specific assembled sequence
recovered thousands of variants that the canonical pipeline had called
as "no-call" or homozygous reference. The clinical consequences are
real and one-directional. A pathogenic variant common in the patient's
population but absent in GRCh38 maps as a false reference match — the
variant is invisible at the alignment step, before any caller sees it.
Polygenic risk scores derived on European cohorts transfer poorly to
non-European populations partly for the same reason. Ancestry inference
under-counts variation in under-represented populations.

#figure(
  image("../../diagrams/lecture-11/06-reference-bias.svg", width: 95%),
  caption: [
    Three variants at a single locus, each common in a different
    ancestry. GRCh38 carries V1; reads carrying V2 or V3 misalign or
    drop. A pangenome graph that includes all three calls them
    correctly for every ancestry.
  ],
) <fig:bias>

#note[
  Reference bias is a training–test distribution mismatch.
  Population genetics has independently re-discovered every conclusion
  the machine-learning literature reached about domain shift: a model
  trained on one distribution and applied to a different one degrades
  in systematic, predictable ways, with the magnitude of degradation
  proportional to the distance between the two distributions. The fix
  in both fields is the same — broaden the training distribution to
  cover the deployment regime. In genomics the broadening takes the
  form of a pangenome graph.
]

=== The Pangenome Graph

A *pangenome graph* — or *graph genome*, the terms are used
interchangeably — represents variation as branching topology in a
directed graph rather than as differences from one linear backbone.
Each node carries a short DNA segment; each edge joins two segments
that co-occur in at least one known haplotype. A *path* through the
graph is one specific haplotype: pick a sequence of edges, concatenate
the segment sequences along the way, and you have one person's chromosome.

At a biallelic #idx("SNP")SNP, the graph has two nodes — one carrying the
reference base, one carrying the alternate — both bracketed by the
shared flanking sequence. At a small #idx("indel")indel, one branch contains the
inserted sequence and another skips it. At a common structural
variant, the graph has two or more branches representing the SV
haplotypes. Common variation is structure, not annotation.

#figure(
  image("../../diagrams/lecture-11/07-linear-vs-graph.svg", width: 95%),
  caption: [
    Same locus, two data structures. The linear reference has
    variants annotated off to one side; the graph has them as
    first-class branches with explicit alternate paths.
  ],
) <fig:linear-vs-graph>

A useful order-of-magnitude calibration. The HPRC v1 pangenome encodes
roughly 110 Gb of sequence — the 47 input haplotypes laid end to end —
in a graph that compresses to about 4 Gb of #idx("GFA")GFA text. The
compression is not lossy: it works because most of the genome is
shared across haplotypes, and the graph only branches where they
disagree. Where two haplotypes are identical the graph runs as a
single chain of nodes; where they diverge it branches; where they
re-converge it merges. The whole structure is a compact representation
of all 47 sequences simultaneously.

=== EE Framing: From Delta-Encoding to Dictionary Encoding

A linear reference plus a per-sample VCF is delta encoding: one
baseline sequence, every other genome described as differences. A
pangenome graph plus per-sample path assignments is dictionary
encoding: the full set of observed haplotype fragments stored in the
codebook, and each sample's chromosome encoded as the indices of the
fragments it uses. The trade-off is the same as in every entropy-coding
scheme. Delta encoding is compact when the deltas are small and the
baseline is representative; it fails when the baseline is not
representative, because the encoder has to spell out long runs of
non-baseline material as "differences." Dictionary encoding is larger
in raw bytes but robust to baseline shift, because every codeword
appears in the dictionary regardless of which sample is being encoded.
Pangenomes apply the dictionary-encoding shift to reference genomes.


== GFA: The Pangenome File Format <sec:gfa>

The text format that carries a pangenome from disk into memory is
*Graphical Fragment Assembly* — GFA. It is tab-separated, one record
per line, with a one-character line-type code at the start of each
line. The grammar is small enough to fit on one page.

```
H  VN:Z:1.0                          # header, version
S  1  ACGTAAGTTTG                    # segment: node ID=1, sequence
S  2  ACGCAAGAATG                    # segment: alternate allele node
S  3  CGATCGATCGA                    # segment: shared flanking sequence
L  1  +  3  +  0M                    # link: node 1 forward → node 3 forward
L  2  +  3  +  0M                    # link: node 2 forward → node 3 forward
P  hap1  1+,3+  *                    # path: haplotype 1 = 1+ then 3+
P  hap2  2+,3+  *                    # path: haplotype 2 = 2+ then 3+
W  sample1  1  chr1  0  22  >1>3     # walk: sample's chromosome through the graph
```

Five record types, all you really need to know:

- `H` — header, with format version and optional metadata.
- `S` — segment, the graph node: an integer ID plus a DNA sequence.
- `L` — link, the graph edge: from-node, from-orientation, to-node,
  to-orientation, plus an overlap #idx("CIGAR")CIGAR (`0M` for no overlap).
- `P` — path, a named haplotype expressed as a comma-separated list of
  oriented nodes.
- `W` — walk, the GFA 1.1 successor to `P`, expressing a sample's
  chromosome walk in a Bandage-compatible angle-bracket notation.

#figure(
  image("../../diagrams/lecture-11/08-gfa-anatomy.svg", width: 95%),
  caption: [
    GFA text and the graph it encodes, side by side. Each line type
    corresponds to one element of the graph: header, node, edge, path,
    walk.
  ],
) <fig:gfa>

#tip[
  The orientation flags on segments and edges (the `+` and `-` after
  each node ID) matter. DNA is double-stranded, and a graph node has
  a forward sequence and its reverse complement. A path that
  traverses `1+,3+` is not the same as `1-,3+`: the first emits the
  forward sequence of node 1, the second emits its reverse
  complement. Most early GFA parsing bugs trace to ignoring orientation.
]

The toolchain around GFA is reasonably mature. *#idx("vg")vg* (the variation
graph toolkit; Garrison et al., 2018) is the reference implementation,
with subcommands for construction, alignment, calling, and projection
back to linear. *minigraph* (Li, 2021) is faster and uses a simpler
graph model that prohibits nested variation. *pggb* — the Pangenome
Graph Builder — constructs a graph from a collection of input genomes
without requiring a designated reference. *odgi* provides
graph-manipulation and analysis primitives, in the same role `samtools`
plays for BAMs.

#warn[
  A pangenome graph is only as diverse as the haplotypes it was
  built from. HPRC v1 includes 47 haplotypes from a deliberately
  diverse panel, but a variant carried by no HPRC sample remains
  unrepresented in the graph; reads carrying that variant will still
  misalign or get called as reference. Graph genomes shrink
  reference bias substantially; they do not eliminate it. The right
  intervention for a specific clinical or population-genetic study
  is to build a study-specific pangenome that includes samples from
  the target population.
]


== Graph-Based Seed-and-Extend <sec:seed-extend>

The algorithmic move that lets reads align to a graph is a careful
generalisation of Chapter 2's seed-and-extend. Every step of the
short-read alignment recipe has to be rethought when the reference is
no longer a line.

Chapter 2's linear seed-and-extend:

+ Pick k-mers from the read.
+ Look up each k-mer in the reference's k-mer index to find candidate
  alignment positions.
+ For each candidate, extend the alignment base-by-base using dynamic
  programming (Smith–Waterman or banded variants).
+ Report the best-scoring alignment.

On a graph reference, every step changes shape:

+ *K-mer indexing on a graph.* Each node carries a sequence, but k-mers
  that span node boundaries — k-mers that overlap an edge — have to be
  indexed as "k-mer plus path-prefix" rather than "k-mer plus position".
  Practical implementations enumerate all length-k paths through the
  graph and index each one. The index size grows with the branching
  factor, but for human pangenomes the growth is manageable.
+ *Seed placement* gives a starting node (not just a position) and an
  offset within the node's sequence.
+ *Extension* is no longer along a single chain. At each step, the
  extension can continue into any out-edge of the current node, and
  the algorithm has to consider all of them. The local DP table
  acquires a third dimension: read position, node, and offset within
  the node.
+ *Scoring* becomes the maximum over all graph paths of length matching
  the read, not just over the linear neighbourhood.

#figure(
  image("../../diagrams/lecture-11/09-graph-seed-extend.svg", width: 95%),
  caption: [
    Seed-and-extend on a DAG. The seed locates a starting node;
    extension branches whenever the graph branches and picks the
    best-scoring continuation.
  ],
) <fig:graph-seed>

=== The vg / minigraph / GraphAligner Toolchain

Three tools dominate the space, separated by their algorithmic choices
more than by their feature sets.

*vg* (Garrison et al., 2018) is the reference implementation, with
full #idx("graph alignment")graph alignment under principled scoring. Its main subcommands
trace the entire pipeline: `vg construct` builds a graph from a
reference FASTA plus a VCF of known variants; `vg giraffe` is the
fast short-read graph aligner; `vg map` is the slower but more
accurate general-purpose aligner; `vg deconstruct` pulls variants from
graph alignments; `vg surject` projects graph alignments back onto a
linear reference for downstream tools that expect linear BAMs. vg is
about 5–10× slower than `bwa-mem` on the same data; tractable, not
free.

*minigraph* (Li, 2021) trades expressiveness for speed. It uses a
simpler graph model that prohibits nested variation (variants inside
variants — for example, an SNP inside an inserted sequence) and
borrows its alignment core from `minimap2`. The result is much faster
construction and alignment, useful when building pangenome scaffolds
from large numbers of input genomes.

*GraphAligner* (Rautiainen & Marschall, 2020) is the long-read-specific
counterpart. It uses a `minimap2`-derived seed–chain–extend with
long-read-aware gap penalties and is the standard choice for aligning
HiFi or ONT reads to a pangenome graph.

The 2024 pipeline shape:

+ Build the graph with `pggb` or `minigraph cactus` from a collection of
  reference-quality whole genomes.
+ Align reads with `vg giraffe` for short reads or `GraphAligner` for
  long reads.
+ Call variants with `vg call` (or a graph-aware caller).
+ Project to linear with `vg surject` if downstream tools expect
  linear BAMs.

#figure(
  image("../../diagrams/lecture-11/10-vg-toolchain.svg", width: 95%),
  caption: [
    The vg toolchain end-to-end. Build, align, call, optionally project
    to linear. File formats are labelled on each arrow.
  ],
) <fig:toolchain>

The practical trade-offs reduce to four numbers. Recall on
under-represented populations rises by 5–20 % on SVs and on common
non-reference variants — the bigger the gap between the linear donor and
the target sample, the bigger the gain. Runtime rises by 5–10×.
Tooling friction is non-trivial: many downstream tools still want
linear BAMs, and `vg surject` bridges that with edge cases that bite at
locus boundaries. And graph quality is bounded by the input panel:
generic pangenomes are a better default than a single linear reference;
study-specific pangenomes are better still for specific cohorts.

#tip[
  A graph alignment's output is a *path* through the graph — a list of
  (node, offset) tuples. Two ways to consume it: read variants
  directly off the path (each branch choice implicitly encodes which
  allele was traversed; `vg call` does this), or project the graph
  alignment onto a chosen linear reference so downstream tools can
  treat it as a BAM. Each choice has its own failure modes; the
  graph-native variant call is more honest, the projected BAM is
  more compatible.
]


== Viterbi on a DAG <sec:viterbi>

The cleanest piece of the whole long-reads-and-pangenome story is the
algorithm underneath graph alignment. It is not a new algorithm. It is
Viterbi dynamic programming on a state graph, exactly the algorithm an
EE student met in their first communications, control, or speech-
recognition course. The linear-chain alignment of Chapter 2 is the
special case where the state graph is a single line; graph alignment
is the general case where the state graph is a directed acyclic graph.

=== From Smith–Waterman to Viterbi

Chapter 2 introduced dynamic programming for sequence alignment as
Smith–Waterman and Needleman–Wunsch — a two-dimensional score matrix
$M[i, j]$, filled cell-by-cell with the recurrence

$ M[i, j] = max cases(M[i-1, j-1] + s(a_i, b_j), M[i-1, j] - g, M[i, j-1] - g, 0) $

where $a_i$ is the read base at position $i$, $b_j$ is the reference
base at position $j$, $s$ is the match/mismatch score, and $g$ is the
gap penalty. The traceback walks the maximum-scoring path back from
the highest-scoring cell.

The cleaner framing for the same algorithm is *Viterbi decoding* on a
trellis. The reference is a chain of states, one per reference base.
The read is observed over time, one character per time step. The
alignment is the maximum-likelihood state path through the trellis
given the noisy observation sequence. The score matrix is the Viterbi
forward accumulator; the traceback is the Viterbi traceback. Every
piece of Chapter 2's alignment machinery transfers directly.

The reframing matters because it immediately generalises. Viterbi
applies to any trellis whose state graph is a DAG, not just to linear
chains. Pangenome alignment is the natural generalisation: the state
graph stops being a line and becomes the pangenome DAG.

=== The Update Equation, Once

Formally, define $V(s, t)$ as the best alignment score ending at graph
state $s$ after consuming $t$ bases of the read. The update is

$ V(s, t) = max_(s' in "pred"(s)) [V(s', t - 1) + "score"(s, "read"[t])] $

with the maximum taken over all predecessors $s'$ of $s$ in the state
graph. At the end of the read, take the highest $V(s, T)$ over all
terminal states $s$ and trace back along the back-pointers to recover
the best path.

For a linear chain, $"pred"(s)$ is a single state — the state
immediately to the left of $s$ — and the algorithm reduces to standard
Smith–Waterman. For a DAG, $"pred"(s)$ is the set of nodes with
in-edges into $s$, and the algorithm handles branches naturally. The
update equation does not change. Only the structure of the
predecessor set does.

#figure(
  image("../figures/ch11/f2-viterbi-update.svg", width: 95%),
  caption: [
    Viterbi on a linear chain (left) versus on a DAG (right). The
    update equation is the same in both cases; only the predecessor
    set of each state differs.
  ],
) <fig:viterbi>

#note[
  The single cleanest genomics-as-signal-processing translation in
  the entire book lives here. Alignment is Viterbi decoding;
  pangenomes simply replace the state graph. The back-pointer logic,
  the optimality argument, the asymptotic complexity, the role of
  match/mismatch scores as log-likelihood-ratios — all of it transfers
  from the channel-decoding literature unchanged. If you understand
  the Viterbi algorithm on a trellis, you understand graph alignment.
]

=== Complexity, Correctness, and When the Graph Wins

The computational accounting follows immediately. Linear-chain Viterbi
runs in $O(N times T)$ time, where $N$ is the reference length and $T$
is the read length. DAG Viterbi runs in $O(|E| times T)$, where $|E|$
is the number of edges in the graph. For an HPRC-scale graph
(~10⁹ edges), this is 5–10× slower than the linear case on the same
input — meaningful but tractable on modern hardware.

Correctness is preserved. Viterbi on a DAG returns the true optimal
path through the graph, in the same sense Smith–Waterman returns the
true optimal alignment on a chain. The optimality property — that the
best path ending at state $s$ at time $t$ extends from the best path
ending at some predecessor at time $t - 1$ — holds for any DAG of
states. Graph alignment is not a heuristic approximation of linear
alignment; it is a strict generalisation that reduces exactly to the
linear case when the reference graph is a line.

The interesting question is when graph DP beats linear DP. For a read
that carries no non-reference content, the best graph path runs
through the reference-identity nodes and graph Viterbi collapses to
linear Viterbi: no penalty, no improvement, just extra runtime. For a
read that carries a variant *not in the linear reference but present
in the graph*, the difference is sharp. The linear alignment
accumulates a mismatch (or gap) penalty at the variant site and
misaligns the read or soft-clips the flanking sequence. The graph
alignment traverses the alternate-allele node at the same locus
without penalty. The score difference is the mismatch cost of the
variant — small for an SNP, large for a multi-kilobase SV. Graph
alignment's benefit concentrates at SV-rich loci and at population-
specific variants. Most of the runtime cost is paid on reads that
gain nothing; the gain is concentrated on the reads that matter.

#warn[
  A graph aligner that returns a path through the graph is harder to
  interpret than a linear aligner that returns a (chromosome,
  position, CIGAR) tuple. Downstream consumers of the alignment have
  to know how to read graph coordinates. Tooling that pre-dates the
  pangenome era handles this poorly; `vg surject` exists for exactly
  this reason, but it introduces edge cases at branch boundaries. A
  graph-aware pipeline that surjects at the end is doing more work
  than a linear pipeline; budget the engineering accordingly.
]


== T2T and HPRC: The State of the Art <sec:t2t-hprc>

The two concrete artifacts of the long-reads-and-pangenome transition
are *T2T-CHM13* (the complete linear reference) and *HPRC v1* (the
47-haplotype pangenome). The first is a single genome assembled to
the last base; the second is a graph of haplotypes that together
span enough human diversity to make graph alignment meaningfully
worthwhile.

=== T2T-CHM13 Revisited

@sec:unlock introduced T2T-CHM13 in passing. The breakdown of the
200 Mb of added sequence relative to GRCh38 is worth seeing in detail:

- *Centromeres* (~70 Mb). Alpha-satellite arrays of 171 bp repeats
  organised into higher-order repeats (HORs). CHM13 resolved the
  centromeric array on every chromosome.
- *Acrocentric short arms* (chr 13, 14, 15, 21, 22; ~40 Mb). Carry
  rDNA gene clusters and satellite DNA. Previously unresolvable in
  any assembly.
- *Segmental duplications* (~50 Mb). Multi-kilobase blocks duplicated
  at greater than 90 % identity. Includes regions important for
  immune function (`NBPF` genes, olfactory receptors) and several
  medically critical loci (`SMN1`/`SMN2`).
- *The Y chromosome* — not in CHM13 itself (CHM13 is XX), but the
  T2T-Y assembly published in 2023 used a separate male donor and the
  current production reference incorporates it as T2T-CHM13 v2.0.

The biological gain showed up immediately. The first per-chromosome
catalogue of centromere dynamics — alpha-satellite HOR variation,
centromere drift, satellite-array breakpoints. Complete rDNA arrays
that finally allowed rRNA-biology research to escape gene-copy
counting. Resolved sequence for `SMN1` and `SMN2` such that copy-number
genotyping of spinal-muscular-atrophy patients no longer required a
custom assay. Production pipelines have begun migrating from GRCh38 to
T2T-CHM13 v2.0; many clinical pipelines remain on GRCh38 for
compatibility, and the migration will take years.

=== HPRC v1 — Forty-Seven Haplotypes

The *Human Pangenome Reference Consortium* released its v1 pangenome
in May 2023 (Liao et al., _Nature_). The release covers 47
haplotype-resolved human genome assemblies — most diploid, from 24
individuals — chosen for ancestral diversity across six broad
groupings: African, European, East Asian, South Asian, Admixed
American, and Oceanian. Each haplotype was assembled with `hifiasm` in
trio mode, with parental short reads providing the k-mer signatures
that disambiguate haplotype origin. Most haplotypes are near-T2T for
most chromosomes.

The pangenome graph itself, built with `pggb`, encodes ~110 Gb of
underlying sequence (47 haplotype assemblies laid end to end) in a
graph of approximately 100 million nodes and 130 million edges. The
graph captures variation absent from GRCh38:

- ~50,000 structural variants (≥ 50 bp) common in at least one
  included population.
- Thousands of smaller-variant haplotype combinations.
- Explicit copy-number variation at segmental duplications and
  tandem repeats.

Downstream impact, measured on graph-aware variant callers:
approximately 10 % more variant calls recovered on diverse samples,
with the largest gain on under-represented ancestries.

#figure(
  image("../figures/ch11/f4-pangenome-recall-gain.svg", width: 95%),
  caption: [
    Recall gain of HPRC v1 pangenome alignment over linear GRCh38, by
    variant class and ancestry. Gains are largest for SVs and for
    populations farthest from the GRCh38 donor.
  ],
) <fig:recall-gain>

#figure(
  image("../../diagrams/lecture-11/12-hprc-pangenome.svg", width: 95%),
  caption: [
    HPRC v1: 47 haplotypes across six ancestry groups, compiled into a
    graph that explicitly represents variation that GRCh38 misses.
  ],
) <fig:hprc>

=== What's Still Missing

The reference is more complete than it has ever been, and it is still
incomplete. Forty-seven haplotypes covers a fraction of human
diversity; HPRC v2 (targeting 350 haplotypes) is expected mid-2025, v3
will likely push past 700. Rare variants specific to individual
families or small populations are still not in the graph and remain
invisible to graph-based calling. The HPRC pangenome is human-only;
pangenome graphs for model organisms — mouse, *Drosophila*,
*Arabidopsis*, yeast — are at varied stages of construction. And a
large fraction of downstream tools (variant annotators, structural-
variant filters, clinical-report generators) still expect linear
inputs, so the `surject` bridge stays load-bearing.

#note[
  The HPRC v1 release was widely framed as "the completion of the
  Human Genome Project". The framing is dramatised but substantially
  accurate. For the first time, the reference explicitly represented
  human variation as topology rather than as annotated differences
  against a single donor. The first move was the 2003 HGP draft. The
  second was T2T-CHM13's gap closure in 2022. The third is the
  pangenome — the recognition that even a perfect single reference is
  the wrong data structure.
]


== Phasing as Source Separation <sec:phasing>

@sec:unlock introduced phasing as one of the four problems long reads
made easy; this section treats it on its own terms and frames it as
the cleanest machine-learning analogy in the chapter.

=== The Problem Setup

A diploid genome carries two copies of every autosome. The two copies
are chemically indistinguishable in the lab — they separate in space
inside the #idx("nucleus")nucleus, but they go into the sequencing tube as a single
homogenised solution. At every heterozygous site, the sample carries
two different alleles: one on the maternal chromosome, one on the
paternal. Phasing is the problem of deciding which maternal allele
sits with which other maternal allele on the same molecule — that is,
of recovering the two parental haplotype sequences from the
post-mixing read pool.

A single heterozygous site cannot be phased without external
information; there is no "before" against which to align the maternal
and paternal allele. Two heterozygous sites on the same read *can* be
linked — the read directly shows which allele at site 1 co-occurs with
which allele at site 2. But short reads are short; typical haplotype
blocks from short-read phasing span only a few kilobases before
phase information runs out and a new block has to start.

Long reads extend the block length dramatically. A 20 kb HiFi read
links 5–10 heterozygous sites at typical human heterozygosity rates.
An ultra-long ONT read in the 100 kb–1 Mb range links dozens. The
phasing block is set by the joint distribution of read length and
coverage; for 30× HiFi, chromosome-scale phase blocks are routine.

=== Tools and the Outcome

Three tools dominate phasing in 2024:

- *hifiasm* (Cheng et al., 2021; 2024). The default HiFi assembler,
  with haplotype-resolved output: the assembly itself produces two
  contigs per chromosome, one per parental haplotype. In trio mode it
  uses parental short reads to label the child's HiFi reads by
  haplotype origin; in #idx("Hi-C")Hi-C mode it uses Hi-C chromatin-conformation
  scaffolding when trio data is unavailable.
- *WhatsHap* (Martin et al., 2016; v2 Garg et al., 2024). Takes a VCF
  of variant calls and a BAM of long reads and emits a phased VCF
  with each heterozygous variant labelled with its haplotype block.
  Standard for clinical data.
- *HapCUT2* (Edge et al., 2017). Graph-cut alternative to WhatsHap's
  min-cost-flow formulation. Similar performance characteristics;
  different internals.

The 2024 outcome is unremarkable: a phased VCF for any diploid human
sequenced to 30× HiFi is routine, with chromosome-scale blocks as the
default expectation. The barrier is no longer the algorithm.

=== EE Framing: Phasing Is Source Separation

The cleanest framing of phasing borrows from audio and communications.
Two underlying source signals — the maternal and paternal haplotypes —
are linearly superimposed in the observed mixture (the read pool). The
job is to recover each source.

In audio, this problem is the cocktail-party / blind-source-separation
problem: two speakers talking simultaneously, one microphone, recover
each voice. In communications, it is joint demodulation of a
multi-access channel. The standard solution toolbox — independent-
component analysis, beamforming, factorial HMMs — assumes you do not
know in advance which sample of the mixture came from which source.

Long-read phasing turns this from unsupervised into *supervised*
source separation. Each long read is a coherent block of the mixture
that has not been pre-mixed — within one read, every allele came from
the same parental chromosome. The label that an unsupervised separation
would have to infer is *carried inside the read itself* as the
combination of heterozygous-site genotypes. Reads with the same
combination of alleles are from the same source; reads with the
mirror combination are from the other source. The remaining residual
problem is to chain coherent blocks across reads whose haplotype
membership has to be resolved by the small overlap of heterozygous
sites between them.

#figure(
  image("../figures/ch11/f3-phasing-source-separation.svg", width: 95%),
  caption: [
    Phasing as supervised source separation. Two haplotype sources are
    superimposed in the read pool; each long read carries its source
    label inside itself by virtue of carrying multiple heterozygous
    sites from one parental chromosome.
  ],
) <fig:source-sep>

Short-read phasing is the *unsupervised* version of the same problem.
With only one or two het sites per read, no single read uniquely
identifies its source; the analysis falls back on the structural prior
that haplotypes drawn from the same population are correlated through
linkage disequilibrium, and that correlation lets a probabilistic model
recover most blocks most of the time, accurately for common variants
and progressively worse for rare ones. Long-read phasing is the
supervised version; short-read phasing is the unsupervised version
with LD as the structural prior. The genomics and ML literatures
agree on which is easier.

#warn[
  Phasing errors *cascade*. A phase switch — an error at heterozygous
  site $k$ in a chromosome-length block — inherits to every variant
  downstream of $k$ inside the same block. Modern long-read phasers
  report phase-switch error rates of 0.1–1 % per block, which sounds
  small until you remember that the downstream half of a long block
  is then sourced from the wrong parent. For clinical use, validate
  long-range phase against trio data or Hi-C; report unvalidated
  long-range phase with caveats.
]

== A Brief Note on Where the SV Section of Chapter 4 Stands Now <sec:sv-recap>

Chapter 4 §4.7 covered structural-variant calling under short reads —
the three-channel framework (discordant pairs, split reads, depth),
the tool zoo (Manta, DELLY, GRIDSS, LUMPY), the accuracy ceiling
around 70 % recall on mid-size SVs. Most of that section is still
correct as a description of short-read SV calling. None of it is
correct as a description of the SV calling that production pipelines
do in 2024 if they have access to long reads.

The replacement workflow is the one this chapter has built up to:
align HiFi or ONT reads with `minimap2` or `GraphAligner`; call SVs
with Sniffles2, CuteSV, or Severus; phase the calls with WhatsHap or
hifiasm trio. Recall jumps from ~70 % to above 95 % at mid-size SVs,
and the events that get caught include the events the short-read
pipeline never had a chance against — multi-kilobase repeats, complex
nested SVs, mobile-element insertions, segmental duplications.

The cost is incremental sequencing dollars. Short reads remain the
workhorse for general clinical WGS at around five hundred dollars per
genome; long reads add a thousand to twenty-five hundred dollars per
genome. For research-grade SV catalogues, rare-disease cases where
short-read analysis has stalled, and any cohort serious about
complete variant ascertainment, the cost difference is increasingly
the right side of the trade.


== Summary <sec:summary>

- Long reads in 2024 are *long and accurate*. PacBio HiFi sits around
  Q33 at 20 kb; modern ONT is Q20+ raw, Q25+ polished, also at 20 kb
  with an ultra-long tail. The reference-grade quadrant — length above
  10 kb, accuracy above Q20 — was empty before 2019; both platforms
  now live in it.
- Long reads make SVs, tandem repeats, phasing, and T2T assembly
  routine. Each of these was hard-or-impossible under short reads.
  The detection problem becomes the alignment problem.
- A single linear reference systematically under-represents
  populations distant from its donor. *Reference bias* has measurable
  clinical consequences — missed pathogenic variants in
  under-represented populations and #idx("PRS")PRS scores that transfer poorly
  across ancestries.
- A *pangenome graph* makes variation first-class topology rather
  than annotation. GFA is the standard text format; vg / minigraph /
  pggb is the toolchain. HPRC v1 (47 haplotypes, 2023) is the
  current human pangenome; v2 (350) is in flight.
- Graph alignment is *Viterbi on a DAG*. Chapter 2's linear-chain
  alignment is the special case where the state graph is a line. The
  update equation is unchanged; only the predecessor set of each
  state changes.
- Phasing is *source separation*. Long reads make it supervised
  (each read carries its source label inside itself by virtue of
  spanning multiple heterozygous sites). Short-read phasing is the
  unsupervised version with linkage disequilibrium as the structural
  prior.


== Exercises <sec:exercises>

#strong[1.] #emph[Quadrant calibration.]
Plot the position of a hypothetical sequencer that emits 5 kb reads at
Q22 onto the read-length × accuracy plane of @fig:quadrants. Does it
sit in the reference-grade quadrant? Justify in one sentence, then
list one analytical capability that this device unlocks compared with
a 5 kb read at Q12 and one that it still cannot match against 20 kb
HiFi.

#strong[2.] #emph[Repeat genotyping.]
A clinical lab has a sample with HTT CAG repeats. Short-read coverage
at 30× gives an estimate of "≥ 40 units" with an upper-bound
uncertainty. The same sample sequenced to 15× HiFi yields ten reads
that each fully span the repeat. The repeat-unit counts on the ten
reads are: 38, 39, 40, 41, 38, 39, 39, 40, 40, 39 (homozygous estimate)
and a second cluster of 76, 75, 77, 76, 77, 76, 77, 76 (heterozygous
estimate). What genotype do you report, and what is the most
informative summary statistic from this dataset for a clinical report?

#strong[3.] #emph[Viterbi update by hand.]
On a small DAG with five states $s_1, ..., s_5$ and edges
$s_1 -> s_2, s_1 -> s_3, s_2 -> s_4, s_3 -> s_4, s_4 -> s_5$, and a
read of length 5, write out the Viterbi update for $V(s_4, 4)$
explicitly as a `max` of two terms. Make explicit which terms come
from `pred`($s_4$). For a uniform match score of $+1$ and uniform
gap penalty of $0$, compute $V(s_5, 5)$ under any deterministic
choice you like for the read sequence.

#strong[4.] #emph[GFA decoding.]
Given the GFA file in @sec:gfa, list the two haplotypes encoded by
`P hap1` and `P hap2`. Re-write the same graph so that node 2 is
traversed in reverse complement orientation on `hap2`; show the
modified `P` line and explain why the resulting sequence is different
from the original.

#strong[5.] #emph[Recall-gain arithmetic.]
@fig:recall-gain reports an absolute recall gain of about 17
percentage points for mid-size SVs on African ancestry samples
(58 % → 75 %). For a clinical pipeline serving 10,000 African-ancestry
patients with average 4,000 structural variants in the 50 bp–10 kb
range per genome, estimate the expected number of additional SVs that
the pangenome pipeline would identify across the cohort relative to
the linear pipeline. Show the arithmetic.

#strong[6.] #emph[Phase-switch cascade.]
A chromosome-length phase block on chromosome 1 has 5,000 heterozygous
variants in a single 200 Mb block. The reported phase-switch error
rate is 0.5 % per block. (a) What fraction of variants downstream of
a single mid-block phase switch are mis-phased? (b) How does this
compare to per-variant Mendelian-error rates from short-read trio
phasing of common SNPs (~ 0.01–0.05 % per variant, accumulating across
the block independently)? Discuss in one paragraph which regime is
more dangerous for an analysis that joins multiple variants — for
example, a compound-heterozygous analysis.

#strong[7.] #emph[(Open-ended.)] Pick one of the three pangenome
tools (vg, minigraph, pggb) and read its primary publication. In one
paragraph, describe the single most surprising design choice the
authors made and the empirical observation that justifies it.


== Further Reading <sec:further-reading>

- *Wenger, A. M., Peluso, P., Rowell, W. J., et al.* (2019). "Accurate
  Circular Consensus Long-Read Sequencing Improves Variant Detection
  and Assembly of a Human Genome." _Nature Biotechnology_ 37:
  1155–1162. The HiFi launch paper.
- *Nurk, S., Koren, S., Rhie, A., et al.* (2022). "The Complete
  Sequence of a Human Genome." _Science_ 376: 44–53. T2T-CHM13.
- *Liao, W.-W., Asri, M., Ebler, J., et al.* (2023). "A Draft Human
  Pangenome Reference." _Nature_ 617: 312–324. HPRC v1.
- *Garrison, E., Sirén, J., Novak, A. M., et al.* (2018). "Variation
  Graph Toolkit Improves Read Mapping by Representing Genetic
  Variation in the Reference." _Nature Biotechnology_ 36: 875–879.
  The vg paper.
- *Li, H.* (2021). "minigraph as a Potential Reference-Free Pangenome
  Construction Approach." _bioRxiv_. The minigraph paper.
- *Garrison, E., Guarracino, A., Heumos, S., et al.* (2024). "Building
  Pangenome Graphs." _Nature Methods_ 21: 2008–2012. pggb.
- *Cheng, H., Concepcion, G. T., Feng, X., Zhang, H., and Li, H.*
  (2021). "Haplotype-Resolved #idx("de novo assembly")De Novo Assembly Using Phased Assembly
  Graphs with hifiasm." _Nature Methods_ 18: 170–175.
- *Sherman, R. M., Forman, J., Antonescu, V., et al.* (2019).
  "Assembly of a Pan-Genome from Deep Sequencing of 910 Humans of
  African Descent." _Nature Genetics_ 51: 30–35. The reference-bias
  paper that motivated much of the pangenome work.
- *Rautiainen, M., and Marschall, T.* (2020). "GraphAligner: Rapid and
  Versatile Sequence-to-Graph Alignment." _Genome Biology_ 21: 253.
- *Martin, M., Patterson, M., Garg, S., et al.* (2016). "WhatsHap:
  Fast and Accurate Read-Based Phasing." _bioRxiv_. The original
  WhatsHap paper; the 2024 Garg et al. update is the current
  reference.
- *GFA specification.* `github.com/GFA-spec/GFA-spec` The format. Read
  it before writing any pangenome-processing code.
