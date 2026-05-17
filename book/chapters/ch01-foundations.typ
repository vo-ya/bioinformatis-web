#import "../theme/book-theme.typ": *

= Foundations: From Cells to Sequences to FASTQ <ch:foundations>

#matters[
  Bioinformatics begins with a translation problem. A cell is a chemical
  machine that stores its operating manual in a polymer; an engineer's job
  is to read that polymer with an instrument and turn its output into a
  text file. Everything else in this book — alignment, assembly, variant
  calling, expression analysis, protein structure — assumes you understand
  what the polymer is, how the instrument reads it, and what the text file
  contains. We spend a chapter on it because most analysis bugs trace back
  to a misunderstanding of one of those three things.
]

A genome is not, strictly speaking, an information object. It is roughly
three picograms of phosphate-backboned nucleic acid coiled inside a
membrane-bound organelle, replicating itself once per cell cycle and
occasionally being transcribed by a polymerase that walks along it at
about fifty bases per second. Calling it _information_ is already a
modelling choice. It is the choice that makes bioinformatics possible.

This chapter walks the path that choice opens up. We start in the cell,
where the polymer physically lives. We work outward through the chemistry
of the double helix and the central dogma, then jump abruptly into
engineering: the sequencing instruments that take a tube of dissolved
DNA as input and emit a text file as output. By the end of the chapter
you will know what is inside that text file, why each character is
there, and what kind of errors you are signing up for when you open it.

The remaining twenty-six chapters of this book all start from that
text file. So it is worth getting right.


== Cells: Where the Information Lives <sec:cells>

The genome of every organism you will encounter in this book lives inside
a cell. Cells come in two architectures.

*Prokaryotes* — bacteria and archaea — keep their DNA loose in the
cytoplasm, organised as a single circular chromosome with no surrounding
membrane. The genome is small, typically one to ten million base pairs.
There is no separation between the place where DNA lives and the place
where it is translated into protein: a ribosome can begin reading an
mRNA before the polymerase has finished writing it.

*Eukaryotes* — everything else, including yeast, plants, fungi, and
animals — wrap their DNA in a nuclear membrane. The genome is large
(typically hundreds of millions to billions of base pairs), organised
across multiple linear chromosomes, and packaged at several levels:
DNA wraps around histone octamers to form nucleosomes; nucleosomes
fold into chromatin fibres; chromatin folds into chromosome territories
inside the nucleus. There is also a second, smaller genome inside every
eukaryotic cell — the mitochondrial genome — a 16-kilobase circular
chromosome that is a evolutionary souvenir from when an
ancestral cell engulfed an aerobic bacterium and never let go.

#figure(
  image("../../diagrams/lecture-01/01-cell-organelles.svg", width: 90%),
  caption: [
    Eukaryotic cell architecture and the locations where DNA, RNA, and
    protein synthesis happen. The nucleus holds the primary genome;
    mitochondria carry their own small circular chromosome.
  ],
) <fig:cell>

You will spend most of your bioinformatics career working with eukaryotic
genomes — human, mouse, plant, yeast. Bacterial genomes show up in
microbiome work (Chapter 23) and in basic-research contexts. When we say
_genome_ without qualification in this book, we mean the diploid nuclear
human genome: 3.05 billion base pairs across 22 autosomes, X, Y, and the
mitochondrion, in two near-identical copies per somatic cell.

#note[
  The "near-identical" matters. The two copies of each autosome differ
  at roughly four to five million positions — single nucleotide
  variants, small insertions and deletions, and occasional structural
  changes. These differences _are_ the variants that Chapter 4 will
  teach you to call.
]

=== A Short Detour: How We Learned the Polymer Was the Carrier

The recognition that DNA — and not protein — is the molecule that
carries genetic information took roughly seventy years to crystallise.
Friedrich Miescher isolated "nuclein" from pus-soaked bandages in 1869
and recognised it as a phosphorus-rich molecule from cell nuclei, but
neither he nor his contemporaries suspected its function. Through the
first half of the twentieth century the consensus assumption was that
chromosomes carried hereditary information in their _protein_ component;
DNA was assumed to be a structural scaffold, too chemically monotonous
to encode anything interesting.

Three experiments changed that. Avery, MacLeod, and McCarty in 1944
showed that the "transforming principle" that converted harmless
pneumococcal bacteria into virulent ones was DNA — not protein, not RNA.
Hershey and Chase in 1952 separated phage DNA from phage protein with
radioactive labels and showed that only the DNA entered infected
bacteria. And Watson, Crick, Franklin, and Wilkins in 1953 published
the double-helix structure, immediately suggesting a mechanism by which
the molecule could replicate itself.

The history matters for one reason: it is hard to overstate how
counter-intuitive the answer felt at the time. DNA has four monomers.
Proteins have twenty. From any nineteenth-century chemist's perspective
the molecule with more letters in its alphabet must be the one doing
the storing. The fact that the molecule with four letters won is the
first piece of bioinformatics intuition you should internalise: an
alphabet's size matters less than how many symbols you string together.


== The Double Helix <sec:dna>

DNA is a polymer of four nucleotides, each of which has three pieces:
a five-carbon sugar (deoxyribose), a phosphate group, and one of four
nitrogen-containing bases — adenine (A), thymine (T), guanine (G), or
cytosine (C). The sugar and phosphate of consecutive nucleotides form
the backbone of the polymer through phosphodiester bonds. The bases
project inward.

Two strands of this polymer wrap around each other into a right-handed
double helix, joined at the bases by hydrogen bonds. The pairing is
specific: A pairs with T through two hydrogen bonds, G pairs with C
through three. This *Watson–Crick base pairing* gives the molecule
its most important property — each strand is a template for the other.
Given one strand, the sequence of the other follows mechanically.

The two strands run in opposite directions. Each strand has a
"5' end" and a "3' end" — these labels refer to which carbon of the
sugar carries the next phosphate. By convention we always write
sequences from 5' to 3', the direction in which DNA polymerases read
and synthesise. When we say a read is "AGCTAA", we mean the 5' end
is the A on the left.

#figure(
  image("../../diagrams/lecture-01/02-dna-double-helix.svg", width: 80%),
  caption: [
    The double helix in cartoon form. Two antiparallel strands; bases
    pair internally; the sugar-phosphate backbone runs along the
    outside in both directions.
  ],
) <fig:helix>

The geometry of the helix is surprisingly regular. B-form DNA (the
common physiological form) has a diameter of about 20 Å, makes one
complete turn every 34 Å, and packs roughly 10.5 base pairs per turn.
The bases stack on top of each other inside the helix at 3.4 Å vertical
spacing. The stacking generates a strong UV absorbance at 260 nm that
gives bioinformaticians their favourite spectrophotometric assay —
"OD-260" for nucleic-acid concentration.

Two grooves run along the outside of the helix. The *major groove* is
about 22 Å wide; the *minor groove*, about 12 Å. The base pairs are
exposed differently in each groove — in particular, the major groove
displays edges of the four base pairs that a protein can distinguish
without unwinding the DNA. Almost every sequence-specific protein you
will hear about — transcription factors, restriction enzymes,
CRISPR-Cas9 — reads DNA through the major groove.

#figure(
  image("../figures/ch01/f1-dna-grooves.svg", width: 75%),
  caption: [
    Cross-section of B-form DNA showing the major and minor grooves.
    Sequence-specific proteins read the bases through the major groove
    without unwinding the helix.
  ],
) <fig:grooves>

#tip[
  The asymmetry of the grooves explains why most protein-DNA contact
  patterns are written as "5'-NNNNNN-3'" rather than as base
  identities — the protein is reading hydrogen-bond donors and
  acceptors at the floor of the major groove, not "an A" or "a T".
]

=== Topology

A long thin polymer like DNA cannot just float around inside a cell:
2 metres of human DNA has to fit inside a nucleus 5 microns across.
Cells solve this by *supercoiling* the molecule. Imagine taking a
telephone cord (long, flexible) and twisting it in your hand — it
coils on itself, knots locally, stores energy in the twist. Real DNA
does exactly this. Topoisomerases, a family of enzymes, manage the
twist by introducing or removing supercoils.

The topology matters for sequencing only obliquely — every protocol
starts by chemically and physically denaturing the DNA — but it
matters constantly in the cell. Replication, transcription, and DNA
repair all change the molecule's topology, and every cell maintains
an active topology budget.


== The Central Dogma <sec:dogma>

The information stored in DNA is read out in two steps. *Transcription*
copies a stretch of DNA into a single-stranded RNA molecule called
messenger RNA (mRNA). *Translation* reads the mRNA in three-nucleotide
groups called codons, each of which specifies one of twenty amino acids
(or a stop signal), and synthesises a protein.

This is the *central dogma of molecular biology*, named by Francis Crick
in 1957. It is a directionality claim: information flows from nucleic
acid to protein, not the other way around. The dogma has exceptions —
RNA can be reverse-transcribed back to DNA (HIV, retrotransposons),
and some viruses use RNA directly as their genome — but the prohibition
that proteins cannot dictate nucleic-acid sequence has held up.

#figure(
  image("../../diagrams/lecture-01/03-central-dogma-flow.svg", width: 90%),
  caption: [
    DNA is transcribed into mRNA; mRNA is translated into protein.
    Reverse transcription (RNA → DNA) is the most important documented
    exception to the unidirectional flow.
  ],
) <fig:dogma>

=== The Genetic Code

The mapping from codons to amino acids is the genetic code. There are
$4^3 = 64$ possible codons but only 20 amino acids plus a stop signal,
so the code is *degenerate*: most amino acids have multiple codons.
The degeneracy is not random — synonymous codons usually differ only
in the third position, a fact that lets the translation machinery
tolerate single-base errors.

#figure(
  image("../../diagrams/lecture-01/04-genetic-code-table.svg", width: 100%),
  caption: [
    The standard genetic code. Each three-nucleotide codon maps to
    one amino acid, with three stop codons (UAA, UAG, UGA). The third
    position is "wobble" — most synonymous codons differ only there.
  ],
) <fig:code>

The code is nearly universal — the same in bacteria, plants, and
humans, with a handful of small variants in mitochondria and some
protozoa. Universality is itself evidence of common ancestry: the code
was fixed before the last universal common ancestor diversified, more
than three billion years ago.

#note[
  When we talk about a *gene* in this book we usually mean a stretch of
  genomic DNA that, after transcription and (in eukaryotes) splicing,
  produces a single mature mRNA encoding a protein. Eukaryotic genes
  are interrupted by *introns* — sequences that are transcribed but
  spliced out before translation. Bacterial genes are not. The
  details of splicing show up again in Chapters 5 and 6.
]


== Reading the Genome: A Brief Historical Arc <sec:history>

Modern DNA sequencing began in the late 1970s. Frederick Sanger and
Allan Maxam, working independently, each published a way to determine
the sequence of bases in a piece of DNA. Sanger's method — using
chain-terminating dideoxynucleotides — proved more scalable and
dominated for the next thirty years.

#figure(
  image("../../diagrams/lecture-01/06-sanger-termination.svg", width: 85%),
  caption: [
    Sanger sequencing by dideoxy chain termination. Four reactions,
    each spiked with a different fluorescently labelled ddNTP. The
    fragments are separated by capillary electrophoresis; the
    sequence is read off the electropherogram.
  ],
) <fig:sanger>

The first revolution came in the late 1980s with automated capillary
electrophoresis, which replaced gel slabs with thin fused-silica
capillaries and let the throughput grow from a few hundred bases per day
per machine to a few hundred thousand. The Human Genome Project, which
ran from 1990 to 2003 and ultimately cost about three billion dollars,
was a capillary-sequencing project. Sanger-method bases dominate the
reference genome you still use today.

The second revolution — *next-generation sequencing* — arrived in 2005
with the 454 platform from Roche, which read DNA by detecting
pyrophosphate release after each base addition. 454 was followed by
Solexa (later Illumina) in 2007, SOLiD from Applied Biosystems shortly
after, and a long tail of competing platforms. The defining feature of
NGS was *parallelism*: instead of running a few dozen sequencing
reactions side by side in a 96-well plate, an NGS flow cell held
hundreds of millions of separate clonal clusters and read them all at
once.

The cost per base fell by six orders of magnitude in fifteen years —
faster than Moore's Law over the same period. The shape of the curve
is one of the most famous graphs in modern biology.

#figure(
  image("../figures/ch01/f2-cost-timeline.svg", width: 95%),
  caption: [
    The NHGRI "Cost per Genome" curve, plotted log-scale from 2001 to
    2025. The discontinuity in 2008 marks the displacement of capillary
    Sanger by short-read NGS; subsequent improvements track better
    chemistry, longer reads, and higher cluster density.
  ],
) <fig:cost>

#matters[
  The cost curve is the single most important context for everything
  in this book. Algorithms that were impossible in 2003 because the
  data was scarce, and algorithms that were impossible in 2013 because
  the data was too plentiful, are both routine today. Where the
  bottleneck moves next — and it _will_ move — is the question to keep
  asking as you read.
]

The third era is now well underway. Long-read platforms — PacBio's
single-molecule real-time (SMRT) sequencing from 2011 and Oxford
Nanopore's protein-pore platform from 2014 — read individual molecules
ten to a hundred kilobases long, at the cost of higher per-base error
rates. The first complete, telomere-to-telomere human genome assembly
was published in 2022, twenty years after the "complete" draft, using
long-read data to fill in the centromeres, ribosomal arrays, and other
repetitive regions that short reads could never resolve. Chapter 11
covers long reads and the pangenome.


== Illumina: Sequencing by Synthesis <sec:illumina>

Illumina sequencing dominates the field by volume and will be the
default platform for the rest of this book unless we say otherwise.
It is worth understanding the physics of the instrument because most
of the systematic errors you will encounter in FASTQ files come
directly from limitations of this physics.

The protocol has three stages: *library preparation* (turning a tube
of genomic DNA into a tube of short, adapter-flanked fragments),
*cluster generation* (anchoring those fragments to a glass flow cell
and clonally amplifying each one in place), and *sequencing*
(extending each cluster's DNA one base at a time, photographing the
fluorescent signal, and decoding the base from the color).

=== Library Preparation

The input to a sequencer is never a long strand of genomic DNA — it
is a population of short fragments, each capped at both ends with
synthetic *adapter* sequences. The adapters serve two purposes: they
let the fragment hybridise to the flow cell, and they carry primer
binding sites for the sequencing reaction itself.

#figure(
  image("../../diagrams/lecture-01/09-adapter-structure.svg", width: 95%),
  caption: [
    Illumina adapter structure. P5 and P7 anneal to the flow cell;
    Rd1 and Rd2 are sequencing primer binding sites; the index
    bar-coded segment lets multiple samples be pooled on one run.
  ],
) <fig:adapter>

Fragmentation is mechanical (sonication) or enzymatic (transposase).
Adapter ligation happens in solution. The library is then size-selected
on a magnetic bead system and PCR-amplified just enough to give
detectable signal — too many PCR cycles introduces duplicate reads
and amplification artefacts, both of which complicate downstream
analysis. The exact protocol varies by application but the shape is
universal: fragment, ligate, size-select, amplify lightly.

=== Cluster Generation by Bridge Amplification

The fragmented library is loaded onto a flow cell — a glass slide
patterned with a high-density lawn of P5 and P7 oligonucleotide
"capture probes." Each library fragment, by virtue of its adapters,
hybridises to one capture probe. A polymerase extends across the
fragment, generating a copy attached to the slide. The fragment then
"bridges" — bends over and hybridises to a nearby capture probe of the
opposite type, where a new round of extension makes a second copy.
Repeating bridge amplification for thirty-five cycles produces a
clonal cluster of about a thousand identical molecules in roughly a
one-micron spot.

#figure(
  image("../../diagrams/lecture-01/08-bridge-amplification.svg", width: 90%),
  caption: [
    Bridge amplification. A single library molecule hybridises to one
    end of a flow-cell capture probe, bends over to the other,
    extends, and cycles. After roughly 35 rounds, each starting
    molecule has produced a clonal cluster of about 1000 copies.
  ],
) <fig:bridge>

Why bother with the cluster? Because the optical readout that follows
needs photons. A single molecule emits very few photons per imaging
cycle, and dyes bleach quickly. Amplifying the molecule into a
thousand-copy cluster gives a thousand times more signal per imaging
event, which is what makes the basecall reliable.

=== Sequencing by Synthesis

With the flow cell loaded, the actual sequencing is conceptually
simple. The instrument floods the flow cell with sequencing primer,
four fluorescently labelled "reversible terminator" nucleotides, and
DNA polymerase. The polymerase adds exactly one terminator nucleotide
to each cluster's growing strand and stops. The instrument images the
flow cell in four channels (one per base), records the dominant color
at each cluster, then chemically removes the terminator and the
fluorescent label. The cluster is now ready to accept the next base.
Repeat for 75, 100, 150, or 300 cycles depending on the instrument and
read mode.

#figure(
  image("../figures/ch01/f3-intensity-trace.svg", width: 95%),
  caption: [
    Four-color intensity trace for a single Illumina cluster across
    20 cycles. At each cycle, one channel dominates; the base call is
    the channel with the highest intensity. The asymmetry between
    channel intensities is real and shapes per-base quality scores.
  ],
) <fig:intensity>

The most important thing to internalise about sequencing-by-synthesis
is that the basecall is a *signal-processing decision*. The raw data
at each cycle is a vector of four fluorescence intensities. Calling a
base means choosing the largest of the four. The quality of that
decision — captured by the Phred quality score we will define in a
moment — depends on how confidently the largest intensity dominates
the others.

#note[
  The two dominant sources of error in Illumina sequencing are
  *phasing* and *prephasing*. Phasing happens when a polymerase fails
  to extend during one cycle and stays behind; prephasing happens when
  it adds two bases in one cycle. After many cycles, an unsynchronised
  cluster looks like a smeared mixture of two adjacent positions in
  the read. This is why Illumina error rates climb sharply at the end
  of long reads.
]


== Long-Read Platforms <sec:long-reads>

Short Illumina reads — 75 to 300 bases, typically — are excellent for
counting and for variant detection. They are not enough to span the
many regions of a real genome that contain repetitive elements longer
than the read itself. Centromeres, ribosomal-RNA arrays, segmental
duplications, and many disease-associated tandem repeats are simply
invisible to short-read sequencing alone.

*PacBio* and *Oxford Nanopore* solve this differently. PacBio's
single-molecule real-time (SMRT) platform watches a single DNA
polymerase work on a single template inside a tiny illuminated well —
a "zero-mode waveguide" — and detects fluorescent nucleotides as they
are incorporated. Read lengths reach tens of kilobases; circular
consensus reads ("HiFi") deliver Q30+ accuracy at ~20 kb reads by
sequencing the same molecule multiple times.

#figure(
  image("../../diagrams/lecture-01/10-pacbio-zmw.svg", width: 80%),
  caption: [
    A PacBio zero-mode waveguide. The aluminium well is 70 nm wide —
    smaller than the wavelength of visible light — so excitation only
    reaches the polymerase at the bottom. Each fluorescent
    nucleotide is detected only when it is incorporated.
  ],
) <fig:zmw>

*Oxford Nanopore* works without a polymerase at all. A motor protein
threads a single DNA strand through a biological pore embedded in a
membrane held at a fixed voltage. As the strand moves through, it
modulates the ionic current across the pore, and the current pattern
is decoded into a base sequence by a recurrent neural network. Read
lengths can exceed a megabase. The raw signal is a continuous current
trace — a "squiggle" — that gets translated to bases by software, not
by chemistry.

#figure(
  image("../../diagrams/lecture-01/11-nanopore-squiggle.svg", width: 95%),
  caption: [
    Nanopore "squiggle": picoampere-scale ionic current versus time as
    a single DNA strand translocates through a protein pore. The
    basecaller is a neural network that maps current segments back
    to bases.
  ],
) <fig:nanopore>

Long-read platforms shift the engineering tradeoff: read length goes
up by a factor of a hundred, but per-base error rates were
historically higher (5–15% for ONT in 2018; 0.1–1% for HiFi). Modern
ONT chemistry has narrowed the gap dramatically. Chapter 11 returns
to long reads when we cover assembly.


== The FASTQ Format <sec:fastq>

Whatever platform produced the reads, the output that lands in your
working directory is a FASTQ file. FASTQ is a four-line-per-record
text format that pairs each read with a per-base quality score.

```
@SRR000001.1 read_id_metadata
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTC
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>
```

The four lines, in order, are:

1. A header line that starts with `@` and contains the read identifier
   (and, in some formats, the run, lane, and tile information).
2. The base sequence itself, in ACGT (and occasionally N for ambiguous).
3. A separator line that is always `+` (occasionally followed by a
   repeat of the header for redundancy).
4. The per-base quality string, the same length as the sequence.

#figure(
  image("../../diagrams/lecture-01/12-fastq-anatomy.svg", width: 95%),
  caption: [
    Anatomy of a FASTQ record. Identifier, sequence, separator,
    quality string. The quality string is the same length as the
    sequence and encodes each base's Phred score as an ASCII character.
  ],
) <fig:fastq>

=== Phred Quality Scores

The quality string encodes a *Phred quality score* per base. A Phred
score $Q$ is a log-transformed probability that the base is wrong:

$ Q = -10 log_10 P_("err") $

A score of Q10 means a 1 in 10 chance of error; Q20 is 1 in 100; Q30
is 1 in 1000; Q40 is 1 in 10,000. Q30 is the de-facto "good base"
threshold for Illumina; anything above Q40 from a single read is
suspicious for that platform.

To fit the score into one ASCII character per base, the score is
offset and ASCII-encoded. Modern FASTQ uses *Phred+33* encoding: the
character is `chr(Q + 33)`. So Q0 is `!`, Q10 is `+`, Q30 is `?`, Q40
is `I`.

#figure(
  image("../figures/ch01/f4-phred-encoding.svg", width: 95%),
  caption: [
    Phred quality conversion. The Q-score → P-error transform is
    log-scaled; the ASCII character is `chr(Q + 33)`. Below: a
    25-base read with each base coloured by its Q-score (red = low,
    green = high). Low-quality tails are the typical Illumina
    signature.
  ],
) <fig:phred>

#warn[
  Historical Phred+64 encoding (Illumina 1.3–1.7, deprecated since
  2011) is still occasionally found in archived data. The two
  encodings differ by 31 in the ASCII value of every quality
  character. Misidentifying the encoding will silently shift every
  base by 31 Q-points, with predictably catastrophic downstream
  consequences. Tools like `fastqc` autodetect; for archival data
  always confirm.
]

=== Paired Reads, Index Reads, and Other Lies

In practice you rarely have just _one_ FASTQ file. Paired-end Illumina
sequencing reads each cluster from both ends, producing two FASTQ
files where each record in file 1 has a matching record at the same
position in file 2. The two files are named with `_R1` and `_R2`
suffixes. Read pairing is enforced by record-position in the file —
never reorder a paired-end FASTQ without a tool that knows.

Multiplexed runs further split each lane by *sample index* — a
short barcode written into the adapter. The index is read as a
separate cycle on the instrument and written into the read header,
or sometimes into separate index FASTQ files. Demultiplexing happens
either on the instrument or as a first step in the analysis pipeline;
either way, by the time a FASTQ reaches you it has usually been split
per-sample already.

#tip[
  The single most useful early-pipeline reality check is `seqkit stats`
  (or `samtools fqimport` + `samtools stats`) on every FASTQ before
  you run anything else. Read count, read length distribution, mean
  Q, and percent N catch a remarkable fraction of upstream problems
  — wrong file pairing, truncated downloads, encoding misidentification —
  before they propagate.
]


== The Pipeline Ahead <sec:pipeline>

#figure(
  image("../../diagrams/lecture-01/05-bioinformatics-pipeline.svg", width: 100%),
  caption: [
    The canonical short-read bioinformatics pipeline. Each box is a
    later chapter: alignment (Chapter 2), assembly (Chapter 3),
    variant calling (Chapter 4), expression analysis (Chapters 5–8),
    and so on through the rest of the book.
  ],
) <fig:pipeline>

Everything from here onward starts with a FASTQ file. Chapter 2 covers
the most common first step: aligning short reads back to a reference
genome. Chapter 3 covers the alternative when you have no reference —
_de novo_ assembly. Chapter 4 turns aligned reads into a list of
variants. Chapters 5 through 8 work with RNA instead of DNA, asking
which genes are expressed and at what level. The remaining chapters
go deeper into each subspeciality.

In every case the basic flow is the same: a sequencing instrument
emits FASTQ, an analysis tool turns FASTQ into something more
interpretable, and a statistical model turns the interpretable thing
into a biological claim. The art of bioinformatics is knowing the
failure modes of each link in that chain. The next 26 chapters are
a tour of those failure modes, structured by the data type and the
biological question.


== Summary <sec:summary>

- A genome is a polymer of four nucleotides organised into a double
  helix. The two strands are antiparallel; each is the complement of
  the other; the sequence is conventionally written 5' to 3'.
- Information flows from DNA to RNA to protein. The genetic code is
  degenerate, near-universal, and almost three billion years old.
- A modern sequencer converts a tube of DNA into a text file. Illumina
  reads short (75–300 bp), at high accuracy, in massive parallel.
  PacBio and Nanopore read long (10kb–1Mb), at lower per-pass accuracy.
- The output text file is FASTQ: four lines per read, with a per-base
  Phred quality score. A Q-score of 30 means a 1-in-1000 chance the
  base is wrong. The encoding offset is `+33` for modern data.
- Every subsequent chapter of this book begins with a FASTQ.


== Exercises <sec:exercises>

#strong[1. ASCII to Phred.] Given the quality string
`!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>`, compute the mean
Q-score, the minimum Q-score, and the position of the minimum. (The
encoding is Phred+33. Hint: in Python, `ord(c) - 33` for each
character.)

#strong[2. Reverse complement.] Write a function that takes a DNA
sequence and returns its reverse complement. Test on
`AGCTTGCAATGCATGAATAG`. Make sure your function handles lowercase
input and that complementing `N` returns `N`.

#strong[3. Codon table check.] Using the genetic code in
@fig:code, translate the sequence `ATGGCCTGAAGCAGATGA` into a protein.
Where does translation start? Where does it stop? What is the
resulting peptide?

#strong[4. Reading the cost curve.] @fig:cost shows the NHGRI
cost-per-genome curve on a log scale. Estimate the doubling time of
sequencing _output_ (bases per dollar) over (a) 2003–2008 and
(b) 2008–2015. Compare with Moore's Law (~18-month doubling). What
biological or engineering factors explain the difference between
the two intervals?

#strong[5. Phred discrimination.] Two reads have mean Q-scores of 25
and 35 respectively. What is the expected number of erroneous bases
per 100 sequenced bases for each? Express your answer to two
significant figures.

#strong[6. (Open-ended.)] Pick one of the platform-specific artefacts
mentioned in this chapter — Illumina phasing, ONT homopolymer slippage,
PCR duplicates — and find one published bioinformatics tool that
explicitly accounts for it. Briefly describe (one paragraph) the
correction it applies.


== Further Reading <sec:further-reading>

- *Watson, J. D., and Crick, F. H.* (1953). "Molecular Structure of
  Nucleic Acids." _Nature_ 171: 737–738. Two pages. Still readable.
- *Sanger, F., Nicklen, S., and Coulson, A. R.* (1977). "DNA
  Sequencing with Chain-Terminating Inhibitors." _PNAS_ 74: 5463–5467.
  The paper that started everything.
- *Bentley, D. R., et al.* (2008). "Accurate Whole Human Genome
  Sequencing using Reversible Terminator Chemistry." _Nature_ 456:
  53–59. The Illumina platform paper.
- *Eisenstein, M.* (2022). "The Complete Human Genome." _Nature
  Methods_ 19: 521–524. A short news piece covering the T2T consortium.
- *Wendisch, A., et al.* (2020). "Lessons Learned from
  High-Throughput Sequencing Technologies." _Annual Review of
  Genomics and Human Genetics_ 21: 109–130. A condensed review of the
  platform space.
