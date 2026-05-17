#import "../theme/book-theme.typ": *

= ChIP-seq, ATAC-seq, and Peak Calling: Detection Theory on the Genome <ch:chip-atac>

#matters[
  DNA sequencing tells you what the letters are; RNA sequencing tells
  you which letters get transcribed; neither tells you _why_. The cell
  decides which genes to turn on by physically marking and unwrapping
  the regions of chromatin where transcription factors land. Reading
  that decision requires a different kind of assay — one that produces
  a coverage signal stacked over regulatory elements rather than over
  genes or exons. This chapter walks the analysis pipeline for two such
  assays, ChIP-seq and ATAC-seq, from the wet-lab chemistry through to
  modern sequence-to-regulation deep-learning models. The mathematical
  spine is detection theory. Once you see the peak-calling problem as
  constant-false-alarm-rate detection against a local Poisson null, the
  rest of the chapter falls into place.
]

A short coverage track does not look like much. A few thousand reads
piled into a 200-base-pair window above background, flanked by a
near-empty desert; another smaller pile-up a few kilobases over; the
genome scrolling on either side. The data type is identical to anything
this book has covered so far — aligned BAM, integer counts per genomic
interval — and the algorithms turn out to be ports of techniques that
EE students have already seen in radar, communications, and classical
DSP. The reason the chapter is worth your time is that the _biological
question_ is new. Regulation lives in a layer between sequence and
transcript that neither Chapter 1's DNA assays nor Chapter 5's RNA
assays can see. ChIP-seq and ATAC-seq are the two workhorse experiments
that open up that layer, and the bioinformatics pipeline they share is
a self-contained tour of detection theory on the genome.

The chapter has six sections. @sec:assays introduces the two assays,
the biological questions each answers, and the kind of signal each
deposits in a BAM. @sec:libraries works through library preparation —
crosslinking, sonication, Tn5 transposition — and the QC signatures
each protocol leaves in the fragment-length distribution. @sec:peaks
is the long middle of the chapter: the MACS2 peak caller as a
constant-false-alarm-rate detector against a local Poisson null, with
narrow versus broad modes and the practical choices around each. @sec:diffbind
ports the negative-binomial GLM machinery from Chapter 6's
differential-expression analysis to peaks-versus-samples count matrices.
@sec:motifs introduces position weight matrices and the matched-filter
view of motif scanning, then covers ATAC footprinting as a $sqrt(N)$
signal-averaging problem. @sec:enformer closes with the sequence-to-regulation
deep-learning models that have started to predict the assay output
directly from DNA.

Every step in this chapter assumes you have a coordinate-sorted,
duplicate-marked BAM aligned to a current human or mouse reference.
The alignment step is identical to Chapter 2; the variant-calling
machinery is not needed; the read-counting and statistical machinery
will be familiar from Chapter 6.


== What ChIP-seq and ATAC-seq Measure <sec:assays>

A typical human cell carries roughly two metres of DNA folded into a
five-micron nucleus. That DNA is not bare — it is wrapped around
histone octamers at 147 base-pair intervals to form _nucleosomes_, and
those nucleosomes are packaged into higher-order structures whose
density varies dramatically along the chromosome. Some regions are
tightly compacted into heterochromatin where no transcription factor
can physically reach the DNA; others sit in loose, accessible
configurations where the regulatory machinery can land. The cell uses
this physical accessibility, plus a layer of covalent histone
modifications and DNA methylation marks, to control which genes are
transcribed in which cell types. Two strands of identical DNA in two
sister cells can be expressed in radically different ways because the
chromatin packaging differs.

*ChIP-seq* (Chromatin ImmunoPrecipitation sequencing) measures where a
specific protein sits on the DNA. You start with a population of cells,
treat them with formaldehyde so that any protein currently touching the
DNA is covalently locked in place, fragment the chromatin, and use an
antibody to pull down fragments carrying the protein of interest.
Reverse the crosslinks, purify the surviving DNA, and sequence it. The
reads pile up at the genomic positions where the protein was bound at
the moment of crosslinking. Targets are most often _transcription
factors_ — CTCF, p53, the oestrogen receptor — or _histone marks_, the
covalent modifications on histone tails that index the state of the
underlying chromatin (H3K4me3 at active promoters, H3K27ac at active
enhancers, H3K27me3 at silenced regions, H3K9me3 at heterochromatin).

*ATAC-seq* (Assay for Transposase-Accessible Chromatin sequencing)
measures where chromatin is open, without committing to a specific
target. The trick is a single enzyme, _Tn5 transposase_, pre-loaded
with sequencing adapters as its payload. Apply Tn5 to permeabilised
cells and it cuts the DNA wherever it can physically reach — which is
only the open, nucleosome-free regions — and inserts adapters at the
cut sites in a single chemical step. The library that emerges is
already adapter-flanked and ready to sequence; the only DNA that gets
amplified is DNA that was accessible to Tn5 when the experiment ran.
ATAC was introduced by Buenrostro and colleagues in 2013 and within
three years had largely displaced its older cousin DNase-seq (which
used DNaseI digestion to do the same job, with worse signal-to-noise
and ten times the input cells).

The two assays are complementary, not competing. ChIP tells you _which
protein_ is bound at a specific position; ATAC tells you _whether the
position is open at all_. A typical project runs ATAC first to find the
accessible regulatory regions, then runs targeted ChIP for the
transcription factors of interest within those regions. Both produce
coverage tracks dominated by sharp pile-ups at regulatory sites,
flanked by near-empty background.

#figure(
  image("../../diagrams/lecture-09/01-chip-vs-atac-overview.svg", width: 95%),
  caption: [
    The two assays in side-by-side schematic form. ChIP-seq pulls down
    DNA bound by a specific antibody target; ATAC-seq maps all
    Tn5-accessible regions in one reaction. Below: the shared cartoon
    chromatin with nucleosomes and TF binding sites, aligned to the
    coverage signal each assay produces.
  ],
) <fig:chip-vs-atac>

#note[
  A gene is not "on" because of its sequence alone. The sequence is the
  contract; the regulatory state of the chromatin is the enforcement.
  ChIP-seq asks which proteins are reading the contract at this
  moment; ATAC-seq asks which parts of the contract are open for
  reading. The combination tells you which regulatory elements are
  active in the cells you sampled — a view neither whole-genome
  sequencing nor RNA-seq can give you.
]

The modern ChIP-seq era began in 2007 with two near-simultaneous
papers, Barski and colleagues in _Cell_ and Robertson and colleagues in
_Nature Methods_, that combined the older microarray-readout assay
ChIP-chip with the freshly available Solexa/Illumina short-read
sequencer. The follow-on project, ENCODE — originally a 2003 pilot
using microarrays, relaunched in 2007 as a sequencing-era effort —
profiled hundreds of transcription factors and histone marks across
dozens of human cell lines and released a single coordinated dataset in
September 2012. The ENCODE release remains the largest single-event
data drop in functional genomics, and the default reference set against
which everything downstream is aligned. If you scan a peak list from a
random ATAC-seq paper, the odds that its consensus accessibility map
borrows from ENCODE are close to one.

=== The Single-Cell Variant

The same Tn5 chemistry runs at per-cell resolution in droplet scATAC-seq
and in the 10x Multiome platform (paired RNA + ATAC on the same
nucleus). The output is dramatically sparser — roughly ten thousand
fragments per cell versus thirty million reads per bulk sample — and
the dimensionality-reduction step uses Latent Semantic Indexing (an
inverse-document-frequency-weighted SVD borrowed from text retrieval)
instead of the PCA familiar from bulk RNA-seq. But the peak-calling
machinery in this chapter transfers cleanly: aggregate per-cell
fragments into a pseudo-bulk track per cluster, then call peaks as
though the pseudo-bulks were bulk. Differential accessibility, motif
scanning, and footprinting all run on the pseudo-bulks. Chapter 8
covers the single-cell side in detail; for the rest of this chapter we
will treat ATAC and ChIP as bulk assays.


== Library Preparation, Biases, and Alignment <sec:libraries>

The wet-lab mechanics matter because each step shapes what the
analysis sees. Skipping the chemistry produces a pipeline that works on
clean reference datasets and breaks silently on real biology.

=== Crosslinking, Sonication, and Tn5 Transposition

*Formaldehyde crosslinking* is the time-stopping move that makes
ChIP-seq possible. Formaldehyde reacts with primary amines on proteins
and with the N7 of adenine and guanine in DNA, forming short covalent
bonds with a reach of about 2 Å. One per cent formaldehyde for ten
minutes at room temperature is standard. Too little crosslink and
proteins wander off their binding sites during the subsequent hours of
fragmentation and pull-down; you end up sequencing where the protein
_was at some point in the protocol_, not where it was at the moment you
added formaldehyde. Too much crosslink and you weld everything to
everything — the antibody pulls down a fog of indirectly attached DNA
and binding specificity drops out. Published protocols report the
concentration and time explicitly; deviations are usually deliberate.

*Sonication* is the physical follow-up that breaks crosslinked
chromatin into ~300 bp fragments using ultrasound. The fragment-length
distribution is roughly log-normal with a long tail. Sonication has its
own subtle biases — it cuts preferentially at nucleosome-free regions,
slightly depleting fragments from heavily packaged heterochromatin — but
for peak calling these biases average out over a typical 200-bp test
window and rarely cause downstream problems.

*Tn5 transposition* replaces all of the above for ATAC-seq. Tn5 is a
prokaryotic transposase whose native role is to cut DNA at its
recognition sequence and insert a transposon in a single reaction. The
ATAC trick is to pre-load Tn5 with sequencing adapters in place of the
natural transposon payload, then apply the loaded enzyme to lightly
permeabilised cells. Tn5 cuts wherever it can physically reach the
DNA, leaving the adapters covalently attached at the cut sites. The
result is a sequencing library that is both adapter-flanked and
accessibility-selected in one chemical step, from a starting input of
~50,000 cells, in under four hours of bench time. The protocol's
combination of simplicity and signal-to-noise is what made ATAC-seq
displace DNase-seq as the default accessibility assay.

The catch is that Tn5 has _sequence preference_. The enzyme cuts more
often at certain ~9-bp local motifs than at others, with a typical
preferred motif resembling `NNNYMNNHN`. Even in a perfectly open
region of DNA, the raw per-base Tn5 cut counts are non-uniform — and
that non-uniformity is precisely what footprint analysis (Section 9.5)
has to deconvolve before it can identify protein-protection signals.
For peak calling the bias is small enough to average out over a 200-bp
window; for footprinting it is not.

#figure(
  image("../../diagrams/lecture-09/03-atac-tn5.svg", width: 90%),
  caption: [
    Tn5 inserts adapters only in nucleosome-free regions; the fragment
    length records how many nucleosomes sit between the two adapter
    insertion points. Sub-nucleosomal fragments span no nucleosome;
    mono-nucleosomal fragments span exactly one.
  ],
) <fig:tn5>

=== The Fragment-Length Signature

Both assays produce libraries with characteristic fragment-length
distributions, and looking at the distribution is the first quality
check you should run on any ATAC or ChIP sample.

ChIP-seq fragments are roughly log-normal around the sonication target
of ~300 bp, with no internal structure because sonication shears
randomly. ATAC-seq fragments are far more interesting: the distribution
is _multi-modal_, with peaks at sub-nucleosomal (< 100 bp,
nucleosome-free linker DNA), mono-nucleosomal (150–180 bp, one
nucleosome between the two Tn5 cuts), di-nucleosomal (~300 bp), and a
faint tri-nucleosomal mode (~450 bp). The 147 bp of DNA wrapped around
a nucleosome's histone octamer plus a ~20 bp linker on each side
explains the spacing exactly.

#figure(
  image("../figures/ch09/f1-fragment-ladder.svg", width: 95%),
  caption: [
    ATAC-seq fragment-length distribution on a log-density axis,
    showing the characteristic sub / mono / di / tri-nucleosomal
    ladder. Overlaid in dashed amber: a typical sonicated ChIP-seq
    library with no laddering. Absence of the ATAC ladder is the
    first hard sign of a failed experiment.
  ],
) <fig:fragment-ladder>

The ATAC ladder is the first QC gate: a sample whose fragment-length
distribution shows a clean nucleosomal pattern is almost certainly
salvageable, and a sample that does not show the pattern is almost
certainly not. Common failure modes include too little Tn5 (everything
ends up at the long tail), too much Tn5 (the library over-fragments
into a single sub-nucleosomal mode), degraded chromatin (the ladder
flattens), and over-aggressive size selection (the longer modes get
clipped). No downstream analytical fix recovers a sample that fails
this QC.

#warn[
  Tn5 sequence preference matters in two stages: it _doesn't_ matter
  much for peak calling, because the bias averages out over the 200-bp
  windows MACS2 tests; it _does_ matter for footprinting, where the
  analysis works at single-base resolution and the cut-site bias
  appears as a systematic confound. Always run TOBIAS, HINT-ATAC, or a
  similar Tn5-bias correction before interpreting a footprint as a
  protein-protection signal. Raw cut counts are not an unbiased
  accessibility measurement.
]

=== Alignment Quirks

Alignment is mostly a callback to Chapter 2: BWA-MEM or Bowtie2 on
paired-end reads, expecting most reads to map uniquely, with the same
mark-duplicates pass that variant-calling pipelines run. Two
ChIP-and-ATAC-specific wrinkles modify the picture.

*The Tn5 offset correction.* Tn5 inserts adapters as a dimer, leaving a
9 bp staggered cut on the two strands. The 5' end of a forward-strand
read maps 4 bp downstream of the cut site; the 5' end of a reverse-strand
read maps 5 bp upstream. To recover the exact Tn5 cut position you
shift forward-strand reads by +4 bp and reverse-strand reads by −5 bp
before building the cut-site pile-up. The ENCODE ATAC-seq pipeline does
this routinely. Skipping the shift blurs footprint analyses by about
9 bp and pushes the footprint trough off-centre relative to the
underlying motif.

*Blacklist regions.* Some genomic regions produce artefactual high
coverage no matter what assay you run: rRNA gene clusters, telomeric
and centromeric repeats, regions with annotation errors. These
positions systematically produce false-positive peaks. The ENCODE
Blacklist (Amemiya, Kundaje, & Boyle 2019) is a curated set of ~400
GRCh38 intervals that every accessibility analysis should exclude
before peak calling.

*PCR duplicates.* Illumina libraries undergo PCR amplification, and
duplicate fragments from the same pre-amplification template should be
collapsed before peak calling. Unlike single-cell RNA-seq, ATAC and
ChIP libraries usually lack UMIs; deduplication uses
mapping-position + fragment-length as a duplicate key. Picard
MarkDuplicates or `samtools markdup` does the job.

#tip[
  Run the fragment-length plot, the read-depth histogram per
  chromosome, and the cross-correlation analysis (NSC and RSC, the
  ENCODE quality metrics) before you touch a peak caller. A bad
  experiment is much cheaper to detect at the BAM stage than after a
  week of downstream analysis on a peakset that should not have
  existed. The ENCODE chromatin-accessibility working group publishes
  reference thresholds for each metric per assay.
]


== Peak Calling as Detection Theory <sec:peaks>

The peak-calling problem is the central computational task of this
chapter. You align thirty million paired-end reads to the human
reference, get a coordinate-sorted BAM in hand, and need to produce a
BED file: a list of genomic intervals at which the read pile-up is
significantly above the local background. Some regions of the genome
will have many tens of reads stacked into a 200 bp window; most regions
will have one read or zero. The task is to decide, at every position on
the three-gigabase reference, whether the local count is signal or
noise.

This is textbook detection theory. At every candidate location, choose
between two hypotheses:

- $H_0$ (null): no binding event; the local coverage is background
  noise drawn from a Poisson process with rate $lambda_(text("local"))$.
- $H_1$ (alternative): a binding event is present; the local coverage
  is background plus signal.

The signal-to-noise ratio in a typical ChIP-seq experiment is
generous — a real TF binding site routinely shows ~5–20× the local
background — but the genome is three billion positions long, so even a
modest false-alarm rate of one in $10^5$ produces tens of thousands of
spurious peaks across a whole-genome scan. Multiple testing is not
optional. The shape of the answer is a detector whose false-alarm rate
is controlled across the entire genome, and whose threshold adapts to
local conditions so that real peaks above local hot regions are not
masked and noise above local cold regions is not falsely amplified.

#figure(
  image("../../diagrams/lecture-09/05-read-pileup.svg", width: 95%),
  caption: [
    A typical 5 kb genome-browser view. A sharp central pile-up sits
    above the local background; a smaller secondary rise sits one
    kilobase downstream. The peak-calling task is to flag the central
    region as a peak while leaving the surrounding background
    untouched.
  ],
) <fig:pileup>

#note[
  The peak-calling problem is structurally identical to _target
  detection in radar_. The genome coordinate plays the role of the
  range-Doppler bin; the per-window count plays the role of the echo
  amplitude; the local-window background plays the role of the clutter
  estimate. The radar tradition has a name for the family of
  detectors that adapt their threshold to local noise: _constant
  false-alarm rate_ (CFAR) detection. The classical reference is the
  cell-averaging CFAR detector of Finn & Johnson (1968), refined into
  order-statistic CFAR by Rohling (1983). MACS2 is one of these
  detectors, with the genome's coverage track standing in for the
  radar's range-bin time series.
]

=== MACS2: The Algorithm

The defining peak caller of the ChIP-seq era is MACS, _Model-based
Analysis of ChIP-Seq_, published by Zhang and colleagues in _Genome
Biology_ in 2008. MACS2 is the 2012 refinement (Feng et al. 2012) that
became the field's default and stayed there. MACS3, released in 2021,
adds ATAC-specific modes and improves memory use but keeps the core
detection formulation untouched. The fact that the formulation has
survived three versions and fifteen years of refinement is unusual; it
is a signal that the level of abstraction is correct.

The algorithm proceeds in six steps.

1. *Build the signal track.* Extend each aligned read by the estimated
   fragment length to recover an approximate fragment position, then
   compute the per-base fragment coverage across the genome. For
   paired-end data the fragment length is observed directly; for
   single-end data MACS2 estimates it from cross-correlation of the
   forward and reverse strand profiles.
2. *Slide a test window.* For each 200-bp window, count fragments
   whose summits fall inside the window. Call this count $c$. The
   200 bp width is chosen to match the resolution of a typical TF
   ChIP-seq peak; broader marks use larger windows.
3. *Estimate a local background rate $lambda_(text("local"))$.* This
   is the central trick. MACS2 computes several local Poisson rate
   estimates — the genome-wide average $lambda_(text("bg"))$, a
   1 kb window $lambda_(1000)$, a 5 kb window $lambda_(5000)$,
   a 10 kb window $lambda_(10000)$, and (if an input or IgG control is
   available) a matched-control count $lambda_(text("input"))$ — and
   takes the maximum:
   $ lambda_(text("local")) = max(lambda_(text("bg")), lambda_(1000),
     lambda_(5000), lambda_(10000), lambda_(text("input"))). $
4. *Compute a Poisson p-value.* Under $H_0$ the observed count $c$ is
   distributed as $"Poisson"(lambda_(text("local")))$, and the p-value
   is $P(X >= c | lambda_(text("local")))$ from the standard Poisson
   tail.
5. *Multi-test correction.* Apply Benjamini–Hochberg across all
   candidate windows in the genome (the same FDR machinery from
   Chapter 6, §6.4) to control the genome-wide false-discovery rate
   at 5 % by default.
6. *Merge adjacent significant windows* into a single peak interval,
   and emit the position of maximum coverage inside the merged region
   as the peak summit.

The whole pipeline runs in a few CPU-minutes on a 30 M-read sample.

#figure(
  image("../figures/ch09/f2-macs2-worked.svg", width: 95%),
  caption: [
    A worked MACS2 call at one candidate window. Multiple local
    background estimates are computed in nested windows around the
    candidate; the maximum becomes $lambda_(text("local"))$; the
    Poisson tail at the observed count produces a p-value that
    survives Benjamini–Hochberg correction.
  ],
) <fig:macs2-worked>

The cleverness lives in step 3. Real genomes violate the
uniform-background assumption in two distinct ways: _biological_ — open
chromatin generally pulls more reads than closed chromatin, regardless
of the specific target — and _technical_ — mappability varies, GC
content affects amplification efficiency, blacklist regions have
pathological pile-ups even after blacklist filtering. The local
estimate $lambda_(1000), lambda_(5000), lambda_(10000)$ absorbs both
sources of non-uniformity into the null distribution. A candidate peak
has to rise above _local_ noise, not just global average noise.

To make the numbers concrete, suppose a candidate 200 bp window shows
$c = 24$ fragments. The genome-wide rate is
$lambda_(text("bg")) = 0.5$ fragments per 200 bp (about 30 million
fragments spread across 3 Gb), the 1 kb estimate is $lambda_(1000) = 2.1$,
the 5 kb estimate is $lambda_(5000) = 3.8$, the 10 kb estimate is
$lambda_(10000) = 2.9$, and the matched input control gives
$lambda_(text("input")) = 1.2$. The local estimate is the maximum,
$lambda_(text("local")) = 3.8$. The Poisson tail at $c = 24$,
$lambda = 3.8$ is roughly $2.4 times 10^(-14)$ — well past
Benjamini–Hochberg threshold at any reasonable cohort size. The
candidate is called a peak.

#note[
  MACS2's max-of-several-windows trick is a _greatest-of CFAR_
  detector in the radar nomenclature. Classical CFAR splits the
  clutter estimate into leading and lagging guard cells around the
  cell-under-test; the greatest-of variant uses the larger of the two
  guard-cell averages to be robust against clutter edges where one
  side is biologically meaningful (a real adjacent peak, the start of
  a heterochromatin block). MACS2's 1 / 5 / 10 kb windows play the
  guard-cell role and the `max` is the greatest-of statistic. Radar
  CFAR has decades of analysis literature on the tradeoffs between
  cell-averaging, order-statistic, and greatest-of variants; the
  greatest-of choice is robust to outliers on the low side at the
  cost of slightly elevated miss rate near real peaks adjacent to
  noise — exactly the tradeoff MACS2 inherits.
]

=== Narrow Versus Broad Peaks

The default MACS2 mode produces narrow peaks: sharp intervals 100–500 bp
wide, suited to most TFs and to the punctate active-promoter and
active-enhancer histone marks. But not every regulatory signal looks
like a sharp peak. Some histone marks cover domains tens to hundreds of
kilobases wide — H3K27me3 over Polycomb-silenced regions, H3K9me3 over
heterochromatin, H3K36me3 over the bodies of actively transcribed
genes. Running default-narrow MACS2 on H3K27me3 fragments these
domains into a long string of sub-peaks, each of which is a real
fragment of one underlying biological signal.

#figure(
  image("../../diagrams/lecture-09/07-narrow-vs-broad.svg", width: 95%),
  caption: [
    Three signal shapes at the same coordinate scale. A typical TF
    ChIP track shows narrow, isolated peaks (CTCF, top). An
    active-promoter mark like H3K4me3 shows slightly wider peaks
    (middle). A silenced-domain mark like H3K27me3 covers broad
    plateaus (bottom) that narrow-peak calling fragments into
    spurious sub-peaks.
  ],
) <fig:narrow-vs-broad>

MACS2's `--broad` flag changes the strategy. The caller uses two
thresholds — a strict primary one to find core enriched regions, then a
more permissive secondary threshold to extend and merge those cores
into larger domains. The output is domain-scale peaks, possibly tens of
kilobases wide. The correct rule of thumb is: narrow-peak mode for TFs,
H3K4me3, H3K27ac, and ATAC-seq; broad-peak mode for H3K27me3, H3K9me3,
and H3K36me3. Specialised broad callers exist — SICER, PePr — that
outperform MACS2 on the very broadest marks but at the cost of an
extra dependency.

A post-processing step on H3K27ac peaks deserves separate mention.
*Super-enhancer analysis* (Whyte et al. 2013) merges nearby H3K27ac
peaks that together span more than ~10 kb of the genome and classifies
the merged region as a "super-enhancer" — a cluster of constituent
enhancers thought to mark master regulators of cell identity. The tool
is ROSE (Rank Ordering of Super-Enhancers). The biological interpretation
of super-enhancers as functionally distinct from ordinary enhancers
remains debated; the operational definition is straightforward.

#warn[
  Default-narrow MACS2 on a broad mark is the most common silent
  failure in chromatin pipelines. The caller will happily emit
  fifty thousand narrow "peaks" on an H3K27me3 dataset, the
  downstream tools will happily count reads inside them, the
  differential-binding GLM will happily produce a list of significant
  hits, and only weeks later will you notice that the underlying
  biology is a few dozen broad domains that have been chopped into
  thousands of fragments. Always pick the mode by the expected signal
  shape, not by the default.
]

=== When Input Controls Matter

A ChIP-seq experiment usually includes an _input control_ — DNA
sonicated from the same starting chromatin but processed without the
antibody pull-down step. The input pile-up captures all the non-specific
biases: mappability, GC content, amplification artefacts, blacklist
regions. Including the input as one of MACS2's $lambda$ estimates lets
the caller subtract those biases out of the background. For TF ChIP
the input typically does not change the called peaks much (the signal
is strong enough that any sensible background works), but for broad
histone marks and for low-enrichment TFs the input becomes essential.

ATAC-seq has no equivalent of input control; the closest analog is a
naked-DNA Tn5 reaction that exposes Tn5 sequence bias but not
chromatin-accessibility bias, and it is rarely run in practice. ATAC's
peak calling therefore relies entirely on the local-window
$lambda$ estimates, which is part of why the local-adaptive trick
matters so much.


== Differential Binding and Accessibility <sec:diffbind>

You have a peak list for one condition. Now your collaborator runs the
same experiment in a knockout, a drug-treated sample, or a different
tissue, and asks: _which peaks changed_? This is the chromatin
equivalent of differential expression, and the statistical machinery
ported almost verbatim from Chapter 6.

The analysis template:

1. Call peaks separately in each sample (MACS2 per sample), or on
   pooled reads across all samples.
2. Build a _consensus peakset_ — the union of peaks across all
   samples, or peaks seen in at least $N$ samples.
3. For each peak × each sample, count reads that overlap the peak.
   The output is a peaks-by-samples count matrix — the same shape as
   Chapter 6's genes-by-samples matrix, with peaks in the row position
   instead of genes.
4. Fit a negative-binomial GLM (DESeq2, edgeR), test each peak for
   differential count between conditions, apply Benjamini–Hochberg
   correction.

Step 4 is literally Chapter 6 again. The count-based DE toolkit works
here because the noise model is identical: integer counts with
overdispersion, per-feature, across a few dozen samples at most. The
input features have changed — genomic intervals instead of transcripts —
but the estimator has not.

#figure(
  image("../../diagrams/lecture-09/08-differential-accessibility.svg", width: 95%),
  caption: [
    Differential accessibility produces the same visual language as
    differential expression. The MA plot has identical axes
    (mean count, log-fold-change between conditions); the
    significance threshold is the same Benjamini–Hochberg-adjusted
    p-value; the underlying statistical machine is the same
    negative-binomial GLM.
  ],
) <fig:diff-access>

=== DiffBind, csaw, and the Choice of Features

Two R packages dominate the differential-binding ecosystem.

*DiffBind* (Stark & Brown 2011, Ross-Innes et al. 2012) takes a set of
peak files and BAM files, builds a consensus peakset, counts reads per
peak per sample, and hands off to DESeq2 or edgeR for the statistical
test. The wrapper handles the chromatin-specific concerns —
consensus-peakset construction, normalisation choice, blacklist
intersection — that the underlying DE engines do not know about. It is
the default choice in most labs.

*csaw* (Lun & Smyth 2014, 2016) takes a different approach: tile the
genome into overlapping windows, count reads per window per sample,
run the negative-binomial GLM on every window, and merge significant
adjacent windows into differential regions _after_ the test. This
avoids the consensus-peakset choice entirely, at the cost of a higher
multiple-testing burden. csaw is the better choice for broad marks
where the consensus-peakset decision is fragile.

The normalisation question deserves a paragraph of its own. RNA-seq
normalisation (TMM, median-of-ratios) assumes that most genes are _not_
differentially expressed, and uses that assumption to find a scaling
factor that aligns the bulk of the distribution across samples. ChIP
and ATAC routinely violate the assumption — a global loss of an active
histone mark in a knockout is a biologically meaningful whole-genome
shift, and median-of-ratios will normalise it away. DiffBind exposes
three options: library-size scaling (simplest, preserves whole-genome
shifts), TMM (assumes most peaks unchanged, normalises away global
shifts), or spike-in normalisation (requires you to spike in another
species' chromatin as a reference, then normalise to the spike-in
counts). The right choice depends on whether you expect global shifts
in your experiment.

#note[
  Differential binding inherits both the strengths and weaknesses of
  count-based DE. Strength: the statistical framework is mature, the
  multiple-testing machinery is well-understood, the per-feature
  shrinkage stabilises low-count estimates. Weakness: the design has
  to fit into a GLM, which means continuous covariates and
  batch effects need explicit modelling, and the peakset has to be
  defined before the test (csaw's window approach evades the second
  weakness but not the first).
]

=== Where Pipelines Silently Fail

Three gotchas show up repeatedly in differential-binding analyses, and
each is a silent failure — no error, just a misleading result.

*Consensus peakset choice biases the test.* If you define the consensus
as "peaks called in $>= 50%$ of samples," any peak specific to one
treatment condition is dropped before the DE test even runs, and you
under-power exactly the comparison you cared about. If you define the
consensus as "peaks in any sample," noisy single-sample false-positive
peaks contaminate the test. The compromise is something like "peaks
called in $>= N$ samples where $N$ is about half the per-condition
sample count," or csaw's window approach, which sidesteps the choice.

*Replicate count is the binding constraint.* ChIP and ATAC experiments
typically use two to three biological replicates per condition. DESeq2
and edgeR work down to two replicates per group, but statistical power
falls off sharply below five replicates per arm. If the biological
effect is subtle, no statistical sophistication will recover power that
the experimental design did not provide. Use five replicates per arm if
you can.

*Batch effects dominate.* Two samples processed on different days,
with different Tn5 lots, or by different operators, will differ
systematically — sometimes more than they differ between conditions.
The DE design has to include a batch covariate; failure to do so
attributes batch variance to the biological condition and produces
inflated significance everywhere. Always include batch in the design,
even if it is not statistically significant on its own.

#warn[
  Differential binding with an altered consensus peakset is the
  silent-failure trap. If the control has 20 000 peaks and the
  treatment has 40 000, using a "peaks in $>= 50%$ of samples"
  consensus drops half the treatment-specific peaks before the test
  runs. No downstream tool can recover what you dropped — the
  treatment-specific signal is gone before the GLM ever sees it.
  Always report how the consensus peakset was constructed and what
  fraction of single-condition peaks made it into the test.
]


== Motifs, Matched Filters, and Footprints <sec:motifs>

Peak calling gives you a list of regulatory intervals; it does not tell
you _which transcription factor_ is bound at each one. ChIP-seq with an
antibody answers that question for one TF at a time, but a typical
project has hundreds of accessible regions and dozens of candidate TFs.
The bridge is motif analysis: scan the underlying DNA sequence under
each peak for the binding-preference patterns of every candidate TF,
and identify likely binding events from sequence alone.

=== Position Weight Matrices

A transcription factor binds DNA with _fractional_ sequence specificity.
It prefers some 8-to-20 bp sequences over others, but real binding is a
distribution, not a single consensus. A TF whose preferred binding site
is `TGACTCA` (the AP-1 family core) will also bind `TGAGTCA` and
`TGACTCT` with somewhat lower affinity, and the relative affinities are
the biologically meaningful quantity.

The standard representation is a *position weight matrix* (PWM): a
$4 times L$ matrix with rows indexed by base (A, C, G, T) and columns
indexed by position in the motif. The entry $w_{b,i}$ is the log of the
ratio between the observed frequency of base $b$ at position $i$ in a
set of known binding sites and the background frequency of $b$ in the
genome,

$ w_{b,i} = log (p_(b,i)) / (q_b), $

where $p_(b,i)$ is the position-specific frequency and $q_b$ is the
background. The matrix is typically estimated from hundreds of known
binding sites — historically from SELEX or protein-binding microarray
experiments, more recently from ChIP-seq peak summits.

Scoring a candidate sequence $s = s_1 s_2 ... s_L$ against the PWM is
the sum of the per-position log-odds contributions,

$ "score"(s) = sum_(i=1)^L w_(s_i, i). $

Higher score means better match. The score is essentially a log-odds
ratio between the TF-binding model and the null background model — how
much more probable is this sequence under the TF's distribution than
under random DNA.

The PWM has a natural visualisation as a _sequence logo_: a stacked
column plot in which the height of each column is the information
content of the position (in bits) and the height of each letter inside
the column is proportional to its frequency at that position. The
information content $H_i$ at position $i$ for a four-letter alphabet
with uniform background is

$ H_i = 2 - sum_b p_(b,i) log_2 (1 / p_(b,i)) , $

which ranges from 0 (no preference; column height vanishes) to 2 bits
(only one base allowed; tall column with a single letter). Stormo and
colleagues introduced the information-content view of PWMs in 1986; the
logo presentation is Schneider & Stephens 1990.

#figure(
  image("../../diagrams/lecture-09/09-pwm-motif.svg", width: 95%),
  caption: [
    The PWM in two equivalent representations: a $4 times L$
    log-odds matrix on the left, the same information rendered as a
    sequence logo on the right. The logo's letter-stack heights are
    the per-position information content; the matrix's entries are
    the log-odds weights.
  ],
) <fig:pwm>

The major curated databases of PWMs are JASPAR (~2 000 motifs across
hundreds of vertebrates; the default first-pass resource), CIS-BP
(~5 000 motifs across 300+ species, denser because it predicts motifs
computationally for TFs without direct experimental data), and
HOCOMOCO (~800 carefully curated human and mouse motifs). The choice
between them is largely about coverage versus quality: HOCOMOCO has
the highest per-motif confidence, JASPAR has broad coverage and the
best maintained tooling, CIS-BP has the broadest coverage at the cost
of including predictions.

=== Motif Scanning as Matched Filtering

Given a PWM and a DNA region, find the positions where the motif's
score is high. The naïve algorithm is the obvious one: slide the PWM
along the sequence one position at a time, compute the score at each
position, and emit positions whose score exceeds a threshold calibrated
to a per-base p-value of $10^(-4)$ or so. Scan both strands (or
equivalently, reverse-complement the PWM and scan once). Tools — FIMO
in the MEME suite, the matchPWM family in Bioconductor, MOODS — differ
in threshold calibration and speed but compute the same underlying
statistic.

#figure(
  image("../../diagrams/lecture-09/10-matched-filter.svg", width: 95%),
  caption: [
    Motif scanning rendered as discrete cross-correlation of a
    4-channel template (the PWM) against a 4-channel one-hot encoding
    of the DNA sequence. The output track is the matched-filter
    response; peaks in the response mark candidate binding sites.
  ],
) <fig:matched-filter>

#note[
  A PWM scan is a _matched filter_ on a 4-channel signal. One-hot
  encode the DNA into a four-row binary signal (one row per base, one
  hot per column). The PWM is a $4 times L$ template. Sliding the
  template along the signal and computing the inner product at each
  offset is, definitionally, discrete cross-correlation. The output
  is high at positions where the local 4-channel signal vector
  matches the template — exactly the optimal-detection statistic for
  a known signal in additive noise that every introductory
  signals-and-systems course covers. The "log-odds PWM on one-hot
  sequence equals matched filter against a 4-channel template" is
  one of the cleanest EE-to-biology isomorphisms in the entire
  pipeline.
]

The threshold calibration is worth a moment. A typical PWM scan over
the human genome ($3 times 10^9$ positions, both strands) produces
many millions of matches at any modest score. Calibrating the threshold
to a per-base p-value of $10^(-4)$ — the standard FIMO default —
restricts hits to roughly one in 10 000 bases on average, or several
hundred thousand hits per genome per TF. That is still far more hits
than a TF actually binds in vivo. The typical ratio of PWM hits to
actual ChIP-confirmed binding sites is 10:1 or 100:1, depending on the
TF and the cell type. In vivo binding is gated by chromatin context,
cooperative TF-TF interactions, post-translational modifications, and
many other factors a PWM does not encode.

#warn[
  Motif hit does not equal TF binding. A PWM scan over open chromatin
  returns a list of positions where this TF _could_ bind given the
  sequence; whether it _does_ bind in your specific cell type and
  condition is a separate question. The most common
  over-interpretation in functional-genomics papers is treating PWM
  hits as if they were binding events. They are not. Binding is what
  ChIP-seq measures; motif hits are what PWMs predict.
]

=== ATAC Footprinting

A TF physically bound to DNA protects a ~20 bp stretch from Tn5
insertion — the enzyme cannot reach the DNA because a protein is in the
way. In ChIP-seq you only know the TF is somewhere inside a ~200 bp
peak. In ATAC-seq, with deep enough coverage, you can see the
_footprint_: a local drop in Tn5 cuts at the exact binding site,
flanked by normal open-chromatin cut density on either side. The
footprint is a coverage trough flanked by shoulders — the inverse
shape of a peak.

The problem is signal-to-noise. At any single motif instance the
footprint is buried in Poisson noise on a handful of cut counts. The
solution is signal averaging across many instances of the same motif.
Take every position relative to the motif centre — say, $-200$ to $+200$
bp — pool the Tn5 cut counts at that position across thousands of motif
instances, and the systematic footprint signal accumulates linearly
while the noise accumulates as $sqrt(N)$.

#figure(
  image("../figures/ch09/f3-footprint-averaging.svg", width: 95%),
  caption: [
    Per-site Tn5 cut profiles are too noisy to read at any individual
    motif instance; the systematic footprint emerges as a coverage
    trough only after averaging over thousands of instances. The SNR
    scales as $sqrt(N)$, which is why footprinting needs both deep
    sequencing and many motif copies.
  ],
) <fig:footprint-averaging>

#note[
  ATAC footprinting is the classic signal-averaging argument from
  every undergraduate measurements course, applied to a genome-scale
  signal. At a single motif instance, the Tn5 cut counts are buried in
  Poisson shot noise; the SNR is poor. Averaging $N$ instances of the
  same motif boosts SNR by $sqrt(N)$, which is what makes the
  footprint visible at $N = 5 000$ but not at $N = 50$. The same
  argument explains why footprinting needs deep coverage (≥ 50 M
  reads per sample) _and_ many motif instances (typically several
  thousand), not just one or the other.
]

@fig:footprint-original shows the same idea on a single TF (CTCF, the
canonical footprint example), with the per-instance profile averaged
over five thousand motif sites. The aggregated trough at the motif
centre is sharp; individual traces show essentially no structure.

#figure(
  image("../../diagrams/lecture-09/11-atac-footprint.svg", width: 95%),
  caption: [
    Aggregated Tn5 cut density relative to the centre of a CTCF
    motif, averaged over thousands of motif instances. The cut
    density drops sharply at the binding site (the footprint) and
    rebounds on either side (the flanking shoulders).
  ],
) <fig:footprint-original>

Tools — TOBIAS (Bentsen et al. 2020), HINT-ATAC (Li et al. 2019),
BaGFoot (Baek et al. 2017) — all implement the same aggregation idea
with different choices of Tn5-bias correction and statistical
thresholding. TOBIAS in particular has become the modern default
because of its explicit Tn5-bias deconvolution step and its built-in
differential-footprint analysis across conditions.

#tip[
  Footprinting is the only analysis in this chapter that detects
  signal by its _absence_ rather than its presence. Everything else
  in the chapter detects pile-ups — regions with more reads than the
  background. Footprinting detects a local deficit of reads, flanked
  by normal coverage. Finding something by what isn't there is a
  different task than finding it by what is there, which is why
  footprinting tools are separate from peak callers and why their
  statistical thresholds need separate calibration.
]


== Sequence-to-Regulation Models <sec:enformer>

Everything up to this point in the chapter is _measurement-driven_:
run an assay, observe the data, call peaks or footprints from the
observation. A different paradigm has emerged that bypasses the
measurement entirely: predict the regulatory landscape directly from
DNA sequence, using a neural network trained on the ENCODE corpus.

The flagship model is *Enformer* (Avsec et al. 2021, DeepMind in
collaboration with Calico), a CNN-plus-transformer hybrid that takes a
~100 kb DNA window as input and predicts thousands of regulatory
outputs simultaneously: CAGE signal (a proxy for promoter activity) at
every position, ATAC and DNase accessibility, ChIP-seq signal for many
TFs and many histone marks. The model is trained on the union of
ENCODE, Roadmap Epigenomics, and FANTOM5 — essentially every
high-quality public regulatory dataset that exists. Performance on
held-out genomic regions reaches within ~15 % of the experimental
replicate reproducibility for many marks, which is to say "about as
good as re-running the experiment" for a significant fraction of the
output features.

#figure(
  image("../../diagrams/lecture-09/12-enformer-architecture.svg", width: 95%),
  caption: [
    Enformer architecture. The CNN block extracts local sequence
    features from a 100 kb one-hot input, downsampling to 128 bp
    bins; the transformer block captures long-range interactions
    across the entire 100 kb window; many output heads predict
    thousands of regulatory tracks in parallel.
  ],
) <fig:enformer>

The architectural choice is worth understanding. The CNN block does
local feature extraction at sub-kilobase resolution — the same kind of
local sequence-pattern detection that a PWM scan does, but with
learned filters and many composition layers. The transformer block on
top captures the long-range interactions that a CNN alone cannot
represent: distal enhancers up to 50 kb from their target promoter,
chromatin-loop endpoints, larger regulatory domains. The
attention mechanism in the transformer learns which distant positions
matter for predicting the output at the central bin. Without the
transformer block the model would underperform on long-range
regulation; without the CNN block the model would lack the local-motif
detectors that read sequence at single-base resolution.

A successor model, *Borzoi* (Linder et al. 2023), extends the context
window to ~500 kb and adds explicit RNA-output heads (so the model
predicts not just chromatin states but actual gene expression as well).
Earlier CNN-only models — DeepSEA, ExPecto, Sei — explore the same
problem with simpler architectures.

What the sequence-to-regulation models do that MACS2 cannot:

- Predict assay output on _sequences that have never been run_. This
  includes other species, synthetic sequences, and personal-genome
  variants. The model takes any 100 kb of DNA and emits a prediction;
  it does not require that the sequence have been sequenced in a
  matching cell line.
- Score _variant effects_ by comparing model outputs for reference and
  variant sequences. This is the most-used application: take a
  candidate regulatory variant from a GWAS or a rare-disease genome,
  feed in the reference and variant 100 kb windows, and compute the
  difference in predicted ATAC or ChIP signal. The result is an
  in-silico estimate of the variant's regulatory impact.
- Capture _long-range interactions_ that the per-base analyses of the
  rest of this chapter cannot see directly. A distal enhancer's
  effect on its target promoter is part of what the model has learned.

What the models still cannot do well: predict accurately in cell types
underrepresented in training data, predict rare cell-type-specific
regulation that does not appear at scale in ENCODE, and model
cooperative-binding nonlinearities beyond what the architecture
captures. The training distribution is large but biased toward
well-studied cell types and well-studied marks; predictions in the
long tail of the regulatory landscape are correspondingly weaker.

#figure(
  image("../figures/ch09/f4-cnn-as-pwm-bank.svg", width: 95%),
  caption: [
    The first convolutional layer of Enformer is a bank of ~256
    learned 15 bp filters operating on the 4-channel one-hot DNA
    input — structurally identical to a bank of PWMs scanning for
    motifs in parallel. Many filters specialise during training onto
    specific TF motifs, with the resulting filter weights
    interpretable as learned PWMs.
  ],
) <fig:cnn-as-pwm>

#note[
  The first layer of Enformer's CNN is a bank of ~256 learned
  convolutional kernels, each 15 bp wide, operating on 4-channel
  one-hot DNA. This is structurally identical to the PWM matched-filter
  framing from @sec:motifs — except the kernels are _learned from
  data_ instead of being hand-curated from known binding sites. Each
  filter specialises during training to fire on a specific motif-like
  pattern; across the 256 filters, the bank rediscovers a large
  subset of known TF motifs (and presumably some new ones). The
  subsequent convolution layers compose these into higher-order
  patterns — spacing between motifs, cooperative binding signatures.
  The architecture sits exactly at the intersection of classical
  signal processing (matched-filter banks) and learned representations
  (data-trained kernels), and it is the cleanest example in genomics
  of how DSP intuition survives intact into the deep-learning era.
]


== Summary <sec:summary>

- ChIP-seq measures protein-DNA binding; ATAC-seq measures chromatin
  accessibility. Both produce read pile-ups at regulatory regions, and
  the analysis pipeline converts pile-ups into a list of regulatory
  intervals plus per-interval count statistics across conditions.
- Library chemistry leaves a signature in the data. The ATAC
  fragment-length distribution shows a sub/mono/di/tri-nucleosomal
  ladder when the experiment works and fails the ladder when it does
  not — the first QC gate of any ATAC pipeline.
- Peak calling is constant-false-alarm-rate detection. MACS2 estimates
  a local Poisson rate from nested windows around each candidate, takes
  the maximum as a greatest-of-CFAR statistic, and tests the observed
  count against the local null. Narrow-peak mode fits TFs and punctate
  marks; broad-peak mode fits domain-scale marks.
- Differential binding is Chapter 6 on a different feature set. The
  same negative-binomial GLM, the same dispersion shrinkage, the same
  Benjamini–Hochberg correction — applied to peaks-by-samples instead
  of genes-by-samples. DiffBind and csaw are the two standard
  wrappers.
- Motif scanning is matched filtering. A PWM is a 4-channel template;
  scanning is discrete cross-correlation of the template against a
  one-hot encoding of the DNA. Hits are far more numerous than actual
  binding events; in-vivo binding is gated by chromatin context that
  PWMs do not model.
- ATAC footprinting is signal averaging over thousands of motif
  instances. SNR scales as $sqrt(N)$; individual instances are too
  noisy to read, but the aggregate footprint emerges cleanly with
  enough depth and enough motif copies.
- Sequence-to-regulation deep learning predicts assay output directly
  from DNA. Enformer is a CNN-plus-transformer trained on the entire
  ENCODE corpus; the first CNN layer is a learned PWM bank, and the
  transformer captures long-range interactions. The architecture is
  classical DSP intuition surviving intact into the deep-learning era.


== Exercises <sec:exercises>

#strong[1.] #emph[Fragment-length QC.]
Download an ATAC-seq BAM file from ENCODE (for example, the GM12878
ATAC-seq experiment, accession ENCSR095QNB) and plot the fragment-length
distribution from positions $0$ to $1 000$ bp on a log y-axis.
Annotate the sub-nucleosomal, mono-nucleosomal, di-nucleosomal, and
tri-nucleosomal modes. Report the ratio of sub-nucleosomal to
mono-nucleosomal fragments — this ratio varies by cell type and is
informative about how decondensed the chromatin is.

#strong[2.] #emph[Local Poisson by hand.]
A candidate 200-bp window contains $c = 18$ fragments. Surrounding
windows give $lambda_(text("bg")) = 0.4$, $lambda_(1000) = 1.8$,
$lambda_(5000) = 2.5$, $lambda_(10000) = 2.2$, $lambda_(text("input")) = 0.9$
(all in fragments per 200 bp). Compute $lambda_(text("local"))$ and
the Poisson tail probability $P(X >= 18 | lambda_(text("local")))$.
Assuming a Bonferroni-style cutoff of $10^(-8)$ for a whole-genome scan
of $1.5 times 10^7$ windows, is this candidate a peak?

#strong[3.] #emph[Narrow versus broad on the same data.]
Run MACS2 in default narrow-peak mode and in `--broad` mode on the
same ENCODE H3K27me3 ChIP-seq dataset (for example, K562 H3K27me3).
Compare the number of peaks, the median peak width, and the median
peak score across the two outputs. Which mode produces output that
matches the underlying biology — single broad Polycomb-silenced
domains versus many fragments of the same domain? Justify briefly.

#strong[4.] #emph[Differential accessibility cascade.]
Download two ATAC-seq conditions (for example, control versus
LPS-stimulated macrophages from GEO accession GSE100383). Call peaks
per sample, build a consensus peakset using two different rules ("any
sample" and "$>= 3$ of 6 samples"), build the counts matrices for both
rules, and run DiffBind. How many peaks are called as differential
under each consensus rule? Where does the difference come from —
under-power on one rule, false positives on the other, or both?

#strong[5.] #emph[PWM scan.]
Pick a TF with a JASPAR motif (CTCF, MA0139.1 is the canonical
example). Scan a 100 kb region around an arbitrary gene of interest
for hits at a per-base p-value of $10^(-4)$. How many hits sit in
ATAC-peak regions versus closed-chromatin regions? Compare the hit
count to an actual CTCF ChIP-seq peak track from ENCODE for the same
region — what fraction of PWM hits are confirmed binders, and what
fraction of ChIP peaks contain at least one PWM hit?

#strong[6.] #emph[Footprint SNR.]
A footprint analysis on a single TF motif instance has a baseline
Tn5-cut density of $lambda = 0.5$ cuts per bp and a footprint depth of
$0.4 lambda$ (i.e. the bound site shows 40 % fewer cuts than the
shoulders). Assuming Poisson noise, how many motif instances must you
average to detect the footprint at a $3 sigma$ level above background?
Show the $sqrt(N)$ scaling explicitly.

#strong[7.] #emph[Enformer variant scoring.]
Pick a published regulatory variant — for example, a fine-mapped GWAS
SNP associated with a phenotype of your choice. Run Enformer
(via the public Hub model) on both the reference and the alternate
sequence centred on the variant. Compare the predicted ATAC and
H3K27ac tracks. Does the model predict a regulatory effect at the
variant position? Compare against the published functional
annotation.

#strong[8.] #emph[(Open-ended.)]
Pick one tool from the chapter — MACS3, TOBIAS, DiffBind, Enformer,
csaw, ROSE, HOCOMOCO — and read its primary publication. In one
paragraph, describe the single most surprising design choice the
authors made and why it works on the empirical data they show.


== Further Reading <sec:further-reading>

- #strong[Zhang, Y., Liu, T., Meyer, C. A., et al.] (2008).
  "Model-Based Analysis of ChIP-Seq (MACS)." #emph[Genome Biology] 9:
  R137. The original MACS paper. Reads cleanly fifteen years later;
  the detection formulation has aged well.

- #strong[Feng, J., Liu, T., Qin, B., Zhang, Y., & Liu, X. S.] (2012).
  "Identifying ChIP-Seq Enrichment Using MACS." #emph[Nature Protocols]
  7: 1728–1740. The MACS2 paper, with the local-$lambda$ refinement and
  the broad-peak mode.

- #strong[Buenrostro, J. D., Giresi, P. G., Zaba, L. C., Chang, H. Y.,
  & Greenleaf, W. J.] (2013). "Transposition of Native Chromatin for
  Fast and Sensitive Epigenomic Profiling of Open Chromatin,
  DNA-Binding Proteins and Nucleosome Position." #emph[Nature Methods]
  10: 1213–1218. The ATAC-seq paper.

- #strong[ENCODE Project Consortium.] (2012). "An Integrated
  Encyclopedia of DNA Elements in the Human Genome." #emph[Nature] 489:
  57–74. The ENCODE 2012 release. The default reference dataset for
  every method in this chapter.

- #strong[Amemiya, H. M., Kundaje, A., & Boyle, A. P.] (2019). "The
  ENCODE Blacklist: Identification of Problematic Regions of the
  Genome." #emph[Scientific Reports] 9: 9354. The blacklist paper.
  Always intersect your peak calls with the matching reference build's
  blacklist before downstream analysis.

- #strong[Ross-Innes, C. S., Stark, R., Teschendorff, A. E., et al.]
  (2012). "Differential Oestrogen Receptor Binding Is Associated with
  Clinical Outcome in Breast Cancer." #emph[Nature] 481: 389–393. The
  DiffBind paper.

- #strong[Lun, A. T. L., & Smyth, G. K.] (2016). "csaw: a Bioconductor
  Package for Differential Binding Analysis of ChIP-seq Data Using
  Sliding Windows." #emph[Nucleic Acids Research] 44: e45. The csaw
  paper.

- #strong[Stormo, G. D.] (1986). "Identifying Coding Sequences." In
  #emph[Nucleic Acid and Protein Sequence Analysis, a Practical
  Approach]. The information-theoretic foundation of PWMs.

- #strong[Castro-Mondragon, J. A., Riudavets-Puig, R., Rauluseviciute,
  I., et al.] (2022). "JASPAR 2022: The 9th Release of the Open-Access
  Database of Transcription Factor Binding Profiles."
  #emph[Nucleic Acids Research] 50: D165–D173. The current JASPAR
  release paper.

- #strong[Bentsen, M., Goymann, P., Schultheis, H., et al.] (2020).
  "ATAC-Seq Footprinting Unravels Kinetics of Transcription Factor
  Binding during Zygotic Genome Activation."
  #emph[Nature Communications] 11: 4267. The TOBIAS paper.

- #strong[Avsec, Ž., Agarwal, V., Visentin, D., et al.] (2021).
  "Effective Gene Expression Prediction from Sequence by Integrating
  Long-Range Interactions." #emph[Nature Methods] 18: 1196–1203. The
  Enformer paper.
