#import "../theme/book-theme.typ": *

= Single-Cell RNA-seq: From Droplets to Cell-Type Atlases <ch:scrna-seq>

#matters[
  Bulk RNA-seq, the subject of Chapter 5, treats a tissue sample as a
  bag of identical cells and reports an average. For a homogeneous
  culture the average is what you want. For everything else — a
  developing embryo, a tumour, a blood draw, a piece of liver — the
  average is a value that no single cell has. Cellular heterogeneity is
  the rule and not the exception, and most biology that matters lives
  in the differences between cells, not in their mean. Single-cell
  RNA-seq is the experimental and computational machinery that lets you
  read those differences. The data are sparser, the pipelines have more
  moving parts, and almost every step is statistically harder than its
  bulk equivalent. This chapter walks the whole stack: how the droplets
  are made, why a unique molecular identifier survives an amplification
  channel, how to draw a defensible line between a real cell and an
  empty droplet, why the standard normalisation tricks of bulk RNA-seq
  fail on a 95 %-zero matrix, and how PCA and UMAP and Leiden combine
  into the canonical visualisation that ends with a labelled UMAP
  hanging over a poster.
]

The technical story of single-cell RNA-seq is shorter than its
biological consequences. The first true single-cell RNA-seq paper,
Tang et al. 2009, hand-picked individual mouse blastomeres into PCR
tubes and processed them one at a time. Throughput was below a hundred
cells per experiment. Six years later the Macosko and Klein groups
independently published droplet-based methods — Drop-seq and inDrops —
that pushed throughput into the thousands of cells per run. In 2017,
10x Genomics' Chromium platform commercialised the same idea with
better chemistry and reproducibility, and within three years it had
become the dominant single-cell instrument in the world. The
computational stack matured in parallel: Cell Ranger, STARsolo, and
alevin-fry on the quantification side; Seurat (R) and Scanpy (Python)
on the analysis side. By the early 2020s a single sample of ten
thousand cells could be processed end-to-end in an afternoon, and the
Human Cell Atlas had begun building a tissue-by-tissue reference of
what every cell type in the human body looks like at the transcript
level.

The five moves of single-cell analysis are: stack the droplets,
deduplicate on UMIs, filter the obvious junk, cluster what is left,
name the cell types. Section 7.1 motivates why single-cell resolution
matters and lays out the 10x Chromium chemistry. Section 7.2 covers
the cell barcodes and unique molecular identifiers that turn a noisy
amplification channel into an invertible one. Section 7.3 walks the
quantification pipeline from BCL files to a sparse count matrix.
Section 7.4 covers the four QC filters that every pipeline applies —
empty droplets, doublets, mitochondrial fraction, ambient RNA —
because every one of them is a known way to make a downstream
analysis silently wrong. Section 7.5 covers normalisation, where the
bulk-RNA-seq toolkit breaks and a single-cell-specific toolkit takes
over. Section 7.6 walks the PCA-then-UMAP dimensionality cascade and
the graph-based Leiden clustering that runs on top of it. Section 7.7
covers cell-type annotation, where a cluster of cells becomes a
biological label.

The chapter assumes you understand bulk RNA-seq from Chapter 5 —
specifically, what a count matrix is and why library size matters.
Everything downstream of the single-cell pipeline — trajectory
inference, RNA velocity, batch integration, spatial transcriptomics —
builds on the cell-by-gene matrix this chapter produces.


== Why Single-Cell, and How the Droplets Are Made <sec:why-droplets>

Bulk RNA-seq measures expression in a population. The output is a
single expression profile per sample, the weighted mean across however
many cells were lysed into the tube. For a tissue like liver — which
contains hepatocytes, Kupffer cells, sinusoidal endothelial cells,
stellate cells, and at least half a dozen rarer types — the bulk
average is a weighted sum that loses every interesting question
about composition. *Which cell type expresses my gene of interest?*
The bulk number does not say. *Does a rare population exist that
expresses something distinct?* The bulk average dilutes a 1 %
population's signal by a factor of a hundred. *Are these cells
discrete states or a continuous gradient?* Bulk averaging cannot
distinguish.

Single-cell RNA-seq dissolves the tissue into a suspension of
individual cells, measures each cell's mRNA pool independently, and
returns a count matrix indexed by (cell, gene). For a sample of ten
thousand cells and twenty thousand genes the matrix is enormous but
sparse: a typical cell yields a few thousand non-zero gene entries
out of twenty thousand. The information that bulk loses sits in the
*differences* between rows of that matrix.

#figure(
  image("../../diagrams/lecture-07/01-bulk-vs-single-cell.svg", width: 92%),
  caption: [
    A tissue containing four distinct cell types. Bulk RNA-seq returns
    a single weighted-mean expression vector. Single-cell RNA-seq
    returns a row per cell and recovers the composition that the bulk
    average had hidden.
  ],
) <fig:bulk-vs-sc>

The empirical claim that drove the field, repeatedly confirmed since
2015, is that *almost no tissue is homogeneous*. Even cell populations
sorted to high purity by flow cytometry — CD4+ T cells, peripheral
blood mononuclear cells, sorted liver hepatocytes — turn out to
contain five or ten functionally distinct sub-states once you look
cell-by-cell. Treating the bulk average as if it described any
individual cell is the ergodic fallacy in transcriptomics:
substituting an ensemble mean for a per-sample observation. In
practice that substitution costs you almost every interesting
biological signal.

#note[
  Bulk RNA-seq is a low-resolution integrator that sums over all cells
  and reports the mean. Single-cell is a high-resolution sampler that
  resolves individual signal sources. The same trade-off appears in
  radio astronomy (beamforming averages many sky pixels versus
  interferometry that resolves individual sources), in microscopy
  (bulk fluorescence versus single-molecule imaging), and in signal
  processing (a power-spectrum estimate versus time-domain sample
  inspection). Going single-cell trades per-sample SNR for resolution.
  Most biological questions, it turns out, need the resolution.
]

=== The 10x Chromium Droplet Workflow

A modern single-cell RNA-seq experiment is fundamentally a microfluidic
encapsulation problem. The 10x Genomics Chromium platform, which has
dominated the market since 2017, isolates individual cells inside
oil-water droplets and runs reverse transcription per-droplet. The
chemistry is short to describe and consequential in its details:

1. *Load.* A cell suspension at one to fifty thousand cells per sample
   is loaded onto a microfluidic chip alongside a large excess of
   barcoded gel beads. The chip merges the three streams — cells,
   beads, oil — at a Y-junction that emits one (cell, bead, oil)
   droplet at a time. Roughly 100 000 droplets form per channel; the
   intended cell-loading rate is set so that the probability of any
   single droplet containing both a cell and a bead is around
   five to ten per cent.
2. *Lyse and prime.* Inside each droplet, the cell is chemically lysed
   and its mRNAs released into solution. Each gel bead carries
   millions of oligonucleotide primers, each consisting of three
   pieces: a *cell barcode* (16 bp; the same sequence across all
   primers on one bead, different across beads), a *unique molecular
   identifier* or UMI (12 bp; random per primer), and a poly-T tail
   (which hybridises to the poly-A of mature mRNA).
3. *Reverse transcribe.* Reverse transcriptase extends each primer
   into a cDNA copy of the captured mRNA. The cell barcode and UMI
   are now physically attached to the cDNA molecule, at its 5' end,
   *before* any amplification has happened.
4. *Pool and sequence.* The emulsion is broken, all the cDNA is
   pooled into a single tube, PCR-amplifies the pool, and the
   amplified library goes onto an Illumina sequencer. Every resulting
   read carries (cell barcode, UMI, transcript sequence). A typical
   v3 Chromium run targets ten thousand cells, fifty thousand reads
   per cell, and roughly three thousand detected genes per cell.

#figure(
  image("../../diagrams/lecture-07/02-droplet-chemistry.svg", width: 92%),
  caption: [
    The 10x Chromium chip. Three input streams (cells, beads, oil)
    converge to produce one (cell, bead, oil) encapsulation per
    droplet; reverse transcription inside the droplet labels every
    cDNA with cell-specific and molecule-specific tags.
  ],
) <fig:droplet>

The Poisson statistics of droplet loading set a hard trade-off. Load
too few cells and most droplets are empty; load too many and a
significant fraction of droplets contain two or more cells — *doublets*,
which appear in the data as fictitious mixed-cell-type entities.
Section 7.4 returns to this.

#warn[
  Chromium is not the only single-cell platform. Smart-seq2 and
  Smart-seq3 (Picelli et al. 2013; Hagemann-Jensen et al. 2020) are
  plate-based — each well holds one cell, each well is sequenced
  independently. Smart-seq protocols recover three to four times as
  many genes per cell as droplet methods, at the cost of two orders of
  magnitude less throughput and no UMIs. Drop-seq (Macosko et al. 2015)
  and inDrops (Klein et al. 2015) are the open-source droplet
  ancestors of Chromium. Microwell-seq (Han et al. 2018) and
  combinatorial-indexing methods (sci-RNA-seq, SPLiT-seq) use
  multi-round barcoding instead of droplets to scale to millions of
  cells per experiment. When reading a single-cell paper, identify the
  platform before comparing any numbers — sensitivity, gene
  detection, and noise floors are not interchangeable.
]


== Cell Barcodes and UMIs: Tags Before the Amplifier <sec:umi>

The two pieces of nucleotide tag attached to every cDNA do
fundamentally different jobs. The cell barcode says *which cell* the
molecule came from; the UMI says *which mRNA molecule within that
cell*. Both are required because the pooling step at the end of the
droplet workflow destroys cell-of-origin information, and because PCR
amplification destroys molecule-count information.

A *cell barcode* is 16 bp drawn from a fixed list of about 6 million
designed sequences that are mutually separated by at least one
Hamming-distance unit. Two consequences. First, all reads from a
single cell share the same 16 bp prefix on the read-1 side of the
library, and grouping reads by cell barcode recovers per-cell
collections. Second, single-base sequencing errors in the barcode are
correctable: any observed barcode at Hamming-1 distance from a valid
barcode can be relabelled to its nearest neighbour. (Hamming-2 errors
are usually left as-is and discarded — they cannot be unambiguously
attributed.)

A *UMI* is a 10-to-12 bp random sequence stamped on each oligo-dT
primer at manufacturing time, before the cell is even captured. The
key word is *random*. Every primer on a bead carries a different UMI;
when reverse transcription runs, every captured mRNA in the droplet
gets a different UMI as its tag. After reverse transcription the cDNA
is amplified by PCR through ten to fifteen cycles, producing tens to
hundreds of copies of every initial cDNA. After sequencing, every
read carries (cell barcode, UMI, transcript). The deduplication step
*collapses on (cell barcode, UMI)*: any group of reads sharing both
fields represents PCR copies of a single original molecule and is
counted once, not many times.

#figure(
  image("../../diagrams/lecture-07/03-barcode-umi.svg", width: 92%),
  caption: [
    Cell barcode and UMI structure. The barcode is the same across
    every primer on one bead and distinct across beads; the UMI is
    different for every primer. PCR amplifies each tagged cDNA into
    many copies that share the original tag.
  ],
) <fig:barcode-umi>

The dedupe arithmetic is worth doing carefully because it justifies
everything that follows. @fig:dedup walks one cell barcode through
three captured molecules, an uncontrolled PCR amplification factor,
and the final collapse on (cell barcode, UMI). Two hundred and nine
sequencing reads, collected because the amplification factor was
about seventy, recover the same answer as three reads would have
under perfect efficiency: this cell has two molecules of gene G1 and
one of gene G2.

#figure(
  image("../figures/ch07/f1-umi-dedup-walkthrough.svg", width: 95%),
  caption: [
    UMI deduplication on a single cell barcode. Three captured mRNAs
    receive three distinct UMIs, PCR expands each clone to a random
    multiple, and the final molecule count is the number of unique
    (cell barcode, UMI) keys — independent of how many duplicates the
    amplifier produced.
  ],
) <fig:dedup>

#note[
  A UMI is an error-correcting tag placed on every molecule *before*
  the amplification channel. The PCR step is a noisy multiplicative
  amplifier (Chapter 1's PCR-as-positive-feedback-amplifier framing),
  and without a pre-amplification tag you have no way to distinguish
  "I sequenced one molecule fifty times" from "I sequenced fifty
  distinct molecules once." The UMI is the serial number stamped on
  each input before it enters the amplifier; counting distinct serial
  numbers at the output recovers the original input count regardless
  of amplification factor. The same idea drives MAC addresses on a
  network, transaction IDs in a distributed database, and isotope
  labelling in mass spectrometry — unique identifiers stable across a
  lossy transport.
]

=== UMI Collisions and Birthday-Paradox Math

The UMI scheme works only if the UMI space is large enough that two
different molecules in the same (cell, gene) rarely receive the same
random tag. The 10x v3 chemistry uses a 12 bp UMI, giving
$4^12 approx 1.7 times 10^7$ distinct UMIs. A typical cell expresses
a few thousand mRNAs of a typical gene at most, and any one gene
sees only one or two molecules per cell. The birthday-paradox
collision probability for $n$ molecules drawn from $N = 4^12$
possibilities is approximately $n^2 / (2N)$ for small $n$; at $n = 10$
the expected collision count is about $3 times 10^(-6)$. Collisions
are negligible at the per-gene level, very rare at the per-cell level,
and noticeable but tractable at the per-dataset level.

A separate concern is *UMI sequencing error*. A 1 % per-base error
rate over a 12 bp UMI gives roughly a 12 % chance that any given UMI
read is one Hamming step away from the true UMI. If two reads from the
same cDNA pick up errors in different positions they will appear as
two different UMIs and over-count by one. Modern dedup tools (UMI-tools,
`alevin-fry`, Cell Ranger 7+) cluster UMIs within a small Hamming radius
before counting, using either directional or adjacency clustering to
collapse error neighbours back onto a single canonical UMI. The
correction recovers about 1–3 % of molecules that naive exact-match
counting would over-report.


== From BCL to a Count Matrix <sec:quantification>

The sequencer emits BCL (Binary Base-Call) files, which `bcl2fastq`
converts into FASTQ-formatted reads — one FASTQ for each end of the
paired-end read structure. The Chromium read layout is asymmetric:
read 1 carries the cell barcode plus UMI (28 bp total for v3), and
read 2 carries the actual cDNA sequence. The quantifier's job is to
align read 2 to the transcriptome (or pseudoalign it; see Chapter 5),
group the alignments by (cell barcode, UMI), deduplicate, and emit a
sparse cell-by-gene count matrix.

Three quantifiers dominate practice in 2024:

- *Cell Ranger* (10x Genomics' own). The reference pipeline; uses STAR
  under the hood with 10x-specific barcode handling and a tightly
  managed output schema. Slower than alternatives but produces the de
  facto standard outputs against which everything else is benchmarked.
- *STARsolo* (Dobin and Gingeras, part of STAR). Adds 10x / Drop-seq /
  Smart-seq barcode handling to vanilla STAR. Two to four times faster
  than Cell Ranger at similar accuracy. The default for anyone not
  locked into the 10x output schema.
- *alevin-fry* (Salmon team; Srivastava et al. 2019, He et al. 2022).
  Pseudoalignment-based — the same trick from Chapter 5's
  bulk-RNA-seq quantifiers, adapted to droplet data. Fastest of the
  three; benchmarks within 2 % of STARsolo on standard datasets.

The output is a sparse cell-by-gene count matrix, written either as
an HDF5 file (`.h5`) or as the Matrix Market triplet format (`.mtx`
plus `barcodes.tsv` plus `features.tsv`). Both formats store only the
non-zero entries, which is the right choice when 95 % of the matrix
is zero. The downstream object that analysis tools consume is the
*AnnData* object (Python, used by Scanpy) or the *SeuratObject* (R,
used by Seurat) — both are essentially "the sparse count matrix plus
per-cell and per-gene metadata in a single container."

#figure(
  image("../../diagrams/lecture-07/04-quant-pipeline.svg", width: 92%),
  caption: [
    The scRNA-seq quantification pipeline. BCL to FASTQ via
    `bcl2fastq`; FASTQ to a sparse cell-by-gene matrix via one of
    Cell Ranger, STARsolo, or alevin-fry; the matrix loaded into an
    AnnData or SeuratObject for downstream analysis.
  ],
) <fig:quant>

#tip[
  Save intermediate artifacts at every stage of the pipeline. Disk is
  cheap; re-sequencing is not. A clean pipeline records BCL, FASTQ,
  BAM (if produced), and the sparse matrix at distinct checkpoints, so
  that any downstream step can be re-run without re-aligning or
  re-sequencing. The single most common pipeline-time waste in
  single-cell analysis is re-running the full quantification stage to
  recover from a downstream bug.
]


== Quality Control: Four Filters Every Pipeline Runs <sec:qc>

The raw cell-by-gene matrix out of any quantifier overstates the
number of cells in the experiment, often by a factor of ten or more,
because the quantifier emits a row for *every barcode it sees*, not
just for every real cell. The QC stage cuts down to a defensible
working set. Four filters run on essentially every single-cell
pipeline, in this order: empty droplets, doublets, dying or stressed
cells (via mitochondrial fraction), and ambient-RNA contamination.

=== Empty Droplets and the Knee Plot

Most droplets in a Chromium run never received a cell. They still
generate reads, because they contain *ambient RNA* — fragments of
mRNA released by other cells during dissociation and floating in the
buffer — and those fragments hybridise to the droplet's bead and
produce barcoded sequencing reads. An empty-droplet barcode is not
silent. It is *low-count*: tens to a few hundred UMIs, versus
thousands to tens of thousands for a real cell.

The canonical visualisation is the *knee plot*: barcodes sorted by
total UMI count in descending order, plotted on log–log axes. Real
cells form a high plateau on the left; empty droplets form a long
low tail on the right; between them sits a steep transition zone
called the *knee*. A simple analysis draws the threshold at the
knee and keeps everything above. The number above varies by run, but
for a Chromium experiment targeting ten thousand cells the knee
typically lands around rank eight to ten thousand and separates a
plateau at $10^4$ UMIs from a tail at $10^1$ to $10^2$ UMIs.

#figure(
  image("../../diagrams/lecture-07/05-knee-plot.svg", width: 92%),
  caption: [
    The canonical knee plot. The plateau on the left is real cells;
    the long low tail on the right is empty droplets; the steep
    transition between them is the knee.
  ],
) <fig:knee>

The naive knee cut works but is wasteful. Some real cells — small
cells, dying cells, cells with naturally low transcriptional output —
have total UMI counts that fall inside the ambiguous transition zone,
and a hard knee threshold throws them out. *EmptyDrops* (Lun et al.
2019) is the standard improvement. It models the ambient RNA profile
as a multinomial distribution over genes, fits the profile from the
deep tail (where every barcode is unambiguously empty), and then
tests each ambiguous-zone barcode for whether its gene distribution
differs significantly from the ambient model. Barcodes whose gene
content is consistent with pure ambient RNA are rejected; barcodes
with significantly cell-like gene content are accepted, even at
moderate total UMI counts. The procedure typically recovers
ten to twenty per cent more real cells than a naive knee cut on the
same data.

#figure(
  image("../figures/ch07/f2-knee-plot-emptydrops.svg", width: 95%),
  caption: [
    Left: the knee plot, with the naive knee cut and the more
    permissive EmptyDrops boundary marked. Right: the multinomial
    deviation density used by EmptyDrops to test ambiguous-zone
    barcodes against the ambient profile fitted from the deep tail.
  ],
) <fig:emptydrops>

#note[
  EmptyDrops is a hypothesis test against a background model fitted
  from the dataset itself. The same idea appears in detector physics
  (signal over background, where the background is estimated from
  off-target regions), in radio astronomy (point sources versus a
  sky-noise spectrum estimated from "empty" pointings), and in
  intrusion detection (anomalous activity versus a baseline profile of
  normal traffic). The hard threshold becomes a calibrated
  false-discovery rate, and the recall gain is real.
]

=== Doublets

When the cell-loading concentration is high enough that two cells
end up in the same droplet, their mRNAs share a single cell barcode
and the resulting "cell" looks like a mixture of two cell types. This
is a *doublet*. The doublet rate scales approximately linearly with
the number of cells loaded per chip, with a 10x Chromium target of
ten thousand cells producing about a 5–8 % doublet fraction. "Load
more cells to amortise the per-cell cost" is a false economy past
roughly six per cent — the cells you lose to doublet filtering wipe
out the savings.

Doublets cluster in a characteristic place in the data. A doublet of
a T cell and a monocyte expresses both T-cell and monocyte markers
and, in the cell-by-cell-distance manifold, sits halfway between the
two pure clusters — a midpoint region that real cells of any type
rarely occupy. Doublet detection algorithms exploit this geometry:

- *Scrublet* (Wolock, Lopez, Klein 2019) generates synthetic doublets
  by summing the count vectors of random pairs of real cells, then
  trains a classifier to distinguish real cells from synthetics, then
  applies the classifier to every real cell to score doublet
  probability.
- *DoubletFinder* (McGinnis et al. 2019) is similar: simulate, project
  into PC space, and score each cell by its proximity to simulated
  doublets in k-NN sense.
- *scDblFinder* (Germain et al. 2021) extends the idea with iterative
  refinement and a gradient-boosting classifier.

All three produce a per-cell doublet score, and a hard threshold
(typically 0.3–0.4) marks cells as likely doublets. Real-world doublet
fractions match the loading-rate prediction within a percent or two
when the algorithms are calibrated against the experimental loading
concentration.

#figure(
  image("../../diagrams/lecture-07/06-doublet-detection.svg", width: 92%),
  caption: [
    Doublet detection. Synthetic doublets are generated by averaging
    pairs of real cells and projected into PC space; they land in the
    between-cluster gaps that real cells avoid. A classifier trained
    on the real-versus-synthetic boundary flags cells whose data
    geometry looks doublet-like.
  ],
) <fig:doublet>

#note[
  Doublet detection is outlier detection in a high-dimensional feature
  space. Real cells lie on a small number of low-dimensional
  manifolds (the cell-type clusters); doublets project to the
  midpoints between manifolds, a region that real cells rarely
  occupy. A classifier trained on the simulated midpoint distribution
  learns to recognise the geometry. The same structural move
  underlies anomaly detection in network traffic (where normal flows
  cluster and intrusion attempts lie in the gaps) and in image-based
  defect detection.
]

=== Dying Cells: Mitochondrial Fraction

A cell whose plasma membrane has been compromised during dissociation
loses most of its cytoplasmic mRNA into the buffer but retains its
mitochondrial RNA, which sits inside the protected mitochondrial
membrane. The signature is unmistakable: a cell with above-normal
mitochondrial-gene expression as a fraction of its total UMIs is a
*dying or stressed* cell, and including it in downstream analysis
contaminates the dataset with cells that no longer represent any
real biological state.

Healthy mammalian cells typically show 5–15 % mitochondrial-gene
expression. The standard QC threshold is anywhere from 10 to 25 %,
depending on the tissue — but the threshold is *dataset-specific*,
not universal. Cardiac myocytes, fast-twitch muscle fibres, and some
neurons have naturally elevated mitochondrial content and tolerate
20–30 % under healthy conditions. Always inspect the per-cell MT
distribution before setting a threshold; never copy a number from a
different tissue.

=== Ambient RNA

Ambient mRNA in the droplet buffer hybridises to every bead's
primers, not just those of the captured cell, so every cell's count
vector contains a contribution from the *ambient profile* — a
weighted mixture of expression from all the other cells in the
sample. The result is a low-level "everyone expresses everything"
floor that obscures real signal at modest expression levels.

Two correction tools:

- *SoupX* (Young and Behjati 2020) estimates the ambient profile from
  the empty-droplet barcodes (the same tail used for EmptyDrops),
  estimates a per-cell contamination fraction from a small set of
  manually picked, cell-type-specific marker genes, and subtracts the
  ambient contribution from each cell's counts.
- *CellBender* (Fleming et al. 2023) uses a variational autoencoder
  to jointly decompose each cell's observed counts into "true cell
  signal", "ambient background", and "technical noise" — a deep
  learning model trained on millions of cells.

In production, SoupX is the lightweight default; CellBender is run
when contamination is high (older samples, low-quality dissociations)
or when the downstream analysis is sensitive to small absolute
expression differences.

#figure(
  image("../../diagrams/lecture-07/07-qc-metrics.svg", width: 92%),
  caption: [
    Per-cell QC dashboard. Every pipeline inspects these four
    distributions — total UMIs, genes per cell, mitochondrial
    fraction, and ambient contamination — and sets per-dataset
    thresholds before proceeding to normalisation.
  ],
) <fig:qc-metrics>

#warn[
  Over-filtering is the failure mode that nobody talks about. A
  pipeline tuned for very high QC stringency throws away cells from
  underrepresented states — small cells, transitional states,
  legitimately stressed but biologically interesting cells. Always
  visualise QC distributions before applying thresholds, and always
  inspect the discarded cells to see whether the dataset's most
  interesting population just got filtered out. A clinical-trial
  analysis lost a treatment-responder signal in 2023 because the QC
  step removed a state-transition population it should have kept.
]


== Normalisation: Where the Bulk Toolkit Breaks <sec:normalisation>

Bulk RNA-seq normalisation (Chapter 5) assumes that most genes are
not differentially expressed and that library size is a smooth
function of biology. Single-cell data breaks both assumptions. A
typical scRNA-seq cell has a few thousand non-zero gene entries out
of twenty thousand total, and library sizes vary by an order of
magnitude across cells in the same dataset for purely technical
reasons. The median-of-ratios estimator that DESeq2 uses on bulk data
becomes unstable on mostly-zero count vectors: the "median" is taken
over whichever genes happen to be non-zero in *both* samples, which is
a biased subset and a noisy estimate.

#figure(
  image("../../diagrams/lecture-07/08-sparse-matrix.svg", width: 92%),
  caption: [
    A slice of a real single-cell count matrix. About 95 % of entries
    are zero; the visible structure is sparse block-diagonal patches
    where groups of cells express their type-defining genes. Most
    bulk-RNA-seq normalisation tricks fail on a matrix this sparse.
  ],
) <fig:sparse>

Two single-cell-specific normalisation methods dominate practice.

*Log-normalisation* — the Seurat `NormalizeData` default and the
Scanpy `pp.normalize_total` + `pp.log1p` pair. Two short steps:

1. For each cell, divide each gene's count by the cell's total UMI
   count, then multiply by a fixed scaling factor (commonly 10 000).
   This produces a counts-per-ten-thousand value, the single-cell
   analogue of CPM.
2. Apply the natural log of one plus the value. The pseudo-count of
   one prevents $log(0)$ and pulls zeros into the same continuous
   range as nonzero observations.

The result is a per-cell normalised expression matrix on a log scale.
Fast, simple, well-behaved for most downstream operations. The price
is a known limitation: the variance of the log-normalised values
still grows with the mean for low-count genes, because log-of-Poisson
is not a variance-stabilising transform. PCA and any other
variance-driven step will over-weight high-expression genes.

*SCTransform* (Hafemeister and Satija 2019, refined in Choudhary and
Satija 2022) is the more principled alternative. For each gene, fit a
regularised negative-binomial generalised linear model that predicts
the count as a function of total UMIs per cell, with the regularisation
sharing information across genes of similar mean expression. The
*Pearson residuals* of the model — the standardised differences
between observed counts and the model's predicted counts — are used
as the normalised values. Two consequences: variance is approximately
stabilised across expression levels, and the residual scale is
zero-centred (negative values are possible, unlike log-norm).

In 2024, log-norm is still the default in both Seurat and Scanpy
because it is fast and good enough for routine clustering and
visualisation. SCTransform is the right choice when variance
stabilisation matters: batch integration, fine-grained sub-cluster
analysis, and any case where the comparison rests on small expression
differences across cells of similar type.

#warn[
  A zero count in a single-cell matrix is not "this gene is not
  expressed." It is "this gene was not captured in this cell," which
  could mean unexpressed *or* expressed-but-undersampled. The two are
  indistinguishable from a single count column. Every single-cell
  inference downstream of normalisation has to respect this — for
  example, marker-gene tests should compare *relative* expression
  patterns across cell groups, not raw counts.
]

=== Highly Variable Gene Selection

After normalisation, the next reduction step is *feature selection*.
A typical scRNA-seq dataset has twenty thousand gene rows, but most
of them are not discriminative: they are uniformly expressed across
every cell type, or uniformly noisy, contributing only nuisance
dimensions to downstream clustering. Keeping all twenty thousand is
not just wasteful — it actively hurts clustering, because in a
high-dimensional space dominated by uninformative dimensions
nearest-neighbour relations stop being meaningful (the curse of
dimensionality, returned to in @sec:dimreduction).

The standard fix is to keep only the top *highly variable genes*
(HVGs): the two to five thousand genes whose variance across cells
exceeds what a mean-variance trend would predict. Two methods
dominate:

- *Variance-mean trend selection* (Seurat's `vst` method, Scanpy's
  default). Fit a smooth mean-variance trend across all genes (a
  trick borrowed from bulk DESeq2 dispersion estimation, Chapter 5).
  Keep the genes with the largest residual variance — those whose
  variance is much larger than the trend predicts for their mean.
- *Pearson-residual variance* (Lause, Berens, Kobak 2021;
  SCTransform's default). The genes whose Pearson residuals have the
  largest variance across cells. Mathematically the cleanest variant
  of HVG selection; computationally heavier.

The top two thousand HVGs typically capture most of the cell-type
discrimination signal. The remaining genes contribute negligible
extra information for clustering and embedding, and removing them
makes everything downstream both faster and more reliable.


== Dimensionality Reduction: PCA, UMAP, Leiden <sec:dimreduction>

After HVG selection, each cell is a point in a roughly two-thousand-
dimensional Euclidean space. This is not yet a workable representation.
Two thousand dimensions is far beyond the point at which Euclidean
geometry behaves usefully: pairwise distances concentrate, nearest
neighbours become indistinguishable from far neighbours, and any
clustering algorithm that depends on distance comparisons fails. The
single-cell analysis pipeline solves this with two successive
dimensionality reductions — a linear one (PCA) followed by a
nonlinear one (UMAP) — and runs clustering on the linear stage.

=== Why Two Reductions

A first reduction by *principal component analysis* (PCA) goes from
roughly 2000 HVGs down to roughly 30 to 50 principal components. PCA
is the singular value decomposition (SVD) of the cells-by-HVGs matrix:
its principal components are linear combinations of HVGs ordered by
how much variance each combination explains. In the canonical
single-cell pipeline, PCA does three jobs at once: it suppresses
noise (later PCs carry mostly noise), it reduces dimensionality to
something distance-comparable, and it produces a representation
on which the rest of the pipeline (k-NN graph, Leiden clustering,
UMAP) is fast.

A second reduction by *UMAP* (or t-SNE; UMAP has been the default
since around 2018) goes from 30–50 PCA components down to two
dimensions, for the sole purpose of visualisation. UMAP is a
nonlinear manifold-learning algorithm that preserves *local*
neighbourhood structure — cells that were nearest neighbours in PC
space remain nearest neighbours in the UMAP embedding — at the cost
of distorting global distances. A UMAP plot is reliable for showing
that two clusters are clusters, that one cluster contains visible
sub-structure, and that some cells sit on a continuous gradient. It
is *not* reliable for reading off how similar two clusters are by
measuring the distance between them.

#figure(
  image("../figures/ch07/f3-pca-umap-geometry.svg", width: 95%),
  caption: [
    The two-reduction cascade. PCA produces a 30-dimensional linear
    projection whose scree plot shows where the variance flattens
    out; UMAP takes those PC coordinates and embeds them in 2D for
    visualisation. Clustering happens on the PC coordinates, never on
    the UMAP plane.
  ],
) <fig:pca-umap>

The diagram in @fig:pca-umap shows the three states. Panel A is the
2000-dimensional HVG matrix where pairwise distances are
uninformative. Panel B is the 30-dimensional PCA representation, with
a scree plot whose elbow at around PC 30 motivates the cutoff and a
PC1-versus-PC2 scatter that already partly separates cell types.
Panel C is the 2-dimensional UMAP embedding where the cell-type
structure is visible at a glance.

=== PCA on Sparse Count Data

The PCA computation on a single-cell HVG matrix is straightforward:
mean-centre each gene to zero mean (sometimes also scale to unit
variance — the right choice depends on the normalisation), then
compute the SVD of the resulting matrix. The top $k$ left-singular
vectors are the cells projected onto the top $k$ principal components.

How many components to keep? The standard heuristic is *enough to
capture 50–80 % of the variance*, typically 30–50. The scree plot —
eigenvalues plotted against component index — shows a steep early
drop followed by a long flat plateau, and the elbow between the two
is the natural cutoff. A more principled approach is *parallel
analysis*: compare the scree plot to one generated by shuffling each
gene's expression across cells (destroying the cell-cell signal but
preserving each gene's marginal distribution), and keep only the
components whose eigenvalues exceed the shuffled baseline. Parallel
analysis is the right tool when the dataset has unusual structure
and the "30–50" heuristic fails.

#figure(
  image("../../diagrams/lecture-07/09-pca-projection.svg", width: 92%),
  caption: [
    PCA output on a single-cell dataset. The scree plot identifies
    the elbow at PC30; the PC1-versus-PC2 projection already separates
    the major cell types as visible blobs in 2D.
  ],
) <fig:pca>

#note[
  PCA is the SVD of the normalised expression matrix, and it produces
  the optimal linear rank-$k$ approximation in the Frobenius-norm
  sense. The same operation underlies JPEG (DCT-basis truncation),
  MP3 (modified-DCT truncation), eigenfaces (face-image principal
  components), and almost every dimensionality reduction in classical
  signal processing. The "top $k$ PCs capture most of the variance"
  story is the matrix-norm analogue of "the first $k$ DFT
  coefficients carry most of the energy of a smooth signal."
]

=== UMAP and t-SNE for Visualisation

Once each cell is a point in PC space, the visualisation problem is
how to lay out the cells on a 2D page so that meaningful relationships
survive. *t-SNE* (van der Maaten and Hinton 2008) was the field
default for a decade. It defines pairwise similarities in the
high-dimensional space using Gaussian kernels and pairwise similarities
in the low-dimensional space using $t$-distributions, and minimises a
Kullback-Leibler divergence between the two distributions. The
$t$-distribution's heavy tails are a feature: it allows distant
high-dimensional points to be placed far apart in 2D without paying
an enormous penalty.

*UMAP* (McInnes, Healy, Melville 2018) is the more recent and now
dominant alternative. UMAP treats the data as a fuzzy topological
structure (a weighted $k$-nearest-neighbours graph in the
high-dimensional space) and finds a low-dimensional layout that
preserves the same topology. In practice UMAP is two to five times
faster than t-SNE, preserves global structure noticeably better
(distances between distant clusters are at least roughly
informative), and has become the de facto single-cell visualisation
standard.

#figure(
  image("../../diagrams/lecture-07/10-umap-embedding.svg", width: 92%),
  caption: [
    A UMAP embedding of a PBMC-like dataset. Clusters, sizes, and
    putative cell-type assignments are all visible at once — which
    is why UMAP, not anything more rigorous, is the figure that
    appears on every single-cell paper.
  ],
) <fig:umap>

#warn[
  UMAP distances above the local-neighbour scale are not metric.
  Two clusters that sit far apart in a UMAP plot are not necessarily
  far apart biologically; UMAP is allowed to push apart any two
  groups that lack a continuous bridge of intermediate cells. Run
  clustering on the PCA coordinates, not on the UMAP coordinates,
  and use UMAP as a presentation layer, not as a metric space.
]

=== Graph-Based Clustering: Leiden

Clustering on 30-dimensional PC coordinates is done graph-theoretically,
not with $k$-means or hierarchical clustering. The reason is that
single-cell data lives on low-dimensional manifolds embedded in the
PC space, and Euclidean-distance-based clustering algorithms cut
across the manifolds in arbitrary ways. Graph-based methods respect
the local geometry by construction.

The pipeline is three steps:

1. Build the *k-nearest-neighbours graph* (kNN, typically $k = 10$ to
   $30$). Each cell is a node; each cell is connected to its $k$
   nearest neighbours in PC space.
2. Optionally refine to a *shared-nearest-neighbours graph* (SNN), in
   which edge weights are the Jaccard overlap of $k$-NN neighbour
   sets. SNN denoises the raw kNN graph by giving more weight to
   pairs of cells that share many neighbours.
3. Run *Leiden* or *Louvain* community detection on the graph. Both
   algorithms maximise modularity — a measure of how much the graph's
   edges concentrate within communities versus crossing between them.

Modularity is defined as
$ Q = (1 / (2m)) sum_(i,j) [A_(i j) - (k_i k_j) / (2m)] delta(c_i, c_j) $
where $A$ is the adjacency matrix, $k_i$ is the degree of node $i$,
$m$ is the total edge count, and $delta(c_i, c_j)$ is 1 when $i$ and
$j$ are in the same community and 0 otherwise. The first term in the
bracket counts within-community edges; the second is the expected
within-community edge count under a degree-preserving random model.
Modularity is therefore the *excess* of within-community edges over a
null random model.

Louvain (Blondel et al. 2008) maximises $Q$ by iteratively moving
nodes between communities to whichever move yields the largest local
gain in $Q$, then aggregating communities into super-nodes and
repeating. It is fast and produces good partitions on most graphs.
Its known failure mode is that some communities it produces are
internally disconnected — Louvain happily places two unconnected
subgraphs into the same community as long as the within-community
modularity is high.

Leiden (Traag, Waltman, van Eck 2019) fixes the disconnection
problem with a refinement step that splits any disconnected community
before aggregation. The result is provably well-connected
communities, generally higher final modularity, and slightly higher
runtime. By 2024 Leiden has replaced Louvain as the single-cell
default in both Scanpy and Seurat.

#figure(
  image("../figures/ch07/f4-leiden-modularity.svg", width: 95%),
  caption: [
    Leiden community detection. Panel A is the raw $k$-nearest-neighbour
    graph in PC space. Panel B is a Louvain partition in which one
    community is internally disconnected (a known failure mode).
    Panel C is the Leiden refinement, which splits the disconnected
    community into two well-connected sub-communities.
  ],
) <fig:leiden>

A critical hyperparameter is the *resolution* $gamma$. Higher
$gamma$ produces more, smaller clusters; lower $gamma$ produces
fewer, larger clusters. The choice is not automatic. The right
practice is to scan resolutions from 0.1 to 2.0, plot the resulting
clusterings on the UMAP, and pick the resolution at which the
cluster boundaries align with known biology (e.g. cluster splits
that reveal distinct marker-gene signatures rather than splits
between cells that share all their markers). Resolution should
follow the biology, not the other way around.

#note[
  Graph-based community detection on a $k$-NN graph is structurally
  identical to modularity-based community detection in social
  networks, citation networks, and bibliometric graphs. The original
  Louvain algorithm was developed for online social networks; its
  use in single-cell is a direct import from that literature.
  Modularity itself, defined by Newman and Girvan in 2004, is the
  same quantity used for everything from Wikipedia article-link
  communities to gene co-expression networks.
]

=== Marker Genes

After clustering, the final step before annotation is *marker gene
identification* — finding genes that are differentially expressed in
one cluster versus the rest. For each cluster $c$ and each gene $g$,
run a differential-expression test comparing cells in $c$ against
all cells not in $c$. The output is a log-fold-change and an adjusted
$p$-value per (cluster, gene); the top genes by log-fold-change are
the cluster's markers.

Two DE tests dominate:

- *Wilcoxon rank-sum* (Seurat and Scanpy defaults). Nonparametric,
  fast, makes no assumption about the count distribution. Robust to
  the weird mean-variance relationship of single-cell counts.
- *MAST* (Finak et al. 2015) is a hurdle model: it jointly models
  the *fraction* of cells expressing the gene at all and the
  *expression level* in the cells that do express it. Slower and
  occasionally more sensitive to low-count markers.

Wilcoxon is the default for routine analysis. MAST is the right
choice when low-count markers matter — for example, when the
distinguishing feature of a sub-population is the absence of a gene
that the surrounding clusters express.

#figure(
  image("../../diagrams/lecture-07/11-marker-heatmap.svg", width: 92%),
  caption: [
    A marker-gene heatmap. Rows are cells grouped by cluster
    assignment; columns are the top five marker genes per cluster.
    The block-diagonal pattern is what a well-resolved clustering
    looks like — each cluster's cells light up its own markers and
    are dim elsewhere.
  ],
) <fig:markers>

#tip[
  If no gene reliably distinguishes cluster A from cluster B, those
  two clusters are probably over-split at the current Leiden
  resolution. Marker-gene inspection is the most reliable empirical
  check on cluster validity — better than any internal cluster
  quality metric. The single most useful early sanity check after a
  scRNA-seq clustering is the markers heatmap.
]


== Cell-Type Annotation <sec:annotation>

The output of the clustering and marker-finding stages is a partition
of cells into clusters plus a ranked list of markers per cluster.
Cell-type annotation is the step where each cluster gets a biological
label — "CD4 T cell", "intermediate monocyte", "endothelial cell" —
that other researchers can use to compare your dataset to theirs.
Two paradigms dominate.

*Marker-based annotation* (manual). For each cluster, look at the top
markers, compare them to a reference of known cell-type-defining
genes, and assign the cluster the label whose markers best match.
Canonical human peripheral-blood markers include MS4A1 (CD20),
CD79A, and CD19 for B cells; CD8A, CD8B, GZMB, and PRF1 for cytotoxic
T cells; CD14, LYZ, and FCN1 for classical monocytes; GNLY, NKG7, and
KLRD1 for NK cells. Marker-based annotation works well when the cell
types are well-characterised and badly when the dataset contains
rare or novel populations that no published reference covers.

*Reference-based annotation* (automated). Given a pre-annotated
reference atlas — typically Tabula Sapiens, an Azimuth-curated
reference for a specific tissue, or a published study's deposited
labels — compute the most-similar reference cell for each query cell
and transfer the reference label. Three tools dominate:

- *SingleR* (Aran et al. 2019) correlates each query cell against the
  mean expression profile of every reference cell type and assigns
  the type with the highest correlation. Fast, stable; the R
  ecosystem default.
- *Azimuth* (Hao et al. 2021, Stuart et al. 2022) is a web
  application running pre-trained references for the most common
  tissues (PBMC, bone marrow, heart). Upload a dataset, get
  annotations in twenty minutes.
- *scArches* / *scPoli* (Lotfollahi et al. 2022) use deep learning to
  map a query dataset into the latent space of a pre-trained
  reference, then transfer labels. Slower; better at handling cell
  states not present in the reference.

Reference-based methods are faster and require less manual labour.
They fail silently when the query contains cell types that the
reference does not cover — every query cell gets assigned to *some*
reference type, even if the right answer is "this cell type is not in
your reference." Validate every reference-based annotation against
marker-gene checks. If a cluster labelled "activated NK cell" has
neutrophil markers as its top differential genes, the annotation is
wrong and the right answer is that the reference did not include
neutrophils.

=== Cell Ontologies and Atlases

For inter-study consistency, cell-type annotations should reference
a shared controlled vocabulary. Three resources matter in 2024:

- *Cell Ontology* (CL; Bard et al. 2005, ongoing). The formal
  ontology of animal cell types. Every type has a CL identifier
  (e.g. CL:0000236 for "B cell"). Hierarchical: parent-child
  relationships between general and specific types let downstream
  tools query at the resolution of interest.
- *Human Cell Atlas* (HCA Consortium 2017 onwards). A growing
  community effort to map every cell type in every human tissue. By
  2024 the HCA had published roughly thirty tissue-scale atlases and
  begun integrating spatial and multi-omic data.
- *Tabula Sapiens* / *Tabula Muris* (Tabula Sapiens Consortium 2022,
  Tabula Muris Consortium 2018). Single-institution, pan-tissue
  atlases with internally consistent annotation. Smaller in scope
  than HCA but more uniform.

Use CL identifiers in your dataset metadata. Vocabulary consistency
is what makes cross-study comparison and meta-analysis tractable.

#note[
  The Human Cell Atlas, launched in 2017, has changed how the field
  defines what a cell type *is*. For most of the twentieth century the
  reference was histology: a cell was named by its morphology, its
  staining pattern, and its anatomical location. The HCA-era
  reference is molecular: a cell is what its transcriptional
  signature says it is. The two definitions agree on the major
  classes — hepatocytes, T cells, neurons — but disagree at the
  fine-grained level, where transcriptomically distinct sub-states
  often share morphology and where morphologically distinct cells
  sometimes share transcriptional programs. For working bioinformatics
  engineers in 2024, the molecular reference is the one to build
  against.
]


== Summary <sec:summary>

- Single-cell RNA-seq measures per-cell transcript abundance. Bulk
  averages hide cell-type composition; single-cell resolves it. The
  cost is per-cell signal-to-noise; the payoff is that heterogeneity
  becomes visible.
- UMIs invert the PCR amplification channel. Each mRNA is tagged
  with a random sequence before amplification; final counts are the
  number of distinct (cell barcode, UMI) keys per gene, not the
  number of reads. The trick is essential at the amplification levels
  single-cell library prep uses.
- Quality control is not optional. Four filters — empty droplets,
  doublets, mitochondrial fraction, ambient RNA — run on every
  serious pipeline. They typically remove twenty to forty per cent of
  observed barcodes before downstream analysis even begins.
- Normalisation must respect sparsity. Bulk methods (median-of-ratios)
  fail on mostly-zero matrices. Log-normalisation or Pearson
  residuals, followed by highly variable gene selection, is the
  modern default.
- Dimensionality reduction is two-stage. PCA produces a linear
  30-to-50-dimensional representation that suppresses noise and
  enables fast nearest-neighbour computation. UMAP takes the PC
  coordinates to two dimensions for visualisation only. Cluster on
  PCs, never on UMAP coordinates.
- Leiden community detection on a $k$-NN graph is the canonical
  clustering algorithm. It maximises modularity and guarantees
  well-connected communities; the resolution parameter is tuned to
  biology, not to a target cluster count.
- Annotation is marker-based or reference-based. Reference-based is
  fast and automated; marker-based is manual but catches reference
  gaps. The right practice is to do both and check that they agree.


== Exercises <sec:exercises>

#strong[1.] #emph[UMI deduplication by hand.]
Three reads share cell barcode `CB42` and gene `GAPDH` with UMIs
`AAACGTTGCAGT`, `AAACGTTGCAGT`, and `AAACGTAGCAGT`. Two of the UMIs
are exactly equal; the third differs from them by one base. (a)
Under exact-match deduplication, how many distinct UMI groups do
the three reads collapse into? (b) Under directional UMI clustering
that allows Hamming-1 collapse onto the more-supported neighbour,
how many distinct UMI groups remain? (c) Argue in one sentence which
answer is more likely to reflect the truth.

#strong[2.] #emph[UMI collision probability.]
A 10x v3 chemistry uses 12 bp UMIs. (a) Compute the size of the UMI
space $N = 4^12$. (b) For a cell expressing 1000 distinct mRNAs of
one gene, compute the expected number of UMI collisions using the
birthday-paradox approximation $n^2 / (2N)$. (c) Why is the
collision count per *gene per cell* the right unit, rather than per
dataset?

#strong[3.] #emph[Reading a knee plot.]
A knee plot of a Chromium run shows a plateau of 7 500 barcodes at
about 12 000 UMIs each, a knee around rank 8 000, and a long tail
extending to rank 100 000 at counts near 30 UMIs each. (a) Estimate
how many cells you would call under a naive knee cut. (b) Estimate
how many additional ambiguous-zone barcodes EmptyDrops typically
recovers (use the 10–20 % rule of thumb). (c) For each of three
known failure modes — a high-doublet sample, a stressed sample with
many low-UMI cells, and a sample with very heterogeneous cell sizes —
predict in one sentence what the knee plot would look like.

#strong[4.] #emph[Doublet rate from loading concentration.]
The 10x Chromium specification states that loading 16 500 cells
yields roughly 10 000 captured cells with a 7.6 % doublet rate, while
loading 33 000 cells yields 20 000 captured cells with a 15.4 %
doublet rate. (a) Verify that the doublet rate scales approximately
linearly with the loading concentration. (b) If you set a doublet
rate budget of 5 %, what loading concentration should you target?
(c) What is the cost trade-off in terms of cells captured per dollar
of sequencing?

#strong[5.] #emph[Why bulk normalisation fails.]
A toy two-cell dataset has six genes. Cell 1 counts are
`[10, 0, 0, 0, 5, 0]`; cell 2 counts are `[0, 8, 0, 0, 0, 12]`.
(a) Compute the DESeq2 median-of-ratios size factor for each cell.
(b) Comment on the result: is the estimator well-defined? (c)
Repeat with log-normalisation (divide each cell by its total UMI,
multiply by 10 000, add 1, take log). Is the answer well-defined now?

#strong[6.] #emph[PCA on a single cell.]
Suppose a single-cell dataset has 10 000 cells and 2 000 HVGs. PCA
produces eigenvalues for the first ten components of
$[2.3, 1.8, 1.4, 1.0, 0.7, 0.45, 0.30, 0.22, 0.18, 0.15]$. (a) What
fraction of the total variance do the first five PCs explain (assume
the remaining PCs sum to about 5 in total)? (b) Where does the elbow
appear to be? (c) Why might you keep PCs past the elbow even if they
individually capture less than 1 % of variance?

#strong[7.] #emph[Modularity by hand.]
A six-node graph has the following undirected edges:
$(1, 2), (1, 3), (2, 3), (3, 4), (4, 5), (4, 6), (5, 6)$. Total
edges $m = 7$. (a) Compute the modularity $Q$ of the partition
${{1, 2, 3}, {4, 5, 6}}$. (b) Repeat for the partition
${{1, 2, 3, 4}, {5, 6}}$. (c) Which partition does Leiden prefer,
and why?

#strong[8.] #emph[Annotation conflict.]
SingleR labels a cluster as "CD4 T cell" with 92 % confidence; the
manual marker check shows the cluster's top markers are CD68, LYZ,
and S100A8 (classical monocyte markers). (a) In one sentence,
explain what is going on. (b) Which annotation should you trust?
(c) What pipeline check would have caught this earlier?

#strong[9.] #emph[(Open-ended.)]
Pick one tool from this chapter — Cell Ranger, alevin-fry,
EmptyDrops, Scrublet, SCTransform, Leiden, Azimuth — and read its
primary publication. In one paragraph, identify the single most
consequential design choice the authors made and the empirical
evidence they offer that it works.


== Further Reading <sec:further-reading>

- *Macosko, E. Z., Basu, A., Satija, R., et al.* (2015). "Highly
  parallel genome-wide expression profiling of individual cells using
  nanoliter droplets." _Cell_ 161: 1202–1214. The Drop-seq paper —
  the open-source predecessor of every commercial droplet-based
  single-cell platform.
- *Klein, A. M., Mazutis, L., Akartuna, I., et al.* (2015). "Droplet
  barcoding for single-cell transcriptomics applied to embryonic
  stem cells." _Cell_ 161: 1187–1201. inDrops; published the same
  week as Macosko et al. by a competing group.
- *Zheng, G. X. Y., Terry, J. M., Belgrader, P., et al.* (2017).
  "Massively parallel digital transcriptional profiling of single
  cells." _Nature Communications_ 8: 14049. The 10x Chromium paper
  and the commercial breakthrough.
- *Lun, A. T. L., Riesenfeld, S., Andrews, T., et al.* (2019).
  "EmptyDrops: distinguishing cells from empty droplets in
  droplet-based single-cell RNA sequencing data." _Genome Biology_
  20: 63. The multinomial-deviation test for the cell-empty
  boundary.
- *Hafemeister, C., & Satija, R.* (2019). "Normalization and variance
  stabilization of single-cell RNA-seq data using regularized
  negative binomial regression." _Genome Biology_ 20: 296. The
  SCTransform paper.
- *Traag, V. A., Waltman, L., & van Eck, N. J.* (2019). "From Louvain
  to Leiden: guaranteeing well-connected communities." _Scientific
  Reports_ 9: 5233. The Leiden refinement of Louvain modularity
  optimisation.
- *McInnes, L., Healy, J., & Melville, J.* (2018). "UMAP: Uniform
  Manifold Approximation and Projection for dimension reduction."
  _arXiv:1802.03426_. The UMAP paper.
- *Wolf, F. A., Angerer, P., & Theis, F. J.* (2018). "Scanpy: large-
  scale single-cell gene expression data analysis." _Genome Biology_
  19: 15. The Python single-cell ecosystem.
- *Stuart, T., Butler, A., Hoffman, P., et al.* (2019). "Comprehensive
  integration of single-cell data." _Cell_ 177: 1888–1902. Seurat v3.
- *The Tabula Sapiens Consortium* (2022). "The Tabula Sapiens: a
  multiple-organ, single-cell transcriptomic atlas of humans."
  _Science_ 376: eabl4896. The reference atlas most modern
  annotation tools build against.
