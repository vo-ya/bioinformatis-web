#import "../theme/book-theme.typ": *

= Data Engineering, File Formats, and Reproducibility <ch:data-engineering>

#matters[
  Every algorithm in the preceding thirteen chapters assumes a piece of
  infrastructure: the read is in a file, the file lives somewhere, the
  reference is the same one the caller was trained on, the tool is
  installed, the pipeline restarts after the cluster fails, and another
  group can rerun the analysis a year from now and get the same answer.
  None of that is automatic. The bytes have to be laid out so common
  queries are cheap, the reference has to be pinned by hash, the tool
  has to be packaged so it survives the next operating-system upgrade,
  and the pipeline has to be expressed in a language that can resume
  from a half-finished state. Skipping any one of those moves is how a
  good analysis becomes an unreproducible one. This chapter is the
  engineering layer; it carries less mathematics than its neighbours
  and more file format.
]

The preceding chapters were about algorithms — how to align reads
(Chapter 2), assemble them (Chapter 3), call variants (Chapter 4),
quantify expression (Chapters 5 through 8), and so on. Every one of
those chapters quietly assumed the data was already on disk, in a
sensible format, with the right reference, accessible to a tool
installed on a machine that could find it. This chapter is about
those assumptions. It is the bioinformatics equivalent of the
electrical engineer's lab-bench discipline: what file formats hold
the bytes, where the bytes live, how to fetch them, how to package the
tools that read them, how to wire those tools into a pipeline that
restarts after failure, and how to prove that the result is the same
when somebody else reruns it on different hardware a year later.

Seven moves carry the workflow forward. Pick a format whose layout
supports the queries you intend to run; fetch the data from the
repository that hosts it; pin the reference genome by hash, not by
filename; ship every tool in a container; express the pipeline as a
directed acyclic graph; benchmark the output against a community
truth set; and archive the whole stack so a future reader can rerun
it. The chapter takes each move in turn. By the end you should be
able to read a published methods paragraph and reconstruct which
formats, repositories, references, containers, workflow language, and
benchmarks the authors implicitly relied on — and which they got
wrong.


== Why File Formats Matter <sec:why-formats>

A bioinformatics file format is a contract between a producer and a
consumer about how bytes are laid out. The interesting contracts say
more than "here are the data": they specify which queries are
cheap and which are expensive. A format that can stream from the
sequencer to the next tool but cannot answer "what reads cover this
500-base window" is fine as a transport format and useless as a
working one. A format that can answer that question in microseconds
on a 100-gigabyte file is the format you keep around for years.

Four engineering dimensions distinguish genomics formats from each
other.

- *Random access.* Given a 100 GB file, can you jump to the region
  `chr7:140753336-140753436` without reading the whole file? Cheap
  random access is the difference between a file you keep cold and a
  file you actively use.
- *Compression ratio.* How many bytes does the format spend per
  biological observation? The naive baseline is the gzipped text
  representation. Good formats beat that by an order of magnitude.
- *Schema evolution.* Can you add a new annotation field without
  invalidating every tool that reads the older variant? The formats
  that survived two decades all have a good story here; the ones
  that didn't, don't.
- *Tool ecosystem.* How many independent libraries can read it, and
  how aggressively are they maintained? A format with one
  authoritative reader is a format that dies when its maintainer
  retires.

#figure(
  image("../../diagrams/lecture-14/01-format-comparison.svg", width: 95%),
  caption: [
    File-format footprints for one 30× human WGS sample. Compression
    is only part of the story — random access via `bgzip` plus index
    is what makes #idx("CRAM")CRAM and #idx("BCF")BCF usable where flat gzip is not.
  ],
) <fig:format-comparison>

The bioinformatics field has accumulated about thirty years of
formats. Some are elegant — #idx("VCF")VCF, #idx("HDF5")HDF5, #idx("Zarr")Zarr, #idx("Parquet")Parquet. Some are
regrettable but entrenched — the #idx("GFF")GFF/GTF split, the #idx("FASTQ")FASTQ header
zoo. Most working time goes into moving data between them, so it
pays to know what each one costs.


== FASTQ, #idx("BAM")BAM, and CRAM <sec:read-formats>

The first three formats every bioinformatician learns are the three
formats that carry sequencing reads from the instrument to the
caller. They form a natural progression from "transport" through
"working" to "archival".

=== FASTQ: the transport format

FASTQ (Chapter 1) is what comes off the sequencer: four lines per
read, identifier, sequence, separator, quality #idx("STRING")string. Plain ASCII,
human-readable, near-universal as an input. It is almost always
shipped as `.fastq.gz` — generic gzip cuts the size by about a
factor of three over raw text. A 30× human whole-genome sample at
2×150 bp produces roughly 90 gigabases of reads, which is around
45 GB of `.fastq.gz`.

Two limitations matter. FASTQ has no index — you cannot ask "show me
the read pair with this identifier" without scanning the file. And
FASTQ carries no alignment information; the reads have not been
placed against a reference yet. FASTQ is a *transport* format,
optimal for moving raw reads from the sequencer to the first tool
that needs them.

#tip[
  #idx("paired-end")Paired-end #idx("Illumina")Illumina data ships as a pair of FASTQ files, `_R1.fq.gz`
  and `_R2.fq.gz`. The pairing is positional — record _n_ in `R1` is
  mated to record _n_ in `R2`. Tools that reorder one file without the
  other silently scramble the pairing. `seqkit pair` and `fastp` are
  the canonical paired-aware repair tools.
]

=== BAM: the working format

BAM (the binary form of #idx("SAM")SAM, the Sequence Alignment / Map format,
specified by Heng Li and colleagues in 2009) is the working
representation of aligned reads. Each record carries the read's
identifier, its alignment position on the reference, the #idx("CIGAR")CIGAR
string describing matches and gaps, the mate's position, the
base sequence, the quality string, the mapping quality, and a list
of optional tags. The records are sorted by genomic coordinate and
the file is *bgzip*-compressed — block-gzip with 64-kilobyte blocks,
each block independently decompressible. An adjacent `.bai` index
maps genomic regions to the byte offsets of the blocks that contain
them, so a tool can jump from `chr7:140753336` to the right block
without reading anything in between.

BAM was the de-facto standard for fifteen years. Its weakness is
that the compression is generic zlib — the same thing you would get
on any random binary file, ignoring everything the encoder knows
about the data. Reads from a sequenced genome are massively
redundant with each other and with the reference: at 30× #idx("coverage")coverage,
the same true base is observed thirty times, and the read mostly
matches the reference it came from. BAM compresses none of that
redundancy explicitly.

=== CRAM: content-aware compression

CRAM (Hsi-Yang Fritz, Leinonen, Cochrane, and Birney, 2011; refined
into the EBI standard by 2013) exploits exactly the redundancy BAM
ignores. The encoder knows what an aligned read *is* — a position on
a reference, plus a small list of differences from that reference,
plus a quality string. So it stores exactly that, column by column,
rather than storing the full base sequence per record.

A single CRAM record stores three things: per-read metadata (name,
flags, position) in columnar compressed blocks; the *differences*
from the reference at this read's mapped position; and the quality
string, which is now the dominant cost. The encoder splits each of
these into its own stream and applies a codec tuned for its
statistics — `gzip` for the metadata, a custom diff encoder for the
substitutions, and rANS (range asymmetric numeral systems) for the
quality scores.

#figure(
  image("../../diagrams/lecture-14/02-cram-scheme.svg", width: 95%),
  caption: [
    CRAM as content-aware lossless compression. The encoder knows
    what a read is — a position plus diffs from the reference — and
    stores exactly that, column by column. Lose the reference, lose
    the data.
  ],
) <fig:cram-scheme>

#note[
  CRAM is to BAM what JPEG is to BMP. BAM treats its input as an
  opaque byte stream and finds whatever patterns it can. CRAM is
  *content-aware*: the encoder knows the data is a stream of aligned
  reads against a known reference, and it stores only what is new.
  The same move separates H.264 from uncompressed video, FLAC from
  raw PCM audio, and Parquet from CSV. Wherever the data has known
  structure, the codec that exploits it beats the codec that does
  not by an order of magnitude.
]

For a 30× human sample, lossless CRAM is typically three to four
times smaller than BAM. Turn on *quality binning* — collapse the 40
distinct Phred scores into a coarser set of four or eight bins
(Illumina's Q-binning, recommended in the EBI guidelines) — and the
size drops by another factor of two. The arithmetic is concrete in
@fig:cram-arithmetic: 80 GB of BAM compresses to 25 GB of lossless
CRAM and 14 GB of Q-binned CRAM, with the savings coming almost
entirely from two streams (the sequence stream collapses by a factor
of twelve; the quality stream by a factor of three with binning).

#figure(
  image("../figures/ch14/f1-cram-compression-arithmetic.svg", width: 95%),
  caption: [
    CRAM arithmetic for a single 30× human whole-genome sample.
    Reference-based encoding collapses the sequence stream; a
    specialised codec collapses the quality stream. The savings are
    real money at biobank scale: 100,000 samples at 50 GB BAM cost
    5 PB; the same samples at 14 GB Q-binned CRAM cost 1.4 PB.
  ],
) <fig:cram-arithmetic>

The cost is real and worth understanding. CRAM is *reference-
dependent*: the encoder records the MD5 hash of the reference used
to compress the file, and decompression requires the same reference
to be present. With the correct reference, decompression is fast and
exact. Without it, a CRAM file is unreadable bytes.

#warn[
  In 2018 the Wellcome Sanger Institute lost or temporarily corrupted
  more than 50 TB of archival CRAM data due to reference-management
  failures. The CRAM files themselves were intact; the chain of
  custody from compressor MD5 to decompression-time reference had
  broken. "We compressed it with CRAM" is not an archival strategy
  unless you also archived the exact reference. Tools like
  `samtools` consult `REF_CACHE` and `REF_PATH` environment
  variables and the EBI reference-fetch server to recover the
  reference at decode time, but the system is fragile against
  network outages and registry changes. Always co-archive the
  reference with the CRAM, and pin both by MD5.
]


== VCF, BCF, and the Tabular Variant Formats <sec:vcf-bcf>

Once reads have been aligned and called against a reference, the
output is variant data. Variants are sparse — only a few million
positions in a 3.1-gigabase genome differ from the reference per
individual — so the natural representation is a table of differences,
not a per-position rewrite. That is what VCF gives you.

VCF (the Variant Call Format, introduced by Petr Daneček and
colleagues in 2011 as part of the 1000 Genomes Project) is a
tab-delimited text file with eight fixed columns —
`CHROM POS ID REF ALT QUAL FILTER INFO` — followed by a FORMAT
column and one column per sample. The semantics are dense and were
covered in Chapter 4; what matters here is the engineering
discipline that makes VCF usable at scale.

Working VCF is always *bgzipped* and *tabix-indexed*. `bgzip` is the
same block-gzip used by BAM: it cuts the file into independently
decompressible 64-kilobyte blocks. `tabix` builds an index keyed on
the `CHROM` and `POS` columns, allowing region queries — "give me
every variant on `chr17` between 41,200,000 and 41,400,000" — to
seek directly to the relevant blocks. Without the index, the file is
still readable, but the random-access guarantee disappears.

#figure(
  image("../../diagrams/lecture-14/03-vcf-structure.svg", width: 95%),
  caption: [
    A VCF record unpacked, with its BCF binary counterpart below.
    Identical information, different ergonomics: VCF is text for
    inspection; BCF is binary for scale.
  ],
) <fig:vcf-structure>

BCF is the binary companion. Same data, fixed-width integer encoding
for the numerical fields, about three times smaller and five to ten
times faster to parse. The trade-off is exactly the same as BAM
versus SAM: BCF is what tools talk to each other in, VCF is what
humans read. For cohorts above a thousand samples, BCF is
effectively mandatory; reading 500,000 sample columns of text per
row burns more CPU on parsing than on analysis.

The pain point is that VCF/BCF scales badly along the sample axis. A
UK-Biobank-scale joint-called VCF with 500,000 samples and even a
modest set of common variants comes to hundreds of gigabytes for
common variants alone, and the per-row format assumes that reading
one row is cheap. Modern extensions handle the biobank regime by
reformatting: *GVCF* (genome VCF, #idx("GATK")GATK's emission format) records
calls at every position including reference; *sparse VCF* and
*project-VCF* formats elide the predominant `0/0` calls; Hail's
*MatrixTable* drops the row-oriented model entirely in favour of a
chunked sample-by-variant matrix on disk. Pick the layout that
matches your dominant query axis — per-variant or per-sample — not
the one the upstream tool defaults to.

#note[
  A VCF is a sparse differential encoding of a genome. Instead of
  transmitting 3 Gb of sequence per sample, you transmit only the
  approximately 4 million differences from a shared reference — a
  compression ratio around 750×. This is entropy coding at the
  file-format level: the reference is the codebook, the VCF is the
  per-sample delta. Indexable compressed VCFs (`bgzip` plus `tabix`)
  use the same design as indexed compressed time-series databases —
  block-level compression for random access by region without
  decompressing the whole file.
]


== Multi-Dimensional Data: HDF5, Zarr, Parquet <sec:hdf5-zarr>

The flat-file-with-index model breaks down for genuinely
multi-dimensional data. A single-cell experiment (Chapters 7 and 8)
produces a matrix of counts indexed by cell, gene, modality, and
batch. A #idx("methylation")methylation experiment produces per-CpG-per-sample matrices
at gigabase scale. A #idx("GWAS")GWAS produces summary-statistic tables across
millions of variants and dozens of traits. These are not
sparse-deltas-from-a-reference; they are dense tensors, and they
need *chunked, random-access, compressible multi-dimensional*
containers.

*HDF5* (the Hierarchical Data Format, version 5, from the HDF
Group, dominant since the late 1990s) is the long-standing answer.
A single HDF5 file is a self-describing binary tree of *groups*
(like directories) containing *datasets* (like files) with named
attributes. Datasets are stored in chunks, each chunk independently
compressed, each chunk independently addressable. The single-cell
community standardised on it via AnnData (`.h5ad`), with conventions
for `X` (the main expression matrix), `obs` (cell metadata),
`var` (gene metadata), and `obsm` / `varm` (reduced-dimensional
embeddings).

HDF5's weakness is that a single HDF5 file is, structurally, one
file. Parallel writes are possible but fiddly. Reading a slice from
the cloud requires either downloading the whole file or making
partial-range HTTP requests, which work but are clunky.

*Zarr* (Alistair Miles and others, 2018+) is the cloud-native
answer. Same chunked multi-dim array, but stored as *a directory of
small files*, one per chunk, each with a deterministic path. Reading
chunk `(7, 12, 3)` of a 3-D array is a single HTTP GET against
`s3://bucket/path/array/7.12.3`. Schema metadata lives in
`.zmetadata`, also a single file. Cloud object stores (S3, GCS,
Azure Blob) handle this layout natively. AnnData's Zarr backend has
become standard for new single-cell datasets above a million cells.

*Parquet* (from the Apache data-lake world; columnar; broadly used
outside bioinformatics) has entered the field via summary-statistic
tables. The columnar layout means that aggregations over a single
column (sum, filter, group-by) touch only the bytes of that column,
not the whole row. UK Biobank's GWAS summary stats ship as Parquet;
#idx("gnomAD")gnomAD's per-variant frequency tables ship as Parquet. Every
modern data-processing ecosystem — Spark, Dask, DuckDB, Polars,
Apache Arrow — reads Parquet natively, which is most of the reason
it has won outside genomics.

A practical scale anchor for the formats in this section:

- Single-sample WGS at 30× coverage: ~45 GB FASTQ.gz → ~80 GB BAM →
  ~25 GB CRAM (~14 GB with Q-binning).
- Cohort WGS at biobank scale: 100,000 samples × 25 GB CRAM ≈ 2.5 PB.
- Single-cell experiment: 1 M cells × 30,000 genes × 4 bytes ≈ 120 GB;
  multi-modal extensions multiply by three to five.
- GWAS summary statistics: 10 M variants × 10 fields × 8 bytes
  ≈ 800 MB per trait; thousands of traits run into multi-TB tables.

Assume from the outset that the data will not fit in RAM. Every
design choice in the field is shaped by this assumption.

#tip[
  Every format in this chapter is solving one engineering problem —
  "how do I get partial access without reading the whole thing?"
  FASTQ/BAM solve it with byte offsets and per-block indices.
  VCF/BCF solve it with tabix-style coordinate indices. HDF5 and
  Zarr solve it with chunk indices on an _n_-dimensional grid.
  Parquet solves it with column-major layout plus row-group indices.
  Same pattern every time: arrange the bytes so common queries land
  on contiguous ranges, and keep an index pointing to them.
]


== Where the Data Lives <sec:repositories>

Bioinformatics has unusual public-data norms relative to other
empirical fields. Major journals require raw sequencing data to be
deposited in one of a small set of repositories before publication.
The accession ID is the canonical citation; the URL is not. A paper
that reports new sequencing data and gives no accession is
unpublishable in most genomics journals.

The major open-access repositories form one tier:

- *SRA* (the Sequence Read Archive, NIH-hosted): the United States
  canonical home for raw sequencing data, free, with accession
  prefixes `SRP` (project), `SRR` (run), `SRS` (sample),
  `SRX` (experiment), and `PRJNA` (BioProject).
- *ENA* (the European #idx("nucleotide")Nucleotide Archive, EBI-hosted): the European
  counterpart; mirrors SRA daily. Prefixes `ERP`, `ERR`, `ERS`, `ERX`,
  `PRJEB`.
- *DDBJ Sequence Read Archive* (Japan): the third leg of the
  International Nucleotide Sequence Database Collaboration; mirrors
  SRA and ENA.
- *GEO* (the Gene Expression Omnibus, NIH-hosted) and its European
  counterpart *ArrayExpress* (EBI): processed-data repositories for
  microarray and #idx("RNA-seq")RNA-seq, typically count matrices and metadata.
  Most modern RNA-seq studies deposit raw reads in SRA and processed
  counts in GEO.

Controlled-access repositories form a second tier, gated by Data Use
Agreements and Institutional Review Board approvals:

- *dbGaP* (database of Genotypes and Phenotypes, NIH-hosted) holds
  individual-level US-cohort data behind an application process.
- *EGA* (the European Genome-phenome Archive, jointly EBI/CRG) is the
  European controlled-access counterpart. UK Biobank, FinnGen, and
  most European clinical cohorts deposit here.
- *GDC* (the Genomic Data Commons, NCI-hosted) holds cancer
  genomics: TCGA, TARGET, and most NCI-funded projects.

Every clinical cohort with identifiable genomic data sits in
dbGaP, EGA, or GDC — never in SRA or ENA. Knowing which tier a
dataset is in matters before you write a grant: open-access data is
free to fetch, controlled-access data takes weeks to gain access
and may require analysis inside a *trusted research environment*
where you cannot move bytes out.

#figure(
  image("../../diagrams/lecture-14/04-repositories.svg", width: 95%),
  caption: [
    The public-data landscape. Open-access on the left mirrors
    internationally; controlled-access on the right gates
    identifiable human data behind data-use committees. Every modern
    paper should cite one or more of these accession IDs.
  ],
) <fig:repositories>

=== Accession grammar

Every public dataset has a deposition record with three things: a
stable accession ID, structured metadata (sample tissue, sequencer
model, library protocol), and a list of associated publications.
The accession is the citable artefact — not the URL, not the file
hash, not the FTP path. A paper's "Data Availability" statement
points to accessions; tools like `pysradb` and `ffq` resolve
accessions back to download URLs.

The prefix grammar is worth memorising because it tells you which
repository to query and what kind of record to expect.
@fig:accession-grammar lays the prefixes out as a table and shows
the containment hierarchy that links them. A `BioProject` contains
one or more `Study` records; each `Study` contains multiple
`Sample` records; each `Sample` is sequenced one or more times in
distinct `Experiment` records; each `Experiment` produces one or
more `Run` records, which are the actual FASTQ files. The same
biological sample resequenced in two libraries is two experiments
under one sample; the same library run twice is two runs under one
experiment.

#figure(
  image("../figures/ch14/f2-accession-grammar.svg", width: 95%),
  caption: [
    Public-data accession grammar. Each repository uses a different
    prefix system, but the underlying containment hierarchy —
    BioProject → Study → Sample → Experiment → Run — is shared
    across the INSDC partners.
  ],
) <fig:accession-grammar>

=== Cloud-hosted copies

For large public datasets, direct download from FTP is slow and
awkward. Cloud-hosted copies have become the de facto standard. The
1000 Genomes Project is hosted on AWS (`s3://1000genomes`), GCS, and
Azure simultaneously. gnomAD lives in Google Cloud Storage at
`gs://gcp-public-data--gnomad`. GTEx and the All of Us research
cohort live on AnVIL and the Researcher Workbench, both Terra-based.
The Human #idx("pangenome")Pangenome Reference Consortium publishes its assemblies on
AWS and GitHub. UK Biobank ships its bulk genomic data via the RAP
(Research Analysis Platform, DNAnexus-hosted), where compute is
billed to the researcher.

Cost models split into two camps. *Free-egress* mirrors charge
nothing for downloads inside the same cloud region but charge for
cross-region transfer. *Requester-pays* buckets bill the person
running the query for byte egress, typically at around US\$0.01 per
gigabyte. The politics of cloud hosting have shifted who pays. Before
cloud mirrors, NIH-funded repositories absorbed storage and
bandwidth costs. Now, the agency pays for baseline storage and the
researcher pays per-analysis. For a 100-TB cohort analysis, the
egress bill alone runs into the thousands of dollars; for a
long-tail small-lab project, this is non-trivial.

The standard fetch tools are mostly mature. `sra-tools` (`prefetch`
plus `fasterq-dump`) is NIH's SRA client. `enaBrowserTool` is the
ENA-side equivalent. `aria2c` with 16 parallel HTTP connections
saturates most network links faster than the official clients;
`awscli` and `gcloud storage cp` handle the cloud mirrors. The
recurring pain points are not the tools but the *versioning*:
repositories rarely delete records, they re-version them, and a
download labelled "the data" can mean different bytes on different
days unless you pin the accession and submission timestamp.


== Reference Data Management <sec:references>

A reference genome is not a single file. "GRCh38" — the human
reference assembly released in 2013 — names a family of related
sequences that differ in which contigs they include and whether the
caller is expected to handle alternate haplotypes.

The major variants in active use:

- *Primary assembly only.* Chromosomes 1–22, X, Y, and the
  mitochondrion. ~3.0 GB FASTA. What most introductory tutorials
  reach for.
- *Primary + unplaced contigs.* The chromosomes plus roughly 300
  unassembled fragments. Adds about 100 MB. Some aligners use these
  to prevent reads from misplacing onto chromosomes.
- *Primary + ALT loci.* The chromosomes plus roughly 250 *alternate
  loci* — sequences for genomic regions with common large-scale
  variation, most famously the HLA region on chromosome 6. ALT-aware
  aligners distribute reads across these loci; ALT-unaware aligners
  treat them as separate chromosomes and fragment coverage.
- *Primary + unplaced + ALT + decoys.* Adds about 2,300 *decoy
  sequences* — artificial sequences that capture known reads not
  matching the main reference (centromere repeats, viral and
  bacterial contamination). Decoys give those reads a sink to flow
  into; without decoys, the reads scatter as misalignments onto
  real chromosomes.
- *Patch releases.* GRCh38 has been updated periodically;
  GRCh38.p14 (the fourteenth patch, 2022) is the current release,
  adding and fixing a modest number of sequences relative to the
  2013 original.

Different downstream tools expect different variants. BWA-MEM with
default settings wants the primary-only assembly for simplicity.
The ALT-aware mode (`bwa-mem -K`, plus a `.alt` file alongside the
FASTA) wants the primary plus ALT variant. DeepVariant's recommended
pipelines expect a specific variant —
`GRCh38_full_analysis_set_plus_decoy_hla.fa` — with a specific
decoy set. Mixing — running DeepVariant against a primary-only
reference, say — produces *garbage silently*: most reads still map,
but the ones that should have flowed into decoys instead land on
chromosomes as false low-quality alignments. There is no error
message. The downstream variant calls are wrong by a few percent in
ways that look like noise.

#figure(
  image("../../diagrams/lecture-14/05-which-grch38.svg", width: 90%),
  caption: [
    Which GRCh38? The same assembly name covers four or five
    compatibility tiers. The right answer depends on what your
    downstream tool expects — pick deliberately, pin by MD5.
  ],
) <fig:which-grch38>

#warn[
  GRCh38 reference mismatches are one of the most common silent
  bug sources in bioinformatics. Symptoms range from "a few percent
  of reads misplaced" (easy to miss) to "wrong variants called in
  HLA" (potentially clinical). Always record the *MD5 hash* of the
  reference used for alignment, not just its filename — filenames
  vary across tool distributions and over time. BAM and CRAM headers
  carry the reference MD5 in the `@SQ M5:` field; check it whenever
  you hand off data between pipeline steps.
]

=== Reference-management tools

Ad-hoc reference management is a minor nightmare. A typical
genomics team accumulates dozens of "GRCh38" files across users,
years, and projects, each subtly different. Four community tools
have grown up around the problem.

- *Ensembl FASTA releases* (EBI) ship versioned reference FASTAs
  with consistent naming and predictable URLs.
- *iGenomes* (Illumina) bundles references plus pre-built aligner
  indices for the major species. Convenient because "grab the
  tarball" is fast; heavyweight at about 100 GB per genome with all
  indices.
- *refgenie* (Stolarczyk and colleagues, 2020) is a reference-asset
  management server. It stores references and derived artefacts —
  BWA index, STAR index, salmon index, dbSNP VCF — under
  stable hash-based IDs. A CLI fetches a specific asset on demand;
  multi-user team servers are supported.
- *NCBI Datasets CLI* is NCBI's modern interface for reference
  downloads, replacing the older `efetch` toolkit for most cases.

Team-scale best practice converges on a small set of rules. Store
references in a shared path, not per-user. Pin by MD5, not by
filename. Co-locate every derived index (`.fai`, `.dict`, BWA, STAR,
salmon) next to the FASTA, so a downstream tool pointed at the
FASTA finds the index it needs. Record the reference used in every
analysis artefact — the BAM header, the VCF header, the pipeline's
provenance log. None of these are technical innovations; they are
laboratory discipline. The discipline matters because reference
mismatches do not crash.

=== GTF, GFF, BED — the annotation companions

A reference assembly is only sequence. Gene models — the locations
of transcripts, exons, coding regions — live in separate annotation
files alongside the FASTA.

*GFF3* (General Feature Format, version 3) is the modern standard.
Tab-delimited, nine columns, with a hierarchical `Parent=` attribute
that links exon records to their parent transcript and transcripts
to their parent gene. *GTF* (Gene Transfer Format) is the older
sibling — every line independent, no hierarchy — still in heavy use
because many older tools never moved on. *BED* (Browser Extensible
Data) is the minimalist tab-delimited format with three to twelve
columns; it is used for everything that does not need GFF's
hierarchy: ChIP peak calls, CpG-island annotations, callable-region
masks, blacklist files. Genome browsers (UCSC, IGV) consume BED
natively.

The pain point is that gene-model sources disagree. Ensembl,
RefSeq, and UCSC's `knownGene` track all annotate the human genome,
and they disagree about where some transcripts start and end. Cross-
referencing annotations from different sources requires careful
namespace conversion. Tools like `gffcompare` and `agat` help
reconcile them.


== Dependency Hell and the Container Era <sec:containers>

Every scientific-computing discipline has dependency problems.
Bioinformatics has them especially badly, and there are concrete
reasons why.

A typical bioinformatics pipeline mixes R (for statistics), Python
(for ML and glue), C and C++ tools (for aligners and variant
callers), and shell scripts (for orchestration). Each ecosystem
brings its own dependency manager, none of them designed to
cooperate. Many tools link against specific versions of system
libraries — `libhts`, `zlib`, `openssl`, `glibc` — and a tool that
worked on Ubuntu 20.04 may not work on Ubuntu 24.04 because
`libcrypto` was bumped. Academic software ages badly: a paper from
2014 with a public GitHub repo often has zero maintenance since.
The README says "run `./install.sh`"; the `install.sh` assumes
Python 2.7. Using it in 2026 is a multi-day project. Bioconductor
releases update every six months, and packages routinely break on
major releases. GPU-dependent tools — AlphaFold, Enformer, scVI,
deep-learning variant callers — demand specific CUDA + driver
combinations that mismatch local hardware silently.

Before containers became standard, "getting the damn thing to run"
frequently dominated actual analysis time. Stories of spending a
week on a single tool install were commonplace and not a sign of
weakness.

#figure(
  image("../../diagrams/lecture-14/06-dependency-hell.svg", width: 95%),
  caption: [
    Bioinformatics dependency hell. A "pipeline" is a stack of
    layers from shell glue down to glibc; a version bump anywhere
    can cascade through everything above it. Containers freeze the
    whole stack.
  ],
) <fig:dependency-hell>

=== Containers

A *container* is an OS-level package that bundles an application
plus all its userland dependencies — libraries, binaries, data —
sharing only the host operating system's kernel. It is not a virtual
machine: there is no hypervisor, no emulated hardware, no
guest-kernel overhead. Containers were enabled by Linux kernel
features developed in the late 2000s (`cgroups` and `namespaces`),
popularised by Docker in 2013, and adopted en masse by bioinformatics
around 2016–2017.

*Docker* is the industry standard. Images are built as layered
filesystems described by a `Dockerfile`; each instruction (`RUN
apt-get install …`, `COPY my-script /opt/`) adds a layer. Images
live in registries — Docker Hub, `quay.io`, the GitHub Container
Registry — and pulling and running an image is one command. The
catch for scientific computing is that Docker requires root-level
privileges on the host to manage its daemon; most academic clusters
forbid this for security reasons.

*Singularity* and its community fork *Apptainer* (after the
Singularity project went commercial in 2021) are the HPC and
academic answer. They run containers as the invoking user, read
Docker images natively (converting them to Singularity's `.sif`
format), and require no root daemon. `apptainer exec image.sif bwa
mem …` is the canonical academic-cluster invocation. Apptainer is
where most academic deployments are now; most published nf-core
pipelines list both a Docker and a Singularity invocation in their
documentation.

#figure(
  image("../../diagrams/lecture-14/07-container-vm.svg", width: 90%),
  caption: [
    Container versus VM. A container shares the kernel with the host
    and isolates everything else via namespaces; a VM emulates
    hardware. For packaging one tool per image, containers are
    vastly cheaper.
  ],
) <fig:container-vm>

#note[
  A container is *lightweight virtualisation*: instead of
  virtualising the hardware, it namespace-isolates processes,
  filesystems, and network. The kernel is shared; everything else is
  sandboxed. This is the software-engineering analogue of a
  system-on-chip IP block — treat "BWA with its exact `libhts`
  version plus its config files" as a black box with well-defined
  inputs and outputs, regardless of what else co-exists on the host.
  Nothing about the container is novel — it is a packaging and
  isolation pattern that pre-dates Docker by decades — but
  bioinformatics benefits disproportionately because per-tool
  complexity is unusually high.
]

=== Conda, mamba, and pixi

Before containers dominated, the userland answer was *conda*. Conda
installs packages into per-project environments (isolated
directories); packages can be prebuilt binaries (not just Python —
also C libraries, R, command-line tools); and conda solves a
SAT-like dependency-resolution problem to find a mutually consistent
version set. The problem is that conda's default solver is slow —
minutes for moderate environments — and sometimes fails to find a
solution at all. *Mamba* (a drop-in replacement using `libsolv`, Red
Hat's fast SAT solver) is now the community default. *Pixi*, newer,
unifies conda packages with project-based lockfiles (analogous to
JavaScript's `package-lock.json`); it is gaining adoption for
reproducible project setup.

Conda versus containers is a question of where the isolation
boundary sits. Conda is lighter — no root, no container runtime,
faster to iterate. Containers are more isolated — the whole OS
layer, not just userland — and reproducible across machines. Modern
practice uses conda locally for development and packages the final
pipeline as a container for deployment. *BioContainers* (the
community project that automatically converts bioconda packages to
Docker and Singularity images, publishing them on `quay.io`) covers
roughly 8,000 tools. A working `bioconda` recipe gets a
BioContainer image for free.

=== Pinning is not the same as containerising

Containers solve the *system-library* and *tool-version* problems
but they do not by themselves make an analysis reproducible. A
fully reproducible run requires pinning at six layers:

#figure(
  image("../figures/ch14/f3-reproducibility-stack.svg", width: 95%),
  caption: [
    The six layers of a reproducible bioinformatics analysis, with
    the canonical pinning artefact for each. A container covers
    layers 2–4; layers 1, 5, and 6 must be pinned independently —
    in the pipeline manifest, the reference manifest, and the
    script source.
  ],
) <fig:reproducibility-stack>

The pipeline code itself, tagged with a Git commit SHA or release
tag. The tool dependencies, locked by a `conda-lock.yml` or
equivalent. The language runtime version. The operating system and
container image, pinned by image *digest* (the immutable
`@sha256:…` content hash), not just a moving tag like `:latest`.
The reference data — every FASTA, dbSNP VCF, GTF, and trained model
— pinned by MD5. And finally the random seeds and hardware
non-determinism: explicit `numpy.random.seed` constants, fixed CUDA
devices, and acknowledgement that floating-point sums on different
GPU models may differ at the last few bits.

#warn[
  "I used a container, so it's reproducible" is not true by itself.
  A container pins layers 2 through 4 of the stack and nothing else.
  Reference data, random seeds, and pipeline version drift independent
  of the container. GPU-based tools such as AlphaFold can produce
  different outputs on different GPU models even with the same
  container, because mixed-precision arithmetic is non-deterministic
  across CUDA versions. A fully reproducible pipeline pins every
  layer explicitly, not just the tool layer.
]


== Workflow Languages <sec:workflows>

A real analysis pipeline has dozens of steps: FASTQ QC, adapter
trimming, alignment, sorting, duplicate-marking, base-quality
recalibration, variant calling, annotation, filtering, summary
generation. Each step uses a different tool, produces intermediate
files, has different CPU and memory requirements, and can be run
in parallel across samples. Shell scripts get you to step five
before they break.

The failure modes are predictable. If step ten fails, you re-run
everything from step one. Parallelising 100 concurrent samples
manually means coordinating the cluster scheduler by hand for every
sample. Asking for 200 GB of RAM per sample on a Slurm cluster
requires knowing the Slurm API. Reproducing the exact input-output
mapping six months later is impossible because the shell history
is gone. A script written for Slurm does not run on AWS Batch and
does not run on Google Cloud.

*Workflow managers* solve these by expressing the pipeline as a
*DAG* — a directed acyclic graph where each node is a step and each
edge is a data dependency. The manager handles scheduling, retry,
caching, parallel execution, and cloud-vs-cluster portability.

#figure(
  image("../../diagrams/lecture-14/08-workflow-dag.svg", width: 95%),
  caption: [
    A variant-calling pipeline as a DAG. Nodes are tools with
    resource specifications; edges are data dependencies. The
    workflow manager turns this graph into a scheduled execution
    plan across the available compute.
  ],
) <fig:workflow-dag>

=== Four workflow languages

Four workflow languages dominate, in roughly the order of their
publication.

*Snakemake* (Köster and Rahmann, 2012) is the oldest. It models the
pipeline as a set of Make-style rules with input and output files;
the DSL is Python, which makes ad-hoc pipeline development quick.
Adoption is academic-heavy and excellent for rapid prototyping.
Cloud support has historically been less mature than Nextflow's but
is catching up via executors for Slurm, Kubernetes, and the major
cloud platforms.

*Nextflow* (Paolo Di Tommaso, Cedric Notredame, and colleagues,
2017) is dataflow-oriented: pipelines are written as *channels* and
*processes*, with channels carrying typed streams of files or
metadata between processes. The DSL is built on Groovy. Nextflow
has strong native support for AWS Batch, Azure, Google Cloud Life
Sciences, and the major schedulers. Resume-on-failure is invoked by
a single `-resume` flag. Industry adoption is heavy.

*WDL* (the Workflow Description Language, Broad Institute, 2015)
was designed for the Cromwell execution engine. It uses a task and
workflow structure with static typing. WDL dominates in the Broad
and Cromwell ecosystems — Terra, the Genomic Data Commons, the
All of Us workbench — and is verbose but clear.

*CWL* (the Common Workflow Language, multi-organisation standard
from 2016 onward) standardises workflow descriptions across multiple
execution engines. Adoption is concentrated in European
bioinformatics infrastructure (ELIXIR) and less common in
GitHub-native projects.

A practical comparison:

#figure(
  table(
    columns: (auto, 1fr, 1fr, 1fr, 1fr),
    align: (left, center, center, center, center),
    inset: 6pt,
    stroke: 0.5pt + rgb("#c4c1b6"),
    table.header(
      [], [*Nextflow*], [*Snakemake*], [*WDL*], [*CWL*],
    ),
    [Paradigm], [Dataflow], [Make-rules], [Task-based], [Task-based],
    [DSL base], [Groovy], [Python], [Custom], [YAML/JSON],
    [Cloud-native], [Excellent], [Growing], [Cromwell-centric], [Moderate],
    [Community], [nf-core], [Academic], [Broad / Terra], [ELIXIR],
    [Learning curve], [Moderate], [Easy], [Easy-moderate], [Steep],
  ),
  caption: [
    The four major workflow languages, compared on the axes that
    drive adoption decisions. Pick by what the rest of the team
    uses; if no team, pick by deployment target.
  ],
) <fig:workflow-table>

#note[
  A workflow manager is a specialised *dataflow runtime*. Processes
  are nodes producing outputs as a function of inputs; the runtime
  schedules execution by topological sort of the dependency DAG,
  handles retries, caches intermediates for resume, and parallelises
  independent branches. This is the same abstraction as Apache
  Spark's DAG scheduler, TensorFlow's computation graph, or the
  Unix `make` utility scaled out to a cluster. Bioinformatics
  re-invented the idea for domain-specific reasons — data size,
  container packaging, cloud deployment — but the core is a
  textbook dataflow runtime.
]

=== Resume, retry, resource declaration

The mechanical features a workflow manager gives you, in
descending order of how often they save a researcher's day:

*Resume on failure.* Every process output is cached by a hash of
its input files, the tool version, and the command line. If step
ten fails and you fix the bug, rerunning skips everything up
through step nine because the cache hits. Nextflow's `.nextflow/cache`
and Snakemake's output-timestamp tracking implement this in slightly
different ways. In practice: a six-hour pipeline that fails at hour
five restarts at hour five plus one minute, not at hour zero.

*Per-process retry.* If a tool fails spuriously — out-of-memory
killed, transient network error, race condition on a shared
filesystem — the workflow retries with backoff. Typical configuration:
three retries with exponential backoff, optionally doubling memory on
each retry.

*Per-process resources.* CPU, memory, walltime, queue, and container
image declared per-process in the pipeline definition, respected by
whichever scheduler the run targets. A typical Nextflow process
declaration:

```
process bwa_mem {
    cpus 16
    memory '32 GB'
    time '4 h'
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'
    ...
}
```

*Execution backends.* Slurm, PBS, SGE, AWS Batch, Azure, Google
Cloud Life Sciences, or local execution — specified in a separate
configuration file, not in the pipeline definition. The same
pipeline runs on any backend by switching the config.

*Provenance reports.* Every run produces an HTML or JSON report
documenting which steps ran, their inputs and outputs, their
runtimes, and their resource usage. Nextflow Tower (commercial
Seqera Platform) productionises this for teams; Snakemake's
`--report` flag produces the equivalent for academic use.

#figure(
  image("../../diagrams/lecture-14/09-resume-logic.svg", width: 95%),
  caption: [
    Resume-on-failure. Every process output is cached by a hash of
    its inputs and command line; a retried run skips cached steps
    and restarts at the first miss. The same cache structure makes
    iterative pipeline development tolerable.
  ],
) <fig:resume-logic>

=== nf-core and the curated-pipeline community

*nf-core* is a curated collection of Nextflow pipelines maintained
by a community of about a thousand contributors (Phil Ewels and
colleagues, 2020). Every nf-core pipeline follows a common style
guide: standardised configuration layout, mandatory testing and
linting, complete documentation, BioContainer-wrapped tools, and
CI/CD that runs the pipeline on test data for every pull request.
Releases are tagged in GitHub and archived to Zenodo with a DOI;
container images are published to `quay.io`.

The flagship pipelines cover most of what a working genomics lab
needs: `nf-core/rnaseq` for differential expression, `nf-core/sarek`
for germline and somatic variant calling, `nf-core/atacseq` and
`nf-core/chipseq` for chromatin assays, `nf-core/scrnaseq` for
single-cell, `nf-core/proteomics` for mass spectrometry,
`nf-core/methylseq` for bisulfite sequencing, `nf-core/fetchngs`
for downloading public data. The list grows by a few pipelines a
year.

Not every lab uses nf-core, but most modern production
bioinformatics runs on something nf-core-shaped: containerised,
DAG-expressed, parameterised by a single config file, archived with
a DOI. If your lab has a "which pipeline for X" problem, check
nf-core first; you will probably find it.

#figure(
  image("../../diagrams/lecture-14/12-nf-core-flow.svg", width: 95%),
  caption: [
    nf-core as a pipeline-curation factory. Every release is
    CI-tested, containerised, and DOI-archived. Community-curated
    is not the same as vendor-curated; nf-core pipelines stay alive
    because dozens of users care.
  ],
) <fig:nf-core-flow>


== Benchmarking and the GIAB Truth Set <sec:benchmarks>

A method that has never been benchmarked is a method nobody trusts.
Bioinformatics has internalised this and built a culture of
community-curated benchmarks. The flagship is the Genome in a
Bottle Consortium.

*GIAB* (NIST-led, Justin Zook and colleagues, 2014 onward) produces
highly validated *truth sets* for human variant calling. The
canonical GIAB sample is *HG002*, an Ashkenazi-Jewish trio member,
sequenced with every major sequencing technology (Illumina short-
read, PacBio HiFi, Oxford Nanopore, Strand-seq, BioNano optical
maps), variant-called by every major caller, and manually curated.
The resulting high-confidence variant VCF plus its accompanying
*confident-region BED* (regions where the truth is reliable) is the
field's gold standard.

Successive releases have widened the confident region and added
more samples:

- GIAB v3.2.2 (2016): the early version cited in most
  variant-caller papers from 2016 through 2019.
- GIAB v4.2.1 (2021): expanded samples (HG001 through HG007),
  improved structural-variant coverage.
- GIAB v5.0 (2024+): T2T-based; covers pericentromeric and
  repetitive regions that were previously excluded as unreliable.

Each release ships three artefacts: a truth VCF (the set of real
variants), a confident-region BED (the regions where the truth set
is reliable), and per-region callability metrics. Outside the
confident region, *the ground truth is undefined* — the regions are
too repetitive for any current technology to fix the truth, and
reporting precision or recall there is meaningless.

#figure(
  image("../../diagrams/lecture-14/10-giab-schema.svg", width: 95%),
  caption: [
    The GIAB truth-set schema. Confident regions (bands) plus truth
    variants (ticks) plus callability metadata. Outside confident
    regions, "accuracy" is not a defined quantity.
  ],
) <fig:giab-schema>

=== hap.py and the arithmetic of precision and recall

*hap.py* (Peter Krusche and colleagues at Illumina, 2019; the
Google fork has been the active version since 2021) is the standard
variant-caller benchmarker. Given a test VCF and the GIAB truth
plus confident region, it does five things: normalises both VCFs
(left-aligns indels, splits multi-allelics, deduplicates);
restricts to the confident region; matches test variants to truth
variants — non-trivial because the same deletion can be written as
`REF=AC ALT=A` or `REF=CA ALT=C`; classifies each variant as a true
positive, false positive, or false negative; and reports precision,
recall, and F1, with per-variant-type breakdowns.

The arithmetic itself is textbook detection theory. Given true
positives (TP), false positives (FP), and false negatives (FN):

$ "precision" = "TP" / ("TP" + "FP") $
$ "recall" = "TP" / ("TP" + "FN") $
$ "F1" = (2 dot "precision" dot "recall") / ("precision" + "recall") $

Precision is the fraction of called variants that are real. Recall
is the fraction of real variants that were called. F1 is the
harmonic mean — pessimistic about whichever of the two is smaller.
A modern caller on HG002 against #idx("GIAB")GIAB v4.2.1 looks like
@fig:happy-matrix: about 3.45 million true-positive SNV calls,
about 2,800 false positives, about 14,500 false negatives, for an
SNV F1 around 0.9975.

#figure(
  image("../figures/ch14/f4-happy-confusion-matrix.svg", width: 95%),
  caption: [
    The #idx("hap.py")hap.py confusion matrix for a representative modern caller
    against GIAB HG002 v4.2.1. Single-number F1 hides the
    per-variant-class variance — SNV F1 sits above 0.99, but #idx("indel")indel
    F1 above 16 bp collapses to around 0.8.
  ],
) <fig:happy-matrix>

The per-variant-type breakdown is the reason hap.py exists. A
caller can be excellent on SNVs and middling on indels; reporting
only the aggregate F1 hides this. hap.py's output table has
separate rows for SNV, INDEL of 1–5 bp, INDEL of 6–15 bp, and
INDEL of 16 bp and longer, each with its own precision and recall.
Modern callers cluster above F1 = 0.999 on SNVs and somewhere
between 0.8 and 0.99 on indels depending on size, with long-read
methods winning the long-indel regime decisively.

#figure(
  image("../../diagrams/lecture-14/11-precision-recall.svg", width: 95%),
  caption: [
    Precision and recall on GIAB HG002 v4.2.1 for several callers.
    Modern short-read callers cluster above F1 = 0.999 for SNVs;
    indel performance is materially worse and spreads more across
    methods. Per-variant-type reporting is not optional.
  ],
) <fig:precision-recall>

#note[
  GIAB is bioinformatics's ImageNet. A community-curated,
  field-standard benchmark that defines a precise evaluation
  protocol, is maintained by a neutral party (NIST), has
  per-category breakdowns so methods can be compared where they
  actually differ, and is expected of every new method paper. The
  same function as ImageNet for image classification, GLUE for
  natural-language processing, and #idx("CASP")CASP for #idx("protein structure")protein structure (which
  the next chapter covers): give the field a way to compare methods
  objectively and drive iterative improvement. Without benchmarks
  of this kind, method rankings are anecdotal and progress stalls.
]

=== What benchmarks do and do not guarantee

A benchmark measures performance on the specific data it evaluates.
A handful of caveats matter when reading "F1 = 0.997 on GIAB HG002"
in a paper.

*Overfitting.* Published methods can be implicitly tuned on GIAB
HG002 — the sample is sequenced, called, and re-called more than any
other in the world. A 0.001 F1 improvement on the benchmark may not
generalise to other samples.

*Coverage.* GIAB covers about 95% of the autosomes for SNVs but
much less for SVs, repeat regions, and HLA. Method rankings outside
the confident region are undefined.

*Population.* HG002 is Ashkenazi-Jewish. GIAB has smaller trio sets
for African and East Asian samples, but their coverage is less
mature. A method that ranks well on GIAB does not necessarily rank
the same way in samples from under-represented populations.

*Technology.* GIAB grew up around Illumina short reads; long-read
benchmarks are newer and less mature. A 2024 reader should interpret
"F1 = 0.997 on GIAB HG002" as "this method is competitive on
standard human WGS in European-ancestry samples." Deployment in
clinical tumour samples, non-human species, or rare conditions
requires additional validation.

Every major sub-field has a benchmark of similar shape. SVs use the
GIAB #idx("SV")SV benchmark plus *Truvari*. RNA-seq quantification has *SQANTI*
and the SRA test sets. Single-cell integration has the *sccloud*
benchmark and the cellxgene "nice tissues" benchmarks. Genome
assembly has *QUAST*. #idx("metagenomics")Metagenomics has *CAMI*. Protein structure
has *CASP*. The common pattern: a well-characterised ground truth,
a reference evaluation script, and a community norm that "if you
publish a new tool, you benchmark on the standard set."


== The Culture: Preprints, Open Source, Archival <sec:culture>

The technical layers of this chapter — formats, repositories,
references, containers, workflows, benchmarks — sit on top of a
cultural layer that is unusual among empirical fields and worth
naming explicitly.

*Preprints first.* Bioinformatics has embraced preprints harder
than most biomedical fields. *bioRxiv* (launched 2013) is where
most methods papers land before peer review, and the community
reads them routinely. *medRxiv* (launched 2019) is the clinical
and epidemiology counterpart. ML-heavy work — #idx("AlphaFold")AlphaFold,
#idx("Enformer")Enformer, #idx("scVI")scVI — also appears on arXiv, often cross-listed.
The typical flow is preprint, community feedback (in GitHub
issues, on social media, in citation comments), revision, journal
submission, peer review, publication. The time between preprint
and formal publication is 6–18 months, and citations begin
accumulating on the preprint version. The practical consequence:
your job as a reader is to evaluate preprints on technical merit
and not wait for peer review to consider them.

*Open-source by default.* Almost all new tools are GitHub-hosted
with permissive licenses (MIT, BSD, Apache). Tools that refuse to
share code are viewed with suspicion. Raw data must be deposited
in a public repository before publication (Part 2). Benchmarks,
pipelines, and documentation are open and community-curated.
Bioconductor packages must ship a vignette; modern Python tools
auto-generate Sphinx or MkDocs documentation. The combination is
unusual: a field where the expected output of new work is a paper
*plus* working code *plus* deposited data *plus* a documented
pipeline. For an EE student entering the field, this is a
material expectation — your published work is expected to come
with working code, and "my tool is proprietary, I share by request"
is read as a weak signal.

*Institutional anchors.* The field runs on a surprisingly small set
of core institutions. NIH (specifically NHGRI and NCBI) funds and
hosts SRA, GenBank, and dbGaP. EBI (the European Bioinformatics
Institute at Wellcome Sanger) operates ENA, Ensembl, and a
significant fraction of European bioinformatics infrastructure. The
Broad Institute publishes GATK and much of the Cromwell/WDL
ecosystem. The Sanger Institute publishes #idx("nf-core")nf-core, #idx("Nextflow")Nextflow, and
the pangenome consortium tooling. UCSC anchors the Genome Browser
and historically much of reference-assembly curation. Commercial
players matter less than in most fields; Illumina dominates
sequencing hardware, but the analysis ecosystem is community-driven.

*The reproducibility crisis.* Bioinformatics has its own version of
the broader reproducibility crisis. Old papers frequently cannot be
rerun because data access has changed (dbGaP permissions expired,
SRA records re-versioned), tool versions were never pinned,
reference builds were mismatched, or the original authors left for
industry and their GitHub repos went stale. Surveys suggest that
30% or more of published bioinformatics pipelines cannot be
re-executed five years later. The modern norms covered in this
chapter — containerised pipelines, MD5-pinned references, code
archived in Zenodo with a DOI, container images archived at a
stable registry — are the field's response. Whether they hold up at
10-year timescales is an open question.

#warn[
  A pipeline that "worked last year" may not work this year even
  with no user changes. Upstream packages update. CDN URLs change.
  GitHub repositories are renamed or deleted. Container registries
  rotate their retention policies. The only durable archival is the
  *exact container image* deposited in a stable registry with an
  immutable identifier — Zenodo's DOI system, the Software Heritage
  Archive, S3 Glacier. "My code is on GitHub" is necessary but not
  sufficient.
]


== Summary <sec:summary>

- *Formats are partial-access engineering.* Every good bioinformatics
  format — BAM, BCF, HDF5, Zarr, Parquet — lays out data so common
  queries touch contiguous bytes, with an index for seek. CRAM goes
  further with content-aware compression, storing reads as positions
  plus diffs against a reference.
- *Reference management is its own discipline.* "#idx("GRCh38")GRCh38" names
  multiple references; know which variant your tool expects; pin by
  MD5, not by filename. Reference mismatches are the silent-failure
  king and produce wrong calls with no error message.
- *Public data lives in tiered repositories.* SRA, ENA, GEO for
  open access; dbGaP, EGA, GDC for clinical and identifiable data.
  Cloud-hosted copies are the working standard at scale;
  requester-pays buckets shift cost from the agency to the analyst.
- *Dependency hell is worse in bioinformatics than in most fields*
  because of mixed-language ecosystems plus ageing academic
  software plus deep system-level dependencies. Containers (#idx("Docker")Docker,
  Apptainer) plus environment managers (conda, mamba, pixi) plus
  BioContainers solve the packaging problem; pinning the rest is on
  you.
- *Workflow languages are dataflow runtimes.* Nextflow, #idx("Snakemake")Snakemake,
  WDL, and CWL all express the pipeline as a DAG of tools with
  resource specs. Resume-on-failure, parallelism, and cloud
  execution are first-class. Nextflow + nf-core is the
  production-cloud standard; Snakemake is the research-cluster
  standard.
- *Benchmarks anchor the field.* GIAB plus hap.py is the variant-
  calling standard. Per-variant-type breakdowns matter; aggregate
  F1 hides large per-class variance.
- *Reproducibility requires pinning at six layers.* Pipeline code,
  tool versions, language runtime, container image digest, reference
  data, and random seeds. Containers cover only the middle three.


== Exercises <sec:exercises>

#strong[1.] #emph[CRAM arithmetic.] Given a 30× human BAM at 80 GB,
estimate the on-disk size of (a) lossless CRAM, (b) CRAM with 4-bin
quality binning, (c) lossless CRAM if the same data were 60×
coverage instead of 30×. State which streams in the CRAM are
linear in coverage and which are constant.

#strong[2.] #emph[VCF compression.] A VCF for a single sample with
4.8 million variants takes 1.5 GB compressed. Estimate the size of
the same file in BCF, and the size of a joint-called VCF for 1,000
samples assuming variant-site overlap follows the typical
common/rare distribution (about 50% of sites in any single sample
are private to that sample, 50% are shared). Sketch the arithmetic.

#strong[3.] #emph[Accession parsing.] You are handed the
accessions `PRJEB12345`, `SRR9876543`, `GSE54321`, `phs000007`,
and `EGAS00001000456`. For each one, identify (a) the repository
that hosts it, (b) the kind of record it points to, (c) whether
the data is open-access or controlled-access. Two sentences each.

#strong[4.] #emph[Reference choice.] You have been given a
#idx("DeepVariant")DeepVariant pipeline that expects the recommended
`GRCh38_full_analysis_set_plus_decoy_hla.fa`. A collaborator hands
you BAMs aligned against `hg38.fa` from UCSC. Describe what will go
wrong, how you would detect it from the BAM header alone, and what
the fix is. (Hint: the `@SQ` `M5:` field.)

#strong[5.] #emph[Container vs conda.] Sketch two reproducibility
recipes for a pipeline that runs `bwa-mem`, `samtools`, and a
Python script that uses `pysam`: one with a single Docker image,
one with a conda environment. Identify three failure modes that
the Docker recipe handles and the conda recipe does not, and one
failure mode that affects both equally.

#strong[6.] #emph[Workflow DAG sizing.] Sketch the DAG for an
`nf-core/sarek`-style germline variant-calling pipeline running 100
samples in parallel. Identify (a) which steps fan out per sample,
(b) which steps fan in across samples (joint calling), (c) which
steps are per-chromosome, (d) the longest critical path.

#strong[7.] #emph[hap.py arithmetic.] A caller produces 3,452,200
true positives, 2,800 false positives, and 14,500 false negatives
on SNVs. Compute precision, recall, and F1 to four significant
figures. Then assume an aggressive filter doubles the FP count to
5,600 and halves the FN count to 7,250 — what does this do to F1?
Express the change in absolute and relative terms.

#strong[8.] #emph[(Open-ended.)] Pick one published bioinformatics
tool from the last five years that includes a docker image, a
workflow-language pipeline (Snakemake or Nextflow), and a GIAB
benchmark in its release. Read its methods paper and the README.
In one paragraph, evaluate whether the release would still run on
your local machine in 2030: which of the six reproducibility layers
in @fig:reproducibility-stack is the weakest?


== Further Reading <sec:further-reading>

- *Li, H., et al.* (2009). "The Sequence Alignment / Map Format and
  SAMtools." _Bioinformatics_ 25: 2078–2079. SAM and BAM in one
  short paper; still the spec the rest of the field is built on.
- *Hsi-Yang Fritz, M., Leinonen, R., Cochrane, G., Birney, E.*
  (2011). "Efficient Storage of High Throughput #idx("DNA")DNA Sequencing
  Data using Reference-Based Compression." _Genome Research_ 21:
  734–740. The CRAM paper.
- *Danecek, P., et al.* (2011). "The Variant Call Format and
  VCFtools." _Bioinformatics_ 27: 2156–2158. The VCF specification
  paper; pairs with the more recent SAMtools / BCFtools review.
- *Köster, J., and Rahmann, S.* (2012). "Snakemake — A Scalable
  Bioinformatics Workflow Engine." _Bioinformatics_ 28: 2520–2522.
- *Di Tommaso, P., Chatzou, M., Floden, E. W., et al.* (2017).
  "Nextflow Enables Reproducible Computational Workflows."
  _Nature Biotechnology_ 35: 316–319.
- *Ewels, P. A., Peltzer, A., Fillinger, S., et al.* (2020). "The
  nf-core Framework for Community-Curated Bioinformatics Pipelines."
  _Nature Biotechnology_ 38: 276–278.
- *Zook, J. M., Chapman, B., Wang, J., et al.* (2014). "Integrating
  Human Sequence Data Sets Provides a Resource of Benchmark #idx("SNP")SNP
  and Indel Genotype Calls." _Nature Biotechnology_ 32: 246–251. GIAB.
- *Krusche, P., Trigg, L., Boutros, P. C., et al.* (2019). "Best
  Practices for Benchmarking Germline Small-Variant Calls in
  Human Genomes." _Nature Biotechnology_ 37: 555–560. The hap.py
  paper and the field's small-variant-benchmark protocol.
- *Stolarczyk, M., Reuter, V. P., Smith, J. P., et al.* (2020).
  "Refgenie: A Reference Genome Resource Manager." _GigaScience_ 9:
  giz149.
- *Kurtzer, G. M., Sochat, V., and Bauer, M. W.* (2017).
  "#idx("Singularity")Singularity: Scientific Containers for Mobility of Compute."
  _PLoS ONE_ 12: e0177459.
- *da Veiga Leprevost, F., Grüning, B. A., et al.* (2017).
  "BioContainers: An Open-Source and Community-Driven Framework
  for Software Standardization." _Bioinformatics_ 33: 2580–2582.
- *GA4GH file-format specifications.* `samtools.github.io/hts-specs/`
  The hub for SAM, BAM, CRAM, VCF, and BCF spec documents. Read
  the spec before writing any code that parses these formats.
