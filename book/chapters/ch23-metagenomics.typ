#import "../theme/book-theme.typ": *

= Metagenomics and the Microbiome: Community Signal Demixing <ch:metagenomics>

#matters[
  The microbiome is the first biological object in this book that is not
  a single organism. It is a community — hundreds of bacterial species,
  a few dozen archaea, fungi, viruses, and protozoans, all in one tube
  — and almost everything you learned from Chapters 2 through 22 has to
  be rebuilt for the case where the sequenced sample contains many
  genomes in unknown proportion. Read alignment becomes read
  classification. Variant calling becomes strain demixing. Differential
  expression becomes differential abundance under a sum-to-one
  constraint. Diversity replaces variance as the central summary
  statistic. The community is treated as a population and its sampling
  becomes the assay. Get the constraints right and the field is
  enormously rewarding; the gut microbiome predicts immunotherapy
  response, the oral microbiome anticipates pancreatic cancer, and the
  hospital sewage stream tracks the spread of antibiotic resistance
  genes between species in real time. Get them wrong and your
  significant taxa are artefacts of compositional bias, your replication
  rate is twenty per cent, and the field finds another reason to
  distrust its own literature. This chapter is the discipline for
  staying on the right side of that line.
]

A microbiome sample is a tube of DNA from hundreds of organisms in
proportions you do not know. The sequencer reads from that mixture
without preference, and the analysis problem is to recover who is
present, in what proportion, doing what, and how that maps to the
host. None of those questions admit the clean
one-genome-one-FASTQ-one-alignment framing of the preceding chapters.
A read might match six bacterial species at the genus level and one
specific strain at the species level; the same species might be a
commensal in one person and a pathogen in another; the abundance of
species A you report is a function of species B's abundance because
the total number of reads is constant. The biology is more interesting
than that of a single human genome and the statistics are harder.

The chapter is organised around the engineer's question: given a
microbiome sample, what sequence of decisions turns it into a defended
claim? The first section walks the scale of the microbial world and
names the two main assays — amplicon sequencing of the 16S rRNA
ribosomal-RNA gene, and shotgun sequencing of total DNA. The second
covers the modern 16S pipeline, from DADA2 amplicon-sequence-variant
inference to QIIME 2 visualisation and taxonomic assignment against
SILVA or the Genome Taxonomy Database. The third walks the shotgun
side: Kraken2 k-mer classification, MetaPhlAn marker-gene profiling,
HUMAnN functional profiling, and metagenome-assembled-genome recovery.
The fourth introduces diversity as the field's central summary
statistic — Shannon entropy, Simpson's index, Faith's phylogenetic
diversity, Bray-Curtis dissimilarity, UniFrac — and the fifth lays
down the compositional-data trap that every microbiome paper from
2008 onward had to learn to escape. The sixth surveys host-microbiome
connections, with the immunotherapy story as the most clinically
load-bearing example. A closing section walks the engineering pitfalls
and the 2024 frontier.

The chapter is the longest single departure from the
one-genome-one-FASTQ paradigm in this book. Every primitive you have
to relearn pays off as a primitive you keep. Compositional thinking
shows up again in single-cell deconvolution; demixing shows up again
in cell-free DNA and microbial cell-free fragments in plasma; the
diversity-as-divergence framing appears in immune-repertoire
sequencing and clonal-population analysis. The microbiome is the
domain where these moves are most exposed; learning them here makes
the rest of the modern bioinformatics literature more readable.


== Why Microbes <sec:ch23-why-microbes>

The human body is a co-organism. Roughly $10^13$ microbial cells live
on and in a typical adult — about as many as human cells, depending on
which census you trust and how you count, with the older one-to-ten
ratio (Savage, 1977) revised down by Sender, Fuchs, and Milo
(*PLoS Biology* 2016) to almost-exactly one-to-one. The microbial
gene count is more lopsided: about three million unique genes across
the human microbiome (Qin et al., *Nature* 2010), versus the roughly
twenty thousand protein-coding genes in the human genome. Most of
those microbial genes encode functions our own cells cannot do —
digesting plant cell walls, synthesising vitamins B12 and K, modifying
drug metabolites, training the immune system, occluding pathogens
at mucosal surfaces.

Ninety-nine per cent of the microbes are bacteria. The remainder is a
diverse minority — archaea, fungi, viruses, protozoans — that
nonetheless matters disproportionately in some niches. *Candida albicans*
overgrowth is a problem in immunocompromised patients;
bacteriophages outnumber bacteria in many sites and shape bacterial
community structure through predation; archaeal *Methanobrevibacter*
species are common gut residents whose hydrogen-consuming metabolism
matters for energy harvest. The site matters too. Gut microbiomes are
the densest and most studied (about $10^14$ cells per gram of faecal
material); oral microbiomes are the second densest and connect to
cardiovascular and pancreatic disease; skin microbiomes are
niche-specific (sebaceous, moist, dry) with distinct *Cutibacterium*,
*Staphylococcus*, and *Corynebacterium* populations; vaginal
microbiomes are usually dominated by one or two *Lactobacillus* species
and matter for preterm birth risk; respiratory microbiomes are small,
hard to sample without host-DNA contamination, and matter in asthma
and chronic obstructive pulmonary disease.

#note[
  The "we are more bacteria than human" line you will encounter in
  popular writing originated with Savage's 1977 estimate of a ten-to-one
  cell ratio and survived for decades. Sender et al.'s 2016 reanalysis
  revised the ratio down to approximately one-to-one, with substantial
  inter-individual variation. The headline number changed; the
  biological point — that the microbiome is a quantitatively serious
  part of any human's cell budget — did not.
]

Microbiome research aims at three questions. The first is taxonomic:
*who is there?* What species are present, and at what abundance?
The second is functional: *what can they do?* Which metabolic
pathways, virulence factors, antibiotic resistance genes are encoded?
The third is host-facing: *what does that mean for the human?* Does
the community predict disease, response to a drug, response to a
dietary intervention? Each question maps to a different sequencing
strategy and a different analysis pipeline.

The strategy that answers question one is *16S ribosomal-RNA amplicon
sequencing*. The 16S gene is a roughly 1500-base-pair bacterial
ribosomal RNA gene with the unusual property of alternating between
strongly conserved regions (where ribosome function constrains every
base) and hypervariable regions (where the constraint is relaxed and
species-level differences accumulate). The conserved regions anchor
PCR primers; the variable regions between them carry the taxonomic
signal. A typical 16S protocol amplifies the V3 — V4 segment (about
460 bp) using universal primers, sequences the product on Illumina at
about ${\$}30$ per sample, and assigns taxonomy by matching the
amplicon to a database of known 16S sequences.

The strategy that answers questions two and three is *shotgun
metagenomics*. Extract total DNA from the sample, fragment it,
sequence randomly, and accept whatever reads come out — including
reads from archaea, fungi, viruses, protozoans, and the host. A
typical shotgun protocol generates about 5 gigabases per sample at a
cost around ${\$}300$. Classification, function, and assembly are all
downstream-computational rather than upstream-experimental.

#figure(
  image("../../diagrams/lecture-23/01-16s-vs-shotgun.svg", width: 95%),
  caption: [
    16S amplicon versus shotgun metagenomics. The choice is constraint-driven:
    cost, biomass, kingdom coverage, and whether functional content
    matters together pick the assay.
  ],
) <fig:assay-compare>

The choice between them is constraint-driven. @fig:design-decision
walks four constraints — budget, biomass, kingdom coverage, functional
content — and converges on one of three platforms.

#figure(
  image("../figures/ch23/f1-design-decision.svg", width: 100%),
  caption: [
    A four-question decision tree for picking the assay. 16S is cheap,
    bacteria-only, and genus-resolution; shotgun is expensive and
    multi-kingdom with strain-resolution; full-length 16S on PacBio
    HiFi sits between them.
  ],
) <fig:design-decision>

A common hybrid design uses 16S on every sample in a large cohort —
hundreds of patients, cheap per-sample — to characterise broad
patterns, then runs shotgun on a stratified subset (typically ten to
twenty per cent) selected to follow up on a hypothesis. The hybrid
keeps cohort-level statistical power and gets the functional content
where it matters most.

#note[
  The microbiome maps cleanly onto population-genetics intuition. The
  relative abundance of a taxon plays the role of an allele frequency.
  The community's richness plays the role of an effective population
  size. Stochastic abundance fluctuations between samples are drift;
  selection by host diet, antibiotic exposure, or immune state is
  selection; faecal microbiota transplant or skin contact is
  migration. Population-genetics metrics like Wright's $F_("st")$ and
  Nei's distance have direct microbial-community analogues in UniFrac
  and Bray-Curtis dissimilarity. Many of the tools in @sec:ch23-diversity
  generalise the same population-comparison statistics that Chapter 13
  used for human ancestry.
]

The deep engineering frame is that a metagenomic sample is a
*mixture of DNA from many genomes*, and recovering the per-taxon
abundances is a *demixing problem*. The reference-genome database is
a dictionary; the reads are an observation; the abundances are the
mixing coefficients. Kraken2 solves the demixing approximately by
exact k-mer matching and least-common-ancestor aggregation;
MetaPhlAn solves it with a curated set of clade-specific marker
genes; Salmon-meta and Bracken use expectation-maximisation to
distribute read mass across reference genomes when no k-mer hit is
unique. Every method in this chapter is a point on a precision —
recall — speed surface for that same demixing problem.


== The 16S Amplicon Workflow <sec:ch23-amplicon>

The bacterial 16S ribosomal-RNA gene is approximately 1500 base pairs
long, encodes a structural RNA component of the 30S ribosomal subunit,
and is present in every bacterial genome (often in multiple copies).
Carl Woese realised in the late 1970s that its alternating
conservation pattern made it an ideal phylogenetic marker: regions
under strong functional constraint stay essentially identical across
all bacteria, while loops and bulges between them tolerate enough
variation to encode hundreds of millions of years of divergence.
Woese's 1977 reorganisation of the tree of life — splitting prokaryotes
into Bacteria and Archaea on the basis of 16S sequence — was a direct
consequence of taking the gene seriously as a clock. Modern 16S
amplicon sequencing is the same idea industrialised.

#figure(
  image("../../diagrams/lecture-23/07-16s-structure.svg", width: 95%),
  caption: [
    The 16S rRNA gene structure. Conserved regions anchor primers;
    variable regions resolve taxonomy. V3 — V4 (~460 bp) is the
    short-read default; full-length sequencing on PacBio HiFi
    promotes the assay to species and strain resolution.
  ],
) <fig:16s-structure>

The variable regions are conventionally numbered V1 through V9. V3 — V4
is the most commonly amplified segment for short-read sequencing
because it fits 2×250 bp paired-end Illumina with overlap, carries
enough variation for genus-level taxonomy, and works well across the
broad range of bacterial taxa that interest most studies. V4 alone is
shorter and amplifies a slightly different range of taxa; V1 — V2
captures different lineages and misses *Bifidobacterium*; full-length
V1 — V9 is now possible on PacBio HiFi and delivers species or
sometimes strain-level resolution. Pick the region with care: the
choice of primer is the dominant taxonomic bias in any 16S study, and
results from different primer sets are not freely comparable.

=== Pace, 1985 and the Culture-Independent Move

The trick of using ribosomal-RNA sequencing to characterise microbial
communities without first culturing the organisms was pioneered by
Norman Pace and colleagues in the mid-1980s. The earlier paradigm —
isolate, culture, identify — was bottlenecked at culturing: more than
ninety per cent of the bacteria in any environmental sample do not
grow on standard laboratory media. Pace, Stahl, Lane, and Olsen's
1985 *PNAS* paper showed that the 5S rRNA gene could be sequenced
directly from environmental DNA, and the obvious extension to 16S
followed within a few years. The line from Pace's work runs through
the rise of "molecular microbial ecology" in the 1990s, through QIIME
in 2010, into the modern industry. The single most important move was
to give up on culturing.

=== DADA2 and Amplicon Sequence Variants

Once the amplicon is sequenced, a 16S workflow has to decide which
reads represent the same biological sequence and which are sequencer
errors. The historical answer was *Operational Taxonomic Units* — read
clusters at 97 % nucleotide identity, on the rough theory that
typical sequencer error rates are at the percent level and
within-species variation in 16S is usually below it. UCLUST and CD-HIT
ran the clustering; mothur and QIIME packaged the pipeline. The
problem with 97 % OTUs is that they erase real, biologically
meaningful single-base differences between closely-related taxa, and
they merge two species that genuinely sit at 98 % 16S identity into
one OTU.

The modern answer is *DADA2* (Callahan, McMurdie, Rosen, Han, Johnson,
and Holmes, *Nature Methods* 2016). DADA2 models the Illumina error
profile explicitly — a substitution rate per cycle, per quality score,
per substitution type — and uses the model to decide which reads are
sequencer-corrupted versions of a true sequence and which are genuinely
new. The output is a set of *Amplicon Sequence Variants*: exact
sequences distinguishable from sequencer noise at single-base
resolution. ASVs replace OTUs in any new analysis and resolve
single-base differences that 97 % clustering used to merge.

#figure(
  image("../../diagrams/lecture-23/02-dada2.svg", width: 95%),
  caption: [
    The DADA2 workflow. Quality-filter, learn per-sample error rates,
    denoise to ASVs, merge paired-end reads, remove chimeras, assign
    taxonomy, build a tree. The error model is what makes single-base
    resolution possible.
  ],
) <fig:dada2>

The DADA2 pipeline in practice runs roughly as follows. Quality-filter
reads at a per-base score (typical: trim where Q drops below 25 — 30)
and at maximum-expected-error thresholds. Learn per-sample error rates
by parametric estimation against a candidate set of high-abundance
sequences. Apply the *core denoising step*: cluster reads into ASVs
under the learned error model, with each candidate ASV's
likelihood ratio against being a sequencer-corrupted version of a more
abundant sequence as the deciding criterion. Merge the paired-end
reads (the overlapping region resolves single-base
ambiguities). Remove chimeras — PCR artefacts where one extension
template-switched during amplification — by checking for a
two-parent-sequence concatenation. Assign taxonomy by naive Bayes or
k-mer-based classification against a reference 16S database. Build a
phylogenetic tree of the ASVs (FastTree is the default) so that
phylogeny-aware diversity metrics can be computed downstream.

#note[
  The denoising step is the part of DADA2 that earns its keep. The
  algorithm models each ASV as the centre of a cloud of variant reads
  generated by Illumina error, and the likelihood of a particular
  read being a member of that cloud rather than a separate ASV is
  computed under the per-sample error model. This is parametric
  empirical Bayes — a step up from the heuristic identity-threshold
  clustering of 97 % OTUs — and it is what lets ASVs distinguish two
  *Lactobacillus* strains that differ by one nucleotide in V4. The
  trade-off is that DADA2 is more sensitive to actual sequencer
  problems: a sample with high error gets fewer ASVs, sometimes
  drastically fewer, which is itself useful diagnostic information.
]

=== QIIME 2 and Reference Databases

QIIME 2 (Quantitative Insights Into Microbial Ecology, Bolyen et al.,
*Nature Biotechnology* 2019) is the dominant umbrella platform for
16S analysis. It descends from the original QIIME (Caporaso et al.,
*Nature Methods* 2010), which packaged the early 2010s pipeline. QIIME 2
adds a workflow-engine flavour: every artefact is versioned, every
analysis step is traceable, and the platform integrates DADA2 as the
default denoiser, vsearch and Deblur as alternatives, and a suite of
diversity and ordination tools. The Galaxy web frontend makes it
runnable without command-line skills.

Taxonomy in 16S is downstream of a curated reference database. Three
databases compete. *SILVA* (Quast et al., 2013, updated quarterly) is
the most comprehensive and the most-used default. *GTDB* (the Genome
Taxonomy Database; Parks et al., *Nature Biotechnology* 2018, 2020)
imposes a phylogenetically coherent taxonomy on whole-genome data
and resolves several long-standing paraphyly problems in older
taxonomies (the genus *Clostridium*, for instance, was split across
twenty-nine distinct groups in the 2018 update). *Greengenes* was the
QIIME 1 default but was effectively frozen at 2013; the Greengenes2
release in 2023 revived it with GTDB-aligned taxonomy. For new
analyses SILVA is the safe default; GTDB is the choice when
phylogenetic rigour matters or when the analysis combines amplicon
data with whole-genome data from the same project.

The output of a DADA2/QIIME 2 run is a set of artefacts that should
feel familiar to anyone who has worked with single-cell RNA-seq. The
*ASV table* is a samples-by-ASVs count matrix. The *taxonomy table*
maps each ASV to a kingdom/phylum/class/order/family/genus/species
assignment with a confidence score. The *phylogenetic tree* of the
ASVs supports phylogeny-aware diversity. The *sample metadata* links
each sample to clinical or experimental variables. With these four
artefacts in hand, every downstream analysis in the chapter is some
function of them.

=== The 2024 Frontier

Three threads matter at the 2024 frontier of amplicon work. The first
is *long-read 16S*. Full-length sequencing on PacBio HiFi promotes the
resolution from genus to species and sometimes strain because the
extra variable regions distinguish lineages that V3 — V4 cannot. The
catch is the cost: about ${\$}80$ per sample, versus ${\$}30$ for short-read
V3 — V4. The trade-off is favourable when the study has a small
cohort and needs strain-level resolution. The second thread is
*multi-region amplicons* — sequencing V1 — V2, V3 — V4, and V6 — V8 in
parallel and combining the calls — which captures lineages that any
single primer pair misses. The third is *learning the error model*:
some recent work uses the read-level confidence outputs of nanopore
basecallers as input to a refined denoiser, with the goal of bringing
nanopore-amplicon 16S into the practical range.


== Shotgun Metagenomics <sec:ch23-shotgun>

Shotgun metagenomics dispenses with the PCR amplification step and
sequences total DNA from a sample. Reads come from every organism
present — bacteria, archaea, fungi, viruses, the host — in proportion
to their DNA abundance. The data are richer and the analyses are
correspondingly heavier. The standard pipeline has three branches off
a common QC-and-host-removal stem.

#figure(
  image("../../diagrams/lecture-23/03-shotgun-pipeline.svg", width: 95%),
  caption: [
    The shotgun metagenomics pipeline. Three analysis branches share a
    common QC and host-DNA-removal stem: taxonomic classification,
    functional profiling, and metagenome-assembled-genome recovery.
  ],
) <fig:shotgun-pipeline>

Step one is read quality control and adapter trimming, usually with
`fastp` or `BBDuk`. Step two is *host-DNA removal*: map the reads
against the host reference (human, mouse, plant, whatever applies)
with `bowtie2` or `bwa-mem2` and discard the mappers. KneadData
packages the bioBakery flavour of this step. Host contamination is
the dominant cost driver in shotgun work — skin samples can be
eighty to ninety-nine per cent human DNA, vaginal samples thirty to
seventy per cent, stool typically under one per cent — and the cost
of usable depth scales accordingly. Step three is taxonomic
classification, which has two competing approaches. Step four is
functional profiling, which depends on having a taxonomic backbone in
place. Step five, optional but valuable when biomass allows, is
metagenome-assembled-genome (MAG) recovery.

=== Kraken2: k-mer Classification

*Kraken2* (Wood, Lu, Langmead, *Genome Biology* 2019) is a k-mer-based
exact classifier. The original Kraken (Wood and Salzberg, 2014) built
a hash table mapping every k-mer in a reference database to the
lowest common ancestor (LCA) of the taxa that contained it; Kraken2
replaced the dense k-mer index with a more compact minimiser-based
scheme and reduced the memory footprint by about an order of
magnitude. The classification itself is simple. For each read,
enumerate its k-mers (default $k = 35$); for each k-mer, look up the
LCA of all reference genomes that contain it; aggregate the
per-k-mer votes into a single read-level call by taking the LCA of
the votes themselves. @fig:kraken-lca walks the steps for a single
sixty-base read.

#figure(
  image("../figures/ch23/f2-kraken-kmer-lca.svg", width: 100%),
  caption: [
    Kraken2 classifies a sixty-base read by decomposing it into
    overlapping thirty-five-mers, looking up the lowest common
    ancestor of each k-mer's reference matches, and aggregating the
    per-k-mer votes into a single read-level call.
  ],
) <fig:kraken-lca>

The engineering geometry is approximate-nearest-neighbour search in
$4^k$-dimensional Hamming space. The reference database is the
indexed dictionary; a query read is decomposed into k-mers; each
k-mer's nearest neighbour in the database provides a taxonomic
identity; the LCA aggregation handles the case where one k-mer
matches several genomes. The same geometry — hashed k-mer lookup with
LCA fallback — appears in MinHash-based sketching (Mash, sourmash),
in DNA-LM tokenisation, and in older locality-sensitive-hashing
literature. Kraken2 is fast: typical runtimes are minutes per sample
on a modern compute node, and the rate-limiting step is database
build, not query.

*Bracken* (Bayesian Reestimation of Abundance with KrakEN; Lu et al.,
2017) post-processes Kraken2 output. Kraken2 emits per-read calls;
many reads end up assigned at the genus or family level because their
k-mers were ambiguous. Bracken redistributes that "stuck-at-the-LCA"
mass down the taxonomy by Bayesian inference: it estimates the
genome-length-corrected abundance of each species under a probabilistic
model of where the ambiguous reads "really" came from. The output is
a clean species-level abundance table that matches what most downstream
tools expect.

#tip[
  When you compare Kraken2/Bracken output to MetaPhlAn output and they
  disagree at the species level — they will — the disagreement is
  almost always informative. MetaPhlAn is conservative; it calls a
  species only when the marker genes hit, and reports zero otherwise.
  Kraken2 is aggressive; every k-mer in the read contributes. The
  truth on a clean sample usually sits between them, and the
  comparison is one of the best free QC checks a shotgun pipeline
  provides.
]

=== MetaPhlAn: Marker-Gene Classification

*MetaPhlAn* (Segata et al., *Nature Methods* 2012; Beghini et al.,
*eLife* 2021 for MetaPhlAn 4; Blanco-Míguez et al., *Nature
Biotechnology* 2023) takes a different cut at the same problem.
Instead of indexing every k-mer of every reference genome, it
precomputes a curated set of *clade-specific marker genes* — about five
genes per species, chosen to be ubiquitous in members of the clade and
absent from non-members. Classification is a Bowtie2 alignment of
reads against the marker database; relative abundance is the per-clade
marker coverage. MetaPhlAn 4 covers about 26,970 species, including
recently discovered uncultured species recovered from MAGs.

MetaPhlAn is slower than Kraken2 — Bowtie2 alignment is roughly an
order of magnitude slower than hashed k-mer lookup — but it is more
accurate at the species level for well-characterised organisms, and
its outputs are easier to interpret biologically because they map to
specific, named genes. The standard 2024 workflow runs both: Kraken2
+ Bracken for breadth (catch everything, including uncultured taxa)
and MetaPhlAn for clinical-grade species-level reporting on the
known clades.

=== HUMAnN: Functional Profiling

*HUMAnN* (the HMP Unified Metabolic Analysis Network; Franzosa et al.,
*Nature Methods* 2018; current release HUMAnN 3) answers the *what
can they do* question. The pipeline maps reads to two reference
resources: *ChocoPhlAn*, a per-species pangenome of clade-specific
genes assembled from sequenced genomes, and *UniRef90*, a clustered
non-redundant protein database covering essentially the
characterisable bacterial protein universe. Reads that map to a
ChocoPhlAn species hit a per-species gene; reads that fail to map to
ChocoPhlAn fall back to a translated nucleotide search against
UniRef90 to catch genes from species not in ChocoPhlAn. The output is
a gene-family-by-sample matrix, then a pathway-by-sample matrix
produced by mapping gene families to MetaCyc pathways.

The HUMAnN output supports a kind of analysis that the taxonomic
pipelines cannot do alone. A microbiome dominated by *Bacteroides*
species can be functionally different from one dominated by *Prevotella*
species even when the relative abundances at the species level look
similar — the pathway repertoires differ. Functional profiling is
the move that lets a microbiome paper claim more than presence and
absence. The bioBakery 4 stack (KneadData → MetaPhlAn 4 → HUMAnN 3)
is the canonical pipeline, and the same operators run it on
everything from infant-microbiome studies to environmental sewage.

=== MAG Recovery

When the biomass is high enough — gut, soil, ocean — shotgun data
support the assembly of near-complete bacterial genomes directly from
the metagenome, with no isolation step. A *Metagenome-Assembled
Genome* (MAG) is the result. The recovery pipeline runs in three
stages. *Assemble* the reads with `megahit` or `SPAdes-meta`, both
of which use multi-k de Bruijn graphs tuned for the
variable-coverage problem of mixed communities (some organisms are
deep, others are shallow, and the assembler has to handle both).
*Bin* the resulting contigs by composition and coverage: contigs from
the same genome share tetranucleotide-frequency profiles and have
correlated read-depth across samples. *MetaBAT 2* (Kang et al.,
*PeerJ* 2019), *CONCOCT*, and *MaxBin 2* are the canonical binners;
each clusters contigs in a feature space that combines composition
and coverage, with cluster cardinality determined empirically.

#figure(
  image("../../diagrams/lecture-23/09-mag-recovery.svg", width: 95%),
  caption: [
    MAG recovery clusters assembled contigs by tetranucleotide
    composition and read coverage; CheckM then assesses each bin's
    completeness and contamination before downstream use.
  ],
) <fig:mag-recovery>

The third stage is quality assessment. *CheckM* (Parks et al.,
*Genome Research* 2015) evaluates each bin by counting the presence
and number of copies of a panel of single-copy marker genes specific
to the lineage. Completeness is the fraction of expected markers
present; contamination is the fraction of markers present in more
than one copy. A *high-quality MAG* by community convention is more
than 90 % complete with less than 5 % contamination; a *medium-quality
MAG* is more than 50 % complete with less than 10 % contamination
(Bowers et al., *Nature Biotechnology* 2017). The *Genomes from Earth's
Microbiomes* (GEM) catalogue (Nayfach et al., *Nature Biotechnology*
2021) provided roughly fifty thousand MAGs from about ten thousand
environmental samples, more than doubling the cultivated bacterial
genome catalogue at the time. The *Unified Human Gastrointestinal
Genome* collection (Almeida et al., 2021) added similar coverage of
the gut. Taxonomy on MAGs runs through GTDB-Tk, which classifies a
bin against the GTDB reference tree by placement on a concatenated
marker-gene phylogeny.

#warn[
  A MAG is not a sequenced genome. It is a consensus assembly from a
  mixed-organism short-read library, and the consensus suppresses
  intra-species variation that may matter biologically. Closely related
  strains in the same sample collapse into a single MAG; rare strains
  fail to bin and disappear into the assembly graph as unclassified
  contigs. The number to internalise is roughly that a "high-quality"
  MAG with 95 % completeness is missing about 250 kb of a typical
  3 Mb genome, and the missing pieces are not random — they are the
  hard-to-assemble repeat regions where lateral gene transfer most
  often lives. MAG-based ARG calls are therefore systematically
  conservative; the ARGs that move most are the ones MAGs miss most.
]


== Diversity Metrics <sec:ch23-diversity>

Diversity is the field's central summary statistic. Population
genetics has $F_("st")$ and Nei's distance; microbiome science has
Shannon entropy, Simpson's index, Bray-Curtis, and UniFrac. The
generalisation is direct: replace the allele with the taxon, the
population with the sample, and the metrics carry over.
Diversity comes in two kinds. *Alpha diversity* is within-sample —
how varied is this one community? *Beta diversity* is between-sample
— how different are these two communities? The third kind, *gamma
diversity*, is across the entire pooled dataset and is rarely used
day to day.

#figure(
  image("../../diagrams/lecture-23/04-diversity.svg", width: 95%),
  caption: [
    Diversity metrics. Alpha measures within-sample richness and
    evenness; beta measures between-sample dissimilarity. Shannon is
    information-theoretic, UniFrac is phylogenetic, Bray-Curtis is
    compositional.
  ],
) <fig:diversity>

=== Alpha Diversity

Three alpha-diversity metrics matter. *Richness* $S$ is the count of
distinct taxa observed in the sample, full stop. It is the simplest
diversity statistic and the most sensitive to sequencing depth —
deeper sequencing finds more rare taxa, so richness comparisons
across samples with unequal depth are misleading unless the data are
rarefied (subsampled to the lowest common depth) or unless
depth-adjusted estimators (Chao 1, ACE) are used.

*Shannon entropy* $H'$ — the same Shannon entropy from information
theory — combines richness and evenness in a single number:

$ H' = -sum_(i=1)^S p_i log p_i $

where $p_i$ is the relative abundance of taxon $i$. The logarithm is
conventionally base $e$ (units: nats) in ecology, base 2 (bits) in
information theory; mixing them silently is a frequent error. Shannon
entropy is maximised when all taxa are equally abundant ($H' = log S$)
and minimises at zero when one taxon dominates entirely. It is
the matched-filter version of the diversity question: an information-
theoretic quantity for an information-theoretic problem.

*Simpson's index* $D$ — strictly, Gini-Simpson — has a more direct
probabilistic interpretation:

$ D = 1 - sum_(i=1)^S p_i^2 $

It is the probability that two reads drawn uniformly at random from
the sample come from different taxa. Like Shannon, it combines
richness and evenness, but it weights evenness more heavily — adding
a rare taxon barely moves Simpson, while it moves Shannon
substantially.

*Faith's phylogenetic diversity* (Faith, 1992) generalises richness by
summing the branch lengths of the phylogenetic tree spanned by the
observed taxa: $"PD" = sum_(b in "tree"(S)) ell_b$. Two samples with the
same Shannon score can differ in PD if one has taxa from many distinct
lineages while the other has many cousins from a single clade. PD is
the metric of choice when biological function correlates with
phylogenetic depth — when "two distinct phyla" matters more than "two
distinct ASVs."

@fig:diversity-math walks the math on two example communities.

#figure(
  image("../figures/ch23/f3-diversity-math.svg", width: 100%),
  caption: [
    Alpha and beta diversity metrics worked numerically on two
    contrasting communities. The first has ten taxa at near-uniform
    abundance; the second has two taxa at 80/20. The formulas and
    their per-community values are tabulated side by side for hand
    verification.
  ],
) <fig:diversity-math>

=== Beta Diversity

Beta-diversity metrics quantify pairwise dissimilarity between
samples. *Bray-Curtis dissimilarity* (Bray and Curtis, 1957) is the
compositional workhorse:

$ d_("BC")(a, b) = (sum_i |a_i - b_i|) / (sum_i (a_i + b_i)) $

It ranges from 0 (identical communities) to 1 (no shared taxa) and
treats abundance differences linearly. *Jaccard* is the
presence-absence cousin — one minus the size of the intersection
divided by the size of the union — and ignores abundance entirely.

*UniFrac* (Lozupone and Knight, *Applied and Environmental
Microbiology* 2005) is phylogeny-aware. Two samples are compared on
the phylogenetic tree of their ASVs. The *unweighted UniFrac*
distance is the fraction of branch length present in only one of the
two samples:

$ d_U(a, b) = (sum_(b in "unique"(a, b)) ell_b) / (sum_(b in "tree"(a union b)) ell_b) $

*Weighted UniFrac* additionally weights each branch by the abundance
difference. Unweighted is the right choice when the question is
about presence and absence; weighted is the right choice when
abundance carries the signal.

#figure(
  image("../../diagrams/lecture-23/10-unifrac.svg", width: 90%),
  caption: [
    UniFrac distance between two samples is the fraction of the
    phylogenetic-tree branch length present in only one. Bray-Curtis
    is composition-only; UniFrac credits shared taxa less when they
    are phylogenetically close.
  ],
) <fig:unifrac>

UniFrac and Bray-Curtis often disagree, and the disagreement is
diagnostic. Two samples that share their taxa with cousins (close on
the tree) will look more similar in UniFrac than in Bray-Curtis; two
samples that share their taxa with strangers (distant on the tree)
will look more similar in Bray-Curtis than in UniFrac. The choice of
metric encodes a biological assumption about what kind of similarity
counts.

Beta-diversity matrices feed two main downstream tools. *Principal
Coordinates Analysis* (PCoA) is a metric scaling of the
sample-by-sample distance matrix into low-dimensional Euclidean space,
where the eigenvectors of the centred distance matrix give the
coordinates. *PERMANOVA* (Anderson, 2001) is the
permutation-based ANOVA on the distance matrix: it partitions
variance among groups defined by sample metadata and tests
significance by shuffling group labels. PCoA gives the picture;
PERMANOVA gives the $p$-value.

#note[
  Diversity is the microbiome field's preferred summary because the
  raw count matrix is too high-dimensional to interpret directly. A
  typical 16S dataset has hundreds to thousands of ASVs, most of them
  rare, most of them with no individually interpretable name. Reducing
  the matrix to a few alpha-diversity scalars per sample and a single
  pairwise distance matrix is a deliberate dimensionality reduction
  that lets the analyst reason about cohort-level patterns.
  Whole-community summaries hide species-level detail; species-level
  analyses lose the gestalt. Both views are needed and neither is
  enough.
]

=== Compositional Data Analysis

Here is the trap that every microbiome paper from 2008 onward had to
learn. *Relative abundances are not absolute abundances.* The total
DNA in a sample is normalised by sequencing depth, so the relative
abundance of species A is a function of the absolute abundance of
every other species in the sample. If species A doubles in absolute
abundance and species B is unchanged, A's relative abundance does
*not* double — it shifts upward, but so do B's relative abundance and
every other taxon's. If species A is unchanged in absolute abundance
and species B halves, A's relative abundance *increases* — and a naive
t-test will report A as significantly enriched even though nothing
happened to A biologically.

@fig:comp-trap walks three scenarios with identical relative-abundance
vectors but completely different absolute-abundance biology. A naive
test calls them all the same and is wrong in at least two of three.

#figure(
  image("../figures/ch23/f4-compositional-trap.svg", width: 100%),
  caption: [
    Three biological scenarios that produce identical relative-abundance
    vectors. A naive t-test treats them identically; only a
    compositional-aware test, anchored to a stable reference taxon,
    can separate them.
  ],
) <fig:comp-trap>

The mathematical framing is *compositional data analysis* (CoDa),
formalised by John Aitchison in 1982 (*The Statistical Analysis of
Compositional Data*). Compositions live on the simplex
$Delta^(D-1) = \{(x_1, dots, x_D) : x_i > 0, sum_i x_i = 1\}$, and the
appropriate algebra on the simplex is not the ordinary vector algebra
of Euclidean space. The *centred log-ratio* (CLR) transform maps a
composition into a $(D-1)$-dimensional Euclidean space where ordinary
statistics apply:

$ "clr"(x)_i = log(x_i) - (1/D) sum_(j=1)^D log(x_j) $

CLR distances are scale-invariant: multiplying every component by a
constant leaves all CLR values unchanged. That is exactly the
property a microbiome statistic should have — sequencing depth varies
sample to sample, and a metric that depends on depth confounds
biology with library size.

*ANCOM-BC* (Analysis of Compositions of Microbiomes with Bias
Correction; Lin and Peddada, *Nature Communications* 2020) is the
current standard differential-abundance test under compositional
constraint. It models log-transformed abundances with a per-sample
bias term estimated from the data, applies a multiple-testing
correction (Benjamini-Hochberg), and reports per-taxon effects in
log-fold-change units relative to the estimated reference frame.
*ALDEx2* (Fernandes et al., 2014) takes a related approach using
Dirichlet-multinomial Bayesian sampling. *songbird* (Morton et al.,
*Nature Communications* 2019) frames differential abundance as a
multinomial regression with explicit reference frames. None of these
tools recovers absolute-abundance change without an anchor; the
reference-frame heuristic — pick a taxon or composite you have reason
to believe is biologically stable, anchor the log-ratio there — is
the practical escape.

#warn[
  Many of the early microbiome literature's findings were artefacts of
  compositional bias. The 2008 — 2015 vintage of papers that used
  t-tests or Wilcoxon tests directly on relative abundances has a
  systematic problem: a substantial fraction of the reported
  "differentially abundant" taxa were not actually changed in
  absolute terms. When ANCOM-BC was applied retrospectively to
  several large public datasets (Lin and Peddada, 2020), between
  twenty and fifty per cent of the original calls did not replicate.
  Reading the older literature requires keeping that revision rate
  in mind.
]


== Host-Microbiome Connections <sec:ch23-host>

The microbiome matters clinically because it changes with disease,
predicts response to therapy, and in some cases causally drives
phenotype. The strongest examples are not metaphors; they are
biomarkers that have been replicated across cohorts. The cleanest one,
and the one most likely to land in clinical practice within the next
decade, is the gut microbiome's effect on cancer immunotherapy
response.

#figure(
  image("../../diagrams/lecture-23/05-microbiome-disease.svg", width: 95%),
  caption: [
    The immunotherapy-response case study. Pre-treatment gut
    microbiomes cluster differently between responders and
    non-responders; *Akkermansia muciniphila* is enriched in
    responders; high-abundance patients have meaningfully better
    survival on anti-PD-1.
  ],
) <fig:microbiome-disease>

=== Immunotherapy Response

Anti-PD-1 and anti-PD-L1 checkpoint inhibitors are the modern
backbone of cancer immunotherapy. Roughly thirty to fifty per cent of
patients respond, depending on tumour type, and the response is often
durable. The variation in who responds is partly explained by
tumour-side factors (mutational burden, neoantigens, PD-L1 expression)
and, since about 2018, partly by gut microbiome composition.

Three papers landed in *Science* on the same day in early 2018. Routy
and colleagues showed in a large French cohort of non-small-cell lung
cancer and renal-cell carcinoma patients that responders to anti-PD-1
had a gut microbiome enriched in *Akkermansia muciniphila*, and that
antibiotic exposure within two months before treatment reduced
response rates. Gopalakrishnan and colleagues showed a comparable
pattern in melanoma patients at MD Anderson. Matson and colleagues
showed the same in melanoma patients at the University of Chicago.
The replication across institutions and tumour types was striking; the
specific taxa varied — *Akkermansia* in some cohorts, *Faecalibacterium*
in others, *Bifidobacterium* in mouse models — but the broad pattern
held.

The mechanism is now reasonably well understood. *Akkermansia* and
similar mucin-degrading commensals modulate the gut immune environment
in ways that prime systemic anti-tumour immunity; mouse experiments
where germ-free animals are colonised with responder-derived microbiota
show better tumour control than animals colonised with non-responder
microbiota. The current clinical wave is *fecal microbiota
transplantation* (FMT) as an adjuvant to checkpoint inhibitors:
Phase I/II trials (Davar et al., *Science* 2021; Baruch et al.,
*Science* 2021) have shown response in previously refractory
melanoma patients. The 2030s clinical picture is likely to include
microbiome modulation as a routine companion to immunotherapy.

#matters[
  The immunotherapy story is the moment microbiome research transitioned
  from "associative" to "clinically actionable." The chain runs
  biomarker → mechanism → intervention, and each link has been
  closed in the laboratory and in early-phase trials. The microbiome
  is now a tractable, modifiable disease-modifying organ — a
  twenty-first-century version of the hospital pharmacy that
  Hippocrates would have recognised in spirit but not in technique.
  The chain breaks if the underlying compositional and statistical
  analyses are not done correctly. The earlier sections of this
  chapter are the prerequisite for taking these clinical claims
  seriously.
]

=== Inflammatory Bowel Disease

Crohn's disease and ulcerative colitis show reproducible gut-microbiome
signatures: reduced diversity (Shannon entropy roughly 1.5 nats lower
than healthy controls), loss of obligate anaerobes —
*Faecalibacterium prausnitzii* in particular, the most-abundant
bacterium in healthy guts — and bloom of pathobionts (*Escherichia*,
*Klebsiella*, *Ruminococcus gnavus*). The *Integrative Human
Microbiome Project* IBD cohort (Lloyd-Price et al., *Nature* 2019)
provided about 1,700 longitudinal samples; the *MetaCardis* and
*HMP1/HMP2* consortia added thousands more. Machine-learning
classifiers reach about 85 % accuracy distinguishing IBD from healthy
controls; whether the signature is causal or consequence of
inflammation is still debated.

=== Other Clinical Connections

The clinical microbiome literature is broader than the immunotherapy
and IBD cases. *Porphyromonas gingivalis*, an oral pathogen of
periodontitis, is enriched in the saliva of pancreatic-cancer
patients up to a decade before diagnosis (Fan et al., *Gut* 2018);
the prospective implication is that oral-microbiome surveillance could
be early detection. *Fusobacterium nucleatum* is enriched in
colorectal-cancer tissue (Castellarin et al., 2012) and promotes
tumour growth in mouse models. *Bacterial drug metabolism* converts,
activates, or inactivates roughly thirty per cent of orally
administered drugs (Spanogiannopoulos et al., *Nature Reviews
Microbiology* 2016) — including the cardiac glycoside digoxin, the
antiviral brivudine, and metformin in type-2 diabetes. *Microbiome-GWAS*
(mGWAS) studies have identified host genetic variants associated
with microbial abundance — the lactose-tolerance variant predicts
*Bifidobacterium* abundance; ABO blood group predicts *Bacteroides*
abundance — suggesting that the microbiome itself is heritable to a
modest degree.

=== Cohort-Scale Resources

Public microbiome data scale into the hundreds of thousands of
samples. The *Human Microbiome Project* (HMP and the follow-on
iHMP) sampled multiple body sites in about 300 and then 1,700
individuals respectively. *MetaCardis* (Forslund et al., *Nature*
2015 onward) covers about 2,000 individuals with cardiometabolic
phenotypes. *TwinsUK Microbiome* provides identical-twin pairs for
heritability analysis. *MGnify* at the EBI aggregates raw amplicon
and shotgun data across studies; the SRA holds the long tail. The
combination — public data plus open pipelines plus a methodological
canon — has made microbiome research one of the most
data-rich subfields of computational biology.


== Frontier Topics <sec:ch23-frontier>

Three frontier topics matter enough to flag.

*Antibiotic resistance gene tracking.* Shotgun metagenomics maps the
mobile-resistome of clinical samples and environmental reservoirs.
The *Comprehensive Antibiotic Resistance Database* (CARD) and the
NCBI *AMRFinderPlus* tool catalogue known resistance genes; reads
that hit them, or assembled contigs that carry them, identify the
presence of resistance markers. The interesting time-series is
hospital sewage: ARGs move between species through lateral gene
transfer at rates that surveillance cultures cannot capture, and
metagenomic monitoring of sewage detects ARG abundance shifts in
response to local antibiotic-stewardship policies. The signal is
strong enough that several European cities have run prospective
sewage-metagenomic surveillance for ARG dynamics.

#figure(
  image("../../diagrams/lecture-23/12-arg-tracking.svg", width: 90%),
  caption: [
    Antibiotic resistance gene abundance traced in hospital sewage
    metagenomes over two years. A policy change at month 8 leaves
    visible inflections in two ARG classes.
  ],
) <fig:arg>

*Strain-level resolution.* Within-species variation matters — pathogenic
*E. coli* and commensal *E. coli* differ in a handful of accessory
genes (Shiga toxin, intimin, fimbriae) that determine clinical
behaviour. *StrainPhlAn* (Truong et al., 2017), *MIDAS* (Nayfach et
al., 2016), and *inStrain* (Olm et al., *Nature Biotechnology* 2021)
infer strain composition from shotgun data by genotyping marker genes
across samples. Long-read shotgun (PacBio HiFi, ONT) has begun to
deliver full-length strain genomes directly from environmental
samples. The 2024 frontier is whether strain-level signatures can
predict clinical phenotype where species-level signatures cannot.

#figure(
  image("../../diagrams/lecture-23/11-strain-resolution.svg", width: 85%),
  caption: [
    Short versus long-read 16S resolution. Genus-level becomes
    species- and strain-level when the full-length 1,500-bp gene is
    sequenced.
  ],
) <fig:strain>

*Foundation models for microbial proteins.* Protein language models
trained on UniRef include enormous microbial subsets — *ESM-2*,
trained on UniRef50, is roughly 60 % microbial by composition. The
embeddings transfer to microbiome-specific tasks: functional
annotation of MAG-recovered genes, antimicrobial-resistance gene
prediction, and metabolite-pathway inference. *Evo* (Nguyen et al.,
*Science* 2024) is a DNA-LM trained on millions of microbial genomes
at the megabase scale; its generative outputs include plausible
designed CRISPR systems and operon structures. The microbiome's data
abundance — public sequence on the petabyte scale — makes it a
natural home for foundation-model approaches that data-limited
single-species domains cannot support.


== Tools, Pitfalls, and Practice <sec:ch23-practice>

The 2024 standard 16S stack is QIIME 2 (full pipeline) or DADA2-in-R
(more flexible at the cost of more manual integration). Diversity and
ordinations run through `phyloseq` (R), `scikit-bio` (Python), or
QIIME 2 itself. Differential abundance runs through ANCOM-BC,
MaAsLin 2 (a linear-model framework that handles longitudinal designs),
or songbird.

The 2024 standard shotgun stack is *bioBakery 4*: KneadData for QC and
host removal, MetaPhlAn 4 for taxonomy, HUMAnN 3 for function. For
reproducibility-critical work, *nf-core/mag* and *Sunbeam* are
Nextflow- and Snakemake-based pipelines that pin tool versions,
parameter sets, and database releases. MAG recovery runs through
*nf-core/mag* or *Atlas*. CheckM and GTDB-Tk gate downstream use.

#figure(
  image("../../diagrams/lecture-23/06-workflow.svg", width: 95%),
  caption: [
    The microbiome analysis workflow from sample to publication. Each
    band carries reproducibility hazards; database-version pinning
    and primer specification should appear in every methods section.
  ],
) <fig:ch23-workflow>

The pitfalls are familiar by now but worth re-listing in one place.

*Sample collection* matters more than any subsequent step. Storage
at $-80$ °C is the gold standard; OMNIgene-GUT and similar stabilising
buffers are acceptable for ambient-temperature shipping;
room-temperature storage shifts composition measurably within twenty-four
hours. Cross-contamination from sample handling is a real source of
batch-confounded signal.

*Sequencing* introduces three biases. Low-biomass samples are
dominated by host DNA, which sets the effective per-sample cost in
the shotgun branch. PCR amplification in the 16S branch preferentially
amplifies some taxa over others — the bias is reproducible per primer
set but not per study. Sequencing depth requirements are about ten
thousand reads per sample for 16S (DADA2 needs enough reads to learn
the error model) and roughly five gigabases per sample for shotgun.

*Bioinformatics* introduces version-dependence. ASVs are objectively
better than 97 % OTUs but legacy data uses OTUs. Reference databases
change taxonomic boundaries between releases — *Clostridium* was
split into many genera in GTDB R03-RS86, and the names of the
resulting clades differ between SILVA and GTDB. Pin database
versions in the methods section; do not assume that re-running an
analysis a year later will reproduce the original calls without
explicit version control.

*Statistical* pitfalls are mostly the compositional-data problem
discussed in @sec:ch23-diversity, plus the usual issues with multiple
testing across thousands of taxa, subject-level confounders (age,
diet, antibiotic exposure, geography), and sample-size requirements.
A working rule of thumb is at least thirty samples per arm for
typical effect sizes, with formal power calculations using real
effect-size priors rather than naive Cohen's d.

#tip[
  The most useful single QC check on any new microbiome dataset is to
  run both Kraken2/Bracken and MetaPhlAn on the same samples, then
  inspect the disagreement matrix at the species level. Disagreements
  are usually informative: a species called confidently by both is
  almost certainly real; a species called only by Kraken2 is often a
  database-induced false positive; a species called only by MetaPhlAn
  is often a marker-gene match for a closely related taxon. The
  comparison is cheap, runs once per project, and catches most upstream
  data problems before they propagate into downstream statistics.
]

=== Reproducibility, Honestly

Microbiome reproducibility is famously poor. Cross-cohort replication
of named taxa as biomarkers is roughly thirty to fifty per cent — much
worse than, say, GWAS replication for common variants. The drivers
are cohort effects (geography, diet, ethnicity, age, sampling
protocol), technical variation (extraction kit, primer choice,
sequencer, library prep), and methodological choices (compositional
versus non-compositional analysis, ASV versus OTU, database version).
The replication rate improves substantially when these are controlled
— same kit, same primers, same database version, compositional-aware
test — but no single study can control for all of them, and the
honest expected replication rate is around fifty per cent for typical
effect sizes through the late 2020s.

The community has begun to standardise. The *MicrobiomeQC* effort
publishes shared positive controls (mock communities of known
composition) that pipelines should recover at known relative
abundance. The *Earth Microbiome Project* shipped a single primer set
and protocol for amplifying environmental samples across thousands of
sites worldwide, and its standardisation is the reason the resulting
catalogue is comparable across studies. *cwl-microbiome* and
*nf-core/mag* commit to a single pipeline definition that runs
identically on different infrastructures. None of these efforts has
solved reproducibility, but each narrows the legitimate sources of
disagreement.


== Summary <sec:ch23-summary>

- A microbiome sample is a mixture of DNA from hundreds of organisms in
  unknown proportions; recovering per-taxon abundances is a community
  signal demixing problem. The two main assays are 16S rRNA amplicon
  sequencing (cheap, genus-resolution, bacteria-only) and shotgun
  metagenomics (expensive, species-strain-resolution, all kingdoms,
  with functional content).
- DADA2 + QIIME 2 is the modern 16S stack. DADA2 infers
  *amplicon sequence variants* — exact sequences distinguishable from
  sequencer error — and replaces the 97 % OTU clusters of the older
  literature. Taxonomy assignment runs against SILVA (default) or
  GTDB (phylogenetically coherent).
- The shotgun stack is bioBakery 4: KneadData (host removal) → MetaPhlAn 4
  (marker-gene taxonomy) or Kraken2 + Bracken (k-mer-LCA classification)
  → HUMAnN 3 (functional profiling). MAG recovery via MetaBAT 2 or
  similar binners, with CheckM as the quality gate.
- Diversity is the field's central summary. Shannon entropy and
  Simpson's index measure alpha (within-sample) diversity;
  Bray-Curtis and UniFrac measure beta (between-sample)
  diversity. UniFrac is phylogeny-aware; Bray-Curtis is composition-only.
  PCoA visualises the distance matrix; PERMANOVA tests group
  differences.
- *Compositional bias is the central statistical hazard of the field.*
  Three completely different absolute-abundance scenarios can produce
  identical relative-abundance vectors, and naive t-tests cannot
  separate them. Centred log-ratio transforms (CLR) and tools that
  build on them (ANCOM-BC, ALDEx2, songbird) are the correct toolkit.
  Reference frames anchored to known-stable taxa or spike-ins recover
  absolute change.
- The gut microbiome predicts immunotherapy response. *Akkermansia
  muciniphila* enrichment in pre-treatment microbiomes correlates with
  anti-PD-1 response across non-small-cell lung cancer, renal-cell
  carcinoma, and melanoma cohorts; faecal microbiota transplantation
  has rescued response in previously refractory patients in
  Phase I/II trials.
- Replication of microbiome biomarkers across cohorts is roughly
  thirty to fifty per cent. The drivers are cohort effects, technical
  variation, and methodological choices. Pipeline standardisation
  (nf-core/mag, Earth Microbiome Project protocols) is closing the
  gap slowly.
- EE framings: metagenomics is community signal demixing; Kraken2's
  k-mer LCA is approximate nearest-neighbour search in $4^k$-dimensional
  Hamming space; Shannon entropy of relative abundances is the
  information-theoretic version of richness; UniFrac is a
  branch-length-weighted set distance on a phylogenetic tree.


== Exercises <sec:ch23-exercises>

#strong[1.] #emph[Diversity by hand.] A faecal sample has the
following ASV read counts: ASV1 = 5000, ASV2 = 3000, ASV3 = 1500,
ASV4 = 400, ASV5 = 80, ASV6 = 20. Compute (a) richness, (b) Shannon
entropy in both nats and bits, (c) Simpson's index, and (d) Pielou's
evenness $J' = H' \/ log S$. Repeat for a rarefied sub-sample of 500
total reads drawn proportionally. Comment on how each metric responds
to the depth change.

#strong[2.] #emph[Bray-Curtis and UniFrac disagree.] Construct two
communities $A$ and $B$ over five taxa $T_1, dots, T_5$ such that
$d_("BC")(A, B) = 0.6$ but $d_U(A, B) < 0.2$, using a small phylogenetic
tree of your design. Show your work. Explain in one paragraph why the
two metrics give such different answers.

#strong[3.] #emph[Compositional rescue.] Generate three synthetic
two-sample datasets that match the three scenarios in @fig:comp-trap.
Apply (a) a naive t-test on relative abundances, (b) a t-test on the
CLR-transformed values, and (c) ANCOM-BC. Report which test flags
which scenario and explain the pattern. A spike-in of a known
constant biomass should be added to one of the three scenarios to
demonstrate the reference-frame rescue.

#strong[4.] #emph[Kraken2 walk-through.] Implement the k-mer LCA
classification algorithm in fifty lines of Python on a toy reference
database of three bacterial genomes ($k = 7$ for tractability). Run
your implementation on a set of simulated 100-base reads with
introduced sequencer errors at 1 % per-base rate; report per-read
classification accuracy. Compare your call to the call from a real
Kraken2 run on the same input.

#strong[5.] #emph[DADA2 versus OTU.] Take a small public 16S dataset
(QIIME 2 ships the Moving-Pictures tutorial as a worked example).
Run the data through (a) DADA2 to produce ASVs and (b) closed-reference
97 % OTU clustering. Compare the resulting tables on the basis of
taxon count, taxon-level concordance with the published call, and the
sensitivity of Shannon entropy across the two representations. Write
one paragraph on which differences are biological and which are
methodological.

#strong[6.] #emph[MAG recovery.] Assemble a small subsampled shotgun
dataset (e.g., one HMP gut sample, downsampled to 500 Mb) with
`megahit`; bin the contigs with MetaBAT 2; assess each bin with CheckM.
Report the completeness/contamination for each bin and identify the
likely taxa via GTDB-Tk. Discuss which bins meet the "high-quality
MAG" threshold and why the others do not.

#strong[7.] #emph[Differential abundance and reference frame.] Take a
public IBD vs healthy 16S cohort (the iHMP-IBD pilot dataset is a
reasonable starting point). Run ANCOM-BC. Then repeat the analysis
with the reference frame anchored to *Faecalibacterium prausnitzii*
explicitly. Report which calls change and which do not. Discuss why
the reference frame matters for biological interpretation when the
disease itself is suspected to alter community-wide biomass.

#strong[8.] #emph[(Open-ended.)] Pick one of the immunotherapy
microbiome papers (Routy et al. 2018, Gopalakrishnan et al. 2018, or
a follow-up cohort study). Audit the analysis pipeline: which
denoiser, which database, which differential-abundance test, which
reference frame, which split protocol. Identify the one
methodological choice you would change in a 2026 reanalysis and
predict in one paragraph how the headline number would shift.


== Further Reading <sec:ch23-further-reading>

- *Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W.,
  Johnson, A. J. A., and Holmes, S. P.* (2016). "DADA2:
  high-resolution sample inference from Illumina amplicon data."
  _Nature Methods_ 13: 581 — 583. The ASV paper. Read with the
  Callahan, McMurdie, Holmes (2017) ASV-versus-OTU follow-up for the
  full argument.
- *Wood, D. E., Lu, J., and Langmead, B.* (2019). "Improved
  metagenomic analysis with Kraken 2." _Genome Biology_ 20: 257. The
  current Kraken2 paper; pair with Lu et al. (2017) on Bracken for
  the abundance-reestimation step.
- *Beghini, F., McIver, L. J., Blanco-Míguez, A., et al.* (2021).
  "Integrating taxonomic, functional, and strain-level profiling of
  diverse microbial communities with bioBakery 3." _eLife_ 10:
  e65088. The bioBakery 3 paper. MetaPhlAn 4 and HUMAnN 3 references
  follow from this.
- *Lozupone, C., and Knight, R.* (2005). "UniFrac: a new
  phylogenetic method for comparing microbial communities." _Applied
  and Environmental Microbiology_ 71: 8228 — 8235. The UniFrac paper;
  short, readable, foundational.
- *Lin, H., and Peddada, S. D.* (2020). "Analysis of compositions of
  microbiomes with bias correction." _Nature Communications_ 11:
  3514. The ANCOM-BC paper. Pair with Gloor et al. (2017) on
  compositional data analysis for the broader framing.
- *Nayfach, S., Roux, S., Seshadri, R., et al.* (2021). "A genomic
  catalog of Earth's microbiomes." _Nature Biotechnology_ 39:
  499 — 509. The GEM catalogue; the moment MAG-based discovery
  surpassed cultivated-genome catalogues in coverage.
- *Routy, B., Le Chatelier, E., Derosa, L., et al.* (2018). "Gut
  microbiome influences efficacy of PD-1-based immunotherapy against
  epithelial tumors." _Science_ 359: 91 — 97. The immunotherapy
  paper. Read with Gopalakrishnan et al. (2018) and Matson et al.
  (2018) for the same-day melanoma replications.
- *Spanogiannopoulos, P., Bess, E. N., Carmody, R. N., and Turnbaugh,
  P. J.* (2016). "The microbial pharmacists within us: a metagenomic
  view of xenobiotic metabolism." _Nature Reviews Microbiology_ 14:
  273 — 287. The drug-metabolism review. The thirty-per-cent number
  in the chapter comes from here.
