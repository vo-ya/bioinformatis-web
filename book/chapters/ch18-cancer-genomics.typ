#import "../theme/book-theme.typ": *

= Cancer Genomics: Integrated Capstone <ch:cancer-genomics>

#matters[
  Cancer is the integrating problem for the first seventeen chapters of
  this book. Every layer of the pipeline you have spent the course
  building — the FASTQ from Chapter 1, the aligner from Chapter 2, the
  variant caller from Chapter 4, the RNA quantifier from Chapters 5 and
  6, the long-read assembler from Chapter 11, the clinical-tiering
  framework from Chapter 14 — meets a problem with stakes attached when
  you point it at a tumour. The biology is the same Darwinian dynamics
  you read about in pre-clinical population genetics, sped up by a
  factor of a thousand and run inside the lifetime of one person. The
  engineering is harder than any single previous chapter because every
  upstream failure mode now has a downstream therapeutic consequence,
  and the same VCF row might decide between immunotherapy, a kinase
  inhibitor, and a clinical trial.
]

A cancer is an evolutionary experiment that the patient is conducting
unwillingly. Somewhere between conception and middle age a single cell
acquires a mutation in a gene that controls its own proliferation, then
acquires a second mutation that disables a check on that proliferation,
then a third that allows it to evade the immune surveillance the body
ordinarily uses to clean up such cells, then a fourth that lets it
ignore the geometric signals telling cells to stay where they are. Each
mutation is a Bernoulli trial against the genome's repair machinery;
each cell division is a chance for one of the trials to succeed. By the
time a lump is palpable the cell has descendants in the billions, those
descendants have diverged from one another in their mutation catalogues,
and the lump is a population of related but non-identical clones with
their own subselection happening inside.

This chapter turns that biology into a bioinformatics workload. We
start with the population-genetic view — cancer as fast-forward
evolution inside a body — and translate the textbook hallmarks
(Hanahan and Weinberg, 2000, 2011, 2022) into the molecular alterations
a sequencer can actually measure. We then walk the cancer-specific
sequencing landscape: tumour-normal pairs, FFPE artefacts, ctDNA and
liquid biopsy, and the panel-versus-exome-versus-genome trade-off that
defines clinical assay choice. The third part walks somatic variant
calling at depth, with Mutect2 as the canonical caller and mutational
signatures as the linear-algebra interlude. Structural variants and
fusions occupy a part of their own because cancer breaks chromosomes
in ways no other field has to handle at scale. Subclonal reconstruction
returns the population genetics; clinical interpretation lands the
sequence in AMP/ASCO/CAP tiers; and the final part walks the path from
a VCF to a drug-gene-cancer-triple decision in a tumour board.

The chapter is the integrated capstone of this book's first arc —
sequencing in Chapter 1, alignment in Chapter 2, variant calling in
Chapter 4, RNA quantification in Chapters 5 and 6, long reads in
Chapter 11, clinical pipelines in Chapter 14, and clinical
interpretation in Chapter 17 — applied to one disease at the depth
where engineering decisions become medical ones.


== Cancer as Applied Evolution <sec:applied-evolution>

The most useful single move in this chapter is to read cancer through
the lens of population genetics. A normal tissue is a population of
cells with shared ancestry, identical-by-descent genomes, and a
recombination-free clonal-only mode of inheritance. A tumour is the
same population running at a higher mutation rate, under a stronger
selection pressure, in a confined geographical space, with no gene
flow in or out. Every concept from Wright-Fisher dynamics carries over.
Variant allele fraction is the allele frequency. Clonal fixation is the
case where a variant reaches frequency one. Selection coefficients
become driver advantages. Drift produces minor subclones the same way
it produces rare alleles in a small population. The metaphor is
exact, and it explains why nearly every method in the rest of the
chapter has a recognisable population-genetic shape.

The Hanahan-Weinberg *hallmarks of cancer* — six in the 2000 paper,
ten in the 2011 update, fourteen in the 2022 revision — are an attempt
to organise that evolutionary endpoint into a list of capabilities the
cancer cell has to acquire. Sustained proliferative signalling. Evasion
of growth suppressors. Resistance to cell death. Replicative immortality.
Induction of angiogenesis. Activation of invasion and metastasis. The
2011 paper added reprogramming of energy metabolism and immune evasion,
plus two enabling characteristics — genome instability and
tumour-promoting inflammation. The 2022 update added phenotypic
plasticity, the role of senescent cells, non-mutational epigenetic
reprogramming, and the polymorphic microbiome. Each hallmark is a
capability; each capability is enabled by specific molecular alterations
that a sequencer can measure.

#figure(
  image("../../diagrams/lecture-18/01-hallmarks.svg", width: 95%),
  caption: [
    The 2022 hallmarks of cancer with the original 2000 capabilities
    on the inner ring and the 2011 and 2022 additions on outer rings.
    Each capability is annotated with an example molecular alteration
    the sequencing pipeline can detect.
  ],
) <fig:hallmarks>

The hallmarks framework is the biological vocabulary. The
bioinformatics vocabulary, the one you will use day-to-day, is a
shorter list: *driver* mutations that confer a selective advantage on
the cell line; *passenger* mutations that happen to ride along; clonal
expansion that fixes a driver in a lineage; subclonal divergence that
produces the heterogeneity the next round of selection acts on. A
typical solid tumour carries between one thousand and one hundred
thousand somatic mutations depending on cancer type — pancreatic
adenocarcinoma sits near the bottom, ultraviolet-driven melanoma near
the top — and only five to fifteen of those are drivers in the
classical sense.

Distinguishing the drivers from the passengers is one of the field's
core technical problems and has produced an industry of methods. The
simplest is recurrence: a variant seen in many patients' tumours, in
the same gene, in the same hot-spot, is overwhelmingly likely to be a
driver. COSMIC's curated Cancer Gene Census collects these; cBioPortal
exposes the per-patient counts; TCGA and PCAWG provide the
pan-cancer denominator. The next step up is statistical: MutSig
(Lawrence et al., *Nature* 2013) and its successor MutSigCV
estimate the per-gene background mutation rate from local context and
flag genes mutated more often than that null. MuSiC, OncodriveFML, and
dN/dS-based methods address the same problem from different
angles. Functional annotation closes the loop — a missense in BRAF's
V600 position is a driver because the biochemistry is settled, even
on first observation. The hottest of the hot-spots are essentially
single-base prior beliefs.

#figure(
  image("../../diagrams/lecture-18/02-drivers-clonal.svg", width: 95%),
  caption: [
    Drivers, passengers, and the clonal tree they imply. Ten drivers
    out of ten thousand mutations on the left; a phylogenetic tree
    with the drivers placed on branches by the order in which they
    rose to fixation on the right. Subclonal branches carry later
    mutations — the canonical site of treatment-resistant recurrence.
  ],
) <fig:drivers>

#note[
  Variant allele fraction is the cancer-genomics analogue of allele
  frequency, and the same intuition applies. A clonal variant in a
  pure diploid tumour has VAF $approx 0.5$; a clonal variant in a
  60 %-pure tumour has VAF $approx 0.3$; a variant present in half
  the tumour cells of a 60 %-pure tumour has VAF $approx 0.15$. A VAF
  of $0.05$ in a 20 %-pure sample could be either clonal or a small
  subclone — you cannot tell without an independent purity estimate.
  This ambiguity is why every modern somatic pipeline produces a
  purity / ploidy estimate before it tries to interpret the VAF
  distribution.
]

Genome instability sits at the centre of the hallmark wheel because it
enables all the others. Mechanisms come in several flavours. Point
mutations from replication error, ultraviolet photodimers, tobacco
adducts, oxidative damage, and spontaneous deamination drive most of
the SNV load and produce the per-process spectra you will meet again
in Section 1.4 as mutational signatures. *Microsatellite instability*
follows from mismatch-repair loss — typically inactivation of MLH1,
MSH2, MSH6, or PMS2 — and produces a characteristic insertion-deletion
spectrum at short repeats. *Chromosomal instability* manifests as
whole-chromosome gains and losses, copy-number changes that scale
genome-wide. *Homologous recombination deficiency*, from BRCA1, BRCA2,
PALB2, or related loss, leaves a specific signature of large
deletions, tandem duplications, and a characteristic single-base
substitution pattern (SBS3) that becomes a clinical biomarker for PARP
inhibitor therapy in Section 1.6.


== The Cancer Sequencing Landscape <sec:landscape>

Cancer sequencing differs from the germline analyses of Chapter 4 in
three operational ways that drive every downstream design choice. The
sample is impure — only a fraction of the input DNA is from tumour
cells, with the rest from stroma, immune infiltrate, and adjacent
normal tissue. The sample is often archival — fixed in formalin and
embedded in paraffin years before the sequencing technology used to
read it existed. And the sample is heterogeneous — a single biopsy
contains cells from multiple subclones at different abundances,
each contributing its own mutation set to the bulk read pool. The
combination forces the field into a small set of well-defined
sequencing designs.

The gold standard is *tumour-normal paired sequencing*. A tumour
sample and a matched normal — peripheral blood, saliva, or adjacent
histologically normal tissue — are sequenced separately. Somatic
variants are by definition present in the tumour and absent from the
normal; germline variants are present in both. The subtraction is
clean. The cost is roughly double a single-sample assay because the
patient pays for two libraries, two flow-cell lanes, and two
alignments. Research consortia — TCGA, PCAWG, ICGC, every modern
discovery-cohort study — almost always run tumour-normal because the
rigour is required for publication and the resulting calls feed the
community knowledge bases.

The clinical reality is more often *tumour-only sequencing*. A single
FFPE block is all that exists; no peripheral-blood draw was banked
when the patient was originally biopsied; the assay's commercial
operator wants to minimise the per-sample cost. Tumour-only calling
filters germline variants statistically using gnomAD frequencies and
ancestry priors instead of subtracting against a matched sample. The
loss in specificity is typically five to ten percent. Commercial
panel assays such as FoundationOne CDx and MSK-IMPACT are run
tumour-only with elaborate germline filtering; some newer designs
add *very-shallow* normal sequencing (5× coverage of the matched
blood) as a cheap germline filter that recovers most of the
matched-pair specificity at a fraction of the cost.

#warn[
  Tumour-only false-positive drivers are a documented hazard. A
  TP53 variant called from a tumour-only sample that turns out to be
  a rare germline variant misclassified by the population frequency
  filter is the canonical bug. Every tumour-only pipeline that
  reports actionable variants has to include a manual review step for
  the strongest hits, and a re-sequencing of the matched normal when
  the variant is consequential. The savings of tumour-only are real
  but they are not free.
]

The second design axis is the input material. Most clinical samples
are *FFPE* — formalin-fixed, paraffin-embedded — because that is what
pathology archives have been producing for a century. Formalin
cross-links DNA to itself and to protein, fragments it to a typical
median size near 150 base pairs, and most importantly deaminates
cytosine bases at a rate that scales with the age of the block. The
result, on a sequencer, is an excess of C-to-T transitions — and the
excess is not random but strand-biased, because the deamination
preferentially affects one strand at a time. Modern callers
(Mutect2's `FilterByOrientationBias`, the equivalent step in
Strelka2) detect the strand asymmetry and reject calls that look like
deamination artefacts. Even so, low-allele-fraction C-to-T calls from
FFPE samples are an artefact-versus-subclone judgment call that
clinical pipelines lean toward filtering aggressively.

The third axis is the assay scale. *Targeted panels* sequence 100 to
500 cancer-relevant genes at very high depth (typically $tilde 500
times$) and dominate clinical care because they catch the actionable
variants at a fraction of the cost of broader sequencing. The major
commercial assays — FoundationOne CDx (324 genes), MSK-IMPACT (468
genes), TSO500 (523 genes) — are panel-shaped. *Whole-exome
sequencing* (WES) covers the 1 % of the genome that codes for protein
at moderate depth ($tilde 100 times$) and dominates research-cohort
work. *Whole-genome sequencing* (WGS) covers the entire 3 Gb at
60 to 100$times$ and captures non-coding drivers, structural
variants, and full mutational-signature inference. WGS is heavy in
research and is gaining clinical adoption for cases where SVs or
mutational-signature biomarkers (HRD score from SBS3) drive therapy.

#figure(
  image("../../diagrams/lecture-18/03-sequencing-designs.svg", width: 100%),
  caption: [
    The three dominant cancer-sequencing designs. Tumour-normal pairs
    on the left, FFPE tumour-only in the middle, ctDNA liquid biopsy
    on the right. Coverage targets, library quirks, and intended
    failure modes differ across the three.
  ],
) <fig:designs>

The newer and most operationally interesting design is *liquid biopsy*
from circulating tumour DNA. ctDNA is the small fraction of cell-free
DNA in plasma that originates from dying tumour cells. A blood draw is
non-invasive, can be repeated arbitrarily often, and captures DNA
from every metastatic site simultaneously rather than the single
geographic location a needle biopsy samples. The price is that ctDNA
is often less than $0.1 %$ of total cell-free DNA in early disease,
which forces coverage in the $30,000 times$ range and the use of
unique molecular identifiers (UMIs) to suppress sequencer noise to
the $10^(-5)$ per-base error floor that subclonal detection requires.
ctDNA is now FDA-approved for several use cases — Guardant360,
FoundationOne Liquid — and is the dominant tool for minimal
residual disease (MRD) monitoring after surgery, where the question is
whether any tumour cells survived.

#note[
  An interpretation hazard in ctDNA is *clonal hematopoiesis of
  indeterminate potential* (CHIP) — clonal mutations in
  hematopoietic stem cells, present in healthy people at increasing
  prevalence with age, that contribute DNA to the cell-free pool.
  A DNMT3A or TET2 variant at low VAF in plasma is more likely to be
  CHIP than tumour-derived, and the standard filter is to require the
  same variant in matched buffy-coat DNA. Without that control, MRD
  pipelines produce false-positive "tumour" calls in older patients.
]

The 2025 clinical workflow tends to combine these designs across a
patient's trajectory. Diagnosis uses a high-coverage tumour-only or
tumour-normal panel for fast actionable-variant identification.
Research-grade WES or WGS may be added in academic centres or for
trial-eligibility evaluation. Longitudinal ctDNA monitors residual
disease and the emergence of resistance subclones. Each modality
answers a different question; the integration happens at the molecular
tumour board (Section 1.7).


== Somatic Variant Calling at Depth <sec:somatic-calling>

The germline framework from Chapter 4 assumes that, at every site,
the truth is one of three genotypes: homozygous reference,
heterozygous, or homozygous alternate. Variant allele fraction in a
diploid sample takes one of three values — $0$, $0.5$, $1$ — and the
caller's job is to choose among them given a stack of reads.
HaplotypeCaller's HMM and GATK's variant-quality recalibration are
both built on top of that three-state assumption.

Somatic calling discards every part of that framework. Tumour purity
shifts the expected VAF of a clonal variant from $0.5$ to $0.5 times
"purity"$. Subclonal variants sit at VAFs below the clonal level. Copy
number changes — common in cancer — change the expected VAF further
because a variant on a chromosome present in three copies has a
different expected fraction than diploid. And the noise floor that
matters is the sequencer's, not biology's, because a $5 %$ VAF call
might be real if the underlying subclone is real and pure noise if
the underlying read errors happen to cluster the right way.

#figure(
  image("../../diagrams/lecture-18/04-mutect2.svg", width: 100%),
  caption: [
    The Mutect2 workflow. Reads from tumour and matched normal feed a
    local-assembly module inherited from HaplotypeCaller; a somatic
    likelihood model emits raw calls; a sequence of cancer-specific
    filters — Panel of Normals, orientation bias, contamination,
    tumour-in-normal — removes the dominant false-positive classes.
    The cancer additions are where most real-versus-artefact decisions
    are made.
  ],
) <fig:mutect2>

*Mutect2* (Cibulskis et al., *Nature Biotechnology* 2013 for the
original MuTect; Benjamin et al., *bioRxiv* 2019 for the GATK 4
version) is the field's most-used somatic caller and a good entry
point for understanding the architecture. The first stage is local
reassembly — exactly as in HaplotypeCaller — to produce a set of
candidate haplotypes from the reads in a window around each putative
variant. The second stage replaces the germline likelihood with a
*somatic likelihood* that scores each read against each haplotype
under the assumption that the alternate allele could be at any
fraction between zero and one, rather than the discrete $0, 0.5, 1$
of the germline case. The third stage is a sequence of
cancer-specific filters: a Panel of Normals (PoN) built from at
least forty normal samples processed identically catches recurrent
alignment artefacts; orientation-bias filtering catches FFPE-style
strand-biased deamination; `CalculateContamination` flags cross-sample
contamination and downweights its evidence; the matched-normal filter
rejects any call with a non-trivial allele fraction in the normal.
The output VCF is conservative, well-annotated, and clinically usable.

*Strelka2* (Kim et al., *Nature Methods* 2018) and the older
*VarScan2* are the main alternatives, each with somewhat different
trade-offs. Strelka2 dominates the speed-per-sample metric and tends
to win on indel sensitivity in benchmark studies; VarScan2 has the
longest history of clinical use and is still common in legacy
pipelines. *LoFreq*, with its strong noise modelling, handles
tumour-only and low-VAF cases competitively. The choice among them is
usually driven by the validation work the calling lab has already
done, not by abstract benchmark performance.

The hardest somatic problem is *low-VAF subclonal detection*. At
Illumina's baseline per-base error rate near $0.1 %$, a $1 %$ VAF
call is at the edge of the noise floor; at $0.1 %$ VAF it is
indistinguishable from sequencer noise without extra machinery. The
extra machinery comes in a few forms. *Depth scaling* — increasing
coverage from $100 times$ to $500 times$ to $10,000 times$ — buys
sensitivity at low VAF by averaging out the per-read noise. *UMIs*
(unique molecular identifiers) tag each original DNA molecule before
amplification, so PCR duplicates can be collapsed into a consensus
that is more accurate than any of its members. *Duplex sequencing*
requires both strands of a single original molecule to agree, which
drops the error floor to $tilde 10^(-7)$ per base at the cost of
deeper coverage still. ctDNA assays use UMIs and duplex chemistry
routinely; the same techniques have started spreading back into
solid-tumour subclonal work.

#tip[
  Clinical actionability has its own VAF cutoffs that are distinct
  from statistical detectability. A variant at VAF $>= 5 %$ is
  reported and acted on; a variant at VAF between $1 %$ and $5 %$ is
  flagged as suspicious and confirmed with an orthogonal assay before
  driving therapy; a variant below $1 %$ is typically not reported
  outside of MRD contexts. The thresholds are not statistical —
  they reflect the operating point of the molecular tumour board.
]

=== Mutational Signatures as NMF

Somatic mutations carry the fingerprint of the process that caused
them. Ultraviolet photodimers preferentially produce C-to-T transitions
at TCN trinucleotides. Tobacco smoke produces C-to-A at specific
flanking contexts. APOBEC deaminases produce C-to-T and C-to-G at TCW
sites. Mismatch-repair deficiency produces an indel-heavy spectrum at
microsatellites and a particular SNV pattern. Spontaneous deamination
of methylated cytosine produces C-to-T at CpG sites and accumulates
clock-like with age.

The framework that made these patterns quantitative is Alexandrov,
Nik-Zainal, Wedge, and Stratton's 2013 *Nature* paper. Every
single-base substitution can be classified by its substitution type
(six options: C-to-A, C-to-G, C-to-T, T-to-A, T-to-C, T-to-G,
collapsed across strand) and its immediate flanking bases (four
options each side), producing $6 times 4 times 4 = 96$
*trinucleotide contexts*. Every tumour's SNV catalogue can be reduced
to a 96-dimensional count vector — its *mutation spectrum*. Distinct
mutagenic processes produce distinct 96-dimensional signatures, and a
tumour's spectrum is a non-negative linear combination of those
signatures weighted by their *exposure* in that tumour.

#figure(
  image("../../diagrams/lecture-18/05-cosmic-signatures.svg", width: 100%),
  caption: [
    Six canonical COSMIC SBS signatures shown as 96-bin spectra.
    SBS1 (clock-like CpG deamination), SBS3 (HRD-associated), SBS4
    (tobacco), SBS6 (mismatch-repair deficiency), SBS7a (ultraviolet),
    and SBS13 (APOBEC C-to-G at TCW). A tumour's mutation spectrum
    is a non-negative weighted sum of signatures like these.
  ],
) <fig:signatures>

The decomposition is the linear-algebra interlude of the chapter. Let
$bold(V) in bb(R)^(96 times N)$ be the observed mutation spectra of
$N$ tumours, stored as columns. Let $bold(W) in bb(R)^(96 times K)$
hold the $K$ signatures, also as columns, with non-negative entries.
Let $bold(H) in bb(R)^(K times N)$ hold the per-tumour exposures of
each signature, also non-negative. The model is

$ bold(V) approx bold(W) bold(H) $

with $bold(W) >= 0$, $bold(H) >= 0$. This is *non-negative matrix
factorisation* (Lee and Seung, *Nature* 1999), the same problem that
appears in hyperspectral unmixing, speech-source separation, and
document-topic modelling. The cancer-genomics literature rediscovered
NMF for this purpose in 2013 and has built a substantial software
ecosystem on top of it: *SignatureAnalyzer* (Broad, on a Bayesian NMF
backbone), *SigProfilerExtractor* (Alexandrov lab; the COSMIC catalogue
curator), *sigminer* in R, and the *mmsig* fitting approach for using
pre-computed signatures on small cohorts.

#note[
  Choice of $K$ — the number of signatures to extract — is the
  classical model-order-selection problem in disguise. Standard
  approaches use a stability-based criterion: factorise repeatedly
  for each candidate $K$ with different random starts, measure the
  cosine similarity of the resulting signatures across runs, and
  pick the largest $K$ for which the signatures remain stable. The
  COSMIC catalogue has $tilde 80$ SBS signatures as of v3.4; not
  all are present in any given cohort.
]

The COSMIC catalogue assigns each signature an etiology where one is
known: SBS1 is spontaneous deamination of 5-methylcytosine at CpG and
accumulates with age in every tumour; SBS2 and SBS13 are APOBEC; SBS3
is HRD-associated and a clinical biomarker for PARP-inhibitor
sensitivity; SBS4 is tobacco; SBS6, SBS15, SBS20, and SBS26 are MMR
deficiency; SBS7a-d are ultraviolet; SBS10a and SBS10b come from POLE
polymerase mutation. The interpretation is straightforward: decompose
the tumour, look at which signatures are present, and read the
mutagenic history.

#figure(
  image("../figures/ch18/f1-signature-nmf.svg", width: 95%),
  caption: [
    Mutational-signature decomposition as NMF. A tumour's 96-bin
    spectrum (left) is approximated as a non-negative weighted sum of
    signatures from a fixed catalogue (middle). The non-negative
    exposure vector (right) reads off the mutagenic processes active
    in the tumour. The same linear-algebra move solves speech source
    separation and document-topic modelling.
  ],
) <fig:nmf>

Signature-based actionability is now clinical reality. SBS3 (HRD) is
the genomic basis for PARP inhibitor (olaparib, niraparib, rucaparib)
sensitivity in ovarian, breast, pancreatic, and prostate cancers. The
MMR-deficient signature combination (SBS6 plus 15 plus 20 plus 26)
flags eligibility for immune checkpoint inhibitor therapy
(pembrolizumab is FDA-approved tumour-agnostic for MSI-H since 2017).
APOBEC signatures are an active research area for therapeutic
vulnerability. Ultraviolet signatures inform melanoma biology but do
not currently change therapy. Signatures have crossed the
research-to-clinic line, and the AMP/ASCO/CAP framework in
Section 1.6 now incorporates them.


== Structural Variants and Gene Fusions <sec:sv-fusions>

Cancer breaks chromosomes in ways no other field has to deal with at
scale. The short-read alignment framework from Chapter 2 produces
useful evidence — discordant pairs, split reads, soft-clipped tails —
but cannot resolve the resulting topology unambiguously. Long reads
do that better; RNA-seq adds an orthogonal view of the resulting
transcripts; the modern cancer SV pipeline integrates all three.

The catalogue of cancer-specific SVs is rich enough to deserve its own
vocabulary. *Chromothripsis* (Stephens et al., *Cell* 2011) is the
catastrophic shattering of one chromosome into dozens of fragments
that reassemble in random order, producing a localised cluster of
SVs that look like a single event when reconstructed correctly. It is
common in osteosarcoma and Li-Fraumeni-associated tumours.
*Chromoplexy* is a multi-chromosome rearrangement common in prostate
cancer. *Chromoanasynthesis* is chromosomal amplification with
variable copy-number. *Extracellular circular DNA* (ecDNA), recognised
as a major mode of oncogene amplification around 2017, carries MYC or
EGFR copies on circular extrachromosomal elements that segregate
randomly during mitosis and reach 50 to 100 copies per cell — a
mechanism for drug resistance and rapid phenotypic adaptation.

#figure(
  image("../../diagrams/lecture-18/06-sv-detection.svg", width: 100%),
  caption: [
    Short reads versus long reads on a complex chromothripsis event.
    Short reads catch the breakpoint signals through discordant
    pairs and split reads but reconstruct the topology only with
    heavy inference; long reads span the breakpoints directly and
    resolve the rearrangement unambiguously.
  ],
) <fig:sv>

SV detection from short reads relies on *Manta* (Chen et al., 2016),
*GRIDSS* (Cameron et al., 2017), and the older *DELLY*. Manta is the
Illumina default and handles deletions, small insertions, and tandem
duplications well; GRIDSS is more sensitive on balanced rearrangements
(inversions, translocations) and is the basis for the LINX
structural-annotation tool used in several clinical pipelines.
Short-read sensitivity on cancer SVs hovers around 60 to 80 % depending
on the SV class. Long-read tools — *Sniffles2* (Smolka et al., 2024) as
the field standard; *Severus* as a cancer-specific successor;
*CuteSV* — push that to 95 % or higher, particularly on the complex
events that short reads cannot reach.

=== Gene Fusions from RNA-seq

A *gene fusion* joins two genes into a single transcript through a
structural rearrangement. The fusion protein has properties neither
parent possessed alone — typically a kinase domain placed under a
constitutively-active promoter or an oligomerisation partner that
forces dimerisation and downstream signalling. Fusions are the
canonical targeted-therapy substrate because the resulting protein is
absent from normal tissue and depends entirely on the rearrangement
for its function. A drug that binds the fusion's active site has
near-zero off-target effect, which is the precision-oncology dream.

The historical first was *BCR-ABL1*, the Philadelphia chromosome
translocation $t(9;22)$ in chronic myeloid leukaemia, and the first
matched drug *imatinib* (FDA approval 2001) inaugurated the modern
era of targeted cancer therapy. *EML4-ALK* in lung adenocarcinoma
(roughly $5 %$ of NSCLC, treatable with crizotinib, alectinib,
lorlatinib), *TMPRSS2-ERG* in prostate cancer (common but not yet
directly targetable), *NTRK* fusions across many cancers (treatable
tumour-agnostic with larotrectinib), and *FGFR* fusions in
cholangiocarcinoma and bladder cancer round out the textbook examples.

#figure(
  image("../../diagrams/lecture-18/07-gene-fusion.svg", width: 95%),
  caption: [
    A gene fusion in three steps. Two parent genes on different
    chromosomes; a chromosomal rearrangement at a breakpoint inside
    intronic regions; a fusion transcript and protein that combine
    domains from the two parents. EML4-ALK in lung adenocarcinoma is
    the canonical example: the ALK kinase domain placed under
    EML4 promoter control becomes a driver.
  ],
) <fig:fusion>

Detection methods follow two complementary routes. *RNA-seq fusion
callers* (STAR-Fusion, *Arriba*, FusionCatcher, JAFFA) look for
chimeric reads that span the junction between two genes' transcripts;
they catch only fusions that produce abundant mRNA, but they catch
those at high specificity. *DNA-based SV callers* detect the
underlying breakpoint; they are more sensitive overall but include
many breakpoints that do not produce a functional transcript. The
modern clinical pipeline runs both and integrates. Arriba (Uhrig et
al., *Genome Research* 2021) is currently the standard clinical
RNA-seq caller because it has a curated cancer-fusion allowlist and
separates high-confidence calls from speculative ones.

Copy-number variation is the third leg of cancer structural change.
Whole-arm losses, whole-chromosome gains, focal amplifications to ten
to a hundred copies, and whole-genome doublings are all common.
Callers include *CNVkit*, *Sequenza*, *PureCN*, *ASCAT*, *TitanCNA*,
and *Facets*; each estimates tumour purity, ploidy, and per-segment
absolute copy number jointly. Clinically critical CNVs include
*HER2* (ERBB2) amplification in breast cancer, which drives
trastuzumab therapy and was the first precision-oncology CNV
biomarker; *MYC* amplification, often ecDNA-resident; *CDKN2A* loss,
which disables a tumour-suppressor pathway across many cancers; and
*EGFR* amplification in glioblastoma and lung cancer.

#warn[
  CNV calling on FFPE samples is famously noisy because the per-segment
  variance is dominated by the input-material variability rather than
  the sequencing. Fresh-frozen samples give cleaner calls; FFPE
  samples need higher coverage and tighter quality controls.
  Trying to call subclonal CNV from a $100 times$ FFPE sample
  produces calls that no clinical assay would report.
]


== Tumour Heterogeneity and Clonal Evolution <sec:heterogeneity>

Every late-stage tumour is a population of subclones. A single bulk
biopsy contains DNA contributed by each subclone at a fraction
proportional to that subclone's prevalence, plus normal DNA from
stroma and immune infiltrate. Recovering the subclones from the bulk
mixture is the *subclonal deconvolution* problem, and it is —
structurally — a blind source separation problem of the kind any EE
student has met in independent-component-analysis or
mixture-model contexts.

The setup is as follows. Let a tumour contain $K$ subclones with
unknown fractions $f_1, ..., f_K$ summing (with the normal fraction)
to one. Each variant either belongs to one subclone or to a common
ancestor of several. The observed VAF of a variant is a linear
combination of the contributions from the subclones that carry it,
weighted by their fractions and modified by tumour purity, ploidy,
and the local copy-number state. The *cancer cell fraction* (CCF) is
the variant's fraction expressed as a fraction of tumour cells rather
than total cells — VAF adjusted for purity, ploidy, and CN.
Reconstruction means clustering variants by their CCF, identifying
the clusters as subclones, and (in the multi-sample case) inferring
the phylogenetic tree relating them.

#figure(
  image("../../diagrams/lecture-18/08-subclonal-tree.svg", width: 100%),
  caption: [
    Subclonal reconstruction in three steps. The empirical
    distribution of VAFs across hundreds of variants (left) shows
    multiple modes; clustering on CCF (middle) identifies the
    subclones; the phylogenetic tree (right) places the subclones
    on a tree rooted at the most recent common ancestor, with trunk
    mutations shared by all leaves.
  ],
) <fig:subclone>

The dominant tools are *PyClone* (Roth et al., *Nature Methods* 2014;
Dirichlet-process clustering with explicit purity and CN priors),
*SciClone* (Miller et al., *PLOS Computational Biology* 2014; variational
inference), *PhyloWGS* (Deshwar et al., 2015; joint clustering and
phylogeny), and the more recent *CliP* for high-throughput clustering.
PyClone-VI is the maintained current version. Their inputs are a
variant set with VAF, CN, and purity annotations; their outputs are
clusters of variants at shared CCFs, optionally with a tree relating
them when multiple samples (multi-region or longitudinal) are available.

The phylogenetic step uses methods that resemble the species-tree
work in Chapter 12 of this book applied to a within-patient sample.
*LICHeE*, *MEDICC2* (which works on copy-number profiles rather than
SNVs), *SPRUCE*, and *CONIPHER* build trees relating subclones from
multi-sample data. The trunk of the tree carries mutations shared by
every cell — early drivers, the foundational events of the tumour's
lineage. Branches carry later mutations, including the resistance
variants that emerge under treatment selection and define the
recurrence subclone. Trunk drivers are the best therapy targets;
branch drivers are the source of treatment resistance.

#note[
  Subclonal deconvolution is identifiability-limited from a single
  bulk sample. Two distinct subclone configurations can produce
  identical bulk VAF distributions, and only multi-region biopsies or
  longitudinal samples break the ambiguity. The PCAWG project's
  pan-cancer subclonal-reconstruction effort (Dentro et al., *Cell*
  2021) is the largest systematic application; it found that a
  majority of tumours show evidence of subclonal structure but that
  reliable phylogeny needs multi-sample data.
]

Single-cell sequencing (Chapter 8) avoids the deconvolution problem
entirely by reading individual cells directly. *scRNA-seq* (10×
Chromium, Smart-seq) reads cell-by-cell expression and reveals the
tumour-microenvironment composition — T cells, macrophages,
fibroblasts, endothelial cells — that matters for immunotherapy
response. *scDNA-seq* (10× Chromium CNV) profiles per-cell copy
number and resolves CNV heterogeneity that bulk averages out.
*scATAC-seq* (Chapter 10) reveals epigenetic heterogeneity. Spatial
transcriptomics — Visium, Xenium, MERFISH — adds geography to the
cell-type catalogue and is starting to enter clinical research
pipelines.

#figure(
  image("../figures/ch18/f2-vaf-mixture.svg", width: 95%),
  caption: [
    Why VAF is a mixture observation. The same biological subclone
    produces different VAF distributions at three tumour purities
    and two copy-number states. Without an independent purity and
    CN estimate, the VAF distribution alone is ambiguous — a small
    minor mode could be a small subclone in a pure sample or a major
    subclone in an impure sample.
  ],
) <fig:vaf-mix>


== Clinical Interpretation: AMP/ASCO/CAP and the Biomarker Trio <sec:clinical>

Chapter 17 walked the ACMG/AMP framework for *germline* variant
classification: a probability scale from benign to pathogenic,
arrived at by aggregating evidence codes. Cancer somatic variants use
a different framework — the *AMP/ASCO/CAP* joint guidelines
(Li, Datto, Duncavage, et al., *Journal of Molecular Diagnostics* 2017)
— and the difference is essential to understand.

ACMG/AMP asks "what is the probability this variant is causing
disease?" — a hypothesis-testing question. AMP/ASCO/CAP asks "what
should I do about this variant?" — a decision-theoretic question. The
output is not a pathogenicity probability but a *tier* that maps to a
recommended action. The two systems coexist because they answer
different questions, and the same clinical laboratory typically
maintains both — one for hereditary risk, one for somatic findings.

The four tiers are as follows. *Tier I* (variants of strong clinical
significance) is FDA-approved-therapy territory: an EGFR L858R in
NSCLC maps to osimertinib, a BRAF V600E in melanoma to dabrafenib
plus trametinib, a HER2 amplification in breast cancer to
trastuzumab. *Tier II* (potential clinical significance) covers
variants with FDA approval in a different cancer type ("off-label"
in the regulatory sense) or evidence from clinical trials or
strong preclinical data. *Tier III* (unknown clinical significance)
sits in a cancer gene but has no current therapeutic implication.
*Tier IV* (benign or likely benign) catches confirmed germline,
known polymorphism, or otherwise non-actionable findings.

#figure(
  image("../../diagrams/lecture-18/09-amp-asco-cap.svg", width: 100%),
  caption: [
    The AMP/ASCO/CAP four-tier system. Each tier is action-oriented:
    Tier I prescribes a specific FDA-approved drug; Tier II points to
    clinical trials or off-label use; Tier III holds variants in
    cancer genes without current therapeutic implications;
    Tier IV catches benign or passenger findings.
  ],
) <fig:tiers>

Each tier has level-of-evidence subdivisions. Tier I-A is FDA-approved
therapy in this cancer type for this variant; Tier I-B is a
professional-guideline (NCCN) recommendation; Tier II-C is
FDA-approved in a different cancer; Tier II-D is clinical-trial or
case-report evidence. The four-by-letter structure produces a
fine-grained action mapping that the clinical reporter can defend
in a tumour board.

The knowledge bases that drive tier assignment are *OncoKB* (Memorial
Sloan Kettering; FDA-recognised for variant annotation), *CIViC*
(crowd-curated, open source; Griffith et al., *Nature Genetics* 2017),
and the newer *ClinGen Somatic Cancer Variant Curation*. A clinical
variant scientist queries OncoKB and CIViC, checks FDA labels for the
matching drug-gene-cancer triple, reviews the NCCN guideline, and
applies the AMP/ASCO/CAP tier with level of evidence. Several
commercial systems (VarSome, QIAGEN QCI) automate the lookup;
proprietary clinical-lab aggregators reproduce most of OncoKB and
CIViC internally with custom curation layers on top.

#note[
  The contrast with ACMG/AMP germline (Chapter 17) is worth
  internalising. Germline aggregates many independent evidence codes
  (population frequency, computational prediction, segregation,
  functional studies) into a pathogenicity score on a Bayes-factor
  scale. Somatic, by contrast, is largely a deterministic lookup
  keyed by (variant, cancer type, drug evidence) into a curated
  knowledge base. The two frameworks operate on different
  computational structures because they answer fundamentally
  different questions.
]

=== TMB, MSI, and HRD — the Biomarker Trio

Beyond single-variant tiers, three aggregate biomarkers drive
precision-oncology decisions. Each is a statistic over the somatic
catalogue rather than a single variant, and each connects to a
specific drug class.

*Tumour mutational burden* (TMB) is the number of non-synonymous
somatic mutations per megabase of coding sequence. A high TMB implies
many neoantigens — mutated protein fragments that the immune system
can recognise as foreign — and predicts response to immune checkpoint
inhibitors (pembrolizumab, nivolumab). The FDA approved a tumour-agnostic
pembrolizumab indication in 2020 for $"TMB" >= 10$ mutations per Mb
on the FoundationOne CDx assay.

*Microsatellite instability* (MSI) is the consequence of
mismatch-repair deficiency, whether from germline Lynch syndrome
(MLH1, MSH2, MSH6, PMS2 loss) or somatic MMR inactivation. MSI-H
tumours have characteristic short-repeat length variability,
characteristic mutational signatures (the SBS6, SBS15, SBS20, SBS26
quartet), elevated TMB, and strong response to immune checkpoint
inhibitors. The FDA's tumour-agnostic pembrolizumab approval for
MSI-H / dMMR tumours (2017) was a landmark — the first time a single
biomarker rather than a single-organ indication drove an approval.

*Homologous recombination deficiency* (HRD) is loss of double-strand
break repair through homologous recombination, most often from
BRCA1, BRCA2, or PALB2 loss but also from other HR-pathway disruption.
HRD-positive tumours show the SBS3 mutational signature, characteristic
indel patterns (ID6), large-scale genomic scarring (loss of heterozygosity,
telomeric allelic imbalance, large state transitions), and exquisite
sensitivity to PARP inhibitors (olaparib, niraparib, rucaparib).
PARP inhibitors are FDA-approved in ovarian, breast, pancreatic, and
prostate cancers contingent on HRD status.

#figure(
  image("../../diagrams/lecture-18/10-biomarker-trio.svg", width: 100%),
  caption: [
    The three precision-oncology biomarkers beyond single variants.
    TMB is a mutation rate; MSI is an instability metric at short
    repeats; HRD is a multi-feature signature combining SBS3, indel
    patterns, and large-scale genomic scarring. Each has an
    assay-specific threshold and each connects to a distinct drug
    class.
  ],
) <fig:ch18-biomarkers>

#warn[
  TMB, MSI, and HRD numbers are *assay-calibrated*, not universal.
  The FDA-approved $"TMB" >= 10$ mut/Mb threshold for tumour-agnostic
  pembrolizumab is specific to FoundationOne CDx; the same TMB
  numeric value from MSK-IMPACT or from whole-exome sequencing means
  something different because the assays sample different gene sets
  with different selection bias. Borrowing thresholds across assays
  is a recurring source of clinical confusion. Every clinical lab
  must validate its own threshold against its own assay's empirical
  distribution.
]

The biomarker trio has a coherent signal-processing structure.
Each is a detection statistic over the somatic catalogue. TMB is an
event-rate estimator. MSI is a noise-level measurement in a specific
frequency band (the short-repeat compartment of the genome). HRD is a
multi-feature detection statistic with a trained decision threshold.
Each has an operating point set by clinical-trial data and an
assay-specific calibration that matters for deployment in exactly
the same way detector thresholds in EE systems are calibrated per
hardware. The decision-theoretic framing is uniform across the trio,
and the trio plus single-variant tiers covers nearly all of modern
precision oncology's actionable surface.


== From Sequencing to Therapy <sec:therapy>

The payoff for the chapter's machinery is a decision in a tumour
board. The workflow follows a recognisable shape: a stage-IV patient
is biopsied; the FFPE block goes to a clinical lab; sequencing on a
tumour-normal or tumour-only panel produces a somatic VCF plus
fusion calls plus copy-number profile plus TMB, MSI, and HRD; a
variant scientist applies AMP/ASCO/CAP and writes a report with
Tier I to IV findings; the molecular tumour board — oncologist,
pathologist, variant scientist, researcher, clinical-trial-unit
representative — reviews the report; a therapy decision is made.

#figure(
  image("../../diagrams/lecture-18/11-therapy-tree.svg", width: 100%),
  caption: [
    The precision-oncology decision tree. Tier I variants drive
    matched targeted therapy; failing that, the TMB / MSI / HRD
    biomarker trio routes the patient to immune checkpoint inhibitor
    or PARP-inhibitor therapy; failing that, standard-of-care
    chemotherapy is the fallback.
  ],
) <fig:therapy>

Targeted therapy decisions follow the *drug-gene-cancer triple*
catalogue: EGFR L858R or exon-19 deletion in NSCLC pairs with
erlotinib, gefitinib, or osimertinib; EML4-ALK in NSCLC with crizotinib,
alectinib, or lorlatinib; BRAF V600E in melanoma with vemurafenib or
the dabrafenib-plus-trametinib combination; HER2 amplification in
breast cancer with trastuzumab; BRCA1/2 loss in ovarian, breast,
pancreatic, or prostate cancer with olaparib or niraparib; KRAS
G12C in NSCLC with sotorasib (the 2021 FDA approval that broke the
long-standing "undruggable" KRAS reputation); BCR-ABL1 in CML with
imatinib, dasatinib, or nilotinib (the original precision-oncology
success story).

#note[
  *Keys fit locks.* The unifying mechanism of every targeted-therapy
  triple is that the drug is designed to bind a specific protein
  shape produced by a specific driver alteration. Imatinib wedges
  into the ATP pocket of BCR-ABL1's constitutively-active kinase
  domain; when the fusion is absent the drug has nothing to bind
  productively; when it is present the drug shuts down the
  downstream signal. This is why targeted therapy is more selective
  than chemotherapy — narrow efficacy, narrow toxicity — and why
  the genome is the determinant of eligibility. Without the
  rearrangement, no key; without the key, no therapy.
]

When no Tier I variant is present, the biomarker trio carries the
load. High TMB or MSI-H routes the patient to pembrolizumab. HRD
routes to PARP inhibitor. PD-L1 expression measured by
immunohistochemistry, complementary to the genomic biomarkers,
informs the choice among checkpoint-inhibitor options. Failing those,
standard-of-care chemotherapy is the fallback, and a Tier II clinical
trial may be sought for any near-actionable findings.

Emerging therapy classes extend the framework rather than replacing
it. *Antibody-drug conjugates* (ADCs) — trastuzumab-deruxtecan for
HER2-low breast cancer is the canonical recent example — pair an
antibody that targets a tumour antigen with a cytotoxic payload that
is delivered once the antibody internalises. *Bispecific antibodies*
(blinatumomab in B-ALL, teclistamab in multiple myeloma) physically
bridge T cells to tumour cells. *CAR-T therapy* engineers a patient's
own T cells to express a tumour-targeting receptor (CD19 for
B-cell malignancies; an expanding catalogue for solid tumours).
*Therapeutic neoantigen vaccines*, often built on the mRNA platforms
that the COVID-19 vaccines popularised, use the patient-specific
somatic catalogue from sequencing to design a tailored vaccine.
Each of these expands the integration surface between genomics and
oncology, and each pulls more of the work the chapter has described
into the clinic.

Drug discovery itself is increasingly structure-based, in the lineage
of AlphaFold (Chapter 15). The 2021 FDA approval of sotorasib, the
first KRAS G12C inhibitor, illustrates the pattern: structure-based
design identified a cryptic binding pocket unique to the G12C
mutant form, enabling selective targeting after decades of KRAS
being called undruggable. AlphaFold-era tools accelerate this style of
discovery by making the input structure cheap.

#figure(
  image("../figures/ch18/f3-tumour-board-flow.svg", width: 95%),
  caption: [
    The tumour-board integration. Sequencing produces a per-patient
    variant package; the variant scientist applies AMP/ASCO/CAP and
    writes a tiered report; the molecular tumour board integrates
    the report with pathology, imaging, and patient-factor inputs
    and produces a therapy decision. The pipeline this book has
    been building feeds one input into a multi-disciplinary
    decision process.
  ],
) <fig:tumour-board>

=== Resistance and Combinations

Single-drug targeted therapy fails over time. The mechanism is the
clonal evolution of Section 1.5 in action: a resistance variant exists
at low VAF in the pre-treatment tumour, the targeted drug removes the
sensitive cells, the resistance subclone expands. The textbook example
is osimertinib resistance in EGFR-mutant NSCLC through EGFR T790M
emergence; longitudinal ctDNA monitoring catches the resistance
mutation before clinical progression and triggers a treatment switch.

Modern strategies push against resistance in three ways. *Upfront
combinations* — dabrafenib plus trametinib for BRAF-mutant melanoma,
encorafenib plus binimetinib — hit two pathway nodes simultaneously
so that no single resistance mutation rescues the cell. *Sequential
therapy* uses the first-line drug until resistance, then switches to
a mechanistically distinct second-line. *Adaptive trials* use ctDNA
monitoring to switch therapy at the molecular signal of resistance
rather than waiting for radiographic progression.

#figure(
  image("../../diagrams/lecture-18/12-industry-ecosystem.svg", width: 100%),
  caption: [
    The cancer-genomics industry ecosystem. Four overlapping
    quadrants — clinical laboratories, pharma and biotech,
    sequencing platforms, and academic cancer centres — and the
    bioinformatics roles inside each. Hybrid clinical, research,
    and ML career paths are increasingly the norm for new graduates.
  ],
) <fig:ecosystem>

#matters[
  The integration is the point. Cancer genomics is not a sequence of
  isolated computational problems — alignment in Chapter 2, variant
  calling in Chapter 4, RNA in Chapter 5, long reads in Chapter 11,
  ML in Chapter 16, clinical genomics in Chapter 17 — but a single
  workflow in which every prior chapter's tools meet a patient. The
  same pipeline that produces a FASTQ for a research cohort feeds a
  tumour board's decision about whether to give pembrolizumab or
  olaparib. The same NMF that recovers topics from a document corpus
  recovers tobacco and ultraviolet signatures from a melanoma
  spectrum. The same blind source separation that handles
  multi-speaker audio reconstructs subclones from bulk DNA. The
  technical depth this book has accumulated lands here, and the
  decisions it informs are the kind no other applied field rivals
  in stakes.
]


== Summary <sec:ch18-summary>

- Cancer is applied evolution: each tumour is a population of cells
  evolving under selection inside a single body, with VAF playing the
  role of allele frequency and drivers playing the role of selected
  alleles. The population-genetic framework from Chapter 12 maps onto
  cancer genomics with very little modification.
- Sequencing design is driven by the clinical question. Tumour-normal
  is the gold standard; tumour-only is common in clinical reality;
  FFPE artefacts (strand-biased C-to-T deamination) need
  orientation-bias filtering; ctDNA enables non-invasive longitudinal
  monitoring at the cost of $30,000 times$ coverage and UMI-based
  noise control.
- Somatic calling discards the germline three-state assumption.
  Mutect2 (and Strelka2, VarScan2) inherit local-assembly from
  HaplotypeCaller and add somatic-specific filters: Panel of Normals,
  orientation bias, contamination, tumour-in-normal. Low-VAF
  subclonal detection requires depth, UMIs, or duplex chemistry.
- Mutational signatures are non-negative matrix factorisation of the
  tumour mutation spectrum. Six substitution types times sixteen
  trinucleotide contexts give a 96-dimensional spectrum; signatures
  are basis vectors; exposures are non-negative coefficients. The
  COSMIC catalogue ($tilde 80$ SBS signatures as of v3.4) carries
  etiology labels for most, and HRD (SBS3) and MMR (SBS6 plus 15
  plus 20 plus 26) are now clinically actionable.
- Structural variants and gene fusions break short reads and need
  long reads or RNA-seq to resolve. BCR-ABL1 (CML, imatinib),
  EML4-ALK (NSCLC, crizotinib), NTRK fusions (tumour-agnostic,
  larotrectinib), and HER2 amplification (breast, trastuzumab) are
  the textbook fusion / CNV examples.
- Subclonal reconstruction is blind source separation. PyClone,
  SciClone, PhyloWGS cluster variants by cancer cell fraction;
  multi-sample data breaks identifiability and enables phylogenetic
  reconstruction. Trunk drivers are the targeted-therapy candidates;
  branch drivers are the source of resistance.
- AMP/ASCO/CAP is the somatic counterpart to ACMG/AMP. Tier I is
  FDA-approved therapy; Tier II is potential / trial / off-label;
  Tier III is unknown; Tier IV is benign. The framework is
  action-centric rather than probability-centric.
- TMB, MSI, and HRD are the three aggregate biomarkers beyond single
  variants. Each is assay-calibrated and connects to a drug class —
  TMB and MSI to checkpoint inhibitors, HRD to PARP inhibitors. Borrowing
  thresholds across assays is a recurring source of clinical error.
- The tumour board is the integrating venue. The pipeline this book
  has built feeds one input into a multi-disciplinary therapy
  decision. The chapter is the capstone of the sequencing,
  alignment, variant-calling, clinical-pipeline, and
  clinical-interpretation arc.


== Exercises <sec:ch18-exercises>

#strong[1.] #emph[VAF arithmetic.] A tumour sample is estimated to be
$45 %$ pure, with diploid background across the locus of interest. A
somatic variant is called at VAF $0.18$. (a) Compute the cancer cell
fraction assuming the variant is present at one copy per affected
tumour cell. (b) If the local copy number is three instead of two
(one allele duplicated), how does the implied CCF change? (c) The
caller reports a second variant at VAF $0.05$ in the same diploid
region. Is it consistent with a subclone present in $25 %$ of tumour
cells? Show the arithmetic.

#strong[2.] #emph[FFPE filter design.] You are designing an
orientation-bias filter for a tumour-only FFPE pipeline. C-to-T
calls account for $42 %$ of all SNV calls in your validation cohort
versus an expected $tilde 25 %$ from fresh-frozen reference data,
suggesting substantial deamination contamination. Propose a
quantitative filter: which per-read summaries do you compute, what
threshold do you set, and how do you validate the threshold without
a matched normal? Sketch the false-positive / false-negative
trade-off.

#strong[3.] #emph[NMF decomposition by hand.] A small toy tumour
has the following counts across four trinucleotide contexts:
T[C>T]G = 80, T[C>T]C = 40, C[C>A]A = 30, T[C>A]A = 20. Decompose
this 4-dimensional spectrum into a non-negative combination of two
"signatures": $bold(W)_1 = (0.5, 0.5, 0, 0)$ (a clock-like CpG
process) and $bold(W)_2 = (0, 0, 0.6, 0.4)$ (a tobacco-like process).
Compute the non-negative exposure weights $h_1, h_2$ that best fit
the observed spectrum under squared-error loss. Verify the
reconstruction.

#strong[4.] #emph[Fusion topology.] @fig:fusion shows EML4-ALK as the
canonical lung-adenocarcinoma fusion. The breakpoint is typically in
EML4 intron 13 and ALK intron 19. Sketch the fusion transcript at
the exon level. Which protein domains from EML4 and ALK survive in
the fusion? Why is the resulting kinase constitutively active even
though the wild-type ALK kinase is not? Cite at least one of the
canonical EML4-ALK papers.

#strong[5.] #emph[Biomarker calibration.] A clinical laboratory
running whole-exome sequencing on solid tumours reports a TMB of 12
mutations per Mb for a patient with NSCLC. Their reference panel
(FoundationOne CDx) reports the cutoff for tumour-agnostic
pembrolizumab as $"TMB" >= 10$ mut/Mb. Should this patient be
considered TMB-high? Discuss the assay-calibration issue. Propose
two empirical experiments the lab could run to translate its WES
TMB onto the panel-calibrated scale.

#strong[6.] #emph[Tier assignment.] A FoundationOne CDx report on a
metastatic colorectal-cancer patient lists: KRAS G12V (canonical
hotspot, no FDA-approved drug in colorectal at the time of writing);
BRAF V600E (canonical hotspot, FDA-approved encorafenib-cetuximab
combination); MLH1 deletion plus MSI-H by NGS; TMB of $32$ mut/Mb;
TP53 R273H (canonical tumour-suppressor hotspot, no direct
therapy). Assign AMP/ASCO/CAP tiers and levels of evidence to each
finding. Which biomarker drives the front-line therapy decision?

#strong[7.] #emph[Subclonal phylogeny.] A patient with metastatic
melanoma has three biopsies: the primary, a lymph-node metastasis,
and a brain metastasis. Each is sequenced; PyClone identifies four
variant clusters with CCFs reported in the following table.

#table(
  columns: 5,
  align: center,
  table.header[Cluster][Primary CCF][Lymph CCF][Brain CCF][Likely subclone],
  [A], [$1.00$], [$1.00$], [$1.00$], [trunk],
  [B], [$0.60$], [$0.10$], [$0.00$], [?],
  [C], [$0.05$], [$0.95$], [$0.40$], [?],
  [D], [$0.00$], [$0.05$], [$0.95$], [?],
)

Sketch a phylogenetic tree consistent with this CCF pattern. Place
B, C, D on branches and justify each placement. Which cluster is
the candidate brain-metastasis-specific driver and what would you
do with that information clinically?

#strong[8.] #emph[(Open-ended.)] Pick one cancer type with a
publicly available molecular sub-typing study — for example, the
TCGA Pan-Cancer paper on lung adenocarcinoma, or the
Pan-Cancer Atlas of Whole Genomes (Campbell et al., *Nature* 2020).
Identify one molecular subtype defined primarily by a structural
or copy-number alteration rather than a single SNV (chromothripsis,
ecDNA amplification, focal amplification, complex rearrangement).
Describe in one paragraph what the alteration is, how the original
study detected and validated it, and what targeted therapy or
clinical trial is currently aligned with it.


== Further Reading <sec:ch18-further-reading>

- #strong[Hanahan, D., and Weinberg, R. A.] (2011). "Hallmarks of
  cancer: the next generation." #emph[Cell] 144: 646–674. The most
  cited paper in cancer biology. Pair with Hanahan (2022),
  #emph[Cancer Discovery] 12: 31–46, for the four 2022 additions.
- #strong[Alexandrov, L. B., Nik-Zainal, S., Wedge, D. C., et al.]
  (2013). "Signatures of mutational processes in human cancer."
  #emph[Nature] 500: 415–421. The paper that introduced the
  NMF-based signature framework. Read with Alexandrov et al. (2020),
  #emph[Nature] 578: 94–101, for the modern catalogue.
- #strong[Cibulskis, K., Lawrence, M. S., Carter, S. L., et al.]
  (2013). "Sensitive detection of somatic point mutations in
  impure and heterogeneous cancer samples." #emph[Nature
  Biotechnology] 31: 213–219. The original MuTect paper. Benjamin et
  al. (2019), #emph[bioRxiv], extends to GATK 4 Mutect2.
- #strong[Roth, A., Khattra, J., Yap, D., et al.] (2014). "PyClone:
  statistical inference of clonal population structure in cancer."
  #emph[Nature Methods] 11: 396–398. The standard reference for
  subclonal deconvolution.
- #strong[Li, M. M., Datto, M., Duncavage, E. J., et al.] (2017).
  "Standards and guidelines for the interpretation and reporting of
  sequence variants in cancer: a joint consensus recommendation of
  the AMP, ASCO, and CAP." #emph[Journal of Molecular Diagnostics]
  19: 4–23. The AMP/ASCO/CAP tier system in full.
- #strong[Campbell, P. J., Getz, G., Korbel, J. O., et al.] (2020).
  "Pan-cancer analysis of whole genomes." #emph[Nature] 578:
  82–93. The PCAWG consortium overview. Companion papers cover
  signatures, structural variants, and subclonal architecture across
  2,800 WGS tumours.
- #strong[Stephens, P. J., Greenman, C. D., Fu, B., et al.] (2011).
  "Massive genomic rearrangement acquired in a single catastrophic
  event during cancer development." #emph[Cell] 144: 27–40. The
  paper that named chromothripsis.
- #strong[Turajlic, S., Sottoriva, A., Graham, T., and Swanton, C.]
  (2019). "Resolving genetic heterogeneity in cancer." #emph[Nature
  Reviews Genetics] 20: 404–416. A modern review of the subclonal
  evolution literature; useful pairing with PyClone for the
  methodological landscape.
