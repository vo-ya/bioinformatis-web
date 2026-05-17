#import "../theme/book-theme.typ": *

= Variant Calling: From Aligned Reads to Called Differences <ch:variant-calling>

#matters[
  Almost every clinical decision a modern genomics workflow informs —
  whether to treat with a particular tyrosine kinase inhibitor, whether
  to put a patient on enhanced cancer surveillance, whether to call a
  rare-disease case solved — rests on a single line of text in a single
  file. That line is a variant call. The pipeline that produces it
  starts with a tube of DNA on an instrument and ends with a string that
  reads `chr7	140753336	.	A	T	612	PASS …`. The four billion
  base-pairs in between have been compressed to a sparse list of
  disagreements with a reference, each annotated with a statistical
  argument for why the disagreement is real. Most analysis bugs in
  human genomics live somewhere in that compression. This chapter walks
  through every step.
]

Two people's genomes differ at roughly one base in a thousand. Multiplied
across the 3.1 gigabases of a human nuclear genome, that comes to three
to four million differences per individual relative to the GRCh38
reference — a few hundred thousand of them inside protein-coding genes,
a few thousand inside the canonical disease-associated set, perhaps a
single digit number with directly clinical consequence. The job of
variant calling is to find those millions of differences in a
billion-read dataset, decide which of them are real, and pass the survivors
on to a human or a software pipeline that turns them into clinical or
biological claims.

The four moves of variant calling are: stack the reads, call by
likelihood, filter honestly, annotate to action. This chapter expands
each of them. Sections 4.1 and 4.2 lay out what a variant is and how
the pileup organises the evidence. Sections 4.3 and 4.4 cover the
statistical machinery — Bayesian genotyping, the math of allele balance,
and how PHRED-encoded scores fall out of the posterior. Section 4.5
turns to the pre-calling pipeline (mark duplicates, INDEL realignment,
base-quality recalibration) and the three generations of callers, from
`samtools mpileup` (2008) through GATK HaplotypeCaller (2010) to
DeepVariant (2017). Section 4.6 takes the special case of somatic
variant calling in tumours, and Section 4.7 covers structural variants,
which need a fundamentally different family of algorithms. Section 4.8
returns to the output: VCF, filtering, and annotation against the
reference databases that turn calls into decisions.

The chapter assumes you have a coordinate-sorted, duplicate-marked BAM
in hand — that is, you understand the alignment step from Chapter 2.
Everything downstream of the VCF — population genetics, expression
analysis, structural biology — assumes you understand what is and is
not in the variant calls.


== What Is a Variant? <sec:what-is-a-variant>

A *variant* is a position in the genome where the sample sequence
differs from a reference sequence. A *mutation* is a change in DNA that
happened during replication, repair, or environmental damage — a
dynamic event in time. The two words are used interchangeably in
practice, but the distinction matters. A variant at the same position
in millions of people is still a variant, but hardly a mutation
anymore: it was a mutation once, in some unknown distant ancestor, and
has been faithfully copied through germlines for tens or hundreds of
thousands of years. A variant unique to one patient's tumour, by
contrast, is both. The variant is what you observe in a dataset; the
mutation is the biological process that created it.

Variants arise in three regimes:

- *Inherited (germline).* Present in the sample from the time of
  conception. Found in essentially every nucleated cell of the body, on
  one or both chromosomes. The substrate of inheritance and population
  genetics; the target of clinical-genetics laboratories.
- *Somatic.* Arose after conception in a subset of cells — most
  infamously in cancer, where a lineage of cells accumulates driver
  mutations. Variant-supporting reads are a _minority_ of reads at the
  position, sometimes a small minority.
- *De novo.* Present in a child but in neither parent — arose in
  the parental germline or in the embryo. A specific, important
  sub-class of germline variants, and the source of much of the
  high-penetrance rare-disease load.

A useful sanity check: at a normal human position the expected VAF
(variant allele frequency) under germline diploidy is 0, 0.5, or 1.0.
If you see 0.15 in a constitutional sample, the most likely
explanations are sample contamination, mosaicism, or a mapping
artifact. If you see 0.15 in a tumour sample at 30 % purity, it might
be a real somatic variant in half of the cancer cells. Same VAF, very
different inference.

=== Classification by Size and Structure

Variants form a continuous size spectrum, but operational pipelines
chop the spectrum into four bins because different tools work on
different bins.

- *Single-nucleotide variants (SNVs).* One-base substitution. A:T → G:C
  or any other flip. The most common variant class by a factor of ten.
- *Insertions and deletions (INDELs).* A few bases inserted or
  deleted, typically 1–50 bp. An order of magnitude rarer than SNVs but
  disproportionately impactful: a 1-bp deletion inside a coding region
  shifts the reading frame and usually knocks out the protein.
- *Structural variants (SVs).* Large changes, typically ≥ 50 bp:
  deletions, insertions, duplications, inversions, translocations.
  Rarer per individual than SNVs by another order of magnitude, but
  implicated in a large fraction of rare disease and cancer.
- *Copy-number variants (CNVs).* The sub-class of SV where the dosage of
  a large region changes (usually duplication or deletion). Often
  measured as a coverage shift rather than as discrete breakpoints.

#figure(
  image("../../diagrams/lecture-04/01-variant-taxonomy.svg", width: 95%),
  caption: [
    The variant taxonomy laid out along a log-scale size axis with
    typical per-genome counts. SNVs dominate by count; structural
    variants dominate in total affected bases.
  ],
) <fig:taxonomy>

The other axis people care about is _inheritance_: germline versus
somatic. Germline callers assume every variant sits at one of three
allele frequencies (0, 0.5, 1.0); somatic callers have to handle a
variant present in 5 % of reads because only 5 % of cells in the tumour
have it. Completely different statistical problem on the same pileup.

=== SNVs and Their Consequences

SNVs are the dominant class — about 3 million per human genome — but
the vast majority are benign. The clinical consequence of a particular
SNV depends almost entirely on where in the genome it falls.

If the SNV lies in a *non-coding region* (about 98 % of the genome) it
may do nothing, or it may alter a regulatory element — a promoter, an
enhancer, a splice site — with effects that range from undetectable to
severe. Most non-coding SNVs are treated as neutral by default; a
small minority, especially those at canonical splice positions, are
known to matter.

If the SNV lies in a *coding region* the consequence depends on how it
changes the codon. The four categories — synonymous, missense,
nonsense, splice-region — partition the space of coding SNVs.

- *Synonymous (silent).* The new codon still codes for the same amino
  acid. No protein change. The default assumption is neutrality;
  modern work has caught exceptions involving codon-usage effects on
  translation speed, RNA folding, and exonic splicing enhancers, but
  the assumption mostly holds.
- *Missense.* The new codon codes for a different amino acid. Impact
  ranges from undetectable to complete loss of function depending on
  which residue and how the chemistry shifts.
- *Nonsense.* The new codon is a stop codon. The protein is truncated
  — usually a loss-of-function allele, because truncated proteins are
  degraded or non-functional.
- *Splice-region.* The SNV disrupts a splice donor or acceptor site at
  an exon boundary. Consequences include exon skipping, intron
  retention, or a frameshift that ripples through the rest of the
  protein.

#figure(
  image("../../diagrams/lecture-04/02-snv-consequences.svg", width: 95%),
  caption: [
    Four consequence classes for a coding-region SNV. The codon context
    determines whether a single-base change is silent, changes an amino
    acid, creates a stop, or disrupts splicing.
  ],
) <fig:snv-consequences>

#note[
  An SNV is a single-bit error in a long genome string. What makes
  DNA different from a random bit stream is its non-uniform cost
  landscape — most positions tolerate substitutions, but a sparse set
  (active-site residues, splice boundaries, transcription-factor
  binding cores) have catastrophic costs. Variant _calling_ is the
  detection problem; variant _interpretation_ is reading out the cost
  landscape at each called position. The same error-found / error-cost
  split appears in convolutional decoding and in any error-correcting
  system with positional weighting.
]

=== INDELs and Frameshifts

INDELs are insertions and deletions of small numbers of bases. They are
roughly an order of magnitude rarer than SNVs but cause disproportionate
harm because of a single arithmetic fact: DNA is read in triplets.

An INDEL whose length is a multiple of three adds or removes whole
codons. The protein is one or a few amino acids longer or shorter, but
the reading frame downstream is preserved. Impact is comparable to a
missense variant: sometimes benign, sometimes severe, never categorically
worse than a substitution.

An INDEL whose length is _not_ a multiple of three shifts the reading
frame. Every codon downstream is now wrong. The translation machinery
reads random-looking nonsense until it hits a premature stop codon —
the protein is truncated, and the truncated transcript is usually
degraded by nonsense-mediated decay. This is a *frameshift* variant,
and it is almost always a complete loss-of-function of the affected
allele.

#figure(
  image("../../diagrams/lecture-04/03-indel-frameshift.svg", width: 95%),
  caption: [
    Frameshift versus in-frame INDEL. Same size range, radically
    different protein consequences, driven entirely by whether the
    length is a multiple of three.
  ],
) <fig:frameshift>

#warn[
  INDELs are harder to align and harder to call than SNVs. A single
  2-bp deletion near the end of a read is often misrepresented as a
  cluster of ten false SNVs if the aligner refuses to open a gap. This
  is why the pre-calling pipeline does INDEL realignment as a separate
  step (@sec:precalling) — and why most false-positive SNV calls in a
  typical pipeline are actually mishandled INDELs.
]

=== Why It Matters

The stakes for variant calling are not abstract. A pathogenic variant
in `BRCA1` shifts a 25-year-old woman's recommended cancer-surveillance
schedule from "screening at 50" to "annual MRI starting at 25." A
`SERPINA1` Z/Z genotype changes pulmonary care for life. A `BRAF` V600E
call in a melanoma biopsy decides whether the patient receives a BRAF
inhibitor or chemotherapy. Every week in every large hospital, a
radiologist's or oncologist's decision hinges on a variant call from a
pipeline like the one in this chapter.

#note[
  The transition from Sanger-era single-gene variant testing to
  genome-wide variant calling happened between 2008 and 2015 and was
  infrastructure-driven. GATK (Broad Institute, first release 2010) and
  `bcftools` (Heng Li, 2011) made calling 3 million variants feasible
  and reproducible. ClinVar (NCBI, 2012) made the clinical-interpretation
  layer shareable. By 2020, a clinical-grade germline VCF for a whole
  exome cost about \$200 to produce, with turnaround under 24 hours —
  cheaper than the office visit that returned the result.
]


== The Pileup: Aligned Reads, Column by Column <sec:pileup>

Before a caller runs, the information at each genomic position is
organised into a *pileup*: at every reference position, list every
read overlapping it together with the base that read carries at that
position and the PHRED quality score of that base. The pileup is the
caller's view of the world.

The `samtools mpileup` format is the canonical text representation.
A single line looks like this:

```
chr1    1000    A    8    ,,...,.,^]   IIIII!III
```

Six tab-separated columns, in order:

1. *Chromosome* — `chr1`.
2. *Position* — `1000`, 1-based.
3. *Reference base* at that position — `A`.
4. *Depth* — the number of reads overlapping the position. Here, 8.
5. *Read bases* — one character per read. `.` means "same as
   reference, on forward strand," `,` means "same as reference, on
   reverse strand," an uppercase letter (`ACGT`) is a mismatch on the
   forward strand, lowercase (`acgt`) is a mismatch on the reverse
   strand. `^]` marks the start of a read (with the second character
   encoding the read's mapping quality), `$` marks the end. `+3AAA` and
   `-2TT` encode small INDELs at the next position.
6. *Base qualities* — one ASCII character per read, PHRED+33 encoded
   (Chapter 1).

The grammar is dense but lossless: every column can be re-derived from
the underlying alignment. A column reading `,,.,,.,` shows seven
reference bases. A column reading `,,A,A,A,` shows seven reads, three
calling `A` and four calling reference — a candidate heterozygous SNV
at allele frequency 3/7.

#figure(
  image("../../diagrams/lecture-04/04-pileup-anatomy.svg", width: 95%),
  caption: [
    A `samtools mpileup` line, decomposed. The format is terse but
    every column can be re-derived from the underlying alignment.
  ],
) <fig:pileup-anatomy>

To make the move from pileup to call concrete, @fig:pileup-cases shows
three columns side by side: a clean heterozygous SNV, a clean
homozygous alternate, and a column of "reference plus a sequencing
error or two" that an inexperienced eye might mistake for a variant.
The same shape of column — a small count of bases with a small count
of qualities — admits three completely different inferences depending
on the distribution of the counts and the qualities.

#figure(
  image("../figures/ch04/f1-pileup-worked-example.svg", width: 95%),
  caption: [
    Three pileup columns at 20× coverage. Column A is a clean
    heterozygous site with balanced REF/ALT and uniformly high
    quality; column B is homozygous alternate; column C looks
    superficially like a low-frequency variant but the two
    non-reference bases are both at the end of low-quality reads on
    the same strand — the signature of noise, not signal.
  ],
) <fig:pileup-cases>

#note[
  A pileup is a histogram per genome column. For every column of the
  reference you have a small count of bases and a matching set of
  quality scores. Variant calling is the task of deciding, column by
  column, whether the column's distribution looks like "all reference"
  or like one of the non-reference genotypes — and how confidently.
]

=== The Vocabulary of a Variant Record

Every caller, and the VCF format itself, speaks the same vocabulary.
Learning it is half the battle.

- *REF* — the reference allele at the site. Usually one base for an
  SNV, longer for an INDEL.
- *ALT* — the alternate allele. A site can be multiallelic
  (`REF=A, ALT=C,G`).
- *DP* — depth of coverage at the site. The number of reads usable
  for calling, after duplicate-marking and mapping-quality filters.
- *AD* — allelic depths. Reads-per-allele as a comma-separated list;
  `AD=12,5` means 12 reads support REF and 5 support ALT.
- *VAF* or *AF* — variant allele frequency. ALT-reads divided by
  total reads. For germline diploids, expected VAFs are very close to
  0, 0.5, or 1.0.
- *GT* — genotype. For diploids, written as `0/0`, `0/1`, `1/1`, or
  `./.`. The digits index into `[REF, ALT1, ALT2, …]`. `0/1` is
  heterozygous, `1/1` is homozygous alternate, `1/2` is compound
  heterozygous with two different ALT alleles, and `./.` is missing.
- *GQ* — genotype quality, PHRED-scaled confidence that the called
  genotype is correct.
- *QUAL* — site-level variant quality. PHRED-scaled probability that
  there is _any_ variant at the site at all.
- *FILTER* — a string label summarising whether the caller thinks the
  call is real. `PASS` means the call survived all filters; anything
  else is a named failure (`LowQual`, `HighDP`, `StrandBias`, …).

#figure(
  image("../../diagrams/lecture-04/05-variant-record-anatomy.svg", width: 95%),
  caption: [
    One position, from reads to VCF line. Each field in the variant
    record is a specific statistic over the pileup column.
  ],
) <fig:record-anatomy>

#warn[
  VCFs are 1-based and inclusive on both ends. BED files are 0-based
  and half-open. `samtools` internally uses both. A variant reported
  in a VCF at `chr1:1000` is the same position as `chr1:999-1000` in
  a BED file. Silently mixing the two produces off-by-one errors in
  every downstream analysis — one of the most common real-world bugs
  in genomic pipelines.
]

=== The Ideal Case and Why It Doesn't Exist

If everything worked perfectly, variant calling would be trivial.
Imagine a 100× sample with zero sequencing errors and perfect
alignment. At a true heterozygous SNV you would see 50 reads carrying
the reference base and 50 reads carrying the variant base. At a true
homozygous site 100 reads would all agree. A simple threshold — call a
variant if the ALT count is at least some number — would work
perfectly.

Real data looks nothing like that. Seven things break the ideal:

1. *Sequencing errors.* Illumina reads carry per-base error rates of
   0.1–1 %. At 100× coverage, a true reference site shows 0–2 spurious
   non-reference bases on average — enough to create false positives
   under naive thresholding.
2. *Non-uniform coverage.* Some positions sit at 10×, some at 150×.
   Thresholds that work at 30× either miss variants at 10× or over-call
   at 150×.
3. *Mapping errors.* Reads from repeats, paralogs, or pseudogenes land
   on the wrong homologous region and show up as fake variants.
4. *Nearby INDELs.* If a nearby indel is not realigned, it gets
   expressed as a cluster of fake SNVs — a well-known source of noise.
5. *Strand bias.* True variants typically show approximately equal
   support from forward and reverse strands. A variant called entirely
   from one strand is almost always an artifact of library prep or PCR.
6. *Systematic base-quality errors.* Certain sequence contexts have
   elevated error rates that the machine's quality scores underestimate.
   Uncalibrated quality scores over-confidently endorse real errors.
7. *Low allele frequency in somatic calling.* A tumour at 30 % purity
   harbouring a clonal variant shows the variant at VAF 0.15 — well
   within noise for a germline-style threshold.

#tip[
  A back-of-envelope number worth knowing: at 30× coverage and 1 %
  per-base error, the expected number of erroneous reads at a given
  position is $30 times 0.01 = 0.3$. The Poisson tail at three or more
  errors is about 0.5 %. Multiply by $3 times 10^9$ positions and you
  get $1.5 times 10^7$ false-positive variant sites — five orders of
  magnitude beyond the truth. Calling cannot be a threshold; it must
  be a model.
]


== The Math: Bayesian Genotyping <sec:bayes>

Generation-two variant callers — `bcftools`, `GATK UnifiedGenotyper`,
and most modern tools inside their assembled-haplotype windows — treat
calling as posterior estimation. For each candidate position they
compute

$ P(G | D) prop P(D | G) P(G) $

where $G$ ranges over the three possible diploid genotypes
($0/0$, $0/1$, $1/1$), $D$ is the observed pileup (bases and
qualities), and $P(G)$ is a prior over genotypes given population
allele frequencies. The caller emits the maximum-posterior genotype as
its call, and the second-best posterior as its uncertainty.

The likelihood $P(D | G)$ factors across independent reads. For a
single read carrying base $b$ at the position with PHRED quality $Q$
and corresponding per-base error probability $epsilon = 10^(-Q/10)$,
the contribution to the likelihood depends on the genotype hypothesis:

- Under $G = 0/0$, the truth is REF on both chromosomes. A read showing
  ALT must be an error: probability $epsilon$. A read showing REF is
  correct: probability $1 - epsilon$.
- Under $G = 1/1$, mirror image: a read showing REF is an error
  ($epsilon$), a read showing ALT is correct ($1 - epsilon$).
- Under $G = 0/1$, the read picked one chromosome at random. With
  probability $1/2$ it sampled the REF chromosome; with probability
  $1/2$ it sampled the ALT chromosome. Then it either reported that
  allele correctly ($1 - epsilon$) or in error ($epsilon$). The
  combined probability of observing REF is
  $ 0.5 (1 - epsilon) + 0.5 epsilon = 0.5 $
  and likewise for ALT.

For $n_R$ reads supporting REF and $n_A$ reads supporting ALT at
quality $Q$ each, the per-genotype likelihoods are

$ P(D | 0/0) = (1 - epsilon)^(n_R) epsilon^(n_A) $
$ P(D | 1/1) = epsilon^(n_R) (1 - epsilon)^(n_A) $
$ P(D | 0/1) = (1/2)^(n_R + n_A) $

In log space, each read contributes a fixed quantum to the
log-likelihood: a read agreeing with a "homozygous for what I saw"
hypothesis adds $log(1 - epsilon)$, a disagreeing read adds
$log(epsilon)$, and under the heterozygous hypothesis every read adds
$log(1/2)$ regardless of which allele it carried.

#figure(
  image("../../diagrams/lecture-04/07-bayesian-genotyping.svg", width: 95%),
  caption: [
    Bayesian genotype calling at one position. Posterior equals
    likelihood times prior, taken over three candidate genotypes; the
    maximum-posterior genotype becomes the call.
  ],
) <fig:bayes>

To make the numbers concrete, suppose a pileup shows eight reads
supporting REF and four reads supporting ALT, all at Q30 (error rate
$epsilon = 10^(-3)$):

$ log P(D | 0/0) = 8 log(0.999) + 4 log(0.001) approx -27.6 $
$ log P(D | 0/1) = 12 log(0.5) approx -8.3 $
$ log P(D | 1/1) = 8 log(0.001) + 4 log(0.999) approx -55.3 $

The heterozygous hypothesis dominates by more than 19 nats over both
homozygous alternatives — a factor of $e^19 approx 1.8 times 10^8$ in
likelihood ratio. Multiply each by a prior — for a random human
position $P(0/0) approx 0.999$, $P(0/1) approx 7 times 10^(-4)$,
$P(1/1) approx 2.5 times 10^(-4)$ — and the prior tilts the homozygous
reference up by about three orders of magnitude. The het hypothesis
still wins comfortably; the call is $0/1$.

#figure(
  image("../figures/ch04/f2-genotype-posterior-math.svg", width: 95%),
  caption: [
    A worked Bayesian genotype call from a 12-read pileup (8 REF, 4
    ALT, all Q30). Log-likelihoods, prior, posterior, and the
    resulting PHRED-scaled genotype quality fall out of the same
    sum of per-read evidence.
  ],
) <fig:posterior-math>

The PHRED-scaled *genotype quality* (GQ) you see in a VCF is the
posterior probability of the second-best genotype, transformed by

$ "GQ" = -10 log_(10) P(G != "MAP")  =  -10 log_(10) (1 - P("MAP" | D)) $

A GQ of 30 says the caller assigns at most $10^(-3)$ probability to any
genotype other than the one it called. GQ of 99 — the default ceiling
for most callers — says the called genotype is at least
$10^(-9.9)$-confident in posterior terms, which is more confidence than
any caller's noise model can actually distinguish from $10^(-12)$;
GQ ≥ 99 is saturated.

#note[
  Bayesian variant calling is posterior estimation with a discrete
  hypothesis set. The three genotype hypotheses are mutually exclusive;
  the likelihood factors across independent-read evidence; the prior is
  learned from population data. This is the same structure as optimal
  detection with a finite-alphabet hypothesis set — MAP decoding of
  convolutional codes, digital-modulation symbol decoding, and
  classifier-combination in ensemble learning. The PHRED score is
  literally the log-likelihood-ratio contribution of each read; Q30
  means "1 in 1000" at the likelihood-ratio level.
]

=== Allele Balance: When the Het Likelihood Lies

The per-read likelihood for a het assumes the read samples either
chromosome with probability $1/2$ independently. Real Illumina data
breaks the assumption at low coverage: at 10× a true het can split
$2/8$ or $8/2$ by chance alone with probability about 9 %, indistinguishable
at low depth from a noisy hom-ref. Worse, reference-bias in alignment
makes the ALT chromosome marginally less likely to map than the REF
chromosome, so true hets at certain regions skew systematically toward
$"AD" = (n_R, n_A)$ with $n_R > n_A$.

Modern callers fold in an *allele-balance* term that penalises the
heterozygous likelihood when the observed REF/ALT ratio is too far
from 0.5. A common form is a beta-binomial likelihood with sequencing
error baked in:

$ L_("het", "AB") prop "Beta-Binomial"(n_A; n_R + n_A, alpha, beta) $

with $alpha, beta$ chosen so the mean is $0.5$ and the variance reflects
real cluster-level skew (typical fits give $alpha = beta approx 50$
for Illumina). The effect is to make $1/9$ or $9/1$ pileups much less
het-favoring than the naive binomial would suggest, while keeping
balanced $4/6$ or $5/5$ pileups firmly heterozygous.

=== From Likelihoods to PHRED-Scaled PLs

VCFs report not the posterior but the *PHRED-scaled likelihoods* (PL
field) for each genotype, normalised so the MAP genotype is 0:

$ "PL"(G) = -10 log_(10) ( P(D | G) / max_(G') P(D | G') ) $

A PL triplet `0,33,99` says the heterozygous likelihood is 1, the
hom-ref likelihood is $10^(-3.3)$ of het, and the hom-alt likelihood
is $10^(-9.9)$ of het. The reader of a VCF reconstructs the
posterior by multiplying through the prior. This separation —
likelihood in PL, prior left to the consumer — is what makes joint
calling across cohorts work: the per-sample likelihoods can be combined
across samples without rerunning per-sample inference.


== Haplotype Callers and the Local Re-Assembly Trick <sec:haplotype>

The per-position model in @sec:bayes treats every column as
independent. That assumption breaks down whenever two variants are
close enough that a single read sees both: the two variants are no
longer independent observations of two columns, they are joint
observations of a haplotype. INDELs make this worse — a single INDEL
can change the read-to-reference alignment for tens of bases around
it, scrambling the per-column analysis.

The fix is *haplotype calling*: in a sliding window around any
candidate variant, locally re-assemble the reads into a small number
of candidate haplotypes, then ask which two haplotypes (with their
specific combination of variants) best explain the reads. GATK
HaplotypeCaller (Broad, 2010) and DeepVariant (Google, 2017) both work
this way, although they differ in how they score haplotypes.

The procedure inside HaplotypeCaller is:

1. *Define an active region* — a stretch of genome where at least
   one read suggests a possible variant (mismatches, soft-clips,
   indels).
2. *Build a De Bruijn graph* (Chapter 3) from the reads in the active
   region, enumerate non-redundant paths, and trim weakly-supported
   branches. The result is a small set of candidate haplotypes —
   typically two to six.
3. *Realign every read against every haplotype* with a pair-HMM,
   producing per-read per-haplotype likelihoods that account for indels
   and gap penalties without the per-position model's blind spots.
4. *Marginalise* over haplotype pairs to recover per-genotype
   likelihoods at every variant inside the region, and emit them in
   the VCF.

The result is a caller that handles INDELs and clustered variants
robustly without needing a separate INDEL-realignment pass. (The
old `GATK IndelRealigner` step was deprecated in GATK4 for exactly
this reason.)

#figure(
  image("../../diagrams/lecture-04/06-indel-realignment.svg", width: 95%),
  caption: [
    INDEL misalignment and the local-realignment fix. A single true
    INDEL manifests as a cluster of false SNVs under naive
    per-position alignment; haplotype-aware calling converges the
    reads onto a single gap.
  ],
) <fig:realignment>

=== DeepVariant: The CNN Replaces the Likelihood

DeepVariant (Poplin et al., 2018) dispenses with the explicit
statistical model. The pipeline is:

1. Run a candidate-finding pass over the BAM to identify positions
   with any non-reference signal.
2. At each candidate position, convert a window of reads into a small
   RGB image. Columns are genome positions, rows are reads, pixel
   channels encode base identity, base quality, mapping quality,
   strand, and mismatch-to-reference.
3. Feed the image through a convolutional neural network trained on
   millions of benchmarked variant examples — initially Inception v3,
   now a custom architecture.
4. Output: a three-way classification (`hom-ref`, `het`, `hom-alt`)
   with a confidence score.

The training data is curated and standardised: Genome in a Bottle's
HG001–HG007 samples each come with an orthogonal-truth VCF derived
from combining many sequencing technologies and assemblies, and the
network learns context-specific features that a Bayesian model would
have to encode explicitly (homopolymer error context, strand-bias
patterns, repeat-region mapping artifacts).

DeepVariant beats every Bayesian caller on benchmark F1 score by a
small but consistent margin, at the cost of roughly 2–3× runtime and a
GPU dependency. The shift from likelihood-based to CNN-based variant
calling tracks a wider pattern in ML-heavy bioinformatics: a
well-defined inference task with millions of labelled examples becomes
a supervised-learning problem the moment someone builds the training
set. The Bayesian formulation is still useful for understanding what
the network is implicitly learning, and for calibrating on corner
cases the training distribution missed.

#note[
  The history is worth one paragraph. `samtools mpileup` ships in
  Heng Li's 2008 SAMtools release — the first widely used likelihood
  variant caller, with a per-position binomial model. `bcftools call`
  inherits and refines the approach. GATK UnifiedGenotyper (2010) adds
  a more principled prior and joint cohort calling. GATK
  HaplotypeCaller (2014) introduces local re-assembly. DeepVariant
  arrives at the end of 2017 and immediately tops the precisionFDA
  Truth Challenge leaderboard. Every modern caller is one of these
  three families.
]


== The Pre-Calling Pipeline <sec:precalling>

Before a caller runs, three pre-processing steps clean up the BAM.
Each removes a specific class of artifact, and each was discovered the
hard way — by chasing false positives back through a pipeline until
the source revealed itself.

*1. Mark duplicates.* PCR amplification during library preparation
creates multiple reads from the same original template fragment.
Duplicate reads carry identical sequence and identical errors;
counting them as independent evidence overstates the support for any
artifact they happen to share. `Picard MarkDuplicates` or `samtools
markdup` flags duplicate reads (defined as reads with identical
start position and orientation; pairs require identical mate position
too) so the caller can ignore them. Cost: a few CPU-minutes. Saves:
dozens to hundreds of spurious variants per genome, concentrated in
high-duplication regions.

*2. INDEL realignment.* A single INDEL near the end of a read is often
misaligned as many small substitutions — the aligner preferred to call
mismatches rather than open a gap because mismatches were locally
cheaper. Local realignment around known or candidate INDEL sites
revisits the region with a gap-aware aligner and produces a cleaner BAM.
Historically, GATK's `IndelRealigner` was a required step. Modern
callers (HaplotypeCaller, DeepVariant) do local re-assembly internally
and don't need a separate realignment pass. If you are running an
older pipeline or a simpler caller (`bcftools`), explicit INDEL
realignment still matters.

*3. Base-quality score recalibration (BQSR).* The sequencer assigns a
PHRED score to every base, but the raw scores are systematically off:
the machine uses a calibration learned at manufacturing time, and
sequence context (homopolymers, GC-extreme regions) plus sample-specific
chemistry produce error rates that deviate from the nominal curve.
BQSR re-estimates the quality scores empirically by looking at _known_
variant sites (dbSNP, known polymorphisms) and treating any mismatch
not at a known site as an error. Mismatches are bucketed by context
(read position, dinucleotide context, original quality), the empirical
error rate per bucket is computed, and a new quality string is written
back to the BAM.

#note[
  BQSR is systematic-error calibration of a measurement instrument.
  The sequencer is a noisy channel with a context-dependent
  systematic distortion. Calibration against known-good reference
  positions fits an error model per context, then corrects
  measurements at unknown positions. The same procedure exists in
  every serious measurement pipeline: dark-frame subtraction in
  astronomy, bias-current compensation in ADCs, NTP clock-drift
  correction in network timing.
]

BQSR matters most when the caller is Bayesian and uses quality scores
as likelihoods. If the caller treats Q30 as $10^(-3)$ error probability
but the real error rate in a particular sequence context is $10^(-2)$,
the caller under-weights real errors and over-calls variants in that
context. DeepVariant absorbs the context dependence into the CNN and
is therefore much less sensitive to whether BQSR ran — one of several
reasons it requires less pipeline orchestration.


== Filtering: Where Pipelines Silently Fail <sec:filtering>

A raw VCF from any caller contains both high-confidence calls and
noise. Filtering splits the two. The standard quality signals are:

- *QUAL* — site-level variant quality (@sec:pileup). Usually a
  per-caller empirical scale, not directly comparable across callers.
- *DP* — too low means too little evidence; too high usually means a
  mapping-error region. Both tails should be filtered.
- *GQ* — genotype quality. Low GQ means the caller could not cleanly
  decide between two genotypes.
- *FS* / *StrandOddsRatio (SOR)* — strand-bias metrics. A variant
  supported only by forward-strand reads is almost always an artifact.
- *MQ* / *MQRankSum* — mapping-quality metrics. Variants in
  low-mapping-quality regions are suspect.
- *ReadPosRankSum* — whether variant-supporting reads cluster near
  read ends (where quality is worst). Clustering signals an artifact.

Two filtering paradigms dominate practice:

- *Hard filtering.* Apply per-metric thresholds: `QUAL < 30` reject,
  `DP < 10` reject, and so on. Simple, interpretable, tunable. Used by
  `bcftools` pipelines and recommended whenever the cohort is too small
  to train a model.
- *VQSR or ML-based filtering.* Fit a Gaussian mixture model (GATK's
  VQSR) or a neural network to a labelled subset of variants from
  known-truth sites; use the model to classify all sites. Strictly
  better than hard filtering when there are enough training data;
  worse when there are not.

#warn[
  Over-filtering is the failure mode nobody talks about. A pipeline
  tuned to maximise precision by aggressive filtering throws away true
  variants in difficult contexts. A variant you filtered out does not
  raise an error — it produces a missed diagnosis. Always report both
  precision and recall on a held-out truth set; 99 % precision at 85 %
  recall is usually worse than 95 % precision at 98 % recall for a
  clinical pipeline.
]

The cascade of filters from raw VCF to actionable shortlist is large.
@fig:filter-cascade walks one canonical example: a whole-genome germline
VCF with about 4.8 million raw calls is whittled down through quality,
region, allele-frequency, and consequence filters to roughly 50
candidate variants suitable for clinical review.

#figure(
  image("../figures/ch04/f3-vcf-filter-cascade.svg", width: 95%),
  caption: [
    The filtering cascade for a typical whole-genome germline VCF.
    Each stage cuts the variant list by a factor of ten or more;
    the survivors at the bottom are a human-manageable shortlist
    for clinical review.
  ],
) <fig:filter-cascade>


== Somatic Calling: A Different Statistical Problem <sec:somatic>

Germline and somatic callers share the pipeline shape — aligned BAMs in,
VCF out — but solve different statistical problems.

A *germline caller* assumes a diploid genome. At every position the
true genotype is one of three (`0/0, 0/1, 1/1`), and the VAF of any
variant should be very close to 0, 0.5, or 1.0. Evidence of a variant
at VAF 0.15 means something is _wrong_ — probably contamination,
mosaicism, or mapping artifact. Germline callers include `GATK
HaplotypeCaller`, DeepVariant, `bcftools call`, and `Strelka2` (germline
mode).

A *somatic caller* expects a mixture. The sample is a tumour with
unknown purity (say 30 % tumour cells, 70 % infiltrating normal tissue)
harbouring a subset of variants each present in a subset of cancer
cells. A variant present in 50 % of tumour cells in a 30 %-pure sample
manifests at VAF $0.5 times 0.3 = 0.15$ — a real, important variant
hiding at the same apparent frequency as sequencing noise. The caller
has to distinguish them.

The standard approach is *matched tumour/normal calling.* Sequence
both the tumour and a matched normal tissue (usually peripheral blood)
from the same patient. At every candidate site, compare the tumour
pileup to the normal pileup: a true somatic variant is present in the
tumour and _absent from the normal_. Callers: `Mutect2` (GATK),
`Strelka2` somatic mode, `VarScan2`, DeepSomatic.

#figure(
  image("../figures/ch04/f4-somatic-vs-germline-vaf.svg", width: 95%),
  caption: [
    Expected VAF distributions under germline diploidy and tumour
    mixture models. Germline VAFs concentrate at 0, 0.5, 1.0; tumour
    VAFs spread continuously and depend on purity, ploidy, and
    sub-clonal structure.
  ],
) <fig:vaf>

#note[
  Somatic calling is signal detection with a prior subtraction. The
  "signal" is the tumour-specific variant allele; the "background" is
  the normal-tissue pileup carrying germline variants you do _not_
  want to call as somatic. Subtracting (or, more honestly, modelling
  jointly) the two distributions is the direct analogue of heterodyne
  detection or differential amplification — measure the difference,
  not the individual levels. The noise floor drops accordingly.
]

Somatic callers also have to handle *sub-clonal structure*. A single
tumour is often a mixture of sub-populations with different variant
complements. A variant in 80 % of cells is a clonal driver, one in 5 %
of cells is a late sub-clonal event. Both are real, and both carry
different implications — clonal variants inform what the tumour _is_,
sub-clonal variants inform resistance and heterogeneity. Sub-clonal
calling requires either very deep sequencing (≥ 500×) or a model that
explicitly enumerates sub-clones.


== Structural Variants <sec:sv>

SNVs and small INDELs are within-read events: a single 150 bp read can
carry the entire signal of the variant. A 5 bp deletion sits inside the
read's alignment as a 5 bp gap; the next read at the same position
reproduces the signal independently.

Structural variants span distances much larger than a read. A 5 kb
deletion cannot be seen inside a 150 bp read — no single read carries
the break. Instead the evidence appears as a pattern across _multiple_
reads: some read pairs land much further apart than expected, some
reads split across the breakpoint and half-map to each side, the depth
dips to zero inside the deleted region. A SV caller reconstructs the
event from this multi-read evidence. Different scales, different
algorithms.

=== The Three Signal Channels

Every SV caller exploits some combination of three signal channels.

*Discordant read pairs.* A paired-end Illumina library has a known
insert-size distribution — mean about 400 bp, standard deviation about
100 bp, for a typical short-fragment library. If a read pair lands with
5000 bp between its ends, the pair has probably straddled a $approx$
4600 bp deletion. The signal at any one pair is weak, but cluster
enough pairs around the same two positions and the deletion localises
to a few hundred bases.

The probability that any single pair from a 400 ± 100 bp library
straddles a 2 kb deletion in a heterozygous sample is roughly equal to
the probability that the fragment was drawn from across the deletion
breakpoints. For a depth of 30× and an average fragment length of 400
bp, the expected number of discordant pairs supporting a het deletion
is about $30 times (L_"del" + 2L_"frag") / (2 L_"frag")$ — a dozen
or two for a 2 kb event, falling to single digits for events near the
fragment length.

*Split reads.* A single read spanning a breakpoint aligns with its
first half matching one location and its second half matching a
distant location. The aligner reports both alignment segments (the
"primary" plus one or more "supplementary" alignments in the BAM);
the boundary between them is the breakpoint itself. Split reads give
single-base-resolution breakpoints — the gold standard for SV calls.
At 30× coverage with 150 bp reads, the expected number of split reads
supporting a het breakpoint is roughly the coverage scaled by the
fraction of read length that needs to overhang each side for the
aligner to commit. For a 20 bp minimum overhang on each side, the
expected count is $30 times (150 - 40) / 150 approx 22$ — comparable
to the discordant-pair count, but with much better localisation.

*Read depth.* A homozygous deletion drops coverage inside the deleted
region to zero. A heterozygous deletion halves it. A homozygous
duplication doubles it. Depth-based calling alone is good for large
CNVs (≥ 10 kb) but loses sensitivity at small events, where the depth
fluctuation lies inside normal coverage noise. Combined with discordant
pairs and split reads, it disambiguates tough cases.

#figure(
  image("../../diagrams/lecture-04/08-sv-detection-signals.svg", width: 95%),
  caption: [
    The three signal channels SV callers exploit, all reporting on the
    same deletion event. Individually each signal is weak; clustered
    and combined, they localise the event to single-base resolution.
  ],
) <fig:sv-signals>

=== The SV Types

Structural variants are classified by what they do to the genome, not
by their size (though size correlates). The six canonical types:

- *Deletion (DEL).* A contiguous region is absent from the sample.
  Can be heterozygous (one copy lost) or homozygous (both lost).
- *Insertion (INS).* A sequence is inserted at a position. Can be a
  novel sequence, a transposable-element insertion (SINE, LINE, SVA),
  or a tandem duplication of adjacent sequence.
- *Inversion (INV).* A segment is reversed in orientation — the
  sequence is the reverse complement of the reference in that region.
  No net loss or gain, but breakpoints can disrupt genes crossed at
  the boundary.
- *Duplication (DUP).* A region is present in multiple copies. Can be
  tandem (adjacent copies), dispersed (copies elsewhere), or segmental
  (large historically shared blocks).
- *Copy-number variant (CNV).* The specific case of a DEL or DUP
  expressed as a copy-number shift rather than as two discrete
  breakpoints — usually reported by depth-based callers operating on a
  single sample without a normal.
- *Translocation (TRA or BND).* A piece of one chromosome is joined to
  a different chromosome. Reciprocal translocations swap arms;
  non-reciprocal events insert one chromosome's material into another
  unidirectionally. Represented in VCFs as breakend ("BND") records
  with mate pointers, not as ordinary REF/ALT lines.

#figure(
  image("../../diagrams/lecture-04/09-sv-types.svg", width: 95%),
  caption: [
    The five canonical SV types, each altering the reference in a
    distinct structural way. Size and mechanism vary; the breakpoint
    bookkeeping varies accordingly.
  ],
) <fig:sv-types>

Real genomes routinely contain *complex SVs* — chains of primitive
operations: an inversion flanked by deletions, a duplication-inversion,
chromothripsis (a region shattered and re-assembled in a single
catastrophe). Modern SV callers recognise some of these; most do not,
and complex SVs are over-represented in the unsolved end of rare-disease
sequencing.

=== SV Callers and Their Tradeoffs

SV callers specialise by signal channel and by input data. The
dominant short-read options in 2024:

- *Manta* (Illumina, 2016). Read-pair plus split-read. Fast,
  conservative. The workhorse for germline SV calling on clinical
  whole-genome sequencing.
- *DELLY* (Tobias Rausch, 2012). Read-pair plus split-read. More
  sensitive than Manta, higher false-positive rate. Good for tumour
  work with aggressive filtering.
- *GRIDSS* (Daniel Cameron, 2017). Read-pair plus split-read plus
  local assembly. State-of-the-art for complex SVs but slower.
- *LUMPY* (Aaron Quinlan group, 2014). Multi-signal integration
  framework that takes pre-computed split-read and discordant-pair
  inputs and runs probabilistic combination.
- *CNVkit*, `cn.MOPS`, `cnvpytor`. Depth-based CNV calling for large
  (≥ 10 kb) copy-number events; complementary to the breakpoint callers.

Long reads change the picture. A single 20 kb PacBio HiFi read spans
most of an SV, so read-based detection becomes direct: if the read's
alignment shows a 5 kb gap, that is the deletion. Long-read-focused
callers — *Sniffles* (Sedlazeck et al., 2018), *pbsv* (PacBio),
*CuteSV* — routinely detect SVs that short-read callers miss,
especially events in repeat-rich regions.

#note[
  Split-read detection is a discontinuity-detection problem on a
  spatially-aligned signal. The read's alignment score-versus-position
  profile has a step at the breakpoint — the aligner either confidently
  maps through the break or it does not. The discontinuity localises
  the event to single-base resolution. The same problem appears in
  edge detection, change-point detection in time series, and glitch
  detection in telemetry; SV callers differ mostly in how aggressively
  they consolidate multiple weak discontinuities into a single strong
  event.
]

The current state of SV discovery looks like this. Short reads are the
workhorse for general clinical WGS because of cost — about \$500 per
genome versus \$1500–3000 per genome for HiFi — and they catch roughly
60–70 % of SVs accurately. Long reads are the gold standard for
research-grade SV catalogues and for any rare-disease case where
short-read analysis has stalled; they catch over 95 % of SVs, including
the hard ones in repeats. Ensemble methods — calling with multiple
tools and intersecting — beat any single caller in published
benchmarks. The first SV-comprehensive human dataset was the
Human Genome Structural Variation Consortium's 2019 release: three
family trios sequenced with short reads plus HiFi plus ONT plus
Strand-seq plus BioNano optical maps plus Hi-C. Processing produced
about 27,000 SVs per haplotype — four to five times what short-read-only
analyses had reported for the same samples. It was the moment the
field accepted that short-read-only SV catalogues are systematically
incomplete.


== The VCF File Format and Annotation <sec:vcf>

VCF (Variant Call Format) is the text format for variant-call output.
Every caller emits VCF; every downstream tool reads it. Understanding
its grammar is non-optional.

A VCF file has three sections:

- *Header.* Lines starting with `##` are meta-information: reference
  genome, caller version, contig lengths, INFO-field and FORMAT-field
  definitions, FILTER definitions. One line starting with `#CHROM`
  defines the column names for the data rows.
- *Data rows.* One per variant site. Fixed first 8 columns (CHROM,
  POS, ID, REF, ALT, QUAL, FILTER, INFO), then a FORMAT column listing
  the per-sample field grammar, then one column per sample.
- *Compression and indexing.* VCF is plain text; compressed to
  `.vcf.gz` with `bgzip` (block-gzip, indexable) and indexed with a
  `.tbi` (tabix) or `.csi` file for random access by genomic range.

A minimal single-sample example:

```
##fileformat=VCFv4.3
##reference=GRCh38
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Per-sample depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO            FORMAT       SAMPLE1
chr1    100234  .       A    G    420   PASS    DP=45;AF=0.51   GT:AD:DP:GQ  0/1:22,23:45:99
chr1    100501  rs12345 C    T    612   PASS    DP=52;AF=1.0    GT:AD:DP:GQ  1/1:0,52:52:99
```

Two rows: a heterozygous SNV at position 100,234 and a homozygous
alternate SNV at 100,501 (the second carries a dbSNP rsID in the ID
column). Every column has a defined meaning; every piece of software
reading the file relies on the header declarations to interpret the
INFO and FORMAT fields.

#figure(
  image("../../diagrams/lecture-04/10-vcf-anatomy.svg", width: 95%),
  caption: [
    The anatomy of a VCF file. The header declares the schema; data
    rows carry one variant per line; per-sample columns follow the
    FORMAT-string grammar.
  ],
) <fig:vcf>

#note[
  A VCF is a sparse differential encoding of a genome. Instead of
  transmitting 3 Gb of sequence per sample, you transmit only the
  $approx$ 4 M differences from a shared reference — a compression
  ratio around 750×. This is entropy coding at the file-format level:
  the reference is the codebook, the VCF is the per-sample delta.
  Indexable compressed VCFs (`bgzip` + `tabix`) are the same design
  as indexed compressed time-series databases — block-level
  compression for random access by region without decompressing the
  whole file.
]

A single VCF can hold many samples — a cohort, a trio, a tumour/normal
pair. Each sample contributes its own column with its own GT, AD, DP,
GQ. Sites are included if any sample has a non-reference call; samples
that don't call get `0/0` or `./.`.

The genotype field deserves one more piece of vocabulary. `0|1` (pipe,
phased) means the caller knows which allele sits on which of the two
parental chromosomes. `0/1` (slash, unphased) means the
allele-to-chromosome assignment is unknown. Phasing matters for
compound-heterozygous analysis and for trio inheritance: unphased data
fundamentally cannot distinguish two heterozygous variants in trans
(on different chromosomes) from in cis (on the same chromosome). Long
reads and parental trios are the two reliable paths to phasing in
practice.

=== Annotation: Turning Calls into Action

A VCF row that says "there is a variant at `chr7:140753336`, REF=A,
ALT=T" tells a clinician nothing on its own. *Annotation* is the step
where every variant is enriched with the information needed to
interpret it:

- *Gene and transcript context.* Which gene does this variant fall in?
  Which exon? Which codon? What amino-acid change?
- *Consequence classification.* Missense, nonsense, splice-region,
  intronic, upstream, intergenic.
- *Population frequency.* How common is the variant in gnomAD, in
  1000 Genomes, in population-specific sub-cohorts?
- *Clinical significance.* Is it in ClinVar? With what pathogenicity
  classification?
- *Functional predictions.* SIFT, PolyPhen, REVEL, CADD score
  estimates of how deleterious a missense variant is.
- *Splicing predictions.* SpliceAI scores for variants near splice
  sites.

Two annotators dominate the ecosystem: *Ensembl VEP* (Variant Effect
Predictor), the gold standard with multi-species support and the
slowest but most complete reading; and *snpEff / snpSift*, faster,
Java-based, widely used in research pipelines. Both read a VCF and
output an annotated VCF with appended INFO fields. A typical VEP
output adds a `CSQ` field with pipe-separated values, one per
transcript the variant intersects:

```
CSQ=T|missense_variant|MODERATE|BRAF|ENSG00000157764|Transcript|
ENST00000288602|protein_coding|15/18|...|p.Val600Glu|...|
REVEL=0.932|SpliceAI=0.01
```

That string says: the alternate allele T causes a missense variant of
moderate impact in BRAF, on transcript ENST00000288602, in exon 15 of
18, changing protein position 600 from valine to glutamate. The REVEL
score is 0.93 (highly deleterious); SpliceAI is 0.01 (no splice
effect). Layered on gnomAD ("not seen in 250,000 exomes") and ClinVar
("Pathogenic, drug-response, MULTIPLE submitters"), the variant becomes
clinically actionable: BRAF V600E, the canonical melanoma driver.

#figure(
  image("../../diagrams/lecture-04/11-annotation-pipeline.svg", width: 95%),
  caption: [
    The variant-annotation pipeline. Each stage joins the variant to
    a reference database and appends fields; the final output is
    clinically actionable.
  ],
) <fig:annotation>

#tip[
  Keep annotation separate from variant calling. A VCF called in 2020
  and annotated with 2024 databases gives a materially different
  clinical interpretation than one annotated with 2020-era databases —
  even though the called variants are identical. Annotation is the
  volatile layer of the pipeline; re-run it whenever the underlying
  databases update, especially ClinVar (refreshed monthly) and gnomAD
  (released yearly).
]


== Summary <sec:summary>

- Variant calling bridges alignment and interpretation. Aligned reads
  in, a filtered annotated VCF out. Every step between exists because
  real data breaks the idealised thresholding model.
- Calling is Bayesian, not thresholding. Per-position posterior genotype
  probabilities are computed from per-read evidence weighted by PHRED
  quality. DeepVariant replaces the explicit Bayesian layer with a CNN;
  the underlying inference structure is the same.
- The pre-calling pipeline is not optional. Mark duplicates, INDEL
  realignment (or haplotype-aware calling), base-quality recalibration.
  Each removes a specific class of artifact; skipping any one leaves a
  systematic bias the caller cannot recover from.
- Structural variants need different algorithms. Read pairs, split
  reads, and depth are the three signal channels. Manta, DELLY, GRIDSS,
  and the long-read callers (Sniffles, pbsv) serve different read
  regimes; ensembling wins on benchmarks.
- VCF is a sparse encoding. Understand the field grammar and you can
  reason about what any variant call means. Annotation is the volatile
  layer that turns calls into clinical interpretations.
- Filtering is where pipelines silently fail. Over-filtering removes
  true variants without errors. Always benchmark against a truth set;
  a 1 % recall drop on the benchmark predicts a 1 % miss rate on real
  samples.


== Exercises <sec:exercises>

#strong[1.] #emph[ASCII pileup parsing.]
Given the `samtools mpileup` line

```
chr1  10500  C  10  .,.,A,A,Aa  IIIIIIIIDB
```

decode the read bases and base qualities, then compute (a) the depth,
(b) the REF count, (c) the ALT count and identity, (d) the strand
distribution of the ALT-supporting reads, and (e) the candidate VAF.
Is this column more likely a real het or an artifact? Justify in one
sentence.

#strong[2.] #emph[Genotype likelihood by hand.]
Set up a pileup of 20 reads at one position: 14 supporting REF and 6
supporting ALT, all at Q30 ($epsilon = 10^(-3)$). Compute
$log P(D | 0/0)$, $log P(D | 0/1)$, and $log P(D | 1/1)$. Which
genotype wins, and by what natural-log likelihood ratio over the
runner-up? Convert the runner-up margin to a PHRED-scaled GQ.

#strong[3.] #emph[Allele balance.]
A heterozygous site in a 20× sample shows AD = (16, 4). Compute the
binomial probability of observing this split or one more extreme
under a fair 50/50 het. Would you trust this call? What additional
evidence would change your mind?

#strong[4.] #emph[Filter cascade arithmetic.]
A whole-genome germline VCF has 4.8 million raw calls. Filters applied
in order: PASS-only retains 88 %; remove ENCODE-blacklist regions
retains 95 % of those; gnomAD AF < 1 % retains 5 %; consequence is
missense / nonsense / frameshift retains 4 %; in ClinVar Pathogenic or
in a curated disease-gene panel retains 8 %. How many variants survive?
At what stage is the largest absolute drop?

#strong[5.] #emph[SV signal counts.]
For a heterozygous 2 kb deletion in a sample sequenced to 30× short-read
coverage with 150 bp reads and a fragment-size distribution of 400 ±
100 bp: estimate (a) the number of discordant read pairs straddling
the breakpoint, (b) the number of split reads supporting the
breakpoint, (c) the expected depth inside the deletion. Show the
arithmetic.

#strong[6.] #emph[Reading a multi-allelic VCF.]
The VCFv4.3 specification handles multi-allelic sites by listing
multiple ALT alleles on a single line. Consider a site at `chr3:5000`
with `REF=AT, ALT=A,ATT`. Describe in one paragraph what each ALT
represents (in terms of insertion / deletion / substitution), and how
a downstream tool should split this line into its constituent variants
for annotation.

#strong[7.] #emph[Germline-versus-somatic at the same VAF.]
A pileup at one site shows 60 reads supporting REF and 12 reads
supporting ALT (VAF = 0.17), all Q30. (a) Under a germline-diploid
prior, what is the maximum-posterior genotype? (b) Under a tumour
prior with 30 % purity and uniform sub-clonal architecture, would the
same pileup be called somatic? Explain how the prior structure
differs.

#strong[8.] #emph[(Open-ended.)]
Pick one tool from the chapter — `Mutect2`, DeepVariant, Manta,
Sniffles, or another — and read its primary publication. In one
paragraph, describe the single most surprising design choice the
authors made and why it works on the empirical data they show.


== Further Reading <sec:further-reading>

- *Li, H., et al.* (2009). "The Sequence Alignment / Map Format and
  SAMtools." _Bioinformatics_ 25: 2078–2079. The original `samtools
  mpileup` paper and the file-format spec the rest of the field is
  built on.
- *McKenna, A., Hanna, M., Banks, E., et al.* (2010). "The Genome
  Analysis Toolkit: A MapReduce Framework for Analyzing
  Next-Generation DNA Sequencing Data." _Genome Research_ 20:
  1297–1303. The GATK paper.
- *Li, H.* (2011). "A Statistical Framework for SNP Calling, Mutation
  Discovery, Association Mapping and Population Genetical Parameter
  Estimation from Sequencing Data." _Bioinformatics_ 27: 2987–2993.
  The `bcftools` Bayesian model in one tight paper.
- *Poplin, R., Chang, P.-C., Alexander, D., et al.* (2018). "A
  Universal SNP and Small-Indel Variant Caller Using Deep Neural
  Networks." _Nature Biotechnology_ 36: 983–987. DeepVariant.
- *Chen, X., Schulz-Trieglaff, O., Shaw, R., et al.* (2016). "Manta:
  Rapid Detection of Structural Variants and INDELs for Germline and
  Cancer Sequencing Applications." _Bioinformatics_ 32: 1220–1222.
- *Sedlazeck, F. J., Rescheneder, P., Smolka, M., et al.* (2018).
  "Accurate Detection of Complex Structural Variations Using
  Single-Molecule Sequencing." _Nature Methods_ 15: 461–468. Sniffles
  and the case for long reads in SV discovery.
- *VCFv4.3 specification.* `samtools.github.io/hts-specs/VCFv4.3.pdf`
  The file format. Read it before writing any VCF-processing code.
- *GATK Best Practices.* `gatk.broadinstitute.org` The canonical
  pipeline guide for short-read germline and somatic calling.
