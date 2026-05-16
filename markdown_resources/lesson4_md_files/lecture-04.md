# Lecture 4 — Variant Calling: From Aligned Reads to Called Differences

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-04.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Classify a DNA variant by its size and structure (SNV, INDEL, SV, CNV) and name the clinical/functional consequence of each.
2. Read a samtools pileup and describe what each column means.
3. Compute a genotype likelihood from a small pileup by hand, and explain why variant calling is Bayesian rather than thresholding.
4. Name the three pre-calling pipeline steps (mark duplicates, INDEL realignment, base-quality recalibration) and explain what each removes.
5. Pick the right caller (GATK, bcftools, DeepVariant, Strelka2, Manta, Delly) given a problem: germline SNVs, somatic variants, structural variants, or small-sample WES.
6. Read a VCF line and extract REF, ALT, QUAL, GT, AD, DP, GQ, VAF.
7. Explain why SVs need a different algorithm family from SNV/INDEL calling, and name the three signal channels SV callers exploit (read-pair, split-read, depth).
8. Annotate a VCF with functional consequences and identify the canonical public databases used at each annotation step.

---

## Part 1 — What is a Variant? (≈30 min)

### 1.1 Variants and mutations (≈7 min)

Two people's genomes differ. On average, at one base out of every thousand. Multiplied across a 3.1 Gb human genome, that's roughly three to four million differences per individual. Every clinical diagnostic genomics workflow and every population-genetics study in the last twenty years exists to find and interpret those differences.

A **variant** is a position in a genome where the sample sequence differs from a reference sequence. A **mutation** is a change in DNA sequence that happened during replication, repair, or environmental damage — a dynamic event in time. In practice the two words are used interchangeably, but the distinction matters: "variant" is what you observe in a dataset; "mutation" is the biological process that created it. A variant at the same position in millions of people is still a variant but hardly a mutation anymore; a variant unique to one patient's tumor is both.

Variants arise in three regimes:

- **Inherited (germline).** Present in the sample from the time of conception. Found in roughly every cell of the body. The substrate of inheritance and population genetics.
- **Somatic.** Arose after conception in a subset of cells — most infamously in cancer, where a lineage of cells accumulates driver mutations. Found in only a fraction of the sample's cells, which means the variant reads are a *minority* of reads at that position.
- **De novo.** Present in a child but in neither parent — arose in the parental germline or in the embryo. A specific, important subset of germline variants.

> **Intuition box**: A variant is a diff between two genome strings. Two genomes are 99.9% identical, and variant calling is the process of finding that last 0.1%. The inputs are aligned reads (from Lecture 2 or an assembly from Lecture 3); the output is a list of positions where the sample disagrees with the reference, with evidence attached.

### 1.2 Classification by size and structure (≈6 min)

Variants form a continuous size spectrum, but we chop that spectrum into four operational bins because different tools work on different bins.

- **Single-nucleotide variants (SNVs)** — one-base substitution. A:T → G:C, or any other flip. The most common variant class by a factor of ten.
- **Insertions / deletions (INDELs)** — a few bases inserted or deleted, typically 1–50 bp. Orders of magnitude rarer than SNVs but disproportionately impactful: a 1-bp deletion inside a coding region creates a frameshift that can knock a protein out entirely.
- **Structural variants (SVs)** — large changes, typically ≥50 bp: deletions, insertions, duplications, inversions, translocations. Rarer still per individual but implicated in a large fraction of rare disease and cancer.
- **Copy-number variants (CNVs)** — a sub-class of SV where the dosage of a large region changes (usually duplication or loss). Can be a single extra copy or many.

**FIGURE — Figure #1: Variant taxonomy and size spectrum** → `diagrams/lecture-04/01-variant-taxonomy.svg`
*A log-scale size axis with SNV, INDEL, SV, and CNV bins laid out along it; annotated with typical counts per human genome.*

The other axis people care about is **inheritance**: germline (present at birth, in every cell) vs somatic (acquired, subset of cells). Germline callers expect every variant to be at ~50% or ~100% allele frequency (one copy of two, or two of two). Somatic callers have to handle a variant present in 5% of reads because only 5% of cells in the tumor have it. Completely different statistical problem on the same pileup.

### 1.3 SNVs and their consequences (≈7 min)

SNVs are the dominant variant class — about 3 million per human genome, out of the 3–4 million total variants. The vast majority are neutral; a small minority are clinically or evolutionarily consequential. The consequence depends on where the SNV falls.

If the SNV is in a **non-coding region** (about 98% of the genome), it may do nothing, or it may alter a regulatory element (promoter, enhancer, splice site) with effects that range from subtle to dramatic. If the SNV is in a **coding region**, the effect depends on how it changes the codon at that position:

- **Synonymous / silent.** The new codon still codes for the same amino acid. No protein change. Historically treated as neutral; modern genetics has caught a few exceptions involving splicing or codon-usage effects, but the default assumption holds.
- **Missense.** The new codon codes for a different amino acid. Protein sequence changes. Impact ranges from "no detectable effect" to "complete loss of function," depending on which residue and how the chemistry shifts.
- **Nonsense.** The new codon is a stop codon. Protein is truncated — usually a loss-of-function allele, because truncated proteins are degraded or non-functional.
- **Splice-region.** The SNV disrupts a splice donor or acceptor site at an exon boundary. Can produce a frameshift, skip an exon, or introduce an intron — consequences ripple through the resulting protein.

**FIGURE — Figure #2: SNV consequence in coding regions** → `diagrams/lecture-04/02-snv-consequences.svg`
*A segment of mRNA with three adjacent codons; four panels below show the four SNV-consequence classes with the resulting codon and amino acid change.*

> **EE framing**: A SNV is a single-bit error in a long genome string. What makes DNA different from a random bit stream is the *non-uniform cost landscape* — most positions tolerate substitutions (synonymous SNVs, non-coding SNVs), but a sparse set of positions have catastrophic costs (active-site residues, splice boundaries). Variant *calling* is detection; variant *interpretation* is reading out the cost landscape at each called position. The same "error-found / error-cost" split applies to convolutional-code decoding and to any error-correcting system with positional weighting.

### 1.4 INDELs and frameshifts (≈6 min)

INDELs are insertions and deletions of small numbers of bases. They are much rarer than SNVs — roughly 1:10 ratio of INDELs:SNVs per genome — but they cause disproportionate harm because of a single arithmetic fact: **DNA is read in triplets.**

An INDEL of length that is a multiple of 3 adds or removes whole codons. The protein is one or a few amino acids longer or shorter, but the **reading frame** downstream is preserved. Impact is like a missense variant: sometimes benign, sometimes severe.

An INDEL of length that is *not* a multiple of 3 shifts the reading frame. Every codon downstream of the INDEL is now wrong. The translation machinery reads random-looking nonsense until it hits a premature stop codon — the protein is truncated, and the truncated transcript is usually degraded by nonsense-mediated decay. This is a **frameshift** variant, and it is almost always a complete loss-of-function of the affected allele.

**FIGURE — Figure #3: Frameshift vs in-frame INDEL** → `diagrams/lecture-04/03-indel-frameshift.svg`
*Two coding sequences, one with a 3-bp deletion (in-frame, one amino acid removed, rest intact) and one with a 1-bp deletion (frameshift, everything downstream changes).*

> **Warning box**: INDELs are harder to align and harder to call than SNVs. A single 2-bp deletion near the end of a read can be misrepresented as a cluster of ten false SNVs if the aligner doesn't open a gap. This is why the pre-calling pipeline does **INDEL realignment** (§3.1) as a separate step — and why most false-positive SNV calls in a typical pipeline are actually mishandled INDELs.

### 1.5 Why we care (≈4 min)

The stakes for variant calling are not abstract:

- **Clinical diagnostics.** A single pathogenic variant in BRCA1 changes cancer surveillance recommendations. A SERPINA1 Z/Z genotype changes pulmonary care. Every week a radiologist's or oncologist's decision hinges on a variant call from a pipeline like the one in this lecture.
- **Rare disease.** Whole-genome / whole-exome sequencing identifies causal variants for otherwise-undiagnosed conditions. The clinical labs running this in 2024 produce ~40% diagnostic yields.
- **Population genetics.** gnomAD's ~1 million exomes / ~100k genomes catalog lets every new variant be looked up for population frequency, which is the strongest single predictor of clinical significance.
- **Cancer genomics.** Somatic variant calling from tumor/normal pairs identifies driver mutations and informs targeted therapy choices (e.g., EGFR variants → specific TKI choices in NSCLC).

> **Historical pointer**: The move from Sanger-era single-gene variant testing to genome-wide variant calling happened between 2008 and 2015, and it was infrastructure-driven. GATK (Broad Institute, first release 2010) and bcftools (Heng Li, 2011) made calling 3 million variants feasible and reproducible. ClinVar (NCBI, 2012) made the clinical-interpretation layer sharable. By 2020 a full clinical-grade germline VCF for a whole exome cost ~$200 to produce, with turnaround under 24 hours.

---

## Part 2 — Variant Calling Mechanics (≈40 min)

### 2.1 What variant calling is, start to finish (≈8 min)

Variant calling is the bridge between a BAM file (reads aligned to a reference) and a VCF file (a list of positions where the sample differs from the reference, with evidence). In between sits the caller.

The canonical inputs:

- A **reference FASTA** (the reference genome sequence, e.g. GRCh38).
- An **aligned BAM** of reads from the sample, coordinate-sorted and indexed.
- Optionally, a **known-variants VCF** (dbSNP, ClinVar) for base-quality recalibration and joint calling.

The canonical output:

- A **VCF file** — one row per variant site per sample, with reference allele, alternate allele(s), per-sample genotype, quality scores, and caller-specific fields.

The canonical pipeline:

```
[reads] → align → [BAM] → mark duplicates → realign INDELs → recalibrate base qualities → [clean BAM]
        → caller → [raw VCF] → filter → [filtered VCF] → annotate → [final VCF]
```

Every stage matters. Skipping any one produces a VCF that looks normal until someone runs downstream analysis on it and the numbers stop making sense.

### 2.2 The pileup: aligned reads, column by column (≈8 min)

Before a caller runs, the information at each position is organised into a **pileup**: at every reference position, list every read that overlaps it, along with the base that read carries at that position and the quality score of that base.

The `samtools mpileup` format is the canonical text representation:

```
chr1    1000    A    8    ,,...,.,^]   IIIII!III
```

Column-by-column:

- **Chromosome** (`chr1`).
- **Position** (`1000`, 1-based).
- **Reference base** at that position (`A`).
- **Depth** — number of reads overlapping the position (`8`).
- **Read bases** — one character per read. `.` means "same as reference, on forward strand." `,` means "same as reference, on reverse strand." Upper-case letter (`ACGT`) means "mismatch, forward strand." Lower-case (`acgt`) means "mismatch, reverse strand." `^]` starts a new read; `$` ends one; `+3AAA` / `-2TT` encode INDELs at the next position.
- **Base qualities** — one ASCII character per read, Phred-encoded (§Lecture 1).

The compressed grammar is dense but lets you eyeball a variant: a column with `,,.,,.,` shows seven reference bases; a column with `,,A,A,A,` shows seven reads where three call `A` and four call reference. That's a candidate heterozygous SNV at frequency 3/7.

**FIGURE — Figure #4: Anatomy of a samtools pileup** → `diagrams/lecture-04/04-pileup-anatomy.svg`
*Six annotated columns of a pileup with callouts for each field, plus an inset showing what the `+N`, `-N`, `^]`, and `$` extensions encode.*

**EMBED — Artifact #1: Pileup Viewer** → `artifacts/lecture-04/01-pileup-viewer.html`
*Paste a small BAM-equivalent alignment (or use a preset); the viewer renders the pileup column-by-column and highlights candidate variant columns in real time as you scroll.*

> **Intuition box**: A pileup is a histogram per genome column. For every column of the reference, you have a small count of bases seen and a small count of quality scores attached. Variant calling is the task of deciding, column by column, whether the column's distribution of bases looks like "all reference" or "one of the non-reference genotypes."

### 2.3 The vocabulary of a variant record (≈8 min)

Every caller (and the VCF format itself) speaks the same vocabulary. Learning it is half the battle.

- **REF** — the reference allele. Usually one base for a SNV, longer for an INDEL.
- **ALT** — the alternate allele. Can be multiple if the site is multiallelic (e.g., `REF=A, ALT=C,G`).
- **DP** — depth of coverage at the site. Number of reads usable for calling.
- **AD** — allelic depths. Reads-per-allele as a comma-separated list. `AD=12,5` means 12 reads support REF, 5 support ALT.
- **VAF / AF** — variant allele frequency. `ALT_reads / total_reads` = 5/17 = 0.294 in the example above. For germline diploid samples, expected VAFs are approximately 0 (homozygous reference), 0.5 (heterozygous), or 1.0 (homozygous alternate).
- **GT** — genotype. For diploids, written as `0/0`, `0/1`, `1/1`, or `./.`. The digits are indices into `[REF, ALT1, ALT2, …]`. `0/0` = homozygous reference. `0/1` = heterozygous. `1/1` = homozygous alternate. `1/2` = compound heterozygous with two different ALT alleles. `./.` = missing, caller couldn't decide.
- **GQ** — genotype quality. Phred-scaled confidence that the called genotype is correct. `GQ=30` means probability of error ≤ 10⁻³.
- **QUAL** — site-level variant quality. Phred-scaled probability that there is *any* variant at this site.
- **FILTER** — a string label summarising whether the caller thinks this call is real. `PASS` means the call survived all filters; anything else is a named failure (e.g., `LowQual`, `HighDP`, `StrandBias`).

**FIGURE — Figure #5: Anatomy of a variant call at one position** → `diagrams/lecture-04/05-variant-record-anatomy.svg`
*A pileup column with 17 reads (12 A, 5 G); below it, the resulting fields populated: REF=A, ALT=G, DP=17, AD=12,5, VAF=0.294, GT=0/1, GQ=95, QUAL=420.*

> **Warning box**: VCFs are 1-based and inclusive on both ends. BED files are 0-based and half-open. samtools internally uses both. A variant reported in a VCF at `chr1:1000` is the same position as `chr1:999-1000` in a BED file. Silently mixing the two produces off-by-one errors in every downstream analysis — one of the most common real-world bugs in genomic pipelines.

### 2.4 The ideal case — and why it doesn't exist (≈8 min)

If everything worked perfectly, variant calling would be trivial.

Imagine a 100× coverage sample with zero sequencing errors and perfect alignment. At a true heterozygous SNV, you'd see 50 reads carrying the reference base and 50 reads carrying the variant base. At a true homozygous site, 100 reads all agree. At a true reference site, 100 reads all match reference. A threshold — "call a variant if ALT count ≥ some number" — would work perfectly.

Real data looks nothing like that. Seven things break the ideal:

1. **Sequencing errors.** Illumina reads have ~0.1–1% per-base error. At 100× coverage, a true reference site will show 0–2 spurious non-reference bases on average — enough to create false positives if you just threshold.
2. **Non-uniform coverage.** Some positions have 10×, some have 150×. Thresholds that work at 30× either miss variants at 10× or over-call at 150×.
3. **Mapping errors.** Reads from repeats, paralogs, or pseudogenes can land on the wrong homologous region and show up as fake variants.
4. **INDELs near the position of interest.** If a nearby indel isn't realigned (§3.1), it gets expressed as a cluster of fake SNVs — a well-known source of noise.
5. **Strand bias.** True variants typically show approximately equal support from forward and reverse strands. A variant called entirely from one strand is almost always an artifact of library prep or PCR.
6. **Systematic base-quality errors.** Certain sequence contexts have elevated error rates that the machine's quality-score assignment underestimates. Uncalibrated quality scores over-confidently endorse real errors.
7. **Low allele frequency (somatic calling).** A tumor at 30% purity harboring a variant in 50% of cancer cells shows the variant at VAF = 0.15 — well within noise for germline-style thresholds. Somatic callers have to distinguish 0.15-VAF true calls from noise at the same VAF.

Every one of those problems has a corresponding fix in modern pipelines. The next three subsections cover the fixes.

> **Discussion prompt**: Why not just call every position where ≥3 reads disagree with reference as a variant? (You'd get on the order of millions of false positives per sample. At 30× coverage and 1% per-base error, each position has roughly 30 × 0.01 = 0.3 expected erroneous reads; Poisson probability of ≥3 errors at a single position is ~0.5%; across 3 × 10⁹ positions, that's ~15 million false positives. A threshold approach fails by five orders of magnitude; you need a model.)

---

## Part 3 — Making Variant Calling Work (≈45 min)

### 3.1 The pre-calling pipeline (≈12 min)

Before a caller runs, three pre-processing steps clean up the BAM. Each removes a specific class of artifact.

**1. Mark duplicates.** PCR amplification during library prep creates multiple reads from the same original template fragment. Those duplicates carry identical sequence and identical errors; counting them as independent evidence overstates the evidence for any artifact they contain. `Picard MarkDuplicates` or `samtools markdup` flags duplicate reads so the caller ignores them. Cost: a few CPU-minutes. Saves: dozens of spurious variants per genome.

**2. INDEL realignment.** A single INDEL near the end of a read is often misaligned as many small substitutions — the aligner preferred to call mismatches rather than open a gap because mismatches were locally cheaper. Local realignment around known or candidate INDEL sites revisits the region with a gap-aware aligner and produces a cleaner BAM.

Historically, GATK's `IndelRealigner` was a required step. Modern callers (GATK4 HaplotypeCaller, DeepVariant) do local re-assembly internally and don't need a separate realignment pass — they look at the reads in windows and construct candidate haplotypes on the fly, effectively doing the realignment inside the caller. If you're running an older pipeline or a simpler caller (bcftools), explicit INDEL realignment still matters.

**FIGURE — Figure #6: INDEL misalignment and local realignment** → `diagrams/lecture-04/06-indel-realignment.svg`
*Top: reads aligned around a true 2-bp deletion, with most reads expressing the deletion as a cluster of 5–6 fake SNVs. Bottom: after local realignment, reads show a clean 2-bp gap and no false SNVs.*

**EMBED — Artifact #2: INDEL Realignment Demo** → `artifacts/lecture-04/02-indel-realignment.html`
*A small region of reads, each with a 2-bp deletion. Toggle "naive alignment" vs "local realignment" and watch the fake SNV cluster disappear, replaced by a clean gap.*

**3. Base-quality score recalibration (BQSR).** The sequencer assigns a Phred quality score to every base. These raw scores are systematically off: the machine uses a calibration learned at manufacturing time, but sequence context (homopolymers, GC-extreme regions) and sample-specific chemistry produce error rates that deviate from the nominal curve. BQSR re-estimates the quality scores empirically by looking at *known* variant sites (dbSNP, known polymorphisms) and assuming any mismatch not at a known site is an error — then bucketing mismatches by context and computing the empirical error rate per bucket. The output is a new per-base quality string that is closer to the empirical truth.

> **EE framing**: BQSR is systematic-error calibration of a measurement instrument. The sequencer is a noisy channel with a systematic distortion — context-dependent transition rates among bases. Calibration against known-good reference positions fits an error model per context, then corrects measurements at unknown positions. The same procedure exists in every serious measurement pipeline: dark-frame subtraction in astronomy, bias-current compensation in ADCs, NTP clock-drift correction in network timing. Without it, your downstream statistics (genotype likelihoods) are computed against an overconfident noise model.

BQSR matters most when your caller is Bayesian and uses quality scores as likelihoods. If the caller treats Q30 as "10⁻³ error probability" but the real error rate at that context is 10⁻², the caller under-weights errors and over-calls variants.

### 3.2 Caller families — from thresholds to deep learning (≈15 min)

Variant callers have evolved through three generations:

**Generation 1 — Heuristic.** Count reads. Apply thresholds. SOAPsnp, early MAQ. Fast, brittle, obsolete.

**Generation 2 — Bayesian.** Treat genotype calling as posterior-probability estimation. For every position, given the observed bases and qualities, compute

```
P(G | D)  ∝  P(D | G) · P(G)
```

where `G` ranges over possible genotypes (`0/0, 0/1, 1/1`), `D` is the observed reads-and-qualities at the position, and `P(G)` is a prior over genotypes given population allele frequencies. The likelihood `P(D | G)` factors across reads: each read's contribution is `P(observed base | true genotype, quality)` computed from the Phred score. Then the caller picks the maximum-posterior genotype.

bcftools and GATK UnifiedGenotyper are second-generation callers. GATK HaplotypeCaller is the modern workhorse — local de-novo assembly of reads in a window, candidate haplotype scoring, and a pair-HMM likelihood for each read against each haplotype. Gives up pure per-position independence for local context.

**FIGURE — Figure #7: Bayesian genotyping at a position** → `diagrams/lecture-04/07-bayesian-genotyping.svg`
*A pileup column; three likelihood curves P(D|0/0), P(D|0/1), P(D|1/1) plotted as functions of read depth; a prior bar; and the resulting posterior with the MAP genotype highlighted.*

**EMBED — Artifact #3: Genotype Likelihood Calculator** → `artifacts/lecture-04/03-genotype-likelihoods.html`
*Enter bases + qualities at a position (or load a preset). The artifact computes P(D|0/0), P(D|0/1), P(D|1/1) using per-read Phred-derived error probabilities, applies a prior, and shows the posterior genotype with the chosen likelihood ratio as the evidence strength.*

> **EE framing**: Bayesian variant calling is posterior estimation with a discrete hypothesis set. The three genotype hypotheses are mutually exclusive, the likelihood factors across independent-read evidence, and the prior is learned from population data. This is the same structure as optimal detection with a finite-alphabet hypothesis set — MAP decoding of convolutional codes, digital-modulation symbol decoding, and classifier-combination in ensemble learning. The Phred score serves as the log-likelihood-ratio contribution of each read, which is why Q30 literally means "one error in 1000" at the likelihood-ratio level.

**Generation 3 — Deep learning.** DeepVariant (Google, 2018) dispenses with an explicit statistical model. Instead:

1. Convert the pileup at every candidate position into a small RGB image — columns are genome positions, rows are reads, pixel colours encode base identity, quality, strand, and mismatch-to-reference.
2. Feed the image through a convolutional neural network trained on millions of benchmarked variant examples.
3. Output: a three-way classification — homozygous reference, heterozygous, homozygous alternate — and a confidence score.

The training data is standardized (Genome in a Bottle's HG001–HG007 samples with orthogonal-truth VCFs), and the CNN — initially Inception-v3, now a custom architecture — learns context-specific features that a Bayesian model would have to encode explicitly. DeepVariant beats every Bayesian caller on benchmark accuracy, at the cost of ~2–3× runtime and GPU dependence.

> **Historical pointer**: The shift from Bayesian to CNN-based calling tracks a wider pattern in ML-heavy bioinformatics. AlphaFold did it for protein structure; DeepVariant did it for variant calling; Inception-like CNNs showed up in pathology and radiology before either. The common structure: a well-defined biological inference task with millions of labeled training examples becomes a supervised-learning problem the moment someone is willing to build the training set. The Bayesian formulation is still useful for understanding what the network is implicitly learning and for calibrating on corner cases the training distribution missed.

### 3.3 Variant-call quality and filtering (≈10 min)

A raw VCF from any caller contains both high-confidence calls and noise. Filtering splits the two. The standard quality signals:

- **QUAL** — site-level variant quality (see §2.3). Usually a per-caller empirical scale, not directly comparable across callers.
- **DP** — too low means too little evidence; too high usually means a mapping-error region (repeat, paralog). Both tails should be filtered.
- **GQ** — genotype quality. Low GQ means the caller couldn't cleanly decide between two genotypes.
- **FS / StrandOddsRatio (SOR)** — strand-bias metrics. A variant supported only by forward-strand reads is probably an artifact.
- **MQ / MQRankSum** — mapping-quality metrics. Variants in low-mapping-quality regions (repeats) are suspect.
- **ReadPosRankSum** — whether variant-supporting reads cluster near read ends (where quality is worst). Clustering = artifact.

Two filtering paradigms:

- **Hard filtering.** Apply per-metric thresholds. `QUAL < 30` → reject. `DP < 10` → reject. Simple, interpretable, tunable. Used by bcftools pipelines, and recommended for small-sample studies where you can't train a model.
- **VQSR / ML-based.** Fit a Gaussian mixture model (GATK's VQSR) or a neural network (hap.py's filter) to a labeled subset of variants from known-truth sites; use the model to classify all sites. Strictly better than hard filtering when you have the data to train it; worse when you don't.

> **Warning box**: Over-filtering is the failure mode nobody talks about. A pipeline tuned to maximise precision (low false-positive rate) by aggressive filtering throws away true variants at low depth or in difficult contexts. A variant you filtered out doesn't generate an error — it generates a missed diagnosis. Always report both precision and recall on a held-out truth set; accepting that 99% precision at 85% recall is worse than 95% precision at 98% recall for most clinical applications.

### 3.4 Somatic vs germline calling (≈8 min)

Germline and somatic callers share the pipeline (aligned BAMs → VCF) but solve different statistical problems.

**Germline callers** assume a diploid genome. At every position, the true genotype is one of three (`0/0, 0/1, 1/1`), and the VAF of a variant should be very close to 0, 0.5, or 1.0. Evidence of a variant at VAF = 0.15 means *something is wrong* — probably a mapping artifact or sample contamination. Germline callers: GATK HaplotypeCaller, DeepVariant, bcftools call, Strelka2 (germline mode).

**Somatic callers** expect a mixture. The sample is a tumor with unknown purity (say, 30% tumor, 70% surrounding normal tissue) harboring a subset of variants each present in a subset of cancer cells. A variant in 50% of tumor cells in a 30% pure sample shows at VAF = 0.15 — that's a real, important variant. The caller has to distinguish it from sequencing noise at the same apparent frequency.

The standard approach: **matched tumor/normal calling.** Sequence both the tumor and a matched normal tissue (usually blood) from the same patient. At every candidate position, compare the tumor pileup to the normal pileup: a true somatic variant is present in the tumor *and absent from the normal*. Callers: Mutect2 (GATK), Strelka2, VarScan2, DeepVariant's somatic mode.

> **EE framing**: Somatic calling is signal detection with a prior subtraction. The "signal" is the tumor-specific variant allele; the "background" is the normal-tissue pileup carrying germline variants that you do *not* want to call as somatic. Subtracting the two distributions (or modelling them jointly, which is what Mutect2 does with a paired-sample likelihood) is the direct analogue of heterodyne detection or differential amplification — measure the difference, not the individual levels. The noise floor drops accordingly.

Somatic callers also have to handle **subclonal structure**: a single tumor is often a mixture of sub-populations (clones) with different variant complements. A variant in 80% of cells is a "clonal" driver; one in 5% of cells is a late "sub-clonal" event. Both are real, and both carry different implications — clonal variants inform what the tumor "is," sub-clonal variants inform resistance and heterogeneity.

---

## Part 4 — Structural Variants (≈45 min)

### 4.1 Why SVs need a different algorithm (≈8 min)

SNVs and small INDELs are within-read events — a single read can carry the complete signal of the variant. A 150 bp read has a SNV visible as one mismatched base; it has a 5 bp deletion visible as a 5 bp gap inside the read's alignment.

Structural variants span distances much larger than a read. A 5 kb deletion cannot be seen inside a 150 bp read — no single read carries the break. Instead, the evidence shows up as a pattern across *multiple* reads: some read pairs land much further apart than expected, some reads split across the breakpoint and half-map to each side, the coverage dips to zero inside the deleted region. A SV caller reconstructs the event from this multi-read evidence.

The three signal channels every SV caller uses:

- **Discordant read pairs.** A paired-end read library has known insert-size distribution (say, mean 400 bp, SD 100 bp). If a read pair lands with 5000 bp between ends, that pair has probably straddled a ~4600 bp deletion. Collect enough discordant pairs clustering around the same two positions and you've localised a deletion.
- **Split reads.** A single read spanning a breakpoint aligns to the reference with the first half matching one location and the second half matching a distant location. The aligner reports both alignment segments; the break in between is the breakpoint itself. Split reads give **single-base-resolution** breakpoints — the gold standard for SV calls.
- **Read depth.** A homozygous deletion drops the coverage inside the deleted region to zero. A heterozygous deletion halves it. A duplication doubles it. Depth-based calling alone can detect large CNVs from a single sample; combined with discordant pairs and split reads, it disambiguates tough cases (inversions vs translocations).

**FIGURE — Figure #8: The three SV detection signals** → `diagrams/lecture-04/08-sv-detection-signals.svg`
*A deletion breakpoint shown with: (a) discordant pairs with larger-than-expected insert size clustering at the breakpoint, (b) split reads mapping half-left, half-right, (c) read depth dropping to zero between the breakpoints.*

> **Intuition box**: SNV calling looks at one column of the pileup at a time. SV calling looks at patterns that span thousands or millions of columns. Different scales, different algorithms. The algorithmic bottleneck is clustering: a thousand discordant pairs near the same ~500 bp region collectively imply a single SV event, and the caller has to group them correctly.

### 4.2 Genome structure context (≈7 min)

SVs reshape the genome at scales that matter for its biology. To understand why a deletion or inversion is consequential, you need a minimum mental model of genome organisation:

- The human genome is about 3.1 Gb of DNA, distributed across 22 autosomes plus the X, Y, and mitochondrial chromosomes.
- Each chromosome is a single linear DNA molecule, with **centromeres** (dense, repeat-rich regions used for chromosome segregation during mitosis) and **telomeres** (protective caps at the ends).
- Genes — the protein-coding and RNA-coding units — are spaced across chromosomes at densities of roughly 1 per 50 kb, with huge regional variation. About 2% of the genome codes for protein.
- The remaining ~98% is non-coding: regulatory regions (promoters, enhancers, insulators), transcribed-but-untranslated RNAs, transposable elements, pseudogenes, and large tracts of "junk" whose function is contested.

Structural variants affect function in two regimes:

- **Local effects** — a deletion removes a gene; an insertion disrupts a splice site; a small inversion flips a regulatory element. Effect is confined to the immediate genomic neighborhood.
- **Global effects** — a translocation joins two chromosomes in a way that fuses two genes (BCR-ABL in CML is the canonical example: a t(9;22) fusion creates a constitutively active tyrosine kinase that drives the cancer). A megabase-scale inversion can disrupt large chromatin domains and affect many genes at once.

### 4.3 SV types and naming (≈10 min)

Structural variants are classified by what they do to the genome, not by their size (though size is correlated). The six canonical types:

- **Deletion (DEL).** A contiguous region is absent from the sample. Loss of sequence. Can be heterozygous (one copy lost, one retained) or homozygous (both copies lost).
- **Insertion (INS).** A sequence is inserted at a position. Can be a new novel sequence, a transposable-element insertion (SINE, LINE), or a tandem duplication of an adjacent sequence.
- **Inversion (INV).** A segment is reversed in orientation — the sequence is the reverse complement of the reference in that region. No loss or gain of sequence, but breakpoints can disrupt genes and regulatory units crossed at the boundaries.
- **Duplication (DUP).** A contiguous region is present in multiple copies. Can be a simple tandem duplication (adjacent copies), a dispersed duplication (copies at other locations), or a segmental duplication (large, historically shared between genomic regions).
- **Copy-number variant (CNV).** The specific case of a DEL or DUP where the number of copies changes. Not a distinct biological mechanism — a CNV is a deletion or duplication measured as a copy-number shift rather than as two breakpoints.
- **Translocation (TRA or BND).** A piece of one chromosome is joined to a different chromosome. Reciprocal translocations swap arms; non-reciprocal translocations insert one chromosome's material into another unidirectionally. Represented in VCFs as breakend ("BND") records rather than as simple deletions/duplications.

**FIGURE — Figure #9: The five canonical SV types** → `diagrams/lecture-04/09-sv-types.svg`
*Five side-by-side panels showing a reference region and the SV-altered version: deletion, insertion, inversion, duplication (tandem), translocation between two chromosomes.*

Two details worth knowing:

- **Complex SVs** — real genomes often contain events that chain multiple primitive operations: an inversion flanked by deletions ("inverted deletion"), a duplication-inversion, chromothripsis (a region shattered and re-assembled). Modern SV callers recognise some of these; most don't.
- **Mobile-element insertions** are a distinct sub-class of INS that deserve their own callers (MELT, xTea). The inserted sequence is a known transposable element (Alu, LINE-1, SVA); the insertion can be typed and annotated even when the breakpoint uncertainty is high.

### 4.4 SV callers and their tradeoffs (≈10 min)

SV callers specialise by signal channel and by input data. The dominant options in 2024:

- **Manta (Illumina / Strelka).** Read-pair + split-read. Short-read focused. Fast, conservative. The workhorse for germline SV calling on clinical WGS.
- **Delly.** Read-pair + split-read. Short-read. More sensitive than Manta but higher false-positive rate. Good for tumor work with aggressive filtering.
- **GRIDSS.** Read-pair + split-read + local assembly. Short-read. State-of-the-art accuracy for complex SVs but slower.
- **SvABA.** Split-read + local assembly. Short-read; supports tumor/normal.
- **LUMPY.** Multi-signal integration framework. Open-source, flexible; requires a separate input from a structural-variant-aware aligner (BWA-MEM + LUMPY or with discordant-pair extractors).
- **CNVkit / cnvpytor.** Depth-based CNV calling. For large (≥ 10 kb) copy-number events; complementary to the breakpoint callers above.

Long reads (PacBio HiFi, ONT) change the picture. A single 20 kb read spans most of an SV, so read-based detection becomes direct: if the read's alignment shows a 5 kb gap, that's the deletion. Long-read-focused callers — **Sniffles**, **pbsv**, **CuteSV** — routinely detect SVs that short-read callers miss, especially those in repeat-rich regions.

**EMBED — Artifact #4: SV Signature Explorer** → `artifacts/lecture-04/04-sv-signatures.html`
*Pick an SV type (deletion, insertion, inversion, duplication, translocation) and see the discordant-pair pattern, split-read pattern, and depth track the event produces. Comparison mode toggles between short reads and long reads to show how long reads simplify detection.*

> **EE framing**: Split-read detection is discontinuity detection on a spatially-aligned signal. The read's alignment score-vs-position profile has a step at the breakpoint — the aligner either confidently maps through the break or it doesn't. Either way, the discontinuity localises the event to single-base resolution. This is the same problem signal-processing tools solve for edge detection, change-point detection in time series, and glitch detection in telemetry. SV callers differ mostly in how aggressively they consolidate multiple weak discontinuities into a single strong event.

### 4.5 Long reads and the future of SV discovery (≈10 min)

Between 2010 and 2020, short-read SV calling was a discipline of workarounds. Every caller hedged against the fundamental limitation that a 150 bp read can't span a 500 bp repeat, and the bulk of SVs live in exactly those repeat-rich regions.

Long reads removed the limitation. PacBio HiFi (15–25 kb, < 0.1% error) and ONT (up to 100 kb+, 5–15% raw error with error correction) span individual SVs directly. The same Telomere-to-Telomere Consortium whose T2T-CHM13 assembly closed the last gaps in the human reference also reported thousands of SVs that were invisible to the short-read-only 1000 Genomes Project dataset.

The current state of SV discovery:

- **Short reads are the workhorse** for general clinical WGS because of cost (~$500/genome vs $1500–3000/genome for HiFi). They catch ~60–70% of SVs accurately.
- **Long reads are the gold standard** for research-grade SV catalogues and for any rare-disease case where short-read analysis has gone nowhere. They catch >95% of SVs, including the hard ones in repeats.
- **Ensemble methods** (calling with multiple tools and intersecting) beat any single caller in benchmarks. The nf-core/structural-variants pipeline and HGSVC's ensemble workflow are the clinical-grade templates.

> **Historical pointer**: The first SV-comprehensive human dataset was HGSVC's 2019 release: three family trios sequenced with short-reads + HiFi + ONT + Strand-seq + BioNano optical maps + Hi-C. Processing produced ~27,000 SVs per haplotype — four to five times as many as short-read-only analyses had reported for the same samples. It was the moment the field accepted that short-read-only SV catalogues were systematically incomplete.

---

## Part 5 — Outputs: VCF and Annotation (≈40 min)

### 5.1 The VCF file format (≈15 min)

VCF (Variant Call Format) is the text format for variant-call output. Every caller emits VCF; every downstream tool reads it. Understanding it is non-optional.

A VCF file has three sections:

- **Header.** Lines starting with `##` are meta-information: reference genome, caller version, contig lengths, INFO-field and FORMAT-field definitions, FILTER definitions. One line starting with `#CHROM` defines the column names for the data rows.
- **Data rows.** One per variant site. Fixed first 8 columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO), then per-sample columns.
- **Compression + indexing.** VCF is plain text; compressed to `.vcf.gz` (bgzip, block-gzip, indexable) and indexed with a `.tbi` file for random access by genomic range.

A simple single-sample example:

```
##fileformat=VCFv4.3
##reference=GRCh38
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Per-sample depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	100234	.	A	G	420	PASS	DP=45;AF=0.51	GT:AD:DP:GQ	0/1:22:23:45:99
chr1	100501	rs12345	C	T	612	PASS	DP=52;AF=1.0	GT:AD:DP:GQ	1/1:0:52:52:99
```

Two rows: a heterozygous SNV at 100234 and a homozygous alternate SNV at 100501 (the second one has a dbSNP rsID in the ID column). Every column has a defined meaning; every piece of software reading the file relies on the header declarations to interpret the INFO and FORMAT fields.

**FIGURE — Figure #10: Anatomy of a VCF file** → `diagrams/lecture-04/10-vcf-anatomy.svg`
*A VCF file with header, column-header, and two data rows; each field annotated with its meaning; a callout showing the FORMAT string and how it maps to per-sample values.*

**EMBED — Artifact #5: VCF Parser and Filter** → `artifacts/lecture-04/05-vcf-parser.html`
*Paste a VCF (or load a preset), see every record parsed into a table, apply QUAL/DP/FILTER filters, and download the filtered output. Handles multi-allelic sites, phased genotypes, and missing data.*

> **EE framing**: VCF is a sparse differential encoding of a genome. Instead of transmitting 3 Gb of sequence per sample, you transmit only the ~4M differences from reference — a compression ratio of ~750×. This is Huffman/entropy coding at the file-format level: the reference is the shared codebook, the VCF is the per-sample delta. Indexable compressed VCFs (bgzip + tabix) are the same design as indexed compressed time-series databases: block-level compression to allow random access into specific regions without decompressing the whole file.

**Multi-sample VCFs.** A single VCF can hold many samples (a cohort, a trio, a tumor/normal pair). Each sample gets its own data column with its own GT, AD, DP, GQ. Sites are included if *any* sample has a non-reference call; samples at sites they don't call get a `0/0` or `./.` entry.

**Phased vs unphased genotypes.** `0|1` (pipe, phased) means the caller knows which allele sits on which of the two parental chromosomes. `0/1` (slash, unphased) means the allele-to-chromosome assignment is unknown. Phasing matters for compound-heterozygous analysis and for trio inheritance — unphased data fundamentally cannot distinguish two heterozygous variants in trans (on different chromosomes) from in cis (on the same chromosome). Long-read sequencing and parental trios are the two reliable paths to phasing.

### 5.2 Filtering VCFs in practice (≈8 min)

A raw VCF straight off a caller is not a final product. Filtering workflow typically:

- **Quality hard filters.** Drop QUAL < 30, DP < 10 or DP > 2×median, GQ < 20. These alone remove most noise without sophisticated training.
- **Region filters.** Drop variants in low-complexity regions, segmental duplications, known-problematic regions (the "decoy" regions and centromeric repeats in GRCh38). The `ENCODE blacklist` is the canonical source.
- **Allele-frequency filters.** For germline research analyses, drop common variants (MAF > 1% in gnomAD) if you're looking for rare-disease variants. For population analyses, keep them.
- **Consequence filters.** After annotation (§5.3), keep only protein-coding consequences — remove intergenic and intronic variants unless you have a specific reason to analyse them.

The filters compose: hard-filter → region-filter → annotate → consequence-filter → output. Every step is a line in a shell script; every step removes orders of magnitude of variants to leave a human-manageable list at the end.

> **Warning box**: Filtering is where pipelines go wrong silently. A too-aggressive filter removes real variants without raising an alarm. Always compute concordance against a benchmark truth set (Genome in a Bottle's HG002 or similar) before trusting a pipeline's output. A 1% drop in recall on the benchmark predicts a 1% miss rate on real samples — unacceptable for a diagnostic pipeline, acceptable for a population-scale study.

### 5.3 Variant annotation (≈12 min)

A VCF says "there's a variant at chr7:140753336, REF=A, ALT=T." That doesn't tell a clinician anything. **Annotation** is the step where each variant is enriched with the information needed to interpret it:

- **Gene and transcript context.** Which gene does this variant fall in? Which exon? Which codon? What amino acid change?
- **Consequence classification.** Missense, nonsense, splice-region, intronic, upstream, intergenic.
- **Population frequency.** How common is this variant in gnomAD? In 1000 Genomes? In population-specific sub-cohorts?
- **Clinical significance.** Is this variant in ClinVar? With what pathogenicity classification?
- **Functional predictions.** Computational scores estimating whether a missense variant is deleterious (SIFT, PolyPhen, REVEL, CADD).
- **Splicing predictions.** Computational scores for variants near splice sites (SpliceAI).

Two annotators dominate the ecosystem:

- **Ensembl VEP** (Variant Effect Predictor). Gold standard, supports every species, maintained by Ensembl. Slowest but most complete. Reads a VCF, outputs an annotated VCF with CSQ fields.
- **snpEff / snpSift.** Faster, Java-based, widely used in research pipelines. Similar capabilities to VEP with different database cadence.

A typical annotation VCF adds INFO fields like:

```
CSQ=T|missense_variant|MODERATE|BRAF|ENSG00000157764|Transcript|ENST00000288602|
    protein_coding|15/18||...|p.Val600Glu|...|REVEL=0.932|SpliceAI=0.01
```

**FIGURE — Figure #11: Variant annotation pipeline** → `diagrams/lecture-04/11-annotation-pipeline.svg`
*A single SNV flowing through a pipeline: input VCF → VEP/snpEff (gene + consequence) → gnomAD join (population frequency) → ClinVar join (clinical significance) → REVEL/CADD (functional prediction) → annotated VCF with all fields present.*

**EMBED — Artifact #6: Variant Annotator** → `artifacts/lecture-04/06-variant-annotator.html`
*Paste a small VCF (or use a preset BRAF / TP53 / BRCA1 example); the annotator shows the gene, consequence, population frequency from a mock gnomAD lookup, a clinical-significance verdict from mock ClinVar entries, and a final interpretation banner (Pathogenic / Likely benign / Uncertain).*

> **Intuition box**: Annotation is joining a variant table to several reference tables. For every variant, look up the gene it's in, the population frequency in gnomAD, the clinical classification in ClinVar, the functional prediction from whatever tool applies. The same variant at `chr7:140753336 A>T` becomes, after annotation, "BRAF V600E, MODERATE impact, gnomAD AF < 1e-5, ClinVar Pathogenic, REVEL 0.932" — which is actionable.

> **Discussion prompt**: Why keep annotation separate from variant calling rather than baking it into the caller? (Keeping them separate means the same caller output can be re-annotated when databases update — new ClinVar classifications, new gnomAD frequencies, new functional predictions. A VCF called in 2020 and annotated with 2024 databases gives a materially different clinical interpretation than one with 2020-era annotation, even though the called variants are identical. Annotation is intentionally the volatile layer.)

---

## Wrap-up (≈10 min)

### What you should take away

- **Variant calling bridges alignment and interpretation.** Aligned reads in; a filtered, annotated VCF out. Every step in between exists because real data breaks the idealised thresholding model.
- **Calling is Bayesian, not thresholding.** Per-position posterior genotype probabilities, computed from per-read evidence weighted by Phred quality. DeepVariant replaces the Bayesian layer with a CNN; the underlying inference structure is the same.
- **The pre-calling pipeline is not optional.** Mark duplicates, INDEL realignment (or haplotype-aware calling), BQSR. Each removes a specific class of artifact; skipping any one leaves a systematic bias that a good caller can't recover from.
- **SVs need different algorithms.** Read pairs, split reads, and depth are the three signal channels. Short-read SV callers (Manta, Delly, GRIDSS) and long-read callers (Sniffles, pbsv) serve different read regimes; ensembling wins on benchmarks.
- **VCF is a sparse encoding.** Understand the field grammar and you can reason about what any variant call means. Annotation is the volatile layer that turns calls into clinical interpretations — and should be re-run whenever databases update.
- **Filtering is where pipelines silently fail.** Over-filtering removes true variants without errors. Always benchmark against a truth set; a 1% recall drop on the benchmark predicts a 1% miss rate on real samples.

### Next lecture

Gene expression analysis: from RNA-seq reads to quantified transcript abundances. Pseudo-alignment vs full alignment, count-to-TPM, differential expression, and the statistical models that underlie DESeq2 and limma.

### Homework

1. Download a public tumor/normal BAM pair (e.g., HCC1143 from the SEQC2 project, or the Genome in a Bottle HG002 reference). Run a germline caller (bcftools or DeepVariant) on the normal. Report: total variants called, PASS fraction, heterozygous/homozygous ratio (should be ~2:1 for a healthy human).
2. Take any VCF from step 1. Annotate with Ensembl VEP (web interface is fine). Report: how many missense variants, how many synonymous, how many nonsense. Does the missense:synonymous ratio match expectations from gnomAD summary statistics?
3. Using the Genotype Likelihood Calculator in §3.2, set up a pileup with 20 reads: 14 REF, 6 ALT, all Q30. Compute P(D|0/0), P(D|0/1), P(D|1/1) by hand (or with the artifact). Which genotype wins, and by what log-likelihood margin?
4. For a hypothetical 2 kb deletion in a heterozygous sample sequenced at 30× short-read coverage: estimate (a) how many discordant read pairs you expect to straddle the breakpoint, (b) how many split reads, (c) the expected depth inside the deletion. Show the arithmetic.
5. Read the VCFv4.3 spec for multi-allelic sites. Explain in one paragraph how a site with `REF=AT, ALT=A,ATT` (a deletion and an insertion at the same position) is represented, and how a downstream tool should handle the two alternates.

### Recommended reading

- McKenna, A., Hanna, M., Banks, E., et al. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research* 20, 1297–1303. (The GATK paper.)
- Poplin, R., Chang, P.-C., Alexander, D., et al. (2018). A universal SNP and small-indel variant caller using deep neural networks. *Nature Biotechnology* 36, 983–987. (The DeepVariant paper.)
- Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. *Bioinformatics* 27, 2987–2993. (The bcftools / samtools calling framework.)
- Chen, X., Schulz-Trieglaff, O., Shaw, R., et al. (2016). Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. *Bioinformatics* 32, 1220–1222.
- Ebler, J., Ebert, P., Clarke, W. E., et al. (2022). Pangenome-based genome inference allows efficient and accurate genotyping across a wide spectrum of variant classes. *Nature Genetics* 54, 518–525. (HGSVC pangenome.)
- VCFv4.3 specification: <https://samtools.github.io/hts-specs/VCFv4.3.pdf>
- GATK Best Practices: <https://gatk.broadinstitute.org/hc/en-us/sections/360007226651>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — What is a Variant?                 | 30 min | 0:30 |
| Part 2 — Variant Calling Mechanics          | 40 min | 1:10 |
| Part 3 — Making Variant Calling Work        | 45 min | 1:55 |
| Part 4 — Structural Variants                | 45 min | 2:40 |
| Part 5 — Outputs: VCF and Annotation        | 40 min | 3:20 |
| Wrap-up                                     | 10 min | 3:30 |

**Total:** ~3h 30min of content.
