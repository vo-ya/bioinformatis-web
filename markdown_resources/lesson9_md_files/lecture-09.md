# Lecture 9 — ChIP-seq, ATAC-seq, and Peak Calling

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-09.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Explain what ChIP-seq and ATAC-seq measure, and name the biological question each answers that DNA- and RNA-level assays cannot.
2. Sketch the ChIP-seq chemistry end-to-end (crosslink → sonicate → immunoprecipitate → sequence) and the ATAC-seq chemistry end-to-end (Tn5 transposition → size-selected fragments → sequence).
3. Read an ATAC-seq fragment-length distribution and identify the sub-nucleosomal, mono-nucleosomal, and di-nucleosomal peaks.
4. Describe MACS2's peak-calling algorithm as local-Poisson detection, derive the local-background lambda estimator, and relate it to CFAR detection in radar.
5. Distinguish narrow-peak from broad-peak calling strategies, and pick the right one for TF ChIP, histone-mark ChIP, and ATAC-seq.
6. Test a peakset for differentially bound or accessible regions using DESeq2/edgeR-style negative-binomial counts — a direct callback to Lecture 6.
7. Read a position weight matrix (PWM), scan a DNA sequence for its hits as a matched-filter problem, and describe the JASPAR / CIS-BP reference databases.
8. Explain ATAC-seq footprinting as an inverse problem from coverage troughs, aggregated over many TF-motif instances to recover the binding shadow.
9. Name two ways modern deep learning predicts regulatory landscapes from DNA sequence (Enformer, Borzoi) and state one thing they do that MACS2 cannot.

---

## Part 1 — What ChIP-seq and ATAC-seq Measure (≈25 min)

### 1.1 The regulation question (≈8 min)

Lectures 1–4 covered DNA sequence: where the letters are, where they differ between samples, and which differences are real. Lectures 5–8 covered RNA: which parts of the genome are transcribed, how much, and in which cell types.

Neither answers a third question: **why** is a gene transcribed in one cell type and not another? The DNA sequence is identical across almost every cell in a body; the expression program is not. Something between genome and transcriptome is controlling which genes are on, and that something is **regulation**: transcription factors binding to DNA, chromatin opening and closing, epigenetic marks being written and erased.

Regulation lives in a different data modality than sequence or transcription. To see it, we need assays that measure binding and accessibility directly:

- **ChIP-seq** (Chromatin ImmunoPrecipitation sequencing) measures where a specific protein binds DNA — any protein with a good antibody.
- **ATAC-seq** (Assay for Transposase-Accessible Chromatin sequencing) measures where DNA is accessible — unwrapped from nucleosomes, available for TFs to read.
- **DNase-seq** (older, similar spirit) measures accessibility via DNaseI digestion instead of transposition.

Both ChIP and ATAC produce sequencing reads that **pile up** at regulatory sites. The analysis task: find the pile-ups, compare them across conditions, and interpret them mechanistically.

> **Intuition box**: A gene is not "on" or "off" because of its sequence alone. The sequence is the contract; the regulatory state is the enforcement. ChIP-seq asks which proteins are reading the contract at this moment; ATAC-seq asks which parts of the contract are open for reading. The combination says which regulatory elements are active in the cells you sampled — a view neither DNA sequencing nor RNA-seq provides.

### 1.2 ChIP-seq — where proteins bind DNA (≈10 min)

A ChIP-seq experiment starts with a cell population and an **antibody** that recognises a specific protein target. The target is most commonly:

- A **transcription factor** (TF) — e.g. CTCF, TP53, FOXA1. The readout: TF binding locations across the genome.
- A **histone mark** — specific chemical modifications to histone proteins, e.g. H3K4me3 (trimethylation of lysine 4 on histone H3 → active promoters), H3K27ac (active enhancers), H3K27me3 (silenced regions), H3K9me3 (heterochromatin). The readout: where that mark decorates the genome.
- A **chromatin remodeler** or **RNA polymerase II** or anything else with a good antibody.

The chemistry:

1. **Crosslink.** Treat cells with formaldehyde. This covalently locks proteins to the DNA they are currently touching — a snapshot of the protein-DNA contact map in that cell state.
2. **Fragment.** Break the cross-linked chromatin into ~200–500 bp fragments by sonication (or by MNase digestion for an alternative called MNase-ChIP).
3. **Immunoprecipitate.** Mix with the specific antibody; pull down antibody-protein-DNA complexes; discard everything else.
4. **Reverse crosslinks, purify DNA, sequence.** The fragments that survive the pull-down came from genomic regions bound by your target protein. Sequencing them shows where.

**FIGURE — Figure #1: What ChIP-seq and ATAC-seq measure** → `diagrams/lecture-09/01-chip-vs-atac-overview.svg`
*Left: a ChIP-seq schematic — crosslinked chromatin is fragmented and pulled down by an antibody specific to a TF or histone mark; surviving fragments are sequenced; reads pile up at bound sites. Right: ATAC-seq — Tn5 transposase inserts sequencing adapters preferentially into open chromatin; reads pile up at accessible sites. Bottom: a cartoon genomic region with TF binding sites, nucleosomes, and open stretches, aligned to the signal each assay produces.*

**FIGURE — Figure #2: ChIP-seq workflow** → `diagrams/lecture-09/02-chip-seq-workflow.svg`
*Four-step schematic: (1) crosslink cells with formaldehyde, (2) sonicate to ~300 bp fragments, (3) pull down with a TF- or histone-specific antibody, (4) reverse crosslinks and sequence the resulting DNA. Each stage labelled with key reagents and expected outputs.*

> **Historical pointer**: The modern ChIP-seq era started with Barski et al. 2007 (*Cell*) and Robertson et al. 2007 (*Nature Methods*), which combined the old ChIP-chip (microarray readout) with the emerging Solexa/Illumina short-read sequencer. The ENCODE Consortium (launched 2003 as ChIP-chip; re-launched in the sequencing era 2007–2012 as "ENCODE 2") then mapped hundreds of TFs and histone marks across dozens of human cell lines. The ENCODE data release (Sept 2012, *Nature* "encyclopedia of DNA elements") remains one of the largest single genomics data releases ever — the default reference that everything downstream is aligned against.

### 1.3 ATAC-seq and DNase-seq — where chromatin is open (≈7 min)

ChIP-seq needs a target and an antibody. What if you want a *map of all regulatory activity* in a tissue without knowing which TFs to ask about?

**ATAC-seq** (Buenrostro et al. 2013) solved this with a clever chemistry trick. The **Tn5 transposase** (from the prokaryotic Tn5 transposon) cuts DNA and inserts a payload in a single step. If you load Tn5 with sequencing adapters as its payload and put it on live chromatin, it inserts adapters into wherever it can physically access the DNA — which is only the **open, nucleosome-free** regions.

The result: short library fragments come out pre-tagged with sequencing adapters, from exactly the open regions. One assay, one enzyme, no antibody needed. The protocol runs from ~50,000 cells in under 4 hours. It transformed epigenomics.

The older analog is **DNase-seq** — instead of Tn5, DNaseI endonuclease preferentially cuts open chromatin. Same output, older, more input DNA required. Most 2024 work uses ATAC.

What ATAC tells you:

- **Where open regulatory regions live** — promoters, enhancers, insulator sites.
- **Where nucleosomes are positioned** — short sub-nucleosomal fragments (< 100 bp) come from accessible linker DNA; mono-nucleosomal (~150 bp) from single-nucleosome-bounded fragments; di-nucleosomal (~300 bp) from two-nucleosome spans. The fragment-length distribution is its own signal (Part 2).
- **Regulatory landscape in rare cell types** — via single-cell ATAC (introduced in L8 §4.2).

What ATAC doesn't tell you: **which TF** is bound in an open region. Only that the region is open. Motif analysis (Part 5) bridges the gap.

**Single-cell ATAC** (callback to L8 §4.2). The same Tn5 chemistry runs at per-cell resolution in droplet scATAC-seq and 10x Multiome (RNA + ATAC on the same nucleus). The data is much sparser — roughly 10,000 fragments per cell vs. ~30 million reads per bulk sample — and dimensionality reduction uses LSI (a TF-IDF-weighted SVD borrowed from text retrieval) instead of the PCA from Lecture 7. But the peak-calling framing in Part 3 largely transfers: you aggregate per-cell fragments into a pseudo-bulk track per cluster, then call peaks as though it were bulk. Everything downstream — differential accessibility (Part 4), motif scanning (Part 5), footprinting (Part 5 with deep-enough aggregate) — runs on the pseudo-bulks.

**EMBED — Artifact #6: ChIP vs ATAC Coverage Viewer** → `artifacts/lecture-09/06-coverage-viewer.html`
*Walk a ~20 kb region with four tracks (CTCF ChIP, H3K4me3, H3K27me3, ATAC). Switch between scenario presets — active promoter, silenced domain, enhancer cluster, active gene body — and watch how the signal across tracks identifies each state. Target aha: no single track resolves the regulatory state; the combination does.*

---

## Part 2 — Library Prep, Biases, and Alignment (≈35 min)

### 2.1 Crosslinking, sonication, Tn5 transposition (≈12 min)

The wet-lab mechanics matter because each step shapes what the analysis sees.

**Formaldehyde crosslinking (ChIP-seq)**. Formaldehyde (HCHO) reacts with primary amines on proteins and with N7 on adenine/guanine in DNA, forming short covalent bonds (~2 Å range). 1% formaldehyde for 10 min at room temperature is standard. Too little: crosslinks don't hold and the IP loses binding information. Too much: you crosslink everything to everything and lose specificity. The concentration and time are part of the protocol; papers report them explicitly.

**Sonication**. Physical breakage of cross-linked chromatin into small fragments using ultrasound. Typical settings produce a distribution centred ~300 bp, with a long tail. Sonication is biased — it preferentially breaks at nucleosome-free regions, subtly depleting fragments from highly packaged heterochromatin. Bioanalyzer traces show the size distribution; the protocol is tuned until the mode sits at the right length.

**Tn5 transposition (ATAC-seq)**. Tn5 is a prokaryotic transposase. In its native role, it cleaves DNA at its recognition sequence and inserts a transposon. For ATAC-seq, Tn5 is pre-loaded with sequencing adapters and applied to a cell's chromatin. It inserts adapters wherever it can physically access the DNA.

Tn5 has a subtle **sequence preference** — it cuts slightly more often at certain ~9-bp motifs than at others. This means the raw cut-site distribution is not uniform even in fully open chromatin. Correction tools (TOBIAS — Part 5, TFBSTools, HINT-ATAC) deconvolve the Tn5 sequence bias before calling footprints.

> **Intuition box**: Crosslinking is a time-stopping trick. Without it, proteins would wander off their sites during the subsequent hours of fragmentation and pull-down, and you'd be sequencing where proteins *were at some point during the experiment*, not where they were **at the moment you added the formaldehyde.** The chemistry freezes the cell state. The quality of the snapshot is set at this step — a bad crosslink produces a blurry experiment no downstream software can fix.

### 2.2 Fragment-length distributions and sequence biases (≈12 min)

Both assays produce libraries with characteristic fragment-length distributions. Looking at the distribution is the first QC step.

**ChIP-seq fragments**: roughly log-normal around the sonication target (~300 bp), with no nucleosomal structure because sonication randomises the breakage points.

**ATAC-seq fragments**: strikingly **multi-modal**, with peaks at:

- **~< 100 bp**: sub-nucleosomal fragments from accessible linker DNA. These are the "nucleosome-free" signal.
- **~150–180 bp**: mono-nucleosomal fragments — the two Tn5 cut sites span exactly one nucleosome (nucleosomes are ~147 bp of DNA wrapped around a histone octamer, plus ~20 bp linker on each side).
- **~300 bp**: di-nucleosomal — two adjacent nucleosomes.
- **~450 bp**: tri-nucleosomal, fainter.

This laddering is a signature of well-performed ATAC-seq. Absence of the pattern = bad experiment (too much or too little Tn5, degraded cells, wrong size selection).

**FIGURE — Figure #3: ATAC-seq Tn5 chemistry** → `diagrams/lecture-09/03-atac-tn5.svg`
*A chromatin cartoon with nucleosomes (histone octamers wrapped with ~147 bp DNA) spaced along a genomic region; Tn5 transposase enzymes shown inserting sequencing adapters in the accessible linker regions; output is a library of short fragments spanning zero, one, or two nucleosomes. Labels indicate the sub-nucleosomal / mono-nucleosomal / di-nucleosomal fragment classes.*

**FIGURE — Figure #4: Fragment-length distribution — nucleosome laddering** → `diagrams/lecture-09/04-fragment-length.svg`
*ATAC-seq fragment-size histogram (log-scale x-axis, linear y). Clear peaks at 50 bp (nucleosome-free), ~180 bp (mono), ~330 bp (di), ~500 bp (tri). Overlay a ChIP-seq fragment distribution (no nucleosomal structure, broad peak at 300 bp) for contrast. Annotations point out each peak's biological origin.*

**EMBED — Artifact #1: Fragment-Size Distribution Explorer** → `artifacts/lecture-09/01-fragment-size.html`
*Simulated ATAC-seq and ChIP-seq fragment distributions. Adjust Tn5 insertion rate, chromatin state, and size-selection; see how the nucleosomal laddering appears, degrades, or disappears. Target aha: the ladder is a signature of healthy, open-chromatin biology — it's the first thing you check before peak calling.*

> **Warning box**: Tn5 has sequence preference. Its preferred insertion motif (roughly 9 bp, approximately `NNNYMNNHN` per Li et al. 2019 bias analyses) means the raw per-base cut counts are not an unbiased measure of accessibility — even in perfectly open DNA, some bases get more cuts than others because Tn5 likes them more. For peak calling the bias averages out over a 200 bp window; for footprinting (Part 5) it does not. TOBIAS, HINT-ATAC, and similar tools explicitly correct for it.

### 2.3 Alignment quirks (≈11 min)

Alignment is mostly a callback to Lecture 2: BWA or Bowtie2 on the paired-end reads, expecting most reads to map uniquely. Two ChIP-and-ATAC-specific wrinkles:

**ATAC-specific Tn5 offset correction.** Tn5 inserts adapters as a dimer, producing a **+4 / −5 staggered cut** on the two strands. The 5' end of a + strand read maps 4 bp after the cut site; the 5' end of a − strand read maps 5 bp before the cut site. To recover the exact Tn5 cut position, shift forward-strand reads +4 bp and reverse-strand reads −5 bp before building the cut-site pile-up. Most ATAC pipelines (e.g. ENCODE) apply this shift routinely. Skipping it blurs footprint analyses by ~9 bp.

**Blacklist regions.** Some parts of the genome produce artifactual high coverage no matter what assay you run — rRNA gene clusters, telomeric / centromeric repeats, regions with annotation errors. These regions systematically produce false-positive peaks in ChIP-seq and ATAC-seq. The **ENCODE Blacklist** (Amemiya et al. 2019) is a curated set of ~400 regions in GRCh38 that every analysis should exclude before peak calling.

**PCR duplicates.** Illumina libraries involve PCR amplification; duplicate fragments from the same pre-amplification template should be collapsed (same as the UMI logic from L7, but for ChIP/ATAC you usually don't have UMIs — use position + fragment length as a duplicate key).

**Shifting of paired-end vs single-end reads.** Some peak callers (MACS2, MACS3) expect you to report fragment length for proper peak modeling. Paired-end reads give it directly; single-end reads require either a separate fragment-size estimation step (cross-correlation analysis, e.g. from SPP) or a configured default.

> **Warning box**: The ENCODE blacklist is version-specific. Using the GRCh37 (hg19) blacklist on a GRCh38-aligned BAM will miss real blacklist regions and fail to filter the right ones. Always match the blacklist to your reference build; nf-core's ATAC-seq and ChIP-seq pipelines pin the correct file per build.

---

## Part 3 — Peak Calling as Detection (≈60 min)

### 3.1 The detection problem (≈10 min)

Suppose you align a ChIP-seq experiment for CTCF. You now have ~30 million reads distributed across the 3-gigabase human genome. Some regions have a strong pile-up (dozens of reads in a 200-bp window); most of the genome has ~0 reads per window.

The **peak-calling task**: given the per-position read count across the genome, declare which windows contain a real signal (a genuine CTCF binding event) and which are background noise.

This is textbook **detection theory** — decide between two hypotheses at every location:

- $H_0$ (null): no binding; the coverage at this position is background noise.
- $H_1$ (alternative): binding event present; the coverage is background + signal.

The signal-to-noise ratio is high in ChIP-seq (a real binding site often has 5–20× the local background count), but the genome is 3 billion positions long, so even rare false alarms — 1 in 10⁵ — produce tens of thousands of spurious peaks. Multiple testing is not optional.

**FIGURE — Figure #5: Read pile-up at a TF binding site** → `diagrams/lecture-09/05-read-pileup.svg`
*A genome browser track style view: a 5 kb window with per-base read coverage plotted on top; a central 200 bp region shows a sharp peak rising well above the baseline coverage; background coverage visible elsewhere. Annotations mark "peak summit", "peak boundaries", and "local background".*

> **EE framing — peak calling as detection**: The peak-calling problem is formally identical to **target detection in radar**: sweep the genome coordinate the way radar sweeps time-delay bins, at each position decide whether the echo amplitude is above a threshold set by the local clutter. The radar tradition distinguishes between **constant-false-alarm-rate (CFAR)** detectors that adapt the threshold to the local noise power, and fixed-threshold detectors that don't. MACS2 is a CFAR detector — it re-estimates the background Poisson rate from a local window around each candidate peak, then thresholds the peak count against that local-adaptive null. Same algorithm your undergraduate detection course would sketch for a phased-array radar.

### 3.2 MACS2 — the algorithm (≈18 min)

**MACS2** (Model-based Analysis of ChIP-Seq, v2 — Zhang et al. 2008, Feng et al. 2012) is the de facto standard. Its core algorithm:

1. **Build a signal track.** Extend each read by the estimated fragment length; compute per-position fragment coverage.
2. **Slide a test window.** For each 200-bp window, count the fragments with summits inside it. Call this the observed count $c$.
3. **Estimate a local background rate $\lambda_{\text{local}}$.** Use the **maximum** of several local background estimates:
   - $\lambda_{\text{bg}}$: genome-wide average count per 200 bp window.
   - $\lambda_{1000}$: average count in a surrounding 1 kb window.
   - $\lambda_{5000}$: 5 kb.
   - $\lambda_{10000}$: 10 kb.
   - (If input/control data is available, $\lambda_{\text{input}}$: equivalent count in matched control.)

   $\lambda_{\text{local}} = \max(\lambda_{\text{bg}}, \lambda_{1000}, \lambda_{5000}, \lambda_{10000}, \lambda_{\text{input}})$.

4. **Compute a Poisson p-value.** Under $H_0$, the count $c$ follows $\text{Poisson}(\lambda_{\text{local}})$. P-value is $P(X \geq c \mid \lambda_{\text{local}})$.
5. **Multi-test correction.** Apply Benjamini-Hochberg across all candidate peaks (callback to Lecture 6 §4) to control FDR at 5% (the default).
6. **Merge adjacent significant windows** into a single peak; estimate the peak summit as the position of maximum coverage within the merged region.

The cleverness is in step 3 — using the **maximum of several local averages** is robust to both local depletion (real biology) and local enrichment (artefacts like blacklisted regions). It's also implicitly CFAR-like: the threshold adapts to local conditions.

**FIGURE — Figure #6: MACS2 local-Poisson detection** → `diagrams/lecture-09/06-macs2-local-poisson.svg`
*A genomic coverage track with three vertical columns: (a) the candidate peak position with count c; (b) local windows at 1 kb, 5 kb, 10 kb on either side with their respective lambda estimates; (c) the Poisson null distribution at lambda_local with the observed c marked on its right tail as "reject H_0". Annotation: "λ_local = max of the local estimates".*

**EMBED — Artifact #2: MACS2 Peak Caller** → `artifacts/lecture-09/02-macs2-peak-caller.html`
*Simulated ChIP-seq coverage track with plantable peaks and tunable background noise. Run the local-Poisson test window-by-window; see which windows are called; compare to ground truth. Adjustable local-window sizes, p-value threshold, BH FDR cutoff. Target aha: raising the p-value threshold inflates false positives; using only global (not local) lambda misses peaks in locally noisy regions.*

> **Historical pointer**: The original MACS paper (Zhang et al., 2008, *Genome Biology*) was one of the first ChIP-seq-specific algorithms, built when the field was six papers old. It replaced naïve fold-enrichment calls that plagued ChIP-seq's first year. MACS2 refined the background model in 2012; MACS3 (2021–) adds ATAC-seq-specific modes and improved memory use. Across three versions and fifteen years, the core detection formulation — observed count vs local Poisson null — has not changed. That is unusual and correct: the detection framing is the right level of abstraction.

### 3.3 Local vs global background (≈12 min)

Why the local-adaptive $\lambda$? Because biology and technology both violate the uniform-background assumption.

**Biological non-uniformity.** Some genomic regions have systematically higher coverage than others. Open chromatin is generally more accessible to any process, so ATAC-seq has extra reads there. Highly expressed genes are more accessible and often more bound. Without local adaptation, the global $\lambda$ under-represents the true background around real hotspots → false positive inflation near real signal.

**Technical non-uniformity.** Mappability varies (regions with many repeats have fewer uniquely-mappable reads → artifactually lower coverage); GC content affects library efficiency; blacklist regions (Part 2) have pathological pile-ups.

**The CFAR solution.** Re-estimate the noise floor locally. The local $\lambda$ averages 1 kb — 10 kb around each candidate, so whatever's elevating the baseline (real biology or artefact) is absorbed into the null distribution. A peak has to rise above the *local* noise, not just the global average.

This is the central idea. Once you see it, peak calling is no longer mysterious.

> **EE framing — adaptive thresholding**: MACS2's `max` across several local-window sizes is conceptually a **guard-cell CFAR** detector. Classical radar CFAR splits the clutter estimate into a leading and lagging guard-cell window around the cell under test, averages the noise estimates from both, and thresholds against that adaptive background. MACS2's 1 kb, 5 kb, 10 kb windows serve exactly this role — they're guard cells centred on the peak candidate, averaged (via `max`, which is a conservative order-statistic choice), and the candidate count is thresholded against that. Radar CFAR has decades of analysis literature on the tradeoffs between mean-level, order-statistic (OS-CFAR), and greatest-of variants; MACS2's choice to take the `max` is a greatest-of-CFAR variant — robust to outliers on the low side, at the cost of slightly elevated miss rate near real peaks adjacent to noise.

**FIGURE — Figure #7: Narrow vs broad peaks** → `diagrams/lecture-09/07-narrow-vs-broad.svg`
*Three coverage tracks aligned to the same x-axis: (a) a TF ChIP (CTCF) — narrow, sharp peaks ~200 bp wide; (b) a histone ChIP for H3K4me3 (active promoter) — broader peaks ~1 kb; (c) a histone ChIP for H3K27me3 (silenced) — very broad, domain-scale peaks spanning 10–100 kb. Each track labelled with its appropriate caller mode (narrow vs broad).*

> **Discussion prompt**: You run MACS2 in default (narrow-peak) mode on a H3K27me3 dataset. The caller finds 50,000 narrow peaks. Your collaborator says "but H3K27me3 forms broad domains of tens to hundreds of kilobases — your narrow calls are fragmenting real domains into pieces." What's the methodological fix? (Use MACS2's `--broad` flag, which merges adjacent windows under a more permissive secondary threshold. For domain-scale marks like H3K27me3, H3K9me3, and H3K36me3, broad-peak calling is the correct mode; narrow-peak output on these marks is an artefact of mismatched algorithm choice. Other callers — SICER, PePr — are designed specifically for broad marks.)

### 3.4 Narrow vs broad peak modes (≈10 min)

**Narrow-peak mode** (MACS2 default). The algorithm described above. Output: sharp, ~100–500 bp peaks. Right for:

- TF ChIP-seq (CTCF, p53, ER, most TFs).
- H3K4me3 (active promoter mark — narrowish).
- H3K27ac (active enhancer mark — narrowish).
- ATAC-seq.

**Broad-peak mode** (MACS2 `--broad`). Uses two thresholds: a strict one to find core regions, a more permissive secondary one to extend/merge into larger domains. Output: domain-scale peaks, possibly tens of kilobases wide. Right for:

- H3K27me3 (Polycomb silenced — domain-scale).
- H3K9me3 (heterochromatin — domain-scale).
- H3K36me3 (gene-body mark — covers whole transcribed regions).

**Superpeak / super-enhancer analysis** (Whyte et al. 2013). A post-processing step on H3K27ac peaks — merge nearby peaks that together span >10 kb of the genome and call the merged region a "super-enhancer". Marks master-regulatory enhancers of cell-type identity. Tools: ROSE (Rank Ordering of Super-Enhancers).

Choosing the right mode matters. Default narrow on a broad mark fragments real domains; default broad on a TF-ChIP over-merges distinct binding events into spurious megapeaks.

---

## Part 4 — Differential Binding and Accessibility (≈30 min)

### 4.1 The DB / DA problem (≈8 min)

You ran ChIP-seq (or ATAC-seq) on two conditions — control vs treated, wild-type vs knockout, normal vs cancer. Same antibody, same protocol, two biological conditions. You want to know **which peaks changed** between conditions.

The analysis template, identical in structure to Lecture 6's differential expression:

1. Call peaks separately in each sample (MACS2 per sample), or on pooled reads across all samples.
2. Build a **consensus peakset**: the union of peaks across all samples, or peaks seen in at least N samples.
3. For each peak × each sample, count reads that overlap. Produces a peaks × samples count matrix — the same shape as a genes × samples matrix from L6.
4. Fit a negative-binomial GLM (DESeq2, edgeR). Test each peak for differential count between conditions. Apply BH correction.

Step 4 is literally a callback to Lecture 6. The count-based DE toolkit works here because the noise model is the same: integer counts with overdispersion, per-feature, across a few dozen samples at most.

**FIGURE — Figure #8: Differential accessibility — MA plot** → `diagrams/lecture-09/08-differential-accessibility.svg`
*MA plot for a differential-accessibility test: x-axis = mean log count across samples; y-axis = log2 fold change between conditions. Points above/below the significance threshold coloured in red/blue. A few highlighted peaks with gene labels (promoters of condition-specific genes). Callback to Lecture 6 §3 volcano / MA conventions.*

> **EE framing — differential accessibility as L6 again**: The count matrix for peaks has the same statistical structure as the count matrix for genes. Same negative-binomial GLM. Same empirical-Bayes dispersion shrinkage. Same Benjamini-Hochberg. The *input features* have changed (genomic intervals instead of transcripts), but the estimator is unchanged. This is a feature of the field: bulk-count DE is a portable tool — once you have counts per feature per sample, the machinery works the same way no matter what the features are.

### 4.2 DiffBind, csaw, and counts across conditions (≈12 min)

Two R packages dominate:

**DiffBind** (Stark & Brown 2011, Ross-Innes et al. 2012). Takes a set of peak files and BAM files; builds the consensus peakset; counts reads per peak per sample; hands off to DESeq2 or edgeR for the statistical test. The wrapper handles the ChIP/ATAC-specific concerns (consensus peakset construction, normalisation choice). Default choice for most labs.

**csaw** (Lun & Smyth 2014, 2016). Alternative approach — window-based rather than peak-based. Tile the genome into overlapping windows; count reads per window per sample; run the DE test on windows; merge significant adjacent windows into differential regions after the test. Avoids the consensus-peakset decision that DiffBind has to make. Higher resolution, more sensitive for broad marks.

**The normalisation question.** Unlike RNA-seq (where you can assume most genes aren't differentially expressed, justifying median-of-ratios), ChIP/ATAC doesn't have that assumption — a whole-genome shift between conditions (e.g. global loss of H3K27ac in a knockout) is biologically real and should *not* be normalised away. DiffBind offers several options: library-size scaling (simplest), TMM (edgeR's default, assumes most peaks unchanged), or spike-in normalisation (requires an added reference species' chromatin). Choose based on whether you expect global shifts.

### 4.3 Gotchas — peakset choice and replicates (≈10 min)

**Consensus peakset choice biases the test.** If you define the consensus as "peaks called in ≥50% of samples," peaks specific to one treatment condition get dropped before the DE test — a severe under-power. If you define as "peaks in any sample," you include noisy single-sample calls and add false positives. Best practice: consensus = peaks called in ≥N samples where N is about half the sample count per condition; or use csaw's window approach to avoid the issue entirely.

**Replicate count.** ChIP/ATAC experiments typically use 2–3 biological replicates per condition. DESeq2/edgeR work down to 2 replicates per group but with much less power than the 5–10 replicate comparisons seen in transcriptomics. Statistics cannot compensate for biology; if the effect needs 5 replicates per arm, use 5.

**EMBED — Artifact #3: Differential Accessibility Explorer** → `artifacts/lecture-09/03-differential-accessibility.html`
*A simulated ATAC-seq experiment — two conditions × 3 replicates × ~500 peaks. Some peaks are truly differential; most are not. Fit a negative-binomial GLM, apply BH, plot MA and volcano. Student adjusts effect size, dispersion, replicate count; sees how power responds.*

> **Warning box**: Differential binding with an altered consensus peakset is a silent-failure trap. If your control has 20,000 peaks and your treatment has 40,000, using "peaks in ≥50% of samples" as the consensus drops half of the treatment-specific peaks entirely — no statistical test can recover what you dropped before the test. Always report how the consensus peakset was constructed and what fraction of single-condition peaks made it in.

---

## Part 5 — Motifs, Matched Filters, and Footprinting (≈35 min)

### 5.1 TF binding motifs and the PWM (≈10 min)

A transcription factor binds DNA with **sequence specificity** — it prefers some 8–20 bp sequences over others. But real binding is not a single consensus sequence; it's a **distribution** over sequences. A TF that mostly likes `TGACTCA` (AP-1) will also bind `TGAGTCA` or `TGACTCT` with lower affinity. The binding specificity is fractional, not binary.

The standard representation is a **Position Weight Matrix (PWM)**:

- 4 rows (A, C, G, T) × *L* columns (motif length, typically 6–20 bp).
- Entry $w_{b, i}$ = log of (observed frequency of base $b$ at position $i$) / (background frequency of $b$). Typically computed from hundreds of known binding sites.

**Scoring** a candidate sequence $s_1 s_2 \ldots s_L$ against a PWM:

$$\text{score}(s) = \sum_{i=1}^{L} w_{s_i, i}$$

Higher score = better match to the TF's binding preference. The score is essentially a log-odds ratio — how much more likely is this sequence to be generated by the TF's binding model than by random background.

**FIGURE — Figure #9: TF motif as PWM and sequence logo** → `diagrams/lecture-09/09-pwm-motif.svg`
*Left: a PWM matrix for CTCF shown as a 4 × 19 numerical table, with weights coloured by sign (positive = accent, negative = subtle). Right: the same motif rendered as a sequence logo (standard IC-height bars; tall CCCTC at the core, variable flanks). Annotations: "each column = one base position", "information content = conservation".*

### 5.2 Motif databases (≈7 min)

Hand-curating PWMs for every TF is expensive. Two community databases collect them:

**JASPAR** (Castro-Mondragon et al. 2022, current version JASPAR 2024). ~2000 curated PWMs across hundreds of vertebrates. Each motif tagged with its TF family, validation source (ChIP-seq, PBM, SELEX), and provenance. The default first-pass database for human and mouse motif scanning.

**CIS-BP** (Weirauch et al. 2014, expanded 2022). ~5000 TF motifs across 300+ species. Denser than JASPAR because it predicts motifs computationally for TFs without direct experimental data (via homology to TFs with known motifs in the same DNA-binding family). Useful when JASPAR has no motif for your TF.

**HOCOMOCO** (Kulakovskiy et al. 2024). ~800 human / mouse motifs, each carefully curated from ChIP-seq data across many experiments. Higher-quality but fewer TFs than JASPAR.

Usage pattern: pick a database, get a PWM for the TF of interest, scan your sequence for hits.

### 5.3 Motif scanning as a matched filter (≈10 min)

Given a PWM and a DNA sequence (say, a 10 kb region around an ATAC peak), find the positions where the motif is likely to be bound.

**Scanning algorithm**:

1. Slide the PWM along the sequence.
2. At each position $i$, compute $\text{score}(s_i s_{i+1} \ldots s_{i+L-1})$ against the PWM.
3. Report positions where the score exceeds a threshold (typically calibrated to a per-base p-value of $10^{-4}$).

Scan both strands (reverse-complement the sequence, or equivalently flip and complement the PWM).

Typical tools: **FIMO** (in the MEME suite), **pwm_scan** (in several bioconductor packages), **MOODS**. They differ in threshold calibration and speed; algorithmically all do the same thing.

**FIGURE — Figure #10: Motif scanning as matched filtering** → `diagrams/lecture-09/10-matched-filter.svg`
*Top: a PWM rendered as a 4-channel template — 4 rows (A, C, G, T) × L columns, with each cell's colour representing the weight. Middle: a DNA sequence rendered as a 4-channel one-hot signal (each position has one row hot). Bottom: the PWM slid along the sequence, producing a correlation-like output track; spikes in the output mark motif hits. Emphasises the signal-processing identity: PWM scan is discrete correlation of the 4-channel sequence against the 4-channel template.*

**EMBED — Artifact #4: PWM Motif Scanner** → `artifacts/lecture-09/04-pwm-motif-scanner.html`
*Input a DNA sequence (or use a preset around a CTCF peak); pick a PWM from JASPAR-style presets; see the scan output as a track aligned under the sequence. Target aha: the hits concentrate in functionally annotated regions; the scan threshold trades sensitivity against false positives.*

> **EE framing — PWM as matched filter**: A PWM scan is a **matched filter** on a 4-channel (A/C/G/T one-hot) signal. Each column of the PWM is a vector of 4 weights; sliding it along the sequence computes, at each position, the inner product between the local 4-channel sequence vector and the template. The output is high at positions where the sequence matches the template, exactly the detection statistic a matched filter produces in any signal-processing course. The "optimal detector under additive Gaussian noise" intuition from EE 101 transfers almost verbatim: log-odds PWM scores on one-hot sequence = matched filter against a 4-channel template. Tools differ in threshold calibration, but the inner operation is this.

> **Warning box**: Motif hit ≠ TF binding. A PWM scan over open chromatin will return far more hits than a TF actually binds in vivo; typical ratios are 10:1 or 100:1 depending on the TF. Binding in real cells is gated by chromatin context, cooperative TF-TF interactions, post-translational modifications of the TF, and many other factors the PWM does not model. Motif hits should be read as "locations this TF *could* bind given the sequence," not "locations it is bound."

### 5.4 ATAC footprinting (≈8 min)

A TF physically bound to DNA protects a ~20 bp stretch from Tn5 insertion — the transposase can't reach the DNA because a protein is in the way. In a ChIP-seq experiment you only know the TF is somewhere in a ~200 bp peak. In an ATAC-seq experiment, if coverage is deep enough, you can see the **footprint**: a local drop in Tn5 cuts at the exact binding site, flanked by normal open-chromatin cuts on either side.

Per-individual-TF-site the footprint is usually too noisy to see. The trick: **aggregate over many instances of the same motif**. Pool the Tn5 cut counts at every position relative to the motif centre, across thousands of motif instances. Random noise averages out; the systematic footprint signal accumulates.

**FIGURE — Figure #11: ATAC footprint at a TF motif** → `diagrams/lecture-09/11-atac-footprint.svg`
*Centre panel: a position-relative-to-motif-centre plot showing Tn5 cut count aggregated over 5,000 CTCF motif instances. Clear dip at the 19-bp motif centre; elevated flanking shoulders. Side panels: individual noisy examples at 5 random instances, then the aggregated clean footprint.*

**EMBED — Artifact #5: ATAC Footprint Analyser** → `artifacts/lecture-09/05-atac-footprint.html`
*Simulate ATAC cut data at N motif instances of a synthetic TF; aggregate to reveal the footprint; adjust coverage depth and motif count; see when aggregation succeeds or fails. Target aha: footprinting is essentially matched-filter averaging over many weak signals — the SNR scales as √N.*

Tools: **TOBIAS** (Bentsen et al. 2020), **HINT-ATAC** (Li et al. 2019), **BaGFoot** (Baek et al. 2017). All implement the same aggregation idea with different bias corrections and statistical thresholds.

> **EE framing — footprinting as inverse problem with averaging**: A single TF footprint is an inverse problem — given a noisy observation at a candidate site, infer whether a protein was bound. Signal-to-noise at a single site is usually poor. Aggregating across thousands of motif instances improves SNR by ~√N (classic signal-averaging argument from any undergraduate measurements course). The footprint only becomes visible in the ensemble, never in the individual. This is why ATAC footprinting needs deep coverage (≥50 million reads per sample) and many motif instances.

> **Intuition box**: The footprint is a **coverage trough**, not a peak. Every other analysis in this lecture detects pile-ups (regions with MORE reads); footprinting detects a local DEFICIT of reads, flanked by normal coverage. Finding something by its absence is a different task than finding it by its presence — which is why footprinting tools are separate from peak callers.

---

## Part 6 — Sequence-to-Regulation Models (≈15 min)

### 6.1 Enformer and its siblings (≈10 min)

Everything up to this point in Lecture 9 is **measurement-driven**: run an assay, observe the data, call peaks or footprints. A different paradigm has emerged: **predict regulatory landscapes directly from DNA sequence**, without running the experiment.

**Enformer** (Avsec, Agarwal et al. 2021, DeepMind + Calico). A CNN + transformer hybrid that takes a ~100 kb DNA sequence and predicts thousands of regulatory outputs: CAGE signal at every position (proxy for TSS activity), ATAC / DNase accessibility, ChIP-seq signal for many TFs and histone marks. Trained on all publicly available ENCODE / Roadmap / FANTOM regulatory data. Performance: within ~15% of experimental replicate reproducibility for many marks — "as good as re-running the experiment" for a significant fraction of features.

**Borzoi** (Linder, Srivastava, Theis, Kelley 2023). Enformer's successor. Extended context (~500 kb), trained on more data including single-cell ATAC.

**DeepSEA / ExPecto / Sei** — earlier CNN-based models with the same goal but less architectural sophistication.

What the sequence-to-function models do that MACS2 can't:

- Predict assay output on **sequences that have never been run** (different species, synthetic sequences, personal genome variants).
- Score **variant effects** — compare model output for reference vs variant sequence to estimate the regulatory impact of a mutation.
- Capture **long-range interactions** — Enformer's 100 kb context lets it learn distal enhancer → promoter relationships from sequence alone.

What they still can't do well: predict accurately in cell types underrepresented in training data; predict rare cell-type-specific regulation; model cooperative-binding nonlinearities beyond what the architecture captures.

**FIGURE — Figure #12: Enformer architecture sketch** → `diagrams/lecture-09/12-enformer-architecture.svg`
*A schematic: DNA sequence → CNN layers (extract local features) → transformer layers (capture long-range interactions) → output heads (separate prediction tracks for CAGE, ATAC, TF ChIPs). Annotate the input context (~100 kb) and the output resolution (~128 bp bins). Forward pointer to Lecture 16's ML-in-genomics synthesis.*

**EMBED — Artifact #7: CNN Filter Bank Visualiser** → `artifacts/lecture-09/07-cnn-filter-bank.html`
*A toy Enformer-first-layer demo: 8 pre-built 8-bp CNN filters initialised to resemble canonical motifs (CTCF, AP1, TATA, GC-box, E-box, etc.). Feed in a DNA sequence; see per-filter activation tracks along the sequence; compare to the PWM-scan output of Artifact #4. Target aha: a learned convolutional filter bank is the same operation as the hand-curated PWM scan, generalised.*

### 6.2 Forward pointer to Lecture 16 (≈5 min)

Sequence-to-function models are one example of a broader pattern that Lecture 16 will synthesize: deep learning architectures matched to genomic problem structure. The pattern so far:

- DeepVariant (L4): CNN on pileups — pileups are image-like; CNNs win.
- scVI / totalVI (L8): VAE with NB likelihood — sparse counts; VAE + count-appropriate likelihood wins.
- Enformer / Borzoi (L9): CNN + transformer on sequence — long-range interactions matter; transformer wins for the global structure, CNN for the local motifs.

Each architecture is *matched* to the inductive biases of the problem. Lecture 16 pulls these patterns out explicitly.

> **EE framing — CNN as learned matched-filter bank**: The first layer of Enformer's CNN is a set of ~256 learned 15-bp convolutional kernels operating on 4-channel one-hot DNA. Structurally this is exactly the PWM matched-filter framing from §5.3 — except the kernels are **learned from data** instead of being hand-curated from known binding sites. Each filter specialises during training to fire on a specific motif-like pattern; across the 256 filters, the bank effectively rediscovers a large subset of known TF motifs (and presumably some unknown ones). Subsequent convolution layers compose these into higher-order patterns — spacing between motifs, cooperative binding signatures. The architecture sits precisely at the intersection of signal processing (matched-filter banks) and learned representations (data-trained kernels), and is the cleanest example in genomics of how classical-DSP intuition survives into the deep-learning era.

---

## Wrap-up (≈10 min)

### What you should take away

- **ChIP-seq measures protein-DNA binding.** ATAC-seq measures chromatin accessibility. Both produce read pile-ups at regulatory regions; peak calling turns pile-ups into a list of regulatory intervals.
- **Library chemistry shapes the downstream signal.** ATAC's Tn5 gives you a fragment-length distribution with characteristic nucleosome laddering — your first QC gate. ChIP's sonication gives you a flat ~300 bp mode. Bias correction matters for footprinting; it matters less for peak calling.
- **Peak calling is detection theory.** MACS2 is a CFAR detector: count at candidate position vs Poisson null estimated from local windows. The local-adaptive background is what makes it work. Narrow peaks for TFs; broad peaks for H3K27me3 / H3K9me3-class marks.
- **Differential binding is Lecture 6 on a different feature set.** Same negative-binomial GLM, same dispersion shrinkage, same BH — on peaks × samples instead of genes × samples. DiffBind / csaw are the standard wrappers.
- **Motif scanning is matched filtering.** PWMs are 4-channel templates; scanning is discrete correlation of template against one-hot sequence. Thresholds trade sensitivity against false positives. Motif hit ≠ binding.
- **ATAC footprinting is signal-averaging over many TF-motif instances.** Aggregation over thousands of sites reveals the coverage trough at the protected binding centre. SNR scales as √N.
- **Deep learning predicts regulatory output from sequence.** Enformer, Borzoi — CNN + transformer architectures trained on ENCODE-scale data. Useful for variant effect prediction and cross-species transfer; less useful for uncommon cell types.

### Next lecture

Methylation, Hi-C, and 3D genome organisation. Bisulfite sequencing as alignment-with-a-known-channel-flip; Hi-C contact maps as correlation matrices; TADs and A/B compartments as eigendecomposition + change-point detection.

### Homework

1. Run MACS2 in narrow and broad mode on the ENCODE K562 CTCF ChIP-seq dataset. Compare the number, median width, and top-peak locations of the two outputs. Which mode fits CTCF biology?
2. Compute the fragment-length distribution of an ATAC-seq sample from ENCODE (e.g. GM12878). Annotate the sub / mono / di / tri nucleosomal peaks. Report the ratio of sub-nucleosomal to mono-nucleosomal fragments; this ratio varies by cell type and is informative.
3. Download two matched ATAC-seq conditions (e.g. control vs LPS-stimulated macrophages from GEO). Call a consensus peakset; build a counts matrix; run DiffBind or csaw; report: how many peaks change, what the top-10 gained / lost peaks are, which genes they fall near.
4. Pick a TF with a JASPAR motif (e.g. CTCF, MA0139.1). Scan a 100 kb region around a gene you're interested in for hits at p < 10⁻⁴. How many hits in open-chromatin regions vs closed? Hit count vs actual ChIP peaks — what's the false-positive rate?
5. Pick a TF with a clear footprint (CTCF is the canonical example). Run TOBIAS or HINT-ATAC on an ATAC-seq sample. Verify: is there a visible footprint at the motif centre? How does footprint depth correlate with local ATAC signal?

### Recommended reading

- Barski, A., Cuddapah, S., Cui, K., et al. (2007). High-resolution profiling of histone methylations in the human genome. *Cell* 129, 823–837. (Early ChIP-seq.)
- Buenrostro, J. D., Giresi, P. G., Zaba, L. C., Chang, H. Y., &amp; Greenleaf, W. J. (2013). Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position. *Nature Methods* 10, 1213–1218. (The ATAC-seq paper.)
- Zhang, Y., Liu, T., Meyer, C. A., et al. (2008). Model-based analysis of ChIP-Seq (MACS). *Genome Biology* 9, R137.
- Feng, J., Liu, T., Qin, B., Zhang, Y., &amp; Liu, X. S. (2012). Identifying ChIP-seq enrichment using MACS. *Nature Protocols* 7, 1728–1740. (The MACS2 protocol paper.)
- ENCODE Project Consortium (2012). An integrated encyclopedia of DNA elements in the human genome. *Nature* 489, 57–74.
- Amemiya, H. M., Kundaje, A., &amp; Boyle, A. P. (2019). The ENCODE blacklist: identification of problematic regions of the genome. *Scientific Reports* 9, 9354.
- Ross-Innes, C. S., Stark, R., Teschendorff, A. E., et al. (2012). Differential oestrogen receptor binding is associated with clinical outcome in breast cancer. *Nature* 481, 389–393. (Early DiffBind / differential ChIP-seq application.)
- Lun, A. T. L., &amp; Smyth, G. K. (2016). csaw: a Bioconductor package for differential binding analysis of ChIP-seq data using sliding windows. *Nucleic Acids Research* 44, e45.
- Castro-Mondragon, J. A., Riudavets-Puig, R., Rauluseviciute, I., et al. (2022). JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles. *Nucleic Acids Research* 50, D165–D173.
- Bentsen, M., Goymann, P., Schultheis, H., et al. (2020). ATAC-seq footprinting unravels kinetics of transcription factor binding during zygotic genome activation. *Nature Communications* 11, 4267. (TOBIAS paper.)
- Avsec, Ž., Agarwal, V., Visentin, D., et al. (2021). Effective gene expression prediction from sequence by integrating long-range interactions. *Nature Methods* 18, 1196–1203. (The Enformer paper.)
- MACS3 documentation: <https://macs3-project.github.io/MACS/>
- TOBIAS tutorial: <https://github.molgen.mpg.de/loosolab/TOBIAS/>
- JASPAR database: <https://jaspar.genereg.net/>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — What ChIP-seq and ATAC-seq Measure          | 25&nbsp;min | 0:25 |
| Part 2 — Library Prep, Biases, and Alignment         | 35&nbsp;min | 1:00 |
| Part 3 — Peak Calling as Detection                    | 60&nbsp;min | 2:00 |
| Part 4 — Differential Binding and Accessibility       | 30&nbsp;min | 2:30 |
| Part 5 — Motifs, Matched Filters, and Footprinting    | 35&nbsp;min | 3:05 |
| Part 6 — Sequence-to-Regulation Models                | 15&nbsp;min | 3:20 |
| Wrap-up                                                 | 10&nbsp;min | 3:30 |

**Total:** ~3h 30min of content.
