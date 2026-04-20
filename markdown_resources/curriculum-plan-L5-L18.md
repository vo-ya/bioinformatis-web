# Bioinformatics for EE Students — Curriculum Plan (L5–L18)

> **Status**: Arc locked in. L1–L4 already drafted. This document captures the plan for L5–L18, including block structure, callback threads, per-lecture notes, and open decisions.
> **Audience**: EE undergraduates / graduates, minimal biology background assumed.
> **Format**: Each lecture targets ~3h–3h 50min of content, matching the L1–L4 pattern.

---

## Curriculum-wide context (recap of L1–L4)

The first four lectures cover DNA end-to-end, one sample at a time:

| # | Lecture | Duration |
|---|---|---|
| L1 | Foundations: From Cells to Sequences to FASTQ | 3h 35min |
| L2 | Read Alignment: From Brute Force to FM Index and Back | 3h 50min |
| L3 | DNA Sequence Assembly: From Reads to Reconstructed Genomes | 3h 25min |
| L4 | Variant Calling: From Aligned Reads to Called Differences | 3h 30min |

**Key concepts and tools already introduced:**
- L1: cell biology, central dogma, sequencing tech (Sanger → NGS → long reads), FASTQ, Phred. Tools: `fastp`, `trim_galore`, FastQC.
- L2: brute force → suffix arrays → k-mer hash → BWT → FM-index, Smith-Waterman, CIGAR, seed-and-extend. Tools: BWA, Bowtie2, Minimap2.
- L3: de novo vs resequencing, Poisson coverage, k-mer error correction, de Bruijn vs overlap graphs, scaffolding, N50. Tools: SPAdes, Flye, Canu, hifiasm, QUAST.
- L4: SNV/INDEL/SV/CNV, pileup, Bayesian genotype calling, DeepVariant (CNN on pileups), BQSR, somatic vs germline, VCF. Tools: GATK HaplotypeCaller, Mutect2, DeepVariant, bcftools, Manta, Delly, GRIDSS, Sniffles, VEP.

**EE framings established:** genetic code ↔ Hamming distance, PCR ↔ positive-feedback amplifier, Phred ↔ dB, alignment ↔ matched filter, BWT ↔ reversible transform, FM-index ↔ compressed lookup, de Bruijn ↔ shift register / Viterbi trellis, genotype calling ↔ MAP decoding, BQSR ↔ instrument calibration, somatic calling ↔ heterodyne detection, VCF ↔ sparse differential encoding.

---

## The post-L4 arc — design rationale

The DNA arc (L1–L4) covered "physical substrate → measurement → alignment → assembly → variants" for one sample. Three axes the DNA arc didn't touch motivate the rest of the course:

1. **DNA → RNA → protein → regulation** — the central dogma in full, plus regulation
2. **One sample → many samples** — populations, cohorts, statistics at scale
3. **Sequence → structure and function** — what the molecules actually *do*

The arc weaves these across 6 blocks:

| Block | Lectures | Focus |
|---|---|---|
| Transcriptomics | L5–L6 | RNA quantification and DE |
| Single-cell | L7–L8 | Single-cell genomics |
| Epigenomics | L9–L10 | Regulation, chromatin, 3D structure |
| Reference & populations | L11–L13 | Long reads, pangenome, popgen, GWAS |
| Infrastructure | L14 | Data engineering, reproducibility |
| Function & structure | L15–L16 | Protein structure, ML synthesis |
| Clinical | L17 | Clinical genomics, ethics |
| Capstone | L18 | Cancer genomics integration |

**Decisions made (with rationale):**
- **Single-cell gets 2 lectures, not 1**: RNA velocity, scVI, batch integration, and spatial all deserve real treatment; cramming kills the EE-flavored content (state estimation, variational inference).
- **Long reads + pangenome gets its own lecture**: T2T (2022) and HPRC (2023) have shifted what "the reference" means; treating long reads as side notes hides this.
- **ML synthesis is late (Option B)**: L16 pulls together architectural patterns from DeepVariant (L4), scVI (L8), Enformer (L9), AlphaFold (L15). Late synthesis lets the lecture be sophisticated about inductive biases rather than abstract.
- **Phylogenetics cut as standalone**: folds into L12 as a coalescent / tree-inference section. Most genomics-relevant face is the coalescent anyway.
- **Data engineering (L14) is its own lecture**: file formats, cloud, containers, workflow languages-as-a-class. Sits mid-course so students appreciate why it matters but before the clinical/capstone where they need it.
- **Clinical/ethics (L17) is its own lecture**: ACMG/AMP, pharmacogenomics, regulatory landscape, ethics. Protected from being squeezed out by the cancer technical integration in L18.
- **Nextflow taught through exercises, not lectures**: workflow languages are syntax, not concepts. L14 covers workflow languages *as a class*; the actual Nextflow teaching happens in homework, culminating in a final pipeline project.

---

## Cross-lecture callback threads

These threads run across multiple lectures. When drafting each lecture, check what to set up and what to pay off.

### Poisson / count-noise thread
- **L3**: introduced as coverage model
- **L6**: returns as negative binomial (overdispersed Poisson) for RNA-seq counts
- **L7**: again with sparse single-cell counts; zero-inflation discussion
- **L9**: again as peak-calling background model (CFAR analog)

### EM / iterative soft-assignment thread
- **L4**: foreshadowed in genotype posteriors
- **L5**: full treatment in Salmon/Kallisto (read-to-transcript assignment under a probabilistic model)
- **L8**: returns in scVI as variational EM
- **EE framing**: iterative soft-decision decoding

### Multiple testing thread
- **L4**: introduced informally in variant filtering
- **L6**: full treatment with BH at 20K-gene scale
- **L13**: scaled up to 10M-SNP GWAS scale; genome-wide significance threshold

### Seed-and-extend thread
- **L2**: established
- **L5**: generalized to splice-aware (long-skip transitions)
- **L11**: generalized to graph alignment on a DAG

### Reference-as-string-is-a-lie thread
- **L4**: hinted at via reference bias in variant calling
- **L11**: payoff in pangenome graphs

### ML-in-genomics thread
- Scattered: DeepVariant (L4), scVI (L8), Enformer (L9), AlphaFold (L15)
- **L16**: synthesis lecture pulls patterns out

### Cancer thread (deliberately seeded for L18 capstone)
- **L4**: somatic vs germline calling
- **L7**: tumor heterogeneity from single-cell perspective
- **L4 + L11**: CNVs from short reads / SVs from long reads
- **L13/L17**: cancer drivers, germline cancer risk
- **L18**: integrated capstone

---

## Per-lecture notes

### L5 — Bulk RNA-seq: splice-aware alignment + quantification (~3h 30min)

**Block**: Transcriptomics

**Proposed sections:**
1. **Biology and why RNA ≠ DNA counts** (~25 min). Transcription, splicing, isoforms, why naive read-per-gene counting fails. Sets up the rest of the lecture.
2. **Splice-aware alignment** (~50 min). STAR, HISAT2. Natural extension of L2: seed-and-extend where reads can span introns. EE framing: HMM with explicit skip states; trellis with one transition edge of enormous length.
3. **Pseudoalignment and quantification** (~60 min). Kallisto, Salmon, RSEM. **Intellectual high point of the lecture** — EM for read-to-transcript assignment ↔ iterative soft-decision decoding. K-mer compatibility classes connect to L2's hash indexing but used for set membership.
4. **Count models and normalization** (~30 min). TPM/CPM/size factors, intro to negative binomial. Callback to L3 Poisson.
5. **Bridge to L6**: state the DE problem, motivate count statistics.

**Pedagogical choice**: keep historical order (splice-aware first, then pseudoalignment). The punch of "you don't actually need base-level alignment for quantification" lands harder after students have felt the splice-aware cost.

**Tools introduced**: STAR, HISAT2, Salmon, Kallisto, RSEM, featureCounts, HTSeq.

**Key EE framings**:
- Splice-aware alignment ↔ HMM with long-skip transitions
- EM for read assignment ↔ iterative soft-decision decoding
- Pseudoalignment ↔ sketch-based set membership / Bloom-filter-adjacent

---

### L6 — Differential expression and count statistics (~3h)

**Block**: Transcriptomics

**Proposed sections:**
1. **The DE problem framed properly** (~20 min). Why a t-test on counts is wrong; what the noise model needs to capture.
2. **Negative binomial as overdispersed Poisson** (~40 min). Direct callback to L3. Biological vs technical variance. Why overdispersion is non-negotiable for RNA-seq.
3. **DESeq2 and edgeR mechanics** (~50 min). Size factor estimation, dispersion shrinkage (empirical Bayes), Wald and LRT tests. Limma-voom as the alternative.
4. **Multiple testing at gene scale** (~30 min). Bonferroni vs BH; FDR control intuition; q-values vs p-values.
5. **What to do with a gene list** (~20 min). Brief intro to GSEA, ORA, pathway databases — sets up what would have been L14 functional annotation, now distributed.

**Tools introduced**: DESeq2, edgeR, limma-voom, GSEA, fgsea, clusterProfiler.

**Key EE framings**:
- Negative binomial ↔ overdispersed Poisson (callback to L3)
- Empirical Bayes shrinkage ↔ regularization in low-data regimes
- BH/FDR ↔ multi-channel detection with controlled false-alarm rate
- Dispersion estimation ↔ noise floor characterization

---

### L7 — Single-cell RNA-seq fundamentals (~3h 30min)

**Block**: Single-cell

**Proposed sections:**
1. **Why single-cell** (~15 min). What bulk hides; cell-type heterogeneity; the questions you can only ask at single-cell resolution.
2. **Droplet protocols and chemistry** (~30 min). 10x Genomics workflow; cell barcodes; UMIs and the PCR duplicate problem (clean callback to L1's PCR ↔ amplifier framing).
3. **From BCL to count matrix** (~30 min). Cell Ranger / STARsolo / alevin pipelines; the cell-by-gene matrix as the canonical output.
4. **QC and preprocessing** (~30 min). Empty droplets (EmptyDrops), doublets (Scrublet, DoubletFinder), mitochondrial fraction filtering, ambient RNA (SoupX, CellBender).
5. **Normalization and feature selection** (~25 min). Why bulk methods fail on sparse data; SCTransform, log-normalization, HVG selection.
6. **Dimensionality reduction and clustering** (~40 min). PCA (linear) → UMAP/t-SNE (nonlinear); Leiden/Louvain clustering; marker gene identification (Wilcoxon, MAST).
7. **The cell-type annotation problem** (~20 min). Reference-based (SingleR, Azimuth) vs marker-based; ontologies (CL, HCA).

**Tools introduced**: Cell Ranger, STARsolo, alevin-fry, Seurat, Scanpy, EmptyDrops, Scrublet, SoupX, SingleR.

**Key EE framings**:
- UMIs ↔ error-correcting tags solving the PCR duplicate problem from L1
- Sparse count matrix ↔ undersampled signal per cell
- UMAP ↔ nonlinear manifold learning; PCA as prior linear projection
- Doublet detection ↔ outlier detection in sparse high-dim space

---

### L8 — Advanced single-cell: trajectories, integration, multi-modal (~3h 30min)

**Block**: Single-cell

**Proposed sections:**
1. **Pseudotime and trajectory inference** (~40 min). Monocle, Slingshot, PAGA. The intuition: order cells along a developmental axis from snapshot data.
2. **RNA velocity** (~40 min). **Strong EE-flavored content.** Spliced/unspliced ratios as derivative estimation; velocyto and scVelo; the dynamical model. EE framing: state-space estimator; literally a derivative estimated from a ratio of measurements.
3. **Batch integration** (~45 min). Why you can't just concatenate experiments; Harmony (linear correction in PCA space), scVI (deep generative model — VAE on counts). EE framing: source separation across batches.
4. **Multi-modal single-cell** (~30 min). CITE-seq (RNA + protein), single-cell ATAC, single-cell multi-ome. Joint embeddings (WNN, MOFA+, totalVI).
5. **Spatial transcriptomics** (~40 min). Visium, MERFISH, Xenium. Image + transcriptome joint modeling; deconvolution of mixed spots. Brief tour, not exhaustive.
6. **Cell-cell communication** (~15 min). CellChat, NicheNet — ligand-receptor inference. Caveats around interpretation.

**Tools introduced**: Monocle3, Slingshot, scVelo, Harmony, scVI, totalVI, Squidpy, CellChat, NicheNet.

**Key EE framings**:
- RNA velocity ↔ state-space estimator from spliced/unspliced ratios
- scVI ↔ variational EM (callback to L5 Salmon EM thread)
- Batch integration ↔ source separation
- Spatial deconvolution ↔ mixed-pixel unmixing in remote sensing

---

### L9 — ChIP-seq, ATAC-seq, and peak calling (~3h 15min)

**Block**: Epigenomics

**Proposed sections:**
1. **What we're measuring** (~25 min). Where proteins bind DNA (ChIP-seq); where chromatin is open (ATAC-seq, DNase-seq). Why these are biology's main tools for studying regulation.
2. **Library prep and biases** (~25 min). Crosslinking, sonication, Tn5 transposase mechanics. GC bias, fragment-length distributions, Tn5 sequence preference.
3. **Alignment with quirks** (~20 min). Brief — mostly a callback to L2. ATAC-specific Tn5 offset correction (+4/-5 shift). Blacklist regions (ENCODE blacklist).
4. **Peak calling as detection** (~50 min). **Major EE-framed section.** MACS2 algorithm; signal vs background; local Poisson vs global Poisson; CFAR analog. Narrow vs broad peaks.
5. **Differential binding / accessibility** (~30 min). DiffBind, csaw. Returns to count statistics from L6 in a new domain.
6. **Footprinting and motif analysis** (~30 min). TF binding site discovery (MEME, HOMER); footprinting in ATAC-seq (TOBIAS); the JASPAR/CIS-BP databases.
7. **Forward pointer to Enformer** (~15 min). Modern deep learning approach to predicting regulatory landscapes from sequence. Sets up L16.

**Tools introduced**: MACS2, MACS3, DiffBind, csaw, TOBIAS, MEME, HOMER, deepTools, ChIPseeker.

**Key EE framings**:
- Peak calling ↔ CFAR detection in radar
- Local vs global background ↔ adaptive vs fixed threshold detection
- Motif PWM ↔ matched filter on a 4-channel signal
- Footprinting ↔ inverse problem from coverage troughs

---

### L10 — Methylation, Hi-C, and 3D genome organization (~3h 15min)

**Block**: Epigenomics

**Proposed sections:**
1. **DNA methylation biology** (~20 min). 5mC, CpG islands, why methylation matters for gene regulation and imprinting.
2. **Bisulfite sequencing and its alignment problem** (~40 min). Chemistry: unmethylated C → U → T after PCR. Why this breaks L2 aligners. Bismark, methylpy, BSMAP — three-letter alignment, in-silico conversion strategies. Recent: long-read direct methylation calling (ONT, PacBio HiFi 5mC).
3. **Differential methylation** (~25 min). DMR calling, methylKit, BSmooth. Beta distribution for proportion data.
4. **Hi-C and chromosome conformation capture** (~40 min). 3C → 4C → 5C → Hi-C → Micro-C lineage. Library prep mechanics; the contact matrix as raw output.
5. **From contacts to structure** (~40 min). Iterative correction (ICE), KR normalization. TADs, A/B compartments, loops. Tools: HiC-Pro, cooler, HiCExplorer, Juicer.
6. **EE framing block** (~20 min). Hi-C contact maps ↔ symmetric correlation matrices; TAD calling ↔ block-diagonal structure detection ↔ image segmentation; A/B compartments ↔ first eigenvector of normalized contact matrix.

**Tools introduced**: Bismark, methylKit, HiC-Pro, cooler, Juicer, HiCExplorer.

**Key EE framings**:
- Bisulfite alignment ↔ alignment with a known channel-induced character flip
- Hi-C contact matrix ↔ correlation/covariance matrix
- A/B compartments ↔ first eigenvector decomposition (PCA on contacts)
- TAD calling ↔ block-diagonal structure detection / change-point detection
- 3D reconstruction from Hi-C ↔ multidimensional scaling on a distance matrix

---

### L11 — Long reads and the pangenome (~3h 30min)

**Block**: Reference & populations

**Why this is its own lecture**: T2T (2022) and HPRC (2023) have shifted what "the reference" is. Treating long reads as side notes hides the field-wide transition. PacBio HiFi is now ~99.9% accurate; ONT is past Q20.

**Proposed sections:**
1. **Long-read tech recap** (~20 min). Quick callback to L1. PacBio HiFi vs ONT; what each is good at in 2024–2026; why accuracy improvements changed the calculus.
2. **What long reads unlock** (~40 min). SVs you can resolve directly (callback to L4); repeats you can span; haplotype phasing without statistical inference; complete telomere-to-telomere assemblies.
3. **The pangenome shift** (~45 min). Why a single linear reference is a lie — reference bias in non-European populations as a concrete, ethically-loaded example. The graph genome concept. GFA format.
4. **Graph-based alignment** (~50 min). **Major EE-framed section.** How seed-and-extend (callback to L2) generalizes when the reference is a DAG. Tools: vg, minigraph, GraphAligner. EE framing: Viterbi on an arbitrary DAG vs a linear chain.
5. **T2T and HPRC** (~30 min). The state of the art. The CHM13 telomere-to-telomere assembly. The 47-haplotype HPRC v1 release. What's still missing.
6. **Phasing and haplotype-resolved assembly** (~25 min). hifiasm in trio mode; HapCUT2; WhatsHap. EE framing: source separation.

**Tools introduced**: vg, minigraph, GraphAligner, pggb, hifiasm (revisited), HapCUT2, WhatsHap.

**Key EE framings**:
- Graph alignment ↔ Viterbi on arbitrary DAG (vs L2's linear chain)
- Pangenome ↔ population-as-codebook representation problem
- Phasing ↔ source separation
- Reference bias ↔ training/test distribution mismatch

---

### L12 — Population genetics fundamentals (~3h 30min)

**Block**: Reference & populations

**Note**: phylogenetics folds in here as a section on coalescent / tree inference (its most genomics-relevant face).

**Proposed sections:**
1. **Allele frequencies and Hardy-Weinberg** (~30 min). The basic accounting. HWE as a null hypothesis; departures as biological signal.
2. **Forces shaping variation** (~45 min). Drift, selection, mutation, migration. Effective population size. The Wright-Fisher model.
3. **Linkage disequilibrium** (~40 min). LD measures (r², D'); LD decay with recombination; LD blocks; haplotype structure. EE framing: LD as autocorrelation along the genome.
4. **The coalescent** (~50 min). Backward-in-time genealogy; coalescent trees; relationship to forward Wright-Fisher; ARGs (ancestral recombination graphs). This is where phylogenetics lives in the modern workflow.
5. **Inferring population history** (~40 min). PSMC, SMC++, Relate. Demographic inference from a single genome; admixture inference (ADMIXTURE, STRUCTURE).
6. **Selection scans** (~25 min). iHS, XP-EHH, Tajima's D. Detecting positive selection in modern populations.

**Tools introduced**: PLINK, ADMIXTURE, PSMC, SMC++, Relate, RAxML (brief), IQ-TREE (brief).

**Key EE framings**:
- LD ↔ autocorrelation along the genome
- Coalescent ↔ branching process / inverse-time stochastic process
- Selection scans ↔ detection of structured signal in a stochastic background
- Admixture inference ↔ mixture-model decomposition

---

### L13 — GWAS and statistical genetics (~3h 30min)

**Block**: Reference & populations

**Proposed sections:**
1. **The GWAS framework** (~30 min). Cohort design; cases vs controls; quantitative traits. The basic regression at every SNP.
2. **The Manhattan plot and what it shows** (~25 min). Genome-wide significance (5×10⁻⁸); QQ plots and inflation; LD-driven peak structure.
3. **Multiple testing at SNP scale** (~30 min). Why 5×10⁻⁸ specifically (LD-corrected Bonferroni for ~10⁶ independent tests). Callback to L4 and L6 thread.
4. **Confounders and corrections** (~40 min). Population stratification (PCA correction, genomic control); cryptic relatedness (kinship matrices); mixed models (BOLT-LMM, SAIGE, REGENIE).
5. **Beyond single-SNP** (~30 min). Heritability estimation (LDSC, GCTA); polygenic risk scores (PRS) and their portability problems across populations.
6. **Fine-mapping and colocalization** (~30 min). From a peak to a causal variant: SuSiE, FINEMAP, CAVIAR. Colocalizing GWAS signal with eQTL signal (coloc).
7. **Rare-variant association** (~25 min). Burden tests (SKAT, SKAT-O); the gnomAD-scale rare variant landscape.

**Tools introduced**: PLINK (revisited), BOLT-LMM, SAIGE, REGENIE, LDSC, GCTA, SuSiE, coloc, SKAT.

**Key EE framings**:
- GWAS ↔ massively-multichannel detection with controlled FAR
- Population structure correction ↔ noise whitening / decorrelation
- Mixed models ↔ generalized least squares with structured covariance
- Fine-mapping ↔ sparse inverse problem (which of N correlated signals is causal)
- LD-based heritability ↔ spectral decomposition of the kinship matrix

---

### L14 — Data engineering, file formats, and reproducibility (~3h)

**Block**: Infrastructure

**Why this lecture exists**: Industry-relevant content that doesn't fit elsewhere. Sits mid-course so students appreciate it but before clinical/capstone where they need it.

**Proposed sections:**
1. **File formats and their tradeoffs** (~40 min). FASTQ → BAM/CRAM (CRAM as content-aware compression); VCF/BCF; HDF5/Zarr for single-cell (AnnData, MuData); Parquet for tabular omics; GFA for graphs (callback to L11). What makes a format good (random access, compression, schema evolution).
2. **Where data lives** (~30 min). SRA, ENA, GEO, dbGaP (controlled access), EGA. Cloud-hosted public data (1000G on AWS, gnomAD). Requester-pays buckets and the politics of who pays for compute.
3. **Reference data management** (~20 min). iGenomes, refgenie, the perils of "which GRCh38" (decoy contigs, alt loci, ALT-aware vs ALT-unaware).
4. **Containerization** (~30 min). Why bioinformatics has *especially bad* dependency hell (R + Python + C tools + system libs). Docker, Singularity/Apptainer, conda/mamba/pixi. BioContainers.
5. **Workflow languages as a class** (~40 min). Nextflow vs Snakemake vs WDL vs CWL. Resource management; resume on failure; cloud execution. **Nextflow itself is taught in homework, not lectures** — this is the conceptual framing.
6. **Benchmarking and validation** (~25 min). GIAB truth sets; precision-recall on variant calls (hap.py); how the field actually evaluates new methods. Reproducibility crises and what they look like.
7. **The culture** (~15 min). Preprints (bioRxiv), the role of nf-core, community-curated pipelines, why bioinformatics has a particularly strong open-source norm.

**Tools introduced/discussed**: samtools, bcftools (revisited), CRAM, AnnData, refgenie, Docker, Singularity, conda, mamba, Nextflow, Snakemake, nf-core, GIAB, hap.py.

**Key EE framings**:
- CRAM ↔ content-aware lossless compression (vs BAM's generic gzip)
- Workflow DAGs ↔ dataflow programming
- Containerization ↔ environment isolation as a hardware-virtualization analog
- GIAB truth sets ↔ ground-truth datasets for ML benchmarking

---

### L15 — Protein structure prediction (AlphaFold-era) (~3h 30min)

**Block**: Function & structure

**Note**: This is mostly an ML lecture in disguise — and that's fine, it's exactly what EE students will engage with. Sets up L16's synthesis.

**Proposed sections:**
1. **Protein structure basics** (~25 min). Primary → secondary → tertiary → quaternary. Domains, folds, the Anfinsen experiment. Why structure matters for function.
2. **The classical era** (~25 min). Homology modeling (Modeller, SWISS-MODEL); threading (I-TASSER); CASP as the field's benchmark. Why this plateaued.
3. **MSA-based methods and coevolution** (~30 min). The pre-AlphaFold breakthrough: contact prediction from MSA covariation (DCA, EVfold). Why coevolving residues are spatially close.
4. **AlphaFold2 architecture** (~60 min). **Major content block.** The Evoformer (axial attention on MSA + pair representation); the structure module (equivariant transformer with IPA); recycling. Confidence metrics: pLDDT, PAE.
5. **AlphaFold3 and successors** (~25 min). Diffusion-based generation; multi-chain and ligand prediction; ESMFold (single-sequence). RoseTTAFold and the open-source alternatives.
6. **Inverse folding and protein design** (~25 min). ProteinMPNN, RFDiffusion. Designing sequences that fold to a target; designing structures de novo.
7. **What this changes** (~20 min). Drug discovery, enzyme engineering, the AlphaFold database (200M predicted structures). What's still hard: dynamics, disorder, multi-conformation states.

**Tools introduced**: AlphaFold2/3, ColabFold, ESMFold, RoseTTAFold, ProteinMPNN, RFDiffusion, PyMOL/ChimeraX (visualization).

**Key EE framings**:
- MSA covariation ↔ inverse covariance estimation (graphical lasso)
- Evoformer ↔ axial attention (efficient 2D attention pattern)
- Structure module IPA ↔ SE(3)-equivariant networks
- Diffusion models ↔ denoising score matching / iterative refinement
- pLDDT ↔ calibrated confidence estimation
- Inverse folding ↔ inverse problem with strong structural priors

---

### L16 — ML in genomics: architectures, pitfalls, frontiers (~3h 30min)

**Block**: Function & structure

**Why late synthesis (Option B)**: Students have already seen DeepVariant (L4), scVI (L8), Enformer (L9), AlphaFold (L15). Pulling the patterns out is more sophisticated when grounded in concrete examples.

**Proposed sections:**
1. **Architectural pattern survey** (~50 min). For each pattern, name the genomic problem and the matching architecture: pileups → CNN (DeepVariant); long-range sequence regulation → dilated conv + transformer (Enformer, Borzoi); count matrices → VAE (scVI, totalVI); molecular graphs → GNN; protein 3D → equivariant transformer (AlphaFold).
2. **Inductive biases — why does each work?** (~40 min). What about each architecture matches what about each problem. Why CNNs work for pileups (translation invariance) but transformers won for protein structure (long-range pair interactions). Why count VAEs need the negative binomial likelihood (callback to L6).
3. **Training data is the bottleneck** (~30 min). Where labeled data comes from in genomics; why it's never enough; weak supervision; self-supervised pretraining.
4. **The data leakage problem** (~30 min). **Critically important and under-taught.** Train/test splits in genomics are *hard*: homologous sequences, related individuals, batch effects across studies. Splitting by chromosome, by species, by family. Why naive random splits are catastrophically optimistic. Real published papers that fell into this trap.
5. **DNA language models** (~30 min). DNABERT, Nucleotide Transformer, HyenaDNA, Evo. What pretraining on DNA gives you; what it doesn't. Honest assessment of the current evidence.
6. **Foundation models for cells** (~25 min). scGPT, Geneformer, scFoundation. Genuinely useful or premature? An honest tour.
7. **What's coming** (~15 min). Multimodal foundation models; protein-DNA co-design; the convergence with AI/ML at large.

**Tools/methods discussed**: DeepVariant, Enformer, Borzoi, scVI, AlphaFold (all revisited), DNABERT, Nucleotide Transformer, HyenaDNA, Evo, scGPT, Geneformer.

**Key EE framings**:
- Architecture choice ↔ inductive bias matching the problem structure
- Self-supervised pretraining ↔ representation learning from unlabeled data
- Data leakage ↔ train/test contamination in correlated data
- DNA LMs ↔ language models on a 4-letter alphabet with very different statistics than text

---

### L17 — Clinical genomics, variant interpretation, and ethics (~3h 30min)

**Block**: Clinical

**Why this lecture exists**: Graduating EE students into the genomics industry without ACMG/AMP, regulatory context, and ethics is a real disservice — they'll be writing code that affects patient care.

**Proposed sections:**
1. **From research to clinic** (~25 min). What changes when a pipeline is clinical: validation, reproducibility, turnaround time, audit trails. CLIA labs, CAP accreditation. LDTs vs IVDs vs RUO assays.
2. **ACMG/AMP variant classification** (~50 min). **Core content.** The 2015 guidelines: pathogenic, likely pathogenic, VUS, likely benign, benign. The PVS1/PS/PM/PP and BA1/BS/BP evidence codes. How the rules combine. Sherloc, ClinGen refinements. ClinVar as the community's classification database.
3. **Variant evidence ecosystem** (~30 min). Population frequency (gnomAD); functional impact predictors (REVEL, AlphaMissense); splice predictors (SpliceAI); ClinVar history; literature mining. How a variant scientist actually works.
4. **Pharmacogenomics** (~30 min). Real worked examples: CYP2D6 and codeine metabolism; warfarin dosing (VKORC1, CYP2C9); HLA-B*57:01 and abacavir; TPMT and thiopurines. PharmGKB, CPIC guidelines, the actual delivered-to-patient examples.
5. **Incidental findings** (~25 min). The ACMG SF list (currently SF v3.x — 80+ genes). What to do when you find something unrelated to the indication. The right-not-to-know debate.
6. **Regulatory landscape** (~25 min). FDA approval of sequencing assays (FoundationOne, MSK-IMPACT). The 2023 FDA LDT rule. EU IVDR. Direct-to-consumer regulation (23andMe FDA history).
7. **Ethics, equity, and data sovereignty** (~25 min). GINA and what it does/doesn't protect. Genetic discrimination in life and disability insurance (GINA carveout). Ancestry bias in clinical databases. The Havasupai case. Indigenous data sovereignty (CARE principles). All of Us, Our Future Health, and what diverse cohorts are trying to fix.

**Tools/databases introduced**: ClinVar, gnomAD (revisited), VarSome, Franklin, InterVar, REVEL, AlphaMissense, SpliceAI, PharmGKB.

**Key EE framings**:
- ACMG/AMP rules ↔ structured rule-based classification (vs end-to-end ML); auditability requirement
- Clinical pipeline validation ↔ FDA/regulatory validation in EE (medical devices, RF certification)
- Variant pathogenicity prediction ↔ multi-evidence Bayesian aggregation
- Ancestry bias ↔ training distribution coverage problem from L16

---

### L18 — Cancer genomics: integrated capstone (~3h 30min)

**Block**: Capstone

**Purpose**: Walk students through what the industry looks like by integrating L1–L17 on a real, important problem. Cancer is ideal because it touches *every* prior lecture.

**Lectures it integrates** (worth being explicit so the lecture earns its capstone status):
- L1: sample prep matters (FFPE artifacts)
- L2–L4: alignment and somatic variant calling (Mutect2, Strelka)
- L5–L6: tumor expression profiling (RNA-seq for fusions, expression subtypes)
- L7–L8: tumor heterogeneity (single-cell tumor + microenvironment)
- L9–L10: epigenomic dysregulation (cancer methylation, accessibility)
- L11: structural variants (long reads for cancer SVs)
- L12–L13: germline cancer risk (BRCA, Lynch syndrome)
- L14: clinical pipeline reproducibility
- L15: drug-target structure (kinase inhibitor design)
- L16: ML in pathology and genomics
- L17: clinical interpretation, ACMG/AMP for cancer (AMP/ASCO/CAP guidelines for somatic)

**Proposed sections:**
1. **What cancer is, genomically** (~25 min). Hallmarks of cancer; driver vs passenger mutations; clonal evolution; tumor heterogeneity. The cancer genome as accumulated damage.
2. **The cancer sequencing landscape** (~30 min). Tumor-only vs tumor-normal; FFPE artifacts; ctDNA / liquid biopsy; clinical panels (FoundationOne, MSK-IMPACT) vs WES vs WGS.
3. **Somatic variant calling at depth** (~40 min). Mutect2/Strelka2 revisited from L4 with cancer-specific concerns: low VAF subclones, contamination, FFPE deamination signatures. Mutational signatures (COSMIC SBS catalog) — what they reveal about etiology.
4. **Structural variants and fusions** (~30 min). Cancer SVs from short reads (Manta, GRIDSS) and long reads (Sniffles, Severus). Gene fusions from RNA-seq (Arriba, STAR-Fusion); BCR-ABL, EML4-ALK as canonical examples.
5. **Tumor heterogeneity and clonal evolution** (~30 min). Subclonal reconstruction (PyClone, SciClone); phylogenetic trees of tumor evolution (LICHeE, MEDICC2); single-cell views (callback to L8).
6. **Clinical interpretation in oncology** (~30 min). AMP/ASCO/CAP tier system (different from L17 ACMG/AMP — somatic-specific); OncoKB, CIViC; tumor mutational burden; microsatellite instability; HRD scoring.
7. **From sequencing to therapy** (~25 min). Targeted therapies and biomarkers; immunotherapy and TMB/MSI; PARP inhibitors and HRD; the actual decision flow in a tumor board.
8. **Wrap-up: what the field looks like in industry** (~20 min). Where graduates work (clinical labs, pharma, biotech, sequencing companies, academic cores); the role of bioinformatics engineers vs computational biologists vs data scientists; final pointers to keep learning.

**Tools introduced**: Mutect2 (revisited), Strelka2, Manta (revisited), GRIDSS (revisited), PyClone, Arriba, STAR-Fusion, SigProfiler, OncoKB, CIViC.

**Capstone exercise**: by this point students will have built a Nextflow pipeline through homework that handles the FASTQ → annotated VCF flow. Final project extends it to a tumor-normal somatic pipeline run on a small public cancer dataset.

---

## Drafting order recommendation

Draft in lecture order (L5 → L18) rather than out of order. The callback threads matter, and drafting L11 before L5 would mean rewriting L11 once L5's framing settles. Each lecture should be drafted at the same level of detail as L2 (with figure callouts, artifact callouts, EE framings, and section timings).

## Open items for later

- **Homework / exercise plan**: Nextflow taught through homework, culminating in the L18 capstone project. To be designed once lectures are drafted.
- **Project tracks**: should there be alternative final-project tracks for students with different interests (e.g. a single-cell track, a structural biology track)? Decide once lectures land.
- **Slide / figure / artifact production**: each lecture needs SVGs and interactive HTML artifacts at the same density as L2. Scope and tooling decision pending.
- **Assessment**: graded homework, exam, final project weighting. To be designed.
