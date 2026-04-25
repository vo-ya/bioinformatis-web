# Lecture 23 (proposed L18) — Metagenomics and the Microbiome

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology assumed
> **File**: provisional `lectures/lecture-23.html` — renumber to `lecture-18.html` when curriculum is reordered.

> **Proposed placement**: insert after L12 (population genetics, becomes new L17) — natural fit because microbial communities are populations, with allele frequencies, drift, selection, and admixture all transferring directly. The metagenome is a population assayed by sequencing rather than by sampling individuals. Placement before GWAS lets the population-genetics machinery generalise from one species (humans) to many (microbiomes) before pivoting to disease association.

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Distinguish 16S rRNA amplicon sequencing from shotgun metagenomics; describe the data each produces and when each is appropriate.
2. Run a 16S workflow with DADA2 / QIIME2: from FASTQ → ASVs (amplicon sequence variants) → taxonomy → diversity statistics.
3. Run a shotgun metagenomics workflow: read QC → taxonomic classification (Kraken2, MetaPhlAn) → functional profiling (HUMAnN) → MAG (metagenome-assembled genome) recovery.
4. Compute and interpret alpha diversity (within-sample) and beta diversity (between-sample) metrics; explain Shannon, Simpson, Bray-Curtis, UniFrac as information-theoretic / phylogenetic measures.
5. Apply differential-abundance testing (ANCOM-BC, MaAsLin2, songbird) with appropriate compositional data corrections.
6. Describe the major host-microbiome connections: gut microbiome and IBD, immunotherapy response, metabolomics, oral microbiome and cancer.
7. Frame metagenomics as detection-and-classification under compositional constraint: amplicon analysis as fingerprint-matched filter, taxonomic classification as approximate-nearest-neighbour search, diversity as information-theoretic divergence.

---

## Part 1 — Why Microbes (≈25 min)

### 1.1 The microbial scale (≈5 min)

The human body harbours ~10¹³ microbial cells (about as many as human cells), encoding ~3 million microbial genes (vs ~25k human genes). 99% of these microbes are bacteria; the rest are archaea, fungi, viruses, and protozoans.

Most microbiome research focuses on:

- **Gut**: largest community (~10¹⁴ cells/g feces), most clinical impact.
- **Oral**: second-largest, links to cardiovascular and pancreatic cancer.
- **Skin**: niche-specific (sebaceous, moist, dry).
- **Vaginal**: dominated by *Lactobacillus*; clinical relevance for preterm birth.
- **Respiratory**: small but important in asthma and COPD.

### 1.2 What we want to know (≈5 min)

Three core questions:

1. **Who's there?** — taxonomic composition (which species, at what abundance).
2. **What can they do?** — functional potential (which metabolic pathways, virulence factors, antibiotic resistance).
3. **What does that mean for the host?** — phenotype association (disease, drug response, diet effect).

Each maps to a different sequencing strategy + analysis pipeline.

### 1.3 The two main assays (≈8 min)

**16S rRNA amplicon sequencing**:

- PCR-amplify a hypervariable region of the bacterial 16S rRNA gene (V3–V4 typical).
- Sequence the amplicon (~250 bp paired-end Illumina, ~$30/sample).
- Cluster reads to taxonomic units (ASVs); assign taxonomy via reference DB (SILVA, GTDB, Greengenes).
- **Advantages**: cheap, well-validated.
- **Limitations**: bacteria only (not archaea/fungi/viruses); genus-level resolution typically; primer biases.

**Shotgun metagenomics**:

- Random sequencing of all DNA in the sample (~5 Gb/sample, ~$300).
- No PCR amplification of marker genes.
- Direct taxonomic classification of reads (Kraken2, MetaPhlAn) + functional profiling (HUMAnN).
- **Advantages**: species- to strain-level resolution; covers all kingdoms; functional content; viable for MAG assembly.
- **Limitations**: 10× cost; bioinformatics-heavy; host-DNA contamination in low-biomass samples.

A hybrid pattern: 16S for screening hundreds of samples, shotgun on selected subset for deep characterisation.

### 1.4 The microbial population analogy (≈5 min)

Treat a microbiome like a population from L12:

- **Allele frequency** ↔ taxon relative abundance.
- **Effective population size** ↔ community richness.
- **Drift** ↔ stochastic abundance fluctuations.
- **Selection** ↔ host or environmental selection (diet, antibiotics).
- **Migration** ↔ external inoculation (FMT, diet, contact).

Population-genetics metrics (Fst, Nei's distance) generalise to microbial-community comparison via UniFrac and Bray-Curtis.

### 1.5 The deep dive (≈2 min)

> **EE framing — metagenomics as community signal demixing**: A metagenomic sample is a **mixture** of DNA from hundreds of organisms. Recovering taxon abundances is a **demixing problem** — the same family as ICA / NMF / source separation that has come up repeatedly in this course (mutational signatures L18, scRNA-seq L8). The reference-genome database is the **dictionary**; reads are the **mixture**; abundances are the **mixing coefficients**. Modern classifiers (Kraken2, Bracken) solve this approximately via k-mer matching; statistical methods (Salmon-meta) solve it via expectation-maximisation on read assignments. The pattern is identical.

**FIGURE — Figure #1: 16S vs shotgun metagenomics decision** → `diagrams/lecture-23/01-16s-vs-shotgun.svg`

---

## Part 2 — 16S Amplicon Workflow (≈40 min)

### 2.1 The 16S rRNA gene (≈4 min)

Bacterial 16S rRNA is a ~1500 bp ribosomal-RNA gene with **alternating conserved and variable regions** (V1 through V9). PCR primers anneal to conserved regions; the variable region between is amplified.

V3–V4 (~460 bp) is the most-used region — long enough for genus-level resolution, fits Illumina 2×250 bp paired-end.

### 2.2 The DADA2 pipeline (≈10 min)

**DADA2** (Callahan et al. 2016) is the modern 16S workflow:

1. Quality-filter and trim reads.
2. **Learn error rates**: model per-sample sequencing error using the Illumina error profile.
3. **Denoise**: cluster reads into **Amplicon Sequence Variants (ASVs)** — exact sequences distinguishable from sequencing error. ASVs replace OTUs (Operational Taxonomic Units) which clustered at 97% identity.
4. **Merge** paired-end reads.
5. **Remove chimeras** (PCR artefacts that combine two parent sequences).
6. **Assign taxonomy** via naive Bayes classifier against SILVA/GTDB/Greengenes.

ASVs are exact sequences, not 97%-clusters — DADA2 resolves single-base differences that legacy OTU methods missed.

### 2.3 QIIME2 (≈4 min)

**QIIME2** is the dominant umbrella platform: end-to-end pipeline with DADA2 as the default denoiser, plus diversity, differential-abundance, and visualisation modules. Integrates with the Galaxy web interface for non-bioinformaticians.

### 2.4 Reference databases (≈4 min)

- **SILVA**: most comprehensive, regularly updated.
- **GTDB**: phylogenetically-coherent taxonomy (no paraphyletic groups).
- **Greengenes 13_8**: legacy default, frozen at 2013.

For new analyses use SILVA (default) or GTDB (for phylogenetic rigour).

### 2.5 Outputs and downstream analysis (≈8 min)

DADA2/QIIME2 produces:

- **ASV table**: rows = samples, columns = ASVs, cells = read counts.
- **Taxonomy table**: ASV → kingdom/phylum/class/order/family/genus/species assignment.
- **Phylogenetic tree** of ASVs (built via FastTree).
- **Sample metadata** linkage.

This is the equivalent of a count matrix + annotation in scRNA-seq (L7) — the same downstream tools (PCA, clustering, differential abundance) apply with appropriate statistical adjustments.

### 2.6 The 2024 frontier (≈5 min)

- **Long-read 16S**: full-length 16S via PacBio HiFi gives strain-level resolution (vs genus-level from short-read V3–V4).
- **Multi-region 16S**: V1V9 amplicons covering all variable regions improve taxonomic resolution.
- **DADA2 + AlphaFold-derived priors**: emerging, marginal improvement.

### 2.7 Worked example (≈5 min)

**FIGURE — Figure #2: DADA2 16S workflow** → `diagrams/lecture-23/02-dada2.svg`

---

## Part 3 — Shotgun Metagenomics (≈40 min)

### 3.1 The shotgun design (≈4 min)

Shotgun metagenomics: extract total DNA, fragment it, sequence randomly. No PCR step. Reads come from any organism in the sample, plus host contamination if present.

Workflow:

1. Read QC (fastp, BBDuk).
2. Host-DNA removal (map to human reference, discard mappers).
3. Taxonomic classification.
4. Functional profiling.
5. Optional: MAG assembly.

### 3.2 Taxonomic classification: Kraken2 (≈10 min)

**Kraken2** (Wood, Lu, Langmead 2019): k-mer-based exact classifier.

- Build a database of k-mers (k=35 default) from all reference genomes (NCBI RefSeq + GenBank).
- For each read, find all matching k-mers in the database.
- Assign read to the **lowest common ancestor (LCA)** of all matched k-mers' taxonomic origins.
- Run time: ~minutes per sample on a typical compute node.

**Bracken** post-processing converts Kraken2 read counts to Bayesian-corrected abundance estimates.

### 3.3 Taxonomic classification: MetaPhlAn (≈6 min)

**MetaPhlAn** (Beghini et al. 2021): marker-gene-based classifier.

- Pre-computed database of clade-specific marker genes (~5 markers per species).
- Reads aligned (Bowtie2) to markers; relative abundance inferred from marker coverage.
- More accurate than Kraken2 at species level for well-characterised organisms.
- Slower (Bowtie2 alignment) but more interpretable.

For a typical workflow: Kraken2 for breadth (catch everything), MetaPhlAn for clinical-grade species calls.

### 3.4 Functional profiling: HUMAnN (≈6 min)

**HUMAnN** (HMP Unified Metabolic Analysis Network, Franzosa 2018): map reads to UniRef90 + ChocoPhlAn (per-species pangenome) → quantify pathway abundance.

Output:

- Per-species, per-gene-family abundances.
- Pathway-level abundances (MetaCyc).
- Stratified output: each pathway abundance broken down by contributing species.

This is the **functional analog of expression analysis** — instead of "which genes are expressed?", it's "which microbial functions are present?".

### 3.5 MAG recovery (≈8 min)

A **Metagenome-Assembled Genome (MAG)** is a near-complete bacterial genome assembled directly from metagenomic reads — no isolation required.

Pipeline:

1. Assemble reads (megahit, SPAdes-meta).
2. **Binning**: cluster contigs by composition (k-mer spectra) + coverage (depth across samples) using MetaBAT2, MaxBin, CONCOCT.
3. Quality assessment (CheckM): completeness + contamination.
4. Taxonomy (GTDB-Tk).

A "high-quality MAG" has > 90% completeness and < 5% contamination. The **Genomes from Earth's Microbiomes (GEM)** catalogue (Nayfach et al. 2021) contains ~50,000 MAGs from ~10,000 environmental samples.

### 3.6 Worked example flow (≈3 min)

**FIGURE — Figure #3: Shotgun metagenomics pipeline** → `diagrams/lecture-23/03-shotgun-pipeline.svg`

> **EE framing — taxonomic classification as approximate-nearest-neighbour search**: Kraken2's k-mer LCA is approximate-nearest-neighbour search in $4^k$-dimensional Hamming space. The reference database is the indexed dictionary; query reads are decomposed into k-mers; the nearest neighbour gives taxonomic identity. Hashing-based exact ANN, the same family as locality-sensitive hashing in EE / signal processing. The LCA aggregation handles ambiguity gracefully — when k-mers match multiple genomes, fall back to common ancestor rather than picking arbitrarily.

---

## Part 4 — Diversity Metrics (≈25 min)

### 4.1 Alpha diversity (≈8 min)

**Alpha diversity** = within-sample diversity. Three flavours:

- **Richness**: number of distinct taxa observed. Easy to compute, sensitive to sequencing depth.
- **Shannon entropy**: $H = -\sum_i p_i \log p_i$. Information-theoretic. Combines richness and evenness.
- **Simpson's index**: $1 - \sum_i p_i^2$. Probability that two reads from the sample are different taxa.

For statistical inference, normalise depth by rarefaction (subsample all samples to same depth) or use depth-adjusted metrics (Faith's PD).

### 4.2 Beta diversity (≈8 min)

**Beta diversity** = between-sample diversity:

- **Bray-Curtis dissimilarity**: $\frac{\sum_i |a_i - b_i|}{\sum_i (a_i + b_i)}$. Compositional difference.
- **Jaccard**: presence/absence-based.
- **UniFrac (unweighted / weighted)**: phylogenetic — accounts for tree distance between observed taxa.

Beta-diversity matrices feed PCoA / PCA visualisations and PERMANOVA hypothesis testing.

### 4.3 The compositional data problem (≈5 min)

**Critical**: microbial abundances are **relative**, not absolute. The total DNA in a sample is normalised by sequencing depth. If species A doubles in absolute abundance, but species B triples, A's relative abundance **decreases** even though A is "growing".

**Compositional data analysis (CoDa)**:

- Aitchison transforms (ALR, CLR, ILR).
- Constraint: components sum to a constant.
- Many "differential abundance" results in the early literature were artefacts of compositional bias.

Modern tools (ANCOM-BC, songbird) explicitly account for compositionality.

### 4.4 Differential abundance: ANCOM-BC (≈4 min)

**ANCOM-BC** (Lin & Peddada 2020): differential-abundance test under compositional constraint.

- Models log-transformed abundances with sample-specific bias.
- Multiple-testing-corrected via Benjamini-Hochberg.
- Recommended over t-test, DESeq2, edgeR for microbiome data.

**FIGURE — Figure #4: Diversity metrics — alpha and beta** → `diagrams/lecture-23/04-diversity.svg`

---

## Part 5 — Host-Microbiome Connections (≈30 min)

### 5.1 Gut microbiome and IBD (≈6 min)

Inflammatory Bowel Disease (Crohn's, ulcerative colitis) shows reproducible gut-microbiome signatures:

- Reduced diversity (Shannon ~ 1.5 lower than healthy).
- Loss of obligate anaerobes (e.g., *Faecalibacterium prausnitzii*).
- Increase of pathobionts (*Escherichia*, *Klebsiella*).

Diagnostic biomarker development is active; the **MetaCardis** and **iHMP-IBD** consortia have produced ~5000-sample datasets that enable machine-learning classifiers (~85% accuracy IBD vs control).

### 5.2 Immunotherapy response (≈8 min)

Cancer-immunotherapy response (anti-PD-1) correlates with gut microbiome composition. Key findings:

- *Akkermansia muciniphila* abundance predicts response in multiple solid tumours (Routy et al. 2018).
- *Bacteroides fragilis* enrichment in non-responders.
- Antibiotics shortly before checkpoint blockade reduce response rate.

Several Phase I/II trials are testing **fecal microbiota transplantation (FMT)** as adjuvant to immunotherapy. ~10 small studies show response improvement in previously refractory patients.

### 5.3 Other clinical connections (≈6 min)

- **Oral microbiome → pancreatic cancer**: *Porphyromonas gingivalis* abundance is elevated in PDAC patients.
- **Colorectal cancer**: *Fusobacterium nucleatum* enrichment in tumour tissue (Castellarin 2012).
- **Type 2 diabetes**: gut microbiome shifts predict response to metformin.
- **Drug metabolism**: bacterial enzymes activate / deactivate ~30% of orally-administered drugs (Spanogiannopoulos 2016).

### 5.4 Cohort-scale population studies (≈4 min)

Major resources:

- **HMP / iHMP** (Human Microbiome Project): 300 + 1700 individuals, multi-site.
- **MetaCardis**: 2000 individuals, cardiometabolic disease focus.
- **TwinsUK Microbiome**: heritability estimation.
- **Fragile Microbiome Atlas**: longitudinal sampling in pregnancy.

Statistical scale: ~200,000 microbiome samples publicly available.

### 5.5 The Hippocrates 2.0 question (≈2 min)

Treating disease via microbiome modulation:

- Direct: FMT, probiotics, prebiotics.
- Indirect: dietary intervention, drug repurposing.
- Engineered: synthetic microbes (live biotherapeutics) currently in Phase II–III for *C. difficile*, IBD, allergy.

The 2030s answer to "what does microbiome do" is "it's a manipulable disease-modifying organ".

### 5.6 Worked example (≈4 min)

**FIGURE — Figure #5: Microbiome-disease association case study** → `diagrams/lecture-23/05-microbiome-disease.svg`

**EMBED — Artifact #5: Microbiome Differential Abundance Explorer** → `artifacts/lecture-23/05-diff-abundance.html`

---

## Part 6 — Tools, Pitfalls, and Practice (≈25 min)

### 6.1 The standard 2024 stack (≈6 min)

For 16S workflows:

- **QIIME2** (full pipeline) or **DADA2 in R**.
- Diversity + ordinations via QIIME2 or `phyloseq` (R) or `scikit-bio` (Python).
- Differential abundance: ANCOM-BC, MaAsLin2, songbird.

For shotgun:

- **bioBakery 4** stack: KneadData (host removal) → MetaPhlAn4 → HUMAnN3.
- Or **Sunbeam / nf-core/mag** for Nextflow-based reproducible workflows.
- MAG recovery: nf-core/mag, Atlas, ATLAS.

### 6.2 Pitfalls (≈8 min)

**Sample collection**:

- Storage temperature matters: -80°C ideal; lyophilisation acceptable.
- Stabilising buffers (RNAlater, OMNIgene) for ambient shipping.
- Cross-contamination from sample handling.

**Sequencing**:

- Low-biomass samples (skin, low-density gut) → high host-DNA contamination.
- PCR amplification bias for 16S — primers preferentially amplify some taxa.
- Sequencing depth: 16S requires 10k reads/sample; shotgun ideally 5 Gb/sample.

**Bioinformatics**:

- ASV vs OTU debate: ASVs are objectively better; legacy data uses OTUs.
- Reference DB updates change taxonomic assignments — pin DB versions.
- Compositional bias as discussed.

**Statistical**:

- Sample-size requirements: for typical effect sizes, 30+ per group.
- Multiple-testing across thousands of taxa requires FDR control.
- Subject-level confounders (age, diet, antibiotics) often dominate the signal.

### 6.3 Reproducibility (≈4 min)

Microbiome reproducibility is famously poor (~30% of associations replicate across cohorts). Drivers:

- Cohort effects (geography, diet, ethnicity).
- Technical variation (DNA extraction, primer choice, sequencer).
- Statistical methodology (compositional vs non-compositional).

The **Strain Resolved Microbiome Wiki** and standardised pipeline efforts (cwl-microbiome) aim to improve reproducibility. Until 2030, expect signature replication around 50%.

### 6.4 The 2024 frontier (≈4 min)

- **Strain-level resolution**: PacBio HiFi for full-length 16S; long-read shotgun.
- **Multi-omics integration**: 16S + metagenomics + metabolomics + host transcriptomics.
- **Engineered microbiome**: live biotherapeutic products entering clinical trials.
- **Foundation models for microbiome**: protein-language-model-style embeddings for microbial proteins (ESM-2 has microbial sub-sets).

### 6.5 Hands-on (≈3 min)

A reasonable first project:

- Download a 16S dataset from MGnify (e.g., a 30-sample IBD cohort).
- Run QIIME2 pipeline.
- Compute alpha and beta diversity.
- Test for IBD vs control differential abundance with ANCOM-BC.
- Visualise with PCoA.

Total: ~1 day of compute, easily fits on a laptop.

**FIGURE — Figure #6: Microbiome analysis workflow** → `diagrams/lecture-23/06-workflow.svg`

---

## Part 7 — Frontier Topics (≈15 min)

### 7.1 Antibiotic resistance gene tracking (≈4 min)

Shotgun metagenomics can map antibiotic resistance genes (ARGs) in clinical samples and environmental reservoirs. Tools: CARD-RGI, AMRFinderPlus. Critical for tracking resistance spread between species (horizontal gene transfer in real time).

### 7.2 Phylogenetic placement (≈3 min)

For low-resolution amplicon data, **phylogenetic placement** (EPA-ng, pplacer) inserts queries into a reference tree without rebuilding it — useful for stable taxonomic comparison across studies.

### 7.3 Strain-level dynamics (≈3 min)

Within-species strain variation matters: different *E. coli* strains range from commensal to pathogen. **StrainPhlAn**, **MIDAS**, **inStrain** infer strain composition from shotgun data.

### 7.4 Microbiome-host genetic interaction (≈3 min)

GWAS for microbiome composition: **microbiome-GWAS (mGWAS)** identifies host genetic variants associated with microbial abundance. Examples: lactose-tolerance variants → *Bifidobacterium*; ABO blood group → *Bacteroides* abundance.

### 7.5 The 2024 frontier (≈2 min)

- **Microbiome foundation models** (Genomic-LM-style): pre-trained on millions of microbial genomes.
- **Synthetic ecology**: design defined consortia for therapeutic delivery.
- **Single-cell metagenomics**: Drop-seq for microbial cells.

---

## Wrap-up (≈10 min)

### What you should take away

- **Two main assays**: 16S (cheap, genus-level, bacteria-only); shotgun (expensive, species-strain-level, multi-kingdom).
- **DADA2 + QIIME2** is the 16S standard; **bioBakery + Kraken2/MetaPhlAn/HUMAnN** is the shotgun standard.
- **Diversity metrics** generalise from population genetics: Shannon entropy = information-theoretic richness; UniFrac = phylogenetic distance; Bray-Curtis = compositional dissimilarity.
- **Compositional bias is critical**: relative abundances are constrained; modern differential-abundance tools (ANCOM-BC) account for this.
- **Host-microbiome connections** are real but small; replication rates ~50%; gut-immunotherapy-response is the most clinically validated.
- **EE framings**: metagenomics as community signal demixing; taxonomic classification as approximate-nearest-neighbour search; diversity as information-theoretic divergence.

### Next lecture

GWAS + Mendelian randomisation (existing L13, becomes new L19). The microbiome population framework here connects to the host-genetics population framework there.

### Homework

1. Download a 16S amplicon dataset from MGnify or the SRA (any IBD or healthy-control cohort). Run DADA2; report ASV count, taxonomic assignments, alpha-diversity Shannon. ~3h compute.
2. From the same dataset, compute beta-diversity Bray-Curtis matrix. PCoA visualisation. Interpret the first two principal coordinates.
3. Run ANCOM-BC for differential-abundance testing case vs control. Report taxa significant at FDR < 0.05.
4. From a public shotgun dataset (e.g., HMP), run Kraken2 and Bracken; compare species call vs the published MetaPhlAn results. Where do they disagree?
5. For one MAG of interest, run CheckM. Report completeness + contamination. Annotate with GTDB-Tk; identify the closest reference.

### Recommended reading

- Callahan, B. J., McMurdie, P. J., Rosen, M. J., et al. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. *Nature Methods* 13, 581–583.
- Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology* 20, 257.
- Beghini, F., McIver, L. J., Blanco-Míguez, A., et al. (2021). Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife* 10, e65088.
- Nayfach, S., Roux, S., Seshadri, R., et al. (2021). A genomic catalog of Earth's microbiomes. *Nature Biotechnology* 39, 499–509.
- Routy, B., Le Chatelier, E., Derosa, L., et al. (2018). Gut microbiome influences efficacy of PD-1-based immunotherapy against epithelial tumors. *Science* 359, 91–97.
- Lin, H., & Peddada, S. D. (2020). Analysis of compositions of microbiomes with bias correction. *Nature Communications* 11, 3514. (ANCOM-BC.)
- Spanogiannopoulos, P., et al. (2016). The microbial pharmacists within us. *Nature Reviews Microbiology* 14, 273–287.
- QIIME2: <https://qiime2.org/>
- bioBakery: <https://huttenhower.sph.harvard.edu/biobakery>
- MGnify: <https://www.ebi.ac.uk/metagenomics/>
- HMP / iHMP: <https://www.ihmpdcc.org/>
- GTDB: <https://gtdb.ecogenomic.org/>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — Why Microbes                                  | 25 min | 0:25 |
| Part 2 — 16S Amplicon Workflow                          | 40 min | 1:05 |
| Part 3 — Shotgun Metagenomics                            | 40 min | 1:45 |
| Part 4 — Diversity Metrics                                | 25 min | 2:10 |
| Part 5 — Host-Microbiome Connections                       | 30 min | 2:40 |
| Part 6 — Tools, Pitfalls, and Practice                       | 25 min | 3:05 |
| Part 7 — Frontier Topics                                       | 15 min | 3:20 |
| Wrap-up                                                          | 10 min | 3:30 |

**Total:** ~3h 30min of content.
