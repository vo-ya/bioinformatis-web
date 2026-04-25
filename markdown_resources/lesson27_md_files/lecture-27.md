# Lecture 27 (proposed L11) — Mass-Spectrometry Proteomics Primer

> **Duration**: ≈3h 0min content
> **Audience**: EE undergraduates / graduates, minimal biology / chemistry assumed
> **File**: provisional `lectures/lecture-27.html` — renumber to `lecture-11.html` when curriculum is reordered.

> **Proposed placement**: insert after L6 (differential expression, becomes new L9) and before L7 (scRNA-seq, becomes new L11). Natural arc: bulk RNA-seq gives transcript abundances → DE analysis identifies regulated genes → MS proteomics gives the protein-level readout — often the more relevant biological output. Position complements L9 / new L13 (ChIP-seq) by giving students the protein-side measurement before single-cell methods. Lecture is ~3h (vs the standard 3.5h) — proteomics is broad and this is a primer, not a deep dive.

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Describe the principles of LC-MS/MS for proteomic measurement: protein digestion, peptide separation, ionisation, fragmentation, mass analysis.
2. Distinguish bottom-up, top-down, and middle-down proteomics; describe DDA (data-dependent) vs DIA (data-independent) acquisition.
3. Run a basic peptide identification: from raw MS spectra → database search (Mascot, MS-Fragger, MaxQuant) → protein inference → false-discovery-rate filtering.
4. Apply quantification methods: label-free, SILAC, TMT, iTRAQ; recognise their statistical properties and sample-size implications.
5. Distinguish proteomics from RNA-seq biologically: protein vs mRNA correlation, post-translational modifications, splice isoforms, secretome.
6. Describe key proteomics applications: biomarker discovery, drug target validation (chemoproteomics), structural proteomics (cross-linking, hydrogen-deuterium exchange).
7. Frame MS proteomics in EE terms: spectra as 1D signals; peptide identification as constrained search; quantification as signal-to-noise estimation; DIA as compressed sensing on a peptide-coverage matrix.

---

## Part 1 — From Mass Spectrum to Peptide (≈25 min)

### 1.1 Why proteins, not just transcripts (≈4 min)

RNA-seq measures **transcript abundance**. But protein abundance is what does the work in cells. The correlation between mRNA and protein abundance is surprisingly low (typically r = 0.4-0.6 across genes), driven by:

- **Translation efficiency**: some transcripts are heavily translated; others sit unused.
- **Post-translational modifications**: phosphorylation, acetylation, ubiquitination radically change protein function without changing transcript count.
- **Protein degradation**: some proteins are short-lived; others persist for days.
- **Secretion**: secreted proteins (cytokines, hormones) are biology-defining but not present in the cell that made them.
- **Splice isoforms** at protein level (which protein isoform is actually produced?).

Many drug targets and disease biomarkers are inferred from protein-level signals that RNA-seq alone misses.

### 1.2 The LC-MS/MS instrument (≈8 min)

The dominant proteomics workflow: **liquid chromatography coupled to tandem mass spectrometry**.

1. **Protein extraction**: from cell, tissue, plasma, etc.
2. **Tryptic digestion**: trypsin (a protease) cleaves proteins after K and R into peptides ~5-25 residues.
3. **LC separation**: peptides flowed through a reverse-phase chromatography column over ~1-3 hours; elute at characteristic times.
4. **Ionisation (ESI)**: peptides aerosolized + ionised in electrospray.
5. **MS1 scan**: parent peptide masses measured.
6. **MS/MS**: selected peptides fragmented (typically by HCD); fragment ions measured (MS2).
7. **Detector**: ion-counting, m/z resolved.

The result: thousands of MS2 spectra per sample, each a fragment-ion pattern from one peptide (typically).

### 1.3 The fragmentation spectrum (≈6 min)

In MS2, peptides fragment along the peptide backbone, producing **b-ions** (N-terminal fragments) and **y-ions** (C-terminal fragments).

For a peptide `LSDPYHRGSP`:

- y1: P (singleton C-terminal).
- y2: SP.
- y3: GSP.
- ...
- b1: L.
- b2: LS.
- ...

Each fragment has a calculable mass. The full b/y series gives the sequence.

A real MS2 spectrum has these masses + many low-intensity noise peaks. Identification = match the observed fragment masses to a candidate peptide's predicted spectrum.

### 1.4 Database search (≈4 min)

Given an observed MS2 spectrum and a protein database (e.g., UniProt human reference):

1. In silico digest all proteins → candidate peptides.
2. For each candidate peptide, predict the theoretical b/y fragment masses.
3. Score against observed spectrum (cross-correlation or related).
4. Best-scoring candidate = peptide identification.

Tools: **Mascot** (commercial, classic), **MaxQuant** (free, dominant), **MS-Fragger** (fast, modern), **Comet**.

### 1.5 The deep dive (≈3 min)

> **EE framing — peptide identification as constrained matched-filter search**: An MS2 spectrum is a 1D signal: m/z (x-axis) × intensity (y-axis). A candidate peptide's theoretical spectrum is a discrete set of expected peaks. Identification = matched-filter detection — does the candidate's predicted-peak pattern correlate with the observed signal? The constraint is **enzyme-specific**: only peptides ending in K or R (after trypsin) are candidates. Modern tools combine matched-filter scoring with statistical models (target-decoy FDR estimation).

**FIGURE — Figure #1: LC-MS/MS workflow** → `diagrams/lecture-27/01-lc-ms-ms.svg`

---

## Part 2 — Acquisition Modes (DDA vs DIA) (≈25 min)

### 2.1 Data-dependent acquisition (DDA) (≈6 min)

Classical mode:

- MS1 scan: identify the most abundant peptide masses.
- Pick the top-N (typically N = 10-20) most intense MS1 peaks.
- Fragment each in MS2 sequentially.
- Repeat throughout the LC gradient.

**Advantages**: well-validated, mature pipelines.
**Limitations**: stochastic — if a peptide isn't in the top-N at some moment, it's missed. Peptide coverage varies between samples.

### 2.2 Data-independent acquisition (DIA) (≈8 min)

Modern alternative:

- MS1 scan as before.
- Instead of picking specific peaks, **systematically fragment all peptides in defined m/z windows** (e.g., 20-Da windows from 400 to 1200 m/z = 40 windows).
- Each MS2 spectrum is a "co-fragmented mixture" — multiple peptides fragmenting simultaneously in the window.
- Identification: deconvolve the mixture using a peptide library or de novo.

**Advantages**: deterministic peptide coverage; reproducible across samples.
**Limitations**: chimeric spectra need deconvolution; needs a spectral library or DIA-NN (DIA-aware neural net).

DIA is becoming dominant for clinical / cohort-scale proteomics where reproducibility matters.

### 2.3 Targeted acquisition (≈4 min)

For absolute quantification of a few specific peptides:

- **MRM (Multiple Reaction Monitoring)**: triple-quadrupole MS; measure specific precursor → fragment pairs.
- **PRM (Parallel Reaction Monitoring)**: high-resolution variant.

Used for clinical biomarker quantification (e.g., troponin I, BNP) where precision matters more than coverage.

### 2.4 The deep dive (≈4 min)

> **EE framing — DIA as compressed sensing on peptides**: In DIA, each MS2 spectrum is a **superposition** of peptide fragmentation patterns. The acquisition matrix is sparse — only peptides in the m/z window contribute. Identification = solving a sparse-recovery problem: find the smallest set of library peptides whose theoretical spectra explain the observed spectrum. This is structurally identical to **compressed sensing**: a sparse signal, an under-determined linear system, recovered via L1 / sparse-coding approaches. Modern DIA tools (DIA-NN, Spectronaut, OpenSWATH) implement this with neural-network-trained spectral libraries.

### 2.5 Worked example (≈3 min)

**FIGURE — Figure #2: DDA vs DIA acquisition** → `diagrams/lecture-27/02-dda-vs-dia.svg`

---

## Part 3 — Identification and FDR (≈30 min)

### 3.1 Target-decoy strategy (≈5 min)

How do we know how many false identifications we have?

**Target-decoy approach**: search not only against the real protein database (target), but also against a reverse / shuffled database (decoy). Decoy hits are by definition false. Compute FDR = decoy hits / target hits at any score threshold.

Classical threshold: **FDR < 1%** at the peptide level.

This generalises to FDR-based ML scorers like Percolator (Käll 2007) — train a classifier on target/decoy features; use predicted score for ranking.

### 3.2 Protein inference (≈6 min)

Peptides identify proteins, but:

- A peptide can be shared across multiple proteins (peptide ambiguity).
- A protein with no unique peptides cannot be unambiguously identified.

**Razor peptide assignment**: shared peptides assigned to the protein with most unique peptides.

**Protein FDR**: separately controlled (typically also < 1%); often more stringent than peptide FDR.

For typical workflows: report peptide-level FDR + protein-level FDR.

### 3.3 Modifications and PTMs (≈8 min)

Protein post-translational modifications (PTMs) change peptide masses:

- **Phosphorylation** (phospho-S, T, Y): +79.97 Da.
- **Acetylation** (acetyl-K, N-terminus): +42.01 Da.
- **Methylation**: +14.02 Da per methyl.
- **Ubiquitination**: complex, GG remnant after trypsin = +114.04 Da on K.
- **Glycosylation**: complex, multi-glycoform.

Database search with **variable modifications** allows modified peptides to match. But each modification doubles the search space; many mods → combinatorial explosion.

Modern tools (MaxQuant, FragPipe) handle ~10 common variable modifications routinely.

### 3.4 Specialised modifications: phosphoproteomics (≈5 min)

Phosphorylation regulates kinase signalling. **Phosphoproteomics**:

- Enrich phosphopeptides via TiO2 or IMAC affinity columns.
- LC-MS/MS as standard.
- Identify phosphosites + their stoichiometry.

Output: ~20,000 phosphosites in a typical experiment. Bridges with kinase signalling networks (L22 / new L10) and CRISPR screens (L24 / new L21).

### 3.5 Worked example (≈3 min)

**FIGURE — Figure #3: MS2 spectrum identification + FDR** → `diagrams/lecture-27/03-ms2-id.svg`

### 3.6 The deep dive (≈3 min)

> **EE framing — PTMs as discrete random codes added to a base signal**: A peptide's mass spectrum is a base signal; a PTM adds a fixed offset to specific b-/y-fragments. Searching for PTMs = expanded matched-filter dictionary that includes modified variants. The combinatorial explosion (PTM combinations on multi-residue peptides) is exactly the **codebook expansion** problem in EE detection theory — search depth grows with the modification dictionary size.

---

## Part 4 — Quantification (≈25 min)

### 4.1 Label-free quantification (LFQ) (≈6 min)

Without isotopic labels, quantify by MS1 ion intensity:

- For each peptide, integrate its MS1 elution profile across the LC-MS run.
- Compare across samples.
- Aggregate per-protein.

**MaxQuant LFQ** (Cox 2014): the standard.

**Advantages**: cheap, widely-applicable, scalable.
**Limitations**: ~10-15% CV; missing values across samples (peptide not detected in all runs).

### 4.2 Stable isotope labelling (≈8 min)

**SILAC (Stable Isotope Labeling by Amino acids in Cell culture)** (Mann 2002):

- Grow cells in media with isotopically labelled amino acids (e.g., heavy K, R).
- Mix with control (light) sample.
- LC-MS/MS — peptide pairs differ in mass by predictable Δ.
- Quantify ratio = relative abundance.

**Advantages**: very accurate (~5% CV); paired comparison eliminates batch effects.
**Limitations**: limited to cell culture (can't do tissue, plasma).

### 4.3 Tandem Mass Tagging (TMT) and iTRAQ (≈6 min)

Modern multiplexed quantification:

- **TMT (Tandem Mass Tag)**: chemical labels (typically 6-, 11-, 16-plex) added to peptides post-digestion.
- All samples mixed; LC-MS/MS as one.
- MS2 spectra contain reporter ions at characteristic m/z; intensity = relative abundance.

**Advantages**: 6-16 samples in one run; reduced LC-MS time per sample.
**Limitations**: ratio compression (signal smearing); reagent cost.

For typical cohort proteomics: TMT-16plex is the standard. ~$5,000 for reagents per 16 samples.

### 4.4 The DIA quantification advantage (≈3 min)

DIA's deterministic peptide coverage gives much better LFQ-style quantification than DDA — typically 5% CV vs 15% for DDA-LFQ. For cohort work where reproducibility matters, DIA is now preferred.

### 4.5 Worked example (≈2 min)

**FIGURE — Figure #4: Quantification methods comparison** → `diagrams/lecture-27/04-quantification.svg`

---

## Part 5 — Statistical Analysis and Differential Abundance (≈25 min)

### 5.1 The protein-level matrix (≈3 min)

Output of a proteomics experiment: protein × sample matrix of relative abundances. Typical scale: ~5,000-8,000 proteins × 30-100 samples.

This is the proteomics analog of the L8/L9 RNA-seq count matrix.

### 5.2 Differential expression analysis (≈8 min)

Adapted from RNA-seq:

- **t-test** (with multiple testing correction): the simplest.
- **limma** (Smyth 2004): linear models with empirical-Bayes shrinkage. Originally for microarrays; adapted for proteomics.
- **DEqMS** (Zhu 2020): proteomics-specific extension of limma, handling missing values.
- **ProDA** (Ahlmann-Eltze 2020): probabilistic dropout-aware proteomics differential expression.

For a typical proteomics study: limma or DEqMS at FDR < 0.05.

### 5.3 Missing values: a proteomics-specific challenge (≈6 min)

Proteomics data has lots of missing values:

- A peptide not detected in some samples might be **missing-not-at-random** (MNAR) — really not present.
- Or **missing-at-random** (MAR) — present but missed by stochastic acquisition.

Imputation strategies:

- **Half-min**: impute as half the lowest observed value (assumes MNAR).
- **KNN**: nearest-neighbour imputation (assumes MAR).
- **Multiple imputation**: average across multiple imputation runs.
- **DIA-based**: simply skip imputation; DIA's coverage is consistent.

The choice matters: MNAR imputation creates bias if the missing actually are MAR.

### 5.4 Pathway analysis (≈5 min)

Same machinery as RNA-seq (L22 / new L10):

- ORA / GSEA on protein abundance changes.
- KEGG / Reactome / GO enrichment.
- Network propagation (STRING-based).

A typical proteomics paper reports DE proteins + pathway enrichment + (often) PPI network analysis.

### 5.5 Worked example (≈3 min)

**EMBED — Artifact #5: Proteomics differential abundance explorer** → `artifacts/lecture-27/05-diff-abundance.html`

---

## Part 6 — Applications (≈20 min)

### 6.1 Biomarker discovery (≈5 min)

Plasma / serum proteomics for clinical biomarkers:

- **Olink** and **SomaScan**: ~5000-protein affinity-based panels — measure protein abundances per blood draw.
- **Simoa**: single-molecule sensitivity for low-abundance proteins.
- **Targeted LC-MS/MS** (PRM): for absolute quantification of validated biomarkers.

Examples:

- **Cardiac biomarkers**: troponin I, BNP, ST2.
- **Cancer biomarkers**: PSA (prostate), CA-125 (ovarian), CA 19-9 (pancreatic).
- **Alzheimer's**: amyloid-β, tau, neurofilament-light.

The plasma proteome has ~5,000 detectable proteins; ~50 are widely-used clinical biomarkers.

### 6.2 Chemoproteomics and target validation (≈5 min)

**Chemoproteomics**: identify the proteins a drug binds in cells.

Workflow:

- Drug + photoaffinity / chemical-probe linker.
- Apply to live cells.
- Pull down crosslinked proteins.
- Identify by MS.
- → list of drug targets (and off-targets).

Tools: **CETSA** (Cellular Thermal Shift Assay), **TPP** (Thermal Proteome Profiling), activity-based protein profiling (ABPP).

This is the experimental complement to AlphaFold-based docking — find what a drug actually binds in cells, not what we predict it should bind.

### 6.3 Structural proteomics (≈5 min)

Identify protein-protein and protein-ligand interactions at near-atomic resolution:

- **Cross-linking MS (XL-MS)**: chemical crosslinker captures interacting residues; MS identifies them.
- **Hydrogen-deuterium exchange MS (HDX-MS)**: probe protein conformation changes upon ligand binding.
- **Native MS**: measure intact protein complexes' masses.

These complement AlphaFold + cryo-EM for understanding drug-target structural details.

### 6.4 Single-cell proteomics (≈3 min)

Emerging frontier (scProteomics):

- **CITE-seq** (with antibodies): not strictly MS, but gives ~100 proteins + scRNA-seq.
- **PiMMS / TRiC**: emerging single-cell MS workflows.
- ~$1k/sample currently; expected to drop with throughput improvements.

Currently mostly research; clinical adoption is years away.

### 6.5 The 2024 frontier (≈2 min)

- **PiMMS / Sonar**: single-cell proteomics via miniaturised LC-MS.
- **Spatial proteomics**: imaging mass cytometry, IMC, Codex; ~50 proteins per cell × spatial context.
- **Proteoforms**: characterising specific protein isoforms (different splice + PTM combinations).
- **PROTAC effects**: chemoproteomics for proximity-inducing degraders.

**FIGURE — Figure #5: Plasma proteome biomarker landscape** → `diagrams/lecture-27/05-biomarkers.svg`
**FIGURE — Figure #6: Chemoproteomics workflow** → `diagrams/lecture-27/06-chemoproteomics.svg`

---

## Part 7 — Tools, Pitfalls, and Workflows (≈15 min)

### 7.1 Standard 2024 stack (≈4 min)

For a typical proteomics analysis:

- **Search engine**: MaxQuant (DDA, free), MS-Fragger / FragPipe (fast DDA, free), DIA-NN (DIA), Spectronaut (DIA, commercial).
- **Quantification**: built-in to search engine.
- **Statistical analysis**: limma (R), DEqMS, MSstats.
- **Visualisation**: Perseus (the MaxQuant ecosystem; visual + stats).

For a clinical Olink / SomaScan study: vendor-specific QC, then standard limma / nlme.

### 7.2 Pitfalls (≈5 min)

**Sample preparation**:

- Contaminants (keratin from skin, blood from sample handling): trace amounts dominate.
- Tryptic digestion variability: incomplete digest produces missed-cleavage peptides.

**Acquisition**:

- DDA stochastic coverage: not all peptides identified in all samples.
- Mass accuracy drift: instrument calibration matters; weekly QC.
- Carryover between samples: solvent gradients matter.

**Identification**:

- Wrong PTM searches → many false-positive PTM IDs.
- Database mismatch (wrong species / version) → systematic mis-identification.
- Multi-protein peptide ambiguity: razor assignment is heuristic.

**Quantification**:

- Missing values: handle explicitly; don't naively impute.
- Ratio compression in TMT: real differences are dampened; fold-change estimates biased.
- Batch effects: even within-experiment, run-to-run variation matters.

### 7.3 Reproducibility (≈3 min)

Proteomics reproducibility is moderate:

- DIA workflows: ~80-85% protein identification overlap across replicates.
- DDA: ~70% identification overlap.
- DEG (differential expression of proteins) ~60-70% replication across cohorts.

Better than microbiome (~50%); worse than RNA-seq (~85%).

### 7.4 The 2024 frontier (≈3 min)

- **DIA + ML libraries**: DIA-NN's neural-net-trained spectral libraries enable database-free DIA analysis.
- **AlphaPeptDeep / Pulsar**: deep-learning-based spectrum prediction for closer matched-filter detection.
- **Foundation models for proteomics**: Pre-trained on millions of MS2 spectra.

**FIGURE — Figure #7: Proteomics analysis workflow** → `diagrams/lecture-27/07-workflow.svg`

---

## Wrap-up (≈10 min)

### What you should take away

- **MS proteomics measures protein abundance**: the post-transcriptional readout. Correlates only modestly with mRNA (~r = 0.5).
- **LC-MS/MS** (bottom-up): protein → trypsin → peptides → LC separation → MS1 / MS2 → identification → quantification.
- **DDA vs DIA**: DDA is stochastic; DIA is deterministic. DIA is becoming the cohort-scale standard.
- **Identification by database search** + target-decoy FDR. Modern tools: MaxQuant, MS-Fragger / FragPipe.
- **Quantification**: LFQ (cheap), SILAC (cell culture), TMT (multiplexed). Choose by question + budget.
- **Applications**: biomarker discovery (Olink, SomaScan, plasma MS), chemoproteomics (drug target validation), structural proteomics (XL-MS, HDX-MS).
- **EE framings**: spectra as 1D signals; peptide identification as constrained matched-filter search; PTMs as codebook expansion; DIA as compressed sensing on peptide-coverage matrix.

### Next lecture

Single-cell RNA-seq fundamentals (existing L7, becomes new L11 in proposed). The protein-level signals you learned here complement single-cell mRNA — modern multimodal experiments combine both.

### Homework

1. Run MaxQuant on a public proteomics dataset (e.g., an HEK293 LC-MS/MS run). Identify peptides at 1% FDR. Report protein count + identification overlap with the published list.
2. Pick a specific PTM (phosphorylation, acetylation). Run a dedicated search; identify modified peptides. Compare to a non-PTM search.
3. Implement label-free quantification on simulated MS1 elution profiles. Add 30% missing values; compare imputation strategies (half-min, KNN, no imputation).
4. From an Olink panel of 100 plasma proteins for a cohort (use a public Bioshare dataset), run differential abundance test for case vs control. Interpret the hits in terms of biology.
5. Walk through a published chemoproteomics study (e.g., kinase-inhibitor target identification via Kinobeads). Re-derive the target identification logic.

### Recommended reading

- Aebersold, R., & Mann, M. (2003). Mass spectrometry-based proteomics. *Nature* 422, 198–207.
- Cox, J., et al. (2014). Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. *Molecular & Cellular Proteomics* 13, 2513–2526.
- Tyanova, S., et al. (2016). The MaxQuant computational platform for mass spectrometry-based shotgun proteomics. *Nature Protocols* 11, 2301–2319.
- Käll, L., et al. (2007). Semi-supervised learning for peptide identification from shotgun proteomics datasets. *Nature Methods* 4, 923–925. (Percolator.)
- Demichev, V., et al. (2020). DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. *Nature Methods* 17, 41–44.
- Mann, M. (2002). The rise of mass spectrometry and the fall of Edman degradation. *Trends in Biochemical Sciences* 27, 233–235. (SILAC origins.)
- Thompson, A., et al. (2003). Tandem Mass Tags: A novel quantification strategy for comparative analysis of complex protein mixtures. *Analytical Chemistry* 75, 1895–1904.
- Savitski, M. M., et al. (2014). Tracking cancer drugs in living cells by thermal profiling of the proteome. *Science* 346, 1255784. (TPP.)
- MaxQuant: <https://www.maxquant.org/>
- DIA-NN: <https://github.com/vdemichev/DiaNN>
- Olink: <https://www.olink.com/>
- ProteomicsDB: <https://www.proteomicsdb.org/>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — From Mass Spectrum to Peptide               | 25 min | 0:25 |
| Part 2 — Acquisition Modes (DDA vs DIA)                | 25 min | 0:50 |
| Part 3 — Identification and FDR                          | 30 min | 1:20 |
| Part 4 — Quantification                                   | 25 min | 1:45 |
| Part 5 — Statistical Analysis and Differential Abundance | 25 min | 2:10 |
| Part 6 — Applications                                       | 20 min | 2:30 |
| Part 7 — Tools, Pitfalls, and Workflows                     | 15 min | 2:45 |
| Wrap-up                                                       | 10 min | 2:55 |

**Total:** ~3h 0min of content.
