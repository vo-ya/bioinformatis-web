# Lecture 24 (proposed L21) — CRISPR Functional Screens and DepMap

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology assumed
> **File**: provisional `lectures/lecture-24.html` — renumber to `lecture-21.html` when curriculum is reordered.

> **Proposed placement**: insert after L19 (causal-inference / GWAS, becomes new L20) and before L22 (data engineering, becomes new L22). The natural arc: GWAS gives candidate loci → MR establishes causal direction → CRISPR screens experimentally test gene-level dependencies → DepMap turns the result into therapy-guidance. Placing CRISPR before clinical genomics gives students the experimental-genomics machinery before the clinical-interpretation lecture.

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Describe the CRISPR-Cas9 mechanism in one paragraph; explain how it's repurposed for genome-scale functional screens.
2. Distinguish library types: CRISPR knockout (CRISPRko, Cas9), CRISPR interference (CRISPRi, dCas9-KRAB), CRISPR activation (CRISPRa, dCas9-VP64), base editors, prime editors.
3. Run a pooled-screen analysis end-to-end: from FASTQ → sgRNA counts → fitness scores via MAGeCK / DrugZ.
4. Interpret a screen as a positive-selection or dropout (negative-selection) experiment; compute log fold change of sgRNA abundance and aggregate to gene-level.
5. Use the DepMap (Cancer Cell Line Encyclopedia) for cancer-gene-dependency lookup; identify lineage-selective dependencies; prioritise drug targets.
6. Apply MAVEs (multiplexed assays of variant effects) and saturation mutagenesis screens for variant interpretation.
7. Frame pooled screens as compressed sensing on a perturbation matrix; sgRNA abundance as a noisy linear measurement; gene scoring as inverse problem.

---

## Part 1 — CRISPR for Functional Genomics (≈30 min)

### 1.1 The CRISPR mechanism (≈6 min)

**CRISPR-Cas9** (Doudna & Charpentier 2012, Zhang 2013): a programmable nuclease.

- **sgRNA (single-guide RNA)**: ~100 nt, 20 nt of which (the "spacer") is complementary to the target DNA.
- **Cas9**: a DNA endonuclease that binds the sgRNA and cleaves DNA where the sgRNA finds a match (with a downstream PAM sequence: NGG for SpCas9).
- **Result**: a double-strand break at the target site → cellular repair (NHEJ usually) → small indels → frameshift → loss-of-function.

For pooled screens, build a library of ~80,000 sgRNAs (4-6 per gene × ~20,000 genes); deliver via lentivirus; each cell gets ~1 sgRNA at a time; selection / dropout reveals which gene knockouts affect fitness.

### 1.2 The screen designs (≈10 min)

**Positive-selection screen**:

- Apply selection (drug, virus, hostile environment) that kills most cells.
- Surviving cells are enriched for sgRNAs targeting genes whose loss confers resistance.
- Read out: pre-selection vs post-selection sgRNA counts.
- High-fold-change sgRNAs → resistance genes.
- Examples: drug-resistance screens (BCR-ABL + imatinib → identify resistance pathway), viral entry (HIV co-receptor screens), tumour-suppressor identification.

**Negative-selection / dropout screen**:

- Grow cells over time without selection.
- sgRNAs targeting essential genes → cells with that knockout die / fail to proliferate → sgRNA counts drop.
- Read out: T0 (initial library) vs T-final (after 14-21 days).
- Drop-out sgRNAs → essential genes.
- Examples: cancer-essentiality (DepMap), context-dependent essentiality (essential in cancer cell lines but not normal — drug targets).

**FIGURE — Figure #1: Pooled screen workflow** → `diagrams/lecture-24/01-pooled-screen.svg`

### 1.3 CRISPRi and CRISPRa (≈6 min)

For non-cutting screens:

**CRISPRi (interference)**:

- **dCas9-KRAB**: catalytically-dead Cas9 fused to a KRAB transcriptional repressor.
- Targets transcription start sites of genes; suppresses transcription without permanent DNA edits.
- **Reversible** (vs CRISPRko which permanently edits).
- Better for essential genes (don't kill cells with single-cell-fatal knockouts).

**CRISPRa (activation)**:

- **dCas9-VP64** (or SunTag, SAM): transcriptional activator.
- Drives expression of normally-silent genes.
- Used for gain-of-function screens (overexpression confers resistance to drug?).

CRISPRko vs CRISPRi vs CRISPRa give complementary information: knockout (gene's removed), interference (gene's reduced), activation (gene's amplified). Modern screens often run all three.

### 1.4 Base editors and prime editors (≈4 min)

**Base editors** (Komor 2016): create specific point mutations without DSBs. Useful for:

- Saturation mutagenesis at a specific gene (every codon changed).
- MAVEs at single-base resolution.
- Discovery of disease-causing variants.

**Prime editors** (Anzalone 2019): even more flexible — write arbitrary edits. Slower per cell but precise.

### 1.5 Pooled vs arrayed screens (≈4 min)

- **Pooled**: all sgRNAs in one tube; phenotype-cell linkage by sgRNA-amplicon sequencing. ~$20k for genome-scale screen.
- **Arrayed**: each sgRNA in a separate well; per-well phenotyping (e.g., imaging, RNA-seq). 10–100× more expensive but allows complex phenotypes.

Most discovery screens are pooled; arrayed screens dominate for downstream characterisation.

---

## Part 2 — Library Design and Quality (≈25 min)

### 2.1 Library catalogues (≈5 min)

Standard human CRISPRko libraries:

- **Brunello** (Doench 2016): 76k sgRNAs, 4 per gene, optimised activity scores.
- **GeCKO v2** (Sanjana 2014): legacy ~120k library.
- **TKOv3** (Mair 2019): 70k sgRNAs, optimised for dropout screens, anti-correlation control.
- **CRISPick** library generator: design custom libraries.

For CRISPRi: **Dolcetto** (Dolcetto-Cas13) or **Horlbeck** library.

### 2.2 sgRNA scoring (≈8 min)

Not all sgRNAs work equally. Two scoring axes:

**On-target activity**: which sgRNAs cut efficiently?

- **Doench Rule Set 2** (Doench 2016): linear regression on sgRNA features → activity score.
- **Azimuth 2.0**: deep-learning extension.
- High-activity sgRNAs reduce noise in screens.

**Off-target effects**: which sgRNAs have unintended targets?

- **CFD score**: per-position mismatch tolerance.
- **MIT score**: alternative.
- High-specificity sgRNAs reduce false positives.

Library design balances activity (high) and specificity (high) against catalogue coverage (4-6 sgRNAs per gene).

### 2.3 Multiplicity of infection (≈5 min)

**MOI** = average number of viral particles per cell. For a clean screen, MOI ~0.3 → ~70% of cells get exactly one sgRNA, ~5% get two or more.

Higher MOI → many cells get multiple sgRNAs → confounded phenotypes.
Lower MOI → fewer infected cells → less data.

The MOI = 0.3 sweet spot is standard.

### 2.4 Coverage and replicates (≈4 min)

- **Cell coverage**: ~500 cells per sgRNA at infection (so each sgRNA has enough founders).
- **Read coverage**: ~500 reads per sgRNA at sequencing.
- **Biological replicates**: typically 2-3 independent infections.

Genome-scale screens require ~100M cells per replicate (75k sgRNAs × 500 cells × MOI 0.3 / 0.3).

### 2.5 The deep dive (≈3 min)

> **EE framing — pooled screens as compressed sensing**: Each cell gets a (typically) single sgRNA, but the population is a mixture of ~80k different perturbations at known fractional abundances. Phenotype readout is the pooled effect of all perturbations. Recovering per-gene fitness from pooled-cell phenotypes is **inverse-problem solving** under known mixing matrix — directly analogous to compressed sensing. The MAGeCK / DrugZ algorithms are statistical solvers for this; modern variants (MAUDE, BAGEL2) explicitly invoke compressed-sensing machinery.

**FIGURE — Figure #2: sgRNA library quality characteristics** → `diagrams/lecture-24/02-library-quality.svg`

---

## Part 3 — Screen Analysis with MAGeCK (≈40 min)

### 3.1 The data (≈4 min)

Output of a screen: **count tables**.

- Rows: sgRNAs.
- Columns: samples (typically pre-selection + post-selection × 2-3 replicates).
- Cells: sgRNA read counts.

For a typical screen, ~80k sgRNAs × 6-8 samples → ~600k cells in the matrix. Each sgRNA has 4-6 instances per gene, so the underlying gene-level matrix is ~20k × 6-8.

### 3.2 sgRNA-level fold change (≈6 min)

Per-sgRNA log fold change (LFC):

$$\text{LFC}_{i} = \log_2 \frac{\text{count post}_i + 1}{\text{count pre}_i + 1}$$

- Negative LFC: sgRNA dropped out → its target gene is essential.
- Positive LFC: sgRNA enriched → its target gene's loss provides selection advantage.

Normalisation: median normalisation, or DESeq-style size factors (count data, after all).

### 3.3 Gene-level aggregation (≈10 min)

Per gene: combine 4-6 sgRNA LFCs into one gene-level score + p-value.

**MAGeCK** (Li 2014, 2015): the standard.

- Computes per-sgRNA p-values (negative-binomial test on counts).
- Aggregates per-gene via **modified RRA (robust rank aggregation)**: median-rank-based test that's robust to a few non-functional sgRNAs.
- Output: per-gene log fold change, p-value, FDR.

**DrugZ** (Wang 2017, 2018): alternative for drug-vs-control screens.

- Z-score-based; specifically designed for resistance-vs-sensitisation.
- More sensitive to gene-level effects than MAGeCK in moderate signals.

**BAGEL2** (Kim 2019): Bayesian classifier for essentiality. Probabilistic essentiality prediction.

For a typical screen: run all three; cross-validate.

### 3.4 Quality control (≈8 min)

Critical QC metrics:

- **sgRNA recovery**: fraction of library sgRNAs detected with ≥ 5 reads. Target ≥ 90%.
- **Replicate concordance**: Pearson correlation between replicates at sgRNA level. Target ≥ 0.9.
- **Essential gene recovery**: fraction of CEGs (core essential genes) showing significant dropout. Target ≥ 80%.
- **Non-essential gene null distribution**: proper centring at LFC = 0 and adequate tail behaviour.

If QC fails, the analysis is uninterpretable. Re-run, don't bandage.

### 3.5 Worked example (≈8 min)

**FIGURE — Figure #3: MAGeCK analysis flow** → `diagrams/lecture-24/03-mageck.svg`

A classic dropout screen on cancer cell lines:

- Library: TKOv3, ~70k sgRNAs.
- 6 cell lines, 2 replicates each, T0 + T21 days.
- MAGeCK output: ~1000 essential genes per cell line, with cell-line-specific subsets ("contextual essentialities").

> **EE framing — gene-level p-value as inverse Z-test on aggregated rank statistic**: MAGeCK's RRA is mathematically the **modified Beta distribution** for the minimum order statistic of a uniform sample. Each gene has $k$ sgRNAs producing $k$ rank-percentile values; under the null (no effect), these are uniform on [0,1]. The MAGeCK score = best (smallest) percentile. Under the null, this is Beta(1, $k$)-distributed; significance is computed analytically. This is the same Mann-Whitney-style rank aggregation used in GSEA (L22-new) — non-parametric, robust to outliers, uniformly minimum-variance unbiased.

### 3.6 Modern variants (≈4 min)

- **MAUDE** (Doench 2018): integrates with the Brunello library design.
- **MELANIE** / **CRISPRcleanR**: corrections for copy-number-variation effects (high-CN regions can fail screens because Cas9 cuts every copy → huge fitness penalty unrelated to gene function).
- **CASA**: cell-line-aware model (handles polymorphic libraries).

---

## Part 4 — DepMap and Cancer Dependencies (≈30 min)

### 4.1 The DepMap project (≈4 min)

**Cancer Dependency Map** (Broad / Sanger / Sanger-DepMap collaboration): genome-scale CRISPR knockout screens in ~1000+ cancer cell lines spanning ~30 cancer types.

Public release (depmap.org): per-gene fitness score (CERES, then Chronos) per cell line, plus matched omics (RNA-seq, methylation, mutations, copy number).

This is the **Atlas of Druggable Cancer Dependencies** — for any candidate drug target, you can ask: "in which cancer cell lines is this gene essential?".

### 4.2 The Chronos algorithm (≈5 min)

**CERES** (Meyers 2017): the original DepMap scoring algorithm; corrects for copy-number bias.

**Chronos** (Dempster 2021): the current scoring algorithm.

- Models per-cell-line, per-sgRNA log fold change.
- Corrects for sgRNA-level activity, copy-number effect, library-batch effect.
- Output: per-gene "gene effect" score for each cell line.
- Score < -1: strong essential. Score = 0: non-essential.

The DepMap ~21,000 genes × ~1000 cell lines is the result.

### 4.3 Lineage-selective dependencies (≈8 min)

For drug-target prioritisation, the key signal is **lineage-selective essentiality**:

- Gene is essential in some cancers but not others.
- → Targeting it could kill specific cancer types without harming normal cells.

Major examples found via DepMap:

- **MCL1** in haematological cancers vs solid tumours.
- **MTAP-deletion synthetic lethality**: MTAP loss in 15% of cancers → PRMT5 dependency. Drug development target.
- **EZH2** in lymphomas with EZH2 mutations.
- **WRN helicase** in MSI-high cancers (synthetic lethal).

The 2024 frontier: ~50 lineage-selective dependencies identified; ~10 in clinical trials; ~3 in advanced clinical use.

### 4.4 Drug-screen integration (≈6 min)

DepMap also contains **drug-screen data** (PRISM, GDSC, CTRPv2): IC50 / AUC values for ~1000 drugs across ~500 cell lines.

Cross-referencing CRISPR + drug screens:

- Drug X kills cell lines where gene Y knockout also kills.
- → Drug X likely targets gene Y or its pathway.
- Used to deduce mechanism of action for novel drugs.

### 4.5 Synthetic lethality (≈4 min)

**Synthetic lethality**: gene A is non-essential, gene B is non-essential; but A + B knockout is lethal. For cancer therapy: if cancer has loss of A (driver mutation), targeting B becomes a cancer-selective killer.

Discovery pipeline:

- Identify cancer-specific gene losses (driver deletions).
- Run CRISPR screens in cancer cell lines vs normal.
- Look for cancer-selective dependencies → synthetic-lethal candidates.

PARP inhibitors in BRCA1/2-deficient cancers (L17, L18) are the canonical clinical example.

### 4.6 Worked example (≈3 min)

**FIGURE — Figure #4: DepMap dependency landscape** → `diagrams/lecture-24/04-depmap.svg`

**EMBED — Artifact #4: DepMap Dependency Browser** → `artifacts/lecture-24/04-depmap-browser.html`

---

## Part 5 — MAVEs and Variant Effect Maps (≈25 min)

### 5.1 The MAVE concept (≈4 min)

**MAVEs (Multiplexed Assays of Variant Effects)**: experimental measurements of all possible variants in a region. Most commonly: deep mutational scanning (DMS) on a protein.

- Generate library of ~all single-amino-acid changes (~20× protein length).
- Couple to a phenotypic assay (binding, growth, expression).
- Count library members before vs after selection.
- Score each variant.

Result: a **per-position-per-AA fitness map** for a protein region.

### 5.2 Saturation mutagenesis with base editors (≈6 min)

Base editors (Cas9-cytosine-deaminase): introduce point mutations at high efficiency without DSBs.

Workflow:

- Library of ~all possible base-edit sgRNAs targeting a gene.
- Each edit produces a specific point mutation.
- Pooled introduction; phenotype selection; sequencing of edited reads.
- Score each variant.

Example: **BRCA1 saturation mutagenesis** (Findlay 2018): tested ~4000 BRCA1 single-base variants for HDR function → published variant-effect map → directly informs ACMG/AMP variant classification (L17).

### 5.3 Variant effect maps in clinical genomics (≈5 min)

The bridge to L17 (clinical genomics):

- ACMG/AMP PS3 (functional studies show damaging effect): MAVEs provide systematic functional data.
- For genes with published MAVEs (BRCA1, MSH2, TP53, etc.), variant interpretation gets a quantitative basis.
- The **MAVEdb** repository aggregates published MAVEs across ~50 genes.

### 5.4 Computational complement: deep-learning MAVE (≈5 min)

For genes without experimental MAVEs:

- **AlphaMissense** (L17): deep-learning predictor trained partly on MAVE data.
- **EVE** (Frazer 2021): unsupervised deep-learning model that learns variant effects from sequence alone.
- Both feed PP3/BP4 evidence in clinical interpretation.

### 5.5 The 2024 frontier (≈3 min)

- **Combinatorial MAVEs**: paired mutations (epistasis maps).
- **Single-cell MAVEs**: link variant to phenotype at single-cell resolution.
- **AI-MAVE integration**: train deep models on existing MAVEs + apply to genes without experimental data.

### 5.6 Worked example (≈2 min)

**FIGURE — Figure #5: MAVE variant effect map** → `diagrams/lecture-24/05-mave-map.svg`

---

## Part 6 — Tools, Pitfalls, and Workflows (≈25 min)

### 6.1 Standard 2024 stack (≈5 min)

For a typical CRISPR-screen analysis:

- **MAGeCK** (Python/R): canonical tool.
- **DrugZ** (R): drug screens.
- **BAGEL2** (R): essentiality classifier.
- **CRISPRcleanR** (R): CN bias correction.
- **MELANIE / CASA**: advanced library-aware tools.

For visualisation: **R/Bioconductor**, **MAGeCKFlute**, or in-house Python/Plotly.

### 6.2 Pitfalls (≈8 min)

**Library**:

- Use a recent, well-validated library (Brunello, TKOv3, Dolcetto). Legacy libraries (older GeCKO) have many low-activity sgRNAs.
- 4-6 sgRNAs/gene minimum; fewer means high noise.

**Experimental**:

- MOI control: > 0.3 → confounded multi-perturbations; enforce strict MOI < 0.3.
- Cell coverage: < 500/sgRNA → noisy results.
- Replicate concordance: < 0.9 means major batch effect; investigate before interpreting.

**Bioinformatics**:

- Copy-number-variation bias: cancer cell lines often have high CN regions; sgRNAs in those regions cut more → inflated fitness penalty unrelated to gene function. CRISPRcleanR or CERES handles this.
- Library DNA contamination: library plasmid DNA contaminates the screen; hash check the count distribution.
- Outlier essential genes: a few highly-essential genes always dominate; trim or rank-aggregate.

**Statistical**:

- Multiple testing: ~20k genes tested; FDR control essential.
- Power: at typical signal sizes, n=2 replicates is barely enough; n=3 is standard; n=4 desirable for weak phenotypes.

### 6.3 Reproducibility (≈4 min)

CRISPR screens have notably better reproducibility than microbiome (~80% replication), driven by:

- Standardised libraries.
- Cell-line uniformity (vs human-cohort heterogeneity).
- Fewer biological covariates (lab-adapted cell lines vs human variation).

But cell-line-specific results don't always generalise to primary tissue — keep cell-line vs tissue distinction in mind.

### 6.4 The 2024 frontier (≈4 min)

- **In vivo CRISPR screens**: pooled screens in transplanted tumours / animal models.
- **Single-cell CRISPR screens**: Perturb-seq, ECCITE-seq → per-cell phenotype + perturbation linkage.
- **Combinatorial perturbations**: paired sgRNAs → epistasis screens.
- **CRISPR-driven RNA targeting**: Cas13a; targeting RNA rather than DNA.

### 6.5 Hands-on (≈4 min)

A reasonable first project:

- Download a published DepMap screen (raw count tables on the DepMap portal).
- Run MAGeCK / DrugZ.
- Identify essential genes.
- Compare to the published list.

Total: ~1 day, fits on a laptop.

**FIGURE — Figure #6: CRISPR screen analysis workflow** → `diagrams/lecture-24/06-workflow.svg`

---

## Part 7 — From Screen to Therapy (≈15 min)

### 7.1 The discovery → drug pipeline (≈4 min)

A typical drug-discovery program informed by CRISPR screens:

1. **Disease model**: cancer cell line / iPSC / patient-derived organoid.
2. **CRISPR screen**: identify essential genes / pathway dependencies.
3. **Hit triage**: filter for druggable targets (small-molecule-binding pockets, surface antigens).
4. **Validation**: arrayed CRISPR knockout, isogenic-pair experiments.
5. **Drug discovery**: virtual screening, hit-to-lead, medicinal chemistry.
6. **Clinical trials**: based on validated target.

The 2010s saw ~5 clinical-grade targets emerge from this pipeline; 2020s expects 10x more.

### 7.2 Pharma adoption (≈4 min)

CRISPR screens are now standard in:

- **Roche / Genentech / Pfizer / Novartis** internal target ID.
- **Recursion / Tempus / In8 / Nimbus**: AI-driven CRISPR-screen-based drug discovery.
- **Vertex Pharmaceuticals**: pioneered MAVE-driven gene-therapy product (Casgevy, FDA-approved 2023, sickle cell disease).

### 7.3 Cancer therapy (≈4 min)

In oncology, DepMap-driven targets in clinical development:

- **PRMT5 inhibitors**: MTAP-deleted cancers (TANGO, GSK).
- **Verismum WRN inhibitors**: MSI-high cancers.
- **PARP inhibitors**: BRCA-deficient cancers (L17, L18).
- **MCL1 inhibitors**: hematologic malignancies.

Each emerged from CRISPR-screen evidence + structural biology + clinical translation.

### 7.4 Beyond cancer (≈3 min)

- **Anti-virals**: CRISPR screens identify host factors required for viral replication. SARS-CoV-2 entry-factor screens (Daniloski 2021) found ACE2 + furin + cathepsin L.
- **Inflammation**: identify cytokine-pathway nodes for autoimmune drug development.
- **Neurodegeneration**: iPSC-derived neuron screens for Alzheimer's, Parkinson's targets.

**EMBED — Artifact #7: From Screen to Drug Target** → `artifacts/lecture-24/07-screen-to-target.html`

---

## Wrap-up (≈10 min)

### What you should take away

- **CRISPR-Cas9 + lentiviral library = genome-scale functional screen.** Knockout, interference, activation each give distinct information.
- **Pooled screens are compressed sensing** on a perturbation matrix; gene-level scoring is rank aggregation across replicate sgRNAs.
- **MAGeCK + DrugZ + BAGEL2** are the analytical workhorses. Critical QC: sgRNA recovery, replicate concordance, CEG dropout.
- **DepMap** is the cancer-dependency atlas: ~1000 cell lines × ~21k genes. Lineage-selective dependencies → drug targets.
- **MAVEs** (saturation mutagenesis) provide systematic variant-effect data → directly informs clinical variant interpretation (L17, L18).
- **CRISPR screens drive modern target discovery**: most emerging cancer drugs trace to CRISPR + structure + clinical translation.
- **EE framings**: pooled screen as compressed sensing; sgRNA abundance as noisy linear measurement; gene scoring as inverse problem; rank aggregation as modified Beta distribution test.

### Next lecture

Data engineering (existing L14, becomes new L22). The reproducibility infrastructure you'll learn there is what makes large-scale CRISPR + DepMap analysis possible.

### Homework

1. Download a public CRISPR-screen dataset (e.g., DepMap Public 23Q4). Run MAGeCK on raw counts; report top 10 essential genes.
2. Compare the same screen analysed by MAGeCK vs DrugZ vs BAGEL2. Where do they agree? Where differ?
3. Pick a cancer type from DepMap; identify the top 5 lineage-selective essential genes. Cross-reference with PRISM drug-screen data; identify drugs targeting these genes.
4. From the BRCA1 MAVE (Findlay 2018), tabulate variant effects for the top 100 ClinVar Pathogenic variants. Compute concordance with experimental MAVE scores. Discuss disagreements.
5. Design a CRISPRko library targeting a 50-gene custom set. Pick optimal sgRNAs using CRISPick / Doench Rule Set 2. Report on-target activity and off-target predicted scores.

### Recommended reading

- Doudna, J. A., & Charpentier, E. (2014). The new frontier of genome engineering with CRISPR-Cas9. *Science* 346, 1258096.
- Doench, J. G., et al. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nature Biotechnology* 34, 184–191.
- Sanjana, N. E., et al. (2014). Improved vectors and genome-wide libraries for CRISPR screening. *Nature Methods* 11, 783–784.
- Mair, B., et al. (2019). Essential gene profiles for human pluripotent stem cells identify uncharacterized genes. *Cell Reports* 27, 599–615.e4. (TKOv3.)
- Li, W., et al. (2014). MAGeCK enables robust identification of essential genes from CRISPR/Cas9 screens. *Genome Biology* 15, 554.
- Wang, B., et al. (2017). Robustly improved CRISPR screening data filtering. *Genome Biology* 18, 165. (DrugZ.)
- Meyers, R. M., et al. (2017). Computational correction of CN bias in CRISPR screens. *Nature Genetics* 49, 1779–1784. (CERES.)
- Dempster, J. M., et al. (2021). Chronos: a CRISPR cell fitness time-series model. *Genome Biology* 22, 343.
- Findlay, G. M., et al. (2018). Accurate classification of BRCA1 variants with saturation genome editing. *Nature* 562, 217–222.
- Tsherniak, A., et al. (2017). Defining a Cancer Dependency Map. *Cell* 170, 564–576.e16.
- DepMap: <https://depmap.org/>
- MAGeCK: <https://sourceforge.net/p/mageck/wiki/Home/>
- MAVEdb: <https://www.mavedb.org/>
- CRISPick: <https://portals.broadinstitute.org/gpp/public/analysis-tools/sgrna-design>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — CRISPR for Functional Genomics              | 30 min | 0:30 |
| Part 2 — Library Design and Quality                    | 25 min | 0:55 |
| Part 3 — Screen Analysis with MAGeCK                   | 40 min | 1:35 |
| Part 4 — DepMap and Cancer Dependencies                | 30 min | 2:05 |
| Part 5 — MAVEs and Variant Effect Maps                  | 25 min | 2:30 |
| Part 6 — Tools, Pitfalls, and Workflows                | 25 min | 2:55 |
| Part 7 — From Screen to Therapy                          | 15 min | 3:10 |
| Wrap-up                                                    | 10 min | 3:20 |

**Total:** ~3h 20min of content.
