# Lecture 17 — Clinical Genomics, Variant Interpretation, and Ethics

> **Duration**: ≈3h 30min content
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-17.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. Distinguish research and clinical bioinformatics pipelines; explain what CLIA / CAP certification, LDT vs IVD, and clinical validation require.
2. Apply the ACMG/AMP 2015 variant-classification framework: identify evidence codes (PVS1, PS1–4, PM1–6, PP1–5, BA1, BS1–4, BP1–7), combine them via the rule set, and arrive at a pathogenic / likely pathogenic / VUS / likely benign / benign call.
3. Assemble a variant evidence dossier using the standard resources — gnomAD for population frequency, REVEL / AlphaMissense for missense-impact prediction, SpliceAI for splice effects, ClinVar for prior classifications.
4. Describe the current state of pharmacogenomics: CYP2D6/codeine, VKORC1/CYP2C9/warfarin, HLA-B*57:01/abacavir, TPMT/thiopurines; use PharmGKB and CPIC guidelines.
5. Explain the ACMG secondary-findings list (v3.x), the right-not-to-know debate, and the operational workflow when incidental findings appear.
6. Outline the US and EU regulatory landscape: FDA-approved sequencing assays (FoundationOne, MSK-IMPACT), the 2023 FDA LDT rule, EU IVDR, and DTC genomics (23andMe history).
7. Articulate the ethics and data-sovereignty issues: GINA protections and their limits, ancestry bias in clinical databases, the Havasupai case, indigenous data sovereignty (CARE principles), and initiatives to improve diversity (All of Us, Our Future Health, H3Africa).

---

## Part 1 — From Research to Clinic (≈20 min)

### 1.1 Why this lecture (≈3 min)

Previous lectures treated bioinformatics as research: pipelines, methods, open questions. This lecture covers the transition from "research result" to "clinical deliverable." When your pipeline's output informs a patient's diagnosis, dosing, or surgical decision, the rules change. Every EE student entering the genomics industry will, sooner or later, write code that touches this boundary.

Clinical bioinformatics differs from research bioinformatics in five fundamental ways:

- **Accuracy at the individual level**. A research pipeline reporting 99.8% recall is excellent; a clinical pipeline that misses 0.2% of pathogenic variants has missed a patient.
- **Reproducibility and auditability**. Every call must be re-derivable years later; every step logged; every tool version pinned.
- **Turnaround time (TAT)**. Research: weeks. Clinical: 3–14 days for most tests.
- **Regulation**. Clinical labs are regulated (CLIA in the US, equivalents elsewhere); research labs are not.
- **Sign-off**. Every clinical report has a named physician signatory. Research papers have authors; clinical reports have a doctor whose license is at stake.

### 1.2 CLIA, CAP, and the lab ecosystem (≈6 min)

In the United States:

- **CLIA** (Clinical Laboratory Improvement Amendments, 1988). Federal regulation of all clinical labs. Any lab testing human samples for health-related results must be CLIA-certified. Genomic tests are "high-complexity" CLIA tests.
- **CAP accreditation** (College of American Pathologists). Voluntary-but-expected peer accreditation; more rigorous than CLIA minimums. Most respected clinical labs carry CAP.
- **NYS CLEP** (New York State Clinical Laboratory Evaluation Program). Separate state-level certification. NYS has the strictest clinical-lab regulation in the US. Getting "NY approval" is a meaningful bar.

Outside the US:

- **ISO 15189**: international medical-laboratory quality standard. Essentially the global equivalent of CLIA + CAP.
- **UKAS**, **DAkkS**: UK and German national accreditation bodies, respectively.

Bioinformatics in CLIA labs:

- Pipeline must have a written **validation document** (Vali-doc) demonstrating accuracy on known samples.
- Every release has a **standard operating procedure (SOP)** with version tracking.
- Any change — even a tool version bump — requires re-validation or at minimum documentation.
- **Internal Quality Control (IQC)**: known positive + negative controls run with each batch.
- **External Quality Assessment (EQA)**: periodic blind proficiency tests (e.g. GenQA).

### 1.3 LDTs, IVDs, RUO (≈6 min)

Assay regulatory classification:

- **LDT** (Laboratory-Developed Test). A test developed, validated, and performed within a single lab. Historically FDA's hands-off position; 2023–2024 rule changes (Part 6) now bring LDTs under FDA oversight.
- **IVD** (In-Vitro Diagnostic device). An FDA-cleared or FDA-approved assay sold as a kit (reagents, instruments, software). Much higher regulatory bar. Examples: FoundationOne CDx, Illumina TSO500.
- **RUO** (Research Use Only). Reagent or tool labelled for research, **not** clinical use. Many academic tools are RUO. Using RUO reagents for clinical reporting is a regulatory violation unless the lab performs their own clinical validation (turning the RUO input into a validated LDT output).

Most modern clinical genomic tests are LDTs built on a mix of FDA-cleared hardware (Illumina sequencers, IVDR-certified library-prep kits) + in-house analysis pipelines.

**FIGURE — Figure #11: Clinical lab regulatory landscape** → `diagrams/lecture-17/11-regulatory-landscape.svg`
*A layered diagram. Bottom: physical infrastructure (sequencers, reagents) with FDA-cleared and IVD-certified items labelled. Middle: the lab itself, with three overlapping certification rings: CLIA (required), CAP (expected), NYS CLEP (strictest). Top: the test output with three classification paths: RUO (non-clinical), LDT (lab-developed, now under FDA enforcement post-2024), IVD (FDA 510(k) or PMA approved). Arrows from physical infrastructure → lab → output, with regulatory constraints labelled at each level.*

**FIGURE — Figure #1: Research vs clinical pipeline comparison** → `diagrams/lecture-17/01-research-vs-clinical.svg`
*Two parallel workflow tracks. Top (research): sample → sequence → align → variant call → analysis → manuscript. Annotations: "1-off analyses", "version floats", "months TAT", "publication standard". Bottom (clinical): sample (with chain-of-custody) → sequence → aligned (against pinned reference) → variant call (validated pipeline) → curation (certified variant scientist) → signed clinical report → EHR. Annotations: "every step SOP-documented", "CLIA / CAP audited", "3–14 day TAT", "physician-signed".*

### 1.4 What changes in the code (≈5 min)

> **Intuition box**: The simplest way to feel the research-to-clinic shift is to picture rerunning last month's analysis. In research, you pull the latest tool version, re-align, see the results wiggle by a few percent, and write a note in the methods section. In a clinical pipeline, rerunning last month's patient sample must produce **bit-identical** output to the original report — or a documented re-validation event is required. Every knob is frozen; every file is hashed; every run has a UUID traceable back through the chain. This isn't bureaucratic overreach — it's what makes the downstream clinical decision legally defensible years later when a patient's lawyer (or the state medical board) asks "how did you reach this call?"

Concrete differences a clinical software engineer notices:

- **No bleeding-edge tool versions**. Production is typically 1–2 years behind the research frontier. Stability > novelty.
- **Containerised everything**. The container image used is part of the validation; updating it requires re-validation.
- **Explicit version pinning** of reference genome, databases (ClinVar, gnomAD, HGMD), every auxiliary file. Pin by MD5 hash.
- **Run-level audit trails**. Every run has a UUID; each file produced has a provenance chain back to raw reads.
- **Automated detection of validation drift**. Regular reruns on reference samples (e.g. HG002 / GIAB) to catch any regression.
- **Fail-closed, not fail-open**. Ambiguous outputs must be reviewed, not auto-reported. A caller that silently outputs "no calls" on a low-coverage region should throw a loud warning, not a quiet empty VCF.

> **EE framing — clinical pipeline validation as FDA / regulatory validation in EE**: Running a clinical bioinformatics pipeline is structurally the same as running an FDA-regulated medical device or RF-certified radio. Every change requires re-validation; every deployment has a version-locked release; every failure mode is characterised in a documented limitations section. The analogy in electrical engineering is the difference between a hobby radio circuit and an FCC-certified broadcast transmitter, or between a research-grade ECG and a 510(k)-cleared clinical one. The engineering is similar in kind — the regulatory and documentation burden is the differentiator. If you've worked in a safety-critical EE domain (automotive, aviation, medical devices), most of the culture here will feel familiar.

---

## Part 2 — ACMG/AMP Variant Classification (≈50 min)

### 2.1 Why a classification standard exists (≈5 min)

Before 2015, clinical labs classified variants independently; the same variant could be called "pathogenic" in one lab and "VUS" in another. This was a problem for patients transferring between providers, for meta-analyses, and for multidisciplinary teams.

In 2015 the **American College of Medical Genetics and Genomics (ACMG)** and the **Association for Molecular Pathology (AMP)** jointly published a 28-page document (Richards et al. 2015, *Genetics in Medicine*) defining the standardised framework. The framework:

- **Five-class output**: pathogenic (P), likely pathogenic (LP), variant of uncertain significance (VUS), likely benign (LB), benign (B).
- **Evidence codes**: individual lines of evidence (e.g. "extremely rare in gnomAD" or "predicted to disrupt protein structure"), each with a specific strength label.
- **Rules for combination**: how evidence codes combine to produce the five-class call.

The goal was to make variant classification **reproducible, auditable, and interpretable** across labs. It has largely succeeded; ACMG/AMP is now the operating standard in almost every clinical genomics lab worldwide.

> **Historical pointer**: The ACMG/AMP 2015 guidelines (Richards, Aziz, Bale et al.) emerged from working groups that had spent three years trying to reconcile the classification systems used at major US labs (GeneDx, Ambry, ARUP, Partners / Laboratory for Molecular Medicine, Myriad, Invitae, and others). The 2015 document codifies not a new classification scheme, but a **common language** across existing ones. The related **ClinGen** initiative (<https://clinicalgenome.org/>) maintains gene- and disease-specific refinements of ACMG/AMP — for some high-interest genes (BRCA1/2, MMR genes, cardiac channelopathies), ClinGen has published expert-curated modifications that supersede the generic rules. ClinVar (the NIH-hosted community classification database) absorbs these as they're published.

### 2.2 The five-class output (≈4 min)

Each class carries a probabilistic interpretation:

- **Pathogenic (P)**: ≥ 99% probability of disease-causing for the patient's phenotype.
- **Likely pathogenic (LP)**: 90–99% probability.
- **Variant of Uncertain Significance (VUS)**: 10–90% probability. **Reporting threshold**: VUS results are returned to physicians but without clinical action.
- **Likely benign (LB)**: 0.1–10% probability of pathogenic.
- **Benign (B)**: ≤ 0.1% probability.

Clinical actionability:

- **P / LP**: acted on (medication choice, surveillance, risk-reducing surgery, family testing).
- **VUS**: not acted on. Reported with a note. Reclassified as evidence accumulates.
- **LB / B**: not acted on. Usually not reported in a clinical context at all.

The **VUS problem**: the majority of variants in most genes are VUS. For rare genes with little population frequency data, VUS rates approach 50%. Reclassifying VUS → P/LP/B as evidence accumulates is a full-time job for variant scientists.

**FIGURE — Figure #2: ACMG/AMP five-class distribution** → `diagrams/lecture-17/02-acmg-classes.svg`
*A horizontal spectrum with five tiers labelled B → LB → VUS → LP → P. Under each, a probability range and an action label. Below the spectrum, a sketch of the distribution of clinical-exome variant calls: most variants are LB/B (common, known benign); a thick VUS tail; a narrow P/LP tip. Annotation: "VUS rate per gene: 15–50%".*

### 2.3 Evidence codes and their strengths (≈14 min)

Twenty-eight evidence codes, organised by direction (pathogenic / benign) and strength. Abbreviated nomenclature:

**Pathogenic evidence** (stronger → weaker):

- **PVS1** (Very Strong): **loss-of-function (LoF)** variant in a gene where LoF is a known disease mechanism. Examples: nonsense, frameshift, canonical splice site, single-exon deletion. Specific caveats apply — e.g. a stop-gain in the last exon may escape NMD and thus not be truly LoF.
- **PS1–PS4** (Strong): PS1 — same amino-acid change as a previously-established pathogenic variant; PS2 — de novo confirmed by parental testing + confirmed parentage; PS3 — well-established functional studies showing damaging effect; PS4 — statistically significant enrichment in cases vs controls.
- **PM1–PM6** (Moderate): PM1 — in a hotspot / critical domain; PM2 — absent from gnomAD (or extremely rare); PM3 — in trans with pathogenic variant in recessive disease; PM4 — protein-length change due to in-frame indel / stop-loss; PM5 — novel missense where different amino-acid change at same residue is known pathogenic; PM6 — assumed de novo without confirmed parentage.
- **PP1–PP5** (Supporting): PP1 — co-segregation with disease in family; PP2 — missense in a gene with low benign-missense rate; PP3 — multiple in-silico predictors agree on damaging effect; PP4 — patient phenotype highly specific for the disease; PP5 — reported as pathogenic in a reputable source (now **deprecated** as of 2018 ClinGen update).

**Benign evidence**:

- **BA1** (Stand-alone): allele frequency ≥ 5% in gnomAD. Stand-alone = enough by itself to classify as B.
- **BS1–BS4** (Strong): BS1 — allele frequency inconsistent with disease prevalence; BS2 — observed in healthy adults at an allele count inconsistent with disease; BS3 — well-established functional studies showing no damaging effect; BS4 — lack of segregation in affected family members.
- **BP1–BP7** (Supporting): BP1 — missense in a gene where truncating variants (not missense) cause disease; BP2 — observed in trans with a pathogenic variant in a recessive gene, in a phenotypically healthy individual; BP3 — in-frame indel in repetitive region without known function; BP4 — multiple predictors agree benign; BP5 — observed in a case with an alternative molecular diagnosis; BP6 — reported benign in a reputable source (**deprecated** 2018); BP7 — silent variant not at splice consensus with no predicted splice effect.

**FIGURE — Figure #3: ACMG/AMP evidence code table** → `diagrams/lecture-17/03-evidence-codes.svg`
*A structured table. Rows: Pathogenic (PVS1, PS1–4, PM1–6, PP1–5); Benign (BA1, BS1–4, BP1–7). Columns: Code, Description, Strength (Very Strong / Strong / Moderate / Supporting / Stand-alone). Colour-coded by strength: red deep → red light for pathogenic, cobalt deep → cobalt light for benign. BA1 highlighted separately in amber as stand-alone.*

### 2.4 Rule combination (≈12 min)

Once evidence codes are assigned, the **combinatorial rules** determine classification:

**Pathogenic (P)**:

1. 1× PVS1 + ≥ 1× (PS or PM or PP with specific weighting), OR
2. ≥ 2× PS, OR
3. 1× PS + ≥ 3× PM, OR
4. 1× PS + 2× PM + 2× PP, OR
5. 1× PS + ≥ 4× PP, OR
6. ≥ 3× PM + 2× PP (etc.)

**Likely pathogenic (LP)**:

- 1× PVS1 + 1× PM, OR
- 1× PS + 1–2× PM, OR
- 1× PS + ≥ 2× PP, OR
- ≥ 3× PM, OR
- etc.

**Benign (B)**:

- 1× BA1 alone, OR
- ≥ 2× BS

**Likely benign (LB)**:

- 1× BS + 1× BP, OR
- ≥ 2× BP

**VUS**: default when none of the above rules are met, or when pathogenic and benign evidence codes conflict.

The specific rules are a lookup table that's essentially memorised by variant scientists. Modern tools (InterVar, VarSome, Franklin) automate the lookup — but understanding what the rules *say* is still the core skill.

> **EE framing — ACMG/AMP rules as structured rule-based classification**: ACMG/AMP classification is an explicit, auditable **rule-based classifier**. Each evidence code is a Boolean (applies / doesn't apply); the rules combine them via a lookup table to produce the class label. This is deliberately **not** end-to-end ML. Reasons: (a) auditability — every call can be traced back to specific evidence; (b) regulatory acceptance — FDA / CE-IVDR approval processes require explainable outputs; (c) expert agreement — the ACMG/AMP framework encodes decades of clinical-genetics consensus that a neural network would have to rediscover. An EE analogy: it's like the difference between a hand-coded FSM for safety-critical signal processing and a learned RNN. The FSM is less powerful but vastly more auditable. In domains where mistakes harm people, the FSM wins. ML-based pathogenicity predictors exist (Part 3) and *inform* evidence codes like PP3/BP4; they don't replace the rule-based framework. This will likely stay true for the foreseeable future.

### 2.5 ClinGen refinements (≈8 min)

The base ACMG/AMP rules are gene-agnostic. For high-interest genes, the **ClinGen** expert panels have published gene-specific refinements:

- **BRCA1 / BRCA2 (breast/ovarian cancer)**: specific rules for when PVS1 applies based on exon location; frequency thresholds tuned for founder populations.
- **TP53 (Li-Fraumeni syndrome)**: modifications for hotspot definitions.
- **MMR genes (Lynch syndrome)**: specific thresholds + additional evidence codes.
- **RASopathy genes**: modifications for de novo handling.
- **Cardiac channelopathies** (SCN5A, KCNQ1): modifications reflecting the high VUS burden in these genes.

When a gene-specific ClinGen refinement exists, **use it**. When it doesn't, fall back to generic ACMG/AMP.

### 2.6 Sherloc and other frameworks (≈3 min)

**Sherloc** (Nykamp et al. 2017, Invitae) is a refined variant-classification framework that extends ACMG/AMP with additional granularity. Specifics:

- More evidence-code subcategories.
- Explicit quantitative thresholds (points per code, summed to a score).
- Designed for high-throughput reclassification at scale.

Most Invitae / Ambry clinical reports follow Sherloc internally but report results using ACMG/AMP five-class labels.

### 2.7 ClinVar as community truth (≈4 min)

**ClinVar** (NIH, freely accessible) is the community's aggregate variant-classification database:

- ~2.5M unique variants as of 2024.
- Each variant has one or more "submissions" from clinical labs, research groups, expert panels.
- Submissions are classified using ACMG/AMP (or Sherloc); ClinVar aggregates them.
- Aggregation result: single-star (one submission) to 4-star (multiple concordant submissions + expert panel + practice guideline).

How clinical variant scientists use ClinVar:

- **Starting point**: what has been said about this variant?
- **Confidence check**: is this variant widely agreed on or controversial?
- **Historical tracker**: was a variant previously LP but reclassified to VUS? (This happens regularly and is important.)

Not all ClinVar submissions are equal. Single-star, no-citation submissions can disagree with 4-star expert-panel classifications; the expert panel is usually right.

> **EE framing — variant classification as multi-evidence Bayesian aggregation**: Under the hood, the ACMG/AMP rule table is an approximation to a Bayesian posterior over the latent "is this variant truly pathogenic?" hypothesis. Each evidence code is a likelihood ratio: PVS1 ≈ LR +350 (strong evidence for pathogenicity), PS ≈ LR +19, PM ≈ LR +4.3, PP ≈ LR +2.1, and on the benign side BA1 is an absolute threshold while BS, BP are the mirror downweights. Combining independent evidence = multiplying likelihood ratios; the rule table encodes which combinations cross the 99% / 90% / 10% / 0.1% posterior-probability thresholds corresponding to the five classes. Tavtigian et al. (2018) published the explicit Bayesian reformulation; it reproduces the 2015 rule table almost exactly. The EE analogy is **Neyman–Pearson multi-sensor detection**: independent sensors (evidence codes) with known likelihood ratios combine into a posterior; thresholds map posteriors to decisions. The rule-based form exists because that's what clinicians can audit; the Bayesian form exists because that's what's mathematically true.

**FIGURE — Figure #4: ClinVar growth + classification distribution** → `diagrams/lecture-17/04-clinvar.svg`
*Top half: growth curve of ClinVar submissions 2013 (launch) → 2024 (~5M submissions, ~2.5M unique variants). Log-y axis. Bottom half: pie / donut chart of classification distribution across unique variants: ~15% P/LP, ~40% VUS, ~40% LB/B, ~5% conflicting. Annotation: "VUS is the plurality — reclassification is the field's ongoing work".*

**EMBED — Artifact #1: ACMG/AMP Classifier** → `artifacts/lecture-17/01-acmg-classifier.html`
*Select evidence codes applying to a variant (PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7). Classifier applies the combinatorial rules; returns the five-class call with a reasoning trace. Target aha: classification is deterministic given evidence codes — the hard part is assigning codes, not combining them.*

---

## Part 3 — The Variant Evidence Ecosystem (≈30 min)

### 3.1 Population frequency: gnomAD (≈6 min)

The single most-used resource in variant classification is **gnomAD** (Karczewski et al. 2020). Allele frequencies across ~800k exomes + ~150k whole genomes, stratified by 8 continental ancestry groups.

Use cases in ACMG/AMP:

- **BA1**: variant present at > 5% in gnomAD → stand-alone benign.
- **BS1**: variant more common than disease prevalence (gene-specific threshold) → strong benign.
- **PM2_Supporting**: absent or extremely rare in gnomAD (various thresholds) → moderate-to-supporting pathogenic.

Specific gene-level frequency thresholds are important:

- For fully-penetrant autosomal dominant diseases: PM2 threshold usually 1/200k alleles.
- For recessive: PM2 threshold usually 1/20k.
- For high-penetrance oncogenes (BRCA1): specific ClinGen-recommended thresholds.

**Ancestry stratification** matters: a variant that's rare in Europeans but common in Africans is still benign; don't use European-only frequencies for a patient of African ancestry.

**FIGURE — Figure #5: gnomAD allele frequency and PM2/BA1 thresholds** → `diagrams/lecture-17/05-gnomad-thresholds.svg`
*A horizontal frequency axis, log-scale 10⁻⁷ to 10⁰. Annotated thresholds: BA1 at 5×10⁻² (shaded amber region to the right), BS1 at disease-prevalence-derived threshold, PM2 at 1/200k for dominant / 1/20k for recessive (shaded cobalt regions). Two example variants plotted: one rare (PM2 applicable), one common (BA1 applicable).*

### 3.2 Functional impact: REVEL, AlphaMissense (≈8 min)

For missense variants, **multiple in-silico predictors** score the likelihood of functional impact. Agreement among predictors supports PP3 (supporting pathogenic) or BP4 (supporting benign).

Major predictors:

- **REVEL** (Ioannidis et al. 2016). Ensemble of 13 older predictors (SIFT, PolyPhen-2, MutationTaster, GERP, etc.). Produces a 0–1 score. Widely used; thresholds: ≥ 0.7 for PP3, ≤ 0.15 for BP4.
- **AlphaMissense** (Cheng et al. 2023). AlphaFold2-derived; predicts likely-pathogenic / likely-benign / ambiguous for every possible missense substitution. State of the art in 2024–2025. Thresholds similar to REVEL.
- **CADD** (Kircher et al. 2014). Ensemble score trained on evolutionary conservation + functional features. Produces a "phred-like" scaled score (typically 0–40). Still widely used but now mostly superseded.
- **SpliceAI** (Jaganathan et al. 2019). Deep-learning predictor of splice-altering variants. Essential for evaluating variants near splice sites but not in canonical sites. Thresholds: ≥ 0.5 for splice disruption evidence (PP3 or PM5 depending on context).

Key caveats:

- These are *decision supports*, not deciders. PP3 is supporting evidence; don't call pathogenic on PP3 alone.
- Predictor **ancestry bias**: predictors trained predominantly on European data may underperform on African variants. This is a real and active research concern.
- **Agreement required**. Using just one predictor is discouraged; ACMG/AMP asks for multiple to agree. ClinGen panels often specify exactly which.

**FIGURE — Figure #6: REVEL / AlphaMissense / SpliceAI comparison** → `diagrams/lecture-17/06-predictors.svg`
*Three side-by-side panels, one per predictor. Each shows: input type (missense / missense / splice), typical score range, recommended thresholds for PP3 / BP4, and a sketch of the predictor's architecture (REVEL: ensemble of simpler models; AlphaMissense: AlphaFold-derived; SpliceAI: deep CNN on sequence). Bottom note: "consensus among predictors is required for PP3 / BP4, not any single one".*

### 3.3 Splicing: SpliceAI deep dive (≈5 min)

Splice-altering variants are an under-appreciated pathogenicity class. A variant 50 bp inside an intron can create a new splice site, skipping an exon and producing a truncated protein.

**SpliceAI** is the workhorse:

- Takes a 10 kb window of genomic sequence centred on the variant.
- Outputs per-position "acceptor gain", "donor gain", "acceptor loss", "donor loss" probabilities.
- If any of these exceed a threshold (0.5 typical), the variant has a predicted splice effect.

Use in ACMG/AMP:

- SpliceAI ≥ 0.5: splice effect evidence, applied as PP3 (supporting) or in specific cases PVS1_strong (very strong) if the variant disrupts a canonical splice site.
- SpliceAI ≥ 0.8: strong evidence of splice disruption; PS3_moderate or PVS1 depending on context.

### 3.4 Literature mining and ClinGen (≈5 min)

**Literature** — case reports, functional studies, family studies — contributes evidence codes PS3 (functional), PS4 (case-control statistics), PP1 (family segregation), PP4 (phenotype).

Tools to find relevant literature:

- **LitVar** (NCBI). Links variants to PubMed articles automatically. Saves hours of hand-searching.
- **Mastermind** (Genomenon). Commercial; broader coverage; subscription-based.
- **ClinVar submissions** themselves frequently include literature citations.

A typical clinical variant-scientist workflow for one variant:

1. Query gnomAD → check frequency across ancestries.
2. Query ClinVar → check prior classifications.
3. Query functional predictors → REVEL / AM / SpliceAI scores.
4. Literature search via LitVar or Mastermind → case reports, functional studies.
5. Query phenotype databases (HGMD, OMIM) → is this gene associated with the patient's disease?
6. Assemble evidence codes.
7. Apply ACMG/AMP rules.
8. Write a **classification rationale** (plain-text justification).
9. Enter the classification into the lab's LIMS; generate the patient-facing report.

Total time per variant: 30 minutes to 4 hours depending on complexity. A clinical variant scientist handles 10–30 variants per day.

### 3.5 Variant interpretation is a profession (≈6 min)

"Variant scientist" (or "variant analyst" / "clinical molecular geneticist") is a distinct professional role. Typical training: M.Sc. or Ph.D. in genetics / molecular biology, followed by 2–3 years of certification-track supervised practice, culminating in ACMG / board certification.

This is not a role an AI replaces today. ML predictors inform individual evidence codes (PP3/BP4); a variant scientist synthesises everything into a clinically-defensible classification with a human-readable rationale. Regulatory frameworks (CLIA, FDA) explicitly require this human-in-the-loop workflow.

The ML-augmented future: variant scientists' productivity doubles as AI pre-fills evidence codes and literature synthesis; human judgment still signs off.

**EMBED — Artifact #2: Variant Evidence Dossier Compiler** → `artifacts/lecture-17/02-evidence-dossier.html`
*Paste a variant (HGVS or chr:pos:ref:alt). Artifact compiles mock query results from gnomAD, ClinVar, REVEL, AlphaMissense, SpliceAI. Suggests applicable ACMG/AMP evidence codes. Shows the rationale. Target aha: variant classification is a structured information-synthesis task — the tooling automates the lookup; the synthesis is what a variant scientist does.*

---

## Part 4 — Pharmacogenomics (≈30 min)

### 4.1 The basics (≈5 min)

**Pharmacogenomics (PGx)**: how genetic variants affect drug response. A person's CYP450 enzyme genotypes affect how fast they metabolise many drugs. Standard dosing assumes average metabolism; for poor or ultra-rapid metabolisers, standard dosing can be dangerous.

Scope:

- Currently ~20 drug-gene pairs have **Clinical Pharmacogenetics Implementation Consortium (CPIC)** Level A evidence — "use this test result to change dosing."
- Another ~40 pairs at CPIC Level B — actionable recommendations but less universally adopted.
- Roughly 50% of FDA-approved drugs have some PGx label, though fewer are routinely genotyped clinically.

The underlying science: cytochrome P450 enzymes (CYP2D6, CYP2C9, CYP2C19, CYP3A4) and other pharmacokinetic enzymes have common loss-of-function or gain-of-function variants. Different ethnic groups carry these variants at different frequencies.

### 4.2 The canonical cases (≈13 min)

**CYP2D6 and codeine**:

- Codeine is a prodrug. CYP2D6 metabolises it to **morphine**, which provides analgesia.
- Poor metabolisers (~7% of Europeans): get no pain relief from codeine.
- Ultra-rapid metabolisers (~3% of Europeans, up to 30% in some North African populations): metabolise codeine so fast that morphine levels spike → respiratory depression. A few paediatric deaths from codeine-containing cough syrups in ultra-rapid metabolisers led the FDA to restrict codeine use in children.
- PGx action: know the patient's CYP2D6 metaboliser status; choose an alternative opioid (hydromorphone, oxycodone for normal metabolism via CYP3A4) if PM or UM.

**Warfarin (VKORC1 + CYP2C9)**:

- Warfarin is the classical anticoagulant. Dosing is notoriously difficult — too little and the patient clots; too much and the patient bleeds.
- **VKORC1** -1639G>A variant: enzyme-sensitivity locus. Affects warfarin requirement.
- **CYP2C9** *2 and *3 alleles: metabolism locus. Slower metabolism → lower dose needed.
- Genotype-guided dosing algorithms (Gage et al. 2008, IWPC 2009) predict steady-state dose from genotype + age + weight + INR. Reduces time-to-stable-INR by ~15%. Not yet universally adopted; DOACs (rivaroxaban, apixaban) have displaced warfarin for most indications.

**HLA-B*57:01 and abacavir**:

- Abacavir is an HIV reverse-transcriptase inhibitor.
- Patients carrying HLA-B*57:01 have a ~50% rate of severe hypersensitivity reaction. Without testing, ~5% of patients experience life-threatening reactions.
- **Pre-prescription HLA-B*57:01 genotyping is standard of care** in every HIV clinic globally. Reduced hypersensitivity incidence by > 95%.

**TPMT / NUDT15 and thiopurines**:

- Thiopurines (azathioprine, mercaptopurine) are used for leukaemia, autoimmune disease, and inflammatory bowel disease.
- TPMT enzyme metabolises them. TPMT deficiency (homozygous *3A, *3B, *3C) causes severe myelosuppression at standard doses.
- NUDT15 variants, common in East Asian populations, have a similar effect.
- Pre-treatment TPMT (and NUDT15 for Asian ancestry) testing reduces severe toxicity events.

### 4.3 Star-allele nomenclature (≈6 min)

PGx uses **star-allele** nomenclature to describe haplotypes:

- **CYP2D6*1**: the reference (wild-type) allele.
- **CYP2D6*4**: a common loss-of-function haplotype (carries a splicing variant + SNPs).
- **CYP2D6*10**: a reduced-function haplotype common in East Asians.
- **CYP2D6*17**: a reduced-function haplotype common in Africans.

Each star allele is a defined **combination of variants** at a gene — a haplotype. The PGx community maintains these definitions at **PharmVar** (<https://pharmvar.org/>).

Metaboliser phenotypes derive from star-allele combinations:

- **CYP2D6 *1/*1**: normal metaboliser.
- ***1/*4**: intermediate metaboliser.
- ***4/*4**: poor metaboliser.
- ***1/*2 ×N** (gene duplication): ultra-rapid metaboliser.

> **Intuition box**: Star alleles are "labels for common haplotypes". Just as "pathogenic / benign" is a clinical label on a variant, "*1 / *4" is a clinical label on a haplotype that combines multiple variants. For each star allele, a lookup table maps it to an activity score (typical: 0, 0.25, 0.5, 1, 2). Summing the two alleles' activity scores produces a metaboliser phenotype ("poor" = 0, "intermediate" = 0.25–1, "normal" = 1–2, "rapid" = 2–3, "ultra-rapid" = 3+). The star-allele abstraction hides the underlying variants from the clinician — they just need "which stars does this patient carry?" and the rest is lookup.

**FIGURE — Figure #7: Pharmacogenomics worked examples** → `diagrams/lecture-17/07-pgx-examples.svg`
*Four-panel grid, one per canonical PGx example. Each panel: drug / gene / variants / metaboliser phenotype / clinical action. CYP2D6 × codeine (poor / normal / ultra-rapid); VKORC1 + CYP2C9 × warfarin (dosing nomogram); HLA-B*57:01 × abacavir (avoid / proceed); TPMT × thiopurines (dose reduce / standard).*

### 4.4 PharmGKB and CPIC (≈3 min)

**PharmGKB** (<https://www.pharmgkb.org/>): the community-maintained PGx knowledge base. For each drug-gene pair:

- Known variants and their effects.
- Metaboliser phenotype definitions.
- Levels of evidence for clinical relevance.
- Links to dosing guidelines.

**CPIC** (Clinical Pharmacogenetics Implementation Consortium): drafts and publishes **dosing guidelines** per drug-gene pair. Format: "if patient genotype is X, give dose Y" tables.

Typical clinical workflow:

1. Order test → genotype the PGx panel (Illumina targeted panel, microarray, or WGS).
2. Software translates genotype → star alleles → metaboliser phenotype.
3. Clinician sees a report with dosing recommendations pulled from CPIC.
4. Clinician adjusts prescription.

Several health systems (Vanderbilt PREDICT, St. Jude PG4KDS, Mayo Clinic) pre-emptively genotype patients and have the PGx results available for any future prescribing decision. The model is working; scaling is slow.

> **EE framing — PGx as discrete-state control policy lookup**: Pharmacogenomics is a canonical **rule-based control policy**. Measure the patient state (genotype → star alleles → metaboliser phenotype), look up the recommended action (standard dose / reduced dose / avoid / alternative drug) in a fixed table (the CPIC guideline), apply. This is the simplest possible closed-loop control: one measurement, discrete state space, table lookup. No ML, no learning — because the evidence-generation machinery (clinical trials, mechanistic studies) has already done the optimisation; CPIC is the compiled answer. The EE analogue is a mode-based controller where each mode has pre-optimised parameters; the runtime's only job is to identify which mode applies. It's also why PGx is "easy" to productionise once the tables exist and "impossible" to productionise for drugs where the tables don't yet exist — the whole game is in building the lookup table.

### 4.5 PGx is under-used (≈3 min)

Despite strong evidence, PGx is **under-utilised** in routine care:

- Reimbursement (insurance) is patchy; patients often pay out of pocket.
- Electronic Health Record (EHR) integration is unevenly implemented; many physicians don't see the genotype when prescribing.
- Physicians are unfamiliar with metaboliser phenotypes beyond a handful of well-known cases.
- The list of actionable PGx drug-gene pairs is growing, but each addition requires reimbursement and EHR plumbing.

**EMBED — Artifact #3: PharmGKB Star-Allele Translator** → `artifacts/lecture-17/03-pharmgkb-star-allele.html`
*Pick a gene (CYP2D6 / CYP2C9 / TPMT / HLA-B) and a genotype (per-allele star-allele choices). Artifact translates to metaboliser phenotype and the corresponding CPIC dosing recommendation for a relevant drug. Target aha: PGx is "genotype → star alleles → phenotype → dosing" lookup, very mechanical once the tables exist.*

---

## Part 5 — Incidental Findings (≈20 min)

### 5.1 The problem (≈5 min)

A patient comes in for **cardiomyopathy testing**. The lab sequences their exome. In addition to the indication-relevant variant, the lab finds a pathogenic **BRCA1** mutation — unrelated to the test indication, but with serious implications (high lifetime breast / ovarian cancer risk).

**Report it? Don't report it? Depends on the patient's consent model?**

This is the **incidental finding** (or "secondary finding") problem. It arises naturally whenever broad sequencing (exome, whole genome) is performed.

### 5.2 The ACMG SF list (≈6 min)

In 2013, ACMG published a **recommended minimum list** of genes for which pathogenic variants should be reported as incidental findings, if the patient consents (or under an opt-out model, unless they opt out). This has become the **ACMG Secondary Findings (SF) list**.

- **SF v1** (2013): 56 genes, primarily cardiovascular + cancer predisposition.
- **SF v2** (2017): 59 genes.
- **SF v3.0 / 3.1 / 3.2** (2021–2024): ~81 genes currently.

Categories covered:

- **Hereditary cancer syndromes**: BRCA1/2, MLH1/MSH2/MSH6/PMS2/EPCAM (Lynch), TP53, APC, MUTYH, NF1/2, PTEN, RET, VHL, etc.
- **Familial hypercholesterolemia**: LDLR, APOB, PCSK9.
- **Inherited arrhythmia syndromes**: SCN5A, KCNQ1, KCNH2, RYR2.
- **Cardiomyopathies**: MYH7, MYBPC3, TNNT2, etc.
- **Malignant hyperthermia susceptibility**: RYR1, CACNA1S.
- **Several others**: Marfan, Ehlers-Danlos vascular type, Wilson disease, etc.

Key criteria for inclusion on the SF list: **actionable** — there's a clear medical intervention (surveillance, prophylactic surgery, lifestyle change, drug avoidance) that changes outcome if the variant is known.

**FIGURE — Figure #8: ACMG SF list v3.x categories** → `diagrams/lecture-17/08-acmg-sf.svg`
*A pie chart / treemap of SF v3.x genes by category. Hereditary cancer syndromes: largest slice (~30 genes). Cardiomyopathies + channelopathies: second largest (~25 genes). Familial hypercholesterolemia: ~3. Others: ~20. Accompanying table lists top 10 genes + their associated diseases.*

### 5.3 The right not to know (≈4 min)

Patients have a right to **not know** incidental findings if they prefer. The ACMG's evolving position:

- **2013**: opt-out model; labs must report unless patient declines.
- **2014**: opt-in model; labs must ask.
- **2021+**: explicit patient-facing consent; specific acknowledgment of SF-list findings.

The **arguments for reporting**:

- Actionable findings save lives.
- Patients benefit from information about their own genome.
- The incremental cost is near-zero once sequencing is done.

The **arguments against universal reporting**:

- Psychological distress in the absence of clear intervention.
- Downstream costs (surveillance, surgery) may not be reimbursed.
- Family members may be implicated without consent.
- Not all "actionable" findings clearly benefit patients; the evidence is genuinely mixed for some genes.

> **Discussion prompt**: A 65-year-old patient consents to exome sequencing for an adult-onset movement disorder indication. The exome reveals a pathogenic BRCA1 variant (hereditary breast / ovarian cancer risk). Should the lab return the finding? What if the patient is 85? What if the patient is a 15-year-old with a rare seizure disorder? What if the patient explicitly opted out of secondary findings at consent? (Starting points: patient age, clinical actionability at that age, consent scope, family implications, state / jurisdictional law.)

### 5.4 Incidental findings in practice (≈5 min)

Practical workflow when an incidental finding is suspected:

1. Confirm the variant (re-sequence with a different chemistry, ideally).
2. Re-classify the variant carefully using ACMG/AMP + ClinGen gene-specific guidelines.
3. Only report if:
   - Pathogenic or LP (not VUS).
   - Patient consented to receive SF-list findings.
   - Gene is on the current ACMG SF list (or otherwise clearly actionable).
4. Route the finding through the appropriate clinical workflow (usually a genetic counsellor before a physician).
5. Document in the clinical report under a distinct "Secondary Findings" section.

Notable failure modes:

- A VUS in a cancer-predisposition gene returned as actionable → patient has unnecessary prophylactic surgery.
- An SF-list variant missed because it wasn't queried → patient dies of a preventable cancer.
- An SF-list variant reported without consent → consent violation with legal liability.

**EMBED — Artifact #4: ACMG SF Incidental-Findings Checker** → `artifacts/lecture-17/04-sf-checker.html`
*Given a candidate variant (preset list), check: is the gene on the SF list? Does the variant meet the reporting threshold (P or LP only)? What was the patient's consent status? Produce a reporting recommendation. Target aha: reporting an incidental finding is a structured decision with multiple gates — gene-on-list + variant-pathogenic + patient-consented — all of which must be true.*

---

## Part 6 — Regulatory Landscape (≈25 min)

### 6.1 FDA-cleared sequencing assays (≈6 min)

The FDA clears sequencing assays through two main pathways:

- **510(k)** clearance: "substantially equivalent" to a predicate device. Faster, less rigorous.
- **PMA** (Premarket Approval): rigorous; required for high-risk Class III devices. Includes companion-diagnostic (CDx) approvals.

Major FDA-cleared / approved genomic assays:

- **FoundationOne CDx** (Foundation Medicine, Roche; 2017 FDA PMA). 324-gene tumour-profiling NGS panel. Approved as CDx for multiple targeted therapies. Heavy industry use.
- **MSK-IMPACT** (Memorial Sloan Kettering; 2017 FDA authorization). 468-gene tumour panel. Unusual: FDA "authorization" for an LDT, paving the way for the 2023 rule change.
- **Illumina TruSight Oncology 500**: large cancer panel.
- **Oncomine Dx Target Test**: ThermoFisher's NGS-based CDx for lung cancer.

For germline (inherited) testing, FDA-cleared kits are fewer: much of the market is LDTs (see next).

### 6.2 The 2023 FDA LDT rule (≈6 min)

Historically, the FDA claimed regulatory authority over LDTs but exercised **enforcement discretion** — in effect, not regulating them. This is how most clinical genomic testing operated.

In **April 2024**, the FDA finalised a rule ending enforcement discretion. Phased implementation over 4 years (2024–2028):

- **Year 1**: adverse-event reporting requirements.
- **Year 2**: quality system regulation compliance.
- **Year 3**: registration + listing requirements.
- **Year 4**: premarket review for high-risk LDTs.

Industry response has been mixed — some large labs (Mayo Clinic, Quest, LabCorp) are preparing; smaller labs are concerned about the regulatory burden. A Congressional review of the rule is ongoing as of 2025.

The practical effect for bioinformatics software: clinical pipelines will need more explicit regulatory documentation. Existing validation practices (Part 1) will need to formalise further.

### 6.3 EU IVDR (≈5 min)

The **In-Vitro Diagnostic Medical Devices Regulation (IVDR)** replaced the previous IVD Directive in May 2022 (after multiple delays).

Classifies IVDs into A/B/C/D risk classes. Genomic tests for severe hereditary disease are typically Class C; tumour profiling may be Class C or D.

Requirements:

- **Notified Body** certification for Class C/D devices.
- **Clinical performance studies** demonstrating accuracy on the intended population.
- **Post-market surveillance** with active follow-up.
- **CE mark** for market access.

Transition challenges: many labs that operated under the previous IVD Directive have struggled to complete IVDR certification in time. The EU has repeatedly extended deadlines. As of 2025, IVDR certification is a significant pain point for European clinical labs.

### 6.4 Direct-to-consumer testing (≈4 min)

DTC genomics is regulated differently from clinical testing.

- **23andMe** (founded 2007). Started returning health-related reports without FDA involvement. In **November 2013**, FDA issued a warning letter demanding cessation. 23andMe complied; shifted to ancestry-only reports.
- **April 2017**: FDA authorised 23andMe to return **limited health-related information** (e.g. BRCA variants for the three Ashkenazi-founder mutations — not general BRCA, a tight scope). The ancestry-plus-light-health model resumed.
- **March 2018**: 23andMe authorised for BRCA (three specific variants) DTC reports.
- **2022+**: 23andMe has broadened its clinical partnerships and diagnostic offerings but remains DTC-first.

**Regulatory question**: should a consumer directly learn they carry a BRCA variant without a physician interpreter? The 23andMe compromise (limited variants, clear disclaimers, optional genetic counselling) is the current answer, but remains controversial.

**FIGURE — Figure #9: FDA approval pathway** → `diagrams/lecture-17/09-fda-pathway.svg`
*A flowchart: starting point "assay concept" → fork into "Research Use Only" (right), "LDT" (center, formerly minimal oversight, now FDA-regulated post-2024), "FDA 510(k) clearance" (left), "FDA PMA approval" (far left for CDx). Each path annotated with typical use cases and examples. Arrow showing the 2024+ LDT rule's convergence toward the regulated pathway.*

### 6.5 What this means for bioinformatics software (≈4 min)

As regulation tightens, bioinformatics software is increasingly treated as part of the medical device:

- **Software as a Medical Device (SaMD)** FDA framework applies to clinical pipelines.
- **Change-management plans** required: every software update has a predetermined scope and validation requirements.
- **Design controls**: documented requirements, design decisions, testing, traceability.
- **Cybersecurity requirements**: especially if the software is cloud-deployed.

In practice: a bioinformatician writing clinical pipelines is writing regulated software. The skills overlap with medical-device firmware engineering (ISO 62304, IEC 62366-1) more than with research scientific computing.

> **Warning box**: In-silico missense predictors (REVEL, AlphaMissense, etc.) are **decision supports**, not **deciders**. They produce probability-like scores with thresholds, which inform a single supporting-evidence code (PP3 or BP4). A common misuse pattern — seen in both automated variant-classification tools and hastily-reviewed pipelines — is to treat the predictor score as directly actionable: "AlphaMissense says pathogenic, report as LP." This is wrong, regulated software or not. A single supporting-level code cannot produce LP by itself under ACMG/AMP; and predictors are known to over-call pathogenic for rare benign variants. Always combine predictors with population frequency, functional evidence, and ClinVar lookups; never sign a report on a predictor score alone.

---

## Part 7 — Ethics, Equity, and Data Sovereignty (≈25 min)

### 7.1 GINA and its limits (≈8 min)

The **Genetic Information Nondiscrimination Act (GINA)**, US federal law signed 2008. Covers:

- **Employment**: employers cannot use genetic information in hiring, firing, or promotion decisions.
- **Health insurance**: group and individual health insurers cannot use genetic information for coverage or pricing decisions.

Does not cover:

- **Life insurance**. A life insurer can ask about (and use) genetic test results.
- **Long-term care insurance**. Same.
- **Disability insurance**. Same.
- **Small employers** (<15 employees).
- **US military** (separate regulatory regime).

The carve-outs are a significant issue. A BRCA1 carrier testing positive may be unable to obtain affordable life insurance. Some states (California, Florida) have stricter state-level protections against insurance discrimination that fill the gaps.

> **Historical pointer**: GINA was signed into law by President George W. Bush on May 21, 2008. It passed the Senate 95–0 and the House 414–1 — extremely rare bipartisan support. Its long legislative history (13 years from first introduction to passage) reflected active industry resistance (insurance, employer groups) that eventually yielded to patient-advocacy coalition pressure after the 2003 Human Genome Project completion made the discrimination potential concrete. The carve-outs for life, disability, and long-term-care insurance were left in as deliberate compromises; bills to close them have been introduced periodically but not passed.

**FIGURE — Figure #10: GINA coverage map** → `diagrams/lecture-17/10-gina-coverage.svg`
*A Venn-diagram-like map of insurance / employment domains. Covered by GINA: health insurance (group + individual), employment (large employers). Not covered: life insurance, disability insurance, long-term care, small employers, military. Shaded regions indicate coverage status. State-law patchwork (California, Florida, Vermont, etc.) shown as small inserts.*

### 7.2 Ancestry bias in clinical databases (≈6 min)

**The problem**: clinical-genomics resources are dominated by people of European ancestry. Consequences cascade:

- **ClinVar**: ~70% of submissions are about variants seen in European patients.
- **gnomAD**: 56% European; African-ancestry subset is growing but remains smaller.
- **Functional predictors**: REVEL / AlphaMissense trained on datasets where European variants dominate.
- **PRS training cohorts**: overwhelmingly European (L13 §5.4).

Consequence for non-European patients:

- **Higher VUS rates**. Variants common in African populations are more often "rare and unknown" in European-dominated ClinVar — they get flagged as potential pathogens.
- **Lower predictor accuracy**. Missense predictors may underperform on populations under-represented in training.
- **Worse clinical actionability**. Fewer firm classifications → fewer actionable findings → worse care.

Initiatives to close the gap:

- **All of Us Research Program** (US): target 1M+ participants with strong minority representation.
- **Our Future Health** (UK): similar goal.
- **H3Africa**: African-focused genomic research consortium.
- **Biobank Japan**, **Biobank Korea**: East Asian coverage.

> **EE framing — ancestry bias as training-distribution coverage problem**: Ancestry bias in clinical databases is the training-distribution problem from L16 in its starkest form. A classifier (ACMG/AMP rules + predictors + ClinVar lookups) trained predominantly on European data doesn't generalise to African patients. The effect is analogous to deploying a model trained on one operating condition to a new condition where distributional shift is known and measurable — but no compensatory correction is applied. Closing the gap requires **diverse training data**, not just "fixing the algorithm"; the algorithm is working as designed on the data it has. This is why initiatives like All of Us, H3Africa, Our Future Health are essential infrastructure — they change the input distribution that every downstream tool builds on.

### 7.3 Havasupai and data sovereignty (≈4 min)

**The Havasupai Tribe v. Arizona Board of Regents** (case filed 2004, settled 2010):

- In 1989, Arizona State University researchers collected blood samples from the Havasupai Tribe (Grand Canyon region) for a diabetes-genetics study.
- Over the following decade, the samples were used for unauthorised research on schizophrenia, inbreeding, and population migration — research the Tribe had not consented to and considered culturally harmful.
- In 2004, Tribal members learned about the additional research.
- The Tribe sued; settled 2010; ASU returned the blood samples (a significant symbolic act); paid $700k; built a health clinic.

Consequences:

- The case became a defining precedent for **indigenous data sovereignty** in US research.
- NIH strengthened tribal-consultation requirements.
- Several large cohort studies now include explicit tribal governance provisions.

### 7.4 CARE principles and indigenous data sovereignty (≈3 min)

The **CARE Principles for Indigenous Data Governance** (Global Indigenous Data Alliance, 2019):

- **Collective benefit**: data ecosystems facilitate positive outcomes for Indigenous Peoples.
- **Authority to control**: Indigenous Peoples' rights to determine how their data are used.
- **Responsibility**: those working with Indigenous data have a duty to share benefits and acknowledge the relationship.
- **Ethics**: Indigenous Peoples' rights and wellbeing are the primary concern.

Complementary to the **FAIR Principles** (Findable, Accessible, Interoperable, Reusable) from the open-data movement; CARE adds the governance layer.

In practice: modern genomic studies involving Indigenous communities increasingly incorporate both FAIR + CARE, with community-level consent, data access restrictions, and benefit-sharing agreements.

### 7.5 Diverse cohorts and the path forward (≈4 min)

Several large-scale initiatives are working to change the landscape:

- **All of Us** (US, 2018+). Target 1M+ US participants with explicit mandate for 50%+ minority enrolment. Enrolment as of 2024: ~700k; data releases ongoing.
- **Our Future Health** (UK, 2020+). Target 5M UK participants.
- **H3Africa** (African Genomic Diversity Consortium, 2012+). Pan-African cohort with local infrastructure building.
- **TOPMed** (NHLBI, 2014+). Multi-ethnic WGS cohort (~150k genomes).
- **PAGE** (Population Architecture using Genomics and Epidemiology, 2008+). US minority-focused consortium.
- **Million Veteran Program** (VA, 2011+). US veterans, more diverse than UK Biobank equivalents.

The data these cohorts generate feeds back into ClinVar, gnomAD, AlphaMissense training sets, and PRS portability models. Progress is measurable: gnomAD v4 (2024) is materially more diverse than v2 (2020). Full equity remains years away.

> **Warning box**: Any clinical pipeline you deploy today has systematic accuracy differences by ancestry. This is not a hypothetical — it's measurable and published. When your pipeline signs out a VUS for a patient of African ancestry, the VUS rate is higher (in absolute terms) than it would be for a matched European patient. When you report ancestry bias as "known limitation" in your test's standard documentation, you are being honest. When you fail to acknowledge it, you are failing the patient.

**EMBED — Artifact #5: GINA Coverage Explorer** → `artifacts/lecture-17/05-gina-coverage.html`
*Scenario-picker: pick an insurance/employment/discrimination scenario (e.g. "BRCA+ patient applies for life insurance"). Artifact returns: GINA covered / not covered / state-law variations. Target aha: GINA's reach is narrower than people assume; life / disability / long-term care are the biggest gaps.*

**EMBED — Artifact #6: FDA / Regulatory Pathway Picker** → `artifacts/lecture-17/06-fda-pathway.html`
*Describe an assay (test type, clinical intent, sample-size for validation). Artifact maps it to the likely regulatory pathway (RUO / LDT / FDA 510(k) / FDA PMA / IVDR), explains what that pathway requires. Target aha: regulatory classification follows the assay's intended use, not a lab's preference; picking the wrong pathway can cost years of delay.*

**FIGURE — Figure #12: Ancestry bias in ClinVar and predictors** → `diagrams/lecture-17/12-ancestry-bias.svg`
*Three stacked bar charts. Top: ClinVar submissions by ancestry of source cohort (~70% European, ~10% East Asian, ~6% African, ~3% South Asian, ~2% Admixed American, ~9% unknown / mixed). Middle: gnomAD v4 ancestry composition, slightly more diverse (~56% European, ~16% African, ~15% East Asian, etc.). Bottom: REVEL / AlphaMissense validation F1 by ancestry, showing measurable gap (European ~0.88 vs African ~0.78). Annotation: "accuracy follows training-distribution coverage — the fix is diverse data, not algorithm tweaks".*

**EMBED — Artifact #7: Ancestry Bias in Variant Predictors** → `artifacts/lecture-17/07-ancestry-bias.html`
*Simulate / visualise REVEL and AlphaMissense performance on variants stratified by ancestry (European vs African). Show accuracy gap. Allow student to toggle "retrained on diverse corpus" and see the gap narrow. Target aha: ancestry bias is a training-distribution coverage problem; it's fixable by expanding the training distribution — gnomAD-scale effort, not an algorithm tweak.*

> **Intuition box**: The ACMG/AMP five-class framework is **not** two-valued "pathogenic vs benign" with a VUS tolerance band. It's a genuine five-way decision, because "we don't know" has a different clinical meaning than "probably not" or "probably yes". Clinicians act on P/LP; wait-and-watch on VUS; effectively ignore LB/B. The VUS bucket carries real information ("the lab looked and couldn't decide") that's lost if you collapse it to either extreme. When reading a clinical report, the distinction between LP (act) and VUS (don't act yet) drives the entire downstream care plan.

---

## Wrap-up (≈10 min)

### What you should take away

- **Clinical bioinformatics is regulated software**. CLIA/CAP/IVDR/FDA compliance shapes every design decision; validation and version-pinning are structural, not optional.
- **ACMG/AMP is the classification standard**. Evidence codes + combinatorial rules → P/LP/VUS/LB/B. Auditable, explainable, explicitly not end-to-end ML.
- **The evidence ecosystem is structured**. gnomAD (frequency), REVEL/AlphaMissense (functional impact), SpliceAI (splicing), ClinVar (prior classifications), literature (PS3/PS4/PP1/PP4). A variant scientist synthesises across these.
- **Pharmacogenomics is operational**. CYP2D6/codeine, VKORC1+CYP2C9/warfarin, HLA-B*57:01/abacavir, TPMT/thiopurines. CPIC guidelines drive dosing. PharmGKB / PharmVar are the reference tables.
- **Incidental findings require explicit consent + a gene list**. ACMG SF v3.x (~81 genes, actionable). Right-not-to-know is real; opt-in consent is current practice.
- **Regulation is tightening**. 2024 FDA LDT rule; EU IVDR; FDA's Software as a Medical Device framework. Clinical bioinformatics is increasingly regulated software engineering.
- **GINA protects narrowly**. Health insurance + employment yes; life / disability / long-term care no. The gaps matter.
- **Ancestry bias is systematic**. ClinVar / gnomAD / predictors dominated by European data; non-European patients get higher VUS rates, lower predictor accuracy, worse care. Cohort diversity (All of Us, H3Africa, Biobank Japan, etc.) is the structural fix.
- **Indigenous data sovereignty matters**. Havasupai, CARE principles. Research consent is community-level, not just individual-level, in these contexts.

### Next lecture

Cancer genomics: integrated capstone. Somatic vs germline; tumour heterogeneity; clonal evolution; MSK-IMPACT / FoundationOne; liquid biopsy (ctDNA); actionable mutations and targeted therapies; immunogenomics; cancer as a synthesis across L1–L17.

### Homework

1. Read the ACMG/AMP 2015 paper (Richards et al. 2015). Pick five evidence codes and give a concrete example for each — a published variant where that code applies.
2. For a synthetic variant you construct (missense, e.g., BRCA1 p.Arg1443Ter), run a classification workflow: pull gnomAD frequency, query ClinVar, query REVEL / AlphaMissense, apply ACMG/AMP rules. Document your reasoning at each step.
3. For one PGx drug-gene pair of your choice (not from the lecture's four canonical cases), read the CPIC guideline and summarise: which variants / star alleles drive which metaboliser phenotypes; what dosing recommendation follows.
4. Read the All of Us Research Program overview. Summarise its target cohort composition vs UK Biobank; what gaps it's designed to close; what the first major data releases showed.
5. Given the tightening LDT regulatory environment (2024 FDA rule; EU IVDR), write a one-page recommendation for how a small academic clinical lab should prepare for the next 3 years. What operational changes do they need to make?

### Recommended reading

- Richards, S., Aziz, N., Bale, S., et al. (2015). Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genetics in Medicine* 17, 405–424.
- Nykamp, K., Anderson, M., Powers, M., et al. (2017). Sherloc: a comprehensive refinement of the ACMG-AMP variant classification criteria. *Genetics in Medicine* 19, 1105–1117.
- Miller, D. T., Lee, K., Chung, W. K., et al. (2022). ACMG SF v3.1 list for reporting of secondary findings in clinical exome and genome sequencing: a policy statement of the American College of Medical Genetics and Genomics. *Genetics in Medicine* 24, 1407–1414.
- Karczewski, K. J., Francioli, L. C., Tiao, G., et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature* 581, 434–443. (gnomAD.)
- Ioannidis, N. M., Rothstein, J. H., Pejaver, V., et al. (2016). REVEL: An ensemble method for predicting the pathogenicity of rare missense variants. *American Journal of Human Genetics* 99, 877–885.
- Cheng, J., Novati, G., Pan, J., et al. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science* 381, eadg7492.
- Jaganathan, K., Panagiotopoulou, S. K., McRae, J. F., et al. (2019). Predicting splicing from primary sequence with deep learning. *Cell* 176, 535–548. (SpliceAI.)
- Caudle, K. E., Klein, T. E., Hoffman, J. M., et al. (2014). Incorporation of pharmacogenomics into routine clinical practice: the Clinical Pharmacogenetics Implementation Consortium (CPIC) guideline development process. *Current Drug Metabolism* 15, 209–217.
- Manolio, T. A., Collins, F. S., Cox, N. J., et al. (2009). Finding the missing heritability of complex diseases. *Nature* 461, 747–753.
- Havasupai Tribe v. Arizona Board of Regents: <https://digitalcommons.law.yale.edu/cbio/2/>
- All of Us Research Program: <https://www.researchallofus.org/>
- CARE Principles: <https://www.gida-global.org/care>
- ClinVar: <https://www.ncbi.nlm.nih.gov/clinvar/>
- gnomAD: <https://gnomad.broadinstitute.org/>
- PharmGKB: <https://www.pharmgkb.org/>
- CPIC: <https://cpicpgx.org/>
- ClinGen: <https://clinicalgenome.org/>
- ACMG SF list: <https://www.acmg.net/ACMG/Medical-Genetics-Practice-Resources/>
- FDA LDT final rule (2024): <https://www.federalregister.gov/documents/2024/05/06/2024-08935>
- EU IVDR: <https://ec.europa.eu/health/md_sector/overview_en>

---

## Appendix — Timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — From Research to Clinic                                 | 20&nbsp;min | 0:20 |
| Part 2 — ACMG/AMP Variant Classification                          | 50&nbsp;min | 1:10 |
| Part 3 — The Variant Evidence Ecosystem                            | 30&nbsp;min | 1:40 |
| Part 4 — Pharmacogenomics                                           | 30&nbsp;min | 2:10 |
| Part 5 — Incidental Findings                                         | 20&nbsp;min | 2:30 |
| Part 6 — Regulatory Landscape                                        | 25&nbsp;min | 2:55 |
| Part 7 — Ethics, Equity, and Data Sovereignty                         | 25&nbsp;min | 3:20 |
| Wrap-up                                                                | 10&nbsp;min | 3:30 |

**Total:** ~3h 30min of content.
