# Lecture 17 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 17 (Clinical Genomics, Variant Interpretation, and Ethics).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-17.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-17/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript.
- Must render standalone.
- Embedded in the lecture via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` exactly once near the end of `<body>`.** C6 smoke gate.

### 1.2 Visual design

- Tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- Classification colours: P red, LP orange, VUS grey, LB muted cobalt, B cobalt.
- Ancestry palette consistent with L12/L13/L17 lecture: African cobalt, European amber, East Asian red, South Asian green, Admixed violet.
- Regulatory colours: RUO grey, LDT amber, IVD cobalt-dark.
- Typography: Inter for UI chrome; JetBrains Mono for HGVS variant names, gene names, star alleles, allele-frequency values, predictor scores.
- Default state is instructive: opens with a pre-computed example.
- Controls grouped in a panel above or to the left of the visualisation.
- Animations ≤ 400 ms.

### 1.3 Interaction model

- **Checkboxes / radio / sliders / dropdowns** — editable inputs.
- **Apply / Run / Reset** — trigger-based.
- Illegal input → quiet inline message (`--fg-muted`); never an `alert()`.

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question:

- ACMG/AMP Classifier → five-class call + reasoning trace listing which rule combination fired.
- Variant Evidence Dossier Compiler → summary of gnomAD / ClinVar / predictor / splice evidence for a variant; suggested evidence codes.
- PharmGKB Star-Allele Translator → metaboliser phenotype + CPIC-recommended dosing action.
- ACMG SF Incidental-Findings Checker → reporting recommendation (report / don't report) with reason.
- GINA Coverage Explorer → coverage verdict for the scenario + applicable state-law additions.
- FDA / Regulatory Pathway Picker → recommended regulatory path + summary of requirements.
- Ancestry Bias in Predictors → performance gap between ancestries; simulation of "after training-diversity fix".

### 1.5 Feasibility gate on user input (required where input is free-form)

- ACMG/AMP Classifier and SF Checker take structured selections (checkboxes / dropdowns).
- Variant Evidence Dossier accepts an HGVS or chr:pos:ref:alt string; validates format before compiling.
- Other artifacts use dropdowns / radio choices only.

### 1.6 Pedagogical constraint

Every artifact produces its named aha moment. If the student plays with the controls and doesn't land on it, the artifact has failed.

### 1.7 Out of scope

- No accounts, no telemetry, no network calls beyond declared CDN libraries (KaTeX permitted for inline math; otherwise none).
- No external data files > 100 KB; all reference data (simulated gnomAD, mock ClinVar submissions, star-allele tables, GINA state-law data) hardcoded in-browser.
- Artifacts are **educational only** — they explicitly state they are not clinical-grade tools.

---

## 2. Artifact #1 — ACMG/AMP Classifier

**File**: `artifacts/lecture-17/01-acmg-classifier.html`
**Lecture anchor**: §2.4 Rule combination
**EE framing reinforced**: rule-based classification; auditable by design; determinstic given evidence codes.

### Teaching purpose

Student toggles which ACMG/AMP evidence codes apply to a hypothetical variant. The classifier applies the rule table and returns the five-class call plus a reasoning trace showing which rule combination fired.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Evidence codes (toggle each that applies):                  │
│                                                             │
│ Pathogenic:                                                 │
│ [ ] PVS1 (Very Strong)                                      │
│ [ ] PS1 [ ] PS2 [ ] PS3 [ ] PS4 (Strong)                    │
│ [ ] PM1 [ ] PM2 [ ] PM3 [ ] PM4 [ ] PM5 [ ] PM6 (Moderate)  │
│ [ ] PP1 [ ] PP2 [ ] PP3 [ ] PP4 (Supporting)                │
│                                                             │
│ Benign:                                                     │
│ [ ] BA1 (Stand-alone)                                       │
│ [ ] BS1 [ ] BS2 [ ] BS3 [ ] BS4 (Strong)                    │
│ [ ] BP1 [ ] BP2 [ ] BP3 [ ] BP4 [ ] BP5 [ ] BP7 (Supporting)│
├─────────────────────────────────────────────────────────────┤
│ Verdict card:                                               │
│   Classification: LIKELY PATHOGENIC                         │
│   Rule fired: PS + PM + 1× PP                                │
│   Posterior-probability estimate: 92% pathogenic             │
├─────────────────────────────────────────────────────────────┤
│ Reasoning trace:                                            │
│   1. Evaluated "BA1 alone" → false.                         │
│   2. Evaluated "2× BS" → false.                             │
│   3. Evaluated "1× PVS1 + ≥ 1× strong/moderate/supp" → true? No.│
│   4. Evaluated "≥ 2× PS" → false.                           │
│   5. Evaluated "1× PS + ≥ 3× PM" → false.                   │
│   6. Evaluated "1× PS + 1–2× PM" (LP rule) → TRUE.          │
│   → Classification = LP.                                    │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Evidence: PS3 + PM2 + PP3 → LIKELY PATHOGENIC.         │   │
│ │ Add another PM or PS → promotes to PATHOGENIC.         │   │
│ │ Add BA1 (frequency > 5%) → rule conflict → VUS.         │  │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- 28 checkboxes, one per evidence code, grouped by direction and strength.
- Reset button.

### What they see

- Classification verdict card.
- Reasoning trace showing which rules were evaluated and which fired.
- Outcome banner with interpretation.

### Target aha

Classification is deterministic given evidence codes. The rule table is explainable; no magic. The hard problem is assigning codes (professional judgment); combining them is mechanical.

### Technical notes

- Pure JS.
- Implement the ACMG/AMP rule table exactly (from Richards et al. 2015 Table 5).
- Handle rule conflicts (P evidence + BA1, etc.) → return VUS with conflict note.
- Reasoning trace: step through candidate rules in descending-strength order.

### Acceptance criteria

- [ ] Default (no codes) → VUS.
- [ ] PVS1 + PM2 + PP3 → PATHOGENIC.
- [ ] PS3 + PM2 + PP3 → LIKELY PATHOGENIC.
- [ ] BA1 alone → BENIGN.
- [ ] Conflicting codes (e.g. PVS1 + BA1) → VUS with conflict note.
- [ ] Reasoning trace is correct and matches the verdict.
- [ ] Opens pre-rendered with a default mid-strength combination.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — Variant Evidence Dossier Compiler

**File**: `artifacts/lecture-17/02-evidence-dossier.html`
**Lecture anchor**: §3.5 Variant interpretation is a profession
**EE framing reinforced**: multi-evidence aggregation; structured information synthesis.

### Teaching purpose

Paste an HGVS variant or pick a preset. The artifact compiles simulated query results from gnomAD, ClinVar, REVEL, AlphaMissense, SpliceAI, and suggests applicable ACMG/AMP evidence codes.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Variant: [ BRCA1:c.5266dupC (p.Gln1756ProfsTer74) ▾ ]       │
│ Presets: [ MLH1 missense / TP53 hotspot / BRCA1 LoF /       │
│            common benign / VUS splice region ]              │
├─────────────────────────────────────────────────────────────┤
│ Evidence sources (simulated):                               │
│                                                             │
│ 1. gnomAD: allele frequency = 2 × 10⁻⁶ (ultra-rare, 1/500k) │
│    → Suggests: PM2_Supporting                               │
│                                                             │
│ 2. ClinVar: 47 submissions, consensus Pathogenic, 3-star     │
│    → Prior strong evidence                                   │
│                                                             │
│ 3. REVEL: 0.87 → PP3 threshold met                          │
│    AlphaMissense: likely pathogenic (0.72) → PP3 agrees     │
│                                                             │
│ 4. SpliceAI: max score 0.03 (no splice effect)              │
│    → BP7 does not apply (coding change); no splice evidence │
│                                                             │
│ 5. Variant is a canonical-splice or LoF in a gene           │
│    where LoF is known mechanism                              │
│    → PVS1 applicable                                        │
├─────────────────────────────────────────────────────────────┤
│ Suggested evidence codes:                                   │
│   PVS1, PM2_Supporting, PP3                                  │
│ Predicted classification (if these apply): PATHOGENIC       │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ This variant has concordant strong + moderate + supp  │  │
│ │ evidence → meets PATHOGENIC rule. ClinVar 3-star       │   │
│ │ expert-panel agreement corroborates.                    │   │
│ │ Try a ambiguous VUS preset to see conflicting evidence. │  │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Variant preset dropdown (6 presets covering the classification spectrum).
- Or paste-box for HGVS / chr:pos:ref:alt; artifact validates format and picks the closest preset.

### What they see

- Multi-source evidence panel (5 sources).
- Suggested ACMG/AMP codes.
- Predicted classification.
- Outcome banner.

### Target aha

Variant classification is structured information synthesis. Multiple sources contribute independent evidence; consistent signals across sources drive confident calls; conflicting signals drive VUS.

### Technical notes

- Pure JS.
- Six hardcoded variant presets with pre-computed evidence queries.
- Paste input validated against HGVS / chr:pos regex; on failure → quiet inline message.
- Evidence-code suggestion mapping: standard threshold table (REVEL ≥ 0.7 → PP3, etc.).

### Acceptance criteria

- [ ] Each preset produces a plausible evidence panel.
- [ ] BRCA1 LoF preset → PVS1 + PM2 → PATHOGENIC.
- [ ] Common benign preset → BA1 → BENIGN.
- [ ] VUS preset → conflicting or insufficient → VUS.
- [ ] Paste box accepts valid HGVS; rejects invalid.
- [ ] Opens pre-rendered with BRCA1 LoF default.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — PharmGKB Star-Allele Translator

**File**: `artifacts/lecture-17/03-pharmgkb-star-allele.html`
**Lecture anchor**: §4.3 Star-allele nomenclature
**EE framing reinforced**: PGx as genotype → star-allele → phenotype → dosing lookup.

### Teaching purpose

Pick a gene and a genotype (two star alleles). The artifact translates to metaboliser phenotype and surfaces the CPIC-recommended dosing action for a relevant drug.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Gene: [ CYP2D6 ▾ / CYP2C9 / TPMT / HLA-B ]                   │
│ Allele 1: [ *1 ▾ ]   Allele 2: [ *4 ▾ ]                     │
│ Drug (one of several relevant): [ codeine ▾ ]                │
├─────────────────────────────────────────────────────────────┤
│ Translation chain:                                           │
│   CYP2D6 *1/*4 → activity score: 1.0 + 0.0 = 1.0             │
│   Metaboliser phenotype: INTERMEDIATE                       │
├─────────────────────────────────────────────────────────────┤
│ CPIC recommendation for codeine:                            │
│   Use codeine at standard label dose; monitor for analgesic │
│   efficacy and side effects.                                │
│   Alternative opioids (morphine, hydromorphone) acceptable. │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Intermediate metaboliser: standard dose acceptable.    │   │
│ │ Change allele 1 to *5 (duplication of *1) → UM         │  │
│ │   → CPIC: AVOID codeine; ultra-rapid morphine formation │  │
│ │ Change to *4/*4 → PM → CPIC: AVOID codeine; no          │   │
│ │   analgesic benefit, use non-CYP2D6-dependent opioid.   │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Gene dropdown: CYP2D6 / CYP2C9 / TPMT / HLA-B.
- Allele 1 + Allele 2 dropdowns: populated with star alleles valid for the chosen gene.
- Drug dropdown: populated with drugs relevant to the chosen gene (codeine for CYP2D6; warfarin for CYP2C9; azathioprine for TPMT; abacavir for HLA-B).

### What they see

- Translation chain: genotype → activity scores → phenotype.
- CPIC dosing recommendation.
- Outcome banner with alternative-allele exploration hints.

### Target aha

PGx is a mechanical lookup once the knowledge bases are populated: genotype → star allele → phenotype → dose. The ML / inference is already done in the evidence aggregation; production is table lookup.

### Technical notes

- Pure JS.
- Hardcoded star-allele activity score tables for CYP2D6, CYP2C9, TPMT (from PharmVar).
- HLA-B handled separately: *57:01 positive or negative; corresponding abacavir action.
- CPIC recommendations hardcoded per (phenotype × drug).

### Acceptance criteria

- [ ] CYP2D6 *1/*1 → normal → standard codeine dose.
- [ ] CYP2D6 *4/*4 → poor metaboliser → avoid codeine.
- [ ] CYP2D6 *1 ×2 (duplication) → ultra-rapid → avoid codeine.
- [ ] HLA-B *57:01 positive → avoid abacavir.
- [ ] TPMT poor metaboliser → dose reduce thiopurine 10×.
- [ ] Opens pre-rendered with CYP2D6 *1/*1 + codeine.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — ACMG SF Incidental-Findings Checker

**File**: `artifacts/lecture-17/04-sf-checker.html`
**Lecture anchor**: §5.4 Incidental findings in practice
**EE framing reinforced**: multi-gate decision (gene-on-list + variant-P/LP + consented).

### Teaching purpose

Simulate the incidental-findings reporting decision. Given a candidate variant (gene + classification + patient consent status), check the three gates and produce a reporting recommendation.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Candidate variant:                                          │
│   Gene: [ BRCA1 ▾ / TP53 / MYH7 / CFTR / ... ]               │
│   Variant: [ preset list specific to the selected gene ]    │
│   Classification: [ P ▾ / LP / VUS / LB / B ]                │
│                                                             │
│ Patient-consent status:                                     │
│   [ opted-in to SF ▾ / opted-out / not asked ]              │
├─────────────────────────────────────────────────────────────┤
│ Gate evaluation:                                            │
│   Gate 1 — Gene on ACMG SF list: BRCA1 ✓                    │
│   Gate 2 — Variant classified P or LP: P ✓                  │
│   Gate 3 — Patient consented: opted-in ✓                    │
│   → ALL THREE GATES PASS                                    │
├─────────────────────────────────────────────────────────────┤
│ Recommendation:                                             │
│   REPORT as Secondary Finding.                              │
│   Route through genetic counsellor before clinical          │
│   release; include in "Secondary Findings" section of       │
│   the clinical report.                                      │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Change classification to VUS → Gate 2 fails →          │   │
│ │   DO NOT REPORT (variant isn't actionable).            │   │
│ │ Change consent to "opted-out" → Gate 3 fails →          │  │
│ │   DO NOT REPORT (respects patient preference).          │  │
│ │ Change gene to a non-SF gene → Gate 1 fails.            │  │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Gene dropdown: selection from ACMG SF v3.x gene list + a few non-SF examples.
- Variant preset per gene.
- Classification dropdown (five classes).
- Consent dropdown (opted-in / opted-out / not asked).

### What they see

- Three-gate evaluation with pass/fail per gate.
- Reporting recommendation.
- Outcome banner.

### Target aha

Incidental-findings reporting is a structured multi-gate decision. All three gates must pass. Explicit rule-based decision process; defensible in audit.

### Technical notes

- Pure JS.
- ACMG SF v3.x gene list hardcoded (81 genes).
- Classification gates: P and LP pass; VUS/LB/B fail.
- Consent gates: opted-in passes; opted-out fails; not-asked produces a "requires consent clarification" warning.

### Acceptance criteria

- [ ] BRCA1 + P + opted-in → REPORT.
- [ ] BRCA1 + VUS + opted-in → DO NOT REPORT (Gate 2 fails).
- [ ] BRCA1 + P + opted-out → DO NOT REPORT (Gate 3 fails).
- [ ] Non-SF gene + P + opted-in → DO NOT REPORT (Gate 1 fails).
- [ ] Default opens with BRCA1 + P + opted-in (the canonical YES case).
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — GINA Coverage Explorer

**File**: `artifacts/lecture-17/05-gina-coverage.html`
**Lecture anchor**: §7.1 GINA and its limits
**EE framing reinforced**: GINA is narrower than commonly assumed; life / disability / LTC are the big gaps.

### Teaching purpose

Student picks a scenario (combination of domain, entity size, and state). Artifact returns whether GINA covers, whether state law fills the gap, and the residual risk.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Scenario:                                                   │
│   Genetic finding: [ BRCA+ ▾ / HTT expansion (Huntington) / │
│                      APOE ε4 / HLA-B*57:01 / ... ]           │
│   Domain: [ health insurance ▾ / life insurance /            │
│             disability / long-term care / employment ]      │
│   Entity size (if employment): [ N/A ▾ / ≥15 / <15 ]         │
│   State (if US): [ California ▾ / Florida / Vermont /       │
│                    most other states ]                       │
├─────────────────────────────────────────────────────────────┤
│ Evaluation:                                                 │
│   Federal GINA: DOES NOT COVER (life insurance is a         │
│     GINA carve-out)                                         │
│   State law (California): CALGINA partial protection on     │
│     life + disability insurance                              │
│   → Patient has some protection in California; none in      │
│     most other states.                                      │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Life insurance + BRCA+ in most states:                 │   │
│ │   insurer may consider genetic test → higher premium   │   │
│ │   or denial legal under federal law.                    │   │
│ │ In California/Florida: state law restricts.             │  │
│ │ Change domain to "health insurance" → GINA protects.    │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Genetic-finding dropdown.
- Domain dropdown.
- Entity size (when applicable).
- State (US) dropdown.

### What they see

- Coverage verdict (federal GINA + state law).
- Outcome banner explaining implications.

### Target aha

GINA is narrower than commonly assumed. Health insurance and large-employer employment yes; life, disability, long-term care, small employers no. State-law patchwork partially fills gaps in a handful of states.

### Technical notes

- Pure JS.
- Hardcoded GINA scope table.
- State-law additions for California, Florida, Vermont, Oregon, Washington, Massachusetts, Connecticut (summarised).

### Acceptance criteria

- [ ] Health insurance + any state → GINA protects.
- [ ] Life insurance + most states → no federal protection, minimal state coverage.
- [ ] Life insurance + California → partial state protection.
- [ ] Employment + small employer → GINA does not protect (< 15 employees).
- [ ] Default opens on life insurance + BRCA + most-states (illustrates the biggest gap).
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — FDA / Regulatory Pathway Picker

**File**: `artifacts/lecture-17/06-fda-pathway.html`
**Lecture anchor**: §6.1 FDA-cleared sequencing assays; §6.2 LDT rule; §6.3 IVDR
**EE framing reinforced**: regulatory pathway follows intended use; picking wrong adds years of delay.

### Teaching purpose

Describe a proposed assay (test type, clinical intent, market, sample size for validation). Artifact maps to the likely regulatory path: RUO / LDT / 510(k) / PMA / IVDR Class C/D.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Assay description:                                          │
│   Test type: [ targeted panel (oncology) ▾ / exome /        │
│                 WGS / PGx panel / NIPT / ... ]               │
│   Clinical intent: [ research ▾ / LDT (single-lab) /        │
│                      commercial kit / companion diagnostic ] │
│   Market: [ US only ▾ / EU / both / global ]                 │
│   Target patient population size: [──●── 1000 ]              │
├─────────────────────────────────────────────────────────────┤
│ Recommended pathway:                                        │
│   US: FDA PMA (Class III; CDx for targeted therapy)          │
│   EU: IVDR Class C                                          │
│                                                             │
│ Summary of requirements:                                    │
│   • Clinical performance studies (≥ 1000 specimens)         │
│   • Quality system (21 CFR 820 / ISO 13485)                  │
│   • Notified Body review (EU)                               │
│   • Labeling + cybersecurity documentation                   │
│   • Estimated timeline: 2–4 years                            │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Companion-diagnostic oncology panel → full PMA / IVDR  │  │
│ │   class-C. Heavy documentation and clinical studies.   │   │
│ │ Change "companion diagnostic" → "LDT single-lab" →     │   │
│ │   FDA LDT pathway (phased 2024+); much faster.         │   │
│ │ Change "research only" → RUO labelling; no FDA review. │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Test-type dropdown.
- Clinical-intent dropdown.
- Market dropdown.
- Patient-population-size slider.

### What they see

- Recommended US + EU regulatory pathways.
- Summary of requirements + timeline estimate.
- Outcome banner with path-switching illustrations.

### Target aha

Regulatory pathway follows the assay's intended use, not the lab's preference. Companion diagnostics are at the highest regulatory bar; LDTs are becoming regulated (2024 FDA rule); RUO is research-only with explicit labelling.

### Technical notes

- Pure JS.
- Hardcoded decision table: (test-type × clinical-intent × market) → (pathway, requirements, timeline).
- Requirements text generated from preset strings.

### Acceptance criteria

- [ ] Oncology CDx + US + EU → PMA + IVDR Class C.
- [ ] Single-lab LDT + US → FDA LDT (phased 2024+).
- [ ] Research-only → RUO.
- [ ] Exome test + commercial kit + US → 510(k) or de novo.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Artifact #7 — Ancestry Bias in Variant Predictors

**File**: `artifacts/lecture-17/07-ancestry-bias.html`
**Lecture anchor**: §7.2 Ancestry bias in clinical databases
**EE framing reinforced**: bias is a training-distribution coverage problem; fix is data, not algorithm.

### Teaching purpose

Simulate REVEL / AlphaMissense performance on variants stratified by ancestry. Show the accuracy gap. Toggle a "retrained on diverse corpus" switch to see how the gap narrows with better training.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Training corpus composition:                                │
│   Eur. fraction: [──●── 0.80 ]                              │
│   African:       [──●── 0.06 ]                              │
│   E. Asian:      [──●── 0.10 ]                              │
│   S. Asian:      [──●── 0.02 ]                              │
│   Admixed Am.:   [──●── 0.02 ]                              │
│ [Re-simulate evaluation]                                    │
├─────────────────────────────────────────────────────────────┤
│ Predictor: [ REVEL ▾ / AlphaMissense ]                      │
├─────────────────────────────────────────────────────────────┤
│ Per-ancestry F1 on held-out variants:                       │
│   European:         0.88                                    │
│   East Asian:       0.84                                    │
│   South Asian:      0.82                                    │
│   Admixed American: 0.82                                    │
│   African:          0.78                                    │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Current Euro-dominated training: 10pp F1 gap vs African.│  │
│ │ Set training corpus to equal 20%-per-ancestry →         │   │
│ │   F1 gap shrinks to ~2pp across populations.            │   │
│ │ Bias follows coverage: balance the data, close the gap. │  │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Five ancestry-fraction sliders (must sum to 1; UI auto-normalises).
- Predictor dropdown (REVEL / AlphaMissense).
- Re-simulate button.

### What they see

- Per-ancestry F1 bar chart.
- Outcome banner showing gap size and the effect of balanced training.

### Target aha

Ancestry bias is quantifiable and follows training-distribution coverage. Balancing the training corpus closes the gap mechanistically — it's a data problem, not a modelling flaw.

### Technical notes

- Pure JS.
- Simulate per-ancestry F1 as a function of training-corpus fraction using a calibrated curve (approximates published disparity magnitudes).
- Bar chart rendered in SVG.

### Acceptance criteria

- [ ] Default (Euro-dominated) → 10pp gap between European and African F1.
- [ ] Equal-fraction corpus → gap < 3pp.
- [ ] Extreme Euro-only corpus → gap > 15pp.
- [ ] Sliders auto-normalise to sum to 1.
- [ ] Opens pre-rendered with Euro-dominated default.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 9. Cross-Artifact Consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Classification colour conventions identical across #1, #2, #4.
- Ancestry palette identical across #7 and L12/L13/L17 lecture figures.
- Regulatory colour conventions identical across #6 and L17 figures.
- Every artifact emits an **outcome banner** per convention §1.4.
- All artifacts explicitly state they are educational and not clinical-grade.

## 10. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls.
- [ ] `<script src="../_shared/resize.js" defer></script>` embedded near `</body>`.
- [ ] Outcome banner or equivalent verdict line visible at end of any user interaction.
- [ ] Educational disclaimer visible somewhere in the artifact.
