# Lecture 4 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 4 (Variant Calling).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-04.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

These conventions apply to every artifact in this lecture. Per-artifact sections below override them only when they need to.

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-04/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript. External libraries only if justified in the per-artifact section.
- The file must render standalone when opened directly in a browser.
- Artifact is embedded in the lecture page via `<iframe>` loaded lazily.

### 1.2 Visual design

- Use the design tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- DNA bases, when rendered, use the `--base-*` palette exactly.
- Typography: **Inter** for UI chrome; **JetBrains Mono** for sequences, pileups, VCF fields, genotypes.
- Default state is instructive: the artifact opens showing a meaningful example, no user input required.
- Controls grouped in a panel above or to the left of the visualization.
- No animations longer than ~400 ms. Motion only when it carries information.

### 1.3 Interaction model

- **Input / sequences / thresholds** — editable text fields, sliders, or dropdowns, validated against the artifact's input alphabet.
- **Step / Play / Pause / Reset** — where the artifact shows an algorithm with discrete steps.
- **Speed** — optional slider, 0.25×–4×. Default 1×.
- Illegal input shows a quiet inline message (`--fg-muted`), not a modal.

### 1.4 Explicit outcome reporting (required)

Every artifact in this lecture answers its own question at the end. The student should never be left to infer the result from the final animation state. Concretely:

- If the artifact calls a variant, it shows a **banner**: "Variant called: chr1:1000 A>G, GT 0/1, QUAL 420" or "No variant — read support below threshold".
- If the artifact scores a genotype likelihood, it shows **the MAP genotype and log-likelihood ratio**.
- If the artifact parses / filters a VCF, it shows **input-row count vs output-row count** and why rows were dropped.
- If the artifact annotates a variant, it shows a **verdict line**: "BRAF V600E · Pathogenic · Rare in population".

### 1.5 Feasibility gate on user input (required where input is free-form)

Artifacts that accept user input (pileup columns, VCF text, pileup files) must pre-flight the input and report *why* it is or isn't runnable before the calculation runs. Rejected inputs get an inline explanation, not a silent blank state.

### 1.6 Pedagogical constraint

Every artifact must produce a **specific realization** — the "target aha moment" named in its section. If the student plays with the artifact and doesn't land on that realization, the artifact has failed and should be revised.

### 1.7 Out of scope

- No logins, accounts, or persistence between sessions.
- No telemetry or analytics.
- No external data files larger than ~50 KB.

---

## 2. Artifact #1 — Pileup Viewer

**File**: `artifacts/lecture-04/01-pileup-viewer.html`
**Lecture anchor**: §2.2 The pileup
**EE framing reinforced**: pileup as a per-column histogram; reading alignment evidence by eye.

### Teaching purpose

Turn the `samtools mpileup` format from opaque text into a visual, column-by-column representation. The student should leave able to eyeball a pileup line and tell whether the column shows a candidate variant.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Mode: [Preset ▾] [Custom pileup]                            │
│ Preset: [ clean region · heterozygous SNV · homozygous SNV ·│
│          INDEL cluster · low-depth region ▾ ]              │
├─────────────────────────────────────────────────────────────┤
│ Raw pileup text (scrollable, monospace):                    │
│   chr1  1000  A  12  ,,..,A,..,aA,     IIIII!IIIIII         │
│   chr1  1001  C  12  ,,...,.,,,,.     IIIII!IIIII           │
│   chr1  1002  G  12  ,,..aA,Aa,.,     IIIII!IIIII           │
│   ...                                                       │
├─────────────────────────────────────────────────────────────┤
│ Visual pileup — each position as a vertical stack of reads: │
│   [ grid: x = position, y = stacked reads ]                 │
│   columns auto-highlighted where VAF > user-set threshold   │
│   threshold slider: [ VAF ≥ 0.20 ]                          │
├─────────────────────────────────────────────────────────────┤
│ Candidate variant columns (auto-detected):                  │
│   · chr1:1000  A→G  DP=12  VAF=0.167  ← marginal, low cov   │
│   · chr1:1002  G→A  DP=12  VAF=0.333  ← likely het          │
│   · chr1:1005  -  INDEL  DP=10                              │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Preset dropdown**: 5 preset pileup excerpts representing different biological scenarios.
- **Custom pileup textarea**: paste a samtools mpileup output fragment.
- **VAF threshold slider**: 0–1, default 0.2.
- **Strand visualization toggle**: colour reads by forward/reverse strand.

### What they see

- Left pane: raw pileup text, colour-coded (matches in neutral, mismatches in `--base-*`, INDELs in `--warning`).
- Right pane: visual stacked bars — one bar per read per position, stacked vertically, coloured by base identity.
- Columns where VAF > threshold get a subtle `--accent-bg` column-wide background highlight.
- Bottom: a table of auto-detected candidate variant columns, each with position, REF→ALT, DP, VAF, and a one-line qualitative label ("marginal", "likely het", "likely hom").

### Target aha moment

Load the "INDEL cluster" preset. The naive pileup shows a dense run of mismatches at five adjacent positions. The student should see from the pattern alone — five consecutive low-depth single-base mismatches — that this is a misaligned INDEL rather than five real SNVs. Connects to Figure 7 and Artifact #2.

### Technical notes

- Pure JS. Parse mpileup tokens per-column: strip `^X` and `$`, handle `+N` / `-N` insertion/deletion prefixes.
- Quality string parsed as Phred ASCII (char code − 33).
- Visual rendering: SVG grid, 10 px per read-row, 14 px per position-column.

### Acceptance criteria

- [ ] Default preset shows at least one clearly heterozygous SNV (6 REF / 6 ALT) at a labeled position.
- [ ] "INDEL cluster" preset visibly breaks the "every mismatch is a SNV" assumption.
- [ ] VAF threshold slider updates candidate-column highlights in real time.
- [ ] Custom input: invalid pileup lines are flagged inline with line number + error.
- [ ] Legible at 720 px artifact width.

---

## 3. Artifact #2 — INDEL Realignment Demo

**File**: `artifacts/lecture-04/02-indel-realignment.html`
**Lecture anchor**: §3.1 Pre-calling pipeline — INDEL realignment
**EE framing reinforced**: local optimisation vs per-read greedy alignment.

### Teaching purpose

Show viscerally why a true 2 bp deletion can look like a cluster of 5 false SNVs under naive alignment, and how local realignment fixes it in one pass.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Scenario: [ 2-bp deletion ▾ ] [ 1-bp insertion ] [ 3-bp del │
│            near read end ]                                  │
│ [ Naive alignment ]  [ Local realignment ]                  │
│ Realignment window: [ ±k bp  k = 20 ]                       │
├─────────────────────────────────────────────────────────────┤
│ Reference:   ACGTACGTATTGCCATAGCATGCATAGCATAGC              │
│                       ↑ true deletion position              │
├─────────────────────────────────────────────────────────────┤
│ Reads (with alignment visualization):                       │
│   r1:  ACGTACGTATT  CCATAGCATGC      (after del)            │
│   r2:  ACGTACGTATT   CCATAGCATGCAT                          │
│   r3:   CGTACGTATTG CCATAGCATGCATA                          │
│   ...                                                       │
│                                                             │
│   [toggle shows naive-alignment mismatches vs re-aligned]   │
├─────────────────────────────────────────────────────────────┤
│ Pileup at the deletion window:                              │
│   Naive:      mismatch cluster at positions 12-16           │
│   Realigned:  clean gap at positions 12-13                  │
│                                                             │
│ Caller verdict:                                             │
│   Naive → 5 false SNV calls                                 │
│   Realigned → 1 true DEL call                               │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Scenario dropdown**: 3 scenarios (2-bp deletion mid-read, 1-bp insertion, 3-bp deletion near read end).
- **Mode toggle**: Naive alignment vs Local realignment.
- **Realignment window size**: ±5, ±10, ±20 bp.

### What they see

- Reads drawn as horizontal bars aligned to the reference; base mismatches coloured `--error`, gap characters `--fg-muted`.
- Toggling between naive and realigned alignment animates read re-positioning with a gap opening where the deletion should be.
- A simulated "caller output" panel at the bottom counting SNVs and INDELs under each alignment mode.

### Target aha moment

Toggle from "Naive" to "Local realignment" on the "3-bp deletion near read end" scenario. The naive alignment shows dispersed mismatches because each read expressed the deletion differently. Realignment converges them to a single 3-bp gap, and the caller goes from "12 false SNVs" to "1 true DEL." The student should connect this to why modern haplotype-aware callers (HaplotypeCaller, DeepVariant) don't need a separate realignment step — they do this internally.

### Technical notes

- Precompute naive alignments (mis-gapped) and the canonical realigned solution for each scenario — don't run a real Smith-Waterman at runtime; just toggle between two baked representations.
- Reads simulated from a reference sequence with known-true variants.

### Acceptance criteria

- [ ] Each scenario shows a clearly different pattern under the two alignment modes.
- [ ] "Caller verdict" counts update on toggle with correct numbers.
- [ ] Window-size slider changes the region displayed without breaking the visualization.

---

## 4. Artifact #3 — Genotype Likelihood Calculator

**File**: `artifacts/lecture-04/03-genotype-likelihoods.html`
**Lecture anchor**: §3.2 Caller families — Bayesian
**EE framing reinforced**: posterior estimation over a discrete hypothesis set.

### Teaching purpose

Make the Bayesian calculation concrete. Given a small pileup (bases and Phred qualities), compute `P(D|0/0)`, `P(D|0/1)`, `P(D|1/1)`, apply a prior, and show the MAP genotype.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Pileup at one position (editable):                          │
│   REF = A                                                   │
│   Reads:  A A A A A A A G G G  (10 reads)                   │
│   Quals:  30 30 30 30 30 30 30 30 30 30  (all Q30)          │
│   [+ add read] [randomize] [preset: clean het · mixed het ·│
│                  low-depth · somatic 15% VAF]               │
├─────────────────────────────────────────────────────────────┤
│ Likelihoods (log scale):                                    │
│   P(D | 0/0) = 10^-10.1                                     │
│   P(D | 0/1) = 10^-3.0   ← MAP                              │
│   P(D | 1/1) = 10^-12.3                                     │
├─────────────────────────────────────────────────────────────┤
│ Prior: [ diploid human · uniform · somatic low-freq ▾ ]     │
│   P(0/0) = 0.999                                            │
│   P(0/1) = 7e-4                                             │
│   P(1/1) = 2e-4                                             │
├─────────────────────────────────────────────────────────────┤
│ Posterior:                                                  │
│   P(0/0 | D) = 5e-12                                        │
│   P(0/1 | D) = 1.0  ← MAP                                   │
│   P(1/1 | D) = 3e-15                                        │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ ✓ Called: 0/1 (heterozygous)                         │   │
│ │   GQ = 95 · log-likelihood margin = 9.3              │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Read bases** — editable list, typed or via preset.
- **Read qualities** — editable list, defaults to Q30, range Q10–Q40.
- **Prior selector** — three options: diploid human (gnomAD-derived), uniform, somatic low-frequency.
- **REF base** — dropdown A/C/G/T.
- **+ Add / − Remove read** — manipulate the pileup size.
- **Presets** — four canned scenarios for quick exploration.

### What they see

- Pileup column rendered visually as stacked coloured bases.
- Likelihood row showing three numerical values on a log scale.
- Prior row showing the three prior probabilities as a small bar chart.
- Posterior row showing the three posterior probabilities, MAP highlighted.
- Outcome banner at the bottom stating the called genotype, genotype quality, and log-likelihood margin vs runner-up.

### Target aha moment

Load the "somatic 15% VAF" preset. With a diploid human prior, the MAP is `0/0` — the caller treats it as a reference call and misses the somatic variant. Switch the prior to "somatic low-frequency." The MAP flips to `0/1`. The student should see that *the prior matters* — somatic callers succeed by adjusting the prior to accommodate low-VAF events that germline callers correctly dismiss.

### Technical notes

- Compute per-read log-likelihoods: for a read with base `b` and quality `q` given true genotype `G`:
  - If G is homozygous with allele `a`: `P(b|G) = 1−ε` if b=a, `ε/3` otherwise, where `ε = 10^(-q/10)`.
  - If G is heterozygous with alleles `a,b`: `P(b|G) = 0.5·P(b|homozygous a) + 0.5·P(b|homozygous b)`.
- Sum log-likelihoods across reads; multiply (add in log) the prior; renormalise for posterior.

### Acceptance criteria

- [ ] Default pileup (7 REF, 3 ALT, all Q30) produces MAP = 0/1 with GQ ≥ 30.
- [ ] All-reference pileup (10 REF, 0 ALT, Q30) produces MAP = 0/0 with high GQ.
- [ ] Prior switch visibly changes the MAP for borderline cases.
- [ ] Invalid input (non-ACGT base, non-numeric quality) shows an inline error.

---

## 5. Artifact #4 — SV Signature Explorer

**File**: `artifacts/lecture-04/04-sv-signatures.html`
**Lecture anchor**: §4.4 SV callers and their tradeoffs
**EE framing reinforced**: anomaly detection across three signal channels.

### Teaching purpose

Let the student toggle through SV types and see how each one produces a distinct pattern in the three signal channels (discordant pairs, split reads, depth). The target aha: SV calling is signal-detection pattern matching.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ SV type: [ Deletion ▾ Insertion · Inversion · Duplication · │
│             Translocation ]                                 │
│ Size: [──●── 3000 bp ]                                      │
│ Zygosity: [ het ▾ hom ]                                     │
│ Read tech: [ Short paired-end ▾ Long reads (HiFi/ONT) ]     │
├─────────────────────────────────────────────────────────────┤
│ Reference:   ═══════════════|====deleted====|═════════════  │
│              ×   x   x   x      (deleted region greyed)     │
│                                                             │
│ Signal tracks:                                              │
│                                                             │
│  Discordant pairs:                                          │
│     ↶    ↶    ↶   (arcs spanning breakpoint, cluster shown) │
│                                                             │
│  Split reads:                                               │
│     [r]  [r]      (single reads mapping half-left, half-   │
│                    right, with visible break)               │
│                                                             │
│  Read depth:                                                │
│     ─────────╲________________╱─────────                   │
│     30×      10%         30×                                │
├─────────────────────────────────────────────────────────────┤
│ Caller verdict:                                             │
│   Discordant-pair clustering: 14 pairs · ✓ signal           │
│   Split-read support: 6 reads · ✓ single-base breakpoint    │
│   Depth drop: 50% · ✓ heterozygous deletion                 │
│   → SV call: DEL chr1:12500-15500 · homozygous              │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **SV type**: 5 options (Deletion, Insertion, Inversion, Duplication, Translocation).
- **SV size**: slider, 100 bp to 100 kb.
- **Zygosity**: het or hom.
- **Read tech**: short paired-end or long reads (toggles the figure to show the long-read equivalent patterns: a single long read spans the entire SV).
- **Depth**: slider, 10× to 60×.

### What they see

- Reference region with the SV marked as a greyed/striped band.
- Three stacked signal tracks: discordant-pair arcs, split-read half-alignments, depth-of-coverage curve.
- Track counts update with SV size and zygosity — a 200 bp deletion produces fewer discordant pairs than a 5 kb deletion.
- Caller verdict at the bottom composing the three signals into a called SV with confidence.

### Target aha moment

Load a 500 bp heterozygous deletion in short-read mode. Only 2–3 discordant pairs, 1 split read, ~50% depth drop over ~500 bp — marginal signal on all three channels. Switch to long-read mode: one long read spans the whole deletion and shows a clean 500 bp gap. The student realizes short-read SV calling for small events is hard because all three channels are weak; long reads replace the multi-read pattern with a single direct observation.

### Technical notes

- All signal patterns precomputed per (SV type, size, zygosity, read tech) combination; no real alignment at runtime.
- Depth curve is a simple step function.
- Discordant pairs and split reads are positioned stochastically within the breakpoint region to give a "real data" feel.

### Acceptance criteria

- [ ] Each SV type shows a distinct pattern across the three channels.
- [ ] Zygosity change visibly halves the depth drop (hom→het).
- [ ] Long-read toggle replaces the three-channel pattern with a single-read direct observation.
- [ ] Caller verdict at the bottom corresponds to the rendered signals.

---

## 6. Artifact #5 — VCF Parser and Filter

**File**: `artifacts/lecture-04/05-vcf-parser.html`
**Lecture anchor**: §5.1 The VCF file format
**EE framing reinforced**: VCF as a sparse differential encoding; filtering as signal thresholding.

### Teaching purpose

Take raw VCF text, parse it into a structured table, apply filters, and report the kept-vs-dropped counts. Student should leave able to read a real VCF and sketch a filtering pipeline.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ VCF source: [ Preset ▾ ] [ Paste your own ]                 │
│ Preset: [ Small germline · Multi-allelic · Multi-sample trio│
│           · Somatic with low-VAF variants ▾ ]               │
├─────────────────────────────────────────────────────────────┤
│ Parsed table (sortable by QUAL, DP, GT):                    │
│  CHROM  POS   REF  ALT  QUAL   FILTER  DP  GT   VAF   Gene  │
│  chr1   1000  A    G    420    PASS    17  0/1  0.29  BRCA1 │
│  chr1   1500  C    T    612    LowQual 8   0/0  0.12  BRCA1 │
│  chr1   2100  G    A,C  85     PASS    12  1/2  —     BRCA2 │
│  ...                                                        │
├─────────────────────────────────────────────────────────────┤
│ Filters:                                                    │
│   QUAL ≥ [──●── 30 ]                                        │
│   DP ≥ [──●── 10 ]  DP ≤ [──●── 100 ]                        │
│   FILTER = PASS only: [✓]                                   │
│   Genotype: [ any ▾ het only · hom alt only ]                │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Input: 47 variants · Kept: 31 · Dropped: 16           │   │
│ │ Dropped by: QUAL=3 · DP=7 · FILTER=4 · GT=2           │   │
│ └──────────────────────────────────────────────────────┘   │
│ [ Download filtered VCF ]                                   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Source**: preset dropdown or paste-your-own textarea.
- **Filters**: QUAL ≥ slider, DP range, FILTER PASS-only toggle, genotype dropdown.
- **Sort**: click any table column header.
- **Download**: exports the kept rows as a VCF.

### What they see

- Parsed table with one row per variant, columns from the mandatory 8 + selected FORMAT/INFO fields.
- Multi-allelic rows shown expanded (one row per ALT allele where meaningful).
- Kept-vs-dropped outcome banner updating in real time as filters change.
- Dropped-by breakdown so the student sees which filter removed what.

### Target aha moment

Start from the "Somatic with low-VAF variants" preset. Default filters (QUAL ≥ 30, DP ≥ 10) reject 40% of input. Drop the QUAL threshold to 15 and watch the banner: kept variants jump, but the dropped breakdown shows it's mostly PASS-labeled FILTER column — filters from the caller were more informative than raw QUAL. The student should see filtering is a composition of multiple signals, not a single threshold.

### Technical notes

- Pure JS VCF parser: split on tabs, parse INFO as `;`-separated `KEY=VALUE`, parse FORMAT + per-sample columns.
- Handle multi-allelic sites and missing data (`.` / `./.`).
- Download via Blob + object URL.

### Acceptance criteria

- [ ] All four presets parse without errors; table renders with ≥ 20 rows each.
- [ ] Filter sliders update the table in real time.
- [ ] Dropped-by breakdown sums to total-dropped.
- [ ] Multi-allelic preset shows at least one row with comma-separated ALT and GT values like `1/2`.
- [ ] Download produces a valid VCF that the same artifact can re-load.

---

## 7. Artifact #6 — Variant Annotator

**File**: `artifacts/lecture-04/06-variant-annotator.html`
**Lecture anchor**: §5.3 Variant annotation
**EE framing reinforced**: joining a variant table to reference databases.

### Teaching purpose

Run a small variant through the full annotation pipeline — gene context → population frequency → clinical significance → functional prediction — and see the final verdict emerge. Student should understand each annotation layer's role.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Input VCF (small, editable):                                │
│   chr7   140753336  A  T   520  PASS  DP=65   GT=0/1        │
│   chr17  43094692   G  A   410  PASS  DP=52   GT=0/1        │
│   [+ preset: BRAF V600E · BRCA1 pathogenic · TP53 LoF ·     │
│     intergenic · synonymous]                                │
├─────────────────────────────────────────────────────────────┤
│ Annotation pipeline (click a stage to see its output):      │
│                                                             │
│  [ VEP: gene + consequence ] → [ gnomAD: popfreq ] →        │
│  [ ClinVar: significance ] → [ REVEL/CADD: functional ]     │
│                                                             │
│  At each stage, new annotations appended.                   │
├─────────────────────────────────────────────────────────────┤
│ Annotated output:                                           │
│   chr7:140753336 A>T                                        │
│     Gene: BRAF                                              │
│     Consequence: missense_variant                           │
│     Protein change: p.Val600Glu                             │
│     gnomAD_AF: 1.2e-5                                       │
│     ClinVar: Pathogenic (melanoma, colorectal cancer)       │
│     REVEL: 0.932  CADD: 26.5                                │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Verdict banner ──────────────────────────────────────┐   │
│ │ ⚠ BRAF V600E · Pathogenic · Rare in population        │   │
│ │   · High functional impact                            │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Paste VCF** or choose from 5 preset variants representing common clinical scenarios.
- **Step-through**: click each annotation stage to see its specific output appended.
- **Run all**: applies all stages at once.

### What they see

- Input VCF rendered as plain text.
- Four pipeline stages shown as clickable boxes; each adds its annotation to the output panel.
- Annotated output builds up stage by stage.
- Verdict banner at the bottom composes into a human-readable interpretation.

### Target aha moment

Run the "intergenic" preset. VEP assigns `Consequence=intergenic_variant`. gnomAD shows it's common. ClinVar has no entry. REVEL says "not applicable (not in coding)." Verdict: "Common intergenic variant · likely benign." Now compare to the "TP53 LoF" preset: missense → rare → Pathogenic → REVEL 0.95 → Verdict: "High-impact tumor suppressor LoF." Same pipeline, radically different verdicts. The student realizes annotation is where a variant becomes clinically actionable — and that the pipeline composes layered evidence, not just one signal.

### Technical notes

- Pure JS with a small hard-coded mini-database of ~20 known variants (covering BRAF V600E, BRCA1/2 common pathogenics, TP53 hotspots, a few common benign SNPs, intergenic examples, synonymous examples). User variants that don't match a known entry fall back to "not annotated" for clinical stages but still get a consequence prediction.
- No external API calls; all annotation is local mock data.

### Acceptance criteria

- [ ] Each of the 5 preset variants produces a distinct verdict banner.
- [ ] Unknown variants pasted in show "not in mock database" for ClinVar / gnomAD but still get a consequence from the VEP-stage heuristic.
- [ ] Pipeline is click-through one stage at a time; output grows with each click.

---

## 8. Cross-Artifact Consistency

- All six artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- The Pileup Viewer (#1) and Genotype Likelihood Calculator (#3) use consistent pileup rendering (same stacked-bar style, same base colour palette).
- The VCF Parser (#5) and Variant Annotator (#6) use identical VCF-line rendering (monospace, tab-aligned).
- Every artifact emits an **outcome banner** per convention §1.4 — no artifact leaves the student to infer the result from the final animation state.

## 9. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls.
- [ ] Outcome banner or equivalent verdict line visible at the end of any user interaction.
- [ ] User-input artifacts pre-flight inputs with explicit pass/fail messaging.
