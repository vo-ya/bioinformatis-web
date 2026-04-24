# Lecture 14 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 14 (Data Engineering, File Formats, and Reproducibility).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-14.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-14/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript.
- Must render standalone.
- Embedded in the lecture via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` exactly once near the end of `<body>`.** C6 smoke gate.

### 1.2 Visual design

- Tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- File-format colours: FASTQ grey, BAM amber-muted, CRAM cobalt, VCF teal, BCF teal-dark, HDF5 violet.
- Status colours: success green `#2d7a3e`, warning amber `#b45309`, error red `#c4342c`.
- Container / env / tool blocks: `--accent-soft` rounded rectangles.
- Workflow DAGs: nodes are rounded rectangles; edges directed with arrow markers.
- Typography: Inter for UI chrome; JetBrains Mono for accession IDs, file names, commands, byte sizes, hash digests.
- Default state is instructive: opens with a pre-computed example.
- Controls grouped in a panel above or to the left of the visualisation.
- Animations ≤ 400 ms.

### 1.3 Interaction model

- **Sliders / toggles / dropdowns / text-paste boxes** — editable parameters with sensible ranges and presets.
- **Run / Simulate / Reset** — for trigger-based computations.
- **Step** (DAG walker) — for iterative animations.
- Illegal input → quiet inline message (`--fg-muted`); never an `alert()`.

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question:

- File-Format Size and Compression Explorer → estimated sizes per format at the current coverage / sample count; cohort-total size.
- SRA / GEO / EGA Accession Parser → repository / resource type / download URL / metadata URL.
- Workflow DAG Visualiser → which nodes are cached vs re-run under a given failure scenario; cache hit/miss counts.
- Container vs Conda Isolation Demo → which tool version wins in each scenario; collision vs isolation verdict.
- GIAB Precision/Recall Simulator → P/R/F1 per variant type; confusion matrix; threshold-scan curve.
- VCF/BCF Record Explorer → parsed fields with types; size comparison VCF vs BCF for the same content.
- Reference Build Chooser → recommended reference variant + justification; expected byte size; compatible-tools list.

### 1.5 Feasibility gate on user input (required where input is free-form)

- Accession Parser + VCF Record Explorer accept pasted text → validate structure; report line-level errors in a quiet inline box.
- Other artifacts use slider / dropdown / preset controls only.

### 1.6 Pedagogical constraint

Every artifact produces its named aha moment. If the student plays with the controls and doesn't land on it, the artifact has failed.

### 1.7 Out of scope

- No accounts, no telemetry, no network calls beyond declared CDN libraries (none here — all logic and presets are local).
- No external data files > 100 KB (per-artifact hardcoded reference tables, preset VCF fragments, accession patterns).

---

## 2. Artifact #1 — File-Format Size and Compression Explorer

**File**: `artifacts/lecture-14/01-format-sizes.html`
**Lecture anchor**: §1.6 A quick note on scale
**EE framing reinforced**: content-aware compression (CRAM) vs generic (BAM); storage cost grows linearly with cohort size.

### Teaching purpose

Slide sequencing parameters (coverage, read length, samples, quality-binning toggle). See the resulting size per format and cumulative cohort total. Show that CRAM + quality binning is ~4× smaller than BAM, and that a 100k-sample biobank stores petabytes differently in different formats.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Coverage:          [──●── 30× ]                             │
│ Read length:       [──●── 150 bp ]                          │
│ Samples:           [──●── 1000 ]  (log)                     │
│ Quality binning:   [ ✓ 8-bin QV ]                           │
│ CRAM reference:    [ ✓ matched ]                            │
├─────────────────────────────────────────────────────────────┤
│ Per-sample size bar chart:                                  │
│   FASTQ.gz  ████████████████████   45 GB                    │
│   BAM       █████████████████████  80 GB                    │
│   CRAM      █████                   25 GB                    │
│   CRAM + QV ███                     18 GB                    │
│   VCF.gz                            1.5 GB                   │
│   BCF                               0.5 GB                   │
├─────────────────────────────────────────────────────────────┤
│ Cohort total (1000 samples):                                │
│   FASTQ.gz   45 TB                                          │
│   BAM        80 TB                                          │
│   CRAM       25 TB                                          │
│   CRAM + QV  18 TB                                          │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ CRAM + QV vs BAM at 30× × 1k samples: saves 62 TB       │   │
│ │ Estimated S3 storage: $1400/mo (BAM) → $310/mo (CRAM+QV)│   │
│ │ Turn coverage to 60× → all sizes double.                │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Coverage (×): 10 – 200, default 30.
- Read length (bp): 50 – 300, default 150.
- Sample count: 1 – 1,000,000 (log slider), default 1000.
- Quality binning toggle (default on).
- CRAM reference matched toggle (default on; off → warn that data is unreadable).

### What they see

- Per-sample size bar chart across formats.
- Cohort total card.
- Outcome banner with comparative savings and an illustrative cloud-storage cost line.

### Target aha

At 30× and 1000 samples, CRAM + QV saves ~62 TB over BAM. At 100k samples this becomes 6.2 PB. Student sees: **format choice is a real cost decision**, not aesthetics.

### Technical notes

- Pure JS.
- Formulas approximate:
  - FASTQ raw ≈ coverage × 3 Gb × 2 (two strands for paired-end) × 1 byte/base + quality bytes → compressed ratio ~3×.
  - BAM ≈ 1.8 × FASTQ.gz (aligned, richer metadata).
  - CRAM ≈ 0.55 × BAM (default), 0.35 × BAM (quality binning).
  - VCF.gz ≈ 4 Mb per sample per chromosome at typical variant density.
  - BCF ≈ 0.3 × VCF.gz.
- Use hard-coded constants per format; no external network.
- Cloud-cost line: use a flat price per TB per month.

### Acceptance criteria

- [ ] Default (30×, 150 bp, 1k samples) shows CRAM at ~25 GB/sample and cohort total ~25 TB.
- [ ] Toggling QV binning shifts CRAM down ~30%.
- [ ] Toggling CRAM reference off → warning banner: "CRAM without matching reference = unreadable".
- [ ] Changing samples to 100k scales cohort total 100×.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — SRA / GEO / EGA Accession Parser

**File**: `artifacts/lecture-14/02-accession-parser.html`
**Lecture anchor**: §2.3 Accessions, metadata, and citation
**EE framing reinforced**: accessions as URIs; repository prefixes as a namespace.

### Teaching purpose

Paste any bioinformatics accession → identify repository, parse the type (project / run / sample / study), generate canonical URLs for metadata and data access. Teaches the accession-prefix conventions.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Accession: [ SRP123456                          ]           │
│ [Parse]  [Try: SRP... / GSE... / PRJNA... / phs.v...]        │
├─────────────────────────────────────────────────────────────┤
│ Parsed:                                                     │
│   Repository:      SRA (NIH)                                 │
│   Type:            Study / Project                            │
│   Mirror in ENA:   ERP counterpart (if matched)              │
│   Access model:    Open                                      │
├─────────────────────────────────────────────────────────────┤
│ Generated URLs:                                              │
│   Metadata (runinfo CSV):    https://trace.ncbi.nlm.nih.gov/… │
│   Data (via sra-tools):      prefetch SRP123456              │
│   ENA mirror:                https://www.ebi.ac.uk/ena/…    │
│   Cloud mirror:              s3://sra-pub-src-1/…            │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ SRP123456 → SRA (NIH) project; 42 runs; 1.2 TB total   │   │
│ │ Mirrored in ENA as ERP789012                           │   │
│ │ Bulk download via `prefetch` (CLI) or AWS s3 sync       │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Text input: accession (any prefix).
- Preset-button row: SRP / SRR / SRS / GSE / GSM / PRJNA / phs / EGAD.
- Parse button (also parses on Enter).

### What they see

- Parsed accession table.
- Generated URL list.
- Outcome banner with interpretation.

### Target aha

Every accession is a URI with a meaningful prefix: SRP = SRA project, SRR = run, GSE = GEO study, PRJNA = BioProject, phs = dbGaP, EGAD = EGA dataset. Student can identify each type at a glance after playing with several presets.

### Technical notes

- Pure JS.
- Hardcoded prefix table with regex-based parsing.
- URLs generated deterministically (no network calls); NCBI, ENA, and EBI URL structures are stable.
- Validation: if the accession doesn't match any known prefix, display an inline "unrecognised — known prefixes: …" message.

### Acceptance criteria

- [ ] Preset buttons auto-populate a valid accession and parse correctly.
- [ ] Pasted valid accessions parse to correct repository and type.
- [ ] Invalid accession produces a quiet inline error (no alert).
- [ ] Generated URLs click-through to actual resources (spot-check a few).
- [ ] Opens pre-rendered with a default SRP example.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — Workflow DAG Visualiser

**File**: `artifacts/lecture-14/03-workflow-dag.html`
**Lecture anchor**: §5.3 Resume, retry, resource management
**EE framing reinforced**: workflow = dataflow DAG; cache hits / misses governed by input-hash × tool-version × command-line.

### Teaching purpose

Load a preset (variant-calling, RNA-seq, or simple toy) DAG. Simulate a failure at a chosen node. Show which nodes are cached, which must re-run. Illustrate what a workflow manager actually automates.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Preset pipeline: [ variant-call ▾ / rna-seq / toy ]         │
│ Nodes / samples: [──●── 5 samples ]                         │
│ Simulate failure at: [ markdup ▾ / call / annot / none ]    │
│ Fix failure: [ ☐ apply fix and resume ]                     │
│ [Run]  [Reset caches]                                       │
├─────────────────────────────────────────────────────────────┤
│ DAG visualisation (SVG):                                    │
│   nodes coloured by state (green = ran, grey = cached,      │
│     red = failed, amber = pending)                          │
│   edges labelled with data product                          │
│   per-sample fan-out visible                                 │
├─────────────────────────────────────────────────────────────┤
│ Cache state table:                                          │
│   step    | hash     | state                                │
│   fastp   | 0xA1B2…  | cached ✓                             │
│   bwa-mem | 0xC3D4…  | cached ✓                             │
│   markdup | —        | failed ✗                             │
│   …                                                         │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Run 1: failed at markdup after 5 samples × (fastp+bwa). │   │
│ │  Cache has 10 entries (5 fastp + 5 bwa-mem).            │   │
│ │ Run 2 with -resume: 10 cache hits, 5 markdup + downstream│   │
│ │  re-runs. Savings: 120 min → 30 min.                    │   │
│ │ Change tool version on markdup → invalidates cache      │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Preset dropdown (variant-call / rna-seq / toy).
- Sample count: 1 – 20, default 5.
- Failure-injection dropdown (node name or none).
- "Apply fix and resume" toggle.
- Run / Reset buttons.

### What they see

- DAG with coloured node states.
- Cache-state table.
- Outcome banner with time savings estimate.

### Target aha

On resume, only the failed step + downstream re-run; upstream cache hits are free. Change a tool version on the failed step → cache invalidates and upstream still hits. Student sees: caching is keyed on inputs + tool, not on time.

### Technical notes

- Pure JS.
- SVG-based DAG renderer; simple layered layout.
- Preset DAGs hardcoded with 6–10 nodes each.
- Cache state: dictionary keyed by node-id; simulated hash is `djb2(name + version + sample)`.
- Time estimates: hardcoded per-node.

### Acceptance criteria

- [ ] Default variant-call preset renders correctly.
- [ ] Simulating a failure at `markdup` and resuming shows correct cache hits.
- [ ] Changing sample count re-renders parallel fan-out.
- [ ] Cache-invalidation on tool-version change demonstrated.
- [ ] Opens pre-rendered with a default fail-and-resume scenario.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — Container vs Conda Environment Isolation

**File**: `artifacts/lecture-14/04-containers.html`
**Lecture anchor**: §4.2 Containers; §4.3 Conda, mamba, pixi
**EE framing reinforced**: containerisation = OS-level isolation; conda = userland-level isolation; host-install = no isolation.

### Teaching purpose

Simulate a shared Linux host trying to run three tools with conflicting dependencies. Compare three isolation strategies: (a) all on host, (b) each in a conda env, (c) each in a container. See which resolves the conflicts.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Tools to run (preset collisions):                           │
│   ☑ BWA 0.7.17 → needs libhts 1.15                          │
│   ☑ bcftools 1.19 → needs libhts 1.19                       │
│   ☑ custom-script.py → needs Python 3.8                      │
│   (host has Python 3.11, libhts 1.17)                        │
├─────────────────────────────────────────────────────────────┤
│ Isolation strategy:                                          │
│ [ host ▾ / conda / container ]                              │
├─────────────────────────────────────────────────────────────┤
│ Simulation:                                                 │
│   Three running processes shown as blocks with their        │
│   effective environment (Python version, libhts version,    │
│   tool versions).                                           │
│   Collisions marked with red stripes.                       │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ host: BWA and bcftools both fail (libhts mismatch).    │   │
│ │ conda: three envs; each tool uses its own libhts; works.│   │
│ │ container: three containers; each with its own full stack│   │
│ │   incl. pinned kernel-level support; works.             │   │
│ │ Kernel shared across containers; host userland unchanged │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Tool-set checkboxes (three preset tools).
- Isolation strategy dropdown: host / conda / container.
- Toggle: "attempt GPU-dependent tool" (adds AlphaFold-style CUDA dependency; increases difficulty).

### What they see

- Three tool blocks with their effective environments.
- Collision indicators when conflicts exist.
- Outcome banner with pass/fail per strategy.

### Target aha

Host strategy fails; conda resolves userland collisions; containers resolve userland + some system-library collisions. Student sees the progression of isolation granularity.

### Technical notes

- Pure JS.
- Hardcoded dependency graph per tool.
- Simulate strategy by resolving whether each tool's dependency graph is satisfiable under the current isolation scheme.
- Conda = separate Python + userland libraries but shared system binaries.
- Container = full userland + OS libraries; shared kernel only.

### Acceptance criteria

- [ ] Default (host strategy, all three tools) → collisions flagged.
- [ ] Switching to conda → collisions resolved.
- [ ] Switching to container → same resolution; plus GPU tool now works.
- [ ] Opens pre-rendered showing the host-failure case.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — GIAB Precision/Recall Simulator

**File**: `artifacts/lecture-14/05-giab-precision-recall.html`
**Lecture anchor**: §6.2 hap.py and precision/recall
**EE framing reinforced**: variant-calling benchmarks as detection theory; per-variant-type calibration matters at scale.

### Teaching purpose

Simulate a caller's output against a GIAB-like truth set. Slide a quality threshold. Show per-variant-type precision, recall, F1; confusion matrix; PR curve across thresholds. Illustrate how a single F1 figure hides per-variant-type failure modes.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Caller profile:                                             │
│   SNV sensitivity (recall ceiling):  [──●── 0.999 ]          │
│   SNV false-positive rate:            [──●── 0.001 ]          │
│   INDEL sensitivity:                  [──●── 0.97  ]          │
│   INDEL false-positive rate:          [──●── 0.02  ]          │
│   Quality threshold:                  [──●── 20    ]          │
│ [Re-simulate 10k variants]                                  │
├─────────────────────────────────────────────────────────────┤
│ Confusion matrix (per variant type):                        │
│   SNV:   TP=9800  FP=10    FN=200   P=0.999  R=0.98  F1=0.99│
│   INDEL: TP=900   FP=20    FN=100   P=0.978  R=0.90  F1=0.94│
├─────────────────────────────────────────────────────────────┤
│ PR curve (threshold scan):                                  │
│   SNV + INDEL curves overlaid                                │
│   Current-threshold marker                                   │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Aggregate F1 = 0.985; SNV F1 = 0.99; INDEL F1 = 0.94.   │   │
│ │ The aggregate hides an indel-recall gap.                 │   │
│ │ Lower threshold to 10 → recall up, precision down.      │   │
│ │ Raise to 40 → recall down, precision up.                 │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- SNV sensitivity / FPR sliders.
- INDEL sensitivity / FPR sliders.
- Quality threshold slider.
- Re-simulate button.

### What they see

- Per-variant-type confusion matrix.
- PR curve across thresholds, per variant type.
- Outcome banner with aggregate + breakdown.

### Target aha

SNV and INDEL metrics diverge; a strong aggregate F1 can hide indel weakness. Changing the threshold trades precision for recall along the PR curve. Student sees: **per-type breakdown is non-negotiable** in variant-calling benchmarks.

### Technical notes

- Pure JS, seeded PRNG.
- Simulate 10k truth variants (90% SNV, 10% INDEL).
- Each caller output: draw TP with probability sensitivity; inject FPs at specified rate.
- Assign per-variant "quality" from a mixture distribution (TPs high, FPs low).
- Confusion matrix recomputes instantly on threshold change.
- PR curve computed across 20–30 threshold values.

### Acceptance criteria

- [ ] Default → SNV F1 > 0.99, INDEL F1 ~ 0.94.
- [ ] Raising quality threshold → precision up, recall down (both types).
- [ ] Lowering INDEL sensitivity → INDEL recall drops; SNV unchanged.
- [ ] PR curves match the thresholded confusion-matrix outputs.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — VCF/BCF Record Explorer

**File**: `artifacts/lecture-14/06-vcf-explorer.html`
**Lecture anchor**: §1.3 VCF and BCF
**EE framing reinforced**: VCF is self-describing text; BCF is the same schema as binary; partial access requires bgzip + tabix indexing.

### Teaching purpose

Paste or select a small VCF text. Parse it row-by-row. Display each field with its type and meaning. Show the same data re-encoded as BCF (conceptually — size comparison, not bit-level output).

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ VCF input:                                                   │
│   Presets: [ tiny ▾ / realistic / multi-sample / malformed ]  │
│   [paste or edit VCF text here]                              │
│ [Parse]                                                      │
├─────────────────────────────────────────────────────────────┤
│ Parsed header:                                              │
│   ##fileformat=VCFv4.3                                      │
│   2 INFO fields, 3 FORMAT fields, 5 samples                 │
├─────────────────────────────────────────────────────────────┤
│ Per-record view (scrollable):                                │
│   Row 1: chr7 12345 rs123 A G 45 PASS AF=0.13 ...           │
│     Each field boxed + typed                                 │
│     INFO key-values expanded                                │
│     Per-sample GT/DP/GQ rendered inline                     │
├─────────────────────────────────────────────────────────────┤
│ Size comparison:                                             │
│   This text VCF:   2.4 KB                                    │
│   Gzipped:         0.8 KB                                    │
│   Equivalent BCF:  0.3 KB (estimated)                        │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ 4 records parsed. 2 INFO fields, 3 FORMAT fields.      │   │
│ │ 1 malformed record flagged at line 12 (missing FORMAT). │   │
│ │ BCF would encode this ~8× smaller than raw text.        │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Preset dropdown (tiny / realistic / multi-sample / malformed).
- Editable textarea.
- Parse button (also on Enter).

### What they see

- Parsed header summary.
- Per-record boxed rendering with typed fields.
- Size-comparison card.
- Outcome banner + any line-level errors.

### Target aha

VCF is human-readable text with a structured schema. BCF is the same schema in binary. Random access via tabix is an indexing story, not a format story. Student recognises that **schema = contract**; VCF vs BCF are two serializations of the same content.

### Technical notes

- Pure JS.
- Implement a minimal VCF parser: header lines (starting with `##`), column line (starting with `#CHROM`), data lines (tab-delimited).
- For INFO: split on `;`, then split each on `=` to key/value.
- For FORMAT + samples: split on `:` in lockstep.
- Size estimates: text length (character count); gzip ~3×; BCF ~ text/8.
- Malformed records: any row with field-count mismatching the header.

### Acceptance criteria

- [ ] Default preset renders correctly with 4+ records parsed.
- [ ] Per-field boxes appear with correct types.
- [ ] Malformed preset flags the specific line.
- [ ] Size-comparison card updates on text edits.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Artifact #7 — Reference Build Chooser

**File**: `artifacts/lecture-14/07-reference-chooser.html`
**Lecture anchor**: §3.1 Why "which GRCh38" is a question
**EE framing reinforced**: reference-variant choice is a deliberate tool-driven decision; MD5 pinning is mandatory.

### Teaching purpose

Given a tool + pipeline goal, recommend which GRCh38 variant to use. Display byte size, compatible tools, and known incompatibilities. Let student pick a tool (BWA, DeepVariant, STAR, etc.) and see the recommendation update.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Pipeline goal: [ germline WGS ▾ / somatic tumour /          │
│                  RNA-seq / ChIP-seq / pangenome build ]      │
│ Primary aligner: [ BWA-MEM ▾ / STAR / minimap2 / Bowtie2 ]   │
│ Variant caller:  [ DeepVariant ▾ / GATK / Sniffles /  — ]    │
│ ALT-aware:       [ ☑ ]                                       │
├─────────────────────────────────────────────────────────────┤
│ Recommendation panel:                                       │
│   Use: GRCh38_full_analysis_set_plus_decoy_hla.fa           │
│   Size: ~3.4 GB                                             │
│   Reason: DeepVariant's recommended pipeline expects         │
│     decoys + HLA alts; without them, variants in HLA region │
│     are miscalled.                                          │
│   MD5: 64b32de2fc934679c16e83a2bc072064                     │
├─────────────────────────────────────────────────────────────┤
│ Compatible tools: ✓                                         │
│   BWA-MEM (with .alt file), DeepVariant, GATK               │
│ Incompatible: ✗                                             │
│   STAR (wants Ensembl release, different chromosome names)   │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Germline WGS + DeepVariant → GRCh38 + decoys + HLA      │   │
│ │ For RNA-seq with STAR, switch to Ensembl release N       │   │
│ │   (different naming convention).                        │   │
│ │ Always pin by MD5, not by filename.                      │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Pipeline goal dropdown.
- Aligner dropdown.
- Variant caller dropdown.
- ALT-aware checkbox.

### What they see

- Recommended reference variant with filename, size, MD5.
- Compatible / incompatible tool lists.
- Outcome banner with reasoning.

### Target aha

Reference choice is not arbitrary. Germline WGS with DeepVariant → decoys + HLA alts. RNA-seq with STAR → Ensembl-style release. Somatic-tumour workflows often want primary-only to keep alt contigs from confusing the caller. Student sees: **pipeline goal + tools drive reference choice**, and mismatches produce silent bugs.

### Technical notes

- Pure JS.
- Hardcoded decision table: (goal × aligner × caller) → recommended reference + reasoning.
- Reference metadata (size, MD5, filename) from the NCBI / EBI / iGenomes canonical bundles.
- Compatible / incompatible lists derived from each tool's documentation.

### Acceptance criteria

- [ ] Default (germline WGS + BWA + DeepVariant) → recommends decoy+HLA bundle.
- [ ] Switching to RNA-seq + STAR → recommends Ensembl release.
- [ ] Switching to somatic tumour + GATK → recommends primary-only variant.
- [ ] Toggling ALT-aware off → warning about losing HLA accuracy.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 9. Cross-Artifact Consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Format / repository / tool colour conventions consistent across artifacts (§1.2 above).
- JetBrains Mono used for all accession IDs, file names, hash digests, byte counts.
- Every artifact emits an **outcome banner** per convention §1.4.
- Preset dropdowns standardised with consistent labelling.

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
