# Lecture 14 — Figures Specification

> **Scope**: Static diagrams for Lecture 14 (Data Engineering, File Formats, and Reproducibility).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow `diagram-style-guide.md` for visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-14.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG; content is stack / diagram / table heavy rather than plot-heavy.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure legible at 720 px; scales to 1200 px.
- Schematic diagrams (layered stacks, DAGs, tables) dominate; use strong use of `--bg-subtle` boxes with `--fg` borders.
- Colour usage: `--accent` cobalt for "modern / recommended"; `--fg-muted` grey for "legacy / deprecated"; amber `#b45309` for "warning / depends on reference"; teal `#0d7377` for "community standard".
- Typography: Inter for UI labels and annotations; JetBrains Mono for file extensions, tool commands, accession IDs, byte sizes, version numbers.
- Escape `&`, `<`, `>` as XML entities (`&amp;`, `&lt;`, `&gt;`).

## Figure Budget

Twelve figures for a ~3h 30min lecture:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | File format compression comparison | Part 1 | Custom SVG |
| 2 | CRAM compression scheme | Part 1 | Custom SVG |
| 3 | VCF/BCF record structure | Part 1 | Custom SVG |
| 4 | Bioinformatics data repository landscape | Part 2 | Custom SVG |
| 5 | "Which GRCh38" decoy / alt variants | Part 3 | Custom SVG |
| 6 | Bioinformatics dependency hell | Part 4 | Custom SVG |
| 7 | Container vs VM stack | Part 4 | Custom SVG |
| 8 | Workflow DAG example | Part 5 | Custom SVG |
| 9 | Resume-on-failure flow | Part 5 | Custom SVG |
| 10 | GIAB truth set schema | Part 6 | Custom SVG |
| 11 | Precision/recall on GIAB benchmark | Part 6 | Custom SVG |
| 12 | nf-core pipeline curation flow | Part 7 | Custom SVG |

---

## Figure 1 — File format compression comparison

**File**: `diagrams/lecture-14/01-format-comparison.svg`
**Lecture anchor**: §1.1 Why file formats matter
**ViewBox**: `0 0 1080 440`

### Purpose

Visual comparison of file sizes and random-access capability across the major bioinformatics formats for a typical 30× human WGS sample.

### Content

**Top panel — horizontal bar chart of file sizes.** Five rows, each a horizontal bar:

- FASTQ.gz — bar length proportional to 45 GB; colour `--fg-muted`.
- BAM — 80 GB; `--accent-soft`.
- CRAM (default) — 25 GB; `--accent`.
- CRAM (+ quality binning) — 18 GB; `--accent` with amber highlight.
- VCF.gz — 1.5 GB; teal.
- BCF — 0.5 GB; teal dark.

Label each bar with its format name (left), its byte size (right), and its "random-access" icon (✓ / ✗ / index).

**Bottom caption strip.** Four-column mini-table:

| Format | Compression ratio vs raw | Random access | Tool ecosystem |
| ------ | ------------------------ | ------------- | -------------- |
| FASTQ.gz | 3× | ✗ | universal |
| BAM | 2× | ✓ (`.bai`) | mature |
| CRAM | 5–7× | ✓ (`.crai`) | mature |
| VCF.gz | 4× | ✓ (`tabix`) | universal |
| BCF | 10× | ✓ (`.csi`) | growing |

### Style notes

- Bars aligned on left edge; length proportional.
- Row labels in JetBrains Mono.
- Random-access icons in a narrow trailing column.

---

## Figure 2 — CRAM compression scheme

**File**: `diagrams/lecture-14/02-cram-scheme.svg`
**Lecture anchor**: §1.2 FASTQ, BAM, CRAM
**ViewBox**: `0 0 1080 480`

### Purpose

Show how CRAM encodes a read as (position, diffs from reference, compressed quality) rather than storing the full read sequence. Contrast with BAM's generic encoding.

### Content

**Top row — BAM record.** A full-width rectangle labelled "BAM record (~300 bytes)". Inside, a stylised byte dump showing: read-name, flags, position, CIGAR, 150-char sequence "ACGTACGT…", 150-char quality "####IIIIII…". Annotation at right: "generic zlib compression on the full block".

**Middle row — reference strip.** A horizontal bar representing the reference genome with the aligned region highlighted. Annotation: "reference.fa (pre-loaded at decode time)".

**Bottom row — CRAM record.** Same-width rectangle labelled "CRAM record (~60 bytes)". Inside, four sub-blocks:

- `pos: chr7:12345` (ref position).
- `diffs: [80:T>C, 122:-AGG, ...]` (sparse, variable-length).
- `quality: [compressed, optionally binned]`.
- `flags, name: columnar codec`.

Arrow from reference strip into "diffs block" showing: "decode: fetch reference slice → apply diffs → recover read sequence".

**Footer warning bar.** `--warning-strong` horizontal band: "⚠ Lose the reference → data is unreadable. Record the reference MD5 in the CRAM header."

### Style notes

- BAM record: `--fg-muted` fill.
- CRAM record: `--accent-soft` fill with distinct blocks shaded.
- Reference arrow: amber dashed.
- Monospace throughout for the byte-level content.

---

## Figure 3 — VCF/BCF record structure

**File**: `diagrams/lecture-14/03-vcf-structure.svg`
**Lecture anchor**: §1.3 VCF and BCF
**ViewBox**: `0 0 1200 480`

### Purpose

Show the anatomy of a VCF line and the equivalent BCF encoding side-by-side.

### Content

**Top band — VCF header.** A three-line stylised header:

- `##fileformat=VCFv4.3`
- `##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">`
- `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`
- `#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE1  SAMPLE2  …`

**Middle band — one annotated VCF row** (tab-delimited). Each field boxed and labelled:

```
chr7    12345    rs123    A    G    45    PASS    AF=0.13;AN=1000    GT:DP:GQ    0/1:28:99    1/1:30:99
```

Annotations pointing to each field from below: "chromosome", "position", "variant ID", "reference allele", "alt allele", "phred quality", "filter status", "INFO key-value", "FORMAT schema", "per-sample fields".

**Bottom band — BCF binary encoding**. Same row represented as typed binary fields with widths labelled: `chrom_id: uint8`, `pos: int32`, `alt_count: uint8`, `ref_len: uint8`, `alt_seq: packed`, etc. Annotation: "same data, ~3× smaller, ~5× faster parse".

**Sidebar annotation**. "bgzipped + tabix = seek by coordinate in O(1). Essential at cohort scale."

### Style notes

- Header lines in muted grey.
- VCF row fields in alternating tinted backgrounds for readability.
- BCF binary fields in small monospace with typed annotations.
- Connecting arrows from VCF row → BCF binary row pointing to "BCF is the same schema, encoded".

---

## Figure 4 — Bioinformatics data repository landscape

**File**: `diagrams/lecture-14/04-repositories.svg`
**Lecture anchor**: §2.1 The public data ecosystem
**ViewBox**: `0 0 1200 520`

### Purpose

A map-style overview of the major public data repositories, split into open-access and controlled-access.

### Content

**Left half — "open-access" band.** Four boxes in a row:

- **SRA** (NIH) — accession `SRP*` / `SRR*` / `SRS*`. "Raw reads."
- **ENA** (EBI) — accession `ERP*` / `ERR*` / `ERS*`. "Mirror of SRA."
- **GEO** (NIH) — accession `GSE*` / `GSM*`. "Processed data (count matrices, microarray)."
- **ArrayExpress** (EBI) — "GEO counterpart; weaker curation."

Arrows between SRA ↔ ENA showing mirroring.

**Right half — "controlled-access" band.** Three boxes:

- **dbGaP** (NIH) — "US individual-level genotype + phenotype."
- **EGA** (EBI/CRG) — "European controlled; UK Biobank, FinnGen."
- **GDC** (NCI) — "Cancer; TCGA, TARGET."

Small locked-padlock icons on the controlled-access boxes.

**Top band — "cloud-hosted copies"**. A thin horizontal band above both halves showing an AWS/GCS icon row with: `s3://1000genomes`, `gs://gcp-public-data--gnomad`, `s3://nasa-gcs-mirror` (etc.). Label: "cloud mirrors = faster bulk access; cost model per-bucket (open vs requester-pays)".

**Bottom caption**. "Every 2025 paper cites an accession ID from one of these. Know which your data is in before you write the data-availability statement."

### Style notes

- Open-access boxes: `--accent-soft` fill.
- Controlled-access boxes: `--warning-soft` fill with padlock icons.
- Cloud band: teal, narrower.
- Accession-prefix labels in monospace.

---

## Figure 5 — "Which GRCh38" decoy / alt variants

**File**: `diagrams/lecture-14/05-which-grch38.svg`
**Lecture anchor**: §3.1 Why "which GRCh38" is a question
**ViewBox**: `0 0 1080 440`

### Purpose

Show how the "GRCh38" label hides four progressively-larger bundle variants with different tool expectations.

### Content

**Four horizontal stacked layers**, each progressively wider to represent larger size:

1. **Primary** — 3.0 GB. Label: "chrs 1–22, X, Y, MT only".
2. **+ Unplaced contigs** — 3.05 GB. Label: "plus ~300 unassembled fragments".
3. **+ ALT loci** — 3.15 GB. Label: "plus ~250 alternative haplotypes (e.g. HLA)".
4. **+ Decoys + HLA alts** — 3.40 GB. Label: "plus ~2300 artificial decoys + HLA ALTs. *GRCh38_full_analysis_set_plus_decoy_hla.fa*".

Each layer annotated on the right with "Compatible tools":

1. "BWA default, STAR default, most academic tools."
2. "Any ALT-unaware tool."
3. "BWA-MEM with `.alt` file, ALT-aware callers."
4. "DeepVariant recommended, Broad pipelines, most 2024+ clinical."

**Bottom band — amber warning strip.** "Mixing these produces silent errors. Record reference MD5 in BAM/CRAM `@SQ M5:` fields; check on every pipeline handoff."

### Style notes

- Each layer in successively darker cobalt.
- Right-side tool-compatibility column in JetBrains Mono.
- MD5 / warning band: full-width amber.

---

## Figure 6 — Bioinformatics dependency hell

**File**: `diagrams/lecture-14/06-dependency-hell.svg`
**Lecture anchor**: §4.1 Why dependency hell is bad in bioinformatics
**ViewBox**: `0 0 1080 520`

### Purpose

Visualise the multi-layer, mixed-language dependency stack that makes a "simple" bioinformatics pipeline painful to install.

### Content

**Exploded-view stack**, bottom to top:

- Kernel layer (linux-6.5, glibc-2.39).
- OS library layer (libcurl, libssl, zlib, libhts).
- Language runtimes (Python 3.11, R 4.3, Perl 5.36).
- Domain libraries (NumPy, pysam, Rsamtools, Bioconductor, AnnData).
- Tools (BWA-MEM 0.7.17, samtools 1.19, GATK 4.5, DeepVariant 1.6, Nextflow 23.10).
- User pipeline code (Snakemake / Nextflow config).

**Red-slash collision indicators** at three points:

- "libhts 1.15 vs 1.17" (tool A wants one, tool B wants the other).
- "Bioconductor 3.18 breaks Bioconductor 3.15 packages."
- "CUDA 11.8 vs 12.1" (GPU tool vs default).

**Right-side annotation box.** Title "Container isolation". Text: "Docker / Apptainer pin the whole stack kernel-upward. Each tool in its own container → collisions cannot happen across tools."

### Style notes

- Stack layers in successive tints.
- Collision indicators: red `/` diagonal slashes across affected layers.
- Right box: `--accent-soft` border.

---

## Figure 7 — Container vs VM stack

**File**: `diagrams/lecture-14/07-container-vm.svg`
**Lecture anchor**: §4.2 Containers: Docker, Singularity, Apptainer
**ViewBox**: `0 0 1080 440`

### Purpose

Classic side-by-side stack comparison. Highlight that containers share the host kernel while VMs virtualise all the way down to hardware.

### Content

**Left column — VM stack.** Bottom to top: "Hardware" → "Hypervisor" → "Guest OS kernel" → "Guest OS userland" → "Application (BWA)". Width: wider, labelled "~1 GB overhead, multi-second boot".

**Right column — Container stack.** Bottom to top: "Hardware" → "Host OS kernel" → "Container runtime (Docker / Apptainer)" → "Application (BWA)". Width: narrower, labelled "<10 MB overhead, instant start".

**Bottom strip — bioinformatics use cases.** Three small container boxes: "BWA + samtools + GATK", "AlphaFold + CUDA 12.1", "scVI + PyTorch". Each labelled with its image size: "~800 MB", "~6 GB", "~4 GB". Annotation: "one image per tool = independent updates".

### Style notes

- VM column: `--fg-muted` tinted blocks.
- Container column: `--accent` tinted blocks.
- Kernel layer highlighted in teal (shared in containers).
- Size labels in JetBrains Mono.

---

## Figure 8 — Workflow DAG example

**File**: `diagrams/lecture-14/08-workflow-dag.svg`
**Lecture anchor**: §5.1 Why workflow languages exist
**ViewBox**: `0 0 1200 560`

### Purpose

Show a realistic pipeline as a DAG of processes, with parallelism annotations and resource requirements per node.

### Content

**DAG with ~10 nodes**, arranged roughly left-to-right:

- `fastp` (QC, trim) — per-sample, 4 cpus / 4 GB.
- `bwa-mem` (align) — per-sample, 16 cpus / 32 GB / 4 h.
- `samtools sort` — per-sample, 8 cpus / 16 GB.
- `samtools markdup` — per-sample, 4 cpus / 8 GB.
- `samtools calmd` — per-sample, 4 cpus / 4 GB.
- `deepvariant` (call) — per-sample × 22 shards → **parallel fan-out annotation** (22 copies of the node).
- `bcftools concat` (merge shards) — per-sample, 2 cpus / 4 GB.
- `bcftools filter` — per-sample, 2 cpus / 2 GB.
- `annotation` (VEP) — per-sample, 4 cpus / 8 GB.
- `multiqc` (report) — cohort-level, 2 cpus / 4 GB.

**Node visuals.** Rounded rectangles; node name in Inter; resource badge in monospace at bottom-right of each node.

**Edges.** Directed arrows; labelled with data product (`reads.fq.gz`, `aligned.bam`, `dedup.bam`, etc.).

**Annotation boxes**:

- Top-left: "Parallelism: 22 shards × 100 samples = 2200 concurrent tasks maximum."
- Bottom: "Workflow manager schedules against cluster / cloud; resume on failure via intermediate caching."

### Style notes

- Per-sample nodes in `--accent-soft`.
- Cohort-level nodes in teal.
- Parallel-fan-out node: stacked shadow (3 rectangles offset) to suggest multiplicity.
- Resource badges in monospace with CPU/mem/time.

---

## Figure 9 — Resume-on-failure flow

**File**: `diagrams/lecture-14/09-resume-logic.svg`
**Lecture anchor**: §5.3 Resume, retry, resource management
**ViewBox**: `0 0 1080 480`

### Purpose

Show how a workflow manager uses per-step caches to skip completed work on resume.

### Content

**Two side-by-side timelines**.

**Left — Run 1 (fails).**

- Step 1 (`fastp`) — green, 5 min, caches its output.
- Step 2 (`bwa-mem`) — green, 120 min, caches.
- Step 3 (`samtools markdup`) — **red X**, fails at 30 min. Cache update skipped.
- Steps 4, 5 — greyed-out, never run.

Cache state snapshot below: green checkmarks for step 1 + step 2 hashes.

**Right — Run 2 with `-resume`.**

- Step 1 — grey "CACHED"; 0 s elapsed.
- Step 2 — grey "CACHED"; 0 s elapsed.
- Step 3 — green, re-runs, succeeds (user fixed the bug), caches.
- Step 4 — green, runs.
- Step 5 — green, runs.

**Footer caption.** "Cache keyed on (input file hash × tool version × command line × resource spec). Bit-for-bit identical invocations re-use the cached output."

### Style notes

- Red X for failure.
- Green boxes for successful runs.
- Grey boxes for cache hits.
- Dashed arrows connecting to the per-step cache state representation.

---

## Figure 10 — GIAB truth set schema

**File**: `diagrams/lecture-14/10-giab-schema.svg`
**Lecture anchor**: §6.1 GIAB: the ground-truth datasets
**ViewBox**: `0 0 1200 420`

### Purpose

Show what a GIAB "truth set" physically is: a reference-coordinate track of confident regions plus a VCF of truth variants within them.

### Content

**Reference coordinate track (top)**. A horizontal bar representing a ~5 Mb region. Below the bar, annotated sub-regions:

- Confident region BED (smooth grey bands covering ~95% of the length).
- Excluded regions (hatched, ~5%): segmental duplications, centromeric repeats, low-complexity.

**Truth VCF track (middle)**. Tick marks within the confident bands representing individual truth variants. Colour per variant type: SNV (cobalt), small INDEL (amber), large INDEL / SV (violet).

**Legend / stats panel (right)**. Small table:

- Truth variants in this region: 12,400 SNV / 1,250 INDEL / 45 SV.
- Confident region: 4,758,000 bp of 5,000,000 bp (95.2%).
- Excluded: 242,000 bp (4.8%) — mostly repeats.

**Annotation at bottom.** "A variant call outside the confident region is **not evaluated** by GIAB — precision/recall is undefined there. These regions are where current technology is weakest."

### Style notes

- Bar: light grey fill.
- Confident sub-regions: filled pale cobalt.
- Excluded regions: cross-hatched.
- Variant ticks: thin vertical lines colour-coded by variant type.

---

## Figure 11 — Precision/recall on GIAB benchmark

**File**: `diagrams/lecture-14/11-precision-recall.svg`
**Lecture anchor**: §6.2 hap.py and precision/recall
**ViewBox**: `0 0 1080 520`

### Purpose

A precision-vs-recall scatter comparing several variant callers on GIAB HG002, with per-variant-type breakdown.

### Content

**Main plot.** X-axis: recall, 0.9–1.0. Y-axis: precision, 0.9–1.0. Iso-F1 contour lines drawn lightly at F1 = 0.95, 0.97, 0.99, 0.999.

**Dots** (each a caller × variant-type pair):

- DeepVariant (short reads) — SNV (cobalt, top-right, F1 ≈ 0.999), INDEL (amber, lower at F1 ≈ 0.98).
- GATK HaplotypeCaller — SNV (slightly lower than DV), INDEL (further below).
- DeepVariant PacBio — SNV (~F1 0.998), INDEL (~F1 0.996, catches up because HiFi reads help indel calling).
- Sniffles (SVs) — well below in both axes, F1 ≈ 0.90 (SVs are harder).

Per-dot annotations with caller name + variant type.

**Inset table (bottom)**. Top-5 rows showing caller / variant type / precision / recall / F1. Monospace.

### Style notes

- Iso-F1 contours: faint grey curves.
- Dots: 5 px filled circles; colour by variant type.
- Legend at top-right.
- Annotations in small italic.

---

## Figure 12 — nf-core pipeline curation flow

**File**: `diagrams/lecture-14/12-nf-core-flow.svg`
**Lecture anchor**: §5.4 Nextflow's nf-core community
**ViewBox**: `0 0 1080 440`

### Purpose

Show how an nf-core pipeline goes from contributor PR to release, highlighting the CI/CD + community-review machinery that distinguishes it from ad-hoc pipelines.

### Content

**Horizontal flow with 5 stations**:

1. **Contributor PR** — user opens a pull request on GitHub.
2. **Lint + test CI** — `nf-core lint` + pipeline CI runs the pipeline on test data; checks style, config, containers.
3. **Community review** — reviewers in the pipeline's Slack channel.
4. **Merge + tag** — PR merges, version tag created.
5. **Release artefact** — Docker/Singularity image built; DOI assigned via Zenodo; announcement in nf-core Slack.

Each station a rounded rectangle with an icon (PR / CI / review / tag / release).

**Below** the flow, a thin band representing infrastructure: GitHub Actions, BioContainers registry, Zenodo archive, Slack community.

**Right-side summary box**. Top-4 nf-core pipeline names (`rnaseq`, `sarek`, `atacseq`, `scrnaseq`) with their star counts / download estimates.

**Annotation**. "Community-curated ≠ vendor-curated. Pipelines stay alive because dozens of users care."

### Style notes

- Flow arrows between stations.
- Infrastructure band: `--bg-subtle` fill.
- Icons drawn inline with clean geometric shapes.
