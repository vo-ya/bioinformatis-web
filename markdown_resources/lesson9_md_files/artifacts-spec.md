# Lecture 9 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 9 (ChIP-seq, ATAC-seq, and Peak Calling).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-09.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

These conventions apply to every artifact in this lecture. Per-artifact sections below override them only when they need to.

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-09/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript. External libraries only if justified per-artifact.
- Must render standalone when opened directly in a browser.
- Embedded in the lecture page via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` near the end of `<body>`.** C6 smoke gate: exactly one occurrence.

### 1.2 Visual design

- Use tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- Coverage tracks: `--accent` filled area on `--bg-muted` plot background; baseline in `--fg-muted`.
- Annotations / thresholds: `--warning` for thresholds, `--error` for blacklisted/failed, `--success` for passed.
- Typography: **Inter** for UI chrome; **JetBrains Mono** for counts, coordinates, gene/TF names, p-values, PWM weights, nucleotide sequences.
- Default state is instructive: the artifact opens showing a meaningful example, no user interaction required.
- Controls grouped in a panel above or to the left of the visualisation.
- No animations longer than ~400 ms.

### 1.3 Interaction model

- **Sliders / dropdowns / inputs** — editable parameters validated against sensible ranges.
- **Play / Step / Reset** — for iterative processes (peak-caller sweep, footprint accumulation).
- **Re-simulate** — where stochastic data is involved.
- Illegal input shows a quiet inline message (`--fg-muted`), not a modal.

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question at the end:

- Fragment-size inspector → reports the identified peaks (sub / mono / di / tri) and their ratios.
- MACS2 peak caller → reports precision / recall vs planted peaks, and the number of false positives.
- Differential accessibility explorer → reports # of up / down / unchanged peaks at the current FDR.
- PWM motif scanner → reports # of hits above the current threshold and the top scoring position.
- ATAC footprint analyser → reports footprint depth and its SNR as a function of number of aggregated instances.

### 1.5 Feasibility gate on user input (required where input is free-form)

Artifacts accepting user input validate before running: check row formats, matching dimensions, sensible bounds; report rejections inline with line numbers. PWM Motif Scanner (which can accept pasted DNA) must check: only ACGTN characters, length ≥ motif length, not empty.

### 1.6 Pedagogical constraint

Every artifact produces a **specific realization** — the target aha moment named in its section. If the student plays with the artifact and doesn't land on that realization, the artifact has failed.

### 1.7 Out of scope

- No logins, accounts, or persistence between sessions.
- No telemetry.
- No external data files larger than ~100 KB.

---

## 2. Artifact #1 — Fragment-Size Distribution Explorer

**File**: `artifacts/lecture-09/01-fragment-size.html`
**Lecture anchor**: §2.2 Fragment-length distributions
**EE framing reinforced**: histogram shape as QC signal; distribution is the measurement, not just summary.

### Teaching purpose

Visualise how ATAC-seq fragment-size distributions produce the characteristic nucleosomal ladder under healthy chromatin conditions, and how the ladder degrades under bad protocol conditions (over-tagmentation, closed chromatin, wrong size selection). Student sees the ladder as the first QC gate.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Assay: [ ATAC ▾ / ChIP ]                                    │
│ Chromatin state:                                            │
│   Open fraction:    [──●── 0.25 ]                           │
│   Nucleosome positioning regularity: [──●── 0.8 ]           │
│ Protocol parameters:                                        │
│   Tn5 insertion rate: [──●── 1.0 ]                          │
│   Size selection min: [──●── 40 ]   max: [──●── 700 ]        │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Fragment-size histogram (SVG, log-x):                       │
│   clear peaks at sub / mono / di / tri nucleosomal sizes    │
│   ChIP mode shows single smooth mode                        │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ ATAC: 4 peaks detected · ratios 3.1 : 2.0 : 0.8 : 0.3 │   │
│ │ Healthy open-chromatin signature ✓                    │   │
│ │ Mode switch to ChIP → single 300 bp mode, no ladder    │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Assay: ATAC / ChIP (default ATAC).
- Open fraction: 0.0–1.0 (default 0.25) — fraction of genome that is accessible.
- Nucleosome positioning regularity: 0.0–1.0 (default 0.8) — how periodic nucleosome spacing is.
- Tn5 insertion rate: 0.1–3.0 (default 1.0) — over-tagmentation knob.
- Size selection: min 0–300 (default 40), max 200–1000 (default 700).
- Re-simulate.

### What they see

- Histogram in an SVG panel: fragment size on log-x, read density on linear y.
- Peak detection: identified peaks marked with labelled arrows (sub, mono, di, tri).
- Outcome banner states # of detected peaks + inter-peak ratios + QC verdict.

### Target aha moment

Default ATAC → clean laddering visible. Set nucleosome regularity = 0.1 → ladder collapses into a single smooth distribution. Set Tn5 rate = 2.5 → sub-nucleosomal peak dominates, nucleosomal peaks shrink (over-tagmentation). Switch to ChIP → single smooth mode at ~300 bp; no laddering. Student sees the distribution shape as diagnostic.

### Technical notes

- Pure JS, seeded mulberry32.
- ATAC simulation: generate N fragments. For each, place Tn5 insertions on a simulated chromatin model (regular spacing with jitter controlled by regularity parameter); compute the resulting fragment size from distance between insertions. Mixture of sub / mono / di / tri based on geometric drop-off.
- ChIP simulation: fragments ~ log-normal(log(300), 0.4).
- Size selection: trim distribution to [min, max] range.
- Peak detection on the histogram: smooth with Gaussian kernel; find local maxima.

### Acceptance criteria

- [ ] Default config produces visible 4-peak ladder.
- [ ] Setting regularity → 0.1 collapses the ladder.
- [ ] Setting Tn5 rate high skews to sub-nucleosomal.
- [ ] Switching to ChIP shows single-mode distribution.
- [ ] Re-simulate produces visually similar but distinct draws.
- [ ] Opens with default pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — MACS2 Peak Caller

**File**: `artifacts/lecture-09/02-macs2-peak-caller.html`
**Lecture anchor**: §3.2 MACS2 algorithm
**EE framing reinforced**: CFAR detection — local-adaptive threshold.

### Teaching purpose

The central algorithmic artifact of the lecture. Walk a window across a simulated ChIP-seq coverage track; compute the local-Poisson λ and p-value; declare peaks; compare to ground truth (planted peaks).

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulated data:                                             │
│   Planted peak count:       [──●── 20 ]                     │
│   Background rate (per bp): [──●── 0.2 ]                    │
│   Signal strength (peak / bg): [──●── 8.0 ]                 │
│   Include noisy region (artifact-like): [ off ▾ / on ]      │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Algorithm parameters:                                       │
│   Use local lambda: [ on ▾ / off (global only) ]            │
│   Local-window size: [──●── 5000 bp ]                       │
│   FDR threshold:    [──●── 0.05 ]                           │
│ [Run MACS2-style caller]                                    │
├─────────────────────────────────────────────────────────────┤
│ Genomic track (SVG):                                        │
│   coverage visualised in accent                             │
│   planted peak positions marked as small green ticks below  │
│   called peaks above track highlighted in warning           │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Called: 23 peaks · True hits: 19 · Missed: 1 · FP: 4  │   │
│ │ Precision 82.6% · Recall 95.0%                        │   │
│ │ Disabling local λ → FP inflates to 18 in noisy region │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Planted peak count: 5–50 (default 20).
- Background rate: 0.05–1.0 reads/bp (default 0.2).
- Signal strength (peak : bg): 2–20 (default 8.0).
- Include noisy region toggle: adds an artifact region with elevated local coverage but no true signal.
- Use local lambda: on / off — toggling off uses only genome-wide $\lambda_{\text{bg}}$.
- Local-window size: 1000–20000 (default 5000).
- FDR threshold: 0.001–0.1 (default 0.05).
- Re-simulate + Run buttons.

### What they see

- Coverage track on SVG, ~100 kb visualised horizontally.
- Planted peaks marked below the track in `--success` green ticks.
- Called peaks marked above the track as spans in `--warning` amber (true hits) or `--error` red (false positives).
- Outcome banner with precision / recall / comparison numbers.

### Target aha moment

Default → high precision and recall; algorithm works. Enable the noisy region; keep local λ on → clean result. Disable local λ → many false positives in the noisy region because the global $\lambda_{\text{bg}}$ is much lower than the local elevated baseline. Student sees why the local-adaptive lambda is the key CFAR trick.

### Technical notes

- Pure JS, seeded mulberry32.
- Simulate a ~100 kb coverage track. Background: Poisson(λ_bg) per 10 bp bin. Planted peaks: at random positions with Gaussian-shaped coverage enhancement scaled by signal strength. Noisy region: elevated λ on a ~5 kb segment.
- MACS2-style algorithm: slide a 200 bp window; at each candidate, compute local lambda via max of several windows if local is on, else just global. Compute Poisson p-value. BH correction. Merge adjacent significant bins.
- Ground truth: the planted peak positions.
- Precision / recall computed at candidate-position level (±200 bp tolerance).
- Feasibility gate: warn if signal strength < 2 (peaks unresolvable).

### Acceptance criteria

- [ ] Default config yields precision ≥ 0.8 and recall ≥ 0.85.
- [ ] Disabling local λ with noisy region on visibly inflates FP count.
- [ ] Tightening FDR raises precision but drops recall.
- [ ] Planted and called peaks rendered on the same SVG coordinate axis.
- [ ] Opens with default pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — Differential Accessibility Explorer

**File**: `artifacts/lecture-09/03-differential-accessibility.html`
**Lecture anchor**: §4.2 DiffBind, csaw, and counts
**EE framing reinforced**: the L6 NB/GLM toolkit applied to peaks × samples.

### Teaching purpose

Apply a DESeq2-style NB GLM + BH correction to a simulated peaks × samples count matrix. Student explores how replicate count, dispersion, and effect size influence power to detect differential accessibility.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulated experiment:                                       │
│   # peaks:            [──●── 500 ]                          │
│   # replicates per condition: [──●── 3 ]                    │
│   Fraction truly differential: [──●── 0.10 ]                │
│   Effect size (log2FC for DE peaks): [──●── 1.5 ]           │
│   Dispersion (φ): [──●── 0.15 ]                             │
│   Library-size variance across samples: [──●── 0.3 ]        │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Analysis params:                                            │
│   FDR cutoff: [──●── 0.05 ]                                 │
│   log2FC cutoff: [──●── 1.0 ]                               │
├─────────────────────────────────────────────────────────────┤
│ Two-panel output:                                           │
│   (a) MA plot — mean log count × log2 FC                     │
│   (b) Volcano — log2 FC × −log10(adjusted p)                 │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Called: 48 peaks · True DE: 50 · Missed: 7 · FP: 5    │   │
│ │ Power: 0.86 · FDR observed: 10.4%                     │   │
│ │ Drop replicates to 2 → power drops to 0.58            │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- # peaks: 100–2000 (default 500).
- Replicates per condition: 2–8 (default 3).
- Fraction truly differential: 0.0–0.5 (default 0.10).
- Effect size log2FC: 0.0–4.0 (default 1.5).
- Dispersion φ: 0.01–1.0 (default 0.15).
- Library-size variance: 0.0–1.0 (default 0.3).
- FDR cutoff: 0.001–0.2 (default 0.05).
- log2FC cutoff: 0.0–3.0 (default 1.0).
- Re-simulate.

### What they see

- MA plot on left, volcano on right — both SVG.
- Significant points highlighted per cutoffs.
- Outcome banner with power, observed FDR, and a "what if" line like "drop replicates to 2 → power X".

### Target aha moment

Default config → reasonable power (~80%). Drop replicates to 2 → power drops to ~50%. Raise dispersion φ to 0.5 → power drops further. Student sees that replicate count and dispersion are the two knobs that most affect power.

### Technical notes

- Pure JS.
- Simulate per peak per sample: negative-binomial counts with mean determined by peak's base rate × condition multiplier × per-sample library-size factor. Dispersion = φ.
- GLM fit per peak: simple pooled MLE for β_1 (condition effect) using Newton-Raphson on the NB likelihood; or approximate with the Wald statistic (β / SE).
- BH correction across all peaks.
- Call counts (precision, recall, FDR) against the planted truth.

### Acceptance criteria

- [ ] Default config yields power ≥ 0.75, observed FDR within 2× nominal.
- [ ] Lowering replicate count visibly reduces power.
- [ ] Raising dispersion reduces power.
- [ ] Effect size of 0 produces a uniform p-value distribution.
- [ ] MA and volcano plots both update on each run.
- [ ] Opens with default pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — PWM Motif Scanner

**File**: `artifacts/lecture-09/04-pwm-motif-scanner.html`
**Lecture anchor**: §5.3 Motif scanning as a matched filter
**EE framing reinforced**: PWM scan = discrete correlation of 4-channel template against 4-channel one-hot signal.

### Teaching purpose

Interactive PWM scan. Student provides a DNA sequence (or uses a preset), selects a PWM, and sees the scan output as a 1D score track aligned under the sequence. The score spikes mark motif hits.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Sequence input:                                             │
│   [ preset: CTCF peak region ▾ / NFκB peak / random / paste ]│
│   Pasted sequence (if "paste" chosen): [ textarea… ]         │
├─────────────────────────────────────────────────────────────┤
│ PWM selection:                                              │
│   Motif: [ CTCF (JASPAR MA0139.1) ▾ / NFkB / AP1 / custom ] │
├─────────────────────────────────────────────────────────────┤
│ Threshold: [──●── 10.0 (log2 score) ]                       │
├─────────────────────────────────────────────────────────────┤
│ Three stacked bands:                                        │
│   (a) sequence rendered as coloured per-base letters         │
│   (b) score track (sliding-window PWM score)                 │
│   (c) hit markers: small coloured boxes at positions above   │
│       threshold                                             │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Hits above threshold: 3 · top score 17.2 at pos 147   │   │
│ │ At threshold 10.0 : 3 hits on forward, 1 on reverse    │   │
│ │ Random-sequence null: ≈1.2 expected hits per kb        │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Sequence preset dropdown + paste textarea.
- PWM preset dropdown (CTCF / NFκB / AP1 / custom-editable).
- Threshold slider: 0–25 (default 10).
- Scan both strands toggle (default on).

### What they see

- Per-base-coloured sequence band (rendered in JetBrains Mono with base-palette colours).
- Score track below, filled area.
- Hit markers at positions above threshold, with tooltip showing score.
- Outcome banner: hit count, top score + position, random-null expectation.

### Target aha moment

Load CTCF preset → a few hits where expected. Load random sequence → almost no hits. Lower threshold to 5 → many spurious hits everywhere. Student sees that threshold is the knob and the random-null rate tells you when your threshold is too permissive.

### Technical notes

- Pure JS.
- PWM presets: hand-curated approximate JASPAR matrices for CTCF, NFκB (RELA), AP1 (FOS::JUN).
- Scoring: standard log-odds PWM scoring. Background base frequencies = 0.25 each.
- Reverse-complement scan: flip-and-complement the PWM.
- Feasibility gate: check input is ACGTN only, length ≥ L, not empty.
- Random null: estimated from 1000 randomly-permuted sequence samples.

### Acceptance criteria

- [ ] CTCF preset + default threshold produces 2–5 hits.
- [ ] Lowering threshold monotonically raises hit count.
- [ ] Pasted random sequence yields near-zero hits at default threshold.
- [ ] Strand toggle visibly changes reverse-strand hits.
- [ ] Opens with default pre-scanned.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — ATAC Footprint Analyser

**File**: `artifacts/lecture-09/05-atac-footprint.html`
**Lecture anchor**: §5.4 ATAC footprinting
**EE framing reinforced**: signal-averaging over N instances → SNR × √N.

### Teaching purpose

Simulate Tn5 cut-count data at many motif instances; aggregate; show that the footprint emerges from the mean of many noisy individual traces. Student sees the √N scaling empirically.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Simulation:                                                 │
│   True footprint depth: [──●── 0.4 ] (0 = no TF bound)      │
│   Flank elevation: [──●── 1.8 ]                             │
│   Per-site coverage (mean cuts / bp): [──●── 0.3 ]          │
│   Number of motif instances: [──●── 1000 ]                  │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Main panel (SVG):                                           │
│   aggregated footprint — x = pos vs motif centre (-200..200)│
│   y = mean cut count per bp                                 │
├─────────────────────────────────────────────────────────────┤
│ Side panel (6 small plots):                                 │
│   6 randomly-selected individual-instance traces            │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Footprint depth (aggregated): 0.42 · SNR: 14.3        │   │
│ │ Signal visible at N ≥ 250; noisy at N = 50            │   │
│ │ Individual traces: no footprint visible               │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- True footprint depth: 0.0–1.0 (default 0.4; 0 = no binding).
- Flank elevation: 1.0–3.0 (default 1.8).
- Per-site coverage (mean cuts/bp): 0.1–2.0 (default 0.3).
- Number of motif instances: 50–5000 (default 1000).
- Re-simulate.

### What they see

- Main aggregated curve showing clear dip at motif centre + flanking elevation.
- Side panel with 6 individual traces (noisy, typically invisible footprint).
- Outcome banner showing aggregated depth, SNR, and a readout of how SNR scales with N.

### Target aha moment

Default (N=1000) → clean footprint with clear dip. Reduce N to 50 → footprint is lost in noise; aggregated curve is too noisy to see it. Student sees that footprinting is fundamentally about **summing many weak signals** — a classic √N-SNR averaging problem.

### Technical notes

- Pure JS.
- Per-site simulation: for each of N instances, generate a 400 bp cut-count trace as Poisson(per-site coverage × (1 + flank_elevation) for flank positions, per-site coverage × (1 − depth) for centre positions).
- Aggregation: mean across sites → aggregated curve.
- SNR: (aggregated-flank − aggregated-centre) / std of centre region across sites.
- Background null: aggregated curve with footprint depth = 0 for comparison.

### Acceptance criteria

- [ ] Default config produces visible footprint with depth ≈ 0.4 and SNR ≥ 10.
- [ ] Reducing N visibly degrades footprint SNR; at N=50 it's near-invisible.
- [ ] Setting depth to 0 produces flat aggregated curve (no footprint).
- [ ] Individual traces are visibly noisier than aggregate.
- [ ] Opens with default pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — ChIP vs ATAC Coverage Viewer

**File**: `artifacts/lecture-09/06-coverage-viewer.html`
**Lecture anchor**: §1.3 ATAC-seq and DNase-seq — where chromatin is open
**EE framing reinforced**: signal-shape differences across assays.

### Teaching purpose

Walk a ~20 kb genomic window; view simulated ChIP-seq (TF or histone) and ATAC-seq coverage tracks side-by-side; annotate with nucleosome positions, TFBS positions, and promoter/enhancer annotations. Student builds intuition for how signal shape varies by assay.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Region: Chr X : 10,000–30,000 (toy coordinates)              │
│ Assay tracks (toggle on/off):                               │
│   [✓] CTCF ChIP-seq                                         │
│   [✓] H3K4me3 ChIP-seq                                      │
│   [✓] H3K27me3 ChIP-seq                                     │
│   [✓] ATAC-seq                                              │
│ Scenario preset: [ active promoter + TF ▾ / silenced domain │
│   / enhancer cluster / active gene body ]                    │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ SVG browser (500px wide):                                   │
│   - Coordinate axis (kb)                                    │
│   - Each selected track with its own colour and label       │
│   - Annotation row: TSS markers, TFBS markers, nucleosome   │
│     positions                                               │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Scenario: active promoter + TF                        │   │
│ │ CTCF ChIP: 2 sharp narrow peaks at TFBS positions     │   │
│ │ H3K4me3: broad peak spanning promoter (~1.5 kb)       │   │
│ │ H3K27me3: flat background — not silenced here          │   │
│ │ ATAC: strong peak at promoter + secondary at enhancer │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Per-track toggles (4 checkboxes).
- Scenario preset dropdown: "active promoter + TF" / "silenced domain (H3K27me3 high)" / "enhancer cluster (H3K27ac peaks)" / "active gene body".
- Re-simulate.

### What they see

- 4 stacked coverage tracks on a shared genomic axis.
- Annotation row with TSS / TFBS / nucleosome positions.
- Outcome banner with per-scenario interpretation of what each track is showing.

### Target aha moment

Switch between scenarios and observe: in "silenced domain" CTCF and ATAC are flat while H3K27me3 is elevated over a broad region; in "active promoter + TF" the opposite pattern. Student sees how the combination of tracks tells you the regulatory state — no single track alone is sufficient.

### Technical notes

- Pure JS.
- 4 simulated tracks with distinct signal shapes per scenario.
- Each track's coverage is built from pre-planted peaks (narrow for TF, medium for H3K4me3 / H3K27ac, broad for H3K27me3, narrow for ATAC) + background Poisson noise.
- SVG browser: horizontal scroll optional but default view shows full 20 kb window.

### Acceptance criteria

- [ ] Each scenario produces visually distinct pattern across the 4 tracks.
- [ ] Toggling tracks on/off works.
- [ ] Opens with default scenario pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Artifact #7 — CNN Filter Bank Visualiser

**File**: `artifacts/lecture-09/07-cnn-filter-bank.html`
**Lecture anchor**: §6.1 Enformer and its siblings
**EE framing reinforced**: CNN first layer as a bank of learned matched filters — same operation as PWM scanning (Artifact #4), generalised.

### Teaching purpose

A toy model of Enformer's first convolutional layer: a bank of 8 pre-built 8-bp CNN filters, each initialised to resemble a canonical TF motif (CTCF, AP1, TATA, GC-box, E-box, NFκB, SP1, poly-A). Student feeds in a DNA sequence; the artifact computes per-filter activation along the sequence and renders each as a horizontal track. Side-by-side with the PWM-scan output from Artifact #4, the student sees that a CNN filter bank is the same discrete-correlation operation as a PWM scan — just with 8 filters instead of 1 and with non-negative ReLU output.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Sequence input:                                             │
│   Preset: [ active promoter (TATA + GC-box) ▾               │
│             enhancer (AP1 + NFκB) / random / paste ]         │
│   Pasted (if chosen): [ textarea… ]                          │
├─────────────────────────────────────────────────────────────┤
│ Filter bank (8 filters, editable weights):                  │
│   [ CTCF ] [ AP1 ] [ TATA ] [ GC-box ]                       │
│   [ E-box ] [ NFkB ] [ SP1 ] [ poly-A ]                      │
│   click a filter to see its weight matrix as a mini logo     │
├─────────────────────────────────────────────────────────────┤
│ Output (SVG):                                               │
│   Top: sequence rendered as coloured per-base letters        │
│   Middle: 8 activation tracks (one per filter), stacked      │
│           post-ReLU; height = filter response at position    │
│   Annotations: strongest-firing filter at each pos           │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Preset: active promoter                              │   │
│ │ Top-firing filters: TATA (pos 120), GC-box (pos 165) │   │
│ │ Filter bank acts as a parallel PWM-scan battery      │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Sequence preset dropdown (active promoter / enhancer / silenced region / random / custom paste).
- Pasted-sequence textarea.
- Filter click to inspect a given filter's weight matrix (rendered as a small sequence logo).
- Threshold slider for the activation display: 0–5 (default 0; ReLU output).

### What they see

- Coloured per-base sequence band.
- 8 stacked activation tracks below, each labelled with the filter name, each showing the ReLU-output of that filter at each position.
- When a filter fires strongly somewhere, its track spikes at that position.
- Annotation row at top of the output showing the strongest-firing filter at each position (i.e. argmax across the filter bank per-position).

### Target aha moment

Feed in the "active promoter" preset. The TATA and GC-box filters spike at the expected positions. Now feed the same sequence into Artifact #4 (PWM Motif Scanner) with the matching JASPAR motifs — the hits align. Student sees: a CNN filter bank is a bank of PWMs evaluated in parallel. The only difference is that CNN filters can be *learned from data* rather than hand-curated from known motifs, which is what Enformer does across 256+ filters.

### Technical notes

- Pure JS, no training needed — filter weights are hand-designed to mimic canonical motif logos.
- 8 filters × 8 bp × 4 channels (A/C/G/T) = 256 weights total. Shipped as constants.
- Activation: per position, compute inner product of filter weights with the local 4-channel one-hot sequence; apply ReLU (max(0, x)).
- Filter logos: render each filter's weight matrix as a sequence logo (same helper as Artifact #4's PWM logo).
- Feasibility gate: pasted sequence must be ACGTN only, length ≥ 8.

### Acceptance criteria

- [ ] Active-promoter preset causes TATA and GC-box filters to spike at distinct positions.
- [ ] Enhancer preset causes AP1 and NFκB filters to spike.
- [ ] Random-sequence preset produces mostly flat tracks with sporadic noise firing.
- [ ] Pasted short sequence (< 8 bp) rejects with inline message.
- [ ] Filter logos render correctly for all 8 filters.
- [ ] Opens with default preset pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 9. Cross-Artifact Consistency

- All seven artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Coverage tracks use the same visual language (filled area in `--accent`, baseline in `--fg-muted`).
- Genomic coordinates rendered in JetBrains Mono; p-values and counts likewise.
- Annotations / thresholds use consistent colours: `--warning` for threshold lines, `--error` for fails, `--success` for passes.
- Every artifact emits an **outcome banner** per convention §1.4.

## 10. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls.
- [ ] `<script src="../_shared/resize.js" defer></script>` embedded near `</body>`.
- [ ] Outcome banner or equivalent verdict line visible at the end of any user interaction.
- [ ] User-input artifacts pre-flight inputs with explicit pass/fail messaging.
