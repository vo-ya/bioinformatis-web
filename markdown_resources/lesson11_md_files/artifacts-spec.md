# Lecture 11 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 11 (Long Reads and the Pangenome).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-11.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-11/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript.
- Must render standalone.
- Embedded in the lecture via `<iframe>` loaded lazily.
- **Every artifact must include `<script src="../_shared/resize.js" defer></script>` exactly once near the end of `<body>`.** C6 smoke gate.

### 1.2 Visual design

- Tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- Graph elements: nodes as rounded rectangles with sequence/ID labels; directed edges with arrow markers; paths as coloured traces overlaying the graph.
- DNA letters: base-palette colours (A `#c4342c`, C `#1e3a8a`, G `#b45309`, T `#2d7a3e`).
- Platform colours: PacBio HiFi `--accent` cobalt; ONT `--base-t` green; Illumina `--fg-muted` grey.
- Typography: Inter for UI chrome; JetBrains Mono for DNA sequences, node IDs, numerical values, GFA fields.
- Default state is instructive: opens with a meaningful example pre-rendered.
- Controls grouped in a panel above or to the left of the visualisation.
- Animations ≤ 400 ms.

### 1.3 Interaction model

- **Sliders / toggles / dropdowns** — editable parameters with sensible ranges.
- **Step / Run / Reset** — for iterative algorithms (seed-and-extend, Viterbi fill).
- **Re-simulate** — for stochastic simulations.
- Illegal input → quiet inline message (`--fg-muted`).

### 1.4 Explicit outcome reporting (required)

Every artifact answers its own question:

- Long-read vs short-read comparator → which features each tech resolves given the current setting.
- SV resolution demo → precision / recall per tech per SV size.
- Pangenome graph viewer → which haplotype is currently traced, which variants each haplotype carries.
- GFA decoder → parsed segment / link / path counts; any syntax errors identified.
- Graph seed-and-extend walker → final aligned path as a sequence of (node, offset); score + comparison to alternate paths.
- Viterbi-on-graph path finder → optimal path + best score; comparison to linear "flattened" alignment.
- Phasing demo → phase accuracy (% of SNP pairs correctly phased); comparison short vs long reads.

### 1.5 Feasibility gate on user input (required where input is free-form)

- GFA decoder: validate file on parse; report line number of any syntax errors.
- Other artifacts use preset / slider inputs only.

### 1.6 Pedagogical constraint

Every artifact produces its named aha moment. If the student plays with the controls and doesn't land on it, the artifact has failed.

### 1.7 Out of scope

- No accounts, no telemetry, no network calls beyond declared CDN libraries (none here).
- No external data files > 100 KB.

---

## 2. Artifact #1 — Long-Read vs Short-Read Comparator

**File**: `artifacts/lecture-11/01-longread-comparator.html`
**Lecture anchor**: §1.3 Why accuracy crossing ~Q20 changes the field
**EE framing reinforced**: observation-window-vs-event-duration tradeoff, per feature type.

### Teaching purpose

Let student pick a sequencing technology and a genomic feature type; see which features are resolved by each tech and which aren't. Teaches that no single tech is universally superior — the right choice depends on the target features.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Technology: [ Illumina ▾ / PacBio HiFi / ONT R10.4 / ONT UL ]│
│ Feature:    [ SNP ▾ / small indel / 5kb SV / 100kb SV /     │
│               short tandem repeat / centromere / phasing ]   │
├─────────────────────────────────────────────────────────────┤
│ Metric cards (per selected tech on selected feature):       │
│   Read length:     150 bp / 20 kb / 20 kb / 300 kb          │
│   Accuracy:        Q30+ / Q30+ / Q22 / Q22                  │
│   Coverage required:  10× / 30× / 30× / 20×                 │
│   Expected precision / recall for this feature              │
├─────────────────────────────────────────────────────────────┤
│ Tech-vs-feature heatmap matrix:                             │
│   rows: techs, cols: features, cell fill = performance      │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Illumina × 5kb SV: precision 88%, recall 58% (inference│   │
│ │   from split + discordant reads is fragile)           │   │
│ │ HiFi × 5kb SV: precision 97%, recall 96% (direct span)│   │
│ │ Switching from Illumina to HiFi boosts recall by ~40pp │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Technology dropdown (Illumina / PacBio HiFi / ONT R10.4 / ONT ultra-long).
- Feature dropdown (SNP / small indel / 5 kb SV / 100 kb SV / tandem repeat / centromere / haplotype phasing).

### What they see

- Metric cards for the currently-selected combination: read length, accuracy, coverage required, expected precision/recall.
- A tech × feature heatmap matrix: rows = techs, columns = features, cell fill = performance (0 to 100% recall with colour gradient).
- Outcome banner describing the current combination and the comparison to alternatives.

### Target aha

Select Illumina × 5 kb SV → precision/recall both modest. Switch to PacBio HiFi × same feature → both >95%. Switch to Illumina × SNP → ~100%/100%; HiFi × SNP → ~100%/100% (same — no advantage for SNP). Student sees: no tech dominates; the right tool depends on the feature.

### Technical notes

- Pure JS.
- Hardcoded precision/recall matrix (4 techs × 7 features) with values approximating 2024 literature.
- Slider-free: all discrete choices.
- Visualisation: SVG with a small dashboard layout.

### Acceptance criteria

- [ ] Default shows Illumina × SNP with high performance.
- [ ] Switching to long-read tech on SV feature dramatically raises recall.
- [ ] Heatmap cells coloured consistently by performance.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 3. Artifact #2 — SV Resolution Demo

**File**: `artifacts/lecture-11/02-sv-resolution.html`
**Lecture anchor**: §2.1 SVs resolved directly
**EE framing reinforced**: read length as observation-window size.

### Teaching purpose

Plant a structural variant; simulate short + long read coverage; see which reads span it and which don't. Report SV-calling precision/recall as a function of SV size and sequencing tech.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ SV configuration:                                           │
│   Type:     [ insertion ▾ / deletion / tandem dup ]         │
│   Size:     [──●── 5000 bp ]                                │
│   Locus:    [──●── 50000 bp in a 100 kb window ]            │
├─────────────────────────────────────────────────────────────┤
│ Read simulation:                                            │
│   Tech:     [ both ▾ / Illumina only / HiFi only ]          │
│   Coverage: [──●── 30× ]                                    │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Visualisation panel:                                        │
│   Top:    reference with SV marked                          │
│   Middle: short-read coverage (spans SV? split? discordant?)│
│   Bottom: long-read coverage (spans / doesn't span)         │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Illumina: 3 reads span (of ~1500); 12 discordant pairs│   │
│ │ HiFi: 23 reads fully span; SV directly resolved        │   │
│ │ Callers on this data: Sniffles2 recall 95%; short-    │   │
│ │   read SV callers 63% for this size                   │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- SV type dropdown: insertion / deletion / tandem duplication (default deletion).
- SV size slider: 50 bp – 50 kb (default 5000 bp).
- Locus slider: where in the 100 kb window (default centred).
- Tech: both / Illumina only / HiFi only (default both).
- Coverage slider: 5× – 60× (default 30×).
- Re-simulate button.

### What they see

- Reference track at top with the SV marked.
- Short-read coverage track: individual 150 bp reads plotted as rectangles; reads spanning the SV highlighted; split reads rendered as two halves with a connecting line; coverage curve showing SV-induced depth changes.
- Long-read coverage track: 20 kb reads plotted; reads spanning the SV highlighted.
- Outcome banner summarising per-tech resolution.

### Target aha

With SV = 5 kb, short reads have very few (or zero) direct spans and rely on discordant/split evidence — fragile. HiFi reads span the event directly. Increasing SV size to 25 kb → even HiFi struggles (would need ultra-long ONT). Student sees that the right observation window depends on event size.

### Technical notes

- Pure JS, seeded mulberry32.
- Simulate read lengths: Illumina ~150 bp constant; HiFi log-normal mean ~18 kb.
- Randomly distribute read start positions uniformly over the reference.
- A read "spans" the SV if its coordinates fully enclose the SV ± some flank (say 100 bp).
- Precision/recall computed approximately from simulation; calibrated against the matrix in Artifact #1 at matching SV sizes.

### Acceptance criteria

- [ ] Default (5 kb deletion, 30×, both techs) → HiFi has many spans, Illumina has few.
- [ ] Increasing SV size to 25 kb reduces HiFi spans.
- [ ] Re-simulate produces distinct draws.
- [ ] Opens pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 4. Artifact #3 — Pangenome Graph Viewer

**File**: `artifacts/lecture-11/03-pangenome-graph.html`
**Lecture anchor**: §3.3 The graph genome concept
**EE framing reinforced**: graph representation as explicit variation encoding.

### Teaching purpose

A small pangenome graph with 3–5 haplotypes. Student toggles between "linear view" (each haplotype as a separate linear sequence) and "graph view" (unified branching graph). Highlight which variants each haplotype carries.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Number of haplotypes:    [──●── 4 ]                         │
│ Number of variant sites: [──●── 5 ]                         │
│ View mode: [ linear ▾ / graph / both side-by-side ]         │
│ Highlight haplotype: [ none ▾ / hap1 / hap2 / hap3 / hap4 ] │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Visualisation (SVG):                                        │
│   linear mode: 4 horizontal strings (hap1..hap4)            │
│   graph mode:  branching DAG with coloured paths            │
├─────────────────────────────────────────────────────────────┤
│ Variant table:                                              │
│   site 1: SNP A/G  · hap1,3 = A; hap2,4 = G                 │
│   site 2: indel    · hap1,2 = insert; hap3,4 = reference    │
│   …                                                         │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Graph nodes: 14 · edges: 20 · paths: 4               │   │
│ │ Unique variant content: 5 variant sites across 4 haps │   │
│ │ Highlighted haplotype 2 traverses nodes {N1, N4, N7…} │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Haplotype count: 2–5 (default 4).
- Variant-site count: 3–8 (default 5).
- View-mode dropdown.
- Highlight-haplotype dropdown.
- Re-simulate.

### What they see

- SVG renderer showing either linear strings or a graph with coloured path traces.
- Variant table listing each variant and which haplotypes carry which allele.
- Summary banner with graph statistics.

### Target aha

Linear view shows 4 sequences side-by-side — redundant everywhere except at variant sites. Graph view collapses the shared sequence and shows variation as branches — same information, much more compact. Student sees why pangenomes represent variation better than a stack of linear references.

### Technical notes

- Pure JS, seeded PRNG.
- Build a random pangenome: shared "backbone" nodes + branching nodes at variant sites. Each haplotype's path is a set of choices at each branch.
- Render linear view: 4 horizontal strings with bases coloured by base palette; variant sites highlighted.
- Render graph view: nodes as boxes, edges as arrows; paths drawn as coloured traces over the graph.
- Highlight mode: make the selected haplotype's path bold; fade others.

### Acceptance criteria

- [ ] Default opens in graph view with 4 haplotypes and 5 variants.
- [ ] Linear view shows 4 separate strings.
- [ ] Highlighting a haplotype makes its path stand out.
- [ ] Variant table stays consistent with the rendered graph.
- [ ] Re-simulate produces new random graph.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 5. Artifact #4 — GFA Format Decoder

**File**: `artifacts/lecture-11/04-gfa-decoder.html`
**Lecture anchor**: §3.4 GFA — Graphical Fragment Assembly format
**EE framing reinforced**: text-to-graph deserialisation as a parse problem.

### Teaching purpose

Paste or pick a GFA text and see the rendered graph. Teaches the GFA line types (H, S, L, P, W) and their graph correspondence.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Preset: [ 2-hap small ▾ / 3-hap with SV / custom paste ]    │
│ GFA text area (read-only for presets, editable on custom):  │
│   [textarea showing the GFA content]                         │
├─────────────────────────────────────────────────────────────┤
│ Parse status: [ 0 errors, 5 segments, 6 links, 3 paths ]     │
├─────────────────────────────────────────────────────────────┤
│ Rendered graph (SVG):                                        │
│   segments as nodes with sequence labels                     │
│   links as directed arrows                                   │
│   paths as coloured trace overlays                           │
├─────────────────────────────────────────────────────────────┤
│ Hover / click to inspect:                                    │
│   - hover segment → show sequence in detail popup           │
│   - hover path → highlight its nodes                        │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Parsed 3 paths across 5 segments                       │   │
│ │ Path 'sample_A' visits: S1+ → S3+ → S5+               │   │
│ │ Error at line 7: unknown segment ID 'S99' in link     │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Preset dropdown (2-hap small / 3-hap with SV / 4-hap diverse / custom paste).
- Custom GFA text area (read-only unless "custom paste" selected).

### What they see

- Parsed GFA rendered as a graph.
- Parse statistics.
- On hover, segment sequence tooltip; path highlight.
- Inline error message with line number if parsing fails.

### Target aha

Paste a GFA file and see the corresponding graph. Remove one `L` line or corrupt one `S` line → observe the parse error at the right line number. Student sees that GFA is a human-readable text serialisation of a graph structure.

### Technical notes

- Pure JS.
- GFA parser: line-by-line, tokenise by tab, dispatch on first field.
- Segment storage: map of ID → sequence.
- Link storage: adjacency list.
- Path storage: list of (segment_id, orientation) tuples.
- Validation: links reference existing segments; paths reference existing segments; no duplicate IDs.
- SVG renderer: force-directed-ish layout for small graphs (≤20 nodes); assign x by topological order, y by branch depth.

### Acceptance criteria

- [ ] Default preset renders correctly.
- [ ] Pasting a malformed GFA produces an inline error at the right line.
- [ ] Hover on segment shows its sequence.
- [ ] Path trace overlay changes when path is highlighted.
- [ ] Opens with default preset pre-rendered.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 6. Artifact #5 — Graph Seed-and-Extend Walker

**File**: `artifacts/lecture-11/05-graph-seed-extend.html`
**Lecture anchor**: §4.1 Seed-and-extend generalised to DAGs
**EE framing reinforced**: seed placement + extension on a graph; branch-awareness of alignment.

### Teaching purpose

Step through seed placement and path-extension on a graph reference. Reader sees how extension branches at graph branches and picks the best-scoring path.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Graph: [ preset small (8 nodes) ▾ / medium (15) / large ]   │
│ Read:  [ preset spanning 2 branches ▾ / custom paste ]      │
│ Custom read: [ ACGTACGT... ]                                │
│ Seed k-mer length: [──●── 6 ]                               │
│ [Place seeds] [Step extension] [Run to completion] [Reset]  │
├─────────────────────────────────────────────────────────────┤
│ Graph visualisation:                                        │
│   Nodes as rectangles with sequences; edges directed        │
│   Seeds highlighted in warning; extension path grown in     │
│     accent; candidate extensions drawn as dashed            │
├─────────────────────────────────────────────────────────────┤
│ Score tracker:                                              │
│   Candidate A (through node 4): score 32                    │
│   Candidate B (through node 5): score 38 ✓ winner           │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Best alignment: path [N1+, N3+, N5+, N7+]             │   │
│ │ Score 38 · read fully aligned at offset 12            │   │
│ │ Losing candidate through N4 scored 32 — matched 6 of  │   │
│ │   the last 10 bases incorrectly                       │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Graph preset dropdown (small / medium / large).
- Read preset dropdown or custom paste.
- Seed k-mer length slider (4–10, default 6).
- Step / Run / Reset buttons.

### What they see

- Graph rendered on an SVG canvas.
- Seeds marked with highlight colour as they're placed.
- Extension animation: each step shows the current extension frontier; at branches, both candidates shown with scores.
- Winning path highlighted at the end; losing candidate faded.

### Target aha

Seed falls on a pre-branch node. Extension runs forward; at the branch, both out-edges are tried in parallel; the higher-scoring branch wins. Student sees graph alignment = linear alignment + branching.

### Technical notes

- Pure JS.
- Preset graphs: 8–20 nodes, 2–4 branch points.
- Seed placement: exact k-mer hash lookup over all k-length paths in the graph.
- Extension: simple banded SW-like scoring at each step; at branches, fork into two candidates; prune candidates whose running score drops more than a threshold below the best.
- Step-through mode: at each frame, show current extensions; after completion, select winner.

### Acceptance criteria

- [ ] Default preset shows seed placement in node 3, then branching extension.
- [ ] Winner is reliably selected.
- [ ] Step mode visibly advances extension.
- [ ] Losing candidate clearly rendered as dashed / faded.
- [ ] Opens pre-computed with seeds placed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 7. Artifact #6 — Viterbi-on-Graph Path Finder

**File**: `artifacts/lecture-11/06-viterbi-graph.html`
**Lecture anchor**: §5.2 Generalising Viterbi to a DAG
**EE framing reinforced**: Viterbi DP on a DAG — the core EE framing.

### Teaching purpose

Compute the optimal Viterbi alignment path through a small graph reference. Student steps through the DP fill, sees how predecessor scores combine, and traces back the optimal path.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Graph + read: [ preset small ▾ / diverged / identical ]     │
│ Mode: [ show DP fill ▾ / show traceback only ]              │
│ Read time step: [──●── 0 ]  (0..T)                           │
│ [Step] [Run to T] [Reset]                                    │
├─────────────────────────────────────────────────────────────┤
│ Main view: trellis-like grid                                │
│   rows: graph states (DAG layout)                           │
│   cols: read time steps                                     │
│   each cell: V(s, t) score                                  │
│   optimal path so far: coloured overlay                     │
├─────────────────────────────────────────────────────────────┤
│ Linear-baseline comparison (side panel):                    │
│   Show what a naive linear Viterbi would compute if the     │
│   graph were flattened into a chain                          │
│   Compare final scores: graph 38 vs linear 29               │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Optimal graph path: N1 → N3 → N5 → N8                 │   │
│ │ Best score: 38 (naive linear: 29, miss)               │   │
│ │ Graph alignment branches at node 3 to exploit variant │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Graph preset dropdown (small / diverged / identical).
- Mode (DP fill animation / traceback only).
- Time step slider.
- Step / Run / Reset.

### What they see

- Trellis grid rendering V(s, t) scores at each iteration.
- Optimal path traced as a coloured overlay.
- Side-panel showing the linear baseline's result for the same read.
- Outcome banner comparing graph vs linear scores.

### Target aha

For a read that crosses a variant site, graph Viterbi scores higher than flattened-linear Viterbi because the graph's branch allows choosing the right allele. Student sees graph alignment is more general than linear and strictly higher (or equal) score.

### Technical notes

- Pure JS.
- Graph presets: 6–10 nodes with 1–2 branches.
- DP fill: standard Viterbi update, V(s, t) = max over predecessors s' of (V(s', t-1) + match_score(s, read[t])).
- Naive linear baseline: "flatten" the graph by picking one specific path (e.g. the reference path) and run linear Viterbi on it; show inferior score.

### Acceptance criteria

- [ ] Default preset runs Viterbi fill to completion.
- [ ] Linear baseline scores lower on read that requires the alternate path.
- [ ] Step-through shows V(s, t) updates.
- [ ] Traceback produces the correct path.
- [ ] Opens with preset pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 8. Artifact #7 — Haplotype Phasing Demo

**File**: `artifacts/lecture-11/07-phasing-demo.html`
**Lecture anchor**: §6.3 Phasing as source separation
**EE framing reinforced**: phasing as supervised vs unsupervised source separation.

### Teaching purpose

Plant a diploid region with heterozygous variants. Simulate short- and long-read fragments. Run a phasing pass. Compute phase accuracy; contrast the two technologies.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Region:                                                     │
│   # heterozygous variants: [──●── 8 ]                       │
│   Region length (kb): [──●── 30 ]                            │
├─────────────────────────────────────────────────────────────┤
│ Read simulation:                                            │
│   Tech: [ compare both ▾ / short only / long only ]         │
│   Coverage: [──●── 30× ]                                    │
│ [Re-simulate]                                               │
├─────────────────────────────────────────────────────────────┤
│ Two visual panels side-by-side:                             │
│   Left: short-read phasing                                  │
│     fragments rendered; each covers 1-2 hets                 │
│     inferred phase block(s) shown                           │
│   Right: long-read phasing                                  │
│     fragments rendered; each spans multiple hets             │
│     single phase block spanning all variants                 │
├─────────────────────────────────────────────────────────────┤
│ Phase-accuracy metrics:                                     │
│   Short reads: 4 phase blocks · 92% pairwise accuracy        │
│   Long reads: 1 phase block · 100% pairwise accuracy         │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ Long reads resolve full phase from reads directly      │   │
│ │ Short reads have multiple phase-block breaks          │   │
│ │ Adding LD prior helps short reads but introduces ~5%  │   │
│ │   phase-switch errors at block boundaries              │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Heterozygous variant count: 3–15 (default 8).
- Region length (kb): 5–50 (default 30).
- Tech mode: compare both / short only / long only.
- Coverage slider.
- Re-simulate.

### What they see

- Two parallel panels (left = short reads, right = long reads).
- Each panel: fragment layout above, inferred phase block(s) below.
- Phase accuracy metrics per tech.
- Outcome banner comparing the two.

### Target aha

Short reads break into multiple phase blocks (each ~few kb); long reads give one continuous phase block spanning all variants. Student sees the direct observational nature of long-read phasing vs the statistical nature of short-read phasing.

### Technical notes

- Pure JS, seeded PRNG.
- Plant a diploid region with 8 heterozygous SNPs and known maternal/paternal haplotypes.
- Simulate short-read fragments (150 bp) covering 1–2 hets each; long-read fragments (20 kb) spanning most or all hets.
- Short-read phasing: greedy connected-component phasing across fragments that cover ≥ 2 hets jointly (like a mini-WhatsHap).
- Long-read phasing: each fragment directly determines phase of the hets it covers.
- Phase accuracy: per-pair of adjacent hets, fraction correctly phased.

### Acceptance criteria

- [ ] Default simulation: long reads produce 1 block with high accuracy; short reads produce multiple blocks.
- [ ] Reducing coverage fragments both approaches.
- [ ] Re-simulate produces new random planted haplotypes.
- [ ] Metrics update correctly.
- [ ] Opens pre-computed.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 9. Artifact #8 — T2T Gap Closer

**File**: `artifacts/lecture-11/08-t2t-gap-closer.html`
**Lecture anchor**: §6.1 T2T-CHM13 revisited
**EE framing reinforced**: visualisation of what ultra-long-read + HiFi assembly added to the canonical reference.

### Teaching purpose

Concretely show what the T2T-CHM13 assembly added to GRCh38 per chromosome, broken out by feature type (centromere, acrocentric short arm, segmental duplication, rDNA). Reader leaves able to name which chromosomes gained the most sequence and what kind of sequence it was.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Compare: [ GRCh38 ▾ / T2T-CHM13 / HPRC pangenome (overlay) ]│
│ Chromosome: [ all ▾ / chr1 / chr9 / chr13 / chr21 / …]      │
│ Feature filter: [ all ▾ / centromere / acrocentric /        │
│                   segmental dup / rDNA / telomere ]          │
├─────────────────────────────────────────────────────────────┤
│ Left panel — chromosome ideogram(s):                        │
│   GRCh38: gaps shown as hatched regions                     │
│   T2T-CHM13: same ideogram, gaps filled with colour-coded   │
│     feature types                                            │
├─────────────────────────────────────────────────────────────┤
│ Right panel — added-sequence bar chart:                     │
│   stacked bars per chromosome by feature type (Mb)          │
│   total "+Mb" readout for the current selection              │
├─────────────────────────────────────────────────────────────┤
│ ┌─ Outcome banner ──────────────────────────────────────┐   │
│ │ All chromosomes: +200 Mb total                         │   │
│ │ chr 13 (acrocentric): +15 Mb, mostly rDNA arrays        │   │
│ │ chr 9: +12 Mb, dominated by alpha-satellite HOR         │   │
│ │ Acrocentric chroms (13/14/15/21/22) gained ~2.5×        │   │
│ │   proportionally more than metacentric chromosomes     │   │
│ └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Compare dropdown: GRCh38 / T2T-CHM13 / HPRC pangenome overlay.
- Chromosome dropdown: "all" summary view or a specific chromosome (1–22, X, Y).
- Feature filter dropdown: all / centromere / acrocentric short arm / segmental duplication / rDNA / telomere.

### What they see

- Ideogram panel (left) showing the current chromosome's structure with gaps hatched (GRCh38) or filled (T2T-CHM13). If "all" is selected, a small-multiples panel shows mini-ideograms for all 24 chromosomes.
- Added-sequence bar chart (right) showing per-chromosome breakdown by feature type.
- Outcome banner with totals and a comparative observation.

### Target aha

View chr 13 in T2T mode. The acrocentric short arm — previously blank in GRCh38 — is now filled, mostly with rDNA and satellite. Compare to chr 1 (metacentric): smaller absolute gain, concentrated in centromere. Student sees that acrocentric chromosomes contributed disproportionately to the T2T gain — consistent with them having had the largest unresolved regions.

### Technical notes

- Pure JS.
- Hardcoded data table: per chromosome, the Mb of each feature type added in T2T vs GRCh38. Values drawn from Nurk et al. 2022 and supplementary tables.
- Ideogram rendering: standard chromosome-bar style, centromere as filled circle / hatched band, gaps as cross-hatch, newly-resolved regions as solid fill colour-coded by feature type.
- Bar chart: simple SVG stacked bars, per chromosome.
- HPRC overlay mode adds per-chromosome variant-count annotations (not a structural change to the ideogram).

### Acceptance criteria

- [ ] Default "all / T2T / all" shows +200 Mb total.
- [ ] chr 13 / 14 / 15 / 21 / 22 visibly gain more than metacentric chromosomes.
- [ ] Feature filter correctly restricts the bar chart.
- [ ] Chromosome dropdown updates the ideogram view.
- [ ] Opens pre-rendered with default selections.
- [ ] HTML parses; JS passes `node --check`; 1× `_shared/resize.js`.

---

## 10. Cross-Artifact Consistency

- All eight artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- Graph visualisations use the same node/edge/path conventions across artifacts #3, #4, #5, #6.
- Platform colours consistent across #1 and #2 (PacBio cobalt, ONT green, Illumina grey).
- Every artifact emits an **outcome banner** per convention §1.4.

## 11. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on colour alone for meaning.
- [ ] No `alert()`, no console spam, no external calls.
- [ ] `<script src="../_shared/resize.js" defer></script>` embedded near `</body>`.
- [ ] Outcome banner or equivalent verdict line visible at end of any user interaction.
- [ ] User-input artifacts pre-flight inputs with explicit pass/fail messaging.
