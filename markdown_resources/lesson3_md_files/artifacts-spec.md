# Lecture 3 — Artifacts Specification

> **Scope**: Interactive artifacts for Lecture 3 (DNA Sequence Assembly).
> **How to use**: hand this file to whoever implements the artifact; each section is self-contained.
> **Companion files**: `lecture-style-guide.md`, `diagram-style-guide.md`, `website-spec.md`, `lecture-03.md`.

---

## 1. Artifact Conventions (Lecture-Wide)

These conventions apply to every artifact in this lecture. Per-artifact sections below only override them when they need to.

### 1.1 Files and layout

- Each artifact is a single self-contained HTML file in `artifacts/lecture-03/NN-name.html`.
- No build step. Vanilla HTML + CSS + JavaScript. External libraries only if justified in the per-artifact section.
- The file must render standalone when opened directly in a browser.
- Artifact is embedded in the lecture page via `<iframe>` loaded lazily.

### 1.2 Visual design

- Use the design tokens from `diagram-style-guide.md` §3 via `../_shared/artifact-theme.css`.
- DNA bases, when rendered, use the `--base-*` palette exactly.
- Typography: **Inter** for UI chrome; **JetBrains Mono** for sequences, k-mers, graph node/edge labels, contig lengths.
- Default state is instructive: the artifact opens showing a meaningful example, no user input required.
- Controls grouped in a panel above or to the left of the visualization.
- No animations longer than ~400 ms. Motion only when it carries information.

### 1.3 Interaction model

- **Input / sequences / thresholds** — editable text fields or sliders, validated against the alphabet (`ACGT` unless noted).
- **Step** — advance one algorithmic step.
- **Play / Pause** — run through steps at ~1 step/sec.
- **Reset** — return to the initial example.
- **Speed** — optional slider, 0.25×–4×. Default 1×.

Illegal input shows a quiet inline message (`--fg-muted`), not a modal.

### 1.4 Pedagogical constraint

Every artifact must produce a **specific realization** — the "target aha moment" named in its section. If the student plays with the artifact and doesn't land on that realization, the artifact has failed and should be revised.

### 1.5 Out of scope

- No logins, accounts, or persistence between sessions.
- No telemetry or analytics.
- No external data files larger than ~50 KB.

---

## 2. Artifact #1 — Coverage Simulator

**File**: `artifacts/lecture-03/01-coverage-simulator.html`
**Lecture anchor**: §2.1 Coverage and average coverage
**EE framing reinforced**: coverage as Poisson sampling; mean vs distribution.

### Teaching purpose

Turn "30× coverage" from a number into a distribution. The student should leave understanding that uniform 30× and Poisson(30) look very similar — and Poisson(5) looks terrible — and that real data is neither.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Genome length L:  [──●── 1000 bp ]                          │
│ Read length r:    [──●── 150 bp  ]                          │
│ Number of reads:  [──●── 200     ]    mean coverage: 30.0×  │
│ Non-uniformity:   [──●── 0%      ]                          │
│ [Drop reads] [Reset]                                        │
├─────────────────────────────────────────────────────────────┤
│ Genome + read pileup (scrollable horizontally):             │
│ ═══════════════════════════════════════════════════════════ │
│     ▬▬▬▬       ▬▬▬▬   ▬▬▬▬                                 │
│        ▬▬▬▬  ▬▬▬▬  ▬▬▬▬                                    │
│  ...                                                        │
├─────────────────────────────────────────────────────────────┤
│ Per-position coverage:                                      │
│  [ area plot — depth over position ]                        │
├─────────────────────────────────────────────────────────────┤
│ Coverage histogram · Poisson(30) overlay                    │
│  [ vertical bars with dashed Poisson curve ]                │
│ Gaps at 0× coverage: 12 positions (1.2% of genome)          │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- Sliders for L, r, number of reads, non-uniformity (0–50%; non-zero values bias reads away from uniform distribution using a simple sinusoidal density).
- **Drop reads**: re-sample reads with current parameters and re-render.
- **Reset**: default to L = 1000, r = 150, N = 200, uniformity = 0%.

### What they see

- Reads animated in as colored bars at their start positions.
- A filled-area coverage-depth plot directly below, updating in real time.
- A histogram of per-position coverage with the theoretical Poisson(λ) curve overlaid.
- A live readout of: observed mean coverage, standard deviation, percent of positions at 0× coverage.

### Target aha moment

Crank the slider to 20% non-uniformity. The mean stays close to 30, but the histogram develops a visible bimodal tail and the "% positions at 0×" readout climbs dramatically. The student should say out loud: "mean is a bad summary of coverage."

### Technical notes

- Use `Math.random()` for read start positions; no dependency on a stats library.
- Non-uniformity implemented by accept/reject sampling against a 1 + α·cos(2πx/L) density.
- Poisson overlay computed in closed form: p(k; λ) = e^(−λ) · λ^k / k!.

### Acceptance criteria

- [ ] Default parameters produce mean ≈ 30, SD ≈ √30 ≈ 5.5, ≤ 1% zero-coverage positions.
- [ ] Non-uniformity slider at 30%+ produces ≥ 5% zero-coverage positions without changing mean by more than 10%.
- [ ] Re-rolling with the same parameters produces different but similar distributions.
- [ ] Legible at 720 px artifact width.

---

## 3. Artifact #2 — k-mer Error Correction

**File**: `artifacts/lecture-03/02-kmer-error-correction.html`
**Lecture anchor**: §2.3 Error correction
**EE framing reinforced**: signal-noise separation in the k-mer frequency domain.

### Teaching purpose

Make the k-mer spectrum concrete. The student should see the two populations (true k-mers, error k-mers) as a clearly bimodal histogram, and watch a threshold separate them.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Reference (hidden, 500 bp): ACGTACGT...                     │
│ Coverage ×:  [──●── 30 ]   k: [──●── 7 ]                    │
│ Error rate:  [──●── 1% ]                                    │
│ [Re-generate reads]                                         │
├─────────────────────────────────────────────────────────────┤
│ k-mer frequency histogram (log-y):                          │
│   [ tall bar at freq 1 | small tail | peak around freq 30 ] │
│   threshold ┈┈┈┈┈┈┈┈                                        │
│ Threshold:  [──●── 3 ]                                      │
├─────────────────────────────────────────────────────────────┤
│ Reads panel (showing 10 reads):                             │
│  read_1:  ACGTACGT·A·CGT... (2 errors corrected)            │
│  read_2:  CGTACGTA (clean)                                  │
│  ...                                                        │
│                                                             │
│ Total corrections:  47 / 512 k-mers flagged                 │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Coverage**: 5–60. Determines how separable the spectrum is.
- **k**: 5–15.
- **Error rate**: 0–5%.
- **Threshold**: 1–20. Determines which k-mers count as errors.
- **Re-generate reads**: resample the reads with a hidden reference and the current parameters.

### What they see

- A bimodal histogram with the true-kmer peak moving as coverage changes, and the error peak staying at frequency 1.
- The threshold line moving as the slider changes; the number of flagged k-mers updating.
- A scrollable list of reads; errors are shown inline with the corrected base highlighted.
- A running count of corrections made.

### Target aha moment

Drop coverage to 5×. The two peaks in the spectrum merge together and threshold selection becomes impossible without misclassifying true k-mers. The student should say: "low coverage breaks error correction because you can't distinguish true k-mers from noise."

### Technical notes

- Hidden reference is a pre-generated 500-bp sequence. Reads are sampled uniformly at random along it; each base is corrupted independently with probability equal to the error rate.
- K-mer frequencies computed directly; no approximation.
- Correction is a single-substitution Hamming-neighbor search: for each error k-mer, check all 3k neighbors at Hamming distance 1, pick any that's above threshold.
- Log-y axis: use `Math.log10` and fixed y-axis range [10⁰, 10⁴].

### Acceptance criteria

- [ ] Default (30× / k = 7 / 1% error) produces a clearly bimodal histogram.
- [ ] Dropping coverage to 5× produces visibly overlapping peaks.
- [ ] Increasing threshold above the true peak correctly flags many true k-mers as "errors" — the artifact should let this happen so the student sees the failure mode.
- [ ] Correction list shows at least 10 corrections at default parameters.

---

## 4. Artifact #3 — De Bruijn Graph Builder

**File**: `artifacts/lecture-03/03-debruijn-builder.html`
**Lecture anchor**: §3.2 De Bruijn graphs
**EE framing reinforced**: graph as shift-register state-transition diagram.

### Teaching purpose

Build a de Bruijn graph from reads interactively. The student should see the graph emerge node-by-node and edge-by-edge, and understand that each edge is one k-mer while each node is one (k−1)-mer.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Reads:                                                      │
│   [ ACGTCC ]  [ CGTCCA ]  [ GTCCAT ]     [+ Add read]       │
│ k:  [──●── 4 ]                                              │
│ [Rebuild] [Reset]                                           │
├─────────────────────────────────────────────────────────────┤
│ k-mer decomposition:                                        │
│   ACGT → ACG ─ACGT─▶ CGT                                    │
│   CGTC → CGT ─CGTC─▶ GTC                                    │
│   GTCC → GTC ─GTCC─▶ TCC                                    │
│   ...                                                       │
├─────────────────────────────────────────────────────────────┤
│ De Bruijn graph:                                            │
│   [ ACG ] → [ CGT ] → [ GTC ] → [ TCC ] → [ CCA ] → [ CAT ] │
│                                                             │
│   (hover a node to see incoming/outgoing edges with reads)  │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Reads input**: up to 8 reads, each 4–15 characters of `ACGT`.
- **k slider**: 2–6 (node labels become illegible above k = 7 in this UI).
- **Add read**: append another read input box.
- **Rebuild**: rebuild the graph from current reads and k.
- **Reset**: return to default 3 reads and k = 4.

### What they see

- The k-mer decomposition panel: for each read, each k-mer split into prefix (k−1)-mer and suffix (k−1)-mer.
- The graph rendered below using a simple force-directed or linear layout (since default reads produce a linear graph, linear is fine for small k; use a force-directed layout for branching cases).
- Hovering a node highlights all edges touching it and, in a small tooltip, lists the (read, offset) pairs that produced each edge.

### Target aha moment

Add a read that introduces branching — e.g. `CGTAC` in addition to the defaults. The graph now has a `CGT` node with two out-edges (one via CGTC, one via CGTA). The student should connect this branch to the "fork in the graph" concept that gets cleaned up in §4.

### Technical notes

- Graph rendering can be plain SVG with manual positioning for the default small cases; force-directed via `d3-force` if that's too limiting. Keep it dependency-free if possible.
- Node labels in JetBrains Mono 12.
- Edge labels in JetBrains Mono 10 `--fg-muted`.

### Acceptance criteria

- [ ] Default (3 reads, k = 4) produces a linear 6-node graph matching Figure 6.
- [ ] Adding a branching read produces a visibly branching graph.
- [ ] Hover on any node displays incoming/outgoing edges with source reads.
- [ ] Changing k to 2 produces a very small graph with many self-loops, to illustrate the degenerate case.

---

## 5. Artifact #4 — Eulerian Path Finder

**File**: `artifacts/lecture-03/04-eulerian-path.html`
**Lecture anchor**: §3.3 Solving de Bruijn graphs: Eulerian paths
**EE framing reinforced**: Hierholzer's linear-time algorithm; assembly as edge-walking.

### Teaching purpose

Step through Hierholzer's algorithm on a pre-built de Bruijn graph. The student should see the algorithm start with a partial cycle, get stuck, splice in a subpath, and terminate with every edge visited exactly once.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Pre-loaded de Bruijn graph · 8 nodes · 12 edges             │
│                                                             │
│ [Step] [Play] [Reset]          unvisited edges: 7 / 12      │
├─────────────────────────────────────────────────────────────┤
│ Graph (visited edges in cobalt, unvisited in gray):         │
│      [ ACG ] ──▶ [ CGT ] ══▶ [ GTA ]                        │
│         ║            ║         ║                            │
│      [ GAC ] ══▶ [ CTA ] ──▶ [ TAC ]                        │
├─────────────────────────────────────────────────────────────┤
│ Current path:  ACG → CGT → GTA → ...                        │
│ Assembled contig: "ACGTA..."                                │
│                                                             │
│ Algorithm trace:                                            │
│   Step 4: at GTA, 2 unvisited edges, taking GTA→TAC         │
│   Step 5: at TAC, 1 unvisited edge, taking TAC→...          │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Step**: advance the algorithm by one edge.
- **Play**: run at ~1 edge/second.
- **Reset**: return to the initial state with no edges visited.
- **Scenario**: dropdown with 3–4 preset graphs (a linear chain, a graph with a single cycle, a graph with multiple cycles requiring splicing, a disconnected graph that fails).

### What they see

- A small pre-built de Bruijn graph with nodes as circles and edges as arrows.
- As the algorithm runs, visited edges turn cobalt; the current position is highlighted; the current path is shown linearized at the bottom with the assembled contig spelled out.
- The trace panel logs each decision: "at node X, chose edge Y because Z unvisited out-edges remain."
- On completion, all edges are cobalt and the assembled contig is shown in full.

### Target aha moment

Use the "requires splicing" preset. The naive walk gets stuck mid-path; the algorithm identifies an intermediate node with unvisited edges and splices a new subpath in. The student should understand that Hierholzer's isn't greedy — it's "greedy plus splicing."

### Technical notes

- Graphs stored as adjacency lists with edge-consumption tracking.
- The "disconnected graph" preset is included specifically to show the failure case: the algorithm halts with unvisited edges and reports that no Eulerian path exists.

### Acceptance criteria

- [ ] All presets complete correctly (or report failure correctly).
- [ ] The "splicing" preset visibly splices a new subpath rather than taking a linear walk.
- [ ] The assembled contig matches the expected genome for each preset.
- [ ] Trace logs each decision with the count of unvisited out-edges.

---

## 6. Artifact #5 — Assembly Topology Inspector

**File**: `artifacts/lecture-03/05-topology-inspector.html`
**Lecture anchor**: §4.2 Graph topology
**EE framing reinforced**: topology patterns as signatures of specific error types.

### Teaching purpose

Show the three canonical topology patterns — tips, bubbles, tangles — and let the student apply cleanup operations. Each cleanup should simplify the graph toward a linear contig.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Preset: [Clean assembly ▾]                                  │
│   • Clean assembly                                          │
│   • One tip (sequencing error)                              │
│   • One bubble (SNP)                                        │
│   • One tangle (repeat)                                     │
│   • Everything (all three on one graph)                     │
│ [Apply: Remove tips] [Apply: Collapse bubbles]              │
│ [Apply: Flag tangles] [Reset]                               │
├─────────────────────────────────────────────────────────────┤
│ Graph (nodes + edges + topology highlights):                │
│   [ main linear path with a tip drawn in warning amber ]    │
│   [ a bubble highlighted in accent-bg ]                     │
│   [ a tangle highlighted in error pink ]                    │
├─────────────────────────────────────────────────────────────┤
│ Contig output:                                              │
│   CONTIG 1: ACGTACGTACGT... (length 128 bp)                 │
│   CONTIG 2: CCATGCAT...      (length 47 bp)                 │
│   Number of contigs: 2                                      │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Preset selector**: 5 canned assembly graphs, each showing 0 to 3 of the topology patterns.
- **Remove tips**: runs the cleanup operation; tips shorter than 2·k disappear.
- **Collapse bubbles**: identifies bubble sub-paths and keeps the one with higher total coverage.
- **Flag tangles**: marks tangle regions as unresolvable; contigs break at tangle boundaries.
- **Reset**: return to the loaded preset with no cleanups applied.

### What they see

- The graph with topology patterns highlighted (tips amber, bubbles accent-bg, tangles error-pink).
- After each cleanup, the graph simplifies and the contig list updates.
- The contig list shows the maximal non-branching paths remaining after cleanup.

### Target aha moment

Load the "Everything" preset. Apply cleanups in the wrong order (flag tangles first), observe that tip/bubble cleanup removes some amount afterward. Then reload and apply in the canonical order (tips → bubbles → tangles), observe that the final result is the same but with a cleaner intermediate state. The student should realize the order of cleanups is somewhat forgiving, and the end state is what matters.

### Technical notes

- Graphs hand-designed per preset, not generated.
- Cleanup operations implemented as pattern-matching on graph topology (detect dead-end branches of length ≤ 2k; detect two-path parallel segments; detect in/out-degree mismatch subgraphs).

### Acceptance criteria

- [ ] Each preset shows the expected topology patterns on load.
- [ ] Each cleanup operation simplifies the graph visibly when applicable.
- [ ] The final contig list is correct for each preset after all appropriate cleanups.
- [ ] Running the same cleanup twice is a no-op (idempotent).

---

## 7. Artifact #6 — Assembly Metrics Calculator

**File**: `artifacts/lecture-03/06-n50-calculator.html`
**Lecture anchor**: §5.2 Assembly metrics
**EE framing reinforced**: N50 as length-weighted median; Nx curve as cumulative distribution.

### Teaching purpose

Make N50, NG50, and the Nx curve concrete with a small worked example. The student should leave able to compute N50 by hand and understand that two assemblies with the same N50 can have very different length distributions.

### UI layout

```
┌─────────────────────────────────────────────────────────────┐
│ Contig lengths (one per line or comma-separated):           │
│   4500, 3200, 2800, 2100, 1500, 900, 700, 600, 400, 300     │
│ Expected genome size:  [ 18000 bp ]                         │
│ [Compute]  [Load preset ▾]                                  │
│   • Small bacterial (good assembly)                         │
│   • Small bacterial (fragmented)                            │
│   • Two assemblies to compare                               │
├─────────────────────────────────────────────────────────────┤
│ Results:                                                    │
│   Total length:       17,000 bp                             │
│   Number of contigs:  10                                    │
│   Largest contig:     4,500 bp                              │
│   N50:  2,800 bp       NG50: 2,100 bp                       │
│   N90:  400 bp         NG90: 300 bp                         │
├─────────────────────────────────────────────────────────────┤
│ Nx curve:                                                   │
│   [ horizontal axis 0-100%; vertical axis contig length ]   │
│   Cursor at 50% highlights N50 value                        │
├─────────────────────────────────────────────────────────────┤
│ (Compare mode:) side-by-side Nx curves for two assemblies   │
└─────────────────────────────────────────────────────────────┘
```

### Student controls

- **Contig lengths input**: free-text area for a list of integers.
- **Expected genome size**: single input, used for NGx computations.
- **Compute**: recomputes metrics.
- **Load preset**: three scenarios with realistic-looking contig distributions.
- **Compare**: toggles a second input set and plots both Nx curves side by side.

### What they see

- All the summary metrics updated live.
- An Nx curve showing contig length at every x from 0 to 100 (step-function curve, since Nx only changes at each contig boundary).
- Crosshair at x = 50% highlighting the N50 value.
- In compare mode, two step curves overlaid in different colours.

### Target aha moment

In compare mode, load two assemblies with identical N50 but very different shapes — one with a single large contig plus many small ones, and one with several medium-sized contigs. The student should see that N50 alone does not distinguish them, and should internalize "always report N50 alongside total length, contig count, and largest contig."

### Technical notes

- Pure JS. No chart libraries needed (simple SVG for the curve works fine).
- Input parsing tolerant of whitespace, commas, newlines.

### Acceptance criteria

- [ ] Default preset produces N50 = 2,800 and NG50 = 2,100 exactly.
- [ ] The compare preset has two assemblies with the same N50 but visibly different shapes.
- [ ] Invalid input (non-numeric, negative) shows a quiet inline error.
- [ ] Nx curve is a step function (not an interpolation).

---

## 8. Cross-Artifact Consistency

- All six artifacts share the same base CSS via `../_shared/artifact-theme.css`.
- The De Bruijn Graph Builder (#3) and the Eulerian Path Finder (#4) use the same graph rendering approach (small SVG, nodes as circles, edges as arrows with JetBrains Mono labels).
- The k-mer Error Correction artifact (#2) and the Coverage Simulator (#1) use consistent histogram styling.

## 9. Testing Checklist (Per Artifact)

- [ ] Opens standalone in the browser, no server, no console errors.
- [ ] Default state demonstrates the teaching point without interaction.
- [ ] All listed controls function.
- [ ] Listed acceptance criteria pass.
- [ ] Legible at 720 px width; degrades gracefully at 1200 px.
- [ ] No reliance on color alone for meaning.
- [ ] No `alert()`, no console spam, no external calls.
