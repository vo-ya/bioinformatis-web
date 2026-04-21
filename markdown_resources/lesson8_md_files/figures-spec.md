# Lecture 8 — Figures Specification

> **Scope**: Static diagrams for Lecture 8 (Advanced Single-Cell: Trajectories, Integration, Multi-Modal).
> **How to use**: hand each figure spec to whoever is drawing the SVG; follow the parent `diagram-style-guide.md` for all visual defaults.
> **Companion files**: `diagram-style-guide.md`, `lecture-style-guide.md`, `artifacts-spec.md`, `lecture-08.md`.

---

## 0. Conventions for This Lecture

- All figures are custom SVG. Content is algorithm-/architecture-heavy; no photographs.
- Filenames use `NN-name-kebab.svg` with zero-padded numbering.
- Each figure must be legible at 720 px and scale cleanly up to 1200 px.
- Cell-type colours should stay consistent with Lecture 7 conventions (cobalt / red / amber / violet / green / teal / pink).
- Monospace (JetBrains Mono) for: gene symbols, equations, rate constants, numeric values. Inter for everything else.
- Source Serif 4 for inset equations where the lecture renders math (only if equations are baked into the SVG rather than layered KaTeX).
- Arrows follow the shared `<marker id="arrow-accent">` and `<marker id="arrow-muted">` pattern.
- Escape `&`, `<`, `>` as XML entities in text content (`&amp;`, `&lt;`, `&gt;`).

## Figure Budget

Ten figures for a ~3h 30min lecture. Placement by part:

| # | Title | Part | Type |
|---|---|---|---|
| 1 | Pseudotime — from snapshot to ordering | Part 1 | Custom SVG |
| 2 | Three trajectory methods side-by-side | Part 1 | Custom SVG |
| 3 | Spliced and unspliced reads | Part 2 | Custom SVG |
| 4 | RNA velocity dynamical model | Part 2 | Custom SVG |
| 5 | Batch effect before and after integration | Part 3 | Custom SVG |
| 6 | scVI VAE architecture | Part 3 | Custom SVG |
| 7 | CITE-seq schematic and joint matrices | Part 4 | Custom SVG |
| 8 | Spatial platforms — Visium vs MERFISH vs Xenium | Part 5 | Custom SVG |
| 9 | Visium spot deconvolution | Part 5 | Custom SVG |
| 10 | Ligand-receptor communication network | Part 6 | Custom SVG |

---

## Figure 1 — Pseudotime, from snapshot to ordering

**File**: `diagrams/lecture-08/01-pseudotime-intuition.svg`
**Lecture anchor**: §1.1 Why pseudotime from a static snapshot
**ViewBox**: `0 0 960 460`

### Purpose

Visual core of the pseudotime intuition. Show that a static snapshot of an asynchronous developmental population implicitly contains the whole trajectory. Students should leave able to explain pseudotime without reading further prose.

### Content

Two stacked panels:

**Top panel — Snapshot of cells on a UMAP.** ~200 cells laid out along a curving trajectory (an S-curve or gentle arc). Cells coloured by a single hue — emphasising that the trajectory is *not* visible in the raw sampling. A few cell icons labelled "naïve", "intermediate", "terminal" at three spots along the curve to anchor the biology.

**Bottom panel — Same cells, coloured by pseudotime.** Same layout; cells now coloured by a 0 → 1 gradient along the curve (e.g. pale `--accent-bg` at start, deep `--accent` at end). A smooth principal curve overlaid in `--fg` 2px. Axis label removed from UMAP (coordinates are arbitrary; only ordering matters).

Annotations: the word "one snapshot" points at the top, "continuous ordering" points at the bottom, "principal curve" labels the overlaid spline.

Caption: "Pseudotime turns a single asynchronous snapshot into a trajectory position per cell."

### Style notes

- Cells: 4 px filled circles.
- Gradient: interpolate `--accent-bg` (`#eef2ff`) → `--accent` (`#1e3a8a`) in 10 discrete steps.
- Principal curve: `--fg` 2px solid with small tick marks perpendicular at four positions.
- Labels: Inter 10 `--fg`.

---

## Figure 2 — Three trajectory methods side-by-side

**File**: `diagrams/lecture-08/02-trajectory-methods.svg`
**Lecture anchor**: §1.2 Monocle, Slingshot, and PAGA
**ViewBox**: `0 0 1080 440`

### Purpose

Comparison chart. Each panel shows the same branching 3-lineage dataset analyzed by a different trajectory method. Students should leave able to pick a method based on their data's topology.

### Content

Three panels laid out horizontally, each labeled at the top:

**(a) Monocle3 — principal graph.** Cells as small dots in a UMAP-like scatter; overlaid graph skeleton (nodes at skeleton junctions, edges along the backbone). The skeleton has a central trunk and two branches. Each cell is connected to its nearest graph node with a light grey line.

**(b) Slingshot — principal curves.** Same UMAP. Instead of a graph, two smooth principal curves are drawn: one curve from the root cluster through the trunk into branch-A, another from root through branch-B. Cluster centres marked with coloured discs.

**(c) PAGA — cluster abstraction.** Same UMAP but the cells faded out; a graph at the cluster level: nodes are coloured discs sized by cluster size, edges are weighted by fraction of kNN edges between clusters. A cleaner topology summary, not per-cell pseudotime.

Bottom annotation row: one line per panel stating the topology assumption ("tree / DAG"; "tree with smooth curves"; "cluster-level graph, topology-agnostic").

### Style notes

- Each panel boxed in `--border` 1px with `--bg-muted` subtle fill.
- Cells: 3 px dots, faded (`opacity 0.4`) in panels where the graph/curves are the focus.
- Graph edges: `--accent` 1.5px.
- Principal curves: `--accent` 2.5px solid with tangent-tick marks.
- Cluster centres in panel c: filled discs in 3 consistent accent-family hues matching the Lecture 7 palette.

---

## Figure 3 — Spliced and unspliced reads

**File**: `diagrams/lecture-08/03-spliced-unspliced.svg`
**Lecture anchor**: §2.1 Spliced vs unspliced reads
**ViewBox**: `0 0 960 440`

### Purpose

Explain the basic observation that drives RNA velocity: intron-containing reads mark newly-transcribed (pre-mRNA) molecules; intron-free reads mark mature (post-splicing) mRNA.

### Content

Three stacked bands:

**Top — Gene model.** A horizontal gene with 4 exons (thick coloured bars) separated by 3 introns (thin lines). Exons labelled E1–E4; introns labeled I1–I3.

**Middle — Three transcript states.**

1. **Unspliced (pre-mRNA)**: full gene shape preserved — all 4 exons plus 3 introns. Tag: "just transcribed, introns still present".
2. **Partially spliced**: 2 exons joined, 1 intron remains. Tag: "in-progress".
3. **Fully spliced (mature mRNA)**: 4 exons joined together, no introns. Tag: "mature, ready for translation".

**Bottom — Read pile-up.** A schematic pile-up of scRNA-seq reads below the gene: most reads (say 80%) fall entirely within exons (green = spliced); some reads (say 20%) span intron-exon boundaries (amber = unspliced). Arrows from each transcript-state to its matching read class.

### Style notes

- Exons: `--base-g` (green family) filled blocks.
- Introns: thin `--fg-muted` lines.
- Spliced reads: `--base-g` mini-bars aligned to exons.
- Unspliced reads: `--warning` mini-bars spanning intron boundaries.
- All labels in Inter 10 `--fg-muted`.

---

## Figure 4 — RNA velocity dynamical model

**File**: `diagrams/lecture-08/04-velocity-model.svg`
**Lecture anchor**: §2.2 The dynamical model and state-space estimator
**ViewBox**: `0 0 1080 440`

### Purpose

The central figure tying the dynamical ODE, the u–s phase diagram, and the UMAP-space velocity arrows together. Students should leave able to re-derive the steady-state relation and explain the arrow-on-UMAP plot.

### Content

Three panels:

**(a) Kinetic diagram.** A small flowchart: α (transcription rate) → u (unspliced) → s (spliced) → ∅ (degradation). Arrow labels: α, β (splicing rate), γ (degradation rate). Below, the two ODEs written in Source Serif italics or SVG text with math glyphs:

- du/dt = α − β·u
- ds/dt = β·u − γ·s

**(b) u–s phase diagram.** Axes: x = spliced count s, y = unspliced count u. A dashed line u = (γ/β) · s represents the steady-state locus. Cells plotted as points on this plane: some cluster near the line (steady state), some are above (inducing: u/s above line), some below (repressing: u/s below line). Colour-code steady/inducing/repressing.

**(c) Velocity-on-UMAP.** A small UMAP of ~100 cells; for each cell, a short arrow in `--accent` pointing in the projected velocity direction. Arrows form a smooth flow field from "immature" cluster to "mature" cluster.

### Style notes

- Kinetic diagram: boxes `--bg-muted` with `--border-strong`; arrows `--accent` 1.5 with standard marker.
- Steady-state line: `--fg-muted` 1px dashed.
- Cells in panel b: `--accent` (steady), `--warning` (inducing), `--base-u` (repressing).
- UMAP arrows: 6 px long, `--accent` with `arrow-accent` marker.

---

## Figure 5 — Batch effect before and after integration

**File**: `diagrams/lecture-08/05-batch-integration.svg`
**Lecture anchor**: §3.1 The batch effect problem
**ViewBox**: `0 0 1000 460`

### Purpose

The canonical "why batch integration matters" visualisation. Students should recognise the shape of the failure and the shape of the fix.

### Content

Two side-by-side UMAPs of the same underlying data:

**(a) Before integration — cells coloured by batch.** Two clear lobes visible: batch-1 cells form 3 clusters on one side; batch-2 cells form 3 parallel clusters on the other side. Each cluster pair should obviously represent the same biology (same cell type from two batches) but they don't mix. Colour by batch label (two distinct accent-family hues).

**(b) After integration — same cells, coloured by cell type.** Three mixed clusters; each cluster contains cells from both batches, intermingled. Colour by cell-type label (three distinct cell-type hues). A small inset legend shows "cell type" and "batch" — the batch colours are now sprinkled uniformly within each cell-type cluster.

Arrow between the two panels labeled "scVI / Harmony / any integration method".

Annotations: "clusters by batch (fail)" under (a); "clusters by biology (success)" under (b).

### Style notes

- Panel backgrounds: `--bg` plain.
- Panels boxed in `--border` 1px.
- Cell dots: 3.5 px.
- Arrow between panels: `--accent` 2 px with marker; label above in `--fg-muted`.

---

## Figure 6 — scVI VAE architecture

**File**: `diagrams/lecture-08/06-scvi-architecture.svg`
**Lecture anchor**: §3.3 scVI — variational autoencoders on counts
**ViewBox**: `0 0 1080 480`

### Purpose

Show scVI's encoder → latent → decoder flow with explicit batch conditioning and NB likelihood. Students should leave able to relate each block to the ELBO equation.

### Content

Left-to-right flow with five blocks:

1. **Input** — a small count-vector visualisation (a row with integer values) labelled "x_i (counts)". Below it, a small one-hot bar labelled "s_i (batch)".

2. **Encoder** — an MLP block (three grey rectangles stacked, arrows between) outputting μ_z and σ_z. Labelled "q_φ(z | x, s)".

3. **Latent z** — a small multivariate gaussian cloud; arrow label "reparameterise". Prior bubble off to the side: N(0, I) with a KL-divergence arrow pointing back to the latent.

4. **Decoder** — an MLP block similar to the encoder but in reverse, taking z and s as inputs, producing ρ (expression mean profile). Labelled "f_θ(z, s)".

5. **NB likelihood** — the final block draws from NB(ρ · l, φ) where l is the per-cell library-size latent; label emphasises "callback to Lecture 6: negative-binomial counts".

Below the flow, the ELBO equation written out in one line:

ELBO = E_q [log p(x | z, s)] − KL(q || p)

with coloured brackets tying each term to its corresponding block in the diagram.

### Style notes

- Blocks: rounded rectangles in `--bg-muted` with `--border-strong` 1.5 outline.
- Arrows in `--accent` 1.5 with standard marker.
- Latent cloud: a soft `--accent-bg` blob with a few dots scattered inside.
- NB likelihood block highlighted with `--warning` outline to flag the non-standard choice (vs Gaussian in vanilla VAE).
- ELBO annotation: colour-match the reconstruction term and KL term to their blocks.

---

## Figure 7 — CITE-seq schematic and joint matrices

**File**: `diagrams/lecture-08/07-cite-seq.svg`
**Lecture anchor**: §4.1 CITE-seq — RNA plus surface protein
**ViewBox**: `0 0 960 460`

### Purpose

Show the CITE-seq chemistry trick (antibodies with DNA tags captured alongside mRNA) and the resulting two-matrix data structure. Students should be able to describe both the chemistry and the data format.

### Content

**Top band — Cell with antibody tags.** A cell cross-section with:
- mRNAs inside the cytoplasm (squiggly lines labelled "mRNA").
- Antibodies bound to surface proteins (Y-shaped icons on the membrane), each labelled with an ADT barcode tag.
- Arrow: "capture in 10x droplet" pointing to the right.

**Bottom band — Two matrices.** Side-by-side block matrices sharing the same rows (cells):

- Left: **RNA matrix** — cells × genes (labelled "20,000 genes"). Visualised as a small heatmap with ~20 visible columns and ~15 rows.
- Right: **Protein matrix** — cells × antibodies (labelled "~100 proteins"). Same 15 rows; ~10 visible columns. Different colour palette than RNA.

A brace above both matrices annotated: "same cells, two modalities".

### Style notes

- Cell cross-section: `--bg-muted` fill with `--border-strong` outline.
- mRNAs: wavy lines in `--base-a`.
- Antibodies: small Y shapes in `--accent`.
- DNA barcode tags on antibodies: small coloured boxes in 3–4 distinct accent hues.
- RNA matrix colour ramp: `--bg` → `--accent`.
- Protein matrix colour ramp: `--bg` → `--base-u` (violet) — different hue emphasises the different modality.

---

## Figure 8 — Spatial platforms, Visium vs MERFISH vs Xenium

**File**: `diagrams/lecture-08/08-spatial-platforms.svg`
**Lecture anchor**: §5.1 Platforms and resolution
**ViewBox**: `0 0 1080 480`

### Purpose

Comparison grid showing the three major spatial transcriptomics platforms side-by-side on resolution, panel size, and throughput. Students should be able to pick a platform given a study design.

### Content

Three panels horizontally at top:

**(a) Visium.** A schematic tissue square with a 5×5 grid of capture spots overlaid. Each spot circular, 55 µm labeled. Below each spot, a small "mix of 5-20 cells" annotation. Stand-out spot enlarged on the side showing multiple cell types underneath.

**(b) MERFISH.** The same tissue, but with individual cells resolved as dots; each cell shown with a few fluorescent-probe signals (small dots) superimposed. Emphasise "single-cell resolution, ~300 gene panel".

**(c) Xenium.** Similar to MERFISH but with a commercial-instrument icon. Note "~400 gene panel, automated, commercial".

Bottom table (a proper SVG `<g>` with text rows): compares three rows × three columns:

| | Resolution | Panel size | Throughput |
|---|---|---|---|
| Visium | spot (55 µm, ~10 cells) | whole transcriptome (~20k) | ~5k spots/slide |
| MERFISH | sub-cellular | 100–500 genes | ~500k cells/slide |
| Xenium | sub-cellular | ~400 genes | comparable to MERFISH |

### Style notes

- Each top panel boxed in `--border` 1px.
- Visium spots: `--accent` circles with `opacity 0.4`.
- MERFISH cells: small `--fg` dots with `--warning` probe signals on top.
- Table: Inter 10 with `--fg-subtle` row dividers.

---

## Figure 9 — Visium spot deconvolution

**File**: `diagrams/lecture-08/09-spot-deconvolution.svg`
**Lecture anchor**: §5.2 The mixed-pixel problem and deconvolution
**ViewBox**: `0 0 960 440`

### Purpose

Show the mixed-pixel problem and its deconvolution output side-by-side. Readers should understand that deconvolution produces per-cell-type *proportions* at each spot, not cell assignments.

### Content

**Left panel — One Visium spot with cells under it.** A circular spot (55 µm) shown as a dashed outline. Inside it, 5 cells of 3 different types drawn as coloured discs: 2 blue, 2 red, 1 amber. Annotation: "spot $i$ contains 5 cells: 40% T cell, 40% B cell, 20% macrophage".

Below the spot, the "observed count vector" x_i rendered as a small horizontal bar: an aggregate signal that is a weighted sum of the three cell types' profiles.

**Right panel — Deconvolution output.** The inferred proportion bar chart:

- Height of bar per cell type: T cell ~40%, B cell ~40%, macrophage ~20%, others near 0.
- Ground-truth proportions overlaid as ticks on the same chart (so reader can see inferred vs true).

Equation between panels: `x_i ≈ Σ_k w_{ik} · µ_k, w_i ≥ 0, Σ w_i = 1`

### Style notes

- Spot circle: dashed `--fg-muted` 1.5.
- Cells under spot: filled discs in 3 distinct cell-type colours matching Figure 8 palette.
- Observed-count bar: `--accent` solid.
- Proportion bar chart: filled bars in matching cell-type colours.
- Ground-truth ticks: `--fg` short horizontal marks over the bar tops.
- Equation: Source Serif italic.

---

## Figure 10 — Ligand-receptor communication network

**File**: `diagrams/lecture-08/10-ligand-receptor.svg`
**Lecture anchor**: §6.1 Ligand-receptor inference
**ViewBox**: `0 0 1000 460`

### Purpose

Canonical "cell-cell communication graph" visualisation. Students should be able to read the edge weights and ligand families.

### Content

A directed graph:

- **Nodes**: 6 cell-type circles arranged around the centre. Label each: "CD4 T", "CD8 T", "B cell", "macrophage", "fibroblast", "endothelial". Circle size proportional to cell count.
- **Edges**: directed arrows between pairs of cell types, representing ligand-receptor signalling.

Edge thickness = aggregate communication score. Edge colour = ligand family. Example axes drawn:
- T cell → macrophage (thick, `--accent`): IFN-γ family.
- Macrophage → endothelial (medium, `--warning`): VEGF family.
- Fibroblast → B cell (thin, `--base-g`): CXCL family.

Side legend: ligand-family colours with 3-4 entries (IFN family, VEGF family, CXCL/chemokine family, TGFβ family).

Under the graph: a small bar chart of "top 10 ligand-receptor pairs by score" with each pair labeled.

### Style notes

- Nodes: filled circles with `--fg` 1.5 outline, labelled in Inter 11.
- Edges: curved arrows with `arrow-accent` or `arrow-muted` markers; thickness ranges 1–4 px.
- Edge colours: 4 distinct accent-family hues for the 4 families.
- Legend on right side, Inter 10.
- Bar chart: simple horizontal bars, `--accent` gradient.
