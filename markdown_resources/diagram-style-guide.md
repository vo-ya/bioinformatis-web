# Diagram Style Guide

> **Purpose**: House style for every custom diagram produced for this bioinformatics course.
> **How to use**: reference this file when creating new diagrams, or paste into a new chat along with a per-lecture `figures-spec.md` entry.
> **Companion files**: `website-spec.md` (site design system), `lecture-style-guide.md` (lecture content), `artifacts-spec.md` (interactive artifacts).

---

## 1. Philosophy

Diagrams in this course exist for three reasons:

1. **Visual rest** — to break up long prose so the eye has somewhere to land.
2. **Visual landmarks** — to anchor the reader's position in a long scrollable page.
3. **Explanation** — to carry information prose can't (spatial structure, flow, simultaneity).

They are not decoration. Every diagram earns its place by making a specific teaching point clearer than text alone.

The style is **editorial and minimal, not textbook-cartoon.** Think:
- **Bartosz Ciechanowski's physics essays** (ciechanow.ski) — clean SVG schematics, flat color, precise lines.
- **Edward Tufte** — minimum ink, no decoration.
- **Molecular Biology of the Cell (Alberts et al.)** when its figures are at their best — biologically accurate at the appropriate level of abstraction, no unnecessary detail.

Not this:
- Overly stylized "explainer video" cartoons.
- 3D-rendered glossy molecular visualizations.
- Clip-art arrows, gradient fills, drop shadows, skeuomorphic paper textures.

---

## 2. Hybrid Strategy — SVG vs Photograph

Every figure is one of two types. Never mix the two within a single figure.

### 2.1 Custom SVG (default)

Use custom SVG for **concepts, flows, abstractions, schematics, and data relationships**:

- Pipeline / workflow diagrams
- Process sequences (transcription, PCR cycle, bridge amplification steps)
- Format layouts (FASTQ record, adapter-flanked fragment)
- Structural schematics (double helix, nanopore cross-section, ZMW cutaway)
- Genetic code tables
- Signal diagrams (nanopore squiggle, Phred/dB intuition)

SVG is editable, lightweight, scales crisply, and uses the site's design tokens directly.

### 2.2 Curated photograph / real image

Use a real image **only** when a photograph conveys something an SVG can't:

- A real Sanger electropherogram trace (not an illustration)
- A NovaSeq / MiSeq / MinION device photograph
- An SEM of a ZMW or flow cell surface
- Historical photographs (e.g., the Sanger lab, the HGP announcement)

**Limit**: no more than 1–2 photographs per lecture. They should be clearly framed as "here's what the actual thing looks like" — a break from the schematic register, not a mixing of it.

**Licensing**: only use images under permissive licenses (CC0, CC BY, public domain, NIH/NHGRI, Wikimedia Commons). Always include attribution in the figure caption when required.

---

## 3. Design Tokens

All custom SVG uses the site's CSS variables. For static SVG files, the hex values are used directly (SVG doesn't inherit CSS vars from parent HTML when loaded via `<img>`). Copy the values below exactly.

```
Surface:
  --bg:             #fcfcfa   (page background — rarely used in SVG)
  --bg-muted:       #f4f3ee   (subtle fills, container backgrounds)
  --bg-inset:       #ebeae3   (slight press-in for labels)
  --border:         #e5e3dc
  --border-strong:  #c4c1b6

Text / Lines:
  --fg:             #0a0a0a   (primary lines, body text)
  --fg-muted:       #525252   (secondary lines, annotations)
  --fg-subtle:      #8a8a85   (tertiary: axis ticks, fine labels)

Accent (used for highlighted / called-out elements):
  --accent:         #1e3a8a   (deep cobalt — primary accent)
  --accent-bright:  #2563eb   (brighter blue — for small pops only)
  --accent-bg:      #eef2ff   (pale cobalt fill — highlighted regions)

DNA bases (must be used whenever bases are shown):
  --base-a:         #c4342c   (Adenine,  muted red)
  --base-t:         #0d7377   (Thymine,  deep teal)
  --base-g:         #b45309   (Guanine,  burnt amber)
  --base-c:         #2d7a3e   (Cytosine, forest green)
  --base-u:         #7c3aed   (Uracil,   violet — RNA only)
  --base-n:         #6b7280   (Unknown,  gray)

Feedback (rare in diagrams; use only when semantically relevant):
  --success:        #15803d
  --warning:        #b45309
  --error:          #b91c1c
```

> **Critical**: the base colors are used in every artifact that renders sequences and must match exactly. Never substitute a "close enough" shade.

---

## 4. Line Style

- **Default stroke width**: `1.5` units (in viewBox space, assuming typical ~800-wide viewBox)
- **Primary / structural lines**: `2.5` units
- **Fine detail lines** (axis ticks, hatching, grid): `1` unit
- **Stroke color**: `#0a0a0a` (--fg) by default; `#525252` (--fg-muted) for secondary lines
- **Stroke linecap**: `round` for open-ended lines, `butt` for lines that terminate at a joint
- **Stroke linejoin**: `round` for organic shapes, `miter` with `stroke-miterlimit="4"` for rectilinear diagrams
- **Dashed lines**: `stroke-dasharray="4 3"` for inferred / optional paths; `stroke-dasharray="2 2"` for fine construction lines
- **No outlines on text** — text is fill-only

## 5. Fill Style

- **Default fill**: `none` — most diagrams are line-art with selective fills
- **Subtle area fills**: `#f4f3ee` (--bg-muted)
- **Highlighted region fills**: `#eef2ff` (--accent-bg)
- **Solid element fills** (circles, dots, badges): `#0a0a0a`, `#1e3a8a`, or a base color as appropriate
- **No gradients. No drop shadows. No blur.**
- **No 3D / pseudo-3D effects** — everything is flat

## 6. Typography in Diagrams

SVG embeds its own typography; use the site's font stack directly via `font-family`.

- **Labels (primary)**: `font-family="Inter, system-ui, sans-serif"` — weight 500–600, size 12–14
- **Labels (secondary / captions)**: Inter, weight 400, size 10–11, color `#525252`
- **Numerical axis labels, tick labels**: Inter, weight 400, size 10, color `#8a8a85`
- **Sequences, code, formulas**: `font-family="JetBrains Mono, SF Mono, Menlo, monospace"` — weight 500, size 12–14
- **Figure titles inside the SVG** (rare — usually handled by the surrounding HTML `<figcaption>`): Source Serif 4, weight 500, size 16–20
- **Never use italic serif inside diagrams** — reserve serif italic for the surrounding prose
- **Uppercase kicker labels** (e.g., "ZOOM", "STEP 1"): Inter, weight 600, `letter-spacing="0.08em"`, color `#1e3a8a`, size 10–11

## 7. Arrows

Arrows are one of the most-abused elements in technical diagrams. Rules:

- **Default arrow**: thin, slightly tapered triangle
- Use an SVG `<marker>` defined once per document (see proof file for the pattern):
  ```xml
  <marker id="arrow-accent" viewBox="0 0 10 10" refX="9" refY="5"
          markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M 0 0 L 10 5 L 0 10 z" fill="#1e3a8a"/>
  </marker>
  ```
- **Arrow color** matches the line it terminates: cobalt for accent lines, `#525252` for secondary
- **Never use open-arrowhead style** (V-shaped) or ornate multi-segment arrows
- **Never use curved arrows** except for genuine rotational/cyclical flows (PCR cycle, Krebs-style loops)
- **Straight or right-angled** paths by default. Diagonals only when they carry meaning.

## 8. Layout Conventions

- **ViewBox**: normalize to a coordinate space where 1 unit ≈ 1 CSS pixel at intended display size. Typical: `viewBox="0 0 860 420"` for a figure displayed at ~720–860px wide.
- **Internal padding**: at least 24–32 units from viewBox edge to content
- **Whitespace**: generous — don't cram. If it looks cramped, make the viewBox bigger.
- **Alignment**: snap to an 8-unit grid when possible. Labels align to gridlines.
- **Reading order**: left-to-right, top-to-bottom, matching reading direction
- **Grouping**: use `<g>` for logical groupings; add `role` and `aria-label` for accessibility on groups that are semantically meaningful

## 9. Accessibility

Every SVG figure must include:

- `role="img"` on the root `<svg>`
- `<title>` and `<desc>` child elements as the first children, for screen readers
- `aria-labelledby="title-id desc-id"` on the root referring to them
- `alt` text on the surrounding `<img>` tag (if the SVG is loaded by reference) that summarizes the figure

Color contrast: all text-on-background combinations must meet WCAG AA (contrast ≥ 4.5 for body text, ≥ 3 for large text). The tokens above are pre-checked for `--fg` on `--bg`, and `--accent` on `--bg` and `--bg-muted`.

**Don't rely on color alone** to convey information. If a highlight is shown in accent color, also use a different stroke weight, an annotation label, or a position change.

## 10. Inline SVG vs Referenced `<img>`

Two ways to include a diagram in a lecture page:

### Use `<img src="...svg">` (default) when:
- The diagram is static (no animation, no interaction)
- It doesn't need to respond to the lecture's CSS vars dynamically
- You want caching benefits

### Inline the SVG (`<svg>...</svg>` in the HTML) when:
- The diagram has hover / click interactions
- It uses CSS animations that need to inherit from the parent page
- Part of the diagram needs to be scripted (rare — that's what artifacts are for)

**Default to `<img>`.** Only inline if there's a specific reason.

## 11. File Organization

```
diagrams/
├── lecture-01/
│   ├── 01-dna-double-helix.svg
│   ├── 02-central-dogma-flow.svg
│   ├── 03-genetic-code-table.svg
│   ├── 04-sanger-chain-termination.svg
│   ├── 05-pcr-thermal-cycle.svg
│   ├── 06-bridge-amplification.svg
│   ├── 07-flow-cell-cluster.svg
│   ├── 08-adapter-structure.svg
│   ├── 09-zmw-cross-section.svg
│   ├── 10-nanopore-squiggle-schematic.svg
│   ├── 11-fastq-anatomy.svg
│   └── photos/
│       ├── novaseq-device.jpg
│       └── minion-device.jpg
├── lecture-02/
│   └── ...
└── _shared/
    └── diagram-theme.css     # (Not strictly needed; SVGs use inline values)
```

- Filename pattern: `NN-name-kebab.svg` with zero-padded numbering
- Numbering is per-lecture, starting at `01`
- Photos live in a `photos/` subdirectory per lecture with appropriate licensing noted in a sibling `LICENSE.md`
- SVG files are authored standalone — open one in a browser and it should render correctly with no external dependencies

## 12. Authoring Workflow

When creating a new diagram for a lecture:

1. **Check `figures-spec.md`** for the lecture to find the spec for this figure.
2. **Open the proof file** (`../fastq-anatomy-proof.html`) to copy the marker/typography/color patterns.
3. **Build the SVG** following this style guide.
4. **Save as** `diagrams/lecture-NN/NN-name.svg` with just the `<svg>...</svg>` (no wrapping HTML).
5. **Verify**: open in a browser standalone — does it render cleanly? Is text legible? Are colors from the palette only?
6. **Embed in the lecture** with an `<img>` tag + `<figcaption>`.
7. **Update the lecture's figure count** and position references if applicable.

## 13. Don'ts — Quick Reference

- ❌ No gradients, drop shadows, blurs, glows, or any decorative effect.
- ❌ No 3D or pseudo-3D (isometric is OK if genuinely informational).
- ❌ No clip-art or stock illustrations imported as SVG.
- ❌ No emoji in diagrams — they're raster, they look cartoony, they don't respect the palette.
- ❌ No gratuitous color. If a line could be black, it's black.
- ❌ No text baked into curves/paths (harder to edit, harder to accessibility-read).
- ❌ No mixed fonts beyond the three listed (Inter, JetBrains Mono, Source Serif 4).
- ❌ No pixel-perfect mimicry of copyrighted textbook figures — redraw concepts in the house style.
- ❌ No 2-pixel-thick arrows, no ornate arrows, no oversized arrowheads.
- ❌ No fill without reason — default to stroke-only line art.

## 14. Pre-Submission Checklist

Before saving a diagram as final:

- [ ] Uses only colors from the design tokens (§3)
- [ ] All text in Inter, JetBrains Mono, or Source Serif 4 only (§6)
- [ ] Default stroke weight 1.5; primary structures 2.5 (§4)
- [ ] No gradients, shadows, 3D, glow effects (§13)
- [ ] Arrows use the standard marker pattern, not default SVG markers (§7)
- [ ] `role="img"`, `<title>`, and `<desc>` present (§9)
- [ ] All base colors (if bases shown) match `--base-*` exactly (§3)
- [ ] Legible at ~720px wide display (test by resizing browser)
- [ ] Renders cleanly when opened standalone in browser (no external deps)
- [ ] Filename follows `NN-name-kebab.svg` pattern (§11)
- [ ] No pixel-perfect copying from copyrighted sources

---

## 15. Canonical Proof File

The file `../fastq-anatomy-proof.html` is the reference implementation of this style. When in doubt about arrow style, typography, kicker labels, callout boxes, or monospace integration, **copy patterns from that file directly**. It is the ground truth.
