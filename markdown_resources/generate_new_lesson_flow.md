# Generating a New Lesson — End-to-End Flow

> **Purpose.** Canonical, step-by-step recipe for bringing a new lecture into the course site. Follow top-to-bottom; nothing earlier in the flow depends on anything later.
> **Audience.** The operator driving the build (in practice: Claude Code + a human reviewer).
> **Companion files.** `website-spec.md`, `lecture-style-guide.md`, `diagram-style-guide.md`, `homepage-spec.md`, `logo-spec.md`.

---

## 0. Prerequisites

Before a new lecture is built, the following must already exist and be current:

- `markdown_resources/website-spec.md`
- `markdown_resources/lecture-style-guide.md`
- `markdown_resources/diagram-style-guide.md`
- `markdown_resources/homepage-spec.md`
- `markdown_resources/logo-spec.md`
- `assets/styles.css`, `assets/lecture.css`, `assets/logo/*.svg`
- `artifacts/_shared/artifact-theme.css`, `artifacts/_shared/resize.js`
- `proof_html_resources/fastq-anatomy-proof.html`, `proof_html_resources/logo-proof-v2.html`

If any of these changes, all subsequent lectures inherit the change on the next deploy — individual lectures never fork the shared design system.

---

## Three phases

| Phase | Who | Output |
|---|---|---|
| **A — Inputs** | Author (a separate Claude chat, or hand-written) | Three per-lecture spec files |
| **B — Build**  | Claude Code | Lecture HTML · figures · artifacts · homepage promotion |
| **C — Quality** | Claude Code + reviewer | Smoke test · review pass · commits · push |

The AI never starts Phase B until Phase A is complete.

---

## Phase A — Inputs (author-side)

The author opens a separate Claude chat, provides a **theme** (one-line lecture topic) and a **list of things to cover** (topic outline), and has that chat produce the three per-lecture spec files. Those files are then dropped into this repo at:

```
markdown_resources/lessonN_md_files/
├── lecture-NN.md           # Lecture content, following lecture-style-guide.md
├── figures-spec.md         # Per-figure build spec for every figure in the lecture
└── artifacts-spec.md       # Per-artifact behaviour + UI spec for every artifact
```

Naming conventions:
- **Folder**: `lesson{N}_md_files/` with single-digit `N` (e.g. `lesson1_md_files`, `lesson2_md_files`, … `lesson16_md_files`).
- **Lecture .md filename**: `lecture-{NN}.md` with zero-padded two-digit `NN` (e.g. `lecture-01.md`, `lecture-02.md`).
- The folder number and file number always match: `lesson3_md_files/lecture-03.md`.

### What each spec file must contain

- **`lecture-NN.md`** — follows `lecture-style-guide.md` §3 skeleton exactly: Duration/Audience/File top-matter blockquote, Learning Objectives (5–8 items, verb-first), numbered Parts with duration, hierarchically numbered subsections, `**EMBED — Artifact #k: Name**` markers, `**FIGURE — Figure #k: Name**` markers, Wrap-up block, Appendix timing table. Every callout uses one of the approved labels (Intuition box, EE framing, Historical pointer, Discussion prompt, Warning box).
- **`figures-spec.md`** — one section per figure, in the order they appear in the lecture. Each section: file path, lecture anchor (§x.y), viewBox, purpose, content description, style notes, style guide §references. Follows the format of `lesson2_md_files/figures-spec.md`.
- **`artifacts-spec.md`** — one section per artifact. Each: file path, lecture anchor, teaching purpose, UI layout (text sketch), student controls, what they see, target "aha moment", technical notes, acceptance criteria. Follows the format of `lesson2_md_files/artifacts-spec.md`.

### Phase A exit check — embed-marker consistency

Before handing off to Phase B, verify that the three files agree with each other. This is the single most common source of downstream rework.

**1. Artifact count matches.**

- Count `**EMBED — Artifact #k: …**` markers in `lecture-NN.md`. Call this `N_a`.
- Count top-level artifact sections (`## Artifact #k — …`) in `artifacts-spec.md`. Call this `M_a`.
- **Require `N_a == M_a`.**

**2. Figure count matches.**

- Count `**FIGURE — Figure #k: …**` markers in `lecture-NN.md`. Call this `N_f`.
- Count top-level figure sections (`## Figure k — …`) in `figures-spec.md`. Call this `M_f`.
- **Require `N_f == M_f`.**

**3. Numbering is contiguous from 1.**

Both artifacts and figures are numbered `1..N` with no gaps. If the lecture text says "Artifact #3" but only Artifacts #1, #2, #4 exist in the spec, fix the spec or the lecture before proceeding.

**4. Filename alignment.**

Each EMBED marker in `lecture-NN.md` names a target path like `artifacts/lecture-NN/NN-name.html`. The corresponding section in `artifacts-spec.md` must declare the same filename under its `File:` field. Same rule for figures.

**5. Section-anchor alignment.**

Each artifact and figure names a `Lecture anchor:` (e.g. `§2.2 Suffix arrays`). That subsection must actually exist in `lecture-NN.md`. If the spec says `§3.4` but the lecture has no 3.4, fix before proceeding.

**6. Callout sanity.**

Skim `lecture-NN.md` once. Every callout uses one of the six labels in `lecture-style-guide.md` §4 — no invented types, no missing colons. The soft count targets (§B1 below) are not enforced here, but clearly off-target distributions (e.g. 20 EE framings or zero Intuition boxes) are cheap to flag now and expensive to rebalance later.

**7. Top-matter exists.**

`lecture-NN.md` has the Duration/Audience/File blockquote at the very top. File value matches `lectures/lecture-NN.html`.

**If any of 1–7 fails, Phase A is not done.** Send the author back to the generator chat with the specific mismatch; do not start Phase B. A 5-minute fix here saves hours of rework in B2-B5.

---

## Phase B — Build (AI-side)

Eight steps, executed in order, one commit per Phase B step or at natural boundaries.

### B1 · Read specs and report an inventory

**Inputs:** the three new per-lecture specs + all five shared specs.
**Outputs:** a short plain-text report (no files modified). Contents:

- Number of Parts, subsections, minutes per Part, total duration.
- Number of figures, per-Part distribution, filename list.
- Number of artifacts, per-Part distribution, filename list.
- Callouts by type (target distribution below).
- Anything unclear, missing, or inconsistent between the three specs.

If the inventory flags anything blocking (missing spec, contradictions), stop and fix the specs before moving on.

**Soft target for callout balance per lecture** (from `lecture-style-guide.md` §4):

| Label | Target count |
|---|---|
| Intuition box | 3–4 |
| EE framing | 5–7 |
| Historical pointer | 1–3 |
| Discussion prompt | 1–3 |
| Warning box | 2–5 |
| Section emphasis (§4.7 bolded lead, no label) | 0–1 |

These are guidance, not gates. The review pass (step C7) checks the distribution and rebalances if needed.

### B2 · Build lecture HTML (scaffold + prose, one pass)

**Inputs:** `lecture-NN.md`, `lecture-style-guide.md`, `website-spec.md`, `homepage-spec.md` (for brand block pattern), existing `lectures/lecture-01.html` or `lecture-02.html` as a template.
**Output:** `lectures/lecture-NN.html` with:

- `<head>`: charset, viewport, title, meta description, og:title, og:description, og:type, favicon link, Google Fonts preconnect + stylesheet, site CSS links, KaTeX CSS.
- Scroll-progress bar div.
- Lecture page wrapper (desktop sticky TOC aside + mobile collapsible TOC details + main column).
- Brand link back to course home.
- Lecture header (kicker, title, subtitle, duration + audience meta).
- Italic-serif tagline (the "N moves" single-line summary — cut freshly from the lecture's themes).
- Learning objectives box.
- All Parts and subsections with prose filled in from the .md.
- Figure embeds as `<figure class="figure figure--wide"><img><figcaption></figcaption></figure>`.
- Artifact embeds as `<figure class="artifact-embed artifact-embed--wide"><figcaption class="artifact-caption">…</figcaption><iframe class="artifact-frame" loading="lazy"></iframe></figure>`.
- Wrap-up section with takeaways, next lecture, homework, references, timing appendix.
- Bottom scripts: scroll-progress updater, TOC IntersectionObserver highlight, iframe postMessage auto-resize, heading-anchor copy, KaTeX auto-render.
- Bottom of body: KaTeX CDN scripts (`katex.min.js` + auto-render).

At this step, figure and artifact paths resolve to files that do not yet exist — the iframes and images will 404 in the browser until B3 and B4 complete. That's expected.

**Callout rendering rules** when translating from markdown:

| Markdown label | HTML class | Label text |
|---|---|---|
| `> **Intuition box**: …` | `callout callout--info` | `Intuition box` |
| `> **EE framing**: …` | `callout callout--info` | `EE framing` (subtitles OK: `EE framing — dB`, etc.) |
| `> **Historical pointer**: …` | `callout callout--note` | `Historical pointer` |
| `> **Discussion prompt**: …` | `callout callout--discussion` | `Discussion prompt` |
| `> **Warning box**: …` | `callout callout--warning` | `Warning box` |
| `> **This is the deepest technical block…**` | `callout callout--note` | *(no label, bold sentence is the lead)* |

Commit: `Lecture NN scaffold + prose`.

### B3 · Build static figures

**Inputs:** `figures-spec.md`, `diagram-style-guide.md`, proof files.
**Output:** `diagrams/lecture-NN/NN-*.svg` — one per figure.

Rules (from `diagram-style-guide.md`):

- Every SVG has `role="img"`, `<title>`, `<desc>`.
- Arrow markers use the `arrow-accent` + `arrow-muted` defs pattern from the proof file.
- Colors from the locked palette only: `--fg`, `--fg-muted`, `--fg-subtle`, `--accent`, `--accent-bg`, `--bg-muted`, `--bg-inset`, `--border`, `--border-strong`. Base colors (`--base-a/t/g/c/u/n`) only when actual nucleotides are shown.
- Fonts: Inter / JetBrains Mono / Source Serif 4 only.
- Default stroke width 1.5, primary structural 2.5.
- No gradients, drop shadows, 3D, glows.
- viewBox as specified per-figure in `figures-spec.md`.

Validation: each SVG parses with `python3 -c "import xml.etree.ElementTree as ET; ET.parse('...')"`.

Commit: `Lecture NN figures`.

### B4 · Build interactive artifacts

**Inputs:** `artifacts-spec.md`, the lecture's `lecture-NN.html` (for placement context), shared `artifact-theme.css` + `resize.js`.
**Output:** `artifacts/lecture-NN/NN-*.html` — one per artifact.

Rules (from `artifacts-spec.md` §1 and `website-spec.md` §7):

- Single self-contained `.html` file per artifact.
- Imports `../_shared/artifact-theme.css`.
- Includes `<script src="../_shared/resize.js" defer></script>` near the end of body.
- Vanilla JS. No framework. No build step. CDN libraries (Chart.js, KaTeX) only if justified per-artifact.
- `<main class="artifact" data-artifact="slug">` as root.
- `<h1>` + one-sentence caption at top (so it reads standalone).
- **Default state is instructive** — artifact opens with a meaningful example, no input required.
- Controls grouped in one panel; outputs in another.
- Font-family: Inter for UI chrome, JetBrains Mono for sequences/algorithmic state.
- All tokens come from the imported theme CSS — no inline hex colours.
- No `alert()`, no `localStorage`, no analytics, no external network calls beyond declared CDN libraries.

Validation: inline JS passes `node --check`; HTML parses well-formed; artifact opens standalone in a browser with no console errors; each acceptance criterion in its `artifacts-spec.md` section passes.

Commit: `Lecture NN artifacts`.

### B5 · Promote homepage card

**Inputs:** `index.html`, running totals from already-shipped lectures.
**Output:** `index.html` with the lecture-N card flipped from placeholder to linked, and the stats strip updated.

Changes:

1. Replace the `<li><div class="lecture-card lecture-card--placeholder">…</div></li>` for lecture N with a linked `<a class="lecture-card" href="lectures/lecture-NN.html">…</a>` carrying real metadata: `≈ Hh MMmin · K figures · M interactive tools`.
2. Update `.homepage-stats`: the running totals (lectures, interactive tools, figures, contact hours).
3. `aria-label` on the linked card follows the pattern: `Lecture NN: Title. Duration approximately H hours M minutes.`

Commit: `Lecture NN homepage promotion`.

---

## Phase C — Quality

### C6 · Smoke test

Start `python3 -m http.server 8000` at the repo root. Verify:

- `http://localhost:8000/index.html` returns 200.
- `http://localhost:8000/lectures/lecture-NN.html` returns 200.
- All `artifacts/lecture-NN/*.html` return 200.
- All `diagrams/lecture-NN/*.svg` return 200.
- `lecture-NN.html` parses well-formed (no unclosed tags).
- Every inline `<script>` block in the lecture and in each artifact passes `node --check`.
- Every new SVG parses as XML.
- Opening the lecture in a browser: no console errors; iframes auto-size; scroll-progress bar tracks scroll; TOC highlights current section.

### C7 · Review pass

Open the rendered lecture and check:

1. **Callout balance** against the soft targets in §B1. If any type is far outside the range, rebalance (remove the weakest, add where a section lacks one).
2. **Figure placement** matches `figures-spec.md` anchors. Every figure has an alt text and a caption that describes what it shows, not what it's called.
3. **Artifact placement** matches `artifacts-spec.md` anchors. Every artifact's default state is visible without user interaction.
4. **Equations** render (KaTeX auto-render fired). No $ leakage.
5. **Cross-references** — every "see Part N" / "Figure K" pointer is correct.
6. **Responsive** — check 720 px (prose column), 1024 px (lecture content), and 1400 px (artifact-wide). Mobile TOC collapses to `<details>`.
7. **Accessibility** — focus-visible outlines; all interactive elements Tab-reachable; no colour-alone signals.

### C8 · Commit cadence and push

One commit per Phase B step (B2, B3, B4, B5) plus any review-pass polish commits from C7. Messages follow the repo's existing style: imperative subject (under 70 chars), body explaining the why, `Co-Authored-By:` trailer. Push once the review pass is clean.

---

## Per-step cheat sheet

| Step | Reads | Writes | Validates | Commit |
|---|---|---|---|---|
| B1 | 3 per-lecture specs, 5 shared specs | — (report only) | — | — |
| B2 | `lecture-NN.md`, style guides, existing lecture as template | `lectures/lecture-NN.html` | HTML well-formed; inline JS `node --check` | yes |
| B3 | `figures-spec.md`, `diagram-style-guide.md`, proof files | `diagrams/lecture-NN/NN-*.svg` | XML parse; title/desc present | yes |
| B4 | `artifacts-spec.md`, shared theme/resize | `artifacts/lecture-NN/NN-*.html` | `node --check`; opens standalone | yes |
| B5 | `index.html`, running totals | `index.html` | card flip complete; stats updated | yes |
| C6 | everything new | — | 200s, parse checks | — |
| C7 | rendered lecture | (possibly) edits | callout balance; placement; rendering | yes if polish |
| C8 | — | — | `git log`/`git push` | push |

---

## Common pitfalls

- **Missing Phase A cross-consistency.** The `lecture-NN.md` embed markers name N artifacts and M figures; `artifacts-spec.md` must spec exactly N, and `figures-spec.md` must spec exactly M. Mismatches discovered in B2 cost a round-trip to the author.
- **Inline hex in artifacts.** Don't. Artifacts import `artifact-theme.css` and use CSS variables. The only place bare hex is allowed is inside SVGs loaded via `<img>` (no CSS inheritance possible) — and only using the documented values from `diagram-style-guide.md` §3.
- **Forgetting KaTeX in `<head>`.** If the lecture has any `$$…$$` or `$…$` math, the lecture HTML must load both the KaTeX CSS in `<head>` and the JS + auto-render scripts near the bottom of `<body>`.
- **Broken iframe paths.** The embed contract says `src="../artifacts/lecture-NN/NN-name.html"`. Both the path prefix and the trailing filename must match exactly — a mismatch silently 404s the iframe.
- **Callout labels drift.** Use the six labels from `lecture-style-guide.md` §4 verbatim — no invented types.
- **Figure alt text identical to caption.** Write alt text that describes the visual content (what a screen reader should hear); write captions that state the teaching point.

---

## Quick reference: relative paths from inside each output file

| Output | Path to shared CSS | Path to shared JS | Path to logo |
|---|---|---|---|
| `lectures/lecture-NN.html` | `../assets/styles.css`, `../assets/lecture.css` | inline | `../assets/logo/favicon.svg` |
| `artifacts/lecture-NN/xx.html` | `../_shared/artifact-theme.css` | `../_shared/resize.js` | — |
| `diagrams/lecture-NN/xx.svg` | colours inlined as hex (no CSS) | — | — |
| `index.html` | `assets/styles.css` | inline | `assets/logo/favicon.svg` |

---

## Version history of this flow

- **v1 (2026-04-19).** Initial draft, derived from the Lecture 1 and Lecture 2 builds. Lecture 1 was built before this flow was written; Lecture 2 was the first lecture to approximately follow it. Lecture 3 onward should follow it exactly.
