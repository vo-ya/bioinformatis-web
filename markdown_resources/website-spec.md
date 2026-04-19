# Bioinformatics Course Website — Technical & Design Spec

> **Audience for this doc**: Claude Code / a developer implementing the site.
> **Audience for the site itself**: Electrical Engineering students taking a bioinformatics course. Lectures are 3–4 hours each.

---

## 1. Goals

1. Host lecture material for an EE-flavored bioinformatics course.
2. Each lecture page contains prose, figures, equations, and **embedded interactive artifacts**.
3. Artifacts must also work standalone (projected full-screen during class, opened in a new tab by students).
4. Zero build step. Static HTML/CSS/JS only. Hostable on GitHub Pages, Netlify, or any Apache/Nginx server.
5. Fast to iterate on: editing a lecture should not require touching artifacts, and vice versa.

## 2. Format Decision: Scrollable Lectures, Not Slides

**Each lecture is a single long scrollable HTML page, not a slide deck.**

Rationale:

- A 3–4 hour lecture is a coherent narrative. Slides would force arbitrary chunking of content that reads naturally as connected prose with inline diagrams and artifacts.
- The same page serves two audiences: the lecturer projecting it during class, and students reviewing it afterwards. Slides are poor review material; a scrollable page with prose reasoning is good study material.
- Interactive artifacts embed naturally into reading flow as iframes; forcing each onto its own slide breaks the rhythm.
- This format matches the visual references chosen for the design system (Blue Bottle, Papers with Code — both scrollable editorial layouts).

Consequences for implementation:

- No presentation framework (no reveal.js, no impress.js, no Slidev).
- No per-slide page breaks, no `[class="slide"]` containers.
- The lecturer teaches live by scrolling through the projected page. Artifacts can be opened standalone in a new tab for full-screen demos.
- **No speaker notes feature.** The lecture prose is written to be read by the student directly; there is no separate authoring layer for the lecturer. Keeping a single audience removes a whole class of edge cases (forgotten toggles, leaked content).
- A future `.presenter-mode` CSS class may scale up type sizes for projection. That is a styling concern, not a format concern.

## 3. Architecture

**Static site. Plain HTML, CSS, vanilla JS.** No framework, no bundler, no Node build.

- Each lecture is a single HTML file that embeds artifacts via `<iframe>`.
- Each artifact is a single self-contained HTML file with inline CSS and JS. External libraries loaded from CDN where needed (Chart.js, KaTeX, etc.).
- A shared `styles.css` handles site-wide design; artifacts import `artifact-theme.css` for design-token consistency but otherwise stand alone.

### Why iframes for artifacts

- **Isolation**: an artifact's JS and CSS can't break the lecture page.
- **Reuse**: the same artifact HTML can be embedded in multiple lectures or opened standalone without modification.
- **Portability**: if the site is later migrated to a different stack (Next.js, Hugo, etc.), the artifacts move unchanged.
- **Simplicity**: no component framework, no hydration, no build.

Trade-off: iframes require explicit height management. Each artifact declares its natural height; the embedding lecture sets iframe height accordingly. A small `postMessage`-based height auto-resize script (`resize.js`) is included so artifacts can tell their parent page how tall they need to be.

## 4. Folder Structure

Lectures and artifacts are kept as separate top-level directories. Inside `artifacts/`, files are grouped by lecture so ownership is obvious at a glance.

```
bioinformatics-course/
├── index.html                          # Course landing page
├── lectures/
│   ├── lecture-01.html
│   ├── lecture-02.html
│   └── ...
├── artifacts/
│   ├── lecture-01/
│   │   ├── 01-dna-explorer.html
│   │   ├── 02-central-dogma.html
│   │   ├── 03-cost-explorer.html
│   │   ├── 04-illumina-basecaller.html
│   │   ├── 05-nanopore-squiggle.html
│   │   └── 06-fastq-inspector.html
│   ├── lecture-02/
│   │   └── ...
│   └── _shared/
│       ├── artifact-theme.css          # Design tokens used by all artifacts
│       └── resize.js                   # postMessage height auto-resize
├── diagrams/
│   ├── lecture-01/
│   │   ├── 01-dna-double-helix.svg
│   │   ├── 02-central-dogma-flow.svg
│   │   ├── ...
│   │   └── photos/                     # Curated real photographs (if any)
│   │       └── novaseq-device.jpg
│   ├── lecture-02/
│   │   └── ...
│   └── _shared/
│       └── diagram-theme.css           # SVG design tokens (mirrors site vars)
├── assets/
│   ├── styles.css                      # Site-wide styles
│   ├── lecture.css                     # Lecture-page styles
│   ├── images/                         # Site-level images (logo, etc.)
│   └── data/                           # e.g. genome-cost.json if externalized
└── README.md
```

### Path conventions

Embedding an artifact from a lecture page (`lectures/lecture-01.html`):

```html
<iframe src="../artifacts/lecture-01/02-central-dogma.html" ...></iframe>
```

Artifact files import shared resources (from `artifacts/lecture-01/<file>.html`, go up one level into `_shared/`):

```html
<link rel="stylesheet" href="../_shared/artifact-theme.css">
<script defer src="../_shared/resize.js"></script>
```

Lecture pages link site CSS:

```html
<link rel="stylesheet" href="../assets/styles.css">
<link rel="stylesheet" href="../assets/lecture.css">
```

Embedding a diagram from a lecture page:

```html
<figure class="figure">
  <img src="../diagrams/lecture-01/02-central-dogma-flow.svg" alt="...">
  <figcaption>...</figcaption>
</figure>
```

SVG diagrams are referenced as `<img src="...svg">` for simplicity (cacheable, no DOM bloat) unless a diagram needs to be interactive, in which case it's inlined. See `diagram-style-guide.md` for when to inline vs reference.

All paths are relative — no absolute URLs — so the site works on any host subpath.

## 5. Design System — "Lecture Notes"

### 5.1 Design references & principles

The aesthetic is synthesized from two references:

- **Blue Bottle Coffee** — restrained, editorial, warm off-white backgrounds, confident cobalt-blue accent, serif headlines paired with clean sans-serif body. The "premium, minimal, confidence-through-restraint" side.
- **Papers with Code (legacy site)** — utilitarian academic, dense information calmly presented, blue links, monospace for technical content. The "no-nonsense, respects your time" side.

Principles:

- **Editorial feel.** Serif headings give the lecture pages the gravitas of a university course or scientific journal; sans body keeps it readable for long sessions.
- **One accent, used sparingly.** Deep cobalt carries the identity. Nothing competes with it.
- **Content-first.** No decoration, no marketing flourishes. Typography and whitespace do the heavy lifting.
- **Readability over density.** Generous line-height, max content width ~720px for prose.
- **Monospace for biology.** Any DNA sequence, FASTQ line, or code snippet renders in a monospace face. Non-negotiable — biological sequences are visually unusable in proportional fonts.
- **Projector-aware.** Warm off-white is readable under lecture-hall projection; pure white washes out.

### 5.2 Color Tokens

Define as CSS custom properties on `:root` in `assets/styles.css`. Artifacts import the same tokens via `artifact-theme.css`.

```css
:root {
  /* ─── Surface ───────────────────────────────────────── */
  --bg:             #fcfcfa;   /* Warm near-white, page background */
  --bg-elevated:    #ffffff;   /* Cards, artifact frames — pure white for contrast */
  --bg-muted:       #f4f3ee;   /* Subtle section tints, code blocks, sequence blocks */
  --bg-inset:       #ebeae3;   /* Slight press-in for table headers, kbd, etc. */

  --border:         #e5e3dc;
  --border-strong:  #c4c1b6;

  /* ─── Text ──────────────────────────────────────────── */
  --fg:             #0a0a0a;   /* Near-black body, editorial weight */
  --fg-muted:       #525252;   /* Captions, metadata */
  --fg-subtle:      #8a8a85;   /* De-emphasized helper text */

  /* ─── Accent ─ Blue Bottle cobalt + PwC link blue ───── */
  --accent:         #1e3a8a;   /* Deep cobalt — signature identity color */
  --accent-hover:   #1e40af;
  --accent-bright:  #2563eb;   /* Brighter blue for inline links */
  --accent-bg:      #eef2ff;   /* Pale blue tint for callouts and highlights */

  /* ─── Semantic — DNA bases (used across every artifact) ─ */
  --base-a:         #c4342c;   /* Adenine   — muted red */
  --base-t:         #0d7377;   /* Thymine   — deep teal, distinct from cobalt accent */
  --base-g:         #b45309;   /* Guanine   — burnt amber */
  --base-c:         #2d7a3e;   /* Cytosine  — forest green */
  --base-u:         #7c3aed;   /* Uracil    — violet (RNA) */
  --base-n:         #6b7280;   /* Unknown   — gray */

  /* ─── Feedback ──────────────────────────────────────── */
  --success:        #15803d;
  --warning:        #b45309;
  --error:          #b91c1c;
  --info:           var(--accent-bright);
}
```

> **Important**: once the base-color palette (`--base-a/t/g/c/u/n`) is set, do not change it. Students will build visual memory across lectures — consistency across artifacts matters more than aesthetic tweaks.

### 5.3 Typography

```css
:root {
  /* Serif for headings — editorial gravitas, Blue Bottle-adjacent */
  --font-serif: "Source Serif 4", "Charter", "Iowan Old Style", Georgia, serif;

  /* Sans for body and UI — tight, neutral, high-density reading */
  --font-sans:  "Inter", system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;

  /* Mono for sequences, FASTQ, code */
  --font-mono:  "JetBrains Mono", "IBM Plex Mono", "SF Mono", Menlo, Consolas, monospace;

  /* Scale */
  --text-xs:    0.75rem;   /* 12px */
  --text-sm:    0.875rem;  /* 14px */
  --text-base:  1rem;      /* 16px */
  --text-lg:    1.125rem;  /* 18px */
  --text-xl:    1.25rem;   /* 20px */
  --text-2xl:   1.5rem;    /* 24px */
  --text-3xl:   2rem;      /* 32px */
  --text-4xl:   2.75rem;   /* 44px — lecture title */

  /* Line-height */
  --leading-tight:   1.2;   /* Large serif headings */
  --leading-normal:  1.6;   /* Body */
  --leading-relaxed: 1.75;  /* Long-form prose */
}
```

**Font loading**: Source Serif 4, Inter, and JetBrains Mono are all available on Google Fonts. Preferred: self-host in `assets/fonts/` for reliability (university networks occasionally block CDNs).

**Usage rules**:

- `h1`, `h2`, `h3` → serif, `--leading-tight`, slightly negative letter-spacing on the largest sizes.
- `h4`, `h5`, `h6` → sans, optionally uppercase or small-caps at smaller sizes for section dividers.
- Body paragraphs, lists, captions → sans.
- Any DNA/RNA sequence, FASTQ line, inline code, code block → mono.
- Links → `--accent-bright`, underline on hover. No underline at rest for body links within prose (cleaner editorial feel); TOC/nav links can be underlined always.

### 5.4 Spacing & Layout

- **Prose max-width**: 720px.
- **Lecture page max-width**: 960px for content, with artifacts allowed to break out to wider widths via an `.artifact-embed--wide` modifier class.
- **Landing page max-width**: 1100px to let the lecture-card grid breathe.
- **Vertical rhythm**: 1.6 line-height body. Heading margins: `2em` above / `0.5em` below for `h2`, proportionally scaled for other levels.
- **Section spacing**: `4rem` between major lecture parts.
- **Grid**: CSS Grid for the lecture-card listing on the landing page; plain flow for lecture pages.

### 5.5 Components

Patterns used throughout. Styled in `assets/styles.css` / `assets/lecture.css`.

| Component | Usage & style notes |
|---|---|
| **Lecture title** | Serif, `--text-4xl`, tight leading, `--fg`. Duration badge (e.g., "≈210 min") in sans, `--text-sm`, `--fg-muted`, aligned to the right of the title. |
| **Part header (`h2`)** | Serif, `--text-3xl`. Optional small-caps kicker above it ("PART 3") in sans, `--fg-subtle`. |
| **Section header (`h3`)** | Serif, `--text-xl`. |
| **Body prose** | Sans, `--text-base`, `--leading-relaxed`, `--fg`. |
| **Callout (info)** | "Intuition box" — analogies, EE perspective framings. Cobalt left border (3px), `--accent-bg` background, subtle rounded corners. Optional small icon. |
| **Callout (note)** | Side remarks, historical context. Gray left border, `--bg-muted` background. |
| **Callout (warning)** | Common misconceptions, pitfalls. Amber left border (`--warning`), faint amber-tinted background. |
| **Callout (discussion)** | Student-facing discussion prompts — "questions to think about." Cobalt left border with a slightly different icon or label than info callouts to distinguish them. |
| **Scroll-progress bar** | 2–3px cobalt bar fixed to the top of the viewport, filling as the user scrolls. Width tracks `scrollY / (scrollHeight - innerHeight)`. |
| **Heading anchor link** | Subtle "¶" or link icon that appears on hover next to each section heading. Clicking copies the deep URL (`#section-id`) to the clipboard. |
| **Equation block** | Centered, `--bg-muted` background, generous padding. KaTeX rendering. |
| **Sequence block** | Monospace, `--bg-muted` background, colored per-base using `--base-*` tokens. Horizontally scrollable on overflow. |
| **Code block** | Monospace, `--bg-muted` background, `--text-sm`, minimal or no syntax highlighting — keep it calm. |
| **Figure** | Image + caption. Caption in `--fg-muted`, `--text-sm`, centered below image. |
| **Artifact embed** | Iframe wrapped in `<figure class="artifact-embed">` with a caption bar showing label, title, and "Open standalone ↗" link. Thin `--border` around the frame; `--bg-elevated` background. |
| **Table** | Clean. `--border` lines, `--bg-inset` header row, optional alternating row shading with `--bg-muted`. Sans, `--text-sm`. |
| **Lecture card (landing)** | `--bg-elevated` background, thin border, serif title, sans meta, hover raises slightly. |
| **Table of contents** | Sticky left sidebar on desktop; collapsed `<details>` at top on mobile. Sans, `--text-sm`. Current section highlighted in cobalt. |

### 5.6 Responsive Behavior

- **Desktop (≥1024px)**: Sticky TOC sidebar left, main content centered. Artifacts embed at their natural width.
- **Tablet (640–1023px)**: TOC collapses to a button that opens an overlay. Content full-width of the container.
- **Mobile (<640px)**: Single column, TOC at top as a collapsible `<details>`. Some artifacts (e.g., flow-cell visualization) may require horizontal scroll — don't fight it.

### 5.7 No Dark Mode

**Deliberate decision: this site is light-mode only.** Rationale:

- The "Lecture Notes" palette is a single coherent aesthetic choice. A dark-mode version would be a different design, not a toggle of the same one — the warm off-white paper tone, the cobalt accent, the base-color harmonies all were picked for this surface.
- Classroom projectors wash out dark backgrounds. Light mode reads better from the back row.
- Supporting two themes doubles the surface area of every visual decision (every new artifact, every figure, every color check). Not worth it for a course site.

**Implementation note for developers**: do not scaffold a `[data-theme="dark"]` attribute, dark-mode CSS variables, or a theme toggle. Do not use `prefers-color-scheme` media queries. Users' OS-level dark-mode preferences should have no effect on this site. If a future need emerges, it becomes a redesign task, not a feature flag.

## 6. Page Templates

### 6.1 Landing Page (`index.html`)

- Course title (serif, `--text-4xl`), instructor, semester.
- 2–3 sentence course description.
- Grid of lecture cards: number, title, duration, short topic summary, "Open lecture →" link in cobalt.
- Optional footer section: syllabus link, GitHub repo link, contact.

### 6.2 Lecture Page (`lectures/lecture-N.html`)

**Format: single long scrollable page, not slides.** The lecture is taught live by scrolling through the page (projected to the classroom) and is used by students afterwards as the authoritative study material. The page is the lecture — there is no separate slide deck or lecture-notes document.

Implications of this decision:

- Content flows as connected prose and inline artifacts. No slide-like chunking.
- All interactive artifacts are embedded inline via iframes exactly where they belong in the reading flow.
- The same page works for live teaching and after-class review.
- Progress through the lecture is indicated visually (scroll-progress bar) so students — and the lecturer — always know where they are.
- For live projection, the lecturer can invoke a `.presenter-mode` CSS class (future enhancement) that scales up type; it does not change the document structure.

Structure:

1. **Scroll-progress bar**: thin cobalt bar fixed to the top of the viewport, filling left-to-right as the reader progresses through the lecture. 2–3 px tall.
2. **Header**: lecture number (kicker), title, estimated duration.
3. **Table of contents** (sticky sidebar on desktop; collapsed on mobile).
4. **Learning objectives**: short bullet list at the top, in a bordered box.
5. **Body**: parts and sections. Artifacts embedded inline.
6. **Wrap-up**: summary, next-lecture pointer, homework, references.

Each **Part** is a `<section>` with an `id` matching the TOC anchor. Each section heading has a subtle anchor-link copy button that appears on hover ("¶" or a link icon) — clicking copies a deep link to that section to the clipboard, so students can share references to specific passages.

Each **artifact embed** looks like:

```html
<figure class="artifact-embed">
  <figcaption class="artifact-caption">
    <span class="artifact-label">Interactive</span>
    <span class="artifact-title">Central Dogma Translator</span>
    <a class="artifact-standalone"
       href="../artifacts/lecture-01/02-central-dogma.html"
       target="_blank" rel="noopener">Open standalone ↗</a>
  </figcaption>
  <iframe src="../artifacts/lecture-01/02-central-dogma.html"
          class="artifact-frame"
          loading="lazy"
          title="Central Dogma Translator"></iframe>
</figure>
```

## 7. Artifact Embedding Contract

Every artifact file must:

1. Be a fully self-contained `.html` file (one file, no external relative dependencies except `../_shared/*`).
2. Import `../_shared/artifact-theme.css` for design tokens.
3. Include `../_shared/resize.js` at the bottom so the parent page can auto-size the iframe.
4. Render sensibly at widths from 600px up to ~1400px.
5. Have a readable `<title>` and a visible heading inside (for standalone viewing).
6. Never make external network calls beyond CDN library fetches (KaTeX, Chart.js). No analytics, no trackers.

### `artifact-theme.css` content

This file should contain:

- The full `:root` custom-property block from §5.2 (color tokens).
- The font-family variables from §5.3.
- A minimal reset + base styles (`body { font-family: var(--font-sans); color: var(--fg); background: var(--bg); }`).
- No component-level styles — those belong in each artifact.

## 8. Hosting

Recommended: **GitHub Pages**. Free, simple, no config beyond enabling Pages on the repo. URL: `https://<user>.github.io/bioinformatics-course/`.

Alternatives: Netlify (nicer preview deployments) or university Apache — rsync the folder.

All paths relative (no absolute URLs) so the site works on any subpath.

## 9. Deliverable Order (Suggested)

1. **Scaffold**: `index.html`, `styles.css`, `lecture.css`, `artifact-theme.css`, `resize.js`, empty `lecture-01.html` skeleton with the TOC and part structure wired up.
2. **Design demo**: Build artifact `artifacts/lecture-01/02-central-dogma.html` first as the design-language proof. Review with stakeholder before proceeding.
3. **Lecture 1 prose**: Fill `lecture-01.html` from the content in `lecture-01.md`.
4. **Remaining artifacts**: Build `01`, `03`, `04`, `05`, `06` in order of rising complexity. `04` (Illumina basecaller) and `05` (Nanopore squiggle) are the most involved.

## 10. Open Questions / Decisions Deferred

- **Math rendering**: KaTeX (faster, smaller) vs MathJax (more complete). **Recommend KaTeX** — we won't hit its limits in this course.
- **Search**: not needed for v1. If added, a static index (Lunr/Fuse) works without a backend.
- **Analytics**: none by default. If the university requires it, use Plausible or a self-hosted equivalent — not GA.
- **Accessibility**: aim for WCAG AA on prose pages. Artifacts will vary; where possible, provide a text-based fallback description in the surrounding lecture prose.
- **Font hosting**: Google Fonts works out of the box; self-host before first production deploy for reliability.
