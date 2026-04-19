# bioinformatics-web

Lecture material and companion website for a bioinformatics course taught at a School of Electrical Engineering. Each lecture is a long scrollable HTML page with inline interactive artifacts and static figures — it's the same document the lecturer projects in class and the students read afterwards.

## Run it locally

No build step. Any static file server works:

```
cd bioinformatics-web
python3 -m http.server 8000
```

Open <http://localhost:8000/>. Alternatives: `npx serve .` if you have Node, or the VS Code "Live Server" extension. Opening files directly from disk (`file://`) mostly works but iframe auto-resize and some browser security policies are more reliable over `http://`.

## Repository layout

```
.
├── index.html                            # Landing page
├── lectures/
│   ├── lecture-01.html                   # Long-form lecture pages
│   └── lecture-02.html
├── artifacts/
│   ├── _shared/                          # Shared design tokens + postMessage resizer
│   ├── lecture-01/                       # 6 interactive HTML artifacts per lecture
│   └── lecture-02/
├── diagrams/
│   ├── lecture-01/                       # Static SVG figures (10-12 per lecture)
│   └── lecture-02/
├── assets/
│   ├── styles.css                        # Site-wide tokens + homepage
│   ├── lecture.css                       # Lecture-page layout, callouts, figures
│   └── logo/                             # Mark, wordmark, favicon variants
├── markdown_resources/                   # Specs + source-of-truth content
│   ├── website-spec.md  homepage-spec.md  logo-spec.md        # Course-wide specs
│   ├── lecture-style-guide.md  diagram-style-guide.md
│   ├── lesson1_md_files/                                      # Lecture 1 sources
│   │   ├── artifacts-spec.md  figures-spec.md  lecture-01.md
│   └── lesson2_md_files/                                      # Lecture 2 sources
│       ├── artifacts-spec.md  figures-spec.md  lecture-02.md
└── proof_html_resources/                 # Standalone-viewable design references
    ├── fastq-anatomy-proof.html
    └── logo-proof-v2.html
```

Each lecture embeds its artifacts via `<iframe>` from `artifacts/lecture-NN/` and its figures via `<img>` from `diagrams/lecture-NN/`. All paths are relative — the site runs under any subpath.

## Specs & authoring guides

These are the sources of truth when adding content. Read the relevant one before editing code.

| File | Use for |
|---|---|
| `markdown_resources/website-spec.md` | Site architecture, folder layout, design tokens, embedding contracts |
| `markdown_resources/homepage-spec.md` | Landing page layout and 16-card lecture list |
| `markdown_resources/logo-spec.md` | Course mark and wordmark |
| `markdown_resources/lecture-style-guide.md` | Voice, structure, callout types, and formatting for lecture prose |
| `markdown_resources/diagram-style-guide.md` | House style for SVG diagrams (colour, stroke, typography, arrows) |
| `markdown_resources/generate_new_lesson_flow.md` | End-to-end recipe for generating a new lecture (phases, steps, validations) |
| `markdown_resources/lesson1_md_files/` | Lecture 1 — behaviour specs and source content |
| `markdown_resources/lesson2_md_files/` | Lecture 2 — behaviour specs and source content |
| `proof_html_resources/fastq-anatomy-proof.html` | Canonical SVG reference — copy arrow markers and typography patterns from here |
| `proof_html_resources/logo-proof-v2.html` | Logo visual reference; course mark is locked to Variant A3 |

## Design system in one glance

- **Light mode only.** No dark-mode scaffolding, no `prefers-color-scheme` queries. See `website-spec.md` §5.7 for why.
- **Three fonts.** Source Serif 4 for headings, Inter for body/UI, JetBrains Mono for any biological sequence or code.
- **One accent.** Deep cobalt (`--accent: #1e3a8a`). Used sparingly.
- **Base colours are locked.** `--base-a`/`--base-t`/`--base-g`/`--base-c`/`--base-u`/`--base-n` — students build visual memory across artifacts.
- **Zero build step.** Plain HTML, CSS, vanilla JS. Chart.js and KaTeX via CDN where needed.

## Current status

- **Lecture 1 — Foundations: From Cells to Sequences to FASTQ.** 3h 35min · 12 figures · 6 interactive tools. Shipped.
- **Lecture 2 — Read Alignment: From Brute Force to FM Index and Back.** 3h 50min · 10 figures · 6 interactive tools. Shipped.
- **Lectures 3–16** — placeholders on the homepage; not yet started.

Lectures shipped so far: **2 of 16.** Interactive tools: **12.** Figures: **22.**

## License

Apache 2.0 (course code and tooling). See `LICENSE` for the full text.

The course mark and wordmark carry their own identity restrictions — see `assets/logo/LICENSE.md`.
