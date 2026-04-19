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
├── index.html                      # Landing page
├── lectures/
│   └── lecture-01.html             # Long-form lecture pages
├── artifacts/
│   ├── _shared/                    # Shared design tokens + postMessage resizer
│   └── lecture-01/                 # One interactive artifact per .html file
├── diagrams/
│   └── lecture-01/                 # Static SVG figures (12 per lecture)
└── assets/
    ├── styles.css                  # Site-wide tokens and landing page
    └── lecture.css                 # Lecture-page layout, callouts, figures
```

Each lecture embeds its artifacts via `<iframe>` from `artifacts/lecture-NN/` and its figures via `<img>` from `diagrams/lecture-NN/`. All paths are relative — the site runs under any subpath.

## Specs & authoring guides

These are the sources of truth when adding content. Read the relevant one before editing code. All spec and style-guide markdown lives under `markdown_resources/`; the two HTML proof files stay at the repo root because they double as standalone-viewable design references.

| File | Use for |
|---|---|
| `markdown_resources/website-spec.md` | Site architecture, folder layout, design tokens, embedding contracts |
| `markdown_resources/homepage-spec.md` | Landing page layout and 16-card lecture list |
| `markdown_resources/logo-spec.md` | Course mark and wordmark |
| `markdown_resources/lecture-style-guide.md` | Voice, structure, callout types, and formatting for lecture prose |
| `markdown_resources/artifacts-spec.md` | Per-artifact behaviour + UI specifications |
| `markdown_resources/diagram-style-guide.md` | House style for SVG diagrams (colour, stroke, typography, arrows) |
| `markdown_resources/figures-spec.md` | Per-lecture figure list with build-order recommendations |
| `markdown_resources/lecture-01.md` | Source markdown for Lecture 1 content |
| `proof_html_resources/fastq-anatomy-proof.html` | Canonical SVG reference — copy arrow markers and typography patterns from here |
| `proof_html_resources/logo-proof-v2.html` | Logo visual reference; course mark is locked to Variant A3 |

## Design system in one glance

- **Light mode only.** No dark-mode scaffolding, no `prefers-color-scheme` queries. See `website-spec.md` §5.7 for why.
- **Three fonts.** Source Serif 4 for headings, Inter for body/UI, JetBrains Mono for any biological sequence or code.
- **One accent.** Deep cobalt (`--accent: #1e3a8a`). Used sparingly.
- **Base colours are locked.** `--base-a`/`--base-t`/`--base-g`/`--base-c`/`--base-u`/`--base-n` — students build visual memory across artifacts.
- **Zero build step.** Plain HTML, CSS, vanilla JS. Chart.js and KaTeX via CDN where needed.

## Current status

- **Lecture 1 — Foundations: From Cells to Sequences to FASTQ.** Prose, six interactive artifacts, and twelve figures all shipped.
- **Lecture 2+ —** not yet started.

## License

Apache 2.0. See `LICENSE`.
