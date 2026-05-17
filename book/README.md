# Book

Print edition of *Bioinformatics for Engineers* — A4, trade-paperback technical-book layout, built from Typst sources.

The book expands on the lectures (longer prose, history, derivations, side-notes, exercises, references); it is **not** a one-to-one transcript. Interactive HTML artifacts from the website are deliberately omitted — they are replaced by additional static figures designed specifically for print.

## Build

```bash
# One-time install
brew install typst

# Compile to PDF
typst compile --root . book/book.typ book/build/book.pdf

# Watch + recompile on save (great for iterating on a chapter)
typst watch --root . book/book.typ book/build/book.pdf
```

`--root .` runs from the repo root so the chapter files can reference figures both under `book/figures/` (new, print-only) and `diagrams/lecture-NN/` (existing, shared with the website).

## Layout

```
book/
├── book.typ                       # Top-level: frontmatter, TOC, chapter includes
├── theme/
│   └── book-theme.typ             # Page geometry, type, headings, admonitions, captions
├── chapters/
│   ├── ch01-foundations.typ
│   └── …                          # Added incrementally, one chapter at a time
├── figures/
│   └── ch01/                      # New SVG figures unique to the print edition
│       ├── f1-dna-grooves.svg
│       ├── f2-cost-timeline.svg
│       ├── f3-intensity-trace.svg
│       └── f4-phred-encoding.svg
└── build/                         # Generated PDFs (gitignored)
```

## Authoring conventions

- Chapter files import `theme/book-theme.typ` to get the admonition functions: `#note[…]`, `#tip[…]`, `#warn[…]`, `#danger[…]`, and `#matters[…]` (the book-specific "Why this matters" pre-chapter callout).
- Reference existing site figures by relative path: `image("../../diagrams/lecture-NN/<file>.svg")`.
- Add new figures under `book/figures/chNN/` and reference with `image("../figures/chNN/<file>.svg")`.
- All figures match the website's house style (see `markdown_resources/diagram-style-guide.md`).

## Fonts

The theme prefers Source Serif 4 / Inter / JetBrains Mono — the same fonts the website uses. The fallback chain is Charter / Helvetica Neue / Menlo (all preinstalled on macOS), so the build works out of the box without installing anything beyond Typst. The font-fallback warnings during compile are expected and harmless.

To use the preferred fonts:
```bash
brew install --cask font-source-serif-4 font-inter font-jetbrains-mono
```

Then re-run `typst compile` and the warnings disappear.

## Status

| Chapter | Source lecture | Word count (approx.) | Figures | Status |
|---------|----------------|----------------------|---------|--------|
| 1. Foundations | L01 | ~9,000 | 12 existing + 4 new = 16 | First draft |
| 2 – 27         | L02 – L27 | —                    | —       | Not started |

The plan is one chapter at a time: write, build, iterate on layout, lock the format, then move on. Compiling the full book is the very last step.
