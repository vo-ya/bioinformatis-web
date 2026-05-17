# Book Authoring Style Guide

> **Purpose**: Reference this file when drafting any chapter of the print edition of *Bioinformatics for Engineers*.
> **How to use in a new chat**: paste this file as context at the start, along with which lecture you are converting and the chapter number. Tell the assistant: "Draft chapter NN following this style guide."
> **Companion files**: `lecture-style-guide.md` (source-of-truth for voice), `diagram-style-guide.md` (figure conventions; print figures inherit the same house style).

---

## 1. Project Identity — What the Book Is and Isn't

The book is the **print edition** of the 27-lecture course already shipped at <https://vo-ya.github.io/bioinformatis-web/>. It is:

- **An expansion**, not a transcript. The website lecture is the spine; the book chapter expands it with history, deeper math, sidebars, exercises, and further-reading references. Aim for ~1.5–2× the word count of the lecture's prose.
- **Print-only typography**. A4 trim, serif body, mirrored running heads, justified prose with first-line indents — the classic trade-paperback technical-book aesthetic. Not the website. Don't try to make the book look like the website.
- **Self-contained per chapter**. A reader who picks up the book mid-stack at Chapter 11 should be able to follow it. Reference earlier chapters by name and number; don't assume continuity of context.

What the book is *not*:

- A transcript or copy-paste of the website's HTML. The website's interactive artifacts are deliberately omitted in print and **replaced with additional static figures** specific to the book (`book/figures/chNN/`).
- A re-skinning of an O'Reilly / Packt / Manning template. The look is *inspired by* trade-paperback conventions but uses our own typography and colour palette.

---

## 2. Toolchain

- **Format**: Typst 0.14+ (https://typst.app/). Source files end in `.typ`.
- **Trim**: A4 (210 × 297 mm).
- **Build**: `typst compile --root . book/book.typ book/build/book.pdf` from the repo root. `--root .` is required so chapter files can reference shared figures under `diagrams/lecture-NN/`.
- **Live preview**: `typst watch --root . book/book.typ book/build/book.pdf` re-renders on save.
- **Output**: a single PDF in `book/build/` (gitignored).

Install Typst with `brew install typst`. No LaTeX, no Ruby, no Pandoc.

---

## 3. Directory Layout

```
book/
├── book.typ                       Top-level: frontmatter, TOC, chapter includes
├── theme/
│   └── book-theme.typ             Theme module (page geometry, type, admonitions)
├── chapters/
│   └── chNN-<slug>.typ            One file per chapter, slug from lecture title
├── figures/
│   └── chNN/                      Print-only SVGs (no equivalent on the website)
└── build/                         Generated PDFs — gitignored
```

Existing site figures stay where they are: chapters reference them with
`image("../../diagrams/lecture-NN/<file>.svg")`. New, print-only figures live
under `book/figures/chNN/` and are referenced as
`image("../figures/chNN/<file>.svg")`.

---

## 4. Typography (codified in `book/theme/book-theme.typ`)

The theme module is the executable source of truth for everything below. Match it; don't try to override it inside a chapter file unless there's a specific reason.

### 4.1 Page geometry

- Trim: A4 (210 × 297 mm).
- Margins: 25 mm top/bottom, 30 mm inside (binding), 22 mm outside.
- Running heads: mirrored across the spine. Book title sits on the outside edge; chapter title on the inside edge. Both in small-caps Inter, muted grey. A 0.5 pt hairline rule sits below the heads.
- Folios (page numbers): outside edge of the footer, Inter 9 pt.

### 4.2 Font stack

| Role | Preferred | First fallback | OS fallback |
|---|---|---|---|
| Body serif | Source Serif 4 | Charter | Georgia |
| Heading / UI sans | Inter | Helvetica Neue | Arial |
| Code / mono | JetBrains Mono | Menlo | Monaco |

The theme declares these as fallback arrays so the build works on any default macOS install without extra font packages. To use the preferred fonts: `brew install --cask font-source-serif-4 font-inter font-jetbrains-mono`.

### 4.3 Body type

- 10.5 pt Source Serif 4, leading 0.75 em, justified.
- Paragraph spacing: 0.65 em.
- First-line indent: **1.2 em, but only on paragraphs that follow another paragraph** (`first-line-indent: (amount: 1.2em, all: false)`). The first paragraph after every heading, figure, list, code block, or admonition stays flush left. This is the standard book-typography convention.
- No widow/orphan tuning yet — revisit once a few chapters land.

### 4.4 Heading hierarchy

| Level | Use | Style | Numbering |
|---|---|---|---|
| `=` (h1) | Chapter | New page; uppercase cobalt kicker "CHAPTER N"; 30 pt serif title; 60 pt cobalt rule below; 36 pt clear space before body | "Chapter 1", "Chapter 2", …  |
| `==` (h2) | Section | 18 pt serif, weight medium; inline cobalt **N.M** prefix; 32 pt above / 14 pt below | 1.1, 1.2, 1.3 …  |
| `===` (h3) | Subsection | 13 pt italic serif; no number; 20 pt above / 8 pt below | none |

Auto-numbering is wired via `set heading(numbering: "1.1")` in the theme; the show rules render the counter manually so the number lands exactly where we want it (kicker for h1, inline prefix for h2, hidden for h3).

### 4.5 Figures

- Captioned with **Figure N.M** prefix in serif italic, Inter 9.5 pt for the rest.
- Indented 2 em on each side of the column for visual breathing room.
- Use `<fig:slug>` labels for cross-references.
- Width: most figures sit at 90–95 % column width; smaller ones at 75–80 %.

### 4.6 Code listings

- Block code: light-grey background (`#f4f3ee`), 0.5 pt hairline border, 10 pt internal padding, JetBrains Mono 9 pt.
- Inline code: same background, 0.92 em mono, tight horizontal padding.

### 4.7 Lists

- Bullets default to `•` (level 1) and en-dash `–` (level 2).
- Numbered lists use `1., 2., 3.` (no parens, no hash).

---

## 5. Admonitions (sidebar callouts)

The theme exposes five callout functions. Use them sparingly — they're a tool for *real* asides, not for highlighting every paragraph. As a rule of thumb, no more than three admonitions per ten pages of prose.

| Function | Colour | Semantic | When to use |
|---|---|---|---|
| `#note[…]` | cobalt | factual aside | A piece of information the reader will want but that interrupts the main flow if inlined. |
| `#tip[…]` | green | practical tip | A specific actionable suggestion ("the single most useful early-pipeline reality check is `seqkit stats`…"). |
| `#warn[…]` | amber | pitfall warning | A common mistake or trap ("Historical Phred+64 encoding is still occasionally found in archived data…"). |
| `#danger[…]` | red | strong caution | A data-destroying or correctness-fatal class of mistake. Use sparingly; `#warn` is usually enough. |
| `#matters[…]` | violet | "Why this matters" | A book-only callout at the top of a chapter or major section explaining stakes and motivation. Typically appears at most twice per chapter. |

Each callout renders as a coloured-left-rule box with an uppercase label. Body text is 10 pt, slightly tighter leading than body prose.

---

## 6. Chapter Structure

Every chapter follows the same backbone:

1. **Title (h1)**: short, descriptive. Match the lecture's `<h1>` title verbatim. Format: `Topic: Tagline`.
2. **Why this matters** (`#matters[…]`, immediately under the chapter title): 100–200 words. Frame the stakes. No filler.
3. **Opening paragraph**: a single non-section paragraph that sets up the chapter before §1. Often a vivid concrete grounding ("A genome is not, strictly speaking, an information object…").
4. **Sections** (h2): typically 5–8 per chapter. Each section is a self-contained topic. Numbered 1.1, 1.2, ….
5. **Subsections** (h3): use sparingly — only when a section genuinely contains two distinct sub-topics. Don't decorate every section with subsections.
6. **Summary** (h2, last functional section): bullet list of 4–6 takeaways. Short sentences.
7. **Exercises** (h2): 5–8 mixed computational and conceptual problems. Number with `**1.** *Title.*` in the body — Typst won't auto-number these as headings.
8. **Further Reading** (h2): 4–6 references. Format: `*Authors* (Year). "Title." *Journal* vol: pages. One-sentence note.`

Length target: **5,000–12,000 words** of prose per chapter, plus figures. Lecture 1 is on the long end (≈9,000 words) because it has to onboard biology; later chapters land closer to 6,000.

---

## 7. Voice

Inherit everything in `lecture-style-guide.md` §2 (Voice & Tone). The book version is the same voice but:

- **Slightly more formal in register** — written prose, not spoken. Fewer parenthetical asides, more careful syntax.
- **More historical depth allowed** — the book has space for "Friedrich Miescher isolated 'nuclein' from pus-soaked bandages in 1869…" passages that the lecture would compress.
- **More explicit derivations** — show the math when there's space. "Phred score Q = -10 log₁₀ P_err" earns a full paragraph in the book; the lecture would compress it to one line.
- **Address the reader, not the audience** — "you will spend most of your bioinformatics career working with eukaryotic genomes" not "we will see in this chapter that…".

Don't preface chapters with "In this chapter you will learn…" boilerplate. The chapter title, *Why this matters* box, and opening paragraph carry that load.

---

## 8. Figures: Existing + Print-Only

Two sources:

- **Existing site figures** (`diagrams/lecture-NN/*.svg`) — already in the house style, reuse them. Reference with `image("../../diagrams/lecture-NN/<file>.svg")` from within `chapters/`.
- **Print-only new figures** (`book/figures/chNN/*.svg`) — designed specifically to replace the website's interactive artifacts. Use the same house style as `diagram-style-guide.md`: warm off-white background (`#fafaf9`), cool-grey default stroke (`#374151`), Inter for labels, JetBrains Mono for sequences and numerics. Project palette only — base colours, `#1e3a8a` cobalt, `#b45309` amber, `#065f46` green, `#dc2626` red, `#7c3aed` violet. No Microsoft-Office defaults.

A typical chapter has 10–13 existing figures plus 3–6 new ones (one per displaced artifact). Always check `diagrams/lecture-NN/` before drafting new figures — there may already be one.

### Caption authoring

- One sentence; two only if the second clarifies a non-obvious detail.
- Active voice ("The polymerase adds…" not "Bases are added by…").
- No "above" / "below" / "left" / "right" — the figure may reflow on print.
- Always reference figures with a label (`<fig:slug>`) and cross-reference with `@fig:slug`.

---

## 9. Exercises

Exercise authoring follows the lecture's pedagogy but with print conventions:

- Number with **boldface ordinal** (`**1.** *Computational.*`) — these aren't real headings.
- Mix **computational** (code-something problems) and **conceptual** (interpret-something problems).
- Include one **open-ended** problem at the end ("Pick a published tool that explicitly accounts for X. Describe in one paragraph the correction it applies.").
- No answer key in the book itself — answers live in the website's Colab notebooks (linked once at the end of the chapter).

---

## 10. Cross-References

- Within a chapter: `@sec:slug` or `@fig:slug` (auto-rendered as "Section 1.3" / "Figure 1.4").
- Across chapters: spell out the chapter number — "Chapter 4 covers variant calling".
- To the website: footnote-style reference at first mention only — "the live Colab exercise (`vo-ya.github.io/bioinformatis-web/`)".

Avoid forward-pointers ("we'll see this in Chapter 22") unless they're load-bearing — they age poorly when chapters get reordered.

---

## 11. Glossary, Index, Notation

Deferred for the first pass. Once 5–6 chapters land we'll know what to factor out:

- **Glossary**: high-frequency technical terms. Build incrementally.
- **Index**: page references for key concepts. Typst supports indexing via labels.
- **Notation table**: standard symbols at the front of the book.

Don't add a glossary entry while drafting a chapter — flag candidate terms with a brief inline definition on first use and let the editorial pass collect them.

---

## 12. Pre-Submission Checklist

Before merging a new chapter to `main`:

- [ ] Chapter file compiles cleanly: `typst compile --root . book/book.typ book/build/book.pdf`
- [ ] All figures resolve (no broken image paths in the output)
- [ ] All `@fig:*` and `@sec:*` cross-references resolve
- [ ] Word count in the 5,000–12,000 range
- [ ] At least one **Why this matters** callout
- [ ] No `#danger` callouts (escalate to `#warn` unless truly data-destroying)
- [ ] Exercises section has ≥ 5 problems and one open-ended
- [ ] Further Reading section has ≥ 4 entries
- [ ] No O'Reilly / Packt / Manning trademark references in prose, comments, or filenames

---

## 13. Build & Commit

One chapter at a time. Each chapter is its own PR. Do **not** compile the full book until the structural conventions are locked across 3+ chapters — until then, the per-chapter PDF is what gets reviewed.
