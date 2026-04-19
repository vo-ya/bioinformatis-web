# Homepage Specification

> **Purpose**: Detailed build spec for the course homepage (`index.html`).
> **Companion files**: `website-spec.md` (site-wide design system), `logo-spec.md` (mark and wordmark), `lecture-style-guide.md` (for lecture page links).

---

## 1. Philosophy

The homepage has two audiences whose needs mostly align but sometimes pull apart:

- **Enrolled students** want to get to Lecture N as fast as possible. They've been here before.
- **Prospective students, colleagues, and curious visitors** want enough context to understand what the course is.

The layout satisfies both by being **content-primary and vertically stacked**: a short identity block establishes context for new visitors without blocking returning students from reaching the lecture grid, which is always the primary content on the page.

The page is **editorial, not promotional**. No marketing language, no carousel, no hero image, no "Welcome!" paragraph. It reads like the front matter of a serious book.

---

## 2. Layout — Direction A (Editorial, Single Column)

The page is centered, single-column, with vertical rhythm as the only structuring mechanism. Max content width 720px, matching lecture prose. The layout flows:

```
  [Logo mark + wordmark]
  
  COURSE TITLE (large serif)
  Subtitle / tagline (regular)
  
  Course description — 3–4 sentences of prose
  
  Course metadata — institution, semester, instructor
  
  ——————————————————— divider ———————————————————
  
  LECTURES (small kicker label)
  
  [Lecture 01 card]
  [Lecture 02 card]
  [Lecture 03 card]
  ...
  [Lecture 16 card]
  
  ——————————————————— divider ———————————————————
  
  Secondary nav: Syllabus · GitHub · Contact
  
  Footer line
```

---

## 3. Content Blocks (in order)

### 3.1 Logo mark + wordmark

Placed at the very top, left-aligned to the content column.

- The helix-as-waveform mark at **28px** height (header size from `logo-spec.md`)
- Wordmark to the right of the mark:
  - Line 1: `Bioinformatics for Engineers` — Source Serif 4, weight 500, size 16px
  - Line 2: `University of Belgrade · School of Electrical Engineering` — Inter, weight 400, size 12px, color `--fg-muted`
- 12px horizontal gap between mark and wordmark
- 48–64px of vertical space below this block before the title

### 3.2 Course title

Large, serif, confident. This is the visual anchor of the page.

- Text: `Bioinformatics for Engineers`
- Font: Source Serif 4, weight 500
- Size: `--text-4xl` (44px at base scale)
- Line height: `--leading-tight` (1.2)
- Letter-spacing: slight negative (−0.01em) for optical tightness at this size
- Color: `--fg`
- Margin below: 8px

### 3.3 Subtitle

A one-line tagline. Optional but recommended — it gives the prospective visitor a two-second read of what the course is about.

- Text: `A one-semester course on how biology became a data problem.`
- Font: Inter, weight 400, size 18px
- Color: `--fg-muted`
- Margin below: 32px

### 3.4 Course description

The 2–4 sentence orientation block.

- Text (placeholder, as approved):
  > A one-semester course on how biology became a data problem. Taught at the School of Electrical Engineering, it covers sequencing technology, the bioinformatics pipeline from raw reads to biological insight, and the algorithms and instruments that made modern genomics possible — framed throughout from an engineer's perspective.
- Font: Inter, weight 400, size `--text-base` (16px)
- Line height: `--leading-relaxed` (1.75)
- Color: `--fg`
- Margin below: 32px

### 3.5 Course metadata

Three short lines, each a key/value pair.

```
Institution  School of Electrical Engineering, University of Belgrade
Semester     Spring 2026
Instructor   Vojislav Varjačić
Course code  [TBD]
```

- Rendered as a two-column grid-like layout (or a plain `<dl>`)
- Labels (left column): Inter, weight 500, size 13px, color `--fg-subtle`, uppercase, letter-spacing 0.08em
- Values (right column): Inter, weight 400, size 14px, color `--fg`
- Rows separated by 8px vertical spacing
- Column gap: 32px
- Labels left-aligned; values left-aligned, starting at a consistent x-position so they form a clean value column
- Note: "Course code" shows `[TBD]` in `--fg-subtle` until a real code is provided

### 3.6 Divider

A horizontal rule that separates the front matter from the lecture list.

- Thin 1px line in `--border-strong`
- 48px margin top, 48px margin bottom
- Full width of the content column (up to max-width 720px)

### 3.7 Lectures kicker + heading

- Kicker: `LECTURES` in Inter, weight 600, size 11px, letter-spacing 0.12em, color `--accent`, uppercase
- Margin below kicker: 8px
- No visible heading below it — the kicker plus the cards is enough. Do not add a redundant "Lecture List" h2.
- Margin below kicker group: 16px

### 3.8 Lecture cards — the 16-card grid

Sixteen cards, stacked vertically, one per row. Each card links to its lecture page.

**Card states**:

- **Published (Lecture 01)**: card is fully linked, shows real metadata
- **Placeholder (Lectures 02–16)**: card is visually muted, non-interactive, shows "Coming soon" badge instead of duration

**Card structure** (published):

```
┌─────────────────────────────────────────────────────┐
│  01   Foundations — From Cells to Sequences to FASTQ│
│       DNA and the central dogma, sequencing         │
│       history, NGS, long-read platforms, FASTQ.     │
│       ≈ 3h 35min · 12 figures · 6 interactive tools │
└─────────────────────────────────────────────────────┘
```

**Card structure** (placeholder):

```
┌─────────────────────────────────────────────────────┐
│  02   [Lecture 02 — TBD]                             │
│       Topic to be announced.                         │
│       Coming soon                                    │
└─────────────────────────────────────────────────────┘
```

**Card anatomy**:
- **Number** (left): Inter, weight 500, size 20px, color `--fg-muted` (published) or `--fg-subtle` (placeholder). Zero-padded ("01", not "1"). Fixed width ~44px so all numbers align.
- **Title** (main): Source Serif 4, weight 500, size 18px, color `--fg` (published) or `--fg-muted` (placeholder). Two-line max with overflow ellipsis.
- **Summary** (below title): Inter, weight 400, size 14px, color `--fg-muted`. One or two lines. Descriptive, not promotional.
- **Metadata line** (bottom): Inter, weight 400, size 12px, color `--fg-subtle`. Shows duration, figure count, artifact count separated by `·`. For placeholders: shows "Coming soon" only.

**Card styling**:
- Background: `--bg-elevated` (pure white)
- Border: 1px solid `--border`
- Border radius: 6px
- Padding: 20px 24px
- Hover (published only): border becomes `--accent`, card shifts up 1px, transition 120ms ease-out
- Placeholder cards: same base styling, but background `--bg-muted`, no hover effect, cursor `default`
- Link wraps the whole card (entire card is clickable), not just the title

**Card spacing**: 12px vertical gap between cards.

**Link target**:
- Published: `lectures/lecture-01.html` (relative path)
- Placeholder: no link, `<div>` instead of `<a>`, or `<a href="#" aria-disabled="true" tabindex="-1">` for structure consistency

### 3.9 Divider (bottom)

Same as §3.6. Separates lecture grid from secondary nav.

### 3.10 Secondary nav

Three inline links, separated by bullet characters.

```
Syllabus  ·  GitHub  ·  Contact
```

- Font: Inter, weight 400, size 14px
- Link color: `--accent-bright`, no underline at rest, underline on hover
- Separator bullet: `--fg-subtle`
- Centered horizontally
- Links:
  - `Syllabus` → to be filled in later; for now link to `#` with a `data-state="placeholder"` attribute
  - `GitHub` → repo URL; for now link to `#`
  - `Contact` → `mailto:` link or `#`; for now link to `#`
- Margin: 32px top, 16px bottom

### 3.11 Footer

A single muted line at the very bottom.

```
© 2026 Vojislav Varjačić · School of Electrical Engineering, University of Belgrade · Last updated [auto]
```

- Font: Inter, weight 400, size 12px
- Color: `--fg-subtle`
- Centered horizontally
- The "Last updated" timestamp is a static placeholder for v1; in v2 can be auto-generated at build time

---

## 4. Responsive Behavior

- **Desktop (≥1024px)**: Content column centered, max-width 720px, 80px horizontal padding on the page body.
- **Tablet (640–1023px)**: Same max-width; horizontal page padding reduced to 48px.
- **Mobile (<640px)**:
  - Horizontal padding: 24px
  - Course title drops to `--text-3xl` (32px)
  - Lecture card number drops to size 18px; card padding reduced to 16px 20px
  - Metadata block: labels and values stack vertically instead of two columns
  - Logo block: wordmark second line may truncate to just "School of EE"

---

## 5. Accessibility

- Every card link has proper focus-visible styling: 2px `--accent` outline with 2px offset.
- Card links have an `aria-label` that combines number, title, and duration: `Lecture 01: Foundations — From Cells to Sequences to FASTQ. Duration approximately 3 hours 35 minutes.`
- Placeholder cards use `aria-disabled="true"` and include `aria-label="Lecture 02 — coming soon"`.
- The lecture list is wrapped in a `<nav aria-label="Course lectures">` or a semantic `<ol>` so screen readers announce it as a structured list.
- All text contrast combinations meet WCAG AA (already ensured by the design tokens).
- Skip-to-content link at the top of the page for keyboard users, visible on focus.

---

## 6. Lecture Card Data — Source of Truth

The 16 cards are populated from a static array inside `index.html`. The first entry is populated with real data; entries 2–16 are placeholder objects. No external data file — the array lives inline in a `<script>` block or, more simply, as hand-written HTML cards since there are only 16 of them and they don't change often.

```js
const lectures = [
  {
    num: "01",
    title: "Foundations — From Cells to Sequences to FASTQ",
    summary: "DNA and the central dogma, sequencing history, NGS, long-read platforms, FASTQ.",
    duration: "≈ 3h 35min",
    figures: 12,
    artifacts: 6,
    href: "lectures/lecture-01.html",
    state: "published",
  },
  // Placeholder entries 02 through 16
  { num: "02", title: "[Lecture 02 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "03", title: "[Lecture 03 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "04", title: "[Lecture 04 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "05", title: "[Lecture 05 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "06", title: "[Lecture 06 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "07", title: "[Lecture 07 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "08", title: "[Lecture 08 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "09", title: "[Lecture 09 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "10", title: "[Lecture 10 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "11", title: "[Lecture 11 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "12", title: "[Lecture 12 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "13", title: "[Lecture 13 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "14", title: "[Lecture 14 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "15", title: "[Lecture 15 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
  { num: "16", title: "[Lecture 16 — TBD]", summary: "Topic to be announced.", state: "placeholder" },
];
```

Either hand-write the 16 cards directly in the HTML (preferred — simpler, no JS on the homepage) or generate them from this array in a `<script>` block. Either approach is fine. If generated, the JS must run without any dependencies.

---

## 7. What NOT to include on the homepage

- No hero illustration or large decorative image.
- No "Welcome to my course!" paragraph.
- No social media links or share buttons.
- No announcement banner, news ticker, or "what's new" section. If news matters later, add it as a separate page.
- No testimonials, no "Featured in..." strip.
- No course rating, stars, or reviews.
- No ads, analytics beacons, or trackers.
- No auto-playing anything.
- No popup, modal, or cookie banner (if the university legally requires a cookie banner, add only what's legally required — no more).

---

## 8. Implementation notes

- Single `index.html` file, roughly 250–300 lines of HTML including inline SVG logo mark.
- All styles via `assets/styles.css` (shared) — no page-specific stylesheet needed.
- The logo mark SVG is inlined in the HTML so it can use currentColor / CSS variables directly if we ever want to.
- Favicon: use the favicon variant from `logo-spec.md`, placed at `/favicon.svg` with a PNG fallback at `/favicon.png` (32×32). Reference both in the `<head>`.
- `<meta>` tags: `description`, `og:title`, `og:description`, `og:image` (use a PNG export of the logo display variant at 1200×630). Charset UTF-8, viewport responsive.
- Lang attribute on `<html>`: `lang="en"`.
