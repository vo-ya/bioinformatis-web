// ─── Book theme ──────────────────────────────────────────────────────
// Trade-paperback technical-book look for the Bioinformatics for
// Engineers print edition. A4 trim, serif body, sans heads, mono code,
// generous margins, mirrored running heads.

#let accent = rgb("#1e3a8a")           // deep cobalt (matches the website)
#let fg     = rgb("#0a0a0a")
#let muted  = rgb("#525252")
#let subtle = rgb("#8a8a85")
#let rule   = rgb("#e5e3dc")
#let code-bg = rgb("#f4f3ee")

#let book-doc(
  title: "Bioinformatics for Engineers",
  subtitle: "How biology became a data problem.",
  author: "Vojislav Varjačić",
  body
) = {
  set document(title: title, author: author)

  set page(
    paper: "a4",
    margin: (top: 25mm, bottom: 25mm, inside: 30mm, outside: 22mm),
    header: context {
      let n = counter(page).get().first()
      if n <= 1 { return }
      let chapters = query(selector(heading.where(level: 1)).before(here()))
      if chapters.len() == 0 { return }
      let chap-title = chapters.last().body
      let book-side  = text(size: 9pt, font: ("Inter", "Helvetica Neue", "Arial"), fill: subtle, smallcaps(title))
      let chap-side  = text(size: 9pt, font: ("Inter", "Helvetica Neue", "Arial"), fill: subtle, smallcaps(chap-title))
      if calc.even(n) {
        grid(columns: (1fr, 1fr), align: (left, right), book-side, chap-side)
      } else {
        grid(columns: (1fr, 1fr), align: (left, right), chap-side, book-side)
      }
      v(-2pt)
      line(length: 100%, stroke: 0.5pt + rule)
    },
    footer: context {
      let n = counter(page).get().first()
      let label = text(size: 9pt, font: ("Inter", "Helvetica Neue", "Arial"), fill: subtle, str(n))
      align(if calc.even(n) { left } else { right }, label)
    },
  )

  // Body type — serif for prose; size tuned for A4 trade-paperback feel.
  set text(font: ("Source Serif 4", "Charter", "Georgia"), size: 10.5pt, lang: "en")
  // Book convention for indentation: indent every paragraph *except* the
  // first one after a heading, figure, list, code block, or callout
  // ("all: false"). This is the classic O-style novel/textbook indent.
  set par(
    leading: 0.75em,
    justify: true,
    first-line-indent: (amount: 1.2em, all: false),
    spacing: 0.65em,
  )

  // Auto-number headings; show rules below choose where each level's
  // number appears (kicker for h1, inline cobalt prefix for h2, hidden
  // for h3).
  set heading(numbering: "1.1")

  // (No per-element indent overrides needed — `first-line-indent.all =
  // false` above already exempts first paragraphs after every block-level
  // element, which is the book convention.)

  // ── Headings ────────────────────────────────────────────────────────
  // h1 = chapter, h2 = section (N.M, inline-numbered), h3 = subsection
  // (no number, italic). Matches the on-website lecture pages: small
  // accent kicker above the chapter title, inline cobalt "N.M" prefix on
  // sections, plain italic subsection heads.
  show heading.where(level: 1): it => {
    pagebreak(weak: true)
    v(40pt)
    block(width: 100%)[
      #text(font: ("Inter", "Helvetica Neue", "Arial"),
            size: 11pt, weight: "semibold", fill: accent,
            tracking: 0.18em,
            upper[Chapter #counter(heading).display("1")])
    ]
    v(12pt)
    block(width: 100%)[
      #text(font: ("Source Serif 4", "Charter", "Georgia"),
            size: 30pt, weight: "medium", it.body)
    ]
    v(10pt)
    block(width: 60pt, height: 2pt, fill: accent)
    v(36pt)
  }

  show heading.where(level: 2): it => {
    v(32pt, weak: true)
    block(width: 100%)[
      #text(font: ("Source Serif 4", "Charter", "Georgia"),
            size: 18pt, weight: "medium")[
        #text(fill: accent)[#counter(heading).display("1.1")]
        #h(0.6em)
        #it.body
      ]
    ]
    v(14pt)
  }

  show heading.where(level: 3): it => {
    v(20pt, weak: true)
    block(width: 100%)[
      #text(font: ("Source Serif 4", "Charter", "Georgia"),
            size: 13pt, weight: "medium", style: "italic", it.body)
    ]
    v(8pt)
  }

  // ── Inline code ────────────────────────────────────────────────────
  show raw.where(block: false): it => {
    box(fill: code-bg, outset: (y: 2pt), inset: (x: 3pt), radius: 2pt,
        text(font: ("JetBrains Mono", "Menlo", "Monaco"), size: 0.92em, it))
  }

  // ── Code blocks ────────────────────────────────────────────────────
  show raw.where(block: true): it => {
    block(width: 100%, fill: code-bg, inset: 10pt, radius: 3pt,
          stroke: 0.5pt + rule,
          text(font: ("JetBrains Mono", "Menlo", "Monaco"), size: 9pt, it))
  }

  // ── Figures ────────────────────────────────────────────────────────
  show figure: set block(breakable: false)
  show figure.caption: it => {
    set text(size: 9.5pt, fill: muted, style: "italic", font: ("Inter", "Helvetica Neue", "Arial"))
    set par(first-line-indent: 0pt, justify: false, leading: 0.55em)
    pad(left: 2em, right: 2em)[
      #text(weight: "semibold", fill: accent)[Figure #it.counter.display("1.1")] \
      #it.body
    ]
  }

  // ── Lists ──────────────────────────────────────────────────────────
  set list(marker: ([•], [–]))
  set enum(numbering: "1.")

  // ── Links ──────────────────────────────────────────────────────────
  show link: it => text(fill: accent, it)

  body
}


// ─── Admonitions ─────────────────────────────────────────────────────
// Trade-paperback side-callouts. Four flavours: note / tip / warning /
// "why this matters" (book-specific). Each is a fenced box with a
// coloured left rule and a small uppercase label.

#let admonition(label, color, body) = {
  block(
    width: 100%,
    stroke: (left: 3pt + color),
    inset: (left: 14pt, right: 12pt, top: 10pt, bottom: 10pt),
    fill: color.lighten(92%),
    breakable: true,
  )[
    #text(font: ("Inter", "Helvetica Neue", "Arial"), size: 8.5pt, weight: "bold",
          fill: color, smallcaps(label))
    #v(4pt, weak: true)
    #set par(first-line-indent: 0pt, leading: 0.7em)
    #set text(size: 10pt)
    #body
  ]
}

#let note(body)    = admonition("Note",          rgb("#1e3a8a"), body)
#let tip(body)     = admonition("Tip",           rgb("#065f46"), body)
#let warn(body)    = admonition("Warning",       rgb("#b45309"), body)
#let danger(body)  = admonition("Caution",       rgb("#b91c1c"), body)
#let matters(body) = admonition("Why this matters", rgb("#6b21a8"), body)


// ─── Indexing ────────────────────────────────────────────────────────
// Mark a term for the back-of-book index. The marker is invisible —
// only its page location matters; the visible prose is unchanged. Use
// `sort:` to override the alphabetisation key (e.g. for "α-synuclein"
// pass sort: "alpha-synuclein").
#let idx(term, sort: none) = {
  let key = if sort == none { term } else { sort }
  [#metadata((term: term, sort: key))<idx>]
}

// Render the alphabetised back-of-book index. Groups entries by term,
// dedupes pages, lists in `term  pp. 12, 47, 102` form, two-column.
#let render-index() = context {
  let entries = query(<idx>)
  let groups = (:)
  for e in entries {
    let v = e.value
    let term = v.term
    let sort-key = v.sort
    let p = counter(page).at(e.location()).first()
    if term in groups {
      if not groups.at(term).pages.contains(p) {
        groups.at(term).pages.push(p)
      }
    } else {
      groups.insert(term, (term: term, sort: sort-key, pages: (p,)))
    }
  }
  let sorted = groups.values().sorted(key: g => lower(g.sort))

  // Render in two columns for compactness.
  set par(first-line-indent: 0pt, leading: 0.55em, justify: false)
  set text(size: 9.5pt)
  show: columns.with(2, gutter: 24pt)

  // Optional grouping by first letter — small header per letter.
  let current-letter = none
  for g in sorted {
    let l = upper(g.sort.at(0))
    if l != current-letter {
      current-letter = l
      v(8pt, weak: true)
      text(font: ("Inter", "Helvetica Neue", "Arial"),
           weight: "semibold", size: 10pt, fill: accent, l)
      v(2pt, weak: true)
    }
    let pages = g.pages.map(str).join(", ")
    block(below: 2pt)[
      #text(weight: "regular")[#g.term]  #h(0.4em)
      #text(fill: muted, size: 9pt)[#pages]
    ]
  }
}


// ─── Cover pages ─────────────────────────────────────────────────────
// Front and back cover, full-bleed A4, no running heads or folios.

#let book-cover(title, subtitle, author, edition, year) = {
  page(
    paper: "a4",
    margin: (top: 50mm, bottom: 30mm, x: 30mm),
    header: none,
    footer: none,
    background: none,
  )[
    #set align(center)
    #set par(first-line-indent: 0pt, justify: false, leading: 0.6em)

    // Hero mark
    #image("../figures/cover-mark.svg", height: 90mm)
    #v(16mm)

    // Kicker
    #block[
      #set text(font: ("Inter", "Helvetica Neue", "Arial"),
                size: 11pt, fill: subtle, tracking: 0.22em)
      #upper("A course in print")
    ]
    #v(14pt)

    // Title — supports multi-line content via line breaks in the
    // caller (e.g. [Bioinformatics for \ Engineers]).
    #block[
      #set text(font: ("Source Serif 4", "Charter", "Georgia"),
                size: 52pt, weight: "medium")
      #set par(leading: 0.32em)
      #title
    ]
    #v(12pt)

    // Cobalt rule
    #block(width: 80pt, height: 3pt, fill: accent)
    #v(16pt)

    // Subtitle
    #block[
      #set text(font: ("Source Serif 4", "Charter", "Georgia"),
                size: 22pt, style: "italic", fill: muted)
      #subtitle
    ]

    // Push author + edition to the bottom band
    #v(1fr)

    #block[
      #set text(font: ("Source Serif 4", "Charter", "Georgia"), size: 16pt)
      #author
    ]
    #v(10pt)
    #block[
      #set text(font: ("Inter", "Helvetica Neue", "Arial"),
                size: 10pt, fill: subtle, tracking: 0.18em)
      #upper(edition + "  ·  " + str(year))
    ]
  ]
}

#let book-back-cover(blurb, site-url, institution, semester, author) = {
  page(
    paper: "a4",
    margin: (top: 30mm, bottom: 30mm, x: 30mm),
    header: none,
    footer: none,
    background: none,
  )[
    #v(10mm)
    // A small mark in the top-right to echo the front cover
    #align(right)[
      #image("../figures/cover-mark.svg", height: 24mm)
    ]
    #v(12pt)
    // The blurb, in body serif
    #set par(justify: true, first-line-indent: 0pt, leading: 0.8em)
    #text(font: ("Source Serif 4", "Charter", "Georgia"), size: 11.5pt, blurb)
    #v(1fr)
    // Footer band — course meta + URL
    #line(length: 100%, stroke: 1pt + accent)
    #v(10pt)
    #set par(first-line-indent: 0pt, leading: 0.7em, justify: false)
    #set text(font: ("Inter", "Helvetica Neue", "Arial"), size: 9.5pt, fill: muted)
    #grid(
      columns: (1fr, auto),
      align: (left, right),
      [
        #text(weight: "semibold", fill: fg)[#author] \
        #institution \
        #semester
      ],
      [
        #text(tracking: 0.12em, upper("Live site"))\
        #text(fill: accent, size: 11pt, site-url)
      ],
    )
  ]
}


// ─── Frontmatter helpers ─────────────────────────────────────────────
#let book-title-page(title, subtitle, author, edition, year) = {
  align(center + horizon)[
    #set par(first-line-indent: 0pt, justify: false, leading: 0.6em)
    #v(-40pt)
    #block[
      #set text(font: ("Inter", "Helvetica Neue", "Arial"),
                size: 11pt, fill: subtle, tracking: 0.18em)
      #upper("A course in print")
    ]
    #v(12pt)
    #block[
      #set text(font: ("Source Serif 4", "Charter", "Georgia"),
                size: 40pt, weight: "medium")
      #set par(leading: 0.32em)
      #title
    ]
    #v(8pt)
    #block(width: 60pt, height: 2pt, fill: muted)
    #v(14pt)
    #block[
      #set text(font: ("Source Serif 4", "Charter", "Georgia"),
                size: 18pt, style: "italic", fill: muted)
      #subtitle
    ]
    #v(60pt)
    #block[
      #set text(font: ("Source Serif 4", "Charter", "Georgia"), size: 14pt)
      #author
    ]
    #v(80pt)
    #block[
      #set text(font: ("Inter", "Helvetica Neue", "Arial"),
                size: 10pt, fill: subtle, tracking: 0.18em)
      #upper(edition + "  ·  " + str(year))
    ]
  ]
  pagebreak()
}

#let copyright-page(year, author, license) = {
  set page(margin: (top: 30mm, bottom: 30mm, inside: 30mm, outside: 22mm))
  set text(size: 9pt, font: ("Inter", "Helvetica Neue", "Arial"), fill: muted)
  set par(justify: false, first-line-indent: 0pt, leading: 0.6em)
  v(1fr)
  [Copyright © #year #author. All rights reserved.\ ]
  v(8pt)
  [Released under #license. Distributed as part of the Bioinformatics for Engineers course at the School of Electrical Engineering, University of Belgrade.]
  v(20pt)
  [Set in Source Serif 4 (body), Inter (sans), and JetBrains Mono (code).]
  v(8pt)
  [Typeset with Typst.]
  v(1fr)
  pagebreak()
}
