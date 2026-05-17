// ─── Book theme ──────────────────────────────────────────────────────
// Trade-paperback technical-book look for the Bioinformatics for
// Engineers print edition. A4 trim, serif body, sans heads, mono code,
// generous margins, mirrored running heads.

#let accent = rgb("#1e3a8a")           // deep cobalt (matches the website)
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
  set par(leading: 0.75em, justify: true, first-line-indent: 1.2em)

  // Auto-number headings; show rules below choose where each level's
  // number appears (kicker for h1, inline cobalt prefix for h2, hidden
  // for h3).
  set heading(numbering: "1.1")

  // No first-line indent right after a heading, blockquote, figure, code, or list.
  show heading: it => { set par(first-line-indent: 0pt); it }
  show figure: it => { set par(first-line-indent: 0pt); it }
  show raw.where(block: true): it => { set par(first-line-indent: 0pt); it }
  show quote: it => { set par(first-line-indent: 0pt); it }

  // ── Headings ────────────────────────────────────────────────────────
  // h1 = chapter, h2 = section (N.M, inline-numbered), h3 = subsection
  // (no number, italic). Matches the on-website lecture pages: small
  // accent kicker above the chapter title, inline cobalt "N.M" prefix on
  // sections, plain italic subsection heads.
  show heading.where(level: 1): it => {
    pagebreak(weak: true)
    v(36pt)
    text(font: ("Inter", "Helvetica Neue", "Arial"),
         size: 11pt, weight: "semibold", fill: accent,
         tracking: 0.18em,
         upper[Chapter #counter(heading).display("1")])
    v(10pt)
    text(font: ("Source Serif 4", "Charter", "Georgia"),
         size: 30pt, weight: "medium", it.body)
    v(8pt)
    block(width: 60pt, height: 2pt, fill: accent)
    v(28pt)
  }

  show heading.where(level: 2): it => {
    v(22pt, weak: true)
    block[
      #text(font: ("Source Serif 4", "Charter", "Georgia"),
            size: 18pt, weight: "medium")[
        #text(fill: accent)[#counter(heading).display("1.1")]
        #h(0.6em)
        #it.body
      ]
    ]
    v(4pt, weak: true)
  }

  show heading.where(level: 3): it => {
    v(14pt, weak: true)
    text(font: ("Source Serif 4", "Charter", "Georgia"),
         size: 13pt, weight: "medium", style: "italic", it.body)
    v(2pt, weak: true)
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
    set par(first-line-indent: 0pt, justify: false)
    pad(left: 2em, right: 2em)[
      *Figure #it.counter.display("1.1")*  #h(0.5em) #it.body
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


// ─── Frontmatter helpers ─────────────────────────────────────────────
#let book-title-page(title, subtitle, author, edition, year) = {
  align(center + horizon)[
    #v(-40pt)
    #text(font: ("Inter", "Helvetica Neue", "Arial"), size: 11pt, fill: subtle, smallcaps("a course in print"))
    #v(12pt)
    #text(font: ("Source Serif 4", "Charter", "Georgia"), size: 44pt, weight: "medium", title)
    #v(4pt)
    #block(width: 60pt, height: 2pt, fill: muted)
    #v(12pt)
    #text(font: ("Source Serif 4", "Charter", "Georgia"), size: 18pt, style: "italic", fill: muted, subtitle)
    #v(60pt)
    #text(font: ("Source Serif 4", "Charter", "Georgia"), size: 14pt, author)
    #v(80pt)
    #text(font: ("Inter", "Helvetica Neue", "Arial"), size: 10pt, fill: subtle, edition + "  ·  " + str(year))
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
