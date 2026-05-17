// Bioinformatics for Engineers — print edition.
// Builds: typst compile book.typ build/book.pdf
//
// To compile a single chapter for review, use chapter-only.typ.

#import "theme/book-theme.typ": *

#show: book-doc.with(
  title: "Bioinformatics for Engineers",
  subtitle: "How biology became a data problem.",
  author: "Vojislav Varjačić",
)

// ─── Front matter ────────────────────────────────────────────────────
#book-title-page(
  "Bioinformatics for Engineers",
  "How biology became a data problem.",
  "Vojislav Varjačić",
  "First edition",
  2026,
)
#copyright-page(2026, "Vojislav Varjačić", "the Apache License 2.0")

// ─── Table of contents ───────────────────────────────────────────────
#{
  set page(margin: (top: 30mm, bottom: 25mm, inside: 30mm, outside: 22mm))
  v(20pt)
  text(font: ("Source Serif 4", "Charter", "Georgia"), size: 24pt, weight: "medium")[Contents]
  v(8pt)
  block(width: 60pt, height: 1pt, fill: rgb("#374151"))
  v(20pt)
  outline(title: none, depth: 2, indent: auto)
  pagebreak()
}

// ─── Chapters ────────────────────────────────────────────────────────
#counter(page).update(1)
#counter(heading).update(0)

#include "chapters/ch01-foundations.typ"
#include "chapters/ch02-read-alignment.typ"
#include "chapters/ch03-assembly.typ"
#include "chapters/ch04-variant-calling.typ"
#include "chapters/ch05-bulk-rna-seq.typ"
#include "chapters/ch06-differential-expression.typ"
#include "chapters/ch07-scrna-seq.typ"
#include "chapters/ch08-advanced-single-cell.typ"
#include "chapters/ch09-chip-atac-peak-calling.typ"
#include "chapters/ch10-methylation-hic.typ"
#include "chapters/ch11-long-reads-pangenome.typ"
#include "chapters/ch12-population-genetics.typ"
#include "chapters/ch13-gwas.typ"
#include "chapters/ch14-data-engineering.typ"
#include "chapters/ch15-protein-structure.typ"
#include "chapters/ch16-ml-in-genomics.typ"

// Future chapters get included here as they are authored.
// #include "chapters/ch17-clinical-genomics.typ"
// ...
