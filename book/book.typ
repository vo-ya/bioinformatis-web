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

// Title is set as content (not a string) so the same explicit line
// break renders identically on the cover and the title page.
#let TITLE_CONTENT = [Bioinformatics for \ Engineers]

// ─── Cover ───────────────────────────────────────────────────────────
#book-cover(
  TITLE_CONTENT,
  "How biology became a data problem.",
  "Vojislav Varjačić",
  "First edition",
  2026,
)

// ─── Front matter ────────────────────────────────────────────────────
#book-title-page(
  TITLE_CONTENT,
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
#include "chapters/ch17-clinical-genomics.typ"
#include "chapters/ch18-cancer-genomics.typ"
#include "chapters/ch19-blast.typ"
#include "chapters/ch20-msa-phylogenetics.typ"
#include "chapters/ch21-hmms.typ"
#include "chapters/ch22-network-biology.typ"
#include "chapters/ch23-metagenomics.typ"
#include "chapters/ch24-crispr-screens.typ"
#include "chapters/ch25-causal-inference.typ"
#include "chapters/ch26-drug-discovery.typ"
#include "chapters/ch27-proteomics-metabolomics.typ"


// ─── Index ───────────────────────────────────────────────────────────
// Back-of-book alphabetised index. Markers are inserted into chapter
// prose via `book/index/tag.py`; the renderer below collects them.
#pagebreak()
#{
  set page(margin: (top: 30mm, bottom: 25mm, inside: 30mm, outside: 22mm))
  v(20pt)
  text(font: ("Source Serif 4", "Charter", "Georgia"),
       size: 28pt, weight: "medium")[Index]
  v(6pt)
  block(width: 60pt, height: 2pt, fill: rgb("#1e3a8a"))
  v(20pt)
  render-index()
}

// ─── Back cover ──────────────────────────────────────────────────────
#book-back-cover(
  [
    Twenty-seven chapters of bioinformatics taught from an engineer's
    point of view, expanded from the live lecture series at the
    School of Electrical Engineering, University of Belgrade. The
    book starts at the cell and the FASTQ file and walks all the way
    out to mass spectrometry, drug discovery, and the AI methods now
    rewriting structural biology — pausing along the way for the
    algorithms (BWT, EM, Viterbi, MR), the statistical machinery
    (NB GLMs, Bayesian genotyping, BH correction), and the
    instruments (Illumina, Nanopore, Orbitrap) that make modern
    genomics possible. Each chapter includes worked exercises, a
    further-reading list, and a companion Colab notebook on the live
    site.

    Written for engineers who never took biology and biologists who
    want the algorithms drawn properly.
  ],
  "vo-ya.github.io/bioinformatis-web",
  "School of Electrical Engineering, University of Belgrade",
  "Spring 2026",
  "Vojislav Varjačić",
)
