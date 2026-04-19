# Lecture Authoring Style Guide

> **Purpose**: Reference this file when drafting any new lecture in this bioinformatics course.
> **How to use in a new chat**: paste this file as context at the start, along with your lecture topic and target duration. Tell the assistant: "Draft the lecture following this style guide."
> **Companion files**: `website-spec.md` (site/visual design), `artifacts-spec.md` (interactive artifact conventions).

---

## 1. Course Identity — What Makes This Course Different

This is a **bioinformatics course taught at a School of Electrical Engineering**. Students have an EE background (signal processing, information theory, control systems, instrumentation, ML) and minimal biology background. The course's identity rests on one core move:

> **Translate biology into things an EE student already understands.**

Not as a gimmick — as the pedagogical spine. Wherever an EE concept genuinely maps to a biological concept, make the map explicit. Where it doesn't, don't force it.

Canonical examples from Lecture 1:
- **Genetic code redundancy** ↔ Hamming distance / error correction.
- **PCR** ↔ positive-feedback amplifier with fixed gain per cycle, noise, saturation regime.
- **Phred quality scores** ↔ dB — "Q30 = 30 dB of base-call confidence."
- **Nanopore basecalling** ↔ sequence-to-sequence ML, CTC, speech recognition.
- **Illumina base calling** ↔ multi-channel signal detection with crosstalk matrix.

Each of these is a real, nontrivial isomorphism. The course earns its identity by finding more of them.

---

## 2. Voice & Tone

**Target voice**: a competent senior engineer explaining something genuinely interesting to a sharp junior engineer. Assumes intelligence. Doesn't dumb down. Doesn't show off.

### Do

- Write in **connected prose**. Paragraphs, not bullet dumps. This is a scrollable lecture, not a slide deck — see §6 for more.
- Be **direct**. "A gene is a stretch of DNA that codes for a functional product" — not "In this section we will explore the concept of what a gene is."
- Use **specific numbers**. "Humans have ~20,000 protein-coding genes" is better than "humans have many genes."
- Use **named people, places, and years** when giving history. "Kary Mullis, 1983" is better than "in the 1980s, a scientist."
- Name **exceptions** honestly. "The genetic code is *nearly* universal — a few exceptions exist in mitochondria and some protists" is better than pretending the rule is perfect.
- Use **em-dashes** for asides when they add rhythm. "Taq polymerase — from Thermus aquaticus, a hot-spring bacterium — is essential because it survives 95°C."
- Admit **uncertainty and open questions**. "Much of the genome's non-coding function is still being mapped."

### Don't

- No hedging AI-voice ("As an AI...", "It's important to note that...", "In the world of biology...").
- No marketing-style exclamation. This is not a textbook trying to seem exciting.
- No padded transitions. "Now that we've covered X, let's move on to Y" — just start Y.
- No condescending simplifications. "Think of DNA as a cookbook" is fine *once* and only if it earns its place; don't keep cooking the metaphor.
- No gratuitous jargon either. First use of every technical term gets a brief definition, a bold on first mention, and ideally the intuition before the term.

### Sentence-level patterns

- Prefer active voice. "RNA polymerase reads the template strand" > "the template strand is read by RNA polymerase."
- Short sentences for emphasis, longer sentences for nuance. Vary.
- Colons before lists, semicolons between related clauses. Respect them.
- Fragments are OK for punch — sparingly. "Substitution errors dominate. Indels are rare."

---

## 3. Document Structure

Every lecture follows this structure. Lengths scale with the topic; the skeleton does not.

```
# Lecture N — Title: Subtitle (optional)

> **Duration**: ≈X min content + breaks
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-NN.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1. [verb + specific outcome]
2. [verb + specific outcome]
3. ... (5–8 objectives total)

---

## Part 1 — [Title] (≈X min)

### 1.1 [Subsection] (≈X min)

[prose...]

### 1.2 [Subsection] (≈X min)

[prose...]

**EMBED — Artifact #N: [Name]** → `artifacts/lecture-NN/NN-name.html`
*[One-sentence description of what the artifact shows.]*

## Part 2 — [Title] (≈X min)

...

---

### ☕ Break (10 min)

---

## Part 3 — [Title] (≈X min)

...

---

## Wrap-up (≈10 min)

### What you should take away

- [bullet, 4–6 items, crisply stated]

### Next lecture

- [1–2 sentence pointer to the following lecture's topic]

### Homework

1. [concrete, specific assignment — something they do, not just read]
2. ...

### Recommended reading

- [Book / chapter]
- [Classic paper, with author, year, title, journal, volume, pages]
- [Another]
- [URL for primary data if relevant]

---

## Appendix — Speaker cues and timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — X | 45 min | 0:45 |
| ... | ... | ... |

**Total**: ~Xh Ymin.
```

### Rules for the skeleton

- **Top matter is a blockquote block** with Duration, Audience, File. Always those three fields, always in that order.
- **Learning Objectives** are a numbered list of 5–8 items. Each starts with a verb (describe, explain, compare, derive, compute, interpret…). Specific outcomes, not vague ones.
- **Parts** are numbered (`Part 1`, `Part 2`, …). Each has an estimated duration in the header.
- **Subsections** are numbered hierarchically (`1.1`, `1.2`, `3.4`, …). Each has its own duration estimate.
- **Breaks** are marked with `☕` and italic `### ☕ Break (N min)` between `---` separators. Typically one break per ~90 min of content.
- **Artifact embeds** use the exact marker shown above — this is how they're located and wired into the HTML later.
- **Wrap-up** always has: takeaways, next-lecture pointer, homework, references. All four.
- **Timing appendix** at the bottom — a table of cumulative time. The lecturer checks this in real time during class.

---

## 4. Callouts — The Only Approved Types

Callouts punctuate prose with a specific rhetorical move. **Never use callouts as speaker notes, and never invent new types.** The site's CSS only styles the ones below.

**All callouts use blockquote syntax: `> **Label**: content`.**

### 4.1 `Intuition box` — analogies and mental models

Use to translate a concept into a student-native mental model. Often EE-flavored.

> **Intuition box**: Think of the cell as a computer. DNA is persistent storage (hard drive). RNA is RAM / message-passing. Proteins are the running programs — they actually do things. The ribosome is the CPU. Mitochondria are the power supply unit, with a surprising amount of its own firmware.

**When to use**: once per major concept, right after the concept is introduced and explained straight. Don't use for every paragraph; lose impact.

### 4.2 `EE framing` — explicit isomorphism with an EE concept

Use when the map to an EE concept is precise enough to be taken seriously, not just metaphorical.

> **EE framing**: PCR is a positive-feedback amplifier with a fixed gain per cycle. Like any amplifier, it has noise (polymerase errors, primer-dimers), a saturation regime (reagents deplete, plateau), and dynamic range limits.

**When to use**: when the analogy carries real explanatory weight — i.e. an EE student can use their existing intuition to predict behavior of the biological system. If you're reaching, leave it out.

**Variant**: `EE framing (dB)`, `EE framing, push harder`, `EE framing — error correction` — subtitles are allowed to point at the specific EE concept being invoked.

### 4.3 `Historical pointer` — context, without disrupting the argument

Use for one-paragraph historical asides that add color and narrative arc.

> **Historical pointer**: Predicting protein structure from amino acid sequence was a grand challenge for 50 years. AlphaFold (2021) largely solved it using transformer-based deep learning. If any single application convinces EE students this field is theirs, it's this.

**When to use**: sparingly — maybe 2–4 times per lecture. Keep them tight.

### 4.4 `Discussion prompt` — questions that should stay open

Use to invite students to think, not to quiz them. Prompts should have interesting answers, not obvious ones.

> **Discussion prompt**: If 98% of your genome doesn't code for protein, what could it possibly be doing? (Answers: regulation, structural chromatin roles, repeat elements, ancient viral DNA, genuinely non-functional sequence…)

**When to use**: once or twice per major part. Optional to include suggested answers in parentheses — useful for the scrollable format because the student is reading alone.

### 4.5 `Warning box` — pitfalls and common mistakes

Use for real gotchas that will trip people up. Not for decoration.

> **Warning box**: Always know which encoding you're dealing with. Mixing them silently produces quality scores that are off by 31 — disastrous and non-obvious. Tools like FastQC auto-detect and warn.

**When to use**: when the consequence of misunderstanding is a real mistake in practice (mis-encoding, mis-interpreted axis, units confusion, file-format gotcha, etc.).

### 4.6 Metadata blockquotes — Duration / Audience / File

Only at the top of the file, in the document header. Not used elsewhere.

### 4.7 `This is the deepest technical block…` — section emphasis

Used sparingly to flag a block that deserves extra attention.

> **This is the deepest technical block of the lecture.** Illumina is ~80% of sequencing output globally — if students understand one platform, make it this one.

**When to use**: max once per lecture, if at all.

### 4.8 ❌ What is NEVER a callout type

- **No Speaker notes.** The lecture is written for the student who will read it directly; the lecturer teaches from the same text projected. There is no hidden authoring layer. If content helps the student, put it in prose; if it only helps the lecturer, delete it.
- No "Fun fact", "Did you know?", "Pro tip" — these collapse the voice into pop-science.
- No "Summary" callouts at the end of sections. If a section needs a summary, the next paragraph can open with one sentence of it; the wrap-up takeaways cover the lecture as a whole.

---

## 5. Writing Patterns

### 5.1 Introducing a technical term

Three-step pattern: **intuition → bolded term → concrete detail or example**.

Good:
> The reaction produces a population of fragments of every possible length, each ending in a known fluorescently labeled base. **Capillary electrophoresis** separates fragments by length with single-nucleotide resolution.

Bad:
> Capillary electrophoresis, which is a technique used in molecular biology, separates fragments...

### 5.2 Lists inside prose

If you have 3+ items that are each ≤1 sentence and parallel, use a bulleted list. If they are each 2+ sentences, use a numbered list with labels.

**Bulleted list**: short, parallel.
```markdown
Functional RNA types worth naming:
- **mRNA** — messenger RNA, the transcript that gets translated into protein.
- **tRNA** — transfer RNA, brings amino acids to the ribosome.
- **rRNA** — ribosomal RNA, structural/catalytic component of ribosomes.
```

**Numbered list**: steps, or long items where numbering aids reference.
```markdown
Thermal cycle (repeated ~30 times):
1. **Denature** (~95°C): strands separate.
2. **Anneal** (~55–65°C): primers bind to their complementary sites.
3. **Extend** (~72°C): polymerase synthesizes new strand.
```

### 5.3 Tables

Use tables for genuinely comparative data — columns must be meaningfully parallel. Keep them compact; no decorative columns.

Good (from Lecture 1):

```markdown
| Platform | Read length | Raw error | Consensus error | Error type |
|---|---|---|---|---|
| Sanger | ~800 bp | <0.1% | — | Random, low |
| Illumina | 150–300 bp | ~0.1–1% | — | Substitutions, quality drops late in read |
| PacBio HiFi (CCS) | 10–25 kb | — | <0.1% | Effectively Sanger-like, at scale |
```

Don't use tables for what would read better as prose. "A gene is made of exons, introns, and UTRs" is prose, not a 3-row table.

### 5.4 Equations

Use KaTeX-compatible LaTeX. Inline as `$Q = -10 \log_{10}(P)$`, display as:

```markdown
```
Q = −10 · log₁₀(P_error)
```
```

For the markdown source file, ASCII-art equations like `Q = −10 · log₁₀(P_error)` are fine — they'll be re-rendered in KaTeX during the HTML build. Keep display equations in fenced code blocks in the markdown so they survive formatting.

### 5.5 ASCII diagrams

Use freely for anything that has a simple schematic structure — pipelines, flow diagrams, sequence layouts. Keep them inside fenced code blocks so indentation survives.

Good (from Lecture 1):

```
5'─[P5 adapter]─[i5 index]─[Read1 primer]─[INSERT]─[Read2 primer]─[i7 index]─[P7 adapter]─3'
```

### 5.6 Specific numbers with `~` for approximations

Always include specific numbers when you have them. Prefix with `~` when approximate. "~1% error" not "a low error rate"; "~20,000 genes" not "many genes".

### 5.7 Bold, italic, monospace

- **Bold** — first mention of a technical term, or emphasis on a consequential word. Not for decoration.
- *Italic* — book titles, or gentle stress. Sparingly.
- `Monospace` — anything that is a literal string: DNA sequences, file names, commands, format names (e.g., `FASTQ`), tool names when referring to the binary (`bwa`, `samtools`).

### 5.8 Numbers and units

- Use thin-style formatting: `3.1 Gb`, `~150 bp`, `2×150 bp`, `$1,000`, `99.9%`.
- Use `%` symbol, not spelled out.
- Scientific notation with `×`: `2³⁰ ≈ 10⁹×`.
- Ranges with en-dash: `10–25 kb`, `55–65°C`.

---

## 6. Scrollable Format — Implications for Writing

Since the lecture is a single long scrollable page (not slides — see `website-spec.md` §2), the writing must support continuous reading, not chunked presentation.

- **No "slide-like" bullet dumps.** A section of 8 one-liner bullets is a failure mode — convert to prose or to a genuinely structured list.
- **No "As we saw on the previous slide…"** — the reader can just scroll up.
- **No per-section summaries** — the takeaways section at the end covers that.
- **Use headings generously** — since the lecture is long, a clear heading hierarchy (H2 for parts, H3 for subsections, H4 if truly needed) gives the reader anchor points. The site auto-builds a TOC from them.
- **Keep paragraphs to 3–6 sentences** — walls of text are hard to read in long form.

---

## 7. Artifacts Integration

When a lecture has interactive artifacts, mark them inline with this exact format:

```markdown
**EMBED — Artifact #N: [Name]** → `artifacts/lecture-NN/NN-name-kebab.html`
*[One-sentence italic description of what the artifact shows.]*
```

Rules:
- **Embed markers go where the concept lands**, immediately after the prose that introduces the concept — not at the end of the section.
- **Number artifacts sequentially** within a lecture (`#1`, `#2`, …).
- **File names** follow `NN-name-kebab.html` with zero-padding.
- **Path** matches the folder structure: `artifacts/lecture-NN/`.

The italic description line is read by students (when the artifact is lazy-loading) and by the developer building the artifact. Keep it concrete.

Artifacts are specified in detail in a separate `artifacts-spec.md` file per lecture. Writing the lecture doesn't require writing the artifact spec — they are separate documents, iterated on separately.

---

## 8. Homework

Homework assignments are mandatory in the wrap-up. They must be:

- **Concrete**: something the student does (writes, computes, plots, derives) — not "read chapter 5".
- **Bounded**: completable in a few hours with the lecture as context.
- **Aligned with EE skills**: scripting in Python or MATLAB, plotting, implementing a small algorithm, analyzing a file format.
- **Answerable**: have a right answer (or a defensible one), not open-ended "reflect on…"

Good (from Lecture 1):

> 1. Download a small FASTQ file.
> 2. Write a Python or MATLAB script that parses it without using Biopython, computes mean quality per read, and plots mean quality as a function of position along the read.
> 3. Answer: at what read position does quality start to drop off, and why?

Bad: "Read about sequencing platforms."

---

## 9. References / Recommended Reading

Always include a short list — 3–6 items. Mix of:

- A **textbook chapter** (e.g., Alberts et al., Molecular Biology of the Cell, Chapter N).
- A **review article** with full citation.
- **Primary sources** when historically important.
- A **URL** for primary data or a tool's documentation if relevant.

Format:

```markdown
- Surname, A. B., & Surname, C. D. (Year). Title of paper. Journal Name volume, pages.
- Book author(s). Book title. (Chapter N)
- https://example.gov/data-page
```

Keep it tight. This is a reading list, not a bibliography.

---

## 10. Prohibited Patterns

Quick list of patterns to avoid, collected from the above:

- ❌ Speaker notes of any kind.
- ❌ "As an AI..." / "I think..." / "It's important to note..."
- ❌ Slide-style bullet dumps without prose connective tissue.
- ❌ Invented callout types ("Fun fact", "Did you know", "Pro tip").
- ❌ Empty transitions ("Now let's move on to...", "In this next section we will explore...").
- ❌ Vague numbers when specific ones are available.
- ❌ "Think of X as Y" metaphors that aren't followed up with the actual mechanism.
- ❌ Tables that would read better as prose.
- ❌ Mid-section summaries.
- ❌ Absolute URLs to the site or artifact files — use relative paths.
- ❌ Calling things "cutting-edge", "revolutionary", "state-of-the-art" without specifics.

---

## 11. Pre-Submission Checklist

Before considering a lecture draft done, verify:

- [ ] Top matter present: Duration, Audience, File — in that order, in a blockquote.
- [ ] Learning Objectives: 5–8 numbered items, each starts with a verb.
- [ ] Parts numbered `Part 1`, `Part 2`, …, each with duration estimate.
- [ ] Subsections numbered hierarchically (1.1, 1.2, …), each with duration.
- [ ] At least one break (`☕ Break (10 min)`) per ~90 min of content.
- [ ] No speaker notes anywhere.
- [ ] All callouts use approved labels only (§4).
- [ ] Artifact embeds use the exact marker format (§7).
- [ ] At least one EE framing callout per major part — or a clear reason why not.
- [ ] Wrap-up has all four: Takeaways, Next lecture, Homework, References.
- [ ] Timing appendix table at the bottom with cumulative times.
- [ ] Total timing matches the course's per-lecture budget (3–4 hours).
- [ ] No absolute URLs.
- [ ] No slide-style bullet dumps replacing prose.
- [ ] Specific numbers with sources where appropriate.

---

## 12. Quick-Start Template

Copy this into a new file and fill it in.

```markdown
# Lecture N — Title

> **Duration**: ≈X min content + breaks
> **Audience**: EE undergraduates / graduates, minimal biology background assumed
> **File**: to be rendered as `lectures/lecture-NN.html`

---

## Learning Objectives

By the end of this lecture, a student should be able to:

1.
2.
3.
4.
5.

---

## Part 1 — [Title] (≈X min)

### 1.1 [Subsection] (≈X min)

[prose]

### 1.2 [Subsection] (≈X min)

[prose]

**EMBED — Artifact #1: [Name]** → `artifacts/lecture-NN/01-name.html`
*[description]*

## Part 2 — [Title] (≈X min)

...

---

### ☕ Break (10 min)

---

## Part 3 — [Title] (≈X min)

...

---

## Wrap-up (≈10 min)

### What you should take away

-

### Next lecture

-

### Homework

1.

### Recommended reading

-

---

## Appendix — Speaker cues and timing summary

| Block | Time | Cumulative |
|---|---|---|
| Part 1 — | X min | 0:X |

**Total**: ~Xh Ymin.
```

---

## 13. How to Prompt in a New Chat

Recommended opening prompt when drafting Lecture N in a fresh chat:

> I'm drafting Lecture N of a bioinformatics course taught at a School of Electrical Engineering. Lectures are 3–4 hours. The topic is [TOPIC]. Draft the lecture following the attached style guide exactly — structure, callout types, voice, EE framings, scrollable format. No speaker notes. Include interactive artifact embed markers where appropriate.
>
> [Attach this style guide + the topic outline]
>
> Start with a short layout (parts + subsection headers + timing) for my review before writing the full prose.

The "layout first, then full draft" step is important. It matches how Lecture 1 was built and catches structural problems before prose is written.
