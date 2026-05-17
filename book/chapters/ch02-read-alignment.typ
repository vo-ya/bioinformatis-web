#import "../theme/book-theme.typ": *

= #idx("read alignment")Read Alignment: From Brute Force to #idx("FM index")FM Index and Back <ch:read-alignment>

#matters[
  Every other chapter in this book starts with a #idx("FASTQ")FASTQ file produced by a
  sequencer. Almost no analysis runs on that file as-is. #idx("variant calling")Variant calling,
  expression quantification, #idx("ChIP-seq")ChIP-seq peak finding, #idx("methylation")methylation analysis,
  structural-variant detection — all of them assume that someone has
  already placed each read onto the genome it came from, with an
  approximate edit distance and a coordinate. That placement is _read
  alignment_, and it is the single most-CPU-burning step in modern
  genomics. Get it right, and the rest of the stack works. Get it wrong,
  and every downstream tool fails in ways that are hard to debug because
  the bug is two pipeline steps upstream. The algorithms in this chapter
  are also some of the most elegant in computer science — invented
  partly for compression, partly for #idx("STRING")string matching, partly for genome
  scale — and they reward studying them as algorithms, not just as
  bioinformatics tools.
]

A sequencing run produces a few hundred million short strings of A, C,
G, and T. Lecture 1 ended there: that is the raw output of the
instrument, and by itself it is almost useless. A FASTQ file is an
unordered pile of 150-base-pair fragments. Nothing biological lives in
that pile; biology lives in the positions on a #idx("chromosome")chromosome that the
fragments came from. Putting them back is the substrate everything
downstream runs on.

The naive version of the problem is one of the oldest in computer
science: given a long text and a short pattern, find all occurrences of
the pattern. Knuth, Morris, and Pratt published their linear-time
algorithm for this in 1977; Boyer and Moore published theirs the same
year. By 1990 a few decades of string-matching theory had accumulated
and the problem looked finished. Then sequencing arrived and the text
got bigger. By 2008 the text was three billion characters and there
were a billion patterns to match against it, each pattern allowed to
differ from the text at a few positions. The classical algorithms broke.
The ones that replaced them are the subject of this chapter.

There is a clean four-move plan in what follows. _Index the reference._
_Shrink the index._ _Score approximate matches._ _Compose the two._ The
first three are individually beautiful. The fourth — composing a fast
exact index with a slow accurate scorer — is the design pattern that
every production aligner is built on. It also turns out to be the same
pattern radar and GPS receivers use, which is not a coincidence.


== The Alignment Problem <sec:problem>

The problem, stated cleanly: given a reference string $R$ over the
alphabet ${A, C, G, T}$ and a read $r$, find all positions $i$ such
that $R[i .. i + |r| - 1]$ matches $r$ within some edit distance $k$.
For the human reference, $|R| approx 3 times 10^9$; for an #idx("Illumina")Illumina
read, $|r| approx 150$; for a typical whole-genome run, you align
about $10^9$ reads against the same reference.

Do the arithmetic once, out loud. Each read has roughly $3 times 10^9$
candidate alignment positions in the reference. At each candidate
position you may need to compare up to 150 characters. That gives
$3 times 10^(11)$ candidate-position-by-character comparisons per read.
Multiply by $10^9$ reads and you have $3 times 10^(20)$ comparisons per
dataset. Even at a billion comparisons per second per CPU core — which
is optimistic — that is roughly ten thousand CPU-years per genome.
Nobody has ten thousand CPU-years per genome.

#figure(
  image("../../diagrams/lecture-02/01-alignment-problem-overview.svg", width: 95%),
  caption: [
    A reference genome with a handful of reads placed against it.
    Exact match, mismatched read, read with a small insertion.
    Alignment is the problem of finding, for every read, the position
    it came from.
  ],
) <fig:overview>

The scale is what forces everything that follows. If the genome were a
kilobyte, brute force would be fine and this chapter would be five
pages long. It is not, so the algorithms used in practice are ones that
were largely invented — or adapted — for exactly this problem. They are
beautiful, and they are also the only way the modern field could exist.

Two complications separate read alignment from plain string matching.
First, reads contain errors. Illumina sequencing has a per-base error
rate of about 0.1 to 1 per cent, dominated by substitutions with a
small tail of indels. A 150-base read therefore has, on average, one
or two errors. Long-read platforms have higher per-base error rates,
sometimes much higher. Second, the genome you sequenced differs from
the reference at about one position in a thousand because of real
biological variation. Together, these two effects make exact-match the
wrong target. The right target is "best approximate match," under some
scoring scheme that quantifies how much each kind of difference costs.

#note[
  The asymmetry between the reference and the reads is everything in
  what follows. The reference is fixed: you index it once, perhaps for
  hours, and then keep the index on disk. The reads are different
  every run, and you have a billion of them. Any work you can move
  from the per-read side to the once-per-reference side is a good
  investment. Any work you do per read had better be cheap.
]

Every production aligner exploits both facts the same way. They
_index_ the reference once, producing a data structure that supports
fast exact lookup of short substrings — twenty to thirty bases,
typically. They use that index to find candidate positions per read —
called _seeds_ — and then run a slower, accurate, error-tolerant
algorithm on small windows around each seed to produce the final
alignment with its errors and indels. This two-phase structure — fast
coarse detection, slow fine estimation — is the _seed-and-extend_
pattern. It is also exactly the same decomposition as coarse
acquisition and fine tracking in a GPS receiver: a cheap correlation
against a short code narrows the candidate delays down to a handful,
and then an expensive tracking loop runs only on those candidates.
Nobody runs the tracking loop over the full search space, because
nobody has the compute. The same is true here.

This chapter is organised by that split. Sections 2.2 and 2.3 build the
fast exact index — first naively, then with suffix arrays, hash tables,
and finally the Burrows-Wheeler / FM index that real aligners actually
use. Section 2.4 builds the slow accurate extension algorithm,
#idx("Smith-Waterman")Smith-Waterman, and shows how its output is serialised into the #idx("CIGAR")CIGAR
string stored in #idx("BAM")BAM records. Section 2.5 composes the two into the
pipeline every aligner — `bwa`, `bowtie2`, `minimap2` — implements.


== Brute Force, Suffix Arrays, and Hash Indices <sec:exact-indices>

We start with the algorithms that no production aligner uses, because
every production aligner has to be measured against them. The unit of
account throughout is the _character comparison_. Every clever index
that follows is, fundamentally, a way to avoid doing them.

=== Brute force

The naive algorithm for finding a read in a reference is to slide the
read along the reference one position at a time and, at each position,
check every character for a match. It is correct, easy to write, and
hopeless at genome scale.

```python
def naive_search(R, r):
    for i in range(len(R) - len(r) + 1):
        if R[i : i + len(r)] == r:
            yield i
```

#figure(
  image("../../diagrams/lecture-02/02-brute-force-sliding.svg", width: 95%),
  caption: [
    Brute force, visualised. At every position the read is compared
    character-by-character against the reference. Matches are dots;
    mismatches are crosses. The algorithm runs across the whole genome,
    one position at a time.
  ],
) <fig:brute-force>

The worst-case cost is $O(|R| dot |r|)$. Real implementations exit
early on mismatch, which helps on random text but is almost worthless
on a genome that is roughly half repetitive sequence. Boyer-Moore and
Knuth-Morris-Pratt, the famous 1977 string-matching speed-ups, do
better — Boyer-Moore by sliding the read in jumps larger than one when
the mismatch tells you it is safe to skip, KMP by pre-computing a
failure table that avoids recomparing characters that already matched.
Both are still $Omega(|R| + |r|)$ per query in their best forms.
Multiplied by $10^9$ reads, they are still nowhere near tractable.

The reason to start here is not historical. It establishes the unit of
account. Every index in the rest of this chapter is a strategy for
doing fewer comparisons, and the count gives you a clean way to
compare them.

=== Suffix arrays

A _suffix_ of a string $R$ is a substring that ends at the last
character. For `R = BANANA$`, the suffixes are `BANANA$`, `ANANA$`,
`NANA$`, `ANA$`, `NA$`, `A$`, and `$`. The trailing `$` is a sentinel
character that compares less than every real character; its purpose is
to force every suffix to be unique and to give every lexicographic
comparison a clean place to terminate.

Here is the first idea that pays for itself. Sort the suffixes
lexicographically. The result is a list of all positions in $R$,
reordered so that if you look at the substring starting at each
position, the substrings are in dictionary order. The _suffix array_
SA stores, for each rank in the sorted order, the position in $R$ where
that suffix begins.

#figure(
  image("../../diagrams/lecture-02/03-suffix-array-structure.svg", width: 95%),
  caption: [
    The #idx("suffix array")suffix array of `BANANA$`. Seven suffixes, listed in sorted
    order. The rows whose suffixes start with `ANA` form a contiguous
    range — which is the reason binary search works on this structure.
  ],
) <fig:suffix-array>

Now the second idea. A read $r$ occurs at position $i$ in $R$ if and
only if $r$ is a prefix of the suffix of $R$ starting at $i$. In the
sorted list, all suffixes that share a prefix are contiguous; they
form a block. So finding all occurrences of $r$ in $R$ reduces to
finding the block of suffixes that begin with $r$. Two binary searches
locate the block's left and right boundaries, each performing at most
$log_2 |R|$ comparisons, each comparison examining at most $|r|$
characters. The total cost per query is $O(|r| log |R|)$.

Put numbers on it. For the human genome, $log_2 (3 times 10^9) approx
32$. With $|r| = 150$, that is roughly 4,800 character comparisons per
read — against the $3 times 10^(11)$ of brute force. Almost eight
orders of magnitude, for the cost of building the suffix array once.

#note[
  The suffix array as a data structure was introduced by Udi Manber
  and Gene Myers in 1990, explicitly pitched as a space-saving
  alternative to the suffix tree, which had dominated stringology
  since Peter Weiner's 1973 paper. The suffix tree is asymptotically
  elegant — linear-time build, linear-time queries — but carries a
  large constant factor in memory. Manber and Myers's move was to
  ask whether you could give up the tree structure, keep only the
  sorted order of the leaves, and still get fast queries. With
  binary search, the answer was yes. Every subsequent genome-scale
  indexing scheme has descended from that reduction.
]

The suffix array is fast at query time but expensive in memory.
Naively, it stores $|R|$ integers of position data. For the human
genome at 3 Gb with four-byte integers, that is 12 GB — before you
store $R$ itself. Sparsified and compressed versions exist, but the
uncompressed form already crosses the line of comfortable.

Construction is not free either. A naive sort of suffixes is
$O(|R|^2 log |R|)$ in the worst case because string comparisons are
themselves $O(|R|)$. The classical linear-time algorithms — DC3,
SA-IS — bring construction to $O(|R|)$, which for the human genome
takes minutes to tens of minutes and is done once.

=== Hash maps and #idx("k-mer")k-mer indices

The second indexing approach drops sortedness and replaces it with
hashing. Pick a length $k$ — typically 12 to 30. Enumerate every
$k$-mer in $R$, and store a map from each distinct $k$-mer to the
list of positions where it occurs.

#figure(
  image("../../diagrams/lecture-02/04-hash-kmer-index.svg", width: 95%),
  caption: [
    A k-mer #idx("hash index")hash index. Each k-mer is a key; the value is the list of
    positions in the reference where it occurs. Lookup is a single
    hash operation; memory is dominated by the lists of positions.
  ],
) <fig:hash-kmer>

To query, hash the first $k$ characters of the read and look up the
position list. Each entry is a candidate — the read might occur there,
or it might just share a $k$-mer there by coincidence. You then verify
by checking the rest of the read. Lookup is $O(1)$ per query; the
verification is the expensive part, and is where the false-positive
rate from short $k$ shows up.

The memory cost is the catch. The hash table has up to $4^k$ buckets,
and across all buckets the position lists sum to $|R|$ entries. For
$k = 12$, the bucket count is $4^(12) approx 1.7 times 10^7$, which
fits easily; for $k = 15$, it is $4^(15) approx 10^9$, which starts to
hurt; beyond that, most buckets are empty and you switch from a dense
array of buckets to a hash table keyed by the bucket itself.

The choice of $k$ is a knob bioinformaticians tune constantly. Larger
$k$ means fewer false candidate positions per lookup, because each
$k$-mer is more specific. But it also means fewer $k$-mers in any one
read survive sequencing errors intact — a read of length 150 with one
error covers $150 - k + 1$ $k$-mers, of which $k$ contain the error.
For $k = 20$ that is 131 $k$-mers per read, 20 of which are destroyed
by a single error, leaving 111 clean seeds; for $k = 100$ you might
have no clean seed at all. The right $k$ balances specificity against
robustness to error, and the balance is one of the things every
aligner exposes as a configurable parameter.

=== Where the memory goes

It is worth tabulating the two indices side-by-side before moving on.
For the human genome:

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, left, left, left),
    stroke: 0.5pt + rule,
    table.header(
      [*Index*], [*Query time*], [*Memory*], [*Build*]
    ),
    [Brute force],     [$O(|R| dot |r|)$ per read],   [$|R|$ bytes],                                 [none],
    [Suffix array],    [$O(|r| log |R|)$],            [$approx 5 |R|$ bytes],                        [$O(|R|)$, minutes],
    [Hash $k$-mer],    [$O(|r|)$ amortised],          [$O(|R| + 4^k)$ bytes],                        [$O(|R|)$],
  ),
  caption: [
    Three candidate indices for genome-scale exact matching. The
    suffix array and hash index both have acceptable query speed and
    unacceptable memory; the question that motivates the next section
    is whether either property can be improved without giving up
    the other.
  ],
) <fig:tradeoff-table>

The suffix array gives logarithmic queries and costs roughly 15 GB for
human. The hash table gives near-constant queries and costs a similar
amount, with the $4^k$ overhead growing with $k$. Neither is crazy by
2020s server standards — a 64 GB workstation handles either — but
both are painful on a laptop and both preclude the common case of
running several alignment jobs at once on a shared machine. The
question that drives the next section is: can we build an index that
is _comparable_ in query speed and _much_ smaller, ideally close to
the size of the compressed reference itself?

#warn[
  When a paper says "we used a hash-based aligner" or "we used a
  suffix-tree aligner," the axis to ask about is not speed — both
  are fast — but _memory_, and specifically memory per reference.
  This is what drove the field to the FM index: not a query-speed
  argument, a memory argument.
]


== The #idx("Burrows-Wheeler transform")Burrows-Wheeler Transform and the FM Index <sec:bwt-fm>

The answer is yes, and the construction is one of the most elegant
algorithms in computer science. The Burrows-Wheeler transform is a
reversible permutation of a string that, paired with a small number of
auxiliary tables, supports exact-match queries in $O(|r|)$ time using
memory close to the size of the original string. It also happens to be
the core of `bzip2`, which is not a coincidence and turns out to be
deeply informative.

=== Construction

The _Burrows-Wheeler transform_ (#idx("BWT")BWT) was introduced by Michael Burrows
and David Wheeler in a 1994 Digital Equipment Corporation technical
report. The report was a compression algorithm, not a string-search
data structure; the search applications would come six years later. The
construction is simple enough to do by hand.

Take the string `ACAACG$`. Write down all cyclic rotations:

```
ACAACG$
CAACG$A
AACG$AC
ACG$ACA
CG$ACAA
G$ACAAC
$ACAACG
```

Sort these rotations lexicographically, treating `$` as smaller than
every letter:

```
$ACAACG
AACG$AC
ACAACG$
ACG$ACA
CAACG$A
CG$ACAA
G$ACAAC
```

The BWT is the _last column_ of this sorted matrix, read top to
bottom: `GC$AAAC`. That is the entire construction. Sort all rotations
of the string with `$` appended, take the last column, you have the
BWT. The first column of the same matrix, called $F$, is just the
characters of the string in sorted order; the last column is called
$L$, and is the transform itself.

#figure(
  image("../../diagrams/lecture-02/05-bwt-rotation-matrix.svg", width: 90%),
  caption: [
    Constructing the BWT of `ACAACG$`. The seven cyclic rotations on
    top are sorted lexicographically; the result on the bottom has
    its first column ($F$) and last column ($L$) labelled. The BWT
    is the $L$ column.
  ],
) <fig:bwt-construction>

Why does this help compression? Because the sort groups rows by the
character that _follows_ them — every rotation whose next character is
an `A` ends up adjacent. The previous character of each of those
rotations, the one that lands in $L$, is drawn from the small set of
characters that precede an `A` in the original string. In natural
language or #idx("DNA")DNA the distribution of "what character precedes which" is
far from uniform; in English text, `h` very often follows `t`, in DNA
the dinucleotide #idx("CpG")CpG is depleted, and so on. That non-uniformity makes
$L$ end up with long runs of the same character — `AAA`, `GG`, `CC` —
which run-length-encode and Huffman-code beautifully. That is why
`bzip2` works.

Why does it help search? Because of a property called the #idx("LF mapping")LF mapping,
which we get to next. The short answer: the BWT is a representation of
$R$ that uses about the same space as $R$ itself and, with a few small
auxiliary arrays, supports exact-match queries of a read in $O(|r|)$
time. No $log |R|$ factor, no $4^k$ bucket table.

#note[
  The BWT is a _reversible transform that moves the string into a
  domain where the operations we care about are cheap_. Compare with
  what the discrete Fourier transform does for a signal: it takes the
  signal from the time domain, where convolution is expensive, to the
  frequency domain, where convolution is pointwise multiplication.
  The discrete cosine and the wavelet transforms do the same trick
  for different operations. The BWT is the string-domain analogue.
  Same pedagogical shape: transform, operate, invert.
]

=== The LF mapping and the inverse transform

The BWT would be useless if it were a one-way function. It is not.
There is a property of the sorted-rotation matrix called the _LF
mapping_, or the _last-first property_, that lets us reconstruct the
original string from $L$ alone — given that we also know $F$, which
is just $L$ sorted.

Here is the property, stated plainly. _The $i$-th occurrence of
character $c$ in $L$ corresponds to the $i$-th occurrence of $c$ in
$F$._ Not the $i$-th occurrence in positional order — the $i$-th
occurrence when you walk down each column top to bottom.

Why is this true? Look at any two rows of the sorted matrix that
happen to end with the same character — say, two rows ending in `A`.
Each of those rows is some cyclic rotation of $R$ that ends in `A`.
The two rotations you get by shifting each of those rows one step to
the right — moving their last character to the front — are themselves
rotations of $R$ that _begin_ with `A`. Because the original matrix
was sorted, and because right-shifting by one character is
order-preserving on rotations of a fixed string, the relative order
of those two right-shifted rotations is the same as the relative order
of the original two. So the order of `A`s down $L$ matches the order
of `A`s down $F$. This is the invariant. Work through a small example
by hand once and it becomes obvious; describe it in words and it
sounds mysterious. The figure below makes the connectors explicit.

#figure(
  image("../../diagrams/lecture-02/06-bwt-lf-mapping.svg", width: 90%),
  caption: [
    The LF mapping. Curved connectors show that the $i$-th occurrence
    of each character in $L$ is the $i$-th occurrence in $F$. This
    invariant is what makes both inversion and #idx("backward search")backward search work.
  ],
) <fig:lf>

In algebra, the LF function on row indices is:

$ "LF"(i) = C[L[i]] + "rank"(L[i], i) $

where $C[c]$ is the count of characters in $R$ strictly less than $c$
(the position in $F$ where the block of $c$ starts), and
$"rank"(c, i)$ is the number of occurrences of $c$ in $L[0 .. i)$.
Notice the algebra is just a careful spelling of "$i$-th occurrence of
the same character" — the $C$ table tells you where each character's
block starts in $F$; the rank tells you which occurrence of that
character you are.

Given LF, reconstructing $R$ from $L$ is mechanical. Start at the row
of the matrix whose $F$ character is `$` — row 0 by construction.
Read off the $L$ character of that row: it is the character that
precedes `$` in the original string, which is the last real character
of $R$. Follow LF from that row to land on a new row whose $F$
character is the second-to-last of $R$. Read its $L$, follow LF
again, and so on. In $|R|$ steps you have read out $R$ in reverse.

#figure(
  image("../figures/ch02/f2-bwt-lf-reconstruction.svg", width: 95%),
  caption: [
    Inverting the BWT of `ACAACG$` by walking LF. Starting at row 0
    (the sentinel row), each step reads $L$ at the current row, jumps
    to the same-rank occurrence in $F$, and prepends the read
    character. After six steps the original `ACAACG$` has been
    reconstructed.
  ],
) <fig:bwt-inverse>

The fact that this works at all is one of those things that feels
mysterious the first time and obvious the third. Burrows and Wheeler's
original paper presents the inverse with about a half-page of prose;
it is worth reading.

=== Backward search: the FM index

The _FM index_, named by its inventors Paolo Ferragina and Giovanni
Manzini in 2000, is the BWT plus a small number of auxiliary arrays,
wired together so that exact-match queries can be answered by a
procedure called _backward search_. The dominant index structure in
modern genome aligners is the FM index, for one reason: it answers
exact-match queries in $O(|r|)$ time using memory proportional to the
compressed reference.

The auxiliary arrays are two. First, the $C$ table we already saw:
for each character $c$ in the alphabet, the number of characters in
$R$ that are lexicographically smaller than $c$. This tells you where
the block of rows in $F$ whose first character is $c$ begins. Second,
the rank function $"rank"(c, i)$: the number of occurrences of $c$ in
$L[0 .. i)$.

The search processes the read _right-to-left_ — hence "backward." At
each step, it maintains a half-open interval $[s p, e p)$ of rows in
the sorted matrix whose $F$-column prefix matches the suffix of the
read seen so far. Initially $s p = 0$ and $e p = |R| + 1$ — the whole
matrix. When you extend the match by one more character $c$ on the
left of what you have already matched, you update:

$ s p_"new" = C[c] + "rank"(c, s p) quad e p_"new" = C[c] + "rank"(c, e p) $

If $s p_"new" < e p_"new"$, the query so far matches somewhere in $R$;
continue. If $s p_"new" gt.eq e p_"new"$, the interval has collapsed
to empty and the read does not occur in $R$ — and you know exactly
which character killed it.

The recurrence is the LF function applied to the boundaries of the
interval. Pause on that for a moment. _Backward search is LF mapping,
performed on a range instead of a single row._ The same algebra that
inverts the BWT also searches it. That is what makes the FM index a
single data structure rather than two.

Walk through a query by hand: search for `ACA` in `ACAACG$`, whose
BWT is $L = $ `GC$AAAC`. The right-to-left trace lives in @fig:fm-search.

#figure(
  image("../figures/ch02/f1-fm-backward-search-walkthrough.svg", width: 100%),
  caption: [
    Backward-search trace for the query `ACA` in `ACAACG$`. The
    half-open interval $[s p, e p)$ starts at $[0, 7)$ — the whole
    matrix — and shrinks to $[2, 3)$ after three rank updates,
    pointing at row 2 of the sorted matrix, whose suffix starts at
    position 0 of the original string.
  ],
) <fig:fm-search>

The whole search was three rank queries and three constant-time
arithmetic updates. Three steps, because $|r| = 3$. The reference was
never looked at during the search — only the BWT and the rank arrays
were consulted. That is the point.

=== Checkpoints and the space-time trade-off

There is a problem hidden in the description above. Computing
$"rank"(c, i)$ on demand requires scanning $L$ from the beginning up
to position $i$ and counting occurrences of $c$, which is $O(|R|)$
per query — worse than brute force. Storing the rank function fully
precomputed for every position of $L$ and every character costs
$O(|R| dot |Sigma|)$. For the human genome and a four-character DNA
alphabet, that is roughly 48 GB of 32-bit integers — four times the
suffix array and firmly in no-can-do territory. Neither extreme works.

The trick every production FM-index implementation uses is to take the
middle path. Precompute and store the rank function at every $d$-th
position of $L$, called a _checkpoint_, and count the residue on the
fly from the nearest stored checkpoint to the query position.
Memory drops by a factor of $d$; query time rises by a factor of $d$,
with a small constant per character of scan, typically implemented
with SIMD `POPCNT` instructions so the scan moves at bytes per
nanosecond.

#figure(
  image("../../diagrams/lecture-02/07-fm-checkpoints-tradeoff.svg", width: 95%),
  caption: [
    Dense versus sparse rank arrays. Storing rank at every position
    gives $O(1)$ queries but costs $O(|R| dot |Sigma|)$ memory;
    storing it at every $d$-th position scales memory by $1/d$ and
    query time by $d$. The trade-off curve is the knob the
    implementer turns.
  ],
) <fig:checkpoints>

For the human genome, $d = 128$ is typical. Memory for the rank
arrays drops from 48 GB to a few hundred megabytes, and the extra
scan adds a nanosecond or two per character of read. The index as a
whole — BWT, $C$ table, sparse suffix array (also checkpointed, for
reverse-mapping rows back to original positions), and sparse rank
arrays — fits in about 4 GB for a human reference. `bwa index`
produces it; `bwa mem` uses it.

#tip[
  Same pattern as keeping a coarse table of transcendental-function
  values in firmware and interpolating between samples — and the same
  trade-off applies. The checkpoint interval $d$ is the knob that
  trades memory for compute. Pick it to match the cache hierarchy:
  a checkpoint stride small enough that the residue scan fits in L1
  is much faster than a stride that forces a trip to L2 every query.
  This is the kind of tuning that distinguishes a good FM-index
  implementation from a great one.
]

The exact-match index is now done. The FM index gives you $O(|r|)$
per query in memory close to the compressed reference. Exact matching
of a billion 150 bp reads against the human genome, on this index, is
a problem you can solve in core hours, not core years. What it does
_not_ do is tolerate errors — and reads have errors. That is the next
problem.


== Score-Based Alignment and CIGAR <sec:scoring>

We now have an index that answers exact-match queries in $O(|r|)$
time using memory close to the compressed reference. But reads are
noisy, and every sequenced read differs from its true position by at
least a few bases — so exact-match is not the problem we actually
need to solve. The rest of the chapter builds the error-tolerant
complement of the index and then composes the two.

=== Why exact matching is not enough

Everything in the previous section was exact matching. Given a read
and a reference, report the positions where the read occurs verbatim.
That is useful for some things — finding known primers, counting
$k$-mer occurrences — but it is not what you want as the final output
of a genome aligner. Reads have errors, reads have real variants,
reads have small insertions and deletions. An exact-match-only
aligner would reject almost every read that has anything biologically
interesting about it.

Concrete numbers. Illumina sequencing has a per-base error rate of
about 0.1 to 1 per cent, dominated by substitutions with a small tail
of indels. A 150 bp read therefore has, on average, between 0.15 and
1.5 sequencing-error-induced mismatches — most reads have at least
one. On top of that, a human individual differs from the #idx("GRCh38")GRCh38
reference at roughly 1 in 1,000 positions due to real biological
variation (single-nucleotide variants), plus roughly 1 in 5,000
positions for small indels, plus occasional larger structural
differences. Any aligner worth using must tolerate at least a handful
of differences per read.

#warn[
  It is tempting to think: "just allow up to $k$ mismatches in the FM
  index and be done." The $k$-mismatch backward-search variant does
  exist. It works for small $k$ (typically 1 or 2) by branching at
  each step into all four possible characters that could close the
  interval, but the branching factor grows as $binom(|r|, k) dot
  |Sigma|^k$, and it does not handle insertions or deletions at all.
  For realistic error budgets, you need a different tool.
]

What we need is a definition of "best approximate match" and an
algorithm that finds it. The definition comes from _scoring_; the
algorithm is _dynamic programming_.

=== Scoring matches, mismatches, and gaps

An alignment between two strings is a way of writing them one above
the other, inserting gap characters so they end up the same length and
each column has either a pair of characters or one character paired
with a gap:

```
Read:      G A - T A C A
Reference: G A T T A C A
```

The read has a gap at position 2: a deletion in the read relative to
the reference, or equivalently an insertion in the reference. Which
one you call it depends on which string you call the query. Columns
with two matching characters are matches; columns with two mismatching
characters are mismatches; columns with a gap are indels.

A _scoring scheme_ assigns a number to each kind of column. The
simplest scheme has three parameters: a match reward (say, $+2$), a
mismatch penalty (say, $-1$), and a gap penalty (say, $-2$). The score
of an alignment is the sum of its column scores. The _best_ alignment
is the one with the highest score.

More realistic schemes use _affine gap penalties_: a larger gap-open
penalty to start a run of gaps, and a smaller gap-extend penalty for
each additional gap. This reflects the biology: when real indels
happen, they are often several bases long, so a run of five gaps
should not cost five times the price of one. Typical numbers for DNA
alignment are match $+1$, mismatch $-4$, gap-open $-6$, gap-extend
$-1$. For protein alignment, the match/mismatch pair is replaced by a
full $20 times 20$ #idx("substitution matrix")substitution matrix — #idx("BLOSUM62")BLOSUM62 is the standard —
because not all amino-acid substitutions are equally plausible.

#note[
  The scoring scheme is a log-likelihood model in disguise. If each
  column of an alignment is an independent observation drawn from a
  probability distribution — match with probability $p_m$, mismatch
  with some lower probability, #idx("indel")indel with some lower probability still
  — then the log-likelihood of the alignment is a weighted sum of
  column counts. The "match reward" is the log-odds of a match
  relative to a mismatch under that model. #idx("affine gap")Affine gap penalties
  correspond to a geometric-length model for indels, which is
  approximately correct for real biological indels. The dynamic
  programming we are about to do is, in every meaningful sense,
  maximum-likelihood decoding under a generative error model — it
  just predates that language by a decade.
]

=== Smith-Waterman: local alignment by #idx("dynamic programming")dynamic programming

The _Smith-Waterman_ algorithm, introduced by Temple Smith and Michael
Waterman in 1981, finds the highest-scoring _local_ alignment between
two sequences — the best-scoring pair of substrings, one from each,
under the chosen scoring scheme. Its global counterpart is
#idx("Needleman-Wunsch")Needleman-Wunsch (1970), which forces both strings to align
end-to-end. For read alignment we almost always want local: the read
should align to _some_ window of the reference, not to the whole
three-billion-base thing.

The algorithm fills a two-dimensional matrix $H$, indexed by the
prefixes of the two sequences. $H[i][j]$ holds the best score
achievable by a local alignment that ends at position $i$ in the
first sequence and position $j$ in the second. The recurrence is:

$ H[i][j] = max cases(
  0,
  H[i-1][j-1] + s(a_i, b_j) quad &"(diagonal: match / mismatch)",
  H[i-1][j] + g                  quad &"(up: gap in sequence 2)",
  H[i][j-1] + g                  quad &"(left: gap in sequence 1)",
) $

where $s(a_i, b_j)$ is the match-or-mismatch score and $g$ is the gap
penalty. The cost is $O(|r| dot |R'|)$ time and $O(|r| dot |R'|)$
space — both products over the lengths of the two strings being
compared. Here $R'$ is the small reference window we are aligning the
read against, not the whole genome; that distinction is what makes the
algorithm fit into a budget.

#note[
  The $0$ at the top of the recurrence is what makes the algorithm
  _local_. Without it, every cell inherits the best score reachable
  from the origin — any prefix that scored negative drags the whole
  alignment down with it. The zero floor lets an alignment forget its
  past and restart from any cell, so the algorithm finds the
  best-scoring substring-pair rather than the best end-to-end
  alignment. Delete the zero and you have Needleman-Wunsch.
]

To find the best alignment, locate the cell with the maximum value in
$H$ — that is the endpoint of the best local alignment. Then _trace
back_: at each cell, you know which of the four candidates won the
max, so you know which neighbour to step to. Keep stepping until you
hit a zero. The path you walked is the alignment.

#figure(
  image("../../diagrams/lecture-02/08-smith-waterman-matrix.svg", width: 95%),
  caption: [
    A filled Smith-Waterman matrix for two short strings. The global
    maximum is the endpoint of the best local alignment; the
    traceback path, walked back to a zero, reconstructs the
    alignment in reverse.
  ],
) <fig:sw>

#note[
  Smith-Waterman is _Viterbi on a 2D trellis_. State is (position in
  sequence 1, position in sequence 2); transitions are the three
  step types (diagonal, up, left); the recurrence is the same
  max-over-predecessors-plus-transition-cost; the traceback is the
  same argmax path. The only difference from a convolutional-code
  #idx("Viterbi")Viterbi decoder is the trellis shape and the zero-floor that makes
  the alignment local. If you have ever implemented Viterbi for an
  error-correcting code, you have implemented Smith-Waterman without
  knowing it.
]

There is a single cardinal rule about Smith-Waterman in practice. You
cannot run it against the whole human genome for every read. A single
read-vs-genome Smith-Waterman fills $150 times 3 times 10^9 approx
5 times 10^(11)$ cells. Multiplied by a billion reads, that is
$5 times 10^(20)$ cell updates per dataset. The brute-force number is
back, with worse constants.

#warn[
  Smith-Waterman is the _wrong algorithm to run on the whole
  reference_. Section 2.5 exists because of this. In practice, you
  only ever run Smith-Waterman on small windows of reference that an
  index has already flagged as candidates — a few hundred base pairs
  per candidate position, not three billion. That brings the cost
  per read down to something that fits in microseconds.
]

=== The band optimisation

There is one more efficiency move worth understanding before we get to
#idx("seed-and-extend")seed-and-extend, because real implementations use it everywhere.

If a seed has already told you _approximately_ where in the reference
window the read should align, the optimal alignment cannot stray far
from the diagonal that the seed defines. A read of length 150 with at
most 5 indels can only have its alignment path drift up or down by 5
cells from the seed diagonal — anywhere further from the diagonal,
the cost of getting there exceeds any plausible match score and the
cell is guaranteed to be sub-optimal.

So you only compute cells inside a narrow _band_ along the seed
diagonal. The band width is a small constant, typically 100 cells for
short reads. Time and space drop from $O(|r| dot |R'|)$ to $O(|r|
dot b)$, with $b$ the band width.

#figure(
  image("../figures/ch02/f3-banded-smith-waterman.svg", width: 100%),
  caption: [
    Full Smith-Waterman versus banded Smith-Waterman. On the left,
    every cell of the DP matrix is computed; on the right, only
    cells inside the cobalt band are filled, the rest skipped. The
    optimal path cannot leave the band as long as the read's error
    budget is bounded.
  ],
) <fig:banded>

Banded alignment is what `bwa mem` and `bowtie2` use during the
extension phase, and what `minimap2` uses to fill in base-level
alignments around chained seeds. Combined with vectorised
implementations that compute many cells per CPU instruction — the
SSE-accelerated stripe-based layout from Michael Farrar's 2007 paper
is the canonical reference — banded Smith-Waterman runs at memory
bandwidth on a modern core, which is the speed limit nobody has
crossed yet.

=== CIGAR strings

The output of Smith-Waterman is a path through a matrix. That path,
serialised as a sequence of edit operations, is what gets stored in
an alignment file. The serialisation format is the _CIGAR string_,
defined as part of the #idx("SAM")SAM/BAM specification.

A CIGAR string is a sequence of `<length><operation>` tokens. The
operations are:

- `M` — alignment match, ambiguous between sequence match and mismatch
  (widely used but imprecise).
- `=` — sequence match (explicit).
- `X` — sequence mismatch (explicit).
- `I` — insertion relative to the reference (a base in the read with
  no counterpart in the reference).
- `D` — deletion relative to the reference (a base in the reference
  with no counterpart in the read).
- `S` — #idx("soft clip")soft clip (read bases that were not aligned but are still
  stored in the BAM record).
- `H` — hard clip (read bases that were trimmed off entirely and not
  retained).
- `N` — skipped region (used for #idx("RNA-seq")RNA-seq, where introns produce large
  reference gaps that are not really deletions).
- `P` — padding (multi-alignment use, rarely seen).

So a CIGAR of `4=1X2=2D3=1I3=` reads as: four matches, one mismatch,
two matches, two deletions from the reference, three matches, one
insertion into the read, three matches. That uniquely encodes the
alignment path, modulo the sequences themselves.

#figure(
  image("../../diagrams/lecture-02/09-cigar-rle-diagram.svg", width: 95%),
  caption: [
    A CIGAR string is run-length encoding of the traceback. The raw
    operation sequence on top is bracketed into runs; the final CIGAR
    string at the bottom is what the BAM record stores.
  ],
) <fig:cigar>

#note[
  The CIGAR string is run-length encoding of the edit-operation
  sequence. The raw traceback is a sequence over a small alphabet
  (`M`, `I`, `D`, sometimes more); many real alignments have long
  runs of `M` interrupted by short `I` or `D` segments. RLE
  compresses this beautifully. The same move is used in raster image
  formats (RLE of scan lines), in fax transmission, and in low-level
  video codecs. CIGAR is to alignment what scan-line RLE is to a
  monochrome image. Nothing more, nothing less.
]

CIGAR strings are what every downstream tool reads. Variant callers
walk them to align pileups across reads. Expression quantifiers count
them to assign reads to transcripts. Visualisation tools use them to
render alignment tracks in IGV or JBrowse. If there is a single
artifact of the alignment step that the rest of the field actually
touches, it is the CIGAR field of the BAM record — column 6, in the
text representation.

#warn[
  The `M` operation is ambiguous: it matches either `=` or `X`. Older
  tools emit `M` and never emit the explicit forms. Newer tools
  sometimes emit `=`/`X`, sometimes `M`, depending on flags. If you
  write a tool that counts mismatches, do not assume that `M` means
  match. Always check both the CIGAR operation and the `MD` tag in
  the BAM record, which encodes mismatches separately. This confusion
  has caused real bugs in published pipelines.
]


== Putting It Together: Seed-and-Extend <sec:seed-extend>

The fast index and the slow scorer combine into the pattern every
production aligner is built on. The composition is mechanical, but
the design choices inside it are where each aligner differs.

=== The two-step pattern

Sections 2.2 and 2.3 gave us an index that answers exact-match
queries in $O(|r|)$ time. Section 2.4 gave us an algorithm that
produces optimal approximate alignments in $O(|r| dot |R'|)$ time
on a window $R'$ of the reference. Neither, alone, is what you want.
The index is fast but intolerant of errors; Smith-Waterman is
accurate but unusably slow against a whole genome. Compose them.

The composition is _seed-and-extend_. In one sentence: use the fast
index to find short exact matches of the read against the reference
(_seeds_), then run Smith-Waterman only on small windows around those
seeds to produce the final approximate alignment. Every production
aligner — `bwa`, `bowtie2`, `minimap2`, `novoalign`, `gem`,
practically any tool worth naming — works this way. Their differences
are entirely in the knobs.

The logic is worth saying out loud. A 150 bp read with a 1 per cent
error rate has 1 or 2 errors on average. For any $k$ shorter than the
expected spacing between errors, _some_ $k$-mer of the read is
error-free. That error-free $k$-mer matches exactly in the reference
at the read's true position. The index finds it in $O(k)$. Smith-
Waterman on a window around that position then produces the full
alignment, including the errors. The argument generalises: as long as
the read's error rate is bounded, _some_ short region of it survives
sequencing intact, and that surviving region is the lever.

#figure(
  image("../../diagrams/lecture-02/10-seed-and-extend-pipeline.svg", width: 100%),
  caption: [
    The whole chapter in one picture. A read enters; the FM index
    returns seeds; Smith-Waterman extends each candidate window; one
    alignment wins and emits a CIGAR string. `bwa mem`, `bowtie2`,
    and `minimap2` differ in the seed type, the extension, and the
    read length they are tuned for, not in the skeleton.
  ],
) <fig:ch02-pipeline>

#note[
  This is the coarse acquisition + fine tracking pattern from
  Section 2.1 again. The FM index is the cheap coarse detector;
  Smith-Waterman is the expensive fine estimator; the aligner pays
  for precision only on the candidates the detector flagged. The
  same pattern shows up in any system that has to localise a signal
  in a large search space under tight latency budgets — radar, GPS,
  speech recognition keyword spotting, network intrusion detection.
  It is one of the canonical engineering moves.
]

=== Three aligners, three sets of knobs

Three short-read and long-read aligners dominate the field. Each is a
different set of choices within the seed-and-extend skeleton.

`bwa mem` (Heng Li, 2013) is the default short-read aligner for human
genomes. It uses the FM index of the reference to find _SMEMs_ —
super-maximal exact matches, exact matches between read and reference
that cannot be extended further in either direction without
introducing a mismatch, and that are not contained in any longer such
match. SMEMs are variable-length: a clean region of the read produces
one long SMEM; a noisy region produces several shorter ones. For each
SMEM of sufficient length, `bwa mem` performs banded Smith-Waterman
extension. The earlier `bwa aln` (Heng Li and Richard Durbin, 2009)
used a different strategy — backward search with at most $k$
mismatches — and was the original #idx("BWA")BWA. `bwa mem` superseded it for
reads longer than 70 bp.

`bowtie2` (Ben Langmead and Steven Salzberg, 2012) also uses an FM
index but with _fixed-length seeds_, typically 22 bp, extracted at
regular intervals along the read. Its extension phase is
SSE-accelerated Smith-Waterman, exploiting the CPU's vector
instructions to compute 8 or 16 cells of the DP matrix in parallel
using stripe-based data layout. The original `bowtie` (Langmead,
Trapnell, Pop, and Salzberg, 2009) was the first widely used
BWT-backed short-read aligner and is what convinced the field that
FM-index-based short-read alignment could run at hundreds of
megabases per CPU-minute in under 2 GB of RAM. `bowtie2` added gapped
alignment to the original's mismatches-only design.

`minimap2` (Heng Li, 2018) is the de facto standard for long reads
— #idx("PacBio")PacBio #idx("HiFi")HiFi, #idx("Oxford Nanopore")Oxford Nanopore — and also handles short reads. Its
seeds are _minimizers_: for each sliding window of length $w$ in the
reference and the read, take the lexicographically smallest $k$-mer
in the window as a representative. This sparsifies the $k$-mer index
by roughly a factor of $w$, cutting memory and I/O while keeping the
property that homologous regions share minimizers with high
probability. After collecting minimizer hits, `minimap2` _chains_
them — groups colinear hits into approximate alignments — and only
then runs banded Smith-Waterman to fill in exact base-level details.

#figure(
  image("../figures/ch02/f4-minimizer-chain-extend.svg", width: 100%),
  caption: [
    The #idx("minimap2")minimap2 skeleton on a long read. Minimizers are sparse
    seeds taken as the lexicographically smallest $k$-mer in each
    sliding window. Hits between read and reference form a dot plot;
    colinear hits chain along the true alignment diagonal, while
    off-diagonal hits are discarded as noise. Banded Smith-Waterman
    runs only on small windows anchored to the chain.
  ],
) <fig:minimap2>

Chaining is essential for long reads because a 10 kb read has
dozens of minimizer hits, most of which should be colinear if the
read is a true match. The chain structure tells you the approximate
alignment before you spend compute on dynamic programming.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, left, left, left, left),
    stroke: 0.5pt + rule,
    table.header(
      [*Aligner*], [*Year*], [*Seed*], [*Extension*], [*Tuned for*]
    ),
    [`bwa mem`],   [2013], [SMEMs (variable)],     [Banded SW],            [Short reads],
    [`bowtie2`],   [2012], [Fixed ~22 bp],         [SSE-accelerated SW],   [Short reads],
    [`minimap2`],  [2018], [Minimizers + chain],   [Banded SW],            [Long reads (and short)],
  ),
  caption: [
    Three production aligners, three sets of choices within the same
    skeleton. Index, seed type, extension, and target read length
    are the four knobs; the structural shape — fast exact seed, slow
    accurate extend — is identical.
  ],
) <fig:aligner-table>

The right aligner for a dataset is usually the one whose tuning
matches the data: `bwa mem` or `bowtie2` for Illumina, `minimap2`
for long reads. But the structure is shared, and once you understand
the pattern, switching between tools is more a matter of file formats
and command-line flags than of fundamentally new ideas.

#note[
  The FM index became practical for genome alignment with the 2009
  #idx("Bowtie")Bowtie paper, and BWA followed the same year with a slightly
  different design. Both demonstrated that BWT-backed short-read
  alignment could fit a human-genome index in roughly 2 GB of RAM,
  more than an order of magnitude smaller than the hash-based
  predecessors (`MAQ`, `ELAND`). The field-wide transition from
  hash-based to BWT-based aligners happened in 2009 to 2011 and has
  not been undone. The exception is long-read alignment, where
  minimap2's hash-based minimizer index won out — partly because
  long reads tolerate denser indices, partly because chaining moved
  the constraint from index size to chain-cost computation.
]


== When the Pattern Breaks <sec:edge-cases>

Seed-and-extend works well for reads that have a single true position
in the reference and a handful of errors. Most reads are like that.
But several cases are not, and a real aligner has to handle them
honestly.

=== Multi-mapping reads

A _multi-mapping_ read comes from a region that occurs more than once
in the genome — a transposable element, a segmental duplication, a
gene with paralogues. About 45 per cent of the human genome is
covered by repeat elements like SINEs and LINEs. A read that falls
entirely inside such an element has many equally good alignments. The
aligner reports one and assigns it a low _mapping quality_ (MAPQ),
which encodes a phred-scaled estimate of the probability that the
reported position is wrong. MAPQ 0 means "I picked one of several
equally good positions essentially at random; do not trust this
position." Most downstream tools filter on MAPQ to avoid drawing
conclusions from reads that could have come from anywhere.

The MAPQ formula varies between aligners and is a fascinating tour of
how aligners model their own uncertainty. `bwa mem`'s mapping quality
is a function of the score difference between the best and second-best
alignment, a regression model fit on simulated data. `bowtie2`'s is
a different function. `minimap2` uses yet another. The numbers are
not strictly comparable across aligners, even though they share a
unit.

=== Structural variants and split alignments

A _structural variant_ produces reads whose alignment is not
contiguous. A read that spans a deletion breakpoint in the sample
matches the reference on one side of the deletion for part of its
length and on the other side for the rest; no single contiguous
Smith-Waterman window aligns the whole thing. Aligners handle this
with _soft-clipping_ — reporting an alignment for the part that fits
and flagging the rest (with `S` in the CIGAR) as unaligned-but-present
— or by reporting a _split alignment_ as a supplementary record. Both
mechanisms exist in SAM/BAM for historical reasons; both are used in
practice. Long-read aligners do this routinely; short-read aligners
do it less, because a 150 bp read rarely spans a #idx("structural variant")structural variant
of interest.

_Chimeric reads_ are the more pathological version of the structural-
variant case: a read that is literally a hybrid of two genomic
regions, usually the result of a library-prep artefact rather than a
real biological event. To the aligner, these look like structural
variants. Distinguishing the two requires signal from many reads and
is a job for the variant caller, not the aligner.

=== Specialised aligners

Bisulfite-converted reads, RNA-seq reads that span splice junctions,
ancient-DNA reads with characteristic damage patterns, and a few
others all need specialised aligners. `bismark` and `methylpy` handle
#idx("bisulfite")bisulfite conversion by allowing a deliberate C-to-T mismatch
asymmetry. `STAR` and `HISAT2` allow reads to align across introns by
extending the CIGAR vocabulary with the `N` operation and by detecting
the splice junctions during seed extension. `mapDamage` models the
post-mortem cytosine deamination that turns ancient-sample bases into
characteristic patterns. The core pattern — index, seed, extend — is
preserved in every one of these tools; the scoring model, the seed
definition, or the handling of clipped bases changes.

#tip[
  When you encounter an unfamiliar aligner, the three questions to
  ask are: what index does it use, what is its seed, and what is its
  extension algorithm. Almost every aligner can be classified along
  those three axes in one paragraph from the methods section. If a
  paper does not answer those three questions clearly, treat it as a
  hint that the design has a hole.
]


== Summary <sec:summary>

- Alignment scale forces indexing. Brute force is correct and
  unusable; every real aligner earns its speed against that baseline.
- The Burrows-Wheeler transform is a reversible permutation that
  moves a string into a domain where exact-match search is cheap —
  the string-domain analogue of moving a signal into the frequency
  domain for convolution.
- The FM index is the BWT plus a rank table, with checkpoints
  trading memory for query time. Under the knob, exact-match
  queries cost $O(|r|)$ time in memory close to the compressed
  reference.
- Smith-Waterman is Viterbi on a 2D trellis. It produces optimal
  local alignments, and its output — serialised as a CIGAR string —
  is the run-length encoding of the traceback path. The band
  optimisation makes it tractable on the small windows the index
  flags.
- Seed-and-extend composes the fast-but-intolerant index with the
  slow-but-accurate dynamic program. Every production aligner is a
  variant of this pattern. The differences are in the seed
  definition, the extension banding, and the tuning for read
  length.


== Exercises <sec:exercises>

*1.* _BWT by hand._ Construct the Burrows-Wheeler transform of the
string `BANANA$`. Show the unsorted rotation matrix, the sorted
rotation matrix, and the final $L$ column. Verify that you get
`ANNB$AA`.

*2.* _Backward search by hand._ Using the BWT you computed in
Exercise 1, build the $C$ table and the rank function (you may write
it as four rows of cumulative counts, one per character). Use
backward search to find all occurrences of `ANA` in `BANANA$`. Show
the $[s p, e p)$ interval at every step. How many occurrences are
there? Where do they start in the original string?

*3.* _Suffix-array baseline._ Implement a naive suffix-array-based
exact matcher in Python. The implementation should: (a) build the
suffix array of a reference string by sorting all suffixes; (b)
answer queries by two binary searches over the array. Compare its
running time against the naive sliding-window brute-force matcher
from Exercise 1 of Chapter 1, on the _E. coli_ K-12 reference genome
(about 4.6 Mb) and 10,000 randomly simulated 100-bp reads. Plot time
per query as a function of reference size by down-sampling the
reference to 10 kb, 100 kb, 1 Mb, and full size.

*4.* _Smith-Waterman by hand._ Fill the full Smith-Waterman matrix for
the read `GCATGCA` against the reference `GATTACA`, using match $+2$,
mismatch $-1$, gap $-2$. Identify the global-max cell and the
traceback path. Write out the optimal local alignment and its CIGAR
string. (Use `=` and `X` for explicit match/mismatch, not the
ambiguous `M`.)

*5.* _Affine gap penalties._ Repeat Exercise 4 with affine gap
penalties: gap-open $-3$, gap-extend $-1$, and the same match and
mismatch scores. (You will need to maintain three matrices: one for
diagonal moves, one for gaps in the reference, one for gaps in the
read. This is the Gotoh (1982) extension to Smith-Waterman.) Does
the optimal alignment change? If so, why?

*6.* _Running `bwa`._ Install `bwa` and `samtools` on a Linux or macOS
machine. Index the _E. coli_ K-12 reference. Simulate 10,000 150-bp
reads from it using `wgsim` or `dwgsim` with default error rates.
Align the simulated reads with `bwa mem`. Use `samtools view` to
inspect the CIGAR strings of the first twenty alignments. Pick one
that contains an indel (`I` or `D`) and one that contains a
soft-clip (`S`). In one sentence each, explain what the alignment is
reporting.

*7.* _Checkpoint stride trade-off._ For a 3 Gb reference, plot the
expected memory #idx("footprint")footprint of the FM-index rank arrays and the
expected query latency, as functions of the checkpoint stride $d in
{1, 4, 16, 64, 256, 1024}$. Assume the residue scan runs at 2 GB/s
(SSE-accelerated `popcnt`) and the rank arrays are 4 bytes per
character per checkpoint. At what stride does the latency cross 100
ns per character? At what stride does memory cross 1 GB?

*8.* _(Open-ended.)_ Pick a specialised aligner from the literature
that is not one of `bwa`, `bowtie2`, or `minimap2`. Read its paper.
Describe in one paragraph: (a) the data it is designed for, (b) the
seed-and-extend choices it makes (index type, seed type, extension
algorithm), and (c) the one design move that distinguishes it from
the canonical short-read or long-read pipeline. Good candidates
include `STAR`, `HISAT2`, `bismark`, `BLASR`, `cuda-seqlib`,
`Strobealign`, or `Winnowmap`.


== Further Reading <sec:further-reading>

- *Burrows, M., and Wheeler, D. J.* (1994). "A Block-Sorting Lossless
  Data Compression Algorithm." _Digital SRC Research Report_ 124. The
  original BWT paper. Compact, readable, and worth reading in
  full — it is twenty pages.
- *Ferragina, P., and Manzini, G.* (2000). "Opportunistic Data
  Structures with Applications." _Proceedings of FOCS 2000_, 390–398.
  The FM-index paper. Dense, but the key construction is in the
  first three pages.
- *Smith, T. F., and Waterman, M. S.* (1981). "Identification of
  Common Molecular Subsequences." _Journal of Molecular Biology_ 147:
  195–197. Two pages. The DP recurrence in its original form.
- *Li, H., and Durbin, R.* (2009). "Fast and Accurate Short Read
  Alignment with Burrows-Wheeler Transform." _Bioinformatics_ 25:
  1754–1760. The original BWA paper.
- *Langmead, B., and Salzberg, S. L.* (2012). "Fast Gapped-Read
  Alignment with Bowtie 2." _Nature Methods_ 9: 357–359. The Bowtie 2
  paper; pair with the 2009 _Genome Biology_ Bowtie paper for the
  full design lineage.
- *Li, H.* (2018). "Minimap2: Pairwise Alignment for #idx("nucleotide")Nucleotide
  Sequences." _Bioinformatics_ 34: 3094–3100. The minimap2 paper.
  The chaining algorithm in its supplementary material is worth the
  detour.
- *Langmead, B.* "Teaching Materials: BWT and FM-Index." Lecture
  slides and worked examples at `langmead-lab.org/teaching-materials`.
  The single best place to develop hand-calculation fluency with
  these structures.
