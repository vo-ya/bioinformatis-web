"""Insert `#idx("term")` markers into chapter Typst sources.

Reads `terms.txt`, walks `book/chapters/ch*.typ`, finds the first
case-insensitive whole-word occurrence of each term in each chapter,
and inserts the index marker immediately before that occurrence.

Skips matches that fall inside:
  - Inline math `$ ... $` (would corrupt math parsing)
  - Inline code `` `...` ``
  - Fenced code blocks ```…```
  - Typst labels `<...>`
  - Existing `#idx(...)` calls (idempotent rerun)

The Typst marker is invisible at render time — only its page
location matters — so inserting it just before the term doesn't
change the visible prose.

Re-runnable. If a chapter is already tagged for a term, the script
moves on.

Usage:
    python3 book/index/tag.py
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
BOOK = HERE.parent
CHAPTERS_DIR = BOOK / "chapters"
TERMS_FILE = HERE / "terms.txt"


def parse_terms(path: Path) -> list[tuple[str, str]]:
    """Return list of (term, sort-key). Empty sort-key means "use term"."""
    out: list[tuple[str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if "=>" in s:
            term, sort = (x.strip() for x in s.split("=>", 1))
        else:
            term, sort = s, ""
        out.append((term, sort))
    return out


def split_safe_spans(text: str) -> list[tuple[int, int]]:
    """Return list of (start, end) character ranges in `text` that are
    safe to insert into — outside math, inline code, fenced code blocks,
    labels, and pre-existing #idx() calls.
    """
    n = len(text)
    safe = [True] * n
    i = 0
    in_fence = False
    while i < n:
        c = text[i]
        # Fenced code block ``` ... ```
        if not in_fence and text[i:i + 3] == "```":
            j = text.find("```", i + 3)
            if j == -1:
                j = n
            else:
                j += 3
            for k in range(i, min(j, n)):
                safe[k] = False
            i = j
            continue
        # Inline code: `...`
        if c == "`":
            j = text.find("`", i + 1)
            if j == -1:
                break
            for k in range(i, j + 1):
                safe[k] = False
            i = j + 1
            continue
        # Inline math: $ ... $
        if c == "$":
            j = text.find("$", i + 1)
            if j == -1:
                break
            for k in range(i, j + 1):
                safe[k] = False
            i = j + 1
            continue
        # Labels: <ch:foo> <sec:foo> <fig:foo> <idx>
        if c == "<":
            j = text.find(">", i + 1)
            if j != -1 and j - i < 80:  # heuristic: real labels are short
                for k in range(i, j + 1):
                    safe[k] = False
                i = j + 1
                continue
        # Existing #idx(...) call — skip the whole call
        if text[i:i + 5] == "#idx(":
            # Find the matching close paren (no nested parens expected)
            j = text.find(")", i + 5)
            if j == -1:
                break
            for k in range(i, j + 1):
                safe[k] = False
            i = j + 1
            continue
        # Cross-references: @sec:foo @fig:foo @ch:foo @something:foo
        # Mask the whole `@xxx:yyy` token so we don't insert a marker
        # between the colon and the label.
        if c == "@":
            m = re.match(r"@\w+:[\w-]+", text[i:])
            if m:
                for k in range(i, i + m.end()):
                    safe[k] = False
                i += m.end()
                continue
        i += 1
    # Compress runs of safe[] into (start, end) ranges.
    ranges: list[tuple[int, int]] = []
    in_range = False
    start = 0
    for k in range(n):
        if safe[k] and not in_range:
            start = k
            in_range = True
        elif not safe[k] and in_range:
            ranges.append((start, k))
            in_range = False
    if in_range:
        ranges.append((start, n))
    return ranges


def find_first_in_safe(text: str, term: str,
                       safe_ranges: list[tuple[int, int]]) -> int | None:
    """Return the start offset of the first whole-word, case-insensitive
    match of `term` that lies entirely inside one of `safe_ranges`.
    Returns None if no such match exists.
    """
    pattern = re.compile(r"(?<![\w-])" + re.escape(term) + r"(?![\w-])",
                         flags=re.IGNORECASE)
    for m in pattern.finditer(text):
        s, e = m.start(), m.end()
        for rs, re_ in safe_ranges:
            if s >= rs and e <= re_:
                return s
        # not safe — keep scanning
    return None


def escape_for_typst_string(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')


def tag_chapter(path: Path, terms: list[tuple[str, str]]) -> tuple[int, int]:
    text = path.read_text(encoding="utf-8")
    tagged = 0
    skipped = 0
    # We apply edits right-to-left so earlier offsets don't shift.
    edits: list[tuple[int, str]] = []
    safe_ranges = split_safe_spans(text)
    for term, sort_key in terms:
        # Skip if this chapter already has an explicit #idx for this term
        # (case-insensitive).
        already = re.search(
            r'#idx\("' + re.escape(term) + r'"', text, re.IGNORECASE)
        if already:
            skipped += 1
            continue
        pos = find_first_in_safe(text, term, safe_ranges)
        if pos is None:
            continue
        if sort_key:
            marker = f'#idx("{escape_for_typst_string(term)}", sort: "{escape_for_typst_string(sort_key)}")'
        else:
            marker = f'#idx("{escape_for_typst_string(term)}")'
        edits.append((pos, marker))
        tagged += 1
    # Apply edits right-to-left.
    edits.sort(key=lambda e: -e[0])
    for pos, marker in edits:
        text = text[:pos] + marker + text[pos:]
    if tagged > 0:
        path.write_text(text, encoding="utf-8")
    return tagged, skipped


def main() -> int:
    terms = parse_terms(TERMS_FILE)
    chapters = sorted(CHAPTERS_DIR.glob("ch*.typ"))
    print(f"  {len(terms)} terms × {len(chapters)} chapters")
    print()
    total = 0
    for ch in chapters:
        tagged, skipped = tag_chapter(ch, terms)
        if tagged or skipped:
            rel = ch.relative_to(BOOK.parent)
            print(f"  {rel}: tagged {tagged}, skipped {skipped}")
        total += tagged
    print()
    print(f"Total markers inserted: {total}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
