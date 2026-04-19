# Logo Specification — Bioinformatics for Engineers

> **Purpose**: Build spec for the course mark and wordmark.
> **Status**: FINAL — locked to Variant A3 (helix with base-colored rungs).
> **Companion files**: `../proof_html_resources/logo-proof-v2.html` (canonical visual reference — section "Variant A3"), `website-spec.md` (design system), `homepage-spec.md` (usage on landing page).

---

## 1. Design Rationale

The mark is a **DNA double helix with base-colored rungs**. Two thin gray sine strands (in `--fg-muted`) support six thicker colored bars — the rungs — rendered in the site's `--base-a`, `--base-t`, `--base-g`, `--base-c` colors. The strands recede; the colored rungs are the hero.

Three things this mark does well:

1. **It's unambiguously a DNA helix.** The biology reading is immediate.
2. **The rungs carry information.** Each rung is an actual base color from the course's design system. Every time a student sees a red/teal/amber/green base inside an artifact, they make an unconscious connection to the logo. The logo teaches its own color vocabulary.
3. **It ties into the site identity at a systems level.** Other marks would be decoration; this one is structurally integrated with how the whole site renders DNA.

The mark is always paired with a serif wordmark in Source Serif 4.

---

## 2. Canonical Geometry

All coordinates are in a `viewBox="0 0 72 72"` coordinate space. The display/header variants use the exact geometry below; favicon variants are simplified per §5.

### 2.1 Recessive strands (two sine curves)

```xml
<path d="M 8 36 Q 18 14, 28 36 T 48 36 T 64 36"
      fill="none" stroke="#525252" stroke-width="1.5"
      stroke-linecap="round"/>
<path d="M 8 36 Q 18 58, 28 36 T 48 36 T 64 36"
      fill="none" stroke="#525252" stroke-width="1.5"
      stroke-linecap="round"/>
```

Two quadratic sine curves, 180° out of phase. The first peaks upward; the second peaks downward. Both use `--fg-muted` (`#525252`) so they recede visually.

### 2.2 Six colored rungs

```xml
<line x1="14" y1="30" x2="14" y2="42" stroke="#c4342c" stroke-width="3" stroke-linecap="round"/>
<line x1="22" y1="26" x2="22" y2="46" stroke="#0d7377" stroke-width="3" stroke-linecap="round"/>
<line x1="30" y1="22" x2="30" y2="50" stroke="#b45309" stroke-width="3" stroke-linecap="round"/>
<line x1="38" y1="22" x2="38" y2="50" stroke="#2d7a3e" stroke-width="3" stroke-linecap="round"/>
<line x1="46" y1="26" x2="46" y2="46" stroke="#c4342c" stroke-width="3" stroke-linecap="round"/>
<line x1="54" y1="30" x2="54" y2="42" stroke="#0d7377" stroke-width="3" stroke-linecap="round"/>
```

Six vertical bars at `x = 14, 22, 30, 38, 46, 54` (spaced 8 units apart). Their vertical extent corresponds to the distance between the two strands at that x-position — longest in the middle where the strands are furthest apart (y 22 ↔ 50, length 28), tapering at the ends (y 30 ↔ 42, length 12).

**Color sequence (left to right)**: A · T · G · C · A · T — using `--base-a #c4342c`, `--base-t #0d7377`, `--base-g #b45309`, `--base-c #2d7a3e`.

Stroke width 3, `stroke-linecap="round"`. The rungs are ~2× the weight of the strands — this is what makes them visually dominant.

### 2.3 Full display-variant SVG

```xml
<svg viewBox="0 0 72 72" xmlns="http://www.w3.org/2000/svg"
     role="img" aria-labelledby="logo-title">
  <title id="logo-title">Bioinformatics for Engineers</title>
  <!-- strands (recessive, muted) -->
  <path d="M 8 36 Q 18 14, 28 36 T 48 36 T 64 36"
        fill="none" stroke="#525252" stroke-width="1.5" stroke-linecap="round"/>
  <path d="M 8 36 Q 18 58, 28 36 T 48 36 T 64 36"
        fill="none" stroke="#525252" stroke-width="1.5" stroke-linecap="round"/>
  <!-- rungs (dominant, colored) -->
  <line x1="14" y1="30" x2="14" y2="42" stroke="#c4342c" stroke-width="3" stroke-linecap="round"/>
  <line x1="22" y1="26" x2="22" y2="46" stroke="#0d7377" stroke-width="3" stroke-linecap="round"/>
  <line x1="30" y1="22" x2="30" y2="50" stroke="#b45309" stroke-width="3" stroke-linecap="round"/>
  <line x1="38" y1="22" x2="38" y2="50" stroke="#2d7a3e" stroke-width="3" stroke-linecap="round"/>
  <line x1="46" y1="26" x2="46" y2="46" stroke="#c4342c" stroke-width="3" stroke-linecap="round"/>
  <line x1="54" y1="30" x2="54" y2="42" stroke="#0d7377" stroke-width="3" stroke-linecap="round"/>
</svg>
```

This is the canonical asset. Save as `assets/logo/mark.svg`.

---

## 3. Size Variants

### 3.1 Display / header (72px, 28px)

Uses the canonical geometry above without modification. SVG scales crisply; the same file is used for display and header contexts.

### 3.2 Medium favicon (32×32)

Uses the canonical geometry but drops to four rungs (center four) to reduce visual noise at this size.

```xml
<svg viewBox="0 0 72 72" ...>
  <path d="M 8 36 Q 18 14, 28 36 T 48 36 T 64 36"
        fill="none" stroke="#525252" stroke-width="2" stroke-linecap="round"/>
  <path d="M 8 36 Q 18 58, 28 36 T 48 36 T 64 36"
        fill="none" stroke="#525252" stroke-width="2" stroke-linecap="round"/>
  <line x1="22" y1="26" x2="22" y2="46" stroke="#c4342c" stroke-width="4" stroke-linecap="round"/>
  <line x1="30" y1="22" x2="30" y2="50" stroke="#b45309" stroke-width="4" stroke-linecap="round"/>
  <line x1="38" y1="22" x2="38" y2="50" stroke="#2d7a3e" stroke-width="4" stroke-linecap="round"/>
  <line x1="46" y1="26" x2="46" y2="46" stroke="#0d7377" stroke-width="4" stroke-linecap="round"/>
</svg>
```

Note the stroke weights increase (strands 2, rungs 4) to preserve perceptual weight at smaller pixel dimensions.

### 3.3 Small favicon (16×16)

Simplified further: strands thicker, only three central rungs (A · G · C).

```xml
<svg viewBox="0 0 72 72" ...>
  <path d="M 8 36 Q 22 12, 36 36 T 64 36"
        fill="none" stroke="#525252" stroke-width="3" stroke-linecap="round"/>
  <path d="M 8 36 Q 22 60, 36 36 T 64 36"
        fill="none" stroke="#525252" stroke-width="3" stroke-linecap="round"/>
  <line x1="22" y1="24" x2="22" y2="48" stroke="#c4342c" stroke-width="5" stroke-linecap="round"/>
  <line x1="36" y1="20" x2="36" y2="52" stroke="#b45309" stroke-width="5" stroke-linecap="round"/>
  <line x1="50" y1="24" x2="50" y2="48" stroke="#2d7a3e" stroke-width="5" stroke-linecap="round"/>
</svg>
```

Three rungs preserves the color-coded character of the mark at the smallest size without becoming mush. T is dropped because its dark teal is closest in perceptual weight to the muted strands; A/G/C cover the chromatic range.

---

## 4. Color Treatment by Background

### 4.1 On page background (`--bg`, `#fcfcfa`) and muted background (`--bg-muted`)

Use the canonical colors exactly as specified. Strands `#525252`, rungs in their respective base colors.

### 4.2 On accent background (`--accent`, `#1e3a8a`)

The mark inverts to remain visible on cobalt:

- **Strands** become white at reduced opacity: `stroke="#ffffff"` with `opacity="0.5"` (or equivalent `stroke-opacity`).
- **Rungs** use brightened versions of the base colors — the saturated base colors lose contrast against cobalt, so we use lighter variants:

| Base | Light background | On cobalt |
|---|---|---|
| A | `#c4342c` | `#ff8a7e` |
| T | `#0d7377` | `#5eb8b3` |
| G | `#b45309` | `#ffb566` |
| C | `#2d7a3e` | `#6fc080` |

These are pre-calculated, not arbitrary — they're the base colors shifted toward higher lightness while preserving hue. Do not substitute.

### 4.3 On photographs or complex backgrounds

Don't do this. If you must place the mark on a photographic background, put it inside a solid white or cobalt container with 24px padding.

---

## 5. Wordmark

### 5.1 Primary wordmark

- **Line 1** — `Bioinformatics for Engineers`
  - Font: Source Serif 4
  - Weight: 500
  - Letter-spacing: −0.005em
  - Color: `--fg` (`#0a0a0a`) on light backgrounds, `#ffffff` on cobalt

- **Line 2 (optional, institutional)** — `University of Belgrade · School of Electrical Engineering`
  - Font: Inter
  - Weight: 400
  - Size: ~70% of line 1
  - Color: `--fg-muted` (`#525252`) on light, `rgba(255,255,255,0.7)` on cobalt

### 5.2 Mark + wordmark composition

Mark on the left, wordmark on the right, aligned on vertical centerline of the two wordmark lines combined. Gap between mark and wordmark: 12px at header size, scales proportionally.

```html
<div class="brand">
  <svg class="brand-mark" width="28" height="28" viewBox="0 0 72 72">
    ...
  </svg>
  <div class="brand-wordmark">
    <div class="brand-title">Bioinformatics for Engineers</div>
    <div class="brand-sub">University of Belgrade · School of Electrical Engineering</div>
  </div>
</div>
```

### 5.3 Sizing table

| Context | Mark | Line 1 | Line 2 |
|---|---|---|---|
| Homepage header block | 28px | 16px | 12px |
| Homepage title (mark only, large hero-like) | Not used — title is typographic only | — | — |
| Lecture page header (future) | 24px | 14px | Not shown |
| Social / OG preview image | 120px | 48px | 24px |
| Email signature / PDF | 20px | 13px | 11px |

---

## 6. File Deliverables

Final files produced from this spec:

```
assets/logo/
├── mark.svg                   # Canonical display/header mark (6 rungs)
├── mark-favicon-32.svg        # Medium favicon (4 rungs)
├── mark-favicon-16.svg        # Small favicon (3 rungs)
├── favicon.svg                # Alias → mark-favicon-32.svg (browsers pick)
├── favicon.ico                # Multi-resolution ICO for legacy browsers
├── favicon-32.png             # 32×32 PNG fallback
├── favicon-16.png             # 16×16 PNG fallback
├── og-image.png               # 1200×630 social preview, mark + wordmark on cobalt
└── LICENSE.md                 # "© Vojislav Varjačić — course mark, all rights reserved"
```

Favicon HTML references in `<head>`:

```html
<link rel="icon" type="image/svg+xml" href="/assets/logo/favicon.svg">
<link rel="icon" type="image/png" sizes="32x32" href="/assets/logo/favicon-32.png">
<link rel="icon" type="image/png" sizes="16x16" href="/assets/logo/favicon-16.png">
<link rel="shortcut icon" href="/assets/logo/favicon.ico">
<meta property="og:image" content="/assets/logo/og-image.png">
```

---

## 7. Clear Space

Protect the mark with clear space equal to 1/6 of the mark's height on all four sides. At 28px display, this is ~5px. Do not place text, rules, or other graphics within this space.

---

## 8. Don'ts

- ❌ Never rotate the mark. The helix is read horizontally.
- ❌ Never stretch or distort. Uniform scale only.
- ❌ Never change the rung color sequence (A-T-G-C-A-T). The pattern is intentional and stable.
- ❌ Never substitute lookalike colors for the base-color palette. The `#c4342c / #0d7377 / #b45309 / #2d7a3e` quartet is locked.
- ❌ Never invert the strand/rung hierarchy (i.e., don't make strands dominant and rungs recessive).
- ❌ Never add outlines, shadows, glows, gradients, or 3D effects.
- ❌ Never animate the mark in brand contexts. Animation is reserved for artifacts where it teaches something.
- ❌ Never place the mark on a low-contrast or busy background without a solid containing fill.
- ❌ Never substitute a sans-serif for the wordmark. Source Serif 4 is load-bearing.
- ❌ Never mix the light-mode and cobalt-mode color treatments in the same render.

---

## 9. Accessibility

- Every SVG includes `role="img"` and a `<title>` element with content `Bioinformatics for Engineers`.
- Wordmark text is real HTML text — selectable, searchable, screen-reader-readable. Never baked into the SVG.
- Color contrast: cobalt accent on near-white background exceeds WCAG AAA for large text. Base-color rungs on white all meet WCAG AA for graphical elements.
- The mark conveys no information not present in the adjacent wordmark — screen readers get the full brand via text, not the SVG.

---

## 10. Construction Summary (for Claude Code)

To build the full suite of logo assets, the following files need to be produced from the geometry in §2–§3:

1. `assets/logo/mark.svg` — paste the §2.3 SVG exactly
2. `assets/logo/mark-favicon-32.svg` — paste the §3.2 SVG
3. `assets/logo/mark-favicon-16.svg` — paste the §3.3 SVG
4. `assets/logo/favicon.svg` — copy of `mark-favicon-32.svg`
5. `assets/logo/favicon-32.png` — render `mark-favicon-32.svg` at 32×32 (via any SVG-to-PNG tool, or skip and let the browser rasterize the SVG)
6. `assets/logo/favicon-16.png` — render `mark-favicon-16.svg` at 16×16
7. `assets/logo/favicon.ico` — compose from 16 + 32 PNGs (optional; modern browsers don't need it)
8. `assets/logo/og-image.png` — 1200×630, cobalt background, inverted mark at 120px + wordmark lines 1 and 2, both in white/muted-white, left-aligned with generous padding. Can be hand-built as an HTML+CSS mockup and screenshotted, or constructed via any SVG-to-raster pipeline.

PNG conversion is only strictly required for the OG image (social previews historically needed raster). Favicon PNGs are helpful but optional — modern browsers happily consume SVG favicons.
