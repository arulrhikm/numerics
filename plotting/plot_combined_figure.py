"""Combine paper Figure 5(right) with both panels of Figure 6 into one figure.

    (a) bound on gate depth, best curve across orders   [from gate_depth.pdf]
    (b) exact operator-norm error, 8-qubit Heisenberg    [from exact-error.pdf]
    (c) the same under a measurement-noise model         [from exact-error.pdf]

This is a *recomposition*, not a rerun: every curve is lifted, as vector art,
out of the PDFs the paper already uses. Nothing is recomputed, so the panels are
guaranteed to agree with the manuscript.

Why composition rather than replotting: Figure 6 has no source in this repo.
``error_analysis/Exact_extrapolation_error.ipynb`` runs at n = 4 (the paper is
n = 8) and has no savefig, and the random-Pauli-X noise model of panel (c)
exists nowhere in code -- only as prose in paper.tex, which never states the
noise rate. Replotting (c) would mean inventing numbers.

REQUIRED INPUT: images/exact-error.pdf, copied out of the Overleaf project.

Panels are placed at equal *plot-frame height*. Font size per unit frame height
is scale-invariant and happens to agree closely between the two sources, so equal
frames also give near-equal type -- which the script verifies and reports rather
than assumes.

Two deliverables come out of this, and they want different things:

    # manuscript figure -- 165mm wide, all three panels, 6.3-7.3pt
    python plotting/plot_combined_figure.py

    # poster block 7 -- 300mm wide, exact-error pair only, 11.5-13.3pt
    python plotting/plot_combined_figure.py --layout row --panels bc \
           --width 300 --out exact-error-poster

Three panels do not fit the poster's slot: at 300mm they measure 7.1-8.3pt, and
at the poster's own 16pt legibility bar nothing here fits at all. The panels are
lifted art and CANNOT be re-typeset, so the only lever is which ones to keep --
hence the poster gets the two exact-error panels, the strongest result, and the
gate-depth bound moves to a caption line.

Panel letters are re-assigned consecutively over whatever ``--panels`` keeps, so
dropping the middle panel yields "(a) (b)", never "(a) (c)".

MEASUREMENT TRAP, if you check the poster fit after embedding this: placing a
clipped PDF leaves zero-height paths in the content stream that ``get_drawings()``
still reports, ~150mm below the visible ink. They render nothing. Filter out
paths with zero width or height before taking a max, or a good board reads as
badly overflowing.
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

import fitz

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

MM = 72.0 / 25.4  # points per mm

# Candidate locations for the one file this script cannot generate itself.
EXACT_ERROR_CANDIDATES = [
    _ROOT.parent / "images" / "exact-error.pdf",
    _ROOT / "plots" / "exact-error.pdf",
    _ROOT.parent / "images" / "exact_error.pdf",
]
GATE_DEPTH = _ROOT / "plots" / "gate_depth.pdf"
OUT_STEM = _ROOT / "plots" / "combined-bounds-and-error"


def find_exact_error(override: str | None = None) -> Path:
    if override:
        p = Path(override)
        if not p.is_file():
            raise SystemExit(f"--exact-error {p} does not exist")
        return p
    for p in EXACT_ERROR_CANDIDATES:
        if p.is_file():
            return p
    raise SystemExit(
        "exact-error.pdf not found. Copy it out of Overleaf to one of:\n  "
        + "\n  ".join(str(p) for p in EXACT_ERROR_CANDIDATES)
        + "\n\nIt cannot be regenerated here: there is no exact-error script in this\n"
          "repo, and the panel (c) noise model exists only as prose in paper.tex."
    )


def _ink_rects(page: fitz.Page) -> list[fitz.Rect]:
    """Every drawing / text span / image rect on the page, page patch excluded."""
    page_area = page.rect.get_area()
    out = []
    for d in page.get_drawings():
        r = d["rect"]
        if r.get_area() > 0.92 * page_area:   # figure / page background patch
            continue
        out.append(r)
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines", []):
            for s in line["spans"]:
                out.append(fitz.Rect(s["bbox"]))
    for img in page.get_images(full=True):
        out.extend(page.get_image_rects(img[0]))
    return out


def panel_groups(page: fitz.Page, *, min_gap_mm: float = 2.5) -> list[fitz.Rect]:
    """Split a 1xN panel figure into panels by projecting ink onto x.

    Identifying the axes rectangle directly is not portable: matplotlib may or
    may not paint an axes background, and when it does not, a rect-based search
    latches onto legend frames instead. Whitespace between panels, on the other
    hand, is present in every such figure. So: mark which x columns carry ink,
    then cut at runs of empty columns wider than ``min_gap_mm``.
    """
    rects = _ink_rects(page)
    if not rects:
        raise SystemExit("page contains no ink")
    step = 0.5
    x_lo = min(r.x0 for r in rects)
    x_hi = max(r.x1 for r in rects)
    n = max(1, int((x_hi - x_lo) / step) + 1)
    occupied = bytearray(n)
    for r in rects:
        i0 = max(0, int((r.x0 - x_lo) / step))
        i1 = min(n - 1, int((r.x1 - x_lo) / step))
        for i in range(i0, i1 + 1):
            occupied[i] = 1

    min_gap = int(min_gap_mm * MM / step)
    groups, start, gap = [], None, 0
    for i in range(n):
        if occupied[i]:
            if start is None:
                start = i
            gap = 0
        else:
            if start is not None:
                gap += 1
                if gap >= min_gap:
                    groups.append((start, i - gap))
                    start, gap = None, 0
    if start is not None:
        groups.append((start, n - 1))

    out = []
    for a, b in groups:
        gx0, gx1 = x_lo + a * step, x_lo + (b + 1) * step
        sel = [r for r in rects if gx0 - 1 <= (r.x0 + r.x1) / 2 <= gx1 + 1]
        box = sel[0]
        for r in sel[1:]:
            box = box | r
        out.append(box)
    return out


def plot_frame_height(page: fitz.Page, clip: fitz.Rect) -> float:
    """Height of the tallest path inside ``clip`` -- i.e. the axes frame.

    Used as the scale reference instead of the clip height, because the clip
    also contains axis titles and tick labels whose share of the total differs
    between the two source figures.
    """
    best = 0.0
    for d in page.get_drawings():
        r = d["rect"]
        if clip.contains(fitz.Rect(r.x0, r.y0, r.x1, r.y1)) or (
            r.x0 >= clip.x0 - 1 and r.x1 <= clip.x1 + 1
        ):
            best = max(best, r.height)
    return best or clip.height


def text_sizes(page: fitz.Page, clip: fitz.Rect | None = None) -> collections.Counter:
    c: collections.Counter = collections.Counter()
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines", []):
            for s in line["spans"]:
                if clip is None or fitz.Rect(s["bbox"]).intersects(clip):
                    c[round(s["size"], 1)] += len(s["text"].strip())
    return c


def dominant_size(counter: collections.Counter) -> float:
    """Most-used font size, ignoring tiny stragglers."""
    return max(counter.items(), key=lambda kv: kv[1])[0] if counter else 0.0


def ink_bbox(page: fitz.Page, x0: float, x1: float) -> fitz.Rect:
    """Bounding box of every drawing/text/image whose centre lies in [x0, x1].

    Used to grow an axes frame back out to include its tick labels, axis title
    and legend, none of which are inside the frame rectangle.
    """
    box = None
    def add(r):
        nonlocal box
        box = r if box is None else box | r
    for d in page.get_drawings():
        r = d["rect"]
        if x0 <= (r.x0 + r.x1) / 2 <= x1:
            add(r)
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines", []):
            for s in line["spans"]:
                r = fitz.Rect(s["bbox"])
                if x0 <= (r.x0 + r.x1) / 2 <= x1:
                    add(r)
    for img in page.get_images(full=True):
        for r in page.get_image_rects(img[0]):
            if x0 <= (r.x0 + r.x1) / 2 <= x1:
                add(r)
    if box is None:
        raise SystemExit(f"no ink found in x range {x0:.0f}-{x1:.0f}")
    return box


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--height", type=float, default=62.0,
                    help="Common axes height for all panels, mm (default 62).")
    ap.add_argument("--gap", type=float, default=5.0,
                    help="Gutter between panels, mm (default 5).")
    ap.add_argument("--exact-error", type=str, default=None,
                    help="Explicit path to exact-error.pdf (default: search images/, plots/).")
    ap.add_argument("--label-pt", type=float, default=0.0,
                    help="Panel-label point size; 0 = match the panels' own type.")
    ap.add_argument("--layout", choices=("row", "grid"), default="grid",
                    help="'row' (default) = panels in one line. 'grid' = (a) centred on "
                         "top with the rest beneath; only worth it if you keep all three, "
                         "since three across at 165mm puts the type at ~4pt.")
    ap.add_argument("--panels", type=str, default="abc",
                    help="Which panels to include (default 'ac'). Three across at "
                         "\\textwidth measures 3.9-4.6pt, which is unreadable; dropping "
                         "the clean exact-error panel gets the remaining two to 5.6-6.5pt. "
                         "Labels are re-lettered consecutively over whatever is kept.")
    ap.add_argument("--out", type=str, default=None,
                    help="Output stem under plots/ (default 'combined-bounds-and-error'). "
                         "Use a second stem for the poster variant so the manuscript "
                         "figure is not overwritten.")
    ap.add_argument("--width", type=float, default=165.0,
                    help="Target output width in mm (default 165 = LaTeX \\textwidth). "
                         "0 keeps the natural size.")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    exact_path = find_exact_error(args.exact_error)

    gd = fitz.open(GATE_DEPTH)
    ee = fitz.open(exact_path)
    gd_pg, ee_pg = gd[0], ee[0]

    print(f"gate_depth.pdf : {gd_pg.rect.width/MM:.0f} x {gd_pg.rect.height/MM:.0f} mm")
    print(f"{exact_path.name:15s}: {ee_pg.rect.width/MM:.0f} x {ee_pg.rect.height/MM:.0f} mm")

    gd_groups = panel_groups(gd_pg)
    ee_groups = panel_groups(ee_pg)
    print(f"panels found   : gate_depth {len(gd_groups)}, exact-error {len(ee_groups)}")
    for nm, gs, pg in (("gate_depth", gd_groups, gd_pg), ("exact-error", ee_groups, ee_pg)):
        for i, r in enumerate(gs):
            print(f"  {nm}[{i}] x {r.x0:6.1f}-{r.x1:6.1f}"
                  f"  ({r.width/MM:.0f} x {r.height/MM:.0f} mm)"
                  f"  frame {plot_frame_height(pg, r)/MM:.0f} mm")

    if len(gd_groups) < 2:
        raise SystemExit("expected 2 panels in gate_depth.pdf (per-order | envelope)")
    if len(ee_groups) < 2:
        raise SystemExit("expected >=2 panels in exact-error.pdf (clean | noisy)")

    # --- source regions -------------------------------------------------
    # (a) the RIGHT panel of gate_depth.pdf = paper Fig 5(right).
    a_clip = gd_groups[-1]
    # (b) first exact-error panel. (c) the second, plus everything to its right:
    # a third group is the shared colorbar, which belongs with (c).
    b_clip = ee_groups[0]
    c_clip = ee_groups[1]
    for extra in ee_groups[2:]:
        c_clip = c_clip | extra

    # Scale so every panel's PLOT FRAME height matches. Font size per unit frame
    # height is scale-invariant, so equal frames also line the type up -- the
    # clip heights cannot be used for this, since axis titles and tick labels
    # take a different share of the total in the two source figures.
    target_h = args.height * MM
    src = {"a": (gd, gd_pg, a_clip), "b": (ee, ee_pg, b_clip), "c": (ee, ee_pg, c_clip)}
    scales = {k: target_h / plot_frame_height(pg, clip) for k, (_, pg, clip) in src.items()}
    print("scales         : " + ", ".join(f"{k}={v:.3f}" for k, v in scales.items()))

    placed = {}
    for k, (_, pg, clip) in src.items():
        dom = dominant_size(text_sizes(pg, clip))
        placed[k] = dom * scales[k]
        print(f"  ({k}) type {dom:.1f}pt native -> {placed[k]:.1f}pt placed")
    spread = max(placed.values()) - min(placed.values())
    print(f"type spread across panels: {spread:.1f}pt"
          + ("  (OK)" if spread <= 1.5 else "  <-- CHECK, panels will look mismatched"))

    # --- compose --------------------------------------------------------
    gap = args.gap * MM
    keys = [k for k in args.panels if k in "abc"]
    if not keys:
        raise SystemExit("--panels must name at least one of a, b, c")
    widths = {k: src[k][2].width * scales[k] for k in keys}
    heights = {k: src[k][2].height * scales[k] for k in keys}

    if args.layout == "row":
        rows = [keys]
    else:
        # (a) alone on top: it comes from a wider figure, so in a single line all
        # three panels need ~270mm and the type collapses at \textwidth.
        top = [k for k in keys if k == "a"]
        bot = [k for k in keys if k != "a"]
        rows = [r for r in (top, bot) if r]

    row_w = [sum(widths[k] for k in r) + gap * (len(r) - 1) for r in rows]
    row_h = [max(heights[k] for k in r) for r in rows]
    total_w = max(row_w)
    total_h = sum(row_h) + gap * (len(rows) - 1)

    out = fitz.open()
    page = out.new_page(width=total_w, height=total_h)

    label_pt = args.label_pt or max(
        dominant_size(text_sizes(pg, clip)) * scales[k] for k, (_, pg, clip) in src.items()
    )

    # Labels run consecutively over the panels actually included: dropping the
    # clean exact-error panel must not leave a figure labelled "(a) (c)".
    letters = {k: chr(ord("a") + i) for i, k in enumerate(keys)}
    if letters != {k: k for k in keys}:
        print("panel labels: " + ", ".join(f"{k}->({v})" for k, v in letters.items()))

    y = 0.0
    for r, rw, rh in zip(rows, row_w, row_h):
        x = (total_w - rw) / 2          # centre short rows
        for k in r:
            doc, pg, clip = src[k]
            w, h = widths[k], heights[k]
            # Bottom-align within the row: x-axis labels share one baseline.
            rect = fitz.Rect(x, y + rh - h, x + w, y + rh)
            page.show_pdf_page(rect, doc, pg.number, clip=clip)
            page.insert_text(
                fitz.Point(x + 2, rect.y0 + label_pt),
                f"({letters[k]})", fontname="hebo", fontsize=label_pt,
            )
            x += w + gap
        y += rh + gap

    if args.width:
        # Uniform rescale to the target width: type shrinks with it, so report
        # what the reader will actually get on the page.
        f = args.width * MM / total_w
        scaled = fitz.open()
        sp = scaled.new_page(width=total_w * f, height=total_h * f)
        sp.show_pdf_page(sp.rect, out, 0)
        out, page = scaled, sp
        total_w, total_h = sp.rect.width, sp.rect.height
        print(f"rescaled to {args.width:.0f} mm wide (x{f:.3f}); "
              f"type on the page: {min(placed.values())*f:.1f}-{max(placed.values())*f:.1f} pt")

    out_stem = OUT_STEM if not args.out else OUT_STEM.parent / args.out
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    out.save(str(out_stem.with_suffix(".pdf")), garbage=3, deflate=True)
    print(f"\nSaved {out_stem.with_suffix('.pdf')}  "
          f"({total_w/MM:.0f} x {total_h/MM:.0f} mm)")

    png = fitz.open(str(out_stem.with_suffix(".pdf")))[0].get_pixmap(matrix=fitz.Matrix(4, 4))
    png.save(str(out_stem.with_suffix(".png")))
    print(f"Saved {out_stem.with_suffix('.png')}")


if __name__ == "__main__":
    main()
