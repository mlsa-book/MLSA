"""Infographic: cause-specific reduction for competing-risks data.

Three-step pipeline (left -> right):
  cause-specific datasets  ->  hazard models  ->  combined quantities

(The original CR dataset is implicit on the left and not drawn; the reader is
expected to bring their own.)

Layout is symmetric about the horizontal mid-line: cause-1 row and cause-2 row
are mirrored above/below, and the combined-quantities box is centred at the
mid-line.
"""
import math
from textwrap import dedent


# ---- palette -------------------------------------------------------------
INK         = "#2E3440"
GREY_LINE   = "#A5ACBA"
HEADER_BG   = "#F0F2F5"
ROW_BG      = "#FFFFFF"
LANE_DIVIDE = "#E1E5EC"

# Softer pastels
C1          = "#C8DCF0"   # cause-1 box fill (soft blue)
C1_HIGH     = "#DAE7F4"   # cause-1 column highlight
C2          = "#F5CDB5"   # cause-2 box fill (soft peach)
C2_HIGH     = "#FAE2D0"   # cause-2 column highlight
COMBINED    = "#D4EBD8"   # combined box fill (soft mint)


# ---- font sizes ----------------------------------------------------------
F_STEP       = 32
F_SUB        = 24
F_BOX_TITLE  = 26
F_TABLE_HDR  = 25
F_TABLE_CELL = 24
F_MATH_BIG   = 29
F_MATH_MED   = 29


# ---- primitives ----------------------------------------------------------

def arrow(x1, y1, x2, y2, head_len=14, head_w=7, stroke_width=1.8):
    dx, dy = x2 - x1, y2 - y1
    L = math.hypot(dx, dy)
    ux, uy = dx / L, dy / L
    px, py = -uy, ux
    bx, by = x2 - head_len * ux, y2 - head_len * uy
    blx, bly = bx + head_w * px, by + head_w * py
    brx, bry = bx - head_w * px, by - head_w * py
    line = (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{bx:.1f}" y2="{by:.1f}" '
            f'stroke="{INK}" stroke-width="{stroke_width}"/>')
    head = (f'<polygon points="{x2:.1f},{y2:.1f} {blx:.1f},{bly:.1f} {brx:.1f},{bry:.1f}" '
            f'fill="{INK}"/>')
    return line + "\n  " + head


def table_centered(cx, cy, cols, rows, title=None, col_widths=None,
                    highlight_col=None, highlight_color=None,
                    row_h=44, header_h=48):
    """Render a table CENTERED on (cx, cy)."""
    n_cols = len(cols)
    if col_widths is None:
        col_widths = [60] * n_cols
    total_w = sum(col_widths)
    total_h = header_h + len(rows) * row_h
    x = cx - total_w / 2
    y = cy - total_h / 2
    starts = [sum(col_widths[:i]) for i in range(n_cols)]
    parts = []

    if title:
        parts.append(
            f'<text x="{cx:.1f}" y="{y - 14:.1f}" text-anchor="middle" '
            f'font-size="{F_BOX_TITLE}" font-weight="600" fill="{INK}">{title}</text>'
        )

    parts.append(
        f'<rect x="{x}" y="{y}" width="{total_w}" height="{header_h}" '
        f'fill="{HEADER_BG}" stroke="{GREY_LINE}" stroke-width="0.8"/>'
    )
    for i, col in enumerate(cols):
        col_cx = x + starts[i] + col_widths[i] / 2
        col_cy = y + header_h * 0.66
        parts.append(
            f'<text x="{col_cx:.1f}" y="{col_cy:.1f}" text-anchor="middle" '
            f'font-size="{F_TABLE_HDR}" font-style="italic" fill="{INK}">{col}</text>'
        )

    for r, row in enumerate(rows):
        ry = y + header_h + r * row_h
        parts.append(
            f'<rect x="{x}" y="{ry}" width="{total_w}" height="{row_h}" '
            f'fill="{ROW_BG}" stroke="{GREY_LINE}" stroke-width="0.8"/>'
        )
        if highlight_col is not None and highlight_color is not None:
            hx = x + starts[highlight_col]
            parts.append(
                f'<rect x="{hx}" y="{ry}" width="{col_widths[highlight_col]}" '
                f'height="{row_h}" fill="{highlight_color}" '
                f'stroke="{GREY_LINE}" stroke-width="0.8"/>'
            )
        for i, cell in enumerate(row):
            col_cx = x + starts[i] + col_widths[i] / 2
            col_cy = ry + row_h * 0.66
            parts.append(
                f'<text x="{col_cx:.1f}" y="{col_cy:.1f}" text-anchor="middle" '
                f'font-size="{F_TABLE_CELL}" fill="{INK}">{cell}</text>'
            )

    parts.append(
        f'<rect x="{x}" y="{y}" width="{total_w}" height="{total_h}" '
        f'fill="none" stroke="{GREY_LINE}" stroke-width="1.2"/>'
    )
    return "\n  ".join(parts), (x, y, total_w, total_h)


def box_centered(cx, cy, w, h, fill, title, body_lines, body_font_size=F_MATH_BIG):
    x = cx - w / 2
    y = cy - h / 2
    parts = [
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" ry="10" '
        f'fill="{fill}" stroke="{GREY_LINE}" stroke-width="1.2"/>',
        f'<text x="{cx:.1f}" y="{y + 32:.1f}" text-anchor="middle" '
        f'font-size="{F_BOX_TITLE}" font-weight="600" fill="{INK}">{title}</text>',
    ]
    spacing = body_font_size + 12
    total = len(body_lines) * spacing
    y0 = y + (h + 32) / 2 - total / 2 + spacing / 2
    for i, line in enumerate(body_lines):
        ly = y0 + i * spacing
        parts.append(
            f'<text x="{cx:.1f}" y="{ly:.1f}" text-anchor="middle" '
            f'font-size="{body_font_size}" fill="{INK}">{line}</text>'
        )
    return "\n  ".join(parts), (x, y, w, h)


# ---- compose -------------------------------------------------------------

W, H = 1520, 920
TOP_MARGIN = 60
MID  = H / 2 + TOP_MARGIN / 2

# Three lanes: cause-split datasets, hazard models, combined quantities
LANE_WIDTHS = [440, 380, 700]
assert sum(LANE_WIDTHS) == W
lane_edges   = [sum(LANE_WIDTHS[:i]) for i in range(len(LANE_WIDTHS) + 1)]
lane_centres = [(lane_edges[i] + lane_edges[i + 1]) / 2 for i in range(3)]

ROW_OFFSET = 200

X1, X2, X3 = lane_centres
Y_MID = MID
Y_TOP = MID - ROW_OFFSET
Y_BOT = MID + ROW_OFFSET


# Vertical separators between lanes
lane_bg_parts = []
for i in range(1, 3):
    lx = lane_edges[i]
    lane_bg_parts.append(
        f'<line x1="{lx:.1f}" y1="40" x2="{lx:.1f}" y2="{H - 40}" '
        f'stroke="{LANE_DIVIDE}" stroke-width="1.2" stroke-dasharray="6 6"/>'
    )


step_labels = []
step_specs = [
    ("Step 1", "split by cause"),
    ("Step 2", "fit one model per cause"),
    ("Step 3", "combine"),
]
for cx, (title, sub) in zip(lane_centres, step_specs):
    step_labels.append(
        f'<text x="{cx:.1f}" y="60" text-anchor="middle" '
        f'font-size="{F_STEP}" font-weight="700" fill="{INK}">{title}</text>'
    )
    step_labels.append(
        f'<text x="{cx:.1f}" y="98" text-anchor="middle" '
        f'font-size="{F_SUB}" font-style="italic" fill="{INK}">{sub}</text>'
    )


# --- Step 1: cause-1 and cause-2 datasets ---
c1_cols = ["t", "x", "δ₁"]
c1_rows = [
    ["3", "x₁", "1"],
    ["4", "x₂", "0"],
    ["2", "x₃", "0"],
]
c1_table, (c1x, c1y, c1w, c1h) = table_centered(
    X1, Y_TOP, c1_cols, c1_rows,
    title="Cause-1 dataset", col_widths=[80, 100, 80],
    highlight_col=2, highlight_color=C1_HIGH,
)

c2_cols = ["t", "x", "δ₂"]
c2_rows = [
    ["3", "x₁", "0"],
    ["4", "x₂", "1"],
    ["2", "x₃", "0"],
]
c2_table, (c2x, c2y, c2w, c2h) = table_centered(
    X1, Y_BOT, c2_cols, c2_rows,
    title="Cause-2 dataset", col_widths=[80, 100, 80],
    highlight_col=2, highlight_color=C2_HIGH,
)


# --- Step 2: hazard model boxes ---
M_W, M_H = 290, 160
m1_box, _ = box_centered(
    X2, Y_TOP, M_W, M_H, fill=C1,
    title="Cause-1 model",
    body_lines=["ĥ₁(τ | x)"], body_font_size=F_MATH_BIG,
)
m2_box, _ = box_centered(
    X2, Y_BOT, M_W, M_H, fill=C2,
    title="Cause-2 model",
    body_lines=["ĥ₂(τ | x)"], body_font_size=F_MATH_BIG,
)


# --- Step 3: combined quantities ---
CB_W, CB_H = 640, 360
cb_box, _ = box_centered(
    X3, Y_MID, CB_W, CB_H, fill=COMBINED,
    title="Combined quantities",
    body_lines=[
        "Ŝ(τ | x) = exp(−Σₑ Ĥₑ(τ | x))",
        "",
        "F̂ₑ(τ | x) = ∫ Ŝ(u⁻| x) ĥₑ(u | x) du",
    ],
    body_font_size=F_MATH_MED,
)


# --- arrows (symmetric pairs) ---
arrows = []

# Step 1 -> Step 2
arrows.append(arrow(c1x + c1w + 6, Y_TOP, X2 - M_W / 2 - 6, Y_TOP))
arrows.append(arrow(c2x + c2w + 6, Y_BOT, X2 - M_W / 2 - 6, Y_BOT))

# Step 2 -> Step 3
arrows.append(arrow(X2 + M_W / 2 + 6, Y_TOP, X3 - CB_W / 2 - 6, Y_MID - 35))
arrows.append(arrow(X2 + M_W / 2 + 6, Y_BOT, X3 - CB_W / 2 - 6, Y_MID + 35))


# ---- assemble -----------------------------------------------------------

svg = dedent(f"""\
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}"
     font-family="DejaVu Sans">
  <rect width="{W}" height="{H}" fill="white"/>

  {chr(10).join("  " + p for p in lane_bg_parts)}

  {chr(10).join("  " + p for p in step_labels)}

  {c1_table}

  {c2_table}

  {m1_box}

  {m2_box}

  {cb_box}

  {chr(10).join("  " + a for a in arrows)}
</svg>
""")

import sys
out = sys.argv[1] if len(sys.argv) > 1 else "/tmp/cr-reduction.svg"
open(out, "w").write(svg)
print(f"Wrote {out}")
