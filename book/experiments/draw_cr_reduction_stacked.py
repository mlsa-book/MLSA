"""Infographic: stacked competing-risks reduction.

Same 4-step pipeline as the separate-models version, but Step 2 stacks
the two cause-specific datasets into a single dataset with an extra
column `c` indicating the cause, and Step 3 fits a single hazard model
with a c x x interaction.
"""
import math
from textwrap import dedent


INK         = "#2E3440"
GREY_LINE   = "#A5ACBA"
HEADER_BG   = "#F0F2F5"
ROW_BG      = "#FFFFFF"
LANE_DIVIDE = "#E1E5EC"

C1_HIGH     = "#DAE7F4"   # rows with c = 1
C2_HIGH     = "#FAE2D0"   # rows with c = 2
COMBINED    = "#D4EBD8"
MODEL_FILL  = "#E5E9F0"   # neutral grey for the single stacked model

F_STEP       = 32
F_SUB        = 24
F_BOX_TITLE  = 26
F_TABLE_HDR  = 25
F_TABLE_CELL = 24
F_MATH_BIG   = 29
F_MATH_MED   = 29


# --- primitives -----------------------------------------------------------

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
                    row_highlight=None, row_h=44, header_h=48):
    """Row-by-row highlight: row_highlight is a list of fill colors (one per data row),
    or None to use the default white."""
    n_cols = len(cols)
    if col_widths is None:
        col_widths = [80] * n_cols
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
        fill = row_highlight[r] if row_highlight is not None else ROW_BG
        parts.append(
            f'<rect x="{x}" y="{ry}" width="{total_w}" height="{row_h}" '
            f'fill="{fill}" stroke="{GREY_LINE}" stroke-width="0.8"/>'
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


def box_centered(cx, cy, w, h, fill, title, body_lines, body_font_size=F_MATH_MED):
    x = cx - w / 2
    y = cy - h / 2
    parts = [
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" ry="10" '
        f'fill="{fill}" stroke="{GREY_LINE}" stroke-width="1.2"/>',
        f'<text x="{cx:.1f}" y="{y + 36:.1f}" text-anchor="middle" '
        f'font-size="{F_BOX_TITLE}" font-weight="600" fill="{INK}">{title}</text>',
    ]
    spacing = body_font_size + 12
    total = len(body_lines) * spacing
    y0 = y + (h + 36) / 2 - total / 2 + spacing / 2
    for i, line in enumerate(body_lines):
        ly = y0 + i * spacing
        parts.append(
            f'<text x="{cx:.1f}" y="{ly:.1f}" text-anchor="middle" '
            f'font-size="{body_font_size}" fill="{INK}">{line}</text>'
        )
    return "\n  ".join(parts), (x, y, w, h)


# --- compose --------------------------------------------------------------

W, H = 1900, 920
TOP_MARGIN = 60
MID = H / 2 + TOP_MARGIN / 2

LANE_WIDTHS = [280, 420, 540, 660]
assert sum(LANE_WIDTHS) == W
lane_edges   = [sum(LANE_WIDTHS[:i]) for i in range(len(LANE_WIDTHS) + 1)]
lane_centres = [(lane_edges[i] + lane_edges[i + 1]) / 2 for i in range(4)]

X1, X2, X3, X4 = lane_centres
Y_MID = MID


# Vertical separators between steps — positioned at the MIDPOINT of each arrow
# (computed after element extents are known below).
lane_bg_parts = []   # populated after element widths are known

# Step headers
step_labels = []
step_specs = [
    ("Step 1", "original data"),
    ("Step 2", "stack causes"),
    ("Step 3", "fit one model"),
    ("Step 4", "combine"),
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


# Step 1: original CR data (same as before)
orig_cols = ["t", "e", "x"]
orig_rows = [
    ["3", "1", "x₁"],
    ["4", "2", "x₂"],
    ["2", "0", "x₃"],
]
orig_table, (ox, oy, ow, oh) = table_centered(
    X1, Y_MID, orig_cols, orig_rows,
    title="Original CR data", col_widths=[80, 80, 100],
)


# Step 2: stacked dataset.
#   Each original observation appears twice — once for cause 1 (c=1)
#   and once for cause 2 (c=2). The status column δ is the indicator
#   that the observation experienced *that specific cause*.
stacked_cols = ["t", "x",  "c", "δ"]
stacked_rows = [
    ["3", "x₁", "1", "1"],
    ["4", "x₂", "1", "0"],
    ["2", "x₃", "1", "0"],
    ["3", "x₁", "2", "0"],
    ["4", "x₂", "2", "1"],
    ["2", "x₃", "2", "0"],
]
row_colors = [C1_HIGH, C1_HIGH, C1_HIGH, C2_HIGH, C2_HIGH, C2_HIGH]
stacked_table, (sx, sy, sw, sh) = table_centered(
    X2, Y_MID, stacked_cols, stacked_rows,
    title="Stacked dataset", col_widths=[70, 90, 70, 70],
    row_highlight=row_colors,
)


# Step 3: single stacked model
M_W, M_H = 420, 260
model_box, _ = box_centered(
    X3, Y_MID, M_W, M_H, fill=MODEL_FILL,
    title="Stacked model",
    body_lines=[
        "ĥ(τ | x, c)",
        "",
        "ĥ₁(τ | x) = ĥ(τ | x, c = 1)",
        "ĥ₂(τ | x) = ĥ(τ | x, c = 2)",
    ],
    body_font_size=F_MATH_BIG,
)


# Step 4: same combined-quantities box
CB_W, CB_H = 540, 360
cb_box, _ = box_centered(
    X4, Y_MID, CB_W, CB_H, fill=COMBINED,
    title="Combined quantities",
    body_lines=[
        "Ŝ(τ | x) = exp(−Σₑ Ĥₑ(τ | x))",
        "",
        "F̂ₑ(τ | x) = ∫ Ŝ(u⁻| x) ĥₑ(u | x) du",
    ],
    body_font_size=F_MATH_MED,
)


# Arrows
arrow_ends = [
    (ox + ow + 6,        sx - 6),
    (sx + sw + 6,        X3 - M_W / 2 - 6),
    (X3 + M_W / 2 + 6,   X4 - CB_W / 2 - 6),
]
arrows = [arrow(x1, Y_MID, x2, Y_MID) for x1, x2 in arrow_ends]

# Vertical dividers at the midpoint of each arrow
for x1, x2 in arrow_ends:
    mid = (x1 + x2) / 2
    lane_bg_parts.append(
        f'<line x1="{mid:.1f}" y1="40" x2="{mid:.1f}" y2="{H - 40}" '
        f'stroke="{LANE_DIVIDE}" stroke-width="1.2" stroke-dasharray="6 6"/>'
    )


# --- assemble -------------------------------------------------------------

svg = dedent(f"""\
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}"
     font-family="DejaVu Sans">
  <rect width="{W}" height="{H}" fill="white"/>

  {chr(10).join("  " + p for p in lane_bg_parts)}

  {chr(10).join("  " + p for p in step_labels)}

  {orig_table}

  {stacked_table}

  {model_box}

  {cb_box}

  {chr(10).join("  " + a for a in arrows)}
</svg>
""")

import sys
out = sys.argv[1] if len(sys.argv) > 1 else "/tmp/cr-reduction-stacked.svg"
open(out, "w").write(svg)
print(f"Wrote {out}")
