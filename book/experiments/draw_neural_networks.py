"""Generate diagrams for the neural networks chapter.

Three figures:
  1. ``ffnn-architecture``  - schematic of a feed-forward network used to
     anchor notation and the layer-wise composition.
  2. ``shared-trunk-cr``    - shared trunk + cause-specific heads for
     competing-risk neural-network architectures.
  3. ``multimodal-fusion``  - modality-specific encoders fused into a shared
     survival head.

Each figure is written as a standalone SVG to the chapter figure folder.
The same SVGs are then rasterised into transparent-background PNGs with
ImageMagick (``convert``) so that the book renders cleanly to both HTML
and PDF.
"""
from __future__ import annotations

import math
import subprocess
from dataclasses import dataclass
from pathlib import Path


FIG_DIR = Path(__file__).resolve().parent.parent / "Figures" / "neuralnetworks"
FIG_DIR.mkdir(parents=True, exist_ok=True)


PAL_INPUT = ("#E8EEF7", "#5B7CA8")
PAL_HIDDEN = ("#FFF1D6", "#C99A2C")
PAL_OUTPUT = ("#E8F5E9", "#3F8A52")
PAL_CR = ("#FDE7E8", "#B45357")
PAL_MOD = ("#F1E6FA", "#7B4FA8")
PAL_FUSION = ("#E8F1F1", "#3E8A8B")
EDGE = "#6E7787"
FONT = "DejaVu Sans"


@dataclass
class Circle:
    cx: float
    cy: float
    r: float

    def boundary_toward(self, x, y):
        dx, dy = x - self.cx, y - self.cy
        d = math.hypot(dx, dy)
        if d == 0:
            return self.cx, self.cy
        return self.cx + self.r * dx / d, self.cy + self.r * dy / d


@dataclass
class Rect:
    x: float
    y: float
    w: float
    h: float

    @property
    def cx(self):
        return self.x + self.w / 2

    @property
    def cy(self):
        return self.y + self.h / 2

    def boundary_toward(self, x, y):
        dx, dy = x - self.cx, y - self.cy
        if dx == 0 and dy == 0:
            return self.cx, self.cy
        ts = []
        if dx > 0:
            ts.append((self.x + self.w - self.cx) / dx)
        if dx < 0:
            ts.append((self.x - self.cx) / dx)
        if dy > 0:
            ts.append((self.y + self.h - self.cy) / dy)
        if dy < 0:
            ts.append((self.y - self.cy) / dy)
        t = min(t for t in ts if t > 0)
        return self.cx + t * dx, self.cy + t * dy


def arrow_svg(x1, y1, x2, y2, stroke=EDGE, sw=1.2, head_len=8, head_w=4,
              opacity=1.0):
    dx, dy = x2 - x1, y2 - y1
    L = math.hypot(dx, dy)
    if L == 0:
        return ""
    ux, uy = dx / L, dy / L
    px, py = -uy, ux
    bx, by = x2 - head_len * ux, y2 - head_len * uy
    bx1, by1 = bx + head_w * px, by + head_w * py
    bx2, by2 = bx - head_w * px, by - head_w * py
    return (
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{bx:.2f}" y2="{by:.2f}" '
        f'stroke="{stroke}" stroke-width="{sw}" stroke-opacity="{opacity}"/>'
        f'<polygon points="{x2:.2f},{y2:.2f} {bx1:.2f},{by1:.2f} '
        f'{bx2:.2f},{by2:.2f}" fill="{stroke}" fill-opacity="{opacity}"/>'
    )


def edge_between(src, dst, sw=0.7, opacity=0.55, *, arrow=False,
                 head_len=6, head_w=3):
    sx, sy = src.boundary_toward(dst.cx, dst.cy)
    tx, ty = dst.boundary_toward(src.cx, src.cy)
    if not arrow:
        return (
            f'<line x1="{sx:.2f}" y1="{sy:.2f}" x2="{tx:.2f}" y2="{ty:.2f}" '
            f'stroke="{EDGE}" stroke-width="{sw}" '
            f'stroke-opacity="{opacity}"/>'
        )
    return arrow_svg(sx, sy, tx, ty, stroke=EDGE, sw=sw,
                     head_len=head_len, head_w=head_w, opacity=opacity)


def circle_node(c, fill, stroke, label="", font_size=14, sw=1.6, label_dy=0):
    parts = [
        f'<circle cx="{c.cx}" cy="{c.cy}" r="{c.r}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{sw}"/>'
    ]
    if label:
        parts.append(
            f'<text x="{c.cx}" y="{c.cy + label_dy + font_size / 3}" '
            f'text-anchor="middle" font-family="{FONT}" '
            f'font-size="{font_size}" fill="#111">{label}</text>'
        )
    return "".join(parts)


def rect_node(r, fill, stroke, label="", font_size=13, sw=1.6, rx=8):
    parts = [
        f'<rect x="{r.x}" y="{r.y}" width="{r.w}" height="{r.h}" rx="{rx}" '
        f'ry="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'
    ]
    if label:
        parts.append(
            f'<text x="{r.cx}" y="{r.cy + font_size / 3}" '
            f'text-anchor="middle" font-family="{FONT}" '
            f'font-size="{font_size}" fill="#111">{label}</text>'
        )
    return "".join(parts)


def text(x, y, label, size=13, weight="normal", anchor="middle", fill="#111"):
    return (
        f'<text x="{x}" y="{y}" text-anchor="{anchor}" '
        f'font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
        f'fill="{fill}">{label}</text>'
    )


def write_svg(name, width, height, body):
    path = FIG_DIR / f"{name}.svg"
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}" '
        f'font-family="{FONT}">{body}</svg>'
    )
    path.write_text(svg)
    return path


# ---- Wiegrebe-style helpers ---------------------------------------------

# Soft palette matching the figures in Wiegrebe et al. (2024).
HDR_BLUE = "#D9E5F4"     # joint / feature subnetwork header fill
HDR_BLUE_LN = "#5B7CA8"
HDR_ORANGE = "#FCE9D6"   # cause-specific / output / reshape header fill
HDR_ORANGE_LN = "#C99A2C"
HDR_GREEN = "#D7EBC5"    # final outcome box fill
HDR_GREEN_LN = "#3F8A52"
LAYER_FILL = "#FFFFFF"
LAYER_LN = "#666666"
DASH = "5,3"             # dashed-border pattern


def dashed_box(rect, header_text, header_fill, header_ln,
               sw=1.2, rx=4):
    """Outer dashed rectangle + solid coloured header bar with bold title."""
    parts = []
    parts.append(
        f'<rect x="{rect.x}" y="{rect.y}" width="{rect.w}" '
        f'height="{rect.h}" rx="{rx}" ry="{rx}" '
        f'fill="white" stroke="{header_ln}" stroke-width="{sw}" '
        f'stroke-dasharray="{DASH}"/>'
    )
    hb_h = 22
    parts.append(
        f'<rect x="{rect.x}" y="{rect.y}" width="{rect.w}" '
        f'height="{hb_h}" rx="{rx}" ry="{rx}" '
        f'fill="{header_fill}" stroke="{header_ln}" stroke-width="{sw}"/>'
    )
    parts.append(text(rect.cx, rect.y + hb_h - 7, header_text,
                      size=12, weight="bold", fill="#222"))
    return "".join(parts)


def layer_box(rect, label, fill=LAYER_FILL, ln=LAYER_LN, sw=1.0):
    parts = []
    parts.append(
        f'<rect x="{rect.x}" y="{rect.y}" width="{rect.w}" '
        f'height="{rect.h}" rx="3" ry="3" fill="{fill}" '
        f'stroke="{ln}" stroke-width="{sw}"/>'
    )
    parts.append(text(rect.cx, rect.cy + 4, label,
                      size=11, fill="#222"))
    return "".join(parts)


def vertical_arrow(x, y0, y1, sw=1.2, head_len=8, head_w=4):
    """Straight downward arrow from (x, y0) to (x, y1)."""
    return arrow_svg(x, y0, x, y1, stroke=EDGE, sw=sw,
                     head_len=head_len, head_w=head_w)


def tabular_icon(cx, cy, cols=4, rows=3, cell=14, ln="#444"):
    """Small grid that represents a table of features."""
    w = cols * cell
    h = rows * cell
    x0 = cx - w / 2
    y0 = cy - h / 2
    parts = []
    for i in range(rows):
        for j in range(cols):
            parts.append(
                f'<rect x="{x0 + j * cell}" y="{y0 + i * cell}" '
                f'width="{cell}" height="{cell}" fill="white" '
                f'stroke="{ln}" stroke-width="1"/>'
            )
    return "".join(parts)


def cylinder_icon(cx, cy, w=46, h=58, fill="white", ln="#444"):
    """Database-cylinder icon (Wiegrebe-style "further modality" symbol)."""
    rx = w / 2
    ry = 7
    top_y = cy - h / 2 + ry
    bot_y = cy + h / 2 - ry
    parts = []
    parts.append(
        f'<path d="M {cx - rx} {top_y} '
        f'L {cx - rx} {bot_y} '
        f'A {rx} {ry} 0 0 0 {cx + rx} {bot_y} '
        f'L {cx + rx} {top_y} '
        f'A {rx} {ry} 0 0 0 {cx - rx} {top_y} Z" '
        f'fill="{fill}" stroke="{ln}" stroke-width="1"/>'
    )
    parts.append(
        f'<ellipse cx="{cx}" cy="{top_y}" rx="{rx}" ry="{ry}" '
        f'fill="{fill}" stroke="{ln}" stroke-width="1"/>'
    )
    for off in (-ry * 0.6, 0, ry * 0.6):
        parts.append(
            f'<ellipse cx="{cx}" cy="{top_y + ry + off + ry}" '
            f'rx="{rx}" ry="{ry}" fill="none" stroke="{ln}" '
            f'stroke-width="0.7" stroke-opacity="0.6"/>'
        )
    return "".join(parts)


def embed_png(x, y, w, h, png_path):
    """Embed an external PNG into the SVG as a base-64 data URI."""
    import base64
    data = Path(png_path).read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    return (
        f'<image x="{x}" y="{y}" width="{w}" height="{h}" '
        f'xlink:href="data:image/png;base64,{b64}" '
        f'href="data:image/png;base64,{b64}" '
        f'preserveAspectRatio="xMidYMid slice"/>'
    )


def rasterise(svg_path, density=300, margin=24):
    out = svg_path.with_suffix(".png")
    subprocess.run(
        [
            "convert",
            "-density",
            str(density),
            "-background",
            "none",
            str(svg_path),
            "-trim",
            "+repage",
            "-bordercolor",
            "none",
            "-border",
            f"{margin}x{margin}",
            str(out),
        ],
        check=True,
    )


# ---- Figure 1: feed-forward neural network ------------------------------

def fig_ffnn():
    """Feed-forward NN diagram: input layer, two hidden layers, output.

    Layout and unit labels mirror the workshop slide
    ``mlcu-dl-intro-workshop/slides/day2-01.qmd`` (3 inputs, 5 + 3 hidden
    units, scalar output ``y``).  Labels use Unicode super/subscript
    characters that render in DejaVu Sans through ImageMagick.  Above each
    column a bold layer caption and the corresponding formula are shown.
    """
    W, H = 920, 480
    parts = []

    x_in, x_h1, x_h2, x_out = 200, 410, 620, 800
    r_node = 24
    spacing = 65
    top_band_y = 30
    centre = H / 2 + 40

    def column(x, labels):
        n = len(labels)
        y0 = centre - (n - 1) * spacing / 2
        ys = [y0 + i * spacing for i in range(n)]
        return [Circle(x, y, r_node) for y in ys]

    in_labels = ["x₁", "x₂", "x₃"]
    h1_labels = ["h₁⁽¹⁾", "h₂⁽¹⁾", "h₃⁽¹⁾", "h₄⁽¹⁾", "h₅⁽¹⁾"]
    h2_labels = ["h₁⁽²⁾", "h₂⁽²⁾", "h₃⁽²⁾"]

    inputs = column(x_in, in_labels)
    h1 = column(x_h1, h1_labels)
    h2 = column(x_h2, h2_labels)
    out = Circle(x_out, centre, r_node + 2)

    # bottom-to-top ordering (h₁ at the bottom, like the slide)
    inputs.reverse()
    h1.reverse()
    h2.reverse()

    for c, lab in zip(inputs, in_labels):
        parts.append(circle_node(c, *PAL_INPUT, label=lab,
                                 font_size=14, sw=1.6, label_dy=-2))
    for c, lab in zip(h1, h1_labels):
        parts.append(circle_node(c, *PAL_HIDDEN, label=lab,
                                 font_size=12, sw=1.4, label_dy=-2))
    for c, lab in zip(h2, h2_labels):
        parts.append(circle_node(c, *PAL_HIDDEN, label=lab,
                                 font_size=12, sw=1.4, label_dy=-2))
    parts.append(circle_node(out, *PAL_OUTPUT, label="y",
                             font_size=15, sw=2, label_dy=-2))

    for a in inputs:
        for b in h1:
            parts.append(edge_between(a, b, sw=0.6, opacity=0.5,
                                      arrow=True))
    for a in h1:
        for b in h2:
            parts.append(edge_between(a, b, sw=0.6, opacity=0.5,
                                      arrow=True))
    for a in h2:
        parts.append(edge_between(a, out, sw=0.8, opacity=0.65,
                                  arrow=True, head_len=8, head_w=4))

    # Top band: bold layer caption + formula for the values in that column.
    parts.append(text(x_in, top_band_y, "input layer",
                      size=13, weight="bold"))
    parts.append(text(x_in, top_band_y + 18,
                      "x ∈ ℝ³",
                      size=12, fill="#444"))

    parts.append(text(x_h1, top_band_y, "hidden layer 1",
                      size=13, weight="bold"))
    parts.append(text(x_h1, top_band_y + 18,
                      "h⁽¹⁾ = σ(W⁽¹⁾x + b⁽¹⁾)",
                      size=12, fill="#444"))

    parts.append(text(x_h2, top_band_y, "hidden layer 2",
                      size=13, weight="bold"))
    parts.append(text(x_h2, top_band_y + 18,
                      "h⁽²⁾ = σ(W⁽²⁾h⁽¹⁾ + b⁽²⁾)",
                      size=12, fill="#444"))

    parts.append(text(x_out, top_band_y, "output layer",
                      size=13, weight="bold"))
    parts.append(text(x_out, top_band_y + 18,
                      "y = W⁽³⁾h⁽²⁾ + b⁽³⁾",
                      size=12, fill="#444"))

    svg_path = write_svg("ffnn-architecture", W, H, "".join(parts))
    rasterise(svg_path)
    return svg_path


# ---- Figure 2: shared-trunk competing-risks network --------------------

def fig_shared_trunk_cr():
    """Wiegrebe-style competing-risks figure.

    Top-to-bottom: tabular data icon -> joint (shared) subnetwork ->
    K cause-specific subnetworks (with output layers) -> concat -> psi()
    transformation -> cumulative incidence functions (CIFs).
    """
    W, H = 920, 980
    parts = []

    cx_fig = W / 2

    # --- 1. Tabular data input ------------------------------------------
    parts.append(text(cx_fig, 28, "Tabular data",
                      size=14, weight="bold", fill="#222"))
    parts.append(tabular_icon(cx_fig, 70, cols=4, rows=3, cell=16))
    parts.append(vertical_arrow(cx_fig, 104, 132))

    # --- 2. Joint subnetwork --------------------------------------------
    joint_w, joint_h = 240, 152
    joint = Rect(cx_fig - joint_w / 2, 138, joint_w, joint_h)
    parts.append(dashed_box(joint, "Joint subnetwork",
                            HDR_BLUE, HDR_BLUE_LN))
    parts.append(layer_box(Rect(joint.x + 50, joint.y + 38,
                                joint.w - 100, 30),
                           "Dense Layer"))
    parts.append(text(joint.cx, joint.y + 92, "⋮", size=18))
    parts.append(layer_box(Rect(joint.x + 50, joint.y + 105,
                                joint.w - 100, 30),
                           "Dense Layer"))

    # --- 3. Fan-out to K=3 cause-specific subnetworks --------------------
    sub_y = 340
    sub_w, sub_h = 200, 130
    cols_x = [80, 360, 640]
    parts.append(vertical_arrow(cx_fig, joint.y + joint.h, sub_y - 6))
    # branching from a centre point above the cause columns
    branch_y = sub_y - 18
    parts.append(
        f'<line x1="{cols_x[0] + sub_w/2}" y1="{branch_y}" '
        f'x2="{cols_x[-1] + sub_w/2}" y2="{branch_y}" '
        f'stroke="{EDGE}" stroke-width="1.2"/>'
    )
    for cx in cols_x:
        parts.append(vertical_arrow(cx + sub_w / 2, branch_y, sub_y - 1,
                                    head_len=7, head_w=3))

    cause_labels = ["Subnetwork (cause 1)", "Subnetwork (cause 2)",
                    "Subnetwork (cause Q)"]
    sub_rects = []
    for x, lab in zip(cols_x, cause_labels):
        r = Rect(x, sub_y, sub_w, sub_h)
        sub_rects.append(r)
        parts.append(dashed_box(r, lab, HDR_ORANGE, HDR_ORANGE_LN))
        parts.append(layer_box(Rect(r.x + 35, r.y + 36,
                                    r.w - 70, 28),
                               "Dense Layer"))
        parts.append(text(r.cx, r.y + 82, "⋮", size=16))
        parts.append(layer_box(Rect(r.x + 35, r.y + 94,
                                    r.w - 70, 28),
                               "Dense Layer"))

    # horizontal dots only between cause 2 and cause Q to suggest variable Q
    mid_y = sub_y + sub_h / 2
    x_mid = (cols_x[1] + sub_w + cols_x[2]) / 2
    parts.append(text(x_mid, mid_y, "⋯", size=22, fill="#555"))

    # --- 4. Output layers ------------------------------------------------
    out_y = 520
    out_h = 130
    output_labels = ["Output Layer 1", "Output Layer 2", "Output Layer Q"]
    out_rects = []
    for r_sub, x, lab in zip(sub_rects, cols_x, output_labels):
        parts.append(vertical_arrow(r_sub.cx, r_sub.y + r_sub.h, out_y - 1))
        r = Rect(x, out_y, sub_w, out_h)
        out_rects.append(r)
        parts.append(dashed_box(r, lab, HDR_ORANGE, HDR_ORANGE_LN))
        parts.append(layer_box(Rect(r.x + 35, r.y + 36,
                                    r.w - 70, 28),
                               "Dense Layer"))
        parts.append(layer_box(Rect(r.x + 35, r.y + 76,
                                    r.w - 70, 28),
                               "Activation"))

    # --- 5. Concat node --------------------------------------------------
    c_cy = 720
    c_r = 18
    parts.append(circle_node(Circle(cx_fig, c_cy, c_r),
                             "white", "#444",
                             label="C", font_size=14, sw=1.4,
                             label_dy=-2))
    for r_out in out_rects:
        # arrow from bottom of each output box, jogging to the C node
        parts.append(arrow_svg(r_out.cx, r_out.y + r_out.h,
                               cx_fig, c_cy - c_r,
                               stroke=EDGE, sw=1.2,
                               head_len=8, head_w=4))

    # --- 6. Transformation* + CIFs ---------------------------------------
    trafo_w, trafo_h = 150, 95
    trafo = Rect(cx_fig - trafo_w / 2, 770, trafo_w, trafo_h)
    parts.append(vertical_arrow(cx_fig, c_cy + c_r, trafo.y - 1))
    parts.append(
        f'<rect x="{trafo.x}" y="{trafo.y}" width="{trafo.w}" '
        f'height="{trafo.h}" rx="4" ry="4" fill="white" '
        f'stroke="#666" stroke-width="1.2" stroke-dasharray="{DASH}"/>'
    )
    parts.append(text(trafo.cx, trafo.y + 18, "Transformation*",
                      size=12, weight="bold", fill="#222"))
    psi_cy = trafo.y + 58
    parts.append(circle_node(Circle(trafo.cx, psi_cy, 20),
                             "white", "#444",
                             label="ψ(·)", font_size=13, sw=1.2,
                             label_dy=-2))

    cif_w, cif_h = 130, 55
    cif = Rect(cx_fig + trafo_w / 2 + 50, trafo.cy - cif_h / 2,
               cif_w, cif_h)
    parts.append(
        f'<rect x="{cif.x}" y="{cif.y}" width="{cif.w}" '
        f'height="{cif.h}" rx="4" ry="4" '
        f'fill="{HDR_GREEN}" stroke="{HDR_GREEN_LN}" stroke-width="1.4"/>'
    )
    parts.append(text(cif.cx, cif.cy + 4, "CIFs",
                      size=13, weight="bold", fill="#222"))
    parts.append(arrow_svg(trafo.x + trafo.w, trafo.cy,
                           cif.x, cif.cy,
                           stroke=EDGE, sw=1.4,
                           head_len=9, head_w=4))

    svg_path = write_svg("shared-trunk-cr", W, H, "".join(parts))
    rasterise(svg_path)
    return svg_path


# ---- Figure 3: multi-modal fusion ---------------------------------------

def fig_multimodal_fusion():
    """Wiegrebe-style multi-modal figure.

    Three columns (Tabular / Image / Further modality) feed into
    modality-specific subnetworks, optionally reshaping their outputs, then
    concatenating to a shared output layer and psi() transformation
    producing the survival function.
    """
    ct_png = FIG_DIR / "synthetic-ct.png"
    if not ct_png.exists():
        subprocess.run(
            ["python3",
             str(Path(__file__).with_name("draw_synthetic_ct.py"))],
            check=True,
        )

    W, H = 980, 1020
    parts = []

    col_xs = [150, 490, 830]
    col_w = 180

    # Column titles
    title_y = 28
    titles = ["Tabular data", "Image data", "Further modality"]
    for x, t in zip(col_xs, titles):
        parts.append(text(x, title_y, t, size=14, weight="bold"))

    # Data icons
    icon_y = 88
    parts.append(tabular_icon(col_xs[0], icon_y, cols=4, rows=3, cell=16))
    img_size = 88
    parts.append(embed_png(col_xs[1] - img_size / 2,
                           icon_y - img_size / 2,
                           img_size, img_size, str(ct_png)))
    parts.append(
        f'<rect x="{col_xs[1] - img_size/2}" y="{icon_y - img_size/2}" '
        f'width="{img_size}" height="{img_size}" '
        f'fill="none" stroke="#444" stroke-width="1"/>'
    )
    parts.append(cylinder_icon(col_xs[2], icon_y, w=50, h=66))

    # Modality-specific subnetworks
    sub_y = 175
    sub_h = 165
    for x in col_xs:
        parts.append(vertical_arrow(x, icon_y + 50, sub_y - 1))

    subnet_specs = [
        ("Tabular subnetwork",
         ["Dense Layer", "⋮", "Dense Layer"]),
        ("Image subnetwork",
         ["Conv2D Block", "⋮", "Conv2D Block"]),
        ("Arbitrary subnetwork",
         ["Arbitrary Layer", "⋮", "Arbitrary Layer"]),
    ]
    sub_rects = []
    for x, (lab, layers) in zip(col_xs, subnet_specs):
        r = Rect(x - col_w / 2, sub_y, col_w, sub_h)
        sub_rects.append(r)
        parts.append(dashed_box(r, lab, HDR_BLUE, HDR_BLUE_LN))
        parts.append(layer_box(Rect(r.x + 30, r.y + 36,
                                    r.w - 60, 28),
                               layers[0]))
        parts.append(text(r.cx, r.y + 92, layers[1], size=16))
        parts.append(layer_box(Rect(r.x + 30, r.y + 108,
                                    r.w - 60, 28),
                               layers[2]))

    # Reshaping for image and further-modality lanes
    resh_y = 380
    resh_h = 100
    img_resh = Rect(col_xs[1] - col_w / 2, resh_y, col_w, resh_h)
    parts.append(vertical_arrow(col_xs[1], sub_rects[1].y + sub_h,
                                img_resh.y - 1))
    parts.append(dashed_box(img_resh, "Reshaping",
                            HDR_ORANGE, HDR_ORANGE_LN))
    parts.append(layer_box(Rect(img_resh.x + 30, img_resh.y + 50,
                                img_resh.w - 60, 28),
                           "Flatten"))
    far_resh = Rect(col_xs[2] - col_w / 2, resh_y, col_w, resh_h)
    parts.append(vertical_arrow(col_xs[2], sub_rects[2].y + sub_h,
                                far_resh.y - 1))
    parts.append(dashed_box(far_resh, "Reshaping (if needed)",
                            HDR_ORANGE, HDR_ORANGE_LN))
    parts.append(layer_box(Rect(far_resh.x + 30, far_resh.y + 50,
                                far_resh.w - 60, 28),
                           "Reshape"))

    # Concat node
    c_cy = 555
    c_r = 18
    c_node = Circle(W / 2, c_cy, c_r)
    parts.append(circle_node(c_node, "white", "#444",
                             label="C", font_size=14, sw=1.4,
                             label_dy=-2))
    tab_bottom = (col_xs[0], sub_rects[0].y + sub_h)
    img_bottom = (col_xs[1], img_resh.y + img_resh.h)
    far_bottom = (col_xs[2], far_resh.y + far_resh.h)
    for bx, by in (tab_bottom, img_bottom, far_bottom):
        parts.append(arrow_svg(bx, by, c_node.cx, c_node.cy - c_r,
                               stroke=EDGE, sw=1.2,
                               head_len=8, head_w=4))

    # Output layer
    out_y = 610
    out_w, out_h = 270, 175
    out = Rect(W / 2 - out_w / 2, out_y, out_w, out_h)
    parts.append(vertical_arrow(W / 2, c_cy + c_r, out.y - 1))
    parts.append(dashed_box(out, "Output Layer",
                            HDR_ORANGE, HDR_ORANGE_LN))
    parts.append(layer_box(Rect(out.x + 50, out.y + 36,
                                out.w - 100, 28), "Dense Layer"))
    parts.append(text(out.cx, out.y + 88, "⋮", size=18))
    parts.append(layer_box(Rect(out.x + 50, out.y + 100,
                                out.w - 100, 28), "Dense Layer"))
    parts.append(layer_box(Rect(out.x + 50, out.y + 140,
                                out.w - 100, 28), "Activation"))

    # Transformation + Survival Function
    trafo_w, trafo_h = 150, 95
    trafo = Rect(W / 2 - trafo_w / 2, out.y + out.h + 20,
                 trafo_w, trafo_h)
    parts.append(vertical_arrow(W / 2, out.y + out.h, trafo.y - 1))
    parts.append(
        f'<rect x="{trafo.x}" y="{trafo.y}" width="{trafo.w}" '
        f'height="{trafo.h}" rx="4" ry="4" fill="white" '
        f'stroke="#666" stroke-width="1.2" stroke-dasharray="{DASH}"/>'
    )
    parts.append(text(trafo.cx, trafo.y + 18, "Transformation*",
                      size=12, weight="bold", fill="#222"))
    parts.append(circle_node(Circle(trafo.cx, trafo.y + 58, 20),
                             "white", "#444",
                             label="ψ(·)", font_size=13, sw=1.2,
                             label_dy=-2))

    surv_w, surv_h = 150, 55
    surv = Rect(trafo.x + trafo.w + 60, trafo.cy - surv_h / 2,
                surv_w, surv_h)
    parts.append(
        f'<rect x="{surv.x}" y="{surv.y}" width="{surv.w}" '
        f'height="{surv.h}" rx="4" ry="4" '
        f'fill="{HDR_GREEN}" stroke="{HDR_GREEN_LN}" stroke-width="1.4"/>'
    )
    parts.append(text(surv.cx, surv.cy - 4, "Survival",
                      size=13, weight="bold", fill="#222"))
    parts.append(text(surv.cx, surv.cy + 14, "Function",
                      size=13, weight="bold", fill="#222"))
    parts.append(arrow_svg(trafo.x + trafo.w, trafo.cy,
                           surv.x, surv.cy,
                           stroke=EDGE, sw=1.4,
                           head_len=9, head_w=4))

    svg_path = write_svg("multimodal-fusion", W, H, "".join(parts))
    rasterise(svg_path)
    return svg_path

def main():
    for p in (fig_ffnn(), fig_shared_trunk_cr(), fig_multimodal_fusion()):
        print(p)


if __name__ == "__main__":
    main()
