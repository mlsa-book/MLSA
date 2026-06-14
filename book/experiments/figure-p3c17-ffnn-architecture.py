"""Generate the ``fig-p3c17-ffnn-architecture`` schematic for the neural
networks chapter.

The figure is written as a standalone SVG to the chapter figure folder
and then rasterised into a transparent-background PNG with ImageMagick
(``convert``) so that the book renders cleanly to both HTML and PDF.
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


def sub_sup_label(cx, cy, base, sub="", sup="", font_size=18,
                  small_ratio=0.62, fill="#111"):
    """Render ``base`` with a subscript and superscript at proper offsets.

    Unicode subscript/superscript glyphs (₁, ⁽¹⁾) don't sit low/high enough
    in DejaVu Sans to read as real sub/super-scripts when rasterised through
    ImageMagick, so we draw three separate text elements with explicit y
    offsets instead.  Returns SVG ``<text>`` markup.
    """
    small = font_size * small_ratio
    parts = []
    # Estimated half-width of the main glyph so the sub/super hang to its
    # right edge rather than to the centre of the unit.
    half_main = font_size * 0.32
    base_y = cy + font_size * 0.34
    parts.append(
        f'<text x="{cx - half_main * 0.4}" y="{base_y:.2f}" '
        f'text-anchor="middle" font-family="{FONT}" '
        f'font-size="{font_size}" fill="{fill}">{base}</text>'
    )
    if sub:
        parts.append(
            f'<text x="{cx + half_main + 1}" y="{base_y + small * 0.45:.2f}" '
            f'text-anchor="start" font-family="{FONT}" '
            f'font-size="{small:.1f}" fill="{fill}">{sub}</text>'
        )
    if sup:
        parts.append(
            f'<text x="{cx + half_main + 1}" '
            f'y="{base_y - font_size * 0.55:.2f}" '
            f'text-anchor="start" font-family="{FONT}" '
            f'font-size="{small:.1f}" fill="{fill}">{sup}</text>'
        )
    return "".join(parts)


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
    """Feed-forward NN diagram: input + two hidden layers + output.

    Layout and unit labels mirror the workshop slide
    ``mlcu-dl-intro-workshop/slides/day2-01.qmd`` (3 inputs, 5 + 3 hidden
    units).  Each layer also carries an explicit bias node labelled ``1``
    that connects to every unit in the next layer.  The final linear
    combination ``z`` is mapped to the output ``y`` through the response
    function ``r(·)``, shown as a yellow box (the link function is
    ``r⁻¹(·)``, consistent with the pseudo-value chapter notation).
    """
    W, H = 1500, 760
    parts = []

    x_in, x_h1, x_h2, x_z = 250, 540, 830, 1060
    r_node = 40
    spacing = 102
    top_band_y = 50
    centre = H / 2 + 35
    # All biases sit on the same horizontal line below the tallest column
    # for visual consistency across layers.
    bias_y = centre + (4 - 1) * spacing / 2 + spacing - 4

    def column(x, labels, bias_below=True):
        n = len(labels)
        y0 = centre - (n - 1) * spacing / 2
        ys = [y0 + i * spacing for i in range(n)]
        nodes = [Circle(x, y, r_node) for y in ys]
        bias = Circle(x, bias_y, r_node - 4) if bias_below else None
        return nodes, bias

    # Labels are stored as (base, subscript, superscript) triples and
    # rendered with sub_sup_label so the indices sit at the correct height.
    in_labels = [("x", "1", ""), ("x", "2", "")]
    h1_labels = [("h", str(j + 1), "(1)") for j in range(4)]
    h2_labels = [("h", str(j + 1), "(2)") for j in range(3)]

    inputs, b_in = column(x_in, in_labels)
    h1, b_h1 = column(x_h1, h1_labels)
    h2, b_h2 = column(x_h2, h2_labels)

    # bottom-to-top visual ordering (index 1 at the bottom)
    inputs.reverse()
    h1.reverse()
    h2.reverse()

    PAL_BIAS = ("#F5F5F5", "#555555")

    # z pre-activation node (sits centred at the height of the layer columns)
    z_node = Circle(x_z, centre, r_node)

    # response-function box (yellow) and output unit
    r_box = Rect(x_z + 88, centre - 36, 106, 72)
    x_out = r_box.x + r_box.w + 88
    out = Circle(x_out, centre, r_node + 2)

    # draw bias nodes first so connections render on top of their faces
    for c, (base, sub, sup) in zip(inputs, in_labels):
        parts.append(circle_node(c, *PAL_INPUT, sw=1.8))
        parts.append(sub_sup_label(c.cx, c.cy, base,
                                   sub=sub, sup=sup, font_size=24))
    parts.append(circle_node(b_in, *PAL_BIAS, label="1",
                             font_size=23, sw=1.6, label_dy=-2))
    for c, (base, sub, sup) in zip(h1, h1_labels):
        parts.append(circle_node(c, *PAL_HIDDEN, sw=1.6))
        parts.append(sub_sup_label(c.cx, c.cy, base,
                                   sub=sub, sup=sup, font_size=21))
    parts.append(circle_node(b_h1, *PAL_BIAS, label="1",
                             font_size=23, sw=1.6, label_dy=-2))
    for c, (base, sub, sup) in zip(h2, h2_labels):
        parts.append(circle_node(c, *PAL_HIDDEN, sw=1.6))
        parts.append(sub_sup_label(c.cx, c.cy, base,
                                   sub=sub, sup=sup, font_size=21))
    parts.append(circle_node(b_h2, *PAL_BIAS, label="1",
                             font_size=23, sw=1.6, label_dy=-2))

    # z, r-box and y
    parts.append(circle_node(z_node, "white", "#444", label="z",
                             font_size=24, sw=1.8, label_dy=-2))
    parts.append(
        f'<rect x="{r_box.x}" y="{r_box.y}" width="{r_box.w}" '
        f'height="{r_box.h}" rx="10" ry="10" '
        f'fill="#FFF7C2" stroke="#C4A000" stroke-width="1.8"/>'
    )
    parts.append(text(r_box.cx, r_box.cy + 8, "r(·)",
                      size=24, fill="#222"))
    parts.append(circle_node(out, *PAL_OUTPUT, label="y",
                             font_size=24, sw=2.2, label_dy=-2))

    # fan-out edges between layers (scaled to match the larger nodes)
    for a in inputs + [b_in]:
        for b in h1:
            parts.append(edge_between(a, b, sw=1.0, opacity=0.55,
                                      arrow=True, head_len=10, head_w=5))
    for a in h1 + [b_h1]:
        for b in h2:
            parts.append(edge_between(a, b, sw=1.0, opacity=0.55,
                                      arrow=True, head_len=10, head_w=5))
    for a in h2 + [b_h2]:
        parts.append(edge_between(a, z_node, sw=1.2, opacity=0.65,
                                  arrow=True, head_len=12, head_w=6))

    # z -> r-box -> y
    parts.append(arrow_svg(z_node.cx + z_node.r, z_node.cy,
                           r_box.x, r_box.cy,
                           stroke=EDGE, sw=2.0, head_len=14, head_w=7))
    parts.append(arrow_svg(r_box.x + r_box.w, r_box.cy,
                           out.cx - out.r, out.cy,
                           stroke=EDGE, sw=2.0, head_len=14, head_w=7))

    # Top band: bold layer caption + formula for the values in that column.
    parts.append(text(x_in, top_band_y, "input layer",
                      size=23, weight="bold"))
    parts.append(text(x_in, top_band_y + 30,
                      "x ∈ ℝ²",
                      size=23, fill="#444"))

    parts.append(text(x_h1, top_band_y, "hidden layer 1",
                      size=23, weight="bold"))
    parts.append(text(x_h1, top_band_y + 30,
                      "h⁽¹⁾ = a(W⁽¹⁾x + b⁽¹⁾)",
                      size=23, fill="#444"))

    parts.append(text(x_h2, top_band_y, "hidden layer 2",
                      size=23, weight="bold"))
    parts.append(text(x_h2, top_band_y + 30,
                      "h⁽²⁾ = a(W⁽²⁾h⁽¹⁾ + b⁽²⁾)",
                      size=23, fill="#444"))

    parts.append(text((z_node.cx + out.cx) / 2, top_band_y,
                      "output layer", size=23, weight="bold"))
    parts.append(text((z_node.cx + out.cx) / 2, top_band_y + 30,
                      "z = W⁽³⁾h⁽²⁾ + b⁽³⁾,   y = r(z)",
                      size=23, fill="#444"))

    svg_path = write_svg("fig-p3c17-ffnn-architecture", W, H, "".join(parts))
    rasterise(svg_path)
    return svg_path



def main():
    print(fig_ffnn())


if __name__ == "__main__":
    main()
