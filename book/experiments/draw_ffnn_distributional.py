"""Emit /tmp/ffnn-distributional.html (auto-refresh preview)."""
from pathlib import Path

CANVAS_W, CANVAS_H = 1140, 780
PANEL_H = 240
PANEL_GAP = 30
PANEL_TOP = [20, 20 + PANEL_H + PANEL_GAP, 20 + 2 * (PANEL_H + PANEL_GAP)]

INPUT_FILL, INPUT_STROKE = "#E8EEF7", "#5B7CA8"
SHARED_FILL, SHARED_STROKE = "#FFF1D6", "#C99A2C"
MU_FILL, MU_STROKE = "#FADBD8", "#B03A2E"
SI_FILL, SI_STROKE = "#D6E4F0", "#2C5582"
Z_FILL, Z_STROKE = "#FFFFFF", "#555555"
R_FILL, R_STROKE = "#FFF7C2", "#C4A000"
OUT_FILL, OUT_STROKE = "#E8F5E9", "#3F8A52"
DIST_FILL, DIST_STROKE = "#F3E5F5", "#6A1B9A"
EDGE = "#444"

X_INPUT = 70
X_SHARED_L, X_SHARED_R = 130, 280       # variant (a) shared hidden
X_SHARED_NARROW_R = 250                  # variant (b) shared hidden ends sooner
X_SUBNET_L, X_SUBNET_R = 295, 415        # variant (b) sub-nets
X_SEPNET_L, X_SEPNET_R = 130, 320        # variant (c) separate networks
X_Z = 470
X_R_L, X_R_R = 510, 555
X_MU = 600
X_DIST_L, X_DIST_R = 700, 920

Y_HEAD_OFFSET = 32

def arrow(x1, y1, x2, y2, head_len=10, head_w=5, color=EDGE, sw=1.4):
    dx, dy = x2 - x1, y2 - y1
    L = (dx*dx + dy*dy) ** 0.5
    if L == 0:
        return ""
    ux, uy = dx / L, dy / L
    px, py = -uy, ux
    bx, by = x2 - head_len * ux, y2 - head_len * uy
    b1x, b1y = bx + head_w * px, by + head_w * py
    b2x, b2y = bx - head_w * px, by - head_w * py
    return (
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{bx:.1f}" y2="{by:.1f}" '
        f'stroke="{color}" stroke-width="{sw}"/>'
        f'<polygon points="{x2:.1f},{y2:.1f} {b1x:.1f},{b1y:.1f} {b2x:.1f},{b2y:.1f}" '
        f'fill="{color}"/>'
    )

def circle(cx, cy, r, fill, stroke, sw=1.6):
    return f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'

def rect(x, y, w, h, fill, stroke, rx=10, sw=1.6):
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" ry="{rx}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')

def text(x, y, s, size=15, weight="normal", anchor="middle", fill="#111"):
    return (f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}" '
            f'font-weight="{weight}" fill="{fill}">{s}</text>')

def label(x, y, s, size=18):
    return text(x, y, s, size=size, weight="bold", anchor="start")

def sub(base, sub_chars):
    """μ̂ / σ̂ with subscript trick using <tspan>."""
    return f'{base}<tspan dy="6" font-size="0.78em">{sub_chars}</tspan><tspan dy="-6"></tspan>'

def panel_a(parts, top):
    cy = top + PANEL_H // 2
    cy_mu = cy - Y_HEAD_OFFSET
    cy_si = cy + Y_HEAD_OFFSET
    parts.append(label(40, top + 24, "(a) shared hidden layers"))
    # input
    parts.append(circle(X_INPUT, cy, 22, INPUT_FILL, INPUT_STROKE))
    parts.append(text(X_INPUT, cy + 6, "x", size=20))
    # shared hidden block
    parts.append(rect(X_SHARED_L, cy - 50, X_SHARED_R - X_SHARED_L, 100, SHARED_FILL, SHARED_STROKE))
    parts.append(text((X_SHARED_L + X_SHARED_R)/2, cy - 6, "shared", size=15, weight="bold"))
    parts.append(text((X_SHARED_L + X_SHARED_R)/2, cy + 14, "hidden layers", size=14))
    # arrow x → shared
    parts.append(arrow(X_INPUT + 22, cy, X_SHARED_L, cy))
    # z_μ, z_σ
    parts.append(circle(X_Z, cy_mu, 18, Z_FILL, Z_STROKE))
    parts.append(text(X_Z, cy_mu + 5, 'z<tspan dy="4" font-size="0.72em">μ</tspan>', size=15))
    parts.append(circle(X_Z, cy_si, 18, Z_FILL, Z_STROKE))
    parts.append(text(X_Z, cy_si + 5, 'z<tspan dy="4" font-size="0.72em">σ</tspan>', size=15))
    # arrows shared → z_μ and z_σ
    parts.append(arrow(X_SHARED_R, cy - 20, X_Z - 18, cy_mu))
    parts.append(arrow(X_SHARED_R, cy + 20, X_Z - 18, cy_si))
    # r boxes
    parts.append(rect(X_R_L, cy_mu - 16, X_R_R - X_R_L, 32, R_FILL, R_STROKE, rx=8))
    parts.append(text((X_R_L + X_R_R)/2, cy_mu + 5, 'r<tspan dy="4" font-size="0.72em">μ</tspan>', size=15))
    parts.append(rect(X_R_L, cy_si - 16, X_R_R - X_R_L, 32, R_FILL, R_STROKE, rx=8))
    parts.append(text((X_R_L + X_R_R)/2, cy_si + 5, 'r<tspan dy="4" font-size="0.72em">σ</tspan>', size=15))
    parts.append(arrow(X_Z + 18, cy_mu, X_R_L, cy_mu))
    parts.append(arrow(X_Z + 18, cy_si, X_R_L, cy_si))
    # μ(x), σ(x)
    parts.append(circle(X_MU, cy_mu, 22, OUT_FILL, OUT_STROKE))
    parts.append(text(X_MU, cy_mu + 5, "μ(x)", size=14))
    parts.append(circle(X_MU, cy_si, 22, OUT_FILL, OUT_STROKE))
    parts.append(text(X_MU, cy_si + 5, "σ(x)", size=14))
    parts.append(arrow(X_R_R, cy_mu, X_MU - 22, cy_mu))
    parts.append(arrow(X_R_R, cy_si, X_MU - 22, cy_si))
    # distribution box
    parts.append(rect(X_DIST_L, cy - 32, X_DIST_R - X_DIST_L, 64, DIST_FILL, DIST_STROKE))
    parts.append(text((X_DIST_L + X_DIST_R)/2, cy + 6,
                      "N(μ(x), σ(x)²)", size=17))
    # arrows μ,σ → dist
    parts.append(arrow(X_MU + 22, cy_mu, X_DIST_L, cy - 14))
    parts.append(arrow(X_MU + 22, cy_si, X_DIST_L, cy + 14))

def panel_b(parts, top):
    cy = top + PANEL_H // 2
    cy_mu = cy - Y_HEAD_OFFSET
    cy_si = cy + Y_HEAD_OFFSET
    parts.append(label(40, top + 24, "(b) shared hidden layers + per-parameter sub-networks"))
    # input
    parts.append(circle(X_INPUT, cy, 22, INPUT_FILL, INPUT_STROKE))
    parts.append(text(X_INPUT, cy + 6, "x", size=20))
    # shared hidden block (narrower)
    parts.append(rect(X_SHARED_L, cy - 50, X_SHARED_NARROW_R - X_SHARED_L, 100, SHARED_FILL, SHARED_STROKE))
    parts.append(text((X_SHARED_L + X_SHARED_NARROW_R)/2, cy - 6, "shared", size=15, weight="bold"))
    parts.append(text((X_SHARED_L + X_SHARED_NARROW_R)/2, cy + 14, "hidden", size=14))
    parts.append(arrow(X_INPUT + 22, cy, X_SHARED_L, cy))
    # sub-net μ and σ
    parts.append(rect(X_SUBNET_L, cy_mu - 25, X_SUBNET_R - X_SUBNET_L, 50, MU_FILL, MU_STROKE))
    parts.append(text((X_SUBNET_L + X_SUBNET_R)/2, cy_mu + 5, "sub-net for μ", size=13))
    parts.append(rect(X_SUBNET_L, cy_si - 25, X_SUBNET_R - X_SUBNET_L, 50, SI_FILL, SI_STROKE))
    parts.append(text((X_SUBNET_L + X_SUBNET_R)/2, cy_si + 5, "sub-net for σ", size=13))
    # arrows shared → sub-nets
    parts.append(arrow(X_SHARED_NARROW_R, cy - 10, X_SUBNET_L, cy_mu))
    parts.append(arrow(X_SHARED_NARROW_R, cy + 10, X_SUBNET_L, cy_si))
    # z
    parts.append(circle(X_Z, cy_mu, 18, Z_FILL, Z_STROKE))
    parts.append(text(X_Z, cy_mu + 5, 'z<tspan dy="4" font-size="0.72em">μ</tspan>', size=15))
    parts.append(circle(X_Z, cy_si, 18, Z_FILL, Z_STROKE))
    parts.append(text(X_Z, cy_si + 5, 'z<tspan dy="4" font-size="0.72em">σ</tspan>', size=15))
    parts.append(arrow(X_SUBNET_R, cy_mu, X_Z - 18, cy_mu))
    parts.append(arrow(X_SUBNET_R, cy_si, X_Z - 18, cy_si))
    # r boxes
    parts.append(rect(X_R_L, cy_mu - 16, X_R_R - X_R_L, 32, R_FILL, R_STROKE, rx=8))
    parts.append(text((X_R_L + X_R_R)/2, cy_mu + 5, 'r<tspan dy="4" font-size="0.72em">μ</tspan>', size=15))
    parts.append(rect(X_R_L, cy_si - 16, X_R_R - X_R_L, 32, R_FILL, R_STROKE, rx=8))
    parts.append(text((X_R_L + X_R_R)/2, cy_si + 5, 'r<tspan dy="4" font-size="0.72em">σ</tspan>', size=15))
    parts.append(arrow(X_Z + 18, cy_mu, X_R_L, cy_mu))
    parts.append(arrow(X_Z + 18, cy_si, X_R_L, cy_si))
    # μ(x), σ(x)
    parts.append(circle(X_MU, cy_mu, 22, OUT_FILL, OUT_STROKE))
    parts.append(text(X_MU, cy_mu + 5, "μ(x)", size=14))
    parts.append(circle(X_MU, cy_si, 22, OUT_FILL, OUT_STROKE))
    parts.append(text(X_MU, cy_si + 5, "σ(x)", size=14))
    parts.append(arrow(X_R_R, cy_mu, X_MU - 22, cy_mu))
    parts.append(arrow(X_R_R, cy_si, X_MU - 22, cy_si))
    # distribution box
    parts.append(rect(X_DIST_L, cy - 32, X_DIST_R - X_DIST_L, 64, DIST_FILL, DIST_STROKE))
    parts.append(text((X_DIST_L + X_DIST_R)/2, cy + 6,
                      "N(μ(x), σ(x)²)", size=17))
    parts.append(arrow(X_MU + 22, cy_mu, X_DIST_L, cy - 14))
    parts.append(arrow(X_MU + 22, cy_si, X_DIST_L, cy + 14))

def panel_c(parts, top):
    cy = top + PANEL_H // 2
    cy_mu = cy - Y_HEAD_OFFSET - 10
    cy_si = cy + Y_HEAD_OFFSET + 10
    parts.append(label(40, top + 24, "(c) separate networks"))
    # input x (single, but feeds both nets)
    parts.append(circle(X_INPUT, cy, 22, INPUT_FILL, INPUT_STROKE))
    parts.append(text(X_INPUT, cy + 6, "x", size=20))
    # network for μ and σ
    parts.append(rect(X_SEPNET_L, cy_mu - 28, X_SEPNET_R - X_SEPNET_L, 56, MU_FILL, MU_STROKE))
    parts.append(text((X_SEPNET_L + X_SEPNET_R)/2, cy_mu + 5, "network for μ", size=14, weight="bold"))
    parts.append(rect(X_SEPNET_L, cy_si - 28, X_SEPNET_R - X_SEPNET_L, 56, SI_FILL, SI_STROKE))
    parts.append(text((X_SEPNET_L + X_SEPNET_R)/2, cy_si + 5, "network for σ", size=14, weight="bold"))
    # arrows x → both networks
    parts.append(arrow(X_INPUT + 22, cy - 8, X_SEPNET_L, cy_mu))
    parts.append(arrow(X_INPUT + 22, cy + 8, X_SEPNET_L, cy_si))
    # z
    parts.append(circle(X_Z, cy_mu, 18, Z_FILL, Z_STROKE))
    parts.append(text(X_Z, cy_mu + 5, 'z<tspan dy="4" font-size="0.72em">μ</tspan>', size=15))
    parts.append(circle(X_Z, cy_si, 18, Z_FILL, Z_STROKE))
    parts.append(text(X_Z, cy_si + 5, 'z<tspan dy="4" font-size="0.72em">σ</tspan>', size=15))
    parts.append(arrow(X_SEPNET_R, cy_mu, X_Z - 18, cy_mu))
    parts.append(arrow(X_SEPNET_R, cy_si, X_Z - 18, cy_si))
    # r boxes
    parts.append(rect(X_R_L, cy_mu - 16, X_R_R - X_R_L, 32, R_FILL, R_STROKE, rx=8))
    parts.append(text((X_R_L + X_R_R)/2, cy_mu + 5, 'r<tspan dy="4" font-size="0.72em">μ</tspan>', size=15))
    parts.append(rect(X_R_L, cy_si - 16, X_R_R - X_R_L, 32, R_FILL, R_STROKE, rx=8))
    parts.append(text((X_R_L + X_R_R)/2, cy_si + 5, 'r<tspan dy="4" font-size="0.72em">σ</tspan>', size=15))
    parts.append(arrow(X_Z + 18, cy_mu, X_R_L, cy_mu))
    parts.append(arrow(X_Z + 18, cy_si, X_R_L, cy_si))
    # μ(x), σ(x)
    parts.append(circle(X_MU, cy_mu, 22, OUT_FILL, OUT_STROKE))
    parts.append(text(X_MU, cy_mu + 5, "μ(x)", size=14))
    parts.append(circle(X_MU, cy_si, 22, OUT_FILL, OUT_STROKE))
    parts.append(text(X_MU, cy_si + 5, "σ(x)", size=14))
    parts.append(arrow(X_R_R, cy_mu, X_MU - 22, cy_mu))
    parts.append(arrow(X_R_R, cy_si, X_MU - 22, cy_si))
    # distribution box
    parts.append(rect(X_DIST_L, cy - 32, X_DIST_R - X_DIST_L, 64, DIST_FILL, DIST_STROKE))
    parts.append(text((X_DIST_L + X_DIST_R)/2, cy + 6,
                      "N(μ(x), σ(x)²)", size=17))
    parts.append(arrow(X_MU + 22, cy_mu, X_DIST_L, cy - 14))
    parts.append(arrow(X_MU + 22, cy_si, X_DIST_L, cy + 14))

parts = []
panel_a(parts, PANEL_TOP[0])
panel_b(parts, PANEL_TOP[1])
panel_c(parts, PANEL_TOP[2])

# separators
for sep_y in (PANEL_TOP[0] + PANEL_H + PANEL_GAP // 2,
              PANEL_TOP[1] + PANEL_H + PANEL_GAP // 2):
    parts.append(f'<line x1="40" y1="{sep_y}" x2="{CANVAS_W-40}" y2="{sep_y}" '
                 f'stroke="#bbb" stroke-width="1.2" stroke-dasharray="6,6"/>')

svg = (
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{CANVAS_W}" height="{CANVAS_H}" '
    f'viewBox="0 0 {CANVAS_W} {CANVAS_H}" font-family="DejaVu Sans">'
    + "".join(parts)
    + '</svg>'
)

html = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <meta http-equiv="refresh" content="2">
  <title>FFNN distributional regression — three architectural variants</title>
  <style>
    body {{ margin: 0; display: flex; justify-content: center; align-items: center; height: 100vh; background: transparent; font-family: 'DejaVu Sans', sans-serif; }}
    svg {{ max-width: 95vw; max-height: 90vh; }}
  </style>
</head>
<body>
{svg}
</body>
</html>
"""
Path("/tmp/ffnn-distributional.html").write_text(html)
print("wrote /tmp/ffnn-distributional.html")
