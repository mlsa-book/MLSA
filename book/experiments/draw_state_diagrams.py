"""Generate two state-diagram SVGs with proper arrowhead geometry.

Conventions
-----------
* Transient state → circle
* Absorbing state → square/rectangle
* Arrowhead tip touches the destination node's boundary.
* Line starts at the source node's boundary.
* Arrowhead base is perpendicular to the line direction.
"""
from __future__ import annotations
import math
from dataclasses import dataclass


# ---- node + arrow primitives --------------------------------------------

@dataclass
class Circle:
    cx: float
    cy: float
    r: float
    def boundary_toward(self, x: float, y: float) -> tuple[float, float]:
        dx, dy = x - self.cx, y - self.cy
        d = math.hypot(dx, dy)
        return (self.cx + self.r * dx / d, self.cy + self.r * dy / d)
    def svg(self, **kw) -> str:
        s = ' '.join(f'{k}="{v}"' for k, v in kw.items())
        return f'<circle cx="{self.cx}" cy="{self.cy}" r="{self.r}" {s}/>'


@dataclass
class Rect:
    x: float
    y: float
    w: float
    h: float
    @property
    def cx(self) -> float: return self.x + self.w / 2
    @property
    def cy(self) -> float: return self.y + self.h / 2
    def boundary_toward(self, x: float, y: float) -> tuple[float, float]:
        # ray from (x,y) toward center; find intersection with rect boundary
        dx, dy = x - self.cx, y - self.cy
        if dx == 0 and dy == 0:
            return (self.cx, self.cy)
        # parametric: pt = center + t*(dx,dy), find smallest t hitting rect edges
        ts = []
        if dx > 0: ts.append((self.x + self.w - self.cx) / dx)
        if dx < 0: ts.append((self.x - self.cx) / dx)
        if dy > 0: ts.append((self.y + self.h - self.cy) / dy)
        if dy < 0: ts.append((self.y - self.cy) / dy)
        t = min(t for t in ts if t > 0)
        return (self.cx + t * dx, self.cy + t * dy)
    def svg(self, **kw) -> str:
        s = ' '.join(f'{k}="{v}"' for k, v in kw.items())
        return f'<rect x="{self.x}" y="{self.y}" width="{self.w}" height="{self.h}" {s}/>'


def arrow(src, dst, perp_offset: float = 0,
          head_len: float = 12, head_w: float = 5,
          stroke_width: float = 1.5) -> str:
    """SVG for an arrow from `src` node to `dst` node.

    perp_offset shifts the whole arrow perpendicular to its direction
    (positive = "below" in screen coords for a left-to-right arrow);
    useful for bidirectional pairs.
    """
    # direction in world coords (source center -> dest center)
    dx0, dy0 = dst.cx - src.cx, dst.cy - src.cy
    L = math.hypot(dx0, dy0)
    ux, uy = dx0 / L, dy0 / L          # unit vector along arrow direction
    px, py = -uy, ux                    # rotate +90° (CCW)

    # perp offset applied to source-center and dest-center
    sx, sy = src.cx + perp_offset * px, src.cy + perp_offset * py
    tx, ty = dst.cx + perp_offset * px, dst.cy + perp_offset * py

    # boundary points using the shifted endpoints
    s_x, s_y = src.boundary_toward(tx, ty)
    s_x += perp_offset * px; s_y += perp_offset * py
    t_x, t_y = dst.boundary_toward(sx, sy)
    t_x += perp_offset * px; t_y += perp_offset * py

    # base of arrowhead
    bx, by = t_x - head_len * ux, t_y - head_len * uy
    blx, bly = bx + head_w * px, by + head_w * py
    brx, bry = bx - head_w * px, by - head_w * py

    def f(v): return f'{v:.2f}'
    line = (f'<line x1="{f(s_x)}" y1="{f(s_y)}" x2="{f(bx)}" y2="{f(by)}" '
            f'stroke="black" stroke-width="{stroke_width}"/>')
    head = (f'<polygon points="{f(t_x)},{f(t_y)} {f(blx)},{f(bly)} {f(brx)},{f(bry)}" '
            f'fill="black"/>')
    return line + '\n  ' + head


def midpoint_offset(src, dst, perp_offset: float, label_offset: float):
    """Return a point near the midpoint of an arrow, offset perpendicular by label_offset."""
    dx0, dy0 = dst.cx - src.cx, dst.cy - src.cy
    L = math.hypot(dx0, dy0)
    ux, uy = dx0 / L, dy0 / L
    px, py = -uy, ux
    mx = (src.cx + dst.cx) / 2 + perp_offset * px
    my = (src.cy + dst.cy) / 2 + perp_offset * py
    return (mx + label_offset * px, my + label_offset * py)


# ---- figure 1: illness-death model --------------------------------------

def illness_death_svg() -> str:
    W, H = 500, 400
    s0 = Circle(90, 200, 42)
    s1 = Circle(400, 90, 42)
    s2 = Rect(358, 268, 84, 84)

    parts = []

    # arrows
    parts.append(arrow(s0, s1))
    parts.append(arrow(s0, s2))
    parts.append(arrow(s1, s2))

    # nodes
    parts.append(s0.svg(fill="white", stroke="black", **{"stroke-width": "1.5"}))
    parts.append('<text x="90" y="208" text-anchor="middle" font-size="24">0</text>')
    parts.append(s1.svg(fill="white", stroke="black", **{"stroke-width": "1.5"}))
    parts.append('<text x="400" y="98" text-anchor="middle" font-size="24">1</text>')
    parts.append(s2.svg(fill="white", stroke="black", **{"stroke-width": "1.5"}))
    parts.append('<text x="400" y="318" text-anchor="middle" font-size="24">2</text>')

    # labels at midpoints, offset perpendicular
    def hlabel(src, dst, name, label_offset=-14):
        mx, my = midpoint_offset(src, dst, perp_offset=0, label_offset=label_offset)
        sub1, sub2 = name
        return (
            f'<text x="{mx-10:.1f}" y="{my+5:.1f}" font-size="18" font-style="italic">h</text>'
            f'<text x="{mx+2:.1f}" y="{my+11:.1f}" font-size="13">{sub1}{sub2}</text>'
            f'<text x="{mx+18:.1f}" y="{my+5:.1f}" font-size="18" font-style="italic">(τ)</text>'
        )

    parts.append(hlabel(s0, s1, ("0", "1"), label_offset=-30))
    parts.append(hlabel(s0, s2, ("0", "2"), label_offset=+34))
    parts.append(hlabel(s1, s2, ("1", "2"), label_offset=-22))

    body = '\n  '.join(parts)
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}"
     font-family="Helvetica,Arial,sans-serif">
  {body}
</svg>
'''


# ---- figure 2: multi-state-prothr-states --------------------------------

def prothr_svg() -> str:
    W, H = 600, 400
    s0 = Circle(90, 200, 78)
    s1 = Circle(430, 90, 78)
    s2 = Rect(370, 240, 120, 120)

    parts = []

    # bidirectional pair (0 <-> 1) with perpendicular offset for separation
    parts.append(arrow(s0, s1, perp_offset=-12))   # 0 -> 1 above
    parts.append(arrow(s1, s0, perp_offset=-12))   # 1 -> 0 also offset, but reversed dir

    # unidirectional 0 -> 2, 1 -> 2
    parts.append(arrow(s0, s2))
    parts.append(arrow(s1, s2))

    # nodes
    parts.append(s0.svg(fill="white", stroke="black", **{"stroke-width": "1.5"}))
    parts.append('<text x="90" y="188" text-anchor="middle" font-size="14">Normal</text>')
    parts.append('<text x="90" y="207" text-anchor="middle" font-size="14">Prothrombin</text>')
    parts.append('<text x="90" y="226" text-anchor="middle" font-size="14">levels (0)</text>')

    parts.append(s1.svg(fill="white", stroke="black", **{"stroke-width": "1.5"}))
    parts.append('<text x="430" y="78"  text-anchor="middle" font-size="14">Abnormal</text>')
    parts.append('<text x="430" y="97"  text-anchor="middle" font-size="14">Prothrombin</text>')
    parts.append('<text x="430" y="116" text-anchor="middle" font-size="14">levels (1)</text>')

    parts.append(s2.svg(fill="white", stroke="black", **{"stroke-width": "1.5"}))
    parts.append('<text x="430" y="295" text-anchor="middle" font-size="16">Death</text>')
    parts.append('<text x="430" y="318" text-anchor="middle" font-size="16">(2)</text>')

    body = '\n  '.join(parts)
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}"
     font-family="Helvetica,Arial,sans-serif">
  {body}
</svg>
'''


# ---- entry point --------------------------------------------------------

import sys, os
out_dir = sys.argv[1] if len(sys.argv) > 1 else "/tmp/mlsa-issue112/book/Figures/survival"
os.makedirs(out_dir, exist_ok=True)
open(os.path.join(out_dir, "illness-death-model.svg"), "w").write(illness_death_svg())
open(os.path.join(out_dir, "multi-state-prothr-states.svg"), "w").write(prothr_svg())
print(f"Wrote SVGs to {out_dir}")
