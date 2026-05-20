"""Generate a small synthetic chest-CT-like image.

Used as the ``Image data`` icon in the multi-modal survival figure.
Produces a recognisable thorax-slice look (body outline, lungs, mediastinum,
spine, ribs) without using any patient data, so the resulting PNG can be
distributed under the same licence as the rest of the book.

Output: ``book/Figures/neuralnetworks/synthetic-ct.png`` (256 x 256, RGB).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFilter


def generate_ct(size: int = 256, seed: int = 7) -> Image.Image:
    rng = np.random.default_rng(seed)

    canvas = Image.new("L", (size, size), 0)
    draw = ImageDraw.Draw(canvas)

    margin = 18
    body_bbox = (margin, margin + 8, size - margin, size - margin)
    draw.ellipse(body_bbox, fill=140)

    cy = size // 2 - 6
    # Lungs: tall, slightly inward-tapering kidney-bean shapes
    # composed of overlapping ellipses (one main + smaller dilations).
    for side in (-1, 1):
        cx = size // 2 + side * 46
        # main lung body, elongated vertically
        draw.ellipse((cx - 26, cy - 60, cx + 26, cy + 60), fill=20)
        # upper apex tapering
        draw.ellipse((cx - 18, cy - 80, cx + 18, cy - 40), fill=20)
        # carve inner edge with a body-coloured ellipse so the lung
        # "cups" toward the mediastinum
        inner_cx = cx - side * 18
        draw.ellipse((inner_cx - 14, cy - 50, inner_cx + 14, cy + 40),
                     fill=140)

    # Mediastinum + heart, off-centre to the patient's left
    heart_cx, heart_cy = size // 2 - 10, cy + 4
    draw.ellipse((heart_cx - 24, heart_cy - 26, heart_cx + 24, heart_cy + 26),
                 fill=180)
    draw.ellipse((heart_cx - 12, heart_cy - 32, heart_cx + 12, heart_cy - 8),
                 fill=170)

    # Spine: at the posterior (bottom in axial slice)
    spine_cx, spine_cy = size // 2, cy + 72
    draw.ellipse((spine_cx - 14, spine_cy - 12, spine_cx + 14, spine_cy + 12),
                 fill=230)
    draw.ellipse((spine_cx - 9, spine_cy - 8, spine_cx + 9, spine_cy + 8),
                 fill=110)

    # Ribs around the body outline
    n_ribs = 7
    for i in range(n_ribs):
        angle = -math.pi / 2 + (i + 1) * math.pi / (n_ribs + 1)
        for side in (-1, 1):
            cx = size / 2 + side * (size / 2 - margin - 4) * math.cos(angle)
            cy_r = size / 2 + (size / 2 - margin - 4) * math.sin(angle) - 8
            draw.ellipse((cx - 5, cy_r - 7, cx + 5, cy_r + 7), fill=230)

    arr = np.array(canvas, dtype=np.float32)
    arr += rng.normal(0.0, 6.0, arr.shape)

    yy, xx = np.indices(arr.shape)
    cx_c, cy_c = size / 2, size / 2 + 4
    r = np.hypot(xx - cx_c, yy - cy_c)
    body_mask = (r < (size / 2 - margin - 2)).astype(np.float32)
    speckle = rng.normal(0.0, 7.0, arr.shape) * body_mask
    arr += speckle

    arr = np.clip(arr, 0, 255).astype(np.uint8)
    img = Image.fromarray(arr, mode="L").convert("RGB")
    img = img.filter(ImageFilter.GaussianBlur(radius=1.0))
    return img


def main() -> None:
    out_dir = Path(__file__).resolve().parent.parent / "Figures" / "neuralnetworks"
    out_dir.mkdir(parents=True, exist_ok=True)
    img = generate_ct()
    out = out_dir / "synthetic-ct.png"
    img.save(out, optimize=True)
    print(out)


if __name__ == "__main__":
    main()
