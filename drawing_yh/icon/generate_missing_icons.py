# -*- coding: utf-8 -*-
"""
Generate 7 missing tissue icons for the bulk RNA-seq B_cell_Aging project:
    cerebellum, frontal, temporal, occipital, lymphnode, testis, duodenum

Each icon is defined as a list of (opacity, SVG path-d) subpaths. We write:
    lib/<name>.svg   — one <path fill="currentColor" fill-opacity=...> per subpath
    lib/<name>.png   — matplotlib-rasterized version (parse same d-string)

Keeping both artefacts in sync via a single source of truth (the path-d strings).

Why this approach:
- Consistent with the library's `change_svg_color` util: it walks <path fill=...>,
  so plain filled paths (no stroke) with `currentColor` play nicely.
- Style roughly matches existing heterogeneous library (brain.svg, liver.svg, etc.):
  single-color iconic silhouettes in 256-viewBox.
- No cairo required (pycairo/cairocffi/rsvg all fail on this box).

Re-run from the repo root:
    python drawing_yh/icon/generate_missing_icons.py
"""
import re
from pathlib import Path as PathLib

import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch

LIB = PathLib(__file__).parent / "lib"
VB = 256  # viewBox side

# ---------------------------------------------------------------------------
# SVG-d parser (absolute M/L/H/V/C/Z only — enough for our hand-drawn shapes)
# ---------------------------------------------------------------------------
_TOKEN_RE = re.compile(r"[MLHVCZmlhvcz]|-?\d+\.?\d*(?:e[+-]?\d+)?")


def parse_d(d: str):
    toks = _TOKEN_RE.findall(d)
    verts, codes = [], []
    cur = (0.0, 0.0)
    start = cur
    i = 0
    while i < len(toks):
        t = toks[i]
        i += 1
        if t in "Mm":
            x, y = float(toks[i]), float(toks[i + 1])
            i += 2
            verts.append((x, y))
            codes.append(MPath.MOVETO)
            cur = (x, y)
            start = cur
            while i < len(toks) and not toks[i][0].isalpha():
                x, y = float(toks[i]), float(toks[i + 1])
                i += 2
                verts.append((x, y))
                codes.append(MPath.LINETO)
                cur = (x, y)
        elif t in "Ll":
            while i < len(toks) and not toks[i][0].isalpha():
                x, y = float(toks[i]), float(toks[i + 1])
                i += 2
                verts.append((x, y))
                codes.append(MPath.LINETO)
                cur = (x, y)
        elif t in "Cc":
            while i < len(toks) and not toks[i][0].isalpha():
                c1 = (float(toks[i]), float(toks[i + 1]))
                c2 = (float(toks[i + 2]), float(toks[i + 3]))
                pt = (float(toks[i + 4]), float(toks[i + 5]))
                i += 6
                verts += [c1, c2, pt]
                codes += [MPath.CURVE4] * 3
                cur = pt
        elif t in "Hh":
            while i < len(toks) and not toks[i][0].isalpha():
                x = float(toks[i])
                i += 1
                verts.append((x, cur[1]))
                codes.append(MPath.LINETO)
                cur = (x, cur[1])
        elif t in "Vv":
            while i < len(toks) and not toks[i][0].isalpha():
                y = float(toks[i])
                i += 1
                verts.append((cur[0], y))
                codes.append(MPath.LINETO)
                cur = (cur[0], y)
        elif t in "Zz":
            verts.append(start)
            codes.append(MPath.CLOSEPOLY)
            cur = start
    return verts, codes


# ---------------------------------------------------------------------------
# Icon definitions — list of (opacity, d) subpaths per icon
# 256×256 viewBox, face-left orientation for brain/lobe icons
# ---------------------------------------------------------------------------
# Lateral-brain silhouette (reused by the three lobe icons as 30%-opacity base)
_BRAIN = (
    "M 40 130 C 40 80 75 50 130 50 "
    "C 180 50 215 75 218 125 "
    "C 220 165 205 190 175 195 "
    "C 145 200 115 200 85 198 "
    "C 55 195 30 180 35 155 "
    "C 32 140 35 135 40 130 Z"
)

ICONS = {
    # bilobed pinched shape (top+bottom notch = vermis pinch)
    "cerebellum": [
        (1.0,
         "M 40 128 C 40 92 80 74 128 86 "
         "C 176 74 216 92 216 128 "
         "C 216 164 176 182 128 170 "
         "C 80 182 40 164 40 128 Z"),
    ],

    # lateral brain + front portion filled solid
    "frontal": [
        (0.3, _BRAIN),
        (1.0,
         "M 40 130 C 40 80 75 50 125 52 "
         "L 118 145 L 60 148 "
         "C 45 145 40 138 40 130 Z"),
    ],

    # lateral brain + lower-middle (temple) filled solid
    "temporal": [
        (0.3, _BRAIN),
        (1.0,
         "M 60 148 L 175 150 "
         "C 172 182 155 200 125 200 "
         "C 90 200 66 180 60 148 Z"),
    ],

    # lateral brain + back-pole filled solid
    "occipital": [
        (0.3, _BRAIN),
        (1.0,
         "M 180 55 C 205 65 218 90 218 125 "
         "C 220 160 208 185 175 195 "
         "L 170 150 "
         "C 185 130 185 90 180 55 Z"),
    ],

    # bean + 3 afferent vessels (left) + 1 efferent (right)
    "lymphnode": [
        (1.0, "M 15 92 L 50 92 L 50 98 L 15 98 Z"),
        (1.0, "M 15 130 L 50 130 L 50 136 L 15 136 Z"),
        (1.0, "M 15 170 L 50 170 L 50 176 L 15 176 Z"),
        (1.0, "M 220 128 L 246 128 L 246 134 L 220 134 Z"),
        (1.0,
         "M 60 105 C 50 70 100 50 150 65 "
         "C 195 75 220 105 220 140 "
         "C 225 175 200 200 160 205 "
         "C 115 210 65 200 58 170 "
         "C 52 145 55 125 60 105 Z"),
    ],

    # testis body + epididymis tube on top + vas deferens stub
    "testis": [
        # epididymis (drawn first, partially tucked behind body)
        (1.0,
         "M 95 75 C 110 55 170 48 200 65 "
         "C 210 75 205 90 185 88 "
         "C 145 82 115 90 95 95 Z"),
        # vas deferens
        (1.0, "M 198 70 L 222 45 L 228 50 L 204 75 Z"),
        # body
        (1.0,
         "M 58 170 C 58 125 95 95 135 95 "
         "C 180 95 205 135 205 170 "
         "C 205 215 175 230 130 230 "
         "C 85 230 58 215 58 170 Z"),
    ],

    # C-shaped duodenal loop + small pylorus nub
    "duodenum": [
        (1.0,
         "M 70 55 C 160 55 215 90 215 135 "
         "C 215 180 160 215 70 215 "
         "L 70 195 "
         "C 145 195 195 165 195 135 "
         "C 195 105 145 75 70 75 Z"),
        (1.0, "M 44 55 L 70 55 L 70 75 L 44 75 Z"),
    ],
}


# ---------------------------------------------------------------------------
def render(name: str, parts):
    # SVG
    body = "".join(
        f'<path fill="currentColor" fill-opacity="{op}" d="{d}"/>'
        for op, d in parts
    )
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="200" height="200" viewBox="0 0 {VB} {VB}">'
        f"{body}</svg>"
    )
    (LIB / f"{name}.svg").write_text(svg, encoding="utf-8")

    # PNG via matplotlib
    fig, ax = plt.subplots(figsize=(2, 2))
    for op, d in parts:
        verts, codes = parse_d(d)
        ax.add_patch(
            PathPatch(
                MPath(verts, codes),
                facecolor="black",
                alpha=op,
                edgecolor="none",
            )
        )
    ax.set_xlim(0, VB)
    ax.set_ylim(0, VB)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(
        LIB / f"{name}.png",
        dpi=150,
        transparent=True,
        bbox_inches="tight",
        pad_inches=0,
    )
    plt.close(fig)
    print(f"  wrote {name}.svg + {name}.png")


def main():
    print(f"generating {len(ICONS)} icons → {LIB}")
    for name, parts in ICONS.items():
        render(name, parts)
    print("done.")


if __name__ == "__main__":
    main()
