# -*- coding: utf-8 -*-
"""
Generate anatomically-detailed tissue icons for the B_cell_Aging project:
    cerebellum, frontal, temporal, occipital, lymphnode, testis,
    duodenum, muscle, fat (wat)

Each icon is a COMPOSITION of multiple SVG primitives (path / ellipse / circle /
rect) to capture real anatomical features without relying on fill-rule tricks.
We write:
    lib/<name>.svg   — one primitive per element, all with fill="currentColor"
    lib/<name>.png   — matplotlib-rasterized from the same primitive spec

Why multiple primitives:
- Single-path designs looked too abstract. Stacking ellipses, rects and paths
  lets us show foliations (cerebellum), vessel valves (lymph node), epididymis
  head+body+tail (testis), mucosal folds (duodenum), etc.
- All primitives use plain `fill` (no stroke), so the library's
  `change_svg_color` util still recolours every element correctly.

Re-run:
    python drawing_yh/icon/generate_missing_icons.py
"""
import math
import re
from pathlib import Path as PathLib

import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Ellipse, Rectangle, Circle

LIB = PathLib(__file__).parent / "lib"
VB = 256


# ---------------------------------------------------------------------------
# SVG path-d parser (absolute M/L/H/V/C/Q/Z; enough for hand-authored shapes)
# ---------------------------------------------------------------------------
_TOK = re.compile(r"[MLHVCQZmlhvcqz]|-?\d+\.?\d*(?:e[+-]?\d+)?")


def parse_d(d: str):
    toks = _TOK.findall(d)
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
        elif t in "Qq":
            while i < len(toks) and not toks[i][0].isalpha():
                c1 = (float(toks[i]), float(toks[i + 1]))
                pt = (float(toks[i + 2]), float(toks[i + 3]))
                i += 4
                verts += [c1, pt]
                codes += [MPath.CURVE3, MPath.CURVE3]
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
# Shared lateral-brain silhouette (face pointing LEFT, 256 viewBox)
# Used as the 30%-opacity base for frontal / temporal / occipital icons
# ---------------------------------------------------------------------------
# Frontal pole (left), parietal dome (top), occipital pole (bottom-right),
# temporal lobe bulge (lower left), slight brainstem notch (bottom-middle-back)
_BRAIN_LATERAL = (
    "M 30 128 "
    "C 28 85 58 45 108 42 "                # frontal/superior frontal
    "C 162 38 210 60 222 108 "             # up to parietal dome
    "C 232 148 228 185 202 200 "           # back-down to parietal-occipital
    "C 185 212 168 213 152 208 "           # occipital pole
    "L 150 202 "
    "C 145 204 138 206 132 205 "
    "C 120 215 100 218 82 212 "            # temporal bulge descending
    "C 58 205 38 188 32 168 "
    "C 26 152 24 138 30 128 Z"
)

# Central-sulcus hint — thin crescent between frontal and parietal
_CENTRAL_SULCUS = (
    "M 132 45 "
    "C 126 75 121 110 118 142 "
    "C 121 141 124 140 127 140 "
    "C 130 108 135 76 140 47 Z"
)

# Sylvian-fissure hint — thin horizontal crescent in the middle
_SYLVIAN = (
    "M 55 152 "
    "C 95 148 140 146 185 142 "
    "C 185 146 185 149 184 152 "
    "C 140 155 95 157 55 158 Z"
)


# ---------------------------------------------------------------------------
# Icon specifications: list of primitives per icon
# Each primitive is a tuple: (kind, params_dict, opacity)
#   kind ∈ {"path", "ellipse", "rect", "circle"}
# ---------------------------------------------------------------------------
ICONS = {
    # Cerebellum: horizontal foliations (ellipses) stacked, tapering at poles,
    # giving the classic banded 'tree of life' appearance. Central vermis is
    # suggested by a very faint vertical band so it doesn't break the stripes.
    "cerebellum": [
        ("ellipse", dict(cx=128, cy=72,  rx=48, ry=9),  1.0),  # top foliation
        ("ellipse", dict(cx=128, cy=90,  rx=70, ry=9),  1.0),
        ("ellipse", dict(cx=128, cy=108, rx=82, ry=10), 1.0),
        ("ellipse", dict(cx=128, cy=127, rx=88, ry=10), 1.0),
        ("ellipse", dict(cx=128, cy=146, rx=88, ry=10), 1.0),  # widest (middle)
        ("ellipse", dict(cx=128, cy=165, rx=82, ry=10), 1.0),
        ("ellipse", dict(cx=128, cy=184, rx=68, ry=10), 1.0),
        ("ellipse", dict(cx=128, cy=202, rx=44, ry=8),  1.0),  # inferior
        # Faint vermis shading (central vertical band, barely visible)
        ("rect",    dict(x=123, y=64, w=10, h=146, rx=5), 0.18),
    ],

    # Lateral brain + FRONTAL lobe highlighted (front third)
    "frontal": [
        ("path",    dict(d=_BRAIN_LATERAL),   0.28),   # base silhouette
        ("path",    dict(d=_CENTRAL_SULCUS),  0.55),   # boundary cue
        ("path",    dict(d=_SYLVIAN),         0.45),
        # Frontal lobe fill
        ("path",    dict(d=(
            "M 30 128 "
            "C 28 85 58 45 108 42 "
            "C 124 41 130 42 134 44 "
            "L 128 95 L 122 148 "
            "C 100 152 72 154 52 152 "
            "C 40 150 30 140 30 128 Z"
        )), 1.0),
    ],

    # Lateral brain + TEMPORAL lobe highlighted.
    # Redesign: strong silhouette + single solid temporal bulge overlay.
    # The temporal lobe is the bulge below the Sylvian fissure (mid-lower side).
    "temporal": [
        # Full brain silhouette at 30% so the lobe pop is unambiguous
        ("path", dict(d=_BRAIN_LATERAL), 0.30),
        # Sylvian fissure hint (upper boundary of temporal)
        ("path", dict(d=_SYLVIAN), 0.55),
        # Temporal lobe fill — prominent bulge under the Sylvian fissure
        ("path", dict(d=(
            "M 52 158 "
            "C 92 160 140 158 186 154 "
            "L 184 172 "
            "C 178 192 156 206 128 208 "
            "C 100 210 76 202 58 190 "
            "C 40 176 38 162 52 158 Z"
        )), 1.0),
        # Small ear-like cue below to anchor "temporal = side of head"
        ("ellipse", dict(cx=52, cy=178, rx=8, ry=12), 1.0),
    ],

    # Lateral brain + OCCIPITAL lobe highlighted.
    # Redesign: full silhouette at low opacity, back pole filled dark,
    # plus an eye-cue on the opposite side makes the posterior orientation read.
    "occipital": [
        ("path", dict(d=_BRAIN_LATERAL), 0.30),
        ("path", dict(d=_SYLVIAN),       0.40),
        # Occipital lobe fill — back third of the brain (upper-back dome pole)
        ("path", dict(d=(
            "M 164 48 "
            "C 194 58 214 80 222 108 "
            "C 230 148 222 184 200 198 "
            "C 182 208 162 209 148 204 "
            "L 148 172 "
            "C 168 164 176 140 176 110 "
            "C 176 86 172 66 164 48 Z"
        )), 1.0),
    ],

    # Lymph node: simple kidney bean + a couple of afferent arrows and one
    # efferent. Smaller detail count → reads cleanly at 24 px.
    "lymphnode": [
        # 2 afferent vessels on the convex side (left)
        ("rect",    dict(x=14, y=102, w=56, h=8, rx=3), 1.0),
        ("rect",    dict(x=14, y=156, w=56, h=8, rx=3), 1.0),
        # Arrow-head valves on the afferents (pointing inward)
        ("path",    dict(d=(
            "M 62 94 L 80 106 L 62 118 Z"
        )), 1.0),
        ("path",    dict(d=(
            "M 62 148 L 80 160 L 62 172 Z"
        )), 1.0),
        # Efferent vessel exiting at the hilum (right side)
        ("rect",    dict(x=186, y=126, w=56, h=10, rx=4), 1.0),
        ("path",    dict(d=(
            "M 232 118 L 250 131 L 232 144 Z"
        )), 1.0),
        # Kidney-bean body with a clear hilum notch on the right
        ("path",    dict(d=(
            "M 82 96 "
            "C 72 68 116 50 152 58 "
            "C 190 66 212 90 214 124 "
            "C 204 124 198 130 198 138 "
            "C 198 148 206 152 212 150 "
            "C 206 180 186 204 150 210 "
            "C 110 216 80 208 70 184 "
            "C 60 156 68 126 82 96 Z"
        )), 1.0),
        # Inner-cortex hint (subtle lighter bean inside)
        ("path",    dict(d=(
            "M 102 112 "
            "C 96 92 128 82 156 90 "
            "C 182 98 194 114 194 132 "
            "C 188 136 186 142 188 148 "
            "C 184 174 168 190 146 194 "
            "C 118 198 96 186 90 168 "
            "C 84 148 92 130 102 112 Z"
        )), 0.22),
    ],

    # Muscle (flexed biceps): reads as "strong arm" at 24 px.
    # Composition: shoulder → upper arm with rounded bicep dome →
    # elbow blend → forearm rising vertically → fist on top.
    # Drawn as a single closed profile path so the silhouette is crisp.
    "muscle": [
        # Shoulder/deltoid blob anchoring the arm at lower-left
        ("ellipse", dict(cx=44, cy=186, rx=34, ry=38), 1.0),
        # Upper arm + flexed bicep (big rounded dome on top, horizontal body)
        ("path", dict(d=(
            "M 30 168 "                    # armpit (lower-left root)
            "C 30 138 58 110 100 104 "     # bicep underside rising
            "C 140 100 170 112 178 138 "   # bicep peak (top of dome)
            "C 182 160 170 180 150 190 "   # top-right of bicep sloping to elbow
            "C 130 200 100 204 72 202 "    # underside of bicep back to armpit
            "C 50 200 34 192 30 168 Z"
        )), 1.0),
        # Elbow joint (round blend between bicep and forearm)
        ("ellipse", dict(cx=166, cy=190, rx=26, ry=28), 1.0),
        # Forearm — vertical, going UP from the elbow toward the fist
        ("path", dict(d=(
            "M 144 60 "
            "L 190 60 "
            "C 198 100 196 144 186 184 "
            "C 182 200 170 206 156 202 "
            "C 142 198 140 184 144 164 "
            "C 150 132 148 96 144 60 Z"
        )), 1.0),
        # Wrist highlight
        ("ellipse", dict(cx=167, cy=58, rx=24, ry=10), 1.0),
        # Fist — large rounded blob at the top of the forearm
        ("path", dict(d=(
            "M 136 54 "
            "C 130 30 146 14 168 14 "
            "C 190 14 206 28 204 52 "
            "C 202 70 194 76 182 76 "
            "L 152 76 "
            "C 140 76 138 66 136 54 Z"
        )), 1.0),
        # Knuckle cuts (white lines on the fist)
        ("rect", dict(x=150, y=22, w=4, h=32, rx=1), 1.0, "white"),
        ("rect", dict(x=162, y=20, w=4, h=36, rx=1), 1.0, "white"),
        ("rect", dict(x=174, y=20, w=4, h=36, rx=1), 1.0, "white"),
        ("rect", dict(x=186, y=22, w=4, h=32, rx=1), 1.0, "white"),
        # Thin cut between bicep dome and shoulder (subtle muscle-definition line)
        ("path", dict(d=(
            "M 66 134 C 74 126 86 122 98 122 "
            "L 98 128 C 88 128 78 132 70 140 Z"
        )), 1.0, "white"),
    ],

    # Adipose tissue / white fat: cluster of 4 big round adipocytes, each a
    # thick dark ring with a bright hollow centre (classic histology look).
    # White-filled inner ellipses punch true holes thanks to the luminance-
    # aware tinting downstream.
    "fat": [
        # 4 large adipocytes in a 2×2 cluster, overlapping slightly
        ("ellipse", dict(cx=80,  cy=80,  rx=60, ry=60), 1.0),
        ("ellipse", dict(cx=80,  cy=80,  rx=46, ry=46), 1.0, "white"),
        ("ellipse", dict(cx=180, cy=84,  rx=58, ry=58), 1.0),
        ("ellipse", dict(cx=180, cy=84,  rx=44, ry=44), 1.0, "white"),
        ("ellipse", dict(cx=82,  cy=186, rx=60, ry=60), 1.0),
        ("ellipse", dict(cx=82,  cy=186, rx=46, ry=46), 1.0, "white"),
        ("ellipse", dict(cx=184, cy=190, rx=58, ry=58), 1.0),
        ("ellipse", dict(cx=184, cy=190, rx=44, ry=44), 1.0, "white"),
        # Small eccentric nucleus on each adipocyte (dark bean pushed to edge)
        ("ellipse", dict(cx=112, cy=50,  rx=8, ry=6), 1.0),
        ("ellipse", dict(cx=214, cy=54,  rx=8, ry=6), 1.0),
        ("ellipse", dict(cx=114, cy=156, rx=8, ry=6), 1.0),
        ("ellipse", dict(cx=216, cy=160, rx=8, ry=6), 1.0),
    ],

    # Testis: ovoid body drawn at 0.45 so the 1.0 epididymis + vas deferens pop
    "testis": [
        # Main ovoid body (testis proper) — lighter so epididymis shows
        ("ellipse", dict(cx=112, cy=160, rx=78, ry=82), 0.45),
        # Tunica albuginea hint — a slightly darker ring inside
        ("ellipse", dict(cx=112, cy=160, rx=66, ry=70), 0.22),
        # Epididymis HEAD — bulb at superior pole (upper right of testis)
        ("ellipse", dict(cx=170, cy=90, rx=30, ry=22), 1.0),
        # Epididymis BODY — substantial curved tube along posterior border
        ("path",    dict(d=(
            "M 192 100 "
            "C 214 128 216 172 204 212 "
            "C 200 228 186 236 174 230 "
            "C 166 222 174 210 182 196 "
            "C 194 160 192 128 178 108 "
            "C 174 100 184 95 192 100 Z"
        )), 1.0),
        # Epididymis TAIL curl at inferior-posterior pole
        ("ellipse", dict(cx=184, cy=228, rx=18, ry=13), 1.0),
        # Vas deferens — substantial tube ascending from tail superiorly
        ("path",    dict(d=(
            "M 168 18 "
            "L 188 18 "
            "C 192 52 186 76 180 90 "
            "L 162 90 "
            "C 158 76 162 52 168 18 Z"
        )), 1.0),
    ],

    # Duodenum: C-loop wall drawn at 0.45 opacity so darker mucosal folds show.
    # Pyloric cap + Ampulla of Vater + DJ flexure at full opacity as anatomical landmarks.
    "duodenum": [
        # C-shape wall ring — light enough that folds on top are visible
        ("path",    dict(d=(
            "M 70 40 "
            "C 160 37 222 78 222 130 "
            "C 222 185 160 225 70 222 "
            "L 70 198 "
            "C 148 200 200 170 200 130 "
            "C 200 95 148 62 70 64 Z"
        )), 0.45),
        # Pyloric cap / duodenal bulb (superior flexure) at top-left, full dark
        ("ellipse", dict(cx=52, cy=52, rx=22, ry=17), 1.0),
        # Ampulla of Vater — pronounced bulge on inner concave curve
        ("ellipse", dict(cx=196, cy=130, rx=13, ry=15), 1.0),
        # Duodenojejunal flexure — kink curving left at bottom, full dark
        ("path",    dict(d=(
            "M 70 222 "
            "C 50 226 32 218 30 200 "
            "C 28 188 36 180 48 184 "
            "C 56 188 64 194 70 198 Z"
        )), 1.0),
        # Transverse mucosal folds (plicae circulares) — dark bars crossing
        # the light wall, spaced along the C perimeter
        ("rect",    dict(x=114, y=38, w=6, h=30, rx=2), 1.0),
        ("rect",    dict(x=148, y=42, w=6, h=30, rx=2), 1.0),
        ("rect",    dict(x=178, y=58, w=6, h=30, rx=2), 1.0),
        ("rect",    dict(x=198, y=82, w=6, h=30, rx=2), 1.0),
        ("rect",    dict(x=206, y=118, w=6, h=26, rx=2), 1.0),
        ("rect",    dict(x=198, y=150, w=6, h=30, rx=2), 1.0),
        ("rect",    dict(x=178, y=176, w=6, h=30, rx=2), 1.0),
        ("rect",    dict(x=148, y=192, w=6, h=30, rx=2), 1.0),
        ("rect",    dict(x=114, y=196, w=6, h=30, rx=2), 1.0),
    ],
}


# ---------------------------------------------------------------------------
def _unpack(part):
    """Accept 3-tuple (kind, params, opacity) or 4-tuple with an extra
    `color` ('black' | 'white'). White primitives act as holes once the
    icon is rendered to PNG (they replace darkness with luminance, and the
    downstream `_tinted_icon` treats luminance as transparency)."""
    if len(part) == 4:
        kind, params, opacity, color = part
    else:
        kind, params, opacity = part
        color = "black"
    return kind, params, opacity, color


def svg_element(kind: str, params: dict, opacity: float,
                color: str = "black") -> str:
    """Serialize one primitive as an SVG element string."""
    # 'currentColor' keeps SVGs recolourable via CSS; 'white' stays literal.
    fill = "currentColor" if color == "black" else color
    attrs = f'fill="{fill}"'
    if opacity != 1.0:
        attrs += f' fill-opacity="{opacity}"'
    if kind == "path":
        return f'<path {attrs} d="{params["d"]}"/>'
    if kind == "ellipse":
        return (f'<ellipse {attrs} cx="{params["cx"]}" cy="{params["cy"]}"'
                f' rx="{params["rx"]}" ry="{params["ry"]}"/>')
    if kind == "circle":
        return f'<circle {attrs} cx="{params["cx"]}" cy="{params["cy"]}" r="{params["r"]}"/>'
    if kind == "rect":
        rx = params.get("rx", 0)
        return (f'<rect {attrs} x="{params["x"]}" y="{params["y"]}"'
                f' width="{params["w"]}" height="{params["h"]}"'
                f' rx="{rx}" ry="{rx}"/>')
    raise ValueError(kind)


def mpl_patch(kind: str, params: dict, opacity: float, color: str = "black"):
    """Create a matplotlib Patch for one primitive."""
    common = dict(facecolor=color, alpha=opacity, edgecolor="none")
    if kind == "path":
        verts, codes = parse_d(params["d"])
        return PathPatch(MPath(verts, codes), **common)
    if kind == "ellipse":
        return Ellipse((params["cx"], params["cy"]),
                       width=params["rx"] * 2, height=params["ry"] * 2, **common)
    if kind == "circle":
        return Circle((params["cx"], params["cy"]), radius=params["r"], **common)
    if kind == "rect":
        # matplotlib Rectangle doesn't do rounded corners; approximate by
        # using FancyBboxPatch if rx>0. For our purposes, sharp corners are fine.
        return Rectangle((params["x"], params["y"]), params["w"], params["h"], **common)
    raise ValueError(kind)


def render(name: str, parts):
    # SVG
    body = "".join(svg_element(*_unpack(p)) for p in parts
                   if _unpack(p)[2] > 0)
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="200" height="200" viewBox="0 0 {VB} {VB}">'
        f"{body}</svg>"
    )
    (LIB / f"{name}.svg").write_text(svg, encoding="utf-8")

    # PNG
    fig, ax = plt.subplots(figsize=(2, 2))
    for part in parts:
        k, p, op, color = _unpack(part)
        if op <= 0:
            continue
        ax.add_patch(mpl_patch(k, p, op, color))
    ax.set_xlim(0, VB)
    ax.set_ylim(0, VB)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(
        LIB / f"{name}.png",
        dpi=150, transparent=True,
        bbox_inches="tight", pad_inches=0,
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
