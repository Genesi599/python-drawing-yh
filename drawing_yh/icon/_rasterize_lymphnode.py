"""
Rasterize lymphnode.svg -> lymphnode.png using svgpathtools + matplotlib.
(Avoids cairo dependency.)
"""
import os
import re
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch

from svgpathtools import parse_path, Line, CubicBezier, QuadraticBezier, Arc

HERE = os.path.dirname(os.path.abspath(__file__))
SVG_PATH = os.path.join(HERE, "lib", "lymphnode.svg")
PNG_PATH = os.path.join(HERE, "lib", "lymphnode.png")

SIZE_PX = 512  # output icon size


def segment_to_mpl(seg, first=False):
    """Return (points, codes) for one segment. If first, prepend MOVETO at start."""
    pts, codes = [], []
    if first:
        pts.append((seg.start.real, seg.start.imag))
        codes.append(MplPath.MOVETO)
    if isinstance(seg, Line):
        pts.append((seg.end.real, seg.end.imag))
        codes.append(MplPath.LINETO)
    elif isinstance(seg, CubicBezier):
        pts += [(seg.control1.real, seg.control1.imag),
                (seg.control2.real, seg.control2.imag),
                (seg.end.real, seg.end.imag)]
        codes += [MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    elif isinstance(seg, QuadraticBezier):
        pts += [(seg.control.real, seg.control.imag),
                (seg.end.real, seg.end.imag)]
        codes += [MplPath.CURVE3, MplPath.CURVE3]
    elif isinstance(seg, Arc):
        # Approximate Arc by sampling.
        samples = 30
        for t in np.linspace(0, 1, samples)[1:]:
            p = seg.point(t)
            pts.append((p.real, p.imag))
            codes.append(MplPath.LINETO)
    else:
        pts.append((seg.end.real, seg.end.imag))
        codes.append(MplPath.LINETO)
    return pts, codes


def svg_path_to_mpl(d_str):
    """Parse SVG path string into a single compound matplotlib Path."""
    p = parse_path(d_str)
    all_pts, all_codes = [], []
    for sub in p.continuous_subpaths():
        for i, seg in enumerate(sub):
            pts, codes = segment_to_mpl(seg, first=(i == 0))
            all_pts += pts
            all_codes += codes
        # close
        all_pts.append(all_pts[-1])
        all_codes.append(MplPath.CLOSEPOLY)
    return MplPath(all_pts, all_codes)


def main():
    tree = ET.parse(SVG_PATH)
    root = tree.getroot()
    # namespace strip
    ns = re.match(r"\{.*\}", root.tag)
    ns = ns.group(0) if ns else ""
    viewBox = root.get("viewBox", "0 0 1024 1024").split()
    vb_x, vb_y, vb_w, vb_h = map(float, viewBox)

    fig, ax = plt.subplots(figsize=(SIZE_PX / 100, SIZE_PX / 100), dpi=100)
    ax.set_xlim(vb_x, vb_x + vb_w)
    ax.set_ylim(vb_y + vb_h, vb_y)  # SVG y grows downward
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1)

    for el in root.findall(f"{ns}path"):
        d = el.get("d")
        fill = el.get("fill", "#000000")
        opacity = float(el.get("fill-opacity", 1.0))
        if not d or fill.lower() == "none":
            continue
        mpl_path = svg_path_to_mpl(d)
        patch = PathPatch(
            mpl_path,
            facecolor=fill,
            edgecolor="none",
            alpha=opacity,
            linewidth=0,
        )
        ax.add_patch(patch)

    fig.savefig(PNG_PATH, dpi=100, transparent=True)
    plt.close(fig)
    print(f"saved: {PNG_PATH}")


if __name__ == "__main__":
    main()
