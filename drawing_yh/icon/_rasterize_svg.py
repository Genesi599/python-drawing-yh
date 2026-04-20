"""
Rasterize <lib>/<name>.svg -> <lib>/<name>.png using svgpathtools + matplotlib.
Avoids cairo dependency.

Usage:
    python _rasterize_svg.py lymphnode duodenum
    python _rasterize_svg.py --all
"""
import argparse
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch

from svgpathtools import parse_path, Line, CubicBezier, QuadraticBezier, Arc

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "lib")
SIZE_PX = 512


def seg_to_mpl(seg, first=False):
    pts, codes = [], []
    if first:
        pts.append((seg.start.real, seg.start.imag))
        codes.append(MplPath.MOVETO)
    if isinstance(seg, Line):
        pts.append((seg.end.real, seg.end.imag)); codes.append(MplPath.LINETO)
    elif isinstance(seg, CubicBezier):
        pts += [(seg.control1.real, seg.control1.imag),
                (seg.control2.real, seg.control2.imag),
                (seg.end.real, seg.end.imag)]
        codes += [MplPath.CURVE4]*3
    elif isinstance(seg, QuadraticBezier):
        pts += [(seg.control.real, seg.control.imag),
                (seg.end.real, seg.end.imag)]
        codes += [MplPath.CURVE3]*2
    elif isinstance(seg, Arc):
        for t in np.linspace(0, 1, 30)[1:]:
            p = seg.point(t)
            pts.append((p.real, p.imag)); codes.append(MplPath.LINETO)
    else:
        pts.append((seg.end.real, seg.end.imag)); codes.append(MplPath.LINETO)
    return pts, codes


def svg_d_to_mpl_path(d_str):
    p = parse_path(d_str)
    all_pts, all_codes = [], []
    for sub in p.continuous_subpaths():
        for i, seg in enumerate(sub):
            pts, codes = seg_to_mpl(seg, first=(i == 0))
            all_pts += pts; all_codes += codes
        all_pts.append(all_pts[-1]); all_codes.append(MplPath.CLOSEPOLY)
    return MplPath(all_pts, all_codes)


def _parse_transform(s):
    """Return affine matrix [[a,c,e],[b,d,f]] for a limited set of SVG transforms."""
    if not s:
        return None
    import numpy as np
    M = np.eye(3)
    for m in re.finditer(r"(matrix|translate|scale|rotate)\s*\(([^)]*)\)", s):
        kind = m.group(1)
        vals = [float(v) for v in re.split(r"[\s,]+", m.group(2).strip())]
        if kind == "matrix":
            a, b, c, d, e, f = vals
            T = np.array([[a, c, e], [b, d, f], [0, 0, 1]])
        elif kind == "translate":
            tx = vals[0]; ty = vals[1] if len(vals) > 1 else 0
            T = np.array([[1, 0, tx], [0, 1, ty], [0, 0, 1]])
        elif kind == "scale":
            sx = vals[0]; sy = vals[1] if len(vals) > 1 else sx
            T = np.array([[sx, 0, 0], [0, sy, 0], [0, 0, 1]])
        elif kind == "rotate":
            th = np.radians(vals[0])
            T = np.array([[np.cos(th), -np.sin(th), 0],
                          [np.sin(th),  np.cos(th), 0],
                          [0, 0, 1]])
        M = M @ T
    return M


def _apply_transform(pts, M):
    import numpy as np
    if M is None:
        return pts
    arr = np.array(pts)
    homo = np.hstack([arr, np.ones((len(arr), 1))])
    out = (M @ homo.T).T[:, :2]
    return [tuple(p) for p in out]


def rasterize(name, size_px=SIZE_PX, verbose=True):
    svg_path = os.path.join(LIB, f"{name}.svg")
    png_path = os.path.join(LIB, f"{name}.png")
    if not os.path.exists(svg_path):
        print(f"[skip] no svg: {svg_path}")
        return False

    tree = ET.parse(svg_path)
    root = tree.getroot()
    ns_m = re.match(r"\{.*\}", root.tag)
    ns = ns_m.group(0) if ns_m else ""

    vb = root.get("viewBox")
    if vb:
        vb_x, vb_y, vb_w, vb_h = map(float, vb.split())
    else:
        vb_x, vb_y = 0.0, 0.0
        vb_w = float(root.get("width", 1024))
        vb_h = float(root.get("height", 1024))

    fig, ax = plt.subplots(figsize=(size_px / 100, size_px / 100), dpi=100)
    ax.set_xlim(vb_x, vb_x + vb_w)
    ax.set_ylim(vb_y + vb_h, vb_y)
    ax.set_aspect("equal")
    ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1)

    def walk(el, parent_M=None, parent_fill=None):
        M = _parse_transform(el.get("transform"))
        if parent_M is not None and M is not None:
            M = parent_M @ M
        elif parent_M is not None:
            M = parent_M
        fill = el.get("fill", parent_fill if parent_fill else "#000000")
        op = float(el.get("fill-opacity", 1.0))
        tag = el.tag.replace(ns, "")

        if tag == "path" and el.get("d") and fill.lower() != "none":
            mpl_path = svg_d_to_mpl_path(el.get("d"))
            if M is not None:
                verts = _apply_transform(mpl_path.vertices, M)
                mpl_path = MplPath(verts, mpl_path.codes)
            ax.add_patch(PathPatch(mpl_path, facecolor=fill, edgecolor="none",
                                   alpha=op, linewidth=0))
        elif tag == "rect":
            x = float(el.get("x", 0)); y = float(el.get("y", 0))
            w = float(el.get("width", 0)); h = float(el.get("height", 0))
            rx = float(el.get("rx", 0));  ry = float(el.get("ry", rx))
            from matplotlib.patches import FancyBboxPatch, Rectangle
            if rx or ry:
                # approximate by drawing rounded rect via Path
                pts = [(x+rx, y), (x+w-rx, y),
                       (x+w-rx, y), (x+w, y), (x+w, y+ry),  # corner
                       (x+w, y+h-ry),
                       (x+w, y+h-ry), (x+w, y+h), (x+w-rx, y+h),
                       (x+rx, y+h),
                       (x+rx, y+h), (x, y+h), (x, y+h-ry),
                       (x, y+ry),
                       (x, y+ry), (x, y), (x+rx, y)]
                codes = [MplPath.MOVETO, MplPath.LINETO,
                         MplPath.CURVE3, MplPath.CURVE3, MplPath.LINETO,
                         MplPath.LINETO,
                         MplPath.CURVE3, MplPath.CURVE3, MplPath.LINETO,
                         MplPath.LINETO,
                         MplPath.CURVE3, MplPath.CURVE3, MplPath.LINETO,
                         MplPath.LINETO,
                         MplPath.CURVE3, MplPath.CURVE3, MplPath.CLOSEPOLY]
                if M is not None:
                    pts = _apply_transform(pts, M)
                ax.add_patch(PathPatch(MplPath(pts, codes), facecolor=fill,
                                       edgecolor="none", alpha=op, linewidth=0))
            else:
                from matplotlib.patches import Rectangle
                ax.add_patch(Rectangle((x, y), w, h, facecolor=fill,
                                       edgecolor="none", alpha=op, linewidth=0))
        elif tag in ("circle", "ellipse"):
            cx = float(el.get("cx", 0)); cy = float(el.get("cy", 0))
            if tag == "circle":
                rx = ry = float(el.get("r", 0))
            else:
                rx = float(el.get("rx", 0)); ry = float(el.get("ry", 0))
            from matplotlib.patches import Ellipse
            ax.add_patch(Ellipse((cx, cy), 2*rx, 2*ry, facecolor=fill,
                                 edgecolor="none", alpha=op, linewidth=0))
        elif tag == "g":
            pass  # children handled below

        for child in el:
            walk(child, M, fill if fill else parent_fill)

    walk(root)

    fig.savefig(png_path, dpi=100, transparent=True)
    plt.close(fig)

    # Auto-crop transparent borders so the icon fills the PNG.
    _crop_transparent(png_path, pad_ratio=0.04)

    if verbose:
        print(f"saved: {png_path}")
    return True


def _crop_transparent(png_path, pad_ratio=0.04):
    """Crop transparent (alpha==0) margins and re-pad with a thin transparent border."""
    from PIL import Image
    im = Image.open(png_path).convert("RGBA")
    arr = np.array(im)
    alpha = arr[..., 3]
    ys, xs = np.where(alpha > 0)
    if ys.size == 0:
        return
    y0, y1 = ys.min(), ys.max() + 1
    x0, x1 = xs.min(), xs.max() + 1
    cropped = im.crop((x0, y0, x1, y1))

    # square canvas with small pad
    w, h = cropped.size
    side = max(w, h)
    pad = int(side * pad_ratio)
    canvas_side = side + 2 * pad
    canvas = Image.new("RGBA", (canvas_side, canvas_side), (0, 0, 0, 0))
    off_x = (canvas_side - w) // 2
    off_y = (canvas_side - h) // 2
    canvas.paste(cropped, (off_x, off_y), cropped)
    canvas.save(png_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()

    if args.all:
        names = [os.path.splitext(f)[0] for f in os.listdir(LIB) if f.endswith(".svg")]
    else:
        names = args.names
    if not names:
        print("Usage: _rasterize_svg.py <name> [<name> ...]  or  --all")
        sys.exit(1)
    for n in names:
        rasterize(n)


if __name__ == "__main__":
    main()
