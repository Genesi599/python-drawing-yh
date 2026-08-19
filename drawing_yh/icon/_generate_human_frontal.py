# -*- coding: utf-8 -*-
"""
Generate a frontal (anterior-view) standing human silhouette icon.

输出: lib/human_frontal.svg + lib/human_frontal.png
风格与 generate_missing_icons.py 一致: 基本图元组合, fill=currentColor。
正面站立, 双臂略外张, 双腿分开 — 供 body map 做背景剪影。
肢体用四角点 polygon(法向偏移), 避免 FancyBboxPatch transform 问题。
"""
from pathlib import Path as PathLib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, Polygon

LIB = PathLib(__file__).parent / "lib"
VB = 256


def limb_poly(p0, p1, w):
    """两端点定义的四肢多边形(矩形, 宽 w), 附加两端半圆由圆补."""
    x0, y0 = p0
    x1, y1 = p1
    dx, dy = x1 - x0, y1 - y0
    L = np.hypot(dx, dy)
    nx, ny = -dy / L * w / 2, dx / L * w / 2  # 法向
    pts = [(x0 + nx, y0 + ny), (x1 + nx, y1 + ny),
           (x1 - nx, y1 - ny), (x0 - nx, y0 - ny)]
    return Polygon(pts, closed=True), pts


def build_prims():
    """返回 (patch, svg_str or None) 列表; None → 由 limb 顶点生成 polygon svg."""
    prims = []

    def add(el, svg):
        prims.append((el, svg))

    # 头 + 颈
    add(Circle((128, 34), 20), '<circle cx="128" cy="34" r="20"/>')
    add(Polygon([(120, 52), (136, 52), (138, 66), (118, 66)], closed=True),
        '<polygon points="120,52 136,52 138,66 118,66"/>')
    # 躯干(肩宽腰窄), 与颈重叠
    add(Polygon([(70, 66), (186, 66), (172, 132), (84, 132)], closed=True),
        '<polygon points="70,66 186,66 172,132 84,132"/>')
    # 肩圆角
    add(Ellipse((74, 72), 24, 20), '<ellipse cx="74" cy="72" rx="12" ry="10"/>')
    add(Ellipse((182, 72), 24, 20), '<ellipse cx="182" cy="72" rx="12" ry="10"/>')
    # 左臂: 肩(72,74)→肘(56,116)→腕(48,156), 与肩圆重叠保证连接
    el, pts = limb_poly((70, 74), (54, 116), 17)
    svg = '<polygon points="' + " ".join(f"{x:.0f},{y:.0f}" for x, y in pts) + '"/>'
    add(el, svg)
    el, pts = limb_poly((54, 116), (47, 156), 14)
    svg = '<polygon points="' + " ".join(f"{x:.0f},{y:.0f}" for x, y in pts) + '"/>'
    add(el, svg)
    # 右臂(镜像 x → 256-x)
    el, pts = limb_poly((186, 74), (202, 116), 17)
    svg = '<polygon points="' + " ".join(f"{x:.0f},{y:.0f}" for x, y in pts) + '"/>'
    add(el, svg)
    el, pts = limb_poly((202, 116), (209, 156), 14)
    svg = '<polygon points="' + " ".join(f"{x:.0f},{y:.0f}" for x, y in pts) + '"/>'
    add(el, svg)
    # 手
    add(Circle((46, 166), 10), '<circle cx="46" cy="166" r="10"/>')
    add(Circle((210, 166), 10), '<circle cx="210" cy="166" r="10"/>')
    # 髋(与躯干重叠)
    add(Polygon([(84, 128), (172, 128), (166, 160), (90, 160)], closed=True),
        '<polygon points="84,128 172,128 166,160 90,160"/>')
    # 双腿(与髋重叠): 左(98,154)→(95,232), 右(158,154)→(161,232)
    el, pts = limb_poly((99, 150), (95, 232), 24)
    svg = '<polygon points="' + " ".join(f"{x:.0f},{y:.0f}" for x, y in pts) + '"/>'
    add(el, svg)
    el, pts = limb_poly((157, 150), (161, 232), 24)
    svg = '<polygon points="' + " ".join(f"{x:.0f},{y:.0f}" for x, y in pts) + '"/>'
    add(el, svg)
    # 足
    add(Ellipse((92, 242), 28, 12), '<ellipse cx="92" cy="242" rx="14" ry="6"/>')
    add(Ellipse((164, 242), 28, 12), '<ellipse cx="164" cy="242" rx="14" ry="6"/>')
    return prims


def render_svg(path):
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {VB} {VB}">',
             '<g fill="currentColor">']
    for el, svg in build_prims():
        parts.append(svg)
    parts.append('</g></svg>')
    path.write_text("\n".join(parts), encoding="utf-8")
    print(f"[SVG] {path}")


def render_png(path, color="#4B5563"):
    fig, ax = plt.subplots(figsize=(4, 4), dpi=128)
    ax.set_xlim(0, VB)
    ax.set_ylim(VB, 0)
    ax.set_aspect("equal")
    ax.axis("off")
    for el, _svg in build_prims():
        el.set_facecolor(color)
        el.set_edgecolor("none")
        ax.add_patch(el)
    fig.savefig(path, transparent=True, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"[PNG] {path}")


if __name__ == "__main__":
    LIB.mkdir(parents=True, exist_ok=True)
    render_svg(LIB / "human_frontal.svg")
    render_png(LIB / "human_frontal.png")
