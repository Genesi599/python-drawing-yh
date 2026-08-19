# -*- coding: utf-8 -*-
"""
Generate pancreas + adrenal icons (black silhouette, 与库内风格统一).

- pancreas: 横置蝌蚪形 — 胰头(右端膨大球)+胰体(渐细长条)+胰尾(左端尖) + 浅锯齿
- adrenal:  三角帽形(金字塔)+顶部小嵴 — 肾上腺典型 shape
输出: lib/<name>.svg + .png (fill=currentColor)
"""
from pathlib import Path as PathLib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, Polygon

LIB = PathLib(__file__).parent / "lib"
VB = 256


def build_pancreas():
    """横置胰腺: 右端胰头(大圆) → 胰体(锥形长条) → 左端胰尾(尖)."""
    prims = []

    def add(el, svg):
        prims.append((el, svg))

    # 胰头: 右端大圆(带十二指肠压迹的近似)
    add(Circle((196, 128), 34), '<circle cx="196" cy="128" r="34"/>')
    # 胰头内侧小圆过渡
    add(Circle((170, 128), 24), '<circle cx="170" cy="128" r="24"/>')
    # 胰体: 锥形长条(右厚左薄), 上缘带浅锯齿
    body_pts = [
        (150, 104), (140, 112), (128, 106), (116, 114), (104, 108),
        (92, 116), (80, 110), (66, 116), (52, 122), (40, 126),  # 上缘锯齿
        (30, 130), (36, 138), (48, 140),                         # 尾尖折返
        (62, 136), (76, 138), (88, 136), (100, 140), (112, 138),
        (124, 142), (136, 140), (150, 152),                      # 下缘
    ]
    add(Polygon(body_pts, closed=True),
        '<polygon points="' + " ".join(f"{x},{y}" for x, y in body_pts) + '"/>')
    # 胰尾小圆补充(使尾端圆润)
    add(Circle((38, 132), 12), '<circle cx="38" cy="132" r="12"/>')
    return prims


def build_adrenal():
    """肾上腺: 三角帽形主体 + 顶部嵴."""
    prims = []

    def add(el, svg):
        prims.append((el, svg))

    # 主体: 倒三角帽(底边在下, 顶部尖)
    tri = [(128, 60), (176, 148), (80, 148)]
    add(Polygon(tri, closed=True),
        '<polygon points="128,60 176,148 80,148"/>')
    # 底部圆角(弧形底边): 三个圆叠出弯曲底缘
    add(Circle((100, 146), 18), '<circle cx="100" cy="146" r="18"/>')
    add(Circle((128, 152), 20), '<circle cx="128" cy="152" r="20"/>')
    add(Circle((156, 146), 18), '<circle cx="156" cy="146" r="18"/>')
    # 顶部嵴(中央小突起)
    ridge = [(120, 66), (128, 44), (136, 66)]
    add(Polygon(ridge, closed=True),
        '<polygon points="120,66 128,44 136,66"/>')
    # 内侧压迹(左右两个小凹, 用白点模拟—略, 保持纯剪影)
    return prims


def render(name, prims, color="#111111"):
    svg_parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {VB} {VB}">',
                 f'<g fill="{color}">']
    for el, svg in prims:
        svg_parts.append(svg)
    svg_parts.append('</g></svg>')
    (LIB / f"{name}.svg").write_text("\n".join(svg_parts), encoding="utf-8")
    print(f"[SVG] {name}")

    fig, ax = plt.subplots(figsize=(2.2, 2.2), dpi=232)
    ax.set_xlim(0, VB)
    ax.set_ylim(VB, 0)
    ax.set_aspect("equal")
    ax.axis("off")
    for el, _svg in prims:
        el.set_facecolor(color)
        el.set_edgecolor("none")
        ax.add_patch(el)
    fig.savefig(LIB / f"{name}.png", transparent=True,
                bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"[PNG] {name}")


if __name__ == "__main__":
    LIB.mkdir(parents=True, exist_ok=True)
    render("pancreas", build_pancreas())
    render("adrenal", build_adrenal())
