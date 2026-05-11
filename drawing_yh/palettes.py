"""
drawing_yh.palettes — 命名色板 + 色盲 / 灰度检查

**主张**:同一类别(组织 / cell type / 年龄 / 物种…)在所有子图中**用同一个色板**。
项目色板集中在本文件,出图时 import,不要每个脚本另起一份。

含两类:
- 通用 categorical 色板(`OKABE_ITO`、`DEEP_20`)
- 项目语义色板(`AGE_GRADIENT`、`SEX`、`SPECIES`、`THREE_GROUP`)

外加色盲 / 灰度可区分性检查工具。
"""
from __future__ import annotations

from typing import Sequence, Tuple


# ============================================================
# 通用 categorical 色板
# ============================================================

# Okabe-Ito —— 色盲友好的 8 色 categorical,科研发表黄金标准。
# Ref: https://jfly.uni-koeln.de/color/
OKABE_ITO = [
    "#E69F00",   # orange
    "#56B4E9",   # sky blue
    "#009E73",   # bluish green
    "#F0E442",   # yellow
    "#0072B2",   # blue
    "#D55E00",   # vermillion
    "#CC79A7",   # reddish purple
    "#000000",   # black
]

# Project deep colors —— 20 色饱和,来自 scatter/dot_chart 历史项目,
# 灰度模式可区分性已经验证过。
DEEP_20 = [
    "#c8001e", "#2a8a3a", "#2a4db5", "#c45e10", "#6a0f8a",
    "#1aa0c4", "#b020b0", "#8a9e00", "#c47090", "#2a7070",
    "#8060c0", "#7a4a10", "#600000", "#30a060", "#606000",
    "#c09060", "#000060", "#606060", "#d04020", "#000000",
]


# ============================================================
# 项目语义色板
# ============================================================

# 年龄(老 → 嫩,viridis 风格 5 档,数值/时间序适用)
AGE_GRADIENT = [
    "#440154",   # 老(深紫)
    "#3B528B",
    "#21908C",
    "#5DC863",
    "#FDE725",   # 嫩(亮黄)
]

# 双性别 F / M(色盲友好的橙紫对比)
SEX = {
    "F": "#E66101",   # 橙
    "M": "#5E3C99",   # 紫
}

# 物种(本项目四物种)
SPECIES = {
    "human":  "#1f77b4",
    "monkey": "#ff7f0e",
    "mouse":  "#2ca02c",
    "rat":    "#d62728",
}

# 三档分组(young / mid / old 之类,蓝→灰→红双向偏差)
THREE_GROUP = ["#2a4db5", "#888888", "#c8001e"]


# ============================================================
# 工具函数
# ============================================================

def _hex_to_rgb(h: str) -> Tuple[float, float, float]:
    """`#RRGGBB` → (r, g, b),分量在 0–1。"""
    h = h.lstrip('#')
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def to_grayscale(color: str) -> float:
    """
    把 hex 颜色映射为灰度亮度(0–1)。用感知加权(Rec. 709)。
    用于灰度可区分性检查 —— 期刊打印有时只用黑白,色板里两个颜色
    灰度值差 < 0.1 通常就太接近、看不出区别。
    """
    r, g, b = _hex_to_rgb(color)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def grayscale_distinct(palette: Sequence[str], threshold: float = 0.1) -> bool:
    """
    色板里任意两色在灰度模式下亮度差是否都 >= `threshold`?
    True 表示色盲 / 灰度打印下也可区分。
    """
    grays = sorted(to_grayscale(c) for c in palette)
    return all(grays[i + 1] - grays[i] >= threshold for i in range(len(grays) - 1))


def preview_palette(palette: Sequence[str], save_as=None, title: str = None):
    """
    画色板预览:上排原色,下排灰度等效。一眼看出色盲 / 灰度下哪些会撞。

    Parameters
    ----------
    palette : list of hex strings
    save_as : str or Path or None
        给路径就保存(用 `drawing_yh.save_fig`,自动 SVG dpi=72 等)。
    title : optional
    """
    import matplotlib.pyplot as plt
    n = len(palette)
    fig, (ax_c, ax_g) = plt.subplots(
        2, 1, figsize=(max(n * 0.6, 4), 1.6), gridspec_kw={'hspace': 0.3}
    )
    for ax, mode in zip((ax_c, ax_g), ('color', 'gray')):
        for i, c in enumerate(palette):
            face = c if mode == 'color' else (to_grayscale(c),) * 3
            ax.bar(i, 1, color=face, edgecolor='black', linewidth=0.5)
            ax.text(i, -0.25, str(i), ha='center', va='top', fontsize=7)
        ax.set_xlim(-0.6, n - 0.4)
        ax.set_ylim(0, 1.2)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_ylabel(mode, fontsize=8)
        for sp in ax.spines.values():
            sp.set_visible(False)
    if title:
        fig.suptitle(title, fontsize=9)
    if save_as is not None:
        from .io import save_fig
        save_fig(fig, save_as)
    return fig
