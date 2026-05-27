"""
hub_spoke — TF hub-spoke 同心圆网络图(衰老课题组标准范式)

参考范本(Liu / Qu lab 出图惯例):
- Y 2023 FOXO3 骨骼肌 Fig4-A: 7 core TFs 内圈 + 14 cell types 外圈
- Huang 2023 WT1 睾丸 Fig6-A: 1 hub TF 中心 + 9 inner TFs + ~150 target genes
- Yang 2023 Liver SREBF2 Fig5-F: hub + core TFs + targets,带 zonation 环

拓扑(3 层同心圆,极坐标布局,不用 networkx spring layout):
- L0 = hub: 中心 (r=0, 单 TF) 或小圈 (r=hub_r, 多 TF 平均铺)
- L1 = mid ring (optional, r=mid_r ≈ 0.5): core/inner TFs
- L2 = outer ring (r=outer_r=1.0): cell type 或 target gene, 角度均匀
- 角度 theta_i = 90 deg - i * 360 deg / N (12 点钟开始顺时针)

视觉编码硬约定:
- 节点大小 ∝ target_count; 3 档 size legend 放图外右下
- 节点颜色: 连续值 -> Blues cmap + colorbar; 离散 -> 用 palette dict
- hub TF: 填蓝 #1f78b4, 白字 italic bold, edgecolor 白
- outer: edgecolor 灰, lw=0.4; label 沿径向放 (x>0 ha='left', x<0 ha='right')
- spoke: hub→outer 灰线 lw=0.3 alpha=0.3, zorder<节点
- TF-TF 全连接 (可选): 蓝线 lw=0.3 alpha=0.3

用法:
>>> from drawing_yh import hub_spoke
>>> outer = {
...     "Fib": {"size": 120, "color": 0.85},
...     "EC":  {"size":  80, "color": 0.62},
...     "Mac": {"size":  45, "color": 0.30},
... }
>>> fig, ax = hub_spoke(hub="FOXO3", outer=outer, cmap="Blues",
...                    color_legend_label="hit score")
"""
from __future__ import annotations

import math
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


HUB_FILL = "#1f78b4"
HUB_EDGE = "#ffffff"
HUB_TEXT = "#ffffff"
MID_FILL = "#1f78b4"
MID_EDGE = "#ffffff"
MID_TEXT = "#ffffff"
OUTER_EDGE = "#7e7e7e"
SPOKE_COLOR = "#9e9e9e"
TF_TF_COLOR = "#5a9bd5"


def _polar_xy(theta_deg: float, r: float) -> tuple[float, float]:
    """角度(deg, 12 点顺时针计) + 半径 -> (x, y) 笛卡尔。"""
    a = math.radians(90.0 - theta_deg)
    return r * math.cos(a), r * math.sin(a)


def _ring_positions(n: int, r: float, start_deg: float = 0.0) -> list[tuple[float, float]]:
    """N 个节点均匀铺在半径 r 的环上,从 12 点顺时针。"""
    if n == 0:
        return []
    return [_polar_xy(start_deg + i * 360.0 / n, r) for i in range(n)]


def _normalize_size(
    values: Iterable[float], smin: float = 30, smax: float = 600
) -> np.ndarray:
    """target_count -> scatter s 参数(area)。简单线性归一到 [smin, smax]。"""
    v = np.asarray(list(values), dtype=float)
    if len(v) == 0 or v.max() == v.min():
        return np.full(len(v), (smin + smax) / 2)
    return smin + (v - v.min()) / (v.max() - v.min()) * (smax - smin)


def _resolve_color(
    val, *, cmap, norm, palette: dict | None
) -> tuple[float, float, float, float]:
    """连续值走 cmap+norm; 离散 key 走 palette dict; None 走灰。"""
    if val is None:
        return (0.75, 0.75, 0.75, 1.0)
    if palette is not None and val in palette:
        c = palette[val]
        if isinstance(c, str):
            from matplotlib.colors import to_rgba
            return to_rgba(c)
        return c
    try:
        v = float(val)
    except (TypeError, ValueError):
        return (0.75, 0.75, 0.75, 1.0)
    if norm is not None and cmap is not None:
        return cmap(norm(v))
    return (0.75, 0.75, 0.75, 1.0)


def hub_spoke(
    hub,                                       # str 或 list[str]
    outer: dict[str, dict],                    # {label: {'size': ..., 'color': ...}}
    mid: dict[str, dict] | None = None,        # {label: {'size': ..., 'color': ...}}
    layout: dict | None = None,                # {'hub_r','mid_r','outer_r'}
    cmap: str | object = "Blues",
    vmin: float | None = None,
    vmax: float | None = None,
    palette: dict | None = None,
    spoke: bool = True,
    spoke_style: dict | None = None,
    tf_tf: bool = False,                       # mid TF 之间是否全连接 (Huang 风格)
    size_legend: tuple = (10, 50, 100),
    color_legend_label: str = "",
    title: str = "",
    label_outer: bool = True,
    label_radius_factor: float = 1.18,
    figsize: tuple | None = None,
    ax: object | None = None,
):
    """画 TF hub-spoke 同心圆 network。

    Parameters
    ----------
    hub : str | list[str]
        中心 hub。单 TF -> 中心点; 多 TF -> 内圈小圆 (r=layout['hub_r'])。
    outer : dict[label, dict]
        外圈节点 (cell type 或 target gene)。每 dict 含:
          - 'size': target_count (int/float), 用于节点大小
          - 'color': float (走 cmap) 或 str key (走 palette)
    mid : dict[label, dict], optional
        中间环 (mid ring) 节点,通常是 core/inner TFs。
    layout : dict, optional
        半径定义。默认 {'hub_r':0.18 (多 hub 时), 'mid_r':0.5, 'outer_r':1.0}。
    cmap : str | Colormap
        颜色映射 (连续值)。默认 'Blues'。
    vmin, vmax : float
        cmap norm 上下界。None 则按 outer color 自动取 min/max。
    palette : dict, optional
        离散 color key -> color hex。优先于 cmap。
    spoke : bool
        是否画 hub -> outer 灰线。
    spoke_style : dict
        spoke 样式覆盖 {'color','lw','alpha'}。
    tf_tf : bool
        mid TFs 之间是否全两两连接 (Huang 2023 WT1 风格)。
    size_legend : tuple
        size legend 3 档参考值(原始 target_count 单位)。
    color_legend_label : str
        colorbar label。
    title : str
        图内 title (短,详细 caption 写到 figcaption)。
    label_outer : bool
        是否标外圈 label (target gene 多时建议 False)。
    label_radius_factor : float
        label 距 outer 节点的径向系数 (默认 r * 1.18)。
    figsize : tuple, optional
        默认 (5.0, 5.0) 单方形。
    ax : matplotlib Axes, optional
        外部 axes (跳过 figsize 创建)。

    Returns
    -------
    (fig, ax) : tuple
    """
    L = {"hub_r": 0.18, "mid_r": 0.5, "outer_r": 1.0}
    if layout:
        L.update(layout)
    spoke_st = {"color": SPOKE_COLOR, "lw": 0.3, "alpha": 0.3}
    if spoke_style:
        spoke_st.update(spoke_style)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize or (5.0, 5.0))
    else:
        fig = ax.figure

    # ---- 颜色 norm: 用 outer 数值统一 ----
    color_vals = []
    for d in outer.values():
        c = d.get("color")
        if c is None:
            continue
        try:
            color_vals.append(float(c))
        except (TypeError, ValueError):
            pass
    if mid is not None:
        for d in mid.values():
            c = d.get("color")
            if c is None:
                continue
            try:
                color_vals.append(float(c))
            except (TypeError, ValueError):
                pass
    if color_vals and palette is None:
        if vmin is None:
            vmin = min(color_vals)
        if vmax is None:
            vmax = max(color_vals)
        if vmin == vmax:
            vmax = vmin + 1e-9
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = None
    if isinstance(cmap, str):
        cmap_obj = plt.get_cmap(cmap)
    else:
        cmap_obj = cmap

    # ---- L2: outer ring ----
    outer_keys = list(outer.keys())
    outer_pos = _ring_positions(len(outer_keys), L["outer_r"])
    outer_sizes = _normalize_size([outer[k].get("size", 1) for k in outer_keys])
    outer_colors = [_resolve_color(outer[k].get("color"), cmap=cmap_obj,
                                    norm=norm, palette=palette) for k in outer_keys]
    outer_xy = dict(zip(outer_keys, outer_pos))

    # ---- L1: mid ring (optional) ----
    mid_keys = list(mid.keys()) if mid else []
    mid_pos = _ring_positions(len(mid_keys), L["mid_r"]) if mid_keys else []
    mid_sizes = (_normalize_size([mid[k].get("size", 1) for k in mid_keys],
                                  smin=200, smax=900)
                 if mid_keys else np.array([]))
    mid_colors = [_resolve_color(mid[k].get("color"), cmap=cmap_obj,
                                  norm=norm, palette=palette) for k in mid_keys]
    mid_xy = dict(zip(mid_keys, mid_pos))

    # ---- L0: hub ----
    if isinstance(hub, str):
        hub_list = [hub]
        hub_pos = [(0.0, 0.0)]
    else:
        hub_list = list(hub)
        hub_pos = _ring_positions(len(hub_list), L["hub_r"])
    hub_xy = dict(zip(hub_list, hub_pos))

    # ---- 画 spokes (hub -> outer) ----
    if spoke:
        for h in hub_list:
            hx, hy = hub_xy[h]
            for k in outer_keys:
                ox, oy = outer_xy[k]
                ax.plot([hx, ox], [hy, oy],
                        color=spoke_st["color"], lw=spoke_st["lw"],
                        alpha=spoke_st["alpha"], zorder=1)
        # 如果有 mid, hub→mid 也画细 spoke
        if mid_keys:
            for h in hub_list:
                hx, hy = hub_xy[h]
                for k in mid_keys:
                    mx, my = mid_xy[k]
                    ax.plot([hx, mx], [hy, my],
                            color=spoke_st["color"], lw=spoke_st["lw"],
                            alpha=spoke_st["alpha"], zorder=1)

    # ---- TF-TF 全连接 (mid 内, Huang 风格) ----
    if tf_tf and len(mid_keys) >= 2:
        for i in range(len(mid_keys)):
            for j in range(i + 1, len(mid_keys)):
                x1, y1 = mid_pos[i]; x2, y2 = mid_pos[j]
                ax.plot([x1, x2], [y1, y2], color=TF_TF_COLOR,
                        lw=0.3, alpha=0.3, zorder=1)
        # hub 多 TF 时 hub 之间也连
        if len(hub_list) >= 2:
            for i in range(len(hub_list)):
                for j in range(i + 1, len(hub_list)):
                    x1, y1 = hub_pos[i]; x2, y2 = hub_pos[j]
                    ax.plot([x1, x2], [y1, y2], color=TF_TF_COLOR,
                            lw=0.3, alpha=0.3, zorder=1)

    # ---- 画节点(scatter,zorder 高) ----
    # outer
    ox_arr = np.array([outer_xy[k][0] for k in outer_keys])
    oy_arr = np.array([outer_xy[k][1] for k in outer_keys])
    ax.scatter(ox_arr, oy_arr, s=outer_sizes, c=outer_colors,
               edgecolor=OUTER_EDGE, lw=0.4, zorder=3)

    # mid
    if mid_keys:
        mx_arr = np.array([mid_xy[k][0] for k in mid_keys])
        my_arr = np.array([mid_xy[k][1] for k in mid_keys])
        ax.scatter(mx_arr, my_arr, s=mid_sizes, c=mid_colors,
                   edgecolor=MID_EDGE, lw=0.6, zorder=4)

    # hub
    hx_arr = np.array([hub_xy[h][0] for h in hub_list])
    hy_arr = np.array([hub_xy[h][1] for h in hub_list])
    hub_size = 1100 if len(hub_list) == 1 else 700
    ax.scatter(hx_arr, hy_arr, s=hub_size, c=[HUB_FILL] * len(hub_list),
               edgecolor=HUB_EDGE, lw=0.8, zorder=5)

    # ---- 节点 label ----
    # hub TF symbol (中心, italic bold 白字)
    for h, (x, y) in zip(hub_list, hub_pos):
        ax.text(x, y, h, ha="center", va="center", fontsize=9,
                fontweight="bold", fontstyle="italic", color=HUB_TEXT,
                zorder=6)
    # mid TF symbol
    for m, (x, y) in zip(mid_keys, mid_pos):
        ax.text(x, y, m, ha="center", va="center", fontsize=8,
                fontweight="bold", fontstyle="italic", color=MID_TEXT,
                zorder=6)
    # outer label (径向放, x>0 ha='left' / x<0 ha='right')
    if label_outer:
        for k, (x, y) in zip(outer_keys, outer_pos):
            r_lbl = L["outer_r"] * label_radius_factor
            scale = r_lbl / max(L["outer_r"], 1e-9)
            lx, ly = x * scale, y * scale
            ha = "left" if x > 0.02 else ("right" if x < -0.02 else "center")
            va = "bottom" if y > 0.02 else ("top" if y < -0.02 else "center")
            ax.text(lx, ly, k, ha=ha, va=va, fontsize=8, color="#222",
                    zorder=6)

    # ---- title ----
    if title:
        ax.set_title(title, fontsize=8)

    # ---- 坐标轴 ----
    pad = L["outer_r"] * (label_radius_factor + 0.35)
    ax.set_xlim(-pad, pad); ax.set_ylim(-pad, pad)
    ax.set_aspect("equal")
    ax.set_axis_off()

    # ---- size legend 3 档 (右下图外) ----
    if size_legend and len(size_legend) >= 1:
        ref_vals = list(size_legend)
        # 同一个 _normalize_size 映射 (用 outer raw counts)
        all_counts = [outer[k].get("size", 1) for k in outer_keys]
        if all_counts and max(all_counts) > min(all_counts):
            cmin, cmax = min(all_counts), max(all_counts)
            ref_sizes = [30 + (v - cmin) / (cmax - cmin) * (600 - 30) for v in ref_vals]
        else:
            ref_sizes = [(30 + 600) / 2] * len(ref_vals)
        size_handles = [
            plt.scatter([], [], s=s, c="#bbbbbb", edgecolor=OUTER_EDGE, lw=0.4,
                        label=f"{int(v)}")
            for v, s in zip(ref_vals, ref_sizes)
        ]
        leg = ax.legend(handles=size_handles, title="Target genes",
                        loc="lower right", bbox_to_anchor=(1.18, -0.05),
                        frameon=False, fontsize=7, title_fontsize=7,
                        labelspacing=1.2, handletextpad=0.6, borderpad=0.4)
        leg.get_title().set_fontsize(7)
        ax.add_artist(leg)

    # ---- colorbar (连续值) ----
    if norm is not None:
        sm = ScalarMappable(norm=norm, cmap=cmap_obj); sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.04,
                            shrink=0.5, aspect=12,
                            location="right", anchor=(0.0, 0.85))
        cbar.set_label(color_legend_label, fontsize=7)
        cbar.ax.tick_params(labelsize=7)

    return fig, ax
