#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
clustermap_annot — 带格内值标注 + 行/列色块的聚类热图(通用版)。

合规 ``drawing_yh.STANDARDS`` 要点:
  - 字号默认走 ``matplotlib.rcParams``(包级 8pt),不硬编码
  - ``figsize=None`` 时按 nrows / ncols 动态计算紧凑尺寸(项目"版面以文字为准")
  - 输出统一通过 ``drawing_yh.save_fig(fig, ..., also=('.pdf', '.svg'))``
  - 颜色 / cmap 复用上游传入,函数内部不写死 hex
  - row_colors / col_colors 支持 ``Series`` 或 ``DataFrame``(多列色块)

典型场景
--------
1. **色 = Z-score(列方向)+ 格内 = 原始量**(本项目最常见):

    >>> from drawing_yh import clustermap_annot, save_fig
    >>> Z = (F - F.mean()) / F.std()                          # 列方向 Z-score
    >>> fig = clustermap_annot(
    ...     Z, F,                                             # 色用 Z,文本用原始 F
    ...     row_colors=group_colors,                          # 功能群色块
    ...     cmap="Reds", vmin=Z.values.min(), vmax=Z.values.max(),
    ...     annot_fmt=lambda v: f"{int(v)}" if float(v).is_integer() else f"{v:.2f}",
    ...     title="...", cbar_label="Z-score (per column)",
    ... )
    >>> save_fig(fig, "out.png", also=(".pdf", ".svg"))

2. **单矩阵 + diverging(相关 r 类)**:

    >>> fig = clustermap_annot(r_mat, cmap="RdBu_r", vmin=-1, vmax=1, center=0,
    ...                        annot_fmt=".2f")

3. **多列 row_colors**(同时标功能群 + celltypist):

    >>> row_cols = pd.DataFrame({
    ...     "group": [COLORS_GRP[g] for g in df.gid],
    ...     "type":  [COLORS_SUB[s] for s in df.subtype],
    ... }, index=Z.index)
    >>> fig = clustermap_annot(Z, F, row_colors=row_cols, ...)

@File    : clustermap_annot.py
@Date    : 2026-06-04
@Author  : yh109
"""
from __future__ import annotations

from itertools import groupby
from typing import Callable, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


def _auto_figsize(n_rows: int, n_cols: int,
                  row_inch: float = 0.18, col_inch: float = 0.45,
                  base_w: float = 1.8, base_h: float = 1.5,
                  min_size: Tuple[float, float] = (3.0, 2.5),
                  max_size: Tuple[float, float] = (10.0, 12.0)) -> Tuple[float, float]:
    """按行列数估算紧凑 figsize(每行 ~0.18in、每列 ~0.45in + 边距)。"""
    w = base_w + n_cols * col_inch
    h = base_h + n_rows * row_inch
    w = min(max(w, min_size[0]), max_size[0])
    h = min(max(h, min_size[1]), max_size[1])
    return (w, h)


def clustermap_annot(
    data_z: pd.DataFrame,
    data_raw: Optional[pd.DataFrame] = None,
    *,
    row_colors: Optional[Union[Sequence, pd.Series, pd.DataFrame]] = None,
    col_colors: Optional[Union[Sequence, pd.Series, pd.DataFrame]] = None,
    row_groups: Optional[Union[Sequence, pd.Series]] = None,
    row_group_labels: Optional[Mapping] = None,
    row_cluster: bool = True,
    col_cluster: bool = False,
    linkage_method: str = "ward",
    metric: str = "euclidean",
    cmap: str = "RdBu_r",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    center: Optional[float] = None,
    col_cmap_groups: Optional[Sequence[Tuple[Sequence[int], str]]] = None,
    col_group_gap: float = 0.0,
    cbar_layout: str = "top-row",
    annot_fmt: Union[str, Callable[[float], str]] = ".2f",
    annot_fontsize: Optional[float] = None,
    annot_dark_threshold: float = 0.55,
    figsize: Optional[Tuple[float, float]] = None,
    cbar_label: str = "",
    title: str = "",
    xtick_rotation: float = 0,
    ytick_right: bool = True,
    xtick_fontsize: Optional[float] = None,
    ytick_fontsize: Optional[float] = None,
    row_colors_label: str = "",
    cbar_pos: Tuple[float, float, float, float] = (0.02, 0.84, 0.04, 0.13),
    dendrogram_ratio: Tuple[float, float] = (0.18, 0.10),
    return_grid: bool = False,
) -> Union[plt.Figure, Tuple[plt.Figure, "sns.matrix.ClusterGrid"]]:
    """带格内值标注 + 行/列色块的聚类热图。

    Parameters
    ----------
    data_z : pd.DataFrame
        驱动**色阶**的矩阵(行 = items,列 = features)。
    data_raw : pd.DataFrame, optional
        驱动**格内文本**的同 shape 矩阵;None 时与 ``data_z`` 同源。
    row_colors, col_colors : array-like / Series / DataFrame, optional
        ``Series`` → 单条色块;``DataFrame`` → 多列色块(seaborn 原生)。
    row_groups : Series-like, optional
        每行的 group id(用于在 row_colors 旁标群名)。dendrogram 重排后,
        函数会找连续同 group 的段并在每段中心写 ``row_group_labels[gid]``。
    row_group_labels : dict {gid: label}, optional
        与 ``row_groups`` 配套。给定时在 row_colors 色块左侧 / dendrogram 区
        标群名(沿 y 轴对齐到段中心,水平方向 ha='right')。
    row_cluster, col_cluster : bool
        默认行聚类、列不聚类(列 = 固定 features)。
    linkage_method, metric : str
        ``scipy.cluster.hierarchy.linkage`` 参数。
    cmap, vmin, vmax, center : matplotlib heatmap 标准参数(单 cmap 模式)。
    col_cmap_groups : list of (col_indices, cmap_name), optional
        **每组列独立 cmap**(覆盖 ``cmap``)。例:
        ``[([0], "Purples"), ([1, 2], "Reds")]`` 第 1 列 Purples、后两列 Reds。
        每组独立计算 vmin/vmax(该组所有 cells 的实际范围)。
    col_group_gap : float, default 0.0
        ``col_cmap_groups`` 模式下列组之间的视觉间隙(以列宽为单位)。0 = 列贴紧,
        ``0.3`` 表示组间留 0.3 个列宽空白(便于一眼区分维度)。
    cbar_layout : {'top-row', 'left-stack', 'bottom-row'}, default 'top-row'
        ``col_cmap_groups`` 模式下多色条的排布:
        - ``top-row``: 顶部横排(默认,不挡 row_colors / yticks)
        - ``left-stack``: 左上堆叠(老行为)
        - ``bottom-row``: 底部横排
    annot_fmt : str or callable
        格内文本 format。字符串 → ``format(v, fmt)``;callable → ``f(v) -> str``。
    annot_fontsize : float, optional
        默认 ``rcParams['font.size'] - 1.5``(主体 8pt → annot 6.5pt)。
    annot_dark_threshold : float
        ``|z/scale| > threshold`` 时白字、否则黑字(scale = max(|vmin|, |vmax|))。
    figsize : (w, h), optional
        None 时按 nrows / ncols 自动算紧凑 figsize(``_auto_figsize``)。
    cbar_label : str
        色条标签。``cbar_pos`` 控制位置(figure 坐标)。
    title : str
        ``fig.suptitle`` 标题。
    xtick_rotation : float
    ytick_right : bool
        True → ytick 文字放右侧、隐藏左侧。
    xtick_fontsize, ytick_fontsize : float, optional
        None 时走 ``rcParams['xtick.labelsize']`` / ``ytick.labelsize``(默认 8pt)。
    row_colors_label : str
        给 row_colors 色带加上小标签(仅 Series 输入有效)。
    cbar_pos : (left, bottom, width, height)
        小色条位置(figure 坐标系)。
    dendrogram_ratio : (row_ratio, col_ratio)
        dendrogram 占比(传给 seaborn.clustermap)。
    return_grid : bool
        True → 返回 (fig, sns.ClusterGrid),便于二次调样。

    Returns
    -------
    fig : matplotlib.figure.Figure
        或 ``(fig, ClusterGrid)`` 当 ``return_grid=True``。

    Notes
    -----
    - 配合 ``drawing_yh.save_fig`` 一次性写 PNG/PDF/SVG。
    - 期刊投稿尺寸用 ``DOUBLE_COL_IN`` / ``SINGLE_COL_IN`` 约束:
      ``figsize=(SINGLE_COL_IN, h)``;否则默认按内容算。
    """
    if data_raw is None:
        data_raw = data_z
    if list(data_z.index) != list(data_raw.index) or list(data_z.columns) != list(data_raw.columns):
        raise ValueError("data_z 与 data_raw 的 index / columns 必须完全一致")

    n_rows, n_cols = data_z.shape

    # 按 STANDARDS:所有文字 = rcParams['font.size'](包级 8pt)
    # "所有文字一致(刻度、轴名、annotation、title)"
    base_font = float(mpl.rcParams.get("font.size", 8))
    if annot_fontsize is None:
        annot_fontsize = base_font
    if xtick_fontsize is None:
        xtick_fontsize = float(mpl.rcParams.get("xtick.labelsize", base_font))
    if ytick_fontsize is None:
        ytick_fontsize = float(mpl.rcParams.get("ytick.labelsize", base_font))

    # figsize 自动
    if figsize is None:
        figsize = _auto_figsize(n_rows, n_cols)

    # row_colors / col_colors 规范化(支持 array / Series / DataFrame)
    def _norm_colors(rc, axis_index, label):
        if rc is None or isinstance(rc, (pd.Series, pd.DataFrame)):
            if isinstance(rc, pd.Series) and label and rc.name is None:
                rc = rc.rename(label)
            return rc
        return pd.Series(list(rc), index=axis_index, name=label or None)

    row_colors = _norm_colors(row_colors, data_z.index, row_colors_label)
    col_colors = _norm_colors(col_colors, data_z.columns, "")

    # ---------- seaborn clustermap ----------
    g = sns.clustermap(
        data_z,
        method=linkage_method, metric=metric,
        cmap=cmap, vmin=vmin, vmax=vmax, center=center,
        row_cluster=row_cluster, col_cluster=col_cluster,
        row_colors=row_colors, col_colors=col_colors,
        figsize=figsize,
        xticklabels=True, yticklabels=True,
        cbar_pos=cbar_pos,
        dendrogram_ratio=dendrogram_ratio,
    )

    # reorder
    row_idx = g.dendrogram_row.reordered_ind if row_cluster else list(range(n_rows))
    col_idx = g.dendrogram_col.reordered_ind if col_cluster else list(range(n_cols))
    z_ord = data_z.iloc[row_idx, col_idx]
    r_ord = data_raw.iloc[row_idx, col_idx]

    ax = g.ax_heatmap

    # ---------- per-group cmap(可选):重画 heatmap 各列组 ----------
    # 每个 cell 的"色阶 scale"用于决定 annot 文字颜色(白/黑)。默认单 cmap 模式
    # scale 是全局 vmin/vmax 的绝对最大;per-group 模式 scale 按列查表。
    col_x_start = None  # col_logical -> 该列 imshow 起点 x(支持组间 gap)
    if col_cmap_groups is not None:
        # 计算每列的 x 起点(允许组间 gap)
        col_x_start = {}
        x_cursor = 0.0
        # 用 col_idx(reorder 后)走索引,但 col_cluster=False 时与 logical 一致
        for k, (cols, _) in enumerate(col_cmap_groups):
            if k > 0:
                x_cursor += col_group_gap
            for c in cols:
                pos = list(col_idx).index(c) if col_cluster else c
                col_x_start[pos] = x_cursor
                x_cursor += 1.0
        total_width = x_cursor

        # 移除 sns.clustermap 默认 QuadMesh + 默认 cbar
        for coll in list(ax.collections):
            coll.remove()
        if g.ax_cbar is not None:
            g.ax_cbar.remove()
        # per-column scale lookup(用于 annot 颜色)
        col_scale = np.ones(n_cols)
        group_specs = []  # (sub_vmin, sub_vmax, cmap_g, cols, label)
        for cols, cmap_g in col_cmap_groups:
            cols = list(cols)
            sub = z_ord.iloc[:, cols].values
            sub_vmin = float(sub.min())
            sub_vmax = float(sub.max())
            for j_logical in cols:
                pos = list(col_idx).index(j_logical) if col_cluster else j_logical
                x0 = col_x_start[pos]
                ax.imshow(z_ord.iloc[:, [pos]].values,
                          cmap=cmap_g, aspect="auto",
                          extent=[x0, x0 + 1, n_rows, 0],
                          vmin=sub_vmin, vmax=sub_vmax,
                          interpolation="nearest", zorder=1)
                col_scale[pos] = max(abs(sub_vmin), abs(sub_vmax)) or 1.0
            group_specs.append((sub_vmin, sub_vmax, cmap_g, cols,
                                ", ".join(str(data_z.columns[c]) for c in cols)))

        # ---------- 多 colorbar 排布 ----------
        n_groups = len(group_specs)
        cb_left, cb_bot, cb_w, cb_h = cbar_pos
        cb_label_fs = base_font   # STANDARDS:cbar label / tick 同 base_font
        for k, (sub_vmin, sub_vmax, cmap_g, cols, label) in enumerate(group_specs):
            if cbar_layout == "top-row":
                # row_dendrogram 之上的空闲区,横向并排
                # (避开 suptitle 在 y=0.985、避开 heatmap 顶端在 y≈0.88)
                bar_w = 0.16
                bar_h_pix = 0.013
                gap = 0.030
                # 从左到右排列;左起留 2% 边距
                left = 0.02 + k * (bar_w + gap)
                bar_ax = g.fig.add_axes([left, 0.935, bar_w, bar_h_pix])
                orient = "horizontal"
            elif cbar_layout == "bottom-row":
                bar_w = 0.11
                bar_h_pix = 0.012
                gap = 0.018
                right = 0.99
                left = right - (k + 1) * (bar_w + gap)
                bar_ax = g.fig.add_axes([left, 0.02, bar_w, bar_h_pix])
                orient = "horizontal"
            else:  # left-stack(老行为,兼容)
                bar_h_each = max(0.012, cb_h / max(n_groups, 1) * 0.55)
                bar_gap = bar_h_each * 0.6
                bar_ax = g.fig.add_axes([cb_left,
                                         cb_bot - k * (bar_h_each + bar_gap),
                                         max(cb_w, 0.10), bar_h_each])
                orient = "horizontal"
            sm = plt.cm.ScalarMappable(cmap=cmap_g,
                                       norm=mpl.colors.Normalize(vmin=sub_vmin, vmax=sub_vmax))
            cb = plt.colorbar(sm, cax=bar_ax, orientation=orient)
            # top-row / bottom-row 横向 cbar 把 label 放上方,避免跟主热图接近
            if cbar_layout in ("top-row", "bottom-row"):
                cb.ax.xaxis.set_label_position("top")
            cb.set_label(label, fontsize=cb_label_fs, labelpad=2)
            cb.ax.tick_params(labelsize=cb_label_fs)

        # 还原 ax_heatmap 范围 + 重设 xtick 到列中心(含 gap)
        ax.set_xlim(0, total_width); ax.set_ylim(n_rows, 0)
        # 把 xtick 移到带 gap 的列中心
        tick_pos = []
        tick_lbl = []
        for j_disp in range(n_cols):
            x0 = col_x_start[j_disp]
            tick_pos.append(x0 + 0.5)
            tick_lbl.append(str(z_ord.columns[j_disp]))
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_lbl)
    else:
        col_scale = None
        total_width = n_cols
        col_x_start = {j: float(j) for j in range(n_cols)}

    # 格内文本
    vmin_eff = vmin if vmin is not None else float(data_z.values.min())
    vmax_eff = vmax if vmax is not None else float(data_z.values.max())
    scale_global = max(abs(vmin_eff), abs(vmax_eff)) or 1.0

    def _fmt(v):
        if callable(annot_fmt):
            return annot_fmt(v)
        if isinstance(v, (int, np.integer)) and isinstance(annot_fmt, str) and ("d" in annot_fmt or "g" in annot_fmt):
            return format(int(v), annot_fmt)
        return format(float(v), annot_fmt)

    for i in range(z_ord.shape[0]):
        for j in range(z_ord.shape[1]):
            z = z_ord.iat[i, j]
            r = r_ord.iat[i, j]
            try:
                if pd.isna(z) or pd.isna(r):
                    continue
            except TypeError:
                pass
            s = col_scale[j] if col_scale is not None else scale_global
            norm = abs(z) / s
            color = "white" if norm > annot_dark_threshold else "black"
            x_text = col_x_start[j] + 0.5
            ax.text(x_text, i + 0.5, _fmt(r),
                    ha="center", va="center",
                    fontsize=annot_fontsize, color=color)

    # 样式
    if ytick_right:
        ax.tick_params(axis="y", labelleft=False, labelright=True, length=0)
        ax.yaxis.set_label_position("right")
    ax.tick_params(axis="x", labelsize=xtick_fontsize)
    ax.tick_params(axis="y", labelsize=ytick_fontsize)
    if xtick_rotation:
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(xtick_rotation)
            lbl.set_ha("right" if xtick_rotation else "center")
    for d_ax in (g.ax_row_dendrogram, g.ax_col_dendrogram):
        if d_ax is None:
            continue
        for coll in d_ax.collections:
            coll.set_linewidth(1.0)

    # ---------- row_groups 标签:dendrogram 重排后每段中心标群名 ----------
    if row_groups is not None and row_group_labels:
        rg_ord = (row_groups.iloc[row_idx].values
                  if isinstance(row_groups, pd.Series)
                  else np.asarray(row_groups)[row_idx])
        y_cursor = 0
        segments = []  # (gid, start_row, end_row_exclusive)
        for gid, gp in groupby(rg_ord):
            L = sum(1 for _ in gp)
            segments.append((gid, y_cursor, y_cursor + L))
            y_cursor += L
        # row_colors 色块的 figure 坐标;若无 row_colors 则用 row_dendrogram 旁
        anchor_ax = g.ax_row_colors if g.ax_row_colors is not None else g.ax_row_dendrogram
        if anchor_ax is not None:
            pos = anchor_ax.get_position()
            n_tot = len(rg_ord)
            for gid, s, e in segments:
                center_row = (s + e) / 2.0
                # seaborn ax_row_colors y 轴默认顶=0 → 底=n_rows
                fig_y = pos.y1 - (center_row / n_tot) * (pos.y1 - pos.y0)
                fig_x = pos.x0 - 0.005  # 略左于 row_colors
                g.fig.text(fig_x, fig_y,
                           str(row_group_labels.get(gid, "")),
                           ha="right", va="center", fontsize=base_font)

    if cbar_label and g.ax_cbar is not None:
        # STANDARDS:cbar label / tick 同 base_font
        g.ax_cbar.set_ylabel(cbar_label, fontsize=base_font)
        g.ax_cbar.tick_params(labelsize=base_font)

    if title:
        # cbar 在 top-row(y≈0.93)时把 title 顶到 0.985 避免重叠;其它布局 y=1.00
        # STANDARDS:title 用 base_font(不放大、不缩小)
        title_y = 0.985 if (col_cmap_groups is not None and cbar_layout == "top-row") else 1.00
        g.fig.suptitle(title, fontsize=base_font, y=title_y)

    if return_grid:
        return g.fig, g
    return g.fig


__all__ = ["clustermap_annot"]
