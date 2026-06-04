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

from typing import Callable, Optional, Sequence, Tuple, Union

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
    row_cluster: bool = True,
    col_cluster: bool = False,
    linkage_method: str = "ward",
    metric: str = "euclidean",
    cmap: str = "RdBu_r",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    center: Optional[float] = None,
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
    row_cluster, col_cluster : bool
        默认行聚类、列不聚类(列 = 固定 features)。
    linkage_method, metric : str
        ``scipy.cluster.hierarchy.linkage`` 参数。
    cmap, vmin, vmax, center : matplotlib heatmap 标准参数。
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

    # 字号默认走 rcParams(包级 8pt)
    base_font = float(mpl.rcParams.get("font.size", 8))
    if annot_fontsize is None:
        annot_fontsize = max(5.0, base_font - 1.5)
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

    # 格内文本
    vmin_eff = vmin if vmin is not None else float(data_z.values.min())
    vmax_eff = vmax if vmax is not None else float(data_z.values.max())
    scale = max(abs(vmin_eff), abs(vmax_eff)) or 1.0

    def _fmt(v):
        if callable(annot_fmt):
            return annot_fmt(v)
        if isinstance(v, (int, np.integer)) and isinstance(annot_fmt, str) and ("d" in annot_fmt or "g" in annot_fmt):
            return format(int(v), annot_fmt)
        return format(float(v), annot_fmt)

    ax = g.ax_heatmap
    for i in range(z_ord.shape[0]):
        for j in range(z_ord.shape[1]):
            z = z_ord.iat[i, j]
            r = r_ord.iat[i, j]
            try:
                if pd.isna(z) or pd.isna(r):
                    continue
            except TypeError:
                pass
            norm = abs(z) / scale
            color = "white" if norm > annot_dark_threshold else "black"
            ax.text(j + 0.5, i + 0.5, _fmt(r),
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

    if cbar_label and g.ax_cbar is not None:
        g.ax_cbar.set_ylabel(cbar_label, fontsize=max(5.0, base_font - 1))
        g.ax_cbar.tick_params(labelsize=max(5.0, base_font - 1))

    if title:
        g.fig.suptitle(title, fontsize=base_font + 1.5, y=1.00)

    if return_grid:
        return g.fig, g
    return g.fig


__all__ = ["clustermap_annot"]
