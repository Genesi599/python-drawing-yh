#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
clustermap_annot — 带格内值标注 + 行色块的聚类热图(通用版)。

基于 ``seaborn.clustermap`` 做底层聚类与渲染,在格内叠加任意"原始值"文本
(允许色阶用 Z-score / 归一化、文本用原始量),并支持行色块(功能群 /
谱系 / 类别注释)。drawing_yh 风格:yticklabel 右侧、紧凑色条、小字号、
彩色 row_colors。

典型场景:
  - 色 = Z-score(列方向),文本 = 原始 hit count / 表达均值
  - 色 = pearson r(范围 -1..1),文本 = 相关系数本身(`data_raw=None`)
  - 色 = -log10 p,文本 = effect size

API
---
    clustermap_annot(data_z, data_raw=None, *, row_colors=None, ...)
returns matplotlib.figure.Figure;`return_grid=True` 时返回 (fig, ClusterGrid)。

@File    : clustermap_annot.py
@Date    : 2026-06-04
@Author  : yh109
"""
from __future__ import annotations

from numbers import Number
from typing import Callable, Iterable, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def clustermap_annot(
    data_z: pd.DataFrame,
    data_raw: Optional[pd.DataFrame] = None,
    *,
    row_colors: Optional[Union[Sequence, pd.Series]] = None,
    col_colors: Optional[Union[Sequence, pd.Series]] = None,
    row_cluster: bool = True,
    col_cluster: bool = False,
    linkage_method: str = "ward",
    metric: str = "euclidean",
    cmap: str = "RdBu_r",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    center: Optional[float] = None,
    annot_fmt: Union[str, Callable[[float], str]] = ".2f",
    annot_fontsize: float = 6.5,
    annot_dark_threshold: float = 0.55,
    figsize: Tuple[float, float] = (7.5, 6.0),
    cbar_label: str = "",
    title: str = "",
    xtick_rotation: float = 0,
    ytick_right: bool = True,
    xtick_fontsize: float = 8,
    ytick_fontsize: float = 7,
    row_colors_label: str = "",
    cbar_pos: Tuple[float, float, float, float] = (0.02, 0.85, 0.045, 0.12),
    return_grid: bool = False,
) -> Union[plt.Figure, Tuple[plt.Figure, "sns.matrix.ClusterGrid"]]:
    """带格内值标注 + 行色块的聚类热图。

    Parameters
    ----------
    data_z : pd.DataFrame
        驱动**色阶**的矩阵(行 = items,列 = features)。常用 Z-score、相关 r、log2FC 等。
    data_raw : pd.DataFrame, optional
        驱动**格内文本**的同 shape 矩阵。若为 None 直接用 ``data_z`` 的数值。
        必须与 ``data_z`` 同 index、同 columns。
    row_colors : array-like or pd.Series, optional
        长度 = ``len(data_z)``,给左侧色块(功能群 / 谱系 / 实验分组等)。
        seaborn 把它和 dendrogram 一起画;若是 pd.Series 会用它的 index 对齐。
    col_colors : 同上,列方向(可选)。
    row_cluster, col_cluster : bool
        默认行聚类、列不聚类(适用于"列 = 固定 features"的场景)。
    linkage_method, metric : str
        传给 ``scipy.cluster.hierarchy.linkage``。
    cmap, vmin, vmax, center : matplotlib heatmap 标准参数。
    annot_fmt : str or callable
        格内文本 format。字符串走 ``format(value, annot_fmt)``;callable
        签名为 ``f(value) -> str``。例:``".2f"`` / ``"{:d}"`` / ``lambda x: f"{int(x)}*"``
        (后者无效,这里只接受 format spec 或函数)。
    annot_fontsize : float
    annot_dark_threshold : float
        ``|normalized_z| > threshold`` 时用白字,否则黑字。normalized = ``|z| /
        max(|vmin|, |vmax|)``,范围 0–1。
    figsize : (w, h)
    cbar_label, title : str
    xtick_rotation : float
    ytick_right : bool
        True → ytick 文字放右侧 + 隐藏左侧。
    xtick_fontsize, ytick_fontsize : float
    row_colors_label : str
        在 row_colors 色块上方加一个轴标题。
    cbar_pos : (left, bottom, width, height) figure 坐标系
        小色条位置;默认放左上角。
    return_grid : bool
        True 时同时返回 seaborn ClusterGrid(便于二次调样式)。

    Returns
    -------
    fig : matplotlib.figure.Figure
    (fig, ClusterGrid) if ``return_grid=True``.

    Notes
    -----
    - 配合 ``drawing_yh.save_fig`` 一次性写 PNG/PDF/SVG。
    - dendrogram 默认线宽 1.0,seaborn 默认偏细;如需更粗在 g.ax_row_dendrogram /
      g.ax_col_dendrogram.collections 上手动调。
    """
    if data_raw is None:
        data_raw = data_z
    if list(data_z.index) != list(data_raw.index) or list(data_z.columns) != list(data_raw.columns):
        raise ValueError("data_z 与 data_raw 的 index / columns 必须完全一致")

    # row_colors: 对齐成与 data_z.index 同序的 Series
    if row_colors is not None and not isinstance(row_colors, pd.Series):
        row_colors = pd.Series(list(row_colors), index=data_z.index, name=row_colors_label or None)
    elif isinstance(row_colors, pd.Series) and row_colors_label:
        row_colors = row_colors.rename(row_colors_label)
    if col_colors is not None and not isinstance(col_colors, pd.Series):
        col_colors = pd.Series(list(col_colors), index=data_z.columns)

    # ---------- seaborn clustermap ----------
    g = sns.clustermap(
        data_z,
        method=linkage_method,
        metric=metric,
        cmap=cmap,
        vmin=vmin, vmax=vmax, center=center,
        row_cluster=row_cluster, col_cluster=col_cluster,
        row_colors=row_colors, col_colors=col_colors,
        figsize=figsize,
        xticklabels=True, yticklabels=True,
        cbar_pos=cbar_pos,
        dendrogram_ratio=(0.18, 0.10),
    )

    # ---------- reorder 后的真实矩阵 ----------
    row_idx = g.dendrogram_row.reordered_ind if row_cluster else list(range(len(data_z)))
    col_idx = g.dendrogram_col.reordered_ind if col_cluster else list(range(len(data_z.columns)))
    z_ord = data_z.iloc[row_idx, col_idx]
    r_ord = data_raw.iloc[row_idx, col_idx]

    # ---------- 格内文本 ----------
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

    # ---------- 样式 ----------
    if ytick_right:
        ax.tick_params(axis="y", labelleft=False, labelright=True, length=0)
        ax.yaxis.set_label_position("right")
    ax.tick_params(axis="x", labelsize=xtick_fontsize)
    ax.tick_params(axis="y", labelsize=ytick_fontsize)
    if xtick_rotation:
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(xtick_rotation)
            lbl.set_ha("right" if xtick_rotation else "center")
    # 加粗 dendrogram 线
    for d_ax in (g.ax_row_dendrogram, g.ax_col_dendrogram):
        if d_ax is None:
            continue
        for coll in d_ax.collections:
            coll.set_linewidth(1.0)

    if cbar_label:
        g.ax_cbar.set_ylabel(cbar_label, fontsize=7)
        g.ax_cbar.tick_params(labelsize=7)

    if title:
        g.fig.suptitle(title, fontsize=9.5, y=1.00)

    if return_grid:
        return g.fig, g
    return g.fig


__all__ = ["clustermap_annot"]
