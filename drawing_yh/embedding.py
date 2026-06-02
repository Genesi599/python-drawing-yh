# -*- coding: utf-8 -*-
"""
Embedding feature plot — unified template (pure matplotlib).

Colour each cell on a 2-D embedding (UMAP / t-SNE / any obsm) by a gene's
expression: one small panel per gene, a shared grey->red colour scale with a
single colour-bar, and a small UMAP1/UMAP2 arrow key in the figure corner
(replacing per-panel ticks).

This is the **drawing layer only** — it takes plain arrays (coords + a value
matrix), NOT an AnnData. The single-cell data prep (read h5ad, pick markers,
pull coords / values) lives in ``single_cell-yh`` and calls ``feature_plot``
here, so ``drawing_yh`` stays free of scanpy/anndata.

Primitives (compose your own grid):

* ``resolve_vlim`` — percentile (``"p2"``/``"p98"``) or absolute colour limits.
* ``scatter_embedding`` — one panel: value-sorted scatter, no ticks.
* ``add_embedding_axes`` — the corner UMAP1/UMAP2 arrow key.
* ``auto_ncols`` — pick a column count for n panels.

One-call: ``feature_plot(coords, values, ...) -> (fig, axes)``.
"""
from __future__ import annotations

from collections.abc import Mapping

import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Colormap, LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

from . import DEFAULT_FONT_SIZE


# grey -> red, the scanpy-style single-channel expression ramp
GREY_RED = ("#F0F0F0", "#B2182B")


# ============================================================
# primitives
# ============================================================
def _resolve_cmap(cmap):
    """Accept a Colormap, a colormap name, or a sequence of colours (turned
    into a LinearSegmentedColormap, e.g. ``("#F0F0F0", "#B2182B")``)."""
    if isinstance(cmap, Colormap):
        return cmap
    if isinstance(cmap, str):
        return plt.get_cmap(cmap)
    return LinearSegmentedColormap.from_list("feature_cmap", list(cmap))


def resolve_vlim(values, vmin="p2", vmax="p98", *, nonzero: bool = True):
    """Resolve colour limits.

    ``vmin`` / ``vmax`` may each be a percentile string (``"p10"``), an
    absolute float, or ``None`` (data min / max). When ``nonzero`` (default),
    percentiles use the > 0 values — sparse snRNA expression is mostly zeros,
    so percentiles over all cells would collapse the scale. Falls back to all
    finite values when everything is 0.
    """
    arr = np.asarray(values, dtype=float).ravel()
    arr = arr[np.isfinite(arr)]
    base = arr[arr > 0] if (nonzero and np.any(arr > 0)) else arr
    if base.size == 0:
        base = np.array([0.0])

    def _num(v):
        if isinstance(v, str):
            return float(np.percentile(base, float(v.lower().lstrip("p"))))
        return float(v)

    lo = _num(vmin) if vmin is not None else float(base.min())
    hi = _num(vmax) if vmax is not None else float(base.max())
    if hi <= lo:
        hi = lo + 1e-9
    return lo, hi


def auto_ncols(n_panels: int, *, preferred_ratio: float = 1.2, max_ncols: int = 6) -> int:
    """Pick a column count keeping the grid close to ``preferred_ratio``
    (ncols / nrows) while minimising empty cells."""
    if n_panels <= 0:
        return 1
    best_ncols, best_score = 1, float("inf")
    for ncols in range(1, max_ncols + 1):
        nrows = int(np.ceil(n_panels / ncols))
        ratio_penalty = abs(ncols / nrows - preferred_ratio)
        empty_penalty = (nrows * ncols - n_panels) * 0.1
        score = ratio_penalty + empty_penalty
        if score < best_score:
            best_score, best_ncols = score, ncols
    return best_ncols


def scatter_embedding(ax, coords, values, *, cmap, vmin, vmax,
                      point_size: float = 6.0, sort_by_value: bool = True,
                      **scatter_kw):
    """Draw one embedding panel coloured by ``values``.

    High-value cells are plotted last (on top) so they are not hidden under the
    grey background; pass ``sort_by_value=False`` to keep input order. Removes
    ticks. Returns the ``PathCollection`` (for a colour-bar).
    """
    coords = np.asarray(coords, dtype=float)
    values = np.asarray(values, dtype=float).ravel()
    order = np.argsort(values, kind="stable") if sort_by_value else np.arange(len(values))
    sc = ax.scatter(
        coords[order, 0], coords[order, 1],
        c=values[order], cmap=cmap, vmin=vmin, vmax=vmax,
        s=point_size, linewidths=0, **scatter_kw,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    return sc


def add_embedding_axes(fig, axis_labels=("UMAP1", "UMAP2"), *,
                       start=(0.08, 0.05), length: float = 0.05,
                       font: int = DEFAULT_FONT_SIZE, lw: float = 1.2):
    """Small corner arrow key (x = ``axis_labels[0]``, y = ``axis_labels[1]``)
    in figure coordinates, replacing per-panel axis ticks for a clean
    multi-panel embedding grid. The y arrow is scaled by the figure aspect so
    both arrows look the same length on screen.
    """
    x0, y0 = start
    fw, fh = fig.get_size_inches()
    aspect = fw / fh
    lx, ly = length, length * aspect
    for dx, dy in ((lx, 0.0), (0.0, ly)):
        fig.patches.append(mpatches.FancyArrowPatch(
            (x0, y0), (x0 + dx, y0 + dy),
            arrowstyle="->", mutation_scale=10, color="black",
            lw=lw, transform=fig.transFigure,
        ))
    fig.text(x0 + lx + 0.004, y0 - 0.006, axis_labels[0],
             fontsize=font, ha="left", va="top")
    fig.text(x0 - 0.006, y0 + ly + 0.004, axis_labels[1],
             fontsize=font, ha="right", va="bottom", rotation=90)


# ============================================================
# one-call convenience
# ============================================================
def _prepare_panels(values, genes, titles):
    """Normalise ``values`` into ``(panel_titles, gene_names, matrix N×G)``.

    Accepts a dict ``{title: (N,) array}``, a DataFrame (columns = genes), or a
    2-D / 1-D array (with optional ``genes`` names).
    """
    if isinstance(values, Mapping):
        keys = list(values.keys())
        mat = np.column_stack([np.asarray(values[k], dtype=float).ravel() for k in keys])
        gene_names = list(genes) if genes is not None else keys
        panel_titles = list(titles) if titles is not None else keys
        return panel_titles, gene_names, mat
    if isinstance(values, pd.DataFrame):
        cols = list(genes) if genes is not None else list(values.columns)
        mat = values[cols].to_numpy(dtype=float)
        panel_titles = list(titles) if titles is not None else cols
        return panel_titles, cols, mat
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 1:
        arr = arr[:, None]
    n_gene = arr.shape[1]
    gene_names = list(genes) if genes is not None else [f"gene{i + 1}" for i in range(n_gene)]
    panel_titles = list(titles) if titles is not None else list(gene_names)
    return panel_titles, gene_names, arr


def feature_plot(
    coords,
    values,
    *,
    genes=None,
    titles=None,
    ncols="auto",
    max_ncols: int = 6,
    preferred_ratio: float = 1.2,
    cmap=GREY_RED,
    vmin="p2",
    vmax="p98",
    share_clim: bool = True,
    nonzero_percentile: bool = True,
    point_size: float = 6.0,
    sort_by_value: bool = True,
    panel_size: float = 1.7,
    figsize: tuple[float, float] | None = None,
    axis_labels=("UMAP1", "UMAP2"),
    show_axes_arrow: bool = True,
    gene_corner_label: bool = True,
    font: int = DEFAULT_FONT_SIZE,
    title_fontsize: int | None = None,
    gene_fontsize: int | None = None,
    show_spines: bool = True,
    cbar_label: str = "Expression",
):
    """Multi-gene embedding feature plot in the unified template style.

    Parameters
    ----------
    coords : (N, 2) array-like
        Embedding coordinates (UMAP / t-SNE / any obsm), one row per cell.
    values : dict | DataFrame | array
        Per-cell expression. ``{title: (N,) values}`` dict, a DataFrame whose
        columns are genes, or a 2-D ``(N, G)`` / 1-D ``(N,)`` array (name the
        columns via ``genes``).
    genes, titles
        Optional gene names (italic corner labels) and panel titles. Default:
        the dict keys / DataFrame columns.
    ncols, max_ncols, preferred_ratio
        Grid layout. ``"auto"`` picks a column count via :func:`auto_ncols`.
    cmap
        Colormap, colormap name, or a ``(low, high)`` colour pair. Default
        grey -> red.
    vmin, vmax, nonzero_percentile
        Colour limits (percentile string / float / None), resolved by
        :func:`resolve_vlim`.
    share_clim
        If True (default), all panels share one colour scale + a single
        colour-bar. If False, each panel is scaled independently with its own
        small colour-bar.
    point_size, sort_by_value
        Scatter marker area and whether to draw high-value cells on top.
    panel_size
        Inches per panel when ``figsize`` is None.
    axis_labels, show_axes_arrow
        Corner UMAP1/UMAP2 arrow key (see :func:`add_embedding_axes`).
    gene_corner_label
        Italic gene name in each panel's top-right corner.

    Returns
    -------
    (fig, axes)
        The figure and the list of per-gene axes (length n genes).
    """
    coords = np.asarray(coords, dtype=float)
    panel_titles, gene_names, mat = _prepare_panels(values, genes, titles)
    n = mat.shape[1]
    if n == 0:
        raise ValueError("no genes / value columns to plot")
    if mat.shape[0] != coords.shape[0]:
        raise ValueError(
            f"coords has {coords.shape[0]} cells but values has {mat.shape[0]}"
        )

    cmap_obj = _resolve_cmap(cmap)
    title_fs = title_fontsize if title_fontsize is not None else font
    gene_fs = gene_fontsize if gene_fontsize is not None else font

    ncols_eff = (auto_ncols(n, preferred_ratio=preferred_ratio, max_ncols=max_ncols)
                 if ncols == "auto" else max(1, int(ncols)))
    nrows = int(np.ceil(n / ncols_eff))

    g_vmin = g_vmax = None
    if share_clim:
        g_vmin, g_vmax = resolve_vlim(mat.ravel(), vmin, vmax, nonzero=nonzero_percentile)

    if figsize is None:
        cbar_w = 0.4 if share_clim else 0.0
        figsize = (panel_size * ncols_eff + cbar_w, panel_size * nrows)

    fig = plt.figure(figsize=figsize)
    n_grid_cols = ncols_eff + (1 if share_clim else 0)
    width_ratios = [1] * ncols_eff + ([0.06] if share_clim else [])
    gs = GridSpec(nrows, n_grid_cols, figure=fig,
                  width_ratios=width_ratios, wspace=0.08, hspace=0.14)

    axes = []
    last_sc = None
    for i in range(n):
        r, c = divmod(i, ncols_eff)
        ax = fig.add_subplot(gs[r, c])
        axes.append(ax)
        if share_clim:
            pv, px = g_vmin, g_vmax
        else:
            pv, px = resolve_vlim(mat[:, i], vmin, vmax, nonzero=nonzero_percentile)
        sc = scatter_embedding(ax, coords, mat[:, i], cmap=cmap_obj,
                               vmin=pv, vmax=px, point_size=point_size,
                               sort_by_value=sort_by_value)
        last_sc = sc
        ax.tick_params(left=False, bottom=False)
        if not show_spines:
            for sp in ax.spines.values():
                sp.set_visible(False)
        if panel_titles[i]:
            ax.set_title(str(panel_titles[i]), fontsize=title_fs, pad=2)
        if gene_corner_label:
            ax.text(0.97, 0.97, str(gene_names[i]), transform=ax.transAxes,
                    fontsize=gene_fs, fontstyle="italic", color="black",
                    va="top", ha="right")
        if not share_clim:
            cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
            cb.ax.tick_params(labelsize=font)

    for j in range(n, nrows * ncols_eff):
        r, c = divmod(j, ncols_eff)
        fig.add_subplot(gs[r, c]).axis("off")

    if share_clim and last_sc is not None:
        cax = fig.add_subplot(gs[:, -1])
        cb = fig.colorbar(last_sc, cax=cax)
        cb.set_label(cbar_label, fontsize=font)
        cb.ax.tick_params(labelsize=font)

    if show_axes_arrow:
        add_embedding_axes(fig, axis_labels, font=font)

    return fig, axes


# ============================================================
# categorical cluster atlas (脑膜 complete-atlas 画法模版)
# ============================================================
def cluster_centroids(coords, labels, *, trim: float = 0.65):
    """每个类别的稳健质心: 先中位, 再取距质心最近 ``trim`` 分位的点重算中位
    (抗离群细胞 / 触角丝须, 标号落在团主体)。返回 ``{label: (x, y, n)}``。
    """
    coords = np.asarray(coords, dtype=float)
    labels = np.asarray(labels)
    out = {}
    for lab in pd.unique(labels):
        m = labels == lab
        x, y = coords[m, 0], coords[m, 1]
        if x.size == 0:
            continue
        cx, cy = np.median(x), np.median(y)
        d = (x - cx) ** 2 + (y - cy) ** 2
        keep = d <= np.quantile(d, trim)
        out[str(lab)] = (float(np.median(x[keep])), float(np.median(y[keep])), int(m.sum()))
    return out


def _atlas_palette(n):
    """n 个类别的定性配色 (tab20 → tab20b → tab20c 循环, 不够再用 hsv 补)。"""
    import matplotlib.cm as cm
    import colorsys
    base = []
    for name in ("tab20", "tab20b", "tab20c"):
        base.extend(list(cm.get_cmap(name).colors))
    if n <= len(base):
        return base[:n]
    extra = [colorsys.hsv_to_rgb((i / max(1, n - len(base))) % 1.0, 0.6, 0.85)
             for i in range(n - len(base))]
    return base + extra


def _draw_atlas_key(ax, order, label_ids, color_of, disp, sublabels, count_of, ncol, fs):
    """右侧独立编号 legend: 2(或 ncol)列, 每项 = 色点 + "编号. 显示名 (计数)" + 可选副标。"""
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    rows = int(np.ceil(len(order) / ncol))
    for idx, lab in enumerate(order):
        col, row = idx // rows, idx % rows
        x0 = 0.02 + col * (1.0 / ncol)
        y = 1.0 - (row + 0.5) / rows
        ax.scatter([x0], [y], s=22, color=color_of[lab], transform=ax.transAxes, clip_on=False)
        txt = f"{label_ids[lab]}. {disp(lab)}"
        if count_of is not None:
            txt += f" ({count_of.get(lab, 0):,})"
        ax.text(x0 + 0.035, y, txt, transform=ax.transAxes, ha="left", va="center", fontsize=fs)
        sub = (sublabels.get(lab) if isinstance(sublabels, Mapping) else None) if sublabels else None
        if sub:
            ax.text(x0 + 0.035, y - 0.45 / rows, str(sub), transform=ax.transAxes,
                    ha="left", va="top", fontsize=fs - 1.2, color="#666666")


def cluster_atlas(
    coords, labels, *,
    order=None, colors=None, excluded=None, label_ids=None,
    title=None, display=None, sublabels=None,
    counts_in_legend=True, centroid="number", min_centroid_cells=0,
    label_offsets=None, legend=True, legend_ncol=2, legend_fontsize=None,
    point_size=0.75, bg_color="#BDBDBD", bg_size=0.06, bg_alpha=0.18,
    alpha=0.78, square=True, axis_labels=("UMAP 1", "UMAP 2"),
    font=DEFAULT_FONT_SIZE, fig_size=None, ax=None,
):
    """分类细胞图谱 UMAP 一键模版 (额叶 final complete-atlas 画法)。

    每类一离散色, **按细胞数 z-order**(大群在底、小群在上不被盖), 团心标编号(白描边),
    右侧 2 列编号 legend(编号→显示名 + 计数), 剔除细胞画淡灰底, 正方形等比(square)
    + UMAP1/2 轴标。分类图谱版一键函数(对应连续值的 ``feature_plot``)。

    coords : (N, 2) x/y 坐标。  labels : (N,) 每个细胞的类别 label。
    order : 编号 / legend 顺序; ``None`` → 按细胞数降序。(画图层序始终大群在底, 与此无关)
    colors : dict ``label->color`` 或与 order 对齐的序列; ``None`` → 自动定性配色。
    excluded : (N,) bool, ``True`` 的细胞画淡灰底, 不参与编号 / 质心 / legend。
    label_ids : dict ``label->编号``; ``None`` → 按 order 自动 1..k。传入则沿用(跨图谱共用全局编号)。
    display : dict 或 callable ``label->显示名``; ``None`` → 原样(可传 celltype_fullnames)。
    sublabels : dict ``label->副标``(如大类 / 谱系), legend 名下加一行灰字。
    centroid : ``"number"``(默认) | ``"name"`` | ``None``。
    min_centroid_cells : 细胞数 < 此值的类不标团心(默认 0 = 全标; 额叶用 80 跳过 tiny)。
    label_offsets : dict ``label->(dx, dy)`` 手动微调团心标号位置(防数字相撞)。
    point_size : 散点大小(默认 0.75 固定; 小群靠 z-order 在上层仍可见)。
    square : ``True`` → data 范围正方形居中 + box_aspect=1。
    返回 ``(fig, ax, label_ids)``, ``label_ids`` = ``{label: 编号}``。
    """
    coords = np.asarray(coords, dtype=float)
    labels = np.asarray(labels).astype(str)
    excluded = (np.zeros(len(labels), bool) if excluded is None
                else np.asarray(excluded, dtype=bool))
    kept = ~excluded
    count_of = pd.Series(labels[kept]).value_counts().to_dict()
    label_offsets = label_offsets or {}

    if order is None:
        order = list(pd.Series(labels[kept]).value_counts().index)
    else:
        order = [str(o) for o in order if (labels[kept] == str(o)).any()]
    # 编号: 默认按 order 顺序 1..k; 传 label_ids 则用之(跨图谱沿用全局编号)
    if label_ids is None:
        label_ids = {lab: i + 1 for i, lab in enumerate(order)}
    else:
        label_ids = {lab: label_ids[lab] for lab in order if lab in label_ids}

    if colors is None:
        pal = _atlas_palette(len(order))
        color_of = {lab: pal[i] for i, lab in enumerate(order)}
    elif isinstance(colors, Mapping):
        color_of = {lab: colors.get(lab, "#9A9A9A") for lab in order}
    else:
        color_of = {lab: colors[i] for i, lab in enumerate(order)}

    def disp(lab):
        if display is None:
            return lab
        return display(lab) if callable(display) else display.get(lab, lab)

    if ax is None and legend:
        nrow = int(np.ceil(len(order) / legend_ncol))
        max_lab = max((len(str(disp(l))) for l in order), default=8)
        key_w = 0.55 + (0.052 * max_lab + 0.45) * legend_ncol
        main_w = fig_size[0] if fig_size else 5.2
        fh = fig_size[1] if fig_size else max(4.4, 0.155 * nrow + 1.7)
        fig, (ax, key_ax) = plt.subplots(
            1, 2, figsize=(main_w + key_w, fh),
            gridspec_kw={"width_ratios": [main_w, key_w], "wspace": 0.02})
    elif ax is None:
        fig, ax = plt.subplots(figsize=fig_size or (5.4, 5.0))
        key_ax = None
    else:
        fig, key_ax = ax.figure, None

    # 剔除细胞: 淡灰背景, 先画
    if excluded.any():
        ax.scatter(coords[excluded, 0], coords[excluded, 1], s=bg_size,
                   c=bg_color, linewidths=0, alpha=bg_alpha, rasterized=True)
    # 类别散点: 大群先画(底层 z 小), 小群在上(不被盖)
    draw_order = sorted(order, key=lambda lab: -count_of.get(lab, 0))
    for z, lab in enumerate(draw_order, start=2):
        m = kept & (labels == lab)
        if not m.any():
            continue
        ax.scatter(coords[m, 0], coords[m, 1], s=point_size, color=color_of[lab],
                   linewidths=0, alpha=alpha, rasterized=True, zorder=z)
    # 团心编号 (白描边; 跳过 tiny; 可手动微调位置)
    if centroid:
        for lab, (cx, cy, n) in cluster_centroids(coords[kept], labels[kept]).items():
            if lab not in label_ids or n < min_centroid_cells:
                continue
            dx, dy = label_offsets.get(lab, (0.0, 0.0))
            txt = str(label_ids[lab]) if centroid == "number" else disp(lab)
            ax.text(cx + dx, cy + dy, txt, ha="center", va="center", fontsize=font,
                    fontweight="bold", zorder=100,
                    path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])
    if square:
        xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
        ymin, ymax = coords[:, 1].min(), coords[:, 1].max()
        xmid, ymid = (xmin + xmax) / 2, (ymin + ymax) / 2
        span = max(xmax - xmin, ymax - ymin) * 1.04
        ax.set_xlim(xmid - span / 2, xmid + span / 2)
        ax.set_ylim(ymid - span / 2, ymid + span / 2)
        ax.set_box_aspect(1)
    ax.set_aspect("equal", adjustable="box")
    if title:
        ax.set_title(title, fontsize=font, pad=4)
    ax.set_xlabel(axis_labels[0], fontsize=font)
    ax.set_ylabel(axis_labels[1], fontsize=font)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    if key_ax is not None:
        _draw_atlas_key(key_ax, order, label_ids, color_of, disp, sublabels,
                        count_of if counts_in_legend else None,
                        legend_ncol, legend_fontsize or (font - 1.2))
    return fig, ax, label_ids


__all__ = [
    "GREY_RED",
    "resolve_vlim",
    "auto_ncols",
    "scatter_embedding",
    "add_embedding_axes",
    "feature_plot",
    "cluster_centroids",
    "cluster_atlas",
]
