#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Marker expression dot plot — unified template.

The default colour ramp is grey->red, matching single_cell-yh (scseq) dot
plots; pass ``cmap="viridis"`` for the brain-aging "Olig panel" look. Per-gene
min-max scaled mean colour (vmin=0 / vmax=1). By default both keys sit at the
right, stacked and vertically aligned (``legend_loc="right"``): a narrow
colour-bar plus a 10/30/60% dot-size key, each with its name as vertical text
aligned in one right-hand column. dot size = pct expressed (thin dark edge
``#303030`` / lw 0.12); ``legend_loc="bottom"`` falls back to a right colour-
bar + bottom-centre horizontal size key. Generous padding so the largest edge
dots are never clipped, and optional vertical block dividers between functional
/ per-subtype gene groups.

Two ways to use the module:

* **Primitives** — call ``dot_sizes`` / ``pad_axes`` / ``style_colorbar`` /
  ``add_size_legend`` / ``add_block_separators`` from inside an existing
  ``plt.subplots`` / ``ax.scatter`` flow. Use this when you need fine control
  over fig size, row highlights, custom colour-bar axes, etc.

* **One-call convenience** — ``marker_dotplot(data, row_order, gene_order,
  ...)`` for simple panels. Returns ``(fig, ax, sc)``.

Both options enforce the same style; the package-level rcParams give Arial 8pt
+ editable fonts on import.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import blended_transform_factory

from . import DEFAULT_FONT_SIZE


# ============================================================
# style constants (single source of truth)
# ============================================================
DOT_CMAP = "viridis"   # 备选连续色板(传 cmap="viridis" 启用)
# 默认 dot 配色:grey -> red,与 single_cell-yh(scseq) 原 dot plot 一致
GREY_RED_CMAP = LinearSegmentedColormap.from_list("grey_to_red", ["#F0F0F0", "#B2182B"])
DOT_EDGE = "#303030"
DOT_EDGE_LW = 0.12
DOT_LEGEND_PCTS = (10, 30, 60)
COLORBAR_LABEL = "Gene-scaled mean"
SEPARATOR_COLOR = "#606060"
SEPARATOR_LW = 0.55
SEPARATOR_ALPHA = 0.65
GRID_COLOR = "#E5E5E5"
GRID_LW = 0.45
GRID_ALPHA = 0.78
HIGHLIGHT_COLOR = "#F2E7C9"  # row band for caveat / sensitivity rows

# layout inch budgets —— 内容驱动 figsize + 边距(8pt Arial)
# 网格 pitch 实算 = max(单行文字含行距, 最大点直径) + 极小 gap → 文字与点都刚好不
# 重叠的最紧值,并随字体 / 数据最大点自适应;边距按标签长度、图例按固定英寸带:全部
# 用英寸算好再换成 figure 分数,随 figsize 与内容自适应、不写死分数(矩形网格版"尽量
# 紧凑";autoshrink_figsize 只缩正方形,不适合本图型,故用 pitch 预算等价实现)。
_CHAR_W_IN = 0.052        # 8pt Arial 单字符近似宽
_GRID_GAP_IN = 0.02       # 相邻格子最小净间隙(pitch = max(文字, 点) + gap)
_TEXT_PITCH_K = 1.30      # 单行文字含行距 ≈ font_pt * k / 72(竖排基因名 / 水平细胞名共用此下限)
_PITCH_MIN_IN = 0.18      # pitch 下限(极小字体 / 极小点时兜底)
_AXIS_LABEL_IN = 0.22     # x/y 轴名(旋转)占位
_ROW_DOT_IN = 0.16        # 左侧 cell-type 彩点占位
_PAD_IN = 0.05
_TITLE_IN = 0.28
_RIGHT_LEGEND_IN = 1.00   # 右侧竖排图例带(colorbar 柱 + 数字 + 竖排名称)
_BOTTOM_LEGEND_IN = 0.80  # 底部横排 size 图例(legend_loc='bottom')


# ============================================================
# primitives
# ============================================================
def dot_sizes(pct, *, scale: float = 2.4, base: float = 4.0):
    """Pct expressed (0-100) -> matplotlib scatter ``s`` (point area).

    Default ``scale=2.4`` / ``base=4.0`` matches the original 10/30/60% size
    key (28 / 76 / 148 pt^2). For sparse-expression panels (e.g. snRNA
    secreted factors where ct-level pct usually < 5%), pass a larger
    ``scale`` so the small dots remain readable, and re-key the legend
    via ``add_size_legend(pcts=...)``.
    """
    return np.clip(np.asarray(pct, dtype=float), 0, 100) * scale + base


def pad_axes(ax, n_cols, n_rows, pad: float = 0.6):
    """Fixed margin around the dot grid so edge dots (top/bottom rows, first/
    last columns) are never clipped, and y is inverted (row 0 at top).

    Use this instead of ``ax.invert_yaxis()`` / matplotlib's auto-margin —
    auto-margin is a fraction of the data range, which collapses to almost
    zero for few-row panels and clips the largest edge dots.
    """
    ax.set_xlim(-pad, (n_cols - 1) + pad)
    ax.set_ylim((n_rows - 1) + pad, -pad)


def style_colorbar(cbar, font: int = DEFAULT_FONT_SIZE, label: str = COLORBAR_LABEL):
    """Label + tick-fontsize the colour-bar in the unified style."""
    cbar.set_label(label, fontsize=font)
    cbar.ax.tick_params(labelsize=font)
    return cbar


def add_size_legend(ax, font: int = DEFAULT_FONT_SIZE, y: float = 0.012,
                    title: str = "Pct expressed",
                    pcts: Sequence[int] | None = None,
                    size_scale: float = 2.4, size_base: float = 4.0,
                    *, orientation: str = "horizontal", bbox=None):
    """Figure-bottom-centre 'Pct expressed' size key (default 10/30/60%) that
    matches the dot sizes from ``dot_sizes``. Uses ``fig.legend`` so the key
    sits at a fixed figure-bottom location regardless of axes height. Leave
    room with ``fig.subplots_adjust(bottom=...)``.

    For sparse panels, override the breakpoints + size mapping, e.g.
    ``add_size_legend(pcts=(1, 5, 15), size_scale=12, size_base=16)``; pair
    with the same ``size_scale`` / ``size_base`` in the actual ``dot_sizes``
    call (or in ``marker_dotplot(size_pcts=..., size_scale=..., size_base=...)``).
    """
    fig = ax.figure
    legend_pcts = list(pcts) if pcts else list(DOT_LEGEND_PCTS)
    handles = [
        ax.scatter([], [], s=float(dot_sizes(p, scale=size_scale, base=size_base)),
                   color="#888888", edgecolors=DOT_EDGE, linewidths=DOT_EDGE_LW,
                   label=f"{p}%")
        for p in legend_pcts
    ]
    if orientation == "vertical":
        anchor = bbox if bbox is not None else (0.86, 0.22)
        return fig.legend(
            handles=handles, title=title, loc="center left",
            bbox_to_anchor=anchor, ncol=1, frameon=False,
            fontsize=font, title_fontsize=font, labelspacing=1.3,
        )
    anchor = bbox if bbox is not None else (0.5, y)
    return fig.legend(
        handles=handles, title=title, loc="lower center",
        bbox_to_anchor=anchor, ncol=len(legend_pcts), frameon=False,
        fontsize=font, title_fontsize=font,
    )


def add_block_separators(ax, blocks, color: str = SEPARATOR_COLOR,
                         lw: float = SEPARATOR_LW, alpha: float = SEPARATOR_ALPHA):
    """Draw vertical dividers between consecutive columns whose block label
    differs. ``blocks`` is a length-n_cols sequence of block labels (one per
    gene column). Pass ``None`` or an all-same sequence to draw nothing.

    Typical block sources:

    * **Per-subtype panels** where each gene is owned by one subtype: pull
      from a per-subtype ``marker_set`` table — `marker_set.drop_duplicates(
      "gene", keep="first")[["gene", "<owner_col>"]]`. The dotplot's own
      data is often a full grid (every row × every gene) where first-
      occurrence collapses to the first row's label, so prefer the
      ``marker_set`` table.

    * **Shared-marker panels** (e.g. identity / QC panels where every row
      shares the same marker list): hand-code a ``gene -> functional_block``
      dict (identity / state / off-lineage / etc.) and feed that.
    """
    if blocks is None:
        return
    previous = None
    for i, block in enumerate(blocks):
        if previous is not None and block != previous:
            ax.axvline(i - 0.5, color=color, linewidth=lw, alpha=alpha)
        previous = block


def add_row_color_dots(ax, colors, *, x: float = -0.022,
                       size: float = 60.0, edgecolor: str = "none"):
    """Colour dot per row just left of the y-axis (e.g. a cell-type colour
    key), aligned with the dot-grid rows.

    ``colors`` is a sequence of length n_rows (row 0 at top, matching
    ``pad_axes``). Use any matplotlib colour; pass ``"#bdbdbd"`` for unknown
    rows. Drawn in a blended transform (x in axes coords, y in data coords) so
    it tracks the rows regardless of axis limits and is never clipped.
    """
    n_rows = len(colors)
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    ax.scatter(
        [x] * n_rows, np.arange(n_rows),
        s=size, c=list(colors), marker="o",
        edgecolors=edgecolor, linewidths=0,
        clip_on=False, zorder=10, transform=trans,
    )


# ============================================================
# one-call convenience
# ============================================================
def _resolve_blocks(block_per_gene, gene_order) -> list | None:
    if block_per_gene is None:
        return None
    if isinstance(block_per_gene, Mapping):
        return [block_per_gene.get(g) for g in gene_order]
    if isinstance(block_per_gene, Sequence) and not isinstance(block_per_gene, str):
        seq = list(block_per_gene)
        if len(seq) != len(gene_order):
            raise ValueError(
                f"block_per_gene length {len(seq)} != len(gene_order) {len(gene_order)}"
            )
        return seq
    raise TypeError(
        "block_per_gene must be a dict (gene -> block) or a sequence aligned "
        f"with gene_order, got {type(block_per_gene).__name__}"
    )


def _grid_layout(n_genes, n_rows, *, font, has_row_colors, has_title,
                 max_row_label, max_gene_label, legend_loc, fig_size,
                 max_pct=100.0, size_scale=2.4, size_base=4.0):
    """按内容算 dot plot 布局(英寸预算 → figure 分数)。

    pitch 实算 = ``max(单行文字含行距, 最大点直径) + gap``,取文字与点都刚好不重叠
    的最紧值,随字体 / 数据最大点自适应(竖排基因名 / 水平细胞名的文字下限一致)。网格
    = ``pitch × 数量``,左 / 下边距按 cell-type / gene 标签长度,右(或下)留固定图例带。
    ``fig_size=None`` 时各英寸预算相加得最紧凑 figsize;给定 figsize 则只据此换算边距
    分数。返回 ``((fig_w, fig_h), {left, right, bottom, top})``,分数已夹紧到合法区间。
    """
    cw = _CHAR_W_IN * (font / DEFAULT_FONT_SIZE)
    text_pitch = font / 72.0 * _TEXT_PITCH_K
    dot_d = 2.0 * (float(dot_sizes(max_pct, scale=size_scale, base=size_base))
                   / np.pi) ** 0.5 / 72.0
    pitch = max(max(text_pitch, dot_d) + _GRID_GAP_IN, _PITCH_MIN_IN)
    left_in = (_AXIS_LABEL_IN + max_row_label * cw
               + (_ROW_DOT_IN if has_row_colors else 0.0) + _PAD_IN)
    bottom_in = _AXIS_LABEL_IN + max_gene_label * cw + _PAD_IN
    top_in = _TITLE_IN if has_title else _PAD_IN + 0.06
    if legend_loc == "right":
        right_in = _RIGHT_LEGEND_IN
    else:
        right_in = _AXIS_LABEL_IN + 0.45
        bottom_in += _BOTTOM_LEGEND_IN
    grid_w = pitch * max(n_genes, 1)
    grid_h = pitch * max(n_rows, 1)
    if fig_size is None:
        fig_size = (left_in + grid_w + right_in, top_in + grid_h + bottom_in)
    fig_w, fig_h = float(fig_size[0]), float(fig_size[1])
    left = min(max(left_in / fig_w, 0.04), 0.6)
    right = max(min(1.0 - right_in / fig_w, 0.97), left + 0.15)
    bottom = min(max(bottom_in / fig_h, 0.04), 0.6)
    top = max(min(1.0 - top_in / fig_h, 0.97), bottom + 0.15)
    return (fig_w, fig_h), {"left": left, "right": right, "bottom": bottom, "top": top}


def _draw_right_legends(fig, ax, sc, *, right_f, fig_w, cbar_label, pcts,
                        size_scale, size_base, font: int = DEFAULT_FONT_SIZE):
    """右侧两个图例:colorbar(柱) + size 点列,严格共用一条中线、纵向对齐;两组
    数字(colorbar 刻度 / pct)左缘对齐,名称都竖排放最右、对齐成列。

    水平位置按「主图右缘 ``right_f`` + 英寸偏移」算,随 figsize 自适应、不写死分数。
    size 点不走 ``fig.legend``(点在图例框内有不可控偏移、对不齐 colorbar),而是用
    figure 坐标(``transFigure`` + ``clip_on=False``)直接画在 ``col_x`` 上;colorbar
    柱中心也设在 ``col_x``,二者精确对齐。colorbar 刻度数字再用 ``tick_params(pad)``
    推到与 pct 标签同一 x,两组数字左对齐。
    """
    legend_pcts = list(pcts) if pcts else list(DOT_LEGEND_PCTS)
    sizes = [float(dot_sizes(p, scale=size_scale, base=size_base)) for p in legend_pcts]
    # 柱宽 = 最大点直径 → colorbar 柱与 size 点列等宽、右缘齐平,数字紧贴两者并左对齐
    r_fx = ((max(sizes) / np.pi) ** 0.5 / 72.0) / fig_w   # 最大点半径(fraction)
    bar_w = 2.0 * r_fx
    col_x = right_f + 0.10 / fig_w + bar_w / 2.0           # 柱中心 = 点列中心(距图 0.10in)
    bar_right = col_x + bar_w / 2.0
    cbar_h = min(0.40, 2.2 / fig.get_size_inches()[1])
    cbar_y0 = 0.90 - cbar_h
    cax = fig.add_axes([col_x - bar_w / 2.0, cbar_y0, bar_w, cbar_h])
    cbar = fig.colorbar(sc, cax=cax)

    y_top, y_bot = 0.34, 0.14
    ys = (np.linspace(y_top, y_bot, len(legend_pcts)) if len(legend_pcts) > 1
          else [(y_top + y_bot) / 2.0])
    # 数字紧贴 柱 / 点 右缘(极小空隙),colorbar 刻度用 pad 推到同一 x、与 pct 左对齐
    label_x = col_x + r_fx + 0.004
    cbar.ax.tick_params(labelsize=font,
                        pad=max((label_x - bar_right) * fig_w * 72.0, 1.0))
    labels = []
    for y, p, s in zip(ys, legend_pcts, sizes):
        ax.scatter([col_x], [y], s=s, color="#888888", edgecolors=DOT_EDGE,
                   linewidths=DOT_EDGE_LW, transform=fig.transFigure,
                   clip_on=False, zorder=5)
        labels.append(ax.text(label_x, y, f"{p}%", transform=fig.transFigure,
                              va="center", ha="left", fontsize=font))

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    cbar_right = inv.transform(cbar.ax.get_tightbbox(renderer).get_points())[1][0]
    size_right = max(inv.transform(t.get_window_extent(renderer).get_points())[1][0]
                     for t in labels)
    title_x = max(cbar_right, size_right) + 0.005
    fig.text(title_x, cbar_y0 + cbar_h / 2.0, cbar_label, rotation=90,
             va="center", ha="left", fontsize=font)
    fig.text(title_x, (y_top + y_bot) / 2.0, "Pct expressed", rotation=90,
             va="center", ha="left", fontsize=font)
    return cbar


def marker_dotplot(
    data: pd.DataFrame,
    row_order: Sequence[str],
    gene_order: Sequence[str],
    *,
    row_col: str = "subtype",
    gene_col: str = "gene",
    avg_col: str = "avg_expr",
    pct_col: str = "pct_expr",
    title: str | None = None,
    xlabel: str = "Marker gene",
    ylabel: str | None = None,
    row_labels: Sequence[str] | None = None,
    fig_size: tuple[float, float] | None = None,
    font: int = DEFAULT_FONT_SIZE,
    block_per_gene=None,
    scale_per_gene: bool = True,
    scale: str | None = None,
    cmap=GREY_RED_CMAP,
    vmin: float | None = None,
    vmax: float | None = None,
    row_colors=None,
    row_color_size: float = 60.0,
    pad: float = 0.6,
    cbar_label: str = COLORBAR_LABEL,
    size_pcts: Sequence[int] | None = None,
    size_scale: float = 2.4,
    size_base: float = 4.0,
    legend_loc: str = "right",
    ax=None,
):
    """One-call marker dot plot in the unified template style.

    Parameters
    ----------
    data : DataFrame
        Long-form table with at least the row label, gene, mean expression
        and percent-expressed columns named by ``row_col`` / ``gene_col`` /
        ``avg_col`` / ``pct_col``.
    row_order, gene_order
        Row (y) and column (x) ordering.
    row_col, gene_col, avg_col, pct_col
        Column names in ``data``.
    title, xlabel, ylabel
        Axis labels and figure title (``ylabel`` defaults to ``row_col``).
    row_labels
        Optional pretty labels to use on the y-axis instead of ``row_order``.
    fig_size
        ``(width, height)`` inches. ``None`` (default) → content-driven inch
        budget: grid sized by per-gene / per-row pitch, margins by label length,
        plus a fixed right legend band. Margins and the right legend also adapt
        to any explicit ``fig_size`` you pass (no hard-coded figure fractions).
    block_per_gene
        Dict ``gene -> block`` or sequence aligned with ``gene_order`` for
        vertical separators. ``None`` (default) draws no separators.
    scale_per_gene
        Legacy switch kept for back-compat. Used only when ``scale`` is None:
        True -> ``scale="gene"``, False -> ``scale="none"``.
    scale
        Colour scaling mode. ``"gene"`` = per-gene (column) min-max to [0, 1]
        (default, the viridis "gene-scaled mean" style). ``"row"`` = per-row
        min-max to [0, 1] (matches scanpy ``standard_scale="group"``).
        ``"none"`` = colour by raw ``avg_col``; colour-bar spans the data range
        (or ``vmin``/``vmax`` if given). When None, falls back to
        ``scale_per_gene``.
    cmap
        Matplotlib colormap name or object. Default grey->red (scseq style,
        ``GREY_RED_CMAP``); pass ``"viridis"`` (``DOT_CMAP``) or any cmap.
    vmin, vmax
        Override the colour limits. Default: 0/1 for ``"gene"``/``"row"``,
        data min/max for ``"none"``.
    row_colors
        Optional per-row colour key drawn just left of the y-axis (cell-type
        colours). Dict ``row -> colour`` (keyed by ``row_order``) or a sequence
        aligned with ``row_order``. ``None`` (default) draws nothing.
    row_color_size
        Marker area (pt^2) for the ``row_colors`` dots.
    pad
        Axis padding (cells of margin around the grid). Default 0.6.
    cbar_label
        Colour-bar label. Default ``"Gene-scaled mean"``.
    ax
        Existing axes to draw into. If ``None``, a new figure is created.

    Returns
    -------
    fig, ax, sc
        The matplotlib figure, axes, and the scatter ``PathCollection`` (for
        attaching a custom colour-bar if needed).
    """
    df = data.loc[data[row_col].isin(row_order) & data[gene_col].isin(gene_order)].copy()
    df[row_col] = pd.Categorical(df[row_col], categories=list(row_order), ordered=True)
    df[gene_col] = pd.Categorical(df[gene_col], categories=list(gene_order), ordered=True)
    df = df.sort_values([row_col, gene_col])

    if scale is None:
        scale = "gene" if scale_per_gene else "none"

    def _minmax(s: pd.Series) -> pd.Series:
        arr = s.to_numpy(dtype=float)
        lo, hi = np.nanmin(arr), np.nanmax(arr)
        if hi > lo:
            return (s - lo) / (hi - lo)
        return pd.Series(np.full(len(s), 0.5), index=s.index)

    if scale == "gene":
        df["__color"] = df.groupby(gene_col, observed=True)[avg_col].transform(_minmax)
        cvmin, cvmax = 0.0, 1.0
    elif scale == "row":
        df["__color"] = df.groupby(row_col, observed=True)[avg_col].transform(_minmax)
        cvmin, cvmax = 0.0, 1.0
    elif scale == "none":
        df["__color"] = df[avg_col].astype(float)
        col = df["__color"].to_numpy(dtype=float)
        cvmin = float(np.nanmin(col)) if col.size else 0.0
        cvmax = float(np.nanmax(col)) if col.size else 1.0
    else:
        raise ValueError("scale must be one of {'gene', 'row', 'none'} or None")
    if vmin is not None:
        cvmin = float(vmin)
    if vmax is not None:
        cvmax = float(vmax)

    x_map = {g: i for i, g in enumerate(gene_order)}
    y_map = {r: i for i, r in enumerate(row_order)}
    x = df[gene_col].astype(str).map(x_map).to_numpy()
    y = df[row_col].astype(str).map(y_map).to_numpy()

    row_disp = list(row_labels) if row_labels is not None else list(row_order)
    max_row_label = max((len(str(r)) for r in row_disp), default=1)
    max_gene_label = max((len(str(g)) for g in gene_order), default=1)
    _pct = df[pct_col].to_numpy(dtype=float)
    max_pct = float(np.clip(_pct.max() if _pct.size else 100.0, 1.0, 100.0))
    layout_size = fig_size if ax is None else tuple(ax.figure.get_size_inches())
    fig_size, frac = _grid_layout(
        len(gene_order), len(row_order), font=font,
        has_row_colors=row_colors is not None, has_title=bool(title),
        max_row_label=max_row_label, max_gene_label=max_gene_label,
        legend_loc=legend_loc, fig_size=layout_size,
        max_pct=max_pct, size_scale=size_scale, size_base=size_base,
    )
    if ax is None:
        fig, ax = plt.subplots(figsize=fig_size)
    else:
        fig = ax.figure
    fig_w = fig.get_size_inches()[0]

    sc = ax.scatter(
        x, y,
        s=dot_sizes(df[pct_col].to_numpy(dtype=float),
                    scale=size_scale, base=size_base),
        c=df["__color"].to_numpy(dtype=float),
        cmap=cmap, vmin=cvmin, vmax=cvmax,
        linewidths=DOT_EDGE_LW, edgecolors=DOT_EDGE,
    )
    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(list(gene_order), rotation=90, ha="center", va="top", fontsize=font)
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(list(row_labels) if row_labels is not None else list(row_order),
                       fontsize=font)
    if title:
        ax.set_title(title, fontsize=font, pad=4)
    ax.set_xlabel(xlabel, fontsize=font)
    ax.set_ylabel(ylabel if ylabel is not None else row_col, fontsize=font)
    ax.grid(color=GRID_COLOR, linewidth=GRID_LW, alpha=GRID_ALPHA)
    ax.tick_params(axis="both", length=0, labelsize=font)
    for spine in ax.spines.values():
        spine.set_visible(False)

    pad_axes(ax, len(gene_order), len(row_order), pad=pad)
    blocks = _resolve_blocks(block_per_gene, gene_order)
    if blocks is not None:
        add_block_separators(ax, blocks)

    if row_colors is not None:
        if isinstance(row_colors, Mapping):
            row_color_seq = [row_colors.get(r, "#bdbdbd") for r in row_order]
        else:
            row_color_seq = list(row_colors)
            if len(row_color_seq) != len(row_order):
                raise ValueError(
                    f"row_colors length {len(row_color_seq)} != "
                    f"len(row_order) {len(row_order)}"
                )
        # 彩点贴近 y 轴(x=-0.010),y 标签 pad 让其右缘与彩点留一点点空隙、不挨着
        ax.tick_params(axis="y", pad=float(np.clip(8.0 + row_color_size / 15.0, 10.0, 18.0)))
        add_row_color_dots(ax, row_color_seq, size=row_color_size, x=-0.010)

    if legend_loc == "right":
        # colorbar 柱 + size 点列共用一条中线、纵向对齐;名称竖排放最右对齐成列
        fig.subplots_adjust(**frac)
        _draw_right_legends(fig, ax, sc, right_f=frac["right"], fig_w=fig_w,
                            cbar_label=cbar_label, pcts=size_pcts,
                            size_scale=size_scale, size_base=size_base, font=font)
    else:  # "bottom": colorbar 右侧 + size 图例 底部横排(旧版式)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.022, pad=0.025)
        style_colorbar(cbar, font=font, label=cbar_label)
        add_size_legend(ax, font=font, pcts=size_pcts,
                        size_scale=size_scale, size_base=size_base)
        fig.subplots_adjust(left=frac["left"], right=0.84,
                            top=frac["top"], bottom=frac["bottom"])
    return fig, ax, sc


__all__ = [
    # constants
    "DOT_CMAP", "GREY_RED_CMAP", "DOT_EDGE", "DOT_EDGE_LW", "DOT_LEGEND_PCTS", "COLORBAR_LABEL",
    "SEPARATOR_COLOR", "SEPARATOR_LW", "SEPARATOR_ALPHA",
    "GRID_COLOR", "GRID_LW", "GRID_ALPHA", "HIGHLIGHT_COLOR",
    # primitives
    "dot_sizes", "pad_axes", "style_colorbar", "add_size_legend",
    "add_block_separators", "add_row_color_dots",
    # one-call
    "marker_dotplot",
]
