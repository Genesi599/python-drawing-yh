# -*- coding: utf-8 -*-
"""
Heatmap with per-row horizontal bars — unified template (pure matplotlib).

Left: a diverging value heatmap (typically per-group mean z-score), genes on x
(optionally split into per-group column blocks by vertical dividers), groups on
y with an optional left colour strip. Right: a few horizontal bars per row
(e.g. that group's top enriched terms), coloured by the group colour, with the
term text annotated at the bar end. The two panels share the y axis so each
bar row lines up with its heatmap row.

This is the **drawing layer only** — it takes an already-computed matrix and an
already-computed per-row bar table (``row -> [(label, value), ...]``). The
single-cell prep (read h5ad, pick markers, mean + z-score, run enrichment)
lives in ``single_cell-yh`` and calls ``heatmap_with_row_bars`` here, so
``drawing_yh`` stays free of scanpy / gseapy.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import transforms
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

from drawing_yh import DEFAULT_FONT_SIZE, DOUBLE_COL_IN


DEFAULT_ROW_COLOR = "#6baed6"
DIVERGING_COLORS = ("#2166ac", "#ffffff", "#b2182b")


def wrap_text(text: str, max_chars: int = 40, max_lines: int = 0) -> str:
    """Word-wrap ``text`` to lines of at most ``max_chars`` (no hyphenation).

    max_lines > 0 时截到该行数、末行尾加省略号(2026-09-20: 密集行距图 term 4 行撞行,
    截断保 8pt 行高, 不缩字号)。"""
    if not text:
        return ""
    words = str(text).split()
    if not words:
        return str(text)
    lines, cur = [], words[0]
    for w in words[1:]:
        if len(cur) + 1 + len(w) <= max_chars:
            cur += " " + w
        else:
            lines.append(cur)
            cur = w
    lines.append(cur)
    if max_lines and len(lines) > max_lines:
        cut = lines[:max_lines]
        cut[-1] = cut[-1][:max(1, max_chars - 1)].rstrip() + "..."
        return "\n".join(cut)
    return "\n".join(lines)


def _resolve_cmap(cmap):
    if isinstance(cmap, str):
        return plt.get_cmap(cmap)
    if hasattr(cmap, "__call__"):  # already a Colormap
        return cmap
    return LinearSegmentedColormap.from_list("diverging", list(cmap))


def _resolve_row_colors(row_colors, row_labels, default):
    if row_colors is None:
        return None
    if isinstance(row_colors, Mapping):
        return [row_colors.get(lbl, default) for lbl in row_labels]
    seq = list(row_colors)
    if len(seq) != len(row_labels):
        raise ValueError(
            f"row_colors length {len(seq)} != len(row_labels) {len(row_labels)}"
        )
    return seq


def _resolve_row_bars(row_bars, row_labels, max_bars_per_row):
    if row_bars is None:
        return None
    if isinstance(row_bars, Mapping):
        return [list(row_bars.get(lbl, []))[:max_bars_per_row] for lbl in row_labels]
    seq = list(row_bars)
    out = []
    for i in range(len(row_labels)):
        out.append(list(seq[i])[:max_bars_per_row] if i < len(seq) else [])
    return out


def heatmap_with_row_bars(
    Z,
    row_labels: Sequence[str],
    *,
    row_bars=None,
    block_sizes: Sequence[int] | None = None,
    row_colors=None,
    gene_labels: Sequence[str] | None = None,
    z_clip: float = 2.0,
    cmap="RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    max_bars_per_row: int = 2,
    bar_xlabel: str = r"$-\log_{10}$(p-value)",
    cbar_label: str = "Expression Z-score",
    left_right_width_ratio: tuple[float, float] = (2.5, 1.0),
    row_color_strip: bool = True,
    strip_width: float = 0.018,
    wrap_chars: int = 40,
    font: int = DEFAULT_FONT_SIZE,
    figsize: tuple[float, float] | None = None,
    fig_height: float | None = None,
    title: str | None = None,
    default_color: str = DEFAULT_ROW_COLOR,
):
    """Diverging heatmap + per-row horizontal bars.

    Parameters
    ----------
    Z : (G, P) array-like
        Values to colour (rows = groups, columns = genes). Typically per-group
        mean z-score.
    row_labels : length-G sequence
        Y labels (already in display form).
    row_bars : dict | sequence | None
        Per-row bars. ``{row_label: [(text, value), ...]}`` or a sequence
        aligned with rows. Each row keeps the first ``max_bars_per_row``. None
        draws no right panel (heatmap only).
    block_sizes : sequence of int | None
        Column block sizes (e.g. each group's marker count). Draws vertical
        dividers between consecutive blocks. None / single block = no dividers.
    row_colors : dict | sequence | None
        Per-row colour, used for the left colour strip and the bar colour.
        ``{row_label: colour}`` or a sequence aligned with rows.
    gene_labels : sequence | None
        Column (x) labels. None (default) hides the x axis (too many genes).
    z_clip, vmin, vmax
        Colour limits. Default symmetric ``[-z_clip, z_clip]``; override with
        explicit ``vmin`` / ``vmax``.
    cmap
        Diverging colormap name / object / colour list (default ``"RdBu_r"``).
    max_bars_per_row, bar_xlabel
        Right-panel bar count per row and its x-axis label.
    cbar_label
        Heatmap colour-bar label.
    left_right_width_ratio
        Width ratio of (heatmap, bars) panels.
    row_color_strip, strip_width
        Draw / size the left per-row colour strip (needs ``row_colors``).

    Returns
    -------
    (fig, (ax_heat, ax_bars))
        ``ax_bars`` is None when ``row_bars`` is None.
    """
    Z = np.asarray(Z, dtype=float)
    if Z.ndim != 2:
        raise ValueError(f"Z must be 2-D (G, P), got shape {Z.shape}")
    G, P = Z.shape
    row_labels = list(row_labels)
    if len(row_labels) != G:
        raise ValueError(f"row_labels length {len(row_labels)} != Z rows {G}")

    cmap_obj = _resolve_cmap(cmap)
    colors = _resolve_row_colors(row_colors, row_labels, default_color)
    bars_by_row = _resolve_row_bars(row_bars, row_labels, max_bars_per_row)
    draw_bars = bars_by_row is not None
    hi = z_clip if vmax is None else vmax
    lo = -z_clip if vmin is None else vmin

    left_w, right_w = left_right_width_ratio
    if figsize is None:
        # 紧凑 default(checklist §5):无 bars 行高 0.18 in;有 bars 时 term wrap 多行
        # 撑行高需 0.28 in 才不撞(实测 wrap_chars≈28 → 通常 2 行 term)
        if fig_height is not None:
            fig_h = float(fig_height)
        elif draw_bars:
            # 行高 0.32 in/row(为 row_bars 多行 term 留 0.22 in 文字 + ~0.1 in 余量)
            fig_h = float(np.clip(0.32 * G + 1.2, 2.4, 11.0))
        else:
            fig_h = float(np.clip(0.18 * G + 0.9, 1.8, 8.0))
        if draw_bars:
            fig_w = min(DOUBLE_COL_IN * 1.4, fig_h * (left_w + right_w))
        else:
            fig_w = min(DOUBLE_COL_IN, fig_h * left_w)
        figsize = (fig_w, fig_h)

    if draw_bars:
        fig, (ax_heat, ax_bars) = plt.subplots(
            1, 2, figsize=figsize,
            gridspec_kw={"width_ratios": [left_w, right_w]},
            sharey=True, constrained_layout=True,
        )
    else:
        fig, ax_heat = plt.subplots(figsize=figsize, constrained_layout=True)
        ax_bars = None

    # ---- left: heatmap ----
    im = ax_heat.imshow(Z, aspect="auto", interpolation="nearest",
                        cmap=cmap_obj, vmin=lo, vmax=hi, origin="upper")
    ax_heat.set_yticks(np.arange(G))
    ax_heat.set_yticklabels(row_labels, fontsize=font)
    ax_heat.set_ylabel("")
    if gene_labels is not None:
        ax_heat.set_xticks(np.arange(P))
        ax_heat.set_xticklabels(list(gene_labels), rotation=90, ha="center",
                                fontsize=font, fontstyle="italic")
    else:
        ax_heat.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_heat.tick_params(axis="y", length=0, pad=4)
    ax_heat.set_xlabel("")

    if block_sizes:
        for x in np.cumsum(list(block_sizes))[:-1]:
            ax_heat.axvline(x - 0.5, color="k", lw=0.5, alpha=0.5)

    # 迷你色条(2026-09-20 杨航): 横窄条贴右下角底缘, 无 label(caption 里说明), 刻度 min/0/max。
    # 位置按英寸换算; 右下角通常是 -logp 刻度区, 无 term(文字在条区上方), 无冲突。
    fw, fh = fig.get_size_inches()
    # 锚在热图区下方右侧(与热图同列, 不进 term 区; 底缘留出 x 轴空间)
    fig.canvas.draw()
    _hb = ax_heat.get_window_extent(fig.canvas.get_renderer())
    left_in = _hb.x1 / fig.dpi - 0.9
    if left_in < _hb.x0 / fig.dpi:
        left_in = _hb.x0 / fig.dpi
    cbax = fig.add_axes([left_in / fw, 0.075, min(0.9, _hb.width / fig.dpi) / fw, 0.05 / fh],
                        label="mini_cbar", in_layout=False)   # 不参与 constrained layout
    cbar = fig.colorbar(im, cax=cbax, orientation="horizontal")
    vmin_, vmax_ = im.get_clim()
    ticks_ = [vmin_, 0.0, vmax_] if vmin_ < 0 < vmax_ else [vmin_, vmax_]
    cbar.set_ticks(ticks_)
    cbar.ax.tick_params(labelsize=font, pad=1, top=True, labeltop=True,
                        bottom=False, labelbottom=False)   # 刻度在条上方, 不垂出底缘
    cbar.outline.set_linewidth(0.5)

    if row_color_strip and colors is not None:
        trans = transforms.blended_transform_factory(ax_heat.transAxes, ax_heat.transData)
        for gi, col in enumerate(colors):
            ax_heat.add_patch(Rectangle(
                (-strip_width, gi - 0.5), strip_width, 1.0,
                transform=trans, facecolor=col, edgecolor="none",
                linewidth=0, zorder=3, clip_on=False,
            ))
        # FIX: y tick label 默认 pad=4 像素易被 strip(axes 左外侧 strip_width × axes 宽度)遮挡。
        # canvas.draw() 后实测 axes 宽度,把 pad 推到 strip 之外 +4 像素余量(checklist §10)。
        fig.canvas.draw()
        ax_w_px = ax_heat.get_window_extent().width
        ax_heat.tick_params(axis="y", length=0, pad=int(strip_width * ax_w_px) + 4)

    # ---- right: per-row bars ----
    if draw_bars:
        all_vals = [v for row in bars_by_row for (_, v) in row]
        xmax = float(max(all_vals)) if all_vals else 1.0
        offsets_two = (-0.2, 0.2)
        trans_bar = transforms.blended_transform_factory(ax_bars.transAxes, ax_bars.transData)
        for gi, row in enumerate(bars_by_row):
            if not row:
                continue
            offs = offsets_two if len(row) >= 2 else (0.0,)
            for j, (term, val) in enumerate(row):
                y = gi + offs[j] if j < len(offs) else gi
                col = colors[gi] if colors is not None else default_color
                ax_bars.barh(y, val, height=0.36, color=col, alpha=0.9, edgecolor="none")
                wrapped = wrap_text(term, max_chars=wrap_chars)   # max_chars>=56 时自然单行不换行不截断
                # 2026-09-20: 文字改放条末端右侧数据坐标(不出 axes, 过 stray 门禁);
                # xlim 预留 55% 右侧给 term 文字(原 annotate 越 axes 被 audit 报 stray)
                ax_bars.text(
                    xmax * 1.06, y, wrapped, ha="left", va="center",   # 统一起点=左对齐(2026-09-20 杨航)
                    fontsize=font, color="black", clip_on=False, linespacing=1, zorder=10,
                )
        bar_x1 = xmax * 1.05  # 条形实际数据范围右界(2026-09-20 杨航: 轴只在柱范围)
        ax_bars.set_xlabel(bar_xlabel, fontsize=font)
        ax_bars.grid(axis="x", linestyle=":", alpha=0.3)
        ax_bars.tick_params(axis="y", left=False, labelleft=False)
        ax_bars.tick_params(axis="x", labelsize=font)
        # 刻度减疏(2026-09-20 杨航: 间距不用那么大): 最多 4 个主刻度
        ax_bars.xaxis.set_major_locator(mticker.MaxNLocator(4))
        ax_bars.set_ylim(ax_heat.get_ylim())
        # 单行 term 完整显示: 迭代量 renderer 实测文字 bbox, 扩 xlim 直至最长 term 右缘入界
        # (固定系数/数据换算在不同 xmax 下都不稳, 2026-09-20 三轮实测后改实测迭代)。
        for _ in range(12):
            fig.canvas.draw()
            rend = fig.canvas.get_renderer()
            ab = ax_bars.get_window_extent(rend)
            txts = [t for t in ax_bars.texts if t.get_text().strip()]
            if not txts:
                break
            right_px = max(t.get_window_extent(rend).x1 for t in txts)
            if right_px <= ab.x1 + 0.5:
                break
            x0, x1 = ax_bars.get_xlim()
            need_frac = (right_px - ab.x1) / max(ab.width, 1.0)
            ax_bars.set_xlim(x0, x1 + need_frac * (x1 - x0) * 1.15)
        # 轴线/刻度只覆盖条形范围(2026-09-20 杨航两轮修正):
        # 根因: MaxNLocator 在 xlim 扩张前设, 扩张后刻度按全 xlim(含 term 文字区)重新生成 ->
        # 刻度被 term 宽度稀释、与柱子对不上。修: xlim 定稿后按 0..bar_x1 数据范围显式 set_xticks
        # (禁 locator 重算), 刻度线+标签只落在柱区内; 轴线用 axes-fraction 线段画到 bar_x1。
        ax_bars.spines["bottom"].set_visible(False)
        ticks = mticker.MaxNLocator(4).tick_values(0.0, bar_x1)
        ticks = [t for t in ticks if 0 <= t <= bar_x1]
        ax_bars.xaxis.set_major_locator(mticker.FixedLocator(ticks))
        ax_bars.set_xticks(ticks)
        x0f = (bar_x1 - ax_bars.get_xlim()[0]) / (ax_bars.get_xlim()[1] - ax_bars.get_xlim()[0])
        ax_bars.add_line(mlines.Line2D([0, x0f], [0, 0], transform=ax_bars.transAxes,
                                       color=ax_bars.spines["left"].get_edgecolor(), linewidth=0.8, clip_on=False))
        for sp in ("top", "right"):
            ax_bars.spines[sp].set_visible(False)


    if title:
        fig.suptitle(title, fontsize=font, fontweight="bold")

    return fig, (ax_heat, ax_bars)


__all__ = [
    "heatmap_with_row_bars",
    "wrap_text",
    "DEFAULT_ROW_COLOR",
    "DIVERGING_COLORS",
]
