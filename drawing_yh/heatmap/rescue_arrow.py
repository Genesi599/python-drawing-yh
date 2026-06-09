"""Rescue-style directional arrow heatmap templates."""
from __future__ import annotations

from typing import Mapping, Sequence

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np
import pandas as pd

DEFAULT_FONT_SIZE = 8
DOUBLE_COL_IN = 6.89

DEFAULT_ARROW_COLORS = {
    "aging_up": "#D9822B",
    "aging_down": "#3E7CB1",
    "rescue_up": "#CC6677",
    "rescue_down": "#1B7F79",
}
DEFAULT_BACKGROUND_COLORS = {
    "aging_up": "#F7E3DD",
    "aging_down": "#DDEBEC",
    "present": "#F1F3F5",
    "missing": "#E6E6E6",
}


def _ordered(values: Sequence, order: Sequence | None) -> list:
    if order is not None:
        return list(order)
    return list(pd.Index(values).drop_duplicates())


def _scale_length(value: float, scale: float, *, min_len: float = 0.05, max_len: float = 0.34) -> float:
    if not np.isfinite(value) or value == 0:
        return 0.0
    length = abs(value) / max(scale, 1e-12) * max_len
    return float(np.clip(length, min_len, max_len))


def _alpha_from_padj(padj: float | None, max_score: float) -> float:
    if padj is None or not np.isfinite(padj) or padj <= 0:
        return 0.95
    score = min(-np.log10(max(padj, 1e-300)), max_score)
    return float(0.35 + 0.60 * score / max(max_score, 1e-12))


def _draw_vertical_arrow(ax, x: float, y0: float, delta: float, color: str, alpha: float) -> None:
    if delta == 0:
        return
    arrow = FancyArrowPatch(
        (x, y0),
        (x, y0 + delta),
        arrowstyle="-|>",
        mutation_scale=7,
        linewidth=0.8,
        color=color,
        alpha=alpha,
        shrinkA=0,
        shrinkB=0,
    )
    ax.add_patch(arrow)


def plot_rescue_arrow_heatmap(
    data: pd.DataFrame,
    *,
    row: str,
    column: str,
    aging_lfc: str,
    rescue_lfc: str,
    aging_padj: str | None = None,
    rescue_padj: str | None = None,
    count: str | None = None,
    row_order: Sequence | None = None,
    column_order: Sequence | None = None,
    row_label_map: Mapping | None = None,
    column_label_map: Mapping | None = None,
    title: str | None = None,
    xlabel: str = "Context",
    ylabel: str = "Gene",
    count_label: str | None = None,
    arrow_colors: Mapping[str, str] | None = None,
    background_colors: Mapping[str, str] | None = None,
    lfc_scale: float | None = None,
    padj_score_cap: float = 12.0,
    show_legend: bool = True,
    show_count: bool = True,
    xtick_rotation: float = 45,
    figsize: tuple[float, float] | None = None,
    width: float | None = DOUBLE_COL_IN,
    cell_width: float = 0.44,
    cell_height: float = 0.24,
):
    """Plot a compact rescue heatmap with paired aging and rescue arrows.

    The template expects one row per feature/context pair. Inside each tile,
    the left arrow encodes the aging effect and the right arrow encodes the
    rescue/treatment effect. Arrow length is proportional to absolute log2FC;
    alpha is proportional to ``-log10(adjusted P)`` when p-value columns are
    provided. Optional cell counts are drawn in the upper-left corner.
    """
    if data.empty:
        raise ValueError("data is empty")

    colors = dict(DEFAULT_ARROW_COLORS)
    if arrow_colors:
        colors.update(arrow_colors)
    bg = dict(DEFAULT_BACKGROUND_COLORS)
    if background_colors:
        bg.update(background_colors)

    df = data.copy()
    rows = _ordered(df[row], row_order)
    cols = _ordered(df[column], column_order)
    n_rows, n_cols = len(rows), len(cols)
    row_pos = {v: i for i, v in enumerate(rows)}
    col_pos = {v: j for j, v in enumerate(cols)}

    if lfc_scale is None:
        vals = pd.concat([df[aging_lfc].abs(), df[rescue_lfc].abs()], ignore_index=True)
        finite = vals[np.isfinite(vals)]
        lfc_scale = float(finite.quantile(0.95)) if len(finite) else 1.0
        if not np.isfinite(lfc_scale) or lfc_scale == 0:
            lfc_scale = 1.0

    if figsize is None:
        fig_w = width if width is not None else float(np.clip(n_cols * cell_width + 2.4, 3.35, 8.5))
        fig_h = float(np.clip(n_rows * cell_height + 1.35, 2.4, 8.5))
        figsize = (fig_w, fig_h)

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)

    present = set(zip(df[row], df[column]))
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            tile_bg = bg["missing"]
            if (r, c) in present:
                sub = df[(df[row] == r) & (df[column] == c)]
                a = float(sub[aging_lfc].iloc[0])
                tile_bg = bg["aging_up"] if a > 0 else bg["aging_down"] if a < 0 else bg["present"]
            ax.add_patch(
                Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=tile_bg, edgecolor="white", linewidth=0.8)
            )

    for _, rec in df.iterrows():
        r = rec[row]
        c = rec[column]
        if r not in row_pos or c not in col_pos:
            continue
        i = row_pos[r]
        j = col_pos[c]
        center_y = i
        aging = float(rec[aging_lfc])
        rescue = float(rec[rescue_lfc])
        aging_len = _scale_length(aging, lfc_scale)
        rescue_len = _scale_length(rescue, lfc_scale)
        aging_delta = -np.sign(aging) * aging_len
        rescue_delta = -np.sign(rescue) * rescue_len

        aging_color = colors["aging_up"] if aging > 0 else colors["aging_down"]
        rescue_color = colors["rescue_up"] if rescue > 0 else colors["rescue_down"]
        aging_alpha = _alpha_from_padj(float(rec[aging_padj]) if aging_padj else None, padj_score_cap)
        rescue_alpha = _alpha_from_padj(float(rec[rescue_padj]) if rescue_padj else None, padj_score_cap)

        x_aging = j - 0.13
        x_rescue = j + 0.13
        _draw_vertical_arrow(ax, x_aging, center_y, aging_delta, aging_color, aging_alpha)
        _draw_vertical_arrow(ax, x_rescue, center_y + aging_delta, rescue_delta, rescue_color, rescue_alpha)

        if show_count and count and pd.notna(rec[count]):
            label = str(int(rec[count])) if float(rec[count]).is_integer() else str(rec[count])
            ax.text(j - 0.41, i - 0.34, label, ha="left", va="top", fontsize=5.8, color="#333333")

    row_labels = [str(row_label_map.get(r, r) if row_label_map else r) for r in rows]
    col_labels = [str(column_label_map.get(c, c) if column_label_map else c) for c in cols]
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(col_labels, rotation=xtick_rotation, ha="right", fontsize=DEFAULT_FONT_SIZE)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels, fontsize=DEFAULT_FONT_SIZE, fontstyle="italic")
    ax.set_xlabel(xlabel, fontsize=DEFAULT_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=DEFAULT_FONT_SIZE)
    if title:
        ax.set_title(title, fontsize=DEFAULT_FONT_SIZE, pad=5)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if count_label:
        ax.text(
            n_cols - 0.5,
            -0.76,
            count_label,
            ha="right",
            va="center",
            fontsize=DEFAULT_FONT_SIZE - 1,
            color="#333333",
            clip_on=False,
        )

    if show_legend:
        handles = [
            Line2D([0], [0], marker=r"$\uparrow$", linestyle="None", color=colors["aging_up"], label="Aging up"),
            Line2D([0], [0], marker=r"$\downarrow$", linestyle="None", color=colors["aging_down"], label="Aging down"),
            Line2D([0], [0], marker=r"$\uparrow$", linestyle="None", color=colors["rescue_up"], label="Rescue up"),
            Line2D([0], [0], marker=r"$\downarrow$", linestyle="None", color=colors["rescue_down"], label="Rescue down"),
        ]
        ax.legend(
            handles=handles,
            frameon=False,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.02),
            ncol=4,
            columnspacing=0.9,
            handletextpad=0.3,
            borderaxespad=0,
        )
    fig.tight_layout(pad=0.5)
    return fig, ax


__all__ = ["plot_rescue_arrow_heatmap"]
