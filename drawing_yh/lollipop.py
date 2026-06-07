"""
lollipop.py — Horizontal lollipop (stem-and-dot) plot for drawing_yh.

Typical use: ranked TF effect sizes, gene rankings, any signed score list.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from drawing_yh import DEFAULT_FONT_SIZE


def lollipop_plot(
    values: pd.Series,
    *,
    sort_by_abs: bool = True,
    top_n: int | None = None,
    color_positive: str = "#d6604d",
    color_negative: str = "#4393c3",
    stem_color: str = "#aaaaaa",
    stem_lw: float = 0.8,
    marker_size: float = 5,
    xlabel: str = "Mean effect (old − young)",
    ylabel: str = "",
    title: str | None = None,
    figsize: tuple[float, float] | None = None,
    show_legend: bool = True,
    legend_labels: tuple[str, str] = ("Up in old", "Down in old"),
) -> tuple:
    """Horizontal lollipop plot for ranked signed effect sizes.

    Parameters
    ----------
    values : pd.Series
        Index = item labels (e.g. TF names), values = signed numeric effect.
        Items are shown as horizontal stems ending in colored dots.
    sort_by_abs : bool
        If True, sort by |value| ascending so the largest bars are at top.
    top_n : int, optional
        Keep only the top-N items by |value| before plotting.
    """
    s = values.dropna().copy()

    if top_n is not None:
        s = s.reindex(s.abs().sort_values(ascending=False).index[:top_n])

    if sort_by_abs:
        # ascending so strongest ends up at the top of the chart
        s = s.reindex(s.abs().sort_values(ascending=True).index)

    n = len(s)
    if n == 0:
        raise ValueError("lollipop_plot: no data to plot after filtering.")

    if figsize is None:
        figsize = (4.0, max(1.5, 0.22 * n + 1.2))

    fig, ax = plt.subplots(figsize=figsize)
    y = np.arange(n)
    colors = [color_positive if v >= 0 else color_negative for v in s.values]

    ax.hlines(y, 0, s.values, colors=stem_color, linewidths=stem_lw)
    ax.scatter(s.values, y, color=colors, s=marker_size ** 2,
               zorder=3, edgecolors="none")
    ax.axvline(0, color="#606060", lw=0.6)

    ax.set_yticks(y)
    ax.set_yticklabels(s.index.tolist(), fontsize=DEFAULT_FONT_SIZE)
    ax.set_xlabel(xlabel, fontsize=DEFAULT_FONT_SIZE)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=DEFAULT_FONT_SIZE)
    if title:
        ax.set_title(title, fontsize=DEFAULT_FONT_SIZE)

    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=DEFAULT_FONT_SIZE)

    if show_legend:
        handles = [
            Patch(facecolor=color_positive, label=legend_labels[0], edgecolor="none"),
            Patch(facecolor=color_negative, label=legend_labels[1], edgecolor="none"),
        ]
        ax.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.12),
            ncol=2,
            frameon=False,
            fontsize=DEFAULT_FONT_SIZE - 1,
        )

    fig.tight_layout()
    return fig, ax


__all__ = ["lollipop_plot"]
