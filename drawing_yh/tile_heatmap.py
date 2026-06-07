"""
tile_heatmap.py — TF-style square-tile heatmap for drawing_yh.

Abstracts the meninges 92_plot_tf_results.py style:
- pcolormesh + set_aspect('equal') → each cell is a square
- White grid lines (edgecolors)
- TwoSlopeNorm for diverging color (RdBu_r)
- Dynamic figsize based on tile_in (inches per cell)
- FDR star annotations in cell centers
- Optional horizontal divider lines (up/down group separator)
- Optional per-column x-tick colors (subtype coloring)
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, Normalize

from drawing_yh import DEFAULT_FONT_SIZE

# Default style constants (meninges-compatible)
_TILE_IN = 0.18       # inch per cell side (square)
_EDGE_COLOR = "white"
_EDGE_LW = 0.5
_W_MARGIN_IN = 2.0    # extra width for y-tick labels + colorbar
_CBAR_H = 0.5         # top headroom for colorbar / title
_DIVIDER_COLOR = "#000000"
_DIVIDER_LW = 0.8


def _compute_figsize(
    n_rows: int,
    n_cols: int,
    col_labels: Sequence[str],
    *,
    tile_in: float = _TILE_IN,
    w_margin_in: float = _W_MARGIN_IN,
    cbar_h: float = _CBAR_H,
) -> tuple[float, float]:
    axes_w = tile_in * n_cols
    axes_h = tile_in * n_rows
    # 45° rotated tick label height estimate
    max_label_len = max((len(c) for c in col_labels), default=1)
    xtick_h = max_label_len * DEFAULT_FONT_SIZE * 0.6 / 72 + 0.3
    w = axes_w + w_margin_in
    h = axes_h + xtick_h + cbar_h + 0.3
    return w, h


def tile_heatmap(
    matrix: pd.DataFrame,
    *,
    row_labels: Sequence[str] | None = None,
    col_labels: Sequence[str] | None = None,
    cmap: str = "RdBu_r",
    diverging: bool = True,
    vmax: float | None = None,
    tile_in: float = _TILE_IN,
    edge_color: str = _EDGE_COLOR,
    edge_lw: float = _EDGE_LW,
    pvalue_matrix: pd.DataFrame | None = None,
    xtick_colors: list[str] | None = None,
    divider_rows: list[int] | None = None,
    divider_color: str = _DIVIDER_COLOR,
    divider_lw: float = _DIVIDER_LW,
    title: str | None = None,
    xlabel: str = "",
    ylabel: str = "TF",
    cbar_label: str | None = None,
    w_margin_in: float = _W_MARGIN_IN,
    cbar_h: float = _CBAR_H,
    figsize: tuple[float, float] | None = None,
) -> tuple:
    """Render a TF-style square-tile heatmap.

    Parameters
    ----------
    matrix : pd.DataFrame
        Rows = TF names, columns = cell-type labels. Values = signed effect.
    row_labels, col_labels : sequence of str, optional
        Display labels for rows/columns. If omitted, ``matrix.index`` and
        ``matrix.columns`` are used.
    pvalue_matrix : pd.DataFrame, optional
        Same index/columns as *matrix*. Cells with padj < 0.05/0.01/0.001
        get *, **, *** annotations.
    xtick_colors : list[str], optional
        One hex color per column (for subtype-colored x-tick labels).
    divider_rows : list[int], optional
        Draw a horizontal line *after* each row index in this list
        (e.g. ``[n_up]`` to separate up/down TF groups).
    """
    n_rows, n_cols = matrix.shape
    if row_labels is None:
        row_labels = [str(r) for r in matrix.index]
    else:
        row_labels = [str(r) for r in row_labels]
    if col_labels is None:
        col_labels = [str(c) for c in matrix.columns]
    else:
        col_labels = [str(c) for c in col_labels]
    if len(row_labels) != n_rows:
        raise ValueError("row_labels length must match matrix rows.")
    if len(col_labels) != n_cols:
        raise ValueError("col_labels length must match matrix columns.")

    if figsize is None:
        figsize = _compute_figsize(
            n_rows, n_cols, col_labels,
            tile_in=tile_in, w_margin_in=w_margin_in, cbar_h=cbar_h,
        )

    fig, ax = plt.subplots(figsize=figsize)

    vals = matrix.values.astype(float)
    finite = vals[np.isfinite(vals)]
    if diverging:
        v = vmax
        if v is None:
            v = float(np.nanmax(np.abs(finite))) if finite.size else 0.5
            v = max(v, 0.5)
        norm = TwoSlopeNorm(vmin=-v, vcenter=0, vmax=v)
    else:
        norm = Normalize(
            vmin=float(np.nanmin(finite)) if finite.size else -1,
            vmax=float(np.nanmax(finite)) if finite.size else 1,
        )

    mesh = ax.pcolormesh(
        np.arange(n_cols + 1),
        np.arange(n_rows + 1),
        vals,
        cmap=cmap,
        norm=norm,
        edgecolors=edge_color,
        linewidth=edge_lw,
    )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)

    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(col_labels, rotation=45, ha="right")
    if xtick_colors is not None:
        for tick_label, color in zip(ax.get_xticklabels(), xtick_colors):
            tick_label.set_color(color)

    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(row_labels)

    # FDR star annotations
    if pvalue_matrix is not None:
        padj_aligned = pvalue_matrix.reindex(
            index=matrix.index, columns=matrix.columns
        ).values
        for i in range(n_rows):
            for j in range(n_cols):
                p = padj_aligned[i, j]
                if pd.isna(p):
                    continue
                mark = (
                    "***" if p < 0.001
                    else "**" if p < 0.01
                    else "*" if p < 0.05
                    else None
                )
                if mark is None:
                    continue
                v_cell = vals[i, j]
                v_ref = vmax if vmax is not None else 0.5
                txt_color = (
                    "#fff"
                    if (np.isfinite(v_cell) and abs(v_cell) > v_ref * 0.6)
                    else "#000"
                )
                ax.text(
                    j + 0.5, i + 0.5, mark,
                    ha="center", va="center",
                    fontsize=DEFAULT_FONT_SIZE - 1,
                    color=txt_color, zorder=4,
                )

    # Horizontal dividers (e.g. up/down group separator)
    if divider_rows:
        for row_idx in divider_rows:
            if 0 < row_idx < n_rows:
                ax.axhline(row_idx, color=divider_color, lw=divider_lw, alpha=0.7)

    ax.set_xlabel(xlabel, fontsize=DEFAULT_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=DEFAULT_FONT_SIZE)
    if title:
        ax.set_title(title, fontsize=DEFAULT_FONT_SIZE)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0, labelsize=DEFAULT_FONT_SIZE)

    cbar = fig.colorbar(mesh, ax=ax, fraction=0.025, pad=0.02)
    if cbar_label:
        cbar.set_label(cbar_label, fontsize=DEFAULT_FONT_SIZE)
    cbar.ax.tick_params(labelsize=DEFAULT_FONT_SIZE)

    fig.tight_layout()
    return fig, ax


__all__ = ["tile_heatmap"]
