"""Standard two-dimensional heatmap helpers for drawing_yh."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, Sequence

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

from drawing_yh import DEFAULT_FONT_SIZE, DOUBLE_COL_IN
from drawing_yh.io import save_fig


DEFAULT_MISSING_COLOR = "#f2f4f7"
DEFAULT_LINE_COLOR = "white"
DEFAULT_DOT_COLOR = "#111827"


def compute_heatmap_figsize(
    n_rows: int,
    n_cols: int,
    *,
    width: float | None = DOUBLE_COL_IN,
    cell_width: float = 0.36,
    cell_height: float = 0.24,
    base_height: float = 1.4,
    min_height: float = 2.2,
    max_height: float = 8.0,
) -> tuple[float, float]:
    """Compute a compact, standards-friendly heatmap size in inches."""
    if width is None:
        width = float(np.clip(n_cols * cell_width + 2.2, 3.35, 10.0))
    height = float(np.clip(n_rows * cell_height + base_height, min_height, max_height))
    return float(width), height


def long_to_matrix(
    data: pd.DataFrame,
    *,
    row: str,
    column: str,
    value: str,
    row_order: Sequence | None = None,
    column_order: Sequence | None = None,
    aggfunc: str = "first",
) -> pd.DataFrame:
    """Pivot a long-form table into the matrix expected by ``plot_tile_heatmap``."""
    matrix = data.pivot_table(index=row, columns=column, values=value, aggfunc=aggfunc)
    if row_order is not None:
        matrix = matrix.reindex(list(row_order))
    if column_order is not None:
        matrix = matrix.reindex(columns=list(column_order))
    return matrix


def significance_stars(
    pvalues: pd.DataFrame,
    thresholds: Sequence[tuple[float, str]] = ((0.001, "***"), (0.01, "**"), (0.05, "*")),
) -> pd.DataFrame:
    """Convert a p-value matrix into standard star labels."""
    labels = pd.DataFrame("", index=pvalues.index, columns=pvalues.columns)
    for threshold, label in sorted(thresholds, key=lambda item: item[0], reverse=True):
        labels = labels.mask(pvalues < threshold, label)
    return labels


def _as_dataframe(data) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        return data.copy()
    return pd.DataFrame(data)


def _display_labels(labels: Iterable, mapping: Mapping | None) -> list[str]:
    if mapping is None:
        return [str(label) for label in labels]
    return [str(mapping.get(label, label)) for label in labels]


def _build_norm(
    values: np.ndarray,
    *,
    center: float | None,
    vmin: float | None,
    vmax: float | None,
):
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        finite = np.array([0.0])

    if center is not None:
        span = float(np.nanmax(np.abs(finite - center)))
        if not np.isfinite(span) or span == 0:
            span = 1.0
        if vmin is None:
            vmin = center - span
        if vmax is None:
            vmax = center + span
        if vmin < center < vmax:
            return mcolors.TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)

    if vmin is None:
        vmin = float(np.nanmin(finite))
    if vmax is None:
        vmax = float(np.nanmax(finite))
    if vmin == vmax:
        vmin -= 1.0
        vmax += 1.0
    return mcolors.Normalize(vmin=vmin, vmax=vmax)


def _add_vector_colorbar(
    fig,
    ax,
    *,
    cmap,
    norm,
    label: str | None,
    ticks: Sequence[float] | None,
    tick_format: str,
):
    fig.canvas.draw()
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + 0.018, pos.y0, 0.018, pos.height])
    cax.set_xlim(0, 1)
    cax.set_ylim(norm.vmin, norm.vmax)

    bins = np.linspace(norm.vmin, norm.vmax, 96)
    for lo, hi in zip(bins[:-1], bins[1:]):
        mid = (lo + hi) / 2
        cax.add_patch(
            Rectangle(
                (0, lo),
                1,
                hi - lo,
                facecolor=cmap(norm(mid)),
                edgecolor=cmap(norm(mid)),
                linewidth=0,
            )
        )

    if ticks is None:
        ticks = np.linspace(norm.vmin, norm.vmax, 5)
    cax.set_xticks([])
    cax.set_yticks(ticks)
    cax.set_yticklabels([tick_format.format(tick) for tick in ticks])
    cax.tick_params(axis="y", labelsize=DEFAULT_FONT_SIZE, length=2, pad=2)
    cax.yaxis.tick_right()
    cax.yaxis.set_label_position("right")
    if label:
        cax.set_ylabel(label, fontsize=DEFAULT_FONT_SIZE)
    for spine in cax.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("#111827")
    return cax


def plot_tile_heatmap(
    data,
    *,
    pvalues=None,
    annotations=None,
    row_order: Sequence | None = None,
    column_order: Sequence | None = None,
    row_label_map: Mapping | None = None,
    column_label_map: Mapping | None = None,
    cmap: str = "RdBu_r",
    center: float | None = 0.0,
    vmin: float | None = None,
    vmax: float | None = None,
    missing_color: str = DEFAULT_MISSING_COLOR,
    linecolor: str = DEFAULT_LINE_COLOR,
    linewidth: float = 0.45,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    cbar_label: str | None = None,
    cbar_ticks: Sequence[float] | None = None,
    cbar_tick_format: str = "{:.1f}",
    show_colorbar: bool = True,
    significance_mode: str | None = None,
    marker_pvalue_threshold: float = 0.05,
    dot_size: float = 8,
    dot_color: str = DEFAULT_DOT_COLOR,
    annotation_color: str = DEFAULT_DOT_COLOR,
    annotation_fontsize: float | None = None,
    xtick_rotation: float = 45,
    ytick_rotation: float = 0,
    figsize: tuple[float, float] | None = None,
    width: float | None = DOUBLE_COL_IN,
    cell_width: float = 0.36,
    cell_height: float = 0.24,
):
    """Render a vector tile heatmap with optional p-value stars or dots."""
    matrix = _as_dataframe(data)
    if row_order is not None:
        matrix = matrix.reindex(list(row_order))
    if column_order is not None:
        matrix = matrix.reindex(columns=list(column_order))

    n_rows, n_cols = matrix.shape
    if figsize is None:
        figsize = compute_heatmap_figsize(
            n_rows,
            n_cols,
            width=width,
            cell_width=cell_width,
            cell_height=cell_height,
        )

    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(right=0.86 if show_colorbar else 0.98)
    cmap_obj = plt.get_cmap(cmap)
    norm = _build_norm(matrix.to_numpy(dtype=float), center=center, vmin=vmin, vmax=vmax)

    for i, row_name in enumerate(matrix.index):
        for j, col_name in enumerate(matrix.columns):
            value = matrix.loc[row_name, col_name]
            facecolor = missing_color if pd.isna(value) else cmap_obj(norm(float(value)))
            ax.add_patch(
                Rectangle(
                    (j, i),
                    1,
                    1,
                    facecolor=facecolor,
                    edgecolor=linecolor,
                    linewidth=linewidth,
                )
            )

    if pvalues is not None and significance_mode is None:
        significance_mode = "star"
    if pvalues is not None:
        pmat = _as_dataframe(pvalues).reindex(index=matrix.index, columns=matrix.columns)
        if significance_mode == "star":
            annotations = significance_stars(pmat)
        elif significance_mode == "dot":
            for i, row_name in enumerate(matrix.index):
                for j, col_name in enumerate(matrix.columns):
                    pval = pmat.loc[row_name, col_name]
                    if pd.notna(pval) and pval < marker_pvalue_threshold:
                        ax.scatter(
                            j + 0.5,
                            i + 0.5,
                            s=dot_size,
                            color=dot_color,
                            linewidths=0,
                            zorder=3,
                        )

    if annotations is not None and significance_mode != "dot":
        amat = _as_dataframe(annotations).reindex(index=matrix.index, columns=matrix.columns)
        fontsize = annotation_fontsize or DEFAULT_FONT_SIZE
        for i, row_name in enumerate(matrix.index):
            for j, col_name in enumerate(matrix.columns):
                label = amat.loc[row_name, col_name]
                if pd.notna(label) and str(label):
                    ax.text(
                        j + 0.5,
                        i + 0.5,
                        str(label),
                        ha="center",
                        va="center",
                        fontsize=fontsize,
                        color=annotation_color,
                        zorder=4,
                    )

    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(
        _display_labels(matrix.columns, column_label_map),
        rotation=xtick_rotation,
        ha="right" if xtick_rotation else "center",
    )
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(_display_labels(matrix.index, row_label_map), rotation=ytick_rotation)
    ax.tick_params(labelsize=DEFAULT_FONT_SIZE, length=0)
    ax.set_xlabel(xlabel or "", fontsize=DEFAULT_FONT_SIZE)
    ax.set_ylabel(ylabel or "", fontsize=DEFAULT_FONT_SIZE)
    if title:
        ax.set_title(title, fontsize=DEFAULT_FONT_SIZE, pad=8)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if show_colorbar:
        _add_vector_colorbar(
            fig,
            ax,
            cmap=cmap_obj,
            norm=norm,
            label=cbar_label,
            ticks=cbar_ticks,
            tick_format=cbar_tick_format,
        )
    return fig, ax


def save_tile_heatmap(
    data,
    path: str | Path,
    *,
    also: Iterable[str] = (".pdf", ".png"),
    close: bool = True,
    **kwargs,
) -> list[Path]:
    """Render and save a standard tile heatmap."""
    fig, _ = plot_tile_heatmap(data, **kwargs)
    written = save_fig(fig, path, also=also)
    if close:
        plt.close(fig)
    return written


def plot_long_heatmap(
    data: pd.DataFrame,
    *,
    row: str,
    column: str,
    value: str,
    pvalue: str | None = None,
    row_order: Sequence | None = None,
    column_order: Sequence | None = None,
    aggfunc: str = "first",
    **kwargs,
):
    """Render a standard heatmap from a long-form table."""
    matrix = long_to_matrix(
        data,
        row=row,
        column=column,
        value=value,
        row_order=row_order,
        column_order=column_order,
        aggfunc=aggfunc,
    )
    pvalues = None
    if pvalue is not None:
        pvalues = long_to_matrix(
            data,
            row=row,
            column=column,
            value=pvalue,
            row_order=matrix.index,
            column_order=matrix.columns,
            aggfunc=aggfunc,
        )
    return plot_tile_heatmap(matrix, pvalues=pvalues, **kwargs)


def save_long_heatmap(
    data: pd.DataFrame,
    path: str | Path,
    *,
    row: str,
    column: str,
    value: str,
    pvalue: str | None = None,
    row_order: Sequence | None = None,
    column_order: Sequence | None = None,
    aggfunc: str = "first",
    also: Iterable[str] = (".pdf", ".png"),
    close: bool = True,
    **kwargs,
) -> list[Path]:
    """Render and save a standard heatmap from a long-form table."""
    fig, _ = plot_long_heatmap(
        data,
        row=row,
        column=column,
        value=value,
        pvalue=pvalue,
        row_order=row_order,
        column_order=column_order,
        aggfunc=aggfunc,
        **kwargs,
    )
    written = save_fig(fig, path, also=also)
    if close:
        plt.close(fig)
    return written


__all__ = [
    "compute_heatmap_figsize",
    "long_to_matrix",
    "significance_stars",
    "plot_tile_heatmap",
    "save_tile_heatmap",
    "plot_long_heatmap",
    "save_long_heatmap",
]
