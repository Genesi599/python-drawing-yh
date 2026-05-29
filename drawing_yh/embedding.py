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


__all__ = [
    "GREY_RED",
    "resolve_vlim",
    "auto_ncols",
    "scatter_embedding",
    "add_embedding_axes",
    "feature_plot",
]
