# -*- coding: utf-8 -*-
"""
UMAP value-overlay: colour cells on an embedding by a per-population
continuous value (e.g. NB log2FC) with diverging colormap, centroid
numbers, and a right-side key table.

Intended for differential-abundance UMAP panels (aging analysis),
but the interface is fully generic — any per-cell float value works.

Key functions
-------------
draw_da_umap_panel(ax, ax_key, x, y, cell_values, cell_labels, **kw)
    Low-level primitive: one scatter + one key, on axes you supply.

umap_value_panels(panels, **kw) -> (fig, cbar_ax)
    High-level: creates the figure, lays out ncols panels, adds a
    shared bottom colorbar.  Returns (fig, cbar_ax).

Example — multi-panel DA UMAP
------------------------------
    panels = []
    for lin in lineage_order:
        df = subcluster_data[lin]
        panels.append(dict(
            x=df.subUMAP1, y=df.subUMAP2,
            cell_values=df.cell_type_final.map(lfc_by_type),
            cell_labels=df.cell_type_final,
            pop_order=subtypes[lin],
            key_lines=build_key_lines(lin, da_within, disp_map),
            title=lin, title_color=lineage_color[lin],
        ))
    fig, cax = umap_value_panels(panels, ncols=2, clip=1.5,
                                  cbar_label="log2FC (within-lineage, old/young)")
"""
from __future__ import annotations

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from . import DEFAULT_FONT_SIZE


# ── default colormap (RdBu 7-stop, matches brain-aging convention) ────────────
_RDBU7 = ["#2166AC", "#4393C3", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#D6604D", "#B2182B"]
DA_CMAP = LinearSegmentedColormap.from_list("da_rdbu", _RDBU7)

GRAY_NA = (0.86, 0.86, 0.86, 1.0)   # "not estimable" fill colour


# ── internal helpers ──────────────────────────────────────────────────────────

def _robust_centroid(x: np.ndarray, y: np.ndarray, trim: float = 0.65):
    """Robust centroid: global median → keep trim-quantile → re-median.

    Identical to ``drawing_yh.embedding.cluster_centroids`` logic, but
    returns a single (cx, cy) pair instead of a dict.
    """
    cx, cy = np.median(x), np.median(y)
    d = (x - cx) ** 2 + (y - cy) ** 2
    keep = d <= np.quantile(d, trim)
    return float(np.median(x[keep])), float(np.median(y[keep]))


def _tight_lims(ax, x: np.ndarray, y: np.ndarray, q: float = 1.0, pad: float = 0.06):
    """Shrink axes to the cell distribution, clipping outlier whiskers."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    xlo, xhi = np.percentile(x, [q, 100 - q])
    ylo, yhi = np.percentile(y, [q, 100 - q])
    mx = pad * (xhi - xlo) or 1.0
    my = pad * (yhi - ylo) or 1.0
    ax.set_xlim(xlo - mx, xhi + mx)
    ax.set_ylim(ylo - my, yhi + my)


def _draw_key_top(
    ax,
    lines: list[tuple[str, str, bool]],
    fs: float,
    gap: float = 1.5,
):
    """Fill *ax* with key rows packed top-to-bottom using physical line height.

    Parameters
    ----------
    lines
        List of ``(text, colour, bold)`` tuples.
    fs
        Font size in points.
    gap
        Line-height multiplier (default 1.5 = 1.5× font size).
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax_h_in = ax.get_position().height * ax.figure.get_figheight()
    lh = (fs * gap / 72.0) / max(ax_h_in, 1e-6)
    y = 1.0 - 0.5 * lh
    for txt, color, bold in lines:
        ax.text(
            0.0 if bold else 0.035, y, txt,
            ha="left", va="center",
            fontsize=fs, color=color,
            fontweight="bold" if bold else "normal",
        )
        y -= lh


# ── primitive ─────────────────────────────────────────────────────────────────

def draw_da_umap_panel(
    ax,
    ax_key,
    x,
    y,
    cell_values,
    cell_labels,
    *,
    norm,
    cmap=None,
    pop_order: list[str] | None = None,
    key_lines: list[tuple[str, str, bool]] | None = None,
    fs: int = DEFAULT_FONT_SIZE,
    title: str | None = None,
    title_color: str = "black",
    point_size: float = 4.5,
    tight_q: float = 1.0,
    tight_pad: float = 0.06,
    key_gap: float = 1.5,
) -> None:
    """Draw one DA-UMAP panel (scatter + centroid numbers) and its right-side key.

    Parameters
    ----------
    ax, ax_key
        Matplotlib axes: scatter plot and key text respectively.
    x, y
        Cell 2-D embedding coordinates (length N).
    cell_values
        Per-cell continuous value; NaN → grey ("not estimable"), length N.
    cell_labels
        Per-cell population label (str), length N.
    norm
        A shared ``TwoSlopeNorm``; typically
        ``TwoSlopeNorm(vmin=-clip, vcenter=0, vmax=clip)``.
        Passing a shared norm keeps colours comparable across panels.
    cmap
        Colormap (default: module-level ``DA_CMAP``, 7-stop RdBu diverging).
    pop_order
        Ordered list of population names for centroid numbering.
        Populations absent in the data are silently skipped.
        If *None*, populations are ordered by first appearance.
    key_lines
        Explicit key content: list of ``(text, colour, bold)`` tuples.
        Build with :func:`build_key_lines` or supply manually.
        If *None*, ``ax_key`` is hidden (no key rendered).
    fs
        Font size in points (default 8).
    title, title_color
        Panel title (top-left, coloured).
    point_size
        Scatter marker area in pt².
    tight_q, tight_pad
        Axis-range trimming parameters passed to :func:`_tight_lims`.
    key_gap
        Line-height multiplier for the key rows.
    """
    _cmap = cmap if cmap is not None else DA_CMAP
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    cv = np.asarray(cell_values, dtype=float)
    cl = np.asarray(cell_labels, dtype=str)

    nan_m = np.isnan(cv)
    clipped = np.clip(cv, norm.vmin, norm.vmax)
    rgba = np.empty((len(cv), 4))
    rgba[nan_m] = GRAY_NA
    if (~nan_m).any():
        rgba[~nan_m] = _cmap(norm(clipped[~nan_m]))

    # gray first (behind), then coloured cells on top
    if nan_m.any():
        ax.scatter(x[nan_m], y[nan_m], s=point_size, c=[GRAY_NA],
                   linewidths=0, alpha=0.6, rasterized=True, zorder=1)
    ax.scatter(x[~nan_m], y[~nan_m], s=point_size, c=rgba[~nan_m],
               linewidths=0, alpha=0.9, rasterized=True, zorder=2)

    _tight_lims(ax, x, y, q=tight_q, pad=tight_pad)

    pops = (
        [p for p in pop_order if (cl == p).any()]
        if pop_order is not None
        else list(pd.unique(cl[~nan_m]))
    )
    for i, pop in enumerate(pops, 1):
        m = cl == pop
        if not m.any():
            continue
        cx, cy = _robust_centroid(x[m], y[m])
        t = ax.text(cx, cy, str(i), ha="center", va="center",
                    fontsize=fs, fontweight="bold", color="black", zorder=5)
        t.set_path_effects([pe.withStroke(linewidth=2.2, foreground="white")])

    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_edgecolor("#bbbbbb")
        sp.set_linewidth(0.8)
    if title:
        ax.set_title(title, fontsize=fs, loc="left", color=title_color, pad=3)

    if key_lines:
        _draw_key_top(ax_key, key_lines, fs, key_gap)
    else:
        ax_key.axis("off")


# ── key-line builder ──────────────────────────────────────────────────────────

def build_key_lines(
    pop_order: list[str],
    da: "pd.DataFrame",
    *,
    header: str | None = None,
    header_color: str = "black",
    display: dict[str, str] | None = None,
    pval_col: str = "p_value",
    fdr_col: str = "fdr",
    yp_col: str = "young_mean_prop",
    op_col: str = "old_mean_prop",
    lfc_col: str = "log2FC",
    idx_col: str | None = None,
    pct_scale: float = 100.0,
    pval_thr: float = 0.1,
    fdr_thr: float = 0.05,
) -> list[tuple[str, str, bool]]:
    """Build ``key_lines`` for :func:`draw_da_umap_panel` from a DA table.

    Parameters
    ----------
    pop_order
        Ordered list of population names matching the UMAP.
    da
        DataFrame with DA results.  Expected columns (overridable):
        ``log2FC``, ``p_value``, ``fdr``, ``young_mean_prop``, ``old_mean_prop``.
        Index (or ``idx_col``) must match the values in ``pop_order``.
    header
        Optional bold first line (e.g. ``"Glia  (within-lineage)"``).
    header_color
        Colour for the header line.
    display
        Optional ``{population: abbreviated_name}`` mapping.
    pval_col, fdr_col, yp_col, op_col, lfc_col
        Column names in *da*.
    idx_col
        Column to use as index; if *None*, uses ``da.index``.
    pct_scale
        Multiply proportions by this to get percentages (default 100).
    pval_thr, fdr_thr
        Raw-p threshold for ``†`` and FDR threshold for ``*``.

    Returns
    -------
    list of (text, colour, bold) tuples ready for ``draw_da_umap_panel``.
    """
    if idx_col is not None:
        da = da.set_index(idx_col)

    lines: list[tuple[str, str, bool]] = []
    if header:
        lines.append((header, header_color, True))

    for i, pop in enumerate(pop_order, 1):
        dname = display.get(pop, pop) if display else pop
        if pop not in da.index:
            lines.append((f"{i}  {dname}  (NE)", "#222222", False))
            continue
        row = da.loc[pop]
        lfc = row.get(lfc_col, np.nan)
        pval = row.get(pval_col, np.nan)
        fdr = row.get(fdr_col, np.nan)
        yp = row.get(yp_col, np.nan)
        op = row.get(op_col, np.nan)

        lfc_str = f"{lfc:+.1f}" if pd.notna(lfc) else "NE"
        pct_str = (
            f"{yp * pct_scale:.0f}→{op * pct_scale:.0f}%"
            if pd.notna(yp) and pd.notna(op) else ""
        )
        mark = (
            " *" if (pd.notna(fdr) and fdr < fdr_thr) else
            (" †" if (pd.notna(pval) and pval < pval_thr) else "")
        )
        extra = f"; {pct_str}" if pct_str else ""
        lines.append((f"{i}  {dname}  ({lfc_str}{extra}){mark}", "#222222", False))

    return lines


# ── high-level figure builder ─────────────────────────────────────────────────

def umap_value_panels(
    panels: list[dict],
    *,
    ncols: int = 2,
    panel_w: float = 5.2,
    panel_h: float = 3.4,
    clip: float = 1.5,
    cmap=None,
    cbar_label: str = "log2FC (old / young)",
    cbar_pos: tuple[float, float, float, float] | None = None,
    figure_title: str | None = None,
    figure_subtitle: str | None = None,
    fs: int = DEFAULT_FONT_SIZE,
):
    """Create a multi-panel (or single-panel) DA UMAP figure.

    Parameters
    ----------
    panels
        List of panel-specification dicts.  Each dict may contain:

        ``x, y``
            Cell coordinates (required).
        ``cell_values``
            Per-cell float; NaN → grey (required).
        ``cell_labels``
            Per-cell population label (required).
        ``pop_order``
            Ordered list of population names for centroid numbering.
        ``key_lines``
            Explicit ``(text, colour, bold)`` tuples; use
            :func:`build_key_lines` to generate from a DA table.
        ``title``
            Panel title (top-left).
        ``title_color``
            Colour for the panel title (default black).

    ncols
        Number of columns in the panel grid.
    panel_w, panel_h
        Width / height per grid cell in inches.
    clip
        Colormap saturation point (±clip log2FC).
    cmap
        Colormap; default ``DA_CMAP`` (7-stop RdBu diverging).
    cbar_label
        Shared colorbar label.
    cbar_pos
        ``[left, bottom, width, height]`` in figure coordinates for the
        shared colorbar.  Auto-positioned at the bottom centre when *None*.
    figure_title, figure_subtitle
        Optional super-title and subtitle rendered above the grid.
    fs
        Font size in points.

    Returns
    -------
    (fig, cbar_ax)
    """
    n = len(panels)
    nrows = int(np.ceil(n / ncols))
    norm = TwoSlopeNorm(vmin=-clip, vcenter=0, vmax=clip)
    _cmap = cmap if cmap is not None else DA_CMAP

    has_header = bool(figure_title or figure_subtitle)
    GX0, GX1 = 0.012, 0.997
    GY0 = 0.085
    GY1 = 0.925 if has_header else 0.970
    cell_w = (GX1 - GX0) / ncols
    cell_h = (GY1 - GY0) / nrows

    fig = plt.figure(figsize=(panel_w * ncols, panel_h * nrows))

    for pi, pspec in enumerate(panels):
        r, c = divmod(pi, ncols)
        x0 = GX0 + c * cell_w
        y0 = GY1 - (r + 1) * cell_h
        ax = fig.add_axes([
            x0 + 0.010 * cell_w, y0 + 0.06 * cell_h,
            0.45 * cell_w,        0.78 * cell_h,
        ])
        ax_key = fig.add_axes([
            x0 + 0.47 * cell_w, y0 + 0.02 * cell_h,
            0.52 * cell_w,       0.92 * cell_h,
        ])
        draw_da_umap_panel(
            ax, ax_key,
            pspec["x"], pspec["y"],
            pspec["cell_values"], pspec["cell_labels"],
            norm=norm, cmap=_cmap,
            pop_order=pspec.get("pop_order"),
            key_lines=pspec.get("key_lines"),
            title=pspec.get("title"),
            title_color=pspec.get("title_color", "black"),
            fs=fs,
        )

    for pi in range(n, nrows * ncols):
        r, c = divmod(pi, ncols)
        x0 = GX0 + c * cell_w
        y0 = GY1 - (r + 1) * cell_h
        fig.add_axes([x0, y0, cell_w, cell_h]).axis("off")

    if cbar_pos is None:
        cbar_pos = [0.40, 0.030, 0.20, 0.013]
    cax = fig.add_axes(cbar_pos)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=_cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", extend="both")
    cb.set_label(cbar_label, fontsize=fs)
    cb.ax.tick_params(labelsize=fs)

    if figure_title:
        fig.text(GX0, GY1 + 0.013, figure_title,
                 ha="left", va="bottom", fontsize=fs + 1, fontweight="bold")
    if figure_subtitle:
        fig.text(GX0, GY1 - 0.005, figure_subtitle,
                 ha="left", va="top", fontsize=fs - 1, color="#555555")

    return fig, cax


__all__ = [
    "DA_CMAP",
    "GRAY_NA",
    "draw_da_umap_panel",
    "build_key_lines",
    "umap_value_panels",
]
