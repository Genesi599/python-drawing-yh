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

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

from drawing_yh import DEFAULT_FONT_SIZE, DOUBLE_COL_IN


DEFAULT_ROW_COLOR = "#6baed6"
DIVERGING_COLORS = ("#2166ac", "#ffffff", "#b2182b")


def wrap_text(text: str, max_chars: int = 40) -> str:
    """Word-wrap ``text`` to lines of at most ``max_chars`` (no hyphenation)."""
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
        fig_h = fig_height if fig_height is not None else float(np.clip(0.22 * G + 1.2, 2.2, 9.0))
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

    cbar = fig.colorbar(im, ax=ax_heat, orientation="horizontal",
                        fraction=0.045, pad=0.06, shrink=0.6)
    cbar.set_label(cbar_label, fontsize=font)
    cbar.ax.tick_params(labelsize=font)

    if row_color_strip and colors is not None:
        trans = transforms.blended_transform_factory(ax_heat.transAxes, ax_heat.transData)
        for gi, col in enumerate(colors):
            ax_heat.add_patch(Rectangle(
                (-strip_width, gi - 0.5), strip_width, 1.0,
                transform=trans, facecolor=col, edgecolor="none",
                linewidth=0, zorder=3, clip_on=False,
            ))

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
                wrapped = wrap_text(term, max_chars=wrap_chars)
                ax_bars.annotate(
                    wrapped, xy=(1, y), xycoords=trans_bar, xytext=(3, 0),
                    textcoords="offset points", ha="left", va="center",
                    fontsize=font, color="black", clip_on=False, linespacing=1,
                )
        ax_bars.set_xlabel(bar_xlabel, fontsize=font)
        ax_bars.grid(axis="x", linestyle=":", alpha=0.3)
        ax_bars.tick_params(axis="y", left=False, labelleft=False)
        ax_bars.tick_params(axis="x", labelsize=font)
        ax_bars.set_ylim(ax_heat.get_ylim())
        ax_bars.set_xlim(0, xmax * 1.25)
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
