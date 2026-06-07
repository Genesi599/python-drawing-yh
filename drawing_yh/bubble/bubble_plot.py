"""Generic bubble (dot) plot for categorical x/y axes.

bubble_plot(items, x_order, y_order, *, size_max, ...)
    Draw a grid of bubbles where x/y are categorical and bubble size encodes
    a numeric value.  Colours encode a third variable (category palette or
    diverging direction).

two_tone_ticklabels(ax, axis, positions, seg_lists, fontsize)
    Replace axis tick labels with multi-coloured AnnotationBbox segments.
    Returns the artist list needed for bbox_extra_artists in save_fig.

Both helpers are independent of CellChat – the caller supplies data.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import (
    AnnotationBbox,
    HPacker,
    TextArea,
    VPacker,
)


def two_tone_ticklabels(
    ax,
    axis: str,
    positions: Sequence,
    seg_lists: Sequence[Sequence[tuple[str, str]]],
    fontsize: int = 8,
):
    """Replace axis tick labels with multi-coloured segments.

    Parameters
    ----------
    ax : Axes
    axis : ``'x'`` or ``'y'``
    positions : sequence of numeric tick positions
    seg_lists : list of [(text, colour), ...] – one inner list per tick
    fontsize : int

    Returns
    -------
    list of Artist  (pass to ``bbox_extra_artists`` in ``savefig`` / ``save_fig``)
    """
    rot = 270 if axis == "x" else 0
    if axis == "y":
        ax.set_yticks(positions)
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
        xycoords = ("axes fraction", "data")
        box_align = (1.0, 0.5)
        xybox = (-5, 0)
    else:
        ax.set_xticks(positions)
        ax.set_xticklabels([])
        ax.tick_params(axis="x", length=0)
        xycoords = ("data", "axes fraction")
        box_align = (0.5, 1.0)
        xybox = (0, -5)

    artists = []
    for pos, segs in zip(positions, seg_lists):
        tas = [
            TextArea(t, textprops=dict(color=c, fontsize=fontsize, rotation=rot))
            for t, c in segs
        ]
        box = (HPacker if axis == "y" else VPacker)(
            children=tas, align="center", pad=0, sep=0
        )
        anchor = (0, pos) if axis == "y" else (pos, 0)
        ab = AnnotationBbox(
            box,
            anchor,
            xycoords=xycoords,
            xybox=xybox,
            boxcoords="offset points",
            box_alignment=box_align,
            frameon=False,
            annotation_clip=False,
            pad=0,
        )
        ax.add_artist(ab)
        artists.append(ab)
    return artists


def bubble_plot(
    items,
    x_order: Sequence,
    y_order: Sequence,
    *,
    size_max: float,
    base_size: float = 5.0,
    max_size: float = 72.0,
    alpha: float = 0.9,
    edgecolor: str = "#333333",
    linewidth: float = 0.3,
    title: str = "",
    x_tick_segs=None,
    y_tick_segs=None,
    x_fontsize: int = 8,
    y_fontsize: int = 8,
    grid_color: str = "#ECECEC",
    cell_w: float = 0.145,
    cell_h: float = 0.145,
    left_in: float = 1.55,
    right_in: float = 1.50,
    bottom_in: float = 2.10,
    top_in: float = 0.42,
    color_legend_handles=None,
    color_legend_title: str = "",
    size_legend_label: str = "prob",
    size_legend_values=None,
):
    """Draw a generic bubble plot with categorical x (columns) and y (rows).

    Parameters
    ----------
    items : iterable of (x_cat, y_cat, size_value, color)
        One bubble per item.  ``x_cat`` / ``y_cat`` must appear in
        ``x_order`` / ``y_order``.
    x_order : sequence
        Ordered column categories (left → right).
    y_order : sequence
        Ordered row categories (bottom → top, i.e. ``y_order[0]`` is the
        bottom-most row).
    size_max : float
        Reference maximum for the size scale.  **Share this value across
        panels** to keep dot sizes comparable.
    base_size : float
        Minimum marker area (scatter *s* units).
    max_size : float
        Maximum marker area increment (scatter *s* units).
    alpha : float
        Bubble transparency.
    edgecolor : str
        Bubble edge colour.
    linewidth : float
        Bubble edge line width.
    title : str
        Axes title.
    x_tick_segs : list[list[(str, str)]] | None
        Per-column two-tone label segments (length == ``len(x_order)``).
        If *None*, plain text tick labels are used.
    y_tick_segs : list[list[(str, str)]] | None
        Per-row two-tone label segments (length == ``len(y_order)``).
        If *None*, plain text tick labels are used.
    x_fontsize, y_fontsize : int
        Font sizes for the respective axis labels.
    grid_color : str
        Background grid line colour.
    cell_w, cell_h : float
        Cell width / height in inches.
    left_in, right_in, bottom_in, top_in : float
        Fixed outer margins in inches.
    color_legend_handles : list[Line2D] | None
        Legend handles for bubble colours.
    color_legend_title : str
    size_legend_label : str
        Title for the size legend.
    size_legend_values : list[float] | None
        Reference size values for the size legend; defaults to
        ``[0.25, 0.5, 1.0] × size_max``.

    Returns
    -------
    fig : Figure
    ax : Axes
    tick_artists : list[Artist]
        Pass to ``bbox_extra_artists`` in ``save_fig`` / ``savefig`` so
        tight-layout includes two-tone labels and legend boxes.
    """
    nC, nR = len(x_order), len(y_order)
    x_idx = {c: j for j, c in enumerate(x_order)}
    y_idx = {r: i for i, r in enumerate(y_order)}

    fig_w = left_in + nC * cell_w + right_in
    fig_h = bottom_in + nR * cell_h + top_in
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(
        left=left_in / fig_w,
        right=1 - right_in / fig_w,
        top=1 - top_in / fig_h,
        bottom=bottom_in / fig_h,
    )

    for x_cat, y_cat, size_val, color in items:
        xi = x_idx.get(x_cat)
        yi = y_idx.get(y_cat)
        if xi is None or yi is None:
            continue
        s = base_size + (float(size_val) / size_max) * max_size
        ax.scatter(
            xi,
            yi,
            s=s,
            color=color,
            edgecolor=edgecolor,
            linewidth=linewidth,
            alpha=alpha,
            zorder=3,
        )

    tick_artists: list = []

    if y_tick_segs is not None:
        tick_artists += two_tone_ticklabels(
            ax, "y", list(range(nR)), y_tick_segs, fontsize=y_fontsize
        )
    else:
        ax.set_yticks(range(nR))
        ax.set_yticklabels(list(y_order), fontsize=y_fontsize)
        ax.tick_params(axis="y", length=0)

    if x_tick_segs is not None:
        tick_artists += two_tone_ticklabels(
            ax, "x", list(range(nC)), x_tick_segs, fontsize=x_fontsize
        )
    else:
        ax.set_xticks(range(nC))
        ax.set_xticklabels(list(x_order), rotation=90, fontsize=x_fontsize)
        ax.tick_params(axis="x", length=0)

    ax.set_xlim(-0.5, nC - 0.5)
    ax.set_ylim(-0.6, nR - 0.4)
    ax.set_title(title, fontsize=8)
    ax.grid(True, color=grid_color, lw=0.4)
    ax.set_axisbelow(True)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)

    if color_legend_handles:
        leg1 = ax.legend(
            handles=color_legend_handles,
            title=color_legend_title,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            fontsize=8,
            title_fontsize=8,
            frameon=False,
            labelspacing=0.3,
            handletextpad=0.4,
        )
        ax.add_artist(leg1)
        tick_artists.append(leg1)

    if size_legend_values is None:
        size_legend_values = [round(size_max * f, 4) for f in (0.25, 0.5, 1.0)]
    sh = [
        Line2D(
            [0],
            [0],
            marker="o",
            ls="",
            mfc="#999999",
            mec="#333333",
            mew=0.3,
            ms=(base_size + (v / size_max) * max_size) ** 0.5 / 2.2,
            label=f"{v:g}",
        )
        for v in size_legend_values
    ]
    leg2 = ax.legend(
        handles=sh,
        title=size_legend_label,
        loc="lower left",
        bbox_to_anchor=(1.01, 0.0),
        fontsize=8,
        title_fontsize=8,
        frameon=False,
        labelspacing=1.0,
        handletextpad=0.6,
    )
    tick_artists.append(leg2)

    return fig, ax, tick_artists
