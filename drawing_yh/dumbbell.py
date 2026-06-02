#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reusable dumbbell chart helpers.
"""
from __future__ import annotations

from collections.abc import Mapping

import matplotlib.pyplot as plt
import pandas as pd

from . import DEFAULT_FONT_SIZE, DOUBLE_COL_IN
from .io import save_fig
from .layout import compute_figsize
from .palettes import OKABE_ITO


YOUNG_COLOR = OKABE_ITO[4]
OLD_COLOR = OKABE_ITO[5]
UP_COLOR = "#c0392b"
DOWN_COLOR = "#2c5fa8"


def _left_margin(labels: pd.Series) -> float:
    max_len = int(labels.astype(str).str.len().max()) if len(labels) else 12
    return min(0.43, max(0.25, 0.18 + 0.0058 * max_len))


def _delta_format(xmax: float) -> str:
    if xmax < 0.02:
        return ".4f"
    if xmax < 0.2:
        return ".3f"
    return ".2f"


def render_mean_dumbbell(
    plot_df: pd.DataFrame,
    *,
    out_base: str,
    title: str,
    xlabel: str,
    subtitle: str | None = None,
    edge_col: str = "edge",
    young_col: str = "sample_mean_young",
    old_col: str = "sample_mean_old",
    delta_col: str = "sample_delta_old_minus_young",
    sig_col: str = "age_sig",
    tick_colors: list | None = None,
    tick_segments: list | None = None,
    width: float = DOUBLE_COL_IN,
    fmt: str | None = None,
    label_pad_frac: float = 0.025,
    xlim_frac: float = 1.30,
    legend_ncol: int = 4,
    note: str | None = None,
    save_kwargs: Mapping | None = None,
) -> None:
    """Render a compact Young-mean to Old-mean dumbbell chart."""
    df = plot_df.reset_index(drop=True).copy()
    if df.empty:
        print(f"Skip empty plot: {out_base}")
        return

    fig_w, fig_h = compute_figsize(
        n_items=len(df),
        width=width,
        per_item_h=0.24,
        base_h=1.55,
        min_h=3.0,
        max_h=8.2,
    )
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    xmax = max(float(df[young_col].max()), float(df[old_col].max()))
    xmax = max(xmax, 1e-9)

    for i, row in df.iterrows():
        delta = row[delta_col]
        color = UP_COLOR if delta > 0 else DOWN_COLOR
        young = float(row[young_col])
        old = float(row[old_col])
        direction = 1 if old >= young else -1
        distance = abs(old - young)
        gap = min(max(xmax * 0.012, 0.001), distance * 0.28)
        start = young + direction * gap
        head = old - direction * gap
        if abs(head - start) > xmax * 0.003:
            ax.plot([start, head], [i, i], color=color, lw=1.55, alpha=0.78, zorder=1)
        ax.scatter(
            head,
            i,
            marker=">" if direction > 0 else "<",
            s=28,
            color=color,
            edgecolor=color,
            linewidth=0,
            alpha=0.90,
            zorder=2,
        )
        ax.scatter(
            row[young_col],
            i,
            s=28,
            color="white",
            edgecolor=YOUNG_COLOR,
            linewidth=1.0,
            zorder=3,
        )
        ax.scatter(
            row[old_col],
            i,
            s=32,
            color=OLD_COLOR,
            edgecolor="black",
            linewidth=0.4,
            zorder=4,
        )

    pad = max(xmax * label_pad_frac, 0.0005)
    fmt = fmt or _delta_format(xmax)

    for i, row in df.iterrows():
        delta = row[delta_col]
        color = UP_COLOR if delta > 0 else DOWN_COLOR
        sig = row[sig_col] if pd.notna(row[sig_col]) else "NA"
        ax.text(
            max(row[young_col], row[old_col]) + pad,
            i,
            f"{delta:+{fmt}}  {sig}",
            va="center",
            ha="left",
            fontsize=DEFAULT_FONT_SIZE,
            color=color,
            clip_on=False,
        )

    ax.set_yticks(range(len(df)))
    if tick_segments is not None:
        # 两段着色的行标签(发送名/接收名各自 lineage 色):用 offsetbox 拼成右对齐的伪 tick label
        from matplotlib.offsetbox import TextArea, HPacker, AnnotationBbox
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
        for _i, _segs in enumerate(tick_segments):
            _boxes = [
                TextArea(_t, textprops=dict(color=_c, fontsize=DEFAULT_FONT_SIZE, fontweight="bold"))
                for _t, _c in _segs
            ]
            _pack = HPacker(children=_boxes, align="baseline", pad=0, sep=0)
            ax.add_artist(
                AnnotationBbox(
                    _pack, (0, _i), xycoords=("axes fraction", "data"),
                    xybox=(-4, 0), boxcoords="offset points", box_alignment=(1.0, 0.5),
                    frameon=False, annotation_clip=False, pad=0,
                )
            )
    else:
        ax.set_yticklabels(df[edge_col], fontsize=DEFAULT_FONT_SIZE)
        if tick_colors is not None:
            for _tick, _c in zip(ax.get_yticklabels(), tick_colors):
                _tick.set_color(_c)
                _tick.set_fontweight("bold")
    ax.set_xlabel(xlabel, fontsize=DEFAULT_FONT_SIZE)
    ax.set_xlim(0, xmax * xlim_frac + pad)
    ax.set_ylim(-0.6, len(df) - 0.4)
    title_text = f"Dumbbell chart — {title}"
    if subtitle:
        title_text += f"\n{subtitle}"
    ax.set_title(title_text, loc="left", fontsize=DEFAULT_FONT_SIZE, pad=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#dddddd", lw=0.5)

    if note:
        ax.text(
            0,
            -0.13,
            note,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=DEFAULT_FONT_SIZE,
            color="#444444",
        )

    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor="white",
            markeredgecolor=YOUNG_COLOR,
            label="Young mean",
            markersize=5,
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=OLD_COLOR,
            markeredgecolor="black",
            label="Old mean",
            markersize=5,
        ),
        plt.Line2D([0, 1], [0, 0], color=UP_COLOR, lw=1.5, marker=">", markevery=[1], label="Old > Young"),
        plt.Line2D([0, 1], [0, 0], color=DOWN_COLOR, lw=1.5, marker="<", markevery=[1], label="Old < Young"),
    ]
    legend_y = -0.24 if len(df) <= 10 else -0.09
    bottom = 0.30 if len(df) <= 10 else 0.14
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, legend_y),
        ncol=legend_ncol,
        frameon=False,
        fontsize=DEFAULT_FONT_SIZE,
    )

    plt.subplots_adjust(
        left=_left_margin(df[edge_col]),
        right=0.95,
        top=0.89,
        bottom=bottom,
    )
    save_fig(fig, out_base + ".pdf", also=(".png", ".svg"), **(save_kwargs or {}))
    plt.close(fig)
    print(f"Wrote: {out_base}.pdf/.png/.svg")
