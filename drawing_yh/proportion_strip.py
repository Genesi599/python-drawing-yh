"""Generic boxplot + per-sample scatter for Young/Old proportion comparisons.

Reusable across tissues (brain, retina, etc.).  Call fraction_by_age_panel()
to draw on a supplied Axes; call fraction_by_age_figure() for a standalone figure.
"""
from __future__ import annotations

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt


DEFAULT_AGE_COLORS: dict[str, str] = {"Young": "#4C78A8", "Old": "#B55364"}


def _lighten(c: str, f: float = 0.78) -> tuple:
    r, g, b = mcolors.to_rgb(c)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


def fraction_by_age_panel(
    ax,
    df,
    labels: list[str],
    title: str,
    ylabel: str,
    *,
    age_col: str = "age_stage",
    val_col: str = "fraction",
    label_col: str = "label",
    age_order: tuple[str, str] = ("Young", "Old"),
    age_colors: dict[str, str] | None = None,
    short_label_fn=None,
    seed: int = 7,
    fs: float = 8,
) -> None:
    """Draw boxplot + per-sample scatter for Young/Old fraction comparisons.

    Parameters
    ----------
    ax : Axes
    df : DataFrame — long format with at least [age_col, val_col, label_col]
    labels : ordered list of label values to display left → right
    title, ylabel : axes title and y-axis label
    age_col : column holding age group strings (default "age_stage")
    val_col : column holding fraction values 0–1 (displayed as %)
    label_col : column holding cell-type / group names
    age_order : (left_age, right_age) — controls box offset direction
    age_colors : mapping age → hex color; defaults to blue/red
    short_label_fn : optional callable that shortens label strings for x-axis ticks
    seed : RNG seed for jitter
    fs : font size
    """
    if age_colors is None:
        age_colors = dict(DEFAULT_AGE_COLORS)
    sub = df[df[label_col].isin(labels)].copy()
    x_pos = {label: i for i, label in enumerate(labels)}
    rng = np.random.default_rng(seed)
    offsets = (-0.18, 0.18)

    for age, offset in zip(age_order, offsets):
        color = age_colors.get(age, "#999999")
        g = sub[sub[age_col] == age]
        data, pos = [], []
        for label in labels:
            vals = (
                g.loc[g[label_col] == label, val_col]
                .dropna()
                .to_numpy(dtype=float)
                * 100
            )
            if len(vals):
                data.append(vals)
                pos.append(x_pos[label] + offset)
        if data:
            ax.boxplot(
                data,
                positions=pos,
                widths=0.30,
                patch_artist=True,
                showfliers=False,
                zorder=2,
                boxprops=dict(facecolor=_lighten(color), edgecolor=color, lw=1.0),
                medianprops=dict(color=color, lw=1.5),
                whiskerprops=dict(color=color, lw=1.0),
                capprops=dict(color=color, lw=1.0),
            )
        xs = (
            np.array([x_pos[x] for x in g[label_col]], dtype=float)
            + offset
            + rng.normal(0, 0.04, len(g))
        )
        ys = g[val_col].to_numpy(dtype=float) * 100
        ax.scatter(
            xs, ys, s=14, color=color, label=age, alpha=0.9,
            edgecolors="white", linewidths=0.3, zorder=3,
        )

    slabels = [short_label_fn(x) if short_label_fn else x for x in labels]
    ax.set_title(title, fontsize=fs)
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(slabels, rotation=35, ha="right", fontsize=fs)
    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.tick_params(labelsize=fs)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.6, alpha=0.7)


def fraction_by_age_figure(
    df,
    labels: list[str],
    title: str,
    ylabel: str,
    outbase,
    *,
    legend_loc: str = "upper center",
    legend_anchor: tuple[float, float] = (0.5, -0.27),
    fig_w: float | None = None,
    fig_h: float = 3.2,
    **panel_kw,
):
    """Standalone figure wrapper around fraction_by_age_panel.

    Parameters
    ----------
    outbase : Path-like without extension; saves .png/.pdf/.svg via save_fig.
    fig_w : override figure width; defaults to auto-scaled from len(labels).
    panel_kw : forwarded to fraction_by_age_panel.

    Returns
    -------
    outbase (Path) — for chaining / copy-to-report logic.
    """
    from drawing_yh import save_fig  # lazy import avoids circular at module load

    if fig_w is None:
        fig_w = max(3.2, min(9.0, 0.68 * len(labels) + 1.4))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fraction_by_age_panel(ax, df, labels, title, ylabel, **panel_kw)
    handles, labs = ax.get_legend_handles_labels()
    seen = dict(zip(labs, handles))
    ax.legend(
        seen.values(), seen.keys(),
        loc=legend_loc, bbox_to_anchor=legend_anchor,
        ncol=2, frameon=False,
        fontsize=panel_kw.get("fs", 8),
    )
    fig.tight_layout()
    from pathlib import Path
    outbase = Path(outbase)
    save_fig(fig, outbase.with_suffix(".png"), also=(".pdf", ".svg"))
    plt.close(fig)
    return outbase
