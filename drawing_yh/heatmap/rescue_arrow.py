"""Rescue-style directional arrow heatmap templates."""
from __future__ import annotations

from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
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

STAIRCASE_PANEL_CONFIG = {
    "Aging Up-Rescued": {
        "aging_color": "#E07020",
        "rescue_color": "#1B7F79",
        "title": "Aging Up-Rescued",
    },
    "Aging Down-Rescued": {
        "aging_color": "#2166AC",
        "rescue_color": "#D0627D",
        "title": "Aging Down-Rescued",
    },
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


def _tint(hex_color: str, alpha: float = 0.30) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[:2], 16), int(h[2:4], 16), int(h[4:], 16)
    r = int(r + (255 - r) * (1 - alpha))
    g = int(g + (255 - g) * (1 - alpha))
    b = int(b + (255 - b) * (1 - alpha))
    return f"#{r:02x}{g:02x}{b:02x}"


def _block_arrow(
    ax,
    x: float,
    y_from: float,
    y_to: float,
    color,
    *,
    shaft_width: float = 0.13,
    head_width: float = 0.30,
    head_frac: float = 0.32,
    head_min: float = 0.06,
    tail_round: float = 0.05,
    alpha: float = 1.0,
    zorder: int = 5,
) -> None:
    dy = y_to - y_from
    if abs(dy) < 1e-6:
        return
    sgn = 1 if dy > 0 else -1
    length = abs(dy)
    head_len = max(min(length * head_frac, length * 0.55), min(head_min, length * 0.9))
    shaft_len = length - head_len
    shaft_half = shaft_width / 2
    head_half = head_width / 2
    radius = min(tail_round, max(shaft_len * 0.45, 0.0), shaft_half * 0.95)
    y_tail = y_from
    y_shaft_top = y_from + sgn * shaft_len
    y_apex = y_from + sgn * length

    if shaft_len <= radius + 1e-6:
        verts = [
            (x - shaft_half, y_tail),
            (x - shaft_half, y_shaft_top),
            (x - head_half, y_shaft_top),
            (x, y_apex),
            (x + head_half, y_shaft_top),
            (x + shaft_half, y_shaft_top),
            (x + shaft_half, y_tail),
            (x - shaft_half, y_tail),
        ]
        codes = [MplPath.MOVETO] + [MplPath.LINETO] * 6 + [MplPath.CLOSEPOLY]
    else:
        verts = [
            (x - shaft_half, y_tail + sgn * radius),
            (x - shaft_half, y_shaft_top),
            (x - head_half, y_shaft_top),
            (x, y_apex),
            (x + head_half, y_shaft_top),
            (x + shaft_half, y_shaft_top),
            (x + shaft_half, y_tail + sgn * radius),
            (x + shaft_half, y_tail),
            (x, y_tail),
            (x - shaft_half, y_tail),
            (x - shaft_half, y_tail + sgn * radius),
        ]
        codes = [
            MplPath.MOVETO,
            MplPath.LINETO,
            MplPath.LINETO,
            MplPath.LINETO,
            MplPath.LINETO,
            MplPath.LINETO,
            MplPath.LINETO,
            MplPath.CURVE3,
            MplPath.CURVE3,
            MplPath.CURVE3,
            MplPath.CURVE3,
        ]

    ax.add_patch(
        PathPatch(
            MplPath(verts, codes),
            facecolor=color,
            edgecolor="white",
            linewidth=0.45,
            joinstyle="miter",
            alpha=alpha,
            zorder=zorder,
        )
    )


class _ArrowProxy:
    def __init__(self, color: str, direction: int):
        self.color = color
        self.direction = direction


class _ArrowHandler(HandlerBase):
    def create_artists(self, legend, handle, xdescent, ydescent, width, height, fontsize, trans):
        center_x = width * 0.5 - xdescent
        pad = height * 0.12
        if handle.direction > 0:
            y_tail = -ydescent + pad
            y_apex = height - ydescent - pad
        else:
            y_tail = height - ydescent - pad
            y_apex = -ydescent + pad
        sgn = 1 if y_apex > y_tail else -1
        length = abs(y_apex - y_tail)
        head_len = length * 0.42
        shaft_len = length - head_len
        shaft_half = 1.1
        head_half = 2.8
        y_shaft = y_tail + sgn * shaft_len
        verts = [
            (center_x - shaft_half, y_tail),
            (center_x - shaft_half, y_shaft),
            (center_x - head_half, y_shaft),
            (center_x, y_apex),
            (center_x + head_half, y_shaft),
            (center_x + shaft_half, y_shaft),
            (center_x + shaft_half, y_tail),
            (center_x - shaft_half, y_tail),
        ]
        codes = [MplPath.MOVETO] + [MplPath.LINETO] * 6 + [MplPath.CLOSEPOLY]
        patch = PathPatch(MplPath(verts, codes), facecolor=handle.color, edgecolor="none", transform=trans)
        return [patch]


def _padj_color(padj: float | None, cmap, *, sig_cutoff: float, cap: float, grey: str = "#cfcfcf"):
    if padj is None or not np.isfinite(padj):
        return grey
    if padj >= sig_cutoff:
        return grey
    score = -np.log10(max(float(padj), 1e-300))
    frac = (score - (-np.log10(sig_cutoff))) / max(cap - (-np.log10(sig_cutoff)), 1e-12)
    frac = float(np.clip(frac, 0.0, 1.0))
    return cmap(0.3 + frac * 0.7)


def _chunk_rows(rows: list, max_rows_per_col: int) -> list[list]:
    if not rows:
        return []
    n_cols = max(1, int(np.ceil(len(rows) / max_rows_per_col)))
    rows_per_col = int(np.ceil(len(rows) / n_cols))
    return [rows[i : i + rows_per_col] for i in range(0, len(rows), rows_per_col)]


def _infer_panel_name(value: float) -> str:
    return "Aging Up-Rescued" if value >= 0 else "Aging Down-Rescued"


def plot_staircase_rescue_heatmap(
    data: pd.DataFrame,
    *,
    row: str,
    column: str,
    aging_lfc: str,
    rescue_lfc: str,
    panel: str | None = None,
    aging_padj: str | None = None,
    rescue_padj: str | None = None,
    count: str | None = None,
    row_count: str | None = None,
    row_order: Sequence | None = None,
    column_order: Sequence | None = None,
    panel_order: Sequence | None = None,
    row_label_map: Mapping | None = None,
    column_label_map: Mapping | None = None,
    panel_config: Mapping | None = None,
    title: str | None = None,
    xlabel: str = "Context",
    ylabel: str = "Gene",
    count_label: str | None = None,
    right_count_label: str = "n Context",
    max_rows_per_col: int = 40,
    cell_in: float | None = None,
    cell_min: float = 0.14,
    cell_max: float = 0.22,
    min_width: float = DOUBLE_COL_IN,
    max_width: float = DOUBLE_COL_IN,
    top_margin: float = 0.86,
    bottom_margin: float = 0.72,
    left_margin: float = 0.86,
    right_margin: float = 0.55,
    panel_gap: float = 0.62,
    subpanel_gap: float = 0.34,
    xtick_rotation: float = 45,
    padj_sig_cutoff: float = 0.05,
    padj_score_cap: float = 10.0,
    show_legend: bool = True,
    show_cell_count: bool = True,
    show_right_counts: bool = True,
    figsize: tuple[float, float] | None = None,
):
    """Draw a staircase rescue heatmap with paired block arrows in each cell.

    This template follows the cross-tissue rescue figure style: genes are split
    into rescued-direction panels, each cell has a baseline, the left arrow
    encodes aging log2FC, and the right arrow starts at the aging tip to encode
    treatment/rescue log2FC. Arrow length follows ``log1p(abs(log2FC))`` and
    color saturation follows ``-log10(adjusted P)``.
    """
    if data.empty:
        raise ValueError("data is empty")

    df = data.copy()
    if panel is None:
        df["_panel"] = df[aging_lfc].map(_infer_panel_name)
        panel = "_panel"
    panels = _ordered(df[panel], panel_order)
    cols = _ordered(df[column], column_order)

    cfg = {k: dict(v) for k, v in STAIRCASE_PANEL_CONFIG.items()}
    if panel_config:
        for key, val in panel_config.items():
            cfg.setdefault(key, {}).update(val)

    if row_order is None:
        row_order_by_panel = {}
        for p in panels:
            sub = df[df[panel] == p]
            counts = sub.groupby(row)[column].nunique()
            row_order_by_panel[p] = counts.sort_values(ascending=False).index.tolist()
    else:
        row_rank = {v: i for i, v in enumerate(row_order)}
        row_order_by_panel = {}
        for p in panels:
            present = pd.Index(df.loc[df[panel] == p, row]).drop_duplicates().tolist()
            row_order_by_panel[p] = sorted(present, key=lambda v: row_rank.get(v, len(row_rank)))

    chunks_by_panel = {p: _chunk_rows(row_order_by_panel[p], max_rows_per_col) for p in panels}
    panels = [p for p in panels if chunks_by_panel.get(p)]
    if not panels:
        raise ValueError("no rows to plot")

    if cell_in is None:
        n_cell_cols = sum(len(chunks_by_panel[p]) * len(cols) for p in panels)
        gap_w = sum(max(0, len(chunks_by_panel[p]) - 1) * subpanel_gap for p in panels)
        gap_w += max(0, len(panels) - 1) * panel_gap
        avail = max_width - left_margin - right_margin - gap_w
        cell_in = float(np.clip(avail / max(n_cell_cols, 1), cell_min, cell_max))
    max_rows = max(len(chunk) for p in panels for chunk in chunks_by_panel[p])
    plot_w = sum(len(chunks_by_panel[p]) * len(cols) * cell_in for p in panels)
    plot_w += sum(max(0, len(chunks_by_panel[p]) - 1) * subpanel_gap for p in panels)
    plot_w += max(0, len(panels) - 1) * panel_gap
    fig_w = left_margin + plot_w + right_margin
    fig_w = max(min_width, min(max_width, fig_w))
    fig_h = top_margin + max_rows * cell_in + bottom_margin
    if figsize is not None:
        fig_w, fig_h = figsize

    vals = pd.concat([df[aging_lfc].abs(), df[rescue_lfc].abs()], ignore_index=True)
    vals = np.log1p(vals[np.isfinite(vals)])
    lfc_scale = float(np.nanpercentile(vals, 95)) if len(vals) else 1.0
    if not np.isfinite(lfc_scale) or lfc_scale <= 0:
        lfc_scale = 1.0

    row_counts = {}
    if row_count:
        row_counts = df.groupby(row)[row_count].max().to_dict()
    elif count:
        row_counts = df.groupby(row)[count].sum().to_dict()

    fig = plt.figure(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")
    if title:
        fig.text(0.5, 1 - 0.03 / fig_h, title, ha="center", va="top", fontsize=DEFAULT_FONT_SIZE)

    side_pad = max(0.0, (fig_w - (left_margin + plot_w + right_margin)) / 2)
    x_cursor = left_margin + side_pad
    y_top = fig_h - top_margin
    axes = []

    def draw_axis(ax, p, genes, *, show_title: bool, show_right_label: bool):
        pcfg = cfg.get(p, {})
        up_cfg = cfg.get("Aging Up-Rescued", STAIRCASE_PANEL_CONFIG["Aging Up-Rescued"])
        down_cfg = cfg.get("Aging Down-Rescued", STAIRCASE_PANEL_CONFIG["Aging Down-Rescued"])
        if p == "Aging Up-Rescued":
            aging_up_color = pcfg.get("aging_color", up_cfg["aging_color"])
            rescue_down_color = pcfg.get("rescue_color", up_cfg["rescue_color"])
            aging_down_color = down_cfg["aging_color"]
            rescue_up_color = down_cfg["rescue_color"]
        elif p == "Aging Down-Rescued":
            aging_down_color = pcfg.get("aging_color", down_cfg["aging_color"])
            rescue_up_color = pcfg.get("rescue_color", down_cfg["rescue_color"])
            aging_up_color = up_cfg["aging_color"]
            rescue_down_color = up_cfg["rescue_color"]
        else:
            aging_up_color = up_cfg["aging_color"]
            aging_down_color = down_cfg["aging_color"]
            rescue_down_color = up_cfg["rescue_color"]
            rescue_up_color = down_cfg["rescue_color"]
        cmap_aging_up = LinearSegmentedColormap.from_list("aging_up_rescue", ["#cfcfcf", aging_up_color])
        cmap_aging_down = LinearSegmentedColormap.from_list("aging_down_rescue", ["#cfcfcf", aging_down_color])
        cmap_rescue_down = LinearSegmentedColormap.from_list("rescue_down", ["#cfcfcf", rescue_down_color])
        cmap_rescue_up = LinearSegmentedColormap.from_list("rescue_up", ["#cfcfcf", rescue_up_color])
        rescue_bg_up = _tint(rescue_up_color, alpha=0.30)
        rescue_bg_down = _tint(rescue_down_color, alpha=0.30)

        n_rows = len(genes)
        n_cols = len(cols)
        ax.set_xlim(0, n_cols)
        ax.set_ylim(n_rows, 0)
        gap = 0.16
        gap2 = gap / 2
        margin_edge = 0.16
        overshoot = 0.18
        max_len = 1.0 - 2 * margin_edge - overshoot
        baseline_color = "#aaaaaa"
        missing_color = "#e0e0e0"
        present_color = "#f2f2f2"

        lookup = {
            (rec[row], rec[column]): rec
            for _, rec in df[(df[panel] == p) & (df[row].isin(genes))].iterrows()
        }
        for j, gene in enumerate(genes):
            for i, col in enumerate(cols):
                rec = lookup.get((gene, col))
                if rec is None:
                    bg = missing_color
                    has_data = False
                else:
                    has_data = True
                    rescue = float(rec[rescue_lfc])
                    aging = float(rec[aging_lfc])
                    is_rescue = np.sign(aging) != np.sign(rescue) and aging != 0 and rescue != 0
                    if is_rescue:
                        bg = rescue_bg_up if rescue > 0 else rescue_bg_down
                    else:
                        bg = present_color
                ax.add_patch(
                    FancyBboxPatch(
                        (i + gap2, j + gap2),
                        1 - gap,
                        1 - gap,
                        boxstyle="round,pad=0,rounding_size=0.06",
                        facecolor=bg,
                        edgecolor="none",
                        zorder=1,
                        mutation_aspect=1,
                    )
                )
                if has_data:
                    aging = float(rec[aging_lfc])
                    base_offset = 1 - margin_edge - overshoot if aging > 0 else margin_edge + overshoot
                    base_y = j + base_offset
                    ax.plot(
                        [i + gap2 + 0.04, i + 1 - gap2 - 0.04],
                        [base_y, base_y],
                        color=baseline_color,
                        linewidth=0.7,
                        zorder=2,
                        solid_capstyle="round",
                        alpha=0.7,
                    )

        for j, gene in enumerate(genes):
            for i, col in enumerate(cols):
                rec = lookup.get((gene, col))
                if rec is None:
                    continue
                aging = float(rec[aging_lfc])
                rescue = float(rec[rescue_lfc])
                base_offset = 1 - margin_edge - overshoot if aging > 0 else margin_edge + overshoot
                base_y = j + base_offset
                cell_top = j + margin_edge
                cell_bottom = j + 1 - margin_edge

                padj_a = float(rec[aging_padj]) if aging_padj else None
                padj_r = float(rec[rescue_padj]) if rescue_padj else None
                aging_dir = -1 if aging > 0 else 1
                rescue_dir = -1 if rescue > 0 else 1
                aging_exp = np.log1p(abs(aging))
                rescue_exp = np.log1p(abs(rescue))
                aging_len_want = min(aging_exp / lfc_scale, 1.0) * max_len
                room_a = (base_y - cell_top) if aging_dir == -1 else (cell_bottom - base_y)
                aging_len = min(aging_len_want, max(0.0, room_a))
                aging_tip = base_y + aging_dir * aging_len

                aging_cmap = cmap_aging_up if aging > 0 else cmap_aging_down
                rescue_cmap = cmap_rescue_up if rescue > 0 else cmap_rescue_down
                aging_color = _padj_color(padj_a, aging_cmap, sig_cutoff=padj_sig_cutoff, cap=padj_score_cap)
                rescue_color = _padj_color(padj_r, rescue_cmap, sig_cutoff=padj_sig_cutoff, cap=padj_score_cap)
                _block_arrow(ax, i + 0.30, base_y, aging_tip, aging_color)

                rescue_len_want = min(rescue_exp / lfc_scale, 1.0) * max_len
                room_r = (aging_tip - cell_top) if rescue_dir == -1 else (cell_bottom - aging_tip)
                rescue_len = min(rescue_len_want, max(0.0, room_r))
                rescue_tip = aging_tip + rescue_dir * rescue_len
                _block_arrow(ax, i + 0.70, aging_tip, rescue_tip, rescue_color)
                if rescue_len_want > rescue_len + 1e-6:
                    cap_y = j + gap2 + 0.045 if rescue_dir == -1 else j + 1 - gap2 - 0.045
                    ax.plot(
                        [i + 0.70 - 0.17, i + 0.70 + 0.17],
                        [cap_y, cap_y],
                        color=rescue_color,
                        linewidth=0.7,
                        solid_capstyle="round",
                        zorder=6,
                    )

                if show_cell_count and count and pd.notna(rec[count]):
                    label = str(int(rec[count])) if float(rec[count]).is_integer() else str(rec[count])
                    ax.text(i + 0.11, j + 0.16, label, ha="left", va="top", fontsize=5.8, color="#333333")

        col_labels = [str(column_label_map.get(c, c) if column_label_map else c) for c in cols]
        row_labels = [str(row_label_map.get(g, g) if row_label_map else g) for g in genes]
        ax.set_xticks(np.arange(n_cols) + 0.5)
        ax.set_xticklabels(col_labels, rotation=xtick_rotation, ha="right", fontsize=DEFAULT_FONT_SIZE)
        ax.set_yticks(np.arange(n_rows) + 0.5)
        ax.set_yticklabels(row_labels, fontsize=DEFAULT_FONT_SIZE, fontstyle="italic")
        ax.tick_params(length=0, pad=1.0)
        if show_title:
            ax.set_title(pcfg.get("title", str(p)), fontsize=DEFAULT_FONT_SIZE, pad=3, fontweight="bold")
        for spine in ax.spines.values():
            spine.set_visible(False)
        if show_right_counts and row_counts:
            for j, gene in enumerate(genes):
                val = row_counts.get(gene)
                if val is None or pd.isna(val):
                    continue
                label = str(int(val)) if float(val).is_integer() else str(val)
                ax.text(n_cols + 0.20, j + 0.5, label, ha="left", va="center", fontsize=DEFAULT_FONT_SIZE - 1)
            if show_right_label:
                ax.text(
                    n_cols + 0.82,
                    n_rows / 2,
                    right_count_label,
                    ha="center",
                    va="center",
                    rotation=90,
                    fontsize=DEFAULT_FONT_SIZE - 1,
                    color="#333333",
                    clip_on=False,
                )

    for p_idx, p in enumerate(panels):
        chunks = chunks_by_panel[p]
        for c_idx, genes in enumerate(chunks):
            ax_w = len(cols) * cell_in
            ax_h = len(genes) * cell_in
            ax = fig.add_axes([x_cursor / fig_w, (y_top - ax_h) / fig_h, ax_w / fig_w, ax_h / fig_h])
            draw_axis(
                ax,
                p,
                genes,
                show_title=c_idx == 0,
                show_right_label=(p_idx == len(panels) - 1 and c_idx == len(chunks) - 1),
            )
            axes.append(ax)
            x_cursor += ax_w
            if c_idx != len(chunks) - 1:
                x_cursor += subpanel_gap
        if p_idx != len(panels) - 1:
            x_cursor += panel_gap

    if ylabel:
        fig.text(0.015, (bottom_margin + max_rows * cell_in / 2) / fig_h, ylabel, rotation=90,
                 ha="left", va="center", fontsize=DEFAULT_FONT_SIZE)
    if xlabel:
        fig.text(0.5, 0.08 / fig_h, xlabel, ha="center", va="bottom", fontsize=DEFAULT_FONT_SIZE)
    if count_label:
        fig.text(0.5, 0.27 / fig_h, count_label, ha="center", va="bottom", fontsize=DEFAULT_FONT_SIZE - 1,
                 color="#333333")

    if show_legend:
        up_cfg = STAIRCASE_PANEL_CONFIG["Aging Up-Rescued"]
        down_cfg = STAIRCASE_PANEL_CONFIG["Aging Down-Rescued"]
        arrow_handles = [
            _ArrowProxy(up_cfg["aging_color"], +1),
            _ArrowProxy(up_cfg["rescue_color"], -1),
            _ArrowProxy(down_cfg["aging_color"], -1),
            _ArrowProxy(down_cfg["rescue_color"], +1),
            mpatches.Patch(facecolor="#cfcfcf", edgecolor="none"),
        ]
        arrow_labels = ["Aging up", "Rescue (down)", "Aging down", "Rescue (up)", "Not sig"]
        bg_handles = [
            mpatches.Patch(facecolor=_tint(down_cfg["rescue_color"], alpha=0.30), edgecolor="#c8c8c8", linewidth=0.6),
            mpatches.Patch(facecolor=_tint(up_cfg["rescue_color"], alpha=0.30), edgecolor="#c8c8c8", linewidth=0.6),
            mpatches.Patch(facecolor="#f2f2f2", edgecolor="#c8c8c8", linewidth=0.6),
            mpatches.Patch(facecolor="#e0e0e0", edgecolor="#c8c8c8", linewidth=0.6),
        ]
        bg_labels = ["Rescue cell (up)", "Rescue cell (down)", "Non-rescue", "No DE data"]
        leg_y = 1 - 0.05 / fig_h
        leg1 = fig.legend(
            handles=arrow_handles,
            labels=arrow_labels,
            title="Arrows - length: |lfc|   saturation: -log10(padj)",
            loc="upper left",
            bbox_to_anchor=(left_margin / fig_w, leg_y),
            frameon=False,
            handlelength=1.0,
            handletextpad=0.4,
            ncol=3,
            columnspacing=0.9,
            handler_map={_ArrowProxy: _ArrowHandler()},
            title_fontproperties={"weight": "bold"},
        )
        leg1._legend_box.align = "left"
        leg2 = fig.legend(
            handles=bg_handles,
            labels=bg_labels,
            title="Cell background",
            loc="upper left",
            bbox_to_anchor=(0.66, leg_y),
            frameon=False,
            handlelength=1.0,
            handletextpad=0.4,
            ncol=2,
            columnspacing=0.9,
            title_fontproperties={"weight": "bold"},
        )
        leg2._legend_box.align = "left"

    return fig, axes


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


__all__ = ["plot_rescue_arrow_heatmap", "plot_staircase_rescue_heatmap"]
