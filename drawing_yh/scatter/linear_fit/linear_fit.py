#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : main.py
@Date    : 2026/2/4 14:22
@Author  : yh109
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import FancyArrowPatch
import matplotlib as mpl


# 1. 让 matplotlib 别把字体子集化
mpl.rcParams['pdf.fonttype'] = 42          # 已经写了，保留
mpl.rcParams['ps.fonttype'] = 42
# 关键：关闭子集化（>=3.6 有效）
mpl.rcParams['pdf.use14corefonts'] = False   # 不用 14 种核心字体
# 2. 指定一个系统里肯定有的字体，避免 mpl 走自带字体
mpl.rcParams['font.family'] = 'Arial'      # 或 'DejaVu Sans', 'SimHei' 等

def _choose_grid(n_panels, target_ratio=16/9):
    best = None
    for cols in range(1, min(n_panels, 8) + 1):
        rows = int(np.ceil(n_panels / cols))
        ratio = (cols * 1.0) / max(1, rows)
        score = abs(ratio - target_ratio)
        if best is None or score < best[0]:
            best = (score, rows, cols)
    return best[1], best[2]


def plot_correlation_barplot(cor_df, output_path, cell_type):
    plt.figure(figsize=(14, 6))
    colors = ["red" if p < 0.05 else "gray"
              for p in cor_df["p_value"]]
    plt.barh(cor_df["tissue_general"], cor_df["pearson_r"],
             color=colors)
    plt.xlabel("Pearson r", fontsize=12)
    plt.ylabel("Tissue", fontsize=12)
    plt.title(f"{cell_type} Ratio vs Age Correlation by Tissue "
              f"(red: p<0.05)", fontsize=13)
    plt.axvline(0, color="black", linewidth=0.8, linestyle="--")
    plt.tight_layout()
    plt.savefig(f"{output_path}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{output_path}.pdf", bbox_inches="tight")
    plt.close()


def plot_all_tissues(cor_df, donor_summary, cell_count_col,
                     cell_ratio_col, age_col, age_unit, xlim_margin,
                     cell_type, output_path, exclude_zero=True,
                     organ_palette=None, tissues_to_plot=None,
                     tissue_colors=None, font_scale=1.8,
                     panel_size=5.0, save_svg=True, dpi=600,
                     scatter_size=25, line_width=1.2, title_pad=2,
                     # Explicit per-element font sizes. When None, fall back
                     # to font_scale × legacy defaults (preserves old callers).
                     title_size=None, tick_size=None, stats_size=None,
                     # Axis-arrow geometry (inches).
                     arrow_len_in=None, arrow_start_in=0.08,
                     arrow_lw=0.8, arrow_head_scale=8,
                     arrow_label_size=None, arrow_label_offset_in=0.14):
    if tissues_to_plot is not None:
        tissues = [t for t in tissues_to_plot if t in cor_df["tissue_general"].values]
        cor_df = cor_df[cor_df["tissue_general"].isin(tissues)]
    else:
        tissues = cor_df["tissue_general"].tolist()

    n_tissues = len(tissues)
    n_rows, n_cols = _choose_grid(n_tissues, target_ratio=8 / 8)

    tissues_to_plot = cor_df['tissue_general'].unique()
    data_for_range = donor_summary[donor_summary['tissue_general'].isin(tissues_to_plot)]
    if exclude_zero and cell_count_col:
        data_for_range = data_for_range[data_for_range[cell_count_col] > 0]
    x_min = data_for_range[age_col].min()
    x_max = data_for_range[age_col].max()

    if tissue_colors is not None:
        palette = [tissue_colors.get(t, (0.5, 0.5, 0.5)) for t in tissues]
    elif organ_palette is not None:
        palette = [organ_palette.get(t, (0.5, 0.5, 0.5)) for t in tissues]
    else:
        palette = sns.color_palette("husl", n_tissues)

    fig_width = n_cols * panel_size
    fig_height = n_rows * panel_size
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(fig_width, fig_height))
    if n_tissues == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for idx, tissue in enumerate(tissues):
        ax = axes[idx]
        row = cor_df[cor_df["tissue_general"] == tissue].iloc[0]
        rho = row["pearson_r"]
        p_val = row["p_value"]

        data = donor_summary[donor_summary["tissue_general"] == tissue]
        if exclude_zero:
            data = data[data[cell_count_col] > 0]

        if len(data) == 0:
            ax.axis("off")
            continue
        color = palette[idx]
        sns.regplot(data=data, x=age_col, y=cell_ratio_col,
                    ax=ax, color=color,
                    scatter_kws={"alpha": 0.6, "s": scatter_size,
                                 "edgecolor": "none"},
                    line_kws={"linewidth": line_width})
        ax.set_xlim(x_min - xlim_margin, x_max + xlim_margin)
        _title_size = title_size if title_size is not None else int(18 * font_scale)
        _tick_size  = tick_size  if tick_size  is not None else int(15 * font_scale)
        _stats_size = stats_size if stats_size is not None else int(13 * font_scale)
        ax.set_title(f"{tissue}", fontsize=_title_size, fontweight="bold",
                     pad=title_pad)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=_tick_size)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        ax.set_aspect("auto")
        ax.text(0.05, 0.95, f"r = {rho:.3f}\np = {p_val:.3f}",
                transform=ax.transAxes, fontsize=_stats_size,
                verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="white",
                          alpha=0.8, edgecolor="gray", linewidth=0.5))

    for idx in range(n_tissues, len(axes)):
        axes[idx].axis("off")

    # Axis arrows: default length = 30% of shortest panel dim (0.3 in for
    # panel_size=1.6), font falls back to legacy 30pt if unset.
    _arrow_len = arrow_len_in if arrow_len_in is not None else max(0.3, panel_size * 0.3)
    _arrow_font = arrow_label_size if arrow_label_size is not None else int(30 * font_scale)
    _add_axis_arrows(fig, fig_width, fig_height, age_unit, cell_type,
                     arrow_len_inch=_arrow_len,
                     arrow_start_inch=arrow_start_in,
                     line_width=arrow_lw,
                     mutation_scale=arrow_head_scale,
                     font_size=_arrow_font,
                     label_offset_inch=arrow_label_offset_in)

    # Tighter rect — axis-arrow space at lower-left was 0.025, 0.01 suffices
    # with the smaller panel_size; also bring panels closer together.
    plt.tight_layout(rect=[0.01, 0.01, 1, 1], w_pad=0.2, h_pad=0.2)
    # Each matplotlib backend only accepts a narrow metadata key set; passing
    # a PDF-shaped dict to SVG raises. Dispatch per-format.
    _common = dict(bbox_inches='tight', pad_inches=0.02)
    plt.savefig(f"{output_path}.png", dpi=dpi,
                metadata={'Software': None}, **_common)
    plt.savefig(f"{output_path}.pdf", dpi=dpi,
                metadata={'Creator': None, 'Producer': None}, **_common)
    if save_svg:
        # SVG must be 72 dpi so 1 pt of font renders as 1 px (so downstream
        # vector editors get correct text metrics).
        plt.savefig(f"{output_path}.svg", dpi=72,
                    metadata={'Creator': None}, **_common)
    plt.close()




def plot_significant_tissues(sig_cor, donor_summary, cell_count_col,
                             cell_ratio_col, age_col, age_unit,
                             xlim_margin, x_min, x_max, cell_type,
                             output_path, exclude_zero=True):
    n_sig = len(sig_cor)
    if n_sig == 0:
        return

    n_rows_sig, n_cols_sig = _choose_grid(n_sig, target_ratio=16/9)
    fig_width_sig = n_cols_sig * 5
    fig_height_sig = n_rows_sig * 5
    fig, axes = plt.subplots(n_rows_sig, n_cols_sig,
                             figsize=(fig_width_sig, fig_height_sig))
    if n_sig == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    sig_palette = sns.color_palette("husl", n_sig)

    for idx, (_, row) in enumerate(sig_cor.iterrows()):
        ax = axes[idx]
        tissue = row["tissue_general"]
        rho = row["pearson_r"]
        p_val = row["p_value"]

        data = donor_summary[
            donor_summary["tissue_general"] == tissue
            ]
        if exclude_zero:
            data = data[data[cell_count_col] > 0]

        color = sig_palette[idx]
        sns.regplot(data=data, x=age_col, y=cell_ratio_col,
                    ax=ax, color=color,
                    scatter_kws={"alpha": 0.6, "s": 200},
                    line_kws={"linewidth": 3})
        ax.set_xlim(x_min - xlim_margin, x_max + xlim_margin)
        ax.set_title(f"{tissue}", fontsize=18, fontweight="bold",
                     pad=10)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=15)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        ax.set_aspect("auto")
        ax.text(0.05, 0.95, f"r = {rho:.3f}\np = {p_val:.3f}",
                transform=ax.transAxes, fontsize=14,
                verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="white",
                          alpha=0.8, edgecolor="gray"))

    for idx in range(n_sig, len(axes)):
        axes[idx].axis("off")

    _add_axis_arrows(fig, fig_width_sig, fig_height_sig, age_unit,
                     cell_type)

    plt.tight_layout(rect=[0.03, 0.03, 1, 1])
    plt.savefig(f"{output_path}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{output_path}.pdf", bbox_inches="tight")
    plt.close()


def _add_axis_arrows(fig, fig_width, fig_height, age_unit, cell_type,
                     arrow_len_inch=1.5, arrow_start_inch=0.2,
                     line_width=4, mutation_scale=35, font_size=30,
                     label_offset_inch=0.3):
    """Lower-left axis-arrow overlay. All geometry in absolute inches so
    arrows look consistent regardless of ``figsize``."""
    arrow_len_x = arrow_len_inch / fig_width
    arrow_len_y = arrow_len_inch / fig_height
    arrow_x_start = arrow_start_inch / fig_width
    arrow_y_start = arrow_start_inch / fig_height

    arrow_x = FancyArrowPatch(
        (arrow_x_start, arrow_y_start),
        (arrow_x_start + arrow_len_x, arrow_y_start),
        transform=fig.transFigure,
        arrowstyle="->", lw=line_width, color="black",
        mutation_scale=mutation_scale)
    arrow_y = FancyArrowPatch(
        (arrow_x_start, arrow_y_start),
        (arrow_x_start, arrow_y_start + arrow_len_y),
        transform=fig.transFigure,
        arrowstyle="->", lw=line_width, color="black",
        mutation_scale=mutation_scale)
    fig.add_artist(arrow_x)
    fig.add_artist(arrow_y)

    # Labels anchored away from the corner so their bboxes never intersect.
    # X-label: left-aligned, starts slightly right of the y-arrow so the
    # rotated y-label never sits on top of its first letter.
    # Y-label: bottom-aligned, rises from slightly above the x-arrow so its
    # first letter (read bottom-up) never touches "Age" below.
    x_label_x = arrow_x_start + label_offset_inch / fig_width
    x_label_y = arrow_y_start - label_offset_inch / fig_height
    y_label_x = arrow_x_start - label_offset_inch / fig_width
    y_label_y = arrow_y_start + label_offset_inch / fig_height

    fig.text(x_label_x, x_label_y,
             f"Age ({age_unit})", ha="left", va="top",
             fontsize=font_size, fontweight="bold")
    fig.text(y_label_x, y_label_y,
             f"{cell_type} Ratio", ha="right", va="bottom",
             fontsize=font_size, fontweight="bold", rotation=90)

# ---------- 新增：器官-物种配色器 ----------
def _make_organ_species_palette(tissues):
    """
    tissues : list, 元素格式为 "Organ_species" 或任意字符串
    返回: dict {tissue: RGB 三元色}, 长度与 tissues 相同
    """
    import matplotlib.colors as mcolors

    # 1. 基础器官色系（Hue 固定）
    organ_hue = {
        "Liver": 0.00,   # 红
        "Lung": 0.45,    # 青绿
        "Kidney": 0.58,  # 蓝
    }

    # 2. 物种在 HSV 里的 (亮度, 饱和度) 微调
    species_vs = {
        "human": (1.0, 1.0),
        "mouse": (0.90, 0.85),
        "monkey": (1.10, 0.90),
        "rat": (0.80, 0.75),
    }

    palette = {}
    for t in tissues:
        parts = t.split("_", 1)
        organ = parts[0]
        species = parts[1] if len(parts) == 2 else "human"

        h = organ_hue.get(organ, 0.15)          # 默认橙
        v_adj, s_adj = species_vs.get(species, (1.0, 1.0))

        # 生成颜色
        rgb = mcolors.hsv_to_rgb((h, min(s_adj, 1.0), min(v_adj, 1.0)))
        palette[t] = rgb
    return palette