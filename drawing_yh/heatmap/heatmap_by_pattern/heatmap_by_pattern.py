#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : heatmap_by_pattern.py
@Date    : 2026/1/22 10:29
@Author  : yh109
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import zscore
from matplotlib.colors import to_rgb, to_hex
from matplotlib.colors import LinearSegmentedColormap

from drawing_yh import DEFAULT_FONT_SIZE, save_fig


def draw_heatmap(expr_df, df_pat, patterns, condition_map, age_series, enrich_dict,
                 feature_meta_df, id_col, name_col, logfc_col,
                 blank_rows=2, figsize=(10, 16), base_fontsize=DEFAULT_FONT_SIZE, fc_cut=0,
                 split_up_down=False, out_dir=None, is_up_regulation=True,
                 enrich_show_n=2):
    LINE_COLOR = 'black'
    LINE_WIDTH = 2
    LINE_ALPHA = 1

    warm_pool = [
        "#FF8888", "#FF9944", "#FFAA66", "#FF88BB",
        "#BB88FF", "#FF7799", "#FFAA44", "#FF8899",
    ]
    cool_pool = [
        "#6699FF", "#66BBFF", "#66FFAA", "#88FFFF",
        "#99CCFF", "#77DDFF", "#99CCFF", "#88BBFF",
    ]

    colors = ["#0042ae", "#f7f7f7", "#fd1a41"]
    custom_seismic = LinearSegmentedColormap.from_list("custom_seismic", colors, N=256)

    pattern_trend = df_pat.groupby('pattern')['overall_trend'].first().to_dict()
    protein_to_gene = dict(zip(feature_meta_df[id_col], feature_meta_df[name_col]))

    big_blocks, pattern_info = [], []
    fc_trend = []
    base = 0

    for pat_idx, pat in enumerate(patterns):
        logfc = pd.to_numeric(df_pat.loc[df_pat['pattern'] == pat, logfc_col],
                              errors="coerce").reindex(expr_df.index)
        genes_pat = logfc.dropna().index
        if genes_pat.empty:
            continue

        up_genes = logfc[logfc > fc_cut].sort_values(ascending=False).index
        down_genes = logfc[logfc < -fc_cut].sort_values(ascending=True).index

        if split_up_down:
            # 只保留上调或下调基因（需要外层循环控制）
            genes_to_use = up_genes if is_up_regulation else down_genes
        else:
            genes_to_use = pd.Index(list(up_genes) + list(down_genes))

        if genes_to_use.empty:
            continue
        mat_pat = expr_df.loc[genes_pat]
        mat_pat = pd.DataFrame(
            zscore(mat_pat, axis=1, nan_policy='omit'),
            index=mat_pat.index,
            columns=mat_pat.columns
        )

        up_genes = logfc[logfc > fc_cut].sort_values(ascending=False).index
        down_genes = logfc[logfc < -fc_cut].sort_values(ascending=True).index
        order = up_genes.tolist() + down_genes.tolist() if len(up_genes) + len(down_genes) > 0 else genes_pat.tolist()
        mat_pat = mat_pat.loc[order]

        pat_start = base
        pat_end = base + len(mat_pat)
        pat_center_y = (pat_start + pat_end) / 2

        enrich_data = enrich_dict.get(pat, [])
        trend = pattern_trend.get(pat, 'Mixed')

        if trend.lower() == 'up':
            fc_color = warm_pool[pat_idx % len(warm_pool)]
        elif trend.lower() == 'down':
            fc_color = cool_pool[pat_idx % len(cool_pool)]
        else:
            fc_color = '#CCCCCC'

        pattern_info.append({
            'pattern': pat,
            'enrich_data': enrich_data,
            'genes_pat': order,  # 添加这行
            'y_start': pat_start,
            'y_end': pat_end,
            'y_center': pat_center_y,
            'trend': trend,
            'fc_color': fc_color
        })

        big_blocks.append(mat_pat)
        fc_trend.extend([fc_color] * len(order))
        base += len(mat_pat)

        blank = pd.DataFrame(np.nan, index=[f"BLANK_{i}" for i in range(blank_rows)],
                             columns=mat_pat.columns)
        big_blocks.append(blank)
        fc_trend.extend(['#FFFFFF'] * blank_rows)
        base += blank_rows

    if not big_blocks:
        print("❌ 没有可画数据")
        return

    big_mat = pd.concat(big_blocks, axis=0)
    sample_order = big_mat.columns.sort_values(key=lambda x: age_series[x])
    big_mat = big_mat[sample_order]

    lut = {"Young": "#1f77b4", "Middle": "#ff7f0e", "Old": "#d62728"}
    sample_colors = sample_order.map(condition_map.map(lut))

    vmin = float(np.nanmin(big_mat.values))
    vmax = float(np.nanmax(big_mat.values))
    norm = plt.Normalize(vmin=vmin, vmax=vmax)

    g = sns.clustermap(
        big_mat,
        method='ward',
        metric='euclidean',
        row_cluster=False,
        col_cluster=False,
        cmap=custom_seismic,
        norm=norm,
        xticklabels=True,
        yticklabels=False,
        figsize=figsize,
        cbar_pos=None,
        dendrogram_ratio=(0.15, 0.15),
        linewidths=0,
    )

    fig = g.fig
    hm_pos = g.ax_heatmap.get_position()

    # 添加FC竖条 (与热图y轴完全对齐)
    hm_ylim = g.ax_heatmap.get_ylim()  # clustermap反转: (n-0.5, -0.5)
    FC_BAR_WIDTH = 0.08
    ax_fc = fig.add_axes([hm_pos.x1, hm_pos.y0, FC_BAR_WIDTH, hm_pos.height])

    for i, color in enumerate(reversed(fc_trend)):
        ax_fc.add_patch(Rectangle((0, i - 0.5), 0.3, 1, facecolor=color, edgecolor='none'))

    ax_fc.set_xlim(0, 1)
    # 使用与热图相同的反转ylim，使得 reversed 后的 y=0(最后一行) 在顶部
    ax_fc.set_ylim(len(fc_trend) - 0.5, -0.5)
    ax_fc.axis('off')

    g.ax_heatmap.tick_params(axis='x', labelsize=base_fontsize)

    bar_h, y0 = 0.02, 1.01
    for i, s in enumerate(sample_order):
        g.ax_heatmap.add_patch(Rectangle((i, y0), 1, bar_h,
                                         facecolor=sample_colors[i],
                                         transform=g.ax_heatmap.get_xaxis_transform(),
                                         clip_on=False))

    age_group_positions = {}
    for i, s in enumerate(sample_order):
        age_group = condition_map[s]
        if age_group not in age_group_positions:
            age_group_positions[age_group] = []
        age_group_positions[age_group].append(i)

    for age_group, indices in age_group_positions.items():
        center_x = (min(indices) + max(indices)) / 2 + 0.5
        g.ax_heatmap.text(center_x, y0 + bar_h + 0.01, age_group,
                          fontsize=base_fontsize, ha='center', va='bottom', fontweight='bold',
                          transform=g.ax_heatmap.get_xaxis_transform(),
                          clip_on=False)

    for info in pattern_info:
        g.ax_heatmap.text(-0.5, info['y_center'], info['pattern'],
                          fontsize=base_fontsize, ha='right', va='center',
                          rotation=0, clip_on=False, fontweight='bold')

    # 准备富集数据
    all_terms, all_pvals, all_genes = [], [], []
    bar_positions, pattern_bounds, bar_colors = [], [], []
    pos = 0
    bar_idx = 0
    pattern_bar_map = {}

    for info in reversed(pattern_info):
        enrich_list = info['enrich_data']
        pat_start = pos

        if not enrich_list:
            # 没有term，显示前6个基因，分成两行，每行3个用分号分隔
            top_genes = [protein_to_gene.get(g, g) for g in info['genes_pat'][:6]]
            top_genes = [g.capitalize() for g in top_genes]

            line1 = '; '.join(top_genes[:3])
            line2 = '; '.join(top_genes[3:])

            for line in [line1, line2]:
                all_terms.append(line)
                all_pvals.append(0)
                all_genes.append('')
                bar_positions.append(pos)
                bar_colors.append(info['fc_color'])
                pos += 1

            pattern_bar_map[info['pattern']] = (pos - 2, pos - 1)
        else:
            df_enrich = pd.DataFrame(enrich_list)
            # 二次筛选：只展示前 3 个
            show_n = enrich_show_n
            df_enrich = df_enrich.sort_values('P-value').head(show_n)
            n_bars = len(df_enrich)
            pat_start_bar = pos
            pat_end_bar = pos + n_bars - 1

            for _, row in df_enrich.iterrows():
                all_terms.append(row['Term'])
                all_pvals.append(row['neg_log10_p'])
                all_genes.append(row['Genes'])
                bar_positions.append(pos)
                bar_colors.append(info['fc_color'])
                pos += 1

            pattern_bar_map[info['pattern']] = (pat_start_bar, pat_end_bar)

        pat_end = pos
        pattern_bounds.append((info['pattern'], pat_start, pat_end))
        pos += 0.5
        bar_idx += 1

    # 添加柱状图axes到右侧
    ENRICH_BAR_WIDTH = 1
    ENRICH_BAR_GAP = 0
    ax_bar = fig.add_axes([hm_pos.x1 + FC_BAR_WIDTH + ENRICH_BAR_GAP, hm_pos.y0, ENRICH_BAR_WIDTH, hm_pos.height])

    ax_bar.barh(bar_positions, all_pvals, color=bar_colors, height=0.6, alpha=0.75)

    def darken_color(hex_color, factor=0.6):
        rgb = to_rgb(hex_color)
        return to_hex([c * factor for c in rgb])
    GENE_FS = base_fontsize - 1
    TERM_FS = base_fontsize
    for idx, (term, genes) in enumerate(zip(all_terms, all_genes)):
        x_pos = all_pvals[idx] * 0.02 if all_pvals[idx] > 0 else 0.01
        term_color = darken_color(bar_colors[idx])

        if genes == '':  # 基因列表行
            ax_bar.text(x_pos, bar_positions[idx], term,
                        va='center', ha='left',
                        fontsize=TERM_FS, color=term_color, fontweight='bold')
        else:  # term行
            term_clean = term.split(' (')[0]
            ax_bar.text(x_pos, bar_positions[idx], term_clean,
                        va='center', ha='left',
                        fontsize=TERM_FS, color='black', fontweight='bold')
            gene_color = darken_color(bar_colors[idx])
            ax_bar.text(x_pos, bar_positions[idx] - 0.25, f'({genes})',
                        va='top', ha='left',
                        fontsize=TERM_FS - 1, color=gene_color, fontweight='bold')

    ax_bar.set_xlabel('-Log10(p)', fontsize=base_fontsize, fontweight='bold')
    ax_bar.spines[['top', 'right', 'left']].set_visible(False)
    ax_bar.tick_params(left=False, labelleft=False, bottom=True, labelsize=base_fontsize - 2)

    ax_bar.set_ylim(-1, pos - 1)

    # 添加pattern标签
    for pat, pat_start, pat_end in pattern_bounds:
        ax_bar.axhline(y=pat_start - 0.8, color=LINE_COLOR, linewidth=LINE_WIDTH, alpha=LINE_ALPHA)
        # pat_mid = (pat_start + pat_end) / 2 - 0.5
        # ax_bar.text(-max(all_pvals) * 0.1, pat_mid, pat,
        #             fontsize=BASE_FONTSIZE, fontweight='bold', ha='right', va='center')

    # 右侧连线
    for info in pattern_info:
        bar_range = pattern_bar_map.get(info['pattern'])
        if bar_range is None:
            continue

        pat_start, pat_end = bar_range
        line_y_bar = pat_start - 0.8

        line_y_fc = len(fc_trend) - info['y_end']

        x1, y1 = ax_fc.transData.transform((0.3, line_y_fc))
        x2, y2 = ax_bar.transData.transform((0, line_y_bar))

        x1, y1 = fig.transFigure.inverted().transform((x1, y1))
        x2, y2 = fig.transFigure.inverted().transform((x2, y2))

        fig.lines.append(plt.Line2D(
            [x1, x2],
            [y1, y2],
            transform=fig.transFigure,
            color=LINE_COLOR, linewidth=LINE_WIDTH, alpha=LINE_ALPHA
        ))

    suffix = "upregulation" if is_up_regulation else "downregulation"
    out_file = out_dir / f"protein_patterns_heatmap_{suffix}.png"
    save_fig(fig, out_file, also=('.pdf', '.svg'))
    plt.close()
    print("✅ 热力图完成：", out_file)
