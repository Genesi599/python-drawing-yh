#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : age_dotmap_by_pattern.py
@Date    : 2026/1/22 11:09
@Author  : yh109
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from statsmodels.nonparametric.smoothers_lowess import lowess



def draw_pattern_pointplot(expr_df, df_pat, patterns, age_series, enrich_dict,
                           feature_meta_df, id_col, name_col, logfc_col,
                           sample_meta, figsize=(14, 16), base_fontsize=20, fc_cut=0,
                           split_up_down=False, out_dir=None, is_up_regulation=True,
                           enrich_show_n=2):
    from matplotlib.gridspec import GridSpec
    from matplotlib.colors import to_rgb, to_hex


    warm_pool = [
        "#FF8888", "#FF9944", "#FFAA66", "#FF88BB",
        "#BB88FF", "#FF7799", "#FFAA44", "#FF8899",
    ]
    cool_pool = [
        "#6699FF", "#66BBFF", "#66FFAA", "#88FFFF",
        "#99CCFF", "#77DDFF", "#99CCFF", "#88BBFF",
    ]

    pattern_trend = df_pat.groupby('pattern')['overall_trend'].first().to_dict()
    protein_to_gene = dict(zip(feature_meta_df[id_col], feature_meta_df[name_col]))

    pattern_info = []
    fc_trend = []
    sample_order = expr_df.columns.sort_values(key=lambda x: age_series[x])

    for pat_idx, pat in enumerate(patterns):
        logfc = pd.to_numeric(df_pat.loc[df_pat['pattern'] == pat, logfc_col],
                              errors="coerce").reindex(expr_df.index)
        genes_pat = logfc.dropna().index
        if genes_pat.empty:
            continue

        up_genes = logfc[logfc > fc_cut].sort_values(ascending=False).index
        down_genes = logfc[logfc < -fc_cut].sort_values(ascending=True).index

        if split_up_down:
            genes_to_use = up_genes if is_up_regulation else down_genes
        else:
            genes_to_use = pd.Index(list(up_genes) + list(down_genes))

        if genes_to_use.empty:
            continue

        order = up_genes.tolist() + down_genes.tolist() if len(up_genes) + len(down_genes) > 0 else genes_pat.tolist()

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
            'genes_pat': order,
            'trend': trend,
            'fc_color': fc_color
        })
        fc_trend.append((pat, fc_color))

    if not pattern_info:
        print("❌ 没有可画数据")
        return

    n_patterns = len(pattern_info)
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(n_patterns, 100, figure=fig, hspace=0.4, wspace=0.3)

    plot_axes = []
    for p_idx, info in enumerate(pattern_info):
        ax = fig.add_subplot(gs[p_idx, :70])
        plot_axes.append(ax)

        # 根据 pattern 名称判断性别
        pat_lower = info['pattern'].lower()
        if 'female' in pat_lower:
            sex_filter = 'F'
        elif 'male' in pat_lower:
            sex_filter = 'M'
        else:
            sex_filter = None

        # 过滤样本
        if sex_filter:
            valid_samples = sample_meta[sample_meta['sex'] == sex_filter].index
            plot_sample_order = [s for s in sample_order if s in valid_samples]
        else:
            plot_sample_order = sample_order

        genes = info['genes_pat']
        gene_colors = plt.cm.tab20(np.linspace(0, 1, len(genes)))

        for g_idx, gene in enumerate(genes):
            expr = expr_df.loc[gene, plot_sample_order].values
            ages = age_series[plot_sample_order].values

            # 标准化表达值
            expr_norm = (expr - expr.mean()) / expr.std()

            ax.scatter(ages, expr_norm, alpha=0.6, s=80, color=gene_colors[g_idx],
                       label=protein_to_gene.get(gene, gene)[:20])

            if len(ages) > 1:
                lowess_result = lowess(expr_norm, ages, frac=0.8)
                fit_x = lowess_result[:, 0]
                fit_y = lowess_result[:, 1]
                ax.plot(fit_x, fit_y, color=gene_colors[g_idx], linewidth=2, alpha=0.8)

        ax.set_ylabel('Expression', fontsize=base_fontsize - 2)
        ax.set_xlabel('Age', fontsize=base_fontsize - 2)
        ax.set_title(f"{info['pattern']}", fontsize=base_fontsize, fontweight='bold')
        ax.tick_params(labelsize=base_fontsize - 4)
        ax.grid(alpha=0.3)

        # 为这个 pattern 创建对应的颜色条
        ax_fc = fig.add_subplot(gs[p_idx, 70:75])
        n_genes = len(info['genes_pat'])
        ax_fc.add_patch(Rectangle((0, 0), 0.3, n_genes, facecolor=info['fc_color'], edgecolor='black', linewidth=1))
        ax_fc.set_xlim(0, 1)
        ax_fc.set_ylim(0, n_genes)
        ax_fc.axis('off')

    FC_BAR_WIDTH = 0.05
    ENRICH_BAR_WIDTH = 0.6
    ENRICH_BAR_GAP = 0.02

    hm_left = 0.1
    hm_width = 0.55
    fc_left = hm_left + hm_width
    enrich_left = fc_left + FC_BAR_WIDTH + ENRICH_BAR_GAP

    # 富集图
    all_terms, all_pvals, all_genes = [], [], []
    bar_positions, bar_colors = [], []
    pos = 0
    pattern_bar_map = {}

    for info in reversed(pattern_info):
        enrich_list = info['enrich_data']
        pat_start = pos

        if not enrich_list:
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

        pos += 0.5

    ax_bar = fig.add_axes([enrich_left, 0.1, ENRICH_BAR_WIDTH, 0.8])
    ax_bar.barh(bar_positions, all_pvals, color=bar_colors, height=0.6, alpha=0.75)

    def darken_color(hex_color, factor=0.6):
        rgb = to_rgb(hex_color)
        return to_hex([c * factor for c in rgb])

    GENE_FS = base_fontsize - 1
    TERM_FS = base_fontsize
    for idx, (term, genes) in enumerate(zip(all_terms, all_genes)):
        x_pos = all_pvals[idx] * 0.02 if all_pvals[idx] > 0 else 0.01
        term_color = darken_color(bar_colors[idx])

        if genes == '':
            ax_bar.text(x_pos, bar_positions[idx], term,
                        va='center', ha='left',
                        fontsize=TERM_FS, color=term_color, fontweight='bold')
        else:
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

    suffix = "upregulation" if is_up_regulation else "downregulation"
    out_file = out_dir / f"protein_patterns_pointplot_{suffix}.png"
    plt.savefig(out_file, dpi=600, bbox_inches='tight')
    plt.savefig(out_file.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()
    print("✅ 点图完成：", out_file)
