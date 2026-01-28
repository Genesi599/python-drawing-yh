#!/usr/bin/env python
# coding: utf-8
"""
无聚版：行方向 z-score + age 排序 + overall_trend 顺序 + 右侧富集柱状图
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import zscore
import gseapy as gp
from tenacity import retry, stop_after_attempt, wait_fixed
import os
from matplotlib.colors import to_rgb, to_hex


@retry(stop=stop_after_attempt(3), wait=wait_fixed(2))
def _enrich_genes(gene_list):
    if len(gene_list) < 5:
        return None
    return gp.enrichr(gene_list=list(gene_list),
                      gene_sets='GO_Biological_Process_2025',
                      organism='Human', outdir=None)


base_path = Path(r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis")
EXPR_CSV = base_path / 'data/abundance_sample_x_protein.csv'
META_CSV = base_path / 'data/sample_meta.csv'
PATTERN_CSV = base_path / 'data/Pattern_Analysis/all_proteins_pattern.csv'
ENRICH_OUT_DIR = base_path / 'data/Pattern_Analysis/enrich_results'

feature_meta = pd.read_csv(base_path / 'data/feature_meta.csv')
protein_to_gene = dict(zip(feature_meta['protein'], feature_meta['Gene name']))
OUT_DIR = base_path / 'figure'
OUT_DIR.mkdir(exist_ok=True)

ENRICH_N_RESULTS = 20
ENRICH_N_DISPLAY = 2

FC_CUT = 0
BLANK_ROWS = 2
FIG_W, FIG_H = 10, 16
BASE_FONTSIZE = 20

colors = ["#0042ae", "#f7f7f7", "#fd1a41"]
custom_seismic = LinearSegmentedColormap.from_list("custom_seismic", colors, N=256)


def enrich_and_save(expr_df, df_pat, patterns, out_dir=ENRICH_OUT_DIR):
    os.makedirs(out_dir, exist_ok=True)
    enrich_dict = {}

    for i, pat in enumerate(patterns, 1):
        print(f"\n[{i}/{len(patterns)}] 处理 pattern: {pat}")

        logfc = pd.to_numeric(df_pat.loc[df_pat['pattern'] == pat, "log2FC_Young_vs_Old"],
                              errors="coerce").reindex(expr_df.index)
        genes_pat = logfc.dropna().index
        print(f"  ├─ 该pattern总蛋白数: {len(genes_pat)}")

        if genes_pat.empty:
            print(f"  └─ ❌ 跳过: 无匹配基因")
            enrich_dict[pat] = []
            continue

        all_genes = expr_df.loc[genes_pat].index.tolist()
        gene_names = [protein_to_gene.get(g, g) for g in all_genes]
        print(f"  ├─ 转换后基因名: {gene_names}")

        print(f"  ├─ 调用Enrichr API...")
        try:
            enrich_res = _enrich_genes(gene_names)
            if enrich_res is None:
                print(f"  ├─ ⚠ API返回None")
                enrich_dict[pat] = []
                continue

            df_enrich = enrich_res.results
            if df_enrich.empty:
                print(f"  ├─ ⚠ 返回空结果")
                enrich_dict[pat] = []
                continue

            df_enrich = df_enrich.sort_values('P-value').head(ENRICH_N_RESULTS)
            df_enrich['neg_log10_p'] = -np.log10(df_enrich['P-value'])

            print(f"  ├─ 取前2个显著项:")
            for idx, row in df_enrich.iterrows():
                print(f"     - {row['Term']} (p={row['P-value']:.2e})")

            enrich_dict[pat] = df_enrich[['Term', 'P-value', 'Genes', 'neg_log10_p']].to_dict('records')

            fpath = f"{out_dir}/{pat.replace(' ', '_')}_enrich.csv"
            df_enrich.to_csv(fpath, index=False)
            print(f"  └─ ✅ 已保存: {fpath}")

        except Exception as e:
            print(f"  ├─ ❌ 错误: {type(e).__name__}: {str(e)}")
            enrich_dict[pat] = []

    print(f"\n=== 富集分析完成 ===")
    return enrich_dict


def draw_heatmap(expr_df, df_pat, patterns, condition_map, age_series, enrich_dict):
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

    pattern_trend = df_pat.groupby('pattern')['overall_trend'].first().to_dict()

    big_blocks, pattern_info = [], []
    fc_trend = []
    base = 0

    for pat_idx, pat in enumerate(patterns):
        logfc = pd.to_numeric(df_pat.loc[df_pat['pattern'] == pat, "log2FC_Young_vs_Old"],
                              errors="coerce").reindex(expr_df.index)
        genes_pat = logfc.dropna().index
        if genes_pat.empty:
            continue
        mat_pat = expr_df.loc[genes_pat]
        mat_pat = pd.DataFrame(
            zscore(mat_pat, axis=1, nan_policy='omit'),
            index=mat_pat.index,
            columns=mat_pat.columns
        )

        up_genes = logfc[logfc > FC_CUT].sort_values(ascending=False).index
        down_genes = logfc[logfc < -FC_CUT].sort_values(ascending=True).index
        order = up_genes.tolist() + down_genes.tolist() if len(up_genes) + len(down_genes) > 0 else genes_pat.tolist()
        mat_pat = mat_pat.loc[order]

        pat_start = base
        pat_end = base + len(mat_pat)
        pat_center_y = (pat_start + pat_end) / 2

        enrich_data = enrich_dict.get(pat, [])
        trend = pattern_trend.get(pat, 'Mixed')

        if trend == 'Up':
            fc_color = warm_pool[pat_idx % len(warm_pool)]
        elif trend == 'Down':
            fc_color = cool_pool[pat_idx % len(cool_pool)]
        else:
            fc_color = '#CCCCCC'

        pattern_info.append({
            'pattern': pat,
            'enrich_data': enrich_data,
            'y_start': pat_start,
            'y_end': pat_end,
            'y_center': pat_center_y,
            'trend': trend,
            'fc_color': fc_color
        })

        big_blocks.append(mat_pat)
        fc_trend.extend([fc_color] * len(order))
        base += len(mat_pat)

        blank = pd.DataFrame(np.nan, index=[f"BLANK_{i}" for i in range(BLANK_ROWS)],
                             columns=mat_pat.columns)
        big_blocks.append(blank)
        fc_trend.extend(['#FFFFFF'] * BLANK_ROWS)
        base += BLANK_ROWS

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
        figsize=(FIG_W, FIG_H),
        cbar_pos=None,
        dendrogram_ratio=(0.15, 0.15),
        linewidths=0,
    )

    fig = g.fig
    hm_pos = g.ax_heatmap.get_position()

    # 添加FC竖条
    FC_BAR_WIDTH = 0.08
    ax_fc = fig.add_axes([hm_pos.x1, hm_pos.y0, FC_BAR_WIDTH, hm_pos.height])

    for i, color in enumerate(reversed(fc_trend)):
        ax_fc.add_patch(Rectangle((0, i), 0.3, 1, facecolor=color, edgecolor='none'))

    ax_fc.set_xlim(0, 1)
    ax_fc.set_ylim(-0.5, len(fc_trend) - 0.5)
    ax_fc.axis('off')

    g.ax_heatmap.tick_params(axis='x', labelsize=BASE_FONTSIZE)

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
                          fontsize=BASE_FONTSIZE, ha='center', va='bottom', fontweight='bold',
                          transform=g.ax_heatmap.get_xaxis_transform(),
                          clip_on=False)

    for info in pattern_info:
        g.ax_heatmap.text(-0.5, info['y_center'], info['pattern'],
                          fontsize=BASE_FONTSIZE, ha='right', va='center',
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
            pattern_bar_map[info['pattern']] = None
        else:
            df_enrich = pd.DataFrame(enrich_list[:ENRICH_N_DISPLAY])
            n_bars = len(df_enrich)
            pat_start_bar = pos
            pat_end_bar = pos + n_bars - 1

            for _, row in df_enrich.iterrows():
                all_terms.append(row['Term'])
                all_pvals.append(row['neg_log10_p'])
                gene_txt = ', '.join(row['Genes'].split(';')[:5]) if pd.notna(row['Genes']) else ''
                all_genes.append(gene_txt)
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
    TERM_FS = BASE_FONTSIZE
    GENE_FS = BASE_FONTSIZE - 1
    for idx, (term, genes) in enumerate(zip(all_terms, all_genes)):
        x_pos = all_pvals[idx] * 0.02
        gene_color = darken_color(bar_colors[idx])
        term_clean = term.split(' (')[0]
        ax_bar.text(x_pos, bar_positions[idx], term_clean,
                    va='center', ha='left',
                    fontsize=TERM_FS, color='black', fontweight='bold')

        ax_bar.text(x_pos, bar_positions[idx] - 0.25, f'({genes})',
                    va='top', ha='left',
                    fontsize=GENE_FS, color=gene_color, fontweight='bold', fontstyle='italic')

    ax_bar.set_xlabel('-Log10(p)', fontsize=BASE_FONTSIZE, fontweight='bold')
    ax_bar.spines[['top', 'right', 'left']].set_visible(False)
    ax_bar.tick_params(left=False, labelleft=False, bottom=True, labelsize=BASE_FONTSIZE - 2)

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

    out_file = OUT_DIR / "protein_patterns_heatmap.png"
    plt.savefig(out_file, dpi=600, bbox_inches='tight')
    plt.savefig(out_file.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()
    print("✅ 热力图完成：", out_file)


def load_or_enrich(expr_df, df_pat, patterns, out_dir=ENRICH_OUT_DIR):
    # 检查是否有完整的富集结果
    if os.path.exists(out_dir):
        csv_files = [f for f in os.listdir(out_dir) if f.endswith('_enrich.csv')]
        if len(csv_files) == len(patterns):
            print("✅ 检测到现存富集结果，直接读取")
            enrich_dict = {}
            for pat in patterns:
                fpath = os.path.join(out_dir, f"{pat.replace(' ', '_')}_enrich.csv")
                if os.path.exists(fpath):
                    df_enrich = pd.read_csv(fpath)
                    enrich_dict[pat] = df_enrich[['Term', 'P-value', 'Genes', 'neg_log10_p']].to_dict('records')
                else:
                    enrich_dict[pat] = []
            return enrich_dict

    print("❌ 未检测到完整富集结果，执行新的富集分析")
    return enrich_and_save(expr_df, df_pat, patterns, out_dir)



def main():
    expr_df = pd.read_csv(EXPR_CSV, index_col=0)
    meta = pd.read_csv(META_CSV, index_col="sample")
    df_pat = pd.read_csv(PATTERN_CSV)
    if 'gene' in df_pat.columns:
        df_pat = df_pat.set_index('gene')
    df_pat.index = df_pat.index.str.replace(r'-\d+$', '', regex=True)
    expr_df.columns = expr_df.columns.str.replace(r'-\d+$', '', regex=True)
    df_pat = df_pat[df_pat['pattern'] != 'Non-significant']

    trend_dict = df_pat.groupby('pattern')['overall_trend'].first().to_dict()

    def sort_key(p):
        return {'Up': 0, 'Down': 1, 'Mixed': 2}.get(trend_dict.get(p, 'Mixed'), 3), p

    patterns = sorted(df_pat['pattern'].unique(), key=sort_key)

    condition_map = meta["condition"].str.capitalize()
    age_series = pd.to_numeric(meta["age"], errors="coerce")
    common_smp = expr_df.index.intersection(condition_map.index)
    expr_df = expr_df.loc[common_smp]
    condition_map = condition_map.loc[common_smp]
    age_series = age_series.reindex(common_smp)
    expr_df = expr_df.T

    print("=== 各 pattern 匹配蛋白数 ===")
    for pat in patterns:
        logfc = pd.to_numeric(df_pat.loc[df_pat['pattern'] == pat, "log2FC_Young_vs_Old"],
                              errors="coerce").reindex(expr_df.index)
        print(f"{pat:30s} {logfc.notna().sum()}")

    # 第一步：富集分析
    print("\n=== 检查富集结果 ===")
    enrich_dict = load_or_enrich(expr_df, df_pat, patterns)

    # 第二步：绘制热力图+富集图
    print("\n=== 生成热力图 ===")
    draw_heatmap(expr_df, df_pat, patterns, condition_map, age_series, enrich_dict)


if __name__ == "__main__":
    main()
