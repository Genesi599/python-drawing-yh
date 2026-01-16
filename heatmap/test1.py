#!/usr/bin/env python
# coding: utf-8
"""
无聚版：行方向 z-score + age 排序 + overall_trend 顺序 + 最小 3 行高
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import zscore

EXPR_CSV    = r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\abundance_sample_x_protein.csv"
META_CSV    = r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\sample_meta.csv"
PATTERN_CSV = r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\Pattern_Analysis\all_proteins_pattern.csv"
OUT_DIR     = Path("outputs")
OUT_DIR.mkdir(exist_ok=True)

FC_CUT      = 0
TOP_ANNOT   = 20
BLANK_ROWS  = 2
FIG_W, FIG_H = 9, 16

colors = ["#0042ae", "#f7f7f7", "#fd1a41"]   # 蓝→白→红
custom_seismic = LinearSegmentedColormap.from_list("custom_seismic", colors, N=256)
# custom_seismic = sns.color_palette("coolwarm", as_cmap=True)


def main():
    expr_df = pd.read_csv(EXPR_CSV, index_col=0)
    meta    = pd.read_csv(META_CSV, index_col="sample")
    df_pat  = pd.read_csv(PATTERN_CSV)
    if 'gene' in df_pat.columns:
        df_pat = df_pat.set_index('gene')
    df_pat.index = df_pat.index.str.replace(r'-\d+$', '', regex=True)
    expr_df.columns = expr_df.columns.str.replace(r'-\d+$', '', regex=True)
    df_pat = df_pat[df_pat['pattern'] != 'Non-significant']

    # 按 overall_trend 排序：Up→Down→Mixed
    trend_dict = df_pat.groupby('pattern')['overall_trend'].first().to_dict()
    def sort_key(p):
        return {'Up': 0, 'Down': 1, 'Mixed': 2}.get(trend_dict.get(p, 'Mixed'), 3), p
    patterns = sorted(df_pat['pattern'].unique(), key=sort_key)

    condition_map = meta["condition"].str.capitalize()
    age_series    = pd.to_numeric(meta["age"], errors="coerce")
    common_smp = expr_df.index.intersection(condition_map.index)
    expr_df    = expr_df.loc[common_smp]
    condition_map = condition_map.loc[common_smp]
    age_series    = age_series.reindex(common_smp)
    expr_df = expr_df.T

    print("=== 各 pattern 匹配蛋白数 ===")
    for pat in patterns:
        logfc = pd.to_numeric(df_pat.loc[df_pat['pattern'] == pat, "log2FC_Young_vs_Old"],
                              errors="coerce").reindex(expr_df.index)
        print(f"{pat:30s} {logfc.notna().sum()}")

    # 拼大矩阵 & z-score（无聚类）
    big_blocks, y_ticks, y_labels = [], [], []
    annotate_genes, annotate_y = [], []
    base = 0

    for pat in patterns:
        logfc = pd.to_numeric(df_pat.loc[df_pat['pattern'] == pat, "log2FC_Young_vs_Old"],
                              errors="coerce").reindex(expr_df.index)
        genes_pat = logfc.dropna().index
        if genes_pat.empty:
            continue
        mat_pat = expr_df.loc[genes_pat]
        # 行方向 z-score
        mat_pat = pd.DataFrame(
            zscore(mat_pat, axis=1, nan_policy='omit'),
            index=mat_pat.index,
            columns=mat_pat.columns
        )

        # ****** 只按 logFC 排序，不聚类 ******
        up_genes   = logfc[logfc >  FC_CUT].sort_values(ascending=False).index
        down_genes = logfc[logfc < -FC_CUT].sort_values(ascending=True).index
        order = up_genes.tolist() + down_genes.tolist() if len(up_genes) + len(down_genes) > 0 else genes_pat.tolist()
        mat_pat = mat_pat.loc[order]

        # 右侧标注
        top_up   = logfc[logfc >  FC_CUT].nlargest(TOP_ANNOT).index
        top_down = logfc[logfc < -FC_CUT].nsmallest(TOP_ANNOT).index
        show = [g for g in top_up if g in mat_pat.index] + \
               [g for g in top_down if g in mat_pat.index]
        for g in show:
            annotate_genes.append(g)
            annotate_y.append(base + mat_pat.index.get_loc(g) + 0.5)

        # 拼块
        big_blocks.append(mat_pat)
        y_labels.extend(mat_pat.index)
        y_ticks.append(base + len(mat_pat) / 2)
        base += len(mat_pat)

        # 空白分隔
        blank = pd.DataFrame(np.nan, index=[f"BLANK_{i}" for i in range(BLANK_ROWS)],
                             columns=mat_pat.columns)
        big_blocks.append(blank)
        y_labels.extend(blank.index)
        base += BLANK_ROWS

    if not big_blocks:
        print("❌ 没有可画数据")
        return

    big_mat = pd.concat(big_blocks, axis=0)

    # 按 age 升序排列样本
    sample_order = big_mat.columns.sort_values(key=lambda x: age_series[x])
    big_mat = big_mat[sample_order]

    # 颜色条
    lut = {"Young": "#1f77b4", "Middle": "#ff7f0e", "Old": "#d62728"}
    sample_colors = sample_order.map(condition_map.map(lut))

    # 画 clustermap（无聚类）
    vmin = float(np.nanmin(big_mat.values))
    vmax = float(np.nanmax(big_mat.values))
    norm = plt.Normalize(vmin=vmin, vmax=vmax)

    g = sns.clustermap(
        big_mat,
        method='ward',
        metric='euclidean',
        row_cluster=False,   # 关键：不再聚类
        col_cluster=False,
        cmap=custom_seismic,
        norm=norm,
        xticklabels=True,
        yticklabels=False,
        figsize=(FIG_W, FIG_H),
        cbar_pos=None,
        dendrogram_ratio=(0.15, 0.15),
        linewidths=0,     # ← 关键：把网格线宽度设为 0
    )

    # 样本颜色条
    bar_h, y0 = 0.02, 1.01
    for i, s in enumerate(sample_order):
        g.ax_heatmap.add_patch(Rectangle((i, y0), 1, bar_h,
                                         facecolor=sample_colors[i],
                                         transform=g.ax_heatmap.get_xaxis_transform(),
                                         clip_on=False))


    # 右侧标注
    if annotate_genes:
        import adjustText as at
        texts = []
        for gene, y in zip(annotate_genes, annotate_y):
            t = g.ax_heatmap.text(big_mat.shape[1] + 0.5, y, gene,
                                  fontsize=9, fontweight='bold',
                                  va='center', ha='left', clip_on=False)
            texts.append(t)
        at.adjust_text(texts, ax=g.ax_heatmap,
                       max_iter=1000, max_move=(0, 5),
                       expand_text=(0, 5), force_text=(0, 500),
                       only_move={'points': 'y', 'text': 'y'},
                       arrowprops=None, autoalign=False)
        # 引线
        def left_edge(t):
            bbox = t.get_window_extent(renderer=g.fig.canvas.get_renderer()) \
                    .transformed(g.ax_heatmap.transData.inverted())
            return bbox.x0, t.get_position()[1]
        gene2y0 = dict(zip(annotate_genes, annotate_y))
        for t in texts:
            x_end, y_end = left_edge(t)
            x_start = big_mat.shape[1]
            y_start = gene2y0[t.get_text()]
            g.ax_heatmap.annotate('', xy=(x_end, y_end),
                                  xytext=(x_start, y_start),
                                  arrowprops=dict(arrowstyle='-', lw=1,
                                                  color='black', shrinkA=0),
                                  annotation_clip=False)

    # 保存
    out_file = OUT_DIR / "protein_patterns_heatmap.png"
    plt.savefig(out_file, dpi=600, bbox_inches='tight')
    plt.savefig(out_file.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()
    print("✅ 无聚类版 z-score + age 排序 + overall_trend 顺序 完成：", out_file)


if __name__ == "__main__":
    main()