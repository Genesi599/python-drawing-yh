#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
热图 + 物种扇形叠加
横=组织（按本图基因集合内出现物种数降序）
纵=基因（上下调各按出现物种数降序取前20，上调在上）
"""

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import matplotlib.patches as mpatches
from matplotlib import cm, colors
import seaborn as sns


def plot_heatmap_fan_species(
        in_csv,
        out_dir=None,
        top_per_group=20,
        font_scale=1.2,
        r_base=0.4,
        figsize=(16, 9),
        dpi_png=600,
        species_colors=None,
        gene_type='all',
):
    """
    绘制热图+物种扇形叠加图

    Parameters
    ----------
    in_csv : str or Path
        输入CSV文件路径
    out_dir : str or Path, optional
        输出目录，默认为in_csv的parent/figure
    top_per_group : int
        上下调各取前N个基因，默认20
    font_scale : float
        字体缩放因子，默认1.2
    r_base : float
        扇形半径，默认0.4
    figsize : tuple
        图表大小，默认(16, 9)
    dpi_png : int
        PNG输出DPI，默认600
    species_colors : ndarray, optional
        自定义物种色板(N, 4)，RGBA格式
    gene_type : str
        基因类型，可选 'up'(上调), 'down'(下调), 'all'(全部)，默认'all'
    """
    in_csv = Path(in_csv)
    if out_dir is None:
        out_dir = in_csv.parent / "figure"
    else:
        out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    # 1. 读表
    df = pd.read_csv(in_csv)
    all_species = sorted(df['species'].unique())

    # 2. 物种色板
    if species_colors is None:
        species_colors = np.array([
            [1.00, 0.65, 0.00, 1.0],
            [0.40, 0.85, 0.20, 1.0],
            [0.95, 0.90, 0.10, 1.0],
            [0.75, 0.20, 0.75, 1.0],
            [0.00, 0.70, 0.70, 1.0],
            [0.90, 0.25, 0.55, 1.0]
        ])
    sp_colors = species_colors[:len(all_species)]


    # 3. 挑基因：上下调各自按「出现物种数」降序取前20
    gene_cnt = df.groupby('names').size().sort_values(ascending=False)
    # 上调基因：logfoldchanges>0
    up_genes = gene_cnt.loc[df[df['logfoldchanges'] > 0]['names'].unique()].nlargest(top_per_group).index
    # 下调基因：logfoldchanges<0
    down_genes = gene_cnt.loc[df[df['logfoldchanges'] < 0]['names'].unique()].nlargest(top_per_group).index

    if gene_type == 'up':
        genes = list(up_genes)
    elif gene_type == 'down':
        genes = list(down_genes)
    else:  # 'all'
        genes = list(up_genes) + list(down_genes)
    # 组织：基因在不同组织的出现频次
    sub_df = df[df['names'].isin(genes)]
    tissue_cnt = sub_df.groupby('tissue')['names'].nunique().sort_values(ascending=False)
    org_order = tissue_cnt.index

    # 4. 组织×基因原始 logFC 矩阵（按新顺序重排）
    mat = (sub_df.pivot_table(index='tissue', columns='names', values='logfoldchanges')
           .fillna(0).reindex(index=org_order, columns=genes))

    # 5. 指数变换矩阵
    exp_mat = np.sign(mat) * np.log1p(np.abs(mat))
    norm = colors.TwoSlopeNorm(vmin=-np.abs(exp_mat.values).max(),
                               vcenter=0,
                               vmax=np.abs(exp_mat.values).max())
    cmap_fill = cm.coolwarm

    # 6. 画布
    n_row, n_col = exp_mat.shape
    fig, ax = plt.subplots(figsize=figsize)

    # 7. 热图（转置 → 横=组织，纵=基因）
    sns.heatmap(
        exp_mat.T,
        cmap=cmap_fill,
        norm=norm,
        linewidths=0.3,
        linecolor='#f2f2f2',
        xticklabels=True,
        yticklabels=True,
        cbar=False,
        ax=ax
    )

    # 8. 叠加固定大小扇形
    for i, tissue in enumerate(exp_mat.T.columns):  # 横轴
        for j, gene in enumerate(exp_mat.T.index):  # 纵轴
            cx = i + 0.5
            cy = j + 0.5
            sub = df[(df['tissue'] == tissue) & (df['names'] == gene)]
            if sub.empty:
                continue
            for k, sp in enumerate(all_species):
                v = sub.loc[sub['species'] == sp, 'logfoldchanges']
                if v.empty or v.isna().all():
                    continue
                val = v.iloc[0]
                r = r_base
                exp_val = np.sign(val) * np.log1p(np.abs(val))
                theta1 = k * 60
                theta2 = (k + 1) * 60
                wedge = mpatches.Wedge(
                    (cx, cy), r, theta1, theta2,
                    facecolor=cmap_fill(norm(exp_val)),
                    edgecolor=sp_colors[k],
                    linewidth=1.2
                )
                ax.add_patch(wedge)

    # 9. 坐标轴
    ax.set_xlabel('Tissue', fontsize=14 * font_scale)
    ax.set_ylabel('Gene', fontsize=14 * font_scale)
    ax.set_title('Cross-Species Aging B Cell DEG', fontsize=16 * font_scale)
    plt.xticks(rotation=45, ha='right', fontsize=10 * font_scale)
    plt.yticks(rotation=0, fontsize=10 * font_scale)

    # 10. 物种图例（带边框扇形）
    legend_handles = []
    for k, sp in enumerate(all_species):
        wedge = mpatches.Wedge((0, 0), 1, 0, 60,
                               facecolor='none',
                               edgecolor=sp_colors[k],
                               linewidth=2)
        legend_handles.append(wedge)
    ax.legend(legend_handles, all_species,
              loc='upper left', bbox_to_anchor=(1.02, 1),
              frameon=False, fontsize=12 * font_scale,
              handlelength=1.5, handletextpad=0.5)

    # 11. colorbar（右下方）
    cbar_ax = fig.add_axes([0.81, 0.15, 0.03, 0.3])
    sm = cm.ScalarMappable(cmap=cmap_fill, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
    cbar.ax.tick_params(labelsize=12 * font_scale)
    cbar.set_label('exp-logFC', fontsize=12 * font_scale)

    plt.tight_layout()
    # 12. 保存
    for ext in ['png', 'pdf']:
        out_path = out_dir / f'heatmap_fan_speciesTop.{ext}'
        fig.savefig(out_path, dpi=dpi_png if ext == 'png' else None, bbox_inches='tight')
        print(f'  ✅ {out_path.name} 完成')
    plt.close()


def main():
    # --------------------------
    in_csv = Path(r"D:\Projects\B_Cell_Aging\B_cell_tissue_DEG_results_filtered.csv")
    plot_heatmap_fan_species(
        in_csv=in_csv,
        top_per_group=20,
        font_scale=1.5,
        r_base=0.4,
        gene_type='up',
    )


if __name__ == '__main__':
    main()
