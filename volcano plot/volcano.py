from math import log10

import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math
from adjustText import adjust_text

def create_volcano_plot(
        input_file, output_file='Volcano_plot.png',
        x_threshold=0.5, y_threshold=-np.log10(0.05),   # 注意用 np.log10
        lfc_col='log2FoldChange', p_col='padj',
        id_col='GeneName'):

    # 读数据
    df = pd.read_csv(input_file)

    # 只要三列，去掉多余
    vol = df[[id_col, lfc_col, p_col]].copy()
    vol['y'] = -np.log10(vol[p_col])

    # 分组：保证每一行都被覆盖
    vol['group'] = 'black'                       # 1. normal
    up_mask = (vol[lfc_col] >= x_threshold) & (vol['y'] >= y_threshold)
    dn_mask = (vol[lfc_col] <= -x_threshold) & (vol['y'] >= y_threshold)
    vol.loc[up_mask, 'group'] = 'tab:red'        # 2. up
    vol.loc[dn_mask, 'group'] = 'tab:blue'       # 3. down

    # 排序（按离原点距离）
    vol['dist'] = vol[lfc_col]**2 + vol['y']**2
    vol = vol.sort_values('dist', ascending=False)

    # 画图
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(vol[lfc_col], vol['y'],
               s=3, c=vol['group'], alpha=0.8)

    # 阈值线
    ax.axvline(-x_threshold, ls='--', color='grey', lw=1)
    ax.axvline(x_threshold, ls='--', color='grey', lw=1)
    ax.axhline(y_threshold, ls='--', color='grey', lw=1)

    # 标签
    ax.set_xlabel('Log2 Fold Change', fontweight='bold')
    ax.set_ylabel('-Log10(adj. p-value)', fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # 给显著下调基因加文字（取前 10 个）
    dn_sig = vol[dn_mask].head(10)
    texts = [ax.text(r[lfc_col]-0.05, r['y'], r[id_col],
                     fontsize=10, style='italic', weight='bold')
             for _, r in dn_sig.iterrows()]
    adjust_text(texts, arrowprops=dict(arrowstyle='->', lw=0.5, color='green'))

    dn_sig = vol[up_mask].head(10)
    texts = [ax.text(r[lfc_col]-0.05, r['y'], r[id_col],
                     fontsize=10, style='italic', weight='bold')
             for _, r in dn_sig.iterrows()]
    adjust_text(texts, arrowprops=dict(arrowstyle='->', lw=0.5, color='green'))

    ax.set_xlabel('Log2 Fold Change', fontweight='bold', fontsize=12)
    ax.set_ylabel('-Log10(adj. p-value)', fontweight='bold', fontsize=12)

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close()

# # 示例使用
# create_volcano_plot(
#     input_file='day7_Gene_Summary.xlsx',
#     output_file='Volcano_plot.png',
#     x_threshold=1,
#     y_threshold=2,
#     pos_lfc_col='pos|lfc',
#     pos_pval_col='pos|p-value',
#     neg_lfc_col='neg|lfc',
#     neg_pval_col='neg|p-value',
#     id_col='id'
# )