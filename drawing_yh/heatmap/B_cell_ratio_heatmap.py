#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : B_cell_ratio_heatmap.py
@Date    : 2026/2/7 17:13
@Author  : yh109
"""
# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : B_cell_ratio_heatmap.py
@Date    : 2026/2/7 16:44
@Author  : yh109
"""
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from collections import defaultdict
import matplotlib as mpl
import seaborn as sns
from pathlib import Path

# 配置字体
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.use14corefonts'] = False
mpl.rcParams['font.family'] = 'Arial'

# 数据路径
data_dir = r"D:\Projects\B_Cell_Aging\B_cell_ratio\correlation_results"
output_dir = r"D:\Projects\B_Cell_Aging\B_cell_ratio\figure\bubble_plot"
os.makedirs(output_dir, exist_ok=True)

# 读取所有文件
data_files = sorted(glob.glob(os.path.join(data_dir, "*.csv")))

corr_data = {}
for fpath in data_files:
    filename = Path(fpath).stem
    species = filename.split('_')[0].capitalize()
    df = pd.read_csv(fpath, index_col=0)
    df['-log10_pval'] = -np.log10(df['Spearman_P_value'])
    corr_data[species] = df


def sort_genes_by_freq_and_strength(corr_data, direction='up', pval_threshold=0.05):
    counts = defaultdict(int)
    sum_abs_r = defaultdict(float)
    for species, df in corr_data.items():
        if direction == 'up':
            subset = df[df['Spearman_Correlation'] > 0]
        else:
            subset = df[df['Spearman_Correlation'] < 0]
        for tissue in subset.index:
            pval = subset.loc[tissue, 'Spearman_P_value']
            if pval < pval_threshold:
                r_val = subset.loc[tissue, 'Spearman_Correlation']
                counts[tissue] += 1
                sum_abs_r[tissue] += abs(r_val)
    avg_abs_r = {t: sum_abs_r[t] / counts[t] for t in counts}
    sorted_tissues = sorted(counts.keys(), key=lambda t: (-counts[t], -avg_abs_r[t]))
    return sorted_tissues


datasets_list = sorted(list(corr_data.keys()))
up_tissues = sort_genes_by_freq_and_strength(corr_data, direction='up', pval_threshold=0.05)[:50]

# 保存选中的tissue
up_tissues_df = pd.DataFrame({
    'Tissue': up_tissues,
    'Direction': ['Up'] * len(up_tissues)
})
up_tissues_df.to_csv(os.path.join(output_dir, 'selected_tissues.csv'), index=False)

# 构建热图数据矩阵和注释矩阵
heatmap_data = np.zeros((len(datasets_list), len(up_tissues)))
annot_data = [['' for _ in range(len(up_tissues))] for _ in range(len(datasets_list))]

for x_idx, tissue in enumerate(up_tissues):
    for y_idx, species in enumerate(datasets_list):
        if tissue in corr_data[species].index:
            r_val = corr_data[species].loc[tissue, 'Spearman_Correlation']
            pval = corr_data[species].loc[tissue, 'Spearman_P_value']
            heatmap_data[y_idx, x_idx] = r_val

            if pval < 0.001:
                annot_data[y_idx][x_idx] = '***'
            elif pval < 0.01:
                annot_data[y_idx][x_idx] = '**'
            elif pval < 0.05:
                annot_data[y_idx][x_idx] = '*'

# 计算配色范围
vmax = max(abs(np.nanmin(heatmap_data)), abs(np.nanmax(heatmap_data)))
vmin = -vmax

# 计算图大小
n_rows, n_cols = heatmap_data.shape
cell_size = 0.3
fig_width = n_cols * cell_size + 10
fig_height = n_rows * cell_size + 1.5
font_scale = 1.5

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# 绘制热图
sns.heatmap(heatmap_data, annot=False, cmap='RdBu_r', center=0, vmin=vmin, vmax=vmax,
            xticklabels=up_tissues, yticklabels=datasets_list,
            cbar_kws={'label': 'Spearman r'}, linewidths=0.5, linecolor='gray',
            ax=ax, square=False)

# 添加显著性标记
for i in range(len(datasets_list)):
    for j in range(len(up_tissues)):
        label = annot_data[i][j]
        if label:
            ax.text(j + 0.5, i + 0.5, label,
                    ha='center', va='center',
                    fontsize=int(15 * font_scale), fontweight='bold', color='black',
                    zorder=4)

# 设置刻度标签
ax.set_xticklabels(up_tissues, fontsize=int(11 * font_scale), rotation=45, ha='right')
ax.set_yticklabels(datasets_list, fontsize=int(12 * font_scale), rotation=0)
# ax.set_title('Up-regulation', fontsize=int(14 * font_scale), fontweight='bold')

# 调整颜色条
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=int(10 * font_scale))
cbar.set_label('Spearman r', fontsize=int(12 * font_scale))

plt.subplots_adjust(bottom=0.25, right=0.9)

# 保存
plt.savefig(os.path.join(output_dir, 'heatmap_bcell.png'), dpi=600, pad_inches=0.1, bbox_inches='tight')
plt.savefig(os.path.join(output_dir, 'heatmap_bcell.pdf'),
            dpi=300, bbox_inches='tight', pad_inches=0.1,
            metadata={'Creator': None, 'Producer': None})
plt.close()

print(f"Heatmap saved to {output_dir}")
print(f"Up-regulated tissues: {len(up_tissues)}")
print(f"Species included: {datasets_list}")