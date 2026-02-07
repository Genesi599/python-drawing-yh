#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : B_cell_ratio_bubble.py
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
    species = filename.split('_')[0]
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

# 计算气泡大小范围（基于-log10(p-value)）
plotted_logp = []
for species, df in corr_data.items():
    for tissue in up_tissues:
        if tissue in df.index:
            plotted_logp.append(df.loc[tissue, '-log10_pval'])

logp_min = np.min(plotted_logp) if plotted_logp else 0
logp_max = np.max(plotted_logp) if plotted_logp else 1
logp_range = logp_max - logp_min if logp_max > logp_min else 1

target_min_size = 100
target_max_size = 600
logp_scale = (target_max_size - target_min_size) / logp_range if logp_range > 0 else 1

# 绘图
fig, ax = plt.subplots(1, 1, figsize=(16, 6))

norm = Normalize(vmin=-1, vmax=1)
cmap = plt.cm.RdBu_r
font_scale = 2

for x_idx, tissue in enumerate(up_tissues):
    df_list = {sp: corr_data[sp] for sp in datasets_list}
    for y_idx, species in enumerate(datasets_list):
        df = df_list[species]
        if tissue in df.index:
            spearman_r = df.loc[tissue, 'Spearman_Correlation']
            logp = df.loc[tissue, '-log10_pval']
            pval = df.loc[tissue, 'Spearman_P_value']

            ax.scatter(x_idx, y_idx, s=target_min_size + (logp - logp_min) * logp_scale,
                       c=[spearman_r], cmap=cmap, norm=norm, edgecolor='black', linewidth=0.5, alpha=0.7, zorder=3)

            if pval < 0.001:
                sig_marker = '***'
            elif pval < 0.01:
                sig_marker = '**'
            elif pval < 0.05:
                sig_marker = '*'
            else:
                sig_marker = ''

            if sig_marker:
                ax.text(x_idx, y_idx, sig_marker, ha='center', va='center',
                        fontsize=int(7 * font_scale), fontweight='bold', color='black', zorder=4, alpha=0.7)

ax.set_xticks(range(len(up_tissues)))
ax.set_xticklabels(up_tissues, fontsize=int(11 * font_scale), rotation=45, ha='right')
ax.set_yticks(range(len(datasets_list)))
ax.set_yticklabels(datasets_list, fontsize=int(12 * font_scale))
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_title('Up-regulation', fontsize=int(14 * font_scale), fontweight='bold')
ax.margins(y=0.3)

# 颜色条（右上）
cbar_ax = fig.add_axes([0.92, 0.55, 0.02, 0.35])
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
cbar.set_label('Spearman r', fontsize=int(12 * font_scale))
cbar.set_ticks([-1, 0, 1])
cbar.ax.tick_params(labelsize=int(10 * font_scale))

# 气泡大小图例（右下）
logp_med = np.median(plotted_logp) if plotted_logp else (logp_min + logp_max) / 2
size_values = [logp_min, logp_med, logp_max]
size_labels = [f'{logp_min:.2f}', f'{logp_med:.2f}', f'{logp_max:.2f}']
size_handles = [plt.scatter([], [], s=target_min_size + (val - logp_min) * logp_scale,
                            color='gray', alpha=0.7, edgecolor='black', linewidth=0.5) for val in size_values]
fig.legend(size_handles, size_labels, loc='upper left', bbox_to_anchor=(0.88, 0.45),
           ncol=1, frameon=True, fontsize=int(10 * font_scale), title='-log₁₀(p-value)',
           title_fontsize=int(11 * font_scale))

plt.subplots_adjust(bottom=0.25, right=0.85)

# 保存
plt.savefig(os.path.join(output_dir, 'bubble_plot_bcell.png'), dpi=600, pad_inches=0.1, bbox_inches='tight')
plt.savefig(os.path.join(output_dir, 'bubble_plot_bcell.pdf'),
            dpi=300,
            bbox_inches='tight',
            pad_inches=0.1,
            metadata={'Creator': None, 'Producer': None})
plt.close()

print(f"Plot saved to {output_dir}")
print(f"Up-regulated tissues: {len(up_tissues)}")
print(f"Species included: {datasets_list}")
