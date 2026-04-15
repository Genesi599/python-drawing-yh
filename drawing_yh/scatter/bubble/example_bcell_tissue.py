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
    species = Path(fpath).stem
    df = pd.read_csv(fpath, index_col=0)
    df['-log2_pval'] = -np.log2(df['Spearman_P_value'])
    corr_data[species] = df


def sort_genes_by_freq_and_strength(corr_data, direction='up'):
    counts = defaultdict(int)
    sum_abs_r = defaultdict(float)
    for species, df in corr_data.items():
        if direction == 'up':
            subset = df[df['Spearman_Correlation'] > 0]
        else:
            subset = df[df['Spearman_Correlation'] < 0]
        for tissue in subset.index:
            r_val = subset.loc[tissue, 'Spearman_Correlation']
            counts[tissue] += 1
            sum_abs_r[tissue] += abs(r_val)
    avg_abs_r = {t: sum_abs_r[t] / counts[t] for t in counts}
    sorted_tissues = sorted(counts.keys(), key=lambda t: (-counts[t], -avg_abs_r[t]))
    return sorted_tissues


datasets_list = sorted(list(corr_data.keys()))
up_tissues = sort_genes_by_freq_and_strength(corr_data, direction='up')[:50]
down_tissues = sort_genes_by_freq_and_strength(corr_data, direction='down')[:50]

# 保存选中的tissue
all_tissues = pd.DataFrame({
    'Tissue': list(up_tissues) + list(down_tissues),
    'Direction': ['Up'] * len(up_tissues) + ['Down'] * len(down_tissues)
})
all_tissues.to_csv(os.path.join(output_dir, 'selected_tissues.csv'), index=False)

# 计算气泡大小范围（基于相关性绝对值）
plotted_corr = []
for species, df in corr_data.items():
    for tissue in up_tissues + down_tissues:
        if tissue in df.index:
            plotted_corr.append(abs(df.loc[tissue, 'Spearman_Correlation']))

corr_min = np.min(plotted_corr) if plotted_corr else 0
corr_max = np.max(plotted_corr) if plotted_corr else 1
corr_range = corr_max - corr_min if corr_max > corr_min else 1

target_min_size = 100
target_max_size = 600
corr_scale = (target_max_size - target_min_size) / corr_range if corr_range > 0 else 1

# 绘图
fig, axes = plt.subplots(2, 1, figsize=(18, 8))

norm = Normalize(vmin=-1, vmax=1)
cmap = plt.cm.RdBu_r
font_scale = 1.4

for row, (tissues, regulation) in enumerate([(up_tissues, 'Up-regulation'), (down_tissues, 'Down-regulation')]):
    ax = axes[row]
    for y_idx, species in enumerate(datasets_list):
        df = corr_data[species]
        for x_idx, tissue in enumerate(tissues):
            if tissue in df.index:
                spearman_r = df.loc[tissue, 'Spearman_Correlation']
                abs_r = abs(spearman_r)
                pval = df.loc[tissue, 'Spearman_P_value']

                ax.scatter(x_idx, y_idx, s=target_min_size + (abs_r - corr_min) * corr_scale,
                           c=[spearman_r], cmap=cmap, norm=norm, edgecolor='black', linewidth=0.5, alpha=0.7, zorder=3)

                # 显著性标记
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

    ax.set_xticks(range(len(tissues)))
    ax.set_xticklabels(tissues, fontsize=int(11 * font_scale), rotation=45, ha='right')
    ax.set_yticks(range(len(datasets_list)))
    ax.set_yticklabels(datasets_list, fontsize=int(12 * font_scale))
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_title(regulation, fontsize=int(14 * font_scale), fontweight='bold')
    ax.margins(y=0.3)
    ax.invert_yaxis()

# 颜色条（左下）
legend_bottom = 0.05 * font_scale
cbar_ax = fig.add_axes([0.12, legend_bottom + 0.08, 0.3, 0.025])
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Spearman r', fontsize=int(12 * font_scale))
cbar.set_ticks([-1, 0, 1])
cbar.ax.tick_params(labelsize=int(10 * font_scale))

# 气泡大小图例（右下）
corr_med = np.median(plotted_corr) if plotted_corr else (corr_min + corr_max) / 2
size_values = [corr_min, corr_med, corr_max]
size_handles = [plt.scatter([], [], s=target_min_size + (val - corr_min) * corr_scale,
                            color='gray', alpha=0.7, edgecolor='black', linewidth=0.5) for val in size_values]
size_labels = ['Low', 'Medium', 'High']
fig.legend(size_handles, size_labels, loc='lower right', bbox_to_anchor=(0.9, legend_bottom + 0.07),
           ncol=3, frameon=False, fontsize=int(11 * font_scale))

fig.text(0.65, legend_bottom + 0.04, 'Bubble size: |Spearman r|', fontsize=int(12 * font_scale))

# 调整布局
hspace = 0.6 * font_scale
bottom = 0.25 * font_scale
plt.subplots_adjust(hspace=hspace, bottom=bottom)

# 保存
plt.savefig(os.path.join(output_dir, 'bubble_plot_bcell.png'), dpi=600, pad_inches=0.1)
plt.savefig(os.path.join(output_dir, 'bubble_plot_bcell.pdf'),
            dpi=300,
            bbox_inches='tight',
            pad_inches=0.1,
            metadata={'Creator': None, 'Producer': None})
plt.close()

print(f"Plot saved to {output_dir}")
print(f"Up-regulated tissues: {len(up_tissues)}")
print(f"Down-regulated tissues: {len(down_tissues)}")
print(f"Species included: {datasets_list}")
