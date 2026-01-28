import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from collections import defaultdict

# 1. 让 matplotlib 别把字体子集化
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42          # 已经写了，保留
mpl.rcParams['ps.fonttype'] = 42
# 关键：关闭子集化（>=3.6 有效）
mpl.rcParams['pdf.use14corefonts'] = False   # 不用 14 种核心字体
# 2. 指定一个系统里肯定有的字体，避免 mpl 走自带字体
mpl.rcParams['font.family'] = 'Arial'      # 或 'DejaVu Sans', 'SimHei' 等

# 数据路径

meta_path = r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\feature_meta.csv"

data_paths = {
    'all': r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\by_sex\protein_age_correlations_significant_all.csv",
    'F': r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\by_sex\protein_age_correlations_significant_F.csv",
    'M': r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\by_sex\protein_age_correlations_significant_M.csv"
}

full_data_paths = {
    'all': r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\by_sex\protein_age_correlations_all.csv",
    'F': r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\by_sex\protein_age_correlations_F.csv",
    'M': r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\by_sex\protein_age_correlations_M.csv"
}

meta = pd.read_csv(meta_path)
protein_to_gene = dict(zip(meta['protein'], meta['Gene name']))

corr_data = {}
full_corr_data = {}
for key, path in data_paths.items():
    df = pd.read_csv(path)
    df['Gene name'] = df['protein'].map(protein_to_gene)
    df['-log2_pval'] = -np.log2(df['pearson_pval'])
    corr_data[key] = df

    full_df = pd.read_csv(full_data_paths[key])
    full_df['Gene name'] = full_df['protein'].map(protein_to_gene)
    full_df['-log2_pval'] = -np.log2(full_df['pearson_pval'])
    full_corr_data[key] = full_df


def sort_genes_by_freq_and_strength(corr_data, direction='up'):
    counts = defaultdict(int)
    sum_abs_r = defaultdict(float)
    for key, df in corr_data.items():
        if direction == 'up':
            subset = df[df['pearson_r'] > 0]
        else:
            subset = df[df['pearson_r'] < 0]
        for gene in subset['Gene name'].unique():
            val = subset.loc[subset['Gene name'] == gene, 'pearson_r'].values[0]
            counts[gene] += 1
            sum_abs_r[gene] += abs(val)
    avg_abs_r = {g: (sum_abs_r[g] / counts[g]) if counts[g] > 0 else 0.0 for g in counts}
    sorted_genes = sorted(list(counts.keys()), key=lambda g: (-counts[g], -avg_abs_r.get(g, 0.0)))
    return sorted_genes


up_genes = sort_genes_by_freq_and_strength(corr_data, direction='up')[:50]
down_genes = sort_genes_by_freq_and_strength(corr_data, direction='down')[:50]

# ↓ 插入保存代码 ↓
gene_to_protein = {v: k for k, v in protein_to_gene.items()}

all_genes = pd.DataFrame({
    'Gene name': list(up_genes) + list(down_genes),
    'protein': [gene_to_protein.get(g, g) for g in list(up_genes) + list(down_genes)],
    'direction': ['Up'] * len(up_genes) + ['Down'] * len(down_genes)
})
all_genes.to_csv('D:\\Projects\\Bone_Marrow_Aging\\proteomics\\analysis\\data\\selected_genes_all.csv', index=False)

# ↑ 插入保存代码 ↑

fig, axes = plt.subplots(2, 1, figsize=(18, 8))
datasets_list = [('all', 'All'), ('F', 'Female'), ('M', 'Male')]

norm = Normalize(vmin=-1, vmax=1)
cmap = plt.cm.RdBu_r
font_scale = 1.4




plotted_abundances = []

for row, (genes, regulation) in enumerate([(up_genes, 'Up-regulation'), (down_genes, 'Down-regulation')]):
    for y_idx, (key, label) in enumerate(datasets_list):
        df = full_corr_data[key]
        for gene in genes:
            match = df[df['Gene name'] == gene]
            if len(match) > 0:
                abundance = match['mean_abundance'].values[0]
                plotted_abundances.append(abundance)

abundance_min = np.min(plotted_abundances)
abundance_max = np.max(plotted_abundances)
abundance_med = np.median(plotted_abundances)
abundance_range = abundance_max - abundance_min

target_min_size = 100
target_max_size = 600
abundance_scale = (target_max_size - target_min_size) / abundance_range if abundance_range > 0 else 1

for row, (genes, regulation) in enumerate([(up_genes, 'Up-regulation'), (down_genes, 'Down-regulation')]):
    ax = axes[row]
    for y_idx, (key, label) in enumerate(datasets_list):
        df = full_corr_data[key]
        for x_idx, gene in enumerate(genes):
            match = df[df['Gene name'] == gene]
            if len(match) > 0:
                pearson_r = match['pearson_r'].values[0]
                abundance = match['mean_abundance'].values[0]
                plotted_abundances.append(abundance)
                ax.scatter(x_idx, y_idx, s=target_min_size + (abundance - abundance_min) * abundance_scale,
                           c=[pearson_r], cmap=cmap, norm=norm, edgecolor='black', linewidth=0.5, alpha=0.7, zorder=3)

                # 添加显著性标记
                pval = match['pearson_pval'].values[0]
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



    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, fontsize=int(11 * font_scale), rotation=45, ha='right')
    ax.set_yticks(range(len(datasets_list)))
    ax.set_yticklabels([label for _, label in datasets_list], fontsize=int(12*font_scale))
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_title(regulation, fontsize=int(14*font_scale), fontweight='bold')
    ax.margins(y=0.3)
    ax.invert_yaxis()


legend_bottom = 0.05 * font_scale
# 颜色条放左下角
cbar_y = 0.05 * font_scale
cbar_ax = fig.add_axes([0.12, legend_bottom + 0.08, 0.3, 0.025])
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Pearson r', fontsize=int(12*font_scale))
cbar.set_ticks([-1, 0, 1])
cbar.ax.tick_params(labelsize=int(10*font_scale))


# 气泡大小图例放右下角

abundance_min = np.min(plotted_abundances)
abundance_max = np.max(plotted_abundances)
abundance_med = np.median(plotted_abundances)

size_values = [abundance_min, abundance_med, abundance_max]
# 图例处
size_handles = [plt.scatter([], [], s=target_min_size + (val - abundance_min) * abundance_scale,
                            color='gray', alpha=0.7, edgecolor='black', linewidth=0.5) for val in size_values]
size_labels = ['Low', 'Medium', 'High']
legend_y = 0.05 * font_scale
fig.legend(size_handles, size_labels, loc='lower right', bbox_to_anchor=(0.9, legend_bottom + 0.07),
           ncol=3, frameon=False, fontsize=int(11*font_scale))

text_y = 0.02 * font_scale
fig.text(0.65, legend_bottom + 0.04, 'Bubble size: mean abundance', fontsize=int(12*font_scale))



hspace = 0.6 * font_scale
bottom = 0.25 * font_scale
plt.subplots_adjust(hspace=hspace, bottom=bottom)
plt.savefig('bubble_plot_adjusted.png', dpi=600, pad_inches=0.1)
plt.savefig('bubble_plot_adjusted.pdf',
            dpi=300,             # 矢量图 dpi 其实无所谓
            bbox_inches='tight',
            pad_inches=0.1,
            metadata={'Creator': None, 'Producer': None})   # 可选：去掉多余元数据


