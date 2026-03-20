import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from collections import defaultdict
import matplotlib.font_manager as fm


# 数据路径
meta_path = r"D:\Projects\Bone_Marrow_Aging\metabonomics\analysis\data\feature_meta.csv"
data_paths = {
    'all': r"D:\Projects\Bone_Marrow_Aging\metabonomics\analysis\data\by_sex\protein_age_correlations_significant_all.csv",
    'F': r"D:\Projects\Bone_Marrow_Aging\metabonomics\analysis\data\by_sex\protein_age_correlations_significant_F.csv",
    'M': r"D:\Projects\Bone_Marrow_Aging\metabonomics\analysis\data\by_sex\protein_age_correlations_significant_M.csv"
}

full_data_paths = {
    'all': r"D:\Projects\Bone_Marrow_Aging\metabonomics\analysis\data\by_sex\protein_age_correlations_all.csv",
    'F': r"D:\Projects\Bone_Marrow_Aging\metabonomics\analysis\data\by_sex\protein_age_correlations_F.csv",
    'M': r"D:\Projects\Bone_Marrow_Aging\metabonomics\analysis\data\by_sex\protein_age_correlations_M.csv"
}


# 列名参数（根据实际列名调整）
meta_columns = {
    'protein_col': 'COMP_ID',
    'gene_name_col': 'COMPOUND_Name',
}

corr_columns = {
    'protein_col': 'protein',
    'pearson_pval_col': 'pearson_pval',
    'pearson_r_col': 'pearson_r',
    'mean_abundance_col': 'mean_abundance',
}

plot_columns = {
    'compound_col': 'COMPOUND_Name',
}


meta = pd.read_csv(meta_path)
protein_to_gene = dict(zip(meta[meta_columns['protein_col']], meta[meta_columns['gene_name_col']]))

corr_data = {}
full_corr_data = {}
for key, path in data_paths.items():
    df = pd.read_csv(path)
    df[plot_columns['compound_col']] = df[corr_columns['protein_col']].map(protein_to_gene)
    df['-log2_pval'] = -np.log2(df[corr_columns['pearson_pval_col']])
    corr_data[key] = df

    full_df = pd.read_csv(full_data_paths[key])
    full_df[plot_columns['compound_col']] = full_df[corr_columns['protein_col']].map(protein_to_gene)
    full_df['-log2_pval'] = -np.log2(full_df[corr_columns['pearson_pval_col']])
    full_corr_data[key] = full_df


# 读取meta数据并建立映射
meta = pd.read_csv(meta_path)
protein_to_gene = dict(zip(meta[meta_columns['protein_col']], meta[meta_columns['gene_name_col']]))

# 读取相关性数据
corr_data = {}
for key, path in data_paths.items():
    df = pd.read_csv(path)
    df['protein'] = df['protein'].astype(str)  # 转换为字符串
    df[plot_columns['compound_col']] = df[corr_columns['protein_col']].map(protein_to_gene)
    df['-log2_pval'] = -np.log2(df[corr_columns['pearson_pval_col']])
    corr_data[key] = df


def sort_genes_by_freq_and_strength(corr_data, corr_columns, compound_col='Gene name', direction='up'):
    counts = defaultdict(int)
    sum_abs_r = defaultdict(float)
    pearson_r_col = corr_columns['pearson_r_col']

    for key, df in corr_data.items():
        if direction == 'up':
            subset = df[df[pearson_r_col] > 0]
        else:
            subset = df[df[pearson_r_col] < 0]
        for gene in subset[compound_col].unique():
            gene_data = subset[subset[compound_col] == gene]
            if len(gene_data) > 0:
                val = gene_data[pearson_r_col].iloc[0]
                counts[gene] += 1
                sum_abs_r[gene] += abs(val)
    avg_abs_r = {g: (sum_abs_r[g] / counts[g]) if counts[g] > 0 else 0.0 for g in counts}
    sorted_genes = sorted(list(counts.keys()), key=lambda g: (-counts[g], -avg_abs_r.get(g, 0.0)))
    return sorted_genes


# 获取前50个基因
up_genes = sort_genes_by_freq_and_strength(corr_data, corr_columns, plot_columns['compound_col'], direction='up')[:40]
down_genes = sort_genes_by_freq_and_strength(corr_data, corr_columns, plot_columns['compound_col'], direction='down')[:40]

# ↓ 插入保存代码 ↓
compound_to_comp_id = {v: k for k, v in protein_to_gene.items()}

all_compounds = pd.DataFrame({
    'COMPOUND_Name': list(up_genes) + list(down_genes),
    'COMP_ID': [compound_to_comp_id.get(c, c) for c in list(up_genes) + list(down_genes)],
    'direction': ['Up'] * len(up_genes) + ['Down'] * len(down_genes)
})
all_compounds.to_csv(r'D:\Projects\Bone_Marrow_Aging\metabonomics\analysis\data\selected_metabolites_all.csv', index=False)
# ↑ 插入保存代码 ↑

# 开始绘图
datasets_list = [('all', 'All'), ('F', 'Female'), ('M', 'Male')]

norm = Normalize(vmin=-1, vmax=1)
cmap = plt.cm.RdBu_r
font_scale = 1.5


def plot_regulation(genes, regulation_label):
    fig, ax = plt.subplots(figsize=(18, 10))
    norm = Normalize(vmin=-1, vmax=1)
    cmap = plt.cm.RdBu_r

    # 底部左侧的颜色条轴
    cbar_ax = fig.add_axes([0.15, 0.08, 0.30, 0.02])

    # 第一遍：收集实际绘制的点的丰度
    plotted_abundances = []
    for y_idx, (key, label) in enumerate(datasets_list):
        df = full_corr_data[key]
        for gene in genes:
            match = df[df[plot_columns['compound_col']] == gene]
            if len(match) > 0:
                abundance = match[corr_columns['mean_abundance_col']].values[0]
                plotted_abundances.append(abundance)

    # 计算丰度范围和缩放
    abundance_min = np.min(plotted_abundances)
    abundance_max = np.max(plotted_abundances)
    abundance_med = np.median(plotted_abundances)
    abundance_range = abundance_max - abundance_min

    target_min_size = 200
    target_max_size = 700
    abundance_scale = (target_max_size - target_min_size) / abundance_range if abundance_range > 0 else 1

    # 第二遍：绘图
    for y_idx, (key, label) in enumerate(datasets_list):
        df = full_corr_data[key]
        for x_idx, gene in enumerate(genes):
            match = df[df[plot_columns['compound_col']] == gene]
            if len(match) > 0:
                pearson_r = match[corr_columns['pearson_r_col']].values[0]
                abundance = match[corr_columns['mean_abundance_col']].values[0]
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
    ax.set_xticklabels(genes, fontsize=int(11 * font_scale), rotation=50, ha='right')
    ax.set_yticks(range(len(datasets_list)))
    ax.set_yticklabels([label for _, label in datasets_list], fontsize=int(12*font_scale))
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_title(regulation_label, fontsize=int(14*font_scale), fontweight='bold')
    ax.margins(y=0.3)
    ax.invert_yaxis()

    # 颜色条
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('Pearson r', fontsize=int(12*font_scale))
    cbar.set_ticks([-1, 0, 1])
    cbar.ax.tick_params(labelsize=int(10*font_scale))

    # 气泡大小图例
    size_values = [abundance_min, abundance_med, abundance_max]
    size_handles = [plt.scatter([], [], s=target_min_size + (val - abundance_min) * abundance_scale,
                                color='gray', alpha=0.7, edgecolor='black', linewidth=0.5) for val in size_values]
    size_labels = ['Low', 'Medium', 'High']
    fig.legend(size_handles, size_labels, loc='lower right',
               bbox_to_anchor=(0.95, 0.10), bbox_transform=fig.transFigure,
               ncol=3, frameon=False, fontsize=int(11 * font_scale))
    fig.text(0.8, 0.06, 'Mean Abundance', ha='center', fontsize=int(11 * font_scale))

    plt.subplots_adjust(bottom=0.6, left=0.2)
    plt.savefig(f'bubble_plot_{regulation_label.lower().replace(" ", "_")}.png', dpi=300, pad_inches=0.1)
    plt.show()

# 分别绘制
plot_regulation(up_genes, 'Up-regulation')
plot_regulation(down_genes, 'Down-regulation')
