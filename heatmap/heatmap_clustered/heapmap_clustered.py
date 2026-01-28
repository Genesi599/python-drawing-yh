#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : heapmap_clustered.py
@Date    : 2026/1/22 16:50
@Author  : yh109
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
import scanpy as sc
from matplotlib.colors import PowerNorm
from collections import defaultdict

# 文件路径
adata_path = "D:/Projects/Bone_Marrow_Aging/sc-seq/data/bone_marrow_aging_processed.h5ad"
protein_corr_paths = {
    'all': "D:/Projects/Bone_Marrow_Aging/proteomics/analysis/data/by_sex/protein_age_correlations_significant_all.csv",
    'F': "D:/Projects/Bone_Marrow_Aging/proteomics/analysis/data/by_sex/protein_age_correlations_significant_F.csv",
    'M': "D:/Projects/Bone_Marrow_Aging/proteomics/analysis/data/by_sex/protein_age_correlations_significant_M.csv"
}
feature_meta_path = "D:/Projects/Bone_Marrow_Aging/proteomics/analysis/data/feature_meta.csv"

# 读取数据
adata = sc.read_h5ad(adata_path)
meta_df = pd.read_csv(feature_meta_path)

# 读取三个数据集并映射基因名
corr_data = {}
protein_to_gene = dict(zip(meta_df['protein'], meta_df['Gene name']))
for key, path in protein_corr_paths.items():
    df = pd.read_csv(path)
    df['Gene name'] = df['protein'].map(protein_to_gene)
    corr_data[key] = df

# 按相同逻辑筛选基因
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

n_genes = 50
genes_pos = sort_genes_by_freq_and_strength(corr_data, direction='up')[:n_genes]
genes_pos = [g for g in genes_pos if g is not None and g in adata.var_names]

genes_neg = sort_genes_by_freq_and_strength(corr_data, direction='down')[:n_genes]
genes_neg = [g for g in genes_neg if g is not None and g in adata.var_names]

# 定义函数用于计算表达矩阵
def compute_expression_matrix(genes, adata, cell_types):
    expr_matrix = []
    for gene in genes:
        gene_idx = adata.var_names.get_loc(gene)
        expr_data = adata.X[:, gene_idx]
        if hasattr(expr_data, 'toarray'):
            expr_data = expr_data.toarray().ravel()
        else:
            expr_data = np.asarray(expr_data).ravel()
        gene_means = [expr_data[adata.obs['cluster_cell_type'] == ct].mean() for ct in cell_types]
        expr_matrix.append(gene_means)
    return pd.DataFrame(expr_matrix, index=genes, columns=cell_types)

# 定义绘图函数
def plot_clustermap(expr_df, cmap, method, output_prefix, height_per_gene=0.2, base_height=3):
    fig_height = base_height + len(expr_df) * height_per_gene

    cluster = sns.clustermap(
        data=expr_df,
        method=method,
        cmap=cmap,
        norm=PowerNorm(gamma=0.5),
        figsize=(12, fig_height),
        xticklabels=True,
        yticklabels=True,
        row_cluster=True,
        col_cluster=True
    )

    left, bottom, width, height = 0.35, 0.15, 0.5, 0.75
    tree_size, tree_gap = 0.08, 0.01

    cluster.ax_heatmap.set_position([left, bottom, width, height])
    cluster.ax_row_dendrogram.set_position([left - tree_size - tree_gap, bottom, tree_size, height])
    cluster.ax_col_dendrogram.set_position([left, bottom + height + tree_gap, width, tree_size])

    for ax in [cluster.ax_row_dendrogram, cluster.ax_col_dendrogram]:
        for line in ax.collections:
            line.set_linewidth(2)

    cluster.ax_heatmap.set_yticklabels(
        cluster.ax_heatmap.get_yticklabels(),
        fontsize=16, fontweight='bold', ha='left'
    )
    cluster.ax_heatmap.tick_params(axis='y', labelleft=False, labelright=True, length=0)
    cluster.ax_heatmap.set_xticklabels(
        cluster.ax_heatmap.get_xticklabels(),
        fontsize=16, fontweight='bold', ha='right', rotation=45
    )

    cluster.cax.remove()
    cbar_ax = plt.axes([left - 0.2, bottom + height + tree_gap + 0.05, 0.15, 0.015])
    cbar = plt.colorbar(cluster.ax_heatmap.collections[0], cax=cbar_ax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=16)

    plt.savefig(f'{output_prefix}_{method}.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

# 提取细胞类型
cluster_col = 'cluster_cell_type'
cell_types = sorted(adata.obs[cluster_col].unique())

# 计算正相关基因表达矩阵（内存中处理）
pos_expr_df = compute_expression_matrix(genes_pos, adata, cell_types)
# 归一化
def z_score(df):
    return (df - df.mean(axis=1).values.reshape(-1, 1)) / df.std(axis=1).values.reshape(-1, 1)

z_score_df_pos = pos_expr_df.loc[pos_expr_df.mean(axis=1).sort_values(ascending=False).index]

# 计算负相关基因表达矩阵（内存中处理）
neg_expr_df = compute_expression_matrix(genes_neg, adata, cell_types)
z_score_df_neg = neg_expr_df.loc[neg_expr_df.mean(axis=1).sort_values(ascending=False).index]

# 定义渐变色
red_cmap = LinearSegmentedColormap.from_list("RedGrad", ["white", "#FF333E"])
blue_cmap = LinearSegmentedColormap.from_list("BlueGrad", ["white", "#1064F3"])

linkage_methods = ['ward']

# 调用绘图函数
for method in linkage_methods:
    plot_clustermap(z_score_df_pos, red_cmap, method, 'D:/Projects/Bone_Marrow_Aging/proteomics/analysis/figure/clustermap_positive')
    plot_clustermap(z_score_df_neg, blue_cmap, method, 'D:/Projects/Bone_Marrow_Aging/proteomics/analysis/figure/clustermap_negative')