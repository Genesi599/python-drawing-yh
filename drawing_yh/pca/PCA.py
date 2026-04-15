# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.patches import Ellipse


def plot_pca_with_ellipses(csv_file, info_path, save_dir,
                           color_dict=None, n_proteins=None, figsize=(10, 8)):
    """
    绘制PCA图并添加分组椭圆

    Parameters:
    -----------
    csv_file : str
        蛋白质丰度文件路径（样本×蛋白质）
    info_path : str
        样本元数据文件路径，需包含'condition'列
    save_dir : str
        输出图片保存目录
    color_dict : dict, optional
        分组配色字典，默认为{'Young': '#1f77b4', 'Middle': '#ff7f0e', 'Old': '#d62728'}
    n_proteins : int, optional
        使用前n个高变蛋白，None表示全部
    figsize : tuple, optional
        图片尺寸
    """

    if color_dict is None:
        color_dict = {'Young': '#1f77b4', 'Middle': '#ff7f0e', 'Old': '#d62728'}

    # 读数据
    meta_df = pd.read_csv(info_path, index_col=0)
    x_df = pd.read_csv(csv_file, index_col=0)

    # 样本名清理
    sample_names = [s.strip() for s in x_df.index.tolist()]
    meta_df.index = meta_df.index.str.strip()

    # 取分组
    group_ser = meta_df.reindex(sample_names)['condition']
    if group_ser.isna().sum() > 0:
        raise KeyError(f"有 {group_ser.isna().sum()} 个样本找不到 condition！")

    group_list = group_ser.tolist()
    print(f"样本数: {len(sample_names)}, group_list 长度: {len(group_list)}")

    # 数据标准化
    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(x_df)

    if n_proteins:
        protein_var = x_scaled.var(axis=0)
        top_indices = np.argsort(protein_var)[-n_proteins:]
        x_scaled = x_scaled[:, top_indices]

    # PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(x_scaled)

    # 构建PCA DataFrame
    pca_df = pd.DataFrame(data=pca_result,
                          columns=['PC1', 'PC2'],
                          index=sample_names)
    pca_df['Group'] = group_list
    pca_df = pca_df.dropna(subset=['Group'])

    # 绘图
    plt.figure(figsize=figsize)
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Group',
                    palette=color_dict, s=300)

    # 添加样本名标签
    for idx, row in pca_df.iterrows():
        plt.annotate(idx, (row['PC1'], row['PC2']), xytext=(10, -10),
                     textcoords='offset points', fontsize=16, ha='left')

    # 添加椭圆
    for group in pca_df['Group'].unique():
        group_data = pca_df[pca_df['Group'] == group]
        points = group_data[['PC1', 'PC2']].values
        mean = points.mean(axis=0)
        cov = np.cov(points.T)
        eigenvalues, eigenvectors = np.linalg.eig(cov)
        angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
        width = 2 * 2 * np.sqrt(eigenvalues[0])
        height = 2 * 2 * np.sqrt(eigenvalues[1])
        ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle,
                          edgecolor=color_dict[group], facecolor='none',
                          linewidth=2, linestyle='--')
        plt.gca().add_patch(ellipse)

    # 调整坐标轴
    pc1_min, pc1_max = pca_df['PC1'].min(), pca_df['PC1'].max()
    pc2_min, pc2_max = pca_df['PC2'].min(), pca_df['PC2'].max()
    margin = 0.1
    plt.xlim(pc1_min - margin * (pc1_max - pc1_min),
             pc1_max + margin * (pc1_max - pc1_min))
    plt.ylim(pc2_min - margin * (pc2_max - pc2_min),
             pc2_max + margin * (pc2_max - pc2_min))

    # 标题和标签
    var1 = pca.explained_variance_ratio_[0] * 100
    var2 = pca.explained_variance_ratio_[1] * 100
    title_str = 'PCA of Protein Abundance with Group Ellipses'
    if n_proteins:
        title_str += f' (Top {n_proteins} Variable Proteins)'

    plt.xlabel(f'Principal Component 1 ({var1:.1f}%)', fontsize=14)
    plt.ylabel(f'Principal Component 2 ({var2:.1f}%)', fontsize=14)
    plt.title(title_str, fontsize=16)
    plt.legend(title='Conditions', fontsize=12, title_fontsize=12)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/pca_plot_with_ellipses.png", dpi=300)
    plt.show()


# 使用示例
if __name__ == '__main__':
    csv_file = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\abundance_sample_x_protein.csv'
    save_dir = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\figure'
    info_path = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\sample_meta.csv'

    plot_pca_with_ellipses(csv_file, info_path, save_dir, n_proteins=500)
