#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : PCA_test.py
@Date    : 2026/1/27 22:20
@Author  : yh109
"""

# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Ellipse


def plot_pca_with_ellipses(csv_file, info_path, save_dir,
                           color_dict=None, n_proteins=None, figsize=(10, 8), filename=None):
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
        分组配色字典
    n_proteins : int, optional
        使用前n个高变蛋白，None表示全部
    figsize : tuple, optional
        图片尺寸
    filename : str, optional
        输出文件名，不含路径
    """

    if color_dict is None:
        color_dict = {'Young': '#1f77b4', 'Middle': '#ff7f0e', 'Old': '#d62728'}

    if filename is None:
        filename = "pca_plot_with_ellipses.png"

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

    # 数据标准化
    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(x_df)
    x_scaled = np.nan_to_num(x_scaled, nan=0.0)

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
    margin = 0.3
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
    plt.savefig(os.path.join(save_dir, filename), dpi=300)
    plt.close()


def jackknife_pca_analysis(csv_file, info_path, output_dir,
                           color_dict=None, n_proteins=None, figsize=(10, 8)):
    """
    留一法交叉验证：去掉每个样本后绘制PCA，用于检测异常样本

    Parameters:
    -----------
    csv_file : str
        蛋白质丰度文件路径
    info_path : str
        样本元数据文件路径
    output_dir : str
        输出图片目录
    color_dict : dict, optional
        分组配色字典
    n_proteins : int, optional
        使用前n个高变蛋白
    figsize : tuple, optional
        图片尺寸
    """

    os.makedirs(output_dir, exist_ok=True)

    # 读数据
    meta_df = pd.read_csv(info_path, index_col=0)
    x_df = pd.read_csv(csv_file, index_col=0)

    sample_names = [s.strip() for s in x_df.index.tolist()]
    meta_df.index = meta_df.index.str.strip()

    print(f"共 {len(sample_names)} 个样本")
    print("开始逐个去掉样本进行PCA...")

    # 绘制全样本PCA
    plot_pca_with_ellipses(csv_file, info_path, output_dir,
                           color_dict=color_dict, n_proteins=n_proteins,
                           figsize=figsize, filename="00_all_samples.png")
    print("已保存全样本PCA图: 00_all_samples.png")

    # 留一法
    for i, excluded_sample in enumerate(sample_names, 1):
        # 去掉该样本
        x_subset = x_df.drop(excluded_sample)

        # 临时保存子集CSV
        temp_csv = f"{output_dir}/.temp_subset.csv"
        x_subset.to_csv(temp_csv)

        try:
            output_name = f"{i:02d}_exclude_{excluded_sample}.png"
            plot_pca_with_ellipses(temp_csv, info_path, output_dir,
                                   color_dict=color_dict, n_proteins=n_proteins,
                                   figsize=figsize, filename=output_name)
            print(f"[{i}/{len(sample_names)}] 已保存: {output_name}")
        finally:
            if os.path.exists(temp_csv):
                os.remove(temp_csv)


def remove_specified_samples_pca(csv_file, info_path, output_dir,
                                  exclude_samples=None, color_dict=None,
                                  n_proteins=None, figsize=(10, 8)):
    """
    去掉指定的多个样本后绘制PCA

    Parameters:
    -----------
    csv_file : str
        蛋白质丰度文件路径
    info_path : str
        样本元数据文件路径
    output_dir : str
        输出图片目录
    exclude_samples : list, optional
        要去掉的样本名列表
    color_dict : dict, optional
        分组配色字典
    n_proteins : int, optional
        使用前n个高变蛋白
    figsize : tuple, optional
        图片尺寸
    """

    if color_dict is None:
        color_dict = {'Young': '#1f77b4', 'Middle': '#ff7f0e', 'Old': '#d62728'}

    if exclude_samples is None:
        exclude_samples = []

    os.makedirs(output_dir, exist_ok=True)

    # 读数据
    x_df = pd.read_csv(csv_file, index_col=0)
    sample_names = [s.strip() for s in x_df.index.tolist()]

    # 清理样本名并去掉指定样本
    exclude_samples = [s.strip() for s in exclude_samples]
    x_subset = x_df.drop([s for s in exclude_samples if s in x_df.index])

    print(f"原样本数: {len(sample_names)}")
    print(f"去掉样本: {exclude_samples}")
    print(f"剩余样本数: {len(x_subset)}")

    # 临时保存子集CSV
    temp_csv = f"{output_dir}/.temp_subset.csv"
    x_subset.to_csv(temp_csv)

    try:
        excluded_str = "_".join(exclude_samples) if exclude_samples else "none"
        filename = f"pca_exclude_{excluded_str}.png"
        plot_pca_with_ellipses(temp_csv, info_path, output_dir,
                               color_dict=color_dict, n_proteins=n_proteins,
                               figsize=figsize, filename=filename)
        print(f"已保存: {filename}")
    finally:
        if os.path.exists(temp_csv):
            os.remove(temp_csv)

if __name__ == '__main__':
    csv_file = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\abundance_sample_x_protein.csv'
    info_path = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\sample_meta.csv'
    output_dir = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\figure\jackknife_pca'

    jackknife_pca_analysis(csv_file, info_path, output_dir, n_proteins=500)
    # remove_specified_samples_pca(csv_file, info_path, output_dir,
    #                              exclude_samples=['F1', 'M3', 'F3'],
    #                              n_proteins=500)
