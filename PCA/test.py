#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : test.py
@Date    : 2026/1/23 13:22
@Author  : yh109
"""
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