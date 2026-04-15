#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==================== 用户输入区（仅需改这里） ============================
csv_file = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\abundance_sample_x_protein.csv'
save_dir = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\figure'
info_path = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\sample_info.csv'

# ---- 读分组信息 ----
# ---- 读分组信息 ----
meta_df = pd.read_csv(info_path, index_col=0)          # 索引=样本名
x_df = pd.read_csv(csv_file, index_col=0)              # 行=蛋白，列=样本
x_idx = x_df.columns                                    # 样本名列表

# 保留全部三组
group_ser = meta_df.reindex(x_idx)['condition']
miss = group_ser.isna().sum()
if miss:
    raise KeyError(f'有 {miss} 个样本找不到 condition！')

# 确保 group_list 的长度与样本数量一致
group_list = group_ser.tolist()

# 根据条件长度打印调试信息
print(f'样本数: {len(x_idx)}, group_list长度: {len(group_list)}')

# 三组配色
color_dict = {'Young': '#1f77b4', 'Middle': '#ff7f0e', 'Old': '#d62728'}
# ======================== 结束输入区 ========================================


from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ==================== PCA分析与可视化 ============================
# 检查输入数据的维度
print(f"x_df shape: {x_df.shape}")  # 输出(蛋白数量, 样本数量)
print(f"meta_df shape: {meta_df.shape}")  # 输出(样本数量, 属性)

# 确认样本数是否一致
x_idx = x_df.columns
if len(x_idx) != len(group_list):
    raise ValueError(f'样本数不匹配: x_df样本数为 {len(x_idx)}, group_list样本数为 {len(group_list)}')

# 数据标准化（注意这里的样本在数据集中为列）
scaler = StandardScaler()
x_scaled = scaler.fit_transform(x_df.T)  # 转置操作使列为样本

# PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(x_scaled)

# 创建 DataFrame 存储 PCA 结果
pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'], index=x_idx)  # 使用样本名作为索引
pca_df['Group'] = group_list

# 绘制 PCA 图
import numpy as np
from matplotlib.patches import Ellipse

# 绘制 PCA 图
plt.figure(figsize=(10, 8))

# 使用 seaborn 绘制散点图
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Group', palette=color_dict, s=100)

# 为每组样本添加圈
unique_groups = pca_df['Group'].unique()

# 根据每组的样本添加椭圆
for group in unique_groups:
    # 获取当前组的样本
    group_data = pca_df[pca_df['Group'] == group]

    # 计算当前组的中心点和标准差
    mean_x = group_data['PC1'].mean()
    mean_y = group_data['PC2'].mean()
    std_x = group_data['PC1'].std()
    std_y = group_data['PC2'].std()

    # 这里使用 2 倍的标准差作为椭圆的半长轴和半短轴
    ellipses = Ellipse((mean_x, mean_y), width=2 * std_x, height=2 * std_y, edgecolor=color_dict[group],
                       facecolor='none', linewidth=2, linestyle='--')  # 设置椭圆的样式

    plt.gca().add_patch(ellipses)

# 设置标题和坐标轴标签
plt.title('PCA of Protein Abundance with Group Ellipses')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='Conditions')
plt.grid()

# 保存和显示图形
plt.savefig(f"{save_dir}/pca_plot_with_ellipses.png")
plt.show()