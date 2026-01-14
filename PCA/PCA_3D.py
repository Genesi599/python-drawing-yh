#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ==================== 用户输入区（仅需改这里） ============================
csv_file = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\abundance_sample_x_protein.csv'
save_dir = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\figure'
info_path = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\sample_info.csv'

# ---- 读分组信息 ----
meta_df = pd.read_csv(info_path, index_col=0)
x_df = pd.read_csv(csv_file, index_col=0)
x_idx = x_df.columns

group_ser = meta_df.reindex(x_idx)['condition']
miss = group_ser.isna().sum()
if miss:
    raise KeyError(f'有 {miss} 个样本找不到 condition！')

group_list = group_ser.tolist()

print(f'样本数: {len(x_idx)}, group_list长度: {len(group_list)}')

color_dict = {'Young': '#1f77b4', 'Middle': '#ff7f0e', 'Old': '#d62728'}

# ==================== PCA分析 ============================
print(f"x_df shape: {x_df.shape}")
print(f"meta_df shape: {meta_df.shape}")

if len(x_idx) != len(group_list):
    raise ValueError(f'样本数不匹配: x_df样本数为 {len(x_idx)}, group_list样本数为 {len(group_list)}')

scaler = StandardScaler()
x_scaled = scaler.fit_transform(x_df.T)

# PCA
pca = PCA(n_components=3)  # 设置为3个主成分
pca_result = pca.fit_transform(x_scaled)

# 创建 DataFrame 存储 PCA 结果
pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2', 'PC3'], index=x_idx)
pca_df['Group'] = group_list

# ==================== 绘制3D交互式PCA图 ============================
# 使用 Plotly 绘制 3D 散点图
fig = px.scatter_3d(
    pca_df,
    x='PC1',
    y='PC2',
    z='PC3',
    color='Group',
    color_discrete_map=color_dict,  # 使用自定义颜色
    labels={'Group': 'Condition'},
    title='3D PCA of Protein Abundance',
    hover_name=pca_df.index
)

# 添加每组的球体
for group in pca_df['Group'].unique():
    group_data = pca_df[pca_df['Group'] == group]
    mean_x = group_data['PC1'].mean()
    mean_y = group_data['PC2'].mean()
    mean_z = group_data['PC3'].mean()

    std_x = group_data['PC1'].std()
    std_y = group_data['PC2'].std()
    std_z = group_data['PC3'].std()

    # 生成球体坐标
    u = np.linspace(0, 2 * np.pi, 100)  # 0到2π生成100个点
    v = np.linspace(0, np.pi, 50)  # 0到π生成50个点
    x_sphere = mean_x + 2 * std_x * np.outer(np.cos(u), np.sin(v))  # 球面方程X坐标
    y_sphere = mean_y + 2 * std_y * np.outer(np.sin(u), np.sin(v))  # 球面方程Y坐标
    z_sphere = mean_z + 2 * std_z * np.outer(np.ones(np.size(u)), np.cos(v))  # 球面方程Z坐标

    # 添加球体到图中
    fig.add_trace(go.Mesh3d(
        x=x_sphere.flatten(),
        y=y_sphere.flatten(),
        z=z_sphere.flatten(),
        color=color_dict[group],
        opacity=0.2,
        showlegend=False  # 不显示图例
    ))

# 显示图形
fig.show()

# 另存为HTML文件
fig.write_html(f"{save_dir}/pca_3d_plot_with_spheres.html")