#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import numpy as np
from matplotlib.patches import Ellipse

# ==================== 用户输入区（仅需改这里） ============================
csv_file = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\abundance_sample_x_protein.csv'
save_dir = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\figure'
info_path = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\sample_info.csv'

# ---- 读分组信息 ----
meta_df = pd.read_csv(info_path, index_col=0)  # 索引=样本名
x_df = pd.read_csv(csv_file, index_col=0)  # 行=蛋白，列=样本
x_idx = x_df.columns  # 样本名列表

# 保留全部三组
group_ser = meta_df.reindex(x_idx)['condition']
miss = group_ser.isna().sum()
if miss:
    raise KeyError(f'有 {miss} 个样本找不到 condition！')

# 确保 group_list 的长度与样本数量一致
group_list = group_ser.tolist()

# 三组配色
color_dict = {'Young': '#1f77b4', 'Middle': '#ff7f0e', 'Old': '#d62728'}
# ======================== 结束输入区 ========================================

# ==================== PLS-DA分析与可视化 ============================
print(f"x_df shape: {x_df.shape}")  # 输出(蛋白数量, 样本数量)
print(f"meta_df shape: {meta_df.shape}")  # 输出(样本数量, 属性)

# 确认样本数是否一致
x_idx = x_df.columns
if len(x_idx) != len(group_list):
    raise ValueError(f'样本数不匹配: x_df样本数为 {len(x_idx)}, group_list样本数为 {len(group_list)}')

# 数据标准化（注意这里的样本在数据集中为列）
scaler = StandardScaler()
x_scaled = scaler.fit_transform(x_df.T)  # 转置操作使列为样本

# 将组别转换为数值标签
y_labels = pd.Series(group_list).astype('category').cat.codes  # 转换为数值型标签

# PLS-DA
pls = PLSRegression(n_components=2)
pls.fit(x_scaled, y_labels)

# 获取模型的结果
X_pls = pls.transform(x_scaled)

# 创建 DataFrame 存储 PLS-DA 结果
pls_df = pd.DataFrame(data=X_pls, columns=['PLS1', 'PLS2'], index=x_idx)  # 使用样本名作为索引
pls_df['Group'] = group_list

# 绘制 PLS-DA 图
plt.figure(figsize=(10, 8))

# 使用 seaborn 绘制散点图
sns.scatterplot(data=pls_df, x='PLS1', y='PLS2', hue='Group', palette=color_dict, s=100)

# 为每组样本添加圈
unique_groups = pls_df['Group'].unique()

for group in unique_groups:
    group_data = pls_df[pls_df['Group'] == group]
    mean_x = group_data['PLS1'].mean()
    mean_y = group_data['PLS2'].mean()

    # 计算标准差
    std_x = group_data['PLS1'].std()
    std_y = group_data['PLS2'].std()

    # 创建椭圆并添加到图中
    ellipses = Ellipse((mean_x, mean_y), width=2 * std_x, height=2 * std_y,
                       edgecolor=color_dict[group], facecolor='none',
                       linewidth=2, linestyle='--')  # 椭圆样式设置
    plt.gca().add_patch(ellipses)

# 设置标题和坐标轴标签
plt.title('PLS-DA of Protein Abundance with Group Ellipses')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='Conditions')
plt.grid()

# 保存和显示图形
plt.savefig(f"{save_dir}/plsda_plot.png")
plt.show()

# 可选：呈现混淆矩阵
# 获取 PLS 预测的分数并转换为类别
predictions = pls.predict(x_scaled)
predicted_labels = np.round(predictions).astype(int).flatten()  # 展平并四舍五入到最近的整数

# 确保预测结果和原标签的一致性
if len(set(predicted_labels)) != len(set(y_labels)):
    print("模型未能识别出足够的类别，可能需要调整预测阈值。")
else:
    # 修正显示标签的处理
    unique_labels = pd.Series(group_list).astype('category').cat.categories
    cm = confusion_matrix(y_labels, predicted_labels)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=unique_labels)
    disp.plot()
    plt.title('Confusion Matrix for PLS-DA')
    plt.savefig(f"{save_dir}/plsda_confusion_matrix.png")
    plt.show()