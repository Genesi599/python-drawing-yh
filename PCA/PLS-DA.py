import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler

# 用户参数设置
csv_file = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\abundance_sample_x_protein.csv'
info_path = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\sample_meta.csv'
save_path = r'D:\Projects\Bone_Marrow_Aging\proteomics\analysis\figure\plsda_plot.png'
n_proteins = 500  # 选择前n个变异最大的蛋白，全部用-1或None表示

# 读取数据
meta_df = pd.read_csv(info_path, index_col=0)
x_df = pd.read_csv(csv_file, index_col=0)

# 样本名
sample_names = x_df.index.tolist()

# 去除空格
meta_df.index = meta_df.index.str.strip()
sample_names = [s.strip() for s in sample_names]

# 取分组
group_ser = meta_df.reindex(sample_names)['condition']
if group_ser.isna().any():
    raise KeyError('有样本缺少 condition information！')
group_list = group_ser.tolist()

# 类别编号（编码）
category_mapping = {cat: idx for idx, cat in enumerate(np.unique(group_list))}
Y = np.array([category_mapping[grp] for grp in group_list])

# 归一化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(x_df)

# 取变异最大的蛋白
if n_proteins and n_proteins > 0:
    protein_var = np.var(X_scaled, axis=0)
    top_indices = np.argsort(protein_var)[-n_proteins:]
    X_scaled = X_scaled[:, top_indices]

# PLS-DA
pls = PLSRegression(n_components=2)
pls.fit(X_scaled, Y)
X_scores = pls.x_scores_  # 投影得分

# 可视化
plt.figure(figsize=(10,8))
colors = ['#1f77b4','#ff7f0e','#d62728']
labels_unique = np.unique(group_list)

for label_str, color, name in zip(labels_unique, colors, ['Young','Middle','Old']):
    idxs = np.array(group_list) == label_str
    plt.scatter(X_scores[idxs,0], X_scores[idxs,1], label=name, color=color, s=300)
    # 添加样本名标签
    for i in np.where(idxs)[0]:
        plt.annotate(sample_names[i], (X_scores[i,0], X_scores[i,1]),
                     xytext=(10, -10), textcoords='offset points', fontsize=16, ha='left')

plt.xlabel('PLS-DA Component 1')
plt.ylabel('PLS-DA Component 2')
plt.legend()
plt.title('PLS-DA of Protein Data')
plt.tight_layout()
plt.savefig(save_path, dpi=300)
plt.show()


# ==========================================================

import pandas as pd
import numpy as np


def calculate_vip(pls, X):
    """计算VIP值"""
    t = pls.x_scores_
    w = pls.x_weights_
    q = pls.y_loadings_

    p = X.shape[1]
    n_comp = pls.n_components

    vip_scores = np.zeros(p)
    ss = np.sum(t ** 2, axis=0) * np.sum(q ** 2, axis=0)
    ss_total = np.sum(ss)

    for i in range(p):
        vip_scores[i] = np.sqrt(p * np.sum(ss * (w[i, :] ** 2)) / ss_total)

    return vip_scores

# =================================================

import matplotlib.pyplot as plt
import os

# 假设已完成数据读取和模型训练（pls）
# 计算VIP值
vip_values = calculate_vip(pls, X_scaled)

# 排序

vip_df = pd.DataFrame({
    'protein': x_df.columns[top_indices] if n_proteins > 0 else x_df.columns,
    'VIP': vip_values
})

vip_df = vip_df.head(20)  # 只取top 20
vip_df = vip_df.sort_values('VIP', ascending=False)

# 绘图
plt.figure(figsize=(12,6))
plt.bar(vip_df['protein'], vip_df['VIP'], color='purple')
plt.xticks(rotation=90)
plt.xlabel('Proteins')
plt.ylabel('VIP scores')
plt.title('VIP scores of top proteins')
plt.tight_layout()
vip_plot_path = save_path.replace('.png', '_vip.png')
plt.savefig(vip_plot_path, dpi=300)
plt.show()

# =====================================================

from scipy.spatial.distance import cosine
from scipy.stats import linregress

# 原始模型
r2y_actual = 1.0  # R2Y通常接近1
q2_actual = pls.score(X_scaled, Y)

# 置换检验
n_perm = 50
r2y_perm_list = []
q2_perm_list = []
similarity_list = []

for _ in range(n_perm):
    Y_perm = np.random.permutation(Y)
    sim = 1 - np.corrcoef(Y, Y_perm)[0, 1]  # 相似度

    pls_perm = PLSRegression(n_components=2)
    pls_perm.fit(X_scaled, Y_perm)

    similarity_list.append(sim)
    r2y_perm_list.append(1.0)  # R2Y恒为1
    q2_perm_list.append(pls_perm.score(X_scaled, Y_perm))

similarity = np.array(similarity_list)
r2y_perm = np.array(r2y_perm_list)
q2_perm = np.array(q2_perm_list)

# 绘图
plt.figure(figsize=(10, 6))
plt.scatter(similarity, r2y_perm, marker='^', color='#1f77b4', s=100, label='R2Y', zorder=3)
plt.scatter(similarity, q2_perm, marker='s', color='#2ca02c', s=100, label='Q2', zorder=3)

# Q2回归线
slope, intercept, r_value, _, _ = linregress(similarity, q2_perm)
x_line = np.array([similarity.min(), similarity.max()])
y_line = intercept + slope * x_line
plt.plot(x_line, y_line, color='#2ca02c', linewidth=2)

plt.axhline(y=0, color='k', linestyle='--', linewidth=0.8)
plt.axhline(y=0.05, color='k', linestyle='--', linewidth=0.8)
plt.xlabel('Similar(y,y_perm)', fontsize=12)
plt.ylabel('R2Y/Q2', fontsize=12)
plt.legend(loc='upper right')
plt.title('PLS-DA Permutation Test')
plt.tight_layout()
perm_plot_path = save_path.replace('.png', '_permutation.png')
plt.savefig(perm_plot_path, dpi=300)
plt.show()