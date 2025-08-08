import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# 模拟数据
np.random.seed(42)
data = np.random.rand(10, 10)
z_score_df = pd.DataFrame(data)

# 自定义的行顺序（例如反转原顺序）
custom_row_order = z_score_df.index[::-1]  # 行顺序调整为从最后一行到第一行
sorted_df = z_score_df.loc[custom_row_order, :]  # 重排数据框的行

# 绘制热图
cluster = sns.clustermap(
    data=sorted_df,
    row_cluster=False,  # 禁用行聚类，固定顺序
    col_cluster=False,  # 禁用列聚类，也固定列顺序
    figsize=(16, 9),
    cmap='vlag',
    xticklabels=True,
    yticklabels=True
)

plt.show()