import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

# 读取CSV文件
df = pd.read_csv('data.csv', index_col=0)


# df = df[df.index.isin(class_df['WBC'])]

# 初始化一个空的DataFrame来存储结果
z_score_df = pd.DataFrame(index=df.index, columns=df.columns)

# 确保所有数据都是数值类型，并处理缺失值
z_score_df = z_score_df.apply(pd.to_numeric, errors='coerce').fillna(0)

# 全局设置字体属性
plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 使用 SimHei 字体以支持中文
plt.rcParams['axes.unicode_minus'] = False  # 确保负号正常显示

# 创建从蓝色到红色的颜色映射
cmap = LinearSegmentedColormap.from_list('blue_red', ['white', 'orange'])

# 创建一个新的图形和轴对象
fig, ax = plt.subplots(dpi=600, figsize=(8, 6))

# 绘制热图
sns.heatmap(data=df,
            cmap=cmap,
            vmin=-0, vmax=100,
            xticklabels=True,
            cbar_kws={'orientation': 'horizontal', 'shrink': 0.2, 'pad': 0.1},
            ax=ax)

# 设置横坐标刻度标签的大小和旋转角度
plt.xticks(rotation=45, fontsize=10, ha='right', va='top', fontweight='bold')
plt.tick_params(axis='x', length=0)  # 隐藏横坐标刻度线

# 将纵坐标标签放到右边并水平显示
ax.yaxis.tick_right()               # 将刻度移至右边
ax.yaxis.set_label_position("right")  # 将标签移至右边
plt.yticks(rotation=0, fontsize=10, fontweight='bold')
plt.tick_params(axis='y', length=0)  # 隐藏横坐标刻度线
# ax.set_yticklabels([])  # 仅隐藏刻度值


# 将图例放到左上方
cbar = ax.collections[0].colorbar
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.xaxis.set_label_position('top')
cbar.ax.set_position([0.1, 0.8, 0.2, 0.15])  # [left, bottom, width, height]
# 设置 Colorbar 刻度值的文字大小
cbar.ax.tick_params(labelsize=15)  # 设置刻度文字大小为 10
for label in cbar.ax.get_xticklabels():  # 设置刻度文字加粗
    label.set_fontweight('bold')
# 添加颜色条注释
cbar.ax.text(0.5, -0.5, 'viable cell', ha='center', va='center', fontsize=18,  fontweight='bold',transform=cbar.ax.transAxes)

# 添加竖条颜色框
height = z_score_df.shape[0]

# 设置横轴标签和纵轴标签的字体大小和格式
plt.xlabel(f"dBET1 (nM)", fontsize=18, fontweight="bold")
plt.ylabel(f"CC-90009 (nM)", fontsize=18, fontweight="bold")


plt.savefig('heatmap_with_colorbar.png', format='png', bbox_inches='tight')