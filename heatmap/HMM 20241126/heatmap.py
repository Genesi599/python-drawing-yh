import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

# 读取CSV文件
df = pd.read_csv('data-1.csv', index_col=0)
# class_df = pd.read_excel('classify.xlsx')

# print(class_df)

# 去掉每个元素中的逗号，然后转换为数值类型
# df = df.applymap(lambda x: pd.to_numeric(str(x).replace(',', ''), errors='coerce'))
# print(df)

# df = df[df.index.isin(class_df['WBC'])]

# 初始化一个空的DataFrame来存储结果
z_score_df = pd.DataFrame(index=df.index, columns=df.columns)



# 对每一行计算Z分数并缩放
for index, row in df.iterrows():
    # 计算每一行的均值和标准差
    mean = row.mean()
    std = row.std()

    # 计算Z分数
    z_score = (row - mean) / std

    # 将结果存储到新的DataFrame中
    z_score_df.loc[index] = z_score

# 确保所有数据都是数值类型，并处理缺失值
z_score_df = z_score_df.apply(pd.to_numeric, errors='coerce').fillna(0)

# 转置DataFrame
# z_score_df = z_score_df.transpose()

print("Z-Score DataFrame:")
print(z_score_df)

# 全局设置字体属性
plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 使用 SimHei 字体以支持中文
plt.rcParams['axes.unicode_minus'] = False  # 确保负号正常显示

# 创建从蓝色到红色的颜色映射
cmap = LinearSegmentedColormap.from_list('blue_red', ['blue', 'white', 'red'])

# 创建一个新的图形和轴对象
fig, ax = plt.subplots(dpi=600, figsize=(8, 20))

# 绘制热图
sns.heatmap(data=z_score_df,
            cmap=cmap,
            vmin=-4, vmax=4,
            xticklabels=True,
            yticklabels=True,
            cbar_kws={'orientation': 'horizontal', 'shrink': 0.2, 'pad': 0.1},
            ax=ax)

# 设置横坐标刻度标签的大小和旋转角度
plt.xticks(rotation=90, fontsize=10, ha='right', va='top', fontweight='bold')
plt.tick_params(axis='x', length=0)  # 隐藏横坐标刻度线

# 将纵坐标标签放到右边并水平显示
ax.yaxis.tick_left()
plt.yticks(rotation=0, fontsize=10, fontweight='bold')
# plt.tick_params(axis='y', length=0)  # 隐藏横坐标刻度线
# ax.set_yticklabels([])  # 仅隐藏刻度值

# ax.yaxis.set_visible(False)  # 隐藏 y 轴

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
cbar.ax.text(0.5, -0.5, 'Column Z-Score', ha='center', va='center', fontsize=18,  fontweight='bold',transform=cbar.ax.transAxes)

# 添加竖条颜色框
height = z_score_df.shape[0]

# 在热图左侧添加颜色条
xposition = -2.5
width1 = 0.6
fontsize1 = 22
# ax.add_patch(Rectangle((xposition, 0), width1, height* 7.7/23, color='black', transform=ax.transData, clip_on=False))
# ax.add_patch(Rectangle((xposition, height/2.8), width1, height* 6.7/23, color='black', transform=ax.transData, clip_on=False))
# ax.add_patch(Rectangle((xposition, 2*height/3), width1, height* 7.7/23, color='black', transform=ax.transData, clip_on=False))
#
# # Add text on the left of each color bar
# ax.text(xposition, height/6, 'Young', va='center', ha='right', fontsize=fontsize1, transform=ax.transData, rotation=90, fontweight='bold')
# ax.text(xposition, height/2, 'Middle', va='center', ha='right', fontsize=fontsize1, transform=ax.transData, rotation=90, fontweight='bold')
# ax.text(xposition, 5*height/6, 'Old', va='center', ha='right', fontsize=fontsize1, transform=ax.transData, rotation=90, fontweight='bold')

plt.savefig('heatmap_with_colorbar.png', format='png', bbox_inches='tight')