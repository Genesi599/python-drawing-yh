import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math
from adjustText import adjust_text

# 导入数据
pd_df = pd.read_excel('O_vs_M.xlsx')
name = 'Old vs Middle'
result = pd.DataFrame()
print(pd_df)

posx = pd_df[pd_df['logFC'] > 0]['logFC']
negx = pd_df[pd_df['logFC'] < 0]['logFC']
x = pd.concat([posx, negx], axis=0)
result['x'] = x

posy = pd_df[pd_df['logFC'] > 0]['P.Value']
negy = pd_df[pd_df['logFC'] < 0]['P.Value']
y = pd.concat([posy, negy], axis=0)
result['y'] = y
result['y'] = -np.log10(result['y'])

# 设置显著性阈值
x_threshold = 0
y_threshold = 1.30103

# 分组为up, normal, down
result['group'] = 'black'
result.loc[(result.x > x_threshold) & (result.y > y_threshold), 'group'] = 'tab:red'  # x=-+x_threshold直接截断
result.loc[(result.x < -x_threshold) & (result.y > y_threshold), 'group'] = 'tab:blue'  # x=-+x_threshold直接截断
result.loc[result.y < y_threshold, 'group'] = 'dimgrey'  # 阈值以下点为灰色

# 添加列
result = result.sort_index()
result = pd.concat([result, pd_df['proteinID']], axis=1)
result['average'] = (result['x'] ** 2 + result['y'] ** 2) ** 0.5
result = result.sort_values(by='average', ascending=False)
print(result.head())
print(result.iloc[0]['proteinID'])

# 确定坐标轴显示范围
xmin = math.floor(result['x'].min())
xmax = math.ceil(result['x'].max())
ymin = 0
ymax = math.ceil(result['y'].max())

# 全局设置字体属性
plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 使用 SimHei 字体以支持中文

# 绘制散点图
fig = plt.figure(figsize=plt.figaspect(9 / 16))  # 确定fig比例（h/w）
fig.suptitle(name, fontsize=16, fontweight='bold')  # 添加标题
ax = fig.add_subplot()
ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')
ax.scatter(result['x'], result['y'], s=20, c=result['group'])
ax.set_ylabel('-Log10(P.Value)', fontweight='bold', fontsize=16)
ax.set_xlabel('Log2 Foldchange', fontweight='bold', fontsize=16)
ax.spines['right'].set_visible(False)  # 去掉右边框
ax.spines['top'].set_visible(False)  # 去掉上边框

# 标注点
# result = result[result['x'] < 0]  # 只取x＜0

texts = [
    ax.text(
        result.iloc[i]['x'] + 0.01,
        result.iloc[i]['y'],
        result.iloc[i]['proteinID'],
        fontsize=8,
        color="black",
        style="italic",
        weight="bold",
        verticalalignment='center',
        horizontalalignment='left'
    )
    for i in range(100)
    if result.iloc[i]['y'] > y_threshold and (result.iloc[i]['x'] > x_threshold or result.iloc[i]['x'] < -x_threshold)
]

adjust_text(texts, arrowprops=dict(arrowstyle='->', lw=0.5, color='green'))


# 水平和竖直线
# ax.vlines(-x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
# ax.vlines(x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
ax.vlines(0, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
# ax.hlines(y_threshold, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖水平线

ax.set_xticks(range(xmin, xmax, 2))  # 设置x轴刻度起点和步长
ax.set_yticks(range(ymin, ymax, 1))  # 设置y轴刻度起点和步长

fig.savefig('Volcano plot.png', dpi=600)  # 保存为eps矢量图
