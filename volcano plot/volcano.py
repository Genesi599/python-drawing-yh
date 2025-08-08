import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math
from adjustText import adjust_text

# 导入数据
pd_df = pd.read_excel('day7_Gene_Summary.xlsx')
result = pd.DataFrame()

posx = pd_df[pd_df['pos|lfc'] > 0]['pos|lfc']
negx = pd_df[pd_df['neg|lfc'] < 0]['neg|lfc']
x = pd.concat([posx, negx], axis=0)
result['x'] = x

posy = pd_df[pd_df['pos|lfc'] > 0]['pos|p-value']
negy = pd_df[pd_df['neg|lfc'] < 0]['neg|p-value']
y = pd.concat([posy, negy], axis=0)
result['y'] = y
result['y'] = -np.log10(result['y'])

# 设置显著性阈值
x_threshold = 1
y_threshold = 2

# 分组为up, normal, down
result['group'] = 'black'
result.loc[(result.x > x_threshold) & (result.y > y_threshold), 'group'] = 'tab:red'  # x=-+x_threshold直接截断
result.loc[(result.x < -x_threshold) & (result.y > y_threshold), 'group'] = 'tab:blue'  # x=-+x_threshold直接截断
result.loc[result.y < y_threshold, 'group'] = 'dimgrey'  # 阈值以下点为灰色

# 添加列
result = result.sort_index()
result = pd.concat([result, pd_df['id']], axis=1)
result['average'] = (result['x'] ** 2 + result['y'] ** 2) ** 0.5
result = result.sort_values(by='average', ascending=False)
print(result.head())
print(result.iloc[0]['id'])

# 确定坐标轴显示范围
xmin = math.floor(result['x'].min())
xmax = math.ceil(result['x'].max())
ymin = 0
ymax = math.ceil(result['y'].max())

# 绘制散点图
fig = plt.figure(figsize=plt.figaspect(9 / 16))  # 确定fig比例（h/w）
ax = fig.add_subplot()
ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')
ax.scatter(result['x'], result['y'], s=2, c=result['group'])
ax.set_ylabel('-Log10(p-value)', fontweight='bold')
ax.set_xlabel('log fold change', fontweight='bold')
ax.spines['right'].set_visible(False)  # 去掉右边框
ax.spines['top'].set_visible(False)  # 去掉上边框

# 标注点
result = result[result['x'] < 0]  # 只取x＜0

texts = [ax.text(result.iloc[i]['x'] + 0.1, result.iloc[i]['y'], result.iloc[i]['id'], fontsize=7, color="black",
                 style="italic", weight="bold", verticalalignment='center', horizontalalignment='left') for i in
         range(10) if result.iloc[i]['y'] > y_threshold and result.iloc[i]['x'] < 0]

adjust_text(texts, arrowprops=dict(arrowstyle='->', lw=0.5, color='green'))


# 水平和竖直线
ax.vlines(-x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
ax.vlines(x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
ax.hlines(y_threshold, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖水平线

ax.set_xticks(range(xmin, xmax, 2))  # 设置x轴刻度起点和步长
ax.set_yticks(range(ymin, ymax, 1))  # 设置y轴刻度起点和步长

fig.savefig('Volcano plot.png', dpi=300)  # 保存为eps矢量图
