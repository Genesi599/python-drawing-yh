import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math
from adjustText import adjust_text
from pathlib import Path

path = Path(r"C:\Users\yh109\OneDrive\桌面\FC.csv")
pd_df = pd.read_csv(path, index_col=0, encoding='utf-8')
pd_df = pd_df[:1000]
pd_df['protein'] = pd_df.index
pd_df = pd_df.sort_values(by='log2FC', ascending=True)
pd_df = pd_df.reset_index(drop=True)
pd_df['ID'] = pd_df.index


result = pd.DataFrame()

posx = pd_df[pd_df['log2FC'] > 0]['ID']
negx = pd_df[pd_df['log2FC'] < 0]['ID']
x = pd.concat([posx, negx], axis=0)
result['x'] = x

posy = pd_df[pd_df['log2FC'] > 0]['log2FC']
negy = pd_df[pd_df['log2FC'] < 0]['log2FC']
y = pd.concat([posy, negy], axis=0)
result['y'] = y
print(result)

# 设置显著性阈值
x_threshold = 0
y_threshold = math.floor(result['y'].min())

# 分组为up, normal, down
result['group'] = 'black'
result.loc[(result.x > x_threshold) & (result.y > y_threshold), 'group'] = 'tab:red'  # x=-+x_threshold直接截断
result.loc[(result.x < -x_threshold) & (result.y > y_threshold), 'group'] = 'tab:blue'  # x=-+x_threshold直接截断
result.loc[result.y < y_threshold, 'group'] = 'dimgrey'  # 阈值以下点为灰色

# 添加列
result = result.sort_index()
result = pd.concat([result, pd_df['GENE']], axis=1)
result['average'] = (result['x'] ** 2 + result['y'] ** 2) ** 0.5
result = result.sort_values(by='average', ascending=False)


# 确定坐标轴显示范围
# xmin = math.floor(result['x'].min())
xmin = 0
xmax = math.ceil(result['x'].max())
ymin = math.floor(result['y'].min())
ymax = math.ceil(result['y'].max())

# 绘制散点图
fig = plt.figure(figsize=plt.figaspect(9 / 16))  # 确定fig比例（h/w）
ax = fig.add_subplot()
ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')
ax.scatter(result['x'], result['y'], s=2, c=result['group'])
ax.set_ylabel('log2 fold change', fontweight='bold')
ax.set_xlabel('rank', fontweight='bold')
ax.spines['right'].set_visible(False)  # 去掉右边框
ax.spines['top'].set_visible(False)  # 去掉上边框

# 标注点
result = result[result['x'] > 0]  # 只取x＜0

texts = [ax.text(result.iloc[i]['x'] + 0.1, result.iloc[i]['y'], result.iloc[i]['GENE'], fontsize=10, color="black",
                 style="italic", weight="bold", verticalalignment='center', horizontalalignment='left') for i in
         range(10) if result.iloc[i]['y'] > y_threshold and result.iloc[i]['x'] > 0]

adjust_text(texts, arrowprops=dict(arrowstyle='->', lw=0.5, color='green'))


# 水平和竖直线
ax.vlines(-x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
ax.vlines(x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
ax.hlines(y_threshold, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖水平线

ax.set_xticks(range(xmin, xmax, 200))  # 设置x轴刻度起点和步长
ax.set_yticks(range(ymin, ymax, 1))  # 设置y轴刻度起点和步长

fig.savefig('Volcano plot.png', dpi=300)  # 保存为eps矢量图
