import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math
from adjustText import adjust_text
from pathlib import Path

path = Path(r"C:\Users\yh109\Nutstore\1\Institute of Zoology\protein\FC.csv")
listpath = Path(r"C:\Users\yh109\Documents\GitHub\enrichment_analysis-yh\data\annolist.txt")
with open(listpath, 'r', encoding='utf-8') as f:
    loaded = [line.rstrip('\n') for line in f]
annolist = loaded
print(annolist)


pd_df = pd.read_csv(path, index_col=0, encoding='utf-8')
pd_df = pd_df[:500]
pd_df['TS13'] = np.log10(pd_df['TS13'])
pd_df['protein'] = pd_df.index
pd_df = pd_df.sort_values(by='log2FC', ascending=True)
pd_df = pd_df.reset_index(drop=True)
pd_df['ID'] = pd_df.index
print(pd_df)



result = pd.DataFrame()

posx = pd_df[pd_df['log2FC'] > 0]['ID']
negx = pd_df[pd_df['log2FC'] < 0]['ID']
x = pd.concat([posx, negx], axis=0)
result['x'] = x

posy = pd_df[pd_df['log2FC'] > 0]['TS13']
negy = pd_df[pd_df['log2FC'] < 0]['TS13']
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
xmax = math.ceil(result['x'].max()*1.1)
ymin = math.floor(result['y'].min())
ymax = math.ceil(result['y'].max())

# 绘制散点图
fig = plt.figure(figsize=plt.figaspect(9 / 9))  # 确定fig比例（h/w）
ax = fig.add_subplot()
ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')
ax.scatter(result['x'], result['y'], s=5, c=result['group'])
ax.set_ylabel(r'log$_\boldsymbol{10}$(Abundance)', fontweight='bold', fontsize=12)
ax.set_xlabel('rank', fontweight='bold', fontsize=12)
ax.spines['right'].set_visible(False)  # 去掉右边框
ax.spines['top'].set_visible(False)  # 去掉上边框

# 标注点
result = result[result['x'] > 0]  # 只取x＜0

# texts = [ax.text(result.iloc[i]['x'] + 0.1, result.iloc[i]['y'], result.iloc[i]['GENE'], fontsize=10, color="black",
#                  style="italic", weight="bold", verticalalignment='center', horizontalalignment='left') for i in
#          range(10) if result.iloc[i]['y'] > y_threshold and result.iloc[i]['x'] > 0]

texts = [ax.text(result.iloc[i]['x'] + 0.1, result.iloc[i]['y'], result.iloc[i]['GENE'], fontsize=9, color="black",
                 style="italic", weight="bold", verticalalignment='center', horizontalalignment='left') for i in
         range(len(result)) if result.iloc[i]['GENE'] in annolist]

adjust_text(texts, arrowprops=dict(arrowstyle='->', lw=0.9, color='green'),
            force_explode = (-10, 2),
            iter_lim = 10000,
            # explode_radius = 10,
            avoid_self=True)


# 水平和竖直线
ax.vlines(-x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
ax.vlines(x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
ax.hlines(y_threshold, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖水平线

ax.set_xticks(range(xmin, xmax, 200))  # 设置x轴刻度起点和步长
ax.set_yticks(range(ymin, ymax, 1))  # 设置y轴刻度起点和步长

fig.savefig('Volcano plot.png', dpi=600)  # 保存为eps矢量图
