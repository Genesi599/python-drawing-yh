import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math
from adjustText import adjust_text

# 导入数据
# Specify the sheet you want to read, e.g., 'Sheet1'
sheet_name = 'Sheet1'
specific_names = [
    'ICAM-1',
    'CXCL16'
]

xlabel_name = 'R'
ylabel_name = 'pvalue'

title_label = 'Correlation'

# Read the specified sheet from the Excel file
pd_df = pd.read_excel('data.xlsx', sheet_name=sheet_name)
# print(pd_df)

x_name = 'spearman_corr'
y_name = 'p_value'
IDname = 'ID'

result = pd.DataFrame()

x = pd_df[x_name]
result['x'] = x
y = pd_df[y_name]
result['y'] = y

result['y'] = -np.log10(result['y'])

# 设置显著性阈值
x_threshold = 1
y_threshold = 2

result = result.dropna()

print(result)

# 添加列, 计算与原点距离
result = result.sort_index()
result = pd.concat([result, pd_df[IDname]], axis=1)
result['average'] = (result['x'] ** 2 + result['y'] ** 2) ** 0.5
result = result.sort_values(by='average', ascending=False)

# 分组为up, normal, down
result['group-edge'] = 'dimgrey'
result['group-face'] = '#d9f5ff'
# result.loc[(result.x > x_threshold) & (result.y > y_threshold), 'group'] = 'dimgrey'  # x=-+x_threshold直接截断
# result.loc[(result.x < -x_threshold) & (result.y > y_threshold), 'group'] = 'dimgrey'  # x=-+x_threshold直接截断
# result.loc[result.y < y_threshold, 'group'] = 'dimgrey'
result.loc[result.x > 0, 'group-face'] = 'pink'
# result.loc[result.y > 0, 'group-edge'] = 'pink'
for i in specific_names:
    if i in result[IDname].tolist():
        result.loc[result[IDname] == i, 'group'] = 'tab:red'  # x=-+x_threshold直接截断

# 确定坐标轴显示范围
xmax = math.ceil(result['x'].max())
xmin = math.floor(result['x'].min())
ymax = math.ceil(result['y'].abs().max()) + 1
ymin = 0

# 全局设置字体属性
plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 使用 SimHei 字体以支持中文
plt.rcParams['axes.unicode_minus'] = False  # 确保负号正常显示

# 绘制散点图
fig = plt.figure(figsize=plt.figaspect(9 / 16))  # 确定fig比例（h/w）
ax = fig.add_subplot()
ax.set(xlim=(xmin - 0.1, xmax + 0.1), ylim=(ymin, ymax + 0.1), title='')
ax.scatter(result['x'], result['y'], s=20, marker='o',
           facecolors=result['group-face'], edgecolors=result['group-edge'])
for a in specific_names:
    if a in result[IDname].tolist():
        ax.scatter(
            result[result[IDname] == a]['x'],
            result[result[IDname] == a]['y'],
            c='red',
            s=20,  # 可以增大点的大小
            zorder=10  # 更高的层级
        )
for a in specific_names:
    if a in result[IDname].tolist():
        ax.scatter(
            result[result[IDname] == a]['x'],
            result[result[IDname] == a]['y'],
            c='red',
            s=20,  # 可以增大点的大小
            zorder=10  # 更高的层级
        )
ax.set_xlabel(xlabel_name, fontweight='bold', fontsize=20)
ax.set_ylabel(ylabel_name, fontweight='bold', fontsize=20)
ax.spines['right'].set_visible(False)  # 去掉右边框
ax.spines['top'].set_visible(False)  # 去掉上边框

# 标注点
# 根据特定条件选择标注点
texts = []
for i in range(len(result)):
    # 添加条件，例如：只标注x值大于某个阈值的点
    for a in specific_names:
        if a in result[IDname].tolist():
            if result.iloc[i][IDname] == a:
                text = ax.text(
                    result.iloc[i]['x'] + 0.2,
                    result.iloc[i]['y'],
                    result.iloc[i][IDname],
                    fontsize=15,
                    color='#430000',
                    style="italic",
                    weight="bold",
                    verticalalignment='center',
                    horizontalalignment='left'
                )
                texts.append(text)

# adjust_text(texts,
#             # arrowprops=dict(arrowstyle='->', lw=1.5, color='skyblue'),
#             force_text=(0.5, 15),      # 增加文字间的排斥力
#             )

# 水平和竖直线
ax.vlines(0, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
# ax.vlines(x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
ax.hlines(0, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖水平线

# ax.set_xticks(range(xmin, xmax, 5))  # 设置x轴刻度起点和步长
ax.set_yticks(range(ymin + 1, ymax, 1))  # 设置y轴刻度起点和步长
# 设置刻度字体大小
plt.xticks(fontsize=20, fontweight='bold')  # x轴刻度字体大小
plt.yticks(fontsize=20, fontweight='bold')  # y轴刻度字体大小

ax.set_title(title_label, fontweight='bold', fontsize=20)
fig.savefig(f"{sheet_name}.png", dpi=600, bbox_inches='tight')  # 保存为eps矢量图
