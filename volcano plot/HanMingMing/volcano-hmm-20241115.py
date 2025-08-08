import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math
from adjustText import adjust_text

# 导入数据
# Specify the sheet you want to read, e.g., 'Sheet1'
sheet_name = 'young  vs middle '
selected_name = [
    # 'hsa-miR-98-5p',
    'hsa-miR-8085',
    # 'hsa-miR-1285-3p',
    'hsa-miR-545-3p',
    'hsa-miR-548e-5p',
    'hsa-miR-500b-3p'
]
xlabel_name = 'log₂(o⁸G IP) in middle'

# Read the specified sheet from the Excel file
pd_df = pd.read_excel('human samples.xlsx', sheet_name=sheet_name)
# print(pd_df)

x_name = 'log2(IP)'
y_name = 'Relative Enrichment 1'

result = pd.DataFrame()

x = pd_df[x_name]
result['x'] = x
y = pd_df[y_name]
result['y'] = y

# result['y'] = -np.log10(result['y'])

# 设置显著性阈值
x_threshold = 1
y_threshold = 2

result = result.dropna()
result = result[result['x'] >= 0]

print(result)


def log_scale(series, base=10):
    """
    保留原始数据正负号的对数缩放

    参数:
    series: 输入的pandas Series
    base: 对数的底数（默认为10）

    返回:
    对数缩放后的Series
    """

    def safe_log(x, base):
        """安全的对数计算"""
        if x == 0:
            return 0
        return np.sign(x) * np.log(abs(x)) / np.log(base)

    return series.apply(lambda x: safe_log(x, base))



# Z分数标准化
result['y'] = log_scale(result['y'])

# 添加列, 计算与原点距离
result = result.sort_index()
result = pd.concat([result, pd_df['miRNAID']], axis=1)
result['average'] = (result['x'] ** 2 + result['y'] ** 2) ** 0.5
result = result.sort_values(by='average', ascending=False)


# 分组为up, normal, down
result['group'] = 'grey'
for i in selected_name:
    if i in result['miRNAID'].tolist():
        result.loc[result['miRNAID'] == i, 'group'] = 'tab:red'  # x=-+x_threshold直接截断


# 确定坐标轴显示范围
xmin = 0
xmax = math.ceil(result['x'].max())

ymax = math.ceil(result['y'].abs().max())+1
ymin = -ymax

# 全局设置字体属性
plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 使用 SimHei 字体以支持中文
plt.rcParams['axes.unicode_minus'] = False  # 确保负号正常显示

# 绘制散点图
fig = plt.figure(figsize=plt.figaspect(9 / 16))  # 确定fig比例（h/w）
ax = fig.add_subplot()
ax.set(xlim=(xmin - 0.1, xmax + 0.1), ylim=(ymin, ymax+0.1), title='')
ax.scatter(result['x'], result['y'], s=20, c=result['group'])
for a in selected_name:
    if a in result['miRNAID'].tolist():
        ax.scatter(
            result[result['miRNAID'] == a]['x'],
            result[result['miRNAID'] == a]['y'],
            c='red',
            s=20,  # 可以增大点的大小
            zorder=10  # 更高的层级
            )
ax.set_xlabel(xlabel_name, fontweight='bold', fontsize=20)
ax.set_ylabel('Relative o⁸G enrichment', fontweight='bold', fontsize=20)
ax.spines['right'].set_visible(False)  # 去掉右边框
ax.spines['top'].set_visible(False)  # 去掉上边框



# 标注点
# 根据特定条件选择标注点
texts = []
for i in range(len(result)):
    # 添加条件，例如：只标注x值大于某个阈值的点
    for a in selected_name:
        if a in result['miRNAID'].tolist():
            if result.iloc[i]['miRNAID'] == a:
                text = ax.text(
                    result.iloc[i]['x'],
                    result.iloc[i]['y']+0.2,
                    result.iloc[i]['miRNAID'],
                    fontsize=15,
                    color="red",
                    style="italic",
                    weight="bold",
                    verticalalignment='bottom',
                    horizontalalignment='left'
                )
                texts.append(text)

adjust_text(texts,
            # arrowprops=dict(arrowstyle='->', lw=1.5, color='skyblue'),
            force_text=(0.5, 15),      # 增加文字间的排斥力
            )

# 水平和竖直线
# ax.vlines(0, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
# ax.vlines(x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
ax.hlines(0, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖水平线

ax.set_xticks(range(xmin, xmax, 1))  # 设置x轴刻度起点和步长
ax.set_yticks(range(ymin+1, ymax, 1))  # 设置y轴刻度起点和步长
# 设置刻度字体大小
plt.xticks(fontsize=20, fontweight='bold')  # x轴刻度字体大小
plt.yticks(fontsize=20, fontweight='bold')  # y轴刻度字体大小

fig.savefig(f"{sheet_name}.png", dpi=600, bbox_inches='tight')  # 保存为eps矢量图
