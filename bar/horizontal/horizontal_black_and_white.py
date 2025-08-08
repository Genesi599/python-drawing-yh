import random
import matplotlib.pyplot as plt
import pandas as pd
from Drawing import color
import time
import numpy as np
from collections import OrderedDict
import matplotlib.ticker as ticker
from matplotlib.ticker import StrMethodFormatter

start = time.time()

# 读取文件
data = pd.read_excel('IC50 summary.xlsx', sheet_name='Sheet2')
data1 = data.set_index(["cancer"])

# 生成颜色列表
color_list = list(map(lambda x: color.color(tuple(x)), color.ncolors(len(data['cancer'].unique()))))
random.shuffle(color_list)
final_list = []
for i in list(range(len(data['cancer'].unique()))):
    list1 = color_list[i:i + 1] * data['cancer'].value_counts(sort=False).iloc[i]

    final_list += list1

datalist = list(range(len(data['cancer'])))
# datalist.reverse()

# 设置全局字体为 Arial
plt.rcParams["font.sans-serif"] = ["Arial"]


plt.close('all')
fig, ax = plt.subplots()
fig.set_size_inches(8, 4)

# # 绘制柱状图时应用不同的颜色
# for i in datalist:
#     plt.bar(data1['cell-line'].iloc[i],
#             data1['IC50'].iloc[i],
#             width=0.7,
#             label=data['cancer'].iloc[i],
#             color=final_list[i],
#             edgecolor='black')

# 定义填充图案列表
hatch_patterns = [
    '///', '//O', 'ooo', 'xxo', '+++', '----', '//oo', '..//', 'OO', 'XX', 'OX'
]

random.shuffle(hatch_patterns)
final_hatch_patterns = hatch_patterns[:len(data['cancer'].unique())]
final_list = []
for i in list(range(len(data['cancer'].unique()))):
    list1 = hatch_patterns[i:i + 1] * data['cancer'].value_counts(sort=False).iloc[i]
    final_list += list1

# 绘制柱状图时应用不同的图案
for i in datalist:
    plt.bar(data1['cell-line'].iloc[i],
            data1['IC50'].iloc[i],
            width=0.7,
            label=data['cancer'].iloc[i],
            color='white',  # 使用白色填充
            edgecolor='black',
            hatch=final_list[i])  # 应用图案

# 横坐标设置
move = -0.5
loc_list = np.arange(move, len(data1['cell-line']) + move, 1)
ax.set_xticks(loc_list)
ax.set_xticklabels(data1['cell-line'])
plt.xticks(fontsize=8, rotation=315, ha='left')
ax.tick_params(axis='x', length=0)  # 隐藏横坐标刻度线


# 添加图例
plt.legend(fontsize=10)

# 设置对数坐标轴标签
plt.ylabel("IC$_{50}$ (μM)")
plt.gca().xaxis.set_label_coords(0.5, 1.08)
plt.ylim((0.001, 100))
plt.yticks(fontsize=9)

# 对数坐标轴设置
plt.yscale("log")
formatter = ticker.ScalarFormatter()
formatter.set_scientific(False)
ax.yaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(StrMethodFormatter('{x:g}'))

# 去掉右边和上边的边框
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# 合并相同图例
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
# by_label = OrderedDict(reversed(list(by_label.items())))  # 倒序排列

# 创建没有边框的图例
plt.legend(
    by_label.values(), by_label.keys(),
    bbox_to_anchor=(1, 1),
    frameon=False  # 去掉图例的边框
)


# plt.show()
fig.savefig('fig.png', dpi=600, bbox_inches='tight')

end = time.time()
print('Running time: %s Seconds' % (end - start))
