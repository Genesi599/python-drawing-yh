import random
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.font_manager import FontProperties

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
unique_number = len(data['cancer'].unique())
color_list = color.ncolors(unique_number, shade='light', alpha=0.7)
# 将颜色列表转换为元组列表
color_list = list(map(tuple, color_list))

random.shuffle(color_list)
final_list = []
for i in list(range(unique_number)):
    list1 = color_list[i:i + 1] * data['cancer'].value_counts(sort=False).iloc[i]

    final_list += list1

datalist = list(range(len(data['cancer'])))
# datalist.reverse()

# 设置全局字体为 Arial
plt.rcParams["font.sans-serif"] = ["Arial"]


plt.close('all')
fig, ax = plt.subplots()
fig.set_size_inches(17, 5)

for i in datalist:
    plt.bar(data1['cell-line'].iloc[i],
            data1['IC50'].iloc[i],
            width=0.7,
            label=data['cancer'].iloc[i],
            color=color.rgba_to_hex(final_list[i]),
            edgecolor='black',
            )

# 横坐标设置
move = 0.5
loc_list = np.arange(move, len(data1['cell-line']) + move, 1)
ax.set_xticks(loc_list)
ax.set_xticklabels(data1['cell-line'])
plt.xticks(fontsize=17, rotation=45, ha='right', fontweight='bold')
ax.tick_params(axis='x', length=0)  # 隐藏横坐标刻度线


# 设置对数坐标轴标签
plt.ylabel("IC$_{50}$ (μM)", fontsize=18, fontweight='bold')
plt.gca().xaxis.set_label_coords(0.5, 1.08)
plt.ylim((0.001, 100))
plt.yticks(fontsize=19, fontweight='bold')

# 对数坐标轴设置
plt.yscale("log")
formatter = ticker.ScalarFormatter()
formatter.set_scientific(False)
ax.yaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(StrMethodFormatter('{x:g}'))

# Set axis line thickness
thickness = 3  # Adjust this value to change the thickness
# ax.spines['top'].set_linewidth(thickness)
ax.spines['bottom'].set_linewidth(thickness)
ax.spines['left'].set_linewidth(thickness)
# ax.spines['right'].set_linewidth(thickness)

# 去掉右边和上边的边框
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# 合并相同图例
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
# by_label = OrderedDict(reversed(list(by_label.items())))  # 倒序排列

# 设置字体属性为粗体
font_properties = FontProperties(weight='bold', size=19)

# 创建没有边框的图例
plt.legend(
    by_label.values(), by_label.keys(),
    bbox_to_anchor=(1, 1),
    frameon=False,  # 去掉图例的边框
    prop=font_properties
)

# Add a horizontal dashed line at y=1
plt.axhline(y=1, color='gray', linestyle='--', linewidth=1.2)

# 创建矩形框

rectangle = plt.Rectangle((30.5,0), 9, 10, linewidth=2, linestyle='--', edgecolor='red', fill=False)

# 添加矩形框

ax.add_patch(rectangle)

# plt.show()
fig.savefig('fig.png', dpi=600, bbox_inches='tight')
plt.savefig('fig.svg', format='svg', bbox_inches='tight')

end = time.time()
print('Running time: %s Seconds' % (end - start))
