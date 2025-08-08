import random
import matplotlib.pyplot as plt
import pandas as pd
from Drawing import color
import time
from collections import OrderedDict
from matplotlib.ticker import ScalarFormatter
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
    list1 = color_list[i:i + 1] * data['cancer'].value_counts(sort=False)[i]
    final_list += list1

# mpl.rcParams["font.sans-serif"] = ['Microsoft YaHei']

plt.close('all')
plt.subplots(len(data['cancer'].unique()), 1,
             sharex='all',
             figsize=(4, 8), dpi=500)

datalist = list(range(len(data['cancer'])))
datalist.reverse()
for i in datalist:
    ax = plt.subplot()
    plt.barh(data1['cell-line'].iloc[i],
             data1['IC50'].iloc[i],
             height=0.65,
             label=data['cancer'].iloc[i],
             color=final_list[i],
             edgecolor='black')
    # 横坐标设置
    ax.xaxis.set_ticks_position('top')
    plt.xscale("log")
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
    # 添加图例
    plt.legend()

# 合并相同图例
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
by_label = OrderedDict(reversed(list(by_label.items())))  # 倒序排列
plt.legend(by_label.values(), by_label.keys(), loc='upper right', bbox_to_anchor=(1.5, 1))
# 设置x轴标签
plt.xlabel("IC$_{50}$ (μM)")
plt.gca().xaxis.set_label_coords(0.5, 1.08)
plt.xlim((0.01, 100))
plt.yticks(fontsize=9)

# 保存图片
plt.savefig('fig.png', dpi=500, bbox_inches='tight')

end = time.time()
print('Running time: %s Seconds' % (end - start))
