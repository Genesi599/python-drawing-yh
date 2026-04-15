import random
import matplotlib.pyplot as plt
import pandas as pd
from Drawing import color
from pylab import mpl

# 读取文件
data = pd.read_excel('IC50 summary.xlsx', sheet_name='Sheet2')

data1 = data.set_index(["cancer"])
print(data1['cell-line'].iloc[1])

# 生成颜色列表
color_list = list(map(lambda x: color.color(tuple(x)), color.ncolors(len(data['cancer'].unique()))))
random.shuffle(color_list)
final_list = []
for i in list(range(len(data['cancer'].unique()))):
    list1 = color_list[i:i + 1] * data['cancer'].value_counts(sort=False)[i]
    print(list1)
    final_list += list1

# df = CRISPRLibData[(CRISPRLibData["cancer"] == '胰腺癌')]

mpl.rcParams["font.sans-serif"] = ['Microsoft YaHei']

plt.close('all')


plt.subplots(len(data['cancer'].unique()), 1,
             sharex='all',
             figsize=(4, 8), dpi=500)

for i in list(range(len(data['cancer'].unique()))):
    ax = plt.subplot(len(data['cancer'].unique()), 1, i + 1)
    plt.barh(data1.at[data['cancer'].unique()[i], 'cell-line'],
             data1.at[data['cancer'].unique()[i], 'IC50'],
             height=0.65,
             label=data['cancer'].unique()[i],
             color=color_list[i],
             edgecolor='black')
    ax.invert_yaxis()
    ax.xaxis.set_ticks_position('top')
    plt.xscale("log")
    plt.legend()

plt.xlabel("IC50 (μM)")
plt.gca().xaxis.set_label_coords(0.5, 1.08)
plt.xlim((0.01, 100))
plt.yticks(fontsize=9)


# 添加颜色标签


# cmap = ListedColormap(colors=final_list)
# Leg = Legend(ax, loc='lower right', labels='1', handles=)
# ax.add_artist(Leg)

# cbar = plt.colorbar(cmap, orientation='vertical')
# cbar.ax.set_yticklabels(CRISPRLibData['cancer'].unique())


# 保存图片
plt.savefig('fig.png', dpi=500, bbox_inches='tight')
