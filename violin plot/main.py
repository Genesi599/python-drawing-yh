import time
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import random

start = time.time()

# random.shuffle(color_list)


# 导入数据，从excel中

df = pd.read_excel('day7_sgRNA_Summary.xlsx')

# df = pd.concat([df1, df2, df3])
pathway = pd.read_excel('pathway.xlsx')
pathway_list = pathway['pathway'].unique().tolist()
print(pathway_list)

# 生成颜色列表
color_list = sns.color_palette("dark", n_colors=len(pathway_list))

for a in list(range(len(pathway_list))):

    one_pathway = pathway[pathway['pathway'] == pathway_list[a]]['GENE'].tolist()
    print(one_pathway)

    # 创建一个包含子图的图表
    fig, axs = plt.subplots(len(one_pathway), 1, figsize=(8, 1.1 * len(one_pathway)), sharex=True)

    for i in list(range(len(one_pathway))):
        print(i)
        df2 = df[df['Gene'] == one_pathway[i]]['LFC']
        # axs[i] = violin
        sns.violinplot(x=df["LFC"], legend='', orient='y', cut=0, linewidth=0.3, bw_method='silverman',
                       bw_adjust=1, inner="quart", width=0.95, fill=True, alpha=0.3, color=color_list[a], ax=axs[i])
        # 加入散点图
        sns.stripplot(x=df2, color=color_list[a], jitter=0.1, ax=axs[i], size=8)
        # 文本框
        text_box = plt.text(0.95, 0.7, one_pathway[i], transform=axs[i].transAxes, fontsize=25, ha='right', va='center')
        axs[i].yaxis.set_visible(False)
        axs[i].set_xlabel('The log fold change of sgRNA', fontsize=25)
        axs[i].spines['bottom'].set_linewidth(2)
        for label in axs[i].get_xticklabels():
            label.set_fontsize(20)
        if i < len(one_pathway) - 1:
            axs[i].axis('off')

        # 去除上下左右的边框（默认该函数会取出右上的边框）
        sns.despine(left=True)  # bottom=True

    # 调整子图之间的间距
    plt.subplots_adjust(wspace=0, hspace=0.1)

    # 纵坐标
    yaxis_ax = fig.add_axes([0.1, 0.15, 0, 0.7])
    yaxis_ax.set_yticks([])
    yaxis_ax.set_xticks([])
    yaxis_ax.set_ylabel(pathway_list[a], fontsize=25)  # labelpad用于调整标签与轴的距离
    yaxis_ax.spines['left'].set_linewidth(2)
    # 保存
    file_name = 'violin' + "-" + str(a) + '.png'
    plt.savefig(file_name, dpi=300)
    print('pathway-'+str(a))
end = time.time()
print('Running time: %s Seconds' % (end - start))
