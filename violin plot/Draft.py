import time
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import random

start = time.time()
# 生成颜色列表
color_list = sns.color_palette("Set1", n_colors=10)
random.shuffle(color_list)
# sns.set(style="whitegrid")


# 导入数据，从excel中

df = pd.read_excel('2-8_sgRNA_Summary.xlsx')
pathway = pd.read_excel('pathway.xlsx')
pathway_list = pathway['pathway'].unique().tolist()

one_pathway = pathway[pathway['pathway'] == pathway_list[1]]['GENE'].tolist()

# 创建一个包含子图的图表
fig, axs = plt.subplots(len(one_pathway), 1, figsize=(8, 1*len(one_pathway)), sharex=True)

print(one_pathway)
for i in list(range(len(one_pathway))):
    print(i)
    print(one_pathway[i])
    df2 = df[df['Gene'] == one_pathway[i]]['LFC']
    # axs[i] = violin
    sns.violinplot(x=df["LFC"], legend='', orient='y', cut=0, linewidth=0.3, bw_method='silverman',
                   bw_adjust=1, inner="quart", width=0.95, fill=True, alpha=0.3, color='red', ax=axs[i])
    # 加入散点图
    sns.stripplot(x=df2, color='red', jitter=0.1, ax=axs[i], size=6)
    # 文本框
    text_box = plt.text(0.95, 0.6, one_pathway[i], transform=axs[i].transAxes, fontsize=18, ha='right', va='center')
    axs[i].yaxis.set_visible(False)
    axs[i].set_xlabel('The log fold change of sgRNA', fontsize=20)
    for label in axs[i].get_xticklabels():
        label.set_fontsize(15)
    if i < len(one_pathway) - 1:
        axs[i].axis('off')
    # 去除上下左右的边框（默认该函数会取出右上的边框）
    sns.despine(left=True)  # bottom=True

# 调整子图之间的间距
plt.subplots_adjust(wspace=0, hspace=0.1)
plt.savefig('violin.png', dpi=300)
end = time.time()
print('Running time: %s Seconds' % (end - start))
