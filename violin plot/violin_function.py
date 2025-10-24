import time
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import random


def load_data(excel_file, pathway_file):
    """
    加载数据文件

    Parameters:
    -----------
    excel_file : str
        sgRNA数据文件路径
    pathway_file : str
        通路数据文件路径

    Returns:
    --------
    df : pd.DataFrame
        sgRNA数据
    pathway_list : list
        通路列表
    """
    df = pd.read_excel(excel_file)
    pathway = pd.read_excel(pathway_file)
    pathway_list = pathway['pathway'].unique().tolist()
    print(pathway_list)
    return df, pathway, pathway_list


def generate_colors(pathway_list):
    """
    生成颜色列表

    Parameters:
    -----------
    pathway_list : list
        通路列表

    Returns:
    --------
    color_list : list
        颜色列表
    """
    color_list = sns.color_palette("dark", n_colors=len(pathway_list))
    return color_list


def plot_single_gene(ax, df, gene_name, color, show_xlabel=True):
    """
    绘制单个基因的小提琴图和散点图

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        子图对象
    df : pd.DataFrame
        数据框
    gene_name : str
        基因名称
    color : tuple
        颜色
    show_xlabel : bool
        是否显示x轴标签
    """
    df2 = df[df['Gene'] == gene_name]['LFC']

    # 绘制小提琴图
    sns.violinplot(x=df["LFC"], legend='', orient='y', cut=0, linewidth=0.3,
                   bw_method='silverman', bw_adjust=1, inner="quart", width=0.95,
                   fill=True, alpha=0.3, color=color, ax=ax)

    # 加入散点图
    sns.stripplot(x=df2, color=color, jitter=0.1, ax=ax, size=8)

    # 添加文本框
    text_box = plt.text(0.95, 0.7, gene_name, transform=ax.transAxes,
                        fontsize=25, ha='right', va='center')

    # 设置轴属性
    ax.yaxis.set_visible(False)
    ax.spines['bottom'].set_linewidth(2)

    # 设置x轴标签
    if show_xlabel:
        ax.set_xlabel('The log fold change of sgRNA', fontsize=25)
    else:
        ax.set_xlabel('')

    for label in ax.get_xticklabels():
        label.set_fontsize(20)


def plot_pathway(df, pathway, pathway_list, color_list, pathway_index, output_file=None):
    """
    绘制单个通路的所有基因图表

    Parameters:
    -----------
    df : pd.DataFrame
        sgRNA数据
    pathway : pd.DataFrame
        通路数据
    pathway_list : list
        通路列表
    color_list : list
        颜色列表
    pathway_index : int
        当前通路索引
    output_file : str, optional
        输出文件名，不指定则使用默认名称
    """
    pathway_name = pathway_list[pathway_index]
    one_pathway = pathway[pathway['pathway'] == pathway_name]['GENE'].tolist()
    print(f"Processing pathway {pathway_index}: {one_pathway}")

    # 创建包含子图的图表
    fig, axs = plt.subplots(len(one_pathway), 1, figsize=(8, 1.1 * len(one_pathway)),
                            sharex=True)

    # 处理单个基因的情况（axs不是数组）
    if len(one_pathway) == 1:
        axs = [axs]

    for i in range(len(one_pathway)):
        is_last = (i == len(one_pathway) - 1)
        plot_single_gene(axs[i], df, one_pathway[i], color_list[pathway_index],
                         show_xlabel=is_last)

        # 隐藏非最后一行的轴
        if i < len(one_pathway) - 1:
            axs[i].axis('off')

    # 去除上下左右的边框
    sns.despine(left=True)

    # 调整子图之间的间距
    plt.subplots_adjust(wspace=0, hspace=0.1)

    # 添加纵坐标标签
    yaxis_ax = fig.add_axes([0.1, 0.15, 0, 0.7])
    yaxis_ax.set_yticks([])
    yaxis_ax.set_xticks([])
    yaxis_ax.set_ylabel(pathway_name, fontsize=25)
    yaxis_ax.spines['left'].set_linewidth(2)

    # 保存图表
    if output_file is None:
        output_file = f'violin-{pathway_index}.png'

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f'Saved: {output_file}')
    plt.close()


def main(excel_file='day7_sgRNA_Summary.xlsx', pathway_file='pathway.xlsx'):
    """
    主函数：执行整个分析流程

    Parameters:
    -----------
    excel_file : str
        sgRNA数据文件路径
    pathway_file : str
        通路数据文件路径
    """
    start = time.time()

    # 加载数据
    df, pathway, pathway_list = load_data(excel_file, pathway_file)

    # 生成颜色列表
    color_list = generate_colors(pathway_list)

    # 绘制每个通路的图表
    for a in range(len(pathway_list)):
        plot_pathway(df, pathway, pathway_list, color_list, a)

    end = time.time()
    print(f'Running time: {end - start:.2f} Seconds')


# 使用示例
if __name__ == '__main__':
    main()