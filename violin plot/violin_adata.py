import time
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import sparse


def load_data(adata, genes, layer=None, group_col='group', gene_col=None):
    """直接在var_names中找基因"""

    # 找出存在的基因
    found_genes = [g for g in genes if g in adata.var_names]

    if not found_genes:
        raise ValueError("指定的基因都未找到")

    # 提示未找到的基因
    missing_genes = [g for g in genes if g not in adata.var_names]
    if missing_genes:
        print(f"Warning: 基因 {missing_genes} 未找到，已跳过")

    # 获取数据
    if layer and layer in adata.layers:
        X = adata[:, found_genes].layers[layer]
    else:
        X = adata[:, found_genes].X

    # 转为密集数组
    if hasattr(X, 'toarray'):
        X = X.toarray()

    # 构建数据框（长格式）
    df_list = []
    for gene_idx, gene_name in enumerate(found_genes):
        gene_data = X[:, gene_idx]
        for cell_idx, expr_value in enumerate(gene_data):
            df_list.append({
                'Gene': gene_name,
                'Group': adata.obs[group_col].values[cell_idx],
                'Expression': expr_value
            })

    df = pd.DataFrame(df_list)
    return df


def generate_colors(groups):
    """生成颜色列表"""
    color_palette = sns.color_palette("Set2", n_colors=len(groups))
    color_dict = {group: color_palette[i] for i, group in enumerate(groups)}
    return color_dict


def plot_single_gene(ax, df, gene_name, cell_types, color_dict,
                     show_ylabel=True,
                     fontsize_gene=16, fontsize_xlabel=12, fontsize_ylabel=13,
                     fontsize_tick=11, point_size=3):
    """
    绘制单个基因在不同细胞类型中的小提琴图和散点图

    Parameters:
    -----------
    fontsize_gene : int
        基因名称的字体大小
    fontsize_xlabel : int
        x轴标签的字体大小
    fontsize_ylabel : int
        y轴标签的字体大小
    fontsize_tick : int
        刻度标签的字体大小
    point_size : int
        散点的大小
    """

    for i, cell_type in enumerate(cell_types):
        group_data = df[(df['Gene'] == gene_name) & (df['CellType'] == cell_type)]['Expression']

        if len(group_data) == 0:
            continue

        # 绘制小提琴图
        sns.violinplot(x=[i] * len(group_data), y=group_data,
                       cut=0, linewidth=0.3,
                       bw_method='silverman', bw_adjust=1, inner="quart",
                       width=0.7,
                       fill=True, alpha=0.3, color=color_dict[cell_type], ax=ax)

        # 加入散点图（调整大小）
        sns.stripplot(x=[i] * len(group_data), y=group_data,
                      color=color_dict[cell_type], jitter=0.1, ax=ax, size=point_size)

    # 添加基因名称文本框
    ax.text(0.02, 0.98, gene_name, transform=ax.transAxes,
            fontsize=fontsize_gene, ha='left', va='top', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # 设置x轴标签为细胞类型名称
    ax.set_xticks(range(len(cell_types)))
    ax.set_xticklabels(cell_types, fontsize=fontsize_xlabel, rotation=45, ha='right')
    ax.set_xlabel('')

    # 设置y轴标签
    if show_ylabel:
        ax.set_ylabel('Expression Level', fontsize=fontsize_ylabel)
    else:
        ax.set_ylabel('')

    # 美化
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    # 设置刻度标签字体大小
    for label in ax.get_yticklabels():
        label.set_fontsize(fontsize_tick)


def plot_genes(df, genes, cell_types, color_dict, output_file=None,
               fontsize_gene=16, fontsize_xlabel=12, fontsize_ylabel=13,
               fontsize_tick=11, point_size=3):
    """
    绘制多个基因的小提琴图

    Parameters:
    -----------
    fontsize_gene : int
        基因名称的字体大小（默认16）
    fontsize_xlabel : int
        x轴标签的字体大小（默认12）
    fontsize_ylabel : int
        y轴标签的字体大小（默认13）
    fontsize_tick : int
        刻度标签的字体大小（默认11）
    point_size : int
        散点的大小（默认3）
    """

    print(f"\n绘制基因: {genes}")

    # 检查基因是否存在
    available_genes = [gene for gene in genes if gene in df['Gene'].unique()]
    if not available_genes:
        print(f"Error: 没有找到指定的基因")
        return False

    genes = available_genes

    # 创建包含子图的图表
    n_genes = len(genes)
    fig, axs = plt.subplots(n_genes, 1,
                            figsize=(10, 3.5 * n_genes),
                            sharex=True)

    # 处理单个基因的情况（axs不是数组）
    if n_genes == 1:
        axs = [axs]

    for i, gene in enumerate(genes):
        is_first = (i == 0)
        plot_single_gene(axs[i], df, gene, cell_types, color_dict,
                         show_ylabel=is_first,
                         fontsize_gene=fontsize_gene,
                         fontsize_xlabel=fontsize_xlabel,
                         fontsize_ylabel=fontsize_ylabel,
                         fontsize_tick=fontsize_tick,
                         point_size=point_size)

    # 调整子图之间的间距
    plt.tight_layout()

    # 保存图表
    if output_file is None:
        output_file = 'violin_plot.png'

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f'Saved: {output_file}')
    plt.close()

    return True


def main(adata, genes, group_col='group', cell_type_col='predicted_labels',
         layer=None, gene_col='symbol', output_dir=None,
         fontsize_gene=16, fontsize_xlabel=12, fontsize_ylabel=13,
         fontsize_tick=11, point_size=3):
    """
    主函数：执行整个分析流程，为每个分组生成独立的图片

    Parameters:
    -----------
    adata : ad.AnnData
        单细胞数据对象
    genes : list
        要绘制的基因列表
        例如：['TP53', 'BRCA1', 'EGFR']
    group_col : str
        adata.obs中的分组列名，默认为'group'
    cell_type_col : str
        adata.obs中的细胞类型列名，默认为'predicted_labels'
    layer : str or None
        adata中的层名称。如果为None，使用adata.raw.X或adata.X
    gene_col : str
        用于基因标识的列名，默认为'symbol'
    output_dir : str, optional
        输出文件目录
    fontsize_gene : int
        基因名称的字体大小（默认16）
    fontsize_xlabel : int
        x轴标签的字体大小（默认12）
    fontsize_ylabel : int
        y轴标签的字体大小（默认13）
    fontsize_tick : int
        刻度标签的字体大小（默认11）
    point_size : int
        散点的大小（默认3，推荐范围1-5）
    """
    import os
    start = time.time()

    # 验证分组列是否存在
    if group_col not in adata.obs.columns:
        available_cols = adata.obs.columns.tolist()
        raise ValueError(f"adata.obs 中不存在列 '{group_col}'。\n可用的列: {available_cols}")

    # 验证细胞类型列是否存在
    if cell_type_col not in adata.obs.columns:
        available_cols = adata.obs.columns.tolist()
        raise ValueError(f"adata.obs 中不存在列 '{cell_type_col}'。\n可用的列: {available_cols}")

    # 获取可用的分组值
    available_groups = sorted(adata.obs[group_col].unique().tolist())
    print(f"列 '{group_col}' 中的所有值: {available_groups}")

    # 获取可用的细胞类型
    available_cell_types = sorted(adata.obs[cell_type_col].unique().tolist())
    print(f"列 '{cell_type_col}' 中的所有值: {available_cell_types}")

    # 创建输出目录
    if output_dir is None:
        output_dir = './violin_plots'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"创建输出目录: {output_dir}")

    # 生成颜色字典
    color_dict = generate_colors(available_cell_types)

    # 为每个分组生成图片
    for group in available_groups:
        print(f"\n{'=' * 60}")
        print(f"处理分组: {group}")
        print(f"{'=' * 60}")

        # 子集化adata，只保留当前分组的数据
        adata_subset = adata[adata.obs[group_col] == group]
        print(f"该分组包含 {adata_subset.n_obs} 个细胞")

        # 加载该分组的数据
        print(f"正在加载数据...")
        df = load_data(adata_subset, genes=genes, layer=layer,
                       group_col=cell_type_col, gene_col=gene_col)

        # 确定该分组中实际存在的细胞类型
        cell_types_in_group = sorted(df['Group'].unique().tolist())
        print(f"该分组中的细胞类型: {cell_types_in_group}")

        # 重命名df中的Group列为CellType
        df.rename(columns={'Group': 'CellType'}, inplace=True)

        # 绘制基因图表
        genes_to_plot = [g for g in genes if g in df['Gene'].unique()]

        if not genes_to_plot:
            print(f"Warning: 该分组中未找到任何指定的基因")
            continue

        output_file = os.path.join(output_dir, f'violin_plot_{group}.png')
        success = plot_genes(df, genes_to_plot, cell_types_in_group,
                             color_dict, output_file=output_file,
                             fontsize_gene=fontsize_gene,
                             fontsize_xlabel=fontsize_xlabel,
                             fontsize_ylabel=fontsize_ylabel,
                             fontsize_tick=fontsize_tick,
                             point_size=point_size)

        if success:
            print(f"✓ {group} 的图片已生成")
        else:
            print(f"✗ {group} 的绘图失败")

    end = time.time()
    print(f"\n{'=' * 60}")
    print(f"所有分组处理完成！")
    print(f'Running time: {end - start:.2f} Seconds')
    print(f"输出目录: {output_dir}")


# 使用示例
if __name__ == '__main__':
    import anndata as ad
    from pathlib import Path
    import scanpy as sc

    # 加载单细胞数据
    path = Path(r"D:\leukocyte_single_cell\Monkey\monkey-B_cell\subluster\data")
    name = "after_Annotation.h5ad"
    adata = sc.read_h5ad(path / name)

    # 查看可用的列
    print("adata.obs 中的可用列:")
    print(adata.obs.columns.tolist())
    print()

    # 为每个分组生成独立的小提琴图
    # 调整字体和散点大小
    main(adata,
         genes=['HOPX', 'LITAF', 'PLEK', 'ZBTB32'],
         group_col='group',
         cell_type_col='cluster_cell_type',
         output_dir='./violin_plots_by_group',
         fontsize_gene=18,  # 基因名称字体大小（推荐16-20）
         fontsize_xlabel=13,  # x轴标签字体大小（推荐11-14）
         fontsize_ylabel=14,  # y轴标签字体大小（推荐12-15）
         fontsize_tick=12,  # 刻度标签字体大小（推荐10-13）
         point_size=1)  # 散点大小（推荐1-4，越小越不显眼）