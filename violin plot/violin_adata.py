import time
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import numpy as np
from scipy import stats
import os


def load_data_by_groups(adata, genes, group_col='group', cell_type_col='predicted_labels',
                        layer=None, remove_zero_cells=True):
    """
    按分组加载数据

    Parameters:
    -----------
    remove_zero_cells : bool
        是否删除表达量为0的细胞（默认True）
        - True: 删除表达量为0的细胞（严格模式）
        - False: 保留表达量为0的细胞（宽松模式）

    Returns:
    --------
    dict : {gene: {group: {cell_type: expression_values}}}
    """

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

    # 构建嵌套字典
    data_dict = {gene: {} for gene in found_genes}

    for gene_idx, gene_name in enumerate(found_genes):
        gene_data = X[:, gene_idx]
        groups = adata.obs[group_col].values
        cell_types = adata.obs[cell_type_col].values

        for group in np.unique(groups):
            data_dict[gene_name][group] = {}
            group_mask = groups == group

            for cell_type in np.unique(cell_types[group_mask]):
                mask = (groups == group) & (cell_types == cell_type)
                values = gene_data[mask]

                # 根据开关决定是否删除表达量为0的细胞
                if remove_zero_cells:
                    values = values[values > 0]

                if len(values) > 0:
                    data_dict[gene_name][group][cell_type] = values

    return data_dict


def normalize_by_cell_count(data_dict, method='equal_sample', enable_resample=True):
    """
    按细胞数进行标准化（重采样）

    Parameters:
    -----------
    data_dict : dict
        原始数据字典
    method : str
        标准化方法：
        - 'equal_sample': 每个细胞类型采样相同数量的细胞
        - 'downsample_to_min': 每个细胞类型采样到最少的那个
        - 'upsample_to_max': 每个细胞类型采样到最多的那个（有重复）
    enable_resample : bool
        是否启用重采样（默认True）
        - True: 进行重采样标准化
        - False: 不进行重采样，保留原始细胞数

    Returns:
    --------
    dict : 标准化后的数据字典（或原始数据）
    """

    if not enable_resample:
        # 不进行重采样，直接返回原始数据
        return data_dict

    normalized_dict = {}

    for gene, groups_data in data_dict.items():
        # 找出最少和最多的细胞数
        all_cell_counts = []
        for group_data in groups_data.values():
            for cell_type_values in group_data.values():
                all_cell_counts.append(len(cell_type_values))

        if not all_cell_counts:
            normalized_dict[gene] = data_dict[gene]
            continue

        if method == 'equal_sample':
            # 使用最少的细胞数
            target_count = min(all_cell_counts)
        elif method == 'downsample_to_min':
            target_count = min(all_cell_counts)
        elif method == 'upsample_to_max':
            target_count = max(all_cell_counts)
        else:
            target_count = min(all_cell_counts)

        normalized_dict[gene] = {}
        for group, group_data in groups_data.items():
            normalized_dict[gene][group] = {}
            for cell_type, values in group_data.items():
                # 重采样
                if len(values) >= target_count:
                    # 下采样：随机选择
                    resampled = np.random.choice(values, size=target_count, replace=False)
                else:
                    # 上采样：有放回地采样
                    resampled = np.random.choice(values, size=target_count, replace=True)

                normalized_dict[gene][group][cell_type] = resampled

    return normalized_dict


def generate_colors(cell_types):
    """为每个细胞类型生成颜色"""
    color_palette = sns.color_palette("husl", n_colors=len(cell_types))
    color_dict = {ct: color_palette[i] for i, ct in enumerate(cell_types)}
    return color_dict


def get_global_xlim(normalized_data, genes_list, groups_list, cell_types_list):
    """
    计算全局的横坐标范围（基于所有基因、分组和细胞类型的表达值）

    Returns:
    --------
    tuple : (x_min, x_max)
    """
    all_values = []

    for gene in genes_list:
        if gene not in normalized_data:
            continue
        for group in groups_list:
            if group not in normalized_data[gene]:
                continue
            for cell_type in cell_types_list:
                if cell_type in normalized_data[gene][group]:
                    values = normalized_data[gene][group][cell_type]
                    all_values.extend(values)

    if not all_values:
        return (0, 1)

    all_values = np.array(all_values)
    x_min = np.min(all_values)
    x_max = np.max(all_values)

    # 添加5%的边距
    margin = (x_max - x_min) * 0.05
    return (x_min - margin, x_max + margin)


def diagnose_cell_counts(adata, genes, group_col='group', cell_type_col='predicted_labels',
                         layer=None):
    """
    诊断函数：查看原始数据中的细胞数分布
    """
    print(f"\n{'=' * 80}")
    print(f"诊断：原始数据中的细胞数分布")
    print(f"{'=' * 80}\n")

    # 获取数据
    if layer and layer in adata.layers:
        X = adata[:, genes].layers[layer]
    else:
        X = adata[:, genes].X

    # 转为密集数组
    if hasattr(X, 'toarray'):
        X = X.toarray()

    groups = adata.obs[group_col].values
    cell_types = adata.obs[cell_type_col].values

    found_genes = [g for g in genes if g in adata.var_names]

    # 统计原始细胞数（未删除0值）
    print("【模式1：保留所有细胞（包含表达量=0）】")
    print("-" * 80)
    for group in sorted(np.unique(groups)):
        print(f"\n{group} 组:")
        group_mask = groups == group
        group_total = 0
        for cell_type in sorted(np.unique(cell_types[group_mask])):
            mask = (groups == group) & (cell_types == cell_type)
            count = np.sum(mask)
            print(f"  {cell_type:20s}: {count:5d} 个细胞")
            group_total += count
        print(f"  {'合计':20s}: {group_total:5d} 个细胞")

    # 统计删除0值后的细胞数
    print(f"\n{'=' * 80}")
    print("【模式2：删除表达量=0的细胞】")
    print("-" * 80)

    for group in sorted(np.unique(groups)):
        print(f"\n{group} 组:")
        group_mask = groups == group

        group_nonzero_counts = {}
        for cell_type in sorted(np.unique(cell_types[group_mask])):
            mask = (groups == group) & (cell_types == cell_type)

            # 对该细胞类型的所有基因统计非零细胞
            nonzero_counts_per_gene = []
            for gene_idx, gene_name in enumerate(found_genes):
                gene_data = X[:, gene_idx]
                values = gene_data[mask]
                nonzero_count = np.sum(values > 0)
                nonzero_counts_per_gene.append(nonzero_count)

            # 显示每个基因的非零细胞数
            min_nonzero = min(nonzero_counts_per_gene) if nonzero_counts_per_gene else 0
            avg_nonzero = int(np.mean(nonzero_counts_per_gene)) if nonzero_counts_per_gene else 0
            max_nonzero = max(nonzero_counts_per_gene) if nonzero_counts_per_gene else 0

            group_nonzero_counts[cell_type] = {
                'min': min_nonzero,
                'avg': avg_nonzero,
                'max': max_nonzero,
                'original': np.sum(mask)
            }

            print(f"  {cell_type:20s}: 原始={group_nonzero_counts[cell_type]['original']:5d} | "
                  f"删除后 min={min_nonzero:5d} avg={avg_nonzero:5d} max={max_nonzero:5d}")

        group_total_original = sum(v['original'] for v in group_nonzero_counts.values())
        group_total_avg = sum(v['avg'] for v in group_nonzero_counts.values())
        print(f"  {'合计':20s}: 原始={group_total_original:5d} | 删除后平均={group_total_avg:5d}")

    print(f"\n{'=' * 80}\n")


def plot_single_gene(gene_name, normalized_data, groups_list,
                     cell_types_list, cell_type_colors,
                     global_xlim=None,
                     fontsize_title=20, fontsize_label=14,
                     fontsize_tick=12, point_size=20,
                     violin_width=0.5, alpha=0.4):
    """
    为单个基因绘制小提琴图并返回axes对象

    Returns:
    --------
    list : axes 列表
    """

    n_groups = len(groups_list)
    n_cell_types = len(cell_types_list)

    # 创建子图：1行 × n_groups列
    fig, axes = plt.subplots(1, n_groups, figsize=(6 * n_groups, 3 + n_cell_types * 0.6),
                             sharey=True)

    # 处理单个分组的情况
    if n_groups == 1:
        axes = [axes]

    # 为每个分组创建小提琴图
    for col_idx, group in enumerate(groups_list):
        ax = axes[col_idx]

        # 准备该分组中所有细胞类型的数据
        plot_data = []
        cell_type_labels = []
        plot_colors = []

        for cell_type in cell_types_list:
            if (group in normalized_data[gene_name] and
                    cell_type in normalized_data[gene_name][group]):

                values = normalized_data[gene_name][group][cell_type]
                if len(values) > 0:
                    plot_data.append(values)
                    cell_type_labels.append(cell_type)
                    plot_colors.append(cell_type_colors.get(cell_type, 'skyblue'))

        if not plot_data:
            ax.text(0.5, 0.5, f'No data', ha='center', va='center',
                    transform=ax.transAxes, fontsize=14)
            ax.set_title(group, fontsize=fontsize_title, fontweight='bold')
            continue

        # 绘制横向小提琴图
        positions = np.arange(len(plot_data))
        parts = ax.violinplot(plot_data, positions=positions, vert=False,
                              widths=violin_width, showmeans=False,
                              showextrema=False)

        # 设置小提琴图颜色
        for pc, color in zip(parts['bodies'], plot_colors):
            pc.set_facecolor(color)
            pc.set_alpha(alpha)
            pc.set_edgecolor('black')
            pc.set_linewidth(0.8)

        # 添加箱线图（显示中位数和四分位数）
        bp = ax.boxplot(plot_data, positions=positions, vert=False, widths=0.15,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(facecolor='white', alpha=0.8),
                        medianprops=dict(color='red', linewidth=2),
                        whiskerprops=dict(color='black', linewidth=1),
                        capprops=dict(color='black', linewidth=1))

        # 添加散点
        for i, (cell_type, values, color) in enumerate(zip(cell_type_labels, plot_data, plot_colors)):
            x_values = np.random.normal(i, 0.02, size=len(values))
            ax.scatter(values, x_values, color=color,
                       s=point_size, alpha=0.5, edgecolors='none')

        # 设置轴标签
        ax.set_yticks(positions)
        ax.set_yticklabels(cell_type_labels, fontsize=fontsize_tick)
        ax.set_xlabel('Expression Level', fontsize=fontsize_label, fontweight='bold')
        ax.set_title(f'{group}\n(n={len(plot_data[0])})',
                     fontsize=fontsize_title, fontweight='bold', pad=15)

        # 添加y轴标签（仅第一列）
        if col_idx == 0:
            ax.set_ylabel('Cell Type', fontsize=fontsize_label, fontweight='bold')

        # 美化
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        ax.tick_params(axis='x', labelsize=fontsize_tick)

        # 只在第一列显示y轴标签
        if col_idx > 0:
            ax.set_yticklabels([])

        # 设置统一的横坐标范围
        if global_xlim is not None:
            ax.set_xlim(global_xlim)

    # 添加细胞类型图例
    legend_elements = [Patch(facecolor=cell_type_colors.get(ct, 'skyblue'),
                             edgecolor='black', label=ct)
                       for ct in cell_types_list
                       if ct in cell_type_colors]

    fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.02),
               ncol=len(cell_types_list), fontsize=fontsize_label - 2, frameon=True)

    # 总标题
    fig.suptitle(f'{gene_name} Expression Across Cell Types and Groups',
                 fontsize=fontsize_title + 2, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0.05, 1, 0.96])

    return fig, axes


def plot_gene_by_celltype_with_normalization(gene_name, normalized_data, groups_list,
                                             cell_types_list, cell_type_colors,
                                             global_xlim=None,
                                             output_file=None,
                                             fontsize_title=20, fontsize_label=14,
                                             fontsize_tick=12, point_size=20,
                                             violin_width=0.5, alpha=0.4):
    """
    为单个基因绘制按细胞类型分列的小提琴图并保存
    """

    fig, axes = plot_single_gene(gene_name, normalized_data, groups_list,
                                 cell_types_list, cell_type_colors,
                                 global_xlim=global_xlim,
                                 fontsize_title=fontsize_title,
                                 fontsize_label=fontsize_label,
                                 fontsize_tick=fontsize_tick,
                                 point_size=point_size,
                                 violin_width=violin_width,
                                 alpha=alpha)

    if output_file is None:
        output_file = f'{gene_name}_celltype_comparison.png'

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f'Saved: {output_file}')
    plt.close()


def plot_combined_genes(genes_list, normalized_data, groups_list,
                        cell_types_list, cell_type_colors,
                        global_xlim=None,
                        output_file=None,
                        fontsize_title=16, fontsize_label=12,
                        fontsize_tick=10, point_size=15,
                        violin_width=0.5, alpha=0.4):
    """
    将多个基因的小提琴图拼接成一个大图

    Parameters:
    -----------
    genes_list : list
        基因列表
    output_file : str
        输出文件名
    """

    n_genes = len(genes_list)
    n_groups = len(groups_list)
    n_cell_types = len(cell_types_list)

    # 创建总的大图：n_genes行 × n_groups列
    fig, all_axes = plt.subplots(n_genes, n_groups,
                                 figsize=(6 * n_groups, 3 * n_genes + n_cell_types * 0.4),
                                 sharey='row')

    # 处理单个行或列的情况
    if n_genes == 1 and n_groups == 1:
        all_axes = [[all_axes]]
    elif n_genes == 1:
        all_axes = [all_axes]
    elif n_groups == 1:
        all_axes = [[ax] for ax in all_axes]

    # 为每个基因绘制小提琴图
    for gene_idx, gene_name in enumerate(genes_list):
        if gene_name not in normalized_data:
            continue

        for group_idx, group in enumerate(groups_list):
            ax = all_axes[gene_idx][group_idx]

            # 准备该分组中所有细胞类型的数据
            plot_data = []
            cell_type_labels = []
            plot_colors = []

            for cell_type in cell_types_list:
                if (group in normalized_data[gene_name] and
                        cell_type in normalized_data[gene_name][group]):

                    values = normalized_data[gene_name][group][cell_type]
                    if len(values) > 0:
                        plot_data.append(values)
                        cell_type_labels.append(cell_type)
                        plot_colors.append(cell_type_colors.get(cell_type, 'skyblue'))

            if not plot_data:
                ax.text(0.5, 0.5, f'No data', ha='center', va='center',
                        transform=ax.transAxes, fontsize=12)
                ax.set_title(f'{group}', fontsize=fontsize_title, fontweight='bold')
                continue

            # 绘制横向小提琴图
            positions = np.arange(len(plot_data))
            parts = ax.violinplot(plot_data, positions=positions, vert=False,
                                  widths=violin_width, showmeans=False,
                                  showextrema=False)

            # 设置小提琴图颜色
            for pc, color in zip(parts['bodies'], plot_colors):
                pc.set_facecolor(color)
                pc.set_alpha(alpha)
                pc.set_edgecolor('black')
                pc.set_linewidth(0.8)

            # 添加箱线图（显示中位数和四分位数）
            bp = ax.boxplot(plot_data, positions=positions, vert=False, widths=0.15,
                            patch_artist=True, showfliers=False,
                            boxprops=dict(facecolor='white', alpha=0.8),
                            medianprops=dict(color='red', linewidth=2),
                            whiskerprops=dict(color='black', linewidth=1),
                            capprops=dict(color='black', linewidth=1))

            # 添加散点
            for i, (cell_type, values, color) in enumerate(zip(cell_type_labels, plot_data, plot_colors)):
                x_values = np.random.normal(i, 0.02, size=len(values))
                ax.scatter(values, x_values, color=color,
                           s=point_size, alpha=0.5, edgecolors='none')

            # 设置轴标签
            ax.set_yticks(positions)
            ax.set_yticklabels(cell_type_labels, fontsize=fontsize_tick)
            ax.set_xlabel('Expression Level', fontsize=fontsize_label, fontweight='bold')

            # 标题：第一行显示分组名
            if gene_idx == 0:
                ax.set_title(f'{group}\n(n={len(plot_data[0])})',
                             fontsize=fontsize_title, fontweight='bold', pad=10)
            else:
                ax.set_title(f'{group}',
                             fontsize=fontsize_title, fontweight='bold', pad=10)

            # y轴标签：第一列显示基因名
            if group_idx == 0:
                ax.set_ylabel(f'{gene_name}\nCell Type', fontsize=fontsize_label, fontweight='bold')
            else:
                ax.set_yticklabels([])

            # 美化
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(axis='x', alpha=0.3, linestyle='--')
            ax.tick_params(axis='x', labelsize=fontsize_tick)

            # 设置统一的横坐标范围
            if global_xlim is not None:
                ax.set_xlim(global_xlim)

    # 添加细胞类型图例
    legend_elements = [Patch(facecolor=cell_type_colors.get(ct, 'skyblue'),
                             edgecolor='black', label=ct)
                       for ct in cell_types_list
                       if ct in cell_type_colors]

    fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.01),
               ncol=len(cell_types_list), fontsize=fontsize_label, frameon=True)

    # 总标题
    fig.suptitle(f'Gene Expression Across Cell Types and Groups',
                 fontsize=fontsize_title + 4, fontweight='bold', y=0.995)

    plt.tight_layout(rect=[0, 0.03, 1, 0.99])

    if output_file is None:
        output_file = 'combined_genes_celltype_comparison.png'

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f'Saved combined plot: {output_file}')
    plt.close()


def main(adata, genes, group_col='group', cell_type_col='predicted_labels',
         groups_list=['Y', 'M', 'O'], layer=None, output_dir=None,
         normalize_method='equal_sample',
         remove_zero_cells=True,
         enable_resample=True,
         fontsize_title=20, fontsize_label=14, fontsize_tick=12,
         point_size=20, violin_width=0.5, alpha=0.4):
    """
    主函数：为每个基因生成按细胞类型分列的对比图

    Parameters:
    -----------
    adata : ad.AnnData
        单细胞数据对象
    genes : list
        要绘制的基因列表
    group_col : str
        分组列名
    cell_type_col : str
        细胞类型列名
    groups_list : list
        分组列表，用于确定绘图顺序
    layer : str or None
        数据层名称
    output_dir : str
        输出文件目录
    normalize_method : str
        细胞数标准化方法：
        - 'equal_sample': 采样到最少细胞数（推荐）
        - 'downsample_to_min': 下采样到最少
        - 'upsample_to_max': 上采样到最多
    remove_zero_cells : bool
        是否删除表达量为0的细胞（默认True）
        - True: 删除表达量为0的细胞（严格模式）
        - False: 保留表达量为0的细胞（宽松模式）
    enable_resample : bool
        是否启用重采样（默认True）
        - True: 进行重采样标准化，使不同细胞类型的细胞数相同
        - False: 不进行重采样，保留原始细胞数
    fontsize_* : int
        各种字体大小参数
    point_size : int
        散点大小
    violin_width : float
        小提琴图宽度
    alpha : float
        透明度
    """

    start = time.time()

    # 验证列是否存在
    if group_col not in adata.obs.columns:
        available_cols = adata.obs.columns.tolist()
        raise ValueError(f"adata.obs 中不存在列 '{group_col}'。\n可用的列: {available_cols}")

    if cell_type_col not in adata.obs.columns:
        available_cols = adata.obs.columns.tolist()
        raise ValueError(f"adata.obs 中不存在列 '{cell_type_col}'。\n可用的列: {available_cols}")

    # 创建输出目录
    if output_dir is None:
        output_dir = './gene_celltype_comparison_plots'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"创建输出目录: {output_dir}\n")

    # 获取实际存在的分组和细胞类型
    available_groups = sorted(adata.obs[group_col].unique().tolist())
    available_cell_types = sorted(adata.obs[cell_type_col].unique().tolist())

    print(f"数据中的所有分组: {available_groups}")
    print(f"数据中的所有细胞类型: {available_cell_types}")

    # 过滤groups_list
    groups_to_plot = [g for g in groups_list if g in available_groups]
    print(f"将要绘制的分组: {groups_to_plot}")
    print(f"将要绘制的细胞类型: {available_cell_types}")
    print(f"细胞数标准化方法: {normalize_method}")

    # 构建删除细胞的描述文本
    remove_cells_desc = "删除了零表达细胞" if remove_zero_cells else "保留了所有细胞"
    resample_desc = "启用了重采样" if enable_resample else "未启用重采样"
    print(f"细胞处理: {remove_cells_desc}，{resample_desc}")
    print()

    # 为每个细胞类型生成颜色
    cell_type_colors = generate_colors(available_cell_types)

    # 加载数据
    print(f"正在加载数据...")
    data_dict = load_data_by_groups(adata, genes, group_col=group_col,
                                    cell_type_col=cell_type_col, layer=layer,
                                    remove_zero_cells=remove_zero_cells)

    # 对细胞数进行标准化（重采样）
    if enable_resample:
        print(f"正在对细胞数进行标准化 ({normalize_method})...")
        normalized_data = normalize_by_cell_count(data_dict, method=normalize_method,
                                                  enable_resample=True)
    else:
        print(f"跳过重采样，使用原始细胞数...")
        normalized_data = normalize_by_cell_count(data_dict, method=normalize_method,
                                                  enable_resample=False)

    # 计算全局的横坐标范围
    print(f"正在计算全局横坐标范围...")
    global_xlim = get_global_xlim(normalized_data, genes, groups_to_plot, available_cell_types)
    print(f"全局横坐标范围: [{global_xlim[0]:.4f}, {global_xlim[1]:.4f}]\n")

    # 为每个基因绘制图表
    print(f"\n{'=' * 80}")
    print(f"开始绘制基因表达量对比图")
    print(f"{'=' * 80}\n")

    valid_genes = []
    for i, gene in enumerate(genes):
        if gene not in normalized_data:
            print(f"[{i + 1}/{len(genes)}] ⚠ {gene}: 未找到\n")
            continue

        print(f"[{i + 1}/{len(genes)}] 正在绘制 {gene}")

        # 统计信息
        gene_has_data = False
        for group in groups_to_plot:
            if group in normalized_data[gene]:
                total_cells = sum(len(v) for v in normalized_data[gene][group].values() if len(v) > 0)
                if total_cells > 0:
                    print(f"    {group}: {total_cells} 个细胞")
                    gene_has_data = True

        if not gene_has_data:
            print(f"    ✗ 该基因在所有分组中没有数据，跳过\n")
            continue

        valid_genes.append(gene)
        output_file = os.path.join(output_dir, f'{gene}_celltype_comparison.png')

        try:
            plot_gene_by_celltype_with_normalization(
                gene, normalized_data, groups_to_plot, available_cell_types,
                cell_type_colors, global_xlim=global_xlim,
                output_file=output_file,
                fontsize_title=fontsize_title,
                fontsize_label=fontsize_label,
                fontsize_tick=fontsize_tick,
                point_size=point_size,
                violin_width=violin_width,
                alpha=alpha
            )
            print(f"    ✓ 生成成功\n")
        except Exception as e:
            print(f"    ✗ 生成失败: {e}\n")

    # 生成拼接图
    if valid_genes:
        print(f"\n{'=' * 80}")
        print(f"正在生成拼接图...")
        print(f"{'=' * 80}\n")

        combined_output = os.path.join(output_dir, 'combined_all_genes_celltype_comparison.png')

        try:
            plot_combined_genes(
                valid_genes, normalized_data, groups_to_plot, available_cell_types,
                cell_type_colors, global_xlim=global_xlim,
                output_file=combined_output,
                fontsize_title=fontsize_title - 4,
                fontsize_label=fontsize_label - 2,
                fontsize_tick=fontsize_tick - 2,
                point_size=point_size - 5,
                violin_width=violin_width,
                alpha=alpha
            )
            print(f"✓ 拼接图生成成功\n")
        except Exception as e:
            print(f"✗ 拼接图生成失败: {e}\n")

    end = time.time()
    print(f"{'=' * 80}")
    print(f"所有基因绘制完成！")
    print(f'耗时: {end - start:.2f} 秒')
    print(f"输出目录: {output_dir}")
    print(f"绘制的基因数: {len(valid_genes)}/{len(genes)}")
    print(f"全局横坐标范围: [{global_xlim[0]:.4f}, {global_xlim[1]:.4f}]")
    print(f"细胞处理: {remove_cells_desc}，{resample_desc}")
    print(f"{'=' * 80}")


# 使用示例
if __name__ == '__main__':
    import scanpy as sc
    from pathlib import Path

    # 加载单细胞数据
    path = Path(r"/dellstorage09/quj_lab/yanghang/leukocyte_single_cell/Human/human-B_cell/data")
    name = "after_Annotation.h5ad"
    adata = sc.read_h5ad(path / name)

    # 查看可用的列
    print("adata.obs 中的可用列:")
    print(adata.obs.columns.tolist())
    print()

    # ============================================================
    # 【第一步】诊断原始细胞数（看看到底有多少细胞）
    # ============================================================
    genes_to_check = ['HOPX', 'LITAF', 'PLEK', 'ZBTB32']
    diagnose_cell_counts(adata,
                         genes=genes_to_check,
                         group_col='group',
                         cell_type_col='cluster_cell_type')

    # ============================================================
    # 【第二步】绘制图表
    # ============================================================
    main(adata,
         genes=genes_to_check,
         group_col='group',
         cell_type_col='cluster_cell_type',
         groups_list=['Y','M','O'],
         output_dir='/dellstorage09/quj_lab/yanghang/leukocyte_single_cell/Human/human-B_cell/figure/gene_celltype_comparison_plots',
         normalize_method='upsample_to_max',
         remove_zero_cells=None,  # ⭐ 是否删除表达量为0的细胞
         enable_resample=False,  # ⭐ 是否进行重采样
         fontsize_title=20,
         fontsize_label=14,
         fontsize_tick=12,
         point_size=20,
         violin_width=0.5,
         alpha=0.4)