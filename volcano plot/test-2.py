import numpy as np  # Scientific computing
import pandas as pd  # Data analysis
import matplotlib.pyplot as plt  # Plotting
from adjustText import adjust_text

def create_volcano_plot(
        input_file,
        output_file='Volcano_plot.png',
        x_threshold=0.5,
        y_threshold=0.5,  # 更新为适当的阈值
        lfc_col='log2_fold_change',  # 仍使用 log2_fold_change
        p_col='adjusted_p_value',  # 使用 adjusted_p_value
        id_col='gene_name',  # 仍使用 gene_name
        selected_genes=None):

    # 读数据
    df = pd.read_csv(input_file)

    # 只保留需要的列
    vol = df[[id_col, lfc_col, 'Mean_Thymus', p_col]].copy()

    # 找到当前最小的非零 p 值
    min_non_zero_pval = vol[vol[p_col] > 0][p_col].min()

    # 替换 p 值为0的情况
    vol[p_col] = vol[p_col].replace(0, min_non_zero_pval)

    # 计算 y 值
    vol['y'] = vol['Mean_Thymus']  # 纵轴为 Mean_Thymus

    # 分组 (normal/up/down)
    vol['group'] = 'grey'  # 默认组
    up_mask = (vol[lfc_col] >= x_threshold) & (vol['y'] >= y_threshold)
    dn_mask = (vol[lfc_col] <= -x_threshold) & (vol['y'] >= y_threshold)
    vol.loc[up_mask, 'group'] = 'tab:red'  # up
    vol.loc[dn_mask, 'group'] = 'tab:blue'  # down

    # 排序 (按距离)
    vol['dist'] = vol[lfc_col] ** 2 + vol['y'] ** 2
    vol = vol.sort_values('dist', ascending=False)

    # 画图
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(vol[lfc_col], vol['y'], s=10, c=vol['group'], alpha=0.8)  # 点更小

    # 阈值线
    ax.axvline(-x_threshold, ls='--', color='grey', lw=1)
    ax.axvline(x_threshold, ls='--', color='grey', lw=1)
    ax.axhline(y_threshold, ls='--', color='grey', lw=1)

    # 设置标签
    ax.set_xlabel('Log2 Fold Change', fontweight='bold', fontsize=12)
    ax.set_ylabel('Mean Thymus Expression', fontweight='bold', fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # 标注选定的基因
    if selected_genes is not None:
        texts = []
        gene_xy = {}  # 记录每个基因的真实数据点坐标
        for i, gene in enumerate(selected_genes):
            if gene in vol[id_col].values:
                row = vol[vol[id_col] == gene].iloc[0]
                x = row[lfc_col]
                y = row['y']
                gene_xy[gene] = (x, y)  # 存真坐标

                # 初始偏移（跟你上一版一致）
                txt = ax.text(x + 0.15 + i*0.08, y + 0.15 + i*0.08, gene,
                              fontsize=9, style='italic', weight='bold',
                              ha='left', va='bottom')
                texts.append(txt)

        # 1. 先让文字按原有策略微调（不碰箭头）
        adjust_text(texts,
                    x=vol[lfc_col].values, y=vol['y'].values,
                    ax=ax,
                    force_points=(0.5, 0.5),
                    force_text=(0.3, 0.3),
                    expand_points=(0, 0),
                    ha='left', va='bottom',
                    only_move={'points': 'y', 'text': 'xy'},
                    force_expand=(0, 0.5),
                    arrowprops=None)  # ← 先不画箭头

        # 2. 再按“最终文字位置→原始数据点”重画箭头
        for gene, (x, y) in gene_xy.items():
            for t in texts:
                if t.get_text() == gene:
                    ax.annotate('',
                                xy=(x, y),                  # 真实数据点
                                xytext=t.get_position(),    # 文字最终位置
                                ha='left', va='bottom',
                                arrowprops=dict(
                                    arrowstyle='->',
                                    color='gray',        # 灰色
                                    lw=0.5,
                                    ls='--',             # 虚线
                                    connectionstyle='arc3,rad=0'))
                    break
    fig.tight_layout()

    # 保存 PNG
    fig.savefig(output_file, dpi=300)
    # 额外保存 PDF
    fig.savefig(output_file.replace('.png', '.pdf'), dpi=300, bbox_inches='tight')

    plt.close()
    print("DataFrame head:\n", df.head())
    print("Log2 Fold Change range:", df[lfc_col].min(), df[lfc_col].max())
    print("Adjusted p-value range:", df[p_col].min(), df[p_col].max())
    # 打印横坐标值的范围
    print("Log2 Fold Change range:", vol[lfc_col].min(), vol[lfc_col].max())


# 示例使用
create_volcano_plot(
    input_file=r"D:\Projects\Thymus_Aging\mouse_atlas\thymus_specific_genes_mapped.csv",
    output_file='Volcano_plot.png',
    x_threshold=1,
    y_threshold=1,  # 根据需要设置合理的阈值
    selected_genes=[
        "STMN1",
        "RGS10",
        "RBM38",
        "FERMT3",
        "BZW2",
        "PRKCB",
        "DOK2",
        "FYB1",
        "CSK",
        "LDHB",
        "PRPS2",
        "GMFG",
        "THY1",
        "CD27"
    ]
)