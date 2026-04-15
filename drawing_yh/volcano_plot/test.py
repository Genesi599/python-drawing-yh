from math import log10
import numpy as np  # Scientific computing
import pandas as pd  # Data analysis
import matplotlib.pyplot as plt  # Plotting
from adjustText import adjust_text


def create_volcano_plot(
        input_file,
        output_file='Volcano_plot.png',
        x_threshold=0.5,
        y_threshold=-np.log10(0.05),
        lfc_col='r',
        p_col='pval_adj',
        id_col='Gene.name',
        selected_genes=None):
    # 读数据
    df = pd.read_csv(input_file)

    # 只保留需要的列
    vol = df[[id_col, lfc_col, p_col]].copy()

    # 找到当前最小的非零 p 值
    min_non_zero_pval = vol[vol[p_col] > 0][p_col].min()

    # 替换 p 值为0的情况
    vol[p_col] = vol[p_col].replace(0, min_non_zero_pval)

    vol['y'] = -np.log10(vol[p_col])

    # 分组 (normal/up/down)
    vol['group'] = 'black'  # 默认组
    up_mask = (vol[lfc_col] >= x_threshold) & (vol['y'] >= y_threshold)
    dn_mask = (vol[lfc_col] <= -x_threshold) & (vol['y'] >= y_threshold)
    vol.loc[up_mask, 'group'] = 'tab:red'  # up
    vol.loc[dn_mask, 'group'] = 'tab:blue'  # down

    # 排序 (按距离)
    vol['dist'] = vol[lfc_col] ** 2 + vol['y'] ** 2
    vol = vol.sort_values('dist', ascending=False)

    # 画图
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(vol[lfc_col], vol['y'], s=30, c=vol['group'], alpha=0.8)

    # 阈值线
    ax.axvline(-x_threshold, ls='--', color='grey', lw=1)
    ax.axvline(x_threshold, ls='--', color='grey', lw=1)
    ax.axhline(y_threshold, ls='--', color='grey', lw=1)

    # 设置标签
    ax.set_xlabel('Correlation Coefficient (r)', fontweight='bold', fontsize=12)
    ax.set_ylabel('-Log10(adj. p-value)', fontweight='bold', fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # 标注选定的基因
    if selected_genes is not None:
        texts = []  # 用于存放所有的文本标签
        for gene in selected_genes:
            if gene in vol[id_col].values:  # 检查基因是否在数据中
                row = vol[vol[id_col] == gene]
                # 获取基因的坐标
                x = row[lfc_col].values[0]
                y = row['y'].values[0]

                # 使用合适的偏移调整文本位置，尽量向下
                text = ax.text(x, y - 0.1, gene,  # 后面数字可以根据需要自己调整
                               fontsize=10, style='italic', weight='bold', color='black')

                texts.append(text)  # 保存文本对象以便后续调整

        # 使用 adjust_text 来处理重叠和文本位置
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='black', lw=0.5),
                    only_move={'points': 'y', 'text': 'xy'})  # 只移动y轴上的文本

    fig.tight_layout()

    # 保存 PNG
    fig.savefig(output_file, dpi=300)
    # 额外保存 PDF
    fig.savefig(output_file.replace('.png', '.pdf'), dpi=300, bbox_inches='tight')

    plt.close()
    print("DataFrame head:\n", df.head())
    print("Correlation Coefficient (r) range:", df[lfc_col].min(), df[lfc_col].max())
    print("Adjusted p-value range:", df[p_col].min(), df[p_col].max())
    # 打印横坐标值的范围
    print("Correlation Coefficient (r) range:", vol[lfc_col].min(), vol[lfc_col].max())


# 示例使用
create_volcano_plot(
    input_file='D:\\Projects\\Thymus_Aging\\proteomics\\impute_spearman_correlation_res.csv',
    output_file='Volcano_plot.png',
    x_threshold=0,
    y_threshold=-np.log10(0.05),
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