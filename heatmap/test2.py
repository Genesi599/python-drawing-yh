from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

root_dir   = Path(r"D:\Projects\Neutrophil_Aging\leukocyte_single_cell\Monkey\monkey-B_cell_celltypist_filter\per_organ_OvsY")
PADJ_CUT   = 1
FC_CUT     = 0
TOP_ANNOT  = 20            # 右侧标注基因数
OUT_FIG    = Path("fig")   # 输出文件夹
OUT_FIG.mkdir(exist_ok=True)
# =================================

# 颜色映射

custom_seismic = sns.color_palette("coolwarm", as_cmap=True)
custom_seismic.set_under("#0c00f3")
custom_seismic.set_over("#ff0008")

linkage_methods = ['single', 'complete', 'ward']

from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# 全局配色
custom_seismic = sns.color_palette("coolwarm", as_cmap=True)
custom_seismic.set_under("#0c00f3")
custom_seismic.set_over("#ff0008")

FC_CUT   = 0
TOP_ANNOT = 20


def plot_sig_clustermap(expr: pd.DataFrame,
                        group_info: pd.Series,
                        out_png: Path,
                        gene_fc: pd.Series | None = None,
                        title: str = ""):
    """
    绘制“基因×样本”热图（无转置版）
    参数
    ----
    expr      : DataFrame,  index=基因, columns=样本
    group_info: Series,   index=样本, values=Young/Middle/Old
    gene_fc   : Series,   index=基因, values=logFC（可选）
    out_png   : Path,     输出png路径
    title     : str,      图标题（暂未绘制，可扩展）
    """
    # 0) 形状检查
    if not expr.columns.equals(group_info.index):
        raise ValueError("expr.columns 必须与 group_info.index 完全一致")

    # 1) 选基因
    TOP_N = 200
    if gene_fc is None or gene_fc.empty:
        var = expr.var(axis=1)                      # 按基因算方差
        sig_genes = var.nlargest(TOP_N).index
        gene_fc = pd.Series(0., index=sig_genes)
    else:
        sig_genes = gene_fc[(gene_fc.abs() > FC_CUT) & gene_fc.notna()].index

    sig_genes = expr.index.intersection(sig_genes)
    if sig_genes.empty:
        print("  无基因可画，跳过"); return

    expr_sub = expr.loc[sig_genes]
    # 删掉样本间无变化的基因
    expr_sub = expr_sub.loc[:, expr_sub.var() > 0]
    if expr_sub.empty or expr_sub.shape[1] < 3:
        print("  剩余基因<3，跳过"); return

    # 2) 上下调分组
    up_genes = gene_fc[gene_fc >  FC_CUT].index
    down_genes = gene_fc[gene_fc < -FC_CUT].index
    avail = []
    if len(up_genes):   avail.append(('up', up_genes))
    if len(down_genes): avail.append(('down', down_genes))
    if not avail:       avail = [('all', expr_sub.index)]

    # 3) 按组聚类
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

    ordered_blocks = []
    for tag, genes in avail:
        block = expr_sub.loc[genes]                # 行=基因，列=样本
        leaves = block.index[leaves_list(
            linkage(pdist(block, metric='euclidean'), method='ward'))]
        ordered_blocks.append((tag, block.loc[leaves]))
    plot_mat = pd.concat([blk for _, blk in ordered_blocks])

    # 4) 样本顺序 & 颜色
    lut = {"Young": "#1f77b4", "Middle": "#ff7f0e", "Old": "#d62728"}
    sample_colors = group_info.map(lut)
    age_order = {'Young': 0, 'Middle': 1, 'Old': 2}
    sample_order = sorted(
        plot_mat.columns,
        key=lambda s: age_order[group_info[s]]
    )
    plot_mat = plot_mat[sample_order]

    # 5) 颜色映射范围
    vmin = float(np.nanmin(plot_mat.values))
    vmax = float(np.nanmax(plot_mat.values))
    norm = plt.Normalize(vmin=vmin, vmax=vmax)

    # 6) 画 clustermap
    g = sns.clustermap(
        plot_mat,
        method='ward',
        metric='euclidean',
        row_cluster=True,
        col_cluster=False,
        cmap=custom_seismic,
        norm=norm,
        xticklabels=True,
        yticklabels=False,
        figsize=(8, 9),
        cbar_pos=None,
        dendrogram_ratio=(0.15, 0.15)
    )

    # 7) 样本颜色条
    bar_height, y_base = 0.02, 1.01
    for i, s in enumerate(sample_order):
        g.ax_heatmap.add_patch(
            Rectangle((i, y_base), 1, bar_height,
                      facecolor=sample_colors[s],
                      transform=g.ax_heatmap.get_xaxis_transform(),
                      clip_on=False)
        )

    # 8) 右侧基因注释
    if gene_fc is not None and (gene_fc.gt(FC_CUT).any() or gene_fc.lt(-FC_CUT).any()):
        top_up   = gene_fc[gene_fc >  FC_CUT].nlargest(TOP_ANNOT).index
        top_down = gene_fc[gene_fc < -FC_CUT].nsmallest(TOP_ANNOT).index
        show_genes = [g for g in top_up if g in plot_mat.index] + \
                     [g for g in top_down if g in plot_mat.index]
        if show_genes:
            import adjustText as at
            all_genes = plot_mat.index.tolist()
            y_coords  = [all_genes.index(g) + 0.5 for g in show_genes]
            x_text    = plot_mat.shape[1] + max(0.1, plot_mat.shape[1] * 0.02)

            texts = []
            np.random.seed(42)
            for gene, y in zip(show_genes, y_coords):
                y_jit = y + np.random.uniform(-0.15, 0.15)
                t = g.ax_heatmap.text(x_text, y_jit, gene,
                                      fontsize=11, fontweight='bold',
                                      va='center', ha='left', clip_on=False)
                texts.append(t)

            # adjust_text 参数保持你原设置
            at.adjust_text(texts, ax=g.ax_heatmap,
                           max_iter=1000, max_move=(0, 5),
                           expand_text=(0, 5), force_text=(0, 500),
                           only_move={'points':'y','text':'y'},
                           arrowprops=None, autoalign=False)

            # 画引线
            def left_edge(t):
                bbox = t.get_window_extent(renderer=g.fig.canvas.get_renderer()) \
                        .transformed(g.ax_heatmap.transData.inverted())
                return bbox.x0, t.get_position()[1]

            gene2y0 = dict(zip(show_genes, y_coords))
            for t in texts:
                x_end, y_end = left_edge(t)
                x_start = plot_mat.shape[1]
                y_start = gene2y0[t.get_text()]
                g.ax_heatmap.annotate('', xy=(x_end, y_end),
                                      xytext=(x_start, y_start),
                                      arrowprops=dict(arrowstyle='-', lw=1,
                                                      color='black', shrinkA=0),
                                      annotation_clip=False)

    # 9) 布局微调（保持你原比例）
    left, bottom, weight, height = 0.2, 0.2, 0.6, 0.7
    tree_size, tree_gap = 0.15, 0.01
    g.ax_heatmap.set_position([left, bottom, weight, height])
    g.ax_row_dendrogram.set_position([left-tree_size-tree_gap, bottom, tree_size, height])
    g.ax_col_dendrogram.set_position([left, bottom+height+tree_gap*8/10, weight, tree_size*8/10])
    for ax in [g.ax_row_dendrogram, g.ax_col_dendrogram]:
        for line in ax.collections:
            line.set_linewidth(2)

    # 10) 保存
    plt.savefig(out_png, dpi=600, bbox_inches='tight')
    plt.savefig(out_png.with_suffix('.pdf'), bbox_inches='tight')
    plt.close()
    print(f"  ✅ {out_png.name} 完成")


def main():
    # 1. 读表达矩阵 → 行样本，列蛋白
    expr_df = pd.read_csv(
        r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\abundance_sample_x_protein.csv",
        index_col=0
    )

    # 2. 读 meta
    meta = pd.read_csv(
        r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\sample_meta.csv",
        index_col="sample"
    )
    condition_map = meta["condition"].str.capitalize()

    # 3. 对齐样本
    common = expr_df.index.intersection(condition_map.index)
    expr_df = expr_df.loc[common]
    condition_map = condition_map.loc[common]

    # 4. 读 logFC
    df = pd.read_csv(
        r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\Pattern_Analysis\all_proteins_pattern.csv",
        index_col="gene"
    )
    logfc = pd.to_numeric(df["log2FC_Young_vs_Old"], errors="coerce")
    logfc = logfc.reindex(expr_df.columns)

    # 5. 转置 → 现在 index=基因，columns=样本
    expr_df = expr_df.T

    # 6. 调用
    out_png = Path("outputs/protein_condition_heatmap.png")
    out_png.parent.mkdir(exist_ok=True)
    plot_sig_clustermap(
        expr=expr_df,
        group_info=condition_map,
        out_png=out_png,
        gene_fc=logfc,
        title="Protein expression by condition"
    )
if __name__ == "__main__":
    main()