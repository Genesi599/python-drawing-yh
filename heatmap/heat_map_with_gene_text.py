#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
横轴=样本  纵轴=显著基因 的热图（clustermap）
完全沿用你已有的颜色条、直方图、左侧色块代码
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.patches import Rectangle
from numba.cpython.new_builtins import max_iterable

# ========== 1. 只改这里 ==========
root_dir   = Path(r"D:\Projects\Neutrophil_Aging\leukocyte_single_cell\Monkey\monkey-B_cell_celltypist_filter\per_organ_OvsY")
PADJ_CUT   = 0.05
FC_CUT     = 0.5
TOP_ANNOT  = 20            # 右侧标注基因数
OUT_FIG    = Path("fig")   # 输出文件夹
OUT_FIG.mkdir(exist_ok=True)
# =================================

# 颜色映射

custom_seismic = sns.color_palette("coolwarm", as_cmap=True)
custom_seismic.set_under("#0c00f3")
custom_seismic.set_over("#ff0008")

linkage_methods = ['single', 'complete', 'ward']

def read_group_info(counts_df):
    """样本名 → old/young"""
    d = {"O": "old", "Y": "young"}
    return pd.Series([d.get(s[0], "unknown") for s in counts_df.index],
                     index=counts_df.index, name="condition")

def plot_sig_clustermap(expr, group_info, gene_fc, out_png, title=""):
    # 1. 取显著基因
    sig_genes = gene_fc[(gene_fc.abs() > FC_CUT) & gene_fc.notna()].index
    sig_genes = sig_genes.intersection(expr.columns)
    if sig_genes.empty:
        print("  无显著基因，跳过")
        return
    expr = expr[sig_genes]
    gene_fc = gene_fc.reindex(expr.columns)

    # 2. 使用原始表达量作为热图颜色；不进行 z-score 标准化
    expr = expr.loc[:, expr.var() > 0]  # 先删掉样本间无变化的基因
    expr_raw = expr.copy()              # 使用原始表达量作为热图颜色
    if expr_raw.empty or expr_raw.shape[1] < 3:
        print("  剩余基因<3，跳过")
        return

    # 3. 上下调分组
    up_genes = gene_fc[gene_fc > FC_CUT].index
    down_genes = gene_fc[gene_fc < -FC_CUT].index
    avail = []
    if len(up_genes) > 0:   avail.append(('up', up_genes))
    if len(down_genes) > 0: avail.append(('down', down_genes))
    if not avail:
        print("  无上调/下调基因，跳过")
        return

    # 4. 分别聚类后合并
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

    ordered_blocks = []
    for tag, genes in avail:
        block_raw = expr_raw[genes].T  # 使用原始表达量，基因为行，样本为列时的转置
        leaves = block_raw.index[leaves_list(linkage(pdist(block_raw, metric='euclidean'), method='ward'))]
        ordered_blocks.append((tag, block_raw.loc[leaves]))
    plot_mat = pd.concat([blk for _, blk in ordered_blocks])

    # 5. 样本年龄颜色映射
    lut = {"young": "#1f77b4", "old": "#d62728"}
    sample_colors = group_info.map(lut)

    # 0. 用数据决定颜色范围（使用原始表达量的最小/最大值）
    all_vals = expr_raw.values
    vmin = float(np.nanmin(all_vals))
    vmax = float(np.nanmax(all_vals))
    norm = plt.Normalize(vmin=vmin, vmax=vmax)

    # 6. 画 clustermap
    # 使用 plot_mat 作为输入，颜色依据原始表达量的范围
    g = sns.clustermap(
        plot_mat,
        method='ward',
        metric='euclidean',
        row_cluster=True,
        col_cluster=False,
        cmap=sns.color_palette("coolwarm", as_cmap=True),
        norm=norm,  # 使用原始表达量进行颜色映射
        xticklabels=True,
        yticklabels=False,
        figsize=(8, 9),
        cbar_pos=None,
        dendrogram_ratio=(0.15, 0.15)
    )

    # 7. 强制 Y 在前 O 在后（不用聚类结果）
    new_sample_order = [
                           s for s in plot_mat.columns if s.startswith('Y')
                       ] + [
                           s for s in plot_mat.columns if s.startswith('O')
                       ]
    plot_mat = plot_mat[new_sample_order]

    # 8. 重画 heatmap 使顺序生效
    g.ax_heatmap.clear()
    sns.heatmap(
        plot_mat,
        cmap=custom_seismic,
        norm=norm,
        ax=g.ax_heatmap,
        xticklabels=True,
        yticklabels=False,
        cbar=False
    )
    g.ax_heatmap.set_ylabel("")  # ← 加这一行
    g.ax_heatmap.set_xticklabels(plot_mat.columns, rotation=90, fontsize=13,
                                 ha='center', va='top', fontweight='bold')
    g.ax_heatmap.tick_params(axis='x', length=0)

    # 9. 画单条样本年龄颜色带（无白条）
    bar_height, y_base = 0.02, 1.01
    for i, s in enumerate(new_sample_order):
        g.ax_heatmap.add_patch(
            Rectangle((i, y_base), 1, bar_height,
                      facecolor=sample_colors[s],
                      transform=g.ax_heatmap.get_xaxis_transform(),
                      clip_on=False)
        )

    # 先拿到整行顺序（按 plot_mat 行序）
    all_genes = plot_mat.index.tolist()

    # TOP_ANNOT 现在用于右侧注释数量
    if len(gene_fc[gene_fc >  FC_CUT].index) > 0 or len(gene_fc[gene_fc < -FC_CUT].index) > 0:
        top_up   = gene_fc[gene_fc >  FC_CUT].nlargest(TOP_ANNOT).index.tolist()
        top_down = gene_fc[gene_fc < -FC_CUT].nsmallest(TOP_ANNOT).index.tolist()
        show_genes = top_up + top_down
        y_coords   = [all_genes.index(g) + 0.5 for g in show_genes]

        # 10. 右侧注释文本（保持原有实现思路）
        import adjustText as at
        fig = g.fig
        ax  = g.ax_heatmap

        # 10.1 物理白板与文本位置准备
        show_len = len(show_genes)
        if show_len > 0:
            # 位置与文本创建
            x_text = plot_mat.shape[1] + max(0.1, plot_mat.shape[1] * 0.02)
            texts = []
            np.random.seed(42)  # 可重复
            for gene, y in zip(show_genes, y_coords):
                y_jit = y + np.random.uniform(-0.15, 0.15)
                t = ax.text(x_text, y_jit, gene,
                            fontsize=11, fontweight='bold',
                            va='center', ha='left',
                            clip_on=False)
                texts.append(t)

            # 调整文本避免重叠
            move_step = 1
            tol = 0.3
            max_round = 100
            for rnd in range(max_round):
                at.adjust_text(texts,
                               ax=ax,
                               max_iter=1000,
                               max_move=(0, 5),
                               expand_text=(0, 5),
                               force_text=(0, 500),
                               force_points=(0, 0),
                               expand_points=(0, 0),
                               autoalign=False,
                               only_move={'points': 'y', 'text': 'y'},
                               lim=3000,
                               arrowprops=None,
                               save_steps=False,
                               add_objects=[])

                gene2y0 = {g: y for g, y in zip(show_genes, y_coords)}
                ys_cur = np.array([t.get_position()[1] for t in texts])
                used = np.zeros(len(ys_cur), bool)
                moved = False

                for i in range(len(ys_cur)):
                    if used[i]:
                        continue
                    group = np.abs(ys_cur - ys_cur[i]) <= tol
                    idx = np.where(group)[0]
                    if len(idx) == 1:
                        used[idx] = True
                        continue
                    used[idx] = True
                    moved = True

                    orig_ys = np.array([gene2y0[t.get_text()] for t in texts])[idx]
                    top_i = idx[orig_ys.argmin()]
                    bot_i = idx[orig_ys.argmax()]

                    x_top, y_top = texts[top_i].get_position()
                    x_bot, y_bot = texts[bot_i].get_position()
                    texts[top_i].set_position((x_top, y_top + move_step))
                    texts[bot_i].set_position((x_bot, y_bot - move_step))

                if not moved:
                    break

            for t in texts:
                t.set_position((x_text, t.get_position()[1]))
                t.set_ha('left')

            def get_left_edge(text_obj):
                bbox = text_obj.get_window_extent(renderer=fig.canvas.get_renderer())
                data_bbox = bbox.transformed(ax.transData.inverted())
                return data_bbox.x0, text_obj.get_position()[1]

            gene2y0 = {g: y for g, y in zip(show_genes, y_coords)}

            for t in texts:
                x_end, y_end = get_left_edge(t)
                x_start = plot_mat.shape[1]
                y_start = gene2y0[t.get_text()]
                ax.annotate('', xy=(x_end, y_end), xytext=(x_start, y_start),
                            arrowprops=dict(arrowstyle='-', lw=1, color='black', shrinkA=0),
                            annotation_clip=False)

        # 11. 在颜色条上方标 Young / Old
        y_text = y_base + bar_height + 0.01
        if any(s.startswith('O') for s in new_sample_order):
            first_o_idx = next(i for i, s in enumerate(new_sample_order) if s.startswith('O'))
            g.ax_heatmap.text((first_o_idx - 1) / 2, y_text, 'Young',
                              transform=g.ax_heatmap.get_xaxis_transform(),
                              ha='center', va='bottom', fontsize=25, fontweight='bold')
            g.ax_heatmap.text((first_o_idx + len(new_sample_order) - 1) / 2, y_text, 'Old',
                              transform=g.ax_heatmap.get_xaxis_transform(),
                              ha='center', va='bottom', fontsize=25, fontweight='bold')
        else:
            g.ax_heatmap.text(len(new_sample_order) / 2, y_text, 'Young',
                              transform=g.ax_heatmap.get_xheatmap_transform(),
                              ha='center', va='bottom', fontsize=25, fontweight='bold')

    # 11. 微调布局
    left, bottom, weight, height = 0.2, 0.2, 0.6, 0.7
    tree_size, tree_gap = 0.15, 0.01
    g.ax_heatmap.set_position([left, bottom, weight, height])
    g.ax_row_dendrogram.set_position([left-tree_size-tree_gap, 0.2, tree_size, 0.7])
    g.ax_col_dendrogram.set_position([left, bottom+height+tree_gap*8/10, weight, tree_size*8/10])
    for ax in [g.ax_row_dendrogram, g.ax_col_dendrogram]:
        for line in ax.collections:
            line.set_linewidth(2)

    plt.savefig(out_png, dpi=600, bbox_inches='tight')
    pdf_path = out_png.with_suffix('.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close()
    print(f"  ✅ {out_png.name} 完成")
    print(f"  ✅ {pdf_path.name} 完成")
    plt.close()
    print(f"  ✅ {out_png.name} 完成")

def main():
    for subdir in list(root_dir.iterdir()):
        if not subdir.is_dir():
            continue
        deg_file  = subdir / "DEG_OvsY.csv"
        expr_file = subdir / "merged_expression_matrix.tsv"
        if not deg_file.exists() or not expr_file.exists():
            continue

        print(f"\n>> {subdir.name}")
        deg = pd.read_csv(deg_file, index_col=0).dropna(subset=["padj"])
        sig = deg[(deg.padj < PADJ_CUT) & (abs(deg.log2FoldChange) > FC_CUT)]
        if sig.shape[0] == 0:
            print("  无显著基因，跳过")
            continue

        counts = pd.read_csv(expr_file, sep='\t', index_col=0).T.dropna(axis=1)
        group_s = read_group_info(counts)
        expr_sig = counts[sig.index.intersection(counts.columns)]
        if expr_sig.shape[1] == 0:
            print("  显著基因与表达矩阵无交集，跳过")
            continue

        # ✅ 直接保存到原目录
        out_png = subdir / f"{subdir.name}_sig_heatmap.png"
        plot_sig_clustermap(expr_sig, group_s, sig.log2FoldChange, out_png,
                            title=f"{subdir.name}  Significant DEG heatmap")


    print("\n=== All done ===")

if __name__ == "__main__":
    main()