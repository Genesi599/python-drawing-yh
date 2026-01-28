#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
热图 + 物种图标叠加
横=组织（按本图基因集合内出现物种数降序）
纵=基因（上下调各按出现物种数降序取前20，上调在上）
"""

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import seaborn as sns
import matplotlib.image as mpimg
import re


def normalize_species_name(species_name):
    """
    将各种物种名称变体统一为标准化名称
    例如：Human_B_cell-cellxgene, Human_B_cell, human → human
    """
    name_lower = species_name.lower()

    # 定义映射规则（按优先级排序，更具体的模式在前）
    patterns = [
        (r'human.*', 'human'),  # 任何包含human的
        (r'mouse.*', 'mouse'),  # 任何包含mouse的
        (r'monkey.*', 'monkey'),  # 任何包含monkey的
        (r'rat.*', 'rat'),  # 任何包含rat的
        (r'macaque.*', 'monkey'),  # 猕猴也归为monkey
        (r'h.*sapiens.*', 'human'),  # 学名
        (r'm.*musculus.*', 'mouse'),  # 学名
        (r'r.*norvegicus.*', 'rat'),  # 学名
    ]

    for pattern, standard_name in patterns:
        if re.match(pattern, name_lower):
            return standard_name

    # 如果没有匹配，返回原始名称的小写版本
    return name_lower


def merge_species_data(*csv_paths, out_csv=None):
    """
    合并多个物种的 DEG 数据

    每个 CSV 应有列：tissue, names, logfoldchanges
    会自动添加 species 列（从文件名提取或参数指定）
    """
    dfs = []
    for csv_path in csv_paths:
        df = pd.read_csv(csv_path)
        if 'species' not in df.columns:
            # 从文件名提取物种名（如 human_deg.csv → human）
            species_name = Path(csv_path).stem.split('_')[0]
            df['species'] = species_name

        # ===== 新增：规范化物种名称 =====
        df['species'] = df['species'].apply(normalize_species_name)
        # 保留原始名称用于调试（可选）
        # df['species_original'] = df['species']

        dfs.append(df)

    merged = pd.concat(dfs, ignore_index=True)

    if out_csv:
        merged.to_csv(out_csv, index=False)
        print(f"✅ 合并数据已保存到 {out_csv}")

    return merged


def load_species_icons(svg_paths):
    species_icon_map = {}
    for path in svg_paths:
        png_path = Path(path).with_suffix('.png')  # 改成 .png
        filename = png_path.stem.lower()
        species_name = filename.split('_')[0].lower()
        species_icon_map[species_name] = str(png_path)
    return species_icon_map


def plot_heatmap_with_icons(
        in_csv,
        proteome_csv_path=None,
        out_dir=None,
        top_per_group=20,
        font_scale=1.2,
        figsize=(16, 9),
        dpi_png=600,
        species_icon_paths=None,
        gene_type='all',
):
    """
    绘制热图+物种图标叠加图
    """
    in_csv = Path(in_csv)
    if out_dir is None:
        out_dir = in_csv.parent / "figure"
    else:
        out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    # 1. 读表
    df = pd.read_csv(in_csv)

    # ===== 新增：确保species列已规范化（如果输入的是未处理的CSV）=====
    if 'species' in df.columns:
        df['species'] = df['species'].apply(normalize_species_name)

    all_species = sorted(df['species'].unique())
    print(f"检测到物种: {all_species}")  # 调试信息

    # 2. 物种图标加载
    if species_icon_paths is None:
        print("请提供物种图标 SVG 文件路径。")
        return
    species_icon_map = load_species_icons(species_icon_paths)

    # ===== 修改：直接使用标准化后的物种名匹配图标 =====
    species_to_icon = {}
    for sp in all_species:
        icon_path = species_icon_map.get(sp)
        if not icon_path:
            # 尝试其他匹配方式
            for icon_sp, path in species_icon_map.items():
                if icon_sp in sp or sp in icon_sp:
                    icon_path = path
                    break
        species_to_icon[sp] = icon_path

    # 调试：打印物种到图标的映射
    print("物种到图标的映射：")
    print(f"可用图标: {species_icon_map}")
    print(f"数据中的物种: {all_species}")
    print(f"最终映射: {species_to_icon}")

    # 3. 挑基因（按物种数排序，合并后的species已经去重）
    gene_cnt = df.groupby('names').size().sort_values(ascending=False)
    up_genes = gene_cnt.loc[df[df['logfoldchanges'] > 0]['names'].unique()].nlargest(top_per_group).index
    down_genes = gene_cnt.loc[df[df['logfoldchanges'] < 0]['names'].unique()].nlargest(top_per_group).index

    if gene_type == 'up':
        genes = list(up_genes)
    elif gene_type == 'down':
        genes = list(down_genes)
    else:
        genes = list(up_genes) + list(down_genes)

    # 4. 组织排序（按本图基因集合内出现物种数降序）
    sub_df = df[df['names'].isin(genes)]
    # ===== 注意：这里计算的是合并后的物种数（去重后）=====
    tissue_cnt = sub_df.groupby('tissue')['species'].nunique().sort_values(ascending=False)
    org_order = tissue_cnt.index
    print(f"组织排序（按物种数）: {list(org_order)}")

    # 5. 矩阵（多个物种同一格子需要特殊处理，这里取平均或最大绝对值）
    # ===== 修改：同一tissue-gene组合多个物种时，取平均logFC =====
    mat = (sub_df.groupby(['tissue', 'names'])['logfoldchanges'].mean()
           .unstack(fill_value=0).reindex(index=org_order, columns=genes))

    exp_mat = np.sign(mat) * np.log1p(np.abs(mat))
    from matplotlib import colors
    norm = colors.TwoSlopeNorm(vmin=-np.abs(exp_mat.values).max(),
                               vcenter=0,
                               vmax=np.abs(exp_mat.values).max())
    cmap_fill = plt.cm.coolwarm

    # 6. 画布 - 添加气泡图子图
    if proteome_csv_path:
        fig, (ax, ax_bubble) = plt.subplots(1, 2, figsize=(figsize[0] + 3, figsize[1]),
                                            gridspec_kw={'width_ratios': [5, 0.4]})
    else:
        fig, ax = plt.subplots(figsize=figsize)
        ax_bubble = None

    # 7. 热图
    sns.heatmap(
        exp_mat.T,
        cmap=cmap_fill,
        norm=norm,
        linewidths=0.3,
        linecolor='#f2f2f2',
        xticklabels=True,
        yticklabels=True,
        cbar=False,
        ax=ax
    )

    # 8. overlay图标（显示所有物种）
    # ===== 修改：显示该tissue-gene组合中所有的物种（去重后）=====
    for i, tissue in enumerate(exp_mat.T.columns):  # 横轴
        for j, gene in enumerate(exp_mat.T.index):  # 纵轴
            cx = i + 0.5
            cy = j + 0.5
            # 获取该组合的所有唯一物种（已合并）
            sub = df[(df['tissue'] == tissue) & (df['names'] == gene)]
            if sub.empty:
                continue

            # 获取唯一的物种列表
            unique_species = sub['species'].unique()

            from matplotlib.offsetbox import OffsetImage, AnnotationBbox

            # 动态调整图标间距
            n_icons = len(unique_species)
            if n_icons > 0:
                offset_range = 0.2
                if n_icons == 1:
                    offsets = [0]
                else:
                    offsets = np.linspace(-offset_range, offset_range, n_icons)

                for idx, sp in enumerate(unique_species):
                    offset_x = offsets[idx]
                    icon_path = species_to_icon.get(sp)
                    if not icon_path or not Path(icon_path).exists():
                        continue
                    try:
                        image = mpimg.imread(icon_path)
                        imagebox = OffsetImage(image, zoom=0.07)
                        ab = AnnotationBbox(imagebox, (cx + offset_x, cy),
                                            frameon=False, xycoords='data',
                                            box_alignment=(0.5, 0.5))
                        ax.add_artist(ab)
                    except Exception as e:
                        print(f"加载图标失败 {icon_path}: {e}")

    # 9. 坐标轴
    ax.set_xlabel('Tissue', fontsize=14 * font_scale)
    ax.set_ylabel('Gene', fontsize=14 * font_scale)
    ax.set_title('Cross-Species Aging B Cell DEG', fontsize=16 * font_scale)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=14 * font_scale)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=14 * font_scale)

    # 10. 颜色条
    if ax_bubble:
        cbar_ax = fig.add_axes([0.81, 0.15, 0.02, 0.3])
    else:
        cbar_ax = fig.add_axes([0.81, 0.15, 0.03, 0.3])
    sm = plt.cm.ScalarMappable(cmap=cmap_fill, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
    cbar.ax.tick_params(labelsize=12 * font_scale)
    cbar.set_label('exp-logFC', fontsize=12 * font_scale)

    plt.tight_layout()

    if proteome_csv_path and ax_bubble:
        prot_df = pd.read_csv(proteome_csv_path)
        prot_df = prot_df.dropna(subset=['gene_name', 'r', 'pval_adj'])

        gene_order = exp_mat.T.index.tolist()
        prot_subset = prot_df[prot_df['gene_name'].isin(gene_order)].copy()
        prot_subset = prot_subset.sort_values('pval_adj').drop_duplicates('gene_name')
        prot_subset['gene_order'] = prot_subset['gene_name'].map(
            {g: i for i, g in enumerate(gene_order)}
        )
        prot_subset = prot_subset.dropna(subset=['gene_order']).sort_values('gene_order')

        print(f"\n气泡绘制前:")
        print(f"prot_subset行数: {len(prot_subset)}")
        print(f"r值范围: {prot_subset['r'].min():.3f} ~ {prot_subset['r'].max():.3f}")
        print(f"pval_adj范围: {prot_subset['pval_adj'].min():.3e} ~ {prot_subset['pval_adj'].max():.3e}")
        print(f"ax_bubble 是否为None: {ax_bubble is None}")

        # 绘制气泡图
        y_pos = np.arange(len(prot_subset))
        sizes = np.abs(prot_subset['r'].values) * 5000
        colors = -np.log10(prot_subset['pval_adj'].values)

        print(f"气泡大小范围: {sizes.min():.1f} ~ {sizes.max():.1f}")
        print(f"气泡颜色值范围: {colors.min():.2f} ~ {colors.max():.2f}")

        scatter = ax_bubble.scatter(np.ones_like(y_pos), y_pos, s=sizes,
                                    c=colors, cmap='coolwarm', alpha=0.7,
                                    edgecolors='black', linewidth=0.5)
        print("气泡已绘制")

    # 10b. 绘制气泡图（如果提供了蛋白质数据）
    if proteome_csv_path and ax_bubble:
        prot_df = pd.read_csv(proteome_csv_path)
        prot_df = prot_df.dropna(subset=['gene_name', 'r', 'pval_adj'])

        # 用热图的基因顺序和名称
        gene_order = exp_mat.T.index.tolist()
        prot_subset = prot_df[prot_df['gene_name'].isin(gene_order)].copy()

        # 为每个基因保留一行数据（如果多行则取p值最小的）
        prot_subset = prot_subset.sort_values('pval_adj').drop_duplicates('gene_name')
        prot_subset['gene_order'] = prot_subset['gene_name'].map(
            {g: i for i, g in enumerate(gene_order)}
        )
        prot_subset = prot_subset.dropna(subset=['gene_order']).sort_values('gene_order')

        # 绘制气泡图
        y_pos = np.arange(len(prot_subset))
        sizes = np.abs(prot_subset['r'].values) * 200  # 气泡大小
        colors = -np.log10(prot_subset['pval_adj'].values)
        colors = np.nan_to_num(colors, posinf=np.nanmax(colors[~np.isinf(colors)]))  # 替换inf

        scatter = ax_bubble.scatter(np.ones_like(y_pos), y_pos, s=sizes,
                                    c=colors, cmap='coolwarm', alpha=0.7,
                                    edgecolors='black', linewidth=0.5)

        ax_bubble.set_yticks(y_pos)
        ax_bubble.set_yticklabels(prot_subset['gene_name'].values, fontsize=9 * font_scale)
        ax_bubble.set_xticks([])
        ax_bubble.set_xlim(-0.5, 1.5)
        ax_bubble.set_ylim(-0.5, len(prot_subset) - 0.5)
        ax_bubble.invert_yaxis()
        ax_bubble.set_xlabel('Proteome', fontsize=12 * font_scale)

        # 颜色条
        cbar2 = plt.colorbar(scatter, ax=ax_bubble)
        cbar2.set_label('-log10(pval_adj)', fontsize=10 * font_scale)

    # 11. 保存
    for ext in ['png', 'pdf']:
        out_path = out_dir / f'heatmap_icons.{ext}'
        fig.savefig(out_path, dpi=dpi_png if ext == 'png' else None, bbox_inches='tight')
        print(f'  ✅ {out_path.name} 完成')
    plt.close()


def convert_svg_to_png(svg_paths):
    """将SVG转换为PNG"""
    import subprocess
    from pathlib import Path
    import shutil

    # 查找 Inkscape
    inkscape_cmd = shutil.which('inkscape')
    if not inkscape_cmd:
        # Windows 默认路径
        possible_paths = [
            r'C:\Program Files\Inkscape\bin\inkscape.exe',
            r'C:\Program Files (x86)\Inkscape\bin\inkscape.exe',
        ]
        for path in possible_paths:
            if Path(path).exists():
                inkscape_cmd = path
                break

    if not inkscape_cmd:
        raise FileNotFoundError('找不到 Inkscape，请确保已安装并添加到 PATH')

    for svg_path in svg_paths:
        png_path = Path(svg_path).with_suffix('.png')
        if png_path.exists():
            continue
        subprocess.run([
            inkscape_cmd,
            str(svg_path),
            f'--export-filename={str(png_path)}',
            '--export-type=png',
            '--export-dpi=96'
        ], check=True)
        print(f'✅ {png_path.name} 转换完成')


# 使用示例
def main():
    in_csv = Path(r"D:\Projects\B_Cell_Aging\B_cell_tissue_DEG_results_filtered.csv")
    svg_paths = [
        r"C:\Users\yh109\Documents\GitHub\python-drawing-yh\icon\lib\human.svg",
        r"C:\Users\yh109\Documents\GitHub\python-drawing-yh\icon\lib\mouse.svg",
        r"C:\Users\yh109\Documents\GitHub\python-drawing-yh\icon\lib\monkey.svg",
        r"C:\Users\yh109\Documents\GitHub\python-drawing-yh\icon\lib\rat.svg"
    ]
    # convert_svg_to_png(svg_paths)
    plot_heatmap_with_icons(
        in_csv=in_csv,
        proteome_csv_path=r"D:\Projects\Thymus_Aging\overlap\thymus_data\proteome_up.csv",
        top_per_group=20,
        font_scale=1.5,
        gene_type='up',
        species_icon_paths=svg_paths,
        figsize=(16, 12)
    )


if __name__ == '__main__':
    main()