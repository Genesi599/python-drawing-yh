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

from drawing_yh import DEFAULT_FONT_SIZE, save_fig


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
        out_dir=None,
        top_per_group=20,
        font_scale=1.0,
        figsize=(16, 9),
        dpi_png=300,
        species_icon_paths=None,
        gene_type='all',
        title='Cross-species B cell DEG',
        output_stem='heatmap_icons',
):
    """
    绘制热图+物种图标叠加图
    """
    in_csv = Path(in_csv)
    if out_dir is None:
        out_dir = in_csv.parent / "figure"
    else:
        out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

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

    # 6. 画布
    fig, ax = plt.subplots(figsize=figsize)

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
    base_fontsize = DEFAULT_FONT_SIZE * font_scale
    ax.set_xlabel('Tissue', fontsize=base_fontsize)
    ax.set_ylabel('Gene', fontsize=base_fontsize)
    if title:
        ax.set_title(title, fontsize=base_fontsize, pad=8)
    plt.xticks(rotation=45, ha='right', fontsize=base_fontsize)
    plt.yticks(rotation=0, fontsize=base_fontsize)

    # 10. 颜色条
    cbar_ax = fig.add_axes([0.81, 0.15, 0.03, 0.3])
    sm = plt.cm.ScalarMappable(cmap=cmap_fill, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
    cbar.ax.tick_params(labelsize=base_fontsize)
    cbar.set_label('exp-logFC', fontsize=base_fontsize)

    plt.tight_layout()
    # 11. 保存
    written = save_fig(fig, out_dir / f'{output_stem}.svg', also=('.pdf', '.png'))
    for out_path in written:
        print(f'  {out_path.name} saved')
    plt.close()
    return written


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
        top_per_group=20,
        font_scale=1.5,
        gene_type='up',
        species_icon_paths=svg_paths
    )


if __name__ == '__main__':
    main()
