import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from adjustText import adjust_text

import matplotlib as mpl


# 1. 让 matplotlib 别把字体子集化
mpl.rcParams['pdf.fonttype'] = 42          # 已经写了，保留
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'        # SVG 文字保持文字对象,不转曲
# 关键：关闭子集化（>=3.6 有效）
mpl.rcParams['pdf.use14corefonts'] = False   # 不用 14 种核心字体
# 2. 指定一个系统里肯定有的字体，避免 mpl 走自带字体
mpl.rcParams['font.family'] = 'Arial'      # 或 'DejaVu Sans', 'SimHei' 等

corr_type = 'Spearman'
font_scale = 3.0
input_dir = r"D:\Projects\B_Cell_Aging\B_cell_ratio\correlation_results"
output_dir = r"D:\Projects\B_Cell_Aging\B_cell_ratio\figure"
color_map_path = r"D:\Projects\B_Cell_Aging\B_cell_ratio\ref\tissue_system_mapping.csv"

color_df = pd.read_csv(color_map_path)
tissue_colors = dict(zip(color_df['Tissue'], color_df['Color']))

csv_files = list(Path(input_dir).glob('*_tissue_age_correlation_results.csv'))

for csv_file in csv_files:
    df = pd.read_csv(csv_file)

    species = csv_file.stem.replace('_tissue_age_correlation_results', '')

    corr_col = f'{corr_type}_Correlation'
    pval_col = f'{corr_type}_P_value'

    df['neg_log_p'] = -np.log10(df[pval_col].clip(lower=1e-300))
    df = df[df['neg_log_p'] <= 30]
    df['color'] = df['Tissue'].map(tissue_colors).fillna('#808080')

    top20_df = df.nlargest(12, 'neg_log_p')

    plt.figure(figsize=(8, 8))

    for idx, row in df.iterrows():
        plt.scatter(row[corr_col], row['neg_log_p'], s=400, alpha=0.6, color=row['color'])

    texts = []
    points = []
    colors = []
    for idx, row in top20_df.iterrows():
        offset_x = np.random.uniform(-0.01, 0.01)
        offset_y = np.random.uniform(-0.1, 0.1)
        text = plt.text(row[corr_col] + offset_x, row['neg_log_p'] + offset_y, row['Tissue'],
                        fontsize=9 * font_scale, color=row['color'])
        texts.append(text)
        points.append((row[corr_col], row['neg_log_p']))
        colors.append(row['color'])

    adjust_text(texts)

    for text, point, color in zip(texts, points, colors):
        text_pos = text.get_position()
        distance = np.sqrt((text_pos[0] - point[0]) ** 2 + (text_pos[1] - point[1]) ** 2)
        if distance > 0.5:
            plt.annotate('', xy=point, xytext=text_pos,
                         arrowprops=dict(arrowstyle='->', color=color, lw=1.5, linestyle=(0, (2, 5))))

    plt.xlabel(f'{corr_type} Correlation', fontsize=12*font_scale)
    plt.ylabel('-log10(P value)', fontsize=12*font_scale)

    plt.xlabel(f'{corr_type} Correlation', fontsize=12 * font_scale)
    plt.ylabel('-log10(P value)', fontsize=12 * font_scale)
    plt.tick_params(axis='both', labelsize=10 * font_scale)

    plt.title(f'{species.capitalize()}', fontsize=14 * font_scale)

    plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', linewidth=3, alpha=0.5)
    plt.axvline(x=0, color='gray', linestyle='--', linewidth=3, alpha=0.5)

    # plt.gca().set_aspect('equal', adjustable='box')

    # plt.tight_layout()

    output_pdf = os.path.join(output_dir, f'{species}_correlation_plot.pdf')
    output_png = os.path.join(output_dir, f'{species}_correlation_plot.png')
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Generated: {output_pdf}")
    print(f"Generated: {output_png}")