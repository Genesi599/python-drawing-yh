import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from adjustText import adjust_text
import matplotlib as mpl
from matplotlib.lines import Line2D

# ============================================================
# 全局设置
# ============================================================
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.use14corefonts'] = False
mpl.rcParams['font.family'] = 'Arial'

corr_type = 'Spearman'
font_scale = 2.0
input_dir = r"D:\Projects\B_Cell_Aging\B_cell_ratio\correlation_results"
output_dir = r"D:\Projects\B_Cell_Aging\B_cell_ratio\figure"
color_map_path = r"D:\Projects\B_Cell_Aging\B_cell_ratio\ref\tissue_system_mapping.csv"

corr_col = f'{corr_type}_Correlation'
pval_col = f'{corr_type}_P_value'

# ============================================================
# 读取颜色映射
# ============================================================
color_df = pd.read_csv(color_map_path)
tissue_colors = dict(zip(color_df['Tissue'], color_df['Color']))

# ============================================================
# 可用的 marker 形状列表（按物种分配）
# ============================================================
marker_list = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', 'h', '<', '>', 'p']

# ============================================================
# 读取所有物种数据
# ============================================================
csv_files = sorted(Path(input_dir).glob('*_tissue_age_correlation_results.csv'))
all_data = []
species_list = []

for csv_file in csv_files:
    df = pd.read_csv(csv_file)
    species = csv_file.stem.replace('_tissue_age_correlation_results', '')
    df['Species'] = species.capitalize()
    df['neg_log_p'] = -np.log10(df[pval_col].clip(lower=1e-300))
    df = df[df['neg_log_p'] <= 30]
    df['color'] = df['Tissue'].map(tissue_colors).fillna('#808080')
    all_data.append(df)
    species_list.append(species.capitalize())

df_all = pd.concat(all_data, ignore_index=True)

# 为每个物种分配 marker
species_markers = {sp: marker_list[i % len(marker_list)] for i, sp in enumerate(species_list)}

# ============================================================
# 绘图
# ============================================================
fig, ax = plt.subplots(figsize=(14, 10))

for species, grp in df_all.groupby('Species'):
    marker = species_markers[species]
    for _, row in grp.iterrows():
        ax.scatter(
            row[corr_col], row['neg_log_p'],
            s=300, alpha=0.6,
            color=row['color'],
            marker=marker,
            edgecolors='white', linewidths=0.5
        )

# ============================================================
# 标注所有显著的点（P < 0.05，即 -log10(P) > 1.301）
# ============================================================
sig_df = df_all[df_all['neg_log_p'] > -np.log10(0.05)]

texts = []
points = []
colors = []
for _, row in sig_df.iterrows():
    text = ax.text(
        row[corr_col], row['neg_log_p'], row['Tissue'],
        fontsize=9 * font_scale, color=row['color']
    )
    texts.append(text)
    points.append((row[corr_col], row['neg_log_p']))
    colors.append(row['color'])

adjust_text(texts, ax=ax)

for text, point, color in zip(texts, points, colors):
    text_pos = text.get_position()
    dist = np.sqrt((text_pos[0] - point[0]) ** 2 + (text_pos[1] - point[1]) ** 2)
    if dist > 0.5:
        ax.annotate(
            '', xy=point, xytext=text_pos,
            arrowprops=dict(arrowstyle='->', color=color, lw=1.2, linestyle=(0, (2, 5)))
        )

# ============================================================
# 参考线
# ============================================================
ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', linewidth=2, alpha=0.5, label='P = 0.05')
ax.axvline(x=0, color='gray', linestyle='--', linewidth=2, alpha=0.5)

# ============================================================
# 坐标轴标签
# ============================================================
ax.set_xlabel(f'{corr_type} Correlation', fontsize=12 * font_scale)
ax.set_ylabel('-log10(P value)', fontsize=12 * font_scale)
ax.tick_params(axis='both', labelsize=10 * font_scale)
ax.set_title('Tissue–Age Correlation Across Species', fontsize=14 * font_scale)

# ============================================================
# 图例：只保留物种 → marker 形状
# ============================================================
species_handles = [
    Line2D([0], [0], marker=species_markers[sp], color='gray', linestyle='None',
           markersize=10, label=sp)
    for sp in species_list
]

ax.legend(
    handles=species_handles, title='Species',
    loc='upper left', bbox_to_anchor=(1.02, 1.0),
    fontsize=8 * font_scale, title_fontsize=9 * font_scale,
    frameon=True
)

# ============================================================
# 保存
# ============================================================
os.makedirs(output_dir, exist_ok=True)
output_pdf = os.path.join(output_dir, 'combined_species_correlation_plot.pdf')
output_png = os.path.join(output_dir, 'combined_species_correlation_plot.png')

plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.savefig(output_png, dpi=300, bbox_inches='tight')
plt.close()

print(f"Generated: {output_pdf}")
print(f"Generated: {output_png}")