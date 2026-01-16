import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pathlib

# ----------------  用户可调参数  ---------------- #
y_col = 'B_cell_ratio_2'          # <--- 只改这里即可换列
# --------------------------------------------- #

# ----------------  样式统一  ---------------- #
sns.set_theme(style="ticks", font_scale=1.3, rc={
    "font.family": "Arial",
    "axes.linewidth": 1.2,
    "axes.edgecolor": "black",
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.transparent": False
})

palette = {"Young": "#0072B2", "Old": "#D55E00"}

# ----------------  数据读取  ---------------- #
file = pathlib.Path(r"C:\Users\yh109\OneDrive\桌面\test.xlsx")
df = pd.read_excel(file)

df = (df
      .assign(Organ=lambda x: x['Tissue'].astype(str).str.extract(r'(^[^_]+CD45)')[0],
              GroupLabel=lambda x: x['group'].map({'Y': 'Young', 'O': 'Old'}))
      .loc[lambda x: x['GroupLabel'].isin(['Young', 'Old'])]
     )
# 确保 y_col 为数值
df[y_col] = pd.to_numeric(df[y_col], errors='coerce')
df = df.dropna(subset=[y_col])

# ----------------  排序  ---------------- #
order = (df.groupby('Organ', observed=False)[y_col].mean()
         .sort_values(ascending=False)
         .index.to_list())
df['Organ'] = pd.Categorical(df['Organ'], categories=order, ordered=True)

# ----------------  绘图  ---------------- #
fig, ax = plt.subplots(figsize=(len(order)*1.2, 5), constrained_layout=True)

# 箱线+散点
sns.boxplot(data=df, x='Organ', y=y_col, hue='GroupLabel',
            palette=palette, width=0.6, linewidth=1.5,
            showcaps=False, boxprops=dict(facecolor='none'), ax=ax)
sns.stripplot(data=df, x='Organ', y=y_col, hue='GroupLabel',
              palette=palette, size=4, dodge=True, jitter=False, alpha=0.8, ax=ax)

# 均值条带
for i, organ in enumerate(order):
    for j, age in enumerate(['Young', 'Old']):
        sub = df[(df['Organ'] == organ) & (df['GroupLabel'] == age)]
        mean_val = sub[y_col].mean()
        x_pos = i - 0.2 + j*0.4
        ax.plot([x_pos, x_pos+0.2], [mean_val, mean_val],
                color=palette[age], lw=4, solid_capstyle='round')

# 坐标轴
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', style='italic')
ax.set_xlabel(None)
ax.set_ylabel("B cell ratio (%)", fontsize=14)
sns.despine(ax=ax, trim=True)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], frameon=False, loc='upper right')
fig.suptitle("B cell ratio in Macaca fascicularis", fontsize=16)
# ----------------  保存  ---------------- #
out_path = pathlib.Path(__file__).with_suffix('').name
fig.savefig(out_path + '.pdf', facecolor='white')
fig.savefig(out_path + '.png', dpi=600, facecolor='white')

plt.show()