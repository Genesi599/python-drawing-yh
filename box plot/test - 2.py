import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib

# ----------------  用户可调参数  ---------------- #
y_col = 'B_cell_ratio_2'          # <--- 想画哪列改这里
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

single_color = "#0072B2"   # Young 单色

# ----------------  数据读取与过滤  ---------------- #
file = pathlib.Path(r"C:\Users\yh109\OneDrive\桌面\test.xlsx")
df = pd.read_excel(file)

df = (df
      .assign(Organ=lambda x: x['Tissue'].astype(str).str.extract(r'(^[^_]+CD45)')[0])
      .query("group == 'Y'")                          # 只留 Young
     )
df[y_col] = pd.to_numeric(df[y_col], errors='coerce')
df = df.dropna(subset=[y_col])

# ----------------  排序  ---------------- #
order = (df.groupby('Organ', observed=False)[y_col].mean()
         .sort_values(ascending=False)
         .index.to_list())
df['Organ'] = pd.Categorical(df['Organ'], categories=order, ordered=True)

# ----------------  绘图  ---------------- #
fig, ax = plt.subplots(figsize=(len(order)*1.2, 5), constrained_layout=True)

# 箱线+散点（单色）
sns.boxplot(data=df, x='Organ', y=y_col, color=single_color,
            width=0.6, linewidth=1.5, showcaps=False,
            boxprops=dict(facecolor='none'), ax=ax)
sns.stripplot(data=df, x='Organ', y=y_col, color=single_color,
              size=4, alpha=0.8, ax=ax)

# 均值横线
for i, organ in enumerate(order):
    mean_val = df[df['Organ'] == organ][y_col].mean()
    ax.plot([i-0.2, i+0.2], [mean_val, mean_val],
            color=single_color, lw=4, solid_capstyle='round')

# 坐标轴
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', style='italic')
ax.set_xlabel(None)
ax.set_ylabel("B cell ratio (%)", fontsize=14)
sns.despine(ax=ax, trim=True)
fig.suptitle("B cell ratio in Macaca fascicularis (young)", fontsize=16)
# ----------------  保存  ---------------- #
out_path = pathlib.Path(__file__).with_suffix('').name + '_Young'
fig.savefig(out_path + '.pdf', facecolor='white')
fig.savefig(out_path + '.png', dpi=600, facecolor='white')

plt.show()