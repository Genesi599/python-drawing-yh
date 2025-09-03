import gseapy as gp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

FC_df = pd.read_csv('FC.csv', index_col=0, encoding='utf-8')

# 设定阈值
thr = 1

# 新增列 'flag'
FC_df['flag'] = np.where(FC_df['log2FC'] >=  thr, 'up',
              np.where(FC_df['log2FC'] <= -thr, 'down', 'stable'))


df = FC_df.dropna(subset=['GENE'])
df_filtered = df[df['flag'] == 'down']

gene_list = df_filtered['GENE'].dropna().tolist()

# 假设你有一个基因列表，使用 EntrezGeneID、GeneSymbol 请确保与数据库一致
# gene_list = ["TP53", "BRCA1", "BRCA2", "EGFR", "VEGFA", "MYC"]

enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets='GO_Biological_Process_2025',   # 也可选 GO_Molecular_Function_2021 等
    organism='Human',                         # 支持 Human / Mouse / ...
    outdir=None                               # 不保存中间文件
)


df = enr.results      # pandas.DataFrame
print(df.columns)

df.to_csv('data/GO_result.csv', index=True, encoding='utf-8')


# ===== 2. 计算 -log10(p-value) 并排序 =====
df['neg_log10_p'] = -np.log10(df['P-value'])
df = df.sort_values('neg_log10_p', ascending=True).tail(10)   # 取显著的前20条

# ===== 3. 绘图 =====

fig, ax = plt.subplots(figsize=(8, 0.8*len(df)))
# 横向柱状图
y_pos = np.arange(len(df))
bar_height = 0.6
ax.barh(y_pos, df['neg_log10_p'], color='#e07c44', height=bar_height)

offset = 0.1  # 距离柱子左边缘的偏移量，可微调
# 把 term 写在柱子里面
for idx, (term, x) in enumerate(zip(df['Term'], df['neg_log10_p'])):
    # 让文字位于柱子中央；x/2 表示柱子中点

    ax.text(offset, idx, term,
            va='center', ha='left',
            color='black', fontsize=12, fontweight='bold')

# 基因标注：每个柱子右侧显示交集基因（用逗号隔开，太长可截断）
for idx, (genes, x) in enumerate(zip(df['Genes'], df['neg_log10_p'])):
    gene_txt = ', '.join(genes.split(';')[:6])   # 取前 6 个基因，防止过长
    gene_txt = '( ' + gene_txt + ' )'
    ax.text(
        offset,
        idx - 0.6*bar_height,
        gene_txt,
        va='top',
        ha='left',
        fontsize=12,
        color='black',
        fontstyle='italic'
    )

# 坐标轴美化
ax.set_yticks([])
ax.set_yticklabels([])

ax.set_xlabel('-Log10(p-value)', fontsize=16)
ax.tick_params(axis='x', labelsize=16)
ax.xaxis.set_major_locator(MultipleLocator(5))   # 2 可换成任意步长

ax.set_title('Go terms and pathways of genes coding protein', pad=20, loc='left', fontsize=18)

plt.tight_layout()
plt.savefig('figure/GO_BP_bar_with_genes.png', dpi=300)