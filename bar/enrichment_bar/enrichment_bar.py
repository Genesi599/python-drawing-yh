import os
import gseapy as gp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# 假设您提供的基因列表
gene_list = [
    "AHCY", "AHSA1", "ALDH16A1", "ANXA5", "ARF3", "BCL2L1",
    "BID", "BLVRA", "CD2AP", "CLDN5", "DBN1", "DPP3",
    "DTYMK", "EEF1D", "ENO2", "FHL2", "FHOD1", "FKBP5",
    "FN3KRP", "GRAP2", "GYG1", "HADH", "KHSRP", "MOB1B",
    "MPP1", "NANS", "NCK2", "NUDT5", "PAICS", "PARVB",
    "PCBP1", "PCYT2", "PDIA4", "PGD", "PGM2L1", "PGP",
    "PPI1", "PTK2B", "PTP4A2", "PYGB", "RBM38", "RGS10",
    "SDC4", "SERPINB9", "SGTA", "SHMT1", "SNTB1", "STIM1",
    "STMN1", "SUS1", "TBC1D10B", "TKT", "UCHL3", "VAT1",
    "VPS26B", "YKT6"
]

# 确保目录存在
output_directory = 'data'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

figure_directory = 'figure'
if not os.path.exists(figure_directory):
    os.makedirs(figure_directory)

# GO 富集分析
go_enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets='GO_Biological_Process_2025',
    organism='Human',
    outdir=None
)

go_df = go_enr.results  # pandas.DataFrame

# 保存 GO 结果
go_df.to_csv(os.path.join(output_directory, 'GO_result.csv'), index=True, encoding='utf-8')

# 计算 -log10(p-value) 并排序 GO
go_df['neg_log10_p'] = -np.log10(go_df['P-value'])
go_df = go_df.sort_values('neg_log10_p', ascending=True).tail(10)

# 绘制 GO 图
plt.figure(figsize=(8, 5))
y_pos_go = np.arange(len(go_df))
plt.barh(y_pos_go, go_df['neg_log10_p'], color='#e07c44', height=0.6)

# 添加基因名标注
for idx, (term, genes, x) in enumerate(zip(go_df['Term'], go_df['Genes'], go_df['neg_log10_p'])):
    plt.text(0.1, idx, term, va='center', ha='left', color='black', fontsize=10, fontweight='bold')
    gene_txt = ', '.join(genes.split(';')[:6])  # 只取前6个基因
    plt.text(0.1, idx - 0.3, f'({gene_txt})', va='top', ha='left', fontsize=10, color='black', fontstyle='italic')

plt.xlabel('-Log10(p-value)', fontsize=16)
plt.title('GO Terms and Pathways', pad=20, loc='left', fontsize=18)
plt.tight_layout()
plt.savefig(os.path.join(figure_directory, 'GO_BP_bar_with_genes.png'), dpi=600)
plt.show()

# KEGG 富集分析
kegg_enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets='KEGG_2021_Human',  # 请确认使用有效的基因集
    organism='Human',
    outdir=None
)

kegg_df = kegg_enr.results  # pandas.DataFrame

# 保存 KEGG 结果
kegg_df.to_csv(os.path.join(output_directory, 'KEGG_result.csv'), index=True, encoding='utf-8')

# 计算 -log10(p-value) 并排序 KEGG
kegg_df['neg_log10_p'] = -np.log10(kegg_df['P-value'])
kegg_df = kegg_df.sort_values('neg_log10_p', ascending=True).tail(10)

# 绘制 KEGG 图
plt.figure(figsize=(8, 5))
y_pos_kegg = np.arange(len(kegg_df))
plt.barh(y_pos_kegg, kegg_df['neg_log10_p'], color='#6b5b95', height=0.6)

# 添加基因名标注
for idx, (term, genes, x) in enumerate(zip(kegg_df['Term'], kegg_df['Genes'], kegg_df['neg_log10_p'])):
    plt.text(0.1, idx, term, va='center', ha='left', color='black', fontsize=10, fontweight='bold')
    gene_txt = ', '.join(genes.split(';')[:6])  # 只取前6个基因
    plt.text(0.1, idx - 0.3, f'({gene_txt})', va='top', ha='left', fontsize=10, color='black', fontstyle='italic')

plt.xlabel('-Log10(p-value)', fontsize=16)
plt.title('KEGG Pathways', pad=20, loc='left', fontsize=18)
plt.tight_layout()
plt.savefig(os.path.join(figure_directory, 'KEGG_BP_bar_with_genes.png'), dpi=600)
plt.show()