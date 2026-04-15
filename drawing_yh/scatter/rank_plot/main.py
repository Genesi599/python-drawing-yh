#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Rank abundance plot - 使用示例

原始文件: sort_fig.py (蛋白丰度排序图)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rank_plot import plot_rank

# ============================================================
# 示例：蛋白丰度 rank plot
# ============================================================
path = Path(r"C:\Users\yh109\Nutstore\1\Institute of Zoology\protein\FC.csv")
listpath = Path(r"C:\Users\yh109\Documents\GitHub\enrichment_analysis-yh\data\annolist.txt")

with open(listpath, 'r', encoding='utf-8') as f:
    annolist = [line.rstrip('\n') for line in f]

pd_df = pd.read_csv(path, index_col=0, encoding='utf-8')
pd_df = pd_df[:500]
pd_df['TS13'] = np.log10(pd_df['TS13'])
pd_df['protein'] = pd_df.index
pd_df = pd_df.sort_values(by='log2FC', ascending=True)
pd_df = pd_df.reset_index(drop=True)
pd_df['ID'] = pd_df.index

fig, ax = plot_rank(
    df=pd_df,
    x_col='ID',
    y_col='TS13',
    label_col='GENE',
    highlight_labels=annolist,
    xlabel='Rank',
    ylabel=r'log$_{10}$(Abundance)',
    output_path='rank_plot_example'
)
