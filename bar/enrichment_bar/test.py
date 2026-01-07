#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
只需改下面两行
"""
CSV_PATH = r"C:\Users\yh109\OneDrive\桌面\bulk_go.csv"   # 原始结果
# ------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def main(csv_file=CSV_PATH):
    # 1. 读入
    df = pd.read_csv(csv_file, sep=None, engine='python', encoding='gbk')

    # 2. 确保 neg_log10_p
    if 'neg_log10_p' not in df.columns:
        df['neg_log10_p'] = -np.log10(df['Adjusted P-value'].clip(1e-300))

    # 3. 四个筛选条件（按最后两列）
    masks = {
        'NET up'            : (df['change'] == 'up')   & (df['class'] == 'NET'),
        'NET up (reversible by DNaseI)'  : (df['change'] == 'down')   & (df['class'] == 'NETaes rescue'),
        'NET down'          : (df['change'] == 'down') & (df['class'] == 'NET'),
        'NET down (reversible by DNaseI)': (df['change'] == 'up') & (df['class'] == 'NETaes rescue'),
    }

    # 4. 画图函数
    def _draw(mask, title, color):
        sub = (df[mask]
               .sort_values('neg_log10_p', ascending=False)
               .head(8)
               .iloc[::-1])          # ← 倒序，让最显著的在上
        fig, ax = plt.subplots(figsize=(7, 4))
        y_pos = np.arange(len(sub))
        ax.barh(y_pos, sub['neg_log10_p'], color=color, height=0.6, alpha=0.75)
        for idx, (term, genes, x) in enumerate(zip(sub['Term'], sub['Genes'], sub['neg_log10_p'])):
            ax.text(0.1, idx, term, va='center', ha='left', fontsize=11, fontweight='bold')
            gene_txt = ', '.join(genes.split(';')[:6])
            ax.text(0.1, idx-0.3, f'({gene_txt})', va='top', ha='left', fontsize=11,
                    color='red' if 'up' in title.lower() else 'blue', fontweight='bold', fontstyle='italic')
        ax.set_xlabel('-Log10(p-value)', fontsize=13, fontweight='bold')
        ax.set_title(title, pad=10, loc='left', fontsize=15, fontweight='bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(left=False, labelleft=False)
        return fig

    # 5. 逐个保存为 PDF
    order = [
        ('NET up',            masks['NET up'],            '#FFB38A'),
        ('NET up (reversible by DNaseI)',  masks['NET up (reversible by DNaseI)'],  '#FF8E72'),
        ('NET down',          masks['NET down'],          '#81C7F4'),
        ('NET down (reversible by DNaseI)',masks['NET down (reversible by DNaseI)'],'#69D4D4'),
    ]
    pdf_names = ['NET_up.pdf', 'NETaes_rescue_up.pdf', 'NET_down.pdf', 'NETaes_rescue_down.pdf']
    png_names = ['NET_up.png', 'NETaes_rescue_up.png', 'NET_down.png', 'NETaes_rescue_down.png']

    for (title, mask, color), pdf, png in zip(order, pdf_names, png_names):
        fig = _draw(mask, title, color)

        # 保存为 PDF
        fig.savefig(pdf, format='pdf', bbox_inches='tight')

        # 保存为 PNG
        fig.savefig(png, format='png', bbox_inches='tight')

        plt.close(fig)
        print(f'saved -> {os.path.abspath(pdf)}')
        print(f'saved -> {os.path.abspath(png)}')

if __name__ == '__main__':
    cmd_csv = sys.argv[1] if len(sys.argv) > 1 else CSV_PATH
    main(cmd_csv)