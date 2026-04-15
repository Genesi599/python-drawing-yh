#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
多项式拟合散点图 - 使用示例

原始文件: scatter.py (CLEC10 数据多项式拟合)
"""

import pandas as pd
import numpy as np
from poly_fit import plot_poly_fit

# ============================================================
# 示例：读取数据并绘制多项式拟合图
# ============================================================
sheet_name = 'Sheet2'
df = pd.read_excel('CLEC10.xlsx', sheet_name=sheet_name)

y_name = 'IFNb'
df = df.dropna()
df = df[df[y_name] != 0]

fig, ax, r2 = plot_poly_fit(
    x=df['age'],
    y=df[y_name],
    degree=4,
    xlabel='Age (year)',
    ylabel='IFNB1',
    output_path='poly_fit_example'
)

print(f'R-squared: {r2:.3f}')
