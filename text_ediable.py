#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : text_ediable.py
@Date    : 2026/1/23 11:11
@Author  : yh109
"""
import matplotlib.pyplot as plt
import matplotlib as mpl


# 1. 让 matplotlib 别把字体子集化
mpl.rcParams['pdf.fonttype'] = 42          # 已经写了，保留
mpl.rcParams['ps.fonttype'] = 42
# 关键：关闭子集化（>=3.6 有效）
mpl.rcParams['pdf.use14corefonts'] = False   # 不用 14 种核心字体
# 2. 指定一个系统里肯定有的字体，避免 mpl 走自带字体
mpl.rcParams['font.family'] = 'Arial'      # 或 'DejaVu Sans', 'SimHei' 等


plt.savefig('bubble_plot_adjusted.pdf',
            dpi=300,             # 矢量图 dpi 其实无所谓
            bbox_inches='tight',
            pad_inches=0.1,
            metadata={'Creator': None, 'Producer': None})   # 可选：去掉多余元数据
