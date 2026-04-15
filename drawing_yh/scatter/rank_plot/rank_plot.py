#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Rank abundance plot - 按丰度排序的散点图，支持差异基因上下调着色与标注
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from adjustText import adjust_text
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.use14corefonts'] = False
mpl.rcParams['font.family'] = 'Arial'


def plot_rank(df, x_col, y_col, label_col,
              highlight_labels=None,
              up_color='tab:red', down_color='tab:blue',
              ns_color='dimgrey', bg_color='black',
              x_threshold=0, y_threshold=None,
              xlabel='Rank', ylabel=r'log$_{10}$(Abundance)',
              n_label=None, point_size=5, font_scale=1.0,
              figsize=None, output_path=None):
    """
    绘制 rank abundance 散点图。

    Parameters
    ----------
    df : DataFrame
        包含 x_col, y_col, label_col 的数据
    x_col, y_col : str
        横纵坐标列名
    label_col : str
        标注列名
    highlight_labels : list, optional
        需要标注的基因/蛋白名列表，若为 None 则不标注
    up_color, down_color, ns_color, bg_color : str
        上调、下调、非显著、背景点颜色
    x_threshold : float
        x 轴阈值（用于区分上下调）
    y_threshold : float, optional
        y 轴阈值（低于此值的点为灰色），默认取 y 最小整数
    xlabel, ylabel : str
        坐标轴标签
    n_label : int, optional
        若无 highlight_labels，标注前 n_label 个最显著点
    point_size : int
        散点大小
    font_scale : float
        字体缩放
    figsize : tuple, optional
        图片大小
    output_path : str, optional
        输出路径（不含扩展名）
    """
    df = df.copy()

    if y_threshold is None:
        y_threshold = math.floor(df[y_col].min())

    # 着色
    df['_color'] = bg_color
    df.loc[(df[x_col] > x_threshold) & (df[y_col] > y_threshold),
           '_color'] = up_color
    df.loc[(df[x_col] < -x_threshold) & (df[y_col] > y_threshold),
           '_color'] = down_color
    df.loc[df[y_col] < y_threshold, '_color'] = ns_color

    xmin = 0
    xmax = math.ceil(df[x_col].max() * 1.1)
    ymin = math.floor(df[y_col].min())
    ymax = math.ceil(df[y_col].max())

    if figsize is None:
        figsize = (9, 9)
    fig, ax = plt.subplots(figsize=figsize)
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    ax.scatter(df[x_col], df[y_col], s=point_size, c=df['_color'])
    ax.set_ylabel(ylabel, fontweight='bold',
                  fontsize=int(12 * font_scale))
    ax.set_xlabel(xlabel, fontweight='bold',
                  fontsize=int(12 * font_scale))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # 标注
    if highlight_labels is not None:
        to_label = df[df[label_col].isin(highlight_labels)]
    elif n_label is not None:
        df['_dist'] = (df[x_col] ** 2 + df[y_col] ** 2) ** 0.5
        to_label = df.nlargest(n_label, '_dist')
    else:
        to_label = pd.DataFrame()

    if len(to_label) > 0:
        texts = [
            ax.text(row[x_col] + 0.1, row[y_col], row[label_col],
                    fontsize=int(9 * font_scale), color='black',
                    style='italic', weight='bold',
                    va='center', ha='left')
            for _, row in to_label.iterrows()
        ]
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle='->', lw=0.9,
                                   color='green'),
                    force_explode=(-10, 2), iter_lim=10000,
                    avoid_self=True)

    # 参考线
    ax.vlines(x_threshold, ymin, ymax, color='dimgrey',
              linestyle='dashed', linewidth=1)
    ax.hlines(y_threshold, xmin, xmax, color='dimgrey',
              linestyle='dashed', linewidth=1)

    if output_path:
        fig.savefig(f'{output_path}.png', dpi=600, bbox_inches='tight')
        fig.savefig(f'{output_path}.pdf', dpi=300, bbox_inches='tight',
                    metadata={'Creator': None, 'Producer': None})
    plt.close(fig)
    return fig, ax
