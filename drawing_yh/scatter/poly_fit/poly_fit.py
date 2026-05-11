#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
多项式拟合散点图
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.use14corefonts'] = False
mpl.rcParams['font.family'] = 'Arial'


def plot_poly_fit(x, y, degree=4, xlabel='Age', ylabel='Value',
                  figsize=(8, 6), ci_scale=0.3, point_size=10,
                  font_scale=1.0, output_path=None):
    """
    绘制多项式拟合散点图，带置信区间。

    Parameters
    ----------
    x, y : array-like
        横纵坐标数据
    degree : int
        多项式阶数
    xlabel, ylabel : str
        坐标轴标签
    figsize : tuple
        图片大小
    ci_scale : float
        置信区间宽度系数（相对于残差标准差的倍数）
    point_size : int
        散点大小
    font_scale : float
        字体缩放
    output_path : str, optional
        输出路径（不含扩展名），同时保存 png 和 pdf
    """
    x = np.asarray(x)
    y = np.asarray(y)

    coeffs = np.polyfit(x, y, degree)
    p = np.poly1d(coeffs)
    x_fit = np.linspace(x.min(), x.max(), 100)
    y_fit = p(x_fit)

    y_pred = p(x)
    residuals = y - y_pred
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(x, y, s=point_size, label='Data Points')
    ax.plot(x_fit, y_fit, color='red',
            label=f'R² = {r_squared:.3f}')
    ax.fill_between(x_fit,
                     y_fit - ci_scale * np.std(residuals),
                     y_fit + ci_scale * np.std(residuals),
                     color='gray', alpha=0.2)
    ax.set_xlabel(xlabel, fontsize=int(20 * font_scale))
    ax.set_ylabel(ylabel, fontsize=int(20 * font_scale))
    ax.tick_params(axis='both', which='major',
                   labelsize=int(15 * font_scale))
    ax.legend(fontsize=int(12 * font_scale))

    if output_path:
        fig.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
        fig.savefig(f'{output_path}.pdf', dpi=300, bbox_inches='tight',
                    metadata={'Creator': None, 'Producer': None})
    plt.close(fig)
    return fig, ax, r_squared
