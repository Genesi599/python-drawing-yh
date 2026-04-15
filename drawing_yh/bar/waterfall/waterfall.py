"""
双栏瀑布图：展示蛋白/基因与某连续变量（如年龄）的相关性
左栏 = 正相关（降序），右栏 = 负相关（升序）
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.family']  = 'Arial'


def plot_waterfall(
    df: pd.DataFrame,
    value_col: str        = 'pearson_r',
    label_col: str        = 'gene_name',
    pval_col: str | None  = 'pearson_pval',
    pval_thresh: float    = 0.05,
    cmap_pos: tuple       = ('#fde0dd', '#c51b8a'),
    cmap_neg: tuple       = ('#deebf7', '#2171b5'),
    figsize: tuple | None = None,
    font_size: float      = 7,
    bar_height: float     = 0.7,
    title: str            = '',
    xlabel: str           = '',
    out_path: str | None  = None,
    dpi: int              = 600,
):
    """
    Parameters
    ----------
    df : DataFrame
        至少含 value_col 和 label_col 两列。
    value_col : str
        用于排序和绘制条形的数值列（如 pearson_r、slope）。
    label_col : str
        Y 轴标签列（如 gene_name）。
    pval_col : str | None
        p 值列名。若提供，则在标签后标星号。
    pval_thresh : float
        标记星号的 p 值阈值。
    cmap_pos / cmap_neg : tuple
        正/负相关的渐变色范围 (浅, 深)。
    figsize : tuple | None
        图像尺寸，None 则自动计算。
    out_path : str | None
        输出路径（自动保存 png/pdf/svg）。
    """
    df = df.copy()
    df[value_col] = pd.to_numeric(df[value_col], errors='coerce')
    df = df.dropna(subset=[value_col])

    pos = df[df[value_col] > 0].sort_values(value_col, ascending=False)
    neg = df[df[value_col] < 0].sort_values(value_col, ascending=True)

    n_pos, n_neg = len(pos), len(neg)
    n_max = max(n_pos, n_neg, 1)

    if figsize is None:
        col_w = 3.2
        row_h = max(font_size * 1.6 / 72, 0.16)
        figsize = (col_w * 2 + 0.6, max(n_max * row_h + 0.8, 2.0))

    fig, (ax_pos, ax_neg) = plt.subplots(
        1, 2, figsize=figsize, facecolor='white',
        gridspec_kw={'wspace': 0.02},
    )

    cm_pos = LinearSegmentedColormap.from_list('pos', cmap_pos, N=256)
    cm_neg = LinearSegmentedColormap.from_list('neg', cmap_neg, N=256)

    def _draw_panel(ax, sub_df, cmap, is_pos):
        n = len(sub_df)
        if n == 0:
            ax.set_visible(False)
            return

        vals   = sub_df[value_col].values
        labels = sub_df[label_col].astype(str).values
        abs_max = max(abs(vals.max()), abs(vals.min()), 1e-9)
        normed = np.abs(vals) / abs_max
        colors = [cmap(v) for v in normed]

        y = np.arange(n)
        ax.barh(y, vals, height=bar_height, color=colors,
                edgecolor='white', linewidth=0.3)

        # 星号标记
        stars = [''] * n
        if pval_col and pval_col in sub_df.columns:
            pvals = sub_df[pval_col].values
            for i, p in enumerate(pvals):
                if pd.notna(p):
                    if p < 0.001:
                        stars[i] = ' ***'
                    elif p < 0.01:
                        stars[i] = ' **'
                    elif p < pval_thresh:
                        stars[i] = ' *'

        label_strs = [f'{l}{s}' for l, s in zip(labels, stars)]

        if is_pos:
            ax.set_yticks(y)
            ax.set_yticklabels(label_strs, fontsize=font_size, fontstyle='italic')
            ax.yaxis.tick_right()
            ax.invert_xaxis()
        else:
            ax.set_yticks(y)
            ax.set_yticklabels(label_strs, fontsize=font_size, fontstyle='italic')
            ax.yaxis.tick_left()

        ax.invert_yaxis()
        ax.set_ylim(n_max - 0.5, -0.5)
        ax.tick_params(axis='y', length=0, pad=3)
        ax.spines[['top', 'right', 'left']].set_visible(False)
        ax.tick_params(axis='x', labelsize=font_size)

    _draw_panel(ax_pos, pos, cm_pos, is_pos=True)
    _draw_panel(ax_neg, neg, cm_neg, is_pos=False)

    # X 轴标签
    if xlabel:
        fig.text(0.5, 0.01, xlabel, ha='center', fontsize=font_size + 1,
                 fontfamily='Arial')

    # 标题
    if n_pos > 0:
        ax_pos.set_title(f'Positive (n={n_pos})', fontsize=font_size + 1,
                         fontweight='bold', pad=6)
    if n_neg > 0:
        ax_neg.set_title(f'Negative (n={n_neg})', fontsize=font_size + 1,
                         fontweight='bold', pad=6)

    if title:
        fig.suptitle(title, fontsize=font_size + 3, fontweight='bold', y=0.99)

    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.05, top=0.93,
                        wspace=0.02)

    if out_path:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        for suffix in ['.png', '.pdf', '.svg']:
            p = out_path.with_suffix(suffix)
            fig.savefig(p, dpi=dpi if suffix == '.png' else None,
                        bbox_inches='tight', facecolor='white')
            print(f'Saved: {p}')

    return fig, (ax_pos, ax_neg)
