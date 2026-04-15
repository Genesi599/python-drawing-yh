# enrichment_combo.py
import os, gseapy as gp, numpy as np, matplotlib.pyplot as plt
from tenacity import retry, stop_after_attempt, wait_fixed
import pandas as pd
import re

# ---------- 重试 ----------
@retry(stop=stop_after_attempt(5), wait=wait_fixed(3))
def _enrich_with_retry(gene_list, gene_set):
    return gp.enrichr(gene_list=gene_list,
                      gene_sets=gene_set,
                      organism='Human',
                      outdir=None)

# ---------- 单幅富集图 ----------


import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'pdf.use14corefonts': False,
    'font.family': 'Arial',
    'font.size': 8,
})

FONT_SIZE = 8
BAR_HEIGHT = 0.5
ROW_HEIGHT_IN = 0.28
FIG_WIDTH_IN = 6.89




def _single_enrich(analysis_type, gene_fontcolor, bar_color,
                   df=None, file=None, n_term=8, alpha=0.75):
    if file is not None:
        df = pd.read_csv(file, sep='\t') if str(file).endswith('.txt') else pd.read_csv(file)
    if df is None:
        raise ValueError("df 和 file 至少提供一个")
    if analysis_type not in ('GO', 'KEGG'):
        raise ValueError("analysis_type must be 'GO' or 'KEGG'")

    df = df.copy()
    df['neg_log10_p'] = -np.log10(df['P-value'])
    df = df.sort_values('neg_log10_p').tail(n_term)

    fig, ax = plt.subplots(figsize=(FIG_WIDTH_IN, len(df) * ROW_HEIGHT_IN + 0.6))
    ax.barh(np.arange(len(df)), df['neg_log10_p'], color=bar_color, height=BAR_HEIGHT, alpha=alpha)

    for idx, (term, genes) in enumerate(zip(df['Term'], df['Genes'])):
        ax.text(0.1, idx, re.sub(r'\s*\(GO:\d+\)', '', term),
                va='center', ha='left', color='black', fontsize=FONT_SIZE)
        ax.text(0.1, idx - 0.28, f'({", ".join(genes.split(";")[:6])})',
                va='top', ha='left', fontsize=FONT_SIZE, color=gene_fontcolor, fontstyle='italic')

    ax.set_xlabel('-Log10(p-value)', fontsize=FONT_SIZE)
    ax.tick_params(axis='x', labelsize=FONT_SIZE, pad=8)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines['bottom'].set_position(('outward', 8))
    ax.set_yticks([])

    fig.canvas.draw()
    fig.set_size_inches(
        max(t.get_window_extent(fig.canvas.get_renderer()).width for t in ax.texts) / fig.dpi + 0.5,
        len(df) * ROW_HEIGHT_IN + 0.6
    )
    fig.tight_layout()
    return fig, ax, df



def enrich_up_down(analysis_type,
                   up_file, down_file,
                   up_color='red', down_color='blue',
                   up_fontcolor='black', down_fontcolor='black',
                   n_term_up=8, n_term_down=8,
                   alpha=0.75,
                   up_title='', down_title=''):
    def _prep(file, n_term):
        df = pd.read_csv(file, sep='\t') if str(file).endswith('.txt') else pd.read_csv(file)
        df = df.copy()
        df['neg_log10_p'] = -np.log10(df['P-value'])
        return df.sort_values('neg_log10_p').tail(n_term)

    df_up = _prep(up_file, n_term_up)
    df_down = _prep(down_file, n_term_down)

    h_up = len(df_up) * ROW_HEIGHT_IN + 0.6
    h_down = len(df_down) * ROW_HEIGHT_IN + 0.6
    fig = plt.figure(figsize=(FIG_WIDTH_IN, h_up + h_down))
    gs = fig.add_gridspec(2, 1, height_ratios=[len(df_up), len(df_down)], hspace=0.4)
    ax_up = fig.add_subplot(gs[0])
    ax_down = fig.add_subplot(gs[1])

    def _draw(ax, df, bar_color, gene_fontcolor, title='', show_xlabel=True):
        ax.barh(np.arange(len(df)), df['neg_log10_p'], color=bar_color, height=BAR_HEIGHT, alpha=alpha)
        for idx, (term, genes) in enumerate(zip(df['Term'], df['Genes'])):
            ax.text(0.1, idx, re.sub(r'\s*\(GO:\d+\)', '', term),
                    va='center', ha='left', color='black', fontsize=FONT_SIZE)
            ax.text(0.1, idx - 0.28, f'({", ".join(genes.split(";")[:6])})',
                    va='top', ha='left', fontsize=FONT_SIZE, color=gene_fontcolor, fontstyle='italic')
        if show_xlabel:
            ax.set_xlabel('-Log10(p-value)', fontsize=FONT_SIZE)
        else:
            ax.tick_params(axis='x', bottom=False, labelbottom=False)
            ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='x', labelsize=FONT_SIZE, pad=8)
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines['bottom'].set_position(('outward', 8))
        ax.set_yticks([])
        if title:
            ax.set_title(title, fontsize=FONT_SIZE, fontweight='normal', loc='left', pad=4)

    _draw(ax_up, df_up, up_color, up_fontcolor, title=up_title, show_xlabel=False)
    _draw(ax_down, df_down, down_color, down_fontcolor, title=down_title, show_xlabel=True)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    max_w = max(t.get_window_extent(renderer).width for t in list(ax_up.texts) + list(ax_down.texts))
    fig.set_size_inches(max_w / fig.dpi + 0.5, h_up + h_down)
    fig.tight_layout()
    return fig, ax_up, ax_down, df_up, df_down



# ---------- 四合一主函数 ----------
def enrichment_combo(up_genes, down_genes,
                     gene_fontcolor_up='red',
                     gene_fontcolor_down='blue',
                     combo_path='enrich.png'):
    """
    输入 up/down 基因 list，直接保存 2×2 四宫格
    """
    color_map = {
        ('GO',  'up'): '#FFB38A',
        ('KEGG','up'): '#FF8E72',
        ('GO',  'dn'): '#81C7F4',
        ('KEGG','dn'): '#69d4d4',
    }
    jobs = [
        (up_genes,   'GO',  'up',  gene_fontcolor_up,  color_map[('GO',  'up')]),
        (up_genes,   'KEGG','up',  gene_fontcolor_up,  color_map[('KEGG','up')]),
        (down_genes, 'GO',  'dn',  gene_fontcolor_down,color_map[('GO',  'dn')]),
        (down_genes, 'KEGG','dn',  gene_fontcolor_down,color_map[('KEGG','dn')]),
    ]

    axes_data = []
    for genes, atype, ud, gfont, bcol in jobs:
        if not genes:
            genes = ['DUMMY']          # 空列表占位
        fig, ax, _ = _single_enrich(genes, atype, gfont, bcol)
        axes_data.append((ax, atype, ud))
        plt.close(fig)

    # 拼 2×2
    combo_fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    pos_map  = {('up', 'GO'): (0, 0), ('up', 'KEGG'): (0, 1),
                ('dn', 'GO'): (1, 0), ('dn', 'KEGG'): (1, 1)}
    title_map= {('up', 'GO'): 'Up-regulated GO',
                ('up', 'KEGG'): 'Up-regulated KEGG',
                ('dn', 'GO'): 'Down-regulated GO',
                ('dn', 'KEGG'): 'Down-regulated KEGG'}

    for (ax_src, atype, ud), (_, _, _, gfont, _) in zip(axes_data, jobs):
        row, col = pos_map[(ud, atype)]
        ax = axes[row, col]

        # bar
        for patch in ax_src.patches:
            ax.barh(patch.get_y() + patch.get_height()/2,
                    patch.get_width(), height=patch.get_height(),
                    color=patch.get_facecolor(), alpha=patch.get_alpha())
        # text 统一字号
        TERM_FS, GENE_FS = 12, 11
        for txt in ax_src.texts:
            t = txt.get_text()
            fs = GENE_FS if t.startswith('(') else TERM_FS
            ax.text(txt.get_position()[0], txt.get_position()[1], t,
                    fontsize=fs, color=txt.get_color(), fontweight='bold',
                    fontstyle='italic' if t.startswith('(') else 'normal',
                    va=txt.get_va(), ha=txt.get_ha())

        ax.set_xlabel('-Log10(p-value)', fontsize=13, fontweight='bold')
        ax.set_title(title_map[(ud, atype)], pad=10, loc='left',
                     fontsize=15, fontweight='bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(left=False, labelleft=False)

    combo_fig.tight_layout()
    combo_fig.savefig(combo_path, dpi=600, bbox_inches='tight')
    combo_fig.savefig(combo_path.replace('.png', '.pdf'),
                      dpi=600, bbox_inches='tight')
    plt.close(combo_fig)