# enrichment_combo.py
import os, gseapy as gp, numpy as np, matplotlib.pyplot as plt
from tenacity import retry, stop_after_attempt, wait_fixed

# ---------- 重试 ----------
@retry(stop=stop_after_attempt(5), wait=wait_fixed(3))
def _enrich_with_retry(gene_list, gene_set):
    return gp.enrichr(gene_list=gene_list,
                      gene_sets=gene_set,
                      organism='Human',
                      outdir=None)

# ---------- 单幅富集图 ----------
def _single_enrich(gene_list, analysis_type, gene_fontcolor,
                   bar_color, n_term=8):
    """返回画好的 fig, ax, df（仅供内部组合）"""
    if analysis_type == 'GO':
        gene_set = 'GO_Biological_Process_2025'
        title_prefix = 'GO Biological Process'
    elif analysis_type == 'KEGG':
        gene_set = 'KEGG_2021_Human'
        title_prefix = 'KEGG Pathway'
    else:
        raise ValueError("analysis_type must be 'GO' or 'KEGG'")

    res = _enrich_with_retry(gene_list, gene_set)
    df = res.results
    df['neg_log10_p'] = -np.log10(df['P-value'])
    df = df.sort_values('neg_log10_p').tail(n_term)  # 前 8

    fig, ax = plt.subplots(figsize=(7, 4))
    y_pos = np.arange(len(df))
    ax.barh(y_pos, df['neg_log10_p'],
            color=bar_color, height=0.6, alpha=0.75)

    # 标注
    for idx, (term, genes, x) in enumerate(
            zip(df['Term'], df['Genes'], df['neg_log10_p'])):
        ax.text(0.1, idx, term,
                va='center', ha='left',
                color='black', fontsize=11, fontweight='bold')

        gene_txt = ', '.join(genes.split(';')[:6])
        ax.text(0.1, idx - 0.3, f'({gene_txt})',
                va='top', ha='left',
                fontsize=11, color=gene_fontcolor,
                fontweight='bold', fontstyle='italic')

    ax.set_xlabel('-Log10(p-value)', fontsize=13, fontweight='bold')
    ax.set_title(title_prefix, pad=10, loc='left', fontsize=15, fontweight='bold')
    ax.spines.top.set_visible(False)
    ax.spines.right.set_visible(False)
    return fig, ax, df


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
    combo_fig, axes = plt.subplots(2, 2, figsize=(16, 9))
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
    plt.close(combo_fig)