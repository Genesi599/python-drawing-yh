import os, gseapy as gp, numpy as np, pandas as pd, matplotlib.pyplot as plt
from tenacity import retry, stop_after_attempt, wait_fixed

# ---------- 通用重试 ----------
@retry(stop=stop_after_attempt(5), wait=wait_fixed(3))
def _enrich_with_retry(gene_list, gene_set):
    return gp.enrichr(gene_list=gene_list,
                      gene_sets=gene_set,
                      organism='Human',
                      outdir=None)

# ---------- 第一步：富集 → CSV ----------
def enrichment_tables(up_genes, down_genes,
                      out_dir='enrich_tables',
                      max_term=50):
    """
    返回 4 个 CSV 路径列表（up/go, up/kegg, down/go, down/kegg）
    每个表按 p 值升序，最多 max_term 行
    """
    os.makedirs(out_dir, exist_ok=True)
    gene_sets = {
        'GO': 'GO_Biological_Process_2025',
        'KEGG': 'KEGG_2021_Human'
    }
    jobs = [
        (up_genes,   'GO',   'up'),
        (up_genes,   'KEGG', 'up'),
        (down_genes, 'GO',   'dn'),
        (down_genes, 'KEGG', 'dn'),
    ]
    paths = []
    for genes, atype, ud in jobs:
        if not genes:
            genes = ['DUMMY']
        res = _enrich_with_retry(genes, gene_sets[atype])
        df = res.results
        df['neg_log10_p'] = -np.log10(df['P-value'])
        df = df.sort_values('P-value').head(max_term)
        fpath = os.path.join(out_dir, f'{ud}_{atype}.csv')
        df.to_csv(fpath, index=False)
        paths.append(fpath)
    return paths

# ---------- 第二步：从 CSV → 2×2 图 ----------
def plot_from_tables(table_dir='enrich_tables',
                     gene_fontcolor_up='red',
                     gene_fontcolor_down='blue',
                     combo_path='enrich.png',
                     top_n=8):
    """
    读取 4 个表，按 neg_log10_p 取前 top_n 行画图
    """
    color_map = {
        ('GO',  'up'): '#FFB38A',
        ('KEGG','up'): '#FF8E72',
        ('GO',  'dn'): '#81C7F4',
        ('KEGG','dn'): '#69d4d4',
    }
    title_map = {
        ('up', 'GO'):   'Up-regulated GO',
        ('up', 'KEGG'): 'Up-regulated KEGG',
        ('dn', 'GO'):   'Down-regulated GO',
        ('dn', 'KEGG'): 'Down-regulated KEGG',
    }
    pos_map = {
        ('up', 'GO'):   (0, 0),
        ('up', 'KEGG'): (0, 1),
        ('dn', 'GO'):   (1, 0),
        ('dn', 'KEGG'): (1, 1),
    }

    # 一次性读 4 个表
    df_dict = {}
    for ud, atype in pos_map.keys():
        fpath = os.path.join(table_dir, f'{ud}_{atype}.csv')
        df = pd.read_csv(fpath)
        df = df.sort_values('neg_log10_p', ascending=False).head(top_n)
        df_dict[(ud, atype)] = df

    # 开始画图
    combo_fig, axes = plt.subplots(2, 2, figsize=(16, 9))
    for (ud, atype), df in df_dict.items():
        row, col = pos_map[(ud, atype)]
        ax = axes[row, col]
        bar_color = color_map[(atype, ud)]
        y_pos = np.arange(len(df))

        # bar
        ax.barh(y_pos, df['neg_log10_p'],
                color=bar_color, height=0.6, alpha=0.75)

        # 文字
        TERM_FS, GENE_FS = 12, 11
        for idx, (term, genes, x) in enumerate(
                zip(df['Term'], df['Genes'], df['neg_log10_p'])):
            ax.text(0.1, idx, term,
                    va='center', ha='left',
                    fontsize=TERM_FS, color='black', fontweight='bold')

            gene_txt = ', '.join(genes.split(';')[:6])
            ax.text(0.1, idx - 0.3, f'({gene_txt})',
                    va='top', ha='left',
                    fontsize=GENE_FS,
                    color=gene_fontcolor_up if ud == 'up' else gene_fontcolor_down,
                    fontweight='bold', fontstyle='italic')

        ax.set_xlabel('-Log10(p-value)', fontsize=13, fontweight='bold')
        ax.set_title(title_map[(ud, atype)], pad=10, loc='left',
                     fontsize=15, fontweight='bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(left=False, labelleft=False)

    combo_fig.tight_layout()
    combo_fig.savefig(combo_path, dpi=600, bbox_inches='tight')
    plt.close(combo_fig)

