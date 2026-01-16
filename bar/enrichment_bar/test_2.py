import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from gseapy import enrichr
from tenacity import retry, stop_after_attempt, wait_fixed

# ---------- 通用重试 ----------
@retry(stop=stop_after_attempt(5), wait=wait_fixed(3))
def _enrich_with_retry(gene_list, gene_set, outdir=None):
    return enrichr(gene_list=gene_list, gene_sets=gene_set, organism='Human', outdir=outdir, cutoff=1e-5)

# ---------- 读取输入文件并按 logfoldchanges 收集上下调基因 ----------
def load_genes_from_table(csv_path, logfc_col='logfoldchanges', gene_col='names'):
    """
    从输入表中返回两组基因名列表：
    - up_genes: logfc > 0 的基因
    - down_genes: logfc < 0 的基因
    仅使用 names 列中的基因名
    """
    df = pd.read_csv(csv_path, sep=',')
    if gene_col not in df.columns or logfc_col not in df.columns:
        raise ValueError(f"输入表需要包含列: {gene_col} 和 {logfc_col}")

    df = df[[gene_col, logfc_col]].rename(columns={gene_col: 'names', logfc_col: 'logfc'})

    up_genes = df.loc[df['logfc'] > 0, 'names'].dropna().astype(str).tolist()
    down_genes = df.loc[df['logfc'] < 0, 'names'].dropna().astype(str).tolist()

    return up_genes, down_genes

# ---------- 第一步：富集 → Excel（输出到同目录，四个 sheet） ----------
def enrichment_to_excel(up_genes, down_genes, base_out_dir,
                        max_term=100, excel_filename='enrich_results.xlsx'):
    """
    将四个富集结果写入同一个 Excel 文件的四个 sheet：
    - up_GO, up_KEGG, dn_GO, dn_KEGG
    每个表按 p 值升序，最多 max_term 行
    """
    os.makedirs(base_out_dir, exist_ok=True)
    gene_sets = {
        'GO': 'GO_Biological_Process_2025',
        'KEGG': 'KEGG_2021_Human'
    }
    jobs = [
        (up_genes,   'GO', 'up'),
        (up_genes,   'KEGG', 'up'),
        (down_genes, 'GO',   'dn'),
        (down_genes, 'KEGG', 'dn'),
    ]

    df_map = {}
    for genes, atype, ud in jobs:
        if not genes:
            genes = ['DUMMY']
        res = _enrich_with_retry(genes, gene_sets[atype], outdir=None)
        df = res.results
        # 找到 P 值列的名称
        if 'P-value' in df.columns:
            pcol = 'P-value'
        else:
            pcol = None
            for cname in ['P-value', 'pvalue', 'p.value', 'pval']:
                if cname in df.columns:
                    pcol = cname
                    break
        if pcol is None:
            raise KeyError("P 值列未在 Enrichr 结果中找到，请检查列名。")
        df['neg_log10_p'] = -np.log10(df[pcol].astype(float))
        df = df.sort_values(pcol, ascending=True).head(max_term)
        df_map[(ud, atype)] = df

    excel_path = os.path.join(base_out_dir, excel_filename)
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df_map[('up', 'GO')].to_excel(writer, sheet_name='up_GO', index=False)
        df_map[('up', 'KEGG')].to_excel(writer, sheet_name='up_KEGG', index=False)
        df_map[('dn', 'GO')].to_excel(writer, sheet_name='dn_GO', index=False)
        df_map[('dn', 'KEGG')].to_excel(writer, sheet_name='dn_KEGG', index=False)

    return excel_path

# ---------- 从 Excel 读取并画图（2x2） ----------
def plot_from_excel(excel_path,
                    top_n=8,
                    combo_path=None,
                    title_prefix=''):
    """
    从 Excel 的四个 sheet 读取数据，按 P 值升序取前 top_n 行画图
    输出到与 Excel 文件同目录的 combo_path
    Sheet 名称：up_GO、up_KEGG、dn_GO、dn_KEGG
    """
    # 读取四个 sheet
    sheets = {
        ('up', 'GO'): 'up_GO',
        ('up', 'KEGG'): 'up_KEGG',
        ('dn', 'GO'): 'dn_GO',
        ('dn', 'KEGG'): 'dn_KEGG',
    }
    df_dict = {}
    for key, sheet in sheets.items():
        df = pd.read_excel(excel_path, sheet_name=sheet)
        # 统一使用存在的 P 值列进行排序
        if 'P-value' in df.columns:
            pcol = 'P-value'
        else:
            pcol = None
            for cname in ['P-value', 'pvalue', 'p.value', 'pval']:
                if cname in df.columns:
                    pcol = cname
                    break

        if pcol is not None:
            df = df.sort_values(pcol, ascending=True).head(top_n)
        else:
            df = df.head(top_n)
        df_dict[key] = df

    # 颜色与文本布局保持原有逻辑
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

    # 开始画图
    combo_fig, axes = plt.subplots(2, 2, figsize=(16, 9))
    for (ud, atype), df in df_dict.items():
        row, col = pos_map[(ud, atype)]
        ax = axes[row, col]
        bar_color = color_map.get((atype, ud), '#888888')
        y_pos = np.arange(len(df))[::-1]

        # bar
        ax.barh(y_pos, df['neg_log10_p'], color=bar_color, height=0.6, alpha=0.75)

        # 文字
        TERM_FS, GENE_FS = 12, 11
        for idx, (term, genes, x) in enumerate(
                zip(df['Term'], df['Genes'], df['neg_log10_p'])):
            ax.text(0.1, idx, term,
                    va='center', ha='left',
                    fontsize=TERM_FS, color='black', fontweight='bold')

            gene_txt = ', '.join(str(g) for g in genes.split(';')[:6])
            ax.text(0.1, idx - 0.3, f'({gene_txt})',
                    va='top', ha='left',
                    fontsize=GENE_FS,
                    color=('#FF0000' if ud == 'up' else '#0000FF'),
                    fontweight='bold', fontstyle='italic')

        ax.set_xlabel('-Log10(p-value)', fontsize=13, fontweight='bold')
        ax.set_title(title_map[(ud, atype)], pad=10, loc='left',
                     fontsize=15, fontweight='bold')
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(left=False, labelleft=False)

    combo_fig.tight_layout()
    if combo_path is None:
        combo_path = os.path.join(os.path.dirname(excel_path), 'enrich.png')
    combo_fig.savefig(combo_path, dpi=600, bbox_inches='tight')
    plt.close(combo_fig)

    return combo_path

# ---------- 主流程 ----------
def main(input_csv_path):
    # 1) 输入文件所在目录，作为输出目录
    input_dir = os.path.dirname(os.path.abspath(input_csv_path))

    # 2) 从输入表中读取 tissue 列的唯一值
    df0 = pd.read_csv(input_csv_path, sep=',')
    if 'tissue' in df0.columns:
        tissue_col = 'tissue'
    elif 'Tissue' in df0.columns:
        tissue_col = 'Tissue'
    else:
        raise ValueError("输入 CSV 中未找到 tissue 列，请确认列名为 tissue/Tissue。")

    tissues = sorted(df0[tissue_col].dropna().unique().tolist())

    results = {}
    for t in tissues:
        # 2a) 过滤当前 tissue 的基因
        df_t = df0[df0[tissue_col] == t]
        if df_t.empty:
            continue

        up_genes = df_t.loc[df_t['logfoldchanges'] > 0, 'names'].dropna().astype(str).tolist()
        down_genes = df_t.loc[df_t['logfoldchanges'] < 0, 'names'].dropna().astype(str).tolist()

        # 2b) 为当前 tissue 富集分析，输出到单独的 Excel（四个 sheet）
        enrich_excel_path = enrichment_to_excel(
            up_genes, down_genes,
            base_out_dir=input_dir,
            max_term=50,
            excel_filename=f'enrich_results_{t}.xlsx'
        )

        # 2c) 生成该 tissue 的图像
        combo_path = os.path.join(input_dir, f'enrich_{t}.png')
        plot_from_excel(enrich_excel_path, top_n=8, combo_path=combo_path)

        results[t] = {
            'enrich_excel': enrich_excel_path,
            'plot_png': combo_path
        }

    return results



# 如果直接作为脚本运行
if __name__ == '__main__':
    # 请替换为你的输入文件路径
    input_csv = r"C:\Users\yh599\Desktop\test.csv"
    result = main(input_csv)
    print(result)