#!/usr/bin/env python
# coding: utf-8
"""
代谢组学（按性别分组）pattern 点图 + KEGG/SMPDB 富集柱状图
富集方法: metabo_enrich_integrated (KEGG + SMPDB 本地通路)
"""
import pandas as pd
import numpy as np
from pathlib import Path
from tenacity import retry, stop_after_attempt, wait_fixed
import os

from drawing_yh.heatmap.heatmap_by_pattern.heatmap_by_pattern import draw_heatmap
from drawing_yh.hybrid_chart.age_dotmap_by_pattern import draw_pattern_pointplot
from metabo_enrich.metabo_enricher_integrated import metabo_enrich_integrated


base_path = Path(r"D:\Projects\Bone_Marrow_Aging\metabonomics\analysis")
feature_meta = pd.read_csv(base_path / 'data/feature_meta.csv', encoding='utf-8-sig', encoding_errors='replace')
protein_to_gene = dict(zip(feature_meta['COMP_ID'], feature_meta['COMPOUND_Name']))

protein_to_kegg = feature_meta.set_index('COMP_ID')['KEGG_ID'].dropna().to_dict()
protein_to_hmdb = feature_meta.set_index('COMP_ID')['HMDB_ID'].dropna().to_dict()


@retry(stop=stop_after_attempt(3), wait=wait_fixed(2))
def _enrich_genes(protein_ids):
    kegg_ids = [protein_to_kegg[p] for p in protein_ids if p in protein_to_kegg]
    hmdb_ids = [protein_to_hmdb[p] for p in protein_ids if p in protein_to_hmdb]

    if len(kegg_ids) + len(hmdb_ids) < 2:
        return None

    kegg_file = Path(r"D:\Projects\ref\enrichment_dataset\KEGG\hsa_pathways_with_compounds.txt")
    smpdb_dir = r"D:\Projects\ref\enrichment_dataset\SMPDB\smpdb_metabolites.csv"

    return metabo_enrich_integrated(
        kegg_c_list=kegg_ids,
        hmdb_list=hmdb_ids,
        kegg_txt_path=kegg_file,
        smpdb_dir=smpdb_dir
    )


EXPR_CSV = base_path / 'data/abundance_sample_x_protein.csv'
META_CSV = base_path / 'data/sample_meta.csv'
data_sub = 'by_sex'
PATTERN_CSV = base_path / f'data/{data_sub}/protein_correlation_summary.csv'
ENRICH_OUT_DIR = base_path / f'data/{data_sub}/enrich_results'

OUT_DIR = Path("outputs")
OUT_DIR.mkdir(exist_ok=True)

ENRICH_N_RESULTS = 20
ENRICH_N_DISPLAY = 2

FC_CUT = 0
BLANK_ROWS = 2
FIG_W, FIG_H = 10, 16
BASE_FONTSIZE = 20


def enrich_and_save(expr_df, df_pat, patterns, out_dir=ENRICH_OUT_DIR):
    os.makedirs(out_dir, exist_ok=True)
    enrich_dict = {}

    for i, pat in enumerate(patterns, 1):
        print(f"\n[{i}/{len(patterns)}] 处理 pattern: {pat}")

        genes_pat = df_pat[df_pat['pattern'] == pat].index.tolist()
        genes_pat = [g for g in genes_pat if g in expr_df.index]
        print(f"  ├─ 该pattern总蛋白数: {len(genes_pat)}")

        if not genes_pat:
            print(f"  └─ ❌ 跳过: 无匹配基因")
            enrich_dict[pat] = []
            continue

        all_genes = expr_df.loc[genes_pat].index.tolist()
        gene_names = [protein_to_gene.get(g, g) for g in all_genes]
        print(f"  ├─ 转换后基因名: {gene_names}")

        print(f"  ├─ 调用Enrichr API...")
        try:
            enrich_res = _enrich_genes(genes_pat)
            if enrich_res is None:
                print(f"  ├─ ⚠ API返回None")
                enrich_dict[pat] = []
                continue

            df_enrich = enrich_res.res2d
            if df_enrich.empty:
                print(f"  ├─ ⚠ 返回空结果")
                enrich_dict[pat] = []
                continue

            df_enrich = df_enrich.sort_values('P-value').head(ENRICH_N_RESULTS)
            df_enrich['neg_log10_p'] = -np.log10(df_enrich['P-value'])

            print(f"  ├─ 取前2个显著项:")
            for idx, row in df_enrich.iterrows():
                print(f"     - {row['Term']} (p={row['P-value']:.2e})")

            enrich_dict[pat] = df_enrich[['Term', 'P-value', 'Genes', 'neg_log10_p']].to_dict('records')

            fpath = f"{out_dir}/{pat.replace(' ', '_')}_enrich.csv"
            df_enrich.to_csv(fpath, index=False)
            print(f"  └─ ✅ 已保存: {fpath}")

        except Exception as e:
            print(f"  ├─ ❌ 错误: {type(e).__name__}: {str(e)}")
            enrich_dict[pat] = []

    print(f"\n=== 富集分析完成 ===")
    return enrich_dict


def load_or_enrich(expr_df, df_pat, patterns, out_dir=ENRICH_OUT_DIR):
    if os.path.exists(out_dir):
        csv_files = [f for f in os.listdir(out_dir) if f.endswith('_enrich.csv')]
        if len(csv_files) == len(patterns):
            print("✅ 检测到现存富集结果，直接读取")
            enrich_dict = {}
            for pat in patterns:
                fpath = os.path.join(out_dir, f"{pat.replace(' ', '_')}_enrich.csv")
                if os.path.exists(fpath):
                    df_enrich = pd.read_csv(fpath)
                    enrich_dict[pat] = df_enrich[['Term', 'P-value', 'Genes', 'neg_log10_p']].to_dict('records')
                else:
                    enrich_dict[pat] = []
            return enrich_dict

    print("❌ 未检测到完整富集结果，执行新的富集分析")
    return enrich_and_save(expr_df, df_pat, patterns, out_dir)


def main():
    expr_df = pd.read_csv(EXPR_CSV, index_col=0)
    meta = pd.read_csv(META_CSV, index_col="sample")
    df_pat = pd.read_csv(PATTERN_CSV)
    if 'gene' in df_pat.columns:
        df_pat = df_pat.set_index('gene')
    df_pat.index = df_pat.index.astype(str).str.replace(r'-\d+$', '', regex=True)
    expr_df.columns = expr_df.columns.astype(str).str.replace(r'-\d+$', '', regex=True)
    df_pat = df_pat[df_pat['pattern'] != 'Non-significant']

    trend_dict = df_pat.groupby('pattern')['overall_trend'].first().to_dict()

    def sort_key(p):
        trend = trend_dict.get(p, 'Mixed')
        trend_lower = trend.lower() if trend else 'mixed'
        return {'up': 0, 'down': 1, 'mixed': 2}.get(trend_lower, 3), p

    patterns = sorted(df_pat['pattern'].unique(), key=sort_key)

    condition_map = meta["condition"].str.capitalize()
    age_series = pd.to_numeric(meta["age"], errors="coerce")
    common_smp = expr_df.index.intersection(condition_map.index)
    expr_df = expr_df.loc[common_smp]
    condition_map = condition_map.loc[common_smp]
    age_series = age_series.reindex(common_smp)
    expr_df = expr_df.T

    print("=== 各 pattern 匹配蛋白数 ===")
    for pat in patterns:
        n_genes = len(df_pat[df_pat['pattern'] == pat])
        print(f"{pat:30s} {n_genes}")

    print("\n=== 检查富集结果 ===")
    enrich_dict = load_or_enrich(expr_df, df_pat, patterns)

    draw_pattern_pointplot(expr_df, df_pat, patterns, age_series, enrich_dict,
                           feature_meta, 'COMP_ID', 'COMPOUND_Name', 'pearson_r',
                           sample_meta=meta, figsize=(10, 16), base_fontsize=20,
                           fc_cut=0, out_dir=base_path / 'figure',
                           split_up_down=True, is_up_regulation=True, enrich_show_n=4)

    draw_pattern_pointplot(expr_df, df_pat, patterns, age_series, enrich_dict,
                           feature_meta, 'COMP_ID', 'COMPOUND_Name', 'pearson_r',
                           sample_meta=meta, figsize=(10, 16), base_fontsize=20,
                           fc_cut=0, out_dir=base_path / 'figure',
                           split_up_down=True, is_up_regulation=False, enrich_show_n=4)


if __name__ == "__main__":
    main()
