# batch_enrich.py
from pathlib import Path
import pandas as pd
from enrichment_bar import enrichment_combo   # 关键函数导入

# ========== 唯一需要手动改的地方 ==========
root_dir = Path(r"D:\Projects\Neutrophil_Aging\NET_bulk\Ctrl_vs_NETs_vs_NETs_DNase_I")
# =========================================

def main():
    for csv_path in root_dir.rglob('*all.csv'):
        print(f'Processing: {csv_path}')
        df = pd.read_csv(csv_path)
        up_genes   = df.loc[df['type'] == 'Up', 'GeneName'].dropna().unique().tolist()
        down_genes = df.loc[df['type'] == 'Down', 'GeneName'].dropna().unique().tolist()

        out_path = csv_path.with_suffix('')  # 去掉 .csv
        enrichment_combo(up_genes, down_genes,
                         gene_fontcolor_up='red',
                         gene_fontcolor_down='blue',
                         combo_path=str(out_path) + '_enrich.png')
        print(f'  -> saved: {out_path}_enrich.png\n')

if __name__ == '__main__':
    main()