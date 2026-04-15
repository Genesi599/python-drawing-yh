from pathlib import Path
import pandas as pd
from enrichment_bar_steps import enrichment_tables   # 第一步函数

# ========== 唯一需要手动改的地方 ==========
root_dir = Path(r"D:\Projects\Neutrophil_Aging\NET_bulk\Ctrl_vs_NETs_vs_NETs_DNase_I")
# =========================================

def main():
    for csv_path in root_dir.rglob('*all.csv'):
        print(f'Processing: {csv_path}')
        df = pd.read_csv(csv_path)
        up_genes   = df.loc[df['type'] == 'Up', 'GeneName'].dropna().unique().tolist()
        down_genes = df.loc[df['type'] == 'Down', 'GeneName'].dropna().unique().tolist()

        # 用原文件名（去掉 .csv）当子目录名
        sub_dir = csv_path.parent / 'enrich_tables' / csv_path.stem
        sub_dir.mkdir(parents=True, exist_ok=True)

        paths = enrichment_tables(up_genes, down_genes,
                                  out_dir=sub_dir,
                                  max_term=300)
        print(f'  -> tables saved to: {sub_dir}\n')

if __name__ == '__main__':
    main()