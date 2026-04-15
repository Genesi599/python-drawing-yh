import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

import pandas as pd
from drawing_yh.bar.waterfall import plot_waterfall

out_dir = os.path.dirname(__file__)

# ── 读取显著蛋白-年龄相关性数据 ─────────────────────────────────
csv_path = r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\by_sex\protein_age_correlations_significant_all.csv"
df = pd.read_csv(csv_path)

plot_waterfall(
    df,
    value_col='pearson_r',
    label_col='gene_name',
    pval_col='pearson_pval',
    title='Protein–Age Correlation (Significant)',
    xlabel='Pearson r',
    font_size=7,
    out_path=os.path.join(out_dir, 'waterfall.png'),
)

print('Done.')
