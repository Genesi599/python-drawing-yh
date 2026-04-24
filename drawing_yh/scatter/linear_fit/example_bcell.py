import pandas as pd
from pathlib import Path
from linear_fit import plot_all_tissues
import numpy as np


def process_age_column(df, age_col='Age'):
    sample_val = df[age_col].dropna().iloc[0]
    if 'Y' in str(sample_val):
        age_unit = 'years'
        df['age'] = pd.to_numeric(df[age_col].str.replace('Y', ''), errors='coerce')
    elif 'M' in str(sample_val):
        age_unit = 'months'
        df['age'] = pd.to_numeric(df[age_col].str.replace('M', ''), errors='coerce')
    else:
        age_unit = 'years'
        df['age'] = pd.to_numeric(df[age_col], errors='coerce')

    df = df.dropna(subset=['age'])
    return df, age_unit


base_dir = Path(r"D:\Projects\B_Cell_Aging\B_cell_ratio")
output_dir = base_dir / "figure" / "linear_fit"
output_dir.mkdir(parents=True, exist_ok=True)

color_map_df = pd.read_csv(base_dir / "ref" / "tissue_system_mapping.csv")
tissue_colors = dict(zip(color_map_df['Tissue'], color_map_df['Color']))

for csv_file in base_dir.glob("*_Bcell_sample_summary_with_counts.csv"):
    species = csv_file.stem.split('_')[0]

    cor_file = base_dir / "correlation_results" / f"{species}_tissue_age_correlation_results.csv"
    if not cor_file.exists():
        print(f"跳过 {species}: 未找到相关性文件")
        continue

    donor_summary = pd.read_csv(csv_file)
    cor_df = pd.read_csv(cor_file)

    donor_summary, age_unit = process_age_column(donor_summary)

    cor_df_plot = cor_df.rename(columns={
        'Tissue': 'tissue_general',
        'Spearman_Correlation': 'spearman_r',
        'Pearson_Correlation': 'pearson_r',
        'Spearman_P_value': 'p_value'
    })

    donor_summary_plot = donor_summary.rename(columns={'Tissue': 'tissue_general'})

    # plot_all_tissues(
    #     cor_df=cor_df_plot,
    #     donor_summary=donor_summary_plot,
    #     exclude_zero=False,
    #     cell_ratio_col="b_cell_ratio",
    #     age_col="age",
    #     age_unit=age_unit,
    #     xlim_margin=2,
    #     cell_type="B Cell",
    #     output_path=str(output_dir / f"{species}_all_tissues"),
    #     cell_count_col=None
    # )

    filtered_cor = cor_df_plot[(cor_df_plot['spearman_r'] > 0) & (cor_df_plot['Sample_count'] >= 20)].copy()
    neg_log_p = -np.log10(filtered_cor['p_value'])
    r_norm = filtered_cor['spearman_r'] / filtered_cor['spearman_r'].max()
    p_norm = neg_log_p / neg_log_p.max()
    filtered_cor['distance'] = np.sqrt(r_norm ** 2 + p_norm ** 2)
    filtered_cor = filtered_cor.sort_values('distance', ascending=False).head(4).drop(columns=['distance'])

    if len(filtered_cor) > 0:
        top_tissues = filtered_cor['tissue_general'].tolist()
        top_colors = {t: tissue_colors.get(t, '#808080') for t in top_tissues}

        plot_all_tissues(
            cor_df=filtered_cor,
            donor_summary=donor_summary_plot,
            exclude_zero=False,
            cell_ratio_col="b_cell_ratio",
            age_col="age",
            age_unit=age_unit,
            xlim_margin=2,
            cell_type="B Cell",
            output_path=str(output_dir / f"{species}_top4_tissues"),
            cell_count_col=None,
            tissue_colors=top_colors,
            font_scale=2.2
        )

    print(f"完成 {species}")