#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : B_cell_linear_fit.py
@Date    : 2026/2/6 12:43
@Author  : yh109
"""
from heatmap_tile_style import *

def prepare_data_from_summary(input_file):
    df = pd.read_csv(input_file)

    df['geneset'] = df['directory']
    df['tissue_category'] = df['file'].str.split('_').str[0]
    df['correlation'] = df['score_mean_effect_size']
    df['pvalue'] = df['score_mean_pvalue']

    return df[['geneset', 'tissue_category', 'correlation', 'pvalue']]


if __name__ == '__main__':
    input_file = "D:/Projects/Neutrophil_Aging/spatial/plots_python/pvalue_summary.csv"
    output_dir = "D:/Projects/Neutrophil_Aging/spatial/plots_python"

    df = prepare_data_from_summary(input_file)
    temp_file = os.path.join(output_dir, "temp_formatted.csv")
    df.to_csv(temp_file, index=False)

    geneset_order = [
        "Neutrophil_surface_protein",
        "Neutrophil_activation",
        "Neutrophil_degranulation",
        "Primary_granule",
        "Secondary_granule",
        "Tertiary_granule",
        "Neutrophil_extracellular_trap",
        "Neutrophil_extracellular_trap_new",
        "Neutrophil_chemotaxis",
        "Neutrophil_migration"
    ]

    plot_correlation_heatmap(temp_file, output_dir, geneset_order=geneset_order, cbar_label="Effect Size")
    os.remove(temp_file)

