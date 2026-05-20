#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : heatmap_tile_style.py
@Date    : 2026/2/3 15:40
@Author  : yh109
"""
import pandas as pd
import os
from pathlib import Path

from drawing_yh.heatmap import save_long_heatmap

def plot_correlation_heatmap(input_file, output_dir, geneset_order=None, pvalue_col='pvalue',
                             font_scale=3.0, cell_size=1.0, cbar_label='Correlation'):
    df = pd.read_csv(input_file)

    if geneset_order is None:
        geneset_order = [
            "Neutrophil_surface_protein",
            "Neutrophil_activation",
            "Neutrophil_degranulation",
            "Primary_granule",
            "Secondary_granule",
            "Tertiary_granule",
            "Neutrophil_extracellular_trap",
            "Neutrophil_chemotaxis",
            "Neutrophil_migration"
        ]

    # 读取对应的tissue文件用于筛选
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    tissue_file = os.path.join(output_dir, f"{base_name}_tissues.txt")

    if os.path.exists(tissue_file):
        with open(tissue_file, 'r') as f:
            tissue_filter = set(line.strip() for line in f if line.strip())
        df = df[df['tissue_category'].isin(tissue_filter)]

    pivot_corr = df.pivot(index='geneset', columns='tissue_category', values='correlation')
    row_order = [g for g in geneset_order if g in pivot_corr.index]
    pivot_corr = pivot_corr.reindex(row_order)
    tissue_mean = pivot_corr.mean(axis=0).sort_values(ascending=False)
    column_order = list(tissue_mean.index)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    scale = max(float(cell_size), 0.1)
    return save_long_heatmap(
        df,
        output_dir / f"{base_name}_heatmap.svg",
        row='geneset',
        column='tissue_category',
        value='correlation',
        pvalue=pvalue_col,
        row_order=row_order,
        column_order=column_order,
        significance_mode='star',
        cbar_label=cbar_label,
        xlabel='Tissue',
        ylabel='Gene set',
        width=None,
        cell_width=0.36 * scale,
        cell_height=0.26 * scale,
        xtick_rotation=90,
        also=(".pdf", ".png"),
    )


def export_tissues_to_txt(input_file, output_dir):
    """Extract tissues from single CSV file and save to txt"""
    df = pd.read_csv(input_file)
    tissues = sorted(df['tissue_category'].unique())

    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}_tissues.txt")

    with open(output_file, 'w') as f:
        for tissue in tissues:
            f.write(tissue + '\n')

    print(f"导出 {len(tissues)} 个器官到: {output_file}")


if __name__ == '__main__':
    input_files = [
        "D:/Projects/Neutrophil_Aging/neutrophil_geneset_score/human_age_correlation_by_tissue.csv",
        "D:/Projects/Neutrophil_Aging/neutrophil_geneset_score/monkey_age_correlation_by_tissue.csv",
        "D:/Projects/Neutrophil_Aging/neutrophil_geneset_score/mouse_age_correlation_by_tissue.csv",
    ]
    output_dir = "D:/Projects/Neutrophil_Aging/neutrophil_geneset_score/figure"

    os.makedirs(output_dir, exist_ok=True)

    for input_file in input_files:
        # export_tissues_to_txt(input_file, output_dir)
        plot_correlation_heatmap(input_file, output_dir)

