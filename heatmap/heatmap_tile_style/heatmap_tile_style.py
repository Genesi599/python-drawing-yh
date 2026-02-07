#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : heatmap_tile_style.py
@Date    : 2026/2/3 15:40
@Author  : yh109
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib as mpl
import os




def plot_correlation_heatmap(input_file, output_dir, geneset_order=None, pvalue_col='pvalue',
                            font_scale=3.0, cell_size=1.0, cbar_label='Correlation'):
    df = pd.read_csv(input_file)

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['pdf.use14corefonts'] = False
    mpl.rcParams['font.family'] = 'Arial'

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

    df['sig_label'] = ''
    df.loc[df[pvalue_col] < 0.001, 'sig_label'] = '***'
    df.loc[(df[pvalue_col] >= 0.001) & (df[pvalue_col] < 0.01), 'sig_label'] = '**'
    df.loc[(df[pvalue_col] >= 0.01) & (df[pvalue_col] < 0.05), 'sig_label'] = '*'

    pivot_corr = df.pivot(index='geneset', columns='tissue_category', values='correlation')
    pivot_sig = df.pivot(index='geneset', columns='tissue_category', values='sig_label')

    pivot_corr = pivot_corr.reindex([g for g in geneset_order if g in pivot_corr.index])

    tissue_mean = pivot_corr.mean(axis=0).sort_values(ascending=False)
    pivot_corr = pivot_corr[tissue_mean.index]
    pivot_sig = pivot_sig[tissue_mean.index]

    vmax = max(abs(pivot_corr.min().min()), abs(pivot_corr.max().max()))
    vmin = -vmax

    n_rows, n_cols = pivot_corr.shape
    fig_width = n_cols * cell_size + 2
    fig_height = n_rows * cell_size + 1

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(pivot_corr, cmap='RdBu_r', center=0, vmin=vmin, vmax=vmax,
                cbar_kws={'label': cbar_label}, linewidths=0.5, linecolor='white',
                ax=ax, square=True)

    for i, geneset in enumerate(pivot_corr.index):
        for j, tissue in enumerate(pivot_corr.columns):
            label = pivot_sig.loc[geneset, tissue]
            if pd.notna(label) and label:
                ax.text(j + 0.5, i + 0.5, label,
                        ha='center', va='center',
                        fontsize=30, fontweight='bold', color='black',
                        verticalalignment='center', horizontalalignment='center')

    plt.xticks(rotation=90, ha='center', fontsize=12 * font_scale)
    plt.yticks(rotation=0, fontsize=12 * font_scale)
    ax.set_xlabel('')
    ax.set_ylabel('')

    cbar = ax.collections[0].colorbar
    ticks = [vmin, vmin / 2, 0, vmax / 2, vmax]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{t:.1f}' for t in ticks])
    cbar.ax.tick_params(labelsize=10 * font_scale)
    cbar.set_label(cbar_label, size=12 * font_scale)

    pos = ax.get_position()
    cbar_pos = cbar.ax.get_position()
    gap_inch = 0.3
    gap_relative = gap_inch / fig_width
    new_x0 = pos.x1 + gap_relative
    cbar.ax.set_position([new_x0, pos.y0, cbar_pos.width, pos.height])

    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_pdf = os.path.join(output_dir, f"{base_name}_heatmap.pdf")
    output_png = os.path.join(output_dir, f"{base_name}_heatmap.png")

    os.makedirs(output_dir, exist_ok=True)

    plt.savefig(output_pdf, dpi=300, pad_inches=0.1, bbox_inches='tight',
                metadata={'Creator': None, 'Producer': None})
    plt.savefig(output_png, dpi=300, pad_inches=0.1, bbox_inches='tight',)
    plt.close()


if __name__ == '__main__':
    input_file = "D:/Projects/Neutrophil_Aging/neutrophil_geneset_score/human_age_correlation_by_tissue.csv"
    output_dir = "D:/Projects/Neutrophil_Aging/neutrophil_geneset_score/figure"
    plot_correlation_heatmap(input_file, output_dir)
    input_file = "D:/Projects/Neutrophil_Aging/neutrophil_geneset_score/monkey_age_correlation_by_tissue.csv"
    plot_correlation_heatmap(input_file, output_dir)
    input_file = "D:/Projects/Neutrophil_Aging/neutrophil_geneset_score/mouse_age_correlation_by_tissue.csv"
    plot_correlation_heatmap(input_file, output_dir)


