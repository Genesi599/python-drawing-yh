from math import log10

import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math
from adjustText import adjust_text


def create_volcano_plot(
        input_file, output_file='Volcano_plot.png',
        x_threshold=0.5, y_threshold=-np.log10(0.05),
        lfc_col='log2FoldChange', p_col='padj',
        id_col='GeneName', break_x=False):
    df = pd.read_csv(input_file)
    vol = df[[id_col, lfc_col, p_col]].copy()
    vol['y'] = -np.log10(vol[p_col])

    vol['group'] = 'black'
    up_mask = (vol[lfc_col] >= x_threshold) & (vol['y'] >= y_threshold)
    dn_mask = (vol[lfc_col] <= -x_threshold) & (vol['y'] >= y_threshold)
    vol.loc[up_mask, 'group'] = 'tab:red'
    vol.loc[dn_mask, 'group'] = 'tab:blue'
    vol['dist'] = vol[lfc_col] ** 2 + vol['y'] ** 2
    vol = vol.sort_values('dist', ascending=False)

    if break_x:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6),
                                       gridspec_kw={'wspace': 0.15})

        dn_vol = vol[vol[lfc_col] <= 0]
        up_vol = vol[vol[lfc_col] >= 0]

        x_min = vol[lfc_col].min()
        x_max = vol[lfc_col].max()
        x_margin = (x_max - x_min) * 0.05

        # 左图：负值范围
        dn_x_min = dn_vol[lfc_col].min()
        dn_x_max = dn_vol[lfc_col].max()
        dn_margin = (dn_x_max - dn_x_min) * 0.05
        ax1.set_xlim(dn_x_min - dn_margin, dn_x_max + dn_margin)

        # 右图：正值范围
        up_x_min = up_vol[lfc_col].min()
        up_x_max = up_vol[lfc_col].max()
        up_margin = (up_x_max - up_x_min) * 0.05
        ax2.set_xlim(up_x_min - up_margin, up_x_max + up_margin)

        ax1.scatter(dn_vol[lfc_col], dn_vol['y'], s=3, c=dn_vol['group'], alpha=0.8)
        ax2.scatter(up_vol[lfc_col], up_vol['y'], s=3, c=up_vol['group'], alpha=0.8)

        y_min, y_max = vol['y'].min(), vol['y'].max()
        y_margin = (y_max - y_min) * 0.05
        ax1.set_ylim(y_min - y_margin, y_max + y_margin)
        ax2.set_ylim(y_min - y_margin, y_max + y_margin)

        ax1.axhline(y_threshold, ls='--', color='grey', lw=1)
        ax1.set_ylabel('-Log10(adj. p-value)', fontweight='bold', fontsize=12)
        ax1.axvline(-x_threshold, ls='--', color='grey', lw=1)
        ax1.set_xlabel('Log2 Fold Change', fontweight='bold', fontsize=12)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        ax2.axhline(y_threshold, ls='--', color='grey', lw=1)
        ax2.axvline(x_threshold, ls='--', color='grey', lw=1)
        ax2.set_xlabel('Log2 Fold Change', fontweight='bold', fontsize=12)
        ax2.spines['top'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.yaxis.set_visible(False)

        d = 0.03
        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=18,
                      linestyle='none', color='k', mec='k', mew=1.5, clip_on=False)
        ax1.plot([1, 1], [1, 0], transform=ax1.transAxes, **kwargs)
        ax2.plot([0, 0], [0, 1], transform=ax2.transAxes, **kwargs)

        dn_sig = vol[dn_mask].head(10)
        texts = [ax1.text(r[lfc_col] + 0.02, r['y'], r[id_col],
                          fontsize=9, style='italic', weight='bold', ha='right')
                 for _, r in dn_sig.iterrows()]
        if texts:
            adjust_text(texts, ax=ax1, arrowprops=dict(arrowstyle='->', lw=0.5))

        up_sig = vol[up_mask].head(10)
        texts = [ax2.text(r[lfc_col] - 0.02, r['y'], r[id_col],
                          fontsize=9, style='italic', weight='bold', ha='left')
                 for _, r in up_sig.iterrows()]
        if texts:
            adjust_text(texts, ax=ax2, arrowprops=dict(arrowstyle='->', lw=0.5))
    else:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.scatter(vol[lfc_col], vol['y'], s=3, c=vol['group'], alpha=0.8)
        ax.axvline(-x_threshold, ls='--', color='grey', lw=1)
        ax.axvline(x_threshold, ls='--', color='grey', lw=1)
        ax.axhline(y_threshold, ls='--', color='grey', lw=1)
        ax.set_xlabel('Log2 Fold Change', fontweight='bold', fontsize=12)
        ax.set_ylabel('-Log10(adj. p-value)', fontweight='bold', fontsize=12)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        dn_sig = vol[dn_mask].head(10)
        texts = [ax.text(r[lfc_col] - 0.05, r['y'], r[id_col],
                         fontsize=10, style='italic', weight='bold')
                 for _, r in dn_sig.iterrows()]
        if texts:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', lw=0.5))

        up_sig = vol[up_mask].head(10)
        texts = [ax.text(r[lfc_col] - 0.05, r['y'], r[id_col],
                         fontsize=10, style='italic', weight='bold')
                 for _, r in up_sig.iterrows()]
        if texts:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', lw=0.5))

    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    fig.savefig(output_file.replace('.png', '.pdf'), dpi=300, bbox_inches='tight')
    plt.close()

# # 示例使用
create_volcano_plot(
    input_file=r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\data\by_sex\protein_age_correlations_significant_all.csv",
    output_file='Volcano_plot.png',
    # x_threshold=0.3,
    y_threshold=-np.log10(0.05),
    lfc_col='pearson_r',
    p_col='pearson_pval',
    id_col='gene_primary',
    break_x=True
)