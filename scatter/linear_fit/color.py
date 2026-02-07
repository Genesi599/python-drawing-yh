#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : color.py
@Date    : 2026/2/4 14:52
@Author  : yh109
"""

import pandas as pd
from matplotlib.colors import to_rgb, to_hex
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

SYSTEM_BASE_COLORS = {
    'Urinary system': '#8B4513',
    'Digestive system': '#90EE90',
    'Muscular system': '#FFB6C1',
    'Nervous system': '#4169E1',
    'Respiratory system': '#E6E6FA',
    'Adipose tissue': '#FFD700',
    'Endocrine system': '#228B22',
    'Integumentary system': '#FFA500',
    'Immune system': '#2F4F4F',
    'Reproductive system': '#87CEEB',
    'Cardiovascular system': '#DC143C'
}


def generate_tissue_colors(csv_path):
    df = pd.read_csv(csv_path)
    tissue_colors = {}

    for system, base_color in SYSTEM_BASE_COLORS.items():
        tissues = df[df['System'] == system]['Tissue'].tolist()
        n = len(tissues)

        if n == 0:
            continue

        if n == 1:
            tissue_colors[tissues[0]] = base_color
        else:
            base_rgb = to_rgb(base_color)
            for i, tissue in enumerate(tissues):
                factor = 0.6 + (i / (n - 1)) * 0.4
                new_rgb = tuple(max(0, min(1, c * factor)) for c in base_rgb)
                tissue_colors[tissue] = to_hex(new_rgb)

    df['Color'] = df['Tissue'].map(tissue_colors)
    df.to_csv(csv_path, index=False)

    return tissue_colors


def save_system_legend(system_colors, save_dir):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis('off')

    patches = [mpatches.Patch(color=color, label=system)
               for system, color in system_colors.items()]

    ax.legend(handles=patches, loc='center', fontsize=12, frameon=False)

    save_path = Path(save_dir)
    fig.savefig(save_path / 'system_colors_legend.png', dpi=300, bbox_inches='tight')
    fig.savefig(save_path / 'system_colors_legend.pdf', bbox_inches='tight')
    plt.close()


csv_path = r"D:\Projects\B_Cell_Aging\B_cell_ratio\ref\tissue_system_mapping.csv"
tissue_colors = generate_tissue_colors(csv_path)

save_dir = Path(csv_path).parent
save_system_legend(SYSTEM_BASE_COLORS, save_dir)