import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from drawing_yh.pie_chart.zoom_pie import plot_zoom_pie

out_dir = os.path.dirname(__file__)

# ── 左侧：大类 ────────────────────────────────────────────────────
major_values = [565, 756, 101, 167, 80, 884, 315, 291]
major_labels = ['Secreted', 'Cell_membrane', 'Cell_periphery',
                'Endomembrane_system', 'Mitochondrion',
                'Cytoplasm', 'Nucleus', 'Unknown']
major_colors = ['#4C72B0', '#DD8452', '#55A868', '#C44E52',
                '#8172B3', '#937860', '#64B5CD']

# ── 右侧：Secreted 的细分（zoom_index=0） ─────────────────────────
detail_values = [251, 173, 65, 33, 30, 12, 1]
detail_labels = ['Secreted', 'Surface film', 'ECM',
                 'Extracellular space', 'Cell surface', 'Exosome', 'Other']

plot_zoom_pie(
    major_values=major_values,
    major_labels=major_labels,
    detail_values=detail_values,
    detail_labels=detail_labels,
    zoom_index=0,
    major_colors=major_colors,
    title_left='Proteome Subcellular Localization',
    title_right='Secreted Proteins',
    show_pct=True,
    show_count=True,
    font_size=8,
    merge_threshold=3.0,
    out_path=os.path.join(out_dir, 'zoom_pie.png'),
)

print('Done.')
