import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from drawing_yh.pie_chart.pie_chart import plot_pie

values = [251, 173, 65, 33, 30, 12, 1]
labels = ['Secreted', 'Surface film', 'ECM', 'Extracellular space',
          'Cell surface', 'Exosome', 'Other']

out_dir = os.path.dirname(__file__)

plot_pie(
    values=values,
    labels=labels,
    title='Secreted Proteins (Donut)',
    donut=True,
    donut_ratio=0.50,
    show_pct=True,
    label_font_size=10,
    pct_font_size=8,
    out_path=os.path.join(out_dir, 'donut.png'),
)
