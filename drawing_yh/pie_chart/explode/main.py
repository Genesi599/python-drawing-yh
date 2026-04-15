import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from drawing_yh.pie_chart.pie_chart import plot_pie

values = [251, 173, 65, 33, 30, 12, 1]
labels = ['Secreted', 'Surface film', 'ECM', 'Extracellular space',
          'Cell surface', 'Exosome', 'Other']

explode = [0.08, 0, 0, 0, 0, 0, 0]   # 突出第一块

out_dir = os.path.dirname(__file__)

plot_pie(
    values=values,
    labels=labels,
    title='Secreted Proteins\n(Secreted highlighted)',
    explode=explode,
    show_pct=True,
    label_font_size=10,
    pct_font_size=8,
    out_path=os.path.join(out_dir, 'explode.png'),
)
