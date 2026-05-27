"""
drawing_yh.network — 网络图

模块化相关性网络(Bokeh):Louvain 切社区 + 两阶段力导向布局 + 密度 halo +
边带 + 高亮环 + 标签去重叠 → HTML / PNG / SVG。

入口
----
>>> import networkx as nx
>>> from drawing_yh.network import plot_modular_network
>>> G = nx.Graph()
>>> # ... 给每个节点设 G.nodes[n]['ntype'] = 'Protein' / 'Metabolite' / ... ──
>>> # ... 每条边设 G.edges[u,v]['r'] = 相关系数 ──
>>> G_drawn, module_id = plot_modular_network(
...     G,
...     node_type_colors={'Protein': '#7aa3d4', 'Metabolite': '#f0b07a'},
...     node_label_map={n: lbl for n, lbl in ...},
...     node_highlight={n: +1 for n in up_set} | {n: -1 for n in down_set},
...     outpath='out/net.png',
... )
"""
from .modular import plot_modular_network
from .layout import layout_meta_then_intra, merge_small_modules
from .halo import compute_density_background, compute_module_halos
from .labels import (
    deoverlap_labels, estimate_label_size, resolve_node_overlaps, fs_px,
)
from .edges import cubic_bezier_samples, quadratic_bezier_samples, bundled_arc
from .palette import MODULE_PALETTE, module_palette, lighten_hex, darken_hex
from .hub_spoke import hub_spoke

__all__ = [
    "plot_modular_network",
    "layout_meta_then_intra",
    "merge_small_modules",
    "compute_density_background",
    "compute_module_halos",
    "deoverlap_labels",
    "estimate_label_size",
    "resolve_node_overlaps",
    "fs_px",
    "cubic_bezier_samples",
    "quadratic_bezier_samples",
    "bundled_arc",
    "MODULE_PALETTE",
    "module_palette",
    "lighten_hex",
    "darken_hex",
    "hub_spoke",
]
