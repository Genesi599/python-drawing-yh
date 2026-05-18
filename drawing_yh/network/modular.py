"""
模块化相关性网络 (Bokeh) — 通用绘图入口。

接收一个已经搭好的 networkx.Graph,每个节点须带 `ntype` 属性(或自定义),
按 Louvain 切社区 → 两阶段力导向布局 → 密度 halo 背景 → 边带 → 节点 +
高亮环 + 标签去重叠 → HTML / PNG / SVG 三件套。

设计原则:
- **数据无关**:不读 corr_*/, 不识别 Protein/Metabolite —— 这些是项目的活,
  调用者负责把 G 搭好(`G.nodes[n]["ntype"]` 给类型,attr 名可配)
- **高亮可选**:age-sig / DE 节点统一走 `node_highlight={node_id: ±1}` 接口,
  +1 显示上调环色,-1 下调环色,缺失则不画环
- **字号统一**:默认 "8px"(与 matplotlib panel 视觉对齐),改 "8pt" 走标准印刷字号
"""
from __future__ import annotations

import math
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import networkx as nx

from bokeh.plotting import figure, save, output_file
from bokeh.models import ColumnDataSource, LabelSet, HoverTool, Range1d, Title
from bokeh.core.properties import value as Value
from bokeh.io import export_png, export_svg

from .palette import MODULE_PALETTE, module_palette, lighten_hex, darken_hex
from .layout import layout_meta_then_intra, merge_small_modules
from .edges import cubic_bezier_samples
from .labels import fs_px, deoverlap_labels, estimate_label_size, resolve_node_overlaps
from .halo import compute_density_background


def _detect_communities(
    G: nx.Graph, seed: int = 1, resolution: float = 1.0,
) -> dict:
    """Louvain 社区检测,返回 node→module_id 映射。孤立节点归为最后一个 orphan 模块。"""
    subnodes = [n for n in G.nodes if G.degree[n] > 0]
    G_sub = G.subgraph(subnodes)

    try:
        communities = nx.community.louvain_communities(
            G_sub, seed=seed, resolution=resolution
        )
    except AttributeError:
        communities = list(nx.community.greedy_modularity_communities(G_sub))

    communities = sorted(communities, key=len, reverse=True)

    module_id = {}
    for i, comm in enumerate(communities):
        for n in comm:
            module_id[n] = i

    orphan_mod = len(communities)
    for n in G.nodes:
        if n not in module_id:
            module_id[n] = orphan_mod
    return module_id


def plot_modular_network(
    G: nx.Graph,
    *,
    # ── 节点类型 / 可视化 ──
    node_type_attr: str = "ntype",
    node_type_colors: dict | None = None,
    node_label_map: dict | None = None,
    node_highlight: dict | None = None,
    highlight_colors: tuple = ("#dc2626", "#1e40af"),
    module_function_map: dict | None = None,

    # ── 社区检测(module_id 为 None 时启用) ──
    module_id: dict | None = None,
    community_seed: int = 1,
    community_resolution: float = 1.0,
    min_module_size: int = 0,

    # ── 节点过滤 ──
    keep_isolated: bool = False,
    min_degree: int = 2,
    max_draw_edges: int = 12000,
    edge_intra_ratio: float = 0.75,

    # ── 布局 ──
    layout_canvas_size: float = 1800.0,
    layout_module_spacing: float = 1.2,
    layout_iterations_meta: int = 200,
    layout_iterations_intra: int = 120,
    layout_intra_radius_base: float = 40.0,
    layout_intra_radius_scale: float = 18.0,
    layout_radial_power: float = 0.55,
    layout_inter_padding: float = 30.0,
    layout_meta_layout: str = "spring",

    # ── halo / 节点 / 边 ──
    halo_buffer_ratio: float = 0.35,
    halo_min_radius: float = 40.0,
    halo_alpha: float = 0.22,
    halo_line_alpha: float = 0.4,
    node_size: float = 5.0,
    hub_size_max: float = 26.0,
    hub_degree_percentile: float = 0.985,
    edge_alpha: float = 0.18,
    edge_intra_mult: float = 1.0,
    edge_inter_mult: float = 0.35,
    edge_width_scale: float = 0.9,
    label_top_n: int = 60,
    label_module_hub_n: int = 0,
    hub_ring_color: str = "#374151",

    # ── 字号 / title ──
    font_family: str = "Arial",
    font_size: str = "8px",
    text_color: str = "#0f172a",
    show_title: bool = False,
    title: str = "Modular network",
    subtitle: str | None = None,

    # ── 输出 ──
    fig_width: int = 850,
    fig_height: int = 850,
    outpath: str | None = None,
    save_module_csv: bool = True,
) -> tuple[nx.Graph, dict]:
    """绘制模块化相关性网络。

    Parameters
    ----------
    G : networkx.Graph
        已构建的图;每个节点须设置 `node_type_attr`(默认 "ntype")属性,
        每条边可带 `r` (相关系数) / `p` / `adjp` 属性(用于边色与边采样)。
    module_id : dict | None
        {node_id: int} 模块归属;None 时本函数内部跑 Louvain。
    node_type_colors : dict
        {ntype_value: hex_color},如 {"Protein": "#7aa3d4", "Metabolite": "#f0b07a"}。
    node_label_map : dict | None
        {str(node_id): display_label};缺失节点 fallback 显示 id。
    node_highlight : dict | None
        {str(node_id): ±1};+1 = 上调高亮环色,-1 = 下调,缺失 = 不画环。
    highlight_colors : (up_color, down_color)
    module_function_map : dict | None
        {module_id: str or list[str]} 模块功能短标签,显示在模块中心 ID 下方。
    label_top_n : int
        全局按 degree 取前 N 个节点标 label(0 = 关掉)。
    label_module_hub_n : int
        **每模块**按 degree 取前 N 个节点标 label(0 = 关掉);同时给这些 hub 节点
        画一圈灰色环(`hub_ring_color`)。与 node_highlight 重叠时高亮色优先。
    hub_ring_color : str
        module hub 节点的环色(默认 `#374151` slate-700)。

    其它参数同原 net_Bokeh_modular.plot_modular_network。

    Returns
    -------
    (G, module_id) : 过滤 + 模块合并后的图与最终模块映射。
    """
    node_type_colors = node_type_colors or {}
    node_label_map = node_label_map or {}
    module_function_map = module_function_map or {}

    # 高亮 dict 归一化:既支持 {id: ±1} 也支持 set(id)
    def _norm_dir(x):
        if x is None:
            return {}
        if isinstance(x, dict):
            return {str(k): (1 if int(v) >= 0 else -1) for k, v in x.items()}
        return {str(k): 1 for k in x}
    node_highlight = _norm_dir(node_highlight)

    # ── 4. 按最低度数过滤节点 ──
    if keep_isolated:
        thr_ = 0
    else:
        thr_ = max(int(min_degree), 1)
    if thr_ > 0:
        low = [n for n in G.nodes if G.degree[n] < thr_]
        if low:
            G = G.copy()
            G.remove_nodes_from(low)
        print(f"  去掉 {len(low)} 个 degree<{thr_} 的节点,剩: {G.number_of_nodes()}")

    # ── 5. Louvain ──
    if module_id is None:
        print(f"▶ Louvain 社区检测 (resolution={community_resolution})")
        module_id = _detect_communities(G, seed=community_seed, resolution=community_resolution)
    raw_mod_count = len(set(module_id.values()))
    if min_module_size and min_module_size > 1:
        module_id = merge_small_modules(G, module_id, min_size=min_module_size)
        merged_count = len(set(module_id.values()))
        print(f"  小模块合并 (min_size={min_module_size}): {raw_mod_count} → {merged_count}")
    n_modules_total = len(set(module_id.values()))
    sizes = Counter(module_id.values())
    print(f"  模块数: {n_modules_total}")
    for m, s in sizes.most_common(15):
        print(f"    模块 {m}: {s} 节点")

    # ── 6-7. 两阶段布局 ──
    print("▶ 计算布局 (meta-graph + per-module spring_layout)")
    pos, module_centers = layout_meta_then_intra(
        G, module_id,
        seed=community_seed,
        iterations_meta=layout_iterations_meta,
        iterations_intra=layout_iterations_intra,
        canvas_size=layout_canvas_size,
        module_spacing=layout_module_spacing,
        intra_radius_base=layout_intra_radius_base,
        intra_radius_scale=layout_intra_radius_scale,
        radial_power=layout_radial_power,
        inter_padding=layout_inter_padding,
        meta_layout=layout_meta_layout,
    )

    # ── 8. 色板 ──
    palette = module_palette(n_modules_total, sizes)

    # ── 8.5 把同一模块内的 highlight (标注) 节点推开,让标签更均匀分布 ──
    if node_highlight:
        labeled_by_mod = defaultdict(list)
        for n in pos:
            if str(n) in node_highlight:
                labeled_by_mod[module_id[n]].append(n)

        for m, lab_nodes in labeled_by_mod.items():
            if len(lab_nodes) < 2:
                continue
            cx, cy = module_centers[m]

            _all_in_mod = [n for n in pos if module_id[n] == m]
            _mod_r_list = [math.hypot(pos[n][0] - cx, pos[n][1] - cy) for n in _all_in_mod]
            if not _mod_r_list:
                continue
            r_actual = float(np.quantile(_mod_r_list, 0.85))
            r_max_in_mod = r_actual

            target_pair_d = max(2 * r_actual / math.sqrt(len(lab_nodes) + 1) * 0.95, 3.0)

            for _ in range(30):
                for i in range(len(lab_nodes)):
                    ni = lab_nodes[i]
                    xi, yi = pos[ni]
                    fx, fy = 0.0, 0.0
                    for j in range(len(lab_nodes)):
                        if i == j:
                            continue
                        nj = lab_nodes[j]
                        xj, yj = pos[nj]
                        dx, dy = xi - xj, yi - yj
                        d = math.sqrt(dx * dx + dy * dy) + 1e-9
                        if d < target_pair_d:
                            mag = (target_pair_d - d) * 0.10
                            fx += dx / d * mag
                            fy += dy / d * mag
                    nx_ = xi + fx
                    ny_ = yi + fy
                    rel_dx = nx_ - cx
                    rel_dy = ny_ - cy
                    rel_r = math.sqrt(rel_dx * rel_dx + rel_dy * rel_dy)
                    if rel_r > r_max_in_mod:
                        nx_ = cx + rel_dx / rel_r * r_max_in_mod
                        ny_ = cy + rel_dy / rel_r * r_max_in_mod
                    pos[ni] = (nx_, ny_)
        print(f"  ► 高亮节点模块内均匀分散完成")

    # ── 11. Hub 阈值 ──
    deg_map = dict(G.degree())
    deg_vals = np.array(list(deg_map.values())) if deg_map else np.array([0])
    hub_thr = max(int(np.quantile(deg_vals, hub_degree_percentile)), 2)
    dmax = max(deg_vals.max(), 1)
    print(f"  Hub 阈值: degree ≥ {hub_thr}")

    if node_highlight:
        sig_in_g = sum(1 for n in G.nodes if str(n) in node_highlight)
        print(f"  高亮节点 (在网络中): {sig_in_g}/{len(node_highlight)}")

    # ── 12. 节点列 ──
    top_labeled = set(
        [n for n, _ in sorted(deg_map.items(), key=lambda kv: -kv[1])[:int(label_top_n)]]
    )
    # 每模块按 global degree 取前 N 个作为 module hub,既标 label 也画 ring
    module_hubs: set = set()
    if label_module_hub_n > 0:
        _by_mod = defaultdict(list)
        for _n, _m in module_id.items():
            if _n in G:
                _by_mod[_m].append(_n)
        for _m, _nodes in _by_mod.items():
            _sorted = sorted(_nodes, key=lambda x: -deg_map.get(x, 0))
            for _n in _sorted[:int(label_module_hub_n)]:
                module_hubs.add(_n)
                top_labeled.add(_n)

    node_ids = list(G.nodes)
    n_x, n_y = [], []
    n_fill, n_size, n_line_c, n_line_w = [], [], [], []
    n_type, n_module, n_degree, n_label, n_full_label = [], [], [], [], []
    n_marker, n_glow_size, n_glow_color = [], [], []
    n_high = []  # "Up"/"Down"/"No"
    n_ring_size, n_ring_color = [], []
    kept_ids = []

    UP_COLOR, DOWN_COLOR = highlight_colors

    # degree → size 连续映射(sqrt 压缩,避免少数超大 hub 吃掉所有尺寸预算)
    d_max_safe = float(max(dmax, 1))

    # 先算每个节点的像素尺寸,换算到 canvas 半径,做节点级碰撞解开
    node_pixel_size = {}
    for n in node_ids:
        if n not in pos:
            continue
        d_ = deg_map[n]
        sf_ = math.sqrt(d_ / d_max_safe) if d_ > 0 else 0.0
        node_pixel_size[n] = node_size + (hub_size_max - node_size) * sf_

    _xs_pre = np.array([pos[n][0] for n in node_ids if n in pos])
    _ys_pre = np.array([pos[n][1] for n in node_ids if n in pos])
    _pad = 180.0
    _x_span = (_xs_pre.max() - _xs_pre.min()) + 2 * _pad
    _y_span = (_ys_pre.max() - _ys_pre.min()) + 2 * _pad
    _canvas_per_px = max(_x_span / fig_width, _y_span / fig_height)
    node_radius_canvas = {
        n: (s / 2.0) * _canvas_per_px for n, s in node_pixel_size.items()
    }
    print(f"▶ 解开节点重叠 (canvas/px ≈ {_canvas_per_px:.2f})")
    pos = resolve_node_overlaps(
        pos, node_radius_canvas,
        iterations=80,
        padding=_canvas_per_px * 0.05,
        shift_frac=0.35,
        min_gap_factor=0.05,
    )

    # ── 计算每个模块标签的 bbox,用于把高亮节点推出 bbox ──
    _xs_real = np.array([pos[n][0] for n in pos])
    _ys_real = np.array([pos[n][1] for n in pos])
    _real_pad = 2.0
    _data_w_real = (_xs_real.max() - _xs_real.min()) + 2 * _real_pad
    _data_h_real = (_ys_real.max() - _ys_real.min()) + 2 * _real_pad
    _target_area = float(fig_width) * float(fig_height)
    _aspect_real = _data_w_real / max(_data_h_real, 1e-6)
    _fw_real = math.sqrt(_target_area * _aspect_real)
    _fh_real = math.sqrt(_target_area / _aspect_real)
    _canvas_per_px_real = max(_data_w_real / _fw_real, _data_h_real / _fh_real)

    _id_px_bbox = fs_px(font_size)
    _fn_px_bbox = fs_px(font_size)
    _line_h_canvas = max(_id_px_bbox, _fn_px_bbox) * 1.15 * _canvas_per_px_real
    _big_thr_bbox = max(20, int(sizes.most_common(1)[0][1] * 0.08))
    _small_thr_bbox = 50
    mod_lbl_rects = []
    for _m, (mcx_, mcy_) in module_centers.items():
        _func = module_function_map.get(_m) or module_function_map.get(int(_m))
        _big_enough = sizes[_m] >= _big_thr_bbox
        if not (_big_enough or _func):
            continue
        if isinstance(_func, str):
            _fn_list = [_func]
        elif isinstance(_func, (list, tuple)):
            _fn_list = [str(x) for x in _func]
        else:
            _fn_list = []
        if sizes[_m] < _small_thr_bbox and len(_fn_list) > 1:
            _fn_list = _fn_list[:1]
        _id_text = f"M{_m}  ·  n={sizes[_m]}"
        _iw, _ih = estimate_label_size(_id_text, _id_px_bbox, _canvas_per_px_real)
        _n_fn = len(_fn_list)
        if _n_fn > 0:
            _max_fn_w = max(estimate_label_size(t, _fn_px_bbox, _canvas_per_px_real)[0] for t in _fn_list)
            _comp_w = max(_iw, _max_fn_w)
            _comp_h = _line_h_canvas * (1 + _n_fn)
        else:
            _comp_w = _iw
            _comp_h = _ih
        _pad_lbl = _canvas_per_px_real * 1.0
        _hw = _comp_w / 2 + _pad_lbl
        _hh = _comp_h / 2 + _pad_lbl
        mod_lbl_rects.append((_m, mcx_ - _hw, mcy_ - _hh, mcx_ + _hw, mcy_ + _hh))

    # ── 把高亮节点从 label bbox 内推出去(沿模块径向)──
    if mod_lbl_rects and node_highlight:
        _shifted = 0
        for n in list(pos):
            if str(n) not in node_highlight:
                continue
            x, y = pos[n]
            m = module_id[n]
            cx, cy = module_centers.get(m, (0.0, 0.0))
            for (mb, xlo, ylo, xhi, yhi) in mod_lbl_rects:
                if not (xlo <= x <= xhi and ylo <= y <= yhi):
                    continue
                dx, dy = x - cx, y - cy
                rr = math.hypot(dx, dy)
                if rr < 0.5:
                    ang = (hash(str(n)) & 0xFF) * 2 * math.pi / 256
                    dx, dy = math.cos(ang), math.sin(ang)
                    rr = 1.0
                ux, uy = dx / rr, dy / rr
                _all_in_mod = [pos[k] for k in pos if module_id[k] == m]
                if _all_in_mod:
                    _dists = [math.hypot(p[0] - cx, p[1] - cy) for p in _all_in_mod]
                    r_max = float(np.quantile(_dists, 0.92))
                else:
                    r_max = 1e9
                step = _canvas_per_px_real * 2.0
                nx_, ny_ = x, y
                for _ in range(20):
                    nx_ += ux * step
                    ny_ += uy * step
                    if not (xlo <= nx_ <= xhi and ylo <= ny_ <= yhi):
                        break
                    if math.hypot(nx_ - cx, ny_ - cy) > r_max:
                        nx_ -= ux * step
                        ny_ -= uy * step
                        break
                pos[n] = (nx_, ny_)
                _shifted += 1
                break
        if _shifted:
            print(f"  ► 高亮节点避让模块标签:微调 {_shifted} 个")

    # 计算每个模块的"散点阈值半径":超过此距离的非高亮节点裁掉
    SCATTER_FACTOR = 0.85
    _mod_scatter_r = {}
    for _m in module_centers:
        _mod_dists = []
        for _nn in pos:
            if module_id.get(_nn) == _m:
                _mod_dists.append(math.hypot(pos[_nn][0] - module_centers[_m][0],
                                             pos[_nn][1] - module_centers[_m][1]))
        if _mod_dists:
            _r_dense = float(np.quantile(_mod_dists, 0.75))
            _mod_scatter_r[_m] = max(_r_dense * 1.5, 5.0)
        else:
            _mod_scatter_r[_m] = 1e9

    _skipped_scatter = 0
    for n in node_ids:
        if n not in pos:
            continue
        x, y = pos[n]
        t = G.nodes[n][node_type_attr]
        m = module_id[n]

        is_highlight = str(n) in node_highlight
        cmx, cmy = module_centers.get(m, (0.0, 0.0))
        _dist_to_mod = math.hypot(x - cmx, y - cmy)
        if _dist_to_mod > _mod_scatter_r.get(m, 1e9) and not is_highlight:
            _skipped_scatter += 1
            continue

        fill = node_type_colors.get(t, "#9ca3af")
        marker = "circle"

        d = deg_map[n]
        size_frac = math.sqrt(d / d_max_safe) if d > 0 else 0.0
        size = node_size + (hub_size_max - node_size) * size_frac

        age_dir = node_highlight.get(str(n), 0)

        # 所有节点用相同的淡描边,不再用 hub 黑边 / 外发光 区分
        line_c = darken_hex(fill, 0.35)
        line_w = 0.4
        glow_size = 0.0
        glow_color = fill

        # 环:高亮节点 = ±方向色;module hub = hub_ring_color;两者重叠以高亮色优先
        if age_dir != 0:
            ring_size = max(size * 1.95, size + 6.0)
            ring_color = UP_COLOR if age_dir > 0 else DOWN_COLOR
        elif n in module_hubs:
            ring_size = max(size * 1.95, size + 6.0)
            ring_color = hub_ring_color
        else:
            ring_size = 0.0
            ring_color = fill

        lab = node_label_map.get(str(n), str(n))

        kept_ids.append(str(n))
        n_x.append(x)
        n_y.append(y)
        n_fill.append(fill)
        n_size.append(size)
        n_line_c.append(line_c)
        n_line_w.append(line_w)
        n_type.append(t)
        n_module.append(m)
        n_degree.append(d)
        n_full_label.append(lab)
        # 高亮节点强制显示标签(覆盖 top_n 限制)
        if age_dir != 0 or n in top_labeled:
            n_label.append(lab)
        else:
            n_label.append("")
        n_marker.append(marker)
        n_glow_size.append(glow_size)
        n_glow_color.append(glow_color)
        n_high.append("Up" if age_dir > 0 else ("Down" if age_dir < 0 else "No"))
        n_ring_size.append(ring_size)
        n_ring_color.append(ring_color)

    if _skipped_scatter > 0:
        print(f"  ► 跳过 {_skipped_scatter} 个离群散点 (>{SCATTER_FACTOR}×模块半径,且非高亮)")

    node_source = ColumnDataSource(data=dict(
        x=n_x, y=n_y, id=kept_ids,
        fill=n_fill, size=n_size,
        line_color=n_line_c, line_width=n_line_w,
        ntype=n_type, module=n_module, degree=n_degree,
        label=n_label, full_label=n_full_label,
        marker=n_marker, glow_size=n_glow_size, glow_color=n_glow_color,
        highlight=n_high,
        ring_size=n_ring_size, ring_color=n_ring_color,
    ))

    # ── 13. 边列(分层采样) ──
    all_edges = list(G.edges(data=True))
    if max_draw_edges and len(all_edges) > max_draw_edges:
        intra_budget = int(max_draw_edges * edge_intra_ratio)
        inter_budget = max_draw_edges - intra_budget
        intra_edges = [e for e in all_edges if module_id[e[0]] == module_id[e[1]]]
        inter_edges = [e for e in all_edges if module_id[e[0]] != module_id[e[1]]]
        intra_edges.sort(key=lambda e: -abs(float(e[2].get("r", 0))))

        # 按源模块分桶
        inter_by_module = defaultdict(list)
        for e in inter_edges:
            inter_by_module[module_id[e[0]]].append(e)
            inter_by_module[module_id[e[1]]].append(e)
        for m in inter_by_module:
            inter_by_module[m].sort(key=lambda e: -abs(float(e[2].get("r", 0))))

        n_mods_with_inter = max(len(inter_by_module), 1)
        per_mod_quota = max(1, inter_budget // n_mods_with_inter)

        kept_keys = set()
        kept_inter = []
        for m, edges in inter_by_module.items():
            for e in edges[:per_mod_quota]:
                k = (min(e[0], e[1]), max(e[0], e[1]))
                if k not in kept_keys:
                    kept_keys.add(k)
                    kept_inter.append(e)

        if len(kept_inter) < inter_budget:
            inter_edges_sorted = sorted(
                inter_edges, key=lambda e: -abs(float(e[2].get("r", 0)))
            )
            for e in inter_edges_sorted:
                if len(kept_inter) >= inter_budget:
                    break
                k = (min(e[0], e[1]), max(e[0], e[1]))
                if k not in kept_keys:
                    kept_keys.add(k)
                    kept_inter.append(e)

        all_edges = intra_edges[:intra_budget] + kept_inter
        print(f"  绘图边: 模块内 {len(intra_edges[:intra_budget])} / "
              f"跨模块 {len(kept_inter)} (按源模块每模块配额 {per_mod_quota})")

        # 强制加入高亮节点的所有边
        if node_highlight:
            already = {(min(u, v), max(u, v)) for u, v, _ in all_edges}
            extra_edges = []
            for n in G.nodes:
                if str(n) not in node_highlight:
                    continue
                for nb in G.neighbors(n):
                    key = (min(n, nb), max(n, nb))
                    if key in already:
                        continue
                    already.add(key)
                    extra_edges.append((n, nb, G[n][nb]))
            all_edges.extend(extra_edges)
            print(f"  + 强制加入高亮节点边: {len(extra_edges)} → 总 {len(all_edges)}")

    edge_xs, edge_ys, edge_colors, edge_widths, edge_alphas = [], [], [], [], []
    # 正相关偏温暖(玫瑰粉),负相关偏冷色(紫蓝) — 更柔和
    POS_COLOR = "#e69ba5"
    NEG_COLOR = "#9eb1ea"

    # ── 跨模块边束:每对模块 (A,B) 在 A/B 边缘各设 trunk 中转点 ──
    pair_count = Counter()
    for uu, vv, _ in all_edges:
        a, b = module_id[uu], module_id[vv]
        if a != b:
            pair_count[(min(a, b), max(a, b))] += 1

    _mod_radius_actual: dict = {}
    for _m in module_centers:
        _cx, _cy = module_centers[_m]
        _rs = [
            math.hypot(pos[_n][0] - _cx, pos[_n][1] - _cy)
            for _n in pos
            if module_id.get(_n) == _m
        ]
        _mod_radius_actual[_m] = max(_rs) if _rs else 0.0

    pair_meta = {}
    trunk_inset = 1.12
    for key, cnt in pair_count.items():
        mA_, mB_ = key
        if mA_ not in module_centers or mB_ not in module_centers:
            continue
        cA = module_centers[mA_]; cB = module_centers[mB_]
        dx = cB[0] - cA[0]; dy = cB[1] - cA[1]
        L = math.hypot(dx, dy) or 1.0
        ux_, uy_ = dx / L, dy / L
        rA_eff = _mod_radius_actual.get(mA_, 0.0) * trunk_inset
        rB_eff = _mod_radius_actual.get(mB_, 0.0) * trunk_inset
        trunk_a = (cA[0] + ux_ * rA_eff, cA[1] + uy_ * rA_eff)
        trunk_b = (cB[0] - ux_ * rB_eff, cB[1] - uy_ * rB_eff)
        pair_meta[key] = dict(trunk_a=trunk_a, trunk_b=trunk_b, total=cnt)

    inter_width_scale = 0.45

    # 高亮节点 ID 集合,决定模块内边的可见度
    _high_set = set()
    for n in G.nodes:
        if str(n) in node_highlight:
            _high_set.add(n)

    # 模块内"无高亮端点"的边的衰减系数
    INTRA_BG_ALPHA_MULT = 0.20
    INTRA_BG_WIDTH_MULT = 0.50

    for u, v, d in all_edges:
        if u not in pos or v not in pos:
            continue
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        r_val = float(d.get("r", 0))
        w = (0.3 + 0.9 * abs(r_val)) * edge_width_scale
        is_cross = module_id[u] != module_id[v]
        if is_cross:
            mA, mB = module_id[u], module_id[v]
            key = (min(mA, mB), max(mA, mB))
            meta = pair_meta.get(key)
            if meta is not None:
                if module_id[u] == key[0]:
                    t_u, t_v = meta["trunk_a"], meta["trunk_b"]
                else:
                    t_u, t_v = meta["trunk_b"], meta["trunk_a"]
                xs_, ys_ = cubic_bezier_samples(
                    (x0, y0), t_u, t_v, (x1, y1), n=22,
                )
                edge_xs.append(xs_); edge_ys.append(ys_)
            else:
                edge_xs.append([x0, x1]); edge_ys.append([y0, y1])
            edge_colors.append(POS_COLOR if r_val >= 0 else NEG_COLOR)
            edge_alphas.append(edge_alpha * edge_inter_mult)
            edge_widths.append(w * inter_width_scale)
        else:
            edge_xs.append([x0, x1])
            edge_ys.append([y0, y1])
            m_col = palette[module_id[u]] if module_id[u] < len(palette) else "#94a3b8"
            edge_colors.append(darken_hex(m_col, 0.1))
            touches_high = (u in _high_set) or (v in _high_set)
            if touches_high:
                edge_alphas.append(edge_alpha * edge_intra_mult)
                edge_widths.append(w)
            else:
                edge_alphas.append(edge_alpha * edge_intra_mult * INTRA_BG_ALPHA_MULT)
                edge_widths.append(w * INTRA_BG_WIDTH_MULT)

    edge_source = ColumnDataSource(dict(
        xs=edge_xs, ys=edge_ys,
        color=edge_colors, width=edge_widths, alpha=edge_alphas,
    ))

    # ── 14. 绘图 ──
    xs_all = np.array(n_x)
    ys_all = np.array(n_y)
    pad = 2
    x_lo, x_hi = xs_all.min() - pad, xs_all.max() + pad
    y_lo, y_hi = ys_all.min() - pad, ys_all.max() + pad
    data_w = x_hi - x_lo
    data_h = y_hi - y_lo
    # 图幅自适应:保持总像素面积 ≈ fig_width × fig_height,但宽高比贴合数据
    target_area = float(fig_width) * float(fig_height)
    aspect = data_w / max(data_h, 1e-6)
    fig_w_eff = int(round(math.sqrt(target_area * aspect)))
    fig_h_eff = int(round(math.sqrt(target_area / aspect)))
    fig_w_eff = max(fig_w_eff, 800)
    fig_h_eff = max(fig_h_eff, 800)
    print(f"  ► 图幅自适应: 数据宽:高={aspect:.2f} → fig {fig_w_eff}×{fig_h_eff} px")

    p = figure(
        width=fig_w_eff, height=fig_h_eff,
        x_range=Range1d(x_lo, x_hi),
        y_range=Range1d(y_lo, y_hi),
        tools="pan,wheel_zoom,box_zoom,reset,save",
        toolbar_location=None,
        background_fill_color="white",
        border_fill_color="white",
        outline_line_color=None,
        match_aspect=True,
        min_border=0,
        min_border_top=0, min_border_bottom=0,
        min_border_left=0, min_border_right=0,
    )
    p.axis.visible = False
    p.grid.visible = False
    p.title.text_font = font_family
    p.title.text_font_size = font_size
    p.title.text_color = text_color

    _canvas_per_px = max(data_w / float(fig_w_eff), data_h / float(fig_h_eff))

    if show_title:
        p.add_layout(Title(
            text=title,
            text_font_size=font_size, text_font_style="normal",
            text_color=text_color, align="center",
            text_font=font_family,
        ), "above")
        if subtitle:
            p.add_layout(Title(
                text=subtitle,
                text_font_size=font_size, text_font_style="normal",
                text_color="#64748b", align="center",
                text_font=font_family,
            ), "above")

    # 密度渐变背景
    density_bg = compute_density_background(
        pos, module_id, palette, sizes,
        x_range=(x_lo, x_hi), y_range=(y_lo, y_hi),
        grid_size=520,
        sigma_pixels=7.0,
        max_alpha=0.06,
        alpha_gamma=1.2,
        density_cutoff=0.08,
        overlap_cap=0.6,
    )
    if density_bg is not None:
        p.image_rgba(
            image=[density_bg],
            x=x_lo, y=y_lo,
            dw=x_hi - x_lo, dh=y_hi - y_lo,
            level="image",
        )

    # 边
    p.multi_line(
        xs="xs", ys="ys",
        line_color="color", line_width="width", line_alpha="alpha",
        line_cap="round",
        source=edge_source,
    )

    # Hub 外发光层(普通节点 glow_size=0,只影响 hub)
    p.scatter(
        x="x", y="y", size="glow_size",
        marker="marker",
        fill_color="glow_color", fill_alpha=0.25,
        line_alpha=0.0,
        source=node_source,
    )

    # 高亮节点环
    p.scatter(
        x="x", y="y", size="ring_size",
        marker="circle",
        fill_color="ring_color", fill_alpha=0.18,
        line_color="ring_color", line_width=1.6, line_alpha=0.95,
        source=node_source,
    )

    # 节点主层
    node_renderer = p.scatter(
        x="x", y="y",
        size="size", marker="marker",
        fill_color="fill", fill_alpha=1.0,
        line_color="line_color", line_width="line_width",
        source=node_source,
    )

    hover = HoverTool(
        renderers=[node_renderer],
        tooltips=[
            ("ID", "@id"),
            ("Type", "@ntype"),
            ("Label", "@full_label"),
            ("Module", "@module"),
            ("Degree", "@degree"),
            ("Highlight", "@highlight"),
        ],
    )
    p.add_tools(hover)

    # ════════ 所有文字标签 → 统一收集 → 同一次 deoverlap 跑 ════════
    _node_px = fs_px(font_size)
    _id_px = fs_px(font_size)
    _fn_px = fs_px(font_size)
    _line_spacing_px = max(_id_px, _fn_px) * 1.15
    _line_spacing_canvas = _line_spacing_px * _canvas_per_px
    _half_offset = _line_spacing_canvas / 2.0

    # 1) 节点标签
    _node_anchors, _node_widths, _node_heights, _node_texts = [], [], [], []
    for _i in range(len(node_source.data["x"])):
        _txt = node_source.data["label"][_i]
        if not _txt:
            continue
        _ax = float(node_source.data["x"][_i])
        _ay = float(node_source.data["y"][_i]) + _canvas_per_px * 6
        _w, _h = estimate_label_size(_txt, _node_px, _canvas_per_px)
        _node_anchors.append((_ax, _ay))
        _node_widths.append(_w); _node_heights.append(_h)
        _node_texts.append(_txt)

    # 2) 模块标签 (ID + 多行功能名 → 复合块)
    big_size_thr = max(20, int(sizes.most_common(1)[0][1] * 0.08))
    SMALL_MODULE_THR = 50
    id_xs, id_ys, id_texts, id_colors = [], [], [], []
    fn_lines_per_mod, fn_colors_per_mod = [], []
    for m, (mcx, mcy) in module_centers.items():
        func = module_function_map.get(m) or module_function_map.get(int(m))
        big_enough = sizes[m] >= big_size_thr
        if not (big_enough or func):
            continue
        col = darken_hex(palette[m] if m < len(palette) else "#64748b", 0.35)
        id_xs.append(mcx); id_ys.append(mcy)
        id_texts.append(f"M{m}  ·  n={sizes[m]}")
        id_colors.append(col)
        if isinstance(func, str):
            fn_list = [func]
        elif isinstance(func, (list, tuple)):
            fn_list = [str(x) for x in func]
        else:
            fn_list = []
        if sizes[m] < SMALL_MODULE_THR and len(fn_list) > 1:
            fn_list = fn_list[:1]
        fn_lines_per_mod.append(fn_list)
        fn_colors_per_mod.append(col)

    _mod_anchors, _mod_widths, _mod_heights = [], [], []
    for k in range(len(id_xs)):
        iw, ih = estimate_label_size(id_texts[k], _id_px, _canvas_per_px)
        fn_lines = fn_lines_per_mod[k]
        n_fn = len(fn_lines)
        if n_fn > 0:
            fn_wh = [estimate_label_size(t, _fn_px, _canvas_per_px) for t in fn_lines]
            max_fn_w = max(w for w, _ in fn_wh)
            comp_w = max(iw, max_fn_w)
            comp_h = _line_spacing_canvas * (1 + n_fn)
        else:
            comp_w = iw
            comp_h = ih
        _mod_anchors.append((id_xs[k], id_ys[k]))
        _mod_widths.append(comp_w); _mod_heights.append(comp_h)

    # 3) 拼成统一列表 → 一起去重叠
    _all_anchors = _node_anchors + _mod_anchors
    _all_widths = _node_widths + _mod_widths
    _all_heights = _node_heights + _mod_heights
    _n_nodes = len(_node_anchors)
    _is_static = [False] * _n_nodes + [True] * len(_mod_anchors)

    if _all_anchors:
        _new_pos = deoverlap_labels(
            _all_anchors, _all_widths, _all_heights,
            extra_pad=_canvas_per_px * 1.6,
            iterations=400, anchor_pull=0.04, push_step=0.5,
            is_static=_is_static,
        )

        # 用 label 实际位置 + bbox 扩展 figure 范围,避免文字被切
        _lbl_x_lo = min(_new_pos[k][0] - _all_widths[k] / 2 for k in range(len(_new_pos)))
        _lbl_x_hi = max(_new_pos[k][0] + _all_widths[k] / 2 for k in range(len(_new_pos)))
        _lbl_y_lo = min(_new_pos[k][1] - _all_heights[k] / 2 for k in range(len(_new_pos)))
        _lbl_y_hi = max(_new_pos[k][1] + _all_heights[k] / 2 for k in range(len(_new_pos)))
        _new_x_lo = min(x_lo, _lbl_x_lo)
        _new_x_hi = max(x_hi, _lbl_x_hi)
        _new_y_lo = min(y_lo, _lbl_y_lo)
        _new_y_hi = max(y_hi, _lbl_y_hi)
        if (_new_x_lo < x_lo - 1e-3 or _new_x_hi > x_hi + 1e-3 or
            _new_y_lo < y_lo - 1e-3 or _new_y_hi > y_hi + 1e-3):
            print(f"  ► 扩展 figure 范围以容纳 label "
                  f"(x:{x_lo:.0f}→{_new_x_lo:.0f},{x_hi:.0f}→{_new_x_hi:.0f}; "
                  f"y:{y_lo:.0f}→{_new_y_lo:.0f},{y_hi:.0f}→{_new_y_hi:.0f})")
            p.x_range.start = _new_x_lo
            p.x_range.end = _new_x_hi
            p.y_range.start = _new_y_lo
            p.y_range.end = _new_y_hi
            _new_dw = _new_x_hi - _new_x_lo
            _new_dh = _new_y_hi - _new_y_lo
            _new_aspect = _new_dw / max(_new_dh, 1e-6)
            _new_fw = int(round(math.sqrt(target_area * _new_aspect)))
            _new_fh = int(round(math.sqrt(target_area / _new_aspect)))
            _new_fw = max(_new_fw, 800)
            _new_fh = max(_new_fh, 800)
            p.width = _new_fw
            p.height = _new_fh

        # 拆分回去
        _new_node_pos = _new_pos[:_n_nodes]
        _new_mod_pos = _new_pos[_n_nodes:]

        # 节点标签层
        if _new_node_pos:
            _node_lbl_src = ColumnDataSource(dict(
                x=[p[0] for p in _new_node_pos],
                y=[p[1] for p in _new_node_pos],
                text=_node_texts,
            ))
            p.add_layout(LabelSet(
                x="x", y="y", text="text",
                source=_node_lbl_src,
                text_font_size=font_size,
                text_align="center", text_baseline="middle",
                text_color=text_color,
                text_font=Value(font_family),
                text_font_style="normal",
            ))

        # 模块标签层
        if _new_mod_pos:
            new_id_xs, new_id_ys = [], []
            new_fn_xs, new_fn_ys, new_fn_texts, new_fn_colors = [], [], [], []
            for k in range(len(id_xs)):
                cx, cy = _new_mod_pos[k]
                fn_lines = fn_lines_per_mod[k]
                n_fn = len(fn_lines)
                if n_fn == 0:
                    new_id_xs.append(cx); new_id_ys.append(cy)
                    continue
                top_y = cy + n_fn * _half_offset
                new_id_xs.append(cx)
                new_id_ys.append(top_y)
                for j, ftxt in enumerate(fn_lines):
                    fy = top_y - (j + 1) * _line_spacing_canvas
                    new_fn_xs.append(cx); new_fn_ys.append(fy)
                    new_fn_texts.append(ftxt)
                    new_fn_colors.append(fn_colors_per_mod[k])

            id_src = ColumnDataSource(dict(
                x=new_id_xs, y=new_id_ys, text=id_texts, color=id_colors,
            ))
            p.add_layout(LabelSet(
                x="x", y="y", text="text",
                source=id_src,
                text_font_size=font_size, text_font_style="normal",
                text_color="color",
                text_align="center", text_baseline="middle",
                text_font=Value(font_family),
                background_fill_color="white",
                background_fill_alpha=0.55,
                border_line_color=None,
            ))
            if new_fn_xs:
                fn_src = ColumnDataSource(dict(
                    x=new_fn_xs, y=new_fn_ys, text=new_fn_texts, color=new_fn_colors,
                ))
                p.add_layout(LabelSet(
                    x="x", y="y", text="text",
                    source=fn_src,
                    text_font_size=font_size, text_font_style="normal",
                    text_color="color",
                    text_align="center", text_baseline="middle",
                    text_font=Value(font_family),
                    background_fill_color="white",
                    background_fill_alpha=0.55,
                    border_line_color=None,
                ))

        print(f"  ► 文字去重叠完成 (统一处理 {_n_nodes} 节点标签 + {len(id_xs)} 模块标签)")

    # ── 保存 ──
    if outpath is not None:
        od = os.path.dirname(outpath)
        if od:
            os.makedirs(od, exist_ok=True)
        base = os.path.splitext(outpath)[0]

        html_path = base + ".html"
        output_file(html_path, title=title if show_title else "Modular network")
        save(p)
        print(f"✓ HTML 已保存: {html_path}")

        png_path = base + ".png"
        svg_path = base + ".svg"
        try:
            from PIL import Image
            Image.MAX_IMAGE_PIXELS = None
        except ImportError:
            pass

        driver = None
        try:
            from selenium import webdriver
            from selenium.webdriver.chrome.options import Options as ChromeOptions
            opts = ChromeOptions()
            opts.add_argument("--headless")
            opts.add_argument("--no-sandbox")
            opts.add_argument("--disable-dev-shm-usage")
            try:
                driver = webdriver.Chrome(options=opts)
            except Exception:
                from selenium.webdriver.firefox.options import Options as FirefoxOptions
                fo = FirefoxOptions()
                fo.add_argument("--headless")
                driver = webdriver.Firefox(options=fo)
        except Exception as e:
            print(f"⚠ 启动 webdriver 失败 ({e}),将尝试无 webdriver 导出")

        # PNG: 隐藏所有 LabelSet。位图作 figure.svg 的栅格背景使用,文字由 SVG
        # 矢量层提供。栅格化的 8pt 文字在缩放后糊成马赛克,且 Bokeh 把每个 text
        # 渲染成 fill + stroke 两份,叠加矢量层后会出现"皇冠"伪影。
        _label_annots = [a for a in p.center if isinstance(a, LabelSet)]
        for a in _label_annots:
            a.visible = False
        try:
            if driver is not None:
                export_png(p, filename=png_path, webdriver=driver)
            else:
                export_png(p, filename=png_path)
            print(f"✓ PNG 已保存: {png_path}  (label-free)")
        except Exception as e:
            print(f"⚠ PNG 导出失败: {e}")
        finally:
            for a in _label_annots:
                a.visible = True

        # SVG:把后端切到 svg 模式,再导出
        try:
            p.output_backend = "svg"
            if driver is not None:
                export_svg(p, filename=svg_path, webdriver=driver)
            else:
                export_svg(p, filename=svg_path)
            print(f"✓ SVG 已保存: {svg_path}")
        except Exception as e:
            print(f"⚠ SVG 导出失败: {e}")
        finally:
            try:
                p.output_backend = "canvas"
            except Exception:
                pass

        if driver is not None:
            try:
                driver.quit()
            except Exception:
                pass

        if save_module_csv:
            module_csv = base + "_modules.csv"
            pd.DataFrame({
                "id": [str(n) for n in G.nodes],
                "ntype": [G.nodes[n][node_type_attr] for n in G.nodes],
                "module": [module_id[n] for n in G.nodes],
                "degree": [deg_map[n] for n in G.nodes],
                "highlight": [int(str(n) in node_highlight) for n in G.nodes],
                "highlight_direction": [
                    "Up" if node_highlight.get(str(n), 0) > 0
                    else ("Down" if node_highlight.get(str(n), 0) < 0 else "")
                    for n in G.nodes
                ],
            }).to_csv(module_csv, index=False, encoding="utf-8")
            print(f"✓ 模块表已保存: {module_csv}")

    return G, module_id
