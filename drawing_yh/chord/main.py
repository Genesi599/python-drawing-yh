"""
drawing_yh.chord — directional chord diagram(基于 mpl_chord_diagram + patches)

适用场景:细胞-细胞通讯网络(sender → receiver,值 = 通讯强度 / event 数)、
任何 "from → to" 的有向加权图(基因调控、迁移流向、引文网络…)。

核心函数 `chord_diagram(matrix, ...)`:
    - matrix: pandas DataFrame,index = sender(源),columns = receiver(目标),值 = 权重
    - 自动按 (row + col 和) 分配每个 node 的 sector 大小
    - 同 sender 内 chord 按 receiver 距离排序(避免交叉)
    - chord 颜色 = sender 端深(sender_color) → receiver 端浅(lightened sender_color)
    - 每条 chord 在 receiver 端 sector arc 内层加一段 sender_color stripe(标识来源)
    - 返回 (fig, ax)

附带 `HEMATOPOIETIC_LINEAGE_COLORS` —— 造血谱系细胞类型的成系配色 dict(BM 项目复用):
    HSPC/progenitor=teal,B/Plasma=blue,DC=purple,T/NK=olive,
    Neutrophil=red,Mono/Mac=orange,Mk=magenta,Erythroid=brown

底层:vendored mpl_chord_diagram 0.4.1 in `_mpl_chord/`(已 patch 5 处:
allow gradient+directed、cend=lighten(sender)、asize×3、intra_gap、receiver inner stripe)
"""
from __future__ import annotations

import matplotlib.colors as _mc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle
import matplotlib.patheffects as pe

# 造血谱系成系配色(全彩色,不用灰/黑;teal=HSPC, blue=B, purple=DC, olive=T/NK,
#                   red=Neu, orange=Mono/Mac, magenta=Mk, brown=Erythroid)
HEMATOPOIETIC_LINEAGE_COLORS = {
    # HSPC / 早期 progenitor — teal/cyan
    'HSC/MPP': '#0d6b6b', 'CLP': '#179a9a', 'CMP': '#33b3b3', 'GMP': '#1a8585',
    'GP': '#26a0a0', 'MDP': '#4dc0c0', 'CDP': '#6bcccc', 'MEMP': '#0f7878',
    # B / Plasma — blue
    'Pro-B': '#1f5fb4', 'Pre-B': '#3679cc', 'ImmB': '#5b9bd9', 'BC': '#7bb3e1', 'Pla': '#0d3a7d',
    # DC — purple
    'pDC': '#7b3da0', 'cDC': '#a06bc0',
    # T / NK — olive / yellow-green
    'TC': '#7a8c1f', 'NKT': '#a3b552',
    # Neutrophil — red
    'Pro-Neu': '#a01818', 'Pre-Neu': '#c42a2a', 'ImmNeu': '#e04547', 'Neu': '#f17c7e',
    # Mono / Mac / EBM — orange
    'Mono': '#cc6600', 'Mac': '#e68a1a', 'EBM': '#ffaa55',
    # Mk — magenta / pink
    'MkP': '#a8186a', 'Mk/Gra': '#d96bb0',
    # Erythroid — brown
    'EryP': '#6b4226', 'Early Ery': '#8c5a3c', 'Late Ery': '#a87156', 'Ery/Gra': '#c79283',
}

# 通用 categorical 备用色(node 没在 color_map 里时按顺序取)
_FALLBACK_PALETTE = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
    '#e377c2', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
]


def lighten(hex_color, factor: float = 0.3):
    """向白色混合;factor=0 原色,factor=1 纯白。返回 (r, g, b) tuple。"""
    rgb = _mc.to_rgb(hex_color)
    return tuple(c + (1.0 - c) * factor for c in rgb)


def _as_square_matrix(mat: pd.DataFrame, order: list | None = None) -> tuple[pd.DataFrame, list]:
    """Align a directed matrix to a square node set."""
    univ = set(mat.index) | set(mat.columns)
    if order is not None:
        nodes = [n for n in order if n in univ] + [n for n in sorted(univ) if n not in set(order)]
    else:
        nodes = sorted(univ)
    return mat.reindex(index=nodes, columns=nodes, fill_value=0.0), nodes


def _add_directed_layout_epsilon(mat: pd.DataFrame) -> pd.DataFrame:
    """Avoid zero in/out degree divisions in the directed chord layout."""
    out = mat.copy()
    positive = out.to_numpy()[out.to_numpy() > 0]
    eps = float(positive.min()) * 1e-6 if len(positive) else 1e-9
    for node in out.index:
        row_sum = out.loc[node].sum()
        col_sum = out[node].sum()
        if (row_sum == 0 or col_sum == 0) and (row_sum + col_sum > 0):
            out.loc[node, node] = eps
    return out


def _slice_interval(start: float, end: float, weights: np.ndarray, gap: float) -> list[tuple[float, float]]:
    """Split a chord interval into weighted sub-intervals."""
    total = float(weights.sum())
    if total <= 0:
        return []
    reverse = end < start
    lo, hi = (end, start) if reverse else (start, end)
    span = max(hi - lo, 0.0)
    n = len(weights)
    use_gap = min(gap, span / max(n * 3, 1))
    usable = max(span - use_gap * max(n - 1, 0), 0.0)
    cur = lo
    parts: list[tuple[float, float]] = []
    for weight in weights:
        width = usable * float(weight) / total
        a, b = cur, cur + width
        parts.append((b, a) if reverse else (a, b))
        cur = b + use_gap
    return parts


def _allocate_role_arcs(
    group_nodes: list[str],
    degree_by_node: dict[str, float],
    span: tuple[float, float],
    pad: float,
) -> dict[str, tuple[float, float]]:
    """Allocate weighted node arcs inside one fixed role span."""
    if not group_nodes:
        return {}
    start, end = float(span[0]), float(span[1])
    direction = 1.0 if end >= start else -1.0
    span_width = abs(end - start)
    use_pad = min(float(pad), span_width / max(len(group_nodes) * 3, 1))
    usable = max(span_width - use_pad * max(len(group_nodes) - 1, 0), span_width * 0.12)
    weights = np.array([max(float(degree_by_node.get(node, 0.0)), 1e-9) for node in group_nodes])
    total = float(weights.sum())
    cur = start
    arcs: dict[str, tuple[float, float]] = {}
    for node, weight in zip(group_nodes, weights):
        arc_width = direction * usable * float(weight) / total
        arcs[node] = (cur, cur + arc_width)
        cur = cur + arc_width + direction * use_pad
    return arcs


def _role_label_position(start: float, end: float, radius: float = 1.09) -> tuple[tuple[float, float, float], bool]:
    mid = (start + end) / 2
    theta = np.deg2rad(mid)
    x = float(np.cos(theta) * radius)
    y = float(np.sin(theta) * radius)
    norm = mid % 360
    if -30 <= norm <= 180:
        return (x, y, mid - 90), False
    return (x, y, mid - 270), True


def _compute_role_split_positions(
    mat: pd.DataFrame,
    nodes: list[str],
    source_nodes: set[str],
    target_nodes: set[str],
    *,
    source_span: tuple[float, float],
    target_span: tuple[float, float],
    pad: float,
    sub_gap: float,
) -> tuple[list[tuple[float, float]], list[bool], list[tuple[float, float, float]], dict[tuple[int, int], tuple[float, float, float, float]]]:
    """Place source/input nodes in the upper span and target/output nodes in the lower span."""
    arr = mat.to_numpy(dtype=float)
    out_deg = arr.sum(axis=1)
    in_deg = arr.sum(axis=0)
    degree = out_deg + in_deg
    degree_by_node = {node: float(degree[i]) for i, node in enumerate(nodes)}

    source_group: list[str] = []
    target_group: list[str] = []
    for i, node in enumerate(nodes):
        is_source = node in source_nodes
        is_target = node in target_nodes
        if is_source and not is_target:
            source_group.append(node)
        elif is_target and not is_source:
            target_group.append(node)
        elif is_source and is_target:
            if out_deg[i] >= in_deg[i]:
                source_group.append(node)
            else:
                target_group.append(node)

    if not source_group or not target_group:
        raise ValueError("role_layout='split' requires at least one source/input node and one target/output node")

    arc_by_node = {
        **_allocate_role_arcs(source_group, degree_by_node, source_span, pad),
        **_allocate_role_arcs(target_group, degree_by_node, target_span, pad),
    }
    arc = [arc_by_node[node] for node in nodes]
    node_pos: list[tuple[float, float, float]] = []
    rotation: list[bool] = []
    for start, end in arc:
        text_pos, should_flip = _role_label_position(start, end)
        node_pos.append(text_pos)
        rotation.append(should_flip)

    angle_by_node = {node: (arc_by_node[node][0] + arc_by_node[node][1]) / 2 for node in nodes}
    source_segments: dict[tuple[int, int], tuple[float, float]] = {}
    target_segments: dict[tuple[int, int], tuple[float, float]] = {}
    for i, src in enumerate(nodes):
        js = [j for j in range(len(nodes)) if arr[i, j] > 0]
        js.sort(key=lambda j: angle_by_node[nodes[j]])
        weights = np.array([arr[i, j] for j in js], dtype=float)
        for j, segment in zip(js, _slice_interval(arc[i][0], arc[i][1], weights, sub_gap)):
            source_segments[(i, j)] = segment
    for j, dst in enumerate(nodes):
        is_ = [i for i in range(len(nodes)) if arr[i, j] > 0]
        is_.sort(key=lambda i: angle_by_node[nodes[i]])
        weights = np.array([arr[i, j] for i in is_], dtype=float)
        for i, segment in zip(is_, _slice_interval(arc[j][0], arc[j][1], weights, sub_gap)):
            target_segments[(i, j)] = segment

    pos: dict[tuple[int, int], tuple[float, float, float, float]] = {}
    for key, src_segment in source_segments.items():
        if key not in target_segments:
            continue
        dst_segment = target_segments[key]
        pos[key] = (*src_segment, *dst_segment)
    return arc, rotation, node_pos, pos


def chord_diagram(
    matrix: pd.DataFrame,
    color_map: dict | None = None,
    *,
    order: list | None = None,      # 显式 sector 顺序(如按谱系分组);None=按 node 名字母序(旧行为)
    top_nodes: list | None = None,  # 这些 node 涉及的 chord 提到最上层(避免被其它 chord 遮挡)
    figsize: tuple = (2.5, 2.5),    # 起始 figsize,autoshrink 会从这缩到 labels 刚好不重叠
    fontsize: int = 9,              # 固定 fontsize(不自动缩,由 figsize 调整防重叠)
    alpha: float = 1.0,             # chord 不透明,深
    pad: float = 3.0,               # sector 间空隙(度)
    chordwidth: float = 0.5,        # chord 弯曲度(0.3 直,0.5 适中,0.7 中段窄)
    intra_gap: float | None = None, # 同 sender 内 chord 子区段空隙(度);None=按 chord 总数自适应
    width: float = 0.1,             # sector arc 厚度(0-1)
    drop_zero_nodes: bool = True,
    radial_labels: bool = True,     # node label 径向(垂直 sector arc),节省版面
    **legacy_kwargs,                 # 静默吃掉所有旧 API 参数(title/title_size/space/lw/ec/direction/height_ratio/...)
):
    """画 directional chord diagram。

    Params
    ------
    matrix
        pandas DataFrame,index = sender,columns = receiver,值 = 权重(>=0)。
        非对称即有向(direction 默认 1,row → col)。
    color_map
        dict {node_name: color}。缺省按 HEMATOPOIETIC_LINEAGE_COLORS 找,
        找不到的 node 按 _FALLBACK_PALETTE 循环取色。
        显式传 {} 则全部用 fallback。
    figsize
        figure 尺寸(英寸)。默认 3.5×3.5(紧凑,单栏)。
    fontsize
        node label 字号。自动 overlap detection 重叠时缩小至最小 5pt。
    alpha
        chord 透明度。0.85 较深,0.5 较透。
    pad
        sector 之间空隙(度)。
    chordwidth
        chord 弯曲度。0.3 偏直,0.5 适中,0.7 中段窄。
    intra_gap
        同 sender 内多条 chord 子区段之间空隙(度)。
    width
        sector arc 在 figure radius 上的厚度。
    drop_zero_nodes
        True 时去掉行和 + 列和都为 0 的 node。

    Returns
    -------
    (fig, ax)
        matplotlib Figure + Axes。Title 不在图内,由 caller 加(如 fig.suptitle)。
        写出用 `drawing_yh.save_fig(fig, 'chord.pdf', also=('.png', '.svg'))`。

    Example
    -------
        import pandas as pd
        from drawing_yh.chord import chord_diagram
        from drawing_yh import save_fig
        fig, _ = chord_diagram(mat)
        save_fig(fig, 'out/chord.pdf', also=('.png', '.svg'))
    """
    from ._mpl_chord.chord_diagram import chord_diagram as _mpl_chord

    # 1. 对齐 + drop zero
    mat = matrix.copy()
    _univ = set(mat.index) | set(mat.columns)
    if order is not None:
        # 显式顺序优先(谱系分组等);order 里没列到的 node 按字母补到末尾
        all_nodes = [n for n in order if n in _univ] + [n for n in sorted(_univ) if n not in set(order)]
    else:
        all_nodes = sorted(_univ)
    mat = mat.reindex(index=all_nodes, columns=all_nodes, fill_value=0.0)
    if drop_zero_nodes:
        keep = [n for n in all_nodes if mat.loc[n].sum() > 0 or mat[n].sum() > 0]
        mat = mat.loc[keep, keep]
        all_nodes = keep

    # 2. 配色
    base_map = dict(HEMATOPOIETIC_LINEAGE_COLORS) if color_map is None else {}
    if color_map:
        base_map.update(color_map)
    cmap = {}
    fb_i = 0
    for n in all_nodes:
        if n in base_map:
            cmap[n] = base_map[n]
        else:
            cmap[n] = _FALLBACK_PALETTE[fb_i % len(_FALLBACK_PALETTE)]
            fb_i += 1
    colors = [cmap[n] for n in all_nodes]

    # 3. 画 chord(vendored mpl_chord_diagram with all patches)
    # alpha floor:legacy production scripts 传 alpha=0.42/0.55 太浅,统一 floor 到 0.85
    eff_alpha = max(alpha, 0.85)

    # intra_gap 自适应:chord 多则 gap 小(否则极小 sub-region 被 gap 吃掉)
    # 公式:总 chord 数 N → eff_gap = clip(8/N, 0.15, 1.5)
    #   5 chord  → 1.5 度;15 → 0.53;50 → 0.16;100+ → 0.15(floor)
    if intra_gap is None:
        import numpy as _np
        n_chords = int((_np.asarray(mat.values) > 0).sum())
        eff_intra_gap = max(0.15, min(1.5, 8.0 / max(n_chords, 1)))
    else:
        eff_intra_gap = intra_gap

    def _render(_fig_sz):
        _fig, _ax = plt.subplots(figsize=(_fig_sz, _fig_sz))
        plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
        _mpl_chord(mat.values, names=all_nodes, colors=colors, ax=_ax,
                   use_gradient=True,
                   chord_colors=colors,
                   directed=True,
                   sort="distance",
                   pad=pad,
                   chordwidth=chordwidth,
                   width=width,
                   gap=0,
                   intra_gap=eff_intra_gap,
                   rotate_names=radial_labels,
                   show=False,
                   fontsize=fontsize,
                   alpha=eff_alpha,
                   top_nodes=set(top_nodes) if top_nodes else None)
        return _fig, _ax

    # 通用 autoshrink — sectors 少时 min_size 大一点防 chord 被 labels 压扁
    from drawing_yh.layout import autoshrink_figsize
    n_sec = len(all_nodes)
    fig, ax = autoshrink_figsize(
        _render,
        initial=figsize[0],
        min_size=max(2.2, 2.2 + (n_sec - 6) * 0.04),
        bbox_shrink=0.5,
    )
    return fig, ax


def lr_role_chord_diagram(
    edges: pd.DataFrame,
    *,
    source: str = "source",
    target: str = "target",
    weight: str = "weight",
    edge_color: str = "edge_color",
    edge_alpha: str | None = None,
    edge_label: str | None = None,
    edge_number: str | None = None,
    edge_number_color: str | None = None,
    edge_number_side: str = "source",
    edge_number_radius: float = 0.90,
    edge_number_jitter: float = 0.0,
    edge_number_fontsize: int | None = None,
    edge_legend_color: str | None = None,
    source_role: str = "send",
    target_role: str = "recept",
    source_role_label: str | None = None,
    target_role_label: str | None = None,
    source_group_label: str | None = None,
    target_group_label: str | None = None,
    group_label_fontsize: int | None = None,
    group_label_radius: float = 1.34,
    label_map: dict | None = None,
    role_colors: dict | None = None,
    show_legend: bool = True,
    legend_title: str = "Cell type",
    legend_order: list | None = None,
    max_legend_items: int | None = None,
    order: list | None = None,
    top_nodes: list | None = None,
    role_layout: str = "split",
    source_span: tuple[float, float] = (8.0, 172.0),
    target_span: tuple[float, float] = (188.0, 352.0),
    role_span_background: bool = True,
    role_span_background_alpha: float = 0.18,
    figsize: tuple = (6.5, 6.5),
    fontsize: int = 8,
    alpha: float = 0.78,
    pad: float = 3.0,
    chordwidth: float = 0.55,
    intra_gap: float | None = None,
    sub_gap: float = 0.04,
    width: float = 0.08,
    radial_labels: bool = True,
    ax=None,
):
    """Draw an LR-style chord diagram with role-colored arcs and edge-colored ribbons.

    This template is for ligand/receptor or sender/receiver gene chords. By
    default source/input nodes occupy the upper arc and target/output nodes
    occupy the lower arc. Node color encodes role, while ribbon color encodes an
    edge-level category such as cell type. Multiple rows with the same
    source-target pair are drawn as parallel weighted ribbons instead of being
    collapsed into one color.
    """
    from ._mpl_chord.chord_diagram import chord_arc, ideogram_arc
    from ._mpl_chord.utilities import compute_positions

    required = {source, target, weight, edge_color}
    if edge_alpha is not None:
        required.add(edge_alpha)
    if edge_label is not None:
        required.add(edge_label)
    if edge_number is not None:
        required.add(edge_number)
    if edge_number_color is not None:
        required.add(edge_number_color)
    if edge_legend_color is not None:
        required.add(edge_legend_color)
    missing = required.difference(edges.columns)
    if missing:
        raise ValueError(f"edges is missing columns: {sorted(missing)}")
    if role_layout not in {"split", "circular"}:
        raise ValueError("role_layout must be 'split' or 'circular'")

    data = edges.copy()
    data = data[data[weight] > 0].copy()
    if data.empty:
        raise ValueError("edges has no positive-weight rows")

    role_colors = {
        source_role: "#F47C7C",
        target_role: "#FFD21A",
        **(role_colors or {}),
    }
    source_role_label = source_role_label or source_role
    target_role_label = target_role_label or target_role
    group_label_fontsize = group_label_fontsize or (fontsize + 2)
    label_map = {str(k): str(v) for k, v in (label_map or {}).items()}

    mat = data.pivot_table(
        index=source, columns=target, values=weight,
        aggfunc="sum", fill_value=0.0,
    )
    mat, nodes = _as_square_matrix(mat, order=order)
    keep = [n for n in nodes if mat.loc[n].sum() > 0 or mat[n].sum() > 0]
    mat = mat.loc[keep, keep]
    nodes = keep
    layout_mat = _add_directed_layout_epsilon(mat)

    source_nodes = set(data[source].astype(str))
    target_nodes = set(data[target].astype(str))
    arc_colors = []
    for node in nodes:
        if node in source_nodes and node not in target_nodes:
            arc_colors.append(role_colors[source_role])
        elif node in target_nodes and node not in source_nodes:
            arc_colors.append(role_colors[target_role])
        else:
            arc_colors.append(role_colors.get("both", "#BDBDBD"))

    def _render(fig_size: float):
        legend_extra = 1.35 if show_legend and ax is None else 0.0
        if ax is None:
            fig, plot_ax = plt.subplots(figsize=(fig_size + legend_extra, fig_size))
            if show_legend:
                plt.subplots_adjust(left=0.02, right=0.70, top=0.98, bottom=0.02)
            else:
                plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
        else:
            fig = ax.figure
            plot_ax = ax

        if intra_gap is None:
            n_chords = int((mat.to_numpy() > 0).sum())
            eff_intra_gap = max(0.15, min(1.5, 8.0 / max(n_chords, 1)))
        else:
            eff_intra_gap = intra_gap

        arc: list[tuple[float, float]] = []
        node_pos: list[tuple[float, float, float]] = []
        rotation: list[bool] = []
        pos: dict[tuple[int, int], tuple[float, float, float, float]] = {}
        if role_layout == "split":
            arc, rotation, node_pos, pos = _compute_role_split_positions(
                mat,
                nodes,
                source_nodes,
                target_nodes,
                source_span=source_span,
                target_span=target_span,
                pad=pad,
                sub_gap=eff_intra_gap,
            )
        else:
            arr = layout_mat.to_numpy(dtype=float)
            out_deg = arr.sum(axis=1)
            in_deg = arr.sum(axis=0)
            degree = out_deg + in_deg
            compute_positions(
                arr, degree, in_deg, out_deg, 0, False,
                {"sort": "distance", "intra_gap": eff_intra_gap},
                True, 360, pad, arc, rotation, node_pos, pos,
            )

        if role_layout == "split" and role_span_background:
            ideogram_arc(
                start=source_span[0], end=source_span[1], radius=1.0,
                color=role_colors[source_role], width=width, alpha=role_span_background_alpha,
                ax=plot_ax,
            )
            ideogram_arc(
                start=target_span[0], end=target_span[1], radius=1.0,
                color=role_colors[target_role], width=width, alpha=role_span_background_alpha,
                ax=plot_ax,
            )

        for i, color in enumerate(arc_colors):
            start, end = arc[i]
            ideogram_arc(start=start, end=end, radius=1.0, color=color,
                         width=width, alpha=1.0, ax=plot_ax)

        node_to_i = {node: i for i, node in enumerate(nodes)}
        top = set(top_nodes or [])
        grouped = data.groupby([source, target], sort=False)
        edge_rank = {str(v): i for i, v in enumerate(legend_order or [])}
        edge_number_fontsize_eff = edge_number_fontsize or max(5, fontsize - 2)
        edge_number_specs: list[dict] = []
        for (src, dst), group in grouped:
            if src not in node_to_i or dst not in node_to_i:
                continue
            if edge_label is not None and legend_order is not None:
                group = group.copy()
                group["_edge_rank"] = (
                    group[edge_label].astype(str).map(edge_rank).fillna(len(edge_rank))
                )
                group = group.sort_values(["_edge_rank", weight], ascending=[True, False])
            i, j = node_to_i[src], node_to_i[dst]
            start1, end1, start2, end2 = pos[(i, j)]
            weights = group[weight].to_numpy(dtype=float)
            src_parts = _slice_interval(start1, end1, weights, sub_gap)
            dst_parts = _slice_interval(start2, end2, weights, sub_gap)
            for (_, row), (s1, e1), (s2, e2) in zip(group.iterrows(), src_parts, dst_parts):
                color = row[edge_color]
                zorder = 6 if (src in top or dst in top) else 2
                row_alpha = float(row[edge_alpha]) if edge_alpha is not None else alpha
                chord_arc(
                    s1, e1, s2, e2,
                    radius=1 - width,
                    gap=0,
                    chordwidth=chordwidth,
                    color=color,
                    cend=color,
                    alpha=row_alpha,
                    ax=plot_ax,
                    use_gradient=False,
                    extent=360,
                    directed=True,
                    zorder=zorder,
                )
                if edge_number is not None:
                    number_value = row[edge_number]
                    if pd.isna(number_value) or str(number_value) in {"", "0", "nan", "None"}:
                        continue
                    if edge_number_side == "target":
                        theta = np.deg2rad((s2 + e2) / 2)
                        number_group = ("target", j)
                    else:
                        theta = np.deg2rad((s1 + e1) / 2)
                        number_group = ("source", i)
                    edge_number_specs.append({
                        "group": number_group,
                        "theta": theta,
                        "text": str(number_value),
                        "color": row[edge_number_color] if edge_number_color is not None else "black",
                        "zorder": zorder + 10,
                    })

        lim = max(1.24, group_label_radius + 0.08 if (source_group_label or target_group_label) else 1.24)
        plot_ax.set_xlim(-lim, lim)
        plot_ax.set_ylim(-lim, lim)
        plot_ax.set_aspect(1)

        number_groups: dict[tuple[str, int], list[dict]] = {}
        for spec in edge_number_specs:
            number_groups.setdefault(spec["group"], []).append(spec)
        def _alt_lane(index: int) -> int:
            """0, +1, -1, +2, -2... for local alternating number lanes."""
            if index == 0:
                return 0
            lane = (index + 1) // 2
            return lane if index % 2 else -lane

        for specs in number_groups.values():
            specs.sort(key=lambda s: s["theta"])
            for k, spec in enumerate(specs):
                spec["lane"] = _alt_lane(k)
                radius = edge_number_radius + spec["lane"] * edge_number_jitter
                spec["radius"] = min(1.08, max(0.74, radius))

        number_specs = [spec for specs in number_groups.values() for spec in specs]
        number_specs.sort(key=lambda s: (s["group"], s["theta"]))
        placed_number_bboxes = []
        radius_offsets = [0.0]
        angle_offsets = [0.0]
        if edge_number_jitter:
            radius_offsets = [edge_number_jitter * lane for lane in [
                0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5,
            ]]
            angle_offsets = np.deg2rad([
                0.0, 1.4, -1.4, 2.8, -2.8, 4.2, -4.2, 5.6, -5.6,
            ])

        def _draw_number(spec: dict, radius: float, theta: float):
            return plot_ax.text(
                np.cos(theta) * radius,
                np.sin(theta) * radius,
                spec["text"],
                ha="center",
                va="center",
                fontsize=edge_number_fontsize_eff,
                fontweight="bold",
                color=spec["color"],
                zorder=spec["zorder"],
                path_effects=[pe.withStroke(linewidth=2.2, foreground="white")],
            )

        for spec in number_specs:
            placed = None
            if edge_number_jitter:
                for dr in radius_offsets:
                    for dt in angle_offsets:
                        radius = min(1.18, max(0.64, spec["radius"] + dr))
                        theta = spec["theta"] + dt
                        txt = _draw_number(spec, radius, theta)
                        fig.canvas.draw()
                        bb = txt.get_window_extent(renderer=fig.canvas.get_renderer()).expanded(1.04, 1.08)
                        if not any(bb.overlaps(prev) for prev in placed_number_bboxes):
                            placed_number_bboxes.append(bb)
                            placed = txt
                            break
                        txt.remove()
                    if placed is not None:
                        break
            if placed is None:
                txt = _draw_number(spec, spec["radius"], spec["theta"])
                if edge_number_jitter:
                    fig.canvas.draw()
                    placed_number_bboxes.append(
                        txt.get_window_extent(renderer=fig.canvas.get_renderer()).expanded(1.04, 1.08)
                    )

        prop = {
            "fontsize": fontsize,
            "ha": "center",
            "va": "center",
            "rotation_mode": "anchor",
            "color": "black",
        }
        if radial_labels:
            for node, text_pos, should_flip in zip(nodes, node_pos, rotation):
                angle = text_pos[2]
                pp = prop.copy()
                rotate = 90
                arc_angle = np.average(arc[nodes.index(node)])
                if 90 < arc_angle < 180 or 270 < arc_angle:
                    rotate = -90
                if 90 < arc_angle < 270:
                    pp["ha"] = "right"
                else:
                    pp["ha"] = "left"
                plot_ax.text(
                    text_pos[0], text_pos[1], label_map.get(str(node), str(node)),
                    rotation=angle + rotate, **pp,
                )
        else:
            for node, text_pos, should_flip in zip(nodes, node_pos, rotation):
                pp = prop.copy()
                if should_flip:
                    pp["va"] = "top"
                else:
                    pp["va"] = "bottom"
                plot_ax.text(
                    text_pos[0], text_pos[1], label_map.get(str(node), str(node)),
                    rotation=text_pos[2], **pp,
                )

        def _add_group_label(group_nodes: set[str], label: str | None):
            if not label:
                return
            idx = [i for i, node in enumerate(nodes) if node in group_nodes]
            if not idx:
                return
            mids = np.deg2rad([np.mean(arc[i]) for i in idx])
            weights = np.array([max(abs(arc[i][1] - arc[i][0]), 0.1) for i in idx])
            x = float(np.sum(np.cos(mids) * weights))
            y = float(np.sum(np.sin(mids) * weights))
            angle = float(np.arctan2(y, x))
            lx = float(np.cos(angle) * group_label_radius)
            ly = float(np.sin(angle) * group_label_radius)
            # Keep the large group labels readable and phrase-like instead of radial.
            if abs(ly) > abs(lx):
                lx = 0.0
            else:
                ly = 0.0
            plot_ax.text(
                lx, ly, label,
                fontsize=group_label_fontsize,
                fontweight="semibold",
                ha="center", va="center",
                color="black",
                zorder=20,
            )

        _add_group_label(source_nodes, source_group_label)
        _add_group_label(target_nodes, target_group_label)

        plot_ax.set_xlim(-lim, lim)
        plot_ax.set_ylim(-lim, lim)
        plot_ax.set_aspect(1)
        plot_ax.axis("off")

        if show_legend:
            role_handles = [
                Patch(facecolor=role_colors[source_role], edgecolor="none", label=source_role_label),
                Patch(facecolor=role_colors[target_role], edgecolor="none", label=target_role_label),
            ]
            leg1 = fig.legend(
                handles=role_handles,
                title="Arc",
                frameon=False,
                loc="upper right",
                bbox_to_anchor=(0.985, 0.93),
                handlelength=1.1,
            )
            if edge_label is not None:
                legend_color_col = edge_legend_color or edge_color
                legend_df = (
                    data.groupby([edge_label, legend_color_col], as_index=False)[weight]
                    .sum()
                    .sort_values(weight, ascending=False)
                )
                if legend_order is not None:
                    rank = {str(v): i for i, v in enumerate(legend_order)}
                    legend_df["_legend_rank"] = (
                        legend_df[edge_label].astype(str).map(rank).fillna(len(rank))
                    )
                    legend_df = legend_df.sort_values(
                        ["_legend_rank", weight],
                        ascending=[True, False],
                    )
                if max_legend_items is not None:
                    legend_df = legend_df.head(max_legend_items)
                edge_handles = [
                    Patch(facecolor=row[legend_color_col], edgecolor="none", label=str(row[edge_label]))
                    for _, row in legend_df.iterrows()
                ]
                fig.legend(
                    handles=edge_handles,
                    title=legend_title,
                    frameon=False,
                    loc="upper right",
                    bbox_to_anchor=(0.985, 0.55),
                    handlelength=1.1,
                )
        return fig, plot_ax

    if ax is not None:
        return _render(figsize[0])

    from drawing_yh.layout import autoshrink_figsize
    fig, ax = autoshrink_figsize(
        _render,
        initial=figsize[0],
        min_size=max(2.2, 2.2 + (len(nodes) - 6) * 0.04),
        max_size=max(figsize[0], 10.5),
        bbox_shrink=0.72,
    )
    return fig, ax


def _first_seen(values: pd.Series | list) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        if pd.isna(value):
            continue
        text = str(value)
        if text not in seen:
            out.append(text)
            seen.add(text)
    return out


def _legend_records(
    data: pd.DataFrame,
    label_col: str,
    color_col: str,
    *,
    order: list | None = None,
    number_col: str | None = None,
) -> list[dict]:
    records = []
    rank = {str(v): i for i, v in enumerate(order or [])}
    work = data[[label_col, color_col] + ([number_col] if number_col else [])].dropna(subset=[label_col])
    work = work.drop_duplicates(subset=[label_col, color_col] + ([number_col] if number_col else []))
    if order is not None:
        work["_rank"] = work[label_col].astype(str).map(rank).fillna(len(rank))
        work = work.sort_values(["_rank", label_col])
    for _, row in work.iterrows():
        rec = {"label": str(row[label_col]), "color": row[color_col]}
        if number_col:
            value = row[number_col]
            rec["number"] = "" if pd.isna(value) else str(value)
        records.append(rec)
    return records


def _draw_lr_role_chord_panel_legend(
    ax,
    data: pd.DataFrame,
    *,
    role_colors: dict,
    source_role: str,
    target_role: str,
    source_role_label: str,
    target_role_label: str,
    edge_label: str | None,
    edge_color: str,
    edge_legend_color: str | None,
    edge_legend_title: str,
    edge_legend_order: list | None,
    edge_number: str | None,
    status: str | None,
    status_palette: dict | None,
    status_title: str,
    number_color_legend: dict | None,
    number_color_title: str,
    title: str,
    fontsize: int,
    heading_fontsize: int,
    title_fontsize: int,
) -> None:
    """Draw a compact shared legend for lr_role_chord_panel."""
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    y = 0.95
    ax.text(0.50, y, title, ha="center", va="center",
            fontsize=title_fontsize, fontweight="bold", transform=ax.transAxes)
    y -= 0.070

    def heading(text: str) -> None:
        nonlocal y
        ax.text(0.10, y, text, fontsize=heading_fontsize, fontweight="bold",
                transform=ax.transAxes)
        y -= 0.040

    def swatch(label: str, color: str, *, edge: str = "none") -> None:
        nonlocal y
        ax.add_patch(Rectangle((0.10, y - 0.014), 0.035, 0.028,
                               facecolor=color, edgecolor=edge, linewidth=0.6,
                               transform=ax.transAxes))
        ax.text(0.16, y, label, fontsize=fontsize, va="center", transform=ax.transAxes)
        y -= 0.038

    heading("Arc")
    swatch(source_role_label, role_colors[source_role])
    swatch(target_role_label, role_colors[target_role])
    y -= 0.024

    if status and status_palette:
        heading(status_title)
        for label, color in status_palette.items():
            edge = "#777777" if str(color).upper() in {"#F0F0F0", "#F2F2F2", "#D0D0D0"} else "none"
            swatch(str(label), color, edge=edge)
        y -= 0.024

    if number_color_legend:
        heading(number_color_title)
        for label, color in number_color_legend.items():
            ax.text(0.10, y, "12", fontsize=fontsize, fontweight="bold",
                    color=color, va="center", transform=ax.transAxes)
            ax.text(0.18, y, str(label), fontsize=fontsize, va="center", transform=ax.transAxes)
            y -= 0.038
        y -= 0.024

    if edge_label:
        heading(edge_legend_title)
        color_col = edge_legend_color or edge_color
        records = _legend_records(
            data,
            edge_label,
            color_col,
            order=edge_legend_order,
            number_col=edge_number,
        )
        if len(records) > 12:
            n_left = (len(records) + 1) // 2
            col_x = [0.10, 0.56]
            row_y0 = y
            row_step = 0.030
            for idx, rec in enumerate(records):
                col = 0 if idx < n_left else 1
                row = idx if idx < n_left else idx - n_left
                yy = row_y0 - row * row_step
                xx = col_x[col]
                if rec.get("number"):
                    ax.text(xx, yy, rec["number"], fontsize=fontsize - 1,
                            fontweight="bold", va="center", transform=ax.transAxes)
                    patch_x = xx + 0.048
                    label_x = xx + 0.088
                else:
                    patch_x = xx
                    label_x = xx + 0.046
                ax.add_patch(Rectangle((patch_x, yy - 0.010), 0.030, 0.021,
                                       facecolor=rec["color"], edgecolor="none",
                                       transform=ax.transAxes))
                ax.text(label_x, yy, rec["label"], fontsize=fontsize - 1,
                        va="center", transform=ax.transAxes)
        else:
            for rec in records:
                if rec.get("number"):
                    ax.text(0.10, y, rec["number"], fontsize=fontsize,
                            fontweight="bold", va="center", transform=ax.transAxes)
                    patch_x = 0.18
                    label_x = 0.23
                else:
                    patch_x = 0.10
                    label_x = 0.16
                ax.add_patch(Rectangle((patch_x, y - 0.014), 0.035, 0.028,
                                       facecolor=rec["color"], edgecolor="none",
                                       transform=ax.transAxes))
                ax.text(label_x, y, rec["label"], fontsize=fontsize,
                        va="center", transform=ax.transAxes)
                y -= 0.038


def lr_role_chord_panel(
    edges: pd.DataFrame,
    *,
    panel: str = "panel",
    panel_order: list | None = None,
    source: str = "source",
    target: str = "target",
    weight: str = "weight",
    edge_color: str = "edge_color",
    edge_alpha: str | None = None,
    edge_label: str | None = None,
    edge_number: str | None = None,
    edge_number_color: str | None = None,
    edge_number_side: str = "source",
    edge_number_radius: float = 0.90,
    edge_number_jitter: float = 0.0,
    edge_legend_color: str | None = None,
    edge_legend_title: str = "Edge group",
    edge_legend_order: list | None = None,
    status: str | None = None,
    status_palette: dict | None = None,
    status_title: str = "Status",
    number_color_legend: dict | None = None,
    number_color_title: str = "Number color",
    source_role: str = "send",
    target_role: str = "recept",
    source_role_label: str | None = None,
    target_role_label: str | None = None,
    source_group_label: str | None = None,
    target_group_label: str | None = None,
    group_label_fontsize: int | None = None,
    group_label_radius: float = 1.50,
    role_colors: dict | None = None,
    label_map: dict | None = None,
    order: list | None = None,
    top_nodes: list | None = None,
    role_layout: str = "split",
    source_span: tuple[float, float] = (8.0, 172.0),
    target_span: tuple[float, float] = (188.0, 352.0),
    role_span_background: bool = True,
    role_span_background_alpha: float = 0.18,
    title: str | None = None,
    legend_title: str = "Legend",
    figsize: tuple | None = None,
    chord_size: float = 5.2,
    legend_width_ratio: float = 0.70,
    fontsize: int = 8,
    panel_title_fontsize: int = 12,
    legend_fontsize: int = 8,
    legend_heading_fontsize: int = 9,
    legend_title_fontsize: int = 13,
    pad: float = 3.0,
    chordwidth: float = 0.55,
    intra_gap: float | None = None,
    sub_gap: float = 0.04,
    width: float = 0.08,
    radial_labels: bool = True,
    wspace: float = 0.0,
    panel_top: float | None = None,
):
    """Draw side-by-side LR chord panels with one shared legend.

    `edges` is a plotted edge table. Each row is one ribbon in one panel.
    The function keeps panel-specific rendering in separate axes while sharing
    arc, status, number-color, and edge-group legends in the rightmost column.
    It is intended for templates such as aging / treatment A / treatment B
    comparisons that use the same LR-cell edge schema.
    """
    required = {panel, source, target, weight, edge_color}
    optional = [
        edge_alpha, edge_label, edge_number, edge_number_color,
        edge_legend_color, status,
    ]
    required.update(c for c in optional if c)
    missing = required.difference(edges.columns)
    if missing:
        raise ValueError(f"edges is missing columns: {sorted(missing)}")

    data = edges.copy()
    data = data[data[weight] > 0].copy()
    if data.empty:
        raise ValueError("edges has no positive-weight rows")
    data[panel] = data[panel].astype(str)

    panels = [str(x) for x in (panel_order or _first_seen(data[panel]))]
    panels = [p for p in panels if p in set(data[panel])]
    if not panels:
        raise ValueError("No panels from panel_order are present in edges")

    role_colors = {
        source_role: "#F47C7C",
        target_role: "#FFD21A",
        **(role_colors or {}),
    }
    source_role_label = source_role_label or source_role
    target_role_label = target_role_label or target_role

    if figsize is None:
        figsize = (chord_size * (len(panels) + legend_width_ratio), chord_size)
    fig = plt.figure(figsize=figsize, constrained_layout=False)
    gs = fig.add_gridspec(
        1,
        len(panels) + 1,
        width_ratios=[1.0] * len(panels) + [legend_width_ratio],
        left=0.012,
        right=0.992,
        top=panel_top if panel_top is not None else (0.82 if title else 0.84),
        bottom=0.060,
        wspace=wspace,
    )
    axes = [fig.add_subplot(gs[0, i]) for i in range(len(panels))]
    legend_ax = fig.add_subplot(gs[0, len(panels)])

    for ax, panel_name in zip(axes, panels):
        panel_edges = data[data[panel].eq(panel_name)].copy()
        lr_role_chord_diagram(
            panel_edges,
            source=source,
            target=target,
            weight=weight,
            edge_color=edge_color,
            edge_alpha=edge_alpha,
            edge_label=edge_label,
            edge_number=edge_number,
            edge_number_color=edge_number_color,
            edge_number_side=edge_number_side,
            edge_number_radius=edge_number_radius,
            edge_number_jitter=edge_number_jitter,
            edge_legend_color=edge_legend_color,
            source_role=source_role,
            target_role=target_role,
            source_role_label=source_role_label,
            target_role_label=target_role_label,
            source_group_label=source_group_label,
            target_group_label=target_group_label,
            group_label_fontsize=group_label_fontsize,
            group_label_radius=group_label_radius,
            role_colors=role_colors,
            label_map=label_map,
            show_legend=False,
            legend_order=edge_legend_order,
            order=order,
            top_nodes=top_nodes,
            role_layout=role_layout,
            source_span=source_span,
            target_span=target_span,
            role_span_background=role_span_background,
            role_span_background_alpha=role_span_background_alpha,
            figsize=(chord_size, chord_size),
            fontsize=fontsize,
            pad=pad,
            chordwidth=chordwidth,
            intra_gap=intra_gap,
            sub_gap=sub_gap,
            width=width,
            radial_labels=radial_labels,
            ax=ax,
        )
        ax.set_title(panel_name, fontsize=panel_title_fontsize, fontweight="bold", pad=10)

    _draw_lr_role_chord_panel_legend(
        legend_ax,
        data,
        role_colors=role_colors,
        source_role=source_role,
        target_role=target_role,
        source_role_label=source_role_label,
        target_role_label=target_role_label,
        edge_label=edge_label,
        edge_color=edge_color,
        edge_legend_color=edge_legend_color,
        edge_legend_title=edge_legend_title,
        edge_legend_order=edge_legend_order,
        edge_number=edge_number,
        status=status,
        status_palette=status_palette,
        status_title=status_title,
        number_color_legend=number_color_legend,
        number_color_title=number_color_title,
        title=legend_title,
        fontsize=legend_fontsize,
        heading_fontsize=legend_heading_fontsize,
        title_fontsize=legend_title_fontsize,
    )
    if title:
        fig.suptitle(title, fontsize=panel_title_fontsize + 2, fontweight="bold", y=0.985)
    return fig, axes, legend_ax
