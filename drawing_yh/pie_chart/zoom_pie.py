
import colorsys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from pathlib import Path

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.use14corefonts'] = False
mpl.rcParams['font.family']  = 'Arial'

DEFAULT_COLORS = [
    '#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3',
    '#937860', '#DA8BC3', '#8C8C8C', '#CCB974', '#64B5CD',
    '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF', '#AEC7E8',
]


def generate_shades(base_color: str, n: int, l_range=(0.55, 0.82)) -> list:
    h = base_color.lstrip('#')
    r, g, b = (int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    hue, l, s = colorsys.rgb_to_hls(r, g, b)
    s = min(s * 1.05, 1.0)
    shades = []
    for i in range(n):
        t = i / max(n - 1, 1)
        l_new = l_range[0] + t * (l_range[1] - l_range[0])
        r2, g2, b2 = colorsys.hls_to_rgb(hue, l_new, s)
        shades.append('#{:02X}{:02X}{:02X}'.format(
            int(r2 * 255), int(g2 * 255), int(b2 * 255)))
    return shades


def _merge_small(values, labels, colors, threshold_pct, keep_index=None, other_color='#BBBBBB'):
    total    = sum(values)
    kept     = []
    other_v  = 0.0
    new_keep = None
    for i, (v, l, c) in enumerate(zip(values, labels, colors)):
        if i == keep_index or v / total * 100 >= threshold_pct:
            if i == keep_index:
                new_keep = len(kept)
            kept.append((v, l, c))
        else:
            other_v += v
    if other_v > 0:
        kept.append((other_v, 'Others', other_color))
    vv, ll, cc = zip(*kept) if kept else ([], [], [])
    return list(vv), list(ll), list(cc), new_keep


def _startangle_so_slice_is_right(values, zoom_index):
    total    = sum(values)
    fracs    = [v / total for v in values]
    mid_frac = sum(fracs[:zoom_index]) + fracs[zoom_index] / 2
    return -mid_frac * 360


def _draw_pie_on_ax(ax, values, labels, colors,
                    show_pct, show_count, pct_distance, label_distance,
                    startangle, font_size, title,
                    min_pct_for_label=3.0,
                    suppress_label_indices=None):
    values   = np.array(values, dtype=float)
    total    = values.sum()
    pcts     = values / total * 100
    suppress = set(suppress_label_indices or [])

    def _make_label(i, l):
        if i in suppress or pcts[i] < min_pct_for_label:
            return ''
        lines = [l]
        if show_pct:
            pct_str = f'{pcts[i]:.1f}%'
            if show_count:
                pct_str += f' (n={int(round(pcts[i]/100*total))})'
            lines.append(pct_str)
        return '\n'.join(lines)

    display_labels = [_make_label(i, l) for i, l in enumerate(labels)]

    result = ax.pie(
        values,
        labels=display_labels,
        colors=colors,
        autopct=None,
        pctdistance=pct_distance,
        labeldistance=label_distance,
        startangle=startangle,
        wedgeprops=dict(linewidth=1.0, edgecolor='white'),
        textprops=dict(fontsize=font_size, fontfamily='Arial'),
    )
    wedges, texts = result[0], result[1]

    # 记录每个标签对应的扇形中心方向（单位向量），用于后续画引导线
    mid_angles = []
    for wedge, t in zip(wedges, texts):
        mid_angle = np.radians((wedge.theta1 + wedge.theta2) / 2)
        mid_angles.append(mid_angle)
        if not t.get_text().strip():
            continue
        x = np.cos(mid_angle) * label_distance
        y = np.sin(mid_angle) * label_distance
        t.set_position((x, y))
        t.set_ha('center')
        t.set_va('center')

    for t in texts:
        t.set_fontsize(font_size)
        t.set_fontweight('normal')
        t.set_multialignment('left')

    if title:
        ax.set_title(title, fontsize=font_size, fontweight='normal',
                     fontfamily='Arial', pad=8)
    ax.set_aspect('equal')
    return wedges, texts, mid_angles


def _ax_to_fig(fig, ax, xy):
    disp = ax.transData.transform(xy)
    return fig.transFigure.inverted().transform(disp)


def _fix_label_overlaps(fig, ax, texts, margin=3.0, max_iter=50):
    texts = [t for t in texts if t.get_text().strip()]
    if len(texts) < 2:
        return

    fig.canvas.draw()

    p0 = ax.transData.transform([0, 0])
    p1 = ax.transData.transform([1, 0])
    px_per_unit_x = abs(p1[0] - p0[0])
    p1y = ax.transData.transform([0, 1])
    px_per_unit_y = abs(p1y[1] - p0[1])

    for _ in range(max_iter):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        bboxes = [t.get_window_extent(renderer) for t in texts]
        moved = False

        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                bi, bj = bboxes[i], bboxes[j]
                if not bi.overlaps(bj):
                    continue

                overlap_x = (min(bi.x1, bj.x1) - max(bi.x0, bj.x0)) + margin
                overlap_y = (min(bi.y1, bj.y1) - max(bi.y0, bj.y0)) + margin

                dx = overlap_x / px_per_unit_x / 2
                dy = overlap_y / px_per_unit_y / 2

                xi, yi = texts[i].get_position()
                xj, yj = texts[j].get_position()

                vx, vy = xi - xj, yi - yj
                norm = max(np.hypot(vx, vy), 1e-9)
                vx, vy = vx / norm, vy / norm

                texts[i].set_position((xi + vx * dx, yi + vy * dy))
                texts[j].set_position((xj - vx * dx, yj - vy * dy))
                moved = True

        if not moved:
            break

    fig.canvas.draw()


def _fix_cross_axes_overlaps(fig, ax_left, left_texts, ax_right, right_texts,
                              margin=3.0, max_iter=30):
    """修正两个 axes 之间的标签重叠。"""
    lt = [t for t in left_texts if t.get_text().strip()]
    rt = [t for t in right_texts if t.get_text().strip()]
    if not lt or not rt:
        return

    fig.canvas.draw()

    for _ in range(max_iter):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        moved = False

        for tl in lt:
            bl = tl.get_window_extent(renderer)
            for tr in rt:
                br = tr.get_window_extent(renderer)
                if not bl.overlaps(br):
                    continue

                overlap_x = (min(bl.x1, br.x1) - max(bl.x0, br.x0)) + margin
                overlap_y = (min(bl.y1, br.y1) - max(bl.y0, br.y0)) + margin

                cl_d = np.array([(bl.x0 + bl.x1) / 2, (bl.y0 + bl.y1) / 2])
                cr_d = np.array([(br.x0 + br.x1) / 2, (br.y0 + br.y1) / 2])
                v = cl_d - cr_d
                norm = max(np.hypot(v[0], v[1]), 1e-9)
                v = v / norm

                push_px = np.array([v[0] * overlap_x, v[1] * overlap_y]) / 2

                pos_l_disp = ax_left.transData.transform(tl.get_position())
                tl.set_position(ax_left.transData.inverted().transform(
                    pos_l_disp + push_px))

                pos_r_disp = ax_right.transData.transform(tr.get_position())
                tr.set_position(ax_right.transData.inverted().transform(
                    pos_r_disp - push_px))

                moved = True

        if not moved:
            break

    fig.canvas.draw()


def _draw_leader_lines(ax, wedges, texts, mid_angles, rim=0.6):
    """在标签与扇形边缘之间画虚线引导线，标签位置以文字当前位置为准。"""
    for wedge, t, angle in zip(wedges, texts, mid_angles):
        if not t.get_text().strip():
            continue
        tx, ty = t.get_position()
        # 扇形边缘点（沿中心角方向）
        ex = np.cos(angle) * rim
        ey = np.sin(angle) * rim
        # 只在标签离扇形边缘有一定距离时才画
        dist = np.hypot(tx - ex, ty - ey)
        if dist < 0.05:
            continue
        ax.plot(
            [ex, tx], [ey, ty],
            color='#555555', linewidth=0.6, linestyle='--', alpha=0.85, zorder=1, clip_on=False,
        )


def plot_zoom_pie(
    major_values: list,
    major_labels: list,
    detail_values: list,
    detail_labels: list,
    zoom_index: int,
    major_colors: list      = None,
    detail_colors: list     = None,
    title_left: str         = '',
    title_right: str        = '',
    show_pct: bool          = True,
    show_count: bool        = False,
    pct_distance: float     = 0.70,
    label_distance: float   = 0.6,
    font_size: float        = 8,
    merge_threshold: float  = 3.0,
    connect_alpha: float    = 0.18,
    figsize: tuple          = None,
    out_path: str           = None,
    show_leader_lines: bool = True,
):
    n_major = len(major_values)
    if major_colors is None:
        major_colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(n_major)]
    base_color = major_colors[zoom_index]

    m_vals, m_labs, m_cols, new_zoom_idx = _merge_small(
        major_values, major_labels, major_colors,
        merge_threshold, keep_index=zoom_index,
    )

    others_idx = next((i for i, l in enumerate(m_labs) if l == 'Others'), None)
    sort_indices = sorted(
        [i for i in range(len(m_vals)) if i != others_idx],
        key=lambda i: m_vals[i], reverse=True
    )
    if others_idx is not None:
        sort_indices.append(others_idx)

    m_vals = [m_vals[i] for i in sort_indices]
    m_labs = [m_labs[i] for i in sort_indices]
    m_cols = [m_cols[i] for i in sort_indices]
    new_zoom_idx = sort_indices.index(new_zoom_idx)

    n_detail_raw = len(detail_values)
    dc_raw = generate_shades(base_color, n_detail_raw) if detail_colors is None \
             else list(detail_colors)
    d_vals, d_labs, d_cols, _ = _merge_small(
        detail_values, detail_labels, dc_raw, merge_threshold,
        other_color=generate_shades(base_color, 1, l_range=(0.88, 0.88))[0],
    )

    left_startangle = _startangle_so_slice_is_right(m_vals, new_zoom_idx)

    # ── 布局度量 ───────────────────────────────────────────────────
    pt2in   = 1 / 72
    cw      = font_size * 0.55 * pt2in
    max_l   = max((len(l) for l in m_labs), default=8)
    max_r   = max((len(l) for l in d_labs), default=8)
    pie_r   = 1.0
    w_l     = 2 * (label_distance * pie_r + max_l * cw + 0.15)
    w_r     = 2 * (label_distance * pie_r + max_r * cw + 0.15) * 0.82

    if figsize is None:
        figsize = (round(w_l + w_r + 0.2, 1),
                   round(max(w_l * 0.82, 2.4), 1))

    fig = plt.figure(figsize=figsize, facecolor='white')

    # ── 初始 axes（先放置再根据实际渲染尺寸校正） ─────────────────
    fig_w    = figsize[0]
    pad      = 0.01
    frac_l   = w_l / fig_w
    frac_r   = w_r / fig_w

    ax_left  = fig.add_axes([pad, 0.06, frac_l, 0.88])
    ax_right = fig.add_axes([pad + frac_l, 0.11, frac_r, 0.78])
    ax_left.set_zorder(2)
    ax_left.patch.set_visible(False)
    ax_right.set_zorder(1)

    left_wedges, left_texts, left_angles = _draw_pie_on_ax(
        ax_left, m_vals, m_labs, m_cols,
        show_pct=show_pct, show_count=show_count,
        pct_distance=pct_distance, label_distance=label_distance,
        startangle=left_startangle, font_size=font_size, title=title_left,
        min_pct_for_label=0,
    )

    right_wedges, right_texts, right_angles = _draw_pie_on_ax(
        ax_right, d_vals, d_labs, d_cols,
        show_pct=show_pct, show_count=show_count,
        pct_distance=pct_distance, label_distance=label_distance,
        startangle=90.0, font_size=font_size, title=title_right,
        min_pct_for_label=0,
    )

    # ── 测量实际饼半径，校正右图位置使间距 = 左饼半径 ─────────────
    fig.canvas.draw()
    _cl = ax_left.transData.transform([0, 0])
    _el = ax_left.transData.transform([1, 0])
    _cr = ax_right.transData.transform([0, 0])
    _er = ax_right.transData.transform([-1, 0])
    r_l_px = _el[0] - _cl[0]
    r_r_px = _cr[0] - _er[0]
    current_gap_px = (_cr[0] - r_r_px) - (_cl[0] + r_l_px)
    shift_fig = (current_gap_px - r_l_px * 0.5) / (fig.dpi * figsize[0])
    pos = ax_right.get_position()
    ax_right.set_position([pos.x0 - shift_fig, pos.y0,
                           pos.width, pos.height])

    fig.canvas.draw()
    _fix_label_overlaps(fig, ax_left,  left_texts,  margin=8.0)
    _fix_label_overlaps(fig, ax_right, right_texts, margin=8.0)
    _fix_cross_axes_overlaps(fig, ax_left, left_texts,
                             ax_right, right_texts, margin=8.0)

    # 引导线在 overlap 修正后画，位置已是最终位置
    if show_leader_lines:
        _draw_leader_lines(ax_left,  left_wedges,  left_texts,  left_angles)
        _draw_leader_lines(ax_right, right_wedges, right_texts, right_angles)

    for t in left_texts:
        t.set_zorder(10)

    wedge  = left_wedges[new_zoom_idx]
    theta1 = np.radians(wedge.theta1)
    theta2 = np.radians(wedge.theta2)
    if theta2 < theta1:
        theta1, theta2 = theta2, theta1

    n_arc = 80

    def _d2d(ax, xy):
        return np.array(ax.transData.transform(xy))
    def _d2f(xy):
        return np.array(fig.transFigure.inverted().transform(xy))

    left_thetas = np.linspace(theta1, theta2, n_arc)
    left_arc_d  = [_d2d(ax_left, [np.cos(t), np.sin(t)]) for t in left_thetas]

    c_r   = _d2d(ax_right, [0, 0])
    rim_r = _d2d(ax_right, [1, 0])
    r_r   = np.linalg.norm(rim_r - c_r)

    def _tangent(P, C, r):
        dv = P - C
        d  = np.linalg.norm(dv)
        if d <= r * 1.01:
            u = dv / max(d, 1e-10)
            return C + r * u, C + r * u
        beta  = np.arctan2(dv[1], dv[0])
        delta = np.arccos(np.clip(r / d, -1, 1))
        T1 = C + r * np.array([np.cos(beta + delta), np.sin(beta + delta)])
        T2 = C + r * np.array([np.cos(beta - delta), np.sin(beta - delta)])
        return T1, T2

    p_lo_d = left_arc_d[0]
    p_hi_d = left_arc_d[-1]

    T1_lo, T2_lo = _tangent(p_lo_d, c_r, r_r)
    t_lo_d = T1_lo if T1_lo[1] < T2_lo[1] else T2_lo

    T1_hi, T2_hi = _tangent(p_hi_d, c_r, r_r)
    t_hi_d = T1_hi if T1_hi[1] > T2_hi[1] else T2_hi

    a_hi = np.arctan2((t_hi_d - c_r)[1], (t_hi_d - c_r)[0])
    a_lo = np.arctan2((t_lo_d - c_r)[1], (t_lo_d - c_r)[0])
    if a_lo < a_hi:
        a_lo += 2 * np.pi
    right_arc_d = [c_r + r_r * np.array([np.cos(a), np.sin(a)])
                   for a in np.linspace(a_hi, a_lo, n_arc)]

    shape_f = [_d2f(p) for p in (left_arc_d + right_arc_d)]

    connector = Polygon(
        shape_f, closed=True,
        facecolor=base_color, alpha=connect_alpha,
        edgecolor=base_color, linewidth=0.8,
        transform=fig.transFigure, zorder=0,
    )
    fig.add_artist(connector)

    if out_path:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        for suffix in ['.png', '.pdf', '.svg']:
            p = out_path.with_suffix(suffix)
            dpi = 600 if suffix == '.png' else None
            fig.savefig(p, dpi=dpi, bbox_inches='tight', facecolor='white')
            print(f"Saved: {p}")

    return fig, ax_left, ax_right
