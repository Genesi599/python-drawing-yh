import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

mpl.rcParams['pdf.fonttype'] = 42       # 字体嵌入 PDF（Illustrator 可编辑）
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['svg.fonttype'] = 'none'   # SVG 保留文字（Illustrator 可编辑）
mpl.rcParams['pdf.use14corefonts'] = False
mpl.rcParams['font.family']  = 'Arial'

DEFAULT_COLORS = [
    '#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3',
    '#937860', '#DA8BC3', '#8C8C8C', '#CCB974', '#64B5CD',
    '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF', '#AEC7E8',
]


def _auto_figsize(
    n_direct_labels: int,
    max_label_chars: int,
    label_font_size: float,
    label_distance: float,
    title: str,
    has_legend: bool,
) -> tuple:
    """
    根据直接标注的标签数量、最长标签字符数和字体大小自动推算 figsize（单位：英寸）。

    思路
    ----
    · 固定饼图半径 pie_r = 3 英寸作为参考基准
    · 标签起点距圆心 = label_distance * pie_r
    · 标签文字宽度 ≈ max_chars * font_size_pt * 0.55 / 72
    · 所需宽度   = 2 * (label_distance * pie_r + text_w + margin)
    · 所需高度：n_direct 个标签均匀分布，最密区域在上下两侧；
                保证相邻标签行距 ≥ 1.5 * font_size_pt / 72 英寸
    · 如果有 legend，右侧额外留 legend_w 的空间
    """
    pt_per_inch = 72.0
    pie_r       = 2.0          # 饼图半径（英寸）
    char_w_inch = label_font_size * 0.55 / pt_per_inch
    text_w      = max_label_chars * char_w_inch
    margin      = 0.5

    legend_w = 1.8 if has_legend else 0.0

    w = 2 * (label_distance * pie_r + text_w + margin) + legend_w

    # 垂直：最坏情况约 n/2 个标签堆在上半或下半
    line_h     = label_font_size * 1.6 / pt_per_inch
    title_h    = (label_font_size + 4) * 2 / pt_per_inch if title else 0.0
    h_labels   = max(1, n_direct_labels // 2 + 1) * line_h / 0.65
    h          = max(w * 0.85, h_labels + title_h + 1.0, 5.0)
    w          = max(w, 5.0)

    return round(w, 2), round(h, 2)


def plot_pie(
    values: list,
    labels: list,
    title: str          = '',
    colors: list        = None,
    explode: list       = None,
    donut: bool         = False,
    donut_ratio: float  = 0.55,
    show_pct: bool      = True,
    show_count: bool    = False,
    pct_distance: float = 0.75,
    label_distance: float = 1.18,
    startangle: float   = 90,
    label_font_size: float = 10,   # pt，外部标签字体大小
    pct_font_size: float   = 8,    # pt，扇区内百分比字体大小
    title_font_size: float = 11,   # pt，标题字体大小
    legend_font_size: float = 9,   # pt，图例字体大小
    small_threshold: float = 5.0,  # 小于此百分比的扇区改用图例，避免标签重叠
    figsize: tuple      = None,    # None = 自动推算
    out_path: str       = None,    # 支持 .png / .pdf / .svg；自动同时保存三种格式
):
    """
    绘制饼图或环形图，支持自动 figsize 和固定 pt 字体，输出 PNG + PDF + SVG。

    Parameters
    ----------
    values            : 各类别数值列表
    labels            : 各类别标签列表
    title             : 图标题
    colors            : 颜色列表，默认使用内置配色
    explode           : 各扇区偏移量，如 [0.08, 0, 0]
    donut             : True = 环形图
    donut_ratio       : 环形内圈半径比例（0~1）
    show_pct          : 扇区内显示百分比
    show_count        : 百分比后附加 n=xxx
    pct_distance      : 百分比距圆心比例
    label_distance    : 外部标签距圆心比例
    label_font_size   : 外部标签字体（pt）
    pct_font_size     : 扇区内百分比字体（pt）
    title_font_size   : 标题字体（pt）
    legend_font_size  : 图例字体（pt）
    small_threshold   : 扇区占比小于此值（%）时改用图例（避免标签重叠）
    figsize           : 图像尺寸，None 自动推算
    out_path          : 输出路径（含扩展名），自动同时保存 PNG / PDF / SVG
    """
    assert len(values) == len(labels), "values 和 labels 长度必须一致"

    values = np.array(values, dtype=float)
    n      = len(values)
    total  = values.sum()

    if colors is None:
        colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(n)]
    if explode is None:
        explode = [0.0] * n

    # ── 区分直接标签 vs 图例标签 ─────────────────────────────────────
    pcts          = values / total * 100
    direct_mask   = pcts >= small_threshold          # True = 显示外部标签
    legend_mask   = ~direct_mask

    display_labels = [l if direct_mask[i] else '' for i, l in enumerate(labels)]
    has_legend     = legend_mask.any()

    # ── 自动 figsize ─────────────────────────────────────────────────
    if figsize is None:
        direct_labels_list = [l for i, l in enumerate(labels) if direct_mask[i]]
        max_chars = max((len(l) for l in direct_labels_list), default=8)
        figsize = _auto_figsize(
            n_direct_labels=int(direct_mask.sum()),
            max_label_chars=max_chars,
            label_font_size=label_font_size,
            label_distance=label_distance,
            title=title,
            has_legend=has_legend,
        )

    # ── 百分比格式 ───────────────────────────────────────────────────
    def autopct_func(pct):
        if not show_pct:
            return ''
        if show_count:
            count = int(round(pct / 100.0 * total))
            return f'{pct:.1f}%\n(n={count})'
        return f'{pct:.1f}%'

    # ── 绘图 ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=figsize)

    wedges, texts, autotexts = ax.pie(
        values,
        labels=display_labels,
        colors=colors,
        explode=explode,
        autopct=autopct_func if show_pct else None,
        pctdistance=pct_distance,
        labeldistance=label_distance,
        startangle=startangle,
        wedgeprops=dict(linewidth=1.5, edgecolor='white'),
        textprops=dict(fontsize=label_font_size, fontfamily='Arial'),
    )

    # 外部标签
    for text in texts:
        text.set_fontsize(label_font_size)
        text.set_fontweight('bold')
        text.set_fontfamily('Arial')

    # 内部百分比
    for autotext in autotexts:
        autotext.set_fontsize(pct_font_size)
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontfamily('Arial')

    # 环形图
    if donut:
        centre_circle = plt.Circle((0, 0), donut_ratio, color='white', zorder=10)
        ax.add_artist(centre_circle)

    # 标题
    if title:
        ax.set_title(title, fontsize=title_font_size, fontweight='bold',
                     fontfamily='Arial', pad=16)

    ax.set_aspect('equal')

    # ── 小扇区图例 ───────────────────────────────────────────────────
    if has_legend:
        legend_handles = []
        legend_texts   = []
        for i, (w, l) in enumerate(zip(wedges, labels)):
            if legend_mask[i]:
                count_str = f' (n={int(values[i])})' if show_count else ''
                legend_texts.append(f'{l}  {pcts[i]:.1f}%{count_str}')
                legend_handles.append(
                    plt.Rectangle((0, 0), 1, 1, fc=colors[i], ec='white', lw=1)
                )
        ax.legend(
            legend_handles, legend_texts,
            loc='lower right',
            bbox_to_anchor=(1.35, 0.0),
            fontsize=legend_font_size,
            frameon=True,
            framealpha=0.85,
            edgecolor='#CCCCCC',
            title='< 5%',
            title_fontsize=legend_font_size,
            prop={'family': 'Arial', 'size': legend_font_size},
        )

    # ── 保存 ─────────────────────────────────────────────────────────
    if out_path:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        for suffix in ['.png', '.pdf', '.svg']:
            p = out_path.with_suffix(suffix)
            dpi = 300 if suffix == '.png' else None
            fig.savefig(p, dpi=dpi, bbox_inches='tight', facecolor='white')
            print(f"Saved: {p}")

    return fig, ax


def nested_pie(
    groups,
    *,
    ax=None,
    figsize: tuple        = (7.0, 7.0),
    inner_radius: float   = 0.42,    # 中心空洞半径(donut hole)
    inner_width: float    = 0.26,    # 内环厚度
    outer_width: float    = 0.30,    # 外环厚度
    gap: float            = 0.0,     # 内外环之间留白
    startangle: float     = 90,
    counterclock: bool    = False,
    edgecolor: str        = 'white',
    edgewidth: float      = 0.6,
    inner_label_min_pct: float = 4.0,   # 内环扇区 ≥ 此占比才标 label
    inner_label_size: float    = 7.5,
    inner_label_color: str     = '#2b2b2b',
    center_text: str      = '',
    center_text_size: float = 11,
):
    """
    两层嵌套环图(sunburst 风格,纯 matplotlib —— 静态导出可控,不依赖 plotly/kaleido)。

    内环 = `groups` 的各大类;外环 = 每个 group 的 `children`,在父扇区角度内细分。
    内外环角度自动对齐(每个 group 的子项数值之和即父值),因此外环每段都落在它所属
    大类的扇区下方,等价于一层 drill-down 的 sunburst。

    Parameters
    ----------
    groups : list[dict]
        有序大类列表,每项:
            {'label': str,                       # 大类名(内环 label)
             'color': hex,                       # 内环扇区颜色
             'value': num,                       # 可省略,默认 = 子项数值之和
             'children': [{'label':str, 'value':num, 'color':hex}, ...]}
        children 为空时该大类外环留一段同色实环。
    inner_radius, inner_width, outer_width, gap : float
        半径 / 环厚控制(数据坐标,饼图坐标系单位 ≈1)。
    inner_label_min_pct : float
        内环扇区占总数比例 ≥ 此值才在环带中点标大类名(避免薄扇区文字挤叠)。
    center_text : str
        中心空洞文字(常用于放总数,如 'n = 2,480')。

    Returns
    -------
    dict : {'fig', 'ax', 'inner_wedges', 'outer_wedges',
            'inner_values', 'outer_segments'}
        outer_segments = [(group_label, child_label, value, color), ...] 与 outer_wedges 对齐。

    Notes
    -----
    · 外环不直接标 label(段数多易糊),调用方用 children 的 color 配一份侧边 legend。
    · 图例 / 标题 / 中心额外注释由调用方在返回的 ax 上叠加,本函数只负责双环本体。
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    inner_values, inner_colors, inner_labels = [], [], []
    outer_values, outer_colors, outer_segments = [], [], []
    for g in groups:
        children = g.get('children', []) or []
        val = g.get('value')
        if val is None:
            val = sum(c['value'] for c in children) if children else 0
        if val <= 0:
            continue
        inner_values.append(val)
        inner_colors.append(g['color'])
        inner_labels.append(g.get('label', ''))
        if children:
            for c in children:
                if c['value'] <= 0:
                    continue
                outer_values.append(c['value'])
                outer_colors.append(c['color'])
                outer_segments.append((g.get('label', ''), c.get('label', ''),
                                       c['value'], c['color']))
        else:
            # 无子项:外环用父颜色补一段,保持角度守恒
            outer_values.append(val)
            outer_colors.append(g['color'])
            outer_segments.append((g.get('label', ''), '', val, g['color']))

    total = float(sum(inner_values)) or 1.0
    wp = dict(edgecolor=edgecolor, linewidth=edgewidth)

    inner_wedges, _ = ax.pie(
        inner_values, colors=inner_colors,
        radius=inner_radius + inner_width,
        startangle=startangle, counterclock=counterclock,
        wedgeprops=dict(width=inner_width, **wp),
    )
    outer_wedges, _ = ax.pie(
        outer_values, colors=outer_colors,
        radius=inner_radius + inner_width + gap + outer_width,
        startangle=startangle, counterclock=counterclock,
        wedgeprops=dict(width=outer_width, **wp),
    )

    # 内环 label:环带中点,水平,仅够大的扇区
    import math as _m
    r_lab = inner_radius + inner_width / 2.0
    for w, lab, val in zip(inner_wedges, inner_labels, inner_values):
        if not lab or val / total * 100 < inner_label_min_pct:
            continue
        ang = _m.radians((w.theta1 + w.theta2) / 2.0)
        ax.text(r_lab * _m.cos(ang), r_lab * _m.sin(ang), lab,
                ha='center', va='center', rotation_mode='anchor',
                fontsize=inner_label_size, color=inner_label_color,
                fontweight='bold')

    if center_text:
        ax.text(0, 0, center_text, ha='center', va='center',
                fontsize=center_text_size, fontweight='bold', color='#2b2b2b')

    ax.set_aspect('equal')
    return dict(fig=fig, ax=ax, inner_wedges=inner_wedges,
                outer_wedges=outer_wedges, inner_values=inner_values,
                outer_segments=outer_segments)
