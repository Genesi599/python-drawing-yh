# -*- coding: utf-8 -*-
"""通用 2/3 集合 Venn 图(符合 drawing-yh 出图标准)。

特性
----
- ``import drawing_yh`` 已设 8pt / Arial / fonttype42 rc;本函数所有文字默认 8pt。
- 配色默认 ``OKABE_ITO``(避开 matplotlib_venn 默认红/绿,色盲友好)。
- 三种圆面积模式(``mode``):
    * ``'proportional'`` 真比例(matplotlib_venn 默认,面积 = 真实子集大小)
    * ``'log'``          log 加权(集合大小悬殊时小圆不消失;面积示意、非真比例)
    * ``'equal'``        等大 schematic(所有圆等大,只靠数字标注)
  无论哪种,**每个区域 label 始终标真实计数**(0 区域留空)。
- ``show_members``:把全集交集成员列在图下方 + 箭头(成员可经 ``member_labels``
  映射成可读名,如 HMDB ID → 代谢物名)。需 ``ax=None`` 由本函数建图。
- 返回 ``(fig, ax)``,出图由调用方 ``save_fig`` 写三格式(不内嵌 savefig)。

用法
----
    import drawing_yh
    from drawing_yh import venn_diagram, save_fig

    fig, ax = venn_diagram(
        [sc, detected, sig],
        [f'sc ({len(sc)})', f'detected ({len(detected)})', f'sig ({len(sig)})'],
        mode='log', show_members=True, member_labels=hmdb2name,
        title='Venn diagram — ...')
    save_fig(fig, 'out/venn.pdf', also=('.png', '.svg'))
"""
from __future__ import annotations

import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

from .. import ONE_HALF_COL_IN
from ..palettes import OKABE_ITO

# matplotlib_venn 的 subset id 顺序
_VENN3_IDS = ('100', '010', '110', '001', '101', '011', '111')
_VENN2_IDS = ('10', '01', '11')


def _region_counts(sets: list[set]) -> dict[str, int]:
    """各区域的真实元素数,key 为 matplotlib_venn 的 subset id。"""
    if len(sets) == 2:
        a, b = sets
        return {'10': len(a - b), '01': len(b - a), '11': len(a & b)}
    a, b, c = sets
    return {
        '100': len(a - b - c), '010': len(b - a - c), '001': len(c - a - b),
        '110': len((a & b) - c), '101': len((a & c) - b), '011': len((b & c) - a),
        '111': len(a & b & c),
    }


def _subset_arg(sets: list[set], mode: str):
    """返回传给 venn2/venn3 的 subsets 实参(决定圆/区域几何大小)。"""
    if mode == 'proportional':
        return list(sets)                                  # 真比例:直接传集合
    if mode == 'equal':
        return tuple(1.0 for _ in (_VENN2_IDS if len(sets) == 2 else _VENN3_IDS))
    if mode == 'log':
        # 圆面积 ~ log10(set size);overlap 区给小常数,保证可见但不喧宾夺主
        w = [math.log10(len(s) + 1) + 1e-3 for s in sets]
        ov = max(min(w) * 0.30, 1e-3)
        if len(sets) == 2:
            return (w[0], w[1], ov)                        # (10, 01, 11)
        return (w[0], w[1], ov, w[2], ov, ov, ov)          # (100,010,110,001,101,011,111)
    raise ValueError(f"mode must be 'proportional'|'log'|'equal', got {mode!r}")


def venn_diagram(sets, labels, *, colors=None, mode='proportional', alpha=0.5,
                 show_members=False, member_labels=None, member_title=None,
                 members_per_line=3, title=None, ax=None, figsize=None,
                 set_label_size=None, subset_label_size=None,
                 edgecolor='white', linewidth=0.6):
    """画一张 2 或 3 集合 Venn 图,返回 ``(fig, ax)``。

    参数
    ----
    sets : list[set] —— 2 或 3 个集合(原始元素,用于算真实区域计数)。
    labels : list[str] —— 集合名(顺序与 sets 对应;通常自带 ``(size)``)。
    colors : 可选,默认 ``OKABE_ITO[:n]``。
    mode : ``'proportional'`` | ``'log'`` | ``'equal'``,见模块 docstring。
    show_members : 在图下方列全集交集成员(需 ``ax=None``)。
    member_labels : dict,把交集元素映射成可读名(如 {HMDB: 代谢物名})。
    member_title : 成员块标题,默认 ``f"{n}-set intersection ({k})"``。
    """
    sets = [set(s) for s in sets]
    n = len(sets)
    if n not in (2, 3):
        raise ValueError(f"venn_diagram 只支持 2 或 3 个集合,收到 {n}")
    if len(labels) != n:
        raise ValueError("labels 长度必须与 sets 一致")

    from matplotlib_venn import venn2, venn3   # 延迟导入(可选依赖,缺了不影响 import drawing_yh)

    if colors is None:
        colors = OKABE_ITO[:n]
    fs_set = set_label_size if set_label_size is not None else mpl.rcParams['font.size']
    fs_sub = subset_label_size if subset_label_size is not None else mpl.rcParams['font.size']

    counts = _region_counts(sets)
    ids = _VENN2_IDS if n == 2 else _VENN3_IDS
    inter = set.intersection(*sets)
    want_members = show_members and ax is None and len(inter) > 0

    # ── figure / axes ──
    if ax is not None:
        fig = ax.figure
        ax_txt = None
    elif want_members:
        fig = plt.figure(figsize=figsize or (ONE_HALF_COL_IN, ONE_HALF_COL_IN * 1.35))
        ax = fig.add_axes([0.04, 0.30, 0.92, 0.66])
        ax_txt = fig.add_axes([0.04, 0.01, 0.92, 0.25])
        ax_txt.axis('off')
    else:
        fig, ax = plt.subplots(figsize=figsize or (ONE_HALF_COL_IN, ONE_HALF_COL_IN))
        ax_txt = None

    draw = venn2 if n == 2 else venn3
    v = draw(_subset_arg(sets, mode), set_labels=labels,
             set_colors=tuple(colors), alpha=alpha, ax=ax)

    # 区域 label 一律覆盖成真实计数(0 区域留空)
    for rid in ids:
        lbl = v.get_label_by_id(rid)
        if lbl is not None:
            c = counts[rid]
            lbl.set_text(str(c) if c else '')
            lbl.set_fontsize(fs_sub)
    # 集合名字号 + 区域描边
    for t in v.set_labels:
        if t is not None:
            t.set_fontsize(fs_set)
    for rid in ids:
        p = v.get_patch_by_id(rid)
        if p is not None:
            p.set_edgecolor(edgecolor)
            p.set_linewidth(linewidth)

    if title:
        ax.set_title(title)

    # 成员列表 + 箭头(原 tuned 特性)
    if want_members:
        names = sorted(str((member_labels or {}).get(x, x)) for x in inter)
        head = member_title or f"{n}-set intersection ({len(names)})"
        body = '\n'.join(', '.join(names[i:i + members_per_line])
                         for i in range(0, len(names), members_per_line))
        center = v.get_patch_by_id('111' if n == 3 else '11')
        if center is not None:
            verts = center.get_path().vertices
            cx, cy = verts[:, 0].mean(), verts[:, 1].mean()
            y0 = ax.get_ylim()[0]
            ax.add_patch(FancyArrowPatch(
                (cx, cy), (cx, y0), arrowstyle='-|>', mutation_scale=7,
                lw=0.7, color='0.35', connectionstyle='arc3,rad=0.08',
                zorder=60, clip_on=False))
        ax_txt.text(0.5, 0.95, head, ha='center', va='top', fontweight='bold',
                    color='#1A5490', transform=ax_txt.transAxes, fontsize=fs_sub)
        ax_txt.text(0.5, 0.60, body, ha='center', va='top',
                    transform=ax_txt.transAxes, fontsize=fs_sub)

    return fig, ax
