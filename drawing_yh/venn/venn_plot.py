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
- **figsize 不预先固定**:默认 ``autoshrink=True`` 用 ``autoshrink_figsize`` 按 label
  bbox 重叠实测,找不重叠的最小尺寸(符合"版面尽量紧凑"标准)。
- **图内 title 求简**:``title`` 直接画在图上,**保持简短**;图型前缀 + 视觉编码写进报告 figcaption,不塞进图(见 preferences 出图规范)。
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

from ..layout import ONE_HALF_COL_IN, autoshrink_figsize
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
                 members_per_line=3, max_members=None, member_order=None, title=None,
                 ax=None, figsize=None, autoshrink=True,
                 set_label_size=None, subset_label_size=None,
                 edgecolor='white', linewidth=0.6):
    """画一张 2 或 3 集合 Venn 图,返回 ``(fig, ax)``。

    参数
    ----
    sets : list[set] —— 2 或 3 个集合(原始元素,用于算真实区域计数)。
    labels : list[str] —— 集合名(顺序与 sets 对应;通常自带 ``(size)``)。
    colors : 可选,默认 ``OKABE_ITO[:n]``。
    mode : ``'proportional'`` | ``'log'`` | ``'equal'``,见模块 docstring。
    title : 图内标题,**保持简短**(图型前缀 + 视觉编码写进报告 figcaption,不塞进图)。
    autoshrink : ``True`` 自动找不重叠的最小 figsize(``figsize`` 给定时忽略)。
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

    full_title = title

    draw = venn2 if n == 2 else venn3
    aspect = 1.2 if want_members else 1.0

    def _draw(axv):
        v = draw(_subset_arg(sets, mode), set_labels=labels,
                 set_colors=tuple(colors), alpha=alpha, ax=axv)
        # 区域 label 一律覆盖成真实计数(0 区域留空)
        for rid in ids:
            lbl = v.get_label_by_id(rid)
            if lbl is not None:
                c = counts[rid]
                lbl.set_text(str(c) if c else '')
                lbl.set_fontsize(fs_sub)
        for t in v.set_labels:
            if t is not None:
                t.set_fontsize(fs_set)
        for rid in ids:
            p = v.get_patch_by_id(rid)
            if p is not None:
                p.set_edgecolor(edgecolor)
                p.set_linewidth(linewidth)
        if n == 3:
            # 默认 2 上 1 下 → 翻转成 1 上 2 下:底部中间空出给箭头 / 成员列表,
            # 单圆标题随之移到顶部,不再与向下箭头重叠(文字仍正立)
            axv.invert_yaxis()
        if full_title:
            axv.set_title(full_title)
        return v

    def _add_members(fig, axv, v):
        # 用 fig.text(非独立 axes)列成员,避免空 axes 被 tight bbox 计入留白
        # member_order 给定 → 按该顺序(如显著性)排;否则按名字字母序。
        if member_order:
            mo = list(dict.fromkeys(member_order))
            ordered = [x for x in mo if x in inter] + sorted(x for x in inter if x not in set(mo))
        else:
            ordered = sorted(inter, key=lambda x: str((member_labels or {}).get(x, x)))
        names = [str((member_labels or {}).get(x, x)) for x in ordered]
        # max_members 给定且超出 → 只列前 N,余者用 "… (+k more)" 省略
        extra = 0
        if max_members is not None and len(names) > max_members:
            extra = len(names) - max_members
            names = names[:max_members]
        head = member_title or f"{n}-set intersection ({len(inter)})"
        body = '\n'.join(', '.join(names[i:i + members_per_line])
                         for i in range(0, len(names), members_per_line))
        if extra:
            body += f'\n… (+{extra} more)'
        center = v.get_patch_by_id('111' if n == 3 else '11')
        if center is not None:
            verts = center.get_path().vertices
            cx, cy = verts[:, 0].mean(), verts[:, 1].mean()
            axv.add_patch(FancyArrowPatch(
                (cx, cy), (cx, axv.get_ylim()[0]), arrowstyle='-|>',
                mutation_scale=7, lw=0.7, color='0.35',
                connectionstyle='arc3,rad=0.08', zorder=60, clip_on=False))
        fig.text(0.5, 0.125, head, ha='center', va='top', fontweight='bold',
                 color='#1A5490', fontsize=fs_sub)
        fig.text(0.5, 0.075, body, ha='center', va='top', fontsize=fs_sub)

    # 调用方给了 ax:直接画进去(成员列表需自建图,这里跳过)
    if ax is not None:
        _draw(ax)
        return ax.figure, ax

    def _build(w, h):
        fig = plt.figure(figsize=(w, h))
        if want_members:
            axv = fig.add_axes([0.04, 0.20, 0.92, 0.68])
            v = _draw(axv)
            _add_members(fig, axv, v)
        else:
            axv = fig.add_axes([0.03, 0.03, 0.94, 0.86])
            _draw(axv)
        return fig, axv

    # figsize 显式给定 → 直接用;否则默认 autoshrink 找最小不重叠尺寸
    if figsize is not None:
        w, h = figsize if isinstance(figsize, (tuple, list)) else (figsize, figsize * aspect)
        return _build(w, h)
    if autoshrink:
        return autoshrink_figsize(
            lambda s: _build(s, s * aspect),
            initial=4.5, min_size=3.0, max_size=7.5,
            get_texts=lambda a: [t for t in (list(a.texts) + [a.title]) if t.get_text()],
            bbox_shrink=0.6,
        )
    return _build(ONE_HALF_COL_IN, ONE_HALF_COL_IN * aspect)
