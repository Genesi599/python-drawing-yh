"""
drawing_yh.layout — figsize 推算 + 文字防重叠工具

- 不预先固定 figsize:按字号 + 内容量动态算
- 文字防重叠靠 `fig.canvas.draw() + get_window_extent()` 读 pixel bbox,不靠经验估
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np


# ============================================================
# 期刊版面英寸换算
# ============================================================
SINGLE_COL_IN = 3.35       # 单栏 ~ 8.5 cm
ONE_HALF_COL_IN = 4.49     # 1.5 栏 ~ 11.4 cm
DOUBLE_COL_IN = 6.89       # 双栏 ~ 17.5 cm


def compute_figsize(
    n_items: int,
    *,
    width: float = DOUBLE_COL_IN,
    per_item_h: float = 0.25,
    base_h: float = 1.5,
    min_h: float = 2.0,
    max_h: float = 8.0,
    aspect: float = None,
) -> Tuple[float, float]:
    """
    按元素数量动态算 figsize,避免预先固定造成虚大 / 虚小。

    高度:`per_item_h * n_items + base_h`,clip 到 [min_h, max_h]。
    如果给了 `aspect`(h/w),按比例算高度,忽略 n_items 公式。

    Parameters
    ----------
    n_items : int
        Y 轴或主排列方向上的元素数(如基因数、样本数、条目数)。
    width : float
        英寸,默认双栏。可换 `SINGLE_COL_IN` / `ONE_HALF_COL_IN`。
    per_item_h, base_h, min_h, max_h : float
        线性公式 + clip 范围。默认值是 dot_chart 实战调出来的。
    aspect : float or None
        若给,直接 `height = width * aspect`,覆盖其它高度逻辑。
    """
    if aspect is not None:
        return (float(width), float(width) * float(aspect))
    h = float(np.clip(n_items * per_item_h + base_h, min_h, max_h))
    return (float(width), h)


# ============================================================
# 文字防重叠
# ============================================================

@dataclass
class TextBBox:
    """单个 Text 在 canvas 上的像素 bbox。"""
    text: object     # matplotlib.text.Text
    x0: float
    y0: float
    x1: float
    y1: float

    @property
    def width(self) -> float:
        return self.x1 - self.x0

    @property
    def height(self) -> float:
        return self.y1 - self.y0


def get_text_bboxes(fig, texts: Sequence) -> list:
    """
    取得每个 Text 的当前像素 bbox。
    内部先 `fig.canvas.draw()`,确保 renderer 状态最新。
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    out = []
    for t in texts:
        bb = t.get_window_extent(renderer=renderer)
        out.append(TextBBox(t, bb.x0, bb.y0, bb.x1, bb.y1))
    return out


def find_overlaps(bboxes: Sequence[TextBBox], gap_px: float = 0.0) -> list:
    """
    返回所有相互重叠的下标对 `[(i, j), ...]`。
    `gap_px`:认为是"重叠"的最小允许间距,> 0 即"靠太近也算"。
    """
    overlaps = []
    n = len(bboxes)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = bboxes[i], bboxes[j]
            if a.x1 + gap_px <= b.x0 or b.x1 + gap_px <= a.x0:
                continue
            if a.y1 + gap_px <= b.y0 or b.y1 + gap_px <= a.y0:
                continue
            overlaps.append((i, j))
    return overlaps


def nudge_no_overlap_y(fig, texts: Sequence, *,
                        step: float = 0.005,
                        max_iter: int = 200,
                        gap_px: float = 0.0) -> int:
    """
    沿 Y 方向(figure 相对坐标)逐步推开重叠文字。

    简单贪心:每轮取所有重叠对,把后者(j)往上推 `step`(figure 坐标系),
    重新 draw + 检测,直到无重叠或达到 `max_iter`。

    Returns
    -------
    int : 实际迭代次数(0 表示一开始就不重叠,`max_iter` 表示放弃)
    """
    for it in range(max_iter):
        bbs = get_text_bboxes(fig, texts)
        overlaps = find_overlaps(bbs, gap_px=gap_px)
        if not overlaps:
            return it
        for _, j in overlaps:
            t = texts[j]
            x, y = t.get_position()
            t.set_position((x, y + step))
    return max_iter


def _shrunk_bbox(bb, factor: float):
    """按中心收缩 bbox 到 factor 比例(0.5 = 缩到 50%)。"""
    from matplotlib.transforms import Bbox
    cx = (bb.x0 + bb.x1) / 2
    cy = (bb.y0 + bb.y1) / 2
    hw = bb.width * factor / 2
    hh = bb.height * factor / 2
    return Bbox.from_extents(cx - hw, cy - hh, cx + hw, cy + hh)


def autoshrink_figsize(
    render_fn,
    initial: float = 2.5,
    min_size: float = 1.8,
    max_size: float = 8.0,
    shrink: float = 0.85,
    grow: float = 1.15,
    bbox_shrink: float = 0.5,
    get_texts=None,
    max_iter: int = 30,
):
    """
    自动找最小 figsize:用 render_fn 在不同 figsize 上 render,检测 label bbox 重叠。
    不重叠就缩小,重叠就放大,直到收敛或触 min/max。

    Why
    ---
    出图标准是"按字号 + 内容量动态算出尽可能小的 figsize"(见 preferences.md / STANDARDS.md)。
    手动调 figsize 不准:sectors 多少、label 长短都影响。autoshrink 用实际 bbox overlap
    决定,跨图型通用。

    Parameters
    ----------
    render_fn : callable(figsize: float) -> (fig, ax)
        每次给个正方形 figsize,返回 (fig, ax)。允许 close 旧 fig。
    initial : float
        起始 figsize(inch)。
    min_size, max_size : float
        figsize 上下限。
    shrink, grow : float
        每次缩 / 放比例。
    bbox_shrink : float
        bbox 收缩比 — `0.5` = 各 label bbox 缩到中心 50% 才检测重叠(允许 labels
        物理接近,因为 matplotlib bbox 含 padding,实际 text 比 bbox 小)。
    get_texts : callable(ax) -> list[Text], 可选
        提取要检测的 text artists。默认 `ax.texts`(适合 chord-style sector labels;
        bar/scatter 之类需要传 `lambda ax: ax.get_xticklabels() + ax.get_yticklabels()`)。
    max_iter : int
        最多迭代次数,防止死循环。

    Returns
    -------
    (fig, ax)
        最终选定 figsize 上的 (fig, ax)。
    """
    import matplotlib.pyplot as _plt
    if get_texts is None:
        get_texts = lambda ax: list(ax.texts)

    def _has_overlap(_fig, _ax):
        _fig.canvas.draw()
        ts = get_texts(_ax)
        if len(ts) < 2:
            return False
        renderer = _fig.canvas.get_renderer()
        bbs = [_shrunk_bbox(t.get_window_extent(renderer), bbox_shrink) for t in ts]
        for i in range(len(bbs)):
            for j in range(i + 1, len(bbs)):
                if bbs[i].overlaps(bbs[j]):
                    return True
        return False

    fig_sz = initial
    fig, ax = render_fn(fig_sz)
    last_good = fig_sz
    iters = 0
    if not _has_overlap(fig, ax):
        # 不重叠,试更小
        while fig_sz > min_size and iters < max_iter:
            iters += 1
            new_sz = max(fig_sz * shrink, min_size)
            _plt.close(fig)
            fig, ax = render_fn(new_sz)
            if _has_overlap(fig, ax):
                _plt.close(fig)
                fig, ax = render_fn(last_good)
                break
            last_good = fig_sz = new_sz
            if fig_sz <= min_size:
                break
    else:
        # 重叠,试更大
        while fig_sz < max_size and iters < max_iter:
            iters += 1
            new_sz = fig_sz * grow
            _plt.close(fig)
            fig, ax = render_fn(new_sz)
            if not _has_overlap(fig, ax):
                break
            fig_sz = new_sz
    return fig, ax
