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


def chord_diagram(
    matrix: pd.DataFrame,
    color_map: dict | None = None,
    *,
    figsize: tuple = (2.5, 2.5),    # 起始 figsize,autoshrink 会从这缩到 labels 刚好不重叠
    fontsize: int = 9,              # 固定 fontsize(不自动缩,由 figsize 调整防重叠)
    alpha: float = 1.0,             # chord 不透明,深
    pad: float = 3.0,               # sector 间空隙(度)
    chordwidth: float = 0.5,        # chord 弯曲度(0.3 直,0.5 适中,0.7 中段窄)
    intra_gap: float = 1.5,         # 同 sender 内 chord 子区段空隙(度)
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
    all_nodes = sorted(set(mat.index) | set(mat.columns))
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

    def _check_overlap(_ax, _fig):
        """Tight overlap check:bbox 各收缩到中心 50% 才算重叠,允许 labels 物理接近"""
        from matplotlib.transforms import Bbox as _Bbox
        _fig.canvas.draw()
        _t = list(_ax.texts)
        if len(_t) < 2:
            return False
        _renderer = _fig.canvas.get_renderer()
        def _shrunk(_b, _f=0.5):
            cx, cy = (_b.x0 + _b.x1) / 2, (_b.y0 + _b.y1) / 2
            hw, hh = _b.width * _f / 2, _b.height * _f / 2
            return _Bbox.from_extents(cx - hw, cy - hh, cx + hw, cy + hh)
        _bb = [_shrunk(t.get_window_extent(_renderer)) for t in _t]
        for ii in range(len(_bb)):
            for jj in range(ii + 1, len(_bb)):
                if _bb[ii].overlaps(_bb[jj]):
                    return True
        return False

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
                   intra_gap=intra_gap,
                   rotate_names=radial_labels,
                   show=False,
                   fontsize=fontsize,
                   alpha=eff_alpha)
        return _fig, _ax

    # Autoshrink figsize:从 figsize 开始往下缩到 labels 临界不重叠
    # 下限 min_figsize:sectors 少时确保 chord 不被 labels 压扁(2.2 base + 每 5 sector +0.2)
    n_sec = len(all_nodes)
    min_figsize = max(2.2, 2.2 + (n_sec - 6) * 0.04)
    fig_sz = figsize[0]
    fig, ax = _render(fig_sz)
    last_good_sz = fig_sz
    if not _check_overlap(ax, fig):
        # 不 overlap,试更小(直到 min_figsize 或 overlap)
        while fig_sz > min_figsize:
            new_sz = max(fig_sz * 0.85, min_figsize)
            plt.close(fig)
            fig, ax = _render(new_sz)
            if _check_overlap(ax, fig):
                # 重叠了 → 回退用 last_good_sz
                plt.close(fig)
                fig, ax = _render(last_good_sz)
                break
            last_good_sz = fig_sz = new_sz
            if fig_sz <= min_figsize:
                break
    else:
        # 初始 figsize 就重叠 → 试更大
        while fig_sz < 8.0:
            new_sz = fig_sz * 1.15
            plt.close(fig)
            fig, ax = _render(new_sz)
            if not _check_overlap(ax, fig):
                break
            fig_sz = new_sz

    return fig, ax
