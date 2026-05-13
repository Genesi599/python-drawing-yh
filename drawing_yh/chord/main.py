"""
drawing_yh.chord — directional chord diagram(基于 pycirclize)

适用场景:细胞-细胞通讯网络(sender → receiver,值 = 通讯强度 / event 数)、
任何 "from → to" 的有向加权图(基因调控、迁移流向、引文网络…)。

核心函数 `chord_diagram(matrix, ...)`:
    - matrix: pandas DataFrame,index = sender(源),columns = receiver(目标),值 = 权重
    - pycirclize 自动按 (行和 + 列和) 分配每个 node 的 sector 长度,ribbon 在 sector 上错开排布
      (标准 circlize chordDiagram 行为 —— 同一 node 的进 / 出 ribbon 不会堆在一点)
    - 返回 (fig, circos),fig 可直接交给 drawing_yh.save_fig 写出

附带 `HEMATOPOIETIC_LINEAGE_COLORS` —— 造血谱系细胞类型的成系配色 dict(BM 项目复用):
    HSPC/progenitor=teal,B/Plasma=blue,DC=purple,T/NK=olive,
    Neutrophil=red,Mono/Mac=orange,Mk=magenta,Erythroid=brown

依赖:pycirclize(已加进 pyproject.toml)
"""
from __future__ import annotations

import matplotlib.colors as _mc
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
    space: int = 3,
    r_lim: tuple = (93, 100),
    alpha: float = 0.42,
    ec: str = "white",
    lw: float = 0.15,
    direction: int = 1,
    label_size: int = 9,
    label_orientation: str = "vertical",
    figsize: tuple = (8, 8),
    title: str | None = None,
    title_size: int = 11,
    drop_zero_nodes: bool = True,
):
    """画 directional chord diagram。

    Params
    ------
    matrix
        pandas DataFrame,index = sender,columns = receiver,值 = 权重(>=0)。
        非对称即有向(direction=1 时 row → col)。
    color_map
        dict {node_name: color}。缺省按 HEMATOPOIETIC_LINEAGE_COLORS 找,
        找不到的 node 按 _FALLBACK_PALETTE 循环取色。
        显式传 {} 则全部用 fallback。
    space
        sector 之间的间隔角度。node 多时调小(如 2),少时可大(如 4)。
    r_lim
        sector 轨道的半径范围(0-100)。
    alpha, ec, lw
        ribbon 的透明度 / 边线颜色 / 边线宽度。
    direction
        1 = ribbon 带方向(sender 一端宽,receiver 一端窄;或反之,按 pycirclize 约定);
        0 = 无方向。
    label_size, label_orientation
        node 标签字号 / 方向("vertical" / "horizontal")。
    figsize
        figure 尺寸(英寸)。
    title, title_size
        总标题(可多行,用 \\n)及字号。
    drop_zero_nodes
        True 时去掉行和 + 列和都为 0 的 node(不参与通讯的细胞类型)。

    Returns
    -------
    (fig, circos)
        matplotlib Figure 与 pycirclize Circos 对象。
        写出建议 `from drawing_yh import save_fig; save_fig(fig, 'chord.pdf', also=('.png',))`。

    Example
    -------
        import pandas as pd
        from drawing_yh.chord import chord_diagram, HEMATOPOIETIC_LINEAGE_COLORS
        from drawing_yh import save_fig

        # mat: index=sender cell type, columns=receiver cell type, values=通讯强度
        fig, _ = chord_diagram(mat, title='Urea-mediated mCCC')
        save_fig(fig, 'out/chord_urea.pdf', also=('.png',))
    """
    try:
        from pycirclize import Circos
    except ImportError as e:  # 边界:外部依赖缺失,给清晰提示
        raise ImportError(
            "chord_diagram 需要 pycirclize:pip install pycirclize"
        ) from e
    import matplotlib.pyplot as plt

    mat = matrix.copy()
    # 对齐 index / columns 的并集(否则 pycirclize 报错)
    all_nodes = sorted(set(mat.index) | set(mat.columns))
    mat = mat.reindex(index=all_nodes, columns=all_nodes, fill_value=0.0)

    if drop_zero_nodes:
        keep = [n for n in all_nodes if mat.loc[n].sum() > 0 or mat[n].sum() > 0]
        mat = mat.loc[keep, keep]
        all_nodes = keep

    # 颜色
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

    circos = Circos.chord_diagram(
        mat,
        space=space,
        cmap=cmap,
        r_lim=r_lim,
        ticks_interval=None,
        label_kws=dict(size=label_size, orientation=label_orientation),
        link_kws=dict(direction=direction, alpha=alpha, ec=ec, lw=lw),
    )
    fig = circos.plotfig(figsize=figsize)
    if title:
        fig.suptitle(title, fontsize=title_size, y=1.02)
    return fig, circos
