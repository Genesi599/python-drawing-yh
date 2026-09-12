"""Reusable icon assets for scientific figures."""
from __future__ import annotations

import re
from pathlib import Path

ICON_DIR = Path(__file__).parent / "lib"

_NUM = r"[-+]?(?:\d*\.\d+|\d+\.?)"
_TRANSLATE = re.compile(rf"translate\(\s*({_NUM})[ ,]+({_NUM})\s*\)")
_SCALE = re.compile(rf"scale\(\s*({_NUM})(?:[ ,]+({_NUM}))?\s*\)")


def tint_items(items, target: str, lo: float = 0.85, hi: float = 1.0):
    """按视觉主色染: 同色 path 面积合并, 最大色组 = 纯目标色, 其余按亮度梯度压暗。

    保留多 path 分层纹理(muscle 红层+暗纹 / lymphnode 深紫管+淡紫球), 平染会糊掉;
    "最亮层=纯色"是错口径(lym 淡紫球小、深紫管占一半墨水, 染出来主色发暗)。
    """
    def _lum(c):
        if not re.fullmatch(r"#[0-9A-Fa-f]{6}", c):
            return 0.0
        return 0.299*int(c[1:3], 16) + 0.587*int(c[3:5], 16) + 0.114*int(c[5:7], 16)
    def _area(p):
        # 渲染像素面积(shoelace 对带洞 path 失效; 洞与轮廓环绕方向相反相互抵消)
        import numpy as np
        fig = plt.figure(figsize=(0.6, 0.6), dpi=100)
        ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
        ax.add_patch(PathPatch(p, facecolor="k", edgecolor="none"))
        fig.canvas.draw()
        arr = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].astype(int)
        a = int((arr.sum(axis=2) < 300).sum())
        plt.close(fig)
        return float(a)
    if "_plt" not in globals():
        import matplotlib.pyplot as plt
    from matplotlib.patches import PathPatch

    # 同色分组合并面积 → 主色组
    from collections import defaultdict
    by_col = defaultdict(float)
    for p, c in items:
        by_col[c] += _area(p)
    main_col = max(by_col.items(), key=lambda kv: kv[1])[0] if by_col else None
    main_lum = _lum(main_col) if main_col else 0.0
    # 主色组与最亮组的亮度差 → 档位梯度(lo..hi); 单色图标 → 全部纯色
    cols = list(by_col)
    if len(cols) <= 1:
        lums = {c: 1.0 for c in cols}
    else:
        lo_l = min(_lum(c) for c in cols)
        hi_l = max(_lum(c) for c in cols)
        span = (hi_l - lo_l) or 1.0
        lums = {c: 0.0 if c == main_col else (_lum(c) - lo_l) / span for c in cols}
    tr, tg, tb = int(target[1:3], 16), int(target[3:5], 16), int(target[5:7], 16)
    out = []
    for p, c in items:
        n = lums.get(c, 0.0)
        f = 1.0 if c == main_col else lo + (hi - lo) * n
        out.append((p, (tr*f/255, tg*f/255, tb*f/255)))
    return out


def icon_path(name: str, ext: str = "png") -> Path:
    """库内资产路径; 不存在返回 None。"""
    p = ICON_DIR / f"{name}.{ext}"
    return p if p.exists() else None


def load_icon_vector(stem: str):
    """库 SVG → ([(matplotlib.Path, fill), ...], aspect) 或 None。

    每个 <path> 独立保留自己的 fill 颜色(多色分层结构完整, 2026-09-12 教训:
    合并单色填充会把多 path 图标糊成一团)。组 transform(translate/scale)已应用。
    坐标归一到单位空间: y∈[0,1](已翻转为 matplotlib y-up), x∈[0,aspect]。
    """
    p = icon_path(stem, "svg")
    if p is None:
        return None
    src = p.read_text(encoding="utf-8")
    from svgpath2mpl import parse_path  # 延迟导入, 无依赖环境不受影响

    # 组 transform(库内最多一层 <g>, 作用于全部 path)
    tr = re.search(r"<g\b[^>]*\btransform\s*=\s*\"([^\"]*)\"", src)
    tx = ty = 0.0
    sx = sy = 1.0
    if tr:
        m = _TRANSLATE.search(tr.group(1))
        if m:
            tx, ty = float(m.group(1)), float(m.group(2))
        m = _SCALE.search(tr.group(1))
        if m:
            sx = float(m.group(1))
            sy = float(m.group(2)) if m.group(2) else sx

    import numpy as np
    from matplotlib.path import Path as MplPath

    raw = []  # [(verts_viewBox, codes, fill)]
    for m in re.finditer(r"<path\b([^>]*)>", src):
        attrs = m.group(1)
        d = re.search(r'\bd\s*=\s*"([^"]*)"', attrs)
        if not d or not d.group(1).strip():
            continue
        fill = re.search(r'\bfill\s*=\s*"([^"]*)"', attrs)
        color = fill.group(1) if fill else "#000000"
        if color in ("none", "transparent"):
            continue
        if not re.fullmatch(r"#[0-9A-Fa-f]{3,8}", color):
            color = "#333333"  # currentColor/命名色/rgb() 等 → 统一深灰(与旧 PNG 版一致)
        try:
            mp = parse_path(d.group(1))
        except Exception:
            continue
        v = mp.vertices.astype(float)
        v = np.column_stack([v[:, 0] * sx + tx, v[:, 1] * sy + ty])
        raw.append((v, mp.codes, color))
    if not raw:
        return None

    # union bbox → fit-to-square: 以最大边为 1, 短边方向两侧居中
    # (2026-09-12: 按高归一会让宽扁图标(liver)与细长图标(spleen)视觉大小差 5 倍;
    #  x 左对齐会让左缘参差, 与旧 PNG 方形画布行为不一致)
    import numpy as np
    allv = np.vstack([v for v, _, _ in raw])
    x0, x1 = allv[:, 0].min(), allv[:, 0].max()
    y0, y1 = allv[:, 1].min(), allv[:, 1].max()
    H = y1 - y0
    Wd = x1 - x0
    if H <= 0 or Wd <= 0:
        return None
    S = max(H, Wd)          # 最大边 → 1
    # 居中位移在"归一化后的单位空间"里算((S-Wd)/2 是 viewBox 原始单位, 直接加会
    # 把图标平移几十个单位 — lymphnode x 偏 87 / liver y 偏 83 的根源)
    x_shift = (1.0 - Wd / S) / 2
    y_shift = (1.0 - H / S) / 2

    items = []
    for v, codes, color in raw:
        vn = np.column_stack([
            (v[:, 0] - x0) / S + x_shift,
            (y1 - v[:, 1]) / S + y_shift,
        ])
        items.append((MplPath(vn, codes), color))
    return items, 1.0  # 方形: DrawingArea 统一 ICON_H_PT 见方, 各行左/右缘共线
