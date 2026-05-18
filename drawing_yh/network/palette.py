"""
模块化网络的配色 + hex 颜色 lighten/darken。

色板按"rank 1+2 不同时是红+绿"挑过(色盲友好);超出 10 个模块的部分按
golden-ratio hue 步长 + 低饱和 pastel 循环。
"""
from __future__ import annotations

import colorsys
from collections import Counter


MODULE_PALETTE: list[str] = [
    "#4E79A7",   # 0 steel blue
    "#F28E2B",   # 1 warm orange
    "#76B7B2",   # 2 teal
    "#E15759",   # 3 coral red
    "#B07AA1",   # 4 plum
    "#59A14F",   # 5 leaf green
    "#EDC948",   # 6 mustard
    "#9C755F",   # 7 umber
    "#FF9DA7",   # 8 pink
    "#4A6FA5",   # 9 ocean blue
]


def module_palette(n_modules: int, sizes_counter: Counter) -> list[str]:
    """按模块大小 rank 分配:前 10 名走 MODULE_PALETTE,其余 golden-ratio pastel 循环。
    返回长度为 n_modules 的 hex 颜色列表,index = module_id。"""
    rank_order = [m for m, _ in sizes_counter.most_common()]
    palette = ["#cbd5e1"] * n_modules
    for rank, m in enumerate(rank_order):
        if m >= n_modules:
            continue
        if rank < len(MODULE_PALETTE):
            palette[m] = MODULE_PALETTE[rank]
        else:
            hue = (rank * 0.61803398875) % 1.0
            r, g, b = colorsys.hsv_to_rgb(hue, 0.35, 0.78)
            palette[m] = f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
    return palette


def lighten_hex(hex_color: str, amount: float = 0.5) -> str:
    """向白色混合 amount∈[0,1] 比例;输入输出都是 '#rrggbb'。
    drawing_yh.chord.lighten 返回 RGB tuple,这里返回 hex 字符串。"""
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    r = int(r + (255 - r) * amount)
    g = int(g + (255 - g) * amount)
    b = int(b + (255 - b) * amount)
    return f"#{r:02x}{g:02x}{b:02x}"


def darken_hex(hex_color: str, amount: float = 0.3) -> str:
    """向黑色混合 amount∈[0,1] 比例;输入输出都是 '#rrggbb'。"""
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    r = int(r * (1 - amount))
    g = int(g * (1 - amount))
    b = int(b * (1 - amount))
    return f"#{r:02x}{g:02x}{b:02x}"
