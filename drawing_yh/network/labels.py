"""
节点/模块标签的去重叠 + bbox 估算 + 节点级碰撞解开。

通用算法,不依赖具体网络数据;调用者负责喂 (x, y) + (width, height) 列表。
"""
from __future__ import annotations

import math
import numpy as np


def fs_px(s) -> float:
    """把 '8px' / '8pt' 字号字符串解析成像素值 (1pt = 96/72 px)。"""
    s = str(s).strip()
    if s.endswith("px"):
        return float(s[:-2])
    if s.endswith("pt"):
        return float(s[:-2]) * (96.0 / 72.0)
    return float(s)


def estimate_label_size(text: str, font_px: float, canvas_per_px: float) -> tuple[float, float]:
    """估算文字 bbox 在 canvas 单位下的宽高。
    Arial 在 96dpi 下经验值:char_w ≈ font_px * 0.55,height ≈ font_px * 1.2。
    canvas-units = px * canvas_per_px (Bokeh 内部 dpi=96)。
    """
    char_w_px = font_px * 0.55
    char_h_px = font_px * 1.20
    n_chars = max(1, len(text))
    return n_chars * char_w_px * canvas_per_px, char_h_px * canvas_per_px


def deoverlap_labels(
    anchors: list,           # [(x, y), ...] 每个标签的锚点 (节点 / 模块中心)
    widths: list,            # [w, ...]      每个标签 bbox 宽 (canvas 单位)
    heights: list,           # [h, ...]
    extra_pad: float = 0.0,  # 标签间最小空隙
    iterations: int = 80,
    anchor_pull: float = 0.05,  # 每帧把标签拉回锚点的强度
    push_step: float = 0.6,
    is_static: list | None = None,  # bool list:True 的项不能移动(障碍物)
):
    """力导向标签去重叠 + 静态障碍支持:
    - 每对相交 bbox 沿短边方向推开
    - is_static[k]=True 的标签锁定在 anchor,其他标签遇到它时只能"绕开"
    - 同时给可动标签回锚点的弹簧力
    返回与 anchors 同长度的 [(x, y), ...] 新位置列表。"""
    n = len(anchors)
    if n == 0:
        return []
    pos = np.array(anchors, dtype=float)
    anc = pos.copy()
    w = np.array(widths, dtype=float) + extra_pad
    h = np.array(heights, dtype=float) + extra_pad
    static = np.array(is_static if is_static is not None else [False] * n, dtype=bool)

    for _ in range(iterations):
        forces = np.zeros_like(pos)
        any_overlap = False
        for i in range(n):
            for j in range(i + 1, n):
                dx = pos[j, 0] - pos[i, 0]
                dy = pos[j, 1] - pos[i, 1]
                gap_x = (w[i] + w[j]) / 2.0
                gap_y = (h[i] + h[j]) / 2.0
                ox = gap_x - abs(dx)
                oy = gap_y - abs(dy)
                if ox > 0 and oy > 0:
                    any_overlap = True
                    # 双方都 static → 跳过(都不能动)
                    if static[i] and static[j]:
                        continue
                    if ox < oy:
                        sign = 1.0 if dx >= 0 else -1.0
                        if dx == 0:
                            sign = 1.0 if (i % 2 == 0) else -1.0
                        push = (ox + 0.5) * sign
                        # 一方 static → 全部由另一方承担;否则平分
                        if static[i]:
                            forces[j, 0] += push
                        elif static[j]:
                            forces[i, 0] -= push
                        else:
                            forces[i, 0] -= push * 0.5
                            forces[j, 0] += push * 0.5
                    else:
                        sign = 1.0 if dy >= 0 else -1.0
                        if dy == 0:
                            sign = 1.0 if (i % 2 == 0) else -1.0
                        push = (oy + 0.5) * sign
                        if static[i]:
                            forces[j, 1] += push
                        elif static[j]:
                            forces[i, 1] -= push
                        else:
                            forces[i, 1] -= push * 0.5
                            forces[j, 1] += push * 0.5

        # 只有非 static 的标签会移动
        for k in range(n):
            if static[k]:
                continue
            pos[k] = pos[k] + forces[k] * push_step
            # 弹回锚点
            pos[k] = pos[k] + (anc[k] - pos[k]) * anchor_pull
        if not any_overlap:
            break
    return [(float(pos[k, 0]), float(pos[k, 1])) for k in range(n)]


def resolve_node_overlaps(
    pos: dict,
    node_radius_canvas: dict,
    iterations: int = 80,
    padding: float = 0.5,
    shift_frac: float = 0.5,
    min_gap_factor: float = 0.0,
):
    """用 cKDTree 邻近查询迭代推开重叠节点。min_gap_factor > 0 时,
    两节点中心间距会被推至 (R_i + R_j) * (1 + min_gap_factor) + padding,
    即"表面间隙 ≥ min_gap_factor × (R_i + R_j)"。"""
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        return pos

    ids = [n for n in pos if n in node_radius_canvas]
    if len(ids) < 2:
        return pos
    P = np.array([pos[n] for n in ids], dtype=np.float64)
    R = np.array([node_radius_canvas[n] for n in ids], dtype=np.float64)
    r_max = float(R.max())
    mult = 1.0 + max(min_gap_factor, 0.0)
    query_r = 2.0 * r_max * mult + padding
    rng = np.random.default_rng(0)

    for it in range(iterations):
        tree = cKDTree(P)
        pairs = tree.query_pairs(r=query_r)
        if not pairs:
            break
        moved = False
        for i, j in pairs:
            dx = P[j, 0] - P[i, 0]
            dy = P[j, 1] - P[i, 1]
            dist = math.hypot(dx, dy)
            min_dist = (R[i] + R[j]) * mult + padding
            if dist >= min_dist:
                continue
            if dist < 1e-9:
                ang = rng.uniform(0, 2 * math.pi)
                dx, dy = math.cos(ang), math.sin(ang)
                dist = 1.0
            overlap = min_dist - dist
            ux, uy = dx / dist, dy / dist
            shift = overlap * shift_frac
            P[i, 0] -= ux * shift
            P[i, 1] -= uy * shift
            P[j, 0] += ux * shift
            P[j, 1] += uy * shift
            moved = True
        if not moved:
            break

    return {ids[i]: (float(P[i, 0]), float(P[i, 1])) for i in range(len(ids))}
