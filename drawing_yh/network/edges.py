"""
跨模块边束的 Bezier 采样 + bundled_arc(u → trunk_A → trunk_B → v)。
"""
from __future__ import annotations

import math
import numpy as np


def cubic_bezier_samples(p0, p1, p2, p3, n: int = 28):
    t = np.linspace(0.0, 1.0, n)
    u = 1.0 - t
    xs = (u ** 3) * p0[0] + 3 * (u ** 2) * t * p1[0] + 3 * u * (t ** 2) * p2[0] + (t ** 3) * p3[0]
    ys = (u ** 3) * p0[1] + 3 * (u ** 2) * t * p1[1] + 3 * u * (t ** 2) * p2[1] + (t ** 3) * p3[1]
    return xs.tolist(), ys.tolist()


def quadratic_bezier_samples(p0, p1, p2, n: int = 32):
    t = np.linspace(0.0, 1.0, n)
    u = 1.0 - t
    xs = (u ** 2) * p0[0] + 2 * u * t * p1[0] + (t ** 2) * p2[0]
    ys = (u ** 2) * p0[1] + 2 * u * t * p1[1] + (t ** 2) * p2[1]
    return xs.tolist(), ys.tolist()


def bundled_arc(node_u, anchor_A, anchor_B, node_v,
                perp=(0.0, 0.0), bow_offset: float = 0.0,
                n_side: int = 12, n_mid: int = 18):
    """u → anchor_A (cubic) + anchor_A → anchor_B (cubic 弧, 沿 perp 方向偏 bow_offset)
    + anchor_B → v (cubic)。不同 bow_offset → 同对模块内不同缕形成鼓形束。"""
    ax, ay = anchor_A
    bx, by = anchor_B
    dx_ab, dy_ab = bx - ax, by - ay
    L_ab = math.hypot(dx_ab, dy_ab) or 1.0
    tan_step = min(L_ab * 0.20, 50.0)

    # 中段:cubic 弧 — 两个控制点沿 perp 偏移 → 对称弧
    p1_m = (ax + dx_ab * 0.33 + perp[0] * bow_offset,
            ay + dy_ab * 0.33 + perp[1] * bow_offset)
    p2_m = (ax + dx_ab * 0.67 + perp[0] * bow_offset,
            ay + dy_ab * 0.67 + perp[1] * bow_offset)
    xs_m, ys_m = cubic_bezier_samples((ax, ay), p1_m, p2_m, (bx, by), n=n_mid)

    # 边界处切线(和弧对齐,保证 side↔middle 平滑)
    tA = (p1_m[0] - ax, p1_m[1] - ay)
    tA_L = math.hypot(*tA) or 1.0
    tgA = (tA[0] / tA_L, tA[1] / tA_L)
    tB = (bx - p2_m[0], by - p2_m[1])
    tB_L = math.hypot(*tB) or 1.0
    tgB = (tB[0] / tB_L, tB[1] / tB_L)

    # side u → anchor_A:进入锚时切线沿 tgA
    p1_u = (node_u[0] + (ax - node_u[0]) * 0.5,
            node_u[1] + (ay - node_u[1]) * 0.5)
    p2_u = (ax - tgA[0] * tan_step, ay - tgA[1] * tan_step)
    xs1, ys1 = cubic_bezier_samples(node_u, p1_u, p2_u, (ax, ay), n=n_side)

    # side anchor_B → v:离开锚时切线沿 tgB
    p1_v = (bx + tgB[0] * tan_step, by + tgB[1] * tan_step)
    p2_v = (node_v[0] + (bx - node_v[0]) * 0.5,
            node_v[1] + (by - node_v[1]) * 0.5)
    xs3, ys3 = cubic_bezier_samples((bx, by), p1_v, p2_v, node_v, n=n_side)

    # 拼接,去掉中段重复端点
    return xs1 + xs_m[1:-1] + xs3, ys1 + ys_m[1:-1] + ys3
