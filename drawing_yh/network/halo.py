"""
模块 halo:密度渐变背景(替代 halo 椭圆)+ 协方差椭圆 凸包 halo。
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
import numpy as np

from .palette import lighten_hex


def compute_density_background(
    pos: dict,
    module_id: dict,
    palette: list,
    sizes: Counter,
    x_range: tuple,
    y_range: tuple,
    grid_size: int = 520,
    sigma_pixels: float = 14.0,
    max_alpha: float = 0.45,
    alpha_gamma: float = 1.6,
    density_cutoff: float = 0.08,
    overlap_cap: float = 1.2,
):
    """对每个模块的节点位置做栅格化 + gaussian_filter,按模块颜色合成 RGBA 背景。
    稀疏区域 (per-module normalized density < density_cutoff) 完全透明 → 白色底。
    只在节点密集处显示模块色云团。"""
    try:
        from scipy.ndimage import gaussian_filter
    except ImportError:
        return None

    x_min, x_max = x_range
    y_min, y_max = y_range
    W = H = grid_size

    # 收集每个模块的节点栅格坐标
    mod_nodes = defaultdict(list)
    for n, m in module_id.items():
        if n in pos:
            mod_nodes[m].append(pos[n])

    # 每模块原始密度(不做 per-module 归一化,保留实际节点密集度差异)
    module_density = {}
    global_peak = 0.0
    for m, pts in mod_nodes.items():
        hit = np.zeros((H, W), dtype=np.float64)
        for x, y in pts:
            j = int((x - x_min) / (x_max - x_min) * (W - 1))
            i = int((y - y_min) / (y_max - y_min) * (H - 1))
            if 0 <= i < H and 0 <= j < W:
                hit[i, j] += 1.0
        if hit.sum() == 0:
            continue
        d = gaussian_filter(hit, sigma=sigma_pixels)
        module_density[m] = d
        if d.max() > global_peak:
            global_peak = d.max()

    # 逐模块归一化:每个模块按自己的峰值归一到 1.0
    # → 所有模块的 halo 深度上限一致,大模块不再"一枝独秀"
    for m in list(module_density.keys()):
        d = module_density[m]
        peak_m = d.max()
        if peak_m <= 0:
            continue
        d = d / peak_m
        d = np.where(d < density_cutoff, 0.0,
                     (d - density_cutoff) / (1.0 - density_cutoff))
        module_density[m] = d

    # 合成:加权颜色平均 + 总密度决定 alpha
    r_c = np.zeros((H, W), dtype=np.float64)
    g_c = np.zeros((H, W), dtype=np.float64)
    b_c = np.zeros((H, W), dtype=np.float64)
    total = np.zeros((H, W), dtype=np.float64)

    for m, d in module_density.items():
        col = palette[m] if m < len(palette) else "#cbd5e1"
        col = lighten_hex(col, 0.12)
        r, g, b = int(col[1:3], 16), int(col[3:5], 16), int(col[5:7], 16)
        r_c += d * r
        g_c += d * g
        b_c += d * b
        total += d

    safe = np.maximum(total, 1e-9)
    R = np.clip(r_c / safe, 0, 255).astype(np.uint8)
    G = np.clip(g_c / safe, 0, 255).astype(np.uint8)
    B = np.clip(b_c / safe, 0, 255).astype(np.uint8)

    # Alpha: gamma > 1 压缩低密度 → 边缘快速消散为白色
    alpha_raw = np.clip(total / overlap_cap, 0, 1.0)
    A = (np.power(alpha_raw, alpha_gamma) * max_alpha * 255).astype(np.uint8)

    rgba = np.empty((H, W, 4), dtype=np.uint8)
    rgba[..., 0] = R
    rgba[..., 1] = G
    rgba[..., 2] = B
    rgba[..., 3] = A
    return rgba.view(dtype=np.uint32).reshape(H, W)


def compute_module_halos(
    pos: dict,
    module_id: dict,
    buffer_ratio: float = 0.35,
    min_radius: float = 40.0,
    std_mult: float = 2.2,
    min_nodes_for_halo: int = 3,
) -> dict:
    """每个模块画协方差椭圆 halo (中心 + 2σ 外扩)。节点数 < min_nodes_for_halo 跳过。
    返回 {module_id: (hx_list, hy_list)},画在底层用 patches/polygon。"""
    grouped = defaultdict(list)
    for n, m in module_id.items():
        if n in pos:
            grouped[m].append(pos[n])

    halos = {}
    for m, pts_list in grouped.items():
        if len(pts_list) < min_nodes_for_halo:
            continue
        pts = np.array(pts_list)
        cx, cy = pts.mean(axis=0)
        try:
            cov = np.cov(pts.T) + np.eye(2) * (min_radius ** 2) * 0.05
            vals, vecs = np.linalg.eigh(cov)
            vals = np.maximum(vals, 0)
            axis_lengths = std_mult * np.sqrt(vals) * (1 + buffer_ratio)
            axis_lengths = np.maximum(axis_lengths, min_radius)
            t = np.linspace(0, 2 * math.pi, 60)
            ellipse_local = np.column_stack([
                axis_lengths[0] * np.cos(t),
                axis_lengths[1] * np.sin(t),
            ])
            ellipse_world = ellipse_local @ vecs.T
            hx = cx + ellipse_world[:, 0]
            hy = cy + ellipse_world[:, 1]
            halos[m] = (hx.tolist(), hy.tolist())
        except Exception:
            continue
    return halos
