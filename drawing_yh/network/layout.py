"""
模块化网络的两阶段布局:
  1. 元图 spring_layout(模块当作超节点)→ 模块中心
  2. 每个模块子图独立 spring_layout → 模块内自然形状,平移到对应模块中心

外加 merge_small_modules:把小于阈值的模块按邻居多数投票并入大模块。
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
import numpy as np
import networkx as nx


def merge_small_modules(
    G: nx.Graph,
    module_id: dict,
    min_size: int = 15,
    max_passes: int = 5,
) -> dict:
    """把小于 min_size 的模块的节点逐一并入它在邻接节点里多数所属的"大模块"。
    迭代直到无变化或达到 max_passes。返回新的 node→module_id 映射(按大小降序重新编号 0..K-1)。"""
    module_id = dict(module_id)
    for _ in range(max_passes):
        sizes = Counter(module_id.values())
        small_mods = {m for m, c in sizes.items() if c < min_size}
        if not small_mods:
            break
        # 在小模块里的节点,看邻居都在哪个模块,挑数量最大的"大模块"作为新归属
        # 若所有邻居都在小模块里,这一轮先不动它,等其他小模块先被合并
        new_assign = {}
        for n, m in module_id.items():
            if m not in small_mods:
                continue
            nb_count = Counter(
                module_id[v] for v in G.neighbors(n)
                if v in module_id and module_id[v] not in small_mods
            )
            if nb_count:
                new_assign[n] = nb_count.most_common(1)[0][0]
        if not new_assign:
            # 仍有小模块但邻居全是小模块;退化处理:并入"最大邻居小模块"
            for n, m in module_id.items():
                if m not in small_mods:
                    continue
                nb_count = Counter(
                    module_id[v] for v in G.neighbors(n) if v in module_id
                )
                nb_count.pop(m, None)
                if nb_count:
                    new_assign[n] = nb_count.most_common(1)[0][0]
            if not new_assign:
                break
        module_id.update(new_assign)
    # 重新编号:0..K-1,按大小降序
    sizes = Counter(module_id.values())
    order = [m for m, _ in sizes.most_common()]
    relabel = {m: i for i, m in enumerate(order)}
    return {n: relabel[m] for n, m in module_id.items()}


def layout_meta_then_intra(
    G: nx.Graph,
    module_id: dict,
    seed: int = 1,
    iterations_meta: int = 200,
    iterations_intra: int = 120,
    canvas_size: float = 1800.0,
    module_spacing: float = 1.0,
    intra_radius_base: float = 40.0,
    intra_radius_scale: float = 18.0,
    radial_power: float = 0.55,
    inter_padding: float = 30.0,
    meta_layout: str = "spring",
    uniform_attract_x: float = 0.40,
    uniform_attract_y: float = 0.95,
) -> tuple[dict, dict]:
    """
    两阶段布局,保持模块间自然分布 + 模块内形状自然:
      1) 元图 spring_layout — 模块当作超节点,跨模块边数 = 超边权重 → 模块中心有机散布
      2) 对每个模块子图独立 spring_layout,直接用其自然形状(不强制圆),
         按节点数缩放,平移到对应模块中心

    meta_layout:
      - "spring"  : 默认,元图弹簧布局 + 碰撞分离
      - "circle"  : 把模块中心摆到圆环上(消除"内部模块挡道")
      - "uniform" : 紧凑打包(碰撞 + 中心吸引,无外框边界,矩形画布友好)

    返回 (pos, module_centers) — pos 是 {node: (x, y)},module_centers 是 {module: (x, y)}。
    """
    sizes = Counter(module_id.values())

    # ── 1. 元图:spring_layout 给拓扑位置,再用碰撞分离消除重叠 ──
    print("  ► 元图布局 (spring + 碰撞分离) …")
    meta = nx.Graph()
    meta.add_nodes_from(sizes.keys())
    cross_w = Counter()
    for u, v in G.edges():
        a, b = module_id[u], module_id[v]
        if a != b:
            cross_w[tuple(sorted((a, b)))] += 1
    # 元边权:log 压缩密度,防止大模块间吸引力过强
    for (a, b), w in cross_w.items():
        density = w / math.sqrt(sizes[a] * sizes[b])
        meta.add_edge(a, b, weight=float(math.log1p(density * 10.0)))

    n_meta = meta.number_of_nodes()
    rng = np.random.default_rng(seed)

    if n_meta >= 2:
        meta_k = module_spacing * (3.0 / math.sqrt(n_meta))
        meta_pos = nx.spring_layout(
            meta, weight="weight", k=meta_k,
            iterations=iterations_meta, seed=seed,
        )
    elif n_meta == 1:
        meta_pos = {next(iter(meta.nodes)): (0.0, 0.0)}
    else:
        meta_pos = {}

    # 元图位置缩放到画布(以 canvas_size 为包围盒)
    if meta_pos:
        mxs = np.array([meta_pos[m][0] for m in meta_pos])
        mys = np.array([meta_pos[m][1] for m in meta_pos])
        mspan = max(mxs.max() - mxs.min(), mys.max() - mys.min(), 1e-6)
        mcx, mcy = (mxs.max() + mxs.min()) / 2, (mys.max() + mys.min()) / 2
        meta_scale = canvas_size / mspan
        module_centers = {
            m: ((meta_pos[m][0] - mcx) * meta_scale,
                (meta_pos[m][1] - mcy) * meta_scale)
            for m in meta_pos
        }
    else:
        module_centers = {}

    # ── 1b. 碰撞分离:按模块半径互相推开,使 halo 不重叠 ──
    def _module_radius(m):
        return intra_radius_base + intra_radius_scale * math.sqrt(max(sizes[m], 1))

    print(f"  ► 碰撞分离 (inter_padding={inter_padding:.0f}) …")
    padding = float(inter_padding)
    mod_keys = list(module_centers.keys())
    for it in range(200):
        moved_any = False
        for i in range(len(mod_keys)):
            m1 = mod_keys[i]
            x1, y1 = module_centers[m1]
            r1 = _module_radius(m1)
            for j in range(i + 1, len(mod_keys)):
                m2 = mod_keys[j]
                x2, y2 = module_centers[m2]
                r2 = _module_radius(m2)
                dx, dy = x2 - x1, y2 - y1
                dist = math.hypot(dx, dy)
                min_dist = r1 + r2 + padding
                if dist < min_dist:
                    if dist < 1e-6:
                        # 完全重合:随机方向推开
                        ang = rng.uniform(0, 2 * math.pi)
                        dx, dy = math.cos(ang), math.sin(ang)
                        dist = 1.0
                    overlap = min_dist - dist
                    ux, uy = dx / dist, dy / dist
                    shift = overlap * 0.55
                    module_centers[m1] = (x1 - ux * shift, y1 - uy * shift)
                    module_centers[m2] = (x2 + ux * shift, y2 + uy * shift)
                    moved_any = True
        if not moved_any:
            print(f"    碰撞分离收敛 (iter={it})")
            break

    # 重新居中
    xs = np.array([module_centers[m][0] for m in module_centers])
    ys = np.array([module_centers[m][1] for m in module_centers])
    cx, cy = (xs.max() + xs.min()) / 2, (ys.max() + ys.min()) / 2
    module_centers = {m: (p[0] - cx, p[1] - cy) for m, p in module_centers.items()}

    # ── 1c. 打破外环对称性:小模块径向随机扰动,再次碰撞 ──
    max_r = max(math.hypot(*module_centers[m]) for m in module_centers)
    for m in list(module_centers):
        x, y = module_centers[m]
        r = math.hypot(x, y)
        # 只对位于外侧一半的小模块扰动,保留大模块位置
        if r > max_r * 0.35 and sizes[m] < max(sizes.values()) * 0.15:
            dr = rng.uniform(-0.25, 0.25) * max_r
            da = rng.uniform(-0.4, 0.4)
            ang = math.atan2(y, x) + da
            new_r = max(r + dr, max_r * 0.25)
            module_centers[m] = (new_r * math.cos(ang), new_r * math.sin(ang))

    # 再一次碰撞分离清理扰动后的重叠
    for it in range(100):
        moved_any = False
        for i in range(len(mod_keys)):
            m1 = mod_keys[i]
            x1, y1 = module_centers[m1]
            r1 = _module_radius(m1)
            for j in range(i + 1, len(mod_keys)):
                m2 = mod_keys[j]
                x2, y2 = module_centers[m2]
                r2 = _module_radius(m2)
                dx, dy = x2 - x1, y2 - y1
                dist = math.hypot(dx, dy)
                min_dist = r1 + r2 + padding
                if dist < min_dist:
                    if dist < 1e-6:
                        ang = rng.uniform(0, 2 * math.pi)
                        dx, dy = math.cos(ang), math.sin(ang)
                        dist = 1.0
                    overlap = min_dist - dist
                    ux, uy = dx / dist, dy / dist
                    shift = overlap * 0.55
                    module_centers[m1] = (x1 - ux * shift, y1 - uy * shift)
                    module_centers[m2] = (x2 + ux * shift, y2 + uy * shift)
                    moved_any = True
        if not moved_any:
            break

    xs = np.array([module_centers[m][0] for m in module_centers])
    ys = np.array([module_centers[m][1] for m in module_centers])
    cx, cy = (xs.max() + xs.min()) / 2, (ys.max() + ys.min()) / 2
    module_centers = {m: (p[0] - cx, p[1] - cy) for m, p in module_centers.items()}

    # ── 1d. 可选:把所有模块中心摆到一个圆环上(消除"内部模块挡道"问题) ──
    if meta_layout == "circle" and len(module_centers) >= 3:
        print("  ► 圆环布局:按 spring 角度顺序均匀分布在圆周上")
        # 按 spring 出来的角度排序,均匀分布到圆上;半径=每模块所需半径之和的合理估计
        entries = [
            (m, math.atan2(c[1], c[0])) for m, c in module_centers.items()
        ]
        entries.sort(key=lambda e: e[1])  # 按角度排序
        n_mod = len(entries)
        # 半径:让最近相邻的两模块刚好不相撞
        radii_list = [_module_radius(m) for m in module_centers]
        max_r_mod = max(radii_list) if radii_list else 1.0
        # 弧长 = 2πR / n_mod ≥ 2*max_r_mod + inter_padding
        R = max(
            (2 * max_r_mod + inter_padding) * n_mod / (2 * math.pi),
            canvas_size * 0.45,
        )
        new_centers = {}
        for i, (m, _) in enumerate(entries):
            ang = 2 * math.pi * i / n_mod
            new_centers[m] = (R * math.cos(ang), R * math.sin(ang))
        module_centers = new_centers

    # ── 1e. 可选:紧凑打包(碰撞 + 中心吸引,无外框边界) ──
    #     做法:模块间只有"碰撞"(避免重叠) + 全局"向心吸引"(全员往原点收)
    #     这样模块只会以"互相贴住"的方式紧密聚拢,而不会被边界推到角和边
    elif meta_layout == "uniform" and len(module_centers) >= 3:
        print("  ► 紧凑打包:碰撞 + 中心吸引")
        mod_keys_list = list(module_centers.keys())
        n_mod = len(mod_keys_list)
        # 初始位置:小范围环形 (随机角度) → 给碰撞一个均匀起点
        rng = np.random.default_rng(seed)
        radii = np.array([_module_radius(m) for m in mod_keys_list], dtype=float)
        max_r = float(radii.max())
        init_R = max_r * 1.5
        arr = np.zeros((n_mod, 2), dtype=float)
        for i in range(n_mod):
            ang = 2 * math.pi * i / n_mod + rng.uniform(-0.2, 0.2)
            arr[i, 0] = init_R * math.cos(ang)
            arr[i, 1] = init_R * math.sin(ang)

        # 跨模块边权 (用于反穿越力)
        mod_idx = {m: i for i, m in enumerate(mod_keys_list)}
        pair_w = np.zeros((n_mod, n_mod), dtype=float)
        for u, v in G.edges():
            a, b = module_id[u], module_id[v]
            if a != b and a in mod_idx and b in mod_idx:
                ia, ib = mod_idx[a], mod_idx[b]
                pair_w[ia, ib] += 1
                pair_w[ib, ia] += 1
        log_w = np.log1p(pair_w)
        heavy_thr = max(np.quantile(pair_w[pair_w > 0], 0.5), 30.0) if pair_w.max() > 0 else 0.0

        n_iter = 600
        for it in range(n_iter):
            diff = arr[None, :, :] - arr[:, None, :]
            dist = np.sqrt((diff ** 2).sum(axis=-1)) + 1e-9
            np.fill_diagonal(dist, np.inf)
            unit = diff / dist[..., None]

            # 1) 尺寸感知碰撞:仅在 r_i + r_j + padding 内推开
            min_d = radii[:, None] + radii[None, :] + inter_padding
            np.fill_diagonal(min_d, 0)
            overlap = np.maximum(min_d - dist, 0.0)
            f_total = (-unit * overlap[..., None]).sum(axis=1) * 0.5

            # 2) 矩形偏好中心吸引:默认 x 弱 / y 强 → 团体水平扁平化,填满矩形画布
            attract = arr.copy()
            attract[:, 0] *= uniform_attract_x
            attract[:, 1] *= uniform_attract_y
            f_total = f_total - attract

            # 3) 反穿越力:对重连模块 (a,b),把走廊里的第三方模块 c 推到一侧
            if it > 80 and heavy_thr > 0:
                for ia in range(n_mod):
                    for ib in range(ia + 1, n_mod):
                        w = pair_w[ia, ib]
                        if w < heavy_thr:
                            continue
                        cA = arr[ia]; cB = arr[ib]
                        ab = cB - cA
                        ab_len = float(np.linalg.norm(ab)) + 1e-9
                        if ab_len < radii[ia] + radii[ib]:
                            continue
                        ab_u = ab / ab_len
                        for ic in range(n_mod):
                            if ic == ia or ic == ib:
                                continue
                            ac = arr[ic] - cA
                            t = float(np.dot(ac, ab_u))
                            if not (0.15 * ab_len < t < 0.85 * ab_len):
                                continue
                            perp_vec = ac - t * ab_u
                            perp_d = float(np.linalg.norm(perp_vec)) + 1e-9
                            danger = radii[ic] + 80.0
                            if perp_d < danger:
                                push_mag = log_w[ia, ib] * 4.0 * (1.0 - perp_d / danger)
                                push_mag = min(push_mag, max_r * 0.05)
                                f_total[ic] += (perp_vec / perp_d) * push_mag

            # 4) 步长冷却
            cooling = max(0.15, 1.0 - it / n_iter)
            arr = arr + f_total * cooling

        # 重新居中
        arr = arr - arr.mean(axis=0)
        module_centers = {mod_keys_list[i]: (float(arr[i, 0]), float(arr[i, 1])) for i in range(n_mod)}

    # ── 2. 每模块独立 spring_layout,保留自然形状 ──
    print("  ► 每个模块独立 spring_layout (自然形状) …")
    pos = {}
    grouped = defaultdict(list)
    for n, m in module_id.items():
        if n in G:
            grouped[m].append(n)

    for m, nodes in grouped.items():
        if m not in module_centers:
            continue
        cx, cy = module_centers[m]
        n_count = len(nodes)
        subG = G.subgraph(nodes)

        if subG.number_of_edges() > 0 and n_count >= 3:
            # 较大 k 让 spring 更分散;iterations_intra 保证收敛
            k = 2.5 / math.sqrt(n_count)
            sub_pos = nx.spring_layout(
                subG, seed=seed + m, k=k, iterations=iterations_intra,
            )
        elif n_count >= 2:
            sub_pos = {n: (rng.normal(0, 0.3), rng.normal(0, 0.3)) for n in nodes}
        else:
            sub_pos = {nodes[0]: (0.0, 0.0)}

        # 中心化
        sx = np.array([sub_pos[n][0] for n in nodes])
        sy = np.array([sub_pos[n][1] for n in nodes])
        sx -= sx.mean()
        sy -= sy.mean()

        target_radius = intra_radius_base + intra_radius_scale * math.sqrt(n_count)
        r_n = np.sqrt(sx ** 2 + sy ** 2)
        max_r = max(r_n.max(), 1e-6)

        # 径向 power 重映射:把堆叠在中心的节点向外推(power < 1 扩散内部)
        # 保留角度,仅拉伸半径;避免"中心黑团 + 空 halo"
        r_norm = r_n / max_r
        r_new = np.power(np.maximum(r_norm, 1e-4), radial_power) * target_radius
        scale_per_node = r_new / np.maximum(r_n, 1e-6)
        for i, n in enumerate(nodes):
            pos[n] = (cx + sx[i] * scale_per_node[i], cy + sy[i] * scale_per_node[i])

    return pos, module_centers
