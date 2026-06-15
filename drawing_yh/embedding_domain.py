# -*- coding: utf-8 -*-
"""Gap-filled domain outlines for categorical embedding atlases.

This module is for the "one readable domain per cell type" UMAP/t-SNE use case:

1. optionally remove isolated low-density cells with a kNN distance filter;
2. keep the main DBSCAN components for each label, while allowing selected
   labels to keep all components;
3. add temporary bridge points between separated main components;
4. draw a concave hull on the real+bridge support;
5. smooth the hull with a closed Catmull-Rom spline and lightly expand it only
   when support points would otherwise fall outside.

The bridge points are only for outline construction. They should not be shown
in the final scatter plot or written back to the source data.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import ConvexHull, cKDTree


@dataclass
class GapfilledOutlineResult:
    """Result returned by :func:`gapfilled_embedding_outline`."""

    curve: np.ndarray
    polygon: object
    label_point: tuple[float, float]
    outside: int
    raw_outside: int
    main_mask: np.ndarray
    component_labels: np.ndarray
    main_component_ids: list[int]
    bridge_edges: list[tuple[int, int, float, int]]
    bridge_points: np.ndarray
    stats: dict


def _polygon_cls():
    from shapely.geometry import Polygon

    return Polygon


def _point_cls():
    from shapely.geometry import Point

    return Point


def _dbscan_cls():
    from sklearn.cluster import DBSCAN

    return DBSCAN


def _polygon_parts(geom):
    if geom.is_empty:
        return []
    if geom.geom_type == "Polygon":
        return [geom]
    return [g for g in getattr(geom, "geoms", []) if g.geom_type == "Polygon"]


def _convex_hull_vertices(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, float)
    if len(points) < 4:
        return points
    hull = ConvexHull(points)
    return points[hull.vertices]


def _dbscan_params(
    n_cells: int,
    *,
    dbscan_eps: float,
    main_cluster_frac: float,
    main_cluster_min: int,
    main_cluster_max: int,
):
    min_samples = max(10, min(80, int(np.sqrt(max(n_cells, 1)) / 4)))
    main_min = max(main_cluster_min, min(main_cluster_max, int(n_cells * main_cluster_frac)))
    return dbscan_eps, min_samples, main_min


def _all_main_components(
    points: np.ndarray,
    *,
    keep_all_components: bool,
    dbscan_eps: float,
    main_cluster_frac: float,
    main_cluster_min: int,
    main_cluster_max: int,
    min_cluster_rel_to_largest: float,
):
    DBSCAN = _dbscan_cls()
    eps, min_samples, main_min = _dbscan_params(
        len(points),
        dbscan_eps=dbscan_eps,
        main_cluster_frac=main_cluster_frac,
        main_cluster_min=main_cluster_min,
        main_cluster_max=main_cluster_max,
    )
    lab = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(points)
    cluster_ids, counts = np.unique(lab[lab >= 0], return_counts=True)
    if len(cluster_ids) == 0:
        return np.ones(len(points), dtype=bool), lab, [], eps, min_samples, main_min

    order = np.argsort(counts)[::-1]
    cluster_ids = cluster_ids[order]
    counts = counts[order]

    if keep_all_components:
        main_ids = [int(cid) for cid in cluster_ids]
        return np.ones(len(points), dtype=bool), lab, main_ids, eps, min_samples, main_min

    main_ids = [int(cid) for cid, count in zip(cluster_ids, counts) if count >= main_min]
    if not main_ids:
        main_ids = [int(cluster_ids[0])]

    largest_count = int(counts[0])
    size_by_id = {int(cid): int(count) for cid, count in zip(cluster_ids, counts)}
    main_ids = [
        cid for cid in main_ids
        if size_by_id[cid] >= largest_count * min_cluster_rel_to_largest
    ]
    if not main_ids:
        main_ids = [int(cluster_ids[0])]
    return np.isin(lab, main_ids), lab, main_ids, eps, min_samples, main_min


def _outline_input(
    points: np.ndarray,
    component_labels: np.ndarray,
    main_component_ids: Sequence[int],
    *,
    max_points: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    hulls = []
    for cid in main_component_ids:
        pts = points[component_labels == cid]
        if len(pts):
            hulls.append(_convex_hull_vertices(pts))
    hull_pts = np.vstack(hulls) if hulls else np.empty((0, 2), float)
    if len(points) <= max_points:
        return np.vstack([points, hull_pts]) if len(hull_pts) else points
    idx = rng.choice(len(points), max_points, replace=False)
    return np.vstack([points[idx], hull_pts])


def _concave_outline_from_points(points: np.ndarray, *, ratio: float):
    from shapely import concave_hull
    from shapely.geometry import MultiPoint, Polygon

    if len(points) < 4:
        support = _convex_hull_vertices(points)
        return Polygon(support), support
    geom = concave_hull(MultiPoint(points), ratio=ratio, allow_holes=False)
    parts = _polygon_parts(geom)
    if not parts:
        support = _convex_hull_vertices(points)
        return Polygon(support), support
    poly = max(parts, key=lambda g: g.area)
    coords = np.asarray(poly.exterior.coords, float)
    return poly, coords[:-1]


def _select_vertices_by_arclength(points: np.ndarray, n: int) -> np.ndarray:
    pts = np.asarray(points, float)
    if len(pts) > 1 and np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]
    if len(pts) <= n:
        return pts

    seg = np.linalg.norm(np.roll(pts, -1, axis=0) - pts, axis=1)
    vertex_s = np.r_[0.0, np.cumsum(seg[:-1])]
    total = float(seg.sum())
    if total <= 0:
        return pts

    targets = np.linspace(0, total, n, endpoint=False)
    chosen = []
    used = set()
    for target in targets:
        delta = np.abs(vertex_s - target)
        delta = np.minimum(delta, total - delta)
        for idx in np.argsort(delta):
            idx = int(idx)
            if idx not in used:
                used.add(idx)
                chosen.append(idx)
                break
    chosen = sorted(chosen, key=lambda i: vertex_s[i])
    return pts[chosen]


def _catmull_rom_closed(
    points: np.ndarray,
    *,
    samples_per_edge: int,
    alpha: float,
) -> np.ndarray:
    pts = np.asarray(points, float)
    if len(pts) < 4:
        return np.vstack([pts, pts[0]])

    curve = []
    n = len(pts)
    for i in range(n):
        p0 = pts[(i - 1) % n]
        p1 = pts[i]
        p2 = pts[(i + 1) % n]
        p3 = pts[(i + 2) % n]

        def tj(ti, pa, pb):
            return ti + np.linalg.norm(pb - pa) ** alpha

        t0 = 0.0
        t1 = tj(t0, p0, p1)
        t2 = tj(t1, p1, p2)
        t3 = tj(t2, p2, p3)
        if min(t1 - t0, t2 - t1, t3 - t2) < 1e-9:
            curve.append(p1)
            continue
        t = np.linspace(t1, t2, samples_per_edge, endpoint=False)[:, None]
        a1 = (t1 - t) / (t1 - t0) * p0 + (t - t0) / (t1 - t0) * p1
        a2 = (t2 - t) / (t2 - t1) * p1 + (t - t1) / (t2 - t1) * p2
        a3 = (t3 - t) / (t3 - t2) * p2 + (t - t2) / (t3 - t2) * p3
        b1 = (t2 - t) / (t2 - t0) * a1 + (t - t0) / (t2 - t0) * a2
        b2 = (t3 - t) / (t3 - t1) * a2 + (t - t1) / (t3 - t1) * a3
        c = (t2 - t) / (t2 - t1) * b1 + (t - t1) / (t2 - t1) * b2
        curve.append(c)
    out = np.vstack(curve)
    return np.vstack([out, out[0]])


def _nearest_hull_pair(h1: np.ndarray, h2: np.ndarray):
    tree = cKDTree(h2)
    dist, idx = tree.query(h1, k=1)
    i = int(np.argmin(dist))
    j = int(idx[i])
    return float(dist[i]), h1[i], h2[j]


def _bridge_tube(
    p0,
    p1,
    *,
    gap_step: float,
    gap_width: float,
    gap_lanes: int,
) -> np.ndarray:
    p0 = np.asarray(p0, float)
    p1 = np.asarray(p1, float)
    vec = p1 - p0
    length = float(np.linalg.norm(vec))
    if length <= gap_step:
        return np.empty((0, 2), float)
    direction = vec / length
    normal = np.array([-direction[1], direction[0]])
    n = max(3, int(np.ceil(length / gap_step)))
    offsets = np.linspace(-gap_width, gap_width, gap_lanes)
    out = []
    for t in np.linspace(0, 1, n + 2)[1:-1]:
        taper = np.sin(np.pi * t) ** 0.55
        center = (1 - t) * p0 + t * p1
        for off in offsets:
            out.append(center + normal * off * taper)
    return np.vstack(out)


def _bridge_points_for_components(
    points: np.ndarray,
    component_labels: np.ndarray,
    main_component_ids: Sequence[int],
    *,
    gap_step: float,
    gap_width: float,
    gap_lanes: int,
):
    if len(main_component_ids) <= 1:
        return np.empty((0, 2), float), []

    hulls = [_convex_hull_vertices(points[component_labels == cid]) for cid in main_component_ids]
    n = len(main_component_ids)
    dmat = np.zeros((n, n), float)
    pairs = {}
    for i in range(n):
        for j in range(i + 1, n):
            dist, p0, p1 = _nearest_hull_pair(hulls[i], hulls[j])
            dmat[i, j] = dmat[j, i] = dist
            pairs[(i, j)] = (p0, p1, dist)

    mst = minimum_spanning_tree(dmat).tocoo()
    bridges = []
    edges = []
    for i, j in zip(mst.row, mst.col):
        a, b = sorted((int(i), int(j)))
        p0, p1, dist = pairs[(a, b)]
        bp = _bridge_tube(p0, p1, gap_step=gap_step, gap_width=gap_width, gap_lanes=gap_lanes)
        if len(bp):
            bridges.append(bp)
        edges.append((int(main_component_ids[a]), int(main_component_ids[b]), float(dist), int(len(bp))))
    if not bridges:
        return np.empty((0, 2), float), edges
    return np.vstack(bridges), edges


def _contains_count(poly, points: np.ndarray):
    exterior = np.asarray(poly.exterior.coords, float)
    inside = MplPath(exterior).contains_points(points, radius=1e-7)
    return int((~inside).sum()), inside


def _outside_count_with_boundary(poly, points: np.ndarray, *, tol: float = 1e-6):
    raw_outside, inside = _contains_count(poly, points)
    if raw_outside == 0:
        return 0, raw_outside
    Point = _point_cls()
    outside_pts = points[~inside]
    true_outside = sum(
        poly.exterior.distance(Point(float(x), float(y))) > tol for x, y in outside_pts
    )
    return int(true_outside), int(raw_outside)


def _outward_buffer_curve(poly, points_for_scale: np.ndarray, *, outward_pad_frac: float):
    scale = max(float(np.ptp(points_for_scale[:, 0])), float(np.ptp(points_for_scale[:, 1])), 1.0)
    pad = max(scale * outward_pad_frac, 0.03)
    expanded = poly.buffer(pad, quad_segs=18)
    if hasattr(expanded, "geoms"):
        expanded = max(expanded.geoms, key=lambda geom: geom.area)
    return np.asarray(expanded.exterior.coords, float), expanded


def label_point_for_outline(curve: np.ndarray) -> tuple[float, float]:
    """Return a label position inside a closed outline curve."""
    Polygon = _polygon_cls()
    Point = _point_cls()
    coords = np.asarray(curve, float)
    poly = Polygon(coords)
    if poly.is_valid and poly.area > 0:
        p = poly.centroid
        if not poly.contains(Point(float(p.x), float(p.y))):
            p = poly.representative_point()
        return float(p.x), float(p.y)
    center = coords.mean(axis=0)
    return float(center[0]), float(center[1])


def gapfilled_embedding_outline(
    points,
    *,
    seed: int = 0,
    keep_all_components: bool = False,
    dbscan_eps: float = 0.55,
    main_cluster_frac: float = 0.010,
    main_cluster_min: int = 120,
    main_cluster_max: int = 800,
    min_cluster_rel_to_largest: float = 0.05,
    concave_ratio: float = 0.24,
    outline_input_max: int = 8000,
    anchor_max: int = 140,
    smooth_anchor_max: int = 70,
    smooth_points_per_edge: int = 16,
    smooth_alpha: float = 0.65,
    gap_step: float = 0.18,
    gap_width: float = 0.42,
    gap_lanes: int = 7,
    outward_pad_frac: float = 0.006,
    log_label: str | None = None,
) -> GapfilledOutlineResult:
    """Build one gap-filled outline around a labelled embedding domain.

    Parameters mirror the final retina atlas defaults. Set
    ``keep_all_components=True`` for genuinely dispersed populations where
    small separated clusters should remain part of the domain (e.g. Neuron in
    the retina atlas). For ordinary subtypes, leave it ``False`` so distant
    small components are not forced into the outline.
    """
    points = np.asarray(points, float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must be a (N, 2) array")
    if len(points) < 4:
        raise ValueError("at least 4 points are required to build a domain outline")

    main_mask, component_labels, main_ids, eps, min_samples, main_min = _all_main_components(
        points,
        keep_all_components=keep_all_components,
        dbscan_eps=dbscan_eps,
        main_cluster_frac=main_cluster_frac,
        main_cluster_min=main_cluster_min,
        main_cluster_max=main_cluster_max,
        min_cluster_rel_to_largest=min_cluster_rel_to_largest,
    )
    main_pts = points[main_mask]
    main_lab = component_labels[main_mask]
    support_mask = np.isin(main_lab, main_ids)
    support_pts = main_pts[support_mask] if support_mask.any() else main_pts
    support_lab = main_lab[support_mask] if support_mask.any() else main_lab

    bridges, edges = _bridge_points_for_components(
        support_pts,
        support_lab,
        main_ids,
        gap_step=gap_step,
        gap_width=gap_width,
        gap_lanes=gap_lanes,
    )
    base_input = (
        support_pts
        if keep_all_components
        else _outline_input(support_pts, support_lab, main_ids, max_points=outline_input_max, seed=seed)
    )
    input_pts = np.vstack([base_input, bridges]) if len(bridges) else base_input

    outline, raw_support = _concave_outline_from_points(input_pts, ratio=concave_ratio)
    anchors = _select_vertices_by_arclength(raw_support, min(len(raw_support), anchor_max))
    smooth_anchors = _select_vertices_by_arclength(raw_support, min(len(raw_support), smooth_anchor_max))

    Polygon = _polygon_cls()
    curve = _catmull_rom_closed(
        smooth_anchors,
        samples_per_edge=smooth_points_per_edge,
        alpha=smooth_alpha,
    )
    poly = Polygon(curve)
    curve_status = "catmull-sparser"
    if (not poly.is_valid) or poly.area < outline.area * 0.60:
        curve = _catmull_rom_closed(anchors, samples_per_edge=10, alpha=smooth_alpha)
        poly = Polygon(curve)
        curve_status = "catmull-fallback"
    if (not poly.is_valid) or poly.area < outline.area * 0.60:
        curve = np.asarray(outline.exterior.coords, float)
        poly = outline
        curve_status = "outline-fallback"

    outside, raw_outside = _outside_count_with_boundary(poly, support_pts)
    if outside:
        curve, poly = _outward_buffer_curve(poly, input_pts, outward_pad_frac=outward_pad_frac)
        outside, raw_outside = _outside_count_with_boundary(poly, support_pts)
        curve_status += "+outward"

    stats = {
        "cells": int(len(points)),
        "main_cells": int(len(main_pts)),
        "outline_support_cells": int(len(support_pts)),
        "excluded_cells": int(len(points) - len(main_pts)),
        "n_components": int(len(main_ids)),
        "n_bridge_points": int(len(bridges)),
        "n_bridge_edges": int(len(edges)),
        "raw_vertices": int(len(raw_support)),
        "n_anchors": int(len(anchors)),
        "curve_status": curve_status,
        "dbscan_eps": float(eps),
        "dbscan_min_samples": int(min_samples),
        "main_min": int(main_min),
    }
    if log_label:
        print(
            f"{log_label}: cells={stats['cells']:,}; main={stats['main_cells']:,}; "
            f"outline_support={stats['outline_support_cells']:,}; "
            f"excluded={stats['excluded_cells']:,}; clusters={stats['n_components']}; "
            f"bridges={stats['n_bridge_points']:,}; edges={stats['n_bridge_edges']}; "
            f"raw_vertices={stats['raw_vertices']} anchors={stats['n_anchors']}; "
            f"curve={curve_status}; outside={outside} raw={raw_outside}; "
            f"dbscan eps={eps:.2f} min_samples={min_samples} main_min={main_min}"
        )

    return GapfilledOutlineResult(
        curve=curve,
        polygon=poly,
        label_point=label_point_for_outline(curve),
        outside=int(outside),
        raw_outside=int(raw_outside),
        main_mask=main_mask,
        component_labels=component_labels,
        main_component_ids=[int(x) for x in main_ids],
        bridge_edges=edges,
        bridge_points=bridges,
        stats=stats,
    )


def embedding_knn_density_mask(
    coords,
    labels=None,
    *,
    q: float = 0.90,
    k: int = 15,
    min_cells: int = 50,
) -> np.ndarray:
    """Return a boolean mask keeping cells in locally dense embedding regions.

    The filter is evaluated per label when ``labels`` is supplied. It is useful
    before domain drawing to avoid isolated speckles being plotted or outlined.
    """
    coords = np.asarray(coords, float)
    if labels is None:
        groups = [(None, np.arange(len(coords)))]
    else:
        labels = np.asarray(labels).astype(str)
        groups = [(lab, np.where(labels == lab)[0]) for lab in pd.unique(labels)]

    keep = np.zeros(len(coords), dtype=bool)
    for _lab, idx in groups:
        pts = coords[idx]
        if len(pts) < min_cells:
            keep[idx] = True
            continue
        kk = min(k + 1, len(pts))
        dist, _ = cKDTree(pts).query(pts, k=kk)
        kth = dist[:, -1]
        keep[idx[kth <= np.quantile(kth, q)]] = True
    return keep


def filter_main_embedding_components(
    coords,
    labels,
    *,
    order=None,
    keep_all_labels: set[str] | None = None,
    min_cells_to_filter: int = 40,
    **outline_kw,
):
    """Keep main DBSCAN components for each label.

    Returns ``(mask, stats_by_label)``. Labels in ``keep_all_labels`` keep all
    DBSCAN components; all other labels drop distant small components using the
    same defaults as :func:`gapfilled_embedding_outline`.
    """
    coords = np.asarray(coords, float)
    labels = np.asarray(labels).astype(str)
    keep_all_labels = set(keep_all_labels or set())
    labels_to_visit = [str(x) for x in (order if order is not None else pd.unique(labels))]
    mask = np.zeros(len(coords), dtype=bool)
    stats_by_label = {}

    defaults = {
        "dbscan_eps": 0.55,
        "main_cluster_frac": 0.010,
        "main_cluster_min": 120,
        "main_cluster_max": 800,
        "min_cluster_rel_to_largest": 0.05,
    }
    defaults.update(outline_kw)

    for label in labels_to_visit:
        idx = np.where(labels == label)[0]
        if len(idx) == 0:
            continue
        if len(idx) < min_cells_to_filter:
            mask[idx] = True
            stats_by_label[label] = {
                "kept": int(len(idx)),
                "total": int(len(idx)),
                "removed": 0,
                "skipped": True,
            }
            continue
        main_mask, lab, main_ids, eps, min_samples, main_min = _all_main_components(
            coords[idx],
            keep_all_components=label in keep_all_labels,
            **defaults,
        )
        mask[idx[main_mask]] = True
        stats_by_label[label] = {
            "kept": int(main_mask.sum()),
            "total": int(len(idx)),
            "removed": int(len(idx) - main_mask.sum()),
            "components": [int(x) for x in main_ids],
            "noise": int((lab < 0).sum()),
            "dbscan_eps": float(eps),
            "dbscan_min_samples": int(min_samples),
            "main_min": int(main_min),
            "skipped": False,
        }
    return mask, stats_by_label


def trim_points_inside_reference_domain(
    coords,
    labels,
    *,
    target_labels: Sequence[str],
    reference_label: str,
    reference_curve=None,
    seed: int = 0,
    **outline_kw,
):
    """Drop target-label points that fall inside a larger reference domain.

    Returns ``(keep_mask, reference_curve, counts)``. This is for cases where a
    small subtype has spillover points inside a dominant subtype domain and the
    atlas semantics are clearer if those points are not displayed.
    """
    coords = np.asarray(coords, float)
    labels = np.asarray(labels).astype(str)
    if reference_curve is None:
        ref_pts = coords[labels == str(reference_label)]
        if len(ref_pts) < 4:
            return np.ones(len(coords), dtype=bool), None, {}
        reference_curve = gapfilled_embedding_outline(ref_pts, seed=seed, **outline_kw).curve

    ref_path = MplPath(np.asarray(reference_curve, float))
    target_mask = np.isin(labels, [str(x) for x in target_labels])
    inside = ref_path.contains_points(coords, radius=1e-6)
    drop = target_mask & inside
    counts = pd.Series(labels[drop]).value_counts().to_dict() if drop.any() else {}
    return ~drop, reference_curve, counts


def draw_gapfilled_embedding_domains(
    ax,
    coords,
    labels,
    *,
    order=None,
    colors=None,
    keep_all_labels: set[str] | None = None,
    fill_alpha: float = 0.24,
    line_alpha: float = 0.92,
    linewidth: float = 0.65,
    linestyle=(0, (4, 3)),
    zorder: float = 0.2,
    seed: int = 0,
    **outline_kw,
) -> dict[str, GapfilledOutlineResult]:
    """Draw filled, dashed gap-filled domains on an existing embedding axis."""
    coords = np.asarray(coords, float)
    labels = np.asarray(labels).astype(str)
    order = [str(x) for x in (order if order is not None else pd.unique(labels))]
    keep_all_labels = set(keep_all_labels or set())
    out: dict[str, GapfilledOutlineResult] = {}

    def _color(label, i):
        if colors is None:
            return f"C{i % 10}"
        if isinstance(colors, Mapping):
            return colors.get(label, f"C{i % 10}")
        return colors[i]

    for i, label in enumerate(order):
        pts = coords[labels == label]
        if len(pts) < 4:
            continue
        res = gapfilled_embedding_outline(
            pts,
            seed=seed + i,
            keep_all_components=label in keep_all_labels,
            **outline_kw,
        )
        col = _color(label, i)
        ax.fill(
            res.curve[:, 0],
            res.curve[:, 1],
            color=col,
            alpha=fill_alpha,
            linewidth=0,
            zorder=zorder,
        )
        ax.plot(
            res.curve[:, 0],
            res.curve[:, 1],
            color=col,
            alpha=line_alpha,
            linewidth=linewidth,
            linestyle=linestyle,
            zorder=zorder + 0.1,
        )
        out[label] = res
    return out


__all__ = [
    "GapfilledOutlineResult",
    "gapfilled_embedding_outline",
    "draw_gapfilled_embedding_domains",
    "embedding_knn_density_mask",
    "filter_main_embedding_components",
    "label_point_for_outline",
    "trim_points_inside_reference_domain",
]
