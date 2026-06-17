# -*- coding: utf-8 -*-
"""PCA helpers for drawing_yh.

Public helpers:
- parallel_centroid_guides: draw N parallel guide lines through group centroids
  on a PCA scatter (or any 2D coord), sharing a single direction that
  maximises the spread of those centroids (= top eigenvector of the centroid
  covariance, i.e. 1D PCA on centroids). Lines are perpendicular to that
  direction, so group separation along the optimal axis is most visible.
"""
from __future__ import annotations

from typing import Mapping, Iterable, Sequence

import numpy as np
import pandas as pd
import matplotlib.axes


def parallel_centroid_guides(
    ax: matplotlib.axes.Axes,
    coords: pd.DataFrame,
    *,
    group_col: str = "Group",
    x_col: str = "PC1",
    y_col: str = "PC2",
    color_map: Mapping[str, str],
    group_order: Sequence[str] | None = None,
    linestyle: str = "--",
    linewidth: float = 0.8,
    alpha: float = 0.6,
    zorder: float = 1,
) -> dict:
    """Draw parallel guide lines through each group's centroid.

    The shared direction maximises spread of group centroids (top eigenvector
    of centroid covariance). Each guide line is perpendicular to that axis
    and passes through one centroid; lines are coloured per group.

    Parameters
    ----------
    ax : matplotlib axes
    coords : DataFrame with at least `x_col`, `y_col`, `group_col`
    color_map : group_name -> color hex/name
    group_order : iterable of group names to include (default = unique in coords)
    linestyle, linewidth, alpha, zorder : passed to ax.axline

    Returns
    -------
    dict with keys:
      'direction' : (dx, dy) unit vector of max-spread axis (None if N<2)
      'centroids' : DataFrame with `group_col`, `x_col`, `y_col`
      'n_lines' : int — number of guide lines actually drawn
    """
    groups = list(group_order) if group_order is not None \
        else list(coords[group_col].dropna().unique())

    cents = []
    for grp in groups:
        sub = coords[coords[group_col] == grp]
        if len(sub) == 0:
            continue
        cents.append({group_col: grp,
                      x_col: float(sub[x_col].mean()),
                      y_col: float(sub[y_col].mean())})
    cents_df = pd.DataFrame(cents)

    if len(cents_df) < 2:
        return {"direction": None, "centroids": cents_df, "n_lines": 0}

    arr = cents_df[[x_col, y_col]].to_numpy()
    centered = arr - arr.mean(axis=0)
    cov = np.cov(centered.T)
    eigvals, eigvecs = np.linalg.eig(cov)
    top = eigvecs[:, np.argmax(eigvals.real)].real
    # Normalise to unit (numerical safety)
    norm = np.linalg.norm(top)
    if norm > 0:
        top = top / norm
    perp = np.array([-top[1], top[0]])  # line direction = perpendicular to top

    n_lines = 0
    for _, row in cents_df.iterrows():
        cx, cy = float(row[x_col]), float(row[y_col])
        col = color_map.get(row[group_col], "#888888")
        ax.axline((cx, cy), (cx + perp[0], cy + perp[1]),
                  color=col, ls=linestyle, lw=linewidth, alpha=alpha,
                  zorder=zorder)
        n_lines += 1

    return {"direction": (float(top[0]), float(top[1])),
            "centroids": cents_df, "n_lines": n_lines}


__all__ = ["parallel_centroid_guides"]
