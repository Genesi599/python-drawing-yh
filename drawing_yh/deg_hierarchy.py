# -*- coding: utf-8 -*-
"""
Hierarchical DEG count bars: two-column, colour-dotted, descending rows.

Intended for "Cell count + DEG hierarchy"-style figures (aging analysis):
bold rows = upstream branch / major type (branch-level test), indented rows =
final subtypes (per-subtype test) -- independent tests, counts NOT additive.
Single-subtype branches collapse to one bold row (same test, no duplicate).

Generic interface: any list of hierarchical rows with (down, up) counts works.

Key function
------------
draw_deg_hierarchy(fig, rows, *, xmax=None, fs=8, dot_s=28, up_c="#B2182B",
                   down_c="#2166AC") -> None
    Draws ONE column (label area + diverging bar area) on the current figure
    via add_axes; caller positions multiple columns and supplies fig size.

Row dict fields
---------------
    label   str            display name
    kind    "major"|"sub"  bold+dark-edge vs indented+lighter
    y       float          row y (assign before calling; see group_rows)
    de      (down, up)     counts; None -> "n.t."
    color   hex            dot colour (UMAP-matched palette recommended)

group_rows(groups, de) -> (rows, n_slots)
    Build row dicts from [(branch, [subtypes])] with descending sort by
    total count (branch-level totals between groups, subtype totals within),
    single-subtype branches collapsed to one row. GAP=0.25 between groups.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.ticker import FuncFormatter, MaxNLocator

from . import DEFAULT_FONT_SIZE


def _tint(c, frac):
    r, g, b = to_rgb(c)
    return (r + (1 - r) * frac, g + (1 - g) * frac, b + (1 - b) * frac)


def group_rows(groups, de, gap=0.25):
    """groups = [(branch, [subtypes])]; de[label]=(down,up).

    Descending by total count (branch totals between groups, subtype totals
    within). Single-subtype branch -> one bold row only. Returns (rows,
    n_slots) with row dicts ready for draw_deg_hierarchy.
    """
    def tot(label):
        dn, up = de.get(label, (0, 0))
        return dn + up

    rows, y = [], 0.0
    for br, subs in sorted(groups, key=lambda g: -tot(g[0])):
        rows.append(dict(label=br, kind="major", y=y, de=de.get(br), tested=br in de,
                         color=None))          # color 由调用方注入(项目色板)
        y += 1
        if len(subs) == 1:
            y += gap
            continue
        for sub in sorted(subs, key=lambda s: -tot(s)):
            rows.append(dict(label=sub, kind="sub", y=y, de=de.get(sub), tested=sub in de,
                             color=None))
            y += 1
        y += gap
    return rows, y


def draw_deg_hierarchy(fig, rows, *, rect_lab, rect_bar, ymax, xmax=None,
                       fs=DEFAULT_FONT_SIZE, dot_s=28,
                       up_c="#B2182B", down_c="#2166AC",
                       xlabel="DEG count  (down | up)"):
    """Draw one column (labels + diverging bars) via fig.add_axes.

    rect_lab/rect_bar: figure-fraction rects. ymax: max row y (top margin).
    xmax shared across columns for cross-column comparability (auto if None).
    """
    ax_lab = fig.add_axes(rect_lab)
    ax_lab.set_xlim(0, 1)
    ax_lab.axis("off")
    ax_lab.set_xticks([])
    ax_lab.set_yticks([])
    ax2 = fig.add_axes(rect_bar)
    for ax in (ax_lab, ax2):
        ax.set_ylim(ymax + 0.6, -0.7)

    tested = [r for r in rows if r["tested"]]
    xmax = xmax or max(max(r["de"]) for r in tested)

    for r in rows:
        col = r.get("color") or ("#333333" if r["kind"] == "sub" else "#222222")
        if r["kind"] == "major":
            ax_lab.scatter([0.030], [r["y"]], s=dot_s, color=col, edgecolors="white",
                           linewidths=0.5, zorder=4, clip_on=False)
            ax_lab.text(0.055, r["y"], r["label"], ha="left", va="center",
                        fontweight="bold", fontsize=fs, color="black")
        else:
            ax_lab.scatter([0.215], [r["y"]], s=dot_s, color=col, edgecolors="white",
                           linewidths=0.5, zorder=4, clip_on=False)
            ax_lab.text(0.240, r["y"], r["label"], ha="left", va="center",
                        fontsize=fs, color="#333333")

    # hierarchy guides: branch row -> each subtype row (same branch label match
    # via caller-supplied 'branch' field when present, else skip)
    for r in rows:
        if r["kind"] != "sub" or "branch" not in r:
            continue
    by_label = {r["label"]: r for r in rows}
    for r in rows:
        if r["kind"] == "sub" and r.get("branch") in by_label:
            my = by_label[r["branch"]]["y"]
            ax_lab.plot([0.145, 0.145], [my + 0.4, r["y"]], color="#bbbbbb", lw=0.8, zorder=1)
            ax_lab.plot([0.145, 0.20], [r["y"], r["y"]], color="#bbbbbb", lw=0.8, zorder=1)
    # vertical spine per branch (connect guide stubs)
    for r in rows:
        if r["kind"] == "major":
            subs = [s for s in rows if s["kind"] == "sub" and s.get("branch") == r["label"]]
            if subs:
                ax_lab.plot([0.145, 0.145], [r["y"] + 0.4, subs[-1]["y"]],
                            color="#bbbbbb", lw=0.8, zorder=1)

    ax2.axvline(0, color="#888888", lw=0.6, zorder=2)
    for r in rows:
        if not r["tested"]:
            ax2.text(xmax * 0.03, r["y"], "n.t.", va="center", ha="left",
                     fontsize=fs, color="#999999", fontstyle="italic")
            continue
        dn, up = r["de"]
        edge, lw, tf = (("#222222", 0.8, 0.0) if r["kind"] == "major"
                        else ("none", 0.0, 0.22))
        fw = "bold" if r["kind"] == "major" else "normal"
        bh = 0.66
        ax2.barh(r["y"], -dn, height=bh, color=_tint(down_c, tf), edgecolor=edge,
                 linewidth=lw, zorder=3)
        ax2.barh(r["y"], up, height=bh, color=_tint(up_c, tf), edgecolor=edge,
                 linewidth=lw, zorder=3)
        if dn:
            ax2.text(-dn - xmax * 0.02, r["y"], f"{dn}", va="center", ha="right",
                     fontsize=fs, color="black", fontweight=fw)
        if up:
            ax2.text(up + xmax * 0.02, r["y"], f"{up}", va="center", ha="left",
                     fontsize=fs, color="black", fontweight=fw)
    ax2.set_xlim(-xmax * 1.30, xmax * 1.30)
    ax2.xaxis.set_major_locator(MaxNLocator(5, integer=True))
    ax2.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(int(v))}"))
    ax2.set_xlabel(xlabel, fontsize=fs)
    ax2.set_yticks([])
    for sp in ("top", "right", "left"):
        ax2.spines[sp].set_visible(False)
    ax2.tick_params(length=0, labelsize=fs)
    ax2.grid(axis="x", color="#eeeeee", lw=0.5, zorder=0)
    return ax_lab, ax2
