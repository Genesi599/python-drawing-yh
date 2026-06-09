"""Volcano plot — drawing-yh checklist-compliant fig-in-fig-out template.

- 自动检测 y 是 raw p (0-1) 还是 -log10 p,raw 时 clip lower=1e-300 自动转
- 三层散点(weak grey + sig_up + sig_dn),配色集中,与 OYdeg 系列 WARM/COOL 兼容
- 阈值线 (x=0 + horizontal nominal p)
- top-N label(可选)+ adjustText,按 p 或距原点排序
- 图例外置底部水平 (`frameon=False`,checklist 8.6)
- 全英文 + unicode minus
"""
from __future__ import annotations
from typing import Optional, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def volcano(
    df: pd.DataFrame,
    x: str,
    y: str,
    *,
    p_cutoff: float = 0.05,
    fc_cutoff: Optional[float] = None,
    up_col: str = "#C0392B",
    down_col: str = "#2166AC",
    weak_col: str = "#bbbbbb",
    label_col: Optional[str] = None,
    label_n: int = 10,
    label_by: str = "p",  # "p" or "dist"
    label_fontsize: float = 6.5,
    label_color_by_dir: bool = True,
    point_size: float = 5.5,
    weak_size: float = 2.2,
    alpha: float = 0.85,
    weak_alpha: float = 0.42,
    rasterize_weak: bool = True,
    figsize: tuple = (3.6, 3.0),
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    add_threshold_lines: bool = True,
    p_for_hline: float = 0.05,
    legend: bool = True,
    legend_ncol: int = 3,
    adjust_text_kwargs: Optional[dict] = None,
    ax=None,
):
    """Volcano scatter.

    Parameters
    ----------
    df : DataFrame  必含 ``x`` 和 ``y`` 列;若指定 ``label_col`` 还需该列。
    x : str         effect (log2FC / score diff) 列名。
    y : str         p-value 列名(0-1 范围)或 -log10 p 列名(自动检测)。
    p_cutoff : float
        显著性 p 阈值(显色 sig up/down)。
    fc_cutoff : float | None
        x effect 阈值(取 |x|≥fc_cutoff 才算 sig);None 不卡 x。
    up_col / down_col / weak_col : color
        三层散点配色。默认与 OYdeg 系列 WARM/COOL 一致(#C0392B / #2166AC)。
    label_col : str | None
        ID 列,用于标注 top-N 点;None 关闭标注。
    label_n : int
        标 top-N 点(`label_by` 排)。
    label_by : "p" | "dist"
        排序依据:"p" 按 y 降序 / "dist" 按距原点距离。
    label_color_by_dir : bool
        True → label 颜色按 effect 方向上色(up=红 / down=蓝),与点同色。
    point_size / weak_size : float
        sig 点 / weak 点散点大小。
    alpha / weak_alpha : float
        sig / weak 点透明度。
    rasterize_weak : bool
        weak 点光栅化(减 PDF/SVG 大小,数据大时必开)。
    figsize : (W, H)
        单栏 3.35-3.6 in 推荐(checklist 5)。
    title : str | None
        左对齐 small title;None 不显示。
    xlabel / ylabel : str | None
        None → 默认 ``f"Effect ({x})"`` / ``"−log10(p)"`` (unicode minus)。
    add_threshold_lines : bool
        画 x=0 垂线 + y=-log10(p_for_hline) 水平虚线。
    legend : bool
        是否显示图例(三类:up / down / weak)。
    legend_ncol : int
        图例列数(默认 3 水平排)。
    adjust_text_kwargs : dict | None
        覆盖默认 adjustText 参数(如 ``force_text`` / ``expand``)。
    ax : matplotlib Axes | None
        在已有 ax 上画;None → 新建 fig+ax。

    Returns
    -------
    (fig, ax) : matplotlib Figure, Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    work = df.copy()
    xv = pd.to_numeric(work[x], errors="coerce")
    yraw = pd.to_numeric(work[y], errors="coerce")

    # 自动检测 y 是 raw p (0-1) 还是 -log10 p
    if yraw.dropna().between(0, 1).all():
        yv = -np.log10(yraw.clip(lower=1e-300))
    else:
        yv = yraw

    work["_x"] = xv
    work["_y"] = yv

    sig_y = yv > -np.log10(p_cutoff)
    sig_mask = sig_y
    if fc_cutoff is not None:
        sig_mask = sig_mask & (xv.abs() >= fc_cutoff)
    sig_up = sig_mask & (xv > 0)
    sig_dn = sig_mask & (xv < 0)
    weak = ~sig_mask

    ax.scatter(xv[weak], yv[weak], s=weak_size, color=weak_col,
               alpha=weak_alpha, edgecolor="none", rasterized=rasterize_weak,
               label=f"p≥{p_cutoff}")
    ax.scatter(xv[sig_up], yv[sig_up], s=point_size, color=up_col,
               alpha=alpha, edgecolor="none", label=f"Up (p<{p_cutoff})")
    ax.scatter(xv[sig_dn], yv[sig_dn], s=point_size, color=down_col,
               alpha=alpha, edgecolor="none", label=f"Down (p<{p_cutoff})")

    if add_threshold_lines:
        ax.axvline(0, color="#888888", lw=0.5, zorder=1)
        ax.axhline(-np.log10(p_for_hline), color="#aaaaaa", lw=0.5, ls="--", zorder=1)

    if label_col is not None and label_n > 0:
        if label_by == "p":
            top = work.sort_values("_y", ascending=False).head(label_n)
        else:
            xr = max(float(xv.abs().max()), 1e-6)
            yr = max(float(yv.max()), 1e-6)
            work["_d"] = np.sqrt((xv / xr) ** 2 + (yv / yr) ** 2)
            top = work.nlargest(label_n, "_d")
        texts = []
        for _, row in top.iterrows():
            col = (up_col if row["_x"] > 0 else down_col) if label_color_by_dir else "#222222"
            texts.append(ax.text(row["_x"], row["_y"], str(row[label_col]),
                                 fontsize=label_fontsize, color=col, zorder=5))
        try:
            from adjustText import adjust_text
            kw = dict(arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
                      expand=(1.4, 1.8), force_text=(0.6, 0.8),
                      force_points=(0.3, 0.4))
            if adjust_text_kwargs:
                kw.update(adjust_text_kwargs)
            adjust_text(texts, ax=ax, **kw)
        except Exception as exc:
            print(f"[volcano] adjustText skipped: {exc}")

    if title is not None:
        ax.set_title(title, loc="left", fontsize=8)
    ax.set_xlabel(xlabel if xlabel is not None else f"Effect ({x})")
    ax.set_ylabel(ylabel if ylabel is not None else "−log10(p)")

    if legend:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
                  ncol=legend_ncol, frameon=False, fontsize=6.5,
                  handletextpad=0.3, columnspacing=0.9, markerscale=1.8)

    fig.tight_layout()
    return fig, ax


# ── 旧 file-based 接口,向后兼容 (deprecated,新代码请用 ``volcano()``) ──
def create_volcano_plot(
        input_file, output_file="Volcano_plot.png",
        x_threshold=0.5, y_threshold=None,
        lfc_col="log2FoldChange", p_col="padj",
        id_col="GeneName"):
    """Deprecated file-based wrapper kept for backward compatibility."""
    import warnings
    warnings.warn("create_volcano_plot(file-based) deprecated; use volcano(df, ...)",
                  DeprecationWarning, stacklevel=2)
    df = pd.read_csv(input_file)
    if y_threshold is None:
        y_threshold = -np.log10(0.05)
    p_cutoff = 10 ** (-y_threshold)
    fig, ax = volcano(df, x=lfc_col, y=p_col,
                      p_cutoff=p_cutoff, fc_cutoff=x_threshold,
                      label_col=id_col, label_n=10, label_by="p",
                      figsize=(9, 6))
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    fig.savefig(output_file.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
