# -*- coding: utf-8 -*-
"""出图后定量门禁: 机械量 overlap / outside / pitch_ratio / marker_ratio。

规则出处: AI-Agent-Hub/hooks/drawing_standards_check.py 的 LAYOUT_GATE_RULE
          + memory/drawing_text_position_audit.md

用法(作图脚本 save_fig 之前):
    from drawing_yh import audit_layout
    print(audit_layout(fig, label="sample_composition"))   # 不达标直接抛 RuntimeError

要点:
  - 内部先 `fig.canvas.draw()` 再量 bbox —— 未渲染时 Text.get_window_extent() 返回过期坐标,
    会产出"假阳性重叠"(图例移位后仍报同一批冲突)。
  - 只统计落在视口内的刻度文字, 否则把不渲染的刻度误报为碰撞。
  - marker_ratio 需要调用方传 row_pitch_in / marker_diam_in; 不传则如实写 N/A。
"""
from __future__ import annotations

import itertools
from typing import Any, Iterable

import numpy as np

__all__ = ["audit_layout", "LayoutAuditError"]

PITCH_MIN, PITCH_MAX = 1.4, 1.9     # 密度目标区间(规则 §5b), 越界只告警不拦
PITCH_HARD_MIN = 1.30                # 低于此值 = 碰撞风险 -> 硬失败
MARKER_RATIO_MAX = 0.90
MIN_FONTSIZE = 8.0                   # 字号锁死 8 pt(硬规则 §5) -> 低于此值硬失败
AREA_TOL = 1.0          # px^2, 忽略亚像素接触
OUTSIDE_TOL = 1.5       # px


class LayoutAuditError(RuntimeError):
    """版式门禁不通过(overlap/outside/pitch/marker 任一越界)。"""


def _ensure_drawable(fig) -> Any:
    """确保 fig 有可用的 Agg renderer 并完成一次 draw。"""
    try:
        fig.canvas.draw()
        return fig.canvas.get_renderer()
    except AttributeError:
        from matplotlib.backends.backend_agg import FigureCanvasAgg

        FigureCanvasAgg(fig)
        fig.canvas.draw()
        return fig.canvas.get_renderer()


def _text_bbox(art, rend):
    """取"文字框本体"bbox。

    ⚠️ 带 arrowprops 的 Annotation, 其 get_window_extent() 会把**引线**算进 bbox —— 三个用引线
    指向不同团的白框明明分开放, 却因引线交叉被判"重叠"(2026-09-16 实测假阳性)。
    有 bbox patch 时改用 patch 的 extent(即可见的文字框)。
    """
    patch = art.get_bbox_patch() if hasattr(art, "get_bbox_patch") else None
    if patch is not None:
        try:
            return patch.get_window_extent(rend)
        except Exception:
            pass
    return art.get_window_extent(rend)


def _in_viewport(labels: Iterable, lo: float, hi: float, axis: str) -> list:
    out = []
    for t in labels:
        if not t.get_text().strip() or not t.get_visible():
            continue
        pos = t.get_position()[0 if axis == "x" else 1]
        if min(lo, hi) - 1e-9 <= pos <= max(lo, hi) + 1e-9:
            out.append(t)
    return out


def _collect_texts(fig) -> list[tuple[str, Any]]:
    """收集参与审计的文字 artist: 视口内刻度 + 轴标题 + 图例文字 + ax.texts 标注。

    返回 (kind, artist); kind 用于字号门禁区分——坐标轴/标题/图例类必须 >=8 pt,
    而 ax.texts 标注允许小字号(规则明确允许格内小缩写与显著性星号)。
    """
    items: list[tuple[str, Any]] = []
    for i, ax in enumerate(fig.get_axes()):
        xlo, xhi = ax.get_xlim()
        ylo, yhi = ax.get_ylim()
        for t in _in_viewport(ax.get_xticklabels(), xlo, xhi, "x"):
            items.append((f"ax{i}.xtick", t))
        for t in _in_viewport(ax.get_yticklabels(), ylo, yhi, "y"):
            items.append((f"ax{i}.ytick", t))
        if ax.title.get_text().strip():
            items.append((f"ax{i}.title", ax.title))
        if ax.xaxis.label.get_text().strip():
            items.append((f"ax{i}.xlabel", ax.xaxis.label))
        if ax.yaxis.label.get_text().strip():
            items.append((f"ax{i}.ylabel", ax.yaxis.label))
        leg = ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                if t.get_text().strip():
                    items.append((f"ax{i}.legend", t))
        for t in ax.texts:
            if t.get_text().strip() and t.get_visible():
                items.append((f"ax{i}.text", t))
    for t in fig.texts:
        if t.get_text().strip() and t.get_visible():
            items.append(("fig.text", t))
    return items


def _pitch_stats(fig, rend) -> dict[str, dict[str, float]]:
    """各轴 8pt 刻度标签的 pitch ratio = 相邻中心距 / 渲染后沿轴尺寸。"""
    stats: dict[str, dict[str, float]] = {}
    for i, ax in enumerate(fig.get_axes()):
        for axis, lims in (("x", ax.get_xlim()), ("y", ax.get_ylim())):
            labels = _in_viewport(ax.get_xticklabels() if axis == "x" else ax.get_yticklabels(),
                                  lims[0], lims[1], axis)
            if len(labels) < 3:
                continue
            boxes = [_text_bbox(t, rend) for t in labels]
            centers = np.sort([(b.x0 + b.x1) / 2 if axis == "x" else (b.y0 + b.y1) / 2
                               for b in boxes])
            size = float(np.median([b.width if axis == "x" else b.height for b in boxes]))
            if size <= 0:
                continue
            ratios = np.diff(centers) / size
            stats[f"ax{i}.{axis}"] = {
                "pitch_min": float(ratios.min()),
                "pitch_median": float(np.median(ratios)),
                "pitch_max": float(ratios.max()),
                "n_labels": len(labels),
            }
    return stats


def audit_layout(
    fig,
    *,
    label: str = "",
    pitch_target: tuple[float, float] = (PITCH_MIN, PITCH_MAX),
    pitch_hard_min: float = PITCH_HARD_MIN,
    min_fontsize: float | None = MIN_FONTSIZE,
    marker_ratio_max: float = MARKER_RATIO_MAX,
    row_pitch_in: float | None = None,
    marker_diam_in: float | None = None,
    raise_on_fail: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """量化版式门禁。返回指标 dict; 不达标且 raise_on_fail=True 时抛 LayoutAuditError。

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        已构建完成的图(函数内部会 draw 一次以取真实 bbox)。
    label : str
        日志前缀(通常传图 stem)。
    row_pitch_in, marker_diam_in : float | None
        散点/点图行高与最大点径(英寸); 两者都给出才算 marker_ratio, 否则记 N/A。
    """
    rend = _ensure_drawable(fig)
    fig_bb = fig.bbox

    texts = _collect_texts(fig)
    overlaps = []
    for (n1, t1), (n2, t2) in itertools.combinations(texts, 2):
        b1, b2 = _text_bbox(t1, rend), _text_bbox(t2, rend)
        if not b1.overlaps(b2):
            continue
        area = ((min(b1.x1, b2.x1) - max(b1.x0, b2.x0))
                * (min(b1.y1, b2.y1) - max(b1.y0, b2.y0)))
        if area > AREA_TOL:
            overlaps.append((round(float(area), 1), f"{n1}:{t1.get_text()}", f"{n2}:{t2.get_text()}"))
    overlaps.sort(reverse=True)

    outside = []
    for n, t in texts:
        b = _text_bbox(t, rend)
        over = max(fig_bb.x0 - b.x0, b.x1 - fig_bb.x1, fig_bb.y0 - b.y0, b.y1 - fig_bb.y1)
        if over > OUTSIDE_TOL:
            outside.append((round(float(over), 1), f"{n}:{t.get_text()}"))

    # 图内标注越出所属坐标区(2026-09-16 加): 质心编号被偏移推到自家面板外(如落在 colorbar 下方),
    # outside 检查只看整幅画布抓不到, 但读者会看到"孤立数字"。
    stray = []
    for ai, ax in enumerate(fig.get_axes()):
        ab = ax.get_window_extent(rend)
        for t in getattr(ax, "texts", []):
            if not t.get_text().strip():
                continue
            b = _text_bbox(t, rend)
            cx, cy = (b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2
            if not (ab.x0 - 2.0 <= cx <= ab.x1 + 2.0 and ab.y0 - 2.0 <= cy <= ab.y1 + 2.0):
                stray.append((ai, t.get_text()[:24]))

    pitch = _pitch_stats(fig, rend)
    # 硬失败: 密度低到有碰撞风险; 软告警: 偏松(>目标上限)或略紧
    hard_pitch = {k: v for k, v in pitch.items() if v["pitch_median"] < pitch_hard_min}
    warn_pitch = {k: v for k, v in pitch.items()
                  if k not in hard_pitch
                  and (v["pitch_median"] < pitch_target[0] or v["pitch_median"] > pitch_target[1])}

    # 字号门禁: 坐标轴/标题/图例类文字 >= 8 pt(标注除外)
    small = []
    if min_fontsize:
        for kind, t in texts:
            txt = t.get_text().strip()
            # 只豁免**符号类标注**(ax.text 里的星号/†‡§/纯数字角标);
            # 刻度·轴标签·标题·图例一律吃 8 pt 硬规则 —— 曾因"≤2 字符"豁免
            # 让 1 字符刻度(Y/O/P)与 5.8 pt 注释白框都溜过门禁(2026-09-16 两次收紧)
            if kind.endswith(".text") and (len(txt) <= 2
                                           or all(ch in "*†‡§#0123456789" for ch in txt)):
                continue
            fs = t.get_fontsize()
            if isinstance(fs, (int, float)) and fs < min_fontsize - 1e-9:
                small.append((round(float(fs), 1), f"{kind}:{t.get_text()[:24]}"))

    if row_pitch_in and marker_diam_in:
        marker_ratio: float | str = round(float(marker_diam_in) / float(row_pitch_in), 3)
    else:
        marker_ratio = "N/A"

    report = {
        "label": label or fig.get_label() or "figure",
        "overlap": len(overlaps),
        "overlap_detail": overlaps[:12],
        "outside": len(outside),
        "outside_detail": outside[:12],
        "stray_text": len(stray),
        "stray_detail": stray[:12],
        "pitch_ratio": pitch,
        "pitch_fail_axes": hard_pitch,
        "pitch_warn_axes": warn_pitch,
        "small_fonts": small,
        "marker_ratio": marker_ratio,
        "passed": True,
        "failures": [],
    }

    fails: list[str] = []
    if overlaps:
        fails.append(f"overlap={len(overlaps)} (须 0): {overlaps[:4]}")
    if outside:
        fails.append(f"outside={len(outside)} (须 0): {outside[:4]}")
    if stray:
        fails.append(f"标签越出坐标区={len(stray)} (须 0): {stray[:4]}")
    if hard_pitch:
        fails.append(f"pitch_ratio < {pitch_hard_min}(碰撞风险): {hard_pitch}")
    if small:
        fails.append(f"字号 < {min_fontsize} pt(硬规则须锁 8 pt): {small[:6]}")
    warns = []
    if warn_pitch:
        warns.append(f"密度偏离目标 {pitch_target}(非阻塞): "
                     + ", ".join(f"{k}:{v['pitch_median']:.2f}" for k, v in warn_pitch.items()))
    report["warnings"] = warns
    if isinstance(marker_ratio, float) and marker_ratio > marker_ratio_max:
        fails.append(f"marker_ratio={marker_ratio} > {marker_ratio_max}")
    report["failures"] = fails
    report["passed"] = not fails

    if verbose:
        print(f"[audit_layout] {report['label']}: overlap={report['overlap']} "
              f"outside={report['outside']} stray={report['stray_text']} marker_ratio={marker_ratio} "
              f"pitch=" + ", ".join(f"{k}:{v['pitch_median']:.2f}" for k, v in pitch.items())
              + ("  OK" if report["passed"] else "  FAIL"))
        for w in warns:
            print(f"[audit_layout]   warn: {w}")

    if fails and raise_on_fail:
        raise LayoutAuditError(f"{report['label']} 版式门禁不通过 -> " + " | ".join(fails))
    return report
