#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Cell-style DA(differential-abundance) dot plot —— 按分类聚合的衰老方向概览。

每行 = 一个 sub 分类(sub_col)。dot size = mean|r|(mean_absr 列),
dot color = DA score(DA 列, [-1,1], 蓝=下调主导 / 红=上调主导),
右侧标 nsig/n,左侧 super 色条(super_col)按大类分组。

data 须已按 super 分组连续排序好(同一 super 的行相邻),并含列:
    sub_col, super_col, 'n', 'nsig', 'mean_absr', 'DA'

代谢(metabolome_subpathways: super=super-pathway, sub=sub-pathway)与
蛋白(proteome_go_summary: super=subcellular location, sub=GO term)共用本函数。
"""
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.use14corefonts'] = False
mpl.rcParams['font.family'] = 'Arial'

PT_PER_IN = 72.0
DEFAULT_DA_CMAP = LinearSegmentedColormap.from_list(
    'da_bwr', ['#2166AC', '#F7F7F7', '#B2182B'])


def _measure_text_width_in(strings, fontsize, family='Arial'):
    fig = plt.figure(figsize=(2, 2), dpi=PT_PER_IN)
    r = fig.canvas.get_renderer()
    max_w = 0.0
    for s in strings:
        t = fig.text(0, 0, str(s), fontsize=fontsize, family=family)
        bb = t.get_window_extent(renderer=r)
        max_w = max(max_w, bb.width / PT_PER_IN)
        t.remove()
    plt.close(fig)
    return max_w


def _check_label_overlap(fig, texts):
    r = fig.canvas.get_renderer()
    bbs = [t.get_window_extent(renderer=r) for t in texts]
    for i in range(len(bbs)):
        for j in range(i + 1, len(bbs)):
            a, b = bbs[i], bbs[j]
            if not (a.y1 <= b.y0 or b.y1 <= a.y0):
                if not (a.x1 <= b.x0 or b.x1 <= a.x0):
                    return True
    return False


def da_dot_summary(
        data, *, sub_col, super_col='super',
        super_colors, title='', out_path=None,
        font_size=8, text_color='black',
        r_min=0.1, r_max=0.6, dot_smin=12, dot_smax=110,
        leg_rs=(0.1, 0.3, 0.5), da_cmap=None,
        pad_in=0.03, bar_w_in=0.09, dot_w_in=0.85, leg_w_in=0.90,
        col_gap_in=0.04, top_pad_in=0.22, bot_pad_in=0.38):
    """画 DA dot summary;若给 out_path 则写 .svg/.png/.pdf 三格式。返回 (fig, axes)。"""
    if da_cmap is None:
        da_cmap = DEFAULT_DA_CMAP
    n = len(data)
    row_pitch_in = font_size * 1.4 / PT_PER_IN

    def size_from_absr(absr):
        a = np.clip(absr, r_min, r_max)
        return dot_smin + (dot_smax - dot_smin) * (a - r_min) / (r_max - r_min)

    # 量一次最长 sub 标签,据此定列宽 + figsize
    label_w = _measure_text_width_in(data[sub_col].tolist(), font_size) + 0.03
    axes_h_in = n * row_pitch_in
    fig_h = top_pad_in + axes_h_in + bot_pad_in
    fig_w = (pad_in + bar_w_in + col_gap_in + label_w + col_gap_in
             + dot_w_in + col_gap_in + leg_w_in + pad_in)

    fig = plt.figure(figsize=(fig_w, fig_h))
    x_bar = pad_in / fig_w
    w_bar = bar_w_in / fig_w
    x_label = x_bar + w_bar + col_gap_in / fig_w
    w_label = label_w / fig_w
    x_dot = x_label + w_label + col_gap_in / fig_w
    w_dot = dot_w_in / fig_w
    x_leg = x_dot + w_dot + col_gap_in / fig_w
    w_leg = leg_w_in / fig_w
    y_ax = bot_pad_in / fig_h
    h_ax = axes_h_in / fig_h

    ax_bar = fig.add_axes([x_bar, y_ax, w_bar, h_ax])
    ax_label = fig.add_axes([x_label, y_ax, w_label, h_ax])
    ax_dot = fig.add_axes([x_dot, y_ax, w_dot, h_ax])
    ax_leg = fig.add_axes([x_leg, y_ax, w_leg, h_ax])

    y = np.arange(n)[::-1]
    ylim = (-0.7, n - 0.3)

    # — super 色条(连续 run 合并成一块) —
    supers = data[super_col].tolist()
    present_supers = []
    for s in supers:
        if s not in present_supers:
            present_supers.append(s)
    runs, start = [], 0
    for i in range(1, n + 1):
        if i == n or supers[i] != supers[start]:
            runs.append((supers[start], start, i - 1))
            start = i
    for s, a, b in runs:
        y_top = y[a] + 0.5
        y_bot = y[b] - 0.5
        ax_bar.add_patch(Rectangle(
            (0, y_bot), 1, y_top - y_bot,
            facecolor=super_colors.get(s, '#AAAAAA'),
            edgecolor='white', linewidth=0.5))
    ax_bar.set_xlim(0, 1); ax_bar.set_ylim(ylim)
    ax_bar.set_xticks([]); ax_bar.set_yticks([])
    for sp in ax_bar.spines.values():
        sp.set_visible(False)

    # — sub 标签(右对齐) —
    ax_label.set_xlim(0, 1); ax_label.set_ylim(ylim)
    label_texts = []
    for yi, name in zip(y, data[sub_col]):
        t = ax_label.text(1.0, yi, str(name), ha='right', va='center',
                          fontsize=font_size, color=text_color)
        label_texts.append(t)
    ax_label.set_xticks([]); ax_label.set_yticks([])
    for sp in ax_label.spines.values():
        sp.set_visible(False)

    # — dot plot —
    sizes = size_from_absr(data['mean_absr'].values)
    norm = Normalize(vmin=-1, vmax=1)
    colors = da_cmap(norm(data['DA'].values))
    ax_dot.set_xlim(0, 1); ax_dot.set_ylim(ylim)
    for yi in y:
        ax_dot.axhline(yi, color='#EFEFEF', lw=0.4, zorder=0)
    ax_dot.axvline(0.3, color='#CFCFCF', lw=0.5, linestyle='--', zorder=0)
    ax_dot.scatter(np.full(n, 0.3), y, s=sizes, c=colors,
                   edgecolor='black', linewidth=0.5, zorder=3)
    for yi, n_i, ns_i in zip(y, data['n'], data['nsig']):
        ax_dot.text(0.62, yi, f'{ns_i}/{n_i}', ha='left', va='center',
                    fontsize=font_size, color=text_color)
    ax_dot.set_xticks([]); ax_dot.set_yticks([])
    for sp in ax_dot.spines.values():
        sp.set_visible(False)

    # — 右侧图例:size + DA colorbar —
    ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1); ax_leg.axis('off')
    leg_sizes = size_from_absr(np.array(list(leg_rs)))
    y0 = 0.97
    ax_leg.text(0.5, y0, 'mean |r|', ha='center', va='top',
                fontsize=font_size, color=text_color)
    for i, (r, s) in enumerate(zip(leg_rs, leg_sizes)):
        yy = y0 - 0.07 - i * 0.06
        ax_leg.scatter(0.25, yy, s=s, facecolor='#BBBBBB',
                       edgecolor='black', linewidth=0.5)
        ax_leg.text(0.42, yy, f'{r:.1f}', ha='left', va='center',
                    fontsize=font_size, color=text_color)
    cb_x0, cb_y0, cb_w, cb_h = 0.12, 0.55, 0.60, 0.035
    cax = ax_leg.inset_axes([cb_x0, cb_y0, cb_w, cb_h])
    cb = plt.colorbar(ScalarMappable(norm=norm, cmap=da_cmap),
                      cax=cax, orientation='horizontal')
    cb.set_ticks([-1, 0, 1]); cb.set_ticklabels(['-1', '0', '1'])
    cb.ax.tick_params(labelsize=font_size, length=2, pad=1, colors=text_color)
    cb.outline.set_linewidth(0.4)
    ax_leg.text(cb_x0 + cb_w / 2, cb_y0 + cb_h + 0.02, 'DA score',
                ha='center', va='bottom', fontsize=font_size,
                color=text_color, transform=ax_leg.transAxes)

    # — 底部 super 图例 —
    handles = [Rectangle((0, 0), 1, 1,
                         facecolor=super_colors.get(s, '#AAAAAA'),
                         edgecolor='none') for s in present_supers]
    fig.legend(handles, present_supers, loc='lower center', ncol=len(handles),
               frameon=False, fontsize=font_size, labelcolor=text_color,
               bbox_to_anchor=(0.5, 0.005), handlelength=0.8, handleheight=0.8,
               columnspacing=0.9, handletextpad=0.3)

    fig.suptitle(title, fontsize=font_size, fontweight='bold',
                 color=text_color, y=1 - 0.01)

    # — 标签重叠则增高重排 —
    for _ in range(6):
        fig.canvas.draw()
        if not _check_label_overlap(fig, label_texts):
            break
        fig_h *= 1.1
        fig.set_size_inches(fig_w, fig_h)
        y_ax = bot_pad_in / fig_h
        h_ax = (fig_h - top_pad_in - bot_pad_in) / fig_h
        for ax in (ax_bar, ax_label, ax_dot, ax_leg):
            pos = ax.get_position()
            ax.set_position([pos.x0, y_ax, pos.width, h_ax])

    if out_path:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        for suf, dpi in [('.svg', PT_PER_IN), ('.png', 220), ('.pdf', None)]:
            fig.savefig(out_path.with_suffix(suf), dpi=dpi, bbox_inches='tight',
                        facecolor='white', metadata={'Creator': None, 'Date': None})
    return fig, (ax_bar, ax_label, ax_dot, ax_leg)
