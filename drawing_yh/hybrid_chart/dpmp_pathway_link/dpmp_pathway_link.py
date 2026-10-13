#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""sig 代谢物(DPMP) x 通路富集关联图 — betaine Fig E/F 构图复刻(草图)。

结构:
- 中间全量窄色带 = 全部检测代谢物(ρ 排序, 60 分箱) → 图的枢纽;
- 左侧 = top10 up / top10 down 放大块(名称 + 色条), 三段式 zoom 指示线 + 朝左小三角喇叭;
- 右侧 = 通路富集棒棒糖(渐进变长竖线束 + 末端空心圆, 长度 = -log10 P);
- 粉/蓝半透明渐变漏斗: 从色带对应段收聚到通路组实际范围。

DIAG=True 时叠加编号诊断线 L1~L5(英文标签, Arial 无中文字形)。
数据: metabolite_age_correlations_significant_all.csv (Spearman p<0.05)
      metabolite_age_correlations_all.csv (全量)
      metabolome_pathway_enrichment/data/enrich_up.csv / enrich_down.csv
"""
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, r'C:/Users/yh109/Documents/GitHub/python-drawing-yh')
import drawing_yh  # noqa: F401  8pt / Arial 全局标准

"""模板化数据入口: 四个 CSV 路径可经 CLI 覆盖。
用法:
  python dpmp_pathway_link.py [--sig S.csv --all A.csv --up U.csv --down D.csv]
默认 = 骨髓微环境项目 BMIF 代谢组(杨航 2026-10-12 定稿参数内置)。
"""
import argparse

_ap = argparse.ArgumentParser()
_ap.add_argument('--sig', default=r'D:/Projects/Bone_Marrow_Aging/metabonomics/analysis/data/by_sex/metabolite_age_correlations_significant_all.csv')
_ap.add_argument('--all', dest='allcsv', default=r'D:/Projects/Bone_Marrow_Aging/metabonomics/analysis/data/by_sex/metabolite_age_correlations_all.csv')
_ap.add_argument('--up', default=r'C:/Users/yh109/Documents/GitHub/projects-yh/Bone_Marrow_Aging/figure/metabolome_pathway_enrichment/data/enrich_up.csv')
_ap.add_argument('--down', default=r'C:/Users/yh109/Documents/GitHub/projects-yh/Bone_Marrow_Aging/figure/metabolome_pathway_enrichment/data/enrich_down.csv')
_ap.add_argument('--outdir', default=None, help='输出目录(默认脚本旁 figure/)')
ARGS, _ = _ap.parse_known_args()

SIG = Path(ARGS.sig)
ALLF = Path(ARGS.allcsv)
UP = Path(ARGS.up)
DN = Path(ARGS.down)
OUTDIR = Path(ARGS.outdir) if ARGS.outdir else Path(__file__).resolve().parent / 'figure'

DIAG = False  # 编号诊断线(已确认, 关闭)

# ——— 数据 ———
sig = pd.read_csv(SIG)
sig = sig[sig.spearman_pval < 0.05]          # 统一 p<0.05 口径(与富集一致)
sig['name'] = sig.COMPOUND_Name.str.replace(r'\*\*', '', regex=True)
sig = sig.sort_values('spearman_rho', ascending=False).reset_index(drop=True)
top_up = sig.head(10)
top_dn = sig.tail(10)

band_df = pd.read_csv(ALLF).sort_values('spearman_rho', ascending=False).reset_index(drop=True)
n_all = len(band_df)
n_up_all = int((band_df.spearman_rho > 0).sum())
n_dn_all = n_all - n_up_all

up = pd.read_csv(UP)
dn = pd.read_csv(DN)
up = up[up['P-value'] < 0.40]
dn = dn[dn['P-value'] < 0.25]


def _term_rows(df):
    cols = list(df.columns)
    return [(re.sub(r'^\[(\w+)\] ', '', r[cols.index('Term')]), float(r[cols.index('P-value')]))
            for r in df.itertuples(index=False, name=None)]


up_terms = _term_rows(up)
dn_terms = _term_rows(dn)
terms = [(t, p, 'up') for t, p in up_terms] + [(t, p, 'dn') for t, p in dn_terms]

n_up_sig = int((sig.spearman_rho > 0).sum())
n_dn_sig = int((sig.spearman_rho < 0).sum())
print('[data] band %d (up %d/down %d); sig %d (up %d/down %d); enrich up %d / down %d'
      % (n_all, n_up_all, n_dn_all, len(sig), n_up_sig, n_dn_sig, len(up_terms), len(dn_terms)))

# ——— 版面 ———
FS = 8
name_w = max(top_up.name.str.len().max(), top_dn.name.str.len().max()) * 0.060 + 0.28
zoom_w = 0.26
band_w = 0.085
blk_w = band_w
funnel_w = 0.34
plot_w = max(len(t) for t, _, _ in terms) * 0.058 + 0.06 + 0.17
tail_w = 0.35
fig_h = 4.25
fig_w = name_w + blk_w + zoom_w + band_w + funnel_w + plot_w + tail_w

fig = plt.figure(figsize=(fig_w, fig_h))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)
ax.axis('off')

ax.add_patch(plt.Rectangle((0, 0), fig_w, fig_h, facecolor='#F2F2F2', edgecolor='none', zorder=-10))
band_ttl = 0.30   # 顶部留白

x_blk1 = name_w
x_band0 = x_blk1 + blk_w + zoom_w
x_band1 = x_band0 + band_w
x_plot0 = x_band1 + funnel_w
x_plot1 = x_plot0 + plot_w

head_y = fig_h - band_ttl - 0.16
ax.text(x_blk1 - 0.02, head_y, 'Metabolites', ha='right', va='center', fontsize=FS)
ax.text(x_plot0, head_y, 'Enriched metabolic pathways', ha='left', va='center', fontsize=FS)



cmap = plt.colormaps['RdBu_r']
vmax = float(band_df.spearman_rho.abs().max())

# ——— 中间: 全量窄色带 ———
band_top = fig_h - band_ttl - 0.50
band_bot = 0.62
band_h = band_top - band_bot
N_BIN = 60
bin_h = band_h / N_BIN
bins = np.array_split(band_df.spearman_rho.to_numpy(), N_BIN)
mean_r = np.array([b.mean() for b in bins])
for i, m in enumerate(mean_r):
    y1 = band_top - i * bin_h
    ax.add_patch(plt.Rectangle((x_band0, y1 - bin_h), band_w, bin_h - 0.004,
                               facecolor=cmap((m + vmax) / (2 * vmax)), edgecolor='none', zorder=2))

up_bin_end = N_BIN * n_up_all / n_all          # 精确分数(不作整箱取整)
rho0_y = band_top - up_bin_end * bin_h          # ρ=0 分界(色带白点)

# ——— 左侧: top10 放大块 ———
row_h = band_h / 21                             # 20 行 + 1 行中缝 = 带高
gap_blk = row_h
up_blk_top = band_top - row_h * 0.5          # 半行偏移: 首行中心低于带顶半行
dn_blk_top = band_top - 11 * row_h - row_h * 0.5
up_ys = [up_blk_top - i * row_h for i in range(10)]
dn_ys = [dn_blk_top - i * row_h for i in range(10)]
for (nm, rho), y in zip(top_up[['name', 'spearman_rho']].itertuples(index=False), up_ys):
    ax.text(x_blk1 - 0.04, y, nm, ha='right', va='center', fontsize=FS, zorder=4)
    ax.add_patch(plt.Rectangle((x_blk1, y - row_h * 0.42), blk_w, row_h * 0.84,
                               facecolor=cmap((rho + vmax) / (2 * vmax)), edgecolor='none', zorder=3))
for (nm, rho), y in zip(top_dn[['name', 'spearman_rho']].itertuples(index=False), dn_ys):
    ax.text(x_blk1 - 0.04, y, nm, ha='right', va='center', fontsize=FS, zorder=4)
    ax.add_patch(plt.Rectangle((x_blk1, y - row_h * 0.42), blk_w, row_h * 0.84,
                               facecolor=cmap((rho + vmax) / (2 * vmax)), edgecolor='none', zorder=3))

# ρ 图例
lg_x, lg_y0, lg_h = 1.05, fig_h - band_ttl - 0.12, 0.55   # 左列文字左侧一点
for i, v in enumerate(np.linspace(vmax, -vmax, 24)):
    ax.add_patch(plt.Rectangle((lg_x, lg_y0 - (i + 1) * lg_h / 24), 0.045, lg_h / 24,
                               facecolor=cmap((v + vmax) / (2 * vmax)), edgecolor='none'))
for v, yy in [(vmax, lg_y0), (0, lg_y0 - lg_h / 2), (-vmax, lg_y0 - lg_h)]:
    ax.text(lg_x - 0.015, yy, '%+.1f' % v, ha='right', va='center', fontsize=FS - 2)
ax.text(lg_x + 0.022, lg_y0 - lg_h - 0.11, 'ρ', ha='center', fontsize=FS - 1)


# ——— zoom 三段式指示(左块 ↔ 色带对应位置) ———
def zoom_pointer(seg_y_center, blk_y_center):
    x_b = x_blk1 + blk_w + 0.015          # 喇叭底边 = 左块右缘
    x_tip = x_b + 0.10
    h = 0.075
    xl0, xl1 = x_tip, x_tip + 0.07
    xr0, xr1 = x_band0 - 0.07, x_band0
    ax.plot([xl0, xl1], [blk_y_center] * 2, color='k', lw=0.6, zorder=6)
    ax.plot([xr0, xr1], [seg_y_center] * 2, color='k', lw=0.6, zorder=6)
    ax.plot([xl1, xr0], [blk_y_center, seg_y_center], color='k', lw=0.6, zorder=6)
    ax.add_patch(plt.Polygon([(x_tip, blk_y_center), (x_b, blk_y_center + h),
                              (x_b, blk_y_center - h)],
                             closed=True, facecolor='white', edgecolor='k', lw=0.8, zorder=7))


up_sig_center = band_top - (n_up_sig / 2) / n_all * band_h
dn_sig_center = band_bot + (n_dn_sig / 2) / n_all * band_h
zoom_pointer(up_sig_center, float(np.mean(up_ys)))
zoom_pointer(dn_sig_center, float(np.mean(dn_ys)))

# ——— 右侧: 通路棒棒糖(竖线束, 左短右长 + 末端空心圆) ———
C_UP, C_DN = '#C51B7D', '#2C7BB6'
col_top = up_blk_top
col_bot = dn_blk_top - 9 * row_h
pitch_all = (col_top - col_bot) / len(terms)
term_ys = [col_top - (i + 0.5) * pitch_all for i in range(len(terms))]
negl = [-math.log10(p) for _, p, _ in terms]
scale = plot_w * 0.78 / max(negl)
HB = min(pitch_all * 0.56, 0.17)
y_name_off = 0.0                   # 通路名居行中心(白描边压线)
y_bar_off = 0.0                    # 竖线束回行中心
for (nm, _p, dirn), y, ng in zip(terms, term_ys, negl):
    col = C_UP if dirn == 'up' else C_DN
    L = max(ng * scale, 0.05)
    NV = max(5, int(round(L / 0.045)))
    yb = y + y_bar_off
    for k in range(NV):
        xk = x_plot0 + (k + 0.5) / NV * L
        hk = HB * (k + 1) / NV
        ax.plot([xk, xk], [yb - hk / 2, yb + hk / 2], color=col, lw=0.65,
                solid_capstyle='butt', zorder=3)
    ax.plot([x_plot0 + L], [yb], marker='o', ms=HB / 0.013, mfc='white', mec=col,
            mew=0.9, zorder=5)
    ax.text(x_plot0 + 0.04, y + y_name_off, nm, ha='left', va='center', fontsize=FS,
            zorder=6, path_effects=[pe.withStroke(linewidth=2.0, foreground='white')])

# 横轴(x = -log10 P): 刻度朝下、数字在轴下
axis_y = min(term_ys) - pitch_all * 0.85
ax.plot([x_plot0, x_plot1], [axis_y] * 2, color='k', lw=0.7)
ax.plot([x_plot0, x_plot0], [axis_y, col_top + 0.10], color='k', lw=0.6)
v = 0.0
while v <= max(negl) + 1e-9:
    xt = x_plot0 + v * scale
    if xt <= x_plot1 - 0.02:
        ax.plot([xt, xt], [axis_y, axis_y - 0.06], color='k', lw=0.6)
        ax.text(xt, axis_y - 0.10, '%g' % v, ha='center', va='top', fontsize=FS)
    v += 0.5
ax.text(x_plot1, axis_y - 0.10, '-log10 P', ha='right', va='top', fontsize=FS)

# ——— 粉/蓝渐变漏斗: 色带段 → 通路组实际范围 ———


def funnel(seg_y0, seg_y1, term_ys_g, color, side, inset_l=0.08, inset_r=0.10):
    """左端 = 色带段(seg_y0..seg_y1), 右端 = 通路组实际跨度;
    靠内侧(接缝 / 相邻通路组一侧)各向内收 inset, 两个衔接块不重叠。"""
    pad_l, pad_r = 0.06, 0.06
    x0 = x_band1 + pad_l
    x1 = x_plot0 - pad_r
    y_l_top, y_l_bot = max(seg_y0, seg_y1), min(seg_y0, seg_y1)
    y_r_top = min(max(term_ys_g) + HB + 0.02, band_top)
    y_r_bot = max(min(term_ys_g) - HB - 0.02, band_bot)
    if side == 'up':          # 内侧 = 下缘
        y_l_bot += inset_l
        y_r_bot += inset_r
    else:                     # 内侧 = 上缘
        y_l_top -= inset_l
        y_r_top -= inset_r
    y_t = max(y_l_top, y_r_top)     # imshow 必须覆盖左右并集, 否则楔形被裁
    y_b = min(y_l_bot, y_r_bot)
    base = mcolors.to_rgb(color)
    cmap_g = mcolors.LinearSegmentedColormap.from_list(
        'g', [(base[0], base[1], base[2], 0.15), (base[0], base[1], base[2], 0.55)])
    grad = np.linspace(0, 1, 64).reshape(1, -1)
    im = ax.imshow(grad, extent=[x0, x1, y_b, y_t], aspect='auto', origin='lower',
                   cmap=cmap_g, zorder=1)
    from matplotlib.patches import Polygon as MPolygon
    quad = MPolygon([(x0, y_l_top), (x0, y_l_bot), (x1, y_r_bot), (x1, y_r_top)],
                    closed=True, transform=ax.transData, facecolor='none', edgecolor='none')
    ax.add_patch(quad)
    im.set_clip_path(quad)


# 漏斗左端 = sig 代谢物在色带上的实际行段(按行 → 分箱边界对齐)
# 左端 = 整个上调段 / 整个下调段(以 rho=0 分界相接)
up_funnel_top = band_top
up_funnel_bot = rho0_y
dn_funnel_top = rho0_y
dn_funnel_bot = band_bot
funnel(up_funnel_top, up_funnel_bot, [t for (_, _, d), t in zip(terms, term_ys) if d == 'up'], C_UP, 'up')
funnel(dn_funnel_top, dn_funnel_bot, [t for (_, _, d), t in zip(terms, term_ys) if d == 'dn'], C_DN, 'dn')

# ——— 编号诊断线(临时) ———
if DIAG:
    _diag = [
        ('L1', band_top - (math.floor((n_up_sig - 1) / n_all * N_BIN) + 1) * bin_h, 'up-sig-end'),
        ('L2', band_top - round(N_BIN * 0.20) * bin_h, 'red-end-20%'),
        ('L3', rho0_y, 'rho=0'),
        ('L4', band_top - round(N_BIN * 0.80) * bin_h, 'blue-start-80%'),
        ('L5', band_top - math.floor((n_all - n_dn_sig) / n_all * N_BIN) * bin_h, 'down-sig-start'),
    ]
    for _id, _y, _lab in _diag:
        ax.plot([0.05, x_plot0 - 0.10], [_y, _y], color='#666666', lw=0.9,
                ls=(0, (4, 3)), zorder=3.5)
        ax.text(0.06, _y + 0.03, '%s %s' % (_id, _lab), ha='left', va='bottom',
                fontsize=FS, color='#666666', zorder=7,
                path_effects=[pe.withStroke(linewidth=2.2, foreground='#F2F2F2')])
    for _c, _y, _dy, _id in [(C_UP, up_funnel_bot, 0.12, 'F1'), (C_DN, dn_funnel_top, -0.16, 'F2')]:
        ax.plot([x_band1, x_plot0 - 0.06], [_y, _y], color=_c, lw=0.9, zorder=3.6)
        ax.text((x_band1 + x_plot0 - 0.06) / 2, _y + _dy, '%s funnel-left' % _id,
                ha='center', va='center', fontsize=FS, color=_c, zorder=7,
                path_effects=[pe.withStroke(linewidth=2.2, foreground='#F2F2F2')])
    print('[diag] ' + ' | '.join('%s@%.0f%%' % (i, (band_top - y) / band_h * 100)
                                 for i, y, _ in _diag))

# ——— 保存(草图快速迭代) ———
out = OUTDIR
out.mkdir(parents=True, exist_ok=True)
out.mkdir(exist_ok=True)
fig.savefig(out / 'figure.svg', bbox_inches='tight')
fig.savefig(out / 'figure.png', dpi=200, bbox_inches='tight')
print('[size] %.2f x %.2f in = %.0f x %.0f mm'
      % (fig_w, fig_h, fig_w * 25.4, fig_h * 25.4))