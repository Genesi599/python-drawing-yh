#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build a static gallery of representative drawing-yh chart examples."""
from __future__ import annotations

import html
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import drawing_yh
from drawing_yh import (
    OKABE_ITO,
    marker_dotplot,
    render_mean_dumbbell,
    save_fig,
    venn_diagram,
)
from drawing_yh.chord import chord_diagram
from drawing_yh.network import hub_spoke


ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
GEN = OUT / "generated"


@dataclass(frozen=True)
class Example:
    slug: str
    category: str
    sub_category: str
    title: str
    image: Path
    source: Path
    use_case: str
    note: str
    tags: tuple[str, ...]


def _rel_url(path: Path) -> str:
    rel = Path(os.path.relpath(path.resolve(), OUT.resolve())).as_posix()
    return quote(rel, safe="/:#?&=%")


def _root_rel(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def _anchor(prefix: str, value: str) -> str:
    chars = []
    last_dash = False
    for ch in value.lower():
        if ch.isalnum():
            chars.append(ch)
            last_dash = False
        elif not last_dash:
            chars.append("-")
            last_dash = True
    return f"{prefix}-{''.join(chars).strip('-') or 'section'}"


def _style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e6e6e6", linewidth=0.5)


def generate_box_plot() -> Path:
    rng = np.random.default_rng(3)
    vals = [
        rng.normal(2.2, 0.35, 28),
        rng.normal(2.9, 0.42, 28),
        rng.normal(3.4, 0.38, 28),
    ]
    out = GEN / "box_plot.png"
    fig, ax = plt.subplots(figsize=(3.35, 2.3))
    bp = ax.boxplot(vals, patch_artist=True, widths=0.55, showfliers=False)
    for patch, color in zip(bp["boxes"], [OKABE_ITO[0], OKABE_ITO[2], OKABE_ITO[5]]):
        patch.set(facecolor=color, alpha=0.65, edgecolor="#333333", linewidth=0.8)
    for key in ("whiskers", "caps", "medians"):
        for line in bp[key]:
            line.set(color="#333333", linewidth=0.8)
    jitter_x = np.repeat(np.arange(1, 4), 28) + rng.normal(0, 0.045, 84)
    ax.scatter(jitter_x, np.concatenate(vals), s=8, c="#222222", alpha=0.45, linewidths=0)
    ax.set_xticklabels(["Young", "Middle", "Old"])
    ax.set_ylabel("Score")
    ax.set_title("Group distribution")
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def generate_scatter_linear_fit() -> Path:
    rng = np.random.default_rng(5)
    x = np.linspace(20, 85, 46)
    y = 0.035 * x + rng.normal(0, 0.28, len(x)) + 1.0
    coef = np.polyfit(x, y, 1)
    out = GEN / "scatter_linear_fit.png"
    fig, ax = plt.subplots(figsize=(3.35, 2.25))
    ax.scatter(x, y, s=18, color=OKABE_ITO[0], alpha=0.78, edgecolor="white", linewidth=0.35)
    ax.plot(x, np.polyval(coef, x), color=OKABE_ITO[1], linewidth=1.2)
    ax.set_xlabel("Age")
    ax.set_ylabel("Module score")
    ax.set_title("Linear trend")
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def generate_bubble_plot() -> Path:
    pathways = ["Complement", "ECM", "OXPHOS", "Cytokine", "Proteostasis"]
    groups = ["Plasma", "BMIF", "Tissue"]
    rng = np.random.default_rng(7)
    rows = []
    for i, p in enumerate(pathways):
        for j, g in enumerate(groups):
            rows.append((p, g, rng.uniform(-1, 1), rng.uniform(10, 95)))
    df = pd.DataFrame(rows, columns=["pathway", "group", "effect", "size"])
    out = GEN / "bubble_plot.png"
    fig, ax = plt.subplots(figsize=(4.25, 2.45))
    x = df["group"].map({g: i for i, g in enumerate(groups)})
    y = df["pathway"].map({p: i for i, p in enumerate(pathways)})
    sc = ax.scatter(
        x,
        y,
        s=df["size"] * 2.2,
        c=df["effect"],
        cmap="RdBu_r",
        vmin=-1,
        vmax=1,
        edgecolor="#303030",
        linewidth=0.25,
    )
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups)
    ax.set_yticks(range(len(pathways)))
    ax.set_yticklabels(pathways)
    ax.set_title("Bubble summary")
    ax.tick_params(length=0)
    ax.grid(color="#e6e6e6", linewidth=0.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label("Effect")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def generate_rank_plot() -> Path:
    rng = np.random.default_rng(11)
    genes = ["VSIG4", "PTX3", "PZP", "LILRB5", "LGALS3", "S100A4", "IL18BP", "SDC1"]
    vals = np.sort(rng.normal(0.36, 0.13, len(genes)))[::-1]
    out = GEN / "rank_plot.png"
    fig, ax = plt.subplots(figsize=(3.7, 2.35))
    colors = [OKABE_ITO[1] if v >= np.median(vals) else OKABE_ITO[0] for v in vals]
    ax.vlines(range(len(genes)), 0, vals, color="#8b949e", linewidth=0.8)
    ax.scatter(range(len(genes)), vals, s=26, color=colors, edgecolor="#222222", linewidth=0.35)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_ylabel("Rank score")
    ax.set_title("Ranked candidates")
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def generate_heatmap_tile() -> Path:
    rng = np.random.default_rng(13)
    mat = rng.normal(0, 1, (6, 7))
    out = GEN / "heatmap_tile.png"
    fig, ax = plt.subplots(figsize=(4.0, 2.65))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(7))
    ax.set_xticklabels(["EC", "Fib", "Mac", "VSMC", "T", "B", "Neu"], rotation=45, ha="right")
    ax.set_yticks(range(6))
    ax.set_yticklabels(["GO", "KEGG", "Reactome", "LR", "Marker", "Module"])
    ax.set_title("Tile heatmap")
    for i in range(mat.shape[0] + 1):
        ax.axhline(i - 0.5, color="white", linewidth=0.7)
    for j in range(mat.shape[1] + 1):
        ax.axvline(j - 0.5, color="white", linewidth=0.7)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("Z-score")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def generate_pca_plot() -> Path:
    rng = np.random.default_rng(17)
    centers = np.array([[-1.4, -0.2], [0.7, 0.9], [1.4, -0.65]])
    out = GEN / "pca_score.png"
    fig, ax = plt.subplots(figsize=(3.35, 2.35))
    for label, center, color in zip(["Young", "Middle", "Old"], centers, [OKABE_ITO[2], OKABE_ITO[0], OKABE_ITO[1]]):
        pts = center + rng.normal(0, 0.27, (14, 2))
        ax.scatter(pts[:, 0], pts[:, 1], s=22, color=color, alpha=0.82, label=label, edgecolor="white", linewidth=0.4)
    ax.set_xlabel("PC1 (38%)")
    ax.set_ylabel("PC2 (17%)")
    ax.set_title("PCA score plot")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    _style_axes(ax)
    fig.subplots_adjust(bottom=0.30)
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def generate_marker_dotplot() -> Path:
    rows = ["Fib", "EC", "Mac", "VSMC", "T cell", "B cell"]
    genes = ["COL1A1", "PECAM1", "LYZ", "ACTA2", "CD3D", "MS4A1", "PTX3", "VSIG4"]
    rng = np.random.default_rng(19)
    data = []
    for r_i, r in enumerate(rows):
        for g_i, g in enumerate(genes):
            base = 0.22 + 0.55 * (r_i == (g_i % len(rows)))
            avg = base + rng.normal(0, 0.08)
            pct = np.clip(8 + 72 * base + rng.normal(0, 7), 0, 100)
            data.append((r, g, avg, pct))
    df = pd.DataFrame(data, columns=["subtype", "gene", "avg_expr", "pct_expr"])
    out = GEN / "marker_dotplot.png"
    fig, _, _ = marker_dotplot(
        df,
        rows,
        genes,
        title="Marker panel",
        xlabel="Marker gene",
        ylabel="Cell type",
        fig_size=(5.2, 3.25),
        block_per_gene={
            "COL1A1": "Lineage",
            "PECAM1": "Lineage",
            "LYZ": "Lineage",
            "ACTA2": "Lineage",
            "CD3D": "Lineage",
            "MS4A1": "Lineage",
            "PTX3": "State",
            "VSIG4": "State",
        },
    )
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def generate_dumbbell() -> Path:
    df = pd.DataFrame(
        {
            "edge": ["Fib -> EC", "Mac -> Fib", "VSMC -> EC", "T cell -> Mac", "B cell -> Mac"],
            "sample_mean_young": [0.18, 0.22, 0.12, 0.09, 0.07],
            "sample_mean_old": [0.34, 0.17, 0.23, 0.05, 0.12],
            "sample_delta_old_minus_young": [0.16, -0.05, 0.11, -0.04, 0.05],
            "age_sig": ["FDR<0.05", "nominal", "FDR<0.05", "nominal", "nominal"],
        }
    )
    base = GEN / "dumbbell"
    render_mean_dumbbell(
        df,
        out_base=str(base),
        title="Young vs Old mean",
        xlabel="Mean communication strength",
        subtitle=None,
        xlim_frac=1.45,
        legend_ncol=2,
    )
    base.with_suffix(".pdf").unlink(missing_ok=True)
    base.with_suffix(".svg").unlink(missing_ok=True)
    return base.with_suffix(".png")


def generate_venn() -> Path:
    out = GEN / "venn_diagram.png"
    a = {"PTX3", "VSIG4", "PZP", "LILRB5", "LGALS3", "IL18BP"}
    b = {"PTX3", "VSIG4", "PZP", "LILRB5", "CRYAB", "CTSG"}
    c = {"PTX3", "VSIG4", "PZP", "MZB1", "BPI", "CD9"}
    fig, _ = venn_diagram(
        [a, b, c],
        ["BMIF", "Plasma", "Array"],
        mode="equal",
        title="Shared candidates",
    )
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def generate_chord() -> Path:
    out = GEN / "chord_diagram.png"
    mat = pd.DataFrame(
        [[0, 5, 1, 3], [2, 0, 4, 1], [3, 1, 0, 5], [1, 3, 2, 0]],
        index=["HSC/MPP", "Mac", "BC", "Neu"],
        columns=["HSC/MPP", "Mac", "BC", "Neu"],
    )
    fig, _ = chord_diagram(mat, figsize=(3.8, 3.8), fontsize=8, pad=4.0)
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def generate_hub_spoke() -> Path:
    out = GEN / "hub_spoke.png"
    outer = {
        "EC": {"size": 92, "color": 0.92},
        "Fib": {"size": 110, "color": 0.85},
        "Mac": {"size": 75, "color": 0.70},
        "VSMC": {"size": 66, "color": 0.55},
        "T cell": {"size": 42, "color": 0.38},
        "B cell": {"size": 38, "color": 0.30},
        "Neuron": {"size": 58, "color": 0.48},
    }
    mid = {
        "JUN": {"size": 80, "color": 0.72},
        "FOS": {"size": 68, "color": 0.62},
        "NFKB1": {"size": 60, "color": 0.52},
    }
    fig, _ = hub_spoke(
        hub="FOXO3",
        mid=mid,
        outer=outer,
        title="TF hub network",
        color_legend_label="Activity",
        size_legend=(40, 80, 110),
        figsize=(4.2, 4.2),
    )
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def generate_all() -> dict[str, Path]:
    GEN.mkdir(parents=True, exist_ok=True)
    return {
        "box_plot": generate_box_plot(),
        "scatter_linear_fit": generate_scatter_linear_fit(),
        "bubble_plot": generate_bubble_plot(),
        "rank_plot": generate_rank_plot(),
        "heatmap_tile": generate_heatmap_tile(),
        "pca_score": generate_pca_plot(),
        "marker_dotplot": generate_marker_dotplot(),
        "dumbbell": generate_dumbbell(),
        "venn_diagram": generate_venn(),
        "chord_diagram": generate_chord(),
        "hub_spoke": generate_hub_spoke(),
    }


def build_examples(generated: dict[str, Path]) -> list[Example]:
    p = ROOT / "drawing_yh"
    return [
        Example("bar-horizontal", "Quantitative", "Bar", "Horizontal bar", p / "bar/horizontal/fig.png", p / "bar/horizontal/main.py", "排序比较、top-N 指标、富集结果横向展示。", "适合标签较长的类别。", ("bar", "ranking")),
        Example("bar-circular", "Quantitative", "Bar", "Circular bar", p / "bar/circle/test.png", p / "bar/circle/main.py", "周期性或环形排列的类别值。", "视觉更强，适合少量类别。", ("bar", "radial")),
        Example("bar-enrichment", "Quantitative", "Bar", "Enrichment bar", p / "bar/enrichment_bar/figure/My_GO_KEGG_combo.png", p / "bar/enrichment_bar/enrichment_bar.py", "GO / KEGG / Reactome 富集结果。", "通常用 -log10(FDR) 或 enrichment score 排序。", ("bar", "enrichment")),
        Example("bar-waterfall", "Quantitative", "Bar", "Waterfall", p / "bar/waterfall/waterfall.png", p / "bar/waterfall/waterfall.py", "按方向和幅度排序的连续变化。", "适合展示个体、基因或通路的正负变化。", ("bar", "waterfall")),
        Example("box-plot", "Quantitative", "Distribution", "Box plot", generated["box_plot"], OUT / "build_gallery.py", "组间分布与离散程度。", "可叠加散点显示样本量。", ("box", "distribution")),
        Example("violin-plot", "Quantitative", "Distribution", "Violin plot", p / "violin_plot/violin.png", p / "violin_plot/violin_function.py", "连续变量的密度分布。", "比 box plot 更突出分布形状。", ("violin", "distribution")),
        Example("scatter-linear", "Quantitative", "Scatter", "Scatter + linear fit", generated["scatter_linear_fit"], OUT / "build_gallery.py", "年龄相关、相关性、回归趋势。", "点表示样本，线表示拟合趋势。", ("scatter", "correlation")),
        Example("bubble-plot", "Quantitative", "Scatter", "Bubble plot", generated["bubble_plot"], OUT / "build_gallery.py", "二维类别矩阵 + 第三变量大小。", "颜色和气泡大小可同时编码。", ("bubble", "matrix")),
        Example("rank-plot", "Quantitative", "Scatter", "Rank plot", generated["rank_plot"], OUT / "build_gallery.py", "候选基因、通路、feature 的排序。", "适合强调 top candidates。", ("rank", "candidate")),
        Example("dot-chart", "Quantitative", "Scatter", "Dot chart", p / "scatter/dot_chart/output/protein_dotplot/TNF_p0.png", p / "scatter/dot_chart/dot_chart.py", "蛋白/代谢物跨组织或分组点图。", "可用颜色、位置、显著性共同编码。", ("dot", "summary")),
        Example("marker-dotplot", "Omics", "Marker / QC", "Marker dot plot", generated["marker_dotplot"], p / "dotplot.py", "单细胞 marker 表达矩阵。", "颜色=gene-scaled mean，大小=pct expressed。", ("single-cell", "marker", "dotplot")),
        Example("dumbbell", "Omics", "Comparison", "Dumbbell chart", generated["dumbbell"], p / "dumbbell.py", "Young vs Old 或处理前后均值比较。", "箭头方向编码升降，适合通讯强度或 pathway score。", ("dumbbell", "comparison")),
        Example("volcano", "Omics", "Differential", "Volcano plot", p / "volcano_plot/Volcano_plot.png", p / "volcano_plot/volcano.py", "差异分析结果，logFC × p-value。", "适合快速筛选显著上/下调。", ("volcano", "DEG")),
        Example("pca", "Omics", "Marker / QC", "PCA score plot", generated["pca_score"], p / "pca/PCA.py", "样本整体结构、批次、分组分离。", "用于 QC 和组间结构展示。", ("PCA", "QC")),
        Example("heatmap-clustered", "Omics", "Heatmap", "Clustered heatmap", p / "heatmap/heatmap_clustered/HMM 20241126/heatmap_with_colorbar_cluster.png", p / "heatmap/heatmap_clustered/heapmap_clustered.py", "矩阵聚类、表达模式、样本/基因分组。", "适合 feature 数中等的全局模式。", ("heatmap", "cluster")),
        Example("heatmap-tile", "Omics", "Heatmap", "Tile heatmap", generated["heatmap_tile"], p / "heatmap/heatmap_tile_style/heatmap_tile_style.py", "通路 × 细胞类型、组织 × feature 的紧凑矩阵。", "适合报告里快速比较方向。", ("heatmap", "tile")),
        Example("dose-heatmap", "Omics", "Heatmap", "Dose-response heatmap", p / "dose_response/heatmap_with_colorbar.png", p / "dose_response/heatmap_with_colorbar.py", "药物浓度 × 组合条件矩阵。", "适合筛选敏感窗口。", ("dose-response", "heatmap")),
        Example("pie", "Composition", "Pie / Donut", "Pie chart", p / "pie_chart/pie/pie.png", p / "pie_chart/pie/main.py", "少量类别的组成比例。", "类别过多时优先换 bar plot。", ("pie", "composition")),
        Example("donut", "Composition", "Pie / Donut", "Donut chart", p / "pie_chart/donut/donut.png", p / "pie_chart/donut/main.py", "组成比例 + 中心注释。", "比普通 pie 更适合放总数或标签。", ("donut", "composition")),
        Example("explode-pie", "Composition", "Pie / Donut", "Exploded pie", p / "pie_chart/explode/explode.png", p / "pie_chart/explode/main.py", "强调某一个组成部分。", "只适合非常明确的 highlight。", ("pie", "highlight")),
        Example("zoom-pie", "Composition", "Pie / Donut", "Zoom pie", p / "pie_chart/zoom_pie/zoom_pie.png", p / "pie_chart/zoom_pie/main.py", "大类 + 小类局部放大。", "适合比例悬殊但又要展示小类。", ("pie", "zoom")),
        Example("venn", "Composition", "Set overlap", "Venn diagram", generated["venn_diagram"], p / "venn/venn_plot.py", "2-3 个集合交集。", "集合大小悬殊时可用 equal/log 模式。", ("venn", "set")),
        Example("dose-curve", "Curves", "Dose / enzyme", "Dose-response curve", p / "dose_response/fig.png", p / "dose_response/main.py", "药物浓度-效应曲线。", "适合 IC50/viability 展示。", ("curve", "dose-response")),
        Example("michaelis-menten", "Curves", "Dose / enzyme", "Michaelis-Menten curve", p / "michaelis_menten/fig.png", p / "michaelis_menten/MM_curve.py", "酶动力学 Vmax/Km 曲线。", "适合反应速率拟合。", ("curve", "enzyme")),
        Example("chord", "Network", "Flow / relation", "Chord diagram", generated["chord_diagram"], p / "chord/main.py", "有向加权网络，如细胞通讯 sender → receiver。", "ribbon 宽度表示权重，颜色跟 sender 对齐。", ("network", "chord")),
        Example("hub-spoke", "Network", "Flow / relation", "Hub-spoke network", generated["hub_spoke"], p / "network/hub_spoke.py", "中心 TF / ligand 与 target/cell type 的辐射网络。", "适合展示一个核心节点的调控范围。", ("network", "TF")),
        Example("tissue-icons", "Utility", "Assets", "Tissue / species icons", p / "icon/lib/brain.png", p / "icon/generate_missing_icons.py", "报告和组合图里的组织/物种图标。", "用于 legend、schematic 或 atlas overview。", ("icon", "utility")),
    ]


def render_html(examples: list[Example]) -> str:
    categories: list[str] = []
    sub_categories: dict[str, list[str]] = {}
    examples_by_sub: dict[tuple[str, str], list[Example]] = {}
    for ex in examples:
        if ex.category not in categories:
            categories.append(ex.category)
        sub_categories.setdefault(ex.category, [])
        if ex.sub_category not in sub_categories[ex.category]:
            sub_categories[ex.category].append(ex.sub_category)
        examples_by_sub.setdefault((ex.category, ex.sub_category), []).append(ex)

    top_anchor = {cat: _anchor("top", cat) for cat in categories}
    sub_anchor = {
        (cat, sub): _anchor("sub", f"{cat}-{sub}")
        for cat, subs in sub_categories.items()
        for sub in subs
    }
    all_subs = [(cat, sub) for cat in categories for sub in sub_categories[cat]]

    top_buttons = "\n".join(
        f'<button class="tab top-tab" data-top="{html.escape(cat)}" data-anchor="{html.escape(top_anchor[cat])}">'
        f'{html.escape(cat)}<span>{sum(1 for ex in examples if ex.category == cat)}</span></button>'
        for cat in categories
    )
    sub_buttons = "\n".join(
        f'<button class="tab sub-tab" data-top="{html.escape(cat)}" data-sub-id="{html.escape(sub_anchor[(cat, sub)])}" data-anchor="{html.escape(sub_anchor[(cat, sub)])}">'
        f'{html.escape(sub)}<span>{len(examples_by_sub[(cat, sub)])}</span></button>'
        for cat, sub in all_subs
    )

    sidebar_parts: list[str] = []
    for cat in categories:
        sidebar_parts.append(
            f'<div class="sidebar-group" data-top="{html.escape(cat)}">'
            f'<a class="sidebar-title" href="#{html.escape(top_anchor[cat])}" data-action="top" data-top="{html.escape(cat)}">{html.escape(cat)}</a>'
        )
        for sub in sub_categories[cat]:
            sid = sub_anchor[(cat, sub)]
            sidebar_parts.append(
                f'<div class="sidebar-sub" data-top="{html.escape(cat)}" data-sub-id="{html.escape(sid)}">'
                f'<a class="sidebar-subtitle" href="#{html.escape(sid)}" data-action="sub" data-top="{html.escape(cat)}" data-sub-id="{html.escape(sid)}">{html.escape(sub)}</a>'
            )
            for ex in examples_by_sub[(cat, sub)]:
                search_text = html.escape((ex.title + " " + ex.use_case + " " + " ".join(ex.tags)).lower())
                sidebar_parts.append(
                    f'<a class="sidebar-leaf" href="#{html.escape(ex.slug)}" data-action="card" data-top="{html.escape(cat)}" '
                    f'data-sub-id="{html.escape(sid)}" data-card="{html.escape(ex.slug)}" data-search="{search_text}">{html.escape(ex.title)}</a>'
                )
            sidebar_parts.append("</div>")
        sidebar_parts.append("</div>")

    def render_card(ex: Example, sid: str) -> str:
        img_url = _rel_url(ex.image)
        source_url = _rel_url(ex.source)
        root_rel = html.escape(_root_rel(ex.source))
        tags = " ".join(f"<span>{html.escape(t)}</span>" for t in ex.tags)
        search_text = html.escape((ex.title + " " + ex.use_case + " " + ex.note + " " + " ".join(ex.tags)).lower())
        return f"""
<article class="card" id="{html.escape(ex.slug)}" data-top="{html.escape(ex.category)}" data-sub-id="{html.escape(sid)}" data-search="{search_text}">
  <a class="thumb" href="{img_url}" title="Open image">
    <img src="{img_url}" alt="{html.escape(ex.title)}">
  </a>
  <div class="card-body">
    <div class="meta">
      <span class="badge">{html.escape(ex.category)}</span>
      <span class="badge sub">{html.escape(ex.sub_category)}</span>
      <span class="source">{root_rel}</span>
    </div>
    <h4>{html.escape(ex.title)}</h4>
    <p>{html.escape(ex.use_case)}</p>
    <p class="note">{html.escape(ex.note)}</p>
    <div class="tags">{tags}</div>
    <div class="links"><a href="{img_url}">原图</a><a href="{source_url}">源码/示例</a></div>
  </div>
</article>"""

    grouped: list[str] = []
    for cat in categories:
        grouped.append(
            f'<section class="top-section" id="{html.escape(top_anchor[cat])}" data-top="{html.escape(cat)}">'
            f'<div class="section-head"><h2>{html.escape(cat)}</h2><span>{sum(1 for ex in examples if ex.category == cat)} examples</span></div>'
        )
        for sub in sub_categories[cat]:
            sid = sub_anchor[(cat, sub)]
            grouped.append(
                f'<section class="sub-section" id="{html.escape(sid)}" data-top="{html.escape(cat)}" data-sub-id="{html.escape(sid)}">'
                f'<div class="sub-head"><h3>{html.escape(sub)}</h3><span>{len(examples_by_sub[(cat, sub)])}</span></div>'
                '<div class="grid">'
            )
            grouped.extend(render_card(ex, sid) for ex in examples_by_sub[(cat, sub)])
            grouped.append("</div></section>")
        grouped.append("</section>")

    legacy_hash = {
        "Bar": {"top": "Quantitative", "subId": sub_anchor[("Quantitative", "Bar")], "target": sub_anchor[("Quantitative", "Bar")]},
        "Distribution": {"top": "Quantitative", "subId": sub_anchor[("Quantitative", "Distribution")], "target": sub_anchor[("Quantitative", "Distribution")]},
        "Scatter": {"top": "Quantitative", "subId": sub_anchor[("Quantitative", "Scatter")], "target": sub_anchor[("Quantitative", "Scatter")]},
        "Omics": {"top": "Omics", "subId": "all", "target": top_anchor["Omics"]},
        "Heatmap": {"top": "Omics", "subId": sub_anchor[("Omics", "Heatmap")], "target": sub_anchor[("Omics", "Heatmap")]},
        "Composition": {"top": "Composition", "subId": "all", "target": top_anchor["Composition"]},
        "Curve": {"top": "Curves", "subId": "all", "target": top_anchor["Curves"]},
        "Curves": {"top": "Curves", "subId": "all", "target": top_anchor["Curves"]},
        "Network": {"top": "Network", "subId": "all", "target": top_anchor["Network"]},
        "Utility": {"top": "Utility", "subId": "all", "target": top_anchor["Utility"]},
    }
    legacy_json = json.dumps(legacy_hash, ensure_ascii=False)

    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>drawing-yh 图型示例库</title>
  <style>
    :root {{
      --bg: #f6f7f4;
      --panel: #ffffff;
      --ink: #1d2528;
      --muted: #667174;
      --line: #d9ded8;
      --teal: #237b73;
      --teal-soft: #dcebe7;
      --amber: #a76517;
      --blue: #315f8c;
      --violet: #6b5b95;
    }}
    * {{ box-sizing: border-box; }}
    html {{ scroll-behavior: smooth; }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font: 15px/1.55 -apple-system, BlinkMacSystemFont, "Segoe UI", "Microsoft YaHei", Arial, sans-serif;
      letter-spacing: 0;
    }}
    .layout {{ display: grid; grid-template-columns: 260px minmax(0, 1fr); min-height: 100vh; }}
    aside {{
      position: sticky;
      top: 0;
      height: 100vh;
      padding: 22px 16px;
      background: #eef1ec;
      border-right: 1px solid var(--line);
      overflow-y: auto;
    }}
    aside h1 {{ font-size: 18px; line-height: 1.25; margin: 0 0 14px; }}
    aside .count {{ color: var(--muted); font-size: 13px; margin-bottom: 18px; }}
    aside a {{
      display: block;
      color: var(--ink);
      text-decoration: none;
    }}
    .sidebar-title {{
      margin-top: 12px;
      padding: 7px 8px;
      border-radius: 6px;
      font-weight: 700;
    }}
    .sidebar-subtitle {{
      margin: 4px 0 2px 10px;
      padding: 5px 8px;
      color: var(--blue);
      font-weight: 600;
      font-size: 13px;
    }}
    .sidebar-leaf {{
      margin-left: 22px;
      padding: 4px 8px;
      font-size: 13px;
      color: var(--ink);
      border-left: 3px solid transparent;
    }}
    aside a:hover {{ color: var(--teal); background: #e6ebe5; }}
    .sidebar-leaf:hover {{ border-left-color: var(--teal); }}
    main {{ padding: 24px 34px 56px; max-width: 1480px; }}
    .toolbar {{
      position: sticky;
      top: 0;
      z-index: 5;
      display: grid;
      gap: 10px;
      padding: 10px 0 18px;
      background: linear-gradient(var(--bg) 82%, rgba(246,247,244,0));
    }}
    .toolbar-top {{
      display: flex;
      gap: 12px;
      align-items: center;
      justify-content: space-between;
      flex-wrap: wrap;
    }}
    .nav-row {{ display: flex; gap: 8px; align-items: center; flex-wrap: wrap; min-height: 36px; }}
    .summary {{ color: var(--muted); font-size: 13px; }}
    input {{
      width: min(360px, 100%);
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 9px 11px;
      background: white;
      color: var(--ink);
      font: inherit;
    }}
    .tab {{
      border: 1px solid var(--line);
      background: white;
      color: var(--ink);
      border-radius: 6px;
      padding: 8px 11px;
      cursor: pointer;
      font: inherit;
      min-height: 34px;
    }}
    .tab span {{ margin-left: 6px; color: var(--muted); font-size: 12px; }}
    .tab.active {{ background: var(--teal); border-color: var(--teal); color: white; }}
    .tab.active span {{ color: rgba(255,255,255,0.78); }}
    .sub-tab.active {{ background: var(--blue); border-color: var(--blue); }}
    .top-section {{ scroll-margin-top: 130px; margin-top: 22px; }}
    .sub-section {{ scroll-margin-top: 130px; margin-top: 16px; }}
    .section-head {{
      display: flex;
      justify-content: space-between;
      gap: 12px;
      align-items: baseline;
      border-bottom: 1px solid var(--line);
      padding-bottom: 8px;
      margin: 28px 0 14px;
    }}
    .section-head h2 {{
      margin: 0;
      font-size: 20px;
    }}
    .section-head span, .sub-head span {{ color: var(--muted); font-size: 13px; }}
    .sub-head {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      margin: 8px 0 10px;
    }}
    .sub-head h3 {{ margin: 0; font-size: 16px; color: var(--blue); }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 16px; }}
    .card {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      overflow: hidden;
      min-height: 100%;
      display: flex;
      flex-direction: column;
      box-shadow: 0 1px 2px rgba(20, 34, 36, 0.05);
    }}
    .thumb {{
      display: grid;
      place-items: center;
      height: 245px;
      padding: 10px;
      background: #fbfbf8;
      border-bottom: 1px solid var(--line);
    }}
    .thumb img {{
      max-width: 100%;
      max-height: 225px;
      object-fit: contain;
      background: white;
    }}
    .card-body {{ padding: 13px 14px 15px; }}
    .meta {{ display: flex; gap: 8px; align-items: center; min-width: 0; }}
    .badge {{
      display: inline-block;
      color: white;
      background: var(--blue);
      border-radius: 5px;
      padding: 2px 7px;
      font-size: 12px;
      line-height: 1.45;
      flex: 0 0 auto;
    }}
    .badge.sub {{ background: var(--violet); }}
    .source {{
      color: var(--muted);
      font-size: 12px;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }}
    h4 {{ margin: 9px 0 6px; font-size: 17px; }}
    p {{ margin: 6px 0; }}
    .note {{ color: var(--muted); }}
    .tags {{ display: flex; flex-wrap: wrap; gap: 5px; margin: 10px 0; }}
    .tags span {{
      border: 1px solid #e1d2c2;
      color: var(--amber);
      border-radius: 4px;
      padding: 1px 6px;
      font-size: 12px;
      background: #fff9f1;
    }}
    .links {{ display: flex; gap: 14px; margin-top: 8px; }}
    .links a {{ color: var(--teal); text-decoration: none; font-weight: 600; }}
    .links a:hover {{ text-decoration: underline; }}
    .no-results {{
      border: 1px dashed var(--line);
      color: var(--muted);
      padding: 18px;
      border-radius: 8px;
      background: rgba(255,255,255,0.55);
      margin-top: 24px;
    }}
    .hidden {{ display: none !important; }}
    @media (max-width: 820px) {{
      .layout {{ grid-template-columns: 1fr; }}
      aside {{ position: relative; height: auto; border-right: 0; border-bottom: 1px solid var(--line); }}
      aside nav {{ display: flex; flex-wrap: wrap; gap: 4px; }}
      .sidebar-sub, .sidebar-group {{ display: contents; }}
      .sidebar-title, .sidebar-subtitle, .sidebar-leaf {{ margin: 0; border-left: 0; border-bottom: 2px solid transparent; }}
      .sidebar-leaf {{ display: none; }}
      main {{ padding: 18px 16px 42px; }}
      .toolbar {{ position: relative; }}
      .toolbar-top {{ align-items: stretch; }}
      .thumb {{ height: 220px; }}
    }}
  </style>
</head>
<body>
  <div class="layout">
    <aside>
      <h1>drawing-yh 图型示例库</h1>
      <div class="count">{len(examples)} 种示例 · 2026-05-28</div>
      <nav>{"".join(sidebar_parts)}</nav>
    </aside>
    <main>
      <div class="toolbar">
        <div class="toolbar-top">
          <input id="search" type="search" placeholder="搜索图型 / 场景 / tag">
          <div class="summary"><span id="visibleCount">{len(examples)}</span> / {len(examples)} examples</div>
        </div>
        <div class="nav-row top-nav">
          <button class="tab top-tab active" data-top="all" data-anchor="">全部</button>
          {top_buttons}
        </div>
        <div class="nav-row sub-nav">
          <button class="tab sub-tab active" data-top="all" data-sub-id="all" data-anchor="">全部二级</button>
          {sub_buttons}
        </div>
      </div>
      {"".join(grouped)}
      <div id="noResults" class="no-results hidden">没有匹配图型</div>
    </main>
  </div>
  <script>
    const search = document.querySelector("#search");
    const topTabs = Array.from(document.querySelectorAll(".top-tab"));
    const subTabs = Array.from(document.querySelectorAll(".sub-tab"));
    const cards = Array.from(document.querySelectorAll(".card"));
    const topSections = Array.from(document.querySelectorAll(".top-section"));
    const subSections = Array.from(document.querySelectorAll(".sub-section"));
    const sidebarGroups = Array.from(document.querySelectorAll(".sidebar-group"));
    const sidebarSubs = Array.from(document.querySelectorAll(".sidebar-sub"));
    const sidebarLeaves = Array.from(document.querySelectorAll(".sidebar-leaf"));
    const visibleCount = document.querySelector("#visibleCount");
    const noResults = document.querySelector("#noResults");
    const legacyMap = {legacy_json};
    let activeTop = "all";
    let activeSub = "all";

    function setHash(anchor) {{
      if (anchor) {{
        history.replaceState(null, "", "#" + anchor);
      }} else {{
        history.replaceState(null, "", location.pathname + location.search);
      }}
    }}

    function syncTabs() {{
      topTabs.forEach(btn => btn.classList.toggle("active", btn.dataset.top === activeTop));
      subTabs.forEach(btn => {{
        const isAll = btn.dataset.subId === "all";
        const belongs = activeTop === "all" || btn.dataset.top === activeTop;
        btn.classList.toggle("hidden", !isAll && !belongs);
        btn.classList.toggle("active", btn.dataset.subId === activeSub);
      }});
    }}

    function applyFilter() {{
      const q = search.value.trim().toLowerCase();
      let visible = 0;
      cards.forEach(card => {{
        const okTop = activeTop === "all" || card.dataset.top === activeTop;
        const okSub = activeSub === "all" || card.dataset.subId === activeSub;
        const okText = !q || card.dataset.search.includes(q);
        const show = okTop && okSub && okText;
        card.classList.toggle("hidden", !show);
        if (show) visible += 1;
      }});
      subSections.forEach(sec => {{
        const visible = Array.from(sec.querySelectorAll(".card")).some(c => !c.classList.contains("hidden"));
        sec.classList.toggle("hidden", !visible);
      }});
      topSections.forEach(sec => {{
        const visible = Array.from(sec.querySelectorAll(".card")).some(c => !c.classList.contains("hidden"));
        sec.classList.toggle("hidden", !visible);
      }});
      sidebarLeaves.forEach(leaf => {{
        const okTop = activeTop === "all" || leaf.dataset.top === activeTop;
        const okSub = activeSub === "all" || leaf.dataset.subId === activeSub;
        const okText = !q || leaf.dataset.search.includes(q);
        leaf.classList.toggle("hidden", !(okTop && okSub && okText));
      }});
      sidebarSubs.forEach(sec => {{
        const anyLeaf = Array.from(sec.querySelectorAll(".sidebar-leaf")).some(x => !x.classList.contains("hidden"));
        sec.classList.toggle("hidden", !anyLeaf);
      }});
      sidebarGroups.forEach(sec => {{
        const anyLeaf = Array.from(sec.querySelectorAll(".sidebar-leaf")).some(x => !x.classList.contains("hidden"));
        sec.classList.toggle("hidden", !anyLeaf);
      }});
      visibleCount.textContent = visible;
      noResults.classList.toggle("hidden", visible !== 0);
      syncTabs();
    }}

    function scrollToAnchor(anchor) {{
      if (!anchor) return;
      const target = document.getElementById(anchor);
      if (target) requestAnimationFrame(() => target.scrollIntoView({{ block: "start" }}));
    }}

    function setState(top, subId, anchor) {{
      activeTop = top || "all";
      activeSub = subId || "all";
      applyFilter();
      setHash(anchor);
      scrollToAnchor(anchor);
    }}

    topTabs.forEach(btn => btn.addEventListener("click", () => {{
      setState(btn.dataset.top, "all", btn.dataset.anchor);
    }}));
    subTabs.forEach(btn => btn.addEventListener("click", () => {{
      const top = btn.dataset.subId === "all" ? activeTop : btn.dataset.top;
      setState(top, btn.dataset.subId, btn.dataset.anchor);
    }}));
    document.querySelectorAll("[data-action]").forEach(link => link.addEventListener("click", event => {{
      event.preventDefault();
      if (link.dataset.action === "top") {{
        setState(link.dataset.top, "all", link.getAttribute("href").slice(1));
      }} else if (link.dataset.action === "sub") {{
        setState(link.dataset.top, link.dataset.subId, link.getAttribute("href").slice(1));
      }} else {{
        setState(link.dataset.top, link.dataset.subId, link.dataset.card);
      }}
    }}));
    search.addEventListener("input", applyFilter);

    function restoreHash() {{
      const raw = decodeURIComponent(location.hash.replace(/^#/, ""));
      if (!raw) {{
        applyFilter();
        return;
      }}
      if (legacyMap[raw]) {{
        const item = legacyMap[raw];
        activeTop = item.top;
        activeSub = item.subId;
        applyFilter();
        scrollToAnchor(item.target);
        return;
      }}
      const target = document.getElementById(raw);
      if (target && target.classList.contains("card")) {{
        activeTop = target.dataset.top;
        activeSub = target.dataset.subId;
      }} else if (target && target.classList.contains("sub-section")) {{
        activeTop = target.dataset.top;
        activeSub = target.dataset.subId;
      }} else if (target && target.classList.contains("top-section")) {{
        activeTop = target.dataset.top;
        activeSub = "all";
      }}
      applyFilter();
      scrollToAnchor(raw);
    }}

    restoreHash();
    window.addEventListener("hashchange", restoreHash);
  </script>
</body>
</html>
"""


def main() -> None:
    generated = generate_all()
    examples = build_examples(generated)
    missing = [ex for ex in examples if not ex.image.exists()]
    if missing:
        miss = "\n".join(f"- {ex.slug}: {ex.image}" for ex in missing)
        raise FileNotFoundError(f"Missing gallery images:\n{miss}")
    html_text = render_html(examples)
    (OUT / "index.html").write_text(html_text, encoding="utf-8")
    manifest = pd.DataFrame(
        [
            {
                "slug": ex.slug,
                "category": ex.category,
                "sub_category": ex.sub_category,
                "title": ex.title,
                "image": _root_rel(ex.image),
                "source": _root_rel(ex.source),
                "use_case": ex.use_case,
                "note": ex.note,
                "tags": ",".join(ex.tags),
            }
            for ex in examples
        ]
    )
    manifest.to_csv(OUT / "gallery_manifest.csv", index=False, encoding="utf-8-sig")
    print(f"Wrote {OUT / 'index.html'}")
    print(f"Wrote {OUT / 'gallery_manifest.csv'}")
    print(f"Examples: {len(examples)}")


if __name__ == "__main__":
    main()
