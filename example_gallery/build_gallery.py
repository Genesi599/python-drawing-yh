#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build a static gallery of representative drawing-yh chart examples."""
from __future__ import annotations

import html
import json
import math
import os
import shutil
from dataclasses import dataclass, replace
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
    report,
    render_mean_dumbbell,
    save_fig,
    venn_diagram,
)
from drawing_yh.chord import chord_diagram, lr_role_chord_panel
from drawing_yh.network import hub_spoke


ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
GEN = OUT / "generated"
PAGES = OUT / "pages"
IMAGE_ASSETS = OUT / "assets" / "gallery_images"
REPO_BLOB_BASE = "https://github.com/Genesi599/python-drawing-yh/blob/main"


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
        row_colors=dict(zip(rows, OKABE_ITO)),
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


def generate_lr_chord_panel() -> Path:
    out = GEN / "lr_chord_panel.png"
    cell_colors = {
        "Neuron": "#4C6EDB",
        "Glia": "#58B368",
        "Fib": "#E07A5F",
        "Myeloid": "#9D4EDD",
        "VSMC": "#A97142",
    }
    base_edges = [
        ("NRXN3", "NLGN1", "Neuron", 1, 0.18, "Old > Young", "Rescued", "Rescued"),
        ("COL4A3", "ITGAV", "Fib", 3, 0.14, "Old > Young", "Rescued", "Not rescued"),
        ("NRG3", "ERBB4", "Myeloid", 4, 0.11, "Old < Young", "Rescued", "Not rescued"),
        ("LAMA2", "ITGB8", "Glia", 2, 0.10, "Old > Young", "Not rescued", "Not rescued"),
        ("TENM2", "ADGRL3", "VSMC", 5, 0.08, "Old > Young", "Not rescued", "Rescued"),
    ]
    rows = []
    for source, target, other_cell, cell_no, weight, direction, vc_status, met_status in base_edges:
        rows.extend([
            ("Aging change", source, target, other_cell, cell_no, weight, direction, "Aging"),
            ("VC rescue", source, target, other_cell, cell_no, weight, direction, vc_status),
            ("Met rescue", source, target, other_cell, cell_no, weight, direction, met_status),
        ])
    df = pd.DataFrame(
        rows,
        columns=["panel", "source", "target", "other_cell", "cell_no", "weight", "direction", "status"],
    )
    df["edge_color_base"] = df["other_cell"].map(cell_colors)
    df["edge_color"] = np.where(df["status"].eq("Not rescued"), "#D0D0D0", df["edge_color_base"])
    df["edge_alpha"] = np.where(df["status"].eq("Not rescued"), 0.28, 0.86)
    df["number"] = np.where(df["status"].eq("Not rescued"), "", df["cell_no"].astype(str))
    df["number_color"] = np.where(df["direction"].eq("Old > Young"), "#E67E22", "#1F77B4")
    fig, _, _ = lr_role_chord_panel(
        df,
        panel_order=["Aging change", "VC rescue", "Met rescue"],
        edge_label="other_cell",
        edge_color="edge_color",
        edge_alpha="edge_alpha",
        edge_number="number",
        edge_number_color="number_color",
        edge_legend_color="edge_color_base",
        edge_legend_title="Other cell number",
        edge_legend_order=list(cell_colors),
        status="status",
        status_palette={"Rescued": "#3A8E3A", "Not rescued": "#D0D0D0"},
        number_color_legend={"Old > Young": "#E67E22", "Old < Young": "#1F77B4"},
        source_role_label="Ligand / sender",
        target_role_label="Receptor / receiver",
        source_group_label="Ligands (sender cells)",
        target_group_label="Receptors (receiver cells)",
        group_label_fontsize=8,
        group_label_radius=1.54,
        legend_title="Legend",
        figsize=(12.4, 3.8),
        chord_size=3.8,
        legend_width_ratio=0.58,
        fontsize=7,
        panel_title_fontsize=9,
        legend_fontsize=7,
        legend_heading_fontsize=8,
        legend_title_fontsize=10,
        edge_number_jitter=0.030,
        width=0.075,
        pad=2.4,
        wspace=0.0,
    )
    fig.savefig(out, dpi=220, bbox_inches="tight")
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


def generate_feature_plot() -> Path:
    from drawing_yh import feature_plot

    rng = np.random.default_rng(23)
    per = 200
    centers = np.array([[-2.6, 0.1], [2.1, 2.0], [1.4, -2.3]])
    coords = np.vstack([c + rng.normal(0, 0.7, (per, 2)) for c in centers])

    def marker(cluster_idx: int, strength: float = 3.0):
        v = rng.gamma(0.3, size=coords.shape[0])
        lo, hi = cluster_idx * per, (cluster_idx + 1) * per
        v[lo:hi] += rng.gamma(strength, size=per)
        return v

    values = {"PECAM1": marker(0), "COL1A1": marker(1), "LYZ": marker(2)}
    out = GEN / "feature_plot.png"
    fig, _ = feature_plot(coords, values, panel_size=1.95, point_size=7,
                          vmin="p2", vmax="p98")
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return out


def generate_heatmap_row_bars() -> Path:
    from drawing_yh import heatmap_with_row_bars

    rng = np.random.default_rng(29)
    rows = ["Fib", "EC", "Mac", "VSMC", "T cell", "B cell"]
    per = 3
    G = len(rows)
    Z = rng.normal(0, 0.6, (G, G * per))
    for i in range(G):
        Z[i, i * per:(i + 1) * per] += 2.3
    row_bars = {
        "Fib": [("extracellular matrix organization", 6.1), ("collagen fibril", 3.2)],
        "EC": [("vasculature development", 5.4), ("angiogenesis", 2.8)],
        "Mac": [("immune response", 4.9), ("phagocytosis", 2.4)],
        "VSMC": [("muscle contraction", 4.1)],
        "T cell": [("T cell activation", 3.8)],
        "B cell": [("B cell receptor signaling", 3.3)],
    }
    colors = dict(zip(rows, OKABE_ITO))
    out = GEN / "heatmap_row_bars.png"
    fig, _ = heatmap_with_row_bars(
        Z, row_labels=rows, row_bars=row_bars,
        block_sizes=[per] * G, row_colors=colors,
        z_clip=2.5, figsize=(6.2, 3.0),
    )
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return out


def generate_oydeg_enrichment_heatmap() -> Path:
    from drawing_yh import plot_oydeg_heatmap_enrichment

    columns = ["Fib", "EC", "Mac"]
    colors = dict(zip(columns, [OKABE_ITO[0], OKABE_ITO[2], OKABE_ITO[1]]))
    genes = [
        "COL1A1", "COL3A1", "POSTN", "MMP2", "CXCL12",
        "PECAM1", "VWF", "KDR", "SOX18", "CLDN5",
        "C1QA", "LST1", "TYROBP", "APOE", "MERTK",
        "DCN", "LUM", "RGS5", "FLT1", "CSF1R",
    ]
    values = pd.DataFrame(0.0, index=genes, columns=columns)
    values.loc[genes[0:5], "Fib"] = [1.8, 1.5, 1.2, -1.1, -0.8]
    values.loc[genes[5:10], "EC"] = [1.6, 1.4, 1.1, -1.2, -0.9]
    values.loc[genes[10:15], "Mac"] = [1.7, 1.3, 1.1, -1.0, -0.8]
    values.loc["DCN", ["Fib", "EC"]] = [0.9, 0.4]
    values.loc["FLT1", ["EC", "Mac"]] = [0.8, 0.3]
    values.loc["CSF1R", "Mac"] = 0.9

    heat_blocks = [
        {"label": "Fib", "color": colors["Fib"], "s": 0, "e": 5},
        {"label": "EC", "color": colors["EC"], "s": 5, "e": 10},
        {"label": "Mac", "color": colors["Mac"], "s": 10, "e": 15},
        {"label": "Shared", "color": "#8b949e", "s": 15, "e": 20},
    ]
    gene_status = {}
    for gene in genes:
        gene_status[gene] = {}
        for col in columns:
            val = values.loc[gene, col]
            if val > 0:
                gene_status[gene][col] = "up"
            elif val < 0:
                gene_status[gene][col] = "down"
    right_blocks = [
        {
            "label": "Fib",
            "color": colors["Fib"],
            "groups": [
                ("Up in old", "#C0392B", [
                    {"term": "Extracellular Matrix Organization", "nlp": 3.4, "db": "G", "genes": ["COL1A1", "COL3A1", "POSTN", "MMP2"]},
                ]),
                ("Down in old", "#2166AC", [
                    {"term": "Chemokine Signaling", "nlp": 1.8, "db": "K", "genes": ["CXCL12"]},
                ]),
            ],
        },
        {
            "label": "EC",
            "color": colors["EC"],
            "groups": [
                ("Up in old", "#C0392B", [
                    {"term": "Vasculature Development", "nlp": 2.8, "db": "G", "genes": ["PECAM1", "VWF", "KDR", "SOX18"]},
                ]),
            ],
        },
        {
            "label": "Mac",
            "color": colors["Mac"],
            "groups": [
                ("Up in old", "#C0392B", [
                    {"term": "Phagosome", "nlp": 2.5, "db": "K", "genes": ["C1QA", "TYROBP", "APOE", "MERTK"]},
                ]),
            ],
        },
    ]
    out = GEN / "oydeg_enrichment_heatmap.png"
    fig, _ = plot_oydeg_heatmap_enrichment(
        values,
        heat_blocks=heat_blocks,
        label_genes=["COL1A1", "POSTN", "CXCL12", "PECAM1", "KDR", "CLDN5", "C1QA", "TYROBP", "APOE", "FLT1"],
        gene_status=gene_status,
        column_colors=colors,
        right_blocks=right_blocks,
        heatmap_title="Top aging DEGs (log2FC)",
        right_title="KEGG + GO enrichment",
        figsize=(6.9, 5.2),
    )
    fig.savefig(out, dpi=220, bbox_inches="tight")
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
        "feature_plot": generate_feature_plot(),
        "heatmap_row_bars": generate_heatmap_row_bars(),
        "oydeg_enrichment_heatmap": generate_oydeg_enrichment_heatmap(),
        "dumbbell": generate_dumbbell(),
        "venn_diagram": generate_venn(),
        "chord_diagram": generate_chord(),
        "lr_chord_panel": generate_lr_chord_panel(),
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
        Example("dot-chart", "Quantitative", "Summary dot", "Dot chart", p / "scatter/dot_chart/output/protein_dotplot/TNF_p0.png", p / "scatter/dot_chart/dot_chart.py", "蛋白/代谢物跨组织或分组点图。", "可用颜色、位置、显著性共同编码。", ("dot", "summary")),
        Example("marker-dotplot", "Single-cell", "Marker expression", "Marker dot plot", generated["marker_dotplot"], p / "dotplot.py", "单细胞 marker 表达矩阵。", "颜色=gene-scaled mean(grey-red，默认 scseq 配色)，大小=pct expressed，左侧彩点=cell-type 颜色。", ("single-cell", "marker", "dotplot")),
        Example("feature-plot", "Single-cell", "Embedding", "Embedding feature plot", generated["feature_plot"], p / "embedding.py", "基因表达投影到 UMAP/tSNE 的多基因网格。", "共享 grey-red colorbar + 角落 UMAP 箭头轴，按非零百分位裁剪。", ("single-cell", "umap", "feature")),
        Example("dumbbell", "Quantitative", "Comparison", "Dumbbell chart", generated["dumbbell"], p / "dumbbell.py", "Young vs Old 或处理前后均值比较。", "箭头方向编码升降，适合通讯强度或 pathway score。", ("dumbbell", "comparison")),
        Example("volcano", "Omics", "Differential", "Volcano plot", p / "volcano_plot/Volcano_plot.png", p / "volcano_plot/volcano.py", "差异分析结果，logFC × p-value。", "适合快速筛选显著上/下调。", ("volcano", "DEG")),
        Example("pca", "QC", "Dimension reduction", "PCA score plot", generated["pca_score"], p / "pca/PCA.py", "样本整体结构、批次、分组分离。", "用于 QC 和组间结构展示。", ("PCA", "QC")),
        Example("heatmap-clustered", "Omics", "Heatmap", "Clustered heatmap", p / "heatmap/heatmap_clustered/HMM 20241126/heatmap_with_colorbar_cluster.png", p / "heatmap/heatmap_clustered/heapmap_clustered.py", "矩阵聚类、表达模式、样本/基因分组。", "适合 feature 数中等的全局模式。", ("heatmap", "cluster")),
        Example("heatmap-tile", "Omics", "Heatmap", "Tile heatmap", generated["heatmap_tile"], p / "heatmap/heatmap_tile_style/heatmap_tile_style.py", "通路 × 细胞类型、组织 × feature 的紧凑矩阵。", "适合报告里快速比较方向。", ("heatmap", "tile")),
        Example("heatmap-row-bars", "Omics", "Heatmap", "Heatmap + per-row bars", generated["heatmap_row_bars"], p / "heatmap/row_bars.py", "z-score 热图 + 每行 top 富集条(如各 cell type GO)。", "左侧行彩条 + 列块分隔，右侧横向条标注通路。", ("single-cell", "heatmap", "enrichment")),
        Example("oydeg-enrichment-heatmap", "Omics", "Heatmap", "OY-DEG heatmap + enrichment", generated["oydeg_enrichment_heatmap"], p / "heatmap/oydeg_enrichment.py", "Old vs Young DEG log2FC 热图 + 动态 gene labels + Up/Down 富集块。", "输入为已整理矩阵、block、label genes 和 enrichment blocks；DEG/enrichment 计算留在分析流程。", ("single-cell", "DEG", "enrichment")),
        Example("dose-heatmap", "Curves", "Dose response", "Dose-response heatmap", p / "dose_response/heatmap_with_colorbar.png", p / "dose_response/heatmap_with_colorbar.py", "药物浓度 × 组合条件矩阵。", "适合筛选敏感窗口。", ("dose-response", "heatmap")),
        Example("pie", "Composition", "Pie / Donut", "Pie chart", p / "pie_chart/pie/pie.png", p / "pie_chart/pie/main.py", "少量类别的组成比例。", "类别过多时优先换 bar plot。", ("pie", "composition")),
        Example("donut", "Composition", "Pie / Donut", "Donut chart", p / "pie_chart/donut/donut.png", p / "pie_chart/donut/main.py", "组成比例 + 中心注释。", "比普通 pie 更适合放总数或标签。", ("donut", "composition")),
        Example("explode-pie", "Composition", "Pie / Donut", "Exploded pie", p / "pie_chart/explode/explode.png", p / "pie_chart/explode/main.py", "强调某一个组成部分。", "只适合非常明确的 highlight。", ("pie", "highlight")),
        Example("zoom-pie", "Composition", "Pie / Donut", "Zoom pie", p / "pie_chart/zoom_pie/zoom_pie.png", p / "pie_chart/zoom_pie/main.py", "大类 + 小类局部放大。", "适合比例悬殊但又要展示小类。", ("pie", "zoom")),
        Example("venn", "Composition", "Set overlap", "Venn diagram", generated["venn_diagram"], p / "venn/venn_plot.py", "2-3 个集合交集。", "集合大小悬殊时可用 equal/log 模式。", ("venn", "set")),
        Example("dose-curve", "Curves", "Dose / enzyme", "Dose-response curve", p / "dose_response/fig.png", p / "dose_response/main.py", "药物浓度-效应曲线。", "适合 IC50/viability 展示。", ("curve", "dose-response")),
        Example("michaelis-menten", "Curves", "Dose / enzyme", "Michaelis-Menten curve", p / "michaelis_menten/fig.png", p / "michaelis_menten/MM_curve.py", "酶动力学 Vmax/Km 曲线。", "适合反应速率拟合。", ("curve", "enzyme")),
        Example("chord", "Network", "Flow / relation", "Chord diagram", generated["chord_diagram"], p / "chord/main.py", "有向加权网络，如细胞通讯 sender → receiver。", "ribbon 宽度表示权重，颜色跟 sender 对齐。", ("network", "chord")),
        Example("lr-chord-panel", "Network", "Flow / relation", "LR chord panel", generated["lr_chord_panel"], p / "chord/main.py", "配体-受体弦图的横向多面板比较。", "共享 legend、可选 rescue 灰边、边编号和方向色，适合 aging / intervention 对照。", ("network", "chord", "ligand-receptor", "panel")),
        Example("hub-spoke", "Network", "Flow / relation", "Hub-spoke network", generated["hub_spoke"], p / "network/hub_spoke.py", "中心 TF / ligand 与 target/cell type 的辐射网络。", "适合展示一个核心节点的调控范围。", ("network", "TF")),
        Example("tissue-icons", "Utility", "Assets", "Tissue / species icons", p / "icon/lib/brain.png", p / "icon/generate_missing_icons.py", "报告和组合图里的组织/物种图标。", "用于 legend、schematic 或 atlas overview。", ("icon", "utility")),
        Example("flat-cell-schematic", "Utility", "Assets", "Flat cell schematic", p / "cell_schematic/templates/flat_eukaryotic_cell_preview.png", p / "cell_schematic/__init__.py", "亚细胞定位、细胞器分布和蛋白定位总览。", "可编辑 SVG + placement 列表；可控制细胞器数量、位置、大小和旋转，并附锚点与素材许可。", ("cell", "organelle", "localization", "schematic")),
    ]


def materialize_gallery_images(examples: list[Example]) -> list[Example]:
    """Copy example images that live outside the gallery into gallery assets."""
    if IMAGE_ASSETS.exists():
        shutil.rmtree(IMAGE_ASSETS)
    IMAGE_ASSETS.mkdir(parents=True, exist_ok=True)

    out_root = OUT.resolve()
    materialized: list[Example] = []
    for ex in examples:
        image = ex.image.resolve()
        if image.is_relative_to(out_root):
            materialized.append(ex)
            continue
        dest = IMAGE_ASSETS / f"{ex.slug}{image.suffix.lower()}"
        shutil.copy2(image, dest)
        materialized.append(replace(ex, image=dest))
    return materialized


def _categories(examples: list[Example]) -> list[str]:
    categories: list[str] = []
    for ex in examples:
        if ex.category not in categories:
            categories.append(ex.category)
    return categories


def _sub_categories(examples: list[Example]) -> dict[str, list[str]]:
    sub_categories: dict[str, list[str]] = {}
    for ex in examples:
        sub_categories.setdefault(ex.category, [])
        if ex.sub_category not in sub_categories[ex.category]:
            sub_categories[ex.category].append(ex.sub_category)
    return sub_categories


def _page_slug(category: str, sub_category: str | None = None) -> str:
    if sub_category is None:
        return _anchor("top", category).removeprefix("top-")
    return _anchor("sub", f"{category}-{sub_category}").removeprefix("sub-")


def _page_href(slug: str) -> str:
    return f"{slug}.html"


def _image_href(path: Path, *, from_pages: bool = True) -> str:
    rel = Path(os.path.relpath(path.resolve(), OUT.resolve())).as_posix()
    if from_pages:
        rel = f"../{rel}"
    return quote(rel, safe="/:#?&=%")


def _source_href(path: Path, *, from_pages: bool = True) -> str:
    return f"{REPO_BLOB_BASE}/{quote(_root_rel(path), safe='/')}"


def _super_switch(categories: list[str], current: str) -> str:
    out = ['<div class="super-switch"><span class="super-label">图类</span>']
    out.append(f'<a class="super-btn{" current" if current == "index" else ""}" href="index.html">Overview</a>')
    for cat in categories:
        slug = _page_slug(cat)
        cls = " current" if cat == current else ""
        out.append(f'<a class="super-btn{cls}" href="{_page_href(slug)}">{html.escape(cat)}</a>')
    out.append("</div>")
    return "".join(out)


def _top_nav(categories: list[str], current: str) -> str:
    items = [(_page_slug(cat), cat, _page_href(_page_slug(cat))) for cat in categories]
    return report.build_nav(items, current=_page_slug(current), label=None)


def _sidebar(
    examples: list[Example],
    current: str,
    current_sub: str | None = None,
    *,
    overview: bool = False,
) -> str:
    items: list[tuple] = [("index", "Overview", "index.html")]
    sub_categories = _sub_categories(examples)
    categories = _categories(examples)
    if overview:
        for cat in categories:
            items.append(("section", cat))
            items.append((_page_slug(cat), f"{cat} overview", _page_href(_page_slug(cat))))
            for sub in sub_categories[cat]:
                items.append((_page_slug(cat, sub), sub, _page_href(_page_slug(cat, sub))))
        return report.build_sidebar(items, current="index")

    items.append(("section", current))
    items.append((_page_slug(current), f"{current} overview", _page_href(_page_slug(current))))
    for sub in sub_categories[current]:
        items.append((_page_slug(current, sub), sub, _page_href(_page_slug(current, sub))))
        if current_sub == sub:
            for ex in examples:
                if ex.category == current and ex.sub_category == sub:
                    items.append((f"{ex.slug}", f"  {ex.title}", f"#{ex.slug}"))
    current_slug = _page_slug(current, current_sub) if current_sub else _page_slug(current)
    return report.build_sidebar(items, current=current_slug)


def _gallery_extra_css() -> str:
    return """
<style>
.super-switch {
  display: flex; gap: 8px; justify-content: center; align-items: center;
  margin: 0 0 14px 0; padding: 8px 10px;
  background: var(--bg-elev); border: 1px solid var(--border); border-radius: 12px;
  flex-wrap: wrap;
}
.super-switch .super-label {
  color: var(--muted); font-size: .78em; font-weight: 600;
  letter-spacing: .08em; margin-right: 6px; text-transform: uppercase;
}
.super-switch a.super-btn {
  padding: 6px 18px; border-radius: 999px; color: var(--fg); font-weight: 600;
  font-size: .95em; border: 1px solid var(--border); background: var(--bg-soft);
}
.super-switch a.super-btn:hover { background: var(--bg); color: var(--link-hover); text-decoration: none; }
.super-switch a.super-btn.current { background: var(--accent); color: #0d1117; border-color: var(--accent-strong); }
.gallery-grid {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 260px), 1fr));
  gap: 18px; margin: 1.2em 0 2em;
}
.gallery-search {
  display: grid; grid-template-columns: auto minmax(0, 1fr) auto;
  gap: 10px; align-items: center; margin: 1em 0 1.2em; padding: 10px 12px;
  border: 1px solid var(--border); border-radius: 8px; background: var(--bg-elev);
}
.gallery-search label {
  color: var(--muted); font-size: .86em; font-weight: 700; letter-spacing: .04em;
  text-transform: uppercase;
}
.gallery-search input {
  min-width: 0; width: 100%; box-sizing: border-box; padding: 8px 10px;
  border: 1px solid var(--border); border-radius: 6px;
  color: var(--fg); background: var(--bg-soft);
}
.gallery-search input:focus {
  outline: 2px solid var(--accent); outline-offset: 1px; border-color: var(--accent-strong);
}
.gallery-search-count { color: var(--muted); font-size: .86em; white-space: nowrap; }
.gallery-empty {
  margin: 1em 0 2em; padding: 18px; text-align: center;
  color: var(--muted); border: 1px dashed var(--border); border-radius: 8px;
  background: var(--bg-soft);
}
.gallery-card[hidden], .gallery-empty[hidden] { display: none !important; }
.gallery-card {
  border: 1px solid var(--border); border-radius: 8px; background: var(--bg-soft);
  overflow: hidden;
}
.gallery-card figure { margin: 0; }
.gallery-card img {
  height: 230px; width: 100%; object-fit: contain; margin: 0;
  border: 0; border-bottom: 1px solid var(--border); border-radius: 0;
}
.gallery-card .body { padding: 12px 14px 14px; }
.gallery-card h3 { margin: 0 0 8px; }
.gallery-card p { margin: .45em 0; }
.gallery-meta {
  color: var(--muted); font-size: .86em; word-break: break-all;
}
.gallery-tags { display: flex; flex-wrap: wrap; gap: 5px; margin-top: 8px; }
.gallery-tags span {
  border: 1px solid var(--border); border-radius: 4px; padding: 1px 6px;
  color: var(--accent-strong); background: var(--bg-elev); font-size: .82em;
}
.gallery-links { display: flex; gap: 14px; margin-top: 10px; font-weight: 600; }
@media (max-width: 640px) {
  .gallery-search { grid-template-columns: 1fr; }
  .gallery-search-count { white-space: normal; }
}
</style>
"""


def _header(examples: list[Example], current: str, current_sub: str | None = None, *, overview: bool = False) -> str:
    categories = _categories(examples)
    if overview:
        nav = report.build_nav(
            [(_page_slug(cat), cat, _page_href(_page_slug(cat))) for cat in categories],
            current=None,
            label=None,
        )
        sidebar = _sidebar(examples, current, overview=True)
        super_nav = _super_switch(categories, "index")
    else:
        nav = _top_nav(categories, current)
        sidebar = _sidebar(examples, current, current_sub)
        super_nav = _super_switch(categories, current)
    return _gallery_extra_css() + "\n" + sidebar + "\n" + super_nav + "\n" + nav + '\n<main class="report-main">'


def _example_card(ex: Example) -> str:
    img_url = _image_href(ex.image)
    source_url = _source_href(ex.source)
    tags = " ".join(f"<span>{html.escape(tag)}</span>" for tag in ex.tags)
    search_text = " ".join(
        [
            ex.slug,
            ex.title,
            ex.category,
            ex.sub_category,
            ex.use_case,
            ex.note,
            _root_rel(ex.source),
            *ex.tags,
        ]
    )
    return f"""
<article class="gallery-card" id="{html.escape(ex.slug)}" data-search="{html.escape(search_text, quote=True)}">
<figure>
<a href="{img_url}"><img src="{img_url}" alt="{html.escape(ex.title)}"></a>
</figure>
<div class="body">
<h3>{html.escape(ex.title)}</h3>
<p>{html.escape(ex.use_case)}</p>
<p><em>{html.escape(ex.note)}</em></p>
<div class="gallery-meta">{html.escape(ex.category)} / {html.escape(ex.sub_category)} · {html.escape(_root_rel(ex.source))}</div>
<div class="gallery-tags">{tags}</div>
<div class="gallery-links"><a href="{img_url}">原图</a><a href="{source_url}">源码/示例</a></div>
</div>
</article>
"""


def _gallery_search(total: int) -> str:
    return f"""
<div class="gallery-search" role="search" data-gallery-search>
<label>Search</label>
<input type="search" placeholder="Title, section, text" autocomplete="off" aria-label="Search examples">
<div class="gallery-search-count" data-gallery-count>{total} examples</div>
</div>
"""


def _gallery_empty_state() -> str:
    return '<p class="gallery-empty" data-gallery-empty hidden>No results</p>'


def _gallery_search_script() -> str:
    return """
<script>
document.addEventListener('DOMContentLoaded', function () {
  document.querySelectorAll('[data-gallery-search]').forEach(function (box) {
    var input = box.querySelector('input[type="search"]');
    var main = box.closest('main') || document;
    var cards = Array.prototype.slice.call(main.querySelectorAll('.gallery-card'));
    var empty = main.querySelector('[data-gallery-empty]');
    var count = box.querySelector('[data-gallery-count]');
    function normalize(value) {
      return (value || '').toLowerCase().trim();
    }
    function applyFilter() {
      var query = normalize(input.value);
      var shown = 0;
      cards.forEach(function (card) {
        var text = normalize(card.getAttribute('data-search') || card.textContent);
        var matched = !query || text.indexOf(query) !== -1;
        card.hidden = !matched;
        if (matched) shown += 1;
      });
      if (empty) empty.hidden = shown !== 0;
      if (count) count.textContent = query ? shown + ' / ' + cards.length + ' examples' : cards.length + ' examples';
    }
    input.addEventListener('input', applyFilter);
    input.addEventListener('search', applyFilter);
    applyFilter();
  });
});
</script>
"""


def _render_report_page(md: str, out_html: Path, title: str, header_html: str) -> None:
    report.render(
        md,
        out_html=out_html,
        css=report.css_path("dark"),
        toc=False,
        title=title,
        favicon_emoji="🎨",
        header_html=header_html,
        footer_html="</main>" + _gallery_search_script(),
        extra_args=["--resource-path", str(OUT)],
    )


def render_report_site(examples: list[Example]) -> None:
    if PAGES.exists():
        shutil.rmtree(PAGES)
    PAGES.mkdir(parents=True, exist_ok=True)
    categories = _categories(examples)
    sub_categories = _sub_categories(examples)

    overview_lines = [
        "# drawing-yh 图型示例库",
        "",
        f"{len(examples)} 种示例图，按标准三层目录组织：图类 → 子类 → 具体图型。",
        "",
        "## 图类",
        "",
    ]
    for cat in categories:
        n = sum(1 for ex in examples if ex.category == cat)
        overview_lines.append(f"- [{cat}]({_page_href(_page_slug(cat))}) — {n} examples")
    overview_lines += ["", "## 全部示例", "", _gallery_search(len(examples)), '<div class="gallery-grid">']
    overview_lines += [_example_card(ex) for ex in examples]
    overview_lines += ["</div>", _gallery_empty_state()]
    _render_report_page(
        "\n".join(overview_lines) + "\n",
        PAGES / "index.html",
        "drawing-yh 图型示例库",
        _header(examples, categories[0], overview=True),
    )

    for cat in categories:
        cat_examples = [ex for ex in examples if ex.category == cat]
        lines = [f"# {cat}", "", f"{len(cat_examples)} 种示例图。", "", "## 子类", ""]
        for sub in sub_categories[cat]:
            sub_examples = [ex for ex in cat_examples if ex.sub_category == sub]
            lines.append(f"- [{sub}]({_page_href(_page_slug(cat, sub))}) — {len(sub_examples)} examples")
        lines += ["", "## 全部示例", "", _gallery_search(len(cat_examples)), '<div class="gallery-grid">']
        lines += [_example_card(ex) for ex in cat_examples]
        lines += ["</div>", _gallery_empty_state()]
        _render_report_page(
            "\n".join(lines) + "\n",
            PAGES / _page_href(_page_slug(cat)),
            f"{cat} gallery",
            _header(examples, cat),
        )

        for sub in sub_categories[cat]:
            sub_examples = [ex for ex in cat_examples if ex.sub_category == sub]
            lines = [
                f"# {sub}",
                "",
                f"{cat} / {sub} — {len(sub_examples)} 种示例图。",
                "",
                _gallery_search(len(sub_examples)),
                '<div class="gallery-grid">',
            ]
            lines += [_example_card(ex) for ex in sub_examples]
            lines += ["</div>", _gallery_empty_state()]
            _render_report_page(
                "\n".join(lines) + "\n",
                PAGES / _page_href(_page_slug(cat, sub)),
                f"{cat} / {sub}",
                _header(examples, cat, sub),
            )

    legacy_routes = {
        "": "pages/index.html",
        "top-quantitative": "pages/quantitative.html",
        "Quantitative": "pages/quantitative.html",
        "Bar": "pages/quantitative-bar.html",
        "Distribution": "pages/quantitative-distribution.html",
        "Scatter": "pages/quantitative-scatter.html",
        "Omics": "pages/omics.html",
        "Heatmap": "pages/omics-heatmap.html",
        "Composition": "pages/composition.html",
        "Curve": "pages/curves.html",
        "Curves": "pages/curves.html",
        "Network": "pages/network.html",
        "Utility": "pages/utility.html",
    }
    for cat in categories:
        legacy_routes[_anchor("top", cat)] = f"pages/{_page_href(_page_slug(cat))}"
        for sub in sub_categories[cat]:
            legacy_routes[_anchor("sub", f"{cat}-{sub}")] = f"pages/{_page_href(_page_slug(cat, sub))}"
    for ex in examples:
        legacy_routes[ex.slug] = f"pages/{_page_href(_page_slug(ex.category, ex.sub_category))}#{ex.slug}"

    routes_json = json.dumps(legacy_routes, ensure_ascii=False)
    (OUT / "index.html").write_text(
        f"""<!doctype html>
<meta charset="utf-8">
<title>redirect</title>
<script>
const routes = {routes_json};
const key = decodeURIComponent(location.hash.replace(/^#/, ""));
location.replace(routes[key] || routes[""]);
</script>
<p>Redirect to <a href="pages/index.html">report</a>.</p>
""",
        encoding="utf-8",
    )


def main() -> None:
    generated = generate_all()
    examples = build_examples(generated)
    missing = [ex for ex in examples if not ex.image.exists()]
    if missing:
        miss = "\n".join(f"- {ex.slug}: {ex.image}" for ex in missing)
        raise FileNotFoundError(f"Missing gallery images:\n{miss}")
    examples = materialize_gallery_images(examples)
    render_report_site(examples)
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
