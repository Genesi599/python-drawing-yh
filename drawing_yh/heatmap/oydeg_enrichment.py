"""OY-DEG heatmap with dynamic gene labels and enrichment blocks.

This module is a pure matplotlib drawing layer. It expects upstream code to
prepare DEG rows, selected label genes, gene direction status, and enrichment
terms. Filtering, enrichment, and term deduplication belong in the analysis
pipeline or single_cell-yh.
"""
from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
import textwrap

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, to_hex, to_rgb
import numpy as np
import pandas as pd


WARM = "#C0392B"
COOL = "#2166AC"
DEFAULT_DIVERGING = ("#2166AC", "#4393C3", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#D6604D", "#B2182B")


def darken(hex_color: str, factor: float = 0.55) -> str:
    """Return a darker version of a hex color."""
    return to_hex([c * factor for c in to_rgb(hex_color)])


def cap_gene(gene: str) -> str:
    """Display gene names with a leading capital unless already all-uppercase."""
    gene = str(gene)
    return gene if gene.isupper() else gene.capitalize()


def _resolve_values(values, columns=None, row_labels=None):
    if isinstance(values, pd.DataFrame):
        df = values.copy()
        if columns is not None:
            df = df.reindex(columns=list(columns))
        if row_labels is not None:
            df = df.reindex(index=list(row_labels))
        cols = list(df.columns)
        rows = [str(x) for x in df.index]
        mat = df.values.astype(float)
        return mat, cols, rows

    mat = np.asarray(values, dtype=float)
    if mat.ndim != 2:
        raise ValueError(f"values must be 2-D, got shape {mat.shape}")
    n_rows, n_cols = mat.shape
    cols = list(columns) if columns is not None else [str(i) for i in range(n_cols)]
    rows = [str(x) for x in row_labels] if row_labels is not None else [str(i) for i in range(n_rows)]
    if len(cols) != n_cols:
        raise ValueError(f"columns length {len(cols)} != values columns {n_cols}")
    if len(rows) != n_rows:
        raise ValueError(f"row_labels length {len(rows)} != values rows {n_rows}")
    return mat, cols, rows


def _resolve_cmap(cmap, missing_color: str):
    if cmap is None:
        cmap_obj = LinearSegmentedColormap.from_list("oydeg_rdbu", list(DEFAULT_DIVERGING))
    elif isinstance(cmap, str):
        cmap_obj = plt.get_cmap(cmap).copy()
    else:
        cmap_obj = cmap.copy() if hasattr(cmap, "copy") else cmap
    if hasattr(cmap_obj, "set_bad"):
        cmap_obj.set_bad(missing_color)
    return cmap_obj


def _display(label, display_labels=None) -> str:
    if display_labels is None:
        return str(label)
    if isinstance(display_labels, Mapping):
        return str(display_labels.get(label, label))
    if callable(display_labels):
        return str(display_labels(label))
    return str(label)


def _as_heat_blocks(heat_blocks, n_rows: int, default_color: str):
    if heat_blocks is None:
        return [dict(label="", color=default_color, s=0, e=n_rows)]
    out = []
    for block in heat_blocks:
        b = dict(block)
        b.setdefault("label", "")
        b.setdefault("color", default_color)
        if "s" not in b or "e" not in b:
            raise ValueError("Each heat block must contain 's' and 'e' row indices")
        out.append(b)
    return out


def _add_block_gaps(mat: np.ndarray, heat_blocks: list[dict], gap_frac: float):
    n, ncol = mat.shape
    gap = max(1, round(n * gap_frac)) if len(heat_blocks) > 1 else 0
    disp_idx = np.arange(n, dtype=int)
    cum = 0
    for bi, block in enumerate(heat_blocks):
        block["s_disp"], block["e_disp"] = block["s"] + cum, block["e"] + cum
        disp_idx[block["s"]:block["e"]] = np.arange(block["s"], block["e"], dtype=int) + cum
        if bi < len(heat_blocks) - 1:
            cum += gap
    n_disp = n + gap * max(len(heat_blocks) - 1, 0)
    mat_disp = np.full((n_disp, ncol), np.nan)
    mat_disp[disp_idx] = mat
    return mat_disp, disp_idx, n_disp


def _resolve_rep_col(rep_col, mat: np.ndarray):
    n, _ = mat.shape
    if rep_col is None:
        rep = []
        for i in range(n):
            row = np.nan_to_num(mat[i], nan=0.0)
            rep.append(int(np.argmax(np.abs(row))))
        return rep
    if isinstance(rep_col, Mapping):
        raise TypeError("Mapping rep_col is not supported; pass a sequence aligned to rows")
    rep = list(rep_col)
    if len(rep) != n:
        raise ValueError(f"rep_col length {len(rep)} != values rows {n}")
    return [int(x) for x in rep]


def _status_for(gene_status, gene: str, column) -> str | None:
    if gene_status is None:
        return None
    if isinstance(gene_status, Mapping):
        item = gene_status.get(gene)
        if isinstance(item, Mapping):
            return item.get(column)
        return gene_status.get((gene, column))
    return None


def _normalise_right_blocks(right_blocks, columns, column_colors, display_labels, default_color: str):
    if right_blocks is None:
        return []
    if isinstance(right_blocks, Mapping):
        blocks = []
        for col in columns:
            groups = right_blocks.get(col, [])
            blocks.append(dict(
                label=_display(col, display_labels),
                color=column_colors.get(col, default_color) if column_colors else default_color,
                groups=groups,
            ))
        return blocks
    out = []
    for block in right_blocks:
        b = dict(block)
        b.setdefault("label", "")
        b.setdefault("color", default_color)
        b.setdefault("groups", [])
        out.append(b)
    return out


def _term_dict(term) -> dict:
    if isinstance(term, Mapping):
        return dict(term)
    if len(term) == 2:
        label, value = term
        return {"term": label, "nlp": value, "db": "", "genes": []}
    raise ValueError("Term entries must be dicts or (term, value) tuples")


def _build_enrichment_layout(
    right_blocks,
    *,
    max_genes_per_term: int,
    gene_formatter: Callable[[str], str],
    gene_wrap: int,
):
    head_h, sub_h, term_h, gene_h, term_gap, gap = 1.0, 0.7, 0.85, 0.65, 0.25, 0.55
    layout, slot = [], 0.0
    for block in right_blocks:
        header_y = slot
        slot += head_h
        groups_layout, any_term = [], False
        for group in block.get("groups", []):
            dlabel, dcolor, terms_raw = group
            terms = [_term_dict(t) for t in terms_raw]
            any_term = any_term or bool(terms)
            sub_y = slot
            slot += sub_h
            items = []
            for term in terms:
                term_y = slot
                slot += term_h
                genes = [str(x) for x in term.get("genes", []) if str(x)]
                shown = genes[:max_genes_per_term]
                extra = len(genes) - len(shown)
                gene_text = ", ".join(gene_formatter(x) for x in shown)
                if extra:
                    gene_text += f" (+{extra})"
                glines = []
                for line in (textwrap.wrap(gene_text, gene_wrap) or [""]):
                    glines.append((slot, line))
                    slot += gene_h
                slot += term_gap
                items.append(dict(term_y=term_y, t=term, glines=glines))
            groups_layout.append(dict(dlabel=dlabel, dcolor=dcolor, sub_y=sub_y, items=items))
        note_y = None
        if not any_term:
            note_y = slot
            slot += 0.85
        layout.append(dict(
            label=block.get("label", ""),
            color=block.get("color", "#777777"),
            header_y=header_y,
            glayout=groups_layout,
            note_y=note_y,
        ))
        slot += gap
    total_slots = max(slot, 1.0)
    max_nlp = max([float(item["t"].get("nlp", 0.0))
                   for block in layout for group in block["glayout"] for item in group["items"]] + [1.5])
    return layout, total_slots, max_nlp


def plot_oydeg_heatmap_enrichment(
    values,
    *,
    columns: Sequence[str] | None = None,
    row_labels: Sequence[str] | None = None,
    heat_blocks: Sequence[Mapping] | None = None,
    label_genes: Sequence[str] | None = None,
    rep_col: Sequence[int] | None = None,
    gene_status: Mapping | None = None,
    column_colors: Mapping | None = None,
    display_labels: Mapping | Callable[[str], str] | None = None,
    right_blocks: Sequence[Mapping] | Mapping | None = None,
    fc_clip: float = 2.0,
    cmap=None,
    missing_color: str = "#ffffff",
    heatmap_title: str | None = None,
    right_title: str = "KEGG + GO enrichment",
    right_xlabel: str = "-log10(padj)",
    colorbar_label: str = "log2FC",
    no_terms_label: str = "no significant pathway",
    max_genes_per_term: int = 8,
    gene_wrap: int = 70,
    gene_formatter: Callable[[str], str] = cap_gene,
    fig_width: float | None = None,
    figsize: tuple[float, float] | None = None,
    slot_height: float = 0.16,
    block_gap_frac: float = 0.015,
    dot_dx: float = 0.05,
    font_sizes: Mapping[str, float] | None = None,
):
    """Draw an OY-DEG heatmap with labeled genes and enrichment blocks.

    Parameters are intentionally presentation-level. Upstream code should pass
    ordered rows, block boundaries, selected labels, and enrichment term lists.

    ``right_blocks`` is a sequence of dictionaries:
    ``{"label": str, "color": hex, "groups": [(group_label, group_color, terms)]}``.
    Each term is a dict with ``term``, ``nlp``, optional ``db`` and ``genes``.
    """
    mat, cols, rows = _resolve_values(values, columns, row_labels)
    n, ncol = mat.shape
    if n == 0 or ncol == 0:
        raise ValueError("values must contain at least one row and one column")

    column_colors = dict(column_colors or {})
    heat_blocks = _as_heat_blocks(heat_blocks, n, "#777777")
    rep_col = _resolve_rep_col(rep_col, mat)
    if label_genes is None:
        label_genes = rows[:min(len(rows), 40)]
    label_genes = [str(g) for g in label_genes]
    row_of = {g: i for i, g in enumerate(rows)}
    label_genes = [g for g in label_genes if g in row_of]
    if not label_genes:
        label_genes = rows[:min(len(rows), 40)]

    cmap_obj = _resolve_cmap(cmap, missing_color)
    norm = TwoSlopeNorm(vmin=-fc_clip, vcenter=0.0, vmax=fc_clip)
    mat_disp, disp_idx, n_disp = _add_block_gaps(np.clip(mat, -fc_clip, fc_clip), heat_blocks, block_gap_frac)

    right_blocks_norm = _normalise_right_blocks(right_blocks, cols, column_colors, display_labels, "#777777")
    layout, total_slots, max_nlp = _build_enrichment_layout(
        right_blocks_norm,
        max_genes_per_term=max_genes_per_term,
        gene_formatter=gene_formatter,
        gene_wrap=gene_wrap,
    )

    fs = {
        "title": 8.0, "xtick": 7.0, "block": 6.5, "cbar": 7.0, "ctick": 6.0,
        "label": 6.0, "head": 8.0, "sub": 7.0, "term": 6.5, "gene": 5.5,
        "note": 6.5, "axis": 7.0, "btick": 6.0,
    }
    if font_sizes:
        fs.update(font_sizes)

    if figsize is None:
        fig_h = max(4.5, total_slots * slot_height, n_disp * slot_height * 0.95)
        fig_w = fig_width if fig_width is not None else 11.0 + ncol * 0.3
        figsize = (fig_w, fig_h)
    fig = plt.figure(figsize=figsize)

    hm_w = max(0.10, ncol * 0.012)
    ax_heat = fig.add_axes([0.05, 0.055, hm_w, 0.90])
    ax_labels = fig.add_axes([0.05 + hm_w + 0.015, 0.055, 0.17, 0.90])
    bar_x0 = 0.05 + hm_w + 0.015 + 0.19
    ax_enrich = fig.add_axes([bar_x0, 0.055, 0.99 - bar_x0 - 0.01, 0.90])
    cax = fig.add_axes([0.05, 0.022, 0.10, 0.010])

    ax_heat.imshow(mat_disp, aspect="auto", cmap=cmap_obj, norm=norm,
                   extent=[0, ncol, n_disp, 0], interpolation="nearest")
    ax_heat.set_xticks([i + 0.5 for i in range(ncol)])
    ax_heat.set_xticklabels([_display(c, display_labels) for c in cols], rotation=45,
                            ha="right", fontsize=fs["xtick"])
    for tick, col in zip(ax_heat.get_xticklabels(), cols):
        tick.set_color(darken(column_colors.get(col, "#000000"), 0.62))
    ax_heat.set_yticks([])
    for i in range(1, ncol):
        ax_heat.axvline(i, color="#ffffff", lw=0.4, zorder=3)
    for block in heat_blocks:
        bc = darken(block.get("color", "#777777"), 0.62)
        s, e = block["s_disp"], block["e_disp"]
        ax_heat.plot([0, 0], [s, e], color=bc, lw=0.9, zorder=4, clip_on=False)
        ax_heat.plot([ncol, ncol], [s, e], color=bc, lw=0.9, zorder=4, clip_on=False)
        if block.get("label"):
            ax_heat.text(-0.08 * ncol, (s + e) / 2, str(block["label"]),
                         va="center", ha="right", fontsize=fs["block"], fontweight="bold",
                         color=bc, clip_on=False)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)
    ax_heat.set_xlim(0, ncol)
    ax_heat.set_ylim(n_disp, 0)
    if heatmap_title:
        ax_heat.set_title(heatmap_title, loc="left", fontsize=fs["title"])

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", extend="both")
    cb.set_label(colorbar_label, fontsize=fs["cbar"])
    cb.ax.tick_params(labelsize=fs["ctick"])

    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(n_disp, 0)
    ax_labels.axis("off")
    nl = max(len(label_genes), 1)
    for i, gene in enumerate(label_genes):
        ly = (i + 0.5) * n_disp / nl
        ri = row_of[gene]
        fc0 = mat[ri, rep_col[ri]]
        color = cmap_obj(norm(float(np.clip(fc0, -fc_clip, fc_clip)))) if not np.isnan(fc0) else "#aaaaaa"
        ax_labels.plot([0.02, 0.09, 0.17, 0.22],
                       [disp_idx[ri] + 0.5, disp_idx[ri] + 0.5, ly, ly],
                       color=color, lw=0.45, alpha=0.85, zorder=1, solid_capstyle="round")
        for j, col in enumerate(cols):
            x = 0.22 + j * dot_dx
            ax_labels.scatter(x, ly, s=12, marker="s", color="#efefef",
                              edgecolors="none", zorder=2)
            status = _status_for(gene_status, gene, col)
            if status:
                ax_labels.scatter(x, ly, s=8, marker=("^" if status == "up" else "v"),
                                  color=column_colors.get(col, "#777777"),
                                  edgecolors="none", zorder=3)
        name_x = 0.22 + ncol * dot_dx + 0.02
        ax_labels.text(name_x, ly, gene_formatter(gene), va="center", ha="left",
                       fontsize=fs["label"], color="#222222")

    ax_enrich.set_xlim(0, max_nlp * 1.05)
    ax_enrich.set_ylim(total_slots, -1.0)
    xt = max_nlp * 0.03
    for block in layout:
        ax_enrich.axhline(block["header_y"] - 0.6, color="#cccccc", lw=0.5, zorder=1)
        ax_enrich.text(0, block["header_y"], str(block["label"]),
                       va="center", ha="left", fontsize=fs["head"], fontweight="bold",
                       color=darken(block["color"], 0.62), zorder=5)
        if block["note_y"] is not None:
            ax_enrich.text(xt, block["note_y"], no_terms_label, va="center", ha="left",
                           fontsize=fs["note"], color="#999999", fontstyle="italic")
        for group in block["glayout"]:
            ax_enrich.text(xt, group["sub_y"], group["dlabel"],
                           va="center", ha="left", fontsize=fs["sub"], fontweight="bold",
                           color=darken(group["dcolor"], 0.7))
            for item in group["items"]:
                term = item["t"]
                y = item["term_y"]
                nlp = float(term.get("nlp", 0.0))
                ax_enrich.barh(y, nlp, height=0.60, color=group["dcolor"],
                               alpha=0.28, edgecolor="none")
                db = str(term.get("db", "")).strip()
                label = f"{term.get('term', '')} ({db})" if db else str(term.get("term", ""))
                ax_enrich.text(xt, y, label, va="center", ha="left",
                               fontsize=fs["term"], color="black", fontweight="bold")
                for gy, line in item["glines"]:
                    ax_enrich.text(xt, gy, line, va="center", ha="left",
                                   fontsize=fs["gene"], color=darken(group["dcolor"]),
                                   fontstyle="italic")
    ax_enrich.set_xlabel(right_xlabel, fontsize=fs["axis"])
    ax_enrich.spines[["top", "right", "left"]].set_visible(False)
    ax_enrich.tick_params(left=False, labelleft=False, labelsize=fs["btick"])
    if right_title:
        ax_enrich.set_title(right_title, loc="left", fontsize=fs["title"])

    axes = {"heatmap": ax_heat, "labels": ax_labels, "enrichment": ax_enrich, "colorbar": cax}
    return fig, axes


__all__ = [
    "plot_oydeg_heatmap_enrichment",
    "darken",
    "cap_gene",
    "WARM",
    "COOL",
]
