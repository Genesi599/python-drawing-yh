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


def _wrap_label(label: str, width: int = 14) -> str:
    label = str(label).replace("_", " ")
    return "\n".join(textwrap.wrap(label, width)) or label


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
    head_h, sub_h, term_h, gene_h, term_gap, gap = 1.05, 0.78, 0.88, 0.70, 0.28, 0.60
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


def _build_enrichment_columns(right_blocks, ncols, **kw):
    """把 enrichment blocks 顺序拆成 ncols 列(按累计高度均衡,保留 block 原序),
    每列独立从 y=0 排版。返回 [(layout, total_slots), ...] + 全局 max_nlp。
    ncols=1 时与单列 _build_enrichment_layout 完全等价。"""
    if ncols <= 1 or len(right_blocks) <= 1:
        layout, total, max_nlp = _build_enrichment_layout(right_blocks, **kw)
        return [(layout, total)], max_nlp
    heights = [_build_enrichment_layout([b], **kw)[1] for b in right_blocks]
    target = sum(heights) / ncols
    cols, cur, cum = [[]], 0.0, 0
    for b, h in zip(right_blocks, heights):
        if len(cols) < ncols and cur >= target - 1e-9 and cols[-1]:
            cols.append([]); cur = 0.0
        cols[-1].append(b); cur += h
    out, max_nlp = [], 1.5
    for cb in cols:
        layout, total, mnlp = _build_enrichment_layout(cb, **kw)
        out.append((layout, total)); max_nlp = max(max_nlp, mnlp)
    return out, max_nlp


def _lay_terms_from(terms_raw, start_y, *, max_genes_per_term, gene_formatter,
                    gene_wrap, term_wrap):
    """从 start_y 往下排一组 term(可换行的 term 标签 + 基因行),返回 (items, end_y)。"""
    term_h, gene_h, term_gap, lbl_line_h = 0.88, 0.70, 0.28, 0.66
    slot = start_y
    items = []
    for t in (_term_dict(x) for x in terms_raw):
        db = str(t.get("db", "")).strip()
        lbl = f"{t.get('term', '')} ({db})" if db else str(t.get("term", ""))
        lbl_lines = textwrap.wrap(lbl, term_wrap) or [lbl]
        term_y = slot
        slot += term_h + (len(lbl_lines) - 1) * lbl_line_h
        genes = [str(x) for x in t.get("genes", []) if str(x)]
        shown = genes[:max_genes_per_term]
        extra = len(genes) - len(shown)
        gtext = ", ".join(gene_formatter(x) for x in shown) + (f" (+{extra})" if extra else "")
        glines = []
        for line in (textwrap.wrap(gtext, gene_wrap) or [""]):
            if line.strip():
                glines.append((slot, line)); slot += gene_h
        slot += term_gap
        items.append(dict(term_y=term_y, lbl_lines=lbl_lines, glines=glines, t=t))
    return items, slot


def _build_updown_layout(right_blocks, *, max_genes_per_term, gene_formatter,
                         gene_wrap, term_wrap):
    """每个 block 的 up / down 两组并排放(左 up / 右 down),块高 = head + max(up,down)。"""
    head_h, sub_h, gap = 1.05, 0.78, 0.55
    layout, slot = [], 0.0
    for block in right_blocks:
        header_y = slot
        slot += head_h
        sides, any_term = [], False
        for group in block.get("groups", [])[:2]:
            dlabel, dcolor, terms_raw = group
            any_term = any_term or bool(terms_raw)
            items, end_y = _lay_terms_from(
                terms_raw, slot + sub_h, max_genes_per_term=max_genes_per_term,
                gene_formatter=gene_formatter, gene_wrap=gene_wrap, term_wrap=term_wrap)
            sides.append(dict(dlabel=dlabel, dcolor=dcolor, sub_y=slot, items=items, end_y=end_y))
        block_end = max([s["end_y"] for s in sides], default=slot)
        note_y = None
        if not any_term:
            note_y = slot; block_end = slot + 0.85
        layout.append(dict(label=block.get("label", ""), color=block.get("color", "#777777"),
                           header_y=header_y, sides=sides, note_y=note_y))
        slot = block_end + gap
    total = max(slot, 1.0)
    max_nlp = max([float(it["t"].get("nlp", 0.0))
                   for b in layout for s in b["sides"] for it in s["items"]] + [1.5])
    return layout, total, max_nlp


def _render_updown(ax_up, ax_down, layout, total_slots, max_nlp, *, fs, no_terms_label):
    """up 组渲染到 ax_up,down 组到 ax_down,header 跨在 ax_up 左缘,两轴 y 同步。"""
    for ax in (ax_up, ax_down):
        ax.set_xlim(0, max_nlp * 1.05)
        ax.set_ylim(total_slots, -1.0)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(left=False, labelleft=False, labelsize=fs["btick"])
    xt = max_nlp * 0.03
    for block in layout:
        ax_up.axhline(block["header_y"] - 0.6, color="#cccccc", lw=0.5, zorder=1)
        ax_down.axhline(block["header_y"] - 0.6, color="#cccccc", lw=0.5, zorder=1)
        ax_up.text(0, block["header_y"], str(block["label"]), va="center", ha="left",
                   fontsize=fs["head"], fontweight="bold",
                   color=darken(block["color"], 0.62), zorder=5, clip_on=False)
        if block["note_y"] is not None:
            ax_up.text(xt, block["note_y"], no_terms_label, va="center", ha="left",
                       fontsize=fs["note"], color="#999999", fontstyle="italic")
        for k, side in enumerate(block["sides"]):
            ax = ax_up if k == 0 else ax_down
            ax.text(xt, side["sub_y"], side["dlabel"], va="center", ha="left",
                    fontsize=fs["sub"], fontweight="bold",
                    color=darken(side["dcolor"], 0.7), clip_on=False)
            for it in side["items"]:
                nlp = float(it["t"].get("nlp", 0.0))
                ax.barh(it["term_y"], nlp, height=0.60, color=side["dcolor"],
                        alpha=0.28, edgecolor="none")
                for li, line in enumerate(it["lbl_lines"]):
                    ax.text(xt, it["term_y"] + li * 0.66, line, va="center", ha="left",
                            fontsize=fs["term"], color="black", fontweight="bold", clip_on=False)
                for gy, line in it["glines"]:
                    ax.text(xt, gy, line, va="center", ha="left", fontsize=fs["gene"],
                            color=darken(side["dcolor"]), fontstyle="italic", clip_on=False)


def _render_enrich_column(ax, layout, total_slots, max_nlp, *, fs,
                          no_terms_label):
    """在单个 enrich 轴上渲染一列 blocks(从 _build_enrichment_layout 出来的 layout)。"""
    ax.set_xlim(0, max_nlp * 1.05)
    ax.set_ylim(total_slots, -1.0)
    xt = max_nlp * 0.03
    for block in layout:
        ax.axhline(block["header_y"] - 0.6, color="#cccccc", lw=0.5, zorder=1)
        ax.text(0, block["header_y"], str(block["label"]),
                va="center", ha="left", fontsize=fs["head"], fontweight="bold",
                color=darken(block["color"], 0.62), zorder=5, clip_on=False)
        if block["note_y"] is not None:
            ax.text(xt, block["note_y"], no_terms_label, va="center", ha="left",
                    fontsize=fs["note"], color="#999999", fontstyle="italic")
        for group in block["glayout"]:
            ax.text(xt, group["sub_y"], group["dlabel"],
                    va="center", ha="left", fontsize=fs["sub"], fontweight="bold",
                    color=darken(group["dcolor"], 0.7), clip_on=False)
            for item in group["items"]:
                term = item["t"]
                y = item["term_y"]
                nlp = float(term.get("nlp", 0.0))
                ax.barh(y, nlp, height=0.60, color=group["dcolor"],
                        alpha=0.28, edgecolor="none")
                db = str(term.get("db", "")).strip()
                label = f"{term.get('term', '')} ({db})" if db else str(term.get("term", ""))
                ax.text(xt, y, label, va="center", ha="left",
                        fontsize=fs["term"], color="black", fontweight="bold", clip_on=False)
                for gy, line in item["glines"]:
                    ax.text(xt, gy, line, va="center", ha="left",
                            fontsize=fs["gene"], color=darken(group["dcolor"]),
                            fontstyle="italic", clip_on=False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(left=False, labelleft=False, labelsize=fs["btick"])


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
    gene_wrap: int = 85,
    gene_formatter: Callable[[str], str] = cap_gene,
    fig_width: float | None = None,
    figsize: tuple[float, float] | None = None,
    slot_height: float = 0.185,
    block_gap_frac: float = 0.015,
    dot_dx: float = 0.060,
    font_sizes: Mapping[str, float] | None = None,
    compact_height: bool = True,
    enrich_ncols: int = 1,
    enrich_updown_cols: bool = False,
    updown_gap: float = 0.012,
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
    enrich_ncols = max(1, int(enrich_ncols))
    if enrich_updown_cols:
        # 每个细胞类型块的 up / down term 左右并排(块高 = head + max(up,down))
        updown_layout, total_slots, max_nlp = _build_updown_layout(
            right_blocks_norm,
            max_genes_per_term=max_genes_per_term,
            gene_formatter=gene_formatter,
            gene_wrap=max(24, gene_wrap // 2),
            term_wrap=38,
        )
        enrich_cols = None
    else:
        # 多列时按列宽缩小 gene_wrap,避免左列文字串到右列
        col_gene_wrap = gene_wrap if enrich_ncols == 1 else max(28, gene_wrap // enrich_ncols)
        enrich_cols, max_nlp = _build_enrichment_columns(
            right_blocks_norm, enrich_ncols,
            max_genes_per_term=max_genes_per_term,
            gene_formatter=gene_formatter,
            gene_wrap=col_gene_wrap,
        )
        total_slots = max((t for _, t in enrich_cols), default=1.0)

    fs = {
        "title": 8.0, "xtick": 8.0, "block": 8.0, "cbar": 8.0, "ctick": 8.0,
        "label": 8.0, "head": 8.0, "sub": 8.0, "term": 8.0, "gene": 8.0,
        "note": 8.0, "axis": 8.0, "btick": 8.0,
    }
    if font_sizes:
        fs.update(font_sizes)

    if figsize is None:
        label_h = max(len(label_genes), 1) * slot_height * 0.72
        row_h = 0.0 if compact_height else n_disp * slot_height * 0.95
        fig_h = max(4.5, total_slots * slot_height, label_h, row_h)
        fig_w = fig_width if fig_width is not None else 13.0
        figsize = (fig_w, fig_h)
    fig = plt.figure(figsize=figsize)

    ax_heat = fig.add_axes([0.06, 0.055, 0.13, 0.90])
    ax_labels = fig.add_axes([0.19, 0.055, 0.18, 0.90])
    dot_x0 = 0.25
    name_x = dot_x0 + ncol * dot_dx + 0.02
    max_gene_chars = max((len(gene_formatter(g)) for g in label_genes), default=6)
    bar_x0 = 0.19 + name_x * 0.18 + max_gene_chars * (fs["label"] * 0.62 / 72) / figsize[0] + 0.015
    bar_x0 = min(max(bar_x0, 0.36), 0.62)
    enrich_region = max(0.985 - bar_x0, 0.13)
    if enrich_updown_cols:
        # up 左 / down 右,两轴贴近(小间隙 updown_gap)
        half = (enrich_region - updown_gap) / 2
        ax_up = fig.add_axes([bar_x0, 0.055, half, 0.90])
        ax_down = fig.add_axes([bar_x0 + half + updown_gap, 0.055, half, 0.90])
        ax_enrich_list = [ax_up, ax_down]
    else:
        col_region = enrich_region / enrich_ncols
        enrich_axw = min(0.13, col_region * 0.92)
        ax_enrich_list = [
            fig.add_axes([bar_x0 + i * col_region, 0.055, enrich_axw, 0.90])
            for i in range(enrich_ncols)
        ]
    ax_enrich = ax_enrich_list[0]
    cax = fig.add_axes([0.205, 0.022, 0.11, 0.010])

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
    for block_i, block in enumerate(heat_blocks):
        bc = darken(block.get("color", "#777777"), 0.62)
        s, e = block["s_disp"], block["e_disp"]
        ax_heat.plot([0, 0], [s, e], color=bc, lw=0.9, zorder=4, clip_on=False)
        ax_heat.plot([ncol, ncol], [s, e], color=bc, lw=0.9, zorder=4, clip_on=False)
        if block_i == 0:
            ax_heat.plot([0, ncol], [s, s], color=bc, lw=0.9, zorder=4, clip_on=False)
        if block_i == len(heat_blocks) - 1:
            ax_heat.plot([0, ncol], [e, e], color=bc, lw=0.9, zorder=4, clip_on=False)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    centers = [(block["s_disp"] + block["e_disp"]) / 2 for block in heat_blocks]
    labels = [_wrap_label(block.get("label", "")) for block in heat_blocks]
    line_row = max((fs["block"] * 1.5 / 72) / (0.90 * figsize[1]) * n_disp, 1e-6)
    adjusted = list(centers)
    for i in range(1, len(adjusted)):
        sep = line_row * ((labels[i - 1].count("\n") + 1) + (labels[i].count("\n") + 1)) / 2 + 0.4 * line_row
        if adjusted[i] < adjusted[i - 1] + sep:
            adjusted[i] = adjusted[i - 1] + sep
    if adjusted and adjusted[-1] > n_disp:
        adjusted = [y - (adjusted[-1] - n_disp) for y in adjusted]
    col_per_frac = ncol / 0.13
    stub_x, diag_x, label_x = -0.010 * col_per_frac, -0.030 * col_per_frac, -0.035 * col_per_frac
    for block, center, y, label in zip(heat_blocks, centers, adjusted, labels):
        if not label:
            continue
        bc = darken(block.get("color", "#777777"), 0.62)
        if abs(y - center) > 0.3 * line_row:
            ax_heat.plot([0, stub_x, diag_x], [center, center, y], color=bc, lw=0.6,
                         alpha=0.85, clip_on=False, zorder=4, solid_capstyle="round",
                         solid_joinstyle="round")
        ax_heat.text(label_x, y, label, va="center", ha="right", fontsize=fs["block"],
                     fontweight="bold", color=bc, clip_on=False, linespacing=0.9)
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
    header_marker_size, grid_marker_size, tri_marker_size, leader_lw = 26, 16, 34, 0.9
    for j, col in enumerate(cols):
        x = dot_x0 + j * dot_dx
        ax_labels.scatter(x, -0.006 * n_disp, s=header_marker_size, marker="s",
                          color=column_colors.get(col, "#777777"), edgecolors="#ffffff",
                          linewidths=0.3, clip_on=False, zorder=3)
    for i, gene in enumerate(label_genes):
        ly = (i + 0.5) * n_disp / nl
        ri = row_of[gene]
        fc0 = mat[ri, rep_col[ri]]
        color = cmap_obj(norm(float(np.clip(fc0, -fc_clip, fc_clip)))) if not np.isnan(fc0) else "#aaaaaa"
        ax_labels.plot([0.00, 0.06, 0.16, 0.22],
                       [disp_idx[ri] + 0.5, disp_idx[ri] + 0.5, ly, ly],
                       color=color, lw=leader_lw, alpha=0.85, zorder=1,
                       solid_capstyle="round", solid_joinstyle="round")
        for j, col in enumerate(cols):
            x = dot_x0 + j * dot_dx
            ax_labels.scatter(x, ly, s=grid_marker_size, marker="s", color="#efefef",
                              edgecolors="none", zorder=2, clip_on=False)
            status = _status_for(gene_status, gene, col)
            if status:
                ax_labels.scatter(x, ly, s=tri_marker_size, marker=("^" if status == "up" else "v"),
                                  color=column_colors.get(col, "#777777"),
                                  edgecolors="none", zorder=3, clip_on=False)
        ax_labels.text(name_x, ly, gene_formatter(gene), va="center", ha="left",
                       fontsize=fs["label"], color="#222222", clip_on=False)

    if enrich_updown_cols:
        _render_updown(ax_enrich_list[0], ax_enrich_list[1], updown_layout,
                       total_slots, max_nlp, fs=fs, no_terms_label=no_terms_label)
        ax_enrich_list[0].set_xlabel(right_xlabel, fontsize=fs["axis"])
        ax_enrich_list[1].set_xlabel(right_xlabel, fontsize=fs["axis"])
        if right_title:
            ax_enrich_list[0].set_title(right_title, loc="left", fontsize=fs["title"])
    else:
        for ci, (ax_e, (col_layout, col_total)) in enumerate(zip(ax_enrich_list, enrich_cols)):
            _render_enrich_column(ax_e, col_layout, total_slots, max_nlp,
                                  fs=fs, no_terms_label=no_terms_label)
            ax_e.set_xlabel(right_xlabel, fontsize=fs["axis"])
            if right_title and ci == 0:
                ax_e.set_title(right_title, loc="left", fontsize=fs["title"])

    axes = {"heatmap": ax_heat, "labels": ax_labels,
            "enrichment": ax_enrich_list[0], "enrichment_cols": ax_enrich_list,
            "colorbar": cax}
    return fig, axes


__all__ = [
    "plot_oydeg_heatmap_enrichment",
    "darken",
    "cap_gene",
    "WARM",
    "COOL",
]
