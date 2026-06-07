"""Enrichment dotplot — scatter matrix: rows=pathway terms, cols=cell types.

Colour encodes -log10(padj); size encodes gene-overlap count.
Designed for KEGG/GO over-representation results in the OYdeg analysis
pipeline but fully generic: pass any dataframe with the right columns.

Public API
----------
plot_enrichment_dotplot(df, *, columns, ...)
    Returns (fig, ax).
"""
from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

import drawing_yh  # noqa: F401  — apply rc on import
from drawing_yh.io import save_fig  # noqa: F401 (re-used by callers)

__all__ = ["plot_enrichment_dotplot", "pick_enrich_terms"]

# ── noise-term exclusion list shared with retina / brain scripts ──────────────
DEFAULT_EXCLUDE_KW = [
    "coronavirus", "covid", "epstein-barr", "herpes", "hepatitis", "influenza",
    "measles", "tuberculosis", "amoebiasis", "trypanosomiasis", "leishmaniasis",
    "malaria", "toxoplasmosis", "staphylococcus", "salmonella", "shigellosis",
    "pertussis", "legionellosis", "yersinia", "papillomavirus",
    "carcinoma", "cancer", "melanoma", "glioma", "leukemia",
    "cardiomyopathy", "myocarditis", "alcoholism", "cocaine", "amphetamine",
    "morphine", "nicotine", "prion disease", "systemic lupus",
    "graft-versus-host", "allograft rejection", "asthma",
    "diabetes", "diabetic", "arthritis", "inflammatory bowel",
    "lipid and atherosclerosis", "bacterial invasion",
    "aldosterone", "cortisol", "cushing", "renin", "insulin secretion",
    "pancreatic secretion", "protein digestion", "parathyroid",
    "prolactin", "growth hormone", "thyroid hormone", "gnrh", "oxytocin",
    "gastric acid", "salivary", "carbohydrate digestion", "proximal tubule",
    "sodium reabsorption", "cardiomyocyte", "cardiac muscle",
    "cardiac conduction", "heart rate", "neuromuscular", "osteoclast",
    "taste transduction", "olfactory transduction", "spinocerebellar",
    "immunodeficiency virus", "trp channels",
]


def _keep_term(term: str, exclude: list[str]) -> bool:
    tl = term.lower()
    return not any(k in tl for k in exclude)


def _dedup_jaccard(df: pd.DataFrame, genes_col: str, thresh: float) -> pd.DataFrame:
    """Greedy Jaccard dedup: keep rows whose gene set is < thresh similar to all kept rows."""
    kept_idx, kept_sets = [], []
    for i, row in df.iterrows():
        gs = set(str(row[genes_col]).split(";"))
        if all(
            (len(gs & ks) / len(gs | ks) if (gs | ks) else 0) < thresh
            for ks in kept_sets
        ):
            kept_idx.append(i)
            kept_sets.append(gs)
    return df.loc[kept_idx]


def pick_enrich_terms(
    df: pd.DataFrame,
    *,
    subtype_col: str = "subtype",
    term_col: str = "Term",
    padj_col: str = "Adjusted P-value",
    hits_col: str = "hits",
    genes_col: str = "Genes",
    term_size_col: str = "term_size",
    padj_thresh: float = 0.05,
    max_term_size: int = 300,
    min_hits: int = 3,
    top_terms: int = 16,
    jaccard_thresh: float = 0.6,
    exclude_keywords: list[str] | None = None,
) -> tuple[list[str], list[str]]:
    """Select top enrichment terms and column order.

    Returns (terms, columns) — both in display order.
    ``columns`` is sorted by descending total enrichment score.
    """
    kw = exclude_keywords if exclude_keywords is not None else DEFAULT_EXCLUDE_KW
    d = df[
        (df[padj_col] < padj_thresh)
        & (df[term_size_col] <= max_term_size)
        & (df[hits_col] >= min_hits)
        & df[term_col].map(lambda t: _keep_term(t, kw))
    ].copy()
    if d.empty:
        return [], []
    d["_score"] = -np.log10(d[padj_col].clip(lower=1e-300))
    agg = (
        d.groupby(term_col)
        .agg(
            _n=(subtype_col, "nunique"),
            _s=("_score", "sum"),
            **{genes_col: (genes_col, lambda x: ";".join(x.astype(str)))},
        )
        .sort_values(["_n", "_s"], ascending=False)
        .reset_index()
    )
    top_df = _dedup_jaccard(agg, genes_col, jaccard_thresh).head(top_terms)
    terms = top_df[term_col].tolist()
    cols = (
        d[d[term_col].isin(terms)]
        .groupby(subtype_col)["_score"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    )
    return terms, cols


def plot_enrichment_dotplot(
    df: pd.DataFrame,
    *,
    columns: list[str],
    terms: list[str] | None = None,
    subtype_col: str = "subtype",
    term_col: str = "Term",
    padj_col: str = "Adjusted P-value",
    hits_col: str = "hits",
    genes_col: str = "Genes",
    term_size_col: str = "term_size",
    direction_color: str = "#B2182B",
    padj_thresh: float = 0.05,
    max_term_size: int = 300,
    min_hits: int = 3,
    top_terms: int = 16,
    jaccard_thresh: float = 0.6,
    exclude_keywords: list[str] | None = None,
    padj_max_log: float = 6.0,
    title: str | None = None,
    term_max_len: int = 46,
    column_display: dict[str, str] | None = None,
) -> tuple[plt.Figure, plt.Axes] | None:
    """Draw an enrichment bubble/scatter dotplot.

    Parameters
    ----------
    df :
        Enrichment results, one row per (term, cell-type). Must contain at
        minimum: *subtype_col*, *term_col*, *padj_col*, *hits_col*,
        *genes_col*, *term_size_col*.
    columns :
        Ordered list of cell-type / subtype names for the x-axis.
        Rows for cell types not in *columns* are ignored.
    terms :
        Pre-selected term list (rows). If None, auto-selected via
        :func:`pick_enrich_terms`.
    direction_color :
        Hex colour for the high end of the -log10(padj) colour ramp
        (low end is light grey ``#e9e9e9``).
    column_display :
        Optional ``{raw_name: display_label}`` mapping applied to x-tick labels.

    Returns
    -------
    (fig, ax) or None if no terms pass the filter.
    """
    kw = exclude_keywords if exclude_keywords is not None else DEFAULT_EXCLUDE_KW
    d = df[df[subtype_col].isin(columns)].copy()

    if terms is None:
        terms, _ = pick_enrich_terms(
            d,
            subtype_col=subtype_col,
            term_col=term_col,
            padj_col=padj_col,
            hits_col=hits_col,
            genes_col=genes_col,
            term_size_col=term_size_col,
            padj_thresh=padj_thresh,
            max_term_size=max_term_size,
            min_hits=min_hits,
            top_terms=top_terms,
            jaccard_thresh=jaccard_thresh,
            exclude_keywords=kw,
        )

    if not terms:
        return None

    # Only keep columns that actually have at least one sig term
    cts = [c for c in columns if c in d[d[term_col].isin(terms) & (d[padj_col] < padj_thresh)][subtype_col].unique()]
    if not cts:
        return None

    cmap = LinearSegmentedColormap.from_list("enr", ["#e9e9e9", direction_color])
    norm = Normalize(-np.log10(padj_thresh), padj_max_log)

    nr, nc = len(terms), len(cts)
    fig_w = max(4.5, 2.9 + nc * 0.34)
    fig_h = max(2.6, 0.7 + nr * 0.24)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sub = d[d[term_col].isin(terms) & d[subtype_col].isin(cts)]
    hmax = max(float(sub[hits_col].max()), 1.0)

    for _, r in sub.iterrows():
        if float(r[padj_col]) >= padj_thresh:
            continue
        ct = r[subtype_col]
        if ct not in cts:
            continue
        i = cts.index(ct)
        j = terms.index(r[term_col])
        lp = -np.log10(max(float(r[padj_col]), 1e-300))
        ax.scatter(
            i, j,
            s=12 + 70 * (float(r[hits_col]) / hmax),
            color=cmap(norm(lp)),
            edgecolor="white",
            linewidth=0.4,
            zorder=3,
        )

    def _short(t: str) -> str:
        t = t.strip()
        return t if len(t) <= term_max_len else t[: term_max_len - 1] + "…"

    xt_labels = [column_display.get(c, c) if column_display else c for c in cts]
    ax.set_xticks(range(nc))
    ax.set_xticklabels(xt_labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(nr))
    ax.set_yticklabels([_short(t).title() for t in terms], fontsize=7)
    ax.set_ylim(-0.6, nr - 0.4)
    ax.set_xlim(-0.6, nc - 0.4)
    ax.invert_yaxis()
    ax.tick_params(length=0)
    ax.grid(True, color="#f0f0f0", lw=0.5, zorder=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    if title:
        ax.set_title(title, loc="left", fontsize=8, pad=6)

    # Colorbar
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.022, pad=0.02, extend="max")
    cb.set_label("-log10(padj)", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    # Size legend
    for h in [3, max(3, int(hmax * 0.5)), int(hmax)]:
        ax.scatter(
            [], [],
            s=12 + 70 * (h / hmax),
            color="#888888",
            edgecolor="white",
            linewidth=0.4,
            label=str(h),
        )
    ax.legend(
        title="gene hits",
        loc="upper left",
        bbox_to_anchor=(1.0, 1.0),
        fontsize=6,
        title_fontsize=6,
        frameon=False,
        labelspacing=1.0,
    )

    fig.tight_layout()
    return fig, ax
