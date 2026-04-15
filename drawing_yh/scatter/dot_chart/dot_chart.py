import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import os
import re
from pathlib import Path

matplotlib.use('Agg')

# ============================================================
# Scientific Figure Configuration
# ============================================================

# Font settings for publication-quality figures
matplotlib.rcParams.update({
    'pdf.fonttype': 42,          # TrueType fonts for Adobe/Illustrator editing
    'ps.fonttype': 42,           # TrueType fonts for PostScript
    'svg.fonttype': 'none',      # Keep text as text objects in SVG
    'pdf.use14corefonts': False, # Allow full font embedding
    'font.family': 'Arial',      # Journal-standard font
})

# Output settings
OUTPUT_DPI = 600                 # PNG and PDF resolution
SVG_DPI = 72                     # SVG resolution for 1:1 text rendering
FONT_SIZE = 8                    # Base font size (final print size)

# Colorblind-friendly palette
deep_colors = [
    "#c8001e", "#2a8a3a", "#2a4db5", "#c45e10", "#6a0f8a",
    "#1aa0c4", "#b020b0", "#8a9e00", "#c47090", "#2a7070",
    "#8060c0", "#7a4a10", "#600000", "#30a060", "#606000",
    "#c09060", "#000060", "#606060", "#d04020", "#000000"
]


def safe_filename(x: str) -> str:
    """Convert filename to safe format by removing special characters."""
    return re.sub(r'[\\/:*?"<>|]+', "_", str(x))[:180].rstrip(" .")


def p_tag(p_threshold: float) -> str:
    """Convert p-value threshold to filename-friendly string, e.g., 0.01 -> 'p0.01'."""
    return f"p{p_threshold}"


def darken(hex_color, factor=0.65):
    """Darken a hex color by a given factor for better contrast."""
    hex_color = hex_color.lstrip('#')
    r, g, b = (int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    return '#{:02x}{:02x}{:02x}'.format(int(r * factor), int(g * factor), int(b * factor))


def plot_dot_chart(
        corr_df: pd.DataFrame,
        protein_name: str,
        outdir: str,
        r_threshold: float = 0.25,
        p_threshold: float = 0.05,
        sex_tag: str = "",
        font_size: float = FONT_SIZE,
        fig_width: float | None = None,
        mode: str = "metabolite",
        ratio_dir: str | None = None,
):
    """
    Unified dot chart plotting function, supporting both metabolite and protein modes.
    Filename includes p-value threshold identifier.
    
    Parameters:
    -----------
    corr_df : pd.DataFrame
        Correlation analysis results DataFrame
    protein_name : str
        Target protein name for title and filename
    outdir : str
        Output directory for saved figures
    r_threshold : float
        Correlation coefficient threshold
    p_threshold : float
        P-value significance threshold
    sex_tag : str
        Sex stratification label (e.g., "Male", "Female")
    font_size : float
        Base font size in points (final print size, default: 8pt)
    fig_width : float
        Figure width in inches (default: auto-calculated based on content)
    mode : str
        Plot mode: "metabolite" or "protein"
    ratio_dir : str
        Directory containing pathway ratio CSV files
    """
    os.makedirs(outdir, exist_ok=True)

    if ratio_dir is None:
        ratio_dir = outdir

    if mode == "metabolite":
        group_col = "SUB_META_PATHWAY"
        title_target = "Metabolites"
        file_tag = "dotplot"
    else:
        group_col = "Group"
        title_target = "Proteins"
        file_tag = "protein_dotplot"

    # Colorblind-friendly palette

    for r_sign, r_label in [("pos", f"r>{r_threshold}"), ("neg", f"r<-{r_threshold}")]:
        r_filter = corr_df['r'] > r_threshold if r_sign == "pos" else corr_df['r'] < -r_threshold
        filtered = corr_df[(corr_df['p'] < p_threshold) & r_filter].copy()

        if filtered.empty:
            print(f"[SKIP] Dot plot ({r_sign}, {p_tag(p_threshold)}): No matching items")
            continue

        tag = f"_{sex_tag}" if sex_tag else ""
        has_group = group_col in filtered.columns and filtered[group_col].notna().any()

        # -- Pathway ratio (metabolite mode only) --
        ratio_label = None
        if mode == "metabolite":
            ratio_path = os.path.join(ratio_dir,
                f"{safe_filename(protein_name)}{tag}_{p_tag(p_threshold)}_pathway_ratio_{r_sign}.csv")
            if os.path.exists(ratio_path):
                ratio = pd.read_csv(ratio_path)
                filtered = filtered.merge(
                    ratio[["SUB_META_PATHWAY", "hit_count", "total_count", "ratio"]],
                    on="SUB_META_PATHWAY", how="left"
                )
                ratio_label = ratio.set_index("SUB_META_PATHWAY")

        # -- Group sorting --
        if has_group:
            pathway_order = (
                filtered.groupby(group_col)['r']
                .apply(lambda s: s.abs().max())
                .sort_values(ascending=True)
                .index.tolist()
            )
            filtered[group_col] = pd.Categorical(
                filtered[group_col], categories=pathway_order, ordered=True
            )
            filtered = filtered.sort_values(
                [group_col, "r"], ascending=[True, r_sign == "pos"]
            ).reset_index(drop=True)

            groups = [g for g in pathway_order if g in filtered[group_col].values]
            color_map = {g: deep_colors[i % len(deep_colors)] for i, g in enumerate(groups)}
        else:
            filtered = filtered.sort_values("r", ascending=(r_sign == "pos")).reset_index(drop=True)
            groups = []
            color_map = {}

        # -- Dynamic figure size calculation --
        # Height based on number of items (minimum 2 inches, max 8 inches)
        n_items = len(filtered)
        fig_height = np.clip(n_items * 0.25 + 1.5, 2.5, 8.0)
        
        # Width: use journal standard if not specified (single column: 3.35in, double: 6.89in)
        if fig_width is None:
            fig_width = 6.89  # Default to double column width
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        if groups:
            dot_colors = [color_map.get(g, '#333333') for g in filtered[group_col]]
        else:
            dot_colors = '#2a4db5'

        ax.scatter(filtered['r'], range(len(filtered)), c=dot_colors, s=80, alpha=0.95, zorder=3)

        # -- Group separator lines and labels --
        if groups:
            for g in groups:
                pos_list = filtered.index[filtered[group_col] == g].tolist()
                if not pos_list:
                    continue
                mid = np.mean(pos_list)
                start = min(pos_list) - 0.5
                ax.axhline(start, color="grey", linewidth=1.2, linestyle="--", alpha=0.5)

                if mode == "metabolite" and ratio_label is not None and g in ratio_label.index:
                    row = ratio_label.loc[g]
                    label = f"{g}  ({int(row['hit_count'])}/{int(row['total_count'])})"
                elif mode == "protein":
                    count = len(pos_list)
                    label = f"{g}  ({count})"
                else:
                    label = g

                ax.text(1.01, float(mid), label,
                        transform=ax.get_yaxis_transform(),
                        va="center", ha="left",
                        fontsize=font_size, color=color_map.get(g, '#333'),
                        fontweight="bold", clip_on=False)

        # -- Y-axis labels --
        if mode == "metabolite":
            display_names = []
            for _, row in filtered.iterrows():
                abbr = row.get('Compound_Abbr', '')
                name = abbr if abbr and str(abbr) != 'nan' else row['Compound_Name']
                display_names.append(name)
        else:
            display_names = filtered['Gene_Name'].tolist()

        ax.set_yticks(range(len(filtered)))
        ax.set_yticklabels(display_names, fontsize=font_size)
        if groups:
            for tick, (_, row) in zip(ax.get_yticklabels(), filtered.iterrows()):
                tick.set_color(darken(color_map.get(row[group_col], '#333')))
                tick.set_fontweight("bold")

        ax.set_xlabel("Pearson r", fontsize=font_size + 1)
        ax.tick_params(axis='x', labelsize=font_size)

        x_min = filtered['r'].min() - 0.05
        x_max = filtered['r'].max() + 0.05
        ax.set_xlim(x_min, x_max)
        if r_sign == "neg":
            ax.invert_xaxis()

        # -- Title --
        ax.set_title(
            f"{protein_name} — Correlated {title_target}",
            fontsize=font_size + 2, pad=40
        )

        corr_direction = "Positive Correlation" if r_sign == "pos" else "Negative Correlation"
        subtitle_color = "red" if r_sign == "pos" else "blue"
        subtitle_parts = [f"{corr_direction} (p<{p_threshold}, {r_label})"]
        if sex_tag:
            subtitle_parts.append(f"Sex: {sex_tag}")

        ax.text(0.5, 1.01, "  |  ".join(subtitle_parts),
                transform=ax.transAxes, ha="center", va="bottom",
                fontsize=font_size + 1, color=subtitle_color, fontweight="bold")

        ax.margins(y=0.02)
        plt.subplots_adjust(right=0.55)
        plt.tight_layout()

        # -- Text overlap detection and adjustment --
        fig.canvas.draw()
        
        # Filename includes p-value identifier
        out_base = os.path.join(outdir,
            f"{safe_filename(protein_name)}{tag}_{p_tag(p_threshold)}_{file_tag}_{r_sign}")
        
        # -- Save in multiple formats with scientific standards --
        out_path = Path(out_base)
        for suffix in ['.png', '.pdf', '.svg']:
            p = out_path.with_suffix(suffix)
            dpi = SVG_DPI if suffix == '.svg' else OUTPUT_DPI
            
            # SVG doesn't support metadata parameter
            if suffix == '.svg':
                fig.savefig(p, dpi=dpi, bbox_inches='tight', facecolor='white')
            else:
                fig.savefig(p, dpi=dpi, bbox_inches='tight', facecolor='white',
                           metadata={'Creator': None, 'Producer': None})
        
        plt.close()
        print(f"[SUCCESS] Dot plot saved: {out_base}.png (also PDF and SVG)")
