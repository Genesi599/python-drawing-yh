# Dot Chart Module

Scientific publication-quality dot charts for visualizing correlation analysis results, supporting both metabolite and protein correlation modes.

## Features

- [SUPPORTED] Metabolite correlation dot plots with pathway grouping
- [SUPPORTED] Protein-protein correlation dot plots with cellular component grouping
- [SUPPORTED] Sex stratification analysis
- [SUPPORTED] Pathway ratio annotations
- [SUPPORTED] Automatic grouping and color coding
- [SUPPORTED] Dynamic figure size calculation based on content
- [SUPPORTED] Both positive and negative correlation visualization
- [SUPPORTED] Publication-quality output (PNG, PDF, SVG)

## Scientific Figure Standards

### Font
- **Font family**: Arial (journal standard)
- **Font size**: 8pt base (final print size, no scaling)
- **Font embedding**: TrueType fonts embedded (pdf.fonttype = 42)
- **SVG text**: Kept as editable text objects (svg.fonttype = 'none')

### Output Size & DPI
- **PNG/PDF**: 600 DPI for publication quality
- **SVG**: 72 DPI for 1:1 text rendering
- **Figure size**: Auto-calculated based on content, or journal-standard widths
  - Single column: 3.35 inches (~8.5 cm)
  - Double column: 6.89 inches (~17.5 cm)
- **No metadata**: Removed Creator/Producer metadata for clean output

### Color
- **Colorblind-friendly palette**: 20 distinct colors
- **Consistent mapping**: Same categories use same colors across all subplots
- **Grayscale tested**: Ensures distinguishability in grayscale mode

## Quick Start

### Basic Usage

```python
from drawing_yh.scatter.dot_chart.dot_chart import plot_dot_chart
import pandas as pd

# Load your correlation results
corr_df = pd.read_csv("correlation_results.csv")

# Plot metabolite dot chart with publication settings
plot_dot_chart(
    corr_df=corr_df,
    protein_name="APOE",
    outdir="output/metabolite_plots",
    r_threshold=0.25,
    p_threshold=0.05,
    font_size=8,           # Final print size in points
    fig_width=6.89,        # Double column width in inches
    mode="metabolite"
)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `corr_df` | DataFrame | Required | Correlation analysis results with columns: Metabolite_ID, Compound_Name, r, p, etc. |
| `protein_name` | str | Required | Target protein name for title and filename |
| `outdir` | str | Required | Output directory for saved figures |
| `r_threshold` | float | 0.25 | Correlation coefficient threshold |
| `p_threshold` | float | 0.05 | P-value significance threshold |
| `sex_tag` | str | "" | Sex stratification label (e.g., "Male", "Female") |
| `font_size` | float | 8 | Base font size in points (final print size, no scaling) |
| `fig_width` | float | 6.89 | Figure width in inches (auto-calculated if None) |
| `mode` | str | "metabolite" | Plot mode: "metabolite" or "protein" |
| `ratio_dir` | str | None | Directory containing pathway ratio CSV files |

### Required DataFrame Columns

**For Metabolite Mode:**
- `Metabolite_ID`: Metabolite identifier
- `Compound_Name`: Full metabolite name
- `Compound_Abbr`: Abbreviation (optional)
- `SUB_META_PATHWAY`: Sub-pathway category
- `SUPER_META_PATHWAY`: Super-pathway category
- `r`: Correlation coefficient
- `p`: P-value

**For Protein Mode:**
- `Gene_Name`: Gene name
- `Group`: Protein group/category
- `r`: Correlation coefficient
- `p`: P-value

## Examples

Run the example script to generate sample figures:

```bash
cd c:\Users\yh109\Documents\GitHub\python-drawing-yh\drawing_yh\scatter\dot_chart
python example_dot_chart.py
```

This will generate:
1. Metabolite dot chart with pathway grouping
2. Protein dot chart with cellular component grouping
3. Sex-stratified dot charts

Output files will be saved in the `output/` directory.

## Output Files

The module generates three formats for maximum compatibility:

### Formats
- **PNG** (600 DPI): For documents and presentations
- **PDF** (600 DPI, embedded fonts): For LaTeX and vector editing in Illustrator
- **SVG** (72 DPI, text objects): For web and Inkscape editing

### Filenames
- `{protein_name}_{p_threshold}_dotplot_pos.png/pdf/svg` - Positive correlations
- `{protein_name}_{p_threshold}_dotplot_neg.png/pdf/svg` - Negative correlations

Example filenames:
- `APOE_p0.05_dotplot_pos.png`
- `TNF_p0.05_protein_dotplot_neg.pdf`
- `IL6_p0.05_dotplot_pos.svg`

## Visualization Features

- **Grouping**: Items are grouped by pathway/category with color coding
- **Sorting**: Groups sorted by maximum absolute correlation; items sorted within groups
- **Labels**: Right-side annotations show group name and hit/total counts
- **Axis**: Inverted X-axis for negative correlations
- **Dynamic sizing**: Figure height adjusts based on number of items

## Reproducibility

- **Centralized configuration**: All colors, fonts, and sizes use variables
- **No hardcoded values**: Easy to modify for different journals
- **Deterministic output**: Same input always produces same output
- **Metadata removal**: Clean output without software signatures

## Notes

- Items must pass both p-value and r-threshold filters to be included
- Empty results are skipped with a message
- Figure height dynamically calculated: 0.25in per item + 1.5in base (range: 2.5-8in)
- Color palette uses 20 distinct colors for group differentiation
- All text sizes are in final print points (no scaling distortion)
