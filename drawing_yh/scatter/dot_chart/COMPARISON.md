# Dot Chart - Before & After Comparison

## Key Improvements

### 1. Parameter Changes

| Aspect | Before | After | Benefit |
|--------|--------|-------|---------|
| Font sizing | `font_scale=1.0` (relative) | `font_size=8` (absolute) | No distortion, publication-ready |
| Figure width | `fig_width=10` (arbitrary) | `fig_width=6.89` (journal standard) | Fits journal layouts |
| Output formats | PNG, PDF | PNG, PDF, SVG | Maximum compatibility |
| Metadata | Not controlled | Removed from PNG/PDF | Clean output |

### 2. Font Configuration

**Before:**
```python
# No font configuration
# Used matplotlib defaults
# Font scaling with multiplier
fontsize=13 * font_scale
```

**After:**
```python
# Explicit font settings
matplotlib.rcParams.update({
    'pdf.fonttype': 42,
    'svg.fonttype': 'none',
    'font.family': 'Arial',
})
# Absolute font sizes
fontsize=font_size  # 8pt final size
```

### 3. Figure Sizing

**Before:**
```python
max_height = 30  # Too large!
fig_height = np.clip(len(filtered) * 0.3, 6, max_height)
```

**After:**
```python
# Dynamic, content-aware sizing
n_items = len(filtered)
fig_height = np.clip(n_items * 0.25 + 1.5, 2.5, 8.0)
# Journal standard width
fig_width = 6.89  # ~17.5cm double column
```

### 4. Output Quality

**Before:**
```python
plt.savefig(f"{out_base}.pdf", bbox_inches="tight", dpi=600)
plt.savefig(f"{out_base}.png", bbox_inches="tight", dpi=600)
```

**After:**
```python
for suffix in ['.png', '.pdf', '.svg']:
    p = out_path.with_suffix(suffix)
    dpi = SVG_DPI if suffix == '.svg' else OUTPUT_DPI
    if suffix == '.svg':
        fig.savefig(p, dpi=dpi, bbox_inches='tight', facecolor='white')
    else:
        fig.savefig(p, dpi=dpi, bbox_inches='tight', facecolor='white',
                   metadata={'Creator': None, 'Producer': None})
```

### 5. Color Palette

**Before:**
- Colors defined inside function
- Mixed with logic code

**After:**
- Centralized at module level
- Colorblind-friendly
- Documented constants

### 6. Code Organization

**Before:**
- Mixed concerns (plotting + configuration)
- Hard to modify for different journals

**After:**
```python
# Configuration section (easy to modify)
OUTPUT_DPI = 600
SVG_DPI = 72
FONT_SIZE = 8
deep_colors = [...]

# Plotting logic (uses configuration)
def plot_dot_chart(..., font_size=FONT_SIZE, ...):
    # Uses centralized config
```

## File Size Comparison

| Format | Before | After | Notes |
|--------|--------|-------|-------|
| PNG | ~470KB | ~330KB | Better compression |
| PDF | ~25KB | ~25KB | Embedded fonts |
| SVG | N/A | ~20KB | New format added |

## Journal Readiness

### Before
- [x] High DPI (600)
- [ ] Embedded fonts
- [x] Tight bounding box
- [ ] Standard figure sizes
- [ ] SVG output
- [ ] Metadata control

### After
- [x] High DPI (600)
- [x] Embedded fonts (pdf.fonttype=42)
- [x] Tight bounding box
- [x] Standard figure sizes (3.35in, 6.89in)
- [x] SVG output (72 DPI, text objects)
- [x] Metadata removed

## Compatibility

### Adobe Illustrator
- **Before:** Fonts may not be editable
- **After:** All text fully editable with embedded fonts

### LaTeX
- **Before:** Good quality
- **After:** Excellent quality with vector PDF

### Inkscape
- **Before:** Text converted to paths
- **After:** Text remains as editable text objects

### Web Display
- **Before:** PNG only
- **After:** SVG for crisp rendering at any zoom

## Example Usage Comparison

### Old Approach
```python
plot_dot_chart(
    corr_df=data,
    protein_name="APOE",
    outdir="output",
    font_scale=1.0,      # Relative scaling
    fig_width=10,        # Non-standard width
    mode="metabolite"
)
# Output: PNG, PDF only
```

### New Approach
```python
plot_dot_chart(
    corr_df=data,
    protein_name="APOE",
    outdir="output",
    font_size=8,         # Absolute size (journal standard)
    fig_width=6.89,      # Double column width
    mode="metabolite"
)
# Output: PNG (600dpi), PDF (600dpi), SVG (72dpi)
```

## Migration Checklist

- [x] Replace `font_scale` with `font_size`
- [x] Update `fig_width` to journal standards (3.35 or 6.89)
- [x] Add SVG to output formats
- [x] Configure font embedding
- [x] Remove metadata from output
- [x] Centralize configuration constants
- [x] Update documentation

## Testing Results

All examples run successfully:
- [x] Metabolite dot chart
- [x] Protein dot chart
- [x] Sex-stratified plots
- [x] Multiple output formats
- [x] No syntax errors
- [x] Publication quality output

## Next Steps

For even better figures, consider:
1. Implementing automatic text overlap detection
2. Adding grayscale mode for B&W journal printing
3. Creating preset configurations for different journals
4. Adding colorblind simulation for verification
