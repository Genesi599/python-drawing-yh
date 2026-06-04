"""
drawing_yh — Python 科研作图工具包

**import 即应用项目科研出图标准**(可编辑字体、Arial、SVG 不转曲、统一字号),
并导出常用工具:`save_fig`、命名色板、`compute_figsize`、文字防重叠 helper。

完整标准见包根 `STANDARDS.md`,或 `LLM-yh/LLM_textbook/sci_plot_subplot_notes.md`。

Quick start
-----------
    from drawing_yh import save_fig, OKABE_ITO, compute_figsize, DOUBLE_COL_IN
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=compute_figsize(n_items=len(df),
                                                   width=DOUBLE_COL_IN))
    # ... 用 OKABE_ITO 或 DEEP_20 上色 ...
    save_fig(fig, 'out/fig1.pdf', also=('.png', '.svg'))   # 一次写三份
"""
from importlib.metadata import version, PackageNotFoundError
import matplotlib as _mpl

try:
    __version__ = version("drawing-yh")
except PackageNotFoundError:
    __version__ = "0.0.0"


# ============================================================
# 项目常量
# ============================================================
DEFAULT_FONT_SIZE = 8        # 印刷字号(5–8pt 范围,统一 8)
SVG_DPI = 72                 # SVG 必须 72 才能 1:1 渲染文字
PNG_DPI = 300                # PNG 默认印刷质量

SINGLE_COL_IN = 3.35         # 单栏 ~ 8.5 cm
ONE_HALF_COL_IN = 4.49       # 1.5 栏 ~ 11.4 cm
DOUBLE_COL_IN = 6.89         # 双栏 ~ 17.5 cm


# ============================================================
# rcParams —— 包级一次性设置(import 时生效)
# ============================================================
RC_DEFAULTS = {
    # —— 字体可编辑性(投稿必需)——
    'pdf.fonttype':         42,         # 嵌入 TrueType,Illustrator 可编辑
    'ps.fonttype':          42,
    'svg.fonttype':         'none',     # SVG 文字保持文字对象,不转曲
    'pdf.use14corefonts':   False,      # 不替换为不可编辑的核心字体
    # —— 字体 ——
    'font.family':          'Arial',    # 期刊通用;装不上会回退到默认 sans-serif
    'font.size':            DEFAULT_FONT_SIZE,
    'axes.titlesize':       DEFAULT_FONT_SIZE,
    'axes.labelsize':       DEFAULT_FONT_SIZE,
    'xtick.labelsize':      DEFAULT_FONT_SIZE,
    'ytick.labelsize':      DEFAULT_FONT_SIZE,
    'legend.fontsize':      DEFAULT_FONT_SIZE,
    'figure.titlesize':     DEFAULT_FONT_SIZE,
    # —— 保存默认 ——
    'savefig.bbox':         'tight',
    'savefig.pad_inches':   0.02,
    'savefig.facecolor':    'white',
}


def set_rc(overrides: dict = None) -> None:
    """
    应用项目标准 rcParams。

    包 import 时已自动调用一次。如果脚本里临时改了 rc(比如调字号),
    想恢复到标准,可以再 `drawing_yh.set_rc()` 一下。
    """
    rc = dict(RC_DEFAULTS)
    if overrides:
        rc.update(overrides)
    _mpl.rcParams.update(rc)


# 包加载时立即应用一次;所有 `import drawing_yh` / `from drawing_yh import …`
# 的脚本都自动享受标准 rc,不必每个文件再写一遍。
set_rc()


# ============================================================
# 顶层导出
# ============================================================
from .io import save_fig                                # noqa: E402
from .palettes import (                                 # noqa: E402
    OKABE_ITO, DEEP_20, AGE_GRADIENT, SEX, SPECIES, THREE_GROUP,
    to_grayscale, grayscale_distinct, preview_palette, family_palette,
)
from .layout import (                                   # noqa: E402
    compute_figsize, get_text_bboxes, find_overlaps, nudge_no_overlap_y,
    autoshrink_figsize,
    TextBBox,
)
from .dumbbell import render_mean_dumbbell              # noqa: E402
from .venn import venn_diagram                          # noqa: E402
from .dotplot import (                                  # noqa: E402
    marker_dotplot,
    dot_sizes, pad_axes, style_colorbar,
    add_size_legend, add_block_separators, add_row_color_dots,
)
from . import dotplot                                   # noqa: E402,F401
from .embedding import (                                # noqa: E402
    feature_plot,
    resolve_vlim, scatter_embedding, add_embedding_axes, auto_ncols,
    cluster_atlas, cluster_centroids,
)
from . import embedding                                 # noqa: E402,F401
from .heatmap import heatmap_with_row_bars, clustermap_annot              # noqa: E402
from . import report                                    # noqa: E402,F401
from .network import hub_spoke                          # noqa: E402

__all__ = [
    # version
    "__version__",
    # constants
    "DEFAULT_FONT_SIZE", "SVG_DPI", "PNG_DPI",
    "SINGLE_COL_IN", "ONE_HALF_COL_IN", "DOUBLE_COL_IN",
    "RC_DEFAULTS",
    # rc
    "set_rc",
    # io
    "save_fig",
    # palettes
    "OKABE_ITO", "DEEP_20", "AGE_GRADIENT", "SEX", "SPECIES", "THREE_GROUP",
    "to_grayscale", "grayscale_distinct", "preview_palette", "family_palette",
    # layout
    "compute_figsize", "get_text_bboxes", "find_overlaps",
    "nudge_no_overlap_y", "autoshrink_figsize", "TextBBox",
    # charts
    "render_mean_dumbbell",
    "venn_diagram",
    "marker_dotplot",
    "feature_plot",
    "heatmap_with_row_bars",
    "clustermap_annot",
    # dot-plot style primitives (re-exported from drawing_yh.dotplot)
    "dot_sizes", "pad_axes", "style_colorbar",
    "add_size_legend", "add_block_separators", "add_row_color_dots",
    # embedding feature-plot primitives (re-exported from drawing_yh.embedding)
    "resolve_vlim", "scatter_embedding", "add_embedding_axes", "auto_ncols",
    "cluster_atlas", "cluster_centroids",
    # subpackages
    "dotplot",
    "embedding",
    "report",
    # network
    "hub_spoke",
]
