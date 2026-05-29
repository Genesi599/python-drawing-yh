# `drawing-yh` 科研出图标准

本包的 `__init__.py` 在 `import drawing_yh` 时**自动应用**以下标准。
源约定见 `LLM-yh/LLM_textbook/sci_plot_subplot_notes.md`,本文件是包级落地说明。

## 一句话用法

```python
from drawing_yh import save_fig, OKABE_ITO, compute_figsize, DOUBLE_COL_IN
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=compute_figsize(n_items=len(df), width=DOUBLE_COL_IN))
# ... 用 OKABE_ITO / DEEP_20 上色 ...
save_fig(fig, 'out/fig1.pdf', also=('.png', '.svg'))   # 一次写三份,各自正确 dpi
```

---

## 字体(可编辑性,投稿必需)

`__init__.py` 包级生效:

```python
'pdf.fonttype':       42,         # TrueType,Illustrator 可编辑
'ps.fonttype':        42,
'svg.fonttype':       'none',     # SVG 文字不转曲
'pdf.use14corefonts': False,      # 不替换为不可编辑核心字体
'font.family':        'Arial',    # 期刊通用
```

## 字号

- 项目统一 **`font_size = 8` pt**(印刷字号 5–8pt 范围;本项目取 8)
- `font.size` / `axes.titlesize` / `axes.labelsize` / `xtick.labelsize` / `ytick.labelsize` / `legend.fontsize` / `figure.titlesize` 全部默认 8
- **不要靠缩放调字号** —— 缩放会失真
- 所有文字一致(刻度、轴名、annotation、title)

## 输出尺寸 / DPI

| 后缀     | dpi    | 备注                                                       |
|--------|------|----------------------------------------------------------|
| `.svg` | **72** | SVG 字号 = `font_size × 72/dpi`,**只有 72 才能 1:1** |
| `.pdf` | _none_ | 矢量,dpi 没意义                                              |
| `.png` | 300+   | 期刊印刷质量普遍要求 ≥ 300                                       |

统一通过 `save_fig(fig, path, also=(…))` 写出 —— 内部按后缀自动选 dpi,
默认带 `bbox_inches='tight'` + `facecolor='white'` + 去 metadata
(`Creator` / `Producer` 都清掉,期刊提交不会带 matplotlib / 系统签名)。

## Figure size

| 期刊版面               | 常量                  | 英寸     |
|------------------------|----------------------|--------|
| 单栏 ~ 8.5 cm          | `SINGLE_COL_IN`      | 3.35   |
| 1.5 栏 ~ 11.4 cm       | `ONE_HALF_COL_IN`    | 4.49   |
| 双栏 ~ 17.5 cm         | `DOUBLE_COL_IN`      | 6.89   |

- `figsize` 直接 = 目标印刷尺寸,**不后期缩放**
- **不要预先固定**:用 `compute_figsize(n_items, width=…)` 按字号 + 内容动态算
- 高度公式默认 `per_item_h * n + base_h`,clip 到 `[min_h, max_h]`,
  来自 `scatter/dot_chart` 实战调出来的值

## Autoshrink figsize(2026-05-16 起,所有图型推荐)

"版面尽量紧凑,只要文字不重叠就好"是统一的出图原则。`autoshrink_figsize` 自动找最小 figsize,
从初始大小开始缩,实测 label bbox 重叠就回退 / 放大。不需要手算 figsize。

```python
from drawing_yh import autoshrink_figsize

def render(figsize):
    fig, ax = plt.subplots(figsize=(figsize, figsize))
    ax.bar(...)
    ax.set_xticks(...); ax.set_xticklabels(labels, rotation=45)
    return fig, ax

fig, ax = autoshrink_figsize(
    render,
    initial=4.0, min_size=2.0,
    get_texts=lambda ax: ax.get_xticklabels(),   # 跟图型相关:bar / scatter / heatmap 用 ticklabels;chord / network 用 ax.texts(default)
    bbox_shrink=0.5,                              # bbox 收缩 50% 才算重叠(允许 labels 物理接近,matplotlib bbox 含 padding)
)
```

**适用场景**:
- chord diagram (`drawing_yh.chord.chord_diagram` 已用),sector label 重叠
- bar chart 长 x 轴 label 容易重叠
- heatmap 多行 / 多列 label
- scatter 多 annotation
- network node label

**Why**:手动调 figsize 不准 — 同一 `figsize=4` 在 6 个 cell 跟 24 个 cell 上视觉完全不同。
让代码实测 label bbox 决定最小可用 figsize。

## 文字防重叠

不靠经验估,用 pixel bbox 测:

```python
from drawing_yh import nudge_no_overlap_y
nudge_no_overlap_y(fig, label_texts, step=0.005)
```

底层:`get_text_bboxes(fig, texts)` 取每个 `Text` 的 `get_window_extent()`,
`find_overlaps(bboxes)` 返回所有重叠对。

## 颜色 —— 同类别跨子图保持一致

| 名字              | 用途                                                |
|------------------|---------------------------------------------------|
| `OKABE_ITO`       | 8 色色盲友好通用 categorical(Okabe-Ito 标准)         |
| `DEEP_20`         | 20 色深色饱和(项目历史色板,灰度可区分性已验证)        |
| `AGE_GRADIENT`    | 5 色 老 → 嫩 渐变(viridis 风格)                  |
| `SEX`             | `{F: 橙, M: 紫}`(色盲友好二色对比)                  |
| `SPECIES`         | `{human, monkey, mouse, rat}` 四色                  |
| `THREE_GROUP`     | young / mid / old 三色(蓝 → 灰 → 红)              |

色盲 / 灰度检查工具:

```python
from drawing_yh import grayscale_distinct, preview_palette
grayscale_distinct(my_palette)                   # bool,灰度下两两亮度差 ≥ 0.1?
preview_palette(my_palette, save_as='preview.png')  # 双行预览:原色 vs 灰度
```

## 单细胞作图模板(dot plot / feature plot / heatmap + 富集)

三个单细胞常用图的**纯 matplotlib 模板**(8pt / 配合 `save_fig` 三格式 / 不依赖 scanpy·gseapy)。
**分层约定**:drawing-yh 只管画(输入已算好的矩阵 / 坐标 / 表),h5ad 读取、选基因、
算均值 / z-score、跑富集等数据准备留在 `single_cell-yh`,由它调这里的模板。

| 图型 | 入口 | 输入 |
|---|---|---|
| Marker dot plot | `marker_dotplot(data, row_order, gene_order, …)` | long 表(row / gene / avg / pct) |
| Embedding feature plot | `feature_plot(coords, values, …)` | coords (N,2) + 值矩阵 / dict / DataFrame |
| 热图 + 每行富集条 | `heatmap_with_row_bars(Z, row_labels, …)` | Z 矩阵 + 每行 `[(term, value)]` |

**Dot plot**(默认 grey-red,scseq 配色;`cmap="viridis"` 切 viridis;图例默认右侧竖排——colorbar 右上 + size 右下,`legend_loc="bottom"` 回到底部横排):

```python
from drawing_yh import marker_dotplot
fig, ax, sc = marker_dotplot(
    long_df, row_order=celltypes, gene_order=genes,
    scale='row',                 # 'gene'(默认,per-gene) / 'row'(=scanpy standard_scale='group') / 'none'
    # cmap 默认 grey->red(scseq);要 viridis 传 cmap="viridis"
    row_colors=celltype_colors,  # 左侧 cell-type 彩点(dict / seq)
    block_per_gene=gene_block,   # 列块竖线分隔(可选)
)
```

**Feature plot**(多基因 UMAP/tSNE 网格,共享 grey-red colorbar + 角落箭头轴):

```python
from drawing_yh import feature_plot
fig, axes = feature_plot(
    coords,                      # adata.obsm['X_umap']
    {'MZB1': v1, 'PTX3': v2},    # 或 DataFrame(列=基因) / (N,G) array + genes=
    vmin='p2', vmax='p98',       # 百分位(默认按非零值)或固定值
    share_clim=True,             # 所有 panel 共用一个 colorbar
    axis_labels=('UMAP1', 'UMAP2'),
)
```

**Heatmap + 富集条**(左 z-score 热图 + 列块分隔 + 左侧行彩条;右每行 top 富集条):

```python
from drawing_yh import heatmap_with_row_bars
fig, (ax_heat, ax_bars) = heatmap_with_row_bars(
    Z, row_labels=celltypes,
    row_bars={ct: [(term, neg_log10_p), …] for ct in celltypes},
    block_sizes=per_group_marker_counts,   # 列块竖线
    row_colors=celltype_colors,            # 左侧彩条 + 柱色
    z_clip=2.0,
)
```

原语(自拼布局时用):dot plot `dot_sizes` / `add_row_color_dots` / `add_block_separators`;
feature plot `scatter_embedding` / `add_embedding_axes` / `resolve_vlim` / `auto_ncols`。

## HTML 报告(`drawing_yh.report`)

Markdown 研究报告 → 标准 GitHub-flavored 单页 / 多页 HTML,
样式规则统一,所有项目复用同一套 CSS。**底层仍走 pandoc**,本模块只
封一层薄壳 + ship 一份默认 CSS。

```python
from drawing_yh import report

# 1) 单页:REPORT.md → REPORT.html(同目录,自动复制 default.css)
report.render('REPORT.md')

# 2) 多页拆分 + 跨页 nav(适合长报告分发 / 截图)
md = open('REPORT.md', encoding='utf-8').read()
chunks = report.split_by_h2(md, [
    ('00_summary',  ['摘要', '项目目标', 'Layer 1']),
    ('01_findings', ['Finding A', 'Finding B', '数字总览']),
    ('02_appendix', ['附录 A', '附录 B', '致谢']),
])
PAGES = [('00_summary','摘要'), ('01_findings','Findings'), ('02_appendix','附录')]
for name, body in chunks.items():
    report.render(body,
                  out_html=f'pages/{name}.html',
                  title=name,
                  header_html=report.build_nav(PAGES, current=name))

# 3) 仅取 CSS 路径(配合 Makefile / 手写 pandoc 命令)
print(report.css_path())   # → /…/drawing_yh/report/templates/default.css
```

**约定**

- 默认 CSS:`templates/default.css`(浅色)/ `templates/dark.css`(暗色 GitHub Dark)——
  980 px 居中列、PingFang/YaHei 中文回退、含 `.report-nav` 顶部 sticky 导航 +
  `aside.report-sidebar` 左侧 sidebar 样式(窄屏 ≤ 1180 折成 48 px 细栏 + hover 展开)
- `render()` 默认 `copy_css=True`:把 CSS 复制到输出旁,链接走相对路径
  —— 否则 `file://` 跨源会被浏览器挡(打开 HTML 全是裸文本)
- `render()` 内部用 pandoc `--metadata pagetitle=` 而不是 `title=`,避免在正文重复渲染 `<h1 class="title">`
- `split_by_h2()` 支持模糊匹配(`'Layer 1'` 命中 `'## Layer 1 · 蛋白 LR…'`),
  默认把 H1 标题前置到每个分页(可关 `h1_passthrough=False`)
- `build_nav` / `build_sidebar` 都接受 `(slug, label)` 或 `(slug, label, href)`,
  以及 `('section', 'group_name')` 作 sidebar 分组头

## Slide-style 多页报告(默认结果报告形式)

3 个模板文件配套使用,从 `drawing_yh.report` 取:

| Helper | 模板文件 | 用途 |
|---|---|---|
| `report.config_template_path()` | `_report_config_template.py` | **共享配置**(SLIDES / CATEGORIES / split_by_heading / load_chunks),纯 import 无副作用 |
| `report.slide_template_path()` | `build_pages_template.py` | **HTML 渲染**(REPORT.md → 多页 dark 网页 + sticky nav + sidebar) |
| `report.combined_config_template_path()` | `_combined_report_config_template.py` | **三层 combined 配置**(多个已有报告 → 顶层 super switch + 原 category nav + sidebar) |
| `report.combined_template_path()` | `build_combined_report_template.py` | **三层 combined 渲染**(非破坏式读取多个子报告,统一生成 combined/pages) |
| `report.pptx_template_path()` | `md_to_pptx_template.py` | **PPTX 渲染**(REPORT.md → native PPTX,继承一份母版 .pptx 的 theme) |

**典型工作流**:

```bash
# 1) 拷三份模板进项目
cp $(python -c "from drawing_yh import report; print(report.config_template_path())") _report_config.py
cp $(python -c "from drawing_yh import report; print(report.slide_template_path())")  build_pages.py
cp $(python -c "from drawing_yh import report; print(report.pptx_template_path())")   md_to_pptx.py

# 2) 改 _report_config.py 的 SLIDES + CATEGORIES;改 md_to_pptx.py 顶端 4 个常量
#    (PROJECT_LABEL / TEMPLATE / OUT / 标题页文案)

# 3) 跑
python build_pages.py     # → pages/*.html(只跑 HTML)
python md_to_pptx.py      # → REPORT.pptx(只跑 PPT,不连带刷 HTML)
```

**三层目录 combined 工作流**(多个已存在子报告合成一个网页):

```bash
mkdir combined
cd combined
cp $(python -c "from drawing_yh import report; print(report.combined_config_template_path())") _combined_report_config.py
cp $(python -c "from drawing_yh import report; print(report.combined_template_path())")        build_combined_report.py
cp $(python -c "from drawing_yh import report; print(report.serve_nocache_template_path())")    serve_nocache.py

# 改 _combined_report_config.py:
# SOURCE_REPORTS = [{'key':'m','label':'脑膜','path':HERE.parent/'meninges'}, ...]
# CROSS_PAGES 按需放综合/跨组织页面

python build_combined_report.py
python serve_nocache.py 8775
```

三层结构 = 顶层 super switch(如组织/物种/模块) + 每个子报告原有 category nav + 左侧 sidebar。
构建是非破坏式,不改各子报告 `REPORT.md`;只读取各子报告 `_report_config.py` / `REPORT.md`,
并把各自 `report_figs/` 加前缀同步到 combined/report_figs/。

**产物布局**:

```
analysis/<project>/
├── REPORT.md             # H2 / H3 划章节
├── _report_config.py     # SLIDES / CATEGORIES,两个脚本共用
├── build_pages.py        # HTML 渲染入口
├── md_to_pptx.py         # PPT 渲染入口
├── index.html            # redirect → pages/index.html
├── pages/
│   ├── index.html        # landing
│   ├── <slug>.html × N   # 每张 ≤ 一页 PPT 量
│   └── dark.css          # 自动复制
└── (D:/Projects/<project>/REPORT.pptx)   # PPT 落 D 盘(代码层不存图/数据)
```

**HTML 报告特征**:dark GitHub 风 / sticky 顶 4-tab(无前缀文字)/ 左 sidebar 按当前 cat 过滤
(landing 例外)/ 窄屏 sidebar 折 48 px 细栏不挪顶 / 每页只一层标题(用 pandoc `pagetitle`)。

**PPTX 特征**:继承母版 master / theme / 字体配色;每张 slide 顶部 = 项目标签 + 标题 + 橙色分隔线;
有图 → 左大图 + 右 ➤ bullets(垂直居中);无图 → 全宽;表格列宽按内容长短分配;字号 14pt
header / 11–14pt body 自动降级;**有图的 slide 不放表;表 + bullets 不共存**;放不下静默丢。

详细约定见 `preferences.md` 的 "结果报告(slide-style HTML)的标准形式" 节。

**前置依赖**

- `pandoc` 必须在 PATH 上 —— Windows: `winget install JohnMacFarlane.Pandoc`,
  macOS: `brew install pandoc`

**发布到 reports.yhtiddly.fun**

组内共享的 slide-style HTML 报告发布到 `https://reports.yhtiddly.fun/<project>/`,
服务器目录为 `todo:/var/www/reports/<project>/`。站点由 nginx basic auth 保护,
账号密码不写入仓库或记忆库。

```python
from drawing_yh import report

# 在含 pages/、report_figs/、index.html 的报告目录运行
report.publish("cellchat_x_BMIF")
```

`publish()` 会上传 `pages/`、`report_figs/` 和根 `index.html`,并把远程目录/文件
权限修成 nginx 可读。新项目发布后用无账号 401、有账号 200 的 curl 检查再发链接。

## 复用性

- 颜色列表、字号、figsize **变量集中管理**,不硬编码在多处
- 代码里多处出现魔法数(字号、宽高、颜色字面量)就抽常量
- **项目级共识应该沉淀进 `drawing_yh`**(rc / palettes / save_fig / layout / report),
  而不是每个项目脚本各写一遍

## 何时绕开标准

包级 rc 是默认值,脚本里可以局部覆盖(例如 poster 字号大一些):

```python
import matplotlib as mpl
mpl.rcParams['font.size'] = 12     # 这张图我要大一点
# ...
from drawing_yh import set_rc
set_rc()                            # 画完别人时恢复
```

或一次性:

```python
import drawing_yh
drawing_yh.set_rc({'font.size': 12})   # 临时整体覆盖
```
