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

## 复用性

- 颜色列表、字号、figsize **变量集中管理**,不硬编码在多处
- 代码里多处出现魔法数(字号、宽高、颜色字面量)就抽常量
- **项目级共识应该沉淀进 `drawing_yh`**(rc / palettes / save_fig / layout),
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
