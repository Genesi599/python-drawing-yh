"""drawing_yh.assembly — 统一排图工具(2026-09-02 合并:
原 figure_assembly/build_layout.py 的 flow 排版 + MRI-yh/assemble_age_compare.py 的固定网格算法)。

模式(每页 page 级 `mode` 键):
  flow(默认) — 子图保持原始物理尺寸(不缩放, 2026-09-02 用户确认规则)逐行流式排版,
               画布按内容自适应; 超宽面板单独成行。
  grid       — 固定网格(原 assemble_age_compare 算法): cell 尺寸 + 行分组标签 +
               顶部标题, 每格内 thumbnail 式等比 fit 居中(默认只缩不放, enlarge=true 才放大)。
  abs        — 面板显式 x_mm/y_mm 绝对定位(旧布局表兼容)。
  有意缩放: 任何模式下面板显式 width_mm/height_mm 即覆盖自然尺寸(须在 YAML 注释原因)。

用法:
    python -m drawing_yh.assembly.build_layout <layout.yaml> [--ink <inkscape.exe>]

layout YAML:
    name / out
    pages:
      - page: {width_mm, height_mm, background}
        mode: flow | grid | abs
        grid: {cell_w_mm, cell_h_mm, cols?, label_w_mm, title, title_h_mm, enlarge}
        gutters: {label_gap_mm, caption_gap_mm, col_gap_mm, row_gap_mm, margin_mm}
        panels:
          - {id, src, label?, caption?, row?(grid 行分组), x_mm/y_mm(abs), width_mm/height_mm(有意缩放)}

实现要点(踩坑修正):
    - 嵌入不用嵌套 <svg>(Inkscape 1.4 导 PDF 时每个嵌套 svg 变独立页),
      改深度合并: 子 SVG children 搬进 <g translate+scale>, id/url(#id)/href 加 "<pid>__" 前缀。
    - 子 SVG 的 inkscape:page / sodipodi:namedview 多页标记必须丢弃。
    - 根元素勿手写 xmlns(ET.register_namespace 已注册, 手写重复属性)。
"""
