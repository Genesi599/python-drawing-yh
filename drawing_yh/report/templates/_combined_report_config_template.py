#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_combined_report_config_template.py — 三层目录 combined 报告配置。

复制到 combined 报告目录后改名 `_combined_report_config.py`。它和
`build_combined_report.py` 配套使用,把多个已有 slide-style 报告拼成一个
统一网页:

  1. super 层:例如 组织 / 队列 / 物种 / 模块
  2. 每个子报告自己的 CATEGORIES 顶部 nav
  3. 每个子报告自己的左侧 sidebar 页面目录

前提:每个子报告目录下已经有普通两层报告的 `_report_config.py` / REPORT.md /
report_figs/。
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "pages"
REPORT_FIGS = HERE / "report_figs"
CSS_NAME = "dark"

# ==================================================================
# 1) 子报告来源 ------------------------------------------------------
# key 会作为合并后页面前缀,只能用英文/数字/下划线。
# path 指向已有普通报告目录,该目录必须有 `_report_config.py`。
# ==================================================================
SOURCE_REPORTS = [
    {
        "key": "a",
        "label": "模块 A",
        "path": HERE.parent / "module_a",
        "landing_title": "Module A report",
        "summary": "atlas + subtype QC + aging screens",
    },
    {
        "key": "b",
        "label": "模块 B",
        "path": HERE.parent / "module_b",
        "landing_title": "Module B report",
        "summary": "mechanism + validation",
    },
]

# ==================================================================
# 2) 综合层 ----------------------------------------------------------
# cross key 对应根 landing(index.html)和跨模块内容页。
# COMBINED_REPORT 只写综合内容;子报告正文仍来自各自 REPORT.md。
# ==================================================================
CROSS_KEY = "cross"
CROSS_LABEL = "综合"
PROJECT_TITLE = "Combined report"
PROJECT_INTRO = (
    "Use the top switch to move across report modules; inside each module "
    "the original section navigation is preserved."
)
COMBINED_REPORT = HERE / "REPORT.md"

# (combined_slug, heading_prefix_in_COMBINED_REPORT, display_label)
CROSS_PAGES = [
    ("x_summary", "Summary", "Summary"),
]

# 可选:把某个子报告里的页面挪到综合层展示。
# 例: {("b", "module_a_vs_b"): CROSS_KEY}
ROUTE_OVERRIDES = {}

# 可选:综合层 sidebar 额外链接。通常只在 ROUTE_OVERRIDES 后补充展示名。
# 例: [("b_module_a_vs_b", "Module A vs B")]
CROSS_EXTRA_LINKS = []

# ==================================================================
# 3) 资源同步 --------------------------------------------------------
# 子报告 report_figs/ 会复制到 combined/report_figs/,文件名自动加
# "<source_key>__" 前缀避免重名。combined/figures/ 会加 "cross__" 前缀。
# ==================================================================
FIGURE_EXTENSIONS = (".png", ".pdf", ".svg", ".csv")
COMBINED_FIGURE_DIR = HERE / "figures"
