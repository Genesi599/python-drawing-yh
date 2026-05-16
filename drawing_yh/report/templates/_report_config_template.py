#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_report_config_template.py — REPORT.md 共享配置 + markdown 切分。

复制到项目目录改名 `_report_config.py`,改这三处即可:
  1. SLIDES        —— 你的 H2/H3 标题清单 + 显示 label + 分类
  2. CATEGORIES    —— 顶部 4-tab 分类
  3. (`SRC` / 图片目录名 默认 'report_figs',按需调整)

`build_pages.py`(出 HTML)和 `md_to_pptx.py`(出 PPTX)都从这里 import,
**纯 import 不会触发任何渲染**,放心用。
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC  = HERE / 'REPORT.md'
# 数据层 figure 目录(production scripts 输出到这里);build_pages.py 跑前自动 sync
# 同名图过来到 HERE/report_figs/。None 禁用 sync(图手动管理)。
FIGURE_SRC = None    # e.g. Path(r'D:\Projects\<project>\<sub>\figure')

# ==================================================================
# SLIDES:每张幻灯 1 行 ----------------------------------------------
# (slug, heading_prefix, display_label, category_key)
# ==================================================================
SLIDES = [
    ('intro_summary',    '摘要',       '摘要',      'intro'),
    ('intro_objective',  '项目目标',   '项目目标',  'intro'),
    # ... 加你自己的
    ('find_A',           'Finding A',  'Finding A', 'findings'),
    ('find_B',           'Finding B',  'Finding B', 'findings'),
    # ... 加你自己的
    ('appx_thanks',      '致谢',       '致谢',      'appendix'),
]

# ==================================================================
# CATEGORIES:顶部 tab(overview 必须排第一,代表 landing)----------
# ==================================================================
CATEGORIES = [
    ('overview', '总览'),
    ('intro',    '摘要 · Layer 1'),
    ('findings', '4 个 Finding'),
    ('appendix', '附录'),
]

# ==================================================================
# 以下框架代码不用改 ------------------------------------------------
# ==================================================================
def cat_of(slug: str) -> str:
    if slug == 'index':
        return 'overview'
    for s, _, _, ck in SLIDES:
        if s == slug:
            return ck
    return 'overview'


def first_in_category(cat_key: str) -> str:
    if cat_key == 'overview':
        return 'index'
    for s, _, _, ck in SLIDES:
        if ck == cat_key:
            return s
    return 'index'


def cat_label_of(ck: str) -> str:
    for k, lbl in CATEGORIES:
        if k == ck:
            return lbl
    return ck


# Splitter:H2 + H3 都当边界,模糊前缀匹配
_H = re.compile(r'^(#{2,3})\s+(.+?)\s*$', re.MULTILINE)


def split_by_heading(md_text: str, splits):
    """splits = list of (slug, heading_prefix). 返回 dict[slug → md chunk]。
    每页只保留自己 chunk 的标题,不 prepend 全文 H1。"""
    matches = list(_H.finditer(md_text))
    bounds = {}
    for i, m in enumerate(matches):
        title = m.group(2)
        start = m.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(md_text)
        bounds[title] = (start, end)

    def resolve(want: str) -> str:
        if want in bounds:
            return want
        cands = [t for t in bounds if t.startswith(want)]
        if len(cands) == 1: return cands[0]
        if len(cands) > 1: raise KeyError(f"prefix {want!r} ambiguous: {cands}")
        raise KeyError(f"{want!r} not found in headings: {list(bounds)[:5]}…")

    out = {}
    for slug, prefix in splits:
        t = resolve(prefix)
        s, e = bounds[t]
        out[slug] = md_text[s:e].rstrip() + '\n'
    return out


_IMG_FIX_PAGES = re.compile(r'(<img\s+[^>]*src=")report_figs/', re.IGNORECASE)


def load_chunks(*, fix_img_for_pages: bool = True) -> dict[str, str]:
    """读 REPORT.md → split_by_heading → 修图相对路径(若产物落在 pages/ 子目录)。
    `fix_img_for_pages=True` 把 'report_figs/' → '../report_figs/';
    `False` 保持原路径(适合 PPTX / 跟 REPORT.md 同根的产物)。
    若你的项目用了别的图片目录名,改本函数里的 _IMG_FIX_PAGES 正则。"""
    md = SRC.read_text(encoding='utf-8')
    chunks = split_by_heading(md, [(slug, prefix) for slug, prefix, _, _ in SLIDES])
    if fix_img_for_pages:
        for k in list(chunks):
            chunks[k] = _IMG_FIX_PAGES.sub(r'\1../report_figs/', chunks[k])
    return chunks
