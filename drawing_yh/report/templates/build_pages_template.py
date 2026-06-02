#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
build_pages_template.py — REPORT.md → 多页 slide-style HTML(暗色主题)。

复制到项目目录改名 `build_pages.py`。配置(SLIDES / CATEGORIES)在
`_report_config.py`,本文件只放 HTML 渲染逻辑。

产物:
  /pages/index.html             — landing 总览
  /pages/<slide>.html × N       — 各 slide
  /pages/dark.css               — 自动复制
  /index.html                   — redirect → pages/index.html

特征:暗色 GitHub 风 / sticky 顶 4-tab(无前缀文字)/ 左 sidebar 按当前 cat 过滤
(landing 例外列全)/ 窄屏 sidebar 折 48 px 细栏 / 每页只一层标题。

约定细节见 `preferences.md` 的 "结果报告(slide-style HTML)" 节。
"""
import shutil
from drawing_yh import report
from _report_config import (
    HERE, SLIDES, CATEGORIES, cat_of, first_in_category, cat_label_of, load_chunks,
    FIGURE_SRC,
)

OUT = HERE / 'pages'
REPORT_TITLE = '项目结果报告'


def sync_report_figs():
    """从 FIGURE_SRC (数据层 figure/) 同步更新过的图到 HERE/report_figs/。
    解决 production scripts 输出到 D:/Projects 但 HTML 引用 report_figs/ 的脱节。
    FIGURE_SRC=None 跳过。"""
    if FIGURE_SRC is None or not FIGURE_SRC.exists():
        return
    report_figs = HERE / 'report_figs'
    if not report_figs.exists():
        return
    updated = 0
    for f in report_figs.iterdir():
        if f.suffix.lower() not in ('.png', '.pdf', '.svg'):
            continue
        src = FIGURE_SRC / f.name
        if src.exists() and src.stat().st_mtime > f.stat().st_mtime:
            shutil.copy2(src, f)
            updated += 1
    if updated:
        print(f'  synced {updated} updated figure(s) from {FIGURE_SRC} → report_figs/')
CSS = report.css_path('dark')        # 'dark' 或 'default'

TOP_NAV = [(ck, lbl, f'{first_in_category(ck)}.html') for ck, lbl in CATEGORIES]


def sidebar_items_for_slide(slug: str):
    """非 landing → 仅当前 cat 的 slide;landing → 列全部按 cat 分组。"""
    cat = cat_of(slug)
    if cat == 'overview':
        items = [('index', '总览页')]
        for ck, lbl in CATEGORIES[1:]:
            items.append(('section', lbl))
            for s, _, dl, c in SLIDES:
                if c == ck:
                    items.append((s, dl))
        return items
    items = [('section', cat_label_of(cat))]
    for s, _, dl, c in SLIDES:
        if c == cat:
            items.append((s, dl))
    return items


def write_slide(slug: str, body_md: str, title: str):
    cat = cat_of(slug)
    nav = report.build_nav(TOP_NAV, current=cat, label=None)
    sidebar = report.build_sidebar(sidebar_items_for_slide(slug), current=slug)
    header = sidebar + '\n' + nav + '\n<main class="report-main">'
    out = report.render(
        body_md,
        out_html=OUT / f'{slug}.html',
        css=CSS,
        title=title,
        toc=False,
        header_html=header,
        footer_html='</main>',
        extra_args=['--resource-path', str(HERE)],
    )
    print(f'  → {out.relative_to(HERE)}  ({out.stat().st_size:,} B)')


def main():
    OUT.mkdir(exist_ok=True)
    sync_report_figs()                  # ★ build 前先同步最新图
    chunks = load_chunks(fix_img_for_pages=True)

    print('rendering slides →', OUT)
    for slug, _, dlabel, _ in SLIDES:
        write_slide(slug, chunks[slug], dlabel)

    # landing(总览)
    landing_md = ['# 总览\n', f'共 **{len(SLIDES)} 张 slide**,按 {len(CATEGORIES)-1} 类组织。每张 ≤ 一页 PPT 量。\n']
    for ck, lbl in CATEGORIES[1:]:
        landing_md.append(f'\n## {lbl}\n')
        for slug, _, dlabel, c in SLIDES:
            if c == ck:
                landing_md.append(f'- [{dlabel}]({slug}.html)')
    write_slide('index', '\n'.join(landing_md), '总览')

    # 根目录 redirect
    redirect = (
        '<!doctype html><meta charset="utf-8">'
        '<meta http-equiv="refresh" content="0; url=pages/index.html">'
        f'<title>{REPORT_TITLE}</title>'
        '<p>跳转到 <a href="pages/index.html">pages/index.html</a> …</p>\n'
    )
    (HERE / 'index.html').write_text(redirect, encoding='utf-8')
    print(f'  → index.html  ({len(redirect):,} B)  [redirect → pages/index.html]')
    print('done.')


if __name__ == '__main__':
    main()
