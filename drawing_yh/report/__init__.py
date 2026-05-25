"""
drawing_yh.report
=================

Render Markdown research reports to GitHub-flavored standalone HTML using a
bundled CSS template.  Thin wrapper around `pandoc`.

Quick start
-----------
    from drawing_yh import report

    # 1) one-call render with bundled CSS
    report.render('REPORT.md')                       # → REPORT.html

    # 2) just need the CSS path (e.g. for a Makefile)
    print(report.css_path())                         # → /…/templates/default.css

    # 3) split a long report into linked pages
    md = open('REPORT.md', encoding='utf-8').read()
    chunks = report.split_by_h2(md, [
        ('00_summary',  ['摘要', '项目目标', 'Layer 1']),
        ('01_findings', ['4 个 finding', '附:gdcat', '数字总览']),
        ('02_appendix', ['基因级整合', '数据物种', '实验验证',
                         '数据互补', '附录 A', '附录 B', '致谢']),
    ])
    pages = [('00_summary', '摘要'),
             ('01_findings', 'Findings'),
             ('02_appendix', '附录')]
    for name, body in chunks.items():
        nav = report.build_nav(pages, current=name)
        report.render(body, out_html=f'pages/{name}.html',
                      title=name, header_html=nav)
"""
from .core import (
    css_path,
    list_templates,
    render,
    split_by_h2,
    build_nav,
    build_sidebar,
    slide_template_path,
    config_template_path,
    combined_config_template_path,
    combined_template_path,
    pptx_template_path,
    serve_nocache_template_path,
)

__all__ = [
    "css_path",
    "list_templates",
    "render",
    "split_by_h2",
    "build_nav",
    "build_sidebar",
    "slide_template_path",
    "config_template_path",
    "combined_config_template_path",
    "combined_template_path",
    "pptx_template_path",
    "serve_nocache_template_path",
]
