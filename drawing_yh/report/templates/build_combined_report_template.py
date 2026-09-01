#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
build_combined_report_template.py — 多个 slide-style 报告 → 三层目录统一网页。

复制到 combined 报告目录后改名 `build_combined_report.py`。配置在
`_combined_report_config.py`。

输出结构:
  pages/index.html          — combined landing
  pages/<key>_index.html    — 每个子报告 landing
  pages/<key>_<slug>.html   — 子报告页面
  pages/x_*.html            — 综合层页面
  report_figs/              — 自动合并图片,加前缀避免重名
  index.html                — redirect → pages/index.html
"""
from __future__ import annotations

import importlib.util
import re
import shutil
from pathlib import Path

from drawing_yh import report
from _combined_report_config import (
    COMBINED_FIGURE_DIR,
    COMBINED_REPORT,
    CROSS_EXTRA_LINKS,
    CROSS_KEY,
    CROSS_LABEL,
    CROSS_PAGES,
    CSS_NAME,
    FIGURE_EXTENSIONS,
    HERE,
    OUT,
    PROJECT_INTRO,
    PROJECT_TITLE,
    REPORT_FIGS,
    ROUTE_OVERRIDES,
    SOURCE_REPORTS,
)
try:
    from _combined_report_config import FAVICON_EMOJI as _FAVICON_EMOJI
except ImportError:
    _FAVICON_EMOJI = None


SUPER_CSS = """
<style>
.super-switch {
  display: flex; flex-wrap: wrap; gap: 6px 8px;
  justify-content: center; align-items: center;
  margin: 0 0 14px 0; padding: 8px 10px;
  background: var(--bg-elev); border: 1px solid var(--border); border-radius: 12px;
}
.super-switch .super-label {
  flex: 0 0 auto;
  color: var(--muted); font-size: .78em; font-weight: 600;
  letter-spacing: .08em; margin-right: 6px; text-transform: uppercase;
}
.super-switch a.super-btn {
  flex: 0 1 auto; min-width: 0; box-sizing: border-box;
  padding: 6px 22px; border-radius: 999px; color: var(--fg); font-weight: 600;
  font-size: .95em; white-space: nowrap; text-align: center;
  border: 1px solid var(--border); background: var(--bg-soft);
  transition: background .15s ease, color .15s ease, border-color .15s ease;
}
.super-switch a.super-btn:hover { background: var(--bg); color: var(--link-hover); text-decoration: none; }
.super-switch a.super-btn.current { background: var(--accent); color: #0d1117; border-color: var(--accent-strong); }
@media (max-width: 900px) {
  .super-switch { justify-content: flex-start; padding: 8px; }
  .super-switch .super-label { flex-basis: 100%; margin-right: 0; }
  .super-switch a.super-btn { padding: 5px 12px; font-size: .9em; }
}
</style>
"""

IMG_MD = re.compile(r"(!\[[^\]]*\]\()(?!(?:https?:|file:|#|/))([^)\s]+)", re.I)
IMG_HTML = re.compile(r'(<img\s+[^>]*src=")(?!(?:https?:|file:|#|/))([^"]+)"', re.I)


def load_py(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def src_key(row) -> str:
    return row["key"]


def src_label(row) -> str:
    return row["label"]


def src_path(row) -> Path:
    return Path(row["path"])


def src_landing_title(row) -> str:
    return row.get("landing_title") or f"{src_label(row)} report"


def src_summary(row) -> str:
    return row.get("summary", "")


def route_for(tag: str, slug: str) -> str:
    return ROUTE_OVERRIDES.get((tag, slug), tag)


def page_slug(tag: str, slug: str) -> str:
    return f"{tag}_{slug.replace('/', '_')}"


def cslug_for(cslug_of, tag: str, slug: str) -> str:
    return cslug_of[tag][slug]


def split_combined_report() -> dict[str, str]:
    if not COMBINED_REPORT.exists():
        return {}
    md = COMBINED_REPORT.read_text(encoding="utf-8")
    if not CROSS_PAGES:
        return {}
    return split_by_heading(md, [(heading, heading) for _cs, heading, _label in CROSS_PAGES])


_H = re.compile(r"^(#{2,3})\s+(.+?)\s*$", re.MULTILINE)


def split_by_heading(md_text: str, splits):
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
        if len(cands) == 1:
            return cands[0]
        if len(cands) > 1:
            raise KeyError(f"prefix {want!r} ambiguous: {cands}")
        raise KeyError(f"{want!r} not found in headings: {list(bounds)[:8]}")

    out = {}
    for slug, prefix in splits:
        t = resolve(prefix)
        s, e = bounds[t]
        out[slug] = md_text[s:e].rstrip() + "\n"
    return out


def build_model():
    cfgs = {}
    chunks = {}
    cslug_of = {}
    pages = []
    sources_by_key = {src_key(row): row for row in SOURCE_REPORTS}

    for row in SOURCE_REPORTS:
        tag = src_key(row)
        cfg = load_py(src_path(row) / "_report_config.py", f"cfg_{tag}")
        cfgs[tag] = cfg
        chunks[tag] = cfg.load_chunks(fix_img_for_pages=False)
        cslug_of[tag] = {"index": f"{tag}_index"}
        for slug, _prefix, label, cat in cfg.SLIDES:
            cslug = page_slug(tag, slug)
            cslug_of[tag][slug] = cslug
            pages.append({
                "cslug": cslug,
                "src": slug,
                "label": label,
                "tag": tag,
                "super": route_for(tag, slug),
                "cat": cat,
            })

    chunks[CROSS_KEY] = split_combined_report()
    cslug_of[CROSS_KEY] = {}
    for cslug, src, label in CROSS_PAGES:
        cslug_of[CROSS_KEY][src] = cslug
        pages.append({
            "cslug": cslug,
            "src": src,
            "label": label,
            "tag": CROSS_KEY,
            "super": CROSS_KEY,
            "cat": CROSS_KEY,
        })
    return sources_by_key, cfgs, chunks, cslug_of, pages


def sync_figs():
    REPORT_FIGS.mkdir(parents=True, exist_ok=True)
    copied = 0
    for row in SOURCE_REPORTS:
        tag = src_key(row)
        rf = src_path(row) / "report_figs"
        if not rf.exists():
            continue
        for f in rf.iterdir():
            if f.suffix.lower() not in FIGURE_EXTENSIONS:
                continue
            dst = REPORT_FIGS / f"{tag}__{f.name}"
            if not dst.exists() or f.stat().st_mtime > dst.stat().st_mtime:
                shutil.copy2(f, dst)
                copied += 1
    if COMBINED_FIGURE_DIR.exists():
        for f in COMBINED_FIGURE_DIR.rglob("*"):
            if f.suffix.lower() not in FIGURE_EXTENSIONS:
                continue
            dst = REPORT_FIGS / f"{CROSS_KEY}__{f.name}"
            if not dst.exists() or f.stat().st_mtime > dst.stat().st_mtime:
                shutil.copy2(f, dst)
                copied += 1
    print(f"  synced {copied} figure file(s) -> {REPORT_FIGS.name}/")


def extra_css_from_source(path: Path) -> str:
    bp = path / "build_pages.py"
    if not bp.exists():
        return ""
    text = bp.read_text(encoding="utf-8")
    m = re.search(r"<style>.*?</style>", text, re.S)
    return m.group(0) if m else ""


def super_items():
    return [(src_key(row), src_label(row), f"{src_key(row)}_index.html") for row in SOURCE_REPORTS] + [
        (CROSS_KEY, CROSS_LABEL, "index.html")
    ]


def super_switch(active: str) -> str:
    out = ['<div class="super-switch"><span class="super-label">切换</span>']
    for key, label, href in super_items():
        cls = " current" if key == active else ""
        out.append(f'<a class="super-btn{cls}" href="{href}">{label}</a>')
    out.append("</div>")
    return "".join(out)


def top_nav_for_source(cfg, tag: str, cslug_of, current_cat: str) -> str:
    items = []
    for ck, label in cfg.CATEGORIES:
        first = "index"
        if ck != "overview":
            for s, _prefix, _display, cat in cfg.SLIDES:
                if cat == ck and route_for(tag, s) == tag:
                    first = s
                    break
        items.append((ck, label, f"{cslug_for(cslug_of, tag, first)}.html"))
    return report.build_nav(items, current=current_cat, label=None)


def sidebar_for_source(cfg, tag: str, cslug_of, src_slug: str, current_cslug: str):
    cat = cfg.cat_of(src_slug)
    if cat == "overview":
        items = [(cslug_for(cslug_of, tag, "index"), "Overview", f"{cslug_for(cslug_of, tag, 'index')}.html")]
        for ck, label in cfg.CATEGORIES[1:]:
            rows = [(cslug_for(cslug_of, tag, s), dl, f"{cslug_for(cslug_of, tag, s)}.html")
                    for s, _p, dl, c in cfg.SLIDES
                    if c == ck and route_for(tag, s) == tag]
            if rows:
                items.append(("section", label))
                items += rows
        return report.build_sidebar(items, current=current_cslug)

    items = [("section", cfg.cat_label_of(cat))]
    items += [(cslug_for(cslug_of, tag, s), dl, f"{cslug_for(cslug_of, tag, s)}.html")
              for s, _p, dl, c in cfg.SLIDES
              if c == cat and route_for(tag, s) == route_for(tag, src_slug)]
    return report.build_sidebar(items, current=current_cslug)


def cross_nav(current: str) -> str:
    return report.build_nav([("overview", "Overview", "index.html"),
                             (CROSS_KEY, CROSS_LABEL, first_cross_href())],
                            current=current, label=None)


def first_cross_href() -> str:
    if CROSS_PAGES:
        return f"{CROSS_PAGES[0][0]}.html"
    if CROSS_EXTRA_LINKS:
        return f"{CROSS_EXTRA_LINKS[0][0]}.html"
    return "index.html"


def cross_sidebar(current_cslug: str):
    items = [("index", "Overview", "index.html")]
    rows = [(cslug, label, f"{cslug}.html") for cslug, _src, label in CROSS_PAGES]
    rows += [(cslug, label, f"{cslug}.html") for cslug, label in CROSS_EXTRA_LINKS]
    if rows:
        items.append(("section", CROSS_LABEL))
        items += rows
    return report.build_sidebar(items, current=current_cslug)


def rewrite_image_paths(md: str, tag: str) -> str:
    def normalized_name(path: str) -> str:
        if path.startswith("../"):
            path = path[3:]
        if path.startswith("report_figs/"):
            path = path[len("report_figs/"):]
        return Path(path).name

    def repl_md(m):
        return f"{m.group(1)}../report_figs/{tag}__{normalized_name(m.group(2))}"

    def repl_html(m):
        return f'{m.group(1)}../report_figs/{tag}__{normalized_name(m.group(2))}"'

    md = IMG_MD.sub(repl_md, md)
    md = IMG_HTML.sub(repl_html, md)
    return md


def rewrite_page_links(md: str, tag: str, pages) -> str:
    if tag == CROSS_KEY:
        return md
    pairs = sorted([(p["src"], p["cslug"]) for p in pages if p["tag"] == tag],
                   key=lambda x: -len(x[0]))
    for src, dst in pairs:
        pat = re.compile(r'(\]\(|href=")(?:\./)?' + re.escape(src) + r"\.html")
        md = pat.sub(lambda m, dst=dst: f"{m.group(1)}{dst}.html", md)
    return md


def render_slide(cslug, body_md, title, super_key, tag, cat, sources, cfgs, cslug_of, pages, css, extra_css):
    if super_key == CROSS_KEY:
        nav = cross_nav("overview" if cslug == "index" else CROSS_KEY)
        sidebar = cross_sidebar(cslug)
    else:
        cfg = cfgs[tag]
        src = "index" if cslug == cslug_for(cslug_of, tag, "index") else next(
            p["src"] for p in pages if p["cslug"] == cslug)
        nav = top_nav_for_source(cfg, tag, cslug_of, cat)
        sidebar = sidebar_for_source(cfg, tag, cslug_of, src, cslug)
    header = extra_css + "\n" + sidebar + "\n" + super_switch(super_key) + "\n" + nav + '\n<main class="report-main">'
    # 含 **PPT 页标题** 字段的 chunk → 结构化卡片(汇报草稿预览将进 PPT 的内容);其余原样
    body_md = report.briefing_card(body_md)
    body = rewrite_page_links(rewrite_image_paths(body_md, tag), tag, pages)
    out = report.render(body, out_html=OUT / f"{cslug}.html", css=css, title=title,
                        toc=False, favicon_emoji=_FAVICON_EMOJI, header_html=header,
                        footer_html="</main>", extra_args=["--resource-path", str(HERE)])
    print(f"  -> {out.relative_to(HERE)} ({out.stat().st_size:,} B)")


def source_landing(row, cfg, tag, cslug_of):
    pages = [(s, dl, c) for s, _p, dl, c in cfg.SLIDES if route_for(tag, s) == tag]
    lines = [f"# {src_landing_title(row)}\n"]
    if src_summary(row):
        lines.append(src_summary(row) + "\n")
    for ck, label in cfg.CATEGORIES[1:]:
        rows = [f"- [{dl}]({cslug_for(cslug_of, tag, s)}.html)" for s, dl, c in pages if c == ck]
        if rows:
            lines.append(f"\n## {label}\n")
            lines += rows
    return "\n".join(lines)


def combined_landing(sources, pages):
    lines = [f"# {PROJECT_TITLE}\n", PROJECT_INTRO + "\n", "\n## Modules\n"]
    for row in SOURCE_REPORTS:
        tag = src_key(row)
        n = sum(1 for p in pages if p["super"] == tag)
        desc = f" — {src_summary(row)}" if src_summary(row) else ""
        lines.append(f"- **{src_label(row)}** — {n} pages{desc}. [Open]({tag}_index.html)")
    if CROSS_PAGES or CROSS_EXTRA_LINKS:
        lines.append(f"\n## {CROSS_LABEL}\n")
        for cslug, _src, label in CROSS_PAGES:
            lines.append(f"- [{label}]({cslug}.html)")
        for cslug, label in CROSS_EXTRA_LINKS:
            lines.append(f"- [{label}]({cslug}.html)")
    return "\n".join(lines)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sync_figs()
    sources, cfgs, chunks, cslug_of, pages = build_model()
    css = report.css_path(CSS_NAME)
    source_css = []
    for row in SOURCE_REPORTS:
        ecss = extra_css_from_source(src_path(row))
        if ecss and ecss not in source_css:
            source_css.append(ecss)
    extra_css = SUPER_CSS + "\n".join(source_css)

    print(f"rendering {len(pages) + len(SOURCE_REPORTS) + 1} pages -> {OUT}")
    for p in pages:
        body = chunks[p["tag"]][p["src"]]
        render_slide(p["cslug"], body, p["label"], p["super"], p["tag"], p["cat"],
                     sources, cfgs, cslug_of, pages, css, extra_css)

    for row in SOURCE_REPORTS:
        tag = src_key(row)
        _entry = row.get("entry_slug", "index")
        _idx = cslug_for(cslug_of, tag, "index")
        if _entry != "index":
            # 配了 entry_slug(入口是某概览 slide)→ <tag>_index 重定向到入口页,
            # 避免再生成一个重复的目录 landing 而成孤儿(2026-06-25)
            _tgt = f"{cslug_for(cslug_of, tag, _entry)}.html"
            (OUT / f"{_idx}.html").write_text(
                '<!doctype html><meta charset="utf-8">'
                f'<meta http-equiv="refresh" content="0; url={_tgt}">'
                f'<title>{PROJECT_TITLE}</title>'
                f'<p>Redirect to <a href="{_tgt}">overview</a>.</p>\n',
                encoding="utf-8")
        else:
            render_slide(_idx,
                         source_landing(row, cfgs[tag], tag, cslug_of),
                         src_landing_title(row), tag, tag, "overview",
                         sources, cfgs, cslug_of, pages, css, extra_css)

    render_slide("index", combined_landing(sources, pages), PROJECT_TITLE,
                 CROSS_KEY, CROSS_KEY, "overview",
                 sources, cfgs, cslug_of, pages, css, extra_css)

    redirect = (
        '<!doctype html><meta charset="utf-8">'
        '<meta http-equiv="refresh" content="0; url=pages/index.html">'
        f'<title>{PROJECT_TITLE}</title><p>Redirect to <a href="pages/index.html">report</a>.</p>\n'
    )
    (HERE / "index.html").write_text(redirect, encoding="utf-8")
    report.check_orphan_pages(OUT)
    print("done.")


if __name__ == "__main__":
    main()
