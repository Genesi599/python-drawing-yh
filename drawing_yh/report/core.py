"""
drawing_yh.report.core
======================

Render Markdown research reports to standalone HTML using a bundled GitHub-flavored
CSS template.  Keeps `pandoc` as the engine (already on most workstations) and
ships only the styling + a thin Python wrapper.

Public API (re-exported from `drawing_yh.report`):

    css_path(name='default')      → Path to bundled CSS
    render(md, out_html=None, ...) → run pandoc, return Path
    split_by_h2(md_text, splits)   → dict[name → markdown chunk]
    build_nav(pages, current)      → HTML <nav> string

Typical project use::

    from drawing_yh import report

    # one-off single-page render
    report.render('REPORT.md')                     # writes index.html beside REPORT.md

    # multi-page split + render
    md = open('REPORT.md', encoding='utf-8').read()
    chunks = report.split_by_h2(md, [
        ('00_summary',  ['摘要', '项目目标', 'Layer 1']),
        ('01_findings', ['4 个 finding', '附:gdcat', '数字总览']),
        ('02_appendix', ['基因级整合', '数据物种', '实验验证',
                         '数据互补', '附录 A', '附录 B', '致谢']),
    ])
    nav_pages = [(n, n.split('_', 1)[1].title()) for n in chunks]
    for name, body_md in chunks.items():
        nav_html = report.build_nav(nav_pages, current=name)
        report.render(body_md, out_html=f'pages/{name}.html',
                      title=name, header_html=nav_html)
"""
from __future__ import annotations

import re
import shlex
import shutil
import subprocess
import tempfile
from importlib import resources
from pathlib import Path
from typing import Iterable, Mapping, Sequence


# ------------------------------------------------------------------
# CSS template lookup
# ------------------------------------------------------------------
def css_path(name: str = 'default') -> Path:
    """Return absolute filesystem path to a bundled CSS template.

    Parameters
    ----------
    name : str
        Template stem.  Currently only ``'default'`` ships.

    Returns
    -------
    pathlib.Path
        Absolute path under the installed package's ``templates/`` directory.
    """
    # importlib.resources.files works for both regular and editable installs
    p = resources.files(__package__).joinpath('templates', f'{name}.css')
    # On editable installs `p` is already a real Path; on zipped installs we
    # materialize a temp copy so pandoc -c can find it on disk.
    if hasattr(p, 'read_bytes') and not Path(str(p)).exists():
        tmp = Path(tempfile.gettempdir()) / f'drawing_yh_{name}.css'
        tmp.write_bytes(p.read_bytes())
        return tmp
    return Path(str(p))


def list_templates() -> list[str]:
    """List available CSS template names."""
    tdir = Path(str(resources.files(__package__).joinpath('templates')))
    return sorted(p.stem for p in tdir.glob('*.css'))


def slide_template_path() -> Path:
    """Return path to the bundled `build_pages_template.py` (HTML 渲染部分)。

    跟 `config_template_path()`(配置)配套使用。Copy 两个到项目目录(分别改名为
    `build_pages.py` / `_report_config.py`),改 `_report_config.py` 里的
    `SLIDES` / `CATEGORIES`,跑 `python build_pages.py` → 得到 sticky 顶 nav +
    filtered 左 sidebar + dark CSS 的多页 slide-style HTML 报告。
    """
    return Path(str(resources.files(__package__).joinpath('templates', 'build_pages_template.py')))


def config_template_path() -> Path:
    """Return path to the bundled `_report_config_template.py`(共享配置)。

    跟 `slide_template_path()` 和 `pptx_template_path()` 都配套使用。
    内容:`SLIDES` / `CATEGORIES` / `split_by_heading()` / `load_chunks()`,
    无副作用,被 build_pages.py / md_to_pptx.py import 不会跑渲染。
    """
    return Path(str(resources.files(__package__).joinpath('templates', '_report_config_template.py')))


def pptx_template_path() -> Path:
    """Return path to the bundled `md_to_pptx_template.py`(PPT 渲染)。

    跟 `config_template_path()` 配套。Copy 到项目目录(改名 `md_to_pptx.py`),
    改 4 处:`PROJECT_LABEL` / `TEMPLATE`(母版 .pptx)/ `OUT` / 标题页文案,
    跑 `python md_to_pptx.py` → 得到继承母版 theme 的 native PPTX。

    布局:每张 slide 顶部有项目标签 + 标题 + 橙色分隔线。有图 → 左大图 + 右
    ➤ bullets;无图 → 全宽。表格列宽按内容长短分配,字号自动 14→11pt。
    """
    return Path(str(resources.files(__package__).joinpath('templates', 'md_to_pptx_template.py')))


def serve_nocache_template_path() -> Path:
    """Return path to the bundled `serve_nocache_template.py`(本地预览 server)。

    `python -m http.server` 套壳,每个响应加 `Cache-Control: no-store` headers,
    避免重生成 HTML / 重做图后浏览器拿缓存看到旧版(尤其 sidebar 跳转时)。

    用法:`cp $(... serve_nocache_template_path())  serve_nocache.py` →
        `.claude/launch.json` 配 `python serve_nocache.py 8766`。
    """
    return Path(str(resources.files(__package__).joinpath('templates', 'serve_nocache_template.py')))


def watch_build_template_path() -> Path:
    """Return path to the bundled `watch_build_template.py`(文件监听自动 rebuild)。

    监听 REPORT.md / _report_config.py / report_figs/ 变动,自动重跑 build_pages.py,
    实现改 md → 浏览器自动看到最新 HTML。配合 serve_nocache.py 使用。

    用法:`cp $(... watch_build_template_path())  watch_build.py` →
        `python watch_build.py`(监听当前目录)。依赖 `pip install watchdog`。
    """
    return Path(str(resources.files(__package__).joinpath('templates', 'watch_build_template.py')))


def combined_config_template_path() -> Path:
    """Return path to `_combined_report_config_template.py`(三层报告配置)。

    跟 `combined_template_path()` 配套。用于把多个已有 slide-style 报告拼成一个
    统一网页:顶层 super switch(如组织/物种/模块) + 子报告 category nav +
    子报告 sidebar。Copy 到 combined 目录后改名 `_combined_report_config.py`。
    """
    return Path(str(resources.files(__package__).joinpath('templates', '_combined_report_config_template.py')))


def combined_template_path() -> Path:
    """Return path to `build_combined_report_template.py`(三层报告渲染入口)。

    Copy 到 combined 目录后改名 `build_combined_report.py`,配置见
    `_combined_report_config.py`。它会非破坏式读取各子报告的 `_report_config.py`
    / `REPORT.md`,同步 report_figs,生成带三层目录的 `pages/*.html`。
    """
    return Path(str(resources.files(__package__).joinpath('templates', 'build_combined_report_template.py')))


# ------------------------------------------------------------------
# Publish helper
# ------------------------------------------------------------------
_PROJECT_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def _validate_project_name(project_name: str) -> str:
    project = str(project_name).strip().strip("/\\")
    if not _PROJECT_NAME_RE.fullmatch(project):
        raise ValueError(
            "project_name must start with an ASCII letter/number and contain "
            "only ASCII letters, numbers, '.', '_' or '-'."
        )
    if project in {".", ".."} or ".." in project.split("."):
        raise ValueError("project_name must not contain parent-directory traversal.")
    return project


def publish(
    project_name: str,
    host: str = "todo",
    *,
    source_dir: str | Path = ".",
    remote_root: str = "/var/www/reports",
    include_index: bool = True,
    fix_remote_permissions: bool = True,
    dry_run: bool = False,
) -> str:
    """Publish a slide-style HTML report directory with ssh/scp.

    The expected local layout is the standard ``drawing_yh.report`` output:
    ``pages/`` plus optional ``report_figs/`` and root ``index.html``.

    Parameters
    ----------
    project_name : str
        Remote project folder name under ``remote_root``. Only conservative
        ASCII names are allowed to avoid accidental path traversal.
    host : str, default "todo"
        SSH host alias from the user's ``~/.ssh/config``.
    source_dir : str | Path, default "."
        Directory containing ``pages/`` and ``report_figs/``.
    remote_root : str, default "/var/www/reports"
        Remote reports root.
    include_index : bool, default True
        Copy root ``index.html`` when present, so ``/<project>/`` redirects to
        ``pages/index.html`` in standard reports.
    fix_remote_permissions : bool, default True
        Make uploaded directories/files readable by nginx-style static servers.
    dry_run : bool, default False
        Print the ssh/scp commands without running them.

    Returns
    -------
    str
        Remote directory path, e.g. ``/var/www/reports/cellchat_x_BMIF``.
    """
    project = _validate_project_name(project_name)
    source = Path(source_dir).resolve()
    remote_dir = f"{remote_root.rstrip('/')}/{project}"

    pages = source / "pages"
    if not pages.is_dir():
        raise FileNotFoundError(f"Required report directory not found: {pages}")

    path_names = ["pages"]
    report_figs = source / "report_figs"
    if report_figs.is_dir():
        path_names.append("report_figs")

    index_html = source / "index.html"
    if include_index and index_html.is_file():
        path_names.append("index.html")

    mkdir_cmd = ["ssh", host, f"mkdir -p {shlex.quote(remote_dir)}"]
    scp_cmd = [
        "scp",
        "-r",
        *path_names,
        f"{host}:{remote_dir}/",
    ]
    chmod_cmd = [
        "ssh",
        host,
        (
            f"find {shlex.quote(remote_dir)} -type d -exec chmod 755 {{}} + && "
            f"find {shlex.quote(remote_dir)} -type f -exec chmod 644 {{}} +"
        ),
    ]

    if dry_run:
        print(f"cd {shlex.quote(str(source))}")
        print(" ".join(mkdir_cmd))
        print(" ".join(shlex.quote(part) for part in scp_cmd))
        if fix_remote_permissions:
            print(" ".join(chmod_cmd))
        return remote_dir

    subprocess.run(mkdir_cmd, check=True)
    subprocess.run(scp_cmd, check=True, cwd=source)
    if fix_remote_permissions:
        subprocess.run(chmod_cmd, check=True)
    return remote_dir


# ------------------------------------------------------------------
# Image toolbar — inline JS injected by render(img_toolbar=True)
# ------------------------------------------------------------------
_IMG_TOOLBAR_SCRIPT = """\
<script>
document.addEventListener('DOMContentLoaded', function () {
  document.querySelectorAll('img').forEach(function (img) {
    if (img.closest('.img-toolbar-wrap')) return;
    var wrap = document.createElement('div');
    wrap.className = 'img-toolbar-wrap';
    img.parentNode.insertBefore(wrap, img);
    wrap.appendChild(img);
    var tb = document.createElement('div');
    tb.className = 'img-toolbar';
    var cpBtn = document.createElement('button');
    cpBtn.textContent = '复制';
    cpBtn.title = 'Copy image to clipboard';
    cpBtn.onclick = function () {
      var c = document.createElement('canvas');
      c.width = img.naturalWidth; c.height = img.naturalHeight;
      c.getContext('2d').drawImage(img, 0, 0);
      c.toBlob(function (blob) {
        navigator.clipboard.write([new ClipboardItem({'image/png': blob})]).then(function () {
          cpBtn.textContent = '✓';
          setTimeout(function () { cpBtn.textContent = '复制'; }, 1200);
        }).catch(function () {
          cpBtn.textContent = '失败';
          setTimeout(function () { cpBtn.textContent = '复制'; }, 1200);
        });
      }, 'image/png');
    };
    var dlBtn = document.createElement('button');
    dlBtn.textContent = '下载';
    dlBtn.title = 'Download image';
    dlBtn.onclick = function () {
      var a = document.createElement('a');
      a.href = img.src;
      a.download = img.src.split('/').pop().split('?')[0] || 'image.png';
      document.body.appendChild(a); a.click(); document.body.removeChild(a);
    };
    tb.appendChild(cpBtn); tb.appendChild(dlBtn);
    wrap.appendChild(tb);
  });
});
</script>
"""

# ------------------------------------------------------------------
# Image lightbox — inline JS injected by render(img_lightbox=True)
# ------------------------------------------------------------------
_IMG_LIGHTBOX_SCRIPT = """\
<script>
document.addEventListener('click', function (e) {
  if (e.target.tagName !== 'IMG' || e.target.closest('.img-toolbar')) return;
  var m = document.createElement('div');
  m.className = 'img-lightbox';
  var img = document.createElement('img');
  img.src = e.target.src;
  img.alt = e.target.alt || '';
  m.appendChild(img);
  m.addEventListener('click', function () { m.remove(); });
  document.body.appendChild(m);
});
document.addEventListener('keydown', function (e) {
  if (e.key === 'Escape') {
    document.querySelectorAll('.img-lightbox').forEach(function (m) { m.remove(); });
  }
});
</script>
"""

# ------------------------------------------------------------------
# pandoc wrapper
# ------------------------------------------------------------------
def _have_pandoc() -> bool:
    return shutil.which('pandoc') is not None


_EMOJI_MAP: list[tuple[re.Pattern, str]] = [
    (re.compile(r'retina|视网膜|眼', re.I), '👁'),
    (re.compile(r'bone.?marrow|骨髓|bmif', re.I), '🦴'),
    (re.compile(r'brain|脑|meninges|脑膜|frontal|额叶', re.I), '🧠'),
    (re.compile(r'b.?cell|b细胞|b_cell', re.I), '🧫'),
    (re.compile(r'thymus|胸腺', re.I), '🫀'),
    (re.compile(r'neutrophil|中性粒', re.I), '🩸'),
    (re.compile(r'aging|衰老', re.I), '⏳'),
]
_EMOJI_DEFAULT = '🔬'


def _emoji_for_title(title: str) -> str:
    for pattern, emoji in _EMOJI_MAP:
        if pattern.search(title):
            return emoji
    return _EMOJI_DEFAULT


def render(
    md,
    out_html: str | Path | None = None,
    css: str | Path | None = None,
    *,
    toc: bool = True,
    toc_depth: int = 3,
    standalone: bool = True,
    title: str | None = None,
    favicon_emoji: str | None = None,
    header_html: str | None = None,
    footer_html: str | None = None,
    copy_css: bool = True,
    img_toolbar: bool = True,
    img_lightbox: bool = True,
    extra_args: Sequence[str] = (),
) -> Path:
    """Render Markdown → standalone HTML via pandoc.

    Parameters
    ----------
    md : str | Path
        Either a path to a `.md` file, or raw Markdown text.  Detection: if it
        looks like a path that exists on disk, treated as a file; otherwise
        as raw text.
    out_html : str | Path, optional
        Output path.  Defaults to ``<md_stem>.html`` next to the source when
        `md` is a file, else ``./report.html``.
    css : str | Path, optional
        Path to a CSS file.  Defaults to the bundled ``default.css``.
    toc : bool
        Insert a table of contents.
    toc_depth : int
        Heading level to descend to in the TOC.
    standalone : bool
        Generate a standalone HTML document (header + body + closing tags).
    title : str, optional
        Document `<title>` and ``<h1>`` (pandoc ``--metadata title=…``).
    header_html : str, optional
        Raw HTML injected into the document via ``--include-before-body``.
        Useful for cross-page nav bars (see :func:`build_nav`).
    footer_html : str, optional
        Raw HTML injected via ``--include-after-body``.  Pair with
        ``header_html`` to wrap the body in a custom container — e.g. open
        ``<main>`` in ``header_html`` and close it in ``footer_html`` so a
        sidebar (see :func:`build_sidebar`) can layout next to the content.
    copy_css : bool, default True
        When ``True``, copy the CSS file next to the output HTML and reference
        it by basename.  Required for the output to render correctly when
        opened over ``http://`` (browsers refuse cross-origin ``file://`` CSS).
        Disable only if you want to share a CSS file across many outputs and
        manage the link path yourself.
    favicon_emoji : str, optional
        Emoji character used as the browser-tab favicon (data-URI SVG).
        When omitted, ``_emoji_for_title(title)`` auto-selects from keywords in
        *title*.  Pass an empty string ``""`` to suppress the favicon entirely.
    img_toolbar : bool, default True
        When ``True``, inject an inline ``<script>`` that adds hover copy /
        download buttons to every ``<img>`` on the page.  Requires the
        ``.img-toolbar-wrap`` / ``.img-toolbar`` CSS classes in the stylesheet.
    img_lightbox : bool, default True
        When ``True``, inject an inline ``<script>`` that opens a fullscreen
        modal overlay when any ``<img>`` is clicked.  Click overlay or press
        ESC to dismiss.  Toolbar buttons are unaffected.
    extra_args : Sequence[str], optional
        Additional pandoc CLI args appended verbatim.

    Returns
    -------
    pathlib.Path
        Path to the generated HTML.
    """
    if not _have_pandoc():
        raise RuntimeError(
            "pandoc not found on PATH. Install from https://pandoc.org/ "
            "or `winget install JohnMacFarlane.Pandoc`."
        )

    # 1. resolve input
    md = Path(md) if isinstance(md, (str, Path)) and Path(str(md)).exists() else md
    src_path: Path
    cleanup_src = False
    if isinstance(md, Path):
        src_path = md
    else:
        # raw text → temp file
        tmp = tempfile.NamedTemporaryFile(
            'w', encoding='utf-8', suffix='.md', delete=False
        )
        tmp.write(md if isinstance(md, str) else str(md))
        tmp.close()
        src_path = Path(tmp.name)
        cleanup_src = True

    # 2. resolve output
    if out_html is None:
        out_html = (
            src_path.with_suffix('.html')
            if not cleanup_src else Path.cwd() / 'report.html'
        )
    out_html = Path(out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)

    # 3. css — by default copy next to output so the link is relative
    css_src = Path(css) if css else css_path('default')
    if copy_css:
        css_dst = out_html.parent / css_src.name
        if not css_dst.exists() or css_dst.read_bytes() != css_src.read_bytes():
            css_dst.write_bytes(css_src.read_bytes())
        css_p = Path(css_dst.name)            # relative ref, same dir as html
    else:
        css_p = css_src

    # 4. optional pre-body / post-body includes
    inc_paths: list[Path] = []
    cleanup_includes: list[Path] = []
    def _spill(html: str | None) -> Path | None:
        if not html:
            return None
        t = tempfile.NamedTemporaryFile(
            'w', encoding='utf-8', suffix='.html', delete=False
        )
        t.write(html); t.close()
        p = Path(t.name)
        cleanup_includes.append(p)
        return p
    inc_before = _spill(header_html)
    inc_after  = _spill(footer_html)

    # 4b. image toolbar script (inline <script> injected after body)
    inc_toolbar = None
    if img_toolbar:
        inc_toolbar = _spill(_IMG_TOOLBAR_SCRIPT)

    # 4b2. image lightbox script (inline <script> injected after body)
    inc_lightbox = None
    if img_lightbox:
        inc_lightbox = _spill(_IMG_LIGHTBOX_SCRIPT)

    # 4c. emoji favicon injected into <head> via --include-in-header
    inc_favicon = None
    _emoji = favicon_emoji if favicon_emoji is not None else (
        _emoji_for_title(title) if title else _EMOJI_DEFAULT
    )
    if _emoji:
        _svg = (
            f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'>"
            f"<text y='.9em' font-size='90'>{_emoji}</text></svg>"
        )
        _favicon_html = f'<link rel="icon" href="data:image/svg+xml,{_svg}">\n'
        inc_favicon = _spill(_favicon_html)

    # 5. build command
    cmd: list[str] = ['pandoc', str(src_path), '-o', str(out_html), '-c', str(css_p)]
    if standalone:
        cmd.append('-s')
    if toc:
        cmd += ['--toc', f'--toc-depth={toc_depth}']
    if title:
        # `pagetitle` 只设 <head><title>(浏览器 tab),不渲染 <h1 class="title">。
        # 如果调用方真的想让 pandoc 把 title 渲染到正文顶部,可以用 extra_args
        # 自己塞 ['--metadata', 'title=…']。
        cmd += ['--metadata', f'pagetitle={title}']
    if inc_favicon:
        cmd += ['--include-in-header', str(inc_favicon)]
    if inc_before:
        cmd += ['--include-before-body', str(inc_before)]
    if inc_after:
        cmd += ['--include-after-body', str(inc_after)]
    if inc_toolbar:
        cmd += ['--include-after-body', str(inc_toolbar)]
    if inc_lightbox:
        cmd += ['--include-after-body', str(inc_lightbox)]
    cmd += list(extra_args)

    try:
        subprocess.run(cmd, check=True)
    finally:
        if cleanup_src:
            src_path.unlink(missing_ok=True)
        for p in cleanup_includes:
            p.unlink(missing_ok=True)

    return out_html


# ------------------------------------------------------------------
# Markdown splitter (by H2 sections)
# ------------------------------------------------------------------
_H2 = re.compile(r'^##\s+(.+?)\s*$', re.MULTILINE)


def _h2_index(md_text: str) -> list[tuple[int, int, str]]:
    """Return list of (line_idx, char_offset, heading_text) for every `## ` line."""
    out = []
    for m in _H2.finditer(md_text):
        out.append((md_text[:m.start()].count('\n'), m.start(), m.group(1)))
    return out


def split_by_h2(
    md_text: str,
    splits: Sequence[tuple[str, Sequence[str]]],
    *,
    h1_passthrough: bool = True,
    fuzzy: bool = True,
) -> dict[str, str]:
    """Cut a Markdown document into named chunks along ``##`` boundaries.

    Parameters
    ----------
    md_text : str
        Full Markdown source.
    splits : sequence of (page_name, list_of_h2_titles)
        Each entry produces one chunk containing those H2 sections (and their
        descendants until the next H2).  Order of `splits` defines page order.
    h1_passthrough : bool
        If True, the document's first ``# `` line (if any) is prepended to
        every chunk so each page has a top-level title.
    fuzzy : bool
        Match a requested heading if any of the headings in the document
        *startswith* it (case-sensitive).  Lets you write ``'Layer 1'`` and
        match ``'## Layer 1 · 蛋白 LR 通讯网络全景(...)'``.

    Returns
    -------
    dict[str, str]
        ``{page_name: markdown_chunk}`` in input order.

    Raises
    ------
    KeyError
        If a requested heading cannot be located.
    """
    h2s = _h2_index(md_text)  # (line, char_offset, title)
    if not h2s:
        raise ValueError("No `## ` headings found in document.")

    # H1 prefix (optional)
    h1_prefix = ''
    if h1_passthrough:
        m = re.search(r'^#\s+.+?$', md_text, flags=re.MULTILINE)
        if m:
            h1_prefix = m.group(0) + '\n\n'

    # Map heading text → (start_offset, end_offset)
    bounds: dict[str, tuple[int, int]] = {}
    for i, (_, start, title) in enumerate(h2s):
        end = h2s[i + 1][1] if i + 1 < len(h2s) else len(md_text)
        bounds[title] = (start, end)

    def resolve(want: str) -> str:
        if want in bounds:
            return want
        if fuzzy:
            cands = [t for t in bounds if t.startswith(want)]
            if len(cands) == 1:
                return cands[0]
            if len(cands) > 1:
                raise KeyError(
                    f"Heading prefix {want!r} is ambiguous, matches: {cands}"
                )
        raise KeyError(
            f"Heading {want!r} not found. Available: {list(bounds.keys())}"
        )

    out: dict[str, str] = {}
    for name, wants in splits:
        pieces = []
        for w in wants:
            t = resolve(w)
            s, e = bounds[t]
            pieces.append(md_text[s:e].rstrip() + '\n')
        out[name] = h1_prefix + '\n'.join(pieces)
    return out


# ------------------------------------------------------------------
# Cross-page nav bar
# ------------------------------------------------------------------
def build_nav(
    pages: Iterable[tuple],
    current: str | None = None,
    *,
    label: str = '页面',
    suffix: str = '.html',
) -> str:
    """Build a small ``<nav class="report-nav">`` HTML block linking sibling pages.

    Parameters
    ----------
    pages : iterable of tuples
        Each entry is ``(slug, display_label)`` or
        ``(slug, display_label, href)``.  When a 3-tuple is given, ``href`` is
        used verbatim (no suffix munging) — useful when nav links straddle
        directory levels (e.g. some links go to siblings, others to ``../``).
    current : str, optional
        The current page slug.  Marked with ``class="current"``.
    label : str
        Leading text in the nav bar.
    suffix : str
        File suffix appended to each 2-tuple slug to form the href.

    Returns
    -------
    str
        HTML string.  Pair with :func:`render` via ``header_html=…``.
    """
    items = ['<nav class="report-nav">']
    if label:
        items.append(f'  <span class="label">{label}</span>')
    for entry in pages:
        if len(entry) == 3:
            slug, lbl, href = entry
        else:
            slug, lbl = entry
            href = f'{slug}{suffix}'
        cls = ' class="current"' if slug == current else ''
        items.append(f'  <a href="{href}"{cls}>{lbl}</a>')
    items.append('</nav>')
    return '\n'.join(items)


# ------------------------------------------------------------------
# Left sidebar (vertical nav, supports section headers)
# ------------------------------------------------------------------
def build_sidebar(
    items: Iterable,
    current: str | None = None,
    *,
    label: str | None = None,
    suffix: str = '.html',
) -> str:
    """Build a vertical ``<aside class="report-sidebar">`` block.

    Items can be regular page links or section headers (rendered as group
    titles, not clickable). Use this together with :func:`render` and a
    ``<main>`` wrapper to produce a docs-style two-column layout::

        sidebar = build_sidebar([
            ('index', '总览'),
            ('section', '4 个 Finding'),
            ('find_A', 'Finding A · 13 candidates'),
            ('find_B', 'Finding B · Urea'),
            ('section', '附录'),
            ('appx_pipeline', '9 步流程'),
        ], current='find_A')

        report.render(md, css=css, ...,
                      header_html=sidebar + '<main class="report-main">',
                      footer_html='</main>')

    Parameters
    ----------
    items : iterable of tuples
        Each entry is one of:

        - ``(slug, display_label)`` — regular link, href = ``slug + suffix``
        - ``(slug, display_label, href)`` — link with explicit href
        - ``('section', section_label)`` — group header (no link)
    current : str, optional
        Slug of the current page.  Marked with ``class="current"``.
    label : str, optional
        Optional small label rendered above the item list (e.g. ``'章节'``).
    suffix : str
        File suffix appended to each 2-tuple slug to form the href.

    Returns
    -------
    str
        HTML string.
    """
    parts = ['<aside class="report-sidebar">']
    if label:
        parts.append(f'  <div class="sidebar-label">{label}</div>')
    for entry in items:
        if len(entry) >= 2 and entry[0] == 'section':
            parts.append(f'  <div class="sidebar-section">{entry[1]}</div>')
            continue
        if len(entry) == 3:
            slug, lbl, href = entry
        else:
            slug, lbl = entry
            href = f'{slug}{suffix}'
        cls = ' class="current"' if slug == current else ''
        parts.append(f'  <a href="{href}"{cls}>{lbl}</a>')
    parts.append('</aside>')
    return '\n'.join(parts)
