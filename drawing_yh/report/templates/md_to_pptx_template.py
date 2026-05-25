#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
md_to_pptx_template.py — REPORT.md → native PPTX,继承一份 PPTX 母版的 theme。

复制到项目目录改名 `md_to_pptx.py`。配置(SLIDES)从 `_report_config.py` import。
**改这 4 处即可**:
  1. PROJECT_LABEL —— 顶部 header 左上的小字(类似 "Retinal (Macaca, bulk)")
  2. TEMPLATE      —— 母版 .pptx 路径(用它的 master / 字体 / 主题色;原文件不动)
  3. OUT           —— 输出 .pptx 路径
  4. main() 里的 title_slide 标题 + 副标题文案

布局:每张 slide 顶部 = 项目标签 + slide 标题 + 橙色分隔线;body 自动选择
  - 有图(`<img>`):左大图 + 右 ➤ bullets(垂直居中)
  - 无图:全宽 ➤ bullets / 表
约定:有图 slide 不放表;表格列宽按内容长短分配;字号 16pt(表 14pt header / 11–14pt body);
      内容超 → 自动截字号 / 截行 / 丢次要 list,**不**显式提示。

注:本模板假设母版尺寸 16:9(13.33 × 7.5 in)。其他比例改下方 SLIDE_W / SLIDE_H。
"""
import re, sys
from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.enum.shapes import MSO_SHAPE
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from PIL import Image

# ----- 改这里 -----
PROJECT_LABEL  = 'My Project × Subtopic'                   # 改成你项目的标签
DIVIDER_COLOR  = RGBColor(0xBE, 0x76, 0x00)                # header 下方分隔线颜色
TEMPLATE = Path(r'D:/path/to/your/template.pptx')          # 母版 PPTX
OUT      = Path(r'D:/path/to/output/REPORT.pptx')          # 输出
# -------------------

# 共享配置(纯 import 不触发任何渲染副作用)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _report_config import HERE, SLIDES, load_chunks   # noqa: E402

# 母版 16:9(13.33 × 7.5 in),hardcode 避免 module-level 创建 Presentation
SLIDE_W = Inches(13.333)
SLIDE_H = Inches(7.500)

# prs / 各 layout 在 main() 里加载,通过 _PRS_CTX 给 add_content_slide 用
_PRS_CTX: dict = {}     # {'blank_layout': layout 对象}


def _load_template() -> Presentation:
    """读母版 + 真删原 slide(连 part + rels + 图片一起删,否则 PowerPoint
    严格校验会报"内容有问题")。继承 master / theme / 字体配色。"""
    prs = Presentation(str(TEMPLATE))
    sld_id_lst = prs.slides._sldIdLst
    for sld_id in list(sld_id_lst):
        prs.part.drop_rel(sld_id.rId)
        sld_id_lst.remove(sld_id)
    return prs

# ------------------------------------------------------------------
# 2. Markdown chunk → 块结构(para / list / table / image / code)
# ------------------------------------------------------------------
def parse_chunk(md: str):
    """返回 (heading, [block, ...])。heading 用 chunk 第一个 H2/H3。"""
    h_m = re.search(r'^(#{2,3})\s+(.+?)\s*$', md, re.MULTILINE)
    heading = h_m.group(2).strip() if h_m else ''
    body = (md[h_m.end():] if h_m else md).strip()

    blocks = []
    for part in re.split(r'\n\s*\n', body):
        part = part.strip()
        if not part:
            continue
        # 跳过 markdown 横线(---、***、___)
        if re.fullmatch(r'[-*_]{3,}', part):
            continue
        # 跳过 image 后面的 *italic caption*(我们的左大图布局不需要再放图注)
        if re.fullmatch(r'\*[^*]+\*', part):
            continue
        if part.startswith('|') and re.search(r'\n\|', part):
            blocks.append({'type': 'table', 'md': part})
        elif re.match(r'<img\s', part, re.I):
            m = re.search(r'src="([^"]+)"', part)
            if m:
                blocks.append({'type': 'image', 'src': m.group(1)})
        elif part.startswith('```'):
            inner = re.sub(r'^```\w*\n?', '', part)
            inner = re.sub(r'\n?```$', '', inner)
            blocks.append({'type': 'code', 'text': inner})
        elif re.match(r'[-*]\s', part):
            items = [re.sub(r'^[-*]\s+', '', l).strip()
                     for l in part.splitlines() if l.strip().startswith(('-', '*'))]
            blocks.append({'type': 'list', 'items': items})
        else:
            blocks.append({'type': 'para', 'text': part})
    return heading, blocks


# ------------------------------------------------------------------
# 3. Helper:在 slide 上的 (cur_y) 位置追加各种块
# ------------------------------------------------------------------
def resolve_image(src: str) -> Path | None:
    """REPORT.md 里的 'report_figs/foo.png' / '../report_figs/foo.png' 等都解析。"""
    name = Path(src).name
    candidates = [
        HERE / 'report_figs' / name,
        Path(r'D:/Projects/Bone_Marrow_Aging/cellchat_x_BMIF/figure') / name,
        Path(r'D:/Projects/Bone_Marrow_Aging/cellchat_x_BMIF/report_figs') / name,
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


BODY_LEFT   = Inches(0.30)
BODY_RIGHT  = Inches(0.30)
BODY_TOP    = Inches(1.00)              # 在 divider(y=0.81)之下留 gap
BODY_W      = SLIDE_W - BODY_LEFT - BODY_RIGHT
BODY_BOTTOM = SLIDE_H - Inches(0.25)


# ------------------------------------------------------------------
# Header layout 仿 retina-bulk:左上项目标签 + 中右 slide 标题 + 横分隔线
# ------------------------------------------------------------------
def add_header(slide, slide_title: str):
    # 左上 项目标签
    pl = slide.shapes.add_textbox(Inches(0.06), Inches(0.41), Inches(5.0), Inches(0.40))
    pl.text_frame.text = PROJECT_LABEL
    pl.text_frame.margin_top = pl.text_frame.margin_bottom = 0
    p = pl.text_frame.paragraphs[0]
    for r in p.runs:
        r.font.name = '微软雅黑'
        r.font.size = Pt(18)
    # 右半 slide 标题(bigger + bold)
    st = slide.shapes.add_textbox(Inches(5.20), Inches(0.13), Inches(8.00), Inches(0.55))
    st.text_frame.text = slide_title
    st.text_frame.margin_top = st.text_frame.margin_bottom = 0
    p = st.text_frame.paragraphs[0]
    for r in p.runs:
        r.font.name = '微软雅黑'
        r.font.size = Pt(28)
        r.font.bold = True
    # 横分隔线 #BE7600,full width × 0.04 in
    div = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, Inches(0.81), SLIDE_W, Inches(0.04))
    div.fill.solid()
    div.fill.fore_color.rgb = DIVIDER_COLOR
    div.line.fill.background()  # 无边框


def add_paragraph(slide, txt: str, x, y, w, *, font_size=16, bold=False):
    # 16pt 中文一行约 25–30 字宽 5 in;每行 ~0.32 in 高
    chars_per_line = max(15, int(w / Emu(914400) * 5))   # rough
    n_lines = max(1, len(txt) // chars_per_line + 1)
    h = Inches(max(0.4, 0.34 * n_lines + 0.18))
    tb = slide.shapes.add_textbox(x, y, w, h)
    tf = tb.text_frame
    tf.word_wrap = True
    tf.text = txt
    for p in tf.paragraphs:
        for r in p.runs:
            r.font.name = '微软雅黑'
            r.font.size = Pt(font_size)
            r.font.bold = bold
    return h


def add_list(slide, items, x, y, w, *, font_size=14):
    n = len(items)
    h = Inches(max(0.45, 0.32 * n + 0.2))
    tb = slide.shapes.add_textbox(x, y, w, h)
    tf = tb.text_frame
    tf.word_wrap = True
    for i, it in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = '• ' + it
        for r in p.runs:
            r.font.name = '微软雅黑'
            r.font.size = Pt(font_size)
    return h


def _bullet_lines(text: str, w_in: float, font_size: int) -> int:
    """估算文字在宽度 w_in (英寸) 的 textbox 里能占多少行。
    标定:微软雅黑 16pt,5 in 宽 ≈ 19 中文字/行 → 3.8 字/(in × 16pt)
    其他字号线性外推:chars_per_line = w_in × 3.8 × 16/font_size。"""
    chars_per_line = max(2, int(w_in * 3.8 * 16 / font_size))
    return max(1, (len(text) + chars_per_line - 1) // chars_per_line)


def _parse_md_table(md_block: str):
    """返回 (cells, n_cols)。cells 是 list of list of str。"""
    rows = [r for r in md_block.splitlines() if r.strip().startswith('|')
            and not re.match(r'^\|\s*[-:|\s]+\|\s*$', r)]
    cells = [[c.strip() for c in r.strip().strip('|').split('|')] for r in rows]
    n_cols = max((len(r) for r in cells), default=0)
    return cells, n_cols


def _compute_col_widths_in(cells, total_w_in: float, n_cols: int,
                            min_col_w_in: float = 0.6) -> list[float]:
    """按每列 max(len(cell)) 加权分配列宽,带 min 兜底。
    充分利用版面 → 短列窄、长列宽,减少不必要的 wrap。"""
    if n_cols == 0: return []
    max_chars = [1] * n_cols
    for row in cells:
        for ci in range(min(len(row), n_cols)):
            max_chars[ci] = max(max_chars[ci], len(row[ci]))
    # min 兜底,但不能加起来超过总宽
    min_alloc = min(min_col_w_in, total_w_in / n_cols)
    total_min = min_alloc * n_cols
    if total_min >= total_w_in:
        return [total_w_in / n_cols] * n_cols
    remaining = total_w_in - total_min
    total_wt = sum(max_chars)
    return [min_alloc + (max_chars[i] / total_wt) * remaining for i in range(n_cols)]


def _cells_h_in(cells, col_widths_in: list[float], font_pt: int) -> float:
    """按每列实际宽度算 wrap 高度(英寸)。"""
    if not cells: return 0.0
    line_h_in = (font_pt + 6) / 72.0
    total_h = 0.05
    for row in cells:
        max_lines = 1
        for ci in range(min(len(row), len(col_widths_in))):
            cw_in = col_widths_in[ci]
            chars_per_line = max(2, int(cw_in * 3.8 * 16 / font_pt))
            n = max(1, (len(row[ci]) + chars_per_line - 1) // chars_per_line)
            max_lines = max(max_lines, n)
        total_h += max_lines * line_h_in + 0.06
    return total_h


def estimate_table_h_in(md_block: str, w_in: float, font_pt: int = 14) -> float:
    cells, n_cols = _parse_md_table(md_block)
    if not cells: return 0.0
    col_widths = _compute_col_widths_in(cells, w_in, n_cols)
    return _cells_h_in(cells, col_widths, font_pt)


def can_fit_table(md_block: str, w_in: float, max_h_in: float,
                   min_font_pt: int = 11) -> bool:
    """min_font 也装不下 → False。"""
    return estimate_table_h_in(md_block, w_in, font_pt=min_font_pt) <= max_h_in


def fit_bullets(items, w_in: float, max_h_in: float, font_size: int = 16):
    """在 (w × max_h) 的盒子里能塞下几条箭头要点 + 每条多长。
    超出的:先尝试缩字数,实在装不下就丢弃。返回 (items_fit, est_h_in)。"""
    line_h = font_size / 72 * 1.45                  # 行高 in
    space_before = 6 / 72                            # 段前空 in
    used = 0.18                                      # 上下 padding
    out = []
    for it in items:
        s = it.strip()
        n_lines = _bullet_lines(s, w_in, font_size)
        item_h = line_h * n_lines + space_before
        if used + item_h > max_h_in:
            # 试试压短到能放进去
            rem_h = max_h_in - used - space_before
            if rem_h < line_h:
                break                                 # 没空间,直接停
            allowed_lines = max(1, int(rem_h / line_h))
            chars_per_line = max(2, int(w_in * 3.8 * 16 / font_size))
            allowed_chars = allowed_lines * chars_per_line - 2   # 留省略号
            if allowed_chars >= 8:
                out.append(s[:allowed_chars].rstrip() + '…')
                used += line_h * allowed_lines + space_before
            break
        out.append(s)
        used += item_h
    return out, used


def add_arrow_bullets_centered(slide, items, x, y, w, h, *, font_size=16):
    """➤ bullets,在 (w × h) 盒子内**垂直居中**(用 PPT 的 vertical_anchor)。"""
    if not items:
        return Inches(0)
    tb = slide.shapes.add_textbox(x, y, w, h)
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    for i, it in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = '➤ ' + it.strip()
        p.space_before = Pt(6)
        for r in p.runs:
            r.font.name = '微软雅黑'
            r.font.size = Pt(font_size)
    return h


def add_arrow_bullets(slide, items, x, y, w, *, font_size=16, max_h=None):
    """➤ 箭头要点。可选 max_h 自动裁切(防溢出)。"""
    w_in = w / 914400  # EMU → in
    if max_h is not None:
        items, est_used = fit_bullets(items, w_in, max_h / 914400, font_size=font_size)
        h = Inches(min(est_used + 0.05, max_h / 914400))
    else:
        # 无限制 → 老式估高
        h = Inches(max(0.45, 0.42 * len(items) + 0.20))
    if not items:
        return Inches(0)
    tb = slide.shapes.add_textbox(x, y, w, h)
    tf = tb.text_frame
    tf.word_wrap = True
    for i, it in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = '➤ ' + it.strip()
        p.space_before = Pt(6)
        for r in p.runs:
            r.font.name = '微软雅黑'
            r.font.size = Pt(font_size)
    return h


def add_md_table(slide, md_block: str, x, y, w, *, max_h=None, font_pt=14, min_font_pt=11):
    """画表,wrap-aware。
    流程:
      ① 在 [min_font, font_pt] 找最大可读字号能装下全表(含 wrap)
      ② 装不下 → 用 min_font 一行一行加,直到再加就溢出 → 末行 "(共 N 行,省略 M)"
      ③ header 至少要装下,否则返回 0(放弃这张表)
    """
    cells, n_cols = _parse_md_table(md_block)
    total_rows = len(cells)
    if not (total_rows and n_cols):
        return Inches(0)

    w_in = w / 914400
    # 关键:列宽按内容长短分配(短列窄、长列宽,充分利用版面)
    col_widths_in = _compute_col_widths_in(cells, w_in, n_cols)

    chosen_font = font_pt
    final_cells = cells

    if max_h is not None:
        max_h_in = max_h / 914400

        # ① 在大 → 小字号里找能装下全表的
        fits_full = False
        for fp in range(font_pt, min_font_pt - 1, -1):
            if _cells_h_in(cells, col_widths_in, fp) <= max_h_in:
                chosen_font = fp
                fits_full = True
                break

        # ② 全表装不下 → min font 截行
        if not fits_full:
            chosen_font = min_font_pt
            kept = [cells[0]]                       # 必带 header
            ellipsis_h = (chosen_font + 6) / 72.0 + 0.06
            for row in cells[1:]:
                trial = kept + [row]
                trial_h = _cells_h_in(trial, col_widths_in, chosen_font)
                if trial_h + ellipsis_h > max_h_in:
                    break
                kept = trial
            n_kept_data = len(kept) - 1
            n_dropped = (total_rows - 1) - n_kept_data
            if n_dropped > 0:
                ell = ['…'] * n_cols
                ell[1 if n_cols > 1 else 0] = f'(共 {total_rows - 1} 行,省略 {n_dropped})'
                kept.append(ell[:n_cols])
            # ③ 还是装不下(光 header 都超):整个放弃
            if _cells_h_in(kept, col_widths_in, chosen_font) > max_h_in:
                return Inches(0)
            final_cells = kept

    n_rows = len(final_cells)
    h = Inches(_cells_h_in(final_cells, col_widths_in, chosen_font))

    tbl = slide.shapes.add_table(n_rows, n_cols, x, y, w, h).table
    # 应用每列宽度(short cell narrow,long cell wide)
    for ci, cw_in in enumerate(col_widths_in):
        if ci < n_cols:
            tbl.columns[ci].width = Inches(cw_in)
    for ri, row in enumerate(final_cells):
        for ci in range(n_cols):
            val = row[ci] if ci < len(row) else ''
            cell = tbl.cell(ri, ci)
            cell.text = val
            for p in cell.text_frame.paragraphs:
                for r in p.runs:
                    r.font.name = '微软雅黑'
                    r.font.size = Pt(chosen_font + 1 if ri == 0 else chosen_font)
                    if ri == 0:
                        r.font.bold = True
    return h


def add_image(slide, img_path: Path, x, y, max_w, max_h):
    """嵌图,先用 Pillow 读真实像素,等比 fit 到 (max_w × max_h),
    一次 add_picture 给定确切宽高 —— 不留孤儿 relationship。"""
    with Image.open(img_path) as im:
        iw, ih = im.size
    aspect = ih / iw if iw else 1
    tw = int(max_w)
    th = int(tw * aspect)
    if th > int(max_h):
        th = int(max_h)
        tw = int(th / aspect)
    pic = slide.shapes.add_picture(str(img_path), x, y, width=tw, height=th)
    return pic.height


def add_code(slide, txt: str, x, y, w):
    n_lines = txt.count('\n') + 1
    h = Inches(max(0.4, 0.22 * n_lines + 0.2))
    tb = slide.shapes.add_textbox(x, y, w, h)
    tf = tb.text_frame
    tf.word_wrap = True
    tf.text = txt
    for p in tf.paragraphs:
        for r in p.runs:
            r.font.name = 'Consolas'
            r.font.size = Pt(11)
    return h


# ------------------------------------------------------------------
# 4. 一张 slide 的完整布局
# ------------------------------------------------------------------
# ------------------------------------------------------------------
# 双栏布局参数(有图时用):左大图 + 右箭头要点
# ------------------------------------------------------------------
LEFT_IMG_X   = Inches(0.30)
LEFT_IMG_W   = Inches(7.50)         # ~56% slide width
RIGHT_TXT_X  = Inches(7.95)
RIGHT_TXT_W  = Inches(5.10)         # ~38% slide width


def layout_image_left(slide, images, others):
    """有图 → 左 7.5 in 大图,右 5.1 in 箭头要点 + 表格。"""
    # 1. 左主图
    main_img = resolve_image(images[0]['src'])
    if main_img:
        avail_h = BODY_BOTTOM - BODY_TOP
        add_image(slide, main_img, LEFT_IMG_X, BODY_TOP, LEFT_IMG_W, avail_h)

    # 2. 右栏:把 para + list 全收编成 ➤ bullets,table 紧跟其后
    bullets, tables, codes = [], [], []
    for blk in others:
        if blk['type'] == 'para':   bullets.append(blk['text'])
        elif blk['type'] == 'list': bullets.extend(blk['items'])
        elif blk['type'] == 'table': tables.append(blk)
        elif blk['type'] == 'code':  codes.append(blk)

    # 右栏内容:**有图的 slide 一律不放表**(易冲突 / 拥挤)。
    # 表里的数据用 bullets 凝练,或读者去网页 / REPORT.md 看
    avail = BODY_BOTTOM - BODY_TOP                       # ~6.25 in
    avail_in = avail / 914400
    rt_w_in  = RIGHT_TXT_W / 914400

    fit_table = None
    dropped_n = len(tables)

    # 2a) 只有 bullets,没表 → 垂直居中
    if bullets and fit_table is None:
        # 留底部 ~0.3 in 给"省略"注脚
        anchor_h = avail - (Inches(0.35) if dropped_n else Inches(0))
        # 截断超长 bullet
        items_fit, _ = fit_bullets(bullets, rt_w_in, anchor_h / 914400, font_size=16)
        add_arrow_bullets_centered(slide, items_fit, RIGHT_TXT_X, BODY_TOP,
                                    RIGHT_TXT_W, anchor_h, font_size=16)
        cur_y = BODY_BOTTOM - (Inches(0.35) if dropped_n else Inches(0))

    # 2b) bullets + 1 张表 → bullets 顶上,表跟下面
    elif fit_table is not None:
        cur_y = BODY_TOP
        if bullets:
            cur_y += add_arrow_bullets(slide, bullets, RIGHT_TXT_X, cur_y, RIGHT_TXT_W,
                                        font_size=16, max_h=avail * 0.35) + Inches(0.1)
        rem = BODY_BOTTOM - cur_y - Inches(0.4 if dropped_n else 0.05)
        cur_y += add_md_table(slide, fit_table['md'], RIGHT_TXT_X, cur_y, RIGHT_TXT_W,
                               max_h=rem, font_pt=14, min_font_pt=11) + Inches(0.05)

    # 2c) 没 bullets 也没能放的表
    else:
        cur_y = BODY_BOTTOM

    # 4) 代码块(如果还有空间,塞右下)
    for cb in codes:
        if cur_y >= BODY_BOTTOM - Inches(0.4): break
        cur_y += add_code(slide, cb['text'], RIGHT_TXT_X, cur_y, RIGHT_TXT_W) + Inches(0.05)

    # 3. 第 2+ 张图(若有):塞到左下方主图之下,或者干脆跳过
    # 大多数 slide 1 张图就够;find_A 这种 2 张的暂时只显第 1 张,正文里有 venn 图就给 venn 让位
    # (后续要的话可以做"主图换 venn / 加缩略图")


def _estimate_height_in(blk, w_in: float = 12.73) -> float:
    """估高(英寸),用 wrap-aware 算法预判会不会溢出。"""
    if blk['type'] == 'list':
        return max(0.45, 0.40 * len(blk['items']) + 0.25)
    if blk['type'] == 'para':
        return max(0.4, 0.04 * (len(blk['text']) // 70 + 1) + 0.32)
    if blk['type'] == 'table':
        return estimate_table_h_in(blk['md'], w_in, font_pt=14)
    if blk['type'] == 'code':
        return max(0.4, 0.22 * (blk['text'].count('\n') + 1) + 0.2)
    return 0.5


def layout_full_width(slide, blocks):
    """无图 → body 全宽。**有表的 slide 表为主**:
       - bullet list 整个丢(那些"观察 / 评论"放网页/md 看)
       - paragraph 截到 60 字(基本 1 行作 label)
       - 表用 wrap-aware 自动选字号,装不下截行
    """
    n_tables = sum(1 for b in blocks if b['type'] == 'table')
    table_dominated = n_tables > 0
    para_max_chars = 60 if table_dominated else 1000
    body_w_in = BODY_W / 914400

    # 表为主时 → 直接过滤 list 块(整段丢)
    n_dropped_lists = 0
    if table_dominated:
        n_dropped_lists = sum(1 for b in blocks if b['type'] == 'list')
        blocks = [b for b in blocks if b['type'] != 'list']

    seen_tbl = 0
    dropped_tbl = 0
    cur_y = BODY_TOP

    for blk in blocks:
        rem = BODY_BOTTOM - cur_y - Inches(0.05)
        if rem < Inches(0.4):
            if blk['type'] == 'table': dropped_tbl += 1
            continue
        if blk['type'] == 'list':
            cur_y += add_arrow_bullets(slide, blk['items'], BODY_LEFT, cur_y, BODY_W,
                                        font_size=16, max_h=rem) + Inches(0.06)
        elif blk['type'] == 'para':
            text = blk['text']
            if len(text) > para_max_chars:
                text = text[:para_max_chars - 1].rstrip() + '…'
            items, _ = fit_bullets([text], body_w_in, rem / 914400, font_size=16)
            if items:
                cur_y += add_arrow_bullets(slide, items, BODY_LEFT, cur_y, BODY_W,
                                            font_size=16, max_h=rem) + Inches(0.04)
        elif blk['type'] == 'table':
            if seen_tbl >= 2:
                dropped_tbl += 1
                continue
            h_used = add_md_table(slide, blk['md'], BODY_LEFT, cur_y, BODY_W,
                                   max_h=rem, font_pt=14, min_font_pt=10)
            if h_used == 0:
                dropped_tbl += 1
            else:
                cur_y += h_used + Inches(0.08)
                seen_tbl += 1
        elif blk['type'] == 'code':
            cur_y += add_code(slide, blk['text'], BODY_LEFT, cur_y, BODY_W) + Inches(0.05)



def add_content_slide(prs, label: str, md_chunk: str):
    slide = prs.slides.add_slide(_PRS_CTX['blank_layout'])
    add_header(slide, label)

    _, blocks = parse_chunk(md_chunk)
    images = [b for b in blocks if b['type'] == 'image']
    others = [b for b in blocks if b['type'] != 'image']

    if images:
        layout_image_left(slide, images, others)
    else:
        layout_full_width(slide, others)


# ------------------------------------------------------------------
# 5. main()
# ------------------------------------------------------------------
def main():
    prs = _load_template()
    _PRS_CTX['blank_layout'] = prs.slide_layouts[6]   # 空白
    title_layout              = prs.slide_layouts[0]   # 标题幻灯片
    chunks = load_chunks(fix_img_for_pages=False)      # PPTX 跟 REPORT.md 同根,不要改图路径

    # 5a. 标题页 —— 改这里
    title_slide = prs.slides.add_slide(title_layout)
    title_slide.shapes.title.text = '<项目名> × <子主题>'
    if len(title_slide.placeholders) > 1:
        sub = title_slide.placeholders[1]
        sub.text = '<作者> · <单位>\n<日期>'

    # 5b. 16 张 content slide
    for slug, _, dlabel, _ in SLIDES:
        add_content_slide(prs, dlabel, chunks[slug])

    # 5c. 保存(被 PowerPoint 锁住时落到带时间戳的备用名)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    try:
        prs.save(str(OUT))
        target = OUT
    except PermissionError:
        from datetime import datetime
        target = OUT.with_stem(OUT.stem + '_' + datetime.now().strftime('%H%M%S'))
        prs.save(str(target))
        print(f'WARN: {OUT} is locked (PowerPoint 还开着?), saved to alt path:')
    print(f'Wrote: {target}')
    print(f'  slides: {len(prs.slides)}')
    print(f'  size:   {target.stat().st_size:,} B')


if __name__ == '__main__':
    main()
