"""统一排图工具 CLI(见包 docstring)。python -m drawing_yh.assembly.build_layout <layout.yaml>"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

MM = 3.7795275591  # mm → px @96dpi
SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)
ET.register_namespace("inkscape", "http://www.inkscape.org/namespaces/inkscape")
ET.register_namespace("sodipodi", "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd")
DEFAULT_INK = r"C:/Program Files/Inkscape/bin/inkscape.exe"

URL_RE = re.compile(r"url\(#([^)]+)\)")


def load_yaml(p: Path) -> dict:
    import yaml
    return yaml.safe_load(p.read_text(encoding="utf-8"))


def parse_len(v: str | None) -> float:
    if not v:
        return 0.0
    v = v.strip()
    for suf, fac in (("mm", MM), ("px", 1.0), ("pt", 96 / 72), ("cm", 10 * MM), ("in", 96)):
        if v.endswith(suf):
            return float(v[: -len(suf)]) * fac
    try:
        return float(v)
    except ValueError:
        return 0.0


def svg_box(src: Path) -> tuple[float, float, float, float]:
    """返回 (minx, miny, width, height) px。"""
    root = ET.parse(src).getroot()
    vb = root.get("viewBox")
    if vb:
        parts = vb.replace(",", " ").split()
        x, y, w, h = (float(v) for v in parts[:4])
        if w > 0 and h > 0:
            return x, y, w, h
    return 0.0, 0.0, parse_len(root.get("width")), parse_len(root.get("height"))


def ensure_svg(src: Path, cache: Path, ink: Path) -> Path:
    """PDF/位图统一转成可嵌 SVG;SVG 原样返回。"""
    if src.suffix.lower() == ".svg":
        return src
    cache.mkdir(parents=True, exist_ok=True)
    dst = cache / (src.stem + ".svg")
    if src.suffix.lower() == ".pdf":
        if not dst.exists() or dst.stat().st_mtime < src.stat().st_mtime:
            subprocess.run([str(ink), "--export-type=svg", "--export-filename", str(dst), str(src)],
                           check=True, capture_output=True)
        return dst
    w, h = _png_size(src)
    href = src.resolve().as_uri()
    root = ET.Element("svg", {
        "width": f"{w}px", "height": f"{h}px", "viewBox": f"0 0 {w} {h}"})
    ET.SubElement(root, "image", {"x": "0", "y": "0", "width": str(w), "height": str(h),
                                  f"{{{XLINK_NS}}}href": href, "href": href})
    ET.ElementTree(root).write(dst, encoding="unicode", xml_declaration=True)
    return dst


def _png_size(p: Path) -> tuple[int, int]:
    import struct
    with open(p, "rb") as f:
        head = f.read(24)
    if head[:8] == b"\x89PNG\r\n\x1a\n":
        return struct.unpack(">II", head[16:24])
    return 1000, 700


def embed_panel(dest: ET.Element, src_svg: Path, pid: str,
                x_px: float, y_px: float, w_px: float, h_px: float,
                vw: float, vh: float, minx: float, miny: float) -> float:
    """深度合并: 子 SVG children → <g translate+scale>,id 全部加前缀防冲突。返回 scale。"""
    root = ET.parse(src_svg).getroot()
    scale = w_px / vw
    idmap: dict[str, str] = {}
    for el in root.iter():
        old = el.get("id")
        if old:
            new = f"{pid}__{old}"
            idmap[old] = new
            el.set("id", new)
    for el in root.iter():
        for attr, val in list(el.attrib.items()):
            if "url(#" in val:
                el.set(attr, URL_RE.sub(lambda m: f"url(#{idmap.get(m.group(1), m.group(1))})", val))
            if (attr.endswith("href") or attr.endswith("}href")) and val.startswith("#"):
                el.set(attr, "#" + idmap.get(val[1:], val[1:]))
    tf = f"translate({x_px - minx * scale:.3f},{y_px - miny * scale:.3f})"
    if abs(scale - 1.0) > 1e-6:
        tf += f" scale({scale:.6f})"
    g = ET.SubElement(dest, "g", {"id": pid, "transform": tf})
    DROP = ("namedview", "page", "pageeffect", "pagegroup")
    for child in list(root):
        tag = child.tag.rsplit("}", 1)[-1]
        if tag in DROP:
            continue  # Inkscape 多页标记不能进主文档(会被当独立页导出)
        g.append(child)
    return scale


def build_page(cfg_page: dict, yaml_dir: Path, ink: Path, out_dir: Path, name: str) -> list[Path]:
    """构建一个 YAML page 配置。flow 模式超页自动拆临时下一页(2026-09-03 用户指示),
    返回全部页的 svg 路径(可能 >1)。"""
    page = cfg_page["page"]
    gut = cfg_page.get("gutters", {})
    margin = gut.get("margin_mm", 10) * MM
    label_gap = gut.get("label_gap_mm", 4) * MM
    cap_gap = gut.get("caption_gap_mm", 3.5) * MM
    col_gap = gut.get("col_gap_mm", 6) * MM
    row_gap = gut.get("row_gap_mm", 10) * MM
    mode = cfg_page.get("mode", "flow" if cfg_page.get("flow", False) else "abs")
    cache = out_dir / ".cache"

    # 1) 源尺寸 + 有意缩放覆盖
    items = []
    for p in cfg_page.get("panels", []):
        src = yaml_dir / p["src"]
        if not src.exists():
            raise SystemExit(f"panel src 不存在: {src}")
        sv = ensure_svg(src, cache, ink)
        minx, miny, vw, vh = svg_box(sv)
        if vw <= 0 or vh <= 0:
            raise SystemExit(f"SVG 尺寸读不到: {sv}")
        if "height_mm" in p:
            h_px = p["height_mm"] * MM
            w_px = h_px * vw / vh
        elif "width_mm" in p:
            w_px = p["width_mm"] * MM
            h_px = w_px * vh / vw
        else:
            w_px, h_px = vw, vh  # 自然尺寸
        items.append({"p": p, "sv": sv, "vw": vw, "vh": vh, "minx": minx, "miny": miny,
                      "w_px": w_px, "h_px": h_px})

    grid_rows: list[tuple[str, list[dict]]] = []
    grid_title = None

    # 2) 定位
    page_chunks: list[tuple[float, float, list[dict], list, str | None]] = []
    if mode == "flow":
        target_w = max(page["width_mm"] * MM, max(it["w_px"] for it in items) + 2 * margin)
        if target_w > page["width_mm"] * MM + 0.5:
            wide = [it["p"].get("id", it["sv"].stem) for it in items
                    if it["w_px"] + 2 * margin > page["width_mm"] * MM]
            print(f"[warn] flow 页超宽: {target_w / MM:.0f}mm > 页宽 {page['width_mm']}mm, "
                  f"超宽 panel: {', '.join(wide)}", file=sys.stderr)
        rows, cur, cur_x = [], [], margin
        for it in items:
            if cur and cur_x + col_gap + it["w_px"] > target_w - margin:
                rows.append(cur)
                cur, cur_x = [], margin
            it["x_px"], it["y_px"] = cur_x, 0.0
            cur.append(it)
            cur_x += (col_gap if len(cur) > 1 else 0) + it["w_px"]
        rows.append(cur)
        W = target_w
        page_h = page["height_mm"] * MM
        usable_bottom = page_h - margin
        max_body = page_h - 2 * margin
        # 逐行装页: 放不下 → 临时下一页(2026-09-03 用户指示)
        chunk_rows: list[list[list[dict]]] = []
        cur_rows: list[list[dict]] = []
        cur_y = margin
        for row in rows:
            body = max(it["h_px"] for it in row)
            rh = body + (label_gap + 10) + cap_gap + row_gap
            if body > max_body:
                ids = ", ".join(it["p"].get("id", it["sv"].stem) for it in row)
                print(f"[warn] 单行超页高: {ids} 图体 {body / MM:.0f}mm > 可用 {max_body / MM:.0f}mm, "
                      f"仍独占一页, 底部可能裁切", file=sys.stderr)
            if cur_rows and cur_y + rh > usable_bottom:
                chunk_rows.append(cur_rows)
                cur_rows, cur_y = [], margin
            cur_rows.append(row)
            cur_y += rh
        chunk_rows.append(cur_rows)
        if len(chunk_rows) > 1:
            print(f"[info] flow 自动分页: {len(items)} panel → {len(chunk_rows)} 页", file=sys.stderr)
        for prows in chunk_rows:
            y = margin
            for row in prows:
                for it in row:
                    it["y_px"] = y + (label_gap + 10)
                y += max(it["h_px"] for it in row) + (label_gap + 10) + cap_gap + row_gap
            page_chunks.append((W, page_h, [it for row in prows for it in row], [], None))
    elif mode == "grid":
        g = cfg_page.get("grid", {})
        cell_w = g.get("cell_w_mm", 60) * MM
        cell_h = g.get("cell_h_mm", 60) * MM
        label_w = g.get("label_w_mm", 15) * MM
        grid_title = g.get("title")
        title_h = (g.get("title_h_mm", 12) if grid_title else 0) * MM
        enlarge = g.get("enlarge", False)
        for it in items:  # 按声明顺序分行(row 值变化即新行)
            rlab = str(it["p"].get("row", ""))
            if not grid_rows or grid_rows[-1][0] != rlab:
                grid_rows.append((rlab, []))
            grid_rows[-1][1].append(it)
        cols = g.get("cols") or max(len(r[1]) for r in grid_rows)
        W = 2 * margin + label_w + cols * cell_w
        H = 2 * margin + title_h + len(grid_rows) * cell_h
        for ri, (_rlab, row_items) in enumerate(grid_rows):
            for ci, it in enumerate(row_items):
                sc = min(cell_w / it["vw"], cell_h / it["vh"])
                if not enlarge:
                    sc = min(sc, 1.0)  # thumbnail 语义: 只缩不放
                w2, h2 = it["vw"] * sc, it["vh"] * sc
                it["w_px"], it["h_px"] = w2, h2
                it["x_px"] = margin + label_w + ci * cell_w + (cell_w - w2) / 2
                it["y_px"] = margin + title_h + ri * cell_h + (cell_h - h2) / 2
        page_chunks.append((W, H, items, grid_rows, grid_title))
    else:  # abs
        W = page["width_mm"] * MM
        H = page["height_mm"] * MM
        for it in items:
            p = it["p"]
            if "x_mm" not in p or "y_mm" not in p:
                raise SystemExit(f"abs 模式 panel 需 x_mm/y_mm: {p.get('id')}")
            it["x_px"], it["y_px"] = p["x_mm"] * MM, p["y_mm"] * MM
        page_chunks.append((W, H, items, [], None))

    # 3) 渲染每页
    svg_paths = []
    for ci, (W, H, chunk_items, chunk_grid_rows, chunk_title) in enumerate(page_chunks, 1):
        svg = ET.Element("svg", {
            "width": f"{W / MM:.1f}mm", "height": f"{H / MM:.1f}mm",
            "viewBox": f"0 0 {W:.2f} {H:.2f}"})
        ET.SubElement(svg, "rect", {"x": "0", "y": "0", "width": f"{W:.2f}", "height": f"{H:.2f}",
                                    "fill": page.get("background", "white")})
        if chunk_title:
            ET.SubElement(svg, "text", {"x": f"{W / 2:.1f}", "y": f"{margin + 8:.1f}",
                                        "text-anchor": "middle", "font-family": "Arial",
                                        "font-weight": "bold", "font-size": "28"}).text = chunk_title
        for rlab, row_items in chunk_grid_rows:
            if rlab and row_items:
                ymid = row_items[0]["y_px"] + row_items[0]["h_px"] / 2
                ET.SubElement(svg, "text", {"x": f"{margin:.1f}", "y": f"{ymid:.1f}",
                                            "font-family": "Arial", "font-weight": "bold",
                                            "font-size": "24", "fill": "rgb(120,0,0)"}).text = rlab
        for it in chunk_items:
            p = it["p"]
            pid = p.get("id", it["sv"].stem)
            embed_panel(svg, it["sv"], pid, it["x_px"], it["y_px"], it["w_px"], it["h_px"],
                        it["vw"], it["vh"], it["minx"], it["miny"])
            if p.get("label"):
                ET.SubElement(svg, "text", {"x": f"{it['x_px']:.1f}",
                                            "y": f"{it['y_px'] - label_gap:.1f}",
                                            "font-family": "Arial", "font-weight": "bold",
                                            "font-size": "28"}).text = str(p["label"])
            if p.get("caption"):
                ET.SubElement(svg, "text", {"x": f"{it['x_px']:.1f}",
                                            "y": f"{it['y_px'] + it['h_px'] + cap_gap:.1f}",
                                            "font-family": "Arial", "font-size": "13"}).text = str(p["caption"])
            if mode == "abs":
                if it["y_px"] + it["h_px"] > H or it["x_px"] + it["w_px"] > W:
                    print(f"[warn] {pid} 越界: 底 {(it['y_px'] + it['h_px']) / MM:.0f}mm / 页高 {H/MM:.0f}mm,"
                          f" 右 {(it['x_px'] + it['w_px']) / MM:.0f}mm / 页宽 {W/MM:.0f}mm", file=sys.stderr)
        out_dir.mkdir(exist_ok=True)
        svg_path = out_dir / f"{name}_raw{ci}.svg"
        ET.ElementTree(svg).write(svg_path, encoding="unicode", xml_declaration=True)
        svg_paths.append(svg_path)
    return svg_paths


def export(svg_path: Path, ink: Path, formats=("pdf", "png")) -> list[Path]:
    outs = []
    for ft in formats:
        dst = svg_path.with_suffix("." + ft)
        r = subprocess.run([str(ink), "--export-type", ft, "--export-area-page",
                            "--export-filename", str(dst), str(svg_path)],
                           capture_output=True, text=True)
        if not dst.exists():
            print(f"[warn] {ft} 导出失败: {r.stderr[-200:]}", file=sys.stderr)
        else:
            outs.append(dst)
    return outs


def main() -> int:
    ap = argparse.ArgumentParser(prog="drawing_yh.assembly.build_layout")
    ap.add_argument("layout", help="布局 YAML 路径")
    ap.add_argument("--ink", default=DEFAULT_INK, help="inkscape.exe 路径")
    args = ap.parse_args()
    yaml_p = Path(args.layout).resolve()
    ink = Path(args.ink)
    if not ink.exists():
        raise SystemExit(f"找不到 Inkscape: {ink}")
    cfg = load_yaml(yaml_p)
    out_dir = yaml_p.parent / cfg.get("out", "out_layout")
    name = cfg.get("name", "layout")
    # 清旧输出(保留 .cache)
    if out_dir.exists():
        for old in out_dir.glob(f"{name}*"):
            if old.is_file():
                old.unlink()
    pages = cfg.get("pages") or [{"page": cfg["page"], "gutters": cfg.get("gutters", {}),
                                  "mode": "flow" if cfg.get("flow", False) else "abs",
                                  "panels": cfg.get("panels", [])}]
    raw_paths = []
    for pg in pages:
        # 每 YAML page 给唯一前缀, 防多 page 的 rawN 重名覆盖
        raw_paths.extend(build_page(pg, yaml_p.parent, ink, out_dir, f"{name}__pg{len(raw_paths)}"))
    final = []
    for k, sp in enumerate(raw_paths, 1):
        dst = out_dir / ((name if len(raw_paths) == 1 else f"{name}_p{k}") + ".svg")
        sp.replace(dst)
        final.append(dst)
    for sp in final:
        outs = export(sp, ink)
        print("SVG:", sp)
        for o in outs:
            print(o.suffix.lstrip(".").upper() + ":", o)
    return 0


if __name__ == "__main__":
    sys.exit(main())
