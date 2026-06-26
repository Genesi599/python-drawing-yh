"""
drawing_yh.io — 标准化的图片保存接口

按项目科研出图标准(详见包根 STANDARDS.md):
- `.svg` 自动用 dpi=72(SVG 字号 1:1 渲染必须)
- `.png` 默认 dpi=300(印刷质量)
- `.pdf` 矢量,不带 dpi
- 默认 `bbox_inches='tight'` + `facecolor='white'` + 去除 metadata
- 一次调用可输出多种后缀
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Union


SVG_DPI = 72                # SVG 必须 72 才能 font_size × 72/dpi = 1:1
PNG_DPI = 300               # PNG 印刷质量;期刊普遍要求 ≥ 300

# 这套 metadata 让导出的 PDF/PNG 不带 matplotlib + 系统签名,期刊提交友好
_CLEAN_METADATA = {'Creator': None, 'Producer': None}


def _dpi_for_suffix(suffix: str, dpi: Optional[int]) -> Optional[int]:
    """根据扩展名挑 dpi:`.svg` → 72,`.pdf` → None(矢量),其它 → 用户给的或 300。"""
    if dpi is not None:
        return dpi
    s = suffix.lower()
    if s == '.svg':
        return SVG_DPI
    if s == '.pdf':
        return None
    return PNG_DPI


def save_fig(
    fig,
    path: Union[str, Path],
    *,
    also: Iterable[str] = (),
    dpi: Optional[int] = None,
    transparent: bool = False,
    facecolor: Optional[str] = 'white',
    bbox_inches: Optional[str] = 'tight',
    pad_inches: Optional[float] = 0.02,
    clean_metadata: bool = True,
    save_data=None,
    save_code: bool = False,
    save_description=None,
    **savefig_kwargs,
) -> list:
    """
    按项目标准保存 figure。

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    path : str or Path
        输出路径,含扩展名(如 `out/fig1.pdf`)。如果还要再输出别的格式见 `also`。
    also : iterable of str
        额外要输出的格式扩展名,例如 `('.png', '.svg')`。**一次调用写多份。**
    dpi : int or None
        手动指定 dpi。None 时按扩展名自动选(`.svg`→72,`.pdf`→None,其它→300)。
    transparent : bool
        透明背景。True 时 `facecolor` 失效。
    facecolor : str or None
        默认白底,适合期刊投稿。
    bbox_inches, pad_inches : 控制裁白边。
    clean_metadata : bool
        是否给 PDF/PNG 加 `metadata={'Creator': None, 'Producer': None}` 去签名。
    save_data : DataFrame / dict / list-of-dict / callable / None
        作图标准 §11:图上**实际渲染的 plotted data**。给了就写 `<stem>.csv`(index=False)。
        数据密集图(bubble/heatmap/dotplot/多系列)应传,接力分析直接读、不必重跑脚本复现。
    save_code : bool
        True 时把**调用 save_fig 的脚本**拷成 `<stem>.snapshot.py`(READ-ONLY 横幅 + 源路径 + 时间戳);
        交互式(python -c / REPL)无源脚本时自动跳过。
    save_description : str / callable / None
        给了就写 `<stem>.md`(图的文字描述,让接力的 AI 不读图也知道内容)。
    **savefig_kwargs : 其它透传给 `fig.savefig`。

    Returns
    -------
    list[Path] : 实际写出的所有路径

    Examples
    --------
    >>> save_fig(fig, 'out/fig1.pdf', also=('.png', '.svg'))
    [PosixPath('out/fig1.pdf'), PosixPath('out/fig1.png'), PosixPath('out/fig1.svg')]
    """
    path = Path(path)
    base = path.with_suffix('')
    suffixes = list(dict.fromkeys([path.suffix] + list(also)))   # 去重保序
    if not any(suffixes):
        suffixes = ['.png']

    written = []
    for s in suffixes:
        if not s:
            continue
        p = base.with_suffix(s)
        kwargs = dict(savefig_kwargs)
        effective_dpi = _dpi_for_suffix(s, dpi)
        if effective_dpi is not None:
            kwargs.setdefault('dpi', effective_dpi)
        if bbox_inches is not None:
            kwargs.setdefault('bbox_inches', bbox_inches)
        if pad_inches is not None:
            kwargs.setdefault('pad_inches', pad_inches)
        if transparent:
            kwargs['transparent'] = True
        elif facecolor is not None:
            kwargs.setdefault('facecolor', facecolor)
        # SVG 不嵌 Creator/Producer,只对 PDF/PNG 加 metadata
        if clean_metadata and s.lower() in ('.pdf', '.png'):
            kwargs.setdefault('metadata', _CLEAN_METADATA)
        fig.savefig(p, **kwargs)
        written.append(p)

    # ── sidecar(作图标准 §11「5–6 文件一组」):plotted data / 代码副本 / 文字描述 ──
    if save_data is not None:
        import pandas as pd
        data = save_data() if callable(save_data) else save_data
        df = data if isinstance(data, pd.DataFrame) else pd.DataFrame(data)
        csv_p = base.with_suffix('.csv')
        df.to_csv(csv_p, index=False)
        written.append(csv_p)
    if save_code:
        snap = _write_code_snapshot(base)
        if snap is not None:
            written.append(snap)
    if save_description is not None:
        desc = save_description() if callable(save_description) else save_description
        md_p = base.with_suffix('.md')
        md_p.write_text(str(desc), encoding='utf-8')
        written.append(md_p)
    return written


def _write_code_snapshot(base: Path) -> Optional[Path]:
    """把调用 save_fig 的脚本拷成 <stem>.snapshot.py + READ-ONLY 横幅。

    交互式调用(无 .py 源,如 python -c / REPL)返回 None、不写。
    """
    import inspect
    import datetime
    caller = None
    for fr in inspect.stack():
        fn = (fr.filename or '').replace('\\', '/')
        if fn.endswith('.py') and '/drawing_yh/' not in fn:
            caller = Path(fr.filename)
            break
    if caller is None or not caller.exists():
        return None
    ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    banner = (
        '# ============================================================\n'
        '# READ-ONLY SNAPSHOT - 作图代码副本,只读\n'
        f'# 改图请改源脚本: {caller.resolve()}\n'
        '# 直接改本文件不会更新原脚本、也不会重新出图\n'
        f'# snapshot at: {ts}\n'
        '# ============================================================\n\n'
    )
    snap = base.with_name(base.name + '.snapshot.py')
    snap.write_text(banner + caller.read_text(encoding='utf-8'), encoding='utf-8')
    return snap
