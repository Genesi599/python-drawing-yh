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
    return written
