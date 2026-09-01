"""drawing_yh.report — 转发 shim(2026-08-28 起)

渲染引擎已唯一化到 `local-report-builder/report_toolkit`(独立公开 repo, 零绘图依赖).
本包只做转发, 保证存量项目脚本的 `from drawing_yh import report` 等引用零改动.
改引擎请去 report_toolkit 改, 不要在这里写实现.
"""
from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path

_TK = Path(os.environ.get(
    "REPORT_TOOLKIT_DIR",
    Path(__file__).resolve().parents[3] / "local-report-builder" / "report_toolkit",
))

try:
    _parent = str(_TK.parent)
    if _parent not in sys.path:
        sys.path.insert(0, _parent)
    from report_toolkit import core as _core  # noqa: E402
except Exception as _e:  # pragma: no cover
    raise ImportError(
        "drawing_yh.report 现为转发 shim, 需要 report_toolkit 可用. "
        f"查找位置: {_TK}. "
        "修复: 设置 REPORT_TOOLKIT_DIR 指向 local-report-builder/report_toolkit, "
        f"或 pip install -e local-report-builder. 原始错误: {_e}"
    ) from _e

# 原 __init__ 的 core 公开 API 全量转发
from report_toolkit.core import (  # noqa: E402,F401
    briefing_card,
    build_nav,
    build_sidebar,
    check_orphan_pages,
    combined_config_template_path,
    combined_template_path,
    config_template_path,
    css_path,
    list_templates,
    pptx_template_path,
    publish,
    render,
    serve_nocache_template_path,
    split_by_h2,
    watch_build_template_path,
)

# 子模块兼容: drawing_yh.report.ppt_engine / .legacy_ppt / .core → report_toolkit 对应模块
# (延迟: PPT 子模块只在被显式 import 时才加载, HTML 构建环境无 python-pptx 也能用)
sys.modules[__name__ + ".core"] = _core


def __getattr__(name):  # PEP 562
    if name in ("ppt_engine", "legacy_ppt"):
        import importlib
        mod = importlib.import_module(f"report_toolkit.{name}")
        sys.modules[__name__ + "." + name] = mod
        return mod
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
