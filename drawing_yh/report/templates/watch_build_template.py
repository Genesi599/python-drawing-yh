#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""watch_build.py — 监听 REPORT.md / _report_config.py / report_figs/ 变动,
自动重跑 build_pages.py,实现改 md → 浏览器自动看到最新 HTML。

配合 serve_nocache.py 使用(它加了 no-cache headers,浏览器刷新即拿新文件)。

用法:
    python watch_build.py                     # 监听当前目录
    python watch_build.py /path/to/report     # 监听指定目录
    python watch_build.py . --debounce 3      # 去抖间隔 3 秒(默认 2)

依赖: pip install watchdog
"""
import os
import subprocess
import sys
import time
import threading
from pathlib import Path

from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

# ── 配置 ──────────────────────────────────────────────
# 监听哪些文件的变动触发 rebuild
WATCH_PATTERNS = {
    'REPORT.md',
    '_report_config.py',
    '_combined_report_config.py',
}
# 监听哪些目录下的图片变动触发 rebuild
WATCH_DIRS = {'report_figs', 'figures'}
# 图片后缀
IMG_SUFFIXES = {'.png', '.jpg', '.jpeg', '.svg', '.pdf'}


class RebuildHandler(FileSystemEventHandler):
    """文件变动 → 去抖 → 自动 build_pages.py"""

    def __init__(self, root: Path, debounce: float = 2.0):
        super().__init__()
        self.root = root
        self.debounce = debounce
        self._timer = None
        self._lock = threading.Lock()
        # 找 build 脚本:优先 build_combined.py(combined 报告),否则 build_pages.py
        combined = root / 'build_combined.py'
        single = root / 'build_pages.py'
        if combined.exists():
            self.build_script = combined
        elif single.exists():
            self.build_script = single
        else:
            print(f'[watch] 警告: {root} 下找不到 build_pages.py 或 build_combined.py')
            self.build_script = single  # 兜底,跑时会报错

    def _should_trigger(self, path: str) -> bool:
        """判断这个文件变动是否应该触发 rebuild"""
        p = Path(path)
        name = p.name
        # 直接匹配的文件名
        if name in WATCH_PATTERNS:
            return True
        # report_figs / figures 目录下的图片
        parts = p.parts
        for wd in WATCH_DIRS:
            if wd in parts and p.suffix.lower() in IMG_SUFFIXES:
                return True
        return False

    def _schedule_rebuild(self):
        """去抖:连续变动只触发一次 rebuild"""
        with self._lock:
            if self._timer is not None:
                self._timer.cancel()
            self._timer = threading.Timer(self.debounce, self._do_rebuild)
            self._timer.start()

    def _do_rebuild(self):
        ts = time.strftime('%H:%M:%S')
        print(f'[watch {ts}] 检测到变动,重建 → python {self.build_script.name}')
        try:
            result = subprocess.run(
                [sys.executable, str(self.build_script)],
                cwd=str(self.root),
                capture_output=True,
                text=True,
                timeout=120,
            )
            if result.returncode == 0:
                print(f'[watch {ts}] 重建完成 ✓')
            else:
                print(f'[watch {ts}] 重建失败 (exit {result.returncode})')
                if result.stderr:
                    # 只打最后 10 行,不刷屏
                    for line in result.stderr.strip().splitlines()[-10:]:
                        print(f'  {line}')
        except subprocess.TimeoutExpired:
            print(f'[watch {ts}] 重建超时 (>120s)')
        except Exception as e:
            print(f'[watch {ts}] 重建异常: {e}')

    def on_modified(self, event):
        if not event.is_directory and self._should_trigger(event.src_path):
            self._schedule_rebuild()

    def on_created(self, event):
        if not event.is_directory and self._should_trigger(event.src_path):
            self._schedule_rebuild()


def main():
    import argparse
    ap = argparse.ArgumentParser(description='监听报告文件变动,自动 rebuild HTML')
    ap.add_argument('root', nargs='?', default='.',
                    help='报告根目录(含 REPORT.md + build_pages.py)')
    ap.add_argument('--debounce', type=float, default=2.0,
                    help='去抖间隔秒数(默认 2)')
    args = ap.parse_args()

    root = Path(args.root).resolve()
    handler = RebuildHandler(root, debounce=args.debounce)
    observer = Observer()
    observer.schedule(handler, str(root), recursive=True)
    observer.start()

    print(f'[watch] 监听 {root}')
    print(f'[watch] 触发文件: {WATCH_PATTERNS}')
    print(f'[watch] 触发目录: {WATCH_DIRS} (图片)')
    print(f'[watch] 去抖: {args.debounce}s')
    print(f'[watch] build: {handler.build_script.name}')
    print(f'[watch] Ctrl+C 停止')
    print()

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
        print('\n[watch] 已停止')
    observer.join()


if __name__ == '__main__':
    main()
