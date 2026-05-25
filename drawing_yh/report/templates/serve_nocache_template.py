#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""serve_nocache.py — `python -m http.server` 套壳,每个响应加 no-cache headers。

为什么:`python -m http.server` 默认无 Cache-Control header,build_pages.py 重生成
HTML / 重做图后浏览器拿缓存看到旧版(尤其 sidebar 跳转时)。这里强制 no-store
让浏览器**永远从磁盘拉新文件**。

用法:`python serve_nocache.py [port]`(默认 8766)
"""
import sys
from functools import partial
from pathlib import Path
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer


class NoCacheHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate, max-age=0')
        self.send_header('Pragma', 'no-cache')
        self.send_header('Expires', '0')
        super().end_headers()


if __name__ == '__main__':
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8766
    # 默认服务脚本所在目录(报告所在 dir),避免依赖 cwd
    root = sys.argv[2] if len(sys.argv) > 2 else str(Path(__file__).resolve().parent)
    handler = partial(NoCacheHandler, directory=root)
    print(f'serve_nocache.py listening on :{port}  root={root}')
    ThreadingHTTPServer(('', port), handler).serve_forever()
