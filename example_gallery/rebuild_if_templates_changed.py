#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Rebuild the example gallery when drawing templates change.

Intended for the local git pre-commit hook. The hook checks staged files; if a
template/module/gallery source changed, this script rebuilds the gallery and
stages the generated HTML, manifest, and preview assets.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GALLERY_DIR = ROOT / "example_gallery"

WATCH_PREFIXES = (
    "drawing_yh/",
    "example_gallery/build_gallery.py",
    "example_gallery/rebuild_if_templates_changed.py",
)
WATCH_SUFFIXES = (".py", ".css", ".csv")

OUTPUT_PREFIXES = (
    "example_gallery/generated/",
    "example_gallery/assets/",
    "example_gallery/pages/",
)
OUTPUT_FILES = (
    "example_gallery/index.html",
    "example_gallery/gallery_manifest.csv",
)


def run_git(args: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", "-C", str(ROOT), *args],
        check=check,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def staged_files() -> list[str]:
    out = run_git(["diff", "--cached", "--name-only", "--diff-filter=ACMRT"]).stdout
    return [line.strip().replace("\\", "/") for line in out.splitlines() if line.strip()]


def should_rebuild(files: list[str]) -> bool:
    for path in files:
        if not path.endswith(WATCH_SUFFIXES):
            continue
        if path.startswith(OUTPUT_PREFIXES) or path in OUTPUT_FILES:
            continue
        if any(path.startswith(prefix) for prefix in WATCH_PREFIXES):
            return True
    return False


def stage_gallery_outputs() -> None:
    paths = [str(GALLERY_DIR / "index.html"), str(GALLERY_DIR / "gallery_manifest.csv")]
    for folder in ["generated", "assets", "pages"]:
        paths.append(str(GALLERY_DIR / folder))
    subprocess.run(["git", "-C", str(ROOT), "add", "--", *paths], check=True)


def main() -> int:
    files = staged_files()
    if not should_rebuild(files):
        print("gallery hook: no template/gallery source change; skip rebuild")
        return 0

    print("gallery hook: template/gallery source changed; rebuilding example gallery")
    subprocess.run([sys.executable, str(GALLERY_DIR / "build_gallery.py")], check=True, cwd=str(ROOT))
    stage_gallery_outputs()
    print("gallery hook: rebuilt and staged example_gallery outputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
