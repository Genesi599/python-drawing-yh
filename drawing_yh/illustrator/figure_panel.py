#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : figure_panel.py
@Date    : 2026/4/8
@Author  : yh109

在 Adobe Illustrator 中创建论文整页组图模板
"""

import os
import time
import win32com.client

# ─────────────── 布局参数 ───────────────
ROWS = 3
COLS = 2
MARGIN_PT = 20
GAP_H_PT = 16
GAP_V_PT = 20
LABEL_FONT_SIZE = 10
SAVE_PATH = r"C:\Users\yh109\Desktop\figure_panel.ai"
# ────────────────────────────────────────

PAGE_W = 595
PAGE_H = 842
USABLE_W = PAGE_W - MARGIN_PT * 2
USABLE_H = PAGE_H - MARGIN_PT * 2
PANEL_W = (USABLE_W - GAP_H_PT * (COLS - 1)) / COLS
PANEL_H = (USABLE_H - GAP_V_PT * (ROWS - 1)) / ROWS


def make_rgb(r, g, b):
    c = win32com.client.Dispatch("Illustrator.RGBColor")
    c.Red = r
    c.Green = g
    c.Blue = b
    return c


def main():
    print("Connecting to Illustrator...")
    ai = win32com.client.Dispatch("Illustrator.Application")

    doc = ai.Documents.Add(2, PAGE_W, PAGE_H)
    time.sleep(1)

    layer = doc.Layers[0]
    layer.Name = "panels"

    labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    gray = make_rgb(230, 230, 230)
    stroke_gray = make_rgb(180, 180, 180)
    black = make_rgb(0, 0, 0)

    for i in range(ROWS * COLS):
        row = i // COLS
        col = i % COLS
        label = labels[i]

        left = MARGIN_PT + col * (PANEL_W + GAP_H_PT)
        top = PAGE_H - MARGIN_PT - row * (PANEL_H + GAP_V_PT)

        # 画占位框（和 main.py 一样的写法）
        rect = layer.PathItems.Rectangle(top, left, PANEL_W, PANEL_H)
        rect.FillColor = gray
        rect.StrokeColor = stroke_gray
        rect.StrokeWidth = 0.5

        # 添加标签
        tf = doc.TextFrames.Add()
        tf.Contents = label
        tf.Position = [left + 3, top - 2]
        tr = tf.TextRange
        ca = tr.CharacterAttributes
        ca.Size = LABEL_FONT_SIZE
        ca.FillColor = black

        print(f"  Panel {label} done")
        time.sleep(0.2)

    # 打印信息
    print(f"\nPanel size: {PANEL_W:.1f} x {PANEL_H:.1f} pt "
          f"({PANEL_W / 72 * 2.54:.1f} x {PANEL_H / 72 * 2.54:.1f} cm)")

    # 保存
    if SAVE_PATH:
        options = win32com.client.Dispatch("Illustrator.IllustratorSaveOptions")
        doc.SaveAs(SAVE_PATH, options)
        print(f"Saved: {SAVE_PATH}")


if __name__ == "__main__":
    main()
