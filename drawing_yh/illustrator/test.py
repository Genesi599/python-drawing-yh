#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : test.py
@Date    : 2026/4/8 09:52
@Author  : yh109
"""

import win32com.client
import os

PDF_PATH = r"D:\Projects\Bone_Marrow_Aging\proteomics\analysis\figure\protein_patterns_heatmap_downregulation.pdf"  # 修改为实际 PDF 路径

ai = win32com.client.Dispatch("Illustrator.Application")

# 新建 A4 文档 (RGB, 595x842 pt)
doc = ai.Documents.Add(2, 595, 842)
layer = doc.Layers[0]
layer.Name = "pdf_layer"

# 将 PDF 作为链接文件嵌入到文档中
# PlaceFile(FileName, Position, ColorModel, Template, Replace)
placed = doc.PlacedItems.Add()
placed.File = PDF_PATH

# 读取 PDF 原始尺寸（PlaceFile 后自动获得）
orig_w = placed.Width
orig_h = placed.Height

# 若超出页面则等比缩放到页面内，否则保持原尺寸
page_w, page_h = 595, 842
scale = min(page_w / orig_w, page_h / orig_h, 1.0)  # 1.0 保证不放大
new_w = orig_w * scale
new_h = orig_h * scale

# 居中放置
left = (page_w - new_w) / 2
top  = page_h - (page_h - new_h) / 2   # AI 坐标 top = 页面高 - 上边距

placed.Width = new_w
placed.Height = new_h
placed.Position = [left, top]

# 嵌入（Embed）：将链接转为内嵌对象，断开与原文件的链接
placed.Embed()

print(f"原始尺寸: {orig_w:.1f} x {orig_h:.1f} pt")
print(f"缩放比例: {scale:.4f}")
print(f"嵌入后尺寸: {new_w:.1f} x {new_h:.1f} pt")
