#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@File    : test.py
@Date    : 2026/1/25 16:03
@Author  : yh109
"""

import xml.etree.ElementTree as ET


def change_svg_colors(svg_path, color_map, output_path=None):
    """
    color_map: {'old_color': 'new_color'} 的字典
    例如: {'#FF0000': '#00FF00', 'red': 'blue'}
    """
    tree = ET.parse(svg_path)
    root = tree.getroot()

    # 定义 SVG 命名空间
    ns = {'svg': 'http://www.w3.org/2000/svg'}

    # 遍历所有元素
    for elem in root.iter():
        for attr in ['fill', 'stroke', 'stop-color', 'color']:
            if attr in elem.attrib:
                for old_color, new_color in color_map.items():
                    if elem.attrib[attr] == old_color:
                        elem.attrib[attr] = new_color

    output = output_path or svg_path
    tree.write(output, encoding='utf-8', xml_declaration=True)


    change_svg_colors('input.svg', {'#FF0000': '#0000FF'}, 'output.svg')