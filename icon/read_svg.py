import xml.etree.ElementTree as ET


def change_svg_color(input_path, output_path, color="#FF0000"):
    """读取SVG，修改颜色，保存新SVG"""
    with open(input_path, 'r', encoding='utf-8') as f:
        svg_content = f.read()

    # 解析SVG
    root = ET.fromstring(svg_content)

    # 修改所有path的fill颜色
    for path in root.iter('{http://www.w3.org/2000/svg}path'):
        path.set('fill', color)
        # 移除可能冲突的style
        if 'style' in path.attrib:
            del path.attrib['style']

    # 保存
    tree = ET.ElementTree(root)
    tree.write(output_path, encoding='utf-8', xml_declaration=True)
    print(f"已保存: {output_path}")


# 使用示例
svg_path = r"C:\Users\yh109\OneDrive\桌面\mouse.svg"
output_path = r"C:\Users\yh109\OneDrive\桌面\human_red.svg"

change_svg_color(svg_path, output_path, color="#FF0000")