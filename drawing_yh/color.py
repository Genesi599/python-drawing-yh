import colorsys
import random

from matplotlib import pyplot as plt


def get_n_hls_colors(num, shade='random', random_seed=None):
    hls_colors = []
    step = 360.0 / num  # Ensure hues are evenly distributed
    if shade == 'light':
        l_range = (60, 80)
    elif shade == 'dark':
        l_range = (30, 50)
    else:
        l_range = (40, 80)

    random.seed(random_seed)
    for i in range(num):
        h = (i * step) % 360
        s = 90 + random.random() * 10
        l = l_range[0] + random.random() * (l_range[1] - l_range[0])
        _hlsc = [h / 360.0, l / 100.0, s / 100.0]
        hls_colors.append(_hlsc)

    return hls_colors


def rgba_to_hex(rgba, include_alpha=True):
    """
    将 RGBA 元组转换为十六进制 RGB 或 RGBA 字符串。

    :param rgba: 包含三个整数和一个小数的 RGBA 元组
    :param include_alpha: 是否包含透明度分量
    :return: 十六进制 RGB 或 RGBA 字符串
    """
    r, g, b, a = rgba
    if include_alpha:
        # 将小数透明度转换为整数
        alpha_int = round(a * 255)
        return "#{:02X}{:02X}{:02X}{:02X}".format(r, g, b, alpha_int)
    else:
        return "#{:02X}{:02X}{:02X}".format(r, g, b)

# Visualize color palette and save to file
def visualize_and_save_palette(colors, save_as):
    # Normalize the RGB components

    # Create a figure for the palette
    plt.figure(figsize=(max(len(colors) * 0.5, 8), 2))
    for i in range(len(colors)):
        plt.bar(i, 1, color=colors[i], edgecolor='black')  # Each color as a bar
        plt.text(i, -0.2, str(i), ha='center', va='center', fontsize=8)  # Add index

    # Remove extra axes and save image
    plt.axis('off')
    if save_as:
        plt.savefig(save_as, bbox_inches='tight', dpi=300)

def ncolors(num, shade='random', format='hex', alpha=1.0, save_as='colormap', randon_seed=None):

    if num < 1 or format not in ['rgb', 'hex', 'rgba']:
        return []

    hls_colors = get_n_hls_colors(num, shade, randon_seed)
    colors = []

    for hlsc in hls_colors:
        _r, _g, _b = colorsys.hls_to_rgb(hlsc[0], hlsc[1], hlsc[2])
        r, g, b = [int(x * 255.0) for x in (_r, _g, _b)]

        if format == 'rgb':
            colors.append([r, g, b])
        elif format == 'hex':
            colors.append(rgba_to_hex((r, g, b, alpha)))
        elif format == 'rgba':
            colors.append([r, g, b, alpha])

    # Visualize and save color palette
    # 设置固定的随机种子
    random.seed(randon_seed)
    # 使用 random.shuffle 打乱列表
    random.shuffle(colors)
    visualize_and_save_palette(colors, save_as)
    return colors



