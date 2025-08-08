import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

def get_heapmap(drug_names, value_names):

    # 读取CSV文件
    df = pd.read_csv('predicted_viability.csv', index_col=0)
    # 方法2：同时修改多列
    columns_to_round = [drug_names[0], drug_names[1]]
    df[columns_to_round] = df[columns_to_round].round(2)
    print(df)
    # 使用pivot_table方法重新组织数据
    pivot_df = df.pivot_table(index=drug_names[1], columns=drug_names[0], values=value_names, aggfunc='mean')

    print(pivot_df)

    # 全局设置字体属性
    plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 使用 SimHei 字体以支持中文
    plt.rcParams['axes.unicode_minus'] = False  # 确保负号正常显示

    # 创建从蓝色到红色的颜色映射
    cmap = LinearSegmentedColormap.from_list('blue_red', ['blue', 'white', 'red'])

    # 创建一个新的图形和轴对象
    fig, ax = plt.subplots(dpi=600, figsize=(8, 6))

    # 绘制热图
    sns.heatmap(data=pivot_df,
                cmap=cmap,
                vmin=-1, vmax=1,
                xticklabels=True,
                cbar_kws={'orientation': 'horizontal', 'shrink': 0.2, 'pad': 0.1},
                ax=ax)

    # 设置横坐标刻度标签的大小和旋转角度
    plt.xticks(rotation=45, fontsize=10, ha='right', va='top', fontweight='bold')
    plt.tick_params(axis='x', length=0)  # 隐藏横坐标刻度线

    # 将纵坐标标签放到右边并水平显示
    ax.yaxis.tick_right()               # 将刻度移至右边
    ax.yaxis.set_label_position("right")  # 将标签移至右边
    plt.yticks(rotation=0, fontsize=10, fontweight='bold')
    plt.tick_params(axis='y', length=0)  # 隐藏横坐标刻度线
    # ax.set_yticklabels([])  # 仅隐藏刻度值


    # 将图例放到左上方
    cbar = ax.collections[0].colorbar
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.set_position([0.1, 0.8, 0.2, 0.15])  # [left, bottom, width, height]
    # 设置 Colorbar 刻度值的文字大小
    cbar.ax.tick_params(labelsize=15)  # 设置刻度文字大小为 10
    for label in cbar.ax.get_xticklabels():  # 设置刻度文字加粗
        label.set_fontweight('bold')
    # 添加颜色条注释
    cbar.ax.text(0.5, -1.5, value_names, ha='center', va='center', fontsize=18,  fontweight='bold',transform=cbar.ax.transAxes)


    # 设置横轴标签和纵轴标签的字体大小和格式
    plt.xlabel(f"{drug_names[0]} (nM)", fontsize=18, fontweight="bold")
    plt.ylabel(f"{drug_names[1]}", fontsize=18, fontweight="bold")


    plt.savefig('heatmap_with_colorbar.png', format='png', bbox_inches='tight')