import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math
from adjustText import adjust_text

# 导入数据
sheet_names = [
    'young',
    'middle',
    'old'
]

for sheet_name in sheet_names:
    pd_df = pd.read_excel('20241202.xlsx', sheet_name=sheet_name)
    result = pd.DataFrame()
    xname = 'log2(IP)'
    yname = 'log2(IP/Input)'
    IDname = 'miRNAID'

    # Use pd.to_numeric with errors='coerce' and drop rows with non-numeric 'col1'
    pd_df = pd_df[pd.to_numeric(pd_df[xname], errors='coerce').notna()]
    pd_df = pd_df[pd.to_numeric(pd_df[yname], errors='coerce').notna()]

    posx = pd_df[pd_df[xname] > 0][xname]
    negx = pd_df[pd_df[xname] < 0][xname]
    x = pd.concat([posx, negx], axis=0)
    result['x'] = x
    # result['x'] = np.log2(result['x'])

    posy = pd_df[pd_df[yname] > 0][yname]
    negy = pd_df[pd_df[yname] < 0][yname]
    y = pd.concat([posy, negy], axis=0)
    result['y'] = y
    # result['y'] = -np.log10(result['y'])

    print(y)

    # 设置显著性阈值
    x_threshold = 0
    y_threshold = -np.log10(0.1)

    # 添加列
    result = result.sort_index()
    result = pd.concat([result, pd_df[IDname]], axis=1)
    result['average'] = (result['x'] ** 2 + result['y'] ** 2) ** 0.5
    result = result.sort_values(by='average', ascending=False)

    # 分组为up, normal, down
    specific_names = ["hsa-miR-98-5p",
                      "hsa-miR-545-3p",
                      "hsa-miR-8085",
                      "hsa-miR-1285-3p",
                      "hsa-miR-500b-3P",
                      "hsa-miR-548e-5p"]  # Replace with the names you want to annotate specially

    result['group'] = 'black'
    result.loc[(result.x > x_threshold) & (result.y > y_threshold), 'group'] = 'dimgrey'  # x=-+x_threshold直接截断
    result.loc[(result.x < -x_threshold) & (result.y > y_threshold), 'group'] = 'dimgrey'  # x=-+x_threshold直接截断
    result.loc[result.y < y_threshold, 'group'] = 'dimgrey'  # 阈值以下点为灰色
    result.loc[result[IDname].isin(specific_names), 'group'] = 'red'

    # 确定坐标轴显示范围
    xmin = math.floor(result['x'].min())
    xmax = math.ceil(result['x'].max())
    ymin = math.floor(result['y'].min())
    ymax = math.ceil(result['y'].max())

    # 全局设置字体属性
    plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 使用 SimHei 字体以支持中文

    # 绘制散点图
    fig = plt.figure(figsize=plt.figaspect(9 / 16))  # 确定fig比例（h/w）
    ax = fig.add_subplot()
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')

    ax.scatter(result['x'], result['y'], s=20, marker='o', facecolors='lightgray', edgecolors=result['group'])
    for a in specific_names:
        if a in result[IDname].tolist():
            ax.scatter(
                result[result[IDname] == a]['x'],
                result[result[IDname] == a]['y'],
                c='red',
                s=20,  # 可以增大点的大小
                zorder=10  # 更高的层级
            )

    ax.set_ylabel('o⁸G enrichment\n[log₂(IP/input)]', fontweight='bold', fontsize=16)
    ax.set_xlabel('log₂(o⁸G IP)', fontweight='bold', fontsize=16)
    # 添加标题
    ax.set_title(sheet_name, fontweight='bold', fontsize=18, color='black', pad=30)  # 设置标题
    # 设置横轴和纵轴刻度文字大小
    ax.tick_params(axis='x', labelsize=14)  # 横轴刻度文字大小设置为12
    ax.tick_params(axis='y', labelsize=12)  # 纵轴刻度文字大小设置为12

    ax.spines['right'].set_visible(False)  # 去掉右边框
    ax.spines['top'].set_visible(False)  # 去掉上边框

    # 标注点
    # result = result[result['x'] < 0]  # 只取x＜0

    # Create the annotations
    texts = [
        ax.text(
            result.iloc[i]['x'],  # Adjusting x-offset for the label
            result.iloc[i]['y'],  # y-coordinate
            result.iloc[i][IDname],  # Text to display (e.g., IDname column value)
            fontsize=18 if result.iloc[i][IDname] in specific_names else 14,  # Larger fontsize for special names
            color="red" if result.iloc[i][IDname] in specific_names else "black",
            # Red for special names, black otherwise
            style="italic" if result.iloc[i][IDname] in specific_names else "normal",  # Italic for both but can change
            weight="bold",
            verticalalignment='top',
            horizontalalignment='left'
        )
        for i in range(len(result))  # Loop over all rows in 'result'
        if result.iloc[i][IDname] in specific_names
    ]

    adjust_text(texts, arrowprops=dict(arrowstyle='->', lw=0.5, color='green'))

    # 水平和竖直线
    # ax.vlines(-x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
    # ax.vlines(x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
    # ax.vlines(0, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
    ax.hlines(0, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖水平线

    ax.set_xticks(range(xmin, xmax, 1))  # 设置x轴刻度起点和步长
    ax.set_yticks(range(ymin, ymax, 1))  # 设置y轴刻度起点和步长

    fig.savefig(f'{sheet_name}.png', dpi=300)  # 保存为eps矢量图
