import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math
from adjustText import adjust_text

# 导入数据
sheet_names = ['y vs m', 'm vs o', 'o vs y']

for sheet_name in sheet_names:
    pd_df = pd.read_excel('1127-2.xlsx', sheet_name=sheet_name)
    result = pd.DataFrame()

    xname = 'fc'
    yname = 'p'
    IDname = 'proteinID'

    posx = pd_df[pd_df[xname] > 0][xname]
    negx = pd_df[pd_df[xname] < 0][xname]
    x = pd.concat([posx, negx], axis=0)
    result['x'] = x
    result['x'] = np.log2(result['x'])

    posy = pd_df[pd_df[yname] > 0][yname]
    negy = pd_df[pd_df[yname] < 0][yname]
    y = pd.concat([posy, negy], axis=0)
    result['y'] = y
    result['y'] = -np.log10(result['y'])

    # 设置显著性阈值
    x_threshold = 0
    y_threshold = -np.log10(0.1)

    # 添加列
    result = result.sort_index()
    result = pd.concat([result, pd_df[IDname]], axis=1)
    result['average'] = (result['x'] ** 2 + result['y'] ** 2) ** 0.5
    result = result.sort_values(by='average', ascending=False)

    # 分组为up, normal, down
    specific_names = ["ICAM-1", "CXCL16"]  # Replace with the names you want to annotate specially

    result['group'] = 'black'
    result.loc[(result.x > x_threshold) & (result.y > y_threshold), 'group'] = 'tab:pink'  # x=-+x_threshold直接截断
    result.loc[result[IDname].isin(specific_names), 'group'] = 'tab:red'
    result.loc[(result.x < -x_threshold) & (result.y > y_threshold), 'group'] = 'tab:blue'  # x=-+x_threshold直接截断
    result.loc[result.y < y_threshold, 'group'] = 'dimgrey'  # 阈值以下点为灰色





    # 确定坐标轴显示范围
    xmin = math.floor(result['x'].min())
    xmax = math.ceil(result['x'].max())
    ymin = 0
    ymax = math.ceil(result['y'].max())

    # 全局设置字体属性
    plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 使用 SimHei 字体以支持中文

    # 绘制散点图
    fig = plt.figure(figsize=plt.figaspect(9 / 16))  # 确定fig比例（h/w）
    ax = fig.add_subplot()
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')
    ax.scatter(result['x'], result['y'], s=20, c=result['group'])
    ax.set_ylabel('-Log10(P.Value)', fontweight='bold', fontsize=16)
    ax.set_xlabel('Log2 Foldchange', fontweight='bold', fontsize=16)
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
            fontsize=12 if result.iloc[i][IDname] in specific_names else 12,  # Larger fontsize for special names
            color="black" if result.iloc[i][IDname] in specific_names else "black",
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
    ax.vlines(0, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
    ax.hlines(y_threshold, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖水平线

    ax.set_xticks(range(xmin, xmax, 1))  # 设置x轴刻度起点和步长
    ax.set_yticks(range(ymin, ymax, 1))  # 设置y轴刻度起点和步长

    fig.savefig(f'{sheet_name}.png', dpi=300)  # 保存为eps矢量图
