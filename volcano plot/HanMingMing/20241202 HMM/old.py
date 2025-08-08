import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math
from adjustText import adjust_text

# 导入数据
sheet_names = ['young', 'middle', 'old']

for sheet_name in sheet_names:
    pd_df = pd.read_excel('20241202.xlsx', sheet_name=sheet_name)
    result = pd.DataFrame()
    xname = 'log2(input)'
    yname = 'log2(IP)'
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
                      # "hsa-miR-1285-3p",
                      "hsa-miR-500b-3P",
                      "hsa-miR-548e-5p"]  # Replace with the names you want to annotate specially

    result['group-edge'] = 'dimgrey'
    result['group-face'] = '#d9f5ff'
    # result.loc[(result.x > x_threshold) & (result.y > y_threshold), 'group'] = 'dimgrey'  # x=-+x_threshold直接截断
    # result.loc[(result.x < -x_threshold) & (result.y > y_threshold), 'group'] = 'dimgrey'  # x=-+x_threshold直接截断
    # result.loc[result.y < y_threshold, 'group'] = 'dimgrey'
    result.loc[result.y > result.x, 'group-face'] = 'pink'
    result.loc[result[IDname].isin(specific_names), 'group-edge'] = 'red'


    # 确定坐标轴显示范围
    xmin = math.floor(result['x'].min())
    xmax = math.ceil(result['x'].max())
    ymin = math.floor(result['y'].min())
    ymax = math.ceil(result['y'].max())

    # 全局设置字体属性
    plt.rcParams['font.sans-serif'] = ['Arial']  # 使用 SimHei 字体以支持中文

    # 绘制散点图
    fig = plt.figure(figsize=plt.figaspect(9 / 16))  # 确定fig比例（h/w）
    ax = fig.add_subplot()
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')

    ax.scatter(result['x'], result['y'], s=50, marker='o',
               facecolors=result['group-face'], edgecolors=result['group-edge'])
    for a in specific_names:
        if a in result[IDname].tolist():
            ax.scatter(
                result[result[IDname] == a]['x'],
                result[result[IDname] == a]['y'],
                c='red',
                s=20,  # 可以增大点的大小
                zorder=10  # 更高的层级
            )

    # 添加 y = x 的虚线
    ax.plot([xmin, xmax], [xmin, xmax], linestyle='--', color='black', linewidth=2, label='y=x')  # 确保起点和终点覆盖范围

    ax.set_ylabel(r'log$_\mathbf{2}$(o⁸G IP)', fontweight='bold', fontsize=18)
    ax.set_xlabel(r'log$_\mathbf{2}$(input)', fontweight='bold', fontsize=18)

    # ax.set_ylabel(r'log₂(o⁸G IP)', fontweight='bold', fontsize=18)
    # ax.set_xlabel(r'log₂(input)', fontweight='bold', fontsize=18)
    # 添加标题
    ax.set_title(sheet_name, fontweight='bold', fontsize=18, color='black', pad=30)  # 设置标题
    # 设置横轴和纵轴刻度文字大小
    ax.tick_params(axis='x', labelsize=14)  # 横轴刻度文字大小设置为12
    ax.tick_params(axis='y', labelsize=14)  # 纵轴刻度文字大小设置为12
    labels = ax.get_xticklabels()
    for label in labels:
        label.set_fontweight('bold')
    labels = ax.get_yticklabels()
    for label in labels:
        label.set_fontweight('bold')

    ax.spines['right'].set_visible(False)  # 去掉右边框
    ax.spines['top'].set_visible(False)  # 去掉上边框

    # 标注点
    # result = result[result['x'] < 0]  # 只取x＜0

    # Create the annotations
    # 筛选符合条件的数据
    filtered_result = result[result['miRNAID'].isin(specific_names)]  # 筛选特定名称
    texts = []

    # 创建标注
    for _, row in filtered_result.iterrows():
        arrow_width = 1
        text_color = '#430000'
        arrow_color = 'black'

        yposition = row['y']+3
        text = ax.text(
            row['x']-2.5,
            yposition,  # 添加y方向的偏移
            row['miRNAID'][4:],
            fontsize=20,  # 标注大小
            color=text_color,
            style="italic",
            weight="bold",
            # verticalalignment='bottom',
            horizontalalignment='center'
        )
        texts.append(text)
        ax.annotate(
            "",  # Empty text because we're only drawing the arrow
            xy=(row['x'], row['y']),  # Point (arrow tip)
            xytext=(row['x']-1.7, yposition),  # Text position
            arrowprops=dict(
                arrowstyle="->",  # Arrow type
                color=arrow_color,     # Arrow color
                lw=arrow_width              # Line width for the arrow
            )
        )


    # adjust_text(texts, arrowprops=dict(arrowstyle='->', lw=0.5, color='green'))

    # 水平和竖直线
    # ax.vlines(-x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
    # ax.vlines(x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
    # ax.vlines(0, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
    # ax.hlines(0, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖水平线

    ax.set_xticks(range(xmin, xmax, 1))  # 设置x轴刻度起点和步长
    ax.set_yticks(range(ymin, ymax, 1))  # 设置y轴刻度起点和步长

    fig.savefig(f'{sheet_name}.png', dpi=300, bbox_inches='tight')  # 保存为eps矢量图
