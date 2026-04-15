import pandas as pd  # Data analysis
import numpy as np  # Scientific computing
import matplotlib.pyplot as plt  # Plotting
import math

# 导入数据
# Specify the sheet you want to read, e.g., 'Sheet1'
sheet_names = ['young  vs middle ', 'young  vs  old', 'middle vs old']
for sheet_name in sheet_names:
    words = sheet_name.split()  # 按空格分割字符串
    last_word = words[-1]  # 获取最后一个单词
    specific_names = [
        'hsa-miR-98-5p',
        'hsa-miR-8085',
        'hsa-miR-1285-3p',
        'hsa-miR-545-3p',
        'hsa-miR-548e-5p',
        'hsa-miR-500b-3p'
    ]
    xlabel_name = r'log$_\mathbf{2}$(o⁸G IP/input)' + f'in {last_word} group'
    # xlabel_name = r'log$_\mathbf{2}$(o⁸G IP/input)' + f'in late middle-aged group'
    # xlabel_name = r'log₂(o⁸G IP/input)' + f'in {last_word} group'
    # xlabel_name = r'log₂(o⁸G IP/input)' + f'in late middle-aged group'

    # Read the specified sheet from the Excel file
    pd_df = pd.read_excel('human samples.xlsx', sheet_name=sheet_name)
    # print(pd_df)

    x_name = 'log2(IP/Input)-2'
    y_name = 'Relative Enrichment 2'
    IDname = 'miRNAID'

    result = pd.DataFrame()

    x = pd_df[x_name]
    result['x'] = x
    y = pd_df[y_name]
    result['y'] = y

    # result['y'] = -np.log10(result['y'])

    # 设置显著性阈值
    x_threshold = 1
    y_threshold = 2

    result = result.dropna()
    result = result[result['x'] >= 0]

    print(result)


    def log_scale(series, base=10):
        """
        保留原始数据正负号的对数缩放

        参数:
        series: 输入的pandas Series
        base: 对数的底数（默认为10）

        返回:
        对数缩放后的Series
        """

        def safe_log(x, base):
            """安全的对数计算"""
            if x == 0:
                return 0
            return np.sign(x) * np.log(abs(x)) / np.log(base)

        return series.apply(lambda x: safe_log(x, base))



    # Z分数标准化
    result['y'] = log_scale(result['y'])

    # 添加列, 计算与原点距离
    result = result.sort_index()
    result = pd.concat([result, pd_df[IDname]], axis=1)
    result['average'] = (result['x'] ** 2 + result['y'] ** 2) ** 0.5
    result = result.sort_values(by='average', ascending=False)


    # 分组为up, normal, down
    result['group-edge'] = 'dimgrey'
    result['group-face'] = '#d9f5ff'
    # result.loc[(result.x > x_threshold) & (result.y > y_threshold), 'group'] = 'dimgrey'  # x=-+x_threshold直接截断
    # result.loc[(result.x < -x_threshold) & (result.y > y_threshold), 'group'] = 'dimgrey'  # x=-+x_threshold直接截断
    # result.loc[result.y < y_threshold, 'group'] = 'dimgrey'
    result.loc[result.y > 0, 'group-face'] = 'pink'
    # result.loc[result.y > 0, 'group-edge'] = 'pink'
    for i in specific_names:
        if i in result['miRNAID'].tolist():
            result.loc[result['miRNAID'] == i, 'group'] = 'tab:red'  # x=-+x_threshold直接截断


    # 确定坐标轴显示范围
    xmin = 0
    xmax = math.ceil(result['x'].max())

    ymax = math.ceil(result['y'].abs().max())+1
    ymin = -ymax

    # 全局设置字体属性
    plt.rcParams['font.sans-serif'] = ['Arial']  # 使用 SimHei 字体以支持中文
    plt.rcParams['axes.unicode_minus'] = False  # 确保负号正常显示

    # 绘制散点图
    fig = plt.figure(figsize=plt.figaspect(9 / 16))  # 确定fig比例（h/w）
    ax = fig.add_subplot()
    ax.set(xlim=(xmin - 0.1, xmax + 0.1), ylim=(ymin, ymax+0.1), title='')
    ax.scatter(result['x'], result['y'], s=35, marker='o',
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
    for a in specific_names:
        if a in result['miRNAID'].tolist():
            ax.scatter(
                result[result['miRNAID'] == a]['x'],
                result[result['miRNAID'] == a]['y'],
                c='red',
                s=20,  # 可以增大点的大小
                zorder=10  # 更高的层级
                )
    ax.set_xlabel(xlabel_name, fontweight='bold', fontsize=20)
    ax.set_ylabel('Relative o⁸G enrichment', fontweight='bold', fontsize=20)
    ax.spines['right'].set_visible(False)  # 去掉右边框
    ax.spines['top'].set_visible(False)  # 去掉上边框


    # 筛选符合条件的数据
    filtered_result = result[result['miRNAID'].isin(specific_names)]  # 筛选特定名称
    texts = []

    # 创建标注
    for _, row in filtered_result.iterrows():
        text = ax.text(
            row['x']+0.6,
            row['y']+1.1,  # 添加y方向的偏移
            row['miRNAID'][4:],
            fontsize=20,  # 标注大小
            color='#430000',
            style="italic",
            weight="bold",
            # verticalalignment='bottom',
            horizontalalignment='center'
        )
        texts.append(text)
        ax.annotate(
            "",  # Empty text because we're only drawing the arrow
            xy=(row['x'], row['y']),  # Point (arrow tip)
            xytext=(row['x']+0.6, row['y']+1.1),  # Text position
            arrowprops=dict(
                arrowstyle="->",  # Arrow type
                color="black",     # Arrow color
                lw=1             # Line width for the arrow
            )
        )

        # Adjust the text positions and add arrows to connect annotations to their scatter points
        # adjust_text(
        #     texts,
        #     ax=ax,
        #     expand_text=(1.2, 1.5),  # Avoid overlapping text
        #     arrowprops=dict(
        #         arrowstyle="->",  # Arrow style
        #         color="gray",      # Arrow color
        #         lw=1               # Line width for the arrow
        #     )
        # )

    # 水平和竖直线
    # ax.vlines(0, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
    # ax.vlines(x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖直线
    ax.hlines(0, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)  # 画竖水平线

    ax.set_xticks(range(xmin, xmax, 1))  # 设置x轴刻度起点和步长
    ax.set_yticks(range(ymin+1, ymax, 1))  # 设置y轴刻度起点和步长
    # 设置刻度字体大小
    plt.xticks(fontsize=20, fontweight='bold')  # x轴刻度字体大小
    plt.yticks(fontsize=20, fontweight='bold')  # y轴刻度字体大小

    fig.savefig(f"{sheet_name}.png", dpi=600, bbox_inches='tight')  # 保存为eps矢量图
