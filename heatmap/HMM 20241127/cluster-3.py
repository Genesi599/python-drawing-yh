import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator
from seaborn.external.kde import gaussian_kde


# 读取文件
sheet_names = ['Sheet1', 'Sheet2', 'Sheet3']
namelist = ['Young', 'Old']
df = pd.read_excel('1127-4.xlsx', index_col=0, sheet_name='Sheet3')

# 去掉每个元素中的逗号，然后转换为数值类型
# df = df.replace(',', '', regex=True).astype(float)

# 读取分类数据
# class_df = pd.read_excel('classify.xlsx')
# classname = 'O vs Y'
# df = df[df.index.isin(class_df[classname])]

# 获取 DataFrame 的列名列表
column_names = df.columns.tolist()
# 打印列名列表
print(column_names)



selected_column_names = ['FPKM_C25-M1-1', 'FPKM_C25-M1-2', 'FPKM_C25-M1-3', 'FPKM_NC-1', 'FPKM_NC-2', 'FPKM_NC-3']
# df = df[selected_column_names]

# 初始化一个空的DataFrame来存储结果
z_score_df = pd.DataFrame(index=df.index, columns=df.columns)

# 对每一行计算Z分数并缩放
for index, row in df.iterrows():
    # 计算每一行的均值和标准差
    mean = row.mean()
    std = row.std()

    # 计算Z分数
    z_score = (row - mean) / std

    # 将结果存储到新的DataFrame中
    z_score_df.loc[index] = z_score

# 确保所有数据都是数值类型，并处理缺失值
z_score_df = z_score_df.apply(pd.to_numeric, errors='coerce').fillna(0)

# 转置DataFrame
# z_score_df = z_score_df.transpose()
# 根据每行的均值排序（从大到小）
z_score_df = z_score_df.loc[z_score_df.mean(axis=1).sort_values(ascending=False).index]
print(z_score_df)
# 全局设置字体属性
plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 使用 SimHei 字体以支持中文
plt.rcParams['axes.unicode_minus'] = False  # 确保负号正常显示
# 创建从蓝色到红色的颜色映射
# cmap = LinearSegmentedColormap.from_list('blue_red', ['#517eb9', 'white', '#e97f64'])

# 定义自定义颜色梯度（以 seismic 为基础，可以自定义更多颜色间隔）

colors = [
    (0.0, "blue"),  # -2: 蓝色
    (0.34, "#8c8cff"),  # -1: 浅蓝色
    (0.45, "white"),  # 0: 白色
    (0.52, "#fff4f4"),  # 1: 橙色
    (0.7, "#ff7c7c"),  # 1: 橙色
    # (0.9, "#ff7c7c"),  # 1: 橙色
    (1.0, "red")  # 2: 红色
]

# colors = [
#     (0.0, "blue"),  # -2: 蓝色
#     (0.5, "#fff4f4"),  # 0: 白色
#     (1.0, "red")  # 2: 红色
# ]

# 创建一个自定义的 LinearSegmentedColormap
custom_seismic = LinearSegmentedColormap.from_list("CustomSeismic", colors, N=256)

cmap = custom_seismic
norm = plt.Normalize(vmin=-4, vmax=4)  # 标准化到 z-score 范围

linkage_methods = [
    'single',       # 最短距离法
    'complete',     # 最长距离法
    'average',      # 平均距离法
    'weighted',     # 加权平均法
    'centroid',     # 质心法
    'median',       # 中值法
    'ward'          # 最小方差法（只适用于欧几里得距离）
]

for i in linkage_methods:
    # 绘制热图
    figsize_w = 8
    figsize_h = 8
    cluster = sns.clustermap(data=z_score_df,
                             method=i,
                             figsize=(figsize_w, figsize_h),
                             cbar_pos=None,  # 通过此参数隐藏颜色条
                             cmap=cmap,
                             vmin=-4, vmax=4,
                             xticklabels=True,
                             yticklabels=True,
                             row_colors=None,
                             dendrogram_ratio=(0.1, 0.15),  # 调整树状图和热图比例
                             row_cluster=True,  # 禁用行聚类
                             col_cluster=False
                             )

    # 仅隐藏纵坐标的 tick 标注（刻度标签）
    cluster.ax_heatmap.xaxis.set_visible(False)  # 隐藏 x 轴
    # cluster.ax_heatmap.yaxis.set_visible(False)  # 隐藏 y 轴
    # 获取生成的对象并更新其子图位置
    left, bottom, weight, height = 0.2, 0.2, 0.6, 0.7
    tree_size = 0.1
    tree_gap = 0.01
    cluster.ax_heatmap.set_position([left, bottom, weight, height])  # [左, 下, 宽度, 高度]
    cluster.ax_row_dendrogram.set_position([left - tree_size - tree_gap, 0.2, tree_size, 0.7])  # 左侧树状图位置（手动增加左侧距离）
    cluster.ax_col_dendrogram.set_position(
        [left, bottom + height + tree_gap * figsize_w / figsize_h, weight, tree_size * figsize_w / figsize_h])  # 顶部树状图位置

    # 修改聚类树状图线条的宽度
    for ax in [cluster.ax_row_dendrogram, cluster.ax_col_dendrogram]:
        for line in ax.collections:  # 遍历所有绘图线
            line.set_linewidth(2)  # 设置线条宽度为 2（可以根据需求调整）

    # 设置横坐标刻度标签的大小和旋转角度
    cluster.ax_heatmap.set_xticklabels(
        cluster.ax_heatmap.get_xticklabels(),
        rotation=90, fontsize=15, ha='center', va='top', fontweight='bold'
    )
    cluster.ax_heatmap.tick_params(axis='x', length=0)  # 隐藏横坐标刻度线

    # 设置纵坐标刻度标签的大小和旋转角度
    # 首先确保显示所有刻度
    # cluster.ax_heatmap.yaxis.set_ticks(range(len(cluster.ax_heatmap.get_yticklabels())))

    # 然后设置纵标签的显示属性
    cluster.ax_heatmap.set_yticklabels(
        cluster.ax_heatmap.get_yticklabels(),
        fontsize=11,
        ha='left',
        va='center',
        fontweight='bold'
    )
    cluster.ax_heatmap.tick_params(axis='y', length=0)  # 隐藏横坐标刻度线

    # 绘制颜色条和直方图（叠加效果）
    position1 = [0.2, 0.94, 0.2, 0.05 * figsize_w / figsize_h]  #(left, bottom, width, height)
    # 添加颜色条轴
    cbar_ax = cluster.fig.add_axes(position1)  # 新建轴(left, bottom, width, height)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # 绘制颜色条 (作为背景)
    cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal', ticks=[])
    cbar.ax.xaxis.set_ticks_position('bottom')  # 刻度放到顶部
    cbar.ax.xaxis.set_label_position('top')
    cbar.set_label('Column Z-Score', fontsize=13, fontweight='bold')
    # 调整刻度文字大小
    cbar.ax.set_xticks(np.arange(-3, 4, 3))  # 定义刻度位置
    cbar.ax.tick_params(labelsize=13)  # 设置刻度字体大小为12（你可以根据需要调整数值）
    cbar.ax.tick_params(axis='x', bottom=False, labelbottom=False)  # 隐藏横轴
    # 设置刻度字体为粗体
    bold_font = font_manager.FontProperties(weight='bold')
    for label in cbar.ax.get_yticklabels():
        label.set_fontproperties(bold_font)

        # 新建直方图轴，确保横纵坐标与颜色条对齐
    hist_ax = cluster.fig.add_axes(position1)  # 直方图位置，[left, bottom, width, height]

    # 将热图数据展平为一维数组（适合直方图）
    data = z_score_df.values.flatten()  # 从 DataFrame 中提取数据，展平为一维数组
    data = np.clip(data, -4, 4)  # 限制数据范围在 [-4, 4]，匹配颜色条和直方图范围

    # 修改直方图背景透明度
    hist_ax.set_facecolor((1, 1, 1, 0))  # 将背景设置为完全透明 (RGBA: 白+透明)

    # 绘制直方图
    hist_data, bin_edges = np.histogram(data, bins=500, range=(-4, 4), density=False)  # 使用计数模式 (density=False)
    hist_ax.bar(bin_edges[:-1], hist_data, width=0.16, color='cyan', alpha=1)  # 通过 alpha 参数设置背景透明度

    # 核密度估计 (拟合曲线，用于叠加在直方图上)
    # kde = gaussian_kde(data)
    # x_values = np.linspace(-3, 3, 500)  # 曲线横坐标范围
    # kde_values = kde(x_values) * len(data) * (bin_edges[1] - bin_edges[0])  # 将 KDE 放缩为计数值
    # hist_ax.plot(x_values, kde_values, color='blue', lw=3.0)  # 绘制拟合曲线

    # 调整直方图样式
    hist_ax.set_xlim(-4, 4)  # 横轴范围，与颜色条一致
    hist_ax.set_ylim(0, hist_data.max() * 1.2)  # 为纵轴的计数留出一定的空间
    hist_ax.set_ylabel('Count', fontsize=13, fontweight='bold')  # 将纵轴从密度改为计数
    # hist_ax.tick_params(axis='x', bottom=False, labelbottom=False)  # 隐藏横轴
    # 调整轴外观
    hist_ax.set_xticks(np.arange(-3, 4, 3))  # 定义刻度位置
    hist_ax.tick_params(axis='x', labelsize=16)  # 调整横轴刻度字体大小
    hist_ax.tick_params(axis='y', labelsize=16)  # 调整纵轴刻度字体大小

    # 设置字体为粗体 (通过遍历设置刻度标签属性)
    bold_font = font_manager.FontProperties(weight='bold')  # 定义粗体属性
    for label in hist_ax.get_xticklabels():  # 遍历 X 轴刻度
        label.set_fontproperties(bold_font)
    for label in hist_ax.get_yticklabels():  # 遍历 Y 轴刻度
        label.set_fontproperties(bold_font)

        # 添加竖条颜色框
    height = z_score_df.shape[0]
    weight = z_score_df.shape[1]
    note_number = 2

    # 添加颜色条
    xposition = 0.15
    width1 = weight / note_number - 0.3
    fontsize1 = 22
    ax = cluster.ax_heatmap
    # 添加颜色条
    height1 = 0.2
    for a in range(0, note_number):
        print(a)
        ax.add_patch(Rectangle((xposition + a/note_number*weight, height + 0.5), width1, height1, color='black', transform=ax.transData, clip_on=False))

        # Add text on the left of each color bar
        ax.text(xposition + a/note_number*weight + width1/2, height + 1, namelist[a], va='top', ha='center', fontsize=fontsize1, transform=ax.transData,
                rotation=0, fontweight='bold')


    plt.savefig(f'cluster-{i}.png', format='png', bbox_inches='tight', dpi=600)
