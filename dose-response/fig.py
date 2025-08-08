import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import pandas as pd

plt.rcParams['font.sans-serif'] = ['Arial']  # 设置全局字体为SimHei
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# 指定文件夹路径
folder_path = 'data/viability'

# 列出文件夹中的所有文件
files = os.listdir(folder_path)

# 过滤出所有以 .csv 结尾的文件
csv_files = [file for file in files if file.endswith('.csv')]

# 创建一个字典来存储每个 CSV 文件的 DataFrame
dfs = {}

# 读取每个 CSV 文件并存储为 DataFrame
for csv_file in csv_files:
    file_path = os.path.join(folder_path, csv_file)
    dfs[csv_file] = pd.read_csv(file_path)

total_df = pd.concat(dfs, ignore_index=True)

dBET57_df = total_df[total_df['CC-90009'] == 0]

result_df = dBET57_df.groupby('dBET1')['viability'].apply(list).reset_index()


# 定义计算平均值和标准误的函数
def calculate_mean_sem(lst):
    if not lst:  # 如果列表为空，返回None
        return None, None
    mean = np.mean(lst)
    sem = np.std(lst) / np.sqrt(len(lst))
    return mean, sem



# 使用 apply 方法计算每个列表的平均值和标准误
result_df[['mean_column', 'sem_column']] = result_df['viability'].apply(calculate_mean_sem).apply(pd.Series)

pd.set_option('display.max_columns', None)  # Display all columns
# 查看结果

result_df = result_df[result_df['dBET1'] != 0]
x_data = result_df['dBET1'].values
y_data = result_df['mean_column'].values


# 定义4PL模型函数
def four_param_logistic(x, bottom, top, ic50, hillslope):
    return bottom + (top - bottom) / (1 + (ic50 / x) ** hillslope)


# 拟合曲线
popt, pcov = curve_fit(four_param_logistic, x_data, y_data, p0=[0, 1, 300, 1])

# # 添加拟合参数信息
param_names = ['Bottom', 'Top', 'IC50', 'HillSlope']

# 打印拟合参数及其95%置信区间
print("\n拟合参数值：")
for name, value, var in zip(param_names, popt, np.diag(pcov)):
    print(f"{name}: {value:.3f} ± {np.sqrt(var):.3f}")
# 生成拟合参数字典
fit_params = {name: {'value': value, 'stderr': np.sqrt(var)} for name, value, var in
              zip(param_names, popt, np.diag(pcov))}




# 绘图
plt.figure(
    figsize=(12, 8)
)

custom_legend = []

# 计算均值和标准误（SE）

# 生成一个对数等间距的数组S_fine，用于后续的精细处理或绘图
# 使用np.logspace函数在x_data的最小值和最大值之间生成100个对数等间距的点
# 这里的对数是基于10的对数，以适应广泛的数据范围并提供更细致的数据处理或绘制需求
print(x_data)

S_fine = np.logspace(np.log10(x_data.min()), np.log10(x_data.max()), 100)
v_fine = four_param_logistic(S_fine, popt[0], popt[1], popt[2], popt[3])


S_fine_max = S_fine.max()
V_fine_max = v_fine.max()
print(S_fine)

color = 'red'
shape = 'o'

# 绘制实验数据点和误差条
plt.errorbar(x_data, result_df['mean_column'].values,
             yerr=result_df['sem_column'].values,
             fmt=shape, capsize=3, markersize=9, elinewidth=3,
             color=color)

# 绘制 Michaelis-Menten 拟合曲线
plt.plot(S_fine, v_fine, color=color, linewidth=3)
# 定制图例，将点形状与曲线样式结合
custom_legend.append(
    Line2D([0], [0], color=color, marker=shape, linestyle='-', markersize=8, linewidth=2)
)

# # 获取 X 轴和 Y 轴的范围，并设置从零开始
# plt.xlim(0, S_fine_max * 1.1)  # 横轴 （从 0 到 `S_fine` 的最大值）
# plt.ylim(0, V_fine_max * 1.2)  # 纵轴 （从 0 到 `v_fine` 的最大值）

# 图表格式设置
plt.xlabel('Substrate Concentration [S]')
plt.ylabel('Reaction Rate v')
plt.legend(handles=custom_legend, fontsize=25, loc="upper left",
           bbox_to_anchor=(1, 1), frameon=False
           )

# 获取主刻度定位器
ax = plt.gca()  # 获取当前坐标轴对象
major_locator = ax.yaxis.get_major_locator()  # 获取 y 轴主刻度定位器
major_ticks = ax.yaxis.get_majorticklocs()  # 获取主刻度值

# 设置横坐标刻度间隔
plt.xscale('log')
# ax.xaxis.set_major_locator(MultipleLocator(S_fine_max / 5))  # 设置主刻度间隔
# 设置纵坐标刻度间隔
# ax.yaxis.set_major_locator(MultipleLocator(V_fine_max)  # 设置主刻度间隔

# 设置刻度线的粗细 (x 轴和 y 轴)
ax.tick_params(axis='both', width=2, length=6)  # 刻度线宽度为 2，长度为 6
ax.spines['bottom'].set_linewidth(2)  # 下边框线
ax.spines['left'].set_linewidth(2)  # 左边框线

# 隐藏右边和上边框线
ax = plt.gca()  # 获取当前的坐标轴对象
ax.spines['right'].set_visible(False)  # 隐藏右边框
ax.spines['top'].set_visible(False)  # 隐藏上边框
ax.spines['left'].set_bounds(0, major_ticks[-2])  # 手动设置 Y 轴框线范围
# ax.spines['bottom'].set_bounds(0, S_fine_max)  # 手动设置 X 轴框线范围

# 设置刻度的字体大小和字体格式
ax.tick_params(axis='x', labelsize=25)  # 横坐标刻度
ax.tick_params(axis='y', labelsize=25)  # 纵坐标刻度
# 设置横轴标签和纵轴标签的字体大小和格式
# plt.xlabel(f"{drug.capitalize()} (μM)", fontsize=28, fontweight="bold")
# plt.ylabel("V (pmol/min/pmol CYP2C9)", fontsize=28, fontweight="bold")
# # 设置标题及字体样式
# plt.title(group, fontsize=26, fontweight="bold")

# 保存图像

# plt.show()
plt.savefig('fig.png', dpi=300, bbox_inches='tight')
plt.savefig('fig.svg', format='svg', bbox_inches='tight')
