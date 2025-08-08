import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
import pandas as pd
from matplotlib import rcParams
from matplotlib.lines import Line2D
import json

# 从 JSON 文件读取列表数据
with open("increased.json", "r") as file:
    increased = json.load(file)
with open("decreased.json", "r") as file:
    decreased = json.load(file)

# 读取Excel文件
file_path = 'diclofenac.xlsx'  # 替换为你的Excel文件路径
df = pd.read_excel(file_path)

# 将数据转换为字典，其中键是列名，值是对应的NumPy数组
data_dict = {}

# 假设第一列是底物浓度，其余列是不同组的数据
S = df.iloc[:, 0].values  # 底物浓度
data_dict['S'] = S

# 从第二列开始，每三列一组，分别对应WT, CYP2C9*3, N218A, T229A
for i in range(1, len(df.columns), 3):
    group_name = df.columns[i]  # 组名
    data_dict[group_name] = df.iloc[:, i:i + 3].values.T


# 定义 Michaelis-Menten 方程
def michaelis_menten(S, Vmax, Km):
    return (Vmax * S) / (Km + S)


def generate_colors(n):
    # 使用 matplotlib 的内置 colormap (tab10, tab20 等)
    cmap = plt.get_cmap('tab20')  # 'tab10' 最多支持 10 种颜色
    colors = [cmap(i % 10) for i in range(n)]  # 取余以循环颜色
    return colors


def generate_shapes(n):
    # 定义固定的点形状集合
    markers = ['o', 's', '^', 'v', '<', '>', 'd', 'p', '*', 'x', '+', ',', '1', '2', '3', '4', '|']
    # 循环取点形状，直到满足 n 个
    shapes = [markers[i % len(markers)] for i in range(n)]
    return shapes


# 生成固定 8 种颜色，每次调用生成的颜色都一样


rcParams['font.family'] = 'Arial'
rcParams['font.weight'] = 'bold'  # 全局字体为粗体


def merge_duplicates(list1, list2):
    # 找到两个列表的交集部分
    duplicates = list(set(list1) & set(list2))
    return duplicates


# 三次重复实验的反应速率 [v] （假设的实验数据，可以替换为实际数据）

keys = list(data_dict.keys())
keys.remove('S')

df = pd.DataFrame(index=keys)
df["Vmax"] = np.nan
df['Km'] = np.nan
df['v_mean'] = [np.array([]) for _ in range(len(df))]
df['v_se'] = [np.array([]) for _ in range(len(df))]
df['color'] = generate_colors(len(df))
df['shape'] = generate_shapes(len(df))

for i in list(df.index):
    print(i)
    v_repeats = data_dict[i]
    # 计算均值和标准误（SE）
    a= v_repeats.mean(axis=0)
    df.at[i, 'v_mean'] = np.array(a)  # 每行取均值
    b = v_repeats.std(axis=0, ddof=1) / np.sqrt(v_repeats.shape[1])
    print(b)
    df.at[i, 'v_se'] = v_repeats.std(axis=0, ddof=1) / np.sqrt(v_repeats.shape[1])  # 标准误 = 标准差 / sqrt(重复次数)
    # 使用 Michaelis-Menten 方程进行拟合
    initial_guess = [3, 0.5]  # 初始值 [Vmax, Km]
    result = curve_fit(michaelis_menten, S, df.loc[i, 'v_mean'], p0=initial_guess)
    # 从 result 中分别提取参数和协方差矩阵
    params = result[0]  # 第一个值是拟合参数 (params)
    covariance = result[1]  # 第二个值是协方差矩阵 (covariance)
    df.loc[i, 'Vmax'] = params[0]
    df.loc[i, 'Km'] = params[1]

keys = merge_duplicates(keys, increased)
keys.insert(0, 'WT')
keys.insert(1, '*3')

df_filtered = df[df.index.isin(keys)]
sorted_df = df_filtered.sort_values(by='Vmax', ascending=False)
print(sorted_df)


# 绘图
plt.figure(
    figsize=(10, 6)
)
S_fine_max = 0
V_fine_max = 0
custom_legend = []

colors = generate_colors(len(keys))
shapes = generate_shapes(len(keys))

for i in list(sorted_df.index):
    # 计算均值和标准误（SE）
    # 生成拟合曲线的数据点
    S_fine = np.linspace(0, S.max(), 100)
    v_fine = michaelis_menten(S_fine, sorted_df.loc[i, 'Vmax'], sorted_df.loc[i, 'Km'])
    if S_fine.max() > S_fine_max:
        S_fine_max = S_fine.max()
    if v_fine.max() > V_fine_max:
        V_fine_max = v_fine.max()

    color = sorted_df.loc[i, 'color']
    shape = sorted_df.loc[i, 'shape']

    # 绘制实验数据点和误差条
    plt.errorbar(S, sorted_df.loc[i, 'v_mean'],
                 yerr=sorted_df.loc[i, 'v_se'],
                 fmt=shape, capsize=3, markersize=7, elinewidth=3,
                 color = color)
    # 绘制 Michaelis-Menten 拟合曲线
    plt.plot(S_fine, v_fine, color=color, linewidth=2)
    # 定制图例，将点形状与曲线样式结合
    custom_legend.append(
        Line2D([0], [0], color=color, marker=shape, linestyle='-', label=i, markersize=8)
    )

# 获取 X 轴和 Y 轴的范围，并设置从零开始
plt.xlim(0, S_fine_max * 1.1)  # 横轴 （从 0 到 `S_fine` 的最大值）
plt.ylim(0, V_fine_max * 1.2)  # 纵轴 （从 0 到 `v_fine` 的最大值）

# 图表格式设置
plt.xlabel('Substrate Concentration [S]')
plt.ylabel('Reaction Rate v')
plt.legend(handles=custom_legend, fontsize=18, loc="upper left",
           bbox_to_anchor=(1, 1), frameon=False
           )


# 设置横轴标签和纵轴标签的字体大小和格式
plt.xlabel("X-axis (Time)", fontsize=20, fontweight="bold")
plt.ylabel("Y-axis (Amplitude)", fontsize=20, fontweight="bold")

# 获取主刻度定位器
ax = plt.gca()  # 获取当前坐标轴对象
major_locator = ax.yaxis.get_major_locator()  # 获取 y 轴主刻度定位器
major_ticks = ax.yaxis.get_majorticklocs()  # 获取主刻度值

# 设置横坐标刻度间隔
ax.xaxis.set_major_locator(MultipleLocator(S_fine_max / 5))  # 设置主刻度间隔
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
print(major_ticks)
ax.spines['left'].set_bounds(0, major_ticks[-2])  # 手动设置 Y 轴框线范围
ax.spines['bottom'].set_bounds(0, S_fine_max)  # 手动设置 X 轴框线范围

# 设置刻度的字体大小和字体格式
ax.tick_params(axis='x', labelsize=20)  # 横坐标刻度
ax.tick_params(axis='y', labelsize=20)  # 纵坐标刻度

# 设置标题及字体样式
# plt.title("Customized Axis Font Example", fontsize=16, fontfamily="sans-serif", fontweight="bold")

# 保存图像

# plt.show()
plt.savefig('fig.png', dpi=600, bbox_inches='tight')
plt.savefig('fig.svg', format='svg', bbox_inches='tight')