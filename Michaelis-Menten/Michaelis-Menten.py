import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import sem
import pandas as pd
from pathlib import Path

# 读取Excel文件
file_path = Path(r"tolbutamide.xlsx")  # 替换为你的Excel文件路径
df = pd.read_excel(file_path)
print(df)

# 将数据转换为字典，其中键是列名，值是对应的NumPy数组
data_dict = {}

# 假设第一列是底物浓度，其余列是不同组的数据
S = df.iloc[:, 0].values  # 底物浓度
data_dict['S'] = S

# 从第二列开始，每三列一组，分别对应WT, CYP2C9*3, N218A, T229A
for i in range(1, len(df.columns), 3):
    group_name = df.columns[i]  # 组名
    data_dict[group_name] = df.iloc[:, i:i+3].values.T


# 定义Michaelis-Menten方程
def michaelis_menten(S, Vmax, Km):
    return Vmax * S / (Km + S)


# S = np.array([10, 25, 50, 100, 250, 500, 1000, 2000])
v1 = data_dict['WT'][0]
v2 = data_dict['WT'][1]
v3 = data_dict['WT'][2]

def compute_clint(S, v1, v2, v3):

    # 计算每组数据的平均值和标准误
    v_avg = np.mean([v1, v2, v3], axis=0)
    v_sem = sem([v1, v2, v3], axis=0)

    # 拟合每组数据并计算参数的平均值和标准误
    params = []
    for v in [v1, v2, v3]:
        p, p_cov = curve_fit(michaelis_menten, S, v, p0=[1, 1])
        params.append(p)

    # 计算Vmax和Km的平均值和标准误
    Vmax_avg = np.mean([p[0] for p in params])
    Vmax_sem = sem([p[0] for p in params])

    Km_avg = np.mean([p[1] for p in params])
    Km_sem = sem([p[1] for p in params])

    # 计算CLint的平均值和标准误
    CLint_avg = Vmax_avg / Km_avg * 1000
    CLint_sem = CLint_avg * np.sqrt((Vmax_sem / Vmax_avg)**2 + (Km_sem / Km_avg)**2)

    # 将结果存储在字典中
    results = {
        'Vmax_avg': Vmax_avg,
        'Vmax_sem': Vmax_sem,
        'Km_avg': Km_avg,
        'Km_sem': Km_sem,
        'CLint_avg': CLint_avg,
        'CLint_sem': CLint_sem
    }
    return results


WT_clint = compute_clint(S, data_dict['WT'][0], data_dict['WT'][1], data_dict['WT'][2])['CLint_avg']
print(WT_clint)

# 打印转换后的数据
for key, value in data_dict.items():
    print(f"{key}:")
    result = compute_clint(S, value[0], value[1], value[2])
    print(f"{result['Vmax_avg']:.2f}±{result['Vmax_sem']:.3f}")
    print(f"{result['Km_avg']:.2f}±{result['Km_sem']:.2f}")
    print(f"{result['CLint_avg']:.2f}±{result['CLint_sem']:.3f}")
    print(f"{result['CLint_avg']/WT_clint:.2%}")

# 创建一个空的DataFrame来存储结果
results_df = pd.DataFrame(columns=['Group', 'Vmax', 'Km', 'CLint', 'Relative clearance'])

# 计算每个group的结果并添加到DataFrame
for key, value in data_dict.items():
    result = compute_clint(S, value[0], value[1], value[2])
    row = {
        'Group': key,
        'Vmax': f"{result['Vmax_avg']:.2f}±{result['Vmax_sem']:.3f}",
        'Km': f"{result['Km_avg']:.2f}±{result['Km_sem']:.2f}",
        'CLint': f"{result['CLint_avg']:.2f}±{result['CLint_sem']:.3f}",
        'Relative clearance': f"{result['CLint_avg']/WT_clint:.2%}"
    }
    results_df = results_df._append(row, ignore_index=True)

# 将结果保存为Excel文件
results_df.to_excel('results.xlsx', index=False)
