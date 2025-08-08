import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def sigmoid(x, bottom, top, ic50, hill):
    """
    四参数Logistic sigmoid函数

    参数:
    x: 浓度值
    bottom: 最低响应值
    top: 最高响应值
    ic50: 半数抑制浓度
    hill: Hill斜率系数
    """
    return bottom + (top - bottom) / (1 + 10 ** ((np.log10(ic50) - np.log10(x)) * hill))


def shift_drug_sensitivity_curve(original_x, original_data, new_x, target_ic50):
    """
    平移药敏曲线到新的横坐标系列

    参数:
    original_x: 原始横坐标序列
    original_data: 原始响应值
    new_x: 新的横坐标序列
    target_ic50: 目标IC50值

    返回:
    new_y: 新横坐标对应的响应值
    params: 拟合参数
    """
    # 拟合原始数据
    try:
        popt, *_ = curve_fit(
            sigmoid,
            original_x,
            original_data,
            p0=[min(original_data), max(original_data), np.median(original_x), 1],
            bounds=([0, 0, 0, 0], [np.inf, np.inf, np.inf, 10])
        )
        bottom, top, current_ic50, hill = popt

        # 计算平移后的IC50
        shift_factor = target_ic50 / current_ic50
        new_ic50 = current_ic50 * shift_factor

        # 计算新的y值
        new_y = sigmoid(new_x, bottom, top, new_ic50, hill)

        return new_y, (bottom, top, new_ic50, hill)

    except RuntimeError as e:
        print(f"拟合错误: {e}")
        return None, None


def process_drug_sensitivity_data(original_x, data1, data2, data3, new_x, target_ic50):
    """
    处理三组平行数据

    参数:
    original_x: 原始横坐标序列
    data1, data2, data3: 三组平行数据
    new_x: 新的横坐标序列
    target_ic50: 目标IC50值
    """
    # 数据验证
    if not (len(original_x) == len(data1) == len(data2) == len(data3)):
        raise ValueError("所有输入数据长度必须相同")

    if not all(x > 0 for x in original_x + new_x):
        raise ValueError("浓度值必须为正数")

    if target_ic50 <= 0:
        raise ValueError("目标IC50必须为正数")

    # 处理三组数据
    results = []
    for i, data in enumerate([data1, data2, data3], 1):
        new_y, params = shift_drug_sensitivity_curve(original_x, data, new_x, target_ic50)
        if new_y is not None:
            results.append((new_y, params))
        else:
            raise RuntimeError(f"组 {i} 数据处理失败")

    # 绘图
    plt.figure(figsize=(15, 5))

    # 返回结果
    return [new_y for new_y, _ in results]


# 使用示例
if __name__ == "__main__":

    # 原始横坐标（浓度）
    original_x = [16000, 4000, 1000, 250, 62.5, 15.625, 3.90625, 0.9765625, 0.244140625]



    # 原始数据
    data1 = [6.21758, 12.4096, 79.0351, 97.45105, 106.659, 95.70471, 106.1298, 99.99118, 96.92186]
    data2 = [6.21758, 14.94973, 71.94391, 106.7119, 105.2302, 104.9127, 113.4327, 109.305, 104.5952]
    data3 = [6.555038, 13.09755, 59.56077, 105.0185, 105.389, 94.96384, 117.19, 98.72111, 98.93279]

    # 新横坐标和目标IC50
    new_x = [4000, 1000, 250, 62.5, 15.625, 3.90625, 0.9765625, 0.244140625]
    target_ic50 = 819.5

    try:
        # 处理数据并获取结果
        new_data = process_drug_sensitivity_data(
            original_x, data1, data2, data3,
            new_x, target_ic50
        )

        # 打印结果
        print("\n平移后的数据:")
        print("新横坐标:", new_x)
        for i, data in enumerate(new_data, 1):
            print(f"组{i}:", data.tolist())

    except Exception as e:
        print(f"错误: {e}")
