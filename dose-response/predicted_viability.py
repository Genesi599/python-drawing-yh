import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
from scipy.stats import ttest_rel



# 创建一个字典来存储每个 CSV 文件的 DataFrame
def predect_viability(file_path, drug_names):

    total_df = pd.read_csv(file_path)
    def calculate_mean_sem(lst):
        if not lst:  # 如果列表为空，返回None
            return None, None
        mean = np.mean(lst)
        sem = np.std(lst) / np.sqrt(len(lst))
        return mean, sem
    def four_param_logistic(x, bottom, top, ic50, hillslope):
        return bottom + (top - bottom) / (1 + (ic50 / x) ** hillslope)

    # 选择 drug1 数据

    dBET57_df = total_df[total_df[drug_names[1]] == 0]

    result_df = dBET57_df.groupby(drug_names[0])['viability'].apply(list).reset_index()
    # 使用 apply 方法计算每个列表的平均值和标准误
    result_df[['mean_column', 'sem_column']] = result_df['viability'].apply(calculate_mean_sem).apply(pd.Series)
    pd.set_option('display.max_columns', None)  # Display all columns
    # 去除X为0的值

    result_df = result_df[result_df[drug_names[0]] != 0]
    x_data = result_df[drug_names[0]].values
    y_data = result_df['mean_column'].values
    # 拟合曲线
    a = curve_fit(four_param_logistic, x_data, y_data, p0=[0, 1, 300, 1])
    popt = a[0]
    pcov = a[1]
    # # 添加拟合参数信息
    param_names = ['Bottom', 'Top', 'IC50', 'HillSlope']
    # 打印拟合参数及其95%置信区间
    # print("\n拟合参数值：")
    # for name, value, var in zip(param_names, popt, np.diag(pcov)):
    #     print(f"{name}: {value:.3f} ± {np.sqrt(var):.3f}")
    # 生成拟合参数字典
    dBET57_fit_params = {name: {'value': value, 'stderr': np.sqrt(var)} for name, value, var in
                  zip(param_names, popt, np.diag(pcov))}

    total_df['predicted_viability'] = np.nan


    def viability(dBET57_conc, fit_params):
        predicted_viability = four_param_logistic(dBET57_conc, fit_params['Bottom']['value'],
                                                  fit_params['Top']['value'],
                                                  fit_params['IC50']['value'],
                                                  fit_params['HillSlope']['value'])
        return predicted_viability


    # 选择 drug2 数据
    dBET57_df = total_df[total_df[drug_names[0]] == 0]
    result_df = dBET57_df.groupby(drug_names[1])['viability'].apply(list).reset_index()
    # 使用 apply 方法计算每个列表的平均值和标准误
    result_df[['mean_column', 'sem_column']] = result_df['viability'].apply(calculate_mean_sem).apply(pd.Series)
    pd.set_option('display.max_columns', None)  # Display all columns
    # 去除X为0的值
    print(result_df)
    result_df = result_df[result_df[drug_names[1]] != 0]
    x_data = result_df[drug_names[1]].values
    print(x_data)
    y_data = result_df['mean_column'].values
    print(y_data)
    # 拟合曲线
    a = curve_fit(four_param_logistic, x_data, y_data,
                  p0=[0, 1, 300, 1],
                  maxfev=5000)
    popt = a[0]
    pcov = a[1]
    # # 添加拟合参数信息
    # param_names = ['Bottom', 'Top', 'IC50', 'HillSlope']
    # 打印拟合参数及其95%置信区间
    # print("\n拟合参数值：")
    # for name, value, var in zip(param_names, popt, np.diag(pcov)):
    #     print(f"{name}: {value:.3f} ± {np.sqrt(var):.3f}")
    # 生成拟合参数字典
    drug2_fit_params = {name: {'value': value, 'stderr': np.sqrt(var)} for name, value, var in
                  zip(param_names, popt, np.diag(pcov))}
    a = viability(1000, drug2_fit_params)
    b = viability(1000, dBET57_fit_params)

    predicted = 1 - ((1-a) + (1-b) - (1-a)*(1-b))

    for i in range(len(total_df)):
        drug1 = total_df.loc[i, drug_names[0]]
        drug2 = total_df.loc[i, drug_names[1]]
        a = viability(drug2, drug2_fit_params)
        b = viability(drug1, dBET57_fit_params)
        predicted = 1 - ((1 - a) + (1 - b) - (1 - a) * (1 - b))
        total_df.loc[i, 'predicted_viability'] = predicted
    total_df['Difference'] = total_df['viability'] - total_df['predicted_viability']
    # 保存为 CSV 文件
    total_df.to_csv(f'predicted_viability.csv', index=False)


    # column1 = total_df['viability']
    # column2 = total_df['predicted_viability']
    # # 进行配对样本t检验
    # t_stat, p_value = ttest_rel(column1, column2)
    # # 输出结果
    # print(f't-statistic: {t_stat}')
    # print(f'p-value: {p_value}')