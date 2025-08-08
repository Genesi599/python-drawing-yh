import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib
from brokenaxes import brokenaxes

matplotlib.use('Agg')  # or 'TkAgg'

# 生成示例数据
sheet_name = 'Sheet1'
df = pd.read_excel('CLEC10.xlsx', sheet_name=sheet_name)
print(df)
y_name = 'CLEC10A'
df = df.dropna()
print(df)
df = df[df[y_name] != 0]
data = df
print(df)
# 使用多项式拟合
degree = 3
coeffs = np.polyfit(data['age'], data[y_name], degree)
p = np.poly1d(coeffs)
x_fit = np.linspace(data['age'].min(), data['age'].max(), 100)
y_fit = p(x_fit)
y_fit = np.log(y_fit)

# 计算置信区间
y_fit_pred = p(data['age'])
residuals = data[y_name] - y_fit_pred
ss_res = np.sum(residuals ** 2)
ss_tot = np.sum((data[y_name] - data[y_name].mean()) ** 2)
r_squared = 1 - (ss_res / ss_tot)
# 绘制图形


fig = plt.figure(figsize=(8, 6))  # Create a figure
# Create broken axes by specifying the ranges for the y-axis
bax = brokenaxes(ylims=((1, 200), (2000, 10000)), hspace=0.05)
# 设置y轴为对数刻度
bax.set_yscale('log')
# First subplot
bax.scatter(data['age'], data[y_name], s=10, label='Data Points')
bax.plot(x_fit, y_fit, color='red', label=f'R-squared: {r_squared:.3f}')
bax.fill_between(x_fit, y_fit - 0.5 * np.std(residuals), y_fit + 0.5 * np.std(residuals), color='gray', alpha=0.2)

bax.set_ylabel(y_name)
bax.legend()
# Optionally: ax1.set_title(f'{y_name} vs Age with Polynomial Fit', fontsize=20)


# 生成示例数据
sheet_name = 'Sheet2'
df = pd.read_excel('CLEC10.xlsx', sheet_name=sheet_name)
print(df)
y_name = 'IFNb'
df = df.dropna()
print(df)
df = df[df[y_name] != 0]
data = df
print(df)
# 使用多项式拟合
degree = 3
coeffs = np.polyfit(data['age'], data[y_name], degree)
p = np.poly1d(coeffs)
x_fit = np.linspace(data['age'].min(), data['age'].max(), 100)
y_fit = p(x_fit)
y_fit = np.log(y_fit)
print(x_fit)
print(y_fit)
# 计算置信区间
y_fit_pred = p(data['age'])
residuals = data[y_name] - y_fit_pred
ss_res = np.sum(residuals ** 2)
ss_tot = np.sum((data[y_name] - data[y_name].mean()) ** 2)
r_squared = 1 - (ss_res / ss_tot)



# Second subplot
bax.scatter(data['age'], data[y_name], s=10, label='Data Points')
bax.plot(x_fit, y_fit, color='red', label=f'R-squared: {r_squared:.3f}')
bax.fill_between(x_fit, y_fit - 0.5 * np.std(residuals), y_fit + 0.5 * np.std(residuals), color='gray', alpha=0.2)

bax.set_xlabel('Age (year)')
# No need to set ylabel again
bax.legend()
# Optionally: ax2.set_title(f'{y_name} vs Age with Polynomial Fit', fontsize=20)




plt.xlabel('Age', fontsize=20)  # 设置X轴标签的字体大小为12
plt.ylabel(
    f'{y_name}',
    # 'IFNB1',
    fontsize=20)  # 设置Y轴标签的字体大小为12


# 设置横纵坐标数字的字体大小
plt.tick_params(axis='both', which='major', labelsize=15)  # 设置主要刻度标签
print(f'R-squared: {r_squared:.3f}')
plt.savefig(f'test.png', dpi=300, bbox_inches='tight')
