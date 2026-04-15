import pandas as pd

# 读取 CSV 文件
df = pd.read_csv('data.csv')

# 查看前几行
print(df.head())

# 访问某一列
print(df['列名'])