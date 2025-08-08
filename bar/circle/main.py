import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch
import matplotlib
matplotlib.use('Agg')  # or 'TkAgg'

sheet_name = 'Sheet1'
df = pd.read_excel('data.xlsx', sheet_name=sheet_name)
# 根据Spearman相关性分组
df['group'] = pd.cut(df['spearman_corr'],
                     bins=[-1, 0, 0.3, 0.5, 0.7, 1],
                     labels=['Negative', 'Very Weak', 'Weak', 'Moderate', 'Strong'],
                     include_lowest=True)

# 定义颜色映射
color_map = {'Very Weak': 'green', 'Weak': 'yellow', 'Moderate': 'orange', 'Strong': 'red', 'Negative': 'grey'}
df['color'] = df['group'].map(color_map)

print(df)

# 绘制柱状图
fig, ax = plt.subplots(figsize=(10, 8))
bars = ax.bar(df['ID'], df['spearman_corr'], color=df['color'])

# 添加图例
legend_elements = [Patch(facecolor=color_map[name], label=name) for name in color_map]
ax.legend(handles=legend_elements, title='Correlation Group')

# 添加标题和标签
ax.set_title('Spearman correlation analysis')
ax.set_xlabel('ID')
ax.set_ylabel('Spearman correlation')

# 旋转x轴标签，防止重叠
plt.xticks(rotation=45, ha='right')

# 显示图表
plt.tight_layout()
plt.savefig(f'test.png', dpi=300)