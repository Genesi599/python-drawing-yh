import matplotlib.pyplot as plt
import numpy as np

# 创建数据
x = np.linspace(0, 10, 100)
y1 = np.sin(x)
y2 = np.cos(x)

# 创建一个包含两个子图的画布
fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharey=True)

# 在第一个子图中绘制数据
axs[0].plot(x, y1, label='sin(x)')
axs[0].set_title('Subplot 1')
axs[0].legend()

# 在第二个子图中绘制相同的数据
axs[1].plot(x, y2, label='cos(x)')
axs[1].set_title('Subplot 2')
axs[1].legend()

# 在主图中添加共享的y轴标签
fig.text(0.5, 0.04, 'Shared Y-Axis', ha='center', va='center', fontsize=14)

# 显示图形
plt.show()