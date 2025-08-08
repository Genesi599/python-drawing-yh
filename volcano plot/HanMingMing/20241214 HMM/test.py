import matplotlib.pyplot as plt
from adjustText import adjust_text

# 数据点
x = [1, 2, 3, 4, 5]
y = [5, 7, 6, 8, 7]

# 要标注的点及对应标注文字
annotations = {
    (2, 7): "Key Point 1",
    (4, 8): "Key Point 2",
    (5, 7): "Key Point 3"
}

# 创建图形
plt.figure(figsize=(8, 6))

# 绘制点
plt.scatter(x, y, color='blue', label='Data points')

# 初始化文本集合
texts = []
for (x_coord, y_coord), text in annotations.items():
    t = plt.text(x_coord, y_coord, text, fontsize=10, color='black')
    texts.append(t)

# 自动调整文字位置以避免重叠
adjust_text(
    texts,
    arrowprops=dict(
        arrowstyle="->",  # 定义箭头样式
        color='red',      # 箭头颜色
        lw=1              # 箭头线宽
    )
)

# 添加图例
plt.legend()

# 添加标题和坐标轴标签
plt.title("Scatter Plot with Auto-adjusted Annotations")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")

# 显示网格和图形
plt.grid(True)
plt.tight_layout()
plt.show()