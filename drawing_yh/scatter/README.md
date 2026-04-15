# scatter - 散点图/气泡图绘制模块

## 目录结构

```
scatter/
├── simple_scatter/     # 基础散点图（带显著性标注）
├── linear_fit/         # 散点图 + 线性回归拟合
├── poly_fit/           # 散点图 + 多项式拟合
├── rank_plot/          # Rank abundance plot（排序散点图）
├── bubble/             # 气泡图（上调/下调基因可视化）
├── dot_chart/          # 点图（分组相关性展示）
└── README.md           # 本文件
```

## 模块说明

| 模块 | 功能 |
|------|------|
| **simple_scatter** | 基础散点图，支持按组织着色、显著性标注、相关性分析 |
| **linear_fit** | 散点图 + Seaborn线性回归拟合，显示r/p值，支持多组织子图 |
| **poly_fit** | 散点图 + 多项式拟合，显示R²和置信区间 |
| **rank_plot** | 按丰度排序的散点图，支持差异基因上下调着色与标注 |
| **bubble** | 气泡图，支持上调/下调基因分类、显著性标记、丰度映射 |
| **dot_chart** | 点图，支持分组蛋白/代谢物相关性展示 |
