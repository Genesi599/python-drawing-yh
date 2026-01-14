library(ropls)  # 提供PCA, PLS-DA, OPLS-DA功能
library(ggplot2)  # 提供绘图功能

# ==================== 用户输入区（仅需改这里） ============================
csv_file <- "D:/Projects/Bone_Marrow_Aging/proteomics/analysis/abundance_sample_x_protein.csv"
save_dir <- "D:/Projects/Bone_Marrow_Aging/proteomics/analysis/figure"
info_path <- "D:/Projects/Bone_Marrow_Aging/proteomics/analysis/sample_info.csv"

# ---- 读分组信息 ----
meta_df <- read.csv(info_path, row.names=1)  # 读取分组信息

# 确保样本名格式一致，去掉空格并替换为点
rownames(meta_df) <- gsub(" ", ".", rownames(meta_df))
x_df <- read.csv(csv_file, row.names=1)  # 读取蛋白质丰度数据
x_idx <- colnames(x_df)
# 保留全部三组
group_ser <- meta_df[match(x_idx, rownames(meta_df)), 'condition']
# 检查缺失 condition 的样本
if (any(is.na(group_ser))) {
    missing_samples <- x_idx[is.na(group_ser)]
    stop(paste('缺失 condition 的样本:', paste(missing_samples, collapse = ", ")))
}

# 保留全部三组
group_ser <- meta_df[match(x_idx, rownames(meta_df)), 'condition']
if (any(is.na(group_ser))) {
    missing_samples <- x_idx[is.na(group_ser)]
    stop(paste('缺失 condition 的样本:', paste(missing_samples, collapse = ", ")))
}

# 确保 group_list 的长度与样本数量一致
group_list <- as.character(group_ser)

# 三组配色
color_dict <- c(Young = "#1f77b4", Middle = "#ff7f0e", Old = "#d62728")
# ======================== 结束输入区 ========================================

# ==================== PCA分析与可视化 ============================
# 数据标准化
x_scaled <- scale(t(x_df))  # 转置操作使行变为样本

# PCA
pca <- opls(x_scaled, predI = 2)  # 使用ropls，取两个主成分
scores <- pca@scoreMN  # 提取样本坐标（得分）

score_df <- data.frame(
  Sample = rownames(scores),
  PC1 = scores[, 1],
  PC2 = scores[, 2],
  Group = factor(group_list, levels = unique(group_list))
)

# 创建绘图
pca_plot <- ggplot(score_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = Group)) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = color_dict) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Score (PCA)",
    x = paste0("PC 1 (", round(pca@modelDF["R2X"][1,] * 100, 1), "%)"),
    y = paste0("PC 2 (", round(pca@modelDF["R2X"][2,] * 100, 1), "%)")
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"))

# 保存和显示图形
ggsave(paste0(save_dir, "/pca_plot_with_ellipses.png"), plot = pca_plot, height = 5.2, width = 5)
print(pca_plot)