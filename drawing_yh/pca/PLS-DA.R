PLS_DA = function(dataMatrix, group_list, color, savepath='', filename='') {
  plsda <- opls(t(dataMatrix), group_list)  # 拟合 PLS-DA 模型
  scores <- plsda@scoreMN

  score_df <- data.frame(
    Sample = rownames(scores),
    Component_1 = scores[, 1],
    Component_2 = scores[, 2],
    Group = factor(group_list, levels = unique(group_list))
  )

  # 可视化
  ggplot(score_df, aes(x = Component_1, y = Component_2, color = Group)) +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(aes(group = Group)) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = color) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Score (PLS-DA)",
      x = paste0("Component 1 (", round(plsda@modelDF["R2X"][1,] * 100, 1), "%)"),
      y = paste0("Component 2 (", round(plsda@modelDF["R2X"][2,] * 100, 1), "%)")
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          legend.position = c(0.98, 0.98),
          legend.justification = c("right", "top"))

  ggsave(paste0(savepath, filename), height = 5.2, width = 5)
}

#调用函数
PLS_DA(dataMatrix, group_list, biocolor, 'data/3_metabolomics/', 'Figure 1. PLS-DA_All_Group.png')
PLS_DA(dataMatrix[,1:12], group_list[1:12], biocolor[1:2], 'data/3_metabolomics/', 'Figure 1. PLS-DA_Normal vs. OGD_R.png')
PLS_DA(dataMatrix[,7:18], group_list[7:18], biocolor[2:3], 'data/3_metabolomics/', 'Figure 1. PLS-DA_HSYA75 vs. OGD_R.png')