OPLS_DA = function(dataMatrix, group_list, color, savepath='', filename='') {
  oplsda <- opls(t(dataMatrix[,1:12]), group_list[1:12], predI = 1, orthoI = NA)
  VIP <- oplsda@vipVn  # 提取变量重要性

  score_df <- data.frame(
    Sample = rownames(oplsda@scoreMN),
    Component_1 = oplsda@scoreMN[, 1],
    Component_2 = oplsda@orthoScoreMN[, 1],
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
      title = "Score (OPLS-DA)",
      x = paste0("Component 1 (", round(oplsda@modelDF["R2X"][1,] * 100, 1), "%)"),
      y = paste0("Orthogonal Component 1 (", round(oplsda@modelDF["R2Y"][2,] * 100, 1), "%)")
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          legend.position = c(0.98, 0.98),
          legend.justification = c("right", "top"))

  ggsave(paste0(savepath, filename), height = 5.2, width = 5)
  return(VIP)
}

#调用函数
VIP1 = OPLS_DA(dataMatrix[,1:12], group_list[1:12], biocolor[1:2], 'data/3_metabolomics/', 'Figure 1. OPLS-DA_Normal vs. OGD_R.png')
VIP2 = OPLS_DA(dataMatrix[,7:18], group_list[7:18], biocolor[2:3], 'data/3_metabolomics/', 'Figure 1. OPLS-DA_HSYA75 vs. OGD_R.png')