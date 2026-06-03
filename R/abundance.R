plot_relative_abundance <- function(top_mirnas, n){
  p <- ggplot(top_mirnas, aes(x = Percentage, y = miRNA)) +
    geom_col(fill = "indianred3", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = paste0(round(Percentage, 2), "%")), 
              hjust = -0.1, size = 1.5, color = "black") +
    labs(
      title = paste0("Top ",n," Most Abundant miRNAs"),
      x = "% of Total Reads in Dataset",
      y = NULL
    ) +
    theme_classic(base_size = 6)+
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(top_mirnas$Percentage) + 2))
  return(p)
}

plot_boxplot_abundance <- function(df_boxplot, color_var, color_palette = NULL){
  color_var_capital <- set_names(color_var)
  p <- ggplot(df_boxplot, aes(x = variable, y = log2_cpm, fill = variable))+
    facet_wrap(~miRNA)+
    geom_boxplot(outlier.shape = NA, linewidth = 0.1, alpha = 0.8)+
    theme_classic(base_size = 5)+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )+
    labs(title = paste0("Expression Distribution of Top miRNAs by ", color_var_capital), x = NULL, y = "Expression (Log2 CPM)")
  if (!is.null(color_palette)) {
    p <- p + scale_fill_manual(values = color_palette)
  } else {
    p <- p + scale_fill_viridis_d(option = "plasma",
                                  begin = 0.2,
                                  end = 0.8)
  }
  return(p)
}


process_abundance <- function(counts_matrix, log2_cpm, metadata, all_variables, color_palette = NULL, n = 20){
  total_reads <- sum(counts_matrix)
  top_mirnas <- as.data.frame(counts_matrix) %>%
    rownames_to_column("miRNA") %>%
    pivot_longer(-miRNA, names_to = "id", values_to = "counts") %>%
    group_by(miRNA) %>%
    summarise(Total_Counts = sum(counts)) %>%
    mutate(Percentage = (Total_Counts / total_reads) * 100) %>%
    arrange(desc(Percentage)) %>%
    slice_head(n = n)
  top_mirnas$miRNA <- factor(top_mirnas$miRNA, levels = rev(top_mirnas$miRNA))
  p <- plot_relative_abundance(top_mirnas, n)
  ggsave(filename = "results/exploratory_analysis/abundance/relative_abundance.png", plot = p,
         height = 7, width = 10, units = "cm")
  
  top_names <- as.character(top_mirnas$miRNA)
  df_boxplot <- as.data.frame(log2_cpm[top_names, ]) %>%
    rownames_to_column("miRNA") %>%
    pivot_longer(-miRNA, names_to = "id", values_to = "log2_cpm") %>%
    left_join(metadata, by = "id")
  
  for (var in all_variables){
    df_boxplot_var <- df_boxplot[,c("miRNA", "log2_cpm", var)]
    colnames(df_boxplot_var)[3] <- "variable"
    p <- plot_boxplot_abundance(df_boxplot_var, var, color_palette[[var]])
    ggsave(filename = paste0("results/exploratory_analysis/abundance/",var,"_abundance.png"), plot = p,
           height = 7, width = 10, units = "cm")
  }
}
