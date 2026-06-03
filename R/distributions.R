plot_density <- function(log_2_cpm_df_var, var, color_palette = NULL){
  var_name <- set_names(var)
  
  p <- ggplot(log2_cpm_df_var, aes(x = log2_cpm))+
    geom_density(aes(color = variable, group = id), linewidth = 0.2, show.legend=FALSE)+
    stat_density(aes(x=log2_cpm, color=variable),
                 geom="line",position="identity")+
    theme_classic(base_size = 7)+
    labs(y = "Density", color = var_name, x = "Log2 CPM", title = paste0("Expression distribution by ", var_name))
  
  if (!is.null(color_palette)) {
    p <- p + scale_color_manual(values = color_palette)
  } else {
    p <- p + scale_color_viridis_d(option = "plasma",
                                  begin = 0.2,
                                  end = 0.8)
  }
  return(p)
}


get_distribution <- function(log2_cpm,
                             metadata,
                             all_variables,
                             color_palette = NULL) {
  log2_cpm_df <- as.data.frame(log2_cpm) %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene,
                 names_to = "id",
                 values_to = "log2_cpm") %>%
    left_join(metadata)
  
  sample_order <- log2_cpm_df %>% group_by(id) %>% summarise(m = mean(log2_cpm)) %>% arrange(desc(m)) %>% pull(id)
  log2_cpm_df <- log2_cpm_df %>% mutate(id = factor(id, levels = sample_order))
  
  for (var in all_variables) {
    log2_cpm_df_var <- log2_cpm_df[,c("id", "log2_cpm", var)]
    colnames(log2_cpm_df_var) <- c("id", "log2_cpm", "variable")
    p <- plot_density(log_2_cpm_df_var, var, color_palette[[var]])
    filename <- paste0("results/exploratory_analysis/distributions/",var,".png")
    ggsave(
      plot = p,
      filename = filename,
      width = 10,
      height = 6,
      units = "cm"
    )
  }
}
  