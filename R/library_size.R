plot_library_sizes <- function(df, size_var, var, color_palette = NULL) {
  
  size_var_name <- set_names(size_var)
  color_var_name <- set_names(var)
  
  group_order <- df %>%
    group_by(variable) %>%
    summarise(mean = mean(size_var)) %>%
    arrange(desc(mean)) %>%
    mutate(variable = factor(variable, levels = unique(variable))) %>%
    pull(variable)
  
  df <- df %>% mutate(variable = factor(variable, levels = group_order))
  
  p <- ggplot(df, aes(x = size_var, y = variable))+
    geom_boxplot(outlier.shape = NA, aes(fill = variable), alpha = 0.6, linewidth = 0.3)+
    geom_jitter(height = 0.3, aes(fill = variable), shape = 21, color = "black", stroke = 0.2, size = 1.5) + # modify
    theme_minimal()+
    labs(y = color_var_name, 
         x = size_var_name)+
    theme_minimal()+
    theme(
      axis.text = element_text(size = 5),
      axis.title = element_text(size = 5.5),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks.x = element_line(color = "black", linewidth = 0.3)
    )+
    scale_x_continuous(labels = scales::number_format(scale = 1e-6, suffix = "M"))
  
  if (!is.null(color_palette)){
    p <- p + scale_fill_manual(values = color_palette)
  }
  return(p)
}


library_sizes <- function(counts, dge, metadata, all_variables, color_palette = NULL) {
  library_sizes_df <- tibble(
    id = colnames(counts),
    raw_library_size = colSums(counts),
    filtered_library_size = dge$samples$lib.size
  ) %>%
    left_join(metadata)

  dirname_files <- "results/exploratory_analysis/library_sizes/"
  for (var in all_variables) {
    for (size_var in c("raw_library_size", "filtered_library_size")){
      basename_file <- paste0(var, "-", size_var, ".png")
      df <- library_sizes_df[,c("id", size_var, var)]
      colnames(df) <- c("id", "size_var", "variable")
      p <- plot_library_sizes(df, size_var, var, color_palette = color_palette[[var]])
      filename <- paste0(dirname_files, basename_file)
      ggsave(filename = filename, plot = p,
             height = 7, width = 10, units = "cm")
    }
  }
}