calculate_pca <- function(exprs_matrix, metadata){
  pca_res <- prcomp(t(exprs_matrix), scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x) %>%
    rownames_to_column("id") %>%
    left_join(metadata)
  
  var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
  names(var_explained) <- colnames(pca_res$x)
  return(list("pca_df" = pca_df, "var_explained" = var_explained))
}

plot_scree <- function(var_explained, n = 10){
  var_explained_df <- as.data.frame(var_explained) %>% rownames_to_column("PC") %>% head(n = n) %>% arrange(desc(var_explained)) %>% mutate(PC = factor(PC, levels = PC))
  
  p <- ggplot(var_explained_df, aes(x = PC, y = var_explained))+
    geom_col(fill = "steelblue2")+
    geom_text(aes(label = paste0(round(var_explained, 1), " %"), y = var_explained + 1)) +
    theme_classic(base_size = 6)+
    labs(y = "Percentage of Variance Explained (%)")+
    theme(
      axis.title.x = element_blank()
    )
  return(p)
}

plot_pca <- function(pca_df_plot, x, y, var_explained, color_var, color_palette = NULL){
  pcx_lab <- paste0(x, " ( ", round(var_explained[x], 1), "% )")
  pcy_lab <- paste0(y, " ( ", round(var_explained[y], 1), "% )")
  color_var_capital <- set_names(color_var)
  
  p <- ggplot(pca_df_plot, aes(x = pcx, y  = pcy, fill = variable))+
    theme_classic(base_size = 6)+
    labs(x = pcx_lab, y = pcy_lab, fill = color_var_capital)+
    geom_point(shape = 21, alpha = 0.8, color = "black", size = 2)+
    theme(
      legend.key.size = unit(0.3, "cm")
    )
  
  if (!is.null(color_palette)) {
    p <- p + scale_fill_manual(values = color_palette)
  } else {
    p <- p + scale_fill_viridis_d(option = "plasma",
                                  begin = 0.2,
                                  end = 0.8)
  }
  return(p)
}


process_pca <- function(exprs_matrix, metadata, all_variables, interest_pcs = c("PC1", "PC2"), color_palette = NULL){
  pca_data <- calculate_pca(exprs_matrix, metadata)
  var_explained <- pca_data$var_explained
  pca_df <- pca_data$pca_df
  
  p <- plot_scree(var_explained)
  ggsave(
    plot = p,
    filename = "results/exploratory_analysis/pca/scree.png",
    width = 10,
    height = 6,
    units = "cm"
  )
  
  combs_matrix <- as.data.frame(t(combn(interest_pcs, 2)))
  for (i in 1:nrow(combs_matrix)) {
    x <- combs_matrix$V1[i]
    y <- combs_matrix$V2[i]
    for (var in all_variables){
      pca_df_plot <- pca_df[,c("id", x, y, var)]
      colnames(pca_df_plot)[c(2, 3, 4)] <- c("pcx", "pcy", "variable")
      p <- plot_pca(pca_df_plot, x, y, var_explained, var, color_palette[[var]])
      ggsave(
        plot = p,
        filename = paste0(
          "results/exploratory_analysis/pca/",
          var, "_", x, "_", y,".png"
        ),
        width = 12,
        height = 7,
        units = "cm"
      )
    }
  }
  
}