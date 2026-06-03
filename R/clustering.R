calculate_clustering <- function(exprs_matrix, distances, methods){
  hcs <- list()
  for (distance in distances){
    if (distance == "corr"){
      sample_cor_matrix <- cor(exprs_matrix, method = "pearson")
      dist_matrix <- as.dist(1 - sample_cor_matrix)
      for (method in methods){
        hc <- hclust(d = dist_matrix, method = method)
        hcs[[paste0(distance,"_",method)]] <- hc
      }
    } else if (distance == "euclidean"){
      dist_matrix <- dist(t(exprs_matrix))
      for (method in methods){
        hc <- hclust(d = dist_matrix, method = method)
        hcs[[paste0(distance,"_",method)]] <- hc
      }
    } else{
      stop("distance should be corr or euclidean")
    }
  }
  return(hcs)
}

prepare_plot_clustering <- function(model, metadata){
  n <- length(model$order)
  dg <- as.dendrogram(model)
  ddata_pts <- dendro_data(dg, type = "rectangle")
  
  ddata_pts$labels <- ddata_pts$labels %>% mutate(id = label) %>% left_join(metadata)
  
  label_df <- ddata_pts$labels
  segment_df <- segment(ddata_pts)
  
  return(list("label_df" = label_df, "segment_df" = segment_df))
}

plot_clustering <- function(label_df, segment_df, color_var, color_palette = NULL){
  color_var_capital <- set_names(color_var)
  p <- ggplot()+
    geom_segment(data = segment_df, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.2)+
    geom_point(data = label_df, aes(x = x, y = y, fill = variable), alpha = 0.8, shape = 21, color = "black", stroke = 0.2, size = 1.5) +
    theme_classic(base_size = 6)+
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.key.size = unit(0.3, "cm")
    )+
    labs(fill = color_var_capital)
  if (!is.null(color_palette)) {
    p <- p + scale_fill_manual(values = color_palette)
  } else {
    p <- p + scale_fill_viridis_d(option = "plasma",
                                  begin = 0.2,
                                  end = 0.8)
  }
  return(p)
}


process_clustering <- function(exprs_matrix, metadata, all_variables, distances, methods, color_palette = NULL){
  basefilename <- "results/exploratory_analysis/clustering/"
  hcs <- calculate_clustering(exprs_matrix, distances, methods)
  for (i in 1:length(hcs)){
    hc <- hcs[[i]]
    hc_name <- names(hcs)[i]
    hc_name <- unlist(strsplit(hc_name, "_"))
    distance <- hc_name[1]
    method <- hc_name[2]
    clust_elements <- prepare_plot_clustering(hc, metadata)
    label_df <- clust_elements$label_df
    segment_df <- clust_elements$segment_df
    for (var in all_variables){
      label_df_var <- label_df[,c("x", "y", "id", var)]
      colnames(label_df_var)[4] <- "variable"
      p <- plot_clustering(label_df_var, segment_df, var, color_palette[[var]])
      filename <- paste0(basefilename, var, "_", distance, "_", method, ".png")
      ggsave(
        plot = p,
        filename = filename,
        width = 10,
        height = 6,
        units = "cm"
      )
    }
  }
}
