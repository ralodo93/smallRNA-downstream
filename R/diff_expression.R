plot_body_heatmap <- function(heatmap_df){
  # scale by miRNA forcing a raw numeric vector instead of a matrix object
  heatmap_df <- heatmap_df %>% 
    group_by(miRNA) %>% 
    mutate(scaled = as.numeric(scale(expression))) %>% 
    ungroup()
  
  # Calculate symmetric limits for Z-score visualization
  max_abs <- max(abs(heatmap_df$scaled), na.rm = TRUE)
  
  p <- ggplot(heatmap_df, aes(x = id, y = miRNA, fill = scaled))+
    geom_tile()+
    scale_fill_gradient2(limits = c(-max_abs, max_abs), low = "brown4", mid = "white", high = "midnightblue")+
    scale_x_discrete(expand = c(0, 0)) + # CRITICAL: Remove blank padding on X axis
    theme_classic(base_size = 6)+
    theme(
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      axis.ticks = element_line(linewidth = 0.05),
      axis.ticks.x = element_blank(), # Cleaner look without bottom ticks if sample names are hidden
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt"),
      axis.line = element_line(linewidth = 0.2)
    )+
    labs(fill = "Z-Score")
  return(p)
}

plot_header_heatmap <- function(heatmap_df, main_var, color_palette = NULL){
  main_var_name <- set_names(main_var)
  heatmap_df <- heatmap_df %>% select(id, variable) %>% unique()
  
  p <- ggplot(heatmap_df, aes(x = id, y = 1, fill = variable))+
    geom_tile(color = "white", linewidth = 0.1)+
    theme_classic(base_size = 6)+
    theme(
      axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt")
    )+
    labs(fill = main_var_name)+
    scale_x_discrete(expand = c(0, 0))
  
  if (!is.null(color_palette)) {
    p <- p + scale_fill_manual(values = color_palette)
  } else {
    p <- p + scale_fill_viridis_d(option = "plasma",
                                  begin = 0.2,
                                  end = 0.8)
  }
  return(p)
}

plot_heatmap <- function(deg_table, design_metadata, v, main_var, color_palette = NULL, top = 50){
  main_var_name <- set_names(main_var)
  
  norm_matrix <- v$E
  deg_genes <- deg_table %>% arrange(adj.P.Val) %>% head(n = top) %>% arrange(logFC) %>% pull(miRNA)
  heatmap_mat <- norm_matrix[deg_genes, , drop = FALSE]
  heatmap_df <- heatmap_mat %>% as.data.frame() %>% rownames_to_column("miRNA") %>% pivot_longer(cols = colnames(heatmap_mat), values_to = "expression", names_to = "id") %>% left_join(design_metadata)
  
  heatmap_df <- heatmap_df[,c("miRNA","id", "expression", main_var)]
  colnames(heatmap_df)[4] <- "variable"
  
  heatmap_df <- heatmap_df %>% arrange(variable) %>% mutate(id = factor(id, levels = unique(id)),
                                                            miRNA = factor(miRNA, levels = deg_genes))
  
  body <- plot_body_heatmap(heatmap_df)
  header <- plot_header_heatmap(heatmap_df, main_var, color_palette)
  
  leg_header <- as_ggplot(get_legend(header))
  leg_body   <- as_ggplot(get_legend(body))
  legends    <- plot_grid(leg_header, leg_body, ncol = 1, rel_heights = c(1, 1), align = "v")
  
  header <- header + theme(legend.position = "none")
  body <- body + theme(legend.position = "none")
  
  heat <- plot_grid(header, body, align = "v", ncol = 1, rel_heights = c(0.1, 2))
  
  f_heat <- plot_grid(heat, legends, nrow = 1, rel_widths = c(1, 0.25))
  return(f_heat)
}

plot_volcano <- function(de_table, contrast_name){
  volcano_colors <- c(
    "Up"              = "#D32F2F", # Rojo para los sobreexpresados
    "Down"            = "#1976D2", # Azul para los subexpresados
    "Not Significant" = "#9E9E9E"  # Gris para el ruido de fondo
  )
  p <- ggplot(de_table, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Significance), size = 0.8, alpha = 0.6) +
    
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.15, color = "black") +
    geom_vline(xintercept = -1, linetype = "dashed", linewidth = 0.15, color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.15, color = "black") +
    
    scale_color_manual(values = volcano_colors) +
    
    labs(
      title = paste("Volcano Plot:", contrast_name),
      x = "log2 (Fold Change)",
      y = "-log10 (Adjusted P-Value)",
      color = "Status"
    ) +
    theme_classic(base_size = 6)+
    theme(
      legend.key.size = unit(0.2, "cm")
    )
}


diff_expression <- function(formula, formula_name, design_metadata, dge, color_palette = NULL){
  
  design <- model.matrix(as.formula(formula), data = design_metadata)
  colnames(design) <- make.names(colnames(design))
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  
  formula_parts <- str_split(formula, "\\+")[[1]]
  formula_parts <- trimws(formula_parts)
  formula_parts <- formula_parts[formula_parts != "~ 0"]
  main_var <- formula_parts[1]
  all_var_names <- all.vars(as.formula(formula))
  other_vars <- all_var_names[all_var_names != main_var]
  target_cols <- colnames(design)[str_detect(colnames(design), main_var)]
  if(length(other_vars) > 0) {
    for(ov in other_vars) {
      target_cols <- target_cols[!str_detect(target_cols, ov)]
    }
  }
  
  combs <- as.data.frame(t(combn(target_cols, 2))) %>% 
    mutate(contrast_str = paste0(V2, " - ", V1))
  
  for(i in 1:nrow(combs)) {
    contrast_name <- combs$contrast_str[i]
    
    # Limpiar nombres para los archivos (ej: "trimesterThird - trimesterFirst" -> "Third-First")
    clean_name <- str_replace_all(contrast_name, " ", "")
    clean_name <- gsub(main_var, "", clean_name)
    
    file_prefix <- paste0(formula_name, "_", clean_name)
    
    # Hacer el contraste matemático controlando el sexo fetal
    contrast_matrix <- makeContrasts(contrasts = contrast_name, levels = design)
    fit_contrast <- contrasts.fit(fit, contrast_matrix)
    fit_contrast <- eBayes(fit_contrast)
    
    # Extraer la tabla completa de resultados
    de_table <- topTable(fit_contrast, coef = 1, number = Inf, sort.by = "P") %>%
      rownames_to_column("miRNA") %>%
      as_tibble() %>%
      mutate(
        Significance = case_when(
          adj.P.Val < 0.05 & logFC > 1  ~ "Up",
          adj.P.Val < 0.05 & logFC < -1 ~ "Down",
          TRUE                          ~ "Not Significant"
        )
      )
    
    deg_table <- de_table %>% filter(Significance != "Not Significant")
    
    if (nrow(deg_table) > 1){
      output_file <- paste0("results/differential_expression/DEG_", file_prefix, ".csv")
      write_csv(deg_table, output_file)
      p <- plot_heatmap(deg_table, design_metadata, v, main_var, color_palette[[main_var]])
      ggsave(filename = paste0("results/differential_expression/DEG_", file_prefix, "_heatmap.png"), plot = p,
             height = 7, width = 10, units = "cm", bg = "white")
      p <- plot_volcano(de_table, contrast_name)
      ggsave(filename = paste0("results/differential_expression/DEG_", file_prefix, "_volcano.png"), plot = p,
             height = 7, width = 10, units = "cm", bg = "white")
    }
    output_file <- paste0("results/differential_expression/DE_", file_prefix, ".csv")
    write_csv(de_table, output_file)
  }
}
