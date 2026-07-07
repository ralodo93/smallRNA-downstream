library(yaml)
library(tidyverse)
library(edgeR)
library(ggdendro)
library(pheatmap)

config <- read_yaml("config.yml")
load("results/rds/raw.rdata")


volcano_plot <- function(de_table, contrast){
  mfc <- max(abs(de_table$logFC))
  p <- ggplot(de_table, aes(x = logFC, y = -log10(P.Value)))+
    geom_point(aes(fill = Significance), size = 2, alpha = 0.8, shape = 21, color = "gray20")+
    scale_x_continuous(limits = c(-mfc, mfc))+
    labs(x = "Log Fold Change", y = "-log10 (pvalue)", title = contrast)+
    theme_classic(base_size = 8) +
    scale_fill_manual(values = c("Not Significant" = "gray80", "Up" = "deepskyblue3", "Down" = "salmon"))+
    guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1)))
  return(p)
}

heatmap_plot <- function(matrix_subset, metadata, contrast, deg_table, config){
  matrix_scaled <- t(scale(t(matrix_subset)))

  # clustering cols
  hc <- hclust(dist(t(matrix_scaled)))
  hcdata <- dendro_data(hc)
  hc_labels <- label(hcdata) |> mutate(id = label) |> inner_join(metadata |> rownames_to_column("id"))
  hc_segments <- segment(hcdata)

  n_samples <- length(hc$order)
  x_limits <- c(0.5, n_samples + 0.5)

  # dendogram
  headers <- list()
  p1 <- ggplot()+
    geom_segment(data = hc_segments, aes(x = x, xend = xend, y = y, yend = yend))+
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    theme_void()+
    theme(plot.margin = margin(0.3, 0, 0, 0, "cm"))
  headers[[1]] <- p1

  # headers
  i <- 2
  for (var in config$qc_variables){
    p2 <- ggplot()+
      geom_tile(data = hc_labels, aes(x = x, y = 1, fill = !!sym(var)), color = "black") +
      scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
      theme_void()+
      theme(plot.margin = margin(0, 0.4, 0.2, 0.4, "cm"),
            legend.position = "right",
            legend.margin = margin(0,0,0,0),
            legend.box.margin = margin(0,0,0,0))
    if (var %in% names(config$color_variables)){
      color <- unlist(config$color_variables[[var]])
      p2 <- p2 + scale_fill_manual(values = color)
    }
    headers[[i]] <- p2
    i <- i + 1
  }

  # tranform matrix to pivot_longer
  long_matrix <- matrix_scaled |> as.data.frame() |> rownames_to_column("miRNA") |>
    pivot_longer(cols = colnames(matrix_scaled), names_to = "id", values_to = "zscore") |>
    inner_join(hc_labels) |> mutate(miRNA = factor(miRNA, levels = deg_table$miRNA))

  lim <- max(abs(long_matrix$zscore), na.rm = TRUE)
  p3 <- ggplot()+
    geom_tile(data = long_matrix, aes(x = x, y = miRNA, fill = zscore))+
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    theme_void(base_size = 8)+
    theme(axis.text.y = element_text(hjust = 1), plot.margin = margin(0, 0, 0.3, 0, "cm"),
          legend.position = "right",
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(0,0,0,0))+
    scale_fill_gradient2(
      low = "#2166AC",      # azul
      mid = "white",
      high = "#B2182B",     # rojo
      midpoint = 0,
      name = "zscore",
      limits = c(-lim, lim)
    )

  headers[[i]] <- p3

  body_p <-
    patchwork::wrap_plots(
      headers,
      ncol = 1,
      heights = c(
        1,
        rep(0.5, length(config$qc_variables)),
        8
      )
    ) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = contrast
    ) &
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(2, "mm")
    )
  return(body_p)
}

get_top_table <- function(fit, coef){
  de_table <- topTable(fit, coef = coef, number = Inf, sort.by = "P") %>% 
    rownames_to_column("miRNA") %>%
    as_tibble() %>%
    mutate(
      Significance = case_when(
        adj.P.Val < 0.05 & logFC > 1  ~ "Up",
        adj.P.Val < 0.05 & logFC < -1 ~ "Down",
        TRUE                           ~ "Not Significant"
      )
    )
  return(de_table)
}

prepare_single_scenario <- function(variable, d){
  formula <- paste0("~ 0 + contrast")
  group <- d$samples[[variable]]
  d_copy <- d
  d_copy$samples$contrast <- group
  scenario <- list(d = d_copy, formula = formula)
  return(scenario)
}

prepare_adjusted_scenario <- function(variables, d){
  main_var <- variables[1]
  others <- variables[-1]
  formula <- paste0("~ 0 + contrast + ", paste(others, collapse = " + "))
  group <- d$samples[[main_var]]
  d_copy <- d
  d_copy$samples$contrast <- group
  scenario <- list(d = d_copy, formula = formula)
  return(scenario)
}

prepare_interaction_single_scenario <- function(variables, d){
  formula <- paste0("~ 0 + contrast")
  group <- interaction(d$samples[,variables], sep = "_")
  d_copy <- d
  d_copy$samples$contrast <- group
  scenario <- list(d = d_copy, formula = formula)
  return(scenario)
}


prepare_interaction_adjusted_scenario <- function(formula, d){
  formula_split <- unlist(strsplit(formula, "\\*"))
  var1 <- all.vars(as.formula(formula_split[1]))
  others <- all.vars(as.formula(paste0("~ ", formula_split[2])))
  var2 <- others[1]
  others <- others[-1]
  formula <- paste0("~ 0 + contrast + ", paste(others, collapse = " + "))
  group <- interaction(d$samples[,c(var1, var2)], sep = "_")
  d_copy <- d
  d_copy$samples$contrast <- group
  scenario <- list(d = d_copy, formula = formula)
  return(scenario)
}

de_analysis <- function(form_scenario, d_scenario, contrasts){
  data <- d_scenario$samples
  mm <- model.matrix(as.formula(form_scenario), data = data)
  colnames(mm) <- sub("^contrast", "", colnames(mm))
  y <- voom(d_scenario, mm)
  fit <- lmFit(y, mm)
  contrast_matrix <- makeContrasts(contrasts = contrasts, levels = mm)
  fit <- contrasts.fit(fit, contrast_matrix)
  fit <- eBayes(fit)
  return(fit)
}

extract_groups <- function(contrast){
  
  groups <- strsplit(
    gsub("[()]", "", contrast),
    "[-+]"
  )[[1]]
  
  groups <- trimws(groups)
  groups <- groups[groups != ""]
  
  unique(groups)
}

de_results <- function(fit, contrasts, d_scenario, config){
  de_tables <- list()
  volcano_plots <- list()
  heatmap_plots <- list()
  for (contrast in contrasts){
    de_table <- get_top_table(fit, contrast)
    de_tables[[contrast]] <- de_table
    vp <- volcano_plot(de_table, contrast)
    volcano_plots[[contrast]] <- vp
    
    deg_table <- de_table |> dplyr::filter(Significance != "Not Significant") |> arrange(adj.P.Val) |> head(n = 30) |> arrange(logFC)
    
    if (nrow(deg_table) > 5){
      selected_groups <- extract_groups(contrast)
      samples <- rownames(d_scenario$samples[d_scenario$samples$contrast %in% selected_groups,])
      log2_cpm <- cpm(d_scenario, prior.count=2, log=TRUE)[deg_table$miRNA,samples]
      hp <- heatmap_plot(log2_cpm, d_scenario$samples, contrast, deg_table, config)
      heatmap_plots[[contrast]] <- hp
    }
    
  }
  return(list(de_tables = de_tables, volcano_plots = volcano_plots, heatmap_plots = heatmap_plots))
}

make_contrasts <- function(formula, contrasts, d, config){
  variables <- all.vars(as.formula(formula))
  
  # check if all variables are columns of metadata; if not launch an error message
  if (!all(variables %in% colnames(d$samples))){
    stop("variables inside formula must be columns of metadata")
  }
  
  # check if there are duplicated variables
  if (anyDuplicated(variables)) {
    stop("Variables in the formula must be unique.")
  }
  
  if (length(variables) == 1){ # ~ genotype
    scenario <- prepare_single_scenario(variables, d)
  } else{
    if (!grepl("\\*", formula)){ # ~ genotype + batch
      scenario <- prepare_adjusted_scenario(variables, d)
    } else{
      
      if (!grepl("\\+", formula)){ # ~ genotype * condition
        scenario <- prepare_interaction_single_scenario(variables, d)
      } else{ # ~ genotype * condition + batch
        scenario <- prepare_interaction_adjusted_scenario(formula, d)
      }
    }
  }
  
  form_scenario <- scenario$formula
  d_scenario <- scenario$d
  
  fit <- de_analysis(form_scenario, d_scenario, contrasts)
  de_tables <- de_results(fit, contrasts, d_scenario, config)
  
  return(de_tables)
}

results <- list()
plot_volcano <- list()
plot_heatmap <- list()

for (formula in names(config$formulas)){
  contrasts <- config$formulas[[formula]]
  results_formula <- make_contrasts(formula, contrasts, d, config)
  results[[formula]] <- results_formula$de_tables
  plot_volcano[[formula]] <- results_formula$volcano_plots
  plot_heatmap[[formula]] <- results_formula$heatmap_plots
}
