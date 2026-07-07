library(yaml)
library(tidyverse)
library(edgeR)

config <- read_yaml("config.yml")
load("results/rds/raw.rdata")

log2_cpm_df <- log2_cpm |> as.data.frame() |> rownames_to_column("gene") |>
  pivot_longer(cols = colnames(log2_cpm), values_to = "log_cpm", names_to = "id") |>
  inner_join(metadata |> rownames_to_column("id"))

## density plot (it must be a normal distribution)
density_plots <- list()

## boxplot of expression values
boxplot_plots <- list()

for (var in config$qc_variables){
  var_sym <- sym(var)
  
  df_plot <- log2_cpm_df [order(log2_cpm_df[[var]]), ]
  df_plot$id <- factor(df_plot$id, levels = unique(df_plot$id))
  
  pdensity <- ggplot(df_plot, aes(x = log_cpm)) +
    geom_density(aes(color = !!var_sym, group = id), linewidth = 0.2, show.legend = FALSE) +
    stat_density(aes(x = log_cpm, color = !!var_sym),
                 geom = "line", position = "identity", linewidth = 0.7) +
    theme_classic(base_size = 8) +
    theme(legend.position = "none", plot.title = element_text(face = "bold")) + 
    labs(x = "Log CPM", y = "Density") +
    facet_grid(vars(!!var_sym))+
    ggtitle(tools::toTitleCase(gsub("_", " ", var)))
  
  pboxplot <- ggplot(df_plot, aes(x = id, y = log_cpm)) +
    geom_boxplot(aes(color = !!var_sym))+
    theme_classic(base_size = 8)+
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.title.x = element_blank(),
      plot.margin = margin(0.2, 0.2, 0.4, 0.2, "cm")
    )+
    labs(y = "Log CPM", color = tools::toTitleCase(gsub("_", " ", var)))+
    ggtitle(tools::toTitleCase(gsub("_", " ", var)))
  
  if (var %in% names(config$color_variables)){
    colors <- unlist(config$color_variables[[var]])
    pdensity <- pdensity + scale_color_manual(values = colors)
    pboxplot <- pboxplot + scale_color_manual(values = colors)
  }
  density_plots[[var]] <- pdensity
  boxplot_plots[[var]] <- pboxplot
}

cowplot::plot_grid(plotlist = density_plots)
cowplot::plot_grid(plotlist = boxplot_plots)
