library(yaml)
library(tidyverse)
library(edgeR)
library(ggdendro)

config <- read_yaml("config.yml")
load("results/rds/raw.rdata")


## PCA plots (samples from same group must be together in 2D map)
pca_res <- prcomp(t(log2_cpm))
var_explained <- round((pca_res$sdev^2) / sum(pca_res$sdev^2) * 100, 2)

pca_coords <- pca_res$x[,c(1,2)] |> as.data.frame() |> rownames_to_column("id") |>
  inner_join(metadata |> rownames_to_column("id"))

pca_plots <- list()

## Clustering plots (samples from same group must be near and in the same brach)
hc <- hclust(dist(t(log2_cpm)))
hcdata <- dendro_data(hc)
hc_labels <- label(hcdata) |> mutate(id = label) |> inner_join(metadata |> rownames_to_column("id"))
hc_segments <- segment(hcdata)

clust_plots <- list()

for (var in config$qc_variables){
  var_sym <- sym(var)
  
  # pca
  xlabel <- paste0("PC 1 ( ", var_explained[1], " % )")
  ylabel <- paste0("PC 2 ( ", var_explained[2], " % )")
  ppca <- ggplot(pca_coords, aes(x = PC1, y = PC2))+
    geom_point(size = 3, alpha = 0.75, aes(color = !!var_sym))+
    theme_classic(base_size = 8)+
    theme(plot.title = element_text(face = "bold")) + 
    labs(x = xlabel, y = ylabel, color = tools::toTitleCase(gsub("_", " ", var))) +
    ggtitle(tools::toTitleCase(gsub("_", " ", var)))+
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
  
  # clustering
  maxy <- max(hc_segments$y)
  limits <- c(-maxy/3, maxy*1.1)
  pclust <- ggplot()+
    geom_segment(data = hc_segments, aes(x = x, xend = xend, y = y, yend = yend))+
    geom_text(data = hc_labels, aes(x = x, y = y, label = label), angle = 90, hjust = 1.2, vjust = 0.5)+
    geom_point(data = hc_labels, aes(x = x, y = y, color = !!var_sym), size = 4)+
    scale_y_continuous(limits = limits)+
    theme_minimal(base_size = 8)+
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title = element_text(face = "bold"))+
    labs(y = "Heigth", color = tools::toTitleCase(gsub("_", " ", var)))+
    ggtitle(tools::toTitleCase(gsub("_", " ", var)))+
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
  
  if (var %in% names(config$color_variables)){
    colors <- unlist(config$color_variables[[var]])
    ppca <- ppca + scale_color_manual(values = colors)
    pclust <- pclust + scale_color_manual(values = colors)
  }
  
  pca_plots[[var]] <- ppca
  clust_plots[[var]] <- pclust
  
}

cowplot::plot_grid(plotlist = pca_plots)
cowplot::plot_grid(plotlist = clust_plots)


