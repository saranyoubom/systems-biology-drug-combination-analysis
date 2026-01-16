# Principal Component Analysis - Figure 3f

rm(list = ls())

library(ggplot2)
library(data.table)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)
library(factoextra)
library(ggfortify)

data <- read_excel("Fig_3f_PCA_Relative_mRNA_expression.xlsx", sheet = 1)
attach(data)

# Remove non-numeric columns for PCA
pca_data <- data %>% select(-Group)

x_order <- c("CTRL", "DAPT", "DKK-1")
data$Group <- factor(data$Group, levels = x_order)

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

summary(pca_result)

biplot <- fviz_pca_biplot(pca_result,
                          col.ind = data$Group, palette = "npg",
                          addEllipses = TRUE, label = "var",
                          col.var = "grey20", repel = TRUE,
                          arrowsize = 0.4,
                          legend.title = "Treatment") +
  theme_minimal(base_size = 8) +
  lims(y = c(-5, 4), x = c(-6, 8)) +
  theme(legend.position = "top", title = element_blank(),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = margin(5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8))

biplot

ggsave("Fig_3f_PCA_Relative_mRNA_expression.png", plot = biplot,
       scale = 1, width = 7.6, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_3f_PCA_Relative_mRNA_expression.svg", plot = biplot,
       scale = 1, width = 7.6, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)
