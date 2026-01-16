# mRNA Expression Heatmap - Figure 1d

rm(list = ls())

library(ggplot2)
library(pheatmap)
library(dplyr)
library(readxl)
library(scales)
library(ggsci)

dat <- read_excel("Fig_1d_Fold_change_mRNA_expression.xlsx", sheet = 1)
attach(dat)

data <- filter(dat, Marker != "Pancreatic functional markers")

heatmap_plot <- ggplot(data = data, aes(y = Treatment, x = Gene, fill = Fold_change)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#3c5488", "grey98", "#e64b35"),
                       values = scales::rescale(c(0, 1, 2.2)),
                       limits = c(0, 2.2),
                       name = "Fold change in\nmRNA expression (AU)") +
  labs(y = "Treatment", x = "Gene") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = margin(5),
        legend.title = element_text(size = 8, colour = "black", angle = 0, vjust = 0.95, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5, face = "italic"),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5)) +
  scale_x_discrete(limits = c("Hes1", "Hey1", "Wnt2", "Wnt2b", "Wnt5a", "Wnt5b", "Wnt9a",
                               "Lef1", "Tcf7", "Tcf7l2"), position = "bottom") +
  scale_y_discrete(limits = c("CTRL", "DAPT", "DKK-1"), position = "left") +
  geom_text(aes(label = sprintf("%.2f", Fold_change), vjust = 0.5), color = "black", size = 3)

print(heatmap_plot)

ggsave("Fig_1d_Heatmap_Fold_change_mRNA_expression.png", plot = heatmap_plot,
       scale = 1, width = 16.7, height = 7.4, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_1d_Heatmap_Fold_change_mRNA_expression.svg", plot = heatmap_plot,
       scale = 1, width = 16.7, height = 7.4, units = "cm",
       dpi = 600, limitsize = TRUE)
