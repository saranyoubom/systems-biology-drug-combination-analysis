# C-peptide Secretion Profile - Figure 3a

rm(list = ls())

library(ggplot2)
library(readr)
library(ggsci)
library(readxl)

dat2 <- read_excel("Fig_3a_Scatterplot_C_peptide.xlsx", sheet = 1)
attach(dat2)

dat2$Glucose <- as.numeric(dat2$Glucose)
dat2 <- filter(dat2, Glucose != "44")

plot_loess <- ggplot(dat2, aes(x = Glucose, y = C_peptide, color = Treatment)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  geom_point(size = 1.5, alpha = 0.5) +
  facet_grid(cols = vars(Assay)) +
  scale_color_npg() +
  labs(x = "Glucose concentration (mM)",
       y = "\nRelative C-peptide release (AU)") +
  lims(y = c(0, 320)) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        strip.text = element_text(size = 8, colour = "black", face = "bold"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = margin(5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8))

plot_loess

ggsave("Fig_3a_Scatterplot_C_peptide.png", plot = plot_loess,
       scale = 1, width = 8, height = 6.5, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_3a_Scatterplot_C_peptide.svg", plot = plot_loess,
       scale = 1, width = 8, height = 6.5, units = "cm",
       dpi = 600, limitsize = TRUE)
