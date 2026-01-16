# Calcium Signaling Gene Expression - Figure 3e

rm(list = ls())

library(ggplot2)
library(dplyr)
library(ggsci)
library(readxl)

dat2 <- read_excel("Fig_3e_Barplot_Log_Fold_Change_mRNA_expression.xlsx", sheet = 1)
attach(dat2)

dat3 <- filter(dat2, Marker == "Intracellular markers")
dat3 <- filter(dat3, Treatment != "CTRL")

# Calculate mean and SEM
dat_summary <- dat3 %>%
  group_by(Gene, Treatment) %>%
  summarise(
    mean_LogFC = mean(LogFC),
    sem_LogFC = sd(LogFC) / sqrt(n()),
    .groups = "drop"
  )

# Get NPG palette colors (2nd and 3rd)
npg_colors <- pal_npg()(10)
custom_colors <- c(npg_colors[2], npg_colors[3])

dat_summary$Gene <- factor(dat_summary$Gene, levels = c("Ptbp1", "Itpr1",
                                                         "Ryr1", "Ryr2", "Ryr3"))

plot_bar <- ggplot(dat_summary, aes(x = Gene, y = mean_LogFC, fill = Treatment, color = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_LogFC - sem_LogFC, ymax = mean_LogFC + sem_LogFC),
                width = 0.2, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(x = "Gene", y = "Log2(Fold change)\nof mRNA expression") +
  lims(y = c(-1.5, 0.5)) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = margin(5),
        axis.title.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5, face = "italic"),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8))

plot_bar

ggsave("Fig_3e_Barplot_Log_Fold_Change_mRNA_expression.png", plot = plot_bar,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_3e_Barplot_Log_Fold_Change_mRNA_expression.svg", plot = plot_bar,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)
