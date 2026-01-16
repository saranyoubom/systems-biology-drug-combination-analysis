# C-peptide at Physiological Glucose - Figure 3b

rm(list = ls())

library(ggplot2)
library(data.table)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)

dat <- read_excel("Fig_3b_Boxplot_C_peptide.xlsx", sheet = 1)
attach(dat)

my_comparison1 <- list(c("CTRL", "DAPT"), c("CTRL", "DKK-1"), c("DAPT", "DKK-1"))

dat1 <- dat
dat1 <- filter(dat1, Glucose == "5.5 mM glucose")

glucose.labs <- c("5.5 mM glucose (Physiological level)")
names(glucose.labs) <- c("5.5 mM glucose")

shapiro.test(dat1$C_peptide) # W = 0.6920, p = 0.0007 (non-normal)

dat1$Glucose <- factor(dat1$Glucose, levels = c("Basal medium", "2.8 mM glucose",
                                                 "5.5 mM glucose", "22 mM glucose",
                                                 "44 mM glucose"))

plot1 <- ggboxplot(dat1, x = "Treatment", y = "C_peptide",
                   palette = "npg",
                   color = "Treatment", fill = "Treatment", alpha = 0.5,
                   add = "jitter",
                   ylim = c(0, 450)) +
  labs(y = "\nC-peptide release (ng/µg DNA)") +
  facet_grid(cols = vars(Glucose), labeller = labeller(Glucose = glucose.labs)) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = margin(5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        strip.text = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8)) +
  stat_compare_means(label.y = 448, label.x.npc = "center", size = 2.5,
                     method = "kruskal.test", hjust = 0.5) +
  stat_compare_means(comparisons = my_comparison1, method = "wilcox.test",
                     tip.length = 0.05, step.increase = 0.25,
                     aes(label = after_stat(p.signif)), paired = FALSE, size = 2.5) +
  stat_summary(fun.data = function(x) {
    mean_value <- mean(x)
    sem_value <- sd(x) / sqrt(length(x))
    data.frame(y = 0, label = paste0(round(mean_value, 2), " ± ", round(sem_value, 2)))
  }, geom = "text", vjust = 1, size = 2.5)

plot1

ggsave("Fig_3b_Boxplot_C_peptide.png", plot = plot1,
       scale = 1, width = 8, height = 6.5, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_3b_Boxplot_C_peptide.svg", plot = plot1,
       scale = 1, width = 8, height = 6.5, units = "cm",
       dpi = 600, limitsize = TRUE)
