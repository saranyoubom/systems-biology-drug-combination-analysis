# INS1 Cell Proliferation Analysis - Figure 1c

rm(list = ls())

library(ggplot2)
library(data.table)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)

# Load data
dat1 <- read_excel("Fig_1c_Boxplot_INS1_proliferation.xlsx", sheet = 1)
attach(dat1)

my_comparison1 <- list(c("CTRL", "DAPT"), c("CTRL", "DKK-1"), c("DAPT", "DKK-1"))

shapiro.test(dat1$DNA) # W = 0.8916, p = 0.1236 (normal)

plot1 <- ggboxplot(dat1, x = "Treatment", y = "DNA",
                   palette = "npg",
                   color = "Treatment",
                   fill = "Treatment",
                   alpha = 0.5,
                   add = "jitter",
                   ylim = c(0, 0.65)) +
  labs(y = "\nTotal DNA concent (µg)") +
  facet_wrap(vars(Assay)) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, colour = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8)) +
  stat_compare_means(label.y = 0.62, label.x.npc = "center", size = 2.5,
                     method = "anova", hjust = 0.5) +
  stat_compare_means(comparisons = my_comparison1, method = "t.test",
                     tip.length = 0.05, step.increase = 0.25,
                     aes(label = after_stat(p.signif)), paired = FALSE, size = 2.5) +
  stat_summary(fun.data = function(x) {
    mean_value <- mean(x)
    sem_value <- sd(x) / sqrt(length(x))
    data.frame(y = 0, label = paste0(round(mean_value, 2), " ± ", round(sem_value, 2)))
  }, geom = "text", vjust = 0, size = 2.5)

plot1

ggsave("Fig_1c_Boxplot_INS1_proliferation.png", plot = plot1,
       scale = 1, width = 8.08, height = 4.57, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_1c_Boxplot_INS1_proliferation.svg", plot = plot1,
       scale = 1, width = 8.08, height = 4.57, units = "cm",
       dpi = 600, limitsize = TRUE)
