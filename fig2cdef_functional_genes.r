# Functional Gene Expression - Figure 2c-f

rm(list = ls())

library(ggplot2)
library(data.table)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)

dat1 <- read_excel("Fig_2cdef_Fold_change_mRNA_expression.xlsx", sheet = 1)
attach(dat1)

my_comparison1 <- list(c("CTRL", "DAPT"), c("CTRL", "DKK-1"), c("DAPT", "DKK-1"))

dat.Glut2 <- filter(dat1, Gene == "Glut2")
dat.Kcnj11 <- filter(dat1, Gene == "Kcnj11")
dat.Cacna1c <- filter(dat1, Gene == "Cacna1c")
dat.Cacna1d <- filter(dat1, Gene == "Cacna1d")

# Glut2 (Glucose uptake and sensing)
dat1.Glut2 <- dat.Glut2
shapiro.test(dat1.Glut2$Fold_change) # W = 0.9017, p = 0.1667 (normal)

Glut2.labs <- c("Glut2 (Glucose uptake and sensing)")
names(Glut2.labs) <- c("Glut2")

plot.Glut2 <- ggboxplot(dat1.Glut2, x = "Treatment", y = "Fold_change",
                        palette = "npg",
                        color = "Treatment", fill = "Treatment", alpha = 0.5,
                        add = "jitter",
                        ylim = c(-0.1, 2)) +
  labs(y = "Fold change in\nmRNA expression (AU)") +
  facet_wrap(vars(Gene), labeller = labeller(Gene = Glut2.labs)) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, colour = "black", face = "bold.italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8)) +
  stat_compare_means(label.y = 2, label.x.npc = "center", size = 2.5,
                     method = "anova", hjust = 0.5) +
  stat_compare_means(comparisons = my_comparison1, method = "t.test",
                     tip.length = 0.1, step.increase = 0.5,
                     aes(label = after_stat(p.signif)), paired = FALSE, size = 2.5) +
  stat_summary(fun.data = function(x) {
    mean_value <- mean(x)
    sem_value <- sd(x) / sqrt(length(x))
    data.frame(y = 0, label = paste0(round(mean_value, 2), " ± ", round(sem_value, 2)))
  }, geom = "text", vjust = 1, size = 2.5)

plot.Glut2

ggsave("Fig_2c_Boxplot_Fold_change_Glut2_expression.png", plot = plot.Glut2,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_2c_Boxplot_Fold_change_Glut2_expression.svg", plot = plot.Glut2,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

# Kcnj11 (ATP-sensitive potassium channel)
dat1.Kcnj11 <- dat.Kcnj11
dat1.Kcnj11$Fold_change <- as.numeric(dat1.Kcnj11$Fold_change)
shapiro.test(dat1.Kcnj11$Fold_change) # W = 0.8584, p = 0.0467 (normal)

Kcnj11.labs <- c("Kcnj11 (ATP-sensitive potassium channel)")
names(Kcnj11.labs) <- c("Kcnj11")

plot.Kcnj11 <- ggboxplot(dat1.Kcnj11, x = "Treatment", y = "Fold_change",
                         palette = "npg",
                         color = "Treatment", fill = "Treatment", alpha = 0.5,
                         add = "jitter",
                         ylim = c(-0.1, 2)) +
  labs(y = "Fold change in\nmRNA expression (AU)") +
  facet_wrap(vars(Gene), labeller = labeller(Gene = Kcnj11.labs)) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, colour = "black", face = "bold.italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8)) +
  stat_compare_means(label.y = 2, label.x.npc = "center", size = 2.5,
                     method = "anova", hjust = 0.5) +
  stat_compare_means(comparisons = my_comparison1, method = "t.test",
                     tip.length = 0.05, step.increase = 0.4,
                     aes(label = after_stat(p.signif)), paired = FALSE, size = 2.5) +
  stat_summary(fun.data = function(x) {
    mean_value <- mean(x)
    sem_value <- sd(x) / sqrt(length(x))
    data.frame(y = 0, label = paste0(round(mean_value, 2), " ± ", round(sem_value, 2)))
  }, geom = "text", vjust = 1, size = 2.5)

plot.Kcnj11

ggsave("Fig_2d_Boxplot_Fold_change_Kcnj11_expression.png", plot = plot.Kcnj11,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_2d_Boxplot_Fold_change_Kcnj11_expression.svg", plot = plot.Kcnj11,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

# Cacna1c (Calcium voltage-gated channel)
dat1.Cacna1c <- dat.Cacna1c
dat1.Cacna1c$Fold_change <- as.numeric(dat1.Cacna1c$Fold_change)
shapiro.test(dat1.Cacna1c$Fold_change) # W = 0.8395, p = 0.0273 (non-normal)

Cacna1c.labs <- c("Cacna1c (Calcium voltage-gated channel)")
names(Cacna1c.labs) <- c("Cacna1c")

plot.Cacna1c <- ggboxplot(dat1.Cacna1c, x = "Treatment", y = "Fold_change",
                          palette = "npg",
                          color = "Treatment", fill = "Treatment", alpha = 0.5,
                          add = "jitter",
                          ylim = c(-0.1, 2)) +
  labs(y = "Fold change in\nmRNA expression (AU)") +
  facet_wrap(vars(Gene), labeller = labeller(Gene = Cacna1c.labs)) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, colour = "black", face = "bold.italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8)) +
  stat_compare_means(label.y = 2, label.x.npc = "center", size = 2.5,
                     method = "kruskal.test", hjust = 0.5) +
  stat_compare_means(comparisons = my_comparison1, method = "wilcox.test",
                     tip.length = 0.1, step.increase = 0.8,
                     aes(label = after_stat(p.signif)), paired = FALSE, size = 2.5) +
  stat_summary(fun.data = function(x) {
    mean_value <- mean(x)
    sem_value <- sd(x) / sqrt(length(x))
    data.frame(y = 0, label = paste0(round(mean_value, 2), " ± ", round(sem_value, 2)))
  }, geom = "text", vjust = 1, size = 2.5)

plot.Cacna1c

ggsave("Fig_2e_Boxplot_Fold_change_Cacna1c_expression.png", plot = plot.Cacna1c,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_2e_Boxplot_Fold_change_Cacna1c_expression.svg", plot = plot.Cacna1c,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

# Cacna1d (Calcium voltage-gated channel)
dat1.Cacna1d <- dat.Cacna1d
dat1.Cacna1d$Fold_change <- as.numeric(dat1.Cacna1d$Fold_change)
shapiro.test(dat1.Cacna1d$Fold_change) # W = 0.8832, p = 0.0964 (normal)

Cacna1d.labs <- c("Cacna1d (Calcium voltage-gated channel)")
names(Cacna1d.labs) <- c("Cacna1d")

plot.Cacna1d <- ggboxplot(dat1.Cacna1d, x = "Treatment", y = "Fold_change",
                          palette = "npg",
                          color = "Treatment", fill = "Treatment", alpha = 0.5,
                          add = "jitter",
                          ylim = c(-0.1, 2)) +
  labs(y = "Fold change in\nmRNA expression (AU)") +
  facet_wrap(vars(Gene), labeller = labeller(Gene = Cacna1d.labs)) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, colour = "black", face = "bold.italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8)) +
  stat_compare_means(label.y = 2, label.x.npc = "center", size = 2.5,
                     method = "anova", hjust = 0.5) +
  stat_compare_means(comparisons = my_comparison1, method = "t.test",
                     tip.length = 0.1, step.increase = 0.8,
                     aes(label = after_stat(p.signif)), paired = FALSE, size = 2.5) +
  stat_summary(fun.data = function(x) {
    mean_value <- mean(x)
    sem_value <- sd(x) / sqrt(length(x))
    data.frame(y = 0, label = paste0(round(mean_value, 2), " ± ", round(sem_value, 2)))
  }, geom = "text", vjust = 1, size = 2.5)

plot.Cacna1d

ggsave("Fig_2f_Boxplot_Fold_change_Cacna1d_expression.png", plot = plot.Cacna1d,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_2f_Boxplot_Fold_change_Cacna1d_expression.svg", plot = plot.Cacna1d,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)
