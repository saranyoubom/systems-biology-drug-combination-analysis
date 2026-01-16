# MLSS v4.0 Validation - Sensitivity and Robustness Analysis

rm(list = ls())
set.seed(42)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(reshape2)
library(gridExtra)

npg_colors <- pal_npg("nrc")(10)

# Load MLSS output
mlss_data <- read.csv("fig6c_mlss_all_combinations.csv", stringsAsFactors = FALSE)
top_15 <- mlss_data %>% head(15)

# Define weight scenarios
weight_scenarios <- data.frame(
  Scenario = c("MLSS (Standard)", "C-Heavy (60:25:15)", "C-Heavy (55:20:25)",
               "B-Heavy (35:50:15)", "V-Heavy (35:30:35)", "P-Heavy (30:20:10:40)",
               "Balanced (40:40:10:10)", "Network-Focus (30:20:40:10)"),
  C_Weight = c(0.45, 0.60, 0.55, 0.35, 0.35, 0.30, 0.40, 0.30),
  B_Weight = c(0.30, 0.25, 0.20, 0.50, 0.30, 0.20, 0.40, 0.20),
  V_Weight = c(0.10, 0.10, 0.15, 0.10, 0.30, 0.10, 0.10, 0.40),
  P_Weight = c(0.15, 0.05, 0.10, 0.05, 0.05, 0.40, 0.10, 0.10),
  NPG_Color = c(npg_colors[1:8]),
  stringsAsFactors = FALSE
)

# Recalculate scores for each scenario
top_15_data <- mlss_data %>% head(15)
scenario_results <- data.frame(Combination = top_15_data$Combination,
                                Original_MLSS = top_15_data$MLSS_v4.0)

for (i in 1:nrow(weight_scenarios)) {
  scenario_name <- weight_scenarios$Scenario[i]
  c_w <- weight_scenarios$C_Weight[i]
  b_w <- weight_scenarios$B_Weight[i]
  v_w <- weight_scenarios$V_Weight[i]
  p_w <- weight_scenarios$P_Weight[i]

  new_base_score <- (c_w * top_15_data$Complementarity_C +
                     b_w * top_15_data$Balance_B +
                     v_w * top_15_data$Network_Coverage_V +
                     p_w * top_15_data$Potency_P) * 9

  new_mlss <- new_base_score + top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty
  scenario_results[[scenario_name]] <- new_mlss
}

# Ranking correlation analysis
rank_data <- apply(scenario_results[, -1], 2, rank)
rank_correlation <- cor(rank_data, method = "spearman")

# Ablation study
ablation_results <- data.frame(Combination = top_15_data$Combination,
                               Original_MLSS = top_15_data$MLSS_v4.0,
                               stringsAsFactors = FALSE)

ablation_results$Without_C <- ((0 * top_15_data$Complementarity_C +
                                0.30/0.55 * top_15_data$Balance_B +
                                0.10/0.55 * top_15_data$Network_Coverage_V +
                                0.15/0.55 * top_15_data$Potency_P) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Without_B <- ((0.45/0.70 * top_15_data$Complementarity_C +
                                0 * top_15_data$Balance_B +
                                0.10/0.70 * top_15_data$Network_Coverage_V +
                                0.15/0.70 * top_15_data$Potency_P) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Without_V <- ((0.45/0.90 * top_15_data$Complementarity_C +
                                0.30/0.90 * top_15_data$Balance_B +
                                0 * top_15_data$Network_Coverage_V +
                                0.15/0.90 * top_15_data$Potency_P) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Without_P <- ((0.45/0.85 * top_15_data$Complementarity_C +
                                0.30/0.85 * top_15_data$Balance_B +
                                0.10/0.85 * top_15_data$Network_Coverage_V +
                                0 * top_15_data$Potency_P) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Base_Only <- ((0.45 * top_15_data$Complementarity_C +
                                0.30 * top_15_data$Balance_B +
                                0.10 * top_15_data$Network_Coverage_V +
                                0.15 * top_15_data$Potency_P) * 9)

# Calculate ablation impact
ablation_impact <- data.frame(
  Component = c("Complementarity (C)", "Balance (B)", "Coverage (V)", "Potency (P)", "Synergy/Antagonism"),
  Avg_Rank_Shift = NA,
  Avg_Score_Drop = NA,
  Top5_Stability = NA,
  NPG_Color = npg_colors[1:5],
  stringsAsFactors = FALSE
)

ablation_cols <- c("Without_C", "Without_B", "Without_V", "Without_P", "Base_Only")
for (i in 1:length(ablation_cols)) {
  ablation_col <- ablation_cols[i]
  rank_diff <- abs(rank(ablation_results$Original_MLSS) - rank(ablation_results[[ablation_col]]))
  ablation_impact$Avg_Rank_Shift[i] <- mean(rank_diff, na.rm = TRUE)
  score_drop <- ablation_results$Original_MLSS - ablation_results[[ablation_col]]
  ablation_impact$Avg_Score_Drop[i] <- mean(score_drop, na.rm = TRUE)
  top5_original <- which(rank(ablation_results$Original_MLSS, ties.method = "first") <= 5)
  top5_ablated <- which(rank(ablation_results[[ablation_col]], ties.method = "first") <= 5)
  overlap <- length(intersect(top5_original, top5_ablated))
  ablation_impact$Top5_Stability[i] <- overlap / 5
}

# Visualization 1: Sensitivity heatmap
heatmap_data <- scenario_results %>%
  select(-Original_MLSS) %>%
  arrange(match(Combination, top_15$Combination))
rownames(heatmap_data) <- heatmap_data$Combination
heatmap_data <- heatmap_data[, -1]

heatmap_normalized <- apply(heatmap_data, 2, function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
})

heatmap_long <- reshape2::melt(as.matrix(heatmap_normalized))
colnames(heatmap_long) <- c("Combination", "Scenario", "Normalized_Score")

p_heatmap <- ggplot(heatmap_long, aes(x = Scenario, y = Combination, fill = Normalized_Score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = npg_colors[6], high = npg_colors[1],
                      name = "Normalized\nMLSS Score") +
  labs(title = "Weight Scenario Sensitivity (Top 15)",
       x = "Weight Scenario", y = "Drug Combination") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black", face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave("figs2a_mlss_sensitivity_top15.png", p_heatmap, width = 7.7, height = 6,
       units = "cm", scale = 2.3, dpi = 600, bg = "white")
ggsave("figs2a_mlss_sensitivity_top15.svg", p_heatmap, width = 7.7, height = 6,
       units = "cm", scale = 2.3, dpi = 600, bg = "white")

# Visualization 2: Component importance
importance_plot_data <- ablation_impact %>%
  pivot_longer(cols = c(Avg_Rank_Shift, Avg_Score_Drop),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Avg_Rank_Shift", "Avg_Score_Drop"),
                         labels = c("Avg Rank Shift", "Avg Score Drop")),
         Component = factor(Component, levels = ablation_impact$Component))

p_importance <- ggplot(importance_plot_data, aes(x = reorder(Component, -abs(Value)),
                                                  y = abs(Value), fill = NPG_Color)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3, width = 0.7) +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  scale_fill_identity() +
  labs(title = "Component Importance in MLSS (Top 15)",
       x = "Component",
       y = "Impact Metric (Left: Rank positions; Right: Score units)") +
  coord_flip() +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave("figs2b_component_importance_top15.png", p_importance, width = 7.7, height = 6,
       units = "cm", scale = 2.3, dpi = 600, bg = "white")
ggsave("figs2b_component_importance_top15.svg", p_importance, width = 7.7, height = 6,
       units = "cm", scale = 2.3, dpi = 600, bg = "white")

# Export statistics
write.csv(rank_correlation, "fig6d_ranking_correlations.csv", row.names = TRUE)
write.csv(ablation_impact, "fig6d_ablation_impact.csv", row.names = FALSE)
