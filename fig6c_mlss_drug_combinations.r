# MLSS v4.0 - Multi-Layered Synergy Score for Drug Combinations

rm(list = ls())
set.seed(42)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(tibble)

npg_colors <- pal_npg("nrc")(10)

# Load data
hub_data <- read.csv("fig6b_gene_centrality_summary.csv", stringsAsFactors = FALSE)
drug_data <- read.csv("drug.target.interaction.csv", stringsAsFactors = FALSE)
corr_mat_raw <- read.csv("Fig_4e_Hes1_significant_correlations_only.csv",
                         stringsAsFactors = FALSE, check.names = FALSE)

gene_names <- corr_mat_raw[[1]]
corr_mat <- corr_mat_raw[, -1]
rownames(corr_mat) <- toupper(gene_names)
colnames(corr_mat) <- toupper(gene_names)

# Extract potency data (pIC50)
potency_raw <- drug_data %>%
  filter(ACT_TYPE %in% c("IC50", "Ki", "EC50"),
         ORGANISM == "Homo sapiens",
         RELATION == "=",
         !is.na(ACT_VALUE)) %>%
  select(DRUG_NAME, GENE, ACT_VALUE, ACT_TYPE) %>%
  mutate(Gene_Standardized = toupper(GENE))

potency_lookup <- potency_raw %>%
  group_by(DRUG_NAME, Gene_Standardized) %>%
  summarise(pIC50 = max(ACT_VALUE, na.rm = TRUE), .groups = "drop") %>%
  distinct(DRUG_NAME, Gene_Standardized, pIC50)

# Filter for hub genes (degree >= 6)
hub_genes_df <- hub_data %>% filter(Degree >= 6)
hub_genes <- unique(toupper(hub_genes_df$Gene))

# Convert correlation matrix to pairwise format
correlations <- as.data.frame(as.matrix(corr_mat)) %>%
  rownames_to_column(var = "Gene1") %>%
  pivot_longer(cols = -Gene1, names_to = "Gene2", values_to = "r") %>%
  filter(Gene1 != Gene2, !is.na(r)) %>%
  mutate(pair_key = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")) %>%
  distinct(pair_key, .keep_all = TRUE)

correlations_slim <- correlations %>% select(pair_key, r)

# Filter drug-target interactions for hub genes
drug_targets <- drug_data %>%
  filter(toupper(GENE) %in% hub_genes, ORGANISM == "Homo sapiens") %>%
  select(DRUG_NAME, GENE, TARGET_NAME, ACTION_TYPE) %>%
  distinct() %>%
  mutate(Gene_Standardized = toupper(GENE))

unique_drugs <- sort(unique(drug_targets$DRUG_NAME))

# Generate all drug combinations
all_drug_combos <- expand.grid(
  Drug_A = unique_drugs,
  Drug_B = unique_drugs,
  stringsAsFactors = FALSE
) %>%
  filter(Drug_A < Drug_B) %>%
  mutate(Combination = paste(Drug_A, "+", Drug_B),
         Combo_ID = row_number())

# Helper functions
get_drug_targets <- function(drug_name) {
  drug_targets %>%
    filter(DRUG_NAME == drug_name) %>%
    pull(Gene_Standardized) %>%
    unique()
}

get_pIC50 <- function(drug_name, gene) {
  gene_upper <- toupper(gene)
  pIC50 <- potency_lookup %>%
    filter(DRUG_NAME == drug_name, Gene_Standardized == gene_upper) %>%
    pull(pIC50)
  if (length(pIC50) == 0) return(NA)
  return(pIC50)
}

calculate_potency_multiplier <- function(drug_a, gene_i, drug_b, gene_j) {
  pIC50_i <- get_pIC50(drug_a, gene_i)
  pIC50_j <- get_pIC50(drug_b, gene_j)
  if (is.na(pIC50_i) || is.na(pIC50_j)) return(0)
  max_pIC50 <- max(pIC50_i, pIC50_j)
  potency_mult <- max_pIC50 / 9
  return(min(potency_mult, 1.0))
}

get_mean_pIC50 <- function(drug_name, targets) {
  targets_vec <- unique(toupper(targets))
  pIC50_values <- sapply(targets_vec, function(g) get_pIC50(drug_name, g))
  mean_pIC50 <- mean(pIC50_values[!is.na(pIC50_values)], na.rm = TRUE)
  return(if (is.nan(mean_pIC50)) 0 else mean_pIC50)
}

# Enrich combinations with targets and pathways
all_drug_combos <- all_drug_combos %>%
  mutate(
    Targets_A = sapply(Drug_A, function(d) paste(get_drug_targets(d), collapse = "|")),
    Targets_B = sapply(Drug_B, function(d) paste(get_drug_targets(d), collapse = "|")),
    n_targets_A = sapply(Drug_A, function(d) length(get_drug_targets(d))),
    n_targets_B = sapply(Drug_B, function(d) length(get_drug_targets(d))),
    n_targets_combined = n_targets_A + n_targets_B
  )

# Calculate correlation statistics and complementarity
get_corr_stats <- function(genes_A, genes_B, correlations_df) {
  genes_A <- unique(toupper(genes_A))
  genes_B <- unique(toupper(genes_B))

  if (length(genes_A) == 0 || length(genes_B) == 0) {
    return(data.frame(Mean_Pos_r = 0, Max_Pos_r = 0, Mean_Neg_r = 0,
                      Sum_Neg_Abs_r = 0, Mean_Abs_r = 0))
  }

  pair_df <- expand.grid(Gene1 = genes_A, Gene2 = genes_B, stringsAsFactors = FALSE) %>%
    mutate(pair_key = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")) %>%
    left_join(correlations_df, by = "pair_key") %>%
    filter(!is.na(r))

  if (nrow(pair_df) == 0) {
    return(data.frame(Mean_Pos_r = 0, Max_Pos_r = 0, Mean_Neg_r = 0,
                      Sum_Neg_Abs_r = 0, Mean_Abs_r = 0))
  }

  pos_r <- pair_df$r[pair_df$r > 0]
  neg_r <- pair_df$r[pair_df$r < 0]

  return(data.frame(
    Mean_Pos_r = if (length(pos_r) > 0) mean(pos_r) else 0,
    Max_Pos_r = if (length(pos_r) > 0) max(pos_r) else 0,
    Mean_Neg_r = if (length(neg_r) > 0) mean(neg_r) else 0,
    Sum_Neg_Abs_r = if (length(neg_r) > 0) sum(abs(neg_r)) else 0,
    Mean_Abs_r = mean(abs(pair_df$r))
  ))
}

all_drug_combos <- all_drug_combos %>%
  rowwise() %>%
  mutate(
    targets_A_set = list(if (nzchar(Targets_A)) unlist(strsplit(Targets_A, "\\|")) else character(0)),
    targets_B_set = list(if (nzchar(Targets_B)) unlist(strsplit(Targets_B, "\\|")) else character(0)),
    intersection_size = length(intersect(targets_A_set, targets_B_set)),
    union_size = length(union(targets_A_set, targets_B_set)),
    Complementarity_C = if (union_size > 0) 1 - (intersection_size / union_size) else 0,
    Balance_B = if ((n_targets_A + n_targets_B) > 0) {
      1 - abs(n_targets_A - n_targets_B) / (n_targets_A + n_targets_B)
    } else {
      0
    },
    corr_stats = list(get_corr_stats(targets_A_set, targets_B_set, correlations_slim))
  ) %>%
  unnest_wider(corr_stats) %>%
  ungroup()

# Network coverage normalization
mean_abs_r_min <- min(all_drug_combos$Mean_Abs_r)
mean_abs_r_max <- max(all_drug_combos$Mean_Abs_r)

all_drug_combos <- all_drug_combos %>%
  mutate(Network_Coverage_V = if (mean_abs_r_max > mean_abs_r_min) {
    (Mean_Abs_r - mean_abs_r_min) / (mean_abs_r_max - mean_abs_r_min)
  } else {
    0
  })

# Calculate potency component
all_drug_combos <- all_drug_combos %>%
  rowwise() %>%
  mutate(
    Mean_pIC50_A = get_mean_pIC50(Drug_A, targets_A_set),
    Mean_pIC50_B = get_mean_pIC50(Drug_B, targets_B_set),
    Mean_pIC50 = mean(c(Mean_pIC50_A, Mean_pIC50_B), na.rm = TRUE),
    Potency_P = min(max(Mean_pIC50 - 4, 0), 5) / 5
  ) %>%
  ungroup()

# Synergy bonus with potency multiplier
synergy_bonus_with_potency <- function(drug_a, targets_a, drug_b, targets_b, correlations_df) {
  targets_a <- unique(toupper(targets_a))
  targets_b <- unique(toupper(targets_b))
  synergy_total <- 0

  for (i in seq_along(targets_a)) {
    for (j in seq_along(targets_b)) {
      gene_i <- targets_a[i]
      gene_j <- targets_b[j]
      pair_key <- paste(pmin(gene_i, gene_j), pmax(gene_i, gene_j), sep = "_")

      r_ij <- correlations_df %>%
        filter(pair_key == !!pair_key) %>%
        pull(r)

      if (length(r_ij) == 0 || is.na(r_ij)) next

      if (r_ij > 0) {
        potency_mult <- calculate_potency_multiplier(drug_a, gene_i, drug_b, gene_j)
        synergy_ij <- r_ij * potency_mult * 1
        synergy_total <- synergy_total + synergy_ij
      }
    }
  }

  return(synergy_total)
}

# Antagonism penalty with potency multiplier
antagonism_penalty_with_potency <- function(drug_a, targets_a, drug_b, targets_b, correlations_df) {
  targets_a <- unique(toupper(targets_a))
  targets_b <- unique(toupper(targets_b))
  penalty_total <- 0

  for (i in seq_along(targets_a)) {
    for (j in seq_along(targets_b)) {
      gene_i <- targets_a[i]
      gene_j <- targets_b[j]
      pair_key <- paste(pmin(gene_i, gene_j), pmax(gene_i, gene_j), sep = "_")

      r_ij <- correlations_df %>%
        filter(pair_key == !!pair_key) %>%
        pull(r)

      if (length(r_ij) == 0 || is.na(r_ij)) next

      if (r_ij < 0) {
        potency_mult <- calculate_potency_multiplier(drug_a, gene_i, drug_b, gene_j)
        penalty_ij <- abs(r_ij) * potency_mult * 1
        penalty_total <- penalty_total + penalty_ij
      }
    }
  }

  return(penalty_total)
}

all_drug_combos <- all_drug_combos %>%
  rowwise() %>%
  mutate(
    Synergy_Bonus = synergy_bonus_with_potency(Drug_A, targets_A_set, Drug_B, targets_B_set, correlations_slim),
    Antagonism_Penalty = antagonism_penalty_with_potency(Drug_A, targets_A_set, Drug_B, targets_B_set, correlations_slim)
  ) %>%
  ungroup()

# Final MLSS v4.0 calculation
all_drug_combos <- all_drug_combos %>%
  mutate(
    Base_Score = (0.45 * Complementarity_C +
                  0.30 * Balance_B +
                  0.10 * Network_Coverage_V +
                  0.15 * Potency_P) * 9,
    MLSS_v4.0 = Base_Score + Synergy_Bonus - Antagonism_Penalty
  ) %>%
  arrange(desc(MLSS_v4.0)) %>%
  mutate(Rank = row_number())

top_15_combos <- all_drug_combos %>% head(15)

# Visualization: Top 15 component breakdown
npg_mlss_colors <- c(
  "Complementarity (45%)" = npg_colors[1],
  "Balance (30%)" = npg_colors[4],
  "Coverage (10%)" = npg_colors[3],
  "Potency (15%)" = npg_colors[2],
  "Synergy Bonus" = npg_colors[5],
  "Antagonism Penalty" = npg_colors[6]
)

component_data <- top_15_combos %>%
  arrange(desc(MLSS_v4.0)) %>%
  mutate(
    C_Contribution = Complementarity_C * 0.45 * 9,
    B_Contribution = Balance_B * 0.30 * 9,
    V_Contribution = Network_Coverage_V * 0.10 * 9,
    P_Contribution = Potency_P * 0.15 * 9,
    Synergy_Component = Synergy_Bonus,
    Penalty_Component = Antagonism_Penalty,
    Combination = factor(Combination, levels = rev(Combination))
  ) %>%
  select(Combination, C_Contribution, B_Contribution, V_Contribution, P_Contribution,
         Synergy_Component, Penalty_Component, MLSS_v4.0) %>%
  pivot_longer(cols = c(C_Contribution, B_Contribution, V_Contribution, P_Contribution,
                        Synergy_Component, Penalty_Component),
               names_to = "Component", values_to = "Score") %>%
  mutate(Component = factor(Component,
                            levels = c("C_Contribution", "B_Contribution", "V_Contribution",
                                       "P_Contribution", "Synergy_Component", "Penalty_Component"),
                            labels = c("Complementarity (45%)", "Balance (30%)",
                                       "Coverage (10%)", "Potency (15%)",
                                       "Synergy Bonus", "Antagonism Penalty")))

mlss_labels <- component_data %>%
  group_by(Combination, MLSS_v4.0) %>%
  slice(1) %>%
  ungroup()

p <- ggplot(component_data, aes(x = Combination, y = Score, fill = Component)) +
  geom_bar(stat = "identity", width = 0.85, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.5, linetype = "solid") +
  geom_text(data = mlss_labels, aes(x = Combination, y = MLSS_v4.0,
                                     label = sprintf("%.2f", MLSS_v4.0)),
            hjust = -0.3, size = 3.5, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(name = "Score\nComponent", values = npg_mlss_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  guides(fill = guide_legend(nrow = 3)) +
  coord_flip() +
  labs(title = "MLSS Component Breakdown of Top 15 Drug Combinations",
       x = NULL, y = "Weighted Score Contribution") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 1, color = "black", family = "Arial"),
    axis.text.y = element_text(size = 12, color = "black", family = "Arial"),
    axis.text.x = element_text(size = 12, color = "black", family = "Arial"),
    axis.title.x = element_text(size = 12, color = "black", family = "Arial"),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12, family = "Arial", hjust = 1, face = "bold"),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.box = "vertical"
  )

ggsave("fig6c_mlss_top15.png", p, width = 7.7, height = 6, units = "cm",
       scale = 2.3, dpi = 600, bg = "white")
ggsave("fig6c_mlss_top15.svg", p, width = 7.7, height = 6, units = "cm",
       scale = 2.3, dpi = 600, bg = "white")

# Export results
write.csv(all_drug_combos %>%
            select(Rank, Combination, Drug_A, Drug_B,
                   n_targets_combined, Mean_pIC50,
                   Complementarity_C, Balance_B, Network_Coverage_V, Potency_P,
                   Mean_Pos_r, Mean_Neg_r,
                   Base_Score, Synergy_Bonus, Antagonism_Penalty, MLSS_v4.0),
          "fig6c_mlss_all_combinations.csv", row.names = FALSE)

write.csv(top_15_combos %>%
            select(Rank, Combination, Drug_A, Drug_B,
                   n_targets_combined, Mean_pIC50,
                   Complementarity_C, Balance_B, Network_Coverage_V, Potency_P,
                   Mean_Pos_r, Mean_Neg_r,
                   Base_Score, Synergy_Bonus, Antagonism_Penalty, MLSS_v4.0),
          "fig6c_mlss_top15.csv", row.names = FALSE)
