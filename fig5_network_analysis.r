# Gene Interaction Network Analysis - Figure 5

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(igraph)
library(ggraph)
library(patchwork)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# Create output directories
dir.create("./figures_publication", showWarnings = FALSE, recursive = TRUE)
dir.create("./figures_publication/data_export", showWarnings = FALSE, recursive = TRUE)

# Helper function: create canonical pair names
create_pair_name <- function(gene1, gene2) {
  pairs <- mapply(function(g1, g2) {
    sorted <- sort(c(g1, g2))
    paste(sorted, collapse = "-")
  }, gene1, gene2, SIMPLIFY = TRUE)
  return(pairs)
}

# Load STRING interactions
string_data <- read.csv("string_interactions.csv", stringsAsFactors = FALSE) %>%
  mutate(Pair = create_pair_name(node1, node2)) %>%
  select(Pair, node1, node2, combined_score) %>%
  rename(STRING_Score = combined_score) %>%
  distinct(Pair, .keep_all = TRUE)

# Load consensus pairs from Fig 4d
validated_pairs_raw <- read.csv("Fig_4d_consensus_pairs.csv", stringsAsFactors = FALSE)
validated_pairs <- validated_pairs_raw %>%
  mutate(
    Gene1 = sapply(strsplit(Consensus_Pairs, "_"), `[`, 1),
    Gene2 = sapply(strsplit(Consensus_Pairs, "_"), `[`, 2),
    Pair = create_pair_name(Gene1, Gene2)
  ) %>%
  select(Pair, Gene1, Gene2) %>%
  distinct(Pair, .keep_all = TRUE)

# Load Hes1-normalized correlations from Fig 4e
correlation_data <- read.csv("Fig_4e_Hes1_significant_correlations_only.csv",
                             stringsAsFactors = FALSE, row.names = 1) %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene1") %>%
  pivot_longer(-Gene1, names_to = "Gene2", values_to = "Correlation") %>%
  filter(!is.na(Correlation)) %>%
  filter(Gene1 < Gene2) %>%
  mutate(
    Pair = create_pair_name(Gene1, Gene2),
    r = Correlation
  ) %>%
  select(Pair, Gene1, Gene2, r)

# Classify pairs
string_set <- string_data$Pair
validated_set <- validated_pairs$Pair

all_pairs <- correlation_data %>%
  filter(abs(r) >= 0.7) %>%
  mutate(
    In_STRING = Pair %in% string_set,
    Four_Method_Validated = Pair %in% validated_set,
    Category = case_when(
      In_STRING & Four_Method_Validated ~ "Known Validated",
      !In_STRING & Four_Method_Validated ~ "Novel Discovery",
      In_STRING & !Four_Method_Validated ~ "STRING Only",
      TRUE ~ "Not Validated"
    )
  )

novel_discoveries <- all_pairs %>%
  filter(Category == "Novel Discovery") %>%
  arrange(desc(abs(r)))

# Summary statistics
summary_stats <- all_pairs %>%
  group_by(Category) %>%
  summarise(
    Count = n(),
    Mean_r = mean(abs(r)),
    Median_r = median(abs(r)),
    .groups = "drop"
  ) %>%
  arrange(desc(Count))

# Custom theme
theme_fig5 <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, family = "Arial"),
      axis.title = element_text(size = 10, face = "bold", color = "black", family = "Arial"),
      axis.text = element_text(size = 10, color = "black", family = "Arial"),
      legend.text = element_text(size = 10, family = "Arial"),
      legend.title = element_text(size = 10, face = "bold", family = "Arial"),
      panel.grid.major = element_line(color = "gray70", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

# Panel A: Gene Interaction Network
string_validated <- all_pairs %>%
  filter(Category == "Known Validated") %>%
  arrange(desc(abs(r))) %>%
  head(30)

novel_top <- novel_discoveries %>% head(20)

remaining_needed <- 50 - nrow(string_validated) - nrow(novel_top)
if (remaining_needed > 0) {
  supplement <- all_pairs %>%
    filter(Category %in% c("STRING Only", "Not Validated")) %>%
    arrange(desc(abs(r))) %>%
    head(remaining_needed)
} else {
  supplement <- data.frame()
}

# Create edge list
edges_list <- list()
if (nrow(string_validated) > 0) {
  edges_list[[1]] <- data.frame(
    from = string_validated$Gene1,
    to = string_validated$Gene2,
    weight = abs(string_validated$r),
    EdgeType = "STRING Validated",
    stringsAsFactors = FALSE
  )
}
if (nrow(novel_top) > 0) {
  edges_list[[2]] <- data.frame(
    from = novel_top$Gene1,
    to = novel_top$Gene2,
    weight = abs(novel_top$r),
    EdgeType = "Novel Discovery",
    stringsAsFactors = FALSE
  )
}
if (nrow(supplement) > 0) {
  edges_list[[3]] <- data.frame(
    from = supplement$Gene1,
    to = supplement$Gene2,
    weight = abs(supplement$r),
    EdgeType = supplement$Category,
    stringsAsFactors = FALSE
  )
}

edges <- do.call(rbind, edges_list)

# Create graph
g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, weights = E(g)$weight, directed = FALSE)

novel_genes <- unique(c(novel_top$Gene1, novel_top$Gene2))
V(g)$node_type <- ifelse(V(g)$name %in% novel_genes, "Novel", "Known")

edge_colors <- case_when(
  E(g)$EdgeType == "Novel Discovery" ~ npg_colors[1],
  E(g)$EdgeType == "STRING Validated" ~ npg_colors[4],
  E(g)$EdgeType == "STRING Only" ~ npg_colors[2],
  TRUE ~ "gray70"
)

E(g)$color <- edge_colors
E(g)$width <- E(g)$weight * 1.5

set.seed(42)
layout_net <- layout_with_fr(g, weights = E(g)$weight)

panel_a <- ggraph(g, layout = layout_net) +
  geom_edge_link(aes(edge_color = EdgeType, edge_width = weight),
                 alpha = 0.6, show.legend = TRUE) +
  geom_node_point(aes(size = degree, color = node_type, fill = node_type),
                  alpha = 0.8, shape = 21) +
  geom_node_text(aes(label = name, color = node_type),
                 size = 2.5, repel = TRUE, max.overlaps = 20,
                 fontface = "italic", family = "Arial") +
  scale_edge_color_manual(
    name = "Interaction",
    values = c("STRING Validated" = npg_colors[4],
               "Novel Discovery" = npg_colors[1],
               "STRING Only" = npg_colors[2],
               "Not Validated" = "gray70")
  ) +
  scale_edge_width_continuous(name = "Abs. Correlation (|r|)", range = c(0.3, 2)) +
  scale_color_manual(
    name = "Gene Status",
    values = c("Novel" = npg_colors[1], "Known" = npg_colors[4])
  ) +
  scale_fill_manual(
    name = "Gene Status",
    values = c("Novel" = npg_colors[1], "Known" = npg_colors[4])
  ) +
  scale_size_continuous(name = "Degree", range = c(3, 10)) +
  labs(title = "Gene Interaction Network") +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = "Arial"),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8, family = "Arial"),
    legend.title = element_text(size = 8, face = "bold", family = "Arial")
  )

ggsave("./figures_publication/Fig5a_Network.png",
       plot = panel_a, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

ggsave("./figures_publication/Fig5a_Network.svg",
       plot = panel_a, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

# Panel B: STRING Coverage Distribution
panel_b_data <- all_pairs %>%
  group_by(Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(
    Percentage = 100 * Count / sum(Count),
    Category = factor(Category,
                      levels = c("Known Validated", "Novel Discovery",
                                 "STRING Only", "Not Validated"))
  )

panel_b <- ggplot(panel_b_data, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
            position = position_stack(vjust = 0.5),
            size = 3, fontface = "bold", color = "white", family = "Arial") +
  scale_fill_manual(
    name = NULL,
    values = c("Known Validated" = npg_colors[4],
               "Novel Discovery" = npg_colors[1],
               "STRING Only" = npg_colors[2],
               "Not Validated" = "gray70")
  ) +
  labs(title = "STRING Database Coverage") +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = "Arial"),
    legend.position = "bottom",
    legend.text = element_text(size = 8, family = "Arial")
  )

ggsave("./figures_publication/Fig5b_Coverage.png",
       plot = panel_b, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

ggsave("./figures_publication/Fig5b_Coverage.svg",
       plot = panel_b, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

write.csv(panel_b_data,
          "./figures_publication/data_export/Fig5B_Coverage_Summary.csv",
          row.names = FALSE)

# Panel C: Correlation Distribution
panel_c <- ggplot(all_pairs, aes(x = abs(r), fill = Category)) +
  geom_histogram(bins = 15, alpha = 0.7, position = "identity",
                 color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "black",
             linewidth = 0.4, alpha = 0.5) +
  facet_wrap(~factor(Category,
                     levels = c("Known Validated", "Novel Discovery",
                                "STRING Only", "Not Validated")),
             ncol = 2, scales = "free_y") +
  scale_fill_manual(
    values = c("Known Validated" = npg_colors[4],
               "Novel Discovery" = npg_colors[1],
               "STRING Only" = npg_colors[2],
               "Not Validated" = "gray70")
  ) +
  labs(title = "Correlation Strength Distribution",
       x = "Absolute Correlation (|r|)",
       y = "Frequency") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10, family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial", color = "black"),
    axis.title = element_text(size = 10, family = "Arial", color = "black")
  )

ggsave("./figures_publication/Fig5c_Distribution.png",
       plot = panel_c, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

ggsave("./figures_publication/Fig5c_Distribution.svg",
       plot = panel_c, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

# Panel D: Top Novel Discoveries
top_novel <- novel_discoveries %>%
  head(15) %>%
  mutate(
    Pair_Label = paste(Gene1, Gene2, sep = "-"),
    Rank = row_number()
  ) %>%
  arrange(desc(abs(r)))

panel_d <- ggplot(top_novel, aes(x = reorder(Pair_Label, abs(r)), y = abs(r))) +
  geom_col(fill = npg_colors[1], alpha = 0.8) +
  geom_text(aes(label = sprintf("%.2f", r)),
            hjust = -0.2, size = 2.5, family = "Arial") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Top 15 Novel Gene Pair Discoveries",
       x = NULL,
       y = "Absolute Correlation (|r|)") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = "Arial"),
    axis.text.y = element_text(size = 8, face = "italic", family = "Arial"),
    axis.text.x = element_text(size = 10, family = "Arial"),
    axis.title = element_text(size = 10, family = "Arial"),
    panel.grid.major.y = element_blank()
  )

ggsave("./figures_publication/Fig5d_TopNovel.png",
       plot = panel_d, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

ggsave("./figures_publication/Fig5d_TopNovel.svg",
       plot = panel_d, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

write.csv(top_novel,
          "./figures_publication/data_export/Fig5D_Top_Novel_Discoveries.csv",
          row.names = FALSE)

# Export summary
write.csv(summary_stats,
          "./figures_publication/data_export/Fig5_Summary_Statistics.csv",
          row.names = FALSE)

write.csv(all_pairs,
          "./figures_publication/data_export/Fig5_All_Classified_Pairs.csv",
          row.names = FALSE)
