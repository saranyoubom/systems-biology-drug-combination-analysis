# Four-way Venn Diagram - Figure 4d

rm(list = ls())

library(VennDiagram)
library(grid)
library(Hmisc)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# Load normalized data
hes1_norm <- read.csv("Fig_4c_Hes1_normalized_data.csv")
kcnj11_norm <- read.csv("Fig_4c_Kcnj11_normalized_data.csv")
ins1_norm <- read.csv("Fig_4c_Ins1_normalized_data.csv")
raw_data <- read.csv("Fig_4c_Raw_expression_data.csv")

# Function to extract significant pairs
get_pairs <- function(data, method_name) {
  data_clean <- na.omit(data)
  cor_results <- rcorr(as.matrix(data_clean), type = "pearson")
  cor_matrix <- cor_results$r
  p_matrix <- cor_results$P

  genes <- colnames(data_clean)
  n_genes <- length(genes)
  pairs <- character()

  for (i in 1:(n_genes-1)) {
    for (j in (i+1):n_genes) {
      if (!is.na(p_matrix[i, j]) && p_matrix[i, j] < 0.05) {
        pairs <- c(pairs, paste(genes[i], genes[j], sep = "_"))
      }
    }
  }
  return(pairs)
}

# Extract pairs from each method
pairs_list <- list(
  "Hes1-\nnormalized" = get_pairs(hes1_norm, "Hes1"),
  "Kcnj11-\nnormalized" = get_pairs(kcnj11_norm, "Kcnj11"),
  "Ins1-\nnormalized" = get_pairs(ins1_norm, "Ins1"),
  "Raw\nexpression" = get_pairs(raw_data, "Raw")
)

# Identify consensus pairs
consensus_pairs <- Reduce(intersect, pairs_list)

# Export
write.csv(data.frame(Pairs = pairs_list[[1]]),
          "Fig_4d_Hes1_normalized_pairs.csv", row.names = FALSE)
write.csv(data.frame(Pairs = pairs_list[[2]]),
          "Fig_4d_Kcnj11_normalized_pairs.csv", row.names = FALSE)
write.csv(data.frame(Pairs = pairs_list[[3]]),
          "Fig_4d_Ins1_normalized_pairs.csv", row.names = FALSE)
write.csv(data.frame(Pairs = pairs_list[[4]]),
          "Fig_4d_Raw_expression_pairs.csv", row.names = FALSE)
write.csv(data.frame(Consensus_Pairs = consensus_pairs),
          "Fig_4d_consensus_pairs.csv", row.names = FALSE)

# Create Venn diagram
png("Fig_4d_Venn_diagram.png", width = 9.16, height = 9.16, units = "cm", res = 600)
venn.plot <- venn.diagram(
  x = pairs_list,
  category.names = names(pairs_list),
  filename = NULL,
  fill = c(npg_colors[4], npg_colors[2], npg_colors[3], npg_colors[7]),
  col = c(npg_colors[4], npg_colors[2], npg_colors[3], npg_colors[7]),
  lwd = 1.5,
  alpha = 0.5,
  cex = 0.8,
  fontface = "bold",
  fontfamily = "Arial",
  cat.cex = 0.7,
  cat.fontface = "bold",
  cat.fontfamily = "Arial",
  cat.default.pos = "outer",
  cat.dist = c(0.08, 0.08, 0.08, 0.08)
)
grid.draw(venn.plot)
dev.off()

svg("Fig_4d_Venn_diagram.svg", width = 9.16/2.54, height = 9.16/2.54)
grid.draw(venn.plot)
dev.off()
