# Consensus Hub Gene Network - Figure 6b

Network analysis of consensus gene pairs with correlation strength visualization using NPG color palette.

## Description

Circular network graph showing:
- **Nodes**: Genes from 4-method consensus (size = degree centrality)
- **Edges**: Correlations (width = |r|, linetype = sign)
- **Colors**: Functional categories (Calcium, Incretin, Wnt, Metabolic)

## Requirements

```r
install.packages(c("tidyverse", "igraph", "ggraph", "scales", "ggsci"))
```

## Usage

```r
source("fig6b_consensus_hub_network.R")
```

## Dependencies

- `Fig_4d_consensus_pairs.csv` - Consensus pairs from Fig 4d
- `Fig_4e_Hes1_significant_correlations_only.csv` - Correlations from Fig 4e

## Output

**Figures**:
- `fig6b_consensus_hub_network.png/svg` (600 dpi, 17.71 × 13.8 cm)

**Data exports**:
- `fig6b_gene_centrality_summary.csv` - Node metrics (degree, betweenness, closeness, eigenvector centrality)
- `fig6b_edge_correlation_statistics.csv` - Edge properties (correlation, strength, type)

## Network Features

### Node Properties
- **Size**: Degree centrality (number of connections)
- **Color**: Functional category
  - Dark blue: Calcium Signaling
  - Red/Orange: Incretin Signaling
  - Green: Wnt Signaling
  - Light blue: Metabolic Regulation
- **Label**: Bold italic for hubs (≥6 connections), italic for others

### Edge Properties
- **Width**: Absolute correlation strength (0.5-3.0 pt)
- **Linetype**: Solid (positive r), Dotted (negative r)
- **Color**: Grey70 (uniform)

### Hub Definition
Genes with degree ≥ 6 connections are classified as hubs

## Network Metrics

The script calculates:
- Degree centrality (number of direct connections)
- Betweenness centrality (control over information flow)
- Closeness centrality (proximity to all nodes)
- Eigenvector centrality (influence based on connections)

## Layout

Fruchterman-Reingold circular layout optimizes node positions for clarity
