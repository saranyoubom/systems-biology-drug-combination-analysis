# Drug Combination Network - Figure 6e

Circular network visualization of top MLSS drug combinations.

## Description

Network graph showing drug-drug relationships:
- **Nodes**: Individual drugs (size = appearance frequency)
- **Edges**: Drug pairs from top MLSS combinations
- **Edge width**: MLSS score magnitude
- **Edge linetype**: Cross-pathway vs same-pathway
- **Node colors**: Primary functional pathway

## Requirements

```r
install.packages(c("tidyverse", "igraph", "ggraph", "scales", "ggsci"))
```

## Usage

```r
source("fig6e_drug_combination_network.R")
```

## Dependencies

- `fig6c_mlss_all_combinations.csv` - MLSS scores from Fig 6c

## Output

**Figures**:
- `fig6e_drug_combination_network.png/svg` (600 dpi, 17.71 × 13.8 cm)

**Data exports**:
- `fig6e_drug_node_statistics.csv` - Drug node properties
- `fig6e_drug_edge_statistics.csv` - Combination edge properties

## Network Properties

### Nodes (Drugs)
- **Size**: Network appearances (1 to max)
- **Color**: Primary pathway
  - Dark blue: Calcium Signaling
  - Red/Orange: Incretin Signaling  
  - Green: Wnt Signaling
  - Light blue: Metabolic Regulation
  - Purple: Multi-pathway
  - Grey: Unknown
- **Label**: Bold for hubs (≥4 appearances)

### Edges (Combinations)
- **Width**: MLSS score (0.5-3.5 pt scale)
  - Thin: Lower scores
  - Medium: Intermediate scores
  - Thick: Higher scores
- **Linetype**: 
  - Solid: Cross-pathway combination
  - Dotted: Same-pathway combination
- **Color**: Grey70 (uniform)

## Hub Definition

Drugs appearing ≥4 times in top 15 combinations are classified as hub drugs

## Layout

Circular layout positions nodes around a circle for balanced visualization

## Interpretation

### Hub Drugs
- Central to multiple high-scoring combinations
- Likely key candidates for repositioning
- Cross multiple pathways

### Edge Patterns
- **Thick solid edges**: Strong cross-pathway combinations
- **Thick dotted edges**: Strong within-pathway combinations
- **Thin edges**: Weaker combinations

## Customization

To analyze more/fewer combinations, change `top_n` parameter (default: 15)
