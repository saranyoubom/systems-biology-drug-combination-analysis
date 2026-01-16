# Gene Interaction Network Analysis

Comprehensive network analysis integrating correlation data with STRING database validation.

## Description

Four-panel figure showing:
- **Panel A**: Gene interaction network with novel discoveries highlighted
- **Panel B**: Pie chart of STRING database coverage
- **Panel C**: Correlation strength distributions by category
- **Panel D**: Top 15 novel gene pair discoveries

## Requirements

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "ggrepel", "igraph", "ggraph", "patchwork", "ggsci"))
```

## Usage

```r
source("fig5_network_analysis.R")
```

## Dependencies

Requires outputs from previous analyses:
- `string_interactions.csv` - STRING database interactions
- `Fig_4d_consensus_pairs.csv` - Consensus pairs from Fig 4d
- `Fig_4e_Hes1_significant_correlations_only.csv` - Correlations from Fig 4e

## Output

**Figures** (all in `./figures_publication/`):
- `Fig5a_Network.png/svg` - Interaction network (600 dpi, 12 × 11.25 cm)
- `Fig5b_Coverage.png/svg` - STRING coverage pie chart
- `Fig5c_Distribution.png/svg` - Correlation distributions
- `Fig5d_TopNovel.png/svg` - Top novel discoveries

**Data exports** (in `./figures_publication/data_export/`):
- `Fig5B_Coverage_Summary.csv`
- `Fig5D_Top_Novel_Discoveries.csv`
- `Fig5_Summary_Statistics.csv`
- `Fig5_All_Classified_Pairs.csv`

## Gene Pair Classification

Pairs with |r| ≥ 0.7 classified as:

1. **Known Validated**: In STRING database + 4-method consensus
2. **Novel Discovery**: Not in STRING + 4-method consensus
3. **STRING Only**: In STRING but not in 4-method consensus
4. **Not Validated**: Neither STRING nor 4-method consensus

## Network Visualization

- Node size: Degree centrality
- Edge width: Correlation strength (|r|)
- Colors: NPG palette
  - Red: Novel discoveries
  - Green: STRING validated
  - Yellow: STRING only
  - Grey: Not validated

## Novel Discoveries

Gene pairs with strong correlations (|r| ≥ 0.7) validated across all 4 normalization methods but absent from STRING database, representing potential new biological relationships.
