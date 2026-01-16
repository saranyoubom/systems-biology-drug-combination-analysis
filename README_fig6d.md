# MLSS Validation Analysis - Figure 6d & Supplementary S2

Sensitivity and robustness testing of MLSS v4.0 scoring algorithm.

## Description

Validates MLSS through:
1. **Weight sensitivity**: Tests 8 alternative weight schemes
2. **Ablation study**: Leave-one-out component analysis
3. **Ranking correlation**: Spearman correlation across scenarios
4. **Stability metrics**: Top 5 overlap consistency

## Requirements

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "ggsci", "reshape2", "gridExtra"))
```

## Usage

```r
source("fig6d_mlss_validation.R")
```

## Dependencies

- `fig6c_mlss_all_combinations.csv` - MLSS output from Fig 6c

## Output

**Figures**:
- `figs2a_mlss_sensitivity_top15.png/svg` - Weight sensitivity heatmap
- `figs2b_component_importance_top15.png/svg` - Ablation analysis

**Data exports**:
- `fig6d_ranking_correlations.csv` - Spearman correlations between scenarios
- `fig6d_ablation_impact.csv` - Component removal impact metrics

## Weight Scenarios Tested

1. **MLSS (Standard)**: 45:30:10:15 (C:B:V:P)
2. **C-Heavy (60:25:15)**: Emphasizes complementarity
3. **C-Heavy (55:20:25)**: Moderate complementarity emphasis
4. **B-Heavy (35:50:15)**: Emphasizes balance
5. **V-Heavy (35:30:35)**: Emphasizes network coverage
6. **P-Heavy (30:20:10:40)**: Emphasizes potency
7. **Balanced (40:40:10:10)**: Equal C and B
8. **Network-Focus (30:20:40:10)**: Emphasizes network coverage

## Ablation Conditions

Tests MLSS with one component removed:
- **Without C**: C=0, renormalize B:V:P
- **Without B**: B=0, renormalize C:V:P
- **Without V**: V=0, renormalize C:B:P
- **Without P**: P=0, renormalize C:B:V
- **Base Only**: Remove synergy/antagonism

## Metrics

### Sensitivity
- **Spearman correlation**: Ranking consistency (0-1)
- **Normalized scores**: Heatmap visualization

### Ablation Impact
- **Avg Rank Shift**: Mean position change (0-15)
- **Avg Score Drop**: Mean MLSS reduction
- **Top 5 Stability**: Proportion maintained in top 5 (0-1)

## Interpretation

### High Robustness
- Spearman ρ > 0.95 (excellent)
- Spearman ρ > 0.85 (good)
- Rank shift < 2 positions
- Score drop < 1 unit

### Component Importance
Higher values indicate more critical components:
- **Complementarity**: Typically highest impact
- **Balance**: Moderate impact
- **Synergy/Antagonism**: Variable impact
- **Coverage, Potency**: Lower but non-negligible

## Visualization Details

### Figure S2a: Sensitivity Heatmap
- **Rows**: Top 15 drug combinations
- **Columns**: 8 weight scenarios
- **Color**: Normalized MLSS (purple=low, red=high)

### Figure S2b: Component Importance
- **Left panel**: Average rank shift (positions)
- **Right panel**: Average score drop (MLSS units)
- **Colors**: NPG palette per component
