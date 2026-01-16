# MLSS v4.0 - Multi-Layered Synergy Score - Figure 6c

Comprehensive drug combination scoring algorithm integrating complementarity, balance, network coverage, potency, synergy, and antagonism.

## Description

MLSS v4.0 formula:
```
MLSS = [C×0.45 + B×0.30 + V×0.10 + P×0.15] × 9
       + Σ(r_pos × Potency_Multiplier) × 1
       - Σ(|r_neg| × Potency_Multiplier) × 1
```

Where:
- **C** = Complementarity (target non-overlap)
- **B** = Balance (target count similarity)
- **V** = Network Coverage (correlation density)
- **P** = Potency (pIC50-derived binding strength)
- **Synergy** = Positive correlations weighted by potency
- **Antagonism** = Negative correlations weighted by potency

## Requirements

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "ggsci", "tibble"))
```

## Usage

```r
source("fig6c_mlss_drug_combinations.R")
```

## Dependencies

- `fig6b_gene_centrality_summary.csv` - Hub genes (from Fig 6b)
- `Fig_4e_Hes1_significant_correlations_only.csv` - Correlations (from Fig 4e)
- `drug.target.interaction.csv` - Drug-target database with pIC50 values

## Output

**Figures**:
- `fig6c_mlss_top15.png/svg` (600 dpi, 17.71 × 13.8 cm)

**Data exports**:
- `fig6c_mlss_all_combinations.csv` - All drug pair scores
- `fig6c_mlss_top15.csv` - Top 15 combinations with detailed metrics

## Scoring Components

### Base Score Components (0-9 scale)
1. **Complementarity (45%)**: Measures target diversity
   - Range: 0 (complete overlap) to 1 (no overlap)
2. **Balance (30%)**: Penalizes unequal target counts
   - Range: 0 (highly imbalanced) to 1 (equal targets)
3. **Network Coverage (10%)**: Normalized correlation density
   - Range: 0 (isolated) to 1 (highly connected)
4. **Potency (15%)**: pIC50-based binding strength
   - Range: 0 (weak, pIC50 ≤4) to 1 (strong, pIC50 ≥9)

### Network Modifiers
5. **Synergy Bonus**: Sum of positive gene correlations × potency multipliers
6. **Antagonism Penalty**: Sum of negative gene correlations × potency multipliers

## Potency Multiplier

```
Potency_Multiplier = min(max(pIC50_i, pIC50_j) / 9, 1.0)
```

Scales synergy/antagonism by strongest binding interaction

## Algorithm Flow

1. Load hub genes (degree ≥ 6)
2. Extract drug-target interactions for hub genes
3. Generate all pairwise drug combinations
4. Calculate base components (C, B, V, P)
5. Compute synergy bonus and antagonism penalty
6. Aggregate into final MLSS score
7. Rank combinations

## Interpretation

- **High MLSS (>6)**: Strong combination potential
  - High complementarity, balanced targets, potent binding
  - Positive gene synergy, minimal antagonism
- **Low MLSS (<3)**: Poor combination
  - Overlapping targets, imbalanced, weak potency
  - Antagonistic gene interactions

## Visualization

Stacked horizontal bar chart showing component contributions for top 15:
- Bottom (base): C, B, V, P contributions
- Top: Synergy bonus
- Red (negative): Antagonism penalty
- Total bar height = MLSS score

## Color Scheme

NPG palette:
- Red: Complementarity
- Dark blue: Balance
- Light blue: Coverage
- Green: Potency
- Yellow: Synergy
- Purple: Antagonism
