# Consensus Gene Pairs - Venn Diagram

Four-way Venn diagram showing overlap of significant gene correlations across normalization methods.

## Description

Identifies gene pairs with significant correlations across all four normalization methods (consensus pairs).

## Requirements

```r
install.packages(c("VennDiagram", "grid", "Hmisc", "ggsci"))
```

## Usage

```r
source("fig4d_venn_diagram.R")
```

## Dependencies

Requires normalized data from `fig4c_normalization_methods.R`:
- `Fig_4c_Hes1_normalized_data.csv`
- `Fig_4c_Kcnj11_normalized_data.csv`
- `Fig_4c_Ins1_normalized_data.csv`
- `Fig_4c_Raw_expression_data.csv`

## Output

**Figures**:
- `Fig_4d_Venn_diagram.png` (600 dpi, 9.16 × 9.16 cm)
- `Fig_4d_Venn_diagram.svg`

**Data exports**:
- `Fig_4d_Hes1_normalized_pairs.csv`
- `Fig_4d_Kcnj11_normalized_pairs.csv`
- `Fig_4d_Ins1_normalized_pairs.csv`
- `Fig_4d_Raw_expression_pairs.csv`
- `Fig_4d_consensus_pairs.csv` - Pairs significant in ALL methods
- `Fig_4d_venn_summary.csv`

## Consensus Pairs

Gene pairs significant (p < 0.05) across all four methods represent high-confidence correlations for downstream analysis.

## Notes

Consensus pairs list used in Figure 4e and Figure 5 analyses.
