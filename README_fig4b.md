# Gene Correlation Matrix

Correlation matrix showing significant Pearson correlations (p < 0.05) between genes.

## Description

Visualizes pairwise gene correlations with NPG color gradient. Non-significant correlations shown as grey.

## Requirements

```r
install.packages(c("ggplot2", "readxl", "ggcorrplot", "dplyr", "ggsci"))
```

## Usage

```r
source("fig4b_correlation_matrix.R")
```

## Input Data

**File**: `Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx`

**Format**: Gene expression matrix with genes as columns, samples as rows

## Output

**Figures**:
- `Fig_4b_Correlation_matrix.png` (600 dpi, 9.16 × 9.16 cm)
- `Fig_4b_Correlation_matrix.svg`

**Data exports**:
- `Fig_4b_full_correlation_matrix.csv` - All correlations
- `Fig_4b_pvalue_matrix.csv` - P-values for each correlation
- `Fig_4b_significant_correlations.csv` - Significant pairs only

## Analysis Method

- **Correlation**: Pearson's r
- **Significance**: p < 0.05 threshold
- **Correction**: None (exploratory analysis)

## Visualization

- Blue: Negative correlations
- White: Zero correlation
- Red: Positive correlations
- Grey: Non-significant (p ≥ 0.05)
