# Normalization Methods Comparison

Comparison of four normalization approaches for gene correlation analysis.

## Description

Compares correlation patterns detected using different normalization strategies.

## Requirements

```r
install.packages(c("ggplot2", "readxl", "reshape2", "Hmisc", "ggsci"))
```

## Usage

```r
source("fig4c_normalization_methods.R")
```

## Input Data

**File**: `Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx`

## Output

**Figures**:
- `Fig_4c_Normalization_comparison.png` (600 dpi, 8.1 × 8.7 cm)
- `Fig_4c_Normalization_comparison.svg`

**Data exports** (4 normalization methods + summary):
- `Fig_4c_Hes1_normalized_data.csv`
- `Fig_4c_Kcnj11_normalized_data.csv`
- `Fig_4c_Ins1_normalized_data.csv`
- `Fig_4c_Raw_expression_data.csv`
- `Fig_4c_normalization_comparison_data.csv`

## Normalization Methods

1. **Hes1-normalized**: Expression ratios relative to Hes1
2. **Kcnj11-normalized**: Expression ratios relative to Kcnj11
3. **Ins1-normalized**: Expression ratios relative to Ins1
4. **Raw expression**: Absolute expression values

## Correlation Categories

- **Significant**: p < 0.05
- **Strong**: |r| > 0.7 and p < 0.05
- **Very Strong**: |r| > 0.8 and p < 0.05

## Notes

Normalized data exported for downstream analyses (Fig 4d, 4e, 5).
