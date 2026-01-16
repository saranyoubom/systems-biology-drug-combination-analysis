# Principal Component Analysis

PCA biplot showing treatment clustering based on mRNA expression profiles.

## Description

Performs dimensionality reduction on relative mRNA expression data to visualize treatment effects.

## Requirements

```r
install.packages(c("ggplot2", "data.table", "tidyr", "ggpubr", "dplyr", "ggsci", "readxl", "factoextra", "ggfortify"))
```

## Usage

```r
source("fig3f_pca_analysis.R")
```

## Input Data

**File**: `Fig_3f_PCA_Relative_mRNA_expression.xlsx`

**Columns**:
- `Group`: CTRL, DAPT, DKK-1
- Multiple gene expression columns (numeric values)

## Output

- `Fig_3f_PCA_Relative_mRNA_expression.png` (600 dpi, 7.6 × 5.3 cm)
- `Fig_3f_PCA_Relative_mRNA_expression.svg`

## Analysis Method

- **Scaling**: Data scaled prior to PCA
- **Confidence ellipses**: 95% confidence intervals
- **Biplot**: Both samples and variables displayed
- **Variable arrows**: Grey with repelling labels

## Visualization

- NPG color palette for treatment groups
- Ellipses show within-group variation
- Arrow length indicates variable contribution
- Arrow direction shows correlation with PCs
