# Glucose Uptake Analysis

Measurement of glucose uptake using 2-NBDG fluorescence intensity across treatment groups.

## Description

Generates a boxplot comparing fold change in 2-NBDG (fluorescent glucose analog) uptake.

## Requirements

```r
install.packages(c("ggplot2", "data.table", "tidyr", "ggpubr", "dplyr", "ggsci", "readxl"))
```

## Usage

```r
source("fig2b_glucose_uptake.R")
```

## Input Data

**File**: `Fig_2b_Boxplot_Glucose_uptake.xlsx`

**Columns**:
- `Treatment`: CTRL, DAPT, DKK-1
- `Fold_change`: Fold change in 2-NBDG intensity
- `Assay`: Assay identifier

## Output

- `Fig_2b_Boxplot_Glucose_uptake.png` (600 dpi, 8 × 5.3 cm)
- `Fig_2b_Boxplot_Glucose_uptake.svg`

## Statistical Analysis

- **Normality test**: Shapiro-Wilk (W = 0.7301, p = 0.0017)
- **Overall comparison**: Kruskal-Wallis test (non-parametric)
- **Pairwise comparisons**: Wilcoxon test
- **Visualization**: Mean ± SEM annotations

## Notes

Non-parametric tests used due to non-normal distribution of data.
