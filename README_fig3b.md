# C-peptide at Physiological Glucose

C-peptide secretion measurement at 5.5 mM glucose (physiological fasting level).

## Description

Compares C-peptide release normalized to DNA content at physiological glucose concentration.

## Requirements

```r
install.packages(c("ggplot2", "data.table", "tidyr", "ggpubr", "dplyr", "ggsci", "readxl"))
```

## Usage

```r
source("fig3b_cpeptide_5.5mm.R")
```

## Input Data

**File**: `Fig_3b_Boxplot_C_peptide.xlsx`

**Columns**:
- `Treatment`: CTRL, DAPT, DKK-1
- `Glucose`: Glucose concentration categories
- `C_peptide`: C-peptide release (ng/µg DNA)

## Output

- `Fig_3b_Boxplot_C_peptide.png` (600 dpi, 8 × 6.5 cm)
- `Fig_3b_Boxplot_C_peptide.svg`

## Statistical Analysis

- **Normality test**: Shapiro-Wilk (W = 0.6920, p = 0.0007)
- **Overall comparison**: Kruskal-Wallis test (non-parametric)
- **Pairwise comparisons**: Wilcoxon test
- **Visualization**: Mean ± SEM annotations

## Notes

Data filtered to show only 5.5 mM glucose condition. C-peptide values normalized to total DNA content.
