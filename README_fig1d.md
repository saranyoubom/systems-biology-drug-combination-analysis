# mRNA Expression Heatmap

Heatmap visualization of fold change in mRNA expression for Notch and Wnt pathway genes.

## Description

Creates a tile heatmap with numerical annotations showing expression changes across treatments.

## Requirements

```r
install.packages(c("ggplot2", "pheatmap", "dplyr", "readxl", "scales", "ggsci"))
```

## Usage

```r
source("fig1d_mrna_heatmap.R")
```

## Input Data

**File**: `Fig_1d_Fold_change_mRNA_expression.xlsx`

**Columns**:
- `Treatment`: CTRL, DAPT, DKK-1
- `Gene`: Gene symbols (Hes1, Hey1, Wnt2, Wnt2b, Wnt5a, Wnt5b, Wnt9a, Lef1, Tcf7, Tcf7l2)
- `Fold_change`: Fold change values
- `Marker`: Marker category (filters out "Pancreatic functional markers")

## Output

- `Fig_1d_Heatmap_Fold_change_mRNA_expression.png` (600 dpi, 16.7 × 7.4 cm)
- `Fig_1d_Heatmap_Fold_change_mRNA_expression.svg`

## Visualization

- Three-color gradient: blue (#3c5488) to grey to red (#e64b35)
- Range: 0 to 2.2 fold change
- Numerical values displayed on tiles
- Gene names in italic
