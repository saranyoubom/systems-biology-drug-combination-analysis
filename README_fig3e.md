# Calcium Signaling Gene Expression

Log2 fold change barplot for intracellular calcium handling genes.

## Description

Displays log-transformed expression changes for genes involved in calcium signaling (Ptbp1, Itpr1, Ryr1, Ryr2, Ryr3).

## Requirements

```r
install.packages(c("ggplot2", "dplyr", "ggsci", "readxl"))
```

## Usage

```r
source("fig3e_calcium_genes.R")
```

## Input Data

**File**: `Fig_3e_Barplot_Log_Fold_Change_mRNA_expression.xlsx`

**Columns**:
- `Treatment`: DAPT, DKK-1 (CTRL excluded)
- `Gene`: Ptbp1, Itpr1, Ryr1, Ryr2, Ryr3
- `LogFC`: Log2 fold change
- `Marker`: "Intracellular markers" (used for filtering)

## Output

- `Fig_3e_Barplot_Log_Fold_Change_mRNA_expression.png` (600 dpi, 8 × 5.3 cm)
- `Fig_3e_Barplot_Log_Fold_Change_mRNA_expression.svg`

## Visualization

- Mean ± SEM error bars
- Custom NPG color palette (2nd and 3rd colors)
- Dodged bars for treatment comparison
- Gene names in italic

## Gene Functions

- **Ptbp1**: Polypyrimidine tract binding protein (RNA splicing)
- **Itpr1**: Inositol 1,4,5-trisphosphate receptor type 1
- **Ryr1/2/3**: Ryanodine receptors (calcium release channels)
