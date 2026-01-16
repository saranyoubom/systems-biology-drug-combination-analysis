# Insulin-related Gene Expression

Expression of Glp1r (GLP-1 receptor) and Ins1 (insulin 1) genes.

## Description

Generates two boxplots for key genes in insulin synthesis and secretion pathways.

## Requirements

```r
install.packages(c("ggplot2", "data.table", "tidyr", "ggpubr", "dplyr", "ggsci", "readxl"))
```

## Usage

```r
source("fig3cd_insulin_genes.R")
```

## Input Data

**File**: `Fig_3cd_Boxplot_Fold_change_mRNA_expression.xlsx`

**Columns**:
- `Treatment`: CTRL, DAPT, DKK-1
- `Gene`: Glp1r, Ins1
- `Fold_change`: Fold change in mRNA expression

## Output

Two separate figure sets (PNG + SVG):
- Fig_3c: Glp1r expression
- Fig_3d: Ins1 expression

## Statistical Analysis

**Per-gene normality tests**:
- Glp1r: W = 0.8817, p = 0.0922 (normal) → ANOVA + t-test
- Ins1: W = 0.9034, p = 0.1754 (normal) → ANOVA + t-test

## Gene Functions

- **Glp1r**: GLP-1 receptor, regulates insulin synthesis and glucose-dependent secretion
- **Ins1**: Insulin 1 gene, encodes proinsulin
