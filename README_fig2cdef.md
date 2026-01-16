# Functional Gene Expression

Expression analysis of key pancreatic beta cell functional genes (Glut2, Kcnj11, Cacna1c, Cacna1d).

## Description

Creates four separate boxplots for genes involved in glucose sensing and ion channel function.

## Requirements

```r
install.packages(c("ggplot2", "data.table", "tidyr", "ggpubr", "dplyr", "ggsci", "readxl"))
```

## Usage

```r
source("fig2cdef_functional_genes.R")
```

## Input Data

**File**: `Fig_2cdef_Fold_change_mRNA_expression.xlsx`

**Columns**:
- `Treatment`: CTRL, DAPT, DKK-1
- `Gene`: Glut2, Kcnj11, Cacna1c, Cacna1d
- `Fold_change`: Fold change in mRNA expression

## Output

Four separate figure sets (PNG + SVG):
- Fig_2c: Glut2 expression
- Fig_2d: Kcnj11 expression
- Fig_2e: Cacna1c expression
- Fig_2f: Cacna1d expression

## Statistical Analysis

**Per-gene normality tests**:
- Glut2: W = 0.9017, p = 0.1667 (normal) → ANOVA + t-test
- Kcnj11: W = 0.8584, p = 0.0467 (normal) → ANOVA + t-test
- Cacna1c: W = 0.8395, p = 0.0273 (non-normal) → Kruskal-Wallis + Wilcoxon
- Cacna1d: W = 0.8832, p = 0.0964 (normal) → ANOVA + t-test

## Gene Functions

- **Glut2**: Glucose transporter and sensor
- **Kcnj11**: ATP-sensitive potassium channel subunit
- **Cacna1c/Cacna1d**: Voltage-gated calcium channels (insulin secretion)
