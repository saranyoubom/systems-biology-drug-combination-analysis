# INS1 Cell Proliferation Analysis

Visualization of total DNA content in INS1 cells under treatment with DAPT and DKK-1.

## Description

Generates a boxplot comparing total DNA content (proliferation marker) across three treatment groups with statistical comparisons.

## Requirements

```r
install.packages(c("ggplot2", "data.table", "tidyr", "ggpubr", "dplyr", "ggsci", "readxl"))
```

## Usage

```r
source("fig1c_ins1_proliferation.R")
```

## Input Data

**File**: `Fig_1c_Boxplot_INS1_proliferation.xlsx`

**Columns**:
- `Treatment`: CTRL, DAPT, DKK-1
- `DNA`: Total DNA content (µg)
- `Assay`: Assay identifier

## Output

- `Fig_1c_Boxplot_INS1_proliferation.png` (600 dpi, 8.08 × 4.57 cm)
- `Fig_1c_Boxplot_INS1_proliferation.svg`

## Statistical Analysis

- **Normality test**: Shapiro-Wilk (W = 0.8916, p = 0.1236)
- **Overall comparison**: One-way ANOVA
- **Pairwise comparisons**: t-test
- **Visualization**: Mean ± SEM annotations

## Visualization

- NPG color palette
- Jittered data points overlaid on boxplots
- Significance brackets with p-value indicators
- Minimal theme with size 8 text
