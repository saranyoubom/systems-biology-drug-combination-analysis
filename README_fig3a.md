# C-peptide Secretion Profile

LOESS-smoothed dose-response curves of C-peptide secretion across glucose concentrations.

## Description

Visualizes glucose-stimulated C-peptide secretion using local regression smoothing.

## Requirements

```r
install.packages(c("ggplot2", "readr", "ggsci", "readxl"))
```

## Usage

```r
source("fig3a_cpeptide_curve.R")
```

## Input Data

**File**: `Fig_3a_Scatterplot_C_peptide.xlsx`

**Columns**:
- `Treatment`: CTRL, DAPT, DKK-1
- `Glucose`: Glucose concentration (mM)
- `C_peptide`: Relative C-peptide release (AU)
- `Assay`: Assay identifier

## Output

- `Fig_3a_Scatterplot_C_peptide.png` (600 dpi, 8 × 6.5 cm)
- `Fig_3a_Scatterplot_C_peptide.svg`

## Analysis Method

- **Smoothing**: LOESS (locally estimated scatterplot smoothing)
- **Standard error bands**: 95% confidence intervals
- **Note**: 44 mM glucose concentration excluded from analysis

## Visualization

- Individual data points overlaid with LOESS curves
- NPG color palette
- Standard error shading
