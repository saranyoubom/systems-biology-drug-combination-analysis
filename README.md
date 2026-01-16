# COMPLETE ANALYSIS PIPELINE - FIGURES 1-6

## Notch/Wnt Signaling in β-Cell Function & Type 2 Diabetes Drug Repositioning via Network Pharmacology

**Project:** Integrative Systems Biology Analysis of Notch and Wnt Signaling Pathways in INS1 β-Cells with Drug Repositioning Strategy  
**Version:** 1.0  
**Date:** 16 January 2026  
**Platform:** R 4.0+  
**License:** MIT  

---

## TABLE OF CONTENTS

1. [Project Overview](#1-project-overview)
2. [Experimental Design](#2-experimental-design)
3. [System Requirements](#3-system-requirements)
4. [Installation Guide](#4-installation-guide)
5. [Complete Analysis Workflow (Figures 1-6)](#5-complete-analysis-workflow-figures-1-6)
6. [File Organization](#6-file-organization)
7. [Data Dependencies](#7-data-dependencies)
8. [Usage Instructions](#8-usage-instructions)
9. [Troubleshooting](#9-troubleshooting)
10. [Citation](#10-citation)
11. [Contact Information](#11-contact-information)

---

## 1. PROJECT OVERVIEW

This repository contains the complete computational and statistical analysis pipeline for investigating Notch and Wnt signaling crosstalk in pancreatic β-cell function and identifying drug combination candidates for Type 2 Diabetes (T2D) therapy.

### STUDY OBJECTIVES:

- Characterize Notch/Wnt signaling modulation effects on β-cell proliferation
- Profile glucose uptake and insulin secretion responses
- Map gene co-expression networks using multi-method validation
- Identify hub genes critical for β-cell function
- Develop Multi-Layered Synergy Score (MLSS) algorithm for drug combinations
- Integrate regulatory approval data (FDA/EMA/PMDA)
- Predict optimal drug repositioning candidates

### BIOLOGICAL CONTEXT:

Type 2 Diabetes is characterized by progressive β-cell dysfunction and loss. This study investigates the role of Notch (DAPT inhibition) and Wnt (DKK-1 antagonism) signaling in β-cell proliferation, insulin secretion, and glucose homeostasis, with implications for therapeutic intervention.

---

## 2. EXPERIMENTAL DESIGN

### CELL MODEL:

- INS1 rat insulinoma cells (β-cell line)
- Treatments: CTRL, DAPT (Notch inhibitor), DKK-1 (Wnt antagonist)
- Biological replicates: n=3-6 per group

### ASSAYS PERFORMED:

1. **Cell proliferation** - DNA quantification
2. **Glucose uptake** - 2-NBDG fluorescence
3. **Insulin secretion** - C-peptide ELISA
4. **Gene expression** - qRT-PCR
   - Notch targets: Hes1, Hey1
   - Wnt pathway: Wnt2, Wnt2b, Wnt5a, Wnt5b, Wnt9a, Lef1, Tcf7, Tcf7l2
   - Insulin synthesis: Ins1, Glp1r
   - Calcium signaling: Cacna1c, Cacna1d, Ryr2, Itpr1

### COMPUTATIONAL ANALYSES:

1. Statistical testing (ANOVA, t-test, Kruskal-Wallis)
2. Gene co-expression network construction (4 methods)
3. Network validation (consensus approach)
4. Hub gene identification (centrality metrics)
5. Drug-target interaction mapping
6. MLSS v4.0 scoring algorithm
7. Regulatory approval integration

---

## 3. SYSTEM REQUIREMENTS

### HARDWARE:

- **Processor:** Intel Core i5 or equivalent (2.5 GHz+)
- **RAM:** Minimum 8 GB (16 GB recommended for Fig 4-6)
- **Storage:** 10 GB free space
- **Graphics:** Integrated graphics sufficient

### SOFTWARE:

- **Operating System:** Windows 10+, macOS 10.14+, or Linux (Ubuntu 18.04+)
- **R version:** 4.0.0 or higher (4.2+ recommended)
- **RStudio:** 1.4+ (optional but recommended)
- **Excel/LibreOffice:** For viewing .xlsx data files

### R PACKAGE DEPENDENCIES:

#### Core Data Manipulation:
- tidyverse (1.3.0+)
- dplyr, tidyr, tibble
- data.table
- readxl, readr

#### Statistical Analysis:
- ggpubr (statistical comparisons)
- Hmisc (correlation analysis)
- psych (partial correlation)

#### Visualization:
- ggplot2 (3.3.0+)
- ggsci (NPG color palettes)
- scales
- pheatmap (heatmaps)
- factoextra (PCA visualization)
- ggfortify
- ggrepel
- patchwork (multi-panel figures)

#### Network Analysis:
- igraph (1.2.6+)
- ggraph
- WGCNA (Bioconductor)

#### Specialized:
- ggalluvial (Sankey diagrams)
- VennDiagram
- gridExtra

---

## 4. INSTALLATION GUIDE

### STEP 1: INSTALL R AND RSTUDIO

```bash
# Download R from CRAN
https://cran.r-project.org/

# Download RStudio (optional)
https://www.rstudio.com/products/rstudio/download/
```

### STEP 2: INSTALL REQUIRED R PACKAGES

```r
# Install CRAN packages
install.packages(c(
  "tidyverse", "dplyr", "tidyr", "tibble", "data.table",
  "readxl", "readr", "ggplot2", "ggsci", "scales", "pheatmap",
  "ggpubr", "factoextra", "ggfortify", "ggrepel", "patchwork",
  "igraph", "ggraph", "ggalluvial", "VennDiagram", "gridExtra",
  "Hmisc", "psych"
))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("WGCNA")
```

### STEP 3: VERIFY INSTALLATION

```r
# Test core packages
library(tidyverse)
library(ggplot2)
library(igraph)
library(ggsci)

# Check R version
R.version.string
```

### STEP 4: CLONE/DOWNLOAD REPOSITORY

```bash
# Option A: Git clone
git clone https://github.com/[username]/[repository-name].git

# Option B: Download ZIP
# Extract to working directory
```

### STEP 5: SET WORKING DIRECTORY

```r
# In R/RStudio:
setwd("path/to/repository")

# Verify:
getwd()
list.files()
```

---

## 5. COMPLETE ANALYSIS WORKFLOW (FIGURES 1-6)

### FIGURE 1: IN VITRO VALIDATION

**Notch/Wnt Modulation Effects on β-Cells**

#### Figure 1c: INS1 Cell Proliferation

- **Script:** `fig1c_ins1_proliferation.r`
- **Purpose:** DNA content quantification (proliferation assay)
- **Input:** `Fig_1c_Boxplot_INS1_proliferation.xlsx`
- **Assay:** PicoGreen DNA quantification
- **Statistics:** ANOVA + pairwise t-tests
- **Output:** Boxplot with jitter (8.08 × 4.57 cm, 600 DPI)

#### Figure 1d: Gene Expression Heatmap

- **Script:** `fig1d_mrna_heatmap.r`
- **Purpose:** mRNA fold change visualization
- **Input:** `Fig_1d_Fold_change_mRNA_expression.xlsx`
- **Genes:** Hes1, Hey1, Wnt family, Tcf/Lef transcription factors
- **Visualization:** Tile heatmap with gradient (blue-white-red)
- **Output:** Heatmap (16.7 × 7.4 cm, 600 DPI)

---

### FIGURE 2: FUNCTIONAL VALIDATION

**Glucose Uptake and Gene Expression Profiling**

#### Figure 2b: Glucose Uptake Analysis

- **Script:** `fig2b_glucose_uptake.r`
- **Purpose:** 2-NBDG fluorescence quantification
- **Input:** `Fig_2b_Boxplot_Glucose_uptake.xlsx`
- **Statistics:** Kruskal-Wallis + Wilcoxon tests (non-parametric)
- **Output:** Boxplot (8 × 5.3 cm, 600 DPI)

#### Figure 2c-f: Functional Gene Expression

- **Script:** `fig2cdef_functional_genes.r`
- **Purpose:** Multi-gene expression profiling
- **Input:** `Fig_2cdef_Boxplot_Fold_change_mRNA_expression.xlsx`
- **Genes:** Slc2a2 (GLUT2), Pdx1, Neurod1, Nkx6.1
- **Statistics:** Mixed ANOVA/Kruskal-Wallis (normality-dependent)
- **Output:** 4 separate boxplots (8 × 5.3 cm each, 600 DPI)

---

### FIGURE 3: INSULIN SECRETION ASSAY

**C-peptide Release and Gene Expression Profiling**

#### Figure 3a: C-peptide Secretion Curve

- **Script:** `fig3a_cpeptide_curve.r`
- **Purpose:** Glucose-stimulated insulin secretion (GSIS)
- **Input:** `Fig_3a_Scatterplot_C_peptide.xlsx`
- **Glucose:** 2.8, 5.5, 11, 22 mM concentrations
- **Visualization:** LOESS smoothing curves with confidence bands
- **Output:** Scatter + LOESS (8 × 6.5 cm, 600 DPI)

#### Figure 3b: C-peptide at 5.5 mM Glucose

- **Script:** `fig3b_cpeptide_5.5mm.r`
- **Purpose:** Basal insulin secretion
- **Input:** `Fig_3b_Boxplot_C_peptide.xlsx`
- **Statistics:** ANOVA + t-tests
- **Output:** Boxplot (8 × 5.3 cm, 600 DPI)

#### Figure 3c-d: Insulin Pathway Genes

- **Script:** `fig3cd_insulin_genes.r`
- **Purpose:** Insulin synthesis and secretion markers
- **Input:** `Fig_3cd_Boxplot_Fold_change_mRNA_expression.xlsx`
- **Genes:** Glp1r (incretin receptor), Ins1 (insulin)
- **Output:** 2 boxplots (8 × 5.3 cm each, 600 DPI)

#### Figure 3e: Calcium Signaling Genes

- **Script:** `fig3e_calcium_genes.r`
- **Purpose:** Ca²⁺ homeostasis markers
- **Genes:** Cacna1c, Cacna1d (L-type Ca²⁺ channels)
- **Output:** Bar plot with log fold change (8 × 5.3 cm, 600 DPI)

#### Figure 3f: Principal Component Analysis

- **Script:** `fig3f_pca_analysis.r`
- **Purpose:** Dimensionality reduction and treatment clustering
- **Input:** `Fig_3f_PCA_Relative_mRNA_expression.xlsx`
- **Method:** PCA biplot with 95% confidence ellipses
- **Output:** PCA biplot (7.6 × 5.3 cm, 600 DPI)

---

### FIGURE 4: MULTI-METHOD NETWORK CONSTRUCTION

**Gene Co-Expression Networks with Consensus Validation**

#### Figure 4b: Correlation Matrix Heatmap

- **Script:** `fig4b_correlation_matrix.r`
- **Purpose:** Pearson correlation visualization
- **Genes:** 21 genes (Wnt, Notch, Ca²⁺, insulin pathways)
- **Visualization:** Heatmap with hierarchical clustering
- **Output:** Correlation heatmap (600 DPI)

#### Figure 4c: Normalization Method Comparison

- **Script:** `fig4c_normalization_methods.r`
- **Purpose:** Compare correlation methods

**Methods:**
1. Pearson (parametric)
2. Spearman (rank-based)
3. Partial correlation (adjusted)
4. WGCNA (weighted network)

- **Output:** Method comparison plots (600 DPI)

#### Figure 4d: Venn Diagram - Consensus Pairs

- **Script:** `fig4d_venn_diagram.r`
- **Purpose:** 4-method consensus identification
- **Criteria:** Gene pairs validated by all 4 methods
- **Output:** `Fig_4d_consensus_pairs.csv` (consensus list)
- **Visualization:** 4-way Venn diagram (600 DPI)

#### Figure 4e: Hes1-Normalized Correlation Matrix

- **Script:** (part of fig4d workflow)
- **Purpose:** Residual correlation after Hes1 adjustment
- **Output:** `Fig_4e_Hes1_significant_correlations_only.csv`
- **Use:** Input for Figure 5 and Figure 6b-f

---

### FIGURE 5: NETWORK VALIDATION & HUB GENES

**STRING Database Integration and Novel Discovery Analysis**

#### Figure 5: Gene Interaction Network Analysis

- **Script:** `fig5_network_analysis.r`
- **Purpose:** Validate consensus pairs with STRING database

##### METHODOLOGY:

1. Load consensus pairs from Figure 4d
2. Load STRING interaction database
3. Classify pairs:
   - **Known Validated:** In both consensus + STRING
   - **Novel Discovery:** Consensus only (not in STRING)
   - **STRING Only:** Database only (not validated)
4. Calculate centrality metrics:
   - Degree centrality (connectivity)
   - Betweenness centrality (information flow)
   - Closeness centrality (network proximity)
   - Eigenvector centrality (influence)

##### HUB GENE CRITERIA:

- Degree ≥ 6 connections
- High betweenness (>median)
- Present in consensus pairs

##### OUTPUT:

- `Fig5a_Network.png/svg` - Force-directed network graph
- Node color: Novel (red) vs Known (blue)
- Edge color: Validation status
- Edge width: Correlation strength
- CSV: Novel discoveries with statistics

---

### FIGURE 6: DRUG REPOSITIONING ANALYSIS

**MLSS v4.0 Algorithm & Regulatory Approval Integration**

#### Figure 6b: Consensus Hub Gene Network

- **Script:** `fig6b_consensus_hub_network.r`
- **Purpose:** Visualize gene-gene correlation network

##### INPUTS:

- `Fig_4d_consensus_pairs.csv`
- `Fig_4e_Hes1_significant_correlations_only.csv`

##### FEATURES:

- Circular layout (Fruchterman-Reingold algorithm)
- Node size = degree centrality
- Edge width = |correlation strength|
- Edge linetype = sign (solid=positive, dotted=negative)
- NPG colors by functional category:
  - Dark blue: Calcium Signaling
  - Red/Orange: Incretin Signaling
  - Green: Wnt Signaling
  - Light blue: Metabolic Regulation

##### OUTPUT:

- `fig6b_consensus_hub_network.png/svg` (600 DPI)
- `fig6b_gene_centrality_summary.csv` (network metrics)
- `fig6b_edge_correlation_statistics.csv`

---

#### Figure 6c: MLSS v4.0 Drug Combination Scoring

- **Script:** `fig6c_mlss_drug_combinations.r`
- **Purpose:** Rank drug combinations using MLSS algorithm

##### MLSS v4.0 FORMULA:

```
MLSS = [C×0.45 + B×0.30 + V×0.10 + P×0.15] × 9
     + Σ(r_pos × Potency_Multiplier) × 1
     - Σ(|r_neg| × Potency_Multiplier) × 1
```

##### COMPONENTS:

- **C = Complementarity** (target non-overlap)
  - Range: 0 (complete overlap) to 1 (disjoint)
- **B = Balance** (target count similarity)
  - Range: 0 (highly imbalanced) to 1 (equal)
- **V = Network Coverage** (correlation density, normalized)
- **P = Potency** (pIC50-derived binding strength)
  - P = min(max(pIC50 - 4, 0), 5) / 5

##### POTENCY MULTIPLIER:

```
Potency_Multiplier = min(max(pIC50_i, pIC50_j) / 9, 1.0)
```

##### INPUTS:

- `fig6b_gene_centrality_summary.csv` (hub genes)
- `Fig_4e_Hes1_significant_correlations_only.csv`
- `drug.target.interaction.csv` (pIC50 data)

##### ALGORITHM FLOW:

1. Filter hub genes (degree ≥ 6)
2. Extract drug-target interactions for hub genes
3. Generate all pairwise drug combinations
4. Calculate base components (C, B, V, P)
5. Compute synergy bonus (positive correlations × potency)
6. Compute antagonism penalty (negative correlations × potency)
7. Aggregate into final MLSS score
8. Rank combinations

##### OUTPUT:

- `fig6c_mlss_top15.png/svg` (stacked bar chart, 600 DPI)
- `fig6c_mlss_all_combinations.csv` (full ranking)
- `fig6c_mlss_top15.csv` (top 15 with detailed metrics)

---

#### Figure 6d + Supplementary S2: MLSS Validation

- **Script:** `fig6d_figs2_mlss_validation.r`
- **Purpose:** Sensitivity and robustness analysis

##### ANALYSES PERFORMED:

**1. Weight Sensitivity (8 scenarios)**
- MLSS Standard: 45:30:10:15 (C:B:V:P)
- C-Heavy variants: 60:25:15, 55:20:25
- B-Heavy: 35:50:15
- V-Heavy: 35:30:35
- P-Heavy: 30:20:10:40
- Balanced: 40:40:10:10
- Network-Focus: 30:20:40:10

**2. Ablation Study (leave-one-out)**
- Without C, B, V, P (renormalized)
- Base only (no synergy/antagonism)

**3. Ranking Correlation**
- Spearman ρ across all scenarios
- Expected: ρ > 0.85 (robust)

**4. Stability Metrics**
- Top 5 overlap consistency
- Mean rank shift (positions)
- Mean score drop (MLSS units)

##### OUTPUT:

- `figs2a_mlss_sensitivity_top15.png/svg` (heatmap)
- `figs2b_component_importance_top15.png/svg` (bar chart)
- `fig6d_ranking_correlations.csv`
- `fig6d_ablation_impact.csv`

---

#### Figure 6e: Drug Combination Network

- **Script:** `fig6e_drug_combination_network.r`
- **Purpose:** Visualize drug-drug relationships from top MLSS
- **Input:** `fig6c_mlss_all_combinations.csv`

##### NETWORK PROPERTIES:

- **Nodes:** Individual drugs
  - Size = appearance frequency in top 15
  - Color = primary functional pathway
  - Label = bold for hubs (≥4 appearances)
- **Edges:** Drug pairs from top 15 combinations
  - Width = MLSS score magnitude
  - Linetype = cross-pathway (solid) vs same-pathway (dotted)
  - Color = grey70 (uniform)
- **Layout:** Circular (balanced visualization)

##### HUB DRUGS (≥4 appearances):

- **Calcium pathway:** Isradipine, Nitrendipine, (S)-nitrendipine
- **Incretin pathway:** Albiglutide, Liraglutide

##### OUTPUT:

- `fig6e_drug_combination_network.png/svg` (600 DPI)
- `fig6e_drug_node_statistics.csv`
- `fig6e_drug_edge_statistics.csv`

---

#### Figure 6f: Regulatory Approval & Mechanistic Sankey

- **Script:** `fig6f_regulatory_mechanistic_sankey.r`
- **Purpose:** 6-layer alluvial diagram integrating FDA/EMA/PMDA data

##### INPUTS:

- `fig6c_mlss_top15.csv`
- `fig6b_gene_centrality_summary.csv`
- `FDA_Approved.csv`
- `EMA_Approved.csv`
- `PMDA_Approved.csv`
- `FDA-EMA-PMDA_Approved.csv` (triple-approved)

##### 6-LAYER STRUCTURE:

**Layers 1-4:**
- Layer 1: Approval Status (FDA-EMA-PMDA, Not Approved)
- Layer 2: Hub Drugs (Albiglutide, Liraglutide, Isradipine, Nitrendipine)
- Layer 3: Action Type (Agonist, Antagonist)
- Layer 4: Hub Genes (Glp1r, Tcf7, Ptbp1, Cacna1d, Cacna1c)

**Layers 5-6:**
- Layer 5: Pathways (Incretin Signaling, Calcium Signaling)
- Layer 6: Clinical Outcomes (Enhanced insulin secretion, β-cell preservation, Metabolic homeostasis)

##### FLOW CALCULATION:

```
Flow = Frequency × Connection_Strength × Confidence
```

##### CONFIDENCE WEIGHTS:

- Tier 1 (FDA+EMA+PMDA): 3.0
- Tier 2 (FDA+EMA): 2.0
- Tier 3 (FDA only): 1.0
- Not approved: 0.3

##### OUTPUT:

- `fig6f_regulatory_mechanistic_sankey.png/svg` (600 DPI)
- `fig6f_hub_drug_approval_profile.csv`
- `fig6f_six_layer_complete_mapping.csv`
- `fig6f_pathway_gene_connectivity.csv`

---

### SUPPLEMENTARY FIGURE S1: VALIDATION

- **Script:** `figS1_validation_analysis.r`
- **Purpose:** Cross-validation and reproducibility testing
- **Output:** Validation plots and statistics (600 DPI)

---

## 6. FILE ORGANIZATION

### RECOMMENDED DIRECTORY STRUCTURE:

```
project_root/
│
├── README.md                                # This comprehensive guide
│
├── scripts/                                 # All R analysis scripts
│   ├── fig1c_ins1_proliferation.r
│   ├── fig1d_mrna_heatmap.r
│   ├── fig2b_glucose_uptake.r
│   ├── fig2cdef_functional_genes.r
│   ├── fig3a_cpeptide_curve.r
│   ├── fig3b_cpeptide_5.5mm.r
│   ├── fig3cd_insulin_genes.r
│   ├── fig3e_calcium_genes.r
│   ├── fig3f_pca_analysis.r
│   ├── fig4b_correlation_matrix.r
│   ├── fig4c_normalization_methods.r
│   ├── fig4d_venn_diagram.r
│   ├── fig5_network_analysis.r
│   ├── fig6b_consensus_hub_network.r
│   ├── fig6c_mlss_drug_combinations.r
│   ├── fig6d_figs2_mlss_validation.r
│   ├── fig6e_drug_combination_network.r
│   ├── fig6f_regulatory_mechanistic_sankey.r
│   └── figS1_validation_analysis.r
│
├── data/                                    # Input data files
│   ├── experimental/                        # Raw experimental data (.xlsx)
│   │   ├── Fig_1c_Boxplot_INS1_proliferation.xlsx
│   │   ├── Fig_1d_Fold_change_mRNA_expression.xlsx
│   │   ├── Fig_2b_Boxplot_Glucose_uptake.xlsx
│   │   ├── Fig_2cdef_Boxplot_Fold_change_mRNA_expression.xlsx
│   │   ├── Fig_3a_Scatterplot_C_peptide.xlsx
│   │   ├── Fig_3b_Boxplot_C_peptide.xlsx
│   │   ├── Fig_3cd_Boxplot_Fold_change_mRNA_expression.xlsx
│   │   └── Fig_3f_PCA_Relative_mRNA_expression.xlsx
│   │
│   ├── processed/                           # Intermediate analysis outputs
│   │   ├── Fig_4d_consensus_pairs.csv
│   │   └── Fig_4e_Hes1_significant_correlations_only.csv
│   │
│   └── databases/                           # External databases
│       ├── drug.target.interaction.csv
│       ├── string_interactions.csv
│       ├── FDA_Approved.csv
│       ├── EMA_Approved.csv
│       ├── PMDA_Approved.csv
│       └── FDA-EMA-PMDA_Approved.csv
│
├── results/                                 # Analysis outputs
│   ├── figures/                             # Publication figures
│   │   ├── fig1c_*.png, fig1c_*.svg
│   │   ├── fig1d_*.png, fig1d_*.svg
│   │   ├── ...
│   │   └── fig6f_*.png, fig6f_*.svg
│   │
│   └── tables/                              # Statistical tables
│       ├── fig6b_gene_centrality_summary.csv
│       ├── fig6c_mlss_all_combinations.csv
│       ├── fig6c_mlss_top15.csv
│       ├── fig6d_ranking_correlations.csv
│       ├── fig6e_drug_node_statistics.csv
│       └── fig6f_hub_drug_approval_profile.csv
│
└── docs/                                    # Additional documentation
    └── methods_details.md
```

---

## 7. DATA DEPENDENCIES

### INPUT DATA FILES REQUIRED:

**Figure 1:**
- `Fig_1c_Boxplot_INS1_proliferation.xlsx`
- `Fig_1d_Fold_change_mRNA_expression.xlsx`

**Figure 2:**
- `Fig_2b_Boxplot_Glucose_uptake.xlsx`
- `Fig_2cdef_Boxplot_Fold_change_mRNA_expression.xlsx`

**Figure 3:**
- `Fig_3a_Scatterplot_C_peptide.xlsx`
- `Fig_3b_Boxplot_C_peptide.xlsx`
- `Fig_3cd_Boxplot_Fold_change_mRNA_expression.xlsx`
- `Fig_3f_PCA_Relative_mRNA_expression.xlsx`

**Figure 4:**
- (Generated from raw expression data)
- **Outputs:** 
  - `Fig_4d_consensus_pairs.csv`
  - `Fig_4e_Hes1_significant_correlations_only.csv`

**Figure 5:**
- `Fig_4d_consensus_pairs.csv` (from Figure 4)
- `Fig_4e_Hes1_significant_correlations_only.csv` (from Figure 4)
- `string_interactions.csv` (STRING database)

**Figure 6:**
- `fig6b_gene_centrality_summary.csv` (from Figure 6b)
- `Fig_4e_Hes1_significant_correlations_only.csv` (from Figure 4)
- `drug.target.interaction.csv` (DrugCentral)
- `FDA_Approved.csv`
- `EMA_Approved.csv`
- `PMDA_Approved.csv`
- `FDA-EMA-PMDA_Approved.csv`

### EXTERNAL DATABASE SOURCES:

**1. DrugCentral (2023)**
- **URL:** https://drugcentral.org/
- **File:** `drug.target.interaction.csv`
- **Contains:** Drug names, targets, pIC50/Ki/EC50 values, organisms

**2. STRING v11.5**
- **URL:** https://string-db.org/
- **File:** `string_interactions.csv`
- **Contains:** Protein-protein interactions, confidence scores

**3. FDA Orange Book**
- **URL:** https://www.fda.gov/drugs/drug-approvals-and-databases/
- **File:** `FDA_Approved.csv`

**4. EMA Product Database**
- **URL:** https://www.ema.europa.eu/
- **File:** `EMA_Approved.csv`

**5. PMDA Approved Drugs**
- **URL:** https://www.pmda.go.jp/english/
- **File:** `PMDA_Approved.csv`

---

## 8. USAGE INSTRUCTIONS

### OPTION 1: RUN COMPLETE PIPELINE

```r
# Set working directory
setwd("path/to/project_root")

# Run all scripts in order
source("scripts/fig1c_ins1_proliferation.r")
source("scripts/fig1d_mrna_heatmap.r")
source("scripts/fig2b_glucose_uptake.r")
source("scripts/fig2cdef_functional_genes.r")
source("scripts/fig3a_cpeptide_curve.r")
source("scripts/fig3b_cpeptide_5.5mm.r")
source("scripts/fig3cd_insulin_genes.r")
source("scripts/fig3e_calcium_genes.r")
source("scripts/fig3f_pca_analysis.r")
source("scripts/fig4b_correlation_matrix.r")
source("scripts/fig4c_normalization_methods.r")
source("scripts/fig4d_venn_diagram.r")
source("scripts/fig5_network_analysis.r")
source("scripts/fig6b_consensus_hub_network.r")
source("scripts/fig6c_mlss_drug_combinations.r")
source("scripts/fig6d_figs2_mlss_validation.r")
source("scripts/fig6e_drug_combination_network.r")
source("scripts/fig6f_regulatory_mechanistic_sankey.r")
source("scripts/figS1_validation_analysis.r")
```

### OPTION 2: RUN INDIVIDUAL FIGURES

```r
# Example: Run only Figure 3 (insulin secretion analysis)
setwd("path/to/project_root")
source("scripts/fig3a_cpeptide_curve.r")
source("scripts/fig3b_cpeptide_5.5mm.r")
source("scripts/fig3cd_insulin_genes.r")
source("scripts/fig3e_calcium_genes.r")
source("scripts/fig3f_pca_analysis.r")
```

### OPTION 3: RUN ONLY NETWORK PHARMACOLOGY (FIG 4-6)

```r
# Requires Fig 4 outputs as dependencies
setwd("path/to/project_root")
source("scripts/fig4d_venn_diagram.r")         # Generates consensus
source("scripts/fig5_network_analysis.r")
source("scripts/fig6b_consensus_hub_network.r")
source("scripts/fig6c_mlss_drug_combinations.r")
source("scripts/fig6d_figs2_mlss_validation.r")
source("scripts/fig6e_drug_combination_network.r")
source("scripts/fig6f_regulatory_mechanistic_sankey.r")
```

### CUSTOM PARAMETER ADJUSTMENT:

```r
# Modify MLSS weights (in fig6c script)
C_weight <- 0.45  # Complementarity (default)
B_weight <- 0.30  # Balance
V_weight <- 0.10  # Network Coverage
P_weight <- 0.15  # Potency

# Change top N combinations
top_n <- 15  # Increase to 30, 50, etc.

# Adjust hub gene threshold
hub_threshold <- 6  # Default: degree >= 6

# Modify significance thresholds
alpha <- 0.05          # P-value cutoff
correlation_cutoff <- 0.7  # |r| threshold
```

---

## 9. TROUBLESHOOTING

### Issue 1: Package installation fails

**Solution:**
```r
# Update R to latest version
# Try installing from source
install.packages("package_name", type = "source")

# For WGCNA on Windows:
BiocManager::install("WGCNA", type = "binary")
```

### Issue 2: File not found errors

**Solution:**
```r
# Check working directory
getwd()

# Verify file paths are relative
list.files("data/experimental")

# Ensure all data files are present
# Re-run upstream scripts (Fig 4 before Fig 5-6)
```

### Issue 3: Memory errors (large datasets)

**Solution:**
```r
# Increase memory limit (Windows)
memory.limit(size = 16000)

# Close unnecessary programs
# Use data.table for large files
# Process data in chunks
```

### Issue 4: Graphics display issues

**Solution:**
```r
# Reset graphics device
dev.off()

# Specify device explicitly
png("output.png", width = 800, height = 600)
print(plot_object)
dev.off()
```

### Issue 5: MLSS scores seem incorrect

**Solution:**
```r
# Verify input files:
# - fig6b_gene_centrality_summary.csv exists
# - drug.target.interaction.csv has pIC50 data
# - Check hub gene threshold (degree >= 6)

# Inspect intermediate values
print(head(all_drug_combos))
summary(all_drug_combos$MLSS_v4.0)
```

### Issue 6: Network visualization crashes

**Solution:**
```r
# Reduce network size
top_n <- 10  # Instead of 50

# Simplify layout
layout <- layout_with_fr(g)  # Instead of complex algorithms

# Save plot to file instead of displaying
ggsave("network.png", plot = p, dpi = 300)
```

### Issue 7: Statistical tests give warnings

**Solution:**
- Check sample sizes (n ≥ 3 per group)
- Verify normality assumptions (Shapiro-Wilk test)
- Use non-parametric tests if needed (Kruskal-Wallis, Wilcoxon)

### Issue 8: NPG colors not displaying correctly

**Solution:**
```r
# Reinstall ggsci package
install.packages("ggsci")
library(ggsci)

# Verify palette
pal_npg("nrc")(10)
```

### PERFORMANCE OPTIMIZATION:

**For faster execution:**
- Reduce DPI (600 → 300) during testing
- Limit top_n parameter (50 → 15)
- Skip validation analyses (Fig 6d)
- Use parallel processing:

```r
library(parallel)
cl <- makeCluster(detectCores() - 1)
```

**For large datasets (>10,000 genes):**
- Filter low-variance genes before correlation
- Use data.table instead of data.frame
- Implement chunked processing

---

## LICENSE INFORMATION

**MIT License**

Copyright (c) 2026 [Saranyou Oontawee, Chulalongkorn University]

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

---

## END OF README

**Document Version:** 1.0  
**Last Updated:** January 16, 2026  
**Total Analysis Scripts:** 20  
**Total Output Figures:** 25+  
**Estimated Analysis Time:** 2-4 hours (complete pipeline)  

For questions or support, please contact: [saranyou.o@chula.ac.th]

---
