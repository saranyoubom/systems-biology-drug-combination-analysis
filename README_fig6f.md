# Regulatory Approval & Mechanistic Pathway Sankey - Figure 6f

Six-layer alluvial diagram integrating regulatory status with mechanistic pathways.

## Description

Flow diagram showing:
1. **Approval Status**: FDA/EMA/PMDA regulatory tiers
2. **Hub Drugs**: High-frequency drugs from top MLSS combinations
3. **Action Type**: Agonist vs Antagonist mechanism
4. **Hub Genes**: Target genes from network analysis
5. **Pathways**: Signaling cascades (Calcium, Incretin)
6. **Clinical Outcomes**: Therapeutic effects

## Requirements

```r
install.packages(c("tidyverse", "ggalluvial", "scales", "ggsci"))
```

## Usage

```r
source("fig6f_regulatory_mechanistic_sankey.R")
```

## Dependencies

- `fig6c_mlss_top15.csv` - Top MLSS combinations
- `fig6b_gene_centrality_summary.csv` - Hub genes from Fig 6b
- `FDA_Approved.csv` - FDA-approved drugs database
- `EMA_Approved.csv` - EMA-approved drugs database
- `PMDA_Approved.csv` - PMDA-approved drugs database
- `FDA-EMA-PMDA_Approved.csv` - Triple-approved drugs

## Output

**Figures**:
- `fig6f_regulatory_mechanistic_sankey.png/svg` (600 dpi, 39.1 × 10.3 cm)

**Data exports**:
- `fig6f_hub_drug_approval_profile.csv` - Hub drug regulatory status
- `fig6f_six_layer_complete_mapping.csv` - Complete flow data
- `fig6f_pathway_gene_connectivity.csv` - Pathway-gene relationships

## Regulatory Tiers

Based on DrugCentral 2023 databases:
- **Tier 1 (FDA-EMA-PMDA)**: Approved by all 3 major agencies (Confidence = 3.0)
- **Tier 2 (FDA-EMA)**: Approved by 2 agencies (Confidence = 2.0)
- **Tier 3 (FDA-Only)**: Single agency approval (Confidence = 1.0)
- **Not Approved**: No regulatory approval (Confidence = 0.3)

## Hub Drug Definition

Drugs appearing ≥4 times in top 15 MLSS combinations

## Flow Calculation

```
Flow = Frequency × Connection_Strength × Confidence
```

Where:
- **Frequency**: Appearances in top 15
- **Connection_Strength**: Gene connectivity (0-1)
- **Confidence**: Regulatory tier weight

## Color Scheme

### Pathways
- **Dark Blue**: Calcium Signaling
- **Red/Orange**: Incretin Signaling

### Approval Status
- **Red**: FDA-EMA-PMDA (Tier 1)
- **Light Blue**: FDA-EMA (Tier 2)
- **Grey**: FDA-Only (Tier 3)
- **Medium Grey**: Not Approved

## Interpretation

### Strong Flows (Thick)
- High drug frequency
- Strong gene connectivity
- High regulatory confidence

### Weak Flows (Thin)
- Lower frequency or connectivity
- Lower/no regulatory approval

### Cross-Pathway Connections
Reveal multi-pathway drug mechanisms

## Clinical Outcomes

1. **Enhanced insulin secretion**: Improved β-cell function
2. **β-cell preservation**: Reduced apoptosis
3. **Metabolic homeostasis**: Glucose regulation

## Data-Driven Validation

All regulatory assignments verified against:
- FDA Orange Book
- EMA Product Database
- PMDA Approved Drugs List
- DrugCentral 2023 (Jan 2026 snapshot)
