# CCA RNA-seq Analysis

Differential expression analysis of cholangiocarcinoma (CCA) tumor vs adjacent normal tissue using TCGA-CHOL RNA-seq data.

## Biological Question
Which genes and pathways are dysregulated in cholangiocarcinoma? With focus on:
- CDK/cell cycle pathway dysregulation
- Immune microenvironment characterization

## Dataset
- Primary: TCGA-CHOL (38 tumor, 9 normal)
- Validation: GSE107943 (GEO)

## Tools
R · DESeq2 · clusterProfiler · EnhancedVolcano · TIMER2.0

## Project Structure
data/      → raw TCGA-CHOL count files
scripts/   → R analysis scripts
results/   → DESeq2 output tables
figures/   → volcano plots, heatmaps, pathway plots

## Author
Shamita | github.com/shamita-bioinfo
