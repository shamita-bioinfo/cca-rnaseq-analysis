
# ================================================
# Script: 02_qc_preprocessing.R
# Author: Shamita
# Purpose: QC and preprocessing of TCGA-CHOL data
# ================================================
# Step 1: Load required libraries
# Libraries are like toolboxes — each one adds new functions to R
library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
# Step 2: Load the saved data
# readRDS loads our saved object back into R memory
data <- readRDS("/mnt/d/cca_rnaseq/data/tcga_chol_rnaseq.rds")
counts <- assay(data, "unstranded")

# Step 4: Calculate total counts per sample
# colSums() adds up all 60,660 gene rows for each patient column
# Low total counts = poor quality sample = like a contaminated tube
sample_counts <- colSums(counts)
# Let's see the total counts for each sample
print(sample_counts)
summary(sample_counts)
# QC DECISION:
# Input: 44 patient samples, 60660 genes
# We calculated total RNA counts per sample using colSums()
# Summary showed:
#   Min    = 25 million
#   Median = 43 million
#   Max    = 63 million
# Rule: remove samples with counts less than 50% of median
# 50% of 43 million = 21 million
# Our minimum (25 million) is above 21 million
# CONCLUSION: All 44 samples pass QC — none removed!
