# ================================================
# Script: 03_deseq2_analysis.R
# Author: Shamita
# Purpose: Compare gene expression between 35 tumor
#          and 9 normal samples in TCGA-CHOL to find
#          differentially expressed genes
# ================================================

# Load required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)

# Load saved TCGA-CHOL data
data <- readRDS("/mnt/d/cca_rnaseq/data/tcga_chol_rnaseq.rds")

# Extract sample information (tumor vs normal labels)
sample_info <- colData(data)

# Check sample type labels — tumor vs normal
print(sample_info$sample_type)
# Create DESeq2 object
# dds = DESeq2 Dataset — this is the object DESeq2 works with
# design = ~ sample_type tells DESeq2 to compare based on tumor vs normal
dds <- DESeqDataSet(data, design = ~ sample_type)
# Filter low expression genes
# Remove genes with less than 10 total counts across all samples
# This reduces noise and speeds up analysis
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
# Run DESeq2 — this does all the statistical heavy lifting
# Normalizes counts, fits model, tests every gene
dds <- DESeq(dds)
# Extract results
# res contains fold change and statistics for every gene
res <- results(dds)
# Summarize results — shows how many genes are up/down regulated
summary(res)
# Save results — so we never need to rerun DESeq2 again!
saveRDS(dds, file = "/mnt/d/cca_rnaseq/results/dds.rds")
saveRDS(res, file = "/mnt/d/cca_rnaseq/results/res.rds")
message("DESeq2 analysis complete!")
