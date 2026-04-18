# ================================================
# Script: 04_volcano_plot.R
# Author: Shamita
# Purpose: Filter 16433 genes to find most meaningful
#          using fold change and statistical significance
#          then visualize using volcano plot
# ================================================
# Load required libraries
library(TCGAbiolinks) 
library(SummarizedExperiment) 
library(DESeq2)  
library(ggplot2)
library(EnhancedVolcano)
# Load DESeq2 results from Script 03
res <- readRDS("/mnt/d/cca_rnaseq/results/res.rds")

# Filter significant genes
# padj < 0.05 = statistically significant
# abs(log2FoldChange) > 1 = at least 2x change
sig_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
# How many significant genes did we find?
nrow(sig_genes)
# Quick summary of results
summary(res)
# Create volcano plot
# Each dot = one gene
# Red dots = significant and highly changed = our candidates!
# Save volcano plot

# Convert res to dataframe first — WSL needs this!
res_df <- as.data.frame(res)


# Save volcano plot
# Store plot in variable p
p <- EnhancedVolcano(res_df,
    lab = rownames(res_df),
    x = 'log2FoldChange',
    y = 'padj',
    title = 'CCA Tumor vs Normal',
    pCutoff = 0.05,
    FCcutoff = 1)

# Save using ggsave — more reliable than png() in WSL
ggsave("/mnt/d/cca_rnaseq/figures/volcano_plot.png",
       plot = p,
       width = 12,
       height = 8,
       dpi = 150)
