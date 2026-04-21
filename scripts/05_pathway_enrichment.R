# ================================================
# Script: 05_pathway_enrichment.R
# Author: Shamita
# Purpose: Identify which biological processes are
#          significantly disrupted in CCA using
#          pathway enrichment analysis (GO + KEGG)
# ================================================
# Load required libraries
library(SummarizedExperiment)
library(DESeq2)
library(httr)
library(clusterProfiler)
library(org.Hs.eg.db)
# Load DESeq2 results from Script 03
res<- readRDS("/mnt/d/cca_rnaseq/results/res.rds")
# Filter significant genes
# padj < 0.05 AND fold change > 2x
sig_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
# Remove version numbers from Ensembl IDs
# e.g. ENSG00000052802.13 → ENSG00000052802
gene_ids <- gsub("\\..*", "", rownames(sig_genes))
# Convert Ensembl IDs to Entrez IDs
# clusterProfiler needs Entrez IDs for pathway analysis
gene_entrez <- bitr(gene_ids,
                    fromType = "ENSEMBL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
# GO enrichment analysis
# Finding which biological processes are disrupted
go_results <- enrichGO(gene = gene_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)
# KEGG pathway analysis
# hsa = Homo sapiens (human)
# This will show us cell cycle, immune pathways etc.
set_config(config(ssl_verifypeer = 0L))
kegg_results <- enrichKEGG(gene = gene_entrez$ENTREZID,
                            organism = "hsa",
                            pvalueCutoff = 0.05)
# Save results and print summary
saveRDS(go_results, file = "/mnt/d/cca_rnaseq/results/go_results.rds")
saveRDS(kegg_results, file = "/mnt/d/cca_rnaseq/results/kegg_results.rds")
# Print top GO results
print(head(go_results))

# Print top KEGG results  
print(head(kegg_results))

message("Pathway Enrichment Analysis is Complete!")
