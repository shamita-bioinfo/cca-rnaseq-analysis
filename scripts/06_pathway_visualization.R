# ================================================
# Script: 06_pathway_visualization.R
# Author: Shamita
# Purpose: Visualize GO and KEGG pathway results
# ================================================

library(clusterProfiler)
library(ggplot2)
library(enrichplot)
# Load saved results from Script 05
go_results <- readRDS("/mnt/d/cca_rnaseq/results/go_results.rds")
kegg_results <- readRDS("/mnt/d/cca_rnaseq/results/kegg_results.rds")
# --- PLOT 1: GO Dot Plot ---
go_plot <- dotplot(go_results, showCategory = 20) +
  ggtitle("Top 20 GO Biological Processes — CCA vs Normal") +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

ggsave("/mnt/d/cca_rnaseq/figures/06_go_dotplot.png",
       plot = go_plot,
       width = 10, height = 9, dpi = 300)
message("GO dot plot saved!")
# --- PLOT 2: KEGG Dot Plot ---
kegg_plot <- dotplot(kegg_results, showCategory = 20) +
  ggtitle("Top 20 KEGG Pathways — CCA vs Normal") +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

ggsave("/mnt/d/cca_rnaseq/figures/06_kegg_dotplot.png",
       plot = kegg_plot,
       width = 10, height = 9, dpi = 300)

message("KEGG dot plot saved!")

# --- PLOT 3: GO Enrichment Map ---
# First we calculate which GO terms share genes with each other
go_results_sim <- pairwise_termsim(go_results)

# Then we draw the network
emap_plot <- emapplot(go_results_sim, showCategory = 30) +
  ggtitle("GO Term Network — CCA vs Normal") +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

ggsave("/mnt/d/cca_rnaseq/figures/06_go_emapplot.png",
       plot = emap_plot,
       width = 12, height = 10, dpi = 300)

message("GO enrichment map saved!")
# --- PLOT 4: Immune GO Pathways ---
# Filter for immune related terms
immune_terms <- c("immune", "inflammatory", "defense",
                  "cytokine", "leukocyte", "lymphocyte",
                  "interferon", "interleukin", "complement")

immune_filter <- grepl(paste(immune_terms, collapse = "|"),
                       go_results@result$Description,
                       ignore.case = TRUE)

go_immune <- go_results
go_immune@result <- go_results@result[immune_filter, ]

message(paste("Immune pathways found:", nrow(go_immune@result)))

# Immune dot plot
immune_dot <- dotplot(go_immune, showCategory = 20) +
  ggtitle("Immune GO Pathways — CCA vs Normal") +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

ggsave("/mnt/d/cca_rnaseq/figures/06_immune_dotplot.png",
       plot = immune_dot,
       width = 10, height = 9, dpi = 300)

message("Immune dot plot saved!")

# Immune enrichment map
go_immune_sim <- pairwise_termsim(go_immune)

immune_emap <- emapplot(go_immune_sim, showCategory = 20) +
  ggtitle("Immune GO Term Network — CCA vs Normal") +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

ggsave("/mnt/d/cca_rnaseq/figures/06_immune_emapplot.png",
       plot = immune_emap,
       width = 12, height = 10, dpi = 300)

message("Immune enrichment map saved!")
message("Script 06 complete!")
