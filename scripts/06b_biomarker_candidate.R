# ================================================
# Script: 06b_biomarker_candidates.R
# Author: Shamita
# Purpose: Identify top biomarker candidate genes
#          from DESeq2 results based on statistical
#          strength and fold change
# ================================================

library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load DESeq2 results from Script 03
res <- readRDS("/mnt/d/cca_rnaseq/results/res.rds")

# Convert to a regular dataframe so we can work with it easily
res_df <- as.data.frame(res)

# Add the Ensembl gene IDs as their own column
# Right now they are row names — we want them as a proper column
res_df$ensembl_id <- rownames(res_df)

# Remove version numbers from Ensembl IDs
# e.g. ENSG00000052802.13 → ENSG00000052802
res_df$ensembl_clean <- gsub("\\..*", "", res_df$ensembl_id)

message(paste("Total genes in results:", nrow(res_df)))
# --- STEP 2: Filter for strong candidates ---

# Remove genes with missing values
# Some genes have NA for padj — these failed statistical testing
res_clean <- res_df[!is.na(res_df$padj), ]
message(paste("After removing NA:", nrow(res_clean)))

# Filter 1 — statistical confidence
# padj < 0.05 means less than 5% chance this is a false finding
res_sig <- res_clean[res_clean$padj < 0.05, ]
message(paste("After significance filter:", nrow(res_sig)))

# Filter 2 — fold change strength
# log2FoldChange > 2 means 4x higher in tumor
# log2FoldChange < -2 means 4x lower in tumor
res_strong <- res_sig[abs(res_sig$log2FoldChange) > 2, ]
message(paste("After fold change filter:", nrow(res_strong)))

# Label each gene as UP or DOWN in tumor
res_strong$direction <- ifelse(res_strong$log2FoldChange > 0, 
                               "Upregulated", 
                               "Downregulated")

# How many up vs down?

table(res_strong$direction)
table(res_strong$direction)
# --- STEP 3: Add real gene names ---
# Right now we only have Ensembl IDs like ENSG00000012048
# We need gene symbols like BRCA1, CDK1, TP53
# mapIds() looks up the gene symbol for each Ensembl ID

res_strong$gene_symbol <- mapIds(org.Hs.eg.db,
                                  keys = res_strong$ensembl_clean,
                                  column = "SYMBOL",
                                  keytype = "ENSEMBL",
                                  multiVals = "first")

# Remove any genes where symbol lookup failed
res_named <- res_strong[!is.na(res_strong$gene_symbol), ]

message(paste("Genes with known symbols:", nrow(res_named)))

# --- STEP 4: Rank by significance ---
# Sort by padj smallest to largest
# Smallest padj = most statistically confident finding
res_ranked <- res_named[order(res_named$padj), ]

# Extract top 25 upregulated genes
top_up <- res_ranked[res_ranked$direction == "Upregulated", ]
top_up <- head(top_up, 25)

# Extract top 25 downregulated genes
top_down <- res_ranked[res_ranked$direction == "Downregulated", ]
top_down <- head(top_down, 25)

# Combine into one table
top_candidates <- rbind(top_up, top_down)

# Print a clean summary to screen
print(top_candidates[, c("gene_symbol", 
                          "log2FoldChange", 
                          "padj", 
                          "direction")])
# --- STEP 5: Save results table ---
# Install writexl if you don't have it
if (!requireNamespace("writexl", quietly = TRUE)) {
    install.packages("writexl")
}
library(writexl)

# Save top candidates as Excel file
write_xlsx(top_candidates[, c("gene_symbol",
                               "log2FoldChange",
                               "padj",
                               "direction")],
           "/mnt/d/cca_rnaseq/results/top_biomarker_candidates.xlsx")

message("Excel file saved!")

# --- STEP 6: Bar plot of top candidates ---
# Take top 15 up and top 15 down for cleaner visual
plot_up <- head(top_up, 15)
plot_down <- head(top_down, 15)
plot_data <- rbind(plot_up, plot_down)

# Sort by fold change for visual clarity
plot_data <- plot_data[order(plot_data$log2FoldChange), ]
plot_data$gene_symbol <- factor(plot_data$gene_symbol, 
                                 levels = plot_data$gene_symbol)

# Draw the bar plot
candidate_plot <- ggplot(plot_data, 
                         aes(x = gene_symbol, 
                             y = log2FoldChange, 
                             fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Upregulated" = "#E64B35",
                                "Downregulated" = "#4DBBD5")) +
  labs(title = "Top Biomarker Candidates — CCA vs Normal",
       x = "Gene",
       y = "Log2 Fold Change",
       fill = "Direction") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

ggsave("/mnt/d/cca_rnaseq/figures/06b_biomarker_candidates.png",
       plot = candidate_plot,
       width = 10, height = 9, dpi = 300)

message("Bar plot saved!")
message("Script 06b complete!")
