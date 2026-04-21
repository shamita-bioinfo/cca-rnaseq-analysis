#=========================================
#script 07 
#Author - Shamita
#purpose - extract the raw file and then filter the cdk genes.
#          normalize the value 
#          plot heatmap with tumor vs normal separated.
#=========================================

library(TCGAbiolinks)   
library(DESeq2)         
library(pheatmap)
# this is the heatmap package
library(org.Hs.eg.db)
library(AnnotationDbi)
       
library(SummarizedExperiment)  
# needed to extract data from TCGA object
#extarcting the data from the script the raw one 
tcga_data <- readRDS("/mnt/d/cca_rnaseq/data/tcga_chol_rnaseq.rds")


# Load DESeq2 object
dds <- readRDS("/mnt/d/cca_rnaseq/results/dds.rds")

# Normalize using VST
vst_data <- vst(dds, blind = FALSE)

# Extract normalized matrix
vst_matrix <- assay(vst_data)

# Clean row names
ensembl_clean <- gsub("\\..*", "", rownames(vst_matrix))

# Convert to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_clean,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")


#which row have the cdk genes
#This gives you 60,660 TRUE/FALSE values
#TRUE where the gene symbol matches one of your CDK genes
#FALSE everywhere else

cdk_genes <- c("CDK1", "CDK4", "CCNB1", "CCNB2",
               "CCNA2", "CDKN2A", "CDKN2C",
               "CDC20", "CDC25A", "CDC25C")

# Filter for CDK genes
cdk_filter <- gene_symbols %in% cdk_genes
cdk_vst <- vst_matrix[cdk_filter, ]

#extracting data 
sample_info <- as.data.frame(colData(tcga_data))

#create annotation 
annotation_col <- data.frame(
  SampleType = sample_info$sample_type,
  row.names = rownames(sample_info)
)

#heatmap plotting
rownames(cdk_vst) <- gene_symbols[cdk_filter]
png("/mnt/d/cca_rnaseq/figures/07_cdk_heatmap.png",
    width = 10, height = 8, units = "in", res = 300)


pheatmap(
  cdk_vst,
  annotation_col = annotation_col,
  scale = "row",
  main = "CDK Cell Cycle Genes — CCA Tumor vs Normal",
  fontsize = 10,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100)
)


dev.off()
message("Script 07 complete — CDK heatmap saved!")

