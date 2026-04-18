# ============================================================
# Script: 01_download_data.R
# Project: CCA RNA-seq (TCGA-CHOL)
# Author: Shamita
# Purpose: Download TCGA-CHOL RNA-seq data robustly
# ============================================================

library(TCGAbiolinks)

# --- Define where data will be saved ---
data_dir <- "~/projects/cca_rnaseq/data"

# --- Build the query ---
# This tells GDC exactly what data we want:
# TCGA-CHOL = cholangiocarcinoma project
# STAR - Counts = raw read counts (what DESeq2 needs)
query <- GDCquery(
  project          = "TCGA-CHOL",
  data.category    = "Transcriptome Profiling",
  data.type        = "Gene Expression Quantification",
  workflow.type    = "STAR - Counts"
)

# --- Download in small chunks (more stable on slow connections) ---
# files.per.chunk = 6 means: download 6 files at a time
# This way, if connection drops, you don't lose everything
GDCdownload(
  query,
  method         = "api",
  directory      = data_dir,
  files.per.chunk = 6
)

# --- Prepare: merge all files into one R object ---
# GDCprepare reads all downloaded files and combines them
# into a SummarizedExperiment object (a standard Bioconductor format)
data <- GDCprepare(
  query,
  directory = data_dir
)

# --- Save the object so you never need to download again ---
saveRDS(data, file = paste0(data_dir, "/tcga_chol_rnaseq.rds"))

message("✅ Download and preparation complete!")
message(paste("Samples downloaded:", ncol(data)))
message(paste("Genes in dataset:", nrow(data)))
