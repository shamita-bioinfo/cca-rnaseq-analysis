# CCA RNA-seq Project Context — Shamita

## My Setup
- Dell laptop, Windows 11, WSL Ubuntu 22.04
- conda environment: cca_rnaseq
- R version 4.5.3
- Working directory: /mnt/d/cca_rnaseq/

## Project Structure
/mnt/d/cca_rnaseq/
├── data/      → TCGA-CHOL RNA-seq files (DONE)
├── scripts/   → R scripts
├── results/   → DESeq2 output
└── figures/   → plots

## GitHub
Username: shamita-bioinfo
Repo: github.com/shamita-bioinfo/cca-rnaseq-analysis

## Project 1 — CCA RNA-seq
Dataset: TCGA-CHOL (44 samples, 60660 genes)
Data saved: /mnt/d/cca_rnaseq/data/tcga_chol_rnaseq.rds

## Progress
- Step 1 DONE: Data downloaded and prepared
- Next: Step 2 — QC and preprocessing (Script 02)

## How to resume
1. Open WSL
2. conda activate cca_rnaseq
3. cd /mnt/d/cca_rnaseq## Detailed Progress Log
4. R
5. data <- readRDS("data/tcga_chol_rnaseq.rds")

### Day 1 — April 18, 2026
PROBLEMS SOLVED:
- Keyboard stuck period key (physical hardware issue)
- C: drive 95% full — moved everything to D: drive
- WSL kept shutting down — fixed with powercfg sleep settings

COMPLETED:
- Downloaded all 44 TCGA-CHOL files (185.94 MB) to /mnt/d/cca_rnaseq/data/
- GDCprepare successful:
  * 44 samples (35 Primary Tumor, 9 Solid Tissue Normal)
  * 60,660 genes
  * Saved as: /mnt/d/cca_rnaseq/data/tcga_chol_rnaseq.rds
  * Assays available: unstranded, stranded_first, stranded_second,
    tpm_unstrand, fpkm_unstrand, fpkm_uq_unstrand
  * Clinical + molecular subtype info added automatically

R SKILLS LEARNED TODAY:
- library() — loads a package
- <- — stores a value into a variable
- c() — combines values into a vector
- GDCquery() — tells GDC what data we want
- GDCdownload() — downloads the files
- GDCprepare() — merges all files into one R object
- saveRDS() — saves R object to disk
- dim() — shows dimensions (rows x columns)
- table() — counts categories
- assay() — extracts count matrix

NEXT IMMEDIATE STEP:
- Start Script 02: QC and preprocessing
- File: /mnt/d/cca_rnaseq/scripts/02_qc_preprocessing.R
- First thing: load the saved RDS file
  data <- readRDS("/mnt/d/cca_rnaseq/data/tcga_chol_rnaseq.rds")
- Then: extract count matrix using assay(data, "unstranded")
- Then: check for low quality samples
- Shamita learns by doing — explain every line, she types herself
- She is a beginner in R but strong in biology

### Day 1 continued — Script 02 complete
COMPLETED:
- QC check on all 44 samples
- All samples pass QC (min 25M > cutoff 21M)
- Keeping all 44 samples for analysis

R SKILLS LEARNED:
- source() — runs an entire script
- assay() — extracts counts from data container
- colSums() — adds up all rows in each column
- summary() — shows min, median, mean, max
- class() — tells you what type of object something is

NEXT: Script 03 — DESeq2 differential expression
- Compare 35 tumor vs 9 normal samples
- Find which genes are ON/OFF in cancer
