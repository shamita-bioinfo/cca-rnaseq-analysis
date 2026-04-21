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

File: scripts/01_data_download.R
What it does: Downloads all 44 TCGA-CHOL RNA-seq samples
from GDC portal and merges them into one R object.

Key functions:
- GDCquery() — filters what data we want from TCGA
  (project=TCGA-CHOL, Transcriptome Profiling, STAR-Counts)
- GDCdownload() — physically downloads 44 files (185MB)
- GDCprepare() — merges all 44 files into one table
  (genes as rows, patients as columns)
- saveRDS() — saves merged object to disk as .rds file

Output: data/tcga_chol_rnaseq.rds
- 44 samples (35 Primary Tumor + 9 Solid Tissue Normal)
- 60,660 genes
- Assays: unstranded, stranded_first, stranded_second,
  tpm_unstrand, fpkm_unstrand, fpkm_uq_unstrand

---

NEXT IMMEDIATE STEP:
- Start Script 02: QC and preprocessing
- File: /mnt/d/cca_rnaseq/scripts/02_qc_preprocessing.R
- First thing: load the saved RDS file
  data <- readRDS("/mnt/d/cca_rnaseq/data/tcga_chol_rnaseq.rds")
- Then: extract count matrix using assay(data, "unstranded")
- Then: check for low quality samples


### Day 1 continued — Script 02 complete
COMPLETED:
- QC check on all 44 samples
- All samples pass QC (min 25M > cutoff 21M)
- Keeping all 44 samples for analysis

## Script 02 — Quality Control
File: scripts/02_qc_preprocessing.R
What it does: Checks all 44 samples have enough sequencing
reads to be reliable. Removes bad quality samples.

Key functions:
- assay() — extracts raw count matrix from data object
- colSums() — adds up all RNA counts per patient sample
  (total reads = library size)
- summary() — shows min, median, mean, max of library sizes

Logic: Any sample below 21 million reads is removed.
Low reads = sequencing failed = unreliable gene counts.

Result: ALL 44 samples passed (minimum was 25M reads)
No samples removed. All 44 kept for analysis.

---

## Script 03 — DESeq2 Differential Expression
File: scripts/03_deseq2_analysis.R
What it does: Compares gene expression between 35 tumor
vs 9 normal samples. Finds which genes are significantly
ON or OFF in cancer.

Key concepts:
- Normalization: DESeq2 corrects for different library
  sizes so samples are comparable
- log2FoldChange: positive = higher in tumor,
  negative = lower in tumor
- padj: adjusted p-value. < 0.05 = less than 5% chance
  this result happened by random chance
- Benjamini-Hochberg correction applied to avoid false
  positives when testing 60,660 genes simultaneously

Filters used: padj < 0.05 AND abs(log2FoldChange) > 1

Results:
- Total genes tested: 38,846
- Genes with valid results: 33,522
- Significant genes (padj<0.05): 14,429
- Upregulated in tumor: 4,418
- Downregulated in tumor: 6,747
- More downregulated = CCA cells losing bile duct identity

Output: results/dds.rds and results/res.rds

---

## Script 04 — Volcano Plot
File: scripts/04_volcano_plot.R
What it does: Visualizes all 33,522 tested genes as a
volcano plot. Each dot is one gene.

How to read the plot:
- X axis = log2FoldChange (left=down, right=up in tumor)
- Y axis = -log10(padj) (higher = more significant)
- Top right = upregulated AND significant (CDK1, CCNB1)
- Top left = downregulated AND significant (CDKN2A)
- Grey = not significant
- Red = significant in both fold change AND p-value

Key functions:
- EnhancedVolcano() — draws the volcano plot
- ggsave() — saves as PNG (more reliable than png() in WSL)

Output: figures/volcano_plot.png

---

## Script 05 — Pathway Enrichment Analysis
File: scripts/05_pathway_enrichment.R
What it does: Takes 11,165 significant genes and groups
them into biological pathways. Asks which processes are
most disrupted in CCA.

Steps in script:
1. Load res.rds from Script 03
2. Filter significant genes (padj<0.05, log2FC>1)
3. Remove version numbers from Ensembl IDs
   (ENSG00000052802.13 → ENSG00000052802)
4. Convert Ensembl IDs → Entrez IDs using bitr()
   (clusterProfiler needs Entrez IDs to query databases)
5. Run GO enrichment (enrichGO) — broad biological functions
6. Run KEGG enrichment (enrichKEGG) — specific named pathways
   Note: SSL fix required — library(httr) +
   set_config(config(ssl_verifypeer=0L)) before KEGG
7. Save both results

25.54% of IDs failed to map — normal, these are ncRNAs
and pseudogenes without Entrez IDs. 75% mapped = sufficient.

Top GO Results:
- Small molecule catabolic process (264 genes)
- Fatty acid metabolic process (247 genes)
- Alpha-amino acid metabolic process (144 genes)
Meaning: Normal bile duct metabolism completely disrupted

Top KEGG Results:
- Complement and coagulation cascades (71 genes) → immune evasion
- Bile secretion (65 genes) → CCA signature, bile duct identity lost
- Fatty acid degradation (37 genes) → metabolic reprogramming
- Tryptophan metabolism/IDO1/IDO2 (36 genes) → immune suppression
- Cell cycle hsa04110 (64 genes) → CDK1, CCNB1, CDKN2A found

Output: results/go_results.rds and results/kegg_results.rds

---

## Script 06 — Pathway Visualization
File: scripts/06_pathway_visualization.R
What it does: Creates 5 publication-quality figures from
GO and KEGG pathway results.

Plots created:
1. GO dotplot — top 20 GO biological processes
   X axis = GeneRatio (fraction of pathway covered)
   Dot size = number of genes
   Dot colour = significance (red=most significant)

2. KEGG dotplot — top 20 KEGG pathways
   Same format as GO dotplot

3. GO enrichment map (emapplot) — network diagram
   pairwise_termsim() calculates gene overlap between terms
   Connected nodes = pathways sharing genes
   Clusters = related biological processes

4. Immune GO dotplot — filtered for immune-related terms
   grepl() used to search GO descriptions for keywords:
   immune, inflammatory, cytokine, complement etc.
   381 immune pathways found in CCA

5. Immune GO network — emapplot of immune terms only

Key functions:
- dotplot() — draws enrichment dot plots
- emapplot() — draws pathway network diagrams
- pairwise_termsim() — calculates similarity between terms
- grepl() — searches text for pattern matching

Output: figures/06_go_dotplot.png
        figures/06_kegg_dotplot.png
        figures/06_go_emapplot.png
        figures/06_immune_dotplot.png
        figures/06_immune_emapplot.png

---

## Script 06b — Biomarker Candidate Discovery
File: scripts/06b_biomarker_candidate.R
What it does: Narrows 38,846 genes down to top 50
strongest biomarker candidates using strict filters.

Filtering pipeline:
- All genes: 38,846
- Remove NA padj: 33,522
- padj < 0.05: 14,429
- abs(log2FC) > 2 (4x change): 6,524
- Known gene symbols only: 4,745
- Top 25 up + top 25 down = 50 candidates

Direction split:
- 2,570 upregulated
- 3,954 downregulated

Top upregulated candidates:
- USH2A (log2FC +6.5) — highest fold change
- KCNN2 (log2FC +5.5) — ion channel, cancer survival
- MSMO1 (log2FC +3.8) — cholesterol biosynthesis
- DHODH (log2FC +3.5) — active drug target in cancer
- HDAC6 (log2FC +2.1) — epigenetic regulator

Top downregulated candidates:
- EPS8L3 (log2FC -7.1) — biggest loss in dataset
- PRSS16 (log2FC -5.4) — immune related protease
- ESRP1 (log2FC -5.2) — RNA splicing, loss = invasion
- EPCAM (log2FC -4.7) — epithelial identity marker
- PKM (log2FC -4.1) — Warburg effect confirmed

Key functions used:
- as.data.frame() — converts DESeq2 results to table
- !is.na() — removes rows with missing values
- abs() — absolute value for both + and - fold changes
- order() — sorts table by a column
- ifelse() — if/else logic to label UP or DOWN
- mapIds() — converts Ensembl IDs to gene symbols
- rbind() — stacks two tables together
- write_xlsx() — saves table as Excel file
- factor() — fixes gene order on plot axis
- coord_flip() — flips axes for horizontal bar chart

Output: results/top_biomarker_candidates.xlsx
        figures/06b_biomarker_candidates.png

---

## Overall


## Scripts Completed
- Script 01 ✅ Data download (TCGA-CHOL, 44 samples, 60660 genes)
- Script 02 ✅ QC (all 44 passed, min 25M reads)
- Script 03 ✅ DESeq2 (4418 up, 6747 down, 14429 significant)
- Script 04 ✅ Volcano plot (figures/volcano_plot.png)
- Script 05 ✅ GO + KEGG pathway enrichment (with SSL fix)
- Script 06 ✅ Pathway figures (GO dotplot, KEGG dotplot, emapplot, immune plots)
- Script 06b ✅ Biomarker candidates (top 50 genes, Excel saved)

## Key Findings
CDK1, CCNB1, CDK4 upregulated (cell cycle drivers)
CDKN2A, CDKN2C downregulated (tumor suppressors lost)
Bile secretion disrupted (CCA signature pathway)
Complement/coagulation = immune evasion
EPCAM downregulated = epithelial identity lost
USH2A highest upregulated (log2FC 6.5)
EPS8L3 highest downregulated (log2FC -7.1)

## Figures Saved
figures/volcano_plot.png
figures/06_go_dotplot.png
figures/06_kegg_dotplot.png
figures/06_go_emapplot.png (if generated)
figures/06_immune_dotplot.png
figures/06_immune_emapplot.png

## Results Saved
results/res.rds
results/dds.rds
results/go_results.rds
results/kegg_results.rds
results/cdk_genes.csv
results/biomarker_candidates.xlsx (top 50 genes)

## Remaining Steps
Script 07 → CDK focused heatmap
Script 08 → Immune deconvolution (TIMER2.0 web tool)
Script 09 → Survival analysis
Script 10 → Validation in GSE107943
Paper writing → bioRxiv preprint
