Overview
Mouse_RNAseq_pipeline.R, is a script to process and analyze RNA sequencing data from mouse study. No published RNA-seq software is required. The script includes several functions for handling, analyzing, and visualizing the data to identify differentially expressed genes under Control vs. Knockout.

Script Operations
1. Data Loading
* Purpose: Load RNA-seq data files including original gene expression counts.
2. Data Preprocessing
* Purpose: Prepare the data for downstream analysis.
* Description: This includes a PCA for assessing the original 8 samples, and remove the abnormal sample; filtering low expressed genes; normalizing data, and transforming counts for statistical analysis.
3. Differential Expression Analysis
* Purpose: Identify genes with expression levels differ significantly between control and knockout.
* Description: The script computes log2 fold changes to quantify expression differences and assesses the statistical significance of these changes using p-values.
4. Sorting Results
* Purpose: Organize and refine analysis results based on statistical metrics.
* Description: Genes are ordered according to their differential expression metrics, such as log2 fold change or P_value.
5. Visualization
* Purpose: Visualize the results for interpretation and presentation.
* Description: The script generate plots such as volcano plots, heatmaps and correlation matrix of Upk genes to display patterns of gene expression across samples.

Usage
* Environment: Ensure that R and necessary packages (e.g., dplyr, ggplot2) are installed.
* R Packages: List of R packages required to run the script
library(stats)library(data.table)library(dplyr)library(ggplot2)library(ggrepel)library(pheatmap) 