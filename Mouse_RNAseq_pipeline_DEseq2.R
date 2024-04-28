#setwd("/Users/Qiushi/Documents/github/mouse/qiushi")
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# read the RNA-seq data
expression_data <- read.table("RNAseqExpressionProfile.txt", header = TRUE, row.names = 1, sep="\t", check.names = FALSE)

# remove genes with all zero counts
expression_data <- expression_data[rowSums(expression_data) > 0, ]

# PCA to assess dataset
dds_pca <- DESeqDataSetFromMatrix(countData = expression_data, colData = data.frame(condition = rep(c("Control1", "Control2", "Control3", "Control4", "KO1", "KO2", "KO3", "KO4"), each=1), row.names = colnames(expression_data)), design = ~ 1)
dds_pca <- estimateSizeFactors(dds_pca)
normalized_counts <- counts(dds_pca, normalized = TRUE)
pca_data <- prcomp(t(normalized_counts), scale. = TRUE)
pca_results <- as.data.frame(pca_data$x)

# PCAplot
pca_plot <- ggplot(pca_results, aes(x = PC1, y = PC2, color = rep(c("Control1", "Control2", "Control3", "Control4", "KO1", "KO2", "KO3", "KO4"), each=1))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of RNA-seq samples", x = "Principal Component 1", y = "Principal Component 2", colour = "Category") +
  scale_color_manual(values = rep(c("blue", "red"), each = 4))
print(pca_plot)

png(filename = "PCA.png", width = 800, height = 600, res = 100)
pca_plot
dev.off()

################################################# Question 1 significant genes #################################################

# remove Control2 from the dataset
expression_data <- expression_data[, -which(colnames(expression_data) == "Control2")]

# filter out genes with zero counts across all samples again
filtered_data <- expression_data[rowSums(expression_data) > 0, ]

# differential expression analysis, normalization with DESeq2
condition <- factor(c(rep("Control", 3), rep("KO", 4)))
coldata <- data.frame(row.names = colnames(data), condition)
dds <- DESeqDataSetFromMatrix(countData = filtered_data,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)

# Wilcoxon test and adjust p-values
res <- results(dds, contrast = c("condition", "KO", "Control"), pAdjustMethod = "fdr")
res <- na.omit(res)

# select significant genes with abs(log2FC >= 1) and adj P-value <= 0.05
significant_genes <- res[abs(res$log2FoldChange) >= 1 & res$padj <= 0.05, ]

# select only required columns
output_genes <- as.data.frame(significant_genes)[, c("log2FoldChange", "padj")]
output_genes$GeneID <- rownames(output_genes)

# write to a CSV file
write.csv(output_genes, file = "DE_genes.csv", row.names = FALSE)

################################################# Question 2 volcano #################################################
# data prep

volcano_data <- as.data.frame(res)
volcano_data$GeneID <- rownames(volcano_data)
volcano_data$Significant <- as.factor(ifelse(volcano_data$padj <= 0.05 & abs(volcano_data$log2FoldChange) >= 1, ifelse(volcano_data$log2FoldChange > 0, "Up regulated", "Down regulated"), "Not significant"))

# top 10 genes: 5 up-regulated and 5 down-regulated by adjusted P-value
top_genes <- volcano_data %>%
  filter(Significant != "Not significant") %>%
  arrange(padj) %>%
  group_by(Significant) %>%
  slice_head(n = 5) %>%
  ungroup()

# plot
volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Up regulated" = "red", "Down regulated" = "blue", "Not significant" = "gray")) +
  labs(title = "Volcano Plot of RNA-seq Data",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  geom_label_repel(data = top_genes,
                   aes(label = GeneID),
                   box.padding = 0.35,
                   point.padding = 0.5,
                   size = 3)
print(volcano_plot)

png(filename = "volcano.png", width = 800, height = 600, res = 100)
volcano_plot
dev.off()

################################################# Question 3 heatmap #################################################
# convert res to a regular data frame, add a column of gene names
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# select the top 10 up-regulated/down-regulated genes
top_up <- res_df[order(res_df$log2FoldChange, decreasing = TRUE),][1:10,]
top_down <- res_df[order(res_df$log2FoldChange),][1:10,]
top_genes <- rbind(top_up, top_down)

# subset the normalized counts data for these top genes
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
top_genes_counts <- normalized_counts[rownames(normalized_counts) %in% top_genes$gene,]

# gene order matches the order of the top genes list！
top_genes_counts <- top_genes_counts[top_genes$gene,]

# log-transforme data
log_transformed_counts <- log2(top_genes_counts + 1)
expression_heatmap = pheatmap(log_transformed_counts, 
         cluster_rows = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Expression Heatmap of Top 20 Differentially Expressed Genes",
         labRow = top_genes$gene,
         labCol = colnames(log_transformed_counts)) 
print(expression_heatmap)

png(filename = "Expression heatmap.png", width = 800, height = 600, res = 100)
expression_heatmap
dev.off()

################################################# Question 4 correlation #################################################
# get Upk genes
upk_genes <- grep("^Upk", rownames(normalized_counts), value = TRUE)
upk_counts <- normalized_counts[upk_genes, ]
  
# calculate correlation matrix
correlation_matrix <- cor(t(upk_counts), use = "pairwise.complete.obs", method = "pearson")
rownames(correlation_matrix) <- upk_genes
colnames(correlation_matrix) <- upk_genes

write.csv(correlation_matrix, file = "correlation_matrixs.csv", )

# correlation plot
correlation_heatmap = pheatmap(correlation_matrix, 
         main = "Correlation Matrix of Upk Gene Expression", 
         color = colorRampPalette(c("blue", "white", "red"))(n = 50),
         show_rownames = TRUE,
         show_colnames = TRUE)
print(correlation_heatmap)

png(filename = "Upk_genes_correlation.png", width = 800, height = 600, res = 100)
correlation_heatmap 
dev.off()
