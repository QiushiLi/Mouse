###written by Qiushi Li###
#setwd("/Users/Qiushi/Documents/github/mouse")

library(stats)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# Read the raw data
data <- fread("RNAseqExpressionProfile.txt", data.table = F, header = TRUE)

# Set the first column as row names
rownames(data) <- data[[1]]
data <- data[,-1]

################################################# PCA assessment #################################################

# remove genes with all zero counts
data <- data[rowSums(data) > 0, ]

# Perform PCA on the none all zero data
pca_result <- prcomp(t(data), scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x)
pca_data <- data.frame(Sample = rownames(pca_scores), PC1 = pca_scores[,1], PC2 = pca_scores[,2])

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
  geom_text(vjust = 1.5, color = "blue") +
  geom_point(aes(color = Sample), size = 3) +
  ggtitle("PCA of RNA-seq Data (All Samples)") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2") +
  theme_minimal()
print(pca_plot)

png(filename = "PCA.png", width = 800, height = 600, res = 100)
pca_plot
dev.off()

# remove Control2 from the dataset as per Dr. Wang suggested
expression_data <- data[, -which(colnames(data) == "Control2")]
expression_data[, 1:7] <- lapply(expression_data[, 1:7], as.numeric)

# filter out genes based on counts across all samples
filtered_expression_data <- expression_data[rowSums(expression_data) > 100, ]
hist(log1p(rowSums(filtered_expression_data)), breaks=50, main="Histogram of Log-transformed Sum of Counts",
     xlab="Log-transformed sum of counts", col="blue")
dim(filtered_expression_data)

# normalize the expression data with size factors calculated from geometric mean for each gene
geo_means <- exp(rowMeans(log(filtered_expression_data + 1)))
size_factors <- apply(filtered_expression_data, 2, function(column) median(column / geo_means))
normalized_expression_data <- sweep(filtered_expression_data, 2, size_factors, "/")

# set thresholds
mean_normalized_counts <- rowMeans(normalized_expression_data)
variance_normalized_counts <- apply(normalized_expression_data, 1, var)
final_filtered_data <- normalized_expression_data[mean_normalized_counts > 10 & variance_normalized_counts > 0.5, ]

# make comparison sets
control <- final_filtered_data[, 1:3]
knockout <- final_filtered_data[, 4:7]

# calculate log2 fold change
mean_control <- rowMeans(final_filtered_data[, 1:3])
mean_knockout <- rowMeans(final_filtered_data[, 4:7])
fold_change <- (mean_knockout + 1) / (mean_control + 1)
log2_fold_change <- log2(fold_change)


# Wilcoxon Rank-Sum Test using approximation due to ties
wilcox_results <- apply(normalized_expression_data, 1, function(x) {
  tryCatch({
    wilcox.test(x[1:3], x[4:7], exact = FALSE)
  }, error = function(e) NA
  )
})

# Extract p-values
p_values <- sapply(wilcox_results, function(x) x$p.value)

# Adjust p-values using the Benjamini-Hochberg method.The lowest adjusted P value is 0.2714.
#No significant genes can be identified with the original requirement.
#Switch to use P value as suggested by Dr. Wang.
#p_adjusted <- p.adjust(p_values, method = "BH")
#summary(p_adjusted)

diff_expr_genes <- which(abs(log2_fold_change) >= 1 & p_values <= 0.05)
diff_genes_data <- data.frame(
  Gene=rownames(final_filtered_data)[diff_expr_genes],
  Log2FoldChange=log2_fold_change[diff_expr_genes],
  PValue=p_values[diff_expr_genes]
)

# write to a CSV file
write.csv(diff_genes_data, file = "diff_genes_data.csv", row.names = FALSE)

################################################# Question 2 volcano #################################################

# data prep for volcano plot
volcano_data <- data.frame(GeneID=rownames(final_filtered_data), Log2FoldChange=log2_fold_change, PValue=p_values)
volcano_data$Significant <- as.factor(ifelse(volcano_data$PValue <= 0.05 & abs(volcano_data$Log2FoldChange) >= 1, ifelse(volcano_data$Log2FoldChange > 0, "Up regulated", "Down regulated"), "Not significant"))

# top 10 genes: 5 up-regulated and 5 down-regulated by P-value
top_genes <- volcano_data %>%
  filter(Significant != "Not significant") %>%
  arrange(PValue) %>%
  group_by(Significant) %>%
  slice_head(n = 5) %>%
  ungroup()

# volcano_plot
volcano_plot <- ggplot(volcano_data, aes(x = Log2FoldChange, y = -log10(PValue), color = Significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Up regulated" = "red", "Down regulated" = "blue", "Not significant" = "gray")) +
  labs(title = "Volcano Plot of RNA-seq Data",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(shape = 16))) +
  geom_label_repel(data = top_genes,
                   aes(label = GeneID),
                   box.padding = 0.35,
                   point.padding = 0.5,
                   size = 3)
print(volcano_plot)

png(filename = "volcano_plot.png", width = 800, height = 600, res = 100)
volcano_plot
dev.off()

################################################# Question 3 heatmap #################################################

# select the top 10 up-regulated/down-regulated genes
#top_up <- diff_genes_data[order(diff_genes_data$log2FoldChange, decreasing = TRUE),][1:10,]
top_up <- diff_genes_data[order(-diff_genes_data$Log2FoldChange), ][1:10, ]
top_down <- diff_genes_data[order(diff_genes_data$Log2FoldChange), ][1:10, ]
top_genes <- rbind(top_up, top_down)

# subset the normalized counts data for these top genes
top_genes_counts <- final_filtered_data[rownames(final_filtered_data) %in% top_genes$Gene,]

# Extract gene names from top_genes, match these genes with the row names of top_gene_counts
gene_order <- top_genes$Gene
ordered_indices <- match(gene_order, rownames(top_genes_counts))
ordered_top_genes_counts <- top_genes_counts[ordered_indices, ]

# log-transforme data
log_transformed_counts <- log2(ordered_top_genes_counts + 1)
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
upk_genes <- grep("^Upk", rownames(normalized_expression_data), value = TRUE)
upk_counts <- normalized_expression_data[upk_genes, ]

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

