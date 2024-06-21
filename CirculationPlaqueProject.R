# README

  # AI was used to assist in the coding of this script.
  # Sample files are omitted and will not be available publicly.

# Packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
install.packages("gplots")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("dplyr")

library(gplots)
library(limma)
library(ggplot2)
library(reshape2)
library(dplyr)

# Batches
file_pathC1 <- ""
file_pathC2 <- ""
file_pathP1 <- ""
file_pathP2 <- ""

circSet1 <- read.delim(file_pathC1)
circSet1DF <- data.frame(circSet1)

circSet2 <- read.delim(file_pathC2)
circSet2DF <- data.frame(circSet2)

plaqSet1 <- read.delim(file_pathP1)
plaqSet1DF <- data.frame(plaqSet1)

plaqSet2 <- read.delim(file_pathP2)
plaqSet2DF <- data.frame(plaqSet2)

# Removing column X from cs1 and cs2
all(is.na(circSet1DF$X))
all(is.na(circSet2DF$X))

circSet1DF <- circSet1DF[, !names(circSet1DF) %in% "X"]
circSet2DF <- circSet2DF[, !names(circSet2DF) %in% "X"]

tail(colnames(circSet1DF))
tail(colnames(circSet2DF))

# Removing .sam.counts from cs1, setting genes as rownames for cs1, cs2, ps1, capitalizing aeN in ps1 and ps2
colnames(circSet1DF) <- gsub(".sam.counts", "", colnames(circSet1DF))

rownames(circSet1DF) <- circSet1DF[, 1]
rownames(circSet2DF) <- circSet2DF[, 1]
rownames(plaqSet1DF) <- plaqSet1DF[, 1]

colnames(plaqSet1DF) <- toupper(colnames(plaqSet1DF))
colnames(plaqSet2DF) <- toupper(colnames(plaqSet2DF))

circSet1DF <- subset(circSet1DF, select = -c(gene))
circSet2DF <- subset(circSet2DF, select = -c(gene))
plaqSet1DF <- subset(plaqSet1DF, select = -c(GENE))

# Checking for any missing values, conversion into matrices
any(is.na(circSet1DF))
any(is.na(circSet2DF))
any(is.na(plaqSet1DF))
any(is.na(plaqSet2DF))

matCirc1 <- as.matrix(circSet1DF)
matCirc2 <- as.matrix(circSet2DF)
matPlaq1 <- as.matrix(plaqSet1DF)
matPlaq2 <- as.matrix(plaqSet2DF)

dim(matCirc1)
dim(matCirc2)
dim(matPlaq1)
dim(matPlaq2)

# Filtering out patients and genes with high % of 0
filterPatientsGenes <- function(unFmatrix, thresholdPatients = 0.70, thresholdGenes = 0.90) {

  zero_proportionP <- colMeans(unFmatrix == 0)
  samples_to_remove <- names(zero_proportionP[zero_proportionP > thresholdPatients])
  filtered_patients_matrix <- unFmatrix[, !(colnames(unFmatrix) %in% samples_to_remove)]
  
  zero_proportionG <- rowMeans(filtered_patients_matrix == 0)
  genes_to_remove <- rownames(filtered_patients_matrix)[zero_proportionG > thresholdGenes]
  filtered_matrix <- filtered_patients_matrix[!(rownames(filtered_patients_matrix) %in% genes_to_remove), , drop = FALSE]
  
  return(filtered_matrix)
}

matCirc1F <- filterPatientsGenes(matCirc1)
matCirc2F <- filterPatientsGenes(matCirc2)
matPlaq1F <- filterPatientsGenes(matPlaq1)
matPlaq2F <- filterPatientsGenes(matPlaq2)

dim(matCirc1F)
dim(matCirc2F)
dim(matPlaq1F)
dim(matPlaq2F)

# Counts per million and quantile normalization of matrices
normalization <- function(counts_matrix) {
  
  total_counts <- colSums(counts_matrix)
  scaling_factor <- total_counts / 1e6
  normalized_counts <- t(t(counts_matrix) / scaling_factor)
  
  normalized_counts2 <- normalizeQuantiles(normalized_counts)
  
  return(normalized_counts2)
}

normCirc1 <- normalization(matCirc1F)
normCirc2 <- normalization(matCirc2F)
normPlaq1 <- normalization(matPlaq1F)
normPlaq2 <- normalization(matPlaq2F)

# Finding common genes and patients across all matrices
geneNamesC1 <- rownames(normCirc1)
geneNamesC2 <- rownames(normCirc2)
geneNamesP1 <- rownames(normPlaq1)
geneNamesP2 <- rownames(normPlaq2)

common_gene_names <- Reduce(intersect, list(geneNamesC1, geneNamesC2, geneNamesP1, geneNamesP2))
length(common_gene_names)

patientNamesCirculation <- c(colnames(normCirc1), colnames(normCirc2))
patientNamesPlaque <- c(colnames(normPlaq1), colnames(normPlaq2))

common_patient_names <- Reduce(intersect, list(patientNamesCirculation, patientNamesPlaque))
length(common_patient_names)

# Remove all genes and patients not in common
notCommonRemove <- function(normalMatrix) {
  
  notCommonMatrix <- normalMatrix[common_gene_names, , drop = FALSE]
  notCommonMatrix2 <- notCommonMatrix[, which(colnames(notCommonMatrix) %in% common_patient_names), drop = FALSE]
  
  return(notCommonMatrix2)
}

cleanC1 <- notCommonRemove(normCirc1)
dim(cleanC1)
cleanC2 <- notCommonRemove(normCirc2)
dim(cleanC2)
cleanP1 <- notCommonRemove(normPlaq1)
dim(cleanP1)
cleanP2 <- notCommonRemove(normPlaq2)
dim(cleanP2)

# Adding sets together, then ordering rows and columns of each matrix 
circulatingAll <- cbind(cleanC1, cleanC2)
ordered_column_namesC <- colnames(circulatingAll)[order(colnames(circulatingAll))]
ordered_circulating <- circulatingAll[, ordered_column_namesC, drop = FALSE]

plaqueAll <- cbind(cleanP1, cleanP2)
ordered_column_namesP <- colnames(plaqueAll)[order(colnames(plaqueAll))]
ordered_plaque <- plaqueAll[, ordered_column_namesP, drop = FALSE]

dim(ordered_circulating)
dim(ordered_plaque)

# --- Analysis ---

# Spearmann distribution analysis
spearman_correlations <- sapply(1:nrow(ordered_circulating), function(i) {
  cor(ordered_circulating[i, ], ordered_plaque[i, ], method = "spearman")
})

hist(spearman_correlations, breaks = 70, main = "", 
     xlab = "Spearman Correlation Coefficient", col = "darkorchid1", border = "black", ylim = c(0, 1200), xlim = c(-0.4, 0.8))
highest_value <- max(spearman_correlations)
lowest_value <- min(spearman_correlations)
abline(v = highest_value, col = "red", lwd = 3)
abline(v = lowest_value, col = "green", lwd = 3)
legend("topright", 
       legend = c(paste("Lowest Value:", round(lowest_value, 2), "; H1FOO"), paste("Highest Value:", round(highest_value, 3), "; DAZ4")), 
       col = c("green", "red"), 
       lty = 1, lwd = 2)

  # Positively correlated: 
    # DAZ4, DAZ1, DAZ2, DAZ3, ZFY
    # All Y chromosome linked..??

  # Negatively correlated:
    # H1FOO, PCDH7, GAS7, OR5L1, CCNH

# Differential expression analysis

t_test_result <- apply(cbind(ordered_circulating, ordered_plaque), 1, function(x) t.test(x[1:240], x[241:480]))
p_values <- sapply(t_test_result, function(x) x$p.value)
fold_change <- log2(rowMeans(ordered_circulating) / rowMeans(ordered_plaque))
neg_log_p_values <- -log10(p_values)

fold_change_threshold <- 1
p_value_threshold <- 0.05

significant_genes <- which(abs(fold_change) >= fold_change_threshold & p_values <= p_value_threshold)
upregulated_genes <- which(fold_change > fold_change_threshold & p_values <= p_value_threshold)
downregulated_genes <- which(fold_change < -fold_change_threshold & p_values <= p_value_threshold)
non_significant_genes <- setdiff(1:length(fold_change), significant_genes)

plot(fold_change, neg_log_p_values, 
     xlab = "Log Fold Change", ylab = "-log(P)", 
     main = "", col = "grey", pch = 20, cex = 0.7)
points(fold_change[upregulated_genes], neg_log_p_values[upregulated_genes], 
       col = ifelse(fold_change[upregulated_genes] > fold_change_threshold, "green", "grey"), pch = 20, cex = 0.7)
points(fold_change[downregulated_genes], neg_log_p_values[downregulated_genes], 
       col = ifelse(fold_change[downregulated_genes] < -fold_change_threshold, "orange", "grey"), pch = 20, cex = 0.7)
abline(h = -log10(0.05), col = "red", lty = 3)
abline(v = fold_change_threshold, col = "blue", lty = 3)
abline(v = -fold_change_threshold, col = "blue", lty = 3)

# PCA
combined_matrix <- cbind(ordered_circulating, ordered_plaque)
group <- factor(c(rep("Circulation", ncol(ordered_circulating)), rep("Plaque", ncol(ordered_plaque))))
pca_res <- prcomp(t(combined_matrix), scale. = TRUE)
pca_data <- data.frame(Sample = colnames(combined_matrix),
                       PC1 = pca_res$x[, 1],
                       PC2 = pca_res$x[, 2],
                       Group = group)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "", x = "Principal Component 1", y = "Principal Component 2") +
  theme(plot.title = element_text(hjust = 0.5))

loadings <- pca_res$rotation
pc1_loadings <- loadings[, 1]
top_genes_pc1 <- names(sort(abs(pc1_loadings), decreasing = TRUE)[1:10])
print(top_genes_pc1)

# Violin
circulationplus <- ordered_circulating + 1
plaqueplus <- ordered_plaque + 1

log_avg_expr1 <- apply(circulationplus, 1, function(x) log(mean(x)))
log_avg_expr2 <- apply(plaqueplus, 1, function(x) log(mean(x)))

ene_log_avg_expr1 <- data.frame(Gene = rownames(ordered_circulating), LogAvgExpression = log_avg_expr1, Condition = "Circulation")
ene_log_avg_expr2 <- data.frame(Gene = rownames(ordered_plaque), LogAvgExpression = log_avg_expr2, Condition = "Plaque")

combined_data <- rbind(ene_log_avg_expr1, ene_log_avg_expr2)

y_min <- floor(min(combined_data$LogAvgExpression))
y_max <- ceiling(max(combined_data$LogAvgExpression))

ggplot(combined_data, aes(x = Condition, y = LogAvgExpression, fill = Condition)) +
  geom_violin(scale = "width", trim = FALSE, width = 0.7, ) +
  scale_y_continuous(breaks = seq(y_min, y_max, by = 1)) +
  labs(title = "",
       x = "", y = "Log Average Expression") +
  theme_minimal() 

genes_in_condition1 <- combined_data %>%
  filter(Condition == "Circulation") %>%
  pull(LogAvgExpression)

genes_in_condition2 <- combined_data %>%
  filter(Condition == "Plaque") %>%
  pull(LogAvgExpression)

length(genes_in_condition1)
length(genes_in_condition2)
