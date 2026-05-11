library(dplyr)

split_means <- lapply(1:10, function(i) {
  file.list <- list.files(
    paste0("pseudobulk_update/split_", i, "/data.splits/"),
    full.names = TRUE,
    pattern = "X_test"
  )
  
  # If no files are found for a split, return NA
  if (length(file.list) == 0) {
    return(data.frame(
      Split = i,
      Mean_Prop_Above_0.7_pct = NA_real_
    ))
  }
  
  matrix.list <- lapply(file.list, function(x) {
    read.csv(x, row.names = 1)
  })
  names(matrix.list) <- gsub("X_test_|.csv", "", basename(file.list))
  
  multicollinearity_results <- lapply(names(matrix.list), function(mat_name) {
    mat <- matrix.list[[mat_name]]
    
    # Handle cases with <2 features (no pairs)
    if (ncol(mat) < 2) {
      return(data.frame(
        Matrix_ID = mat_name,
        Num_Features = ncol(mat),
        Total_Feature_Pairs = 0,
        Max_Abs_Cor = NA_real_,
        Median_Abs_Cor = NA_real_,
        Pairs_Above_0.7 = 0,
        Prop_Above_0.7 = NA_real_
      ))
    }
    
    cor_matrix <- cor(mat, method = "spearman")
    upper_tri <- cor_matrix[upper.tri(cor_matrix)]
    abs_cors <- abs(upper_tri)
    
    data.frame(
      Matrix_ID = mat_name,
      Num_Features = ncol(mat),
      Total_Feature_Pairs = length(abs_cors),
      Max_Abs_Cor = max(abs_cors, na.rm = TRUE),
      Median_Abs_Cor = median(abs_cors, na.rm = TRUE),
      Pairs_Above_0.7 = sum(abs_cors > 0.7, na.rm = TRUE),
      Prop_Above_0.7 = sum(abs_cors > 0.7, na.rm = TRUE) / length(abs_cors)
    )
  })
  
  final_summary_df <- bind_rows(multicollinearity_results)
  
  data.frame(
    Split = i,
    Mean_Prop_Above_0.7_pct = mean(final_summary_df$Prop_Above_0.7, na.rm = TRUE) * 100
  )
})

all_split_means <- bind_rows(split_means)

# Save results
write.csv(
  all_split_means,
  "mean_prop_above_0.7_by_split.csv",
  row.names = FALSE
)

all_split_means


##### Checking if genes are coexpressed with confounder genes
library(Seurat)
library(edgeR)

features <- read.csv('results_update/all_features.csv')


pbmc <- readRDS('pbmc_female.control_managed.RDS')

CD4.subset <- subset(pbmc, cell_type_detailed == "CD4-positive, alpha-beta T cell" & disease == 'disease')
CD4.bulk <- AggregateExpression(CD4.subset, slot='counts', group.by='ind_cov')$COMBAT_LogNorm
non_zero_genes <- rowSums(CD4.bulk) > 0
CD4.bulk <- CD4.bulk[non_zero_genes, ]

CD4_cor_matrix <- cor(as.matrix(t(CD4.bulk)), method = "spearman")
upper_tri <- CD4_cor_matrix[upper.tri(CD4_cor_matrix)]

CD4_features <- subset(features, celltype == 'CD4_positive,_alpha_beta_T_cell.chrX')$feature


cd4_cor_subset <- CD4_cor_matrix[rownames(CD4_cor_matrix) %in% CD4_features, , drop = FALSE]

high_cor_indices <- which(abs(cd4_cor_subset) > 0.7, arr.ind = TRUE)

# 3. Build a clean data frame with the results
high_cor_table <- data.frame(
CD4_Predictor = rownames(cd4_cor_subset)[high_cor_indices[, 1]],
Correlated_Gene = colnames(cd4_cor_subset)[high_cor_indices[, 2]],
Abs_Correlation_Score = cd4_cor_subset[high_cor_indices]
)

# 4. Remove self-correlations (e.g., IL2RG correlating 1.0 with IL2RG)
high_cor_table <- high_cor_table[high_cor_table$CD4_Predictor != high_cor_table$Correlated_Gene, ]

# 5. Sort the table from highest correlation to lowest
high_cor_table <- high_cor_table[order(-abs(high_cor_table$Abs_Correlation_Score)), ]


CD8.subset <- subset(pbmc, cell_type_detailed == "CD8-positive, alpha-beta T cell" & disease == 'disease')
CD8.bulk <- AggregateExpression(CD8.subset, slot='counts', group.by='ind_cov')$COMBAT_LogNorm
non_zero_genes <- rowSums(CD8.bulk) > 0
CD8.bulk <- CD8.bulk[non_zero_genes, ]
CD8_cor_matrix <- cor(as.matrix(t(CD8.bulk)), method = "spearman")
upper_tri <- CD8_cor_matrix[upper.tri(CD8_cor_matrix)]
CD8_features <- subset(features, celltype == 'CD8_positive,_alpha_beta_T_cell.chrX')$feature
cd8_cor_subset <- CD8_cor_matrix[rownames(CD8_cor_matrix) %in% CD8_features, , drop = FALSE]
high_cor_indices <- which(abs(cd8_cor_subset) > 0.7, arr.ind = TRUE)
high_cor_table_cd8 <- data.frame(
  CD8_Predictor = rownames(cd8_cor_subset)[high_cor_indices[, 1]],
  Correlated_Gene = colnames(cd8_cor_subset)[high_cor_indices[, 2]],
  Abs_Correlation_Score = cd8_cor_subset[high_cor_indices]
)
high_cor_table_cd8 <- high_cor_table_cd8[high_cor_table_cd8$CD8_Predictor != high_cor_table_cd8$Correlated_Gene, ]
high_cor_table_cd8 <- high_cor_table_cd8[order(-abs(high_cor_table_cd8$Abs_Correlation_Score)), ]


#####
model_results <- read.csv('results_update/Supplementary_Table_1.csv')
features <- read.csv('results_update/all_features.csv')

cMono.chrX <- list()
for(i in 1:10){
    cMono.chrX[[i]] <- read.csv(paste0('pseudobulk_update/split_', i, '/combined/ensemble/metrics_Classical_monocyte.chrX.csv'))
}
cMono.chrX <- do.call(rbind, cMono.chrX)
cMono.chrX <- data.frame(split=1:10, cMono.chrX)

summary(cMono.chrX$MCC)

cMono.SLE <- list()
for(i in 1:10){
    cMono.SLE[[i]] <- read.csv(paste0('pseudobulk_update/split_', i, '/combined/ensemble/metrics_Classical_monocyte.SLE.csv'))
}
cMono.SLE <- do.call(rbind, cMono.SLE)
cMono.SLE <- data.frame(split=1:10, cMono.SLE)

# Report mean and standard deviation of MCC for both cell types
mean_cMono.chrX <- mean(cMono.chrX$MCC)
sd_cMono.chrX <- sd(cMono.chrX$MCC)
mean_cMono.SLE <- mean(cMono.SLE$MCC)
sd_cMono.SLE <- sd(cMono.SLE$MCC)
cat('Classical_monocyte.chrX - Mean MCC:', round(mean_cMono.chrX, 2), 'SD:', round(sd_cMono.chrX, 2), '\n')
cat('Classical_monocyte.SLE - Mean MCC:', round(mean_cMono.SLE, 2), 'SD:', round(sd_cMono.SLE, 2), '\n')


subset(features, celltype == 'Classical_monocyte.chrX')


# Repeat for Non_classical_monocyes
ncMono.chrX <- list()
for(i in 1:10){
    ncMono.chrX[[i]] <- read.csv(paste0('pseudobulk_update/split_', i, '/combined/ensemble/metrics_Non_classical_monocyte.chrX.csv'))
}
ncMono.chrX <- do.call(rbind, ncMono.chrX)
ncMono.chrX <- data.frame(split=1:10, ncMono.chrX)
ncMono.SLE <- list()
for(i in 1:10){
    ncMono.SLE[[i]] <- read.csv(paste0('pseudobulk_update/split_', i, '/combined/ensemble/metrics_Non_classical_monocyte.SLE.csv'))
}
ncMono.SLE <- do.call(rbind, ncMono.SLE)
ncMono.SLE <- data.frame(split=1:10, ncMono.SLE)
mean_ncMono.chrX <- mean(ncMono.chrX$MCC)
sd_ncMono.chrX <- sd(ncMono.chrX$MCC)
mean_ncMono.SLE <- mean(ncMono.SLE$MCC)
sd_ncMono.SLE <- sd(ncMono.SLE$MCC)
cat('Non_classical_monocyte.chrX - Mean MCC:', round(mean_ncMono.chrX, 2), 'SD:', round(sd_ncMono.chrX, 2), '\n')
cat('Non_classical_monocyte.SLE - Mean MCC:', round(mean_ncMono.SLE, 2), 'SD:', round(sd_ncMono.SLE, 2), '\n')
subset(features, celltype == 'Non_classical_monocyte.chrX') 


subset(edgeR[['CD8_positive,_alpha_beta_T_cell']], gene == 'TKTL1')