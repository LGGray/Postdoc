library(tidyverse)
library(data.table)

#############################
# Read in aged bodymap data #
#############################
bodymap_dirs <- list.dirs('adult_aged_bodymap', recursive = FALSE, full.names = FALSE)
bodymap_dirs <- grep('_', bodymap_dirs, value = TRUE)

# No pred 70:30
bodymap_no_pred_70 <- lapply(bodymap_dirs, function(dir){
    tmp <- fread(list.files(paste0('adult_aged_bodymap/', dir), pattern = 'no_predicted_70_links_full_table.txt', full.names = TRUE, recursive = TRUE), data.table = FALSE)
    return(tmp)
})
names(bodymap_no_pred_70) <- bodymap_dirs

# No pred 65:35
bodymap_no_pred_65 <- lapply(bodymap_dirs, function(dir){
    tmp <- fread(list.files(paste0('adult_aged_bodymap/', dir), pattern = 'no_predicted_65_links_full_table.txt', full.names = TRUE, recursive = TRUE), data.table = FALSE)
    return(tmp)
})
names(bodymap_no_pred_65) <- bodymap_dirs

# Pred 70:30
bodymap_pred_70 <- lapply(bodymap_dirs, function(dir){
    tmp <- fread(list.files(paste0('adult_aged_bodymap/', dir), pattern = '^predicted_70_links_full_table.txt', full.names = TRUE, recursive = TRUE), data.table = FALSE)
    return(tmp)
})
names(bodymap_pred_70) <- bodymap_dirs

# Pred 65:35
bodymap_pred_65 <- lapply(bodymap_dirs, function(dir){
    tmp <- fread(list.files(paste0('adult_aged_bodymap/', dir), pattern = '^predicted_65_links_full_table.txt', full.names = TRUE, recursive = TRUE), data.table = FALSE)
    return(tmp)
})
names(bodymap_pred_65) <- bodymap_dirs

####################
# Read in TAC data #
####################
TAC_dirs <- list.dirs('F1_TAC_Sarah', recursive = FALSE, full.names = FALSE)
TAC_dirs <- grep('He', TAC_dirs, value = TRUE)

# No pred 70:30
TAC_no_pred_70 <- lapply(TAC_dirs, function(dir){
    tmp <- fread(list.files(paste0('F1_TAC_Sarah/', dir), pattern = 'no_predicted_70_links_full_table.txt', full.names = TRUE, recursive = TRUE), data.table = FALSE)
    return(tmp)
})
names(TAC_no_pred_70) <- TAC_dirs

# No pred 65:35
TAC_no_pred_65 <- lapply(TAC_dirs, function(dir){
    tmp <- fread(list.files(paste0('F1_TAC_Sarah/', dir), pattern = 'no_predicted_65_links_full_table.txt', full.names = TRUE, recursive = TRUE), data.table = FALSE)
    return(tmp)
})
names(TAC_no_pred_65) <- TAC_dirs

# Pred 70:30
TAC_pred_70 <- lapply(TAC_dirs, function(dir){
    tmp <- fread(list.files(paste0('F1_TAC_Sarah/', dir), pattern = '^predicted_70_links_full_table.txt', full.names = TRUE, recursive = TRUE), data.table = FALSE)
    return(tmp)
})
names(TAC_pred_70) <- TAC_dirs

# Pred 65:35
TAC_pred_65 <- lapply(TAC_dirs, function(dir){
    tmp <- fread(list.files(paste0('F1_TAC_Sarah/', dir), pattern = '^predicted_65_links_full_table.txt', full.names = TRUE, recursive = TRUE), data.table = FALSE)
    return(tmp)
})
names(TAC_pred_65) <- TAC_dirs


# Subset for heart
no_pred_70_heart_9w <- bodymap_no_pred_70[['He_9w']]
no_pred_65_heart_9w <- bodymap_no_pred_65[['He_9w']]
pred_70_heart_9w <- bodymap_pred_70[['He_9w']]
pred_65_heart_9w <- bodymap_pred_65[['He_9w']]
no_pred_70_heart_78w <- bodymap_no_pred_70[['He_78w']]
no_pred_65_heart_78w <- bodymap_no_pred_65[['He_78w']]
pred_70_heart_78w <- bodymap_pred_70[['He_78w']]
pred_65_heart_78w <- bodymap_pred_65[['He_78w']]
no_pred_70_TAC <- TAC_no_pred_70[['He_TAC28d_RNA_XistBxC_-+_XX']]
no_pred_65_TAC <- TAC_no_pred_65[['He_TAC28d_RNA_XistBxC_-+_XX']]
pred_70_TAC <- TAC_pred_70[['He_TAC28d_RNA_XistBxC_-+_XX']]
pred_65_TAC <- TAC_pred_65[['He_TAC28d_RNA_XistBxC_-+_XX']]

# Plot the nrow of each dataframe grouped by reference and ratio
plot_data <- data.frame(
    Sample = c('Heart_9w', 'Heart_9w', 'Heart_9w', 'Heart_9w',
               'Heart_78w', 'Heart_78w', 'Heart_78w', 'Heart_78w',
               'TAC', 'TAC', 'TAC', 'TAC'),
    reference = c('NoPred_70', 'NoPred_65', 'Pred_70', 'Pred_65',
                  'NoPred_70', 'NoPred_65', 'Pred_70', 'Pred_65',
                  'NoPred_70', 'NoPred_65', 'Pred_70', 'Pred_65'),
    Num_Links = c(nrow(no_pred_70_heart_9w), nrow(no_pred_65_heart_9w), nrow(pred_70_heart_9w), nrow(pred_65_heart_9w),
                  nrow(no_pred_70_heart_78w), nrow(no_pred_65_heart_78w), nrow(pred_70_heart_78w), nrow(pred_65_heart_78w),
                  nrow(no_pred_70_TAC), nrow(no_pred_65_TAC), nrow(pred_70_TAC), nrow(pred_65_TAC))
)
plot_data$reference <- factor(plot_data$reference, levels = c('NoPred_70', 'NoPred_65', 'Pred_70', 'Pred_65'))

pdf('LINKS_study/compare_refseq.pdf')
ggplot(plot_data, aes(x = reference, y = Num_Links, fill = Sample)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = 'Number of ASE Links by Reference and Ratio', x = 'Sample', y = 'Number of ASE Links')
dev.off()


# Read in snRNAseq data

### Adult files ###
adult.files <- list.files('adult_aged_heart_snRNAseq/9w/Allelome.LINK', pattern = 'links_full_table.txt', full.names = TRUE, recursive = TRUE)
adult.files <- grep('2025_10_30', adult.files, value=TRUE)

adult_files_nopredict <- grep('no_predicted', adult.files, value=T)
adult_predict_link <- lapply(adult_files_predict, function(f) {
    tmp <- read.delim(f)
    return(tmp)
})
names(adult_predict_link) <- lapply(adult_files_nopredict, function(f) {
    ct <- strsplit(f, '/')[[1]][5]
    return(ct)
})
adult_nopredict_df <- data.frame(
    Celltype = gsub('_.+', '', names(adult_predict_link)),
    reference_ratio = sub('^[^_]+_', '', names(adult_predict_link)),
    Num_Links = sapply(adult_predict_link, nrow)
)


adult_files_predict <- grep('no_predicted', adult.files, value=T, invert=TRUE)
adult_predict_link <- lapply(adult_files_predict, function(f) {
    tmp <- read.delim(f)
    return(tmp)
})
names(adult_predict_link) <- lapply(adult_files_predict, function(f) {
    ct <- strsplit(f, '/')[[1]][5]
    return(ct)
})
adult_predict_df <- data.frame(
    Celltype = gsub('_.+', '', names(adult_predict_link)),
    reference_ratio = sub('^[^_]+_', '', names(adult_predict_link)),
    Num_Links = sapply(adult_predict_link, nrow)
)
adult_predict_df$reference_ratio <- paste0('predicted_', adult_predict_df$reference_ratio)

### Aged files ###
aged.files <- list.files('adult_aged_heart_snRNAseq/78w/Allelome.LINK', pattern = 'links_full_table.txt', full.names = TRUE, recursive = TRUE)
aged.files <- grep('2025_10_30', aged.files, value=TRUE)

aged_files_nopredict <- grep('no_predicted', aged.files, value=T)
aged_nopredict_link <- lapply(aged_files_nopredict, function(f) {
    tmp <- read.delim(f)
    return(tmp)
})
names(aged_nopredict_link) <- lapply(aged_files_nopredict, function(f) {
    ct <- strsplit(f, '/')[[1]][5]
    return(ct)
})
aged_nopredict_df <- data.frame(
    Celltype = gsub('_.+', '', names(aged_nopredict_link)),
    reference_ratio = sub('^[^_]+_', '', names(aged_nopredict_link)),
    Num_Links = sapply(aged_nopredict_link, nrow)
)


aged_files_predict <- grep('predicted', aged.files, value=T, invert=TRUE)
aged_predict_link <- lapply(aged_files_predict, function(f) {
    tmp <- read.delim(f)
    return(tmp)
})
names(aged_predict_link) <- lapply(aged_files_predict, function(f) {
    ct <- strsplit(f, '/')[[1]][5]
    return(ct)
})
aged_predict_df <- data.frame(
    Celltype = gsub('_.+', '', names(aged_predict_link)),
    reference_ratio = sub('^[^_]+_', '', names(aged_predict_link)),
    Num_Links = sapply(aged_predict_link, nrow)
)
aged_predict_df$reference_ratio <- paste0('predicted_', aged_predict_df$reference_ratio)

# Combine adult and aged dataframes
combined_dataframe <- rbind(adult_nopredict_df, adult_predict_df, aged_nopredict_df, aged_predict_df)
combined_dataframe$reference_ratio <- factor(combined_dataframe$reference_ratio, levels = c('no_predicted_70', 'no_predicted_65', 'predicted_70', 'predicted_65'))

pdf('LINKS_study/compare_refseq_snRNAseq.pdf')
ggplot(combined_dataframe, aes(x = reference_ratio, y = Num_Links, fill = Celltype)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = 'Number of ASE Links by Reference and Ratio in snRNAseq', x = 'o', y = 'Number of ASE Links')
dev.off()   