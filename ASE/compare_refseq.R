library(tidyverse)
library(data.table)

# Read in aged bodymap data
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
names(no_pred_65) <- bodymap_dirs

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

# Read in TAC data
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


