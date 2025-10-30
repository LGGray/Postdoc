library(tidyverse)
library(data.table)

### ------ merge replicates for adult bodymap ------- ###

adult_metadata <- read.csv('adult_metadata.csv', sep=';')

sample_dirs <- list.dirs('adult', recursive=FALSE)
sample_dirs <- gsub('adult/', '', sample_dirs)
adult_metadata <- subset(adult_metadata, Sample_ID %in% sample_dirs)
adult_metadata$organ_age <- lapply(strsplit(adult_metadata$FID_comment, '_'), function(x) paste(x[1], x[2], sep="_"))

# No Predicted
lapply(unique(adult_metadata$organ_age), function(x) {
  samples <- subset(adult_metadata, organ_age == x)$Sample_ID
  rep1 <- fread(paste0('adult/', samples[1], '/annotation_no_predicted_merged_locus_table.txt'))
  rep2 <- fread(paste0('adult/', samples[2], '/annotation_no_predicted_merged_locus_table.txt'))
  rep3 <- fread(paste0('adult/', samples[3], '/annotation_no_predicted_merged_locus_table.txt'))

  # Filter for min total reads
  min_total_reads <- 20
  rep1 <- rep1 %>% filter(total_reads >= min_total_reads)
  rep2 <- rep2 %>% filter(total_reads >= min_total_reads)
  rep3 <- rep3 %>% filter(total_reads >= min_total_reads)

   # get overlapping genes and filter these from each locus table
  overlap <- rep1 %>% merge(rep2,by=c("chr","start","end","name")) %>% merge(rep3,by=c("chr","start","end","name"))
  rep1 <- rep1 %>% dplyr::filter(name %in% overlap$name)
  rep2 <- rep2 %>% dplyr::filter(name %in% overlap$name)
  rep3 <- rep3 %>% dplyr::filter(name %in% overlap$name)
  # bind all replicates and select median AR, min total reads and min AS
  merge <- rep1 %>% rbind(rep2) %>% rbind(rep3)
  locus_tab <- merge %>% group_by(name) %>% 
    dplyr::summarise(chr = unique(chr),start = unique(start), end = unique(end), name = unique(name), 
    A1_reads = round(median(A1_reads),0), A2_reads = round(median(A2_reads),0), 
    total_reads = round(median(total_reads),0), allelic_score = min(abs(allelic_score)), 
    allelic_ratio = median(allelic_ratio))
  locus_tab <- locus_tab %>% dplyr::select("chr","start","end","name","A1_reads",
                                           "A2_reads","total_reads","allelic_score",
                                           "allelic_ratio")
  locus_tab <- locus_tab %>% arrange(chr,start)
  # output summarized locus table
  if (!dir.exists(x)) {
  dir.create(x)
  }
  write.table(locus_tab, paste0(x, '/', 'annotation_no_predicted_locus_table_reps.txt'), row.names=FALSE, quote=FALSE, sep='\t')
})

# With predicted genes
lapply(unique(adult_metadata$organ_age), function(x) {
  samples <- subset(adult_metadata, organ_age == x)$Sample_ID
  rep1 <- fread(paste0('adult/', samples[1], '/annotation_merged_locus_table.txt'))
  rep2 <- fread(paste0('adult/', samples[2], '/annotation_merged_locus_table.txt'))
  rep3 <- fread(paste0('adult/', samples[3], '/annotation_merged_locus_table.txt'))

  # Filter for min total reads
  min_total_reads <- 20
  rep1 <- rep1 %>% filter(total_reads >= min_total_reads)
  rep2 <- rep2 %>% filter(total_reads >= min_total_reads)
  rep3 <- rep3 %>% filter(total_reads >= min_total_reads)

   # get overlapping genes and filter these from each locus table
  overlap <- rep1 %>% merge(rep2,by=c("chr","start","end","name")) %>% merge(rep3,by=c("chr","start","end","name"))
  rep1 <- rep1 %>% dplyr::filter(name %in% overlap$name)
  rep2 <- rep2 %>% dplyr::filter(name %in% overlap$name)
  rep3 <- rep3 %>% dplyr::filter(name %in% overlap$name)
  # bind all replicates and select median AR, min total reads and min AS
  merge <- rep1 %>% rbind(rep2) %>% rbind(rep3)
  locus_tab <- merge %>% group_by(name) %>% 
    dplyr::summarise(chr = unique(chr),start = unique(start), end = unique(end), name = unique(name), 
    A1_reads = round(median(A1_reads),0), A2_reads = round(median(A2_reads),0), 
    total_reads = round(median(total_reads),0), allelic_score = min(abs(allelic_score)), 
    allelic_ratio = median(allelic_ratio))
  locus_tab <- locus_tab %>% dplyr::select("chr","start","end","name","A1_reads",
                                           "A2_reads","total_reads","allelic_score",
                                           "allelic_ratio")
  locus_tab <- locus_tab %>% arrange(chr,start)
  # output summarized locus table
  if (!dir.exists(x)) {
  dir.create(x)
  }
  write.table(locus_tab, paste0(x, '/', 'annotation_locus_table_reps.txt'), row.names=FALSE, quote=FALSE, sep='\t')
})

### ------ merge replicates for aged bodymap ------- ###
aged_metadata <- read.csv('aged_metadata.csv', sep=',')

sample_dirs <- list.dirs('aged', recursive=FALSE)
sample_dirs <- gsub('aged/', '', sample_dirs)
aged_metadata <- subset(aged_metadata, Sample_ID %in% sample_dirs)
aged_metadata <- subset(aged_metadata, FID_comment != '')
aged_metadata$organ_age <- lapply(strsplit(aged_metadata$FID_comment, '_'), function(x) paste(x[1], x[2], sep="_"))

# No Predicted genes
lapply(unique(aged_metadata$organ_age), function(x) {
  samples <- subset(aged_metadata, organ_age == x)$Sample_ID
  rep1 <- fread(paste0('aged/', samples[1], '/annotation_no_predicted_merged_locus_table.txt'))
  rep2 <- fread(paste0('aged/', samples[2], '/annotation_no_predicted_merged_locus_table.txt'))
  rep3 <- fread(paste0('aged/', samples[3], '/annotation_no_predicted_merged_locus_table.txt'))

  # Filter for min total reads
  min_total_reads <- 20
  rep1 <- rep1 %>% filter(total_reads >= min_total_reads)
  rep2 <- rep2 %>% filter(total_reads >= min_total_reads)
  rep3 <- rep3 %>% filter(total_reads >= min_total_reads)

   # get overlapping genes and filter these from each locus table
  overlap <- rep1 %>% merge(rep2,by=c("chr","start","end","name")) %>% merge(rep3,by=c("chr","start","end","name"))
  rep1 <- rep1 %>% dplyr::filter(name %in% overlap$name)
  rep2 <- rep2 %>% dplyr::filter(name %in% overlap$name)
  rep3 <- rep3 %>% dplyr::filter(name %in% overlap$name)
  # bind all replicates and select median AR, min total reads and min AS
  merge <- rep1 %>% rbind(rep2) %>% rbind(rep3)
  locus_tab <- merge %>% group_by(name) %>% 
    dplyr::summarise(chr = unique(chr),start = unique(start), end = unique(end), name = unique(name), 
    A1_reads = round(median(A1_reads),0), A2_reads = round(median(A2_reads),0), 
    total_reads = round(median(total_reads),0), allelic_score = min(abs(allelic_score)), 
    allelic_ratio = median(allelic_ratio))
  locus_tab <- locus_tab %>% dplyr::select("chr","start","end","name","A1_reads",
                                           "A2_reads","total_reads","allelic_score",
                                           "allelic_ratio")
  locus_tab <- locus_tab %>% arrange(chr,start)
  # output summarized locus table
  if (!dir.exists(x)) {
  dir.create(x)
  }
  write.table(locus_tab, paste0(x, '/', 'annotation_no_predicted_locus_table_reps.txt'), row.names=FALSE, quote=FALSE, sep='\t')
})

# With predicted genes
lapply(unique(aged_metadata$organ_age), function(x) {
  samples <- subset(aged_metadata, organ_age == x)$Sample_ID
  rep1 <- fread(paste0('aged/', samples[1], '/annotation_merged_locus_table.txt'))
  rep2 <- fread(paste0('aged/', samples[2], '/annotation_merged_locus_table.txt'))
  rep3 <- fread(paste0('aged/', samples[3], '/annotation_merged_locus_table.txt'))

  # Filter for min total reads
  min_total_reads <- 20
  rep1 <- rep1 %>% filter(total_reads >= min_total_reads)
  rep2 <- rep2 %>% filter(total_reads >= min_total_reads)
  rep3 <- rep3 %>% filter(total_reads >= min_total_reads)

   # get overlapping genes and filter these from each locus table
  overlap <- rep1 %>% merge(rep2,by=c("chr","start","end","name")) %>% merge(rep3,by=c("chr","start","end","name"))
  rep1 <- rep1 %>% dplyr::filter(name %in% overlap$name)
  rep2 <- rep2 %>% dplyr::filter(name %in% overlap$name)
  rep3 <- rep3 %>% dplyr::filter(name %in% overlap$name)
  # bind all replicates and select median AR, min total reads and min AS
  merge <- rep1 %>% rbind(rep2) %>% rbind(rep3)
  locus_tab <- merge %>% group_by(name) %>% 
    dplyr::summarise(chr = unique(chr),start = unique(start), end = unique(end), name = unique(name), 
    A1_reads = round(median(A1_reads),0), A2_reads = round(median(A2_reads),0), 
    total_reads = round(median(total_reads),0), allelic_score = min(abs(allelic_score)), 
    allelic_ratio = median(allelic_ratio))
  locus_tab <- locus_tab %>% dplyr::select("chr","start","end","name","A1_reads",
                                           "A2_reads","total_reads","allelic_score",
                                           "allelic_ratio")
  locus_tab <- locus_tab %>% arrange(chr,start)
  # output summarized locus table
  if (!dir.exists(x)) {
  dir.create(x)
  }
  write.table(locus_tab, paste0(x, '/', 'annotation_locus_table_reps.txt'), row.names=FALSE, quote=FALSE, sep='\t')
})

##################################
# Merge replicates for TAC model #
##################################

tac_metadata <- read.csv('metadata.csv', sep=',')
tac_metadata <- subset(tac_metadata, Sample_ID %in% gsub('\\.\\/', '', list.dirs(recursive=FALSE)))
tac_metadata$condition <- gsub('_\\d', '', tac_metadata$FID_comment)

# No Predicted genes
lapply(unique(tac_metadata$condition), function(x) {
  samples <- subset(tac_metadata, condition == x)$Sample_ID
  rep1 <- fread(paste0(samples[1], '/annotation_no_predicted_merged_locus_table.txt'))
  rep2 <- fread(paste0(samples[2], '/annotation_no_predicted_merged_locus_table.txt'))
  rep3 <- fread(paste0(samples[3], '/annotation_no_predicted_merged_locus_table.txt'))

  # Filter for min total reads
  min_total_reads <- 20
  rep1 <- rep1 %>% filter(total_reads >= min_total_reads)
  rep2 <- rep2 %>% filter(total_reads >= min_total_reads)
  rep3 <- rep3 %>% filter(total_reads >= min_total_reads)

  # get overlapping genes and filter these from each locus table
  overlap <- rep1 %>% merge(rep2,by=c("chr","start","end","name")) %>% merge(rep3,by=c("chr","start","end","name"))
  rep1 <- rep1 %>% dplyr::filter(name %in% overlap$name)
  rep2 <- rep2 %>% dplyr::filter(name %in% overlap$name)
  rep3 <- rep3 %>% dplyr::filter(name %in% overlap$name)
  # bind all replicates and select median AR, min total reads and min AS
  merge <- rep1 %>% rbind(rep2) %>% rbind(rep3)
  locus_tab <- merge %>% group_by(name) %>% 
    dplyr::summarise(chr = unique(chr),start = unique(start), end = unique(end), name = unique(name), 
    A1_reads = round(median(A1_reads),0), A2_reads = round(median(A2_reads),0), 
    total_reads = round(median(total_reads),0), allelic_score = min(abs(allelic_score)), 
    allelic_ratio = median(allelic_ratio))
  locus_tab <- locus_tab %>% dplyr::select("chr","start","end","name","A1_reads",
                                           "A2_reads","total_reads","allelic_score",
                                           "allelic_ratio")
  locus_tab <- locus_tab %>% arrange(chr,start)
  # output summarized locus table
  if (!dir.exists(x)) {
  dir.create(x)
  }
  write.table(locus_tab, paste0(x, '/', 'annotation_no_predicted_locus_table_reps.txt'), row.names=FALSE, quote=FALSE, sep='\t')
})

# With predicted genes
lapply(unique(tac_metadata$condition), function(x) {
  samples <- subset(tac_metadata, condition == x)$Sample_ID
  rep1 <- fread(paste0(samples[1], '/annotation_merged_locus_table.txt'))
  rep2 <- fread(paste0(samples[2], '/annotation_merged_locus_table.txt'))
  rep3 <- fread(paste0(samples[3], '/annotation_merged_locus_table.txt'))

  # Filter for min total reads
  min_total_reads <- 20
  rep1 <- rep1 %>% filter(total_reads >= min_total_reads)
  rep2 <- rep2 %>% filter(total_reads >= min_total_reads)
  rep3 <- rep3 %>% filter(total_reads >= min_total_reads)

  # get overlapping genes and filter these from each locus table
  overlap <- rep1 %>% merge(rep2,by=c("chr","start","end","name")) %>% merge(rep3,by=c("chr","start","end","name"))
  rep1 <- rep1 %>% dplyr::filter(name %in% overlap$name)
  rep2 <- rep2 %>% dplyr::filter(name %in% overlap$name)
  rep3 <- rep3 %>% dplyr::filter(name %in% overlap$name)
  # bind all replicates and select median AR, min total reads and min AS
  merge <- rep1 %>% rbind(rep2) %>% rbind(rep3)
  locus_tab <- merge %>% group_by(name) %>% 
    dplyr::summarise(chr = unique(chr),start = unique(start), end = unique(end), name = unique(name), 
    A1_reads = round(median(A1_reads),0), A2_reads = round(median(A2_reads),0), 
    total_reads = round(median(total_reads),0), allelic_score = min(abs(allelic_score)), 
    allelic_ratio = median(allelic_ratio))
  locus_tab <- locus_tab %>% dplyr::select("chr","start","end","name","A1_reads",
                                           "A2_reads","total_reads","allelic_score",
                                           "allelic_ratio")
  locus_tab <- locus_tab %>% arrange(chr,start)
  # output summarized locus table
  if (!dir.exists(x)) {
  dir.create(x)
  }
  write.table(locus_tab, paste0(x, '/', 'annotation_locus_table_reps.txt'), row.names=FALSE, quote=FALSE, sep='\t')
})