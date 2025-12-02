library(tidyverse)
library(data.table)

### ------ merge replicates for adult bodymap ------- ###

adult_metadata <- read.csv('adult_metadata.csv', sep=',')

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

## Merge replicates for cardiac celltypes ##
metadata <- read.csv('SraRunTable.csv')
metadata$sex_tissue <- paste0(metadata$sex, '_', metadata$tissue)

date <- '2025_11_20'
lapply(unique(metadata$sex_tissue), function(x){
  samples <- subset(metadata, sex_tissue == x)$Run
  rep1 <- fread(paste0(samples[1], '/', date, '_', samples[1], '_Aligned.sortedByCoord.out_annotation_us.bed_1/locus_table.txt'))
  rep2 <- fread(paste0(samples[2], '/', date, '_', samples[2], '_Aligned.sortedByCoord.out_annotation_us.bed_1/locus_table.txt'))
  rep3 <- fread(paste0(samples[3], '/', date, '_', samples[3], '_Aligned.sortedByCoord.out_annotation_us.bed_1/locus_table.txt'))

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
  dir.create(gsub(" ", "_", x))
  }
  write.table(locus_tab, paste0(gsub(" ", "_", x), '/', 'annotation_us_locus_table_reps.txt'), row.names=FALSE, quote=FALSE, sep='\t')
})

# Merge replicates for H3K4me3 chipseq
rep_1 <- read.delim('2025_11_25_rep_1.clean_H3K4me3_consensus_peaks.bed_1/locus_table.txt')
rep_2 <- read.delim('2025_11_25_rep_2.clean_H3K4me3_consensus_peaks.bed_1/locus_table.txt')

# # # Filter for min total reads
# min_total_reads <- 20
# rep_1 <- rep_1 %>% filter(total_reads >= min_total_reads)
# rep_2 <- rep_2 %>% filter(total_reads >= min_total_reads)

# # Filter for allelic ratio >= 0.7 and <= 0.3
# rep_1 <- rep_1 %>% filter(allelic_ratio >= 0.65 | allelic_ratio <= 0.35)
# rep_2 <- rep_2 %>% filter(allelic_ratio >= 0.65 | allelic_ratio <= 0.35)

# get overlapping genes and filter these from each locus table
overlap <- rep_1 %>% merge(rep_2,by=c("chr","start","end","name"))
rep_1 <- rep_1 %>% dplyr::filter(name %in% overlap$name)
rep_2 <- rep_2 %>% dplyr::filter(name %in% overlap$name)
# bind all replicates and select median AR, min total reads and min AS
merge <- rep_1 %>% rbind(rep_2)
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
write.table(locus_tab, paste0('H3K4me3_chipseq_annotation_locus_table_reps.txt'), row.names=FALSE, quote=FALSE, sep='\t')

# Write BED file
ap <- locus_tab
bed <- ap[c(1:3,7,9)]
bed$name <- paste0(ap$name,"_reads:",ap$A1_reads,"/",ap$A2_reads,"_score:",ap$allelic_score,"_ratio:",ap$allelic_ratio)
bed$score <- "1000"
bed$strand <- "."
bed$thickStart <-  paste0(bed$start)
bed$thickEnd <- paste0(bed$end)
bed$itemRgb[bed$allelic_ratio >= 0.6 & 
            bed$allelic_ratio < 0.8] <- "207,136,131"
bed$itemRgb[bed$allelic_ratio >= 0.8 & 
            bed$allelic_ratio <= 1] <- "139,25,19"
bed$itemRgb[bed$allelic_ratio <= 0.4 & 
            bed$allelic_ratio > 0.2] <- "140,164,202"
bed$itemRgb[bed$allelic_ratio <= 0.2 & 
            bed$allelic_ratio >=0] <- "48,74,153"
bed$itemRgb[bed$allelic_ratio > 0.4 & 
            bed$allelic_ratio < 0.6] <- "28,100,45"
bed$itemRgb[bed$total_reads < 20] <- "80,80,80"
bed <- bed[,c(1:3,6:11)]
colnames(bed) <- c("track",paste0("name= H3K4me3"),paste0("description= H3K4me3"),"visibility=2","itemRgb=On","useScore=1","","","")
write.table(bed, 'H3K4me3_ENCODE.bed', quote=FALSE, sep = "\t", col.names = FALSE, row.names =FALSE)

# Merge replicates for H3K27ac chipseq
rep_1 <- read.delim('2025_11_28_rep_1.clean_H3K27ac_consensus_peaks.bed_1/locus_table.txt')
rep_2 <- read.delim('2025_11_28_rep_2.clean_H3K27ac_consensus_peaks.bed_1/locus_table.txt')

# # # Filter for min total reads
# min_total_reads <- 20
# rep_1 <- rep_1 %>% filter(total_reads >= min_total_reads)
# rep_2 <- rep_2 %>% filter(total_reads >= min_total_reads)

# # Filter for allelic ratio >= 0.7 and <= 0.3
# rep_1 <- rep_1 %>% filter(allelic_ratio >= 0.65 | allelic_ratio <= 0.35)
# rep_2 <- rep_2 %>% filter(allelic_ratio >= 0.65 | allelic_ratio <= 0.35)

# get overlapping genes and filter these from each locus table
overlap <- rep_1 %>% merge(rep_2,by=c("chr","start","end","name"))
rep_1 <- rep_1 %>% dplyr::filter(name %in% overlap$name)
rep_2 <- rep_2 %>% dplyr::filter(name %in% overlap$name)
# bind all replicates and select median AR, min total reads and min AS
merge <- rep_1 %>% rbind(rep_2)
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
write.table(locus_tab, paste0('H3K27ac_chipseq_annotation_locus_table_reps.txt'), row.names=FALSE, quote=FALSE, sep='\t')

# Write BED file
ap <- locus_tab
bed <- ap[c(1:3,7,9)]
bed$name <- paste0(ap$name,"_reads:",ap$A1_reads,"/",ap$A2_reads,"_score:",ap$allelic_score,"_ratio:",ap$allelic_ratio)
bed$score <- "1000"
bed$strand <- "."
bed$thickStart <-  paste0(bed$start)
bed$thickEnd <- paste0(bed$end)
bed$itemRgb[bed$allelic_ratio >= 0.6 & 
            bed$allelic_ratio < 0.8] <- "207,136,131"
bed$itemRgb[bed$allelic_ratio >= 0.8 &  
            bed$allelic_ratio <= 1] <- "139,25,19"
bed$itemRgb[bed$allelic_ratio <= 0.4 & 
            bed$allelic_ratio > 0.2] <- "140,164,202"
bed$itemRgb[bed$allelic_ratio <= 0.2 & 
            bed$allelic_ratio >=0] <- "48,74,153"
bed$itemRgb[bed$allelic_ratio > 0.4 & 
            bed$allelic_ratio < 0.6] <- "28,100,45"
bed$itemRgb[bed$total_reads < 20] <- "80,80,80"
bed <- bed[,c(1:3,6:11)]
colnames(bed) <- c("track",paste0("name= H3K27ac"),paste0("description= H3K27ac"),"visibility=2","itemRgb=On","useScore=1","","","")
write.table(bed, 'H3K27ac_ENCODE.bed', quote=FALSE, sep = "\t", col.names = FALSE, row.names =FALSE)