##***********************************##
## Allelome.PRO summarize replicates ## 
##***********************************##
## Tim Hasenbein
## Last modification 01.2023
## Creation: 06.2022
## Summarize the replicates from allelome.pro and take overlap
## allelic ratio is summarized by median, everything else with min
## Update 10.03.2023: changed the selection of min total reads to median total reads and also median for the reads of A1 and A2
library(tidyverse)
library(ggplot2)
library(data.table)


######------ Set environment -----###### 
input <- "/Users/timhasenbein/Desktop/mnt_3/ge26jec2/02_project_allelome/04_bodymap_9w/02_allelome_pro/02_strands_merged/"
output <- "/Users/timhasenbein/Desktop/mnt_2/andergassen_lab/03_project_allelome/03_bodymap_9w/02_allelome_pro/03_replicates_merged/"
samples <- c("Ao_9w","Br_9w","He_9w","Ki_9w","Li_9w","Lu_9w","Mu_9w","Sp_9w")
min_total_reads <- 20

sample <- "Br_9w"
######------  merge replicates  ------###### 
for (sample in samples){
  # get merged locus tables for each replicate
  s1 <- fread(paste0(input,sample,"_RNA_BxC_XX_1_merged_locus_table.txt"), header = T)
  s2 <- fread(paste0(input,sample,"_RNA_BxC_XX_2_merged_locus_table.txt"), header = T)
  s3 <- fread(paste0(input,sample,"_RNA_BxC_XX_3_merged_locus_table.txt"), header = T)
  # filter total reads
  s1 <- s1 %>% dplyr::filter(total_reads >= min_total_reads)
  s2 <- s2 %>% dplyr::filter(total_reads >= min_total_reads)
  s3 <- s3 %>% dplyr::filter(total_reads >= min_total_reads)
  # get overlapping genes and filter these from each locus table
  overlap <- s1 %>% merge(s2,by=c("chr","start","end","name")) %>% merge(s3,by=c("chr","start","end","name"))
  s1 <- s1 %>% dplyr::filter(name %in% overlap$name)
  s2 <- s2 %>% dplyr::filter(name %in% overlap$name)
  s3 <- s3 %>% dplyr::filter(name %in% overlap$name)
  # bind all replicates and select median AR, min total reads and min AS
  merge <- s1 %>% rbind(s2) %>% rbind(s3)
  locus_tab <- merge %>% group_by(name) %>% 
    dplyr::summarise(chr = unique(chr),start = unique(start), end = unique(end), name = unique(name), A1_reads = round(median(A1_reads),0), A2_reads = round(median(A2_reads),0), total_reads = round(median(total_reads),0), allelic_score = min(abs(allelic_score)), allelic_ratio = median(allelic_ratio))
  locus_tab <- locus_tab %>% dplyr::select("chr","start","end","name","A1_reads",
                                           "A2_reads","total_reads","allelic_score",
                                           "allelic_ratio")
  locus_tab <- locus_tab %>% arrange(chr,start)
  # output summarized locus table
  write.table(locus_tab,paste0(output,sample,"_locus_table_reps_update.txt"),col.names=T,row.names=F,sep="\t",quote=F)
  # make bedfile
  bed <- data.frame(V1=locus_tab$chr,V2=locus_tab$start,V3=locus_tab$end,V4=paste0(locus_tab$name,"_median_reads:",locus_tab$A1_reads,"/",locus_tab$A2_reads,"_min_score:",locus_tab$allelic_score,"_median_ratio:",locus_tab$allelic_ratio),V5=1000,V6=".",V7=locus_tab$start,V8=locus_tab$end,V9=locus_tab$allelic_ratio)
  # assign color
  bed$itemRgb[bed$V9 >= 0.65 & 
                bed$V9 < 0.8] <- "207,136,131"
  bed$itemRgb[bed$V9 >= 0.8 & 
                bed$V9 <= 1] <- "139,25,19"
  bed$itemRgb[bed$V9 <= 0.35 & 
                bed$V9 > 0.2] <- "140,164,202"
  bed$itemRgb[bed$V9 <= 0.2 & 
                bed$V9 >=0] <- "48,74,153"
  bed$itemRgb[bed$V9 > 0.35 & 
                bed$V9 < 0.65] <- "28,100,45"
  bed <- bed%>%dplyr::select(-V9)
  # output summarized bed file
  write.table(bed,paste0(output,sample,"_reps_colupdate.bed"),
              col.names=F,row.names=F,sep="\t",quote=F)
}
  
