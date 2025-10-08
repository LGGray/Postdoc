##*******************************##
##  Allelome.PRO  merge strands  ##
##*******************************##
## Project: Allelome.Link
## Tim Hasenbein
## Last modification 01.2023
## Creation: 12.2022
## Merge fwd and reverse strand
library(tidyverse)
library(data.table)


######------  set environment  ------###### 
input <- "/Users/timhasenbein/Desktop/mnt_3/ge26jec2/02_project_allelome/04_bodymap_9w/02_allelome_pro/01_stranded_APruns/"
output <- "/Users/timhasenbein/Desktop/mnt_3/ge26jec2/02_project_allelome/04_bodymap_9w/02_allelome_pro/"
samples <- c("2023_01_04_Ao_9w_RNA_BxC_XX_1", "2023_01_04_Ao_9w_RNA_BxC_XX_2", "2023_01_04_Ao_9w_RNA_BxC_XX_3","2023_01_04_Br_9w_RNA_BxC_XX_1", "2023_01_04_Br_9w_RNA_BxC_XX_2", "2023_01_04_Br_9w_RNA_BxC_XX_3","2023_01_04_He_9w_RNA_BxC_XX_1", "2023_01_04_He_9w_RNA_BxC_XX_2", "2023_01_04_He_9w_RNA_BxC_XX_3","2023_01_04_Ki_9w_RNA_BxC_XX_1", "2023_01_04_Ki_9w_RNA_BxC_XX_2", "2023_01_04_Ki_9w_RNA_BxC_XX_3","2023_01_04_Li_9w_RNA_BxC_XX_1", "2023_01_04_Li_9w_RNA_BxC_XX_2", "2023_01_04_Li_9w_RNA_BxC_XX_3","2023_01_04_Lu_9w_RNA_BxC_XX_1", "2023_01_04_Lu_9w_RNA_BxC_XX_2", "2023_01_04_Lu_9w_RNA_BxC_XX_3","2023_01_04_Mu_9w_RNA_BxC_XX_1", "2023_01_04_Mu_9w_RNA_BxC_XX_2", "2023_01_04_Mu_9w_RNA_BxC_XX_3","2023_01_04_Sp_9w_RNA_BxC_XX_1", "2023_01_04_Sp_9w_RNA_BxC_XX_2", "2023_01_04_Sp_9w_RNA_BxC_XX_3")


######------  merge strands  ------###### 
for (sample in samples){
  sample_name <- gsub("2023_01_04_","",sample)
  # merge locus table
  fwd_locus <- fread(paste0(input,sample,"_fwd_annotation_f.bed_1/locus_table.txt"), header = T)
  rev_locus <- fread(paste0(input,sample,"_rev_annotation_r.bed_1/locus_table.txt"), header = T)
  merge_locus <- fwd_locus%>%rbind(rev_locus)  
  merge_locus <- merge_locus%>%arrange(chr,start)
  # merge bed files
  fwd_bed <- fread(paste0(input,sample,"_fwd_annotation_f.bed_1/BED_files/",sample_name,"_fwd.bed"), header = F)
  rev_bed <- fread(paste0(input,sample,"_rev_annotation_r.bed_1/BED_files/",sample_name,"_rev.bed"), header = F)
  merge_bed <- fwd_bed%>%rbind(rev_bed)  
  merge_bed <- merge_bed%>%arrange(V1,V2)
 # write output
  write.table(merge_locus,paste0(output,sample_name,"_merged_locus_table.txt"),col.names=T, row.names=F,sep="\t",quote=F)
  write.table(merge_bed,paste0(output,sample_name,"_merged.bed"),col.names=F, row.names=F,sep="\t",quote=F)
}
  