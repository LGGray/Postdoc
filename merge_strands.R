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

forward <- commandArgs(trailingOnly = TRUE)[1]
reverse <- commandArgs(trailingOnly = TRUE)[2]
output <- commandArgs(trailingOnly = TRUE)[3]

# Merge locus table
fwd_locus <- fread(forward, header = T)
rev_locus <- fread(reverse, header = T)
merge_locus <- fwd_locus%>%rbind(rev_locus)  
merge_locus <- merge_locus%>%arrange(chr,start)

write.table(merge_locus, paste0(output, '/merged_locus_table.txt'), col.names=T, row.names=F, sep="\t", quote=F)
