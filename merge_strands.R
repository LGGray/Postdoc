##*******************************##
##  Allelome.PRO  merge strands  ##
##*******************************##

library(data.table)

sample <- commandArgs(trailingOnly = TRUE)[1]

forward <- list.files(paste0('adult/', sample), pattern='locus_table.txt', recursive=T, full.names = TRUE)[1]
reverse <- list.files(paste0('adult/', sample), pattern='locus_table.txt', recursive=T, full.names = TRUE)[2]


# Merge locus table
fwd_locus <- fread(forward, header = T)
rev_locus <- fread(reverse, header = T)
merge_locus <- rbind(fwd_locus, rev_locus)
merge_locus <- merge_locus[order(merge_locus$chr, merge_locus$start), ]

write.table(merge_locus, paste0('adult/', sample, '/merged_locus_table.txt'), col.names=T, row.names=F, sep="\t", quote=F)
