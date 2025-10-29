library(data.table)
library(dplyr)


# #### Match RefSeq transcript ID to gene name
# gtf <- read.delim('/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/GCF_000001635.27_GRCm39_genomic.gtf', comment.char = "#", header = FALSE)
# gtf_transcripts <- subset(gtf, V3 == "transcript")

# attrs <- strsplit(as.character(gtf_transcripts$V9), "; ", fixed = TRUE)
# gene_ids <- vapply(attrs, function(x) {
#   if (length(x) >= 1) sub("^gene_id\\s*", "", x[1]) else NA_character_
# }, character(1))
# transcript_ids <- vapply(attrs, function(x) {
#   ti <- grep("^transcript_id", x, value = TRUE)
#   if (length(ti)) sub("^transcript_id\\s*", "", ti[1]) else NA_character_
# }, character(1))
# biotype_ids <- vapply(attrs, function(x) {
#   bt <- grep("^transcript_biotype", x, value = TRUE)
#   if (length(bt)) sub("^transcript_biotype\\s*", "", bt[1]) else NA_character_
# }, character(1))

# # make gtf data.tables
# dt_gtf <- data.table(gtf_transcripts)[, .(transcript_id = transcript_ids, gene_ids, biotype_ids)]

# Read in links table
link_file <- list.files(pattern = "links_full_table.txt", full.names = TRUE, recursive = TRUE)
link_file <- grep("2025-10-27", link_file, value=TRUE)
links_table <- fread(link_file, data.table=FALSE)

# Add gene name and biotype
# dt_x <- data.table(links_table)
# dt_x[, gene_base := dt_gtf[.SD, on = .(transcript_id = name_base), gene_ids]]
# dt_x[, biotype_base := dt_gtf[.SD, on = .(transcript_id = name_base), biotype_ids]]
# dt_x[, gene_target := dt_gtf[.SD, on = .(transcript_id = name_target), gene_ids]]
# dt_x[, biotype_target := dt_gtf[.SD, on = .(transcript_id = name_target), biotype_ids]]
# links_table <- as.data.frame(dt_x)


bedpe <- links_table %>%
  mutate(
    name = paste(name_base, name_target, mechanism, sep = "_"),
    score = ifelse(mechanism == "enhancing", 1, -1),
    strand1 = ".",
    strand2 = ".",
    color = ifelse(mechanism == "enhancing", "0,128,0", "188,34,13")
  ) %>%
  select(chr_base, start_base, end_base, chr_target, start_target, end_target,
         name, score, strand1, strand2, color)

write.table(bedpe, "links.bedpe", sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE)

# Create IGV-compatible bedpe with header, then append the bedpe contents
header <- 'track type=interact name="lncRNA-pcGene_links" description="Enhancing and Repressive" useScore=1'
writeLines(header, con = "links_IGV.bedpe")
if (!file.exists("links.bedpe")) stop("links.bedpe not found")
cat(readLines("links.bedpe"), file = "links_IGV.bedpe", sep = "\n", append = TRUE)


locus_tab <- fread('locus_table_reps.txt', data.table=FALSE)
# make bedfile
bed <- data.frame(V1=locus_tab$chr,V2=locus_tab$start,V3=locus_tab$end,V4=paste0(locus_tab$name,"_median_reads:",locus_tab$A1_reads,"/",locus_tab$A2_reads,"_min_score:",round(locus_tab$allelic_score, 2),"_median_ratio:",locus_tab$allelic_ratio),V5=1000,V6=".",V7=locus_tab$start,V8=locus_tab$end,V9=locus_tab$allelic_ratio)
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

write.table(bed, 'locus_table_reps.bed', row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
