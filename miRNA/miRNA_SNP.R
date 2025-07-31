### Code to overlap CAST or BL6 SNPs with miRNA primary transcripts on chrX ###

# Read in miRNA coordinates
miRNA_gff <- read.delim('miRNA/mmu.gff3', comment.char = "#", header = FALSE)
miRNA_gff <- subset(miRNA_gff, V3 == 'miRNA_primary_transcript' & V1 == 'chrX')

# Read in BL6 SNP
BL6 <- read.delim('mouse_SNP/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz', comment.char = "#", header = FALSE)
BL6_X <- subset(BL6, V1 == 'X')

hits_BL6 <- lapply(BL6_X$V2, function(x){
    # Find the miRNA that overlaps with the SNP position
    subset(miRNA, V4 <= x & V5 >= x)
})
hits_BL6[sapply(hits, nrow) > 0]


# Read in CAST SNPs
CAST <- read.delim('mouse_SNP/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz', comment.char = "#", header = FALSE)
CAST_X <- subset(CAST, V1 == 'X')

hits_cast <- lapply(CAST_X$V2, function(x){
    # Find the miRNA that overlaps with the SNP position
    subset(miRNA, V4 <= x & V5 >= x)
})
names(hits_cast) <- CAST_X$V2

result <- hits_cast[sapply(hits_cast, nrow) > 0]

save(result, file = 'miRNA/CAST_miRNA_SNPs.RData')

