library(dplyr)
library(tidyr)
library(DESeq2)


gencode <- read.delim('/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/gencode.vM37.primary_assembly.annotation.gtf', header=FALSE, comment.char='#')
gencode <- subset(gencode, V3 == 'gene')
gene_name <- lapply(gencode$V9, function(x) {
    ENSG <- gsub('gene_id ', '', strsplit(x, '; ')[[1]][1])
    gene <- gsub('gene_name ', '', strsplit(x, '; ')[[1]][3])
    data.frame(ENSG=ENSG, gene=gene)
})
gene_name <- bind_rows(gene_name)

########################
# Read in bodymap data #
########################
samples <- c(list.dirs('adult', recursive=FALSE), list.dirs('aged', recursive=FALSE))
samples <- gsub('^(adult|aged)/', '', samples)
adult_metadata <- read.delim('adult_metadata.csv', sep=';')
aged_metadata  <- read.delim('aged_metadata.csv', sep=',')
metadata <- rbind(adult_metadata, aged_metadata)
metadata <- metadata[metadata$Sample_ID %in% samples & metadata$FID_comment != "", ]

metadata$tissue <- unlist(lapply(metadata$FID_comment, function(x) strsplit(x, '_')[[1]][1]))
metadata$age  <- lapply(metadata$FID_comment, function(x) strsplit(x, '_')[[1]][2])
metadata$tissue_age <- paste0(metadata$tissue, '_', metadata$age)

tissue_GEM <- lapply(unique(metadata$tissue), function(x) {
    tmp <- subset(metadata, tissue == x)

    adult_rep1 <- read.delim(paste0('adult/', tmp$Sample_ID[1], '/', tmp$Sample_Name[1], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[1], '_rep1')))
    adult_rep2 <- read.delim(paste0('adult/', tmp$Sample_ID[2], '/', tmp$Sample_Name[2], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[2], '_rep2')))
    adult_rep3 <- read.delim(paste0('adult/', tmp$Sample_ID[3], '/', tmp$Sample_Name[3], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[3], '_rep3')))
    aged_rep1  <- read.delim(paste0('aged/', tmp$Sample_ID[4], '/', tmp$Sample_Name[4], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[4], '_rep1')))
    aged_rep2  <- read.delim(paste0('aged/', tmp$Sample_ID[5], '/', tmp$Sample_Name[5], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[5], '_rep2')))
    aged_rep3  <- read.delim(paste0('aged/', tmp$Sample_ID[6], '/', tmp$Sample_Name[6], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[6], '_rep3')))
    merged <- Reduce(function(x, y) merge(x, y, by='geneID', all=FALSE), list(adult_rep1, adult_rep2, adult_rep3, aged_rep1, aged_rep2, aged_rep3))
    merged <- merged[grep('^ENS', merged$geneID), ]
    # # add gene names
    # merged <- merge(merged, gene_name, by.x='geneID', by.y='ENSG', all.x=TRUE)
    rownames(merged) <- merged$geneID
    merged <- merged[, 2:(ncol(merged))]
    return(merged)
})
save(tissue_GEM, file='bodymap_tissue_GEM.Rdata')

write.table(merged, file=paste0('DEG/', x, '_GEM.txt'), sep='\t', quote=FALSE)

# Set Up BiomaRt
library(biomaRt)
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

###########################################
# DESeq2 differential expression analysis #
###########################################


exprMat <- read.delim('DEG/He_GEM.txt', header=TRUE)

metadata <- data.frame(
  sample    = colnames(exprMat),
  condition = c('adult', 'adult', 'adult', 'aged', 'aged', 'aged')
)
metadata$condition <- factor(metadata$condition, levels=c('adult', 'aged'))

# 2️⃣ Construct DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = exprMat,
  colData   = metadata,
  design    = ~ condition
)

# Pre-filter low-count genes
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Extract results
res <- results(dds, alpha=0.05)
res <- lfcShrink(dds, coef="condition_aged_vs_adult", type="apeglm")

# Order by adjusted p-value
resOrdered <- as.data.frame(res[order(res$padj), ])

# Map ENSMUSG to gene names with biomaRt
ensembl_values <- rownames(resOrdered)
gene_map <- getBM(
    attributes = c('ensembl_gene_id_version', 'mgi_symbol'),
    filters    = 'ensembl_gene_id_version',
    values     = ensembl_values,
    mart       = mart
)
# rename to match downstream merge by.y='ensembl_gene_id_version'
resOrdered <- merge(resOrdered, gene_map, by.x='row.names', by.y='ensembl_gene_id_version', all.x=TRUE)
colnames(resOrdered)[1] <- 'ensembl_gene_id_version'
colnames(resOrdered)[ncol(resOrdered)] <- 'gene_name'

# Write results to file
write.table(resOrdered, file="DESeq2/He_9w_vs_78w.txt", sep="\t", quote=FALSE, row.names=FALSE)




####################
# Read in TAC data #
####################


foo <- read.delim('adult/Sample_22L011981/22L011981_stranded_reverse.count', header=FALSE,
                         col.names=c('geneID', 'adult_1'),
                         stringsAsFactors=FALSE)
foo <- foo[grep('^ENS', foo$V1), ]
sum(foo$V2)

conda activate RNAseq
cd adult/Sample_22L011981
htseq-count -s no -r pos -f bam 22L011981_Aligned.sortedByCoord.out.bam ../../../GRCm39/gencode.vM37.primary_assembly.annotation.gtf > test_noStrand.count