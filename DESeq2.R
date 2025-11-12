library(dplyr)
library(tidyr)
library(DESeq2)
library(apeglm)
library(org.Mm.eg.db)
library(AnnotationDbi)


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
adult_metadata <- read.delim('adult_metadata.csv', sep=',')
aged_metadata  <- read.delim('aged_metadata.csv', sep=',')
metadata <- rbind(adult_metadata, aged_metadata)
metadata <- metadata[metadata$Sample_ID %in% samples & metadata$FID_comment != "", ]

metadata$tissue <- unlist(lapply(metadata$FID_comment, function(x) strsplit(x, '_')[[1]][1]))
metadata$age  <- lapply(metadata$FID_comment, function(x) strsplit(x, '_')[[1]][2])
metadata$tissue_age <- paste0(metadata$tissue, '_', metadata$age)

lapply(unique(metadata$tissue), function(x) {
    tmp <- subset(metadata, tissue == x)

    adult_rep1 <- read.delim(paste0('adult/', tmp$Sample_ID[1], '/', tmp$Sample_Name[1], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[1], '_rep1')))
    adult_rep2 <- read.delim(paste0('adult/', tmp$Sample_ID[2], '/', tmp$Sample_Name[2], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[2], '_rep2')))
    adult_rep3 <- read.delim(paste0('adult/', tmp$Sample_ID[3], '/', tmp$Sample_Name[3], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[3], '_rep3')))
    aged_rep1  <- read.delim(paste0('aged/', tmp$Sample_ID[4], '/', tmp$Sample_Name[4], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[4], '_rep1')))
    aged_rep2  <- read.delim(paste0('aged/', tmp$Sample_ID[5], '/', tmp$Sample_Name[5], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[5], '_rep2')))
    aged_rep3  <- read.delim(paste0('aged/', tmp$Sample_ID[6], '/', tmp$Sample_Name[6], '_stranded_reverse.count'), header=FALSE, col.names=c('geneID', paste0(tmp$tissue_age[6], '_rep3')))
    merged <- Reduce(function(x, y) merge(x, y, by='geneID', all=FALSE), list(adult_rep1, adult_rep2, adult_rep3, aged_rep1, aged_rep2, aged_rep3))
    merged <- merged[grep('^_', merged$geneID, invert=TRUE), ]
    rownames(merged) <- merged$geneID
    merged <- merged[, 2:(ncol(merged))]
    write.table(merged, file=paste0('DEG/', x, '_GEM.txt'), sep='\t', quote=FALSE)
})
names(tissue_GEM) <- unique(metadata$tissue)
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

# Construct DESeq2 object
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

# Map REFSEQ to gene names with org.Mm.eg.db
refseq_values <- sub("\\.\\d+$", "", rownames(resOrdered))

sym <- mapIds(org.Mm.eg.db, keys = refseq_values, keytype = "REFSEQ",
              column = "SYMBOL", multiVals = "first")
resOrdered$gene_symbol <- unname(sym[refseq_values])

# Write results to file
write.table(resOrdered, file="DESeq2/He_9w_vs_78w.txt", sep="\t", quote=FALSE, row.names=FALSE)




####################
# Read in TAC data #
####################
metadata <- read.delim('metadata.csv', sep=',' )

count_files <- list.files('.', pattern='_stranded_reverse.count', recursive=TRUE, full.names=TRUE)

GEM <- lapply(count_files, function(x) {
  sample_name <- gsub('_stranded_reverse.count', '', basename(x))
  df <- read.delim(x, header=FALSE, col.names=c('geneID', paste0('Sample_', sample_name)))
  return(df)
})
GEM <- Reduce(function(x, y) merge(x, y, by='geneID', all=TRUE), GEM)
GEM <- GEM[grep('^_', GEM$geneID, invert=TRUE), ]
metadata <- metadata[grep('XistBxC_-\\+_XX', metadata$FID_comment), ]

GEM <- GEM[, c(1, which(colnames(GEM) %in% metadata$Sample_ID))]

metadata <- metadata[match(colnames(GEM[,-1]), metadata$Sample_ID), ]
rownames(GEM) <- GEM$geneID
GEM <- GEM[, -1]
write.table(GEM, file='DEG/TAC_GEM.txt', sep='\t', quote=FALSE, row.names=FALSE)



metadata2 <- data.frame(
  sample    = colnames(GEM),
  condition = c('TAC', 'TAC', 'TAC', 'sham', 'sham', 'sham')
)
metadata2$condition <- factor(metadata2$condition, levels=c('TAC', 'sham'))

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = GEM,
  colData   = metadata2,
  design    = ~ condition
)

# Pre-filter low-count genes
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
# Run DESeq2 pipeline
dds <- DESeq(dds)
# Extract results
res <- results(dds, alpha=0.05)
res <- lfcShrink(dds, coef="condition_sham_vs_TAC", type="apeglm")
# Order by adjusted p-value
resOrdered <- as.data.frame(res[order(res$padj), ])
# Map REFSEQ to gene names with org.Mm.eg.db
refseq_values <- sub("\\.\\d+$", "", rownames(resOrdered))

sym <- mapIds(org.Mm.eg.db, keys = refseq_values, keytype = "REFSEQ",
              column = "SYMBOL", multiVals = "first")
resOrdered$gene_symbol <- unname(sym[refseq_values])
subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)

# Write results to file
write.table(resOrdered, file="DEG/TAC_vs_sham.txt", sep="\t", quote=FALSE, row.names=FALSE)


adult_aged <- read.delim('../adult_aged_bodymap/DESeq2/He_9w_vs_78w.txt', header=TRUE)

combined <- merge(adult_aged, resOrdered, by='gene_name', suffixes=c('_adult_aged', '_TAC_sham'))

combined_sig <- subset(combined, padj_adult_aged < 0.05 & padj_TAC_sham < 0.05)

# Correlate on -log10(padj)
pdf('../LINKS_study/figures/DEG_correlation_adult_aged_vs_TAC_sham.pdf')
plot(-log10(combined_sig$padj_adult_aged), -log10(combined_sig$padj_TAC_sham),
     xlab='-log10(padj) Adult Aged', ylab='-log10(padj) TAC Sham',
     main='DEG Correlation: Adult Aged vs TAC Sham')
abline(0, 1, col='red', lty=2)
dev.off()
