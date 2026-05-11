library(tidyverse)
library(data.table)
library(VennDiagram)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(ggplot2)
library(ggrepel)

## GTF file for TPM calculation
gtf_file <- '/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39/20250512_RefSeq/annotation.gtf'

# Parse gene lengths from GTF (union of exon lengths per gene)
get_gene_lengths <- function(gtf_file) {
    gtf <- read.table(gtf_file, header = FALSE, sep = '\t', comment.char = '#',
                      col.names = c('seqname','source','feature','start','end',
                                    'score','strand','frame','attribute'))
    exons <- gtf[gtf$feature == 'exon', ]

    # Extract gene_id from attribute field (handles both quoted and unquoted values)
    exons$gene_id <- sub('.*gene_id "?([^";]+)"?;.*', '\\1', exons$attribute)

    # For each gene, compute union of exon lengths (non-overlapping)
    gene_lengths <- tapply(seq_len(nrow(exons)), exons$gene_id, function(idx) {
        starts <- exons$start[idx]
        ends   <- exons$end[idx]
        # reduce to non-overlapping intervals
        ord    <- order(starts)
        starts <- starts[ord]; ends <- ends[ord]
        merged_len <- 0L
        cur_start  <- starts[1]; cur_end <- ends[1]
        for (i in seq_along(starts)) {
            if (starts[i] <= cur_end) {
                cur_end <- max(cur_end, ends[i])
            } else {
                merged_len <- merged_len + (cur_end - cur_start + 1L)
                cur_start  <- starts[i]; cur_end <- ends[i]
            }
        }
        merged_len + (cur_end - cur_start + 1L)
    })
    setNames(as.integer(gene_lengths), names(gene_lengths))
}

# Calculate TPM from a count matrix (genes x samples)
# counts: matrix or data.frame with gene names as rownames
# gene_lengths: named integer vector of gene lengths in bp
calculate_tpm <- function(counts, gene_lengths) {
    counts <- as.matrix(counts)
    shared <- intersect(rownames(counts), names(gene_lengths))
    counts  <- counts[shared, , drop = FALSE]
    lengths <- gene_lengths[shared]

    rpk <- counts / (lengths / 1e3)           # reads per kilobase
    tpm <- sweep(rpk, 2, colSums(rpk) / 1e6, '/') # normalise per sample
    tpm
}

# Usage:
gene_lengths <- get_gene_lengths(gtf_file)

### Read in Adult MCL metadata and counts
adult_MCL_metadata <- read.csv('cardiac_RNAseq/SraRunTable.csv')
colnames(adult_MCL_metadata)[1] <- 'Sample_Name'
adult_MCL_metadata$condition <- ifelse(grepl('Xist', adult_MCL_metadata$genotype), 'Xist', 'WT')
adult_MCL_metadata$AGE <- 'adult'
# new column where if tissue mapps to Macrophages = MP, Cardiomyocytes = CM, Cardiac fibroblasts = CF, Endothelial cells = EC
adult_MCL_metadata$cell_type <- ifelse(grepl('Macrophages', adult_MCL_metadata$tissue), 'MP',
                                ifelse(grepl('Cardiomyocytes', adult_MCL_metadata$tissue), 'CM',
                                       ifelse(grepl('Cardiac fibroblasts', adult_MCL_metadata$tissue), 'CF',
                                              ifelse(grepl('Endothelial cells', adult_MCL_metadata$tissue), 'EC', NA))))
# list all count files and read into a list
facs_files <- list.files('cardiac_RNAseq', pattern = 'stranded_no.count', full.names = TRUE, recursive = TRUE)
count.list <- lapply(facs_files, function(f) {
    sample <- strsplit(f, '/')[[1]][2]
    tmp <- read.delim(f, header = FALSE)
    rownames(tmp) <- tmp$V1
    colnames(tmp) <- c('gene_id', sample)
    return(tmp)
})
# merge all count files into one data frame
facs_counts <- Reduce(function(x, y) merge(x, y, by = 'gene_id', all = TRUE), count.list)
rownames(facs_counts) <- facs_counts$gene_id
facs_counts <- facs_counts[grep('^__', rownames(facs_counts), invert=TRUE), ]
facs_counts[is.na(facs_counts)] <- 0
# Write raw counts to file for later use
write.table(facs_counts, 'cardiac_RNAseq/ExprMat.txt', sep = '\t', quote = FALSE, row.names = TRUE)
# Calculate TPM
facs_counts_tpm <- calculate_tpm(facs_counts[,-1], gene_lengths)
write.table(facs_counts_tpm, 'cardiac_RNAseq/ExprMat_TPM.txt', sep = '\t', quote = FALSE, row.names = TRUE)

# subset for Xist and Create long data for plotting
xist_counts_adult <- as.data.frame(facs_counts_tpm[rownames(facs_counts_tpm) == 'Xist', , drop = FALSE])
xist_counts_adult_long <- pivot_longer(xist_counts_adult, cols = everything(), names_to = 'sample', values_to = 'tpm')
xist_counts_adult_long <- merge(xist_counts_adult_long, adult_MCL_metadata[,c('Sample_Name', 'cell_type', 'AGE', 'condition', 'sex')], by.x = 'sample', by.y = 'Sample_Name')
xist_counts_adult_long <- subset(xist_counts_adult_long, sex == 'female' & condition == 'Xist')
xist_counts_adult_long$tpm <- as.numeric(xist_counts_adult_long$tpm)


### Read in Aged MCL metadata and counts
aged_MCL_metadata <- read.csv('aged_cardiac_RNAseq/Project_699_lims.csv')
aged_MCL_metadata$cell_type <- strsplit(aged_MCL_metadata$FID_comment, '_') |> sapply(`[`, 1)
aged_MCL_metadata$sex <- strsplit(aged_MCL_metadata$FID_comment, '_') |> sapply(`[`, 6)
aged_MCL_metadata$genotype <- strsplit(aged_MCL_metadata$FID_comment, '_') |> sapply(`[`, 4)
aged_MCL_metadata$AGE <- 'aged'

# list all count files and read into a list
aged_facs_files <- list.files('aged_cardiac_RNAseq', pattern = 'stranded_no.count', full.names = TRUE, recursive = TRUE)
count.list <- lapply(aged_facs_files, function(f) {
    sample <- strsplit(f, '/')[[1]][2]
    tmp <- read.delim(f, header = FALSE)
    rownames(tmp) <- tmp$V1
    colnames(tmp) <- c('gene_id', sample)
    return(tmp)
})
# merge all count files into one data frame
aged_facs_counts <- Reduce(function(x, y) merge(x, y, by = 'gene_id', all = TRUE), count.list)
rownames(aged_facs_counts) <- aged_facs_counts$gene_id
aged_facs_counts <- aged_facs_counts[grep('^__', rownames(aged_facs_counts), invert=TRUE), ]
aged_facs_counts[is.na(aged_facs_counts)] <- 0
# Write raw counts to file for later use
write.table(aged_facs_counts, 'aged_cardiac_RNAseq/ExprMat.txt', sep = '\t', quote = FALSE, row.names = TRUE)
# Calculate TPM
aged_facs_counts_tpm <- calculate_tpm(aged_facs_counts[,-1], gene_lengths)
write.table(aged_facs_counts_tpm, 'aged_cardiac_RNAseq/ExprMat_TPM.txt', sep = '\t', quote = FALSE, row.names = TRUE)

# subset for Xist and Create long data for plotting
xist_counts_aged <- as.data.frame(aged_facs_counts_tpm[rownames(aged_facs_counts_tpm) == 'Xist', , drop = FALSE])
xist_counts_aged_long <- pivot_longer(xist_counts_aged, cols = everything(), names_to = 'sample', values_to = 'tpm')
xist_counts_aged_long <- merge(xist_counts_aged_long, aged_MCL_metadata[,c('Sample_Name', 'cell_type', 'sex', 'genotype', 'AGE')], by.x = 'sample', by.y = 'Sample_Name')
xist_counts_aged_long <- subset(xist_counts_aged_long, sex == 'XX' & genotype == 'XBxC')
xist_counts_aged_long$tpm <- as.numeric(xist_counts_aged_long$tpm)


### Read in TAC MCL metadata and counts
TAC_MCL_metadata <- read.csv('TAC_cardiac_RNAseq/Project_1050_lims_fixed.csv')
TAC_MCL_metadata$cell_type <- strsplit(TAC_MCL_metadata$FID_comment, '_') |> sapply(`[`, 1)
TAC_MCL_metadata$sex <- ifelse(grepl('XX', TAC_MCL_metadata$FID_comment), 'XX', 'XY')
TAC_MCL_metadata$AGE <- ifelse(grepl('Sham', TAC_MCL_metadata$FID_comment), 'Sham', 'TAC')

# list all count files and read into a list
TAC_facs_files <- list.files('TAC_cardiac_RNAseq', pattern = 'stranded_no.count', full.names = TRUE, recursive = TRUE)
count.list <- lapply(TAC_facs_files, function(f) {
    sample <- strsplit(f, '/')[[1]][2]
    tmp <- read.delim(f, header = FALSE)
    rownames(tmp) <- tmp$V1
    colnames(tmp) <- c('gene_id', sample)
    return(tmp)
})
# merge all count files into one data frame
TAC_facs_counts <- Reduce(function(x, y) merge(x, y, by = 'gene_id', all = TRUE), count.list)
rownames(TAC_facs_counts) <- TAC_facs_counts$gene_id
TAC_facs_counts <- TAC_facs_counts[grep('^__', rownames(TAC_facs_counts), invert=TRUE), ]
TAC_facs_counts[is.na(TAC_facs_counts)] <- 0
# Write raw counts to file for later use
write.table(TAC_facs_counts, 'TAC_cardiac_RNAseq/ExprMat.txt', sep = '\t', quote = FALSE, row.names = TRUE)
# Calculate TPM
TAC_facs_counts_tpm <- calculate_tpm(TAC_facs_counts[,-1], gene_lengths)
write.table(TAC_facs_counts_tpm, 'TAC_cardiac_RNAseq/ExprMat_TPM.txt', sep = '\t', quote = FALSE, row.names = TRUE)

# subset for Xist and Create long data for plotting
xist_counts_TAC <- as.data.frame(TAC_facs_counts_tpm[rownames(TAC_facs_counts_tpm) == 'Xist', , drop = FALSE])
xist_counts_TAC_long <- pivot_longer(xist_counts_TAC, cols = everything(), names_to = 'sample', values_to = 'tpm')
xist_counts_TAC_long <- merge(xist_counts_TAC_long, TAC_MCL_metadata[,c('Sample_Name', 'cell_type', 'sex', 'AGE')], by.x = 'sample', by.y = 'Sample_Name')
xist_counts_TAC_long <- subset(xist_counts_TAC_long, sex == 'XX')

# Combine adult, aged, and TAC data for plotting
combined_xist_counts <- rbind(
    xist_counts_adult_long[, c('tpm', 'sample', 'cell_type', 'AGE')], 
    xist_counts_aged_long[, c('tpm', 'sample', 'cell_type', 'AGE')], 
    xist_counts_TAC_long[, c('tpm', 'sample', 'cell_type', 'AGE')])
combined_xist_counts$cell_type <- factor(combined_xist_counts$cell_type, levels = c('MP', 'CM', 'CF', 'EC'))
combined_xist_counts$AGE <- factor(combined_xist_counts$AGE, levels = c('adult', 'aged', 'Sham', 'TAC'))


library(ggpubr)
pdf('combined_Xist_boxplot.pdf')
my_comparisons <- list( c("adult", "aged"), c("Sham", "TAC") )
ggboxplot(combined_xist_counts, x = "AGE", y = "tpm",
          color = "AGE", palette = "jco",
          add = "jitter", add.params = list(color = "black", size = 1.5, alpha = 0.7)) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+ 
  stat_compare_means(label.y = 50) +
  ylab("Xist TPM") +
  facet_wrap(~cell_type) +
  theme(legend.position = "none") +
  xlab('')
dev.off()


pdf('combined_Xist_boxplot_no_celltype.pdf')
my_comparisons <- list( c("adult", "aged"), c("Sham", "TAC") )
ggboxplot(combined_xist_counts, x = "AGE", y = "tpm",
          color = "AGE", palette = "jco",
          add = "jitter", add.params = list(color = "black", size = 1.5, alpha = 0.7)) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+ 
  stat_compare_means(label.y = 50) +
  ylab("Xist TPM") +
  theme(legend.position = "none") +
  xlab('')
dev.off()

## Repeat for a housekeeping gene - Gapdh
gapdh_counts_adult <- as.data.frame(facs_counts_tpm[rownames(facs_counts_tpm) == 'Gapdh', , drop = FALSE])
gapdh_counts_adult_long <- pivot_longer(gapdh_counts_adult, cols = everything(), names_to = 'sample', values_to = 'tpm')
gapdh_counts_adult_long <- merge(gapdh_counts_adult_long, adult_MCL_metadata[,c('Sample_Name', 'cell_type', 'AGE', 'condition', 'sex')], by.x = 'sample', by.y = 'Sample_Name')
gapdh_counts_adult_long <- subset(gapdh_counts_adult_long, sex == 'female' & condition == 'Xist')
gapdh_counts_adult_long$tpm <- as.numeric(gapdh_counts_adult_long$tpm)

gapdh_counts_aged <- as.data.frame(aged_facs_counts_tpm[rownames(aged_facs_counts_tpm) == 'Gapdh', , drop = FALSE])
gapdh_counts_aged_long <- pivot_longer(gapdh_counts_aged, cols = everything(), names_to = 'sample', values_to = 'tpm')
gapdh_counts_aged_long <- merge(gapdh_counts_aged_long, aged_MCL_metadata[,c('Sample_Name', 'cell_type', 'sex', 'genotype', 'AGE')], by.x = 'sample', by.y = 'Sample_Name')
gapdh_counts_aged_long <- subset(gapdh_counts_aged_long, sex == 'XX' & genotype == 'XBxC')
gapdh_counts_aged_long$tpm <- as.numeric(gapdh_counts_aged_long$tpm)

gapdh_counts_TAC <- as.data.frame(TAC_facs_counts_tpm[rownames(TAC_facs_counts_tpm) == 'Gapdh', , drop = FALSE])
gapdh_counts_TAC_long <- pivot_longer(gapdh_counts_TAC, cols = everything(), names_to = 'sample', values_to = 'tpm')
gapdh_counts_TAC_long <- merge(gapdh_counts_TAC_long, TAC_MCL_metadata[,c('Sample_Name', 'cell_type', 'sex', 'AGE')], by.x = 'sample', by.y = 'Sample_Name')
gapdh_counts_TAC_long <- subset(gapdh_counts_TAC_long, sex == 'XX')
gapdh_counts_TAC_long$tpm <- as.numeric(gapdh_counts_TAC_long$tpm) 

combined_gapdh_counts <- rbind(
    gapdh_counts_adult_long[, c('tpm', 'sample', 'cell_type', 'AGE')], 
    gapdh_counts_aged_long[, c('tpm', 'sample', 'cell_type', 'AGE')], 
    gapdh_counts_TAC_long[, c('tpm', 'sample', 'cell_type', 'AGE')])
combined_gapdh_counts$cell_type <- factor(combined_gapdh_counts$cell_type, levels = c('MP', 'CM', 'CF', 'EC'))
combined_gapdh_counts$AGE <- factor(combined_gapdh_counts$AGE, levels = c('adult', 'aged', 'Sham', 'TAC'))

pdf('combined_Gapdh_boxplot.pdf')
my_comparisons <- list( c("adult", "aged"), c("Sham", "TAC") )
ggboxplot(combined_gapdh_counts, x = "AGE", y = "tpm",
          color = "AGE", palette = "jco",
          add = "jitter", add.params = list(color = "black", size = 1.5, alpha = 0.7)) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+ 
  stat_compare_means(label.y = 50) +
  ylab("Gapdh TPM") +
  facet_wrap(~cell_type) +
  theme(legend.position = "none") +
  xlab('')
dev.off()

library(DESeq2)