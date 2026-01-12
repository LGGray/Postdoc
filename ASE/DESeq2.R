library(dplyr)
library(tidyr)
library(DESeq2)
library(apeglm)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)


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
# Add gene
resOrdered$gene <- rownames(resOrdered)

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
write.table(GEM, file='DEG/TAC_GEM.txt', sep='\t', quote=FALSE, row.names=TRUE)



metadata2 <- data.frame(
  sample    = colnames(GEM),
  condition = c('TAC', 'TAC', 'TAC', 'sham', 'sham', 'sham')
)
metadata2$condition <- factor(metadata2$condition, levels=c('sham', 'TAC'))

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
res <- lfcShrink(dds, coef="condition_TAC_vs_sham", type="apeglm")
# Order by adjusted p-value
resOrdered <- as.data.frame(res[order(res$padj), ])
# Add gene
resOrdered$gene <- rownames(resOrdered)

# Write results to file
write.table(resOrdered, file="DEG/TAC_vs_sham.txt", sep="\t", quote=FALSE, row.names=FALSE)

##################################################################
# Correlation of sig genes between aged vs adult and TAC vs sham #
##################################################################

adult_aged <- read.delim('../adult_aged_bodymap/DESeq2/He_9w_vs_78w.txt', header=TRUE)
adult_aged_sig <- subset(adult_aged, padj < 0.05)
TAC_sham  <- read.delim('../F1_TAC_Sarah/DEG/TAC_vs_sham.txt', header=TRUE)
TAC_sham_sig <- subset(TAC_sham, padj < 0.05)

combined <- merge(
  adult_aged_sig[, c('gene', 'log2FoldChange', 'padj')],
  TAC_sham_sig[, c('gene', 'log2FoldChange', 'padj')],
  by='gene',
  suffixes=c('_adult_aged', '_TAC_sham')
)
correlation_res <- cor.test(combined$log2FoldChange_adult_aged, combined$log2FoldChange_TAC_sham, method='spearman')

# Correlate on -log10(padj)
pdf('../LINKS_study/figures/DEG_correlation_adult_aged_vs_TAC_sham.pdf')
plot(combined$log2FoldChange_adult_aged, combined$log2FoldChange_TAC_sham,
     xlab='Aged vs Adult log2FC', ylab='TAC vs Sham log2FC',
     main='')
legend('topright', legend=paste0('Spearman rho = ', round(correlation_res$estimate, 3), '\nP-value = ***'), bty='n')
abline(0, 1, col='red', lty=2)
dev.off()

set.seed(42)
clusters <- kmeans(combined[, c('log2FoldChange_adult_aged', 'log2FoldChange_TAC_sham')], centers=4)$cluster

pdf('../LINKS_study/figures/DEG_adult_aged_vs_TAC_sham_clusters.pdf')
ggplot(combined, aes(x=log2FoldChange_adult_aged, y=log2FoldChange_TAC_sham, colour = factor(clusters))) +
  geom_point(alpha = 0.6) +
  labs(colour = "Cluster")
dev.off()

genes_in_cluster <- split(combined$gene, clusters)

lapply(genes_in_cluster, function(x){
  x[x %in% c('Spaar', 'Lgr6', 'Lmh1', 'Hmgb1', 'Katnal1')]
})

library(fgsea)

pathways <- gmtPathways('../pathways/m5.go.bp.v2025.1.Mm.symbols.gmt')
overrepresentation <- lapply(genes_in_cluster, function(genes) {
  fora(genes=genes, universe=combined$gene, pathways=pathways)
})

lapply(overrepresentation, function(x) subset(x, padj < 0.1)[, c('pathway', 'padj', 'overlapGenes')])



# Visualise DEGs

adult_age <- read.delim('adult_aged_bodymap/DESeq2/He_9w_vs_78w.txt', header=TRUE)

pdf('LINKS_study/figures/DEG_adult_aged_volcano.pdf')
adult_age$color <- 'grey'
adult_age$color[adult_age$log2FoldChange > 1 & adult_age$padj < 0.05] <- 'red'
adult_age$color[adult_age$log2FoldChange < -1 & adult_age$padj < 0.05] <- 'blue'

top_up <- adult_age %>%
  filter(log2FoldChange > 1 & padj < 0.05) %>%
  top_n(5, wt=log2FoldChange)
top_down <- adult_age %>%
  filter(log2FoldChange < -1 & padj < 0.05) %>%
  top_n(-5, wt=log2FoldChange)

adult_age$label <- ''
adult_age$label[adult_age$gene %in% top_up$gene] <- top_up$gene
adult_age$label[adult_age$gene %in% top_down$gene] <- top_down$gene
if ('Pah' %in% adult_age$gene) {
  adult_age$label[adult_age$gene == 'Pah'] <- 'Pah'
}

ggplot(adult_age, aes(x=log2FoldChange, y=-log10(padj), color=color)) +
  geom_point(alpha=0.4) +
  scale_color_identity() +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3,
    color = 'black',
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  theme_minimal() +
  xlab('Log2 Fold Change (Aged vs Adult)') +
  ylab('-Log10 Adjusted P-value') +
  geom_hline(yintercept=-log10(0.05), linetype='dashed', color='black') +
  geom_vline(xintercept=c(-1, 1), linetype='dashed', color='black')
dev.off()

##################
sham_TAC <- read.delim('F1_TAC_Sarah/DEG/TAC_vs_sham.txt', header=TRUE)
pdf('LINKS_study/figures/DEG_TAC_vs_sham_volcano.pdf')
sham_TAC$color <- 'grey'
sham_TAC$color[sham_TAC$log2FoldChange > 1 & sham_TAC$padj < 0.05] <- 'red'
sham_TAC$color[sham_TAC$log2FoldChange < -1 & sham_TAC$padj < 0.05] <- 'blue'
top_up <- sham_TAC %>%
  filter(log2FoldChange > 1 & padj < 0.05) %>%
  top_n(5, wt=log2FoldChange)
top_down <- sham_TAC %>%
  filter(log2FoldChange < -1 & padj < 0.05) %>%
  top_n(-5, wt=log2FoldChange)



sham_TAC$label <- ''
sham_TAC$label[sham_TAC$gene %in% top_up$gene] <- top_up$gene
sham_TAC$label[sham_TAC$gene %in% top_down$gene] <- top_down$gene
sham_TAC[sham_TAC$gene %in% c('Nppa','Nppb','Myh7','Postn','Col1a1'), ]$label <- c('Nppa','Nppb','Myh7','Postn','Col1a1')
ggplot(sham_TAC, aes(x=log2FoldChange, y=-log10(padj), color=color)) +
  geom_point(alpha=0.4) +
  scale_color_identity() +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3,
    color = 'black',
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  theme_minimal() +
  xlab('Log2 Fold Change (TAC vs Sham)') +
  ylab('-Log10 Adjusted P-value') +
  geom_hline(yintercept=-log10(0.05), linetype='dashed', color='black') +
  geom_vline(xintercept=c(-1, 1), linetype='dashed', color='black')
dev.off()


noncoding_coding <- read.delim('LINKS_study/Refseq_coding_noncoding.txt', header = TRUE)
nrow(subset(sham_TAC, gene %in% subset(noncoding_coding, coding_noncoding == "coding")$name & log2FoldChange > 1 & padj < 0.05))

#########################################################
# Plot correlation of adult vs age and sham vs TAC DEGs # 
#########################################################
combined <- merge(
  adult_age[, c('gene', 'log2FoldChange', 'padj')],
  sham_TAC[, c('gene', 'log2FoldChange', 'padj')],
  by='gene',
  suffixes=c('_adult_age', '_TAC_sham')
)
combined_deg <- subset(combined, padj_adult_age < 0.05 & padj_TAC_sham < 0.05)

cor_res <- cor.test(combined_deg$log2FoldChange_adult_age, combined_deg$log2FoldChange_TAC_sham, method='spearman')
combined_deg$label <- dplyr::case_when(
  combined_deg$log2FoldChange_adult_age > 2 & combined_deg$log2FoldChange_TAC_sham > 1  ~ combined_deg$gene,
  combined_deg$log2FoldChange_adult_age > 1.5 & combined_deg$log2FoldChange_TAC_sham < -1.5 ~ combined_deg$gene,
  TRUE ~ ''
)

pdf('LINKS_study/figures/scatterplot_adult_age_vs_TAC_sham_all_DEGs.pdf')
ggplot(combined_deg, aes(x=log2FoldChange_adult_age, y=log2FoldChange_TAC_sham)) +
  geom_point(alpha=0.6) +
  theme_minimal() +
  xlab('Log2 Fold Change (Aged vs Adult)') +
  ylab('Log2 Fold Change (TAC vs Sham)') +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3,
    color = 'black',
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_vline(xintercept=0, linetype='dashed', color='black') +
  annotate('text', x=Inf, y=Inf, 
           label=paste0('Spearman rho = ', round(cor_res$estimate, 3), '\nP-value = ***'), 
           hjust=1.1, vjust=1.1, size=3) +
  # Colour by same direction or not
  aes(color = ifelse((log2FoldChange_adult_age > 0 & log2FoldChange_TAC_sham > 0) | (log2FoldChange_adult_age < 0 & log2FoldChange_TAC_sham < 0), 'Same direction', 'Opposite direction')) +
  scale_color_manual(values = c('Same direction' = 'blue', 'Opposite direction' = 'red')) +
  labs(color = 'Direction') +
  theme(legend.position = 'none')
dev.off()

##############################
# Combine adult and age heart samples and TAC and sham heart samples to visualise in PCA #
##############################
adult_age <- read.delim('adult_aged_bodymap/DEG/He_GEM.txt', header=TRUE)
adult_age$geneID <- rownames(adult_age)
TAC_sham  <- read.delim('F1_TAC_Sarah/DEG/TAC_GEM.txt', header=TRUE)
colnames(TAC_sham) <- c('TAC_rep1', 'TAC_rep2', 'TAC_rep3', 'sham_rep1', 'sham_rep2', 'sham_rep3')
TAC_sham$geneID <- rownames(TAC_sham)
combined <- merge(adult_age, TAC_sham, by='geneID', all=TRUE)




# Build sample metadata
sample_names <- colnames(combined)[-1]
study <- ifelse(grepl('_9w_', sample_names) | grepl('_78w_', sample_names), 'adult_aged_bodymap', 'F1_TAC_Sarah')
condition <- dplyr::case_when(
  grepl('_9w_', sample_names) ~ 'adult',
  grepl('_78w_', sample_names)  ~ 'aged',
  grepl('^TAC_', sample_names)   ~ 'TAC',
  grepl('^sham_', sample_names)  ~ 'sham',
  TRUE ~ 'unknown'
)
coldata <- data.frame(row.names = sample_names, study=factor(study), condition=factor(condition))

# Prepare count matrix
counts_mat <- as.matrix(combined[, -1])
rownames(counts_mat) <- combined$geneID
counts_mat[is.na(counts_mat)] <- 0
storage.mode(counts_mat) <- "integer"

# DESeq2 object with batch term (not for DE testing here, just for normalization metadata)
dds_combined <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = coldata,
  design    = ~ 1
)

# Normalize and variance-stabilize
dds_combined <- estimateSizeFactors(dds_combined)
vst_mat <- assay(vst(dds_combined, blind=TRUE))

# Batch correction for PCA (protect condition)
library(limma)
vst_bc <- removeBatchEffect(vst_mat, batch=coldata$study,
                            design = model.matrix(~ condition, data=coldata))

# Write batch-corrected VST matrix
write.table(vst_bc, file='LINKS_study/batch_corrected_VST_adult_aged_TAC_sham.txt', sep='\t', quote=FALSE, row.names=TRUE)

# PCA
pca <- prcomp(t(vst_bc), scale.=FALSE)
pca_df <- data.frame(pca$x[,1:2], condition=coldata$condition, study=coldata$study, sample=rownames(coldata))

pdf('LINKS_study/figures/PCA_adult_aged_TAC_sham_batch_corrected.pdf')
ggplot(pca_df, aes(PC1, PC2, color=condition)) +
  geom_point(size=3, alpha=0.8) +
  theme_minimal() +
  labs(
       x=paste0('PC1 (', round(100*summary(pca)$importance[2,1],1), '%)'),
       y=paste0('PC2 (', round(100*summary(pca)$importance[2,2],1), '%)'))
dev.off()

head(sort(abs(pca$rotation[,1]), decreasing=TRUE), 20)
head(sort(abs(pca$rotation[,2]), decreasing=TRUE), 20)



library(ComplexHeatmap)
library(circlize)
pdf('LINKS_study/figures/heatmap_adult_aged_TAC_sham_batch_corrected.pdf')
Heatmap(vst_bc,
        name = "VST (batch-corrected)",
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_split = coldata$condition,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)
dev.off()


# Plot a boxplot of Hmgb1 expression across conditions
# Use batch-corrected VST values for cross-study comparison
vst_bc <- read.delim('LINKS_study/batch_corrected_VST_adult_aged_TAC_sham.txt', header=TRUE)
hmgb1_expr <- vst_bc[rownames(vst_bc) == 'Hmgb1', ]

# Convert to long format
hmgb1_expr_long <- data.frame(
  sample = colnames(hmgb1_expr),
  expression = as.numeric(hmgb1_expr[1, ])
)

hmgb1_expr_long$condition <- dplyr::case_when(
  grepl('_9w_', hmgb1_expr_long$sample) ~ 'adult',
  grepl('_78w_', hmgb1_expr_long$sample)  ~ 'aged',
  grepl('^TAC_', hmgb1_expr_long$sample)   ~ 'TAC',
  grepl('^sham_', hmgb1_expr_long$sample)  ~ 'sham',
  TRUE ~ 'unknown'
)

# Get DESeq2 results for annotation
adult_age_deseq <- read.delim('adult_aged_bodymap/DESeq2/He_9w_vs_78w.txt', header=TRUE)
tac_sham_deseq <- read.delim('F1_TAC_Sarah/DEG/TAC_vs_sham.txt', header=TRUE)
hmgb1_adult_age <- subset(adult_age_deseq, gene == 'Hmgb1')
hmgb1_tac_sham <- subset(tac_sham_deseq, gene == 'Hmgb1')

# Set condition factor levels for proper ordering
hmgb1_expr_long$condition <- factor(hmgb1_expr_long$condition, levels=c('adult', 'aged', 'sham', 'TAC'))

# Calculate y positions for significance bars
y_max <- max(hmgb1_expr_long$expression, na.rm=TRUE)
y_range <- diff(range(hmgb1_expr_long$expression, na.rm=TRUE))
bar_height <- y_max + 0.1 * y_range

sample_colours <- c("adult" = "#E69F00", "aged" = "#56B4E9", "TAC" = "#009E73", "sham" = "#CC79A7", "FACS" = "#A67C8C")
pdf('LINKS_study/figures/Hmgb1_expression_boxplot.pdf', width=10, height=6)
ggplot(hmgb1_expr_long, aes(x=condition, y=expression, fill=condition)) +
  geom_boxplot() +
  geom_jitter(width=0.2, size=2, alpha=0.8) +
  theme_minimal() +
  ylab('Hmgb1 Expression (VST batch-corrected)') +
  xlab('') +
  scale_fill_manual(values = sample_colours) +
  theme(legend.position = 'none') +
  # Adult vs Aged significance bar
  geom_segment(aes(x=1, xend=2, y=bar_height, yend=bar_height), color='black') +
  geom_segment(aes(x=1, xend=1, y=bar_height, yend=bar_height - 0.02*y_range), color='black') +
  geom_segment(aes(x=2, xend=2, y=bar_height, yend=bar_height - 0.02*y_range), color='black') +
  annotate('text', x=1.5, y=bar_height + 0.03*y_range, 
           label=paste0('padj = ', round(hmgb1_adult_age$padj, 3)), size=3.5) +
  # Sham vs TAC significance bar
  geom_segment(aes(x=3, xend=4, y=bar_height, yend=bar_height), color='black') +
  geom_segment(aes(x=3, xend=3, y=bar_height, yend=bar_height - 0.02*y_range), color='black') +
  geom_segment(aes(x=4, xend=4, y=bar_height, yend=bar_height - 0.02*y_range), color='black') +
  annotate('text', x=3.5, y=bar_height + 0.03*y_range, 
           label=paste0('padj = ', round(hmgb1_tac_sham$padj, 3)), size=3.5) +
  coord_cartesian(ylim=c(min(hmgb1_expr_long$expression, na.rm=TRUE) - 0.1*y_range, 
                          bar_height + 0.1*y_range))
dev.off()
