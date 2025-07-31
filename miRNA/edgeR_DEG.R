library(edgeR)
library(fgsea)

# Combine HTseq count into one file
count_files <- list.files('HTseq_counts', pattern = "*.count", full.names = TRUE)

counts <- do.call(cbind, lapply(count_files, function(file) {
  read.delim(file, header = FALSE, row.names = 1)
}))

# Set column names to sample names
colnames(counts) <- gsub('_stranded_reverse.count', '', basename(count_files))

# Write counts file
write.csv(counts, 'combined_counts_ENSMUSG.csv', row.names = TRUE)

# 1. Remove version numbers from rownames to use as a filter for biomaRt
ensembl_ids_no_version <- gsub("\\..*$", "", rownames(counts))

# Convert from ENSMUSG to gene symbols
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_symbols <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'ensembl_gene_id',
                      values = unique(ensembl_ids_no_version),
                      mart = mart)

counts <- counts[ensembl_ids_no_version %in% gene_symbols$ensembl_gene_id, ]
ensembl_ids_no_version_mapped <- ensembl_ids_no_version[ensembl_ids_no_version %in% gene_symbols$ensembl_gene_id]
counts <- rowsum(as.matrix(counts), group = gene_symbols$external_gene_name[match(ensembl_ids_no_version_mapped, gene_symbols$ensembl_gene_id)])
# Write counts with gene symbols
write.csv(counts, 'combined_counts_gene_symbols.csv', row.names = TRUE)

# Read in metadata
metadata <- read.csv('SraRunTable.csv')
metadata <- metadata[match(metadata$Run, colnames(counts)),]
metadata <- metadata[, c('Run', 'AGE', 'genotype', 'sex', 'tissue')]

### Split data by cell type ###
metadata_split <- split(metadata, metadata$tissue)

## Run DGE analysis for each tissue type ##
dge_results <- list()

for (tissue in names(metadata_split)) {
  message(paste("Processing tissue:", tissue))
  
  tissue_metadata <- metadata_split[[tissue]]
  
  # Ensure there are at least two groups to compare for the model
  if (length(unique(tissue_metadata$sex)) < 2) {
    message(paste("Skipping", tissue, "- not enough groups for comparison."))
    next
  }
  tissue_metadata$sex <- relevel(factor(tissue_metadata$sex), ref = "male")
  y <- DGEList(counts = counts[, tissue_metadata$Run], group = tissue_metadata$sex)
  
  keep <- filterByExpr(y, group = tissue_metadata$sex)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  y <- calcNormFactors(y)
  
  design <- model.matrix(~sex, data = tissue_metadata)
  y <- estimateDisp(y, design)
  
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  dge_results[[tissue]] <- as.data.frame(topTags(qlf, n = Inf))
}

save(dge_results, file = 'dge_results.RData')

#### Enrichment of miRNA targets ####
pathways <-  gmtPathways('miRNA/m3.mirdb.v2025.1.Mm.symbols.gmt')
names(pathways) <- gsub('MIR_', 'mmu-miR-', names(pathways))

fgseaRes_list <- list()
for (tissue in names(dge_results)) {
  ranks <- sign(dge_results[[tissue]]$logFC) * -log10(dge_results[[tissue]]$FDR)
  names(ranks) <- rownames(dge_results[[tissue]])
  
  fgseaRes_list[[tissue]] <- fgsea(pathways = pathways, stats = ranks)
}

lapply(fgseaRes_list, function(x) {
  subset(x, padj < 0.05)
})

miRNA_SNP <- unlist(lapply(result, function(x){
    tmp <- gsub('Name=', '', strsplit(x$V9, ";")[[1]][3])
    gsub('-mir-', '-miR-', tmp)
}))