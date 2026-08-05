library(Seurat)
library(ggplot2)
library(patchwork)

heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

# QC metrics and cell type composition, grouped by sample
group_var <- 'sample'

if (!'percent.mt' %in% colnames(heart[[]])) {
  heart[['percent.mt']] <- PercentageFeatureSet(heart, pattern = '^mt-', assay = 'RNA')
}

qc_features <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt')
qc_titles <- c('UMI count', 'Number of features', 'Mitochondrial %')

violins <- lapply(seq_along(qc_features), function(i) {
  VlnPlot(heart, features = qc_features[i], group.by = group_var, pt.size = 0) +
    ggtitle(qc_titles[i]) +
    labs(x = NULL, y = NULL) +
    NoLegend() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

prop <- as.data.frame(prop.table(
  table(heart$celltype, heart[[group_var, drop = TRUE]]), margin = 2))
colnames(prop) <- c('celltype', 'sample', 'proportion')

stacked <- ggplot(prop, aes(x = sample, y = proportion, fill = celltype)) +
  geom_col(width = 0.7) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  labs(x = NULL, y = NULL, title = 'Cell type proportion', fill = 'Cell type') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

qc_panel <- wrap_plots(c(violins, list(stacked)), ncol = 2)
ggsave('QC_celltype_panels.pdf', qc_panel, width = 11, height = 9)

##################################################################
# Quick and dirty Wilcoxon DE: 9w vs 78w and Sham vs TAC         #
# n=1 animal per condition, so cells are pseudoreplicates and    #
# p-values are strongly anti-conservative. Ranking only.         #
##################################################################
library(data.table)

DefaultAssay(heart) <- 'RNA'
if (inherits(heart[['RNA']], 'Assay5')) {
  heart[['RNA']] <- JoinLayers(heart[['RNA']])
}
heart <- NormalizeData(heart, verbose = FALSE)

Idents(heart) <- 'sample'

run_wilcox <- function(obj, ident.1, ident.2) {
  res <- FindMarkers(obj, ident.1 = ident.1, ident.2 = ident.2,
                     test.use = 'wilcox', logfc.threshold = 0, min.pct = 0.1)
  res$gene <- rownames(res)
  res$contrast <- paste0(ident.1, '_vs_', ident.2)
  res
}

de_age <- run_wilcox(heart, '78w', '9w')
de_tac <- run_wilcox(heart, 'TAC', 'Sham')

# chrX genes from the mm39 RefSeq annotation (same source as core_escape_SNPs.R)
annotation <- fread('../GRCm39/annotation_us.bed')
chrX_genes <- unique(annotation[V1 == 'chrX', V4])

top_chrX <- function(res, n = 50, padj_cutoff = 0.05) {
  res <- res[res$gene %in% chrX_genes & res$p_val_adj < padj_cutoff, ]
  res[order(-abs(res$avg_log2FC)), ][seq_len(min(n, nrow(res))), ]
}

top_age <- top_chrX(de_age)
top_tac <- top_chrX(de_tac)

write.csv(de_age, 'DE_wilcox_78w_vs_9w_all_genes.csv', row.names = FALSE)
write.csv(de_tac, 'DE_wilcox_TAC_vs_Sham_all_genes.csv', row.names = FALSE)
write.csv(top_age, 'DE_wilcox_78w_vs_9w_top_chrX.csv', row.names = FALSE)
write.csv(top_tac, 'DE_wilcox_TAC_vs_Sham_top_chrX.csv', row.names = FALSE)

print(top_age)
print(top_tac)

Idents(heart) <- 'celltype'

Col4a6
Mid1ip1

