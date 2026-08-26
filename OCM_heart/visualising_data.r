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

##################################################################
# Hypertrophic markers in cardiomyocytes: Sham vs TAC            #
# Restricted to myocytes on purpose - these transcripts are so   #
# abundant that they show up in the ambient background of every  #
# droplet, so a per-celltype view reads as if lymphocytes were   #
# expressing them. Only the myocyte comparison is interpretable. #
##################################################################

DefaultAssay(heart) <- 'RNA'

cm_types <- c('Ventricular Cardiomyocytes', 'Cardiomyocytes (stressed)')
cm <- subset(heart, subset = celltype %in% cm_types & sample %in% c('Sham', 'TAC'))
cm$sample   <- factor(cm$sample, levels = c('Sham', 'TAC'))
cm$celltype <- droplevels(factor(cm$celltype))

# The stressed cluster is TAC-enriched, so pooling it in would fold a shift in
# cell-type composition into what should be a per-cell expression comparison.
# Violins and stats use ventricular myocytes only; the dot plot keeps both
# subtypes, where they are shown as separate rows rather than merged.
vcm <- subset(cm, subset = celltype == 'Ventricular Cardiomyocytes')

# Ctgf is Ccn2 in newer annotations - keep whichever the object carries
hyper_markers <- c('Nppa', 'Nppb', 'Acta1', 'Ankrd1', 'Myh7', 'Ctgf', 'Ccn2',
                   'Xirp2', 'Myh6', 'Atp2a2')
absent <- setdiff(hyper_markers, rownames(cm))
if (length(absent)) message('Not in object: ', paste(absent, collapse = ', '))
hyper_markers <- intersect(hyper_markers, rownames(cm))

# per-gene violins, Sham vs TAC, ventricular cardiomyocytes only
cm_vln <- VlnPlot(vcm, features = hyper_markers, group.by = 'sample',
                  pt.size = 0, ncol = 3, cols = c('#f4a582', '#2b8a8a')) &
  labs(x = NULL) &
  NoLegend()
cm_vln <- cm_vln + plot_annotation(title = 'Ventricular cardiomyocytes: Sham vs TAC')
ggsave('CM_hypertrophic_markers_violin.pdf', cm_vln,
       width = 9, height = 3 * ceiling(length(hyper_markers) / 3))

# compact summary: both myocyte subtypes x condition in one panel
cm$group <- factor(
  paste0(ifelse(as.character(cm$celltype) == 'Cardiomyocytes (stressed)',
                'Stressed CM', 'Ventricular CM'), ' - ', cm$sample),
  levels = c('Ventricular CM - Sham', 'Ventricular CM - TAC',
             'Stressed CM - Sham', 'Stressed CM - TAC'))
cm$group <- droplevels(cm$group)

cm_dot <- DotPlot(cm, features = hyper_markers, group.by = 'group') +
  scale_colour_gradient2(low = '#2166ac', mid = 'grey90', high = '#b2182b',
                         midpoint = 0, name = 'Scaled\nexpression') +
  labs(x = NULL, y = NULL, title = 'Hypertrophic markers in cardiomyocytes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('CM_hypertrophic_markers_dotplot.pdf', cm_dot, width = 8, height = 4)

# numbers behind the plots: mean log-normalised expression and % detected
expr_cm <- GetAssayData(vcm, assay = 'RNA', layer = 'data')
cm_stats <- do.call(rbind, lapply(hyper_markers, function(g) {
  v <- expr_cm[g, ]
  s <- vcm$sample == 'Sham'
  t <- vcm$sample == 'TAC'
  data.frame(gene = g,
             mean_Sham = mean(v[s]), mean_TAC = mean(v[t]),
             delta_mean = mean(v[t]) - mean(v[s]),
             pct_Sham = mean(v[s] > 0), pct_TAC = mean(v[t] > 0))
}))
cm_stats <- cm_stats[order(-cm_stats$delta_mean), ]
write.csv(cm_stats, 'CM_hypertrophic_markers_stats.csv', row.names = FALSE)

cat('\n--- Hypertrophic markers, ventricular cardiomyocytes (n = 1 animal per group) ---\n')
print(cm_stats, row.names = FALSE)
cat('\nCells per group:\n'); print(table(cm$celltype, cm$sample))

##################################################################
# Fibrosis markers in fibroblasts: Sham vs TAC                   #
# Same reasoning as the myocyte section - collagens are abundant #
# enough in a fibrotic ventricle to leak into the background of  #
# every droplet, so the per-cell comparison is restricted to the #
# cells that actually make them.                                 #
##################################################################

# activation / myofibroblast markers plus core ECM
fibro_markers <- c('Postn', 'Cthrc1', 'Meox1', 'Acta2', 'Col1a1', 'Col3a1',
                   'Col8a1', 'Fn1', 'Timp1', 'Comp', 'Thbs4', 'Ctgf', 'Ccn2')
absent_fib <- setdiff(fibro_markers, rownames(heart))
if (length(absent_fib)) message('Not in object: ', paste(absent_fib, collapse = ', '))
fibro_markers <- intersect(fibro_markers, rownames(heart))

fib <- subset(heart, subset = celltype == 'Fibroblasts' & sample %in% c('Sham', 'TAC'))
fib$sample <- factor(fib$sample, levels = c('Sham', 'TAC'))

fib_vln <- VlnPlot(fib, features = fibro_markers, group.by = 'sample',
                   pt.size = 0, ncol = 3, cols = c('#f4a582', '#2b8a8a')) &
  labs(x = NULL) &
  NoLegend()
fib_vln <- fib_vln + plot_annotation(title = 'Fibroblasts: Sham vs TAC')
ggsave('Fibrosis_markers_violin.pdf', fib_vln,
       width = 9, height = 3 * ceiling(length(fibro_markers) / 3))

# is activation fibroblast-specific, or shared with the other stromal cells?
# matched with grepl because the pericyte/epicardial labels contain a dash
# whose exact character is easy to get wrong
stroma <- subset(heart, subset = sample %in% c('Sham', 'TAC'))
stroma_ct <- as.character(stroma$celltype)
stroma <- subset(stroma, cells = colnames(stroma)[
  grepl('Fibroblast|Pericyte|Epicard', stroma_ct)])

ct <- as.character(stroma$celltype)
stroma$group <- factor(paste0(
  ifelse(grepl('Fibroblast', ct), 'Fibroblasts',
         ifelse(grepl('Pericyte', ct), 'Pericyte/SMC', 'Epicardial/Meso')),
  ' - ', stroma$sample))
stroma$group <- droplevels(stroma$group)

fib_dot <- DotPlot(stroma, features = fibro_markers, group.by = 'group') +
  scale_colour_gradient2(low = '#2166ac', mid = 'grey90', high = '#b2182b',
                         midpoint = 0, name = 'Scaled\nexpression') +
  labs(x = NULL, y = NULL, title = 'Fibrosis markers across stromal cells') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Fibrosis_markers_dotplot.pdf', fib_dot, width = 9, height = 4.5)

expr_fib <- GetAssayData(fib, assay = 'RNA', layer = 'data')
fib_stats <- do.call(rbind, lapply(fibro_markers, function(g) {
  v <- expr_fib[g, ]
  s <- fib$sample == 'Sham'
  t <- fib$sample == 'TAC'
  data.frame(gene = g,
             mean_Sham = mean(v[s]), mean_TAC = mean(v[t]),
             delta_mean = mean(v[t]) - mean(v[s]),
             pct_Sham = mean(v[s] > 0), pct_TAC = mean(v[t] > 0))
}))
fib_stats <- fib_stats[order(-fib_stats$delta_mean), ]
write.csv(fib_stats, 'Fibrosis_markers_stats.csv', row.names = FALSE)

cat('\n--- Fibrosis markers, fibroblasts (n = 1 animal per group) ---\n')
print(fib_stats, row.names = FALSE)

# fibrosis also shows up as fibroblast expansion, not only as expression
st <- heart$sample %in% c('Sham', 'TAC')
comp <- prop.table(table(heart$celltype[st],
                         factor(heart$sample[st], c('Sham', 'TAC'))), margin = 2)
cat('\nCell type composition (% of cells per sample):\n')
print(round(comp * 100, 2))
write.csv(as.data.frame.matrix(round(comp * 100, 2)),
          'Celltype_composition_Sham_TAC.csv')
