# 04_chrX_feature_analysis.R
#
# In-depth analysis of X chromosome features: escape gene annotation,
# DEG overlap, Jaccard similarity between gene sets, consistently
# selected features, and XIST/TSIX correlation.
#
# Requires:
#   - ../results_update/selected_features.RData  (from 02_feature_concordance.R)
#   - ../results_update/top_celltypes.RData       (from 01_model_comparison.R)
#   - ../edgeR.list.R  (helper function)
#
# Outputs:
#   - Supplementary_Table_2.csv
#   - Multiple heatmap and volcano PDFs in ../results_update/

library(dplyr)
library(ggplot2)
library(reshape2)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(ggrepel)

# ---------- Load shared data ----------
load('../results_update/selected_features.RData')
load('../results_update/top_celltypes.RData')
SLE <- read.delim('../SLE_DisGeNet.tsv')$Gene
load('../escapees.Rdata')
escape <- rownames(escape)

source('../edgeR.list.R')
edgeR <- deg.list('../edgeR', filter = FALSE)

# ============================================================
# SECTION 1 — chrX feature heatmap with escape gene annotation
# ============================================================

chrX.list <- selected_features$all_features[grep('.chrX',
  names(selected_features$all_features), value = TRUE)]
names(chrX.list) <- gsub('.chrX', '', names(chrX.list))

names(chrX.list)[9] <- "CD8 positive, alpha-beta T cell"
names(chrX.list)[5] <- "CD4 positive, alpha-beta T cell"
names(chrX.list) <- gsub('_', ' ', names(chrX.list))

top_list <- chrX.list[top_celltypes]
names(top_list)[1] <- "CD4 positive, alpha-beta T cell"
names(top_list)[2] <- "CD8 positive, alpha-beta T cell"
names(top_list)[3] <- "Progenitor cell"

mat <- fromList(top_list)
rownames(mat) <- unique(unlist(top_list))

annotation <- ifelse(rownames(mat) %in% escape, 'Escapee', 'Non-escapee')

mat        <- mat[order(annotation, decreasing = FALSE), ]
annotation <- annotation[order(annotation, decreasing = FALSE)]

row_ha <- rowAnnotation(
  Feature_type = annotation,
  col = list(Feature_type = c("Escapee" = "#2284F5", "Non-escapee" = "#F56022")),
  show_annotation_name = FALSE
)

pdf('../results_update/heatmap_chrX_features_clusters.pdf', width = 10, height = 10)
col_fun <- colorRamp2(c(0, 1), c("white", "red"))
p <- Heatmap(as.matrix(mat),
             name = 'Features', col = col_fun,
             show_heatmap_legend = FALSE,
             cluster_rows = FALSE,
             right_annotation = row_ha)
print(p)
dev.off()

# ============================================================
# SECTION 2 — Jaccard index between gene sets
# ============================================================

jaccard_index <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}

if (!dir.exists('../results_update/jaccard_heatmaps')) {
  dir.create('../results_update/jaccard_heatmaps')
}

celltypes <- unique(gsub('.HVG.autosome|.HVG|.autosome|.SLE|.chrX', '',
                          names(selected_features$all_features)))

for (celltype in celltypes) {
  tmp <- selected_features$all_features[grep(paste0('^', celltype),
                                              names(selected_features$all_features))]
  combinations <- combn(names(tmp), 2)
  mtx <- matrix(0, ncol = 5, nrow = 5, dimnames = list(names(tmp), names(tmp)))
  for (i in 1:ncol(combinations)) {
    mtx[combinations[1, i], combinations[2, i]] <-
      jaccard_index(tmp[[combinations[1, i]]], tmp[[combinations[2, i]]])
    mtx[combinations[2, i], combinations[1, i]] <-
      mtx[combinations[1, i], combinations[2, i]]
  }
  diag(mtx) <- 1
  colnames(mtx) <- gsub(celltype, '', colnames(mtx)) %>% gsub('^\\.', '', .)
  rownames(mtx) <- gsub(celltype, '', rownames(mtx)) %>% gsub('^\\.', '', .)
  pdf(paste0('../results_update/jaccard_heatmaps/', celltype, '.pdf'))
  col <- circlize::colorRamp2(c(0, 1), c("white", "red"))
  print(Heatmap(mtx, col = col, name = 'Jaccard Index',
                column_title = replace.names('CD16+.NK.cells'),
                cluster_columns = FALSE, cluster_rows = FALSE))
  dev.off()
}

# --- chrX vs HVG/SLE Jaccard matrix ---
result_list <- list()
for (celltype in celltypes) {
  result <- lapply(c('.HVG', '.SLE'), function(x) {
    jaccard_index(selected_features$all_features[[paste0(celltype, x)]],
                  selected_features$all_features[[paste0(celltype, '.chrX')]])
  })
  names(result) <- c('HVG', 'SLE')
  tmp <- do.call(rbind, result)
  colnames(tmp) <- celltype
  result_list[[celltype]] <- tmp
}
jaccard_matrix <- do.call(cbind, result_list)
colnames(jaccard_matrix) <- gsub('_', ' ', colnames(jaccard_matrix))

pdf('../results_update/chrX.geneset.jaccard.heatmap.pdf')
col <- circlize::colorRamp2(c(0, 1), c("white", "red"))
Heatmap(t(jaccard_matrix), col = col, name = 'Jaccard Index',
        cluster_columns = FALSE, cluster_rows = FALSE)
dev.off()

# --- Overlap statistics ---
overlap_list <- list()
for (celltype in celltypes) {
  X_genes     <- selected_features$all_features[[paste0(celltype, '.chrX')]]
  num_chrX    <- length(X_genes)
  HVG_overlap <- length(intersect(selected_features$all_features[[paste0(celltype, '.HVG')]],
                                   X_genes))
  HVG_Jaccard <- jaccard_index(selected_features$all_features[[paste0(celltype, '.HVG')]],
                                X_genes)
  SLE_overlap <- length(intersect(selected_features$all_features[[paste0(celltype, '.SLE')]],
                                   X_genes))
  SLE_Jaccard <- jaccard_index(selected_features$all_features[[paste0(celltype, '.SLE')]],
                                X_genes)
  DEG_overlap <- length(intersect(X_genes,
                   subset(edgeR[[celltype]], abs(logFC) > 0.1 & FDR < 0.05)$gene))
  DEG_jaccard <- jaccard_index(X_genes,
                   subset(edgeR[[celltype]], abs(logFC) > 0.1 & FDR < 0.05)$gene)
  overlap_list[[celltype]] <- data.frame(
    celltype = gsub('_', ' ', celltype), num_chrX = num_chrX,
    HVG_overlap = HVG_overlap, HVG_Jaccard = round(HVG_Jaccard, 4),
    SLE_overlap = SLE_overlap, SLE_Jaccard = round(SLE_Jaccard, 4),
    DEG_overlap = DEG_overlap, DEG_jaccard = round(DEG_jaccard, 4),
    row.names = NULL
  )
}
overlap_df <- bind_rows(overlap_list)
write.csv(overlap_df, '../results_update/chrX_model_overlap_df.csv', row.names = FALSE)

# --- Intersecting features across gene sets ---
intersect_list <- list()
for (celltype in celltypes) {
  result <- unlist(lapply(c('.HVG', '.SLE'), function(x) {
    intersect(selected_features$all_features[[paste0(celltype, x)]],
              selected_features$all_features[[paste0(celltype, '.chrX')]])
  }))
  intersect_list[[celltype]] <- result
}
sort(table(unlist(intersect_list)))

# ============================================================
# SECTION 3 — DEG heatmaps for selected features
# ============================================================

deg_features_list <- list()
for (i in seq_along(edgeR)) {
  tmp <- lapply(grep(names(edgeR)[i], names(selected_features$all_features), value = TRUE),
                function(x) subset(edgeR[[i]], gene %in% selected_features$all_features[[x]]))
  names(tmp) <- names(selected_features$all_features)[
    grep(names(edgeR)[i], names(selected_features$all_features))
  ]
  deg_features_list <- c(deg_features_list, tmp)
}

deg_top_features <- deg_features_list[
  names(deg_features_list) %in% paste0(gsub(' ', '_', top_celltypes), '.chrX')
]

# Heatmap of cell type × gene coloured by logFC
combined_deg_top_features <- bind_rows(deg_top_features, .id = 'celltype')
combined_deg_top_features <- combined_deg_top_features[-2, ]
degs_mtx <- reshape2::dcast(combined_deg_top_features, gene ~ celltype, value.var = 'logFC')
rownames(degs_mtx) <- degs_mtx$gene
degs_mtx <- degs_mtx[, -1]
degs_mtx[is.na(degs_mtx)] <- 0
colnames(degs_mtx) <- replace.names(colnames(degs_mtx))

pdf('../results_update/chrX_DEG_heatmap.pdf')
col <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
Heatmap(as.matrix(degs_mtx), col = col, name = 'logFC',
        cluster_columns = TRUE, cluster_rows = TRUE)
dev.off()

# --- Shared features DEG heatmap ---
degs <- lapply(celltypes, function(x) subset(edgeR[[x]], gene %in% intersect_list[[x]]))
names(degs) <- celltypes

combined_degs <- bind_rows(degs, .id = 'celltype')
degs_mtx <- reshape2::dcast(combined_degs, gene ~ celltype, value.var = 'logFC')
rownames(degs_mtx) <- degs_mtx$gene
degs_mtx <- degs_mtx[, -1]
degs_mtx[is.na(degs_mtx)] <- 0
colnames(degs_mtx) <- replace.names(colnames(degs_mtx))

fdr_matrix <- reshape2::dcast(combined_degs, gene ~ celltype, value.var = 'FDR')
rownames(fdr_matrix) <- fdr_matrix$gene
colnames(fdr_matrix) <- replace.names(colnames(fdr_matrix))
fdr_matrix <- fdr_matrix[, -1]
significance_matrix <- ifelse(fdr_matrix < 0.05, "*", "ns")
significance_matrix <- ifelse(fdr_matrix < 0.01, "**", significance_matrix)
significance_matrix <- ifelse(fdr_matrix < 0.001, "***", significance_matrix)
significance_matrix[is.na(significance_matrix)] <- ""

pdf('results_update/shared_features_deg_heatmap.pdf')
col <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
Heatmap(as.matrix(degs_mtx), col = col, name = 'logFC',
        cluster_columns = TRUE, cluster_rows = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(significance_matrix[i, j], x, y, gp = gpar(fontsize = 10))
        })
dev.off()

# ============================================================
# SECTION 4 — chrX Jaccard index between top cell types
# ============================================================

chrX_features <- selected_features$all_features[grep('.chrX',
                   names(selected_features$all_features))]
names(chrX_features) <- gsub('.chrX', '', names(chrX_features))

top_celltypes <- c('CD16+.NK.cells', 'Memory.B.cells', 'Tcm.Naive.helper.T.cells',
                   'Regulatory.T.cells', 'DC1', 'DC2', 'pDC',
                   'Tem.Effector.helper.T.cells', 'Non.classical.monocytes')
chrX_features <- chrX_features[top_celltypes]

sort(table(unlist(chrX_features)))

combinations <- combn(names(chrX_features), 2)
jaccard_mtx <- matrix(0, ncol = 9, nrow = 9,
                      dimnames = list(names(chrX_features), names(chrX_features)))
for (i in 1:ncol(combinations)) {
  jaccard_mtx[combinations[1, i], combinations[2, i]] <-
    jaccard_index(chrX_features[[combinations[1, i]]],
                  chrX_features[[combinations[2, i]]])
  jaccard_mtx[combinations[2, i], combinations[1, i]] <-
    jaccard_index(chrX_features[[combinations[1, i]]],
                  chrX_features[[combinations[2, i]]])
}
diag(jaccard_mtx) <- 1

colnames(jaccard_mtx) <- replace.names(colnames(jaccard_mtx))
rownames(jaccard_mtx) <- replace.names(rownames(jaccard_mtx))
jaccard_mtx[lower.tri(jaccard_mtx)] <- 0

pdf('results_update/chrX_jaccard_heatmap.pdf')
col <- circlize::colorRamp2(c(0, 1), c("white", "red"))
Heatmap(jaccard_mtx, col = col, name = 'Jaccard Index')
dev.off()

# --- DEGs for top chrX features: Supplementary Table 2 ---
degs <- lapply(top_celltypes, function(x) subset(edgeR[[x]], gene %in% chrX_features[[x]]))
names(degs) <- top_celltypes

combined_degs <- bind_rows(degs, .id = 'celltype')
combined_degs$is_escape <- ifelse(combined_degs$gene %in% rownames(escape), 'True', 'False')
write.csv(combined_degs, 'results_update/Supplementary_Table_2.csv', row.names = FALSE)

degs_mtx <- reshape2::dcast(combined_degs, gene ~ celltype, value.var = 'logFC')
rownames(degs_mtx) <- degs_mtx$gene
degs_mtx <- degs_mtx[, -1]
degs_mtx[is.na(degs_mtx)] <- 0
colnames(degs_mtx) <- replace.names(colnames(degs_mtx))

fdr_matrix <- reshape2::dcast(combined_degs, gene ~ celltype, value.var = 'FDR')
rownames(fdr_matrix) <- fdr_matrix$gene
colnames(fdr_matrix) <- replace.names(colnames(fdr_matrix))
fdr_matrix <- fdr_matrix[, -1]
significance_matrix <- ifelse(fdr_matrix < 0.05, "*", "ns")
significance_matrix <- ifelse(fdr_matrix < 0.01, "**", significance_matrix)
significance_matrix <- ifelse(fdr_matrix < 0.001, "***", significance_matrix)
significance_matrix[is.na(significance_matrix)] <- ""

top_genes <- names(sort(table(combined_degs$gene), decreasing = TRUE))[1:20]

pdf('results_update/chrX_heatmap.pdf')
col <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
ann <- rowAnnotation(foo = anno_mark(
  at = which(rownames(degs_mtx) %in% top_genes),
  labels = rownames(degs_mtx)[rownames(degs_mtx) %in% top_genes]
))
Heatmap(as.matrix(degs_mtx), col = col, name = 'logFC',
        cluster_columns = TRUE, cluster_rows = TRUE,
        show_row_names = FALSE, right_annotation = ann)
dev.off()

# ============================================================
# SECTION 5 — Potential escapee volcano plots
# ============================================================

potential_escapees <- subset(combined_degs, abs(logFC) >= 0.5 & FDR < 0.05)
potential_escapees <- potential_escapees[, c('celltype', 'gene', 'logFC', 'FDR')]
unique(potential_escapees$gene)

combined_degs$celltype <- replace.names(combined_degs$celltype)

pdf('../results_update/potential_escapees.pdf')
ggplot(combined_degs, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(abs(logFC) > 0.5 & FDR < 0.05, 'red', 'grey'))) +
  scale_color_identity() +
  geom_text_repel(data = subset(combined_degs, abs(logFC) > 0.5 & FDR < 0.05),
                  aes(label = gene), color = 'black', max.overlaps = Inf, size = 2) +
  theme_minimal() +
  theme(legend.position = 'none') +
  labs(x = 'logFC', y = '-log10(FDR)') +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dotted', color = 'black') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', color = 'black') +
  facet_wrap(~celltype, ncol = 4)
dev.off()

# --- Escape gene heatmap ---
combined_degs_escape <- subset(combined_degs, gene %in% rownames(escape))
escape_mtx <- reshape2::dcast(combined_degs_escape, gene ~ celltype, value.var = 'logFC')
rownames(escape_mtx) <- escape_mtx$gene
escape_mtx <- escape_mtx[, -1]
escape_mtx[is.na(escape_mtx)] <- 0
colnames(escape_mtx) <- replace.names(colnames(escape_mtx))

escape_fdr_matrix <- reshape2::dcast(combined_degs_escape, gene ~ celltype, value.var = 'FDR')
rownames(escape_fdr_matrix) <- escape_fdr_matrix$gene
colnames(escape_fdr_matrix) <- replace.names(colnames(escape_fdr_matrix))
escape_fdr_matrix <- escape_fdr_matrix[, -1]
escape_significance_matrix <- ifelse(escape_fdr_matrix < 0.05, "*", "ns")
escape_significance_matrix <- ifelse(escape_fdr_matrix < 0.01, "**", escape_significance_matrix)
escape_significance_matrix <- ifelse(escape_fdr_matrix < 0.001, "***", escape_significance_matrix)
escape_significance_matrix[is.na(escape_significance_matrix)] <- ""

pdf('results_update/escape_heatmap.pdf')
col <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
Heatmap(as.matrix(escape_mtx), col = col, name = 'logFC',
        cluster_columns = TRUE, cluster_rows = TRUE,
        show_row_names = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(escape_significance_matrix[i, j], x, y, gp = gpar(fontsize = 2))
        })
dev.off()

subset(edgeR[['Regulatory.T.cells']],
       gene %in% c('PDCD1', 'CTLA4', 'TIM3', 'LAG3', 'BTLA', 'TIGIT')
       )[, c('gene', 'logFC', 'FDR')]

# ============================================================
# SECTION 6 — Consistently selected features (all 10 splits)
# ============================================================

chrX_features <- selected_features$top_features[grep('.chrX',
                   names(selected_features$top_features))]
names(chrX_features) <- gsub('.chrX', '', names(chrX_features))
chrX_features <- chrX_features[top_celltypes]

degs <- lapply(top_celltypes, function(x) {
  tmp <- merge(edgeR[[x]], data.frame('gene' = chrX_features[[x]]), by = 'gene', all.y = TRUE)
  return(tmp[, c('gene', 'logFC', 'FDR')])
})
names(degs) <- top_celltypes

combined_degs <- dplyr::bind_rows(degs, .id = 'celltype')
combined_degs$celltype <- replace.names(combined_degs$celltype)
write.csv(combined_degs, 'results_update/top_chrX.consistent.csv', row.names = FALSE)

plots_list <- list()
for (cell in names(chrX_features)) {
  features <- chrX_features[[cell]]
  mtx <- readRDS(paste0(cell, '.chrX.RDS'))
  mtx <- mtx[, c('class', features)]
  class <- mtx$class
  mtx_scaled <- scale(as.matrix(mtx[, -1]))
  rownames(mtx_scaled) <- class
  mtx_melt <- reshape2::melt(mtx_scaled)
  mtx_melt$Var1 <- factor(mtx_melt$Var1, levels = c('control', 'disease'))

  plot <- ggplot(mtx_melt, aes(x = Var2, y = value, color = Var1)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(), alpha = 0.5) +
    theme_minimal() +
    labs(x = '', y = 'z-score', color = 'Condition') +
    ggtitle(replace.names(cell)) +
    theme(plot.title = element_text(size = 18)) +
    theme(axis.text.x = element_text(size = 12)) +
    scale_color_manual(values = c('blue', 'red'))
  plots_list[[cell]] <- plot
}

pdf('results_update/consistent_features_boxplot.pdf', width = 15, height = 20)
grid.arrange(grobs = plots_list, ncol = 2, nrow = 5)
dev.off()

# ============================================================
# SECTION 7 — XIST / TSIX correlation in CD16+ NK cells
# ============================================================

cell <- 'CD16+.NK.cells'
mtx <- readRDS(paste0(cell, '.chrX.RDS'))
mtx_scaled <- data.frame(scale(mtx[, c('XIST', 'TSIX')]))
mtx_scaled$class <- mtx$class

coorelation_result <- lapply(split(mtx_scaled, mtx_scaled$class), function(x) {
  cor.test(x$XIST, x$TSIX, method = 'spearman')
})
cor_df <- data.frame(
  class     = names(coorelation_result),
  cor_value = sapply(coorelation_result, function(x) x$estimate),
  p_value   = sapply(coorelation_result, function(x) x$p.value)
)

pdf('results_update/CD16NK_XIST_TSIX.pdf')
ggplot(data.frame(mtx_scaled), aes(x = XIST, y = TSIX)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_minimal() +
  labs(x = 'XIST', y = 'TSIX') +
  facet_wrap(~class) +
  geom_text(data = cor_df, aes(x = Inf, y = Inf, label = paste0(
    'r = ', round(cor_value, 2),
    ', p = ', round(p_value, 2)
  )),
  hjust = 1.1, vjust = 1.1, size = 3, color = 'black')
dev.off()
