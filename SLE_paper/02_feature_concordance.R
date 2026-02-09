# 02_feature_concordance.R
#
# Analyses feature selection consistency across cross-validation splits,
# generates per-celltype heatmaps, and runs cell composition (propeller)
# analysis.
#
# Requires:
#   - ../results_update/top_celltypes.RData  (from 01_model_comparison.R)
#
# Outputs:
#   - ../results_update/selected_features.RData
#   - ../results_update/all_features.csv
#   - ../results_update/top_features.csv
#   - ../results_update/top_celltypes_all_features.csv
#   - ../results_update/propellor_results.RData
#   - Per-celltype feature heatmap PDFs

library(dplyr)
library(ggplot2)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(speckle)

# ---------- Load shared data ----------
load('../escapees.Rdata')
escape <- rownames(escape)
chrX <- read.delim('../chrX_biomaRt.txt')$Gene.name
chrX <- unique(chrX[!chrX %in% ''])

load('../results_update/top_celltypes.RData')

# ============================================================
# SECTION 1 — Feature concordance across CV splits
# ============================================================

feature_list <- list()
for (i in 1:10) {
  feature.files <- list.files(paste0('split_', i, '/features'),
                               pattern = 'combined', full.names = TRUE)
  features <- lapply(feature.files, function(x) read.csv(x)$Feature)
  names(features) <- gsub('combined_|.csv', '', basename(feature.files))
  feature_list[[i]] <- features
}
celltype.geneset <- names(feature_list[[1]])
# Subset for chrX and SLE
celltype.geneset <- celltype.geneset[grep('.chrX|.SLE', celltype.geneset)]

result_list <- list()
for (i in celltype.geneset) {
  tmp <- lapply(1:10, function(x) feature_list[[x]][[i]])
  result_list[[i]] <- tmp
}
names(result_list) <- celltype.geneset

if (!dir.exists('../results_update/feature_heatmap')) {
  dir.create('../results_update/feature_heatmap')
}

# --- Plot heatmap of selected features per celltype/geneset ---
for (file in celltype.geneset) {
  mtx <- fromList(result_list[[file]])
  rownames(mtx) <- unique(unlist(result_list[[file]]))
  colnames(mtx) <- paste('split', 1:10, sep = '_')
  mtx <- mtx[order(rowSums(mtx), decreasing = TRUE), ]

  col_fun <- colorRamp2(c(0, 1), c("white", "black"))
  pdf(paste0('../results_update/feature_heatmap/', file, '.pdf'))
  p <- Heatmap(as.matrix(mtx),
               name = 'Features',
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               col = col_fun,
               row_names_gp = gpar(fontsize = 5),
               show_heatmap_legend = FALSE)
  print(p)
  dev.off()
}

all_features <- lapply(result_list, function(x) unique(unlist(x)))
top_features <- lapply(result_list, function(x) {
  tmp <- table(unlist(x))
  top <- names(tmp[tmp == 10])
  if (length(top) == 0) list() else top
})

# ============================================================
# SECTION 2 — Save selected features
# ============================================================

selected_features <- list(all_features = all_features, top_features = top_features)
save(selected_features, file = '../results_update/selected_features.RData')

# Write as CSV
all_features_mtx <- bind_rows(lapply(names(selected_features$all_features), function(x) {
  data.frame(celltype = x, feature = selected_features$all_features[[x]])
}))
write.csv(all_features_mtx, '../results_update/all_features.csv', row.names = FALSE)

selected_features$top_features[['B_cell_Naive.chrX']]    <- "NA"
selected_features$top_features[['Conventional_DC.chrX']]  <- "NA"

top_features_mtx <- bind_rows(lapply(names(selected_features$top_features), function(x) {
  data.frame(celltype = x, feature = selected_features$top_features[[x]])
}))
write.csv(top_features_mtx, '../results_update/top_features.csv', row.names = FALSE)

# ============================================================
# SECTION 3 — Top cell type feature heatmaps
# ============================================================

# Edit celltype names to match filenames
top_celltypes <- gsub("\\+|-| ", "_", top_celltypes)

# -- Top features (selected in all 10 splits) --
top_celltypes_top_features <- selected_features$top_features[
  names(selected_features$top_features) %in% paste0(gsub(' ', '_', top_celltypes), '.chrX')
]
top_celltypes_top_features_mtx <- fromList(top_celltypes_top_features)
rownames(top_celltypes_top_features_mtx) <- unique(unlist(top_celltypes_top_features))
colnames(top_celltypes_top_features_mtx) <- gsub('.chrX', '',
                                                  names(top_celltypes_top_features_mtx)) %>%
  gsub('_', ' ', .)

pdf('../results_update/top_celltypes_top_features_heatmap.pdf')
col_fun <- colorRamp2(c(0, 1), c("white", "red"))
p <- Heatmap(as.matrix(top_celltypes_top_features_mtx),
             name = 'Features', col = col_fun, show_heatmap_legend = FALSE)
print(p)
dev.off()

# -- All features (selected in any split) --
top_celltypes_all_features <- selected_features$all_features[
  names(selected_features$all_features) %in% paste0(gsub(' ', '_', top_celltypes), '.chrX')
]
top_celltypes_all_features_mtx <- fromList(top_celltypes_all_features)
rownames(top_celltypes_all_features_mtx) <- unique(unlist(top_celltypes_all_features))
colnames(top_celltypes_all_features_mtx) <- gsub('.chrX', '',
                                                  names(top_celltypes_all_features_mtx)) %>%
  gsub('_', ' ', .)

pdf('../results_update/top_celltypes_all_features_heatmap.pdf', width = 10, height = 10)
col_fun <- colorRamp2(c(0, 1), c("white", "red"))
p <- Heatmap(as.matrix(top_celltypes_all_features_mtx),
             name = 'Features', col = col_fun, show_heatmap_legend = FALSE)
print(p)
dev.off()

lapply(names(top_celltypes_all_features), function(x) {
  data.frame(celltype = x, feature = top_celltypes_all_features[[x]])
}) %>%
  bind_rows() %>%
  write.csv('../results_update/top_celltypes_all_features.csv', row.names = FALSE)

# Euclidean clustering on the top celltypes feature matrix
clustering <- hclust(dist(t(top_celltypes_all_features_mtx)), method = 'complete')
groups <- cutree(clustering, k = 3)

# ============================================================
# SECTION 4 — Cell composition analysis (propeller)
# ============================================================

metadata <- read.delim('../metadata.tsv')

index <- c(read.csv('split_1/train_index.csv')$rownames,
           read.csv('split_1/test_index.csv')$rownames)
metadata <- metadata[metadata$ind_cov %in% index, ]
metadata$disease <- ifelse(metadata$disease == 'systemic lupus erythematosus',
                           'disease', 'control')
metadata$disease <- factor(metadata$disease, levels = c('disease', 'control'))

# Create detailed cell type column
metadata$cell_type_detailed <- case_when(
  metadata$cell_type == "CD4-positive, alpha-beta T cell" & metadata$ct_cov == "T4_naive" ~ "CD4 T cell Naive",
  metadata$cell_type == "CD4-positive, alpha-beta T cell" & metadata$ct_cov == "T4_em"    ~ "CD4 T cell Effector Memory",
  metadata$cell_type == "CD4-positive, alpha-beta T cell" & metadata$ct_cov == "T4_reg"   ~ "CD4 T cell Treg",
  metadata$cell_type == "CD8-positive, alpha-beta T cell" & metadata$ct_cov == "T8_naive" ~ "CD8 T cell Naive",
  metadata$cell_type == "CD8-positive, alpha-beta T cell" & metadata$ct_cov == "CytoT_GZMH+" ~ "CD8 T cell Cytotoxic GZMH+",
  metadata$cell_type == "CD8-positive, alpha-beta T cell" & metadata$ct_cov == "CytoT_GZMK+" ~ "CD8 T cell Cytotoxic GZMK+",
  metadata$cell_type == "CD8-positive, alpha-beta T cell" & metadata$ct_cov == "T_mait"   ~ "CD8 T cell MAIT",
  metadata$cell_type == "natural killer cell" & metadata$ct_cov == "NK_bright" ~ "NK cell Bright",
  metadata$cell_type == "natural killer cell" & metadata$ct_cov == "NK_dim"    ~ "NK cell Dim",
  metadata$cell_type == "B cell" & metadata$ct_cov == "B_naive"    ~ "B cell Naive",
  metadata$cell_type == "B cell" & metadata$ct_cov == "B_mem"      ~ "B cell Memory",
  metadata$cell_type == "B cell" & metadata$ct_cov == "B_plasma"   ~ "B cell Plasma",
  metadata$cell_type == "B cell" & metadata$ct_cov == "B_atypical" ~ "B cell Atypical",
  metadata$cell_type == "progenitor cell"            ~ "Progenitor cell",
  metadata$cell_type == "conventional dendritic cell" ~ "Conventional DC",
  metadata$cell_type == "plasmacytoid dendritic cell" ~ "Plasmacytoid DC",
  metadata$cell_type == "classical monocyte"          ~ "Classical monocyte",
  metadata$cell_type == "non-classical monocyte"      ~ "Non-classical monocyte",
  metadata$cell_type == "lymphocyte"                  ~ "Lymphocyte",
  metadata$cell_type == "plasmablast"                 ~ "Plasmablast",
  TRUE ~ metadata$cell_type
)

output.asin <- propeller(clusters = metadata$cell_type_detailed,
                          sample = metadata$ind_cov,
                          group = metadata$disease,
                          transform = 'asin')
output.asin$BaselineProp.clusters <- gsub("-| ", "_", output.asin$BaselineProp.clusters)
foo <- subset(output.asin,
              BaselineProp.clusters %in% top_celltypes & FDR < 0.05
              )[, c('BaselineProp.clusters', 'Tstatistic', 'FDR')]
foo$Tstatistic <- round(foo$Tstatistic, 2)
foo$FDR        <- round(foo$FDR, 6)

save(output.asin, file = '../results_update/propellor_results.RData')

# ============================================================
# SECTION 5 — Escape gene enrichment test
# ============================================================

hits <- rownames(top_celltypes_all_features_mtx)[
  rownames(top_celltypes_all_features_mtx) %in% escape
]

selected <- selected_features$all_features[["CD8_positive,_alpha_beta_T_cell.chrX"]]

a <- sum(selected %in% escape)
b <- sum(!(selected %in% escape))
c <- sum(chrX[!chrX %in% selected] %in% escape)
d <- sum(!(chrX[!chrX %in% selected] %in% escape))
chisq.test(matrix(c(a, b, c, d), nrow = 2))

# --- edgeR DEGs for CD4/CD8 T cell features ---
source('../edgeR.list.R')
edgeR <- deg.list('../edgeR', filter = FALSE)

top_celltypes_all_features <- selected_features$all_features[c(
  'CD4_positive,_alpha_beta_T_cell.chrX',
  'CD8_positive,_alpha_beta_T_cell.chrX'
)]
edgeR_top_celltypes_all_features <- lapply(names(top_celltypes_all_features), function(x) {
  tmp <- edgeR[[gsub('.chrX', '', x)]]
  subset(tmp, gene %in% top_celltypes_all_features[[x]])[, c('gene', 'logFC', 'FDR')]
})
names(edgeR_top_celltypes_all_features) <- names(top_celltypes_all_features)
