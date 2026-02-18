library(Seurat)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
library(scDEED)
options(future.globals.maxSize = 8 * 1024^3)

adult.file <- '9w/outs/cellbender_output_seurat.h5'
data.data <- Read10X_h5(filename = adult.file, use.names = TRUE)

adult <- CreateSeuratObject(counts = data.data)

# Preprocessing and Normalisation
adult <- PercentageFeatureSet(adult, pattern = "^MT-", col.name = "percent.mt")
adult <- subset(adult, subset = nCount_RNA > 500 & nCount_RNA < 20000 & nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
adult <- subset(adult, features = rownames(adult)[Matrix::rowSums(GetAssayData(adult, layer = "counts") > 0) > 10])
adult <- SCTransform(adult, vars.to.regress = "percent.mt", verbose = FALSE)

# Dimension reduction
adult <- RunPCA(adult)

pdf('9w/elbow_plot.pdf')
ElbowPlot(adult)
dev.off()

adult <- RunUMAP(adult, dims = 1:10)

# Doublet detection
sweep.res.list_adult <- paramSweep(adult, PCs = 1:10, sct = TRUE)
sweep.stats_adult <- summarizeSweep(sweep.res.list_adult, GT = FALSE)
bcmvn_adult <- find.pK(sweep.stats_adult)

    bcmvn_adult <- bcmvn_adult %>%
    mutate(pK_num = as.numeric(as.character(pK)))

pdf('9w/pK.pdf')
ggplot(bcmvn_adult, aes(x = pK_num, y = BCmetric)) + 
geom_line() +
geom_point() +
labs(x = "pK", y = "BCmvn (BCmetric)", title = "DoubletFinder pK sweep") +
theme_minimal()
dev.off()

## pK = 0.05
annotations <- adult@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*nrow(adult@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
adult <- doubletFinder(adult, PCs = 1:10, pN = 0.25, pK = 0.05, nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)

# Remove doublets
adult <- subset(adult, subset = DF.classifications_0.25_0.05_865 == "Singlet")

cell_index <- as.data.frame(colnames(adult))
rownames(cell_index) <- cell_index[,1]
write.table(cell_index, file = '9w/cell_index.txt', col.names = FALSE, row.names = TRUE, quote = FALSE)

saveRDS(adult, file = '9w/adult.RDS')

###################
aged.file <- '78w/outs/cellbender_output_seurat.h5'
aged.data <- Read10X_h5(filename = aged.file, use.names = TRUE)

aged <- CreateSeuratObject(counts = aged.data)

# Preprocessing and normalisation
aged <- PercentageFeatureSet(aged, pattern = "^MT-", col.name = "percent.mt")
aged <- subset(aged, subset = nCount_RNA > 500 & nCount_RNA < 20000 & nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
aged <- subset(aged, features = rownames(aged)[Matrix::rowSums(GetAssayData(aged, layer = "counts") > 0) > 10])
aged <- SCTransform(aged, vars.to.regress = "percent.mt", verbose = FALSE)

# Dimension reduction
aged <- RunPCA(aged)

pdf('78w/elbow_plot.pdf')
ElbowPlot(aged)
dev.off()

aged <- RunUMAP(aged, dims = 1:10)

# Doublet detection
sweep.res.list_aged <- paramSweep(aged, PCs = 1:10, sct = TRUE)
sweep.stats_aged <- summarizeSweep(sweep.res.list_aged, GT = FALSE)
bcmvn_aged <- find.pK(sweep.stats_aged)

bcmvn_aged <- bcmvn_aged %>%
  mutate(pK_num = as.numeric(as.character(pK)))

pdf('78w/pK.pdf')
ggplot(bcmvn_aged, aes(x = pK_num, y = BCmetric)) + 
  geom_line() +
  geom_point() +
  labs(x = "pK", y = "BCmvn (BCmetric)", title = "DoubletFinder pK sweep") +
  theme_minimal()
dev.off()

## pK = 0.19
annotations <- aged@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*nrow(aged@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
aged <- doubletFinder(aged, PCs = 1:10, pN = 0.25, pK = 0.19, nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)

# Remove doublets
aged <- subset(aged, subset = DF.classifications_0.25_0.19_804 == "Singlet")

cell_index <- as.data.frame(colnames(aged))
rownames(cell_index) <- cell_index[,1]
write.table(cell_index, file = '78w/cell_index.txt', col.names = FALSE, row.names = TRUE, quote = FALSE)

saveRDS(aged, file = '78w/aged.RDS')

# Integration - https://satijalab.org/seurat/articles/seurat5_integration
adult <- readRDS('9w/adult.RDS')
aged <- readRDS('78w/aged.RDS')

adult$age <- 'adult'
aged$age  <- 'aged'

merged <- merge(adult, aged, add.cell.ids = c("adult","aged"), project = "combined")
DefaultAssay(merged) <- "RNA"
merged <- SCTransform(merged, verbose = FALSE)

merged <- RunPCA(merged)
merged[["SCT"]] <- split(merged[["SCT"]], f = merged$age)
merged <- IntegrateLayers(object = merged, method = RPCAIntegration, assay='SCT')
merged <- JoinLayers(merged)


pdf('elbow_plot.pdf')
ElbowPlot(merged, reduction = "pca", ndims=30)
dev.off()


merged <- FindNeighbors(merged, reduction = "integrated.dr", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.1, cluster.name = "seurat_clusters")
Idents(merged) <- 'seurat_clusters'

merged <- RunUMAP(merged, reduction = "integrated.dr", dims = 1:30)


# ### scDEED for integrated data - https://github.com/JSB-UCLA/scDEED
# K=
# pre_embedding = 'integrated.dr'
# rpca_space = Embeddings(merged, pre_embedding)
# permuted_rpca_space = rpca_space
# ## We will now permute the integrated RPCA space.
# ##This is the same as the Permuted function, but we have swapped rows/columns because embeddings are cell x features while the scale.data slot is features x cells
# set.seed(1000)
# for (i in 1:dim(permuted_rpca_space)[2])
# {
#   row = pracma::randperm(dim(permuted_rpca_space)[1])
#   permuted_rpca_space[,i]=rpca_space[row,i]

# }
# ##Now, we will add our permuted RPCA space as a reduced dimension object. 
# data.permuted = merged
# data.permuted[[pre_embedding]] <- CreateDimReducObject(embeddings =permuted_rpca_space , key = "integratedrpca_", assay = DefaultAssay(data.permuted))

# ##Now we can call scDEED and 1) provide our own permuted object 2) specify that the pre_embedding space is integrated.cca
# result = scDEED(merged, K = K, pre_embedding = pre_embedding, permuted = data.permuted, reduction.method = 'umap')
# result$num_dubious

# # Minimum number of dubious cells at n_neighbors, min_dist
# m <- 0.1
# n <- 5

# merged <- RunUMAP(merged, dims = 1:5, min.dist = m, n.neighbors = n, seed.use = 100)

# ### Find neighbours and clusters
# merged <- FindNeighbors(merged, reduction = "pca", dims = 1:K)
# merged <- FindClusters(merged, resolution = 0.1, cluster.name = "seurat_clusters")


Idents(merged) <- 'seurat_clusters'
pdf('adult_aged_umap.pdf')
DimPlot(merged, reduction = "umap", label = TRUE)
dev.off()

saveRDS(merged, file = 'adult_aged_merged.RDS')

# Load in marker genes
marker_genes <- read.table('marker.genes.txt', header = FALSE)

# Plot marker genes
pdf('adult_aged_marker_genes.pdf', width = 10, height = 8)
DotPlot(merged, features = marker_genes$V1) + RotatedAxis() +
theme(axis.text.x = element_text(angle = 90, vjust=1, size=8))
dev.off()

pdf('adult_aged_marker_genes_heatmap.pdf', width = 10, height = 8)
DoHeatmap(merged, features = marker_genes$V1, assay = "SCT") + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

pdf('adult_aged_ventricular_CM_markers.pdf')
VlnPlot(merged, features = c("Nppa","Nppb","Pln","Tnnt2","Myh6"), group.by = "seurat_clusters")
dev.off()
table(Idents(merged), merged$sample)

merged$celltype <- plyr::mapvalues(Idents(merged),
from = 0:11,
to = c("FB","ventCM","ventCM","EC","EC","SMC",
"Macro","atrialCM","LEC","Epicard","Adipo","Glial"))

Idents(merged) <- 'celltype'

## Set cluster colours
cluster_colors <- list(Glial="#EC68A3", FB="#459AD5", LEC="#11B6EB", Adipo="#9188C0", 
                      SMC="#2EB7BE", Epicard="#AD7AB3", Macro="#36B28F", atrialCM="#59B031", 
                      EC="#9AA921", ventCM="#EE766F")


pdf('adult_aged_umap_celltype.pdf')
DimPlot(merged, reduction = "umap", label = FALSE, cols = cluster_colors)
dev.off()

saveRDS(merged, file = 'adult_aged_merged.RDS')
merged <- readRDS('adult_aged_merged.RDS')

pdf('batch_umap.pdf')
DimPlot(merged, reduction = "umap", group.by = "age")
dev.off()

# Export cell IDs for sinto - Allelome.PRO2
cell_ids <- data.frame(cell_id = colnames(merged), cell_type = merged$celltype, age = merged$age)
cell_ids$cell_type <- factor(cell_ids$cell_type, levels = names(sort(table(cell_ids$cell_type))))
cell_ids$age <- factor(cell_ids$age, levels = c("adult", "aged"))

cell_ids$cell_id <- gsub('adult_|aged_', '', cell_ids$cell_id)
cell_ids$age <- factor(cell_ids$age, levels = c('aged', 'adult'))

adult <- cell_ids[cell_ids$age == 'adult', 1:2]
aged <- cell_ids[cell_ids$age == 'aged', 1:2]

write.table(adult, file = '9w/cell_index.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(aged, file = '78w/cell_index.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

# Clustered bar chart of cell type percentage split by age
pdf('celltype_percentage_by_age.pdf')
ggplot(cell_ids, aes(x = age, fill = cell_type)) +
  geom_bar(position = "fill") +
  labs(x = "", y = "Cell composition", title = "") +
  theme_minimal() +
  scale_fill_manual(values = cluster_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, by = 0.25)) +
  coord_flip() +
  theme(legend.position = "bottom", legend.title = element_blank())
dev.off()

cell_ids <- data.frame(cell_id = colnames(merged), age = merged$age)
adult <- cell_ids[cell_ids$age == 'adult', , drop = FALSE]
aged <- cell_ids[cell_ids$age == 'aged', , drop = FALSE]
write.table(adult, file = 'pseudobulk/adult_cell_index.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(aged, file = 'pseudobulk/aged_cell_index.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

### Compare gene expression between adult and aged cells
DefaultAssay(merged) <- 'RNA'
merged <- JoinLayers(merged)

pseudobulk_fc <- function(
  obj,
  cell_ident,                 
  group_by,                   
  assay = "RNA",
  layer = "counts",
  contrast = NULL,            
  min_total_count = 10,       
  prior.count = 1             
) {
  stopifnot(group_by %in% colnames(obj[[]]))
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Install edgeR.")

  # 1) Subset to the cell type (by ident label or by explicit cells)
  if (is.character(cell_ident) && all(cell_ident %in% Idents(obj))) {
    x <- subset(obj, idents = cell_ident)
  } else {
    x <- subset(obj, cells = cell_ident)
  }

  # 2) Get raw counts (not scale.data)
  m <- LayerData(x, assay = "RNA", layer = "counts")
  if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")

  # 3) Build pseudobulk (sum counts per group)
  grp <- x[[group_by]][, 1, drop = TRUE]
  groups <- split(colnames(x), grp)
  pb_counts <- do.call(cbind, lapply(groups, function(cells)
    Matrix::rowSums(m[, cells, drop = FALSE])
  ))
  colnames(pb_counts) <- names(groups)

  # 4) Filter low-count genes (optional but sensible)
  keep <- Matrix::rowSums(pb_counts) >= min_total_count
  pb_counts <- pb_counts[keep, , drop = FALSE]

  # 5) Normalise to logCPM
  logCPM <- edgeR::cpm(pb_counts, prior.count = prior.count, log = TRUE)

  # 6) Work out contrast
  grps <- colnames(pb_counts)
  if (is.null(contrast)) {
    if (length(grps) != 2) {
      stop("Provide 'contrast = c(target, reference)' when there are not exactly 2 groups.")
    }
    contrast <- grps
  } else {
    stopifnot(length(contrast) == 2, all(contrast %in% grps))
  }
  fc <- logCPM[, contrast[1]] - logCPM[, contrast[2]]

  # 7) Tidy result
  res <- data.frame(
    gene = rownames(pb_counts),
    count_ref = pb_counts[, contrast[2]],
    count_tar = pb_counts[, contrast[1]],
    logCPM_ref = logCPM[, contrast[2]],
    logCPM_tar = logCPM[, contrast[1]],
    log2FC = as.numeric(fc),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  res <- res[order(-abs(res$log2FC)), ]

  list(
    counts = pb_counts,
    logCPM = logCPM,
    contrast = paste0(contrast[1], " vs ", contrast[2]),
    results = res
  )
}

result_list <- list()
for(celltype in levels(merged)){
  result_list[[celltype]] <- pseudobulk_fc(merged, cell_ident = celltype, group_by = "age", contrast = c("aged", "adult"))$res
}

saveRDS(result_list, file = 'adult_aged_log2FC_results.RDS')


library(ComplexHeatmap)
library(circlize)

common_genes <- Reduce(intersect, lapply(result_list, function(x) x$gene))
log2FC_matrix <- do.call(cbind, lapply(result_list, function(x) {
  x <- x[x$gene %in% common_genes, ]
  x <- x[match(common_genes, x$gene), ]
  x$log2FC
}))
rownames(log2FC_matrix) <- common_genes

pdf('adult_aged_log2FC_heatmap.pdf')
Heatmap(log2FC_matrix, name = "log2FC", col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        show_row_names = FALSE, show_column_names = TRUE,
        clustering_distance_columns = 'spearman', clustering_method_columns = 'ward.D2',
        cluster_rows = TRUE, cluster_columns = TRUE)
dev.off()


#### Read in Allelome.PRO2 output
load("adult.Allelome.PRO2.RData")
load("aged.Allelome.PRO2.RData")



# Create empty matrix with same dimensions as merged
allelic_matrix <- matrix(NA, nrow = nrow(merged), ncol = ncol(merged))
rownames(allelic_matrix) <- rownames(merged)
colnames(allelic_matrix) <- colnames(merged)

lapply(names(adult.Allelome), function(ct) {
  adult.Allelome_mean <- adult.Allelome[[ct]] %>%
    group_by(gene) %>%
    summarise(allelic_ratio = mean(allelic_ratio)) %>%
    ungroup() %>%
    data.frame()

  for(i in 1:nrow(adult.Allelome_mean)) {
    gene <- adult.Allelome_mean[i, "gene"]
    if(gene %in% rownames(allelic_matrix)) {
      allelic_matrix[gene, subset(cell_ids, cell_type == ct & age == "adult")$cell_id] <- adult.Allelome_mean[i, "allelic_ratio"]
    } else {
      next
    }
  }

  aged.Allelome_mean <- aged.Allelome[[ct]] %>%
    group_by(gene) %>%
    summarise(allelic_ratio = mean(allelic_ratio)) %>%
    ungroup() %>%
    data.frame()

  for(i in 1:nrow(aged.Allelome_mean)) {
    gene <- aged.Allelome_mean[i, "gene"]
    if(gene %in% rownames(allelic_matrix)) {
      allelic_matrix[gene, subset(cell_ids, cell_type == ct & age == "aged")$cell_id] <- aged.Allelome_mean[i, "allelic_ratio"]
    } else {
      next
    }
  }
})


library(data.table)
library(Matrix)

# 0) Prep labels
cell_ids$cell_id   <- as.character(cell_ids$cell_id)
cell_ids$cell_type <- as.character(cell_ids$cell_type)
cell_ids$age       <- as.character(cell_ids$age)
stopifnot(identical(colnames(merged), cell_ids$cell_id))

# 1) Make a group per cell: (cell_type × age)
cell_ids <- as.data.table(cell_ids)
cell_ids[, grp := paste(cell_type, age, sep = "||")]
grp_levels <- unique(cell_ids$grp)

# 2) Bind adult/aged Allelome into one long DT with a 'grp' column
#    (assumes each list element [[ct]] has columns: gene, allelic_ratio)
stack_one <- function(L, age_label) {
  if (length(L) == 0) return(data.table())
  rbindlist(lapply(names(L), function(ct) {
    if (is.null(L[[ct]]) || !nrow(L[[ct]])) return(NULL)
    data.table(
      gene = as.character(L[[ct]]$gene),
      allelic_ratio = as.numeric(L[[ct]]$allelic_ratio),
      grp  = paste(ct, age_label, sep = "||")
    )
  }), use.names = TRUE, fill = TRUE)
}

DT <- rbind(
  stack_one(adult.Allelome, "adult"),
  stack_one(aged.Allelome,  "aged"),
  use.names = TRUE, fill = TRUE
)

# 3) Per-gene mean within each group
DT <- DT[!is.na(gene) & !is.na(allelic_ratio)]
DTm <- DT[, .(allelic_ratio = mean(allelic_ratio, na.rm = TRUE)), by = .(gene, grp)]

# 4) Build the **genes × groups** matrix A, aligned to rownames(merged) and grp_levels
#    (fast dcast + reindex)
A <- dcast(DTm, gene ~ grp, value.var = "allelic_ratio", fun.aggregate = mean)
gene_order <- match(rownames(merged), A$gene)
# Create a full matrix with NA, then fill matching rows
A_full <- matrix(NA_real_, nrow = nrow(merged), ncol = length(grp_levels),
                 dimnames = list(rownames(merged), grp_levels))
ok <- !is.na(gene_order)
if (any(ok)) {
  # Reorder A's columns to grp_levels, dropping the 'gene' column
  A_mat <- as.matrix(A[, grp_levels, with = FALSE])  # columns may be missing -> handled below
  # If some groups absent in A, add them as NA columns
  missing_cols <- setdiff(grp_levels, colnames(A_mat))
  if (length(missing_cols)) {
    A_mat <- cbind(A_mat, matrix(NA_real_, nrow = nrow(A_mat), ncol = length(missing_cols),
                                 dimnames = list(NULL, missing_cols)))
    A_mat <- A_mat[, grp_levels, drop = FALSE]
  }
  A_full[ok, ] <- A_mat[gene_order[ok], , drop = FALSE]
}

# 5) Build **groups × cells** one-hot sparse matrix P
g_idx <- match(cell_ids$grp, grp_levels)
P <- sparseMatrix(i = g_idx,
                  j = seq_len(nrow(cell_ids)),
                  x = 1, dims = c(length(grp_levels), nrow(cell_ids)),
                  dimnames = list(grp_levels, cell_ids$cell_id))

# 6) Broadcast in one go: (genes × groups) %*% (groups × cells) = (genes × cells)
allelic_matrix <- A_full %*% P
# Convert to base matrix (optional)
allelic_matrix <- as.matrix(allelic_matrix)

# Add allelic_matrix to merged assay
merged <- CreateAssayObject(counts = allelic_matrix)

# Plot UMAP for Med14 and colour by allelic ratio
# From 0.5-0.9 shades of green
# 0.9 = Orange
# 1 = Red 
library(scales)

# Define color breaks
breaks <- seq(0.5, 1, by = 0.1)
colors <- c("green", "yellow", "orange", "red")

# Create a color palette
col_fun <- colorRamp2(breaks, colors)

pdf('allelic_ratio_Med14.pdf')
FeaturePlot(merged, features = "Med14", reduction = "umap", cols = col_fun(breaks), min.cutoff = 0.5, max.cutoff = 1) +
  theme_minimal() +
  labs(title = "Allelic Ratio of Med14") +
  theme(plot.title = element_text(hjust = 0.5))


############################################
# Visualising per cell expression of Hmgb1 #
############################################

# Feature Plot
pdf('Hmgb1_expression.pdf')
FeaturePlot(merged, features = "Hmgb1", reduction = "umap", cols = c("lightblue", "darkblue"), min.cutoff = "q10", max.cutoff = "q90")
dev.off()

# Violin Plot split by age and cell type
pdf('Hmgb1_violin.pdf', width = 12, height = 6)
VlnPlot(merged, features = "Hmgb1", group.by = "age", split.by = "celltype", cols = cluster_colors[levels(merged)])
dev.off()

# Percent of cells expressing Hmgb1 in each cell type and age group
hmgb1_expr <- merged@assays$SCT$counts["Hmgb1", ]
cell_ids$Hmgb1_expr <- hmgb1_expr[match(rownames(cell_ids), colnames(merged))]
hmgb1_summary <- cell_ids %>%
  group_by(cell_type, age) %>%
  summarise(percent_expressing = mean(Hmgb1_expr > 0) * 100) %>%
  ungroup()
hmgb1_summary$age <- factor(hmgb1_summary$age, levels = c("adult", "aged"))
pdf('Hmgb1_percent_expressing.pdf', width = 10, height = 6)
ggplot(hmgb1_summary, aes(x = cell_type, y = percent_expressing, fill = age)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Cell Type", y = "Percent Expressing Hmgb1", title = "Percentage of Cells Expressing Hmgb1 by Cell Type and Age") +
  theme_minimal() +
  scale_fill_manual(values = c("adult" = "lightblue", "aged" = "darkblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())
dev.off() 

####################################################
# Visualising per cell expression of 5930430L01Rik #
####################################################
DefaultAssay(merged) <- "SCT"
pdf('5930430L01Rik_expression.pdf')
FeaturePlot(merged, features = "5930430L01Rik", reduction = "umap", cols = c("lightblue", "darkblue"), min.cutoff = "q10", max.cutoff = "q90")
dev.off()

pdf('5930430L01Rik_violin.pdf', width = 12, height = 6)
VlnPlot(merged, features = "5930430L01Rik", group.by = "age", split.by = "celltype", cols = cluster_colors[levels(merged)])
dev.off()
  
Rik_expr <- merged@assays$SCT$counts["5930430L01Rik", ]
cell_ids$Rik_expr <- Rik_expr[match(rownames(cell_ids), colnames(merged))]
rik_summary <- cell_ids %>%
  group_by(cell_type, age) %>%
  summarise(percent_expressing = mean(Rik_expr > 0) * 100) %>%
  ungroup() 
rik_summary$age <- factor(rik_summary$age, levels = c("adult", "aged"))
pdf('5930430L01Rik_percent_expressing.pdf', width = 10, height = 6)
ggplot(rik_summary, aes(x = cell_type, y = percent_expressing, fill = age)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Cell Type", y = "Percent Expressing 5930430L01Rik", title = "Percentage of Cells Expressing 5930430L01Rik by Cell Type and Age") +
  theme_minimal() +
  scale_fill_manual(values = c("adult" = "lightblue", "aged" = "darkblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())
dev.off()


##############################################
# Visualising per cell expression of Gm15408 #
##############################################
pdf('Gm15408_expression.pdf')
FeaturePlot(merged, features = "Gm15408", reduction = "umap", cols = c("lightblue", "darkblue"), min.cutoff = "q10", max.cutoff = "q90")
dev.off()
pdf('Gm15408_violin.pdf', width = 12, height = 6)
VlnPlot(merged, features = "Gm15408", group.by = "age", split.by = "celltype", cols = cluster_colors[levels(merged)])
dev.off() 

Gm15408_expr <- merged@assays$SCT$counts["Gm15408", ]
cell_ids$Gm15408_expr <- Gm15408_expr[match(rownames(cell_ids), colnames(merged))]
Gm15408_summary <- cell_ids %>%
  group_by(cell_type, age) %>%
  summarise(percent_expressing = mean(Gm15408_expr > 0) * 100) %>%
  ungroup()
Gm15408_summary$age <- factor(Gm15408_summary$age, levels = c("adult", "aged"))
pdf('Gm15408_percent_expressing.pdf', width = 10, height = 6)
ggplot(Gm15408_summary, aes(x = cell_type, y = percent_expressing, fill = age)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Cell Type", y = "Percent Expressing Gm15408", title = "Percentage of Cells Expressing Gm15408 by Cell Type and Age") +
  theme_minimal() +
  scale_fill_manual(values = c("adult" = "lightblue", "aged" = "darkblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())
dev.off()