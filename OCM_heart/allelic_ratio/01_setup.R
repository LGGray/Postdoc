# ---------------------------------------------------------------------------
# 01 - Setup: load the Seurat object, subcluster the ventricular
# cardiomyocytes, and produce the basic QC / Xist expression plots.
#
# Run from the OCM_heart/ directory:  Rscript allelic_ratio/01_setup.R
# Writes: heart_seurat_object_SCT.rds (adds the celltype_sub column)
#
# Run this FIRST -- 02-04 expect the saved object to carry celltype_sub.
# ---------------------------------------------------------------------------
source("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/OCM_heart/allelic_ratio/00_functions.R")

############################
# Read in single cell data #
############################
heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

# Subcluster the Ventricular Cardiomyocytes to check for the two suspected
# subpopulations. Kept as a separate metadata column (celltype_sub) rather
# than overwriting `celltype`, since downstream code filters on the exact
# string "Ventricular Cardiomyocytes" and shouldn't silently lose those cells.
vcm <- subset(heart, subset = celltype == "Ventricular Cardiomyocytes")

DefaultAssay(vcm) <- "RNA"
vcm <- SCTransform(vcm, verbose = FALSE)
vcm <- RunPCA(vcm, features = VariableFeatures(vcm), npcs = 30)

pdf("Allelic_ratio_results/VCM_subcluster_ElbowPlot.pdf")
ElbowPlot(vcm, ndims = 30)
dev.off()

# NB: check the elbow plot output and adjust `dims` below if needed before
# trusting the clustering
vcm <- RunUMAP(vcm, dims = 1:15)
vcm <- FindNeighbors(vcm, dims = 1:15)
vcm <- FindClusters(vcm, resolution = 0.05)

pdf("Allelic_ratio_results/VCM_subcluster_UMAP.pdf")
DimPlot(vcm, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

# Check whether the split tracks sample/condition rather than a real
# biological subpopulation (e.g. batch effect)
pdf("Allelic_ratio_results/VCM_subcluster_UMAP_by_sample.pdf")
DimPlot(vcm, reduction = "umap", group.by = "sample", raster = FALSE)
dev.off()

# Markers distinguishing the subclusters
vcm <- PrepSCTFindMarkers(vcm)
vcm.markers <- FindAllMarkers(vcm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.table(vcm.markers, file = "Allelic_ratio_results/VCM_subcluster_markers.txt")

pdf("Allelic_ratio_results/VCM_subcluster_markers_dotplot.pdf", width = 10, height = 8)
DotPlot(vcm, features = unique(unlist(lapply(split(vcm.markers, vcm.markers$cluster), function(x) {
  x[order(x$avg_log2FC, decreasing = TRUE), "gene"][1:10]
})))) + RotatedAxis() +
theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8))
dev.off()

# depth + complexity per subcluster — is this a technical split?
pdf("Allelic_ratio_results/VCM_subcluster_QC_violinplots.pdf", width = 10, height = 8)
VlnPlot(vcm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "seurat_clusters", pt.size = 0)
dev.off()
table(vcm$seurat_clusters, vcm$sample)



# Push subcluster labels back onto the full heart object as a new column
heart$celltype_sub <- as.character(heart$celltype)
match_idx <- match(Cells(vcm), colnames(heart))
heart$celltype_sub[match_idx] <- paste0("Ventricular Cardiomyocytes_", Idents(vcm))
heart$celltype_sub <- factor(heart$celltype_sub)

heart_no_cluster1 <- subset(heart, subset = celltype_sub != "Ventricular Cardiomyocytes_1")

# Plot UMAP with subcluster labels
pdf("Allelic_ratio_results/heart_UMAP_with_VCM_subclusters.pdf")
DimPlot(heart_no_cluster1, reduction = "umap", group.by = "celltype_sub", label = TRUE, raster = FALSE) +
  theme(legend.position = "none")
dev.off()

saveRDS(heart, file = "heart_seurat_object_SCT.rds")

heart_metadata <- heart@meta.data
# box plot of nCount_RNA per sample, split by celltype
pdf('Allelic_ratio_results/nCount_RNA_boxplot_by_celltype.pdf')
ggplot(heart_metadata, aes(x = celltype, y = nCount_RNA, fill = sample)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~sample) +
  labs(y = "nCount_RNA", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.off()

# Boxplot of Xist expression (approximate TPM) per sample, split by celltype
# Xist gene length from mm10 chrX:102,503,972-102,537,573
xist_length_kb <- (102537573 - 102503972 + 1) / 1000

xist_counts <- GetAssayData(heart, assay = "RNA", layer = "counts")["Xist", ]
total_counts <- heart$nCount_RNA

xist_df <- data.frame(
  sample = heart$sample,
  celltype = heart$celltype,
  Xist_TPM = (xist_counts / xist_length_kb) / (total_counts / 1e6)
)

pdf('Allelic_ratio_results/xist_expression_TPM_boxplot_by_celltype.pdf')
ggplot(xist_df, aes(x = sample, y = Xist_TPM, fill = sample)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~celltype) +
  labs(y = "Xist TPM (approx.)", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.off()

