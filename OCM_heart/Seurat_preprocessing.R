library(Seurat)
library(ddqcR)
library(glmGamPoi)
options(future.globals.maxSize = 8 * 1024^3)

# ! important ! run this command below to modify input file for Seurat
# ptrepack --complevel 5 tiny_output_filtered.h5:/matrix tiny_output_filtered_seurat.h5:/matrix

# Local override of ddqcR::initialQC with corrected feature matching
initialQC <- function(
	data,
	basic.n.genes = 100,
	basic.percent.mt = 80,
	mt.prefix = "MT-",
	rb.prefix = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"
) {
	mt.features <- grep(mt.prefix, rownames(data), ignore.case = TRUE, value = TRUE)
	rb.features <- grep(rb.prefix, rownames(data), ignore.case = TRUE, value = TRUE)

	data[["percent.mt"]] <- if (length(mt.features) > 0) {
		PercentageFeatureSet(data, features = mt.features)
	} else {
		rep(0, ncol(data))
	}

	data[["percent.rb"]] <- if (length(rb.features) > 0) {
		PercentageFeatureSet(data, features = rb.features)
	} else {
		rep(0, ncol(data))
	}

	data <- subset(data, subset = nFeature_RNA >= basic.n.genes & percent.mt <= basic.percent.mt)
	data
}
mt.prefix <- "mt-"
rb.prefix <- "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa"


# 9w mouse heart
heart_9w <- Read10X_h5('Hearts_OCM/outs/per_sample_outs/He_9w_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered_seurat.h5')
heart_9w <- CreateSeuratObject(heart_9w)
heart_9w
# Filter with ddqcR
heart_9w <- initialQC(heart_9w, mt.prefix = mt.prefix, rb.prefix = rb.prefix)
summary(heart_9w$percent.mt)
summary(heart_9w$percent.rb)
pdf("ddqcR_QC_metrics_heart_9w.pdf")
df.qc <- ddqc.metrics(heart_9w)
dev.off()
heart_9w <- filterData(heart_9w, df.qc)
# Add sample to metadata
heart_9w$sample <- "9w"
# SCTransform normalisation
heart_9w <- SCTransform(heart_9w, vars.to.regress = c("percent.mt"), verbose = FALSE)

# 78w mouse heart
heart_78w <- Read10X_h5('Hearts_OCM/outs/per_sample_outs/He_78w_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered_seurat.h5')
heart_78w <- CreateSeuratObject(heart_78w)
heart_78w
# Filter with ddqcR
heart_78w <- initialQC(heart_78w, mt.prefix = mt.prefix, rb.prefix = rb.prefix)
summary(heart_78w$percent.mt)
summary(heart_78w$percent.rb)
pdf("ddqcR_QC_metrics_heart_78w.pdf")
df.qc <- ddqc.metrics(heart_78w)
dev.off()
heart_78w <- filterData(heart_78w, df.qc)
heart_78w$sample <- "78w"
# SCTransform normalisation
heart_78w <- SCTransform(heart_78w, vars.to.regress = c("percent.mt"), verbose = FALSE)

# TAC mouse heart
heart_TAC <- Read10X_h5('Hearts_OCM/outs/per_sample_outs/He_TAC_28d_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered_seurat.h5')
heart_TAC <- CreateSeuratObject(heart_TAC)
heart_TAC
# Filter with ddqcR
heart_TAC <- initialQC(heart_TAC, mt.prefix = mt.prefix, rb.prefix = rb.prefix)
summary(heart_TAC$percent.mt)
summary(heart_TAC$percent.rb)
pdf("ddqcR_QC_metrics_heart_TAC.pdf")
df.qc <- ddqc.metrics(heart_TAC)
dev.off()
heart_TAC <- filterData(heart_TAC, df.qc)
heart_TAC$sample <- "TAC"
# SCTransform normalisation
heart_TAC <- SCTransform(heart_TAC, vars.to.regress = c("percent.mt"), verbose = FALSE)

# Sham mouse heart
heart_Sham <- Read10X_h5('Hearts_OCM/outs/per_sample_outs/He_Sham_28d_snRNA_XBxC_het_XX/count/sample_raw_feature_bc_matrix.cellbender_filtered_seurat.h5')
heart_Sham <- CreateSeuratObject(heart_Sham)
heart_Sham
# Filter with ddqcR
heart_Sham <- initialQC(heart_Sham, mt.prefix = mt.prefix, rb.prefix = rb.prefix)
summary(heart_Sham$percent.mt)
summary(heart_Sham$percent.rb)
pdf("ddqcR_QC_metrics_heart_Sham.pdf")
df.qc <- ddqc.metrics(heart_Sham)
dev.off()
heart_Sham <- filterData(heart_Sham, df.qc)
heart_Sham$sample <- "Sham"
# SCTransform normalisation
heart_Sham <- SCTransform(heart_Sham, vars.to.regress = c("percent.mt"), verbose = FALSE)


# Merge all samples into one Seurat object
heart <- merge(heart_9w, y = c(heart_78w, heart_TAC, heart_Sham), add.cell.ids = c("9w", "78w", "TAC", "Sham"))
heart@meta.data$sample <- factor(heart@meta.data$sample, levels = c("9w", "78w", "Sham", "TAC"))

# Violin plots of key metrics per samples
pdf("QC_violin_plots_heart.pdf")
VlnPlot(heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), group.by = "sample", ncol = 2, pt.size = 0)
dev.off()

# SCTransform normalisation on merged object (raw counts) and regressing out percent.mt
DefaultAssay(heart) <- "RNA"
heart <- SCTransform(heart, vars.to.regress = "percent.mt", verbose = FALSE)

# PCA and UMAP
heart <- RunPCA(heart, features = VariableFeatures(heart), npcs = 50)

pdf("Heart_ElbowPlot.pdf")
ElbowPlot(heart, ndims = 50)
dev.off()

# Construct UMAP embedding
heart <- RunUMAP(heart, dims = 1:20)

heart <- readRDS('heart_seurat_object_SCT.rds')


# Clustering
heart <- FindNeighbors(heart, dims = 1:20)
heart <- FindClusters(heart, resolution = 0.1)

pdf("Heart_UMAP_samples.pdf")
DimPlot(heart, reduction = "umap", group.by = "sample", raster=FALSE)
dev.off()

pdf("Heart_UMAP_splitsample.pdf")
DimPlot(heart, reduction = "umap", group.by = "sample", split.by = "sample", raster=FALSE)
dev.off()

pdf("Heart_UMAP_clusters_0.1.pdf")
DimPlot(heart, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

pdf('Heart_UMAP_markers.pdf')
FeaturePlot(heart, features = c("Tnnt2", "Myh6", "Actc1", "Col1a1", "Pecam1", "Ptprc"))
dev.off()

saveRDS(heart, file = "heart_seurat_object_SCT.rds")
heart <- readRDS("heart_seurat_object_SCT.rds")

# Find markers for each cluster
heart <- PrepSCTFindMarkers(heart)
heart.markers <- FindAllMarkers(heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

# Save markers to file
write.table(heart.markers, file = "heart_cluster_markers_logfc_0.5_resolution_0.1.RData")

lapply(split(heart.markers, heart.markers$cluster), function(x){
    x[order(x$avg_log2FC, decreasing = TRUE),][1:5,'gene']
})

pdf('fibroblast_markers.pdf')
FeaturePlot(heart, features = c("Col1a1", "Col1a2", "Dcn", "Postn", "Pdgfra"))
dev.off()

# Using Sarah's marker genes
marker_genes <- read.table('../adult_aged_heart_snRNAseq/marker.genes.txt', header = FALSE)

# Plot marker genes
pdf('marker_genes_dotplot.pdf', width = 10, height = 8)
DotPlot(heart, features = marker_genes$V1) + RotatedAxis() +
theme(axis.text.x = element_text(angle = 90, vjust=1, size=8))
dev.off()

# Downsample to equal cells per cluster so each cluster occupies the same width
cells_equal <- unlist(lapply(levels(Idents(heart)), function(cl) {
  cells <- WhichCells(heart, idents = cl)
  if (length(cells) > 1000) sample(cells, 1000) else cells
}))

marker_genes <- marker_genes[marker_genes$V1 %in% rownames(heart), ]

# 1. Switch the default assay back to RNA
DefaultAssay(heart) <- "RNA"

# 2. Normalize the RNA data (if you haven't already done so previously)
heart <- NormalizeData(heart)

# 3. Scale only the genes you care about for the heatmap
heart <- ScaleData(heart)
pdf('marker_genes_heatmap_0.1.pdf', width = 10, height = 8)
DoHeatmap(
  heart,
  features = marker_genes$V1,
  assay = "RNA",
  cells = cells_equal
) + scale_fill_gradient2(low = "blue", mid = "white", high = "red")
dev.off()


save(heart, file = "heart_seurat_object_SCT.RData")
