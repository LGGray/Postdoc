library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

localdir <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial/9w/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))

DefaultAssay(object) <- "Spatial.008um"

vln.plot <- VlnPlot(object, features = 'nCount_Spatial.008um', pt.size = 0) + 
              theme(axis.text = element_text(size = 4)) + NoLegend()

count.plot <- SpatialFeaturePlot(object, features = 'nCount_Spatial.008um') + 
                theme(legend.position = "right")

# Note that many spots have very few counts, in-part due to low cellular density in certain tissue regions
pdf('test.figures/spatial_008um_count_distribution.pdf')
vln.plot | count.plot
dev.off()

# Normalize both 8um and 16um bins
DefaultAssay(object) <- "Spatial.008um"
object <- NormalizeData(object)

DefaultAssay(object) <- "Spatial.016um"
object <- NormalizeData(object)

# A few circular high-count hotspots (bubbles/debris) sit far above the tissue
# median and saturate the colour scale, so clip every feature plot at q02/q98.
lo <- "q02"; hi <- "q98"

# Switch spatial resolution to 16um from 8um
# Ttn: pan-cardiomyocyte sarcomere gene - confirms myocardium, near-uniform
DefaultAssay(object) <- "Spatial.016um"
p1 <- SpatialFeaturePlot(object, features = "Ttn",
                         min.cutoff = lo, max.cutoff = hi) +
        ggtitle("Ttn expression (16um)")

# Switch back to 8um
# Myl7: atrial cardiomyocytes - marks the atrial region against ventricle
DefaultAssay(object) <- "Spatial.008um"
p2 <- SpatialFeaturePlot(object, features = "Myl7",
                         min.cutoff = lo, max.cutoff = hi) +
        ggtitle("Myl7 expression (8um)")

pdf('test.figures/spatial_008um_16um_expression.pdf')
p1 | p2
dev.off()

# Cardiomyocyte markers are abundant enough to resolve at 8um.
# Nppa  - atrial + stressed myocardium, strongest chamber-level contrast
# Myl7  - atrial cardiomyocytes
# Myl2  - ventricular (compact myocardium) cardiomyocytes
# Myh7  - beta-MHC, subendocardial / stress gradient
# Ttn   - pan-cardiomyocyte, sarcomere (tissue outline / CM density)
cm.markers <- c("Nppa", "Myl7", "Myl2", "Myh7", "Ttn")

# Stromal/vascular markers are far sparser per bin and read as near-empty at
# 8um, where most bins carry 0 counts for them. Bin at 16um to recover signal.
# Postn - activated fibroblasts, fibrotic / perivascular regions
# Col1a1- interstitial and perivascular collagen
# Pecam1- endothelium, vasculature
stromal.markers <- c("Postn", "Col1a1", "Pecam1", "Vwf", "Dcn")

DefaultAssay(object) <- "Spatial.008um"
cm.markers <- cm.markers[cm.markers %in% rownames(object)]
pdf('test.figures/spatial_008um_cardiac_markers.pdf', width = 12, height = 8)
SpatialFeaturePlot(object, features = cm.markers, ncol = 3,
                   min.cutoff = lo, max.cutoff = hi)
dev.off()

DefaultAssay(object) <- "Spatial.016um"
stromal.markers <- stromal.markers[stromal.markers %in% rownames(object)]
pdf('test.figures/spatial_016um_stromal_markers.pdf', width = 12, height = 8)
SpatialFeaturePlot(object, features = stromal.markers, ncol = 3,
                   min.cutoff = lo, max.cutoff = hi)
dev.off()

# Note that data is already normalized
DefaultAssay(object) <- "Spatial.008um"
object <- FindVariableFeatures(object)
object <- ScaleData(object)

# We select 50,0000 cells and create a new 'sketch' assay
object <- SketchData(object = object,
                     ncells = 50000,
                     method = "LeverageScore",
                     sketched.assay = "sketch")
# Switch analysis to sketched cells
DefaultAssay(object) <- "sketch"

# Perform the clustering workflow
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")

# ElbowPlot to determine the number of PCs to use for clustering
pdf('test.figures/sketch_elbow_plot.pdf')
ElbowPlot(object, ndims = 50, reduction = "pca.sketch")
dev.off()

object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:10)

# resolution = 1 gave 67 clusters, far more than heart supports: the elbow puts
# real signal in ~6-8 PCs, and at 8um a bin is a fragment of a single (large)
# cardiomyocyte, so much of that structure is sub-cellular/depth, not identity.
# Sweep instead and pick the resolution where cluster count plateaus.
# The first sweep (0.1-1.0) never plateaued: 16, 25, 40, 45, 48, 61, 68. A count
# that climbs monotonically means the graph is fragmenting a continuum rather
# than resolving discrete populations. Shattering starts at 0.2->0.3 (+15), so
# re-sweep below 0.2 to test whether ~16 clusters is genuine structure.
res.sweep <- c(0.02, 0.05, 0.08, 0.1, 0.15, 0.2)
object <- FindClusters(object, resolution = res.sweep, verbose = FALSE)

n.clust <- sapply(res.sweep, function(r)
  length(unique(object[[paste0("sketch_snn_res.", r)]][, 1])))
print(data.frame(resolution = res.sweep, n_clusters = n.clust))

# Working resolution - revisit once the sweep above is inspected
object$seurat_cluster.sketched <- object[["sketch_snn_res.0.1"]][, 1]
Idents(object) <- "seurat_cluster.sketched"

object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:10)

object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:10,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label=F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "none")

# switch to full dataset
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label=F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "none")

pdf('test.figures/sketch_clustering.pdf', width = 12, height = 6)
p1 | p2
dev.off()

pdf('test.figures/sketch_spatial.pdf', width = 12, height = 6)
SpatialDimPlot(object, label = T, repel = T, label.size = 4)
dev.off()

# ---------------------------------------------------------------------------
# Clustering at 16um
#
# The 8um sweep never plateaued at any resolution (0.02-1.0 gave 4, 10, 14, 16,
# 21, 25, 40, 45, 48, 61, 68 - essentially linear). Cluster count there is a
# free parameter, not a property of the data: an 8um bin is a fragment of a
# single large cardiomyocyte, so bin transcriptomes form a continuum. 16um bins
# are closer to whole cells. 4x fewer bins than 8um, so cluster directly rather
# than sketching.
# ---------------------------------------------------------------------------
DefaultAssay(object) <- "Spatial.016um"
cat("16um bins:", ncol(object[["Spatial.016um"]]), "\n")

object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "Spatial.016um", reduction.name = "pca.016um")

pdf('test.figures/elbow_016um.pdf')
ElbowPlot(object, ndims = 50, reduction = "pca.016um")
dev.off()

object <- FindNeighbors(object, assay = "Spatial.016um",
                        reduction = "pca.016um", dims = 1:15)

# Sweep gave 7, 10, 14, 18, 20, 26, 31 over 0.1-1.5. No true plateau: the flat
# stretch at 0.5-0.8 (~7 clusters per unit of resolution) re-accelerates to ~30
# by 0.8-1.0. But the whole curve is much shallower than 8um (10-40 per unit
# vs 40-150), and the elbow moves from ~7 to ~15 real PCs, so 16um bins are
# genuinely closer to whole cells - a gentler continuum, not a discrete one.
# 0.5 (18 clusters) is a working choice justified by marker annotation below,
# NOT by a plateau. Report it that way.
res.sweep.16 <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5)
# graph.name set explicitly: FindClusters otherwise infers it from DefaultAssay,
# which silently picks the wrong graph if the assay drifted in an interactive run
object <- FindClusters(object, graph.name = "Spatial.016um_snn",
                       resolution = res.sweep.16, verbose = FALSE)

n.clust.16 <- sapply(res.sweep.16, function(r)
  length(unique(object[[paste0("Spatial.016um_snn_res.", r)]][, 1])))
print(data.frame(resolution = res.sweep.16, n_clusters = n.clust.16))

object$cluster.016um <- object[["Spatial.016um_snn_res.0.5"]][, 1]
Idents(object) <- "cluster.016um"

object <- RunUMAP(object, reduction = "pca.016um",
                  reduction.name = "umap.016um", dims = 1:15)

pdf('test.figures/clustering_016um.pdf', width = 14, height = 6)
DimPlot(object, reduction = "umap.016um", label = TRUE) |
  SpatialDimPlot(object, label = TRUE, repel = TRUE, label.size = 3)
dev.off()

# Acceptance test: does the clustering separate atrium from ventricle?
# Nppa/Myl7 mark atrium and Myl2 marks ventricle in the marker panel, so any
# usable resolution must split those into different clusters.
pdf('test.figures/chamber_check_016um.pdf', width = 10, height = 6)
VlnPlot(object, features = c("Nppa", "Myl7", "Myl2"),
        stack = TRUE, flip = TRUE) + NoLegend()
dev.off()

# ---------------------------------------------------------------------------
# Which clusters are real?
#
# At res 0.5 the chamber check passes (Nppa+ in 12/13 only, spatially confined
# to the atrial patch), but clusters 0-10 and 15 are salt-and-pepper across the
# ventricle. That is expected for genuinely dispersed cell types (capillary
# endothelium, fibroblasts, macrophages) and NOT expected for depth/technical
# splits, so markers are what separate the two cases. A cluster with no
# meaningful positive markers is a slice of the cardiomyocyte continuum and
# should be merged.
# ---------------------------------------------------------------------------
markers.016um <- FindAllMarkers(object, only.pos = TRUE,
                                min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers.016um, 'test.figures/markers_016um_res0.5.csv', row.names = FALSE)

top.markers <- markers.016um %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)
print(as.data.frame(top.markers), max = 1000)

# Clusters returning few or no markers are continuum slices, not populations
n.markers <- table(markers.016um$cluster)
print(n.markers)

# Is a cluster driven by sequencing depth rather than identity?
pdf('test.figures/cluster_qc_016um.pdf', width = 12, height = 5)
VlnPlot(object, features = c("nCount_Spatial.016um", "nFeature_Spatial.016um"),
        pt.size = 0, ncol = 2) + NoLegend()
dev.off()

# ---------------------------------------------------------------------------
# Annotation at res 0.5
#
# NB judge clusters by effect size, not marker count: 9 (Acta2/Tagln log2FC ~7),
# 3 (Rgs5), 5 (Lyz2), 6/15 (Hba/Hbb) each return 1-5 markers but are
# unambiguous cell types. Filtering on marker count would discard the entire
# vascular and immune compartment.
#
# 14 and 16 are near-empty bins (~50 and ~20 counts vs ~600 elsewhere), which
# is why their "markers" have pct.1 << pct.2 - they are depleted for
# everything. Low library size also inflates normalised values, so these bins
# produce spurious bright streaks in SpatialFeaturePlot. Drop them.
# ---------------------------------------------------------------------------
cluster.labels <- c(
  "0"  = "Ventricular CM",        # Gja1, Casq2
  "1"  = "Fibroblast",            # Gsn, Dcn
  "2"  = "Endothelium",           # Fabp4
  "3"  = "Pericyte",              # Rgs5
  "4"  = "Stressed CM",           # Xirp2, Ankrd1
  "5"  = "Macrophage",            # Lyz2
  "6"  = "Erythrocyte",           # Hbb-bs, Hba-a1
  "7"  = "Ventricular CM",        # Cox7a1 only - continuum slice
  "8"  = "Atrial CM",             # Myl4, Myl7
  "9"  = "Smooth muscle",         # Acta2, Tagln, Myl9
  "10" = "Epicardium/adipocyte",  # Mgst1, Fabp5
  "11" = "Fibroblast",            # Dcn, Gsn, Mgp, Dpt, Cfh
  "12" = "Stressed CM",           # Acta1, Ankrd1, Myh7, Nppa, Nppb
  "13" = "Atrial CM",             # Myl4, Myl7, Sln, Mybphl
  "14" = "Low-count artefact",    # ~50 counts
  "15" = "Erythrocyte",           # Hba/Hbb
  "16" = "Low-count artefact"     # ~20 counts
)
# unname() is required: the lookup carries cluster IDs as names, and Seurat
# matches names against cell barcodes rather than assigning positionally
celltype <- unname(cluster.labels[as.character(object$cluster.016um)])
names(celltype) <- colnames(object)
stopifnot(!anyNA(celltype))
object$celltype.016um <- celltype

# Drop the near-empty bins before any downstream analysis or figures
keep <- colnames(object)[object$celltype.016um != "Low-count artefact"]
object.clean <- subset(object, cells = keep)

Idents(object.clean) <- "celltype.016um"
pdf('test.figures/annotated_016um.pdf', width = 14, height = 6)
DimPlot(object.clean, reduction = "umap.016um", label = TRUE, repel = TRUE) |
  SpatialDimPlot(object.clean, label = FALSE)
dev.off()