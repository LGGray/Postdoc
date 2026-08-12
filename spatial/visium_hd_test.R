library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


localdir <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial/9w/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = 8)

DefaultAssay(object) <- "Spatial.008um"

vln.plot <- VlnPlot(object, features = 'nCount_Spatial.008um', pt.size = 0) + 
              theme(axis.text = element_text(size = 4)) + NoLegend()

count.plot <- SpatialFeaturePlot(object, features = 'nCount_Spatial.008um') + 
                theme(legend.position = "right")

# Note that many spots have very few counts, in-part due to low cellular density in certain tissue regions
pdf('test.figures/spatial_008um_count_distribution.pdf')
vln.plot | count.plot
dev.off()

# Normalize
object <- NormalizeData(object)

# A few circular high-count hotspots (bubbles/debris) sit far above the tissue
# median and saturate the colour scale, so clip every feature plot at q02/q98.
lo <- "q02"; hi <- "q98"

# Ttn: pan-cardiomyocyte sarcomere gene - confirms myocardium, near-uniform
DefaultAssay(object) <- "Spatial.008um"
p1 <- SpatialFeaturePlot(object, features = "Ttn",
                         min.cutoff = lo, max.cutoff = hi) +
        ggtitle("Ttn expression (8um)")

# Myl7: atrial cardiomyocytes - marks the atrial region against ventricle
p2 <- SpatialFeaturePlot(object, features = "Myl7",
                         min.cutoff = lo, max.cutoff = hi) +
        ggtitle("Myl7 expression (8um)")

pdf('test.figures/spatial_008um_Ttn_Myl7_expression.pdf')
p1 | p2
dev.off()

# Cardiomyocyte markers are abundant enough to resolve at 16um.
# Nppa  - atrial + stressed myocardium, strongest chamber-level contrast
# Myl7  - atrial cardiomyocytes
# Myl2  - ventricular (compact myocardium) cardiomyocytes
# Myh7  - beta-MHC, subendocardial / stress gradient
# Ttn   - pan-cardiomyocyte, sarcomere (tissue outline / CM density)
cm.markers <- c("Nppa", "Myl7", "Myl2", "Myh7", "Ttn")

# Stromal/vascular markers are far sparser per bin and read as near-empty at
# 16um, where most bins carry 0 counts for them. Bin at 16um to recover signal.
# Postn - activated fibroblasts, fibrotic / perivascular regions
# Col1a1- interstitial and perivascular collagen
# Pecam1- endothelium, vasculature
stromal.markers <- c("Postn", "Col1a1", "Pecam1", "Vwf", "Dcn")

# Visualise X chromosome markers
chrX.markers <- c("Xist", "Tsix", "Kdm5c", "Eif2s3x", "Ddx3x")

cm.markers <- cm.markers[cm.markers %in% rownames(object)]
pdf('test.figures/spatial_008um_cardiac_markers.pdf', width = 12, height = 8)
SpatialFeaturePlot(object, features = cm.markers, ncol = 3,
                   min.cutoff = lo, max.cutoff = hi)
dev.off()

stromal.markers <- stromal.markers[stromal.markers %in% rownames(object)]
pdf('test.figures/spatial_008um_stromal_markers.pdf', width = 12, height = 8)
SpatialFeaturePlot(object, features = stromal.markers, ncol = 3,
                   min.cutoff = lo, max.cutoff = hi)
dev.off()

chrX.markers <- chrX.markers[chrX.markers %in% rownames(object)]
pdf('test.figures/spatial_008um_chrX_markers.pdf', width = 12, height = 8)
SpatialFeaturePlot(object, features = chrX.markers, ncol = 3,
                   min.cutoff = lo, max.cutoff = hi)
dev.off()

pdf('test.figures/spatial_008um_DMD_plot.pdf')
SpatialFeaturePlot(object, features = "Dmd",
                   min.cutoff = lo, max.cutoff = hi)
dev.off()


# Note that data is already normalized
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
pdf('test.figures/sketch_elbow_plot_008um.pdf')
ElbowPlot(object, ndims = 50, reduction = "pca.sketch")
dev.off()

object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:15)


res.sweep <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5)
object <- FindClusters(object, resolution = res.sweep, verbose = FALSE)

n.clust <- sapply(res.sweep, function(r)
  length(unique(object[[paste0("sketch_snn_res.", r)]][, 1])))
print(data.frame(resolution = res.sweep, n_clusters = n.clust))

# Working resolution - revisit once the sweep above is inspected
object$seurat_cluster.sketched <- object[["sketch_snn_res.0.1"]][, 1]
Idents(object) <- "seurat_cluster.sketched"

object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:15)

object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:15,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label=F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "none")

# switch to full dataset
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label=F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "none")

pdf('test.figures/sketch_clustering_008um.pdf', width = 12, height = 6)
p1 | p2
dev.off()

pdf('test.figures/sketch_spatial_008um.pdf', width = 12, height = 6)
SpatialDimPlot(object, label = T, repel = T, label.size = 4)
dev.off()

Idents(object) <- "seurat_cluster.projected"
cells <- CellsByIdentities(object, idents=0:length(levels(object))-1)
p <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], 
                    cols.highlight = c("#FFFF00","grey50"), facet.highlight = T, combine=T) + NoLegend()

pdf('test.figures/sketch_spatial_highlight_008um.pdf', width = 12, height = 6)
print(p)
dev.off()

# Create downsampled object to make visualization easier
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
object_subset <- subset(object, cells = Cells(object[['Spatial.008um']]), downsample=1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = 'Spatial.008um', only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
pdf('test.figures/top_gene_heatmap.pdf', width = 12, height = 6)
print(p)
dev.off()

DefaultAssay(object) <- "Spatial.008um"
object <- FindVariableFeatures(object, nfeatures = 3000)
vf <- grep("^(mt-|Hb[ab]-|Rp[sl])", VariableFeatures(object),
           value = TRUE, invert = TRUE)
VariableFeatures(object) <- vf


saveRDS(object, file = "spatial_008um_seurat_object.RDS")

object <- readRDS("spatial_008um_seurat_object.RDS")

.libPaths(c("~/R/matrix-dev", .libPaths()))
library(SeuratWrappers)
library(Banksy)
source("../Postdoc/spatial/banksy_patch.R")


# k_geom : Local neighborhood size. Larger values will yield larger domains
# lambda : Influence of the neighborhood. Larger values yield more spatially coherent domains
object <- RunBanksy(object, lambda = 0.8, verbose = TRUE,
                    assay = 'Spatial.008um', slot = 'data',
                    features = 'variable',
                    k_geom = 200, use_agf = FALSE, chunk_size = 5000)

saveRDS(object, file = "spatial_008um_seurat_object_banksy.RDS")

DefaultAssay(object) <- "BANKSY"
object <- RunPCA(object, assay = 'BANKSY', reduction.name = "pca.banksy", features = rownames(object), npcs = 30)

pdf('test.figures/banksy_elbow_plot.pdf')
ElbowPlot(object, reduction = "pca.banksy", ndims = 30)
dev.off()

object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:15)

# Resolution sweep. res 0.1 gave 13 domains, but only one was a distinct tissue:
# the Myl7+/Myl2- atrial appendage at top-right (45.6% Myl7 detection vs 0.5-2.9%
# elsewhere). The other 12 split homogeneous ventricular myocardium along a
# sequencing-depth gradient - their Myl2 detection rate tracks median nCount
# almost monotonically (101 counts -> 76% detected, 194 counts -> 90%). This
# section realistically holds 3-4 domains, so the useful range sits well below
# the Seurat default. Reuses the existing graph, so the sweep is cheap.
banksy.res <- c(0.02, 0.05, 0.1, 0.2, 0.3)
object <- FindClusters(object, resolution = banksy.res, verbose = FALSE)
n.domains <- sapply(banksy.res, function(r)
  length(unique(object[[paste0("BANKSY_snn_res.", r)]][, 1])))
print(data.frame(resolution = banksy.res, n_domains = n.domains))

# Detection rate, not median: at 8um most bins are zero for any given gene, so a
# median of 0 is uninformative and reads as "absent" when the gene is merely
# sparse - that is what hid the atrium in the first pass. Atrium is Myl7+/Myl2-,
# ventricle Myl2+/Myl7-, Nppa marks trabecular / subendocardial myocardium.
# nCount is carried alongside to show whether a domain is a tissue or an artefact.
DefaultAssay(object) <- "Spatial.008um"
marker.df <- FetchData(object, vars = c("Myl7", "Myl2", "Nppa",
                                        "nCount_Spatial.008um"))

domain.summary <- do.call(rbind, lapply(banksy.res, function(r) {
  d <- marker.df
  d$domain <- object[[paste0("BANKSY_snn_res.", r)]][, 1]
  out <- aggregate(cbind(Myl7, Myl2, Nppa) ~ domain, d,
                   function(x) round(100 * mean(x > 0), 1))
  out$median_nCount <- aggregate(nCount_Spatial.008um ~ domain, d, median)[, 2]
  out$n_bins <- as.vector(table(d$domain))
  cbind(resolution = r, out)
}))
print(domain.summary)
write.csv(domain.summary, 'test.figures/banksy_resolution_sweep_summary.csv',
          row.names = FALSE)

# The domain map decides the resolution, not the domain count - one panel each.
res.plots <- lapply(seq_along(banksy.res), function(i) {
  col <- paste0("BANKSY_snn_res.", banksy.res[i])
  SpatialDimPlot(object, images = "slice1.008um", group.by = col, label = FALSE) +
    ggtitle(paste0("res = ", banksy.res[i], "  (", n.domains[i], " domains)")) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 10))
})
pdf('test.figures/banksy_resolution_sweep_008um.pdf', width = 15, height = 10)
print(wrap_plots(res.plots, ncol = 3))
dev.off()

# Working resolution - revisit once the sweep above is inspected
object$banksy_cluster <- object[["BANKSY_snn_res.0.2"]][, 1]
Idents(object) <- "banksy_cluster"
p <- SpatialDimPlot(object, images = "slice1.008um", group.by = "banksy_cluster", label = T, repel = T, label.size = 4) 
pdf('test.figures/banksy_spatial_008um.pdf')
print(p)
dev.off()

banksy_cells <- CellsByIdentities(object)
p <- SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00","grey50"),facet.highlight = T, combine=T) + NoLegend()
pdf('test.figures/banksy_spatial_highlight_008um.pdf')
print(p)
dev.off()

# Integration with snRNA-seq
ref <- readRDS('../OCM/heart_seurat_object_SCT.rds')
ref$celltype <- Idents(ref)
ref <- subset(ref, subset = !celltype %in% c("Cardiomyocytes (stressed)", "Epicardial - Mesothelial cells"))
ref$celltype <- droplevels(factor(ref$celltype))


ref <- NormalizeData(ref)
ref <- FindVariableFeatures(ref); ref <- ScaleData(ref); ref <- RunPCA(ref)

DefaultAssay(object) <- "sketch"
anchors <- FindTransferAnchors(reference = ref, query = object,
                               reference.reduction = "pca", dims = 1:30)
preds <- TransferData(anchorset = anchors, refdata = ref$celltype, dims = 1:30)
object <- AddMetaData(object, preds)

# 1. Working partition - the sweep showed 0.2 is where atrium and background separate
object$banksy_cluster <- object[["BANKSY_snn_res.0.2"]][, 1]

# 2. Project labels from the 50k sketched bins to all 248k
DefaultAssay(object) <- "sketch"
object <- ProjectData(
  object            = object,
  assay             = "Spatial.008um",
  full.reduction    = "full.pca.sketch",
  sketched.assay    = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model        = "umap.sketch",
  dims              = 1:15,
  refdata           = list(full_celltype = "predicted.id")
)

# 3. Composition per domain
DefaultAssay(object) <- "Spatial.008um"
comp <- object[[]] %>%
  dplyr::filter(!is.na(full_celltype)) %>%
  dplyr::count(banksy_cluster, full_celltype) %>%
  dplyr::group_by(banksy_cluster) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

p <- ggplot(comp, aes(x = factor(banksy_cluster), y = prop, fill = full_celltype)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "BANKSY domain", y = "% of bins", fill = "Cell type",
       title = "Cell type composition of BANKSY domains (8um bins)") +
  theme_minimal()

pdf('test.figures/banksy_celltype_composition.pdf', width = 10, height = 6)
print(p)
dev.off()





##############
adult_localdir <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial/9w/outs"
adult_object <- Load10X_Spatial(data.dir = adult_localdir, bin.size = 16)

DefaultAssay(adult_object) <- "Spatial.016um"

vln.plot <- VlnPlot(adult_object, features = 'nCount_Spatial.016um', pt.size = 0) + 
              theme(axis.text = element_text(size = 4)) + NoLegend()

count.plot <- SpatialFeaturePlot(adult_object, features = 'nCount_Spatial.016um') + 
                theme(legend.position = "right")

# Note that many spots have very few counts, in-part due to low cellular density in certain tissue regions
pdf('test.figures/adult_spatial_016um_count_distribution.pdf')
vln.plot | count.plot
dev.off()

# Normalize
adult_object <- NormalizeData(adult_object)

# A few circular high-count hotspots (bubbles/debris) sit far above the tissue
# median and saturate the colour scale, so clip every feature plot at q02/q98.
lo <- "q02"; hi <- "q98"

# Ttn: pan-cardiomyocyte sarcomere gene - confirms myocardium, near-uniform
p1 <- SpatialFeaturePlot(adult_object, features = "Ttn",
                         min.cutoff = lo, max.cutoff = hi) +
        ggtitle("Ttn expression (16um)")

# Myl7: atrial cardiomyocytes - marks the atrial region against ventricle
p2 <- SpatialFeaturePlot(adult_object, features = "Myl7",
                         min.cutoff = lo, max.cutoff = hi) +
        ggtitle("Myl7 expression (16um)")

pdf('test.figures/adult_spatial_016um_cardiac_expression.pdf')
p1 | p2
dev.off()
###
aged_localdir <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial/78w/outs"
aged_object <- Load10X_Spatial(data.dir = aged_localdir, bin.size = 16)

DefaultAssay(aged_object) <- "Spatial.016um"

vln.plot <- VlnPlot(aged_object, features = 'nCount_Spatial.016um', pt.size = 0) + 
              theme(axis.text = element_text(size = 4)) + NoLegend()

count.plot <- SpatialFeaturePlot(aged_object, features = 'nCount_Spatial.016um') + 
                theme(legend.position = "right")

# Note that many spots have very few counts, in-part due to low cellular density in certain tissue regions
pdf('test.figures/aged_spatial_016um_count_distribution.pdf')
vln.plot | count.plot
dev.off()

# Normalize
aged_object <- NormalizeData(aged_object)

# A few circular high-count hotspots (bubbles/debris) sit far above the tissue
# median and saturate the colour scale, so clip every feature plot at q02/q98.
lo <- "q02"; hi <- "q98"

# Ttn: pan-cardiomyocyte sarcomere gene - confirms myocardium, near-uniform
DefaultAssay(aged_object) <- "Spatial.016um"
p1 <- SpatialFeaturePlot(aged_object, features = "Ttn",
                         min.cutoff = lo, max.cutoff = hi) +
        ggtitle("Ttn expression (16um)")

# Myl7: atrial cardiomyocytes - marks the atrial region against ventricle
p2 <- SpatialFeaturePlot(aged_object, features = "Myl7",
                         min.cutoff = lo, max.cutoff = hi) +
        ggtitle("Myl7 expression (16um)")

pdf('test.figures/aged_spatial_016um_expression.pdf')
p1 | p2
dev.off()

###############################################################################
# Adult (9w) vs aged (78w): X-linked genes called DE with age in the bulk
# heart RNA-seq. The question here is not whether they differ - that came from
# the bulk data - but whether the difference is spatially structured or uniform
# across the section, which is exactly what bulk cannot tell you.
#
# Caveats worth carrying into any interpretation of this panel:
#  - One section per age, so this is descriptive. No test is run and none would
#    be valid: the replicate unit is the animal, not the bin.
#  - LogNormalize corrects per-bin depth but not section-level differences in
#    RNA quality or ambient background, so absolute intensity is only
#    semi-quantitative between the two sections.
###############################################################################

chrX.de <- c('Shroom4', 'Tspan7', 'Sh3kbp1', 'Med14', 'Kctd12b', '2210013O21Rik')

age.objects <- list(adult = adult_object, aged = aged_object)
for (nm in names(age.objects)) DefaultAssay(age.objects[[nm]]) <- "Spatial.016um"

# Drop anything absent from either section rather than letting FetchData fail
# halfway through building the figure.
in.both <- Reduce(intersect, lapply(age.objects, rownames))
absent  <- setdiff(chrX.de, in.both)
if (length(absent)) message("Not in both sections, skipped: ", paste(absent, collapse = ", "))
chrX.de <- intersect(chrX.de, in.both)
if (!length(chrX.de)) stop("None of the requested genes are present in both sections")

# Fetch every gene once per section and reuse the vectors for the colour limits,
# the violins and the summary table. FetchData on a full 16um section is not
# cheap, and the earlier per-call version pulled each gene three times over.
# [, gene] would break on 2210013O21Rik-style names, which data.frame() mangles
# into X2210013O21Rik - index positionally instead.
expr <- lapply(age.objects, function(o) {
  d <- FetchData(o, vars = chrX.de)
  setNames(as.list(as.data.frame(d)), chrX.de)
})

# Per-object q02/q98 cutoffs would give each section its own colour scale, which
# is exactly what makes an adult/aged comparison unreadable. Pool the values
# across both sections and clip to one shared range per gene.
shared_limits <- function(gene, probs = c(0.02, 0.98)) {
  v <- unlist(lapply(expr, `[[`, gene), use.names = FALSE)
  lims <- as.numeric(quantile(v, probs, na.rm = TRUE))
  # Sparse genes can be zero at both quantiles, which collapses the scale.
  if (!is.finite(diff(lims)) || diff(lims) <= 0) lims <- range(v, na.rm = TRUE)
  if (diff(lims) <= 0) lims <- c(0, max(1, lims[2]))
  lims
}

# viridis rather than Seurat's Spectral default: sequential and perceptually
# uniform, so the same colour means the same value in both sections.
gene_spatial <- function(obj, gene, lims, title) {
  SpatialFeaturePlot(obj, features = gene,
                     min.cutoff = lims[1], max.cutoff = lims[2]) +
    scale_fill_viridis_c(option = "magma", limits = lims,
                         oob = scales::squish, name = NULL) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
}

# Two gene pairs per row. One gene per row would be a 13x27in strip at six
# genes; this keeps the adult/aged pair adjacent (the only comparison that
# matters) while staying a printable figure.
GENES_PER_ROW <- 2

plots <- unlist(lapply(chrX.de, function(g) {
  lims <- shared_limits(g)
  list(gene_spatial(age.objects$adult, g, lims, paste0(g, " - adult 9w")),
       gene_spatial(age.objects$aged,  g, lims, paste0(g, " - aged 78w")))
}), recursive = FALSE)

n.rows <- ceiling(length(chrX.de) / GENES_PER_ROW)
panel <- wrap_plots(plots, ncol = 2 * GENES_PER_ROW) +
  plot_annotation(title = "X-linked age-DE genes from bulk heart RNA-seq (Visium HD, 16um bins)",
                  subtitle = "Colour scale shared between sections within each gene; n = 1 section per age",
                  theme = theme(plot.title = element_text(face = "bold")))

ggsave('test.figures/adult_vs_aged_chrX_DE.pdf', panel,
       width = 3.4 * 2 * GENES_PER_ROW, height = 3.6 * n.rows, limitsize = FALSE)
ggsave('test.figures/adult_vs_aged_chrX_DE.png', panel, dpi = 300, bg = "white",
       width = 3.4 * 2 * GENES_PER_ROW, height = 3.6 * n.rows, limitsize = FALSE)

# Distributions moved out of the spatial figure and faceted - at six genes a
# violin wedged beside each tissue pair is too small to read, and facetting also
# puts all six side by side, which is the more useful view.
# Detected bins only: most 16um bins are zero for these genes, so a violin over
# all bins is a spike at zero. The zero fraction is in the summary table.
viol.df <- do.call(rbind, lapply(chrX.de, function(g) {
  do.call(rbind, lapply(names(expr), function(nm) {
    data.frame(gene = g, age = nm, value = expr[[nm]][[g]])
  }))
}))
viol.df$gene <- factor(viol.df$gene, levels = chrX.de)
viol.df$age  <- factor(viol.df$age, levels = names(expr))

violins <- ggplot(viol.df[viol.df$value > 0, ], aes(age, value, fill = age)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.8) +
  stat_summary(fun = median, geom = "point", size = 2, colour = "black") +
  facet_wrap(~ gene, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c(adult = "#4C72B0", aged = "#C44E52")) +
  labs(title = "Bin-level expression, detected bins only",
       x = NULL, y = "log-normalised expression") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave('test.figures/adult_vs_aged_chrX_DE_violin.pdf', violins, width = 10, height = 6)
ggsave('test.figures/adult_vs_aged_chrX_DE_violin.png', violins, dpi = 300, bg = "white",
       width = 10, height = 6)

# The numbers behind the panel. Detection rate and mean over detected bins are
# reported separately because a shift in either alone reads very differently:
# fewer bins expressing is not the same result as the same bins expressing less.
chrX.summary <- do.call(rbind, lapply(chrX.de, function(g) {
  do.call(rbind, lapply(names(age.objects), function(nm) {
    v   <- expr[[nm]][[g]]
    raw <- GetAssayData(age.objects[[nm]], assay = "Spatial.016um", layer = "counts")[g, ]
    data.frame(
      gene          = g,
      age           = nm,
      n_bins        = length(v),
      total_counts  = sum(raw),
      bins_detected = sum(v > 0),
      pct_detected  = round(100 * mean(v > 0), 2),
      mean_all      = round(mean(v), 4),
      mean_detected = round(if (any(v > 0)) mean(v[v > 0]) else NA_real_, 4)
    )
  }))
}))
write.csv(chrX.summary, 'test.figures/adult_vs_aged_chrX_DE_summary.csv', row.names = FALSE)
print(chrX.summary)

# At six genes some will be too sparse at 16um for the tissue plot to say
# anything - flag them rather than letting a near-empty panel read as "absent".
too.sparse <- unique(chrX.summary$gene[chrX.summary$pct_detected < 1])
if (length(too.sparse)) {
  message("Detected in <1% of bins in at least one section - tissue panel is not ",
          "informative for: ", paste(too.sparse, collapse = ", "))
}

