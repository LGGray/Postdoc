# ---------------------------------------------------------------------------
# 04 - Regenerate figures from the 03 checkpoint.
#
# 03 writes multiome_signac_joint.rds once the merged LSI, clustering and cell
# typing are done - which is the multi-hour part. This reloads that and redraws
# the figures only, so iterating on presentation details costs minutes instead
# of rerunning FeatureMatrix over both fragment files.
#
# Deliberately a separate script rather than a resume branch inside 03: 03 is
# working and this needs no changes to it.
#
#   conda activate seurat_env
#   Rscript multiome/04_figures_from_checkpoint.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(ggplot2)
  library(dplyr); library(readr); library(tibble); library(patchwork)
})
options(future.globals.maxSize = 32 * 1024^3)

BASE <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2"
WORK <- file.path(BASE, "adult_aged_multiome")
OUT  <- file.path(WORK, "figures_signac")
RDS  <- file.path(WORK, "multiome_signac_joint.rds")
source(file.path(BASE, "Postdoc", "multiome", "00_helpers.R"))
say <- function(...) cat(sprintf(...), "\n", sep = "")

ESCAPE_GENES <- c("Kdm5c","Kdm6a","Ddx3x","Eif2s3x","Ftx","Jpx",
                  "Pbdc1","Utp14a","Akap17a","Sts")

if (!file.exists(RDS)) stop("no checkpoint at ", RDS, " - run 03 first")
say("loading %s", RDS)
obj <- readRDS(RDS)
say("%d nuclei, %d clusters, assays: %s", ncol(obj),
    length(unique(obj$seurat_clusters)), paste(Assays(obj), collapse = ", "))

obj <- set_meta(obj, "celltype_short", short_labels(obj$celltype_provisional))
CT_LEVELS <- sort(unique(obj$celltype_short))
SC_COL <- celltype_scale(CT_LEVELS, "colour")
say("cell types (%d): %s", length(CT_LEVELS), paste(CT_LEVELS, collapse = ", "))
Idents(obj) <- "celltype_short"

# ---- three panels, one shared legend ----
panel <- function(red, ttl) {
  p <- DimPlot(obj, reduction = red, group.by = "celltype_short",
               pt.size = 0.3, shuffle = TRUE) +
    ggtitle(ttl) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.title = element_text(size = 9))
  if (!is.null(SC_COL)) p <- p + SC_COL
  p
}
dev_open(file.path(OUT, "umap_joint.pdf"), width = 15, height = 5.5)
print(
  (panel("rna.umap", "RNA") | panel("atac.umap", "ATAC (common peaks, merged LSI)") |
   panel("wnn.umap", "WNN joint")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 9))
)
dev.off()
say("wrote umap_joint.pdf")

# ---- the panel to present: one large plot, direct labels ----
dev_open(file.path(OUT, "umap_wnn_labelled.pdf"), width = 9, height = 7.5)
p <- DimPlot(obj, reduction = "wnn.umap", group.by = "celltype_short",
             label = TRUE, repel = TRUE, label.size = 4, pt.size = 0.4,
             shuffle = TRUE) +
  NoLegend() + ggtitle("Joint RNA + ATAC (WNN), 9w and 78w") +
  theme(plot.title = element_text(size = 13, face = "bold"))
if (!is.null(SC_COL)) p <- p + SC_COL
print(p)
dev.off()
say("wrote umap_wnn_labelled.pdf")

# ---- by sample ----
dev_open(file.path(OUT, "umap_joint_by_sample.pdf"), width = 12, height = 5)
p1 <- DimPlot(obj, reduction = "wnn.umap", group.by = "sample", pt.size = 0.3,
              shuffle = TRUE) + ggtitle("WNN, by sample")
p2 <- DimPlot(obj, reduction = "wnn.umap", split.by = "sample",
              group.by = "celltype_short", pt.size = 0.3) +
  ggtitle("n=1 per age - descriptive only")
if (!is.null(SC_COL)) p2 <- p2 + SC_COL
print(p1 | p2)
dev.off()
say("wrote umap_joint_by_sample.pdf")

# ---- WNN modality weights ----
wcols <- grep("\\.weight$", colnames(obj@meta.data), value = TRUE)
if (length(wcols)) {
  dev_open(file.path(OUT, "wnn_modality_weights.pdf"), width = 9, height = 6)
  print(VlnPlot(obj, features = wcols, group.by = "celltype_short",
                pt.size = 0, ncol = 1) &
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                axis.title.x = element_blank()))
  dev.off()
  say("wrote wnn_modality_weights.pdf (%s)", paste(wcols, collapse = ", "))
}

# ---- escape gene activity, if 03 built the assay ----
if ("GeneActivity" %in% Assays(obj)) {
  DefaultAssay(obj) <- "GeneActivity"
  esc <- intersect(ESCAPE_GENES, rownames(obj[["GeneActivity"]]))
  if (length(esc)) {
    dev_open(file.path(OUT, "escape_gene_activity.pdf"), width = 12, height = 8)
    print(FeaturePlot(obj, features = esc, reduction = "wnn.umap", ncol = 4) &
            theme(plot.title = element_text(size = 9)))
    print(DotPlot(obj, features = esc, group.by = "celltype_short") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            ggtitle("Core escape genes: ATAC gene activity"))
    dev.off()
    say("wrote escape_gene_activity.pdf (%d genes)", length(esc))
  }
} else say("no GeneActivity assay in the checkpoint - skipping")

# ---- coverage tracks ----
DefaultAssay(obj) <- "ATAC"
ann <- Annotation(obj)
avail <- if (!is.null(ann)) intersect(c("Xist", ESCAPE_GENES), ann$gene_name) else character(0)
if (length(avail)) {
  dev_open(file.path(OUT, "coverage_escape_loci.pdf"), width = 10, height = 7)
  for (g in avail) {
    pp <- try(CoveragePlot(obj, region = g, features = g, expression.assay = "SCT",
                           extend.upstream = 5000, extend.downstream = 5000), silent = TRUE)
    if (!inherits(pp, "try-error")) print(pp + patchwork::plot_annotation(title = g))
    else say("  CoveragePlot failed for %s", g)
  }
  dev.off()
  say("wrote coverage_escape_loci.pdf (%d loci)", length(avail))
} else say("no annotation on the object - skipping coverage tracks")

say("")
say("figures refreshed under %s", OUT)
