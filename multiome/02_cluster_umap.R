# ---------------------------------------------------------------------------
# 02 - Clustering, cell typing and UMAPs for the multiome heart data.
#
# Figures for presentation, and the cluster/celltype assignment that the
# per-CELLTYPE allelic route needs later (sinto wants a barcode -> group map).
#
# Depends on 01_qc_barcodes.R having run: it reads qc/qc_pass_<id>.csv and
# scores only nuclei that passed ddqcR plus the ATAC floors.
#
# NO SIGNAC DEPENDENCY, deliberately - it appears nowhere in this project, so
# assuming it would be a new install on a deadline. The ATAC side instead
# reuses the LSI that cellranger-arc already computed
# (outs/analysis/dimensionality_reduction/atac/lsa_projection.csv), attached as
# a Seurat DimReduc. FindMultiModalNeighbors is Seurat, not Signac, so WNN
# works from that.
#
# WHY THE MERGED UMAP IS RNA-ONLY AND WNN IS PER-SAMPLE. cellranger computed
# each sample's LSI independently, so 9w component 3 and 78w component 3 are
# not the same axis and the two matrices cannot be concatenated into a shared
# space. Doing it anyway would produce a joint embedding whose structure is an
# artefact. So:
#   merged, both samples : RNA only  -> cell types across age
#   per sample           : RNA + ATAC WNN -> the multiome view
# Recomputing LSI over a common peak set would fix this, and needs Signac.
#
# A NOTE FOR THE TALK. There is one animal per age, so any 9w-vs-78w
# difference in these plots is confounded with animal and with data quality -
# 78w has 2.2x the nuclei and better ATAC. The per-cluster sample composition
# table below quantifies how much of the structure is sample-driven. Cell type
# identity is safe to show; age comparisons are descriptive only.
#
# Run in seurat_env:
#   conda activate seurat_env
#   Rscript multiome/02_cluster_umap.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(readr)
  library(tibble); library(patchwork)
})
options(future.globals.maxSize = 16 * 1024^3)
set.seed(1)

BASE <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2"
WORK <- file.path(BASE, "adult_aged_multiome")
QC   <- file.path(WORK, "qc")
OUT  <- file.path(WORK, "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SAMPLES    <- c("9w", "78w")
RESOLUTION <- 0.1     # 0.1 is what Seurat_preprocessing.R used to land 11 clusters
N_PCS      <- 30
LSI_DIMS   <- 2:15    # component 1 of an ATAC LSI tracks depth; conventionally dropped
MARKER_FILE <- file.path(BASE, "adult_aged_heart_snRNAseq", "marker.genes.txt")

# Provisional scoring panels. The LABELS are the ones already used in
# OCM_heart/Seurat_preprocessing.R so the multiome figures line up with the
# earlier snRNA ones. Assignment is by highest mean module score per cluster
# and is PROVISIONAL - check it against the dotplot before believing it.
PANELS <- list(
  "Ventricular Cardiomyocytes"      = c("Tnnt2","Myh6","Actc1","Ttn","Ryr2","Pln","Myl2","Tnni3","Actn2","Trdn"),
  "Cardiomyocytes (stressed)"       = c("Nppa","Nppb","Myh7","Ankrd1","Acta1","Xirp2"),
  "Fibroblasts"                     = c("Col1a1","Col3a1","Dcn","Pdgfra","Gsn","Postn","Lum"),
  "Endothelial cells"               = c("Pecam1","Cdh5","Egfl7","Kdr","Emcn","Ly6c1","Plvap","Tie1"),
  "Endocardium"                     = c("Npr3","Nfatc1","Vwf","Cytl1"),
  "Lymphatic endothelial"           = c("Lyve1","Prox1","Mmrn1","Flt4","Ccl21a"),
  "Pericytes - Smooth muscle cells" = c("Rgs5","Abcc9","Kcnj8","Myh11","Acta2","Tagln","Notch3"),
  "Macrophages"                     = c("Ptprc","Cd68","Adgre1","Csf1r","Lyz2","Mrc1","Fcgr1"),
  "T cells"                         = c("Cd3e","Cd3d","Cd8a","Themis","Skap1"),
  "B cells"                         = c("Cd79a","Cd79b","Ms4a1","Ighm"),
  "Epicardial - Mesothelial cells"  = c("Wt1","Msln","Upk3b","Krt19")
)

say <- function(...) cat(sprintf(...), "\n", sep = "")
pdfout <- function(name, w = 9, h = 7) pdf(file.path(OUT, name), width = w, height = h)

# ---------------------------------------------------------------------------
# load QC-passing nuclei
# ---------------------------------------------------------------------------
load_one <- function(id) {
  passf <- file.path(QC, sprintf("qc_pass_%s.csv", id))
  h5    <- file.path(WORK, id, "outs", "filtered_feature_bc_matrix.h5")
  stopifnot(file.exists(passf), file.exists(h5))
  pass <- read_csv(passf, show_col_types = FALSE)

  mat <- Read10X_h5(h5)
  if (is.list(mat)) mat <- mat[["Gene Expression"]]
  keep <- intersect(colnames(mat), pass$barcode)
  say("%s: %d QC-pass barcodes, %d found in matrix", id, nrow(pass), length(keep))

  o <- CreateSeuratObject(mat[, keep], project = id)
  o$sample     <- id
  o$percent.mt <- PercentageFeatureSet(o, pattern = "^mt-")
  meta <- pass %>% filter(barcode %in% keep) %>%
    select(barcode, atac_fragments, frip, tss_frac) %>%
    column_to_rownames("barcode")
  o <- AddMetaData(o, meta[colnames(o), , drop = FALSE])
  o
}
objs <- lapply(SAMPLES, load_one); names(objs) <- SAMPLES

# ---------------------------------------------------------------------------
# merged RNA analysis
# ---------------------------------------------------------------------------
say("--- merged RNA: SCTransform, PCA, UMAP, clustering ---")
merged <- merge(objs[[1]], y = objs[[2]], add.cell.ids = SAMPLES)
# Plain SCTransform with no vars.to.regress, matching the merged-object call in
# Seurat_preprocessing.R (its per-sample calls that regress percent.mt are all
# commented out). Not regressing is the right choice here rather than just the
# consistent one: median percent.mt is ~10% and is strongly cell-type-linked -
# cardiomyocytes are the most mitochondria-rich cell type in the heart - so
# regressing it out would partly regress out cardiomyocyte identity.
merged <- SCTransform(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = N_PCS, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:N_PCS, verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:N_PCS, verbose = FALSE)
merged <- FindClusters(merged, resolution = RESOLUTION, verbose = FALSE)
say("clusters at resolution %.2f: %d", RESOLUTION, length(levels(merged)))

# ---- provisional cell typing ----
present <- lapply(PANELS, function(g) intersect(g, rownames(merged)))
drop <- names(present)[lengths(present) == 0]
if (length(drop)) say("panels with no genes present, skipped: %s", paste(drop, collapse = ", "))
present <- present[lengths(present) > 0]
merged <- AddModuleScore(merged, features = present, name = "panel", verbose = FALSE)
score_cols <- paste0("panel", seq_along(present))
colnames(merged@meta.data)[match(score_cols, colnames(merged@meta.data))] <- names(present)

per_cluster <- merged@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(names(present)), mean), .groups = "drop")
assign <- names(present)[apply(per_cluster[, names(present)], 1, which.max)]
names(assign) <- as.character(per_cluster$seurat_clusters)
merged$celltype_provisional <- assign[as.character(merged$seurat_clusters)]
say("provisional assignment:")
print(as.data.frame(table(cluster = merged$seurat_clusters,
                          celltype = merged$celltype_provisional)) %>% filter(Freq > 0))
write_csv(per_cluster, file.path(OUT, "cluster_panel_scores.csv"))

# ---- the honesty check: is a cluster just one sample? ----
comp <- merged@meta.data %>%
  count(seurat_clusters, sample) %>%
  tidyr::pivot_wider(names_from = sample, values_from = n, values_fill = 0) %>%
  mutate(total = rowSums(across(all_of(SAMPLES))),
         pct_9w = round(100 * .data[["9w"]] / total, 1))
say("per-cluster sample composition (a cluster near 0%% or 100%% 9w is sample-driven):")
print(as.data.frame(comp), row.names = FALSE)
write_csv(comp, file.path(OUT, "cluster_sample_composition.csv"))

# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
pdfout("umap_rna_clusters.pdf", 10, 8)
print(DimPlot(merged, label = TRUE, repel = TRUE) + ggtitle("Merged RNA, clusters"))
print(DimPlot(merged, group.by = "celltype_provisional", label = TRUE, repel = TRUE,
              label.size = 3) +
        ggtitle("Merged RNA, provisional cell type") +
        theme(legend.position = "bottom"))
print(DimPlot(merged, group.by = "sample") + ggtitle("Merged RNA, by sample"))
print(DimPlot(merged, split.by = "sample", group.by = "celltype_provisional") +
        ggtitle("Split by sample - n=1 per age, so read cautiously"))
dev.off()

pdfout("umap_rna_qc.pdf", 11, 8)
print(FeaturePlot(merged, features = c("nCount_RNA","nFeature_RNA","percent.mt",
                                       "atac_fragments","frip","tss_frac")))
dev.off()

if (file.exists(MARKER_FILE)) {
  mk <- read.table(MARKER_FILE, header = FALSE)$V1
  mk <- unique(mk[mk %in% rownames(merged)])
  say("marker.genes.txt: %d of %d genes present", length(mk),
      length(read.table(MARKER_FILE, header = FALSE)$V1))
  pdfout("marker_dotplot.pdf", 14, 7)
  print(DotPlot(merged, features = mk, group.by = "celltype_provisional") +
          RotatedAxis() + ggtitle("Established marker panel (marker.genes.txt)"))
  dev.off()
} else say("marker.genes.txt not found at %s - skipping dotplot", MARKER_FILE)

pdfout("marker_featureplots.pdf", 12, 9)
key <- intersect(c("Tnnt2","Myh6","Nppa","Col1a1","Pecam1","Rgs5","Ptprc","Msln","Xist"),
                 rownames(merged))
print(FeaturePlot(merged, features = key, ncol = 3))
dev.off()

say("--- FindAllMarkers ---")
DefaultAssay(merged) <- "SCT"
merged <- PrepSCTFindMarkers(merged)
am <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25,
                     logfc.threshold = 0.5, verbose = FALSE)
write_csv(am, file.path(OUT, "cluster_markers.csv"))
say("wrote %d marker rows", nrow(am))

# ---------------------------------------------------------------------------
# per-sample WNN, using cellranger's ATAC LSI
# ---------------------------------------------------------------------------
wnn_one <- function(id) {
  o <- objs[[id]]
  lsif <- file.path(WORK, id, "outs", "analysis", "dimensionality_reduction",
                    "atac", "lsa_projection.csv")
  if (!file.exists(lsif)) { say("%s: no LSI projection, skipping WNN", id); return(invisible(NULL)) }

  lsi <- read.csv(lsif, row.names = 1, check.names = FALSE)
  common <- intersect(colnames(o), rownames(lsi))
  say("%s WNN: %d nuclei with both RNA and ATAC LSI", id, length(common))
  if (length(common) < 50) { say("%s: too few, skipping WNN", id); return(invisible(NULL)) }
  o <- subset(o, cells = common)

  o <- SCTransform(o, verbose = FALSE)   # as above: percent.mt not regressed
  o <- RunPCA(o, npcs = N_PCS, verbose = FALSE)

  emb <- as.matrix(lsi[common, , drop = FALSE])
  colnames(emb) <- paste0("lsi_", seq_len(ncol(emb)))
  o[["lsi"]] <- CreateDimReducObject(embeddings = emb, key = "lsi_", assay = "SCT")

  dims_atac <- intersect(LSI_DIMS, seq_len(ncol(emb)))
  o <- FindMultiModalNeighbors(o, reduction.list = list("pca", "lsi"),
                               dims.list = list(1:N_PCS, dims_atac), verbose = FALSE)
  o <- RunUMAP(o, nn.name = "weighted.nn", reduction.name = "wnn.umap",
               reduction.key = "wnnUMAP_", verbose = FALSE)
  o <- FindClusters(o, graph.name = "wsnn", resolution = RESOLUTION, verbose = FALSE)
  o <- RunUMAP(o, reduction = "lsi", dims = dims_atac,
               reduction.name = "atac.umap", reduction.key = "atacUMAP_", verbose = FALSE)
  o <- RunUMAP(o, dims = 1:N_PCS, reduction.name = "rna.umap",
               reduction.key = "rnaUMAP_", verbose = FALSE)

  pdfout(sprintf("umap_wnn_%s.pdf", id), 13, 5)
  print(
    (DimPlot(o, reduction = "rna.umap",  label = TRUE) + ggtitle(sprintf("%s RNA", id))) |
    (DimPlot(o, reduction = "atac.umap", label = TRUE) + ggtitle(sprintf("%s ATAC (cellranger LSI)", id))) |
    (DimPlot(o, reduction = "wnn.umap",  label = TRUE) + ggtitle(sprintf("%s WNN joint", id)))
  )
  dev.off()
  saveRDS(o, file.path(WORK, sprintf("multiome_%s_wnn.rds", id)))
  invisible(NULL)
}
say("--- per-sample WNN ---")
for (id in SAMPLES) try(wnn_one(id))

# ---------------------------------------------------------------------------
# save, and export the group map the per-celltype allelic route needs
# ---------------------------------------------------------------------------
saveRDS(merged, file.path(WORK, "multiome_merged_rna.rds"))

exp <- merged@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(barcode = sub("^(9w|78w)_", "", cell)) %>%
  select(sample, barcode, seurat_clusters, celltype_provisional)
write_csv(exp, file.path(OUT, "nucleus_celltype.csv"))

# sinto group files, one per sample, grouped by provisional cell type. Uses the
# same two-barcode-space discipline as 01: the group label is what names the
# output BAM, so RNA and ATAC agree.
for (id in SAMPLES) {
  pass <- read_csv(file.path(QC, sprintf("qc_pass_%s.csv", id)), show_col_types = FALSE)
  e <- exp %>% filter(sample == id) %>%
    inner_join(pass %>% select(barcode, gex_barcode, atac_barcode), by = "barcode") %>%
    mutate(group = gsub("[^A-Za-z0-9]+", "_", celltype_provisional))
  write.table(e[, c("gex_barcode","group")],
              file.path(OUT, sprintf("sinto_%s_rna_bycelltype.txt", id)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(e[, c("atac_barcode","group")],
              file.path(OUT, sprintf("sinto_%s_atac_bycelltype.txt", id)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  say("%s: wrote per-celltype sinto files (%d nuclei, %d groups)",
      id, nrow(e), length(unique(e$group)))
}
say("")
say("figures and tables under %s", OUT)
say("PROVISIONAL cell types - verify against marker_dotplot.pdf before presenting.")
