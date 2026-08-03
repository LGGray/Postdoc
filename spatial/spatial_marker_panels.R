# Representative marker panels for the Visium HD cardiac sections.
# Spatial + UMAP plots of endothelial / fibroblast / other cell-type scores,
# plus a detection check for Mex3a. Written for the interview presentation:
# a small number of clean panels, not a full analysis.
#
# Run on the cluster from adult_aged_spatial/, e.g.
#   Rscript spatial_marker_panels.R
#
# Uses spaceranger's own UMAP and graph clustering rather than recomputing,
# so the panels match what the web summary shows.

library(Seurat)
library(ggplot2)
library(Matrix)

##### ---------------------- CONFIG ---------------------- #####
BASE      <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial"
SAMPLE    <- "9w"      # adult = baseline; better complexity than 78w
BIN       <- "square_016um"
OUT_DIR   <- file.path(BASE, "figures", paste0(SAMPLE, "_", BIN))

# Point size for the spatial plots. Tune if bins look too sparse or too merged:
# smaller bins (008um) need a smaller value, 032um a larger one.
POINT_SIZE <- 0.35

# Marker sets. Endothelial and fibroblast are the collaborator's; the rest are
# standard mouse heart markers for context.
MARKERS <- list(
  Endothelial   = c("Cdh5", "Pecam1", "Flt1", "Vwf"),
  Fibroblast    = c("Pdgfra", "Col1a1", "Dcn"),
  Cardiomyocyte = c("Tnnt2", "Myh6", "Actc1", "Ttn"),
  VSMC          = c("Myh11", "Tagln", "Acta2"),
  Myeloid       = c("Ptprc", "Csf1r", "Adgre1", "Cd68"),
  Lymphoid      = c("Cd3e", "Cd79a", "Ms4a1")
)

# A module score built from one surviving marker is just that gene's expression
# under a cell-type label, which is a claim the panel can't support. Require at
# least this many detected markers before scoring a set.
MIN_MARKERS <- 2

# Endothelial subtypes - only plotted if they look detectable.
EC_SUBTYPES <- c("Bmx", "Mecom", "Gja5", "Rgcc", "Ackr1", "Fabp4")

# Single genes worth showing on their own.
SINGLE_GENES <- c("Pecam1", "Cdh5", "Col1a1", "Pdgfra", "Gsn", "Myl2", "Nppa")

# The gene in question, plus reference genes spanning a range of expression
# levels so "not detected" can be interpreted against a known scale.
GENE_OF_INTEREST <- "Mex3a"
REFERENCE_GENES  <- c("Tnnt2", "Pecam1", "Col1a1", "Cd3e", "Mecom", "Gja5")
##### ------------------------------------------------------ #####

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

bin_dir <- file.path(BASE, SAMPLE, "outs", "binned_outputs", BIN)
if (!dir.exists(bin_dir)) stop("No such bin directory: ", bin_dir)

# ---------------------------------------------------------------- load

message("Reading counts ...")
# If Read10X_h5 errors on the compression, repack first (see Seurat_preprocessing.R):
#   ptrepack --complevel 5 in.h5:/matrix out.h5:/matrix
counts <- Read10X_h5(file.path(bin_dir, "filtered_feature_bc_matrix.h5"))

# Spatial coordinates. spaceranger 3.x+ writes parquet; older versions csv.
read_positions <- function(dir) {
  pq <- file.path(dir, "spatial", "tissue_positions.parquet")
  cs <- file.path(dir, "spatial", "tissue_positions.csv")
  if (file.exists(pq)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("tissue_positions.parquet needs the 'arrow' package: install.packages('arrow')")
    }
    as.data.frame(arrow::read_parquet(pq))
  } else if (file.exists(cs)) {
    read.csv(cs)
  } else {
    stop("No tissue_positions file under ", file.path(dir, "spatial"))
  }
}
pos <- read_positions(bin_dir)

# spaceranger's precomputed secondary analysis.
umap_file    <- file.path(bin_dir, "analysis", "umap",
                          "gene_expression_2_components", "projection.csv")
cluster_file <- file.path(bin_dir, "analysis", "clustering",
                          "gene_expression_graphclust", "clusters.csv")
for (f in c(umap_file, cluster_file)) {
  if (!file.exists(f)) stop("Missing secondary analysis output: ", f)
}
umap     <- read.csv(umap_file)
clusters <- read.csv(cluster_file)

# tissue_positions uses "barcode", the analysis CSVs use "Barcode", and read.csv
# turns "UMAP-1" into "UMAP.1". Normalise all of it up front.
norm_bc <- function(df) {
  i <- match("barcode", tolower(names(df)))
  if (is.na(i)) stop("No barcode column in: ", paste(names(df), collapse = ", "))
  names(df)[i] <- "barcode"
  df
}
pos      <- norm_bc(pos)
umap     <- norm_bc(umap)
clusters <- norm_bc(clusters)
# Don't hardcode the coordinate/cluster column names - they vary across
# spaceranger versions ("UMAP-1" vs "UMAP.1" after read.csv, "Cluster" vs
# "graphclust"). Once barcode is normalised, take what's left positionally.
umap_cols <- setdiff(names(umap), "barcode")
clus_cols <- setdiff(names(clusters), "barcode")
if (length(umap_cols) < 2) {
  stop("UMAP projection has no coordinate columns: ", paste(names(umap), collapse = ", "))
}
if (length(clus_cols) < 1) {
  stop("clusters.csv has no cluster column: ", paste(names(clusters), collapse = ", "))
}
message("UMAP columns: ", paste(umap_cols[1:2], collapse = ", "),
        " | cluster column: ", clus_cols[1])

# ---- QC (3' WTA: mito/hemo/ribo are present, unlike probe-based) ----
obj <- CreateSeuratObject(counts, project = SAMPLE, assay = "Spatial",
                          min.cells = 3, min.features = 10)

obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^mt-")
obj[["percent.hb"]]   <- PercentageFeatureSet(obj, pattern = "^Hb[ab]-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")

message("nCount  quantiles: ", paste(round(quantile(obj$nCount_Spatial,  c(.01,.05,.25,.5,.95))), collapse=" "))
message("nFeature quantiles: ", paste(round(quantile(obj$nFeature_Spatial, c(.01,.05,.25,.5,.95))), collapse=" "))
message("percent.mt median: ", round(median(obj$percent.mt), 1))

MIN_COUNT <- 100   # 16um; drop to ~15 for 008um
MIN_FEAT  <- 50
MAX_MT    <- 50   # loose on purpose - cardiomyocytes are legitimately 30-40%

keep <- obj$nCount_Spatial >= MIN_COUNT & obj$nFeature_Spatial >= MIN_FEAT & obj$percent.mt <= MAX_MT
message("Keeping ", sum(keep), "/", ncol(obj), " bins (", round(100*mean(keep),1), "%)")
obj <- obj[, keep]

# Keep a log-normalised Spatial assay alongside SCT. The detection check and the
# gene-of-interest panel have to read the raw scale: SCT stores *corrected*
# counts, and its variance stabilisation is exactly the wrong lens when the
# question is "is this gene seen at all".
obj <- NormalizeData(obj, assay = "Spatial", verbose = FALSE)

# SCTransform for the module scores and marker panels - better behaved than
# LogNormalize when bins are shallow (median ~600 UMI here). vst.flavor "v2"
# needs glmGamPoi; fall back to v1 rather than erroring on a fresh cluster env.
vst_flavor <- if (requireNamespace("glmGamPoi", quietly = TRUE)) "v2" else NULL
if (is.null(vst_flavor)) {
  message("glmGamPoi not found - using sctransform v1 (slower). ",
          "BiocManager::install('glmGamPoi') to speed this up.")
}
message("SCTransform ...")
obj <- SCTransform(obj, assay = "Spatial", new.assay.name = "SCT",
                   vst.flavor = vst_flavor, verbose = FALSE)
DefaultAssay(obj) <- "SCT"


bc <- colnames(obj)
idx <- function(df) match(bc, df$barcode)

# pxl_row is the vertical axis and pxl_col the horizontal one. Image coordinates
# run downwards from the top-left, so y is reversed at plotting time.
meta <- data.frame(
  x       = pos$pxl_col_in_fullres[idx(pos)],
  y       = pos$pxl_row_in_fullres[idx(pos)],
  UMAP_1  = umap[[umap_cols[1]]][idx(umap)],
  UMAP_2  = umap[[umap_cols[2]]][idx(umap)],
  cluster = factor(clusters[[clus_cols[1]]][idx(clusters)]),
  row.names = bc
)
# Check each column explicitly. A NULL column is dropped silently by data.frame,
# so testing meta$UMAP_1 alone can pass on an object that has no UMAP at all.
for (cn in c("x", "y", "UMAP_1", "UMAP_2", "cluster")) {
  if (!cn %in% names(meta)) stop("meta is missing column: ", cn)
  n_na <- sum(is.na(meta[[cn]]))
  if (n_na) warning(n_na, "/", nrow(meta), " bins have no ", cn,
                    " - check bin size consistency between matrix and analysis outputs")
}
obj <- AddMetaData(obj, meta)

# ---------------------------------------------------------------- scores

# Genes now live in two places. The Spatial assay holds everything that survived
# CreateSeuratObject; SCT additionally drops genes too sparse to fit (min_cells =
# 5 by default). Keep the two reasons apart - "absent from the reference" and
# "detected, but too sparse for SCT to model" are very different statements to
# make about a gene.
in_spatial <- rownames(obj[["Spatial"]])
in_sct     <- rownames(obj[["SCT"]])

present <- function(g) intersect(g, in_sct)

report_missing <- function(g, nm) {
  absent <- setdiff(g, in_spatial)
  sparse <- setdiff(intersect(g, in_spatial), in_sct)
  if (length(absent)) message("  ", nm, ": not in reference - ", paste(absent, collapse = ", "))
  if (length(sparse)) message("  ", nm, ": too sparse for SCT - ", paste(sparse, collapse = ", "))
}

MIN_MARKERS <- 2
for (nm in names(MARKERS)) {
  found <- present(MARKERS[[nm]])
  report_missing(MARKERS[[nm]], nm)
  if (length(found) < MIN_MARKERS) {
    message("  ", nm, ": skipped, only ", length(found), " of ",
            length(MARKERS[[nm]]), " markers detected")
    next
  }
  obj <- AddModuleScore(obj, features = list(found), name = paste0(nm, "_"),
                        ctrl = 50, seed = 42)
  # AddModuleScore appends an index; rename to the plain cell-type name.
  names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "_1")] <- nm
}

# ---------------------------------------------------------------- plotting

# Shared look: no axes, fixed aspect, dark background so bright signal carries
# on a projector.
theme_panel <- function() {
  theme_void(base_size = 11) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.key.width = unit(0.3, "cm")
    )
}

# Draws high values last so signal is never buried under zero-valued bins.
spatial_plot <- function(obj, feature, title = feature, limits = NULL) {
  df <- data.frame(
    x = obj$x, y = obj$y,
    value = FetchData(obj, vars = feature)[, 1]
  )
  df <- df[order(df$value), ]
  ggplot(df, aes(x, y, colour = value)) +
    geom_point(size = POINT_SIZE, shape = 16) +
    scale_colour_viridis_c(option = "magma", name = NULL, limits = limits) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = title) +
    theme_panel()
}

umap_plot <- function(obj, feature, title = feature) {
  df <- data.frame(
    x = obj$UMAP_1, y = obj$UMAP_2,
    value = FetchData(obj, vars = feature)[, 1]
  )
  df <- df[order(df$value), ]
  ggplot(df, aes(x, y, colour = value)) +
    geom_point(size = POINT_SIZE, shape = 16) +
    scale_colour_viridis_c(option = "magma", name = NULL) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
}

cluster_plot <- function(obj) {
  df <- data.frame(x = obj$x, y = obj$y, cluster = obj$cluster)
  ggplot(df, aes(x, y, colour = cluster)) +
    geom_point(size = POINT_SIZE, shape = 16) +
    scale_y_reverse() +
    coord_fixed() +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(title = paste0(SAMPLE, " - graph-based clusters (", BIN, ")")) +
    theme_panel()
}

save_panel <- function(p, name, w = 6, h = 6) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
}

message("Plotting ...")

save_panel(cluster_plot(obj), "00_clusters", w = 7)

# The two panels the collaborator actually asked for: score on tissue, and the
# same score on the UMAP.
for (nm in names(MARKERS)) {
  if (!nm %in% names(obj@meta.data)) next
  save_panel(spatial_plot(obj, nm, paste0(nm, " score")), paste0("01_spatial_", nm))
  save_panel(umap_plot(obj, nm, paste0(nm, " score")),    paste0("02_umap_", nm))
}

for (g in present(SINGLE_GENES)) {
  save_panel(spatial_plot(obj, g), paste0("03_gene_", g))
}

for (g in present(EC_SUBTYPES)) {
  save_panel(spatial_plot(obj, g), paste0("04_ECsubtype_", g))
}

# Combined figure: endothelial and fibroblast, spatial beside UMAP.
# Guarded on both scores existing - a marker set skipped above would otherwise
# take the script down here, after every other panel has already been written.
combo_sets <- c("Endothelial", "Fibroblast")
have_combo <- all(combo_sets %in% names(obj@meta.data))
if (!have_combo) {
  message("Skipping combined figure - missing score(s): ",
          paste(setdiff(combo_sets, names(obj@meta.data)), collapse = ", "))
}
if (have_combo && requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combo <- (spatial_plot(obj, "Endothelial", "Endothelial score") |
            umap_plot(obj, "Endothelial", "Endothelial score")) /
           (spatial_plot(obj, "Fibroblast", "Fibroblast score") |
            umap_plot(obj, "Fibroblast", "Fibroblast score"))
  ggsave(file.path(OUT_DIR, "05_combined_EC_FB.png"), combo,
         width = 11, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, "05_combined_EC_FB.pdf"), combo, width = 11, height = 10)
}

# ---------------------------------------------------------------- detection

# Report raw counts, not normalized values - the question is whether the gene
# is seen at all, and normalization can make a handful of counts look like
# meaningful signal.
# assay = "Spatial" is load-bearing: DefaultAssay is SCT by this point, and SCT's
# counts layer holds corrected counts, not the observed UMIs.
raw <- GetAssayData(obj, assay = "Spatial", layer = "counts")
total_umi <- sum(raw)

detection_row <- function(g) {
  if (!g %in% rownames(raw)) {
    return(data.frame(gene = g, in_reference = FALSE, total_counts = NA,
                      bins_detected = NA, pct_bins = NA, per_million = NA))
  }
  v <- raw[g, ]
  data.frame(
    gene          = g,
    in_reference  = TRUE,
    total_counts  = sum(v),
    bins_detected = sum(v > 0),
    pct_bins      = round(100 * sum(v > 0) / length(v), 3),
    per_million   = round(1e6 * sum(v) / total_umi, 3)
  )
}

detection <- do.call(rbind, lapply(c(GENE_OF_INTEREST, REFERENCE_GENES), detection_row))
detection$sample    <- SAMPLE
detection$bin       <- BIN
detection$total_umi <- total_umi

write.csv(detection, file.path(OUT_DIR, "detection_check.csv"), row.names = FALSE)
message("\nDetection check (raw counts, ", format(total_umi, big.mark = ","), " UMIs total):")
print(detection[, c("gene", "total_counts", "bins_detected", "pct_bins", "per_million")])

# Plot the gene of interest regardless, so "not detected" is visible rather
# than asserted.
# Plotted off the Spatial assay on purpose. If the gene is barely there, SCT may
# have dropped it entirely, and where it hasn't, the corrected values would
# flatter a handful of counts into something that reads as signal.
if (GENE_OF_INTEREST %in% in_spatial) {
  DefaultAssay(obj) <- "Spatial"
  save_panel(spatial_plot(obj, GENE_OF_INTEREST), paste0("06_", GENE_OF_INTEREST))
  save_panel(umap_plot(obj, GENE_OF_INTEREST),    paste0("06_", GENE_OF_INTEREST, "_umap"))
  DefaultAssay(obj) <- "SCT"
} else {
  message(GENE_OF_INTEREST, " not in the filtered matrix - no panel drawn")
}

saveRDS(obj, file.path(OUT_DIR, paste0(SAMPLE, "_", BIN, "_marker_panels.rds")))
message("\nWritten to ", OUT_DIR)
