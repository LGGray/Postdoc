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
names(umap)[names(umap) != "barcode"][1:2] <- c("UMAP_1", "UMAP_2")

# ---------------------------------------------------------------- assemble

obj <- CreateSeuratObject(counts, project = SAMPLE, assay = "Spatial")
obj <- NormalizeData(obj, verbose = FALSE)

bc <- colnames(obj)
idx <- function(df) match(bc, df$barcode)

# pxl_row is the vertical axis and pxl_col the horizontal one. Image coordinates
# run downwards from the top-left, so y is reversed at plotting time.
meta <- data.frame(
  x       = pos$pxl_col_in_fullres[idx(pos)],
  y       = pos$pxl_row_in_fullres[idx(pos)],
  UMAP_1  = umap$UMAP.1[idx(umap)],
  UMAP_2  = umap$UMAP.2[idx(umap)],
  cluster = factor(clusters$Cluster[idx(clusters)]),
  row.names = bc
)
if (anyNA(meta$x) || anyNA(meta$UMAP_1)) {
  warning("Some barcodes had no coordinate or UMAP match - check bin size consistency")
}
obj <- AddMetaData(obj, meta)

# ---------------------------------------------------------------- scores

present <- function(g) intersect(g, rownames(obj))

for (nm in names(MARKERS)) {
  found <- present(MARKERS[[nm]])
  missing <- setdiff(MARKERS[[nm]], found)
  if (length(missing)) message("  ", nm, ": not in reference - ", paste(missing, collapse = ", "))
  if (!length(found)) { message("  ", nm, ": skipped, no markers found"); next }
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
if (requireNamespace("patchwork", quietly = TRUE)) {
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
raw <- GetAssayData(obj, layer = "counts")
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
if (GENE_OF_INTEREST %in% rownames(obj)) {
  save_panel(spatial_plot(obj, GENE_OF_INTEREST), paste0("06_", GENE_OF_INTEREST))
  save_panel(umap_plot(obj, GENE_OF_INTEREST),    paste0("06_", GENE_OF_INTEREST, "_umap"))
}

saveRDS(obj, file.path(OUT_DIR, paste0(SAMPLE, "_", BIN, "_marker_panels.rds")))
message("\nWritten to ", OUT_DIR)
