# Read the spaceranger cell-segmentation output into a Seurat object.
#
# WHY THIS EXISTS. Everything in spatial/ works on a square grid - 2um bins
# counted into tiles, labelled at 8um, pooled per tile for the allelic ratio. A
# tile is not a cell, and the seminar feedback was that a fixed grid both splits
# single cardiomyocytes across tiles and mixes a capillary into the myocyte
# beside it. spaceranger 4.1 segments the CytAssist image and assigns every 2um
# bin to a cell, so the same library can be read as cells instead. This script
# is step one of that alternative: load the segmented matrix, attach each cell's
# geometry, and describe what the segmentation actually produced. It does not
# cluster, annotate or touch allelic counts.
#
# NOTHING IS FILTERED OUT HERE. QC is computed and written as a `qc_pass` flag
# on every cell, so the thresholds can be judged from the figures before a later
# script commits to them. Subsetting is that script's decision, not this one's.
#
# WHAT IT READS, under <SPATIAL_DIR>/<SAMPLE>/outs/segmented_outputs/:
#   <MATRIX>_feature_cell_matrix.h5   counts, one column per called cell -
#                                     filtered_feature_cell_matrix.h5 by
#                                     default, MATRIX=raw for the other one.
#                                     Barcodes are "cellid_%09d-1", NOT spatial
#                                     barcodes, and the integer is the cell_id
#                                     used by the geojson below.
#   cell_segmentations.geojson        one Polygon per cell in full-resolution
#                                     image pixels, with cell_centroid and
#                                     nucleus_centroid in `properties`.
#   nucleus_segmentations.geojson     the nucleus polygon for the same cell_id.
#                                     Optional; skipped if absent.
#   spatial/scalefactors_json.json    microns_per_pixel, for areas and microns.
#   analysis/clustering/*/clusters.csv  spaceranger's own cell clusters, joined
#                                     in when present so its view is on the same
#                                     table as ours.
#
# THE COORDINATES ARE NOT IN THE H5. The segmented h5 is a plain 10x matrix -
# matrix/{barcodes,data,indices,indptr,shape,features} and nothing else. Every
# position, area and outline lives in the geojson, which is why this script
# parses a 50 MB json rather than reading a tissue_positions file like the
# binned scripts do.
#
# PIXELS TO MICRONS. `microns_per_pixel` from segmented_outputs/spatial/
# scalefactors_json.json (0.172619 for 9w), NOT the "Microscope Image Pixel
# Size" in metrics_summary.csv (0.173795). They differ by 0.7% and only the
# first is right: it puts a single 2um bin at exactly 2.00um across, and the
# median cell area it gives for 9w is 88.1um^2 against the 88.0um^2 that
# metrics_summary.csv reports. The other scale factor misses by a micron^2.
#
# CELL COUNTS DO NOT MATCH ACROSS THE FILES, and that is expected. For 9w:
#   raw h5        95610 cells      every segmented object
#   geojson       94608 polygons   those with an outline
#   filtered h5   94581 cells      those spaceranger calls cells   <- the default
# The filtered set is a strict subset of both of the others, so every cell it
# holds has a polygon and the join below drops nothing. Under MATRIX=raw, 1002
# barcodes have no polygon and therefore no position; those are dropped, with a
# count printed, because a cell with no coordinates is no use to a spatial
# analysis. The 27 cells that have a polygon but are not in the filtered matrix
# are segmented objects that failed spaceranger's own UMI filter.
#
# READ THE SEGMENTATION BEFORE TRUSTING IT. For 9w the median cell is 88um^2 and
# the median NUCLEUS is 8um^2 - two 2um bins. A real nuclear cross-section is
# 30-50um^2, and an adult cardiomyocyte cross-section is several hundred um^2.
# So these "cells" are small nuclear seeds grown outwards, not traced myocytes,
# and a myocyte is likely split across several of them. fig3 (area) and fig5
# (counts against area) are the panels that say whether that is tolerable for
# the question being asked. This is exactly the comparison the tiling feedback
# asked for, so the numbers are written out rather than assumed.
#
# WHAT IT WRITES, and it writes nowhere else, in <OUT_DIR>
# (default <SPATIAL_DIR>/spatial_V2/<SAMPLE>/ - a fresh directory, so nothing
# from the tile analysis can be overwritten; a re-run does overwrite its own
# previous outputs, which all carry fixed names):
#   <SAMPLE>_segmented_seurat.rds  Seurat object, assay "Segmented", one FOV
#                                  ("cells") of centroids so SpatialFeaturePlot
#                                  and the Seurat v5 niche functions work
#   cells.tsv.gz                   one row per cell: ids, centroid in pixels and
#                                  microns, cell and nucleus area, nCount,
#                                  nFeature, pct_mt, spaceranger cluster, qc_pass
#   qc_summary.tsv                 the counts and medians quoted above, per run
#   provenance.tsv                 inputs, sizes, thresholds, session info
#   fig1  tissue map, cells coloured by nCount                     <- the picture
#   fig2  nCount / nFeature / pct_mt distributions with thresholds
#   fig3  cell and nucleus area distributions
#   fig4  barcode-rank (knee) curve, raw and filtered marked
#   fig5  nCount against cell area - the segmentation sanity check
#   fig6  tissue map of qc_pass, so discards can be seen to be where they should
#   fig7  Xist and a marker panel per cell, if the genes are present
#
# Run on the cluster, one sample per invocation:
#   sbatch slurm/spatial_v2_load_cells.slurm 9w 78w
# or by hand from adult_aged_spatial/:
#   SAMPLE=9w Rscript ~/Postdoc/spatial_V2/01_load_segmented_cells.R

# Must be line 1 of anything that touches Seurat/Matrix on this cluster (see
# spatial/banksy_patch.R for why).
.libPaths(c("~/R/matrix-dev", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})
for (p in c("data.table", "jsonlite", "patchwork")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("R package '", p, "' is required")
}
library(data.table)
library(patchwork)

##### ---------------------- CONFIG ---------------------- #####
SPATIAL_DIR <- Sys.getenv("SPATIAL_DIR",
  "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial")
BASE    <- SPATIAL_DIR
SAMPLE  <- Sys.getenv("SAMPLE", "9w")            # 9w = adult, 78w = aged
MATRIX  <- Sys.getenv("MATRIX", "filtered")      # filtered | raw
OUT_DIR <- Sys.getenv("OUT_DIR", "")
# A directory of its own, named after this code directory, so the V2 outputs
# never land beside the tile analysis in annotation/ or ase*/. Nothing else in
# adult_aged_spatial/ writes here.
if (!nzchar(OUT_DIR)) OUT_DIR <- file.path(BASE, "spatial_V2", SAMPLE)
stopifnot(MATRIX %in% c("raw", "filtered"))

env_flag <- function(name, default) {
  v <- toupper(trimws(Sys.getenv(name, "")))
  if (!nzchar(v)) return(default)
  if (v %in% c("1", "TRUE", "T", "YES")) return(TRUE)
  if (v %in% c("0", "FALSE", "F", "NO"))  return(FALSE)
  stop("Cannot read ", name, " = '", Sys.getenv(name), "' as a logical")
}
env_num <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (!nzchar(v)) return(default)
  x <- suppressWarnings(as.numeric(v))
  if (!is.finite(x)) stop(name, " = '", v, "' is not a number")
  x
}

# QC thresholds. These FLAG, they do not filter - see the header. The defaults
# are deliberately loose: for 9w the median cell carries 199 UMIs over 122
# genes, so a 100-UMI floor would already discard a third of the section, and
# cardiomyocytes are legitimately 30-40% mitochondrial.
MIN_COUNT <- env_num("MIN_COUNT", 50)
MIN_FEAT  <- env_num("MIN_FEAT", 25)
MAX_MT    <- env_num("MAX_MT", 50)
# Area bounds in um^2. The floor drops one- and two-bin specks; the ceiling
# flags likely merges of several cells. A 2um bin is 4um^2, so 20um^2 is five
# bins. Set MAX_AREA to 0 to disable the ceiling.
MIN_AREA  <- env_num("MIN_AREA", 20)
MAX_AREA  <- env_num("MAX_AREA", 1000)

# spaceranger's own clustering of the same cells, joined onto cells.tsv.gz for
# comparison. Empty string skips it.
SR_CLUSTERING <- Sys.getenv("SR_CLUSTERING", "gene_expression_graphclust")

MT_PATTERN <- Sys.getenv("MT_PATTERN", "^mt-")   # GRCm39 gene symbols
GENE_PANEL <- strsplit(Sys.getenv("GENE_PANEL", "Xist,Ttn,Myl2,Nppa,Pecam1,Dcn"),
                       "\\s*[,;]\\s*")[[1]]

SAVE_RDS   <- env_flag("SAVE_RDS", TRUE)
POINT_SIZE <- env_num("POINT_SIZE", 0.15)
FLIP_X <- env_flag("FLIP_X", FALSE)              # as in spatial/tile_ratio_map.R
FLIP_Y <- env_flag("FLIP_Y", FALSE)
RASTERISE <- TRUE

SEG_DIR <- file.path(BASE, SAMPLE, "outs", "segmented_outputs")
H5      <- file.path(SEG_DIR, paste0(MATRIX, "_feature_cell_matrix.h5"))
MTX_DIR <- file.path(SEG_DIR, paste0(MATRIX, "_feature_cell_matrix"))
CELL_GJ <- file.path(SEG_DIR, "cell_segmentations.geojson")
NUC_GJ  <- file.path(SEG_DIR, "nucleus_segmentations.geojson")
SCALEF  <- file.path(SEG_DIR, "spatial", "scalefactors_json.json")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
message("sample     : ", SAMPLE)
message("segmented  : ", SEG_DIR)
message("matrix     : ", MATRIX)
message("out        : ", OUT_DIR)

##### ---------------------- 1. COUNTS ---------------------- #####
# Read10X_h5 needs hdf5r; fall back to the mtx directory beside the h5, which
# holds the same matrix, rather than failing the job over a missing package.
if (file.exists(H5) && requireNamespace("hdf5r", quietly = TRUE)) {
  message("reading ", basename(H5))
  counts <- Read10X_h5(H5)
} else if (dir.exists(MTX_DIR)) {
  if (file.exists(H5))
    message("note: hdf5r not installed, reading ", basename(MTX_DIR), " instead")
  message("reading ", MTX_DIR)
  counts <- Read10X(MTX_DIR)
} else {
  stop("No segmented matrix found. Looked for:\n  ", H5, "\n  ", MTX_DIR)
}
# A multi-modal h5 would come back as a list; these libraries are Gene
# Expression only, so say so rather than silently taking the first element.
if (is.list(counts)) {
  if (!"Gene Expression" %in% names(counts))
    stop("Matrix has no 'Gene Expression' modality: ", paste(names(counts), collapse = ", "))
  counts <- counts[["Gene Expression"]]
}
message("  ", nrow(counts), " genes x ", ncol(counts), " cells")

# "cellid_000000001-1" -> 1. The integer is the join key to both geojsons.
cell_id <- as.integer(sub("-1$", "", sub("^cellid_0*", "", colnames(counts))))
if (anyNA(cell_id))
  stop("Could not parse a cell_id out of every barcode - is this really a ",
       "segmented matrix? First barcode: ", colnames(counts)[1])

##### ---------------------- 2. GEOMETRY ---------------------- #####
if (!file.exists(SCALEF)) stop("Missing ", SCALEF)
MPP <- jsonlite::fromJSON(SCALEF)$microns_per_pixel
if (!is.numeric(MPP) || !is.finite(MPP)) stop("No usable microns_per_pixel in ", SCALEF)
message("microns_per_pixel: ", format(MPP, digits = 9))

# Shoelace formula. The polygons are small, closed and non-self-intersecting, so
# this is exact and avoids a dependency on sf, which is not in seurat_env.
poly_area_px <- function(m) {
  n <- nrow(m)
  if (n < 3) return(0)
  j <- c(2:n, 1)
  abs(sum(m[, 1] * m[j, 2] - m[j, 1] * m[, 2])) / 2
}
# jsonlite hands a Polygon's coordinates back as a list of rings, or as a
# 1 x npoints x 2 array when every ring in the file has the same length. Take
# the outer ring either way; these are simple polygons with no holes.
outer_ring <- function(co) {
  if (is.list(co)) as.matrix(co[[1]]) else matrix(co[1, , ], ncol = 2)
}

read_segmentation <- function(path, what) {
  message("parsing ", basename(path), " (", what, ") - this reads ~50 MB of json")
  gj <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  pr <- gj$features$properties
  ce <- gj$features$geometry$coordinates
  d  <- data.table(cell_id = as.integer(pr$cell_id),
                   area_px = vapply(ce, function(co) poly_area_px(outer_ring(co)), numeric(1)))
  # Both files carry a nucleus_centroid; only the cell file carries a
  # cell_centroid, so take each centroid from the file that defines it.
  cen_col <- if (what == "cell") "cell_centroid" else "nucleus_centroid"
  cen <- do.call(rbind, pr[[cen_col]])
  d[, `:=`(cx_px = cen[, 1], cy_px = cen[, 2])]
  d
}

if (!file.exists(CELL_GJ)) stop("Missing ", CELL_GJ)
cells <- read_segmentation(CELL_GJ, "cell")
setnames(cells, c("area_px", "cx_px", "cy_px"), c("cell_area_px", "x_px", "y_px"))
message("  ", nrow(cells), " cell polygons")

nuc <- NULL
if (file.exists(NUC_GJ)) {
  nuc <- read_segmentation(NUC_GJ, "nucleus")
  setnames(nuc, c("area_px", "cx_px", "cy_px"),
           c("nucleus_area_px", "nucleus_x_px", "nucleus_y_px"))
  message("  ", nrow(nuc), " nucleus polygons")
} else {
  message("note: no nucleus_segmentations.geojson - nucleus columns will be NA")
}

##### ---------------------- 3. JOIN ---------------------- #####
meta <- data.table(barcode = colnames(counts), cell_id = cell_id)
meta <- cells[meta, on = "cell_id"]            # keep every matrix column for now
if (!is.null(nuc)) meta <- nuc[meta, on = "cell_id"]
setkey(meta, NULL)
setcolorder(meta, c("barcode", "cell_id"))
# Restore matrix order: the join above reorders, and the object must line up
# with the columns of `counts`.
meta <- meta[match(colnames(counts), barcode)]
stopifnot(identical(meta$barcode, colnames(counts)))

no_poly <- is.na(meta$x_px)
if (any(no_poly)) {
  message("dropping ", sum(no_poly), " of ", ncol(counts),
          " cells with no polygon in cell_segmentations.geojson (no position)")
  keep   <- !no_poly
  counts <- counts[, keep, drop = FALSE]
  meta   <- meta[keep]
}
if (ncol(counts) == 0) stop("No cell in the matrix has a segmentation polygon")

meta[, `:=`(
  x         = x_px * if (FLIP_X) -1 else 1,
  y         = y_px * if (FLIP_Y) -1 else 1,
  x_um      = x_px * MPP,
  y_um      = y_px * MPP,
  cell_area_um2 = cell_area_px * MPP^2
)]
if (!is.null(nuc)) meta[, nucleus_area_um2 := nucleus_area_px * MPP^2]

##### ---------------------- 4. OBJECT AND QC ---------------------- #####
obj <- CreateSeuratObject(counts = counts, assay = "Segmented", project = SAMPLE)
obj$sample  <- SAMPLE
obj$cell_id <- meta$cell_id
for (cl in c("x", "y", "x_px", "y_px", "x_um", "y_um", "cell_area_um2",
             "nucleus_area_um2", "nucleus_x_px", "nucleus_y_px")) {
  if (cl %in% names(meta)) obj[[cl]] <- meta[[cl]]
}
obj$pct_mt <- PercentageFeatureSet(obj, pattern = MT_PATTERN)
if (all(obj$pct_mt == 0))
  message("note: MT_PATTERN '", MT_PATTERN, "' matched no gene - pct_mt is all zero")

nc <- obj$nCount_Segmented
nf <- obj$nFeature_Segmented
obj$qc_pass <- nc >= MIN_COUNT & nf >= MIN_FEAT & obj$pct_mt <= MAX_MT &
  obj$cell_area_um2 >= MIN_AREA & (MAX_AREA <= 0 | obj$cell_area_um2 <= MAX_AREA)
message("qc_pass: ", sum(obj$qc_pass), " / ", ncol(obj),
        " (", round(100 * mean(obj$qc_pass), 1), "%)")

# spaceranger's own clusters, for the same cells.
if (nzchar(SR_CLUSTERING)) {
  cl_csv <- file.path(SEG_DIR, "analysis", "clustering", SR_CLUSTERING, "clusters.csv")
  if (file.exists(cl_csv)) {
    cl <- fread(cl_csv)
    obj$sr_cluster <- factor(cl$Cluster[match(colnames(obj), cl$Barcode)])
    message("joined ", SR_CLUSTERING, ": ", sum(!is.na(obj$sr_cluster)), " cells labelled")
  } else {
    message("note: no ", SR_CLUSTERING, "/clusters.csv - sr_cluster skipped")
  }
}

# A centroids FOV, so the object is spatial to Seurat itself and not only to the
# ggplot maps below: SpatialFeaturePlot, ImageDimPlot and the v5 niche functions
# all read it. Guarded - the object is still worth writing if this fails.
fov_ok <- tryCatch({
  cf <- data.frame(x = obj$x_px, y = obj$y_px, cell = colnames(obj),
                   stringsAsFactors = FALSE)
  obj[["cells"]] <- CreateFOV(cf, type = "centroids", assay = "Segmented")
  TRUE
}, error = function(e) { message("note: could not build the FOV - ", conditionMessage(e)); FALSE })

##### ---------------------- 5. TABLES ---------------------- #####
out_cols <- intersect(c("barcode", "cell_id", "x_px", "y_px", "x_um", "y_um",
                        "cell_area_um2", "nucleus_area_um2",
                        "nucleus_x_px", "nucleus_y_px"), names(meta))
cells_out <- meta[, ..out_cols]
cells_out[, `:=`(nCount = as.numeric(nc), nFeature = as.numeric(nf),
                 pct_mt = as.numeric(obj$pct_mt), qc_pass = obj$qc_pass)]
if ("sr_cluster" %in% names(obj[[]])) cells_out[, sr_cluster := obj$sr_cluster]
fwrite(cells_out, file.path(OUT_DIR, "cells.tsv.gz"), sep = "\t")

med <- function(x) if (all(is.na(x))) NA_real_ else round(median(x, na.rm = TRUE), 2)
qc <- data.table(
  sample              = SAMPLE,
  matrix              = MATRIX,
  cells_in_matrix     = length(cell_id),
  cells_with_polygon  = ncol(obj),
  cells_dropped_no_polygon = sum(no_poly),
  cells_qc_pass       = sum(obj$qc_pass),
  genes               = nrow(obj),
  median_nCount       = med(nc),
  median_nFeature     = med(nf),
  median_pct_mt       = med(obj$pct_mt),
  median_cell_area_um2    = med(obj$cell_area_um2),
  median_nucleus_area_um2 = if ("nucleus_area_um2" %in% names(obj[[]])) med(obj$nucleus_area_um2) else NA_real_,
  microns_per_pixel   = MPP,
  MIN_COUNT = MIN_COUNT, MIN_FEAT = MIN_FEAT, MAX_MT = MAX_MT,
  MIN_AREA = MIN_AREA, MAX_AREA = MAX_AREA)
fwrite(qc, file.path(OUT_DIR, "qc_summary.tsv"), sep = "\t")
print(t(qc))

##### ---------------------- 6. FIGURES ---------------------- #####
theme_map <- theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank())

ras <- function(g) {
  if (RASTERISE && requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(g, dpi = 300) else g
}
save_fig <- function(name, p, w = 8, h = 7) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}

df <- as.data.table(obj[[]])
df[, `:=`(x = obj$x, y = obj$y)]

# fig1 - the section, cells coloured by depth. Clipped at q99: a handful of
# merged cells otherwise take the whole colour scale.
hi <- quantile(df$nCount_Segmented, 0.99)
p1 <- ggplot(df, aes(x, y, colour = pmin(nCount_Segmented, hi))) +
  geom_point(size = POINT_SIZE, shape = 16) +
  scale_colour_viridis_c(name = "UMIs/cell", option = "magma") +
  scale_y_reverse() + coord_fixed() + theme_map +
  labs(title = paste0(SAMPLE, ": ", format(ncol(obj), big.mark = ","),
                      " segmented cells"),
       subtitle = "colour clipped at the 99th percentile")
save_fig("fig1_cells_nCount", ras(p1))

# fig2 - depth and mitochondrial content, with the flags drawn on.
vln <- function(col, thr, lab, logx = TRUE) {
  g <- ggplot(df, aes(x = 1, y = .data[[col]])) +
    geom_violin(fill = "grey85", colour = NA) +
    geom_hline(yintercept = thr, colour = "firebrick", linetype = 2) +
    labs(x = NULL, y = lab, title = paste0(lab, "  (line: ", thr, ")")) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
  if (logx) g <- g + scale_y_log10()
  g
}
p2 <- vln("nCount_Segmented", MIN_COUNT, "UMIs per cell") |
  vln("nFeature_Segmented", MIN_FEAT, "genes per cell") |
  vln("pct_mt", MAX_MT, "% mitochondrial", logx = FALSE)
save_fig("fig2_qc_distributions", p2, w = 9, h = 4)

# fig3 - geometry. The panel that says whether these are cells: for 9w the
# nucleus median is two 2um bins, well under a real nuclear cross-section.
area_hist <- function(col, lab, vlines) {
  ggplot(df[is.finite(get(col))], aes(.data[[col]])) +
    geom_histogram(bins = 80, fill = "grey40") +
    geom_vline(xintercept = vlines, colour = "firebrick", linetype = 2) +
    scale_x_log10() +
    labs(x = expression(um^2), y = "cells", title = lab) +
    theme_minimal(base_size = 9)
}
p3 <- area_hist("cell_area_um2", "cell area", c(MIN_AREA, if (MAX_AREA > 0) MAX_AREA))
if ("nucleus_area_um2" %in% names(df))
  p3 <- p3 | area_hist("nucleus_area_um2", "nucleus area", numeric(0))
save_fig("fig3_area", p3, w = 8, h = 4)

# fig4 - barcode rank. Where the matrix stops being cells; MATRIX=filtered cuts
# somewhere on this curve, so it is worth seeing before choosing one.
kn <- data.table(rank = seq_len(ncol(obj)), n = sort(nc, decreasing = TRUE))
p4 <- ggplot(kn[n > 0], aes(rank, n)) +
  geom_line() + scale_x_log10() + scale_y_log10() +
  geom_hline(yintercept = MIN_COUNT, colour = "firebrick", linetype = 2) +
  labs(x = "cell rank", y = "UMIs", title = paste0(SAMPLE, " (", MATRIX, "): barcode rank"),
       subtitle = paste0("line: MIN_COUNT = ", MIN_COUNT)) +
  theme_minimal(base_size = 9)
save_fig("fig4_knee", p4, w = 5, h = 4)

# fig5 - counts against area. A segmentation that is capturing cells gives a
# rising cloud; a flat one means area is being assigned independently of where
# the transcripts are.
p5 <- ggplot(df[cell_area_um2 > 0 & nCount_Segmented > 0],
             aes(cell_area_um2, nCount_Segmented)) +
  geom_point(size = 0.1, alpha = 0.15, shape = 16) +
  scale_x_log10() + scale_y_log10() +
  labs(x = expression(cell~area~(um^2)), y = "UMIs per cell",
       title = "depth against size",
       subtitle = sprintf("Spearman rho = %.2f",
                          cor(df$cell_area_um2, df$nCount_Segmented,
                              method = "spearman", use = "complete.obs"))) +
  theme_minimal(base_size = 9)
# mgcv is a recommended package but is not in every conda R, and a missing
# smoother is not worth failing the figure over.
if (requireNamespace("mgcv", quietly = TRUE)) {
  p5 <- p5 + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                         colour = "firebrick", se = FALSE)
}
save_fig("fig5_counts_vs_area", ras(p5), w = 5, h = 4.5)

# fig6 - where the discards are. Scattered is fine; a whole region failing means
# the thresholds are wrong for that tissue, not that the tissue is bad.
p6 <- ggplot(df, aes(x, y, colour = qc_pass)) +
  geom_point(size = POINT_SIZE, shape = 16) +
  scale_colour_manual(values = c(`TRUE` = "grey70", `FALSE` = "firebrick"), name = "qc_pass") +
  scale_y_reverse() + coord_fixed() + theme_map +
  labs(title = paste0("qc_pass: ", sum(df$qc_pass), " / ", nrow(df)))
save_fig("fig6_qc_map", ras(p6))

# fig7 - a few genes per cell, the visual check that the object is the section.
panel <- intersect(GENE_PANEL, rownames(obj))
if (length(panel)) {
  ex <- LayerData(obj, assay = "Segmented", layer = "counts")[panel, , drop = FALSE]
  gl <- lapply(panel, function(g) {
    d <- copy(df)[, e := as.numeric(ex[g, ])]
    d <- d[order(e)]                       # expressing cells drawn on top
    ggplot(d, aes(x, y, colour = e)) +
      geom_point(size = POINT_SIZE, shape = 16) +
      scale_colour_viridis_c(name = "UMIs", option = "viridis",
                             limits = c(0, max(1, quantile(d$e, 0.999)))) +
      scale_y_reverse() + coord_fixed() + theme_map + labs(title = g)
  })
  save_fig("fig7_gene_panel", wrap_plots(lapply(gl, ras), ncol = 3),
           w = 12, h = 4 * ceiling(length(panel) / 3))
} else {
  message("note: none of GENE_PANEL is in the matrix - fig7 skipped")
}

##### ---------------------- 7. OBJECT AND PROVENANCE ---------------------- #####
if (SAVE_RDS) {
  rds <- file.path(OUT_DIR, paste0(SAMPLE, "_segmented_seurat.rds"))
  message("writing ", rds)
  saveRDS(obj, rds)
}

prov <- data.table(
  field = c("script", "run_at", "sample", "matrix", "seg_dir", "h5", "cell_geojson",
          "nucleus_geojson", "microns_per_pixel", "cells_in_matrix", "cells_kept",
          "qc_pass", "fov", "sr_clustering", "R", "Seurat"),
  value = c("spatial_V2/01_load_segmented_cells.R",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), SAMPLE, MATRIX, SEG_DIR,
            if (file.exists(H5)) H5 else MTX_DIR, CELL_GJ,
            if (is.null(nuc)) "absent" else NUC_GJ,
            format(MPP, digits = 9), length(cell_id), ncol(obj), sum(obj$qc_pass),
            if (fov_ok) "cells" else "none", SR_CLUSTERING,
            paste(R.version$major, R.version$minor, sep = "."),
            as.character(packageVersion("Seurat"))))
fwrite(prov, file.path(OUT_DIR, "provenance.tsv"), sep = "\t")

message("done: ", OUT_DIR)
