# Representative marker panels for the Visium HD cardiac sections.
#
# Deliverable: a handful of slide-ready figures for a 6-minute talk, not a full
# analysis. Everything is descriptive - one section, no statistics, no claim
# beyond "the tissue and the assay behave the way they should".
#
#   fig1   chamber orientation (Myl2/Myl7/Nppa/Ttn) - where the LV is
#   fig1b  UMIs per bin - the QC panel, in case anyone asks
#   fig2   cluster annotation: tissue map beside the matching UMAP   <- hero
#   fig2b  raw graph-based clusters, unannotated
#   fig3   endothelial score: tissue beside UMAP                     <- the ask
#   fig4   fibroblast score:  tissue beside UMAP                     <- the ask
#   fig5   both signatures in one 2x2, for the single-slide version
#   fig6   endothelial marker genes, one panel each
#   fig7   fibroblast marker genes (+ Gsn), one panel each
#   fig8   endothelial subtype markers (only if detectable)
#   fig9   Mex3a detection check: positive bins on tissue + a rate comparison
#
# Run on the cluster from adult_aged_spatial/. The sample is an argument to the
# wrapper (default 78w), so both sections use the same script:
#   sbatch ~/Postdoc/slurm/spatial_marker_panels.slurm 78w
#   sbatch ~/Postdoc/slurm/spatial_marker_panels.slurm 9w 78w
#
# Uses spaceranger's own UMAP and graph clustering rather than recomputing, so
# the panels match what the web summary shows and the script stays cheap.

# Must be line 1 of anything that touches Seurat/Matrix on this cluster - the
# stock Matrix in the conda env is too old for the sparse-matrix classes Seurat
# hands around (see banksy_patch.R).
.libPaths(c("~/R/matrix-dev", .libPaths()))

library(Seurat)
library(ggplot2)
library(Matrix)

##### ---------------------- CONFIG ---------------------- #####
# Sample, bin and output directory come from the environment so one script
# serves both sections without an edit; the defaults reproduce the original 9w
# run. SPATIAL_DIR rather than BASE on purpose - the slurm wrapper already uses
# BASE for the project root, one level up, and picking that up silently would
# send every path to the wrong place.
#
#   SAMPLE=78w Rscript spatial_marker_panels.R
SPATIAL_DIR <- Sys.getenv("SPATIAL_DIR",
  "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial")
BASE      <- SPATIAL_DIR
SAMPLE    <- Sys.getenv("SAMPLE", "9w")   # 9w = adult baseline, 78w = aged
BIN       <- Sys.getenv("BIN", "square_016um")
OUT_DIR   <- Sys.getenv("OUT_DIR", "")
if (!nzchar(OUT_DIR)) OUT_DIR <- file.path(BASE, "figures", paste0(SAMPLE, "_", BIN))

# Environment variables arrive as strings; "FALSE" is TRUE to as.logical()'s
# careless cousins and silently NA to some, so parse the flags explicitly.
env_flag <- function(name, default) {
  v <- toupper(trimws(Sys.getenv(name, "")))
  if (!nzchar(v)) return(default)
  if (v %in% c("1", "TRUE", "T", "YES")) return(TRUE)
  if (v %in% c("0", "FALSE", "F", "NO"))  return(FALSE)
  stop("Cannot read ", name, " = '", Sys.getenv(name), "' as a logical")
}

# Point size for the spatial plots. Tune if bins look too sparse or too merged:
# smaller bins (008um) need a smaller value, 032um a larger one.
POINT_SIZE <- as.numeric(Sys.getenv("POINT_SIZE", "0.35"))
if (!is.finite(POINT_SIZE) || POINT_SIZE <= 0) stop("POINT_SIZE must be a positive number")

# Section orientation. Visium images are not filed apex-down; flip until the LV
# free wall sits where an anatomist expects it. y is always reversed on top of
# this, because image coordinates run downwards from the top-left.
# Leave both FALSE to keep these panels in the same frame as the allelic-ratio
# tile maps (tile_gene_ar_maps.R), which never flip - a flip here would make the
# two figures look comparable while being mirrored relative to each other.
FLIP_X <- env_flag("FLIP_X", FALSE)
FLIP_Y <- env_flag("FLIP_Y", FALSE)

# Light panels (grey zeros, dark signal) print and project cleanly on the usual
# white slide template. DARK_MODE swaps to black background + magma, which is
# more striking on a dark template but unusable on paper.
DARK_MODE <- FALSE

# Vector output at 50-100k bins makes a PDF no laptop wants to open in a talk.
# Rasterise the point layer and keep the text/axes as vectors.
RASTERISE <- TRUE

# Normalisation for the score and marker panels: "LogNormalize" or "SCT".
#
# LogNormalize by default, matching Seurat's own Visium HD vignette and
# visium_hd_test.R. SCTransform is what the spot-level Visium vignette
# recommends, on the argument that per-spot depth tracks cell density and
# anatomy - but its real payoff is variance stabilisation for HVG selection and
# clustering, and this script does neither: the UMAP and the clusters come from
# spaceranger. What is left is the cost. SCT drops genes below min_cells = 5,
# which is the wrong filter when the deliverable is "is Mex3a detectable" and
# "plot the EC subtype markers if they are there", and it puts the panels on a
# different scale from the detection check.
NORMALISATION <- "LogNormalize"
stopifnot(NORMALISATION %in% c("LogNormalize", "SCT"))

# Marker sets. Endothelial and fibroblast are the collaborator's; the rest are
# standard mouse heart markers for context, deliberately coarse - 16um bins hold
# several cells, so anything finer than these classes is over-reading the data.
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

# Low end of the colour ramp for module-score panels, as a quantile.
#
# Scores are unimodal and centred near zero with a long right tail. Anchoring
# the ramp at the minimum puts the bulk of the tissue in its middle - which is
# what made the first endothelial panel read as uniform orange noise - and
# wastes the grey end entirely. At 0.5 the lower half of the bins go flat grey
# and the whole ramp is spent on the signal. Raise towards 0.7 for a sparser,
# punchier panel; 0 restores the old min-to-q99 behaviour.
SCORE_QLO <- 0.5

# Neighbourhood radius, in bins, for the smoothed score panels (0 disables).
# Capillary endothelium sits in nearly every bin at a count or two, so the raw
# score is legitimately salt-and-pepper - true, but unreadable from the back of
# a room. Averaging over a (2r+1)^2 box turns it into a density map. At 16um,
# r = 2 is an 80um kernel. The unsmoothed panel is always written too; the
# smoothed one is a visualisation, not a new measurement, and is labelled so.
SMOOTH_RADIUS <- 2

# Cluster annotation thresholds (see the annotation block). A cluster is named
# after its highest mean marker-set z-score, but only when that score is
# positive and clears the runner-up by MIN_MARGIN. Inspect
# cluster_annotation.csv and retune rather than guessing twice.
MIN_Z      <- 0
MIN_MARGIN <- 0.10

# Endothelial subtypes - only plotted if they look detectable.
EC_SUBTYPES <- c("Bmx", "Mecom", "Gja5", "Rgcc", "Ackr1", "Fabp4")

# Chamber orientation. Myl2 = ventricular (so: LV/RV free wall and septum),
# Myl7 = atrial, Nppa = atrial + trabecular/stressed myocardium. These three are
# what lets a viewer point at the LV in the other panels.
ORIENTATION_GENES <- c("Myl2", "Myl7", "Nppa", "Ttn")

# Fibroblast activation marker, shown next to the resting fibroblast markers.
FB_EXTRA <- "Gsn"

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

# Detection of the gene of interest, measured before any filtering. Whether a
# gene is seen at all should not depend on QC thresholds chosen for plotting.
prefilter_detection <- local({
  g <- c(GENE_OF_INTEREST, REFERENCE_GENES)
  tot <- sum(counts)
  do.call(rbind, lapply(g, function(gene) {
    if (!gene %in% rownames(counts)) {
      return(data.frame(gene = gene, in_reference = FALSE, total_counts = NA_real_,
                        bins_detected = NA_integer_, pct_bins = NA_real_,
                        per_million = NA_real_))
    }
    v <- counts[gene, ]
    data.frame(gene = gene, in_reference = TRUE, total_counts = sum(v),
               bins_detected = sum(v > 0),
               pct_bins    = round(100 * sum(v > 0) / length(v), 4),
               per_million = round(1e6 * sum(v) / tot, 4))
  }))
})
prefilter_detection$stage     <- "pre-QC (all bins in filtered matrix)"
prefilter_detection$total_umi <- sum(counts)

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

# PercentageFeatureSet errors outright when the pattern matches nothing, which
# takes the script down over a metric that is only ever printed.
pct_or_zero <- function(o, pattern) {
  if (!any(grepl(pattern, rownames(o)))) {
    message("No features match ", pattern, " - reporting 0")
    return(setNames(rep(0, ncol(o)), colnames(o)))
  }
  # Seurat v4 returns a one-column data.frame here, v5 a named numeric vector,
  # so [, 1] is not portable. v5 also assembles it layer by layer, which need
  # not follow colnames(o) - reorder by cell name whenever names are available.
  p <- PercentageFeatureSet(o, pattern = pattern)
  v <- if (is.data.frame(p) || is.matrix(p)) setNames(p[, 1], rownames(p)) else p
  if (!is.null(names(v)) && all(colnames(o) %in% names(v))) v <- v[colnames(o)]
  v
}
obj[["percent.mt"]]   <- pct_or_zero(obj, "^mt-")
obj[["percent.hb"]]   <- pct_or_zero(obj, "^Hb[ab]-")
obj[["percent.ribo"]] <- pct_or_zero(obj, "^Rp[sl]")

message("nCount  quantiles: ", paste(round(quantile(obj$nCount_Spatial,  c(.01,.05,.25,.5,.95))), collapse=" "))
message("nFeature quantiles: ", paste(round(quantile(obj$nFeature_Spatial, c(.01,.05,.25,.5,.95))), collapse=" "))
message("percent.mt median: ", round(median(obj$percent.mt), 1))

MIN_COUNT <- 100   # 16um; drop to ~15 for 008um
MIN_FEAT  <- 50
MAX_MT    <- 50   # loose on purpose - cardiomyocytes are legitimately 30-40%

keep <- obj$nCount_Spatial >= MIN_COUNT & obj$nFeature_Spatial >= MIN_FEAT & obj$percent.mt <= MAX_MT
message("Keeping ", sum(keep), "/", ncol(obj), " bins (", round(100*mean(keep),1), "%)")
obj <- obj[, keep]

# The Spatial assay is log-normalised either way. The detection check and the
# gene-of-interest panel read its counts layer, never a normalised one: the
# question there is whether the gene is seen at all, and any normalisation is
# the wrong lens for that.
obj <- NormalizeData(obj, assay = "Spatial", verbose = FALSE)

if (NORMALISATION == "SCT") {
  # vst.flavor "v2" needs glmGamPoi; fall back to v1 rather than erroring on a
  # fresh cluster env. v1's parameter regularisation is the part Lause et al.
  # criticised, so v2 is worth having.
  vst_flavor <- if (requireNamespace("glmGamPoi", quietly = TRUE)) "v2" else NULL
  if (is.null(vst_flavor)) {
    message("glmGamPoi not found - using sctransform v1 (slower). ",
            "BiocManager::install('glmGamPoi') to speed this up.")
  }
  message("SCTransform ...")
  obj <- SCTransform(obj, assay = "Spatial", new.assay.name = "SCT",
                     vst.flavor = vst_flavor, verbose = FALSE)
}

# Everything downstream reads PANEL_ASSAY, so the switch is one variable.
PANEL_ASSAY <- if (NORMALISATION == "SCT") "SCT" else "Spatial"
DefaultAssay(obj) <- PANEL_ASSAY
message("Panels drawn from the ", PANEL_ASSAY, " assay (", NORMALISATION, ")")


bc <- colnames(obj)
idx <- function(df) match(bc, df$barcode)

# pxl_row is the vertical axis and pxl_col the horizontal one. Image coordinates
# run downwards from the top-left, so y is reversed at plotting time.
# array_row/array_col are the bin's integer position on the capture grid, which
# is what the smoothing box filter indexes into. Carried alongside the pixel
# coordinates rather than derived from them - the pixel grid is rotated and
# scaled relative to the array, and rounding it back would be guesswork.
meta <- data.frame(
  x       = pos$pxl_col_in_fullres[idx(pos)],
  y       = pos$pxl_row_in_fullres[idx(pos)],
  arow    = pos$array_row[idx(pos)],
  acol    = pos$array_col[idx(pos)],
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
if (SMOOTH_RADIUS > 0 && (all(is.na(meta$arow)) || all(is.na(meta$acol)))) {
  message("tissue_positions has no array_row/array_col - disabling smoothed panels")
  SMOOTH_RADIUS <- 0
}
if (FLIP_X) meta$x <- -meta$x
if (FLIP_Y) meta$y <- -meta$y
obj <- AddMetaData(obj, meta)

# ---------------------------------------------------------------- scores

# Under SCT a gene can go missing for two quite different reasons: absent from
# the reference, or present but below SCT's min_cells = 5 and therefore not
# modelled. Those are different statements to make about a gene, so keep them
# apart. Under LogNormalize the second set is empty by construction, which is
# most of the reason for preferring it here.
in_spatial <- rownames(obj[["Spatial"]])
in_panel   <- rownames(obj[[PANEL_ASSAY]])

present <- function(g) intersect(g, in_panel)

report_missing <- function(g, nm) {
  absent <- setdiff(g, in_spatial)
  sparse <- setdiff(intersect(g, in_spatial), in_panel)
  if (length(absent)) message("  ", nm, ": not in reference - ", paste(absent, collapse = ", "))
  if (length(sparse)) message("  ", nm, ": dropped by ", NORMALISATION, " as too sparse - ",
                              paste(sparse, collapse = ", "))
}

used_markers <- list()
for (nm in names(MARKERS)) {
  found <- present(MARKERS[[nm]])
  report_missing(MARKERS[[nm]], nm)
  if (length(found) < MIN_MARKERS) {
    message("  ", nm, ": skipped, only ", length(found), " of ",
            length(MARKERS[[nm]]), " markers detected")
    next
  }
  used_markers[[nm]] <- found
  obj <- AddModuleScore(obj, features = list(found), name = paste0(nm, "_"),
                        ctrl = 50, seed = 42)
  # AddModuleScore appends an index; rename to the plain cell-type name.
  names(obj@meta.data)[names(obj@meta.data) == paste0(nm, "_1")] <- nm
}
score_cols <- names(used_markers)

# Record which genes actually went into each score. Without this the figure
# legend says "Endothelial score" and nobody can reconstruct what that means.
write.csv(
  data.frame(set     = rep(names(used_markers), lengths(used_markers)),
             gene    = unlist(used_markers, use.names = FALSE),
             row.names = NULL),
  file.path(OUT_DIR, "markers_used.csv"), row.names = FALSE)

# ------------------------------------------------------------ annotation

# Name each spaceranger cluster after the marker set it is most enriched for.
#
# Cluster level, not bin level, and on purpose: a 16um bin spans several cells,
# so a per-bin winner-take-all label would be a hard assignment to soft data.
# Scores are z-scored across bins first, otherwise the comparison is between
# module scores on different scales rather than between cell types.
annotate_clusters <- function(obj, score_cols) {
  z  <- scale(as.matrix(obj@meta.data[, score_cols, drop = FALSE]))
  # droplevels: QC can empty a cluster entirely, and split() would then hand
  # colMeans a zero-row matrix and put a row of NaN in the annotation table.
  cl <- droplevels(factor(obj$cluster))
  m  <- t(vapply(split(seq_len(nrow(z)), cl),
                 function(i) colMeans(z[i, , drop = FALSE]),
                 numeric(length(score_cols))))
  colnames(m) <- score_cols

  ord    <- t(apply(m, 1, sort, decreasing = TRUE))
  top    <- ord[, 1]
  second <- if (ncol(ord) > 1) ord[, 2] else rep(-Inf, nrow(ord))
  best   <- colnames(m)[max.col(m, ties.method = "first")]

  lab <- ifelse(top > MIN_Z & (top - second) >= MIN_MARGIN, best, "Unassigned")
  # cbind an explicit data.frame rather than passing the matrix to data.frame():
  # that mangles the score columns into round.m..3..Endothelial and friends.
  cbind(
    data.frame(cluster = rownames(m), label = lab,
               n_bins = as.vector(table(cl)),
               top_score = round(top, 3), margin = round(top - second, 3),
               row.names = NULL, stringsAsFactors = FALSE),
    as.data.frame(round(m, 3), row.names = NULL))
}

annot <- NULL
if (length(score_cols) >= 2) {
  annot <- annotate_clusters(obj, score_cols)
  write.csv(annot, file.path(OUT_DIR, "cluster_annotation.csv"), row.names = FALSE)
  message("\nCluster annotation (mean z-scored module score per cluster):")
  print(annot[, c("cluster", "label", "n_bins", "top_score", "margin")])
  obj$celltype <- factor(annot$label[match(as.character(obj$cluster), annot$cluster)])
  message("Unassigned clusters: ", sum(annot$label == "Unassigned"), "/", nrow(annot),
          " - retune MIN_Z / MIN_MARGIN against cluster_annotation.csv if that looks wrong")
} else {
  message("Fewer than 2 marker sets scored - skipping cluster annotation")
}

# ---------------------------------------------------------------- plotting

# Shared look. Light by default: grey zeros with dark saturated signal, which
# survives both a projector and a printer. DARK_MODE is the black/magma variant.
BG_COL   <- if (DARK_MODE) "black" else "white"
FG_COL   <- if (DARK_MODE) "grey90" else "black"
EXPR_COLS <- if (DARK_MODE) {
  c("#000004", "#3B0F70", "#8C2981", "#DE4968", "#FE9F6D", "#FCFDBF")
} else {
  c("grey90", "#FEC98D", "#FE9F6D", "#F1605D", "#B63679", "#721F81", "#2D1160")
}

theme_panel <- function() {
  theme_void(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = BG_COL, colour = NA),
      panel.background = element_rect(fill = BG_COL, colour = NA),
      legend.background = element_rect(fill = BG_COL, colour = NA),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13, colour = FG_COL),
      plot.subtitle    = element_text(hjust = 0.5, size = 10, colour = FG_COL),
      legend.title     = element_text(size = 10, colour = FG_COL),
      legend.text      = element_text(size = 9,  colour = FG_COL),
      legend.key.width = unit(0.30, "cm"),
      legend.key.height= unit(0.75, "cm")
    )
}

theme_umap <- function() {
  theme_classic(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = BG_COL, colour = NA),
      panel.background = element_rect(fill = BG_COL, colour = NA),
      legend.background = element_rect(fill = BG_COL, colour = NA),
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13, colour = FG_COL),
      plot.subtitle = element_text(hjust = 0.5, size = 10, colour = FG_COL),
      axis.title    = element_text(size = 10, colour = FG_COL),
      axis.text     = element_blank(), axis.ticks = element_blank(),
      axis.line     = element_line(colour = FG_COL),
      legend.title  = element_text(size = 10, colour = FG_COL),
      legend.text   = element_text(size = 9,  colour = FG_COL)
    )
}

# ~60k bins per panel is a 40MB vector PDF that stalls a presentation. Rasterise
# the points, keep the labels vector.
gpoint <- function(...) {
  g <- geom_point(...)
  if (RASTERISE && requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::rasterise(g, dpi = 300)
  } else {
    g
  }
}
if (RASTERISE && !requireNamespace("ggrastr", quietly = TRUE)) {
  message("ggrastr not installed - PDFs will be large vector files. ",
          "install.packages('ggrastr') to fix.")
}

# A handful of bubbles/debris bins sit far above the tissue and saturate every
# colour scale, so cap the high end at a quantile and squish the rest in. Falls
# back to the observed range when the gene is so sparse that the quantile is
# indistinguishable from the floor.
# qlo anchors the low end: 0 means "the observed minimum" (right for a gene,
# where zero really is the floor), a quantile means "everything below this is
# background" (right for a module score, which has no meaningful zero).
expr_limits <- function(v, qlo = 0, hi = 0.99) {
  v <- v[is.finite(v)]
  lo <- if (qlo <= 0) min(v) else as.numeric(quantile(v, qlo))
  h  <- as.numeric(quantile(v, hi))
  if (!is.finite(h) || h <= lo) h <- max(v)
  if (h <= lo) h <- lo + 1
  c(lo, h)
}

expr_scale <- function(v, name = NULL, qlo = 0) {
  scale_colour_gradientn(colours = EXPR_COLS, name = name,
                         limits = expr_limits(v, qlo = qlo), oob = scales::squish)
}

# Box mean over the capture grid, via an integral image so the cost does not
# grow with the radius. Bins dropped by QC leave holes: dividing by the count
# of present neighbours rather than (2r+1)^2 keeps the tissue edge from fading
# into the background, which a fixed denominator would do.
box_sum <- function(M, r) {
  nr <- nrow(M); nc <- ncol(M)
  I <- matrix(0, nr + 1L, nc + 1L)
  I[-1L, -1L] <- t(apply(apply(M, 2, cumsum), 1, cumsum))
  ii <- seq_len(nr); jj <- seq_len(nc)
  r1 <- pmax(ii - r, 1L); r2 <- pmin(ii + r, nr)
  c1 <- pmax(jj - r, 1L); c2 <- pmin(jj + r, nc)
  I[r2 + 1L, c2 + 1L] - I[r1, c2 + 1L] - I[r2 + 1L, c1] + I[r1, c1]
}

smooth_grid <- function(values, arow, acol, r = SMOOTH_RADIUS) {
  # arrow::read_parquet can hand back integer64 for the array columns, which
  # neither is.finite() nor matrix indexing handles the way you would expect.
  arow <- as.numeric(arow); acol <- as.numeric(acol)
  ok <- is.finite(values) & is.finite(arow) & is.finite(acol)
  out <- rep(NA_real_, length(values))
  if (!any(ok) || r < 1) return(values)
  i <- as.integer(arow[ok] - min(arow[ok]) + 1L)
  j <- as.integer(acol[ok] - min(acol[ok]) + 1L)
  idx <- cbind(i, j)
  S <- matrix(0, max(i), max(j)); S[idx] <- values[ok]
  N <- matrix(0, max(i), max(j)); N[idx] <- 1
  out[ok] <- box_sum(S, r)[idx] / box_sum(N, r)[idx]
  out
}

fetch_value <- function(obj, feature, assay = NULL) {
  old <- DefaultAssay(obj)
  if (!is.null(assay)) DefaultAssay(obj) <- assay
  v <- FetchData(obj, vars = feature)[, 1]
  DefaultAssay(obj) <- old
  v
}

# Draws high values last so signal is never buried under zero-valued bins.
spatial_plot <- function(obj, feature, title = feature, subtitle = NULL,
                         assay = NULL, legend_name = NULL, qlo = 0,
                         smooth = FALSE) {
  value <- fetch_value(obj, feature, assay)
  # Smoothed panels are drawn slightly larger: at POINT_SIZE the bins leave
  # visible gaps, which stipples something meant to read as a continuous map.
  psize <- POINT_SIZE * if (smooth) 1.4 else 1
  if (smooth) value <- smooth_grid(value, obj$arow, obj$acol)
  df <- data.frame(x = obj$x, y = obj$y, value = value)
  df <- df[complete.cases(df), ]
  df <- df[order(df$value), ]
  ggplot(df, aes(x, y, colour = value)) +
    gpoint(size = psize, shape = 16) +
    expr_scale(df$value, legend_name, qlo = qlo) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = title, subtitle = subtitle) +
    theme_panel()
}

umap_plot <- function(obj, feature, title = feature, subtitle = NULL,
                      assay = NULL, legend_name = NULL, qlo = 0) {
  df <- data.frame(x = obj$UMAP_1, y = obj$UMAP_2, value = fetch_value(obj, feature, assay))
  df <- df[complete.cases(df), ]
  df <- df[order(df$value), ]
  ggplot(df, aes(x, y, colour = value)) +
    gpoint(size = POINT_SIZE, shape = 16) +
    expr_scale(df$value, legend_name, qlo = qlo) +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
    theme_umap()
}

# Okabe-Ito: colourblind-safe and legible at projector contrast, unlike the
# default hue wheel once six clusters are on screen.
CT_COLS <- c(
  Cardiomyocyte = "#E69F00", Endothelial = "#0072B2", Fibroblast = "#009E73",
  VSMC          = "#CC79A7", Myeloid     = "#D55E00", Lymphoid   = "#56B4E9",
  Unassigned    = "grey75"
)

discrete_plot <- function(df, title, subtitle = NULL, umap = FALSE, cols = NULL) {
  df <- df[complete.cases(df), ]
  p <- ggplot(df, aes(x, y, colour = group)) +
    gpoint(size = POINT_SIZE, shape = 16) +
    guides(colour = guide_legend(override.aes = list(size = 3.5), title = NULL)) +
    labs(title = title, subtitle = subtitle)
  if (!is.null(cols)) p <- p + scale_colour_manual(values = cols, drop = FALSE)
  if (umap) {
    p + labs(x = "UMAP 1", y = "UMAP 2") + theme_umap()
  } else {
    p + scale_y_reverse() + coord_fixed() + theme_panel()
  }
}

save_panel <- function(p, name, w = 6, h = 6) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = BG_COL)
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, bg = BG_COL)
}

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork) else
  message("patchwork not installed - only single panels will be written. ",
          "install.packages('patchwork') for the combined figures.")

# Every multi-panel figure goes through here, so a missing patchwork degrades to
# "some figures skipped" rather than an error two thirds of the way through.
save_combo <- function(plots, name, ncol, w, h, title = NULL, subtitle = NULL) {
  if (!has_patchwork) return(invisible(NULL))
  p <- patchwork::wrap_plots(plots, ncol = ncol)
  if (!is.null(title)) {
    p <- p + patchwork::plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme(plot.background = element_rect(fill = BG_COL, colour = NA),
                    plot.title    = element_text(face = "bold", size = 15, colour = FG_COL),
                    plot.subtitle = element_text(size = 10, colour = FG_COL)))
  }
  save_panel(p, name, w = w, h = h)
}

CAPTION_BIN <- paste0(SAMPLE, " mouse heart - Visium HD, ",
                      sub("square_0*", "", BIN), " bins")
# Bin edge in microns, for the smoothing kernel caption. Parsed rather than
# hardcoded so switching BIN to square_008um does not silently mislabel it.
BIN_UM <- as.integer(sub("^[^0-9]*0*([0-9]+)um$", "\\1", BIN))
if (is.na(BIN_UM)) stop("Cannot parse a bin size out of BIN = ", BIN)
# Averaging compresses the range, so a smoothed panel and a raw one carry
# different scales. Name them differently or the figure invites a comparison
# that is several-fold wrong.
SMOOTH_LAB <- paste0("score\n(", 2 * SMOOTH_RADIUS + 1, "x",
                     2 * SMOOTH_RADIUS + 1, " mean)")

message("Plotting ...")

# --- fig1: what the section is, before any interpretation of it -------------
orient <- present(ORIENTATION_GENES)
if (length(orient)) {
  save_combo(
    lapply(orient, function(g) spatial_plot(obj, g, g, legend_name = "log expr")),
    "fig1_orientation", ncol = 2,
    w = 5.2 * min(2, length(orient)), h = 5.2 * ceiling(length(orient) / 2),
    title = paste0("Chamber orientation - ", CAPTION_BIN),
    subtitle = "Myl2 ventricular (LV/RV, septum) | Myl7 atrial | Nppa atrial + trabecular | Ttn pan-cardiomyocyte")
}
save_panel(spatial_plot(obj, "nCount_Spatial", "UMIs per bin",
                        subtitle = CAPTION_BIN, legend_name = "nCount"),
           "fig1b_depth", w = 6.5)

# --- fig2: the hero. Annotated tissue map beside the matching UMAP ----------
if (!is.null(annot)) {
  lv <- c(intersect(names(CT_COLS), levels(obj$celltype)),
          setdiff(levels(obj$celltype), names(CT_COLS)))
  obj$celltype <- factor(as.character(obj$celltype), levels = lv)
  sp <- discrete_plot(data.frame(x = obj$x, y = obj$y, group = obj$celltype),
                      "Tissue", cols = CT_COLS)
  um <- discrete_plot(data.frame(x = obj$UMAP_1, y = obj$UMAP_2, group = obj$celltype),
                      "Transcriptome (UMAP)", umap = TRUE, cols = CT_COLS)
  save_combo(list(sp, um), "fig2_annotation", ncol = 2, w = 12.5, h = 6,
             title = paste0("Cell-type composition - ", CAPTION_BIN),
             subtitle = paste("Graph-based clusters labelled by their most enriched marker-set score.",
                              "Bins span several cells, so labels are enrichment, not purity."))
}
# Unannotated cluster map as a fallback / sanity check.
save_panel(discrete_plot(data.frame(x = obj$x, y = obj$y, group = obj$cluster),
                         paste0("Graph-based clusters - ", CAPTION_BIN)),
           "fig2b_clusters", w = 7.5)

# --- fig3 / fig4: the two panels the collaborator asked for -----------------
score_fig <- list(Endothelial = "fig3_endothelial", Fibroblast = "fig4_fibroblast")
for (nm in names(score_fig)) {
  if (!nm %in% score_cols) {
    message("Skipping ", nm, " figure - score not computed")
    next
  }
  genes <- paste(used_markers[[nm]], collapse = " + ")
  save_combo(
    list(spatial_plot(obj, nm, "Tissue", legend_name = "score", qlo = SCORE_QLO),
         umap_plot(obj, nm, "Transcriptome (UMAP)", legend_name = "score", qlo = SCORE_QLO)),
    score_fig[[nm]], ncol = 2, w = 12.5, h = 6,
    title = paste0(nm, " signature - ", CAPTION_BIN),
    subtitle = paste0("Module score over ", genes,
                      "  |  colour from the ", 100 * SCORE_QLO,
                      "th percentile up; grey = lower half of bins"))

  # Smoothed companion. Same score, same figure layout, different rendering -
  # kept as its own file so the raw panel above is always the one on record.
  if (SMOOTH_RADIUS > 0) {
    save_combo(
      list(spatial_plot(obj, nm, "Tissue - smoothed", legend_name = SMOOTH_LAB,
                        qlo = SCORE_QLO, smooth = TRUE),
           umap_plot(obj, nm, "Transcriptome (UMAP)", legend_name = "score",
                     qlo = SCORE_QLO)),
      paste0(score_fig[[nm]], "_smoothed"), ncol = 2, w = 12.5, h = 6,
      title = paste0(nm, " signature - ", CAPTION_BIN),
      subtitle = paste0("Module score over ", genes, "  |  tissue panel averaged over a ",
                        2 * SMOOTH_RADIUS + 1, "x", 2 * SMOOTH_RADIUS + 1,
                        " bin box (", (2 * SMOOTH_RADIUS + 1) * BIN_UM,
                        "um) - display smoothing only"))
  }
}

# Both signatures in one 2x2, for the single-slide version.
combo_sets <- intersect(c("Endothelial", "Fibroblast"), score_cols)
if (length(combo_sets) == 2) {
  save_combo(
    list(spatial_plot(obj, "Endothelial", "Endothelial - tissue",
                      legend_name = if (SMOOTH_RADIUS > 0) SMOOTH_LAB else "score",
                      qlo = SCORE_QLO, smooth = SMOOTH_RADIUS > 0),
         umap_plot(obj, "Endothelial", "Endothelial - UMAP", legend_name = "score",
                   qlo = SCORE_QLO),
         spatial_plot(obj, "Fibroblast", "Fibroblast - tissue",
                      legend_name = if (SMOOTH_RADIUS > 0) SMOOTH_LAB else "score",
                      qlo = SCORE_QLO, smooth = SMOOTH_RADIUS > 0),
         umap_plot(obj, "Fibroblast", "Fibroblast - UMAP", legend_name = "score",
                   qlo = SCORE_QLO)),
    "fig5_EC_FB_combined", ncol = 2, w = 12.5, h = 11.5,
    title = paste0("Endothelial and fibroblast signatures - ", CAPTION_BIN),
    subtitle = if (SMOOTH_RADIUS > 0)
      paste0("Tissue panels averaged over a ", 2 * SMOOTH_RADIUS + 1, "x",
             2 * SMOOTH_RADIUS + 1, " bin box - display smoothing only") else NULL)
} else {
  message("Skipping combined EC/FB figure - missing score(s): ",
          paste(setdiff(c("Endothelial", "Fibroblast"), score_cols), collapse = ", "))
}

# --- fig6 / fig7: the individual marker genes behind the scores -------------
ec_genes <- present(MARKERS$Endothelial)
if (length(ec_genes)) {
  save_combo(lapply(ec_genes, function(g) spatial_plot(obj, g, g, legend_name = "log expr")),
             "fig6_EC_markers", ncol = 2,
             w = 5.2 * min(2, length(ec_genes)), h = 5.2 * ceiling(length(ec_genes) / 2),
             title = paste0("Endothelial marker genes - ", CAPTION_BIN))
}
fb_genes <- present(c(MARKERS$Fibroblast, FB_EXTRA))
if (length(fb_genes)) {
  save_combo(lapply(fb_genes, function(g) spatial_plot(obj, g, g, legend_name = "log expr")),
             "fig7_FB_markers", ncol = 2,
             w = 5.2 * min(2, length(fb_genes)), h = 5.2 * ceiling(length(fb_genes) / 2),
             title = paste0("Fibroblast marker genes - ", CAPTION_BIN),
             # Not "Gsn marks activated fibroblasts": Gsn is broadly expressed
             # and at baseline this panel is retracing the same ECM-rich
             # structures as Dcn, not reporting activation.
             subtitle = paste("Pdgfra and Col1a1 are low-abundance transcripts;",
                              "Dcn and Gsn are abundant secreted ECM genes and are not fibroblast-exclusive"))
}

# --- fig8: endothelial subtypes, only if they are actually there ------------
ec_sub <- present(EC_SUBTYPES)
if (length(ec_sub)) {
  save_combo(lapply(ec_sub, function(g) spatial_plot(obj, g, g, legend_name = "log expr")),
             "fig8_EC_subtypes", ncol = 3,
             w = 5.2 * min(3, length(ec_sub)), h = 5.2 * ceiling(length(ec_sub) / 3),
             title = paste0("Endothelial subtype markers - ", CAPTION_BIN),
             subtitle = "Bmx / Mecom / Gja5 arterial | Rgcc / Ackr1 venular | Fabp4 microvascular")
} else {
  message("No endothelial subtype markers detected - skipping fig8")
}

# ---------------------------------------------------------------- detection

# Report raw counts, not normalized values - the question is whether the gene
# is seen at all, and normalization can make a handful of counts look like
# meaningful signal.
# assay = "Spatial" is load-bearing: under NORMALISATION = "SCT" the default
# assay is SCT by this point, and SCT's counts layer holds corrected counts,
# not the observed UMIs.
raw <- GetAssayData(obj, assay = "Spatial", layer = "counts")
total_umi <- sum(raw)

detection_row <- function(g) {
  if (!g %in% rownames(raw)) {
    return(data.frame(gene = g, in_reference = FALSE, total_counts = NA_real_,
                      bins_detected = NA_integer_, pct_bins = NA_real_,
                      per_million = NA_real_))
  }
  v <- raw[g, ]
  data.frame(
    gene          = g,
    in_reference  = TRUE,
    total_counts  = sum(v),
    bins_detected = sum(v > 0),
    pct_bins      = round(100 * sum(v > 0) / length(v), 4),
    per_million   = round(1e6 * sum(v) / total_umi, 4)
  )
}

detection <- do.call(rbind, lapply(c(GENE_OF_INTEREST, REFERENCE_GENES), detection_row))
detection$stage     <- "post-QC (bins used for the figures)"
detection$total_umi <- total_umi

detection_all <- rbind(prefilter_detection[, names(detection)], detection)
detection_all$sample <- SAMPLE
detection_all$bin    <- BIN

write.csv(detection_all, file.path(OUT_DIR, "detection_check.csv"), row.names = FALSE)
message("\nDetection check (raw counts, ", format(total_umi, big.mark = ","), " UMIs post-QC):")
print(detection[, c("gene", "total_counts", "bins_detected", "pct_bins", "per_million")])

# Plot the gene of interest regardless, so "not detected" is visible rather
# than asserted.
#
# Presence/absence, not a colour ramp: at a handful of counts a continuous scale
# renders as an empty panel, which looks like a failed plot rather than a
# result. Positive bins are drawn on top and oversized so single bins are
# findable from the back of a room.
# Off the Spatial counts on purpose, under either normalisation: nothing here
# should flatter a handful of UMIs into something that reads as signal.
if (GENE_OF_INTEREST %in% rownames(raw)) {
  v   <- as.numeric(raw[GENE_OF_INTEREST, ])
  pct <- round(100 * mean(v > 0), 2)
  df  <- data.frame(x = obj$x, y = obj$y, pos = v > 0)
  df  <- df[complete.cases(df), ]
  df  <- df[order(df$pos), ]

  p_tissue <- ggplot(df, aes(x, y, colour = pos)) +
    gpoint(data = subset(df, !pos), size = POINT_SIZE, shape = 16) +
    gpoint(data = subset(df,  pos), size = POINT_SIZE * 4, shape = 16) +
    # Named labels, not positional: if the gene is detected in zero bins the
    # TRUE level never appears and a length-2 label vector errors out - which
    # is exactly the case this panel exists to draw.
    scale_colour_manual(values = c(`FALSE` = "grey88", `TRUE` = "#D55E00"),
                        labels = c(`FALSE` = "0 counts",
                                   `TRUE`  = paste0(GENE_OF_INTEREST, "+")),
                        name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3.5))) +
    scale_y_reverse() + coord_fixed() +
    labs(title = paste0(GENE_OF_INTEREST, " - detected bins"),
         subtitle = paste0(sum(v > 0), " / ", length(v), " bins (", pct, "%), ",
                           sum(v), " UMIs total")) +
    theme_panel()

  # The comparison is the point: "0.1% of bins" means nothing until it sits
  # beside a gene everyone agrees is expressed and one everyone agrees is rare.
  bar_df <- detection[detection$in_reference, ]
  bar_df$gene <- factor(bar_df$gene, levels = bar_df$gene[order(bar_df$pct_bins)])
  bar_df$is_goi <- bar_df$gene == GENE_OF_INTEREST
  p_bar <- ggplot(bar_df, aes(gene, pmax(pct_bins, 1e-3), fill = is_goi)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0(pct_bins, "%")), hjust = -0.15, size = 3.2, colour = FG_COL) +
    scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#D55E00"), guide = "none") +
    scale_y_log10(expand = expansion(mult = c(0, 0.35))) +
    coord_flip() +
    labs(title = "Detection rate in context", y = "% of bins with >=1 count (log scale)", x = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.background = element_rect(fill = BG_COL, colour = NA),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 13, colour = FG_COL),
          axis.text = element_text(colour = FG_COL), axis.title = element_text(colour = FG_COL))

  save_combo(list(p_tissue, p_bar), paste0("fig9_", GENE_OF_INTEREST), ncol = 2,
             w = 12.5, h = 6,
             title = paste0(GENE_OF_INTEREST, " detection - ", CAPTION_BIN),
             subtitle = "Raw UMI counts, no normalisation - the question is detection, not level")
  save_panel(p_tissue, paste0("fig9b_", GENE_OF_INTEREST, "_tissue"), w = 6.5)
} else {
  message(GENE_OF_INTEREST, " not in the filtered matrix - no panel drawn")
}

saveRDS(obj, file.path(OUT_DIR, paste0(SAMPLE, "_", BIN, "_marker_panels.rds")))
message("\nWritten to ", OUT_DIR)
print(sessionInfo())
