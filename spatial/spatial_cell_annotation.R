# Cell-type annotation of the Visium HD cardiac sections at 8um bins.
#
# WHAT THIS IS FOR. A high-level picture of the biology - which cell types are
# where, how much of the section each one covers, and whether the two ages look
# alike - and a per-bin label that the allelic analysis can pool on. Depth rules
# out per-cell escape (a cardiomyocyte carries ~7 informative chrX molecules), so
# the strategy is: cell IDENTITY at bin resolution, allelic QUANTITY pooled per
# type. This script provides the identity; ase_tile_sweep.R's DOMAIN_TSV hook
# consumes it.
#
# Two annotation routes, both written out so they can be compared:
#
#   marker    Module scores over curated marker sets, z-scored across bins, and
#             each bin labelled by its winning set when it wins clearly. Needs
#             nothing but the section. This is the label the exports use.
#   transfer  Seurat label transfer from the OCM snRNA-seq object (the same
#             animals' nuclei, already typed), on a leverage-score sketch and
#             projected to every bin. Optional: runs when the reference exists
#             and DO_TRANSFER is not 0. Independent of the marker sets, so
#             agreement between the two is the sanity check.
#
# WHY 8um. A 2um bin is subcellular and has too few counts to type; a 16um bin
# holds several cells and blurs capillary endothelium into the myocytes around
# it. 8um is the smallest bin at which a marker set has a few counts to score,
# and it is the finest grid spaceranger runs its own clustering on. Bins are
# still not cells - a label means "this bin is dominated by", not "this bin is".
#
# WHAT IT WRITES, in <SPATIAL_DIR>/annotation/<SAMPLE>_<BIN>/:
#   celltype_8um.tsv        barcode, array_row, array_col, x, y, nCount,
#                           celltype (marker route), celltype_ref (transfer),
#                           one z-score column per marker set
#   domain_2um.tsv          barcode, domain - every in-tissue 2um bin carrying
#                           its 8um parent's label. THIS is the DOMAIN_TSV input
#                           for ase_tile_sweep.R. Uncompressed on purpose: that
#                           script reads it with plain fread(), and the RNAseq
#                           env has no R.utils to decompress in-process.
#   composition.csv         bins, percent and mm^2 per cell type
#   marker_set_means.csv    mean z-score of every set within every label - the
#                           diagonal should dominate; off-diagonal structure
#                           says which sets are confusable at this bin size
#   fig1  tissue map coloured by cell type                         <- the picture
#   fig1b one panel per cell type, that type highlighted
#   fig2  composition bar chart
#   fig3  marker-set x label heatmap (fig for marker_set_means.csv)
#   fig4  neighbourhood fraction per non-myocyte type - where the stroma,
#         vasculature and immune cells concentrate
#   fig5  spaceranger UMAP coloured by label, if the 8um secondary analysis exists
#   fig6  label-transfer map and marker-vs-transfer agreement, if transfer ran
#   composition_all_samples.pdf   in annotation/, once more than one sample has
#                                 been run - the age comparison, descriptive only
#
# Run on the cluster, one sample per invocation:
#   sbatch ~/Postdoc/slurm/spatial_cell_annotation.slurm 9w 78w
# or by hand from adult_aged_spatial/:
#   SAMPLE=9w Rscript ~/Postdoc/spatial/spatial_cell_annotation.R

# Must be line 1 of anything that touches Seurat/Matrix on this cluster (see
# banksy_patch.R for why).
.libPaths(c("~/R/matrix-dev", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})
for (p in c("data.table", "patchwork")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("R package '", p, "' is required")
}
library(data.table)
library(patchwork)

##### ---------------------- CONFIG ---------------------- #####
SPATIAL_DIR <- Sys.getenv("SPATIAL_DIR",
  "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial")
BASE    <- SPATIAL_DIR
SAMPLE  <- Sys.getenv("SAMPLE", "9w")          # 9w = adult, 78w = aged
BIN     <- Sys.getenv("BIN", "square_008um")
OUT_DIR <- Sys.getenv("OUT_DIR", "")
if (!nzchar(OUT_DIR)) OUT_DIR <- file.path(BASE, "annotation", paste0(SAMPLE, "_", BIN))

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

# The OCM snRNA-seq object, typed in OCM_heart/Seurat_preprocessing.R. Lives one
# level above adult_aged_spatial/ in the project tree (visium_hd_test.R reads it
# as ../OCM/...). Transfer is skipped, with a message, if it is not there.
REF_RDS     <- Sys.getenv("REF_RDS", file.path(dirname(BASE), "OCM", "heart_seurat_object_SCT.rds"))
DO_TRANSFER <- env_flag("DO_TRANSFER", TRUE)
# Reference cell types to leave out of the transfer. visium_hd_test.R dropped
# "Cardiomyocytes (stressed)" and "Epicardial - Mesothelial cells"; the default
# here keeps everything, because a label that is absent from the reference can
# never be transferred and its bins are then mislabelled as the nearest thing.
REF_DROP    <- strsplit(Sys.getenv("REF_DROP", ""), "\\s*[,;]\\s*")[[1]]
SKETCH_N    <- as.integer(env_num("SKETCH_N", 50000))
TRANSFER_DIMS <- 1:30

# "bin": every bin labelled on its own scores. "cluster": spaceranger's 8um
# graph clusters are labelled by their mean scores and bins inherit the label -
# smoother, and the right choice if the per-bin map is salt-and-pepper.
ANNOTATION_LEVEL <- Sys.getenv("ANNOTATION_LEVEL", "bin")
stopifnot(ANNOTATION_LEVEL %in% c("bin", "cluster"))

# QC. 8um bins are shallow: the marker-panel script's 16um floor of 100 UMIs
# would discard most of the section here.
MIN_COUNT <- env_num("MIN_COUNT", 20)
MIN_FEAT  <- env_num("MIN_FEAT", 15)
MAX_MT    <- env_num("MAX_MT", 50)     # loose: cardiomyocytes are legitimately 30-40%

# Labelling rule, in two steps.
#
# 1. Non-myocyte sets compete: a bin takes the set with the highest z-scored
#    module score if that score clears MIN_Z, beats the runner-up by
#    MIN_MARGIN, AND the bin carries at least MIN_MARKER_UMIS raw counts on
#    that set's genes. The count gate is not optional. A z-score is spiky for a
#    sparse set - one stray Cd3e transcript in a myocyte bin is several SDs
#    above that set's mean - and without the gate a synthetic section came out
#    24% Unassigned and 6% lymphocyte. Two counts on two different markers is
#    the least that should be called evidence of a cell.
# 2. Bins no non-myocyte set claims fall to the myocyte class (ventricular or
#    atrial, whichever scores higher) provided that score is at least
#    MIN_Z_BG, and to "Unassigned" otherwise. Myocytes are the background of a
#    heart section - most of the area, and their transcripts are in every bin
#    as ambient - so "myocyte unless there is evidence otherwise" is how the
#    tissue is read, and the rule says so rather than pretending the sets are
#    symmetric.
#
# Look at fig3 and marker_set_means.csv before retuning.
MIN_Z           <- env_num("MIN_Z", 0.5)
MIN_MARGIN      <- env_num("MIN_MARGIN", 0.25)
MIN_MARKER_UMIS <- env_num("MIN_MARKER_UMIS", 2)
MIN_Z_BG        <- env_num("MIN_Z_BG", -0.5)
BG_SETS <- c("Ventricular cardiomyocyte", "Atrial cardiomyocyte")

# Marker sets, coarse on purpose: an 8um bin cannot support subtypes. Names are
# chosen to line up with the OCM snRNA labels where the two overlap, so the
# agreement heatmap reads along its diagonal.
#
# Two known ambiguities worth carrying in mind when reading the maps:
#   - Lyve1 is also on resident cardiac macrophages, so a "Lymphatic
#     endothelial" bin next to a macrophage-rich region wants a second look
#     (Prox1 and Flt4 are the discriminating ones).
#   - Ttn, Myh6 and Ryr2 are so abundant that every bin carries some; that is
#     ambient plus adjacency, and it is why the myocyte set is scored against
#     the others by z-score rather than by raw level.
MARKERS <- list(
  `Ventricular cardiomyocyte` = c("Myl2", "Myh6", "Tnnt2", "Actc1", "Ryr2", "Ttn"),
  `Atrial cardiomyocyte`      = c("Myl7", "Myl4", "Nppa", "Sln"),
  Fibroblast                  = c("Pdgfra", "Col1a1", "Col1a2", "Dcn", "Lum", "Gsn"),
  Endothelial                 = c("Cdh5", "Pecam1", "Flt1", "Kdr", "Fabp4", "Cd36"),
  Endocardial                 = c("Npr3", "Vwf", "Cgnl1", "Emcn"),
  `Lymphatic endothelial`     = c("Lyve1", "Prox1", "Flt4", "Mmrn1"),
  `Pericyte / SMC`            = c("Rgs5", "Kcnj8", "Abcc9", "Pdgfrb", "Myh11", "Acta2", "Tagln"),
  Macrophage                  = c("C1qa", "C1qb", "Cd68", "Adgre1", "F13a1", "Mrc1"),
  Lymphocyte                  = c("Cd3e", "Cd3g", "Cd79a", "Ms4a1", "Ptprc"),
  Epicardial                  = c("Msln", "Upk3b", "Wt1", "Lrrn4"),
  Adipocyte                   = c("Adipoq", "Plin1", "Cidec", "Lep"),
  `Schwann / neuronal`        = c("Plp1", "Mpz", "Kcna1", "Nrxn1")
)
MIN_MARKERS <- 2       # a set with one detected gene is that gene, not a type

# Neighbourhood radius in bins for the fig4 fraction maps. At 8um, r = 6 is a
# 104um box, about the width of five myocytes: fine enough to keep vessels and
# scars visible, coarse enough that single mislabelled bins vanish.
SMOOTH_RADIUS <- as.integer(env_num("SMOOTH_RADIUS", 6))

POINT_SIZE <- env_num("POINT_SIZE", 0.2)   # 8um bins; ~0.35 for 016um
FLIP_X <- env_flag("FLIP_X", FALSE)         # keep FALSE to match the tile maps
FLIP_Y <- env_flag("FLIP_Y", FALSE)
RASTERISE <- TRUE

# Okabe-Ito plus a few, fixed by name so a type keeps its colour across samples
# and figures. Unassigned is drawn first and light so the labelled bins sit on top.
CT_COLS <- c(
  `Ventricular cardiomyocyte` = "#E8C27A",
  `Atrial cardiomyocyte`      = "#E69F00",
  Fibroblast                  = "#009E73",
  Endothelial                 = "#0072B2",
  Endocardial                 = "#56B4E9",
  `Lymphatic endothelial`     = "#7FDBD4",
  `Pericyte / SMC`            = "#CC79A7",
  Macrophage                  = "#D55E00",
  Lymphocyte                  = "#8B0000",
  Epicardial                  = "#6A3D9A",
  Adipocyte                   = "#F0E442",
  `Schwann / neuronal`        = "#4B4B4B",
  Unassigned                  = "grey85"
)
##### ------------------------------------------------------ #####

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) message(sprintf(...))

BIN_UM <- as.integer(sub("^[^0-9]*0*([0-9]+)um$", "\\1", BIN))
if (is.na(BIN_UM)) stop("Cannot parse a bin size out of BIN = ", BIN)
bin_dir <- file.path(BASE, SAMPLE, "outs", "binned_outputs", BIN)
if (!dir.exists(bin_dir)) stop("No such bin directory: ", bin_dir)

# ---------------------------------------------------------------- input

read_positions <- function(dir) {
  pq <- file.path(dir, "spatial", "tissue_positions.parquet")
  cs <- file.path(dir, "spatial", "tissue_positions.csv")
  p <- if (file.exists(pq)) {
    if (!requireNamespace("arrow", quietly = TRUE))
      stop("tissue_positions.parquet needs the 'arrow' package")
    as.data.frame(arrow::read_parquet(pq))
  } else if (file.exists(cs)) {
    read.csv(cs)
  } else stop("No tissue_positions file under ", file.path(dir, "spatial"))
  names(p) <- tolower(names(p))
  # arrow can hand back integer64, which %/% and match() mishandle.
  for (cn in c("in_tissue", "array_row", "array_col",
               "pxl_row_in_fullres", "pxl_col_in_fullres")) {
    if (cn %in% names(p)) p[[cn]] <- as.numeric(p[[cn]])
  }
  p[p$in_tissue == 1, ]
}

msg("Reading %s ...", bin_dir)
# spaceranger writes both the h5 and the Matrix Market directory. Prefer the h5
# (one file, faster) and fall back to the directory where hdf5r is absent.
h5  <- file.path(bin_dir, "filtered_feature_bc_matrix.h5")
mtx <- file.path(bin_dir, "filtered_feature_bc_matrix")
counts <- if (file.exists(h5) && requireNamespace("hdf5r", quietly = TRUE)) {
  Read10X_h5(h5)
} else if (dir.exists(mtx)) {
  if (file.exists(h5)) msg("  hdf5r not installed - reading the Matrix Market directory instead")
  Read10X(mtx)
} else stop("No filtered_feature_bc_matrix(.h5) under ", bin_dir)
if (is.list(counts)) counts <- counts[["Gene Expression"]]
pos <- read_positions(bin_dir)

obj <- CreateSeuratObject(counts, project = SAMPLE, assay = "Spatial",
                          min.cells = 3, min.features = 1)
rm(counts)

pct_or_zero <- function(o, pattern) {
  if (!any(grepl(pattern, rownames(o)))) return(setNames(rep(0, ncol(o)), colnames(o)))
  p <- PercentageFeatureSet(o, pattern = pattern)
  v <- if (is.data.frame(p) || is.matrix(p)) setNames(p[, 1], rownames(p)) else p
  if (!is.null(names(v)) && all(colnames(o) %in% names(v))) v <- v[colnames(o)]
  v
}
obj[["percent.mt"]] <- pct_or_zero(obj, "^mt-")

keep <- obj$nCount_Spatial >= MIN_COUNT & obj$nFeature_Spatial >= MIN_FEAT &
        obj$percent.mt <= MAX_MT & colnames(obj) %in% pos$barcode
msg("QC: keeping %d / %d bins (%.1f%%) at >= %g UMIs, >= %g genes, <= %g%% mt",
    sum(keep), ncol(obj), 100 * mean(keep), MIN_COUNT, MIN_FEAT, MAX_MT)
obj <- obj[, keep]

i <- match(colnames(obj), pos$barcode)
obj$array_row <- pos$array_row[i]
obj$array_col <- pos$array_col[i]
obj$x <- pos$pxl_col_in_fullres[i] * if (FLIP_X) -1 else 1
obj$y <- pos$pxl_row_in_fullres[i] * if (FLIP_Y) -1 else 1

obj <- NormalizeData(obj, verbose = FALSE)

# ---------------------------------------------------------------- scores

present <- function(g) intersect(g, rownames(obj))
used_markers <- list()
for (nm in names(MARKERS)) {
  found  <- present(MARKERS[[nm]])
  absent <- setdiff(MARKERS[[nm]], found)
  if (length(absent)) msg("  %s: not in matrix - %s", nm, paste(absent, collapse = ", "))
  if (length(found) < MIN_MARKERS) {
    msg("  %s: SKIPPED, only %d of %d markers detected", nm, length(found), length(MARKERS[[nm]]))
    next
  }
  used_markers[[nm]] <- found
}
if (length(used_markers) < 3) stop("Fewer than three marker sets scored - nothing to annotate with")

# One AddModuleScore call for all sets, so the control-gene bins are shared.
# ctrl controls are drawn per expression bin (24 bins), so on a small feature
# space 50 is unsatisfiable; scale it to what the matrix can supply.
n_ctrl <- min(50L, max(5L, nrow(obj) %/% 24L - max(lengths(used_markers))))
if (n_ctrl < 50L) msg("  %d genes in the matrix - using %d control genes per set instead of 50", nrow(obj), n_ctrl)
obj <- AddModuleScore(obj, features = unname(used_markers), name = "set_",
                      ctrl = n_ctrl, seed = 42)
score_cols <- paste0("set_", seq_along(used_markers))
scores <- as.matrix(obj@meta.data[, score_cols])
colnames(scores) <- names(used_markers)
obj@meta.data[, score_cols] <- NULL

# z-score across bins per set. A set that does not vary at all (every marker
# at zero everywhere) cannot rank anything and is dropped rather than left as
# a column of NaN that would poison max.col().
sds <- apply(scores, 2, sd)
if (any(sds == 0)) {
  msg("  dropping zero-variance set(s): %s", paste(names(sds)[sds == 0], collapse = ", "))
  scores <- scores[, sds > 0, drop = FALSE]
}
z <- scale(scores)
set_names <- colnames(z)

write.csv(data.frame(set = rep(names(used_markers), lengths(used_markers)),
                     gene = unlist(used_markers, use.names = FALSE), row.names = NULL),
          file.path(OUT_DIR, "markers_used.csv"), row.names = FALSE)

# ------------------------------------------------------------ annotation

# Winner and runner-up per row, vectorised: 250k bins x 12 sets is too many
# rows for apply(sort).
top_two <- function(m) {
  n <- nrow(m)
  b1 <- max.col(m, ties.method = "first")
  top <- m[cbind(seq_len(n), b1)]
  m2 <- m; m2[cbind(seq_len(n), b1)] <- -Inf
  b2 <- max.col(m2, ties.method = "first")
  list(best = colnames(m)[b1], top = top, second = m2[cbind(seq_len(n), b2)])
}
label_rows <- function(m) {
  if (!ncol(m)) return(rep("Unassigned", nrow(m)))
  tt <- top_two(m)
  ok <- is.finite(tt$top) & tt$top >= MIN_Z & (tt$top - tt$second) >= MIN_MARGIN
  ifelse(!is.na(ok) & ok, tt$best, "Unassigned")
}
# zm: rows x sets of z-scores. elig: same shape, FALSE where a set has too few
# raw counts in that row to be eligible (NULL = no gate, for cluster means).
assign_labels <- function(zm, elig = NULL) {
  fg <- setdiff(colnames(zm), BG_SETS)
  bg <- intersect(colnames(zm), BG_SETS)
  zf <- zm[, fg, drop = FALSE]
  if (!is.null(elig)) zf[!elig[, fg, drop = FALSE]] <- -Inf
  lab <- label_rows(zf)
  if (length(bg)) {
    # The same count gate applies between the two myocyte classes: the
    # ventricular set is pan-myocyte apart from Myl2 and is eligible almost
    # everywhere, while the atrial set is near zero outside the atrium, so a
    # single stray Myl7 count must not flip a ventricular bin to atrial.
    zb <- zm[, bg, drop = FALSE]
    if (!is.null(elig)) zb[!elig[, bg, drop = FALSE]] <- -Inf
    bi <- max.col(zb, ties.method = "first")
    bt <- zb[cbind(seq_len(nrow(zb)), bi)]
    fall <- lab == "Unassigned" & is.finite(bt) & bt >= MIN_Z_BG
    lab[fall] <- bg[bi[fall]]
  }
  lab
}

# Raw counts on each set's genes per bin, for the eligibility gate.
raw <- GetAssayData(obj, assay = "Spatial", layer = "counts")
set_umis <- vapply(set_names, function(nm) Matrix::colSums(raw[used_markers[[nm]], , drop = FALSE]),
                   numeric(ncol(obj)))
rm(raw)
eligible <- set_umis >= MIN_MARKER_UMIS

if (ANNOTATION_LEVEL == "bin") {
  celltype <- assign_labels(z, eligible)
} else {
  cl_file <- file.path(bin_dir, "analysis", "clustering",
                       "gene_expression_graphclust", "clusters.csv")
  if (!file.exists(cl_file)) stop("ANNOTATION_LEVEL=cluster needs ", cl_file)
  cl <- read.csv(cl_file)
  names(cl) <- tolower(names(cl))
  cl_of <- cl[[setdiff(names(cl), "barcode")[1]]][match(colnames(obj), cl$barcode)]
  cl_of[is.na(cl_of)] <- "none"
  cl_f <- factor(cl_of)
  m <- t(vapply(split(seq_len(nrow(z)), cl_f),
                function(ii) colMeans(z[ii, , drop = FALSE]), numeric(ncol(z))))
  colnames(m) <- set_names
  # No count gate at cluster level: a cluster mean already averages stray
  # counts away, and the gate would only ask whether the average bin has two.
  cl_lab <- setNames(assign_labels(m), levels(cl_f))
  cl_lab["none"] <- "Unassigned"
  celltype <- unname(cl_lab[as.character(cl_f)])
  write.csv(cbind(cluster = rownames(m), label = cl_lab[rownames(m)],
                  n_bins = as.vector(table(cl_f)), round(m, 3)),
            file.path(OUT_DIR, "cluster_annotation.csv"), row.names = FALSE)
}

ct_levels <- c(intersect(names(CT_COLS), set_names), "Unassigned")
obj$celltype <- factor(celltype, levels = ct_levels)
comp <- as.data.frame(table(celltype = obj$celltype))
comp$pct  <- round(100 * comp$Freq / sum(comp$Freq), 2)
comp$mm2  <- round(comp$Freq * (BIN_UM / 1000)^2, 3)
comp$sample <- SAMPLE
names(comp)[names(comp) == "Freq"] <- "n_bins"
write.csv(comp, file.path(OUT_DIR, "composition.csv"), row.names = FALSE)
msg("\nComposition (%s, %s, %s-level labels):", SAMPLE, BIN, ANNOTATION_LEVEL)
print(comp[, c("celltype", "n_bins", "pct", "mm2")], row.names = FALSE)

# Mean z per set within each label. The diagonal is the check that the labels
# mean what they say; a large off-diagonal entry is a set that co-varies with
# another at this bin size and should be read together with it.
set_means <- t(vapply(split(seq_len(nrow(z)), obj$celltype),
                      function(ii) if (length(ii)) colMeans(z[ii, , drop = FALSE])
                                   else rep(NA_real_, ncol(z)),
                      numeric(ncol(z))))
colnames(set_means) <- set_names
write.csv(round(set_means, 3), file.path(OUT_DIR, "marker_set_means.csv"))

# -------------------------------------------------------- label transfer

obj$celltype_ref <- NA_character_
obj$celltype_ref_score <- NA_real_
transfer_ran <- FALSE
if (DO_TRANSFER && file.exists(REF_RDS)) {
  msg("\nLabel transfer from %s ...", REF_RDS)
  transfer_ran <- tryCatch({
    ref <- readRDS(REF_RDS)
    # meta.data$celltype is the annotation; Idents() on heart_seurat_object_SCT.rds
    # is a single level ("all"), so taking Idents() unconditionally transferred
    # one constant label to every bin and the 2026-09-08 run recorded
    # celltype_ref = "all" for all 238844 bins with score 1. Prefer the column,
    # fall back to Idents only when it is absent - the same order
    # OCM_heart/scDblFinder_rates.R uses.
    ref$celltype <- if ("celltype" %in% colnames(ref@meta.data)) {
      as.character(ref$celltype)
    } else {
      msg("  reference: no meta.data$celltype, falling back to Idents()")
      as.character(Idents(ref))
    }
    keep_ref <- !is.na(ref$celltype) & nzchar(ref$celltype)
    if (!all(keep_ref)) {
      msg("  reference: dropping %d nuclei with no label", sum(!keep_ref))
      ref <- ref[, keep_ref]
    }
    if (length(unique(ref$celltype)) < 2)
      stop("reference carries ", length(unique(ref$celltype)), " label level(s) (",
           paste(head(unique(ref$celltype), 3), collapse = ", "),
           "): there is nothing to transfer. Check meta.data$celltype and Idents() on ",
           REF_RDS)
    if (length(REF_DROP) && nzchar(REF_DROP[1])) {
      ref <- subset(ref, subset = !celltype %in% REF_DROP)
      msg("  reference: dropped %s", paste(REF_DROP, collapse = ", "))
    }
    DefaultAssay(ref) <- "RNA"
    ref <- NormalizeData(ref, verbose = FALSE)
    ref <- FindVariableFeatures(ref, nfeatures = 2000, verbose = FALSE)
    ref <- ScaleData(ref, verbose = FALSE)
    ref <- RunPCA(ref, npcs = max(TRANSFER_DIMS), verbose = FALSE)
    msg("  reference: %d nuclei, %d types", ncol(ref), length(unique(ref$celltype)))

    # Sketch, anchor, project - the Seurat Visium HD recipe. Anchoring on all
    # ~250k bins is hours and tens of GB; on 50k it is minutes.
    DefaultAssay(obj) <- "Spatial"
    obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- SketchData(obj, ncells = min(SKETCH_N, ncol(obj)),
                      method = "LeverageScore", sketched.assay = "sketch", verbose = FALSE)
    DefaultAssay(obj) <- "sketch"
    obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = max(TRANSFER_DIMS), reduction.name = "pca.sketch", verbose = FALSE)
    anchors <- FindTransferAnchors(reference = ref, query = obj,
                                   reference.reduction = "pca",
                                   dims = TRANSFER_DIMS, verbose = FALSE)
    preds <- TransferData(anchorset = anchors, refdata = ref$celltype,
                          dims = TRANSFER_DIMS, verbose = FALSE)
    obj <- AddMetaData(obj, preds)
    obj <- ProjectData(obj, assay = "Spatial", sketched.assay = "sketch",
                       sketched.reduction = "pca.sketch",
                       full.reduction = "full.pca.sketch", dims = TRANSFER_DIMS,
                       refdata = list(celltype_ref = "predicted.id"), verbose = FALSE)
    DefaultAssay(obj) <- "Spatial"
    if ("celltype_ref.score" %in% names(obj@meta.data))
      obj$celltype_ref_score <- obj$celltype_ref.score
    rm(ref, anchors, preds); gc(verbose = FALSE)
    # Belt and braces: a transfer that assigns one label to every bin is a
    # failure whatever caused it, and is worse than no transfer because the
    # crosstab and fig6 look like output.
    got <- unique(obj$celltype_ref[!is.na(obj$celltype_ref)])
    if (length(got) < 2)
      stop("transfer produced a single label (", paste(head(got, 1), collapse = ""),
           ") for every bin - refusing to record it")
    msg("  transferred labels: %s",
        paste(names(sort(table(obj$celltype_ref), decreasing = TRUE)), collapse = ", "))
    TRUE
  }, error = function(e) {
    msg("  label transfer FAILED - continuing with the marker route only:\n    %s",
        conditionMessage(e))
    FALSE
  })
} else if (DO_TRANSFER) {
  msg("\nNo reference at %s - skipping label transfer (set REF_RDS to point at it)", REF_RDS)
}

# ---------------------------------------------------------------- exports

meta_out <- data.frame(
  barcode = colnames(obj), array_row = obj$array_row, array_col = obj$array_col,
  x = obj$x, y = obj$y, nCount = obj$nCount_Spatial,
  celltype = as.character(obj$celltype),
  celltype_ref = obj$celltype_ref, celltype_ref_score = obj$celltype_ref_score,
  round(z, 4), check.names = FALSE, stringsAsFactors = FALSE)
fwrite(meta_out, file.path(OUT_DIR, "celltype_8um.tsv"), sep = "\t")

# Every in-tissue 2um bin inherits its parent's label. spaceranger's coarser
# bins are 4x4 blocks of 2um bins with array index = 2um index %/% 4 (the same
# rule spase_common.R uses to rebin). Matched on array index, not on barcode
# text, and checked: if the join misses, the rule is wrong for this spaceranger
# and the export must not be trusted.
k <- BIN_UM %/% 2
pos2_dir <- file.path(BASE, SAMPLE, "outs", "binned_outputs", "square_002um")
if (dir.exists(pos2_dir)) {
  # Column is k8, not `key`: that name is data.table()'s own argument and
  # would set a key instead of making a column (see tile_ratio_map.R).
  pos2 <- as.data.table(read_positions(pos2_dir))
  pos2[, k8 := paste(as.integer(array_row) %/% k, as.integer(array_col) %/% k)]
  lab8 <- data.table(k8 = paste(as.integer(obj$array_row), as.integer(obj$array_col)),
                     domain = as.character(obj$celltype))
  # Which 2um parents exist at all (in tissue at 8um), before QC removed any.
  parents <- paste(as.integer(pos$array_row), as.integer(pos$array_col))
  hit <- mean(pos2$k8 %in% parents)
  if (hit < 0.9) stop(sprintf(
    "Only %.1f%% of 2um tissue bins have an 8um parent by array index - the %%/%% %d rule does not hold here",
    100 * hit, k))
  dom <- lab8[pos2, on = "k8"]
  dom[is.na(domain), domain := "Unannotated"]        # parent failed QC
  fwrite(dom[, .(barcode, domain)], file.path(OUT_DIR, "domain_2um.tsv"), sep = "\t")
  msg("\ndomain_2um.tsv: %d 2um bins, %.1f%% labelled (parents passing QC)",
      nrow(dom), 100 * mean(dom$domain != "Unannotated"))
  msg("  use with:  DOMAIN_TSV=%s DOMAIN_KEEP='Ventricular cardiomyocyte' ... in ase_tile_sweep.R",
      file.path(OUT_DIR, "domain_2um.tsv"))
} else {
  msg("\nNo square_002um directory - domain_2um.tsv not written")
}

# ---------------------------------------------------------------- figures

gpoint <- function(...) {
  g <- geom_point(...)
  if (RASTERISE && requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(g, dpi = 300) else g
}
theme_panel <- function() {
  theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "#52514e"),
          plot.caption = element_text(size = 7, colour = "#52514e", hjust = 0),
          legend.text = element_text(size = 8), legend.title = element_text(size = 9))
}
save_panel <- function(p, name, w = 7, h = 7) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}
CAPTION <- sprintf("%s mouse heart - Visium HD, %dum bins, %s-level marker labels",
                   SAMPLE, BIN_UM, ANNOTATION_LEVEL)

df <- data.frame(x = obj$x, y = obj$y, celltype = obj$celltype)
cols <- CT_COLS[levels(df$celltype)]

# fig1: Unassigned drawn first so the labelled bins are on top.
o <- order(df$celltype == "Unassigned", decreasing = TRUE)
p1 <- ggplot(df[o, ], aes(x, y, colour = celltype)) +
  gpoint(size = POINT_SIZE, shape = 16) +
  scale_colour_manual(values = cols, drop = FALSE, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3.5), ncol = 1)) +
  scale_y_reverse() + coord_fixed() +
  labs(title = "Cell-type map", subtitle = CAPTION,
       caption = "A bin is labelled by the marker set it is most enriched for; bins span parts of several cells, so labels are dominance, not purity.") +
  theme_panel()
save_panel(p1, "fig1_celltype_map", w = 9, h = 7)

# fig1b: one facet per type, everything else grey.
types <- setdiff(levels(df$celltype), "Unassigned")
types <- types[types %in% df$celltype]
facets <- lapply(types, function(t) {
  d <- df; d$hit <- d$celltype == t; d <- d[order(d$hit), ]
  ggplot(d, aes(x, y)) +
    gpoint(data = d[!d$hit, ], colour = "grey88", size = POINT_SIZE, shape = 16) +
    gpoint(data = d[d$hit, ], colour = CT_COLS[[t]], size = POINT_SIZE * 2, shape = 16) +
    scale_y_reverse() + coord_fixed() +
    labs(title = sprintf("%s  (%.1f%%)", t, comp$pct[comp$celltype == t])) +
    theme_panel()
})
if (length(facets)) {
  nc <- min(4, length(facets))
  save_panel(wrap_plots(facets, ncol = nc) + plot_annotation(subtitle = CAPTION),
             "fig1b_celltype_facets", w = 4.2 * nc, h = 4.2 * ceiling(length(facets) / nc))
}

# fig2: composition.
p2 <- ggplot(comp[comp$n_bins > 0, ], aes(x = reorder(celltype, pct), y = pct, fill = celltype)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sprintf("%.1f%%  (%.2f mm2)", pct, mm2)), hjust = -0.1, size = 3) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  coord_flip() +
  labs(title = "Tissue composition", subtitle = CAPTION, x = NULL, y = "% of QC-passing bins") +
  theme_classic(base_size = 11)
save_panel(p2, "fig2_composition", w = 7, h = 5)

# fig3: set x label heatmap.
hm <- as.data.frame(as.table(set_means))
names(hm) <- c("label", "set", "z")
hm$label <- factor(hm$label, levels = rev(rownames(set_means)))
hm$set   <- factor(hm$set, levels = colnames(set_means))
p3 <- ggplot(hm, aes(set, label, fill = z)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.1f", z)), size = 2.6) +
  scale_fill_gradient2(low = "#184f95", mid = "white", high = "#b02a2a", midpoint = 0, name = "mean z") +
  labs(title = "Marker-set score by assigned label", subtitle = CAPTION,
       x = "marker set", y = "assigned label",
       caption = "Rows: bins carrying a label. Columns: mean z-scored module score of each set in those bins. The diagonal should dominate.") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_panel(p3, "fig3_marker_set_heatmap", w = 8, h = 6.5)

# fig4: neighbourhood fraction of each non-myocyte type. Box mean over the
# capture grid via an integral image (after spatial_marker_panels.R), divided
# by the count of present neighbours so the tissue edge does not fade.
box_sum <- function(M, r) {
  nr <- nrow(M); nc <- ncol(M)
  I <- matrix(0, nr + 1L, nc + 1L)
  I[-1L, -1L] <- t(apply(apply(M, 2, cumsum), 1, cumsum))
  ii <- seq_len(nr); jj <- seq_len(nc)
  r1 <- pmax(ii - r, 1L); r2 <- pmin(ii + r, nr)
  c1 <- pmax(jj - r, 1L); c2 <- pmin(jj + r, nc)
  I[r2 + 1L, c2 + 1L] - I[r1, c2 + 1L] - I[r2 + 1L, c1] + I[r1, c1]
}
smooth_grid <- function(values, arow, acol, r) {
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
if (SMOOTH_RADIUS > 0) {
  nbr_types <- setdiff(types, c("Ventricular cardiomyocyte", "Atrial cardiomyocyte"))
  nbr <- lapply(nbr_types, function(t) {
    f <- smooth_grid(as.numeric(obj$celltype == t), obj$array_row, obj$array_col, SMOOTH_RADIUS)
    d <- data.frame(x = obj$x, y = obj$y, f = f)
    d <- d[order(d$f), ]
    hi <- max(quantile(d$f, 0.99, na.rm = TRUE), 1e-3)
    ggplot(d, aes(x, y, colour = f)) +
      gpoint(size = POINT_SIZE * 1.4, shape = 16) +
      scale_colour_gradientn(colours = c("grey92", "#FEC98D", "#F1605D", "#721F81", "#2D1160"),
                             limits = c(0, hi), oob = scales::squish, name = "fraction") +
      scale_y_reverse() + coord_fixed() +
      labs(title = t) + theme_panel()
  })
  if (length(nbr)) {
    nc <- min(4, length(nbr))
    save_panel(wrap_plots(nbr, ncol = nc) +
                 plot_annotation(title = "Where each type concentrates",
                                 subtitle = sprintf("%s  |  fraction of bins carrying the label within a %dx%d bin box (%d um)",
                                                    CAPTION, 2 * SMOOTH_RADIUS + 1, 2 * SMOOTH_RADIUS + 1,
                                                    (2 * SMOOTH_RADIUS + 1) * BIN_UM)),
               "fig4_neighbourhood_fraction", w = 4.6 * nc, h = 4.4 * ceiling(length(nbr) / nc))
  }
}

# fig5: spaceranger's own 8um UMAP, if it ran one.
umap_file <- file.path(bin_dir, "analysis", "umap", "gene_expression_2_components", "projection.csv")
if (file.exists(umap_file)) {
  um <- read.csv(umap_file); names(um) <- tolower(names(um))
  j <- match(colnames(obj), um$barcode)
  ucols <- setdiff(names(um), "barcode")
  d <- data.frame(x = um[[ucols[1]]][j], y = um[[ucols[2]]][j], celltype = obj$celltype)
  d <- d[complete.cases(d), ]
  d <- d[order(d$celltype == "Unassigned", decreasing = TRUE), ]
  p5 <- ggplot(d, aes(x, y, colour = celltype)) +
    gpoint(size = POINT_SIZE, shape = 16) +
    scale_colour_manual(values = cols, drop = FALSE, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3.5), ncol = 1)) +
    labs(title = "spaceranger UMAP by marker label", subtitle = CAPTION, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 11) + theme(axis.text = element_blank(), axis.ticks = element_blank())
  save_panel(p5, "fig5_umap", w = 8, h = 6.5)
}

# fig6: transfer map and agreement.
if (transfer_ran) {
  d <- data.frame(x = obj$x, y = obj$y, ref = factor(obj$celltype_ref))
  d <- d[complete.cases(d), ]
  p6a <- ggplot(d, aes(x, y, colour = ref)) +
    gpoint(size = POINT_SIZE, shape = 16) +
    guides(colour = guide_legend(override.aes = list(size = 3.5), ncol = 1, title = NULL)) +
    scale_y_reverse() + coord_fixed() +
    labs(title = "Label transfer from OCM snRNA-seq", subtitle = basename(REF_RDS)) +
    theme_panel()
  tab <- table(marker = obj$celltype, transfer = obj$celltype_ref, useNA = "no")
  agree <- as.data.frame(prop.table(tab, 1))
  agree$Freq[is.nan(agree$Freq)] <- 0
  p6b <- ggplot(agree, aes(transfer, marker, fill = Freq)) +
    geom_tile() + geom_text(aes(label = ifelse(Freq >= 0.05, sprintf("%.0f", 100 * Freq), "")), size = 2.6) +
    scale_fill_gradient(low = "white", high = "#2D6E5D", limits = c(0, 1), name = "row %") +
    labs(title = "Marker label vs transferred label", x = "transferred (snRNA reference)", y = "marker route",
         caption = "Each row sums to 100 over the transferred labels. Agreement is the check; the two routes share no genes by construction.") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  save_panel(p6a | p6b, "fig6_label_transfer", w = 15, h = 6.5)
  write.csv(as.data.frame.matrix(tab), file.path(OUT_DIR, "marker_vs_transfer_crosstab.csv"))
}

# ------------------------------------------------- all samples together

# Once more than one sample has a composition.csv, draw them side by side.
# Descriptive only: one section per age, so there is no test to run.
all_comp <- Sys.glob(file.path(BASE, "annotation", paste0("*_", BIN), "composition.csv"))
if (length(all_comp) > 1) {
  ac <- rbindlist(lapply(all_comp, fread), fill = TRUE)
  ac[, sample := factor(sample, levels = c("9w", "78w", setdiff(unique(sample), c("9w", "78w"))))]
  ac[, celltype := factor(celltype, levels = names(CT_COLS))]
  pa <- ggplot(ac[n_bins > 0], aes(sample, pct, fill = celltype)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = CT_COLS, drop = FALSE, name = NULL) +
    labs(title = "Composition by sample", x = NULL, y = "% of QC-passing bins",
         caption = "n = 1 section per age. Differences are descriptive and confounded with section plane and depth.") +
    theme_classic(base_size = 11)
  pb <- ggplot(ac[n_bins > 0 & !celltype %in% c("Ventricular cardiomyocyte", "Unassigned")],
               aes(celltype, pct, fill = sample)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c(`9w` = "#4C72B0", `78w` = "#C44E52"), name = NULL) +
    labs(title = "Non-myocyte types by sample", x = NULL, y = "% of bins") +
    theme_classic(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(BASE, "annotation", paste0("composition_all_samples_", BIN, ".pdf")),
         pa | pb, width = 13, height = 5.5, bg = "white")
  fwrite(ac, file.path(BASE, "annotation", paste0("composition_all_samples_", BIN, ".csv")))
  msg("\nWrote the cross-sample composition figure for: %s", paste(unique(ac$sample), collapse = ", "))
}

# ------------------------------------------------------------- provenance

# NOT data.table(key = ..., value = ...): `key` is data.table()'s own argument.
prov <- data.table(
  k = c("script", "run_at", "sample", "bin", "annotation_level", "min_count", "min_feat",
        "max_mt", "min_z", "min_margin", "min_marker_umis", "min_z_bg", "sets_scored",
        "transfer_ran", "ref_rds", "ref_drop", "sketch_n", "smooth_radius", "n_bins_qc",
        "seurat_version", "r_version"),
  v = c("spatial/spatial_cell_annotation.R", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        SAMPLE, BIN, ANNOTATION_LEVEL, MIN_COUNT, MIN_FEAT, MAX_MT, MIN_Z, MIN_MARGIN,
        MIN_MARKER_UMIS, MIN_Z_BG, paste(set_names, collapse = ","), transfer_ran,
        if (transfer_ran) REF_RDS else "", paste(REF_DROP, collapse = ","),
        SKETCH_N, SMOOTH_RADIUS, ncol(obj),
        as.character(packageVersion("Seurat")),
        paste(R.version$major, R.version$minor, sep = ".")))
setnames(prov, c("key", "value"))
fwrite(prov, file.path(OUT_DIR, "provenance.tsv"), sep = "\t")

saveRDS(obj, file.path(OUT_DIR, paste0(SAMPLE, "_", BIN, "_annotated.rds")))
msg("\nWritten to %s", OUT_DIR)
