# ---------------------------------------------------------------------------
# Shared helpers for the multiome scripts, following the 00_functions.R
# convention used in OCM_heart/allelic_ratio/.
# ---------------------------------------------------------------------------

# Mean z-scored expression per marker panel, per cell.
#
# This replaces Seurat's AddModuleScore, which failed on the merged
# Signac object with "No cell overlap between new meta data and Seurat
# object" - its metadata plumbing does not survive a merge of objects that
# carry a ChromatinAssay. Computing the scores directly removes that
# dependency and keeps the cell names under our control.
#
# Each gene is z-scored ACROSS CELLS before averaging, which is what makes the
# panels comparable to each other: without it a panel of highly expressed genes
# (Tnnt2, Myh6) would outscore every other panel in every cluster purely on
# magnitude, and the argmax would be meaningless. This is not the control-gene
# binning of the Tirosh method that AddModuleScore implements, but for "which
# panel scores highest in this cluster" it is sufficient and transparent.
panel_scores <- function(obj, panels, assay = "SCT", layer = "data") {
  dat <- SeuratObject::GetAssayData(obj, assay = assay, layer = layer)
  keep_panels <- vapply(panels, function(g) length(intersect(g, rownames(dat))) > 0, logical(1))
  panels <- panels[keep_panels]
  if (!length(panels)) stop("no panel has any gene present in assay '", assay, "'")

  out <- matrix(NA_real_, nrow = ncol(dat), ncol = length(panels),
                dimnames = list(colnames(dat), names(panels)))
  for (nm in names(panels)) {
    g <- intersect(panels[[nm]], rownames(dat))
    m <- as.matrix(dat[g, , drop = FALSE])
    sds <- apply(m, 1, stats::sd)
    z <- (m - rowMeans(m)) / pmax(sds, 1e-8)
    out[, nm] <- colMeans(z)
  }
  out
}

# Attach panel scores, assign each cluster the highest-scoring panel, and
# report. Returns the object with `celltype_provisional` set.
assign_celltypes <- function(obj, panels, cluster_col = "seurat_clusters",
                             assay = "SCT", layer = "data", verbose = TRUE) {
  ps <- panel_scores(obj, panels, assay = assay, layer = layer)

  # Align explicitly rather than trusting the names to match. This is the exact
  # failure that took out AddModuleScore, so it is checked loudly instead of
  # being allowed to produce an all-NA column.
  common <- intersect(rownames(obj@meta.data), rownames(ps))
  if (!length(common)) {
    stop("cell names disagree between meta.data and assay '", assay, "'.\n",
         "  meta.data e.g.: ", paste(head(rownames(obj@meta.data), 2), collapse = ", "), "\n",
         "  assay     e.g.: ", paste(head(rownames(ps), 2), collapse = ", "))
  }
  if (length(common) < nrow(obj@meta.data) && verbose) {
    message(sprintf("panel_scores: %d of %d cells matched",
                    length(common), nrow(obj@meta.data)))
  }
  for (nm in colnames(ps)) {
    v <- rep(NA_real_, nrow(obj@meta.data)); names(v) <- rownames(obj@meta.data)
    v[common] <- ps[common, nm]
    obj@meta.data[[nm]] <- as.numeric(v)
  }

  cl <- as.character(obj@meta.data[[cluster_col]])
  means <- vapply(split(seq_len(nrow(ps)), cl[match(rownames(ps), rownames(obj@meta.data))]),
                  function(i) colMeans(ps[i, , drop = FALSE], na.rm = TRUE),
                  numeric(ncol(ps)))
  means <- t(means)                      # clusters x panels
  best <- colnames(ps)[apply(means, 1, which.max)]
  names(best) <- rownames(means)

  # as.character() and a DIRECT meta.data assignment, both load-bearing.
  # `best[cl]` inherits its names from `best`, i.e. the CLUSTER LABELS - so
  # `obj$celltype_provisional <- best[cl]` hands AddMetaData a vector named
  # "0","1","2",... which matches no cell, and fails with "No cell overlap
  # between new meta data and Seurat object". Writing into meta.data by
  # position, with names stripped, avoids AddMetaData's name matching
  # altogether. Same reason the score columns above use as.numeric().
  ct <- as.character(best[cl])
  if (all(is.na(ct))) {
    stop("every cluster failed to score; cluster_col '", cluster_col,
         "' may not hold the clusters (values seen: ",
         paste(head(unique(cl), 5), collapse = ", "), ")")
  }
  obj <- set_meta(obj, "celltype_provisional", ct)

  if (verbose) {
    message("provisional cluster -> cell type:")
    print(data.frame(cluster = rownames(means), celltype = best[rownames(means)],
                     n = as.integer(table(cl)[rownames(means)]), row.names = NULL))
  }
  # @misc, not attr(): Seurat objects are S4 and a plain attribute can be
  # dropped silently by any method that reconstructs the object.
  obj@misc[["panel_cluster_means"]] <- means
  obj
}

# Assign a metadata column BY POSITION, with names stripped.
#
# Route every metadata write through this. Seurat's `obj$col <- v` goes through
# AddMetaData, which matches on names whenever v has them - and in R almost
# every lookup produces names you did not ask for:
#   best[cl]                    is named by CLUSTER LABEL
#   CELLTYPE_SHORT[x]           is named by FULL CELL TYPE
#   ifelse(test, yes, no)       inherits the names of `test`
# Each of those looks like a plain character vector and fails identically with
# "No cell overlap between new meta data and Seurat object". That error has now
# cost three job submissions from three different call sites, so the guard
# belongs at the boundary rather than in each function.
set_meta <- function(obj, name, values) {
  n <- nrow(obj@meta.data)
  if (length(values) != n) {
    stop("set_meta('", name, "'): ", length(values), " values for ", n, " cells")
  }
  obj@meta.data[[name]] <- unname(values)
  obj
}

# ---------------------------------------------------------------------------
# Sample order
# ---------------------------------------------------------------------------

# Adult before aged, everywhere. This has to be declared because the default is
# WRONG rather than merely arbitrary: sorted alphabetically "78w" comes before
# "9w", since "7" precedes "9" as a character. So any axis, legend or facet left
# to ggplot's own ordering silently presents aged first and reads as though the
# comparison runs backwards. Route every sample column through as_sample().
SAMPLE_LEVELS <- c("9w", "78w")
as_sample <- function(x) factor(as.character(x), levels = SAMPLE_LEVELS)

# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

# cairo_pdf when available, plain pdf otherwise.
#
# pdf() writes Type-1 fonts in latin1 and cannot measure anything outside that
# encoding, which is where "font width unknown for character 0x09 in encoding
# latin1" comes from and why some labels came out mispositioned. cairo_pdf is
# UTF-8 throughout.
#
# onefile = TRUE is NOT optional: cairo_pdf defaults to FALSE, which treats the
# filename as a per-page pattern and would silently split every multi-page
# figure into separate files.
dev_open <- function(path, width = 9, height = 7) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(path, width = width, height = height, onefile = TRUE)
  } else {
    message("cairo unavailable; falling back to pdf() - expect latin1 font warnings")
    grDevices::pdf(path, width = width, height = height)
  }
}

# Short display labels. Full names stay in `celltype_provisional` and in every
# exported table; these exist only so a panel small enough to sit three-across
# is legible. "Pericytes - Smooth muscle cells" is 31 characters and was
# colliding with three other labels.
CELLTYPE_SHORT <- c(
  "Ventricular Cardiomyocytes"      = "Ventricular CM",
  "Cardiomyocytes (stressed)"       = "CM (stressed)",
  "Fibroblasts"                     = "Fibroblasts",
  "Endothelial cells"               = "Endothelial",
  "Endocardium"                     = "Endocardium",
  "Lymphatic endothelial"           = "Lymphatic EC",
  "Pericytes - Smooth muscle cells" = "Pericyte/SMC",
  "Macrophages"                     = "Macrophages",
  "T cells"                         = "T cells",
  "B cells"                         = "B cells",
  "Epicardial - Mesothelial cells"  = "Epicardial/Meso"
)
short_labels <- function(x) {
  out <- CELLTYPE_SHORT[as.character(x)]
  # unname() is load-bearing: `out` is named by the FULL cell type, and ifelse
  # inherits the names of its `test`. See set_meta above.
  unname(ifelse(is.na(out), as.character(x), out))   # unmapped labels pass through
}

# Okabe-Ito, a published palette designed and tested for colour-vision
# deficiency. Seurat's default is evenly-spaced HCL hues, which are not
# CVD-safe - a real problem for a figure shown to a room. Fixed order, never
# cycled: past 8 categories this returns NULL and the caller keeps Seurat's
# default rather than inventing a 9th hue. Yellow is placed last because it is
# the weakest of the eight against a white scatter background.
#
# NO GREY AMONG THE CATEGORIES. Grey reads as "missing/other/inactive" to
# anyone used to looking at these plots, so a grey cluster is ambiguous no
# matter how it is labelled. Keeping grey OUT of the categorical set is what
# makes it unambiguous as na.value below - the convention is only useful if it
# is reserved. Okabe-Ito's own eighth colour is black, used here in its place.
OKABE_ITO <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
               "#56B4E9", "#D55E00", "#000000", "#F0E442")

celltype_scale <- function(levels, aes = c("colour", "fill")) {
  aes <- match.arg(aes)
  levels <- unique(as.character(levels))
  if (length(levels) > length(OKABE_ITO)) {
    message(sprintf("%d categories exceeds the %d-colour fixed order; keeping the default scale",
                    length(levels), length(OKABE_ITO)))
    return(NULL)
  }
  vals <- setNames(OKABE_ITO[seq_along(levels)], levels)
  if (aes == "colour") ggplot2::scale_colour_manual(values = vals, na.value = "grey80")
  else                 ggplot2::scale_fill_manual(values = vals, na.value = "grey80")
}
