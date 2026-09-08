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
  obj$celltype_provisional <- best[cl]

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
