# ---------------------------------------------------------------------------
# Build the barcode -> cell type group maps used to pseudobulk the OCM hearts,
# with scDblFinder doublets excluded.
#
# Run from the OCM/ directory:
#   Rscript $POSTDOC_ROOT/OCM_heart/pseudobulk_celltype_ids.R
#
# Input:  Allelic_ratio_results/scDblFinder_per_cell.txt
#         (written by OCM_heart/scDblFinder_rates.R; columns
#          cell / sample / celltype / class / score, cell = <sample>_<barcode>)
#
# This table -- not heart_seurat_object_SCT.rds -- is the input on purpose.
# It already carries both the singlet/doublet call and the final Idents()
# cell type for every nucleus, so the group map is built from the same
# doublet call that the rest of the allelic-ratio analysis filters on
# (00_functions.R:drop_doublets), and a 6 GB object never has to be loaded.
#
# Output (all under pseudobulk_celltype/):
#   sinto_celltype_<cohort>.txt   barcode <TAB> group label, no header.
#                                 Doublets are already gone, so neither
#                                 route below can ever put a doublet read
#                                 into a pseudobulk BAM.
#                                 - `sinto filterbarcodes --cells` reads this
#                                   directly (MODE=sinto)
#                                 - the merge route awk's column 2 out of it
#                                   to build per-group BAM lists (MODE=merge)
#   celltype_labels.txt           label <TAB> celltype, for reading figures
#   pseudobulk_cell_counts.txt    cohort/celltype/label/n_cells -- carried
#                                 into the analysis, because a pseudobulk of
#                                 5 nuclei and one of 4000 are not comparable
#                                 measurements of the same thing.
#
# The label is the cell type with every run of non-alphanumeric characters
# collapsed to "_": three of the eleven OCM cell types contain spaces,
# parentheses or hyphens ("Pericytes - Smooth muscle cells",
# "Cardiomyocytes (stressed)"), and the label becomes a BAM *filename* that
# Allelome.PRO2.sh then interpolates unquoted into its job directory name.
# ---------------------------------------------------------------------------

COHORTS  <- c("9w", "78w", "Sham", "TAC")
OUT_DIR  <- "pseudobulk_celltype"
DBL_FILE <- Sys.getenv("DOUBLET_FILE",
                       "Allelic_ratio_results/scDblFinder_per_cell.txt")

# Groups below this are still written out -- dropping them here would hide
# them, and the read depth of a pseudobulk is what actually decides whether
# its allelic ratio is usable. This only controls a warning, and the number
# is carried to the analysis in pseudobulk_cell_counts.txt.
MIN_CELLS <- as.integer(Sys.getenv("MIN_CELLS", "20"))

# Set ALLOW_NO_DOUBLETS=1 only if the input table is already singlet-only.
ALLOW_NO_DOUBLETS <- nzchar(Sys.getenv("ALLOW_NO_DOUBLETS", ""))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(DBL_FILE)) {
  stop("doublet/celltype table not found: ", DBL_FILE,
       "\n  Run OCM_heart/scDblFinder_rates.R first, or point DOUBLET_FILE at it.")
}

cells <- read.delim(DBL_FILE, header = TRUE, stringsAsFactors = FALSE)
need  <- c("cell", "sample", "celltype", "class")
if (!all(need %in% names(cells))) {
  stop(DBL_FILE, " needs columns ", paste(need, collapse = ", "),
       "; has ", paste(names(cells), collapse = ", "))
}

cells$sample   <- as.character(cells$sample)
cells$celltype <- as.character(cells$celltype)
cells$class    <- as.character(cells$class)

# ---- doublet removal ------------------------------------------------------
n_dbl <- sum(cells$class == "doublet")
if (n_dbl == 0 && !ALLOW_NO_DOUBLETS) {
  stop(DBL_FILE, " has no rows with class == 'doublet'.\n",
       "  Pseudobulking this would produce BAMs that still contain doublet ",
       "reads while being named as if they did not.\n",
       "  If the table really is singlet-only already, re-run with ",
       "ALLOW_NO_DOUBLETS=1.")
}

singlets <- cells[cells$class == "singlet", , drop = FALSE]
message(sprintf("doublets excluded: %d of %d nuclei (%.2f%%); %d singlets kept",
                n_dbl, nrow(cells), 100 * n_dbl / nrow(cells), nrow(singlets)))

# per-sample rate, so the numbers in the log match scDblFinder_rates_per_sample.txt
print(as.data.frame(table(sample = cells$sample, class = cells$class)))

missing_ct <- is.na(singlets$celltype) | !nzchar(singlets$celltype)
if (any(missing_ct)) {
  message(sprintf("dropping %d singlets with no cell type label", sum(missing_ct)))
  singlets <- singlets[!missing_ct, , drop = FALSE]
}

# ---- label sanitisation ---------------------------------------------------
sanitise <- function(x) {
  y <- gsub("[^A-Za-z0-9]+", "_", x)
  y <- gsub("^_+|_+$", "", y)
  y
}

ct_levels <- sort(unique(singlets$celltype))
labels <- data.frame(label    = sanitise(ct_levels),
                     celltype = ct_levels,
                     stringsAsFactors = FALSE)

if (anyDuplicated(labels$label)) {
  dup <- labels$label[duplicated(labels$label)]
  stop("cell type labels collide after sanitisation: ",
       paste(unique(dup), collapse = ", "),
       "\n  Two different cell types would be merged into one pseudobulk BAM.")
}
if (any(!nzchar(labels$label))) {
  stop("a cell type sanitises to an empty label: ",
       paste(labels$celltype[!nzchar(labels$label)], collapse = ", "))
}

singlets$label <- labels$label[match(singlets$celltype, labels$celltype)]

write.table(labels, file.path(OUT_DIR, "celltype_labels.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- per-cohort sinto maps ------------------------------------------------
counts <- list()

for (co in COHORTS) {
  sub <- singlets[singlets$sample == co, , drop = FALSE]
  if (!nrow(sub)) {
    warning("no singlets for cohort ", co, " -- no map written")
    next
  }

  # `cell` is <sample>_<barcode>; sinto and the per-cell BAM names both use
  # the bare cellranger barcode, so strip exactly the leading "<cohort>_".
  prefix  <- paste0("^", co, "_")
  has_pre <- grepl(prefix, sub$cell)
  if (!all(has_pre)) {
    stop(co, ": ", sum(!has_pre), " cells do not start with '", co, "_' (e.g. ",
         sub$cell[!has_pre][1], "). The barcode cannot be recovered safely.")
  }
  barcode <- sub(prefix, "", sub$cell)

  if (anyDuplicated(barcode)) {
    stop(co, ": duplicated barcodes after stripping the cohort prefix")
  }

  map <- data.frame(barcode, sub$label, stringsAsFactors = FALSE)
  # tab separated, no header, no quotes: sinto filterbarcodes' --cells format
  write.table(map, file.path(OUT_DIR, paste0("sinto_celltype_", co, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  tab <- as.data.frame(table(label = sub$label), stringsAsFactors = FALSE)
  counts[[co]] <- data.frame(
    cohort   = co,
    label    = tab$label,
    celltype = labels$celltype[match(tab$label, labels$label)],
    n_cells  = as.integer(tab$Freq),
    stringsAsFactors = FALSE
  )

  message(sprintf("%s: %d singlet nuclei -> %d groups", co, nrow(sub), nrow(tab)))
}

if (!length(counts)) stop("no cohort produced a map")

cell_counts <- do.call(rbind, counts)
cell_counts <- cell_counts[order(cell_counts$cohort, -cell_counts$n_cells), ]
write.table(cell_counts, file.path(OUT_DIR, "pseudobulk_cell_counts.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n--- nuclei per pseudobulk (doublets already removed) ---\n")
print(cell_counts, row.names = FALSE)

thin <- cell_counts[cell_counts$n_cells < MIN_CELLS, , drop = FALSE]
if (nrow(thin)) {
  cat(sprintf("\nNOTE: %d of %d pseudobulks have < %d nuclei. They are still\n",
              nrow(thin), nrow(cell_counts), MIN_CELLS))
  cat("written and still scored; the analysis filters on n_cells and on read\n")
  cat("depth, which is the quantity that actually matters:\n")
  print(thin, row.names = FALSE)
}

cat(sprintf("\nWrote maps for %d cohorts to %s/\n", length(counts), OUT_DIR))
