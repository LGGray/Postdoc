# ---------------------------------------------------------------------------
# Derive spatial marker panels from the snRNA-seq annotation.
#
# The hand-curated sets in spatial/spatial_cell_annotation.R carry 4-7 genes
# each. At 8um a bin holds ~120 UMIs, so a 5-gene panel rarely reaches the
# 2-UMI eligibility gate and small cells are lost: macrophages are 10.5% of
# nuclei but 0.1% of bins. Bigger, data-derived panels raise the expected
# marker UMIs per bin at any bin size, which is the one lever that does not
# trade against spatial resolution.
#
#   DOUBLET_FILE=Allelic_ratio_results/scDblFinder_per_cell.txt \
#     Rscript OCM_heart/derive_marker_panels.R [OUT_DIR] [N_PER_SET]
#
# Writes into OUT_DIR (default marker_panels/):
#   marker_panels_snRNA.csv    set,gene,rank,pct_in,pct_out,pct_vcm,mean_expr
#   marker_panels_summary.txt  genes per set and why genes were dropped
#
# Feed the CSV to the annotation with MARKER_CSV=<path>.
#
# SELECTION. Ranking on log fold change alone picks lowly expressed genes with
# large ratios, which is the opposite of what a count gate needs. So: filter
# for specificity first, then rank the survivors by abundance, so each panel is
# the most-detected genes among the specific ones.
#   1. detected in >= MIN_PCT_IN of the type's nuclei
#   2. detected in <= MAX_PCT_OUT of all other nuclei
#   3. detected in <= MAX_PCT_VCM of ventricular myocytes. This one is the
#      point of the exercise: myocyte transcript is the ambient background of
#      every bin, so a marker that myocytes also express fires everywhere and
#      is worse than useless. Not applied to the myocyte panel itself.
#   4. each gene assigned to at most one set (its best), so the sets stay
#      competitive rather than sharing genes
#   5. no mitochondrial, ribosomal, haemoglobin or predicted-gene models
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

args     <- commandArgs(TRUE)
OUT_DIR  <- if (length(args) >= 1) args[1] else "marker_panels"
N_PER_SET<- if (length(args) >= 2) as.integer(args[2]) else 50L
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_PCT_IN  <- 0.30
MAX_PCT_OUT <- 0.15
MAX_PCT_VCM <- 0.10

# snRNA label -> spatial set name. B and T cells merge: the spatial sets carry
# one Lymphocyte class because 8um cannot separate them.
MAP <- c(
  "Ventricular Cardiomyocytes"      = "Ventricular cardiomyocyte",
  "Fibroblasts"                     = "Fibroblast",
  "Endothelial cells"               = "Endothelial",
  "Macrophages"                     = "Macrophage",
  "Pericytes - Smooth muscle cells" = "Pericyte / SMC",
  "Endocardium"                     = "Endocardial",
  "Lymphatic endothelial"           = "Lymphatic endothelial",
  "Epicardial - Mesothelial cells"  = "Epicardial",
  "B cells"                         = "Lymphocyte",
  "T cells"                         = "Lymphocyte")
# "Cardiomyocytes (stressed)" is deliberately absent: it is not a spatial set.
VCM_SET <- "Ventricular cardiomyocyte"

heart <- readRDS("heart_seurat_object_SCT.rds")
lab <- if ("celltype" %in% colnames(heart@meta.data)) as.character(heart$celltype) else as.character(Idents(heart))
if (length(unique(lab)) < 2) stop("reference has one label level - check meta.data$celltype")
heart$celltype <- lab

DOUBLET_FILE <- Sys.getenv("DOUBLET_FILE", "")
if (nzchar(DOUBLET_FILE)) {
  d <- read.delim(DOUBLET_FILE, stringsAsFactors = FALSE)
  dbl <- intersect(unique(d$cell[d$class == "doublet"]), colnames(heart))
  if (!length(dbl)) stop("no doublet keys matched colnames(heart)")
  message(sprintf("excluding %d doublets of %d nuclei (%.1f%%)",
                  length(dbl), ncol(heart), 100 * length(dbl) / ncol(heart)))
  heart <- heart[, setdiff(colnames(heart), dbl)]
} else warning("DOUBLET_FILE not set - doublets are NOT excluded", immediate. = TRUE)

heart$set <- unname(MAP[heart$celltype])
heart <- heart[, !is.na(heart$set)]
message("nuclei per set:"); print(table(heart$set))

DefaultAssay(heart) <- "RNA"
heart <- NormalizeData(heart, verbose = FALSE)
Idents(heart) <- "set"

cnt  <- GetAssayData(heart, assay = "RNA", layer = "counts")
dat  <- GetAssayData(heart, assay = "RNA", layer = "data")
sets <- sort(unique(heart$set))
is_vcm <- heart$set == VCM_SET
pct_vcm <- Matrix::rowMeans(cnt[, is_vcm, drop = FALSE] > 0)

DROP <- "^mt-|^Rp[sl]|^Hb[ab]-|^Gm[0-9]+$|Rik$|^AC[0-9]|^AY[0-9]"

# Per set: detection in and out, mean normalised expression in.
tab <- do.call(rbind, lapply(sets, function(s) {
  i <- heart$set == s
  data.frame(set = s, gene = rownames(cnt),
             pct_in    = Matrix::rowMeans(cnt[, i,  drop = FALSE] > 0),
             pct_out   = Matrix::rowMeans(cnt[, !i, drop = FALSE] > 0),
             pct_vcm   = pct_vcm,
             mean_expr = Matrix::rowMeans(dat[, i, drop = FALSE]),
             row.names = NULL, stringsAsFactors = FALSE)
}))

keep <- tab$pct_in >= MIN_PCT_IN & tab$pct_out <= MAX_PCT_OUT & !grepl(DROP, tab$gene)
keep <- keep & (tab$set == VCM_SET | tab$pct_vcm <= MAX_PCT_VCM)
cand <- tab[keep, ]

# One gene, one set: keep it where it is most specific (largest pct_in - pct_out).
cand$spec <- cand$pct_in - cand$pct_out
cand <- cand[order(-cand$spec), ]
cand <- cand[!duplicated(cand$gene), ]

# Rank the survivors by abundance: for a count gate, expected UMIs is what matters.
cand <- cand[order(cand$set, -cand$mean_expr), ]
panels <- do.call(rbind, lapply(split(cand, cand$set), function(d) {
  d <- head(d, N_PER_SET); d$rank <- seq_len(nrow(d)); d
}))
panels <- panels[order(panels$set, panels$rank),
                 c("set", "gene", "rank", "pct_in", "pct_out", "pct_vcm", "mean_expr")]
write.csv(panels, file.path(OUT_DIR, "marker_panels_snRNA.csv"), row.names = FALSE)

con <- file(file.path(OUT_DIR, "marker_panels_summary.txt"), "w")
wr <- function(...) { cat(sprintf(...), file = con); cat(sprintf(...)) }
wr("derived %d-gene panels, %d sets, from %d nuclei\n", N_PER_SET, length(unique(panels$set)), ncol(heart))
wr("filters: pct_in >= %.2f, pct_out <= %.2f, pct_vcm <= %.2f (myocyte panel exempt)\n\n",
   MIN_PCT_IN, MAX_PCT_OUT, MAX_PCT_VCM)
wr("%-28s %6s %10s %10s\n", "set", "genes", "med pct_in", "med expr")
for (s in sort(unique(panels$set))) {
  d <- panels[panels$set == s, ]
  wr("%-28s %6d %10.2f %10.3f\n", s, nrow(d), median(d$pct_in), median(d$mean_expr))
}
short <- setdiff(sets, unique(panels$set))
if (length(short)) wr("\nNO PANEL (nothing cleared the filters): %s\n", paste(short, collapse = ", "))
thin <- sapply(split(panels$gene, panels$set), length)
if (any(thin < N_PER_SET))
  wr("\nfewer than %d genes: %s\n", N_PER_SET,
     paste(sprintf("%s (%d)", names(thin)[thin < N_PER_SET], thin[thin < N_PER_SET]), collapse = ", "))
wr("\nNOT DERIVED - these have no snRNA counterpart and keep their curated sets:\n")
wr("  Atrial cardiomyocyte, Adipocyte, Schwann / neuronal\n")
close(con)
message("wrote ", file.path(OUT_DIR, "marker_panels_snRNA.csv"))
