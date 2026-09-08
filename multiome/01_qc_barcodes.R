# ---------------------------------------------------------------------------
# 01 - Multiome QC and per-nucleus barcode export.
#
# Produces the nucleus set that everything allelic downstream is scored on, and
# the two sinto barcode files needed to split the BAMs per nucleus.
#
# THE BARCODE TRAP THIS SOLVES. A multiome nucleus has TWO barcode sequences,
# one per modality, linked 1:1 by a 10x translation table. Verified on this
# data: in per_barcode_metrics.csv, `barcode` == `gex_barcode` in 100% of rows
# and `barcode` == `atac_barcode` in 0% of them. So the GEX BAM's CB tags are
# in gex_barcode space and the ATAC BAM's are in atac_barcode space. Feeding
# one barcode list to sinto for both BAMs makes the ATAC pass match NOTHING -
# a silent complete miss, not a partial one.
#
# Both exported files therefore use the CANONICAL barcode as the sinto GROUP,
# and differ only in the lookup column:
#     sinto_<id>_rna.txt   gex_barcode  <tab>  barcode
#     sinto_<id>_atac.txt  atac_barcode <tab>  barcode
# sinto names each output BAM after the group, so both modalities yield
# <barcode>.bam and a nucleus joins across them on filename alone.
#
# WHY FIXED THRESHOLDS AND NOT ddqcR. OCM_heart/Seurat_preprocessing.R defines
# a local initialQC override but every call site is commented out; the QC that
# actually runs there is the fixed subset() below. This follows the working
# path. ddqcR is still installed if we want to revisit it.
#
# WHY THE THRESHOLDS NEED REVIEWING RATHER THAN REUSING. The OCM numbers were
# set on CellBender-filtered snRNA. This data is shallower per nucleus - median
# 1353 UMIs at 9w and 863 at 78w - so a floor tuned there can bite much harder
# here, especially at 78w. The attrition table below reports what each single
# criterion removes so the cost is visible before anything is committed to.
#
# Run in seurat_env (NOT the RNAseq env the allelic scripts use):
#   conda activate seurat_env
#   Rscript multiome/01_qc_barcodes.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

BASE <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2"
WORK <- file.path(BASE, "adult_aged_multiome")
OUT  <- file.path(WORK, "qc")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SAMPLES <- c("9w", "78w")

# ---- thresholds: REVIEW against qc_distributions_*.pdf before trusting ----
TH <- list(
  nCount_RNA_min     = 500,
  nCount_RNA_max     = 20000,
  nFeature_RNA_min   = 200,
  nFeature_RNA_max   = 5000,
  percent_mt_max     = 5,
  # ATAC floors are deliberately permissive. Sample-level FRiP is only 0.14
  # (9w) and 0.20 (78w), so a per-nucleus FRiP floor set by scRNA intuition
  # would discard most of the data.
  atac_fragments_min = 1000,
  frip_min           = 0.05
)

say <- function(...) cat(sprintf(...), "\n", sep = "")

qc_one <- function(id) {
  say("=========== %s ===========", id)
  h5  <- file.path(WORK, id, "outs", "filtered_feature_bc_matrix.h5")
  pbm <- file.path(WORK, id, "outs", "per_barcode_metrics.csv")
  stopifnot(file.exists(h5), file.exists(pbm))

  # cellranger-arc h5 is multimodal; take the GEX assay only. ATAC QC comes
  # from per_barcode_metrics.csv instead, which avoids a Signac dependency.
  mat <- Read10X_h5(h5)
  if (is.list(mat)) {
    stopifnot("Gene Expression" %in% names(mat))
    mat <- mat[["Gene Expression"]]
  }
  obj <- CreateSeuratObject(mat, project = id)
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^mt-")

  m <- read_csv(pbm, show_col_types = FALSE) %>%
    filter(is_cell == 1) %>%
    mutate(
      frip           = ifelse(atac_fragments > 0, atac_peak_region_fragments / atac_fragments, NA_real_),
      tss_frac       = ifelse(atac_fragments > 0, atac_TSS_fragments / atac_fragments, NA_real_),
      atac_mito_frac = ifelse(atac_raw_reads > 0, atac_mitochondrial_reads / atac_raw_reads, NA_real_)
    )

  # The h5 is already the filtered matrix, so its barcodes should be exactly
  # the is_cell set. Assert rather than assume - a mismatch means the two files
  # came from different runs.
  common <- intersect(colnames(obj), m$barcode)
  say("barcodes: h5 %d, is_cell %d, shared %d", ncol(obj), nrow(m), length(common))
  if (length(common) != ncol(obj)) {
    say("WARNING: %d h5 barcodes absent from the is_cell set", ncol(obj) - length(common))
  }

  md <- obj@meta.data %>%
    tibble::rownames_to_column("barcode") %>%
    inner_join(m, by = "barcode")

  # ---- distributions, so thresholds are chosen and not inherited ----
  long <- md %>%
    select(barcode, nCount_RNA, nFeature_RNA, percent.mt,
           atac_fragments, frip, tss_frac) %>%
    pivot_longer(-barcode, names_to = "metric", values_to = "value")
  pdf(file.path(OUT, sprintf("qc_distributions_%s.pdf", id)), width = 10, height = 6)
  print(
    ggplot(long, aes(value)) +
      geom_histogram(bins = 80) +
      facet_wrap(~metric, scales = "free") +
      scale_x_continuous(trans = "log1p") +
      labs(title = sprintf("%s: per-nucleus QC metrics (log1p x)", id),
           subtitle = sprintf("n = %d nuclei called by cellranger", nrow(md))) +
      theme_minimal()
  )
  dev.off()

  # ---- attrition, one criterion at a time ----
  crit <- list(
    nCount_RNA_min     = md$nCount_RNA   >= TH$nCount_RNA_min,
    nCount_RNA_max     = md$nCount_RNA   <= TH$nCount_RNA_max,
    nFeature_RNA_min   = md$nFeature_RNA >= TH$nFeature_RNA_min,
    nFeature_RNA_max   = md$nFeature_RNA <= TH$nFeature_RNA_max,
    percent_mt_max     = md$percent.mt   <= TH$percent_mt_max,
    atac_fragments_min = md$atac_fragments >= TH$atac_fragments_min,
    frip_min           = !is.na(md$frip) & md$frip >= TH$frip_min
  )
  att <- tibble(
    criterion = names(crit),
    threshold = unlist(TH[names(crit)]),
    fails     = vapply(crit, function(k) sum(!k), integer(1)),
    pct_fail  = round(100 * vapply(crit, function(k) mean(!k), numeric(1)), 2)
  )
  say("attrition per criterion (each in isolation):")
  print(as.data.frame(att), row.names = FALSE)

  keep <- Reduce(`&`, crit)
  say("PASS %d / %d nuclei (%.1f%%)", sum(keep), length(keep), 100 * mean(keep))

  pass <- md[keep, ]
  write_csv(att,  file.path(OUT, sprintf("attrition_%s.csv", id)))
  write_csv(pass, file.path(OUT, sprintf("qc_pass_%s.csv", id)))

  # ---- sinto barcode files: group is the CANONICAL barcode in both ----
  write.table(pass[, c("gex_barcode", "barcode")],
              file.path(OUT, sprintf("sinto_%s_rna.txt", id)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(pass[, c("atac_barcode", "barcode")],
              file.path(OUT, sprintf("sinto_%s_atac.txt", id)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  say("wrote sinto_%s_rna.txt and sinto_%s_atac.txt (%d nuclei each)", id, id, nrow(pass))

  tibble(sample = id, called = nrow(md), passed = nrow(pass),
         median_umi = median(md$nCount_RNA),
         median_frags = median(md$atac_fragments))
}

summ <- bind_rows(lapply(SAMPLES, qc_one))
say("")
say("=========== summary ===========")
print(as.data.frame(summ), row.names = FALSE)
write_csv(summ, file.path(OUT, "qc_summary.csv"))
say("")
say("Review qc_distributions_*.pdf and attrition_*.csv, adjust TH at the top,")
say("and re-run before splitting BAMs - the nucleus set fixes everything after it.")
