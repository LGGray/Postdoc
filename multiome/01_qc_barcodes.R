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
# RNA cutoffs come from ddqcR and ATAC floors are applied separately; see the
# "QC method" section below for why, and for what the inherited fixed cutoffs
# did to this data.
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

# ---- QC method ----
#
# RNA cutoffs come from ddqcR, which is the right tool here for a specific
# reason rather than a general preference. Cardiomyocytes are the most
# mitochondria-rich cell type in the heart, so ANY global mito cutoff - the
# fixed 5% or a global MAD - preferentially deletes the cell type this project
# is actually about. ddqcR clusters first and takes MAD-based cutoffs WITHIN
# each cluster, so a legitimately high-mito cardiomyocyte population is judged
# against other cardiomyocytes instead of against fibroblasts.
#
# History worth not repeating: the fixed `percent.mt < 5` inherited from
# OCM_heart/Seurat_preprocessing.R removed 84.8% of 9w and 83.9% of 78w nuclei
# on its own, more than every other criterion combined. That number was
# calibrated on CellBender-FILTERED input, and cardiac ambient RNA is heavily
# mitochondrial, so CellBender was stripping most of that signal before the 5%
# test ever applied. It was never transferable to un-denoised counts.
#
# ddqcR is RNA-only, so the ATAC floors below are applied independently and
# intersected. They also act as the backstop for ddqcR's known failure mode:
# per-cluster MAD is relative, so a cluster consisting ENTIRELY of poor
# droplets can survive it. atac_fragments and FRiP measure quality absolutely.
ATAC_TH <- list(
  atac_fragments_min = 1000,
  frip_min           = 0.05
)
# Reference-only, to quantify what ddqcR buys against what we had. Never
# applied - reported in the log so the choice of method is defensible.
REF <- list(percent_mt_fixed = 5, mads = 3)

say <- function(...) cat(sprintf(...), "\n", sep = "")

if (!requireNamespace("ddqcR", quietly = TRUE)) {
  stop("ddqcR is not installed in this environment.\n",
       "  It is used by pSS_preprocessing.R, so it exists somewhere - check\n",
       "  which env that runs in, or install into seurat_env. Do NOT silently\n",
       "  fall back to a global cutoff: that is the failure this script exists\n",
       "  to avoid.")
}
suppressPackageStartupMessages(library(ddqcR))

# Local override of ddqcR::initialQC, carried over from
# OCM_heart/Seurat_preprocessing.R. The packaged version hardcodes the HUMAN
# "MT-" prefix; against a mouse object it matches nothing, percent.mt comes
# back 0 for every nucleus, and the mito criterion then passes everything
# silently. pSS_preprocessing.R can use the stock function because its data is
# human.
initialQC <- function(data,
                      basic.n.genes = 100,
                      basic.percent.mt = 80,
                      mt.prefix = "^mt-",
                      rb.prefix = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa") {
  mt.features <- grep(mt.prefix, rownames(data), ignore.case = TRUE, value = TRUE)
  rb.features <- grep(rb.prefix, rownames(data), ignore.case = TRUE, value = TRUE)
  data[["percent.mt"]] <- if (length(mt.features) > 0) {
    PercentageFeatureSet(data, features = mt.features)
  } else rep(0, ncol(data))
  data[["percent.rb"]] <- if (length(rb.features) > 0) {
    PercentageFeatureSet(data, features = rb.features)
  } else rep(0, ncol(data))
  subset(data, subset = nFeature_RNA >= basic.n.genes & percent.mt <= basic.percent.mt)
}

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

  # Confirm the pattern matches real genes. If it matched nothing, percent.mt
  # would be 0 everywhere and the criterion would silently pass all nuclei -
  # the opposite failure to the one seen, but just as quiet.
  mt_genes <- grep("^mt-", rownames(obj), value = TRUE)
  say("mito genes matched by '^mt-': %d [%s]", length(mt_genes),
      paste(head(mt_genes, 13), collapse = ", "))
  if (length(mt_genes) == 0) say("WARNING: no mito genes matched; percent.mt is meaningless")

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

  # ---- quantiles, so a threshold can be picked as a number ----
  QS <- c(.01, .05, .10, .25, .50, .75, .90, .95, .99)
  metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt",
               "atac_fragments", "frip", "tss_frac")
  qtab <- t(vapply(metrics, function(v) quantile(md[[v]], QS, na.rm = TRUE),
                   numeric(length(QS))))
  colnames(qtab) <- paste0("p", round(100 * QS))
  say("quantiles across the %d called nuclei:", nrow(md))
  print(round(qtab, 3))
  write.csv(qtab, file.path(OUT, sprintf("quantiles_%s.csv", id)))

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

  # ---- RNA: ddqcR, per-cluster data-driven cutoffs ----
  n_called <- ncol(obj)
  obj <- initialQC(obj)
  say("initialQC (n.genes >= 100, percent.mt <= 80): %d -> %d nuclei",
      n_called, ncol(obj))

  # Called bare, exactly as pSS_preprocessing.R does. ddqc.metrics takes
  # tuning arguments (MAD threshold, clustering resolution) but the package is
  # not installed locally to check their names against, and a wrong argument
  # name here is an error rather than a silent default. Verify with
  # ?ddqc.metrics on the cluster before tuning; the defaults are the published
  # ones and are a reasonable starting point.
  pdf(file.path(OUT, sprintf("ddqc_%s.pdf", id)), width = 10, height = 7)
  df.qc <- ddqc.metrics(obj)
  dev.off()
  obj <- filterData(obj, df.qc)
  say("ddqcR (package defaults): -> %d nuclei", ncol(obj))
  # write.csv, not write_csv: df.qc may carry its cluster ids as rownames and
  # as_tibble(rownames=) errors when they are absent. This handles both.
  write.csv(df.qc, file.path(OUT, sprintf("ddqc_metrics_%s.csv", id)))

  rna_pass <- colnames(obj)

  # ---- ATAC: absolute floors, independent of ddqcR ----
  atac_pass <- md$barcode[
    md$atac_fragments >= ATAC_TH$atac_fragments_min &
    !is.na(md$frip) & md$frip >= ATAC_TH$frip_min
  ]
  say("ATAC floors (fragments >= %d, FRiP >= %.2f): %d / %d nuclei pass",
      ATAC_TH$atac_fragments_min, ATAC_TH$frip_min, length(atac_pass), nrow(md))

  keep_bc <- intersect(rna_pass, atac_pass)
  say("BOTH: %d / %d nuclei (%.1f%%)", length(keep_bc), n_called,
      100 * length(keep_bc) / n_called)

  # ---- reference comparison, reported not applied ----
  ref_fixed <- sum(md$percent.mt <= REF$percent_mt_fixed)
  ref_mad   <- sum(md$percent.mt <= median(md$percent.mt, na.rm = TRUE) +
                     REF$mads * mad(md$percent.mt, na.rm = TRUE), na.rm = TRUE)
  say("for reference, mito criterion ALONE would keep: fixed %g%% -> %d, global %d-MAD -> %d, of %d",
      REF$percent_mt_fixed, ref_fixed, REF$mads, ref_mad, nrow(md))

  pass <- md[md$barcode %in% keep_bc, ]
  write_csv(pass, file.path(OUT, sprintf("qc_pass_%s.csv", id)))

  # ---- sinto barcode files: group is the CANONICAL barcode in both ----
  write.table(pass[, c("gex_barcode", "barcode")],
              file.path(OUT, sprintf("sinto_%s_rna.txt", id)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(pass[, c("atac_barcode", "barcode")],
              file.path(OUT, sprintf("sinto_%s_atac.txt", id)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  say("wrote sinto_%s_rna.txt and sinto_%s_atac.txt (%d nuclei each)", id, id, nrow(pass))

  tibble(sample = id, called = n_called, passed = nrow(pass),
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
