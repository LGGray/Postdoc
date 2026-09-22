# ---------------------------------------------------------------------------
# 08 - Allelic ATAC: accessibility of the inactive X, and its concordance with
#      RNA escape in the same nuclei.
#
# Consumes what slurm/multiome_allelome_atac.slurm produced. Step 7 of
# multiome/ANALYSIS_PLAN.md - the axis multiome adds and RNA alone cannot give.
#
# ---------------------------------------------------------------------------
# THIS IS NOT THE SAME QUANTITY AS 05, AND MUST NOT BE CALLED "ESCAPE"
# ---------------------------------------------------------------------------
# A1 = C57BL/6 (B6), A2 = CAST/EiJ, the same orientation 05_allelome_plots.R
# documents. Xist is deleted on B6, so CAST is the inactive X in every nucleus.
#
#   on the GEX BAM   A2/(A1+A2) on chrX is escape TRANSCRIPTION   (~3-4%)
#   on the ATAC BAM  A2/(A1+A2) on chrX is Xi ACCESSIBILITY       (~27-47%)
#
# Those are different quantities and the gap between them is the result, not a
# discrepancy: the inactive X is far more accessible than it is transcribed.
# Nothing in this script should be read as, or renamed to, escape. The column
# is `xi_access` throughout.
#
# ---------------------------------------------------------------------------
# WHY THE MAPPING-BIAS CORRECTION IS RE-ESTIMATED HERE RATHER THAN BORROWED
# ---------------------------------------------------------------------------
# ANALYSIS_PLAN.md is explicit that the bias needs re-estimating separately for
# ATAC because the read length differs. It is not a detail: GEX is 101bp STAR,
# ATAC is 50/49bp BWA, and a shorter read carries fewer SNPs and tolerates
# fewer mismatches before it fails to place, so CAST reads are lost at a
# different rate. The lambda model is the same as 06_gene_level_escape.R and
# the estimate comes from THIS table's own autosomes, under the same -q 30
# filter the chrX numbers were measured under.
#
# ---------------------------------------------------------------------------
# THE PER-NUCLEUS JOIN IS THE POINT OF HAVING DONE A MULTIOME AT ALL
# ---------------------------------------------------------------------------
# cellranger-arc writes the SAME CB tag on both BAMs - the gex_barcode - which
# slurm/multiome_allelome_atac.slurm established by sampling (25/25 CB tags
# matched gex_barcode, 0/25 matched atac_barcode). So an ATAC barcode and a GEX
# barcode with the same string are the same nucleus, and RNA output can be
# plotted against chromatin state per nucleus rather than per cell type.
#
# That also makes a zero-overlap join the one failure this script treats as
# fatal. sinto writes an empty BAM rather than erroring on a barcode that
# matches nothing, so a barcode-space mistake upstream arrives here as a silent
# absence of correlation, which looks exactly like a negative biological
# result. See the barcode section of the ATAC slurm header.
#
#   conda activate seurat_env
#   Rscript multiome/08_atac_allelome.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr); library(tibble)
})

BASE <- Sys.getenv("CLUSTER_MOUNT", "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2")
WORK <- file.path(BASE, "adult_aged_multiome")
OUT  <- Sys.getenv("FIG_OUT", file.path(WORK, "figures_atac"))
# Default to DSS, not SCRATCH, which is the opposite of 05_allelome_plots.R and
# deliberate: multiome_allelome_atac.slurm has no copy-back stage, so this tree
# only exists on DSS once it has been rsynced off scratch by hand. Pointing at
# the copy is what makes the job reproducible after the purge. ATAC_TREE
# overrides it back to scratch for a run against a tree still in flight.
# Note --export=NONE: set that INSIDE the job script, not in the submitting
# shell, or it will not arrive.
TREE <- Sys.getenv("ATAC_TREE", file.path(WORK, "multiome_allelome_atac"))
# The GEX results this is compared against, written by 05_allelome_plots.R.
RNA_DIR <- Sys.getenv("RNA_FIG", file.path(WORK, "figures_allelome"))
source(file.path(BASE, "Postdoc", "multiome", "00_helpers.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SAMPLES    <- c("9w", "78w")
ANNOT_BASE <- "chr_annotation_mm39.bed"
AUTOSOMES  <- paste0("chr", 1:19)
A1_IS      <- "B6"
# Per-nucleus gate. Higher than the 20 used for GEX in 05 because ATAC is the
# deeper modality per nucleus here - median informative chrX fragments is ~112
# against ~20-25 informative chrX molecules for GEX - so 50 costs few nuclei
# and buys a materially tighter per-nucleus proportion.
MIN_INFORMATIVE <- as.integer(Sys.getenv("MIN_INFORMATIVE", "50"))

say <- function(...) cat(sprintf(...), "\n", sep = "")

if (!dir.exists(TREE)) {
  message("no ATAC allelome tree at ", TREE,
          " - relying on cached TSVs; set ATAC_TREE to rescan")
}

# ---------------------------------------------------------------------------
# consolidate - same directory grammar as 05_allelome_plots.R
# ---------------------------------------------------------------------------
parse_group <- function(dirs) {
  x <- sub("^[0-9]{4}_[0-9]{2}_[0-9]{2}_", "", dirs)
  sub(paste0("_", ANNOT_BASE, "_[0-9]+$"), "", x)
}

load_tree <- function(tree) {
  cache <- file.path(OUT, sprintf("allelome_atac_%s.tsv", tree))
  if (file.exists(cache) && !identical(Sys.getenv("REBUILD"), "1")) {
    say("using cached %s (REBUILD=1 to rescan)", basename(cache))
    return(read_tsv(cache, show_col_types = FALSE) %>% mutate(sample = as_sample(sample)))
  }
  rows <- list()
  partial <- FALSE
  for (id in SAMPLES) {
    root  <- file.path(TREE, tree, "out", id)
    paths <- list.files(root, pattern = "^locus_table\\.txt$",
                        recursive = TRUE, full.names = TRUE)
    # Count the run directories too, and compare. An Allelome.PRO2 run directory
    # without a locus_table.txt in it is either a job that died partway or - far
    # more likely here - a tree that is still being rsynced off scratch, since
    # multiome_allelome_atac.slurm has no copy-back stage and the copy is done
    # by hand after the job. That distinction matters because of the CACHE
    # below: scanning mid-copy would write a partial consolidation to disk and
    # every later run would read it silently, so this has to be loud.
    dirs <- list.files(root, pattern = paste0("_", ANNOT_BASE, "_[0-9]+$"))
    say("%s / %s: %d locus tables in %d run directories", tree, id, length(paths), length(dirs))
    if (length(dirs) && length(paths) < 0.98 * length(dirs)) {
      partial <- TRUE
      say("  WARNING: %d of %d run directories have no locus_table.txt.",
          length(dirs) - length(paths), length(dirs))
      say("  If the rsync from scratch is still running, WAIT - and once it has")
      say("  finished re-run with REBUILD=1, or the cache written below keeps")
      say("  this partial result for good.")
    }
    if (!length(paths)) next
    grp <- parse_group(basename(dirname(paths)))
    for (i in seq_along(paths)) {
      d <- tryCatch(read.delim(paths[i], stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(d) || !nrow(d) || !all(c("chr","A1_reads","A2_reads") %in% names(d))) next
      d <- aggregate(cbind(A1_reads, A2_reads) ~ chr, data = d, FUN = sum)
      d$sample <- id; d$group <- grp[i]
      rows[[length(rows) + 1]] <- d
    }
  }
  if (!length(rows)) return(NULL)
  out <- bind_rows(rows) %>% mutate(total = A1_reads + A2_reads)
  if (partial) {
    # Still return the data so the run is useful, but do NOT persist it. A
    # cache is only safe when it is a complete scan, and an incomplete one is
    # indistinguishable from a complete one once written.
    say("NOT caching %s - the scan was incomplete (see the warning above)",
        basename(cache))
  } else {
    write_tsv(out, cache)
    say("wrote %s (%d rows)", basename(cache), nrow(out))
  }
  mutate(out, sample = as_sample(sample))
}

ct <- load_tree("celltype")
pn <- load_tree("pernucleus")
if (is.null(ct) && is.null(pn)) stop("no ATAC Allelome.PRO2 output found under ", TREE)

# A per-nucleus tree that is present for one sample only is the expected state
# partway through, and is worth naming rather than letting it turn into a
# one-sample figure without comment.
if (!is.null(pn)) {
  have <- sort(unique(as.character(pn$sample)))
  if (length(have) < length(SAMPLES)) {
    say("NOTE: per-nucleus ATAC present for %s only (missing %s).",
        paste(have, collapse = ", "), paste(setdiff(SAMPLES, have), collapse = ", "))
    say("  multiome_allelome_atac.slurm writes to SCRATCH and has no copy-back")
    say("  stage, so a finished sample can still be absent here. Check scratch")
    say("  before concluding the run is incomplete.")
  }
}

# ---------------------------------------------------------------------------
# mapping bias - measured on THIS modality, then applied
# ---------------------------------------------------------------------------
# p_obs = p / (p + (1-p)(1-lambda)); autosomal truth is p = 0.5, so
# lambda = 2 - 1/p_auto and p = p_obs(1-lambda) / (1 - p_obs*lambda).
# Identical algebra to 06_gene_level_escape.R, different lambda.
pool  <- bind_rows(ct, pn)
bias_src <- if (!is.null(ct)) ct else pn   # cell types are the same fragments, partitioned
auto <- bias_src %>% filter(chr %in% AUTOSOMES) %>%
  summarise(A1 = sum(A1_reads), A2 = sum(A2_reads))
p_auto <- auto$A1 / (auto$A1 + auto$A2)
if (!is.finite(p_auto) || p_auto <= 0 || p_auto >= 1) {
  stop("autosomal B6 fraction is ", p_auto, " - no usable autosomal ATAC reads. ",
       "Check the run used chr_annotation_mm39.bed and -q 30 (NOT -q 255, which ",
       "returns zero reads from a BWA BAM without erroring).")
}
LAMBDA <- 2 - 1 / p_auto
debias <- function(p_obs) p_obs * (1 - LAMBDA) / (1 - p_obs * LAMBDA)
# Accessibility of Xi, corrected. Input is the raw B6 fraction.
xi_of <- function(a1, a2) 1 - debias(a1 / (a1 + a2))

xw <- bias_src %>% filter(chr == "chrX") %>% summarise(A1 = sum(A1_reads), A2 = sum(A2_reads))
p_x <- xw$A1 / (xw$A1 + xw$A2)

say("")
say("============ ATAC orientation and mapping bias (A1 = %s) ============", A1_IS)
say("autosomal B6 fraction        %.4f  (%d / %d fragments)", p_auto, auto$A1, auto$A1 + auto$A2)
say("  -> CAST read loss lambda   %.4f", LAMBDA)
say("pooled chrX B6 fraction      %.4f  (%d / %d fragments)", p_x, xw$A1, xw$A1 + xw$A2)
say("pooled Xi accessibility      %.4f raw, %.4f bias-corrected", 1 - p_x, xi_of(xw$A1, xw$A2))
say("  (0.5 = Xi as accessible as Xa; 0 = fully closed)")
say("====================================================================")
say("")

# Sanity checks, as WARNINGS not stops. Unlike the GEX run there is no strong
# prior on where this number should land - that is the open question - so a
# hard floor like 05's CHRX_A1_MIN would be asserting the answer. What can be
# checked is the ordering: Xi cannot be more accessible than a balanced
# autosome, and chromatin that is closed cannot be transcribed.
xi_pool <- xi_of(xw$A1, xw$A2)
if (xi_pool > 0.5) {
  say("WARNING: pooled Xi accessibility is %.3f, above the 0.5 an autosome gives.", xi_pool)
  say("  That says the inactive X is MORE accessible than the active one, which")
  say("  is not a result - look first at the Xic mask extent (is _no_Xist right?)")
  say("  and at whether -q 30 was actually applied.")
}

# ---------------------------------------------------------------------------
# per cell type
# ---------------------------------------------------------------------------
wilson <- function(k, n) {
  if (n == 0) return(c(NA, NA))
  p <- k / n; z <- 1.96; d <- 1 + z^2 / n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / d,
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / d)
}

atac_ct <- NULL
if (!is.null(ct)) {
  atac_ct <- ct %>%
    mutate(region = case_when(chr == "chrX" ~ "chrX",
                              chr %in% AUTOSOMES ~ "autosomal", TRUE ~ NA_character_)) %>%
    filter(!is.na(region)) %>%
    group_by(sample, group, region) %>%
    summarise(A1 = sum(A1_reads), A2 = sum(A2_reads), .groups = "drop") %>%
    mutate(n = A1 + A2,
           xi_raw    = A2 / n,
           xi_access = xi_of(A1, A2))
  ci <- t(mapply(function(a1, a2) wilson(a2, a1 + a2), atac_ct$A1, atac_ct$A2))
  # The interval is computed on the RAW proportion then pushed through the same
  # monotone correction as the point estimate, so the bars stay consistent with
  # the point rather than being a raw interval around a corrected value.
  atac_ct$lo <- 1 - debias(1 - ci[, 1])
  atac_ct$hi <- 1 - debias(1 - ci[, 2])
  # Recover the full cell type the same way 06 does, by inverting the exact
  # gsub 03_signac_joint.R applied, rather than guessing at the punctuation.
  m <- setNames(names(CELLTYPE_SHORT), gsub("[^A-Za-z0-9]+", "_", names(CELLTYPE_SHORT)))
  full <- unname(m[atac_ct$group])                  # unname: see set_meta in 00_helpers.R
  atac_ct$celltype <- ifelse(is.na(full), atac_ct$group, full)
  atac_ct$label <- short_labels(atac_ct$celltype)

  say("per-celltype chrX Xi accessibility (bias-corrected), Wilson 95%% CI:")
  print(as.data.frame(atac_ct %>% filter(region == "chrX") %>%
    select(sample, label, n, xi_raw, xi_access, lo, hi) %>% arrange(sample, desc(n))),
    row.names = FALSE, digits = 3)
  write_tsv(atac_ct, file.path(OUT, "atac_celltype_xi.tsv"))

  dev_open(file.path(OUT, "atac_celltype_xi.pdf"), width = 11, height = 6.5)
  print(
    ggplot(atac_ct %>% filter(region == "chrX"),
           aes(reorder(label, xi_access), xi_access, colour = sample)) +
      geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey60") +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, position = position_dodge(0.4)) +
      geom_point(size = 2.6, position = position_dodge(0.4)) +
      coord_flip() + ylim(0, 0.55) +
      sample_scale("colour") +
      labs(x = NULL, y = "Xi accessibility (CAST fragment fraction, bias-corrected)",
           title = "Accessibility of the inactive X, per cell type",
           subtitle = "dotted = 0.5, i.e. Xi as accessible as Xa; n=1 animal per age, so the age contrast is descriptive") +
      theme_minimal()
  )
  print(
    ggplot(atac_ct, aes(reorder(label, xi_access), xi_access, colour = sample, shape = region)) +
      geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey60") +
      geom_point(size = 2.6, position = position_dodge(0.4)) + coord_flip() +
      sample_scale("colour") +
      labs(x = NULL, y = "CAST fragment fraction (bias-corrected)",
           title = "chrX against the autosomal control, ATAC",
           subtitle = "autosomes sit at 0.5 by construction after correction; the chrX deficit is the silencing") +
      theme_minimal()
  )
  dev.off()
  say("wrote atac_celltype_xi.pdf")
}

# ---------------------------------------------------------------------------
# RNA vs ATAC, per cell type - the headline
# ---------------------------------------------------------------------------
# The RNA side is read from 05_allelome_plots.R's export and re-corrected HERE
# with its own lambda. 05 reports raw numbers, and comparing a corrected ATAC
# value against a raw RNA one would put a few points of pure methodology into
# the gap that this figure exists to measure.
rna_f <- file.path(RNA_DIR, "celltype_escape.tsv")
if (!is.null(atac_ct) && file.exists(rna_f)) {
  rna <- read_tsv(rna_f, show_col_types = FALSE)
  rauto <- rna %>% filter(region == "autosomal") %>% summarise(A1 = sum(A1), A2 = sum(A2))
  p_auto_rna <- rauto$A1 / (rauto$A1 + rauto$A2)
  LAMBDA_RNA <- 2 - 1 / p_auto_rna
  debias_rna <- function(p) p * (1 - LAMBDA_RNA) / (1 - p * LAMBDA_RNA)

  say("")
  say("RNA lambda %.4f (from %s) vs ATAC lambda %.4f - re-estimated, not borrowed",
      LAMBDA_RNA, basename(rna_f), LAMBDA)

  cmp <- rna %>% filter(region == "chrX") %>%
    transmute(sample = as_sample(sample), group,
              rna_n = n, rna_escape = 1 - debias_rna(A1 / (A1 + A2))) %>%
    inner_join(atac_ct %>% filter(region == "chrX") %>%
                 transmute(sample, group, atac_n = n, xi_access, label),
               by = c("sample", "group"))

  if (!nrow(cmp)) {
    say("WARNING: no cell type matched between the RNA and ATAC tables. The group")
    say("  names come from the same sinto export, so a zero-row join means one")
    say("  side was built from a different tree. Skipping the concordance figure.")
  } else {
    cmp <- cmp %>% mutate(ratio = xi_access / pmax(rna_escape, 1e-9))
    write_tsv(cmp, file.path(OUT, "rna_atac_celltype.tsv"))
    say("")
    say("RNA escape vs Xi accessibility, per cell type:")
    print(as.data.frame(cmp %>% select(sample, label, rna_escape, xi_access, ratio) %>%
                          arrange(sample, desc(xi_access))), row.names = FALSE, digits = 3)

    dev_open(file.path(OUT, "rna_atac_concordance.pdf"), width = 10, height = 7)
    print(
      ggplot(cmp, aes(rna_escape, xi_access, colour = sample)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
        geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey70") +
        geom_point(aes(size = atac_n), alpha = 0.85) +
        geom_text(aes(label = label), size = 3, vjust = -1.1, show.legend = FALSE) +
        scale_size_continuous(name = "informative\nATAC fragments") +
        sample_scale("colour") +
        expand_limits(x = 0, y = 0) +
        labs(x = "RNA escape (CAST molecule fraction, bias-corrected)",
             y = "Xi accessibility (CAST fragment fraction, bias-corrected)",
             title = "The inactive X is far more accessible than it is transcribed",
             subtitle = "dashed = equality; every point above it is chromatin open without output. Both axes bias-corrected with their own lambda.") +
        theme_minimal()
    )
    print(
      ggplot(cmp, aes(reorder(label, ratio), ratio, colour = sample)) +
        geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
        geom_point(size = 2.8, position = position_dodge(0.4)) + coord_flip() +
        sample_scale("colour") +
        labs(x = NULL, y = "Xi accessibility / RNA escape",
             title = "Accessibility-to-output ratio on the inactive X",
             subtitle = "1 = chromatin state fully explains transcriptional output") +
        theme_minimal()
    )
    dev.off()
    say("wrote rna_atac_concordance.pdf")
  }
} else if (is.null(atac_ct)) {
  say("no per-celltype ATAC tree; skipping the RNA/ATAC concordance")
} else {
  say("no %s; run 05_allelome_plots.R first for the RNA/ATAC concordance", rna_f)
}

# ---------------------------------------------------------------------------
# per nucleus
# ---------------------------------------------------------------------------
if (!is.null(pn)) {
  lab <- file.path(WORK, "figures_signac", "nucleus_celltype_joint.csv")
  meta <- if (file.exists(lab)) read_csv(lab, show_col_types = FALSE) else NULL

  pnx <- pn %>% filter(chr == "chrX") %>%
    transmute(sample, barcode = group, A1 = A1_reads, A2 = A2_reads,
              n_inf = A1_reads + A2_reads, xi_access = xi_of(A1_reads, A2_reads))
  if (!is.null(meta)) {
    # as_sample() again AFTER the join: meta$sample is character from read_csv,
    # and a join between a factor and a character key resolves to character,
    # dropping the levels. facet_wrap reads the column's own order and never
    # consults sample_scale(). See SAMPLE_LEVELS in 00_helpers.R.
    pnx <- pnx %>% left_join(meta %>% select(sample, barcode, celltype_provisional),
                             by = c("sample", "barcode")) %>%
      mutate(sample = as_sample(sample),
             label = short_labels(ifelse(is.na(celltype_provisional),
                                         "unlabelled", celltype_provisional)))
  } else pnx$label <- "all nuclei"

  say("")
  say("per-nucleus ATAC chrX: %d nuclei, median informative fragments %.0f",
      nrow(pnx), median(pnx$n_inf))
  say("  nuclei at or above the %d-fragment gate: %d (%.1f%%)", MIN_INFORMATIVE,
      sum(pnx$n_inf >= MIN_INFORMATIVE), 100 * mean(pnx$n_inf >= MIN_INFORMATIVE))
  write_tsv(pnx, file.path(OUT, "atac_pernucleus_xi.tsv"))

  keep <- pnx %>% filter(n_inf >= MIN_INFORMATIVE)
  dev_open(file.path(OUT, "atac_pernucleus_xi.pdf"), width = 11, height = 7)
  print(
    ggplot(pnx, aes(n_inf, xi_access)) +
      geom_point(alpha = 0.25, size = 0.7) +
      geom_vline(xintercept = MIN_INFORMATIVE, linetype = "dashed", colour = OKABE_ITO[6]) +
      geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey60") +
      scale_x_log10() + facet_wrap(~sample) +
      labs(x = "informative chrX fragments per nucleus (log)", y = "Xi accessibility",
           title = "Per-nucleus Xi accessibility against the evidence behind it",
           subtitle = "the funnel IS the noise: low-count nuclei spread widely with no biology involved") +
      theme_minimal()
  )
  p <- ggplot(keep, aes(xi_access, colour = sample)) +
    geom_density(linewidth = 0.9) +
    geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey60") +
    sample_scale("colour") +
    labs(x = "Xi accessibility (CAST fragment fraction, bias-corrected)", y = "density",
         title = sprintf("Per-nucleus Xi accessibility (>= %d informative fragments)", MIN_INFORMATIVE),
         subtitle = "n=1 animal per age - the age contrast is descriptive only") +
    theme_minimal()
  print(p)
  print(p + facet_wrap(~label, scales = "free_y") +
          labs(title = "Per-nucleus Xi accessibility by cell type"))
  # autosomal control at nucleus level: must centre on 0.5 after correction
  pna <- pn %>% filter(chr %in% AUTOSOMES) %>%
    group_by(sample, group) %>%
    summarise(A1 = sum(A1_reads), A2 = sum(A2_reads), .groups = "drop") %>%
    mutate(n_inf = A1 + A2, ratio = 1 - debias(A1 / (A1 + A2))) %>%
    filter(n_inf >= MIN_INFORMATIVE)
  print(
    ggplot(pna, aes(ratio, colour = sample)) + geom_density(linewidth = 0.9) +
      geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey60") +
      sample_scale("colour") +
      labs(x = "autosomal CAST fragment fraction (bias-corrected)", y = "density",
           title = "Autosomal control, per nucleus, ATAC",
           subtitle = "centres on 0.5 by construction - a shifted or bimodal curve means the correction is not capturing the bias") +
      theme_minimal()
  )
  dev.off()
  say("wrote atac_pernucleus_xi.pdf")

  # -------------------------------------------------------------------------
  # the same nucleus, both modalities
  # -------------------------------------------------------------------------
  rna_pn_f <- file.path(RNA_DIR, "pernucleus_escape.tsv")
  if (file.exists(rna_pn_f)) {
    rna_pn <- read_tsv(rna_pn_f, show_col_types = FALSE)
    # Re-correct the RNA side with its own lambda, as above. LAMBDA_RNA may
    # already exist from the per-celltype block; recover it here when that block
    # was skipped, so the per-nucleus join does not silently fall back to raw.
    if (!exists("LAMBDA_RNA")) {
      rf <- file.path(RNA_DIR, "celltype_escape.tsv")
      if (file.exists(rf)) {
        rr <- read_tsv(rf, show_col_types = FALSE) %>% filter(region == "autosomal") %>%
          summarise(A1 = sum(A1), A2 = sum(A2))
        LAMBDA_RNA <- 2 - 1 / (rr$A1 / (rr$A1 + rr$A2))
      } else LAMBDA_RNA <- 0
    }
    debias_rna <- function(p) p * (1 - LAMBDA_RNA) / (1 - p * LAMBDA_RNA)

    both <- rna_pn %>%
      transmute(sample = as.character(sample), barcode,
                rna_n = n_inf, rna_escape = 1 - debias_rna(A1 / (A1 + A2))) %>%
      inner_join(pnx %>% transmute(sample = as.character(sample), barcode,
                                   atac_n = n_inf, xi_access, label),
                 by = c("sample", "barcode")) %>%
      mutate(sample = as_sample(sample))

    # THE CHECK THAT MUST BE FATAL. Both trees are keyed on the gex_barcode, so
    # a nucleus present in both files must join. Zero overlap means a
    # barcode-space error upstream - and sinto writes an EMPTY BAM rather than
    # failing on a barcode that matches nothing, so that error has no message of
    # its own and arrives here looking like an absence of correlation.
    if (!nrow(both)) {
      stop("zero barcode overlap between the RNA and ATAC per-nucleus tables.\n",
           "  Both are keyed on cellranger-arc's CB tag, which is the gex_barcode\n",
           "  for BOTH modalities - see the barcode section of\n",
           "  slurm/multiome_allelome_atac.slurm. A zero-row join therefore means\n",
           "  one side was split on sinto_<id>_atac_bycelltype.txt (the atac_barcode\n",
           "  space), which produces empty BAMs silently.\n",
           "  RNA e.g.: ", paste(head(rna_pn$barcode, 2), collapse = ", "), "\n",
           "  ATAC e.g.: ", paste(head(pnx$barcode, 2), collapse = ", "))
    }
    say("")
    say("nuclei with BOTH modalities: %d (%.0f%% of ATAC nuclei, %.0f%% of RNA nuclei)",
        nrow(both), 100 * nrow(both) / nrow(pnx), 100 * nrow(both) / nrow(rna_pn))

    bk <- both %>% filter(atac_n >= MIN_INFORMATIVE, rna_n >= 20)
    say("  passing both gates (ATAC >= %d, RNA >= 20): %d", MIN_INFORMATIVE, nrow(bk))
    write_tsv(both, file.path(OUT, "rna_atac_pernucleus.tsv"))

    if (nrow(bk) >= 30) {
      # Spearman, not Pearson. Both axes are bounded proportions on very
      # different scales with a hard floor at 0, and the RNA side is near that
      # floor for most nuclei - a rank correlation says whether nuclei ORDER
      # together without assuming the relationship is linear.
      cors <- bk %>% group_by(sample) %>%
        summarise(n = n(),
                  rho = suppressWarnings(cor(rna_escape, xi_access, method = "spearman")),
                  p = suppressWarnings(cor.test(rna_escape, xi_access,
                                                method = "spearman")$p.value),
                  .groups = "drop")
      say("")
      say("per-nucleus RNA/ATAC rank correlation:")
      print(as.data.frame(cors), row.names = FALSE, digits = 3)
      write_tsv(cors, file.path(OUT, "rna_atac_pernucleus_cor.tsv"))

      dev_open(file.path(OUT, "rna_atac_pernucleus.pdf"), width = 11, height = 6.5)
      print(
        ggplot(bk, aes(rna_escape, xi_access)) +
          geom_point(alpha = 0.3, size = 0.9) +
          geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = OKABE_ITO[2]) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
          facet_wrap(~sample) +
          labs(x = "RNA escape (CAST molecule fraction)",
               y = "Xi accessibility (CAST fragment fraction)",
               title = "RNA output against chromatin state in the SAME nucleus",
               subtitle = sprintf("nuclei passing both gates; Spearman rho = %s",
                                  paste(sprintf("%s %.2f", cors$sample, cors$rho), collapse = ", "))) +
          theme_minimal()
      )
      print(
        ggplot(bk, aes(rna_escape, xi_access, colour = label)) +
          geom_point(alpha = 0.5, size = 0.9) + facet_wrap(~sample) +
          labs(x = "RNA escape (CAST molecule fraction)", y = "Xi accessibility",
               colour = NULL, title = "The same, split by cell type") +
          theme_minimal()
      )
      dev.off()
      say("wrote rna_atac_pernucleus.pdf")
    } else {
      say("only %d nuclei pass both gates - too few for a per-nucleus correlation,", nrow(bk))
      say("  so the table is written but the figure is skipped.")
    }
  } else {
    say("no %s; skipping the per-nucleus RNA/ATAC join", rna_pn_f)
  }
}

say("")
say("output under %s", OUT)
say("A1 = %s. chrX CAST fraction here is ACCESSIBILITY of the inactive X,", A1_IS)
say("not escape transcription - see the header before quoting any number.")
