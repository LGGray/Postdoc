# ---------------------------------------------------------------------------
# Shared libraries and helper functions for the allelic ratio analysis.
# Split out of Allelic_ratio_testing.R, which is kept intact as the original.
# Sourced by 01-05; not meant to be run on its own.
# ---------------------------------------------------------------------------

library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(broom.mixed)


run_dispersion_lrt <- function(df, ref_level, comparison_level) {
  df$sample <- factor(df$sample, levels = c(ref_level, comparison_level))

  fit_const_disp <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                            family = betabinomial(link = "logit"), data = df)
  fit_var_disp <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                          dispformula = ~ sample,
                          family = betabinomial(link = "logit"), data = df)

  lrt_row <- data.frame(anova(fit_const_disp, fit_var_disp))[2, , drop = FALSE]
  disp_coef_name <- paste0("sample", comparison_level)
  disp_beta <- unname(fixef(fit_var_disp)$disp[disp_coef_name])

  direction <- if (is.na(disp_beta)) {
    NA_character_
  } else if (disp_beta < 0) {
    paste(comparison_level, "more heterogeneous")
  } else {
    paste(comparison_level, "less heterogeneous")
  }

  lrt_row$disp_beta <- disp_beta
  lrt_row$direction <- direction
  lrt_row
}

# Simulation-based false-positive rate for run_dispersion_lrt, quantifying the
# n=1-animal-per-condition caveat: if the two conditions truly share the SAME
# dispersion, and each condition is still just a single animal, how often does
# the LRT call a "significant" difference purely from animal-to-animal noise?
#
# Mechanics: simulate two animals from IDENTICAL true (mu, phi), but let each
# animal's *realized* phi wobble around that shared truth on the log scale
# (animal_disp_sd) -- this represents ordinary biological variability between
# individuals, which a single-animal-per-condition design can never separate
# from a genuine condition effect. Cells are simulated at the real per-cell
# read depths and real group sizes so the simulation matches the actual data
# structure. animal_disp_sd cannot be estimated from n=1 data -- it has to be
# assumed, so report a range of plausible values rather than trusting one.
simulate_dispersion_null_fpr <- function(df, ref_level, comparison_level,
                                         animal_disp_sd = 0.3, n_reps = 1000,
                                         seed = 1) {
  set.seed(seed)
  df$sample <- factor(df$sample, levels = c(ref_level, comparison_level))

  # anchor simulation parameters to this celltype/comparison's real data:
  # real per-group mean allelic proportion and the real shared dispersion
  # under the null (single phi fit across both groups)
  fit_null <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                      family = betabinomial(link = "logit"), data = df)
  mu_hat <- predict(fit_null,
                     newdata = data.frame(sample = factor(c(ref_level, comparison_level),
                                                           levels = c(ref_level, comparison_level))),
                     type = "response")
  mu_ref  <- unname(mu_hat[1])
  mu_comp <- unname(mu_hat[2])
  phi_true <- exp(unname(fixef(fit_null)$disp[["(Intercept)"]]))

  ref_reads  <- df$total_reads[df$sample == ref_level]
  comp_reads <- df$total_reads[df$sample == comparison_level]

  p_values <- rep(NA_real_, n_reps)
  n_failed <- 0

  for (i in seq_len(n_reps)) {
    # each simulated animal gets its own realized dispersion, both drawn from
    # the same distribution since the true dispersion is identical by design
    phi_ref  <- phi_true * exp(rnorm(1, 0, animal_disp_sd))
    phi_comp <- phi_true * exp(rnorm(1, 0, animal_disp_sd))

    p_ref  <- rbeta(length(ref_reads),  mu_ref  * phi_ref,  (1 - mu_ref)  * phi_ref)
    p_comp <- rbeta(length(comp_reads), mu_comp * phi_comp, (1 - mu_comp) * phi_comp)

    A1_ref  <- rbinom(length(ref_reads),  ref_reads,  p_ref)
    A1_comp <- rbinom(length(comp_reads), comp_reads, p_comp)

    df_sim <- data.frame(
      A1_reads = c(A1_ref, A1_comp),
      A2_reads = c(ref_reads - A1_ref, comp_reads - A1_comp),
      sample = c(rep(ref_level, length(ref_reads)), rep(comparison_level, length(comp_reads)))
    )

    lrt_row <- tryCatch(run_dispersion_lrt(df_sim, ref_level, comparison_level),
                         error = function(e) NULL, warning = function(w) NULL)

    if (is.null(lrt_row)) {
      n_failed <- n_failed + 1
    } else {
      p_values[i] <- lrt_row$Pr..Chisq.
    }
  }

  list(
    fpr = mean(p_values < 0.05, na.rm = TRUE),
    n_failed = n_failed,
    n_reps = n_reps,
    p_values = p_values
  )
}

# Logistic regression LRT testing whether the proportion of monoallelic-like
# cells (allelic_ratio <= 0.1 or >= 0.9) differs between two sample levels.
run_monoallelic_lrt <- function(df, ref_level, comparison_level) {
  df$sample <- factor(df$sample, levels = c(ref_level, comparison_level))

  fit_null <- glm(monoallelic ~ 1, family = binomial, data = df)
  fit_full <- glm(monoallelic ~ sample, family = binomial, data = df)

  lrt <- anova(fit_null, fit_full, test = "LRT")

  coef_name <- paste0("sample", comparison_level)
  beta <- unname(coef(fit_full)[coef_name])

  direction <- if (is.na(beta)) {
    NA_character_
  } else if (beta > 0) {
    paste(comparison_level, "more monoallelic")
  } else {
    paste(comparison_level, "less monoallelic")
  }

  data.frame(
    Df = lrt$Df[2],
    Deviance = lrt$Deviance[2],
    p_value = lrt[["Pr(>Chi)"]][2],
    beta = beta,
    direction = direction
  )
}

fdr_to_stars <- function(fdr) {
  ifelse(is.na(fdr), "ns",
         ifelse(fdr <= 0.001, "***",
                ifelse(fdr <= 0.01, "**",
        ifelse(fdr <= 0.05, "*", "ns"))))
}

# ---------------------------------------------------------------------------
# Total-read cutoff and the output directory keyed to it.
#
# The cutoff is not a neutral QC knob: 06_cutoff_sweep.R shows mean AR rising
# from 0.85 at 10-25 chrX reads to 0.97 above 200, a real depth bias rather
# than sampling noise, so results are only comparable within one cutoff.
# Hence one directory per cutoff -- re-running at a different value writes
# somewhere new instead of silently overwriting the previous set.
#
# Defined here rather than in 02 so that 03 and 05, which read 02's per-cell
# metadata handoff, cannot drift onto a different cutoff than the file they
# are reading was built with.
#
# Override for a one-off:  MIN_TOTAL_READS=50 Rscript allelic_ratio/02_whole_chrX.R
MIN_TOTAL_READS <- as.integer(Sys.getenv("MIN_TOTAL_READS", "30"))
stopifnot(!is.na(MIN_TOTAL_READS), MIN_TOTAL_READS >= 1)

# Which Allelome.PRO2 tree the results were built from. Everything below hangs
# off this, so pointing 02-08 at the deduplicated tree is one variable, not a
# fork of 982 lines:
#
#   RESULTS_ROOT=Allelic_ratio_results_dedup Rscript allelic_ratio/02_whole_chrX.R
#
# The two roots must not be mixed. total_reads means reads in the default tree
# and molecules in the deduplicated one, so a cutoff of 30 is a different filter
# in each -- see 09_dedup_comparison.R for how much different.
RESULTS_ROOT <- Sys.getenv("RESULTS_ROOT", "Allelic_ratio_results")
dir.create(RESULTS_ROOT, showWarnings = FALSE, recursive = TRUE)

CUTOFF_DIR <- file.path(RESULTS_ROOT,
                        paste0("cutoff_", MIN_TOTAL_READS))
dir.create(CUTOFF_DIR, showWarnings = FALSE, recursive = TRUE)

# The raw per-cell allelic ratio table is cutoff-INDEPENDENT: it is every cell
# with any chrX coverage, and both 02 and 06 filter it themselves. It stays at
# the top level rather than being duplicated into each cutoff directory.
ALLELIC_RATIOS_FILE <- file.path(RESULTS_ROOT, "whole_chr_allelic_ratios.txt")

# The core escape block (04) needs its OWN cutoff, not MIN_TOTAL_READS. Its
# total_reads counts reads over roughly twenty escape genes, not the whole X,
# so the same number means something very different: the whole-chrX median is
# ~50 reads while the block's is a handful. Applying 30 here would discard
# almost every cell. Hence a separate constant and a separate directory --
# the two cutoffs are not comparable and should not share a folder name.
MIN_CEB_READS <- as.integer(Sys.getenv("MIN_CEB_READS", "5"))
stopifnot(!is.na(MIN_CEB_READS), MIN_CEB_READS >= 1)

# The per-gene block (03) has two thresholds of its own, and neither is
# MIN_TOTAL_READS: that one is a whole-chrX depth per cell, while these are
# per-GENE depths, which are one to two orders of magnitude smaller.
#
# Both were hardcoded as a bare `10` in three places in 03_all_genes.R while a
# name that was never defined anywhere, MIN_READS, appeared in its plot titles
# and comments. That name reached a runtime expression at 03_all_genes.R:325,
# so the block threw as soon as control got that far; every earlier reference to
# it was a comment and so silently described a number that was not there. Named
# here, used there, so the figures now state the threshold they were actually
# built with.
#
#   MIN_GENE_READS  pooled over the cells of one (gene, celltype), across
#                   samples. Selects which genes have any signal at all.
#   MIN_CELL_READS  within a single cell, for the per-cell scatter and violin
#                   panels, where one cell's ratio has to mean something on its
#                   own.
# The core escape gene set, in one place.
#
# 03_all_genes.R flagged these by reading core_escape_genes_gene_df.txt, a file
# that is READ there and WRITTEN nowhere in this repo - so
# is_core_escape was silently all FALSE and the escape overlay in
# VCM_subcluster_per_gene_delta_scatter.pdf was an empty layer with no error
# anywhere. The list itself was never missing: it is the same eleven genes
# OCM_heart/core_escape_SNPs.R uses to build the core escape SNP bed. Naming it
# here removes the phantom file from the loop.
#
# Keep in step with new_core_escape_genes in OCM_heart/core_escape_SNPs.R, which
# is a standalone script on the reference tree and does not source this file.
CORE_ESCAPE_GENES <- c("Kdm5c", "Kdm6a", "Ddx3x", "Eif2s3x", "Ftx",
                       "5530601H04Rik", "Jpx", "Pbdc1", "Utp14a",
                       "Akap17a", "Sts")

MIN_GENE_READS <- as.integer(Sys.getenv("MIN_GENE_READS", "10"))
MIN_CELL_READS <- as.integer(Sys.getenv("MIN_CELL_READS", "10"))
stopifnot(!is.na(MIN_GENE_READS), MIN_GENE_READS >= 1,
          !is.na(MIN_CELL_READS), MIN_CELL_READS >= 1)

CEB_DIR <- file.path(RESULTS_ROOT,
                     paste0("core_escape_cutoff_", MIN_CEB_READS))
dir.create(CEB_DIR, showWarnings = FALSE, recursive = TRUE)

# As with the whole-chrX table above, the core escape block's raw per-cell
# ratios are cutoff-independent -- 04 writes them before applying its filter,
# and 07 sweeps cutoffs over them. Top level, not inside a cutoff directory.
CEB_RATIOS_FILE <- file.path(RESULTS_ROOT, "core_escape_block_new_allelic_ratio_table.txt")

# ---------------------------------------------------------------------------
# Reading an Allelome.PRO2 output tree.
#
# Allelome.PRO2 names each run <date>_<bam>_<annotation>_<run no.>, so the cell
# barcode has to be recovered from the directory name. Do NOT do this by
# splitting the full path on "_" and taking a fixed index, the way the original
# Allelic_ratio_testing.R does -- that index depends on how many underscores are
# in the tree's own directory name, so it is 4/5 under Allelome.PRO2/ and 6/7
# under Allelome.PRO2_all_genes/, and silently wrong under any new tree name.
# ---------------------------------------------------------------------------
parse_allelome_cell <- function(paths) {
  d <- basename(dirname(paths))
  d <- sub("^[0-9]{4}_[0-9]{2}_[0-9]{2}_", "", d)   # drop the date stamp
  sub("_[^_]*\\.bed_[0-9]+$", "", d)                # drop _<annotation>_<run no.>
}

# The cell key the Seurat object uses, <sample>_<barcode>, from a locus table
# path. Trees are laid out as <tree>/<sample>/<job dir>/locus_table.txt.
#
# sinto names each cell's BAM from its cells file, and this project has written
# that file both ways: some trees carry <barcode>.bam and some <sample>_<barcode>.bam.
# Prepending the sample unconditionally is therefore wrong half the time, and
# wrong in the worst way - it produces "9w_9w_AAACCC-1", which is a valid string
# that matches no cell, so the subset comes back empty instead of erroring.
# Testing for the prefix is correct under both conventions.
allelome_cell_barcode <- function(paths, sample_of = NULL) {
  cell <- parse_allelome_cell(paths)
  smp  <- if (is.null(sample_of)) basename(dirname(dirname(paths))) else sample_of
  ifelse(startsWith(cell, paste0(smp, "_")), cell, paste0(smp, "_", cell))
}

# Refuse to carry on with barcodes that do not match the object.
#
# This is the failure mode the whole helper above exists to prevent: a wrong
# barcode parse yields strings Seurat simply does not find, so subset() returns
# few or no cells and every downstream number is computed on that remainder
# without a single error being raised. 04_core_escape.R - the LOX analysis - had
# a hardcoded field index that only held for one tree name.
check_barcode_match <- function(barcodes, obj_cells, what = "cells",
                                min_frac = 0.5) {
  n <- length(unique(barcodes))
  hit <- length(intersect(unique(barcodes), obj_cells))
  frac <- if (n) hit / n else 0
  message(sprintf("%s: %d of %d parsed barcodes match the object (%.1f%%)",
                  what, hit, n, 100 * frac))
  if (hit == 0) {
    stop(what, ": NONE of the ", n, " parsed barcodes match the Seurat object.\n",
         "  parsed:  ", paste(head(unique(barcodes), 3), collapse = ", "), "\n",
         "  object:  ", paste(head(obj_cells, 3), collapse = ", "), "\n",
         "  The barcode parse is wrong - see allelome_cell_barcode(). Nothing ",
         "downstream would\n  have errored, it would just have been computed on ",
         "an empty subset.")
  }
  if (frac < min_frac) {
    warning(what, ": only ", round(100 * frac, 1), "% of parsed barcodes match ",
            "the object. Expected most of them; check the tree and the object ",
            "are the same experiment.")
  }
  invisible(frac)
}

# One row per (cell, chromosome), counts summed over whatever intervals the
# annotation bed defined.
#
# Summing is correct whether the bed carries one interval per chromosome or one
# per gene, but it does NOT make two trees built from different beds comparable:
# a per-gene bed counts only reads inside genes, a whole-chromosome interval
# counts everything on the chromosome. Compare trees only where the annotation
# matches.
#
# allelic_ratio is deliberately not taken from the locus table. Two conventions
# are in play in this project (A1/total, and max(A1,A2)/total) and averaging a
# per-locus ratio column is not the same as a depth-weighted whole-chromosome
# ratio anyway. Both forms are returned so callers state which they mean:
#   ar_a1  = A1 / total, directional, 0-1
#   ar_dom = max(A1, A2) / total, dominant-allele fraction, 0.5-1
load_allelome_tree <- function(tree, samples = c("9w", "78w", "Sham", "TAC")) {
  paths <- unlist(lapply(samples, function(s) {
    list.files(file.path(tree, s), pattern = "locus_table.txt",
               recursive = TRUE, full.names = TRUE)
  }))
  if (!length(paths)) stop("No locus tables under ", tree)

  sample_of <- basename(dirname(dirname(paths)))
  cell_of   <- allelome_cell_barcode(paths, sample_of)

  rows <- lapply(seq_along(paths), function(i) {
    tmp <- tryCatch(read.delim(paths[i], header = TRUE, stringsAsFactors = FALSE),
                    error = function(e) NULL)
    if (is.null(tmp) || !nrow(tmp)) return(NULL)
    if (!all(c("chr", "A1_reads", "A2_reads") %in% names(tmp))) return(NULL)
    tmp$chr <- as.character(tmp$chr)
    agg <- aggregate(cbind(A1_reads, A2_reads) ~ chr, data = tmp, FUN = sum)
    # cell_of is already <sample>_<barcode>; allelome_cell_barcode() added the
    # prefix only where the BAM name lacked it, so nothing is prepended twice.
    agg$cell_barcode <- cell_of[i]
    agg$sample <- sample_of[i]
    agg
  })

  out <- bind_rows(rows)
  if (!nrow(out)) stop("No readable locus tables under ", tree)
  out$total_reads <- out$A1_reads + out$A2_reads
  out <- out[out$total_reads > 0, ]
  out$ar_a1  <- out$A1_reads / out$total_reads
  out$ar_dom <- pmax(out$A1_reads, out$A2_reads) / out$total_reads
  out[, c("cell_barcode", "sample", "chr",
          "A1_reads", "A2_reads", "total_reads", "ar_a1", "ar_dom")]
}

# The autosomal control. Reads over all autosomes in one cell, which carry the
# same B6-reference mapping bias against CAST reads, the same library skew and
# the same duplication structure as chrX, but no monoallelic biology. It is the
# empirical null: chrX is only interesting where it exceeds this.
AUTOSOMES <- paste0("chr", 1:19)

collapse_autosomes <- function(df, autosomes = AUTOSOMES) {
  a <- df[df$chr %in% autosomes, ]
  if (!nrow(a)) return(a[0, ])
  out <- a %>%
    group_by(cell_barcode, sample) %>%
    summarise(A1_reads = sum(A1_reads), A2_reads = sum(A2_reads),
              n_chr = n(), .groups = "drop") %>%
    as.data.frame()
  out$chr <- "autosomal"
  out$total_reads <- out$A1_reads + out$A2_reads
  out <- out[out$total_reads > 0, ]
  out$ar_a1  <- out$A1_reads / out$total_reads
  out$ar_dom <- pmax(out$A1_reads, out$A2_reads) / out$total_reads
  out
}

# Whole-chromosome code assigns one value per cell with match(), which silently
# takes the FIRST matching row. That is correct only when the annotation bed had
# one interval per chromosome (chr_annotation_mm39.bed, 21 intervals). Fed a
# per-gene bed (annotation_us_mm39_chrX.bed, 3227 chrX intervals) the same code
# hands every cell one arbitrary gene while calling it the whole X. 9w was
# scored that way in the original tree, so this is a real failure mode, not a
# hypothetical. 10_build_ratio_table.R sums to one row per (cell, chr); this
# asserts whatever is actually being read satisfies that.
assert_one_row_per_cell <- function(df, what = "allelic ratio table") {
  dup <- df$cell_barcode[duplicated(df$cell_barcode)]
  if (length(dup)) {
    n_per <- sort(table(df$cell_barcode[df$cell_barcode %in% dup]), decreasing = TRUE)
    stop(sprintf(
      paste0("%s has %d cells with more than one row (worst: %s, %d rows).\n",
             "  Whole-chromosome code takes the first row per cell, so these ",
             "cells would silently get one interval, not the chromosome.\n",
             "  Rebuild with 10_build_ratio_table.R, which sums to one row per ",
             "(cell, chr)."),
      what, length(unique(dup)), names(n_per)[1], as.integer(n_per)[1]))
  }
  invisible(df)
}

# Which annotation bed each sample in a tree was scored against. Allelome.PRO2
# puts the bed basename in every run directory name, so this is recoverable
# after the fact - and it has to be, because comparing two trees is only valid
# where the annotation matches. Returns one row per (sample, annotation), so a
# sample scored against two beds shows up as two rows.
allelome_annotations <- function(tree, samples = c("9w", "78w", "Sham", "TAC")) {
  rows <- lapply(samples, function(s) {
    d <- list.dirs(file.path(tree, s), recursive = FALSE, full.names = FALSE)
    a <- unique(sub("^.*_([^_]*\\.bed)_[0-9]+$", "\\1", d[grepl("\\.bed_[0-9]+$", d)]))
    if (!length(a)) return(NULL)
    data.frame(sample = s, annotation = a, stringsAsFactors = FALSE)
  })
  bind_rows(rows)
}
