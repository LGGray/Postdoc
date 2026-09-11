# ---------------------------------------------------------------------------
# 05 - Consolidate Allelome.PRO2 output and plot the escape results.
#
# Reads whatever multiome_allelome_run.slurm managed to finish, from its stable
# SCRATCH path, so a run that hit the wall clock still gets plotted. Writes one
# consolidated TSV per tree to DSS and caches it, so re-plotting is instant.
#
# ALLELE ORIENTATION: A1 = C57BL/6 (B6), A2 = CAST/EiJ. Confirmed directly,
# which matters because the toolkit does not document it - score.R only shows
# allelic_ratio = A1_reads/total_reads, and A1/A2 arrive from the SNP bed's
# allele column through the pileup with no label attached.
#
# So escape = A2/(A1+A2), the CAST fraction: Xist is deleted on the B6 allele,
# B6 therefore cannot be silenced, and CAST is the inactive X in every nucleus.
#
# The pooled chrX skew is kept as an INDEPENDENT CHECK rather than as the
# source of the orientation. With the orientation known from outside the data
# that is a real test - a chrX A1 fraction that is NOT strongly B6 contradicts
# the whole chain and points at the Xic mask, the SNP file or the genotype.
# Deriving the orientation FROM the skew, as an earlier version did, was mildly
# circular: it assumed the biology in order to measure it.
#
# The loader mirrors load_allelome_tree() in
# OCM_heart/allelic_ratio/00_functions.R rather than sourcing it: that file
# calls dir.create() at load time and needs glmmTMB/broom.mixed, neither of
# which belongs in a plotting job.
#
#   conda activate seurat_env
#   Rscript multiome/05_allelome_plots.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr); library(tibble)
})

BASE <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2"
WORK <- file.path(BASE, "adult_aged_multiome")
OUT  <- file.path(WORK, "figures_allelome")
# ALLELOME_TREE overrides the scratch location, for a tree that has been
# copied to DSS ahead of a purge. Same Sys.getenv override convention as
# POSTDOC_ROOT in OCM_heart/allelic_ratio/. Note --export=NONE: set this
# INSIDE the job script, not in the submitting shell, or it will not arrive.
SCR  <- Sys.getenv("ALLELOME_TREE", file.path(Sys.getenv("SCRATCH"), "multiome_allelome"))
source(file.path(BASE, "Postdoc", "multiome", "00_helpers.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SAMPLES    <- c("9w", "78w")
ANNOT_BASE <- "chr_annotation_mm39.bed"
AUTOSOMES  <- paste0("chr", 1:19)
PRIOR_ESCAPE <- 0.127     # snRNA/spatial estimate this should reproduce
A1_IS      <- "B6"        # A1 = C57BL/6, A2 = CAST/EiJ - confirmed, not inferred
CHRX_A1_MIN <- 0.60       # below this, the data contradicts the known orientation
MIN_INFORMATIVE <- 20     # per-nucleus gate, applied at ANALYSIS time on the
                          # actual informative chrX count - not a UMI proxy

say <- function(...) cat(sprintf(...), "\n", sep = "")
# A missing tree is not fatal: the cached TSVs are enough to re-plot, which is
# the common case once scratch has been purged.
if (!dir.exists(SCR)) {
  message("no allelome tree at ", SCR,
          " - relying on cached TSVs; set ALLELOME_TREE to rescan")
}

# ---------------------------------------------------------------------------
# consolidate
# ---------------------------------------------------------------------------
# Allelome.PRO2 names each run <DATE>_<bambase>_<annot>_<minread>; bambase is
# the sinto group, i.e. a barcode (per nucleus) or a cell type (per cell type).
parse_group <- function(dirs) {
  x <- sub("^[0-9]{4}_[0-9]{2}_[0-9]{2}_", "", dirs)
  sub(paste0("_", ANNOT_BASE, "_[0-9]+$"), "", x)
}

load_tree <- function(tree) {
  cache <- file.path(OUT, sprintf("allelome_%s.tsv", tree))
  if (file.exists(cache) && !identical(Sys.getenv("REBUILD"), "1")) {
    say("using cached %s (REBUILD=1 to rescan)", basename(cache))
    return(read_tsv(cache, show_col_types = FALSE) %>% mutate(sample = as_sample(sample)))
  }
  rows <- list()
  for (id in SAMPLES) {
    paths <- list.files(file.path(SCR, tree, "out", id), pattern = "^locus_table\\.txt$",
                        recursive = TRUE, full.names = TRUE)
    say("%s / %s: %d locus tables", tree, id, length(paths))
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
  write_tsv(out, cache)
  say("wrote %s (%d rows)", basename(cache), nrow(out))
  mutate(out, sample = as_sample(sample))
}

ct <- load_tree("celltype")
pn <- load_tree("pernucleus")
if (is.null(ct) && is.null(pn)) stop("no Allelome.PRO2 output found under ", SCR)

# ---------------------------------------------------------------------------
# orientation
# ---------------------------------------------------------------------------
pool <- bind_rows(ct, pn)
px <- pool %>% filter(chr == "chrX") %>% summarise(a1 = sum(A1_reads), a2 = sum(A2_reads))
pa <- pool %>% filter(chr %in% AUTOSOMES) %>% summarise(a1 = sum(A1_reads), a2 = sum(A2_reads))
ar_x    <- px$a1 / (px$a1 + px$a2)
ar_auto <- pa$a1 / (pa$a1 + pa$a2)

say("")
# A1 = B6 is known, so escape is the CAST fraction. Not derived from the data.
escape_of <- function(a1, a2) a2 / (a1 + a2)

say("============ orientation check (A1 = %s, known) ============", A1_IS)
say("pooled chrX      A1(%s) fraction: %.4f  (%d / %d reads)", A1_IS, ar_x, px$a1, px$a1 + px$a2)
say("pooled autosomal A1(%s) fraction: %.4f  (%d / %d reads)", A1_IS, ar_auto, pa$a1, pa$a1 + pa$a2)
say("implied pooled chrX escape (CAST): %.4f", 1 - ar_x)
if (ar_x < CHRX_A1_MIN) {
  stop(sprintf(paste0("pooled chrX A1(B6) fraction is %.3f, below the %.2f floor.\n",
    "  A1 = B6 is known independently, and Xist is deleted on B6, so B6 cannot be\n",
    "  silenced and chrX MUST be strongly B6-skewed. This says it is not, which\n",
    "  contradicts the chain rather than revising the orientation. Candidates, in\n",
    "  rough order of likelihood: the Xic mask (is _no_Xist the right extent?), the\n",
    "  SNP file, the sample genotype, or the dedup/MAPQ filters. Not plotting a\n",
    "  number that cannot be interpreted."), ar_x, CHRX_A1_MIN))
}
say("-> consistent with A1 = B6 and full skewing; escape = A2/(A1+A2) = CAST fraction")
say("autosomal mapping bias: A1 excess of %+.4f from 0.5 (B6-reference bias)",
    ar_auto - 0.5)
say("===========================================================")
say("")

wilson <- function(k, n) {
  if (n == 0) return(c(NA, NA))
  p <- k / n; z <- 1.96; d <- 1 + z^2 / n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / d,
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / d)
}

# ---------------------------------------------------------------------------
# per cell type
# ---------------------------------------------------------------------------
if (!is.null(ct)) {
  ctx <- ct %>%
    mutate(region = case_when(chr == "chrX" ~ "chrX",
                              chr %in% AUTOSOMES ~ "autosomal", TRUE ~ NA_character_)) %>%
    filter(!is.na(region)) %>%
    group_by(sample, group, region) %>%
    summarise(A1 = sum(A1_reads), A2 = sum(A2_reads), .groups = "drop") %>%
    mutate(n = A1 + A2, escape = escape_of(A1, A2))
  ci <- t(mapply(function(a1, a2) {
    k <- if (A1_IS == "B6") a2 else a1
    wilson(k, a1 + a2)
  }, ctx$A1, ctx$A2))
  ctx$lo <- ci[, 1]; ctx$hi <- ci[, 2]
  ctx$label <- short_labels(ctx$group)

  say("per-celltype chrX escape (CAST fraction), Wilson 95%% CI:")
  print(as.data.frame(ctx %>% filter(region == "chrX") %>%
    select(sample, label, n, escape, lo, hi) %>% arrange(sample, desc(n))),
    row.names = FALSE, digits = 3)
  write_tsv(ctx, file.path(OUT, "celltype_escape.tsv"))

  dev_open(file.path(OUT, "celltype_escape.pdf"), width = 11, height = 6.5)
  print(
    ggplot(ctx %>% filter(region == "chrX"),
           aes(reorder(label, escape), escape, colour = sample)) +
      geom_hline(yintercept = PRIOR_ESCAPE, linetype = "dashed", colour = "grey40") +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2,
                    position = position_dodge(0.4)) +
      geom_point(size = 2.6, position = position_dodge(0.4)) +
      coord_flip() +
      sample_scale("colour") +
      labs(x = NULL, y = "chrX escape (CAST fraction)",
           title = "Per-cell-type chrX escape",
           subtitle = sprintf("dashed = %.1f%% from snRNA/spatial; bars are Wilson 95%% CI; n=1 animal per age",
                              100 * PRIOR_ESCAPE)) +
      theme_minimal()
  )
  print(
    ggplot(ctx, aes(reorder(label, escape), escape, colour = sample, shape = region)) +
      geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey60") +
      geom_point(size = 2.6, position = position_dodge(0.4)) + coord_flip() +
      sample_scale("colour") +
      labs(x = NULL, y = "escape / CAST fraction",
           title = "chrX against the autosomal control",
           subtitle = "autosomes should sit at 0.5; their offset is the B6-reference mapping bias") +
      theme_minimal()
  )
  dev.off()
  say("wrote celltype_escape.pdf")
}

# ---------------------------------------------------------------------------
# per nucleus
# ---------------------------------------------------------------------------
if (!is.null(pn)) {
  lab <- file.path(WORK, "figures_signac", "nucleus_celltype_joint.csv")
  meta <- if (file.exists(lab)) read_csv(lab, show_col_types = FALSE) else NULL

  pnx <- pn %>% filter(chr == "chrX") %>%
    transmute(sample, barcode = group, A1 = A1_reads, A2 = A2_reads,
              n_inf = A1_reads + A2_reads, escape = escape_of(A1_reads, A2_reads))
  if (!is.null(meta)) {
    pnx <- pnx %>% left_join(meta %>% select(sample, barcode, celltype_provisional),
                             by = c("sample", "barcode")) %>%
      mutate(label = short_labels(ifelse(is.na(celltype_provisional),
                                         "unlabelled", celltype_provisional)))
  } else pnx$label <- "all nuclei"

  say("per-nucleus chrX: %d nuclei, median informative molecules %.0f",
      nrow(pnx), median(pnx$n_inf))
  say("  nuclei at or above the %d-molecule gate: %d (%.1f%%)", MIN_INFORMATIVE,
      sum(pnx$n_inf >= MIN_INFORMATIVE), 100 * mean(pnx$n_inf >= MIN_INFORMATIVE))
  write_tsv(pnx, file.path(OUT, "pernucleus_escape.tsv"))

  keep <- pnx %>% filter(n_inf >= MIN_INFORMATIVE)
  SC <- celltype_scale(sort(unique(pnx$label)), "colour")

  dev_open(file.path(OUT, "pernucleus_escape.pdf"), width = 11, height = 7)
  print(
    ggplot(pnx, aes(n_inf, escape)) +
      geom_point(alpha = 0.25, size = 0.7) +
      geom_vline(xintercept = MIN_INFORMATIVE, linetype = "dashed", colour = OKABE_ITO[6]) +
      geom_hline(yintercept = PRIOR_ESCAPE, linetype = "dashed", colour = "grey40") +
      scale_x_log10() + facet_wrap(~sample) +
      labs(x = "informative chrX molecules per nucleus (log)", y = "escape (CAST fraction)",
           title = "Per-nucleus escape against the evidence behind it",
           subtitle = "the funnel IS the noise: low-count nuclei spread widely with no biology involved") +
      theme_minimal()
  )
  p <- ggplot(keep, aes(escape, colour = sample)) +
    geom_density(linewidth = 0.9) +
    geom_vline(xintercept = PRIOR_ESCAPE, linetype = "dashed", colour = "grey40") +
    sample_scale("colour") +
    labs(x = "escape (CAST fraction)", y = "density",
         title = sprintf("Per-nucleus escape distribution (>= %d informative molecules)", MIN_INFORMATIVE),
         subtitle = "n=1 animal per age - the age contrast is descriptive only") +
    theme_minimal()
  print(p)
  print(p + facet_wrap(~label, scales = "free_y") +
          labs(title = "Per-nucleus escape by cell type"))
  print(
    ggplot(keep, aes(reorder(label, escape, median), escape, fill = sample)) +
      geom_boxplot(outlier.size = 0.4, position = position_dodge(0.8)) +
      geom_hline(yintercept = PRIOR_ESCAPE, linetype = "dashed", colour = "grey40") +
      coord_flip() +
      sample_scale("fill") +
      labs(x = NULL, y = "escape (CAST fraction)",
           title = "Per-nucleus escape by cell type and age") + theme_minimal()
  )
  # autosomal control at nucleus level: the same plot must centre on 0.5
  pna <- pn %>% filter(chr %in% AUTOSOMES) %>%
    group_by(sample, group) %>%
    summarise(A1 = sum(A1_reads), A2 = sum(A2_reads), .groups = "drop") %>%
    mutate(n_inf = A1 + A2, ratio = escape_of(A1, A2)) %>%
    filter(n_inf >= MIN_INFORMATIVE)
  print(
    ggplot(pna, aes(ratio, colour = sample)) + geom_density(linewidth = 0.9) +
      geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey60") +
      sample_scale("colour") +
      labs(x = sprintf("autosomal %s fraction", if (A1_IS == "B6") "CAST" else "B6"),
           y = "density", title = "Autosomal control, per nucleus",
           subtitle = "should centre on 0.5; any offset is the B6-reference mapping bias, not biology") +
      theme_minimal()
  )
  dev.off()
  say("wrote pernucleus_escape.pdf")
}

say("")
say("output under %s", OUT)
say("A1 = %s (known); escape reported as the CAST fraction", A1_IS)
