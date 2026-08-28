# Per-allele rates under every SNP mask, side by side (task 2 of
# spatial/NEXT_ANALYSIS.md).
#
# The question: how much of the CAST rate difference between the two sections
# (4.79 -> 6.94 informative chrX CAST UMIs per 1000 informative autosomal UMIs,
# FC 1.449) survives masking the Xic properly rather than just the Xist gene
# span? Ftx, Jpx and Tsix sit outside the annotated Xist interval and all
# intergenic Xic reads survive it, and in this design those are effectively
# CAST-only - so a "no-Xist" mask can leave the entire effect behind.
#
# Nothing is re-counted here and nothing is re-run. ase_bin_allele_counts.py
# writes a column pair per --subset-bed in ONE pass over the BAM, so the masked
# and unmasked chrX numbers come from the same UMI calls and are exactly
# comparable. That is the whole reason for doing it as subsets rather than as a
# second counting run against a second SNP bed: two runs would differ by
# whatever else changed between them.
#
# BOTH MASKS ARE ALWAYS REPORTED. Never switch silently - the point of the
# table is the difference between them.
#
# Reported per sample and per set:
#   rate_ref / rate_alt   informative chrX UMIs of ONE allele per 1000
#                         informative autosomal UMIs. Reported separately,
#                         never as the ratio alone: a ratio conflates a
#                         numerator change with a denominator change, which is
#                         exactly what went wrong the first time round.
#   fold change between samples, per allele, with a tile-level bootstrap CI.
#
# The bootstrap resamples TILES, not UMIs. UMIs within a tile are not
# independent - one amplified molecule, one cell's burst - so a UMI-level
# interval would be far too narrow.
#
# n = 1 animal per age. Every fold change here describes two sections.
#
# Run:  sbatch ~/Postdoc/slurm/spatial_mask_comparison.slurm
#  or:  Rscript ~/Postdoc/spatial/ase_mask_comparison.R

.libPaths(c("~/R/matrix-dev", .libPaths()))

msg <- function(...) message(sprintf(...))
msg_table <- function(x) message(paste(capture.output(print(x)), collapse = "\n"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

##### ---------------------- CONFIG ---------------------- #####
BASE    <- Sys.getenv("BASE",
  "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial")
SAMPLES <- strsplit(Sys.getenv("SAMPLES", "9w,78w"), ",")[[1]]
TILE_UM <- as.integer(Sys.getenv("TILE_UM", "64"))
OUT_DIR <- Sys.getenv("OUT_DIR", file.path(BASE, "ase"))
N_BOOT  <- as.integer(Sys.getenv("N_BOOT", "10000"))

# Tiles below this many informative chrX UMIs still count toward the rates -
# the rate is a pooled quantity and dropping shallow tiles would bias it - but
# the per-tile ratio summaries at the bottom need a floor to mean anything.
MIN_TILE_UMI <- 20L
##### ------------------------------------------------------ #####

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
k_bins <- TILE_UM %/% 2L

fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), ...)
  else fread(path, ...)
}

load_tiles <- function(smp) {
  f <- file.path(BASE, "ase", smp, "bin_allele_counts.tsv")
  if (!file.exists(f)) { msg("  no %s - skipping %s", f, smp); return(NULL) }
  b <- fread_any(f)
  setnames(b, tolower(names(b)))
  b[, `:=`(tile_row = array_row %/% k_bins, tile_col = array_col %/% k_bins)]
  cnt <- grep("_(ref|alt)$", names(b), value = TRUE)
  t <- b[, lapply(.SD, sum), by = .(tile_row, tile_col), .SDcols = cnt]
  t[, `:=`(a_n = a_ref + a_alt, x_n = x_ref + x_alt, sample = smp)]
  msg("  %s: %d tiles at %d um, %d chrX and %d autosomal informative UMIs",
      smp, nrow(t), TILE_UM, t[, sum(x_n)], t[, sum(a_n)])
  t
}

msg("Tile size %d um, samples: %s, %d bootstrap draws",
    TILE_UM, paste(SAMPLES, collapse = ", "), N_BOOT)
tl <- rbindlist(lapply(SAMPLES, load_tiles), fill = TRUE)
if (!nrow(tl)) stop("No sample had a count table")
tl[, sample := factor(sample, levels = intersect(SAMPLES, unique(sample)))]

SETS <- sub("_ref$", "", grep("_ref$", names(tl), value = TRUE))
SETS <- setdiff(SETS, "a")
if (!length(SETS)) stop("No <name>_ref columns in the count tables")
msg("Sets: %s", paste(SETS, collapse = ", "))
if (!"noxic" %in% SETS) {
  msg("NOTE: no 'noxic' columns, so there is no Xic-masked comparison to make.")
  msg("      Re-run the counting pass with --subset-bed noxic=<Xic bed>:complement")
  msg("      (slurm/spatial_ase_sweep.slurm already does this once the interval")
  msg("      beds exist - build them with OCM_heart/core_escape_SNPs.R).")
}

# ------------------------------------------------------------- bootstrap
#
# One tile resample per replicate, shared across all sets and both alleles, so
# the sets are compared on the same draws and their differences carry the
# correlation between them rather than being treated as independent.
boot_rates <- function(d, sets, n_boot) {
  n <- nrow(d)
  a <- d$a_n
  num <- lapply(sets, function(s) list(ref = d[[paste0(s, "_ref")]],
                                       alt = d[[paste0(s, "_alt")]]))
  names(num) <- sets
  out <- array(NA_real_, dim = c(n_boot, length(sets), 2),
               dimnames = list(NULL, sets, c("ref", "alt")))
  for (b in seq_len(n_boot)) {
    w <- tabulate(sample.int(n, n, replace = TRUE), nbins = n)
    den <- sum(w * a)
    if (den == 0) next
    for (i in seq_along(sets)) {
      out[b, i, 1] <- 1000 * sum(w * num[[i]]$ref) / den
      out[b, i, 2] <- 1000 * sum(w * num[[i]]$alt) / den
    }
  }
  out
}

set.seed(1)
boots <- list()
point <- rbindlist(lapply(levels(tl$sample), function(smp) {
  d <- tl[sample == smp]
  boots[[smp]] <<- boot_rates(d, SETS, N_BOOT)
  rbindlist(lapply(SETS, function(s) {
    r <- sum(d[[paste0(s, "_ref")]]); al <- sum(d[[paste0(s, "_alt")]])
    data.table(sample = smp, set = s, tiles = nrow(d), auto_umi = sum(d$a_n),
               ref_umi = r, alt_umi = al,
               rate_ref = 1000 * r / sum(d$a_n),
               rate_alt = 1000 * al / sum(d$a_n),
               ratio = if (r + al > 0) r / (r + al) else NA_real_)
  }))
}))

ci <- rbindlist(lapply(names(boots), function(smp) {
  bb <- boots[[smp]]
  rbindlist(lapply(SETS, function(s) rbindlist(lapply(c("ref", "alt"), function(al) {
    v <- bb[, s, al]
    v <- v[is.finite(v)]
    data.table(sample = smp, set = s, allele = al,
               lo = as.numeric(quantile(v, 0.025)),
               hi = as.numeric(quantile(v, 0.975)))
  }))))
}))

rates <- melt(point, id.vars = c("sample", "set", "tiles", "auto_umi", "ratio",
                                 "ref_umi", "alt_umi"),
              measure.vars = c("rate_ref", "rate_alt"),
              variable.name = "allele", value.name = "rate")
rates[, allele := sub("rate_", "", as.character(allele))]
rates <- merge(rates, ci, by = c("sample", "set", "allele"))
rates[, `:=`(strain = fifelse(allele == "ref", "Bl6 (A1)", "CAST (A2)"),
             umi = fifelse(allele == "ref", ref_umi, alt_umi),
             sample = factor(sample, levels = levels(tl$sample)))]

msg("\n--- per-allele rates: chrX UMIs of one allele per 1000 autosomal UMIs ---")
msg_table(rates[order(sample, set, allele),
                .(sample, set, strain, umi, rate = round(rate, 3),
                  ci = sprintf("[%.3f, %.3f]", lo, hi))])

# ------------------------------------------------------- fold changes
#
# Between the two samples, per set and per allele, with a bootstrap CI on the
# ratio itself. The two samples are independent sections, so their replicates
# are paired arbitrarily - which is what an independent bootstrap of a ratio of
# two independent estimates is.
fc <- NULL
if (length(names(boots)) >= 2L) {
  lv <- names(boots)
  a1 <- lv[1]; a2 <- lv[length(lv)]
  fc <- rbindlist(lapply(SETS, function(s) rbindlist(lapply(c("ref", "alt"), function(al) {
    v1 <- boots[[a1]][, s, al]; v2 <- boots[[a2]][, s, al]
    ok <- is.finite(v1) & is.finite(v2) & v1 > 0
    r <- v2[ok] / v1[ok]
    p1 <- point[sample == a1 & set == s][[paste0("rate_", al)]]
    p2 <- point[sample == a2 & set == s][[paste0("rate_", al)]]
    data.table(set = s, allele = al,
               strain = if (al == "ref") "Bl6 (A1)" else "CAST (A2)",
               rate_young = p1, rate_old = p2,
               fc = if (p1 > 0) p2 / p1 else NA_real_,
               lo = as.numeric(quantile(r, 0.025)),
               hi = as.numeric(quantile(r, 0.975)))
  }))))
  msg("\n--- %s -> %s fold change in the per-allele rate ---", a1, a2)
  msg("n = 1 animal per age. This is a comparison of two sections.")
  msg_table(fc[order(allele, set),
               .(set, strain, young = round(rate_young, 3),
                 old = round(rate_old, 3), FC = round(fc, 3),
                 ci = sprintf("[%.3f, %.3f]", lo, hi),
                 sig = fifelse(lo > 1 | hi < 1, "*", ""))])

  # The acceptance criterion, stated.
  if (all(c("x", "noxic") %in% SETS)) {
    for (al in c("alt", "ref")) {
      u <- fc[set == "x" & allele == al]
      m <- fc[set == "noxic" & allele == al]
      xic <- fc[set == "xic" & allele == al]
      msg("\n=== ACCEPTANCE CRITERION (task 2), %s ===", u$strain)
      msg("Xist-span mask only : %.3f -> %.3f, FC %.3f [%.3f, %.3f]",
          u$rate_young, u$rate_old, u$fc, u$lo, u$hi)
      msg("Xist +- 500 kb mask : %.3f -> %.3f, FC %.3f [%.3f, %.3f]",
          m$rate_young, m$rate_old, m$fc, m$lo, m$hi)
      if (nrow(xic)) {
        msg("the masked region alone: %.3f -> %.3f, FC %.3f [%.3f, %.3f]",
            xic$rate_young, xic$rate_old, xic$fc, xic$lo, xic$hi)
        msg("  (the Xic is %.1f%% / %.1f%% of the unmasked rate in the two samples)",
            100 * xic$rate_young / u$rate_young,
            100 * xic$rate_old / u$rate_old)
      }
      if (is.finite(u$fc) && is.finite(m$fc)) {
        # "How much survives" is only a meaningful question if there was an
        # effect to survive. With an unmasked fold change at 1 the denominator
        # is zero and the share is arbitrarily large in either direction, so it
        # is not reported rather than reported as a big number.
        surv <- if (abs(u$fc - 1) < 0.02) NA_real_ else (m$fc - 1) / (u$fc - 1)
        if (!is.finite(surv)) {
          msg("-> The unmasked fold change is %.3f, too close to 1 to apportion.",
              u$fc)
        } else {
          msg("-> %.0f%% of the effect survives Xic masking.", 100 * surv)
          if (surv < 0.5) {
            msg("   Most of it was the Xic. The unmasked number cannot be reported;")
            msg("   rebuild the tile map on the no_Xic500kb SNP bed before using it")
            msg("   (slurm/spatial_sinto_tiles.slurm 64 0 no_Xic500kb).")
          } else {
            msg("   The effect is not the Xic. Report the masked number as primary")
            msg("   and the unmasked one alongside it.")
          }
        }
      }
    }
  }
}

# --------------------------------------------- per-tile ratio summaries
#
# The tile allelic ratio, for the record and for comparison with the published
# map. NOT an age comparison: it is dominated by depth and by library noise, and
# the two sections are not depth-matched (task 6). Reported per mask because the
# mask is the only thing changing here.
tr <- rbindlist(lapply(SETS, function(s) {
  d <- copy(tl)
  d[, `:=`(n = get(paste0(s, "_ref")) + get(paste0(s, "_alt")),
           r = get(paste0(s, "_ref")))]
  d <- d[n >= MIN_TILE_UMI]
  if (!nrow(d)) return(NULL)
  d[, .(set = s, tiles = .N, median_n = as.integer(median(n)),
        mean_ratio = round(mean(r / n), 4),
        pooled_ratio = round(sum(r) / sum(n), 4),
        frac_gt_0.9 = round(mean(r / n > 0.9), 4),
        frac_lt_0.1 = round(mean(r / n < 0.1), 4)), by = sample]
}))
msg("\n--- per-tile chrX allelic ratio, tiles with >= %d UMIs in that set ---",
    MIN_TILE_UMI)
msg("Depth is not matched between samples, so do NOT read this as an age effect.")
msg_table(tr[order(sample, set)])

fwrite(rates, file.path(OUT_DIR, sprintf("mask_comparison_rates_%dum.csv", TILE_UM)))
if (!is.null(fc)) {
  fwrite(fc, file.path(OUT_DIR,
                       sprintf("mask_comparison_fold_change_%dum.csv", TILE_UM)))
}
fwrite(tr, file.path(OUT_DIR, sprintf("mask_comparison_tile_ratio_%dum.csv", TILE_UM)))

##### --------------------------- figures --------------------------- #####

pdf(file.path(OUT_DIR, sprintf("mask_comparison_%dum.pdf", TILE_UM)),
    width = 8.5, height = 5.5)
th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8, colour = "#52514e"),
        plot.caption = element_text(size = 7, colour = "#52514e", hjust = 0))

print(
  ggplot(rates, aes(set, rate, fill = sample)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_errorbar(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.8), width = 0.2,
                  linewidth = 0.3) +
    facet_wrap(~ strain, scales = "free_y") +
    scale_fill_manual(values = c("#184f95", "#b02a2a")[seq_len(uniqueN(rates$sample))],
                      limits = levels(rates$sample), name = NULL) +
    labs(title = "Per-allele chrX rate under each SNP set",
         subtitle = paste("Informative chrX UMIs of ONE allele per 1000",
                          "informative autosomal UMIs, with a tile bootstrap",
                          "95% CI.\n'x' is the whole chrX as counted (Xist span",
                          "already masked in the SNP bed); 'noxic' additionally",
                          "masks Xist +- 500 kb."),
         x = NULL, y = "UMIs per 1000 autosomal UMIs",
         caption = paste("Both alleles are shown separately and always. The",
                         "ratio alone conflates a numerator change with a",
                         "denominator change.")) +
    th + theme(axis.text.x = element_text(angle = 30, hjust = 1))
)

if (!is.null(fc)) {
  print(
    ggplot(fc, aes(set, fc, colour = strain)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
      geom_pointrange(aes(ymin = lo, ymax = hi),
                      position = position_dodge(width = 0.5)) +
      scale_colour_manual(values = c("Bl6 (A1)" = "#b02a2a",
                                     "CAST (A2)" = "#184f95"), name = NULL) +
      coord_flip() +
      labs(title = sprintf("Fold change in the per-allele rate, by SNP set"),
           subtitle = paste("Tile bootstrap 95% CI.",
                            "n = 1 animal per age: this compares two sections,",
                            "not two ages."),
           x = NULL, y = "fold change (older / younger)") +
      th
  )
}
invisible(dev.off())

prov <- data.table(
  k = c("script", "run_at", "tile_um", "samples", "sets", "n_boot",
        "min_tile_umi", "counts_md5s"),
  v = c("spatial/ase_mask_comparison.R",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"), TILE_UM,
        paste(levels(tl$sample), collapse = ","),
        paste(SETS, collapse = ","), N_BOOT, MIN_TILE_UMI,
        paste(vapply(levels(tl$sample), function(s) unname(tools::md5sum(
          file.path(BASE, "ase", s, "bin_allele_counts.tsv"))), ""),
          collapse = ",")))
setnames(prov, c("key", "value"))
fwrite(prov, file.path(OUT_DIR,
                       sprintf("provenance_mask_comparison_%dum.tsv", TILE_UM)),
       sep = "\t")

msg("\nWrote mask_comparison_*_%dum.{csv,pdf} to %s", TILE_UM, OUT_DIR)
