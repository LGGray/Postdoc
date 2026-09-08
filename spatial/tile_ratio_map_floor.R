# Redraw the tile maps on tiles deep enough to carry a ratio, and show first the
# depth distribution that the floor has to be chosen against.
#
#   page 1   informative chrX units per tile, one histogram per sample, with the
#            floor marked and the retained fraction of TILES and of UNITS stated
#            separately. The second number is the one that matters for the
#            pooled statistic; the first is the one that matters for the map.
#   page 2   pooled ratio as a function of the floor. If the line is flat, the
#            floor is cosmetic and changes what the map LOOKS like without
#            changing what it SAYS - which is the claim this script makes.
#   page 3+  the panels of tile_ratio_map.R, on the tiles that clear the floor.
#
# WHY A FLOOR. tile_ratio_map.R applies none: any tile with x_n > 0 gets a ratio
# (its line 333, fifelse(x_n > 0, x_a1 / x_n, NA)). A tile holding two
# informative reads therefore takes a ratio of 0.0, 0.5 or 1.0 and is painted at
# full saturation, because the extreme OCM bins are the loudest colours on the
# ramp. Those tiles carry no evidence and a disproportionate share of the visual
# impression.
#
# WHY 20. It is the number tile_ratio_map.R's own panel_depth already names:
# ~20 INDEPENDENT observations to separate a ratio of 0.8 from 0.5. Read it with
# that adjective attached. At COUNT_UNIT=read a unit is a read, reads of one
# molecule are not independent observations of its allele, and the measured
# duplication on this data is 2.5-3.8x on average and up to 34x at a pile-up
# locus. So a 20-READ floor is 5-8 molecules and removes the indefensible tiles
# rather than delivering the stated power; MIN_X_UNITS=50 is the honest
# read-level stand-in for 20 molecules. Both are one environment variable, and
# page 2 shows what either does to the answer.
#
# USAGE (cluster; seurat_env, NOT the repo-wide RNAseq env)
#   conda activate seurat_env
#   SNP_LABEL=dup Rscript ~/Postdoc/spatial/tile_ratio_map_floor.R
#   SNP_LABEL=dup MIN_X_UNITS=50 Rscript ~/Postdoc/spatial/tile_ratio_map_floor.R
#   COUNT_UNIT=molecule PYSAM tree: set SNP_LABEL to match the directory suffix
#
# Every other knob - BASE, SAMPLES, TILE_UM, SNP_LABEL, COUNT_UNIT, the colour
# ramp, the panels - comes from tile_ratio_map.R and is deliberately NOT
# restated here, so the two scripts cannot drift apart.

.libPaths(c("~/R/matrix-dev", .libPaths()))

##### ------------------- source the loaders ------------------ #####
# tile_ratio_map.R next to this file, not a copy of it. TILE_RATIO_MAP_LIB makes
# it define collect_sample() and every panel_* without drawing its own figures
# or touching its own OUT_PDF.
# SPATIAL_DIR first, because --file= names the OUTER script when this file is
# itself sourced (a test harness, or an interactive session) and would then
# resolve to the wrong directory. Explicit beats inferred; inferred beats a
# hardcoded home path.
this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
SPATIAL_DIR <- Sys.getenv("SPATIAL_DIR")
if (!nzchar(SPATIAL_DIR)) {
  SPATIAL_DIR <- if (length(this_file)) dirname(normalizePath(this_file[1]))
                 else "~/Postdoc/spatial"
}
if (!file.exists(file.path(SPATIAL_DIR, "tile_ratio_map.R")))
  stop(sprintf("no tile_ratio_map.R in %s - set SPATIAL_DIR", SPATIAL_DIR))
Sys.setenv(TILE_RATIO_MAP_LIB = "1")
source(file.path(SPATIAL_DIR, "tile_ratio_map.R"))

##### ---------------------- CONFIG ---------------------- #####
# The floor, in whatever COUNT_UNIT the tree was scored in.
MIN_X_UNITS <- as.integer(Sys.getenv("MIN_X_UNITS", "20"))
if (is.na(MIN_X_UNITS) || MIN_X_UNITS < 1) stop("MIN_X_UNITS must be a positive integer")

# Alongside the unfloored figures, never on top of them: the whole point is to
# hold the two next to each other.
OUT_PDF <- Sys.getenv("OUT_PDF",
                      file.path(BASE, "ase",
                                sprintf("tile_ratio_map_%dum%s_min%d.pdf",
                                        TILE_UM, SUF, MIN_X_UNITS)))

# Floors for the page-2 sweep. 1 rather than 0 because "no floor" already means
# "at least one informative unit" - that is what x_n > 0 does.
SWEEP <- c(1, 2, 5, 10, 20, 50, 100)

# A fourth grey, darker than COL_NA (grey70, pending) and much darker than
# COL_FOOT (never submitted), so the calls panel can say "below floor" as its
# own category instead of borrowing one of the other two.
COL_FLOOR <- "#7d7b73"

##### ------------------- apply the floor ------------------- #####

# Re-derive everything downstream of the ratio on the kept tiles. Not a filter
# on the drawn data: the null band, and therefore every z and every call, is
# estimated FROM the tile set, so a floor that changed which tiles are drawn but
# left auto_sd fitted on all of them would report calls against a null the
# figure no longer shows. Mirrors tile_ratio_map.R's lines 335-372.
apply_floor <- function(d) {
  d <- copy(d)
  d[, shallow := !is.na(x_n) & x_n < MIN_X_UNITS]

  # Pending tiles are grey in the ratio panels and so are below-floor tiles, so
  # on an unfinished run the two are conflated and the map cannot be read. Say
  # so rather than drawing it quietly; use tile_ratio_map.R until the run ends.
  n_pending <- sum(d$submitted & is.na(d$x_n))
  if (n_pending > 0L) {
    msg("  WARNING: %d tiles are submitted but not yet scored. They and the", n_pending)
    msg("           below-floor tiles are BOTH grey in the ratio panels, so this")
    msg("           figure cannot distinguish them. Plot a finished run.")
  }

  n_scored <- sum(!is.na(d$x_n) & d$x_n > 0)
  n_keep   <- sum(!d$shallow & !is.na(d$x_n) & d$x_n > 0)
  u_all    <- sum(d$x_n, na.rm = TRUE)
  u_keep   <- sum(d[shallow == FALSE]$x_n, na.rm = TRUE)
  msg("  floor %d %s: %d of %d scored tiles kept (%.1f%%), holding %.1f%% of the chrX %s",
      MIN_X_UNITS, UNIT_N, n_keep, n_scored, 100 * n_keep / max(n_scored, 1),
      100 * u_keep / max(u_all, 1), UNIT_N)

  d[shallow == TRUE, x_ratio := NA_real_]

  # The null band, refitted. Autosomal spread on the KEPT tiles: a shallow-chrX
  # tile is a shallow tile generally, so leaving them in inflates the observed
  # variance and every |z| that is measured against it.
  sc <- d[shallow == FALSE & !is.na(a_ratio) & a_n > 0]
  if (nrow(sc) >= MIN_TILES_FOR_SD) {
    v_obs   <- var(sc$a_ratio)
    v_binom <- mean(sc$a_ratio * (1 - sc$a_ratio) / sc$a_n)
    auto_sd <- sqrt(max(v_obs - v_binom, 0))
    msg("  refitted null band %.4f on %d kept tiles (was %.4f on all of them)",
        auto_sd, nrow(sc), d$auto_sd[1])
  } else {
    auto_sd <- AUTO_SD_FALLBACK
    msg("  only %d kept tiles - using the %.3f fallback null band", nrow(sc), auto_sd)
  }
  d[, auto_sd := auto_sd]
  d[, se := sqrt(a_ratio * (1 - a_ratio) / x_n + auto_sd^2)]
  d[, z  := (x_ratio - a_ratio) / se]
  d[, x_bin := cut(x_ratio, breaks = OCM_BREAKS, include.lowest = TRUE,
                   right = TRUE, labels = OCM_LABELS)]
  # "below floor" ahead of the pending/not-submitted tests, so a tile that is
  # both shallow and unsubmitted is named by the reason this figure applied.
  d[, call := fcase(
      shallow,                          "below floor",
      is.na(x_ratio) &  submitted,      "pending",
      is.na(x_ratio) & !submitted,      "not submitted",
      z >  Z_CALL,                      "Bl6-skewed",
      z < -Z_CALL,                      "CAST-skewed",
      default = "mixed")]
  d
}

##### -------------------- new panels -------------------- #####

# Per-sample retained fractions, for the page-1 subtitle and the log.
retained <- function(out) {
  out[!is.na(x_n) & x_n > 0,
      .(tiles = .N, kept = sum(x_n >= MIN_X_UNITS),
        units = sum(x_n), units_kept = sum(x_n[x_n >= MIN_X_UNITS]),
        med = as.integer(median(x_n)), q1 = as.integer(quantile(x_n, 0.25)),
        q3 = as.integer(quantile(x_n, 0.75))), by = sample]
}

# Page 1. Log x because tile depth spans three or four orders of magnitude here;
# on a linear axis the shallow tiles that the floor is about are a single bar
# against the origin, which is precisely the detail being decided.
panel_depth_hist <- function(out) {
  r <- retained(out)
  sub <- paste(sprintf("%s: median %d (IQR %d-%d), %d of %d tiles and %.1f%% of %s clear %d",
                       r$sample, r$med, r$q1, r$q3, r$kept, r$tiles,
                       100 * r$units_kept / r$units, UNIT_N, MIN_X_UNITS),
               collapse = "\n")
  ggplot(out[!is.na(x_n) & x_n > 0], aes(x_n)) +
    geom_histogram(bins = 60, fill = "#3B5FB6", colour = NA) +
    geom_vline(xintercept = MIN_X_UNITS, linetype = 2, colour = "#8B1913") +
    annotate("text", x = MIN_X_UNITS, y = Inf,
             label = sprintf(" floor = %d", MIN_X_UNITS),
             hjust = 0, vjust = 1.5, size = 2.8, colour = "#8B1913") +
    scale_x_log10(breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000),
                  labels = scales::comma) +
    facet_wrap(~sample, ncol = 1, scales = "free_y") +
    labs(title = sprintf("informative chrX %s per %d um tile", UNIT_N, TILE_UM),
         subtitle = sub, x = sprintf("chrX %s per tile (log scale)", UNIT_N),
         y = "tiles",
         caption = if (identical(COUNT_UNIT, "read"))
           paste("These are READS. Reads of one molecule are not independent",
                 "observations of its allele, so a 20-read floor is 5-8\nmolecules",
                 "at the 2.5-3.8x duplication measured here. It removes tiles that",
                 "cannot carry a ratio at all;\nit does not deliver the ~20",
                 "independent observations that separate 0.8 from 0.5.") else
           paste("Units are molecules, so the floor is in independent observations",
                 "and the ~20 power figure applies directly.")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8, colour = "#52514e"),
          plot.caption  = element_text(size = 7, colour = "#52514e", hjust = 0))
}

# Page 2. The floor is a decision about what to DRAW; this is the check that it
# is not also a decision about what the data say. Pooled ratio is unit-weighted,
# so dropping tiles that hold a few percent of the units should move it by
# almost nothing - and if it does move, that is the finding, not the map.
panel_floor_sweep <- function(out) {
  sw <- rbindlist(lapply(SWEEP, function(f)
    out[!is.na(x_n) & x_n >= f,
        .(floor = f, tiles = .N, pooled = sum(x_a1) / sum(x_n)), by = sample]))
  ggplot(sw, aes(floor, pooled, colour = sample)) +
    geom_line() + geom_point(size = 1.6) +
    geom_vline(xintercept = MIN_X_UNITS, linetype = 2, colour = "#52514e") +
    geom_text(aes(label = tiles), vjust = -1, size = 2.4, show.legend = FALSE) +
    scale_x_log10(breaks = SWEEP) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = c("#184f95", "#b02a2a"), name = NULL) +
    labs(title = sprintf("does the floor change the answer? pooled chrX ratio vs floor"),
         subtitle = sprintf("dashed line is the floor this figure uses (%d %s); labels are the tiles retained",
                            MIN_X_UNITS, UNIT_N),
         x = sprintf("minimum chrX %s per tile (log scale)", UNIT_N),
         y = sprintf("pooled Bl6 (A1) fraction, %s", UNIT_WEIGHT),
         caption = paste("A flat line means the floor is cosmetic: it removes tiles",
                         "that cannot carry a ratio without moving the\npooled",
                         "statistic, because that statistic is",
                         paste0(UNIT_WEIGHT, ".") ,
                         "A rising or falling line means the shallow tiles\nwere",
                         "not merely noisy but biased, which is a result and",
                         "belongs in the text.")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          plot.subtitle = element_text(size = 8, colour = "#52514e"),
          plot.caption  = element_text(size = 7, colour = "#52514e", hjust = 0))
}

##### ---------------- panels that must be re-said ---------------- #####

# tile_ratio_map.R's wording assumes grey means "submitted, not yet scored".
# With a floor on a finished run grey means "scored, below the floor", and that
# is the sentence the reader needs. Redefining it here rather than editing the
# original keeps the unfloored figures exactly as they were; the panel functions
# resolve it at call time, so they pick this up.
scored_line <- function(n_ok, n_sub, extra = "") {
  sprintf("%d tiles clear the %d-%s floor; grey = scored but below it%s",
          n_ok, MIN_X_UNITS, UNIT_1, extra)
}

# Same panel as the original with "below floor" added as its own key. Copied
# rather than parameterised because the legend is the part that changes, and a
# reader comparing the two figures should be able to see the difference here.
panel_call <- function(d) {
  lv <- c("Bl6-skewed", "CAST-skewed", "mixed", "below floor", "pending",
          "not submitted")
  d2 <- copy(d)[, call := factor(call, levels = lv)]
  base_map(d2) +
    geom_tile(data = d2[call != "not submitted"], aes(fill = call),
              width = d$side[1], height = d$side[1], colour = NA) +
    scale_fill_manual(values = c("Bl6-skewed" = COL_BL6, "CAST-skewed" = COL_CAST,
                                 "mixed" = COL_MID, "below floor" = COL_FLOOR,
                                 "pending" = COL_NA, "not submitted" = COL_FOOT),
                      drop = FALSE, name = NULL) +
    guides(fill = guide_legend(override.aes = list(colour = "#c3c2b7"))) +
    labs(title = sprintf("%s - where there is DEPTH to call a skew, not where skew is",
                         d$sample[1]),
         subtitle = sprintf("|z| > %d vs this sample's own %.3f autosomal null band, refitted on the tiles above the floor: %s\n%s",
                            Z_CALL, d$auto_sd[1],
                            paste(sprintf("%s %d", lv[1:4],
                                          sapply(lv[1:4], function(l) sum(d2$call == l, na.rm = TRUE))),
                                  collapse = ", "),
                            depth_confound(d)),
         caption = paste("STILL A POWER MAP, even with the floor. z = (x_ratio - a_ratio)/se and se shrinks with depth,",
                         "so a tile is\ncoloured when it is DEEP, not when it is skewed. The floor removes the tiles",
                         "with no power at all;\nit does not make the remaining colour a statement about biology. See the",
                         "beta-binomial null of\nNEXT_ANALYSIS.md task 7 for a panel that would answer that."))
}

##### ------------------------- run ------------------------- #####

msg("Tile size %d um, samples: %s, SNP mask: %s, floor %d %s",
    TILE_UM, paste(SAMPLES, collapse = ", "), SNP_LABEL, MIN_X_UNITS, UNIT_N)

all_raw <- list(); all_flt <- list()
for (s in SAMPLES) {
  msg("[%s]", s)
  d <- tryCatch(collect_sample(s),
                error = function(e) { msg("  %s", conditionMessage(e)); NULL })
  if (is.null(d)) next
  if (!any(!is.na(d$x_ratio))) { msg("  no tile has a chrX ratio - skipping"); next }
  all_raw[[s]] <- d
  all_flt[[s]] <- apply_floor(d)
}

if (!length(all_flt)) {
  msg("Nothing to plot.")
} else {
  ord <- function(x) {
    x[, sample := factor(sample, levels = intersect(SAMPLES, unique(sample)))]
    x
  }
  raw <- ord(rbindlist(all_raw, fill = TRUE))
  flt <- ord(rbindlist(all_flt, fill = TRUE))

  pdf(OUT_PDF, width = 9, height = 8)
  print(panel_depth_hist(raw))
  print(panel_floor_sweep(raw))
  for (s in names(all_flt)) {
    d <- all_flt[[s]]
    print(panel_ratio_ocm(d, he = TRUE))
    print(panel_ratio_ocm(d, he = FALSE))
    print(panel_ratio(d, he = FALSE))
    print(panel_auto(d))
    print(panel_auto_zoom(d))
    print(panel_depth(d))
    print(panel_call(d))
  }
  print(panel_hist(flt))
  print(panel_violin(flt))
  invisible(dev.off())

  csv <- sub("\\.pdf$", ".csv", OUT_PDF)
  fwrite(flt[, .(sample, tile, x, y, n_bins, x_a1, x_a2, x_n, a_a1, a_a2, a_n,
                 x_ratio, x_bin, a_ratio, z, call, submitted, shallow)], csv)

  snp_bed <- file.path(dirname(BASE), "GRCm39",
                       sprintf("SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_%s.bed",
                               SNP_LABEL))
  prov <- data.table(
    k = c("script", "run_at", "tile_um", "samples", "snp_label", "snp_bed",
          "snp_bed_md5", "annotation", "count_unit", "min_x_units", "z_call",
          "auto_sd_refitted", "tiles_scored", "tiles_kept",
          "pooled_x_all", "pooled_x_kept"),
    v = c("spatial/tile_ratio_map_floor.R",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          TILE_UM, paste(levels(flt$sample), collapse = ","), SNP_LABEL, snp_bed,
          if (file.exists(snp_bed)) unname(tools::md5sum(snp_bed)) else
            "bed not readable from here",
          ANNOT_BASE, COUNT_UNIT, MIN_X_UNITS, Z_CALL,
          paste(sprintf("%s=%.4f", names(all_flt),
                        vapply(all_flt, function(d) d$auto_sd[1], 0)), collapse = ","),
          paste(sprintf("%s=%d", names(all_raw),
                        vapply(all_raw, function(d) sum(!is.na(d$x_ratio)), 0L)), collapse = ","),
          paste(sprintf("%s=%d", names(all_flt),
                        vapply(all_flt, function(d) sum(!is.na(d$x_ratio)), 0L)), collapse = ","),
          paste(sprintf("%s=%.4f", names(all_raw),
                        vapply(all_raw, function(d) sum(d$x_a1, na.rm = TRUE) /
                                                    sum(d$x_n, na.rm = TRUE), 0)), collapse = ","),
          paste(sprintf("%s=%.4f", names(all_flt),
                        vapply(all_flt, function(d) {
                          k <- d[shallow == FALSE]
                          sum(k$x_a1, na.rm = TRUE) / sum(k$x_n, na.rm = TRUE)
                        }, 0)), collapse = ",")))
  setnames(prov, c("key", "value"))
  fwrite(prov, sub("\\.pdf$", "_provenance.tsv", OUT_PDF), sep = "\t")

  msg("\nWrote %s\n       %s\n       %s", OUT_PDF, csv,
      sub("\\.pdf$", "_provenance.tsv", OUT_PDF))

  # The table the floor has to survive: if pooled_kept differs from pooled_all
  # by more than rounding, the shallow tiles were biased and not just noisy.
  msg("\nPooled chrX ratio, all scored tiles vs tiles clearing the floor:")
  cmp <- merge(
    raw[!is.na(x_ratio), .(tiles_all = .N, units_all = sum(x_n),
                           pooled_all = round(sum(x_a1) / sum(x_n), 4)), by = sample],
    flt[!is.na(x_ratio), .(tiles_kept = .N, units_kept = sum(x_n),
                           pooled_kept = round(sum(x_a1) / sum(x_n), 4)), by = sample],
    by = "sample")
  cmp[, `:=`(d_pooled = round(pooled_kept - pooled_all, 4),
             pct_units_kept = round(100 * units_kept / units_all, 1))]
  print(cmp[order(sample)])
}
