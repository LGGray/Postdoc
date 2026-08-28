# Where on chrX do the excess UMIs sit? (Task 1 of spatial/NEXT_ANALYSIS.md)
#
# The single analysis that settles the headline per-allele-rate result. The chrX
# aggregate everything else is built on is counted against whole-chromosome
# intervals with no gene requirement, so gene identity - and position - is
# discarded before any downstream script sees it. That is why a +45% CAST rate
# carried by ~10% of tiles cannot currently be told apart from a locus artefact.
#
# Inputs, per sample, both written by ase_bin_allele_counts.py:
#   ase/<sample>/bin_allele_counts.tsv                 per-2um-bin totals
#   ase/<sample>/chrX_windows_<TILE_UM>um.tsv.gz       --window-out, allele-split
#                                                      per 100 kb window x tile
#
# What it produces, in ase/chrX_window_profile_<TILE_UM>um.{pdf,csv}:
#
#   1. A positional profile of chrX UMIs along the chromosome, allele-split,
#      for the top X-share decile of tiles against all the others, per sample.
#      Rates, not counts: UMIs per 1000 informative autosomal UMIs in the same
#      tiles, so two groups of very different depth are comparable. This is the
#      per-allele statistic the analysis plan makes primary, resolved by
#      position instead of summed over the chromosome.
#
#   2. The concentration verdict. Excess is defined per window and per allele as
#      rate(top decile) - rate(rest), and the windows are ranked by it. If a
#      handful of windows carry most of the excess, the +45% is a locus artefact
#      and the Xic is the prime suspect (mm39 chrX:100-104 Mb). If the excess is
#      spread in proportion to each window's normal output - i.e. the ratio of
#      the two profiles is flat - it is an Xi-wide effect.
#
#   3. The same decomposition for the 9w -> 78w change, which is the comparison
#      the paper would actually rest on. n = 1 per age, so this is a description
#      of two sections, not an age effect.
#
# Run:  sbatch ~/Postdoc/slurm/spatial_window_profile.slurm 64
#  or:  TILE_UM=64 Rscript ~/Postdoc/spatial/ase_chrX_window_profile.R

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

# The window table must have been written at the SAME tile size, or the tile
# grids do not line up and the join silently drops nearly everything.
WINDOW_TSV <- function(s) file.path(BASE, "ase", s,
                                    sprintf("chrX_windows_%dum.tsv.gz", TILE_UM))
COUNTS_TSV <- function(s) file.path(BASE, "ase", s, "bin_allele_counts.tsv")

# Tiles are stratified by X-share, x_n / (x_n + a_n). The top decile is the
# group under suspicion: in 78w those tiles carry 198 chrX UMIs against 3879
# autosomal, 3.4x the chrX depth of other tiles at LOWER autosomal depth, which
# is not a transcriptional state anything explains.
N_STRATA <- 10L
TOP_STRATUM <- N_STRATA          # which stratum is "suspect"

# A tile needs some autosomal depth for its X-share to mean anything; with
# a_n = 3 the share is quantised to thirds.
MIN_AUTO_UMI <- 200L

# mm39. Xist sits at chrX:102.28 Mb; the Xic and its regulators span roughly
# 100-104 Mb. Drawn as a band on every positional panel and reported separately
# in the concentration table, because it is the specific hypothesis.
XIC_FROM <- 100e6
XIC_TO   <- 104e6

# How many top windows the descriptive tables are stated over. 10 of ~1700 is
# 0.6% of the chromosome.
TOP_K <- 10L

# The verdict is NOT a fold-change threshold on the top-K share. Ranking windows
# by excess and then reporting the top few overstates concentration by
# construction - on the synthetic control, 200 windows of pure Poisson noise with
# no spike at all put 13% of the "excess" in the top 10 windows, 3.1x their share
# of normal output. Any fixed multiple either passes that or fails a real effect.
#
# So the test is exact and per window instead. Condition on each window's total
# across the two tile groups: under the Xi-wide hypothesis the two groups have
# the same positional profile and differ only by one global factor, so a UMI in
# window w falls in the top decile with probability N_top / (N_top + N_rest)
# regardless of w. Enrichment is then the upper tail of a binomial, BH-corrected
# across windows. The global factor - the thing the Xi-wide hypothesis is about -
# is absorbed by the conditioning and cannot make a window look significant.
FDR <- 0.01
# Concentrated if few windows are significant AND they carry most of the excess.
CONC_MAX_WINDOWS <- 5L
CONC_MIN_SHARE   <- 0.5
##### ------------------------------------------------------ #####

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), ...)
  else fread(path, ...)
}

k_bins <- TILE_UM %/% 2L

# ------------------------------------------------------------------- load

load_sample <- function(smp) {
  ct <- COUNTS_TSV(smp); wt <- WINDOW_TSV(smp)
  if (!file.exists(ct)) { msg("  no count table %s - skipping %s", ct, smp); return(NULL) }
  if (!file.exists(wt)) {
    msg("  no window table %s", wt)
    msg("  Re-run the counting pass with --window-out; see slurm/spatial_ase_sweep.slurm.")
    return(NULL)
  }
  bins <- fread_any(ct)
  setnames(bins, tolower(names(bins)))
  bins[, `:=`(tile_row = array_row %/% k_bins, tile_col = array_col %/% k_bins)]
  tiles <- bins[, .(x_ref = sum(x_ref), x_alt = sum(x_alt),
                    a_ref = sum(a_ref), a_alt = sum(a_alt),
                    n_bins = .N), by = .(tile_row, tile_col)]
  tiles[, `:=`(x_n = x_ref + x_alt, a_n = a_ref + a_alt)]
  tiles[, x_share := fifelse(x_n + a_n > 0, x_n / (x_n + a_n), NA_real_)]

  n0 <- nrow(tiles)
  tiles <- tiles[a_n >= MIN_AUTO_UMI & !is.na(x_share)]
  msg("  %s: %d tiles, %d with >= %d autosomal UMIs", smp, n0, nrow(tiles),
      MIN_AUTO_UMI)
  if (!nrow(tiles)) return(NULL)

  # Deciles of X-share. Ties at the boundary are broken by rank so the strata
  # have equal counts even when the share is quantised at low depth.
  tiles[, stratum := as.integer(cut(frank(x_share, ties.method = "first"),
                                    breaks = N_STRATA, labels = FALSE))]
  tiles[, group := fifelse(stratum == TOP_STRATUM, "top decile", "rest")]

  win <- fread_any(wt)
  setnames(win, tolower(names(win)))
  need <- c("chrom", "win_start", "tile_row", "tile_col", "ref", "alt")
  if (!all(need %in% names(win))) {
    stop("window table ", wt, " lacks ", paste(setdiff(need, names(win)),
                                               collapse = ", "))
  }
  win <- merge(win, tiles[, .(tile_row, tile_col, group, stratum, a_n)],
               by = c("tile_row", "tile_col"))
  if (!nrow(win)) {
    stop("No window row joined a tile for ", smp, ".\n",
         "  The window table was almost certainly written at a different ",
         "--window-tile-um\n  than TILE_UM = ", TILE_UM, ".")
  }
  win[, sample := smp]
  tiles[, sample := smp]
  list(tiles = tiles, win = win)
}

msg("Tile size %d um, samples: %s", TILE_UM, paste(SAMPLES, collapse = ", "))
loaded <- lapply(SAMPLES, load_sample)
names(loaded) <- SAMPLES
loaded <- loaded[!vapply(loaded, is.null, TRUE)]
if (!length(loaded)) stop("No sample had both a count table and a window table")

tiles <- rbindlist(lapply(loaded, `[[`, "tiles"))
win   <- rbindlist(lapply(loaded, `[[`, "win"))
tiles[, sample := factor(sample, levels = intersect(SAMPLES, unique(sample)))]
win[, sample := factor(sample, levels = intersect(SAMPLES, unique(sample)))]

# Sanity: reproduce the tile-level numbers the analysis plan quotes, so a
# mismatch between the two tables shows up here rather than as a wrong figure.
msg("\nTile groups (the suspects are the top X-share decile):")
msg_table(tiles[, .(tiles = .N,
                    x_n = as.integer(median(x_n)),
                    a_n = as.integer(median(a_n)),
                    x_share = round(median(x_share), 4),
                    x_ratio = round(sum(x_ref) / sum(x_n), 4),
                    a_ratio = round(sum(a_ref) / sum(a_n), 4),
                    n_bins = as.integer(median(n_bins))),
                by = .(sample, group)][order(sample, -group)])

# The window table is a strict subset of the chrX counts in the main table -
# every chrX UMI with a decided allele lands in exactly one window - so the two
# must agree. They will not if the two files came from different runs.
chk <- merge(
  win[chrom == unique(win$chrom)[1], .(win_umi = sum(ref + alt)), by = sample],
  tiles[, .(tile_umi = sum(x_n)), by = sample], by = "sample")
chk[, agreement := round(win_umi / tile_umi, 4)]
msg("\nWindow table vs main table, chrX UMIs (should be ~1, below 1 only by the")
msg("tiles dropped for autosomal depth):")
msg_table(chk)

##### -------------------- positional profiles -------------------- #####

X_CHROM <- if ("chrX" %in% win$chrom) "chrX" else unique(win$chrom)[1]
if (X_CHROM != "chrX") msg("NOTE: no chrX in the window table, using %s", X_CHROM)

# Autosomal depth per group, the denominator of every rate below. Taken from the
# tiles, not from the window table: it has to be the same denominator the
# per-allele rate uses elsewhere, and the window table only holds chrX.
denom <- tiles[, .(auto_umi = sum(a_n), tiles = .N), by = .(sample, group)]

prof <- win[chrom == X_CHROM,
            .(ref = sum(ref), alt = sum(alt)), by = .(sample, group, win_start)]
prof <- merge(prof, denom, by = c("sample", "group"))
prof[, `:=`(rate_ref = 1000 * ref / auto_umi,
            rate_alt = 1000 * alt / auto_umi,
            rate_n   = 1000 * (ref + alt) / auto_umi)]

# Windows with no UMI in a group are absent, not zero. Fill them, so a window
# that fires only in the top decile is visible as a spike against zero rather
# than dropping out of the comparison.
grid <- CJ(sample = unique(prof$sample), group = unique(prof$group),
           win_start = sort(unique(prof$win_start)), unique = TRUE)
prof <- merge(grid, prof, by = c("sample", "group", "win_start"), all.x = TRUE)
prof <- merge(prof[, !c("auto_umi", "tiles")], denom,
              by = c("sample", "group"))
for (cc in c("ref", "alt")) prof[is.na(get(cc)), (cc) := 0L]
prof[, `:=`(rate_ref = 1000 * ref / auto_umi,
            rate_alt = 1000 * alt / auto_umi,
            rate_n   = 1000 * (ref + alt) / auto_umi,
            in_xic   = win_start >= XIC_FROM & win_start < XIC_TO)]

fwrite(prof, file.path(OUT_DIR, sprintf("chrX_window_profile_%dum.csv", TILE_UM)))

##### ------------------- concentration of excess ------------------- #####

# One row per window per allele: the top decile's rate, the rest's rate, and the
# difference. Ranked by the difference, cumulated, so "the top 10 windows explain
# X% of the excess" is read straight off it.
wide <- dcast(prof, sample + win_start + in_xic ~ group,
              value.var = c("ref", "alt", "rate_ref", "rate_alt"))
setnames(wide, gsub(" ", "_", names(wide)))

concentration <- function(allele) {
  rt <- paste0("rate_", allele, "_top_decile")
  rr <- paste0("rate_", allele, "_rest")
  if (!all(c(rt, rr) %in% names(wide))) return(NULL)
  d <- wide[, .(sample, win_start, in_xic,
                r_top = get(rt), r_rest = get(rr))]
  d[is.na(r_top), r_top := 0]
  d[is.na(r_rest), r_rest := 0]
  d[, excess := r_top - r_rest]
  d[, allele := allele]
  d <- d[order(sample, -excess)]
  d[, rank := seq_len(.N), by = sample]
  # Denominators: the POSITIVE excess only. Negative-excess windows are real
  # information but they cannot be part of "which windows carry the excess", and
  # netting them off would let a big negative window inflate the share of the
  # positives past 100%.
  d[, pos_excess_total := sum(pmax(excess, 0)), by = sample]
  d[, cum_share := cumsum(pmax(excess, 0)) / pos_excess_total, by = sample]
  # The fair baseline: what share of the REST group's output those same windows
  # carry. The busiest windows are expected to carry the most of anything.
  d[, rest_total := sum(r_rest), by = sample]
  d[]
}

conc <- rbindlist(lapply(c("ref", "alt"), concentration))

# --- the exact per-window enrichment test ---------------------------------
#
# Counts, not rates: the test needs the raw numerator. n_tiles_top is carried
# alongside because a window whose excess comes from ONE tile is a different
# claim from one that fires across the whole decile, and the binomial - which
# treats UMIs as independent - cannot tell them apart on its own.
tiles_per_window <- win[chrom == X_CHROM & group == "top decile" & ref + alt > 0,
                        .(n_tiles_top = uniqueN(paste(tile_row, tile_col))),
                        by = .(sample, win_start)]

enrich <- rbindlist(lapply(c("ref", "alt"), function(al) {
  ct <- paste0(al, "_top_decile"); cr <- paste0(al, "_rest")
  if (!all(c(ct, cr) %in% names(wide))) return(NULL)
  d <- wide[, .(sample, win_start, in_xic, allele = al,
                c_top = as.numeric(get(ct)), c_rest = as.numeric(get(cr)))]
  d[is.na(c_top), c_top := 0]
  d[is.na(c_rest), c_rest := 0]
  d[, `:=`(N_top = sum(c_top), N_rest = sum(c_rest)), by = sample]
  d[, pi_top := N_top / (N_top + N_rest)]
  d[, m := c_top + c_rest]
  # Upper tail P(X >= c_top). Enrichment in the suspect group only: a window
  # that is DEPLETED there is not what task 1 is asking about.
  d[, p := fifelse(m > 0, pbinom(c_top - 1, m, pi_top, lower.tail = FALSE), 1)]
  d[, q := p.adjust(p, method = "BH"), by = sample]
  d[, obs_over_exp := fifelse(m > 0, c_top / (m * pi_top), NA_real_)]
  d[]
}))
enrich <- merge(enrich, tiles_per_window, by = c("sample", "win_start"),
                all.x = TRUE)
enrich[is.na(n_tiles_top), n_tiles_top := 0L]
conc <- merge(conc, enrich[, .(sample, win_start, allele, c_top, c_rest, m,
                               obs_over_exp, p, q, n_tiles_top)],
              by = c("sample", "win_start", "allele"), all.x = TRUE)
fwrite(conc[order(sample, allele, -excess)],
       file.path(OUT_DIR, sprintf("chrX_window_excess_%dum.csv", TILE_UM)))

# min()/max() of an empty selection returns +-Inf as a double while a non-empty
# one returns the column's own type, which makes the by-group result types
# disagree and errors. Both are forced to double and to NA when there is nothing
# to take an extremum of.
extremum <- function(f, v) {
  v <- as.numeric(v)
  if (!length(v)) NA_real_ else f(v)
}
verdict <- conc[, {
  sig <- !is.na(q) & q < FDR
  pos <- pmax(excess, 0)
  .(n_windows_tested = .N,
    n_sig = sum(sig),
    sig_excess_share = round(sum(pos[sig]) / sum(pos), 4),
    sig_in_xic = sum(sig & in_xic),
    max_sig_over_exp = round(extremum(max, obs_over_exp[sig]), 2),
    min_tiles_in_sig = extremum(min, n_tiles_top[sig]),
    top_k_excess_share = round(sum(pos[rank <= TOP_K]) / sum(pos), 4),
    total_excess_per_1000 = round(pos_excess_total[1], 3))
}, by = .(sample, allele)]
verdict[, call := fifelse(
  n_sig > 0 & n_sig <= CONC_MAX_WINDOWS & sig_excess_share >= CONC_MIN_SHARE,
  "CONCENTRATED",
  fifelse(n_sig == 0 | sig_excess_share < CONC_MIN_SHARE, "spread", "mixed"))]

msg("\n--- per-window enrichment test, BH q < %.3g over %d windows ---",
    FDR, uniqueN(prof$win_start))
msg("Conditioning on each window's total absorbs the global rate factor, so a")
msg("significant window is one whose SHARE of chrX output differs, not one that")
msg("is merely busy. min_tiles_in_sig is how many distinct tiles the thinnest")
msg("significant window draws on - 1 means a single tile is the whole result.")
msg_table(verdict[order(sample, allele)])

# The Xic specifically, since that is the named hypothesis.
xic <- conc[, .(xic_excess_share = round(sum(pmax(excess, 0) * in_xic) /
                                           pos_excess_total[1], 4),
                xic_rest_share = round(sum(r_rest * in_xic) / rest_total[1], 4),
                xic_windows = sum(in_xic)),
            by = .(sample, allele)]
msg("\nXic (chrX:%.0f-%.0f Mb) share of the excess, against its share of normal output:",
    XIC_FROM / 1e6, XIC_TO / 1e6)
msg_table(xic[order(sample, allele)])

msg("\nTop %d windows by CAST (alt) excess:", TOP_K)
msg_table(conc[allele == "alt" & rank <= TOP_K][order(sample, rank),
               .(sample, rank, Mb = round(win_start / 1e6, 2), in_xic,
                 c_top, c_rest, obs_exp = round(obs_over_exp, 2),
                 q = signif(q, 2), n_tiles = n_tiles_top,
                 excess = round(excess, 4), cum_share = round(cum_share, 3))])

##### ------------------- the 9w -> 78w change ------------------- #####

age <- NULL
if (uniqueN(prof$sample) >= 2L) {
  lv <- levels(prof$sample)
  a1 <- lv[1]; a2 <- lv[length(lv)]
  ag <- dcast(prof[, .(sample, group, win_start, in_xic, rate_ref, rate_alt)],
              group + win_start + in_xic ~ sample,
              value.var = c("rate_ref", "rate_alt"))
  age <- rbindlist(lapply(c("ref", "alt"), function(al) {
    c1 <- paste0("rate_", al, "_", a1); c2 <- paste0("rate_", al, "_", a2)
    if (!all(c(c1, c2) %in% names(ag))) return(NULL)
    d <- ag[, .(group, win_start, in_xic, allele = al,
                r_young = get(c1), r_old = get(c2))]
    d[is.na(r_young), r_young := 0]
    d[is.na(r_old), r_old := 0]
    d[, delta := r_old - r_young]
    d <- d[order(group, -delta)]
    d[, rank := seq_len(.N), by = group]
    d[, cum_share := cumsum(pmax(delta, 0)) / sum(pmax(delta, 0)), by = group]
    d[]
  }))
  fwrite(age, file.path(OUT_DIR,
                        sprintf("chrX_window_age_delta_%dum.csv", TILE_UM)))
  msg("\n--- %s -> %s change, by window (n = 1 per age: two sections, not an age effect) ---",
      a1, a2)
  msg_table(age[rank <= TOP_K & allele == "alt"][order(group, rank),
                .(group, rank, Mb = round(win_start / 1e6, 2), in_xic,
                  r_young = round(r_young, 4), r_old = round(r_old, 4),
                  delta = round(delta, 4), cum_share = round(cum_share, 3))])
  msg("\nShare of the total CAST rate increase carried by the top %d windows:",
      TOP_K)
  msg_table(age[allele == "alt", .(
    top_k_share = round(sum(pmax(delta, 0)[rank <= TOP_K]) /
                          sum(pmax(delta, 0)), 4),
    xic_share = round(sum(pmax(delta, 0) * in_xic) / sum(pmax(delta, 0)), 4),
    total_increase_per_1000 = round(sum(pmax(delta, 0)), 3)), by = group])
}

##### --------------------------- figures --------------------------- #####

pdf(file.path(OUT_DIR, sprintf("chrX_window_profile_%dum.pdf", TILE_UM)),
    width = 10, height = 6.5)

xic_band <- annotate("rect", xmin = XIC_FROM / 1e6, xmax = XIC_TO / 1e6,
                     ymin = -Inf, ymax = Inf, fill = "#f2c9c9", alpha = 0.45)
theme_prof <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8, colour = "#52514e"),
        plot.caption = element_text(size = 7, colour = "#52514e", hjust = 0))

pl <- melt(prof, id.vars = c("sample", "group", "win_start", "in_xic"),
           measure.vars = c("rate_ref", "rate_alt"),
           variable.name = "allele", value.name = "rate")
pl[, allele := factor(fifelse(allele == "rate_ref", "Bl6 (Xa)", "CAST (Xi)"),
                      levels = c("Bl6 (Xa)", "CAST (Xi)"))]

print(
  ggplot(pl, aes(win_start / 1e6, rate, colour = group)) +
    xic_band +
    geom_line(linewidth = 0.35) +
    facet_grid(allele ~ sample, scales = "free_y") +
    scale_colour_manual(values = c("top decile" = "#b02a2a",
                                   "rest" = "#184f95"), name = NULL) +
    labs(title = sprintf("chrX UMIs along the chromosome, %d um tiles", TILE_UM),
         subtitle = paste("Rate, not count: UMIs per 1000 informative autosomal",
                          "UMIs in the same tiles. Shaded band is the Xic",
                          "(chrX:100-104 Mb, mm39)."),
         x = "chrX position (Mb, mm39)", y = "UMIs per 1000 autosomal UMIs",
         caption = paste("If the top decile's excess were an Xi-wide effect the",
                         "red line would sit above the blue one by a constant",
                         "FACTOR across the chromosome.",
                         "\nIf it is a locus artefact it will be a few spikes.")) +
    theme_prof
)

# The same thing as a ratio, which is the form that separates the two
# hypotheses: flat = Xi-wide, spiky = locus.
pl2 <- copy(wide)
for (al in c("ref", "alt")) {
  rt <- paste0("rate_", al, "_top_decile"); rr <- paste0("rate_", al, "_rest")
  if (all(c(rt, rr) %in% names(pl2))) {
    # A pseudocount at the scale of one UMI in the rest group, so a window that
    # is empty there gives a bounded ratio instead of an infinity.
    eps <- pl2[get(rr) > 0, min(get(rr))] / 2
    pl2[, (paste0("lfc_", al)) := log2((get(rt) + eps) / (get(rr) + eps))]
  }
}
pl2m <- melt(pl2, id.vars = c("sample", "win_start", "in_xic"),
             measure.vars = grep("^lfc_", names(pl2), value = TRUE),
             variable.name = "allele", value.name = "lfc")
pl2m[, allele := factor(fifelse(allele == "lfc_ref", "Bl6 (Xa)", "CAST (Xi)"),
                        levels = c("Bl6 (Xa)", "CAST (Xi)"))]
print(
  ggplot(pl2m[is.finite(lfc)], aes(win_start / 1e6, lfc)) +
    xic_band +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_point(aes(colour = in_xic), size = 0.5) +
    geom_smooth(method = "loess", span = 0.2, se = FALSE, colour = "#0b0b0b",
                linewidth = 0.4, formula = y ~ x) +
    facet_grid(allele ~ sample) +
    scale_colour_manual(values = c("FALSE" = "#52514e", "TRUE" = "#b02a2a"),
                        guide = "none") +
    labs(title = "Top X-share decile against the rest, per 100 kb window",
         subtitle = paste("log2 of the rate ratio. FLAT means the top decile is",
                          "the whole Xi turned up by one factor - a real",
                          "Xi-wide effect.\nSPIKES mean a locus artefact, and",
                          "the loess line will not move with them."),
         x = "chrX position (Mb, mm39)", y = "log2 (top decile / rest)") +
    theme_prof
)

print(
  ggplot(conc[allele == "alt"], aes(rank, cum_share, colour = sample)) +
    geom_hline(yintercept = c(0.5, 0.9), linetype = 2, colour = "grey60") +
    geom_line() +
    scale_x_log10() +
    labs(title = "How few windows carry the CAST excess",
         subtitle = paste("Cumulative share of the positive top-decile excess,",
                          "windows ranked by excess. A curve that reaches 0.5",
                          "within a few windows\nis a locus artefact; one that",
                          "climbs slowly over hundreds is an Xi-wide effect."),
         x = sprintf("window rank (of %d)", uniqueN(prof$win_start)),
         y = "cumulative share of excess") +
    theme_prof
)

# Top decile vs rest per window, on log rates. Proportionality is a straight
# line of slope 1 through the origin; a locus artefact is a few points off it.
print(
  ggplot(wide[rate_alt_rest > 0 & rate_alt_top_decile > 0],
         aes(rate_alt_rest, rate_alt_top_decile)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
    geom_point(aes(colour = in_xic), size = 0.6) +
    facet_wrap(~ sample) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c("FALSE" = "#52514e", "TRUE" = "#b02a2a"),
                        name = "in Xic") +
    labs(title = "CAST rate per window: top X-share decile vs the rest",
         subtitle = paste("Dashed line is equality. A uniform Xi-wide increase",
                          "is a cloud parallel to it and above it;",
                          "\na locus artefact is a few points far above it."),
         x = "rest, UMIs per 1000 autosomal", y = "top decile") +
    theme_prof
)

if (!is.null(age)) {
  print(
    ggplot(age[allele == "alt"], aes(win_start / 1e6, delta, colour = group)) +
      xic_band +
      geom_hline(yintercept = 0, colour = "grey40") +
      geom_point(size = 0.5) +
      facet_wrap(~ group, ncol = 1, scales = "free_y") +
      scale_colour_manual(values = c("top decile" = "#b02a2a",
                                     "rest" = "#184f95"), guide = "none") +
      labs(title = "Where the CAST rate increase lives, per window",
           subtitle = paste("Older minus younger, per 1000 autosomal UMIs.",
                            "n = 1 animal per age: this describes two sections,",
                            "not an age effect."),
           x = "chrX position (Mb, mm39)", y = "rate difference") +
      theme_prof
  )
}
invisible(dev.off())

##### --------------------------- verdict --------------------------- #####

msg("\n=== ACCEPTANCE CRITERION (task 1) ===")
for (i in seq_len(nrow(verdict))) {
  v <- verdict[i]
  msg("%s, %s allele: %d of %d windows enriched at q < %.3g, carrying %.1f%% of the excess (thinnest draws on %s tiles) -> %s",
      as.character(v$sample), v$allele, v$n_sig, v$n_windows_tested, FDR,
      100 * v$sig_excess_share,
      if (is.finite(v$min_tiles_in_sig)) as.character(v$min_tiles_in_sig) else "-",
      v$call)
}
alt_calls <- verdict[allele == "alt", call]
if (any(alt_calls == "CONCENTRATED")) {
  msg("\nCONCENTRATED. The excess sits in a handful of windows, so the +45%% CAST")
  msg("  rate is a locus artefact, not an Xi-wide effect.")
  msg("  -> Task 2 (mask the Xic properly) before anything else, and do not run")
  msg("     tasks 5-8 on the unmasked counts. Read the top-window table above:")
  msg("     if the windows are in the Xic band the mask fixes it; if they are")
  msg("     elsewhere, that locus needs its own explanation first. If n_tiles is")
  msg("     1 or 2, one tile is the entire result and no mask will help.")
} else if (any(alt_calls == "mixed")) {
  msg("\nMIXED. Several windows are individually enriched but they do not carry")
  msg("  most of the excess. Mask them (task 2 for the Xic, an explicit bed for")
  msg("  anything else), re-run, and see what is left before going on.")
} else {
  msg("\nSPREAD. No small set of windows carries the excess: it is distributed")
  msg("  across the chromosome roughly in proportion to normal output, which is")
  msg("  what a genuine Xi-wide effect looks like. It becomes the central result")
  msg("  - subject to task 3's imprinted control showing the per-UMI error floor")
  msg("  is not what is being measured, and to replication: n = 1 per age.")
}

# NOT data.table(key = ..., value = ...): `key` is data.table()'s own argument,
# so that form silently tries to set a key instead of making a column.
prov <- data.table(
  k = c("script", "run_at", "tile_um", "samples", "n_strata", "top_stratum",
        "min_auto_umi", "top_k", "fdr", "conc_max_windows", "conc_min_share",
        "xic_from", "xic_to", "window_tables", "window_table_md5s"),
  v = c("spatial/ase_chrX_window_profile.R",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"), TILE_UM,
        paste(names(loaded), collapse = ","), N_STRATA, TOP_STRATUM,
        MIN_AUTO_UMI, TOP_K, FDR, CONC_MAX_WINDOWS, CONC_MIN_SHARE,
        XIC_FROM, XIC_TO,
        paste(vapply(names(loaded), WINDOW_TSV, ""), collapse = ","),
        paste(vapply(names(loaded), function(s)
          unname(tools::md5sum(WINDOW_TSV(s))), ""), collapse = ",")))
setnames(prov, c("key", "value"))
fwrite(prov, file.path(OUT_DIR,
                       sprintf("provenance_chrX_window_profile_%dum.tsv", TILE_UM)),
       sep = "\t")

msg("\nWrote to %s:", OUT_DIR)
for (f in sprintf(c("chrX_window_profile_%dum.pdf", "chrX_window_profile_%dum.csv",
                    "chrX_window_excess_%dum.csv", "chrX_window_age_delta_%dum.csv",
                    "provenance_chrX_window_profile_%dum.tsv"), TILE_UM)) {
  if (file.exists(file.path(OUT_DIR, f))) msg("  %s", f)
}
