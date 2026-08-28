# Choose a tile size for the spatial allele-specific analysis.
#
# Input is the per-2um-bin allelic count table from ase_bin_allele_counts.py.
# Because tile sizes are nested (each is a whole number of 2um bins per side),
# every candidate tiling is a pure summation of that one table - so the sweep
# costs minutes, and sinto only ever gets run once, at the size chosen here.
#
# Three things are measured, and they disagree with each other on purpose:
#
#   1. Coverage      fraction of tiles reaching a usable UMI count. Rises with
#                    tile size. Sets the floor.
#   2. Dispersion    rho(s), the fraction of the maximum possible mosaic
#                    variance in allelic ratio that is real rather than
#                    sampling noise. Falls with tile size as tiles start
#                    averaging over more than one clonal patch. Sets the ceiling.
#   3. C(d)          probability that two informative UMIs d microns apart carry
#                    the same allele. Needs no tiling at all, so it estimates
#                    the patch size without the thing we are trying to choose
#                    contaminating the answer. This is the one to trust.
#
# Everything is computed for chrX and, identically, for the autosomal control.
# Autosomes are biallelic and should have no spatial structure, so autosomal
# rho and C(d) are the empirical null - they already contain the reference
# mapping bias from using the standard 10x B6 reference, plus whatever technical
# spatial structure the slide has. Only chrX above that null means anything.
#
# EVERY EXTRA COLUMN PAIR in the count table is carried through automatically.
# ase_bin_allele_counts.py --subset-bed writes <name>_ref / <name>_alt for any
# interval set, so escape vs non-escape chrX, Xic-masked chrX and the imprinted
# positive controls all appear here as additional curves with no code change.
# What each one is for:
#
#   escape / nonescape   Partitions the C(d) floor. C(4um) = 0.785 on chrX says
#                        ~12% of chrX UMI allele calls disagree with their own
#                        cell's XCI state. If dropping the escapees lifts the
#                        short range toward 1, that floor is escape and the
#                        global skew is ~0.88. If it does not, the floor is
#                        per-UMI assignment error and the chrX numbers built on
#                        it are measuring a noise rate.
#
#   imppat / impmat      The positive control that decides between those two.
#                        An imprinted locus is monoallelic in EVERY cell, so
#                        two UMIs from one should agree at any separation. The
#                        two sets are separate because H19/Igf2, Airn/Igf2r and
#                        Meg3/Dlk1 are reciprocally imprinted - pooling them
#                        would mix UMIs carrying opposite alleles and land C(d)
#                        near 0.5 on perfect data. See OCM_heart/core_escape_SNPs.R.
#
#   noxic / xic          chrX with Xist +- 500 kb masked, and the masked-out
#                        region on its own, for comparison against the
#                        Xist-gene-only mask the SNP bed applies.
#
# Set names come from the --subset-bed flags in slurm/spatial_ase_sweep.slurm.
# Column names are lowercased on read, so name them in lower case there to keep
# the two ends readable together.
#
# NOTE on the imprinted control: it is a test of the per-UMI ERROR FLOOR, not of
# spatial resolution. Imprinting is uniform across the tissue, so C_imprinted(d)
# is expected flat and high, NOT decaying - there is no patch structure for it
# to resolve. A flat curve near 1 is the pass condition; decay would be the
# surprise. The resolution question the analysis plan also asks of it cannot be
# answered this way, because no monoallelic-and-spatially-patchy locus exists to
# calibrate against; the chrX excess over the autosomal null remains the only
# handle on patch structure.
#
# Run on the cluster:  sbatch ~/Postdoc/slurm/spatial_ase_sweep.slurm

.libPaths(c("~/R/matrix-dev", .libPaths()))

# print() writes to stdout and message() to stderr, and SLURM sends those to
# different files. Route tables through message() so the log reads in order.
msg_table <- function(x) message(paste(capture.output(print(x)), collapse = "\n"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

##### ---------------------- CONFIG ---------------------- #####
BASE      <- Sys.getenv("BASE",
  "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial")
# The SLURM array sets SAMPLE; the default is what you get running by hand.
SAMPLE    <- Sys.getenv("SAMPLE", "9w")
IN_TSV    <- file.path(BASE, "ase", SAMPLE, "bin_allele_counts.tsv")
# Every in-tissue 2um bin, including the ones with no informative UMI. Needed
# for the coverage denominator - IN_TSV only holds bins that carry data, and
# using that as the denominator would report near-perfect coverage at 8um.
TISSUE_TSV <- paste0(IN_TSV, ".tissue_bins.tsv.gz")
OUT_DIR   <- file.path(BASE, "ase", SAMPLE)

# Candidate tile sizes, in microns. Every one is a whole number of 2um bins per
# side, so the tilings nest exactly and the sweep is strict aggregation with no
# re-binning artefacts. 2um and 4um are in the list to make the point that they
# are empty, not because they are plausible.
SIZES_UM <- c(2, 4, 8, 16, 32, 64, 96, 128, 192, 256, 384, 512)

# Minimum informative UMIs for a tile to enter the dispersion estimate. Below
# ~5 the per-tile ratio is almost pure noise and the moment estimator becomes
# unstable, though the beta-binomial MLE handles it better.
MIN_UMI <- 5

# Reporting thresholds for the coverage curve. From a two-sided binomial test
# against 0.5 at alpha = 0.05, 80% power: ~10 UMIs detects an AR of 0.9, ~20
# detects 0.8, ~50 detects 0.7. Calling which X is active in a near-clonal
# patch is the 0.9 case; measuring escape is the 0.7 case and will not work
# per-tile at any size, which is what the domain pseudobulk is for.
COVERAGE_TARGETS <- c(10, 20, 50)

# Pair-correlation sampling. Pairs are drawn per distance stratum rather than
# uniformly at random, because uniform sampling over a 2D section is swamped by
# long distances and leaves the short range - the range that matters - noisy.
CD_MAX_UM     <- 2000
CD_N_STRATA   <- 32
CD_PER_STRATUM <- 200000L

# Which column sets get a C(d) curve. Cost is linear in the number of sets, so
# a table with six subsets takes ~4x as long as the original two. "all" is the
# default because the subsets are the whole point of the current round; set it
# to "chrX,autosome" to get the old behaviour and the old runtime.
CD_SETS <- Sys.getenv("CD_SETS", "all")

# Angular sectors for C(d). 0 = one all-angles curve, which is the simple
# reading and the default. Set to 4 to also split C(d) by direction on the
# capture grid, which is worth doing if the isotropic curve looks odd:
# cardiomyocytes run in bundles, so clonal patches are probably elongated, and
# an elongated patch measured isotropically returns a length describing neither
# axis. Costs nothing but a second pass over the same sampled pairs.
CD_N_SECTORS <- 0L

# Optional: restrict to bins in one tissue domain, e.g. cardiomyocyte bins from
# the BANKSY clustering. A TSV with columns barcode, domain, plus the domain
# values to keep. This is the trick that buys the most - a tile restricted to
# one lineage stops mixing clonal histories, so you can push the tile size up
# for depth without paying the usual bias penalty. NULL disables it.
DOMAIN_TSV    <- NULL
DOMAIN_KEEP   <- NULL

# Set once you have looked at the plots, to write the barcode -> tile map that
# sinto consumes. NULL means "sweep only, tell me what you found".
# Also settable from the environment (TILE_UM=64 sbatch ...) so you do not have
# to edit and re-commit the script between the sweep and the sinto run.
TILE_UM_FOR_SINTO <- if (nzchar(Sys.getenv("TILE_UM"))) {
  as.integer(Sys.getenv("TILE_UM"))
} else {
  NULL
}
# Tiles below this many informative chrX UMIs are dropped from the sinto map -
# no point spending a BAM and an Allelome.PRO2 invocation on an empty tile.
# chrX is the limiting signal; autosomes have ~7x more, so a tile that clears
# this on chrX is comfortable on autosomes. Set to 0 to keep every tissue tile.
SINTO_MIN_UMI <- 10
##### ------------------------------------------------------ #####

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# fread() can only decompress .gz in-process via R.utils, which is not in the
# RNAseq env. Piping through gzip is equivalent and needs no extra package.
fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), ...)
  else fread(path, ...)
}

if (!file.exists(IN_TSV)) stop("No count table: ", IN_TSV)
bins <- fread_any(IN_TSV)
setnames(bins, tolower(names(bins)))
stopifnot(all(c("barcode", "array_row", "array_col",
                "x_ref", "x_alt", "a_ref", "a_alt") %in% names(bins)))

if (!is.null(DOMAIN_TSV)) {
  dom <- fread(DOMAIN_TSV)
  setnames(dom, tolower(names(dom)))
  keep <- dom[domain %in% DOMAIN_KEEP, barcode]
  n0 <- nrow(bins)
  bins <- bins[barcode %in% keep]
  message("Domain filter: kept ", nrow(bins), " / ", n0, " bins with data")
}

# Column sets. chrX and the autosomal control are always there; anything else
# the counter emitted as a <name>_ref / <name>_alt pair joins the analysis
# automatically, in the order the columns appear.
SET_COLS <- list(chrX = c("x_ref", "x_alt"), autosome = c("a_ref", "a_alt"))
extra <- setdiff(sub("_ref$", "", grep("_ref$", names(bins), value = TRUE)),
                 c("x", "a"))
for (nm in extra) {
  alt <- paste0(nm, "_alt")
  if (!alt %in% names(bins)) {
    warning("Column ", nm, "_ref has no matching ", alt, " - ignoring it")
    next
  }
  SET_COLS[[nm]] <- c(paste0(nm, "_ref"), alt)
}
if (length(extra)) message("Subset columns found: ", paste(extra, collapse = ", "))

bins[, x_n := x_ref + x_alt]
bins[, a_n := a_ref + a_alt]
if (nrow(bins) == 0L) {
  stop("No bins carry an informative UMI in ", IN_TSV, "\n",
       "  The counting pass produced nothing. Check its log for the ",
       "pre-flight line\n  (CB tags resolving to bins) and the allele-match ",
       "rate, then re-run with FORCE=1.")
}
message("Bins with any informative UMI: ", nrow(bins))
message("chrX informative UMIs: ", sum(bins$x_n),
        "   autosomal: ", sum(bins$a_n))

if (!file.exists(TISSUE_TSV)) stop("No tissue bin list: ", TISSUE_TSV)
tissue <- fread_any(TISSUE_TSV)
setnames(tissue, tolower(names(tissue)))
if (!is.null(DOMAIN_TSV)) {
  # The denominator has to be restricted the same way the numerator is, and
  # from the same source - a domain-restricted numerator over a whole-tissue
  # denominator would make every coverage number look worse than it is.
  tissue <- tissue[barcode %in% keep]
}
message("In-tissue 2um bins: ", nrow(tissue))

# Tiles that overlie tissue at each candidate size. This, not the number of
# tiles that happen to carry a read, is what the coverage fraction divides by.
n_tissue_tiles <- sapply(SIZES_UM, function(s) {
  k <- s / 2
  uniqueN(tissue[, .(array_row %/% k, array_col %/% k)])
})
names(n_tissue_tiles) <- as.character(SIZES_UM)

# ------------------------------------------------------------- dispersion

# Moment estimator. Var(p_hat) under a single shared p is p(1-p)/n, so whatever
# observed variance exceeds that, scaled by the maximum possible p(1-p), is the
# fraction of tiles that genuinely differ from each other.
rho_moment <- function(x, n) {
  p    <- x / n
  pbar <- sum(x) / sum(n)
  denom <- pbar * (1 - pbar)
  if (denom <= 0) return(NA_real_)
  (var(p) - mean(denom / n)) / denom
}

# Beta-binomial MLE. Same quantity, but it handles wildly uneven n properly
# instead of leaning on a mean-of-reciprocals, so it is the headline number.
# Parameterised directly in (mu, rho): Var(p_hat) = mu(1-mu)[1/n + (1-1/n)rho],
# so rho is the intraclass correlation and lives on (0, 1).
rho_betabinom <- function(x, n) {
  nll <- function(par) {
    mu  <- plogis(par[1])
    rho <- plogis(par[2])
    s   <- 1 / rho - 1
    a   <- mu * s
    b   <- (1 - mu) * s
    if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(1e10)
    -sum(lbeta(x + a, n - x + b) - lbeta(a, b))
  }
  fit <- tryCatch(
    optim(c(0, qlogis(0.05)), nll, method = "Nelder-Mead",
          control = list(maxit = 2000)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence != 0) return(NA_real_)
  plogis(fit$par[2])
}

# Median of `v` padded out to n_total with zeros, without building the padding:
# at 2um that would be a 4M-element vector per call, allocated 24 times.
median_with_zeros <- function(v, n_total) {
  n_total <- as.integer(n_total)
  if (n_total <= 0L) return(NA_real_)
  n_zero <- max(0L, n_total - length(v))
  s <- sort(v)
  at <- function(i) if (i <= n_zero) 0 else s[i - n_zero]
  if (n_total %% 2L == 1L) as.numeric(at((n_total + 1L) %/% 2L))
  else (at(n_total %/% 2L) + at(n_total %/% 2L + 1L)) / 2
}

sweep_one <- function(size_um, ref_col, alt_col, label) {
  k <- size_um / 2                       # 2um bins per tile side
  tmp <- data.table(tr = bins$array_row %/% k, tc = bins$array_col %/% k,
                    r  = bins[[ref_col]],      a  = bins[[alt_col]])
  tile <- tmp[, .(x = sum(r), n = sum(r) + sum(a)), by = .(tr, tc)]
  tile <- tile[n > 0]
  n_tissue <- n_tissue_tiles[[as.character(size_um)]]

  cov <- sapply(COVERAGE_TARGETS, function(t) sum(tile$n >= t) / n_tissue)
  names(cov) <- paste0("frac_ge_", COVERAGE_TARGETS)

  use <- tile[n >= MIN_UMI]
  rm_ <- if (nrow(use) >= 20) rho_moment(use$x, use$n) else NA_real_
  rb_ <- if (nrow(use) >= 20) rho_betabinom(use$x, use$n) else NA_real_

  res <- data.table(chrom_set = label, size_um = size_um,
                    n_tissue_tiles = n_tissue,
                    n_tiles_with_data = nrow(tile), n_tiles_used = nrow(use),
                    # Median over tissue tiles, so the empty ones count. The
                    # median over tiles-with-data is the flattering number and
                    # is not the one that decides whether a grid map is
                    # supportable.
                    median_umi = median_with_zeros(tile$n, n_tissue),
                    mean_umi = sum(tile$n) / n_tissue,
                    p_bar = sum(tile$x) / sum(tile$n),
                    rho_moment = rm_, rho_bb = rb_)
  # Assigned explicitly rather than splicing as.list(cov) into the constructor,
  # which is not guaranteed to become one column per element.
  for (nm in names(cov)) set(res, j = nm, value = cov[[nm]])
  res[]
}

message("\nSweeping tile sizes ...")
sweep <- rbindlist(lapply(names(SET_COLS), function(lab) {
  cols <- SET_COLS[[lab]]
  rbindlist(lapply(SIZES_UM, sweep_one, ref_col = cols[1], alt_col = cols[2],
                   label = lab))
}))
fwrite(sweep, file.path(OUT_DIR, "tile_size_sweep.csv"))
msg_table(sweep[chrom_set == "chrX",
                .(size_um, n_tissue_tiles, median_umi, frac_ge_10, frac_ge_20,
                  p_bar = round(p_bar, 4), rho_bb = round(rho_bb, 4))])
message("\nAutosomal control (same tiles, expected rho ~ 0):")
msg_table(sweep[chrom_set == "autosome",
                .(size_um, median_umi, p_bar = round(p_bar, 4),
                  rho_bb = round(rho_bb, 4))])

# ---------------------------------------------------- pair correlation C(d)

# For two bins with counts (r1,a1) and (r2,a2), the chance that one UMI drawn
# from each carries the same allele is (r1*r2 + a1*a2)/(n1*n2). Averaging that
# over many pairs at a given separation estimates C(d) with no binning and no
# within-bin resampling.
pair_correlation <- function(ref_col, alt_col, label) {
  d <- data.table(r = bins$array_row, c = bins$array_col,
                  ref = bins[[ref_col]], alt = bins[[alt_col]])
  d <- d[ref + alt > 0]
  d[, n := ref + alt]
  if (nrow(d) < 1000) return(NULL)
  # Integer key for O(1) partner lookup. Visium HD arrays are ~3350 bins per
  # side at 2um; 100000 leaves room to spare and keeps the key in a double.
  KEYMUL <- 100000
  d[, bkey := r * KEYMUL + c]
  setkey(d, bkey)

  dists <- unique(round(exp(seq(log(4), log(CD_MAX_UM), length.out = CD_N_STRATA))))
  out <- rbindlist(lapply(dists, function(dd) {  # dd is the nominal separation
    nb <- dd / 2                                   # separation in bin units
    th <- runif(CD_PER_STRATUM, 0, 2 * pi)
    dr <- as.integer(round(nb * cos(th)))
    dc <- as.integer(round(nb * sin(th)))
    # Integer rounding of the offset means the realised separation drifts from
    # the nominal one, worst at short range. Report what was actually sampled.
    ok <- !(dr == 0 & dc == 0)
    anchor <- d[sample.int(nrow(d), CD_PER_STRATUM, replace = TRUE)]
    q <- data.table(bkey = (anchor$r + dr) * KEYMUL + (anchor$c + dc),
                    r1 = anchor$ref, a1 = anchor$alt, n1 = anchor$n,
                    d_real = 2 * sqrt(dr^2 + dc^2),
                    # Fold to 0-180: an offset and its negation are the same
                    # separation, so opposite directions share a sector.
                    sector = as.integer(floor(
                      (atan2(dc, dr) %% pi) / pi * max(1L, CD_N_SECTORS))))[ok]
    hit <- d[q, on = "bkey", nomatch = 0]
    if (nrow(hit) < 100) return(NULL)
    hit[, agree := (r1 * ref + a1 * alt) / (n1 * n)]

    # sector -1 is the all-angles curve, kept so the isotropic reading stays
    # available and directly comparable to the per-sector ones.
    #
    # dist_nominal vs dist_um is the integer-rounding drift, and it is reported
    # rather than left implicit because it is worst at exactly the short range
    # the escape and imprinted comparisons depend on: at a nominal 4um the
    # offset is +-1 bin, so the realised separations are 2, 2.83 and 4um and
    # the mean lands near 3, not 4.
    iso <- hit[, .(sector = -1L, dist_nominal = dd, dist_um = mean(d_real),
                   d_sd = sd(d_real), d_min = min(d_real), d_max = max(d_real),
                   n_pairs = .N, C = mean(agree), se = sd(agree) / sqrt(.N))]
    per <- if (CD_N_SECTORS > 0L) {
      hit[, .(dist_nominal = dd, dist_um = mean(d_real), d_sd = sd(d_real),
              d_min = min(d_real), d_max = max(d_real),
              n_pairs = .N, C = mean(agree),
              se = sd(agree) / sqrt(.N)), by = sector][n_pairs >= 100]
    } else NULL
    rbindlist(list(iso, per), use.names = TRUE)[, chrom_set := label][]
  }))
  out
}

cd_sets <- if (CD_SETS == "all") names(SET_COLS) else
  intersect(strsplit(CD_SETS, ",")[[1]], names(SET_COLS))
message("Estimating pair correlation for: ", paste(cd_sets, collapse = ", "))
set.seed(1)
cd <- rbindlist(lapply(cd_sets, function(lab) {
  cols <- SET_COLS[[lab]]
  n_bins_set <- bins[get(cols[1]) + get(cols[2]) > 0, .N]
  n_umi_set  <- bins[, sum(get(cols[1]) + get(cols[2]))]
  message(sprintf("  %-12s %d bins carry a UMI, %d UMIs total%s",
                  lab, n_bins_set, n_umi_set,
                  if (n_bins_set < 1000) "  -- too few, skipped" else ""))
  r <- pair_correlation(cols[1], cols[2], lab)
  # A subset with too little depth returns NULL rather than a noisy curve. Say
  # so here: an absent curve in the PDF is otherwise indistinguishable from a
  # subset bed that matched nothing.
  if (is.null(r)) return(NULL)
  r[, `:=`(n_bins_set = n_bins_set, n_umi_set = n_umi_set)][]
}), use.names = TRUE)
fwrite(cd, file.path(OUT_DIR, "pair_correlation.csv"))
if (!nrow(cd)) stop("No set produced a pair-correlation curve")

# The short-range table, all sets together. This is the whole of tasks 3 and 4
# in one place: the imprinted sets give the per-UMI error floor, escape vs
# non-escape says whether the chrX floor is escape, and chrX and the autosomal
# control are the reference points.
#
# err is the per-UMI allele-call error rate implied by C, on the model that each
# UMI is called wrong independently with probability e:
#   C = (1-e)^2 + e^2  for a truly monoallelic locus
# solved for e, taking the root below 0.5. It is only meaningful for a set that
# SHOULD be monoallelic everywhere - the imprinted ones - which is why it is
# printed for all sets but interpreted only for those.
implied_err <- function(C) {
  ok <- is.finite(C) & C >= 0.5 & C <= 1
  out <- rep(NA_real_, length(C))
  out[ok] <- (1 - sqrt(2 * C[ok] - 1)) / 2
  out
}
short <- cd[sector == -1L][order(chrom_set, dist_um)][
  , .SD[1:min(3L, .N)], by = chrom_set][
  , .(dist_nominal = dist_nominal[1], dist_um = round(dist_um[1], 2),
      d_range = sprintf("%.1f-%.1f", d_min[1], d_max[1]),
      n_pairs = n_pairs[1], C = round(C[1], 4), se = round(se[1], 4),
      C_mean3 = round(mean(C), 4), implied_umi_err = round(implied_err(C[1]), 4)),
  by = chrom_set]
message("\nShort-range pair correlation, every set (realised separations, not nominal):")
msg_table(short)
fwrite(short, file.path(OUT_DIR, "pair_correlation_short_range.csv"))

imp <- short[grepl("^imp", chrom_set)]
if (nrow(imp)) {
  message("\n--- imprinted positive control ---")
  for (i in seq_len(nrow(imp))) {
    message(sprintf(
      "%s: C(%.1fum) = %.4f -> implied per-UMI error %.1f%%",
      imp$chrom_set[i], imp$dist_um[i], imp$C[i],
      100 * imp$implied_umi_err[i]))
  }
  cx0 <- short[chrom_set == "chrX", C]
  if (length(cx0) && all(is.finite(imp$C))) {
    if (max(imp$C) >= 0.97) {
      message("PASS: an imprinted locus is called consistently, so per-UMI ",
              "assignment is not the\n  ceiling, and the chrX floor of ",
              sprintf("%.3f", cx0), " is biology (escape or leak), not error.")
    } else if (max(imp$C) <= cx0 + 0.05) {
      message("FAIL: the imprinted control sits at ", sprintf("%.3f", max(imp$C)),
              ", no better than chrX's ", sprintf("%.3f", cx0), ".")
      message("  Per-UMI allele assignment, not biology, is the ceiling. The ",
              "tile analysis cannot\n  go further on this data and the ",
              "per-allele rate is measuring a noise rate.")
    } else {
      message("PARTIAL: the imprinted control is above chrX but below 0.97. ",
              "Some of the chrX\n  floor is error and some is biology; the ",
              "split is the gap between the two.")
    }
    message("  Before reading any of that: check the per-locus depth and ",
            "ratios in the counter's\n  --locus-out table. A set pooled over ",
            "loci that are not all monoallelic in heart\n  (Grb10 and Ube3a ",
            "are the usual offenders) fails this test for the wrong reason.")
  }
}
if (all(c("escape", "nonescape") %in% short$chrom_set)) {
  message("\n--- escape partition ---")
  message(sprintf("all chrX %.4f   core-escape only %.4f   chrX minus escape %.4f",
                  short[chrom_set == "chrX", C],
                  short[chrom_set == "escape", C],
                  short[chrom_set == "nonescape", C]))
  if (short[chrom_set == "nonescape", C] - short[chrom_set == "chrX", C] > 0.02) {
    message("Removing the escapees lifts the short range, so escape is part of ",
            "the floor.")
  } else {
    message("Removing the escapees changes nothing: the floor is not escape at ",
            "these 11 genes.")
  }
}

# C(d) = Cinf + (C0 - Cinf) * exp(-d / L). L is the patch length; tiles should
# sit at or below it, and below L/2 to resolve patch shape rather than just
# detect that patches exist.
fit_decay <- function(dt) {
  if (is.null(dt) || nrow(dt) < 8) return(NA_real_)
  st <- list(Cinf = min(dt$C), C0 = max(dt$C),
             L = max(10, dt$dist_um[which.min(abs(dt$C - mean(range(dt$C))))]))
  f <- tryCatch(nls(C ~ Cinf + (C0 - Cinf) * exp(-dist_um / L), data = dt,
                    start = st, control = nls.control(warnOnly = TRUE)),
                error = function(e) NULL)
  if (is.null(f)) return(NA_real_)
  L <- as.numeric(coef(f)["L"])
  if (!is.finite(L) || L <= 0 || L > CD_MAX_UM) return(NA_real_)
  L
}
# nls on a curve that barely decays returns a singular gradient rather than an
# answer, which is what happened to chrX. This needs no model: interpolate the
# distance at which C falls halfway from its short-range value to its plateau.
half_decay <- function(dt) {
  if (is.null(dt) || nrow(dt) < 6) return(NA_real_)
  dt <- dt[order(dist_um)]
  C0   <- mean(head(dt$C, 3))
  Cinf <- mean(tail(dt$C, 5))
  if (!is.finite(C0 - Cinf) || (C0 - Cinf) <= 0) return(NA_real_)
  half <- (C0 + Cinf) / 2
  i <- which(dt$C <= half)[1]
  if (is.na(i) || i == 1L) return(NA_real_)
  with(dt, dist_um[i - 1] + (half - C[i - 1]) *
           (dist_um[i] - dist_um[i - 1]) / (C[i] - C[i - 1]))
}

# The decisive comparison. Autosomes are diploid, so any short-range allelic
# correlation there is technical plus transcriptional bursting - one cell's
# burst spreads across the bins that cell covers. chrX carries that same floor
# PLUS clonal X-inactivation, so it is the excess over the autosomal curve, not
# the chrX curve itself, that measures patch structure.
cx <- cd[chrom_set == "chrX"    & sector == -1L][order(dist_um)]
ca <- cd[chrom_set == "autosome" & sector == -1L][order(dist_um)]
excess <- NULL
if (nrow(cx) > 1L && nrow(ca) > 1L) {   # approx() needs two points
  ca_at <- approx(ca$dist_um, ca$C, xout = cx$dist_um, rule = 2)$y
  excess <- data.table(dist_um = cx$dist_um,
                       C_chrX = cx$C, C_auto = ca_at,
                       excess = cx$C - ca_at,
                       se = cx$se)
  fwrite(excess, file.path(OUT_DIR, "pair_correlation_excess.csv"))
  message("\nchrX excess over the autosomal null, short range:")
  msg_table(excess[dist_um <= 200,
                   .(dist_um = round(dist_um), C_chrX = round(C_chrX, 4),
                     C_auto = round(C_auto, 4), excess = round(excess, 4),
                     se = round(se, 4))])
}

L_x_free <- half_decay(cx)
L_a_free <- half_decay(ca)
message(sprintf("\nHalf-decay (model-free): chrX %.0f um, autosome %.0f um",
                L_x_free, L_a_free))

L_x <- fit_decay(cd[chrom_set == "chrX" & sector == -1L])
L_a <- fit_decay(cd[chrom_set == "autosome" & sector == -1L])
message(sprintf("Correlation length (all angles): chrX %.0f um, autosome %.0f um",
                L_x, L_a))

# Per-sector decay lengths. The spread between them is the anisotropy, and the
# SHORT axis is what should set the tile size: a square tile large enough to
# span the long axis will already be mixing patches across the fibre.
sector_L <- NULL
if (CD_N_SECTORS > 0L) {
  secs <- sort(unique(cd[sector >= 0L, sector]))
  sector_L <- rbindlist(lapply(secs, function(k) {
    data.table(sector = k,
               angle_deg = round((k + 0.5) * 180 / CD_N_SECTORS),
               L_chrX = fit_decay(cd[chrom_set == "chrX" & sector == k]),
               L_auto = fit_decay(cd[chrom_set == "autosome" & sector == k]))
  }))
  fwrite(sector_L, file.path(OUT_DIR, "pair_correlation_by_sector.csv"))
  message("Directional correlation length (chrX):")
  msg_table(sector_L)
}

# ------------------------------------------------------------------ figures

pdf(file.path(OUT_DIR, "tile_size_sweep.pdf"), width = 7.5, height = 5.5)

# SAMPLE goes in every title. Both samples write a file of this name into their
# own directory, and the two PDFs were previously distinguishable only by parent
# directory - which has already caused one mix-up.
titled <- function(s) paste0(SAMPLE, " - ", s)

cov_long <- melt(sweep, id.vars = c("chrom_set", "size_um"),
                 measure.vars = paste0("frac_ge_", COVERAGE_TARGETS),
                 variable.name = "target", value.name = "frac")
cov_long[, target := sub("frac_ge_", ">= ", target)]
print(
  ggplot(cov_long[chrom_set == "chrX"], aes(size_um, frac, colour = target)) +
    geom_line() + geom_point() +
    scale_x_log10(breaks = SIZES_UM) +
    labs(title = titled("Coverage: tiles reaching a usable depth"),
         subtitle = "chrX informative UMIs per tile",
         x = "tile size (um)", y = "fraction of tissue tiles",
         colour = "UMIs") +
    theme_bw()
)

print(
  ggplot(sweep, aes(size_um, rho_bb, colour = chrom_set)) +
    geom_line() + geom_point() +
    scale_x_log10(breaks = SIZES_UM) +
    labs(title = titled("Mosaic dispersion rho(s)"),
         subtitle = paste("Beta-binomial intraclass correlation.",
                          "chrX above the autosomal null is real structure;",
                          "the fall-off marks where tiles start averaging",
                          "over more than one patch."),
         x = "tile size (um)", y = "rho") +
    theme_bw() + theme(plot.subtitle = element_text(size = 7))
)

if (nrow(cd)) {
  print(
    ggplot(cd[sector == -1L & chrom_set %in% c("chrX", "autosome")],
           aes(dist_um, C, colour = chrom_set)) +
      geom_ribbon(aes(ymin = C - 2 * se, ymax = C + 2 * se, fill = chrom_set),
                  alpha = 0.15, colour = NA) +
      geom_line() +
      scale_x_log10() +
      labs(title = titled("Pair correlation C(d) - no tiling involved"),
           subtitle = sprintf("Decay length: chrX %.0f um, autosome %.0f um",
                              L_x, L_a),
           x = "separation (um)",
           y = "P(two UMIs share an allele)") +
      theme_bw()
  )

  # Every set on one axis, with the two reference lines that make it readable:
  # 1 is what a monoallelic locus should give, 0.5 what an unstructured
  # biallelic one gives. The imprinted curves belong up against the top; how far
  # short they fall IS the per-UMI error floor, and the same distance is
  # subtracted from every chrX number in this project.
  if (uniqueN(cd$chrom_set) > 2L) {
    print(
      ggplot(cd[sector == -1L], aes(dist_um, C, colour = chrom_set)) +
        geom_hline(yintercept = c(0.5, 1), linetype = 2, colour = "grey60") +
        geom_ribbon(aes(ymin = C - 2 * se, ymax = C + 2 * se, fill = chrom_set),
                    alpha = 0.12, colour = NA) +
        geom_line() +
        scale_x_log10() +
        # coord_cartesian, not scale limits: a curve that leaves the panel
        # should be clipped, not dropped with a warning about missing values.
        coord_cartesian(ylim = c(0.4, 1.02)) +
        labs(title = titled("C(d) by SNP set: the floor, partitioned"),
             subtitle = paste("imp* are imprinted positive controls - monoallelic",
                              "in every cell, so they should sit flat near 1 and",
                              "the shortfall is per-UMI error.\nescape/nonescape",
                              "partition chrX; noXic is chrX with Xist +- 500 kb",
                              "masked. Sets are counted on the same UMIs."),
             x = "separation (um)", y = "P(two UMIs share an allele)",
             colour = "SNP set", fill = "SNP set") +
        theme_bw() + theme(plot.subtitle = element_text(size = 7))
    )
  }

  if (!is.null(excess)) {
    print(
      ggplot(excess, aes(dist_um, excess)) +
        geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
        geom_ribbon(aes(ymin = excess - 2 * se, ymax = excess + 2 * se),
                    alpha = 0.15) +
        geom_line() +
        scale_x_log10() +
        labs(title = titled("chrX clonal signal: C(d) above the autosomal null"),
             subtitle = paste("Autosomes carry the same cell-footprint and",
                              "bursting correlation without clonal XCI, so",
                              "the excess is the mosaicism. Flat at zero",
                              "means no patch structure at this depth."),
             x = "separation (um)", y = "C_chrX(d) - C_autosome(d)") +
        theme_bw() + theme(plot.subtitle = element_text(size = 7))
    )
  }

  if (CD_N_SECTORS > 0L && nrow(cd[sector >= 0L])) {
    print(
      ggplot(cd[sector >= 0L],
             aes(dist_um, C, colour = factor(round((sector + 0.5) * 180 /
                                                   CD_N_SECTORS)))) +
        geom_line() +
        facet_wrap(~ chrom_set) +
        scale_x_log10() +
        labs(title = titled("C(d) by direction on the capture grid"),
             subtitle = paste("Sectors that separate on chrX but not on the",
                              "autosomal control mean the patches are",
                              "elongated - set the tile size from the",
                              "shortest axis."),
             x = "separation (um)", y = "P(two UMIs share an allele)",
             colour = "angle (deg)") +
        theme_bw() + theme(plot.subtitle = element_text(size = 7))
    )
  }
}

# The intuitive version of rho(s): at a well-chosen size the tile ratios are
# dispersed or frankly bimodal; too large and they collapse onto a single spike.
for (s in intersect(c(16, 32, 64, 128, 256), SIZES_UM)) {
  k <- s / 2
  tl <- bins[, .(x = sum(x_ref), n = sum(x_ref) + sum(x_alt)),
             by = .(array_row %/% k, array_col %/% k)][n >= MIN_UMI]
  if (nrow(tl) < 50) next
  print(
    ggplot(tl, aes(x / n)) +
      geom_histogram(bins = 40, fill = "grey30") +
      labs(title = titled(paste0("chrX allelic ratio per tile, ", s, " um")),
           subtitle = paste0(nrow(tl), " of ",
                             n_tissue_tiles[[as.character(s)]],
                             " tissue tiles clear ", MIN_UMI,
                             " UMIs; median ", median(tl$n), " UMIs among those"),
           x = "B6 fraction", y = "tiles") +
      theme_bw()
  )
}
invisible(dev.off())

# --------------------------------------------------------- recommendation

x <- sweep[chrom_set == "chrX"]
best_rho <- if (all(is.na(x$rho_bb))) NULL else x[which.max(rho_bb)]
usable   <- x[frac_ge_10 >= 0.5]
message("\n--- recommendation ---")
if (!is.null(excess)) {
  e0 <- mean(head(excess$excess, 3))
  s0 <- mean(head(excess$se, 3))
  if (is.finite(e0) && e0 > 2 * s0) {
    message(sprintf("chrX carries clonal signal: excess %.4f at short range (%.1f SE).",
                    e0, e0 / s0))
    if (!is.na(L_x_free)) {
      message(sprintf("  Patch half-decay ~%.0f um -> tile at <= %.0f um.",
                      L_x_free, L_x_free / 2))
    }
  } else {
    message(sprintf("No clonal signal: chrX C(d) sits %.4f above the autosomal null at short range (%.1f SE).",
                    e0, if (s0 > 0) e0 / s0 else NA_real_))
    message("  Either patches are smaller than the depth allows resolving, or")
    message("  diffusion has smeared them. A grid map has nothing to show; the")
    message("  per-domain pseudobulk is where the signal will be.")
  }
}
if (!is.na(L_x)) {
  message(sprintf("C(d) puts the chrX patch length at ~%.0f um (autosomal null %.0f um).",
                  L_x, L_a))
  message(sprintf("  -> to resolve patch shape, tile at <= %.0f um.", L_x / 2))
}
if (!is.null(sector_L) && sum(!is.na(sector_L$L_chrX)) >= 2) {
  Lmin <- min(sector_L$L_chrX, na.rm = TRUE)
  Lmax <- max(sector_L$L_chrX, na.rm = TRUE)
  message(sprintf("Directional: %.0f um across the short axis, %.0f um along the long one (%.1fx).",
                  Lmin, Lmax, Lmax / Lmin))
  if (Lmax / Lmin >= 1.5) {
    message(sprintf("  Patches are elongated, so the short axis governs: tile at <= %.0f um.",
                    Lmin / 2))
    message("  Square tiles are a fair first pass, but a tile sized for the ",
            "long axis will be mixing patches across the fibre.")
  }
}
if (is.null(best_rho)) {
  message("rho could not be estimated at any tile size - too few tiles clear ",
          "MIN_UMI. That is itself the answer: the grid is depth-limited.")
} else {
  message(sprintf("rho peaks at %d um (rho = %.3f).",
                  best_rho$size_um, best_rho$rho_bb))
}
if (nrow(usable)) {
  message(sprintf("Smallest tile where >=50%% of tiles clear 10 UMIs: %d um.",
                  min(usable$size_um)))
} else {
  message("No tile size gets 50% of tiles to 10 chrX UMIs. The grid map is not ",
          "supportable at this depth - use the domain pseudobulk instead, or ",
          "restrict to cardiomyocyte bins via DOMAIN_TSV and re-run.")
}
message("Set TILE_UM_FOR_SINTO once you have looked at ",
        file.path(OUT_DIR, "tile_size_sweep.pdf"))

# ------------------------------------------------- barcode -> tile for sinto

if (!is.null(TILE_UM_FOR_SINTO)) {
  k <- TILE_UM_FOR_SINTO / 2
  tile_id <- function(r, cc) sprintf("tile_%dum_r%04d_c%04d",
                                     as.integer(TILE_UM_FOR_SINTO),
                                     as.integer(r %/% k), as.integer(cc %/% k))

  # Depth comes from `bins`, which only holds bins carrying an informative UMI.
  bins[, tile := tile_id(array_row, array_col)]
  depth <- bins[, .(n_chrX = sum(x_n), n_auto = sum(a_n)), by = tile]
  keep  <- depth[n_chrX >= SINTO_MIN_UMI, tile]

  # ...but the MAP must come from every tissue bin. Four fifths of tissue bins
  # carry no informative UMI, and their reads still belong in the tile: they
  # cover SNPs this counter skipped, they carry autosomal signal, and they make
  # up the denominator Allelome.PRO2 works from. Building the map out of `bins`
  # would hand sinto a fifth of the data.
  tissue[, tile := tile_id(array_row, array_col)]
  map <- tissue[tile %in% keep, .(barcode, tile)]

  message(sprintf("\nsinto map at %dum: %d tiles, %d of %d tissue bins (%.0f%%)",
                  TILE_UM_FOR_SINTO, length(keep), nrow(map), nrow(tissue),
                  100 * nrow(map) / nrow(tissue)))
  message(sprintf("  bins per tile: %.0f of a possible %d",
                  nrow(map) / max(1L, length(keep)), as.integer(k)^2))
  message(sprintf("  informative UMIs per tile, median: chrX %d, autosomal %d",
                  as.integer(median(depth[tile %in% keep, n_chrX])),
                  as.integer(median(depth[tile %in% keep, n_auto]))))
  if (length(keep) > 5000) {
    warning("More than 5000 tiles - sinto holds one output BAM open per tile ",
            "and will run into the open-file limit. Raise TILE_UM_FOR_SINTO ",
            "or SINTO_MIN_UMI.")
  }

  fwrite(map, file.path(OUT_DIR, sprintf("sinto_tiles_%dum.tsv",
                                         TILE_UM_FOR_SINTO)),
         sep = "\t", col.names = FALSE)
  fwrite(depth[tile %in% keep][order(-n_chrX)],
         file.path(OUT_DIR, sprintf("tile_depth_%dum.csv", TILE_UM_FOR_SINTO)))
}

# ------------------------------------------------------------- provenance
#
# Which SNP bed, which count table, which thresholds. The count table's own
# sidecar (written by ase_bin_allele_counts.py) carries the bed path and md5, so
# it is folded in here rather than restated - if it is missing, that is itself
# worth recording, because it means the counts predate the sidecar and their bed
# cannot be identified.
# NOT data.table(key = ..., value = ...): `key` is data.table()'s own argument,
# so that form sets a key instead of making a column, and errors.
prov <- data.table(
  k = c("script", "sample", "run_at", "counts_tsv", "counts_tsv_md5",
          "tissue_tsv", "min_umi", "coverage_targets", "sizes_um",
          "cd_max_um", "cd_n_strata", "cd_per_stratum", "cd_sets",
          "cd_n_sectors", "domain_tsv", "tile_um_for_sinto", "sinto_min_umi",
          "sets", "r_version"),
  v = c("spatial/ase_tile_sweep.R", SAMPLE,
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            IN_TSV, unname(tools::md5sum(IN_TSV)), TISSUE_TSV,
            MIN_UMI, paste(COVERAGE_TARGETS, collapse = ","),
            paste(SIZES_UM, collapse = ","), CD_MAX_UM, CD_N_STRATA,
            CD_PER_STRATUM, paste(cd_sets, collapse = ","), CD_N_SECTORS,
            if (is.null(DOMAIN_TSV)) "NULL" else DOMAIN_TSV,
            if (is.null(TILE_UM_FOR_SINTO)) "NULL" else TILE_UM_FOR_SINTO,
            SINTO_MIN_UMI, paste(names(SET_COLS), collapse = ","),
            paste(R.version$major, R.version$minor, sep = ".")))
setnames(prov, c("key", "value"))
counter_prov <- paste0(IN_TSV, ".provenance.tsv")
if (file.exists(counter_prov)) {
  cp <- fread(counter_prov, colClasses = "character")
  prov <- rbindlist(list(prov, cp[, .(key = paste0("counter.", key), value)]))
} else {
  miss <- data.table(k = "counter.provenance",
                     v = paste0("MISSING: ", counter_prov,
                       " - the SNP bed behind these counts cannot be identified"))
  setnames(miss, c("key", "value"))
  prov <- rbindlist(list(prov, miss))
}
fwrite(prov, file.path(OUT_DIR, "provenance_tile_sweep.tsv"), sep = "\t")
message("Provenance in ", file.path(OUT_DIR, "provenance_tile_sweep.tsv"))

message("\nDone. Outputs in ", OUT_DIR)
