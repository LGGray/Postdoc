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
# Run on the cluster:  sbatch ~/Postdoc/slurm/spatial_ase_sweep.slurm

.libPaths(c("~/R/matrix-dev", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

##### ---------------------- CONFIG ---------------------- #####
BASE      <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial"
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
SINTO_MIN_UMI <- 10
##### ------------------------------------------------------ #####

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(IN_TSV)) stop("No count table: ", IN_TSV)
bins <- fread(IN_TSV)
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

bins[, x_n := x_ref + x_alt]
bins[, a_n := a_ref + a_alt]
message("Bins with any informative UMI: ", nrow(bins))
message("chrX informative UMIs: ", sum(bins$x_n),
        "   autosomal: ", sum(bins$a_n))

if (!file.exists(TISSUE_TSV)) stop("No tissue bin list: ", TISSUE_TSV)
tissue <- fread(TISSUE_TSV)
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

  data.table(chrom_set = label, size_um = size_um,
             n_tissue_tiles = n_tissue,
             n_tiles_with_data = nrow(tile), n_tiles_used = nrow(use),
             # Median over tissue tiles, so the empty ones count. The median
             # over tiles-with-data is the flattering number and is not the
             # one that decides whether a grid map is supportable.
             median_umi = as.numeric(median(c(tile$n, rep(0, max(0, n_tissue - nrow(tile)))))),
             mean_umi = sum(tile$n) / n_tissue,
             p_bar = sum(tile$x) / sum(tile$n),
             rho_moment = rm_, rho_bb = rb_,
             as.list(cov))
}

message("\nSweeping tile sizes ...")
sweep <- rbindlist(c(
  lapply(SIZES_UM, sweep_one, ref_col = "x_ref", alt_col = "x_alt", label = "chrX"),
  lapply(SIZES_UM, sweep_one, ref_col = "a_ref", alt_col = "a_alt", label = "autosome")
))
fwrite(sweep, file.path(OUT_DIR, "tile_size_sweep.csv"))
print(sweep[chrom_set == "chrX",
            .(size_um, n_tissue_tiles, median_umi, frac_ge_10, frac_ge_20, rho_bb)])

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
    anchor <- d[sample.int(nrow(d), CD_PER_STRATUM, replace = TRUE)]
    # Integer rounding of the offset means the realised separation drifts from
    # the nominal one, worst at short range. Report what was actually sampled.
    ok <- !(dr == 0 & dc == 0)
    q <- data.table(bkey = (anchor$r + dr) * KEYMUL + (anchor$c + dc),
                    r1 = anchor$ref, a1 = anchor$alt, n1 = anchor$n,
                    d_real = 2 * sqrt(dr^2 + dc^2))[ok]
    hit <- d[q, on = "bkey", nomatch = 0]
    if (nrow(hit) < 100) return(NULL)
    agree <- (hit$r1 * hit$ref + hit$a1 * hit$alt) / (hit$n1 * hit$n)
    data.table(chrom_set = label, dist_um = mean(hit$d_real),
               n_pairs = nrow(hit), C = mean(agree),
               se = sd(agree) / sqrt(nrow(hit)))
  }))
  out
}

message("Estimating pair correlation ...")
set.seed(1)
cd <- rbindlist(list(
  pair_correlation("x_ref", "x_alt", "chrX"),
  pair_correlation("a_ref", "a_alt", "autosome")
))
fwrite(cd, file.path(OUT_DIR, "pair_correlation.csv"))

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
L_x <- fit_decay(cd[chrom_set == "chrX"])
L_a <- fit_decay(cd[chrom_set == "autosome"])
message(sprintf("Correlation length: chrX %.0f um, autosome %.0f um",
                L_x, L_a))

# ------------------------------------------------------------------ figures

pdf(file.path(OUT_DIR, "tile_size_sweep.pdf"), width = 7.5, height = 5.5)

cov_long <- melt(sweep, id.vars = c("chrom_set", "size_um"),
                 measure.vars = paste0("frac_ge_", COVERAGE_TARGETS),
                 variable.name = "target", value.name = "frac")
cov_long[, target := sub("frac_ge_", ">= ", target)]
print(
  ggplot(cov_long[chrom_set == "chrX"], aes(size_um, frac, colour = target)) +
    geom_line() + geom_point() +
    scale_x_log10(breaks = SIZES_UM) +
    labs(title = "Coverage: tiles reaching a usable depth",
         subtitle = "chrX informative UMIs per tile",
         x = "tile size (um)", y = "fraction of tissue tiles",
         colour = "UMIs") +
    theme_bw()
)

print(
  ggplot(sweep, aes(size_um, rho_bb, colour = chrom_set)) +
    geom_line() + geom_point() +
    scale_x_log10(breaks = SIZES_UM) +
    labs(title = "Mosaic dispersion rho(s)",
         subtitle = paste("Beta-binomial intraclass correlation.",
                          "chrX above the autosomal null is real structure;",
                          "the fall-off marks where tiles start averaging",
                          "over more than one patch."),
         x = "tile size (um)", y = "rho") +
    theme_bw() + theme(plot.subtitle = element_text(size = 7))
)

if (nrow(cd)) {
  print(
    ggplot(cd, aes(dist_um, C, colour = chrom_set)) +
      geom_ribbon(aes(ymin = C - 2 * se, ymax = C + 2 * se, fill = chrom_set),
                  alpha = 0.15, colour = NA) +
      geom_line() +
      scale_x_log10() +
      labs(title = "Pair correlation C(d) - no tiling involved",
           subtitle = sprintf("Decay length: chrX %.0f um, autosome %.0f um",
                              L_x, L_a),
           x = "separation (um)",
           y = "P(two UMIs share an allele)") +
      theme_bw()
  )
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
      labs(title = paste0("chrX allelic ratio per tile, ", s, " um"),
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
if (!is.na(L_x)) {
  message(sprintf("C(d) puts the chrX patch length at ~%.0f um (autosomal null %.0f um).",
                  L_x, L_a))
  message(sprintf("  -> to resolve patch shape, tile at <= %.0f um.", L_x / 2))
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
  bins[, tile := sprintf("tile_%dum_r%04d_c%04d", as.integer(TILE_UM_FOR_SINTO),
                         as.integer(array_row %/% k), as.integer(array_col %/% k))]
  depth <- bins[, .(n = sum(x_n)), by = tile]
  keep  <- depth[n >= SINTO_MIN_UMI, tile]
  map   <- bins[tile %in% keep, .(barcode, tile)]

  # sinto filterbarcodes holds one output BAM open per group. A few thousand is
  # already heavy on file descriptors and I/O; tens of thousands will not run.
  message(sprintf("\nsinto map: %d tiles, %d bins, median %d UMIs/tile",
                  length(keep), nrow(map),
                  as.integer(median(depth[tile %in% keep, n]))))
  if (length(keep) > 5000) {
    warning("More than 5000 tiles - sinto will struggle with the open file ",
            "handles. Raise TILE_UM_FOR_SINTO or SINTO_MIN_UMI.")
  }
  fwrite(map, file.path(OUT_DIR, sprintf("sinto_tiles_%dum.tsv",
                                         TILE_UM_FOR_SINTO)),
         sep = "\t", col.names = FALSE)
  fwrite(depth[tile %in% keep][order(-n)],
         file.path(OUT_DIR, sprintf("tile_depth_%dum.csv", TILE_UM_FOR_SINTO)))
}

message("\nDone. Outputs in ", OUT_DIR)
