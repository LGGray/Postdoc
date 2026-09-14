# Quantify spatial structure in the per-tile chrX allelic ratio.
#
# WHAT THIS ASKS, AND WHAT IT CANNOT ASK.
#
# The CAST X is the inactive X in EVERY cell of this cross (B6 mother x CAST
# father, Xist deleted on the B6 X). There is no mosaic. So this script is not
# looking for XCI patches - there are none to find, and ase_tile_sweep.R:777
# already says so on the rho(s) figure. The question here is narrower and is
# the only one the design supports:
#
#     does the LEVEL of escape vary systematically from place to place?
#
# A positive answer is not clonal patchiness. The candidate explanations, in
# the order they should be ruled out, are:
#   1. depth       - coverage is strongly spatially structured (Moran's I on
#                    log depth is 0.35-0.63 on this data) and the ratio is
#                    depth-biased, so depth alone can induce ratio structure.
#   2. composition - cell types differ in escape and are spatially organised.
#                    THIS IS NOT YET TESTED HERE; see the note at the end.
#   3. tissue architecture / sectioning artefacts.
#   4. genuine regional differences in escape.
#
# The autosomal ratio is carried through every statistic as the negative
# control. It has no monoallelic biology, so any structure it shows is
# technical and invalidates the chrX statement over the same tiles.
#
# WHY MORAN'S I AND NOT rho / phi. NEXT_ANALYSIS.md's "Do not do" list already
# rules out reading a patch size off rho(s) (monotone, no plateau) and reporting
# phi from the scase table (unidentifiable at 16um). Both are variance-partition
# statistics: they say how much tiles differ, not whether nearby tiles differ
# LESS than distant ones. Only the second is a spatial question, and Moran's I
# is the direct estimator of it.
#
# Moran's I has one further advantage on a duplicate-inclusive tree. phi and rho
# are measured against a binomial n, and with PCR duplicates kept that n is
# inflated by the duplication factor (2.5-3.8x here, and not uniform across the
# slide), so both are inflated by an amount that is not biology. Moran's I is
# computed on the ratio VALUES and does not reference n at all. Duplication
# makes each tile noisier than its n implies, which DEPRESSES I - so an I that
# is significantly positive on this tree is conservative, not inflated.
#
# USAGE (cluster; seurat_env, NOT the repo-wide RNAseq env)
#   conda activate seurat_env
#   Rscript ~/Postdoc/spatial/tile_spatial_structure.R
#   TREE=ase_pysam_dup_64um COUNT_UNIT=read MIN_X=20 Rscript ...
#
# Every knob is an environment variable with the default stated beside it, and
# the slurm wrapper passes them positionally because of --export=NONE.

.libPaths(c("~/R/matrix-dev", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

##### ---------------------------- CONFIG ---------------------------- #####
BASE   <- Sys.getenv("BASE",
  "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial")
# Which counting tree. The duplicate-inclusive pysam one is the default because
# it is the one that actually contains duplicates: Allelome.PRO2 cannot see them
# (Allelome.PRO2.sh:289 runs mpileup with no --ff, so mpileup's default
# UNMAP,SECONDARY,QCFAIL,DUP mask drops them before scoring, and the AP2 _dup
# tree is byte-identical to the dedup one).
TREE    <- Sys.getenv("TREE", "ase_pysam_dup_64um")
SAMPLES <- strsplit(Sys.getenv("SAMPLES", "9w,78w"), ",")[[1]]
# read or umi. The tile_chrom_counts.tsv table carries both.
COUNT_UNIT <- Sys.getenv("COUNT_UNIT", "read")
# Minimum informative chrX units per tile. A tile below this carries no ratio
# worth correlating and contributes only noise, which depresses I.
MIN_X   <- as.integer(Sys.getenv("MIN_X", "20"))
# Reads per informative molecule, for n_eff. 1 for a umi run. Only the
# beta-binomial null uses it; Moran's I does not reference n.
DUP_FACTOR <- as.numeric(Sys.getenv("DUP_FACTOR",
                                    if (COUNT_UNIT == "read") "3" else "1"))
NPERM   <- as.integer(Sys.getenv("NPERM", "999"))
# Correlogram: Moran's I at Chebyshev ring distance 1..MAX_LAG. Tile side in um
# converts a lag to a physical distance; it is NOT read from the tree name
# because a tree can be re-tiled without being renamed.
TILE_UM   <- as.numeric(Sys.getenv("TILE_UM", "64"))
MAX_LAG   <- as.integer(Sys.getenv("MAX_LAG", "12"))
# Fewer permutations per lag than for the global statistic: MAX_LAG x 3 fields
# of them, and a correlogram is read as a shape, not off one lag's p-value.
NPERM_LAG <- as.integer(Sys.getenv("NPERM_LAG", "199"))
SEED    <- as.integer(Sys.getenv("SEED", "1"))
OUT_DIR <- Sys.getenv("OUT_DIR", file.path(BASE, TREE))
TAG     <- Sys.getenv("TAG", paste0(COUNT_UNIT, "_min", MIN_X))

set.seed(SEED)
AUTOSOMES <- paste0("chr", 1:19)
msg <- function(...) message(sprintf(...))

##### ------------------------- SPATIAL CORE ------------------------- #####
# Queen adjacency on the integer tile grid. The grid is the right object here,
# not pixel centroids: tiles are a regular lattice by construction, and an
# integer neighbour test cannot drift the way a distance threshold on centroids
# can when a tile is partly off-tissue.
build_edges <- function(trow, tcol) {
  key <- paste(trow, tcol, sep = "_")
  idx <- setNames(seq_along(key), key)
  offs <- expand.grid(dr = -1:1, dc = -1:1)
  offs <- offs[!(offs$dr == 0 & offs$dc == 0), ]
  ii <- integer(0); jj <- integer(0)
  for (k in seq_len(nrow(offs))) {
    j <- idx[paste(trow + offs$dr[k], tcol + offs$dc[k], sep = "_")]
    ok <- !is.na(j)
    ii <- c(ii, which(ok)); jj <- c(jj, as.integer(j[ok]))
  }
  list(i = ii, j = jj, W = length(ii))
}

# Edges between tiles at Chebyshev ring distance exactly d, i.e. the ring of
# tiles d steps out. d = 1 reproduces build_edges(). Rings rather than cumulative
# discs, because a correlogram has to show where I FALLS TO the null, and a
# cumulative statistic cannot fall - it only dilutes.
ring_edges <- function(trow, tcol, d) {
  key  <- as.numeric(trow) * 1e5 + as.numeric(tcol)
  ord  <- order(key); skey <- key[ord]
  offs <- expand.grid(dr = -d:d, dc = -d:d)
  offs <- offs[pmax(abs(offs$dr), abs(offs$dc)) == d, , drop = FALSE]
  ii <- integer(0); jj <- integer(0)
  for (k in seq_len(nrow(offs))) {
    m  <- match(as.numeric(trow + offs$dr[k]) * 1e5 +
                as.numeric(tcol + offs$dc[k]), skey)
    ok <- !is.na(m)
    ii <- c(ii, which(ok)); jj <- c(jj, ord[m[ok]])
  }
  list(i = ii, j = jj, W = length(ii))
}

# Global Moran's I with a permutation null. The permutation is the honest null
# here rather than the normal approximation: the tile ratios are bounded,
# heavily skewed towards 1 and have wildly uneven precision, so the analytic
# variance of I is not trustworthy on them.
morans_I <- function(e, v, nperm = NPERM) {
  n <- length(v); dev <- v - mean(v); den <- sum(dev^2)
  if (den == 0 || e$W == 0) return(list(I = NA_real_, p = NA_real_, z = NA_real_))
  I <- (n / e$W) * sum(dev[e$i] * dev[e$j]) / den
  perm <- vapply(seq_len(nperm), function(k) {
    d <- dev[sample.int(n)]
    (n / e$W) * sum(d[e$i] * d[e$j]) / sum(d^2)
  }, numeric(1))
  list(I = I, p = (sum(perm >= I) + 1) / (nperm + 1),
       z = (I - mean(perm)) / stats::sd(perm))
}

# Local Moran (LISA), conditional permutation: tile i is held fixed and the
# other n-1 values are shuffled, which is the null "this tile's neighbourhood is
# no more similar to it than a random neighbourhood would be".
local_moran <- function(e, v, nperm = 499) {
  n <- length(v); dev <- v - mean(v); m2 <- sum(dev^2) / n
  nbsum <- vapply(seq_len(n), function(i) sum(dev[e$j[e$i == i]]), numeric(1))
  Ii <- dev * nbsum / m2
  deg <- tabulate(e$i, nbins = n)
  cnt <- integer(n)
  for (k in seq_len(nperm)) {
    # A permuted neighbour-sum of the right SIZE for each tile. A contiguous run
    # of the shuffled vector is exchangeable with deg[i] values drawn at random,
    # and is O(n) instead of O(sum(deg)). The vector is doubled so a run that
    # starts near the end WRAPS rather than being truncated - truncating would
    # shorten the run, shrink |sim| and make every p-value optimistic.
    d  <- dev[sample.int(n)]
    cs <- cumsum(c(0, d, d))
    st <- sample.int(n, n, replace = TRUE)
    sim <- dev * (cs[st + deg] - cs[st]) / m2
    cnt <- cnt + (abs(sim) >= abs(Ii))
  }
  p <- (cnt + 1) / (nperm + 1)
  data.table(lisa = Ii, lisa_p = p, lisa_q = p.adjust(p, "BH"), n_nb = deg)
}

##### ---------------------------- LOAD ------------------------------ #####
# tile_chrom_counts.tsv names the read columns a1_reads/a2_reads and the
# molecule ones a1_umi/a2_umi - plural on one, singular on the other. COUNT_UNIT
# stays singular ("read"/"umi") to match every other script in spatial/, so the
# mapping is explicit here rather than pasted.
UNIT_COL <- c(read = "reads", umi = "umi")[COUNT_UNIT]
if (is.na(UNIT_COL)) stop("COUNT_UNIT must be 'read' or 'umi', got '", COUNT_UNIT, "'")
a1c <- paste0("a1_", UNIT_COL); a2c <- paste0("a2_", UNIT_COL)
load_sample <- function(smp) {
  f <- file.path(BASE, TREE, smp, "tile_chrom_counts.tsv")
  if (!file.exists(f)) { msg("  no tile_chrom_counts.tsv for %s - skipped", smp); return(NULL) }
  d <- fread(f)
  for (col in c(a1c, a2c, "trow", "tcol", "chrom"))
    if (!col %in% names(d)) stop("column ", col, " missing from ", f)
  d[, arm := fifelse(chrom == "chrX", "x", fifelse(chrom %in% AUTOSOMES, "a", NA_character_))]
  d <- d[!is.na(arm)]
  w <- dcast(d[, .(a1 = sum(get(a1c)), n = sum(get(a1c)) + sum(get(a2c))),
               by = .(trow, tcol, arm)],
             trow + tcol ~ arm, value.var = c("a1", "n"), fill = 0)
  setnames(w, c("a1_x", "n_x", "a1_a", "n_a"), c("x_a1", "x_n", "a_a1", "a_n"),
           skip_absent = TRUE)
  w[, sample := smp]
  w[x_n >= MIN_X & a_n > 0]
}

##### ---------------------------- RUN ------------------------------- #####
summ <- list(); per_tile <- list(); correlo <- list()
for (smp in SAMPLES) {
  d <- load_sample(smp); if (is.null(d) || !nrow(d)) next
  msg("\n=== %s : %d tiles with >= %d chrX %ss ===", smp, nrow(d), MIN_X, COUNT_UNIT)
  d[, `:=`(x_ratio = x_a1 / x_n, a_ratio = a_a1 / a_n)]
  e <- build_edges(d$trow, d$tcol)
  msg("  %d tiles, %d adjacency edges (mean %.1f neighbours)",
      nrow(d), e$W, e$W / nrow(d))

  # The four global statistics. Depth is logged because it spans orders of
  # magnitude and Moran's I on the raw counts would be driven by a few tiles.
  vars <- list(chrX_ratio      = d$x_ratio,
               chrX_log_depth  = log(d$x_n),
               auto_ratio      = d$a_ratio,
               auto_log_depth  = log(d$a_n))
  for (nm in names(vars)) {
    r <- morans_I(e, vars[[nm]])
    msg("  Moran's I  %-16s %+.4f   perm p = %.3f   z = %+.1f", nm, r$I, r$p, r$z)
    summ[[length(summ) + 1]] <- data.table(
      sample = smp, variable = nm, statistic = "morans_I",
      value = r$I, perm_p = r$p, z = r$z, n_tiles = nrow(d))
  }

  # Depth adjustment. If the chrX structure is induced by the (strongly
  # structured) depth field, residualising on depth removes it. If I is
  # unchanged, depth is not the explanation.
  fit1 <- stats::lm(x_ratio ~ log(x_n), data = d)
  r1 <- morans_I(e, stats::residuals(fit1))
  msg("  Moran's I  %-16s %+.4f   perm p = %.3f   (depth slope %+.4f / log unit)",
      "ratio | depth", r1$I, r1$p, stats::coef(fit1)[2])
  # Plus the autosomal ratio, which absorbs per-tile mapping bias.
  fit2 <- stats::lm(x_ratio ~ log(x_n) + log(a_n) + a_ratio, data = d)
  r2 <- morans_I(e, stats::residuals(fit2))
  msg("  Moran's I  %-16s %+.4f   perm p = %.3f", "ratio | depth+auto", r2$I, r2$p)
  for (nm in c("ratio_given_depth", "ratio_given_depth_auto")) {
    r <- if (nm == "ratio_given_depth") r1 else r2
    summ[[length(summ) + 1]] <- data.table(
      sample = smp, variable = nm, statistic = "morans_I",
      value = r$I, perm_p = r$p, z = r$z, n_tiles = nrow(d))
  }

  # The beta-binomial null, NEXT_ANALYSIS task 7. Fitted on the AUTOSOMAL
  # control (no monoallelic biology, same tiles, same library) and centred on
  # the sample's OWN pooled chrX ratio, not on 0.5.
  #
  # n_eff, not n. On a duplicate-inclusive read tree the counts are inflated by
  # the duplication factor and phi absorbs all of it, so phi here is mostly a
  # statement about PCR, not about tissue. It is reported for completeness and
  # must not be quoted as biological overdispersion.
  d[, `:=`(x_neff = x_n / DUP_FACTOR, a_neff = a_n / DUP_FACTOR)]
  v_obs  <- stats::var(d$a_ratio)
  v_bin  <- mean(d$a_ratio * (1 - d$a_ratio) / d$a_neff)
  phi    <- v_obs / v_bin
  tau    <- sqrt(max(v_obs - v_bin, 0))
  xbar   <- sum(d$x_a1) / sum(d$x_n)
  d[, null_sd := sqrt(phi * xbar * (1 - xbar) / x_neff + tau^2)]
  d[, z_bb := (x_ratio - xbar) / null_sd]
  d[, q_bb := p.adjust(2 * stats::pnorm(-abs(z_bb)), "BH")]
  msg("  autosomal control: phi = %.1fx binomial, tau = %.4f (on n_eff = n/%.1f)",
      phi, tau, DUP_FACTOR)
  msg("  pooled chrX ratio %.4f; tiles differing from it at FDR 5%%: %d (%.1f%%)",
      xbar, sum(d$q_bb < 0.05), 100 * mean(d$q_bb < 0.05))
  summ[[length(summ) + 1]] <- data.table(
    sample = smp, variable = "beta_binomial_null", statistic = "phi",
    value = phi, perm_p = NA_real_, z = NA_real_, n_tiles = nrow(d))

  d <- cbind(d, local_moran(e, d$x_ratio))
  # Two different questions, and the second is usually the informative one here.
  # Surviving BH over ~4000 tiles needs a CONCENTRATED signal (hotspots). A
  # global I that is significant while nothing survives BH means the structure
  # is real but DIFFUSE - spread thinly over many tiles rather than pooled into
  # a few regions. So report the enrichment of small p as well as the FDR count,
  # because "0 tiles at FDR 5%" on its own reads as "no structure" and that is
  # not what it means.
  enr <- mean(d$lisa_p < 0.05) / 0.05
  msg("  LISA: %d tiles at FDR 5%%; %.1f%% of tiles at raw p < 0.05 (%.1fx the 5%% null)",
      sum(d$lisa_q < 0.05), 100 * mean(d$lisa_p < 0.05), enr)
  msg("        -> %s", if (sum(d$lisa_q < 0.05) > 0) "concentrated enough to localise"
                       else if (enr > 1.3) "real but DIFFUSE; no tile is individually callable"
                       else "consistent with no local structure")
  summ[[length(summ) + 1]] <- data.table(
    sample = smp, variable = "lisa_p05_enrichment", statistic = "ratio_to_null",
    value = enr, perm_p = NA_real_, z = NA_real_, n_tiles = nrow(d))

  # ---- correlogram: I at ring distance 1..MAX_LAG -------------------------
  #
  # WHAT THE DECAY LENGTH IS AND IS NOT. It is NOT a patch size. The CAST X is
  # inactive in every cell of this cross, so there are no clonal XCI domains for
  # a correlogram to measure, and NEXT_ANALYSIS.md already rules out reading a
  # scale off rho(s) for the same reason. What it does separate is short-range
  # from long-range structure, which have different causes:
  #   lag 1-2 (<= ~130 um)  tile-boundary bleed, segmentation, local depth
  #                         correlation -> technical.
  #   lag 5-10 (0.3-0.6 mm) regional/anatomical -> most plausibly cell-type
  #                         composition, which is not controlled here.
  #
  # Three fields, because the shape is only interpretable against the other two:
  #   x_ratio        the field of interest
  #   resid          the same field with depth and the autosomal ratio removed;
  #                  if the decay survives this it is not a coverage gradient
  #   a_ratio        the negative control; it should sit on the null at EVERY
  #                  lag, and if it does not, no lag of the chrX curve is safe
  #
  # A global gradient across the section produces slow decay with no
  # characteristic scale, which is not patchiness either - the flat C(d) excess
  # already on record is exactly that. Read the shape: a knee is a scale, a
  # straight slow decline is a gradient.
  resid2 <- stats::residuals(fit2)
  cg <- rbindlist(lapply(seq_len(MAX_LAG), function(lag) {
    re <- ring_edges(d$trow, d$tcol, lag)
    if (re$W < 50) return(NULL)
    rbindlist(lapply(list(x_ratio = d$x_ratio, resid = resid2, a_ratio = d$a_ratio),
      function(v) { r <- morans_I(re, v, NPERM_LAG)
                    data.table(I = r$I, perm_p = r$p, z = r$z) }), idcol = "field")[
      , `:=`(sample = smp, lag = lag, dist_um = lag * TILE_UM, n_pairs = re$W)]
  }))
  if (nrow(cg)) {
    msg("  correlogram (chrX ratio), %s:", "I at ring distance")
    for (k in seq_len(nrow(cg[field == "x_ratio"]))) {
      r <- cg[field == "x_ratio"][k]
      msg("    lag %2d (%4.0f um)  I = %+.4f  p = %.3f  (%d pairs)",
          r$lag, r$dist_um, r$I, r$perm_p, r$n_pairs)
    }
    # 12 lags x 3 fields per sample. An isolated lag at p ~ 0.01 is what that many
    # tests produce on a null field, so the shape of the curve is the evidence and
    # a single flagged lag is not. BH within each field makes that explicit.
    cg[, perm_q := p.adjust(perm_p, "BH"), by = field]
    xr_cg <- cg[field == "x_ratio"]
    first_ns <- xr_cg[perm_p >= 0.05, min(lag)]

    # The decay length, reported because it is the number people ask for, and
    # caveated because the curve has no knee. An exponential fit to a monotone
    # decline always returns a length; that length is a description of a
    # GRADIENT, not evidence of a characteristic scale. A real scale would show
    # as a knee - a plateau then a drop - and would survive changing MAX_LAG.
    pos <- xr_cg[I > 0]
    if (nrow(pos) >= 4) {
      ef <- stats::lm(log(I) ~ lag, data = pos)
      b  <- unname(stats::coef(ef)[2]); r2 <- summary(ef)$r.squared
      if (b < 0 && r2 >= 0.5) {
        msg("  e-folding length %.0f um (R2 = %.2f on log I ~ lag). A LENGTH, not a scale:",
            -TILE_UM / b, r2)
        msg("    with no knee in the curve this describes a gradient's steepness only.")
      } else if (b < 0) {
        msg("  no usable decay length: log I ~ lag fits at R2 = %.2f, so the curve", r2)
        msg("    is not an exponential decline and the number would be meaningless.")
      }
    }
    if (is.finite(first_ns)) {
      msg("  -> chrX I first reaches the null at lag %d (%.0f um).",
          first_ns, first_ns * TILE_UM)
      msg("     That is a correlation length, NOT a patch size - there is no")
      msg("     mosaic in this cross. See the note above this block.")
    } else {
      msg("  -> chrX I is still above the null at lag %d (%.0f um): no decay",
          MAX_LAG, MAX_LAG * TILE_UM)
      msg("     within the range tested. A gradient, not a scale. Raise MAX_LAG.")
    }
    ac <- cg[field == "a_ratio" & perm_q < 0.05]
    if (nrow(ac)) {
      msg("  WARNING: the autosomal control is non-null after BH at lag(s) %s.",
          paste(ac$lag, collapse = ", "))
      msg("    The control carries no monoallelic biology, so structure in it is")
      msg("    technical and the chrX curve is not safe to read at those lags.")
    } else if (nrow(cg[field == "a_ratio" & perm_p < 0.05])) {
      msg("  (autosomal control has isolated raw p < 0.05 lag(s) %s, none surviving",
          paste(cg[field == "a_ratio" & perm_p < 0.05]$lag, collapse = ", "))
      msg("   BH across lags - consistent with multiple testing, control is clean.)")
    }
    correlo[[smp]] <- cg
  }
  per_tile[[smp]] <- d
}

if (!length(per_tile)) stop("no samples produced tiles - check TREE and MIN_X")
summ <- rbindlist(summ); tiles <- rbindlist(per_tile, fill = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
fwrite(summ,  file.path(OUT_DIR, sprintf("tile_spatial_structure_%s.csv", TAG)))
fwrite(tiles, file.path(OUT_DIR, sprintf("tile_spatial_structure_tiles_%s.csv", TAG)))
if (length(correlo))
  fwrite(rbindlist(correlo),
         file.path(OUT_DIR, sprintf("tile_spatial_correlogram_%s.csv", TAG)))

##### ---------------------------- PLOTS ----------------------------- #####
pdf(file.path(OUT_DIR, sprintf("tile_spatial_structure_%s.pdf", TAG)),
    width = 11, height = 5.5)
gi <- summ[statistic == "morans_I"]
gi[, variable := factor(variable, unique(variable))]
# Facet in the order SAMPLES was given (9w then 78w), not alphabetically -
# "78w" sorts before "9w" as a string and puts aged on the left.
gi[, sample := factor(sample, SAMPLES[SAMPLES %in% unique(sample)])]
print(
  ggplot(gi, aes(variable, value, fill = perm_p < 0.05)) +
    geom_col() + facet_wrap(~sample) + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#b02a2a", `FALSE` = "grey70"),
                      name = "perm p < 0.05") +
    labs(title = "Global Moran's I per tile field",
         subtitle = paste0("Tree ", TREE, ", ", COUNT_UNIT, " level, tiles with >= ",
                           MIN_X, " chrX ", COUNT_UNIT, "s. ", NPERM, " permutations."),
         caption = paste("auto_ratio is the negative control: it should be ~0.",
                         "\nIf chrX_ratio survives 'ratio | depth' it is not a coverage artefact.",
                         "\nThis is NOT a test for XCI patches - the CAST X is inactive in every cell."),
         x = NULL, y = "Moran's I") + theme_bw())

if (length(correlo)) {
  cgall <- rbindlist(correlo)
  cgall[, sample := factor(sample, SAMPLES[SAMPLES %in% unique(sample)])]
  cgall[, field := factor(field, c("x_ratio", "resid", "a_ratio"),
                          c("chrX ratio", "chrX ratio | depth + autosomal",
                            "autosomal ratio (control)"))]
  print(
    ggplot(cgall, aes(dist_um, I, colour = field)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      geom_line() +
      geom_point(aes(shape = perm_p < 0.05), size = 2) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                         name = "perm p < 0.05") +
      scale_colour_manual(values = c("#b02a2a", "#1b6ca8", "grey55"), name = NULL) +
      facet_wrap(~sample) +
      labs(title = "Correlogram: Moran's I by ring distance",
           subtitle = sprintf("Tile side %.0f um; ring lag 1..%d; %d permutations per lag",
                              TILE_UM, MAX_LAG, NPERM_LAG),
           caption = paste(
             "Where the curve meets the dashed null is a CORRELATION LENGTH, not a patch size:",
             "\nthe CAST X is inactive in every cell, so there are no clonal XCI domains to size.",
             "\nDecay by lag 1-2 is technical (tile bleed, local depth); persistence to lag 5-10 is regional,",
             "\nmost plausibly cell-type composition. A straight slow decline with no knee is a gradient, not a scale.",
             "\nThe control must sit on the null at EVERY lag, or no lag of the chrX curve is safe."),
           x = "ring distance (um)", y = "Moran's I") +
      theme_bw() + theme(plot.caption = element_text(size = 7, hjust = 0)))
}

for (smp in names(per_tile)) {
  d <- per_tile[[smp]]
  print(
    ggplot(d, aes(tcol, -trow, fill = lisa_q < 0.05)) +
      geom_tile() + coord_equal() +
      scale_fill_manual(values = c(`TRUE` = "#b02a2a", `FALSE` = "grey88"),
                        name = "LISA FDR < 5%") +
      labs(title = sprintf("%s - where the local structure is", smp),
           subtitle = sprintf("%d of %d tiles locally coherent",
                              sum(d$lisa_q < 0.05), nrow(d)),
           x = NULL, y = NULL) + theme_bw())
  print(
    ggplot(d, aes(x_n, x_ratio)) +
      geom_point(size = 0.4, alpha = 0.3) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = "#b02a2a") +
      scale_x_log10() +
      labs(title = sprintf("%s - the depth confound, directly", smp),
           subtitle = "Any slope here can induce ratio structure from depth structure alone",
           x = sprintf("chrX %ss per tile", COUNT_UNIT), y = "chrX B6 fraction") +
      theme_bw())
}
invisible(dev.off())

msg("\nWrote tile_spatial_structure_%s.{csv,pdf} to %s", TAG, OUT_DIR)
msg("")
msg("NOT YET CONTROLLED: cell-type composition. Cell types differ in escape and")
msg("are spatially organised, so composition can produce every number above with")
msg("escape itself spatially uniform. To close it, join a per-tile composition")
msg("vector from spatial/spatial_cell_annotation.R and add it to the fit2 model,")
msg("then re-read 'ratio | depth+auto'. Until then a positive I says structure")
msg("EXISTS, not that escape itself varies.")
