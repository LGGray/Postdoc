# Map the per-tile Allelome.PRO2 chrX allelic ratio back onto the slide.
#
# A quick look at the tile run while it is still in flight: draws every tile
# that has a locus table, and a grey box for every tile that was handed to
# sinto but has not been scored yet. The faint footprint underneath is every
# in-tissue tile, including the ones below SINTO_MIN_UMI that were never
# submitted - so "grey" always means "pending", never "no tissue".
#
# Reads whatever exists, in this order:
#   1. ase/<sample>/Allelome.PRO2_tiles_<um>um/<tile>/locus_table.txt  (collected)
#   2. $SCRATCH/spatial_tiles_<sample>_<um>um/allelome/<dated>/locus_table.txt
# so it works before spatial_sinto_tiles.slurm reaches its collect stage.
#
# Run on the cluster:  sbatch ~/Postdoc/slurm/spatial_tile_map.slurm
#              or:     TILE_UM=64 Rscript ~/Postdoc/spatial/tile_ratio_map.R

.libPaths(c("~/R/matrix-dev", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

##### ---------------------- CONFIG ---------------------- #####
BASE    <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial"
SAMPLES <- strsplit(Sys.getenv("SAMPLES", "9w,78w"), ",")[[1]]
TILE_UM <- as.integer(Sys.getenv("TILE_UM", "64"))
BIN     <- "square_002um"
ANNOT_BASE <- "chr_annotation_mm39.bed"
OUT_PDF <- Sys.getenv("OUT_PDF",
                      file.path(BASE, "ase",
                                sprintf("tile_ratio_map_%dum.pdf", TILE_UM)))

# Draw the H&E underneath the ratio panel. Needs the png package; skipped with a
# message if it is not installed rather than failing the whole script.
USE_HE  <- TRUE

# The capture array can sit rotated relative to the image. Axis-aligned squares
# are drawn either way, so a large rotation would show as a staircased edge -
# flip these if the map comes out mirrored against the H&E.
FLIP_X  <- FALSE
FLIP_Y  <- FALSE

# Chromosomes excluded from the autosomal control.
AUTOSOMES <- paste0("chr", 1:19)

# Per-tile autosomal ratios are overdispersed relative to binomial, so the null
# band is the OBSERVED autosomal spread, not the binomial one. Estimated per
# sample: 9w came out at 0.027 and 78w at 0.051 on comparable depth, so a single
# hardcoded value fits one sample and inflates every z in the other. AUTO_SD
# below is only the fallback for a sample with too few scored tiles to estimate.
AUTO_SD_FALLBACK <- 0.032
MIN_TILES_FOR_SD <- 30
Z_CALL  <- 3

# Diverging pair, neutral MID-POINT DELIBERATELY NEAR-WHITE rather than the
# usual grey: grey is spent on NA here, and a grey midpoint would make a
# genuinely biallelic tile indistinguishable from an unprocessed one.
COL_CAST <- "#184f95"   # ratio -> 0, CAST (A2) dominant
COL_MID  <- "#f0efec"   # 0.5
COL_BL6  <- "#b02a2a"   # ratio -> 1, Bl6 (A1) dominant
# ---- the lab allelic-ratio ramp, from OCM_heart/allelic_ratio/ -------------
# Verbatim from 02_whole_chrX.R and 07_core_escape_cutoff_sweep.R, which agree
# on it, so these tile maps bin and colour exactly like the OCM UMAPs.
#
# Blue = A2/CAST-dominant, green ~ biallelic, olive/yellow -> orange -> dark red
# = increasingly monoallelic A1/Bl6. The "red to green" half is the top of it.
#
# Preferred over the ASE/adult_aged_snRNAseq.R ramp (green->yellow->orange->red
# over 0.5-1.0) for three reasons: it covers the whole 0-1 range, so the
# CAST-biased tiles get a real colour instead of being clamped at the bottom of
# a scale that has none; it is binned rather than continuous, which is honest
# about per-tile precision (at ~74 reads the SE is 0.05, so 0.1-wide bins are
# roughly the resolution the data supports and a smooth ramp would imply more);
# and grey-for-NA is already its convention.
#
# One caveat worth knowing: lightness is not monotone along this ramp - it
# peaks at 0.80-0.90 (L* 0.74) and falls to 0.41 at 0.95-1.00, so the darkest
# bins sit at BOTH ends. In a UMAP with few extreme cells that is harmless; on
# a 464-tile map, do not read darkness as "high ratio". Read the legend.
# CVD: poles separate at dE 19.7-30.8, but adjacent bins get as close as 3.2
# (0.30-0.40 vs 0.40-0.50), so neighbouring bins are not reliably separable -
# which is why the calls panel exists alongside this one.
OCM_BREAKS <- c(seq(0, 0.9, by = 0.1), 0.95, 1.0)
OCM_COLORS <- c("#2B3186", "#3B5FB6", "#38749F", "#367373", "#2D6E5D", "#1E652D",
                "#658C2D", "#8D9F25", "#B3B112", "#C97314", "#8B1913")
OCM_LABELS <- c("0.00-0.10", "0.10-0.20", "0.20-0.30", "0.30-0.40", "0.40-0.50",
                "0.50-0.60", "0.60-0.70", "0.70-0.80", "0.80-0.90",
                "0.90-0.95", "0.95-1.00")
names(OCM_COLORS) <- OCM_LABELS

# grey70 is the OCM na.value. Kept for the pending tiles specifically, so the
# grey means the same thing here as it does there; the never-submitted
# footprint stays lighter so the two are still distinguishable.
COL_NA <- "grey70"
COL_FOOT <- "#e6e5e1"   # in-tissue, never submitted
##### ------------------------------------------------------ #####

msg <- function(...) message(sprintf(...))

read_positions <- function(dir) {
  pq <- file.path(dir, "spatial", "tissue_positions.parquet")
  LEGACY <- c("barcode", "in_tissue", "array_row", "array_col",
              "pxl_row_in_fullres", "pxl_col_in_fullres")
  for (cs in file.path(dir, "spatial",
                       c("tissue_positions.tsv", "tissue_positions.csv",
                         "tissue_positions_list.csv"))) {
    if (!file.exists(cs)) next
    p <- as.data.table(fread(cs))
    setnames(p, tolower(names(p)))
    # tissue_positions_list.csv carries no header, so fread invents V1..V6 and
    # every lookup below would miss without ever erroring. Name them.
    if (all(grepl("^v[0-9]+$", names(p))) && ncol(p) >= length(LEGACY)) {
      setnames(p, seq_along(LEGACY), LEGACY)
    }
    return(p)
  }
  if (file.exists(pq)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("tissue_positions.parquet needs the 'arrow' package, and no ",
           ".tsv/.csv fallback is present under ", file.path(dir, "spatial"))
    }
    p <- as.data.table(arrow::read_parquet(pq))
    setnames(p, tolower(names(p)))
    return(p)
  }
  stop("No tissue_positions file under ", file.path(dir, "spatial"))
}

# scalefactors_json.json, without taking a jsonlite dependency for one number.
lowres_scalef <- function(dir) {
  f <- file.path(dir, "spatial", "scalefactors_json.json")
  if (!file.exists(f)) return(NA_real_)
  txt <- paste(readLines(f, warn = FALSE), collapse = " ")
  m <- regmatches(txt, regexpr('"tissue_lowres_scalef"\\s*:\\s*[0-9.eE+-]+', txt))
  if (!length(m)) return(NA_real_)
  as.numeric(sub('.*:\\s*', '', m))
}

# One locus table -> one row: chrX counts, and every autosome pooled.
read_locus <- function(path) {
  lt <- tryCatch(fread(path, showProgress = FALSE), error = function(e) NULL)
  if (is.null(lt) || !nrow(lt)) return(NULL)
  setnames(lt, tolower(names(lt)))
  need <- c("chr", "a1_reads", "a2_reads", "total_reads")
  if (!all(need %in% names(lt))) return(NULL)
  x <- lt[chr == "chrX"]
  a <- lt[chr %in% AUTOSOMES]
  if (!nrow(x)) return(NULL)
  data.table(
    x_a1 = sum(x$a1_reads), x_a2 = sum(x$a2_reads), x_n = sum(x$total_reads),
    a_a1 = sum(a$a1_reads), a_a2 = sum(a$a2_reads), a_n = sum(a$total_reads)
  )
}

collect_sample <- function(smp) {
  ase_dir   <- file.path(BASE, "ase", smp)
  final_dir <- file.path(ase_dir, sprintf("Allelome.PRO2_tiles_%dum", TILE_UM))
  map_tsv   <- file.path(ase_dir, sprintf("sinto_tiles_%dum.tsv", TILE_UM))
  if (!file.exists(map_tsv)) {
    msg("  no sinto map at %s - skipping %s", map_tsv, smp); return(NULL)
  }

  # The submitted universe. This, not the set of locus tables, is the
  # denominator - it is what defines a grey box.
  map <- fread(map_tsv, header = FALSE, col.names = c("barcode", "tile"))
  submitted_tiles <- unique(map$tile)

  # Collected layout first, then the live scratch layout.
  lt_paths <- character(0)
  if (dir.exists(final_dir)) {
    lt_paths <- Sys.glob(file.path(final_dir, "*", "locus_table.txt"))
    names(lt_paths) <- basename(dirname(lt_paths))
  }
  scratch <- Sys.getenv("SCRATCH")
  if (!length(lt_paths) && nzchar(scratch)) {
    ap2 <- file.path(scratch, sprintf("spatial_tiles_%s_%dum", smp, TILE_UM),
                     "allelome")
    p <- Sys.glob(file.path(ap2, paste0("*_", ANNOT_BASE, "_*"), "locus_table.txt"))
    if (length(p)) {
      # 2026_08_27_tile_64um_r0019_c0050_chr_annotation_mm39.bed_1 -> tile name
      nm <- basename(dirname(p))
      nm <- sub("^[0-9]{4}_[0-9]{2}_[0-9]{2}_", "", nm)
      nm <- sub(paste0("_", ANNOT_BASE, "_[0-9]+$"), "", nm)
      lt_paths <- p; names(lt_paths) <- nm
      msg("  reading from scratch: %s", ap2)
    }
  }
  if (!length(lt_paths)) {
    msg("  no locus tables found for %s - skipping", smp); return(NULL)
  }
  # A resubmitted job leaves _1.._n for the same tile; keep the last.
  lt_paths <- lt_paths[!duplicated(names(lt_paths), fromLast = TRUE)]

  scored <- rbindlist(lapply(seq_along(lt_paths), function(i) {
    r <- read_locus(lt_paths[[i]])
    if (is.null(r)) return(NULL)
    r[, tile := names(lt_paths)[i]][]
  }))
  msg("  %d locus tables, %d with a chrX row, of %d submitted tiles",
      length(lt_paths), nrow(scored), length(submitted_tiles))
  if (!nrow(scored)) return(NULL)

  # --- tile geometry, from the bins themselves ---
  bin_dir <- file.path(BASE, smp, "outs", "binned_outputs", BIN)
  pos <- read_positions(bin_dir)
  if (!all(c("array_row", "array_col") %in% names(pos))) {
    stop("tissue_positions for ", smp, " has no array_row/array_col")
  }
  it <- grep("^in_tissue", names(pos), value = TRUE)[1]
  if (!is.na(it)) pos <- pos[get(it) == 1]
  # arrow can hand back integer64, which %/% and matrix indexing mishandle.
  pos[, `:=`(array_row = as.numeric(array_row), array_col = as.numeric(array_col),
             px = as.numeric(pxl_col_in_fullres), py = as.numeric(pxl_row_in_fullres))]

  k <- TILE_UM / 2                      # bins per tile side
  pos[, tile := sprintf("tile_%dum_r%04d_c%04d", TILE_UM,
                        as.integer(array_row %/% k), as.integer(array_col %/% k))]

  he_png <- if (USE_HE) file.path(bin_dir, "spatial", "tissue_lowres_image.png")
            else NA_character_

  sf <- lowres_scalef(bin_dir)
  if (!is.finite(sf)) { sf <- 1; msg("  no lowres scalefactor - plotting in full-res pixels") }
  geom <- pos[, .(x = mean(px) * sf, y = mean(py) * sf, n_bins = .N), by = tile]

  # Tile side in the plotted units. Taken from the array->pixel slopes so it is
  # right even where a tile is only partly covered by tissue.
  sx <- coef(lm(x ~ array_col + array_row,
                data = pos[, .(x = px * sf, array_col, array_row)]))
  sy <- coef(lm(y ~ array_col + array_row,
                data = pos[, .(y = py * sf, array_col, array_row)]))
  side <- k * mean(c(sqrt(sx[["array_col"]]^2 + sy[["array_col"]]^2),
                     sqrt(sx[["array_row"]]^2 + sy[["array_row"]]^2)))
  msg("  tile side %.1f plot units (%d bins), lowres scalef %.4f", side, k, sf)

  d <- merge(geom, scored, by = "tile", all.x = TRUE)
  d[, submitted := tile %in% submitted_tiles]
  # Row/col off the tile name. Needed for adjacency: whether a call clusters in
  # space is the difference between biology and per-tile noise, and pixel
  # centroids cannot answer it as cleanly as the integer grid can.
  d[, `:=`(trow = as.integer(sub(".*_r([0-9]+)_c[0-9]+$", "\\1", tile)),
           tcol = as.integer(sub(".*_r[0-9]+_c([0-9]+)$", "\\1", tile)))]
  d[, `:=`(x_ratio = fifelse(x_n > 0, x_a1 / x_n, NA_real_),
           a_ratio = fifelse(a_n > 0, a_a1 / a_n, NA_real_))]
  sc <- d[!is.na(a_ratio) & a_n > 0]
  if (nrow(sc) >= MIN_TILES_FOR_SD) {
    # Observed variance in the autosomal ratio, less the binomial variance it
    # already contains. What remains is the technical per-tile spread, which is
    # the null band chrX has to clear.
    v_obs  <- var(sc$a_ratio)
    v_binom <- mean(sc$a_ratio * (1 - sc$a_ratio) / sc$a_n)
    auto_sd <- sqrt(max(v_obs - v_binom, 0))
    msg("  autosomal spread: sd %.4f observed, %.4f binomial -> null band %.4f (%.1fx binomial)",
        sqrt(v_obs), sqrt(v_binom), auto_sd, sqrt(v_obs) / sqrt(v_binom))
  } else {
    auto_sd <- AUTO_SD_FALLBACK
    msg("  only %d scored tiles - using the %.3f fallback null band",
        nrow(sc), auto_sd)
  }
  d[, se := sqrt(a_ratio * (1 - a_ratio) / x_n + auto_sd^2)]
  d[, z := (x_ratio - a_ratio) / se]
  # Binned exactly as OCM_heart/allelic_ratio does it: same breaks, same
  # include.lowest/right, same labels, so a tile falls in the same bin a cell
  # with that ratio would.
  d[, x_bin := cut(x_ratio, breaks = OCM_BREAKS, include.lowest = TRUE,
                   right = TRUE, labels = OCM_LABELS)]
  d[, call := fcase(
      is.na(x_ratio) &  submitted, "pending",
      is.na(x_ratio) & !submitted, "not submitted",
      z >  Z_CALL, "Bl6-skewed",
      z < -Z_CALL, "CAST-skewed",
      default = "mixed")]
  d[, `:=`(sample = smp, side = side, FLIPX = FLIP_X, FLIPY = FLIP_Y,
         he = he_png, auto_sd = auto_sd)]
  if (FLIP_X) d[, x := -x]
  if (FLIP_Y) d[, y := -y]

  d
}

##### ------------------------- draw ------------------------ #####

# Recessive frame: the marks carry everything, the axes carry nothing. Pixel
# coordinates are not a quantity anyone reads off a slide.
theme_slide <- theme_void(base_size = 10) +
  theme(legend.position = "right",
        legend.key.height = unit(14, "pt"),
        legend.key.width  = unit(8, "pt"),
        plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "#52514e"),
        plot.caption  = element_text(size = 7, colour = "#52514e", hjust = 0))

he_layer <- function(d) {
  f <- d$he[1]
  if (is.na(f) || !file.exists(f)) {
    msg("  no %s - drawing without the H&E", basename(as.character(f)))
    return(NULL)
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    msg("  png package absent - drawing without the H&E"); return(NULL)
  }
  img <- png::readPNG(f)
  # The lowres PNG is the full-res image scaled by tissue_lowres_scalef, which
  # is exactly the space the tile centroids were put into above, so the raster
  # spans 0..ncol / 0..nrow of the PNG itself.
  w <- dim(img)[2]; h <- dim(img)[1]
  xr <- if (isTRUE(d$FLIPX[1])) c(-w, 0) else c(0, w)
  yr <- if (isTRUE(d$FLIPY[1])) c(-h, 0) else c(0, h)
  annotation_raster(img, xmin = xr[1], xmax = xr[2], ymin = yr[1], ymax = yr[2],
                    interpolate = TRUE)
}

# y is reversed everywhere: image coordinates run downwards from the top-left.
base_map <- function(d, he = FALSE) {
  g <- ggplot(d, aes(x, y))
  if (he) { l <- he_layer(d); if (!is.null(l)) g <- g + l }
  g +
    geom_tile(data = d[submitted == FALSE], width = d$side[1], height = d$side[1],
              fill = COL_FOOT, colour = NA) +
    coord_fixed() + scale_y_reverse() + theme_slide
}

panel_ratio <- function(d, he = FALSE) {
  n_ok <- sum(!is.na(d$x_ratio)); n_sub <- sum(d$submitted)
  base_map(d, he) +
    geom_tile(data = d[submitted & is.na(x_ratio)],
              width = d$side[1], height = d$side[1], fill = COL_NA, colour = NA) +
    geom_tile(data = d[!is.na(x_ratio)], aes(fill = x_ratio),
              width = d$side[1], height = d$side[1], colour = NA) +
    scale_fill_gradient2(low = COL_CAST, mid = COL_MID, high = COL_BL6,
                         midpoint = 0.5, limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),
                         name = "Bl6 (A1)\nfraction") +
    labs(title = sprintf("%s - chrX allelic ratio per %d um tile", d$sample[1], TILE_UM),
         subtitle = sprintf("%d of %d submitted tiles scored; grey = submitted, not yet scored",
                            n_ok, n_sub),
         caption = paste("Faint boxes are in-tissue tiles below SINTO_MIN_UMI that were never submitted.",
                         "\nMidpoint is near-white, not grey, so a biallelic tile stays distinct from an unscored one."))
}

panel_ratio_ocm <- function(d, he = FALSE) {
  n_ok <- sum(!is.na(d$x_ratio)); n_sub <- sum(d$submitted)
  lo <- sum(d$x_ratio < 0.5, na.rm = TRUE)
  base_map(d, he) +
    geom_tile(data = d[submitted & is.na(x_bin)],
              width = d$side[1], height = d$side[1], fill = COL_NA, colour = NA) +
    geom_tile(data = d[!is.na(x_bin)], aes(fill = x_bin),
              width = d$side[1], height = d$side[1], colour = NA) +
    scale_fill_manual(values = OCM_COLORS, limits = OCM_LABELS, drop = FALSE,
                      na.value = COL_NA, name = "Allelic ratio") +
    guides(fill = guide_legend(ncol = 1, reverse = TRUE)) +
    labs(title = sprintf("%s - chrX allelic ratio per %d um tile (OCM bins)",
                         d$sample[1], TILE_UM),
         subtitle = sprintf("%d of %d submitted tiles scored; %d (%.0f%%) below 0.5. Grey = submitted, not yet scored.",
                            n_ok, n_sub, lo, 100 * lo / max(n_ok, 1)),
         caption = paste("Bins and colours from OCM_heart/allelic_ratio (02_whole_chrX.R).",
                         "\nLightness is not monotone along this ramp - the darkest bins are at both ends, so read the legend, not the darkness."))
}

panel_auto <- function(d) {
  base_map(d) +
    geom_tile(data = d[submitted & is.na(a_ratio)],
              width = d$side[1], height = d$side[1], fill = COL_NA, colour = NA) +
    geom_tile(data = d[!is.na(a_ratio)], aes(fill = a_ratio),
              width = d$side[1], height = d$side[1], colour = NA) +
    scale_fill_gradient2(low = COL_CAST, mid = COL_MID, high = COL_BL6,
                         midpoint = 0.5, limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),
                         name = "Bl6 (A1)\nfraction") +
    labs(title = sprintf("%s - autosomal control, same tiles, same scale", d$sample[1]),
         subtitle = sprintf("Should be flat near-white throughout. Observed sd %.3f on a median of %d reads.",
                            sd(d$a_ratio, na.rm = TRUE),
                            as.integer(median(d$a_n, na.rm = TRUE))),
         caption = paste("Any structure here is technical and invalidates the chrX panel over the same tiles.",
                         "\nDeliberately NOT the OCM ramp: autosomal ratios all sit in 0.40-0.60, which is two of its bins, so it would show nothing."))
}

panel_auto_zoom <- function(d) {
  # Same data as panel_auto, on a +-3 null-band scale instead of 0-1. The wide
  # panel answers "is the control flat"; this one answers "and if it is not,
  # is the residual structured or speckle" - which is the question that decides
  # whether the chrX calls over these tiles can be trusted.
  sdv <- d$auto_sd[1]
  lim <- max(3 * sdv, 0.05)
  d2 <- copy(d)[, a_clip := pmin(pmax(a_ratio - 0.5, -lim), lim)]
  base_map(d2) +
    geom_tile(data = d2[submitted & is.na(a_ratio)],
              width = d$side[1], height = d$side[1], fill = COL_NA, colour = NA) +
    geom_tile(data = d2[!is.na(a_clip)], aes(fill = a_clip),
              width = d$side[1], height = d$side[1], colour = NA) +
    scale_fill_gradient2(low = COL_CAST, mid = COL_MID, high = COL_BL6,
                         midpoint = 0, limits = c(-lim, lim),
                         name = "autosomal\ndeviation\nfrom 0.5") +
    labs(title = sprintf("%s - autosomal control, opened up to +-%.3f", d$sample[1], lim),
         subtitle = sprintf("Null band %.4f. Speckle is technical noise; a patch or a gradient is a bias with geometry.",
                            sdv),
         caption = "Clipped at the scale limits, so an extreme tile reads as the pole rather than dominating the ramp.")
}

panel_depth <- function(d) {
  base_map(d) +
    geom_tile(data = d[submitted & is.na(x_n)],
              width = d$side[1], height = d$side[1], fill = COL_NA, colour = NA) +
    geom_tile(data = d[!is.na(x_n)], aes(fill = x_n),
              width = d$side[1], height = d$side[1], colour = NA) +
    scale_fill_gradient(low = "#cde2fb", high = "#0d366b", trans = "sqrt",
                        name = "chrX reads") +
    labs(title = sprintf("%s - informative chrX reads per tile", d$sample[1]),
         subtitle = sprintf("median %d; a ratio needs ~20 reads to separate 0.8 from 0.5",
                            as.integer(median(d$x_n, na.rm = TRUE))))
}

panel_call <- function(d) {
  lv <- c("Bl6-skewed", "CAST-skewed", "mixed", "pending", "not submitted")
  d2 <- copy(d)[, call := factor(call, levels = lv)]
  base_map(d2) +
    geom_tile(data = d2[call != "not submitted"], aes(fill = call),
              width = d$side[1], height = d$side[1], colour = NA) +
    scale_fill_manual(values = c("Bl6-skewed" = COL_BL6, "CAST-skewed" = COL_CAST,
                                 "mixed" = COL_MID, "pending" = COL_NA,
                                 "not submitted" = COL_FOOT),
                      drop = FALSE, name = NULL) +
    # COL_FOOT on a white legend background is an invisible swatch; the border
    # is what makes "not submitted" readable as a key rather than a gap.
    guides(fill = guide_legend(override.aes = list(colour = "#c3c2b7"))) +
    labs(title = sprintf("%s - call vs each tile's own autosomal ratio", d$sample[1]),
         subtitle = sprintf("|z| > %d, against this sample's own %.3f autosomal null band: %s",
                            Z_CALL, d$auto_sd[1],
                            paste(sprintf("%s %d", lv[1:4],
                                          sapply(lv[1:4], function(l) sum(d2$call == l, na.rm = TRUE))),
                                  collapse = ", ")))
}

clustering_test <- function(d, label, n_perm = 2000L) {
  sc <- d[!is.na(x_ratio)]
  if (!nrow(sc)) return(invisible(NULL))
  k <- sum(sc$call == label)
  if (k < 5) { msg("  %s: only %d tiles - too few to test for clustering", label, k); return(invisible(NULL)) }

  # 4-neighbour adjacency among SCORED tiles only. Using every tile would count
  # a pending neighbour as "not a match" and bias the statistic toward
  # dispersion purely because the run is incomplete.
  key <- function(r, cc) r * 100000L + cc
  present <- key(sc$trow, sc$tcol)
  # Observed: adjacent pairs that are both `label`.
  count_pairs <- function(is_lab) {
    kk <- present[is_lab]
    sum(c((key(sc$trow[is_lab] + 1L, sc$tcol[is_lab]) %in% kk),
          (key(sc$trow[is_lab], sc$tcol[is_lab] + 1L) %in% kk)))
  }
  obs <- count_pairs(sc$call == label)
  # Null: the same number of labels placed at random over the scored tiles, so
  # the test is conditioned on the shape of the region actually processed.
  perm <- vapply(seq_len(n_perm), function(i) {
    idx <- rep(FALSE, nrow(sc)); idx[sample.int(nrow(sc), k)] <- TRUE
    count_pairs(idx)
  }, numeric(1))
  pval <- (1 + sum(perm >= obs)) / (1 + n_perm)
  msg("  %s: %d tiles, %d adjacent pairs (null %.1f, p = %.3f) -> %s",
      label, k, obs, mean(perm), pval,
      if (pval < 0.05) "CLUSTERED, consistent with a clonal population"
      else "scattered, consistent with per-tile noise")
  invisible(NULL)
}

msg("Tile size %d um, samples: %s", TILE_UM, paste(SAMPLES, collapse = ", "))
all_d <- list()
pdf(OUT_PDF, width = 9, height = 8)
for (s in SAMPLES) {
  msg("[%s]", s)
  d <- tryCatch(collect_sample(s), error = function(e) { msg("  %s", conditionMessage(e)); NULL })
  if (is.null(d)) next
  if (!any(!is.na(d$x_ratio))) { msg("  no tile has a chrX ratio - skipping plots"); next }
  print(panel_ratio_ocm(d, he = TRUE))
  print(panel_ratio_ocm(d, he = FALSE))
  print(panel_ratio(d, he = FALSE))
  print(panel_auto(d))
  print(panel_auto_zoom(d))
  print(panel_depth(d))
  print(panel_call(d))
  msg("  spatial clustering of the calls:")
  for (lab in c("CAST-skewed", "Bl6-skewed", "mixed")) clustering_test(d, lab)
  all_d[[s]] <- d
}
invisible(dev.off())

##### --------------------- distributions --------------------- #####
# The maps show WHERE; these show WHERE MOST TILES SIT, which a map cannot.
# Worth having because the OCM ramp puts four green/olive bins across 0.40-0.80,
# so a field of strongly Bl6-biased tiles reads as "quite green" and the eye
# under-reads the skew. Counts settle it.

panel_hist <- function(out) {
  # Same bins and colours as the map, so a bar and a tile of the same colour
  # mean the same thing. Counts, not density: the reader wants "how many tiles".
  h <- out[!is.na(x_bin), .N, by = .(sample, x_bin)]
  ggplot(h, aes(x_bin, N, fill = x_bin)) +
    geom_col(width = 0.85) +
    geom_text(aes(label = N), vjust = -0.35, size = 2.6, colour = "#52514e") +
    # 0.5 sits on the boundary between the 5th and 6th bin.
    geom_vline(xintercept = 5.5, linetype = 2, colour = "#52514e") +
    annotate("text", x = 5.5, y = Inf, label = "  0.5", hjust = 0, vjust = 1.4,
             size = 2.6, colour = "#52514e") +
    scale_fill_manual(values = OCM_COLORS, limits = OCM_LABELS, drop = FALSE,
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    facet_wrap(~ sample, ncol = 1, scales = "free_y") +
    labs(title = sprintf("Where the tiles sit: chrX allelic ratio per %d um tile", TILE_UM),
         subtitle = "Same bins and colours as the maps. Dashed line = 0.5.",
         x = "chrX allelic ratio (Bl6 / total)", y = "tiles",
         caption = paste("No legend: the bars are positioned by the quantity they are coloured by,",
                         "so the x axis already carries it.")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          plot.subtitle = element_text(size = 8, colour = "#52514e"),
          plot.caption  = element_text(size = 7, colour = "#52514e", hjust = 0))
}

panel_violin <- function(out) {
  # chrX against its own autosomal control on one axis. This is the figure that
  # makes the point in a single frame: the control is a spike on 0.5, chrX is
  # broad and displaced - and neither of those is visible in a map.
  long <- melt(out[!is.na(x_ratio)], id.vars = "sample",
               measure.vars = c("x_ratio", "a_ratio"),
               variable.name = "chrom_set", value.name = "ratio")
  long[, chrom_set := factor(fifelse(chrom_set == "x_ratio", "chrX", "autosomal"),
                             levels = c("chrX", "autosomal"))]
  long <- long[!is.na(ratio)]
  lab <- out[!is.na(x_ratio), .(
    n = .N,
    med = median(x_ratio),
    pooled = sum(x_a1) / sum(x_n),
    lo = mean(x_ratio < 0.5)), by = sample][order(sample)]

  ggplot(long, aes(sample, ratio, fill = chrom_set)) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "#52514e") +
    geom_violin(position = position_dodge(width = 0.8), width = 0.75,
                colour = NA, alpha = 0.85, scale = "width", trim = TRUE) +
    # The box carries the median and quartiles the violin only implies.
    #
    # group = interaction(...) is load-bearing. Setting fill = "white" as a
    # fixed parameter REMOVES the fill aesthetic, and with it the grouping the
    # violins get from fill = chrom_set - so geom_boxplot fell back to grouping
    # by x alone and drew ONE box per sample over chrX and autosomal pooled
    # together. It looked plausible (a median near 0.54) and was meaningless.
    geom_boxplot(aes(group = interaction(sample, chrom_set)),
                 position = position_dodge(width = 0.8), width = 0.12,
                 outlier.shape = NA, colour = "#0b0b0b",
                 fill = "white", alpha = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = c(chrX = "#2D6E5D", autosomal = "grey70"),
                      name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = "chrX against its own autosomal control, per tile",
         subtitle = paste0(
           paste(sprintf("%s: n=%d, median %.2f, read-weighted %.2f, %.0f%% below 0.5",
                         as.character(lab$sample), lab$n, lab$med, lab$pooled,
                         100 * lab$lo),
                 collapse = "   |   "),
           "\nread-weighted well below the median means the DEEP tiles are the low-ratio ones"),
         x = NULL, y = "allelic ratio (Bl6 / total)",
         caption = paste("Violins scaled to equal width, so shape is comparable between",
                         "groups of different size; the box is the median and quartiles.",
                         "\nThe autosomal spike on 0.5 is the null the chrX distribution",
                         "has to be judged against.")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8, colour = "#52514e"),
          plot.caption  = element_text(size = 7, colour = "#52514e", hjust = 0))
}

if (length(all_d)) {
  out <- rbindlist(all_d, fill = TRUE)
  # Order the panels by SAMPLES (9w, 78w), not alphabetically - "78w" sorts
  # before "9w" as a string, which puts the aged sample first and reads
  # backwards. intersect keeps the SAMPLES order while dropping any sample that
  # produced no tiles, so an absent one cannot leave an empty facet behind.
  out[, sample := factor(sample, levels = intersect(SAMPLES, unique(sample)))]
  pdf(sub("\\.pdf$", "_distribution.pdf", OUT_PDF), width = 7, height = 7)
  print(panel_hist(out))
  print(panel_violin(out))
  invisible(dev.off())
  msg("Wrote %s", sub("\\.pdf$", "_distribution.pdf", OUT_PDF))
  csv <- sub("\\.pdf$", ".csv", OUT_PDF)
  fwrite(out[, .(sample, tile, x, y, n_bins, x_a1, x_a2, x_n, a_a1, a_a2, a_n,
                 x_ratio, x_bin, a_ratio, z, call, submitted)], csv)
  msg("\nWrote %s\n       %s", OUT_PDF, csv)
  msg("Pooled chrX ratio per sample.")
  msg("NOTE: tiles are scored in tile-name order, so the finished set is a")
  msg("      contiguous band of the section, not a random sample of it. These")
  msg("      pooled numbers describe that band until the run completes.")
  print(out[!is.na(x_ratio), .(tiles = .N, chrX_reads = sum(x_n),
                               pooled_x = round(sum(x_a1) / sum(x_n), 3),
                               pooled_a = round(sum(a_a1) / sum(a_n), 3)),
                           by = sample][order(sample)])
} else {
  msg("Nothing to plot.")
}
