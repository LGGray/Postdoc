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

# Per-tile autosomal ratios are overdispersed relative to binomial (measured at
# 4.1x on the first 30 tiles), so the null band is +-0.03, not +-0.008. Folded
# into the z so a tile is only called when it clears the real spread.
AUTO_SD <- 0.032
Z_CALL  <- 3

# Diverging pair, neutral MID-POINT DELIBERATELY NEAR-WHITE rather than the
# usual grey: grey is spent on NA here, and a grey midpoint would make a
# genuinely biallelic tile indistinguishable from an unprocessed one.
COL_CAST <- "#184f95"   # ratio -> 0, CAST (A2) dominant
COL_MID  <- "#f0efec"   # 0.5
COL_BL6  <- "#b02a2a"   # ratio -> 1, Bl6 (A1) dominant
COL_NA   <- "#a8a8a4"   # submitted, not yet scored
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
  d[, `:=`(x_ratio = fifelse(x_n > 0, x_a1 / x_n, NA_real_),
           a_ratio = fifelse(a_n > 0, a_a1 / a_n, NA_real_))]
  d[, se := sqrt(a_ratio * (1 - a_ratio) / x_n + AUTO_SD^2)]
  d[, z := (x_ratio - a_ratio) / se]
  d[, call := fcase(
      is.na(x_ratio) &  submitted, "pending",
      is.na(x_ratio) & !submitted, "not submitted",
      z >  Z_CALL, "Bl6-skewed",
      z < -Z_CALL, "CAST-skewed",
      default = "mixed")]
  d[, `:=`(sample = smp, side = side, FLIPX = FLIP_X, FLIPY = FLIP_Y,
         he = he_png)]
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
         caption = "Any structure here is technical and invalidates the chrX panel over the same tiles.")
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
    labs(title = sprintf("%s - call vs each tile's own autosomal ratio", d$sample[1]),
         subtitle = sprintf("|z| > %d, with the %.3f autosomal overdispersion folded in: %s",
                            Z_CALL, AUTO_SD,
                            paste(sprintf("%s %d", lv[1:4],
                                          sapply(lv[1:4], function(l) sum(d2$call == l, na.rm = TRUE))),
                                  collapse = ", ")))
}

msg("Tile size %d um, samples: %s", TILE_UM, paste(SAMPLES, collapse = ", "))
all_d <- list()
pdf(OUT_PDF, width = 9, height = 8)
for (s in SAMPLES) {
  msg("[%s]", s)
  d <- tryCatch(collect_sample(s), error = function(e) { msg("  %s", conditionMessage(e)); NULL })
  if (is.null(d)) next
  if (!any(!is.na(d$x_ratio))) { msg("  no tile has a chrX ratio - skipping plots"); next }
  print(panel_ratio(d, he = TRUE))
  print(panel_ratio(d, he = FALSE))
  print(panel_auto(d))
  print(panel_depth(d))
  print(panel_call(d))
  all_d[[s]] <- d
}
invisible(dev.off())

if (length(all_d)) {
  out <- rbindlist(all_d, fill = TRUE)
  csv <- sub("\\.pdf$", ".csv", OUT_PDF)
  fwrite(out[, .(sample, tile, x, y, n_bins, x_a1, x_a2, x_n, a_a1, a_a2, a_n,
                 x_ratio, a_ratio, z, call, submitted)], csv)
  msg("\nWrote %s\n       %s", OUT_PDF, csv)
  msg("Pooled chrX ratio per sample:")
  print(out[!is.na(x_ratio), .(tiles = .N, chrX_reads = sum(x_n),
                               pooled_x = round(sum(x_a1) / sum(x_n), 3),
                               pooled_a = round(sum(a_a1) / sum(a_n), 3)), by = sample])
} else {
  msg("Nothing to plot.")
}
