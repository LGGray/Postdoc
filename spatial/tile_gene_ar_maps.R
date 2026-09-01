# Slide maps: every 64 um tile coloured by ONE GENE's allelic ratio, 9w beside 78w.
#
# tile_ratio_map.R maps chrX POOLED - every informative molecule on the
# chromosome, so nearly every tile is coloured and the map reads as a map. This
# maps one gene at a time, which is a different picture with a different failure
# mode, and the point of the layout is to make that failure mode visible rather
# than to hide it. Read the next paragraph before reading a map.
#
# WHAT A COLOURED TILE MEANS HERE. Outside Smpx, a 64 um tile holds about ONE
# informative molecule of any gene in this panel: at 9w, >= 2 UMIs is reached by
# 167/1047 Tspan7 tiles, 10/228 Kdm5c tiles and 0/44 Ftx tiles, and 97-100% of
# tiles read exactly 0 or 1 even on duplicate-inclusive reads, because every
# duplicate of the one molecule votes the same way. So a tile is not usually
# showing you a ratio - it is showing you which allele its single molecule came
# from. Dark red = that molecule was B6, dark blue = it was CAST. The map is a
# map of coin flips whose BIAS is the quantity of interest, and only the density
# of blue across a region means anything, never one blue tile.
#
# Coverage is the other half of it and gets its own row: for most of this panel
# only 0.3-5% of in-tissue tiles carry the gene at all (Ftx 44 of ~3950 tiles,
# Ddx3x 28, Shroom4 16), so an apparent pattern in the AR row is very often the
# expression pattern showing through. Compare the two rows before believing one.
#
# Run on the cluster:
#   conda activate seurat_env
#   Rscript ~/Postdoc/spatial/tile_gene_ar_maps.R
#   GENES=Kdm6a,Ftx,H19 Rscript ~/Postdoc/spatial/tile_gene_ar_maps.R
# or from slurm/spatial_tile_gene_ar.slurm, which runs it after the distributions.
#
# INPUT is the same tile x gene table the distributions use, plus the 2 um
# tissue_positions for geometry:
#   <IN_ROOT>/<sample>/tile_gene_counts.tsv.gz
#   <sample>/outs/binned_outputs/square_002um/spatial/tissue_positions.parquet
#
# WHICH TREE, AND WHAT "reads" MEANS IN IT. spatial_tile_locus_map.slurm has
# three counting modes and each writes its own OUT_ROOT:
#   ase_pysam_64um       umi    duplicates DROPPED. a1_reads = DEDUPLICATED reads
#                               (Allelome.PRO2's statistic), a1_umi = molecules.
#   ase_pysam_dup_64um   dup    duplicates KEPT. a1_reads = duplicate-INCLUSIVE
#                               reads, a1_umi = molecules.   <- the default here
#   ase_pysam_reads_64um reads  deduplicated reads only. NEVER RUN - no such tree.
# So LEVELS=reads against the default IN_ROOT is PCR-duplicate-inclusive, not
# deduplicated. For deduplicated reads point IN_ROOT at ase_pysam_64um; its
# read columns are exactly that, and no recount is needed.
#
# AND THE MOLECULE COLUMNS DIFFER BETWEEN THE TWO TREES - they are not the same
# number counted twice. chrX at 9w: 116,354 informative molecules at AR 0.8728
# in the umi tree, 160,901 at AR 0.8980 in the dup tree. Keeping duplicates
# finds 38% more MOLECULES (26% at 78w), because a (bin, UB) whose unmarked
# representative read misses a SNP is lost outright when duplicates are dropped
# while its marked duplicates would have covered one. That moves the answer:
# 10.2% escape here against the 12.7% the umi tree gives. Neither is wrong, but
# a number quoted from one tree cannot be compared with one from the other.
#
# LEVEL. The default is the duplicate-inclusive read level, because that is what
# the dup tree exists to show - but on THIS figure the two levels are nearly
# identical, and for a reason worth stating: a tile's colour is set by which
# allele its molecules came from, and duplicating a molecule cannot change that.
# The levels differ where a tile holds several molecules with a pile-up on one
# of them, which is Smpx and little else. The COVERAGE row is always drawn in
# MOLECULES whatever LEVEL says, because reads there would draw a PCR pile-up
# map rather than a coverage map.
#
# n = 1 animal per age. 9w beside 78w is a description of two sections, never an
# age effect.

.libPaths(c("~/R/matrix-dev", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

##### ---------------------- CONFIG ---------------------- #####
BASE    <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial"
IN_ROOT <- Sys.getenv("IN_ROOT", file.path(BASE, "ase_pysam_dup_64um"))
SAMPLES <- strsplit(Sys.getenv("SAMPLES", "9w,78w"), ",")[[1]]
TILE_UM <- as.integer(Sys.getenv("TILE_UM", "64"))
BIN     <- "square_002um"
LEVEL   <- Sys.getenv("LEVEL", "reads")          # reads (dup) | umi
OUT_DIR <- Sys.getenv("OUT_DIR", IN_ROOT)
OUT_PDF <- Sys.getenv("OUT_PDF",
                      file.path(OUT_DIR, sprintf("tile_gene_ar_maps_%dum_%s.pdf",
                                                 TILE_UM, LEVEL)))

# Restrict to a few genes, for a quick look: GENES=Kdm6a,Ftx
GENES   <- Sys.getenv("GENES", "")

# Draw the H&E underneath. Needs the png package; skipped with a message
# otherwise, exactly as tile_ratio_map.R does it.
USE_HE  <- !identical(tolower(Sys.getenv("USE_HE", "TRUE")), "false")

# Below this many coloured tiles a gene gets a page anyway, but the title says
# so - a map of nine tiles is a location, not a pattern.
MIN_TILES_MAP <- 25L

# The capture array can sit rotated relative to the image. Flip these if the map
# comes out mirrored against the H&E. Same meaning as in tile_ratio_map.R.
FLIP_X <- FALSE
FLIP_Y <- FALSE

# ---- the lab allelic-ratio ramp, verbatim from tile_ratio_map.R ------------
# Blue = A2/CAST-dominant, green ~ biallelic, olive/yellow -> orange -> dark red
# = increasingly monoallelic A1/Bl6. Binned rather than continuous, which is
# honest about per-tile precision - and on this figure especially, where almost
# every tile lands in the bottom or the top bin because it holds one molecule.
# Lightness is NOT monotone along it: the darkest bins are at BOTH ends, so read
# the legend, never the darkness.
OCM_BREAKS <- c(seq(0, 0.9, by = 0.1), 0.95, 1.0)
OCM_COLORS <- c("#2B3186", "#3B5FB6", "#38749F", "#367373", "#2D6E5D", "#1E652D",
                "#658C2D", "#8D9F25", "#B3B112", "#C97314", "#8B1913")
OCM_LABELS <- c("0.00-0.10", "0.10-0.20", "0.20-0.30", "0.30-0.40", "0.40-0.50",
                "0.50-0.60", "0.60-0.70", "0.70-0.80", "0.80-0.90",
                "0.90-0.95", "0.95-1.00")
names(OCM_COLORS) <- OCM_LABELS
COL_FOOT <- "#e6e5e1"   # in-tissue, but this gene has no informative molecule
##### ------------------------------------------------------ #####

msg <- function(...) message(sprintf(...))

source_panel <- function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else NULL
  for (d in c(here, "spatial", "~/Postdoc/spatial", ".")) {
    f <- file.path(d, "tile_gene_panel.R")
    if (file.exists(f)) { source(f, local = FALSE); return(f) }
  }
  stop("Cannot find tile_gene_panel.R next to this script")
}
msg("panel from %s", source_panel())

# Side-by-side needs a composer. patchwork comes with Seurat, so seurat_env has
# it; gridExtra is the fallback, and one-panel-per-page the fallback to that -
# a degraded layout beats a script that will not run.
# Braced deliberately. At TOP LEVEL R finishes the assignment at the end of the
# first line and then meets a bare `else`, which is a parse error - the same
# idiom reads fine inside a function body, which is why tile_ratio_map.R gets
# away with it and this did not.
COMPOSE <- {
  if (requireNamespace("patchwork", quietly = TRUE)) "patchwork"
  else if (requireNamespace("gridExtra", quietly = TRUE)) "gridExtra"
  else "pages"
}
if (COMPOSE == "pages")
  msg("neither patchwork nor gridExtra - drawing one panel per page")

##### ------------------- GEOMETRY ----------------------- #####
# read_positions() and lowres_scalef() are lifted VERBATIM from
# tile_ratio_map.R, and the tile naming below reproduces its
# `array_row %/% k` construction exactly, so a tile lands in the same place on
# both figures. Kept as a copy rather than a shared source file because
# tile_ratio_map.R is a settled figure script and this one should not be able to
# break it; if either changes, change both.
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

lowres_scalef <- function(dir) {
  f <- file.path(dir, "spatial", "scalefactors_json.json")
  if (!file.exists(f)) return(NA_real_)
  txt <- paste(readLines(f, warn = FALSE), collapse = " ")
  m <- regmatches(txt, regexpr('"tissue_lowres_scalef"\\s*:\\s*[0-9.eE+-]+', txt))
  if (!length(m)) return(NA_real_)
  as.numeric(sub('.*:\\s*', '', m))
}

# One row per in-tissue tile: its centroid in lowres pixel space, the tile side
# in those same units, and where the H&E for that sample lives. This is the
# FOOTPRINT - every tile the tissue covers, whether or not any given gene was
# seen in it - which is what stops an unexpressed gene from reading as no tissue.
tile_geometry <- function(smp) {
  bin_dir <- file.path(BASE, smp, "outs", "binned_outputs", BIN)
  pos <- read_positions(bin_dir)
  if (!all(c("array_row", "array_col") %in% names(pos)))
    stop("tissue_positions for ", smp, " has no array_row/array_col")
  it <- grep("^in_tissue", names(pos), value = TRUE)[1]
  if (!is.na(it)) pos <- pos[get(it) == 1]
  # arrow can hand back integer64, which %/% and matrix indexing mishandle.
  pos[, `:=`(array_row = as.numeric(array_row), array_col = as.numeric(array_col),
             px = as.numeric(pxl_col_in_fullres), py = as.numeric(pxl_row_in_fullres))]

  k <- TILE_UM / 2                      # bins per tile side
  pos[, tile := sprintf("tile_%dum_r%04d_c%04d", TILE_UM,
                        as.integer(array_row %/% k), as.integer(array_col %/% k))]

  sf <- lowres_scalef(bin_dir)
  if (!is.finite(sf)) { sf <- 1; msg("  %s: no lowres scalefactor - full-res pixels", smp) }
  geom <- pos[, .(x = mean(px) * sf, y = mean(py) * sf, n_bins = .N), by = tile]

  # Tile side in the plotted units, from the array->pixel slopes, so it is right
  # even where a tile is only partly covered by tissue.
  sx <- coef(lm(x ~ array_col + array_row, data = pos[, .(x = px * sf, array_col, array_row)]))
  sy <- coef(lm(y ~ array_col + array_row, data = pos[, .(y = py * sf, array_col, array_row)]))
  side <- k * mean(c(sqrt(sx[["array_col"]]^2 + sy[["array_col"]]^2),
                     sqrt(sx[["array_row"]]^2 + sy[["array_row"]]^2)))
  msg("  %s: %d in-tissue tiles, side %.1f plot units, lowres scalef %.4f",
      smp, nrow(geom), side, sf)

  geom[, `:=`(sample = smp, side = side,
              he = if (USE_HE) file.path(bin_dir, "spatial", "tissue_lowres_image.png")
                   else NA_character_)]
  if (FLIP_X) geom[, x := -x]
  if (FLIP_Y) geom[, y := -y]
  geom[]
}

##### --------------------- COUNTS ----------------------- #####
read_genes <- function(s, keep) {
  f <- file.path(IN_ROOT, s, "tile_gene_counts.tsv")
  if (!file.exists(f) && file.exists(paste0(f, ".gz"))) f <- paste0(f, ".gz")
  if (!file.exists(f)) stop("No tile_gene_counts.tsv[.gz] under ", file.path(IN_ROOT, s))
  d <- fread(f, select = c("tile", "chrom", "gene", "a1_umi", "a2_umi", "n_umi",
                           "a1_reads", "a2_reads", "n_reads"), showProgress = FALSE)
  d[gene %chin% keep][, sample := s][]
}
read_chrom <- function(s) {
  f <- file.path(IN_ROOT, s, "tile_chrom_counts.tsv")
  if (!file.exists(f)) stop("No tile_chrom_counts.tsv under ", file.path(IN_ROOT, s))
  fread(f, showProgress = FALSE)[chrom == "chrX",
    .(tile, chrom, gene = "chrX (all genes pooled)", a1_umi, a2_umi, n_umi,
      a1_reads, a2_reads, n_reads)][, sample := s][]
}

want <- if (nzchar(GENES)) trimws(strsplit(GENES, ",")[[1]]) else PANEL_GENES()
msg("Reading %s", IN_ROOT)
cnt <- rbindlist(c(lapply(SAMPLES, read_genes, keep = want),
                   if ("chrX (all genes pooled)" %chin% want) lapply(SAMPLES, read_chrom)),
                 use.names = TRUE)
if (!nrow(cnt)) stop("None of the requested genes are in the counts: ",
                     paste(want, collapse = ", "))

missing <- setdiff(want, cnt[, unique(gene)])
if (length(missing))
  msg("no informative molecule anywhere, so no map: %s", paste(missing, collapse = ", "))

msg("Building tile geometry")
geo <- rbindlist(lapply(SAMPLES, tile_geometry))

SLAB <- c("9w" = "adult (9w)", "78w" = "aged (78w)")
meta <- PANEL_META()

##### ---------------------- DRAW ------------------------ #####
theme_slide <- theme_void(base_size = 10) +
  theme(legend.position = "right",
        legend.key.height = unit(11, "pt"),
        legend.key.width  = unit(8, "pt"),
        plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "#52514e"))

he_layer <- function(d) {
  f <- d$he[1]
  if (is.na(f) || !file.exists(f)) return(NULL)
  if (!requireNamespace("png", quietly = TRUE)) return(NULL)
  img <- png::readPNG(f)
  # The lowres PNG is the full-res image scaled by tissue_lowres_scalef, which is
  # exactly the space the tile centroids were put into, so the raster spans
  # 0..ncol / 0..nrow of the PNG itself.
  w <- dim(img)[2]; h <- dim(img)[1]
  xr <- if (FLIP_X) c(-w, 0) else c(0, w)
  yr <- if (FLIP_Y) c(-h, 0) else c(0, h)
  annotation_raster(img, xmin = xr[1], xmax = xr[2], ymin = yr[1], ymax = yr[2],
                    interpolate = TRUE)
}

# y is reversed everywhere: image coordinates run downwards from the top-left.
base_map <- function(d, he = USE_HE) {
  g <- ggplot(d, aes(x, y))
  if (he) { l <- he_layer(d); if (!is.null(l)) g <- g + l }
  g +
    geom_tile(width = d$side[1], height = d$side[1], fill = COL_FOOT, colour = NA) +
    coord_fixed() + scale_y_reverse() + theme_slide
}

panel_ar <- function(d, g_, smp) {
  hit <- d[!is.na(ar)]
  base_map(d) +
    geom_tile(data = hit, aes(fill = ar_bin),
              width = d$side[1], height = d$side[1], colour = NA) +
    scale_fill_manual(values = OCM_COLORS, limits = OCM_LABELS, drop = FALSE,
                      na.value = COL_FOOT, name = "AR\n(B6 / total)") +
    labs(title = sprintf("%s - %s", SLAB[[smp]], g_),
         subtitle = sprintf("%d of %d in-tissue tiles carry it (%.1f%%); %.0f%% of those hold one molecule, so their colour is one allele call",
                            nrow(hit), nrow(d), 100 * nrow(hit) / nrow(d),
                            100 * mean(hit$mono)))
}

panel_cov <- function(d, g_, smp) {
  hit <- d[!is.na(ar)]
  base_map(d) +
    geom_tile(data = hit, aes(fill = n_umi),
              width = d$side[1], height = d$side[1], colour = NA) +
    scale_fill_gradient(low = "#cde2fb", high = "#0d366b", trans = "sqrt",
                        na.value = COL_FOOT, name = "molecules") +
    labs(title = sprintf("%s - %s coverage", SLAB[[smp]], g_),
         subtitle = sprintf("informative MOLECULES per tile (never reads: duplicates would draw a PCR pile-up map). median %g, max %d",
                            median(hit$n_umi), as.integer(max(hit$n_umi))))
}

# One gene, one page: AR across the top, coverage underneath, ages left to right.
one_gene <- function(g_) {
  a1 <- paste0("a1_", LEVEL); nn <- paste0("n_", LEVEL)
  panes <- list(); cov <- list(); n_hit <- 0L
  for (s in SAMPLES) {
    g <- geo[sample == s]
    c1 <- cnt[sample == s & gene == g_,
              .(tile, ar = fifelse(get(nn) > 0, get(a1) / get(nn), NA_real_),
                n_umi, mono = n_umi > 0 & (a1_umi == 0 | a2_umi == 0))]
    d <- merge(g, c1, by = "tile", all.x = TRUE)
    d[, ar_bin := cut(ar, breaks = OCM_BREAKS, include.lowest = TRUE,
                      right = TRUE, labels = OCM_LABELS)]
    n_hit <- n_hit + sum(!is.na(d$ar))
    if (!sum(!is.na(d$ar))) {   # no data for this gene in this section
      panes[[s]] <- base_map(d) +
        labs(title = sprintf("%s - %s", SLAB[[s]], g_),
             subtitle = "no informative molecule in this section")
      cov[[s]] <- panes[[s]]
      next
    }
    panes[[s]] <- panel_ar(d, g_, s)
    cov[[s]]   <- panel_cov(d, g_, s)
  }
  list(panes = unname(c(panes, cov)), n_hit = n_hit)
}

msg("Drawing %d genes -> %s", length(unique(cnt$gene)), OUT_PDF)
dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)
pdf(OUT_PDF, width = 11, height = 9)
pages <- 0L
for (g_ in PANEL_GENES(have = cnt[, unique(gene)])) {
  r <- one_gene(g_)
  ex <- meta[gene == g_, expect]
  tag <- if (length(ex) && !is.na(ex[1]))
           sprintf("  |  imprinted control, expected AR = %g (%s)", ex[1],
                   if (ex[1] == 1) "maternal = B6" else "paternal = CAST") else ""
  thin <- if (r$n_hit < MIN_TILES_MAP)
            sprintf("  |  ONLY %d COLOURED TILES IN TOTAL - a location, not a pattern",
                    r$n_hit) else ""
  cap <- paste0(g_, tag, thin)
  if (COMPOSE == "patchwork") {
    print(patchwork::wrap_plots(r$panes, ncol = length(SAMPLES)) +
          patchwork::plot_annotation(
            title = cap,
            caption = "Top row: allelic ratio, dark red = B6 (active X), dark blue = CAST (inactive X). Bottom row: how many informative molecules were behind it.\nPale tiles are in-tissue but carry no informative molecule of this gene. n = 1 animal per age.",
            theme = theme(plot.title = element_text(face = "bold", size = 12),
                          plot.caption = element_text(size = 7, colour = "#52514e", hjust = 0))))
  } else if (COMPOSE == "gridExtra") {
    gridExtra::grid.arrange(grobs = r$panes, ncol = length(SAMPLES), top = cap)
  } else {
    for (p in r$panes) print(p)
  }
  pages <- pages + 1L
}
invisible(dev.off())
msg("Wrote %s (%d genes)", OUT_PDF, pages)
