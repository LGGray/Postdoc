# Slide map: every 64 um tile coloured by the allelic ratio of the ELEVEN
# core escape genes POOLED, 9w beside 78w, with the depth behind each tile
# directly underneath it.
#
# WHY POOLED AND NOT PER GENE. tile_gene_ar_maps.R draws one gene per page, and
# for this panel that is 11 pages on which almost every coloured tile holds a
# single molecule (Ftx 44 tiles of ~3950 at 9w, Ddx3x 28). Pooling the panel is
# the only way to get more than one molecule into a tile at all, and it is the
# right pooling for the question being asked of it: LOSS OF THE INACTIVE X.
# These genes are the only chrX genes that report the Xi, so if a tile has lost
# its Xi they all go monoallelic B6 together, and their molecules are
# interchangeable evidence for that one event.
#
# WHAT A COLOURED TILE MEANS, AND READ THIS BEFORE READING THE MAP.
# Even pooled, the panel yields a MEDIAN OF ONE informative molecule per tile:
# 1070 of 3947 in-tissue tiles carry any panel molecule at 9w (799 of 4308 at
# 78w), 80% of those hold exactly one (86% at 78w), and the deepest tile on
# either section holds five. So the top row is mostly NOT showing a ratio - it
# is showing which allele a tile's single molecule came from. Dark red = that
# molecule was B6, dark blue = it was CAST. Only the DENSITY of blue over a
# region carries information; one dark red tile is one coin flip.
#
# WHAT THE MAP CANNOT DO, stated here because the figure invites the inference.
# "Monoallelic for escape genes" means zero CAST molecules seen. Under ordinary
# escape that happens by sampling alone with probability AR^n, and the pooled
# panel AR is 0.714 (9w) / 0.762 (78w), so a 5%-level monoallelic call needs
# n >= 9 molecules at 9w and n >= 12 at 78w. NO TILE ON EITHER SECTION REACHES
# EITHER NUMBER. At depth >= 3, 14 of 34 tiles (9w) and 17 of 21 (78w) carry no
# CAST molecule, against 12.4 and 16.0 expected under the per-gene null - i.e.
# exactly the sampling floor. Coarsening 64 um tiles to 256/320/512 um super-
# tiles (which needs no recount, trow/tcol are in the table) leaves 0 units at
# BH q < 0.05 in either section. This figure is therefore a description of where
# the panel's molecules landed and which allele they carried. It is not, and at
# this depth cannot be, a map of X loss.
#
# THE PANEL IS NOT AR-HOMOGENEOUS, which matters for the pooled ratio. Ftx, Jpx
# and 5530601H04Rik are Xic-linked: with B6 carrying the Xist deletion they are
# expressed FROM the inactive CAST X rather than escaping from it, and they sit
# at AR 0.64 / 0.32-0.59 / 0.17-0.25 against Pbdc1 at 0.87. Pooling them pulls
# the panel AR down, which makes a monoallelic threshold look easier to clear
# than it is. XIC_GENES below names them and DROP_XIC=TRUE removes them, which
# is the conservative panel; the default keeps them so the figure matches
# new_core_escape_genes in OCM_heart/core_escape_SNPs.R as written.
#
# Sts contributes essentially nothing (0 informative molecules at 9w, 7 at 78w)
# and is kept only so the gene list matches the source of truth.
#
# COUNTING LEVEL is MOLECULES, and unlike tile_gene_ar_maps.R that is not
# overridable-by-default here. The tree is the --keep-duplicates one, where the
# read columns are duplicate-inclusive AND the duplication is allele-asymmetric
# (B6 ~17x against CAST ~10x on chrX at 9w, per chrX_escape_qc.R), so a
# read-level pooled ratio biases toward monoallelic B6 - the exact artefact this
# figure would be read as evidence of. LEVEL=reads is accepted and prints a
# warning on the console and on the page.
#
# n = 1 animal per age. 9w beside 78w is a description of two sections, never an
# age effect.
#
# Run on the cluster:
#   conda activate seurat_env
#   Rscript ~/Postdoc/spatial/tile_escape_panel_map.R
#   MIN_UMI=2 Rscript ~/Postdoc/spatial/tile_escape_panel_map.R   # drop singletons
#   DROP_XIC=TRUE Rscript ~/Postdoc/spatial/tile_escape_panel_map.R
#
# INPUT (read-only; nothing is recounted - the pysam run already carries every
# gene, because it was scored against the genome-wide annotation_us.bed):
#   <IN_ROOT>/<sample>/tile_gene_counts.tsv.gz
#   <sample>/outs/binned_outputs/square_002um/spatial/tissue_positions.parquet
#
# read_positions(), lowres_scalef(), tile_geometry(), the OCM ramp and the H&E
# layer are copies from tile_gene_ar_maps.R, which took them from
# tile_ratio_map.R. Copied rather than shared for the reason given there: those
# are settled figure scripts and a new one must not be able to break them. A
# tile lands in the same place on all three figures. If either changes, change
# all three.

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
LEVEL   <- Sys.getenv("LEVEL", "umi")            # umi (default) | reads
OUT_DIR <- Sys.getenv("OUT_DIR", IN_ROOT)

# new_core_escape_genes, verbatim from OCM_heart/core_escape_SNPs.R.
ESCAPE_GENES <- c("Kdm5c", "Kdm6a", "Ddx3x", "Eif2s3x", "Ftx",
                  "5530601H04Rik", "Jpx", "Pbdc1", "Utp14a",
                  "Akap17a", "Sts")

# The Xic-linked members: expressed from the Xi in this genotype rather than
# escaping from it. See the header.
XIC_GENES <- c("Ftx", "Jpx", "5530601H04Rik")
DROP_XIC  <- identical(toupper(Sys.getenv("DROP_XIC", "FALSE")), "TRUE")

# Minimum pooled molecules before a tile is coloured in the AR row. 1 shows
# every tile that has anything, which is the honest default given that the
# depth row sits underneath saying how thin that is. 2 or 3 gives a sparser
# map in which a colour is at least a ratio.
MIN_UMI <- as.integer(Sys.getenv("MIN_UMI", "1"))

USE_HE  <- !identical(tolower(Sys.getenv("USE_HE", "TRUE")), "false")

# Flip if the map comes out mirrored against the H&E. Same meaning as in
# tile_ratio_map.R and tile_gene_ar_maps.R.
FLIP_X <- FALSE
FLIP_Y <- FALSE

TAG <- paste0(if (DROP_XIC) "noXic" else "full", "_", LEVEL,
              if (MIN_UMI > 1) paste0("_min", MIN_UMI) else "")
OUT_PDF <- Sys.getenv("OUT_PDF",
                      file.path(OUT_DIR, sprintf("tile_escape_panel_map_%dum_%s.pdf",
                                                 TILE_UM, TAG)))
OUT_CSV <- Sys.getenv("OUT_CSV", sub("\\.pdf$", ".csv.gz", OUT_PDF))

# ---- the lab allelic-ratio ramp, verbatim from tile_ratio_map.R ------------
# Blue = A2/CAST-dominant, green ~ biallelic, olive/yellow -> orange -> dark red
# = increasingly monoallelic A1/Bl6. Lightness is NOT monotone along it: the
# darkest bins are at BOTH ends, so read the legend, never the darkness.
OCM_BREAKS <- c(seq(0, 0.9, by = 0.1), 0.95, 1.0)
OCM_COLORS <- c("#2B3186", "#3B5FB6", "#38749F", "#367373", "#2D6E5D", "#1E652D",
                "#658C2D", "#8D9F25", "#B3B112", "#C97314", "#8B1913")
OCM_LABELS <- c("0.00-0.10", "0.10-0.20", "0.20-0.30", "0.30-0.40", "0.40-0.50",
                "0.50-0.60", "0.60-0.70", "0.70-0.80", "0.80-0.90",
                "0.90-0.95", "0.95-1.00")
names(OCM_COLORS) <- OCM_LABELS
COL_FOOT   <- "#e6e5e1"    # in-tissue, but no informative panel molecule
FOOT_ALPHA <- if (USE_HE) 0.55 else 1
##### ------------------------------------------------------ #####

msg <- function(...) message(sprintf(...))

if (!LEVEL %in% c("umi", "reads")) stop("LEVEL must be umi or reads, got ", LEVEL)
# Braced deliberately: unbraced, only the FIRST msg() is conditional and the
# other two print on every run - the same trap tile_gene_ar_maps.R documents for
# a top-level if/else.
if (LEVEL == "reads") {
  msg("WARNING: LEVEL=reads on a --keep-duplicates tree is duplicate-inclusive and")
  msg("         allele-asymmetric (B6 ~17x vs CAST ~10x). The pooled ratio will be")
  msg("         biased toward monoallelic B6. Use LEVEL=umi for any claim.")
}

genes <- if (DROP_XIC) setdiff(ESCAPE_GENES, XIC_GENES) else ESCAPE_GENES
msg("Panel: %d genes (%s)", length(genes), paste(genes, collapse = ", "))

COMPOSE <- {
  if (requireNamespace("patchwork", quietly = TRUE)) "patchwork"
  else if (requireNamespace("gridExtra", quietly = TRUE)) "gridExtra"
  else "pages"
}
if (COMPOSE == "pages")
  msg("neither patchwork nor gridExtra - drawing one panel per page")

##### ------------------- GEOMETRY ----------------------- #####
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

# One row per in-tissue tile: centroid in lowres pixel space, tile side in the
# same units, and where the H&E lives. This is the FOOTPRINT - every tile the
# tissue covers, whether or not the panel was seen in it - which is what stops
# a sparse panel from reading as absent tissue.
tile_geometry <- function(smp) {
  bin_dir <- file.path(BASE, smp, "outs", "binned_outputs", BIN)
  pos <- read_positions(bin_dir)
  if (!all(c("array_row", "array_col") %in% names(pos)))
    stop("tissue_positions for ", smp, " has no array_row/array_col")
  it <- grep("^in_tissue", names(pos), value = TRUE)[1]
  if (!is.na(it)) pos <- pos[get(it) == 1]
  pos[, `:=`(array_row = as.numeric(array_row), array_col = as.numeric(array_col),
             px = as.numeric(pxl_col_in_fullres), py = as.numeric(pxl_row_in_fullres))]

  k <- TILE_UM / 2                      # 2 um bins per tile side
  pos[, tile := sprintf("tile_%dum_r%04d_c%04d", TILE_UM,
                        as.integer(array_row %/% k), as.integer(array_col %/% k))]

  sf <- lowres_scalef(bin_dir)
  if (!is.finite(sf)) { sf <- 1; msg("  %s: no lowres scalefactor - full-res pixels", smp) }
  geom <- pos[, .(x = mean(px) * sf, y = mean(py) * sf, n_bins = .N), by = tile]

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
# Pool the panel WITHIN a tile: sum a1 and a2 over the genes, so one tile is one
# row. n_umi = a1_umi + a2_umi holds in this table by construction (ties are
# dropped by the counter's majority vote, never split), verified over the first
# 800k rows of the 9w file, so summing the allele columns and summing n agree.
read_panel <- function(s) {
  f <- file.path(IN_ROOT, s, "tile_gene_counts.tsv")
  if (!file.exists(f) && file.exists(paste0(f, ".gz"))) f <- paste0(f, ".gz")
  if (!file.exists(f)) stop("No tile_gene_counts.tsv[.gz] under ", file.path(IN_ROOT, s))
  d <- fread(f, select = c("tile", "chrom", "gene",
                           "a1_umi", "a2_umi", "n_umi",
                           "a1_reads", "a2_reads", "n_reads"),
             showProgress = FALSE)
  # chrom is checked as well as gene: a panel name could in principle collide
  # with an autosomal gene of the same name in the annotation, and an autosomal
  # molecule in an Xi-loss panel would be silent nonsense.
  d <- d[chrom == "chrX" & gene %chin% genes]
  seen <- d[n_umi > 0, unique(gene)]
  drop <- setdiff(genes, seen)
  if (length(drop))
    msg("  %s: no informative molecule anywhere for %s", s, paste(drop, collapse = ", "))
  out <- d[, .(a1_umi = sum(a1_umi), a2_umi = sum(a2_umi), n_umi = sum(n_umi),
               a1_reads = sum(a1_reads), a2_reads = sum(a2_reads),
               n_reads = sum(n_reads), n_genes = sum(n_umi > 0)), by = tile]
  out[, sample := s][]
}

msg("Reading %s", IN_ROOT)
cnt <- rbindlist(lapply(SAMPLES, read_panel), use.names = TRUE)
if (!nrow(cnt)) stop("No panel molecules in ", IN_ROOT)

msg("Building tile geometry")
geo <- rbindlist(lapply(SAMPLES, tile_geometry))

a1c <- paste0("a1_", LEVEL); nc <- paste0("n_", LEVEL)
D <- list()
for (s in SAMPLES) {
  c1 <- cnt[sample == s, .(tile, a1_umi, a2_umi, n_umi, a1_reads, a2_reads, n_reads,
                           n_genes,
                           depth = get(nc),
                           ar = fifelse(get(nc) > 0, get(a1c) / get(nc), NA_real_))]
  d <- merge(geo[sample == s], c1, by = "tile", all.x = TRUE)
  # A tile below MIN_UMI keeps its depth (the lower row still shows it) and
  # loses only its colour in the AR row, so the two rows are not on different
  # tile sets and the top row cannot imply coverage the bottom row denies.
  d[, ar_shown := fifelse(!is.na(depth) & depth >= MIN_UMI, ar, NA_real_)]
  d[, ar_bin := cut(ar_shown, breaks = OCM_BREAKS, include.lowest = TRUE,
                    right = TRUE, labels = OCM_LABELS)]
  D[[s]] <- d
}

# One coverage scale for both ages, computed over both. Left free, ggplot would
# scale each panel to its own maximum and equal blue would mean unequal depth -
# which would invalidate the side-by-side comparison the layout exists for.
COV_LIM <- c(0, max(unlist(lapply(D, function(d) suppressWarnings(max(d$depth, na.rm = TRUE)))),
                    na.rm = TRUE))
msg("Shared depth scale: 0-%g %s", COV_LIM[2], if (LEVEL == "umi") "molecules" else "reads")

##### ---------------------- STATS ----------------------- #####
# Printed, not just plotted: the numbers a reader needs to size the figure's
# claim, on the data actually loaded rather than from the header.
for (s in SAMPLES) {
  d <- D[[s]]; hit <- d[!is.na(depth) & depth > 0]
  pooled_ar <- sum(hit$a1_umi) / sum(hit$n_umi)
  # Capped: at pooled_ar == 1 - a section in which no panel molecule was CAST -
  # the unguarded loop never terminates, and that is a real possible input, not
  # a hypothetical one.
  n_need <- { n <- 1L; while (pooled_ar^n >= 0.05 && n < 1000L) n <- n + 1L; n }
  msg("%s: %d of %d in-tissue tiles carry a panel molecule (%.1f%%); depth median %g, max %g",
      s, nrow(hit), nrow(d), 100 * nrow(hit) / nrow(d),
      median(hit$depth), max(hit$depth))
  msg("    pooled panel AR (molecules) %.4f over %d molecules; %.0f%% of tiles hold one",
      pooled_ar, sum(hit$n_umi), 100 * mean(hit$n_umi == 1))
  msg("    a 5%%-level monoallelic call needs %d molecules; %d tiles reach it",
      n_need, sum(hit$depth >= n_need))
}

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
    geom_tile(width = d$side[1], height = d$side[1], fill = COL_FOOT,
              colour = NA, alpha = FOOT_ALPHA) +
    coord_fixed() + scale_y_reverse() + theme_slide
}

# Falls back to the sample name itself: SAMPLES is settable, and an unknown
# name should retitle the panel, not abort the figure on a [[ ]] subscript.
SLAB_MAP <- c("9w" = "adult (9w)", "78w" = "aged (78w)")
SLAB <- function(s) if (s %in% names(SLAB_MAP)) SLAB_MAP[[s]] else s

panel_ar <- function(d, smp) {
  hit <- d[!is.na(ar_shown)]
  base_map(d) +
    geom_tile(data = hit, aes(fill = ar_bin),
              width = d$side[1], height = d$side[1], colour = NA) +
    scale_fill_manual(values = OCM_COLORS, limits = OCM_LABELS, drop = FALSE,
                      na.value = COL_FOOT, name = "AR\n(B6 / total)") +
    labs(title = SLAB(smp),
         subtitle = sprintf("%d of %d tiles (%.1f%%)  |  pooled AR %.3f  |  %.0f%% one molecule",
                            nrow(hit), nrow(d), 100 * nrow(hit) / nrow(d),
                            sum(hit$a1_umi) / sum(hit$n_umi),
                            100 * mean(hit$n_umi == 1)))
}

panel_cov <- function(d, smp, lim) {
  hit <- d[!is.na(depth) & depth > 0]
  base_map(d) +
    geom_tile(data = hit, aes(fill = depth),
              width = d$side[1], height = d$side[1], colour = NA) +
    scale_fill_gradient(low = "#cde2fb", high = "#0d366b", trans = "sqrt",
                        limits = lim, na.value = COL_FOOT,
                        name = if (LEVEL == "umi") "molecules" else "reads") +
    labs(title = paste(SLAB(smp), "- depth"),
         subtitle = sprintf("%s per tile: median %g, 90th pct %g, max %g",
                            if (LEVEL == "umi") "molecules" else "reads",
                            median(hit$depth),
                            as.numeric(quantile(hit$depth, 0.9)), max(hit$depth)))
}

panes <- c(lapply(SAMPLES, function(s) panel_ar(D[[s]], s)),
           lapply(SAMPLES, function(s) panel_cov(D[[s]], s, COV_LIM)))

cap <- sprintf("Core escape panel pooled per %d um tile - %d genes%s, %s level",
               TILE_UM, length(genes),
               if (DROP_XIC) " (Xic-linked members dropped)" else "", LEVEL)

foot <- paste0(
  "Top row: allelic ratio of ", paste(genes, collapse = ", "), " pooled within each tile. ",
  "Dark red = B6 (active X), dark blue = CAST (inactive X); pale tiles are in-tissue with no informative panel molecule",
  if (MIN_UMI > 1) sprintf(" or fewer than %d", MIN_UMI) else "", ".\n",
  "Bottom row: ", if (LEVEL == "umi") "informative molecules" else "reads",
  " behind each tile, same scale for both ages, so equal blue is equal depth. Read it before reading the row above it.\n",
  "MOST TILES HOLD ONE MOLECULE, so their colour is which allele that molecule came from, not a ratio. Only the density of blue over a region means anything.\n",
  "A 5%-level monoallelic (Xi-loss) call needs ~9 molecules at 9w and ~12 at 78w given the pooled panel AR; no tile on either section reaches that, and no super-tile up to 512 um survives BH.\n",
  "Genotype: B6 mother x CAST father, B6 carries the Xist deletion, so CAST is the inactive X in every cell. n = 1 animal per age - 9w beside 78w is two sections, not an age effect.",
  if (LEVEL == "reads")
    "\nLEVEL=reads on a duplicate-inclusive tree: duplication is allele-asymmetric (B6 ~17x vs CAST ~10x), so this ratio is biased toward monoallelic B6. Not for any claim."
  else "")

msg("Drawing -> %s", OUT_PDF)
dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)
pdf(OUT_PDF, width = 11, height = 9)
if (COMPOSE == "patchwork") {
  print(patchwork::wrap_plots(panes, ncol = length(SAMPLES), guides = "collect") +
        patchwork::plot_annotation(
          title = cap, caption = foot,
          theme = theme(plot.title = element_text(face = "bold", size = 12),
                        plot.caption = element_text(size = 7, colour = "#52514e", hjust = 0))))
} else if (COMPOSE == "gridExtra") {
  gridExtra::grid.arrange(grobs = panes, ncol = length(SAMPLES), top = cap)
} else {
  for (p in panes) print(p)
}
invisible(dev.off())

# The per-tile table behind the figure, so a number quoted off the map can be
# traced to the tile it came from without re-reading the 4.7 M-row gene table.
out <- rbindlist(lapply(SAMPLES, function(s)
  D[[s]][, .(sample, tile, x, y, n_bins, n_genes,
             a1_umi, a2_umi, n_umi, a1_reads, a2_reads, n_reads,
             ar_umi = fifelse(n_umi > 0, a1_umi / n_umi, NA_real_),
             ar_reads = fifelse(n_reads > 0, a1_reads / n_reads, NA_real_),
             plotted = !is.na(ar_shown))]), use.names = TRUE)
fwrite(out, OUT_CSV)
msg("Wrote %s and %s", OUT_PDF, OUT_CSV)
