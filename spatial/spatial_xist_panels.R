# Xist expression on the Visium HD cardiac sections, 9w beside 78w.
#
# WHAT THIS IS FOR. A quick descriptive look at where Xist is and how much of
# it there is, reusing the annotated objects rather than re-reading the
# spaceranger matrices. Nothing here is allelic: this is total Xist UMIs per
# bin, normalised for depth.
#
# READ THE GENOTYPE BEFORE READING THE NUMBERS. B6 mother x CAST father with
# Xist deleted on the B6 allele, so every Xist molecule in these sections comes
# from the CAST X - the inactive one. Total Xist here is already single-allele
# by construction, and a 9w vs 78w difference is a difference in CAST Xist
# output, not a shift in which allele is expressed.
#
# WHAT IT WRITES, in <SPATIAL_DIR>/annotation/<GENE>_<BIN>/:
#   fig1_<gene>_map            tissue map per sample, shared colour scale
#   fig2_<gene>_by_celltype    detection rate and depth-normalised level per
#                              cell type, both ages side by side
#   fig3_<gene>_level_violin   level among positive bins only, per cell type -
#                              i.e. how much, given any at all
#   <gene>_summary.csv         the numbers behind fig2, per sample x cell type
#
# Depth normalisation is POOLED, not a mean of per-bin ratios: sum(gene UMIs) /
# sum(total UMIs) x 10^4 within each group. At 8um most bins carry a handful of
# UMIs, so per-bin ratios are dominated by the shallow bins and a mean over
# them is not an expression estimate. The log-normalised per-bin value is still
# used for the maps and the violin, where the point is the spatial pattern and
# the spread rather than a group estimate.
#
# Bins are not cells (see spatial_cell_annotation.R): a cell-type label means
# "dominated by", so a per-type Xist level is a neighbourhood average.
#
# At 8um most positive bins carry exactly one gene UMI, so fig1 is closer to a
# detection map than to a level map, and detection is partly a depth map. fig2's
# lower panel and fig3 are the panels to read a level off.
#
# Run on the cluster:
#   sbatch slurm/spatial_xist_panels.slurm            # Xist, 9w and 78w
#   sbatch slurm/spatial_xist_panels.slurm 9w         # one sample
#   sbatch --export=ALL,GENE=Tsix slurm/spatial_xist_panels.slurm
# or by hand from adult_aged_spatial/:
#   GENE=Xist Rscript ~/Postdoc/spatial/spatial_xist_panels.R

# Must be line 1 of anything that touches Seurat/Matrix on this cluster (see
# banksy_patch.R for why).
.libPaths(c("~/R/matrix-dev", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})
for (p in c("data.table", "patchwork")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("R package '", p, "' is required")
}
library(data.table)
library(patchwork)

##### ---------------------- CONFIG ---------------------- #####
SPATIAL_DIR <- Sys.getenv("SPATIAL_DIR",
  "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial")
BASE    <- SPATIAL_DIR
BIN     <- Sys.getenv("BIN", "square_008um")
GENE    <- Sys.getenv("GENE", "Xist")
SAMPLES <- trimws(strsplit(Sys.getenv("SAMPLES", "9w,78w"), "[,[:space:]]+")[[1]])
SAMPLES <- SAMPLES[nzchar(SAMPLES)]
if (!length(SAMPLES)) stop("SAMPLES is empty")

OUT_DIR <- Sys.getenv("OUT_DIR", "")
if (!nzchar(OUT_DIR)) OUT_DIR <- file.path(BASE, "annotation", paste0(GENE, "_", BIN))

POINT_SIZE <- as.numeric(Sys.getenv("POINT_SIZE", "0.35"))
if (!is.finite(POINT_SIZE) || POINT_SIZE <= 0) stop("POINT_SIZE must be a positive number")

# Upper end of the map colour scale, as a quantile of the POSITIVE bins. A few
# bins carry an order of magnitude more Xist than the rest and stretching the
# scale to them flattens everything else to one colour.
CAP_Q <- as.numeric(Sys.getenv("CAP_Q", "0.99"))

# Violin panel only: bins per sample x cell type, sampled for plotting speed.
# Densities from 20k points are visually identical to densities from 150k.
VIOLIN_N <- as.integer(Sys.getenv("VIOLIN_N", "20000"))

RASTERISE <- TRUE

# One colour per sample, in the order SAMPLES gives. Extended rather than
# recycled past two, so a third section does not silently reuse 9w's colour.
SAMPLE_COLS <- c("#4C72B0", "#C44E52", "#55A868", "#8172B2", "#CCB974")

# Same palette as spatial_cell_annotation.R, so a cell type is the same colour
# in every figure of this set.
CT_COLS <- c(
  `Ventricular cardiomyocyte` = "#E8C27A",
  `Atrial cardiomyocyte`      = "#E69F00",
  Fibroblast                  = "#009E73",
  Endothelial                 = "#0072B2",
  Endocardial                 = "#56B4E9",
  `Lymphatic endothelial`     = "#7FDBD4",
  `Pericyte / SMC`            = "#CC79A7",
  Macrophage                  = "#D55E00",
  Lymphocyte                  = "#8B0000",
  Epicardial                  = "#6A3D9A",
  Adipocyte                   = "#F0E442",
  `Schwann / neuronal`        = "#4B4B4B",
  Unassigned                  = "grey85"
)
##### ------------------------------------------------------ #####

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) message(sprintf(...))

BIN_UM <- as.integer(sub("^[^0-9]*0*([0-9]+)um$", "\\1", BIN))
if (is.na(BIN_UM)) stop("Cannot parse a bin size out of BIN = ", BIN)

# ---------------------------------------------------------------- input

# One annotated object at a time, reduced to the columns this script needs and
# then dropped: the 8um objects are ~400 MB each and nothing below wants the
# matrix again.
read_one <- function(sample) {
  d   <- file.path(BASE, "annotation", paste0(sample, "_", BIN))
  rds <- file.path(d, paste0(sample, "_", BIN, "_annotated.rds"))
  if (!file.exists(rds)) stop("No annotated object at ", rds,
                              " - run spatial_cell_annotation.R for ", sample, " first")
  msg("Reading %s ...", rds)
  obj <- readRDS(rds)
  if (!GENE %in% rownames(obj[["Spatial"]])) stop(GENE, " is not in the ", sample,
    " matrix (CreateSeuratObject dropped it at min.cells = 3, or the name differs)")

  cnt <- as.numeric(LayerData(obj, assay = "Spatial", layer = "counts")[GENE, ])
  # The annotated object was NormalizeData()'d, so the data layer is there;
  # recompute LogNormalize by hand if a future object arrives without it,
  # rather than silently plotting raw counts on a log scale.
  nrm <- tryCatch(as.numeric(LayerData(obj, assay = "Spatial", layer = "data")[GENE, ]),
                  error = function(e) NULL)
  if (is.null(nrm)) {
    msg("  no data layer - recomputing log1p(count / nCount * 1e4)")
    nrm <- log1p(cnt / obj$nCount_Spatial * 1e4)
  }

  out <- data.table(
    sample   = sample,
    x        = as.numeric(obj$x),
    y        = as.numeric(obj$y),
    celltype = as.character(obj$celltype),
    nCount   = as.numeric(obj$nCount_Spatial),
    count    = cnt,
    expr     = nrm
  )
  rm(obj); invisible(gc())
  msg("  %d bins, %s+ in %d (%.1f%%), %d %s UMIs",
      nrow(out), GENE, sum(out$count > 0), 100 * mean(out$count > 0),
      as.integer(sum(out$count)), GENE)
  out
}

dt <- rbindlist(lapply(SAMPLES, read_one))
dt[, sample := factor(sample, levels = SAMPLES)]
ct_levels <- c(intersect(names(CT_COLS), unique(dt$celltype)),
               setdiff(unique(dt$celltype), names(CT_COLS)))
dt[, celltype := factor(celltype, levels = ct_levels)]

# ---------------------------------------------------------------- summary

pooled <- function(x, n) if (sum(n) > 0) sum(x) / sum(n) * 1e4 else NA_real_

summ <- dt[, .(
  n_bins       = .N,
  n_pos        = sum(count > 0),
  pct_pos      = round(100 * mean(count > 0), 2),
  gene_umis    = as.integer(sum(count)),
  total_umis   = as.integer(sum(nCount)),
  per_10k      = round(pooled(count, nCount), 3),
  mean_umis    = round(mean(count), 4),
  median_depth = median(nCount)
), by = .(sample, celltype)][order(sample, -n_bins)]

overall <- dt[, .(
  celltype     = "ALL BINS",
  n_bins       = .N,
  n_pos        = sum(count > 0),
  pct_pos      = round(100 * mean(count > 0), 2),
  gene_umis    = as.integer(sum(count)),
  total_umis   = as.integer(sum(nCount)),
  per_10k      = round(pooled(count, nCount), 3),
  mean_umis    = round(mean(count), 4),
  median_depth = median(nCount)
), by = sample]

summ_out <- rbind(overall, summ, fill = TRUE)
fwrite(summ_out, file.path(OUT_DIR, paste0(tolower(GENE), "_summary.csv")))
msg("\n%s per 10k UMIs, all bins:", GENE)
print(overall[, .(sample, n_bins, pct_pos, gene_umis, per_10k)])

# ---------------------------------------------------------------- figures

gpoint <- function(...) {
  g <- geom_point(...)
  if (RASTERISE && requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(g, dpi = 300) else g
}
theme_panel <- function() {
  theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "#52514e"),
          plot.caption = element_text(size = 7, colour = "#52514e", hjust = 0),
          legend.text = element_text(size = 8), legend.title = element_text(size = 9))
}
save_panel <- function(p, name, w = 7, h = 7) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}
if (length(SAMPLES) > length(SAMPLE_COLS))
  stop("Only ", length(SAMPLE_COLS), " sample colours defined - add more to SAMPLE_COLS")
sample_cols <- setNames(SAMPLE_COLS[seq_along(SAMPLES)], SAMPLES)

GENOTYPE <- "B6 x CAST F1, Xist deleted on the B6 allele - all Xist signal is from the CAST (inactive) X."
gl <- tolower(GENE)

# fig1: one map per sample on a shared scale, so the two ages are comparable by
# eye. Zero bins drawn first and in grey - on a continuous scale they would sit
# at the bottom of the viridis ramp and read as low signal rather than none.
# Scale limits come from the POSITIVE bins, not from zero. At 8um a bin holds
# a few dozen UMIs, so one gene UMI already normalises to a large value and
# anchoring the ramp at 0 leaves every positive bin in the top sliver of the
# colour scale - the map then only shows detection. Both ends are quantiles, so
# the ramp spans the middle 98% of the positive bins and the tails squish.
pos <- dt$expr[dt$count > 0]
cap <- if (length(pos)) as.numeric(quantile(pos, CAP_Q)) else 1
lo  <- if (length(pos)) as.numeric(quantile(pos, 1 - CAP_Q)) else 0
if (!is.finite(cap) || cap <= 0) cap <- 1
if (!is.finite(lo) || lo >= cap) lo <- 0
maps <- lapply(SAMPLES, function(s) {
  d <- dt[sample == s][order(count > 0)]
  ggplot(d, aes(x, y)) +
    gpoint(data = d[count == 0], colour = "grey90", size = POINT_SIZE, shape = 16) +
    gpoint(data = d[count > 0], aes(colour = expr), size = POINT_SIZE, shape = 16) +
    scale_colour_viridis_c(option = "magma", direction = -1, limits = c(lo, cap),
                           oob = scales::squish, name = paste0(GENE, "\n(log-norm)")) +
    scale_y_reverse() + coord_fixed() +
    labs(title = s, subtitle = sprintf("%s+ in %.1f%% of bins - %.2f UMIs per 10k",
         GENE, summ_out[sample == s & celltype == "ALL BINS", pct_pos],
         summ_out[sample == s & celltype == "ALL BINS", per_10k])) +
    theme_panel()
})
p1 <- wrap_plots(maps, nrow = 1, guides = "collect") +
  plot_annotation(
    title = paste(GENE, "expression"),
    subtitle = sprintf("Visium HD, %dum bins, log-normalised; grey = no %s UMI; colour spans the %.0fth-%.0fth percentile of positive bins (%.2f-%.2f)",
                       BIN_UM, GENE, 100 * (1 - CAP_Q), 100 * CAP_Q, lo, cap),
    caption = paste(GENOTYPE, "Bins span parts of several cells."),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "#52514e"),
                  plot.caption = element_text(size = 7, colour = "#52514e", hjust = 0)))
save_panel(p1, paste0("fig1_", gl, "_map"), w = 6.5 * length(SAMPLES), h = 7)

# fig2: the two things a detection question actually splits into - how often
# any UMI is seen, and how much there is once depth is accounted for. Cell types
# with too few bins to mean anything are dropped rather than plotted as noise.
MIN_BINS <- 100L
keep_ct <- summ[, .(n = sum(n_bins)), by = celltype][n >= MIN_BINS, as.character(celltype)]
s2 <- summ[as.character(celltype) %in% keep_ct]
s2[, celltype := factor(as.character(celltype), levels = intersect(ct_levels, keep_ct))]
bar <- function(col, ylab, ttl) {
  ggplot(s2, aes(celltype, get(col), fill = sample)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    scale_fill_manual(values = sample_cols, name = NULL) +
    labs(title = ttl, x = NULL, y = ylab) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          plot.title = element_text(face = "bold", size = 11),
          panel.grid.major.x = element_blank())
}
# The cell-type names are long enough that repeating them under both panels
# costs more vertical space than the bars get; only the lower panel keeps them.
p2 <- ((bar("pct_pos", "% of bins", paste0("Bins with at least one ", GENE, " UMI")) +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())) /
       bar("per_10k", paste(GENE, "UMIs per 10k"), "Depth-normalised level (pooled within group)")) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste(GENE, "by cell type"),
    subtitle = sprintf("%dum bins, cell types with >= %d bins; a label is dominance, not purity",
                       BIN_UM, MIN_BINS),
    caption = paste(GENOTYPE,
                    "Detection rate tracks bin depth as well as expression - read it beside the level panel."),
    theme = theme(plot.title = element_text(face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "#52514e"),
                  plot.caption = element_text(size = 7, colour = "#52514e", hjust = 0)))
save_panel(p2, paste0("fig2_", gl, "_by_celltype"), w = 9, h = 8)

# fig3: level among positive bins only. Separates "fewer bins have any" from
# "the bins that have it have less", which fig2's two panels together imply but
# do not show directly.
set.seed(42)
v <- dt[count > 0 & as.character(celltype) %in% keep_ct]
v <- v[, .SD[sample.int(.N, min(.N, VIOLIN_N))], by = .(sample, celltype)]
v[, celltype := factor(as.character(celltype), levels = intersect(ct_levels, keep_ct))]
p3 <- ggplot(v, aes(celltype, expr, fill = sample)) +
  geom_violin(position = position_dodge(width = 0.8), width = 0.75,
              scale = "width", linewidth = 0.2) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.12,
               outlier.shape = NA, linewidth = 0.25, fill = "white") +
  scale_fill_manual(values = sample_cols, name = NULL) +
  labs(title = paste(GENE, "level in positive bins"),
       subtitle = sprintf("log-normalised, %s+ bins only, up to %s bins per group; %dum bins",
                          GENE, format(VIOLIN_N, big.mark = ","), BIN_UM),
       caption = paste(GENOTYPE,
                       "Zero bins excluded on purpose - this panel is level given detection, not abundance."),
       x = NULL, y = paste(GENE, "(log-norm)")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "#52514e"),
        plot.caption = element_text(size = 7, colour = "#52514e", hjust = 0),
        panel.grid.major.x = element_blank())
save_panel(p3, paste0("fig3_", gl, "_level_violin"), w = 9, h = 5.5)

# ---------------------------------------------------------------- provenance

# Column is built as k and renamed: `key` is data.table()'s own argument and
# passing it would set a key instead of making a column (the same trap
# spatial_cell_annotation.R and tile_ratio_map.R flag). The header still has to
# read key/value to match the other provenance.tsv files.
prov <- data.table(
  k = c("script", "run_at", "gene", "bin", "samples", "cap_quantile",
        "violin_n", "min_bins_per_celltype", "seurat_version", "r_version"),
  value = as.character(c("spatial/spatial_xist_panels.R",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            GENE, BIN, paste(SAMPLES, collapse = ","), CAP_Q, VIOLIN_N, MIN_BINS,
            as.character(packageVersion("Seurat")),
            paste0(R.version$major, ".", R.version$minor))))
setnames(prov, "k", "key")
fwrite(prov, file.path(OUT_DIR, "provenance.tsv"), sep = "\t")

msg("\nWrote %s", OUT_DIR)
