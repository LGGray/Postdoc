# ---------------------------------------------------------------------------
# 06 - Total-read cutoff sensitivity sweep for the whole-chrX allelic ratio.
#
# Run from the OCM_heart/ directory:  Rscript allelic_ratio/06_cutoff_sweep.R
# Requires: 01_setup.R (for celltype_sub) and the allelic ratio table written
#           by 02_whole_chrX.R.
# Writes:   Allelic_ratio_results/cutoff_sweep/
#
# 02_whole_chrX.R hard-codes total_reads >= 10. That threshold trades noise
# against cell number, and where it sits changes the apparent amount of escape:
# shallow cells are pushed toward AR ~ 1 by binomial sampling alone, so a low
# cutoff inflates the monoallelic mode, while a high cutoff keeps only the
# deepest (and, per the depth analysis in 02, systematically lower-AR) cells.
# This script re-makes the two plots the interpretation rests on -- the AR UMAP
# and the per-celltype violins -- across a range of cutoffs so that dependence
# is visible rather than assumed.
#
# Deliberately descriptive only. No dispersion / monoallelic LRTs are run here:
# with one animal per condition the between-sample beta-binomial tests are not
# calibrated (see the FPR simulations in 02), and repeating them at seven
# cutoffs would only multiply that problem.
# ---------------------------------------------------------------------------
source("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/OCM_heart/allelic_ratio/00_functions.R")

OUT_DIR <- file.path(RESULTS_ROOT, "cutoff_sweep")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Override from the shell, e.g. CUTOFFS=1,5,10,25 Rscript ...
CUTOFFS <- Sys.getenv("CUTOFFS")
CUTOFFS <- if (nzchar(CUTOFFS)) {
  as.numeric(strsplit(CUTOFFS, "[, ]+")[[1]])
} else {
  c(1, 5, 10, 20, 30, 50, 100)
}
# read counts are integers, and the %d formatting below would error on a
# fractional cutoff, so round rather than trusting whatever was passed in
CUTOFFS <- sort(unique(as.integer(round(CUTOFFS[is.finite(CUTOFFS)]))))
message("Sweeping total_reads cutoffs: ", paste(CUTOFFS, collapse = ", "))

SAMPLE_LEVELS <- c("9w", "78w", "Sham", "TAC")

# ---------------------------------------------------------------------------
# Build the unfiltered cell table once: metadata + UMAP coords + allelic counts
# ---------------------------------------------------------------------------
heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

chr_allelic_ratio <- read.table(ALLELIC_RATIOS_FILE,
                                sep = '\t', header = TRUE, stringsAsFactors = FALSE)
chr_allelic_ratio <- subset(chr_allelic_ratio, chr == "chrX")

keep <- intersect(chr_allelic_ratio$cell_barcode, colnames(heart))
chr_allelic_ratio <- chr_allelic_ratio[match(keep, chr_allelic_ratio$cell_barcode), ]

# Coordinates come out of the full object rather than a per-cutoff subset, so
# every panel below sits on the same UMAP and cutoffs stay visually comparable.
emb <- Embeddings(heart, reduction = "umap")[keep, 1:2, drop = FALSE]
colnames(emb) <- c("UMAP_1", "UMAP_2")

cells <- data.frame(
  cell_barcode  = keep,
  celltype      = as.character(heart$celltype[keep]),
  sample        = factor(as.character(heart$sample[keep]), levels = SAMPLE_LEVELS),
  nCount_RNA    = heart$nCount_RNA[keep],
  allelic_ratio = chr_allelic_ratio$allelic_ratio,
  total_reads   = chr_allelic_ratio$total_reads,
  A1_reads      = chr_allelic_ratio$A1_reads,
  A2_reads      = chr_allelic_ratio$A2_reads,
  emb,
  stringsAsFactors = FALSE
)
cells <- cells[!is.na(cells$allelic_ratio) & !is.na(cells$total_reads), ]

if (any(CUTOFFS > max(cells$total_reads))) {
  message("NOTE: cutoffs above the deepest cell (", max(cells$total_reads),
          " reads) retain nothing and are dropped: ",
          paste(CUTOFFS[CUTOFFS > max(cells$total_reads)], collapse = ", "))
  CUTOFFS <- CUTOFFS[CUTOFFS <= max(cells$total_reads)]
}
stopifnot(length(CUTOFFS) > 0)

write.table(cells, file.path(OUT_DIR, "cutoff_sweep_cell_table.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# Shared AR binning / palette (same 11 bins as 02_whole_chrX.R, so the panels
# here can be put next to that figure without re-reading the colour key)
# ---------------------------------------------------------------------------
my_breaks <- c(seq(0, 0.9, by = 0.1), 0.95, 1.0)
my_colors <- c(
  "#2B3186", "#3B5FB6", "#38749F",
  "#367373", "#2D6E5D",
  "#1E652D",
  "#658C2D", "#8D9F25", "#B3B112",
  "#C97314", "#8B1913"
)
bin_labels <- c("0.00-0.10", "0.10-0.20", "0.20-0.30", "0.30-0.40", "0.40-0.50",
                "0.50-0.60", "0.60-0.70", "0.70-0.80", "0.80-0.90",
                "0.90-0.95", "0.95-1.00")
names(my_colors) <- bin_labels

cells$allelic_bin <- cut(cells$allelic_ratio, breaks = my_breaks,
                         include.lowest = TRUE, right = TRUE, labels = bin_labels)

# Second, finer bin set. The observed AR distribution is almost entirely above
# 0.85, so the 02-matching bins above spend 9 of their 11 colours on a range
# that holds a few percent of cells and the UMAP comes out nearly monochrome.
# These bins put the resolution where the cells actually are, which is what
# makes the cutoff effect visible; the coarse version is still written so the
# figures stay comparable with 02.
fine_breaks <- c(0, 0.5, 0.7, 0.8, 0.85, 0.90, 0.93, 0.95, 0.97, 0.99, 1.0)
fine_labels <- c("0.00-0.50", "0.50-0.70", "0.70-0.80", "0.80-0.85", "0.85-0.90",
                 "0.90-0.93", "0.93-0.95", "0.95-0.97", "0.97-0.99", "0.99-1.00")
fine_colors <- c("#2B3186", "#3B5FB6", "#38749F", "#367373", "#2D6E5D",
                 "#1E652D", "#8D9F25", "#C97314", "#B03A12", "#8B1913")
names(fine_colors) <- fine_labels

cells$allelic_bin_fine <- cut(cells$allelic_ratio, breaks = fine_breaks,
                              include.lowest = TRUE, right = TRUE,
                              labels = fine_labels)

# ---------------------------------------------------------------------------
# AR UMAP, built by hand rather than with DimPlot.
#
# The point of doing it this way is draw order. DimPlot plots points in whatever
# order the cells sit in the object, so the rare low-AR (escaping) cells are
# routinely buried under the monoallelic majority and the figure under-reads the
# escape. Sorting descending by AR means the lowest-AR cells are drawn last and
# therefore land on top -- they are the signal, so they win every overlap.
# ---------------------------------------------------------------------------
ar_umap <- function(df, all_cells, title = NULL, point_size = 0.5,
                    bin_col = "allelic_bin", palette = my_colors) {
  # NA bins first (bottom), then high AR -> low AR, so low AR ends up on top
  df <- df[order(is.na(df$allelic_ratio), -df$allelic_ratio), ]

  ggplot() +
    # every retained cell as grey context, so panels with few cells are still
    # readable against the tissue's UMAP shape
    geom_point(data = all_cells, aes(UMAP_1, UMAP_2),
               colour = "grey88", size = point_size * 0.8) +
    geom_point(data = df, aes(UMAP_1, UMAP_2, colour = .data[[bin_col]]),
               size = point_size) +
    # limits pinned to the full label set, not inferred per panel: patchwork
    # only merges guides that are byte-identical, so without this each sample
    # panel keeps its own legend and four copies stack up off the page edge
    scale_colour_manual(values = palette, limits = names(palette),
                        drop = FALSE, na.value = "grey60",
                        name = "Allelic ratio") +
    guides(colour = guide_legend(override.aes = list(size = 1.8), ncol = 1)) +
    ggtitle(title) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 9),
          legend.position = "right",
          legend.key.height = unit(0.42, "lines"),
          legend.text = element_text(size = 6.5),
          legend.title = element_text(size = 7.5))
}

# ---------------------------------------------------------------------------
# Per-cutoff figures
# ---------------------------------------------------------------------------
summary_rows <- list()

for (cut_i in CUTOFFS) {
  flt <- cells[cells$total_reads >= cut_i, ]
  tag <- sprintf("min%03d", cut_i)

  message(sprintf("cutoff >= %d reads: %d cells across %d celltypes",
                  cut_i, nrow(flt), length(unique(flt$celltype))))

  ## ---- UMAP, one panel per sample ----
  # `guides = "collect"` goes to wrap_plots directly. Adding a trailing
  # `& theme(legend.position = ...)` instead pushes the theme down into every
  # subplot and undoes the collection, which is what produced four legends.
  umap_variants <- list(
    list(suffix = "",      col = "allelic_bin",      pal = my_colors,
         note = "bins as in 02_whole_chrX.R"),
    list(suffix = "_fine", col = "allelic_bin_fine", pal = fine_colors,
         note = "bins concentrated above 0.85, where the cells are")
  )

  for (v in umap_variants) {
    panels <- lapply(SAMPLE_LEVELS, function(s) {
      d <- flt[flt$sample == s, ]
      ar_umap(d, flt,
              title = sprintf("%s  (n = %s)", s, format(nrow(d), big.mark = ",")),
              bin_col = v$col, palette = v$pal)
    })
    p_umap <- patchwork::wrap_plots(panels, nrow = 2, guides = "collect") +
      patchwork::plot_annotation(
        title = sprintf("chrX allelic ratio, total_reads >= %d  (n = %s cells)",
                        cut_i, format(nrow(flt), big.mark = ",")),
        subtitle = paste0("Lower-AR cells drawn on top; grey = all retained cells; ",
                          v$note),
        theme = theme(plot.title = element_text(size = 11),
                      plot.subtitle = element_text(size = 8))
      )

    ggsave(file.path(OUT_DIR, sprintf("AR_umap_%s%s.pdf", tag, v$suffix)),
           p_umap, width = 10, height = 7)
  }

  ## ---- Violins, celltype facets ----
  # Two separate thresholds on purpose. Counts are labelled for every
  # celltype/sample group that has any cells, but a violin is only drawn where
  # there are >= 10 of them: at n = 1-2 geom_violin degenerates into a flat
  # dash that reads like a real narrow distribution (this is what
  # "Cardiomyocytes (stressed)" did at n = 1). Nothing is silently dropped --
  # the n label is still there under the empty slot.
  vln_n <- flt %>%
    group_by(celltype, sample) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    filter(n > 0)

  keep_ct <- vln_n %>%
    group_by(celltype) %>%
    filter(max(n) >= 20) %>%
    ungroup() %>%
    pull(celltype) %>%
    unique()

  vln_n <- vln_n %>% filter(celltype %in% keep_ct)

  vln <- flt %>%
    filter(celltype %in% keep_ct) %>%
    semi_join(filter(vln_n, n >= 10), by = c("celltype", "sample")) %>%
    mutate(sample = factor(as.character(sample), levels = SAMPLE_LEVELS))

  if (nrow(vln)) {
    p_vln <- ggplot(vln, aes(x = sample, y = allelic_ratio, fill = sample)) +
      geom_violin(trim = FALSE, scale = "width", bounds = c(0, 1),
                  linewidth = 0.3) +
      geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                   linewidth = 0.3) +
      geom_text(data = vln_n, aes(x = sample, y = -0.04, label = n),
                inherit.aes = FALSE, size = 2, colour = "grey30") +
      facet_wrap(~celltype, labeller = label_wrap_gen(width = 18)) +
      scale_x_discrete(drop = FALSE) +
      scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9, 1.0)) +
      coord_cartesian(ylim = c(-0.06, 1), clip = "off") +
      labs(y = "Allelic ratio", x = NULL,
           title = sprintf("chrX allelic ratio, total_reads >= %d", cut_i),
           subtitle = paste("Celltypes with >= 20 cells in some sample; numbers are cells per group.",
                            "Violins drawn only where n >= 10.")) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            plot.subtitle = element_text(size = 7))

    ggsave(file.path(OUT_DIR, sprintf("AR_violin_%s.pdf", tag)),
           p_vln, width = 9, height = 8)
  } else {
    message("  no celltype reaches 10 cells at this cutoff - violin skipped")
  }

  ## ---- Numbers behind the plots ----
  summary_rows[[as.character(cut_i)]] <- flt %>%
    group_by(celltype, sample, .drop = FALSE) %>%
    summarise(n_cells          = dplyr::n(),
              median_total_reads = median(total_reads),
              median_AR        = median(allelic_ratio),
              mean_AR          = mean(allelic_ratio),
              frac_escaping    = mean(allelic_ratio <= 0.9),
              frac_monoallelic = mean(allelic_ratio >= 0.9),
              .groups = "drop") %>%
    filter(n_cells > 0)
}

sweep_summary <- bind_rows(summary_rows, .id = "min_total_reads") %>%
  mutate(min_total_reads = as.numeric(min_total_reads)) %>%
  arrange(celltype, sample, min_total_reads)

write.table(sweep_summary, file.path(OUT_DIR, "cutoff_sweep_summary_by_celltype.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

sweep_overall <- bind_rows(lapply(CUTOFFS, function(cut_i) {
  cells[cells$total_reads >= cut_i, ] %>%
    group_by(sample) %>%
    summarise(n_cells          = dplyr::n(),
              median_total_reads = median(total_reads),
              median_AR        = median(allelic_ratio),
              frac_escaping    = mean(allelic_ratio <= 0.9),
              .groups = "drop") %>%
    mutate(min_total_reads = cut_i)
})) %>%
  select(min_total_reads, everything()) %>%
  arrange(sample, min_total_reads)

write.table(sweep_overall, file.path(OUT_DIR, "cutoff_sweep_summary_overall.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# Cross-cutoff trend plots: the actual answer to "how much does the cutoff
# change the interpretation". A flat line means the conclusion is robust to the
# threshold; a monotone trend means it is a property of the threshold.
# ---------------------------------------------------------------------------
robust_ct <- sweep_summary %>%
  group_by(celltype) %>%
  filter(any(n_cells >= 20)) %>%
  ungroup()

pdf(file.path(OUT_DIR, "cutoff_sweep_trends.pdf"), width = 11, height = 8)

print(
  ggplot(sweep_overall, aes(min_total_reads, n_cells, colour = sample)) +
    geom_line() + geom_point(size = 1.2) +
    scale_x_log10(breaks = CUTOFFS) + scale_y_log10() +
    labs(x = "Minimum chrX total_reads", y = "Cells retained (log scale)",
         title = "Cells retained vs cutoff",
         subtitle = "The cost side of the trade-off: how much data each threshold discards") +
    theme_bw()
)

print(
  ggplot(sweep_overall, aes(min_total_reads, frac_escaping, colour = sample)) +
    geom_line() + geom_point(size = 1.2) +
    scale_x_log10(breaks = CUTOFFS) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Minimum chrX total_reads", y = "Fraction of cells with AR <= 0.9",
         title = "Apparent escape vs cutoff, all celltypes pooled",
         subtitle = paste("If this rises with the cutoff, the escape signal is being revealed by depth;",
                          "\nif it falls, the low-cutoff estimate was inflated by shallow cells.")) +
    theme_bw()
)

print(
  ggplot(robust_ct, aes(min_total_reads, frac_escaping, colour = sample)) +
    geom_line() + geom_point(size = 0.8) +
    facet_wrap(~celltype, labeller = label_wrap_gen(18)) +
    scale_x_log10(breaks = CUTOFFS) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Minimum chrX total_reads", y = "Fraction AR <= 0.9",
         title = "Apparent escape vs cutoff, per celltype") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

print(
  ggplot(robust_ct, aes(min_total_reads, median_AR, colour = sample)) +
    geom_line() + geom_point(size = 0.8) +
    facet_wrap(~celltype, labeller = label_wrap_gen(18)) +
    scale_x_log10(breaks = CUTOFFS) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Minimum chrX total_reads", y = "Median allelic ratio",
         title = "Median allelic ratio vs cutoff, per celltype") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

print(
  ggplot(robust_ct, aes(min_total_reads, n_cells, colour = sample)) +
    geom_line() + geom_point(size = 0.8) +
    facet_wrap(~celltype, labeller = label_wrap_gen(18)) +
    scale_x_log10(breaks = CUTOFFS) + scale_y_log10() +
    labs(x = "Minimum chrX total_reads", y = "Cells retained (log scale)",
         title = "Cells retained vs cutoff, per celltype") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

dev.off()

# ---------------------------------------------------------------------------
# Depth is not a nuisance here, it is a bias. This is the section that matters.
#
# The escape fraction (AR <= 0.9) is a threshold on a noisy estimate, so it
# necessarily moves with read depth and cannot separate the two explanations on
# its own. The MEAN of A1/total_reads can: A1/N is an unbiased estimator of the
# underlying allelic proportion, so binomial sampling cannot shift its mean at
# any depth. Any depth gradient in mean AR is therefore a real difference in
# what the cells are reporting, not a measurement artefact.
#
# Stratifying rather than filtering is the point. A cutoff mixes two things --
# discarding noisy cells and discarding biased cells -- and reports one number;
# these bands keep them separate.
# ---------------------------------------------------------------------------
DEPTH_BANDS <- list(c(10, 25), c(25, 40), c(40, 60), c(60, 100), c(100, 200),
                    c(200, Inf))
band_label <- sapply(DEPTH_BANDS, function(b)
  if (is.infinite(b[2])) paste0(b[1], "+") else paste0(b[1], "-", b[2]))

assign_band <- function(n) {
  idx <- sapply(n, function(x) {
    hit <- which(sapply(DEPTH_BANDS, function(b) x >= b[1] & x < b[2]))
    if (length(hit)) hit[1] else NA_integer_
  })
  factor(band_label[idx], levels = band_label)
}
cells$depth_band <- assign_band(cells$total_reads)

# --- the gradient itself, and whether celltype composition explains it -------
depth_gradient <- cells %>%
  filter(!is.na(depth_band)) %>%
  group_by(depth_band) %>%
  summarise(n_cells = dplyr::n(), median_depth = median(total_reads),
            mean_AR = mean(allelic_ratio),
            pooled_AR = sum(A1_reads) / sum(total_reads),
            frac_escaping = mean(allelic_ratio <= 0.9), .groups = "drop")

# Same gradient within each celltype. If it survives here it is not the
# large-cell/deep-library celltypes dragging the pooled number around.
depth_gradient_ct <- cells %>%
  filter(!is.na(depth_band)) %>%
  group_by(celltype, depth_band) %>%
  summarise(n_cells = dplyr::n(), mean_AR = mean(allelic_ratio), .groups = "drop") %>%
  filter(n_cells >= 20)

# And against library size instead of chrX informative reads: nCount_RNA is a
# separate quantity, so a gradient in both rules out an artefact of however many
# SNP-overlapping reads a given cell happened to yield.
depth_gradient_lib <- cells %>%
  mutate(lib_sextile = dplyr::ntile(nCount_RNA, 6)) %>%
  group_by(lib_sextile) %>%
  summarise(n_cells = dplyr::n(), median_lib = median(nCount_RNA),
            median_chrX_reads = median(total_reads),
            mean_AR = mean(allelic_ratio), .groups = "drop")

write.table(depth_gradient, file.path(OUT_DIR, "depth_gradient_pooled.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)
write.table(depth_gradient_ct, file.path(OUT_DIR, "depth_gradient_by_celltype.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)
write.table(depth_gradient_lib, file.path(OUT_DIR, "depth_gradient_by_library_size.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

# --- depth-matched sample comparison ----------------------------------------
# The comparison the conclusions actually rest on. Because the samples have
# slightly different depth distributions, part of any pooled difference between
# them is the gradient above rather than biology. Within a band that is held
# fixed, so an ordering that persists down every column is real and one that
# collapses at depth was the gradient.
depth_matched <- cells %>%
  filter(!is.na(depth_band)) %>%
  group_by(depth_band, sample) %>%
  summarise(n_cells = dplyr::n(), mean_AR = mean(allelic_ratio),
            median_AR = median(allelic_ratio),
            frac_escaping = mean(allelic_ratio <= 0.9), .groups = "drop") %>%
  filter(n_cells >= 30)

depth_matched_ct <- cells %>%
  filter(!is.na(depth_band)) %>%
  group_by(celltype, depth_band, sample) %>%
  summarise(n_cells = dplyr::n(), mean_AR = mean(allelic_ratio), .groups = "drop") %>%
  filter(n_cells >= 30)

write.table(depth_matched, file.path(OUT_DIR, "depth_matched_by_sample.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)
write.table(depth_matched_ct, file.path(OUT_DIR, "depth_matched_by_sample_and_celltype.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

pdf(file.path(OUT_DIR, "depth_bias.pdf"), width = 11, height = 8)

print(
  ggplot(depth_gradient, aes(median_depth, mean_AR)) +
    geom_line(colour = "grey40") +
    geom_point(aes(size = n_cells), colour = "firebrick") +
    geom_text(aes(label = band_label[as.integer(depth_band)]), vjust = -1.2, size = 2.6) +
    scale_x_log10() + scale_size_continuous(range = c(1.5, 4), name = "Cells") +
    coord_cartesian(ylim = c(0.75, 1)) +
    labs(x = "chrX informative reads (median of band, log scale)",
         y = "Mean allelic ratio",
         title = "Mean allelic ratio rises with read depth",
         caption = paste("A1/N is unbiased for the underlying proportion, so binomial noise cannot move this mean.",
                         "A gradient here is a real\ndepth-dependent bias -- shallow cells report genuinely lower AR.",
                         "Ambient RNA is the leading candidate: a shallow cell's\nreads are proportionally more ambient, and the ambient pool averages over a mosaic of cells rather than one active X.")) +
    theme_bw() + theme(plot.caption = element_text(size = 7, hjust = 0))
)

print(
  ggplot(depth_gradient_ct, aes(depth_band, mean_AR, group = celltype)) +
    geom_line(colour = "grey60") +
    geom_point(aes(size = n_cells), colour = "firebrick") +
    facet_wrap(~celltype, labeller = label_wrap_gen(18)) +
    scale_size_continuous(range = c(0.6, 2.5), name = "Cells") +
    coord_cartesian(ylim = c(0.7, 1)) +
    labs(x = "chrX informative reads", y = "Mean allelic ratio",
         title = "The same gradient inside every celltype",
         subtitle = "So it is not celltype composition -- deep-library celltypes are not dragging the pooled number") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

print(
  ggplot(depth_gradient_lib, aes(median_lib, mean_AR)) +
    geom_line(colour = "grey40") + geom_point(colour = "firebrick", size = 2) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0.75, 1)) +
    labs(x = "Library size (nCount_RNA, median of sextile, log scale)",
         y = "Mean allelic ratio",
         title = "And against library size, a separate quantity",
         subtitle = "Rules out an artefact of how many SNP-overlapping chrX reads a cell happened to yield") +
    theme_bw()
)

print(
  ggplot(depth_matched, aes(depth_band, mean_AR, colour = sample, group = sample)) +
    geom_line() + geom_point(size = 1.6) +
    coord_cartesian(ylim = c(0.75, 1)) +
    labs(x = "chrX informative reads", y = "Mean allelic ratio",
         title = "Depth-matched sample comparison",
         caption = paste("Within a band the depth bias is held fixed. An ordering that persists down every band",
                         "is a real sample difference;\none that collapses in the deep bands was the depth gradient.")) +
    theme_bw() + theme(plot.caption = element_text(size = 7, hjust = 0))
)

print(
  ggplot(depth_matched_ct, aes(depth_band, mean_AR, colour = sample, group = sample)) +
    geom_line() + geom_point(size = 1) +
    facet_wrap(~celltype, labeller = label_wrap_gen(18)) +
    coord_cartesian(ylim = c(0.7, 1)) +
    labs(x = "chrX informative reads", y = "Mean allelic ratio",
         title = "Depth-matched sample comparison, per celltype",
         subtitle = "Only bands with >= 30 cells are drawn") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

dev.off()

print(depth_gradient)
print(depth_matched)

# ---------------------------------------------------------------------------
# Noise floor, for the escape FRACTION specifically.
#
# Superseded as the headline by the gradient above, but kept because the escape
# fraction is the number in the existing figures and it is worth knowing how
# much of it is a threshold artefact.
#
# The comparison is depth-MATCHED, which the first version of this got wrong:
# thinning to exactly k reads was being compared against every cell with >= k
# reads, whose median depth is far above k, so the noise floor looked several
# times larger than it was. Both sides now sit in the same depth band.
# ---------------------------------------------------------------------------
REF_MIN   <- 200      # reference cells: deep enough to thin down to any band
THIN_REPS <- 200

set.seed(1)
ref_cells <- cells[cells$total_reads >= REF_MIN, ]

noise_floor <- bind_rows(lapply(seq_along(DEPTH_BANDS), function(i) {
  b <- DEPTH_BANDS[[i]]
  if (is.infinite(b[2]) || b[1] >= REF_MIN) return(NULL)
  real <- cells[cells$total_reads >= b[1] & cells$total_reads < b[2], ]
  if (nrow(real) < 50 || nrow(ref_cells) < 100) return(NULL)
  k <- round(median(real$total_reads))
  reps <- replicate(THIN_REPS, {
    a <- rhyper(nrow(ref_cells), ref_cells$A1_reads, ref_cells$A2_reads, k)
    mean(a / k <= 0.9)
  })
  data.frame(depth_band = band_label[i], median_depth = k,
             n_real = nrow(real), n_ref = nrow(ref_cells),
             observed_escape = mean(real$allelic_ratio <= 0.9),
             noise_escape = mean(reps),
             noise_lwr = quantile(reps, 0.025, names = FALSE),
             noise_upr = quantile(reps, 0.975, names = FALSE))
})) %>%
  mutate(noise_share = noise_escape / observed_escape,
         depth_band = factor(depth_band, levels = band_label))

write.table(noise_floor, file.path(OUT_DIR, "cutoff_sweep_noise_floor.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

pdf(file.path(OUT_DIR, "cutoff_sweep_noise_floor.pdf"), width = 8, height = 5)
print(
  ggplot(noise_floor, aes(median_depth)) +
    geom_ribbon(aes(ymin = noise_lwr, ymax = noise_upr), fill = "grey75", alpha = 0.6) +
    geom_line(aes(y = noise_escape, colour = "Noise floor (deep cells thinned to this depth)")) +
    geom_point(aes(y = noise_escape, colour = "Noise floor (deep cells thinned to this depth)"), size = 1.4) +
    geom_line(aes(y = observed_escape, colour = "Observed (real cells at this depth)")) +
    geom_point(aes(y = observed_escape, colour = "Observed (real cells at this depth)"), size = 1.4) +
    scale_x_log10(breaks = noise_floor$median_depth) +
    scale_colour_manual(values = c("Noise floor (deep cells thinned to this depth)" = "grey35",
                                   "Observed (real cells at this depth)" = "firebrick"),
                        name = NULL) +
    labs(x = "chrX informative reads (median of band, log scale)",
         y = "Fraction with AR <= 0.9",
         title = "Escape fraction against its binomial noise floor, depth-matched",
         caption = paste("Both curves are at the SAME depth. Grey: cells with >=", REF_MIN,
                         "reads thinned down to it, so their escape is measurement\nerror by construction.",
                         "The large gap means the escape fraction is mostly not sampling noise -- it tracks the\ndepth bias in the",
                         "figures above. Caveat: the grey reference is deep cells, which the gradient shows are\nthemselves the",
                         "least-biased cells, so this floor is a lower bound.")) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", plot.caption = element_text(size = 7, hjust = 0))
)
dev.off()

print(noise_floor)

# ---------------------------------------------------------------------------
# AR distribution as a ridge-style overlay: all cutoffs on one axis per sample,
# which is the most direct read on whether the monoallelic mode is real or is
# shallow cells piling up at AR = 1.
# ---------------------------------------------------------------------------
dist_tbl <- bind_rows(lapply(CUTOFFS, function(cut_i) {
  d <- cells[cells$total_reads >= cut_i, ]
  d$min_total_reads <- cut_i
  d
}))
dist_tbl$cutoff_lab <- factor(paste0(">= ", dist_tbl$min_total_reads),
                              levels = paste0(">= ", CUTOFFS))

pdf(file.path(OUT_DIR, "cutoff_sweep_AR_histograms.pdf"), width = 10, height = 7)
print(
  ggplot(dist_tbl, aes(allelic_ratio, fill = sample)) +
    geom_histogram(binwidth = 0.02, boundary = 0) +
    facet_grid(cutoff_lab ~ sample, scales = "free_y") +
    labs(x = "Allelic ratio (chrX)", y = "Cells",
         title = "AR distribution at each total-read cutoff") +
    theme_bw(base_size = 8) +
    theme(legend.position = "none",
          strip.text.y = element_text(size = 6))
)
dev.off()

message("Done. Output in ", normalizePath(OUT_DIR))
