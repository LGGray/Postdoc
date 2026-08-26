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

OUT_DIR <- "Allelic_ratio_results/cutoff_sweep"
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

chr_allelic_ratio <- read.table('Allelic_ratio_results/whole_chr_allelic_ratios.txt',
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

# ---------------------------------------------------------------------------
# AR UMAP, built by hand rather than with DimPlot.
#
# The point of doing it this way is draw order. DimPlot plots points in whatever
# order the cells sit in the object, so the rare low-AR (escaping) cells are
# routinely buried under the monoallelic majority and the figure under-reads the
# escape. Sorting descending by AR means the lowest-AR cells are drawn last and
# therefore land on top -- they are the signal, so they win every overlap.
# ---------------------------------------------------------------------------
ar_umap <- function(df, all_cells, title = NULL, point_size = 0.5) {
  # NA bins first (bottom), then high AR -> low AR, so low AR ends up on top
  df <- df[order(is.na(df$allelic_ratio), -df$allelic_ratio), ]

  ggplot() +
    # every retained cell as grey context, so panels with few cells are still
    # readable against the tissue's UMAP shape
    geom_point(data = all_cells, aes(UMAP_1, UMAP_2),
               colour = "grey88", size = point_size * 0.8) +
    geom_point(data = df, aes(UMAP_1, UMAP_2, colour = allelic_bin),
               size = point_size) +
    scale_colour_manual(values = my_colors, drop = FALSE, na.value = "grey60",
                        name = "Allelic ratio") +
    ggtitle(title) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 9))
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
  panels <- lapply(SAMPLE_LEVELS, function(s) {
    d <- flt[flt$sample == s, ]
    ar_umap(d, flt,
            title = sprintf("%s  (n = %s)", s, format(nrow(d), big.mark = ",")))
  })
  p_umap <- patchwork::wrap_plots(panels, nrow = 2) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("chrX allelic ratio, total_reads >= %d  (n = %s cells)",
                      cut_i, format(nrow(flt), big.mark = ",")),
      subtitle = "Lower-AR cells drawn on top; grey = all retained cells",
      theme = theme(plot.title = element_text(size = 11),
                    plot.subtitle = element_text(size = 8))
    ) &
    theme(legend.position = "right")

  ggsave(file.path(OUT_DIR, sprintf("AR_umap_%s.pdf", tag)),
         p_umap, width = 10, height = 7)

  ## ---- Violins, celltype facets ----
  # Only celltypes that still have cells in at least one sample at this cutoff;
  # empty facets otherwise crowd out the ones carrying data.
  vln <- flt %>%
    group_by(celltype) %>%
    filter(dplyr::n() >= 10) %>%
    ungroup() %>%
    mutate(sample = factor(as.character(sample), levels = SAMPLE_LEVELS))

  if (nrow(vln)) {
    vln_n <- vln %>%
      group_by(celltype, sample) %>%
      summarise(n = dplyr::n(), .groups = "drop")

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
           subtitle = "Celltypes with >= 10 cells at this cutoff; numbers are cells per violin") +
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
