# ---------------------------------------------------------------------------
# 07 - Total-read cutoff sensitivity sweep for the CORE ESCAPE BLOCK.
#
# Run from the OCM_heart/ directory:  Rscript allelic_ratio/07_core_escape_cutoff_sweep.R
# Requires: 01_setup.R, and 04_core_escape.R to have written the raw ratio
#           table (CEB_RATIOS_FILE). Does not need 04's filtered outputs.
# Writes:   Allelic_ratio_results/core_escape_cutoff_sweep/
#
# The counterpart to 06_cutoff_sweep.R, but the core escape block is a much
# harder measurement problem and the sweep has to be read differently.
#
# Whole chrX has a median of ~50 informative reads per cell. This block --
# roughly twenty escape genes -- has a median of THREE, and a maximum of 26.
# At three reads the allelic ratio can only take four values, so a large part
# of what the per-cell AR distribution looks like is set by arithmetic rather
# than biology. That makes the cutoff choice far more consequential here than
# for whole chrX, and it pushes in the opposite direction: raising the cutoff
# does not mainly reduce an ambient-driven bias, it buys the resolution needed
# for the ratio to mean anything at all.
#
# Descriptive only, as in 06. No between-sample beta-binomial tests: n=1 animal
# per condition, and here the per-cell estimates are far noisier besides.
# ---------------------------------------------------------------------------
source("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/OCM_heart/allelic_ratio/00_functions.R")

OUT_DIR <- "Allelic_ratio_results/core_escape_cutoff_sweep"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Cutoffs are an order of magnitude lower than 06's on purpose: the deepest
# cell in this block has 26 reads, so 06's 30-100 range would retain nothing.
CUTOFFS <- Sys.getenv("CUTOFFS")
CUTOFFS <- if (nzchar(CUTOFFS)) {
  as.numeric(strsplit(CUTOFFS, "[, ]+")[[1]])
} else {
  c(1, 2, 3, 4, 5, 8, 10)
}
CUTOFFS <- sort(unique(as.integer(round(CUTOFFS[is.finite(CUTOFFS)]))))

SAMPLE_LEVELS <- c("9w", "78w", "Sham", "TAC")
MIN_CELLS_PER_CUTOFF <- 50   # below this a panel is not worth drawing

# ---------------------------------------------------------------------------
# Cell table: raw (unfiltered) ratios + metadata + a shared UMAP
# ---------------------------------------------------------------------------
heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

ceb <- read.delim(CEB_RATIOS_FILE, header = TRUE, stringsAsFactors = FALSE)
ceb <- subset(ceb, chr == "chrX")
stopifnot(!any(duplicated(ceb$cell_barcode)))   # 04 aggregates to one row per cell

keep <- intersect(ceb$cell_barcode, colnames(heart))
ceb <- ceb[match(keep, ceb$cell_barcode), ]

emb <- Embeddings(heart, reduction = "umap")[keep, 1:2, drop = FALSE]
colnames(emb) <- c("UMAP_1", "UMAP_2")

cells <- data.frame(
  cell_barcode  = keep,
  celltype      = as.character(heart$celltype[keep]),
  sample        = factor(as.character(heart$sample[keep]), levels = SAMPLE_LEVELS),
  nCount_RNA    = heart$nCount_RNA[keep],
  allelic_ratio = ceb$allelic_ratio,
  total_reads   = ceb$total_reads,
  A1_reads      = ceb$A1_reads,
  A2_reads      = ceb$A2_reads,
  emb,
  stringsAsFactors = FALSE
)
cells <- cells[!is.na(cells$allelic_ratio) & !is.na(cells$total_reads), ]

CUTOFFS <- CUTOFFS[sapply(CUTOFFS, function(k) sum(cells$total_reads >= k)) >= MIN_CELLS_PER_CUTOFF]
stopifnot(length(CUTOFFS) > 0)
message("Sweeping core escape cutoffs: ", paste(CUTOFFS, collapse = ", "))
message("Depth: median ", median(cells$total_reads), " reads, max ", max(cells$total_reads))

write.table(cells, file.path(OUT_DIR, "ceb_cutoff_sweep_cell_table.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# Palette. Unlike whole chrX, this block's AR really does span 0-1 (mean ~0.63,
# these are escape genes), so the even 11-bin scheme from 02 is the right one
# here and no compressed high-end variant is needed.
# ---------------------------------------------------------------------------
my_breaks <- c(seq(0, 0.9, by = 0.1), 0.95, 1.0)
my_colors <- c("#2B3186", "#3B5FB6", "#38749F", "#367373", "#2D6E5D", "#1E652D",
               "#658C2D", "#8D9F25", "#B3B112", "#C97314", "#8B1913")
bin_labels <- c("0.00-0.10", "0.10-0.20", "0.20-0.30", "0.30-0.40", "0.40-0.50",
                "0.50-0.60", "0.60-0.70", "0.70-0.80", "0.80-0.90",
                "0.90-0.95", "0.95-1.00")
names(my_colors) <- bin_labels
cells$allelic_bin <- cut(cells$allelic_ratio, breaks = my_breaks,
                         include.lowest = TRUE, right = TRUE, labels = bin_labels)

# Same hand-built UMAP as 06: DimPlot draws in object order, which buries the
# low-AR cells under the rest. Sorting descending by AR puts them on top.
ar_umap <- function(df, all_cells, title = NULL, point_size = 0.5) {
  df <- df[order(is.na(df$allelic_ratio), -df$allelic_ratio), ]
  ggplot() +
    geom_point(data = all_cells, aes(UMAP_1, UMAP_2), colour = "grey88",
               size = point_size * 0.8) +
    geom_point(data = df, aes(UMAP_1, UMAP_2, colour = allelic_bin),
               size = point_size) +
    scale_colour_manual(values = my_colors, limits = names(my_colors),
                        drop = FALSE, na.value = "grey60", name = "Allelic ratio") +
    guides(colour = guide_legend(override.aes = list(size = 1.8), ncol = 1)) +
    ggtitle(title) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(),
          plot.title = element_text(size = 9), legend.position = "right",
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
  tag <- sprintf("min%02d", cut_i)
  message(sprintf("cutoff >= %d reads: %d cells", cut_i, nrow(flt)))

  panels <- lapply(SAMPLE_LEVELS, function(s) {
    d <- flt[flt$sample == s, ]
    ar_umap(d, flt, title = sprintf("%s  (n = %s)", s,
                                    format(nrow(d), big.mark = ",")))
  })
  p_umap <- patchwork::wrap_plots(panels, nrow = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("Core escape block allelic ratio, total_reads >= %d  (n = %s cells)",
                      cut_i, format(nrow(flt), big.mark = ",")),
      subtitle = "Lower-AR cells drawn on top; grey = all retained cells",
      theme = theme(plot.title = element_text(size = 11),
                    plot.subtitle = element_text(size = 8)))
  ggsave(file.path(OUT_DIR, sprintf("ceb_AR_umap_%s.pdf", tag)),
         p_umap, width = 10, height = 7)

  vln_n <- flt %>%
    group_by(celltype, sample) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    filter(n > 0)
  keep_ct <- vln_n %>% group_by(celltype) %>% filter(max(n) >= 20) %>%
    ungroup() %>% pull(celltype) %>% unique()
  vln_n <- vln_n %>% filter(celltype %in% keep_ct)
  vln <- flt %>%
    filter(celltype %in% keep_ct) %>%
    semi_join(filter(vln_n, n >= 10), by = c("celltype", "sample")) %>%
    mutate(sample = factor(as.character(sample), levels = SAMPLE_LEVELS))

  if (nrow(vln)) {
    p_vln <- ggplot(vln, aes(x = sample, y = allelic_ratio, fill = sample)) +
      geom_violin(trim = FALSE, scale = "width", bounds = c(0, 1), linewidth = 0.3) +
      geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
      geom_text(data = vln_n, aes(x = sample, y = -0.04, label = n),
                inherit.aes = FALSE, size = 2, colour = "grey30") +
      facet_wrap(~celltype, labeller = label_wrap_gen(width = 18)) +
      scale_x_discrete(drop = FALSE) +
      scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
      coord_cartesian(ylim = c(-0.06, 1), clip = "off") +
      labs(y = "Allelic ratio (core escape block)", x = NULL,
           title = sprintf("Core escape block AR, total_reads >= %d", cut_i),
           subtitle = paste("Celltypes with >= 20 cells in some sample; numbers are cells per group.",
                            "Violins drawn only where n >= 10.")) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none", plot.subtitle = element_text(size = 7))
    ggsave(file.path(OUT_DIR, sprintf("ceb_AR_violin_%s.pdf", tag)),
           p_vln, width = 9, height = 8)
  }

  summary_rows[[as.character(cut_i)]] <- flt %>%
    group_by(celltype, sample, .drop = FALSE) %>%
    summarise(n_cells = dplyr::n(),
              median_total_reads = median(total_reads),
              mean_AR = mean(allelic_ratio), median_AR = median(allelic_ratio),
              frac_escaping = mean(allelic_ratio <= 0.9),
              frac_exactly_1 = mean(allelic_ratio == 1),
              frac_exactly_0 = mean(allelic_ratio == 0),
              .groups = "drop") %>%
    filter(n_cells > 0)
}

sweep_summary <- bind_rows(summary_rows, .id = "min_total_reads") %>%
  mutate(min_total_reads = as.numeric(min_total_reads)) %>%
  arrange(celltype, sample, min_total_reads)
write.table(sweep_summary, file.path(OUT_DIR, "ceb_cutoff_sweep_summary_by_celltype.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

sweep_overall <- bind_rows(lapply(CUTOFFS, function(cut_i) {
  cells[cells$total_reads >= cut_i, ] %>%
    group_by(sample) %>%
    summarise(n_cells = dplyr::n(), median_total_reads = median(total_reads),
              mean_AR = mean(allelic_ratio), median_AR = median(allelic_ratio),
              frac_escaping = mean(allelic_ratio <= 0.9),
              frac_exactly_1 = mean(allelic_ratio == 1),
              .groups = "drop") %>%
    mutate(min_total_reads = cut_i)
})) %>% select(min_total_reads, everything()) %>% arrange(sample, min_total_reads)
write.table(sweep_overall, file.path(OUT_DIR, "ceb_cutoff_sweep_summary_overall.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# The resolution problem, which dominates everything else in this block.
#
# At N reads the ratio can only land on one of N+1 values, so at the observed
# depths most "AR = 1" cells are not evidence of monoallelic expression -- they
# are cells where no other value was reachable. This is also why the escape
# fraction here moves OPPOSITE to whole chrX as the cutoff rises: deeper cells
# can finally express a ratio below 0.9, so apparent escape goes UP.
# ---------------------------------------------------------------------------
resolution <- cells %>%
  filter(total_reads <= 12) %>%
  group_by(total_reads) %>%
  summarise(n_cells = dplyr::n(),
            attainable_values = dplyr::first(total_reads) + 1,
            distinct_AR_seen = dplyr::n_distinct(round(allelic_ratio, 4)),
            frac_exactly_1 = mean(allelic_ratio == 1),
            frac_exactly_0 = mean(allelic_ratio == 0),
            frac_escaping = mean(allelic_ratio <= 0.9),
            mean_AR = mean(allelic_ratio), .groups = "drop")
write.table(resolution, file.path(OUT_DIR, "ceb_resolution_by_depth.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

pdf(file.path(OUT_DIR, "ceb_resolution_limit.pdf"), width = 10, height = 7)
print(
  ggplot(resolution, aes(total_reads)) +
    geom_col(aes(y = n_cells), fill = "grey80") +
    geom_text(aes(y = n_cells, label = paste0(attainable_values, " vals")),
              vjust = -0.4, size = 2.4, colour = "grey30") +
    labs(x = "chrX informative reads over the core escape block", y = "Cells",
         title = "Most cells sit at a depth where the ratio has almost no resolution",
         subtitle = "Label = number of allelic ratios attainable at that depth (N+1)") +
    theme_bw()
)
print(
  ggplot(resolution, aes(total_reads)) +
    geom_line(aes(y = frac_exactly_1, colour = "AR exactly 1.00")) +
    geom_point(aes(y = frac_exactly_1, colour = "AR exactly 1.00")) +
    geom_line(aes(y = frac_exactly_0, colour = "AR exactly 0.00")) +
    geom_point(aes(y = frac_exactly_0, colour = "AR exactly 0.00")) +
    geom_line(aes(y = frac_escaping, colour = "AR <= 0.90 (\"escaping\")")) +
    geom_point(aes(y = frac_escaping, colour = "AR <= 0.90 (\"escaping\")")) +
    scale_colour_manual(values = c("AR exactly 1.00" = "firebrick",
                                   "AR exactly 0.00" = "steelblue",
                                   "AR <= 0.90 (\"escaping\")" = "grey30"), name = NULL) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "chrX informative reads over the core escape block", y = "Fraction of cells",
         title = "Apparent escape rises with depth here -- the opposite of whole chrX",
         caption = paste("Not a biological gradient. At 1-3 reads a cell often cannot express any ratio below 0.9,",
                         "so it is counted as non-escaping\nby arithmetic. As depth rises those cells become able to",
                         "report an intermediate ratio and move into the escaping group.")) +
    theme_bw() + theme(plot.caption = element_text(size = 7, hjust = 0), legend.position = "bottom")
)
dev.off()

# ---------------------------------------------------------------------------
# Depth-stratified and depth-matched sample comparison, as in 06: within a
# band the resolution limit is held fixed, so a sample ordering that persists
# down every band is not an artefact of differing depth distributions.
# ---------------------------------------------------------------------------
CEB_BANDS <- list(c(1, 2), c(2, 3), c(3, 5), c(5, 8), c(8, Inf))
band_label <- sapply(CEB_BANDS, function(b)
  if (is.infinite(b[2])) paste0(b[1], "+") else if (b[2] - b[1] == 1) as.character(b[1])
  else paste0(b[1], "-", b[2] - 1))
cells$depth_band <- factor(band_label[sapply(cells$total_reads, function(x) {
  hit <- which(sapply(CEB_BANDS, function(b) x >= b[1] & x < b[2]))
  if (length(hit)) hit[1] else NA_integer_
})], levels = band_label)

depth_matched <- cells %>%
  filter(!is.na(depth_band)) %>%
  group_by(depth_band, sample) %>%
  summarise(n_cells = dplyr::n(), mean_AR = mean(allelic_ratio),
            frac_escaping = mean(allelic_ratio <= 0.9), .groups = "drop") %>%
  filter(n_cells >= 30)
write.table(depth_matched, file.path(OUT_DIR, "ceb_depth_matched_by_sample.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

pdf(file.path(OUT_DIR, "ceb_cutoff_sweep_trends.pdf"), width = 11, height = 8)
print(
  ggplot(sweep_overall, aes(min_total_reads, n_cells, colour = sample)) +
    geom_line() + geom_point(size = 1.4) +
    scale_x_continuous(breaks = CUTOFFS) + scale_y_log10() +
    labs(x = "Minimum reads over the core escape block", y = "Cells retained (log scale)",
         title = "Cells retained vs cutoff",
         subtitle = "Far steeper than whole chrX: the block's median depth is only a few reads") +
    theme_bw()
)
print(
  ggplot(sweep_overall, aes(min_total_reads, mean_AR, colour = sample)) +
    geom_line() + geom_point(size = 1.4) +
    scale_x_continuous(breaks = CUTOFFS) + coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Minimum reads over the core escape block", y = "Mean allelic ratio",
         title = "Mean allelic ratio vs cutoff",
         subtitle = "Flat lines mean the sample comparison is not being driven by the threshold") +
    theme_bw()
)
print(
  ggplot(sweep_overall, aes(min_total_reads, frac_escaping, colour = sample)) +
    geom_line() + geom_point(size = 1.4) +
    scale_x_continuous(breaks = CUTOFFS) + coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Minimum reads over the core escape block", y = "Fraction AR <= 0.9",
         title = "Apparent escape vs cutoff",
         subtitle = "Rises with the cutoff: see ceb_resolution_limit.pdf -- this is the resolution limit lifting, not biology") +
    theme_bw()
)
print(
  ggplot(depth_matched, aes(depth_band, mean_AR, colour = sample, group = sample)) +
    geom_line() + geom_point(size = 1.6) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Reads over the core escape block", y = "Mean allelic ratio",
         title = "Depth-matched sample comparison",
         caption = "Within a band the resolution limit is fixed, so an ordering that holds down every band is a real sample difference.") +
    theme_bw() + theme(plot.caption = element_text(size = 7, hjust = 0))
)
dev.off()

pdf(file.path(OUT_DIR, "ceb_cutoff_sweep_AR_histograms.pdf"), width = 10, height = 7)
dist_tbl <- bind_rows(lapply(CUTOFFS, function(k) {
  d <- cells[cells$total_reads >= k, ]; d$cutoff_lab <- paste0(">= ", k); d
}))
dist_tbl$cutoff_lab <- factor(dist_tbl$cutoff_lab, levels = paste0(">= ", CUTOFFS))
print(
  ggplot(dist_tbl, aes(allelic_ratio, fill = sample)) +
    geom_histogram(binwidth = 0.05, boundary = 0) +
    facet_grid(cutoff_lab ~ sample, scales = "free_y") +
    labs(x = "Allelic ratio (core escape block)", y = "Cells",
         title = "AR distribution at each cutoff",
         subtitle = "Spikes at 0 and 1 are the resolution limit, not monoallelic biology") +
    theme_bw(base_size = 8) +
    theme(legend.position = "none", strip.text.y = element_text(size = 6))
)
dev.off()

print(resolution)
print(sweep_overall)
message("Done. Output in ", normalizePath(OUT_DIR))
