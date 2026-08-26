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
# How much of the "escape" at a given cutoff is just binomial noise?
#
# The trend plots show the escape fraction falling as the cutoff rises, but on
# their own they cannot say whether the low-cutoff number was noise or whether
# the high cutoff is selecting a different population. This does separate them.
#
# Take only the deep cells, where AR is measured precisely, and thin each one's
# own reads down to a shallow depth by hypergeometric sampling. That holds the
# cells fixed and changes only the read count, so whatever escape appears is
# noise by construction. If thinned deep cells reproduce the escape seen in real
# cells at that depth, the low-cutoff signal was measurement error; if real
# shallow cells still show more, something depth-associated is left over
# (ambient contamination is the leading candidate, cf. the depth analysis in 02).
# ---------------------------------------------------------------------------
DEEP_MIN   <- 100     # "measured accurately" reference set
THIN_DEPTH <- c(10, 20, 30, 50, 75)
THIN_REPS  <- 200

set.seed(1)
deep <- cells[cells$total_reads >= DEEP_MIN, ]

thin_escape <- bind_rows(lapply(SAMPLE_LEVELS, function(s) {
  d <- deep[deep$sample == s, ]
  if (nrow(d) < 30) return(NULL)
  bind_rows(lapply(THIN_DEPTH, function(k) {
    k_cell <- pmin(k, d$total_reads)      # a cell cannot give up more reads than it has
    reps <- replicate(THIN_REPS, {
      a <- rhyper(nrow(d), d$A1_reads, d$A2_reads, k_cell)
      mean(a / k_cell <= 0.9)
    })
    data.frame(sample = s, depth = k, n_deep = nrow(d),
               thinned_escape = mean(reps),
               thinned_lwr = quantile(reps, 0.025, names = FALSE),
               thinned_upr = quantile(reps, 0.975, names = FALSE))
  }))
}))

# Matched comparison: real cells at the same nominal depth, i.e. the escape
# fraction actually observed once that cutoff is applied.
observed_escape <- bind_rows(lapply(THIN_DEPTH, function(k) {
  cells[cells$total_reads >= k, ] %>%
    group_by(sample) %>%
    summarise(observed_escape = mean(allelic_ratio <= 0.9),
              n_cells = dplyr::n(), .groups = "drop") %>%
    mutate(depth = k)
}))

thin_cmp <- thin_escape %>%
  mutate(sample = factor(sample, levels = SAMPLE_LEVELS)) %>%
  left_join(observed_escape, by = c("sample", "depth")) %>%
  mutate(escape_at_deep = sapply(as.character(sample), function(s)
           mean(deep$allelic_ratio[deep$sample == s] <= 0.9)),
         # what share of the observed escape the noise floor alone accounts for
         noise_share = thinned_escape / observed_escape) %>%
  arrange(sample, depth)

write.table(thin_cmp, file.path(OUT_DIR, "cutoff_sweep_noise_floor.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

pdf(file.path(OUT_DIR, "cutoff_sweep_noise_floor.pdf"), width = 10, height = 5)
print(
  ggplot(thin_cmp, aes(depth)) +
    geom_ribbon(aes(ymin = thinned_lwr, ymax = thinned_upr),
                fill = "grey75", alpha = 0.6) +
    geom_line(aes(y = thinned_escape, colour = "Noise floor (deep cells, reads thinned)")) +
    geom_point(aes(y = thinned_escape, colour = "Noise floor (deep cells, reads thinned)"),
               size = 1.2) +
    geom_line(aes(y = observed_escape, colour = "Observed (real cells at this cutoff)")) +
    geom_point(aes(y = observed_escape, colour = "Observed (real cells at this cutoff)"),
               size = 1.2) +
    geom_hline(aes(yintercept = escape_at_deep), linetype = "dotted") +
    facet_wrap(~sample, nrow = 1) +
    scale_x_log10(breaks = THIN_DEPTH) +
    scale_colour_manual(values = c("Noise floor (deep cells, reads thinned)" = "grey35",
                                   "Observed (real cells at this cutoff)" = "firebrick"),
                        name = NULL) +
    labs(x = "Minimum chrX total_reads", y = "Fraction with AR <= 0.9",
         title = "How much apparent escape survives once binomial noise is accounted for",
         caption = paste("Grey: cells with >=", DEEP_MIN, "reads, thinned to the depth on the x-axis --",
                         "escape that is measurement error by construction.\nRed: escape actually observed",
                         "at that cutoff. Dotted: escape among the deep cells at full depth.",
                         "The gap between grey and red is the part\nthe noise floor does not explain.")) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          plot.caption = element_text(size = 7, hjust = 0))
)
dev.off()

print(thin_cmp)

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
