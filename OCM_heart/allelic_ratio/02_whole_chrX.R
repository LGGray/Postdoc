# ---------------------------------------------------------------------------
# 02 - Whole X chromosome analysis.
#
# Run from the OCM_heart/ directory:  Rscript allelic_ratio/02_whole_chrX.R
# Requires: 01_setup.R to have been run (needs celltype_sub on the object)
# Writes:   Allelic_ratio_results/whole_chr_cell_metadata.txt, which 03 and 05
#           read back so they can run as independent jobs.
# ---------------------------------------------------------------------------
source("allelic_ratio/00_functions.R")

heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

###############################
# Whole X chromosome analysis #
################################
# barcodes <- lapply(c('9w', '78w', 'Sham', 'TAC'), function(x) {
#   list.files(paste0('Allelome.PRO2/', x), pattern = 'locus_table.txt', recursive = TRUE, full.names = TRUE)
# })
# barcodes <- unlist(barcodes)

# condition <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[4])
# cellid <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[5])

# chr_allelic_ratio <- lapply(barcodes, function(x) {
#   if (file.exists(x)) {
#     tmp <- read.delim(x, header = TRUE)
#     if (nrow(tmp) == 0) return(NULL)
#     tmp$chr <- as.character(tmp$chr)
#     # tmp <- subset(tmp, chr=="chrX")
#     if (nrow(tmp) == 0) return(NULL)
#     tmp
#   } else {
#     NULL
#   }
# })
# names(chr_allelic_ratio) <- paste0(condition, "_", cellid)
# chr_allelic_ratio <- chr_allelic_ratio[!sapply(chr_allelic_ratio, is.null)]
# chr_allelic_ratio <- bind_rows(chr_allelic_ratio, .id = "cell_barcode")

# # save table of allelic ratios for whole chrX
# write.table(chr_allelic_ratio, 'Allelic_ratio_results/whole_chr_allelic_ratios.txt', sep = '\t', row.names = FALSE, quote = FALSE)
chr_allelic_ratio <- read.table('Allelic_ratio_results/whole_chr_allelic_ratios.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

#### IMPORTANT TO FILTER FOR chrX #######
chr_allelic_ratio <- subset(chr_allelic_ratio, chr == "chrX")

# Subset seurat object by barcodes
subset_heart <- subset(heart, cells = chr_allelic_ratio$cell_barcode)
chr_allelic_ratio <- chr_allelic_ratio[chr_allelic_ratio$cell_barcode %in% colnames(subset_heart), ]
chr_allelic_ratio <- chr_allelic_ratio[match(colnames(subset_heart), chr_allelic_ratio$cell_barcode), ]
subset_heart$allelic_ratio <- chr_allelic_ratio$allelic_ratio
subset_heart$total_reads <- chr_allelic_ratio$total_reads
subset_heart$A1_reads <- chr_allelic_ratio$A1_reads
subset_heart$A2_reads <- chr_allelic_ratio$A2_reads


# Plot distribution of total reads
pdf('Allelic_ratio_results/whole_chr_total_reads_distribution.pdf')
plot(
  density(subset_heart$total_reads, na.rm = TRUE),
  main = "Total Reads Distribution",
  xlab = "Total Reads"
)
abline(v = 10, col = "red", lty = 2)
dev.off()
#
################################################### 
# filter for cells with at least 10 reads on chrX #
###################################################
subset_heart_flt <- subset(subset_heart, subset = total_reads >= 10)

pdf('Allelic_ratio_results/whole_chr_UMAP_celltypes.pdf')
DimPlot(subset_heart_flt, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3) +
  theme(legend.position = "none")
dev.off()

# UMAP plot coloured by allelic ratio split by sample
my_breaks <- seq(0, 1, by = 0.05)

my_colors <- c(
  "#2B3186", "#374795", "#3B5FB6", "#38749F",
  "#37758B", "#367373", "#2D6E5D", "#2A7050",
  "#1E652D", "#1C642D", "#0F7031", "#357B30", "#4E8330",
  "#658C2D", "#78962A", "#8D9F25", "#A2A71D", "#B3B112",
  "#C97314", "#8B1913"
)
my_breaks  <- seq(0, 1, by = 0.05)
bin_labels <- sprintf("%.2f–%.2f", head(my_breaks, -1), tail(my_breaks, -1))

subset_heart_flt$allelic_bin <- cut(
  subset_heart_flt$allelic_ratio,
  breaks = my_breaks,
  include.lowest = TRUE,
  right = TRUE,
  labels = bin_labels        # ensures all 20 levels exist
)


my_breaks <- c(seq(0, 0.9, by = 0.1), 0.95, 1.0)

my_colors <- c(
  "#2B3186", "#3B5FB6", "#38749F",
  "#367373", "#2D6E5D",
  "#1E652D",
  "#658C2D","#8D9F25", "#B3B112",
  "#C97314", "#8B1913"
)

bin_labels <- c(
  "0.00–0.10",
  "0.10–0.20",
  "0.20–0.30",
  "0.30–0.40",
  "0.40–0.50",
  "0.50–0.60",
  "0.60–0.70",
  "0.70–0.80",
  "0.80–0.90",
  "0.90–0.95",
  "0.95–1.00"
)

subset_heart_flt$allelic_bin <- cut(
  subset_heart_flt$allelic_ratio,
  breaks = my_breaks,
  include.lowest = TRUE,
  right = TRUE,
  labels = bin_labels
)


# plot each sample seperately
p_9w <- DimPlot(
  subset(subset_heart_flt, subset = sample == "9w"),
  reduction = "umap",
  group.by = "allelic_bin"
) + ggtitle("9w")

p_78w <- DimPlot(
  subset(subset_heart_flt, subset = sample == "78w"),
  reduction = "umap",
  group.by = "allelic_bin"
) + ggtitle("78w")

p_Sham <- DimPlot(
  subset(subset_heart_flt, subset = sample == "Sham"),
  reduction = "umap",
  group.by = "allelic_bin"
) + ggtitle("Sham")

p_TAC <- DimPlot(
  subset(subset_heart_flt, subset = sample == "TAC"),
  reduction = "umap",
  group.by = "allelic_bin"
) + ggtitle("TAC")

# Plot all samples together with patchwork
library(patchwork)
p_all_samples <- (p_9w | p_78w) / (p_Sham | p_TAC) &
  scale_color_manual(
    values   = my_colors,
    drop     = FALSE,
    na.value = "grey70"
  ) &
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  filename = file.path("Allelic_ratio_results/allelic_ratio_umap_plot_split_by_sample.pdf"),
  plot = p_all_samples,
  width = 10,
  height = 7
)

# Distribution of allelic ratio split by sample -- a plainer complement to the
# UMAP color scheme above, using the same 0.1-wide bins
ar_dist_df <- data.frame(
  allelic_ratio = subset_heart_flt$allelic_ratio,
  sample = factor(subset_heart_flt$sample, levels = c("9w", "78w", "Sham", "TAC"))
)

pdf('Allelic_ratio_results/allelic_ratio_distribution_by_sample.pdf', width = 12, height = 4)
ggplot(ar_dist_df, aes(x = allelic_ratio, fill = sample)) +
  geom_histogram(binwidth = 0.1, boundary = 0) +
  facet_wrap(~sample, nrow = 1) +
  labs(x = "Allelic ratio", y = "Number of cells", fill = NULL) +
  theme_bw() +
  guides(fill = "none")
dev.off()

pdf('Allelic_ratio_results/allelic_ratio_vs_n_snp_reads.pdf', width = 10, height = 10)
ggplot(subset_heart_flt@meta.data, aes(total_reads, allelic_ratio)) +
  geom_point(size = 0.4, alpha = 0.3) +
  geom_hline(yintercept = 0.5, colour = "#8b1913", linetype = 2) +
  scale_x_log10() + facet_wrap(~ sample)
dev.off()

with(subset_heart_flt@meta.data, table(sample, below_half = allelic_ratio < 0.5))


# Create metadata for statistical testing
metadata_whole_chr <- subset_heart_flt@meta.data

# Handoff to 03_all_genes.R and 05_lox_sensitivity.R. This is the one shared
# object that was not already being written to disk, so without it those two
# scripts would only work inside a single long session.
write.table(metadata_whole_chr,
            'Allelic_ratio_results/whole_chr_cell_metadata.txt',
            sep = '\t', quote = FALSE, col.names = NA)

# number of cells per celltype and condition
cell_counts <- metadata_whole_chr %>%
  group_by(celltype, sample) %>%
  summarise(n_cells = n(), .groups = "drop")
write.table(cell_counts, 'Allelic_ratio_results/whole_chr_cell_counts_per_celltype_and_condition.txt', sep = '\t', row.names = FALSE, quote = FALSE)

cell_ids <- data.frame(barcode = colnames(subset_heart_flt), cell_id = colnames(subset_heart_flt), age = subset_heart_flt$sample, celltype = subset_heart_flt$celltype)
cell_ids$sample <- factor(cell_ids$age, levels = c("TAC", "Sham", "78w", "9w"))
cell_ids$barcode <- gsub('9w_|78w_|Sham_|TAC_', '', cell_ids$barcode)
pdf('Allelic_ratio_results/whole_chr_cell_counts_per_celltype_and_condition.pdf', width = 10, height = 10)
ggplot(cell_ids, aes(x = sample, fill = celltype)) +
  geom_bar(position = "fill") +
  labs(x = "", y = "Cell composition", title = "") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, by = 0.25)) +
  coord_flip() +
  theme(legend.position = "bottom", legend.title = element_blank())
dev.off()

# Statistical testing - Does the variation in allelic ratio differ between conditions for each cell type?
adult_vs_aged_lrt <- split(metadata_whole_chr, metadata_whole_chr$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("9w", "78w"))
    if (length(unique(df$sample)) < 2) return(NULL)
    run_dispersion_lrt(df, ref_level = "9w", comparison_level = "78w")
  })
adult_vs_aged_lrt <- bind_rows(adult_vs_aged_lrt, .id = "celltype")
adult_vs_aged_lrt$FDR <- p.adjust(adult_vs_aged_lrt$Pr..Chisq., method = "fdr")

subset(adult_vs_aged_lrt, FDR < 0.05)

# Quantify the n=1-per-condition caveat for the celltypes flagged above: if
# 9w and 78w truly had the SAME dispersion, how often would this exact LRT
# still call it significant, purely from animal-to-animal noise? Swept across
# a few plausible animal-effect magnitudes since that value can't be
# estimated from n=1 data -- treat this as a range, not a single number.
adult_vs_aged_sig_celltypes <- subset(adult_vs_aged_lrt, FDR < 0.05)$celltype

adult_vs_aged_fpr <- sapply(adult_vs_aged_sig_celltypes, function(ct) {
  df <- subset(metadata_whole_chr, celltype == ct & sample %in% c("9w", "78w"))
  sapply(c(0.15, 0.3, 0.5), function(sd) {
    simulate_dispersion_null_fpr(df, "9w", "78w", animal_disp_sd = sd, n_reps = 1000)$fpr
  })
})
rownames(adult_vs_aged_fpr) <- c("fpr_animal_sd_0.15", "fpr_animal_sd_0.3", "fpr_animal_sd_0.5")
t(adult_vs_aged_fpr)

Sham_vs_TAC_lrt <- split(metadata_whole_chr, metadata_whole_chr$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("Sham", "TAC"))
    if (length(unique(df$sample)) < 2) return(NULL)
    run_dispersion_lrt(df, ref_level = "Sham", comparison_level = "TAC")
  })
Sham_vs_TAC_lrt <- bind_rows(Sham_vs_TAC_lrt, .id = "celltype")
Sham_vs_TAC_lrt$FDR <- p.adjust(Sham_vs_TAC_lrt$Pr..Chisq., method = "fdr")

subset(Sham_vs_TAC_lrt, FDR < 0.05)

# Same false-positive-rate calibration as above, for Sham vs TAC
Sham_vs_TAC_sig_celltypes <- subset(Sham_vs_TAC_lrt, FDR < 0.05)$celltype

Sham_vs_TAC_fpr <- sapply(Sham_vs_TAC_sig_celltypes, function(ct) {
  df <- subset(metadata_whole_chr, celltype == ct & sample %in% c("Sham", "TAC"))
  sapply(c(0.15, 0.3, 0.5), function(sd) {
    simulate_dispersion_null_fpr(df, "Sham", "TAC", animal_disp_sd = sd, n_reps = 1000)$fpr
  })
})
rownames(Sham_vs_TAC_fpr) <- c("fpr_animal_sd_0.15", "fpr_animal_sd_0.3", "fpr_animal_sd_0.5")
t(Sham_vs_TAC_fpr)

# Effect sizes: fraction of cells escaping XCI (AR <= 0.9) and median AR per
# celltype/condition. Complements the dispersion LRT p-values above with an
# interpretable magnitude.
AR_COL <- "allelic_ratio"


frac_escaping <- metadata_whole_chr %>%
  rename(AR = all_of(AR_COL)) %>%
  { if ("total_reads" %in% names(.)) filter(., total_reads >= MIN_READS) else . } %>%
  mutate(sample = factor(sample, levels = c("9w", "78w", "Sham", "TAC"))) %>%
  group_by(celltype, sample) %>%
  summarise(n = n(),
            escaping = mean(AR <= 0.9),
            median_AR = median(AR),
            .groups = "drop") %>%
  filter(n >= 20)

write.table(frac_escaping, 'Allelic_ratio_results/whole_chr_fraction_escaping_per_celltype_and_condition.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

pdf('Allelic_ratio_results/whole_chr_fraction_escaping_barplot.pdf')
ggplot(frac_escaping, aes(sample, 100*escaping, fill = sample)) +
  geom_col() +
  facet_wrap(~celltype, labeller = label_wrap_gen(14)) +
  labs(x = NULL, y = "Escaping cells (%)  [AR ≤ 0.9]", fill = "Condition") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())
dev.off()

violin_tbl <- metadata_whole_chr  %>%
  mutate(
    sample = factor(sample, levels = c("9w", "78w", "Sham", "TAC")),
    sample_idx = as.numeric(sample)
  )

violin_ymax <- violin_tbl %>%
  group_by(celltype) %>%
  summarise(y_top = max(allelic_ratio, na.rm = TRUE), .groups = "drop")

violin_ann <- bind_rows(
  adult_vs_aged_lrt %>%
    transmute(celltype, comp = "aa", x1 = 1, x2 = 2, star = fdr_to_stars(FDR)),
  Sham_vs_TAC_lrt %>%
    transmute(celltype, comp = "st", x1 = 3, x2 = 4, star = fdr_to_stars(FDR))
) %>%
  inner_join(violin_ymax, by = "celltype") %>%
  mutate(comp = factor(comp, levels = c("aa", "st"))) %>%
  group_by(celltype) %>%
  arrange(comp, .by_group = TRUE) %>%
  mutate(
    y = pmin(1.10, y_top + 0.04 + (row_number() - 1) * 0.07),
    y_tick = y - 0.02,
    x_mid = (x1 + x2) / 2
  ) %>%
  ungroup()

pdf('Allelic_ratio_results/whole_chr_allelic_ratio_celltype_violin_plot_facet_wrap.pdf')
ggplot(violin_tbl, aes(x = sample_idx, y = allelic_ratio, fill = sample)) +
  geom_violin(trim = FALSE, scale = "width", bounds = c(0, 1)) +
  geom_segment(data = violin_ann,
               aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE,
               color = "black",
               linewidth = 0.3) +
  geom_segment(data = violin_ann,
               aes(x = x1, xend = x1, y = y_tick, yend = y),
               inherit.aes = FALSE,
               color = "black",
               linewidth = 0.3) +
  geom_segment(data = violin_ann,
               aes(x = x2, xend = x2, y = y_tick, yend = y),
               inherit.aes = FALSE,
               color = "black",
               linewidth = 0.3) +
  geom_text(data = violin_ann,
            aes(x = x_mid, y = y + 0.01, label = star),
            inherit.aes = FALSE,
            color = "black",
            size = 3.2,
            vjust = 0) +
  facet_wrap(~celltype, labeller = label_wrap_gen(width = 18)) +
  scale_x_continuous(breaks = 1:4, labels = levels(violin_tbl$sample)) +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9, 1.0)) +
  coord_cartesian(ylim = c(0, 1.12), clip = "off") +
  labs(y = "Allelic ratio", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(5.5, 18, 5.5, 5.5),
        legend.position = "none")
dev.off()

prop_tbl <- subset_heart_flt@meta.data %>%
  group_by(celltype, sample) %>%
  summarise(
    n = n(),
    n_biallelic = sum(allelic_ratio < 0.9, na.rm = TRUE),
    p_biallelic = n_biallelic / n,
    se = sqrt(p_biallelic * (1 - p_biallelic) / n),
    .groups = "drop"
  )

pdf('Allelic_ratio_results/whole_chr_fraction_biallelic_cells_per_celltype_and_condition.pdf')
ggplot(prop_tbl, aes(x = sample, y = p_biallelic, color = sample)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, p_biallelic - 1.96*se),
                    ymax = pmin(1, p_biallelic + 1.96*se)),
                width = 0.2) +
  facet_wrap(~celltype) +
  labs(y = "Fraction biallelic-like cells", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")
dev.off()

saveRDS(subset_heart_flt, file = "Allelic_ratio_results/whole_chr_subset_heart_flt.RDS")

# Cell IDs for the highest-read-count Ventricular Cardiomyocyte per sample (for IGV inspection)
metadata_whole_chr$cell_barcode <- rownames(metadata_whole_chr)

cm_cells <- metadata_whole_chr %>%
  filter(celltype == "Ventricular Cardiomyocytes")

top_reads_cells <- cm_cells %>%
  group_by(sample) %>%
  slice_max(total_reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(selection = "top_reads")

biallelic_cells <- cm_cells %>%
  filter(allelic_ratio < 0.7) %>%
  group_by(sample) %>%
  slice_max(total_reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(selection = "biallelic")

top_cm_cells <- bind_rows(top_reads_cells, biallelic_cells) %>%
  mutate(barcode_stripped = gsub('9w_|78w_|Sham_|TAC_', '', cell_barcode)) %>%
  select(sample, selection, cell_barcode, barcode_stripped, total_reads, allelic_ratio)

write.table(top_cm_cells, 'Allelic_ratio_results/top_ventricular_cardiomyocyte_per_sample.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

top_cm_cells

