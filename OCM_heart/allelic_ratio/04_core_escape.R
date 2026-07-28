# ---------------------------------------------------------------------------
# 04 - Core escape gene analysis (block-level, whole chromosome).
#
# Run from the OCM_heart/ directory:  Rscript allelic_ratio/04_core_escape.R
# Requires: 01_setup.R to have been run.
# Writes:   Allelic_ratio_results/core_escape_block_cell_metadata.txt for 05.
# ---------------------------------------------------------------------------
source("allelic_ratio/00_functions.R")

heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

#################################################
# Core escape genes analysis (block, whole-chr) #
#################################################
barcodes <- lapply(c('9w', '78w', 'Sham', 'TAC'), function(x) {
  list.files(paste0('Allelome.PRO2_core_escape/', x), pattern = 'locus_table.txt', recursive = TRUE, full.names = TRUE)
})
barcodes <- unlist(barcodes)

condition <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[6])
cellid <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[7])

core_escape_block_ratio <- lapply(barcodes, function(x) {
  if (file.exists(x)) {
    tmp <- read.delim(x, header = TRUE)
    if (nrow(tmp) == 0) return(NULL)
    tmp$chr <- as.character(tmp$chr)
    tmp <- subset(tmp, chr == "chrX")
    if (nrow(tmp) == 0) return(NULL)
    tmp
  } else {
    NULL
  }
})
names(core_escape_block_ratio) <- paste0(condition, "_", cellid)
core_escape_block_ratio <- core_escape_block_ratio[!sapply(core_escape_block_ratio, is.null)]
core_escape_block_ratio <- bind_rows(core_escape_block_ratio, .id = "cell_barcode")

# Write out the core escape block allelic ratio table
write.table(core_escape_block_ratio, 'Allelic_ratio_results/core_escape_block_allelic_ratio_table.txt', sep = '\t', row.names = FALSE, quote = FALSE)

# Subset seurat object by barcodes
subset_heart_ceb <- subset(heart, cells = core_escape_block_ratio$cell_barcode)
core_escape_block_ratio <- core_escape_block_ratio[core_escape_block_ratio$cell_barcode %in% colnames(subset_heart_ceb), ]
core_escape_block_ratio <- core_escape_block_ratio[match(colnames(subset_heart_ceb), core_escape_block_ratio$cell_barcode), ]
subset_heart_ceb$allelic_ratio <- core_escape_block_ratio$allelic_ratio
subset_heart_ceb$total_reads <- core_escape_block_ratio$total_reads
subset_heart_ceb$A1_reads <- core_escape_block_ratio$A1_reads
subset_heart_ceb$A2_reads <- core_escape_block_ratio$A2_reads

# Plot distribution of total reads
pdf('Allelic_ratio_results/core_escape_block_total_reads_distribution.pdf')
plot(
  density(subset_heart_ceb$total_reads, na.rm = TRUE),
  main = "Total Reads Distribution",
  xlab = "Total Reads"
)
dev.off()

###################################################
# filter for cells with at least 5 reads on chrX #
###################################################
subset_heart_ceb_flt <- subset(subset_heart_ceb, subset = total_reads >= 5)

table(subset_heart_ceb_flt$sample)


# UMAP plot coloured by allelic ratio split by sample
samples_ceb <- levels(subset_heart_ceb_flt$sample)
plots_ceb <- lapply(samples_ceb, function(s) {
  FeaturePlot(subset(subset_heart_ceb_flt, subset = sample == s),
              features = "allelic_ratio",
              min.cutoff = 0,
              max.cutoff = 1) +
    scale_color_gradientn(colors = my_colors,
                          breaks = seq(0, 1, by = 0.1),
                          limits = c(0, 1),
                          oob = scales::squish,
                          name = "Allelic ratio") +
    theme(axis.title.y = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9)) +
    ggtitle(s)
})
pdf('Allelic_ratio_results/core_escape_block_allelic_ratio_umap_plot_split_by_sample.pdf', width = 10, height = 10)
wrap_plots(plots_ceb, ncol = 2, nrow = 2) + plot_layout(guides = "collect")
dev.off()

# Create metadata for statistical testing
metadata_ceb <- subset_heart_ceb_flt@meta.data

# sort by total read count
metadata_ceb <- metadata_ceb[order(metadata_ceb$total_reads, decreasing = TRUE), ]

# number of cells per celltype and condition
cell_counts_ceb <- metadata_ceb %>%
  group_by(celltype, sample) %>%
  summarise(n_cells = n(), .groups = "drop")
write.table(cell_counts_ceb, 'Allelic_ratio_results/core_escape_block_cell_counts_per_celltype_and_condition.txt', sep = '\t', row.names = FALSE, quote = FALSE)

# Statistical testing - is the proportion of putative LOX cells
# (allelic_ratio >= 0.9, i.e. no reads from the inactive/paternal allele of a
# gene expected to escape XCI) higher between conditions for each cell type?
# One-sided by design: under skewed XCI (maternal always active), AR <= 0.1
# would imply loss of the active allele, a distinct and much rarer event.
metadata_ceb$monoallelic <- as.integer(metadata_ceb$allelic_ratio >= 0.9)

adult_vs_aged_monoallelic_lrt_ceb <- split(metadata_ceb, metadata_ceb$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("9w", "78w"))
    if (length(unique(df$sample)) < 2) return(NULL)
    run_monoallelic_lrt(df, ref_level = "9w", comparison_level = "78w")
  })
adult_vs_aged_monoallelic_lrt_ceb <- bind_rows(adult_vs_aged_monoallelic_lrt_ceb, .id = "celltype")
adult_vs_aged_monoallelic_lrt_ceb$FDR <- p.adjust(adult_vs_aged_monoallelic_lrt_ceb$p_value, method = "fdr")

subset(adult_vs_aged_monoallelic_lrt_ceb, FDR < 0.05)

Sham_vs_TAC_monoallelic_lrt_ceb <- split(metadata_ceb, metadata_ceb$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("Sham", "TAC"))
    if (length(unique(df$sample)) < 2) return(NULL)
    run_monoallelic_lrt(df, ref_level = "Sham", comparison_level = "TAC")
  })
Sham_vs_TAC_monoallelic_lrt_ceb <- bind_rows(Sham_vs_TAC_monoallelic_lrt_ceb, .id = "celltype")
Sham_vs_TAC_monoallelic_lrt_ceb$FDR <- p.adjust(Sham_vs_TAC_monoallelic_lrt_ceb$p_value, method = "fdr")

subset(Sham_vs_TAC_monoallelic_lrt_ceb, FDR < 0.05)

# Stacked barplot of the number of putative LOX cells (AR >= 0.9) vs other
# cells, per sample and cell type. Raw counts (N=1 per sample, so proportions
# would overstate precision the data doesn't support).
metadata_ceb$LOX_status <- factor(
  case_when(
    metadata_ceb$allelic_ratio <= 0.10 ~ "LOX-like (AR <= 0.10)",
    metadata_ceb$allelic_ratio >= 0.9 ~ "LOX-like (AR >= 0.9)",
    TRUE ~ "Other"
  ),
  levels = c("LOX-like (AR <= 0.10)", "Other", "LOX-like (AR >= 0.9)")
)
metadata_ceb$sample <- factor(metadata_ceb$sample, levels = c("9w", "78w", "Sham", "TAC"))

# Plotting allelic ratio per cell type and condition
violin_tbl_ceb <- metadata_ceb %>%
  mutate(
    sample = factor(sample, levels = c("9w", "78w", "Sham", "TAC")),
    sample_idx = as.numeric(sample)
  )

pdf('Allelic_ratio_results/core_escape_block_allelic_ratio_celltype_violin_plot_facet_wrap.pdf')
ggplot(violin_tbl_ceb, aes(x = sample_idx, y = allelic_ratio, fill = sample)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.3, color = "black") +
  facet_wrap(~celltype) +
  scale_x_continuous(breaks = 1:4, labels = levels(violin_tbl_ceb$sample)) +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9, 1.0)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Allelic ratio", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(5.5, 18, 5.5, 5.5),
        legend.position = "none")
dev.off()

pdf('Allelic_ratio_results/core_escape_block_LOX_cell_counts_stacked_barplot.pdf')
ggplot(metadata_ceb, aes(x = sample, fill = LOX_status)) +
  geom_bar(position = "fill") +
  facet_wrap(~celltype) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Proportion of cells", x = NULL, fill = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
dev.off()

# Save file 
saveRDS(subset_heart_ceb, 'Allelic_ratio_results/subset_heart_core_escape_block.RDS')


# Correlation of Xist expression and AR
AR_Xist_df <- data.frame(
  Xist = GetAssayData(subset_heart_ceb_flt, assay = "SCT", layer = "data")["Xist", ],
  allelic_ratio = subset_heart_ceb_flt$allelic_ratio,
  sample = subset_heart_ceb_flt$sample,
  celltype = subset_heart_ceb_flt$celltype
)

AR_Xist_cor <- split(AR_Xist_df, list(AR_Xist_df$celltype, AR_Xist_df$sample), drop = TRUE) %>%
  lapply(function(x) {
    if (nrow(x) < 3) return(NULL)
    tmp <- cor.test(x$Xist, x$allelic_ratio, method = "spearman")
    data.frame(sample = unique(x$sample), celltype = unique(x$celltype), cor = tmp$estimate, p_value = tmp$p.value)
  }) %>%
  bind_rows()
AR_Xist_cor$FDR <- p.adjust(AR_Xist_cor$p_value, method = "fdr")

# Scatter of Xist expression vs allelic ratio, per celltype and sample
ann_df <- AR_Xist_cor %>%
  mutate(star = fdr_to_stars(FDR),
         label = paste0(sample, ": rho=", round(cor, 2), " ", star)) %>%
  group_by(celltype) %>%
  summarise(label = paste(label, collapse = "\n"), .groups = "drop")

pdf('Allelic_ratio_results/core_escape_block_Xist_vs_AR_scatter.pdf', width = 10, height = 10)
ggplot(AR_Xist_df, aes(x = allelic_ratio, y = Xist, color = sample)) +
  geom_point(size = 0.6, alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  geom_text(data = ann_df, aes(x = 0.02, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1.05, size = 2.2) +
  facet_wrap(~celltype) +
  labs(x = "Allelic ratio (core escape genes)", y = "Xist expression (SCT, normalized)", color = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()

# Same, subset to Ventricular Cardiomyocytes only, faceted by sample
vcm_xist_df <- AR_Xist_df %>%
  filter(celltype == "Ventricular Cardiomyocytes")

vcm_ann <- AR_Xist_cor %>%
  filter(celltype == "Ventricular Cardiomyocytes") %>%
  mutate(label = paste0("rho=", round(cor, 2), " ", fdr_to_stars(FDR)))

pdf('Allelic_ratio_results/core_escape_block_Xist_vs_AR_scatter_VCM.pdf', width = 10, height = 3.5)
ggplot(vcm_xist_df, aes(x = allelic_ratio, y = Xist)) +
  geom_point(size = 0.8, alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  geom_text(data = vcm_ann, aes(x = 0.02, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1.2, size = 3) +
  facet_wrap(~sample, nrow = 1) +
  labs(x = "Allelic ratio (core escape genes)", y = "Xist expression\n(SCT, normalized)",
       title = "Ventricular Cardiomyocytes") +
  theme_bw()
dev.off()











# Per-gene breakdown for cells flagged as putative LOX (block-level AR >= 0.9),
# to eyeball concordance across the 4 core escape genes before committing to a
# formal multi-gene threshold.
lox_barcodes <- rownames(subset(metadata_ceb, monoallelic == 1))

gene_df_percell <- read.table('Allelic_ratio_results/core_escape_genes_gene_df.txt',
                               sep = '\t', header = TRUE, stringsAsFactors = FALSE)

lox_gene_breakdown <- gene_df_percell %>%
  dplyr::filter(cell_barcode %in% lox_barcodes) %>%
  dplyr::select(cell_barcode, sample, celltype, name, allelic_ratio, A1_reads, A2_reads, total_reads) %>%
  dplyr::arrange(cell_barcode, name)

write.table(lox_gene_breakdown,
            'Allelic_ratio_results/core_escape_block_LOX_cells_per_gene_breakdown.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

lox_gene_breakdown_wide <- reshape(
  as.data.frame(lox_gene_breakdown),
  idvar = c("cell_barcode", "sample", "celltype"),
  timevar = "name",
  direction = "wide"
)

write.table(lox_gene_breakdown_wide,
            'Allelic_ratio_results/core_escape_block_LOX_cells_per_gene_breakdown_wide.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

# Confirm LOX status: require AR >= 0.9 independently at all four core escape
# genes (missing coverage at any gene fails this - only cells with full
# 4-gene support are counted as confirmed).
gene_cols <- grep("^allelic_ratio\\.", colnames(lox_gene_breakdown_wide), value = TRUE)

lox_gene_breakdown_wide$confirmed_LOX <- apply(
  lox_gene_breakdown_wide[, gene_cols, drop = FALSE], 1,
  function(x) all(!is.na(x) & x >= 0.9)
)
sum(lox_gene_breakdown_wide$confirmed_LOX)

write.table(lox_gene_breakdown_wide,
            'Allelic_ratio_results/core_escape_block_LOX_cells_per_gene_breakdown_wide.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

confirmed_lox_barcodes <- lox_gene_breakdown_wide$cell_barcode[lox_gene_breakdown_wide$confirmed_LOX]

# Handoff to 05_lox_sensitivity.R. Written at the END of this script rather
# than where metadata_ceb is first built, because the monoallelic column that
# 05 filters on is not added until line ~112 above.
write.table(metadata_ceb,
            'Allelic_ratio_results/core_escape_block_cell_metadata.txt',
            sep = '\t', quote = FALSE, col.names = NA)

