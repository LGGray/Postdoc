# ---------------------------------------------------------------------------
# 02 - Whole X chromosome analysis.
#
# Run from the OCM_heart/ directory:  Rscript allelic_ratio/02_whole_chrX.R
# Requires: 01_setup.R to have been run (needs celltype_sub on the object)
# Writes:   Allelic_ratio_results/whole_chr_cell_metadata.txt, which 03 and 05
#           read back so they can run as independent jobs.
# ---------------------------------------------------------------------------
source("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/OCM_heart/allelic_ratio/00_functions.R")

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
chr_allelic_ratio <- read.table(ALLELIC_RATIOS_FILE, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

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
pdf(file.path(CUTOFF_DIR, 'whole_chr_total_reads_distribution.pdf'))
plot(
  density(subset_heart$total_reads, na.rm = TRUE),
  main = "Total Reads Distribution",
  xlab = "Total Reads"
)
abline(v = MIN_TOTAL_READS, col = "red", lty = 2)
dev.off()
#
#####################################################
# filter for cells with at least MIN_TOTAL_READS reads
# on chrX (see 00_functions.R; one output dir per value)
#####################################################
subset_heart_flt <- subset(subset_heart, subset = total_reads >= MIN_TOTAL_READS)

pdf(file.path(CUTOFF_DIR, 'whole_chr_UMAP_celltypes.pdf'))
DimPlot(subset_heart_flt, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3) +
  theme(legend.position = "none")
dev.off()

# UMAP plot coloured by allelic ratio split by sample
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
  labels = bin_label
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
  filename = file.path(CUTOFF_DIR, "allelic_ratio_umap_plot_split_by_sample.pdf"),
  plot = p_all_samples,
  width = 10,
  height = 7
)


# Create metadata for statistical testing
metadata_whole_chr <- subset_heart_flt@meta.data

# Handoff to 03_all_genes.R and 05_lox_sensitivity.R. This is the one shared
# object that was not already being written to disk, so without it those two
# scripts would only work inside a single long session.
write.table(metadata_whole_chr,
            file.path(CUTOFF_DIR, 'whole_chr_cell_metadata.txt'),
            sep = '\t', quote = FALSE, col.names = NA)

# number of cells per celltype and condition
cell_counts <- metadata_whole_chr %>%
  group_by(celltype, sample) %>%
  summarise(n_cells = n(), .groups = "drop")
write.table(cell_counts, file.path(CUTOFF_DIR, 'whole_chr_cell_counts_per_celltype_and_condition.txt'), sep = '\t', row.names = FALSE, quote = FALSE)

cell_ids <- data.frame(barcode = colnames(subset_heart_flt), cell_id = colnames(subset_heart_flt), age = subset_heart_flt$sample, celltype = subset_heart_flt$celltype)
cell_ids$sample <- factor(cell_ids$age, levels = c("TAC", "Sham", "78w", "9w"))
cell_ids$barcode <- gsub('9w_|78w_|Sham_|TAC_', '', cell_ids$barcode)
pdf(file.path(CUTOFF_DIR, 'whole_chr_cell_counts_per_celltype_and_condition.pdf'), width = 10, height = 10)
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

# Save output 
write.table(adult_vs_aged_lrt, file.path(CUTOFF_DIR, 'whole_chr_adult_vs_aged_dispersion_LRT_results.txt'), sep = '\t', row.names = FALSE, quote = FALSE)

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

write.table(Sham_vs_TAC_lrt, file.path(CUTOFF_DIR, 'whole_chr_Sham_vs_TAC_dispersion_LRT_results.txt'), sep = '\t', row.names = FALSE, quote = FALSE)

# Effect sizes: fraction of cells escaping XCI (AR <= 0.9) and median AR per
# celltype/condition. Complements the dispersion LRT p-values above with an
# interpretable magnitude.
AR_COL <- "allelic_ratio"


frac_escaping <- metadata_whole_chr %>%
  rename(AR = all_of(AR_COL)) %>%
  { if ("total_reads" %in% names(.)) filter(., total_reads >= MIN_TOTAL_READS) else . } %>%
  mutate(sample = factor(sample, levels = c("9w", "78w", "Sham", "TAC"))) %>%
  group_by(celltype, sample) %>%
  summarise(n = n(),
            escaping = mean(AR <= 0.9),
            median_AR = median(AR),
            .groups = "drop") %>%
  filter(n >= 20)

write.table(frac_escaping, file.path(CUTOFF_DIR, 'whole_chr_fraction_escaping_per_celltype_and_condition.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)

pdf(file.path(CUTOFF_DIR, 'whole_chr_fraction_escaping_barplot.pdf'))
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

pdf(file.path(CUTOFF_DIR, 'whole_chr_allelic_ratio_celltype_violin_plot_facet_wrap.pdf'))
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

pdf(file.path(CUTOFF_DIR, 'whole_chr_fraction_biallelic_cells_per_celltype_and_condition.pdf'))
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

saveRDS(subset_heart_flt, file = file.path(CUTOFF_DIR, "whole_chr_subset_heart_flt.RDS"))

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

write.table(top_cm_cells, file.path(CUTOFF_DIR, 'top_ventricular_cardiomyocyte_per_sample.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)

top_cm_cells


# ###########################################################################
# NEGATIVE RESULT - DO NOT INTERPRET THE OUTPUTS OF THIS BLOCK.
#
# Ran 2026-07-29. The conclusion is that this dataset cannot answer the
# cell cycle question, for a reason that is a property of the tissue rather
# than of the code:
#
#   whole_chrX_tricycle_embedding.pdf shows a BLOB, not a ring. Real cell
#   cycle signal makes cells trace a circle in the tricycle embedding. Here
#   almost every cell sits on the origin, so there is no cycling population
#   to detect - as expected for adult myocardium.
#
# Consequences, all of which were visible in the outputs:
#   - theta is atan2() of each cell's position in that embedding. For a cell
#     at the origin the angle is set purely by noise, so theta is essentially
#     uniform on [0, 2*pi] and unrelated to cell cycle. This is why the UMAP
#     coloured by theta is confetti with no spatial structure.
#   - estimate_Schwabe_stage() still assigns all five stages to large numbers
#     of cells, which is not credible in a post-mitotic tissue. Those stage
#     calls track total RNA content, not cell cycle.
#   - AR declines monotonically across the Schwabe stage ordering in all four
#     samples. That trend is real and reproducible but is NOT cell cycle: it
#     is the library-size / ambient RNA confound (see the depth block at the
#     end of this file). A monotone decline is also the wrong shape for
#     replication timing, which predicts a peak in S returning to baseline
#     in G2.
#   - A handful of harmonic fits reach FDR < 0.05, concentrated in Sham,
#     which has the most cells. That is a power signature acting on the depth
#     confound, not evidence of a cell cycle effect.
#
# The code is kept because the negative is informative and the diagnostic
# (embedding ring vs blob) is the right first check if this is ever revisited
# in a proliferative tissue or dataset. Do not re-run and report these
# figures as a cell cycle result.
# ###########################################################################

# ===========================================================================
# Cell cycle position (tricycle) vs allelic ratio.
#
# tricycle projects each cell onto a reference cell cycle embedding and returns
# a continuous position theta in [0, 2*pi]. Per the package docs:
#   0.5*pi          start of S phase
#   pi              start of G2M
#   1.5*pi          middle of M
#   1.75*pi - 0.25*pi   G1/G0
#
# Motivation: the Xi replicates late in S phase, the Xa early. So there is a
# window in mid-to-late S where the Xa is already at two DNA copies and the Xi
# is still at one. If transcription scales with template copy number, that
# transient 2:1 dosage imbalance should push the whole-chrX allelic ratio UP
# during S, and it should return to baseline in G2 once the Xi has replicated
# too. The prediction is therefore specific: a peak in AR at theta between
# 0.5*pi and pi, not a general association with cycling.
#
# Note this is a dosage artefact of replication timing, not a change in escape.
# Escape genes are the only Xi contribution to chrX AR, so if the effect is
# real it means apparent escape is partly a function of where a cell sits in S
# phase - which would be a confounder worth controlling for elsewhere in the
# analysis, not a biological result about escape itself.
#
# READ BEFORE INTERPRETING:
#  - Ventricular cardiomyocytes are post-mitotic. With no S phase there is no
#    replication-timing effect to detect, so a null result there is expected
#    and uninformative. The cell types that can actually answer this are
#    fibroblasts, macrophages and endothelial cells.
#  - The readout is steady-state RNA, not nascent transcription. mRNA
#    half-lives are of the same order as S phase, so a transient doubling of
#    template will be heavily buffered: expect the effect to be attenuated and
#    phase-delayed relative to the DNA-level prediction. If this is snRNA-seq
#    the buffering is weaker and the effect should be easier to see.
#  - theta is a projection onto a reference, not a measurement. In sparse
#    single-cell/nucleus data it is noisy, and noise in a predictor attenuates
#    the association towards null.
#  - n = 1 animal per condition, so these are within-sample cell-level
#    associations only. Do not compare amplitudes between conditions.
# ===========================================================================

# NOTE ON MASKING: SingleCellExperiment pulls in matrixStats/S4Vectors/IRanges
# and tricycle pulls in AnnotationDbi. Between them they mask dplyr's count(),
# select(), filter(), slice(), rename(), first() and last(). Because these are
# attached after 00_functions.R has loaded dplyr, they win. Every dplyr verb
# below is namespace-qualified for that reason - do not "tidy" the dplyr::
# prefixes away, matrixStats::count() on a data frame fails with
# "Argument 'x' is not a vector: list".
library(tricycle)
library(SingleCellExperiment)

# tricycle expects library-size-normalised log counts, not SCT residuals, so
# build the SCE from the RNA assay. JoinLayers is a no-op on Seurat v4 objects
# and on v5 objects whose layers are already joined.
subset_heart_flt <- NormalizeData(subset_heart_flt, assay = "RNA", verbose = FALSE)
cc_sce <- SingleCellExperiment(
  assays = list(logcounts = GetAssayData(subset_heart_flt, assay = "RNA", layer = "data"))
)

cc_sce <- project_cycle_space(cc_sce, species = "mouse", gname.type = "SYMBOL")
cc_sce <- estimate_cycle_position(cc_sce)
cc_sce <- estimate_Schwabe_stage(cc_sce, species = "mouse", gname.type = "SYMBOL")

stopifnot(identical(colnames(cc_sce), colnames(subset_heart_flt)))
subset_heart_flt$tricyclePosition <- cc_sce$tricyclePosition
subset_heart_flt$CCStage          <- cc_sce$CCStage

cc_df <- data.frame(
  theta         = subset_heart_flt$tricyclePosition,
  CCStage       = subset_heart_flt$CCStage,
  allelic_ratio = subset_heart_flt$allelic_ratio,
  A1_reads      = subset_heart_flt$A1_reads,
  A2_reads      = subset_heart_flt$A2_reads,
  total_reads   = subset_heart_flt$total_reads,
  sample        = subset_heart_flt$sample,
  celltype      = subset_heart_flt$celltype,
  UMAP_1        = Embeddings(subset_heart_flt, reduction = "umap")[, 1],
  UMAP_2        = Embeddings(subset_heart_flt, reduction = "umap")[, 2]
) %>% dplyr::filter(!is.na(theta))

write.table(cc_df, file.path(CUTOFF_DIR, 'whole_chrX_tricycle_cell_positions.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)

# How much cell cycle variation is there at all? If a cell type is entirely
# G1/G0 the model below has nothing to detect, and that is a property of the
# tissue rather than a result.
cc_stage_counts <- cc_df %>%
  dplyr::count(sample, celltype, CCStage, name = "n_cells") %>%
  group_by(sample, celltype) %>%
  mutate(pct_of_celltype = round(100 * n_cells / sum(n_cells), 1)) %>%
  ungroup() %>%
  arrange(sample, celltype, CCStage)

write.table(cc_stage_counts, file.path(CUTOFF_DIR, 'whole_chrX_tricycle_stage_counts.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)

# --- Beta-binomial model of AR against cell cycle position -------------------
# theta is circular, so it cannot enter as a linear term - 0 and 2*pi are the
# same point. sin(theta) + cos(theta) is the first circular harmonic: it fits
# one peak and one trough around the cycle. The LRT against an intercept-only
# model is a 2 df test of "does AR depend on cell cycle position at all".
# Models are kept so the stats table and the fitted curves below come from the
# same fit rather than refitting for plotting.
cc_split <- split(cc_df, list(cc_df$celltype, cc_df$sample), drop = TRUE)
cc_split <- cc_split[sapply(cc_split, nrow) >= 30]

cc_models <- lapply(cc_split, function(x) {
  m0 <- try(glmmTMB(cbind(A1_reads, A2_reads) ~ 1,
                    family = betabinomial, data = x), silent = TRUE)
  m1 <- try(glmmTMB(cbind(A1_reads, A2_reads) ~ sin(theta) + cos(theta),
                    family = betabinomial, data = x), silent = TRUE)
  if (inherits(m0, "try-error") || inherits(m1, "try-error")) return(NULL)
  an <- try(anova(m0, m1), silent = TRUE)
  if (inherits(an, "try-error")) return(NULL)
  list(m1 = m1, p_value = data.frame(an)[2, "Pr..Chisq."])
})
keep_cc   <- !sapply(cc_models, is.null)
cc_split  <- cc_split[keep_cc]
cc_models <- cc_models[keep_cc]

cc_ar_bb <- base::Map(function(mm, x) {
  cf <- fixef(mm$m1)$cond
  b_sin <- unname(cf["sin(theta)"])
  b_cos <- unname(cf["cos(theta)"])
  data.frame(
    sample = unique(x$sample), celltype = unique(x$celltype), n_cells = nrow(x),
    # Amplitude of the AR swing around the cycle, on the log-odds scale.
    amplitude = sqrt(b_sin^2 + b_cos^2),
    # Position on the cycle where AR is highest, in radians [0, 2*pi).
    peak_theta = atan2(b_sin, b_cos) %% (2 * pi),
    p_value = mm$p_value
  )
}, cc_models, cc_split) %>% bind_rows()
rownames(cc_ar_bb) <- NULL
if (nrow(cc_ar_bb)) {
  cc_ar_bb$FDR <- p.adjust(cc_ar_bb$p_value, method = "fdr")
  # Amplitude as a fold change in odds between the peak and trough of the cycle.
  cc_ar_bb$peak_trough_OR <- exp(2 * cc_ar_bb$amplitude)
  cc_ar_bb$peak_phase <- cut(cc_ar_bb$peak_theta,
                             breaks = c(0, 0.5 * pi, pi, 1.5 * pi, 2 * pi),
                             labels = c("G1/S", "S", "G2M", "M/G1"),
                             include.lowest = TRUE)
  # The replication-timing prediction is specific about WHERE the peak falls,
  # not just that one exists. A significant harmonic peaking outside S is not
  # support for the hypothesis and needs a different explanation.
  cc_ar_bb$peak_in_S_window <- cc_ar_bb$peak_theta >= 0.5 * pi &
                               cc_ar_bb$peak_theta <= pi
}

write.table(cc_ar_bb, file.path(CUTOFF_DIR, 'whole_chrX_tricycle_AR_betabinomial.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)

# --- Pre-specified contrast: S phase vs G1/G0 --------------------------------
# The omnibus harmonic test above spends 2 df detecting a peak anywhere on the
# cycle. The replication-timing hypothesis predicts one specific place, so a
# targeted 1 df contrast has considerably more power.
#
# Lumping S with G2M would dilute it: by G2 the Xi has replicated too, dosage
# is back to 2:2, and AR should have returned to baseline. G2M is therefore
# kept as its own level and acts as an internal negative control - the
# hypothesis predicts S > G1/G0 with G2M indistinguishable from G1/G0.
# Windows are defined on theta using the mapping documented by tricycle
# (0.5*pi start of S, pi start of G2M, 1.5*pi middle of M, 1.75*pi-0.25*pi
# G1/G0) rather than on the discrete CCStage labels. That keeps the contrast
# independent of tricycle's stage-naming and matches the shaded band in the
# figure exactly. Cells between windows are left NA and excluded.
cc_df <- cc_df %>%
  mutate(cc_phase = factor(dplyr::case_when(
    theta >= 1.75 * pi | theta <= 0.25 * pi ~ "G1/G0",
    theta >= 0.5 * pi  & theta <= pi        ~ "S",
    theta >  pi        & theta <= 1.5 * pi  ~ "G2M",
    TRUE                                    ~ NA_character_
  ), levels = c("G1/G0", "S", "G2M")))

stopifnot(sum(!is.na(cc_df$cc_phase)) > 0)

cc_phase_bb <- split(cc_df, list(cc_df$celltype, cc_df$sample), drop = TRUE) %>%
  lapply(function(x) {
    x <- dplyr::filter(x, !is.na(cc_phase))
    if (nrow(x) < 30) return(NULL)
    tb <- table(droplevels(x$cc_phase))
    if (!all(c("G1/G0", "S") %in% names(tb)) || min(tb[c("G1/G0", "S")]) < 5) return(NULL)
    x$cc_phase <- factor(x$cc_phase, levels = c("G1/G0", "S", "G2M"))
    m <- try(glmmTMB(cbind(A1_reads, A2_reads) ~ cc_phase,
                     family = betabinomial, data = x), silent = TRUE)
    if (inherits(m, "try-error")) return(NULL)
    s <- summary(m)$coefficients$cond
    get_row <- function(nm) if (nm %in% rownames(s)) s[nm, ] else c(NA, NA, NA, NA)
    s_row   <- get_row("cc_phaseS")
    g2m_row <- get_row("cc_phaseG2M")
    data.frame(sample = unique(x$sample), celltype = unique(x$celltype),
               n_G1 = unname(tb["G1/G0"]), n_S = unname(tb["S"]),
               n_G2M = if ("G2M" %in% names(tb)) unname(tb["G2M"]) else 0L,
               # Predicted > 1 if the Xa:Xi template ratio is transiently 2:1.
               # Useful benchmark: a clean doubling of Xa template with the Xi
               # unreplicated gives OR = 2 exactly, whatever the baseline AR
               # (the odds Xa/Xi simply double). RNA buffering and the fact
               # that only part of the S window carries the imbalance should
               # pull the observed value well below 2, so expect 1 < OR < 2.
               # An OR at or above 2 implies something beyond gene dosage.
               OR_S_vs_G1 = exp(s_row[1]), p_S_vs_G1 = s_row[4],
               # Internal negative control: predicted ~ 1.
               OR_G2M_vs_G1 = exp(g2m_row[1]), p_G2M_vs_G1 = g2m_row[4])
  }) %>% bind_rows()
rownames(cc_phase_bb) <- NULL
if (nrow(cc_phase_bb)) {
  cc_phase_bb$FDR_S_vs_G1 <- p.adjust(cc_phase_bb$p_S_vs_G1, method = "fdr")
}

write.table(cc_phase_bb, file.path(CUTOFF_DIR, 'whole_chrX_tricycle_AR_S_vs_G1.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)

# --- Plots -------------------------------------------------------------------
# Cyclic hue palette, so theta = 0 and theta = 2*pi get the same colour and the
# scale does not imply a start and an end where the biology has neither.
theta_cols <- hcl(h = seq(15, 375, length.out = 100), c = 80, l = 60)

pdf(file.path(CUTOFF_DIR, 'whole_chrX_tricycle_umap.pdf'), width = 12, height = 3.5)
ggplot(cc_df, aes(x = UMAP_1, y = UMAP_2, colour = theta)) +
  geom_point(size = 0.6, alpha = 0.9) +
  facet_wrap(~sample, nrow = 1) +
  scale_colour_gradientn(colours = theta_cols, limits = c(0, 2 * pi),
                         breaks = c(0, 0.5, 1, 1.5, 2) * pi,
                         labels = c("0", "0.5pi\nS", "pi\nG2M", "1.5pi\nM", "2pi"),
                         name = "Cell cycle\nposition") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw() +
  theme(panel.grid = element_blank())
dev.off()

# Pooled AR within theta bins, with Wilson intervals, plus the fitted first
# harmonic. Same display logic as the Xist figures in 04.
wilson_ci_cc <- function(x, n, z = 1.96) {
  p <- x / n
  d <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / d
  halfw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  data.frame(est = p, lwr = pmax(0, centre - halfw), upr = pmin(1, centre + halfw))
}

theta_breaks <- seq(0, 2 * pi, length.out = 9)
cc_bins <- cc_df %>%
  mutate(theta_bin = cut(theta, breaks = theta_breaks, include.lowest = TRUE)) %>%
  group_by(sample, celltype, theta_bin) %>%
  summarise(theta = mean(theta), A1 = sum(A1_reads), N = sum(total_reads),
            n_cells = n(), .groups = "drop") %>%
  dplyr::filter(N > 0, n_cells >= 5)
cc_bins <- bind_cols(cc_bins, wilson_ci_cc(cc_bins$A1, cc_bins$N))

# Fitted curve from the harmonic model, reusing the fits from above.
cc_sig <- dplyr::select(cc_ar_bb, sample, celltype, FDR)

cc_pred <- base::Map(function(mm, x) {
  nd <- data.frame(theta = seq(0, 2 * pi, length.out = 200), A1_reads = 1, A2_reads = 1)
  p <- predict(mm$m1, newdata = nd, type = "link", se.fit = TRUE)
  data.frame(sample = unique(x$sample), celltype = unique(x$celltype),
             theta = nd$theta,
             fit = plogis(p$fit),
             lwr = plogis(p$fit - 1.96 * p$se.fit),
             upr = plogis(p$fit + 1.96 * p$se.fit))
}, cc_models, cc_split) %>%
  bind_rows() %>%
  left_join(cc_sig, by = c("sample", "celltype")) %>%
  mutate(sig = FDR < 0.05)

pdf(file.path(CUTOFF_DIR, 'whole_chrX_tricycle_AR_vs_position.pdf'), width = 11, height = 13)
ggplot(cc_df, aes(x = theta, y = allelic_ratio)) +
  # Shaded band is the predicted window: S phase, where the Xa has replicated
  # and the Xi has not. Drawn first so it sits behind the data.
  annotate("rect", xmin = 0.5 * pi, xmax = pi, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.55) +
  geom_point(aes(size = total_reads, colour = theta), alpha = 0.3) +
  geom_ribbon(data = cc_pred, aes(y = fit, ymin = lwr, ymax = upr, fill = sig),
              alpha = 0.2, colour = NA) +
  geom_line(data = cc_pred, aes(y = fit, linetype = sig), linewidth = 0.7) +
  geom_errorbar(data = cc_bins, aes(y = est, ymin = lwr, ymax = upr),
                width = 0, colour = "black", linewidth = 0.4) +
  geom_point(data = cc_bins, aes(y = est), size = 1.5, colour = "black") +
  facet_grid(celltype ~ sample, labeller = labeller(celltype = label_wrap_gen(18))) +
  scale_colour_gradientn(colours = theta_cols, limits = c(0, 2 * pi), guide = "none") +
  scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey50"), guide = "none") +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), guide = "none") +
  scale_size_continuous(range = c(0.3, 2.5), name = "Informative\nreads") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2) * pi,
                     labels = c("0", "0.5pi\nS", "pi\nG2M", "1.5pi\nM", "2pi")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Cell cycle position (tricycle)", y = "Allelic ratio (chrX)",
       caption = paste0("Curve: beta-binomial first harmonic, solid = FDR < 0.05. Black points: pooled AR within theta bins, Wilson 95% CI.\n",
                        "Grey band is S phase, where the Xa is replicated and the late-replicating Xi is not - AR is predicted to rise there and return to baseline in G2.\n",
                        "Ventricular cardiomyocytes are post-mitotic - a flat line there reflects the absence of S-phase cells, not an absence of association.")) +
  theme_bw() +
  theme(strip.text.y = element_text(size = 7),
        plot.caption = element_text(size = 7, hjust = 0))
dev.off()

# Polar view of the same thing: pooled AR around the cycle.
# geom_line under coord_polar errors on zero-row facet panels ("replacement has
# 1 row, data has 0") rather than skipping them, and facet_grid always draws
# the full celltype x sample grid including combinations with no cells. So
# panel on a single dropped factor with facet_wrap instead, and keep only
# groups with enough bins to draw a line through.
cc_bins_polar <- cc_bins %>%
  group_by(sample, celltype) %>%
  dplyr::filter(dplyr::n() >= 3) %>%
  ungroup() %>%
  arrange(celltype, sample) %>%
  mutate(panel = factor(paste0(celltype, " - ", sample),
                        levels = unique(paste0(celltype, " - ", sample))))

pdf(file.path(CUTOFF_DIR, 'whole_chrX_tricycle_AR_polar.pdf'), width = 12, height = 12)
ggplot(cc_bins_polar, aes(x = theta, y = est, group = panel)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, linewidth = 0.4, colour = "grey30") +
  geom_point(aes(colour = theta, size = n_cells)) +
  geom_line(colour = "grey30", linewidth = 0.4) +
  facet_wrap(~panel, ncol = 4, drop = TRUE, labeller = label_wrap_gen(24)) +
  coord_polar(theta = "x", start = 0) +
  scale_x_continuous(limits = c(0, 2 * pi), breaks = c(0, 0.5, 1, 1.5) * pi,
                     labels = c("G1/G0", "S", "G2M", "M")) +
  scale_colour_gradientn(colours = theta_cols, limits = c(0, 2 * pi), guide = "none") +
  scale_size_continuous(range = c(0.8, 3), name = "Cells") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Pooled allelic ratio",
       caption = paste0("Dashed circle marks AR = 0.5 (balanced). Bins with < 5 cells omitted; ",
                        "groups with < 3 bins not shown.\nThe line does not close the circle - ",
                        "the gap at the top spans the wrap from the last bin back to the first.")) +
  theme_bw() +
  theme(strip.text = element_text(size = 7),
        plot.caption = element_text(size = 7, hjust = 0))
dev.off()

# Discrete cross-check: AR by Schwabe stage.
pdf(file.path(CUTOFF_DIR, 'whole_chrX_tricycle_AR_by_stage.pdf'), width = 11, height = 4)
cc_df %>%
  dplyr::filter(!is.na(CCStage)) %>%
  ggplot(aes(x = CCStage, y = allelic_ratio)) +
  geom_jitter(width = 0.25, size = 0.3, alpha = 0.25, colour = "grey40") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "steelblue", width = 0.55) +
  facet_wrap(~sample, nrow = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Schwabe stage (tricycle)", y = "Allelic ratio (chrX)",
       caption = "Descriptive only - cells within a sample, n = 1 animal per condition.") +
  theme_bw() +
  theme(plot.caption = element_text(size = 7, hjust = 0))
dev.off()

cc_ar_bb
cc_phase_bb

# Diagnostic that decides whether any of the above is interpretable: real cell
# cycle signal makes cells trace a ring in the tricycle embedding. A blob means
# the cells sit on the origin, theta is the angle of a near-zero-length vector,
# and every downstream theta is noise. This dataset gives a blob.
pdf(file.path(CUTOFF_DIR, 'whole_chrX_tricycle_embedding.pdf'), width = 6, height = 6)
plot(reducedDim(cc_sce, "tricycleEmbedding"),
     pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.3),
     xlab = "tricycle PC1", ylab = "tricycle PC2",
     main = "tricycle embedding (ring = cell cycle signal, blob = none)")
abline(h = 0, v = 0, col = "grey70", lty = 2)
dev.off()


# ===========================================================================
# Sequencing depth as a confounder of the allelic ratio.
#
# This came out of the cell cycle block above: AR tracked the Schwabe stage
# ordering, and that ordering tracks total RNA content rather than cell cycle.
# The underlying issue matters on its own.
#
# Under a beta-binomial, depth should affect the PRECISION of a cell's AR and
# not its expected value - more reads means a tighter estimate of the same p.
# So a fitted mean that moves with depth is evidence of a genuine bias
# mechanism, not sampling noise.
#
# The likely mechanism: baseline chrX AR is ~0.95, i.e. the Xi contributes
# only ~5% of chrX reads, because only escapees are transcribed from it. A
# roughly constant absolute amount of ambient RNA, dominated by the more
# highly expressed Xa, is then a large proportional contamination in a shallow
# cell and a small one in a deep cell. That inflates AR in shallow cells and
# would show up as a negative coefficient below.
#
# Why it matters: if cells, cell types or samples differ in depth, they differ
# in APPARENT escape for purely technical reasons. Any condition comparison of
# escape needs to be shown not to be depth-driven.
#
# Consistent with the AR ~ 0 excess in 04_core_escape.R, which showed a ~17x
# excess of impossible monoallelic calls over binomial expectation at a stable
# 6-8% rate across all four samples - the same artefact process seen from the
# other end.
# ===========================================================================

if (!"nCount_RNA" %in% colnames(subset_heart_flt@meta.data)) {
  stop("nCount_RNA not found in metadata. Substitute whichever library-size ",
       "column this object carries (e.g. nCount_Spatial, nCount_SCT is NOT a ",
       "substitute - it is post-normalisation).")
}

depth_df <- data.frame(
  allelic_ratio = subset_heart_flt$allelic_ratio,
  A1_reads      = subset_heart_flt$A1_reads,
  A2_reads      = subset_heart_flt$A2_reads,
  total_reads   = subset_heart_flt$total_reads,
  nCount_RNA    = subset_heart_flt$nCount_RNA,
  nFeature_RNA  = subset_heart_flt$nFeature_RNA,
  sample        = subset_heart_flt$sample,
  celltype      = subset_heart_flt$celltype
) %>%
  dplyr::filter(!is.na(nCount_RNA), nCount_RNA > 0) %>%
  mutate(log10_depth = log10(nCount_RNA))

# Library depth and chrX informative reads are not the same quantity, so test
# both. total_reads is the direct driver of the binomial noise; nCount_RNA is
# the library size that ambient contamination scales against.
fit_depth <- function(x, predictor) {
  if (nrow(x) < 30) return(NULL)
  f <- as.formula(paste("cbind(A1_reads, A2_reads) ~", predictor))
  m <- try(glmmTMB(f, family = betabinomial, data = x), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  s <- summary(m)$coefficients$cond
  if (!predictor %in% rownames(s)) return(NULL)
  data.frame(sample = unique(x$sample), celltype = unique(x$celltype),
             n_cells = nrow(x), predictor = predictor,
             beta = s[predictor, 1], se = s[predictor, 2],
             OR = exp(s[predictor, 1]), p_value = s[predictor, 4])
}

depth_df$log10_total_reads <- log10(depth_df$total_reads)
depth_split <- split(depth_df, list(depth_df$celltype, depth_df$sample), drop = TRUE)

depth_bb <- bind_rows(
  bind_rows(lapply(depth_split, fit_depth, predictor = "log10_depth")),
  bind_rows(lapply(depth_split, fit_depth, predictor = "log10_total_reads"))
)
rownames(depth_bb) <- NULL
if (nrow(depth_bb)) {
  depth_bb <- depth_bb %>%
    group_by(predictor) %>%
    mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
    ungroup() %>%
    arrange(predictor, sample, celltype)
}

write.table(depth_bb, file.path(CUTOFF_DIR, 'whole_chrX_AR_vs_depth_betabinomial.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)

# Direction summary: how many groups show AR falling with depth, and how many
# of those are significant. A consistent negative sign across many independent
# groups is the signal, not any single p-value.
depth_direction <- depth_bb %>%
  group_by(predictor) %>%
  summarise(n_groups = dplyr::n(),
            n_negative = sum(beta < 0),
            n_sig = sum(FDR < 0.05),
            n_sig_negative = sum(FDR < 0.05 & beta < 0),
            median_OR = median(OR),
            .groups = "drop")

write.table(depth_direction, file.path(CUTOFF_DIR, 'whole_chrX_AR_vs_depth_direction_summary.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)

# Forest plot of the depth coefficients.
pdf(file.path(CUTOFF_DIR, 'whole_chrX_AR_vs_depth_forest.pdf'), width = 10, height = 6)
depth_bb %>%
  mutate(sig = FDR < 0.05,
         lwr = exp(beta - 1.96 * se),
         upr = exp(beta + 1.96 * se),
         celltype = factor(celltype, levels = rev(sort(unique(celltype))))) %>%
  ggplot(aes(x = OR, y = celltype, colour = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.5) +
  geom_point(aes(size = n_cells)) +
  facet_grid(predictor ~ sample) +
  scale_x_log10() +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40"),
                      name = NULL, labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "ns")) +
  scale_size_continuous(range = c(1, 3.5), name = "Cells") +
  labs(x = "Odds ratio per 10-fold increase in depth (95% CI)", y = NULL,
       caption = paste0("OR < 1: deeper cells show LOWER allelic ratio, i.e. more apparent escape. ",
                        "Under a beta-binomial, depth should not move the mean at all,\nso a consistent ",
                        "departure from 1 indicates a bias mechanism (ambient RNA is the leading ",
                        "candidate) rather than sampling noise.")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7),
        plot.caption = element_text(size = 7, hjust = 0))
dev.off()

# Binned view: pooled AR against library size, which shows the shape of the
# relationship rather than just its slope.
# Local copy rather than reusing wilson_ci_cc from the cell cycle block above,
# so this section still runs if that block is ever removed.
wilson_ci_depth <- function(x, n, z = 1.96) {
  p <- x / n
  d <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / d
  halfw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  data.frame(est = p, lwr = pmax(0, centre - halfw), upr = pmin(1, centre + halfw))
}

depth_bins <- depth_df %>%
  group_by(sample, celltype) %>%
  mutate(depth_bin = dplyr::ntile(nCount_RNA, 6)) %>%
  group_by(sample, celltype, depth_bin) %>%
  summarise(depth = median(nCount_RNA), A1 = sum(A1_reads), N = sum(total_reads),
            n_cells = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(N > 0, n_cells >= 10)
depth_bins <- bind_cols(depth_bins, wilson_ci_depth(depth_bins$A1, depth_bins$N))

pdf(file.path(CUTOFF_DIR, 'whole_chrX_AR_vs_depth_binned.pdf'), width = 11, height = 13)
ggplot(depth_bins, aes(x = depth, y = est)) +
  geom_line(aes(group = interaction(sample, celltype)), colour = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, linewidth = 0.4) +
  geom_point(aes(size = n_cells), colour = "steelblue") +
  facet_grid(celltype ~ sample, labeller = labeller(celltype = label_wrap_gen(18))) +
  scale_x_log10() +
  scale_size_continuous(range = c(0.8, 3), name = "Cells") +
  labs(x = "Library size (nCount_RNA, log scale)", y = "Pooled allelic ratio (chrX)",
       caption = "Pooled counts within library-size sextiles, Wilson 95% CI. A downward slope means shallow cells look more monoallelic.") +
  theme_bw() +
  theme(strip.text.y = element_text(size = 7),
        plot.caption = element_text(size = 7, hjust = 0))
dev.off()

depth_direction
depth_bb