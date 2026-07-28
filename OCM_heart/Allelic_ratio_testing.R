library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(broom.mixed)


run_dispersion_lrt <- function(df, ref_level, comparison_level) {
  df$sample <- factor(df$sample, levels = c(ref_level, comparison_level))

  fit_const_disp <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                            family = betabinomial(link = "logit"), data = df)
  fit_var_disp <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                          dispformula = ~ sample,
                          family = betabinomial(link = "logit"), data = df)

  lrt_row <- data.frame(anova(fit_const_disp, fit_var_disp))[2, , drop = FALSE]
  disp_coef_name <- paste0("sample", comparison_level)
  disp_beta <- unname(fixef(fit_var_disp)$disp[disp_coef_name])

  direction <- if (is.na(disp_beta)) {
    NA_character_
  } else if (disp_beta < 0) {
    paste(comparison_level, "more heterogeneous")
  } else {
    paste(comparison_level, "less heterogeneous")
  }

  lrt_row$disp_beta <- disp_beta
  lrt_row$direction <- direction
  lrt_row
}

# Simulation-based false-positive rate for run_dispersion_lrt, quantifying the
# n=1-animal-per-condition caveat: if the two conditions truly share the SAME
# dispersion, and each condition is still just a single animal, how often does
# the LRT call a "significant" difference purely from animal-to-animal noise?
#
# Mechanics: simulate two animals from IDENTICAL true (mu, phi), but let each
# animal's *realized* phi wobble around that shared truth on the log scale
# (animal_disp_sd) -- this represents ordinary biological variability between
# individuals, which a single-animal-per-condition design can never separate
# from a genuine condition effect. Cells are simulated at the real per-cell
# read depths and real group sizes so the simulation matches the actual data
# structure. animal_disp_sd cannot be estimated from n=1 data -- it has to be
# assumed, so report a range of plausible values rather than trusting one.
simulate_dispersion_null_fpr <- function(df, ref_level, comparison_level,
                                         animal_disp_sd = 0.3, n_reps = 1000,
                                         seed = 1) {
  set.seed(seed)
  df$sample <- factor(df$sample, levels = c(ref_level, comparison_level))

  # anchor simulation parameters to this celltype/comparison's real data:
  # real per-group mean allelic proportion and the real shared dispersion
  # under the null (single phi fit across both groups)
  fit_null <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                      family = betabinomial(link = "logit"), data = df)
  mu_hat <- predict(fit_null,
                     newdata = data.frame(sample = factor(c(ref_level, comparison_level),
                                                           levels = c(ref_level, comparison_level))),
                     type = "response")
  mu_ref  <- unname(mu_hat[1])
  mu_comp <- unname(mu_hat[2])
  phi_true <- exp(unname(fixef(fit_null)$disp[["(Intercept)"]]))

  ref_reads  <- df$total_reads[df$sample == ref_level]
  comp_reads <- df$total_reads[df$sample == comparison_level]

  p_values <- rep(NA_real_, n_reps)
  n_failed <- 0

  for (i in seq_len(n_reps)) {
    # each simulated animal gets its own realized dispersion, both drawn from
    # the same distribution since the true dispersion is identical by design
    phi_ref  <- phi_true * exp(rnorm(1, 0, animal_disp_sd))
    phi_comp <- phi_true * exp(rnorm(1, 0, animal_disp_sd))

    p_ref  <- rbeta(length(ref_reads),  mu_ref  * phi_ref,  (1 - mu_ref)  * phi_ref)
    p_comp <- rbeta(length(comp_reads), mu_comp * phi_comp, (1 - mu_comp) * phi_comp)

    A1_ref  <- rbinom(length(ref_reads),  ref_reads,  p_ref)
    A1_comp <- rbinom(length(comp_reads), comp_reads, p_comp)

    df_sim <- data.frame(
      A1_reads = c(A1_ref, A1_comp),
      A2_reads = c(ref_reads - A1_ref, comp_reads - A1_comp),
      sample = c(rep(ref_level, length(ref_reads)), rep(comparison_level, length(comp_reads)))
    )

    lrt_row <- tryCatch(run_dispersion_lrt(df_sim, ref_level, comparison_level),
                         error = function(e) NULL, warning = function(w) NULL)

    if (is.null(lrt_row)) {
      n_failed <- n_failed + 1
    } else {
      p_values[i] <- lrt_row$Pr..Chisq.
    }
  }

  list(
    fpr = mean(p_values < 0.05, na.rm = TRUE),
    n_failed = n_failed,
    n_reps = n_reps,
    p_values = p_values
  )
}

# Logistic regression LRT testing whether the proportion of monoallelic-like
# cells (allelic_ratio <= 0.1 or >= 0.9) differs between two sample levels.
run_monoallelic_lrt <- function(df, ref_level, comparison_level) {
  df$sample <- factor(df$sample, levels = c(ref_level, comparison_level))

  fit_null <- glm(monoallelic ~ 1, family = binomial, data = df)
  fit_full <- glm(monoallelic ~ sample, family = binomial, data = df)

  lrt <- anova(fit_null, fit_full, test = "LRT")

  coef_name <- paste0("sample", comparison_level)
  beta <- unname(coef(fit_full)[coef_name])

  direction <- if (is.na(beta)) {
    NA_character_
  } else if (beta > 0) {
    paste(comparison_level, "more monoallelic")
  } else {
    paste(comparison_level, "less monoallelic")
  }

  data.frame(
    Df = lrt$Df[2],
    Deviance = lrt$Deviance[2],
    p_value = lrt[["Pr(>Chi)"]][2],
    beta = beta,
    direction = direction
  )
}

fdr_to_stars <- function(fdr) {
  ifelse(is.na(fdr), "ns",
         ifelse(fdr <= 0.001, "***",
                ifelse(fdr <= 0.01, "**",
        ifelse(fdr <= 0.05, "*", "ns"))))
}

############################
# Read in single cell data #
############################
heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

# Subcluster the Ventricular Cardiomyocytes to check for the two suspected
# subpopulations. Kept as a separate metadata column (celltype_sub) rather
# than overwriting `celltype`, since downstream code filters on the exact
# string "Ventricular Cardiomyocytes" and shouldn't silently lose those cells.
vcm <- subset(heart, subset = celltype == "Ventricular Cardiomyocytes")

DefaultAssay(vcm) <- "RNA"
vcm <- SCTransform(vcm, verbose = FALSE)
vcm <- RunPCA(vcm, features = VariableFeatures(vcm), npcs = 30)

pdf("Allelic_ratio_results/VCM_subcluster_ElbowPlot.pdf")
ElbowPlot(vcm, ndims = 30)
dev.off()

# NB: check the elbow plot output and adjust `dims` below if needed before
# trusting the clustering
vcm <- RunUMAP(vcm, dims = 1:15)
vcm <- FindNeighbors(vcm, dims = 1:15)
vcm <- FindClusters(vcm, resolution = 0.05)

pdf("Allelic_ratio_results/VCM_subcluster_UMAP.pdf")
DimPlot(vcm, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

# Check whether the split tracks sample/condition rather than a real
# biological subpopulation (e.g. batch effect)
pdf("Allelic_ratio_results/VCM_subcluster_UMAP_by_sample.pdf")
DimPlot(vcm, reduction = "umap", group.by = "sample", raster = FALSE)
dev.off()

# Markers distinguishing the subclusters
vcm <- PrepSCTFindMarkers(vcm)
vcm.markers <- FindAllMarkers(vcm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.table(vcm.markers, file = "Allelic_ratio_results/VCM_subcluster_markers.txt")

pdf("Allelic_ratio_results/VCM_subcluster_markers_dotplot.pdf", width = 10, height = 8)
DotPlot(vcm, features = unique(unlist(lapply(split(vcm.markers, vcm.markers$cluster), function(x) {
  x[order(x$avg_log2FC, decreasing = TRUE), "gene"][1:10]
})))) + RotatedAxis() +
theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8))
dev.off()

# depth + complexity per subcluster — is this a technical split?
pdf("Allelic_ratio_results/VCM_subcluster_QC_violinplots.pdf", width = 10, height = 8)
VlnPlot(vcm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "seurat_clusters", pt.size = 0)
dev.off()
table(vcm$seurat_clusters, vcm$sample)



# Push subcluster labels back onto the full heart object as a new column
heart$celltype_sub <- as.character(heart$celltype)
match_idx <- match(Cells(vcm), colnames(heart))
heart$celltype_sub[match_idx] <- paste0("Ventricular Cardiomyocytes_", Idents(vcm))
heart$celltype_sub <- factor(heart$celltype_sub)

heart_no_cluster1 <- subset(heart, subset = celltype_sub != "Ventricular Cardiomyocytes_1")

# Plot UMAP with subcluster labels
pdf("Allelic_ratio_results/heart_UMAP_with_VCM_subclusters.pdf")
DimPlot(heart_no_cluster1, reduction = "umap", group.by = "celltype_sub", label = TRUE, raster = FALSE) +
  theme(legend.position = "none")
dev.off()

saveRDS(heart, file = "heart_seurat_object_SCT.rds")

heart_metadata <- heart@meta.data
# box plot of nCount_RNA per sample, split by celltype
pdf('Allelic_ratio_results/nCount_RNA_boxplot_by_celltype.pdf')
ggplot(heart_metadata, aes(x = celltype, y = nCount_RNA, fill = sample)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~sample) +
  labs(y = "nCount_RNA", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.off()

# Boxplot of Xist expression (approximate TPM) per sample, split by celltype
# Xist gene length from mm10 chrX:102,503,972-102,537,573
xist_length_kb <- (102537573 - 102503972 + 1) / 1000

xist_counts <- GetAssayData(heart, assay = "RNA", layer = "counts")["Xist", ]
total_counts <- heart$nCount_RNA

xist_df <- data.frame(
  sample = heart$sample,
  celltype = heart$celltype,
  Xist_TPM = (xist_counts / xist_length_kb) / (total_counts / 1e6)
)

pdf('Allelic_ratio_results/xist_expression_TPM_boxplot_by_celltype.pdf')
ggplot(xist_df, aes(x = sample, y = Xist_TPM, fill = sample)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~celltype) +
  labs(y = "Xist TPM (approx.)", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.off()

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

#####################
# All gene analysis #
#####################
barcodes <- lapply(c('9w', '78w', 'Sham', 'TAC'), function(x) {
  list.files(paste0('Allelome.PRO2_all_genes/', x), pattern = 'locus_table.txt', recursive = TRUE, full.names = TRUE)
})
barcodes <- unlist(barcodes)

condition <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[6])
cellid <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[7])

read_gene_table <- function(path, sample_name, cell_name) {
  if (!file.exists(path)) return(NULL)

  tmp <- tryCatch(
    read.delim(path, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if (is.null(tmp)) return(NULL)

  tmp <- as.data.frame(tmp)
  if (!"chr" %in% names(tmp) || !"name" %in% names(tmp) || nrow(tmp) == 0) return(NULL)

  tmp$chr <- as.character(tmp$chr)
  tmp$name <- as.character(tmp$name)

  tmp %>%
    dplyr::mutate(
      sample = sample_name,
      cellid = cell_name,
      cell_barcode = paste0(sample_name, "_", cell_name)
    )
}

gene_df <- lapply(seq_along(barcodes), function(i) {
  read_gene_table(barcodes[i], condition[i], cellid[i])
}) %>%
  bind_rows()

# Match cell types from Seurat object
gene_df$celltype <- Idents(heart)[gene_df$cell_barcode]

# write out the all genes allelic ratio table
write.table(gene_df, 'Allelic_ratio_results/all_genes_gene_df.txt', sep = '\t', row.names = FALSE, quote = FALSE)

################################################################
# Xist allelic purity as a direct ambient-contamination probe  #
################################################################
# Xist is transcribed ONLY from the inactive X, so in a cell with intact XCI
# its reads should be ~100% one allele -- and the OPPOSITE allele to the rest
# of chrX, which is expressed from the active X. Any minor-allele Xist signal
# is therefore ambient soup + mapping/sequencing error, not biology.
#
# Because the soup pools cells with random XCI, its Xist is ~50:50. A cell with
# contamination fraction c shows a Xist minor-allele fraction of ~c/2, so
# ambient ~= 2 * minor fraction. This is an UPPER bound (mapping/seq error is
# folded in), but the error floor is shared across cells, so the comparison
# BETWEEN VCM subclusters is still clean.
#
# The point of doing it this way: it is measured on the same Allelome.PRO2
# read-counting pipeline as everything else, unlike an expression-based
# estimate (decontX), which cannot see A1_reads/A2_reads at all.
xist_allelic <- gene_df %>%
  filter(name == "Xist", total_reads >= 10) %>%
  mutate(
    xist_ratio      = A1_reads / (A1_reads + A2_reads),
    xist_minor_frac = pmin(xist_ratio, 1 - xist_ratio),
    ambient_est     = pmin(2 * xist_minor_frac, 1)
  )

if (nrow(xist_allelic) == 0) {
  warning("No Xist rows with total_reads >= 10 -- check Xist is in the Allelome.PRO2 annotation")
} else {

xist_allelic$celltype_sub <- heart$celltype_sub[xist_allelic$cell_barcode]
xist_allelic$chrX_ratio <- metadata_whole_chr$allelic_ratio[
  match(xist_allelic$cell_barcode, rownames(metadata_whole_chr))]

# 1. Sanity check: Xist should be near-monoallelic (bimodal at 0 and 1).
#    A unimodal pile at 0.5 would mean the allele assignment is wrong or
#    contamination is severe enough to swamp the signal.
pdf('Allelic_ratio_results/Xist_allelic_ratio_distribution.pdf', width = 10, height = 4)
ggplot(xist_allelic, aes(x = xist_ratio)) +
  geom_histogram(bins = 50) +
  facet_wrap(~sample) +
  labs(x = "Xist allelic ratio (A1 / total)", y = "cells") +
  theme_bw()
dev.off()

# 2. Orientation check: Xist comes from the inactive X, the rest of chrX from
#    the active X, so these should be strongly ANTI-correlated. If they are
#    positively correlated, A1/A2 are not consistently oriented between the
#    per-gene and whole-chromosome tables.
pdf('Allelic_ratio_results/Xist_vs_chrX_allelic_ratio.pdf', width = 6, height = 5)
ggplot(xist_allelic, aes(x = chrX_ratio, y = xist_ratio)) +
  geom_point(size = 0.4, alpha = 0.3) +
  labs(x = "whole-chrX allelic ratio (active X)",
       y = "Xist allelic ratio (inactive X)") +
  theme_bw()
dev.off()

cor.test(xist_allelic$chrX_ratio, xist_allelic$xist_ratio, method = "spearman")

# 3. THE TEST: does the ambient estimate differ between the VCM subclusters?
vcm_xist <- xist_allelic %>%
  filter(grepl("^Ventricular Cardiomyocytes", celltype_sub)) %>%
  mutate(celltype_sub = droplevels(factor(celltype_sub)))

vcm_xist %>%
  group_by(celltype_sub) %>%
  summarise(n = n(),
            median_ambient = median(ambient_est, na.rm = TRUE),
            mean_ambient   = mean(ambient_est, na.rm = TRUE),
            median_reads   = median(total_reads, na.rm = TRUE),
            .groups = "drop")

pdf('Allelic_ratio_results/Xist_ambient_estimate_by_VCM_subcluster.pdf', width = 6, height = 5)
ggplot(vcm_xist, aes(x = celltype_sub, y = ambient_est)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(y = "estimated ambient fraction (2 x Xist minor allele)", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

wilcox.test(ambient_est ~ celltype_sub, data = vcm_xist)

}

######################################################################
# Per-gene discriminator: ambient contamination vs. genuine escape   #
######################################################################
# The pooled minor-allele fraction confounds three things: escape, ambient
# RNA, and mapping/sequencing error. This separates the first two by their
# SHAPE across genes rather than their magnitude:
#
#   ambient RNA  is a per-CELL property -> shifts EVERY chrX gene equally
#   escape       is a per-GENE property -> shifts a HANDFUL of loci
#
# So if subcluster 0 differs from subcluster 1 by a near-uniform offset across
# essentially all chrX genes, that is contamination. If the difference sits in
# a few genes with the rest on the diagonal, that is real differential escape
# and those genes are the result.
#
# Allele orientation: the active-X allele is whichever is in the majority over
# the whole chromosome for that cell. Majority orientation is robust to
# contamination -- a truly monoallelic cell with contamination c still reads
# 1 - c/2, which stays > 0.5 for any c < 1.

# per-cell orientation, taken from the whole-chrX ratio
cell_orient <- data.frame(
  cell_barcode = rownames(metadata_whole_chr),
  celltype_sub = metadata_whole_chr$celltype_sub,
  chrX_ratio   = metadata_whole_chr$allelic_ratio,
  stringsAsFactors = FALSE
)
cell_orient <- cell_orient[!is.na(cell_orient$chrX_ratio), ]
cell_orient$active_allele <- ifelse(cell_orient$chrX_ratio > 0.5, "A1", "A2")

# chrX gene-level reads for VCM cells only, oriented against the active X
gene_x <- gene_df %>%
  filter(chr == "chrX", !is.na(A1_reads), !is.na(A2_reads),
         (A1_reads + A2_reads) > 0) %>%
  inner_join(cell_orient, by = "cell_barcode") %>%
  filter(grepl("^Ventricular Cardiomyocytes", celltype_sub)) %>%
  mutate(
    sub            = sub("^Ventricular Cardiomyocytes_", "c", celltype_sub),
    inactive_reads = ifelse(active_allele == "A1", A2_reads, A1_reads),
    gene_total     = A1_reads + A2_reads
  )

# pool over cells WITHIN a subcluster, after orienting each cell individually.
# (Pooling without orienting first would wash everything to 0.5, since which
# allele is active varies randomly from cell to cell.)
per_gene <- gene_x %>%
  group_by(name, sub) %>%
  summarise(inactive = sum(inactive_reads),
            total    = sum(gene_total),
            n_cells  = n(),
            .groups  = "drop") %>%
  filter(total >= 50) %>%
  mutate(inactive_frac = inactive / total)

g0 <- per_gene %>% filter(sub == "c0") %>%
  select(name, frac0 = inactive_frac, tot0 = total, inact0 = inactive, n0 = n_cells)
g1 <- per_gene %>% filter(sub == "c1") %>%
  select(name, frac1 = inactive_frac, tot1 = total, inact1 = inactive, n1 = n_cells)

per_gene_wide <- inner_join(g0, g1, by = "name") %>%
  mutate(delta = frac0 - frac1)

# flag the known core escape genes so they can be seen against the background
core_escape_names <- character(0)
if (file.exists('Allelic_ratio_results/core_escape_genes_gene_df.txt')) {
  core_escape_names <- unique(read.table(
    'Allelic_ratio_results/core_escape_genes_gene_df.txt',
    sep = '\t', header = TRUE, stringsAsFactors = FALSE)$name)
}
per_gene_wide$is_core_escape <- per_gene_wide$name %in% core_escape_names

write.table(per_gene_wide,
            'Allelic_ratio_results/VCM_subcluster_per_gene_inactive_fraction.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

# ---- READOUT 1: the shape of the per-gene difference ----
# unimodal and clearly off zero -> ambient. centred on zero with tail -> escape.
pdf('Allelic_ratio_results/VCM_subcluster_per_gene_delta_histogram.pdf', width = 6, height = 5)
ggplot(per_gene_wide, aes(x = delta)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = 0, colour = "red", linetype = 2) +
  labs(x = "inactive-allele fraction: subcluster 0 - subcluster 1",
       y = "chrX genes") +
  theme_bw()
dev.off()

# ---- READOUT 2: same thing as a scatter ----
# ambient -> cloud sits parallel to but offset from y=x.
# escape  -> cloud sits ON y=x with a few outliers.
pdf('Allelic_ratio_results/VCM_subcluster_per_gene_scatter.pdf', width = 6, height = 5.5)
ggplot(per_gene_wide, aes(x = frac1, y = frac0)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = 2) +
  geom_point(aes(colour = is_core_escape), size = 0.9, alpha = 0.6) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                      name = "core escape gene") +
  labs(x = "inactive-allele fraction (subcluster 1)",
       y = "inactive-allele fraction (subcluster 0)") +
  theme_bw()
dev.off()

# ---- READOUT 3: the sign test, which is the quantitative version ----
# Under ambient, nearly every gene shifts the SAME way, so the proportion of
# genes with delta > 0 heads for 0 or 1. Under differential escape affecting
# only a few loci, it sits near 0.5.
n_pos <- sum(per_gene_wide$delta > 0, na.rm = TRUE)
n_tot <- sum(!is.na(per_gene_wide$delta))

cat("genes tested:            ", n_tot, "\n")
cat("median delta:            ", round(median(per_gene_wide$delta, na.rm = TRUE), 4), "\n")
cat("genes shifted up in c0:  ", n_pos, "/", n_tot,
    " (", round(100 * n_pos / n_tot, 1), "%)\n", sep = "")
print(binom.test(n_pos, n_tot, p = 0.5))

# ---- per-gene outlier test ----
# NB: reads are pooled across cells, so these p-values ignore cell-level
# overdispersion and are anticonservative. Treat them as a ranking to pick
# candidate genes out of the tail, NOT as calibrated significance -- the
# shape readouts above are the actual evidence.
per_gene_wide$p_raw <- mapply(function(i0, t0, i1, t1) {
  tryCatch(fisher.test(matrix(c(i0, t0 - i0, i1, t1 - i1), nrow = 2))$p.value,
           error = function(e) NA_real_)
}, per_gene_wide$inact0, per_gene_wide$tot0,
   per_gene_wide$inact1, per_gene_wide$tot1)
per_gene_wide$FDR <- p.adjust(per_gene_wide$p_raw, method = "BH")

per_gene_wide %>%
  arrange(FDR) %>%
  select(name, frac0, frac1, delta, tot0, tot1, is_core_escape, FDR) %>%
  head(20)

# Plot distribution of total reads across all genes
pdf('Allelic_ratio_results/all_genes_total_reads_distribution.pdf')
plot(
  density(gene_df$total_reads, na.rm = TRUE),
  main = "Total Reads Distribution",
  xlab = "Total Reads"
)
dev.off()

dim(gene_df)
dim(subset(gene_df, total_reads >= 10))

# per-cell analysis for genes with signal (cell resolution preserved)

gene_df_signal <- gene_df %>%
  filter(!is.na(A1_reads), !is.na(A2_reads), !is.na(total_reads), total_reads > 0)

if (nrow(gene_df_signal) > 0) {
  signal_gene_celltype <- gene_df_signal %>%
    group_by(name, celltype) %>%
    summarise(total_reads_sum = sum(total_reads, na.rm = TRUE), .groups = "drop") %>%
    filter(total_reads_sum >= 10)

  signal_genes <- signal_gene_celltype %>%
    group_by(name) %>%
    summarise(total_reads_sum = sum(total_reads_sum), .groups = "drop") %>%
    arrange(desc(total_reads_sum))

  write.table(signal_gene_celltype,
              'Allelic_ratio_results/all_genes_signal_genes_by_celltype.txt',
              sep = '\t', row.names = FALSE, quote = FALSE)

  # keep only the (gene, celltype) combinations that individually cleared MIN_READS,
  # pooled across samples -- other celltypes for the same gene are dropped, not just shown empty
  per_cell_signal_genes <- gene_df_signal %>%
    semi_join(signal_gene_celltype, by = c("name", "celltype"))

  per_cell_signal_genes$sample <- factor(per_cell_signal_genes$sample, levels = c("9w", "78w", "Sham", "TAC"))
  per_cell_signal_genes$name <- factor(per_cell_signal_genes$name, levels = signal_genes$name)

  write.table(per_cell_signal_genes,
              'Allelic_ratio_results/all_genes_per_cell_signal_genes.txt',
              sep = '\t', row.names = FALSE, quote = FALSE)

  pdf('Allelic_ratio_results/all_genes_signal_total_reads_barplot.pdf')
  ggplot(signal_genes, aes(x = reorder(name, total_reads_sum), y = total_reads_sum)) +
    geom_col(fill = '#3b7a57', width = 0.7) +
    coord_flip() +
    labs(x = NULL, y = 'Total reads (summed across cells in qualifying celltypes)', title = paste('Genes with >=', MIN_READS, 'reads in at least one celltype (pooled across samples)')) +
    theme_bw()
  dev.off()

  pdf('Allelic_ratio_results/all_genes_per_cell_allelic_ratio_by_gene.pdf')
  ggplot(per_cell_signal_genes, aes(x = celltype, y = allelic_ratio)) +
    geom_jitter(aes(size = total_reads, color = sample), width = 0.2, height = 0, alpha = 0.6) +
    facet_wrap(~name) +
    scale_size_continuous(range = c(0.5, 4)) +
    labs(x = NULL, y = 'A1 / total', size = 'Total reads', color = 'Sample') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()
}

# inspect the signal-gene ranking to judge whether MIN_READS is picking out
# a sensible number of "main player" genes, before committing to a single one
head(signal_genes, 20)
quantile(signal_genes$total_reads_sum, probs = c(0.5, 0.75, 0.9, 0.95, 0.99))

# stricter gene selection: require at least one cell with total_reads > 10 in
# every celltype x sample combination, restricted to the 4 largest celltypes
# (requiring every celltype in the object left zero genes -- too sparse)
top4_celltypes <- names(sort(table(Idents(heart)), decreasing = TRUE))[1:4]
n_combos <- length(top4_celltypes) * 4

gene_combo_coverage <- gene_df_signal %>%
  filter(total_reads > 10, celltype %in% top4_celltypes) %>%
  distinct(name, celltype, sample) %>%
  group_by(name) %>%
  summarise(n_combos_covered = n(), .groups = "drop")

full_coverage_genes <- gene_combo_coverage %>%
  filter(n_combos_covered == n_combos) %>%
  pull(name)

length(full_coverage_genes)

# violin plot of allelic ratio by sample, for a single gene of interest
# set target_gene to whichever gene name you want to look at
target_gene <- signal_genes$name[1]

single_gene_df <- per_cell_signal_genes %>%
  filter(name %in% target_gene & total_reads >= 10)
single_gene_df$celltype <- droplevels(single_gene_df$celltype)

pdf('Allelic_ratio_results/all_genes_top6_allelic_ratio_violin_by_sample.pdf', width = 16, height = 12)
ggplot(single_gene_df, aes(x = sample, y = allelic_ratio, fill = sample)) +
  geom_violin(trim = TRUE) +
  geom_jitter(aes(size = total_reads), width = 0.15, height = 0, alpha = 0.5) +
  scale_size_continuous(range = c(0.5, 4)) +
  labs(x = NULL, y = 'Allelic ratio', title = paste(target_gene, collapse = ", "), size = 'Total reads') +
  theme_bw() +
  guides(fill = "none") +
  facet_grid(name ~ celltype, labeller = labeller(celltype = label_wrap_gen(width = 14)))
dev.off()

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

# Repeat the dispersion/heterogeneity test after excluding putative LOX cells
# (AR >= 0.9) - a lost inactive allele collapses a cell's ratio to a fixed
# point with no sampling variance, which can trivially shift/mask the
# between-condition dispersion comparison for an escape gene.
lox_barcodes_ceb <- rownames(subset(metadata_ceb, monoallelic == 1))

metadata_whole_chr <- metadata_whole_chr[!rownames(metadata_whole_chr) %in% lox_barcodes_ceb, ]

adult_vs_aged_lrt_noLOX <- split(metadata_whole_chr, metadata_whole_chr$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("9w", "78w"))
    if (length(unique(df$sample)) < 2) return(NULL)
    run_dispersion_lrt(df, ref_level = "9w", comparison_level = "78w")
  })
adult_vs_aged_lrt_noLOX <- bind_rows(adult_vs_aged_lrt_noLOX, .id = "celltype")
adult_vs_aged_lrt_noLOX$FDR <- p.adjust(adult_vs_aged_lrt_noLOX$Pr..Chisq., method = "fdr")

subset(adult_vs_aged_lrt_noLOX, FDR < 0.05)

Sham_vs_TAC_lrt_noLOX <- split(metadata_whole_chr, metadata_whole_chr$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("Sham", "TAC"))
    if (length(unique(df$sample)) < 2) return(NULL)
    run_dispersion_lrt(df, ref_level = "Sham", comparison_level = "TAC")
  })
Sham_vs_TAC_lrt_noLOX <- bind_rows(Sham_vs_TAC_lrt_noLOX, .id = "celltype")
Sham_vs_TAC_lrt_noLOX$FDR <- p.adjust(Sham_vs_TAC_lrt_noLOX$Pr..Chisq., method = "fdr")

subset(Sham_vs_TAC_lrt_noLOX, FDR < 0.05)


