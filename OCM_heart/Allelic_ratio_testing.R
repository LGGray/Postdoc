library(Seurat)
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(broom.mixed)


# Colour scale
my_colors <- c("#1c642d", "#1c642d", "#0F7031", "#357B30", "#4E8330", "#658C2D", "#78962A" ,"#8D9F25", "#A2A71D", "#B3B112", "#C97314" ,"#8b1913")
my_breaks <- seq(0.5,1,by = 0.05)

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

heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

###############################
# Whole X chromosome analysis #
################################
barcodes <- lapply(c('9w', '78w', 'Sham', 'TAC'), function(x) {
  list.files(paste0('Allelome.PRO2/', x), pattern = 'locus_table.txt', recursive = TRUE, full.names = TRUE)
})
barcodes <- unlist(barcodes)

condition <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[4])
cellid <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[5])

allelic_ratio <- lapply(barcodes, function(x) {
  if (file.exists(x)) {
    tmp <- read.delim(x, header = TRUE)
    subset(tmp, chr=="chrX")$allelic_ratio
  } else {
    NULL
  }
})
names(allelic_ratio) <- paste0(condition, "_", cellid)
allelic_ratio <- allelic_ratio[!sapply(allelic_ratio, is.null)]
allelic_ratio <- unlist(allelic_ratio)

# Subset seurat object by barcodes
subset_heart <- subset(heart, cells = names(allelic_ratio))
allelic_ratio <- allelic_ratio[names(allelic_ratio) %in% colnames(subset_heart)]
allelic_ratio <- allelic_ratio[match(colnames(subset_heart), names(allelic_ratio))]
subset_heart$allelic_ratio <- allelic_ratio

total_reads <- lapply(barcodes, function(x) {
  if (file.exists(x)) {
    tmp <- read.delim(x, header = TRUE)
    subset(tmp, chr=="chrX")$total_reads
  } else {
    NULL
  }
})
total_reads <- unlist(total_reads)
names(total_reads) <- paste0(condition, "_", cellid)
total_reads <- total_reads[names(total_reads) %in% colnames(subset_heart)]
total_reads <- total_reads[match(colnames(subset_heart), names(total_reads))]
subset_heart$total_reads <- total_reads

A1_reads <- lapply(barcodes, function(x) {
  if (file.exists(x)) {
    tmp <- read.delim(x, header = TRUE)
    subset(tmp, chr=="chrX")$A1_reads
  } else {
    NULL
  }
})
A1_reads <- unlist(A1_reads)
names(A1_reads) <- paste0(condition, "_", cellid)
A1_reads <- A1_reads[names(A1_reads) %in% colnames(subset_heart)]
A1_reads <- A1_reads[match(colnames(subset_heart), names(A1_reads))]
subset_heart$A1_reads <- A1_reads

A2_reads <- lapply(barcodes, function(x) {
  if (file.exists(x)) {
    tmp <- read.delim(x, header = TRUE)
    subset(tmp, chr=="chrX")$A2_reads
  } else {
    NULL
  }
})
A2_reads <- unlist(A2_reads)
names(A2_reads) <- paste0(condition, "_", cellid)
A2_reads <- A2_reads[names(A2_reads) %in% colnames(subset_heart)]
A2_reads <- A2_reads[match(colnames(subset_heart), names(A2_reads))]
subset_heart$A2_reads <- A2_reads


# Plot distribution of total reads
pdf('Allelic_ratio_results/whole_chr_total_reads_distribution.pdf')
hist(subset_heart$total_reads, breaks = 10, main = "Total Reads Distribution", xlab = "Total Reads")
dev.off()
#
################################################### 
# filter for cells with at least 20 reads on chrX #
###################################################
subset_heart_flt <- subset(subset_heart, subset = total_reads >= 20)

pdf('Allelic_ratio_results/whole_chr_allelic_ratio_celltype_violin_plot_facet_wrap.pdf')
ggplot(subset_heart_flt@meta.data, aes(x = sample, y = allelic_ratio, fill = after_stat(y))) +
  geom_violin() +
  scale_fill_gradientn(colors = my_colors,
                       breaks = my_breaks,
                       limits = c(0.5, 1),
                       oob = scales::squish,
                       name = "Allelic ratio") +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# UMAP plot coloured by allelic ratio split by sample
pdf('Allelic_ratio_results/allelic_ratio_umap_plot.pdf')
FeaturePlot(subset_heart_flt,
            features = "allelic_ratio",
            split.by = "sample",
            cols = my_colors,
            min.cutoff = 0.5,
            max.cutoff = 1)
dev.off()

# Create metadata for statistical testing
metadata <- subset_heart_flt@meta.data


# number of cells per celltype and condition
cell_counts <- metadata %>%
  group_by(celltype, sample) %>%
  summarise(n_cells = n())
write.table(cell_counts, 'Allelic_ratio_results/whole_chr_cell_counts_per_celltype_and_condition.txt', sep = '\t', row.names = FALSE, quote = FALSE)


# Statistical testing - Does the variation in allelic ratio differ between conditions for each cell type?
adult_vs_aged_lrt <- split(metadata, metadata$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("9w", "78w"))
    if (length(unique(df$sample)) < 2) return(NULL)
    run_dispersion_lrt(df, ref_level = "9w", comparison_level = "78w")
  })
adult_vs_aged_lrt <- bind_rows(adult_vs_aged_lrt, .id = "celltype")
adult_vs_aged_lrt$FDR <- p.adjust(adult_vs_aged_lrt$Pr..Chisq., method = "fdr")

subset(adult_vs_aged_lrt, FDR < 0.05)

Sham_vs_TAC_lrt <- split(metadata, metadata$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("Sham", "TAC"))
    if (length(unique(df$sample)) < 2) return(NULL)
    run_dispersion_lrt(df, ref_level = "Sham", comparison_level = "TAC")
  })
Sham_vs_TAC_lrt <- bind_rows(Sham_vs_TAC_lrt, .id = "celltype")
Sham_vs_TAC_lrt$FDR <- p.adjust(Sham_vs_TAC_lrt$Pr..Chisq., method = "fdr")

subset(Sham_vs_TAC_lrt, FDR < 0.05)


# Plotting Fraction of biallelic-like cells per cell type and condition
fdr_to_stars <- function(fdr) {
  ifelse(is.na(fdr), "ns",
         ifelse(fdr <= 0.001, "***",
                ifelse(fdr <= 0.01, "**",
        ifelse(fdr <= 0.05, "*", "ns"))))
}

violin_tbl <- metadata %>%
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
ggplot(violin_tbl, aes(x = sample_idx, y = allelic_ratio, fill = after_stat(y))) +
  geom_violin(trim = FALSE, scale = "width") +
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
  scale_fill_gradientn(colors = my_colors,
                       breaks = my_breaks,
                       limits = c(0.5, 1),
                       oob = scales::squish,
                       name = "Allelic ratio") +
  facet_wrap(~celltype) +
  scale_x_continuous(breaks = 1:4, labels = levels(violin_tbl$sample)) +
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

save(subset_heart_flt, file = "Allelic_ratio_results/whole_chr_subset_heart_flt.RData")

##############################
# Core escape genes analysis #
##############################
barcodes <- lapply(c('9w', '78w', 'Sham', 'TAC'), function(x) {
  list.files(paste0('Allelome.PRO2_core_escape_genes/', x), pattern = 'locus_table.txt', recursive = TRUE, full.names = TRUE)
})
barcodes <- unlist(barcodes)

condition <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[7])
cellid <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[8])

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

# match cell types from Seurat object
gene_df$celltype <- Idents(heart)[gene_df$cell_barcode]

write.table(gene_df, 'Allelic_ratio_results/core_escape_genes_gene_df.txt', sep = '\t', row.names = FALSE, quote = FALSE)
gene_df <- read.table('Allelic_ratio_results/core_escape_genes_gene_df.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)


pdf('Allelic_ratio_results/core_escape_genes_total_reads_distribution.pdf')
hist(gene_df$total_reads, breaks = 2, main = "Total Reads Distribution", xlab = "Total Reads")
dev.off()


# Plot a PCA of allelic ratios across core escape genes coloured by cell type and shape by sample
pca_input <- gene_df %>%
  dplyr::select(cell_barcode, sample, celltype, name, allelic_ratio) %>%
  dplyr::distinct()

pca_wide <- reshape(
  as.data.frame(pca_input),
  idvar = c("cell_barcode", "sample", "celltype"),
  timevar = "name",
  direction = "wide"
)

gene_cols <- grep("^allelic_ratio\\.", colnames(pca_wide), value = TRUE)

if (length(gene_cols) >= 2 && nrow(pca_wide) >= 3) {
  # Use only complete observations (no imputation).
  complete_rows <- complete.cases(pca_wide[, gene_cols, drop = FALSE])
  pca_wide_cc <- pca_wide[complete_rows, , drop = FALSE]
  pca_mat <- as.matrix(pca_wide_cc[, gene_cols, drop = FALSE])

  # Remove columns with no variance.
  keep_cols <- apply(pca_mat, 2, function(x) stats::sd(x) > 0)
  pca_mat <- pca_mat[, keep_cols, drop = FALSE]

  if (ncol(pca_mat) >= 2 && nrow(pca_mat) >= 3) {
    pca_fit <- prcomp(pca_mat, center = TRUE, scale. = TRUE)
    pca_scores <- as.data.frame(pca_fit$x[, 1:2, drop = FALSE])
    pca_scores$celltype <- pca_wide_cc$celltype
    pca_scores$sample <- pca_wide_cc$sample

    var_exp <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
    pc1_lab <- paste0("PC1 (", round(100 * var_exp[1], 1), "%)")
    pc2_lab <- paste0("PC2 (", round(100 * var_exp[2], 1), "%)")

    pdf('Allelic_ratio_results/core_escape_genes_allelic_ratio_PCA.pdf')
    print(ggplot(pca_scores, aes(x = PC1, y = PC2, color = celltype, shape = sample)) +
      geom_point(alpha = 0.8, size = 1.6) +
      labs(x = pc1_lab, y = pc2_lab) +
      theme_bw())
    dev.off()
  }
}

# Heirachical bayesian approach

# Informative priors from bulk median allelic ratios.
# Aged medians available for Ddx3x, Kdm5c, Kdm6a; Eif2s3x uses adult value.
bulk_prior_tbl <- data.frame(
  name = c("Ddx3x", "Eif2s3x", "Kdm5c", "Kdm6a"),
  adult_median = c(0.68, 0.74, 0.68, 0.62),
  aged_median = c(0.68, 0.74, 0.65, 0.67),
  stringsAsFactors = FALSE
) %>%
  mutate(prior_mean = rowMeans(cbind(adult_median, aged_median), na.rm = TRUE))

# Prior concentration controls how strongly bulk estimates regularize sparse scRNA counts.
prior_kappa <- 20

gene_df_bayes <- gene_df %>%
  filter(name %in% bulk_prior_tbl$name,
         !is.na(A1_reads), !is.na(A2_reads), !is.na(total_reads), total_reads > 0)

if (nrow(gene_df_bayes) > 0) {
  bayes_gene <- gene_df_bayes %>%
    group_by(sample, celltype, name) %>%
    summarise(
      A1_reads = sum(A1_reads, na.rm = TRUE),
      A2_reads = sum(A2_reads, na.rm = TRUE),
      total_reads = sum(total_reads, na.rm = TRUE),
      n_cells = n_distinct(cell_barcode),
      .groups = "drop"
    ) %>%
    left_join(bulk_prior_tbl, by = "name") %>%
    mutate(
      alpha_prior = prior_mean * prior_kappa,
      beta_prior = (1 - prior_mean) * prior_kappa,
      post_alpha = alpha_prior + A1_reads,
      post_beta = beta_prior + A2_reads,
      post_mean = post_alpha / (post_alpha + post_beta),
      post_median = qbeta(0.5, post_alpha, post_beta),
      post_lower95 = qbeta(0.025, post_alpha, post_beta),
      post_upper95 = qbeta(0.975, post_alpha, post_beta),
      post_prob_gt_prior = 1 - pbeta(prior_mean, post_alpha, post_beta),
      post_prob_gt_0.9 = 1 - pbeta(0.9, post_alpha, post_beta)
    )

  bayes_gene$sample <- factor(bayes_gene$sample, levels = c("9w", "78w", "Sham", "TAC"))

  

  write.table(
    bayes_gene,
    'Allelic_ratio_results/core_escape_genes_bayes_posterior_by_gene.txt',
    sep = '\t', row.names = FALSE, quote = FALSE
  )



  bayes_celltype <- bayes_gene %>%
    group_by(sample, celltype) %>%
    summarise(
      A1_reads = sum(A1_reads, na.rm = TRUE),
      A2_reads = sum(A2_reads, na.rm = TRUE),
      total_reads = sum(total_reads, na.rm = TRUE),
      n_cells = sum(n_cells, na.rm = TRUE),
      prior_mean = mean(prior_mean, na.rm = TRUE),
      alpha_prior = sum(alpha_prior, na.rm = TRUE),
      beta_prior = sum(beta_prior, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      post_alpha = alpha_prior + A1_reads,
      post_beta = beta_prior + A2_reads,
      post_mean = post_alpha / (post_alpha + post_beta),
      post_median = qbeta(0.5, post_alpha, post_beta),
      post_lower95 = qbeta(0.025, post_alpha, post_beta),
      post_upper95 = qbeta(0.975, post_alpha, post_beta),
      post_prob_gt_prior = 1 - pbeta(prior_mean, post_alpha, post_beta),
      post_prob_gt_0.9 = 1 - pbeta(0.9, post_alpha, post_beta)
    )

  bayes_celltype$sample <- factor(bayes_celltype$sample, levels = c("9w", "78w", "Sham", "TAC"))

  write.table(
    bayes_celltype,
    'Allelic_ratio_results/core_escape_genes_bayes_posterior_by_celltype.txt',
    sep = '\t', row.names = FALSE, quote = FALSE
  )

  pdf('Allelic_ratio_results/core_escape_genes_bayes_posterior_mean_heatmap_by_gene.pdf')
  ggplot(bayes_gene, aes(x = sample, y = celltype, fill = post_mean)) +
    geom_tile() +
    facet_wrap(~name) +
    scale_fill_gradientn(colors = my_colors,
                         breaks = my_breaks,
                         limits = c(0.5, 1),
                         oob = scales::squish) +
    labs(fill = 'Posterior mean', x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()

  pdf('Allelic_ratio_results/core_escape_genes_bayes_posterior_credible_width_by_gene.pdf')
  ggplot(
    bayes_gene %>% mutate(ci_width = post_upper95 - post_lower95),
    aes(x = sample, y = celltype, fill = ci_width)
  ) +
    geom_tile() +
    facet_wrap(~name) +
    scale_fill_gradient(low = 'white', high = '#ef8a62') +
    labs(fill = '95% CrI width', x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()
}

# PCA for bayes_gene posterior means
if (exists("bayes_gene") && nrow(bayes_gene) > 0) {
  pca_input_bayes <- bayes_gene %>%
    dplyr::select(sample, celltype, name, post_mean) %>%
    dplyr::distinct()

  pca_wide_bayes <- reshape(
    as.data.frame(pca_input_bayes),
    idvar = c("sample", "celltype"),
    timevar = "name",
    direction = "wide"
  )

  gene_cols_bayes <- grep("^post_mean\\.", colnames(pca_wide_bayes), value = TRUE)

  if (length(gene_cols_bayes) >= 2 && nrow(pca_wide_bayes) >= 3) {
    complete_rows_bayes <- complete.cases(pca_wide_bayes[, gene_cols_bayes, drop = FALSE])
    pca_wide_bayes_cc <- pca_wide_bayes[complete_rows_bayes, , drop = FALSE]
    pca_mat_bayes <- as.matrix(pca_wide_bayes_cc[, gene_cols_bayes, drop = FALSE])

    keep_cols_bayes <- apply(pca_mat_bayes, 2, function(x) stats::sd(x) > 0)
    pca_mat_bayes <- pca_mat_bayes[, keep_cols_bayes, drop = FALSE]

    if (ncol(pca_mat_bayes) >= 2 && nrow(pca_mat_bayes) >= 3) {
      pca_fit_bayes <- prcomp(pca_mat_bayes, center = TRUE, scale. = TRUE)
      pca_scores_bayes <- as.data.frame(pca_fit_bayes$x[, 1:2, drop = FALSE])
      pca_scores_bayes$celltype <- pca_wide_bayes_cc$celltype
      pca_scores_bayes$sample <- pca_wide_bayes_cc$sample

      var_exp_bayes <- (pca_fit_bayes$sdev^2) / sum(pca_fit_bayes$sdev^2)
      pc1_lab_bayes <- paste0("PC1 (", round(100 * var_exp_bayes[1], 1), "%)")
      pc2_lab_bayes <- paste0("PC2 (", round(100 * var_exp_bayes[2], 1), "%)")

      pdf('Allelic_ratio_results/core_escape_genes_bayes_posterior_PCA.pdf')
      print(ggplot(pca_scores_bayes, aes(x = PC1, y = PC2, color = celltype, shape = sample)) +
        geom_point(alpha = 0.8, size = 2) +
        labs(x = pc1_lab_bayes, y = pc2_lab_bayes) +
        theme_bw())
      dev.off()
    }
  }
}


##############################
## Pseudobulk and cell-level LOX inference for low-depth data ##
##############################
gene_df <- gene_df %>%
  filter(!is.na(A1_reads), !is.na(A2_reads), !is.na(total_reads), total_reads > 0)

if (nrow(gene_df) > 0) {
  # Pseudobulk at sample x celltype x gene level
  pb_gene <- gene_df %>%
    group_by(sample, celltype, name) %>%
    summarise(
      A1_reads = sum(A1_reads, na.rm = TRUE),
      A2_reads = sum(A2_reads, na.rm = TRUE),
      total_reads = sum(total_reads, na.rm = TRUE),
      allelic_ratio = A1_reads / pmax(total_reads, 1),
      n_cells = n_distinct(cell_barcode),
      .groups = "drop"
    )
  pb_gene$sample <- factor(pb_gene$sample, levels = c("9w", "78w", "Sham", "TAC"))
  write.table(pb_gene,
              'Allelic_ratio_results/core_escape_genes_pseudobulk_by_gene.txt',
              sep = '\t', row.names = FALSE, quote = FALSE)

  # Pseudobulk at sample x celltype level across all selected genes
  pb_celltype <- gene_df %>%
    group_by(sample, celltype) %>%
    summarise(
      A1_reads = sum(A1_reads, na.rm = TRUE),
      A2_reads = sum(A2_reads, na.rm = TRUE),
      total_reads = sum(total_reads, na.rm = TRUE),
      allelic_ratio = A1_reads / pmax(total_reads, 1),
      n_cells = n_distinct(cell_barcode),
      n_gene_obs = n(),
      .groups = "drop"
    )
  pb_celltype$sample <- factor(pb_celltype$sample, levels = c("9w", "78w", "Sham", "TAC"))
  write.table(pb_celltype,
              'Allelic_ratio_results/core_escape_genes_pseudobulk_by_celltype.txt',
              sep = '\t', row.names = FALSE, quote = FALSE)

  compare_pseudobulk <- function(df, g1, g2) {
    split(df, df$celltype) %>%
      lapply(function(x) {
        x <- subset(x, sample %in% c(g1, g2))
        if (length(unique(x$sample)) < 2) return(NULL)
        x <- x[match(c(g1, g2), x$sample), ]
        if (any(is.na(x$sample))) return(NULL)

        mat <- matrix(c(x$A1_reads[1], x$A2_reads[1], x$A1_reads[2], x$A2_reads[2]),
                      nrow = 2, byrow = TRUE)
        ft <- fisher.test(mat)

        data.frame(
          celltype = x$celltype[1],
          group1 = g1,
          group2 = g2,
          A1_group1 = x$A1_reads[1],
          A2_group1 = x$A2_reads[1],
          A1_group2 = x$A1_reads[2],
          A2_group2 = x$A2_reads[2],
          OR_group2_vs_group1 = unname(ft$estimate),
          p.value = ft$p.value,
          stringsAsFactors = FALSE
        )
      }) %>%
      bind_rows()
  }

  pb_9w_78w <- compare_pseudobulk(pb_celltype, "9w", "78w")
  if (nrow(pb_9w_78w) > 0) {
    pb_9w_78w$FDR <- p.adjust(pb_9w_78w$p.value, method = "fdr")
    write.table(pb_9w_78w,
                'Allelic_ratio_results/core_escape_genes_pseudobulk_9w_vs_78w.txt',
                sep = '\t', row.names = FALSE, quote = FALSE)
  }

  pb_sham_tac <- compare_pseudobulk(pb_celltype, "Sham", "TAC")
  if (nrow(pb_sham_tac) > 0) {
    pb_sham_tac$FDR <- p.adjust(pb_sham_tac$p.value, method = "fdr")
    write.table(pb_sham_tac,
                'Allelic_ratio_results/core_escape_genes_pseudobulk_Sham_vs_TAC.txt',
                sep = '\t', row.names = FALSE, quote = FALSE)
  }

  # Pseudobulk visualisation
  pdf('Allelic_ratio_results/core_escape_genes_pseudobulk_heatmap_by_gene.pdf')
  ggplot(pb_gene, aes(x = sample, y = celltype, fill = allelic_ratio)) +
    geom_tile() +
    facet_wrap(~name) +
    scale_fill_gradientn(colors = my_colors,
                         breaks = my_breaks,
                         limits = c(0.5, 1),
                         oob = scales::squish) +
    labs(fill = 'A1 / total', x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()

  pdf('Allelic_ratio_results/core_escape_genes_pseudobulk_heatmap_pooled.pdf')
  ggplot(pb_celltype, aes(x = sample, y = celltype, fill = allelic_ratio)) +
    geom_tile() +
    scale_fill_gradientn(colors = my_colors,
                         breaks = my_breaks,
                         limits = c(0.5, 1),
                         oob = scales::squish) +
    labs(fill = 'A1 / total', x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()

  pb_long <- bind_rows(
    pb_celltype %>% transmute(sample, celltype, allele = 'A1', reads = A1_reads),
    pb_celltype %>% transmute(sample, celltype, allele = 'A2', reads = A2_reads)
  )
  pb_long$sample <- factor(pb_long$sample, levels = c("9w", "78w", "Sham", "TAC"))

  pdf('Allelic_ratio_results/core_escape_genes_pseudobulk_A1_A2_stacked.pdf')
  ggplot(pb_long, aes(x = sample, y = reads, fill = allele)) +
    geom_col(position = 'fill') +
    facet_wrap(~celltype) +
    labs(y = 'Proportion of pooled reads', x = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()

  comp_tbl <- bind_rows(
    if (nrow(pb_9w_78w) > 0) pb_9w_78w %>% mutate(comparison = '9w vs 78w') else NULL,
    if (nrow(pb_sham_tac) > 0) pb_sham_tac %>% mutate(comparison = 'Sham vs TAC') else NULL
  )

  if (nrow(comp_tbl) > 0) {
    comp_tbl <- comp_tbl %>%
      mutate(
        log2OR = log2(OR_group2_vs_group1),
        sig = ifelse(FDR <= 0.05, 'FDR<=0.05', 'ns')
      )

    pdf('Allelic_ratio_results/core_escape_genes_pseudobulk_OR_forest.pdf')
    ggplot(comp_tbl, aes(x = log2OR, y = celltype, color = sig)) +
      geom_vline(xintercept = 0, linetype = 2, color = 'grey50') +
      geom_point(size = 2) +
      facet_wrap(~comparison) +
      labs(x = 'log2(OR group2/group1)', y = NULL) +
      theme_bw()
    dev.off()
  }

  # Cell-level posterior probability of LOX across selected genes
  # LOX proxy: posterior P(true maternal fraction > 0.9) using Beta(1,1) prior.
  cell_agg <- gene_df %>%
    group_by(cell_barcode, sample, celltype) %>%
    summarise(
      A1_reads = sum(A1_reads, na.rm = TRUE),
      A2_reads = sum(A2_reads, na.rm = TRUE),
      total_reads = sum(total_reads, na.rm = TRUE),
      n_genes_detected = n_distinct(name),
      .groups = "drop"
    ) %>%
    mutate(
      post_prob_lox = 1 - pbeta(0.9, A1_reads + 1, A2_reads + 1),
      lox_high_conf = post_prob_lox >= 0.9
    )

  write.table(cell_agg,
              'Allelic_ratio_results/core_escape_genes_cell_level_posterior_LOX.txt',
              sep = '\t', row.names = FALSE, quote = FALSE)

  lox_prop <- cell_agg %>%
    group_by(celltype, sample) %>%
    summarise(
      n = n(),
      n_lox_high_conf = sum(lox_high_conf, na.rm = TRUE),
      p_lox_high_conf = n_lox_high_conf / n,
      se = sqrt(p_lox_high_conf * (1 - p_lox_high_conf) / n),
      .groups = "drop"
    )

  pdf('Allelic_ratio_results/core_escape_genes_fraction_high_conf_LOX_cells.pdf')
  ggplot(lox_prop, aes(x = sample, y = p_lox_high_conf, color = sample)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = pmax(0, p_lox_high_conf - 1.96 * se),
                      ymax = pmin(1, p_lox_high_conf + 1.96 * se)),
                  width = 0.2) +
    facet_wrap(~celltype) +
    labs(y = 'Fraction high-confidence LOX cells', x = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none')
  dev.off()
}

saveRDS(gene_df, file = "Allelic_ratio_results/core_escape_genes_gene_df.RData")