# ---------------------------------------------------------------------------
# 04 - Core escape gene analysis (block-level, whole chromosome).
#
# Run from the OCM_heart/ directory:  Rscript allelic_ratio/04_core_escape.R
# Requires: 01_setup.R to have been run.
# Writes:   Allelic_ratio_results/core_escape_block_cell_metadata.txt for 05.
# ---------------------------------------------------------------------------
source("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/OCM_heart/allelic_ratio/00_functions.R")

heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

#################################################
# Core escape genes analysis (block, whole-chr) #
#################################################
barcodes <- lapply(c('9w', '78w', 'Sham', 'TAC'), function(x) {
  list.files(paste0('Allelome.PRO2_core_escape_new/', x), pattern = 'locus_table.txt', recursive = TRUE, full.names = TRUE)
})
barcodes <- unlist(barcodes)

condition <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[7])
cellid <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[8])

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
write.table(core_escape_block_ratio, 'Allelic_ratio_results/core_escape_block_new_allelic_ratio_table.txt', sep = '\t', row.names = FALSE, quote = FALSE)

# Subset seurat object by barcodes
subset_heart_ceb <- subset(heart, cells = core_escape_block_ratio$cell_barcode)
core_escape_block_ratio <- core_escape_block_ratio[core_escape_block_ratio$cell_barcode %in% colnames(subset_heart_ceb), ]
core_escape_block_ratio <- core_escape_block_ratio[match(colnames(subset_heart_ceb), core_escape_block_ratio$cell_barcode), ]
subset_heart_ceb$allelic_ratio <- core_escape_block_ratio$allelic_ratio
subset_heart_ceb$total_reads <- core_escape_block_ratio$total_reads
subset_heart_ceb$A1_reads <- core_escape_block_ratio$A1_reads
subset_heart_ceb$A2_reads <- core_escape_block_ratio$A2_reads

# Plot distribution of total reads
pdf('Allelic_ratio_results/core_escape_block_new_total_reads_distribution.pdf')
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
my_colors <- c(
  "#2B3186", "#374795", "#3B5FB6", "#38749F",
  "#37758B", "#367373", "#2D6E5D", "#2A7050",
  "#1E652D", "#1C642D", "#0F7031", "#357B30", "#4E8330",
  "#658C2D", "#78962A", "#8D9F25", "#A2A71D", "#B3B112",
  "#C97314", "#8B1913"
)
my_breaks  <- seq(0, 1, by = 0.05)
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
pdf('Allelic_ratio_results/core_escape_block_new_allelic_ratio_umap_plot_split_by_sample.pdf', width = 10, height = 10)
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

pdf('Allelic_ratio_results/core_escape_block_new_allelic_ratio_celltype_violin_plot_facet_wrap.pdf')
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

pdf('Allelic_ratio_results/core_escape_block_new_LOX_cell_counts_stacked_barplot.pdf')
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
saveRDS(subset_heart_ceb, 'Allelic_ratio_results/subset_heart_core_escape_block_new.RDS')


# ==========================================================================
# Beta-binomial model of allelic ratio ~ Xist expression, per sample x celltype.
xist_ar <- data.frame(
  Xist        = GetAssayData(subset_heart_ceb_flt, assay = "SCT", layer = "data")["Xist", ],
  A1_reads    = subset_heart_ceb_flt$A1_reads,
  A2_reads    = subset_heart_ceb_flt$A2_reads,
  total_reads = subset_heart_ceb_flt$total_reads,
  allelic_ratio = subset_heart_ceb_flt$allelic_ratio,
  sample      = subset_heart_ceb_flt$sample,
  celltype    = subset_heart_ceb_flt$celltype
)

# Fit one beta-binomial model per sample x celltype. Keep the fitted models so
# the stats table and the fitted curves in the plots come from the same fit.
bb_split <- split(xist_ar, list(xist_ar$celltype, xist_ar$sample), drop = TRUE)
bb_split <- bb_split[sapply(bb_split, function(x) nrow(x) >= 20 && length(unique(x$Xist)) >= 5)]

bb_models <- lapply(bb_split, function(x) {
  m <- try(glmmTMB(cbind(A1_reads, A2_reads) ~ Xist,
                   family = betabinomial, data = x), silent = TRUE)
  if (inherits(m, "try-error")) NULL else m
})
keep_bb   <- !sapply(bb_models, is.null)
bb_split  <- bb_split[keep_bb]
bb_models <- bb_models[keep_bb]

AR_Xist_bb <- Map(function(m, x) {
  s <- summary(m)$coefficients$cond["Xist", ]
  data.frame(sample = unique(x$sample), celltype = unique(x$celltype),
             beta = s[1], se = s[2], OR = exp(s[1]),
             p_value = s[4], n_cells = nrow(x))
}, bb_models, bb_split) %>% bind_rows()
rownames(AR_Xist_bb) <- NULL
AR_Xist_bb$FDR <- p.adjust(AR_Xist_bb$p_value, method = "fdr")

write.table(AR_Xist_bb, 'Allelic_ratio_results/core_escape_block_new_Xist_vs_AR_betabinomial.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# Plots for the beta-binomial fit.
#
# Xist on x, allelic ratio on y (Xist is the predictor, AR the modelled
# response). Point size = informative read depth, so the pile-up of cells at
# AR = 0 / AR = 1 is visibly the low-depth end. Fitted curve + 95% CI comes
# from the beta-binomial model; the overlaid points with error bars are pooled
# counts within Xist quintiles (sum(A1) / sum(total), Wilson interval), which
# is an assumption-light check on the model curve.
# ---------------------------------------------------------------------------

# Fitted curve + 95% CI on the response scale, per modelled group.
bb_pred <- Map(function(m, x) {
  nd <- data.frame(Xist = seq(min(x$Xist), max(x$Xist), length.out = 100),
                   A1_reads = 1, A2_reads = 1)
  p <- predict(m, newdata = nd, type = "link", se.fit = TRUE)
  data.frame(sample = unique(x$sample), celltype = unique(x$celltype),
             Xist = nd$Xist,
             fit = plogis(p$fit),
             lwr = plogis(p$fit - 1.96 * p$se.fit),
             upr = plogis(p$fit + 1.96 * p$se.fit))
}, bb_models, bb_split) %>% bind_rows()

# Pooled allelic ratio within Xist quintiles, with a Wilson interval.
wilson_ci <- function(x, n, z = 1.96) {
  p <- x / n
  d <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / d
  halfw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  data.frame(est = p, lwr = pmax(0, centre - halfw), upr = pmin(1, centre + halfw))
}

bb_bins <- lapply(bb_split, function(x) {
  br <- unique(quantile(x$Xist, probs = seq(0, 1, 0.2), na.rm = TRUE))
  x$bin <- if (length(br) < 3) 1L else cut(x$Xist, breaks = br, include.lowest = TRUE, labels = FALSE)
  binned <- x %>%
    group_by(bin) %>%
    summarise(Xist = mean(Xist), A1 = sum(A1_reads), N = sum(total_reads), .groups = "drop") %>%
    filter(N > 0)
  binned$sample   <- unique(x$sample)
  binned$celltype <- unique(x$celltype)
  cbind(binned, wilson_ci(binned$A1, binned$N))
}) %>% bind_rows()

# Significance flag drives the colour of the fitted curve, so significant
# panels are readable at a glance in the full grid.
sig_lookup <- AR_Xist_bb %>%
  mutate(sig = FDR < 0.05) %>%
  dplyr::select(sample, celltype, sig)

bb_pred <- left_join(bb_pred, sig_lookup, by = c("sample", "celltype"))
bb_bins <- left_join(bb_bins, sig_lookup, by = c("sample", "celltype"))

# Stats annotation: OR and stars for the significant groups, "ns" otherwise.
bb_ann <- AR_Xist_bb %>%
  mutate(sig   = FDR < 0.05,
         label = ifelse(sig,
                        paste0("OR=", sprintf("%.2f", OR), " ", fdr_to_stars(FDR),
                               "\nFDR=", format.pval(FDR, digits = 2), "  n=", n_cells),
                        paste0("ns  n=", n_cells)))

sig_cols <- c(`TRUE` = "firebrick", `FALSE` = "grey40")

pdf('Allelic_ratio_results/core_escape_block_new_Xist_vs_AR_betabinomial.pdf', width = 11, height = 14)
ggplot(xist_ar, aes(x = Xist, y = allelic_ratio)) +
  geom_point(aes(size = total_reads), alpha = 0.3, colour = "steelblue") +
  geom_ribbon(data = bb_pred, aes(y = fit, ymin = lwr, ymax = upr, fill = sig),
              alpha = 0.2, colour = NA) +
  geom_line(data = bb_pred, aes(y = fit, colour = sig), linewidth = 0.7) +
  geom_errorbar(data = bb_bins, aes(y = est, ymin = lwr, ymax = upr),
                width = 0, colour = "black", linewidth = 0.4) +
  geom_point(data = bb_bins, aes(y = est), size = 1.6, colour = "black") +
  geom_text(data = bb_ann, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.15, size = 2.4) +
  facet_grid(celltype ~ sample, labeller = labeller(celltype = label_wrap_gen(18))) +
  scale_size_continuous(range = c(0.3, 2.5), name = "Informative\nreads") +
  scale_colour_manual(values = sig_cols, guide = "none") +
  scale_fill_manual(values = sig_cols, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Xist expression (SCT, normalized)",
       y = "Allelic ratio (core escape genes)",
       caption = "Curve: beta-binomial fit +/- 95% CI (red = FDR < 0.05). Black points: pooled AR within Xist quintiles, Wilson 95% CI.") +
  theme_bw() +
  theme(strip.text.y = element_text(size = 7),
        plot.caption = element_text(size = 7, hjust = 0))
dev.off()

# Same, Ventricular Cardiomyocytes only, faceted by sample.
vcm_pts  <- filter(xist_ar,    celltype == "Ventricular Cardiomyocytes")
vcm_pred <- filter(bb_pred,    celltype == "Ventricular Cardiomyocytes")
vcm_bins <- filter(bb_bins,    celltype == "Ventricular Cardiomyocytes")
vcm_ann_bb <- filter(bb_ann,   celltype == "Ventricular Cardiomyocytes")

pdf('Allelic_ratio_results/core_escape_block_new_Xist_vs_AR_betabinomial_VCM.pdf', width = 11, height = 3.5)
ggplot(vcm_pts, aes(x = Xist, y = allelic_ratio)) +
  geom_point(aes(size = total_reads), alpha = 0.3, colour = "steelblue") +
  geom_ribbon(data = vcm_pred, aes(y = fit, ymin = lwr, ymax = upr, fill = sig),
              alpha = 0.2, colour = NA) +
  geom_line(data = vcm_pred, aes(y = fit, colour = sig), linewidth = 0.8) +
  geom_errorbar(data = vcm_bins, aes(y = est, ymin = lwr, ymax = upr),
                width = 0, colour = "black", linewidth = 0.4) +
  geom_point(data = vcm_bins, aes(y = est), size = 1.8, colour = "black") +
  geom_text(data = vcm_ann_bb, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.15, size = 3) +
  facet_wrap(~sample, nrow = 1) +
  scale_size_continuous(range = c(0.3, 3), name = "Informative\nreads") +
  scale_colour_manual(values = sig_cols, guide = "none") +
  scale_fill_manual(values = sig_cols, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Xist expression (SCT, normalized)",
       y = "Allelic ratio\n(core escape genes)",
       title = "Ventricular Cardiomyocytes",
       caption = "Curve: beta-binomial fit +/- 95% CI (red = FDR < 0.05). Black points: pooled AR within Xist quintiles, Wilson 95% CI.") +
  theme_bw() +
  theme(plot.caption = element_text(size = 7, hjust = 0))
dev.off()

# Effect-size summary across all groups: OR with 95% CI on the log-odds scale,
# significant groups highlighted. Compact alternative to reading the table.
bb_forest <- AR_Xist_bb %>%
  mutate(sig = FDR < 0.05,
         lwr = exp(beta - 1.96 * se),
         upr = exp(beta + 1.96 * se),
         celltype = factor(celltype, levels = rev(sort(unique(celltype)))))

pdf('Allelic_ratio_results/core_escape_block_new_Xist_vs_AR_betabinomial_forest.pdf', width = 8, height = 5)
ggplot(bb_forest, aes(x = OR, y = celltype, colour = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, linewidth = 0.5) +
  geom_point(aes(size = n_cells)) +
  facet_wrap(~sample, nrow = 1) +
  scale_x_log10() +
  scale_colour_manual(values = sig_cols, name = NULL,
                      labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "ns")) +
  scale_size_continuous(range = c(1, 3.5), name = "Cells") +
  labs(x = "Odds ratio per unit Xist (95% CI)", y = NULL,
       caption = "OR < 1: higher Xist associated with lower allelic ratio, i.e. AR closer to 0.5 and\nmore balanced biallelic expression of the core escape genes (AR = 1 is maternal/active allele only).") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7),
        plot.caption = element_text(size = 7, hjust = 0))
dev.off()


# ===========================================================================
# Diagnostics for the extreme-AR cells.
#
# These mice carry an Xist repA deletion on the maternal allele, so the
# maternal X cannot initiate silencing and is always the active X. AR is
# therefore Xa / (Xa + Xi) with a fixed allele assignment:
#   AR -> 1   all reads from the active X: escape lost, or the Xi is gone.
#   AR -> 0   all reads from the inactive X and none from the active X. Such a
#             cell would be nullisomic for X expression and is not viable, so
#             these are artefact by construction.
#
# The AR ~ 0 cells are therefore a class where the truth is known to be wrong,
# which makes them a calibration set for the false-positive rate on the
# AR >= 0.9 LOX calls.
# ===========================================================================

# --- 1. Observed vs expected frequency of extreme AR under the fitted model --
#
# NOTE ON PARAMETERISATION: glmmTMB's betabinomial uses Beta(mu*phi,
# (1-mu)*phi) with phi = sigma(model). If that ever changes upstream, the
# expected counts below become wrong, so the pmf check underneath is a guard.
dbetabinom <- function(k, n, mu, phi) {
  a <- mu * phi
  b <- (1 - mu) * phi
  exp(lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b))
}
stopifnot(abs(sum(dbetabinom(0:20, 20, 0.65, 3)) - 1) < 1e-8)

# For each cell: probability of observing zero maternal reads (AR = 0) and of
# observing all maternal reads (AR = 1), given its own depth and fitted mean.
# The binomial columns ignore overdispersion and show how much of the extremes
# the beta-binomial has already absorbed.
extreme_expected <- Map(function(m, x) {
  mu  <- fitted(m)
  phi <- sigma(m)
  n   <- x$total_reads
  data.frame(
    sample = unique(x$sample), celltype = unique(x$celltype), n_cells = nrow(x),
    obs_AR0  = sum(x$A1_reads == 0),
    exp_AR0_bb  = sum(dbetabinom(0, n, mu, phi)),
    exp_AR0_bin = sum(dbinom(0, n, mu)),
    obs_AR1  = sum(x$A2_reads == 0),
    exp_AR1_bb  = sum(dbetabinom(n, n, mu, phi)),
    exp_AR1_bin = sum(dbinom(n, n, mu))
  )
}, bb_models, bb_split) %>% bind_rows()
rownames(extreme_expected) <- NULL

extreme_expected <- extreme_expected %>%
  mutate(
    # > 1 means more extreme cells than the noise model can generate.
    excess_AR0 = obs_AR0 / exp_AR0_bb,
    excess_AR1 = obs_AR1 / exp_AR1_bb,
    # Transfer the AR = 0 artefact rate across to AR = 1, scaling by the
    # sampling asymmetry the model expects (AR = 1 is intrinsically more likely
    # than AR = 0 whenever mean AR > 0.5, so the rates are not interchangeable).
    excess_AR0_cells    = pmax(0, obs_AR0 - exp_AR0_bb),
    # Scaling is unstable when the model expects almost no AR = 0 cells, so
    # report NA there rather than a number that looks precise and is not.
    predicted_artefact_AR1 = ifelse(exp_AR0_bb < 0.5, NA_real_,
                                    excess_AR0_cells * exp_AR1_bb / exp_AR0_bb),
    artefact_frac_of_AR1   = ifelse(obs_AR1 == 0, NA_real_,
                                    predicted_artefact_AR1 / obs_AR1)
  )

# Guard against a wrong sigma() parameterisation upstream: an expected count
# can never exceed the number of cells it was summed over.
stopifnot(all(extreme_expected$exp_AR0_bb <= extreme_expected$n_cells + 1e-6),
          all(extreme_expected$exp_AR1_bb <= extreme_expected$n_cells + 1e-6))

write.table(extreme_expected, 'Allelic_ratio_results/core_escape_block_new_extreme_AR_observed_vs_expected.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

# --- 2. Read depth of the AR ~ 0 cells vs the rest ---------------------------
# If the AR ~ 0 cells are the binomial noise floor they should be concentrated
# at low depth. If their depth matches the rest, a depth threshold will not
# remove them and a second artefact process (ambient RNA, misassignment,
# doublets) is implicated.
depth_extreme <- xist_ar %>%
  mutate(ar_class = case_when(allelic_ratio <= 0.10 ~ "AR <= 0.10",
                              allelic_ratio >= 0.90 ~ "AR >= 0.90",
                              TRUE ~ "Intermediate")) %>%
  group_by(sample, ar_class) %>%
  summarise(n_cells = n(),
            median_reads = median(total_reads),
            mean_reads = mean(total_reads),
            q25 = quantile(total_reads, 0.25),
            q75 = quantile(total_reads, 0.75),
            .groups = "drop")

depth_test <- xist_ar %>%
  mutate(is_ar0 = allelic_ratio <= 0.10) %>%
  split(.$sample) %>%
  lapply(function(x) {
    if (length(unique(x$is_ar0)) < 2) return(NULL)
    tt <- wilcox.test(total_reads ~ is_ar0, data = x)
    data.frame(sample = unique(x$sample),
               median_reads_other = median(x$total_reads[!x$is_ar0]),
               median_reads_AR0   = median(x$total_reads[x$is_ar0]),
               n_AR0 = sum(x$is_ar0), p_value = tt$p.value)
  }) %>% bind_rows()
if (nrow(depth_test)) depth_test$FDR <- p.adjust(depth_test$p_value, method = "fdr")

write.table(depth_extreme, 'Allelic_ratio_results/core_escape_block_new_extreme_AR_depth_summary.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)
write.table(depth_test, 'Allelic_ratio_results/core_escape_block_new_extreme_AR_depth_test.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

pdf('Allelic_ratio_results/core_escape_block_new_extreme_AR_depth.pdf', width = 8, height = 4)
xist_ar %>%
  mutate(ar_class = case_when(allelic_ratio <= 0.10 ~ "AR <= 0.10",
                              allelic_ratio >= 0.90 ~ "AR >= 0.90",
                              TRUE ~ "Intermediate")) %>%
  ggplot(aes(x = ar_class, y = total_reads, fill = ar_class)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
  facet_wrap(~sample, nrow = 1) +
  scale_y_log10() +
  labs(x = NULL, y = "Informative reads per cell (log10)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
dev.off()

# --- 3. Sensitivity refit excluding the AR <= 0.10 cells ---------------------
# Those cells sit at logit(-Inf) and drag the curve down wherever they fall on
# the Xist axis. If the VCM result depends on them it is not a real result.
xist_ar_sens <- filter(xist_ar, allelic_ratio > 0.10)

sens_split <- split(xist_ar_sens, list(xist_ar_sens$celltype, xist_ar_sens$sample), drop = TRUE)
sens_split <- sens_split[sapply(sens_split, function(x) nrow(x) >= 20 && length(unique(x$Xist)) >= 5)]

sens_models <- lapply(sens_split, function(x) {
  m <- try(glmmTMB(cbind(A1_reads, A2_reads) ~ Xist,
                   family = betabinomial, data = x), silent = TRUE)
  if (inherits(m, "try-error")) NULL else m
})
keep_sens   <- !sapply(sens_models, is.null)
sens_split  <- sens_split[keep_sens]
sens_models <- sens_models[keep_sens]

AR_Xist_bb_sens <- Map(function(m, x) {
  s <- summary(m)$coefficients$cond["Xist", ]
  data.frame(sample = unique(x$sample), celltype = unique(x$celltype),
             OR_sens = exp(s[1]), p_sens = s[4], n_cells_sens = nrow(x))
}, sens_models, sens_split) %>% bind_rows()
rownames(AR_Xist_bb_sens) <- NULL
AR_Xist_bb_sens$FDR_sens <- p.adjust(AR_Xist_bb_sens$p_sens, method = "fdr")

# Side-by-side comparison: full fit vs AR > 0.10 fit.
bb_sensitivity <- AR_Xist_bb %>%
  dplyr::select(sample, celltype, OR, FDR, n_cells) %>%
  left_join(AR_Xist_bb_sens, by = c("sample", "celltype")) %>%
  mutate(n_dropped   = n_cells - n_cells_sens,
         OR_change   = OR_sens - OR,
         still_sig   = FDR < 0.05 & FDR_sens < 0.05,
         lost_sig    = FDR < 0.05 & FDR_sens >= 0.05)

write.table(bb_sensitivity, 'Allelic_ratio_results/core_escape_block_new_Xist_vs_AR_betabinomial_sensitivity.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

pdf('Allelic_ratio_results/core_escape_block_new_Xist_vs_AR_betabinomial_sensitivity.pdf', width = 7, height = 5)
ggplot(bb_sensitivity, aes(x = OR, y = OR_sens, colour = FDR < 0.05)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 1, colour = "grey85") +
  geom_hline(yintercept = 1, colour = "grey85") +
  geom_point(aes(size = n_cells), alpha = 0.8) +
  facet_wrap(~sample, nrow = 1) +
  scale_colour_manual(values = sig_cols, name = NULL,
                      labels = c(`TRUE` = "FDR < 0.05 (full fit)", `FALSE` = "ns")) +
  scale_size_continuous(range = c(1, 3.5), name = "Cells") +
  labs(x = "OR, all cells", y = "OR, excluding AR <= 0.10 cells",
       caption = "Points on the dashed line are unaffected by the artefactual AR ~ 0 cells.") +
  theme_bw() +
  theme(plot.caption = element_text(size = 7, hjust = 0), legend.position = "bottom")
dev.off()

# --- 4. Xist == 0 cross-tabulated against the AR >= 0.90 LOX calls -----------
# Xist is transcribed from the paternal X, which is the Xi in this model. A
# cell that has lost the Xi loses Xist and loses the Xi-derived escapee reads
# at the same time, so Xist == 0 is a marker for Xi loss that does not use the
# allelic data. Caveat: this is a robustness check, not independent evidence -
# it is a coarsened version of the same two variables as the model above.
xist_lox <- xist_ar %>%
  mutate(xist_zero = ifelse(Xist == 0, "Xist == 0", "Xist > 0"),
         LOX_call  = ifelse(allelic_ratio >= 0.90, "LOX-like (AR >= 0.90)", "Other"))

xist_lox_tab <- xist_lox %>%
  count(sample, celltype, xist_zero, LOX_call, name = "n_cells")

write.table(xist_lox_tab, 'Allelic_ratio_results/core_escape_block_new_Xist_zero_vs_LOX_crosstab.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

# Fisher test per sample x celltype: are LOX-like cells enriched among the
# Xist == 0 cells?
xist_lox_test <- xist_lox %>%
  split(list(.$celltype, .$sample), drop = TRUE) %>%
  lapply(function(x) {
    tb <- table(factor(x$xist_zero, levels = c("Xist > 0", "Xist == 0")),
                factor(x$LOX_call,  levels = c("Other", "LOX-like (AR >= 0.90)")))
    if (any(rowSums(tb) == 0) || any(colSums(tb) == 0)) return(NULL)
    ft <- fisher.test(tb)
    data.frame(sample = unique(x$sample), celltype = unique(x$celltype),
               n_cells = nrow(x),
               n_xist_zero = sum(x$xist_zero == "Xist == 0"),
               pct_LOX_in_xist_zero = 100 * tb[2, 2] / sum(tb[2, ]),
               pct_LOX_in_xist_pos  = 100 * tb[1, 2] / sum(tb[1, ]),
               OR = unname(ft$estimate), p_value = ft$p.value)
  }) %>% bind_rows()
rownames(xist_lox_test) <- NULL
if (nrow(xist_lox_test)) xist_lox_test$FDR <- p.adjust(xist_lox_test$p_value, method = "fdr")

write.table(xist_lox_test, 'Allelic_ratio_results/core_escape_block_new_Xist_zero_vs_LOX_fisher.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

pdf('Allelic_ratio_results/core_escape_block_new_Xist_zero_vs_LOX.pdf', width = 10, height = 6)
xist_lox %>%
  count(sample, celltype, xist_zero, LOX_call, name = "n_cells") %>%
  group_by(sample, celltype, xist_zero) %>%
  mutate(prop = n_cells / sum(n_cells), n_group = sum(n_cells)) %>%
  ungroup() %>%
  filter(n_group >= 10) %>%
  ggplot(aes(x = xist_zero, y = prop, fill = LOX_call)) +
  geom_col() +
  geom_text(aes(label = ifelse(LOX_call == "Other", paste0("n=", n_group), "")),
            y = 1.04, size = 2.2, colour = "grey30") +
  facet_grid(celltype ~ sample, labeller = labeller(celltype = label_wrap_gen(18))) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
  labs(x = NULL, y = "Proportion of cells", fill = NULL,
       caption = "Xist is transcribed from the Xi, so Xist == 0 marks candidate Xi-loss cells without using the allelic data.\nGroups with < 10 cells omitted.") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(size = 7),
        legend.position = "bottom",
        plot.caption = element_text(size = 7, hjust = 0))
dev.off()


# ===========================================================================
# Multi-panel UMAP figure, laid out after the senescence/Xist figure style:
#   a  UMAP coloured by allelic ratio, per sample
#   b  UMAP coloured by Xist expression, per sample
#   c  bivariate UMAP, AR x Xist, with a 2x2 colour key
#   d  Xist expression in LOX-like (high AR) vs other cells, per sample
#
# The red quadrant in (c) is the one that matters: high AR with low Xist is
# the candidate Xi-loss population, since Xist is transcribed from the Xi and
# would be lost along with the Xi-derived escapee reads.
# ===========================================================================

umap_name <- grep("umap", names(subset_heart_ceb_flt@reductions), value = TRUE, ignore.case = TRUE)[1]
stopifnot(!is.na(umap_name))

umap_emb <- Embeddings(subset_heart_ceb_flt, reduction = umap_name)[, 1:2]
colnames(umap_emb) <- c("UMAP_1", "UMAP_2")

umap_df <- data.frame(
  umap_emb,
  Xist          = GetAssayData(subset_heart_ceb_flt, assay = "SCT", layer = "data")["Xist", ],
  allelic_ratio = subset_heart_ceb_flt$allelic_ratio,
  total_reads   = subset_heart_ceb_flt$total_reads,
  sample        = subset_heart_ceb_flt$sample,
  celltype      = subset_heart_ceb_flt$celltype
)

# Thresholds for the bivariate classes. AR uses the same 0.9 cut as LOX_status
# above so the figure and the LOX calls agree. Xist is a median split, which is
# dataset-relative by design - the value is printed into the axis label so the
# figure stays self-documenting.
ar_hi_cut   <- 0.90
xist_hi_cut <- median(umap_df$Xist, na.rm = TRUE)

umap_df <- umap_df %>%
  mutate(ar_class   = ifelse(allelic_ratio >= ar_hi_cut, "hi", "lo"),
         xist_class = ifelse(Xist >= xist_hi_cut, "hi", "lo"),
         bivar = factor(paste0("AR ", ar_class, " / Xist ", xist_class),
                        levels = c("AR lo / Xist lo", "AR hi / Xist lo",
                                   "AR lo / Xist hi", "AR hi / Xist hi")))

bivar_cols <- c("AR lo / Xist lo" = "#C9C9C9",   # neither
                "AR hi / Xist lo" = "#E0342A",   # candidate Xi loss
                "AR lo / Xist hi" = "#5AA14F",   # escaping, Xi present
                "AR hi / Xist hi" = "#F5A623")   # high AR despite Xist

umap_theme <- theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        legend.key.height = unit(0.8, "cm"))

# --- a: allelic ratio -------------------------------------------------------
p_a <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, colour = allelic_ratio)) +
  geom_point(size = 0.5, alpha = 0.9) +
  facet_wrap(~sample, nrow = 1) +
  scale_color_gradientn(colors = my_colors, breaks = seq(0, 1, by = 0.25),
                        limits = c(0, 1), oob = scales::squish,
                        name = "Allelic\nratio") +
  labs(x = "UMAP 1", y = "UMAP 2", tag = "a") +
  umap_theme

# --- b: Xist expression -----------------------------------------------------
p_b <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, colour = Xist)) +
  geom_point(size = 0.5, alpha = 0.9) +
  facet_wrap(~sample, nrow = 1) +
  scale_color_gradient(low = "#D9D9D9", high = "#3F7F32", name = "Xist\nexpression") +
  labs(x = "UMAP 1", y = "UMAP 2", tag = "b") +
  umap_theme

# --- c: bivariate AR x Xist -------------------------------------------------
# Points are drawn grey-first so the informative classes sit on top rather than
# being buried under the bulk of the cells.
umap_df_ordered <- umap_df %>% arrange(desc(bivar == "AR lo / Xist lo"))

p_c <- ggplot(umap_df_ordered, aes(x = UMAP_1, y = UMAP_2, colour = bivar)) +
  geom_point(size = 0.5, alpha = 0.9) +
  facet_wrap(~sample, nrow = 1) +
  scale_colour_manual(values = bivar_cols, guide = "none") +
  labs(x = "UMAP 1", y = "UMAP 2", tag = "c") +
  umap_theme

# 2x2 colour key, drawn as its own plot and placed to the right of panel c.
# The lo/hi levels are set explicitly: a discrete scale would otherwise order
# them alphabetically and flip both axes relative to the low-left / high-right
# convention.
key_df <- expand.grid(ar = c("lo", "hi"), xist = c("lo", "hi"), stringsAsFactors = FALSE)
key_df$bivar <- paste0("AR ", key_df$ar, " / Xist ", key_df$xist)
key_df$ar    <- factor(key_df$ar,   levels = c("lo", "hi"))
key_df$xist  <- factor(key_df$xist, levels = c("lo", "hi"))

p_key <- ggplot(key_df, aes(x = ar, y = xist, fill = bivar)) +
  geom_tile(colour = "white", linewidth = 1) +
  scale_fill_manual(values = bivar_cols, guide = "none") +
  scale_x_discrete(labels = c(lo = paste0("< ", ar_hi_cut), hi = paste0(">= ", ar_hi_cut))) +
  scale_y_discrete(labels = c(lo = paste0("< ", round(xist_hi_cut, 2)),
                              hi = paste0(">= ", round(xist_hi_cut, 2)))) +
  labs(x = "Allelic ratio", y = "Xist") +
  coord_fixed() +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7))

# Key sits in its own column on the right, lining up with the colour bars of
# panels a and b. key_width is the knob if it does not line up exactly - it is
# a fraction of the UMAP row width, and the a/b legend columns are roughly 12%.
key_width <- 0.12
p_c_keyed <- (p_c | p_key) + plot_layout(widths = c(1, key_width))

# --- d: Xist in LOX-like vs other cells -------------------------------------
lox_box_df <- umap_df %>%
  mutate(LOX_call = factor(ifelse(allelic_ratio >= ar_hi_cut, "LOX-like", "Other"),
                           levels = c("Other", "LOX-like")))

lox_box_test <- lox_box_df %>%
  split(.$sample) %>%
  lapply(function(x) {
    if (length(unique(x$LOX_call)) < 2) return(NULL)
    tt <- wilcox.test(Xist ~ LOX_call, data = x)
    data.frame(sample = unique(x$sample), p_value = tt$p.value,
               median_other = median(x$Xist[x$LOX_call == "Other"]),
               median_lox   = median(x$Xist[x$LOX_call == "LOX-like"]),
               n_other = sum(x$LOX_call == "Other"), n_lox = sum(x$LOX_call == "LOX-like"))
  }) %>% bind_rows()
if (nrow(lox_box_test)) {
  lox_box_test$FDR   <- p.adjust(lox_box_test$p_value, method = "fdr")
  lox_box_test$label <- paste0("p=", format.pval(lox_box_test$p_value, digits = 2))
}

write.table(lox_box_test, 'Allelic_ratio_results/core_escape_block_new_Xist_by_LOX_call_test.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

p_d <- ggplot(lox_box_df, aes(x = LOX_call, y = Xist, fill = LOX_call)) +
  geom_jitter(width = 0.25, size = 0.3, alpha = 0.25, colour = "grey30") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55) +
  geom_text(data = lox_box_test, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.4, size = 3) +
  facet_wrap(~sample, nrow = 1) +
  scale_fill_manual(values = c(Other = "#C9C9C9", `LOX-like` = "#E0342A"), guide = "none") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Xist expression\n(SCT, normalized)", tag = "d") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey95"),
        plot.margin = margin(14, 5.5, 5.5, 5.5))

pdf('Allelic_ratio_results/core_escape_block_new_AR_Xist_umap_panel.pdf', width = 13, height = 15)
# Tags are set per panel with labs(tag=) rather than tag_levels, so the colour
# key is not counted as a panel and lettering stays a-d.
(p_a / p_b / p_c_keyed / p_d) +
  plot_annotation(caption = paste0(
                    "a) Allelic ratio on UMAP. b) Xist expression on UMAP. ",
                    "c) Bivariate classes; red = high AR with low Xist, the candidate Xi-loss population.\n",
                    "d) Xist expression in LOX-like (AR >= ", ar_hi_cut, ") vs other cells, Wilcoxon rank-sum. ",
                    "Xist split at the median (", round(xist_hi_cut, 2), ")."),
                  theme = theme(plot.caption = element_text(size = 8, hjust = 0)))
dev.off()

# Same figure restricted to Ventricular Cardiomyocytes, where the Xist-AR
# association is significant in all four samples.
umap_df_vcm <- filter(umap_df_ordered, celltype == "Ventricular Cardiomyocytes")
lox_box_vcm <- filter(lox_box_df, celltype == "Ventricular Cardiomyocytes")

lox_box_test_vcm <- lox_box_vcm %>%
  split(.$sample) %>%
  lapply(function(x) {
    if (length(unique(x$LOX_call)) < 2) return(NULL)
    tt <- wilcox.test(Xist ~ LOX_call, data = x)
    data.frame(sample = unique(x$sample), p_value = tt$p.value,
               label = paste0("p=", format.pval(tt$p.value, digits = 2)))
  }) %>% bind_rows()

write.table(lox_box_test_vcm, 'Allelic_ratio_results/core_escape_block_new_Xist_by_LOX_call_test_VCM.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

p_c_vcm <- ggplot(umap_df_vcm, aes(x = UMAP_1, y = UMAP_2, colour = bivar)) +
  geom_point(size = 0.7, alpha = 0.9) +
  facet_wrap(~sample, nrow = 1) +
  scale_colour_manual(values = bivar_cols, guide = "none") +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Ventricular Cardiomyocytes", tag = "a") +
  umap_theme

p_d_vcm <- ggplot(lox_box_vcm, aes(x = LOX_call, y = Xist, fill = LOX_call)) +
  geom_jitter(width = 0.25, size = 0.4, alpha = 0.3, colour = "grey30") +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55) +
  geom_text(data = lox_box_test_vcm, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.4, size = 3) +
  facet_wrap(~sample, nrow = 1) +
  scale_fill_manual(values = c(Other = "#C9C9C9", `LOX-like` = "#E0342A"), guide = "none") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Xist expression\n(SCT, normalized)", tag = "b") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey95"),
        plot.margin = margin(14, 5.5, 5.5, 5.5))

pdf('Allelic_ratio_results/core_escape_block_new_AR_Xist_umap_panel_VCM.pdf', width = 13, height = 8)
(((p_c_vcm | p_key) + plot_layout(widths = c(1, key_width))) / p_d_vcm)
dev.off()
