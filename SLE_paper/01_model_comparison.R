# 01_model_comparison.R
#
# Compares individual ML model and ensemble performance across
# feature selection methods, gene sets, and cross-validation splits.
#
# Inputs:
#   - split_{1..10}/{method}/metrics/  (individual model metric CSVs)
#   - split_{1..10}/combined/ensemble/ (ensemble metric CSVs)
#
# Outputs:
#   - ../results_update/models_metrics_df.RData
#   - ../results_update/ensemble_metrics.csv
#   - ../results_update/Supplementary_Table_1.csv
#   - ../results_update/average_model_metrics.csv
#   - ../results_update/average_ensemble_metrics.csv
#   - ../results_update/top_celltypes.RData
#   - ../results_update/top_celltypes_average.csv
#   - Multiple comparison PDFs in ../results_update/

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ComplexHeatmap)
library(circlize)

# ---------- Shared colour palettes ----------
gene.set.colours <- c('chrX' = '#8A0798', 'SLE' = '#D90750')
method.colours   <- c('boruta' = '#BFFB00', 'enet' = '#B875B1',
                       'intersection' = '#D2EDF6', 'combined' = '#4DB748')
model.colours    <- c('logit' = '#F4CE03', 'RF' = '#BCEA9D',
                       'SVM' = '#99B2F5', 'GBM' = '#F5B29E',
                       'MLP' = '#26779E', 'ensemble' = '#F5A2F5')

methods <- c('boruta', 'enet', 'intersection', 'combined')
models  <- c('logit', 'RF', 'SVM', 'GBM', 'MLP')

# ---------- Shared helper: format cell type names ----------
format_celltypes <- function(celltypes) {
  formatted_celltypes <- gsub("_", " ", celltypes)
  formatted_celltypes <- gsub("Non classical", "Non-classical", formatted_celltypes)
  formatted_celltypes <- gsub("alpha beta", "alpha-beta", formatted_celltypes)
  formatted_celltypes <- gsub("Effector Memory", "Effector-Memory", formatted_celltypes)
  formatted_celltypes <- gsub("Cytotoxic GZMH ", "Cytotoxic-GZMH", formatted_celltypes)
  formatted_celltypes <- gsub("Cytotoxic GZMK ", "Cytotoxic-GZMK", formatted_celltypes)
  formatted_celltypes <- gsub("natural killer cell", "Natural Killer cell", formatted_celltypes)
  formatted_celltypes <- gsub("CD4-positive", "CD4 positive", formatted_celltypes)
  formatted_celltypes <- gsub("CD8-positive", "CD8 positive", formatted_celltypes)
  return(formatted_celltypes)
}

# ---------- Shared helper: 95 % CI ----------
calc_CI <- function(x) {
  mean_score    <- mean(x)
  std_dev_score <- sd(x)
  n <- length(x)
  alpha <- 0.05
  t_critical   <- qt(alpha / 2, df = n - 1, lower.tail = FALSE)
  margin_error  <- t_critical * (std_dev_score / sqrt(n))
  lower_bound   <- mean_score - margin_error
  upper_bound   <- mean_score + margin_error
  return(c(lower_bound, upper_bound))
}

# ============================================================
# SECTION 1 — Individual model metrics across splits
# ============================================================

models_metrics_list <- list()
for (i in 1:10) {
  metric.files <- unlist(lapply(methods, function(method) {
    list.files(paste0('split_', i, '/', method, '/metrics'),
               pattern = 'metrics_', full.names = TRUE)
  }))
  metric.files <- metric.files[file.size(metric.files) > 0]
  metrics <- lapply(metric.files, read.csv)
  names(metrics) <- gsub('split_\\d+/|metrics/|.csv', '', metric.files) %>%
    gsub('_metrics', '', .) %>%
    gsub('/', '_', .)
  metrics_df <- bind_rows(metrics, .id = 'celltype') %>%
    mutate(
      model    = str_extract(celltype, "logit|RF|SVM|GBM|MLP"),
      gene.set = str_extract(celltype, "HVG\\.autosome|chrX|autosome|HVG|SLE"),
      method   = str_extract(celltype, "boruta|enet|intersection|combined"),
      celltype = gsub("^(boruta_|enet_|intersection_|combined_)?(logit_|RF_|SVM_|GBM_|MLP_)", "", celltype),
      celltype = gsub("\\.HVG\\.autosome|\\.SLE|\\.chrX|\\.autosome|\\.HVG", "", celltype)
    )
  metrics_df <- subset(metrics_df, gene.set %in% c('chrX', 'SLE'))
  models_metrics_list[[i]] <- metrics_df
}
names(models_metrics_list) <- paste('split', 1:10, sep = '_')

# --- Test: method effect on MCC across splits ---
compare_method <- lapply(models_metrics_list,
                         function(x) kruskal.test(MCC ~ method, data = x)$p.value)
compare_method <- t(bind_rows(compare_method, .id = 'split'))
compare_method <- data.frame(p.value = compare_method[, 1],
                              FDR = p.adjust(compare_method[, 1], method = 'fdr'))

# Combine all splits
models_metrics_df <- bind_rows(models_metrics_list, .id = 'split')
models_metrics_df$model    <- factor(models_metrics_df$model,
                                     levels = c('logit', 'RF', 'SVM', 'GBM', 'MLP'))
models_metrics_df$gene.set <- factor(models_metrics_df$gene.set, levels = c('chrX', 'SLE'))
models_metrics_df$split    <- factor(models_metrics_df$split,
                                     levels = paste('split', 1:10, sep = '_'))

save(models_metrics_df, file = '../results_update/models_metrics_df.RData')

# Feature counts per gene set
lapply(split(models_metrics_df, models_metrics_df$gene.set), function(x) {
  c(mean = round(mean(x$n_features), 2), sd = round(sd(x$n_features), 2))
})

# --- Plot: method comparison across splits ---
pdf('../results_update/compare_method_across_splits.pdf')
ggplot(models_metrics_df, aes(x = method, y = MCC, colour = method, group = method)) +
  geom_jitter(width = 0.2, alpha = 1) +
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  theme_minimal() +
  labs(x = 'Method', y = 'MCC') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = method.colours, name = 'Method') +
  facet_wrap(~split, ncol = 5, nrow = 2)
dev.off()

# Pairwise method comparisons
lapply(models_metrics_list, function(x) {
  pairwise.wilcox.test(x$MCC, x$method, p.adjust.method = 'fdr', alternative = 'greater')
})

# Filter for 'combined' method for subsequent model/gene set tests
models_metrics_list <- lapply(models_metrics_list,
                               function(x) subset(x, method == 'combined'))

# --- Test: gene set effect on MCC ---
compare_gene.set <- lapply(models_metrics_list,
                            function(x) kruskal.test(MCC ~ gene.set, data = x)$p.value)
compare_gene.set <- t(bind_rows(compare_gene.set, .id = 'split'))
compare_gene.set <- data.frame(p.value = compare_gene.set[, 1],
                                FDR = p.adjust(compare_gene.set[, 1], method = 'fdr'))

pdf('../results_update/compare_geneset_across_splits.pdf')
comparisons <- list(c('chrX', 'SLE'))
ggplot(models_metrics_df, aes(x = gene.set, y = MCC, colour = gene.set, group = gene.set)) +
  geom_boxplot() +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE,
              test = 'wilcox.test', color = 'black', step_increase = 0.1) +
  theme_minimal() +
  labs(x = 'Gene Set', y = 'MCC') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = gene.set.colours, name = 'Gene Set') +
  facet_wrap(~split, ncol = 5, nrow = 2)
dev.off()

# --- Test: model effect on MCC ---
compare_model <- lapply(models_metrics_list,
                         function(x) kruskal.test(MCC ~ model, data = x)$p.value)
compare_model <- t(bind_rows(compare_model, .id = 'split'))
compare_model <- data.frame(p.value = compare_model[, 1],
                              FDR = p.adjust(compare_model[, 1], method = 'fdr'))

lapply(models_metrics_list, function(x) {
  pairwise.wilcox.test(x$MCC, x$model, p.adjust.method = 'fdr', alternative = 'greater')
})

pdf('../results_update/compare_model_across_splits.pdf')
ggplot(models_metrics_df, aes(x = model, y = MCC, colour = model, group = model)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = 'Model', y = 'MCC') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = model.colours, name = 'Model') +
  facet_wrap(~split, ncol = 5, nrow = 2)
dev.off()

# --- Average individual model metrics across splits ---
average_model_metrics_df <- models_metrics_df %>%
  group_by(celltype, gene.set, model) %>%
  summarise(MCC = mean(MCC),
            n_features = round(mean(n_features), 0)) %>%
  data.frame()
average_model_metrics_df$model    <- factor(average_model_metrics_df$model,
                                            levels = c('logit', 'RF', 'SVM', 'GBM', 'MLP'))
average_model_metrics_df$gene.set <- factor(average_model_metrics_df$gene.set,
                                            levels = c('chrX', 'SLE'))

comparisons <- list(c('chrX', 'SLE'))
pdf('../results_update/avg_model_MCC_geneset.pdf')
ggplot(average_model_metrics_df, aes(x = gene.set, y = MCC, colour = gene.set)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE,
              test = 'wilcox.test', color = 'black', step_increase = 0.1) +
  theme_minimal() +
  labs(x = 'Gene Set', y = 'MCC') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = gene.set.colours, name = 'Gene Set')
dev.off()

# Pairwise model comparisons
kruskal.test(MCC ~ model, data = average_model_metrics_df)
MCC_model <- pairwise.wilcox.test(average_model_metrics_df$MCC,
                                   average_model_metrics_df$model,
                                   p.adjust.method = 'fdr')
MCC_model$p.value[is.na(MCC_model$p.value)] <- 1

pdf('../results_update/pairwise_wilcox_model.heatmap.pdf')
Heatmap(MCC_model$p.value,
        col = circlize::colorRamp2(c(1, 0.05), c("blue", "red")),
        name = 'FDR')
dev.off()

comparisons <- list(c('logit', 'MLP'), c('RF', 'MLP'),
                    c('SVM', 'MLP'), c('GBM', 'MLP'))
pdf('../results_update/avg_model_MCC_model.pdf')
ggplot(average_model_metrics_df, aes(x = model, y = MCC, colour = model)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE,
              test = 'wilcox.test', color = 'black', step_increase = 0.1) +
  theme_minimal() +
  labs(x = 'Model', y = 'MCC') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = model.colours, name = 'Model')
dev.off()

# ============================================================
# SECTION 2 — Ensemble model metrics
# ============================================================

ensemble_metrics_list <- list()
for (i in 1:10) {
  metric.files <- list.files(paste0('split_', i, '/combined/ensemble'),
                              pattern = 'metrics_', full.names = TRUE)
  metrics <- lapply(metric.files, read.csv)
  names(metrics) <- gsub('split_\\d+/|ensemble/metrics_|.csv', '', metric.files) %>%
    gsub('/', '_', .)
  metrics_df <- bind_rows(metrics, .id = 'celltype') %>%
    mutate(
      gene.set = str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
      celltype = gsub('combined_|.HVG.autosome', '', celltype),
      celltype = gsub('.SLE|.chrX|.autosome|.HVG', '', celltype)
    )
  metrics_df <- subset(metrics_df, gene.set %in% c('chrX', 'SLE'))
  ensemble_metrics_list[[i]] <- metrics_df
}
names(ensemble_metrics_list) <- paste('split', 1:10, sep = '_')

# Test: gene set effect on ensemble MCC
compare_gene.set <- lapply(ensemble_metrics_list,
                            function(x) kruskal.test(MCC ~ gene.set, data = x)$p.value)
compare_gene.set <- t(bind_rows(compare_gene.set, .id = 'split'))
compare_gene.set <- data.frame(p.value = compare_gene.set[, 1],
                                FDR = p.adjust(compare_gene.set[, 1], method = 'fdr'))

# Combine ensemble metrics
ensemble_metrics_df <- bind_rows(ensemble_metrics_list, .id = 'split')
ensemble_metrics_df$gene.set <- factor(ensemble_metrics_df$gene.set, levels = c('chrX', 'SLE'))
ensemble_metrics_df$split    <- factor(ensemble_metrics_df$split,
                                        levels = paste('split', 1:10, sep = '_'))
ensemble_metrics_df$method   <- 'ensemble'

# Format cell type names
average_model_metrics_df$celltype <- format_celltypes(average_model_metrics_df$celltype)
ensemble_metrics_df$celltype      <- format_celltypes(ensemble_metrics_df$celltype)

write.csv(ensemble_metrics_df, '../results_update/ensemble_metrics.csv', row.names = FALSE)
ensemble_metrics_df <- read.csv('../results_update/ensemble_metrics.csv')
ensemble_metrics_df$model <- 'ensemble'

# --- Combined all-models comparison ---
all_models <- rbind(
  models_metrics_df[, c('split', 'celltype', 'gene.set', 'model', 'MCC')],
  ensemble_metrics_df[, c('split', 'celltype', 'gene.set', 'model', 'MCC')]
)
all_models$model <- factor(all_models$model,
                           levels = c('logit', 'RF', 'SVM', 'GBM', 'MLP', 'ensemble'))
pairwise.wilcox.test(all_models$MCC, all_models$model, p.adjust.method = 'fdr')

pdf('../results_update/all_models_boxplot.pdf')
ggplot(all_models, aes(x = model, y = MCC, colour = model)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  theme_minimal() +
  labs(x = 'Model', y = 'MCC') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = model.colours, name = 'Model')
dev.off()

range(subset(ensemble_metrics_df, celltype == 'Progenitor cell' & gene.set == 'chrX')$MCC)

write.csv(average_model_metrics_df, '../results_update/average_model_metrics.csv')

# Supplementary Table 1
models_metrics_df <- models_metrics_df[, colnames(ensemble_metrics_df)]
Supplementary_Table_1 <- rbind(models_metrics_df, ensemble_metrics_df)
write.csv(Supplementary_Table_1, '../results_update/Supplementary_Table_1.csv', row.names = FALSE)

# --- Ensemble gene set comparison across splits ---
comparisons <- list(c('chrX', 'SLE'))
pdf('../results_update/ensemble_geneset_across_splits_MCC.pdf')
ggplot(ensemble_metrics_df, aes(x = gene.set, y = MCC, colour = gene.set)) +
  geom_boxplot() +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE,
              test = 'wilcox.test', color = 'black', step_increase = 0.1) +
  theme_minimal() +
  labs(x = 'gene set', y = 'MCC') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = gene.set.colours) +
  facet_wrap(~split, ncol = 5, nrow = 2)
dev.off()

wilcox.test(MCC ~ gene.set, data = ensemble_metrics_df)

pdf('../results_update/ensemble_geneset_MCC.pdf')
comparisons <- list(c('chrX', 'SLE'))
ggplot(ensemble_metrics_df, aes(x = gene.set, y = MCC, colour = gene.set)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, colour = 'black', fill = NA) +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE,
              test = 'wilcox.test', color = 'black', step_increase = 0.1) +
  theme_minimal() +
  labs(x = 'Gene Set', y = 'MCC') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = gene.set.colours, name = 'Gene Set')
dev.off()

# --- Identify top cell types (no significant gene set difference) ---
compare_celltype <- lapply(
  split(ensemble_metrics_df, ensemble_metrics_df$celltype),
  function(x) wilcox.test(MCC ~ gene.set, data = x)$p.value
)
compare_celltype <- t(bind_rows(compare_celltype, .id = 'celltype'))
compare_celltype <- data.frame(p.value = compare_celltype[, 1],
                                FDR = p.adjust(compare_celltype[, 1], method = 'fdr'))
subset(compare_celltype, FDR > 0.05)
top_celltypes <- rownames(subset(compare_celltype, FDR > 0.05))

save(top_celltypes, file = '../results_update/top_celltypes.RData')

p.adjust(unlist(lapply(
  split(ensemble_metrics_df, ensemble_metrics_df$celltype),
  function(x) pairwise.wilcox.test(x$MCC, x$gene.set, p.adjust.method = 'fdr')$p.value
)), method = 'fdr')

# --- Ensemble by cell type ---
comparisons <- list(c('chrX', 'SLE'))
pdf('../results_update/ensemble_geneset_celltype_MCC.pdf', width = 10, height = 15)
ggplot(ensemble_metrics_df, aes(x = gene.set, y = MCC, colour = gene.set)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, colour = 'black', fill = NA) +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE,
              test = 'wilcox.test', color = 'black', step_increase = 0.1) +
  theme_minimal() +
  labs(x = 'Gene Set', y = 'MCC') +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size = 8)) +
  scale_colour_manual(values = gene.set.colours, name = 'Gene Set') +
  facet_wrap(~celltype, ncol = 5, nrow = 5, strip.position = 'bottom')
dev.off()

pdf('../results_update/top_models_geneset.pdf', width = 10, height = 10)
top_models <- subset(ensemble_metrics_df, celltype %in% top_celltypes)
ggplot(top_models, aes(x = gene.set, y = MCC, colour = gene.set)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE,
              test = 'wilcox.test', color = 'black', step_increase = 0.1) +
  theme_minimal() +
  labs(x = '', y = 'MCC') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = gene.set.colours, name = 'Gene Set') +
  facet_wrap(~celltype, strip.position = 'bottom')
dev.off()

sort(table(subset(ensemble_metrics_df,
                  gene.set == 'chrX' & MCC > 0.7 & celltype %in% top_celltypes
                  )[, c('celltype', 'MCC')]$celltype))

# ============================================================
# SECTION 3 — Average ensemble metrics & ensemble vs models
# ============================================================

average_ensemble_metrics <- lapply(
  split(ensemble_metrics_df,
        interaction(ensemble_metrics_df$celltype, ensemble_metrics_df$gene.set, sep = '_')),
  function(x) {
    data.frame(
      F1 = mean(x$F1), F1_lower = calc_CI(x$F1)[1], F1_upper = calc_CI(x$F1)[2],
      AUPRC = mean(x$AUPRC), AUPRC_lower = calc_CI(x$AUPRC)[1],
      AUPRC_upper = calc_CI(x$F1)[2],
      MCC = mean(x$MCC), MCC_lower = calc_CI(x$MCC)[1], MCC_upper = calc_CI(x$MCC)[2],
      n_features = round(mean(x$n_features), 0)
    )
  }
)
average_ensemble_metrics <- bind_rows(average_ensemble_metrics, .id = 'celltype_gene.set') %>%
  separate_wider_delim(celltype_gene.set, delim = "_", names = c('celltype', 'gene.set')) %>%
  data.frame()
average_ensemble_metrics$gene.set <- factor(average_ensemble_metrics$gene.set,
                                            levels = c('chrX', 'SLE'))

write.csv(average_ensemble_metrics, '../results_update/average_ensemble_metrics.csv',
          row.names = FALSE)

pdf('../results_update/avg_ensemble_MCC_geneset.pdf')
comparisons <- list(c('chrX', 'SLE'))
ggplot(average_ensemble_metrics, aes(x = gene.set, y = MCC, colour = gene.set)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, colour = 'black', fill = NA) +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE,
              test = 'wilcox.test', color = 'black', step_increase = 0.1) +
  theme_minimal() +
  labs(x = 'Gene Set', y = 'MCC') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = gene.set.colours, name = 'Gene Set')
dev.off()

pdf('../results_update/avg_ensemble_MCC_forest.pdf', width = 10, height = 5)
ggplot(average_ensemble_metrics, aes(x = MCC, y = celltype)) +
  geom_errorbarh(aes(xmin = MCC_lower, xmax = MCC_upper)) +
  geom_vline(xintercept = 0.7, linetype = 'dotted', color = 'red') +
  geom_point() +
  labs(x = 'MCC', y = '') +
  theme(axis.text.y = element_text(size = 12)) +
  facet_wrap(~gene.set, ncol = 5, nrow = 1)
dev.off()

top_celltypes_average <- subset(average_ensemble_metrics,
                                celltype %in% top_celltypes)[, c('celltype', 'gene.set',
                                                                  'MCC', 'MCC_lower', 'MCC_upper')]
top_celltypes_average$MCC       <- round(top_celltypes_average$MCC, 2)
top_celltypes_average$MCC_lower <- round(top_celltypes_average$MCC_lower, 2)
top_celltypes_average$MCC_upper <- round(top_celltypes_average$MCC_upper, 2)

write.csv(top_celltypes_average, '../results_update/top_celltypes_average.csv')

subset(average_ensemble_metrics,
       gene.set == 'chrX' & MCC_upper > 0.7 & celltype %in% top_celltypes
       )[, c('celltype', 'MCC', 'MCC_lower', 'MCC_upper')]

pdf('../results_update/avg_ensemble_MCC_by_n_features.pdf', width = 10, height = 5)
ggplot(average_ensemble_metrics, aes(x = n_features, y = MCC)) +
  geom_point() +
  theme_minimal() +
  labs(x = 'Number of Features', y = 'MCC') +
  geom_smooth(method = 'lm', se = FALSE) +
  stat_cor(r.digits = 2, cor.coef.name = 'rho', method = "spearman") +
  facet_wrap(~gene.set, ncol = 5, nrow = 1, scales = 'free_x')
dev.off()

comparisons <- list(c('chrX', 'SLE'))
pdf('../results_update/avg_ensemble_n_features_geneset.pdf')
ggplot(average_ensemble_metrics, aes(x = gene.set, y = n_features, colour = gene.set)) +
  geom_boxplot() +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE,
              test = 'wilcox.test', color = 'black', step_increase = 0.2) +
  theme_minimal() +
  labs(x = 'Gene Set', y = 'Number of Features') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = gene.set.colours, name = 'Gene Set')
dev.off()

pairwise.wilcox.test(average_ensemble_metrics$n_features,
                     average_ensemble_metrics$gene.set,
                     p.adjust.method = 'fdr')

# --- Ensemble vs individual models ---
average_ensemble_metrics$model <- 'ensemble'
combined_metrics <- bind_rows(average_model_metrics_df, average_ensemble_metrics)
combined_metrics$model <- factor(combined_metrics$model,
                                 levels = c('logit', 'RF', 'SVM', 'GBM', 'MLP', 'ensemble'))

pairwise.wilcox.test(combined_metrics$MCC, combined_metrics$model, p.adjust.method = 'fdr')

pdf('../results_update/models_ensemble_boxplot.pdf')
combinations <- list(c('logit', 'ensemble'), c('RF', 'ensemble'),
                     c('SVM', 'ensemble'), c('GBM', 'ensemble'),
                     c('MLP', 'ensemble'))
ggplot(combined_metrics, aes(x = model, y = MCC, colour = model)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, colour = 'black', fill = NA) +
  geom_signif(comparisons = combinations, map_signif_level = TRUE,
              test = 'wilcox.test', color = 'black', step_increase = 0.1) +
  theme_minimal() +
  labs(x = 'Model', y = 'MCC') +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = model.colours, name = 'Model')
dev.off()
