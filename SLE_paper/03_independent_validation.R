# 03_independent_validation.R
#
# Analyses classification metrics from the independent validation cohort
# (GSE135779) predicted by the ensemble models.
#
# Inputs:
#   - results_update/metrics_*.csv  (from predict_independent_SLE.py)
#
# Outputs:
#   - results_update/predicting_metrics.txt
#   - results_update/predicting_boxplot.pdf

library(dplyr)
library(ggplot2)
library(stringr)

gene.set.colours <- c('chrX' = '#8A0798', 'SLE' = '#D90750')

# ============================================================
# Read in and format prediction metrics
# ============================================================

metrics.files <- list.files('results_update', pattern = 'csv', full.names = TRUE)

metrics <- lapply(metrics.files, read.csv)
names(metrics) <- gsub('.csv', '', basename(metrics.files))

metrics_df <- bind_rows(metrics, .id = 'filename')

metrics_df$split    <- as.numeric(gsub('.*_split_(.*)', '\\1', metrics_df$filename))
metrics_df$age      <- str_extract(metrics_df$filename, 'adult|child')
metrics_df$age      <- factor(metrics_df$age, levels = c('adult', 'child'))
metrics_df$geneset  <- str_extract(metrics_df$filename, 'chrX|SLE')
metrics_df$geneset  <- factor(metrics_df$geneset, levels = c('chrX', 'SLE'))
metrics_df$celltype <- gsub('^metrics_|\\..*', '', metrics_df$filename)
metrics_df$celltype <- gsub('_', ' ', metrics_df$celltype)

write.table(metrics_df, 'results_update/predicting_metrics.txt', row.names = FALSE)
metrics_df <- read.table('results_update/predicting_metrics.txt', header = TRUE)

# ============================================================
# Plot MCC boxplots by cell type, gene set, and age group
# ============================================================

pdf('results_update/predicting_boxplot.pdf')
ggplot(metrics_df, aes(x = MCC, y = celltype, colour = geneset)) +
  geom_boxplot() +
  scale_colour_manual(values = c('#8A0798', '#D90750')) +
  geom_vline(xintercept = 0.7, linetype = 'dashed', colour = 'grey') +
  facet_wrap(~age) +
  labs(title = '') +
  theme_minimal() +
  theme(legend.position = 'right')
dev.off()

# ============================================================
# Statistical tests: chrX vs SLE within each age group
# ============================================================

# --- Adult ---
adult.test <- lapply(unique(metrics_df$celltype), function(x) {
  tmp <- subset(metrics_df, age == 'adult' & celltype == x)
  wilcox.test(MCC ~ geneset, tmp)
})
names(adult.test) <- unique(metrics_df$celltype)
adult.test <- do.call(rbind, lapply(adult.test, function(x) {
  data.frame(statistic = x$statistic, p.value = x$p.value)
}))
adult.test$FDR <- p.adjust(adult.test$p.value, method = 'fdr')
subset(adult.test, FDR > 0.05)

# --- Child ---
child.test <- lapply(unique(metrics_df$celltype), function(x) {
  tmp <- subset(metrics_df, age == 'child' & celltype == x)
  wilcox.test(MCC ~ geneset, alternative = 'greater', tmp)
})
names(child.test) <- unique(metrics_df$celltype)
child.test <- do.call(rbind, lapply(child.test, function(x) {
  data.frame(statistic = x$statistic, p.value = x$p.value)
}))
child.test$FDR <- p.adjust(child.test$p.value, method = 'fdr')
subset(child.test, FDR < 0.05)

# ============================================================
# Mean MCC with 95 % CI per cell type / gene set / age
# ============================================================

mean_metrics <- metrics_df %>%
  group_by(celltype, geneset, age) %>%
  summarise(mean_MCC  = round(mean(MCC), 2),
            lower_CI  = round(quantile(MCC, 0.025), 2),
            upper_CI  = round(quantile(MCC, 0.975), 2)) %>%
  ungroup() %>%
  data.frame()

subset(mean_metrics, age == 'adult' & mean_MCC > 0.7)
