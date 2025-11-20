library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(reshape2)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(ggsignif)
library(gridExtra)
library(ggpubr)

gene.set.colours <- c('chrX'='#8A0798', 'SLE'='#D90750')
method.colours <- c('boruta'='#BFFB00', 'enet'='#B875B1', 'intersection'='#D2EDF6', 'combined'='#4DB748')
model.colours <- c('logit'='#F4CE03', 'RF'='#BCEA9D', 'SVM'='#99B2F5', 'GBM'='#F5B29E', 'MLP'='#26779E', 'ensemble'='#F5A2F5')

methods <- c('boruta', 'enet', 'intersection', 'combined')
models <- c('logit', 'RF', 'SVM', 'GBM', 'MLP')

load('../escapees.Rdata')
escape <- rownames(escape)
SLE <- read.delim('../SLE_DisGeNet.tsv')$Gene
chrX <- read.delim('../chrX_biomaRt.txt')$Gene.name
chrX <- unique(chrX[!chrX %in% ''])


# Read in metrics for all ML models across gene sets and splits
models_metrics_list <- list()
for(i in 1:10){
    metric.files <- unlist(lapply(methods, function(method){
    list.files(paste0('split_', i, '/', method, '/metrics'), pattern='metrics_', full.names=TRUE)
    }))
    # Check if metric files are empty
    metric.files <- metric.files[file.size(metric.files) > 0]
    metrics <- lapply(metric.files, read.csv)
    names(metrics) <- gsub('split_\\d+/|metrics/|.csv', '', metric.files) %>% 
        gsub('_metrics', '', .) %>%
        gsub('/', '_', .)
    metrics_df <- bind_rows(metrics, .id = 'celltype') %>%
    mutate(
        model = str_extract(celltype, "logit|RF|SVM|GBM|MLP"),
        gene.set = str_extract(celltype, "HVG\\.autosome|chrX|autosome|HVG|SLE"),
        method = str_extract(celltype, "boruta|enet|intersection|combined"),
        celltype = gsub("^(boruta_|enet_|intersection_|combined_)?(logit_|RF_|SVM_|GBM_|MLP_)", "", celltype),
        celltype = gsub("\\.HVG\\.autosome|\\.SLE|\\.chrX|\\.autosome|\\.HVG", "", celltype)
    )
    # Filter out autosome, HVG, and HVG.autosome gene sets
    metrics_df <- subset(metrics_df, gene.set %in% c('chrX', 'SLE'))
    models_metrics_list[[i]] <- metrics_df
}
names(models_metrics_list) <- paste('split', 1:10, sep='_')

### Testing for difference in MCC score between method across splits
compare_method <- lapply(models_metrics_list, function(x){kruskal.test(MCC ~ method, data=x)$p.value})
compare_method <- t(bind_rows(compare_method, .id='split'))
compare_method <- data.frame(p.value=compare_method[,1], FDR=p.adjust(compare_method[,1], method='fdr'))

# Combine metrics from all splits
models_metrics_df <- bind_rows(models_metrics_list, .id='split')
models_metrics_df$model <- factor(models_metrics_df$model, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP'))
# models_metrics_df$method <- factor(models_metrics_df$method, levels=c('boruta', 'enet', 'intersection', 'combined'))
models_metrics_df$gene.set <- factor(models_metrics_df$gene.set, levels=c('chrX', 'SLE'))
models_metrics_df$split <- factor(models_metrics_df$split, levels=paste('split', 1:10, sep='_'))

save(models_metrics_df, file='../results_update/models_metrics_df.RData')
load('../results_update/models_metrics_df.RData')

lapply(split(models_metrics_df, models_metrics_df$gene.set), function(x){
    # Mean and s.d of n_features
    c(mean=round(mean(x$n_features), 2), sd=round(sd(x$n_features), 2))
})

pdf('../results_update/compare_method_across_splits.pdf')
ggplot(models_metrics_df, aes(x=method, y=MCC, colour=method, group=method)) +
    geom_jitter(width=0.2, alpha=1) +
    geom_boxplot(outlier.shape = NA, color = 'black', fill=NA) +
    theme_minimal() +
    labs(x='Method', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=method.colours, name='Method') +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

# Pairwise comparisons between methods for MCC score
lapply(models_metrics_list, function(x){
    pairwise.wilcox.test(x$MCC, x$method, p.adjust.method='fdr', alternative='greater')
})

# Filter for combined method
models_metrics_list <- lapply(models_metrics_list, function(x) subset(x, method == 'combined'))

# Testing for differences in MCC score between gene sets across splits
compare_gene.set <- lapply(models_metrics_list, function(x){kruskal.test(MCC ~ gene.set, data=x)$p.value})
compare_gene.set <- t(bind_rows(compare_gene.set, .id='split'))
compare_gene.set <- data.frame(p.value=compare_gene.set[,1], FDR=p.adjust(compare_gene.set[,1], method='fdr'))

pdf('../results_update/compare_geneset_across_splits.pdf')
comparisons <- list(c('chrX', 'SLE'))
ggplot(models_metrics_df, aes(x=gene.set, y=MCC, colour=gene.set, group=gene.set)) +
    geom_boxplot() +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, 
            test='wilcox.test', color='black', step_increase=0.1) +
    theme_minimal() +
    labs(x='Gene Set', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set') +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

# Testing for difference in MCC score between models across splits
compare_model <- lapply(models_metrics_list, function(x){kruskal.test(MCC ~ model, data=x)$p.value})
compare_model <- t(bind_rows(compare_model, .id='split'))
compare_model <- data.frame(p.value=compare_model[,1], FDR=p.adjust(compare_model[,1], method='fdr'))

lapply(models_metrics_list, function(x){
    pairwise.wilcox.test(x$MCC, x$model, p.adjust.method='fdr', alternative='greater')
})

pdf('../results_update/compare_model_across_splits.pdf')
ggplot(models_metrics_df, aes(x=model, y=MCC, colour=model, group=model)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x='Model', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=model.colours, name='Model') +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

# Average metrics across splits to get a single value for each model
average_model_metrics_df <- models_metrics_df %>%
    group_by(celltype, gene.set, model) %>%
    summarise(
        MCC=mean(MCC),
        n_features=round(mean(n_features),0)) %>%
    data.frame()
average_model_metrics_df$model <- factor(average_model_metrics_df$model, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP'))
average_model_metrics_df$gene.set <- factor(average_model_metrics_df$gene.set, levels=c('chrX','SLE'))

comparisons <- list(c('chrX', 'SLE'))
pdf('../results_update/avg_model_MCC_geneset.pdf')
ggplot(average_model_metrics_df, aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape=NA, color='black', fill=NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, 
            test='wilcox.test', color='black', step_increase=0.1) +
    theme_minimal() +
    labs(x='Gene Set', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set')
dev.off()

# Pairwise comparisons between MCC and model
kruskal.test(MCC ~ model, data=average_model_metrics_df)
MCC_model <- pairwise.wilcox.test(average_model_metrics_df$MCC, average_model_metrics_df$model, 
    p.adjust.method='fdr')
MCC_model$p.value[is.na(MCC_model$p.value)] <- 1

pdf('../results_update/pairwise_wilcox_model.heatmap.pdf')
Heatmap(MCC_model$p.value, col = circlize::colorRamp2(c(1, 0.05), c("blue", "red")), name = 'FDR')
dev.off()

comparisons <- list(c('logit', 'MLP'), c('RF', 'MLP'), c('SVM', 'MLP'), c('GBM', 'MLP'))
pdf('../results_update/avg_model_MCC_model.pdf')
ggplot(average_model_metrics_df, aes(x=model, y=MCC, colour=model)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, 
    test='wilcox.test', color='black', step_increase=0.1) +
    theme_minimal() +
    labs(x='Model', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=model.colours, name='Model')
dev.off()

### Read in metrics file for ensemble models across gene sets and splits ###
ensemble_metrics_list <- list()
for(i in 1:10){
    metric.files <- list.files(paste0('split_', i, '/combined/ensemble'), pattern='metrics_', full.names=TRUE)
    metrics <- lapply(metric.files, read.csv)
    names(metrics) <- gsub('split_\\d+/|ensemble/metrics_|.csv', '', metric.files) %>% gsub('/', '_', .)
    metrics_df <- bind_rows(metrics, .id='celltype') %>%
        mutate(
        gene.set=str_extract(celltype, "HVG.autosome|chrX|autosome|HVG|SLE"),
        celltype=gsub('combined_|.HVG.autosome', '', celltype),
        celltype=gsub('.SLE|.chrX|.autosome|.HVG', '', celltype))
    # Filter out autosome, HVG, and HVG.autosome gene sets
    metrics_df <- subset(metrics_df, gene.set %in% c('chrX', 'SLE'))
    ensemble_metrics_list[[i]] <- metrics_df
}
names(ensemble_metrics_list) <- paste('split', 1:10, sep='_')

### Testing for difference in MCC score between gene sets across splits
compare_gene.set <- lapply(ensemble_metrics_list, function(x){kruskal.test(MCC ~ gene.set, data=x)$p.value})
compare_gene.set <- t(bind_rows(compare_gene.set, .id='split'))
compare_gene.set <- data.frame(p.value=compare_gene.set[,1], FDR=p.adjust(compare_gene.set[,1], method='fdr'))

# Combine metrics from all splits
ensemble_metrics_df <- bind_rows(ensemble_metrics_list, .id='split')
# Format gene.set as factor
ensemble_metrics_df$gene.set <- factor(ensemble_metrics_df$gene.set, levels=c('chrX','SLE'))
# Format split as factor
ensemble_metrics_df$split <- factor(ensemble_metrics_df$split, levels=paste('split', 1:10, sep='_'))
# Add ensemble as method
ensemble_metrics_df$method <- 'ensemble'

format_celltypes <- function(celltypes) {
  # Replace underscores with spaces
  formatted_celltypes <- gsub("_", " ", celltypes)
  
  # Add dashes where appropriate (e.g., "Non classical" -> "Non-classical")
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

# Format celltype names
average_model_metrics_df$celltype <- format_celltypes(average_model_metrics_df$celltype)
ensemble_metrics_df$celltype <- format_celltypes(ensemble_metrics_df$celltype)

write.csv(ensemble_metrics_df, '../results_update/ensemble_metrics.csv', row.names=FALSE)
ensemble_metrics_df <- read.csv('../results_update/ensemble_metrics.csv')
ensemble_metrics_df$model <- 'ensemble'

# Add ensemble model to model_metrics_df
all_models <- rbind(models_metrics_df[,c('split', 'celltype', 'gene.set', 'model', 'MCC')], 
ensemble_metrics_df[,c('split', 'celltype', 'gene.set', 'model', 'MCC')])
all_models$model <- factor(all_models$model, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP', 'ensemble'))
pairwise.wilcox.test(all_models$MCC, all_models$model, p.adjust.method='fdr')

pdf('../results_update/all_models_boxplot.pdf')
ggplot(all_models, aes(x=model, y=MCC, colour=model)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape = NA, color = 'black', fill = NA) +
    theme_minimal() +
    labs(x='Model', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=model.colours, name='Model')
dev.off()

range(subset(ensemble_metrics_df, celltype == 'Progenitor cell' & gene.set == 'chrX')$MCC)

write.csv(average_model_metrics_df,'../results_update/average_model_metrics.csv')
#average_model_metrics_df <- read.csv('../results_update/average_model_metrics.csv')

# order colnames(ensemble_metrics_df) by colnames(models_metrics_df)
models_metrics_df <- models_metrics_df[,colnames(ensemble_metrics_df)]
Supplementary_Table_1 <- rbind(models_metrics_df, ensemble_metrics_df)
write.csv(Supplementary_Table_1, '../results_update/Supplementary_Table_1.csv', row.names=FALSE)

comparisons <- list(c('chrX', 'SLE'))
pdf('../results_update/ensemble_geneset_across_splits_MCC.pdf')
ggplot(ensemble_metrics_df, aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_boxplot() +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, 
    test='wilcox.test', color='black', step_increase=0.1) +
    theme_minimal() +
    labs(x='gene set', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours) +
    facet_wrap(~split, ncol=5, nrow=2)
dev.off()

wilcox.test(MCC ~ gene.set, data=ensemble_metrics_df)

pdf('../results_update/ensemble_geneset_MCC.pdf')
comparisons <- list(c('chrX', 'SLE'))
ggplot(ensemble_metrics_df, aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape=NA, colour='black', fill=NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, test='wilcox.test', 
    color='black', step_increase=0.1) +
    theme_minimal() +
    labs(x='Gene Set', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set')
dev.off()

compare_celltype <- lapply(split(ensemble_metrics_df, ensemble_metrics_df$celltype), function(x){
    wilcox.test(MCC ~ gene.set, data=x)$p.value
})
compare_celltype <- t(bind_rows(compare_celltype, .id='celltype'))
compare_celltype <- data.frame(p.value=compare_celltype[,1], FDR=p.adjust(compare_celltype[,1], method='fdr'))
subset(compare_celltype, FDR > 0.05)
top_celltypes <- rownames(subset(compare_celltype, FDR > 0.05))

tmp <- 'CD8 positive, alpha beta T cell'
# Calculate the 95% CI
calc_CI <- function(x){
    mean_score <- mean(x)
    std_dev_score <- sd(x)
    n <- length(x)
    alpha <- 0.05
    # Calculate the t critical value
    t_critical <- qt(alpha / 2, df = n - 1, lower.tail = FALSE)
    # Calculate the margin of error
    margin_error <- t_critical * (std_dev_score / sqrt(n))
    # Calculate the lower and upper bounds of the 95% CI
    lower_bound <- mean_score - margin_error
    upper_bound <- mean_score + margin_error
    return(c(round(lower_bound, 2), round(upper_bound, 2)))
}
calc_CI(subset(metrics_df, celltype == tmp & geneset == 'chrX')$MCC)

save(top_celltypes, file = '../results_update/top_celltypes.RData')
load('../results_update/top_celltypes.RData')

p.adjust(unlist(lapply(split(ensemble_metrics_df, ensemble_metrics_df$celltype), function(x){
    pairwise.wilcox.test(x$MCC, x$gene.set, p.adjust.method='fdr')$p.value
})), method='fdr')

comparisons <- list(c('chrX', 'SLE'))
pdf('../results_update/ensemble_geneset_celltype_MCC.pdf', width=10, height=15)
ggplot(ensemble_metrics_df, aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape=NA, colour='black', fill=NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, 
    test='wilcox.test', color='black', step_increase=0.1) +
    theme_minimal() +
    labs(x='Gene Set', y='MCC') +
    theme(
        axis.text.x=element_blank(),
        strip.text = element_text(size = 8)) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set') +
    facet_wrap(~celltype, ncol=5, nrow=5, strip.position = 'bottom')
dev.off()

pdf('../results_update/top_models_geneset.pdf', width=10, height=10)
top_models <- subset(ensemble_metrics_df, celltype %in% top_celltypes)
ggplot(top_models, aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape = NA, color='black', fill=NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, 
    test='wilcox.test', color='black', step_increase=0.1) +
    theme_minimal() +
    labs(x='', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set') +
    facet_wrap(~celltype, strip.position = 'bottom')
dev.off()

sort(table(subset(ensemble_metrics_df, gene.set == 'chrX' & MCC > 0.7 & celltype %in% top_celltypes)[,c('celltype', 'MCC')]$celltype))

calc_CI <- function(x){
    mean_score <- mean(x)
    std_dev_score <- sd(x)
    n <- length(x)
    alpha <- 0.05
    # Calculate the t critical value
    t_critical <- qt(alpha / 2, df = n - 1, lower.tail = FALSE)
    # Calculate the margin of error
    margin_error <- t_critical * (std_dev_score / sqrt(n))
    # Calculate the lower and upper bounds of the 95% CI
    lower_bound <- mean_score - margin_error
    upper_bound <- mean_score + margin_error
    return(c(lower_bound, upper_bound))
}

### Average metrics across splits to get a single value for each model ###
average_ensemble_metrics <- lapply(split(ensemble_metrics_df, interaction(ensemble_metrics_df$celltype, ensemble_metrics_df$gene.set, sep='_')), function(x){
    data.frame(F1=mean(x$F1), F1_lower=calc_CI(x$F1)[1], F1_upper=calc_CI(x$F1)[2],
    AUPRC=mean(x$AUPRC), AUPRC_lower=calc_CI(x$AUPRC)[1], AUPRC_upper=calc_CI(x$F1)[2], 
    MCC=mean(x$MCC), MCC_lower=calc_CI(x$MCC)[1], MCC_upper=calc_CI(x$MCC)[2],
    n_features=round(mean(x$n_features), 0))
})
average_ensemble_metrics <- bind_rows(average_ensemble_metrics, .id='celltype_gene.set') %>%
    separate_wider_delim(celltype_gene.set, delim = "_", names = c('celltype', 'gene.set')) %>%
    data.frame()
average_ensemble_metrics$gene.set <- factor(average_ensemble_metrics$gene.set, levels=c('chrX','SLE'))

write.csv(average_ensemble_metrics, '../results_update/average_ensemble_metrics.csv', row.names=FALSE)

pdf('../results_update/avg_ensemble_MCC_geneset.pdf')
comparisons <- list(c('chrX', 'SLE'))
ggplot(average_ensemble_metrics, aes(x=gene.set, y=MCC, colour=gene.set)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape = NA, colour='black', fill=NA) +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, 
    test='wilcox.test', color='black', step_increase=0.1) +
    theme_minimal() +
    labs(x='Gene Set', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set')
dev.off()

pdf('../results_update/avg_ensemble_MCC_forest.pdf', width=10, height=5)
ggplot(average_ensemble_metrics, aes(x=MCC, y=celltype)) +
    geom_errorbarh(aes(xmin = MCC_lower, xmax = MCC_upper)) +
    geom_vline(xintercept = 0.7, linetype = 'dotted', color='red') +
    geom_point() +
    # geom_point(data=best_model_point, color='red') +
    labs(x='MCC', y='') +
    theme(axis.text.y=element_text(size=12)) +
    facet_wrap(~gene.set, ncol=5, nrow=1)
dev.off()

top_celltypes_average <- subset(average_ensemble_metrics, celltype %in% top_celltypes)[,c('celltype', 'gene.set', 'MCC', 'MCC_lower', 'MCC_upper')]
# round 'MCC', 'MCC_lower', 'MCC_upper' to 2 decimal
top_celltypes_average$MCC <- round(top_celltypes_average$MCC, 2)
top_celltypes_average$MCC_lower <- round(top_celltypes_average$MCC_lower, 2)
top_celltypes_average$MCC_upper <- round(top_celltypes_average$MCC_upper, 2)

write.csv(top_celltypes_average, '../results_update/top_celltypes_average.csv')

subset(average_ensemble_metrics, gene.set == 'chrX' & MCC_upper > 0.7 & celltype %in% top_celltypes)[,c('celltype', 'MCC', 'MCC_lower', 'MCC_upper')]

pdf('../results_update/avg_ensemble_MCC_by_n_features.pdf', width=10, height=5)
ggplot(average_ensemble_metrics, aes(x=n_features, y=MCC)) +
    geom_point() +
    theme_minimal() +
    labs(x='Number of Features', y='MCC') +
    geom_smooth(method='lm', se=FALSE) +
    stat_cor(r.digits = 2, cor.coef.name='rho', method = "spearman") +
    facet_wrap(~gene.set, ncol=5, nrow=1, scales='free_x')
dev.off()

comparisons <- list(c('chrX', 'SLE'))
pdf('../results_update/avg_ensemble_n_features_geneset.pdf')
ggplot(average_ensemble_metrics, aes(x=gene.set, y=n_features, colour=gene.set)) +
    geom_boxplot() +
    geom_signif(comparisons=comparisons, map_signif_level=TRUE, 
    test='wilcox.test', color='black', step_increase=0.2) +
    theme_minimal() +
    labs(x='Gene Set', y='Number of Features') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=gene.set.colours, name='Gene Set')
dev.off()

pairwise.wilcox.test(average_ensemble_metrics$n_features, average_ensemble_metrics$gene.set, 
p.adjust.method='fdr')

# Add everage models and average ensemble to see if ensemble is better

average_ensemble_metrics$model <- 'ensemble'
combined_metrics <- bind_rows(average_model_metrics_df, average_ensemble_metrics)
combined_metrics$model <- factor(combined_metrics$model, levels=c('logit', 'RF', 'SVM', 'GBM', 'MLP', 'ensemble'))

pairwise.wilcox.test(combined_metrics$MCC, combined_metrics$model, p.adjust.method='fdr')

pdf('../results_update/models_ensemble_boxplot.pdf')
combinations <- list(c('logit', 'ensemble'), c('RF', 'ensemble'), c('SVM', 'ensemble'), c('GBM', 'ensemble'), c('MLP', 'ensemble'))
ggplot(combined_metrics, aes(x=model, y=MCC, colour=model)) +
    geom_jitter(width=0.2) +
    geom_boxplot(outlier.shape = NA, colour = 'black', fill = NA) +
    geom_signif(comparisons=combinations, map_signif_level=TRUE, 
    test='wilcox.test', color='black', step_increase=0.1) +
    theme_minimal() +
    labs(x='Model', y='MCC') +
    theme(axis.text.x=element_blank()) +
    scale_colour_manual(values=model.colours, name='Model')
dev.off()

### Read in the features selected across splits to visualise concordance ###
feature_list <- list()
for(i in 1:10){
    feature.files <- list.files(paste0('split_', i, '/features'), pattern='combined', full.names=TRUE)
    features <- lapply(feature.files, function(x) read.csv(x)$Feature)
    names(features) <- gsub('combined_|.csv', '', basename(feature.files))
    feature_list[[i]] <- features
}
celltype.geneset <- names(feature_list[[1]])
# Subset for chrX and SLE
celltype.geneset <- celltype.geneset[grep('.chrX|.SLE', celltype.geneset)]

result_list <- list()
for(i in celltype.geneset){
    tmp <- lapply(1:10, function(x){
    feature_list[[x]][[i]]
    })
    result_list[[i]] <- tmp
}
names(result_list) <- celltype.geneset

if (!dir.exists('../results_update/feature_heatmap')) {dir.create('../results_update/feature_heatmap')}

# Plot heatmap of selected features
for(file in celltype.geneset){
    mtx <- fromList(result_list[[file]])
    rownames(mtx) <- unique(unlist(result_list[[file]]))
    colnames(mtx) <- paste('split', 1:10, sep='_')
    # Order rows by total number of features selected
    mtx <- mtx[order(rowSums(mtx), decreasing = TRUE),]

    col_fun <- colorRamp2(c(0, 1), c("white", "black"))
    pdf(paste0('../results_update/feature_heatmap/', file, '.pdf'))
    p <- Heatmap(as.matrix(mtx), 
            name='Features', 
            cluster_columns=FALSE,
            cluster_rows=FALSE,
            col=col_fun,
            row_names_gp = gpar(fontsize = 5),
            show_heatmap_legend=FALSE)
    print(p)
    dev.off()
}

all_features <- lapply(result_list, function(x) unique(unlist(x)))
top_features <- lapply(result_list, function(x) {
    tmp <- table(unlist(x))
    top <- names(tmp[tmp == 10])
    if (length(top) == 0) list() else top
})

### Save selected features ###
selected_features <- list(all_features=all_features, top_features=top_features)
save(selected_features, file='../results_update/selected_features.RData')
load('../results_update/selected_features.RData')

# Write as .csv file
all_features_mtx <- bind_rows(lapply(names(selected_features$all_features), function(x){
    data.frame(celltype=x, feature=selected_features$all_features[[x]])
}))
write.csv(all_features_mtx, '../results_update/all_features.csv', row.names=FALSE)

selected_features$top_features[['B_cell_Naive.chrX']] <- "NA"
selected_features$top_features[['Conventional_DC.chrX']] <- "NA"

top_features_mtx <- bind_rows(lapply(names(selected_features$top_features), function(x){
    data.frame(celltype=x, feature=selected_features$top_features[[x]])
}))
write.csv(top_features_mtx, '../results_update/top_features.csv', row.names=FALSE)

#load('../results_update/selected_features.RData')

## Heatmap of selected top X chromosome genes across top celltypes
# Edit celltype names to match filenames
top_celltypes <- gsub("\\+|-| ", "_", top_celltypes)

top_celltypes_top_features <- selected_features$top_features[names(selected_features$top_features) %in% paste0(gsub(' ', '_', top_celltypes), '.chrX')]
top_celltypes_top_features_mtx <- fromList(top_celltypes_top_features)
rownames(top_celltypes_top_features_mtx) <- unique(unlist(top_celltypes_top_features))
colnames(top_celltypes_top_features_mtx) <- gsub('.chrX', '', names(top_celltypes_top_features_mtx)) %>% gsub('_', ' ', .)

pdf('../results_update/top_celltypes_top_features_heatmap.pdf')
col_fun <- colorRamp2(c(0, 1), c("white", "red"))
p <- Heatmap(as.matrix(top_celltypes_top_features_mtx), 
    name='Features', 
    col=col_fun,
    show_heatmap_legend=FALSE)
print(p)
dev.off()

## Heatmap of all selected X chromosome genes across top celltypes
top_celltypes_all_features <- selected_features$all_features[names(selected_features$all_features) %in% paste0(gsub(' ', '_', top_celltypes), '.chrX')]
top_celltypes_all_features_mtx <- fromList(top_celltypes_all_features)
rownames(top_celltypes_all_features_mtx) <- unique(unlist(top_celltypes_all_features))
colnames(top_celltypes_all_features_mtx) <- gsub('.chrX', '', names(top_celltypes_all_features_mtx)) %>% gsub('_', ' ', .)

pdf('../results_update/top_celltypes_all_features_heatmap.pdf', width=10, height=10)
col_fun <- colorRamp2(c(0, 1), c("white", "red"))
p <- Heatmap(as.matrix(top_celltypes_all_features_mtx), 
    name='Features', 
    col=col_fun,
    show_heatmap_legend=FALSE)
print(p)
dev.off()

lapply(names(top_celltypes_all_features), function(x){
    data.frame(celltype=x, feature=top_celltypes_all_features[[x]])
}) %>%
    bind_rows() %>%
    write.csv('../results_update/top_celltypes_all_features.csv', row.names=FALSE)


# Perfrom Euclidean clustering on the top celltypes all features matrix
clustering <- hclust(dist(t(top_celltypes_all_features_mtx)), method='complete')
# Print the groups
groups <- cutree(clustering, k=3)

# Read in propellor results
library(speckle)
metadata <- read.delim('../metadata.tsv')

index <- c(read.csv('split_1/train_index.csv')$rownames, read.csv('split_1/test_index.csv')$rownames)
metadata <- metadata[metadata$ind_cov %in% index,]
metadata$disease <- ifelse(metadata$disease == 'systemic lupus erythematosus', 'disease', 'control')
metadata$disease <- factor(metadata$disease, levels=c('disease', 'control'))

# Create detailed cell type column
metadata$cell_type_detailed <- case_when(
  # CD4 T cell subsets
  metadata$cell_type == "CD4-positive, alpha-beta T cell" & metadata$ct_cov == "T4_naive" ~ "CD4 T cell Naive",
  metadata$cell_type == "CD4-positive, alpha-beta T cell" & metadata$ct_cov == "T4_em"    ~ "CD4 T cell Effector Memory",
  metadata$cell_type == "CD4-positive, alpha-beta T cell" & metadata$ct_cov == "T4_reg"   ~ "CD4 T cell Treg",

  # CD8 T cell subsets
  metadata$cell_type == "CD8-positive, alpha-beta T cell" & metadata$ct_cov == "T8_naive" ~ "CD8 T cell Naive",
  metadata$cell_type == "CD8-positive, alpha-beta T cell" & metadata$ct_cov == "CytoT_GZMH+" ~ "CD8 T cell Cytotoxic GZMH+",
  metadata$cell_type == "CD8-positive, alpha-beta T cell" & metadata$ct_cov == "CytoT_GZMK+" ~ "CD8 T cell Cytotoxic GZMK+",
  metadata$cell_type == "CD8-positive, alpha-beta T cell" & metadata$ct_cov == "T_mait"   ~ "CD8 T cell MAIT",

  # NK cell subsets
  metadata$cell_type == "natural killer cell" & metadata$ct_cov == "NK_bright" ~ "NK cell Bright",
  metadata$cell_type == "natural killer cell" & metadata$ct_cov == "NK_dim"    ~ "NK cell Dim",

  # B cell subsets
  metadata$cell_type == "B cell" & metadata$ct_cov == "B_naive"    ~ "B cell Naive",
  metadata$cell_type == "B cell" & metadata$ct_cov == "B_mem"      ~ "B cell Memory",
  metadata$cell_type == "B cell" & metadata$ct_cov == "B_plasma"   ~ "B cell Plasma",
  metadata$cell_type == "B cell" & metadata$ct_cov == "B_atypical" ~ "B cell Atypical",

  # Progenitor cells
  metadata$cell_type == "progenitor cell" ~ "Progenitor cell",

  # Dendritic cells
  metadata$cell_type == "conventional dendritic cell" ~ "Conventional DC",
  metadata$cell_type == "plasmacytoid dendritic cell" ~ "Plasmacytoid DC",

  # Monocytes
  metadata$cell_type == "classical monocyte"       ~ "Classical monocyte",
  metadata$cell_type == "non-classical monocyte"   ~ "Non-classical monocyte",

  # Fallback for lymphocyte and plasmablast without ct_cov
  metadata$cell_type == "lymphocyte"               ~ "Lymphocyte",
  metadata$cell_type == "plasmablast"              ~ "Plasmablast",

  # Default to cell_type for anything else
  TRUE ~ metadata$cell_type
)

output.asin <- propeller(clusters=metadata$cell_type_detailed, sample=metadata$ind_cov, group=metadata$disease, transform='asin')
output.asin$BaselineProp.clusters <- gsub("-| ", "_", output.asin$BaselineProp.clusters)
foo <- subset(output.asin, BaselineProp.clusters %in% top_celltypes & FDR < 0.05)[,c('BaselineProp.clusters','Tstatistic', 'FDR')]
foo$Tstatistic <- round(foo$Tstatistic, 2)
foo$FDR <- round(foo$FDR, 6)

save(output.asin, file='../results_update/propellor_results.RData')

# Match rownames to escape genes
hits <- rownames(top_celltypes_all_features_mtx)[rownames(top_celltypes_all_features_mtx) %in% escape]

selected <- selected_features$all_features[["CD8_positive,_alpha_beta_T_cell.chrX"]]

a <- sum(selected %in% escape)
b <- sum(!(selected %in% escape))
c <- sum(chrX[!chrX %in% selected] %in% escape)
d <- sum(!(chrX[!chrX %in% selected] %in% escape))
chisq.test(matrix(c(a, b, c, d), nrow=2))


### Read in edgeR results and subset top_features for DEGs ###
source('../edgeR.list.R')
edgeR <- deg.list('../edgeR', filter=FALSE)

top_celltypes_all_features <- selected_features$all_features[c('CD4_positive,_alpha_beta_T_cell.chrX',
'CD8_positive,_alpha_beta_T_cell.chrX')]
edgeR_top_celltypes_all_features <- lapply(names(top_celltypes_all_features), function(x){
    tmp <- edgeR[[gsub('.chrX', '', x)]]
    subset(tmp, gene %in% top_celltypes_all_features[[x]])[,c('gene', 'logFC', 'FDR')]
})
names(edgeR_top_celltypes_all_features) <- names(top_celltypes_all_features)

### Heatmap of all chrX features with genes labelled ###
library('UpSetR')
library('ComplexHeatmap')
library(circlize)

load('results_update/selected_features.RData')
load('results_update/top_celltypes.RData')
SLE <- read.delim('/directflow/SCCGGroupShare/projects/lacgra/DisGeNet/SLE.tsv')$Gene
load('/directflow/SCCGGroupShare/projects/lacgra/datasets/XCI/escapees.Rdata')


chrX.list <- selected_features$all_features[grep('.chrX',
names(selected_features$all_features), value = TRUE)]
names(chrX.list) <- gsub('.chrX', '', names(chrX.list))

names(chrX.list)[9] <- "CD8 positive, alpha-beta T cell"
names(chrX.list)[5] <- "CD4 positive, alpha-beta T cell"
names(chrX.list) <- gsub('_', ' ', names(chrX.list))

top_list <- chrX.list[top_celltypes]
names(top_list)[1] <- "CD4 positive, alpha-beta T cell"
names(top_list)[2] <- "CD8 positive, alpha-beta T cell"
names(top_list)[3] <- "Progenitor cell"

mat <- fromList(top_list)
rownames(mat) <- unique(unlist(top_list))

annotation <- ifelse(rownames(mat) %in% escape,
                     'Escapee', 
                     'Non-escapee')

mat <- mat[order(annotation, decreasing = FALSE), ]
annotation <- annotation[order(annotation, decreasing = FALSE)]

row_ha <- rowAnnotation(
  Feature_type = annotation,
  col = list(Feature_type = c("Escapee" = "#2284F5", "Non-escapee" = "#F56022")),
  show_annotation_name = FALSE
)

pdf('../results_update/heatmap_chrX_features_clusters.pdf', width=10, height=10)
col_fun <- colorRamp2(c(0, 1), c("white", "red"))
p <- Heatmap(as.matrix(mat), 
    name='Features', 
    col=col_fun,
    show_heatmap_legend=FALSE,
    cluster_rows=FALSE,
    right_annotation = row_ha)
print(p)
dev.off()


### Predicting independent data metrics ###
library(dplyr)
library(ggplot2)
library(stringr)


metrics.files <- list.files('results_update', pattern ='csv', full.names = TRUE)

metrics <- lapply(metrics.files, read.csv)
names(metrics) <- gsub('.csv', '', basename(metrics.files))

metrics_df <- bind_rows(metrics, .id = 'filename')

metrics_df$split <- as.numeric(gsub('.*_split_(.*)', '\\1', metrics_df$filename))
metrics_df$age <- str_extract(metrics_df$filename, 'adult|child')
metrics_df$age <- factor(metrics_df$age, levels = c('adult', 'child'))
metrics_df$geneset <- str_extract(metrics_df$filename, 'chrX|SLE')
metrics_df$geneset <- factor(metrics_df$geneset, levels = c('chrX', 'SLE'))
metrics_df$celltype <- gsub('^metrics_|\\..*', '', metrics_df$filename)
metrics_df$celltype <- gsub('_', ' ', metrics_df$celltype)

write.table(metrics_df, 'results_update/predicting_metrics.txt', row.names = FALSE)
metrics_df <- read.table('results_update/predicting_metrics.txt', header=TRUE)

pdf('results_update/predicting_boxplot.pdf')
ggplot(metrics_df, aes(x=MCC, y=celltype, colour=geneset)) +
  geom_boxplot() +
  scale_colour_manual(values=c('#8A0798', '#D90750')) +
  geom_vline(xintercept=0.7, linetype='dashed', colour='grey') +
  facet_wrap(~age) +
  labs(title='') +
  theme_minimal() +
  theme(legend.position = 'right')
dev.off()

# test for adult chrX vs SLE
adult.test <- lapply(unique(metrics_df$celltype), function(x){
    tmp <- subset(metrics_df, age == 'adult' & celltype == x)
    wilcox.test(MCC ~ geneset, tmp)
})
names(adult.test) <- unique(metrics_df$celltype)
adult.test <- do.call(rbind, lapply(adult.test, function(x) data.frame(statistic = x$statistic, p.value = x$p.value)))
adult.test$FDR <- p.adjust(adult.test$p.value, method = 'fdr')
subset(adult.test, FDR > 0.05)

child.test <- lapply(unique(metrics_df$celltype), function(x){
    tmp <- subset(metrics_df, age == 'child' & celltype == x)
    wilcox.test(MCC ~ geneset, alternative='greater', tmp)
})
names(child.test) <- unique(metrics_df$celltype)
child.test <- do.call(rbind, lapply(child.test, function(x) data.frame(statistic = x$statistic, p.value = x$p.value)))
child.test$FDR <- p.adjust(child.test$p.value, method = 'fdr')
subset(child.test, FDR < 0.05)

# Calculate the mean and 95% CI of MCC for each cell type and geneset
mean_metrics <- metrics_df %>%
  group_by(celltype, geneset, age) %>%
  summarise(mean_MCC = round(mean(MCC), 2), 
            lower_CI = round(quantile(MCC, 0.025), 2), 
            upper_CI = round(quantile(MCC, 0.975), 2)) %>%
  ungroup() %>%
  data.frame()

subset(mean_metrics, age == 'adult' & mean_MCC > 0.7)

#### Analysing FLARE predictions ####
library(dplyr)
library(ggplot2)
library(stringr)

setwd('/directflow/SCCGGroupShare/projects/lacgra/SLE')


metrics.files <- list.files('FLARE/figures', pattern ='csv', full.names = TRUE)

metrics <- lapply(metrics.files, read.csv)
names(metrics) <- gsub('metrics_|.csv', '', basename(metrics.files))

metrics_df <- bind_rows(metrics, .id = 'filename')

metrics_df$split <- as.numeric(gsub('.*_split_(.*)', '\\1', metrics_df$filename))
metrics_df$geneset <- str_extract(metrics_df$filename, 'chrX|SLE')
metrics_df$geneset <- factor(metrics_df$geneset, levels = c('chrX', 'SLE'))
metrics_df$celltype <- gsub('^metrics_|\\..*', '', metrics_df$filename)
metrics_df$celltype <- gsub('_', ' ', metrics_df$celltype)

write.csv(metrics_df, 'FLARE/predicting_FLARE_metrics.csv', row.names = FALSE)

pdf('FLARE/predicting_FLARE_boxplot.pdf')
ggplot(metrics_df, aes(x=MCC, y=celltype, colour=geneset)) +
  geom_boxplot() +
  scale_colour_manual(values=c('#8A0798', '#D90750')) +
  geom_vline(xintercept=0.7, linetype='dashed', colour='grey') +
  facet_wrap(~geneset) +
  labs(title='') +
  theme_minimal() +
  theme(legend.position = 'right')
dev.off()

# test for  chrX vs SLE
testing <- lapply(unique(metrics_df$celltype), function(x){
    tmp <- subset(metrics_df, celltype == x)
    wilcox.test(MCC ~ geneset, tmp)
})
names(testing) <- unique(metrics_df$celltype)
testing <- do.call(rbind, lapply(testing, function(x) data.frame(statistic = x$statistic, p.value = x$p.value)))
testing$FDR <- p.adjust(testing$p.value, method = 'fdr')
subset(testing, FDR > 0.05)


# Calculate the mean and 95% CI of MCC for each cell type and geneset
mean_metrics <- metrics_df %>%
  group_by(celltype, geneset, age) %>%
  summarise(mean_MCC = round(mean(MCC), 2), 
            lower_CI = round(quantile(MCC, 0.025), 2), 
            upper_CI = round(quantile(MCC, 0.975), 2)) %>%
  ungroup() %>%
  data.frame()

subset(mean_metrics, age == 'child' & mean_MCC > 0.7)


#########################

### calculate jaccard index between gene sets ###
# Jaccard index
jaccard_index <- function(x, y){
    intersect <- length(intersect(x, y))
    union <- length(union(x, y))
    return(intersect / union)
}

ifelse(dir.exists('../results_update/jaccard_heatmaps') == F, dir.create('../results_update/jaccard_heatmaps'))

celltypes <- unique(gsub('.HVG.autosome|.HVG|.autosome|.SLE|.chrX', '', names(selected_features$all_features)))
for(celltype in celltypes){
    tmp <- selected_features$all_features[grep(paste0('^',celltype), names(selected_features$all_features))]
    combinations <- combn(names(tmp), 2)
    mtx <- matrix(0, ncol=5, nrow=5, dimnames = list(names(tmp), names(tmp)))
    for(i in 1:ncol(combinations)){
        mtx[combinations[1,i], combinations[2,i]] <- jaccard_index(tmp[[combinations[1,i]]], tmp[[combinations[2,i]]])
        mtx[combinations[2,i], combinations[1,i]] <- mtx[combinations[1,i], combinations[2,i]]
    }
    diag(mtx) <- 1
    colnames(mtx) <- gsub(celltype, '', colnames(mtx)) %>% gsub('^\\.', '', .)
    rownames(mtx) <- gsub(celltype, '', rownames(mtx)) %>% gsub('^\\.', '', .)
    pdf(paste0('../results_update/jaccard_heatmaps/', celltype, '.pdf'))
    col = circlize::colorRamp2(c(0, 1), c("white", "red"))
    print(Heatmap(mtx, col=col, name = 'Jaccard Index', column_title=replace.names('CD16+.NK.cells'),
    cluster_columns=FALSE, cluster_rows=FALSE))
    dev.off()
}

result_list <- list()
for(celltype in celltypes){
    result <- lapply(c('.HVG', '.SLE'), function(x){
    jaccard_index(selected_features$all_features[[paste0(celltype, x)]], selected_features$all_features[[paste0(celltype, '.chrX')]])
    })
    names(result) <- c('HVG', 'SLE')
    tmp <- do.call(rbind, result)
    colnames(tmp) <- celltype
    result_list[[celltype]] <- tmp
}
jaccard_matrix <- do.call(cbind, result_list)
colnames(jaccard_matrix) <- gsub('_', ' ', colnames(jaccard_matrix))

pdf('../results_update/chrX.geneset.jaccard.heatmap.pdf')
col = circlize::colorRamp2(c(0, 1), c("white", "red"))
Heatmap(t(jaccard_matrix), col=col, name = 'Jaccard Index',
cluster_columns=FALSE, cluster_rows=FALSE)
dev.off()

overlap_list <- list()
for(celltype in celltypes){
    X_genes <- selected_features$all_features[[paste0(celltype, '.chrX')]]
    num_chrX <- length(X_genes)
    HVG_overlap <- length(intersect(selected_features$all_features[[paste0(celltype, '.HVG')]], X_genes))
    HVG_Jaccard <- jaccard_index(selected_features$all_features[[paste0(celltype, '.HVG')]], X_genes)
    SLE_overlap <- length(intersect(selected_features$all_features[[paste0(celltype, '.SLE')]], X_genes))
    SLE_Jaccard <- jaccard_index(selected_features$all_features[[paste0(celltype, '.SLE')]], X_genes)
    DEG_overlap <- length(intersect(X_genes, subset(edgeR[[celltype]], abs(logFC) > 0.1 & FDR < 0.05)$gene))
    DEG_jaccard <- jaccard_index(X_genes, subset(edgeR[[celltype]], abs(logFC) > 0.1 & FDR < 0.05)$gene)
    overlap_list[[celltype]] <- data.frame(celltype=gsub('_', ' ', celltype), num_chrX=num_chrX, HVG_overlap=HVG_overlap, 
    HVG_Jaccard=round(HVG_Jaccard,4), SLE_overlap=SLE_overlap, SLE_Jaccard=round(SLE_Jaccard,4), DEG_overlap=DEG_overlap, 
    DEG_jaccard=round(DEG_jaccard,4), row.names=NULL)
}
overlap_df <- bind_rows(overlap_list)
write.csv(overlap_df, '../results_update/chrX_model_overlap_df.csv', row.names=FALSE)

intersect_list <- list()
for(celltype in celltypes){
    result <- unlist(lapply(c('.HVG', '.SLE'), function(x){
    intersect(selected_features$all_features[[paste0(celltype, x)]], selected_features$all_features[[paste0(celltype, '.chrX')]])
    }))
    intersect_list[[celltype]] <- result
}
sort(table(unlist(intersect_list)))

### Read in edgeR results and subset top_features for DEGs ###
source('/directflow/SCCGGroupShare/projects/lacgra/PhD/functions/edgeR.list.R')
edgeR <- deg.list('../edgeR', filter=FALSE)

deg_features_list <- list()
for(i in seq_along(edgeR)){
    tmp <- lapply(grep(names(edgeR)[i], names(selected_features$all_features), value = TRUE), function(x){
        subset(edgeR[[i]], gene %in% selected_features$all_features[[x]])
    })
    names(tmp) <- names(selected_features$all_features)[grep(names(edgeR)[i], names(selected_features$all_features))]
    deg_features_list <- c(deg_features_list, tmp)
}

deg_top_features <- deg_features_list[names(deg_features_list) %in% paste0(gsub(' ', '_', top_celltypes), '.chrX')]

# Create heatmap of cell type by gene coloured by logFC
combined_deg_top_features <- bind_rows(deg_top_features, .id='celltype')
combined_deg_top_features <- combined_deg_top_features[-2,]
degs_mtx <- reshape2::dcast(combined_deg_top_features, gene ~ celltype, value.var='logFC')
rownames(degs_mtx) <- degs_mtx$gene
degs_mtx <- degs_mtx[,-1]
degs_mtx[is.na(degs_mtx)] <- 0
colnames(degs_mtx) <- replace.names(colnames(degs_mtx))

pdf('../results_update/chrX_DEG_heatmap.pdf')
col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
Heatmap(as.matrix(degs_mtx), col=col, name = 'logFC',
        cluster_columns=TRUE, cluster_rows=TRUE)
dev.off()


degs <- lapply(celltypes, function(x){
    subset(edgeR[[x]], gene %in% intersect_list[[x]])
})
names(degs) <- celltypes

combined_degs <- bind_rows(degs, .id='celltype')
degs_mtx <- reshape2::dcast(combined_degs, gene ~ celltype, value.var='logFC')
rownames(degs_mtx) <- degs_mtx$gene
degs_mtx <- degs_mtx[,-1]
degs_mtx[is.na(degs_mtx)] <- 0
colnames(degs_mtx) <- replace.names(colnames(degs_mtx))

fdr_matrix <- reshape2::dcast(combined_degs, gene ~ celltype, value.var='FDR')
rownames(fdr_matrix) <- fdr_matrix$gene
colnames(fdr_matrix) <- replace.names(colnames(fdr_matrix))
fdr_matrix <- fdr_matrix[,-1]
significance_matrix <- ifelse(fdr_matrix < 0.05, "*", "ns")
significance_matrix <- ifelse(fdr_matrix < 0.01, "**", significance_matrix)
significance_matrix <- ifelse(fdr_matrix < 0.001, "***", significance_matrix)
significance_matrix[is.na(significance_matrix)] <- ""

pdf('results_update/shared_features_deg_heatmap.pdf')
col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
Heatmap(as.matrix(degs_mtx), col=col, name = 'logFC',
        cluster_columns=TRUE, cluster_rows=TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(significance_matrix[i, j], x, y, gp = gpar(fontsize = 10))
        })
dev.off()

### Calculate Jaccard index between chrX features ###
chrX_features <- selected_features$all_features[grep('.chrX', names(selected_features$all_features))]
names(chrX_features) <- gsub('.chrX', '', names(chrX_features))

top_celltypes <- c('CD16+.NK.cells', 'Memory.B.cells','Tcm.Naive.helper.T.cells', 'Regulatory.T.cells', 'DC1', 'DC2', 'pDC', 'Tem.Effector.helper.T.cells', 'Non.classical.monocytes')
chrX_features <- chrX_features[top_celltypes]

sort(table(unlist(chrX_features)))

combinations <- combn(names(chrX_features), 2)
jaccard_mtx <- matrix(0, ncol=9, nrow=9, dimnames = list(names(chrX_features), names(chrX_features)))
for(i in 1:ncol(combinations)){
    jaccard_mtx[combinations[1, i], combinations[2, i]] <- jaccard_index(chrX_features[[combinations[1, i]]], chrX_features[[combinations[2, i]]])
    jaccard_mtx[combinations[2, i], combinations[1, i]] <- jaccard_index(chrX_features[[combinations[1, i]]], chrX_features[[combinations[2, i]]])
}
diag(jaccard_mtx) <- 1

colnames(jaccard_mtx) <- replace.names(colnames(jaccard_mtx))
rownames(jaccard_mtx) <- replace.names(rownames(jaccard_mtx))
jaccard_mtx[lower.tri(jaccard_mtx)] <- 0

pdf('results_update/chrX_jaccard_heatmap.pdf')
col = circlize::colorRamp2(c(0, 1), c("white", "red"))
Heatmap(jaccard_mtx, col=col, name = 'Jaccard Index')
dev.off()

degs <- lapply(top_celltypes, function(x){
    subset(edgeR[[x]], gene %in% chrX_features[[x]])
})
names(degs) <- top_celltypes

combined_degs <- bind_rows(degs, .id='celltype')
combined_degs$is_escape <- ifelse(combined_degs$gene %in% rownames(escape), 'True', 'False')
write.csv(combined_degs, 'results_update/Supplementary_Table_2.csv', row.names=FALSE)

degs_mtx <- reshape2::dcast(combined_degs, gene ~ celltype, value.var='logFC')
rownames(degs_mtx) <- degs_mtx$gene
degs_mtx <- degs_mtx[,-1]
degs_mtx[is.na(degs_mtx)] <- 0
colnames(degs_mtx) <- replace.names(colnames(degs_mtx))

fdr_matrix <- reshape2::dcast(combined_degs, gene ~ celltype, value.var='FDR')
rownames(fdr_matrix) <- fdr_matrix$gene
colnames(fdr_matrix) <- replace.names(colnames(fdr_matrix))
fdr_matrix <- fdr_matrix[,-1]
significance_matrix <- ifelse(fdr_matrix < 0.05, "*", "ns")
significance_matrix <- ifelse(fdr_matrix < 0.01, "**", significance_matrix)
significance_matrix <- ifelse(fdr_matrix < 0.001, "***", significance_matrix)
significance_matrix[is.na(significance_matrix)] <- ""

top_genes <- names(sort(table(combined_degs$gene), decreasing=TRUE))[1:20]

pdf('results_update/chrX_heatmap.pdf')
col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
ann <- rowAnnotation(foo = anno_mark(at = which(rownames(degs_mtx) %in% top_genes), 
                                    labels = rownames(degs_mtx)[rownames(degs_mtx) %in% top_genes]))
Heatmap(as.matrix(degs_mtx), col=col, name = 'logFC',
        cluster_columns=TRUE, cluster_rows=TRUE,
        show_row_names=FALSE, right_annotation = ann,
        # cell_fun = function(j, i, x, y, width, height, fill) {
        #     grid.text(significance_matrix[i, j], x, y, gp = gpar(fontsize = 2))
        # })
        )
dev.off()


potential_escapees <- subset(combined_degs, abs(logFC) >= 0.5 & FDR < 0.05)
potential_escapees <- potential_escapees[,c('celltype', 'gene', 'logFC', 'FDR')]
unique(potential_escapees$gene)

combined_degs$celltype <- replace.names(combined_degs$celltype)

library(ggrepel)

pdf('../results_update/potential_escapees.pdf')
ggplot(combined_degs, aes(x=logFC, y=-log10(FDR))) +
    geom_point(aes(color=ifelse(abs(logFC) > 0.5 & FDR < 0.05, 'red', 'grey'))) +
    scale_color_identity() +  # Use colors as provided
    geom_text_repel(data=subset(combined_degs, abs(logFC) > 0.5 & FDR < 0.05), 
                    aes(label=gene), color='black', max.overlaps=Inf, size=2) +
    theme_minimal() +
    theme(legend.position='none') +
    labs(x='logFC', y='-log10(FDR)') +
    geom_vline(xintercept = c(-0.5, 0.5), linetype='dotted', color='black') +
    geom_hline(yintercept = -log10(0.05), linetype='dotted', color='black') +
    facet_wrap(~celltype, ncol=4)
dev.off()


combined_degs_escape <- subset(combined_degs, gene %in% rownames(escape))
escape_mtx <- reshape2::dcast(combined_degs_escape, gene ~ celltype, value.var='logFC')
rownames(escape_mtx) <- escape_mtx$gene
escape_mtx <- escape_mtx[,-1]
escape_mtx[is.na(escape_mtx)] <- 0
colnames(escape_mtx) <- replace.names(colnames(escape_mtx))

escape_fdr_matrix <- reshape2::dcast(combined_degs_escape, gene ~ celltype, value.var='FDR')
rownames(escape_fdr_matrix) <- escape_fdr_matrix$gene
colnames(escape_fdr_matrix) <- replace.names(colnames(escape_fdr_matrix))
escape_fdr_matrix <- escape_fdr_matrix[,-1]
escape_significance_matrix <- ifelse(escape_fdr_matrix < 0.05, "*", "ns")
escape_significance_matrix <- ifelse(escape_fdr_matrix < 0.01, "**", escape_significance_matrix)
escape_significance_matrix <- ifelse(escape_fdr_matrix < 0.001, "***", escape_significance_matrix)
escape_significance_matrix[is.na(escape_significance_matrix)] <- ""

pdf('results_update/escape_heatmap.pdf')
col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
# ann <- rowAnnotation(foo = anno_mark(at = which(rownames(escape_mtx) %in% top_genes), 
#                                     labels = rownames(degs_mtx)[rownames(degs_mtx) %in% top_genes]))
Heatmap(as.matrix(escape_mtx), col=col, name = 'logFC',
        cluster_columns=TRUE, cluster_rows=TRUE,
        show_row_names=FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(escape_significance_matrix[i, j], x, y, gp = gpar(fontsize = 2))
        })
dev.off()

subset(edgeR[['Regulatory.T.cells']], gene %in% c('PDCD1', 'CTLA4', 'TIM3', 'LAG3', 'BTLA','TIGIT'))[,c('gene', 'logFC', 'FDR')]

### Consistently selected features ###
chrX_features <- selected_features$top_features[grep('.chrX', names(selected_features$top_features))]
names(chrX_features) <- gsub('.chrX', '', names(chrX_features))
chrX_features <- chrX_features[top_celltypes]

degs <- lapply(top_celltypes, function(x){
    tmp <- merge(edgeR[[x]], data.frame('gene'=chrX_features[[x]]), by='gene', all.y=TRUE)
    return(tmp[,c('gene', 'logFC', 'FDR')])
})
names(degs) <- top_celltypes

combined_degs <- dplyr::bind_rows(degs, .id='celltype')
combined_degs$celltype <- replace.names(combined_degs$celltype)
write.csv(combined_degs, 'results_update/top_chrX.consistent.csv', row.names=FALSE)

plots_list <- list()
for(cell in names(chrX_features)){
    features <- chrX_features[[cell]]
    mtx <- readRDS(paste0(cell, '.chrX.RDS'))
    mtx <- mtx[, c('class', features)]
    class <- mtx$class
    mtx_scaled <- scale(as.matrix(mtx[,-1]))
    rownames(mtx_scaled) <- class
    mtx_melt <- reshape2::melt(mtx_scaled)
    mtx_melt$Var1 <- factor(mtx_melt$Var1, levels=c('control', 'disease'))
    
    
    plot <- ggplot(mtx_melt, aes(x=Var2, y=value, color=Var1)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(position=position_jitterdodge(), alpha=0.5) +
        theme_minimal() +
        labs(x='', y='z-score', color='Condition') +
        ggtitle(replace.names(cell)) +
        theme(plot.title = element_text(size = 18)) +
        theme(axis.text.x = element_text(size=12)) +
        scale_color_manual(values=c('blue', 'red'), )
    plots_list[[cell]] <- plot
}

pdf('results_update/consistent_features_boxplot.pdf', width=15, height=20)
grid.arrange(grobs = plots_list, ncol = 2, nrow = 5)
dev.off()

cell <- 'CD16+.NK.cells'
mtx <- readRDS(paste0(cell, '.chrX.RDS'))
mtx_scaled <- data.frame(scale(mtx[, c('XIST', 'TSIX')]))
mtx_scaled$class <- mtx$class

coorelation_result <- lapply(split(mtx_scaled, mtx_scaled$class), function(x){
    cor.test(x$XIST, x$TSIX, method='spearman')
})
# Create a data frame for annotations
cor_df <- data.frame(
    class = names(coorelation_result),
    cor_value = sapply(coorelation_result, function(x) x$estimate),
    p_value = sapply(coorelation_result, function(x) x$p.value)
)

pdf('results_update/CD16NK_XIST_TSIX.pdf')
ggplot(data.frame(mtx_scaled), aes(x=XIST, y=TSIX)) +
    geom_point() +
    geom_smooth(method='lm', se=FALSE) +
    theme_minimal() +
    labs(x='XIST', y='TSIX') +
    facet_wrap(~class) +
    geom_text(data=cor_df, aes(x=Inf, y=Inf, label=paste0(
        'r = ', round(cor_value, 2),
        ', p = ', round(p_value, 2))),
        hjust=1.1, vjust=1.1, size=3, color='black')
dev.off()