library(Seurat)
library(dplyr)
library(ggplot2)

heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

barcodes <- lapply(c('9w', '78w', 'Sham', 'TAC'), function(x) {
  list.files(paste0('Allelome.PRO2/', x), pattern = 'locus_table.txt', recursive = TRUE, full.names = TRUE)
})
barcodes <- unlist(barcodes)

condition <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[4])
cellid <- strsplit(dirname(barcodes), '_') %>% sapply(function(x) x[5])

subset_heart <- subset(heart, cells = paste0(condition, "_", cellid))

allelic_ratio <- lapply(barcodes, function(x) {
  if (file.exists(x)) {
    tmp <- read.delim(x, header = TRUE)
    subset(tmp, chr=="chrX")$allelic_ratio
  } else {
    NULL
  }
})
allelic_ratio <- unlist(allelic_ratio)
names(allelic_ratio) <- paste0(condition, "_", cellid)
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


# filter for cells with at least 10 reads on chrX
subset_heart <- subset(subset_heart, subset = total_reads >= 10)


# Violin plot by celltype and condition - facet_wrap so each cell type has its own panel
pdf('allelic_ratio_celltype_violin_plot.pdf')
VlnPlot(subset_heart, features = "allelic_ratio", group.by = "celltype", split.by = "sample")
dev.off()

pdf('allelic_ratio_celltype_violin_plot_facet_wrap.pdf')
ggplot(subset_heart@meta.data, aes(x = sample, y = allelic_ratio, fill = sample)) +
  geom_violin() +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# UMAP plot coloured by allelic ratio split by sample
pdf('allelic_ratio_umap_plot.pdf')
FeaturePlot(subset_heart, features = "allelic_ratio", split.by = "sample")
dev.off()


metadata <- subset_heart@meta.data
metadata$partial_escape <- metadata$allelic_ratio < 0.9

library(broom)
adult_vs_aged <- split(metadata, metadata$celltype) %>%
  lapply(function(x) {

    df <- subset(x, sample %in% c("9w", "78w"))
    df$sample <- factor(df$sample, levels = c("9w", "78w"))

    if (length(unique(df$sample)) < 2) return(NULL)

    tidy(glm(partial_escape ~ sample + log1p(total_reads), family = binomial, data = df))
})
# adjust p-values for multiple testing
adult_vs_aged <- lapply(adult_vs_aged, function(x) {
  if (is.null(x)) return(NULL)
  x$fdr <- p.adjust(x$p.value, method = "fdr")
  x
})

sham_vs_tac <- split(metadata, metadata$celltype) %>%
  lapply(function(x) {

    df <- subset(x, sample %in% c("Sham", "TAC"))
    df$sample <- factor(df$sample, levels = c("Sham", "TAC"))

    if (length(unique(df$sample)) < 2) return(NULL)

    tidy(glm(partial_escape ~ sample + log1p(total_reads), family = binomial, data = df))
})
sham_vs_tac <- lapply(sham_vs_tac, function(x) {
  if (is.null(x)) return(NULL)
  x$fdr <- p.adjust(x$p.value, method = "fdr")
  x
})

# number fo cells per celltype and condition
cell_counts <- metadata %>%
  group_by(celltype, sample) %>%
  summarise(n_cells = n())
write.table(cell_counts, 'cell_counts_per_celltype_and_condition.txt', sep = '\t', row.names = FALSE, quote = FALSE)

saveRDS(subset_heart, 'heart_seurat_object_SCT_AR.rds')


library(glmmTMB)
library(broom.mixed)

fit_bb <- function(df, ref_level) {
  df$sample <- factor(df$sample, levels = ref_level)
  glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
          family = betabinomial(link = "logit"),
          data = df)
}

adult_vs_aged_bb <- split(metadata, metadata$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("9w", "78w"))
    if (length(unique(df$sample)) < 2) return(NULL)
    fit <- fit_bb(df, c("9w", "78w"))
    tidy(fit, effects = "fixed")
  })

adult_vs_aged_bb <- lapply(adult_vs_aged_bb, function(x) {
  if (is.null(x)) return(NULL)
  x$fdr <- p.adjust(x$p.value, method = "fdr")
  x
})

sham_vs_tac_bb <- split(metadata, metadata$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("Sham", "TAC"))
    if (length(unique(df$sample)) < 2) return(NULL)
    fit <- fit_bb(df, c("Sham", "TAC"))
    tidy(fit, effects = "fixed")
  })
sham_vs_tac_bb <- lapply(sham_vs_tac_bb, function(x) {
  if (is.null(x)) return(NULL)
  x$fdr <- p.adjust(x$p.value, method = "fdr")
  x
})

adult_vs_aged_disp <- split(metadata, metadata$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("9w", "78w"))
    if (length(unique(df$sample)) < 2) return(NULL)
    glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                        dispformula = ~ sample,
                        family = betabinomial(link = "logit"),
                        data = df)
  })

Sham_vs_TAC_disp <- split(metadata, metadata$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("Sham", "TAC"))
    if (length(unique(df$sample)) < 2) return(NULL)
    glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                        dispformula = ~ sample,
                        family = betabinomial(link = "logit"),
                        data = df)
  })

adult_vs_aged_lrt <- split(metadata, metadata$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("9w", "78w"))
    if (length(unique(df$sample)) < 2) return(NULL)
  fit_const_disp <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                              family = betabinomial(link = "logit"), data = df)
  fit_var_disp   <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                              dispformula = ~ sample,
                              family = betabinomial(link = "logit"), data = df)
  anova(fit_const_disp, fit_var_disp)
  })
adult_vs_aged_lrt <- lapply(adult_vs_aged_lrt, function(x) {
  if (is.null(x)) return(NULL)
  x$fdr <- p.adjust(x$`Pr(>Chisq)`[2], method = "fdr")
  x
})

Sham_vs_TAC_lrt <- split(metadata, metadata$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("Sham", "TAC"))
    if (length(unique(df$sample)) < 2) return(NULL)
  fit_const_disp <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                              family = betabinomial(link = "logit"), data = df)
  fit_var_disp   <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                              dispformula = ~ sample,
                              family = betabinomial(link = "logit"), data = df)
  anova(fit_const_disp, fit_var_disp)
  })
Sham_vs_TAC_lrt <- lapply(Sham_vs_TAC_lrt, function(x) {
  if (is.null(x)) return(NULL)
  x$fdr <- p.adjust(x$`Pr(>Chisq)`[2], method = "fdr")
  x
})