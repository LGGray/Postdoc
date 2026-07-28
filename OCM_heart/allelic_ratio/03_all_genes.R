# ---------------------------------------------------------------------------
# 03 - Per-gene ("single gene") allelic ratio analysis.
#
# Run from the OCM_heart/ directory:  Rscript allelic_ratio/03_all_genes.R
# Requires: 01_setup.R and 02_whole_chrX.R to have been run.
# ---------------------------------------------------------------------------
source("allelic_ratio/00_functions.R")

heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

# per-cell whole-chrX table, written by 02_whole_chrX.R
metadata_whole_chr <- read.table(
  'Allelic_ratio_results/whole_chr_cell_metadata.txt',
  sep = '\t', header = TRUE, row.names = 1,
  stringsAsFactors = FALSE, check.names = FALSE)

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

