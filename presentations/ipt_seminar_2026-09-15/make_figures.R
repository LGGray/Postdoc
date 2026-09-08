# ---------------------------------------------------------------------------
# Figures for the IPT seminar / lab meeting, 15 Sep 2026.
#
# Reads ONLY summary tables already produced on the cluster (mounted locally at
# /Users/lachlang/cluster) and redraws them in one consistent style for slides.
# Nothing here recomputes an analysis; every number comes from an existing
# OCM_heart / spatial / multiome output. Runs on the laptop:
#
#   Rscript presentations/ipt_seminar_2026-09-15/make_figures.R
#
# The per-gene pseudobulk heatmap (F09) additionally needs the consolidated
# per-nucleus locus tables written by consolidate_gene_tables.py (same folder).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(forcats); library(ggrepel)
})

CL  <- Sys.getenv("CLUSTER_MOUNT", "/Users/lachlang/cluster")
OUT <- Sys.getenv("FIG_OUT", "/Users/lachlang/Downloads/IPT_Seminar_2026-09-15/figures")
DAT <- Sys.getenv("FIG_DATA", "/Users/lachlang/Downloads/IPT_Seminar_2026-09-15/data")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
OCM <- file.path(CL, "OCM"); SPA <- file.path(CL, "adult_aged_spatial"); MUL <- file.path(CL, "adult_aged_multiome")

MONO_AR <- 0.90          # OCM_heart/allelic_ratio/00_functions.R
MIN_TOTAL_READS <- 30    # cutoff_30 directory

# ---- style -----------------------------------------------------------------
SAMPLE_COL <- c("9w" = "#2B7BBA", "78w" = "#E2711D", "Sham" = "#8C8C8C", "TAC" = "#7B3294")
SAMPLE_LAB <- c("9w" = "Adult (9w)", "78w" = "Aged (78w)", "Sham" = "Sham", "TAC" = "TAC")
# Allelic-ratio palette used throughout OCM_heart (blue = biallelic/escape ... red = monoallelic)
AR_BREAKS <- c(seq(0, 0.9, by = 0.1), 0.95, 1.0)
AR_COLS   <- c("#2B3186", "#3B5FB6", "#38749F", "#367373", "#2D6E5D", "#1E652D",
               "#658C2D", "#8D9F25", "#B3B112", "#C97314", "#8B1913")
AR_LABELS <- c("0.00-0.10","0.10-0.20","0.20-0.30","0.30-0.40","0.40-0.50","0.50-0.60",
               "0.60-0.70","0.70-0.80","0.80-0.90","0.90-0.95","0.95-1.00")
ar_bin <- function(x) cut(x, breaks = AR_BREAKS, include.lowest = TRUE, right = TRUE, labels = AR_LABELS)
ar_fill_cont <- function(name = "Allelic ratio\n(B6 / total)") {
  mids <- (head(AR_BREAKS, -1) + tail(AR_BREAKS, -1)) / 2
  scale_fill_gradientn(colours = AR_COLS, values = mids, limits = c(0, 1), name = name, na.value = "white")
}
ESCAPE_GENES <- c("Kdm5c","Kdm6a","Ddx3x","Eif2s3x","Utp14a","Akap17a","Pbdc1","Ftx","Jpx","Sts","5530601H04Rik")

short_ct <- function(x) {
  dplyr::recode(x,
    "Ventricular Cardiomyocytes" = "Ventricular CM",
    "Cardiomyocytes (stressed)" = "CM (stressed)",
    "Pericytes - Smooth muscle cells" = "Pericytes / SMC",
    "Endothelial cells" = "Endothelial",
    "Lymphatic endothelial" = "Lymphatic EC",
    "Epicardial - Mesothelial cells" = "Epicardial")
}
CT_ORDER <- c("Ventricular CM","Fibroblasts","Endothelial","Macrophages","Pericytes / SMC",
              "Endocardium","Lymphatic EC","B cells","T cells","CM (stressed)","Epicardial")
theme_set(theme_classic(base_size = 15) +
          theme(strip.background = element_blank(),
                strip.text = element_text(face = "bold", size = 13),
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12, colour = "grey30"),
                legend.title = element_text(size = 12)))
AGE <- c("9w", "78w")          # the slides show adult vs aged only
STRESS <- c("Sham", "TAC")     # Sham/TAC variants go to a side folder
OUT_ST <- file.path(OUT, "sham_tac"); dir.create(OUT_ST, showWarnings = FALSE)
save_fig <- function(p, name, w, h, dir = OUT) {
  f <- file.path(dir, paste0(name, ".png"))
  ggsave(f, p, width = w, height = h, dpi = 300, bg = "white", device = ragg::agg_png)
  message("wrote ", basename(f))
}
sample_factor <- function(x, levels = c("9w","78w","Sham","TAC")) factor(x, levels = levels)

# ===========================================================================
# PART 1 - snRNA-seq (OCM: 9w, 78w, Sham, TAC)
# ===========================================================================
cells <- read_tsv(file.path(OCM, "Allelic_ratio_results/cutoff_sweep/cutoff_sweep_cell_table.txt"),
                  show_col_types = FALSE) %>%
  mutate(sample = sample_factor(sample), celltype = short_ct(celltype))
counts_all <- read_tsv(file.path(OCM, "cell_counts_per_celltype_and_condition.txt"), show_col_types = FALSE) %>%
  mutate(sample = sample_factor(sample), celltype = factor(short_ct(celltype), CT_ORDER))

# F00 - nuclei per cell type and sample -------------------------------------
p <- ggplot(counts_all %>% filter(!celltype %in% c("CM (stressed)","Epicardial")),
            aes(celltype, n_cells, fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
  labs(x = NULL, y = "Nuclei", title = "Nuclei per cell type after QC",
       subtitle = "One animal per condition") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
save_fig(p, "F00_snRNA_nuclei_per_celltype", 9, 5)

# F01 - informative chrX reads per nucleus ------------------------------------
f01 <- function(samples, dir) {
  d <- cells %>% filter(sample %in% samples) %>% mutate(sample = sample_factor(sample, samples))
  ann <- d %>% group_by(sample) %>%
    summarise(n = n(), kept = sum(total_reads >= MIN_TOTAL_READS), .groups = "drop") %>%
    mutate(lab = sprintf("%s of %s nuclei\n>= %d reads", comma(kept), comma(n), MIN_TOTAL_READS))
  p <- ggplot(d, aes(total_reads, fill = sample)) +
    geom_histogram(bins = 45, colour = NA) +
    geom_vline(xintercept = MIN_TOTAL_READS, linetype = 2) +
    geom_text(data = ann, aes(x = 2000, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = 1, vjust = 1.3, size = 3.8) +
    scale_x_log10(labels = comma) +
    scale_fill_manual(values = SAMPLE_COL, guide = "none") +
    facet_wrap(~sample, labeller = labeller(sample = SAMPLE_LAB), scales = "free_y") +
    labs(x = "SNP-overlapping chrX reads per nucleus (log scale)", y = "Nuclei",
         title = "Allelic depth on chrX per nucleus",
         subtitle = sprintf("Dashed line: %d-read cutoff used for per-nucleus allelic ratios", MIN_TOTAL_READS))
  save_fig(p, "F01_snRNA_chrX_reads_per_nucleus", 9, 4.2, dir = dir)
}
f01(AGE, OUT); f01(STRESS, OUT_ST)

# F02 - autosomes vs chrX: the model works -----------------------------------
wc <- read_tsv(file.path(OCM, "Allelic_ratio_results/whole_chr_allelic_ratios.txt"), show_col_types = FALSE)
wc2 <- wc %>%
  mutate(set = ifelse(chr == "chrX", "chrX", ifelse(chr %in% paste0("chr", 1:19), "Autosomes", NA))) %>%
  filter(!is.na(set)) %>%
  group_by(cell_barcode, set) %>%
  summarise(A1 = sum(A1_reads), A2 = sum(A2_reads), .groups = "drop") %>%
  mutate(total = A1 + A2, ar = A1 / total, sample = sub("_.*", "", cell_barcode)) %>%
  filter(total >= MIN_TOTAL_READS)
f02 <- function(samples, dir) {
  d <- wc2 %>% filter(sample %in% samples) %>% mutate(sample = sample_factor(sample, samples))
  ann2 <- d %>% group_by(sample, set) %>%
    summarise(pooled = sum(A1) / sum(total), n = n(), .groups = "drop") %>%
    mutate(lab = sprintf("pooled %.2f\nn = %s", pooled, comma(n)))
  p <- ggplot(d, aes(sample, ar, fill = sample)) +
    geom_violin(scale = "width", bounds = c(0, 1), colour = NA, alpha = 0.9) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    geom_hline(yintercept = 0.5, linetype = 3, colour = "grey40") +
    geom_text(data = ann2, aes(sample, -0.02, label = lab), inherit.aes = FALSE, size = 3.3, vjust = 1) +
    facet_wrap(~set) +
    scale_fill_manual(values = SAMPLE_COL, guide = "none") +
    scale_x_discrete(labels = SAMPLE_LAB) +
    coord_cartesian(ylim = c(-0.12, 1.02)) +
    labs(x = NULL, y = "B6 fraction per nucleus  (A1 / total)",
         title = "Per-nucleus allelic ratio: autosomes vs chrX",
         subtitle = "Autosomes near 0.5 (slight B6 mapping bias); chrX B6-dominated because the CAST X is inactive")
  save_fig(p, "F02_snRNA_autosome_vs_chrX", 8, 4.2, dir = dir)
}
f02(AGE, OUT); f02(STRESS, OUT_ST)

# F03 - whole-chrX AR by cell type and condition ------------------------------
meta <- read_tsv(file.path(OCM, "Allelic_ratio_results/cutoff_30/whole_chr_cell_metadata.txt"),
                 show_col_types = FALSE) %>%
  mutate(celltype = factor(short_ct(celltype), CT_ORDER)) %>%
  filter(!celltype %in% c("CM (stressed)", "Epicardial"))
f03 <- function(samples, dir) {
  d <- meta %>% filter(sample %in% samples) %>% mutate(sample = sample_factor(sample, samples))
  p <- ggplot(d, aes(sample, allelic_ratio, fill = sample)) +
    geom_violin(scale = "width", bounds = c(0, 1), colour = NA, alpha = 0.9) +
    geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white") +
    geom_hline(yintercept = MONO_AR, linetype = 2, colour = "grey30") +
    facet_wrap(~celltype, nrow = 2) +
    scale_fill_manual(values = SAMPLE_COL, guide = "none") +
    scale_x_discrete(labels = SAMPLE_LAB) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1)) +
    labs(x = NULL, y = "Whole-chrX allelic ratio per nucleus",
         title = "Whole-chromosome allelic ratio per nucleus",
         subtitle = sprintf("Nuclei with >= %d SNP-overlapping chrX reads; dashed line = %.1f monoallelic boundary", MIN_TOTAL_READS, MONO_AR)) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save_fig(p, "F03_snRNA_chrX_AR_violin_by_celltype", 12, 6.5, dir = dir)
}
f03(AGE, OUT); f03(STRESS, OUT_ST)

# F04 - fraction of nuclei below the monoallelic boundary --------------------
fe <- read_tsv(file.path(OCM, "Allelic_ratio_results/cutoff_30/whole_chr_fraction_escaping_per_celltype_and_condition.txt"),
               show_col_types = FALSE) %>%
  mutate(celltype = factor(short_ct(celltype), CT_ORDER)) %>%
  filter(!celltype %in% c("CM (stressed)", "Epicardial")) %>%
  mutate(se = sqrt(escaping * (1 - escaping) / n))
f04 <- function(samples, dir, subtitle) {
  d <- fe %>% filter(sample %in% samples) %>% mutate(sample = sample_factor(sample, samples))
  p <- ggplot(d, aes(celltype, 100 * escaping, fill = sample)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_errorbar(aes(ymin = 100 * (escaping - se), ymax = 100 * (escaping + se)),
                  position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.3) +
    scale_fill_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
    labs(x = NULL, y = sprintf("Nuclei with chrX AR < %.1f (%%)", MONO_AR),
         title = "Nuclei with detectable biallelic chrX expression", subtitle = subtitle) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
  save_fig(p, "F04_snRNA_fraction_biallelic_nuclei", 10, 5.5, dir = dir)
}
f04(AGE, OUT, "Error bars: binomial SE. One animal per age - descriptive only")
f04(STRESS, OUT_ST, "Error bars: binomial SE. One animal per condition - descriptive only")

# F05 - depth bias: mean AR rises with allelic depth --------------------------
bands <- c(10, 25, 40, 60, 100, 200, Inf)
f05 <- function(samples, dir) {
  db <- cells %>% filter(total_reads >= 10, sample %in% samples) %>%
    mutate(sample = sample_factor(sample, samples),
           band = cut(total_reads, bands, right = FALSE,
                      labels = c("10-25","25-40","40-60","60-100","100-200","200+"))) %>%
    group_by(sample, band) %>%
    summarise(n = n(), pooled = sum(A1_reads) / sum(total_reads),
              frac_bi = mean(allelic_ratio < MONO_AR), .groups = "drop")
  p1 <- ggplot(db, aes(band, pooled, colour = sample, group = sample)) +
    geom_line(linewidth = 1) + geom_point(aes(size = n)) +
    scale_colour_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
    scale_size_area(max_size = 6, name = "Nuclei") +
    labs(x = "SNP-overlapping chrX reads per nucleus", y = "Pooled B6 fraction",
         title = "Apparent escape depends on allelic depth")
  p2 <- ggplot(db, aes(band, 100 * frac_bi, colour = sample, group = sample)) +
    geom_line(linewidth = 1) + geom_point(aes(size = n), show.legend = FALSE) +
    scale_colour_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
    scale_size_area(max_size = 6) +
    labs(x = "SNP-overlapping chrX reads per nucleus", y = sprintf("Nuclei with AR < %.1f (%%)", MONO_AR),
         title = "Shallow nuclei look more biallelic") + guides(colour = "none")
  p <- (p1 | p2) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p, "F05_snRNA_depth_bias", 12, 5.5, dir = dir)
}
f05(AGE, OUT); f05(STRESS, OUT_ST)

# F06 - UMAP coloured by whole-chrX AR, per sample ----------------------------
f06 <- function(samples, dir, w) {
  um <- cells %>% filter(total_reads >= MIN_TOTAL_READS, sample %in% samples) %>%
    mutate(bin = ar_bin(allelic_ratio), sample = sample_factor(sample, samples))
  p <- ggplot(um, aes(UMAP_1, UMAP_2, colour = bin)) +
    geom_point(size = 0.5, alpha = 0.9) +
    scale_colour_manual(values = setNames(AR_COLS, AR_LABELS), drop = FALSE, name = "chrX allelic\nratio",
                        guide = guide_legend(override.aes = list(size = 3))) +
    facet_wrap(~sample, labeller = labeller(sample = SAMPLE_LAB), nrow = 1) +
    coord_equal() +
    labs(title = "Whole-chrX allelic ratio per nucleus on the UMAP",
         subtitle = sprintf("Nuclei with >= %d SNP-overlapping chrX reads", MIN_TOTAL_READS)) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  save_fig(p, "F06_snRNA_AR_umap_by_sample", w, 5.2, dir = dir)
}
f06(AGE, OUT, 11); f06(STRESS, OUT_ST, 11)

# F07 - core escape genes: posterior mean AR per cell type ---------------------
bp <- read_tsv(file.path(OCM, "Allelic_ratio_results/core_escape_genes_bayes_posterior_by_gene.txt"),
               show_col_types = FALSE) %>%
  mutate(celltype = factor(short_ct(celltype), rev(CT_ORDER))) %>%
  filter(!celltype %in% c("CM (stressed)", "Epicardial"), total_reads >= 10)
f07 <- function(samples, dir) {
  d <- bp %>% filter(sample %in% samples) %>% mutate(sample = sample_factor(sample, samples))
  p <- ggplot(d, aes(sample, celltype, fill = post_mean, size = total_reads)) +
    geom_point(shape = 21, colour = "grey20") +
    facet_wrap(~name, nrow = 1) +
    ar_fill_cont("Posterior mean\nallelic ratio") +
    scale_size_area(max_size = 10, name = "Informative\nreads") +
    scale_x_discrete(labels = SAMPLE_LAB) +
    labs(x = NULL, y = NULL, title = "Core escape genes: per cell type allelic ratio",
         subtitle = "Beta-binomial posterior (prior from Hoelzl et al. 2025 bulk medians); reads pooled across nuclei of a cell type") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save_fig(p, "F07_snRNA_core_escape_posterior", 12, 5.5, dir = dir)
}
f07(AGE, OUT); f07(STRESS, OUT_ST)

# F20 - Xist vs allelic ratio, beta-binomial odds ratio per cell type ----------
xb <- read_tsv(file.path(OCM, "Allelic_ratio_results/core_escape_cutoff_5/core_escape_block_new_Xist_vs_AR_betabinomial.txt"),
               show_col_types = FALSE) %>%
  mutate(celltype = factor(short_ct(celltype), rev(CT_ORDER)),
         lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se), sig = FDR < 0.05)
f20 <- function(samples, dir) {
  d <- xb %>% filter(sample %in% samples) %>% mutate(sample = sample_factor(sample, samples))
  p <- ggplot(d, aes(OR, celltype, colour = sig)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25) +
    geom_point(aes(size = n_cells)) +
    facet_wrap(~sample, labeller = labeller(sample = SAMPLE_LAB)) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.1, 6)) +
    scale_colour_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "grey45"), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "ns"), name = NULL) +
    scale_size_area(max_size = 6, name = "Nuclei") +
    labs(x = "Odds ratio per unit Xist (95% CI)", y = NULL,
         title = "More Xist, more monoallelic",
         subtitle = "Beta-binomial: core-escape-block allelic ratio ~ Xist. OR < 1: higher Xist, more monoallelic")
  save_fig(p, "F20_snRNA_Xist_vs_AR_forest", 10, 5, dir = dir)
}
f20(AGE, OUT); f20(STRESS, OUT_ST)

# F08 - core escape block, adult vs aged per cell type ------------------------
pb <- read_tsv(file.path(OCM, "Allelic_ratio_results/core_escape_genes_pseudobulk_by_celltype.txt"),
               show_col_types = FALSE) %>%
  mutate(celltype = short_ct(celltype), cast = A2_reads / total_reads,
         se = sqrt(cast * (1 - cast) / total_reads))
ab <- read_tsv(file.path(OCM, "Allelic_ratio_results/core_escape_genes_pseudobulk_9w_vs_78w.txt"),
               show_col_types = FALSE) %>%
  mutate(celltype = short_ct(celltype), star = ifelse(FDR <= 0.001, "***", ifelse(FDR <= 0.01, "**", ifelse(FDR <= 0.05, "*", "ns"))))
pb2 <- pb %>% filter(sample %in% c("9w", "78w"), total_reads >= 30) %>%
  mutate(sample = sample_factor(sample, c("9w","78w")), celltype = factor(celltype, rev(CT_ORDER)))
ab2 <- ab %>% mutate(celltype = factor(celltype, rev(CT_ORDER))) %>%
  inner_join(pb2 %>% group_by(celltype) %>% summarise(xmax = max(cast + se), .groups = "drop"), by = "celltype")
p <- ggplot(pb2, aes(cast, celltype, colour = sample)) +
  geom_line(aes(group = celltype), colour = "grey60", linewidth = 1) +
  geom_errorbarh(aes(xmin = cast - se, xmax = cast + se), height = 0.2) +
  geom_point(size = 4) +
  geom_text(data = ab2, aes(x = xmax + 0.03, y = celltype, label = star), inherit.aes = FALSE, size = 4) +
  scale_colour_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
  scale_x_continuous(labels = percent, limits = c(0, 0.62)) +
  labs(x = "CAST (inactive X) fraction of pooled reads", y = NULL,
       title = "Core escape genes: inactive-X expression, adult vs aged",
       subtitle = "Pooled over Kdm5c, Kdm6a, Ddx3x, Eif2s3x; Fisher test FDR; n = 1 animal per age") +
  theme(legend.position = "top")
save_fig(p, "F08_snRNA_core_escape_adult_vs_aged", 8.5, 5.5)

# F09 - per-gene pseudobulk heatmap (grant Figure 2 style) --------------------
gf <- file.path(DAT, "all_genes_per_cell.tsv")
logf <- file.path(DAT, "consolidate.log")
consolidated <- file.exists(gf) && file.exists(logf) && any(grepl("^FINISHED", readLines(logf)))
if (consolidated) {
  g <- read_tsv(gf, show_col_types = FALSE, col_types = cols(.default = col_character(),
                A1_reads = col_double(), A2_reads = col_double(), total_reads = col_double(), start = col_double()))
  g <- g %>% inner_join(cells %>% select(cell_barcode, celltype), by = "cell_barcode")
  pbg <- g %>% group_by(sample, celltype, gene, start) %>%
    summarise(A1 = sum(A1_reads), A2 = sum(A2_reads), .groups = "drop") %>%
    mutate(total = A1 + A2, ar = A1 / total) %>%
    filter(total >= 20, !celltype %in% c("CM (stressed)", "Epicardial"))
  make_heat <- function(samples, name, title, min_ct = 6, landscape = TRUE) {
    d <- pbg %>% filter(sample %in% samples)
    # slide-sized: genes testable (>= 20 reads) in at least min_ct cell types in EACH sample
    keep <- d %>% count(sample, gene) %>% filter(n >= min_ct) %>% count(gene) %>%
      filter(n == length(samples)) %>% pull(gene)
    d <- d %>% filter(gene %in% keep) %>%
      mutate(gene = fct_reorder(gene, start), celltype = factor(celltype, CT_ORDER),
             sample = sample_factor(sample, samples))
    ng <- n_distinct(d$gene); message(name, ": ", ng, " genes")
    sub <- "Pseudobulk per cell type, >= 20 SNP-overlapping reads; green/blue = expression from both X, red = monoallelic"
    if (landscape) {
      p <- ggplot(d, aes(gene, celltype, fill = ar)) +
        geom_tile(colour = "white", linewidth = 0.3) +
        facet_wrap(~sample, ncol = 1, labeller = labeller(sample = SAMPLE_LAB)) +
        ar_fill_cont("Allelic ratio\n(B6 / total)") +
        scale_y_discrete(limits = rev) +
        labs(x = sprintf("chrX genes testable in >= %d cell types per age (left = distal Xp, ordered by position)", min_ct), y = NULL,
             title = title, subtitle = sub) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7,
                                         face = ifelse(levels(d$gene) %in% ESCAPE_GENES, "bold", "plain")),
              panel.grid = element_blank(), legend.position = "right")
      save_fig(p, name, 13, 6.3)
    } else {
      p <- ggplot(d, aes(celltype, gene, fill = ar)) +
        geom_tile(colour = "white", linewidth = 0.3) +
        facet_wrap(~sample, labeller = labeller(sample = SAMPLE_LAB)) +
        ar_fill_cont("Allelic ratio\n(B6 / total)") +
        scale_y_discrete(limits = rev) +
        labs(x = NULL, y = "chrX genes (ordered by position, top = distal Xp)", title = title, subtitle = sub) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7),
              panel.grid = element_blank())
      save_fig(p, name, 9, min(30, max(6, 0.17 * ng + 2.5)))
    }
    invisible(d)
  }
  make_heat(c("9w","78w"), "F09_snRNA_chrX_gene_heatmap_adult_aged", "chrX allelic ratio per gene and cell type: adult vs aged")
  make_heat(c("9w","78w"), "F09_snRNA_chrX_gene_heatmap_adult_aged_full", "chrX allelic ratio per gene and cell type: adult vs aged (all testable genes)", min_ct = 1, landscape = FALSE)
  make_heat(c("Sham","TAC"), "F09b_snRNA_chrX_gene_heatmap_sham_tac", "chrX allelic ratio per gene and cell type: Sham vs TAC")
  # count escaping genes per cell type (AR < MONO_AR) - simple summary
  esc <- pbg %>% filter(sample %in% AGE) %>% mutate(escape = ar < MONO_AR) %>% group_by(sample, celltype) %>%
    summarise(n_genes = n(), n_escape = sum(escape), .groups = "drop") %>%
    mutate(sample = sample_factor(sample, AGE), celltype = factor(celltype, CT_ORDER))
  write_tsv(esc, file.path(DAT, "escape_gene_counts_by_celltype.tsv"))
  p <- ggplot(esc, aes(celltype, n_escape, fill = sample)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(aes(label = n_genes, y = n_escape + 0.5, group = sample), position = position_dodge(width = 0.8), size = 2.8, colour = "grey30") +
    scale_fill_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
    labs(x = NULL, y = sprintf("Genes with pseudobulk AR < %.1f", MONO_AR),
         title = "Escaping genes per cell type", subtitle = "Grey number: genes testable (>= 20 reads) in that cell type") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
  save_fig(p, "F09c_snRNA_escape_gene_counts", 10, 5.5)
} else message("SKIP F09: consolidation of per-nucleus gene tables not finished yet")

# F21/F22 - cell-type UMAP and QC panels for slide 6 (all nuclei, from the Seurat metadata dump)
# Preferred input: the cluster dump of the Seurat object (all nuclei). Fallback
# until that exists: UMAP + counts from cutoff_sweep_cell_table (every nucleus
# with chrX coverage, i.e. essentially all) and QC metrics from the cutoff_30
# metadata (86-89% of nuclei). The figure subtitle says which was used.
smf <- file.path(DAT, "seurat_metadata_umap.tsv")
if (file.exists(smf)) {
  sm <- read_tsv(smf, show_col_types = FALSE) %>%
    mutate(celltype = short_ct(celltype), sample = as.character(sample))
  qc_note <- "All nuclei after QC"
} else {
  qc_meta <- read_tsv(file.path(OCM, "Allelic_ratio_results/cutoff_30/whole_chr_cell_metadata.txt"), show_col_types = FALSE) %>%
    rename(cell_barcode = 1) %>% select(cell_barcode, nFeature_RNA, percent.mt)
  sm <- cells %>% select(cell_barcode, sample, celltype, nCount_RNA, UMAP_1, UMAP_2) %>%
    mutate(sample = as.character(sample)) %>% left_join(qc_meta, by = "cell_barcode")
  qc_note <- "Features and mito %: nuclei with >= 30 chrX reads (86-89%); UMAP and proportions: all nuclei"
}
{
  CT_COL <- c("Ventricular CM" = "#D9A21B", "Fibroblasts" = "#E8735A", "Endothelial" = "#A8A41F",
              "Macrophages" = "#6BB33B", "Pericytes / SMC" = "#22B573", "Endocardium" = "#1FB7B7",
              "Lymphatic EC" = "#2E9FE0", "B cells" = "#4B6FE0", "T cells" = "#9B6BE0",
              "CM (stressed)" = "#E05BCB", "Epicardial" = "#F06E9C")
  f21 <- function(samples, dir) {
    d <- sm %>% filter(sample %in% samples)
    lab <- d %>% group_by(celltype) %>% summarise(x = median(UMAP_1), y = median(UMAP_2), n = n(), .groups = "drop") %>% filter(n >= 20)
    p <- ggplot(d, aes(UMAP_1, UMAP_2, colour = celltype)) +
      geom_point(size = 0.35, alpha = 0.8) +
      geom_label_repel(data = lab, aes(x, y, label = celltype), inherit.aes = FALSE, size = 3.6,
                       label.size = 0, fill = alpha("white", 0.75), seed = 1) +
      scale_colour_manual(values = CT_COL, guide = "none") + coord_equal() +
      labs(x = "UMAP 1", y = "UMAP 2", title = sprintf("%s nuclei, adult and aged", comma(nrow(d)))) +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
    save_fig(p, "F21_snRNA_UMAP_celltypes", 6.5, 6.5, dir = dir)
  }
  f22 <- function(samples, dir) {
    d <- sm %>% filter(sample %in% samples) %>% mutate(sample = sample_factor(sample, samples))
    vln <- function(col, ttl) ggplot(d, aes(sample, .data[[col]], fill = sample)) +
      geom_violin(scale = "width", colour = NA) + geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
      scale_fill_manual(values = SAMPLE_COL, guide = "none") + scale_x_discrete(labels = SAMPLE_LAB) +
      labs(x = NULL, y = NULL, title = ttl)
    prop <- d %>% count(sample, celltype) %>% group_by(sample) %>% mutate(p = n / sum(n)) %>% ungroup() %>%
      mutate(celltype = factor(celltype, names(CT_COL)))
    pc <- ggplot(prop, aes(sample, p, fill = celltype)) + geom_col(width = 0.7) +
      scale_fill_manual(values = CT_COL, name = NULL) + scale_x_discrete(labels = SAMPLE_LAB) +
      scale_y_continuous(labels = percent, expand = c(0, 0)) +
      labs(x = NULL, y = NULL, title = "Cell type proportion") + theme(legend.text = element_text(size = 9))
    p <- (vln("nCount_RNA", "UMI count") | vln("nFeature_RNA", "Number of features")) /
         (vln("percent.mt", "Mitochondrial %") | pc) +
      plot_annotation(caption = qc_note, theme = theme(plot.caption = element_text(size = 9, colour = "grey40")))
    save_fig(p, "F22_snRNA_QC_panels", 8, 6.5, dir = dir)
  }
  f21(AGE, OUT); f22(AGE, OUT); f21(STRESS, OUT_ST); f22(STRESS, OUT_ST)
}

# ===========================================================================
# PART 2 - Visium HD spatial (9w, 78w)
# ===========================================================================
comp <- read_csv(file.path(SPA, "annotation/composition_all_samples_square_008um.csv"), show_col_types = FALSE) %>%
  mutate(sample = sample_factor(sample, c("9w","78w")))
ct_lv <- comp %>% filter(sample == "9w") %>% arrange(desc(pct)) %>% pull(celltype)
comp <- comp %>% mutate(celltype = factor(celltype, ct_lv))
p1 <- ggplot(comp, aes(sample, pct, fill = celltype)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c(setNames(c("#E3B769","#C98A2E","#1B9E77","#1F78B4","#A6CEE3","#B2DF8A","#E7298A","#D95F02","#7570B3","#66A61E","#E6AB02","#666666","#D9D9D9"), ct_lv)), name = NULL) +
  scale_x_discrete(labels = SAMPLE_LAB) +
  labs(x = NULL, y = "% of QC-passing 8 um bins", title = "Composition per section")
p2 <- ggplot(comp %>% filter(!celltype %in% c("Ventricular cardiomyocyte","Unassigned")),
             aes(celltype, pct, fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
  labs(x = NULL, y = "% of bins", title = "Non-ventricular-myocyte labels") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
save_fig(p1 + p2 + plot_layout(widths = c(1, 1.6)), "F10_spatial_composition", 13, 5.5)

# F11 - tile size vs precision -----------------------------------------------
prec <- bind_rows(lapply(c("9w","78w"), function(s)
  read_csv(file.path(SPA, "ase", s, "escape_precision.csv"), show_col_types = FALSE) %>% mutate(sample = s))) %>%
  mutate(sample = sample_factor(sample, c("9w","78w")))
p1 <- prec %>% select(sample, size_um, frac_ge_10, frac_ge_20) %>%
  pivot_longer(starts_with("frac"), names_to = "thr", values_to = "frac") %>%
  mutate(thr = ifelse(thr == "frac_ge_10", ">= 10 UMIs", ">= 20 UMIs")) %>%
  ggplot(aes(size_um, frac, colour = sample, linetype = thr)) +
  geom_vline(xintercept = 64, colour = "grey60", linetype = 3) +
  geom_line(linewidth = 1) + geom_point() +
  scale_x_log10(breaks = c(2,8,16,32,64,128,256,512)) +
  scale_colour_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
  scale_linetype(name = "chrX informative\nUMIs per tile") +
  labs(x = "Tile size (um)", y = "Fraction of tissue tiles", title = "Coverage: tiles reaching usable allelic depth")
p2 <- prec %>% filter(!is.na(se_escape)) %>%
  select(sample, size_um, se_escape, mde) %>%
  pivot_longer(c(se_escape, mde), names_to = "stat", values_to = "v") %>%
  mutate(stat = ifelse(stat == "mde", "Minimum detectable difference", "SE of escape fraction")) %>%
  ggplot(aes(size_um, v, colour = sample, linetype = stat)) +
  geom_vline(xintercept = 64, colour = "grey60", linetype = 3) +
  geom_line(linewidth = 1) + geom_point() +
  scale_x_log10(breaks = c(16,32,64,128,256,512)) +
  scale_colour_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
  scale_linetype(name = NULL) +
  labs(x = "Tile size (um)", y = "Per-tile precision (escape fraction)", title = "Precision: what a single tile can resolve")
save_fig((p1 | p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(size = 10)),
         "F11_spatial_tile_precision", 13, 6)

# F12 - 64 um tile maps: chrX vs autosomes ------------------------------------
# Tiles are placed on their row/column indices (64 um each), which is exact;
# the x/y columns are image coordinates whose spacing drifts.
tm <- read_csv(file.path(SPA, "ase/tile_ratio_map_64um.csv"), show_col_types = FALSE) %>%
  mutate(sample = sample_factor(sample, c("9w","78w")),
         row = as.integer(sub(".*_r(\\d+)_c\\d+$", "\\1", tile)),
         col = as.integer(sub(".*_c(\\d+)$", "\\1", tile)),
         x_cast = x_a2 / x_n, a_cast = a_a2 / a_n)
tml <- tm %>% filter(!is.na(x_n), x_n >= 10) %>%
  select(sample, row, col, `chrX (inactive-X escape)` = x_cast, `Autosomes (control)` = a_cast) %>%
  pivot_longer(c(`chrX (inactive-X escape)`, `Autosomes (control)`), names_to = "set", values_to = "cast") %>%
  mutate(set = factor(set, c("chrX (inactive-X escape)", "Autosomes (control)")))
tissue <- tm %>% select(sample, row, col)
p <- ggplot(tml) +
  geom_tile(data = tissue, aes(col, row), fill = "grey88", width = 1, height = 1) +
  geom_tile(aes(col, row, fill = cast), width = 1, height = 1) +
  facet_grid(set ~ sample, labeller = labeller(sample = SAMPLE_LAB)) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 0.6), oob = squish, name = "CAST fraction\n(0 = fully silenced)") +
  scale_y_reverse() + coord_equal() +
  labs(title = "Allelic ratio per 64 um tile", x = NULL, y = NULL,
       subtitle = "Tiles with >= 10 informative UMIs; grey = tissue tiles below depth. Autosomes sit near 0.5, chrX near 0.1 everywhere") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
save_fig(p, "F12_spatial_tile_maps_64um", 11, 10)

# F13 - tile distributions ----------------------------------------------------
p <- ggplot(tml, aes(cast, fill = set)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.7, colour = NA) +
  facet_wrap(~sample, labeller = labeller(sample = SAMPLE_LAB), scales = "free_y") +
  scale_fill_manual(values = c("#7B3294", "#B8B8B8"), name = NULL) +
  labs(x = "CAST fraction per 64 um tile", y = "Tiles", title = "Per-tile CAST fraction: chrX vs autosomes") +
  theme(legend.position = "top")
save_fig(p, "F13_spatial_tile_distribution", 10, 4.5)

# F14 - pair correlation vs distance: no spatial structure --------------------
pc <- bind_rows(lapply(c("9w","78w"), function(s)
  read_csv(file.path(SPA, "ase", s, "pair_correlation.csv"), show_col_types = FALSE) %>% mutate(sample = s))) %>%
  filter(chrom_set %in% c("chrX", "autosome", "imppat", "nonescape")) %>%
  mutate(sample = sample_factor(sample, c("9w","78w")),
         chrom_set = recode(chrom_set, chrX = "chrX", autosome = "Autosomes",
                            imppat = "Imprinted (paternal), monoallelic control", nonescape = "chrX minus core escape genes"),
         chrom_set = factor(chrom_set, c("chrX", "chrX minus core escape genes", "Autosomes", "Imprinted (paternal), monoallelic control")))
pnull <- data.frame(yint = c(0.8728^2 + (1 - 0.8728)^2, 0.5),
                    lab = c("no-structure expectation for chrX (p = 0.873)", "biallelic expectation"))
p1 <- ggplot(pc, aes(dist_um, C, colour = chrom_set)) +
  geom_hline(data = pnull, aes(yintercept = yint), linetype = 3, colour = "grey50") +
  geom_ribbon(aes(ymin = C - 1.96 * se, ymax = C + 1.96 * se, fill = chrom_set), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~sample, labeller = labeller(sample = SAMPLE_LAB)) +
  scale_x_log10(breaks = c(4, 16, 64, 256, 1000, 2000)) +
  scale_colour_manual(values = c("#7B3294", "#C2A5CF", "#8C8C8C", "#E2711D"), name = NULL) +
  scale_fill_manual(values = c("#7B3294", "#C2A5CF", "#8C8C8C", "#E2711D"), name = NULL) +
  labs(x = "Distance between UMI pairs (um)", y = "P(same allele)  C(d)",
       title = "Do neighbouring molecules share an allele more than distant ones?",
       subtitle = "Flat C(d) from 4 um to 2 mm = no patches of escape; imprinted loci show the assay resolves monoallelic expression") +
  theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
save_fig(p1, "F14_spatial_pair_correlation", 12, 6)

# F15 - imprinted-locus controls ------------------------------------------------
imp <- bind_rows(lapply(c("9w","78w"), function(s)
  read_tsv(file.path(SPA, "ase", s, "subset_locus_counts.tsv"), show_col_types = FALSE) %>% mutate(sample = s))) %>%
  filter(subset %in% c("imppat", "impmat"), n >= 10) %>%
  mutate(sample = sample_factor(sample, c("9w","78w")),
         expected = ifelse(subset == "imppat", "Paternally expressed (expect CAST)", "Maternally expressed (expect B6)"),
         cast_frac = alt / n,
         excluded = locus %in% c("Cdkn1c", "Mest", "Impact"),
         locus = fct_reorder(locus, cast_frac))
p <- ggplot(imp, aes(locus, cast_frac, fill = expected, alpha = !excluded)) +
  geom_col() +
  geom_text(aes(label = n, y = cast_frac + 0.03), size = 3.2, alpha = 1) +
  facet_wrap(~sample, labeller = labeller(sample = SAMPLE_LAB)) +
  scale_fill_manual(values = c("#D95F02", "#1B9E77"), name = NULL) +
  scale_alpha_manual(values = c(0.35, 1), guide = "none") +
  scale_y_continuous(labels = percent, limits = c(0, 1.08)) +
  labs(x = NULL, y = "CAST fraction of informative UMIs",
       title = "Imprinted loci confirm the cross and the allele calls",
       subtitle = "B6 mother x CAST father: paternal loci read CAST, maternal loci read B6. Faded = loci biallelic in heart, excluded. Numbers = UMIs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_fig(p, "F15_spatial_imprinted_controls", 12, 6)

# F16 - CAST fraction along chrX in 100 kb windows ------------------------------
aw <- bind_rows(lapply(c("9w","78w"), function(s)
  read_tsv(file.path(SPA, "ase", s, sprintf("artifact_windows_%s.tsv", s)), show_col_types = FALSE) %>% mutate(sample = s))) %>%
  filter(is_x, n >= 20) %>%
  mutate(sample = sample_factor(sample, c("9w","78w")), mb = win_start / 1e6, cast = alt / n,
         genes = ifelse(is.na(genes) | genes == "", "no annotated gene", genes))
top <- aw %>% group_by(sample) %>% slice_max(alt, n = 6) %>% ungroup() %>%
  mutate(lab = sprintf("%.1f Mb: %s", mb, sub(",.*", "", genes)))
p <- ggplot(aw, aes(mb, cast)) +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey50") +
  geom_point(aes(size = n), alpha = 0.5, colour = "#2B7BBA") +
  geom_point(data = top, aes(size = n), colour = "#D7301F") +
  geom_text_repel(data = top, aes(label = lab), size = 3, max.overlaps = 20, colour = "#D7301F", box.padding = 0.5) +
  facet_wrap(~sample, ncol = 1, labeller = labeller(sample = SAMPLE_LAB)) +
  scale_size_area(max_size = 8, name = "Molecules\nin window") +
  scale_y_continuous(labels = percent) +
  labs(x = "chrX position (Mb, mm39)", y = "CAST fraction per 100 kb window",
       title = "Where the chrX CAST signal comes from",
       subtitle = "Red = windows carrying the most CAST molecules: mostly unannotated or multicopy loci reading ~100% CAST, i.e. allele-call artefacts, not escape")
save_fig(p, "F16_spatial_chrX_window_scan", 12, 7)

# F17 - partition of the chromosome-wide escape estimate -----------------------
# Values from spatial/NEXT_ANALYSIS.md, status 2026-09-03 (gene-body split of the same UMIs)
part <- tribble(
  ~set, ~`9w`, ~`78w`,
  "chrX, all molecules",               0.1272, 0.1264,
  "chrX, outside gene bodies",         0.4750, 0.4981,
  "chrX, inside gene bodies",          0.0612, 0.0613,
  "chrX genes, minus 7 impossible (>50%)", 0.0284, 0.0299,
  "Autosomes, inside gene bodies",     0.4845, 0.4821) %>%
  pivot_longer(-set, names_to = "sample", values_to = "cast") %>%
  mutate(sample = sample_factor(sample, c("9w","78w")), set = factor(set, rev(unique(set))))
p <- ggplot(part, aes(cast, set, fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = percent(cast, accuracy = 0.1), x = cast + 0.015, group = sample),
            position = position_dodge(width = 0.8), hjust = 0, size = 3.6) +
  scale_fill_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
  scale_x_continuous(labels = percent, limits = c(0, 0.62)) +
  labs(x = "CAST (inactive X) fraction of informative UMIs", y = NULL,
       title = "The 12.7% chromosome-wide 'escape' is mostly artefact",
       subtitle = "Non-genic chrX molecules are 16% of chrX but carry 60% of its CAST signal, at the autosomal value") +
  theme(legend.position = "top")
save_fig(p, "F17_spatial_escape_partition", 11, 5)

# F18 - per-gene escape from the spatial data (spASE scASE) --------------------
sc <- bind_rows(lapply(c("9w","78w"), function(s)
  read_tsv(file.path(SPA, "ase", s, sprintf("spase_scase_%s_16um.tsv", s)), show_col_types = FALSE,
           col_types = cols(.default = col_character())) %>% mutate(sample = s))) %>%
  mutate(across(c(umi, ref, alt, p, ci.low, ci.high), as.numeric)) %>%
  filter(is_x == "TRUE", converged == "TRUE", impossible != "TRUE", umi >= 100) %>%
  mutate(sample = sample_factor(sample, c("9w","78w")),
         core = gene %in% ESCAPE_GENES)
both <- sc %>% count(gene) %>% filter(n == 2) %>% pull(gene)
sc2 <- sc %>% filter(gene %in% both) %>%
  mutate(gene = fct_reorder(gene, p, .fun = mean))
p <- ggplot(sc2, aes(p, gene, colour = sample)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_errorbarh(aes(xmin = ci.low, xmax = ci.high), height = 0, position = position_dodge(width = 0.6)) +
  geom_point(aes(size = umi, shape = core), position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), labels = c(`TRUE` = "known escapee", `FALSE` = "other chrX gene"), name = NULL) +
  scale_size_area(max_size = 6, name = "Informative\nUMIs") +
  scale_x_continuous(labels = percent) +
  labs(x = "Escape fraction  (CAST / total, 95% CI)", y = NULL,
       title = "Per-gene escape from the inactive X in the spatial data",
       subtitle = "chrX genes with >= 100 informative UMIs in both sections; 16 um pixels, binomial CIs") +
  theme(axis.text.y = element_text(face = ifelse(levels(sc2$gene) %in% ESCAPE_GENES, "bold", "plain")),
        legend.position = "right")
save_fig(p, "F18_spatial_per_gene_escape", 10, 8)

# ===========================================================================
# PART 3 - Multiome (9w, 78w)
# ===========================================================================
nj <- read_csv(file.path(MUL, "figures_signac/nucleus_celltype_joint.csv"), show_col_types = FALSE) %>%
  mutate(sample = sample_factor(sample, c("9w","78w")), celltype = factor(short_ct(celltype_provisional), CT_ORDER)) %>%
  count(sample, celltype)
p <- ggplot(nj, aes(celltype, n, fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = n, group = sample), position = position_dodge(width = 0.8), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = SAMPLE_COL, labels = SAMPLE_LAB, name = NULL) +
  labs(x = NULL, y = "Nuclei passing joint QC", title = "Multiome nuclei per cell type",
       subtitle = "WNN clustering on RNA + ATAC; labels transferred from the snRNA-seq marker panels") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
save_fig(p, "F19_multiome_nuclei_per_celltype", 9, 5)

message("all done")
