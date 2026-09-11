# ---------------------------------------------------------------------------
# chrX allelic ratio per gene and cell type, adult (9w) vs aged (78w), from the
# Allelome.PRO2 pseudobulk in OCM/Allelome.PRO2_pseudobulk_celltype.
#
# Allelome.PRO2 already pools one BAM per cell type per group, so
# locus_table.txt IS the pseudobulk and no per-cell summing is needed. The
# genotype, the allele coding and every cutoff are documented in
# ap2_pseudobulk_common.R, which this script sources.
#
# Gene labels are coloured by whether the gene sits in the distal window, and
# each figure carries the distal-enrichment statistics in its caption, so the
# heatmaps and chrX_distal_escape_enrichment.R always report the same numbers.
#
#   Rscript OCM_heart/pseudobulk_celltype_chrX_heatmap.R
#
# Config is by environment variable, see ap2_pseudobulk_common.R. For the
# stress arm: GROUPS=Sham,TAC Rscript OCM_heart/pseudobulk_celltype_chrX_heatmap.R
# ---------------------------------------------------------------------------
.this_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(sub("^--file=", "", f[1])) else "OCM_heart"
}
source(file.path(.this_dir(), "ap2_pseudobulk_common.R"))

BASE <- GROUPS[1]; COMP <- GROUPS[2]
# titles follow the contrast, so the Sham/TAC arm is not labelled "adult vs aged"
LAB <- function(g) sub(" \\(.*", "", SAMPLE_LAB[[g]])          # "Adult (9w)" -> "Adult"
CONTRAST_LAB <- sprintf("%s vs %s", LAB(BASE), LAB(COMP))
# the stress arm has ~3x as many testable genes as the age arm, so the figure
# has to widen with the gene count or the labels collide
fig_width <- function(ng) min(40, max(11, 4 + 0.115 * ng))

pbg <- load_ap2_chrX()
write_tsv(pbg %>% select(sample, celltype, gene, start, distal, A1, A2, total, ar),
          file.path(OUT, "chrX_pseudobulk_celltype_allelic_ratio.tsv"))
message(nrow(pbg), " gene x cell type x group observations with >= ", MIN_READS, " reads")

CAPTION <- distal_tests(pbg)$caption

# Gene labels: bold for known constitutive escapees, and optionally orange for
# genes in the distal window so the shading of the enrichment figure carries
# over. Set LABEL_DISTAL=0 for plain labels; the file names then end in _plain.
LABEL_DISTAL <- Sys.getenv("LABEL_DISTAL", "1") != "0"
SUFFIX       <- if (LABEL_DISTAL) "" else "_plain"
DISTAL_GENES <- unique(pbg$gene[pbg$distal])
label_cols <- function(levels_in_order) {
  if (LABEL_DISTAL) ifelse(levels_in_order %in% DISTAL_GENES, "#B5540F", "grey25") else "grey25"
}
gene_axis <- function(levels_in_order) {
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7,
                                   face = ifelse(levels_in_order %in% ESCAPE_GENES, "bold", "plain"),
                                   colour = label_cols(levels_in_order)),
        panel.grid = element_blank(), legend.position = "right")
}
SUB <- paste0(sprintf("Allelome.PRO2 pseudobulk per cell type, >= %d SNP-overlapping reads; green/blue = expression from both X, red = monoallelic",
                      MIN_READS),
              if (LABEL_DISTAL) sprintf(". Orange gene labels: %s", DISTAL_LAB) else "")

# ---- F1: genes testable in >= MIN_CT cell types in every group -------------
keep <- pbg %>% count(sample, gene) %>% filter(n >= MIN_CT) %>% count(gene) %>%
  filter(n == length(GROUPS)) %>% pull(gene)
d <- pbg %>% filter(gene %in% keep) %>% mutate(gene = fct_reorder(gene, start))
message(n_distinct(d$gene), " genes testable in >= ", MIN_CT, " cell types in every group")
p1 <- ggplot(d, aes(gene, celltype, fill = ar)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_wrap(~sample, ncol = 1, labeller = labeller(sample = SAMPLE_LAB)) +
  ar_fill_cont() + scale_y_discrete(limits = rev, drop = TRUE) +
  labs(x = sprintf("chrX genes testable in >= %d cell types per group (left = distal Xp, ordered by position)", MIN_CT),
       y = NULL, title = sprintf("chrX allelic ratio per gene and cell type: %s", CONTRAST_LAB),
       subtitle = SUB, caption = CAPTION) +
  gene_axis(levels(d$gene))
save_fig(p1, paste0("AP2_chrX_gene_heatmap", SUFFIX), fig_width(n_distinct(d$gene)), 7.4)

# ---- F2: same genes, change with age ---------------------------------------
delta <- paired_by_gene(pbg %>% filter(gene %in% keep)) %>%
  mutate(gene = fct_reorder(as.character(gene), start))
# a handful of large shifts would otherwise flatten everything else, so clamp
# the scale at the 98th percentile of |change| and squish the outliers onto it
lim <- max(0.1, unname(quantile(abs(delta$delta), 0.98, na.rm = TRUE)))
p2 <- ggplot(delta, aes(gene, celltype, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2B3186", mid = "grey95", high = "#8B1913", midpoint = 0,
                       limits = c(-lim, lim), oob = scales::squish,
                       name = sprintf("%s - %s\nallelic ratio", COMP, BASE), na.value = "white") +
  scale_y_discrete(limits = rev, drop = TRUE) +
  labs(x = "chrX genes measured in both groups (left = distal Xp, ordered by position)", y = NULL,
       title = sprintf("Change in chrX allelic ratio, %s, per gene and cell type", CONTRAST_LAB),
       subtitle = sprintf("Red: more monoallelic in %s. Blue: more expression from the inactive (CAST) X in %s",
                          LAB(COMP), LAB(COMP)),
       caption = CAPTION) +
  gene_axis(levels(delta$gene))
save_fig(p2, paste0("AP2_chrX_gene_heatmap_delta", SUFFIX), fig_width(n_distinct(delta$gene)), 5.7)

# ---- F3: every testable gene, portrait, for the record ---------------------
dall <- pbg %>% mutate(gene = fct_reorder(gene, start))
p3 <- ggplot(dall, aes(celltype, gene, fill = ar)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_wrap(~sample, labeller = labeller(sample = SAMPLE_LAB)) +
  ar_fill_cont() + scale_y_discrete(limits = rev) +
  labs(x = NULL, y = "chrX genes (ordered by position, top = distal Xp)",
       title = sprintf("chrX allelic ratio per gene and cell type: %s (all testable genes)", CONTRAST_LAB),
       subtitle = SUB, caption = CAPTION) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 4, colour = label_cols(rev(levels(dall$gene)))),
        panel.grid = element_blank())
save_fig(p3, paste0("AP2_chrX_gene_heatmap_full", SUFFIX), 9,
         min(40, max(6, 0.09 * n_distinct(dall$gene) + 2.5)))
