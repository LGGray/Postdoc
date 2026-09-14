# ---------------------------------------------------------------------------
# Adult, Sham, TAC and aged side by side, EVERY chrX gene, no filtering.
#
# chrX_four_group_comparison.R panel C answers "what do the comparable genes
# do": >= MIN_READS reads, testable in all four groups, and in >= MIN_CT cell
# types. That is the right gene set for a contrast and the wrong one for
# looking at the chromosome. This script draws the other picture - every chrX
# gene Allelome.PRO2 reports at all, in every cell type, in all four groups -
# so nothing is chosen by a cutoff.
#
#   Rscript OCM_heart/chrX_four_group_all_genes_heatmap.R
#
# Config as in ap2_pseudobulk_common.R, plus:
#   MIN_READS         1 here, not 20: every gene with at least one
#                     SNP-overlapping read. See the depth caveat below.
#   MAX_GENE_LABELS   150. Above this the x axis is unlabelled and carries
#                     Mb ticks instead - 883 gene names do not fit on a page.
#   DROP_CT_ALL       0 by default, so CM (stressed) and Epicardial are KEPT
#                     here even though the filtered scripts drop them.
#
# THE DEPTH CAVEAT, and why the bottom panel exists. At total_reads = 1 the
# allelic ratio can only be 0 or 1, at 2 only 0, 0.5 or 1. Removing the depth
# cutoff therefore does not reveal more monoallelic genes, it manufactures
# them: a tile at AR = 1 in the sparse part of the chromosome is a statement
# about coverage, not about escape. The per-gene depth track underneath is
# aligned to the same gene axis so any block of red can be read against the
# reads that produced it, and MIN_READS stays an environment variable so the
# same figure can be redrawn at 5, 10 or 20 to see what survives.
# ---------------------------------------------------------------------------
if (Sys.getenv("GROUPS")    == "") Sys.setenv(GROUPS = "9w,Sham,TAC,78w")
if (Sys.getenv("CONTRAST")  == "") Sys.setenv(CONTRAST = "four_group_all_genes")
if (Sys.getenv("MIN_READS") == "") Sys.setenv(MIN_READS = "1")
.this_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(sub("^--file=", "", f[1])) else "OCM_heart"
}
source(file.path(.this_dir(), "ap2_pseudobulk_common.R"))
suppressPackageStartupMessages(library(patchwork))

MAX_GENE_LABELS <- as.integer(Sys.getenv("MAX_GENE_LABELS", "150"))
KEEP_ALL_CT     <- Sys.getenv("DROP_CT_ALL", "0") == "0"

# No matched gene set, no MIN_CT, and by default no cell type dropped either.
pbg <- load_ap2_chrX(groups = GROUPS,
                     drop_ct = if (KEEP_ALL_CT) character(0) else DROP_CT) %>%
  mutate(sample = factor(sample, GROUPS), gene = fct_reorder(gene, start))
stopifnot(nrow(pbg) > 0)

n_gene <- n_distinct(pbg$gene)
message(n_gene, " chrX genes with >= ", MIN_READS, " read(s) in at least one ",
        "group x cell type; ", nrow(pbg), " tiles over ",
        n_distinct(pbg$celltype), " cell types and ", length(GROUPS), " groups")
write_tsv(pbg %>% select(sample, celltype, gene, start, distal, A1, A2, total, ar),
          file.path(OUT, "chrX_all_genes_allelic_ratio.tsv"))

# How much of the picture is depth rather than biology, stated as a number.
depth <- pbg %>% group_by(gene, start, distal) %>%
  summarise(total = sum(total), n_obs = n(), .groups = "drop")
pct_thin <- 100 * mean(pbg$total < 10)
CAPTION <- sprintf(
  "Every chrX gene Allelome.PRO2 reports at >= %d read(s): %d genes, %d gene x cell type x group observations, no gene or cell type excluded.\n%.0f%% of tiles carry fewer than 10 SNP-overlapping reads, where the ratio is coarse by construction (at 1 read it can only be 0 or 1) - read the tiles against the depth track below.\nDistal window: the first %g Mb plus the last %g Mb of chrX, marked by the dashed lines. One animal per group, so the groups differ descriptively only.",
  MIN_READS, n_gene, nrow(pbg), pct_thin, DISTAL_P_MB, DISTAL_Q_MB)
writeLines(CAPTION, file.path(OUT, "four_group_all_genes_caption.txt"))

# ---- x axis ----------------------------------------------------------------
# Genes are a discrete axis ordered by position, so the tiles stay even width
# and the sparse regions do not collapse. With 883 of them no name fits, so
# label ticks with the Mb coordinate of the gene sitting at each 20 Mb mark.
lev  <- levels(pbg$gene)
pos  <- depth$start[match(lev, as.character(depth$gene))]
label_genes <- n_gene <= MAX_GENE_LABELS
if (label_genes) {
  x_scale <- scale_x_discrete(drop = FALSE)
  x_theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7,
                                              face = ifelse(lev %in% ESCAPE_GENES, "bold", "plain"),
                                              colour = ifelse(lev %in% unique(as.character(pbg$gene[pbg$distal])),
                                                              "#B5540F", "grey25")))
} else {
  marks <- seq(0, floor(CHRX_LEN / 1e6 / 20) * 20, by = 20)
  idx   <- vapply(marks * 1e6, function(m) which.min(abs(pos - m)), integer(1))
  x_scale <- scale_x_discrete(drop = FALSE, breaks = lev[idx], labels = marks)
  x_theme <- theme(axis.text.x = element_text(size = 9))
}
# distal boundaries, drawn on the gene axis rather than in Mb
b_p <- sum(pos <= DISTAL_P_MB * 1e6) + 0.5
b_q <- sum(pos <  CHRX_LEN - DISTAL_Q_MB * 1e6) + 0.5
# Drawn as lines rather than the shaded rects the other scripts use: opaque
# tiles would hide a rect behind them and tint the allelic-ratio colours in
# front of them, and on the log-scaled depth panel an -Inf rect is a NaN.
edges <- geom_vline(xintercept = c(b_p, b_q), colour = "#B5540F", linetype = 2, linewidth = 0.4)

# ---- main heatmap ----------------------------------------------------------
# White tile borders eat a 883-gene axis alive, so they are only drawn when
# the genes are wide enough to have an edge.
tile_border <- if (n_gene <= 300) list(colour = "white", linewidth = 0.25) else list(colour = NA)
pMain <- ggplot(pbg, aes(gene, celltype, fill = ar)) +
  do.call(geom_tile, tile_border) + edges +
  # facet_grid rather than facet_wrap, and drop = TRUE, so a group only shows
  # the cell types it actually has: 9w has no stressed CM and no Epicardial
  # nuclei at all, and those two rows were drawn empty in the adult panel and
  # read as missing data rather than as an absent cell type. space = "free_y"
  # keeps the tiles the same height in every panel despite the row counts
  # differing, which scales = "free_y" alone would not.
  facet_grid(rows = vars(sample), scales = "free_y", space = "free_y",
             switch = "y", labeller = labeller(sample = SAMPLE_LAB)) +
  ar_fill_cont() + x_scale + scale_y_discrete(limits = rev, drop = TRUE) +
  labs(x = NULL, y = NULL,
       title = "chrX allelic ratio per gene and cell type: adult, Sham, TAC, aged",
       subtitle = sprintf("All %d chrX genes, no read, gene-set or cell-type filter. Blank = no SNP-overlapping read for that gene in that cell type; a cell type absent from a group is not drawn at all. Green/blue = expression from both X, red = monoallelic",
                          n_gene)) +
  theme(panel.grid = element_blank(), legend.position = "right",
        axis.text.y = element_text(size = 8), axis.ticks.x = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 13)) +
  x_theme

# ---- depth track -----------------------------------------------------------
pDepth <- ggplot(depth %>% mutate(gene = factor(as.character(gene), lev)),
                 aes(gene, total)) + edges +
  geom_col(aes(fill = total < 10 * n_obs), width = 1) +
  # the flag is a MEAN over the gene's tiles, not a minimum: an orange gene can
  # still hold one well-covered cell type, a grey one a few thin ones
  scale_fill_manual(values = c(`TRUE` = "#C97314", `FALSE` = "grey35"), name = NULL,
                    breaks = c(TRUE, FALSE),
                    labels = c(`TRUE`  = "< 10 reads per tile on average (below the usual MIN_GENE_READS cutoff)",
                               `FALSE` = ">= 10 reads per tile on average")) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  x_scale +
  labs(x = if (label_genes) "chrX genes (ordered by position, left = distal Xp)"
           else "chrX position (Mb, GRCm39), genes ordered by position - dashed: distal window boundaries",
       y = "Reads\n(all cell types,\nall groups)", caption = CAPTION) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(), legend.position = "right") +
  x_theme

save_fig(pMain / pDepth + plot_layout(heights = c(7, 1.2)),
         "AP2_chrX_four_group_all_genes_heatmap",
         min(44, max(13, 0.022 * n_gene + 8)), 14)
