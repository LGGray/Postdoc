# ---------------------------------------------------------------------------
# The escaping genes only, adult / Sham / TAC / aged, on one readable page.
#
# chrX_four_group_all_genes_heatmap.R draws every gene Allelome.PRO2 reports.
# That is the right figure for "what does the chromosome look like" and the
# wrong one for "which genes escape": 883 genes over 9-11 cell types is ~50%
# blank, no gene name fits, and the escaping tiles are scattered too thinly to
# read as a pattern. This script keeps the same tiles and the same colour
# scale and selects the genes.
#
#   Rscript OCM_heart/chrX_four_group_escape_heatmap.R
#
# WHY NOT A BREADTH FILTER. The obvious tightening - keep genes measurable in
# a majority of cell types in every group - selects on coverage rather than on
# escape, and coverage is uneven between the four animals. At >= 10 reads per
# tile it keeps 139 genes but only 6 of the 11 genes in ESCAPE_GENES: Ddx3x
# and Eif2s3x escape in 25/25 and 24/25 measured tiles at a mean AR of 0.63,
# as convincingly as anything that survives, and are dropped purely because
# 78w covers 5 cell types for them instead of 6. It also throws away exactly
# the cell-type-restricted escape the question is about.
#
# So genes are selected on escape evidence instead, in two arms:
#
#   consistent      >= ESC_MIN_TILES tiles measured, escaping (AR < ESCAPE_AR)
#                   in >= ESC_MIN_FRAC of them. Escape seen across cell types.
#   cell-type       in at least one cell type, measured in >= ESC_CT_MIN_GROUPS
#   restricted      groups and escaping in every one of them. A gene that
#                   escapes in one cell type and is monoallelic in the other
#                   eight lands here rather than being filtered away.
#
# The two arms are drawn as separate x facets, so the restricted genes are
# visible as a class instead of being silently mixed into the consistent ones.
# This recovers 10 of the 11 ESCAPE_GENES; the one it does not, Sts, has 6
# measured tiles in total and arrives through the restricted arm.
#
# Config as in ap2_pseudobulk_common.R, plus:
#   ESC_MIN_READS      10. Reads per gene x cell type x group tile. Below ~10
#                      the ratio is coarse by construction, so a tile has to
#                      clear this before it counts as evidence either way.
#   ESC_MIN_TILES      12, of the 36 possible (9 cell types x 4 groups).
#   ESC_MIN_FRAC       0.5
#   ESC_CT_MIN_GROUPS  3
#
# One animal per group, so the four columns differ descriptively only.
# ---------------------------------------------------------------------------
if (Sys.getenv("GROUPS")   == "") Sys.setenv(GROUPS = "9w,Sham,TAC,78w")
if (Sys.getenv("CONTRAST") == "") Sys.setenv(CONTRAST = "four_group_escape")
if (Sys.getenv("MIN_READS")== "") Sys.setenv(MIN_READS = "1")
.this_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(sub("^--file=", "", f[1])) else "OCM_heart"
}
source(file.path(.this_dir(), "ap2_pseudobulk_common.R"))
suppressPackageStartupMessages(library(patchwork))

ESC_MIN_READS     <- as.integer(Sys.getenv("ESC_MIN_READS", "10"))
ESC_MIN_TILES     <- as.integer(Sys.getenv("ESC_MIN_TILES", "12"))
ESC_MIN_FRAC      <- as.numeric(Sys.getenv("ESC_MIN_FRAC", "0.5"))
ESC_CT_MIN_GROUPS <- as.integer(Sys.getenv("ESC_CT_MIN_GROUPS", "3"))

# DROP_CT here, unlike the all-genes figure: CM (stressed) and Epicardial are
# the two cell types that fail the read cutoff in most groups, and leaving them
# in costs two near-empty rows in every panel without adding a gene.
pbg_all <- load_ap2_chrX(groups = GROUPS)
pbg <- pbg_all %>% filter(total >= ESC_MIN_READS) %>%
  mutate(sample = factor(sample, GROUPS))
n_ct <- n_distinct(pbg$celltype)
message(n_distinct(pbg$gene), " chrX genes with at least one tile at >= ",
        ESC_MIN_READS, " reads, over ", n_ct, " cell types and ",
        length(GROUPS), " groups")

# ---- gene selection --------------------------------------------------------
per_gene <- pbg %>%
  group_by(gene) %>%
  summarise(start = first(start),
            n_tiles = n(), n_esc = sum(ar < ESCAPE_AR),
            frac_esc = n_esc / n_tiles, mean_ar = mean(ar),
            n_ct_measured = n_distinct(celltype), .groups = "drop")

consistent <- per_gene %>%
  filter(n_tiles >= ESC_MIN_TILES, frac_esc >= ESC_MIN_FRAC) %>%
  pull(gene) %>% as.character()

# One row per gene x cell type: measured in enough groups, escaping in all.
per_gene_ct <- pbg %>%
  group_by(gene, celltype) %>%
  summarise(n_grp = n(), n_esc = sum(ar < ESCAPE_AR), .groups = "drop") %>%
  filter(n_grp >= ESC_CT_MIN_GROUPS, n_esc == n_grp)
restricted <- setdiff(as.character(unique(per_gene_ct$gene)), consistent)

stopifnot(length(consistent) > 0)
message(length(consistent), " consistent + ", length(restricted),
        " cell-type-restricted = ", length(consistent) + length(restricted),
        " genes; ", length(intersect(c(consistent, restricted), ESCAPE_GENES)),
        "/", length(ESCAPE_GENES), " of ESCAPE_GENES recovered")

# Escape means biallelic, i.e. AR near 0.5. AR < ESCAPE_AR also admits genes
# sitting near ZERO, which is expression from the CAST allele only - the
# inactive X in every nucleus under this genotype - and is an artefact rather
# than escape. Not filtered here, because the boundary that defines escape
# across the repo is MONO_AR and this script must not invent a second one, but
# named on the console and in the caption so it cannot pass unnoticed.
BIALLELIC_FLOOR <- as.numeric(Sys.getenv("BIALLELIC_FLOOR", "0.25"))
cast_only <- per_gene %>%
  filter(as.character(gene) %in% c(consistent, restricted),
         mean_ar < BIALLELIC_FLOOR) %>%
  arrange(mean_ar)
if (nrow(cast_only)) {
  message("WARNING: ", nrow(cast_only), " selected gene(s) sit below AR ",
          BIALLELIC_FLOOR, ", i.e. CAST-only rather than biallelic: ",
          paste(sprintf("%s (AR %.3f over %d tiles)", cast_only$gene,
                        cast_only$mean_ar, cast_only$n_tiles), collapse = "; "),
          ". Expression from the inactive X alone is a mapping or annotation ",
          "artefact, not escape - check before using the figure.")
}

CLASS_LAB <- c(consistent = "Escaping across cell types",
               restricted = "Cell-type restricted")

# Ordered by chrX position. A horizontal gene axis reads as a coordinate
# whether or not it is one, so ordering by mean AR and marking the distal genes
# in orange invited the wrong reading of the whole figure. Position ordering
# makes the distal window a region of the axis instead of a property of the
# labels, so the orange is gone and the boundaries are drawn as lines below.
#
# The axis is still DISCRETE - one column per gene, even width - so spacing is
# rank in position, not distance in Mb. Same choice as the all-genes figure,
# for the same reason: proportional spacing collapses the dense regions.
sel <- pbg %>%
  filter(gene %in% c(consistent, restricted)) %>%
  left_join(per_gene %>% select(gene, mean_ar), by = "gene") %>%
  mutate(class = factor(if_else(as.character(gene) %in% consistent,
                                "consistent", "restricted"),
                        levels = names(CLASS_LAB)),
         gene = fct_reorder(factor(as.character(gene)), start))

write_tsv(sel %>%
            mutate(class = unname(CLASS_LAB[as.character(class)])) %>%
            select(class, sample, celltype, gene, start, distal, A1, A2, total, ar),
          file.path(OUT, "chrX_escape_genes_allelic_ratio.tsv"))
write_tsv(per_gene %>%
            filter(as.character(gene) %in% c(consistent, restricted)) %>%
            mutate(class = if_else(as.character(gene) %in% consistent,
                                   "consistent", "restricted")) %>%
            arrange(class, start),
          file.path(OUT, "chrX_escape_genes_selection.tsv"))

lev <- levels(sel$gene)

# Distal window boundaries, one pair per x facet because scales = "free_x"
# gives each facet its own 1..n discrete axis. A boundary is only drawn when it
# actually falls inside that facet's gene range - a facet holding no distal
# gene would otherwise get a line pinned against its edge.
q_start <- CHRX_LEN - DISTAL_Q_MB * 1e6
gene_pos <- sel %>% distinct(gene, start) %>% arrange(start)
bounds <- tibble::tibble(
  xint = c(sum(gene_pos$start <= DISTAL_P_MB * 1e6),
           sum(gene_pos$start <  q_start)) + 0.5
) %>%
  filter(xint > 0.5, xint < nrow(gene_pos) + 0.5)

# strwrap, because the caption is laid out against the plot width and a long
# line is silently clipped at the device edge rather than wrapped.
wrap <- function(x, width = 155) paste(strwrap(x, width = width), collapse = "\n")
CAPTION <- paste(
  wrap(sprintf("Tiles are one gene x cell type x group pseudobulk at >= %d SNP-overlapping reads; blank = below that cutoff, not necessarily silent. AR is the B6 (active X) fraction, so low = biallelic = escape, and the boundary used throughout the repo is %.2f.",
               ESC_MIN_READS, ESCAPE_AR)),
  wrap(sprintf("Genes are selected on escape rather than on coverage: %d escaping in >= %.0f%% of >= %d measured tiles, plus %d escaping in every measured group of at least one cell type (>= %d groups); which arm each gene came from is in chrX_escape_genes_selection.tsv. %d of the %d genes in ESCAPE_GENES are recovered.",
               length(consistent), 100 * ESC_MIN_FRAC, ESC_MIN_TILES,
               length(restricted), ESC_CT_MIN_GROUPS,
               length(intersect(c(consistent, restricted), ESCAPE_GENES)),
               length(ESCAPE_GENES))),
  if (nrow(cast_only))
    wrap(sprintf("CAUTION: %s below AR %.2f, i.e. expressed from CAST only. CAST is the inactive X in every nucleus here, so that is a mapping or annotation artefact, not escape.",
                 paste(sprintf("%s sits", cast_only$gene), collapse = " and "),
                 BIALLELIC_FLOOR)),
  wrap(sprintf("Genes run left to right by chrX coordinate; dashed lines mark the distal window, the first %g Mb plus the last %g Mb of the chromosome.",
               DISTAL_P_MB, DISTAL_Q_MB)),
  "One animal per group, so the four rows differ descriptively only.",
  sep = "\n")
writeLines(CAPTION, file.path(OUT, "four_group_escape_caption.txt"))

p <- ggplot(sel, aes(gene, celltype, fill = ar)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_vline(data = bounds, aes(xintercept = xint),
             inherit.aes = FALSE, colour = "#B5540F",
             linetype = 2, linewidth = 0.5) +
  # One block, not one panel per selection arm. Splitting the arms broke the
  # gene axis into two independent coordinate runs, which is exactly the
  # misreading the position ordering was meant to fix. The arm is not marked on
  # the figure at all - it is in chrX_escape_genes_selection.tsv - so the only
  # things the gene labels encode are the name and the prior-report flag.
  # free_y so a group only shows the cell types it actually has, space =
  # "free_y" to keep tile height equal across groups.
  facet_grid(rows = vars(sample), scales = "free_y", space = "free_y",
             switch = "y", labeller = labeller(sample = SAMPLE_LAB)) +
  ar_fill_cont() +
  scale_x_discrete(drop = TRUE) +
  scale_y_discrete(limits = rev, drop = TRUE) +
  labs(x = NULL, y = NULL,
       title = "chrX escape per gene and cell type: adult, Sham, TAC, aged",
       subtitle = wrap(sprintf("%d of %d measured chrX genes, ordered by chrX position; one column per gene, so spacing is order, not Mb. Bold = previously reported escaper, dashed = distal window boundary",
                               length(lev), n_distinct(pbg$gene)), 125),
       caption = CAPTION) +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text.y = element_text(size = 9),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 13),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9,
                                   face = ifelse(lev %in% ESCAPE_GENES, "bold", "plain"),
                                   colour = "black"))

save_fig(p, "AP2_chrX_four_group_escape_heatmap",
         max(9, 0.34 * length(lev) + 4), 9)
