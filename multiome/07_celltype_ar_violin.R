# ---------------------------------------------------------------------------
# 07 - Per-nucleus whole-chrX allelic ratio, violin per cell type.
#
# The multiome counterpart of the snRNA figure that
# OCM_heart/allelic_ratio/02_whole_chrX.R writes as
# whole_chr_allelic_ratio_celltype_violin_plot_facet_wrap.pdf. Same quantity,
# same gate, same monoallelic line, so the two can be read side by side.
#
# Y AXIS IS ar_b6 = A1_reads / total_reads, the B6 fraction, 0-1 DIRECTIONAL.
# Named explicitly because this repo carries three different quantities under
# the name `allelic_ratio` - see the note in CLAUDE.md and
# 10_build_ratio_table.R:64-68. B6 is the ACTIVE X here (Xist is deleted on B6,
# so CAST is the inactive X in every nucleus), which means ar_b6 near 1 is
# monoallelic and anything lower is escape. MONO_AR = 0.90 is the boundary,
# taken from OCM_heart/allelic_ratio/00_functions.R rather than restated.
#
# TWO DELIBERATE DIFFERENCES FROM THE snRNA FIGURE.
#
#  * NO SIGNIFICANCE STARS. The snRNA version annotates 9w vs 78w and Sham vs
#    TAC with FDR stars from a dispersion LRT. There is one animal per age here,
#    so that test has no biological replication - the same setting where the
#    OCM calibration put the false positive rate between 8 and 81 percent.
#    ANALYSIS_PLAN.md is explicit that the age contrast is descriptive only, so
#    the panels carry nucleus counts instead of stars.
#
#  * PANEL n IS PRINTED ON EVERY FACET. At roughly 20-25 informative chrX
#    molecules per nucleus, a violin can look perfectly smooth and still be
#    almost entirely binomial noise. The count is the reader's only defence
#    against that, so it is not optional here.
#
#   Rscript multiome/07_celltype_ar_violin.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

CL   <- Sys.getenv("CLUSTER_MOUNT", "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2")
WORK <- file.path(CL, "adult_aged_multiome")
TREE <- Sys.getenv("ALLELOME_TREE", file.path(WORK, "multiome_allelome"))
OUT  <- Sys.getenv("FIG_OUT", file.path(WORK, "figures_allelome"))
CACHE <- file.path(OUT, "allelome_pernucleus.tsv")
LABELS <- file.path(WORK, "figures_signac", "nucleus_celltype_joint.csv")
source(file.path(CL, "Postdoc", "multiome", "00_helpers.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n", sep = "")

SAMPLES   <- c("9w", "78w")
AUTOSOMES <- paste0("chr", 1:19)
MIN_TOTAL_READS <- as.integer(Sys.getenv("MIN_TOTAL_READS", "30"))  # as OCM 00_functions.R
MONO_AR   <- as.numeric(Sys.getenv("MONO_AR", "0.90"))              # ditto
MIN_CELLS <- as.integer(Sys.getenv("MIN_CELLS", "20"))              # as 02_whole_chrX.R

# ---------------------------------------------------------------------------
# load
# ---------------------------------------------------------------------------
# 05_allelome_plots.R already consolidates this tree into a TSV, so prefer the
# cache and fall back to scanning. Scanning 4547 directories is slow enough to
# be worth avoiding but cheap enough not to need a second cache of its own.
load_pernucleus <- function() {
  if (file.exists(CACHE) && !identical(Sys.getenv("REBUILD"), "1")) {
    say("using cached %s (REBUILD=1 to rescan)", basename(CACHE))
    return(read_tsv(CACHE, show_col_types = FALSE) %>% mutate(sample = as_sample(sample)))
  }
  rows <- list()
  for (id in SAMPLES) {
    paths <- list.files(file.path(TREE, "pernucleus", "out", id),
                        pattern = "^locus_table\\.txt$", recursive = TRUE, full.names = TRUE)
    say("%s: %d locus tables", id, length(paths))
    grp <- sub("_chr_annotation_mm39\\.bed_[0-9]+$", "",
               sub("^[0-9]{4}_[0-9]{2}_[0-9]{2}_", "", basename(dirname(paths))))
    for (i in seq_along(paths)) {
      d <- tryCatch(read.delim(paths[i], stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(d) || !nrow(d) || !all(c("chr","A1_reads","A2_reads") %in% names(d))) next
      d <- aggregate(cbind(A1_reads, A2_reads) ~ chr, data = d, FUN = sum)
      d$sample <- id; d$group <- grp[i]
      rows[[length(rows) + 1]] <- d
    }
  }
  if (!length(rows)) stop("no per-nucleus output under ", TREE)
  bind_rows(rows) %>% mutate(total = A1_reads + A2_reads, sample = as_sample(sample))
}

pn <- load_pernucleus()
say("%d rows, %d nuclei", nrow(pn), dplyr::n_distinct(paste(pn$sample, pn$group)))

# The whole point of this figure is the 9w/78w comparison per cell type, so a
# table holding one sample must not quietly produce a one-violin plot. The cache
# is the way that happens: 05_allelome_plots.R writes it from whatever had been
# scored at the time, and an interrupted run or a half-copied tree leaves it
# short. Checked here rather than left to the reader to notice.
missing <- setdiff(SAMPLES, as.character(unique(pn$sample)))
if (length(missing)) {
  stop("no per-nucleus data for: ", paste(missing, collapse = ", "),
       "\n  this figure compares ", paste(SAMPLES, collapse = " and "),
       " per cell type and cannot be drawn from one sample.",
       "\n  the cache at ", CACHE, " may predate the full run;",
       " re-run with REBUILD=1 to rescan the tree.")
}
say("nuclei per sample: %s",
    paste(sprintf("%s=%d", names(table(pn$sample[pn$chr == "chrX"])),
                  table(pn$sample[pn$chr == "chrX"])), collapse = ", "))

# ---- mapping bias, from the autosomes of the same nuclei ----
# Reads align to a B6 reference, so CAST reads align slightly worse and every
# ratio is pulled toward B6. Derivation of the correction is in
# 06_gene_level_escape.R; it is repeated here rather than sourced because that
# script pulls in the whole gene-level tree to run.
pa <- pn %>% filter(chr %in% AUTOSOMES) %>% summarise(A1 = sum(A1_reads), A2 = sum(A2_reads))
p_auto <- pa$A1 / (pa$A1 + pa$A2)
LAMBDA <- 2 - 1 / p_auto
debias <- function(p) p * (1 - LAMBDA) / (1 - p * LAMBDA)
say("autosomal B6 fraction %.4f  ->  CAST read loss lambda %.4f", p_auto, LAMBDA)

# ---- cell type labels ----
if (!file.exists(LABELS)) stop("no cell type labels at ", LABELS, " - run 03_signac_joint.R")
meta <- read_csv(LABELS, show_col_types = FALSE) %>%
  select(sample, barcode, celltype_provisional)

ar <- pn %>% filter(chr == "chrX") %>%
  transmute(sample = as_sample(sample), barcode = group,
            A1 = A1_reads, A2 = A2_reads, total_reads = total,
            ar_b6 = A1_reads / total, ar_b6_corr = debias(A1_reads / total)) %>%
  left_join(meta, by = c("sample", "barcode")) %>%
  filter(!is.na(celltype_provisional)) %>%
  mutate(celltype = short_labels(celltype_provisional))

say("%d nuclei with a cell type label; median informative chrX reads %.0f",
    nrow(ar), median(ar$total_reads))

# ---------------------------------------------------------------------------
# gate
# ---------------------------------------------------------------------------
# The gate is on the ACTUAL informative chrX read count, not on a UMI proxy. At
# this depth it removes most nuclei, and that is the point: below it the ratio
# is determined by sampling rather than by biology.
flt <- ar %>% filter(total_reads >= MIN_TOTAL_READS)
say("nuclei at or above the %d-read gate: %d of %d (%.1f%%)",
    MIN_TOTAL_READS, nrow(flt), nrow(ar), 100 * nrow(flt) / nrow(ar))

keep <- flt %>% count(celltype, sample) %>% filter(n >= MIN_CELLS)
dropped <- setdiff(unique(flt$celltype), unique(keep$celltype))
if (length(dropped)) say("cell types below %d nuclei in every sample, dropped: %s",
                         MIN_CELLS, paste(dropped, collapse = ", "))
flt <- flt %>% semi_join(keep, by = c("celltype", "sample"))

# Facets ordered by how much evidence stands behind them, best first, rather
# than alphabetically. A reader scanning left to right then meets the panels
# that can carry an interpretation before the ones that cannot.
ord <- flt %>% count(celltype) %>% arrange(desc(n)) %>% pull(celltype)
flt <- flt %>% mutate(celltype = factor(celltype, ord))
lab <- flt %>% count(celltype, sample) %>%
  group_by(celltype) %>%
  summarise(txt = paste(sprintf("%s n=%d", sample, n), collapse = "\n"), .groups = "drop")

write_tsv(flt, file.path(OUT, "celltype_ar_pernucleus.tsv"))
say("")
say("nuclei per cell type and sample, after the gate:")
print(as.data.frame(flt %>% count(celltype, sample) %>%
  pivot_wider(names_from = sample, values_from = n, values_fill = 0)), row.names = FALSE)

# ---------------------------------------------------------------------------
# the figure
# ---------------------------------------------------------------------------
# bounds = c(0, 1) matters: an allelic ratio cannot leave [0, 1], and without it
# a violin with trim = FALSE draws a density tail past both ends, which reads as
# nuclei that do not exist. Same call as the snRNA version.
violin <- function(yvar, ttl, sub) {
  ggplot(flt, aes(sample, .data[[yvar]], fill = sample)) +
    geom_violin(trim = FALSE, scale = "width", bounds = c(0, 1), linewidth = 0.3) +
    geom_boxplot(width = 0.12, outlier.size = 0.3, fill = "white", linewidth = 0.3) +
    geom_hline(yintercept = MONO_AR, linetype = "dashed", colour = "grey35") +
    geom_text(data = lab, aes(x = 1.5, y = 0.02, label = txt), inherit.aes = FALSE,
              size = 2.5, colour = "grey25", vjust = 0, lineheight = 0.95) +
    facet_wrap(~celltype, labeller = label_wrap_gen(width = 18)) +
    sample_scale("fill") + sample_x() +
    scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9, 1.0)) +
    coord_cartesian(ylim = c(0, 1.02)) +
    labs(x = NULL, y = "Allelic ratio (B6 / total)", fill = NULL,
         title = ttl, subtitle = sub) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 8.5),
          plot.subtitle = element_text(size = 8.5, colour = "grey30"))
}

SUB <- sprintf(paste0("per nucleus, whole chrX, >= %d informative reads; dashed = monoallelic bound %.2f\n",
                      "n = 1 animal per age, so the 9w/78w difference is descriptive and carries no test"),
               MIN_TOTAL_READS, MONO_AR)

dev_open(file.path(OUT, "celltype_ar_violin.pdf"), width = 10, height = 7.5)
print(violin("ar_b6", "Whole-chrX allelic ratio per cell type", SUB))
print(violin("ar_b6_corr", "Whole-chrX allelic ratio per cell type, mapping-bias corrected",
             paste0(SUB, sprintf("\nautosomal B6 fraction %.4f corrected to 0.5", p_auto))))
dev.off()
say("wrote celltype_ar_violin.pdf")

# ---- companion: fraction of nuclei below the monoallelic bound ----
# The interpretable magnitude behind the violins, and the same summary the
# snRNA script writes as whole_chr_fraction_escaping_barplot.pdf.
esc <- flt %>% group_by(celltype, sample) %>%
  summarise(n = n(), escaping = mean(ar_b6_corr < MONO_AR),
            median_ar = median(ar_b6_corr), .groups = "drop")
write_tsv(esc, file.path(OUT, "celltype_fraction_escaping.tsv"))

dev_open(file.path(OUT, "celltype_fraction_escaping.pdf"), width = 9, height = 6)
print(
  ggplot(esc, aes(sample, 100 * escaping, fill = sample)) +
    geom_col() +
    geom_text(aes(label = sprintf("n=%d", n)), vjust = -0.4, size = 2.6, colour = "grey25") +
    facet_wrap(~celltype, labeller = label_wrap_gen(14)) +
    sample_scale("fill") + sample_x() +
    labs(x = NULL, y = sprintf("Nuclei below the monoallelic bound (%%)  [AR < %.2f]", MONO_AR),
         fill = NULL, title = "Fraction of nuclei escaping, per cell type",
         subtitle = "bias-corrected; n = 1 animal per age, descriptive only") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), strip.background = element_blank())
)
dev.off()
say("wrote celltype_fraction_escaping.pdf")
say("")
say("output under %s", OUT)
