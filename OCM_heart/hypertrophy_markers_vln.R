# ---------------------------------------------------------------------------
# Hypertrophic marker violins, Sham vs TAC, ventricular cardiomyocytes.
#
# Redraws the slide-17 panel in the deck's own palette and adds a per-gene
# Wilcoxon rank-sum test. Doublets are excluded via DOUBLET_FILE, the same
# convention as OCM_heart/allelic_ratio/00_functions.R.
#
#   DOUBLET_FILE=Allelic_ratio_results/scDblFinder_per_cell.txt \
#     Rscript OCM_heart/hypertrophy_markers_vln.R [OUT_DIR]
#
# Writes into OUT_DIR (default hypertrophy_markers/):
#   CM_hypertrophic_markers_violin.pdf / .png   the figure
#   CM_hypertrophic_markers_stats.txt           the test table
#
# On the statistics. The test is over CELLS, and there is one animal per
# condition, so the nuclei are pseudoreplicates: the p value measures how many
# cells were sequenced, not how reproducible the effect is between animals.
# With ~1000 cells even a trivial shift returns p < 1e-10. Read the effect
# size instead - the rank-biserial correlation is the probability that a
# random TAC cell exceeds a random Sham cell, rescaled to [-1, 1], and it does
# not inflate with n. The figure prints the effect size first for that reason.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})

OUT_DIR <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "hypertrophy_markers"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Deck palette (presentations/ipt_seminar_2026-09-15/make_figures.R SAMPLE_COL)
SAMPLE_COL <- c(Sham = "#1B9E77", TAC = "#7B3294")

heart <- readRDS("heart_seurat_object_SCT.rds")
DefaultAssay(heart) <- "RNA"
if (!"celltype" %in% colnames(heart@meta.data)) heart$celltype <- Idents(heart)

# --- doublets -------------------------------------------------------------
# scDblFinder_per_cell.txt is keyed by colnames(heart) (scDblFinder_rates.R),
# so the join is on cell name directly.
DOUBLET_FILE <- Sys.getenv("DOUBLET_FILE", "")
if (nzchar(DOUBLET_FILE)) {
  if (!file.exists(DOUBLET_FILE)) stop("DOUBLET_FILE does not exist: ", DOUBLET_FILE)
  d <- read.delim(DOUBLET_FILE, stringsAsFactors = FALSE)
  if (!all(c("cell", "class") %in% names(d)))
    stop("DOUBLET_FILE needs `cell` and `class` columns: ", DOUBLET_FILE)
  dbl <- unique(d$cell[d$class == "doublet"])
  if (!length(dbl)) stop("DOUBLET_FILE has no class == 'doublet' rows: ", DOUBLET_FILE)
  hit <- intersect(dbl, colnames(heart))
  if (!length(hit))
    stop("None of the ", length(dbl), " doublet keys match colnames(heart). ",
         "Expected keys like ", paste(head(colnames(heart), 2), collapse = ", "))
  message(sprintf("excluding %d doublets of %d nuclei (%.1f%%)",
                  length(hit), ncol(heart), 100 * length(hit) / ncol(heart)))
  heart <- heart[, setdiff(colnames(heart), hit)]
} else {
  warning("DOUBLET_FILE not set - doublets are NOT excluded", immediate. = TRUE)
}

# --- cells ----------------------------------------------------------------
# Ventricular myocytes only. The stressed-CM cluster is TAC-enriched, so
# pooling it in would fold a change in cell-type composition into what should
# be a per-cell expression comparison.
vcm <- subset(heart, subset = celltype == "Ventricular Cardiomyocytes" &
                              sample %in% c("Sham", "TAC"))
vcm$sample <- factor(vcm$sample, levels = c("Sham", "TAC"))
message("cells: ", paste(names(table(vcm$sample)), table(vcm$sample),
                         sep = " = ", collapse = ", "))

# Ctgf is Ccn2 in newer annotations - keep whichever the object carries
MARKERS <- c("Nppa", "Nppb", "Acta1", "Ankrd1", "Myh7", "Ctgf", "Ccn2",
             "Xirp2", "Myh6", "Atp2a2")
absent <- setdiff(MARKERS, rownames(vcm))
if (length(absent)) message("not in object: ", paste(absent, collapse = ", "))
MARKERS <- intersect(MARKERS, rownames(vcm))

# --- per-gene Wilcoxon ----------------------------------------------------
expr <- GetAssayData(vcm, assay = "RNA", layer = "data")
is_t <- vcm$sample == "TAC"; is_s <- vcm$sample == "Sham"
stats <- do.call(rbind, lapply(MARKERS, function(g) {
  v <- as.numeric(expr[g, ])
  w <- suppressWarnings(wilcox.test(v[is_t], v[is_s]))
  # U / (n1 n2) is P(TAC > Sham) with ties at 0.5; 2 * that - 1 is the
  # rank-biserial correlation, an effect size that does not grow with n.
  auc <- unname(w$statistic) / (sum(is_t) * sum(is_s))
  data.frame(gene = g,
             mean_Sham  = mean(v[is_s]), mean_TAC = mean(v[is_t]),
             log2FC     = log2((mean(expm1(v[is_t])) + 1) / (mean(expm1(v[is_s])) + 1)),
             pct_Sham   = mean(v[is_s] > 0), pct_TAC = mean(v[is_t] > 0),
             auc        = auc,
             rank_bis   = 2 * auc - 1,
             p_value    = w$p.value)
}))
stats$FDR <- p.adjust(stats$p_value, method = "BH")
stats <- stats[order(-abs(stats$rank_bis)), ]
write.table(stats, file.path(OUT_DIR, "CM_hypertrophic_markers_stats.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n--- ventricular cardiomyocytes, Sham vs TAC (doublet-free) ---\n")
print(stats, row.names = FALSE, digits = 3)

fmt_p <- function(p) ifelse(p < 2.2e-16, "p < 2e-16", paste0("p = ", format.pval(p, digits = 2)))
lab <- setNames(sprintf("r = %+.2f, %s", stats$rank_bis, fmt_p(stats$FDR)), stats$gene)

# --- figure ---------------------------------------------------------------
pl <- lapply(MARKERS, function(g) {
  VlnPlot(vcm, features = g, group.by = "sample", pt.size = 0,
          cols = SAMPLE_COL) +
    labs(x = NULL, y = "Expression level", title = g, subtitle = lab[[g]]) +
    # sizes are set for the slide, where this is scaled to ~5.6 in wide
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12, colour = "grey30"),
          axis.title.y = element_text(size = 11),
          legend.position = "none")
})
fig <- wrap_plots(pl, ncol = 3) +
  plot_annotation(
    title = "Hypertrophic markers in ventricular cardiomyocytes: Sham vs TAC",
    caption = paste0(
      "Wilcoxon rank-sum over cells, BH-adjusted. r is the rank-biserial correlation: ",
      "+1 means every TAC cell exceeds every Sham cell.\n",
      "One animal per condition, so cells are pseudoreplicates - the p value scales with ",
      "the number of cells sequenced, not with reproducibility. Read r, not p.\n",
      "scDblFinder doublets excluded."),
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.caption = element_text(size = 11, hjust = 0, colour = "grey30")))

h <- 3.1 * ceiling(length(MARKERS) / 3)
ggsave(file.path(OUT_DIR, "CM_hypertrophic_markers_violin.pdf"), fig, width = 10, height = h)
ggsave(file.path(OUT_DIR, "CM_hypertrophic_markers_violin.png"), fig, width = 10, height = h, dpi = 300)
message("wrote ", file.path(OUT_DIR, "CM_hypertrophic_markers_violin.png"))
