# ---------------------------------------------------------------------------
# 09 - Per-nucleus whole-chrX allelic ratio per cell type, RNA and ATAC as two
#      panels of one figure.
#
# The two-modality version of 07_celltype_ar_violin.R, which draws the RNA half
# alone. Reads the consolidated per-nucleus tables that 05_allelome_plots.R and
# 08_atac_allelome.R cache, so it needs no tree scan and runs in seconds.
#
# ---------------------------------------------------------------------------
# Y AXIS IS ar_b6 = A1_reads / total_reads, THE B6 FRACTION, 0-1 DIRECTIONAL
# ---------------------------------------------------------------------------
# Named explicitly because this repo carries three different quantities under
# the name `allelic_ratio` - see CLAUDE.md and 10_build_ratio_table.R:64-68.
# A1 = C57BL/6, A2 = CAST/EiJ. Xist is deleted on B6, so B6 is the ACTIVE X and
# CAST is the inactive X in every nucleus. ar_b6 near 1 is monoallelic; lower
# means signal from the inactive X.
#
# THE SAME NUMBER MEANS TWO DIFFERENT THINGS IN THE TWO PANELS, and the panels
# are stacked rather than overlaid so that is hard to miss:
#
#   RNA panel   1 - ar_b6 is escape TRANSCRIPTION from Xi   (~0.03-0.07)
#   ATAC panel  1 - ar_b6 is ACCESSIBILITY of Xi            (~0.29-0.49)
#
# The gap between the panels is the result - Xi is far more accessible than it
# is transcribed - not a discrepancy between two measurements of one thing.
# MONO_AR is therefore drawn on the RNA panel ONLY: a monoallelic bound is a
# statement about transcription and has no meaning applied to accessibility.
#
# ---------------------------------------------------------------------------
# PAIRED BY DEFAULT: THE SAME NUCLEI IN BOTH PANELS
# ---------------------------------------------------------------------------
# cellranger-arc writes the same CB tag on both BAMs, so a nucleus can be
# matched across modalities by barcode (99-100% join - see 08_atac_allelome.R).
# By default this figure keeps only nuclei that clear BOTH gates, so the two
# panels describe one population of cells and the comparison between them is
# paired rather than a contrast of two different subsets with different depth
# profiles. PAIRED=0 relaxes that to every nucleus passing its own modality's
# gate, which keeps more cells at the cost of that guarantee.
#
#   conda activate seurat_env
#   Rscript multiome/09_celltype_ar_violin_rna_atac.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

CL   <- Sys.getenv("CLUSTER_MOUNT", "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2")
WORK <- file.path(CL, "adult_aged_multiome")
RNA_DIR  <- Sys.getenv("RNA_FIG",  file.path(WORK, "figures_allelome"))
ATAC_DIR <- Sys.getenv("ATAC_FIG", file.path(WORK, "figures_atac"))
OUT      <- Sys.getenv("FIG_OUT",  file.path(WORK, "figures_allelome"))
LABELS   <- file.path(WORK, "figures_signac", "nucleus_celltype_joint.csv")
source(file.path(CL, "Postdoc", "multiome", "00_helpers.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n", sep = "")

SAMPLES   <- c("9w", "78w")
AUTOSOMES <- paste0("chr", 1:19)
# Different gates for the two modalities, because the depth differs by roughly
# 4x per nucleus: median informative chrX is ~20-25 molecules for RNA and ~75
# fragments for ATAC. A single shared gate would either throw away most RNA
# nuclei or let through ATAC ratios built on a handful of fragments.
MIN_RNA   <- as.integer(Sys.getenv("MIN_RNA",  "30"))   # as OCM 00_functions.R / 07
MIN_ATAC  <- as.integer(Sys.getenv("MIN_ATAC", "50"))   # as 08_atac_allelome.R
MONO_AR   <- as.numeric(Sys.getenv("MONO_AR",  "0.90")) # OCM_heart/allelic_ratio/00_functions.R
MIN_CELLS <- as.integer(Sys.getenv("MIN_CELLS", "20"))  # as 02_whole_chrX.R
PAIRED    <- !identical(Sys.getenv("PAIRED", "1"), "0")

MODALITY_LEVELS <- c("RNA (escape transcription)", "ATAC (Xi accessibility)")

# ---------------------------------------------------------------------------
# load - the caches, not the trees
# ---------------------------------------------------------------------------
read_cache <- function(path, what) {
  if (!file.exists(path)) {
    stop("no ", what, " per-nucleus table at ", path, "\n",
         "  run ", if (what == "RNA") "multiome_allelome_plots.slurm"
                   else "multiome_atac_plots.slurm", " first.")
  }
  d <- read_tsv(path, show_col_types = FALSE)
  need <- c("chr", "A1_reads", "A2_reads", "sample", "group")
  if (!all(need %in% names(d))) {
    stop(path, " is missing ", paste(setdiff(need, names(d)), collapse = ", "))
  }
  d
}

rna_raw  <- read_cache(file.path(RNA_DIR,  "allelome_pernucleus.tsv"),      "RNA")
atac_raw <- read_cache(file.path(ATAC_DIR, "allelome_atac_pernucleus.tsv"), "ATAC")

# A cache holding one sample must not quietly become a one-violin figure. Same
# check 07_celltype_ar_violin.R makes, for the same reason: these tables are
# written from whatever had been scored at the time, and an interrupted run or
# a half-copied tree leaves them short. 08_atac_allelome.R refuses to cache an
# incomplete scan, but the RNA cache predates that guard.
for (nm in c("RNA", "ATAC")) {
  d <- if (nm == "RNA") rna_raw else atac_raw
  miss <- setdiff(SAMPLES, as.character(unique(d$sample)))
  if (length(miss)) {
    stop("the ", nm, " per-nucleus cache has no data for: ", paste(miss, collapse = ", "),
         "\n  re-run the consolidation with REBUILD=1 to rescan the tree.")
  }
}

# ---- mapping bias, per modality, from the autosomes of the same nuclei ------
# Reads align to a B6 reference, so CAST reads align slightly worse and every
# ratio is pulled toward B6. Derivation is in 06_gene_level_escape.R.
#
# RE-ESTIMATED PER MODALITY, NOT SHARED. GEX is 101bp STAR, ATAC is 50/49bp
# BWA; a shorter read carries fewer SNPs and tolerates fewer mismatches before
# it fails to place, so CAST reads are lost at a different rate. Measured, the
# two lambdas differ by about 15%, which is larger than several of the
# per-cell-type differences this figure is drawn to show.
lambda_of <- function(d) {
  pa <- d %>% filter(chr %in% AUTOSOMES) %>%
    summarise(A1 = sum(A1_reads), A2 = sum(A2_reads))
  p_auto <- pa$A1 / (pa$A1 + pa$A2)
  if (!is.finite(p_auto) || p_auto <= 0 || p_auto >= 1) {
    stop("autosomal B6 fraction is ", p_auto, " - cannot estimate the mapping bias")
  }
  list(p_auto = p_auto, lambda = 2 - 1 / p_auto)
}
debias_with <- function(p, lambda) p * (1 - lambda) / (1 - p * lambda)

L_RNA  <- lambda_of(rna_raw)
L_ATAC <- lambda_of(atac_raw)
say("RNA  autosomal B6 fraction %.4f  ->  lambda %.4f", L_RNA$p_auto,  L_RNA$lambda)
say("ATAC autosomal B6 fraction %.4f  ->  lambda %.4f", L_ATAC$p_auto, L_ATAC$lambda)

# ---- cell type labels ------------------------------------------------------
if (!file.exists(LABELS)) stop("no cell type labels at ", LABELS, " - run 03_signac_joint.R")
meta <- read_csv(LABELS, show_col_types = FALSE) %>%
  select(sample, barcode, celltype_provisional)

chrx <- function(d, lam, gate, modality) {
  d %>% filter(chr == "chrX") %>%
    transmute(sample = as.character(sample), barcode = group,
              A1 = A1_reads, A2 = A2_reads, total_reads = A1_reads + A2_reads,
              ar_b6 = A1_reads / (A1_reads + A2_reads)) %>%
    mutate(ar_b6_corr = debias_with(ar_b6, lam),
           cast_corr  = 1 - ar_b6_corr,
           modality   = modality,
           gate       = gate)
}

rna  <- chrx(rna_raw,  L_RNA$lambda,  MIN_RNA,  MODALITY_LEVELS[1])
atac <- chrx(atac_raw, L_ATAC$lambda, MIN_ATAC, MODALITY_LEVELS[2])
say("nuclei with a chrX row: RNA %d, ATAC %d", nrow(rna), nrow(atac))

# ---------------------------------------------------------------------------
# gate, and the pairing
# ---------------------------------------------------------------------------
rna_ok  <- rna  %>% filter(total_reads >= MIN_RNA)
atac_ok <- atac %>% filter(total_reads >= MIN_ATAC)
say("passing their own gate: RNA %d (>= %d), ATAC %d (>= %d)",
    nrow(rna_ok), MIN_RNA, nrow(atac_ok), MIN_ATAC)

both_key <- inner_join(rna_ok %>% select(sample, barcode),
                       atac_ok %>% select(sample, barcode),
                       by = c("sample", "barcode"))
say("clearing BOTH gates: %d nuclei", nrow(both_key))
if (PAIRED && !nrow(both_key)) {
  stop("no nucleus clears both gates, so a paired figure cannot be drawn.\n",
       "  Both tables are keyed on cellranger-arc's CB tag, which is the\n",
       "  gex_barcode for BOTH modalities - a zero-row join points at a\n",
       "  barcode-space error upstream rather than at the gates being strict.\n",
       "  Set PAIRED=0 to draw the unpaired version regardless.")
}

if (PAIRED) {
  rna_ok  <- rna_ok  %>% semi_join(both_key, by = c("sample", "barcode"))
  atac_ok <- atac_ok %>% semi_join(both_key, by = c("sample", "barcode"))
  say("PAIRED: both panels drawn on the same %d nuclei", nrow(both_key))
} else {
  say("PAIRED=0: each panel uses every nucleus passing its own gate")
}

ar <- bind_rows(rna_ok, atac_ok) %>%
  left_join(meta, by = c("sample", "barcode")) %>%
  # as_sample() AFTER the join, and it is not redundant: meta$sample arrives
  # from read_csv as character, and a join between a factor and a character key
  # resolves to character with the levels discarded - after which everything
  # reading the column's own order sorts "78w" ahead of "9w". See SAMPLE_LEVELS
  # in 00_helpers.R.
  mutate(sample = as_sample(sample)) %>%
  filter(!is.na(celltype_provisional)) %>%
  mutate(celltype = short_labels(celltype_provisional),
         modality = factor(modality, MODALITY_LEVELS))

# The threshold is applied per cell type AND sample, on the WORSE of the two
# modalities, and removes that combination from BOTH panels.
#
# Two separate decisions, both deliberate. Taking the minimum across modalities
# means a cell type thin in ATAC cannot appear as a confident RNA violin beside
# an empty space - the figure's whole claim is a comparison between the panels,
# so a violin that has no counterpart is worse than no violin. Keeping the
# threshold per SAMPLE rather than per cell type is the 07_celltype_ar_violin.R
# convention: at n = 1 animal per age the 9w and 78w violins are separate
# claims, and one of them being thin is no reason to discard the other.
cnt <- ar %>% count(celltype, sample, modality) %>%
  group_by(celltype, sample) %>% summarise(n = min(n), .groups = "drop")
keep    <- cnt %>% filter(n >= MIN_CELLS) %>% select(celltype, sample)
dropped <- cnt %>% filter(n < MIN_CELLS)
if (nrow(dropped)) {
  say("dropped, under %d nuclei in at least one modality:", MIN_CELLS)
  print(as.data.frame(dropped), row.names = FALSE)
}
ar <- ar %>% semi_join(keep, by = c("celltype", "sample"))
if (!nrow(ar)) stop("every cell type fell below the ", MIN_CELLS, "-nucleus threshold")

# Cell types ordered by evidence, best first, so a reader scanning left to
# right meets the panels that can carry an interpretation before the ones that
# cannot. Ordered on the ATAC count because that is the sparser panel once the
# pairing is applied.
ord <- ar %>% filter(modality == MODALITY_LEVELS[2]) %>%
  count(celltype) %>% arrange(desc(n)) %>% pull(celltype)
ord <- c(ord, setdiff(unique(ar$celltype), ord))
ar <- ar %>% mutate(celltype = factor(celltype, ord))

write_tsv(ar, file.path(OUT, "celltype_ar_rna_atac.tsv"))
say("")
say("nuclei per cell type, sample and modality:")
print(as.data.frame(ar %>% count(modality, celltype, sample) %>%
  pivot_wider(names_from = sample, values_from = n, values_fill = 0)), row.names = FALSE)

# ---------------------------------------------------------------------------
# the figure
# ---------------------------------------------------------------------------
# MONO_AR on the RNA panel ONLY. A monoallelic bound is a claim about
# transcriptional output; applied to an accessibility ratio it would be a line
# with no meaning, and drawing it on both panels would invite exactly the
# misreading the header warns about.
mono_line <- data.frame(modality = factor(MODALITY_LEVELS[1], MODALITY_LEVELS),
                        y = MONO_AR)
# The autosomal expectation, which after correction is 0.5 by construction in
# both modalities - so it is the one reference line that IS common to both.
half_line <- data.frame(modality = factor(MODALITY_LEVELS, MODALITY_LEVELS), y = 0.5)

counts <- ar %>% count(modality, celltype, sample) %>%
  group_by(modality, celltype) %>%
  summarise(txt = paste(n, collapse = "/"), .groups = "drop")

# `zoom` clips the VIEW with coord_cartesian, which leaves the violin's density
# estimated on all the data and merely stops drawing outside the window. Using
# scale_y_continuous(limits=) instead would drop those nuclei BEFORE the density
# is computed and silently redraw the shape. The count outside the window is put
# in the subtitle by the caller, because a zoom that does not say what it hid is
# the figure lying by omission.
#
# NOTE scales = "free_y" was tried here first and does nothing: the RNA ratios
# run the full 0-1 range - a gated nucleus really can sit at 0 - so ggplot picks
# the same limits for both panels and the "free" page came out identical to the
# fixed one. Zooming has to be explicit.
violin <- function(yvar, ylab, ttl, sub, zoom = c(0, 1), show_mono = TRUE) {
  p <- ggplot(ar, aes(celltype, .data[[yvar]], fill = sample)) +
    geom_violin(trim = FALSE, scale = "width", bounds = c(0, 1), linewidth = 0.3,
                position = position_dodge(0.85)) +
    # group = interaction(...) is LOAD-BEARING. Setting fill = "white" outside
    # aes() removes fill from this layer's mapping, and ggplot then derives
    # `group` from the remaining discrete aesthetics - which is x alone. Both
    # samples collapse into ONE box drawn at the cell type's centre, showing
    # the pooled distribution while the violins beside it are per sample. It
    # looks like a deliberate summary rather than a bug, which is why it is
    # pinned here instead of left to the default.
    geom_boxplot(aes(group = interaction(celltype, sample)),
                 width = 0.14, outlier.size = 0.25, fill = "white", linewidth = 0.25,
                 position = position_dodge(0.85), show.legend = FALSE) +
    geom_hline(data = half_line, aes(yintercept = y),
               linetype = "dotted", colour = "grey60", linewidth = 0.3)
  if (show_mono) {
    p <- p + geom_hline(data = mono_line, aes(yintercept = y),
                        linetype = "dashed", colour = "grey35", linewidth = 0.35)
  }
  p +
    geom_text(data = counts, aes(x = celltype, y = -Inf, label = txt),
              inherit.aes = FALSE, vjust = -0.6, size = 2.3, colour = "grey30") +
    facet_wrap(~modality, nrow = 2) +
    coord_cartesian(ylim = zoom) +
    sample_scale("fill") +
    labs(x = NULL, y = ylab, fill = NULL, title = ttl, subtitle = sub) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 9.5),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          legend.position = "top")
}

SUB <- sprintf(paste0(
  "per nucleus, whole chrX, bias-corrected (RNA lambda %.3f, ATAC lambda %.3f, re-estimated per modality)\n",
  "gates: RNA >= %d informative molecules, ATAC >= %d informative fragments%s. ",
  "counts under each cell type are 9w/78w.\n",
  "n = 1 animal per age, so the 9w/78w difference is descriptive and carries no test"),
  L_RNA$lambda, L_ATAC$lambda, MIN_RNA, MIN_ATAC,
  if (PAIRED) "; the SAME nuclei in both panels" else "; each panel independent")

ZOOM_LO <- 0.5
hidden <- ar %>% group_by(modality) %>%
  summarise(n_below = sum(ar_b6_corr < ZOOM_LO), n = n(), .groups = "drop") %>%
  mutate(txt = sprintf("%s %d/%d (%.1f%%)", sub(" .*", "", modality),
                       n_below, n, 100 * n_below / n))

dev_open(file.path(OUT, "celltype_ar_violin_rna_atac.pdf"), width = 10, height = 8)

# Page 1, the full 0-1 range. This is the figure that carries the result: the
# RNA panel sits against the monoallelic bound while the ATAC panel sits far
# below it, on one scale, so the distance between the panels can be read
# directly. It is first because that distance is the finding.
print(violin("ar_b6_corr", "Allelic ratio (B6 / total)",
             "Whole-chrX allelic ratio per cell type, RNA and ATAC",
             paste0(SUB, "\nfull 0-1 range - the distance between the panels is the result: Xi is far more accessible than it is transcribed")))

# Page 2, zoomed to the top half. The RNA panel on page 1 is compressed against
# 1.0 because that is genuinely where the mass is, which makes the differences
# BETWEEN cell types unreadable. This trades the between-panel comparison back
# for that, and names what it hides.
print(violin("ar_b6_corr", "Allelic ratio (B6 / total)",
             sprintf("The same, zoomed to %.1f-1.0", ZOOM_LO),
             paste0(SUB, sprintf("\nZOOMED: densities are computed on all nuclei, but nuclei below %.1f are outside the view - %s",
                                 ZOOM_LO, paste(hidden$txt, collapse = ", "))),
             zoom = c(ZOOM_LO, 1)))

# Page 3, the same data read the other way up. 05_allelome_plots.R and
# 08_atac_allelome.R both report the CAST fraction, so this is the page that
# lines up against those tables without the reader doing 1 - x in their head.
print(violin("cast_corr", "CAST fraction (escape / Xi accessibility)",
             "The same data as the CAST fraction",
             paste0(SUB, "\n1 - allelic ratio: escape transcription on the RNA panel, Xi accessibility on the ATAC panel"),
             show_mono = FALSE))
dev.off()
say("")
say("wrote celltype_ar_violin_rna_atac.pdf (3 pages)")

# ---- the numbers behind the violins ----------------------------------------
summ <- ar %>% group_by(modality, celltype, sample) %>%
  summarise(n = n(),
            median_ar = median(ar_b6_corr),
            median_cast = median(cast_corr),
            mean_cast = mean(cast_corr),
            .groups = "drop")
write_tsv(summ, file.path(OUT, "celltype_ar_rna_atac_summary.tsv"))
say("")
say("median allelic ratio (B6/total), bias-corrected:")
print(as.data.frame(summ %>% select(modality, celltype, sample, n, median_ar) %>%
  pivot_wider(names_from = sample, values_from = c(n, median_ar), values_fill = 0)),
  row.names = FALSE, digits = 3)

say("")
say("output under %s", OUT)
