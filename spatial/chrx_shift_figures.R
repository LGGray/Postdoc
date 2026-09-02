# Figures for the chrX allelic-shift scan and the escape-clustering test.
#
# Reads what the two python scripts wrote - no statistics are recomputed here,
# so the figures cannot disagree with the tables:
#   chrx_allelic_shift_<level>.tsv        chrx_allelic_shift.py
#   chrx_shift_tiles_<level>.tsv.gz       chrx_allelic_shift.py
#   chrx_escape_clustering.tsv            chrx_escape_clustering.py
#
# PAGE 1 IS THE DIAGNOSTIC, AND IT IS FIRST ON PURPOSE. Autosomal genes cannot
# shift allelically between two sections, so their z should lie on the diagonal
# of a normal QQ plot. It does not - that departure is the reason every p-value
# on later pages is the corrected one, and anyone reading the hit list should
# see the size of the correction before the hits.
#
# THE TILE PAGES PLOT trow/tcol INDEX SPACE, not on the H&E. That is
# deliberate: this figure is about adjacency, which is exactly what the index
# grid represents, and it needs neither arrow nor the tissue positions. For the
# same genes on the H&E, use tile_gene_ar_maps.R.
#
#   IN_ROOT=... LEVEL=umi Rscript ~/Postdoc/spatial/chrx_shift_figures.R
#
##### --------------------------- CONFIG --------------------------- #####
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE    <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial"
IN_ROOT <- Sys.getenv("IN_ROOT", file.path(BASE, "ase_pysam_dup_64um"))
OUT_DIR <- Sys.getenv("OUT_DIR", IN_ROOT)
LEVEL   <- Sys.getenv("LEVEL", "umi")
CHROM   <- Sys.getenv("CHROM", "chrX")
Q_CUT   <- as.numeric(Sys.getenv("Q_CUT", "0.05"))
MAX_MAP_GENES <- as.integer(Sys.getenv("MAX_MAP_GENES", "8"))
GENES_PER_PAGE <- as.integer(Sys.getenv("GENES_PER_PAGE", "4"))
# Row index increases downward in the tile grid, as in an image.
FLIP_Y  <- !identical(tolower(Sys.getenv("FLIP_Y", "TRUE")), "false")

OUT_PDF <- Sys.getenv("OUT_PDF",
                      file.path(OUT_DIR, sprintf("chrx_shift_figures_%s.pdf", LEVEL)))

# Sample -> label. Kept as the section labels, never "adult"/"aged" as groups:
# one section per age means these are two sections, not two groups.
SLAB <- c("9w" = "9w section", "78w" = "78w section")
COL_AGE <- c("9w" = "#184f95", "78w" = "#c97314")
##### ------------------------------------------------------------- #####

msg <- function(...) message(sprintf(...))

HAVE_REPEL <- requireNamespace("ggrepel", quietly = TRUE)
if (!HAVE_REPEL) msg("note: ggrepel not installed - gene labels may overlap")

f_shift <- file.path(IN_ROOT, sprintf("chrx_allelic_shift_%s.tsv", LEVEL))
f_clust <- file.path(IN_ROOT, "chrx_escape_clustering.tsv")
f_tiles <- file.path(IN_ROOT, sprintf("chrx_shift_tiles_%s.tsv.gz", LEVEL))
for (f in c(f_shift, f_clust, f_tiles)) {
  if (!file.exists(f)) stop("Missing ", f, " - run the python steps first ",
                            "(slurm/spatial_chrx_shift.slurm)")
}

sh <- fread(f_shift, na.strings = "NA")
cl <- fread(f_clust, na.strings = c("NA", ""))
ti <- fread(f_tiles)
msg("Read %d gene rows, %d clustering rows, %d tile rows", nrow(sh), nrow(cl), nrow(ti))

# The two sample columns are named after the samples themselves, so recover them
# rather than hardcoding 9w/78w - the scan can be run on any pair.
scol <- grep("^esc_", names(sh), value = TRUE)
if (length(scol) != 2) stop("expected two esc_<sample> columns, got: ",
                            paste(scol, collapse = ", "))
S1 <- sub("^esc_", "", scol[1]); S2 <- sub("^esc_", "", scol[2])
msg("Comparison: %s -> %s (d_esc is %s minus %s)", S1, S2, S2, S1)
setnames(sh, c(paste0("esc_", S1), paste0("esc_", S2),
               paste0("n_", S1), paste0("n_", S2)),
         c("esc_1", "esc_2", "n_1", "n_2"))

auto <- sh[testable == 1 & !chrom %chin% c(CHROM, "chrY", "chrM")]
targ <- sh[testable == 1 & chrom == CHROM]
sig  <- targ[!is.na(q_corr) & q_corr < Q_CUT][order(q_corr)]
msg("%d %s genes testable, %d at q_corr < %g; %d autosomal null genes",
    nrow(targ), CHROM, nrow(sig), Q_CUT, nrow(auto))

theme_set(theme_bw(base_size = 9) +
          theme(panel.grid.minor = element_blank(),
                plot.subtitle = element_text(colour = "#555555", size = 7.5)))

pdf(OUT_PDF, width = 11, height = 8.5)

# ---- page 1: the calibration diagnostic ------------------------------------
qq <- function(z, lab) {
  z <- sort(z[is.finite(z)])
  data.table(theo = qnorm(ppoints(length(z))), obs = z, which = lab)
}
qd <- rbind(qq(auto$z_raw, "autosomal, raw"),
            qq(auto$z_corr, "autosomal, corrected"),
            qq(targ$z_raw, sprintf("%s, raw", CHROM)),
            qq(targ$z_corr, sprintf("%s, corrected", CHROM)))
# Named, so the mapping cannot be reassigned by whatever order the factor
# levels happen to sort into.
COL_QQ <- c("#8a8a8a", "#184f95", "#c97314", "#8B1913")
names(COL_QQ) <- c("autosomal, raw", "autosomal, corrected",
                   sprintf("%s, raw", CHROM), sprintf("%s, corrected", CHROM))
lam_raw <- sqrt(median(auto$z_raw^2, na.rm = TRUE) / qchisq(0.5, 1))
lam_cor <- sqrt(median(auto$z_corr^2, na.rm = TRUE) / qchisq(0.5, 1))
print(
  ggplot(qd, aes(theo, obs, colour = which)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#8a8a8a") +
    geom_point(size = 0.7, alpha = 0.7) +
    scale_colour_manual(values = COL_QQ, breaks = names(COL_QQ), name = NULL) +
    labs(title = "Why the p-values are corrected: the test is overdispersed",
         subtitle = sprintf(paste0("Autosomal genes are under no allelic regulation, so their z should follow the dashed line. ",
                                   "Inflation on autosomes: %.2f raw, %.2f corrected.\n",
                                   "The correction is a per-depth-stratum genomic control fitted on the autosomes ",
                                   "and applied to %s. Sections: %s vs %s; n = 1 animal per age."),
                            lam_raw, lam_cor, CHROM, S1, S2),
         x = "expected z (standard normal)", y = "observed z") +
    theme(legend.position = "bottom")
)

# ---- page 2: effect size against significance ------------------------------
vd <- rbind(auto[, .(gene, d_esc, q_corr, set = "autosomal (null)")],
            targ[, .(gene, d_esc, q_corr, set = CHROM)])
vd[, lq := -log10(pmax(q_corr, 1e-300))]
print(
  ggplot(vd[set == "autosomal (null)"], aes(d_esc, lq)) +
    geom_point(colour = "#cccccc", size = 0.6) +
    geom_point(data = vd[set == CHROM], colour = "#184f95", size = 1.4) +
    geom_hline(yintercept = -log10(Q_CUT), linetype = 2, colour = "#c97314") +
    # ggrepel only if it is installed; a plain label is worse but not fatal,
    # exactly as the other spatial scripts treat png and patchwork.
    (if (HAVE_REPEL)
       ggrepel::geom_text_repel(data = vd[set == CHROM & q_corr < Q_CUT],
                                aes(label = gene), size = 2.6, max.overlaps = 30)
     else geom_text(data = vd[set == CHROM & q_corr < Q_CUT], aes(label = gene),
                    size = 2.6, vjust = -0.8)) +
    labs(title = sprintf("%s genes against the autosomal null", CHROM),
         subtitle = sprintf(paste0("x > 0 = more escape in %s. Grey = autosomal genes put through the identical test; ",
                                   "they mark what this comparison\nreturns when there is nothing to find. ",
                                   "Dashed line: q_corr = %g."), S2, Q_CUT),
         x = sprintf("change in escape fraction (%s - %s)", S2, S1),
         y = expression(-log[10](q[corrected])))
)

# ---- page 3: per-gene escape, section against section ----------------------
if (nrow(sig)) {
  dl <- melt(sig[, .(gene, esc_1, esc_2, d_esc)],
             id.vars = c("gene", "d_esc"), value.name = "esc",
             variable.name = "which")
  dl[, sample := ifelse(which == "esc_1", S1, S2)]
  dl[, gene := factor(gene, levels = sig[order(d_esc), gene])]
  print(
    ggplot(dl, aes(esc, gene)) +
      geom_line(aes(group = gene), colour = "#8a8a8a", linewidth = 0.4) +
      geom_point(aes(colour = sample), size = 2.4) +
      scale_colour_manual(values = COL_AGE, breaks = c(S1, S2), name = NULL) +
      labs(title = "Escape fraction per gene, one point per section",
           subtitle = paste0("Genes at q_corr < ", Q_CUT,
                             ", ordered by the size of the change. Both directions appear: ",
                             "a uniform loss of silencing would give one.\n",
                             "With one section per age this is a difference between two sections, not an age effect."),
           x = "escape fraction (inactive-X molecules / informative molecules)", y = NULL)
  )
}

# ---- page 4: clustering, unstratified against depth-stratified -------------
cc <- cl[test %chin% c("escape_unstratified", "escape_depth_stratified")]
if (nrow(cc)) {
  cc[, z := as.numeric(z)]
  # A z of +6 on an expected count of 0.24 is arithmetic, not clustering. Those
  # rows are faded as well as starred, because a star is easy to miss and the
  # point is the same size as a well-powered one.
  cc[, flagged := underpowered == 1 | grepl("uninformative", note)]
  cc[, lab := ifelse(flagged, "*", "")]
  cc[, test := factor(test, levels = c("escape_unstratified", "escape_depth_stratified"),
                      labels = c("null: detection only", "null: detection + depth"))]
  cc[, sample := factor(sample, levels = c(S1, S2))]
  print(
    ggplot(cc[!is.na(z)], aes(z, gene, colour = test, shape = sample)) +
      geom_vline(xintercept = 0, colour = "#8a8a8a") +
      geom_vline(xintercept = 1.96, linetype = 3, colour = "#c97314") +
      geom_point(aes(alpha = ifelse(flagged, 0.3, 1)), size = 2.2,
                 position = position_dodge(width = 0.5)) +
      scale_alpha_identity(guide = "none") +
      geom_text(aes(label = lab), colour = "black", size = 4, vjust = -0.4,
                position = position_dodge(width = 0.5), show.legend = FALSE) +
      scale_colour_manual(values = c("#8a8a8a", "#184f95"), name = NULL) +
      scale_shape_manual(values = c(16, 17), name = NULL) +
      labs(title = "Is escape spatially clustered?",
           subtitle = paste0("Join counts of adjacent escape-positive tiles against a permutation null. ",
                             "Blue is the null that also holds tile depth fixed;\nthe gap between grey and blue is how much ",
                             "apparent clustering was depth. FADED AND STARRED = underpowered or expected\ncount < 2: ",
                             "no evidence either way, not evidence of dispersion. Dotted line: z = 1.96."),
           x = "z (positive = clustered)", y = NULL)
  )
}

# ---- page 5: the tiles behind it ------------------------------------------
gmap <- if (nrow(sig)) head(sig$gene, MAX_MAP_GENES) else character(0)
if (length(gmap)) {
  td <- ti[gene %chin% gmap]
  td[, status := ifelse(a2 > 0, "escape-positive", "detected, no escape")]
  td[, sample := factor(sample, levels = c(S1, S2), labels = SLAB[c(S1, S2)])]
  # A few genes per page, not all of them: coord_fixed is non-negotiable for a
  # spatial map, and eight rows of it on one page squeezes each section into a
  # strip too small to see adjacency in - which is the only thing this page is for.
  chunks <- split(gmap, ceiling(seq_along(gmap) / GENES_PER_PAGE))
  for (i in seq_along(chunks)) {
    dd <- td[gene %chin% chunks[[i]]]
    dd[, gene := factor(gene, levels = chunks[[i]])]
    pg <- ggplot(dd, aes(tcol, trow)) +
      geom_tile(data = dd[status == "detected, no escape"], fill = "#e6e5e1") +
      geom_tile(data = dd[status == "escape-positive"], fill = "#2B3186") +
      facet_grid(gene ~ sample, switch = "y") +
      coord_fixed() +
      labs(title = sprintf("Escape-positive tiles, in tile index space (%d of %d)",
                           i, length(chunks)),
           subtitle = paste0("Blue = at least one inactive-X (CAST) molecule;\n",
                             "pale = gene detected, no escape molecule.\n",
                             "Adjacency here is what page 4 tests.\n",
                             "Index space, not the H&E - for the H&E use tile_gene_ar_maps.R.\n",
                             "Most blue tiles hold ONE escape molecule: this is presence, not degree."),
           x = "tile column", y = "tile row") +
      theme(axis.text = element_text(size = 6),
            strip.text.y.left = element_text(angle = 0),
            panel.grid = element_blank())
    print(if (FLIP_Y) pg + scale_y_reverse() else pg)
  }
}

invisible(dev.off())
msg("Wrote %s", OUT_PDF)
