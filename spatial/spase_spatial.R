# Does a gene's escape from the inactive X vary across the section?
#
# spASE's spase() fits a 2D thin-plate regression spline to each gene's allelic
# ratio and likelihood-ratio-tests it against a flat baseline. spase_scase.R
# gives each gene ONE escape fraction; this asks whether that number is the
# same everywhere on the slide.
#
# WHAT THIS IS NOT LOOKING FOR. spASE was built to find clonal XCI patches in
# random-XCI F1 tissue. This genotype has an Xist deletion on the B6 X, so the
# CAST X is inactive in EVERY cell and there is no mosaic and no patches, by
# construction. A significant spatial fit here would mean the ESCAPE FRACTION
# ITSELF varies across the tissue - regional, not clonal. That is a real
# question, but it is a different one, and nothing in this output should be
# described as an XCI patch.
#
# THE PRIOR IS THAT THIS RETURNS NOTHING. The aggregate answer is already in:
# rho at 64um is 0.006 (9w) and 0.003 (78w), and C(d) is flat from 4um to 2mm
# and fully accounted for by the global escape fraction, so escape does not
# vary spatially at the chromosome level. Per gene it has never been tested,
# which is why this is worth running - but a null result is the expected one
# and is a result.
#
# THE AUTOSOMAL CONTROL IS THE WHOLE POINT. Under the null the LRT p-values
# should be uniform. Autosomal genes are biallelic everywhere and have no
# spatial allelic structure to find, so their p-value distribution measures
# this test's false-positive rate ON THIS DATA - including the overdispersion
# and the spatially smooth depth (Moran's I on a_n is 0.67-0.80) that could
# manufacture a spatial signal on their own. If the autosomal p-values are not
# uniform, no chrX hit can be believed, and the script says so rather than
# ranking the chrX hits anyway.
#
# The control set is DEPTH-MATCHED to the chrX set, one autosomal gene per chrX
# gene at the nearest total UMI count. An unmatched control has different power
# and its p-value distribution is then not the null for the chrX tests - the
# same confound NEXT_ANALYSIS task 6 raises about comparing the two sections.
#
# ONE KNOWN WEAKNESS IN THE TEST ITSELF. spASE computes the p-value as
# 1 - pchisq(abs(lrt.stat), df.diff) (spase.R). The abs() means a spline that
# fits WORSE than the baseline - a negative LRT, which optimiser trouble can
# produce - is scored as evidence FOR the spline instead of against it. That
# inflates false positives, and it is another reason the autosomal calibration
# below is not optional.
#
# Run on the cluster:  sbatch ~/Postdoc/slurm/spatial_spase.slurm 16 "" "" spatial
#              or:     SAMPLE=9w Rscript ~/Postdoc/spatial/spase_spatial.R

SPASE_DIR <- Sys.getenv("SPASE_HOME", "")
if (!nzchar(SPASE_DIR)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  SPASE_DIR <- if (length(f)) dirname(normalizePath(f)) else getwd()
}
source(file.path(SPASE_DIR, "spase_common.R"))

# Degrees of freedom for the thin-plate spline. spASE's own default is 5 and
# its docs say to try several: too few cannot express a real pattern, too many
# spend df on noise and cost power. DF=5,8 is a reasonable pair to compare.
DF <- as.integer(Sys.getenv("DF", "5"))

# spase()'s own min.umi default is 500, applied inside the function, and it
# ERRORS rather than warns if a gene in `genes` was removed by it. It is set
# explicitly here and the same filter is applied to the gene list below, so
# what gets fitted is decided in one place.
SPATIAL_MIN_UMI <- as.integer(Sys.getenv("SPATIAL_MIN_UMI", "500"))

# Cap on chrX genes fitted, deepest first. A spline per gene over every pixel
# is far more expensive than scase(); set this to bound the job, 0 for no cap.
SPATIAL_TOPN <- as.integer(Sys.getenv("SPATIAL_TOPN", "0"))

# Depth-matched autosomal controls to fit PER chrX gene. The controls are the
# calibration, and a calibration is only as good as its size: at 1 per chrX
# gene the null can never be larger than the test set, and below ~50 genes the
# KS check has no power. 3 costs 3x the control fitting time and buys a null
# that can actually detect miscalibration.
CTRL_PER_GENE <- as.integer(Sys.getenv("CTRL_PER_GENE", "3"))

load_spase()
dat <- load_spase_input()

##### ------------------ GENE SELECTION ------------------ #####
# spase() applies min.pixels and min.umi to the matrix itself and then errors
# if a requested gene did not survive, so the same filter is applied here
# first. n_pixels and umi are computed on the same quantities spase() uses:
# a gene x pixel row exists only when a molecule was counted, so n_pixels is
# rowSums((m1>0)|(m2>0)) and umi is rowSums(m1)+rowSums(m2).
elig <- dat$keep[n_pixels >= MIN_PIXELS & umi >= SPATIAL_MIN_UMI]
xg <- elig[chrom == X_CHROM][order(-umi)]
ag <- elig[chrom != X_CHROM]

if (!nrow(xg))
  stop("No chrX gene has >=", MIN_PIXELS, " pixels and >=", SPATIAL_MIN_UMI,
       " UMIs. Lower SPATIAL_MIN_UMI or use a coarser GENE_BIN_UM.",
       call. = FALSE)
if (SPATIAL_TOPN > 0 && nrow(xg) > SPATIAL_TOPN) xg <- xg[seq_len(SPATIAL_TOPN)]

# Depth-matched autosomal control: for each chrX gene, the CTRL_PER_GENE
# unused autosomal genes with the closest total UMI. Greedy from the deepest
# chrX gene down, without replacement, so the two sets have comparable power
# gene for gene rather than merely on average.
pool <- copy(ag)
matched <- character(0)
for (u in xg$umi) {
  for (k in seq_len(CTRL_PER_GENE)) {
    if (!nrow(pool)) break
    i <- which.min(abs(pool$umi - u))
    matched <- c(matched, pool$gene[i])
    pool <- pool[-i]
  }
}
ctrl <- ag[gene %in% matched]

if (!nrow(ctrl))
  stop("No eligible autosomal gene to build the control set from. Without it ",
       "the spatial test has no calibration and its p-values cannot be read.",
       call. = FALSE)
if (nrow(ctrl) < nrow(xg))
  warning("Only ", nrow(ctrl), " autosomal control genes for ", nrow(xg),
          " chrX genes - the null is thinner than the test set.", call. = FALSE)

cat(sprintf(paste0("Fitting %d chrX genes and %d depth-matched autosomal ",
                   "controls (df=%d, min.umi=%d)\n"),
            nrow(xg), nrow(ctrl), DF, SPATIAL_MIN_UMI))
cat(sprintf("  chrX UMI/gene:    median %.0f  range %.0f-%.0f\n",
            median(xg$umi), min(xg$umi), max(xg$umi)))
cat(sprintf("  control UMI/gene: median %.0f  range %.0f-%.0f\n",
            median(ctrl$umi), min(ctrl$umi), max(ctrl$umi)))

fit_genes <- c(xg$gene, ctrl$gene)

##### -------------------- FIT spase -------------------- #####
t0 <- Sys.time()
sp <- spase(matrix1 = dat$m_alt, matrix2 = dat$m_ref,
            covariates = dat$coords, genes = fit_genes,
            df = DF, min.pixels = MIN_PIXELS, min.umi = SPATIAL_MIN_UMI,
            cores = CORES, verbose = VERBOSE)
cat(sprintf("spase() took %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

res <- as.data.table(sp$result)
res[, gene := as.character(gene)]
res <- merge(res, dat$per_gene[, .(gene, chrom, n_pixels, umi, ref, alt)],
             by = "gene", all.x = TRUE)
res[, is_x := chrom == X_CHROM]
res[, converged := is.finite(chisq.p)]
if (any(!res$converged))
  cat(sprintf(paste0("%d of %d spline fits did not converge (%d chrX): %s\n"),
              sum(!res$converged), nrow(res),
              sum(!res$converged & res$is_x),
              paste(head(res[converged == FALSE]$gene, 8), collapse = ", ")))

# spase() puts a BH q-value over EVERY gene it fitted, which here includes the
# autosomal control set - those are the null, not hypotheses, and including
# them changes the chrX q-values. Recomputed within chrX only; spASE's own
# column is renamed so it cannot be quoted by accident.
setnames(res, "qval", "qval.spase.allgenes", skip_absent = TRUE)
res[is_x == TRUE & converged == TRUE,
    qval.chrx := p.adjust(chisq.p, method = "BH")]

##### --------------- THE CALIBRATION CHECK -------------- #####
ctl_p <- res[is_x == FALSE & converged == TRUE]$chisq.p
xp     <- res[is_x == TRUE  & converged == TRUE]$chisq.p
ctl_frac05 <- mean(ctl_p < 0.05)
ks <- suppressWarnings(ks.test(ctl_p, "punif"))

# A uniform null puts 5% of controls below 0.05, and well above that means the
# test is anticonservative on this data, so a chrX hit at the same threshold
# means nothing.
#
# But a KS test on a handful of controls cannot detect even gross
# miscalibration - on 10 controls it fails to reject uniform against almost
# anything - so "did not fail" is not "passed". Below MIN_CTRL_VERDICT controls
# the verdict is UNDETERMINED, not usable: that is the difference between
# having a calibration and merely having no evidence against one.
MIN_CTRL_VERDICT <- as.integer(Sys.getenv("MIN_CTRL_VERDICT", "50"))
VERDICT <- if (length(ctl_p) < MIN_CTRL_VERDICT) "UNDETERMINED" else
           if (ctl_frac05 <= 0.15 && ks$p.value > 0.01) "usable" else
           "NOT CALIBRATED"
CALIBRATED <- VERDICT == "usable"

cat("\n--- calibration on the depth-matched autosomal control ---\n")
cat(sprintf("  controls with p < 0.05: %.1f%% (%d/%d); uniform predicts 5%%\n",
            100 * ctl_frac05, sum(ctl_p < 0.05), length(ctl_p)))
cat(sprintf("  KS test against uniform: D = %.3f, p = %.3g\n",
            unname(ks$statistic), ks$p.value))
cat(sprintf("  VERDICT: %s\n", VERDICT))
if (VERDICT == "UNDETERMINED")
  cat("  Only ", length(ctl_p), " control genes (need ", MIN_CTRL_VERDICT,
      "). The KS test cannot detect\n",
      "  miscalibration at this size, so the chrX p-values below are\n",
      "  PROVISIONAL - they have no working null. Lower SPATIAL_MIN_UMI to\n",
      "  admit more autosomal genes, or use a coarser GENE_BIN_UM.\n", sep = "")
if (VERDICT == "NOT CALIBRATED")
  cat("  Autosomal genes have no spatial allelic structure to find, so an\n",
      "  excess of small p-values there is false positives. Do not report\n",
      "  the chrX hits below as spatial ASE. Try a lower DF, a coarser\n",
      "  GENE_BIN_UM, or treat the chrX-vs-control comparison as the only\n",
      "  readable output.\n", sep = "")

n_sig <- res[is_x == TRUE & converged == TRUE & qval.chrx < 0.05, .N]
cat(sprintf(paste0("\nchrX genes with a significant spatial fit (q < 0.05): ",
                   "%d of %d (controls at p<0.05 imply ~%.1f by chance)\n"),
            n_sig, length(xp), ctl_frac05 * length(xp)))
if (n_sig == 0)
  cat("  A null result, and the expected one: rho at 64um was 0.006/0.003 and\n",
      "  C(d) is flat, so escape does not vary spatially at the chromosome\n",
      "  level either. Per gene, it now also does not.\n", sep = "")

setorder(res, -is_x, chisq.p)
out_tsv <- file.path(OUT_DIR,
  sprintf("spase_spatial_%s_%dum_df%d%s.tsv", SAMPLE, GENE_BIN_UM, DF, SUF))
fwrite(res, out_tsv, sep = "\t")
cat("Wrote", nrow(res), "gene rows to", out_tsv, "\n")

##### ------------------- PROVENANCE --------------------- #####
side <- file.path(OUT_DIR,
  sprintf("spase_spatial_%s_%dum_df%d%s.provenance.tsv",
          SAMPLE, GENE_BIN_UM, DF, SUF))
writeLines(c("key\tvalue",
  spase_provenance_common(dat),
  paste0("stage\tspatial"),
  paste0("df\t", DF),
  paste0("spatial_min_umi\t", SPATIAL_MIN_UMI),
  paste0("spatial_topn\t", SPATIAL_TOPN),
  paste0("ctrl_per_gene\t", CTRL_PER_GENE),
  paste0("n_genes_chrX_fitted\t", length(xp)),
  paste0("n_genes_control_fitted\t", length(ctl_p)),
  paste0("n_genes_nonconverged\t", sum(!res$converged)),
  paste0("control_frac_p_lt_0.05\t", sprintf("%.4f", ctl_frac05)),
  paste0("control_ks_D\t", sprintf("%.4f", unname(ks$statistic))),
  paste0("control_ks_p\t", sprintf("%.4g", ks$p.value)),
  paste0("calibrated\t", CALIBRATED),
  paste0("verdict\t", VERDICT),
  paste0("min_ctrl_verdict\t", MIN_CTRL_VERDICT),
  paste0("n_chrX_significant_q05\t", n_sig)), side)
cat("Provenance written to", side, "\n")

##### --------------------- FIGURES ---------------------- #####
spase_theme()

# Panel A - THE figure. A QQ plot of the LRT p-values against the uniform
# expectation, chrX against its depth-matched control. Two series, so a legend.
# The control tracks the diagonal iff the test is calibrated; chrX above the
# diagonal is the only form a real spatial result can take here.
qq <- rbind(
  data.table(p = sort(xp),   set = "chrX"),
  data.table(p = sort(ctl_p), set = "autosomal control (depth-matched)"))
qq[, expected := (seq_len(.N) - 0.5) / .N, by = set]
qq[, set := factor(set, levels = c("chrX",
                                   "autosomal control (depth-matched)"))]
# A strong spatial fit returns chisq.p of exactly 0, so -log10 p is Inf and
# ggplot silently clips those points against the top of the panel, where they
# read as a row of near-identical results. Capped explicitly instead, and drawn
# as triangles so a capped point is never mistaken for a measured one.
qq[, obs := -log10(p)]
fin <- qq$obs[is.finite(qq$obs)]
CEIL <- if (length(fin)) max(fin) * 1.15 + 0.5 else 5
qq[, capped := !is.finite(obs) | obs > CEIL]
qq[, obs_plot := pmin(fifelse(is.finite(obs), obs, CEIL), CEIL)]
n_capped <- sum(qq$capped)
pA <- ggplot(qq, aes(-log10(expected), obs_plot, colour = set,
                     shape = capped)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55",
              linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 1.3, alpha = 0.8) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  scale_colour_manual(values = c(chrX = C_X,
    `autosomal control (depth-matched)` = C_A)) +
  labs(x = "expected -log10 p (uniform)",
       y = sprintf("observed -log10 p (capped at %.1f)", CEIL),
       colour = NULL,
       title = sprintf("Spatial ASE test, calibrated on autosomes - %s",
                       SAMPLE),
       subtitle = wrap(sprintf(paste0("Thin-plate spline vs flat baseline, ",
         "LRT (spASE spase, df=%d). %d chrX genes, %d depth-matched autosomal ",
         "controls. Controls have no spatial allelic structure to find, so ",
         "they should track the dashed line: %.1f%% sit below p=0.05 against ",
         "the 5%% a uniform null predicts (KS p=%.3g). Verdict: %s. ",
         "%d point(s) shown as triangles are capped (p below the cap). ",
         "%dum pixels, mask %s. n = 1 animal."),
         DF, length(xp), length(ctl_p), 100 * ctl_frac05, ks$p.value,
         VERDICT, n_capped, GENE_BIN_UM, SNP_LABEL))) +
  theme(legend.position = "top")

# Panel B - where the tested genes sit in depth, so a reader can see that the
# control really is matched and that any hit is not simply the deepest gene.
dep <- rbind(
  data.table(umi = res[is_x == TRUE & converged == TRUE]$umi, set = "chrX"),
  data.table(umi = res[is_x == FALSE & converged == TRUE]$umi,
             set = "autosomal control (depth-matched)"))
dep[, set := factor(set, levels = c("chrX",
                                    "autosomal control (depth-matched)"))]
pB <- ggplot(dep, aes(umi, after_stat(scaled), colour = set, fill = set)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  scale_x_log10() +
  scale_colour_manual(values = c(chrX = C_X,
    `autosomal control (depth-matched)` = C_A)) +
  scale_fill_manual(values = c(chrX = C_X,
    `autosomal control (depth-matched)` = C_A)) +
  labs(x = "informative UMIs per gene (log)", y = "density (scaled to 1)",
       colour = NULL, fill = NULL, title = "Depth of the two gene sets",
       subtitle = wrap(paste0("The control is matched gene-for-gene on depth. ",
         "If these two curves separate, the control is not the null for the ",
         "chrX tests and panel A cannot be read."))) +
  theme(legend.position = "top")

pdf_path <- file.path(OUT_DIR,
  sprintf("spase_spatial_%s_%dum_df%d%s.pdf", SAMPLE, GENE_BIN_UM, DF, SUF))
pdf(pdf_path, width = 7.5, height = 6)
print(pA); print(pB)

# Fitted surfaces, but ONLY when the test is calibrated and something is
# actually significant. Plotting a smooth surface for a gene with no spatial
# signal produces a picture that looks like structure regardless, which is
# exactly how a null result gets published as a finding.
if (CALIBRATED && n_sig > 0) {
  top <- res[is_x == TRUE & converged == TRUE & qval.chrx < 0.05][
    order(chisq.p)][seq_len(min(6, .N))]$gene
  cat("Plotting fitted surfaces for:", paste(top, collapse = ", "), "\n")
  ok <- tryCatch({
    plotSpase(matrix1 = dat$m_alt, matrix2 = dat$m_ref,
              covariates = dat$coords, spasefit = sp, genes = top)
    TRUE
  }, error = function(e) {
    message("plotSpase failed: ", conditionMessage(e))
    FALSE
  })
  if (!ok)
    cat("Surfaces skipped; the test results above are unaffected.\n")
} else {
  cat("No fitted surfaces plotted:",
      if (!CALIBRATED) sprintf("calibration verdict is %s.\n", VERDICT)
      else "no chrX gene is significant.\n")
}
dev.off()
cat("Wrote", pdf_path, "\n")

##### ---------------- TOP OF THE TABLE ------------------ #####
cat("\nchrX genes by spatial-fit p-value:\n")
print(res[is_x == TRUE & converged == TRUE][order(chisq.p)][
  seq_len(min(20, .N)),
  .(gene, umi, n_pixels, phi.baseline = round(phi.baseline, 4),
    phi.full = round(phi.full, 4), p = signif(chisq.p, 3),
    q = signif(qval.chrx, 3))])
cat("\nA significant spatial fit here is regional variation in ESCAPE.",
    "It is not\nan XCI patch: this genotype has no mosaic.",
    "n = 1 animal per age.\n")
