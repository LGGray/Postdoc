# Per-gene escape from the inactive X, via spASE's beta-binomial model.
#
# WHAT THIS ADDS OVER THE TILE ANALYSIS. Everything in ase_tile_sweep.R and
# tile_ratio_map.R has already summed over gene identity: the chrX aggregate is
# built against whole-chromosome intervals, so a 12.7% escape fraction is one
# number for the whole chromosome. This script keeps gene identity and asks the
# per-gene question, with two things the tile analysis does not have:
#
#   1. A BETA-BINOMIAL null, fitted per gene. NEXT_ANALYSIS task 7 records that
#      the autosomal ratio here is 15-35x overdispersed relative to binomial,
#      which is why ~62% of tiles were being called skewed against a binomial
#      null that cannot be right. scase() fits the overdispersion (phi) per gene
#      rather than assuming it, so the CIs it reports are the ones to quote.
#   2. Genes as the unit, which is what the escape literature is written in and
#      what makes the result checkable against it. Ddx3x and Jpx near 50%,
#      Kdm6a/Kdm5c in the high 20s and Pbdc1 near 12% were already the best
#      evidence the assay measures escape; this puts CIs on them.
#
# WHY p IS THE ESCAPE FRACTION AND NEEDS NO FLIPPING. These animals carry an
# Xist deletion on the B6 X, so the CAST X is the inactive X in every cell and
# CAST expression from chrX IS escape. scase() models the probability of
# matrix1, so matrix1 is passed the ALT (CAST, Xi) counts and matrix2 the REF
# (B6, Xa) counts - see spase_common.R, which builds both. The `p` column is
# then the escape fraction directly.
#
# WHAT scase()'s OWN pval TESTS, AND WHY IT IS NOT THE ONE YOU WANT. It tests
# H0: p = 0.5, i.e. "is this gene biallelic". On chrX in this genotype that is
# trivially rejected for everything except a full escapee, so it ranks genes by
# depth more than by biology. The columns this script adds - pval.floor and
# pval.chrx - test against the false-escape floor and against the
# chromosome-wide escape fraction, which are the two nulls that mean something
# here. scase()'s pval/qval are passed through unchanged but renamed
# pval.biallelic/qval.biallelic so they cannot be quoted by accident.
#
# For the 2D spatial question - does a gene's escape vary across the section -
# see spase_spatial.R, which runs off the same matrices.
#
# Run on the cluster:  sbatch ~/Postdoc/slurm/spatial_spase.slurm
#              or:     SAMPLE=9w Rscript ~/Postdoc/spatial/spase_scase.R

# Shared config, input loading and matrix build, so this script and
# spase_spatial.R cannot disagree about what matrix1 is.
SPASE_DIR <- Sys.getenv("SPASE_HOME", "")
if (!nzchar(SPASE_DIR)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  SPASE_DIR <- if (length(f)) dirname(normalizePath(f)) else getwd()
}
source(file.path(SPASE_DIR, "spase_common.R"))

# The measured false-escape rate: a B6 molecule called CAST. From the maternally
# expressed imprinted control - 194 B6-expressed UMIs at H19/Rian/Igf2r/Meg3
# with zero on the wrong allele, so <=1.5% by the rule of three. That is the
# null a per-gene escape claim has to beat, not zero.
ESCAPE_FLOOR <- as.numeric(Sys.getenv("ESCAPE_FLOOR", "0.015"))

load_spase()
dat <- load_spase_input()

##### -------------------- FIT scase -------------------- #####
cat("Fitting beta-binomial per gene (scase)...\n")
fit <- scase(matrix1 = dat$m_alt, matrix2 = dat$m_ref,
             min.cells = MIN_PIXELS, cores = CORES, verbose = VERBOSE)
setDT(fit)
fit[, gene := as.character(gene)]
setnames(fit,
         c("pval", "qval"), c("pval.biallelic", "qval.biallelic"),
         skip_absent = TRUE)
res <- merge(fit, dat$per_gene[, .(gene, chrom, n_pixels, umi, ref, alt)],
             by = "gene", all.x = TRUE)
res[, is_x := chrom == X_CHROM]

# scase() sets `flag` and returns NA when the beta-binomial does not converge,
# which happens when a gene's overdispersion sits on the phi = 0 boundary (a
# gene that really is binomial) or when depth is too thin to identify it. Those
# genes have no estimate and are excluded from the calibration, the tests and
# the figures - reported here rather than left to surface as a dropped-rows
# warning from ggplot.
res[, converged := is.finite(logit.p) & is.finite(logit.p.sd)]
if (any(!res$converged))
  cat(sprintf(paste0("%d of %d genes did not converge (%d chrX) and carry ",
                     "no estimate: %s\n"),
              sum(!res$converged), nrow(res),
              sum(!res$converged & res$is_x),
              paste(head(res[converged == FALSE]$gene, 8), collapse = ", ")))

##### --------------- MAPPING-BIAS CALIBRATION ----------- #####
# We map F1 B6xCAST reads to a B6 reference, so CAST reads are slightly less
# mappable and every escape estimate sits BELOW its true value. The autosomal
# genes are biallelic by construction, so their fitted p should be 0.5; however
# far it sits below 0.5 is the bias, measured on the same UMI calls by the same
# model. Correction is applied on the logit scale, which is where the model's
# standard errors live.
auto <- res[is_x == FALSE & converged == TRUE]
bias_logit <- median(auto$logit.p)
# SE of a median, ~1.253 sigma/sqrt(n). This covers the uncertainty in the
# SHIFT, not gene-to-gene variation in bias, which the autosomal MAD below
# reports and which no single shift can absorb.
bias_se <- 1.253 * sd(auto$logit.p) / sqrt(nrow(auto))
auto_mad <- mad(auto$logit.p)

cat(sprintf(paste0("Autosomal calibration on %d genes: alt fraction %.4f ",
                   "(bias %+.4f), logit MAD %.3f\n"),
            nrow(auto), expit(bias_logit), expit(bias_logit) - 0.5, auto_mad))

res[, logit.p.adj  := logit.p - bias_logit]
res[, logit.sd.adj := sqrt(logit.p.sd^2 + bias_se^2)]
res[, p.adj        := expit(logit.p.adj)]
res[, ci.low.adj   := expit(logit.p.adj - 1.96 * logit.sd.adj)]
res[, ci.high.adj  := expit(logit.p.adj + 1.96 * logit.sd.adj)]

##### ------------------- THE TWO TESTS ------------------ #####
# 1. Escape above the false-escape floor. One-sided: the only direction that
#    can be argued from these data is that a gene escapes MORE than the error
#    rate. This is the test behind "gene X escapes".
res[, z.floor := (logit.p.adj - logit(ESCAPE_FLOOR)) / logit.sd.adj]
res[, pval.floor := pnorm(z.floor, lower.tail = FALSE)]
res[is_x == TRUE & converged == TRUE,
    qval.floor := p.adjust(pval.floor, method = "BH")]

# 2. Escape different from the chromosome as a whole, pooled over the chrX
#    genes that passed. Two-sided, and only meaningful on chrX. This is the
#    test behind "gene X escapes more than the average escapee".
# Pooled from the COUNTS, so it includes genes whose fit did not converge:
# the pool is not a fit and does not need one.
x_pool <- res[is_x == TRUE, sum(alt) / sum(ref + alt)]
x_pool_adj <- expit(logit(x_pool) - bias_logit)
res[, z.chrx := (logit.p.adj - logit(x_pool_adj)) / logit.sd.adj]
res[, pval.chrx := 2 * pnorm(abs(z.chrx), lower.tail = FALSE)]
res[is_x == TRUE & converged == TRUE,
    qval.chrx := p.adjust(pval.chrx, method = "BH")]
res[is_x == FALSE, `:=`(z.chrx = NA_real_, pval.chrx = NA_real_)]

cat(sprintf("Pooled chrX escape: %.4f raw, %.4f bias-corrected\n",
            x_pool, x_pool_adj))

setorder(res, -is_x, -umi)
out_tsv <- file.path(OUT_DIR,
  sprintf("spase_scase_%s_%dum%s.tsv", SAMPLE, GENE_BIN_UM, SUF))
fwrite(res, out_tsv, sep = "\t")
cat("Wrote", nrow(res), "gene rows to", out_tsv, "\n")

##### ------------------- PROVENANCE --------------------- #####
side <- file.path(OUT_DIR,
  sprintf("spase_scase_%s_%dum%s.provenance.tsv", SAMPLE, GENE_BIN_UM, SUF))
writeLines(c("key\tvalue",
  spase_provenance_common(dat),
  paste0("stage\tscase"),
  paste0("escape_floor\t", ESCAPE_FLOOR),
  paste0("n_genes_fitted\t", nrow(res)),
  paste0("n_genes_chrX\t", sum(res$is_x)),
  paste0("n_genes_nonconverged\t", sum(!res$converged)),
  paste0("n_genes_chrX_nonconverged\t", sum(!res$converged & res$is_x)),
  paste0("autosomal_alt_fraction\t", sprintf("%.6f", expit(bias_logit))),
  paste0("bias_logit\t", sprintf("%.6f", bias_logit)),
  paste0("bias_logit_se\t", sprintf("%.6f", bias_se)),
  paste0("autosomal_logit_mad\t", sprintf("%.6f", auto_mad)),
  paste0("chrX_pooled_escape_raw\t", sprintf("%.6f", x_pool)),
  paste0("chrX_pooled_escape_adj\t", sprintf("%.6f", x_pool_adj))), side)
cat("Provenance written to", side, "\n")

##### --------------------- FIGURES ---------------------- #####
spase_theme()

TOPN <- as.integer(Sys.getenv("TOPN", "40"))
xg <- res[is_x == TRUE & converged == TRUE][order(-umi)][seq_len(min(TOPN, .N))]
xg[, gene := factor(gene, levels = gene[order(p.adj)])]

# Panel A - the result. Magnitude with uncertainty, one series, so no legend:
# the title names it. Reference lines carry the two nulls that were tested.
pA <- ggplot(xg, aes(p.adj, gene)) +
  geom_vline(xintercept = 0.5, colour = "grey75", linewidth = 0.3) +
  geom_vline(xintercept = x_pool_adj, colour = "grey55",
             linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = ESCAPE_FLOOR, colour = "grey55",
             linetype = "dotted", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci.low.adj, xmax = ci.high.adj),
                 height = 0, colour = C_X, linewidth = 0.5) +
  geom_point(colour = C_X, size = 1.5) +
  scale_x_continuous(limits = c(0, NA)) +
  labs(x = "escape fraction (CAST / total, bias-corrected)", y = NULL,
       title = sprintf("Per-gene escape from the inactive X - %s", SAMPLE),
       subtitle = wrap(sprintf(paste0("beta-binomial (spASE scase), 95%% CI. ",
         "Top %d chrX genes by depth of %d fitted. Dotted = false-escape ",
         "floor %.1f%%; dashed = pooled chrX %.1f%%; grey = 50%% (full ",
         "escape). %dum pixels, mask %s. n = 1 animal."),
         nrow(xg), sum(res$is_x), 100 * ESCAPE_FLOOR, 100 * x_pool_adj,
         GENE_BIN_UM, SNP_LABEL))) +
  theme(axis.text.y = element_text(size = 6.5))

# Panel B - the calibration. Two series, so a legend, and the autosomal set
# should sit on 0.5: that is what makes the chrX values interpretable.
cal <- rbind(
  data.table(p = res[is_x == TRUE & converged == TRUE]$p, set = "chrX"),
  data.table(p = res[is_x == FALSE & converged == TRUE]$p,
             set = "autosomal control"))
cal[, set := factor(set, levels = c("chrX", "autosomal control"))]
# Each density scaled to peak at 1. The question is where the two
# distributions SIT, and the autosomal set is far more concentrated than chrX,
# so on a common density scale its spike flattens chrX into the axis.
pB <- ggplot(cal, aes(p, after_stat(scaled), fill = set, colour = set)) +
  geom_density(alpha = 0.35, linewidth = 0.5, adjust = 1.1) +
  geom_vline(xintercept = 0.5, colour = "grey55", linetype = "dashed",
             linewidth = 0.3) +
  scale_fill_manual(values = c(chrX = C_X, `autosomal control` = C_A)) +
  scale_colour_manual(values = c(chrX = C_X, `autosomal control` = C_A)) +
  labs(x = "alt (CAST) fraction, uncorrected", y = "density (scaled to 1)",
       fill = NULL,
       colour = NULL, title = "Calibration: the autosomal set is the null",
       subtitle = wrap(sprintf(paste0("Autosomal alt fraction %.4f (bias %+.4f ",
         "from 0.5), logit MAD %.3f across %d genes. The correction is a ",
         "single shift and cannot absorb gene-to-gene bias variation - the ",
         "MAD is its size."),
         expit(bias_logit), expit(bias_logit) - 0.5, auto_mad, nrow(auto)))) +
  theme(legend.position = "top")

# Panel C - the standing objection in this project is that a chrX result is
# really a depth result. One series; if escape trends with depth, it shows here.
pC <- ggplot(res[is_x == TRUE & converged == TRUE], aes(umi, p.adj)) +
  geom_hline(yintercept = x_pool_adj, colour = "grey55",
             linetype = "dashed", linewidth = 0.3) +
  geom_point(colour = C_X, size = 1.1, alpha = 0.7) +
  scale_x_log10() +
  labs(x = "informative UMIs per gene (log)",
       y = "escape fraction (corrected)",
       title = "Escape against depth",
       subtitle = wrap(paste0("A trend here would mean the per-gene escape ",
                              "estimate is reading coverage, not biology.")))

pdf_path <- file.path(OUT_DIR,
  sprintf("spase_scase_%s_%dum%s.pdf", SAMPLE, GENE_BIN_UM, SUF))
# Panel A is one row per gene, so the page has to grow with the gene count -
# at a fixed height 10 genes are lost in whitespace and 60 collide.
pdf(pdf_path, width = 7.5,
    height = max(6, min(12, 2.6 + 0.18 * nrow(xg))))
print(pA); print(pB); print(pC)
dev.off()
cat("Wrote", pdf_path, "\n")

##### ---------------- TOP OF THE TABLE ------------------ #####
cat("\nchrX genes by depth (escape, corrected, 95% CI):\n")
print(res[is_x == TRUE & converged == TRUE][order(-umi)][seq_len(min(20, .N)),
     .(gene, umi, n_pixels, escape = round(p.adj, 4),
       lo = round(ci.low.adj, 4), hi = round(ci.high.adj, 4),
       phi = round(phi, 3), q.floor = signif(qval.floor, 3))])
cat("\nDo not report a difference between sections as an age effect:",
    "n = 1 animal per age.\n")
