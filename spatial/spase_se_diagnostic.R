# Does scase()'s standard error know about the overdispersion it fitted?
#
# THE OBSERVATION THAT PROMPTED THIS. In the first run, SE(logit.p) divided by
# the plain binomial SE 1/sqrt(n*p*(1-p)) came out at 1.00-1.05 for every gene
# in both sections, and did not move with phi:
#
#   phi = 0.000   n =  10   median ratio 1.000
#   phi in (0,1)  n = 310   median ratio 1.012
#   phi = 1.000   n = 170   median ratio 1.010
#
# A beta-binomial with real overdispersion has a standard error LARGER than the
# binomial one, by a factor that grows with phi and with depth per unit. A
# ratio that sits at 1.0 across the whole range of phi says the fitted
# overdispersion is not reaching the reported SE. If that is right then the
# CIs in spase_scase.R are binomial CIs, the q-values down to 1e-240 are the
# exact overstatement that script was written to avoid, and NEXT_ANALYSIS task
# 7 - the autosomal ratio is 15-35x overdispersed relative to binomial - is
# unaddressed rather than addressed.
#
# It could also be right that these genes really are near-binomial at the pixel
# level and phi is unidentified, pinned on a boundary, and meaningless. That
# would be a different problem with the same appearance, so this script does
# not try to tell the two apart by reasoning. It measures.
#
# THREE INDEPENDENT CHECKS, in increasing order of authority:
#
#   1. RATIO vs PHI. The observation above, on every gene, with the boundary
#      pile-up quantified. Cheap, and it is the symptom rather than the cause.
#   2. INDEPENDENT REFIT. The same per-gene pixel counts through aod::betabin,
#      which is a different implementation of the same model. If its SE tracks
#      phi and scase()'s does not, the model is fine and the reporting is not.
#   3. PARAMETRIC BOOTSTRAP. The decisive one, because it needs no second
#      implementation to be correct. For each gene, simulate replicate pixel
#      counts from a beta-binomial at the gene's OWN fitted (p, phi) over its
#      OWN pixel depths, re-estimate p from each replicate, and take the SD of
#      logit p-hat across replicates. That is what the standard error of this
#      estimator IS, by definition. Compare it to the reported SE.
#
#      Read it like this: if the bootstrap SD matches the binomial SE, the data
#      really are binomial at this level and phi is noise - the CIs are honest
#      and only phi is meaningless. If the bootstrap SD is several times the
#      reported SE, the reported SE is wrong and every CI and q-value in the
#      scase table is too narrow by that factor.
#
# The script states which of those two the answer is, and prints the inflation
# factor to apply if it is the second.
#
# Run on the cluster:  sbatch ~/Postdoc/slurm/spatial_spase.slurm 16 "" "" sediag
#              or:     SAMPLE=9w Rscript ~/Postdoc/spatial/spase_se_diagnostic.R

SPASE_DIR <- Sys.getenv("SPASE_HOME", "")
if (!nzchar(SPASE_DIR)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  SPASE_DIR <- if (length(f)) dirname(normalizePath(f)) else getwd()
}
source(file.path(SPASE_DIR, "spase_common.R"))

# Genes to put through checks 2 and 3. Both cost real time per gene, and the
# question is about the estimator rather than about any particular gene, so a
# depth-spread sample answers it as well as the full set would.
N_GENES  <- as.integer(Sys.getenv("SE_N_GENES", "60"))
N_BOOT   <- as.integer(Sys.getenv("SE_N_BOOT", "300"))
SEED     <- as.integer(Sys.getenv("SE_SEED", "1"))
# Below this the SE question is swamped by the estimate's own noise.
MIN_UMI_SE <- as.integer(Sys.getenv("SE_MIN_UMI", "200"))

set.seed(SEED)
spase_theme()
load_spase()
dat <- load_spase_input()

SCASE_TSV <- Sys.getenv("SCASE_TSV",
  file.path(OUT_DIR, sprintf("spase_scase_%s_%dum%s.tsv",
                             SAMPLE, GENE_BIN_UM, SUF)))

##### ------------------ CHECK 1: RATIO vs PHI ----------------- #####
if (!file.exists(SCASE_TSV))
  stop("No scase table at ", SCASE_TSV, "\n  Run spase_scase.R first: ",
       "sbatch slurm/spatial_spase.slurm 16 \"\" \"\" scase", call. = FALSE)
fit <- fread(SCASE_TSV)
fit <- fit[is.finite(logit.p) & is.finite(logit.p.sd) & p > 0 & p < 1 & umi > 0]

# The binomial standard error of logit p-hat, by the delta method:
# Var(logit p) = Var(p) / (p(1-p))^2 = p(1-p)/n / (p(1-p))^2 = 1/(n p (1-p)).
fit[, sd_binom := 1 / sqrt(umi * p * (1 - p))]
fit[, ratio := logit.p.sd / sd_binom]
fit[, phi_grp := fifelse(phi < 1e-3, "phi = 0",
                  fifelse(phi > 1 - 1e-3, "phi = 1", "0 < phi < 1"))]

cat(sprintf("\n=== CHECK 1: reported SE / binomial SE, %s ===\n", SAMPLE))
print(fit[umi >= MIN_UMI_SE,
          .(genes = .N, median_ratio = round(median(ratio), 4),
            min = round(min(ratio), 3), max = round(max(ratio), 3)),
          by = phi_grp][order(phi_grp)])
cat(sprintf(paste0(
  "\nBoundary pile-up: phi is at exactly 0 for %d genes and exactly 1 for %d, ",
  "of %d\n  (%.0f%% on a boundary). A Hessian-based SE is not valid at a ",
  "boundary in any case.\n"),
  sum(fit$phi < 1e-3), sum(fit$phi > 1 - 1e-3), nrow(fit),
  100 * mean(fit$phi < 1e-3 | fit$phi > 1 - 1e-3)))

# Spearman rather than Pearson: the claim is only that the SE should RISE with
# phi, not that it should rise linearly.
rho <- suppressWarnings(cor(fit[umi >= MIN_UMI_SE]$phi,
                            fit[umi >= MIN_UMI_SE]$ratio, method = "spearman"))
cat(sprintf("Spearman(phi, SE ratio) = %+.3f. %s\n", rho,
  if (is.na(rho)) "not computable"
  else if (abs(rho) < 0.2)
    "No relationship: the fitted overdispersion is not reaching the SE."
  else "The SE does move with phi."))

##### ------------- THE GENE SAMPLE FOR CHECKS 2 AND 3 ---------- #####
# Spread over depth rather than taken from the top, so a conclusion cannot be
# an artefact of only looking at the deepest genes, and spread over chrX and
# the autosomes so it is not an artefact of only looking at chrX.
cand <- fit[umi >= MIN_UMI_SE & gene %in% rownames(dat$m_alt)]
if (!nrow(cand)) stop("No gene passes SE_MIN_UMI and is in the matrix.",
                      call. = FALSE)
pick <- function(d, k) {
  if (nrow(d) <= k) return(d$gene)
  d[order(umi)][round(seq(1, .N, length.out = k))]$gene
}
sel <- unique(c(pick(cand[is_x == TRUE], ceiling(N_GENES / 2)),
                pick(cand[is_x == FALSE], floor(N_GENES / 2))))
cat(sprintf("\nSample for checks 2-3: %d genes (%d chrX), depth %d-%d UMIs\n",
            length(sel), sum(cand[gene %in% sel]$is_x),
            min(cand[gene %in% sel]$umi), max(cand[gene %in% sel]$umi)))

# Per-gene pixel vectors, from the same matrices scase() was given, so nothing
# here can be fitting a different quantity.
alt_of <- function(g) {
  v <- dat$m_alt[g, ]; n <- v + dat$m_ref[g, ]
  k <- n > 0
  list(y = as.numeric(v[k]), n = as.numeric(n[k]))
}

##### ------------- CHECK 2: INDEPENDENT REFIT (aod) ----------- #####
have_aod <- requireNamespace("aod", quietly = TRUE)
ref <- NULL
if (!have_aod) {
  cat("\n=== CHECK 2 skipped: package 'aod' is not installed ===\n",
      "  install.packages(\"aod\") in the RNAseq env. It is one of spASE's\n",
      "  own undeclared dependencies, so it is usually already there.\n")
} else {
  cat(sprintf("\n=== CHECK 2: refit with aod::betabin, %d genes ===\n",
              length(sel)))
  ref <- rbindlist(lapply(sel, function(g) {
    d <- alt_of(g)
    if (length(d$y) < 10) return(NULL)
    df <- data.frame(y = d$y, n = d$n)
    # No `warnings` argument: aod has carried and dropped one across versions,
    # and passing it to a build that does not take it turns a fittable gene
    # into a silent NULL row.
    m <- tryCatch(suppressWarnings(
                    aod::betabin(cbind(y, n - y) ~ 1, ~ 1, data = df)),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    co <- tryCatch(summary(m)@Coef, error = function(e) NULL)
    if (is.null(co) || !nrow(co)) return(NULL)
    ph <- tryCatch(as.numeric(m@random.param[1]), error = function(e) NA_real_)
    data.table(gene = g, aod_logit_p = co[1, 1], aod_sd = co[1, 2],
               aod_phi = ph)
  }), fill = TRUE)
  if (is.null(ref) || !nrow(ref)) {
    cat("  aod::betabin returned nothing usable on any gene.\n"); ref <- NULL
  } else {
    ref <- merge(ref, fit[, .(gene, umi, p, phi, logit.p, logit.p.sd,
                              sd_binom, ratio, is_x)], by = "gene")
    ref[, sd_ratio_aod := aod_sd / sd_binom]
    ref[, sd_ratio_scase_vs_aod := logit.p.sd / aod_sd]
    cat(sprintf(paste0(
      "  aod SE / binomial SE:   median %.3f (range %.3f-%.3f)\n",
      "  scase SE / aod SE:      median %.3f\n",
      "  Spearman(phi, aod SE ratio) = %+.3f\n"),
      median(ref$sd_ratio_aod), min(ref$sd_ratio_aod), max(ref$sd_ratio_aod),
      median(ref$sd_ratio_scase_vs_aod),
      suppressWarnings(cor(ref$aod_phi, ref$sd_ratio_aod,
                           method = "spearman", use = "complete.obs"))))
    cat("  If aod's ratio rises above 1 where scase's does not, the model is\n",
        "  fine and the SE scase reports is the problem.\n")
  }
}

##### ----------- CHECK 3: PARAMETRIC BOOTSTRAP (decisive) ------ #####
# Beta-binomial draw at the gene's own fitted (p, phi) over its own pixel
# depths. spASE parameterises the overdispersion as phi in [0,1]; the
# corresponding beta shape parameters are a = p*(1/phi - 1), b = (1-p)*(1/phi - 1),
# which is the standard rho <-> (a,b) map with rho = 1/(a+b+1). phi at either
# boundary has no interior beta, so those genes are drawn binomially (phi = 0)
# or as a two-point mixture (phi = 1) - stated in the output rather than
# silently substituted.
rbetabinom <- function(n, prob, phi) {
  if (is.na(phi) || phi <= 1e-6) return(rbinom(length(n), n, prob))
  if (phi >= 1 - 1e-6) return(rbinom(length(n), n, rbinom(length(n), 1, prob)))
  s <- 1 / phi - 1
  rbinom(length(n), n, rbeta(length(n), prob * s, (1 - prob) * s))
}

cat(sprintf(paste0(
  "\n=== CHECK 3: parametric bootstrap, %d genes x %d replicates ===\n",
  "Simulating from each gene's OWN fitted (p, phi) over its OWN pixel depths,\n",
  "re-estimating p by the pooled fraction, and taking the SD of logit p-hat.\n",
  "That SD is the standard error of this estimator by definition.\n"),
  length(sel), N_BOOT))

boot <- rbindlist(lapply(sel, function(g) {
  d <- alt_of(g)
  row <- fit[gene == g]
  if (!nrow(row) || length(d$n) < 10) return(NULL)
  pg <- row$p[1]; phg <- row$phi[1]
  lp <- replicate(N_BOOT, {
    y <- rbetabinom(d$n, pg, phg)
    ph <- sum(y) / sum(d$n)
    # A replicate that lands on 0 or 1 has no logit. Nudged by half a molecule
    # rather than dropped: dropping the extreme replicates is exactly the way
    # to make a bootstrap SD come out too small.
    ph <- min(max(ph, 0.5 / sum(d$n)), 1 - 0.5 / sum(d$n))
    logit(ph)
  })
  data.table(gene = g, boot_sd = sd(lp), boot_mean = mean(lp))
}), fill = TRUE)

if (is.null(boot) || !nrow(boot)) {
  stop("The bootstrap produced no rows - nothing to conclude.", call. = FALSE)
}
boot <- merge(boot, fit[, .(gene, umi, p, phi, logit.p, logit.p.sd, sd_binom,
                            ratio, is_x, phi_grp)], by = "gene")
boot[, `:=`(boot_vs_reported = boot_sd / logit.p.sd,
            boot_vs_binom    = boot_sd / sd_binom)]

cat("\nBy phi group:\n")
print(boot[, .(genes = .N,
               boot_over_reported = round(median(boot_vs_reported), 3),
               boot_over_binomial = round(median(boot_vs_binom), 3)),
           by = phi_grp][order(phi_grp)])

med_rep <- median(boot$boot_vs_reported)
cat(sprintf("\nOverall: the bootstrap SD is %.2fx the SE scase reports.\n",
            med_rep))

##### ------------------------ THE VERDICT ---------------------- #####
verdict <- if (med_rep < 1.25) {
  "SEs_ok"
} else if (median(boot$boot_vs_binom) < 1.25) {
  "phi_meaningless"
} else {
  "SEs_too_narrow"
}
cat("\n=== VERDICT ===\n")
if (verdict == "SEs_ok") {
  cat(sprintf(paste0(
    "The reported SEs hold up: the bootstrap SD is %.2fx them, within noise.\n",
    "The flat ratio against phi in check 1 then means these genes really are\n",
    "near-binomial at the pixel level and phi is unidentified, not that the\n",
    "SE is wrong. The scase CIs and q-values can be quoted as they stand.\n"),
    med_rep))
} else if (verdict == "phi_meaningless") {
  cat(sprintf(paste0(
    "The bootstrap SD is %.2fx the reported SE but only %.2fx the BINOMIAL ",
    "SE.\nThe data are near-binomial at the pixel level and the fitted phi is ",
    "noise;\nthe CIs are roughly right, but nothing in this table is evidence ",
    "of\noverdispersion and phi should not be reported. Note this does not ",
    "contradict\nthe 15-35x overdispersion seen at the TILE level - that is a ",
    "different unit,\nand the tile result is about between-tile variance, not ",
    "within-gene.\n"),
    med_rep, median(boot$boot_vs_binom)))
} else {
  cat(sprintf(paste0(
    "The reported SEs are too narrow by a factor of about %.1f. Every CI in\n",
    "spase_scase_%s_%dum%s.tsv is that much too tight and every q-value is\n",
    "correspondingly overstated - which is what produced values down to ",
    "1e-240.\n\nWhat to do: multiply logit.p.sd by %.2f and recompute, or keep ",
    "using the\n.emp columns spase_scase.R now writes, which widen the SE by ",
    "the autosomal\nlogit MAD instead. Re-check which of the two is wider on ",
    "this data before\nquoting either.\n"),
    med_rep, SAMPLE, GENE_BIN_UM, SUF, med_rep))
}
cat("\nThis is a property of the estimator, not of the animals - it does not",
    "\ndepend on n = 1 per age and applies to both sections equally.\n")

##### ------------------------- OUTPUTS -------------------------- #####
out <- merge(boot, if (is.null(ref)) boot[, .(gene)] else
  ref[, .(gene, aod_logit_p, aod_sd, aod_phi, sd_ratio_aod,
          sd_ratio_scase_vs_aod)], by = "gene", all.x = TRUE)
out_tsv <- file.path(OUT_DIR,
  sprintf("spase_se_diagnostic_%s_%dum%s.tsv", SAMPLE, GENE_BIN_UM, SUF))
fwrite(out[order(-umi)], out_tsv, sep = "\t")
cat("\nWrote", nrow(out), "gene rows to", out_tsv, "\n")

side <- file.path(OUT_DIR,
  sprintf("spase_se_diagnostic_%s_%dum%s.provenance.tsv",
          SAMPLE, GENE_BIN_UM, SUF))
writeLines(c("key\tvalue",
  spase_provenance_common(dat),
  paste0("stage\tse_diagnostic"),
  paste0("scase_tsv\t", SCASE_TSV),
  paste0("scase_tsv_md5\t", tools::md5sum(SCASE_TSV)[[1]]),
  paste0("se_n_genes\t", length(sel)),
  paste0("se_n_boot\t", N_BOOT),
  paste0("se_seed\t", SEED),
  paste0("se_min_umi\t", MIN_UMI_SE),
  paste0("aod_available\t", have_aod),
  paste0("median_reported_over_binomial\t",
         sprintf("%.4f", median(fit[umi >= MIN_UMI_SE]$ratio))),
  paste0("spearman_phi_vs_se_ratio\t", sprintf("%.4f", rho)),
  paste0("frac_phi_on_boundary\t",
         sprintf("%.4f", mean(fit$phi < 1e-3 | fit$phi > 1 - 1e-3))),
  paste0("median_boot_over_reported\t", sprintf("%.4f", med_rep)),
  paste0("median_boot_over_binomial\t",
         sprintf("%.4f", median(boot$boot_vs_binom))),
  paste0("verdict\t", verdict)), side)
cat("Provenance written to", side, "\n")

##### -------------------------- FIGURES ------------------------- #####
pdf_path <- file.path(OUT_DIR,
  sprintf("spase_se_diagnostic_%s_%dum%s.pdf", SAMPLE, GENE_BIN_UM, SUF))
pdf(pdf_path, width = 7.5, height = 5.5)

print(ggplot(fit[umi >= MIN_UMI_SE], aes(phi, ratio)) +
  geom_hline(yintercept = 1, colour = "grey55", linetype = "dashed",
             linewidth = 0.3) +
  geom_point(colour = C_X, size = 1.1, alpha = 0.6) +
  labs(x = "fitted overdispersion phi", y = "reported SE / binomial SE",
       title = sprintf("Does the reported SE respond to phi? - %s", SAMPLE),
       subtitle = wrap(sprintf(paste0(
         "One point per gene with >=%d UMIs. A beta-binomial SE must rise ",
         "with phi. A flat cloud at 1.0 means the fitted overdispersion is ",
         "not reaching the standard error. Spearman %+.3f."),
         MIN_UMI_SE, rho))))

print(ggplot(boot, aes(logit.p.sd, boot_sd,
                       colour = ifelse(is_x, "chrX", "autosomal"))) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55",
              linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 1.3, alpha = 0.75) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(chrX = C_X, autosomal = C_A)) +
  labs(x = "SE reported by scase (log)",
       y = "bootstrap SD of logit p-hat (log)", colour = NULL,
       title = "The decisive comparison",
       subtitle = wrap(sprintf(paste0(
         "%d replicates per gene, simulated from each gene's own fitted ",
         "(p, phi) over its own pixel depths. Points on the dashed line mean ",
         "the reported SE is right. Points above it mean it is too narrow, by ",
         "the vertical distance. Median factor %.2f."),
         N_BOOT, med_rep))) +
  theme(legend.position = "top"))

if (!is.null(ref))
  print(ggplot(ref, aes(phi, sd_ratio_aod)) +
    geom_hline(yintercept = 1, colour = "grey55", linetype = "dashed",
               linewidth = 0.3) +
    geom_point(colour = C_A, size = 1.3, alpha = 0.75) +
    labs(x = "phi fitted by scase", y = "aod::betabin SE / binomial SE",
         title = "An independent implementation of the same model",
         subtitle = wrap(paste0(
           "If this rises with phi where the scase panel above is flat, the ",
           "beta-binomial is identifiable on this data and the problem is in ",
           "what scase reports, not in the model."))))

dev.off()
cat("Wrote", pdf_path, "\n")
