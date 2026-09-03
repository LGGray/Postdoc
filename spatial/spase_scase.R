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
# WHY NON-CONVERGENCE IS NOT A FITTING PROBLEM HERE. In the first run every
# single chrX gene that scase() failed on had EXACTLY ZERO CAST molecules -
# 71 of 122 at 9w, 42 of 85 at 78w, no exceptions. It is not thin depth:
# Hsd17b10 failed with 2611 UMIs and Idh3g with 9363. The likelihood is simply
# maximised on the p = 0 boundary, where logit(p) is -Inf and there is no
# Hessian to invert. Dropping those rows did three bad things, all fixed below:
#
#   1. It conditioned every downstream count on the gene having at least one
#      CAST molecule, so "20 of 51 chrX genes escape" was really 20 of 122.
#   2. It discarded the strongest silencing evidence in the dataset. Those
#      genes are not missing data, they are 25232 B6 molecules with zero CAST.
#   3. It hid the fact that the 1.5% false-escape floor cannot be right. At a
#      uniform 1.5% per-molecule error, a gene with 2611 molecules expects ~39
#      CAST and sees none. FALSE_ESCAPE_FROM_ZEROS below re-derives the floor
#      from the observed number of zero-CAST genes instead.
#
# Zero-CAST genes are now kept, with p = 0 and a one-sided Clopper-Pearson
# upper bound, and are marked fit_status = "zero_alt".
#
# WHY THE MODEL'S OWN SEs ARE NOT USED ON THEIR OWN. In the first run
# SE(logit.p) divided by the plain binomial SE was 1.00-1.05 for every gene,
# and did not move with phi: 1.000 at phi = 0, 1.010 at phi = 1. The fitted
# overdispersion was not reaching the reported standard error, so the CIs were
# binomial CIs and the q-values (down to 1e-240) were the exact overstatement
# this script was written to avoid. spase_se_diagnostic.R tests that directly.
# Until it is resolved, the tests below are ALSO run against an empirical null
# whose width is the autosomal logit MAD - gene-to-gene scatter measured on
# genes whose true value is known - and the .emp columns are the ones to quote.
#
# WHY SOME GENES ARE FLAGGED IMPOSSIBLE. In this genotype the B6 X is active in
# every cell, so escape from the CAST Xi cannot exceed 50% - that would mean
# the inactive X out-expressing the active one. Genes above it are reporting a
# broken allele call, not biology; ase_artifact_scan.R localises the SNPs
# responsible. They are flagged, excluded from the pooled figure, and the
# pooled figure is reported both ways so the exclusion is visible.
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

# The false-escape rate: a B6 molecule called CAST. The historical value comes
# from the maternally expressed imprinted control - 194 B6-expressed UMIs at
# H19/Rian/Igf2r/Meg3 with zero on the wrong allele, so <=1.5% by the rule of
# three. 194 molecules is a thin basis for it, and chrX itself turns out to
# contain far more evidence about the same quantity (see the header), so unless
# ESCAPE_FLOOR is set explicitly the floor is re-derived from the data below
# and this only supplies the fallback and the comparison.
ESCAPE_FLOOR_IMPRINTED <- as.numeric(Sys.getenv("ESCAPE_FLOOR", "0.015"))
ESCAPE_FLOOR_FIXED <- nzchar(Sys.getenv("ESCAPE_FLOOR", ""))

# Escape cannot exceed this: the CAST X is inactive in every cell, so a gene
# reading above half is reporting the inactive X out-expressing the active one.
MAX_ESCAPE <- as.numeric(Sys.getenv("MAX_ESCAPE", "0.5"))

# Genes to drop outright, comma separated - normally the ones whose SNPs
# ase_artifact_scan.R has shown to be mis-declared. Empty by default: the
# impossible-gene flag below is derived from the data and does not need a list.
EXCLUDE_GENES <- Sys.getenv("EXCLUDE_GENES", "")
EXCLUDE_GENES <- trimws(strsplit(EXCLUDE_GENES, ",")[[1]])
EXCLUDE_GENES <- EXCLUDE_GENES[nzchar(EXCLUDE_GENES)]

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

# scase() returns NA when the beta-binomial does not converge. Two very
# different situations produce that, and collapsing them is what made the first
# run's denominators wrong:
#
#   zero_alt   The gene has no CAST molecule at all, so the likelihood is
#              maximised at p = 0 where logit(p) is -Inf. There is nothing
#              wrong with the data - this IS the measurement, and on chrX in
#              this genotype it is the expected result for a silenced gene.
#              Handled below with an exact bound instead of a fitted one.
#   nofit      Everything else: the optimiser failed with CAST molecules
#              present. Genuinely no estimate.
res[, converged := is.finite(logit.p) & is.finite(logit.p.sd)]
res[, fit_status := fifelse(converged, "fitted",
                     fifelse(alt == 0, "zero_alt", "nofit"))]

# One-sided 95% upper bound for zero successes in n trials (Clopper-Pearson,
# which reduces to 1 - alpha^(1/n); the rule of three is its 3/n
# approximation). n is the molecule count: this is molecule-unit data, so the
# molecules are the independent trials and no deflation is needed.
res[, ci.high.exact := NA_real_]
res[fit_status == "zero_alt", `:=`(p = 0, ci.low = 0,
                                   ci.high.exact = 1 - 0.05^(1 / umi))]

cat(sprintf(paste0(
  "Fit status: %d fitted, %d zero_alt (no CAST molecule at all), %d nofit.\n",
  "  On %s: %d fitted, %d zero_alt, %d nofit, of %d genes.\n"),
  sum(res$fit_status == "fitted"), sum(res$fit_status == "zero_alt"),
  sum(res$fit_status == "nofit"), X_CHROM,
  sum(res$is_x & res$fit_status == "fitted"),
  sum(res$is_x & res$fit_status == "zero_alt"),
  sum(res$is_x & res$fit_status == "nofit"), sum(res$is_x)))
if (any(res$fit_status == "nofit"))
  cat("  nofit genes: ",
      paste(head(res[fit_status == "nofit"]$gene, 8), collapse = ", "), "\n")

zx <- res[is_x == TRUE & fit_status == "zero_alt"]
if (nrow(zx))
  cat(sprintf(paste0(
    "  The %d zero-CAST %s genes carry %d B6 molecules between them and zero\n",
    "  CAST. Deepest: %s. These are a silencing result, not missing data.\n"),
    nrow(zx), X_CHROM, sum(zx$ref),
    paste(sprintf("%s (%d)", zx[order(-umi)][seq_len(min(4, .N))]$gene,
                  zx[order(-umi)][seq_len(min(4, .N))]$umi), collapse = ", ")))

##### --------------- MAPPING-BIAS CALIBRATION ----------- #####
# We map F1 B6xCAST reads to a B6 reference, so CAST reads are slightly less
# mappable and every escape estimate sits BELOW its true value. The autosomal
# genes are biallelic by construction, so their fitted p should be 0.5; however
# far it sits below 0.5 is the bias, measured on the same UMI calls by the same
# model. Correction is applied on the logit scale, which is where the model's
# standard errors live.
if (length(EXCLUDE_GENES)) {
  n_before <- nrow(res)
  res <- res[!gene %in% EXCLUDE_GENES]
  cat(sprintf("Dropped %d gene(s) named in EXCLUDE_GENES.\n",
              n_before - nrow(res)))
}

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

##### -------------- GENES THAT CANNOT BE ESCAPE ---------- #####
# Escape is expression from the INACTIVE X, so it is bounded by the active X's
# output: at full escape the two alleles are equal and the fraction is 0.5.
# Above that the inactive X would have to out-express the active one. Such
# genes are an allele call that is the wrong way round; the SNPs responsible
# are what ase_artifact_scan.R goes after.
res[, impossible := is_x & converged & ci.low.adj > MAX_ESCAPE]
res[, suspect := is_x & converged & p.adj > 0.9 * MAX_ESCAPE]
imp <- res[impossible == TRUE][order(-p.adj)]
if (nrow(imp)) {
  cat(sprintf(paste0(
    "\n%d %s gene(s) are above the %.0f%% ceiling with the whole CI clear of ",
    "it -\nnot escape, a broken allele call:\n"),
    nrow(imp), X_CHROM, 100 * MAX_ESCAPE))
  print(imp[, .(gene, umi, ref, alt, escape = round(p.adj, 3),
                lo = round(ci.low.adj, 3), hi = round(ci.high.adj, 3))])
  cat("  Between them they carry", imp[, sum(alt)], "of the",
      res[is_x == TRUE, sum(alt)], "CAST molecules on", X_CHROM, "\n")
}

##### --------------- THE FALSE-ESCAPE FLOOR -------------- #####
# Re-derived from chrX itself. Under a uniform per-molecule false-escape rate
# f, a gene with n molecules shows zero CAST with probability (1-f)^n, so the
# EXPECTED number of zero-CAST genes across the observed depths is
# sum_i (1-f)^n_i. Setting that equal to the number actually observed and
# solving for f uses every chrX gene, not a set selected for being zero, so it
# is not circular.
#
# It is an UPPER bound, and that is the useful direction: genes that genuinely
# escape contribute no zeros whatever f is, so they push the observed count
# down and the solved f up. Whatever comes out, the true floor is at most that.
floor_from_zeros <- function(n, n_zero) {
  if (!length(n) || n_zero >= length(n)) return(0)
  if (n_zero <= 0) return(NA_real_)
  g <- function(f) sum((1 - f)^n) - n_zero
  if (g(1e-9) < 0) return(0)               # fewer zeros than f=0 predicts
  tryCatch(uniroot(g, c(1e-9, 0.5), tol = 1e-10)$root,
           error = function(e) NA_real_)
}
xall <- res[is_x == TRUE]
n_zero_x <- sum(xall$alt == 0)
floor_zeros <- floor_from_zeros(xall$umi, n_zero_x)
exp_zeros_imprinted <- sum((1 - ESCAPE_FLOOR_IMPRINTED)^xall$umi)
n_zero_auto <- res[is_x == FALSE, sum(alt == 0)]

cat(sprintf(paste0(
  "\nFalse-escape floor, re-derived from %s:\n",
  "  %d of %d %s genes have zero CAST molecules.\n",
  "  Uniform rate consistent with that count: %.5f (%.3f%%) - an upper bound.\n",
  "  The imprinted-control floor of %.3f%% predicts %.2f such genes, so it is\n",
  "  %s.\n",
  "  Autosomal genes with zero CAST: %d (expected ~0; they are biallelic).\n"),
  X_CHROM, n_zero_x, nrow(xall), X_CHROM,
  floor_zeros, 100 * floor_zeros, 100 * ESCAPE_FLOOR_IMPRINTED,
  exp_zeros_imprinted,
  if (is.na(floor_zeros)) "not comparable"
  else if (exp_zeros_imprinted < 0.5 * n_zero_x)
    "far too conservative and would be falsified by this data"
  else "consistent with it",
  n_zero_auto))

ESCAPE_FLOOR <- if (ESCAPE_FLOOR_FIXED || is.na(floor_zeros) || floor_zeros <= 0)
  ESCAPE_FLOOR_IMPRINTED else floor_zeros
cat(sprintf("  Testing against %.5f (%s).\n", ESCAPE_FLOOR,
            if (ESCAPE_FLOOR_FIXED) "ESCAPE_FLOOR set explicitly"
            else if (identical(ESCAPE_FLOOR, ESCAPE_FLOOR_IMPRINTED))
              "imprinted control; the re-derivation gave nothing usable"
            else "re-derived above"))

##### -------------------- THE POOLED VALUE --------------- #####
# From the COUNTS, so it includes genes with no fit: the pool is not a fit and
# does not need one. Reported three ways because the exclusions move it by more
# than a factor of two, and quoting one number without the others would hide
# that the headline figure rests on a handful of genes.
pool_of <- function(d) if (!nrow(d)) NA_real_ else
  expit(logit(d[, sum(alt) / sum(ref + alt)]) - bias_logit)
x_pool_all  <- pool_of(xall)
x_pool_poss <- pool_of(xall[impossible == FALSE | is.na(impossible)])
x_pool_clean <- pool_of(xall[(suspect == FALSE | is.na(suspect))])
x_pool_adj <- if (is.na(x_pool_poss)) x_pool_all else x_pool_poss
x_pool <- xall[, sum(alt) / sum(ref + alt)]

cat(sprintf(paste0(
  "\nPooled %s escape (bias-corrected):\n",
  "  every gene                      %.4f\n",
  "  excluding the %2d impossible     %.4f   <- the one to quote\n",
  "  excluding everything over %.0f%%   %.4f\n"),
  X_CHROM, x_pool_all, nrow(imp), x_pool_poss,
  100 * 0.9 * MAX_ESCAPE, x_pool_clean))

##### ------------------- THE TWO TESTS ------------------ #####
# Each is run twice. The model version divides by the beta-binomial's own SE.
# The .emp version widens that by the autosomal logit MAD - the scatter of
# genes whose true value is known to be 0.5 - which absorbs both the
# gene-to-gene mapping bias a single shift cannot and, while the SE problem in
# the header is unresolved, the overdispersion the model is not propagating.
# Quote the .emp columns.
res[, logit.sd.emp := sqrt(logit.sd.adj^2 + auto_mad^2)]
res[, ci.low.emp  := expit(logit.p.adj - 1.96 * logit.sd.emp)]
res[, ci.high.emp := expit(logit.p.adj + 1.96 * logit.sd.emp)]

# 1. Escape above the false-escape floor. One-sided: the only direction these
#    data can argue is that a gene escapes MORE than the error rate. This is
#    the test behind "gene X escapes".
res[, z.floor := (logit.p.adj - logit(ESCAPE_FLOOR)) / logit.sd.adj]
res[, pval.floor := pnorm(z.floor, lower.tail = FALSE)]
res[, z.floor.emp := (logit.p.adj - logit(ESCAPE_FLOOR)) / logit.sd.emp]
res[, pval.floor.emp := pnorm(z.floor.emp, lower.tail = FALSE)]

# 2. Escape different from the chromosome as a whole. Two-sided, chrX only:
#    the test behind "gene X escapes more than the average escapee".
res[, z.chrx := (logit.p.adj - logit(x_pool_adj)) / logit.sd.adj]
res[, pval.chrx := 2 * pnorm(abs(z.chrx), lower.tail = FALSE)]
res[, z.chrx.emp := (logit.p.adj - logit(x_pool_adj)) / logit.sd.emp]
res[, pval.chrx.emp := 2 * pnorm(abs(z.chrx.emp), lower.tail = FALSE)]
res[is_x == FALSE, `:=`(z.chrx = NA_real_, pval.chrx = NA_real_,
                        z.chrx.emp = NA_real_, pval.chrx.emp = NA_real_)]

# BH over the chrX genes that carry a test. The impossible ones are excluded
# from the multiple-testing family as well as from the pool - leaving them in
# would spend the correction on genes already known not to be measuring escape.
testable <- res$is_x & res$converged & !(res$impossible %in% TRUE)
for (v in c("floor", "chrx", "floor.emp", "chrx.emp"))
  res[testable, (paste0("qval.", v)) :=
        p.adjust(get(paste0("pval.", v)), method = "BH")]

# A zero-CAST gene has no z, but it does have a verdict: its exact upper bound
# either clears the floor or does not, and for the deep ones it sits well
# below. Recorded so those rows are not silently blank.
res[fit_status == "zero_alt",
    below_floor_exact := ci.high.exact < ESCAPE_FLOOR]

cat(sprintf(paste0(
  "\nAbove the floor at q<0.05, of %d %s genes in total:\n",
  "  model SE:     %d\n  empirical SE: %d   <- the one to quote\n",
  "  plus %d gene(s) with zero CAST whose exact upper bound is already below\n",
  "  the floor, i.e. measurably silenced rather than merely untested.\n"),
  nrow(xall), X_CHROM,
  sum(res$is_x & res$qval.floor < 0.05, na.rm = TRUE),
  sum(res$is_x & res$qval.floor.emp < 0.05, na.rm = TRUE),
  sum(res$below_floor_exact %in% TRUE)))

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
  paste0("escape_floor_source\t",
         if (ESCAPE_FLOOR_FIXED) "ESCAPE_FLOOR env"
         else if (identical(ESCAPE_FLOOR, ESCAPE_FLOOR_IMPRINTED)) "imprinted"
         else "rederived_from_zero_alt_genes"),
  paste0("escape_floor_imprinted\t", ESCAPE_FLOOR_IMPRINTED),
  paste0("escape_floor_from_zeros\t", sprintf("%.8f", floor_zeros)),
  paste0("max_escape\t", MAX_ESCAPE),
  paste0("exclude_genes\t", paste(EXCLUDE_GENES, collapse = ",")),
  paste0("n_genes_fitted\t", nrow(res)),
  paste0("n_genes_chrX\t", sum(res$is_x)),
  # Kept under their original keys so the two runs can be diffed, but they are
  # now a decomposition rather than one bucket: zero_alt is a measurement and
  # nofit is the only real failure.
  paste0("n_genes_nonconverged\t", sum(!res$converged)),
  paste0("n_genes_chrX_nonconverged\t", sum(!res$converged & res$is_x)),
  paste0("n_genes_chrX_zero_alt\t", sum(res$is_x & res$fit_status == "zero_alt")),
  paste0("n_genes_chrX_nofit\t", sum(res$is_x & res$fit_status == "nofit")),
  paste0("n_genes_chrX_impossible\t", sum(res$impossible %in% TRUE)),
  paste0("chrX_umi_impossible\t", if (nrow(imp)) sum(imp$alt) else 0L),
  paste0("autosomal_alt_fraction\t", sprintf("%.6f", expit(bias_logit))),
  paste0("bias_logit\t", sprintf("%.6f", bias_logit)),
  paste0("bias_logit_se\t", sprintf("%.6f", bias_se)),
  paste0("autosomal_logit_mad\t", sprintf("%.6f", auto_mad)),
  paste0("chrX_pooled_escape_raw\t", sprintf("%.6f", x_pool)),
  paste0("chrX_pooled_escape_adj\t", sprintf("%.6f", x_pool_adj)),
  paste0("chrX_pooled_escape_all_genes\t", sprintf("%.6f", x_pool_all)),
  paste0("chrX_pooled_escape_no_impossible\t", sprintf("%.6f", x_pool_poss)),
  paste0("chrX_pooled_escape_no_suspect\t", sprintf("%.6f", x_pool_clean)),
  paste0("n_chrX_above_floor_model\t",
         sum(res$is_x & res$qval.floor < 0.05, na.rm = TRUE)),
  paste0("n_chrX_above_floor_empirical\t",
         sum(res$is_x & res$qval.floor.emp < 0.05, na.rm = TRUE))), side)
cat("Provenance written to", side, "\n")

##### --------------------- FIGURES ---------------------- #####
spase_theme()

TOPN <- as.integer(Sys.getenv("TOPN", "40"))
# Impossible genes are off the panel: they are not escape and putting them on
# an escape axis invites them to be read as the strongest escapees, which is
# how they got quoted in the first place. Their absence is stated in the
# subtitle rather than left silent.
xg <- res[is_x == TRUE & fit_status %in% c("fitted", "zero_alt") &
          !(impossible %in% TRUE)][order(-umi)][seq_len(min(TOPN, .N))]
# Zero-CAST genes are on the same axis, at 0, with a bar running to their exact
# upper bound. Dropping them is what made the first run's escape distribution
# look like it started at 3%: the genes at the bottom were the ones removed.
xg[, `:=`(pt = fifelse(fit_status == "zero_alt", 0, p.adj),
          lo = fifelse(fit_status == "zero_alt", 0, ci.low.emp),
          hi = fifelse(fit_status == "zero_alt", ci.high.exact, ci.high.emp))]
xg[, gene := factor(gene, levels = gene[order(pt, hi)])]
n_zero_panel <- sum(xg$fit_status == "zero_alt")

# Panel A - the result. Two series only because a zero-CAST gene is a bound and
# a fitted gene is an estimate, and drawing them identically would claim more
# for the bounds than the data supports.
pA <- ggplot(xg, aes(pt, gene)) +
  geom_vline(xintercept = MAX_ESCAPE, colour = "grey75", linewidth = 0.3) +
  geom_vline(xintercept = x_pool_adj, colour = "grey55",
             linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = ESCAPE_FLOOR, colour = "grey55",
             linetype = "dotted", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi,
                     colour = fifelse(fit_status == "zero_alt",
                                      "zero CAST (95% upper bound)",
                                      "fitted (95% CI)")),
                 height = 0, linewidth = 0.5) +
  geom_point(aes(colour = fifelse(fit_status == "zero_alt",
                                  "zero CAST (95% upper bound)",
                                  "fitted (95% CI)")), size = 1.5) +
  scale_colour_manual(values = setNames(c(C_X, C_A),
    c("fitted (95% CI)", "zero CAST (95% upper bound)"))) +
  scale_x_continuous(limits = c(0, NA)) +
  labs(x = "escape fraction (CAST / total, bias-corrected)", y = NULL,
       colour = NULL,
       title = sprintf("Per-gene escape from the inactive X - %s", SAMPLE),
       subtitle = wrap(sprintf(paste0(
         "Top %d of %d %s genes by depth, %d of them with no CAST molecule at ",
         "all. Intervals are widened by the autosomal logit MAD (%.3f), not ",
         "the model SE alone - see the header. Dotted = false-escape floor ",
         "%.2f%%; dashed = pooled %s %.1f%%; grey = %.0f%% (full escape). ",
         "%d gene(s) above the ceiling are excluded. %dum pixels, mask %s. ",
         "n = 1 animal."),
         nrow(xg), sum(res$is_x), X_CHROM, n_zero_panel, auto_mad,
         100 * ESCAPE_FLOOR, X_CHROM, 100 * x_pool_adj, 100 * MAX_ESCAPE,
         nrow(imp), GENE_BIN_UM, SNP_LABEL))) +
  theme(axis.text.y = element_text(size = 6.5), legend.position = "top")

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
# really a depth result. Zero-CAST genes are drawn too, at 0: they are the
# deepest genes on the chromosome, so leaving them out would remove exactly the
# points that decide whether the trend exists.
depth_dt <- res[is_x == TRUE & fit_status %in% c("fitted", "zero_alt") &
                !(impossible %in% TRUE)]
depth_dt[, `:=`(pt = fifelse(fit_status == "zero_alt", 0, p.adj),
                grp = fifelse(fit_status == "zero_alt", "zero CAST", "fitted"))]
pC <- ggplot(depth_dt, aes(umi, pt, colour = grp)) +
  geom_hline(yintercept = x_pool_adj, colour = "grey55",
             linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = ESCAPE_FLOOR, colour = "grey55",
             linetype = "dotted", linewidth = 0.3) +
  geom_point(size = 1.1, alpha = 0.7) +
  scale_x_log10() +
  scale_colour_manual(values = c(fitted = C_X, `zero CAST` = C_A)) +
  labs(x = "informative UMIs per gene (log)",
       y = "escape fraction (corrected)", colour = NULL,
       title = "Escape against depth",
       subtitle = wrap(paste0(
         "A trend here would mean the per-gene escape estimate is reading ",
         "coverage, not biology. The zero-CAST genes run to the right-hand ",
         "end of the depth axis, which is what rules out their zeros being a ",
         "detection limit."))) +
  theme(legend.position = "top")

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
cat("\n", X_CHROM, " genes by depth. `escape` is bias-corrected; `hi` is the ",
    "95% upper\nlimit, from the empirical-null CI where the gene was fitted ",
    "and from the exact\nbound where it had no CAST molecule at all:\n", sep = "")
print(res[is_x == TRUE & fit_status %in% c("fitted", "zero_alt")][
      order(-umi)][seq_len(min(25, .N)),
     .(gene, umi, n_pixels, status = fit_status,
       escape = round(fifelse(fit_status == "zero_alt", 0, p.adj), 4),
       hi = round(fifelse(fit_status == "zero_alt", ci.high.exact,
                          ci.high.emp), 4),
       phi = round(phi, 3),
       q.floor.emp = signif(qval.floor.emp, 3),
       impossible = impossible %in% TRUE)])

cat("\nWHAT TO QUOTE FROM THIS RUN:\n")
cat(sprintf("  pooled %s escape        %.4f  (excluding %d impossible gene(s))\n",
            X_CHROM, x_pool_adj, nrow(imp)))
cat(sprintf("  false-escape floor      %.5f (%s)\n", ESCAPE_FLOOR,
            if (ESCAPE_FLOOR_FIXED) "set by hand" else
            if (identical(ESCAPE_FLOOR, ESCAPE_FLOOR_IMPRINTED)) "imprinted control"
            else "re-derived from the zero-CAST genes"))
cat(sprintf("  genes above it          %d of %d %s genes (empirical null)\n",
            sum(res$is_x & res$qval.floor.emp < 0.05, na.rm = TRUE),
            sum(res$is_x), X_CHROM))
cat("  use the .emp columns, not qval.floor/qval.chrx, until",
    "spase_se_diagnostic.R\n  shows the model SE responds to phi.\n")
cat("\nDo not report a difference between sections as an age effect:",
    "n = 1 animal per age.\n")
