# ---------------------------------------------------------------------------
# Shared libraries and helper functions for the allelic ratio analysis.
# Split out of Allelic_ratio_testing.R, which is kept intact as the original.
# Sourced by 01-05; not meant to be run on its own.
# ---------------------------------------------------------------------------

library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(broom.mixed)


run_dispersion_lrt <- function(df, ref_level, comparison_level) {
  df$sample <- factor(df$sample, levels = c(ref_level, comparison_level))

  fit_const_disp <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                            family = betabinomial(link = "logit"), data = df)
  fit_var_disp <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                          dispformula = ~ sample,
                          family = betabinomial(link = "logit"), data = df)

  lrt_row <- data.frame(anova(fit_const_disp, fit_var_disp))[2, , drop = FALSE]
  disp_coef_name <- paste0("sample", comparison_level)
  disp_beta <- unname(fixef(fit_var_disp)$disp[disp_coef_name])

  direction <- if (is.na(disp_beta)) {
    NA_character_
  } else if (disp_beta < 0) {
    paste(comparison_level, "more heterogeneous")
  } else {
    paste(comparison_level, "less heterogeneous")
  }

  lrt_row$disp_beta <- disp_beta
  lrt_row$direction <- direction
  lrt_row
}

# Simulation-based false-positive rate for run_dispersion_lrt, quantifying the
# n=1-animal-per-condition caveat: if the two conditions truly share the SAME
# dispersion, and each condition is still just a single animal, how often does
# the LRT call a "significant" difference purely from animal-to-animal noise?
#
# Mechanics: simulate two animals from IDENTICAL true (mu, phi), but let each
# animal's *realized* phi wobble around that shared truth on the log scale
# (animal_disp_sd) -- this represents ordinary biological variability between
# individuals, which a single-animal-per-condition design can never separate
# from a genuine condition effect. Cells are simulated at the real per-cell
# read depths and real group sizes so the simulation matches the actual data
# structure. animal_disp_sd cannot be estimated from n=1 data -- it has to be
# assumed, so report a range of plausible values rather than trusting one.
simulate_dispersion_null_fpr <- function(df, ref_level, comparison_level,
                                         animal_disp_sd = 0.3, n_reps = 1000,
                                         seed = 1) {
  set.seed(seed)
  df$sample <- factor(df$sample, levels = c(ref_level, comparison_level))

  # anchor simulation parameters to this celltype/comparison's real data:
  # real per-group mean allelic proportion and the real shared dispersion
  # under the null (single phi fit across both groups)
  fit_null <- glmmTMB(cbind(A1_reads, A2_reads) ~ sample,
                      family = betabinomial(link = "logit"), data = df)
  mu_hat <- predict(fit_null,
                     newdata = data.frame(sample = factor(c(ref_level, comparison_level),
                                                           levels = c(ref_level, comparison_level))),
                     type = "response")
  mu_ref  <- unname(mu_hat[1])
  mu_comp <- unname(mu_hat[2])
  phi_true <- exp(unname(fixef(fit_null)$disp[["(Intercept)"]]))

  ref_reads  <- df$total_reads[df$sample == ref_level]
  comp_reads <- df$total_reads[df$sample == comparison_level]

  p_values <- rep(NA_real_, n_reps)
  n_failed <- 0

  for (i in seq_len(n_reps)) {
    # each simulated animal gets its own realized dispersion, both drawn from
    # the same distribution since the true dispersion is identical by design
    phi_ref  <- phi_true * exp(rnorm(1, 0, animal_disp_sd))
    phi_comp <- phi_true * exp(rnorm(1, 0, animal_disp_sd))

    p_ref  <- rbeta(length(ref_reads),  mu_ref  * phi_ref,  (1 - mu_ref)  * phi_ref)
    p_comp <- rbeta(length(comp_reads), mu_comp * phi_comp, (1 - mu_comp) * phi_comp)

    A1_ref  <- rbinom(length(ref_reads),  ref_reads,  p_ref)
    A1_comp <- rbinom(length(comp_reads), comp_reads, p_comp)

    df_sim <- data.frame(
      A1_reads = c(A1_ref, A1_comp),
      A2_reads = c(ref_reads - A1_ref, comp_reads - A1_comp),
      sample = c(rep(ref_level, length(ref_reads)), rep(comparison_level, length(comp_reads)))
    )

    lrt_row <- tryCatch(run_dispersion_lrt(df_sim, ref_level, comparison_level),
                         error = function(e) NULL, warning = function(w) NULL)

    if (is.null(lrt_row)) {
      n_failed <- n_failed + 1
    } else {
      p_values[i] <- lrt_row$Pr..Chisq.
    }
  }

  list(
    fpr = mean(p_values < 0.05, na.rm = TRUE),
    n_failed = n_failed,
    n_reps = n_reps,
    p_values = p_values
  )
}

# Logistic regression LRT testing whether the proportion of monoallelic-like
# cells (allelic_ratio <= 0.1 or >= 0.9) differs between two sample levels.
run_monoallelic_lrt <- function(df, ref_level, comparison_level) {
  df$sample <- factor(df$sample, levels = c(ref_level, comparison_level))

  fit_null <- glm(monoallelic ~ 1, family = binomial, data = df)
  fit_full <- glm(monoallelic ~ sample, family = binomial, data = df)

  lrt <- anova(fit_null, fit_full, test = "LRT")

  coef_name <- paste0("sample", comparison_level)
  beta <- unname(coef(fit_full)[coef_name])

  direction <- if (is.na(beta)) {
    NA_character_
  } else if (beta > 0) {
    paste(comparison_level, "more monoallelic")
  } else {
    paste(comparison_level, "less monoallelic")
  }

  data.frame(
    Df = lrt$Df[2],
    Deviance = lrt$Deviance[2],
    p_value = lrt[["Pr(>Chi)"]][2],
    beta = beta,
    direction = direction
  )
}

fdr_to_stars <- function(fdr) {
  ifelse(is.na(fdr), "ns",
         ifelse(fdr <= 0.001, "***",
                ifelse(fdr <= 0.01, "**",
        ifelse(fdr <= 0.05, "*", "ns"))))
}

# ---------------------------------------------------------------------------
# Total-read cutoff and the output directory keyed to it.
#
# The cutoff is not a neutral QC knob: 06_cutoff_sweep.R shows mean AR rising
# from 0.85 at 10-25 chrX reads to 0.97 above 200, a real depth bias rather
# than sampling noise, so results are only comparable within one cutoff.
# Hence one directory per cutoff -- re-running at a different value writes
# somewhere new instead of silently overwriting the previous set.
#
# Defined here rather than in 02 so that 03 and 05, which read 02's per-cell
# metadata handoff, cannot drift onto a different cutoff than the file they
# are reading was built with.
#
# Override for a one-off:  MIN_TOTAL_READS=50 Rscript allelic_ratio/02_whole_chrX.R
MIN_TOTAL_READS <- as.integer(Sys.getenv("MIN_TOTAL_READS", "30"))
stopifnot(!is.na(MIN_TOTAL_READS), MIN_TOTAL_READS >= 1)

CUTOFF_DIR <- file.path("Allelic_ratio_results",
                        paste0("cutoff_", MIN_TOTAL_READS))
dir.create(CUTOFF_DIR, showWarnings = FALSE, recursive = TRUE)

# The raw per-cell allelic ratio table is cutoff-INDEPENDENT: it is every cell
# with any chrX coverage, and both 02 and 06 filter it themselves. It stays at
# the top level rather than being duplicated into each cutoff directory.
ALLELIC_RATIOS_FILE <- "Allelic_ratio_results/whole_chr_allelic_ratios.txt"

# The core escape block (04) needs its OWN cutoff, not MIN_TOTAL_READS. Its
# total_reads counts reads over roughly twenty escape genes, not the whole X,
# so the same number means something very different: the whole-chrX median is
# ~50 reads while the block's is a handful. Applying 30 here would discard
# almost every cell. Hence a separate constant and a separate directory --
# the two cutoffs are not comparable and should not share a folder name.
MIN_CEB_READS <- as.integer(Sys.getenv("MIN_CEB_READS", "5"))
stopifnot(!is.na(MIN_CEB_READS), MIN_CEB_READS >= 1)

CEB_DIR <- file.path("Allelic_ratio_results",
                     paste0("core_escape_cutoff_", MIN_CEB_READS))
dir.create(CEB_DIR, showWarnings = FALSE, recursive = TRUE)

# As with the whole-chrX table above, the core escape block's raw per-cell
# ratios are cutoff-independent -- 04 writes them before applying its filter,
# and 07 sweeps cutoffs over them. Top level, not inside a cutoff directory.
CEB_RATIOS_FILE <- "Allelic_ratio_results/core_escape_block_new_allelic_ratio_table.txt"
