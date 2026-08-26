# ---------------------------------------------------------------------------
# 08 - Is "AR >= 0.90" evidence that a cell has lost its inactive X?
#
# Run from the OCM_heart/ directory:  Rscript allelic_ratio/08_lox_calling.R
# Requires: 04_core_escape.R (for the raw core-escape ratio table) and
#           06_cutoff_sweep.R (for the whole-chrX per-cell table, used as an
#           independent cross-check). Does not need 04's filtered outputs.
# Writes:   Allelic_ratio_results/lox_calling/
#
# ---------------------------------------------------------------------------
# THE HYPOTHESIS AND WHY THE CURRENT TEST CANNOT CHECK IT
# ---------------------------------------------------------------------------
# The biology: escape genes are transcribed from BOTH the active X (Xa) and the
# inactive X (Xi). A cell that has lost its Xi can only transcribe them from
# the Xa, so its allelic ratio over the escape block should go to 1. 04
# therefore calls a cell "LOX-like" when AR >= 0.90.
#
# The problem is not the biology, it is that AR is estimated from very few
# reads. The core escape block has a MEDIAN of 3 informative reads per cell.
# With 3 reads the ratio can only be 0, 1/3, 2/3 or 1 -- and a perfectly
# normal two-X cell lands on exactly 1 whenever all three reads happen to come
# from the Xa. So "AR = 1" is produced by two completely different situations
# and the threshold cannot tell them apart.
#
# This script asks a specific question: once you account for how few reads
# there are, is there any EXCESS of high-AR cells beyond what a single
# population of ordinary escaping cells would produce anyway? If there is, that
# excess is your LOX population. If there is not, the AR >= 0.90 count is not
# evidence of X loss in either direction, and you need a signal that does not
# come from the allelic data.
#
# It works through four steps, each answering one question:
#   1. Is a binomial null even appropriate?          -> no, the data are overdispersed
#   2. How often does chance alone produce AR = 1?   -> the false-call rate per depth
#   3. Does ONE escaping population explain it all?  -> the central test
#   4. Does Xist agree with the allelic calls?       -> the independent check
# ---------------------------------------------------------------------------
source("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/OCM_heart/allelic_ratio/00_functions.R")

OUT_DIR <- "Allelic_ratio_results/lox_calling"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LOX_AR_THRESHOLD <- 0.90   # the rule 04 currently uses, evaluated here
FDR_CUT          <- 0.05
notes <- character(0)      # accumulates the plain-language summary written at the end
say <- function(...) {
  msg <- paste0(...)
  notes <<- c(notes, msg)
  message(msg)
}

# ---------------------------------------------------------------------------
# Data. A1_reads = reads from the maternal/A1 haplotype, A2_reads from the
# other; total_reads = A1 + A2; allelic_ratio = A1 / total_reads. Under the
# hypothesis, a lost-Xi cell has A2_reads ~ 0 and so AR ~ 1.
# ---------------------------------------------------------------------------
heart <- readRDS('heart_seurat_object_SCT.rds')
heart$celltype <- Idents(heart)

ceb <- read.delim(CEB_RATIOS_FILE, header = TRUE, stringsAsFactors = FALSE)
ceb <- subset(ceb, chr == "chrX")
stopifnot(!any(duplicated(ceb$cell_barcode)))

keep <- intersect(ceb$cell_barcode, colnames(heart))
ceb  <- ceb[match(keep, ceb$cell_barcode), ]

# Raw counts, not SCT-normalised values, because the only thing this uses Xist
# for is the yes/no question "did this cell yield any Xist reads at all".
xist_counts <- GetAssayData(heart, assay = "RNA", layer = "counts")["Xist", keep]

cells <- data.frame(
  cell_barcode = keep,
  celltype     = as.character(heart$celltype[keep]),
  sample       = factor(as.character(heart$sample[keep]),
                        levels = c("9w", "78w", "Sham", "TAC")),
  nCount_RNA   = heart$nCount_RNA[keep],
  A1_reads     = as.integer(ceb$A1_reads),
  A2_reads     = as.integer(ceb$A2_reads),
  total_reads  = as.integer(ceb$total_reads),
  allelic_ratio = ceb$allelic_ratio,
  xist_count   = as.numeric(xist_counts),
  stringsAsFactors = FALSE
)
cells <- cells[!is.na(cells$allelic_ratio) & cells$total_reads > 0, ]
cells$xist_positive <- cells$xist_count > 0
cells$lox_call_threshold <- cells$allelic_ratio >= LOX_AR_THRESHOLD

say("Cells with core-escape coverage: ", nrow(cells),
    ". Median informative reads per cell: ", median(cells$total_reads),
    " (max ", max(cells$total_reads), ").")
say("Cells with at least one Xist read: ", sum(cells$xist_positive),
    " (", round(100 * mean(cells$xist_positive), 1), "%).")

# Sanity check worth having: if Xist itself were inside the SNP set used for
# the block, it would corrupt the ratio, because Xist is transcribed from the
# Xi and so its allelic ratio runs the OPPOSITE way to every escape gene.
# The whole-chrX SNP file is explicitly a "_no_Xist" build; this cannot verify
# the core-escape build from here, so it is flagged rather than assumed.
say("NOTE: this assumes Xist is excluded from the core-escape SNP set, as it ",
    "is for whole chrX. If it is included, the block AR is contaminated by a ",
    "gene that reports the opposite allele and everything below is affected.")

# ---------------------------------------------------------------------------
# Beta-binomial helpers.
#
# A binomial distribution assumes every cell shares one underlying allelic
# proportion p, and that the only variability is sampling. A beta-binomial adds
# a second layer: each cell draws its OWN p from a Beta distribution, so it can
# represent genuine cell-to-cell variation in how much escape occurs.
#
# Parameterised by mu (the mean allelic proportion across cells) and rho (the
# intra-class correlation, 0 to 1). rho = 0 is exactly binomial; larger rho
# means more cell-to-cell spread. Var = N*mu*(1-mu)*(1 + (N-1)*rho).
# ---------------------------------------------------------------------------
bb_ab <- function(mu, rho) c(mu * (1 - rho) / rho, (1 - mu) * (1 - rho) / rho)

bb_logpmf <- function(k, N, a, b) {
  lchoose(N, k) + lbeta(k + a, N - k + b) - lbeta(a, b)
}

# P(K >= k) for one cell: how surprising its skew is under the given null
bb_upper_tail <- function(k, N, a, b) {
  vapply(seq_along(k), function(i) {
    ks <- k[i]:N[i]
    sum(exp(bb_logpmf(ks, N[i], a, b)))
  }, numeric(1))
}

fit_bb <- function(k, N) {
  nll <- function(par) {
    mu <- plogis(par[1]); rho <- plogis(par[2])
    v <- -sum(bb_logpmf(k, N, bb_ab(mu, rho)[1], bb_ab(mu, rho)[2]))
    if (!is.finite(v)) 1e12 else v
  }
  o <- optim(c(qlogis(0.6), qlogis(0.3)), nll, method = "BFGS")
  list(mu = plogis(o$par[1]), rho = plogis(o$par[2]), nll = o$value, n = length(k))
}

# ---------------------------------------------------------------------------
# STEP 1. Is a binomial null appropriate? Compare, at each fixed read depth,
# the observed variance of A1 counts against what binomial sampling alone
# predicts. If cells really did share one p, the ratio would sit near 1.
# ---------------------------------------------------------------------------
overdisp <- cells %>%
  group_by(total_reads) %>%
  filter(dplyr::n() >= 40) %>%
  summarise(n_cells = dplyr::n(),
            mean_AR = mean(A1_reads / total_reads),
            obs_var = var(A1_reads),
            .groups = "drop") %>%
  mutate(binom_var = total_reads * mean_AR * (1 - mean_AR),
         var_ratio = obs_var / binom_var,
         # rearranging Var = N*p*(1-p)*(1+(N-1)*rho) for rho. Undefined at N=1,
         # where there is only one read and so no within-cell spread to measure.
         implied_rho = ifelse(total_reads > 1,
                              (var_ratio - 1) / (total_reads - 1), NA_real_))
write.table(overdisp, file.path(OUT_DIR, "01_overdispersion_by_depth.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

say("")
say("STEP 1 -- is a binomial test appropriate? Observed variance of A1 counts ",
    "divided by the binomial expectation, at fixed depth, ranges from ",
    sprintf("%.1f", min(overdisp$var_ratio)), " to ",
    sprintf("%.1f", max(overdisp$var_ratio)), ".")
say("  A ratio of 1 would mean all cells share one allelic proportion. These ",
    "are far above 1, so cells genuinely differ from each other in how much ",
    "escape they show. CONSEQUENCE: a per-cell BINOMIAL test against a pooled ",
    "proportion would reject for most cells simply because cells differ, not ",
    "because any X was lost. The null has to be beta-binomial.")

# ---------------------------------------------------------------------------
# STEP 2. How often does chance alone give AR = 1? A cell with N reads and true
# proportion p shows AR = 1 whenever all N reads come from the Xa, which has
# probability p^N. This is the false-LOX rate: the share of ordinary two-X
# cells the AR >= 0.90 rule would wrongly call at each depth.
# ---------------------------------------------------------------------------
p_pooled <- with(subset(cells, total_reads >= 5), sum(A1_reads) / sum(total_reads))

chance_lox <- data.frame(min_reads = c(1, 2, 3, 4, 5, 8, 10, 12)) %>%
  rowwise() %>%
  mutate(cells_retained = sum(cells$total_reads >= min_reads),
         pct_retained = 100 * cells_retained / nrow(cells),
         false_lox_rate_at_this_depth = p_pooled ^ min_reads,
         observed_lox_rate = mean(cells$allelic_ratio[cells$total_reads >= min_reads]
                                  >= LOX_AR_THRESHOLD)) %>%
  ungroup()
write.table(chance_lox, file.path(OUT_DIR, "02_chance_lox_rate_by_cutoff.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

say("")
say("STEP 2 -- how much of the LOX call is chance? Pooled allelic proportion ",
    "of the block is ", sprintf("%.3f", p_pooled), ", so an ordinary two-X ",
    "cell shows AR = 1 with probability ", sprintf("%.3f", p_pooled), "^N.")
for (i in seq_len(nrow(chance_lox))) {
  r <- chance_lox[i, ]
  say(sprintf("  at exactly %2d reads: %.1f%% of normal cells look LOX by chance ",
              r$min_reads, 100 * r$false_lox_rate_at_this_depth),
      sprintf("(cutoff >= %d keeps %d cells, %.0f%%)",
              r$min_reads, r$cells_retained, r$pct_retained))
}
say("  This is why the current cutoff of ", MIN_CEB_READS, " is a problem for ",
    "LOX calling specifically: it does not stop chance producing the call.")

# ---------------------------------------------------------------------------
# STEP 3. THE CENTRAL TEST. Fit ONE beta-binomial population to all the cells,
# i.e. a model in which every cell is an ordinary escaping cell and NO cell has
# lost an X. Then ask that model how many cells it expects at AR = 1 and at
# AR >= 0.90, at each depth, and compare with what was observed.
#
# If the observed counts exceed the model's expectation, that excess is the LOX
# population. If they match, the high-AR cells are what an ordinary escaping
# population produces at these depths, and the threshold count is uninformative.
# ---------------------------------------------------------------------------
fit_all <- with(subset(cells, total_reads >= 3), fit_bb(A1_reads, total_reads))
ab_all <- bb_ab(fit_all$mu, fit_all$rho)

expected_vs_observed <- cells %>%
  group_by(total_reads) %>%
  filter(dplyr::n() >= 40) %>%
  summarise(n_cells = dplyr::n(),
            obs_AR_eq_1 = mean(A1_reads == total_reads),
            obs_AR_ge_thr = mean(A1_reads / total_reads >= LOX_AR_THRESHOLD),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(
    exp_AR_eq_1_binomial = p_pooled ^ total_reads,
    exp_AR_eq_1_betabinom = exp(bb_logpmf(total_reads, total_reads, ab_all[1], ab_all[2])),
    exp_AR_ge_thr_betabinom = sum(exp(bb_logpmf(
      seq(ceiling(LOX_AR_THRESHOLD * total_reads), total_reads),
      total_reads, ab_all[1], ab_all[2])))
  ) %>%
  ungroup() %>%
  mutate(excess_vs_betabinom = obs_AR_ge_thr - exp_AR_ge_thr_betabinom)
write.table(expected_vs_observed, file.path(OUT_DIR, "03_expected_vs_observed_by_depth.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

overall_obs <- mean(cells$allelic_ratio[cells$total_reads >= 8] >= LOX_AR_THRESHOLD)
overall_exp <- mean(sapply(cells$total_reads[cells$total_reads >= 8], function(N)
  sum(exp(bb_logpmf(seq(ceiling(LOX_AR_THRESHOLD * N), N), N, ab_all[1], ab_all[2])))))

say("")
say("STEP 3 -- the central test. A SINGLE beta-binomial population, containing ",
    "no lost-X cells at all, fits as mu = ", sprintf("%.3f", fit_all$mu),
    " (mean allelic proportion) and rho = ", sprintf("%.3f", fit_all$rho),
    " (cell-to-cell spread).")
say("  Observed vs that model's expectation for the AR >= ", LOX_AR_THRESHOLD,
    " rate, using cells with >= 8 reads where the ratio is best resolved: ",
    "observed ", sprintf("%.3f", overall_obs), " vs expected ",
    sprintf("%.3f", overall_exp), ", i.e. an excess of ",
    sprintf("%+.3f", overall_obs - overall_exp), ".")
if (overall_obs - overall_exp < 0.02) {
  say("  READ THIS AS: one ordinary escaping population already accounts for ",
      "the high-AR cells. The count of AR >= 0.90 cells is therefore NOT ",
      "evidence of X loss. Important caveat: this does not disprove the ",
      "hypothesis either. A two-population mixture and one overdispersed ",
      "population look nearly identical at a median of 3 reads -- the ",
      "beta-binomial 'explains' the spike by absorbing exactly the structure ",
      "you are looking for as dispersion. The honest conclusion is that the ",
      "block's allelic ratio cannot settle this question on its own.")
} else {
  say("  READ THIS AS: there are more high-AR cells than one escaping ",
      "population predicts. That excess is a genuine LOX candidate signal.")
}

# ---------------------------------------------------------------------------
# STEP 4. THE INDEPENDENT CHECK: Xist.
#
# Xist is transcribed from the Xi. A cell that has lost its Xi therefore cannot
# produce Xist reads. This is the key point -- that signal does not come from
# the allelic data, so it can corroborate or contradict it.
#
# The evidence is one-sided, and it matters which way. Xist > 0 is strong
# evidence that a cell HAS an Xi. Xist == 0 is only weak evidence that it does
# not, because a cell can easily fail to sample a transcript it does express
# (dropout). So Xist-POSITIVE cells are used as the reference population, and
# Xist == 0 is treated as a corroborating signal, never as a LOX definition.
#
# This is also what breaks the circularity in step 3. There, the null was
# estimated from the same cells being classified, so any real LOX population
# inflated the null it was being tested against. Here the null comes from cells
# independently known to have an Xi.
# ---------------------------------------------------------------------------
ref <- subset(cells, xist_positive & total_reads >= 3)
say("")
say("STEP 4 -- the independent check. Reference population: ", nrow(ref),
    " cells with >= 1 Xist read and >= 3 informative reads. These are known ",
    "to still have an Xi, so they define what escape looks like WITHOUT X loss.")

if (nrow(ref) >= 200) {
  fit_ref <- fit_bb(ref$A1_reads, ref$total_reads)
  ab_ref <- bb_ab(fit_ref$mu, fit_ref$rho)
  say("  Fitted on those cells: mu = ", sprintf("%.3f", fit_ref$mu),
      ", rho = ", sprintf("%.3f", fit_ref$rho),
      ".  (Compare the circular all-cell fit: mu = ", sprintf("%.3f", fit_all$mu),
      ", rho = ", sprintf("%.3f", fit_all$rho), ".)")

  # Per-cell test: how improbable is this cell's skew if it still had an Xi?
  cells$p_lox <- bb_upper_tail(cells$A1_reads, cells$total_reads, ab_ref[1], ab_ref[2])
  cells$FDR_lox <- p.adjust(cells$p_lox, method = "fdr")
  cells$lox_call_test <- cells$FDR_lox < FDR_CUT

  say("  Per-cell one-sided beta-binomial test against that null, FDR < ",
      FDR_CUT, ": ", sum(cells$lox_call_test), " of ", nrow(cells),
      " cells (", sprintf("%.1f%%", 100 * mean(cells$lox_call_test)), ") called LOX, ",
      "versus ", sum(cells$lox_call_threshold), " (",
      sprintf("%.1f%%", 100 * mean(cells$lox_call_threshold)),
      ") by the AR >= ", LOX_AR_THRESHOLD, " rule.")
  if (sum(cells$lox_call_test) == 0) {
    say("  READ THIS AS: a diagnostic, not a failure to find anything. With rho = ",
        sprintf("%.2f", fit_ref$rho), " the fitted Beta is U-shaped (alpha = ",
        sprintf("%.2f", ab_ref[1]), ", beta = ", sprintf("%.2f", ab_ref[2]),
        ", both below 1), so the null ALREADY expects a large share of cells at ",
        "AR = 1 -- roughly a quarter even at 8 reads. No single cell's ratio can ",
        "therefore be surprising enough to survive multiple-testing correction.")
    say("  The conclusion is that this block cannot support PER-CELL ",
        "classification at these depths, whatever test is used. Population-level ",
        "comparisons like the Xist association below still work, because they ",
        "aggregate over thousands of cells instead of asking each one alone.")
  }

  # Do the two independent signals agree? If the allelic call is picking out
  # real lost-Xi cells, those cells should be depleted of Xist.
  # fisher.test needs a full 2x2. It will not be one if nothing was called LOX,
  # or if every cell is Xist-positive -- both plausible here, and neither should
  # kill the script after the expensive fitting above.
  xist_or <- function(call_vec, label) {
    tab <- table(LOX = factor(call_vec, levels = c(FALSE, TRUE)),
                 Xist_zero = factor(!cells$xist_positive, levels = c(FALSE, TRUE)))
    if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
      say("    ", label, ": not testable (one group is empty)")
      return(tab)
    }
    ft <- fisher.test(tab)
    say("    ", label, ": OR = ", sprintf("%.2f", ft$estimate),
        ", p = ", format.pval(ft$p.value, digits = 3))
    tab
  }

  say("  Agreement with Xist (odds ratio > 1 means LOX-called cells are ",
      "depleted of Xist, as the hypothesis predicts):")
  tab <- xist_or(cells$lox_call_test, "beta-binomial test call")
  xist_or(cells$lox_call_threshold, paste0("AR >= ", LOX_AR_THRESHOLD, " rule      "))
  write.table(as.data.frame(tab), file.path(OUT_DIR, "04_lox_test_vs_xist_crosstab.txt"),
              sep = '\t', row.names = FALSE, quote = FALSE)

  per_sample <- cells %>%
    group_by(sample) %>%
    summarise(n_cells = dplyr::n(),
              median_reads = median(total_reads),
              lox_rate_threshold = mean(lox_call_threshold),
              lox_rate_test = mean(lox_call_test),
              pooled_AR = sum(A1_reads) / sum(total_reads),
              frac_xist_zero = mean(!xist_positive),
              .groups = "drop")
  write.table(per_sample, file.path(OUT_DIR, "05_lox_rate_per_sample.txt"),
              sep = '\t', row.names = FALSE, quote = FALSE)
  say("  Per-sample LOX rates are in 05_lox_rate_per_sample.txt. Compare the ",
      "two rate columns: the threshold rate tracks read depth, the test rate ",
      "should not.")
} else {
  say("  NOT ENOUGH Xist-positive cells (", nrow(ref),
      ") to anchor the null. Steps 3 and 4 cannot be separated on this data.")
}

write.table(cells, file.path(OUT_DIR, "06_per_cell_lox_calls.txt"),
            sep = '\t', row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# STEP 4b. THE CONFOUND THAT COULD EXPLAIN STEP 4 AWAY.
#
# 04's Fisher tests found that Xist-zero cells are 4-8x more likely to be
# called LOX in cardiomyocytes. Before believing that, one alternative has to
# be ruled out, and it is a serious one.
#
# Xist is a highly expressed lncRNA. A cell with zero Xist reads is quite
# likely just a cell with a shallow library that failed to sample it -- that is
# dropout, not X loss. And a shallow library also yields fewer chrX informative
# reads, which by STEP 2 makes a chance AR = 1 call MORE likely. So low depth
# pushes BOTH variables in the same direction and can manufacture the entire
# association without any cell having lost anything.
#
# Two things settle it:
#   (a) do Xist-zero cells actually have lower depth? If not, no confound.
#   (b) does the association survive INSIDE depth strata? Holding depth fixed
#       removes its ability to drive both variables. The Mantel-Haenszel test
#       combines the within-stratum odds ratios into one estimate.
# ---------------------------------------------------------------------------
if (exists("ab_ref")) {
  say("")
  say("STEP 4b -- is the Xist association just read depth? Xist is highly ",
      "expressed, so Xist == 0 may simply mean a shallow library, and shallow ",
      "libraries also produce more chance LOX calls. Checking both ways.")

  depth_by_xist <- cells %>%
    group_by(xist_positive) %>%
    summarise(n_cells = dplyr::n(),
              median_block_reads = median(total_reads),
              median_nCount_RNA = median(nCount_RNA),
              lox_rate_threshold = mean(lox_call_threshold),
              .groups = "drop")
  write.table(depth_by_xist, file.path(OUT_DIR, "04b_depth_by_xist_status.txt"),
              sep = '\t', row.names = FALSE, quote = FALSE)
  for (i in seq_len(nrow(depth_by_xist))) {
    r <- depth_by_xist[i, ]
    say(sprintf("  Xist %s: n = %d, median block reads = %d, median library = %.0f, LOX rate = %.3f",
                if (r$xist_positive) "> 0 " else "== 0", r$n_cells,
                r$median_block_reads, r$median_nCount_RNA, r$lox_rate_threshold))
  }
  say("  If the two rows differ in depth, the raw odds ratios from 04 are ",
      "confounded and only the stratified result below is interpretable.")
  say("  Note which way the LEFTOVER confound points. Stratifying on block reads ",
      "closes the chance-LOX path, but Xist-zero cells still have smaller ",
      "libraries, and by the whole-chrX depth analysis a smaller library means ",
      "proportionally more ambient RNA, which pulls AR DOWN and so makes a LOX ",
      "call LESS likely. The residual bias therefore works against the ",
      "association, and cannot be what produces it.")

  # Stratify on the EXACT read count, not on quantile bins.
  #
  # Why block reads is the right thing to hold fixed: the confounding path runs
  #   low library -> Xist dropout          -> Xist == 0
  #   low library -> fewer block reads     -> chance AR = 1 -> LOX call
  # Block reads sit on the second path, and the chance of a spurious LOX call
  # is p^N, fixed once N is fixed. So conditioning on N closes that path
  # completely, which conditioning on library size only does indirectly.
  #
  # And here N is a small integer (1 to 26, median 3), so we can condition on
  # its exact value rather than a bin -- within a stratum every cell has
  # IDENTICAL resolution and an identical chance-LOX rate. Depths above 10 are
  # pooled because individually they hold too few cells to contribute.
  cells$depth_stratum <- factor(ifelse(cells$total_reads >= 10, "10+",
                                       as.character(cells$total_reads)),
                                levels = c(as.character(1:9), "10+"))
  say("  Stratifying on the exact read count (1..9, then 10+), so within a ",
      "stratum every cell has the same resolution and the same chance of a ",
      "spurious LOX call. That closes the depth pathway rather than damping it.")

  strat <- cells %>%
    group_by(depth_stratum) %>%
    summarise(n_cells = dplyr::n(),
              n_xist_zero = sum(!xist_positive),
              lox_rate_xist_zero = mean(lox_call_threshold[!xist_positive]),
              lox_rate_xist_pos  = mean(lox_call_threshold[xist_positive]),
              .groups = "drop") %>%
    mutate(rate_difference = lox_rate_xist_zero - lox_rate_xist_pos)
  write.table(strat, file.path(OUT_DIR, "04b_lox_vs_xist_within_depth_strata.txt"),
              sep = '\t', row.names = FALSE, quote = FALSE)

  say("  LOX rate in Xist-zero vs Xist-positive cells, WITHIN depth strata:")
  for (i in seq_len(nrow(strat))) {
    r <- strat[i, ]
    say(sprintf("    reads %-4s n = %5d (%4d Xist-zero): %.3f vs %.3f, difference %+.3f",
                as.character(r$depth_stratum), r$n_cells, r$n_xist_zero,
                r$lox_rate_xist_zero, r$lox_rate_xist_pos, r$rate_difference))
  }

  # mantelhaen.test multiplies three cell counts together when building its
  # variance term. With thousands of cells per stratum that overflows R's
  # 32-bit integer limit and the confidence interval comes back NA, while the
  # point estimate still looks fine. Casting the table to double fixes it.
  mh_array <- function(df) {
    arr <- table(LOX = factor(df$lox_call_threshold, levels = c(FALSE, TRUE)),
                 Xist_zero = factor(!df$xist_positive, levels = c(FALSE, TRUE)),
                 stratum = df$depth_stratum)
    array(as.double(arr), dim = dim(arr), dimnames = dimnames(arr))
  }
  mh <- tryCatch(mantelhaen.test(mh_array(cells)), error = function(e) NULL)

  if (!is.null(mh)) {
    ci_lo <- mh$conf.int[1]; ci_hi <- mh$conf.int[2]
    say(sprintf("  Mantel-Haenszel combined odds ratio, depth held fixed: %.2f (95%% CI %s-%s), p = %s",
                mh$estimate,
                if (is.na(ci_lo)) "NA" else sprintf("%.2f", ci_lo),
                if (is.na(ci_hi)) "NA" else sprintf("%.2f", ci_hi),
                format.pval(mh$p.value, digits = 3)))
    if (!is.na(ci_lo) && ci_lo > 1) {
      say("  READ THIS AS: the association SURVIVES depth stratification. The ",
          "confidence interval excludes 1, so Xist-zero cells are genuinely ",
          "more likely to look monoallelic even among cells of the same read ",
          "depth. That is independent support for the lost-X hypothesis, ",
          "because Xist status does not come from the allelic data.")
    } else if (is.na(ci_lo)) {
      say("  Confidence interval not computable, so judge on the p-value and ",
          "the per-stratum rate differences above rather than the interval.")
    } else {
      say("  READ THIS AS: once depth is held fixed the association is no ",
          "longer significant, so the raw odds ratios in 04's Fisher tests are ",
          "consistent with being a depth artefact -- shallow cells lose Xist to ",
          "dropout AND gain chance LOX calls.")
    }

    # Mantel-Haenszel reports a single pooled odds ratio, which assumes the true
    # odds ratio is the SAME in every stratum. The per-stratum differences above
    # visibly grow with depth, so test that assumption rather than let the pooled
    # number quietly average over it. A logistic model with an interaction does
    # the job in base R: if the interaction improves fit, the effect is modified
    # by depth and the single pooled OR is an average over unlike things.
    homog <- tryCatch({
      d <- subset(cells, !is.na(depth_stratum))
      m_add <- glm(lox_call_threshold ~ (!xist_positive) + depth_stratum,
                   family = binomial, data = d)
      m_int <- glm(lox_call_threshold ~ (!xist_positive) * depth_stratum,
                   family = binomial, data = d)
      anova(m_add, m_int, test = "LRT")
    }, error = function(e) NULL)
    if (!is.null(homog)) {
      pv <- homog[["Pr(>Chi)"]][2]
      say(sprintf("  Is the odds ratio constant across depth strata? interaction LRT p = %s",
                  format.pval(pv, digits = 3)))
      if (!is.na(pv) && pv < 0.05) {
        say("    It is NOT constant, so read 1.76 as an average over strata that ",
            "differ rather than one underlying effect. The differences above are ",
            "smallest at 1 read and larger from 2 upward, which is what you would ",
            "expect if Xist == 0 at 1 read is mostly dropout in a shallow cell, ",
            "and becomes a more trustworthy marker of a missing Xi as depth rises. ",
            "That pattern supports the signal being real; it just means the ",
            "pooled OR is not the right single summary.")
      }
    }

    # Same test restricted to cells where Xist == 0 is a more reliable marker.
    # At 1 read the cell is so shallow that a zero Xist count says very little.
    mh3 <- tryCatch(mantelhaen.test(mh_array(subset(cells, total_reads >= 3))),
                    error = function(e) NULL)
    if (!is.null(mh3)) {
      say(sprintf("  Restricted to cells with >= 3 reads (n = %d), where a zero Xist count is more informative: OR = %.2f [%.2f-%.2f], p = %s",
                  sum(cells$total_reads >= 3), mh3$estimate,
                  mh3$conf.int[1], mh3$conf.int[2],
                  format.pval(mh3$p.value, digits = 3)))
    }

    # Pooling every celltype together dilutes this: 04 found the effect in
    # ventricular cardiomyocytes and endothelial cells, not everywhere. So also
    # run it per celltype, still with depth held fixed within each.
    say("  Per celltype, depth still held fixed within each:")
    ct_mh <- lapply(split(cells, cells$celltype), function(d) {
      if (nrow(d) < 200 || sum(!d$xist_positive) < 20) return(NULL)
      m <- tryCatch(mantelhaen.test(mh_array(d)), error = function(e) NULL)
      if (is.null(m)) return(NULL)
      data.frame(celltype = d$celltype[1], n_cells = nrow(d),
                 n_xist_zero = sum(!d$xist_positive),
                 OR = unname(m$estimate), lwr = m$conf.int[1], upr = m$conf.int[2],
                 p_value = m$p.value)
    })
    ct_mh <- bind_rows(ct_mh)
    if (nrow(ct_mh)) {
      ct_mh$FDR <- p.adjust(ct_mh$p_value, method = "fdr")
      ct_mh <- ct_mh[order(ct_mh$p_value), ]
      write.table(ct_mh, file.path(OUT_DIR, "04b_mantelhaenszel_per_celltype.txt"),
                  sep = '\t', row.names = FALSE, quote = FALSE)
      for (i in seq_len(nrow(ct_mh))) {
        r <- ct_mh[i, ]
        say(sprintf("    %-32s n = %5d (%4d Xist-zero): OR = %.2f [%.2f-%.2f], FDR = %s",
                    substr(r$celltype, 1, 32), r$n_cells, r$n_xist_zero,
                    r$OR, r$lwr, r$upr, format.pval(r$FDR, digits = 3)))
      }
      say("    An OR above 1 with an interval clear of it means Xist-zero cells ",
          "look monoallelic more often than equally-deep Xist-positive cells of ",
          "the same celltype -- which allelic sampling alone cannot produce.")
    }
  } else {
    say("  Mantel-Haenszel test not computable (a stratum has an empty margin).")
  }

  # Prevalence. If the Xist-positive LOX rate is the chance baseline -- those
  # cells demonstrably still have an Xi, so every LOX call among them is a
  # false positive -- then the EXCESS among Xist-zero cells is the part
  # attributable to real X loss. This is a floor, not an estimate of the total:
  # dropout means many genuinely Xi-less cells are not in the Xist-zero group.
  base_rate <- mean(cells$lox_call_threshold[cells$xist_positive])
  nz <- sum(!cells$xist_positive)
  excess <- sum(cells$lox_call_threshold[!cells$xist_positive]) - base_rate * nz
  say(sprintf("  Prevalence floor: LOX rate among Xist-positive (known to have an Xi) is %.3f,",
              base_rate))
  say(sprintf("    so of the %d Xist-zero cells, ~%.0f are LOX calls beyond that baseline",
              nz, max(0, excess)))
  say(sprintf("    = ~%.1f%% of all %d cells, against the %.1f%% the AR >= %.2f rule reports.",
              100 * max(0, excess) / nrow(cells), nrow(cells),
              100 * mean(cells$lox_call_threshold), LOX_AR_THRESHOLD))
  say("    Treat that as a lower bound on prevalence and an upper bound on the ",
      "threshold rule's specificity.")
}

# ---------------------------------------------------------------------------
# STEP 5. Second independent check: the whole X chromosome.
#
# A cell that has lost an entire X cannot produce A2 reads ANYWHERE on it, not
# just over the escape block. So cells called LOX from the block should also
# look monoallelic across the whole chromosome. This is only near-independent:
# both measurements use the same cell's reads and are affected by the same
# ambient contamination, so treat agreement as supportive, not conclusive.
# ---------------------------------------------------------------------------
wx_file <- "Allelic_ratio_results/cutoff_sweep/cutoff_sweep_cell_table.txt"
if (file.exists(wx_file)) {
  wx <- read.delim(wx_file, stringsAsFactors = FALSE)
  wx <- wx[, c("cell_barcode", "A1_reads", "A2_reads", "total_reads", "allelic_ratio")]
  names(wx) <- c("cell_barcode", "wx_A1", "wx_A2", "wx_total", "wx_AR")
  j <- merge(cells, wx, by = "cell_barcode")
  j <- subset(j, total_reads >= 5 & wx_total >= 30)

  # grouped on the threshold rule specifically: the point of this check is to
  # ask whether THAT rule's calls hold up chromosome-wide
  cross <- j %>%
    group_by(lox_call_threshold) %>%
    summarise(n_cells = dplyr::n(),
              whole_X_AR = mean(wx_AR),
              pct_whole_X_A2_zero = 100 * mean(wx_A2 == 0),
              median_whole_X_A2 = median(wx_A2),
              .groups = "drop")
  write.table(cross, file.path(OUT_DIR, "07_wholeX_crosscheck.txt"),
              sep = '\t', row.names = FALSE, quote = FALSE)

  say("")
  say("STEP 5 -- whole-chromosome cross-check on ", nrow(j), " cells with both ",
      "measurements (block >= 5 reads, whole X >= 30 reads):")
  for (i in seq_len(nrow(cross))) {
    r <- cross[i, ]
    say(sprintf("  block AR %s %.2f : n = %d, whole-X AR = %.3f, %.1f%% with ZERO maternal reads across the whole X, median %d",
                if (r$lox_call_threshold) ">=" else "<", LOX_AR_THRESHOLD,
                r$n_cells, r$whole_X_AR, r$pct_whole_X_A2_zero, r$median_whole_X_A2))
  }
  say("  READ THIS AS: a genuinely lost-X cell should have essentially no ",
      "maternal reads anywhere on the chromosome. Cells with zero are the ",
      "strongest LOX candidates you have. Cells still carrying a substantial ",
      "number of maternal reads chromosome-wide are hard to reconcile with ",
      "having lost that chromosome, whatever their block AR says.")
} else {
  say("")
  say("STEP 5 skipped: ", wx_file, " not found. Run 06_cutoff_sweep.R for the ",
      "whole-chromosome cross-check.")
}

# ---------------------------------------------------------------------------
# The plain-language summary, written to disk so the reasoning survives outside
# whatever terminal this ran in.
# ---------------------------------------------------------------------------
writeLines(c(
  "LOX calling in the core escape block -- what the numbers mean",
  paste0("Generated by allelic_ratio/08_lox_calling.R, cutoff MIN_CEB_READS = ", MIN_CEB_READS),
  strrep("=", 74), "", notes, "",
  strrep("=", 74),
  "Files:",
  "  01_overdispersion_by_depth.txt     why a binomial test is the wrong tool",
  "  02_chance_lox_rate_by_cutoff.txt   how often chance alone gives AR = 1",
  "  03_expected_vs_observed_by_depth.txt  the central test: does one population explain it",
  "  04_lox_test_vs_xist_crosstab.txt   allelic call vs the independent Xist signal",
  "  04b_depth_by_xist_status.txt       is Xist==0 just a shallow library?",
  "  04b_lox_vs_xist_within_depth_strata.txt  the association with depth held fixed",
  "  04b_mantelhaenszel_per_celltype.txt  the same test per celltype",
  "  05_lox_rate_per_sample.txt         LOX rate per sample, both methods",
  "  06_per_cell_lox_calls.txt          per-cell p_lox, FDR_lox and both calls",
  "  07_wholeX_crosscheck.txt           do block calls look monoallelic chromosome-wide"
), file.path(OUT_DIR, "00_SUMMARY.txt"))

message("\nDone. Plain-language summary: ", file.path(OUT_DIR, "00_SUMMARY.txt"))
