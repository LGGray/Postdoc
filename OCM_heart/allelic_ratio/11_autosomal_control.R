# ---------------------------------------------------------------------------
# 11 - Whole-chromosome allelic ratio for chrX AND the autosomes.
#
# The question 02 cannot answer on its own: how much of a cell's chrX allelic
# ratio is X inactivation, and how much is the assay?
#
# We map F1 B6xCAST reads to a B6 reference, so CAST-allele reads are slightly
# less mappable and every ratio sits above its true value in the B6 direction.
# There is also whatever library skew and duplication structure the sample
# carries. All of that applies to the autosomes too, where there is no
# monoallelic biology to explain it -- so the autosomal ratio in the same cell,
# through the same pipeline, is an empirical null that already contains the
# bias. chrX is only interesting where it exceeds that null.
#
# Requires 10_build_ratio_table.R to have been run for the same RESULTS_ROOT,
# with a genome-wide annotation. Run from OCM_heart/, in seurat_env:
#
#   RESULTS_ROOT=Allelic_ratio_results_dedup Rscript allelic_ratio/11_autosomal_control.R
# ---------------------------------------------------------------------------
source("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/OCM_heart/allelic_ratio/00_functions.R")

OUT_DIR <- file.path(CUTOFF_DIR, "autosomal_control")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

tbl <- read.table(ALLELIC_RATIOS_FILE, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE)
stopifnot(all(c("chr", "cell_barcode", "A1_reads", "A2_reads") %in% names(tbl)))
if (!"autosomal" %in% tbl$chr) {
  stop("No `autosomal` rows in ", ALLELIC_RATIOS_FILE,
       " -- rebuild it with 10_build_ratio_table.R from a genome-wide tree.")
}

x <- tbl[tbl$chr == "chrX", ]
a <- tbl[tbl$chr == "autosomal", ]
message("cells with chrX: ", nrow(x), " | with autosomal: ", nrow(a))

# ---------------------------------------------------------------------------
# The cutoff is applied to each arm on its own count. Gating the autosomes at
# MIN_TOTAL_READS chrX reads would select cells on the very quantity being
# controlled for; gating chrX on autosomal depth would do the reverse. Both
# arms need enough reads to estimate a ratio, and that is all the filter says.
# ---------------------------------------------------------------------------
x <- x[x$total_reads >= MIN_TOTAL_READS, ]
a <- a[a$total_reads >= MIN_TOTAL_READS, ]

paired <- inner_join(
  x[, c("cell_barcode", "sample", "A1_reads", "A2_reads", "total_reads", "ar_a1", "ar_dom")],
  a[, c("cell_barcode", "sample", "A1_reads", "A2_reads", "total_reads", "ar_a1", "ar_dom")],
  by = c("cell_barcode", "sample"), suffix = c("_x", "_a"))
stopifnot(nrow(paired) > 0)
message("cells with both arms >= ", MIN_TOTAL_READS, " reads: ", nrow(paired))

# ar_dom is bounded below at 0.5 by construction, so an autosomal ar_dom of,
# say, 0.62 does NOT mean 62:38 skew -- with few reads the max of two roughly
# equal counts is above 0.5 on average purely from sampling. Report the excess
# over the autosomal value in the same cell, which cancels that floor to first
# order, and report the directional ratio separately where sign matters.
paired$ar_excess <- paired$ar_dom_x - paired$ar_dom_a

# The skew-free expectation for ar_dom at this depth, so the autosomal value can
# be read against something. E[max(k, n-k)/n] for k ~ Binom(n, 0.5).
expected_ar_dom <- function(n) {
  vapply(n, function(nn) {
    if (!is.finite(nn) || nn < 1) return(NA_real_)
    k <- 0:nn
    sum(dbinom(k, nn, 0.5) * pmax(k, nn - k) / nn)
  }, numeric(1))
}
paired$ar_dom_a_null <- expected_ar_dom(paired$total_reads_a)
paired$ar_dom_x_null <- expected_ar_dom(paired$total_reads_x)

per_sample <- paired %>%
  group_by(sample) %>%
  summarise(n              = n(),
            med_reads_x    = median(total_reads_x),
            med_reads_a    = median(total_reads_a),
            med_ar_x       = median(ar_dom_x),
            med_ar_a       = median(ar_dom_a),
            med_null_a     = median(ar_dom_a_null),
            med_excess     = median(ar_excess),
            med_dir_x      = median(ar_a1_x),
            med_dir_a      = median(ar_a1_a),
            frac_x_mono    = mean(ar_dom_x >= 0.90),
            frac_a_mono    = mean(ar_dom_a >= 0.90),
            .groups = "drop")
print(as.data.frame(per_sample))

# frac_a_mono is the false monoallelic rate: cells the pipeline would call
# monoallelic on autosomes, where no such thing exists. Whatever it is, chrX
# calls at the same depth carry at least that much of it.
message(sprintf("\nAutosomal cells at ar_dom >= 0.90 (false monoallelic rate): %.1f%%",
                100 * mean(paired$ar_dom_a >= 0.90)))
message(sprintf("chrX cells at ar_dom >= 0.90: %.1f%%",
                100 * mean(paired$ar_dom_x >= 0.90)))

# Paired test per sample: is chrX skewed beyond its own cell's autosomes?
# Wilcoxon signed-rank on the within-cell difference, which needs no assumption
# about the shape of the bias, only that it is shared between the two arms.
tests <- bind_rows(lapply(split(paired, paired$sample), function(d) {
  if (nrow(d) < 10) return(NULL)
  w <- wilcox.test(d$ar_dom_x, d$ar_dom_a, paired = TRUE)
  data.frame(sample = d$sample[1], n = nrow(d),
             median_excess = median(d$ar_excess),
             p = w$p.value, stringsAsFactors = FALSE)
}))
if (nrow(tests)) {
  tests$fdr <- p.adjust(tests$p, method = "BH")
  tests$stars <- fdr_to_stars(tests$fdr)
  print(tests)
  write.table(tests, file.path(OUT_DIR, "chrX_vs_autosomal_paired_test.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

write.table(paired, file.path(OUT_DIR, "paired_chrX_autosomal.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(per_sample, file.path(OUT_DIR, "per_sample_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

long <- bind_rows(
  transform(paired[, c("sample", "cell_barcode")], arm = "chrX",
            ar = paired$ar_dom_x, depth = paired$total_reads_x),
  transform(paired[, c("sample", "cell_barcode")], arm = "autosomal",
            ar = paired$ar_dom_a, depth = paired$total_reads_a))

p_v <- ggplot(long, aes(arm, ar, fill = arm)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, fill = "white") +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "red") +
  facet_wrap(~ sample) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(title = "Whole-chromosome allelic ratio, chrX vs its own autosomal control",
       subtitle = sprintf("dominant-allele fraction, both arms >= %d reads; red = the 0.90 monoallelic line",
                          MIN_TOTAL_READS),
       x = NULL, y = "max(A1, A2) / total") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "chrX_vs_autosomal_violin.pdf"), p_v, width = 9, height = 6)

p_s <- ggplot(paired, aes(ar_dom_a, ar_dom_x)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_point(aes(colour = total_reads_x), alpha = 0.4, size = 0.8) +
  scale_colour_viridis_c(trans = "log10", name = "chrX reads") +
  facet_wrap(~ sample) +
  coord_fixed(xlim = c(0.5, 1), ylim = c(0.5, 1)) +
  labs(title = "Per cell: chrX ratio against the same cell's autosomal ratio",
       subtitle = "points on the diagonal are explained by the technical null",
       x = "autosomal", y = "chrX") +
  theme_bw()
ggsave(file.path(OUT_DIR, "chrX_vs_autosomal_scatter.pdf"), p_s, width = 9, height = 8)

# Depth is the confounder 06 already found on chrX (mean AR rising with read
# count). If the autosomal arm shows the same curve, it is the assay, not the X.
p_d <- ggplot(long, aes(depth, ar, colour = arm)) +
  geom_point(alpha = 0.2, size = 0.6) +
  geom_smooth(method = "loess", se = TRUE, span = 0.9) +
  scale_x_log10() +
  facet_wrap(~ sample) +
  labs(title = "Allelic ratio against depth, both arms",
       subtitle = "a shared downward trend is the sampling floor of a bounded statistic, not biology",
       x = "reads (log10)", y = "max(A1, A2) / total") +
  theme_bw()
ggsave(file.path(OUT_DIR, "AR_vs_depth_both_arms.pdf"), p_d, width = 10, height = 6)

message("Wrote ", OUT_DIR)
