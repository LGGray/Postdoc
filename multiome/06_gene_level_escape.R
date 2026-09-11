# ---------------------------------------------------------------------------
# 06 - Per-gene escape from the gene-level pseudobulk Allelome.PRO2 output.
#
# Reads what slurm/multiome_allelome_genelevel.slurm wrote under
# adult_aged_multiome/allelome_gene_level/<sample>/, which is one row per GENE
# rather than the one row per CHROMOSOME that 05_allelome_plots.R consumes.
#
# WHAT THIS ANSWERS THAT 05 CANNOT. 05 reports a single pooled chrX number per
# cell type. That number is an expression-weighted average over every chrX gene,
# so it cannot say which genes escape, cannot be compared against a published
# escape list, and cannot support anything positional. It is also the same
# pooling artefact as pooling reads across cell types, one level up.
#
# COLUMN CONVENTION, stated because this repo has been bitten by it. chrX
# carries THREE different quantities that all get called an allelic ratio:
#     ar_b6       A1 / total,  0-1, the B6 fraction. B6 is the ACTIVE X here,
#                 so ar_b6 near 1 is monoallelic and lower means escape. This
#                 is the OCM convention (ap2_pseudobulk_common.R, MONO_AR).
#     escape_cast A2 / total = 1 - ar_b6, the CAST fraction. This is what
#                 05_allelome_plots.R calls `escape`.
#     ar_b6_corr  ar_b6 after the mapping-bias correction below.
# Both raw columns are kept in the exported table. Do not collapse them.
#
#   Rscript multiome/06_gene_level_escape.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(tibble)
})

CL     <- Sys.getenv("CLUSTER_MOUNT", "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2")
ROOT   <- Sys.getenv("AP2_ROOT", file.path(CL, "adult_aged_multiome", "allelome_gene_level"))
OUT    <- Sys.getenv("FIG_OUT", file.path(CL, "adult_aged_multiome", "figures_gene_level"))
COMMON <- Sys.getenv("AP2_COMMON", file.path(CL, "Postdoc", "OCM_heart", "ap2_pseudobulk_common.R"))
SAMPLES   <- strsplit(Sys.getenv("GROUPS", "9w,78w"), ",")[[1]]
MIN_READS <- as.integer(Sys.getenv("MIN_READS", "20"))   # SNP-overlapping reads per gene
MIN_SNP   <- as.integer(Sys.getenv("MIN_SNP", "1"))      # informative SNPs per gene
PRIOR_ESCAPE <- 0.127     # snRNA/spatial estimate the pooled result should reproduce
AUTOSOMES <- paste0("chr", 1:19)
ANNOT     <- "annotation_us_mm39_gene_level.bed"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n", sep = "")

source(file.path(CL, "Postdoc", "multiome", "00_helpers.R"))

# ---------------------------------------------------------------------------
# loader - ALL chromosomes, not just chrX
# ---------------------------------------------------------------------------
# The autosomes are not incidental here. They are the only honest null for
# chrX: same B6 reference, same mapping bias, same repeat content, no
# monoallelic biology. Keeping them in the same table means the bias is
# measured under exactly the filters the chrX numbers were measured under.
load_gene_level <- function(root = ROOT, samples = SAMPLES) {
  rows <- list()
  for (s in samples) {
    dirs <- list.files(file.path(root, s), pattern = "bed_1$", full.names = TRUE)
    say("%s: %d Allelome.PRO2 output directories", s, length(dirs))
    for (d in dirs) {
      f <- file.path(d, "locus_table.txt")
      if (!file.exists(f)) next
      grp <- sub(paste0("_", ANNOT, "_1$"), "", sub("^[0-9]{4}_[0-9]{2}_[0-9]{2}_", "", basename(d)))
      tab <- tryCatch(read_tsv(f, show_col_types = FALSE,
                               col_types = cols(chr = col_character(), name = col_character(),
                                                .default = col_double())),
                      error = function(e) NULL)
      if (is.null(tab) || !nrow(tab)) next
      rows[[length(rows) + 1]] <- tab %>%
        transmute(sample = s, group = grp, gene = name, chr, start,
                  A1 = A1_reads, A2 = A2_reads, total = total_reads)
    }
  }
  if (!length(rows)) stop("no gene-level output under ", root,
                          " - run slurm/multiome_allelome_genelevel.slurm first")
  bind_rows(rows)
}

dat <- load_gene_level()

# Recover the full cell type from the sinto group name by INVERTING the exact
# transformation 03_signac_joint.R applied, rather than guessing at it. That
# script writes `gsub("[^A-Za-z0-9]+", "_", celltype_provisional)`, which is not
# a simple underscore-for-space swap: "Cardiomyocytes (stressed)" gains a
# TRAILING underscore from the close bracket, and the " - " in
# "Pericytes - Smooth muscle cells" collapses to one. Applying the same gsub to
# the known cell-type names builds the map that cannot drift from 03.
group_to_full <- function(g) {
  m <- setNames(names(CELLTYPE_SHORT),
                gsub("[^A-Za-z0-9]+", "_", names(CELLTYPE_SHORT)))
  out <- unname(m[g])                       # unname: see set_meta in 00_helpers.R
  ifelse(is.na(out), g, out)                # unmapped groups pass through
}

# The whole-sample pseudobulk BAM is named <id>_gex_dedup, so it arrives as a
# group alongside the cell types. It is the SAME READS as the cell types summed,
# so it must never be analysed in the same table as them - that would double
# count every read. Split it out here and use it only for the pooled gate and
# the bias estimate.
is_whole <- grepl("_gex_dedup$", dat$group)
whole <- dat %>% filter(is_whole)
ct    <- dat %>% filter(!is_whole)
say("%d gene rows across %d cell types, plus %d whole-sample rows",
    nrow(ct), dplyr::n_distinct(ct$group), nrow(whole))

# ---- informative SNPs per gene, the other confounder of a positional trend --
# A gene with few usable SNPs is measured worse, not necessarily more
# biallelically. Written by the slurm job; absent is a warning, not an error.
# V4 is the gene name in these beds (same layout 01_setup.R relies on for Xist)
# and bedtools intersect -c appends its count as the LAST column.
snpf <- file.path(ROOT, "gene_level_snp_count.bed")
if (file.exists(snpf)) {
  sc <- read.delim(snpf, header = FALSE, stringsAsFactors = FALSE)
  sc <- tibble(gene = sc[[4]], n_snp = as.integer(sc[[ncol(sc)]])) %>%
    group_by(gene) %>% summarise(n_snp = max(n_snp), .groups = "drop")
  ct    <- left_join(ct,    sc, by = "gene")
  whole <- left_join(whole, sc, by = "gene")
  say("joined informative SNP counts for %d genes", nrow(sc))
} else {
  say("NOTE: no gene_level_snp_count.bed found - SNP density not controlled for")
  ct$n_snp <- NA_integer_; whole$n_snp <- NA_integer_
}

# ---------------------------------------------------------------------------
# mapping bias, MEASURED AND THEN APPLIED
# ---------------------------------------------------------------------------
# Reads are aligned to a B6 reference, so CAST reads align slightly worse and
# every ratio is pulled toward B6. 05_allelome_plots.R measures this offset and
# prints it, but reports uncorrected chrX numbers. At chromosome level that was
# a rounding matter; at gene level it shifts which genes cross an escape cutoff,
# so it is corrected here.
#
# Model: CAST reads survive at rate (1 - lambda) relative to B6. Then
#     p_obs = p / (p + (1 - p)(1 - lambda))
# On autosomes the truth is p = 0.5, so p_auto = 1 / (2 - lambda), giving
#     lambda = 2 - 1 / p_auto
# and inverting the first equation recovers the corrected ratio
#     p = p_obs (1 - lambda) / (1 - p_obs lambda)
# Check: lambda = 0 leaves p unchanged, and p_obs = p_auto returns exactly 0.5.
auto <- whole %>% filter(chr %in% AUTOSOMES) %>% summarise(A1 = sum(A1), A2 = sum(A2))
p_auto <- auto$A1 / (auto$A1 + auto$A2)
LAMBDA <- 2 - 1 / p_auto
debias <- function(p_obs) p_obs * (1 - LAMBDA) / (1 - p_obs * LAMBDA)

xw <- whole %>% filter(chr == "chrX") %>% summarise(A1 = sum(A1), A2 = sum(A2))
p_x <- xw$A1 / (xw$A1 + xw$A2)

say("")
say("================ pooled gate and mapping bias ================")
say("autosomal B6 fraction      %.4f  (%d / %d reads)", p_auto, auto$A1, auto$A1 + auto$A2)
say("  -> CAST read loss lambda %.4f", LAMBDA)
say("pooled chrX B6 fraction    %.4f  (%d / %d reads)", p_x, xw$A1, xw$A1 + xw$A2)
say("pooled chrX escape (CAST)  %.4f raw, %.4f bias-corrected", 1 - p_x, 1 - debias(p_x))
say("snRNA/spatial prior        %.4f", PRIOR_ESCAPE)
say("==============================================================")
say("")
if (1 - debias(p_x) < 0.5 * PRIOR_ESCAPE || 1 - debias(p_x) > 2 * PRIOR_ESCAPE) {
  say("WARNING: the pooled gene-level escape is more than a factor of two from")
  say("  the snRNA/spatial prior. ANALYSIS_PLAN.md step 4 says to stop and find")
  say("  out why before building on it. Candidates: the Xic mask extent, the SNP")
  say("  build, or the MAPQ filter.")
}

# ---------------------------------------------------------------------------
# per-gene table
# ---------------------------------------------------------------------------
# MIN_READS is a real gate, not a formality. A 2-read gene has an allelic ratio
# of 0, 0.5 or 1 by construction and carries no information about escape, so
# without it the escape-rate denominator fills with noise.
gene <- ct %>%
  filter(total >= MIN_READS, is.na(n_snp) | n_snp >= MIN_SNP) %>%
  mutate(ar_b6       = A1 / total,
         ar_b6_corr  = debias(ar_b6),
         escape_cast = 1 - ar_b6_corr,
         region      = ifelse(chr == "chrX", "chrX",
                              ifelse(chr %in% AUTOSOMES, "autosomal", NA_character_)),
         celltype    = short_labels(group_to_full(group))) %>%
  filter(!is.na(region))

# Per-gene binomial CI. This is defensible where the pooled-chrX Wilson interval
# in 05 was not: within ONE gene the reads are close to exchangeable, whereas
# pooling thousands of genes with wildly different ratios makes a binomial
# interval far too narrow for the quantity it is being read as.
ci <- t(mapply(function(a2, n) {
  if (n == 0) return(c(NA, NA))
  as.numeric(binom.test(a2, n)$conf.int)
}, gene$A2, gene$total))
gene$escape_lo <- debias(1 - ci[, 2]) ; gene$escape_hi <- debias(1 - ci[, 1])

write_tsv(gene, file.path(OUT, "gene_level_escape.tsv"))
say("wrote gene_level_escape.tsv (%d gene x cell type rows, %d on chrX)",
    nrow(gene), sum(gene$region == "chrX"))

gx <- gene %>% filter(region == "chrX")

# Pooled over cell types, which is the level published escape calls are made at
# and the level the chromosome-wide figures below read at. Reads are summed from
# the raw counts and the ratio recomputed - averaging per-cell-type ratios would
# weight a 30-read cell type the same as a 3000-read one.
# THE THRESHOLD IS APPLIED AFTER POOLING, NOT BEFORE. Pooling the already
# filtered per-cell-type table would silently drop any gene whose coverage is
# spread thin but adequate in total - 15 reads in each of five cell types is 75
# informative reads, and is exactly the kind of broadly expressed gene an escape
# analysis must not lose. So this pools the RAW rows and thresholds the sum.
pool_ct <- function(d) {
  d %>% group_by(sample, gene, chr, start) %>%
    summarise(A1 = sum(A1), A2 = sum(A2), n_snp = max(n_snp), n_ct = n(),
              .groups = "drop") %>%
    mutate(total = A1 + A2) %>%
    filter(total >= MIN_READS, is.na(n_snp) | n_snp >= MIN_SNP) %>%
    mutate(ar_b6 = A1 / total, ar_b6_corr = debias(ar_b6),
           escape_cast = 1 - ar_b6_corr,
           region = ifelse(chr == "chrX", "chrX", "autosomal"))
}
ct_typed  <- ct %>% filter(chr == "chrX" | chr %in% AUTOSOMES)
gene_pool <- pool_ct(ct_typed)
gx_pool   <- gene_pool %>% filter(region == "chrX")
write_tsv(gx_pool, file.path(OUT, "gene_level_escape_pooled.tsv"))

say("")
say("chrX genes passing >= %d reads, per sample x cell type:", MIN_READS)
print(as.data.frame(gx %>% count(sample, celltype) %>% arrange(sample, desc(n))), row.names = FALSE)

# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
ESCAPE_GENES <- c("Kdm5c","Kdm6a","Ddx3x","Eif2s3x","Ftx","Jpx",
                  "Pbdc1","Utp14a","Akap17a","Sts")

# ggrepel is not in every env this runs in, and a missing suggested package must
# not cost the whole figure. Falls back to plain labels.
label_layer <- function(d) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = d, aes(label = gene), size = 3,
                             max.overlaps = 20, show.legend = FALSE)
  } else {
    geom_text(data = d, aes(label = gene), size = 3, vjust = -0.8, show.legend = FALSE)
  }
}

dev_open(file.path(OUT, "gene_level_escape.pdf"), width = 12, height = 7.5)

# 1. the figure chromosome-level output could never produce
print(
  ggplot(gx_pool, aes(start / 1e6, escape_cast, colour = sample)) +
    geom_hline(yintercept = PRIOR_ESCAPE, linetype = "dashed", colour = "grey40") +
    geom_point(alpha = 0.6, size = 1.2) +
    label_layer(gx_pool %>% filter(gene %in% ESCAPE_GENES)) +
    scale_colour_manual(values = setNames(OKABE_ITO[1:2], SAMPLES)) +
    labs(x = "position on chrX (Mb)", y = "escape (CAST fraction, bias-corrected)",
         title = "Per-gene chrX escape along the chromosome",
         subtitle = sprintf("reads pooled over cell types; >= %d informative reads per gene; dashed = %.1f%% pooled prior",
                            MIN_READS, 100 * PRIOR_ESCAPE)) +
    theme_minimal()
)

# 2. escape against the evidence behind it - the funnel IS the noise
print(
  ggplot(gx_pool, aes(total, escape_cast)) +
    geom_point(alpha = 0.3, size = 0.8) + scale_x_log10() +
    geom_hline(yintercept = PRIOR_ESCAPE, linetype = "dashed", colour = "grey40") +
    facet_wrap(~sample) +
    labs(x = "informative reads per gene (log)", y = "escape (CAST fraction)",
         title = "Per-gene escape against read depth",
         subtitle = "a gene near 0 or 1 at low depth is undersampled, not necessarily silent") +
    theme_minimal()
)

# 3. the autosomal null, on the same axis
print(
  ggplot(gene_pool, aes(escape_cast, colour = region)) +
    geom_density(linewidth = 0.9) + facet_wrap(~sample) +
    geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey60") +
    labs(x = "CAST fraction (bias-corrected)", y = "density",
         title = "chrX against the autosomal null",
         subtitle = "autosomes centre on 0.5 after correction; chrX should sit far below") +
    theme_minimal()
)

# 4. the named escape genes, per cell type
esc <- gx %>% filter(gene %in% ESCAPE_GENES)
if (nrow(esc)) {
  print(
    ggplot(esc, aes(reorder(gene, escape_cast), escape_cast, colour = sample)) +
      geom_errorbar(aes(ymin = escape_lo, ymax = escape_hi), width = 0.25,
                    position = position_dodge(0.5)) +
      geom_point(size = 2.4, position = position_dodge(0.5)) +
      coord_flip() + facet_wrap(~celltype) +
      scale_colour_manual(values = setNames(OKABE_ITO[1:2], SAMPLES)) +
      labs(x = NULL, y = "escape (CAST fraction, bias-corrected)",
           title = "Core escape genes, per cell type",
           subtitle = "bars are per-gene binomial 95% CI; n=1 animal per age, so the age contrast is descriptive") +
      theme_minimal()
  )
}
dev.off()
say("wrote gene_level_escape.pdf")

# ---------------------------------------------------------------------------
# distal enrichment, reusing the OCM statistics rather than reimplementing them
# ---------------------------------------------------------------------------
# chrX_distal_escape_enrichment.R already carries the Fisher test, the circular
# permutation that keeps escape clusters intact, the cluster-level binomial and
# the window sweep, all calibrated against Hoelzl et al. 2025. Sourcing it keeps
# one implementation of those tests in the project.
if (file.exists(COMMON) && length(SAMPLES) == 2) {
  Sys.setenv(AP2_ROOT = ROOT, GROUPS = paste(SAMPLES, collapse = ","), FIG_OUT = OUT)
  source(COMMON)
  pbg <- gx %>%
    transmute(sample = factor(sample, SAMPLES), celltype, gene, start,
              A1, A2, total, ar = ar_b6_corr,
              distal = is_distal(start),
              region = factor(ifelse(is_distal(start), DISTAL_LAB, "internal"),
                              c(DISTAL_LAB, "internal")))
  # distal_tests() recomputes the ratio from raw A1/A2 and applies ESCAPE_AR to
  # it, so the bias correction above does NOT reach the escape calls there. That
  # is fine for this test and wrong to patch around: the correction is monotone
  # and global, so using it would be exactly equivalent to moving the cutoff, and
  # a distal-vs-internal contrast within one chromosome is unchanged by that.
  # Only the ABSOLUTE escape counts differ, and those come from the table above.
  res <- distal_tests(pbg, groups = SAMPLES)
  write_tsv(res$fisher,      file.path(OUT, "distal_fisher.tsv"))
  write_tsv(res$permutation, file.path(OUT, "distal_permutation.tsv"))
  write_tsv(res$window_sweep, file.path(OUT, "distal_window_sweep.tsv"))
  say("")
  say("%s", res$caption)
} else {
  say("")
  say("skipping distal enrichment: need both samples and %s", COMMON)
}

say("")
say("output under %s", OUT)
