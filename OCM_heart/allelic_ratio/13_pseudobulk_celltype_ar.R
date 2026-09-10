# ---------------------------------------------------------------------------
# 13 - Per-cell-type PSEUDOBULK allelic ratios, all genes, all chromosomes,
#      and the positional question: are genes biallelic at chromosome ends?
#
# Run from the OCM/ directory:
#   Rscript $POSTDOC_ROOT/OCM_heart/allelic_ratio/13_pseudobulk_celltype_ar.R
#
# Reads  Allelome.PRO2_pseudobulk_celltype/<cohort>/<jobdir>/locus_table.txt
#        written by slurm/pseudobulk_celltype_allelome.slurm, one pseudobulk
#        per (cohort, cell type), doublets already excluded upstream by
#        OCM_heart/pseudobulk_celltype_ids.R.
#
# This is a different measurement from 03_all_genes.R, not a rerun of it:
#   03  one Allelome.PRO2 run per NUCLEUS, chrX annotation only. A single
#       nucleus gives a handful of reads per gene, so a per-gene ratio there
#       is mostly 0 or 1 by construction.
#   13  one run per (cohort, cell type) with every nucleus of that type
#       pooled, gene bodies on all 21 chromosomes. Cell-to-cell heterogeneity
#       is gone, per-gene depth is one to three orders of magnitude higher,
#       and a ratio near 0.5 can actually be resolved. That trade is the whole
#       point: the positional question is about genes, not cells.
#
# COLUMN NAMING, per CLAUDE.md: Allelome.PRO2's own `allelic_ratio` column is
# A1_reads/total_reads, directional, 0-1. A1 is B6. It is renamed here to
# `ar_a1` and never called `allelic_ratio`, and the undirected 0.5-1 form is
# `ar_dom`. The two are not interchangeable and the repo has been bitten by
# sharing the name before.
#
# GENOTYPE, per CLAUDE.md: B6 mother x CAST father, Xist deleted on the B6 X.
# The CAST X is the inactive X in every nucleus. So on chrX, A1 (B6) is the
# ACTIVE allele and the CAST fraction 1 - ar_a1 is ESCAPE. On autosomes
# neither allele is privileged and ar_dom is the quantity of interest.
#
# THE n=1 CAVEAT, per the standing note on this dataset: one animal per
# condition. Nothing here compares cohorts inferentially. The positional
# contrast is computed WITHIN a cohort x cell type, comparing terminal genes
# against depth- and SNP-matched interior genes of the same pseudobulk;
# cohort differences are reported descriptively only.
#
# There is deliberately no regression model of ratio on position. A logistic
# fit of the thresholded biallelic call on distance-to-end was tried first and
# removed: it discards the continuous measurement in favour of a MONO_AR
# indicator, and its p-value assumes genes are independent draws when terminal
# genes are spatially clustered by definition -- the hypothesis cannot also be
# the null. A matched contrast with a shuffled-label reference answers the
# actual question without either problem.
# ---------------------------------------------------------------------------
source(file.path(Sys.getenv("POSTDOC_ROOT",
                            "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc"),
                 "OCM_heart/allelic_ratio/00_functions.R"))

library(tidyr)

ALLELOME_TREE <- Sys.getenv("PSEUDOBULK_TREE", "Allelome.PRO2_pseudobulk_celltype")
PB_DIR        <- Sys.getenv("PSEUDOBULK_DIR", "pseudobulk_celltype")
CHROM_SIZES   <- Sys.getenv("CHROM_SIZES",
  "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/chr_annotation_mm39.bed")
SNP_COUNT_BED <- file.path(PB_DIR, "gene_level_snp_count.bed")

COHORTS <- c("9w", "78w", "Sham", "TAC")
# the cohort labels the question was asked in
COHORT_ALIAS <- c("9w" = "adult", "78w" = "aged", "Sham" = "Sham", "TAC" = "TAC")

OUT_DIR <- file.path(RESULTS_ROOT, "pseudobulk_celltype")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Per-gene depth cutoff. MIN_GENE_READS (default 10) is the repo-wide per-gene
# threshold from 00_functions.R -- the same number 03 uses -- so the pseudobulk
# and per-cell per-gene tables are thresholded alike.
MIN_PB_CELLS <- as.integer(Sys.getenv("MIN_PB_CELLS", "20"))
stopifnot(!is.na(MIN_PB_CELLS), MIN_PB_CELLS >= 1)

AUTOSOME_SET <- paste0("chr", 1:19)

# What counts as "escaping" for the chrX summary. 0.10 sits just under the
# ~12.7% whole-chrX CAST fraction measured in this cross, so a gene has to
# carry at least about the chromosome-average escape to be counted.
ESCAPE_CUT <- as.numeric(Sys.getenv("ESCAPE_CUT", "0.10"))

message(sprintf("tree=%s  MIN_GENE_READS=%d  MIN_PB_CELLS=%d  MONO_AR=%.2f",
                ALLELOME_TREE, MIN_GENE_READS, MIN_PB_CELLS, MONO_AR))

# ---------------------------------------------------------------------------
# 1. Ingest
# ---------------------------------------------------------------------------
# The cell type is recovered from the <celltype>.log that Allelome.PRO2 writes
# inside each job directory, NOT by splitting the directory name on "_": the
# labels themselves contain underscores ("Pericytes_Smooth_muscle_cells"), so
# positional splitting -- which is what 03_all_genes.R does for barcodes --
# would silently mis-assign them. The directory-name regex is the fallback.
label_from_jobdir <- function(dir) {
  logs <- list.files(dir, pattern = "\\.log$", full.names = FALSE)
  if (length(logs) == 1L) return(sub("\\.log$", "", logs))
  base <- basename(dir)
  # <YYYY_MM_DD>_<sample>_<annotation.bed>_<minreads>
  lab <- sub("^[0-9]{4}_[0-9]{2}_[0-9]{2}_(.+)_[^_]*\\.bed_[0-9]+$", "\\1", base)
  if (identical(lab, base)) return(NA_character_)
  lab
}

read_pseudobulk <- function(cohort) {
  cohort_dir <- file.path(ALLELOME_TREE, cohort)
  if (!dir.exists(cohort_dir)) {
    warning("no directory for cohort ", cohort, ": ", cohort_dir)
    return(NULL)
  }
  jobdirs <- list.dirs(cohort_dir, recursive = FALSE, full.names = TRUE)
  jobdirs <- jobdirs[file.exists(file.path(jobdirs, "locus_table.txt"))]
  if (!length(jobdirs)) {
    warning("cohort ", cohort, ": no locus_table.txt under ", cohort_dir)
    return(NULL)
  }

  rows <- lapply(jobdirs, function(d) {
    lab <- label_from_jobdir(d)
    if (is.na(lab)) {
      warning("cannot recover a cell type label from ", d, " -- skipped")
      return(NULL)
    }
    tb <- tryCatch(read.delim(file.path(d, "locus_table.txt"), header = TRUE,
                              stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (is.null(tb) || !nrow(tb)) {
      warning("empty locus_table.txt in ", d)
      return(NULL)
    }
    need <- c("chr", "start", "end", "name",
              "A1_reads", "A2_reads", "total_reads", "allelic_ratio")
    if (!all(need %in% names(tb))) {
      warning(d, ": locus_table.txt is missing ",
              paste(setdiff(need, names(tb)), collapse = ", "), " -- skipped")
      return(NULL)
    }
    tb$cohort  <- cohort
    tb$label   <- lab
    tb$jobdir  <- basename(d)
    tb[, c("cohort", "label", "jobdir", need)]
  })

  bind_rows(rows)
}

gene_df <- bind_rows(lapply(COHORTS, read_pseudobulk))
if (!nrow(gene_df)) {
  stop("no pseudobulk locus tables found under ", ALLELOME_TREE,
       ".\n  Run slurm/pseudobulk_celltype_bams.slurm then ",
       "slurm/pseudobulk_celltype_allelome.slurm first.")
}

# One annotation was used for every run; if not, the ratios are not comparable.
# The annotation is what is left of the job directory name after the date
# prefix, the cell type label and the trailing minreads are removed. Done this
# way, not with a regex on "_", because both the label and the annotation
# basename contain underscores -- the pattern in 00_functions.R's
# allelome_annotations() cannot span them and returns just "level.bed".
annots <- unique(mapply(function(jd, lab) {
  a <- sub(paste0("^[0-9]{4}_[0-9]{2}_[0-9]{2}_", lab, "_"), "", jd)
  sub("_[0-9]+$", "", a)
}, gene_df$jobdir, gene_df$label, USE.NAMES = FALSE))
if (length(annots) != 1L) {
  stop("pseudobulks were scored against more than one annotation: ",
       paste(annots, collapse = ", "),
       "\n  Ratios from different beds are not comparable. Re-run the odd one out.")
}
message("annotation: ", annots)

gene_df <- gene_df %>%
  dplyr::rename(gene = name, ar_a1 = allelic_ratio) %>%
  dplyr::mutate(
    cohort      = factor(cohort, levels = COHORTS),
    condition   = unname(COHORT_ALIAS[as.character(cohort)]),
    total_reads = as.integer(total_reads),
    A1_reads    = as.integer(A1_reads),
    A2_reads    = as.integer(A2_reads),
    # recomputed rather than trusted: locus_table's allelic_ratio is rounded
    # to 2 dp, which is coarse enough to move a gene across MONO_AR
    ar_a1  = A1_reads / total_reads,
    ar_dom = pmax(A1_reads, A2_reads) / total_reads
  )

stopifnot(all(abs(gene_df$A1_reads + gene_df$A2_reads - gene_df$total_reads) < 1e-9))

# ---------------------------------------------------------------------------
# 2. Cell type labels, nuclei counts, SNP density
# ---------------------------------------------------------------------------
lab_file <- file.path(PB_DIR, "celltype_labels.txt")
if (file.exists(lab_file)) {
  labs <- read.delim(lab_file, stringsAsFactors = FALSE)
  gene_df$celltype <- labs$celltype[match(gene_df$label, labs$label)]
  miss <- unique(gene_df$label[is.na(gene_df$celltype)])
  if (length(miss)) {
    warning("labels absent from ", lab_file, ": ", paste(miss, collapse = ", "))
    gene_df$celltype[is.na(gene_df$celltype)] <-
      gene_df$label[is.na(gene_df$celltype)]
  }
} else {
  warning(lab_file, " not found -- using the sanitised labels as cell type names")
  gene_df$celltype <- gene_df$label
}

cnt_file <- file.path(PB_DIR, "pseudobulk_cell_counts.txt")
if (file.exists(cnt_file)) {
  cnts <- read.delim(cnt_file, stringsAsFactors = FALSE)
  key <- paste(gene_df$cohort, gene_df$label, sep = "|")
  gene_df$n_cells <- cnts$n_cells[match(key, paste(cnts$cohort, cnts$label, sep = "|"))]
  if (anyNA(gene_df$n_cells)) {
    warning(sum(is.na(gene_df$n_cells)), " rows have no nuclei count in ", cnt_file)
  }
} else {
  warning(cnt_file, " not found -- n_cells will be NA and MIN_PB_CELLS cannot be applied")
  gene_df$n_cells <- NA_integer_
}

# Informative SNPs per gene. A subtelomeric gene with three usable SNPs and a
# mid-chromosome gene with three hundred are not equally well measured, and
# B6/CAST divergence is not uniform along a chromosome -- so this has to be a
# covariate in any positional model, not an afterthought.
if (file.exists(SNP_COUNT_BED)) {
  snpc <- read.delim(SNP_COUNT_BED, header = FALSE, stringsAsFactors = FALSE)
  names(snpc)[c(1, 2, 3, 4, ncol(snpc))] <-
    c("chr", "start", "end", "gene", "n_snps")
  snpc <- snpc[, c("chr", "start", "end", "gene", "n_snps")]
  k1 <- paste(gene_df$chr, gene_df$start, gene_df$end, gene_df$gene, sep = "|")
  k2 <- paste(snpc$chr, snpc$start, snpc$end, snpc$gene, sep = "|")
  gene_df$n_snps <- snpc$n_snps[match(k1, k2)]
  message(sprintf("SNPs per gene joined for %.1f%% of rows",
                  100 * mean(!is.na(gene_df$n_snps))))
} else {
  warning(SNP_COUNT_BED, " not found -- n_snps unavailable, ",
          "positional models will omit the SNP-density covariate")
  gene_df$n_snps <- NA_integer_
}

# ---------------------------------------------------------------------------
# 3. Chromosome geometry
# ---------------------------------------------------------------------------
if (!file.exists(CHROM_SIZES)) {
  stop("chromosome sizes not found: ", CHROM_SIZES)
}
sizes <- read.delim(CHROM_SIZES, header = FALSE, stringsAsFactors = FALSE)
sizes <- data.frame(chr = as.character(sizes[[1]]),
                    chr_len = as.numeric(sizes[[3]]),
                    stringsAsFactors = FALSE)
sizes <- sizes[!duplicated(sizes$chr), ]

gene_df$chr_len <- sizes$chr_len[match(gene_df$chr, sizes$chr)]
if (anyNA(gene_df$chr_len)) {
  bad <- unique(gene_df$chr[is.na(gene_df$chr_len)])
  warning("no length for ", paste(bad, collapse = ", "), " -- those rows dropped")
  gene_df <- gene_df[!is.na(gene_df$chr_len), ]
}

# Mouse chromosomes are TELOCENTRIC: the centromere sits at the very start, so
# the two "ends" are not the same object. Low coordinates are pericentromeric
# (satellite-rich, poorly mappable, and the annotation itself starts ~3 Mb in);
# high coordinates are the true distal telomere. Pooling them into one
# "distance to nearest end" hides that, so all three are carried and the
# proximal/distal split is reported separately.
gene_df <- gene_df %>%
  dplyr::mutate(
    mid          = (start + end) / 2,
    dist_prox    = mid,                       # to the centromeric start
    dist_dist    = chr_len - mid,             # to the distal telomere
    dist_nearest = pmin(mid, chr_len - mid),
    nearest_end  = ifelse(mid <= chr_len - mid, "proximal", "distal"),
    rel_pos      = mid / chr_len,
    chr_class    = ifelse(chr %in% AUTOSOME_SET, "autosome",
                   ifelse(chr == "chrX", "chrX", "other")),
    # escape is only defined on chrX, and only because the CAST X is the Xi
    escape_frac  = ifelse(chr == "chrX", 1 - ar_a1, NA_real_),
    biallelic    = ar_dom < MONO_AR
  )

DIST_BREAKS <- c(0, 1, 2, 5, 10, 20, Inf) * 1e6
DIST_LABELS <- c("0-1 Mb", "1-2 Mb", "2-5 Mb", "5-10 Mb", "10-20 Mb", ">20 Mb")
gene_df$dist_bin <- cut(gene_df$dist_nearest, breaks = DIST_BREAKS,
                        labels = DIST_LABELS, include.lowest = TRUE)

# ---------------------------------------------------------------------------
# 4. The deliverable tables
# ---------------------------------------------------------------------------
KEEP <- c("cohort", "condition", "celltype", "label", "n_cells",
          "chr", "start", "end", "gene", "n_snps",
          "A1_reads", "A2_reads", "total_reads", "ar_a1", "ar_dom",
          "escape_frac", "biallelic",
          "chr_len", "mid", "rel_pos",
          "dist_prox", "dist_dist", "dist_nearest", "nearest_end",
          "dist_bin", "chr_class")
all_rows <- gene_df[, KEEP]

write.table(all_rows, file.path(OUT_DIR, "pseudobulk_celltype_ar_all_genes_raw.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Filtered: the table every figure and model below uses.
ar <- all_rows %>%
  dplyr::filter(total_reads >= MIN_GENE_READS,
                is.na(n_cells) | n_cells >= MIN_PB_CELLS)

if (!nrow(ar)) {
  stop("nothing survives total_reads >= ", MIN_GENE_READS,
       " and n_cells >= ", MIN_PB_CELLS)
}

write.table(ar, file.path(OUT_DIR, sprintf(
  "pseudobulk_celltype_ar_all_genes_min%dreads.txt", MIN_GENE_READS)),
  sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\n%d gene x pseudobulk rows total; %d pass total_reads >= %d and n_cells >= %d\n",
            nrow(all_rows), nrow(ar), MIN_GENE_READS, MIN_PB_CELLS))

# Wide gene x (condition, cell type) matrix of the directional ratio, the
# convenient form for eyeballing one gene across the design.
wide <- ar %>%
  dplyr::mutate(group = paste(condition, label, sep = "__")) %>%
  dplyr::select(chr, start, end, gene, group, ar_a1) %>%
  tidyr::pivot_wider(names_from = group, values_from = ar_a1,
                     values_fn = mean) %>%
  dplyr::arrange(chr, start)
write.table(wide, file.path(OUT_DIR, "pseudobulk_celltype_ar_a1_wide.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Per-pseudobulk QC. Read this before any figure: a pseudobulk with a few
# hundred informative reads cannot support a per-gene ratio at all.
qc <- all_rows %>%
  dplyr::group_by(condition, cohort, celltype, label, n_cells) %>%
  dplyr::summarise(
    genes_any        = dplyr::n(),
    genes_pass       = sum(total_reads >= MIN_GENE_READS),
    genes_pass_chrX  = sum(total_reads >= MIN_GENE_READS & chr == "chrX"),
    snp_reads_total  = sum(total_reads),
    median_gene_reads = stats::median(total_reads),
    .groups = "drop"
  ) %>%
  dplyr::arrange(condition, dplyr::desc(genes_pass))
write.table(qc, file.path(OUT_DIR, "pseudobulk_celltype_qc.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n--- per-pseudobulk QC ---\n"); print(as.data.frame(qc), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Direction check: the genotype has to be visible in the numbers
# ---------------------------------------------------------------------------
# With Xist deleted on the B6 X, the CAST X is the Xi in every nucleus, so
# pooled chrX must come out B6-dominant (ar_a1 well above 0.5) while the
# autosomes sit at ~0.5. If that is not what the table says, the allele
# assignment is wrong and nothing downstream means anything.
dir_check <- ar %>%
  dplyr::group_by(condition, label, chr_class) %>%
  dplyr::summarise(mean_ar_a1 = mean(ar_a1),
                   median_ar_a1 = stats::median(ar_a1),
                   n_genes = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = chr_class,
                     values_from = c(mean_ar_a1, median_ar_a1, n_genes))
write.table(dir_check, file.path(OUT_DIR, "pseudobulk_direction_check.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n--- direction check: chrX should be B6-dominant, autosomes ~0.5 ---\n")
print(as.data.frame(dir_check), row.names = FALSE)

x_mean <- mean(ar$ar_a1[ar$chr_class == "chrX"])
a_mean <- mean(ar$ar_a1[ar$chr_class == "autosome"])
cat(sprintf("\npooled mean ar_a1 (B6 fraction): chrX %.3f, autosomes %.3f\n",
            x_mean, a_mean))
if (!is.na(x_mean) && x_mean < a_mean) {
  warning("chrX is LESS B6-dominant than the autosomes. Under this genotype ",
          "that should be impossible -- check the allele assignment (A1 = B6) ",
          "before interpreting anything below.")
}

# ---------------------------------------------------------------------------
# 6. The positional question
# ---------------------------------------------------------------------------
# Read this before the numbers. "Biallelic at chromosome ends" is two
# different questions on the two chromosome classes:
#
#   chrX      Every gene is expected monoallelic (active B6 X) except escapers.
#             So "biallelic" here means ESCAPE, and the question is whether
#             escape concentrates towards the chromosome ends.
#
#   autosome  Every gene is expected biallelic already, imprinted loci aside.
#             So a positional trend on the autosomes is not biology, it is the
#             technical baseline -- mappability, SNP density and depth all vary
#             towards the ends. It is the control the chrX trend has to beat.
#
# And the confounder that has to be held: apparent monoallelism is manufactured
# by low depth. At total_reads = 10 the smallest non-zero minor fraction is
# 0.1, so a shallow gene is pushed towards ar_dom = 1 whether or not it is
# monoallelic. Depth is therefore stratified in the summary AND carried as a
# covariate in the model. Any trend that survives both is worth a second look.

summarise_bins <- function(df, ...) {
  df %>%
    dplyr::group_by(...) %>%
    dplyr::summarise(
      n_genes       = dplyr::n(),
      frac_bi       = mean(biallelic),
      mean_ar_dom   = mean(ar_dom),
      median_ar_dom = stats::median(ar_dom),
      mean_escape   = mean(escape_frac),
      median_reads  = stats::median(total_reads),
      median_snps   = stats::median(n_snps, na.rm = TRUE),
      .groups = "drop"
    )
}

bin_summary <- summarise_bins(ar, condition, label, celltype, chr_class, dist_bin)
write.table(bin_summary, file.path(OUT_DIR, "pseudobulk_dist_bin_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# proximal (pericentromeric) vs distal (telomeric) ends kept apart
end_summary <- ar %>%
  dplyr::filter(dist_nearest <= 10e6) %>%
  summarise_bins(condition, label, chr_class, nearest_end, dist_bin)
write.table(end_summary, file.path(OUT_DIR, "pseudobulk_end_type_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Depth-stratified version of the same summary: bins compared only against
# genes of comparable depth, then averaged back over strata with equal weight.
# ntile() rather than cut(quantile()) because it is row-aligned inside
# group_by/mutate; building the strata separately and pasting the vector back
# would order them by group and the rows by their original position.
ar <- ar %>%
  dplyr::group_by(condition, label) %>%
  dplyr::mutate(depth_stratum = dplyr::ntile(total_reads, 5)) %>%
  dplyr::ungroup() %>%
  as.data.frame()

depth_matched <- ar %>%
  dplyr::group_by(condition, label, chr_class, depth_stratum, dist_bin) %>%
  dplyr::summarise(n_genes = dplyr::n(), frac_bi = mean(biallelic),
                   mean_ar_dom = mean(ar_dom), .groups = "drop") %>%
  dplyr::group_by(condition, label, chr_class, dist_bin) %>%
  dplyr::summarise(
    n_strata          = dplyr::n(),
    frac_bi_matched   = mean(frac_bi),
    ar_dom_matched    = mean(mean_ar_dom),
    n_genes           = sum(n_genes),
    .groups = "drop"
  )
write.table(depth_matched, file.path(OUT_DIR, "pseudobulk_dist_bin_depth_matched.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n--- fraction biallelic (ar_dom < ", MONO_AR,
    ") by distance to the nearest chromosome end ---\n", sep = "")
print(as.data.frame(
  bin_summary %>%
    dplyr::group_by(chr_class, dist_bin) %>%
    # the weighted means MUST come before n_genes is overwritten: dplyr
    # evaluates summarise() arguments in order and each one shadows the
    # column of the same name for the arguments after it
    dplyr::summarise(frac_bi      = stats::weighted.mean(frac_bi, n_genes),
                     mean_ar_dom  = stats::weighted.mean(mean_ar_dom, n_genes),
                     median_reads = stats::median(median_reads),
                     n_genes      = sum(n_genes),
                     .groups = "drop")
), row.names = FALSE)

# ---- terminal vs matched-interior genes ------------------------------------
# The question is whether genes near a chromosome end look different from
# genes that are not. That needs a reference set, and the only honest one is
# genes matched on the two things that move an allelic ratio on their own:
# per-gene read depth and the number of informative SNPs.
#
# So: stratify every gene in a pseudobulk by depth quintile x SNP-count
# tertile, and inside each stratum compare terminal genes against interior
# ones. The reported statistic is the stratum-weighted difference in the
# CONTINUOUS ratio -- escape fraction on chrX, ar_dom on the autosomes --
# because that is the measurement. The biallelic call at MONO_AR is carried
# alongside as a secondary readout, since it is the repo-wide boundary, but it
# is a thresholded version of the same number and is not the primary result.
#
# The reference distribution comes from shuffling the terminal/interior label
# WITHIN strata, which is exactly the null being claimed ("position is
# unrelated to the ratio, once depth and SNP density are held") and needs no
# distributional assumption.
#
# Its one real limitation, stated once: neighbouring genes share mappability
# and reads, and shuffling within strata breaks that spatial correlation, so
# the permutation p is somewhat optimistic. It is a reference, not a
# certificate. The check that actually matters is whether the SAME SIGN shows
# up in independent cell types -- reported in the consistency table below.
TERM_MB <- as.numeric(Sys.getenv("TERMINAL_MB", "5"))
N_PERM  <- as.integer(Sys.getenv("N_PERM", "2000"))
stopifnot(!is.na(TERM_MB), TERM_MB > 0, !is.na(N_PERM), N_PERM >= 99)

# The window goes in the filename of everything that depends on it. chrX and
# the autosomes do not want the same window -- the autosomes have 360 genes in
# the terminal 1 Mb alone, while chrX has 12, so chrX has to be run wider --
# and without this the second run would silently overwrite the first. The
# per-gene tables are window-independent and keep their plain names.
TERM_TAG <- paste0("_term", sub("\\.0$", "", format(TERM_MB, trim = TRUE)), "Mb")

# The primary metric differs by chromosome class because the biology does:
# on chrX every gene is monoallelic bar escapers, so escape fraction IS the
# quantity; on the autosomes nothing is privileged, so ar_dom is.
ar$pos_metric     <- ifelse(ar$chr_class == "chrX", ar$escape_frac, ar$ar_dom)
ar$biallelic_num  <- as.numeric(ar$biallelic)

PERM_METRICS <- c(pos_metric = "ratio", biallelic_num = "frac_biallelic")

# Stratum-weighted terminal-minus-interior difference, with a within-stratum
# permutation null. Vectorised over metrics and over permutations: for a
# stratum of n genes with k terminal, a permuted statistic needs only the SUM
# of a random k-subset, so one sample.int + one colSums per stratum per draw.
matched_contrast <- function(df, end, n_perm = N_PERM, seed = 42) {
  win <- TERM_MB * 1e6
  near_either <- df$dist_nearest <= win
  terminal <- switch(end,
    nearest  = near_either,
    proximal = df$dist_prox <= win,
    distal   = df$dist_dist <= win,
    stop("unknown end: ", end))

  # the interior is genes far from BOTH ends, so the opposite end never
  # contaminates the reference set in the proximal/distal comparisons
  interior <- !near_either
  keep <- terminal | interior
  df <- df[keep, , drop = FALSE]
  terminal <- terminal[keep]
  if (sum(terminal) < 10 || sum(!terminal) < 20) return(NULL)

  metrics <- names(PERM_METRICS)
  ok_rows <- stats::complete.cases(df[, c(metrics, "total_reads")])
  df <- df[ok_rows, , drop = FALSE]; terminal <- terminal[ok_rows]
  if (sum(terminal) < 10 || sum(!terminal) < 20) return(NULL)

  d_str <- dplyr::ntile(df$total_reads, 5)
  s_str <- if (all(!is.na(df$n_snps))) dplyr::ntile(df$n_snps, 3) else rep(1L, nrow(df))
  stratum <- paste(d_str, s_str, sep = ".")

  V <- as.matrix(df[, metrics, drop = FALSE])
  storage.mode(V) <- "double"

  set.seed(seed)
  parts <- list()
  for (st in unique(stratum)) {
    idx <- which(stratum == st)
    tt  <- terminal[idx]
    n <- length(idx); k <- sum(tt)
    if (k < 1 || (n - k) < 2) next          # a difference is undefined here
    Vs <- V[idx, , drop = FALSE]
    Tot <- colSums(Vs)
    S_obs <- colSums(Vs[tt, , drop = FALSE])
    obs_s <- S_obs / k - (Tot - S_obs) / (n - k)
    null_s <- vapply(seq_len(n_perm), function(i) {
      S <- colSums(Vs[sample.int(n, k), , drop = FALSE])
      S / k - (Tot - S) / (n - k)
    }, numeric(length(metrics)))
    if (length(metrics) == 1L) null_s <- matrix(null_s, nrow = 1L)
    parts[[st]] <- list(w = k, obs = obs_s, null = null_s,
                        mt = S_obs / k, mi = (Tot - S_obs) / (n - k))
  }
  if (!length(parts)) return(NULL)

  w <- vapply(parts, function(x) x$w, 0)
  obs  <- Reduce(`+`, Map(function(x) x$obs  * x$w, parts)) / sum(w)
  null <- Reduce(`+`, Map(function(x) x$null * x$w, parts)) / sum(w)

  p_emp <- (1 + rowSums(abs(null) >= abs(obs))) / (1 + n_perm)
  nm <- rowMeans(null); nsd <- apply(null, 1, stats::sd)
  # a degenerate null (every gene in every usable stratum has the same value,
  # e.g. all biallelic) has zero spread, so z is undefined rather than huge
  z <- ifelse(nsd > 0, (obs - nm) / nsd, NA_real_)

  # stratum-weighted group means, per metric -- same weights as the difference
  mt <- Reduce(`+`, Map(function(x) x$mt * x$w, parts)) / sum(w)
  mi <- Reduce(`+`, Map(function(x) x$mi * x$w, parts)) / sum(w)

  data.frame(
    end            = end,
    metric         = unname(PERM_METRICS[metrics]),
    n_terminal     = sum(terminal),
    n_interior     = sum(!terminal),
    n_strata       = length(parts),
    mean_terminal  = unname(mt),
    mean_interior  = unname(mi),
    diff_matched   = unname(obs),
    null_sd        = unname(nsd),
    z              = unname(z),
    p_perm         = unname(p_emp),
    direction      = ifelse(obs > 0, "higher at the ends", "lower at the ends"),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

contrast_input <- ar[ar$chr_class %in% c("autosome", "chrX"), ]
contrast_keys  <- unique(contrast_input[, c("condition", "label", "celltype", "chr_class")])

end_contrasts <- bind_rows(lapply(seq_len(nrow(contrast_keys)), function(i) {
  k <- contrast_keys[i, ]
  sub <- contrast_input[contrast_input$condition == k$condition &
                        contrast_input$label     == k$label &
                        contrast_input$chr_class == k$chr_class, ]
  bind_rows(lapply(c("nearest", "proximal", "distal"), function(e) {
    r <- matched_contrast(sub, e)
    if (is.null(r)) return(NULL)
    cbind(k, r, row.names = NULL)
  }))
}))

# Which contrasts were possible at all. A 5 Mb window at the end of a 169 Mb
# chromosome holds few genes, and fewer still pass the depth cutoff, so
# matched_contrast() declines any pseudobulk with under 10 terminal or 20
# interior genes. That has to be visible: an absent row is a pseudobulk that
# could not be measured, not one with no effect.
coverage <- bind_rows(lapply(seq_len(nrow(contrast_keys)), function(i) {
  k <- contrast_keys[i, ]
  sub <- contrast_input[contrast_input$condition == k$condition &
                        contrast_input$label     == k$label &
                        contrast_input$chr_class == k$chr_class, ]
  win <- TERM_MB * 1e6
  near <- sub$dist_nearest <= win
  bind_rows(lapply(c("nearest", "proximal", "distal"), function(e) {
    term <- switch(e, nearest = near,
                      proximal = sub$dist_prox <= win,
                      distal   = sub$dist_dist <= win)
    cbind(k, data.frame(end = e,
                        n_terminal = sum(term),
                        n_interior = sum(!near),
                        tested = sum(term) >= 10 & sum(!near) >= 20,
                        row.names = NULL))
  }))
}))
write.table(coverage, file.path(OUT_DIR, paste0("pseudobulk_end_contrast_coverage", TERM_TAG, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)
n_skip <- sum(!coverage$tested)
if (n_skip) {
  cat(sprintf("\nNOTE: %d of %d (pseudobulk x chromosome class x end) contrasts were not testable\n",
              n_skip, nrow(coverage)))
  cat(sprintf("      with a %g Mb terminal window (need >=10 terminal and >=20 interior genes).\n", TERM_MB))
  cat("      Per chromosome class:\n")
  print(as.data.frame(
    coverage %>% dplyr::group_by(chr_class, end) %>%
      # n() - sum(tested), not sum(!tested): summarise() evaluates its
      # arguments in order and `tested = sum(tested)` would shadow the column
      # for every argument after it, silently reporting zero skips
      # `skipped` BEFORE `tested`: summarise() evaluates its arguments in
      # order, so once `tested = sum(tested)` has run the column is a scalar
      # and any later reference to it silently reports zero skips
      dplyr::summarise(n_contrasts = dplyr::n(),
                       skipped  = sum(!tested),
                       tested   = sum(tested),
                       median_n_terminal = stats::median(n_terminal),
                       .groups = "drop")), row.names = FALSE)
  cat("      Widen with TERMINAL_MB=10 if chrX is mostly skipped; see coverage table.\n")
}

if (nrow(end_contrasts)) {
  end_contrasts <- end_contrasts %>%
    dplyr::arrange(chr_class, metric, end, dplyr::desc(abs(z)))
  write.table(end_contrasts, file.path(OUT_DIR, paste0("pseudobulk_end_vs_matched_interior", TERM_TAG, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("\n--- genes within %g Mb of a chromosome end vs depth- and SNP-matched interior genes ---\n",
              TERM_MB))
  cat("    diff_matched > 0 = higher at the ends; p_perm from label shuffling within strata\n")
  print(as.data.frame(
    end_contrasts %>%
      dplyr::filter(metric == "ratio") %>%
      dplyr::select(condition, celltype, chr_class, end, n_terminal,
                    mean_terminal, mean_interior, diff_matched, z, p_perm)
  ), row.names = FALSE, digits = 3)

  # The replication check. One cell type showing a trend is one measurement;
  # the same sign in independent cell types of the same animal is the thing
  # worth acting on. No FDR here on purpose -- these fits share genes and
  # positions, so they are not independent tests and a q-value would be
  # arithmetic without meaning.
  consistency <- end_contrasts %>%
    dplyr::group_by(condition, chr_class, end, metric) %>%
    dplyr::summarise(
      n_celltypes   = dplyr::n(),
      n_positive    = sum(diff_matched > 0),
      n_p_below_05  = sum(p_perm < 0.05),
      median_diff   = stats::median(diff_matched),
      .groups = "drop"
    ) %>%
    dplyr::mutate(sign_agreement = pmax(n_positive, n_celltypes - n_positive) / n_celltypes) %>%
    dplyr::arrange(metric, chr_class, end, dplyr::desc(sign_agreement))
  write.table(consistency, file.path(OUT_DIR, paste0("pseudobulk_end_sign_consistency", TERM_TAG, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n--- sign consistency across cell types (the replication check) ---\n")
  print(as.data.frame(consistency %>% dplyr::filter(metric == "ratio")),
        row.names = FALSE, digits = 3)
} else {
  warning("no terminal-vs-interior contrast could be computed -- ",
          "too few genes inside ", TERM_MB, " Mb of an end")
  end_contrasts <- data.frame()
}

# ---- chrX escape at the ends, on its own terms -----------------------------
# On chrX the interesting quantity is the escape fraction itself, not a
# biallelic/monoallelic call, so it gets a beta-scale summary as well.
if (any(ar$chr_class == "chrX")) {
  x_bins <- ar %>%
    dplyr::filter(chr_class == "chrX") %>%
    dplyr::group_by(condition, label, dist_bin) %>%
    dplyr::summarise(n_genes = dplyr::n(),
                     mean_escape = mean(escape_frac),
                     median_escape = stats::median(escape_frac),
                     frac_escaping = mean(escape_frac > ESCAPE_CUT),
                     median_reads = stats::median(total_reads),
                     .groups = "drop")
  write.table(x_bins, file.path(OUT_DIR, "pseudobulk_chrX_escape_by_dist.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n--- chrX escape fraction (CAST / total) by distance to end, pooled ---\n")
  print(as.data.frame(
    x_bins %>% dplyr::group_by(dist_bin) %>%
      dplyr::summarise(mean_escape = stats::weighted.mean(mean_escape, n_genes),
                       n = sum(n_genes),
                       .groups = "drop")
  ), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# 7. Figures
# ---------------------------------------------------------------------------
pdf(file.path(OUT_DIR, paste0("pseudobulk_celltype_chromosome_ends", TERM_TAG, ".pdf")),
    width = 12, height = 7)

for (cc in c("autosome", "chrX")) {
  d <- bin_summary %>% dplyr::filter(chr_class == cc, !is.na(dist_bin))
  if (!nrow(d)) next
  print(
    ggplot(d, aes(x = dist_bin, y = frac_bi, colour = celltype, group = celltype)) +
      geom_line(alpha = 0.7) + geom_point(size = 1.4) +
      facet_wrap(~ condition, nrow = 1) +
      labs(title = sprintf("%s: fraction biallelic (ar_dom < %.2f) vs distance to nearest chromosome end", cc, MONO_AR),
           subtitle = sprintf("pseudobulk per cell type, doublets excluded, genes with >= %d SNP-overlapping reads", MIN_GENE_READS),
           x = "distance to nearest chromosome end", y = "fraction biallelic") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                         legend.position = "right")
  )
  print(
    ggplot(d, aes(x = dist_bin, y = median_reads, colour = celltype, group = celltype)) +
      geom_line(alpha = 0.7) + geom_point(size = 1.4) +
      facet_wrap(~ condition, nrow = 1) + scale_y_log10() +
      labs(title = sprintf("%s: the confounder -- median per-gene depth by the same bins", cc),
           subtitle = "a depth gradient towards the ends would produce a biallelic gradient on its own",
           x = "distance to nearest chromosome end", y = "median reads per gene (log10)") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
}

dm <- depth_matched %>% dplyr::filter(!is.na(dist_bin),
                                      chr_class %in% c("autosome", "chrX"))
if (nrow(dm)) {
  print(
    ggplot(dm, aes(x = dist_bin, y = frac_bi_matched, colour = chr_class,
                   group = interaction(chr_class, label))) +
      geom_line(alpha = 0.5) +
      facet_wrap(~ condition, nrow = 1) +
      labs(title = "Depth-matched fraction biallelic vs distance to nearest end",
           subtitle = "genes compared only within quintiles of per-gene depth, then averaged over strata; one line per cell type",
           x = "distance to nearest chromosome end", y = "depth-matched fraction biallelic") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
}

# proximal (pericentromeric) vs distal (telomeric): telocentric mouse
# chromosomes make these two different places, not two instances of "an end"
es <- end_summary %>% dplyr::filter(chr_class %in% c("autosome", "chrX"),
                                    !is.na(dist_bin))
if (nrow(es)) {
  print(
    ggplot(es, aes(x = dist_bin, y = frac_bi, colour = nearest_end,
                   group = interaction(nearest_end, label))) +
      geom_line(alpha = 0.5) +
      facet_grid(chr_class ~ condition) +
      labs(title = "Pericentromeric (proximal) vs telomeric (distal) ends kept apart",
           subtitle = "mouse chromosomes are telocentric: low coordinates are centromere, not telomere",
           x = "distance to that end", y = "fraction biallelic") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
}

# along-chromosome profile, cell types pooled within a condition
prof <- ar %>%
  dplyr::filter(chr_class %in% c("autosome", "chrX")) %>%
  dplyr::mutate(pos_bin = cut(rel_pos, breaks = seq(0, 1, 0.05),
                              include.lowest = TRUE)) %>%
  dplyr::group_by(condition, chr_class, pos_bin) %>%
  dplyr::summarise(n = dplyr::n(), frac_bi = mean(biallelic),
                   mean_ar_dom = mean(ar_dom), .groups = "drop") %>%
  dplyr::mutate(pos_mid = (as.integer(pos_bin) - 0.5) * 0.05)
if (nrow(prof)) {
  print(
    ggplot(prof, aes(x = pos_mid, y = frac_bi, colour = condition)) +
      geom_line() + geom_point(size = 1) +
      facet_wrap(~ chr_class, scales = "free_y") +
      labs(title = "Fraction biallelic along the chromosome (0 = centromeric start, 1 = distal telomere)",
           subtitle = "cell types pooled; 5% relative-position bins",
           x = "relative position along chromosome", y = "fraction biallelic") +
      theme_bw()
  )
}

if (nrow(end_contrasts)) {
  ec <- end_contrasts %>% dplyr::filter(metric == "ratio")
  print(
    ggplot(ec, aes(x = diff_matched, y = reorder(celltype, diff_matched),
                   colour = condition)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
      # +-2 permutation SDs of the within-strata null, for scale
      geom_linerange(aes(xmin = diff_matched - 2 * null_sd,
                         xmax = diff_matched + 2 * null_sd),
                     alpha = 0.5) +
      geom_point(size = 2) +
      facet_grid(chr_class ~ end, scales = "free_x") +
      labs(title = sprintf("Genes within %g Mb of a chromosome end vs matched interior genes",
                           TERM_MB),
           subtitle = "difference in escape fraction (chrX) / ar_dom (autosomes), matched on depth and SNP count; bars are +-2 SD of the permutation null",
           x = "terminal minus matched interior", y = NULL) +
      theme_bw()
  )
}

if (any(ar$chr_class == "chrX")) {
  print(
    ggplot(ar %>% dplyr::filter(chr_class == "chrX"),
           aes(x = mid / 1e6, y = escape_frac, colour = condition)) +
      geom_point(alpha = 0.25, size = 0.7) +
      geom_smooth(se = TRUE, method = "loess", formula = y ~ x, span = 0.4) +
      labs(title = "chrX escape fraction (CAST / total) along the X",
           subtitle = "CAST is the inactive X under this genotype, so this is escape; cell types pooled",
           x = "position on chrX (Mb)", y = "escape fraction") +
      theme_bw()
  )
}

dev.off()

# ---------------------------------------------------------------------------
cat("\nWrote to ", OUT_DIR, ":\n", sep = "")
cat("  pseudobulk_celltype_ar_all_genes_raw.txt      every gene x pseudobulk, unfiltered\n")
cat(sprintf("  pseudobulk_celltype_ar_all_genes_min%dreads.txt  the filtered table\n", MIN_GENE_READS))
cat("  pseudobulk_celltype_ar_a1_wide.txt            gene x (condition, cell type)\n")
cat("  pseudobulk_celltype_qc.txt                    depth per pseudobulk -- read first\n")
cat("  pseudobulk_direction_check.txt                chrX vs autosome sanity check\n")
cat("  pseudobulk_dist_bin_summary.txt               by distance to nearest end\n")
cat("  pseudobulk_dist_bin_depth_matched.txt         same, depth-stratified\n")
cat("  pseudobulk_end_type_summary.txt               proximal vs distal ends\n")
cat("  pseudobulk_chrX_escape_by_dist.txt            chrX escape by distance\n")
cat(sprintf("  pseudobulk_end_vs_matched_interior%s.txt   terminal vs matched interior genes\n", TERM_TAG))
cat(sprintf("  pseudobulk_end_contrast_coverage%s.txt     which contrasts were testable at all\n", TERM_TAG))
cat(sprintf("  pseudobulk_end_sign_consistency%s.txt      does the sign repeat across cell types\n", TERM_TAG))
cat(sprintf("  pseudobulk_celltype_chromosome_ends%s.pdf  figures\n", TERM_TAG))
cat("\nn=1 animal per condition: cohort differences above are descriptive.\n")
cat("The terminal-vs-interior contrasts are within-pseudobulk and matched on\n")
cat("depth and SNP count. Read them as effect sizes; the permutation p ignores\n")
cat("the correlation between neighbouring genes, so sign agreement across cell\n")
cat("types is the check that matters.\n")
