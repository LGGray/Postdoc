# QC layer for the per-gene chrX escape screen: SNP support, counting level,
# tile recurrence, and between-age concordance.
#
# This does NOT replace spatial/escape_genes.R. That script is the screen - it
# pools molecules per gene, tests against the chromosome's own escape rate,
# applies BH and the AR < 0.90 effect floor, and already holds out the
# above-biallelic genes. Everything here is the four checks it does not make,
# written so its gene list can be filtered by them:
#
#   1. SNP SUPPORT. escape_genes.R has no idea how many informative SNPs a gene
#      owns, so a gene resting on one SNP is ranked beside a gene resting on a
#      thousand. Llph-ps2 accumulates ~300 informative molecules from a SINGLE
#      chrX SNP; its header calls that multimapping, which is a good guess, but
#      one mis-specified SNP produces the identical signature and neither can be
#      told from the pooled ratio alone. MIN_SNP is a hard floor instead.
#
#   2. COUNTING LEVEL. tile_gene_counts.tsv carries a1_umi/n_umi AND
#      a1_reads/n_reads. escape_genes.R reads only the UMI columns, which is
#      correct and safe by construction - but nothing in the repo records how
#      wrong the read columns are for this purpose, so the next person to reach
#      for "depth > 20" reaches for the wrong one. In the --keep-duplicates run
#      a read-level depth filter of 20 is satisfied by ONE OR TWO molecules
#      (9 of 85,973 gene x tile rows reach 20 molecules in 9w; 22,943 reach 20
#      reads), and duplication is allele-asymmetric - B6 duplicates ~17x against
#      CAST ~10x in 9w - so pooling reads biases escape DOWNWARD. This script
#      quantifies both on the data at hand rather than asserting it.
#
#   3. TILE RECURRENCE. "Escapes in multiple tiles" is a fair question, but a
#      per-tile RATIO needs depth a 64 um tile does not have (median 37 chrX
#      molecules over ~22 genes; 74% of gene x tile rows are a single molecule).
#      A per-tile PRESENCE call needs one molecule, so recurrence is counted as
#      tiles carrying >= 1 CAST molecule out of tiles where the gene is seen at
#      all. That is the strongest per-tile statement the depth supports.
#
#   4. BETWEEN-AGE CONCORDANCE. escape_genes.R requires a gene be testable in
#      both samples, which is the right gate, but it does not flag a gene that
#      reaches escape in one age while sitting at AR >= MAX_AR in the other at
#      comparable depth. With n = 1 animal per age that pattern is a library
#      difference far more often than an ageing effect, and it should be visible
#      in the table rather than inferred from two rows.
#
# Read-only with respect to the tile run. Run on the cluster:
#   sbatch ~/Postdoc/slurm/spatial_gene_escape_qc.slurm 64 pysam_dup
# or interactively:
#   conda activate seurat_env
#   IN_DIR=$BASE/adult_aged_spatial/ase_pysam_dup_64um \
#     Rscript ~/Postdoc/spatial/chrX_escape_qc.R
#
# GENOTYPE, which sets the direction of every number here, as in escape_genes.R:
# B6 carries the Xist deletion so the CAST X is the inactive X in every cell.
# a1 = B6 = active X, a2 = CAST = inactive X, AR = a1/n, escape = 1 - AR.

.libPaths(c("~/R/matrix-dev", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

##### ---------------------- CONFIG ---------------------- #####
BASE      <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2"
SPATIAL   <- file.path(BASE, "adult_aged_spatial")
TILE_UM   <- as.integer(Sys.getenv("TILE_UM", "64"))
SNP_LABEL <- Sys.getenv("SNP_LABEL", "pysam_dup")
SAMPLES   <- strsplit(Sys.getenv("SAMPLES", "9w,78w"), ",")[[1]]

# Directory holding <sample>/tile_gene_counts.tsv[.gz]. Default follows the
# ase_pysam[_dup]_<um>um naming the tile runs use.
IN_DIR <- Sys.getenv("IN_DIR",
                     file.path(SPATIAL, sprintf("ase_%s_%dum", SNP_LABEL, TILE_UM)))

# The SNP bed the tile run was scored against. Only used to COUNT informative
# SNPs per gene, so its exact coordinate convention does not matter here - an
# off-by-one cannot change a count of SNPs inside a multi-kb gene.
SNP_BED <- Sys.getenv("SNP_BED",
                      file.path(BASE, "GRCm39",
                                "SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed"))
# chrX gene intervals. Multiple rows per gene (transcript variants) are expected
# and handled: a SNP is counted once per gene however many transcripts it hits.
ANNOT <- Sys.getenv("ANNOT",
                    file.path(BASE, "GRCm39", "annotation_us_mm39_chrX.bed"))
# The lab's empirical core escape list, used only to report recovery - never to
# filter. A screen that can only find genes already on the list is not a screen.
CORE_BED <- Sys.getenv("CORE_BED",
                       file.path(BASE, "GRCm39", "annotation_us_mm39_core_escape.bed"))

OUT_PDF <- Sys.getenv("OUT_PDF",
                      file.path(IN_DIR, sprintf("chrX_escape_qc_%dum.pdf", TILE_UM)))

# Pooled informative MOLECULES per gene per sample before a gene is reported.
# escape_genes.R uses 30. 20 is used here deliberately, and the difference is
# not cosmetic: Ddx3x (28), Jpx (22) and Eif2s3x (21) are core escape genes that
# clear 20 and not 30, so the stricter floor drops three of the eight core genes
# this data recovers. The cost is wider intervals, which is why the interval is
# printed beside every ratio and why REQUIRE_BOTH_AGES is available. Raise to 30
# to reproduce escape_genes.R's gene set exactly.
MIN_UMI <- as.integer(Sys.getenv("MIN_UMI", "20"))

# Informative SNPs a gene must own. 20 keeps Ddx3x (26 SNPs, real) and drops
# Llph-ps2 (1), Gm6829 (12), Gm14719 (7) and Ndufb11 (10) - the last being the
# gene that tops the read-level list and is not escape at all. There is no
# principled optimum; this is the value at which the known artefacts leave and
# the known escapees stay on this dataset, so re-check it if the SNP bed changes.
MIN_SNP <- as.integer(Sys.getenv("MIN_SNP", "20"))

# Effect-size floor on the AR scale, matching escape_genes.R's MAX_AR: AR < 0.90
# is "at least 10% of molecules come from the inactive X".
#
# WORTH KNOWING BEFORE TRUSTING IT: the chrX pooled molecule-level AR is ~0.90
# itself, so this cut sits essentially ON the chromosome mean rather than above
# it. It is an effect-size floor, not a significance test - escape_genes.R
# supplies the test, against the chromosome's own rate, which is the right null.
MAX_AR <- as.numeric(Sys.getenv("MAX_AR", "0.90"))

# Escape has a ceiling at biallelic expression, so AR below this needs the
# ACTIVE B6 allele to be silent, deleted or unmappable - which conventional
# escape cannot produce. Same idea as escape_genes.R's MAX_ESCAPE (0.60 on the
# escape scale = 0.40 on the AR scale); stated on the AR scale here so every
# threshold in this file lives on one scale.
MIN_AR_ESCAPE <- as.numeric(Sys.getenv("MIN_AR_ESCAPE", "0.40"))

# Require a gene to clear MIN_UMI in both ages. FALSE keeps single-age genes
# and marks them, because at MIN_UMI = 20 that is where Ddx3x and Eif2s3x (9w)
# and Jpx (78w) live. TRUE reproduces escape_genes.R's `ok` gate.
REQUIRE_BOTH_AGES <- as.logical(Sys.getenv("REQUIRE_BOTH_AGES", "FALSE"))

msg <- function(...) cat(sprintf(...), "\n", sep = "")

##### ----------------------- LOAD ----------------------- #####
read_one <- function(s) {
  f <- file.path(IN_DIR, s, "tile_gene_counts.tsv")
  if (!file.exists(f) && file.exists(paste0(f, ".gz"))) f <- paste0(f, ".gz")
  if (!file.exists(f)) stop("no tile_gene_counts.tsv[.gz] for sample ", s, " in ", IN_DIR)
  fread(f, select = c("tile", "chrom", "gene",
                      "a1_umi", "a2_umi", "n_umi",
                      "a1_reads", "a2_reads", "n_reads"))[, sample := s][]
}
msg("Reading %s", IN_DIR)
all <- rbindlist(lapply(SAMPLES, read_one))
all[, sample := factor(sample, levels = SAMPLES)]

x <- all[chrom == "chrX"]
msg("chrX gene x tile rows: %s",
    paste(sprintf("%s=%d", SAMPLES, x[, .N, keyby = sample]$N), collapse = ", "))

# Autosomal ratio per sample: both the biallelic reference and the mapping-bias
# estimate, kept because full escape lands here rather than at 0.50.
auto <- all[chrom %chin% paste0("chr", 1:19),
            .(auto_ar = sum(a1_umi) / sum(n_umi),
              auto_ar_reads = sum(a1_reads) / sum(n_reads)), by = sample]
msg("\nAutosomal ratio (= what full escape looks like):")
print(auto[, .(sample, molecules = round(auto_ar, 4), reads = round(auto_ar_reads, 4))])

##### -------- CHECK 2: what the counting level costs -------- #####
# Done first because it decides whether the read columns may be used at all
# downstream. They are not used, but the number belongs in the log next to the
# data it was measured on rather than in a comment.
DEPTH_TEST <- 20L
dg <- x[, .(rows = .N,
            umi_pass  = sum(n_umi   > DEPTH_TEST),
            read_pass = sum(n_reads > DEPTH_TEST),
            umi_med   = as.numeric(median(n_umi)),
            read_med  = as.numeric(median(n_reads)),
            umi_max   = max(n_umi),
            read_max  = max(n_reads),
            single    = mean(n_umi == 1L),
            # Allele-asymmetric duplication is the mechanism that biases a
            # read-level escape estimate downward, so measure it directly.
            dup_a1    = sum(a1_reads) / sum(a1_umi),
            dup_a2    = sum(a2_reads) / sum(a2_umi),
            ar_umi    = sum(a1_umi) / sum(n_umi),
            ar_reads  = sum(a1_reads) / sum(n_reads)), by = sample]
msg("\n--- CHECK 2: a depth filter of %d means different things per level ---", DEPTH_TEST)
print(dg[, .(sample, rows,
             pass_molecules = umi_pass, pass_reads = read_pass,
             med_molecules = umi_med, med_reads = read_med,
             frac_single_molecule = round(single, 3))])
msg("Per-allele duplication and what it does to the pooled chrX ratio:")
print(dg[, .(sample,
             dup_B6 = round(dup_a1, 1), dup_CAST = round(dup_a2, 1),
             asymmetry = round(dup_a1 / dup_a2, 2),
             chrX_AR_molecules = round(ar_umi, 4),
             chrX_AR_reads = round(ar_reads, 4))])
msg("Read-level pooling shifts the chrX ratio by %s. A POSITIVE shift means",
    paste(sprintf("%s %+.3f", as.character(dg$sample),
                  dg$ar_reads - dg$ar_umi), collapse = ", "))
msg("read-level pooling under-states escape. Either way, every number below")
msg("uses MOLECULES only.")

##### -------- CHECK 1: informative SNPs per gene -------- #####
msg("\n--- CHECK 1: informative SNP support per gene ---")
snp_n <- NULL
if (file.exists(SNP_BED) && file.exists(ANNOT)) {
  bed <- fread(ANNOT, select = 1:4, col.names = c("chrom", "start", "end", "gene"))
  bed <- bed[chrom == "chrX"]
  snp <- fread(SNP_BED, select = 1:2, col.names = c("chrom", "pos"))
  snp <- snp[chrom == "chrX"]
  msg("chrX SNPs in bed: %d over %d annotation intervals (%d genes)",
      nrow(snp), nrow(bed), uniqueN(bed$gene))
  snp[, `:=`(s = pos, e = pos)]
  setkey(bed, chrom, start, end)
  ov <- foverlaps(snp[, .(chrom, s, e, pos)], bed,
                  by.x = c("chrom", "s", "e"), type = "within", nomatch = NULL)
  # unique(gene, pos): a SNP inside three transcript variants of one gene is one
  # informative SNP for that gene, not three.
  snp_n <- unique(ov[, .(gene, pos)])[, .(snps = .N), by = gene]
  print(head(snp_n[order(-snps)], 5))
} else {
  warning("SNP bed or annotation not readable - MIN_SNP cannot be applied.\n",
          "  SNP_BED: ", SNP_BED, "\n  ANNOT:   ", ANNOT)
  msg("SNP bed or annotation missing; SNP filter SKIPPED (snps reported as NA).")
}

core <- character()
if (file.exists(CORE_BED)) {
  core <- unique(fread(CORE_BED, select = 4, col.names = "gene")$gene)
  msg("Core escape list: %d genes (%s)", length(core), paste(core, collapse = ", "))
}

##### --------- POOL, plus CHECK 3: recurrence --------- #####
g <- x[, .(a1 = sum(a1_umi), a2 = sum(a2_umi), n = sum(n_umi),
           # Tiles where the gene is detected at all, and of those, tiles
           # carrying at least one inactive-X molecule. One molecule is the
           # whole point: a per-tile RATIO is not supportable at this depth.
           tiles = .N,
           tiles_cast = sum(a2_umi > 0L),
           # Reads kept only so the table can show what the wrong level says.
           n_reads = sum(n_reads), ar_reads = sum(a1_reads) / sum(n_reads)),
       by = .(sample, gene)]
g[, `:=`(ar = a1 / n, escape = a2 / n, recur = tiles_cast / tiles)]

wilson <- function(k, n, z = 1.96) {
  p <- k / n; d <- 1 + z^2 / n
  c1 <- (p + z^2 / (2 * n)) / d
  h  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(lo = c1 - h, hi = c1 + h)
}
g[, c("ar_lo", "ar_hi") := wilson(a1, n)]
g <- if (is.null(snp_n)) g[, snps := NA_integer_][] else snp_n[g, on = "gene"]
# A gene in tile_gene_counts.tsv but absent from ANNOT gets no SNP count. That
# is an annotation mismatch, not a gene with no SNPs, and binning it silently
# into low_snp_support would hide a real problem - so say how many, and which.
if (!is.null(snp_n)) {
  orphan <- unique(g[is.na(snps), gene])
  if (length(orphan))
    msg("NOTE: %d chrX gene(s) counted by the tile run are absent from\n      %s - reported with snps = 0, so they land in low_snp_support.\n      Check the two annotations match before trusting that bucket: %s",
        length(orphan), ANNOT,
        paste(head(orphan, 8), collapse = ", "))
}
g[is.na(snps), snps := 0L]
g[, core := gene %chin% core]

scored <- g[n >= MIN_UMI]
msg("\nGenes reaching MIN_UMI = %d: %s",
    MIN_UMI, paste(sprintf("%s=%d", SAMPLES, scored[, uniqueN(gene), keyby = sample]$V1),
                   collapse = ", "))

# A gene enters the table if its interval clears MAX_AR in at least one age -
# the same "significant somewhere" logic as escape_genes.R, expressed on the
# interval rather than on a p-value because this script does not test.
cand_genes <- scored[ar_hi < MAX_AR, unique(gene)]
if (REQUIRE_BOTH_AGES) {
  both <- scored[, .N, by = gene][N == length(SAMPLES), gene]
  cand_genes <- intersect(cand_genes, both)
  msg("REQUIRE_BOTH_AGES = TRUE: restricted to genes scored in all %d ages.",
      length(SAMPLES))
}
msg("Genes with 95%% CI upper bound below MAX_AR = %.2f in >=1 age: %d",
    MAX_AR, length(cand_genes))

LEVELS <- c("escape", "inactive_X_exclusive", "age_discordant", "low_snp_support")
LABELS <- c(escape               = "Escape",
            inactive_X_exclusive = "Inactive-X-exclusive (not escape)",
            age_discordant       = "Age-discordant (one age only)",
            low_snp_support      = "Insufficient SNP support")

##### ------------------ CLASSIFY ------------------ #####
# Order matters: a gene is described by the FIRST reason it is not a clean
# escape call, so a 1-SNP gene is reported as unsupported rather than as
# whatever its ratio happens to look like.
cls <- scored[gene %chin% cand_genes,
              .(min_ar = min(ar), max_ar_any = max(ar), ages = .N,
                min_snps = min(snps)), by = gene]
cls[, flag := fifelse(
  min_snps < MIN_SNP,                 "low_snp_support",
  fifelse(min_ar < MIN_AR_ESCAPE,     "inactive_X_exclusive",
  fifelse(max_ar_any >= MAX_AR,       "age_discordant",
                                      "escape")))]
# CHECK 4 guard. A single-age gene cannot be discordant, only unconfirmed. It
# should already be impossible to reach age_discordant with one age - entering
# the table needs ar_hi < MAX_AR, which forces ar < MAX_AR for that same age -
# so this is a belt-and-braces line, not load-bearing. `unconfirmed` is the
# column that actually carries "only one age had the depth to say".
cls[flag == "age_discordant" & ages < length(SAMPLES), flag := "escape"]
cls[, unconfirmed := ages < length(SAMPLES)]
# Ordered by LEVELS, not alphabetically, so the TSV and the plot come out in
# the order the classes are meant to be read rather than in flag-name order.
cls[, flag := factor(flag, levels = LEVELS)]

out <- cls[scored[gene %chin% cand_genes], on = "gene", nomatch = NULL]
setorder(out, flag, -min_ar, gene, sample)

msg("\n--- CLASSIFICATION ---")
for (f in LEVELS) {
  gs <- unique(out[flag == f, gene])
  if (!length(gs)) next
  msg("\n%s: %d gene%s", LABELS[[f]], length(gs), if (length(gs) > 1) "s" else "")
  print(out[flag == f, .(gene, sample, snps, n,
                         AR = round(ar, 3),
                         CI = sprintf("[%.3f-%.3f]", ar_lo, ar_hi),
                         tiles_cast, tiles, recur = round(recur, 2),
                         AR_reads = round(ar_reads, 3),
                         core = fifelse(core, "*", ""))])
}

esc <- unique(out[flag == "escape", gene])
if (length(core)) {
  rec <- intersect(core, esc)
  msg("\nCore escape genes recovered: %d/%d (%s)",
      length(rec), length(core), paste(rec, collapse = ", "))
  miss <- setdiff(core, esc)
  if (length(miss)) {
    msg("Not recovered: %s", paste(miss, collapse = ", "))
    # Depth is the usual reason but not the only one, and a core escape gene
    # dropped for SNP support or for landing in another class is a different
    # problem entirely - so name the actual reason rather than assuming depth.
    print(g[gene %chin% miss,
            .(gene, sample, n, snps, AR = round(ar, 3),
              reason = fcase(n < MIN_UMI,  sprintf("below MIN_UMI (%d)", MIN_UMI),
                             snps < MIN_SNP, sprintf("below MIN_SNP (%d)", MIN_SNP),
                             ar_hi >= MAX_AR, sprintf("CI upper >= MAX_AR (%.2f)", MAX_AR),
                             default = "scored, but classified outside `escape`")
              )][order(gene, sample)])
  }
}

TSV <- sub("\\.pdf$", ".tsv", OUT_PDF)
fwrite(out[, .(gene, sample, flag, core, unconfirmed, snps,
               n, a1, a2, ar, ar_lo, ar_hi, escape,
               tiles, tiles_cast, recur, n_reads, ar_reads)],
       TSV, sep = "\t")

##### -------------------- THE PLOT -------------------- #####
# Three panels, in the order the argument has to be made: what the counting
# level does, then the per-gene ratios, then how many tiles each one recurs in.
pdf(OUT_PDF, width = 11, height = 7.5, onefile = TRUE)

d1 <- melt(dg[, .(sample, molecules = umi_pass, reads = read_pass)],
           id.vars = "sample", variable.name = "level", value.name = "pass")
p1 <- ggplot(d1, aes(level, pass, fill = sample)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(aes(label = format(pass, big.mark = ",", trim = TRUE)),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("#2a78d6", "#eb6834")[seq_along(SAMPLES)],
                    breaks = SAMPLES, name = NULL) +
  labs(title = sprintf("A depth filter of %d, applied to each counting level", DEPTH_TEST),
       subtitle = sprintf(paste("gene x tile observations reaching depth > %d.",
                                "%.0f%%-%.0f%% of rows are a SINGLE molecule, so a read-level",
                                "threshold of %d is\nsatisfied by one or two molecules -",
                                "and duplication is allele-asymmetric (B6 %.1fx vs CAST %.1fx in %s),",
                                "so\npooling reads under-states escape."),
                          DEPTH_TEST, 100 * min(dg$single), 100 * max(dg$single),
                          DEPTH_TEST, dg$dup_a1[1], dg$dup_a2[1], as.character(dg$sample)[1]),
       x = NULL, y = sprintf("gene x tile rows with depth > %d", DEPTH_TEST)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8, colour = "#52514e"))
print(p1)

fp <- copy(out)
fp[, flag := droplevels(factor(flag, levels = LEVELS, labels = LABELS[LEVELS]))]
fp[, glab := paste0(ifelse(core, "* ", ""), gene, ifelse(unconfirmed, " (1 age)", ""))]
setorder(fp, flag, min_ar)
fp[, glab := factor(glab, levels = unique(glab))]
p2 <- ggplot(fp, aes(ar, glab, colour = sample)) +
  geom_vline(xintercept = MIN_AR_ESCAPE, linetype = 3, colour = "#8a8a8a") +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "#8a8a8a") +
  geom_vline(xintercept = MAX_AR, linetype = 2, colour = "#eb6834") +
  geom_errorbarh(aes(xmin = ar_lo, xmax = ar_hi), height = 0, linewidth = 0.7, alpha = 0.55) +
  geom_point(size = 2.1) +
  facet_grid(flag ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = c("#2a78d6", "#eb6834")[seq_along(SAMPLES)],
                      breaks = SAMPLES, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "Pooled allelic ratio per chrX gene, molecules only",
       subtitle = paste0(
         sprintf("AR = B6 / total and B6 is the ACTIVE X, so escape sits BETWEEN 0.5 and %.2f. ", MAX_AR),
         sprintf("Below %.2f the ACTIVE allele is silent,\nwhich escape cannot produce. ", MIN_AR_ESCAPE),
         sprintf("Bars are 95%% Wilson intervals; genes need >= %d molecules and >= %d SNPs. ", MIN_UMI, MIN_SNP),
         "* = core escape gene.\nn = 1 animal per age: a screen, not an age comparison."),
       x = "allelic ratio (B6 / total)", y = NULL,
       caption = sprintf("Dotted %.2f = escape floor. Dashed 0.50 = biallelic. Orange %.2f = the effect-size cut, which sits almost exactly on the chrX pooled mean.",
                         MIN_AR_ESCAPE, MAX_AR)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 0, size = 8),
        plot.subtitle = element_text(size = 8, colour = "#52514e"),
        plot.caption = element_text(size = 7, colour = "#52514e", hjust = 0))
print(p2)

rc <- out[flag == "escape"]
if (nrow(rc)) {
  setorder(rc, -recur)
  rc[, glab := factor(paste0(ifelse(core, "* ", ""), gene),
                      levels = unique(paste0(ifelse(core, "* ", ""), gene)))]
  p3 <- ggplot(rc, aes(recur, glab, fill = sample)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.64) +
    geom_text(aes(label = sprintf("%d/%d", tiles_cast, tiles)),
              position = position_dodge(width = 0.72), hjust = -0.12, size = 2.7) +
    scale_fill_manual(values = c("#2a78d6", "#eb6834")[seq_along(SAMPLES)],
                      breaks = SAMPLES, name = NULL) +
    scale_x_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Tile recurrence of the escape calls",
         subtitle = paste("Share of tiles where the gene is detected that carry >= 1 inactive-X (CAST) molecule.",
                          "\nPresence, not ratio: one molecule is all a per-tile call can be built on at this depth."),
         x = "tiles with a CAST molecule", y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8, colour = "#52514e"))
  print(p3)
}
invisible(dev.off())

##### ------------------ PROVENANCE ------------------ #####
prov <- data.table(
  k = c("script", "run_at", "in_dir", "tile_um", "samples", "snp_label",
        "snp_bed", "snp_bed_md5", "annotation", "core_bed",
        "min_umi", "min_snp", "max_ar", "min_ar_escape", "require_both_ages",
        "chrX_AR_molecules", "chrX_AR_reads", "genes_scored", "genes_escape"),
  v = c("spatial/chrX_escape_qc.R", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        IN_DIR, TILE_UM, paste(SAMPLES, collapse = ","), SNP_LABEL,
        SNP_BED,
        if (file.exists(SNP_BED)) unname(tools::md5sum(SNP_BED)) else
          "bed not readable from here",
        ANNOT, CORE_BED,
        MIN_UMI, MIN_SNP, MAX_AR, MIN_AR_ESCAPE, REQUIRE_BOTH_AGES,
        paste(sprintf("%s=%.4f", as.character(dg$sample), dg$ar_umi), collapse = ","),
        paste(sprintf("%s=%.4f", as.character(dg$sample), dg$ar_reads), collapse = ","),
        uniqueN(scored$gene), length(esc)))
setnames(prov, c("key", "value"))
fwrite(prov, sub("\\.pdf$", "_provenance.tsv", OUT_PDF), sep = "\t")

msg("\nWrote %s\n      %s\n      %s",
    OUT_PDF, TSV, sub("\\.pdf$", "_provenance.tsv", OUT_PDF))
msg("\nTo filter escape_genes.R's ranking by SNP support, join the `snps` column")
msg("from %s on `gene` and drop rows below MIN_SNP = %d.", basename(TSV), MIN_SNP)
