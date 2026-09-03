# Where does the chrX CAST signal actually come from?
#
# WHY THIS EXISTS. The chromosome-wide chrX escape number is 12.7% in both
# sections. Splitting the same molecules by whether they sit inside an
# annotated gene body takes it apart:
#
#                        9w      78w
#   chrX, whole        0.1272  0.1264
#   chrX, gene bodies  0.0612  0.0613
#   chrX, outside      0.4750  0.4981     <- biallelic
#   autosomes, whole   0.4821  0.4841
#
# Non-genic chrX molecules are 16% of chrX but carry 60% of its CAST signal,
# and they sit at the autosomal value. That is not escape from an inactive X -
# in this genotype the B6 X is active in every cell, so a locus reading 50/50,
# let alone 100% CAST, is reporting something other than allelic expression.
# The autosomal genic/non-genic split barely moves (0.482 vs 0.485), so this is
# chrX-specific and not a general property of the annotation.
#
# This script localises that signal and produces a mask, in two passes that can
# be run at different times:
#
#   PASS 1, WINDOWS. Runs off the --window-out table the counting pass already
#     wrote (ase/<sample>/chrX_windows_64um.tsv.gz). No BAM pass, minutes.
#     Ranks 100kb windows by CAST molecules and annotates each with how much of
#     it is gene body. This is enough to name the offending regions.
#   PASS 2, SNPs. Needs the --snp-out ledger, which is a fresh counting pass
#     (slurm/spatial_artifact_scan.slurm, mode `count`). Drops to the SNP and
#     asks which individual SNPs carry a window's or a gene's skew, and writes
#     the mask bed. Skipped with a message if the ledger is absent.
#
# WHAT CALIBRATES THE CALL. The autosomes. They are biallelic by construction
# in this cross, so an autosomal SNP or window sitting near 0 or 1 cannot be
# biology - it is the artefact, measured directly, on the same reads, by the
# same counter. The autosomal rate at a given depth and threshold is therefore
# the false-positive rate for the identical call made on chrX. Nothing here
# thresholds chrX against a number chosen by eye.
#
# THE HYPOTHESIS BEING TESTED, for the flipped SNPs specifically. The SNP bed
# is C57BL/6NJ x CAST/EiJ but the reference is GRCm39, which is C57BL/6J. At
# sites where 6NJ carries a non-reference allele, the bed's "B6" allele is not
# the reference base, so a genuine B6 read scores as CAST. That predicts (a)
# flipped SNPs cluster in blocks rather than scattering, since strain-private
# haplotypes are blocks, and (b) at a flipped SNP the observed base is
# overwhelmingly the declared ALT with almost nothing on the declared REF -
# a clean 1.0, not a noisy 0.7. Both are checked below. The alternative -
# paralogous reads piling in from elsewhere - predicts intermediate fractions
# and enrichment in multicopy families instead. The script reports which
# pattern each region matches rather than assuming one.
#
# Run on the cluster:  sbatch ~/Postdoc/slurm/spatial_artifact_scan.slurm
#              or:     SAMPLE=9w Rscript ~/Postdoc/spatial/ase_artifact_scan.R

SPASE_DIR <- Sys.getenv("SPASE_HOME", "")
if (!nzchar(SPASE_DIR)) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  SPASE_DIR <- if (length(f)) dirname(normalizePath(f)) else getwd()
}
source(file.path(SPASE_DIR, "spase_common.R"))

WINDOWS   <- Sys.getenv("WINDOWS",
  file.path(ASE_DIR, sprintf("chrX_windows_64um%s.tsv.gz", SUF)))
SNP_LEDGER <- Sys.getenv("SNP_LEDGER",
  file.path(ASE_DIR, sprintf("snp_ledger%s.tsv.gz", SUF)))
GENE_BED  <- Sys.getenv("GENE_BED",
  "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/GRCm39/annotation_us_mm39_gene_level.bed")

# 100kb, matching the counting pass's --window-size default. Read from the
# table itself below where possible; this is only the fallback.
WIN_SIZE  <- as.integer(Sys.getenv("WINDOW_SIZE", "100000"))

# A window or SNP needs this many molecules / observations before its allelic
# fraction is worth looking at. Both are deliberately low: the point is to
# catch artefacts, and the autosomal arm reports what that costs in false
# positives at exactly these depths.
MIN_WIN_UMI <- as.integer(Sys.getenv("MIN_WIN_UMI", "50"))
MIN_SNP_OBS <- as.integer(Sys.getenv("MIN_SNP_OBS", "20"))

# The alt fraction above which a SNP is called flipped. Not tuned: 0.90 is far
# outside anything an autosomal SNP can reach honestly, and the autosomal false
# positive rate at this threshold is measured and printed rather than assumed.
FLIP_ALT <- as.numeric(Sys.getenv("FLIP_ALT", "0.90"))

spase_theme()

##### ------------------- GENE BODY COVERAGE ------------------- #####
if (!file.exists(GENE_BED))
  stop("No gene bed at ", GENE_BED, "\n  Rebuild it with OCM_heart/",
       "core_escape_SNPs.R - see slurm/spatial_spase.slurm.", call. = FALSE)
gb <- fread(GENE_BED, header = FALSE,
            select = 1:4, col.names = c("chrom", "start", "end", "gene"))
cat(sprintf("Gene bed: %d intervals, %d on %s\n",
            nrow(gb), sum(gb$chrom == X_CHROM), X_CHROM))

# Genic base pairs per window. Overlapping gene bodies must not be
# double-counted or a window comes out >100% genic, so the intervals are
# flattened per chromosome first.
genic_bp_per_window <- function(bed, win) {
  setorder(bed, chrom, start, end)
  flat <- bed[, {
    s <- start; e <- end; keep_s <- s[1]; keep_e <- e[1]
    os <- integer(0); oe <- integer(0)
    if (.N > 1) for (i in 2:.N) {
      if (s[i] <= keep_e) { keep_e <- max(keep_e, e[i]) }
      else { os <- c(os, keep_s); oe <- c(oe, keep_e); keep_s <- s[i]; keep_e <- e[i] }
    }
    list(start = c(os, keep_s), end = c(oe, keep_e))
  }, by = chrom]
  # One row per (interval, window it touches), clipped to the window.
  flat[, iid := .I]
  parts <- flat[, .(w = (start %/% win):(end %/% win)), by = .(iid, chrom, start, end)]
  parts[, `:=`(a = pmax(start, w * win), b = pmin(end, (w + 1) * win))]
  parts[b > a, .(genic_bp = sum(b - a)), by = .(chrom, w)]
}
gcov <- genic_bp_per_window(copy(gb), WIN_SIZE)

##### ----------------------- PASS 1: WINDOWS ----------------------- #####
if (!file.exists(WINDOWS)) {
  cat("\nNo window table at", WINDOWS, "- skipping pass 1.\n",
      " Produce it with ase_bin_allele_counts.py --window-out.\n")
  win <- NULL
} else {
  w <- fread(WINDOWS)
  setnames(w, tolower(names(w)))
  need <- c("chrom", "win_start", "ref", "alt")
  if (!all(need %in% names(w)))
    stop("Expected ", paste(need, collapse = ", "), " in ", WINDOWS,
         "; got ", paste(names(w), collapse = ", "), call. = FALSE)
  # Trust the file over the env var: a table written at a different
  # --window-size would silently mis-join to the gene coverage otherwise.
  st <- sort(unique(w$win_start))
  if (length(st) > 1) {
    d <- min(diff(st))
    if (d != WIN_SIZE) {
      cat(sprintf("Window table uses %d bp windows, not %d - using %d.\n",
                  d, WIN_SIZE, d))
      WIN_SIZE <- d
      gcov <- genic_bp_per_window(copy(gb), WIN_SIZE)
    }
  }

  win <- w[, .(ref = sum(ref), alt = sum(alt)), by = .(chrom, win_start)]
  win[, w := win_start %/% WIN_SIZE]
  win <- merge(win, gcov, by = c("chrom", "w"), all.x = TRUE)
  win[is.na(genic_bp), genic_bp := 0]
  win[, `:=`(n = ref + alt, genic_frac = genic_bp / WIN_SIZE)]
  win[, alt_frac := alt / n]
  win[, is_x := chrom == X_CHROM]

  # Which genes a window touches, for reading the table.
  gb[, `:=`(w_lo = start %/% WIN_SIZE, w_hi = end %/% WIN_SIZE)]
  gnames <- gb[, .(w = w_lo:w_hi), by = .(chrom, gene)][
    , .(genes = paste(head(unique(gene), 4), collapse = ",")), by = .(chrom, w)]
  win <- merge(win, gnames, by = c("chrom", "w"), all.x = TRUE)
  win[is.na(genes), genes := ""]

  tot <- win[is_x == TRUE, .(ref = sum(ref), alt = sum(alt))]
  cat(sprintf("\n=== PASS 1: %s windows (%d kb) ===\n%s total ref %d alt %d, "
              "alt fraction %.4f\n",
              X_CHROM, WIN_SIZE %/% 1000, X_CHROM, tot$ref, tot$alt,
              tot$alt / (tot$ref + tot$alt)))

  setorder(win, -alt)
  top <- win[is_x == TRUE][seq_len(min(25, .N))]
  top[, cum_pct := 100 * cumsum(alt) / tot$alt]
  cat("\nWindows carrying the most CAST molecules:\n")
  print(top[, .(window = sprintf("%.2f-%.2f Mb", win_start / 1e6,
                                 (win_start + WIN_SIZE) / 1e6),
               ref, alt, alt_frac = round(alt_frac, 3),
               genic = round(genic_frac, 2), cum_pct = round(cum_pct, 1),
               genes)], nrows = 30)

  # The headline decomposition, and the autosomal version of it as the control.
  dec <- win[n >= MIN_WIN_UMI, .(
    windows = .N, ref = sum(ref), alt = sum(alt),
    alt_frac = sum(alt) / sum(ref + alt)),
    by = .(chrom = ifelse(is_x, X_CHROM, "autosomes"),
           mostly_genic = genic_frac >= 0.5)]
  setorder(dec, chrom, -mostly_genic)
  cat("\nGenic vs non-genic, windows with >=", MIN_WIN_UMI, " molecules:\n",
      sep = "")
  print(dec)
  cat("\nThe autosomal rows are the control: they should agree with each other.",
      "\nIf the chrX rows do not, the chrX difference is positional, not allelic.\n")

  # How concentrated is it? A handful of loci accounting for most of the signal
  # is a mask; a smooth spread over the chromosome is biology.
  xw <- win[is_x == TRUE][order(-alt)]
  n_half <- which(cumsum(xw$alt) >= 0.5 * tot$alt)[1]
  n_90 <- which(cumsum(xw$alt) >= 0.9 * tot$alt)[1]
  cat(sprintf("\nConcentration: %d window(s) carry half of all %s CAST "
              "molecules, %d carry 90%%, out of %d windows with any.\n",
              n_half, X_CHROM, n_90, nrow(xw[alt > 0])))

  fwrite(win[order(chrom, win_start)],
         file.path(OUT_DIR, sprintf("artifact_windows_%s%s.tsv", SAMPLE, SUF)),
         sep = "\t")

  # Figure: every window, alt fraction against depth, chrX vs autosome. An
  # artefact is a high-depth point pinned near 1.0; escape is a spread near the
  # bottom. Slots 1/2 as everywhere else in this project.
  pw <- ggplot(win[n >= MIN_WIN_UMI],
               aes(n, alt_frac, colour = ifelse(is_x, "chrX", "autosomal"))) +
    geom_hline(yintercept = 0.5, colour = "grey75", linewidth = 0.3) +
    geom_point(size = 1.1, alpha = 0.6) +
    scale_x_log10() + scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = c(chrX = C_X, autosomal = C_A)) +
    labs(x = "molecules in window (log)", y = "CAST fraction", colour = NULL,
         title = sprintf("Per-window allelic fraction - %s", SAMPLE),
         subtitle = wrap(sprintf(paste0(
           "%d kb windows with >=%d molecules. The autosomal cloud sits at ",
           "0.5 by construction. chrX points near 1.0 are not escape: in this ",
           "genotype the B6 X is active in every cell, so a locus cannot be ",
           "CAST-only. Mask %s."), WIN_SIZE %/% 1000, MIN_WIN_UMI, SNP_LABEL))) +
    theme(legend.position = "top")
}

##### ------------------------ PASS 2: SNPs ------------------------ #####
snp <- NULL
if (!file.exists(SNP_LEDGER)) {
  cat("\n=== PASS 2 skipped: no SNP ledger at", SNP_LEDGER, "===\n",
      " It needs a counting pass with --snp-out:\n",
      "   sbatch ~/Postdoc/slurm/spatial_artifact_scan.slurm count\n",
      " Pass 1 above already names the regions; pass 2 names the SNPs inside\n",
      " them and writes the mask.\n")
} else {
  snp <- fread(SNP_LEDGER)
  setnames(snp, tolower(names(snp)))
  need <- c("chrom", "pos", "ref_allele", "alt_allele", "obs_ref", "obs_alt")
  if (!all(need %in% names(snp)))
    stop("Expected ", paste(need, collapse = ", "), " in ", SNP_LEDGER,
         "; got ", paste(names(snp), collapse = ", "), call. = FALSE)
  snp[, `:=`(obs = obs_ref + obs_alt, is_x = chrom == X_CHROM)]
  snp[, alt_frac := obs_alt / obs]
  if (!"gene" %in% names(snp)) snp[, gene := ""]
  snp[is.na(gene), gene := ""]

  cat(sprintf("\n=== PASS 2: %d informative SNPs (%d on %s, %d autosomal) ===\n",
              nrow(snp), sum(snp$is_x), X_CHROM, sum(!snp$is_x)))

  ok <- snp[obs >= MIN_SNP_OBS]
  # THE CALIBRATION. On an autosome the truth is 0.5, so every autosomal SNP
  # over the threshold is a false positive by definition. That rate is what the
  # chrX count has to be read against.
  fp_auto <- ok[is_x == FALSE, mean(alt_frac >= FLIP_ALT)]
  n_x_flip <- ok[is_x == TRUE, sum(alt_frac >= FLIP_ALT)]
  n_x <- ok[is_x == TRUE, .N]
  cat(sprintf(paste0(
    "SNPs with >=%d observations: %d autosomal, %d on %s.\n",
    "  autosomal SNPs at CAST fraction >=%.2f: %.4f  <- false positive rate\n",
    "  %s SNPs at the same threshold:          %.4f  (%d of %d)\n",
    "  excess over the autosomal rate: %.4f, i.e. ~%.0f of the %d are\n",
    "  more than the calibration explains.\n"),
    MIN_SNP_OBS, ok[is_x == FALSE, .N], n_x, X_CHROM, FLIP_ALT,
    fp_auto, X_CHROM, n_x_flip / n_x, n_x_flip, n_x,
    n_x_flip / n_x - fp_auto, n_x_flip - fp_auto * n_x, n_x_flip))

  # Is a flipped SNP clean or noisy? The 6NJ-vs-6J hypothesis predicts almost
  # nothing on the declared ref allele; paralogous pile-up predicts a mixture.
  fl <- ok[is_x == TRUE & alt_frac >= FLIP_ALT]
  if (nrow(fl)) {
    cat(sprintf(paste0(
      "\nShape of the %s flips: median CAST fraction %.4f, %d of %d at ",
      "exactly 1.0 (no\n  observation at all on the declared B6 allele). ",
      "A clean 1.0 is the signature of a\n  mis-declared allele; an ",
      "intermediate value is the signature of reads arriving\n  from a ",
      "paralogous locus.\n"),
      X_CHROM, median(fl$alt_frac), sum(fl$obs_ref == 0), nrow(fl)))

    # Do they cluster? Strain-private haplotype blocks are blocks; mapping
    # noise is scattered. Compared against the flips' own null - the spacing of
    # the SNPs that were testable, not of the genome.
    setorder(fl, chrom, pos)
    gaps <- fl[, diff(pos), by = chrom]$V1
    base <- ok[is_x == TRUE][order(pos), diff(pos)]
    if (length(gaps) > 5)
      cat(sprintf(paste0(
        "  Spacing: median gap between consecutive flipped SNPs %.0f bp, ",
        "against %.0f bp\n  between consecutive testable SNPs. %s\n"),
        median(gaps), median(base),
        if (median(gaps) < 0.5 * median(base))
          "They cluster - consistent with haplotype blocks."
        else "They do not obviously cluster."))
  }

  ##### -------------- WHICH SNPs CARRY THE SKEWED GENES -------------- #####
  # The genes whose fitted escape came out above 45%, which cannot be escape.
  # Named from the data rather than hard-coded, so this keeps working when the
  # mask changes.
  per_gene <- snp[is_x == TRUE & nzchar(gene),
                  .(snps = .N, obs = sum(obs), obs_alt = sum(obs_alt),
                    mol = sum(mol_ref + mol_alt),
                    alt_frac = sum(obs_alt) / sum(obs)), by = gene]
  susp <- per_gene[obs >= 50 & alt_frac > 0.45][order(-obs_alt)]
  if (nrow(susp)) {
    cat(sprintf("\n%d %s gene(s) above 45%% CAST - impossible as escape:\n",
                nrow(susp), X_CHROM))
    print(susp)
    cat("\nPer-SNP breakdown. `top_snp_share` is the fraction of the gene's",
        "CAST\nobservations carried by its single worst SNP: near 1.0 means",
        "one bad SNP,\nnot a bad gene.\n")
    brk <- snp[is_x == TRUE & gene %in% susp$gene & obs >= 5][
      order(gene, -obs_alt)]
    shr <- brk[, .(n_snps = .N,
                   n_flipped = sum(alt_frac >= FLIP_ALT),
                   top_snp_pos = pos[1],
                   top_snp_alleles = paste0(ref_allele[1], "/", alt_allele[1]),
                   top_snp_share = obs_alt[1] / sum(obs_alt)), by = gene]
    print(shr[order(-top_snp_share)])
  } else {
    cat("\nNo chrX gene above 45% CAST in this ledger.\n")
  }

  ##### ---------------------- THE MASK ---------------------- #####
  # Every SNP the calibration says cannot be honest. Written as a bed in the
  # same convention as the SNP file it subtracts from, so it can be fed
  # straight to the mask-building step rather than re-derived.
  mask <- ok[alt_frac >= FLIP_ALT, .(chrom, pos, obs, obs_ref, obs_alt,
                                     alt_frac, gene, is_x)]
  bed <- mask[, .(chrom, start = pos, end = pos + 1L,
                  name = ifelse(nzchar(gene), gene, "."),
                  score = obs, strand = ".")]
  bed_path <- file.path(OUT_DIR,
    sprintf("artifact_snps_%s%s.bed", SAMPLE, SUF))
  fwrite(bed[order(chrom, start)], bed_path, sep = "\t", col.names = FALSE)
  fwrite(mask[order(chrom, pos)],
         file.path(OUT_DIR, sprintf("artifact_snps_%s%s.tsv", SAMPLE, SUF)),
         sep = "\t")
  cat(sprintf(paste0("\nWrote %d masked SNPs (%d on %s, %d autosomal) to\n  %s\n",
              "Both arms are in the file on purpose: the autosomal ones are ",
              "the false\npositives this threshold costs, and dropping them ",
              "from the mask would hide\nthat cost. Intersect the SNP bed ",
              "against this to build the next mask.\n"),
              nrow(bed), sum(mask$is_x), X_CHROM, sum(!mask$is_x), bed_path))

  # What the chromosome looks like once they are gone. This is the number the
  # per-gene analysis should be re-run against.
  cln <- snp[!(paste(chrom, pos) %in% mask[, paste(chrom, pos)])]
  for (cc in c(TRUE, FALSE)) {
    a <- snp[is_x == cc, sum(mol_alt)]; r <- snp[is_x == cc, sum(mol_ref)]
    a2 <- cln[is_x == cc, sum(mol_alt)]; r2 <- cln[is_x == cc, sum(mol_ref)]
    cat(sprintf("%-10s CAST fraction %.4f before mask, %.4f after (%d of %d "
                "molecules dropped)\n",
                if (cc) X_CHROM else "autosomes",
                a / (a + r), a2 / (a2 + r2), (a + r) - (a2 + r2), a + r))
  }
}

##### ------------------------- PROVENANCE ------------------------- #####
side <- file.path(OUT_DIR, sprintf("artifact_scan_%s%s.provenance.tsv",
                                   SAMPLE, SUF))
writeLines(c("key\tvalue",
  paste0("sample\t", SAMPLE),
  paste0("snp_label\t", SNP_LABEL),
  paste0("windows\t", WINDOWS),
  paste0("windows_md5\t", if (file.exists(WINDOWS))
    tools::md5sum(WINDOWS)[[1]] else ""),
  paste0("window_size\t", WIN_SIZE),
  paste0("snp_ledger\t", SNP_LEDGER),
  paste0("snp_ledger_md5\t", if (file.exists(SNP_LEDGER))
    tools::md5sum(SNP_LEDGER)[[1]] else ""),
  paste0("gene_bed\t", GENE_BED),
  paste0("gene_bed_md5\t", tools::md5sum(GENE_BED)[[1]]),
  paste0("min_win_umi\t", MIN_WIN_UMI),
  paste0("min_snp_obs\t", MIN_SNP_OBS),
  paste0("flip_alt\t", FLIP_ALT),
  paste0("n_masked_snps\t", if (is.null(snp)) "" else
    nrow(snp[obs >= MIN_SNP_OBS & alt_frac >= FLIP_ALT])),
  paste0("R_version\t", paste(R.version$major, R.version$minor, sep = ".")),
  paste0("run_utc\t", format(Sys.time(), tz = "UTC"))), side)
cat("\nProvenance written to", side, "\n")

##### --------------------------- FIGURES --------------------------- #####
pdf_path <- file.path(OUT_DIR,
  sprintf("artifact_scan_%s%s.pdf", SAMPLE, SUF))
pdf(pdf_path, width = 7.5, height = 5.5)
if (!is.null(win)) {
  print(pw)
  # Position along chrX. An artefact is a spike; escape is spread.
  print(ggplot(win[is_x == TRUE & n >= MIN_WIN_UMI],
               aes(win_start / 1e6, alt_frac, size = n)) +
    geom_hline(yintercept = 0.5, colour = "grey75", linewidth = 0.3) +
    geom_point(colour = C_X, alpha = 0.6) +
    scale_size_area(max_size = 5) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = sprintf("%s position (Mb)", X_CHROM), y = "CAST fraction",
         size = "molecules",
         title = sprintf("CAST fraction along %s - %s", X_CHROM, SAMPLE),
         subtitle = wrap(paste0(
           "Escape is spread along the chromosome at a low level. A large ",
           "point near 1.0 is a locus where the allele call itself is wrong, ",
           "and it dominates the chromosome-wide aggregate because it is both ",
           "deep and extreme."))))
}
if (!is.null(snp)) {
  s2 <- snp[obs >= MIN_SNP_OBS]
  print(ggplot(s2, aes(alt_frac, after_stat(scaled),
                       fill = ifelse(is_x, "chrX", "autosomal"),
                       colour = ifelse(is_x, "chrX", "autosomal"))) +
    geom_density(alpha = 0.35, linewidth = 0.5) +
    geom_vline(xintercept = FLIP_ALT, colour = "grey45", linetype = "dotted",
               linewidth = 0.3) +
    scale_fill_manual(values = c(chrX = C_X, autosomal = C_A)) +
    scale_colour_manual(values = c(chrX = C_X, autosomal = C_A)) +
    labs(x = "CAST fraction at the SNP", y = "density (scaled to 1)",
         fill = NULL, colour = NULL,
         title = sprintf("Per-SNP allelic fraction - %s", SAMPLE),
         subtitle = wrap(sprintf(paste0(
           "SNPs with >=%d base observations. Dotted = the flip threshold ",
           "%.2f. Autosomal mass beyond it is the false positive rate; chrX ",
           "mass beyond it in excess of that is the artefact."),
           MIN_SNP_OBS, FLIP_ALT))) +
    theme(legend.position = "top"))
}
dev.off()
cat("Wrote", pdf_path, "\n")
cat("\nn = 1 animal per age: do not read a difference between the two",
    "sections as an age effect.\n")
