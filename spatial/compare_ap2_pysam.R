# Numerical comparison of the two per-tile counting routes.
#
#   route A   sinto split -> Allelome.PRO2 per tile -> locus_table.txt
#             ase/<sample>/Allelome.PRO2_tiles_<um>um[_<label>]/
#   route B   ase_tile_locus_counts.py, one pysam pass over the BAM
#             ase_pysam_<um>um/<sample>/tile_chrom_counts.tsv
#
# The claim route B has to earn is narrow and worth stating exactly:
#
#   AT READ LEVEL, ON THE SAME TILES, THE TWO ROUTES COUNT THE SAME THING.
#
# It is a replication claim, not an agreement-between-methods claim. Both routes
# see the same reads, the same SNP bed and the same duplicate flag, so a
# disagreement is a difference in bookkeeping, not sampling noise. That is why
# nothing below is a hypothesis test: there is no null anyone believes, and a
# p-value against "these two counts of the same molecules differ by chance"
# would be theatre. The honest statistics are (i) how many tiles agree to the
# integer, (ii) how big the disagreements are where they exist, and (iii)
# whether they depend on depth - a constant fractional loss and a depth-driven
# one have completely different causes.
#
# The read-vs-molecule difference is NOT part of that claim. It is the intended
# change, and section 5 quantifies it separately. Do not read section 5 as
# disagreement between the routes; it is disagreement between counting units,
# measured within route B alone on one set of molecules.
#
# USAGE (interactive, on the cluster)
#   conda activate seurat_env
#   Rscript ~/Postdoc/spatial/compare_ap2_pysam.R                  # 64 um
#   TILE_UM=32 Rscript ~/Postdoc/spatial/compare_ap2_pysam.R       # if both exist
#   SAMPLES=9w Rscript ~/Postdoc/spatial/compare_ap2_pysam.R       # one sample
#   MIN_N=50 Rscript ~/Postdoc/spatial/compare_ap2_pysam.R         # stricter
#
# Everything is printed to stderr as it goes, so it is readable live. Two files
# are also written next to the pysam tables: the joined per-tile table, so any
# claim here can be checked by hand, and an optional PDF.

.libPaths(c("~/R/matrix-dev", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
})

##### ---------------------- CONFIG ---------------------- #####
# Same BASE and the same env-var names as tile_ratio_map.R, deliberately: this
# script has to look at exactly the trees that produced the figures.
BASE    <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial"
SAMPLES <- strsplit(Sys.getenv("SAMPLES", "9w,78w"), ",")[[1]]
TILE_UM <- as.integer(Sys.getenv("TILE_UM", "64"))
ANNOT_BASE <- "chr_annotation_mm39.bed"

SNP_LABEL <- Sys.getenv("SNP_LABEL", "no_Xist")
SUF     <- if (SNP_LABEL == "no_Xist") "" else paste0("_", SNP_LABEL)

# Where route B's tidy table lives. Set by slurm/spatial_tile_locus_map.slurm.
PYSAM_ROOT <- Sys.getenv("PYSAM_ROOT",
                         file.path(BASE, sprintf("ase_pysam_%dum", TILE_UM)))

# Ratios on a handful of counts are noise, and a correlation over them measures
# the depth distribution rather than the agreement. Tile-level ratio statistics
# are therefore reported on tiles with at least MIN_N units on BOTH sides; the
# integer-agreement counts in section 2 use every shared tile and need no cutoff.
MIN_N <- as.integer(Sys.getenv("MIN_N", "20"))

# Match tile_ratio_map.R: permanent storage is the source of record and $SCRATCH
# supplements it. If this is FALSE and the tile job was killed before its
# collect stage, tiles present in the figures would show up here as "only in
# pysam", which would be an artefact of where we looked, not a real difference.
SCRATCH_SUPPLEMENT <- !identical(tolower(Sys.getenv("SCRATCH_SUPPLEMENT", "TRUE")),
                                 "false")

OUT_TSV <- Sys.getenv("OUT_TSV",
                      file.path(PYSAM_ROOT,
                                sprintf("compare_ap2_pysam_%dum.tsv", TILE_UM)))
OUT_PDF <- Sys.getenv("OUT_PDF",
                      file.path(PYSAM_ROOT,
                                sprintf("compare_ap2_pysam_%dum.pdf", TILE_UM)))

AUTOSOMES <- paste0("chr", 1:19)

msg  <- function(...) message(sprintf(...))
rule <- function(t) msg("\n===== %s %s", t, strrep("=", max(0, 58 - nchar(t))))
# Percentages of a possibly-zero denominator, without printing NaN at the reader.
pc <- function(num, den) if (den > 0) sprintf("%.1f%%", 100 * num / den) else "n/a"

##### ------------------- ROUTE A: Allelome.PRO2 ------------------- #####

# One locus table -> one row. Deliberately keeps a1, a2 AND total separately,
# where tile_ratio_map.R's read_locus() keeps only a1 and total.
#
# That is not pedantry. tile_ratio_map.R computes its ratio as a1/total_reads,
# which equals a1/(a1+a2) only if Allelome.PRO2's total_reads column holds the
# informative reads and not every read over the locus. Nothing in this repo
# states which it is, and the two differ by the informative fraction - a factor
# of many. Section 2 checks it explicitly rather than assuming, because if
# total > a1+a2 then every published tile ratio is a1 over the wrong
# denominator, which would matter far more than anything else on this page.
read_locus_ap2 <- function(path) {
  lt <- tryCatch(fread(path, showProgress = FALSE), error = function(e) NULL)
  if (is.null(lt) || !nrow(lt)) return(NULL)
  setnames(lt, tolower(names(lt)))
  need <- c("chr", "a1_reads", "a2_reads", "total_reads")
  if (!all(need %in% names(lt))) return(NULL)
  x <- lt[chr == "chrX"]
  a <- lt[chr %in% AUTOSOMES]
  if (!nrow(x)) return(NULL)          # same inclusion rule tile_ratio_map.R uses
  data.table(
    x_a1_ap2 = sum(x$a1_reads), x_a2_ap2 = sum(x$a2_reads),
    x_tot_ap2 = sum(x$total_reads),
    a_a1_ap2 = sum(a$a1_reads), a_a2_ap2 = sum(a$a2_reads),
    a_tot_ap2 = sum(a$total_reads)
  )
}

collect_ap2 <- function(smp) {
  ase_dir   <- file.path(BASE, "ase", smp)
  final_dir <- file.path(ase_dir, sprintf("Allelome.PRO2_tiles_%dum%s",
                                          TILE_UM, SUF))
  lt_paths <- character(0)
  if (dir.exists(final_dir)) {
    lt_paths <- Sys.glob(file.path(final_dir, "*", "locus_table.txt"))
    names(lt_paths) <- basename(dirname(lt_paths))
  }
  msg("  %-5s route A: %d locus tables in %s", smp, length(lt_paths), final_dir)

  scratch <- Sys.getenv("SCRATCH")
  if (SCRATCH_SUPPLEMENT && nzchar(scratch)) {
    ap2 <- file.path(scratch, sprintf("spatial_tiles_%s_%dum", smp, TILE_UM),
                     paste0("allelome", SUF))
    sp <- Sys.glob(file.path(ap2, paste0("*_", ANNOT_BASE, "_*"),
                             "locus_table.txt"))
    if (length(sp)) {
      nm <- basename(dirname(sp))
      nm <- sub("^[0-9]{4}_[0-9]{2}_[0-9]{2}_", "", nm)
      nm <- sub(paste0("_", ANNOT_BASE, "_[0-9]+$"), "", nm)
      names(sp) <- nm
      sp <- sp[!duplicated(names(sp), fromLast = TRUE)]   # resubmitted _1.._n
      sp <- sp[!names(sp) %in% names(lt_paths)]           # permanent wins
      if (length(sp)) {
        msg("          + %d not yet collected, read from $SCRATCH", length(sp))
        lt_paths <- c(lt_paths, sp)
      }
    }
  }
  if (!length(lt_paths)) return(NULL)

  out <- rbindlist(lapply(seq_along(lt_paths), function(i) {
    r <- read_locus_ap2(lt_paths[[i]])
    if (is.null(r)) return(NULL)
    r[, tile := names(lt_paths)[i]][]
  }))
  msg("          %d of them have a chrX row", nrow(out))
  if (!nrow(out)) return(NULL)
  out[, sample := smp][]
}

##### ---------------------- ROUTE B: pysam ---------------------- #####

collect_pysam <- function(smp) {
  f <- file.path(PYSAM_ROOT, smp, "tile_chrom_counts.tsv")
  if (!file.exists(f)) {
    msg("  %-5s route B: MISSING %s", smp, f)
    return(NULL)
  }
  d <- fread(f, showProgress = FALSE)
  msg("  %-5s route B: %d tile x chromosome rows in %s", smp, nrow(d), f)
  # Pool to the same two groups route A is reduced to. The autosomal set is
  # chr1..chr19 on both sides, which is what tile_ratio_map.R pools for its null.
  x <- d[chrom == "chrX", .(x_a1_py = sum(a1_umi), x_a2_py = sum(a2_umi),
                            x_a1_rd = sum(a1_reads), x_a2_rd = sum(a2_reads)),
         by = tile]
  a <- d[chrom %in% AUTOSOMES, .(a_a1_py = sum(a1_umi), a_a2_py = sum(a2_umi),
                                 a_a1_rd = sum(a1_reads), a_a2_rd = sum(a2_reads)),
         by = tile]
  out <- merge(x, a, by = "tile", all = TRUE)
  for (j in names(out)) if (j != "tile") set(out, which(is.na(out[[j]])), j, 0L)
  # Route A only admits a tile that has a chrX row, so apply the same rule here
  # or the tile-set comparison would be counting different things.
  out <- out[(x_a1_rd + x_a2_rd) > 0 | (x_a1_py + x_a2_py) > 0]
  msg("          %d tiles with an informative chrX unit", nrow(out))
  out[, sample := smp][]
}

##### ------------------------ COLLECT ------------------------ #####

rule("0. Inputs")
msg("tile size %d um   SNP label %s   samples %s",
    TILE_UM, SNP_LABEL, paste(SAMPLES, collapse = ", "))
ap2 <- rbindlist(lapply(SAMPLES, collect_ap2), fill = TRUE)
pys <- rbindlist(lapply(SAMPLES, collect_pysam), fill = TRUE)
if (!nrow(ap2)) stop("No Allelome.PRO2 locus tables found - nothing to compare")
if (!nrow(pys)) stop("No pysam tile_chrom_counts.tsv found - run the counting pass first")

# all = TRUE, so a tile scored by only one route survives the join as a row of
# NAs on the other side. That is the point: those tiles are section 1's subject,
# and dropping them would hide exactly the disagreement worth knowing about.
d <- merge(ap2, pys, by = c("sample", "tile"), all = TRUE)
d[, `:=`(in_ap2 = !is.na(x_a1_ap2), in_py = !is.na(x_a1_rd))]

##### -------------------- 1. TILE SETS -------------------- #####

rule("1. Tile sets")
msg("A tile counts as present if it has at least one informative chrX unit,")
msg("which is the rule tile_ratio_map.R uses to accept a tile at all.\n")
for (smp in SAMPLES) {
  s <- d[sample == smp]
  if (!nrow(s)) next
  nb <- s[in_ap2 & in_py, .N]; na <- s[in_ap2 & !in_py, .N]
  np <- s[!in_ap2 & in_py, .N]
  # (in_ap2), not in_ap2: a lone symbol in i is looked up in the calling scope,
  # so the bare form asks for a variable that does not exist rather than
  # filtering on the column.
  msg("%-5s  route A %5d   route B %5d   both %5d   Jaccard %.4f",
      smp, s[(in_ap2), .N], s[(in_py), .N], nb, nb / (nb + na + np))
  msg("       only in A %4d   only in B %4d", na, np)
  if (np > 0) {
    q <- quantile(s[!in_ap2 & in_py, x_a1_rd + x_a2_rd], c(0, .5, .9, 1),
                  na.rm = TRUE)
    msg("       B-only tiles, chrX reads: min %.0f median %.0f p90 %.0f max %.0f",
        q[1], q[2], q[3], q[4])
    msg("       -> expected: these are tiles below SINTO_MIN_UMI that were never")
    msg("          submitted to sinto, so route A never had the chance to score them.")
    msg("          If their median depth is high, that explanation is wrong.")
  }
  if (na > 0) {
    q <- quantile(s[in_ap2 & !in_py, x_a1_ap2 + x_a2_ap2], c(0, .5, 1),
                  na.rm = TRUE)
    msg("       A-only tiles, chrX reads: min %.0f median %.0f max %.0f", q[1], q[2], q[3])
    msg("       -> NOT expected. Route B reads the whole BAM, so a tile route A")
    msg("          scored should always be reachable. Investigate before trusting")
    msg("          anything below: likely a tile-name or off-tissue mismatch.")
  }
}

##### ------------- 2. READ-LEVEL REPLICATION ------------- #####

# The comparison that carries the claim. Route B's read-level columns against
# route A's, on the tiles both scored.
agree <- function(s, a1a, a2a, a1b, a2b, what) {
  t <- s[!is.na(get(a1a)) & !is.na(get(a1b))]
  if (!nrow(t)) { msg("  no shared tiles for %s", what); return(invisible(NULL)) }
  A1 <- t[[a1a]]; A2 <- t[[a2a]]; B1 <- t[[a1b]]; B2 <- t[[a2b]]
  nA <- A1 + A2; nB <- B1 + B2
  same <- A1 == B1 & A2 == B2
  msg("  %s, %d shared tiles", what, nrow(t))
  msg("    identical a1 AND a2 counts        %5d  (%s)", sum(same), pc(sum(same), nrow(t)))
  msg("    pooled: A %d units ratio %.4f | B %d units ratio %.4f",
      sum(nA), sum(A1) / max(sum(nA), 1), sum(nB), sum(B1) / max(sum(nB), 1))
  msg("    total units B/A                   %.4f", sum(nB) / max(sum(nA), 1))

  # Depth agreement, per tile. Relative difference is the right scale: an
  # absolute difference of 5 units is trivial on a deep tile and total on a
  # shallow one.
  ok <- nA > 0 & nB > 0
  if (!any(ok)) { msg("    no tile has units on both sides"); return(invisible(NULL)) }
  rel <- (nB[ok] - nA[ok]) / nA[ok]
  msg("    per-tile depth: Pearson %.5f  Spearman %.5f",
      suppressWarnings(cor(nA[ok], nB[ok])),
      suppressWarnings(cor(nA[ok], nB[ok], method = "spearman")))
  msg("    relative depth difference (B-A)/A: median %+.4f  IQR [%+.4f, %+.4f]",
      median(rel), quantile(rel, .25), quantile(rel, .75))
  msg("    tiles within 1%% / 5%% on depth      %s / %s",
      pc(sum(abs(rel) <= 0.01), length(rel)), pc(sum(abs(rel) <= 0.05), length(rel)))

  # Ratio agreement, on tiles deep enough for a ratio to mean anything.
  deep <- ok & nA >= MIN_N & nB >= MIN_N
  if (sum(deep) >= 3) {
    rA <- A1[deep] / nA[deep]; rB <- B1[deep] / nB[deep]
    dd <- rB - rA
    msg("    per-tile ratio (n >= %d, %d tiles): Pearson %.5f",
        MIN_N, sum(deep), suppressWarnings(cor(rA, rB)))
    msg("      difference B-A: mean %+.5f  median %+.5f  sd %.5f  max|d| %.5f",
        mean(dd), median(dd), sd(dd), max(abs(dd)))
    msg("      95%% of tiles within                %.5f", quantile(abs(dd), 0.95))
  } else {
    msg("    fewer than 3 tiles reach n >= %d - no ratio statistics", MIN_N)
  }
  invisible(NULL)
}

rule("2. Read-level replication: route B reads vs route A")
msg("If the two routes are the same measurement, these agree to the integer.")
msg("They should, because both count reads passing the same filters over the")
msg("same SNP bed. Anything else is a bookkeeping difference to be explained.\n")
for (smp in SAMPLES) {
  s <- d[sample == smp]; if (!nrow(s)) next
  msg("--- %s ---", smp)
  agree(s, "x_a1_ap2", "x_a2_ap2", "x_a1_rd", "x_a2_rd", "chrX")
  agree(s, "a_a1_ap2", "a_a2_ap2", "a_a1_rd", "a_a2_rd", "autosomes (chr1-19)")
}

# The denominator check flagged in read_locus_ap2's comment.
rule("2b. Does Allelome.PRO2's total_reads equal a1+a2?")
msg("tile_ratio_map.R computes its ratio as a1/total_reads. That is a1/(a1+a2)")
msg("only if total_reads holds the informative reads. If it does not, every")
msg("published tile ratio has the wrong denominator - check this first.\n")
for (smp in SAMPLES) {
  s <- d[sample == smp & in_ap2]; if (!nrow(s)) next
  bad_x <- s[x_tot_ap2 != x_a1_ap2 + x_a2_ap2, .N]
  bad_a <- s[a_tot_ap2 != a_a1_ap2 + a_a2_ap2, .N]
  msg("%-5s chrX      %d of %d tiles where total != a1+a2 (%s)",
      smp, bad_x, nrow(s), pc(bad_x, nrow(s)))
  msg("      autosomes %d of %d (%s)", bad_a, nrow(s), pc(bad_a, nrow(s)))
  if (bad_x > 0) {
    ex <- s[x_tot_ap2 != x_a1_ap2 + x_a2_ap2][1]
    msg("      example: a1 %d + a2 %d = %d, total_reads %d  (ratio to total %.4f,",
        ex$x_a1_ap2, ex$x_a2_ap2, ex$x_a1_ap2 + ex$x_a2_ap2, ex$x_tot_ap2,
        ex$x_a1_ap2 / ex$x_tot_ap2)
    msg("               ratio to a1+a2 %.4f) - the two conventions differ here.",
        ex$x_a1_ap2 / (ex$x_a1_ap2 + ex$x_a2_ap2))
  }
}

##### ------------- 3. IS DISAGREEMENT DEPTH-DRIVEN? ------------- #####

rule("3. Does the disagreement depend on depth?")
msg("A constant fractional offset points at a filter one route applies and the")
msg("other does not - e.g. route B requires a UB tag on every read, route A")
msg("does not. A difference that grows or shrinks with depth points at")
msg("something structural instead, such as tiles being split differently.\n")
for (smp in SAMPLES) {
  s <- d[sample == smp & in_ap2 & in_py]
  s <- s[(x_a1_ap2 + x_a2_ap2) > 0]
  if (nrow(s) < 20) next
  s[, nA := x_a1_ap2 + x_a2_ap2][, nB := x_a1_rd + x_a2_rd]
  s[, rel := (nB - nA) / nA]
  # Tied depths collapse the decile edges; with fewer than three distinct edges
  # there are no deciles to speak of and cut() would error rather than say so.
  br <- unique(quantile(s$nA, seq(0, 1, .1)))
  if (length(br) < 3) {
    msg("--- %s: depths too tied for deciles (%d distinct edges) ---", smp, length(br))
    next
  }
  s[, dec := cut(nA, breaks = br, include.lowest = TRUE, labels = FALSE)]
  msg("--- %s, chrX, by depth decile of route A ---", smp)
  msg("   decile   tiles   median A depth   median (B-A)/A")
  print(s[!is.na(dec), .(tiles = .N,
                         median_A = as.integer(median(nA)),
                         median_rel = round(median(rel), 4)), by = dec][order(dec)])
}

##### ------------- 4. READS VS MOLECULES (route B only) ------------- #####

rule("4. Reads vs molecules, within route B")
msg("This is the intended difference, not a discrepancy. Both columns come from")
msg("one pass over one set of molecules, so the comparison is exact: the same")
msg("molecules, counted once each and then once per surviving read.\n")
for (smp in SAMPLES) {
  s <- d[sample == smp & in_py]; if (!nrow(s)) next
  xr <- sum(s$x_a1_rd + s$x_a2_rd); xu <- sum(s$x_a1_py + s$x_a2_py)
  ar <- sum(s$a_a1_rd + s$a_a2_rd); au <- sum(s$a_a1_py + s$a_a2_py)
  msg("--- %s ---", smp)
  msg("  chrX      %9d reads over %8d molecules   %.2f reads/molecule",
      xr, xu, xr / max(xu, 1))
  msg("  autosomes %9d reads over %8d molecules   %.2f reads/molecule",
      ar, au, ar / max(au, 1))
  msg("  chrX pooled Bl6 fraction:  %.4f per read   %.4f per molecule",
      sum(s$x_a1_rd) / max(xr, 1), sum(s$x_a1_py) / max(xu, 1))
  msg("  autosomal pooled (the reference-bias control, expect ~0.5):")
  msg("                             %.4f per read   %.4f per molecule",
      sum(s$a_a1_rd) / max(ar, 1), sum(s$a_a1_py) / max(au, 1))
  t <- s[(x_a1_rd + x_a2_rd) >= MIN_N & (x_a1_py + x_a2_py) >= MIN_N]
  if (nrow(t) >= 3) {
    rr <- t$x_a1_rd / (t$x_a1_rd + t$x_a2_rd)
    ru <- t$x_a1_py / (t$x_a1_py + t$x_a2_py)
    msg("  per-tile chrX ratio, %d tiles with n >= %d both ways:", nrow(t), MIN_N)
    msg("    Pearson %.4f   mean shift (molecule - read) %+.4f   sd of shift %.4f",
        suppressWarnings(cor(rr, ru)), mean(ru - rr), sd(ru - rr))
    msg("    -> a shift this size with a high correlation means the two units")
    msg("       rank tiles alike but put them at different absolute ratios.")
  }
}

##### ------- 5. PER-GENE DUPLICATION, if the gene table exists ------- #####

rule("5. Per-gene duplication multiplier (chrX)")
msg("Where the read-level average gets its weighting from. A gene with a high")
msg("reads/molecule ratio contributes far more to a read-level pooled ratio")
msg("than to a molecule-level one, which is the whole mechanism behind the")
msg("difference in section 4.\n")
for (smp in SAMPLES) {
  gz <- file.path(PYSAM_ROOT, smp, "tile_gene_counts.tsv.gz")
  if (!file.exists(gz)) { msg("%-5s no %s - skipped", smp, gz); next }
  g <- fread(cmd = sprintf("zcat %s", shQuote(gz)), showProgress = FALSE)
  gx <- g[chrom == "chrX", .(reads = sum(a1_reads + a2_reads),
                             umis  = sum(a1_umi + a2_umi),
                             a1_umi = sum(a1_umi)), by = gene]
  gx <- gx[umis >= MIN_N][, mult := reads / umis]
  if (!nrow(gx)) { msg("%-5s no chrX gene reaches %d molecules", smp, MIN_N); next }
  setorder(gx, -mult)
  msg("--- %s: %d chrX genes with >= %d molecules, multiplier %.2f - %.2f ---",
      smp, nrow(gx), MIN_N, min(gx$mult), max(gx$mult))
  msg("  most duplicated (these dominate a read-level pooled ratio):")
  print(head(gx[, .(gene, umis, reads, mult = round(mult, 1),
                    bl6_frac = round(a1_umi / umis, 3))], 8))
  msg("  share of chrX reads carried by the top 1%% of genes by multiplier: %s",
      pc(sum(head(gx, max(1L, nrow(gx) %/% 100))$reads), sum(gx$reads)))
}

##### ------------------------- OUTPUT ------------------------- #####

fwrite(d, OUT_TSV, sep = "\t")
msg("\nWrote the joined per-tile table to %s", OUT_TSV)
msg("  columns *_ap2 = route A reads, *_rd = route B reads, *_py = route B molecules")

ok <- requireNamespace("ggplot2", quietly = TRUE)
if (!ok) {
  msg("ggplot2 not available - skipping the PDF; the numbers above are the point.")
} else {
  suppressPackageStartupMessages(library(ggplot2))
  p <- d[in_ap2 & in_py]
  p[, `:=`(nA = x_a1_ap2 + x_a2_ap2, nB = x_a1_rd + x_a2_rd,
           nU = x_a1_py + x_a2_py)]
  p <- p[nA > 0 & nB > 0]
  pdf(OUT_PDF, width = 7, height = 6)
  print(ggplot(p, aes(nA, nB)) +
          geom_abline(slope = 1, intercept = 0, colour = "#c0392b") +
          geom_point(alpha = 0.25, size = 0.6) +
          scale_x_log10() + scale_y_log10() + facet_wrap(~sample) +
          labs(title = sprintf("chrX depth per %d um tile, read level", TILE_UM),
               subtitle = "Red is y = x. Both axes are reads, so agreement is the claim.",
               x = "Allelome.PRO2 reads", y = "pysam reads") +
          theme_bw(base_size = 10))
  q <- p[nA >= MIN_N & nB >= MIN_N]
  if (nrow(q) >= 3) {
    print(ggplot(q, aes(x_a1_ap2 / nA, x_a1_rd / nB)) +
            geom_abline(slope = 1, intercept = 0, colour = "#c0392b") +
            geom_point(alpha = 0.25, size = 0.6) +
            coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) + facet_wrap(~sample) +
            labs(title = sprintf("chrX Bl6 fraction per tile, read level (n >= %d)", MIN_N),
                 x = "Allelome.PRO2", y = "pysam") +
            theme_bw(base_size = 10))
    print(ggplot(q[nU >= MIN_N], aes(x_a1_rd / nB, x_a1_py / nU)) +
            geom_abline(slope = 1, intercept = 0, colour = "#c0392b") +
            geom_point(alpha = 0.25, size = 0.6) +
            coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) + facet_wrap(~sample) +
            labs(title = "The intended difference: reads vs molecules, pysam only",
                 subtitle = "Off the red line by construction - one molecule, one vote",
                 x = "Bl6 fraction per read", y = "Bl6 fraction per molecule") +
            theme_bw(base_size = 10))
  }
  invisible(dev.off())
  msg("Wrote %s", OUT_PDF)
}
