# Build every SNP and interval bed the allelic analyses consume, from one
# annotation and one SNP file, so the coordinate conventions are handled in one
# place instead of per script.
#
# Outputs (all written to OUT_DIR):
#
#   SNP beds - subsets of the genome-wide SNP file, consumed by Allelome.PRO2
#   and by spatial/ase_bin_allele_counts.py:
#     SNPfile_..._mm39_no_Xist.bed          Xist gene span dropped (the live one)
#     SNPfile_..._mm39_no_Xic500kb.bed      Xist +- XIC_FLANK dropped
#     SNPfile_..._mm39_core_escape_new.bed  the 11-gene core escape set only
#
#   Interval beds - consumed by ase_bin_allele_counts.py --subset-bed, which
#   counts UMIs inside (or outside) them as extra columns of ONE counting pass:
#     annotation_us_mm39_chrX.bed           all chrX genes
#     annotation_us_mm39_core_escape.bed    core escape, one row per gene
#     annotation_us_mm39_Xic500kb.bed       the Xic mask interval itself
#     annotation_us_mm39_imprinted_pat.bed  paternally expressed imprinted loci
#     annotation_us_mm39_imprinted_mat.bed  maternally expressed imprinted loci
#
# Also writes bed_provenance.tsv: every output, its md5, and the inputs it came
# from. Three "no-Xist" builds are named in this repo and only one is live, so a
# result cannot otherwise be attributed to the bed that produced it.
#
# ---------------------------------------------------------------------------
# COORDINATE CONVENTIONS - the one thing to get right here
#
# annotation_us.bed is a real BED: V2 is a 0-based start, V3 an exclusive end,
# so the interval covers 0-based positions [V2, V3).
#
# SNPfile_C57BL6_NJ_CAST_EiJ.bed is NOT. Its column V2 is the 1-based VCF
# POSITION and V3 is V2 + 1 - the `pos, pos+1` form. That is established
# independently by the spatial counter's offset probe, which scored
# --snp-offset -1 at 100.0% against +0 at 49.0% on 200k reads
# (slurm/spatial_ase_sweep.slurm), i.e. the 0-based coordinate pysam reports is
# one LESS than column V2.
#
# So a 1-based SNP position p lies inside a BED interval [V2, V3) iff
#   V2 <= p - 1 < V3   <=>   p in (V2, V3]   <=>   p in seq(V2 + 1, V3)
#
# The previous version of this script used seq(V2, V3), which is one base too
# wide at the LEFT edge and correct at the right. It only ever mattered for a
# single SNP per interval boundary, but "documented or fixed, not ambiguous" is
# the standard here, so it is fixed and the offset is named in one constant.
# ---------------------------------------------------------------------------

library(tidyverse)
library(data.table)

REF_BASE   <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/Y_references/mm39"
ANNOT_FILE <- file.path(REF_BASE, "20250512_RefSeq", "annotation_us.bed")
SNP_FILE   <- file.path(REF_BASE, "SNPs", "SNPfile_C57BL6_NJ_CAST_EiJ.bed")
OUT_DIR    <- Sys.getenv("OUT_DIR", ".")

# Half-open BED [start, end) holds the 1-based SNP positions start+1 .. end.
SNP_TO_BED_SHIFT <- 1L

# How far either side of Xist to mask. 500 kb reaches Ftx, Jpx, Tsix and the
# intergenic Xic, all of which report the inactive X in this design and are
# therefore effectively CAST-only - exactly the reads a "no-Xist" mask leaves
# behind. In this annotation Xist is chrX:102,503,978-102,526,839 (mm39, 22.9 kb),
# so the mask covers chrX:102,003,978-103,026,839. The corresponding mm10 locus
# is at 103.46 Mb - a ~956 kb shift, which is the check that the annotation and
# the BAM are on the same assembly.
XIC_FLANK <- 500000L

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(OUT_DIR, f)

annotation <- fread(ANNOT_FILE)
SNPfile    <- fread(SNP_FILE)

# 1-based SNP positions covered by a set of BED intervals. unlist(Map(seq, ...))
# rather than seq() on the columns directly: seq() is not vectorised, so the
# direct form silently uses only the first row when the annotation has one
# interval per TRANSCRIPT - which it does - and errors outright when it has
# more than one.
bed_positions <- function(dt) {
  if (!nrow(dt)) return(integer(0))
  unique(unlist(Map(seq, dt$V2 + SNP_TO_BED_SHIFT, dt$V3), use.names = FALSE))
}

# Collapse transcript rows to one interval per gene.
gene_level <- function(dt) {
  dt %>%
    dplyr::group_by(V1, V4, V6) %>%
    dplyr::summarise(V2 = min(V2), V3 = max(V3), V5 = 0L, .groups = "drop") %>%
    dplyr::select(V1, V2, V3, V4, V5, V6) %>%
    dplyr::arrange(V1, V2)
}

written <- list()
write_bed <- function(x, file, what) {
  write.table(x, out(file), row.names = FALSE, col.names = FALSE,
              quote = FALSE, sep = "\t")
  written[[file]] <<- data.table(file = file, rows = nrow(x), what = what)
  message(sprintf("  %-46s %8d rows   %s", file, nrow(x), what))
}

# create annotation summed up for all chromosome
annotation_summed <- annotation %>%
  group_by(V1) %>%
  summarise(V2 = min(V2), V3 = max(V3), .groups = "drop")

##### ---------------------- Xist / Xic masks ---------------------- #####

Xist_annotation <- annotation[V4 == "Xist", ]
if (!nrow(Xist_annotation)) stop("No Xist row in ", ANNOT_FILE)
message(sprintf("Xist: %d annotation rows, chrX:%d-%d (BED, 0-based)",
                nrow(Xist_annotation), min(Xist_annotation$V2),
                max(Xist_annotation$V3)))

Xist_span <- bed_positions(Xist_annotation)
SNPfile_no_Xist <- SNPfile %>%
  dplyr::filter(!(V1 == "chrX" & V2 %in% Xist_span))
write_bed(SNPfile_no_Xist, "SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xist.bed",
          sprintf("Xist span masked (%d SNPs dropped)",
                  nrow(SNPfile) - nrow(SNPfile_no_Xist)))

# The Xic mask. One interval spanning Xist +- XIC_FLANK, clamped at 0.
Xic_interval <- data.table(
  V1 = "chrX",
  V2 = max(0L, min(Xist_annotation$V2) - XIC_FLANK),
  V3 = max(Xist_annotation$V3) + XIC_FLANK,
  V4 = sprintf("Xic_Xist_plusminus_%dkb", XIC_FLANK %/% 1000L),
  V5 = 0L,
  V6 = "+"
)
write_bed(Xic_interval, "annotation_us_mm39_Xic500kb.bed",
          sprintf("Xist +- %d kb, one interval", XIC_FLANK %/% 1000L))

# Which annotated genes the mask swallows. Printed, not just written: the whole
# argument for the mask is that Ftx, Jpx and Tsix are inside it, so it is worth
# seeing them named in the log.
in_mask <- annotation[V1 == "chrX" & V2 < Xic_interval$V3 & V3 > Xic_interval$V2]
message("  genes inside the Xic mask: ",
        paste(sort(unique(in_mask$V4)), collapse = ", "))

Xic_span <- bed_positions(Xic_interval)
SNPfile_no_Xic <- SNPfile %>%
  dplyr::filter(!(V1 == "chrX" & V2 %in% Xic_span))
write_bed(SNPfile_no_Xic,
          "SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_no_Xic500kb.bed",
          sprintf("Xist +- %d kb masked (%d SNPs dropped, %d more than the Xist-only mask)",
                  XIC_FLANK %/% 1000L,
                  nrow(SNPfile) - nrow(SNPfile_no_Xic),
                  nrow(SNPfile_no_Xist) - nrow(SNPfile_no_Xic)))

##### -------------------------- chrX genes -------------------------- #####

annotation_chrX <- annotation[V1 == "chrX", ]
write_bed(annotation_chrX, "annotation_us_mm39_chrX.bed",
          "all chrX transcripts")

##### ------------------------- core escape ------------------------- #####

core_escape_genes <- c("Kdm5c", "Kdm6a", "Ddx3x", "Eif2s3x")

new_core_escape_genes <- c("Kdm5c", "Kdm6a", "Ddx3x", "Eif2s3x", "Ftx",
                           "5530601H04Rik", "Jpx", "Pbdc1", "Utp14a",
                           "Akap17a", "Sts")

core_escape <- annotation[V4 %in% new_core_escape_genes, ]
missing <- setdiff(new_core_escape_genes, core_escape$V4)
if (length(missing)) message("  NOT in the annotation: ",
                             paste(missing, collapse = ", "))

core_escape_span <- bed_positions(core_escape)
SNPfile_core_escape <- SNPfile %>%
  dplyr::filter(V1 == "chrX" & V2 %in% core_escape_span)
write_bed(SNPfile_core_escape,
          "SNPfile_C57BL_6NJxCAST_EiJ_sorted_mm39_core_escape_new.bed",
          sprintf("%d core escape genes", uniqueN(core_escape$V4)))

core_escape_gene_level <- gene_level(core_escape)
write_bed(core_escape_gene_level, "annotation_us_mm39_core_escape.bed",
          "core escape, one interval per gene")

##### -------------------------- imprinted -------------------------- #####
#
# Positive control for per-UMI allele assignment, for spatial/ase_tile_sweep.R.
# An imprinted locus is monoallelic in every cell, so two UMIs from it should
# carry the SAME allele at any separation - which is the assumption the chrX
# pair correlation C(d) is being asked to test, and nothing in the pipeline
# currently establishes it. C_imprinted(4um) near 1 says the 12% disagreement
# on chrX is biology; C_imprinted(4um) also ~0.78 says it is the method's
# per-UMI error floor and the chrX number is measuring noise.
#
# SPLIT BY DIRECTION, and this is not optional. H19/Igf2, Airn/Igf2r and
# Meg3/Dlk1 are RECIPROCALLY imprinted pairs: within each pair one gene is
# expressed from the maternal chromosome and the other from the paternal one.
# Pooling them would put UMIs carrying opposite alleles in one set, and C(d)
# over that mixture would sit near 0.5 with perfectly clean data - looking
# exactly like a failed positive control. Each bed below is
# direction-homogeneous, so C(d) within it is a clean read on the error floor.
#
# Which strain is which parent does not matter for C(d) - only that every locus
# in a set agrees. It does show up in the per-locus ref fraction that
# ase_bin_allele_counts.py --locus-out writes: a maternally expressed set
# sitting near ref = 1 means B6 is the mother in this cross, near ref = 0 means
# CAST is. Read it off there rather than assuming.
#
# Imprinting is also tissue-dependent (Grb10 and Ube3a are the usual
# offenders), which is why the downstream script filters on the OBSERVED
# per-locus ratio as well as on depth rather than trusting this list.

imprinted_paternal <- c(
  # canonical paternally expressed
  "Igf2", "Airn", "Snrpn", "Peg3", "Dlk1", "Peg10", "Mest", "Nnat",
  "Plagl1", "Impact", "Nap1l5", "Magel2", "Ndn", "Mkrn3", "Slc38a4",
  "Kcnq1ot1", "Rtl1", "Usp29"
)
imprinted_maternal <- c(
  # canonical maternally expressed. Zim1 belongs HERE, not in the paternal set:
  # it is the maternally expressed reciprocal partner of Peg3 at the same locus.
  "H19", "Igf2r", "Meg3", "Rian", "Mirg", "Cdkn1c", "Zim1",
  "Phlda2", "Osbpl5", "Ascl2", "Slc22a18", "Zrsr1"
)

# Imprinted, but NOT usable in a positive control that has to be clean. Kept as
# a named list because they are still subtracted from both control beds below -
# dropping a gene from the control does not stop its reads landing inside a
# neighbour's collapsed interval.
#
#   Ube3a   maternal in brain, biallelic in most peripheral tissues
#   Gnas    direction depends on the transcript (Nesp maternal, Nespas paternal)
#   Grb10   maternal in most tissues, paternal in brain
#   Tssc4   weak and tissue-specific
imprinted_unreliable <- c("Ube3a", "Gnas", "Grb10", "Tssc4", "Nesp", "Nespas")
# The five loci the analysis plan names, kept as its own bed so the headline
# control can be run on textbook-solid ground before the wider list is used to
# buy depth.
imprinted_core <- c("H19", "Igf2", "Airn", "Igf2r", "Snrpn", "Peg3",
                    "Meg3", "Dlk1")

# Remove from `a` every base covered by `b`, keeping a's name on each surviving
# piece. An interval can be split in two, or vanish entirely.
#
# This is load-bearing, not tidying. gene_level() collapses transcript rows to
# one min-max interval per gene, and for genes with a long non-coding form that
# span is far larger than the coding gene: Snrpn collapses to 468 kb because of
# the Snhg14 transcript, which in mouse runs antisense ACROSS Ube3a. Ube3a is
# maternally expressed. Without this subtraction, 3377 of the paternal set's
# 4115 SNPs sit in that one interval, Ube3a's reads are counted as paternal, and
# the positive control drifts toward 0.5 - failing for a reason that has nothing
# to do with the per-UMI error rate it is meant to measure.
subtract_intervals <- function(a, b) {
  if (!nrow(a) || !nrow(b)) return(a)
  out <- list()
  for (i in seq_len(nrow(a))) {
    lo <- a$V2[i]; hi <- a$V3[i]
    cuts <- b[b$V1 == a$V1[i] & b$V3 > lo & b$V2 < hi, ]
    cuts <- cuts[order(cuts$V2), ]
    pieces <- list(c(lo, hi))
    for (j in seq_len(nrow(cuts))) {
      nxt <- list()
      for (pc in pieces) {
        if (cuts$V3[j] <= pc[1] || cuts$V2[j] >= pc[2]) { nxt[[length(nxt) + 1]] <- pc; next }
        if (cuts$V2[j] > pc[1]) nxt[[length(nxt) + 1]] <- c(pc[1], cuts$V2[j])
        if (cuts$V3[j] < pc[2]) nxt[[length(nxt) + 1]] <- c(cuts$V3[j], pc[2])
      }
      pieces <- nxt
    }
    for (pc in pieces) {
      out[[length(out) + 1]] <- data.table(V1 = a$V1[i], V2 = pc[1], V3 = pc[2],
                                           V4 = a$V4[i], V5 = 0L, V6 = a$V6[i])
    }
  }
  if (!length(out)) return(a[0, ])
  rbindlist(out)[order(V1, V2)]
}

# Intervals of every imprinted gene, whatever direction and whatever its
# reliability. Whichever control bed is being built, everything in here that is
# not part of that set gets cut out of it.
imprinted_all <- unique(c(imprinted_paternal, imprinted_maternal,
                          imprinted_unreliable))
imprinted_annot <- annotation[V4 %in% imprinted_all & V1 != "chrX" & V1 != "chrY", ]
imprinted_gl <- gene_level(imprinted_annot)

imprinted_bed <- function(genes, label, file) {
  keep <- imprinted_gl[imprinted_gl$V4 %in% genes, ]
  miss <- setdiff(genes, keep$V4)
  if (length(miss)) message("  ", label, ": not in the annotation - ",
                            paste(miss, collapse = ", "))
  if (!nrow(keep)) { message("  ", label, ": nothing to write"); return(keep) }
  drop <- imprinted_gl[!imprinted_gl$V4 %in% genes, ]
  g <- subtract_intervals(as.data.table(keep), as.data.table(drop))
  # Report per gene what the subtraction cost, so a control that has been eaten
  # by its neighbours is visible here and not inferred from a flat curve later.
  before <- as.data.table(keep)[, .(kb0 = sum(V3 - V2) / 1000), by = V4]
  after  <- g[, .(kb1 = sum(V3 - V2) / 1000, pieces = .N), by = V4]
  cmp <- merge(before, after, by = "V4", all.x = TRUE)
  cmp[is.na(kb1), `:=`(kb1 = 0, pieces = 0L)]
  cut <- cmp[kb1 < kb0 - 0.001][order(kb1 - kb0)]
  if (nrow(cut)) {
    message("  ", label, ": trimmed where another imprinted gene overlaps -")
    for (i in seq_len(nrow(cut))) {
      message(sprintf("    %-10s %8.1f -> %8.1f kb (%d piece%s)",
                      cut$V4[i], cut$kb0[i], cut$kb1[i], cut$pieces[i],
                      if (cut$pieces[i] == 1L) "" else "s"))
    }
  }
  if (any(cmp$kb1 == 0)) {
    message("  ", label, ": entirely removed - ",
            paste(cmp[kb1 == 0, V4], collapse = ", "))
  }
  write_bed(g, file, sprintf("%s, %d of %d genes, %d intervals after removing overlaps",
                             label, uniqueN(g$V4), length(genes), nrow(g)))
  g
}

imp_pat <- imprinted_bed(imprinted_paternal, "paternally expressed",
                         "annotation_us_mm39_imprinted_pat.bed")
imp_mat <- imprinted_bed(imprinted_maternal, "maternally expressed",
                         "annotation_us_mm39_imprinted_mat.bed")
imprinted_bed(intersect(imprinted_paternal, imprinted_core),
              "core, paternally expressed",
              "annotation_us_mm39_imprinted_core_pat.bed")
imprinted_bed(intersect(imprinted_maternal, imprinted_core),
              "core, maternally expressed",
              "annotation_us_mm39_imprinted_core_mat.bed")

# Chromosomes the imprinted beds reach. The spatial counter fetches these
# interval by interval rather than whole, so the control costs seconds - but it
# only knows to do that because the bed names them, so print them here.
message("  imprinted chromosomes: ",
        paste(sort(unique(c(imp_pat$V1, imp_mat$V1))), collapse = ","))

# How much SNP density each imprinted locus has. A locus with no informative SNP
# cannot contribute to the control however well it is expressed, and finding
# that out here is free.
# One row per GENE, summed over its pieces: after the overlap subtraction a gene
# can be several intervals, and a per-interval table would double-count nothing
# but would also hide how much SNP density the gene has left.
snp_density <- function(g, label) {
  if (!nrow(g)) return(NULL)
  d <- rbindlist(lapply(seq_len(nrow(g)), function(i) {
    n <- SNPfile[V1 == g$V1[i] & V2 > g$V2[i] & V2 <= g$V3[i], .N]
    data.table(set = label, gene = g$V4[i], chr = g$V1[i],
               kb = (g$V3[i] - g$V2[i]) / 1000, snps = n)
  }))
  d <- d[, .(kb = round(sum(kb)), snps = sum(snps), pieces = .N),
         by = .(set, gene, chr)]
  d[order(-snps)]
}
dens <- rbindlist(list(snp_density(imp_pat, "paternal"),
                       snp_density(imp_mat, "maternal")))
if (nrow(dens)) {
  fwrite(dens, out("imprinted_snp_density.tsv"), sep = "\t")
  message("\nInformative SNPs per imprinted locus (before any expression):")
  print(dens)
  if (any(dens$snps == 0)) {
    message("  loci with no SNP at all: ",
            paste(dens[snps == 0, gene], collapse = ", "))
  }
}

##### ------------------------- provenance ------------------------- #####

prov <- rbindlist(written)
prov[, md5 := vapply(file, function(f) unname(tools::md5sum(out(f))), "")]
prov <- rbindlist(list(
  data.table(file = c(ANNOT_FILE, SNP_FILE),
             rows = c(nrow(annotation), nrow(SNPfile)),
             what = c("INPUT annotation", "INPUT SNP file"),
             md5 = vapply(c(ANNOT_FILE, SNP_FILE),
                          function(f) unname(tools::md5sum(f)), "")),
  prov
))
prov[, `:=`(built = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            snp_to_bed_shift = SNP_TO_BED_SHIFT,
            xic_flank = XIC_FLANK)]
fwrite(prov, out("bed_provenance.tsv"), sep = "\t")
message("\nProvenance in ", out("bed_provenance.tsv"))
