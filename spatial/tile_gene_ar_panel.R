# Per-tile allelic ratio for a NAMED PANEL of genes, adult (9w) beside aged (78w).
#
# escape_genes.R screens chrX and plots whatever wins. This plots a fixed list,
# including AUTOSOMAL imprinted loci, which the screen never looks at - and the
# imprinted block is the point: it is the only control that can tell a genuine
# CAST signal on chrX from a B6 molecule miscalled as CAST.
#
# Run on the cluster:
#   conda activate seurat_env
#   Rscript ~/Postdoc/spatial/tile_gene_ar_panel.R
#   IN_ROOT=.../ase_pysam_64um OUT_PREFIX=tile_gene_ar_umi Rscript ...
# or from slurm/spatial_tile_gene_ar.slurm.
#
# INPUT is the tile x gene table the pysam counter writes:
#   <IN_ROOT>/<sample>/tile_gene_counts.tsv.gz
# with a1/a2 at BOTH levels. Default IN_ROOT is the duplicate-inclusive 64 um
# tree (ase_pysam_dup_64um), where a1_reads/a2_reads keep PCR duplicates and
# a1_umi/a2_umi do not - the same molecules, counted twice over.
#
# GENOTYPE, which sets the direction of every number here (same as
# escape_genes.R and tile_ratio_map.R): B6 mother x CAST father, and B6 carries
# the Xist deletion, so the CAST X is the inactive X in EVERY cell.
#   a1 = B6 = REF = maternal = the active X      AR = a1 / (a1 + a2)
#   a2 = CAST = ALT = paternal = the inactive X
# On chrX:  AR ~ 1 is fully silenced, and escape = 1 - AR.
# On the imprinted loci: AR ~ 1 = maternally expressed, AR ~ 0 = paternally
# expressed. Those are the two EXPECTED values, drawn per gene as a short bar,
# so a locus landing away from its expectation is visible without arithmetic.
#
# WHY THE IMPRINTED BLOCK IS WORTH MORE THAN IT LOOKS. Reference mapping bias
# runs B6-ward (autosomal AR minus 0.5: +0.028 at 9w, +0.043 at 78w), so it
# pushes chrX towards "silenced" and makes escape an UNDER-estimate. The only
# error that could FAKE escape is a B6 molecule called CAST, and a maternally
# expressed imprinted gene - H19, Meg3, Rian, Cdkn1c, Igf2r, Grb10, Zim1 - is
# pure B6 by construction. Whatever CAST fraction they show is the false-escape
# floor, measured rather than assumed. The paternal loci (AR ~ 0) test the
# mirror error and also confirm the cross direction.
#
# TWO THINGS IN THE REQUESTED PANEL ARE NOT MEASURABLE, and the script says so
# on the figure rather than quietly dropping them:
#
#   Xist. The SNP mask this tree was counted against is
#   SNPfile_..._mm39_no_Xist.bed - Xist SNPs are deliberately absent, so Xist
#   has no informative positions and never appears in the gene table. That is
#   not a fixable oversight here: B6 carries the Xist DELETION, so every Xist
#   read is CAST by construction and its AR would be ~0 whatever the skew is.
#   Xist cannot report skew in this genotype. The skew panel therefore uses
#   chrX POOLED PER TILE, read straight off tile_chrom_counts.tsv, which is the
#   quantity Xist would have been a proxy for.
#
#   Mcts2, Snurf, Snrpn, Kcnq1ot1, Zrsr1, Peg13, 2210013O21Rik, 4930578C19Rik.
#   Absent from the RefSeq annotation the counter used (annotation_us.bed) under
#   those names, so no read was ever assigned to them. Kcnq1 (the sense gene) and
#   Zrsr2 (chrX) ARE present and are added beside their missing partners, marked
#   as substitutes. Snrpn is the one real loss: it is the cleanest paternal
#   control there is, and the cross direction was originally fixed on it.
#
# WHY DUPLICATE-INCLUSIVE READS NEED A CORRECTED INTERVAL. Keeping duplicates
# multiplies the count without adding an observation - roughly 12 reads per
# distinct molecule on this data. A Wilson interval on n_reads would therefore
# be about 3.5x too narrow. The pooled interval on the read page is computed at
# n_eff = n_umi, the molecules those reads collapse onto, taken per gene per
# sample from the same table. The POINT stays the duplication-weighted read
# ratio, because that is the statistic the dup tree exists to show; only its
# precision is honest about how many independent molecules are behind it.
#
# n = 1 animal per age. 9w beside 78w is a description of two sections, never an
# age effect, and nothing here should be tested as one.

.libPaths(c("~/R/matrix-dev", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

##### ---------------------- CONFIG ---------------------- #####
BASE    <- "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/adult_aged_spatial"
IN_ROOT <- Sys.getenv("IN_ROOT", file.path(BASE, "ase_pysam_dup_64um"))
SAMPLES <- strsplit(Sys.getenv("SAMPLES", "9w,78w"), ",")[[1]]
TILE_UM <- as.integer(Sys.getenv("TILE_UM", "64"))
OUT_PREFIX <- Sys.getenv("OUT_PREFIX", "tile_gene_ar_panel")
OUT_DIR    <- Sys.getenv("OUT_DIR", IN_ROOT)

# Which levels get a page. "reads" is the duplicate-inclusive statistic in the
# dup tree; "umi" is the same molecules deduplicated, from the same file, and is
# drawn as a second page because the two disagree on this data and the disagreement
# is itself worth seeing (chrX pooled 0.898 per molecule vs 0.921 per read at 9w).
LEVELS <- strsplit(Sys.getenv("LEVELS", "reads,umi"), ",")[[1]]

# A tile earns a dot at this many units. Set per level, because a "read" in the
# dup tree is worth ~1/12 of a molecule: 24 dup reads and 2 molecules are
# roughly the same amount of evidence, and one threshold for both would make the
# read page look twelve times better informed than it is.
MIN_TILE_N <- c(reads = 24L, umi = 2L)
MIN_TILES_VIOLIN <- 20L    # fewer tiles than this: dots only, gene marked *

# A violin also needs the tiles to be able to TAKE intermediate values. Outside
# Smpx, a 64 um tile carries about one informative molecule of any gene in this
# panel: at 9w, >= 2 UMIs is reached by 167/1047 Tspan7 tiles, 10/228 Kdm5c
# tiles and 0/44 Ftx tiles, and 97-100% of tiles read exactly 0 or 1 even on
# duplicate-inclusive reads - because every duplicate of the one molecule votes
# the same way. A violin over stacks at 0 and 1 draws an hourglass and invites
# it to be read as bimodal biology; it is arithmetic. So a gene earns a violin
# only if the MEDIAN drawn tile has more than one molecule, and everything else
# gets the composition page instead, which says what the two stacks actually are.
MIN_TILE_MEDIAN_UMI <- 2
##### ------------------------------------------------------ #####

msg <- function(...) message(sprintf(...))

# The gene panel lives in one file, shared with tile_gene_ar_maps.R, so the
# distributions and the maps cannot end up showing different gene sets. Resolved
# from this script's own directory first, so a checkout anywhere works.
source_panel <- function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else NULL
  for (d in c(here, "spatial", "~/Postdoc/spatial", ".")) {
    f <- file.path(d, "tile_gene_panel.R")
    if (file.exists(f)) { source(f, local = FALSE); return(invisible(f)) }
  }
  stop("Cannot find tile_gene_panel.R next to this script")
}
msg("panel from %s", source_panel())

##### ----------------------- LOAD ----------------------- #####
want <- PANEL$gene
read_genes <- function(s) {
  f <- file.path(IN_ROOT, s, "tile_gene_counts.tsv")
  if (!file.exists(f) && file.exists(paste0(f, ".gz"))) f <- paste0(f, ".gz")
  if (!file.exists(f)) stop("No tile_gene_counts.tsv[.gz] under ", file.path(IN_ROOT, s))
  fread(f, select = c("tile", "trow", "tcol", "chrom", "gene",
                      "a1_umi", "a2_umi", "n_umi", "a1_reads", "a2_reads", "n_reads"),
        showProgress = FALSE)[, sample := s][]
}
read_chrom <- function(s) {
  f <- file.path(IN_ROOT, s, "tile_chrom_counts.tsv")
  if (!file.exists(f)) stop("No tile_chrom_counts.tsv under ", file.path(IN_ROOT, s))
  fread(f, showProgress = FALSE)[, sample := s][]
}

msg("Reading %s", IN_ROOT)
all_g <- rbindlist(lapply(SAMPLES, read_genes))
all_c <- rbindlist(lapply(SAMPLES, read_chrom))

SLAB <- c("9w" = "adult (9w)", "78w" = "aged (78w)")
lab_sample <- function(dt) {
  dt[, sample := factor(SLAB[sample], levels = unname(SLAB[SAMPLES]))][]
}
lab_sample(all_g); lab_sample(all_c)

# The autosomal ratio, per sample and per level: the biallelic reference and the
# mapping-bias estimate in one number. From the CHROMOSOME table, so it is every
# informative unit on chr1-19 and not just the ones an annotated gene claimed.
AUTOS <- paste0("chr", 1:19)
base <- all_c[chrom %chin% AUTOS, .(
  auto_ar_umi   = sum(a1_umi) / sum(n_umi),
  auto_ar_reads = sum(a1_reads) / sum(n_reads),
  dup_factor    = sum(n_reads) / sum(n_umi)), by = sample]
cat("\nautosomal ratio (= what full escape looks like) and duplication:\n"); print(base)

# What is actually present. Reported before anything is plotted, because a gene
# silently missing from a control panel is the failure mode that matters here.
present <- all_g[, unique(gene)]
missing <- setdiff(want, c(present, "chrX (all genes pooled)"))
if (length(missing)) {
  cat("\nNOT IN THE GENE TABLE (no informative SNP under this annotation/mask):\n  ",
      paste(missing, collapse = ", "), "\n")
  cat("  Xist is expected here - the mask is SNPfile_..._no_Xist.bed, and B6 carries\n",
      "  the Xist deletion, so Xist could not report skew even with SNPs. The Skew\n",
      "  panel uses chrX pooled per tile instead.\n")
}
sub_ok <- SUBS[gene %chin% present]
if (nrow(sub_ok)) {
  cat("\nsubstitutes plotted for missing panel members:\n")
  print(sub_ok[, .(plotted = gene, in_place_of = for_gene)])
}

##### -------------------- ASSEMBLE ---------------------- #####
# chrX pooled per tile, given a gene's shape so it flows through the same code.
xtile <- all_c[chrom == "chrX", .(sample, tile, trow, tcol, chrom,
                                  gene = "chrX (all genes pooled)",
                                  a1_umi, a2_umi, n_umi, a1_reads, a2_reads, n_reads)]

keep <- c(want, sub_ok$gene)
dat <- rbind(all_g[gene %chin% keep,
                   .(sample, tile, trow, tcol, chrom, gene, a1_umi, a2_umi, n_umi,
                     a1_reads, a2_reads, n_reads)],
             xtile, use.names = TRUE)

# One row per gene per sample, pooled over EVERY tile including the thin ones.
pool <- dat[, .(a1_umi = sum(a1_umi), a2_umi = sum(a2_umi), n_umi = sum(n_umi),
                a1_reads = sum(a1_reads), a2_reads = sum(a2_reads),
                n_reads = sum(n_reads), tiles = .N,
                med_tile_umi = as.numeric(median(n_umi)),
                # Tiles that can only say "B6" or "CAST" because they hold one
                # molecule. Near 100% here for everything but Smpx, and the
                # single most important number for reading the dots.
                frac_tile_mono = sum(n_umi > 0 & (a1_umi == 0 | a2_umi == 0)) /
                                 pmax(sum(n_umi > 0), 1L)),
            by = .(sample, gene, chrom)]
# n_umi == 0 with n_reads > 0 is a REAL row, not a bug: decide() drops a
# molecule whose reads split evenly between the alleles, while those reads each
# still carried their own majority into the read columns. Rare (29 rows in this
# panel) but it divides by zero if left alone, so every molecule-level quantity
# is NA there and the read-level interval falls back to n_eff = 1.
pool[, `:=`(ar_umi     = fifelse(n_umi > 0, a1_umi / n_umi, NA_real_),
            ar_reads   = fifelse(n_reads > 0, a1_reads / n_reads, NA_real_),
            dup_factor = fifelse(n_umi > 0, n_reads / n_umi, NA_real_))]

# Wilson, with n given separately from k/n so the read-level interval can be
# widened to the number of MOLECULES behind those reads. Duplicates inflate k
# and n together and would otherwise buy precision that does not exist.
wilson <- function(p, n, z = 1.96) {
  n <- pmax(n, 1e-9)
  d  <- 1 + z^2 / n
  c1 <- (p + z^2 / (2 * n)) / d
  h  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(lo = pmax(0, c1 - h), hi = pmin(1, c1 + h))
}
pool[, c("lo_umi", "hi_umi")     := wilson(ar_umi, n_umi)]
pool[n_umi == 0, c("lo_umi", "hi_umi") := NA_real_]
pool[, c("lo_reads", "hi_reads") := wilson(ar_reads, pmax(n_umi, 1))]  # n_eff = molecules

meta <- PANEL_META(sub_ok)
pool <- merge(pool, meta, by = "gene", all.x = TRUE)
dat  <- merge(dat,  meta, by = "gene", all.x = TRUE)

ord <- PANEL_GENES(have = pool$gene, subs = sub_ok)

##### ---------------------- OUTPUT ---------------------- #####
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

csv <- pool[gene %chin% ord][order(group, match(gene, ord), sample),
  .(group, gene, chrom, sample, tiles,
    n_umi, ar_umi = round(ar_umi, 4),
    umi_lo = round(lo_umi, 4), umi_hi = round(hi_umi, 4),
    n_reads, ar_reads = round(ar_reads, 4),
    reads_lo = round(lo_reads, 4), reads_hi = round(hi_reads, 4),
    dup_factor = round(dup_factor, 2),
    med_tile_umi, frac_tile_mono = round(frac_tile_mono, 3),
    expect, escape_umi = round(1 - ar_umi, 4))]
f_csv <- file.path(OUT_DIR, paste0(OUT_PREFIX, "_pooled.csv"))
fwrite(csv, f_csv)
msg("Wrote %s", f_csv)

f_tiles <- file.path(OUT_DIR, paste0(OUT_PREFIX, "_per_tile.csv.gz"))
fwrite(dat[gene %chin% ord][order(group, match(gene, ord), sample, tile)], f_tiles)
msg("Wrote %s", f_tiles)

cat("\npooled allelic ratio (a1 = B6 = maternal = active X):\n")
print(csv[, .(group, gene, sample, tiles, n_umi, ar_umi, n_reads, ar_reads,
              med_tile_umi, frac_tile_mono, expect)])

# --- one page per level per group -------------------------------------------
draw <- function(level, grp) {
  a1 <- paste0("a1_", level); nn <- paste0("n_", level)
  min_n <- MIN_TILE_N[[level]]
  unit  <- if (level == "umi") "molecules" else "reads (PCR duplicates kept)"
  unit_s <- if (level == "umi") "molecules" else "dup reads"   # legend title

  gs <- ord[ord %chin% meta[group == grp, gene]]
  gs <- gs[gs %chin% pool[group == grp, gene]]
  if (!length(gs)) return(NULL)

  d <- dat[group == grp & gene %chin% gs & get(nn) >= min_n]
  d[, N := get(nn)]                    # aes() cannot use get(); name it here
  d[, AR := get(a1) / N]
  d[, gene := factor(gene, levels = gs)]

  pl <- pool[group == grp & gene %chin% gs]
  pl[, `:=`(gene = factor(gene, levels = gs),
            AR = get(paste0("ar_", level)),
            lo = get(paste0("lo_", level)), hi = get(paste0("hi_", level)))]
  pl <- pl[is.finite(AR)]

  # Enough tiles AND enough molecules per tile - see MIN_TILE_MEDIAN_UMI.
  nt  <- d[, .(N = .N, med = as.numeric(median(n_umi))), by = .(gene, sample)]
  fatg <- nt[, .(ok = all(N >= MIN_TILES_VIOLIN) && all(med >= MIN_TILE_MEDIAN_UMI) &&
                      uniqueN(sample) == length(SAMPLES)),
             by = gene][ok == TRUE, as.character(gene)]
  fat <- d[as.character(gene) %chin% fatg]

  # Genes with no plottable tile at all still get an axis slot and a diamond -
  # dropping them would hide exactly the loci this panel exists to check.
  labs_v <- setNames(paste0(gs, ifelse(gs %chin% fatg, "", " *"),
                            ifelse(gs %chin% sub_ok$gene,
                                   paste0(" (for ", sub_ok$for_gene[match(gs, sub_ok$gene)], ")"),
                                   "")), gs)

  exp_dt <- unique(pl[!is.na(expect), .(gene, expect)])

  ggplot(d, aes(gene, AR, fill = sample)) +
    geom_hline(yintercept = 1, linetype = 3, colour = "#8a8a8a") +
    geom_hline(yintercept = 0, linetype = 3, colour = "#8a8a8a") +
    geom_hline(yintercept = mean(base[[paste0("auto_ar_", level)]]),
               linetype = 2, colour = "#184f95") +
    annotate("text", x = 0.5, y = mean(base[[paste0("auto_ar_", level)]]),
             label = "autosomal ratio (= biallelic)", hjust = 0, vjust = -0.6,
             size = 2.4, colour = "#184f95") +
    { if (nrow(exp_dt))
        geom_crossbar(data = exp_dt, aes(x = gene, y = expect, ymin = expect, ymax = expect),
                      inherit.aes = FALSE, width = 0.8, colour = "#c97314",
                      linewidth = 0.4, fatten = 0) } +
    { if (nrow(fat)) geom_violin(data = fat, position = position_dodge(width = 0.8),
                                 width = 0.8, colour = NA, alpha = 0.4,
                                 scale = "width", trim = TRUE) } +
    geom_point(aes(size = N),
               position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.8),
               alpha = 0.35, shape = 16, stroke = 0) +
    geom_pointrange(data = pl, aes(gene, AR, ymin = lo, ymax = hi, group = sample),
                    inherit.aes = FALSE, position = position_dodge(width = 0.8),
                    shape = 23, size = 0.35, fill = "white", colour = "#0b0b0b",
                    stroke = 0.6, fatten = 3) +
    scale_fill_manual(values = setNames(c("#2D6E5D", "#b02a2a"), unname(SLAB[SAMPLES])),
                      name = NULL) +
    scale_size_continuous(range = c(0.4, 2.4), name = paste(unit_s, "in tile"),
                          trans = "sqrt") +
    scale_x_discrete(labels = labs_v, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = sprintf("%s: allelic ratio per %d um tile, %s", grp, TILE_UM, unit),
         subtitle = paste0(
           "AR = B6 / total. B6 is maternal AND the active X, so on chrX low AR = escape, ",
           "and on the imprinted loci AR ~ 1 = maternal, AR ~ 0 = paternal.\n",
           "White diamond = pooled ratio over ALL tiles incl. the undrawn ones; its bar is a 95% Wilson interval ",
           if (level == "reads") "at n_eff = MOLECULES, not duplicate reads.\n" else "on molecules.\n",
           sprintf("Dots are tiles with >= %d %s. * = no violin: too few tiles (< %d) or a median tile holding < %d molecules, so its ratio could only be 0 or 1. Read the diamond. ",
                   min_n, unit, MIN_TILES_VIOLIN, MIN_TILE_MEDIAN_UMI),
           "n = 1 animal per age: a description of two sections, not an age effect."),
         x = NULL, y = "allelic ratio (B6 / total)",
         caption = paste(
           "Dotted lines = the two monoallelic poles. Dashed = the autosomal ratio, which is what a fully biallelic locus reads once B6-ward mapping bias is allowed for.",
           if (nrow(exp_dt)) "\nOrange bar = the allele this locus is expected on from its imprinting status; distance from it is the false-call floor, measured rather than assumed." else "",
           if (level == "reads") "\nDuplicate-inclusive reads are not independent observations - the dots inherit that, the interval does not." else "")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.subtitle = element_text(size = 7.5, colour = "#52514e"),
          plot.caption  = element_text(size = 6.5, colour = "#52514e", hjust = 0))
}

# WHY THIS PAGE EXISTS. When almost every tile holds one molecule, "the per-tile
# allelic ratio" is not a distribution - it is a coin, and the only thing the
# tiles can report is HOW MANY of them came up CAST. That is worth seeing
# directly: it is the same information as the pooled diamond, but on the scale
# the question was asked on, and it is the only view in which a gene whose
# escape sits in a few tiles would look different from one whose escape is
# spread evenly. The two ages are drawn as stacked rows over a shared x axis, so
# a gene is read down the column.
draw_comp <- function(grp) {
  gs <- ord[ord %chin% meta[group == grp, gene]]
  gs <- gs[gs %chin% pool[group == grp, gene]]
  if (!length(gs)) return(NULL)

  d <- dat[group == grp & gene %chin% gs & n_umi > 0]
  if (!nrow(d)) return(NULL)
  d[, cls := fifelse(a2_umi == 0, "all B6 (AR = 1)",
              fifelse(a1_umi == 0, "all CAST (AR = 0)", "mixed (0 < AR < 1)"))]
  d[, cls := factor(cls, levels = c("all B6 (AR = 1)", "mixed (0 < AR < 1)",
                                    "all CAST (AR = 0)"))]
  d[, gene := factor(gene, levels = gs)]

  cmp <- d[, .N, by = .(sample, gene, cls)][, frac := N / sum(N), by = .(sample, gene)]
  tot <- d[, .(N = .N), by = .(sample, gene)]

  ggplot(cmp, aes(gene, frac, fill = cls)) +
    geom_col(width = 0.75) +
    geom_text(data = tot, aes(gene, 1.04, label = N), inherit.aes = FALSE,
              size = 2.1, colour = "#52514e") +
    facet_wrap(~ sample, ncol = 1) +
    scale_fill_manual(values = c("all B6 (AR = 1)" = "#b02a2a",
                                 "mixed (0 < AR < 1)" = "#C97314",
                                 "all CAST (AR = 0)" = "#184f95"), name = NULL) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%"),
                       breaks = seq(0, 1, 0.25), limits = c(0, 1.08)) +
    labs(title = sprintf("%s: what each %d um tile could say, %s",
                         grp, TILE_UM, "molecule level"),
         subtitle = paste0(
           "Every tile with at least one informative molecule, classified. Numbers above the bars are tile counts.\n",
           "Outside Smpx a tile holds about ONE molecule of these genes, so its ratio can only be 0 or 1 - the dot plots' two stacks are this bar, and the blue fraction is the escape estimate.\n",
           "n = 1 animal per age: a description of two sections, not an age effect."),
         x = NULL, y = "tiles",
         caption = "Blue = the only molecule in that tile came from the INACTIVE (CAST) X. On the imprinted controls it is the false-call rate; on chrX it is escape.") +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.subtitle = element_text(size = 7.5, colour = "#52514e"),
          plot.caption  = element_text(size = 6.5, colour = "#52514e", hjust = 0))
}

f_pdf <- file.path(OUT_DIR, paste0(OUT_PREFIX, ".pdf"))
pdf(f_pdf, width = 13, height = 6, onefile = TRUE)
pages <- 0L
for (lv in LEVELS) for (grp in levels(PANEL$group)) {
  p <- draw(lv, grp)
  if (!is.null(p)) { print(p); pages <- pages + 1L }
}
for (grp in levels(PANEL$group)) {
  p <- draw_comp(grp)
  if (!is.null(p)) { print(p); pages <- pages + 1L }
}
invisible(dev.off())
msg("Wrote %s (%d pages: %s x %s)", f_pdf, pages,
    paste(LEVELS, collapse = "/"), paste(levels(PANEL$group), collapse = "/"))
