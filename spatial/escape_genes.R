# Identify chrX genes that escape X-inactivation, and show the per-tile spread.
#
# Interactive. On the cluster:
#   conda activate seurat_env
#   cd $BASE/adult_aged_spatial/ase_pysam_64um
#   R
#   source("~/Postdoc/spatial/escape_genes.R")
# Everything worth changing is in the config block.
#
# GENOTYPE, which sets the direction of every number here: B6 carries the Xist
# deletion, so B6 cannot be inactivated and the CAST X is the inactive X in
# EVERY cell. a1 = B6 = active X, a2 = CAST = inactive X. Escape is therefore
#   escape = a2 / (a1 + a2) = 1 - AR
# A fully silenced gene sits at AR ~ 1. A fully escaping gene sits at the
# AUTOSOMAL ratio (~0.52), not at 0.50, because reference mapping is B6-biased -
# which makes escape under-estimated, never over-estimated.
#
# WHY GENES ARE POOLED ACROSS TILES BEFORE ANYTHING IS CALLED. A gene inside one
# 64 um tile carries a handful of molecules: only ~2% of chrX gene x tile rows
# reach 6 UMIs. Filtering per tile therefore selects the most highly expressed
# genes and nothing else - Reps2 is 0.1% escape over 7,154 UMIs and still throws
# tiles reading 5/6, because at n = 6 "AR < 0.9" only means "at least one CAST
# molecule". Meanwhile the real escapees are invisible to it: Kdm6a has 59
# informative UMIs in the whole section, Ftx 36, Ddx3x 12. Calls come from the
# pooled counts; tiles are only ever used for the picture.
#
# AND THE PICTURE HAS A LIMIT WORTH KNOWING BEFORE YOU READ IT. Escape genes are
# lowly expressed, so most of them have almost no per-tile data: at >= 2 UMIs a
# tile, Ndufb11 has 3346 tiles and Kdm5c has 7. A violin over 7 points is a
# drawing, not a distribution, so genes under MIN_TILES_VIOLIN get their dots
# and their pooled estimate but no violin, and are marked * on the axis. The
# pooled point and its interval are the estimate for every gene either way.

.libPaths(c("~/R/matrix-dev", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

##### ---------------------- CONFIG ---------------------- #####
DIR    <- "."            # directory holding 9w/ and 78w/
N_TOP  <- 15             # how many genes to plot
OUT    <- "escape_genes.pdf"

# Pooled informative UMIs per gene per sample before a gene is testable. 30 is
# deliberate, not lax: at 100 the list loses Pbdc1 and Usp9x, both real, and
# BH already controls the error rate over however many genes clear the bar
# (163 genes at 30, 86 at 100). Raise it if you want a shorter, safer list.
MIN_UMI          <- 30
MIN_TILE_UMI     <- 2    # per-tile UMIs before a tile earns a dot
MIN_TILES_VIOLIN <- 20   # tiles needed before a violin is drawn at all
FDR_CUT          <- 0.05

# ABOVE-BIALLELIC GENES. Escape has a ceiling: the Xi allele switching fully on
# gives EQUAL expression from both alleles, which is 0.5 - or the autosomal ratio
# 0.48 once B6-ward mapping bias is allowed for. So the entire escape range is
# ~0.10-0.50, and this threshold sits above it with headroom. It does not filter
# escape; it catches genes where the INACTIVE X supplies most of the transcript,
# which requires the ACTIVE (B6) allele to be silent, deleted or unmappable.
#
# Seven genes do that here, and they are not scattered - which is why they are
# worth looking at rather than deleting:
#   Slc16a2  102.74 Mb  98/98%  |  210-240 kb distal to Xist (102.50-102.53).
#   Gm53059  102.74 Mb  85/96%  |  The B6 chromosome is the one carrying the
#                                  Xist deletion, so a cis effect of the
#                                  deletion allele on its immediate neighbours
#                                  would look exactly like this.
#   Llph-ps2  13.08 Mb  89/94%  |  processed pseudogene lying between Usp9x
#                                  (12.94) and Ddx3x (13.15) - multimapping.
#   Fmr1      67.72 Mb  88/90%  |  neighbours, 680 kb apart, both extreme:
#   Aff2      68.40 Mb  98/97%  |  something regional.
#   Gm14719   48.75 Mb  64/68%  |  predicted gene.
#   Mid1     168.50 Mb  61/83%  |  PAR, but only 33/35 UMIs - wide interval.
#
# One reading would be real: a gene that escapes AND carries a strong cis-eQTL
# favouring CAST. That cannot be separated from an artefact without a reciprocal
# cross, so the default is to hold them out - but they are always PRINTED, and
# KEEP_ABOVE_BIALLELIC brings them into the ranking and the plot if you want to
# look at them directly. They are marked ! on the axis when kept.
MAX_ESCAPE           <- 0.60
KEEP_ABOVE_BIALLELIC <- TRUE

# EFFECT-SIZE FLOOR, on the AR scale, because that is the scale the criterion is
# usually stated on: a gene must reach AR < MAX_AR in at least one sample as
# well as passing FDR. AR = B6 = active X, so AR < 0.90 is "at least 10% of
# molecules come from the inactive X", the working definition of escape.
#
# Significance alone is not enough on its own here. Ndufb11 has 15,880
# informative UMIs, so its 7.4% escape is significant against the 6.1%
# chromosome rate while being nowhere near escape; the floor drops it. This is
# the same effect-size-versus-p-value problem as MIN_ENRICH in tile_ratio_map.R.
MAX_AR               <- 0.90

##### ----------------------- LOAD ----------------------- #####
read_one <- function(s) {
  f <- file.path(DIR, s, "tile_gene_counts.tsv")
  if (!file.exists(f) && file.exists(paste0(f, ".gz"))) f <- paste0(f, ".gz")
  fread(f, select = c("tile", "chrom", "gene", "a1_umi", "a2_umi", "n_umi"))[, sample := s][]
}
all <- rbind(read_one("9w"), read_one("78w"))
all[, sample := factor(sample, levels = c("9w", "78w"),
                       labels = c("adult (9w)", "aged (78w)"))]

# The autosomal ratio is both the biallelic reference and the mapping-bias
# estimate. Per sample, because the two libraries do not share it.
base <- all[chrom %chin% paste0("chr", 1:19),
            .(auto_ar = sum(a1_umi) / sum(n_umi)), by = sample]
cat("autosomal ratio (= what full escape looks like):\n"); print(base)

x <- all[chrom == "chrX"]

##### ------------------- POOLED CALLS ------------------- #####
g <- x[, .(a1 = sum(a1_umi), a2 = sum(a2_umi), n = sum(n_umi), tiles = .N),
       by = .(sample, gene)]
g[, escape := a2 / n]

# The null is the chromosome's own gene-assigned escape rate, not zero. Zero is
# the wrong null - there is a false-escape floor (<=1.5%, from the imprinted
# controls), so testing against it would call most of chrX significant. Testing
# against the chrX rate asks the question actually being asked: does this gene
# escape MORE than chrX does on average.
g[, chrx_rate := sum(a2) / sum(n), by = sample]
cat("\nnull rate (chrX gene-assigned escape):\n")
print(unique(g[, .(sample, chrx_rate = round(chrx_rate, 4))]))

g[, p := NA_real_]          # so the column exists even if nothing clears MIN_UMI
g[n >= MIN_UMI,
  p := mapply(function(k, nn, r) binom.test(k, nn, r, alternative = "greater")$p.value,
              a2, n, chrx_rate)]
g[, fdr := p.adjust(p, "BH"), by = sample]

# Wilson interval, so a gene with 30 UMIs is not shown as if it had 15,000.
wilson <- function(k, n, z = 1.96) {
  p <- k / n; d <- 1 + z^2 / n
  c1 <- (p + z^2 / (2 * n)) / d
  h  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(lo = c1 - h, hi = c1 + h)
}
g[, c("esc_lo", "esc_hi") := wilson(a2, n)]

# Testable in BOTH samples, so adult and aged are comparable gene by gene.
ok <- g[n >= MIN_UMI, .N, by = gene][N == 2L, gene]

above_bi <- g[gene %chin% ok][, .(mx = max(escape)), by = gene][mx > MAX_ESCAPE, gene]
if (length(above_bi)) {
  cat("\n--- ABOVE BIALLELIC (escape >", MAX_ESCAPE, "): the INACTIVE X supplies most of",
      "the transcript,\n    which needs the ACTIVE B6 allele silent, deleted or unmappable.",
      if (KEEP_ABOVE_BIALLELIC) "KEPT (marked ! below).\n" else "HELD OUT.\n")
  print(g[gene %chin% above_bi, .(gene, sample, n, escape = round(escape, 3))][order(gene, sample)])
}

cand <- if (KEEP_ABOVE_BIALLELIC) g[gene %chin% ok] else g[gene %chin% setdiff(ok, above_bi)]
# Significant in at least one animal. n = 1 per age, so this is a screen for
# genes worth looking at - NOT an age comparison. Nothing here supports a
# statement about ageing without replication.
# Both conditions: significant AND big enough to be escape. Either alone gives
# the wrong list - FDR alone lets Ndufb11 in on 15,880 UMIs at 7.4%, and AR
# alone lets in any thinly covered gene that happened to draw a few CAST UMIs.
sig  <- cand[fdr < FDR_CUT, unique(gene)]
big  <- cand[(1 - escape) < MAX_AR, unique(gene)]
hit  <- intersect(sig, big)
rank <- cand[gene %chin% hit, .(m = mean(escape)), by = gene][order(-m)]
top  <- head(rank$gene, N_TOP)

cat("\n---", length(sig), "pass FDR <", FDR_CUT, "|", length(big), "reach AR <", MAX_AR,
    "|", length(hit), "do both; plotting", length(top), "---\n")
print(cand[gene %chin% top][order(match(gene, top), sample),
           .(gene, sample, n, escape = round(escape, 3),
             CI = sprintf("[%.3f-%.3f]", esc_lo, esc_hi), fdr = signif(fdr, 2))])

##### --------------------- THE PLOT --------------------- #####
pl <- x[gene %chin% top & n_umi >= MIN_TILE_UMI]
pl[, AR := a1_umi / n_umi]

nt <- pl[, .(n_tiles = .N), by = .(gene, sample)]
cat("\ntiles with >=", MIN_TILE_UMI, "UMIs (these are the dots):\n")
print(dcast(nt, gene ~ sample, value.var = "n_tiles", fill = 0)[order(match(gene, top))])

# A * marks a gene whose violin was suppressed for want of tiles, so the axis
# itself says which panels are dots-only. Built as "which genes EARN a violin"
# rather than "which are thin", because a gene with no plottable tiles at all is
# missing from nt entirely and would otherwise slip through unmarked.
fat_genes <- nt[, .(ok = all(n_tiles >= MIN_TILES_VIOLIN) && uniqueN(sample) == 2L),
                by = gene][ok == TRUE, gene]
# * = no violin (too few tiles). ! = above biallelic, so not escape whatever
# else it is - only ever present when KEEP_ABOVE_BIALLELIC is TRUE.
labs_v <- setNames(paste0(top,
                          ifelse(top %chin% fat_genes, "", " *"),
                          ifelse(top %chin% above_bi, " !", "")), top)

pl[, gene := factor(gene, levels = top)]
pool <- g[gene %chin% top][, .(gene = factor(gene, levels = top), sample,
                               AR = a1 / n, lo = 1 - esc_hi, hi = 1 - esc_lo)]
fat <- pl[gene %in% fat_genes]   # %in%, not %chin%: gene is a factor by now

p <- ggplot(pl, aes(gene, AR, fill = sample)) +
  geom_hline(yintercept = 1, linetype = 3, colour = "#8a8a8a") +
  geom_hline(yintercept = mean(base$auto_ar), linetype = 2, colour = "#184f95") +
  annotate("text", x = 0.5, y = mean(base$auto_ar), label = "full escape (biallelic)",
           hjust = 0, vjust = -0.6, size = 2.5, colour = "#184f95") +
  { if (nrow(fat)) geom_violin(data = fat, position = position_dodge(width = 0.8),
                               width = 0.8, colour = NA, alpha = 0.4,
                               scale = "width", trim = TRUE) } +
  geom_point(aes(size = n_umi),
             position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.8),
             alpha = 0.35, shape = 16, stroke = 0) +
  # The estimate. Pooled over every molecule, including the tiles too thin to
  # be drawn, with the Wilson interval - which is the part that tells you
  # whether a gene with three dots means anything.
  geom_pointrange(data = pool, aes(gene, AR, ymin = lo, ymax = hi, group = sample),
                  inherit.aes = FALSE, position = position_dodge(width = 0.8),
                  shape = 23, size = 0.35, fill = "white", colour = "#0b0b0b",
                  stroke = 0.6, fatten = 3) +
  scale_fill_manual(values = c("adult (9w)" = "#2D6E5D", "aged (78w)" = "#b02a2a"),
                    name = NULL) +
  scale_size_continuous(range = c(0.4, 2.4), name = "UMIs in tile") +
  scale_x_discrete(labels = labs_v) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "chrX escape candidates: allelic ratio per 64 um tile",
       subtitle = paste0(
         "AR = B6 / total, and B6 is the ACTIVE X, so LOW AR = escape. ",
         "White diamond + bar = pooled ratio over ALL molecules, with 95% CI.\n",
         sprintf("Dots are tiles with >= %d informative UMIs. * = too few tiles for a violin (< %d); read the diamond. ",
                 MIN_TILE_UMI, MIN_TILES_VIOLIN),
         "n = 1 animal per age: a screen, not an age comparison."),
       x = NULL, y = "allelic ratio (B6 / total)",
       caption = paste("Dotted line = fully silenced Xi. Dashed line = the autosomal ratio, which is what full",
                       "escape looks like once B6-ward mapping bias is allowed for.",
                       "\nViolins are lumpy because a tile with few molecules can only take a few AR values -",
                       "that is arithmetic, not biology.")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.subtitle = element_text(size = 8, colour = "#52514e"),
        plot.caption  = element_text(size = 7, colour = "#52514e", hjust = 0))

ggsave(OUT, p, width = 13, height = 6)
cat("\nWrote", normalizePath(OUT), "\n")
