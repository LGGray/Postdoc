# The gene panel, shared by tile_gene_ar_panel.R (distributions) and
# tile_gene_ar_maps.R (slide maps), so the two figures cannot drift apart.
#
# Defines PANEL, SUBS and the PANEL_GENES helper. Sourced, not run.
#
# GENOTYPE, which sets the direction of every number downstream (same as
# escape_genes.R and tile_ratio_map.R): B6 mother x CAST father, and B6 carries
# the Xist deletion, so the CAST X is the inactive X in EVERY cell.
#   a1 = B6 = REF = maternal = the active X      AR = a1 / (a1 + a2)
#   a2 = CAST = ALT = paternal = the inactive X
# On chrX:  AR ~ 1 is fully silenced, and escape = 1 - AR.
# On the imprinted loci: AR ~ 1 = maternally expressed, AR ~ 0 = paternally
# expressed - which is what `expect` below records.

# `expect` is the allelic ratio the locus SHOULD sit at on prior grounds:
# 1 = maternal = B6, 0 = paternal = CAST, NA = no prior (chrX escape is the
# thing being measured). Imprinting direction is the canonical one; Grb10 is
# maternal in peripheral tissue and paternal in CNS, which is heart = maternal.
PANEL <- rbindlist(list(
  data.table(gene = "chrX (all genes pooled)", group = "Skew", expect = NA_real_),

  # BIALLELIC REFERENCE. An autosomal gene under no allelic regulation should sit
  # at the autosomal ratio, which is 0.515 per molecule at 9w - not 0.5, and the
  # gap is the mapping bias toward the B6 reference. 0.5 is recorded as `expect`
  # because that is the biological expectation; the measured autosomal ratio is
  # printed beside it by tile_gene_ar_panel.R, and the distance between the two
  # is the bias. These are the only genes on the panel whose ratio is known in
  # advance without an imprinting or XCI argument, which is what makes them the
  # reference the imprinted controls are read against.
  #
  # MEASURED ON THIS DATA, so the pages are not a surprise:
  #   Malat1  1734 molecules AR 0.498 (9w), 1267 AR 0.499 (78w), 29 SNPs over
  #           7.0 kb - textbook, and the best biallelic control in this dataset.
  #   Gapdh     84 molecules AR 0.940 (9w),   13 AR 1.000 (78w), 11 SNPs over
  #           4.6 kb. NOT biallelic as measured, and on an autosome that is a
  #           method signal rather than biology: 79 B6 against 5 CAST cannot be
  #           allelic regulation of Gapdh. Thin SNP support over a short, highly
  #           expressed gene with many processed pseudogenes is the likely
  #           cause. It is on the panel BECAUSE it fails - a control that only
  #           ever contains passing genes cannot show the reader what failure
  #           looks like - but do not quote it as a biallelic reference.
  # Neat1 (chr19, 140 SNPs, 2604 molecules AR 0.447 at 9w / 5533 AR 0.504 at 78w)
  # is the obvious third if a deeper one is wanted; it is left off so the group
  # stays the two that were asked for.
  data.table(group = "Biallelic (control)", expect = 0.5, gene = c(
    "Malat1", "Gapdh")),

  data.table(group = "Imprinted (control)", expect = 0, gene = c(
    "Mcts2", "Sgce", "Peg3", "Snhg14", "Snurf", "Snrpn", "Igf2",
    "Kcnq1ot1", "Plagl1", "Zrsr1", "Peg13", "Slc38a4", "Airn")),
  data.table(group = "Imprinted (control)", expect = 1, gene = c(
    "Zim1", "H19", "Cdkn1c", "Grb10", "Meg3", "Rian", "Igf2r")),

  data.table(group = "chrX escape (heart, age)", expect = NA_real_, gene = c(
    "Shroom4", "Tspan7", "Sh3kbp1", "Med14", "Kctd12b", "2210013O21Rik",
    "Plp1", "4930578C19Rik", "Smpx", "Slc16a2", "Ftx", "Utp14a", "Pbdc1",
    "Ddx3x", "Kdm5c", "Kdm6a"))
), use.names = TRUE)
# Biallelic first after Skew: it is the reference the other two groups are read
# against, and tile_gene_ar_panel.R pages the groups in this order.
PANEL[, group := factor(group, levels = c("Skew", "Biallelic (control)",
                                          "Imprinted (control)",
                                          "chrX escape (heart, age)"))]

# Stand-ins for panel members the annotation does not carry. Plotted next to the
# gene they replace and labelled, never silently swapped in for it.
#
# Absent under these names from annotation_us.bed, so no read was ever assigned
# to them: Mcts2, Snurf, Snrpn, Zrsr1, Peg13, 2210013O21Rik, 4930578C19Rik.
# Snrpn is the real loss - it is the cleanest paternal control there is, and the
# cross direction was originally fixed on it.
#
# KCNQ1OT1 IS NOT ONE OF THEM, and the earlier version of this comment saying so
# was wrong. It IS in annotation_us.bed, at chr7:142,766,847-142,850,284, and it
# has thousands of informative SNPs. It reads zero molecules for a different
# reason: it sits ENTIRELY INSIDE Kcnq1 (chr7:142,660,613-142,980,787), and
# Intervals.name_at() in ase_bin_allele_counts.py attributes a position to the
# containing interval with the SMALLEST START. Kcnq1 starts 106 kb earlier, so
# it takes every SNP in the Kcnq1ot1 body and Kcnq1ot1 gets none.
#
# That also means KCNQ1 IS NOT A STAND-IN FOR IT - it is the two loci summed.
# Kcnq1 is maternally expressed and Kcnq1ot1 is its paternally expressed
# antisense, over the same coordinates, so the pooled ratio is a mixture of
# opposite directions: AR 0.602 at 9w and 0.494 at 78w, where a clean maternal
# control should read ~1. Plotting it as an imprinted control was showing a
# reciprocal pair averaged to the middle. Kcnq1 is therefore REMOVED from the
# substitutes rather than relabelled.
#
# To get Kcnq1ot1 itself, the SNPs have to be re-attributed, which is a counting
# change and not a plotting one: slurm/spatial_tile_nested_loci.slurm recounts
# chr7 alone against a bed in which Kcnq1ot1 is the only interval over its body,
# and writes a small side tree that EXTRA_ROOT merges into these figures.
#
# Xist is absent for a different and unfixable reason: the mask these trees were
# counted against is SNPfile_..._mm39_no_Xist.bed, and even with SNPs, B6 carries
# the Xist DELETION, so every Xist read is CAST by construction and its ratio
# would read ~0 whatever the skew is. Xist cannot report skew in this genotype.
# "chrX (all genes pooled)" above is the quantity it would have been a proxy for.
SUBS <- data.table(
  gene = "Zrsr2",
  for_gene = "Zrsr1",
  group = "Imprinted (control)",
  # Zrsr2 is the chrX paralogue of Zrsr1 and is NOT imprinted - it is subject to
  # XCI, so its expectation is the chrX one, not an imprinted one.
  expect = NA_real_)

# Panel order: the order it was asked for, with each substitute directly after
# the gene it stands in for, so the reader finds it where the missing one was.
# `have` restricts to what is actually in the counts.
PANEL_GENES <- function(have = NULL, subs = SUBS) {
  ord <- PANEL$gene
  for (i in seq_len(nrow(subs))) {
    at <- match(subs$for_gene[i], ord)
    ord <- append(ord, subs$gene[i], after = if (is.na(at)) length(ord) else at)
  }
  if (is.null(have)) ord else ord[ord %chin% have]
}

# gene -> group/expect, for both scripts, on character throughout so the group
# levels cannot come out reordered by whichever side happened to be first.
PANEL_META <- function(subs = SUBS) {
  m <- rbind(PANEL[, .(gene, group = as.character(group), expect)],
             subs[,  .(gene, group = as.character(group), expect)])
  m[, group := factor(group, levels = levels(PANEL$group))][]
}
