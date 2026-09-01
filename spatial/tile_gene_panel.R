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
PANEL[, group := factor(group, levels = c("Skew", "Imprinted (control)",
                                          "chrX escape (heart, age)"))]

# Stand-ins for panel members the annotation does not carry. Plotted next to the
# gene they replace and labelled, never silently swapped in for it.
#
# Absent under these names from annotation_us.bed, so no read was ever assigned
# to them: Mcts2, Snurf, Snrpn, Kcnq1ot1, Zrsr1, Peg13, 2210013O21Rik,
# 4930578C19Rik. Snrpn is the real loss - it is the cleanest paternal control
# there is, and the cross direction was originally fixed on it.
#
# Xist is absent for a different and unfixable reason: the mask these trees were
# counted against is SNPfile_..._mm39_no_Xist.bed, and even with SNPs, B6 carries
# the Xist DELETION, so every Xist read is CAST by construction and its ratio
# would read ~0 whatever the skew is. Xist cannot report skew in this genotype.
# "chrX (all genes pooled)" above is the quantity it would have been a proxy for.
SUBS <- data.table(
  gene = c("Kcnq1", "Zrsr2"),
  for_gene = c("Kcnq1ot1", "Zrsr1"),
  group = "Imprinted (control)",
  # Kcnq1 is the maternally expressed sense gene of the Kcnq1ot1 locus. Zrsr2 is
  # the chrX paralogue of Zrsr1 and is NOT imprinted - it is subject to XCI, so
  # its expectation is the chrX one, not an imprinted one.
  expect = c(1, NA_real_))

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
