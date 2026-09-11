# ---------------------------------------------------------------------------
# Are escape genes enriched at the distal ends of chrX, and does that get
# stronger with age?
#
# Tests the positional claim of Hoelzl et al., Nat Aging 2025
# (doi:10.1038/s43587-025-00856-8), which reported 84% of age-specific escapees
# in the first 20 Mb or last 40 Mb of chrX, against the Allelome.PRO2
# per-cell-type pseudobulk in OCM/Allelome.PRO2_pseudobulk_celltype.
#
# Five tests, run by distal_tests() in ap2_pseudobulk_common.R:
#   1. Fisher exact  - escapee (AR <= ESCAPE_AR) vs not, distal vs internal, in
#                      each group and for genes that only escape in the older one
#   2. Permutation   - same contrast, rotating the escape-status vector along
#                      the ordered gene list, which preserves the clustering of
#                      escapees that makes genes non-independent
#   3. Wilcoxon      - change in allelic ratio, distal vs internal, no cutoff
#   4. Window sweep  - test 1 across window widths, so the answer cannot be an
#                      artefact of the published 20/40 Mb choice
#   5. Cluster level - collapse newly escaping genes into escape clusters and
#                      test those, since neighbouring escapees are one event
#
#   Rscript OCM_heart/chrX_distal_escape_enrichment.R
#
# Config, all optional environment variables, lives in ap2_pseudobulk_common.R:
# CLUSTER_MOUNT, AP2_ROOT, FIG_OUT, MIN_READS, MIN_CT, GROUPS, ESCAPE_AR,
# CHRX_LEN, DISTAL_P_MB, DISTAL_Q_MB, CLUSTER_MB, N_PERM, SEED. For the stress
# arm instead of the age arm:
#   GROUPS=Sham,TAC Rscript OCM_heart/chrX_distal_escape_enrichment.R
# ---------------------------------------------------------------------------
.this_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(sub("^--file=", "", f[1])) else "OCM_heart"
}
source(file.path(.this_dir(), "ap2_pseudobulk_common.R"))
suppressPackageStartupMessages({ library(patchwork); library(ggrepel) })

BASE <- GROUPS[1]; COMP <- GROUPS[2]
LAB <- function(g) sub(" \\(.*", "", SAMPLE_LAB[[g]])
CONTRAST_LAB <- sprintf("%s vs %s", tolower(LAB(BASE)), tolower(LAB(COMP)))

pbg <- load_ap2_chrX()
message(nrow(pbg), " gene x cell type x group observations at >= ", MIN_READS, " reads")

st <- distal_tests(pbg)
gene_lvl <- st$gene_level; wide <- st$wide; CAPTION <- st$caption
message(nrow(wide), " chrX genes measured in both groups")
for (nm in c("fisher", "permutation", "cluster", "clusters_gained", "wilcoxon", "window_sweep")) {
  write_tsv(st[[nm]], file.path(OUT, paste0("distal_enrichment_", nm, ".tsv")))
  cat("\n=== ", nm, " ===\n", sep = ""); print(as.data.frame(st[[nm]]), digits = 3)
}
writeLines(CAPTION, file.path(OUT, "distal_enrichment_caption.txt"))
cat("\n", CAPTION, "\n", sep = "")

# ---- figure ----------------------------------------------------------------
shade <- list(
  annotate("rect", xmin = -Inf, xmax = DISTAL_P_MB, ymin = -Inf, ymax = Inf, fill = "#E2711D", alpha = 0.10),
  annotate("rect", xmin = (CHRX_LEN / 1e6) - DISTAL_Q_MB, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#E2711D", alpha = 0.10))
lab_genes <- wide %>% filter(gained | delta < -0.1) %>% arrange(delta) %>% head(12)

pA <- ggplot(gene_lvl, aes(start / 1e6, ar)) + shade +
  geom_hline(yintercept = ESCAPE_AR, linetype = 2, colour = "grey40") +
  geom_point(aes(colour = ar <= ESCAPE_AR), size = 1.6, alpha = 0.85) +
  facet_wrap(~sample, ncol = 1, labeller = labeller(sample = SAMPLE_LAB)) +
  scale_colour_manual(values = c(`TRUE` = "#2B3186", `FALSE` = "grey55"),
                      labels = c(`TRUE` = sprintf("escape (AR <= %.2f)", ESCAPE_AR), `FALSE` = "monoallelic"),
                      name = NULL) +
  scale_x_continuous(limits = c(0, CHRX_LEN / 1e6), expand = c(0.01, 0)) +
  labs(x = NULL, y = "Allelic ratio (B6 / total)",
       title = sprintf("Escape from X inactivation along chrX, %s", CONTRAST_LAB),
       subtitle = sprintf("Gene-level pseudobulk, reads pooled over cell types, >= %d SNP-overlapping reads. Shaded: %s",
                          MIN_READS, DISTAL_LAB)) +
  theme(legend.position = "top")

pB <- ggplot(wide, aes(start / 1e6, delta)) + shade +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(aes(colour = gained), size = 1.8, alpha = 0.9) +
  geom_text_repel(data = lab_genes, aes(label = gene), size = 3, min.segment.length = 0,
                  max.overlaps = 20, seed = 1) +
  scale_colour_manual(values = c(`TRUE` = "#8B1913", `FALSE` = "grey60"),
                      labels = c(`TRUE` = "newly escaping", `FALSE` = "no change in call"), name = NULL) +
  scale_x_continuous(limits = c(0, CHRX_LEN / 1e6), expand = c(0.01, 0)) +
  labs(x = "chrX position (Mb, GRCm39)",
       y = sprintf("Change in allelic ratio\n(%s - %s)", COMP, BASE), caption = CAPTION) +
  theme(legend.position = "top")

save_fig(pA / pB + plot_layout(heights = c(2, 1.35)), "AP2_chrX_distal_escape_enrichment", 11, 9)

# escape fraction per 10-Mb bin, the simplest view of the same result
binned <- gene_lvl %>% mutate(bin = floor(start / 1e7) * 10) %>%
  group_by(sample, bin) %>% summarise(n = n(), pct = 100 * mean(escape), .groups = "drop") %>%
  filter(n >= 5)
pC <- ggplot(binned, aes(bin + 5, pct, fill = sample)) + shade +
  geom_col(position = position_dodge(width = 8), width = 7.5) +
  scale_fill_manual(values = c("#2B7BBA", "#E2711D")[seq_along(GROUPS)],
                    labels = SAMPLE_LAB[GROUPS], name = NULL) +
  scale_x_continuous(limits = c(0, CHRX_LEN / 1e6), expand = c(0.01, 0)) +
  labs(x = "chrX position (Mb, GRCm39)", y = sprintf("%% of tested genes with AR <= %.2f", ESCAPE_AR),
       title = "Escape per 10-Mb window along chrX",
       subtitle = sprintf("Bins with at least 5 tested genes. Shaded: %s", DISTAL_LAB), caption = CAPTION) +
  theme(legend.position = "top")
save_fig(pC, "AP2_chrX_escape_by_position_bins", 11, 5)
