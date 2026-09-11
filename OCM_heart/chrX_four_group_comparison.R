# ---------------------------------------------------------------------------
# Adult, Sham, TAC and aged side by side on a matched gene set.
#
# The pairwise scripts test each arm on whatever genes that arm can measure,
# so the age arm (294 testable genes) and the stress arm (400) are not
# comparable as they stand. Here every group is restricted to the genes
# testable in ALL FOUR, which makes the distal-enrichment numbers directly
# comparable and separates three different questions:
#
#   Sham vs adult  - does surgery and batch alone move allelic ratios?
#   TAC vs Sham    - the effect of pressure overload, against its own control
#   Aged vs adult  - the effect of age
#
#   Rscript OCM_heart/chrX_four_group_comparison.R
#
# Config as in ap2_pseudobulk_common.R. GROUPS defaults to the four arms in
# the order they are plotted; the first is the reference for Sham and aged.
# ---------------------------------------------------------------------------
if (Sys.getenv("GROUPS") == "")   Sys.setenv(GROUPS = "9w,Sham,TAC,78w")
if (Sys.getenv("CONTRAST") == "") Sys.setenv(CONTRAST = "four_group")
.this_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(sub("^--file=", "", f[1])) else "OCM_heart"
}
source(file.path(.this_dir(), "ap2_pseudobulk_common.R"))
suppressPackageStartupMessages(library(patchwork))
set.seed(as.integer(Sys.getenv("SEED", "1")))

stopifnot(length(GROUPS) >= 2)
REF <- GROUPS[1]

# ---- matched gene set ------------------------------------------------------
pbg_all <- load_ap2_chrX(groups = GROUPS)
matched <- pbg_all %>% distinct(sample, gene) %>% count(gene) %>%
  filter(n == length(GROUPS)) %>% pull(gene)
pbg <- pbg_all %>% filter(gene %in% matched)
message(length(matched), " genes testable at >= ", MIN_READS, " reads in all ",
        length(GROUPS), " groups, out of ",
        n_distinct(pbg_all$gene), " testable in at least one")
write_tsv(pbg %>% select(sample, celltype, gene, start, distal, A1, A2, total, ar),
          file.path(OUT, "chrX_matched_geneset_allelic_ratio.tsv"))

gene_lvl <- pbg %>%
  group_by(sample, gene, start, distal, region) %>%
  summarise(A1 = sum(A1), A2 = sum(A2), n_ct = n(), .groups = "drop") %>%
  mutate(total = A1 + A2, ar = A1 / total, escape = ar <= ESCAPE_AR)

# ---- per group: how much escape, and how distal is it ----------------------
per_group <- bind_rows(lapply(GROUPS, function(g) {
  d <- gene_lvl %>% filter(sample == g)
  fisher_distal(d$distal, d$escape, g) %>%
    mutate(pct_escape = 100 * mean(d$escape),
           pct_escape_distal = 100 * mean(d$escape[d$distal]),
           pct_escape_internal = 100 * mean(d$escape[!d$distal]))
})) %>% rename(group = test)

# ---- three contrasts, each on the same matched gene set --------------------
CONTRASTS <- list(c(REF, "Sham"), c("Sham", "TAC"), c(REF, "78w"))
CONTRASTS <- Filter(function(x) all(x %in% GROUPS) && x[1] != x[2], CONTRASTS)

contrast_one <- function(base, comp) {
  w <- gene_lvl %>% filter(sample %in% c(base, comp)) %>%
    select(sample, gene, start, distal, ar, escape) %>%
    pivot_wider(names_from = sample, values_from = c(ar, escape)) %>%
    mutate(delta  = .data[[paste0("ar_", comp)]] - .data[[paste0("ar_", base)]],
           gained = !.data[[paste0("escape_", base)]] & .data[[paste0("escape_", comp)]])
  lab <- sprintf("%s -> %s", base, comp)
  clus <- w %>% filter(gained) %>% mutate(cluster = cluster_ids(start)) %>%
    group_by(cluster) %>%
    summarise(n_genes = n(), start = median(start), genes = paste(gene, collapse = ", "), .groups = "drop") %>%
    mutate(distal = is_distal(start))
  p_bg <- mean(w$distal)
  fisher_distal(w$distal, w$gained, lab) %>%
    mutate(median_delta_distal = median(w$delta[w$distal]),
           median_delta_internal = median(w$delta[!w$distal]),
           wilcox_p = wilcox.test(delta ~ distal, data = w)$p.value,
           n_clusters = nrow(clus), n_clusters_distal = sum(clus$distal),
           cluster_p = if (nrow(clus) > 0)
             binom.test(sum(clus$distal), nrow(clus), p = p_bg, alternative = "greater")$p.value else NA_real_,
           cluster_best_possible_p = p_bg ^ nrow(clus))
}
contrasts_tbl <- bind_rows(lapply(CONTRASTS, function(x) contrast_one(x[1], x[2])))
gained_genes <- bind_rows(lapply(CONTRASTS, function(x) {
  w <- gene_lvl %>% filter(sample %in% x) %>% select(sample, gene, start, distal, ar, escape) %>%
    pivot_wider(names_from = sample, values_from = c(ar, escape))
  w %>% filter(!.data[[paste0("escape_", x[1])]], .data[[paste0("escape_", x[2])]]) %>%
    transmute(contrast = sprintf("%s -> %s", x[1], x[2]), gene, Mb = start / 1e6, distal,
              ar_base = .data[[paste0("ar_", x[1])]], ar_comp = .data[[paste0("ar_", x[2])]])
}))

for (nm in c("per_group", "contrasts_tbl", "gained_genes")) {
  write_tsv(get(nm), file.path(OUT, paste0("four_group_", nm, ".tsv")))
  cat("\n=== ", nm, " ===\n", sep = ""); print(as.data.frame(get(nm)), digits = 3)
}

CAPTION <- sprintf(
  "All four arms restricted to the %d chrX genes testable at >= %d reads in every group, so the contrasts are directly comparable.\nDistal window: the first %g Mb plus the last %g Mb of chrX. Escape called at AR <= %.2f. One animal per group.",
  length(matched), MIN_READS, DISTAL_P_MB, DISTAL_Q_MB, ESCAPE_AR)
writeLines(CAPTION, file.path(OUT, "four_group_caption.txt"))

# ---- figures ---------------------------------------------------------------
COLS <- SAMPLE_COL[GROUPS]
shade <- list(
  annotate("rect", xmin = -Inf, xmax = DISTAL_P_MB, ymin = -Inf, ymax = Inf, fill = "#E2711D", alpha = 0.10),
  annotate("rect", xmin = (CHRX_LEN / 1e6) - DISTAL_Q_MB, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#E2711D", alpha = 0.10))

# A. escape rate per group, split by chromosome region. The headline panel:
#    escape rises in both TAC and aged, but only aged pulls the distal bar up.
bars <- per_group %>%
  select(group, distal = pct_escape_distal, internal = pct_escape_internal) %>%
  pivot_longer(-group, names_to = "where", values_to = "pct") %>%
  mutate(group = factor(group, GROUPS),
         where = factor(where, c("distal", "internal"), c(DISTAL_LAB, "internal chrX")))
pA <- ggplot(bars, aes(group, pct, fill = group, alpha = where)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.4, colour = "grey25", show.legend = FALSE) +
  scale_fill_manual(values = COLS, guide = "none") +
  scale_alpha_manual(values = c(1, 0.35), name = NULL) +
  scale_x_discrete(labels = SAMPLE_LAB[GROUPS]) +
  expand_limits(y = max(bars$pct) * 1.15) +
  labs(x = NULL, y = sprintf("Escape rate\n(%% of genes, AR <= %.2f)", ESCAPE_AR),
       title = "Escape rate by group, split by position on chrX",
       subtitle = "Solid: distal window. Pale: internal chromosome") +
  theme(legend.position = "top")

# B. where escape sits along the chromosome, one row per group
pB <- ggplot(gene_lvl %>% mutate(sample = factor(sample, GROUPS)), aes(start / 1e6, ar)) + shade +
  geom_hline(yintercept = ESCAPE_AR, linetype = 2, colour = "grey40") +
  geom_point(aes(colour = escape), size = 1.4, alpha = 0.85) +
  facet_wrap(~sample, ncol = 1, labeller = labeller(sample = SAMPLE_LAB)) +
  scale_colour_manual(values = c(`TRUE` = "#2B3186", `FALSE` = "grey60"),
                      labels = c(`TRUE` = sprintf("escape (AR <= %.2f)", ESCAPE_AR), `FALSE` = "monoallelic"),
                      name = NULL) +
  scale_x_continuous(limits = c(0, CHRX_LEN / 1e6), expand = c(0.01, 0)) +
  labs(x = "chrX position (Mb, GRCm39)", y = "Allelic ratio (B6 / total)",
       title = "Escape along chrX in each group", caption = CAPTION) +
  theme(legend.position = "top")

save_fig(pA / pB + plot_layout(heights = c(1, 2.4)), "AP2_chrX_four_group_summary", 11, 11)

# C. the matched heatmap: genes testable in >= MIN_CT cell types in every group
keep <- pbg %>% count(sample, gene) %>% filter(n >= MIN_CT) %>% count(gene) %>%
  filter(n == length(GROUPS)) %>% pull(gene)
d <- pbg %>% filter(gene %in% keep) %>%
  mutate(gene = fct_reorder(gene, start), sample = factor(sample, GROUPS))
message(n_distinct(d$gene), " genes testable in >= ", MIN_CT, " cell types in every group")
DISTAL_GENES <- unique(pbg$gene[pbg$distal])
pC <- ggplot(d, aes(gene, celltype, fill = ar)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_wrap(~sample, ncol = 1, labeller = labeller(sample = SAMPLE_LAB)) +
  ar_fill_cont() + scale_y_discrete(limits = rev, drop = TRUE) +
  labs(x = sprintf("chrX genes testable in >= %d cell types in every group (left = distal Xp, ordered by position)", MIN_CT),
       y = NULL, title = "chrX allelic ratio per gene and cell type: adult, Sham, TAC, aged",
       subtitle = sprintf("Matched gene set. Orange gene labels: %s", DISTAL_LAB), caption = CAPTION) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7,
                                   face = ifelse(levels(d$gene) %in% ESCAPE_GENES, "bold", "plain"),
                                   colour = ifelse(levels(d$gene) %in% DISTAL_GENES, "#B5540F", "grey25")),
        panel.grid = element_blank(), legend.position = "right")
save_fig(pC, "AP2_chrX_four_group_heatmap", 13, 12)
