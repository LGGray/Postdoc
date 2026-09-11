# ---------------------------------------------------------------------------
# Shared config, palette and loader for the Allelome.PRO2 per-cell-type
# pseudobulk under OCM/Allelome.PRO2_pseudobulk_celltype.
# Sourced by pseudobulk_celltype_chrX_heatmap.R and
# chrX_distal_escape_enrichment.R; not meant to be run on its own.
#
# Genotype: BL6 Xist-deletion mother x CAST father. XCI is fully skewed by
# construction, CAST is the inactive X in every cell, so no sorting of cells by
# active allele is needed before pseudobulking. A1 = BL6 = active X, and
# AR = A1/total near 1 means monoallelic while lower values mean escape.
# The SNP file is the no_Xist one, so Xist itself is never counted.
#
# Every setting below can be overridden with an environment variable, so the
# same scripts work on a re-run, another cutoff, or the Sham/TAC arm.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(forcats); library(scales)
})

CL   <- Sys.getenv("CLUSTER_MOUNT", "/Users/lachlang/cluster")
ROOT <- Sys.getenv("AP2_ROOT", file.path(CL, "OCM", "Allelome.PRO2_pseudobulk_celltype"))
MIN_READS <- as.integer(Sys.getenv("MIN_READS", "20"))   # SNP-overlapping reads per gene per cell type
MIN_CT    <- as.integer(Sys.getenv("MIN_CT", "6"))       # cell types a gene must be testable in, per group
GROUPS    <- strsplit(Sys.getenv("GROUPS", "9w,78w"), ",")[[1]]   # first = baseline, second = comparison

# Each contrast gets its own folder, so running the Sham/TAC arm cannot
# overwrite the results of the age arm. FIG_OUT still overrides the whole path.
CONTRAST <- Sys.getenv("CONTRAST", paste(GROUPS, collapse = "_vs_"))
OUT <- Sys.getenv("FIG_OUT",
                  file.path("/Users/lachlang/Downloads/OCM_AP2_pseudobulk_celltype", CONTRAST))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
message("contrast ", CONTRAST, "; writing to ", OUT)
ESCAPE_AR <- as.numeric(Sys.getenv("ESCAPE_AR", "0.9"))  # AR <= this is escape (Hoelzl 2025; MONO_AR in the repo)

# chrX geometry. The distal window is the one used in Hoelzl et al. 2025:
# the first 20 Mb and the last 40 Mb of the chromosome.
CHRX_LEN    <- as.numeric(Sys.getenv("CHRX_LEN", "169476592"))    # GRCm39
DISTAL_P_MB <- as.numeric(Sys.getenv("DISTAL_P_MB", "20"))
DISTAL_Q_MB <- as.numeric(Sys.getenv("DISTAL_Q_MB", "40"))
DROP_CT     <- c("CM (stressed)", "Epicardial")          # too few chrX reads to plot

is_distal <- function(start) start <= DISTAL_P_MB * 1e6 | start >= CHRX_LEN - DISTAL_Q_MB * 1e6
DISTAL_LAB <- sprintf("distal window (first %g Mb + last %g Mb)", DISTAL_P_MB, DISTAL_Q_MB)

# ---- style ----------------------------------------------------------------
SAMPLE_LAB <- c("9w" = "Adult (9w)", "78w" = "Aged (78w)", "Sham" = "Sham", "TAC" = "TAC")
SAMPLE_COL <- c("9w" = "#2B7BBA", "Sham" = "#1B9E77", "TAC" = "#7B3294", "78w" = "#E2711D")
AR_BREAKS <- c(seq(0, 0.9, by = 0.1), 0.95, 1.0)
AR_COLS   <- c("#2B3186", "#3B5FB6", "#38749F", "#367373", "#2D6E5D", "#1E652D",
               "#658C2D", "#8D9F25", "#B3B112", "#C97314", "#8B1913")
ar_fill_cont <- function(name = "Allelic ratio\n(B6 / total)") {
  # Colour stops sit at the bin midpoints, but the first and last stop must be
  # pinned to 0 and 1: scales::gradient_n_pal() interpolates with rule = 1, so
  # any ratio outside range(values) returns NA and the tile is drawn white.
  # Without this, AR = 1 (A2_reads == 0, the most monoallelic genes) vanishes.
  mids <- (head(AR_BREAKS, -1) + tail(AR_BREAKS, -1)) / 2
  mids[1] <- 0; mids[length(mids)] <- 1
  scale_fill_gradientn(colours = AR_COLS, values = mids, limits = c(0, 1), name = name, na.value = "white")
}
ESCAPE_GENES <- c("Kdm5c","Kdm6a","Ddx3x","Eif2s3x","Utp14a","Akap17a","Pbdc1","Ftx","Jpx","Sts","5530601H04Rik")
CT_ORDER <- c("Ventricular CM","Fibroblasts","Endothelial","Macrophages","Pericytes / SMC",
              "Endocardium","Lymphatic EC","B cells","T cells","CM (stressed)","Epicardial")
CT_MAP <- c("Ventricular_Cardiomyocytes" = "Ventricular CM",
            "Cardiomyocytes_stressed"       = "CM (stressed)",
            "Pericytes_Smooth_muscle_cells" = "Pericytes / SMC",
            "Endothelial_cells"             = "Endothelial",
            "Lymphatic_endothelial"         = "Lymphatic EC",
            "Epicardial_Mesothelial_cells"  = "Epicardial",
            "B_cells" = "B cells", "T_cells" = "T cells",
            "Fibroblasts" = "Fibroblasts", "Macrophages" = "Macrophages",
            "Endocardium" = "Endocardium")
theme_set(theme_classic(base_size = 15) +
          theme(strip.background = element_blank(),
                strip.text = element_text(face = "bold", size = 13),
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12, colour = "grey30"),
                plot.caption = element_text(size = 10, colour = "grey25", hjust = 0),
                legend.title = element_text(size = 12)))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 200)
  message("wrote ", name)
}

# ---- loader ----------------------------------------------------------------
# Returns one row per group x cell type x chrX gene passing MIN_READS.
load_ap2_chrX <- function(groups = GROUPS, root = ROOT, min_reads = MIN_READS, drop_ct = DROP_CT) {
  read_one <- function(grp, dir) {
    ct_raw <- sub("_annotation_us_mm39_gene_level\\.bed_1$", "",
                  sub("^[0-9]{4}_[0-9]{2}_[0-9]{2}_", "", basename(dir)))
    read_tsv(file.path(dir, "locus_table.txt"), show_col_types = FALSE,
             col_types = cols(chr = col_character(), name = col_character(), .default = col_double())) %>%
      filter(chr == "chrX") %>%
      transmute(sample = grp, celltype_raw = ct_raw, gene = name, start,
                A1 = A1_reads, A2 = A2_reads, total = total_reads)
  }
  dirs <- lapply(groups, function(g) list.files(file.path(root, g), pattern = "bed_1$", full.names = TRUE))
  names(dirs) <- groups
  stopifnot(lengths(dirs) > 0)
  bind_rows(lapply(groups, function(g) bind_rows(lapply(dirs[[g]], read_one, grp = g)))) %>%
    mutate(celltype = unname(CT_MAP[celltype_raw]),
           ar = A1 / total,
           distal = is_distal(start),
           region = factor(ifelse(distal, DISTAL_LAB, "internal"), c(DISTAL_LAB, "internal")),
           sample = factor(sample, groups)) %>%
    filter(!is.na(celltype), !celltype %in% drop_ct, total >= min_reads) %>%
    mutate(celltype = factor(celltype, CT_ORDER)) %>%
    droplevels()
}

# Wide, one row per cell type x gene measured in BOTH groups, with the change.
paired_by_gene <- function(pbg, groups = GROUPS) {
  pbg %>% select(sample, celltype, gene, start, distal, region, ar) %>%
    pivot_wider(names_from = sample, values_from = ar) %>%
    filter(!is.na(.data[[groups[1]]]), !is.na(.data[[groups[2]]])) %>%
    mutate(delta = .data[[groups[2]]] - .data[[groups[1]]])
}

# ---------------------------------------------------------------------------
# Distal-enrichment statistics.
# distal_tests() returns every test table plus a ready-made caption, so any
# figure in this folder can carry the same numbers without recomputing them.
# ---------------------------------------------------------------------------
CLUSTER_MB <- as.numeric(Sys.getenv("CLUSTER_MB", "2.5"))  # single-linkage gap defining one escape cluster
N_PERM     <- as.integer(Sys.getenv("N_PERM", "10000"))

# 2x2 Fisher: rows distal/internal, columns hit/not. One-sided, for enrichment.
fisher_distal <- function(distal, hit, label) {
  tab <- matrix(c(sum(distal & hit), sum(distal & !hit),
                  sum(!distal & hit), sum(!distal & !hit)), nrow = 2, byrow = TRUE,
                dimnames = list(c("distal", "internal"), c("hit", "miss")))
  ft <- fisher.test(tab, alternative = "greater")
  tibble(test = label, n_hits = sum(hit), n_tested = length(hit),
         pct_hits_distal = 100 * sum(distal & hit) / max(1, sum(hit)),
         pct_tested_distal = 100 * mean(distal),
         odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}

# Escapees sit in clusters, so genes are not independent and Fisher overstates
# significance. Rotating the escape-status vector along genes ordered by
# position keeps every cluster intact and only moves where the clusters land.
perm_distal <- function(df, hit_col, label, n_perm = N_PERM) {
  d <- df %>% arrange(start)
  hit <- d[[hit_col]]; distal <- d$distal
  if (sum(hit) == 0)
    return(tibble(test = label, observed_pct_distal = NA_real_, null_mean_pct = NA_real_, p_value = NA_real_))
  obs <- sum(distal & hit) / sum(hit)
  null <- vapply(sample.int(length(hit), n_perm, replace = TRUE), function(k) {
    h <- hit[c(seq.int(k, length(hit)), seq_len(k - 1))]
    sum(distal & h) / sum(h)
  }, numeric(1))
  tibble(test = label, observed_pct_distal = 100 * obs, null_mean_pct = 100 * mean(null),
         p_value = (1 + sum(null >= obs)) / (1 + n_perm))
}

# A new cluster starts wherever the gap to the previous gene exceeds gap_mb.
# Hoelzl et al. describe escape as organised in 2.5-Mb regions, so neighbouring
# escapees are one event rather than several independent ones.
cluster_ids <- function(start, gap_mb = CLUSTER_MB) {
  o <- order(start); ids <- integer(length(start))
  ids[o] <- cumsum(c(1, diff(sort(start)) > gap_mb * 1e6))
  ids
}

distal_tests <- function(pbg, groups = GROUPS, seed = as.integer(Sys.getenv("SEED", "1"))) {
  set.seed(seed)
  BASE <- groups[1]; COMP <- groups[2]

  # gene-level pseudobulk: pool reads over cell types, the level the published
  # escape calls are made at. Cell-type resolution is kept for the Wilcoxon.
  gene_lvl <- pbg %>%
    group_by(sample, gene, start, distal, region) %>%
    summarise(A1 = sum(A1), A2 = sum(A2), n_ct = n(), .groups = "drop") %>%
    mutate(total = A1 + A2, ar = A1 / total, escape = ar <= ESCAPE_AR)

  wide <- gene_lvl %>% select(sample, gene, start, distal, region, ar, escape, n_ct) %>%
    pivot_wider(names_from = sample, values_from = c(ar, escape, n_ct)) %>%
    filter(!is.na(.data[[paste0("ar_", BASE)]]), !is.na(.data[[paste0("ar_", COMP)]])) %>%
    mutate(delta  = .data[[paste0("ar_", COMP)]] - .data[[paste0("ar_", BASE)]],
           gained = !.data[[paste0("escape_", BASE)]] & .data[[paste0("escape_", COMP)]])

  gained_lab <- sprintf("escape gained %s -> %s", BASE, COMP)

  # 1. Fisher exact
  fish <- bind_rows(
    fisher_distal(gene_lvl$distal[gene_lvl$sample == BASE], gene_lvl$escape[gene_lvl$sample == BASE],
                  sprintf("escapees in %s", BASE)),
    fisher_distal(gene_lvl$distal[gene_lvl$sample == COMP], gene_lvl$escape[gene_lvl$sample == COMP],
                  sprintf("escapees in %s", COMP)),
    fisher_distal(wide$distal, wide$gained, gained_lab))

  # 2. circular permutation, preserving escape clusters
  perm <- bind_rows(
    perm_distal(gene_lvl %>% filter(sample == BASE) %>% mutate(h = escape), "h", sprintf("escapees in %s", BASE)),
    perm_distal(gene_lvl %>% filter(sample == COMP) %>% mutate(h = escape), "h", sprintf("escapees in %s", COMP)),
    perm_distal(wide %>% mutate(h = gained), "h", gained_lab))

  # 3. Wilcoxon on the change in allelic ratio, no escape cutoff
  pairs <- paired_by_gene(pbg, groups)
  wil <- bind_rows(
    tibble(level = "gene (reads pooled over cell types)",
           n_distal = sum(wide$distal), n_internal = sum(!wide$distal),
           median_delta_distal = median(wide$delta[wide$distal]),
           median_delta_internal = median(wide$delta[!wide$distal]),
           p_value = wilcox.test(delta ~ distal, data = wide)$p.value),
    tibble(level = "gene x cell type",
           n_distal = sum(pairs$distal), n_internal = sum(!pairs$distal),
           median_delta_distal = median(pairs$delta[pairs$distal]),
           median_delta_internal = median(pairs$delta[!pairs$distal]),
           p_value = wilcox.test(delta ~ distal, data = pairs)$p.value))

  # 4. window sweep, keeping the published 1:2 p-arm to q-arm shape
  sweep <- bind_rows(lapply(c(5, 10, 15, 20, 25, 30, 40, 50), function(mb) {
    dis <- wide$start <= mb * 1e6 | wide$start >= CHRX_LEN - 2 * mb * 1e6
    fisher_distal(dis, wide$gained, sprintf("p %g Mb / q %g Mb", mb, 2 * mb)) %>% mutate(window_mb = mb)
  }))

  # 5. cluster level: the honest unit of replication, and the number that says
  # how much power tests 2 and 5 can possibly have.
  clus <- wide %>% filter(gained) %>% mutate(cluster = cluster_ids(start)) %>%
    group_by(cluster) %>%
    summarise(n_genes = n(), start = median(start), genes = paste(gene, collapse = ", "), .groups = "drop") %>%
    mutate(distal = is_distal(start))
  p_bg <- mean(wide$distal)
  cluster_test <- tibble(n_clusters = nrow(clus), n_distal = sum(clus$distal),
                         expected_distal = p_bg * nrow(clus),
                         best_possible_p = p_bg ^ nrow(clus),
                         p_value = binom.test(sum(clus$distal), nrow(clus), p = p_bg,
                                              alternative = "greater")$p.value)

  hit <- fish %>% filter(test == gained_lab)
  caption <- sprintf(
    "Distal enrichment of newly escaping genes: %.0f%% of the %d genes gaining escape lie in the %s, against %.0f%% of all %d tested genes.\nFisher exact one-sided OR = %.1f, p = %s. Wilcoxon on change in allelic ratio, distal vs internal, p = %s.\nThose genes form only %d escape clusters at %g Mb, %d of them distal, binomial p = %s, and with this many clusters no result could fall below p = %.3f.",
    hit$pct_hits_distal, hit$n_hits, DISTAL_LAB, hit$pct_tested_distal, hit$n_tested,
    hit$odds_ratio, format.pval(hit$p_value, digits = 2),
    format.pval(wil$p_value[wil$level == "gene x cell type"], digits = 2),
    cluster_test$n_clusters, CLUSTER_MB, cluster_test$n_distal,
    format.pval(cluster_test$p_value, digits = 2), cluster_test$best_possible_p)

  list(fisher = fish, permutation = perm, cluster = cluster_test, clusters_gained = clus,
       wilcoxon = wil, window_sweep = sweep, gene_level = gene_lvl, wide = wide, caption = caption)
}
