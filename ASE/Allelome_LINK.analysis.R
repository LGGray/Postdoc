library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(UpSetR)
library(patchwork)
library(eulerr)
library(magrittr)
library(ggalluvial)

age.colours <- c("adult" = "#94A2AB", "aged" = "#284456")

#### Match RefSeq transcript ID to gene name
gtf <- read.delim('../GRCm39/GCF_000001635.27_GRCm39_genomic.gtf', comment.char = "#", header = FALSE)
gtf_transcripts <- subset(gtf, V3 == "transcript")

attrs <- strsplit(as.character(gtf_transcripts$V9), "; ", fixed = TRUE)
gene_ids <- vapply(attrs, function(x) {
  if (length(x) >= 1) sub("^gene_id\\s*", "", x[1]) else NA_character_
}, character(1))
transcript_ids <- vapply(attrs, function(x) {
  ti <- grep("^transcript_id", x, value = TRUE)
  if (length(ti)) sub("^transcript_id\\s*", "", ti[1]) else NA_character_
}, character(1))
biotype_ids <- vapply(attrs, function(x) {
  bt <- grep("^transcript_biotype", x, value = TRUE)
  if (length(bt)) sub("^transcript_biotype\\s*", "", bt[1]) else NA_character_
}, character(1))

# make gtf data.tables
dt_gtf <- data.table(gtf_transcripts)[, .(transcript_id = transcript_ids, gene_ids, biotype_ids)]

#### Read in Allelome.PRO2 ouput
adult.files <- list.files('9w', pattern = 'locus_table.txt', full.names = TRUE, recursive = TRUE)
adult.Allelome <- lapply(adult.files, function(file) {
  read.delim(file)
})
names(adult.Allelome) <- unlist(lapply(adult.files, function(x){
  tmp <- strsplit(x, '/')[[1]][[3]]
  strsplit(tmp, '_')[[1]][4]
}))

# Add gene name and biotype
adult.Allelome <- lapply(adult.Allelome, function(x) {
  dt_x <- data.table(x)
  dt_x[, gene := dt_gtf[.SD, on = .(transcript_id = name), gene_ids]]
  dt_x[, biotype := dt_gtf[.SD, on = .(transcript_id = name), biotype_ids]]
  as.data.frame(dt_x)
})

save(adult.Allelome, file = "adult.Allelome.PRO2.RData")
load("adult.Allelome.PRO2.RData")

aged.files <- list.files('78w', pattern = 'locus_table.txt', full.names = TRUE, recursive = TRUE)
aged.Allelome <- lapply(aged.files, function(file) {
  read.delim(file)
})
names(aged.Allelome) <- unlist(lapply(aged.files, function(x){
  tmp <- strsplit(x, '/')[[1]][[3]]
  strsplit(tmp, '_')[[1]][4]
}))

# Add gene name and biotype
aged.Allelome <- lapply(aged.Allelome, function(x) {
  dt_x <- data.table(x)
  dt_x[, gene := dt_gtf[.SD, on = .(transcript_id = name), gene_ids]]
  dt_x[, biotype := dt_gtf[.SD, on = .(transcript_id = name), biotype_ids]]
  as.data.frame(dt_x)
})

save(aged.Allelome, file = "aged.Allelome.PRO2.RData")


# Violin plot of allelic ratios across cell types split by age
combined.adult <- bind_rows(adult.Allelome, .id = "celltype")
combined.adult$age <- 'adult'
combined.aged <- bind_rows(aged.Allelome, .id = "celltype")
combined.aged$age <- 'aged'
all.combined <- bind_rows(combined.adult, combined.aged)
all.combined_autosome <- subset(all.combined, chr != 'X')

pdf("figures/autosomal_allelic_ratios_adult_aged.pdf", width = 8, height = 6)
ggplot(all.combined_autosome, aes(x = celltype, y = allelic_ratio, color = age)) +
  # draw a semi-transparent violin split by age
  geom_violin(trim = TRUE, aes(fill = age), alpha = 0.3,
              position = position_dodge(width = 0.8)) +
  # overlay boxplots at the same dodged positions
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  # dashed reference lines
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "black") +
  labs(x = "", y = "Allelic Ratio") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# Proportion of ASE genes across celltypes
proportion_ase_df <- all.combined_autosome %>%
  group_by(celltype, age) %>%
  summarise(
    total_genes   = n(),
    ase_genes     = sum(allelic_ratio <= 0.3 | allelic_ratio >= 0.7),
    proportion_ase = ase_genes / total_genes,

    # proportion of ALL genes that are ASE AND protein-coding / non-coding
    ase_pc_gene   = sum((allelic_ratio <= 0.3 | allelic_ratio >= 0.7) & biotype == "mRNA") / total_genes,
    ase_nc_gene   = sum((allelic_ratio <= 0.3 | allelic_ratio >= 0.7) & biotype == "lnc_RNA") / total_genes,
    .groups = "drop"
  )

# Stacked bar plot of proportion_ase split by age but filled by ase_pc_gene and ase_nc_gene
stacked_plot_df <- proportion_ase_df %>%
  transmute(
    celltype, age, ase_genes,
    `Pc-genes` = ase_pc_gene * 100,
    `Nc-genes` = ase_nc_gene * 100
  ) %>%
  pivot_longer(
    cols = c(`Pc-genes`, `Nc-genes`),
    names_to = "biotype",
    values_to = "percent"
  ) %>%
  mutate(biotype = factor(biotype, levels = c("Pc-genes", "Nc-genes")))

stacked_label_df <- stacked_plot_df %>%
  group_by(celltype, age) %>%
  summarise(top = sum(percent), ase_genes = first(ase_genes), biotype = first(biotype), .groups = "drop") %>%
  mutate(label = paste0("italic(n) == ", ase_genes),
         y = top + 0.3)

pdf("figures/proportion_ase_genes_adult_aged.pdf", width = 8, height = 7)
ggplot(stacked_plot_df, aes(x = celltype, y = percent, fill = biotype)) +
  geom_col(width = 0.7) +
  geom_text(data = stacked_label_df,
            aes(y = y, label = label),
            angle = 90, vjust = 0.5, size = 3, parse = TRUE) +
  facet_wrap(~ age, nrow = 1) +                   # drop this line if you don’t want facets
  scale_fill_manual(values = c("Nc-genes" = "grey30", "Pc-genes" = "grey80")) +
  labs(x = NULL, y = "Fraction of ASE genes (%)", fill = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 11)
  )
dev.off()

# # Piechart of biotype - split by age
# proportion_biotype_df <- all.combined_autosome %>%
#     filter(biotype != 'mRNA' & biotype != 'transcript') %>%
#     group_by(celltype, age) %>%
#     summarise(
#         lncRNA = sum(biotype == "lnc_RNA") / n() * 100,
#         snRNA = sum(biotype == "snRNA") / n() * 100,
#         snoRNA = sum(biotype == "snoRNA") / n() * 100,
#         antisense_RNA = sum(biotype == "antisense_RNA") / n() * 100,
#         unprocessed_RNA = sum(biotype == "primary_transcript") / n() * 100,
#         .groups = "drop")

# pdf("figures/biotype_piechart_adult_aged.txt")
# ggplot(plot_df, aes(x = "", y = percent, fill = biotype)) +
#   geom_bar(stat = "identity", position = "fill") +
#   coord_polar("y") +
#   facet_wrap(~ age, nrow = 1) +
#   labs(x = NULL, y = "Fraction of ASE genes (%)", fill = NULL) +
#   theme_void()
# dev.off()

### Read in Allelome.LINK output
adult.files <- list.files('9w/Allelome.LINK', recursive=TRUE, pattern='_links_full_table.txt', full.names=TRUE)
adult.LINK <- lapply(adult.files, function(file) {
  read.delim(file)
})
names(adult.LINK) <- unlist(lapply(adult.files, function(x){
  tmp <- strsplit(x, '/')[[1]][[3]]
  strsplit(tmp, '_')[[1]][4]
}))

aged.files <- list.files('78w/Allelome.LINK', recursive=TRUE, pattern='_links_full_table.txt', full.names=TRUE)
aged.LINK <- lapply(aged.files, function(file) {
  read.delim(file)
})
names(aged.LINK) <- unlist(lapply(aged.files, function(x){
  tmp <- strsplit(x, '/')[[1]][[3]]
  strsplit(tmp, '_')[[1]][4]
}))

# Loop over adult and aged to annotate genes
adult.LINK <- lapply(adult.LINK, function(x) {
  dt_x <- data.table(x)
  dt_x[, gene_base := dt_gtf[.SD, on = .(transcript_id = name_base), gene_ids]]
  dt_x[, biotype_base := dt_gtf[.SD, on = .(transcript_id = name_base), biotype_ids]]
  dt_x[, gene_target := dt_gtf[.SD, on = .(transcript_id = name_target), gene_ids]]
  dt_x[, biotype_target := dt_gtf[.SD, on = .(transcript_id = name_target), biotype_ids]]
  as.data.frame(dt_x)
})

save(adult.LINK, file = "adult.Allelome.LINK.RData")
load("adult.Allelome.LINK.RData")

aged.LINK <- lapply(aged.LINK, function(x) {
  dt_x <- data.table(x)
  dt_x[, gene_base := dt_gtf[.SD, on = .(transcript_id = name_base), gene_ids]]
  dt_x[, gene_target := dt_gtf[.SD, on = .(transcript_id = name_target), gene_ids]]
  as.data.frame(dt_x)
})

save(aged.LINK, file = "aged.Allelome.LINK.RData")
load("aged.Allelome.LINK.RData")

#### Visualising the data ###

# Clustered bar chart of the proportion of enhancive and repressive links across cell types, split by age
combined.adult <- bind_rows(adult.LINK, .id = "celltype")
combined.adult$age <- 'adult'
combined.aged <- bind_rows(aged.LINK, .id = "celltype")
combined.aged$age <- 'aged'

combined.LINK <- bind_rows(combined.adult, combined.aged)
combined.LINK <- unique(combined.LINK[,-14])

pdf('figures/combined_links_barplot.pdf')
ggplot(combined.LINK, aes(x = celltype, fill = mechanism)) +
  # show absolute counts (total number of links) stacked by mechanism and faceted by age
  geom_bar(position = "stack") +
  labs(title = "",
     x = "",
     y = "No. of linkages") + 
  scale_fill_manual(values = c("enhancing" = "#8AAAA1", "repressing" = "#C36F6D")) +
  theme_minimal() +
  facet_wrap(~ age) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## Boxplot of linkage scores 
pdf("figures/combined_linkage_score_boxplot.pdf")
ggplot(combined.LINK, aes(x = celltype, y = linkage_score, fill = age)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = "",
         x = "",
         y = "Linkage Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

link_key <- function(df) paste(df$gene_base, df$gene_target, df$mechanism, sep="|")

# Venn diagram of links per cell type - adult and aged
adult_sets <- lapply(adult.LINK, function(x) unique(link_key(x)))
aged_sets  <- lapply(aged.LINK,  function(x) unique(link_key(x)))

plots <- lapply(names(adult_sets), function(ct){
  fit <- euler(c(
    'Adult' = length(adult_sets[[ct]]),
    'Aged'  = length(aged_sets[[ct]]),
    "Adult&Aged" = length(intersect(adult_sets[[ct]], aged_sets[[ct]]))
  ))
  plot(fit, quantities = TRUE, fill = age.colours, main=ct)

})
wrap_plots(plots, ncol = 4)
ggsave("figures/venn_small_multiples_per_celltype.pdf", width = 12, height = 8)


# Links which change during aging
summ_change <- lapply(names(adult_sets), function(ct){
  A <- unique(adult_sets[[ct]]); G <- unique(aged_sets[[ct]])
  tibble(
    celltype = ct,
    status   = c("Retained","Gained_in_Aged","Lost_in_Aged"),
    count    = c(length(intersect(A,G)),
                 length(setdiff(G,A)),
                 length(setdiff(A,G)))
  )
}) |> bind_rows() |>
  group_by(celltype) |>
  mutate(frac = count / sum(count))

ggplot(summ_change, aes(x=celltype, y=frac, fill=status)) +
  geom_col() +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Cell type", y="% of links", title="") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("Retained"="#94A2AB", "Gained_in_Aged"="#284456", "Lost_in_Aged"="#D9D9D6"))
ggsave("figures/retained_gained_lost_stackedbars.pdf", width=10, height=6)

jaccard <- lapply(names(adult_sets), function(ct){
  A <- adult_sets[[ct]]; G <- aged_sets[[ct]]
  tibble(celltype=ct, jaccard=length(intersect(A,G))/length(union(A,G)))
}) |> bind_rows()

ggplot(jaccard, aes(x=reorder(celltype, -jaccard), y=jaccard)) +
  geom_point(size=3) + ylim(0,1) +
  labs(x="Cell type", y="Jaccard (Adult vs Aged)", title="Overall concordance by cell type") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("figures/jaccard_per_celltype.pdf", width=8, height=5)


# Determine how many mechanisms switch
pair_key <- function(df) paste(df$gene_base, df$gene_target, sep="|")

# Bar plot
flip_tab <- lapply(names(adult.LINK), function(ct){
  a <- adult.LINK[[ct]] %>%
    mutate(pair = pair_key(.), mech = mechanism) %>%
    select(pair, mech)
  a <- unique(a)

  g <- aged.LINK[[ct]] %>%
    mutate(pair = pair_key(.), mech = mechanism) %>%
    select(pair, mech)
  g <- unique(g)

  merged <- inner_join(a, g, by = "pair", suffix = c("_adult","_aged"), relationship = "many-to-many")

  summarise(merged,
            celltype = ct,
            flips = sum(mech_adult != mech_aged, na.rm = TRUE),
            same  = sum(mech_adult == mech_aged, na.rm = TRUE))
}) %>% bind_rows()

ggplot(flip_tab, aes(x=reorder(celltype, -flips), y=flips)) +
  geom_col() +
  labs(x="Cell type", y="# mechanism flips", title="") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("figures/mechanism_flips_per_celltype.pdf", width=8, height=5)

# Check ID for coding or non-coding gene
is_lnc_refseq <- function(tx) grepl("^(NR|XR)_", tx)


flip_table <- lapply(names(adult.LINK), function(ct){
  a <- adult.LINK[[ct]] %>%
    mutate(pair = pair_key(.),
           base_tx = name_base,
           base_is_lnc = is_lnc_refseq(name_base)) %>%
    select(pair, mechanism_adult = mechanism, base_tx, base_is_lnc, gene_base, gene_target)
  a <- unique(a)

  g <- aged.LINK[[ct]] %>%
    mutate(pair = pair_key(.),
           base_tx = name_base,
           base_is_lnc = is_lnc_refseq(name_base)) %>%
    select(pair, mechanism_aged = mechanism, base_tx_aged = base_tx, base_is_lnc_aged = base_is_lnc)
  g <- unique(g)

  inner_join(a, g, by = "pair", relationship = "many-to-many") %>%
    mutate(flip = mechanism_adult != mechanism_aged,
           celltype = ct)
}) %>% bind_rows()

# flips only
flips_only <- flip_table %>% filter(flip)


tx_counts <- function(lst) {
  lapply(lst, function(df) df %>% select(transcript = name, total_reads))
}

adult_tx <- tx_counts(adult.Allelome)
aged_tx  <- tx_counts(aged.Allelome)

# Merge adult vs aged, compute log2FC per transcript for each cell type
tx_fc <- lapply(names(adult_tx), function(ct){
  a <- adult_tx[[ct]] %>% group_by(transcript) %>% summarise(count_adult = sum(total_reads), .groups="drop")
  g <- aged_tx[[ct]]  %>% group_by(transcript) %>% summarise(count_aged  = sum(total_reads), .groups="drop")
  full_join(a, g, by="transcript") %>%
    mutate(count_adult = replace_na(count_adult, 0),
           count_aged  = replace_na(count_aged,  0),
           # add a small prior to avoid log(0); tune if you like
           log2FC_tx = log2((count_aged + 1) / (count_adult + 1)),
           is_lnc = is_lnc_refseq(transcript)) %>%
    mutate(celltype = ct)
})
names(tx_fc) <- names(adult_tx)

subset(tx_fc[['ventCM1']], transcript == 'NM_009785.1')
tx_fc <- bind_rows(tx_fc)

### Read in fold change results
adult_aged_fc <- readRDS('adult_aged_log2FC_results.RDS')


### Understanding gain and loss
# Gained in aged (new link appears in aged) 
# Enhancing: base up & target up
# Repressing: base up & target down

# Lost in aged (link present in adult, gone in aged)
# Enhancing: base down & target down
# Repressing: base down & target up

# Tunable thresholds
lfc_up_thr   <-  0.5
lfc_down_thr <- -0.5
min_logCPM   <-  1.0

# Helper to fetch base/target FC + expression in each age
# adult_aged_fc[[ct]] expected to contain: gene, log2FC, logCPM_ref (adult), logCPM_tar (aged)
get_fc_row <- function(fc_table, gene_symbol) {
  out <- subset(fc_table, gene == gene_symbol)
  if (nrow(out) == 0) return(NULL)
  out[1, ]
}

# Evaluate one link under "Gained in aged" rules
score_gained <- function(link, fc_table) {
  parts <- strsplit(link, "\\|")[[1]]
  base <- parts[1]; target <- parts[2]; mech <- parts[3]

  base_fc   <- get_fc_row(fc_table, base)
  target_fc <- get_fc_row(fc_table, target)
  if (is.null(base_fc) || is.null(target_fc)) return(NA_real_)

  # Expression filters
  base_expr_aged   <- !is.na(base_fc$logCPM_tar)   && base_fc$logCPM_tar   >= min_logCPM
  target_expr_aged <- !is.na(target_fc$logCPM_tar) && target_fc$logCPM_tar >= min_logCPM
  if (!(base_expr_aged && target_expr_aged)) return(NA_real_)

  base_up    <- !is.na(base_fc$log2FC)   && base_fc$log2FC   >= lfc_up_thr
  target_up  <- !is.na(target_fc$log2FC) && target_fc$log2FC >= lfc_up_thr
  target_down<- !is.na(target_fc$log2FC) && target_fc$log2FC <= lfc_down_thr

  if (mech == "enhancing" && base_up && target_up) return(1)
  if (mech == "repressing" && base_up && target_down) return(1)
  return(0)
}

# Evaluate one link under "Lost in aged" rules
score_lost <- function(link, fc_table) {
  parts <- strsplit(link, "\\|")[[1]]
  base <- parts[1]; target <- parts[2]; mech <- parts[3]

  base_fc   <- get_fc_row(fc_table, base)
  target_fc <- get_fc_row(fc_table, target)
  if (is.null(base_fc) || is.null(target_fc)) return(NA_real_)

  # Expression filters: require base expressed in ADULT (where link existed)
  base_expr_adult   <- !is.na(base_fc$logCPM_ref)   && base_fc$logCPM_ref   >= min_logCPM
  target_expr_adult <- !is.na(target_fc$logCPM_ref) && target_fc$logCPM_ref >= min_logCPM
  if (!(base_expr_adult && target_expr_adult)) return(NA_real_)

  # For "lost", we expect base to fall in aged OR at least be not expressed in aged
  base_down   <- !is.na(base_fc$log2FC) && base_fc$log2FC <= lfc_down_thr
  base_no_aged<- is.na(base_fc$logCPM_tar) || base_fc$logCPM_tar < min_logCPM
  base_decline <- base_down || base_no_aged

  target_up    <- !is.na(target_fc$log2FC) && target_fc$log2FC >= lfc_up_thr
  target_down  <- !is.na(target_fc$log2FC) && target_fc$log2FC <= lfc_down_thr

  if (mech == "enhancing" && base_decline && target_down) return(1)
  if (mech == "repressing" && base_decline && target_up)  return(1)
  return(0)
}

# Wrapper to score a vector of links
test_concordance_expr <- function(links, fc_table, mode = c("gained","lost")) {
  mode <- match.arg(mode)
  sapply(links, function(x) {
    if (mode == "gained") score_gained(x, fc_table) else score_lost(x, fc_table)
  })
}

# Permutation: shuffle MECHANISMS only (keeps base/target identity & FCs)
permute_mechanisms <- function(links) {
  if (length(links) == 0) return(links)
  mechs <- sapply(links, function(x) strsplit(x, "\\|")[[1]][3])
  shuffled <- sample(mechs, length(mechs), replace = FALSE)
  sapply(seq_along(links), function(i) {
    parts <- strsplit(links[i], "\\|")[[1]]
    parts[3] <- shuffled[i]
    paste(parts, collapse="|")
  })
}

# Main per-celltype computation
set.seed(123)

mechanism_match_expr_perm <- lapply(names(adult_sets), function(ct){
    print(ct)
  fc_table <- adult_aged_fc[[ct]]

  Lost_in_Aged   <- setdiff(adult_sets[[ct]], aged_sets[[ct]])
  Gained_in_Aged <- setdiff(aged_sets[[ct]], adult_sets[[ct]])

  # Observed
  gained_obs_vec <- test_concordance_expr(Gained_in_Aged, fc_table, "gained")
  lost_obs_vec   <- test_concordance_expr(Lost_in_Aged,   fc_table, "lost")

  gained_obs <- mean(gained_obs_vec, na.rm = TRUE) * 100
  lost_obs   <- mean(lost_obs_vec,   na.rm = TRUE) * 100

  # Permutations
  nperm <- 1000
  gained_perm <- replicate(nperm, {
    links_p <- permute_mechanisms(Gained_in_Aged)
    mean(test_concordance_expr(links_p, fc_table, "gained"), na.rm = TRUE) * 100
  })
  lost_perm <- replicate(nperm, {
    links_p <- permute_mechanisms(Lost_in_Aged)
    mean(test_concordance_expr(links_p, fc_table, "lost"), na.rm = TRUE) * 100
  })

  data.frame(
    celltype = ct,
    Gained_obs = gained_obs,
    Lost_obs   = lost_obs,
    Gained_pval = mean(gained_perm >= gained_obs),
    Lost_pval   = mean(lost_perm   >= lost_obs),
    Gained_N = sum(!is.na(gained_obs_vec)),
    Lost_N   = sum(!is.na(lost_obs_vec))
  )
})

mechanism_match <- do.call(rbind, mechanism_match_expr_perm)
write.table(mechanism_match, file = "figures/mechanism_match_perm_results.txt", sep = "\t", row.names = FALSE)


# How many links are shared or unique to each cell type - adult and aged
adult_links <- lapply(names(adult.LINK), function(ct){
    unique(paste(adult.LINK[[ct]]$gene_base, adult.LINK[[ct]]$gene_target, adult.LINK[[ct]]$mechanism, sep='_'))
})
names(adult_links) <- names(adult.LINK)
adult_links_matrix <- fromList(adult_links)
pdf("figures/adult_links_upset.pdf", width = 8, height = 6, onefile = FALSE)
upset(adult_links_matrix, nsets = length(adult_links), order.by = "freq", 
      main.bar.color = "steelblue", sets.bar.color = "darkorange")
dev.off()

aged_links <- lapply(names(aged.LINK), function(ct){
    unique(paste(aged.LINK[[ct]]$gene_base, aged.LINK[[ct]]$gene_target, aged.LINK[[ct]]$mechanism, sep='_'))
})
names(aged_links) <- names(aged.LINK)
aged_links_matrix <- fromList(aged_links)
pdf("figures/aged_links_upset.pdf", width = 8, height = 6, onefile = FALSE)
upset(aged_links_matrix, nsets = length(aged_links), order.by = "freq", 
      main.bar.color = "steelblue", sets.bar.color = "darkorange")
dev.off()


### Subset for non-coding base transcripts
adult.LINK_nc <- lapply(adult.LINK, function(x) {
    x[grepl('NR|XR', x$name_base), ]
})
aged.LINK_nc <- lapply(aged.LINK, function(x) {
    x[grepl('NR|XR', x$name_base), ]
})

unique_links_nc <- lapply(names(aged.LINK_nc), function(ct) {
    setdiff(unique(paste(adult.LINK_nc[[ct]]$gene_base, adult.LINK_nc[[ct]]$gene_target, adult.LINK_nc[[ct]]$mechanism, sep='_')),
            unique(paste(aged.LINK_nc[[ct]]$gene_base, aged.LINK_nc[[ct]]$gene_target, aged.LINK_nc[[ct]]$mechanism, sep='_')))
})
names(unique_links_nc) <- names(aged.LINK_nc)
unique_links_nc_matrix <- fromList(unique_links_nc)

pdf("figures/discordant_links_nc_upset.pdf", width = 8, height = 6, onefile = FALSE)
upset(unique_links_nc_matrix, nsets = length(unique_links_nc), order.by = "freq", 
      main.bar.color = "steelblue", sets.bar.color = "darkorange")
dev.off()

### Read in fold change test 
adult_aged_fc <- readRDS('adult_aged_log2FC_results.RDS')

unique_links_nc_fc <- lapply(names(unique_links_nc), function(ct){
    base <- unlist(lapply(unique_links_nc[[ct]], function(x) strsplit(x, '_')[[1]][1]))
    target <- unlist(lapply(unique_links_nc[[ct]], function(x) strsplit(x, '_')[[1]][2]))

    base_fc <- subset(adult_aged_fc[[ct]], gene %in% base)
    base_fc$type <- 'base'

    target_fc <- subset(adult_aged_fc[[ct]], gene %in% target)
    target_fc$type <- 'target'

    rbind(base_fc, target_fc)
})

unique_links_nc_fc[[1]]

nrow(subset(adult.LINK_nc[[1]], gene_base %in% subset(unique_links_nc_fc[[1]], type == 'base')$gene))
nrow(subset(adult.LINK_nc[[1]], gene_target %in% subset(unique_links_nc_fc[[1]], type == 'target')$gene))

nrow(subset(aged.LINK_nc[[1]], gene_base %in% subset(unique_links_nc_fc[[1]], type == 'base')$gene))
nrow(subset(aged.LINK_nc[[1]], gene_target %in% subset(unique_links_nc_fc[[1]], type == 'target')$gene))
