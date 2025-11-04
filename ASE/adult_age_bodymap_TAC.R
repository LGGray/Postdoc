library(tidyverse)
library(data.table)
library(VennDiagram)
library(eulerr)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)

# Function to extract the link key - base|target|mechanism
link_key <- function(df) paste(df$name_base, df$name_target, df$mechanism, sep="|")

# Jaccard index function
jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  if (union == 0) {
    return(NA)  # Avoid division by zero
  }
  return(intersection / union)
}

sample_colours <- c("adult" = "#E69F00", "aged" = "#56B4E9", "TAC" = "#009E73", "Sham" = "#CC79A7")

# #### Match RefSeq transcript ID to gene name
# gtf <- read.delim('GRCm39/GCF_000001635.27_GRCm39_genomic.gtf', comment.char = "#", header = FALSE)
# gtf_transcripts <- subset(gtf, V3 == "transcript")

# attrs <- strsplit(as.character(gtf_transcripts$V9), "; ", fixed = TRUE)
# gene_ids <- vapply(attrs, function(x) {
#   if (length(x) >= 1) sub("^gene_id\\s*", "", x[1]) else NA_character_
# }, character(1))
# transcript_ids <- vapply(attrs, function(x) {
#   ti <- grep("^transcript_id", x, value = TRUE)
#   if (length(ti)) sub("^transcript_id\\s*", "", ti[1]) else NA_character_
# }, character(1))
# biotype_ids <- vapply(attrs, function(x) {
#   bt <- grep("^transcript_biotype", x, value = TRUE)
#   if (length(bt)) sub("^transcript_biotype\\s*", "", bt[1]) else NA_character_
# }, character(1))

# # make gtf data.tables
# dt_gtf <- data.table(gtf_transcripts)[, .(transcript_id = transcript_ids, gene_ids, biotype_ids)]

########################################
# Read in data from adult aged bodymap #
########################################
bodymap_dirs <- list.dirs('adult_aged_bodymap', recursive = FALSE, full.names = FALSE)
bodymap_dirs <- grep('_', bodymap_dirs, value = TRUE)

bodymap <- lapply(bodymap_dirs, function(dir) {
    tmp <- paste0('adult_aged_bodymap/', dir, '/no_predicted_65/no_predicted_65_links_full_table.txt')
    read.delim(tmp)
})
names(bodymap) <- bodymap_dirs

# # Loop over bodymap to annotate genes
# bodymap <- lapply(bodymap, function(x) {
#   dt_x <- data.table(x)
#   dt_x[, gene_base := dt_gtf[.SD, on = .(transcript_id = name_base), gene_ids]]
#   dt_x[, biotype_base := dt_gtf[.SD, on = .(transcript_id = name_base), biotype_ids]]
#   dt_x[, gene_target := dt_gtf[.SD, on = .(transcript_id = name_target), gene_ids]]
#   dt_x[, biotype_target := dt_gtf[.SD, on = .(transcript_id = name_target), biotype_ids]]
# #   # Filter out same-gene links (different isoforms of same gene)
# #   dt_x <- dt_x[gene_base != gene_target | is.na(gene_base) | is.na(gene_target)]
#   as.data.frame(dt_x)
# })


# Subset for noncoding-coding link
noncoding_coding <- read.delim('GRCm39/Refseq_coding_noncoding.txt', header = TRUE)
head(noncoding_coding)
coding <- subset(noncoding_coding, coding_noncoding == "coding")$name
noncoding <- subset(noncoding_coding, coding_noncoding == "noncoding")$name

bodymap <- lapply(bodymap, function(df) {
  subset(df, name_base %in% noncoding & name_target %in% coding)
})

save(bodymap, file = 'adult_aged_bodymap/bodymap_links_no_predicted_65.RData')
load('adult_aged_bodymap/bodymap_links_no_predicted_65.RData')


#######################################
# Compare link between adult and aged #
#######################################

# Colours for adult and aged
bodymap_link_keys <- lapply(bodymap, function(df) unique(link_key(df)))

links_mtx <- fromList(bodymap_link_keys)
rownames(links_mtx) <- unique(unlist(bodymap_link_keys))
links_mtx <- links_mtx[, c('Ao_9w', 'Br_9w', 'He_9w', 'Ki_9w', 'Li_9w', 'Lu_9w', 'Mu_9w', 'Sp_9w',
                           'Ao_78w', 'Br_78w', 'He_78w', 'Ki_78w', 'Li_78w', 'Lu_78w', 'Mu_78w', 'Sp_78w')]

# Heatmap of links across all adult and aged bodymap samples
pdf('LINKS_study/figures/heatmap_adult_aged_bodymap_links_no_predicted_65.pdf')
ann <- HeatmapAnnotation(
  Age = c(rep("Adult", 8), rep("Aged", 8)),
  col = list(Age = c("Adult" = "#E69F00", "Aged" = "#56B4E9"))
)

Heatmap(links_mtx,
        name = "Link",
        col = c("white", "red"),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = "",
        heatmap_legend_param = list(
          title = "Link Presence",
          at = c(0, 1),
          labels = c("Absent", "Present")
        ),
        top_annotation = ann
)
dev.off()

# UpSet plot of links across all adult and aged bodymap samples
pdf('LINKS_study/figures/UpSet_adult_aged_bodymap_links_no_predicted_65.pdf', onefile = FALSE, width = 10, height = 6)
upset(links_mtx,
      nsets = ncol(links_mtx),
      order.by = "freq",
      main.bar.color = "black",
      sets.bar.color = "black",
      mainbar.y.label = "Number of Links",
      sets.x.label = "Links per Sample"
)
dev.off()


# Venn diagram of adult vs aged links per tissue
tissues <- c('Ao', 'Br', 'He', 'Ki', 'Li', 'Lu', 'Mu', 'Sp')
tissue_name <- list(
  'Ao' = 'Aorta',
  'Br' = 'Brain',
  'He' = 'Heart',
  'Ki' = 'Kidney',
  'Li' = 'Liver',
  'Lu' = 'Lung',
  'Mu' = 'Muscle',
  'Sp' = 'Spleen'
)
plot_list <- lapply(tissues, function(tissue) {
  adult <- bodymap[[paste0(tissue, '_9w')]]
  aged <- bodymap[[paste0(tissue, '_78w')]]
  adult_sets <- unique(link_key(adult))
  aged_sets  <- unique(link_key(aged))
  fit <- euler(list(
      'Adult' = adult_sets,
      'Aged'  = aged_sets
    ))
  plot(fit, quantities = TRUE, fill = sample_colours[c("adult", "aged")], main=tissue_name[[tissue]])
})

pdf('LINKS_study/figures/venn_adult_aged_bodymap_links_no_predicted_65.pdf', width = 12, height = 6)
grid.arrange(grobs = plot_list, ncol = 4, nrow = 2)
dev.off()


# Compare the number of links gained/lost between adult and aged across all tissues
summ_change <- lapply(tissues, function(tissue) {
    adult <- bodymap[[paste0(tissue, '_9w')]]
    aged <- bodymap[[paste0(tissue, '_78w')]]
    A <- unique(link_key(adult))
    G <- unique(link_key(aged))
  tibble(
    organ = tissue_name[[tissue]],
    status   = c("Retained","Gained_in_Aged","Lost_in_Aged"),
    count    = c(length(intersect(A,G)),
                 length(setdiff(G,A)),
                 length(setdiff(A,G)))
  )
}) |> bind_rows() |>
  group_by(organ) |>
  mutate(frac = count / sum(count))

pdf('LINKS_study/figures/barplot_adult_aged_bodymap_link_changes_no_predicted_65.pdf', width = 8, height = 6)
ggplot(summ_change, aes(x=organ, y=frac, fill=status)) +
  geom_col() +
  scale_y_continuous(labels=scales::percent) +
  labs(x="", y="", title="") +
  theme_minimal() +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("Retained"="#94A2AB", "Gained_in_Aged"="#284456", "Lost_in_Aged"="#D9D9D6"))
dev.off()

# Count how many links switch mechanism between adult and aged
pair_key <- function(df) paste(df$name_base, df$name_target, sep="|")

mech_switch_by_tissue <- lapply(tissues, function(tissue) {
  adult <- bodymap[[paste0(tissue, '_9w')]]
  aged  <- bodymap[[paste0(tissue, '_78w')]]

  # mechanisms per base|target pair (within each age)
  adult_map <- lapply(split(adult$mechanism, pair_key(adult)), function(x) unique(na.omit(x)))
  aged_map  <- lapply(split(aged$mechanism,  pair_key(aged)),  function(x) unique(na.omit(x)))

  common_pairs <- intersect(names(adult_map), names(aged_map))

  if (length(common_pairs) == 0) {
    return(tibble(
      organ = tissue_name[[tissue]],
      n_common_pairs = 0,
      n_switch_any = 0,
      n_switch_strict = 0,
      frac_switch_any = NA_real_,
      frac_switch_strict = NA_real_
    ))
  }

  sw <- vapply(common_pairs, function(p) {
    a <- adult_map[[p]]
    g <- aged_map[[p]]
    any_switch <- ("enhancing" %in% a && "repressing" %in% g) ||
                  ("repressing" %in% a && "enhancing" %in% g)
    strict_switch <- (setequal(a, "enhancing") && setequal(g, "repressing")) ||
                     (setequal(a, "repressing") && setequal(g, "enhancing"))
    c(any = any_switch, strict = strict_switch)
  }, FUN.VALUE = c(any = FALSE, strict = FALSE))

  n_any <- sum(sw["any", ])
  n_strict <- sum(sw["strict", ])

  tibble(
    organ = tissue_name[[tissue]],
    n_common_pairs = length(common_pairs),
    n_switch_any = n_any,
    n_switch_strict = n_strict,
    frac_switch_any = n_any / length(common_pairs),
    frac_switch_strict = n_strict / length(common_pairs)
  )
}) |> bind_rows()

# Identify links which do switch
mech_switch_links <- lapply(tissues, function(tissue) {
  adult <- bodymap[[paste0(tissue, '_9w')]]
  aged  <- bodymap[[paste0(tissue, '_78w')]]

  # mechanisms per base|target pair (within each age)
  adult_map <- lapply(split(adult$mechanism, pair_key(adult)), function(x) unique(na.omit(x)))
  aged_map  <- lapply(split(aged$mechanism,  pair_key(aged)),  function(x) unique(na.omit(x)))

  common_pairs <- intersect(names(adult_map), names(aged_map))

  sw_links <- lapply(common_pairs, function(p) {
    a <- adult_map[[p]]
    g <- aged_map[[p]]
    any_switch <- ("enhancing" %in% a && "repressing" %in% g) ||
                  ("repressing" %in% a && "enhancing" %in% g)
    if (any_switch) {
      return(p)
    } else {
      return(NULL)
    }
  })
  sw_links <- unlist(sw_links)
  names(sw_links) <- NULL
  return(sw_links)
})
names(mech_switch_links) <- tissues

##################################################################
# Understanding changes in link mechanism between adult and aged #
##################################################################

deseq <- read.delim('adult_aged_bodymap/DESeq2/He_9w_vs_78w.txt', header = TRUE)
deseq$ensembl_gene_id_version <- gsub('\\..+', '', deseq$ensembl_gene_id_version)

# Plot a venn diagram of adult and aged heart
pdf('LINKS_study/figures/venn_adult_aged_bodymap_He_links_no_predicted_65.pdf')
adult_links <- link_key(bodymap[['He_9w']])
aged_links  <- link_key(bodymap[['He_78w']])
fit <- euler(list(
    'Adult' = adult_links,
    'Aged'  = aged_links
  ))
plot(fit, quantities = TRUE, fill = sample_colours[c("adult", "aged")], main='Heart')
dev.off()

adult_specific <- setdiff(link_key(bodymap[['He_9w']]), link_key(bodymap[['He_78w']]))
aged_specific  <- setdiff(link_key(bodymap[['He_78w']]), link_key(bodymap[['He_9w']]))

adult_specific <- lapply(adult_specific, function(link) {
  data.frame(base = strsplit(link, "\\|")[[1]][1],
              target = strsplit(link, "\\|")[[1]][2])
})
adult_specific <- bind_rows(adult_specific)

library(biomaRt)
# convert gene name to ensembl id
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
base_map <- getBM(attributes = c('mgi_symbol', 'ensembl_gene_id'),
                  filters = 'mgi_symbol',
                  values = adult_specific$base,
                  mart = mart)

target_map <- getBM(attributes = c('mgi_symbol', 'ensembl_gene_id'),
                     filters = 'mgi_symbol',
                     values = adult_specific$target,
                     mart = mart)
base_deg <- subset(deseq, ensembl_gene_id_version %in% base_map$ensembl_gene_id & padj < 0.05)
target_deg <- subset(deseq, ensembl_gene_id_version %in% target_map$ensembl_gene_id & padj < 0.05)
adult_deg <- rbind(
  data.frame(gene = base_deg$gene_name,
             log2FoldChange = base_deg$log2FoldChange,
             padj = base_deg$padj,
             role = 'base'),
  data.frame(gene = target_deg$gene_name,
             log2FoldChange = target_deg$log2FoldChange,
             padj = target_deg$padj,
             role = 'target')
)

aged_specific <- lapply(aged_specific, function(link) {
  data.frame(base = strsplit(link, "\\|")[[1]][1],
              target = strsplit(link, "\\|")[[1]][2])
})
aged_specific <- bind_rows(aged_specific)
# convert gene name to ensembl id
base_map <- getBM(attributes = c('mgi_symbol', 'ensembl_gene_id'),
                  filters = 'mgi_symbol',
                  values = aged_specific$base,
                  mart = mart)  
target_map <- getBM(attributes = c('mgi_symbol', 'ensembl_gene_id'),
                     filters = 'mgi_symbol',
                     values = aged_specific$target,
                     mart = mart)
base_deg <- subset(deseq, ensembl_gene_id_version %in% base_map$ensembl_gene_id & padj < 0.05)
target_deg <- subset(deseq, ensembl_gene_id_version %in% target_map$ensembl_gene_id & padj < 0.05)
aged_deg <- rbind(
  data.frame(gene = base_deg$gene_name,
             log2FoldChange = base_deg$log2FoldChange,
             padj = base_deg$padj,
             role = 'base'),
  data.frame(gene = target_deg$gene_name,
             log2FoldChange = target_deg$log2FoldChange,
             padj = target_deg$padj,
             role = 'target')
)

deg_result <- rbind(adult_deg, data.frame(gene = base_deg$gene_name,
             log2FoldChange = base_deg$log2FoldChange,
             padj = base_deg$padj,
             role = 'base'))
write.table(deg_result, file = 'LINKS_study/DE_genes_in_adult_aged_bodymap_He_links_no_predicted_65.txt', sep = "\t", row.names = FALSE, quote = FALSE)

selected_genes <- unique(c(adult_deg$gene, 'H19'))

subset(deseq, gene_name %in% selected_genes)[,c('gene_name', 'log2FoldChange', 'padj')]

# Plot a volcano plot highlighting DE genes involved in lost/gained links
deseq$highlight <- ifelse(deseq$gene_name %in% selected_genes, deseq$gene_name, NA)
pdf('LINKS_study/figures/volcano_adult_aged_bodymap_He_links_no_predicted_65.pdf', width = 8, height = 6)
ggplot(deseq, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(alpha=0.4) +
  geom_point(data=subset(deseq, !is.na(highlight)), color='red') +
  # geom_text(data=subset(deseq, !is.na(highlight)), aes(label=highlight), vjust=1.5, color='black', position=position_jitter(width=0.2, height=0)) +
  theme_minimal() +
  labs(title='Heart: Adult vs Aged', x='Log2 Fold Change', y='-Log10 Adjusted P-value') +
  theme(legend.position = "none")
dev.off()

#############################################
# Are links changing due to switches in ASE #
#############################################
adult_allelome <- read.delim('adult_aged_bodymap/He_9w/annotation_no_predicted_locus_table_reps.txt', header = TRUE)
aged_allelome  <- read.delim('adult_aged_bodymap/He_78w/annotation_no_predicted_locus_table_reps.txt', header = TRUE)

adult_ase <- apply(adult_specific, 1, function(row) {
  base <- subset(adult_allelome, name == row[['base']])
  base <- base[, c('name','allelic_ratio'), drop = FALSE]
  if (nrow(base) == 0) base <- data.frame(name = row[['base']], allelic_ratio = NA_real_)
  base$gene_type <- 'base'
  target <- subset(adult_allelome, name == row[['target']])
  target <- target[, c('name','allelic_ratio'), drop = FALSE]
  if (nrow(target) == 0) target <- data.frame(name = row[['target']], allelic_ratio = NA_real_)
  target$gene_type <- 'target'
  rbind(base, target)[,c('name', 'allelic_ratio', 'gene_type')]
})
adult_ase <- bind_rows(adult_ase)

aged_ase <- apply(adult_specific, 1, function(row) {
  base <- subset(aged_allelome, name == row[['base']])
  base <- base[, c('name','allelic_ratio'), drop = FALSE]
  if (nrow(base) == 0) base <- data.frame(name = row[['base']], allelic_ratio = NA_real_)
  base$gene_type <- 'base'

  target <- subset(aged_allelome, name == row[['target']])
  target <- target[, c('name','allelic_ratio'), drop = FALSE]
  if (nrow(target) == 0) target <- data.frame(name = row[['target']], allelic_ratio = NA_real_)
  target$gene_type <- 'target'

  rbind(base, target)  # now both have exactly the same columns
})
aged_ase <- bind_rows(aged_ase)

combined_adult <- merge(adult_ase, aged_ase, by = "name", suffixes = c("_adult", "_aged"))

pdf('LINKS_study/figures/scatter_adult_aged_bodymap_He_ASE_adult_specific_links_no_predicted_65.pdf')
ggplot(unique(combined_adult), aes(x=allelic_ratio_adult, y=allelic_ratio_aged, color=gene_type_adult)) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = c(0.35, 0.65), linetype = 2) +
  geom_hline(yintercept = c(0.35, 0.65), linetype = 2) +
  theme_minimal() +
  labs(title='ASE in Adult-specific Links: Adult vs Aged', x='Allelic Ratio (Adult)', y='Allelic Ratio (Aged)') +
  geom_abline(slope=1, intercept=0, linetype="dashed", color = "red")
dev.off()

cor(combined_adult$allelic_ratio_adult, combined_adult$allelic_ratio_aged, method = "spearman")


adult_ase <- apply(aged_specific, 1, function(row) {
  base <- subset(adult_allelome, name == row[['base']])
  base <- base[, c('name','allelic_ratio'), drop = FALSE]
  if (nrow(base) == 0) base <- data.frame(name = row[['base']], allelic_ratio = NA_real_)
  base$gene_type <- 'base'
  target <- subset(adult_allelome, name == row[['target']])
  target <- target[, c('name','allelic_ratio'), drop = FALSE]
  if (nrow(target) == 0) target <- data.frame(name = row[['target']], allelic_ratio = NA_real_)
  target$gene_type <- 'target'
  rbind(base, target)[,c('name', 'allelic_ratio', 'gene_type')]
})
adult_ase <- bind_rows(adult_ase)
aged_ase <- apply(aged_specific, 1, function(row) {
  base <- subset(aged_allelome, name == row[['base']])
  base <- base[, c('name','allelic_ratio'), drop = FALSE]
  if (nrow(base) == 0) base <- data.frame(name = row[['base']], allelic_ratio = NA_real_)
  base$gene_type <- 'base'
  target <- subset(aged_allelome, name == row[['target']])
  target <- target[, c('name','allelic_ratio'), drop = FALSE]
  if (nrow(target) == 0) target <- data.frame(name = row[['target']], allelic_ratio = NA_real_)   
  target$gene_type <- 'target'
  rbind(base, target)[,c('name', 'allelic_ratio', 'gene_type')]
})
aged_ase <- bind_rows(aged_ase) 

combined_aged <- merge(adult_ase, aged_ase, by = "name", suffixes = c("_adult", "_aged"))
pdf('LINKS_study/figures/scatter_adult_aged_bodymap_He_ASE_aged_specific_links_no_predicted_65.pdf')
ggplot(unique(combined_aged), aes(x=allelic_ratio_adult, y=allelic_ratio_aged, color=gene_type_aged)) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = c(0.35, 0.65), linetype = 2) +
  geom_hline(yintercept = c(0.35, 0.65), linetype = 2) +
  theme_minimal() +
  labs(title='ASE in Aged-specific Links: Adult vs Aged', x='Allelic Ratio (Adult)', y='Allelic Ratio (Aged)') +
  geom_abline(slope=1, intercept=0, linetype="dashed", color = "red")
dev.off()


############################################################################
# Overlap with snRNAseq cell-type specific links from adult and aged heart #
############################################################################
adult.files <- list.files('adult_aged_heart_snRNAseq/9w/Allelome.LINK', pattern = 'no_predicted_65_links_full_table.txt', full.names = TRUE, recursive = TRUE)
adult.files <- grep('2025_10_31', adult.files, value = TRUE)
adult_link <- lapply(adult.files, function(f) {
    tmp <- read.delim(f)
    return(tmp)
})
names(adult_link) <- lapply(adult.files, function(f) {
    ct <- strsplit(f, '/')[[1]][5]
    ct <- gsub('_no_predicted_65', '', ct)
    return(ct)
})

# Filter for noncoding-coding links only
adult_link <- lapply(adult_link, function(df) {
  subset(df, name_base %in% noncoding & name_target %in% coding)
})

save(adult_link, file = 'adult_aged_heart_snRNAseq/adult.Allelome.LINK_no_predicted_65.RData')
load('adult_aged_heart_snRNAseq/adult.Allelome.LINK_no_predicted_65.RData')

aged.files <- list.files('adult_aged_heart_snRNAseq/78w/Allelome.LINK', pattern = 'no_predicted_65_links_full_table.txt', full.names = TRUE, recursive = TRUE)
aged.files <- grep('2025_10_31', aged.files, value = TRUE)
aged_link <- lapply(aged.files, function(f) {
    tmp <- read.delim(f)
    return(tmp)
})
names(aged_link) <- lapply(aged.files, function(f) {
    ct <- strsplit(f, '/')[[1]][5]
    ct <- gsub('_no_predicted_65', '', ct)
    return(ct)
})

# Filter for noncoding-coding links only
aged_link <- lapply(aged_link, function(df) {
  subset(df, name_base %in% noncoding & name_target %in% coding)
})

save(aged_link, file = 'adult_aged_heart_snRNAseq/aged.Allelome.LINK_no_predicted_65.RData')
load('adult_aged_heart_snRNAseq/aged.Allelome.LINK_no_predicted_65.RData')


# Visualise cell type specific links
adult_mtx <- fromList(lapply(adult_link, function(x) { unique(link_key(x)) }))
rownames(adult_mtx) <- unique(unlist(lapply(adult_link, function(x) { unique(link_key(x)) })))
colnames(adult_mtx) <- paste0(colnames(adult_mtx), '_adult')
aged_mtx <- fromList(lapply(aged_link, function(x) { unique(link_key(x)) }))
rownames(aged_mtx) <- unique(unlist(lapply(aged_link, function(x) { unique(link_key(x)) })))
colnames(aged_mtx) <- paste0(colnames(aged_mtx), '_aged')
combined_mtx <- merge(adult_mtx, aged_mtx, by = "row.names", all = TRUE)

rownames(combined_mtx) <- combined_mtx$Row.names
combined_mtx <- combined_mtx[, -1]
combined_mtx[is.na(combined_mtx)] <- 0

pdf('LINKS_study/figures/heatmap_adult_aged_heart_snRNAseq_links_no_predicted_65.pdf')
Heatmap(as.matrix(combined_mtx),
        name = "Link",
        col = c("white", "red"),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = "",
        heatmap_legend_param = list(
          title = "Link Presence",
          at = c(0, 1),
          labels = c("Absent", "Present")
        )
)
dev.off()

pdf('LINKS_study/figures/UpSet_adult_aged_heart_snRNAseq_links_no_predicted_65.pdf', onefile = FALSE, width = 10, height = 10)
upset(combined_mtx,
      nsets = ncol(combined_mtx),
      order.by = "freq",
      main.bar.color = "black",
      sets.bar.color = "black",
      mainbar.y.label = "Number of Links",
      sets.x.label = "Links per Sample"
)
dev.off()

# Venn diagram of adult vs aged links per cell type
plot_list <- lapply(names(adult_link), function(ct) {
  adult_ct <- adult_link[[ct]]
  aged_ct  <- aged_link[[ct]]
  adult_sets <- unique(link_key(adult_ct))
  aged_sets  <- unique(link_key(aged_ct))
  fit <- euler(list(
      'Adult' = adult_sets,
      'Aged'  = aged_sets
    ))
  plot(fit, quantities = TRUE, fill = age.colours, main=ct)
})
pdf('LINKS_study/figures/venn_adult_aged_heart_snRNAseq_links_no_predicted_65.pdf', width = 12, height = 6)
grid.arrange(grobs = plot_list, ncol = 6, nrow = 2)
dev.off()

# Plot heatmap of the overlap between adult and aged heart cell-type specific links
common_links <- lapply(names(adult_link), function(ct) {
    adult_ct <- adult_link[[ct]]
    aged_ct  <- aged_link[[ct]]
    adult_sets <- unique(link_key(adult_ct))
    aged_sets  <- unique(link_key(aged_ct))
    intersect(adult_sets, aged_sets)
})
names(common_links) <- names(adult_link)
common_links_mtx <- fromList(common_links)
rownames(common_links_mtx) <- unique(unlist(common_links))

pdf('LINKS_study/figures/heatmap_adult_aged_heart_snRNAseq_common_links.pdf')
Heatmap(as.matrix(common_links_mtx),
        name = "Link",
        col = c("white", "red"),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = "",
        heatmap_legend_param = list(
          title = "Link Presence",
          at = c(0, 1),
          labels = c("Absent", "Present")
        )
)
dev.off()

pdf('LINKS_study/figures/UpSet_adult_aged_heart_snRNAseq_common_links.pdf', onefile = FALSE, width = 10, height = 6)
upset(common_links_mtx,
      nsets = ncol(common_links_mtx),
      order.by = "freq",
      main.bar.color = "black",
      sets.bar.color = "#94A2AB",
      mainbar.y.label = "Number of Links",
      sets.x.label = "Links per Sample"
)
dev.off()

# Overlap with adult heart links and calculate Jaccard index
adult_He <- link_key(bodymap[['He_9w']])
adult_jaccard <- lapply(adult_link, function(ct){
    ct_links <- unique(link_key(ct))
    jaccard_index(adult_He, ct_links)
    # print(length(intersect(adult_He, ct_links)))
})

aged_He <- link_key(bodymap[['He_78w']])
aged_jaccard <- lapply(aged_link, function(ct){
    ct_links <- unique(link_key(ct))
    jaccard_index(aged_He, ct_links)
    # print(length(intersect(aged_He, ct_links)))
})
aged_jaccard['Epicard'] <- 0

aged_jaccard <- aged_jaccard[names(adult_jaccard)]

plot_data <- data.frame(
  Cell_Type = names(adult_jaccard),
  Adult = unlist(adult_jaccard),
  Aged = unlist(aged_jaccard)
)

pdf('LINKS_study/figures/heatmap_jaccard_adult_aged_heart_snRNAseq_links_no_predicted_65.pdf')
Heatmap(as.matrix(plot_data[, -1]),
        name = "Jaccard Index",
        col = colorRamp2(c(0, max(plot_data[, -1], na.rm=TRUE)), c("white", "red")),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = "",
        heatmap_legend_param = list(
          title = "Jaccard Index",
          at = c(0, max(plot_data[, -1], na.rm=TRUE)/2, max(plot_data[, -1], na.rm=TRUE)),
          labels = c("0", round(max(plot_data[, -1], na.rm=TRUE)/2, 2), round(max(plot_data[, -1], na.rm=TRUE), 2))
        )
)
dev.off()


##################################
# Read in data from TAC and Sham #
##################################
TAC_dirs <- list.dirs('F1_TAC_Sarah', recursive = FALSE, full.names = FALSE)
TAC_dirs <- grep('He', TAC_dirs, value = TRUE)

TAC_SHAM <- lapply(TAC_dirs, function(dir) {
    tmp <- paste0('F1_TAC_Sarah/', dir, '/no_predicted_65/no_predicted_65_links_full_table.txt')
    read.delim(tmp)
})
names(TAC_SHAM) <- TAC_dirs

TAC_SHAM <- lapply(TAC_SHAM, function(df) {
  subset(df, name_base %in% noncoding & name_target %in% coding)
})

save(TAC_SHAM, file = 'F1_TAC_Sarah/TAC_SHAM_links_no_predicted_65.RData')
load('F1_TAC_Sarah/TAC_SHAM_links_no_predicted_65.RData')

# # Loop over bodymap to annotate genes
# TAC_SHAM <- lapply(TAC_SHAM, function(x) {
#   dt_x <- data.table(x)
#   dt_x[, gene_base := dt_gtf[.SD, on = .(transcript_id = name_base), gene_ids]]
#   dt_x[, biotype_base := dt_gtf[.SD, on = .(transcript_id = name_base), biotype_ids]]
#   dt_x[, gene_target := dt_gtf[.SD, on = .(transcript_id = name_target), gene_ids]]
#   dt_x[, biotype_target := dt_gtf[.SD, on = .(transcript_id = name_target), biotype_ids]]
# #   # Filter out same-gene links (different isoforms of same gene)
# #   dt_x <- dt_x[gene_base != gene_target | is.na(gene_base) | is.na(gene_target)]
#   as.data.frame(dt_x)
# })


# Venn diagram of TAC vs Sham links in heart split by enhancing and repressive in the one plot
tac_heart_links <- TAC_SHAM[['He_TAC28d_RNA_XistBxC_-+_XX']]
sham_heart_links <- TAC_SHAM[['He_TACSham28d_RNA_XistBxC_-+_XX']]
tac_enhancing <- unique(link_key(subset(tac_heart_links, mechanism == "enhancing")))
tac_repressive <- unique(link_key(subset(tac_heart_links, mechanism == "repressing")))
sham_enhancing <- unique(link_key(subset(sham_heart_links, mechanism == "enhancing")))
sham_repressive <- unique(link_key(subset(sham_heart_links, mechanism == "repressing")))



# Capture the two euler plots as grobs and arrange them side-by-side in one PDF
fit_enhancing <- euler(list(
  TAC  = tac_enhancing,
  Sham = sham_enhancing
))
fit_enhancing <- plot(fit_enhancing, quantities = TRUE, fill = sample_colours[c("TAC", "Sham")], main='Enhancing Links')
fit_repressive <- euler(list(
  TAC  = tac_repressive,
  Sham = sham_repressive
))
fit_repressive <- plot(fit_repressive, quantities=TRUE, fill = sample_colours[c("TAC", "Sham")], main='Repressive Links')
plot_list <- list(fit_enhancing, fit_repressive)

pdf('LINKS_study/figures/venn_TAC_sham_links_no_predicted_65.pdf')
grid.arrange(grobs = plot_list, ncol = 2)
dev.off()

tac_unique_enhancing <- setdiff(tac_enhancing, sham_enhancing)
tac_unique_repressive <- setdiff(tac_repressive, sham_repressive)
sham_unique_enhancing <- setdiff(sham_enhancing, tac_enhancing)
sham_unique_repressive <- setdiff(sham_repressive, tac_repressive)

##################################################
# Compare link between adult, aged and TAC heart #
##################################################

adult_heart_repressive <- unique(link_key(subset(bodymap[['He_9w']], mechanism == "repressing")))
adult_heart_enhancing <- unique(link_key(subset(bodymap[['He_9w']], mechanism == "enhancing")))
aged_heart_repressive  <- unique(link_key(subset(bodymap[['He_78w']], mechanism == "repressing")))
aged_heart_enhancing  <- unique(link_key(subset(bodymap[['He_78w']], mechanism == "enhancing")))

adult_aged_tac_sham_enhancing <- euler(list(
  Adult = adult_heart_enhancing,
  Aged  = aged_heart_enhancing,
  TAC   = tac_unique_enhancing,
  Sham  = sham_unique_enhancing
))
adult_aged_tac_sham_enhancing <- plot(adult_aged_tac_sham_enhancing, quantities = TRUE, fill = sample_colours, main='Enhancing Links in Heart')

adult_aged_tac_sham_repressive <- euler(list(
  Adult = adult_heart_repressive,
  Aged  = aged_heart_repressive,
  TAC   = tac_unique_repressive,
  Sham  = sham_unique_repressive
))
adult_aged_tac_sham_repressive <- plot(adult_aged_tac_sham_repressive, quantities=TRUE, fill=sample_colours, main='Repressive Links in Heart')
plot_list <- list(adult_aged_tac_sham_enhancing, adult_aged_tac_sham_repressive)
pdf('LINKS_study/figures/venn_adult_aged_TAC_Sham_links_no_predicted_65.pdf')
grid.arrange(grobs = plot_list, nrow = 2)
dev.off()

enhancing_list <- list(
  adult_heart_enhancing = adult_heart_enhancing,
  aged_heart_enhancing  = aged_heart_enhancing,
  tac_unique_enhancing  = tac_unique_enhancing,
  sham_unique_enhancing = sham_unique_enhancing
)
enhancing_mtx <- fromList(enhancing_list)

pdf('LINKS_study/figures/upset_adult_aged_TAC_Sham_enhancing_links_no_predicted_65.pdf', onefile = FALSE)
upset(enhancing_mtx,
      nsets = ncol(enhancing_mtx),
      order.by = "freq",
      main.bar.color = "black",
      sets.bar.color = "black",
      mainbar.y.label = "Number of Links",
      sets.x.label = "Links per Sample"
)
dev.off()


repressive_list <- list(
  adult_heart_repressive = adult_heart_repressive,
  aged_heart_repressive  = aged_heart_repressive,
  tac_unique_repressive  = tac_unique_repressive,
  sham_unique_repressive = sham_unique_repressive
)
repressive_mtx <- fromList(repressive_list)
pdf('LINKS_study/figures/upset_adult_aged_TAC_Sham_repressive_links_no_predicted_65.pdf', onefile = FALSE)
upset(repressive_mtx,
      nsets = ncol(repressive_mtx),
      order.by = "freq",
      main.bar.color = "black",
      sets.bar.color = "black",
      mainbar.y.label = "Number of Links",
      sets.x.label = "Links per Sample"
)
dev.off()


# Aged ∩ TAC but not Adult
aged_tac_only_rep <- setdiff(intersect(aged_heart_repressive, tac_unique_repressive), adult_heart_repressive)
aged_tac_only_enh <- setdiff(intersect(aged_heart_enhancing, tac_unique_enhancing), adult_heart_enhancing)
# Adult ∩ Aged but not TAC
adult_aged_only_rep <- setdiff(intersect(adult_heart_repressive, aged_heart_repressive), tac_unique_repressive)
adult_aged_only_enh <- setdiff(intersect(adult_heart_enhancing, aged_heart_enhancing), tac_unique_enhancing)
# Adult ∩ TAC but not Aged
adult_tac_only_rep <- setdiff(intersect(adult_heart_repressive, tac_unique_repressive), aged_heart_repressive)
adult_tac_only_enh <- setdiff(intersect(adult_heart_enhancing, tac_unique_enhancing), aged_heart_enhancing)

overlap_lists <- list(
  aged_tac_only_rep = aged_tac_only_rep,
  aged_tac_only_enh = aged_tac_only_enh,
  adult_aged_only_rep = adult_aged_only_rep,
  adult_aged_only_enh = adult_aged_only_enh,
  adult_tac_only_rep = adult_tac_only_rep,
  adult_tac_only_enh = adult_tac_only_enh
)
save(overlap_lists, file = 'LINKS_study/adult_aged_TAC_link_overlaps_no_predicted_65.RData')
load('LINKS_study/adult_aged_TAC_link_overlaps_no_predicted_65.RData')

##################################################
# Compare link between adult, aged and Sham heart #
adult_aged_sham_enhancing <- euler(list(
  Adult = adult_heart_enhancing,
  Aged  = aged_heart_enhancing,
  Sham  = sham_unique_enhancing
))
adult_aged_sham_enhancing <- plot(adult_aged_sham_enhancing, quantities = TRUE, fill = sample_colours, main='Enhancing Links in Heart') 
adult_aged_sham_repressive <- euler(list(
  Adult = adult_heart_repressive,
  Aged  = aged_heart_repressive,
  Sham  = sham_unique_repressive
))
adult_aged_sham_repressive <- plot(adult_aged_sham_repressive, quantities=TRUE, fill=sample_colours, main='Repressive Links in Heart')
plot_list <- list(adult_aged_sham_enhancing, adult_aged_sham_repressive)
pdf('LINKS_study/figures/venn_adult_aged_Sham_links_no_predicted_65.pdf')
grid.arrange(grobs = plot_list, nrow = 2)
dev.off()

################################
# Analysis of repressive links #
################################

adult_TAC <- intersect(adult_repressive, tac_repressive)
aged_TAC <- intersect(aged_repressive, tac_repressive)

# Overlap with cell-type specific links from adult and aged heart snRNA-seq
load('adult_aged_heart_snRNAseq/adult.Allelome.LINK.RData')
load('adult_aged_heart_snRNAseq/aged.Allelome.LINK.RData')

snRNAseq_adult <- lapply(names(adult.LINK), function(ct) {
    ct_links <- unique(link_key(subset(adult.LINK[[ct]], mechanism == "repressing")))
    if (length(ct_links) == 0) {
        return(data.frame())
    }
    tmp <- data.frame(celltype=ct, link=ct_links)
    subset(tmp, link %in% adult_TAC)
})
names(snRNAseq_adult) <- names(adult.LINK)
snRNAseq_adult_df <- bind_rows(snRNAseq_adult, .id = "celltype")
write.table(snRNAseq_adult_df, file='LINKS_study/snRNAseq_adult_TAC_links.txt', sep='\t', quote=FALSE, row.names=FALSE)


snRNAseq_aged <- lapply(names(aged.LINK), function(ct) {
    ct_links <- unique(link_key(subset(aged.LINK[[ct]], mechanism == "repressing")))
    if (length(ct_links) == 0) {
        return(data.frame())
    }
    tmp <- data.frame(celltype=ct, link=ct_links)
    subset(tmp, link %in% aged_TAC)
})
names(snRNAseq_aged) <- names(aged.LINK)
snRNAseq_aged_df <- bind_rows(snRNAseq_aged, .id = "celltype")
write.table(snRNAseq_aged_df, file='LINKS_study/snRNAseq_aged_TAC_links.txt', sep='\t', quote=FALSE, row.names=FALSE)


############################
# Map mouse genes to human #
############################
library(biomaRt)
library(dplyr)
library(tidyr)
library(readr)

# marts
mmart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")


adult <- read.delim("snRNAseq_adult_TAC_links.txt")
aged  <- read.delim("snRNAseq_aged_TAC_links.txt")

tac_links_genes <- bind_rows(adult, aged) %>%
  separate(link, into = c("gene_base","gene_target","mechanism"), sep = "\\|") %>%
  select(gene_base, gene_target) %>%
  distinct()

mouse_syms <- unique(c(tac_links_genes$gene_base, tac_links_genes$gene_target))

attrs_gene <- c("ensembl_gene_id","mgi_symbol")
gene_tbl <- getBM(
  attributes = attrs_gene,
  filters = "mgi_symbol",
  values  = mouse_syms,
  mart    = mmart
)

attrs_hom <- c("ensembl_gene_id",
               "hsapiens_homolog_ensembl_gene",
               "hsapiens_homolog_associated_gene_name",
               "hsapiens_homolog_orthology_type",
               "hsapiens_homolog_orthology_confidence",
               "hsapiens_homolog_perc_id")
hom_tbl <- getBM(
  attributes = attrs_hom,
  filters = "mgi_symbol",
  values  = mouse_syms,
  mart    = mmart
)

res <- gene_tbl %>%
  left_join(hom_tbl, by = "ensembl_gene_id") %>%
  relocate(mgi_symbol, ensembl_gene_id)

# Genes to follow up in human
human_genes <- unique(res$hsapiens_homolog_associated_gene_name)[-1]

#################################
# Read in Gtex data and overlap #
#################################

gtex_link <- read.delim('LINKS_study/GTEx_ASE_links.txt')

gtex_link_flt <- subset(gtex_link, SEX == 'Female' & Tissue %in% c('Heart - Left Ventricle', 'Heart - Atrial Appendage') & Mechanism == 'repressing')

human_links <- toupper(unlist(overlap_lists))
human_links <- gsub('REPRESSING', 'repressing', human_links)
human_links <- gsub('ENHANCING', 'enhancing', human_links)

paste0(gtex_link_flt$ncRNA, "|", gtex_link_flt$Name_pcGene, "|", gtex_link_flt$Mechanism) %in% human_links