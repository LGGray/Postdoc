library(tidyverse)
library(data.table)
library(VennDiagram)
library(eulerr)


# Function to extract the link key - base|target|mechanism
link_key <- function(df) paste(df$gene_base, df$gene_target, df$mechanism, sep="|")

# Jaccard index function
jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  if (union == 0) {
    return(NA)  # Avoid division by zero
  }
  return(intersection / union)
}

#### Match RefSeq transcript ID to gene name
gtf <- read.delim('GRCm39/GCF_000001635.27_GRCm39_genomic.gtf', comment.char = "#", header = FALSE)
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

########################################
# Read in data from adult aged bodymap #
########################################
bodymap_dirs <- list.dirs('adult_aged_bodymap', recursive = FALSE, full.names = FALSE)
bodymap_dirs <- grep('_', bodymap_dirs, value = TRUE)

bodymap <- lapply(bodymap_dirs, function(dir) {
    tmp <- read.delim(list.files(paste0('adult_aged_bodymap/', dir), pattern = 'links_full_table.txt', full.names = TRUE, recursive = TRUE))
    return(tmp)
})
names(bodymap) <- bodymap_dirs

# Loop over bodymap to annotate genes
bodymap <- lapply(bodymap, function(x) {
  dt_x <- data.table(x)
  dt_x[, gene_base := dt_gtf[.SD, on = .(transcript_id = name_base), gene_ids]]
  dt_x[, biotype_base := dt_gtf[.SD, on = .(transcript_id = name_base), biotype_ids]]
  dt_x[, gene_target := dt_gtf[.SD, on = .(transcript_id = name_target), gene_ids]]
  dt_x[, biotype_target := dt_gtf[.SD, on = .(transcript_id = name_target), biotype_ids]]
#   # Filter out same-gene links (different isoforms of same gene)
#   dt_x <- dt_x[gene_base != gene_target | is.na(gene_base) | is.na(gene_target)]
  as.data.frame(dt_x)
})

save(bodymap, file = 'adult_aged_bodymap/bodymap_links.RData')
load('adult_aged_bodymap/bodymap_links.RData')

#######################################
# Compare link between adult and aged #
#######################################
# Example: compare links in heart tissue between adult and aged
age.colours <- c("adult" = "#94A2AB", "aged" = "#284456")

tissue <- 'He'
adult <- bodymap[[paste0(tissue, '_9w')]]
aged <- bodymap[[paste0(tissue, '_78w')]]

jaccard_index(
  unique(link_key(adult)),
  unique(link_key(aged))
)

adult_sets <- unique(link_key(adult))
aged_sets  <- unique(link_key(aged))

pdf('adult_aged_bodymap/venn_adult_aged_He_links.pdf')
fit <- euler(c(
    'Adult' = length(adult_sets),
    'Aged'  = length(aged_sets),
    "Adult&Aged" = length(intersect(adult_sets, aged_sets))
  ))
plot(fit, quantities = TRUE, fill = age.colours, main='Heart')
dev.off()

# Compatre the number of links gained/lost between adult and aged across all tissues
# Links which change during aging

tissues <- c('Ao', 'Br', 'He', 'Ki', 'Li', 'Lu', 'Mu', 'Sp')
organs <- list(
  'Ao' = 'Aorta',
  'Br' = 'Brain',
  'He' = 'Heart',
  'Ki' = 'Kidney',
  'Li' = 'Liver',
  'Lu' = 'Lung',
  'Mu' = 'Muscle',
  'Sp' = 'Spleen'
)

summ_change <- lapply(tissues, function(tissue) {
    adult <- bodymap[[paste0(tissue, '_9w')]]
    aged <- bodymap[[paste0(tissue, '_78w')]]
    A <- unique(adult_sets[paste0(tissue, '_9w')])
    G <- unique(aged_sets[paste0(tissue, '_78w')])
  tibble(
    organ = organs[[tissue]],
    status   = c("Retained","Gained_in_Aged","Lost_in_Aged"),
    count    = c(length(intersect(A,G)),
                 length(setdiff(G,A)),
                 length(setdiff(A,G)))
  )
}) |> bind_rows() |>
  group_by(organ) |>
  mutate(frac = count / sum(count))

ggplot(summ_change, aes(x=organ, y=frac, fill=status)) +
  geom_col() +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Organ", y="% of links", title="") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("Retained"="#94A2AB", "Gained_in_Aged"="#284456", "Lost_in_Aged"="#D9D9D6"))
ggsave("adult_aged_bodymap/retained_gained_lost_stackedbars.pdf", width=10, height=6)


##################################
# Read in data from TAC and Sham #
##################################
TAC_dirs <- list.dirs('F1_TAC_Sarah', recursive = FALSE, full.names = FALSE)
TAC_dirs <- grep('He', TAC_dirs, value = TRUE)

TAC_SHAM <- lapply(TAC_dirs, function(dir) {
    tmp <- read.delim(list.files(paste0('F1_TAC_Sarah/', dir), pattern = 'links_full_table.txt', full.names = TRUE, recursive = TRUE))
    return(tmp)
})
names(TAC_SHAM) <- TAC_dirs

# Loop over bodymap to annotate genes
TAC_SHAM <- lapply(TAC_SHAM, function(x) {
  dt_x <- data.table(x)
  dt_x[, gene_base := dt_gtf[.SD, on = .(transcript_id = name_base), gene_ids]]
  dt_x[, biotype_base := dt_gtf[.SD, on = .(transcript_id = name_base), biotype_ids]]
  dt_x[, gene_target := dt_gtf[.SD, on = .(transcript_id = name_target), gene_ids]]
  dt_x[, biotype_target := dt_gtf[.SD, on = .(transcript_id = name_target), biotype_ids]]
#   # Filter out same-gene links (different isoforms of same gene)
#   dt_x <- dt_x[gene_base != gene_target | is.na(gene_base) | is.na(gene_target)]
  as.data.frame(dt_x)
})

save(TAC_SHAM, file = 'F1_TAC_Sarah/TAC_SHAM_links.RData')
load('F1_TAC_Sarah/TAC_SHAM_links.RData')

##################################################
# Compare link between adult, aged and TAC heart #
##################################################

adult_heart <- unique(link_key(bodymap[['He_9w']]))
aged_heart  <- unique(link_key(bodymap[['He_78w']]))
sham_heart  <- unique(link_key(TAC_SHAM[['He_TACSham28d_RNA_XistBxC_-+_XX']]))
tac_heart   <- unique(link_key(TAC_SHAM[['He_TAC28d_RNA_XistBxC_-+_XX']]))

sample_colours <- c("Adult_Heart"="#94A2AB", "Aged_Heart"="#284456", "TAC_Heart"="#1F78B4")

pdf('LINKS_study/figures/venn_adult_aged_TAC_links.pdf')
fit <- euler(c(
    'Adult' = length(adult_heart),
    'Aged'  = length(aged_heart),
    'TAC'   = length(tac_heart),
    "Adult&Aged" = length(intersect(adult_heart, aged_heart)),
    "Adult&TAC"  = length(intersect(adult_heart, tac_heart)),
    "Aged&TAC"   = length(intersect(aged_heart, tac_heart)),
    "Adult&Aged&TAC" = length(Reduce(intersect, list(adult_heart, aged_heart, tac_heart)))
  ))
plot(fit, quantities = TRUE, fill = sample_colours, main='Overlap of ncRNA-pcGene Links in Heart')
dev.off()

# Jaccard indices between all conditions
conditions <- list(
  'Adult_Heart' = adult_heart,
  'Aged_Heart'  = aged_heart,
  'TAC_Heart'   = tac_heart
)
jaccard_results <- expand.grid(Condition1 = names(conditions), Condition2 = names(conditions)) %>%
  rowwise() %>%
  mutate(Jaccard_Index = jaccard_index(conditions[[Condition1]], conditions[[Condition2]]))

# Fill as matrix
jaccard_matrix <- matrix(NA, nrow=length(conditions), ncol=length(conditions))
rownames(jaccard_matrix) <- names(conditions)
colnames(jaccard_matrix) <- names(conditions)
for (i in 1:nrow(jaccard_results)) {
  jaccard_matrix[jaccard_results$Condition1[i], jaccard_results$Condition2[i]] <- jaccard_results$Jaccard_Index[i]
}

# Plot heatmap of Jaccard indices - upper triangle
library(ComplexHeatmap)
library(circlize)

pdf('LINKS_study/figures/jaccard_heatmap_adult_aged_TAC_links.pdf')
Heatmap(jaccard_matrix,
        name = "Jaccard Index",
        col = colorRamp2(c(0, 0.5, 1), c("white", "blue", "red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (i < j) {
            grid.text(sprintf("%.2f", jaccard_matrix[i, j]), x, y)
          }
        })
dev.off()

### Repeat Venn diagram but split by repressive and enhancing links
adult_enhancing <- unique(link_key(subset(bodymap[['He_9w']], mechanism == "enhancing")))
adult_repressive <- unique(link_key(subset(bodymap[['He_9w']], mechanism == "repressing")))
aged_enhancing  <- unique(link_key(subset(bodymap[['He_78w']], mechanism == "enhancing")))
aged_repressive  <- unique(link_key(subset(bodymap[['He_78w']], mechanism == "repressing")))
tac_enhancing   <- unique(link_key(subset(TAC_SHAM[['He_TAC28d_RNA_XistBxC_-+_XX']], mechanism == "enhancing")))
tac_repressive   <- unique(link_key(subset(TAC_SHAM[['He_TAC28d_RNA_XistBxC_-+_XX']], mechanism == "repressing")))

pdf('LINKS_study/figures/venn_adult_aged_TAC_links_enhancing.pdf')
fit_enhancing <- euler(c(
    'Adult' = length(adult_enhancing),
    'Aged'  = length(aged_enhancing),
    'TAC'   = length(tac_enhancing),
    "Adult&Aged" = length(intersect(adult_enhancing, aged_enhancing)),
    "Adult&TAC"  = length(intersect(adult_enhancing, tac_enhancing)),
    "Aged&TAC"   = length(intersect(aged_enhancing, tac_enhancing)),
    "Adult&Aged&TAC" = length(Reduce(intersect, list(adult_enhancing, aged_enhancing, tac_enhancing)))
))
plot(fit_enhancing, quantities = TRUE, fill = sample_colours, main='Enhancing Links in Heart')
dev.off()

pdf('LINKS_study/figures/venn_adult_aged_TAC_links_repressive.pdf')
fit_repressive <- euler(c(
    'Adult' = length(adult_repressive),
    'Aged'  = length(aged_repressive),
    'TAC'   = length(tac_repressive),
    "Adult&Aged" = length(intersect(adult_repressive, aged_repressive)),
    "Adult&TAC"  = length(intersect(adult_repressive, tac_repressive)),
    "Aged&TAC"   = length(intersect(aged_repressive, tac_repressive)),
    "Adult&Aged&TAC" = length(Reduce(intersect, list(adult_repressive, aged_repressive, tac_repressive)))
))
plot(fit_repressive, quantities = TRUE, fill = sample_colours, main='Repressive Links in Heart')
dev.off()

# Repeat venn diagram but split by repressive and enhancing links and the base gene biotype is ncRNA
adult_enhancing_ncRNA <- unique(link_key(subset(bodymap[['He_9w']], mechanism == "enhancing" & biotype_base != "mRNA")))
adult_repressive_ncRNA <- unique(link_key(subset(bodymap[['He_9w']], mechanism == "repressing" & biotype_base != "mRNA")))
aged_enhancing_ncRNA  <- unique(link_key(subset(bodymap[['He_78w']], mechanism == "enhancing" & biotype_base != "mRNA")))
aged_repressive_ncRNA  <- unique(link_key(subset(bodymap[['He_78w']], mechanism == "repressing" & biotype_base != "mRNA")))
tac_enhancing_ncRNA   <- unique(link_key(subset(TAC_SHAM[['He_TAC28d_RNA_XistBxC_-+_XX']], mechanism == "enhancing" & biotype_base != "mRNA" & base_gene == 'Magi2'))
tac_repressive_ncRNA   <- unique(link_key(subset(TAC_SHAM[['He_TAC28d_RNA_XistBxC_-+_XX']], mechanism == "repressing" & biotype_base != "mRNA")))

pdf('LINKS_study/figures/venn_adult_aged_TAC_links_enhancing_ncRNA.pdf')
fit_enhancing <- euler(c(
    'Adult' = length(adult_enhancing_ncRNA),
    'Aged'  = length(aged_enhancing_ncRNA),
    'TAC'   = length(tac_enhancing_ncRNA),
    "Adult&Aged" = length(intersect(adult_enhancing_ncRNA, aged_enhancing_ncRNA)),
    "Adult&TAC"  = length(intersect(adult_enhancing_ncRNA, tac_enhancing_ncRNA)),
    "Aged&TAC"   = length(intersect(aged_enhancing_ncRNA, tac_enhancing_ncRNA)),
    "Adult&Aged&TAC" = length(Reduce(intersect, list(adult_enhancing_ncRNA, aged_enhancing_ncRNA, tac_enhancing_ncRNA)))
))
plot(fit_enhancing, quantities = TRUE, fill = sample_colours, main='Enhancing Links in Heart - ncRNA')
dev.off()

pdf('LINKS_study/figures/venn_adult_aged_TAC_links_repressive_ncRNA.pdf')
fit_repressive <- euler(c(
    'Adult' = length(adult_repressive_ncRNA),
    'Aged'  = length(aged_repressive_ncRNA),
    'TAC'   = length(tac_repressive_ncRNA),
    "Adult&Aged" = length(intersect(adult_repressive_ncRNA, aged_repressive_ncRNA)),
    "Adult&TAC"  = length(intersect(adult_repressive_ncRNA, tac_repressive_ncRNA)),
    "Aged&TAC"   = length(intersect(aged_repressive_ncRNA, tac_repressive_ncRNA)),
    "Adult&Aged&TAC" = length(Reduce(intersect, list(adult_repressive_ncRNA, aged_repressive_ncRNA, tac_repressive_ncRNA)))
))
plot(fit_repressive, quantities = TRUE, fill = sample_colours, main='Repressive Links in Heart - ncRNA')
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

human_links <- toupper(intersect(aged_repressive_ncRNA, tac_repressive_ncRNA))
human_links <- gsub('REPRESSING', 'repressing', human_links)

paste0(gtex_link_flt$ncRNA, "|", gtex_link_flt$Name_pcGene, "|", gtex_link_flt$Mechanism) %in% human_links