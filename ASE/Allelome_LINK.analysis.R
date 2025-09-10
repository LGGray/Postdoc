library(data.table)
library(ggplot2)
library(UpSetR)

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

dt_gtf <- data.table(
  seqname = gtf_transcripts$V1,
  start = gtf_transcripts$V4,
  end = gtf_transcripts$V5,
  transcript_id = transcript_ids,
  gene_id = gene_ids
)

# make gtf data.tables
dt_gtf <- data.table(gtf_transcripts)[, .(transcript_id = transcript_ids, gene_ids)]

#### Read in Allelome.PRO2 ouput
adult.files <- list.files('9w', pattern = 'locus_table.txt', full.names = TRUE, recursive = TRUE)
adult.Allelome <- lapply(adult.files, function(file) {
  read.delim(file)
})
names(adult.Allelome) <- unlist(lapply(adult.files, function(x){
  tmp <- strsplit(x, '/')[[1]][[3]]
  strsplit(tmp, '_')[[1]][4]
}))

adult.Allelome <- lapply(adult.Allelome, function(x) {
  dt_x <- data.table(x)
  dt_x[, gene := dt_gtf[.SD, on = .(transcript_id = name), gene_ids]]
  as.data.frame(dt_x)
})

save(adult.Allelome, file = "adult.Allelome.PRO2.RData")

aged.files <- list.files('78w', pattern = 'locus_table.txt', full.names = TRUE, recursive = TRUE)
aged.Allelome <- lapply(aged.files, function(file) {
  read.delim(file)
})
names(aged.Allelome) <- unlist(lapply(aged.files, function(x){
  tmp <- strsplit(x, '/')[[1]][[3]]
  strsplit(tmp, '_')[[1]][4]
}))

aged.Allelome <- lapply(aged.Allelome, function(x) {
  dt_x <- data.table(x)
  dt_x[, gene := dt_gtf[.SD, on = .(transcript_id = name), gene_ids]]
  as.data.frame(dt_x)
})

save(aged.Allelome, file = "aged.Allelome.PRO2.RData")


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
  dt_x[, gene_target := dt_gtf[.SD, on = .(transcript_id = name_target), gene_ids]]
  as.data.frame(dt_x)
})

save(adult.LINK, file = "adult.Allelome.LINK.RData")

aged.LINK <- lapply(aged.LINK, function(x) {
  dt_x <- data.table(x)
  dt_x[, gene_base := dt_gtf[.SD, on = .(transcript_id = name_base), gene_ids]]
  dt_x[, gene_target := dt_gtf[.SD, on = .(transcript_id = name_target), gene_ids]]
  as.data.frame(dt_x)
})

save(aged.LINK, file = "aged.Allelome.LINK.RData")

#### Visualising the data ###
library(ggplot2)

# How many links per cell type for adult and aged. Barplot
adult_links_count <- sapply(adult.LINK, nrow)
aged_links_count <- sapply(aged.LINK, nrow)

df_links <- data.frame(
  celltype = rep(names(adult_links_count), 2),
  age = rep(c("adult", "aged"), each = length(adult_links_count)),
  links = c(adult_links_count, aged_links_count)
)

pdf("figures/links_per_celltype_adult_aged.pdf", width = 8, height = 6)
ggplot(df_links, aes(x = celltype, y = links, fill = age)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of Links per cell type",
       x = "Cell Type",
       y = "Number of Links") +
  theme_minimal()
dev.off()

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


# How many links are concordant between adult and aged
concordant <- lapply(names(aged.LINK), function(ct) {
    intersect(unique(paste(adult.LINK[[ct]]$gene_base, adult.LINK[[ct]]$gene_target, adult.LINK[[ct]]$mechanism, sep='_')),
            unique(paste(aged.LINK[[ct]]$gene_base, aged.LINK[[ct]]$gene_target, aged.LINK[[ct]]$mechanism, sep='_')))
})
names(concordant) <- names(aged.LINK)
concordant_matrix <- fromList(concordant)

pdf("figures/concordant_links_upset.pdf", width = 8, height = 6, onefile = FALSE)
upset(concordant_matrix, nsets = length(concordant), order.by = "freq", 
      main.bar.color = "steelblue", sets.bar.color = "darkorange")
dev.off()

# How many links are unique to each cell type - adult and aged
unique_links <- lapply(names(aged.LINK), function(ct) {
    setdiff(unique(paste(adult.LINK[[ct]]$gene_base, adult.LINK[[ct]]$gene_target, adult.LINK[[ct]]$mechanism, sep='_')),
            unique(paste(aged.LINK[[ct]]$gene_base, aged.LINK[[ct]]$gene_target, aged.LINK[[ct]]$mechanism, sep='_')))
})
names(unique_links) <- names(aged.LINK)
unique_links_matrix <- fromList(unique_links)

pdf("figures/discordant_links_upset.pdf", width = 8, height = 6, onefile = FALSE)
upset(unique_links_matrix, nsets = length(unique_links), order.by = "freq", 
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
