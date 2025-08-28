library(ggplot2)
library(edgeR)

# Read in human and mouse miRNA coordinated and plot a box plot for the frequency of miRNAs per chromosome
human_miRNA <- read.delim('hsa.gff3', sep='\t', comment.char='#', header=FALSE)
human_miRNA <- subset(human_miRNA, V3 == 'miRNA')
human_miRNA$species <- 'Human'
mouse_miRNA <- read.delim('mmu.gff3', sep='\t', comment.char='#', header=FALSE)
mouse_miRNA <- subset(mouse_miRNA, V3 == 'miRNA')
mouse_miRNA$species <- 'Mouse'

combined_miRNA <- rbind(human_miRNA, mouse_miRNA)

# Histogram of miRNA counts per chromosome, ordered by frequency and split by species to human and mouse are separate
pdf('miRNA_histogram.pdf')
ggplot(combined_miRNA, aes(x=V1, fill=species)) +
  geom_bar() +
  labs(title='miRNA Counts per Chromosome', x='Chromosome', y='Count') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  facet_wrap(~species, scales='free_x')
dev.off()

library(ggplot2)
library(dplyr)
library(tidytext)

# Read, keep only miRNA features, tag species
human_miRNA <- read.delim("hsa.gff3", sep = "\t", comment.char = "#", header = FALSE) |>
  subset(V3 == "miRNA") |>
  transform(species = "Human")

mouse_miRNA <- read.delim("mmu.gff3", sep = "\t", comment.char = "#", header = FALSE) |>
  subset(V3 == "miRNA") |>
  transform(species = "Mouse")

combined_miRNA <- rbind(human_miRNA, mouse_miRNA)

# Count miRNAs per chromosome within each species
counts <- combined_miRNA |>
  count(species, V1, name = "count") |>
  mutate(chr_ordered = reorder_within(V1, -count, species),
         highlight = ifelse(V1 %in% c("chrX", "X"), "chrX", "Other"))

pdf("miRNA_histogram.pdf", width = 10, height = 6)
ggplot(counts, aes(x = chr_ordered, y = count, fill = highlight)) +
  geom_col() +
  facet_wrap(~ species, scales = "free_x") +
  scale_x_reordered() +
  scale_fill_manual(values = c("chrX" = "red", "Other" = "grey70")) +
  labs(title = "miRNA Counts per Chromosome",
       x = "Chromosome",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.off()



file.list <- list.files('.', pattern='.csv')
edgeR_list <- lapply(file.list, function(x) {
  read.csv(x, header=TRUE)
})
names(edgeR_list) <- gsub('.edgeR.csv', '', file.list)

degs <- lapply(edgeR_list, function(x){
    subset(x, FDR < 0.05)
})
degs[['Cardiac_fibroblasts']]$gene