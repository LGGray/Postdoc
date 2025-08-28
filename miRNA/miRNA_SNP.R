### Code to overlap CAST or BL6 SNPs with miRNA primary transcripts on chrX ###

# Read in miRNA coordinates
miRNA_gff <- read.delim('miRNA/mmu.gff3', comment.char = "#", header = FALSE)
miRNA_gff <- subset(miRNA_gff, V3 == 'miRNA' & V1 == 'chrX')
miRNA_gff$V10 <- unlist(lapply(miRNA_gff$V9, function(x) {
    # Extract the miRNA name from the GFF file
    gsub('Name=', '', strsplit(x, ';')[[1]][3])
}))

# Read in BL6 SNP
BL6 <- read.delim('mouse_SNP/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz', comment.char = "#", header = FALSE)
BL6_X <- subset(BL6, V1 == 'X')

hits_BL6 <- lapply(BL6_X$V2, function(x){
    # Find the miRNA that overlaps with the SNP position
    subset(miRNA_gff, V4 <= x & V5 >= x)
})
hits_BL6[sapply(hits_BL6, nrow) > 0]


# Read in CAST SNPs
CAST <- read.delim('mouse_SNP/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz', comment.char = "#", header = FALSE)
CAST_X <- subset(CAST, V1 == 'X')

hits_cast <- lapply(CAST_X$V2, function(x){
    # Find the miRNA that overlaps with the SNP position
    subset(miRNA_gff, V4 <= x & V5 >= x)
})
names(hits_cast) <- paste0(CAST_X$V2, ':', CAST_X$V4, '-', CAST_X$V5)

result <- hits_cast[sapply(hits_cast, nrow) > 0]

save(result, file = 'miRNA/CAST_mature_miRNA_SNPs.RData')


foo <- unlist(lapply(result, function(x) {
    # Extract the miRNA name from the GFF file
    gsub('Name=', '', strsplit(x$V9, ';')[[1]][3])
}))

### Read in cardiac miRNAs ###
library(dplyr)
library(tidyr)
cardiac_miRNA <- read.delim('miRNA/cardiac_miRNA/miRNA_counts_matrix.tsv')

cardiac_miRNA$miRNA[grep('miR-374b', cardiac_miRNA$miRNA)]

target_miRNAs <- subset(cardiac_miRNA, miRNA %in% c('mmu-miR-92a-2-5p', 'mmu-miR-743a-5p', 'mmu-miR-741-3p', 'mmu-miR-741-5p',
'mmu-miR-871-5p', 'mmu-miR-7092-5p', 'mmu-miR-223-3p', 'mmu-miR-223-5p', 'mmu-miR-421-3p', 'mmu-miR-421-5p', 
"mmu-miR-374b-3p", "mmu-miR-374b-5p", "mmu-miR-374c-3p", "mmu-miR-374c-5p"))

target_miRNAs_scaled <- data.frame(t(scale(t(target_miRNAs))))

# Long format
target_miRNAs_long <- target_miRNAs %>%
  pivot_longer(
    cols = -miRNA, 
    names_to = "cell", 
    values_to = "count"
  ) %>%
  # Add new column cell type where CF = Fibroblasts, CM = Myocytes, EC = Endothelial Cells, CM = Macrophages
  mutate(celltype = sub("\\..*", "", cell)) %>%
  mutate(celltype = case_when(
    celltype == 'CF' ~ 'Fibroblasts',
    celltype == 'CM' ~ 'Myocytes',
    celltype == 'EC' ~ 'Endothelial Cells',
    celltype == 'MC' ~ 'Macrophages',
    TRUE ~ celltype
  ))

pdf('miRNA/cardiac_miRNA_boxplot.pdf')
ggplot(target_miRNAs_long, aes(x = miRNA, y = count, fill = celltype)) +
  geom_boxplot(
    outlier.shape = NA, alpha = 0.6,
    position = position_dodge2(width = dodge_w, preserve = "single")
  ) +
  geom_jitter(
    size = 2, alpha = 0.8,
    position = position_jitterdodge(
      jitter.width = 0.15, jitter.height = 0, dodge.width = dodge_w
    )
  ) +
  labs(y = "Z-score scaled counts", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  legend.title = element_blank())
dev.off()




library(dplyr)
library(tidyr)

# Example df is called cardiac_miRNA
df_avg <- cardiac_miRNA %>%
  pivot_longer(
    cols = -miRNA, 
    names_to = "cell", 
    values_to = "count"
  ) %>%
  mutate(celltype = sub("\\..*", "", cell)) %>%   # take part before the dot
  group_by(miRNA, celltype) %>%
  summarise(avg_count = mean(count), .groups = "drop") %>%
  pivot_wider(
    names_from = celltype, 
    values_from = avg_count
  ) %>%
  mutate(across(-miRNA, round)) %>%
  data.frame()

subset(df_avg, miRNA %in% miRNA_gff$V10)

## read in small RNA atlas ###
atlas <- read.delim('miRNA/small_RNA_atlas.tsv')

rownames(atlas)[grep('Mir374', rownames(atlas))]

target_miRNAs <- atlas[c(2092, 2112, 2118, 2123, 2143, 2152, 2157, 2158, 2159), grep('Heart', colnames(atlas))]
scaled_target_miRNAs <- t(data.frame(scale(t(target_miRNAs))))
scaled_target_miRNAs <- scaled_target_miRNAs[1:8,]

write.table(target_miRNAs, file = 'miRNA/CAST_mature_miRNA_Atlas.txt', sep = '\t', quote = FALSE, row.names = TRUE)

atlas[c(2092, 2112, 2118, 2123, 2143), grep('Lung', colnames(atlas))]

# Turn target_miRNAs into a long format
library(dplyr)
library(tidyr)
library(stringr)

df_long <- data.frame(scaled_target_miRNAs) %>%
  tibble::rownames_to_column("miRNA") %>%
  pivot_longer(
    cols = -miRNA,
    names_to = "sample",
    values_to = "count"
  ) %>%
  mutate(sex = ifelse(str_detect(sample, "Female"), "Female", "Male"))



# keep a consistent sex order (optional)
df_long$sex <- factor(df_long$sex, levels = c("Female","Male"))

dodge_w <- 0.8

pdf('miRNA/miRNA_expression_atlas.pdf')
ggplot(df_long, aes(x = miRNA, y = count, fill = sex)) +
  geom_boxplot(
    outlier.shape = NA, alpha = 0.6,
    position = position_dodge2(width = dodge_w, preserve = "single")
  ) +
  geom_jitter(
    aes(colour = sex),
    size = 2, alpha = 0.8,
    position = position_jitterdodge(
      jitter.width = 0.15, jitter.height = 0, dodge.width = dodge_w
    )
  ) +
  labs(y = "Z-score scaled counts", x = "") +
  theme_minimal() +
  theme(legend.title = element_blank())
dev.off()



# run Wilcoxon rank-sum test for each miRNA
df_long %>%
  group_by(miRNA) %>%
  wilcox.test(count ~ sex, data=.)


lapply(split(df_long, df_long$miRNA), function(x){
  t.test(count ~ sex, data=x)
})