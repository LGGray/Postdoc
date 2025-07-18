# Load packages
library(biomaRt)
library(ggplot2)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)

### Read in escape data ###
cardiac_cell_types <- read.delim('cardiac_celltypes.txt')
embryo <- read.delim('embryo_organ.txt')
young <- read.delim('young_organ.txt')
adult <- read.delim('adult_organ.txt')
aged <- read.delim('aged_organ.txt')

# Convert median_allelic_ratio to numeric
cardiac_cell_types$median_allelic_ratio <- as.numeric(gsub(',', '.', cardiac_cell_types$median_allelic_ratio))
embryo$median_allelic_ratio<- as.numeric(gsub(',', '.', embryo$median_allelic_ratio))
young$median_allelic_ratio <- as.numeric(gsub(',', '.', young$median_allelic_ratio))
adult$median_allelic_ratio <- as.numeric(gsub(',', '.', adult$median_allelic_ratio))
aged$median_allelic_ratio <- as.numeric(gsub(',', '.', aged$median_allelic_ratio))

### Read in miRNA gff3 file ###
miRNA_gff <- read.delim("mmu.gff3", comment.char = "#", header = FALSE)

file_list <- list(
    cardiac_cell_types = cardiac_cell_types,
    embryo = embryo,
    young = young,
    adult = adult,
    aged = aged
)

#### Iterate over the files

results_list <- list()

for(file in names(file_list)) {
    ### Filter for allelic_ratio =< 0.9 for XCI escape genes ###
    escape_genes <- subset(file_list[[file]], median_allelic_ratio <= 0.9)

    ### Find overlaps between escape genes and miRNA primary transcripts on chrX ###
    overlap_list <- lapply(1:nrow(escape_genes), function(x){
        tmp <- subset(miRNA_gff, V1 == 'chrX' & 
        V3 == 'miRNA_primary_transcript' & 
        V4 >= escape_genes$start[x] & V5 < escape_genes$end[x])
        if (nrow(tmp) > 0) {
            sapply(strsplit(tmp$V9, ';'), function(y) gsub('Name=', '', y[3]))
        } else {
            character(0)
        }
    })
    names(overlap_list) <- escape_genes$name

    escape_genes$miRNA_overlap <- overlap_list

    if (file == "cardiac_cell_types") {
        miRNA_list <- lapply(split(escape_genes, escape_genes$celltype), function(df) {
            unique(unlist(df$miRNA_overlap))
        })
    } else {
        miRNA_list <- lapply(split(escape_genes, escape_genes$organ), function(df) {
            unique(unlist(df$miRNA_overlap))
        })
    }

    result <- subset(escape_genes, sapply(overlap_list, length) > 0)

    results_list[[file]] <- result

    result$miRNA_overlap <- sapply(result$miRNA_overlap, paste, collapse = ",")
    write.table(result, file = paste0(file, '_host_miRNAs.txt'), sep = "\t", row.names = FALSE, quote = FALSE)
}


### calculate the number of mouse miRNAs per chromosome
sort(table(subset(miRNA_gff, V3 == 'miRNA_primary_transcript')$V1))
# calculate the percentage of miRNAs per chromosome
miRNA_counts <- table(subset(miRNA_gff, V3 == 'miRNA_primary_transcript')$V1)
miRNA_percentage <- prop.table(miRNA_counts) * 100

### calculate the number of human miRNAs per chromosome
hsa_miRNA_gff <- read.delim("hsa.gff3", comment.char = "#", header = FALSE)
sort(table(subset(hsa_miRNA_gff, V3 == 'miRNA_primary_transcript')$V1))
hsa_miRNA_counts <- table(subset(hsa_miRNA_gff, V3 == 'miRNA_primary_transcript')$V1)
hsa_miRNA_percentage <- prop.table(hsa_miRNA_counts) * 100


#######################
### Read in miRDB predictions ###
targets <- read.delim('miRDB_v6.0_prediction_result.txt', header=F)

### Filter for mouse miRNAs ###
mouse_targets <- targets[grep('mmu', targets$V1), ]
### Filter for high confidence targets ###
mouse_targets <- mouse_targets[mouse_targets$V3 >= 80, ]

### Read in escape data ###
cardiac_cell_types <- read.delim('cardiac_celltypes.txt')
embryo <- read.delim('embryo_organ.txt')
young <- read.delim('young_organ.txt')
adult <- read.delim('adult_organ.txt')
aged <- read.delim('aged_organ.txt')

### Filter for allelic_ratio < 0.9 for XCI escape genes ###
escape_genes <- cardiac_cell_types[cardiac_cell_types$allelic_ratio <= 0.9, ]

### Get RefSeq mRNA IDs for escape genes ###
# Use Ensembl for mouse
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get RefSeq mRNA IDs (NM_) from gene symbol
results <- getBM(
    attributes = c("external_gene_name", "refseq_mrna"),
    filters = "external_gene_name",
    values = 'Ftx',
    mart = mart
)

refseq <- subset(results, refseq_mrna != "")

### Subset miRDB targets for escape genes ###
escape_targets <- mouse_targets[mouse_targets$V2 %in% refseq$refseq_mrna, ]

# Match refseq$external_gene_name to escape_targets$V2
escape_targets$gene_name <- refseq$external_gene_name[match(escape_targets$V2, refseq$refseq_mrna)]

# Create list of results 
miRNA_gene_list <- split(escape_targets$gene_name, escape_targets$V1)
target_size <- unlist(lapply(miRNA_gene_list, length))

### Plot a histogram of the number of targets per miRNA ###
pdf('figures/celltype_all_miRNA_target_histogram.pdf',6)
ggplot(data = data.frame(target_size), aes(x = target_size)) +
    geom_histogram(binwidth = 1) +
    labs(title = "Number of Targets per miRNA", x = "Number of Targets", y = "Frequency") +
    theme_minimal()
dev.off()

### Read in miRNA gff3 file ###
miRNA_gff <- read.delim("mmu.gff3", comment.char = "#", header = FALSE)
# Subset for X chromosome miRNAs
miRNA_gff_chrX <- miRNA_gff[miRNA_gff$V1 == "chrX", ]
chrX_miRNAs <- lapply(strsplit(miRNA_gff_chrX$V9, ';'), function(x) gsub('Name=', '', x[3]))

# subset miRNA_gene_list for chrX miRNAs
chrX_miRNA_gene_list <- miRNA_gene_list[names(miRNA_gene_list) %in% chrX_miRNAs]
chrX_target_size <- unlist(lapply(chrX_miRNA_gene_list, length))

### Plot a histogram of the number of targets per chrX miRNA ###
pdf('figures/celltype_chrX_miRNA_target_histogram.pdf')
ggplot(data = data.frame(chrX_target_size), aes(x = chrX_target_size)) +
    geom_histogram(binwidth = 1) +
    labs(title = "Number of Targets per chrX miRNA", x = "Number of Targets", y = "Frequency") +
    theme_minimal()
dev.off()   




### Check if the escape gene is a target of host miRNA ###
result <- lapply(1:nrow(escape_genes), function(x) {
    if(escape_genes[x,'miRNA_overlap'] == '') {
        return(character(0))
    }
    escape_genes[x, 'name'] %in% miRNA_gene_list[[escape_genes[x,'miRNA_overlap']]]
})
unlist(result)


result <- sapply(1:nrow(escape_genes), function(x) {
    overlaps <- unlist(strsplit(escape_genes[x, 'miRNA_overlap'], ","))
    any(sapply(overlaps, function(miRNA) {
        escape_genes[x, 'name'] %in% miRNA_gene_list[[miRNA]]
    }))
})### Read in miRDB predictions ###
targets <- read.delim('miRDB_v6.0_prediction_result.txt', header=F)

### Filter for mouse miRNAs ###
mouse_targets <- targets[grep('mmu', targets$V1), ]
### Filter for high confidence targets ###
mouse_targets <- mouse_targets[mouse_targets$V3 >= 80, ]

### Read in escape data ###
cardiac_cell_types <- read.delim('cardiac_celltypes.txt')
embryo <- read.delim('embryo_organ.txt')
young <- read.delim('young_organ.txt')
adult <- read.delim('adult_organ.txt')
aged <- read.delim('aged_organ.txt')

### Filter for allelic_ratio < 0.9 for XCI escape genes ###
escape_genes <- cardiac_cell_types[cardiac_cell_types$allelic_ratio <= 0.9, ]

### Get RefSeq mRNA IDs for escape genes ###
# Use Ensembl for mouse
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get RefSeq mRNA IDs (NM_) from gene symbol
results <- getBM(
    attributes = c("external_gene_name", "refseq_mrna"),
    filters = "external_gene_name",
    values = 'Ftx',
    mart = mart
)

refseq <- subset(results, refseq_mrna != "")

### Subset miRDB targets for escape genes ###
escape_targets <- mouse_targets[mouse_targets$V2 %in% refseq$refseq_mrna, ]

# Match refseq$external_gene_name to escape_targets$V2
escape_targets$gene_name <- refseq$external_gene_name[match(escape_targets$V2, refseq$refseq_mrna)]

# Create list of results 
miRNA_gene_list <- split(escape_targets$gene_name, escape_targets$V1)
target_size <- unlist(lapply(miRNA_gene_list, length))

### Plot a histogram of the number of targets per miRNA ###
pdf('figures/celltype_all_miRNA_target_histogram.pdf',6)
ggplot(data = data.frame(target_size), aes(x = target_size)) +
    geom_histogram(binwidth = 1) +
    labs(title = "Number of Targets per miRNA", x = "Number of Targets", y = "Frequency") +
    theme_minimal()
dev.off()

### Read in miRNA gff3 file ###
miRNA_gff <- read.delim("mmu.gff3", comment.char = "#", header = FALSE)
# Subset for X chromosome miRNAs
miRNA_gff_chrX <- miRNA_gff[miRNA_gff$V1 == "chrX", ]
chrX_miRNAs <- lapply(strsplit(miRNA_gff_chrX$V9, ';'), function(x) gsub('Name=', '', x[3]))

# subset miRNA_gene_list for chrX miRNAs
chrX_miRNA_gene_list <- miRNA_gene_list[names(miRNA_gene_list) %in% chrX_miRNAs]
chrX_target_size <- unlist(lapply(chrX_miRNA_gene_list, length))

### Plot a histogram of the number of targets per chrX miRNA ###
pdf('figures/celltype_chrX_miRNA_target_histogram.pdf')
ggplot(data = data.frame(chrX_target_size), aes(x = chrX_target_size)) +
    geom_histogram(binwidth = 1) +
    labs(title = "Number of Targets per chrX miRNA", x = "Number of Targets", y = "Frequency") +
    theme_minimal()
dev.off()   




### Check if the escape gene is a target of host miRNA ###
result <- lapply(1:nrow(escape_genes), function(x) {
    if(escape_genes[x,'miRNA_overlap'] == '') {
        return(character(0))
    }
    escape_genes[x, 'name'] %in% miRNA_gene_list[[escape_genes[x,'miRNA_overlap']]]
})
unlist(result)


result <- sapply(1:nrow(escape_genes), function(x) {
    overlaps <- unlist(strsplit(escape_genes[x, 'miRNA_overlap'], ","))
    any(sapply(overlaps, function(miRNA) {
        escape_genes[x, 'name'] %in% miRNA_gene_list[[miRNA]]
    }))
})