noncoding_coding <- read.delim('LINKS_study/Refseq_coding_noncoding.txt', header = TRUE)
coding <- subset(noncoding_coding, coding_noncoding == "coding")$name
noncoding <- subset(noncoding_coding, coding_noncoding == "noncoding")$name

# Adult heart
He_9w <- read.delim('adult_aged_bodymap/He_9w/AL_65/AL_65.bedpe', skip=1, header=TRUE)
He_9w$base <- unlist(lapply(He_9w$name, function(x) strsplit(x, '_')[[1]][1]))
He_9w$target <- unlist(lapply(He_9w$name, function(x) strsplit(x, '_')[[1]][2]))

He_9w <- subset(He_9w, base %in% noncoding & target %in% coding)

infile <- 'adult_aged_bodymap/He_9w/AL_65/AL_65.bedpe'
outfile <- 'adult_aged_bodymap/He_9w/AL_65/AL_65.filtered.bedpe'
first_line <- readLines(infile, n = 1)
writeLines(first_line, con = outfile)
write.table(He_9w, file = outfile, sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = TRUE, append = TRUE)

# Aged heart
He_78w <- read.delim('adult_aged_bodymap/He_78w/AL_65/AL_65.bedpe', skip=1, header=TRUE)
He_78w$base <- unlist(lapply(He_78w$name, function(x) strsplit(x, '_')[[1]][1]))
He_78w$target <- unlist(lapply(He_78w$name, function(x) strsplit(x, '_')[[1]][2]))

He_78w <- subset(He_78w, base %in% noncoding & target %in% coding)

infile <- 'adult_aged_bodymap/He_78w/AL_65/AL_65.bedpe'
outfile <- 'adult_aged_bodymap/He_78w/AL_65/AL_65.filtered.bedpe'
first_line <- readLines(infile, n = 1)
writeLines(first_line, con = outfile)
write.table(He_78w, file = outfile, sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = TRUE, append = TRUE)

# TAC
TAC <- read.delim('F1_TAC_Sarah/He_TAC28d_RNA_XistBxC_-+_XX/AL_65/AL_65.bedpe', skip=1, header=TRUE)
TAC$base <- unlist(lapply(TAC$name, function(x) strsplit(x, '_')[[1]][1]))
TAC$target <- unlist(lapply(TAC$name, function(x) strsplit(x, '_')[[1]][2]))
TAC <- subset(TAC, base %in% noncoding & target %in% coding)

infile <- 'F1_TAC_Sarah/He_TAC28d_RNA_XistBxC_-+_XX/AL_65/AL_65.bedpe'
outfile <- 'F1_TAC_Sarah/He_TAC28d_RNA_XistBxC_-+_XX/AL_65/AL_65.filtered.bedpe'
first_line <- readLines(infile, n = 1)
writeLines(first_line, con = outfile)
write.table(TAC, file = outfile, sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = TRUE, append = TRUE)

# Sham
Sham <- read.delim('F1_TAC_Sarah/He_TACSham28d_RNA_XistBxC_-+_XX/AL_65/AL_65.bedpe', skip=1, header=TRUE)
Sham$base <- unlist(lapply(Sham$name, function(x) strsplit(x, '_')[[1]][1]))
Sham$target <- unlist(lapply(Sham$name, function(x) strsplit(x, '_')[[1]][2]))
Sham <- subset(Sham, base %in% noncoding & target %in% coding)
infile <- 'F1_TAC_Sarah/He_TACSham28d_RNA_XistBxC_-+_XX/AL_65/AL_65.bedpe'
outfile <- 'F1_TAC_Sarah/He_TACSham28d_RNA_XistBxC_-+_XX/AL_65/AL_65.filtered.bedpe'
first_line <- readLines(infile, n = 1)
writeLines(first_line, con = outfile)
write.table(Sham, file = outfile, sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = TRUE, append = TRUE)