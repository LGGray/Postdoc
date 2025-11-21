#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript filter_bedpe.R <input_bedpe> <input_links_table> <output_prefix>\n")
  cat("Example: Rscript filter_bedpe.R AL_65.bedpe AL_65_links_full_table.txt AL_65.filtered\n")
  quit(status = 1)
}

infile_bedpe <- args[1]
infile_links <- args[2]
outfile_prefix <- args[3]

# Load coding/noncoding classification
noncoding_coding <- read.delim('/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/LINKS_study/Refseq_coding_noncoding.txt', header = TRUE)
coding <- subset(noncoding_coding, coding_noncoding == "coding")$name
noncoding <- subset(noncoding_coding, coding_noncoding == "noncoding")$name

# Read and filter bedpe
cat("Reading bedpe file:", infile_bedpe, "\n")
bedpe <- read.delim(infile_bedpe, skip=1, header=TRUE)
bedpe$base <- unlist(lapply(bedpe$name, function(x) strsplit(x, '_')[[1]][1]))
bedpe$target <- unlist(lapply(bedpe$name, function(x) strsplit(x, '_')[[1]][2]))
bedpe_filtered <- subset(bedpe, base %in% noncoding & target %in% coding)

# Write filtered bedpe
outfile_bedpe <- paste0(outfile_prefix, ".bedpe")
cat("Writing filtered bedpe to:", outfile_bedpe, "\n")
first_line <- readLines(infile_bedpe, n = 1)
writeLines(first_line, con = outfile_bedpe)
write.table(bedpe_filtered, file = outfile_bedpe, sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = TRUE, append = TRUE)

# Read and filter links table
cat("Reading links table:", infile_links, "\n")
links <- read.delim(infile_links, header = TRUE)
links_filtered <- subset(links, name_base %in% noncoding & name_target %in% coding)

# Write filtered links table
outfile_links <- paste0(outfile_prefix, "_links_full_table.txt")
cat("Writing filtered links table to:", outfile_links, "\n")
write.table(links_filtered, file = outfile_links, sep = '\t', quote = FALSE,
            row.names = FALSE, col.names = TRUE)

cat("Done! Filtered", nrow(bedpe_filtered), "links from", nrow(bedpe), "total.\n")