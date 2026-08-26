# ---------------------------------------------------------------------------
# 05 - LOX sensitivity analysis.
#
# Split out of the core escape section deliberately: this does not analyse the
# core escape block itself, it uses the LOX calls derived there to RE-RUN the
# whole-chrX dispersion tests from 02. Keeping it separate makes it obvious
# that these results supersede/qualify 02's, rather than being new findings.
#
# Run from the OCM_heart/ directory:  Rscript allelic_ratio/05_lox_sensitivity.R
# Requires: 02_whole_chrX.R and 04_core_escape.R to have been run.
# ---------------------------------------------------------------------------
source("allelic_ratio/00_functions.R")

metadata_whole_chr <- read.table(
  file.path(CUTOFF_DIR, 'whole_chr_cell_metadata.txt'),
  sep = '\t', header = TRUE, row.names = 1,
  stringsAsFactors = FALSE, check.names = FALSE)

metadata_ceb <- read.table(
  file.path(CEB_DIR, 'core_escape_block_cell_metadata.txt'),
  sep = '\t', header = TRUE, row.names = 1,
  stringsAsFactors = FALSE, check.names = FALSE)

# Repeat the dispersion/heterogeneity test after excluding putative LOX cells
# (AR >= 0.9) - a lost inactive allele collapses a cell's ratio to a fixed
# point with no sampling variance, which can trivially shift/mask the
# between-condition dispersion comparison for an escape gene.
lox_barcodes_ceb <- rownames(subset(metadata_ceb, monoallelic == 1))

metadata_whole_chr <- metadata_whole_chr[!rownames(metadata_whole_chr) %in% lox_barcodes_ceb, ]

adult_vs_aged_lrt_noLOX <- split(metadata_whole_chr, metadata_whole_chr$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("9w", "78w"))
    if (length(unique(df$sample)) < 2) return(NULL)
    run_dispersion_lrt(df, ref_level = "9w", comparison_level = "78w")
  })
adult_vs_aged_lrt_noLOX <- bind_rows(adult_vs_aged_lrt_noLOX, .id = "celltype")
adult_vs_aged_lrt_noLOX$FDR <- p.adjust(adult_vs_aged_lrt_noLOX$Pr..Chisq., method = "fdr")

subset(adult_vs_aged_lrt_noLOX, FDR < 0.05)

Sham_vs_TAC_lrt_noLOX <- split(metadata_whole_chr, metadata_whole_chr$celltype) %>%
  lapply(function(x) {
    df <- subset(x, sample %in% c("Sham", "TAC"))
    if (length(unique(df$sample)) < 2) return(NULL)
    run_dispersion_lrt(df, ref_level = "Sham", comparison_level = "TAC")
  })
Sham_vs_TAC_lrt_noLOX <- bind_rows(Sham_vs_TAC_lrt_noLOX, .id = "celltype")
Sham_vs_TAC_lrt_noLOX$FDR <- p.adjust(Sham_vs_TAC_lrt_noLOX$Pr..Chisq., method = "fdr")

subset(Sham_vs_TAC_lrt_noLOX, FDR < 0.05)


