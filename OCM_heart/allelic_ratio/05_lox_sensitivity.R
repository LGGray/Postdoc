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
source("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/OCM_heart/allelic_ratio/00_functions.R")

metadata_whole_chr <- read.table(
  file.path(CUTOFF_DIR, 'whole_chr_cell_metadata.txt'),
  sep = '\t', header = TRUE, row.names = 1,
  stringsAsFactors = FALSE, check.names = FALSE)

metadata_ceb <- read.table(
  file.path(CEB_DIR, 'core_escape_block_cell_metadata.txt'),
  sep = '\t', header = TRUE, row.names = 1,
  stringsAsFactors = FALSE, check.names = FALSE)

# TWO RATIO CONVENTIONS MEET HERE, and they used to meet unlabelled.
#
#   whole_chr_cell_metadata.txt (from 02, via 10_build_ratio_table.R) carries
#   `allelic_ratio` = ar_dom = max(A1, A2) / total, bounded 0.5-1.
#
#   core_escape_block_cell_metadata.txt (from 04) carries `ar_a1`, the
#   DIRECTIONAL Allelome.PRO2 ratio A1 / total, bounded 0-1, plus `ar_dom`.
#   It used to call the directional column `allelic_ratio` too, which is what
#   made this join look safe when it was not.
#
# Only `monoallelic` crosses the join, and 04 now builds it from ar_dom, so a
# cell that lost EITHER allele is excluded. Under the old directional test only
# the cells that lost A2 were - half the LOX cells this analysis exists to
# remove were being kept.
stopifnot("monoallelic" %in% names(metadata_ceb))
if ("allelic_ratio" %in% names(metadata_whole_chr) &&
    min(metadata_whole_chr$allelic_ratio, na.rm = TRUE) < 0.5 - 1e-9) {
  warning("whole_chr_cell_metadata.txt has allelic_ratio below 0.5, so it is ",
          "the DIRECTIONAL ratio, not ar_dom as 02 and 06 assume. ",
          "Rebuild it with 10_build_ratio_table.R.")
}

# Repeat the dispersion/heterogeneity test after excluding putative LOX cells
# (ar_dom >= 0.9) - a lost inactive allele collapses a cell's ratio to a fixed
# point with no sampling variance, which can trivially shift/mask the
# between-condition dispersion comparison for an escape gene.
lox_barcodes_ceb <- rownames(subset(metadata_ceb, monoallelic == 1))
message(sprintf("LOX-like cells excluded: %d of %d in the core escape block",
                length(lox_barcodes_ceb), nrow(metadata_ceb)))
n_before <- nrow(metadata_whole_chr)

metadata_whole_chr <- metadata_whole_chr[!rownames(metadata_whole_chr) %in% lox_barcodes_ceb, ]
message(sprintf("whole-chrX cells: %d -> %d after the exclusion",
                n_before, nrow(metadata_whole_chr)))
if (nrow(metadata_whole_chr) == n_before) {
  warning("The exclusion removed nothing. The two metadata files are keyed on ",
          "cell barcodes that do not overlap - check both were built from the ",
          "same tree and the same Seurat object.")
}

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


