# ---------------------------------------------------------------------------
# 10 - Build the per-cell, per-chromosome allelic ratio table from an
#      Allelome.PRO2 output tree.
#
# This is the ingestion step that was commented out at the top of
# 02_whole_chrX.R. It exists as a script because there is now more than one
# tree to ingest: the read-level one and the deduplicated one written by
# scAllelome_dedup.slurm, and because the autosomes are now kept rather than
# discarded at read time.
#
# Writes ALLELIC_RATIOS_FILE (whole_chr_allelic_ratios.txt) under RESULTS_ROOT,
# in the shape 02 and 06 already expect, plus an `autosomal` pseudo-chromosome
# row per cell pooling chr1-19.
#
# Run from the OCM_heart/ directory, in seurat_env:
#
#   # read-level tree, the original
#   Rscript allelic_ratio/10_build_ratio_table.R
#
#   # deduplicated tree
#   TREE=Allelome.PRO2_dedup RESULTS_ROOT=Allelic_ratio_results_dedup \
#     Rscript allelic_ratio/10_build_ratio_table.R
#
# Then everything downstream follows the same RESULTS_ROOT:
#   RESULTS_ROOT=Allelic_ratio_results_dedup Rscript allelic_ratio/02_whole_chrX.R
# ---------------------------------------------------------------------------
source("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/OCM_heart/allelic_ratio/00_functions.R")

TREE    <- Sys.getenv("TREE", "Allelome.PRO2")
SAMPLES <- strsplit(Sys.getenv("SAMPLES", "9w,78w,Sham,TAC"), "[, ]+")[[1]]

message("Tree:         ", TREE)
message("RESULTS_ROOT: ", RESULTS_ROOT)
message("Samples:      ", paste(SAMPLES, collapse = ", "))

df <- load_allelome_tree(TREE, SAMPLES)

# Which samples actually have autosomal rows. scAllelome.slurm gave 9w a
# chrX-only annotation, so in the read-level tree 9w has none and every
# autosomal comparison involving it is empty rather than wrong. Say so loudly
# here rather than letting it surface as a silently missing facet later.
cov <- df %>%
  group_by(sample) %>%
  summarise(n_cells   = n_distinct(cell_barcode),
            has_chrX  = any(chr == "chrX"),
            n_autosom = n_distinct(chr[chr %in% AUTOSOMES]),
            .groups = "drop")
print(as.data.frame(cov))
if (any(cov$n_autosom == 0)) {
  warning("No autosomal rows for: ",
          paste(cov$sample[cov$n_autosom == 0], collapse = ", "),
          " -- that tree was built with a chrX-only annotation, so it carries ",
          "no technical control. Rebuild it genome-wide before using 11.")
}

auto <- collapse_autosomes(df)
out  <- bind_rows(df, auto)

# 02 and 06 read a column literally named allelic_ratio and treat >= 0.9 as
# monoallelic, i.e. the dominant-allele convention. Write that as the default
# and keep the directional ratio alongside it, so nothing downstream changes
# meaning but the direction is still recoverable.
out$allelic_ratio <- out$ar_dom

write.table(out, ALLELIC_RATIOS_FILE, sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote ", ALLELIC_RATIOS_FILE, "  (", nrow(out), " rows, ",
        length(unique(out$cell_barcode)), " cells)")
