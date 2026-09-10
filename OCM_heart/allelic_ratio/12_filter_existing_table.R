# ---------------------------------------------------------------------------
# 12 - Build a doublet-free whole-chrX ratio table from an ALREADY INGESTED
#      table, without walking an Allelome.PRO2 tree.
#
# Why this exists: the read-level Allelome.PRO2 tree that 10_build_ratio_table.R
# ingests is no longer on disk (only Allelome.PRO2_all_genes and
# Allelome.PRO2_core_escape_new remain), so 10 has nothing to read. The
# per-nucleus numbers it produced do survive, in
# Allelic_ratio_results/whole_chr_allelic_ratios.txt, and that is all this step
# needs. The per-cell BAMs under OCM/sinto are intact, so the tree can be
# rebuilt with slurm/scAllelome.slurm if the full ingestion is ever wanted back.
#
# This does what 10 does, minus the tree walk:
#   - drops the scDblFinder doublet nuclei (DOUBLET_FILE, see 00_functions.R)
#   - derives `sample` from the <sample>_<barcode> key
#   - sums to one row per (cell, chr)
#   - computes ar_a1 and ar_dom, and sets allelic_ratio <- ar_dom
#   - adds the `autosomal` pseudo-chromosome row that stage 11 needs
#
# CONVENTION FIX -- read this before comparing against the old results.
#
# The source table's `allelic_ratio` column is ar_a1: directional, A1 / total,
# on 0-1, which is what Allelome.PRO2's locus table gives. 02, 06 and 10 all
# read a column of that name as ar_dom, max(A1, A2) / total on 0.5-1, and cut
# monoallelic at >= MONO_AR. On chrX the two agree for 98.6% of cells because
# B6 (A1) is the active X in this cross -- but they disagree for the ~1.4% with
# ar_a1 < 0.5, the cells that lost the ACTIVE B6 X. Those are exactly the
# cells 08_lox_calling.R is about, and reading ar_a1 as ar_dom put them at the
# biallelic end instead of the monoallelic one.
#
# This step writes ar_dom into allelic_ratio, per the contract documented at
# 10_build_ratio_table.R:64-68. So results under RESULTS_ROOT differ from
# Allelic_ratio_results/ for TWO reasons: the doublets are gone, AND those
# ~1.4% of cells are now read in the right direction. The count is printed
# below so the two effects can be told apart.
#
# Run from the OCM_heart/ directory, in seurat_env:
#
#   DOUBLET_FILE=Allelic_ratio_results/scDblFinder_per_cell.txt \
#     RESULTS_ROOT=Allelic_ratio_results_nodoublet \
#     Rscript allelic_ratio/12_filter_existing_table.R
#
# Then 02, 03, 05, 06 and 11 follow the same RESULTS_ROOT.
# ---------------------------------------------------------------------------
# POSTDOC_ROOT lets these scripts be parsed and syntax-checked off the cluster;
# unset, it is the cluster path these have always used, so job scripts need no change.
source(file.path(Sys.getenv("POSTDOC_ROOT",
                            "/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc"),
                 "OCM_heart/allelic_ratio/00_functions.R"))

SOURCE_TABLE <- Sys.getenv("SOURCE_TABLE",
                           "Allelic_ratio_results/whole_chr_allelic_ratios.txt")
SAMPLES <- strsplit(Sys.getenv("SAMPLES", "9w,78w,Sham,TAC"), "[, ]+")[[1]]

if (!file.exists(SOURCE_TABLE)) {
  stop("SOURCE_TABLE does not exist: ", SOURCE_TABLE)
}
# Refuse to read and write the same path: with RESULTS_ROOT left at its default
# this would overwrite the only surviving copy of the ingested table, and the
# tree it came from is gone.
if (normalizePath(SOURCE_TABLE) ==
      normalizePath(ALLELIC_RATIOS_FILE, mustWork = FALSE)) {
  stop("SOURCE_TABLE and the output are the same file:\n  ", SOURCE_TABLE,
       "\n  Set RESULTS_ROOT to a new directory. The Allelome.PRO2 tree this ",
       "table came from is no longer on disk, so this file is not regenerable ",
       "without rebuilding it.")
}

message("Source:       ", SOURCE_TABLE)
message("RESULTS_ROOT: ", RESULTS_ROOT)
message("Output:       ", ALLELIC_RATIOS_FILE)

raw <- read.delim(SOURCE_TABLE, header = TRUE, stringsAsFactors = FALSE)
need <- c("cell_barcode", "chr", "A1_reads", "A2_reads")
if (!all(need %in% names(raw))) {
  stop("SOURCE_TABLE is missing ", paste(setdiff(need, names(raw)), collapse = ", "),
       "; has: ", paste(names(raw), collapse = ", "))
}
message("Read ", nrow(raw), " rows, ", length(unique(raw$cell_barcode)), " cells")

# `sample` from the <sample>_<barcode> key. The four sample names carry no
# underscore, so the first field is the sample; anything else means the key is
# not what the rest of the pipeline assumes and is worth stopping over.
raw$sample <- sub("_.*$", "", raw$cell_barcode)
unexpected <- setdiff(unique(raw$sample), SAMPLES)
if (length(unexpected)) {
  stop("Unexpected sample prefixes in cell_barcode: ",
       paste(head(unexpected, 5), collapse = ", "),
       "\n  Expected one of: ", paste(SAMPLES, collapse = ", "))
}

raw <- drop_doublets(raw, paste0("12 (", basename(SOURCE_TABLE), ")"))

# One row per (cell, chr), summing whatever intervals the annotation defined.
# Idempotent on a table that already has one row per (cell, chr), and correct
# on one scored against a per-gene bed -- the failure mode
# assert_one_row_per_cell() exists to catch downstream.
df <- raw %>%
  dplyr::group_by(cell_barcode, sample, chr) %>%
  dplyr::summarise(A1_reads = sum(A1_reads), A2_reads = sum(A2_reads),
                   n_interval = dplyr::n(), .groups = "drop") %>%
  as.data.frame()
if (max(df$n_interval) > 1) {
  message("Note: up to ", max(df$n_interval), " intervals per (cell, chr) were ",
          "summed -- the source was scored against a per-gene bed.")
}
df$n_interval <- NULL

df$total_reads <- df$A1_reads + df$A2_reads
df <- df[df$total_reads > 0, ]
df$ar_a1  <- df$A1_reads / df$total_reads
df$ar_dom <- pmax(df$A1_reads, df$A2_reads) / df$total_reads

# How much of the change downstream is the convention fix rather than the
# doublet removal. These are the cells whose allelic_ratio moves by more than
# rounding because ar_a1 < 0.5 is being replaced by 1 - ar_a1.
x <- df[df$chr == "chrX", ]
flipped <- sum(x$ar_a1 < 0.5)
message(sprintf("chrX cells: %d; convention fix changes %d of them (%.2f%%) -- ",
                nrow(x), flipped, 100 * flipped / max(nrow(x), 1)),
        "these had ar_a1 < 0.5 and were read as biallelic.")

auto <- collapse_autosomes(df)
message("Autosomal rows added: ", nrow(auto), " cells")

out <- bind_rows(df, auto)
# 02 and 06 read a column literally named allelic_ratio and treat >= MONO_AR as
# monoallelic, i.e. the dominant-allele convention -- same as
# 10_build_ratio_table.R. Both forms are kept alongside it.
out$allelic_ratio <- out$ar_dom
out <- out[, c("cell_barcode", "sample", "chr", "A1_reads", "A2_reads",
               "total_reads", "ar_a1", "ar_dom", "allelic_ratio")]

write.table(out, ALLELIC_RATIOS_FILE, sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote ", ALLELIC_RATIOS_FILE, "  (", nrow(out), " rows, ",
        length(unique(out$cell_barcode)), " cells)")
