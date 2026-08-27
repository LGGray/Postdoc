# ---------------------------------------------------------------------------
# 09 - Read-level vs molecule-level allelic ratios.
#
# Does deduplication change the chrX allelic ratio, and specifically does it
# dissolve the AR >= 0.90 population that 08 is trying to interpret?
#
# Allelome.PRO2 counts reads. Every read of one UMI carries the same allele by
# construction, so PCR duplicates add weight to a cell's ratio without adding
# evidence about which allele was present. In expectation that is unbiased; what
# it destroys is n. A cell with three real molecules and heavy amplification on
# one of them reads as near-monoallelic. 39.4% of reads in this data are flagged
# duplicate, and a cell carries far fewer informative chrX molecules than a bulk
# library, so this is not a small-print correction.
#
# scAllelome_dedup.slurm rescores the same cells after -F 0x400 -q 255. This
# script joins the two trees on cell barcode and asks what moved.
#
# Run from the OCM_heart/ directory, in seurat_env (NOT RNAseq):
#   Rscript allelic_ratio/09_dedup_comparison.R
# ---------------------------------------------------------------------------
source("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/OCM_heart/allelic_ratio/00_functions.R")

SAMPLES  <- c("9w", "78w", "Sham", "TAC")
TREE_RAW <- "Allelome.PRO2"
TREE_DED <- "Allelome.PRO2_dedup"
OUT_DIR  <- "Allelic_ratio_results/dedup_comparison"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# The cutoff applies to the READ-level tree only, and is the one thing here that
# must not be applied symmetrically: gating the deduplicated tree at the same
# number would silently discard exactly the cells whose read count was inflated,
# which is the effect being measured. Cells are selected on the read-level tree
# and then followed into the deduplicated one whatever their molecule count.
message("Selecting cells on read-level total_reads >= ", MIN_TOTAL_READS)

# ---------------------------------------------------------------------------
# Allelome.PRO2 names its output <date>_<bam>_<annotation>_<run no.>, so the
# cell barcode has to be recovered from the directory name. Do NOT do this by
# splitting the full path on "_" and taking a fixed index, the way 02 and the
# all-genes block in Allelic_ratio_testing.R do - that index is a function of
# how many underscores are in the tree's own directory name, so it is 4/5 for
# Allelome.PRO2/ and 6/7 for Allelome.PRO2_all_genes/. Parse the basename.
# ---------------------------------------------------------------------------
parse_cell <- function(paths) {
  d <- basename(dirname(paths))
  d <- sub("^[0-9]{4}_[0-9]{2}_[0-9]{2}_", "", d)   # drop the date stamp
  sub("_[^_]*\\.bed_[0-9]+$", "", d)                # drop _<annotation>_<run>
}

load_tree <- function(tree) {
  paths <- unlist(lapply(SAMPLES, function(s) {
    list.files(file.path(tree, s), pattern = "locus_table.txt",
               recursive = TRUE, full.names = TRUE)
  }))
  if (!length(paths)) stop("No locus tables under ", tree)

  sample_of <- sub(paste0("^", tree, "/([^/]+)/.*$"), "\\1", paths)
  cell_of   <- parse_cell(paths)

  rows <- lapply(seq_along(paths), function(i) {
    tmp <- tryCatch(read.delim(paths[i], header = TRUE, stringsAsFactors = FALSE),
                    error = function(e) NULL)
    if (is.null(tmp) || !nrow(tmp) || !"chr" %in% names(tmp)) return(NULL)
    tmp <- tmp[as.character(tmp$chr) == "chrX", , drop = FALSE]
    if (!nrow(tmp)) return(NULL)
    data.frame(cell_barcode = paste0(sample_of[i], "_", cell_of[i]),
               sample       = sample_of[i],
               A1_reads     = sum(tmp$A1_reads),
               A2_reads     = sum(tmp$A2_reads),
               stringsAsFactors = FALSE)
  })
  out <- bind_rows(rows)
  # Ratio recomputed from the summed counts rather than averaging the per-locus
  # allelic_ratio column, so it is depth-weighted the same way in both trees.
  out$total_reads <- out$A1_reads + out$A2_reads
  out$AR <- ifelse(out$total_reads > 0, out$A1_reads / out$total_reads, NA_real_)
  out[!is.na(out$AR), ]
}

raw <- load_tree(TREE_RAW)
ded <- load_tree(TREE_DED)
message("read-level cells: ", nrow(raw), " | molecule-level cells: ", nrow(ded))

cmp <- inner_join(raw, ded, by = c("cell_barcode", "sample"),
                  suffix = c("_raw", "_ded"))
cmp <- cmp[cmp$total_reads_raw >= MIN_TOTAL_READS, ]
cmp$retained   <- cmp$total_reads_ded / cmp$total_reads_raw
cmp$AR_shift   <- cmp$AR_ded - cmp$AR_raw
cmp$to_centre  <- abs(cmp$AR_ded - 0.5) - abs(cmp$AR_raw - 0.5)
stopifnot(nrow(cmp) > 0)

message("paired cells: ", nrow(cmp),
        " | median read retention: ", round(median(cmp$retained), 3))

# ---------------------------------------------------------------------------
# The question 08 asked. If the extreme-AR population is duplication-driven, it
# should relax toward 0.5 once each molecule is counted once.
# ---------------------------------------------------------------------------
extreme <- cmp[cmp$AR_raw >= 0.90 | cmp$AR_raw <= 0.10, ]
message(sprintf("cells at AR <= 0.10 or >= 0.90 (read level): %d / %d (%.1f%%)",
                nrow(extreme), nrow(cmp), 100 * nrow(extreme) / nrow(cmp)))
if (nrow(extreme)) {
  still <- sum(extreme$AR_ded >= 0.90 | extreme$AR_ded <= 0.10)
  message(sprintf("  still extreme after deduplication: %d (%.1f%%)",
                  still, 100 * still / nrow(extreme)))
  message(sprintf("  median |AR - 0.5| change: %+.3f", median(extreme$to_centre)))
}

# A cell whose molecule count is too low to say anything is not evidence either
# way, so report how many of the extreme cells simply ran out of depth.
if (nrow(extreme)) {
  message(sprintf("  of those, %d fall below %d molecules after deduplication",
                  sum(extreme$total_reads_ded < MIN_TOTAL_READS), MIN_TOTAL_READS))
}

per_sample <- cmp %>%
  group_by(sample) %>%
  summarise(n            = n(),
            med_retained = median(retained),
            med_AR_raw   = median(AR_raw),
            med_AR_ded   = median(AR_ded),
            frac_ext_raw = mean(AR_raw >= 0.90 | AR_raw <= 0.10),
            frac_ext_ded = mean(AR_ded >= 0.90 | AR_ded <= 0.10),
            .groups = "drop")
print(as.data.frame(per_sample))

write.table(cmp, file.path(OUT_DIR, "paired_cell_ratios.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(per_sample, file.path(OUT_DIR, "per_sample_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

p_pair <- ggplot(cmp, aes(AR_raw, AR_ded)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey70") +
  geom_point(aes(colour = total_reads_ded), alpha = 0.4, size = 0.8) +
  scale_colour_viridis_c(trans = "log10", name = "molecules") +
  facet_wrap(~ sample) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "chrX allelic ratio, read level vs molecule level",
       subtitle = sprintf("cells with >= %d reads before deduplication; points below the diagonal on the right moved toward 0.5",
                          MIN_TOTAL_READS),
       x = "AR from all reads", y = "AR from deduplicated, unique reads") +
  theme_bw()
ggsave(file.path(OUT_DIR, "AR_paired_scatter.pdf"), p_pair, width = 9, height = 8)

p_ret <- ggplot(cmp, aes(retained)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = median(cmp$retained), colour = "red") +
  facet_wrap(~ sample, scales = "free_y") +
  labs(title = "Fraction of chrX reads surviving -F 0x400 -q 255",
       subtitle = "red = median. A long left tail means cells whose ratio was set by a few amplified molecules.",
       x = "molecules / reads", y = "cells") +
  theme_bw()
ggsave(file.path(OUT_DIR, "read_retention.pdf"), p_ret, width = 9, height = 6)

message("Wrote ", OUT_DIR)
