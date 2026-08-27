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
# duplicate and a further 7.4% are not uniquely mapped, and a nucleus carries far
# fewer informative chrX molecules than a bulk library, so this is not small
# print.
#
# scAllelome_dedup.slurm rescores the same cells after -F 0x400 -q 255. This
# joins the two trees on cell barcode and asks what moved.
#
# Run from the OCM_heart/ directory, in seurat_env (NOT RNAseq):
#   Rscript allelic_ratio/09_dedup_comparison.R
# ---------------------------------------------------------------------------
source("/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/OCM_heart/allelic_ratio/00_functions.R")

SAMPLES  <- c("9w", "78w", "Sham", "TAC")
TREE_RAW <- Sys.getenv("TREE_RAW", "Allelome.PRO2")
TREE_DED <- Sys.getenv("TREE_DED", "Allelome.PRO2_dedup")
OUT_DIR  <- file.path(RESULTS_ROOT, "dedup_comparison")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Which samples can be compared at all. Two trees only pair validly where the
# same annotation bed was used, and that has not always been true: scAllelome.slurm
# scored 9w against a chrX-only per-gene bed while the other three got a
# whole-chromosome one. Detect it from the run directory names rather than
# hard-coding a sample list, so this stays right after the read-level 9w tree is
# rebuilt (FILTER=0 sbatch scAllelome_dedup.slurm) and there is nothing to
# remember to update.
annot_raw <- allelome_annotations(TREE_RAW, SAMPLES)
annot_ded <- allelome_annotations(TREE_DED, SAMPLES)
annot <- full_join(annot_raw, annot_ded, by = "sample",
                   suffix = c("_raw", "_ded"))
print(as.data.frame(annot))

multi <- annot$sample[duplicated(annot$sample)]
if (length(multi)) {
  stop("More than one annotation in a tree for: ", paste(unique(multi), collapse = ", "),
       "\n  Move the older runs aside; mixing them gives a cell rows from both beds.")
}
CONFOUNDED <- annot$sample[is.na(annot$annotation_raw) |
                           is.na(annot$annotation_ded) |
                           annot$annotation_raw != annot$annotation_ded]
if (length(CONFOUNDED)) {
  message("Annotation differs between trees for: ", paste(CONFOUNDED, collapse = ", "),
          " -- reported per sample but excluded from pooled statistics, since a ",
          "paired difference there is not attributable to deduplication.")
} else {
  message("All samples share an annotation between trees; pooling all of them.")
}

# The cutoff applies to the READ-level tree only, and this is the one place it
# must not be symmetric: gating the deduplicated tree at the same number would
# discard exactly the cells whose read count was inflated, which is the effect
# being measured. Cells are chosen on reads, then followed into molecules
# whatever their molecule count turns out to be.
message("Selecting cells on read-level total_reads >= ", MIN_TOTAL_READS)

load_chrX <- function(tree) {
  df <- load_allelome_tree(tree, SAMPLES)
  df <- df[df$chr == "chrX", ]
  if (!nrow(df)) stop("No chrX rows under ", tree)
  df$AR <- df$ar_dom       # dominant-allele fraction, the convention 02/08 use
  df[, c("cell_barcode", "sample", "A1_reads", "A2_reads", "total_reads", "AR")]
}

raw <- load_chrX(TREE_RAW)
ded <- load_chrX(TREE_DED)
message("read-level cells: ", nrow(raw), " | molecule-level cells: ", nrow(ded))

cmp <- inner_join(raw, ded, by = c("cell_barcode", "sample"),
                  suffix = c("_raw", "_ded"))
cmp <- cmp[cmp$total_reads_raw >= MIN_TOTAL_READS, ]
stopifnot(nrow(cmp) > 0)
cmp$retained  <- cmp$total_reads_ded / cmp$total_reads_raw
cmp$AR_shift  <- cmp$AR_ded - cmp$AR_raw
cmp$confounded <- cmp$sample %in% CONFOUNDED

message("paired cells: ", nrow(cmp),
        " | median read retention: ", round(median(cmp$retained), 3))

clean <- cmp[!cmp$confounded, ]
message("pooled over ", nrow(clean), " cells",
        if (length(CONFOUNDED)) paste0(" (excluding ", paste(CONFOUNDED, collapse = ", "), ")") else "")
stopifnot(nrow(clean) > 0)

# ---------------------------------------------------------------------------
# The question 08 asked. If the extreme-AR population is duplication-driven, it
# should relax toward 0.5 once each molecule is counted once. AR here is the
# dominant-allele fraction, so it is bounded in [0.5, 1] and there is only one
# extreme tail.
# ---------------------------------------------------------------------------
extreme <- clean[clean$AR_raw >= 0.90, ]
message(sprintf("cells at AR >= 0.90 (read level): %d / %d (%.1f%%)",
                nrow(extreme), nrow(clean), 100 * nrow(extreme) / nrow(clean)))
if (nrow(extreme)) {
  still <- sum(extreme$AR_ded >= 0.90)
  message(sprintf("  still >= 0.90 after deduplication: %d (%.1f%%)",
                  still, 100 * still / nrow(extreme)))
  message(sprintf("  median AR change: %+.3f", median(extreme$AR_shift)))
  # A cell whose molecule count is too low to say anything is not evidence
  # either way, so separate "moved" from "ran out of depth".
  message(sprintf("  of those, %d fall below %d molecules after deduplication",
                  sum(extreme$total_reads_ded < MIN_TOTAL_READS), MIN_TOTAL_READS))
}

per_sample <- cmp %>%
  group_by(sample) %>%
  summarise(n            = n(),
            confounded   = any(confounded),
            med_retained = median(retained),
            med_AR_raw   = median(AR_raw),
            med_AR_ded   = median(AR_ded),
            frac_ext_raw = mean(AR_raw >= 0.90),
            frac_ext_ded = mean(AR_ded >= 0.90),
            .groups = "drop")
print(as.data.frame(per_sample))

write.table(cmp, file.path(OUT_DIR, "paired_cell_ratios.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(per_sample, file.path(OUT_DIR, "per_sample_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

p_pair <- ggplot(cmp, aes(AR_raw, AR_ded)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_point(aes(colour = total_reads_ded), alpha = 0.4, size = 0.8) +
  scale_colour_viridis_c(trans = "log10", name = "molecules") +
  facet_wrap(~ sample) +
  coord_fixed(xlim = c(0.5, 1), ylim = c(0.5, 1)) +
  labs(title = "chrX allelic ratio, read level vs molecule level",
       subtitle = sprintf("cells with >= %d reads before deduplication; points below the diagonal moved toward 0.5%s",
                          MIN_TOTAL_READS,
                          if (length(CONFOUNDED)) paste0(" (", paste(CONFOUNDED, collapse = ", "), ": annotation-confounded)") else ""),
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
